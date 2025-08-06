#include "Wake.h"
#include "Beam.h"

Wake::Wake()
{
  loss = 0;
  lossref="";
  radius=2.5e-3;
  conductivity=0;
  relaxation=0;
  roundpipe=true;
  transient=false;
  ztrans=0;
  gap=0;
  lgap=1;
  hrough=0;
  lrough=1;
  ns=0;
}

Wake::~Wake(){}


void Wake::usage(){

  cout << "List of keywords for Wake" << endl;
  cout << "&wake" << endl;
  cout << " string loss =  0 / reference" << endl;
  cout << " double radius = 2.5e-3" << endl;
  cout << " bool   roundpipe   = true" << endl;
  cout << " string material  = <empty>" << endl;
  cout << " double conductivity = 0e-6" << endl;
  cout << " double relaxation  = 0e-6" << endl;
  cout << " double gap = 0e-6" << endl;
  cout << " double lgap  = 1.0" << endl;
  cout << " double hrough = 0.0" << endl;
  cout << " double lrough = 1.0" << endl;
  cout << " bool transient = false" << endl;
  cout << " double ztrans = 0" << endl;
  cout << "&end" << endl << endl;
  return;
}

// input parameter

bool Wake::init(int rank, int size, map<string,string> *arg,  Time *time, Setup *setup, Beam *beam, Profile *prof)
{

  string material="";
  map<string,string>::iterator end=arg->end();
  map<string,string>::iterator iter=arg->begin();
  
  if (arg->find("loss")!=end    ){this->reference(arg->at("loss"),&loss,&lossref); arg->erase(arg->find("loss"));}
  if (arg->find("radius")!=end) {radius = atof(arg->at("radius").c_str());  arg->erase(arg->find("radius"));}
  if (arg->find("conductivity")!=end) {conductivity= atof(arg->at("conductivity").c_str());  arg->erase(arg->find("conductivity"));}
  if (arg->find("relaxation")!=end) {relaxation = atof(arg->at("relaxation").c_str());  arg->erase(arg->find("relaxation"));}
  if (arg->find("roundpipe")!=end)    {roundpipe    = atob(arg->at("roundpipe").c_str());  arg->erase(arg->find("roundpipe"));}
  if (arg->find("material")!=end){material = arg->at("material"); arg->erase(arg->find("material"));}
  if (arg->find("gap")!=end)  {gap = atof(arg->at("gap").c_str());  arg->erase(arg->find("gap"));}
  if (arg->find("lgap")!=end) {lgap = atof(arg->at("lgap").c_str());  arg->erase(arg->find("lgap"));}
  if (arg->find("hrough")!=end) {hrough = atof(arg->at("hrough").c_str());  arg->erase(arg->find("hrough"));}
  if (arg->find("lrough")!=end) {lrough = atof(arg->at("lrough").c_str());  arg->erase(arg->find("lrough"));}
  if (arg->find("transient")!=end)    {transient    = atob(arg->at("transient").c_str());  arg->erase(arg->find("transient"));}
  if (arg->find("ztrans")!=end) {ztrans = atof(arg->at("ztrans").c_str());  arg->erase(arg->find("ztrans"));}


  if (arg->size()!=0){
    if (rank==0){ cout << "*** Error: Unknown elements in &wake" << endl; this->usage();}
    return false;
  }

  // check for external wake function
  
 string wrongProf="";
  if (prof->check(lossref)== false)  { wrongProf=lossref;}
  if (wrongProf.size() > 0){
    if (rank==0){cout << "*** Error: Unknown profile reference in &wake: " << wrongProf << endl;}
    return false;
  }    

  if (rank==0){cout << "Generating wakefield potentials..." << endl; }

  vector<double> s;
  ns=time->getPosition(&s);  // total number of slices and vector with the s-position
  int nsNode=time->getNodeNSlice();
  
  ns=ns*time->getSampleRate();
  ds=(s[1]-s[0])/time->getSampleRate();

  wakeres = new double [ns];
  wakegeo = new double [ns];
  wakerou = new double [ns];
  wakeext = new double [nsNode];

  for (int j=0; j<nsNode; j++){
    int i=j+time->getNodeOffset();
    wakeext[j]=prof->value(s[i],loss,lossref);  // external wake
  }


  if ((material=="CU") || (material=="Cu") || (material=="cu")){
    conductivity=5.813e7;
    relaxation=8.1e-6;
  }

  if ((material=="AL") || (material=="Al") || (material=="al")){
    conductivity=3.571e7;
    relaxation=2.4e-6;
  }

  // calculate the single particle wakes
  this->singleWakeResistive(rank);
  this->singleWakeGeometric(rank);
  this->singleWakeRoughness(rank);


  // transfer wakes into beam class
  beam->initWake(ns,nsNode, ds,wakeext,wakeres,wakegeo,wakerou,ztrans,radius,transient);

  delete[] wakeres;
  delete[] wakegeo;
  delete[] wakerou;
  delete[] wakeext;

  return true;

}



void Wake::KernelRoughness(vector< complex<double> > *K, complex<double> q1, complex<double> q2)
{
  int N=K->size();
  complex<double> dq=(q2-q1)/static_cast<double>(N-1);
  complex<double> i=complex<double> (0,1);

  for (int j=0; j< N;j++){
    complex<double> q=q1+static_cast<double>(j)*dq;
    complex<double> S=(sqrt(2.*q+1.)-i*sqrt(2.*q-1.))*q/sqrt(4.*q*q-1.);
    K->at(j)=(S+1.)/(1.-i*rrough*q*S)/(1.+i*rrough*q);
  }
  K->at(0)*=0.5;
  K->at(N-1)*=0.5;   // for trapazoidal integration
}



double Wake::TrapIntegrateRoughness(vector< complex<double> > *K, complex<double> q1, complex<double> q2, double tau)
{
  int N=K->size();
  complex<double> dq=(q2-q1)/static_cast<double>(N-1);
  complex<double> i=complex<double> (0,1);  
  complex<double> val=0;

  for (int j=0; j< N;j++){
    complex<double> q=q1+static_cast<double>(j)*dq;
    val+=exp(-i*q*tau)*K->at(j);
  }
  val*=dq;
  return val.real();
}





void Wake::singleWakeRoughness(int rank)
{

  if (hrough <=0){
    for (int i=0; i< ns; i++){
      wakerou[i]=0;
    }
    return;
  }


  double pi=2.*asin(1);
  rrough = pi*pi*pi/lrough/lrough/lrough*hrough*hrough*radius; // aspect ratio parameter. If much smaller than one an inductive mode
             // is expected, if unity or larger it is rather a synchronous mode.

  double tau=0;  // example 
  complex<double> q1 = complex<double> (0,0);
  complex<double> q2 = complex<double> (0,2e-3);
  complex<double> q3 = complex<double> (1,2e-3);
  complex<double> q4 = complex<double> (1,0);
  complex<double> q5 = complex<double> (100,0);

  int N=128;
  vector< complex<double> > K1,K2,K3,K4;
  K1.resize(N);
  K2.resize(8*N);
  K3.resize(N);
  K4.resize(8*N);

  this->KernelRoughness(&K1,q1,q2);
  this->KernelRoughness(&K2,q2,q3);
  this->KernelRoughness(&K3,q3,q4);
  this->KernelRoughness(&K4,q4,q5);


  double res;
  double coef=rrough/pi*4/radius/radius*1.6e-19/4/pi/8.854e-12;   

  for (int i=0; i<ns; i++){
    tau=2*pi*ds*i/lrough;
    res =this->TrapIntegrateRoughness(&K1,q1,q2,tau);
    res+=this->TrapIntegrateRoughness(&K2,q2,q3,tau);
    res+=this->TrapIntegrateRoughness(&K3,q3,q4,tau);
    res+=this->TrapIntegrateRoughness(&K4,q4,q5,tau);
    wakerou[i]=coef*res;
  }
  
  if (rank==0){
       cout << "Roughness Wake calculated..." << endl;   
  }

 
}


void Wake::singleWakeGeometric(int rank)
{


  for (int i=0; i< ns; i++){
      wakegeo[i]=0;
  }

  if (gap <=0){
    return;
  }

  double pi=2.*asin(1.);
  double coef=-vacimp*ce/(pi*pi*radius*lgap)*2*sqrt(0.5*gap); // scaling coefficient
  if (!roundpipe) { coef*=0.956; }     //
  for (int i = 0;i<ns;i++){
    wakegeo[i]=coef*sqrt(ds*i);     
  }
   
  if (rank==0){
       cout << "Geometric Wake calculated..." << endl;   
  }

}




 
/**
 * @brief Compute the longitudinal resistive-wall wake w(s) for round or flat chambers.
 *
 * Physics source for all formulas: Bane & Stupakov, SLAC-PUB-10707 (rev. Oct. 2004).
 *
 * Key relations used here (see cited equation numbers in the paper):
 *   • AC (Drude) conductivity in k-space (Eq. (1)):
 *       σ̃(k) = σ0 / (1 − i k c τ),  with ctau = c τ stored in member `relaxation` [m].
 *
 *   • Surface impedance (good conductor, Leontovich; re-expression of the wall physics):
 *       ζ(k) = (1 − i) * sqrt( k / ( 2 σ̃(k) Z0 ) )   [SI units].
 *
 *   • Round (circular) pipe (radius a), m=0 Leontovich form (equivalent to Eqs. (2)–(3)):
 *       Z(k) = [ Z0 / (2π a) ] / [ 1/ζ(k) − i k a / 2 ]        [Ω/m].
 *
 *   • Flat (parallel plates, half-gap a) kernel integral (Eq. (13), recast via ζ(k) and x=qa):
 *       Z(k) = ∫_0^∞  [ Z0/(2π a) ] /
 *               { cosh x [ cosh x / ζ(k) − i k a sinhc(x) ] }  dx
 *       where  sinhc(x) = sinh(x)/x  (evaluated with a stable series near x=0).
 *
 *   • Wake from impedance (inverse FT; one-sided cosine consistent with W(s<0)=0):
 *       w(s) = (2c/π) ∫_0^{kmax} Re Z(k) cos(k s) dk .
 *
 * Numerics:
 *   • Build Z(k) on a uniform k-grid, k ∈ [0, kmax], kmax = κmax/s0 with κmax≈100.
 *   • Flat case: trapezoid in x ∈ [0, XMAX] with endpoint ½-weights.
 *   • Cosine transform: single block for both geometries; endpoint ½-weights in k.
 *   • Final: multiply by (−e) once (converts V/C to signed energy-kick units for electrons).
 *
 * Implementation notes:
 *   • Prints "k  Re(Zk)  Im(Zk)" on rank==0 for quick cross-checks; comment out later.
 *   • Tunables: KAPPA_MAX, NK (k-grid); XMAX, NQ (flat x-integral).
 */
void Wake::singleWakeResistive(int rank)
{
    // Zero output; early exit if no conductivity
    for (int i = 0; i < ns; ++i) wakeres[i] = 0.0;
    if (conductivity <= 0.0) return;

    // ---- Constants & characteristic scales ----
    const double c_light = 299792458.0;                // m/s
    const double Z0      = vacimp;                     // Ω  (vacuum impedance, SI)
    const double a       = radius;                     // m  (round radius or flat half-gap)
    const double s0      = std::pow(2.0 * a * a / (Z0 * conductivity), 1.0 / 3.0); // (2 a^2 / Z0 σ0)^(1/3)
    const double pi      = 2.0 * std::asin(1.0);

    // k-grid (shared)
    const double       KAPPA_MAX = 100.0;              // κ_max
    const unsigned int NK        = 1000;               // # of k-intervals → NK+1 samples
    const double       kmax      = KAPPA_MAX / s0;     // 1/m
    const double       dk        = kmax / static_cast<double>(NK);

    // Flat Z(k): x-integral controls (∫_0^∞ → [0, XMAX] with trapezoid)
    const double XMAX = 15.0;
    const int    NQ   = 20000;

    // Final −e factor (V/C → signed energy-kick units)
    constexpr double e_charge = 1.602176634e-19;       // C

    using cd = std::complex<double>;

    // ---- Helpers: stable sinhc, Drude σ̃(k), surface impedance ζ(k) ----

    // Numerically stable sinhc(x) = sinh(x)/x with a small-x series
    auto sinhc_stable = [](double x) {
        const double ax = std::abs(x);
        if (ax < 1e-5) {                               // ≈ 1 + x^2/6 + x^4/120
            const double x2 = x * x;
            return 1.0 + x2 * (1.0/6.0 + x2 * (1.0/120.0));
        }
        return std::sinh(x) / x;
    };

    // σ̃(k) = σ0 / (1 − i k c τ)   (SLAC-PUB-10707, Eq. (1))
    auto sigma_k = [&](double k) -> cd {
        return conductivity / (cd(1.0, 0.0) - cd(0.0, 1.0) * k * relaxation);
    };

    // ζ(k) = (1 − i) sqrt( k / ( 2 σ̃(k) Z0 ) )   (Leontovich surface impedance, SI form)
    auto zeta_k = [&](double k) -> cd {
        const cd arg  = cd(k, 0.0) / (cd(2.0, 0.0) * sigma_k(k) * cd(Z0, 0.0));
        const cd root = std::sqrt(arg);
        return cd(1.0, -1.0) * root;
    };

    // ---- Geometry-specific Z(k) ----

    // Flat (parallel plates) integrand recast in ζ(k) and x = q a.
    // Equivalent in structure to SLAC-PUB-10707 Eq. (13) after variable change and SI units.
    auto Zint_flat = [&](double k, double x) -> cd {
        const double C   = std::cosh(x);
        const double Shc = sinhc_stable(x);
        const cd z       = zeta_k(k);

        // denom = cosh(x) * ( cosh(x)/ζ  −  i k a sinhc(x) )
        const cd term1 = cd(C, 0.0) / z;
        const cd term2 = cd(0.0, -k * a * Shc);
        const cd denom = cd(C, 0.0) * (term1 + term2);

        // prefactor = Z0 / (2π a)   [Ω/m]
        const double pref = Z0 / (2.0 * pi * a);
        if (!std::isfinite(denom.real()) || !std::isfinite(denom.imag())) return cd(0.0, 0.0);
        return cd(pref, 0.0) / denom;
    };

    // Flat Z(k) via trapezoid in x ∈ [0, XMAX]
    auto Zk_flat = [&](double k) -> cd {
        if (k == 0.0) return cd(0.0, 0.0);
        const double dx = XMAX / static_cast<double>(NQ - 1);
        cd acc = Zint_flat(k, 0.0) * 0.5;              // endpoint ½-weights
        for (int j = 1; j < NQ - 1; ++j) acc += Zint_flat(k, dx * static_cast<double>(j));
        acc += Zint_flat(k, XMAX) * 0.5;
        return acc * dx;
    };

    // Round (circular) Z(k), Leontovich m=0 (equivalent to the Eqs. (2)–(3) form in 10707)
    // Z(k) = [ Z0/(2π a) ] / [ 1/ζ(k) − i k a / 2 ]
    auto Zk_round = [&](double k) -> cd {
        if (k == 0.0) return cd(0.0, 0.0);
        const cd z   = zeta_k(k);
        const cd den = cd(1.0, 0.0) / z - cd(0.0, 0.5 * k * a);
        const double pref = Z0 / (2.0 * pi * a);
        if (!std::isfinite(den.real()) || !std::isfinite(den.imag())) return cd(0.0, 0.0);
        return cd(pref, 0.0) / den;
    };

    // ---- Build Re Z(k) on a uniform k-grid ----
    std::vector<double> ks(NK + 1), ReZ(NK + 1);
    for (unsigned int j = 0; j <= NK; ++j) ks[j] = dk * static_cast<double>(j);

    if (rank == 0) std::cout << std::setprecision(16);

    if (roundpipe) {
        for (unsigned int j = 0; j <= NK; ++j) {
            const double k = ks[j];
            const cd Z = Zk_round(k);
            ReZ[j] = Z.real();
            if (rank == 0) std::cout << k << " " << Z.real() << " " << Z.imag() << "\n";
        }
    } else {
        for (unsigned int j = 0; j <= NK; ++j) {
            const double k = ks[j];
            const cd Z = Zk_flat(k);
            ReZ[j] = Z.real();
            if (rank == 0) std::cout << k << " " << Z.real() << " " << Z.imag() << "\n";
        }
    }

    // ---- One-sided cosine transform (shared) ----
    // As stated in SLAC-PUB-10707, the wake is the inverse FT of Z(k);
    // with W(s<0)=0 this gives the one-sided cosine form used here.
    ReZ.front() *= 0.5;                                // endpoint ½-weights
    ReZ.back()  *= 0.5;
    const double cos_pref = (2.0 * c_light / pi) * dk; // (2c/π) Δk

    for (int i = 0; i < ns; ++i) {
        const double s_i = i * ds;
        double acc = 0.0;
        for (unsigned int j = 0; j <= NK; ++j) {
            acc += ReZ[j] * std::cos(ks[j] * s_i);
        }
        wakeres[i] = acc * cos_pref;                   // V/C
    }

    // ---- Final conversion: V/C → signed energy kick (electron sign) ----
    for (int i = 0; i < ns; ++i) wakeres[i] *= -e_charge;

    if (rank == 0) {
        std::cout << "Resistive Wake (" << (roundpipe ? "ROUND" : "FLAT")
                  << ") calculated (s0 = " << s0 << ")..." << std::endl;
    }
}

