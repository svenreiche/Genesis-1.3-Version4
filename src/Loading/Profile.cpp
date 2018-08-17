#include "Profile.h"


Profile::Profile()
{
}

Profile::~Profile()
{
  prof.clear();
}

bool Profile::init(int rank, map<string,string> *arg,string element)
{

  ProfileBase *p;
  string label;

  if (element.compare("&profile_const")==0){
    p=(ProfileBase *)new ProfileConst();
    label=p->init(rank,arg);
  } 
  if (element.compare("&profile_gauss")==0){
    p=(ProfileBase *)new ProfileGauss();
    label=p->init(rank,arg);
  } 
  if (element.compare("&profile_polynom")==0){
    p=(ProfileBase *)new ProfilePolynom();
    label=p->init(rank,arg);
  } 
  if (element.compare("&profile_step")==0){
    p=(ProfileBase *)new ProfileStep();
    label=p->init(rank,arg);
  } 
  if (element.compare("&profile_file")==0){
    p=(ProfileBase *)new ProfileFile();
    label=p->init(rank,arg);
  } 

  if (label.size()<1){
    return false;
  } else {
    prof[label]=p;
  }
  if (rank==0) {cout << "Adding profile with label: " << label << endl;}
  return true;
}

double Profile::value(double s, double val, string label)
{
  if ((label.size()<1)||(prof.find(label)==prof.end())){  
    return val;
  } else {
    return  prof[label]->value(s);
  }
}






//------------------------------------
// individual profiles


string ProfileConst::init(int rank, map<string,string>*arg)
{

  string label="";
  c0=0;
  map<string,string>::iterator end=arg->end();

  if (arg->find("label")!=end){label = arg->at("label");  arg->erase(arg->find("label"));}
  if (arg->find("c0")!=end)   {c0    = atof(arg->at("c0").c_str());  arg->erase(arg->find("c0"));}

  if (arg->size()!=0){
    if (rank==0){ cout << "*** Error: Unknown elements in &profile_const" << endl; this->usage();}
    return "";
  }
  if ((label.size()<1)&&(rank==0)){
    cout << "*** Error: Label not defined in &profile_const" << endl; this->usage();
  }
  return label;
}

double ProfileConst::value(double z)
{
  return c0;
}

void ProfileConst::usage(){
  cout << "List of keywords for PROFILE_CONST" << endl;
  cout << "&profile_const" << endl;
  cout << " string label = <empty>" << endl;
  cout << " double c0 = 0" << endl;
  cout << "&end" << endl << endl;
  return;
}

//-----------------------

string ProfilePolynom::init(int rank, map<string,string>*arg)
{
  string label="";
  map<string,string>::iterator end=arg->end();
  c.resize(5);
  for (int i=0; i< c.size();i++){ c[i]=0;}
  

  if (arg->find("label")!=end){label = arg->at("label");  arg->erase(arg->find("label"));}
  if (arg->find("c0")!=end)   {c[0]    = atof(arg->at("c0").c_str());  arg->erase(arg->find("c0"));}
  if (arg->find("c1")!=end)   {c[1]    = atof(arg->at("c1").c_str());  arg->erase(arg->find("c1"));}
  if (arg->find("c2")!=end)   {c[2]    = atof(arg->at("c2").c_str());  arg->erase(arg->find("c2"));}
  if (arg->find("c3")!=end)   {c[3]    = atof(arg->at("c3").c_str());  arg->erase(arg->find("c3"));}
  if (arg->find("c4")!=end)   {c[4]    = atof(arg->at("c4").c_str());  arg->erase(arg->find("c4"));}


  if (arg->size()!=0){
    if (rank==0){ cout << "*** Error: Unknown element in &profile_polynom" << endl; this->usage();}
    return "";
  }
  if ((label.size()<1)&&(rank==0)){
    cout << "*** Error: Label not defined in &profile_polynom" << endl; this->usage();
  }
  return label;
}


double ProfilePolynom::value(double z)
{
  double val=0;
  double zsave=1;
  for (int i=0;i<c.size();i++){
    val+=c[i]*zsave;
    zsave*=z;
  }
  return val;
}

void ProfilePolynom::usage(){
  cout << "List of keywords for PROFILE_POLYNOM" << endl;
  cout << "&profile_polynom" << endl;
  cout << " string label = <empty>" << endl;
  cout << " double c0 = 0" << endl;
  cout << " double c1 = 0" << endl;
  cout << " double c2 = 0" << endl;
  cout << " double c3 = 0" << endl;
  cout << " double c4 = 0" << endl;
  cout << "&end" << endl << endl;
  return;
}

//-----------------------

 string ProfileStep::init(int rank, map<string,string>*arg)
{

  string label="";
  c0=0;
  sstart=0;
  send=0;

  map<string,string>::iterator end=arg->end();

  if (arg->find("label")!=end)  {label = arg->at("label");  arg->erase(arg->find("label"));}
  if (arg->find("c0")!=end)     {c0    = atof(arg->at("c0").c_str());  arg->erase(arg->find("c0"));}
  if (arg->find("s_start")!=end){sstart= atof(arg->at("s_start").c_str());  arg->erase(arg->find("s_start"));}
  if (arg->find("s_end")!=end)  {send  = atof(arg->at("s_end").c_str());  arg->erase(arg->find("s_end"));}

  if (arg->size()!=0){
    if (rank==0){ cout << "*** Error: Unknown elements in &profile_step" << endl; this->usage();}
    return "";
  }
  if ((label.size()<1)&&(rank==0)){
    cout << "*** Error: Label not defined in &profile_step" << endl; this->usage();
  }
  return label;
}

double ProfileStep::value(double z)
{
  if ((z>=sstart) && (z <=send)){
    return c0;
  }
  return 0;
}

void ProfileStep::usage(){

  cout << "List of keywords for PROFILE_STEP" << endl;
  cout << "&profile_step" << endl;
  cout << " string label = <empty>" << endl;
  cout << " double c0 = 0" << endl;
  cout << " double s_start = 0" << endl;
  cout << " double s_end = 0" << endl;
  cout << "&end" << endl << endl;
  return;
}

//-----------------------------------

string ProfileGauss::init(int rank, map<string,string>*arg)
{
  string label="";
  c0=0;
  s0=0;
  sig=1;

  map<string,string>::iterator end=arg->end();

  if (arg->find("label")!=end){label = arg->at("label");  arg->erase(arg->find("label"));}
  if (arg->find("c0")!=end)   {c0    = atof(arg->at("c0").c_str());  arg->erase(arg->find("c0"));}
  if (arg->find("s0")!=end)   {s0    = atof(arg->at("s0").c_str());  arg->erase(arg->find("s0"));}
  if (arg->find("sig")!=end)  {sig   = atof(arg->at("sig").c_str());  arg->erase(arg->find("sig"));}

  if (arg->size()!=0){
    if (rank==0){ cout << "*** Error: Unknown elements in &profile_gauss" << endl; this->usage();}
    return "";
  }
  if ((label.size()<1)&&(rank==0)){
    cout << "*** Error: Label not defined in &profile_gauss" << endl; this->usage();
  }
  return label;
}

double ProfileGauss::value(double z)
{
  return c0*exp(-0.5*(z-s0)*(z-s0)/sig/sig);
}

void ProfileGauss::usage(){ 
  cout << "List of keywords for PROFILE_GAUSS" << endl;
  cout << "&profile_gauss" << endl;
  cout << " string label = <empty>" << endl;
  cout << " double c0 = 0" << endl;
  cout << " double s0 = 0" << endl;
  cout << " double sig= 1" << endl;
  cout << "&end" << endl << endl;
  return;
  return;
}


//-----------------------------------

string ProfileFile::init(int rank, map<string,string>*arg)
{
  string label="";
  xdataset="";
  ydataset="";
  isTime=false;
  revert=false;

  map<string,string>::iterator end=arg->end();

  if (arg->find("label")!=end){label = arg->at("label");  arg->erase(arg->find("label"));}
  if (arg->find("xdata")!=end)   {xdataset = arg->at("xdata");  arg->erase(arg->find("xdata"));}
  if (arg->find("ydata")!=end)   {ydataset = arg->at("ydata");  arg->erase(arg->find("ydata"));}
  if (arg->find("isTime")!=end)  {isTime = atob(arg->at("isTime").c_str()); arg->erase(arg->find("isTime"));}
  if (arg->find("reverse")!=end) {revert = atob(arg->at("reverse").c_str()); arg->erase(arg->find("reverse"));}
  

  if (arg->size()!=0){
    if (rank==0){ cout << "*** Error: Unknown elements in &profile_file" << endl; this->usage();}
    return "";
  }
  if ((label.size()<1)&&(rank==0)){
    cout << "*** Error: Label not defined in &profile_file" << endl; this->usage();
  }

  int ndata=-1;
  bool success;

  success=this->simpleReadDouble1D(xdataset,&xdat);
  if (!success){
    if (rank==0){
      cout << "*** Error: Cannot read the HDF5 dataset: " << xdataset << endl;
    }
    return "";
  }

  success=this->simpleReadDouble1D(ydataset,&ydat);
  if (!success){
    if (rank==0){
      cout << "*** Error: Cannot read the HDF5 dataset: " << ydataset << endl;
    }
    return "";
  }


  if (isTime){ 
    for (int i=0; i<xdat.size();i++){
      xdat[i]*=3e8;         // scale time variable to space varial by multiplying the speed of light
    }  
  }
  
  if (revert){
    double xmin=xdat[0];
    double xmax=xdat[xdat.size()-1];
    reverse(xdat.begin(),xdat.end());
    reverse(ydat.begin(),ydat.end());
    for (int i=0;i<xdat.size();i++){
      xdat[i]=-xdat[i]+xmin+xmax;    // get the correct time window
    }
  }
  return label;
}

double ProfileFile::value(double z)
{
  if (z<xdat[0]){ return ydat[0]; }
  if (z>xdat[xdat.size()-1]){ return ydat[xdat.size()-1]; }
  int idx=0;
  while(z>=xdat[idx]){
    idx++;
  }
  idx--;
  double wei=(z-xdat[idx])/(xdat[idx+1]-xdat[idx]);
  double val=ydat[idx]*(1-wei)+wei*ydat[idx+1];
  return val;
}

void ProfileFile::usage(){ 
  cout << "List of keywords for PROFILE_FILE" << endl;
  cout << "&profile_file" << endl;
  cout << " string label = <empty>" << endl;
  cout << " string xdata = <empty>" << endl;
  cout << " string ydata = <empty>" << endl;
  cout << " bool isTime  = false" << endl;
  cout << " bool reverse = false" << endl;
  cout << "&end" << endl << endl;
  return;
  return;
}
