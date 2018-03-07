#include "LoadField.h"

LoadField::LoadField()
{
  lambda=0;
  lambdaref="";
  power=0;
  powerref="";
  phase=0;
  phaseref="";
  z0=0;
  z0ref="";
  w0=100e-6;
  w0ref="";
  dgrid=1e-3;
  ngrid=151;
  harm=1;
  nx=0;
  ny=0;
  xcen=0;
  ycen=0;
  xangle=0;
  yangle=0;
  add=false;
}

LoadField::~LoadField(){}



void LoadField::usage(){

  cout << "List of keywords for FIELD" << endl;
  cout << "&field" << endl;
  cout << " double lambda = lambdaref" << endl;
  cout << " double power = 0 / reference" << endl;
  cout << " double phase = 0 / reference" << endl;
  cout << " double waist_pos = 0 / reference" << endl;
  cout << " double waist_size = 100e-6 / reference" << endl;
  cout << " double xcenter = 0" << endl;
  cout << " double ycenter = 0" << endl;
  cout << " double xangle = 0" << endl;
  cout << " double yangle = 0" << endl;
  cout << " double dgrid = 1e-3 / from existing field" << endl;
  cout << " int ngrid = 151 / from existing field" << endl;
  cout << " int harm = 1" << endl;
  cout << " int nx = 0" << endl;
  cout << " int ny = 0" << endl;
  cout << " bool accumulate = false " << endl;
  cout << "&end" << endl << endl;
  return;
}

bool LoadField::init(int rank, int size, map<string,string> *arg, vector<Field *> *fieldin,  Setup *setup, Time *time, Profile *prof)
{

  bool dotime=time->isTime();                  // check for time simulation
  double sample=static_cast<double>(time->getSampleRate());         // check slice length
  lambda=setup->getReferenceLength();

  map<string,string>::iterator end=arg->end();

  if (arg->find("lambda")!=end){this->reference(arg->at("lambda"),&lambda,&lambdaref); arg->erase(arg->find("lambda"));}
  if (arg->find("power")!=end) {this->reference(arg->at("power"),&power,&powerref); arg->erase(arg->find("power"));}
  if (arg->find("phase")!=end) {this->reference(arg->at("phase"),&phase,&phaseref); arg->erase(arg->find("phase"));}
  if (arg->find("waist_pos")!=end) {this->reference(arg->at("waist_pos"),&z0,&z0ref); arg->erase(arg->find("waist_pos"));}
  if (arg->find("waist_size")!=end){this->reference(arg->at("waist_size"),&w0,&w0ref); arg->erase(arg->find("waist_size"));}  
  if (arg->find("xcenter")!=end) {xcen = atof(arg->at("xcenter").c_str()); arg->erase(arg->find("xcenter"));}
  if (arg->find("ycenter")!=end) {ycen = atof(arg->at("ycenter").c_str()); arg->erase(arg->find("ycenter"));}
  if (arg->find("xangle")!=end)  {xangle = atof(arg->at("xangle").c_str()); arg->erase(arg->find("xangle"));}
  if (arg->find("yangle")!=end)  {yangle = atof(arg->at("yangle").c_str()); arg->erase(arg->find("yangle"));}
  if (arg->find("dgrid")!=end) {dgrid = atof(arg->at("dgrid").c_str()); arg->erase(arg->find("dgrid"));}
  if (arg->find("ngrid")!=end) {ngrid = atoi(arg->at("ngrid").c_str());  arg->erase(arg->find("ngrid"));}
  if (arg->find("harm")!=end)  {harm  = atoi(arg->at("harm").c_str());  arg->erase(arg->find("harm"));}
  if (arg->find("nx")!=end)    {nx  = atoi(arg->at("nx").c_str());  arg->erase(arg->find("nx"));}
  if (arg->find("ny")!=end)    {ny  = atoi(arg->at("ny").c_str());  arg->erase(arg->find("ny"));}
  if (arg->find("accumulate")!=end) {add  = atob(arg->at("accumulate"));  arg->erase(arg->find("accumulate"));}

  if (arg->size()!=0){
    if (rank==0){ cout << "*** Error: Unknown elements in &field" << endl; this->usage();}
    return false;
  }

  // check for existing field record;
  int idx=-1;
  Field *field;
  for (int i=0; i<fieldin->size();i++){
    if (fieldin->at(i)->harm==harm){
      field=fieldin->at(i);
      idx=i;
      if (add){
         dgrid=field->gridmax;    // take grid size from existing field
         ngrid=field->ngrid; 
      }
    }
  }

  if (idx<0){
    field=new Field;
    add=false;
  }

  // field points now to the record which is filled


  if (rank==0){ 
    if (add) { cout << "Adding "; } else  { cout << "Generating "; }
    cout << "input radiation field for HARM = " << harm <<  " ..." << endl; 
  }

  vector<double> s;
  int nslice=time->getPosition(&s);

  field->init(time->getNodeNSlice(),ngrid,dgrid,lambda,sample*lambda,s[0],harm);
  
  if (idx<0){
    fieldin->push_back(field);
    idx=fieldin->size()-1;
  }

  complex< double >  *fieldslice = new complex<double> [ngrid*ngrid];
  FieldSlice slice;
  GaussHermite gh;

  for (int j=0; j<time->getNodeNSlice(); j++){
    int i=j+time->getNodeOffset();
    slice.lambda=prof->value(s[i],lambda,lambdaref);
    slice.power=prof->value(s[i],power,powerref);
    slice.phase=prof->value(s[i],phase,phaseref);
    slice.z0=prof->value(s[i],z0,z0ref);
    slice.w0=prof->value(s[i],w0,w0ref);
    slice.xcen=xcen;
    slice.ycen=ycen;
    slice.xangle=xangle;
    slice.yangle=yangle;
    slice.nx=nx;
    slice.ny=ny;
    slice.harm=harm;
    gh.loadGauss(fieldslice,&slice,dgrid,ngrid);
    if (add){
      for (int k=0; k<ngrid*ngrid;k++){
        fieldin->at(idx)->field[j].at(k)+=fieldslice[k];
      } 
    } else {
      for (int k=0; k<ngrid*ngrid;k++){
        fieldin->at(idx)->field[j].at(k)=fieldslice[k]; 
      }      
    }

  }


  delete [] fieldslice;  

  // if (idx<0){fieldin->push_back(field);}

  

  return true;

}
