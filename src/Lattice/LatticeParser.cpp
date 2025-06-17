#include <algorithm> // for std::find
#include <utility>
#include "LatticeParser.h"
#include "SeriesManager.h"
#include "SeriesParser.h"

LatticeParser::LatticeParser()
{
}


LatticeParser::~LatticeParser()
{
}

void LatticeParser::addSequence(const std::string& label_in, std::string argument_in, int rank, SeriesManager *sm)
{
    vector <string> par;
    string fld,val;
    std::map <string,string> maparg;
    maparg.insert({"label",label_in});
    std::string seqtype;

    this->chop(std::move(argument_in),&par);
    for (int i=0;i<par.size();i++) {
        if (par[i].empty()) {  // in the case of C: COR = {};  one is getting a zero length string.
            continue;
        }
        size_t pos = par[i].find_first_of("=");
        if (pos == string::npos) {
            if (rank == 0) {
                cout << "*** Warning: Ignoring invalid format: " << par[i] << " for element " << label_in << endl;
            }
            continue;
        }
        fld = par[i].substr(0, pos);
        val = par[i].erase(0, pos + 1);
        this->trim(fld);
        this->trim(val);
        if (fld.compare("type") == 0){
            seqtype = "&sequence_"+val;
        } else {
            maparg.insert({fld, val});
        }
    }
    auto *seriesparser = new SeriesParser;
    seriesparser->init(rank, &maparg, seqtype, sm);
    delete seriesparser;
}

bool LatticeParser::parse(string file, string line, int rank, vector<Element *> &lat, SeriesManager *sm)
{

  istringstream input;
  string instring;
  ostringstream os;
  
    ifstream fin(file.c_str(),ios_base::in);

    if (!fin){
      if (rank==0) {cout << "*** Error: Cannot open magnetic lattice file: " << file << endl;}
      return false;
    }
    while(getline(fin,instring,'\n')){
	os << instring <<endl;
    }
    fin.close();
    input.str(os.str());
 

  //------------------------------------------------------
  // step one - coarse parsing of the input deck

  string comstring="";
  vector<string> content;

  while(getline(input,instring)){    // read line
    //    cout << instring<<endl;
    this->trim(instring);
    if ((!instring.compare(0,1,"#") || instring.length() < 1)){ continue; } // skip comment and empty rows
    comstring.append(" ");
    comstring.append(instring);  // add all content into one string
  }


  for (int i=0; i<comstring.size();i++){ // convert to lower case
    comstring[i]=tolower(comstring[i]);
  }

  size_t pos;
  while ((pos=comstring.find_first_of(";")) !=string::npos){  // split into individual lines
     instring=comstring.substr(0,pos);
     this->trim(instring);
     content.push_back(instring);
     comstring.erase(0,pos+1);
  }
  
  this->trim(comstring);
  if ((comstring.length()>0) && (rank==0)){
    cout << "*** Warning: Ignoring incomplete last line in magnetic lattice" << endl;
  } 

  // -----------------------------------------------------------------------
  // step two - parse each individual line of the lattice according to the format
  // label: type =(content);
  // e.g.   "QF1: Quadrupole = {L=0.2, k1=0.8}; 

  label.clear();
  type.clear();
  argument.clear();

  string inlabel, intype, inargument;
  bool error=false;
 
  for(int i=0;i<content.size();i++){
    
    if ((pos=content[i].find_first_of(':'))==string::npos){
      if (rank==0){ cout<< "*** Error: Invalid Format in lattice file: " << content[i] << endl;}
      error=true;
      continue;
    }
    inlabel=content[i].substr(0,pos);
    this->trim(inlabel);
    content[i].erase(0,pos+1);
    if ((pos=content[i].find_first_of('='))==string::npos){
      if (rank==0){ cout<< "*** Error: Invalid Format in lattice file: " << content[i] << endl;}
      error=true;
      continue;
    }
    intype=content[i].substr(0,pos);
    this->trim(intype);
    content[i].erase(0,pos+1);
    
    if ((pos=content[i].find_first_of('{'))==string::npos){
      if (rank==0){ cout<< "*** Error: Invalid Format in lattice file: " << content[i] << endl;}
      error=true;
      continue;
    }
    content[i].erase(0,pos+1);
    if ((pos=content[i].find_first_of('}'))==string::npos){
      if (rank==0){ cout<< "*** Error: Invalid Format in lattice file: " << content[i] << endl;}
      error=true;
      continue;
    }

    content[i].erase(pos,content[i].size());
    inargument=content[i];
    this->trim(inargument);
    type.push_back(intype.substr(0,4));
    argument.push_back(inargument);

    if(std::find(label.begin(), label.end(), inlabel)!=label.end()) {
      // name of this lattice element was already defined -> error
      if(0==rank) {
        cout << "*** Error in lattice: Label " << inlabel << " redefined" << endl;
        error=true;
        continue;
      }
    }
    label.push_back(inlabel);
  }

  if (error){ return false; }

  // -------------------------------------------------------------
  // step 2.5 - extract seuqence definition and add them to seriesmanager.

  for(std::size_t i = 0; i < type.size(); i++){
      if (type.at(i).compare("sequ") ==0) {
          this->addSequence(label.at(i),argument.at(i),rank,sm);
      }
  }



    // -------------------------------------------------------------
  // step 3 - resolving all references

  
  int recursion = 20;
  for (int i=0;i<label.size();i++){
    if (type[i].compare("line")!=0){
      error=this->resolve(i,recursion-1,rank);
      if (error==false){ return false; }
    }   
  } 


  // --------------------------------------------------------------
  // step 4 - unrolling the line
  
  for (int i=0;i<line.size();i++){
    line[i]=tolower(line[i]);      
  }

  int idx=this->findIndex(&label,line);

  sequence.clear();
  zref.clear();
  zoff.clear();
  refele=-1;
  
  if ((idx>-1) && (type[idx].compare("line")==0)) {
    error=this->unroll(idx, recursion-1, rank);
  } else {
      if (rank==0) {cout << "*** Error: Lattice file does not contain the beamline: " << line << endl;}
      return false;  
  }


  if (error==false){ return false; }

  //---------------------------------------------------------------------------------
  // step 5 - initiate

  double z;
  
  for (int i=0;i<sequence.size();i++){
    if (zref[i]<0){
      z=0+zoff[i];
    } else {
      z=lat[zref[i]]->z+lat[zref[i]]->l+zoff[i];
    } 
    idx=this->findIndex(&label,sequence[i]);
    // //    cout << "Element: " << sequence[i] << " Type: " << type[idx] << " Pos:" << z << " Offset: " << zoff[i] << endl;
    // cout << "Element " << i << " Type: " << type[idx] << endl;
    if (type[idx].compare("quad")==0){ error=false; lat.push_back(this->parseQuad(idx,rank,z));}
    if (type[idx].compare("undu")==0){ error=false; lat.push_back(this->parseID(idx,rank,z,sm)); }
    if (type[idx].compare("drif")==0){ error=false; lat.push_back(this->parseDrift(idx,rank,z));}
    if (type[idx].compare("corr")==0){ error=false; lat.push_back(this->parseCorrector(idx,rank,z));}
    if (type[idx].compare("chic")==0){ error=false; lat.push_back(this->parseChicane(idx,rank,z,sm)); }
    if (type[idx].compare("mark")==0){ error=false; lat.push_back(this->parseMarker(idx,rank,z)); }
    if (type[idx].compare("phas")==0){ error=false; lat.push_back(this->parsePhaseshifter(idx,rank,z,sm)); }
    if (error) { 
      if (rank==0) {cout << "*** Error: Unknown element type " << type[idx] << " in lattice" << endl;}
      return false;
    }
  }
  return true;
}


ID *LatticeParser::parseID(int idx,int rank, double zin, SeriesManager *sm)
{
  ID *ele=new ID;
  
  ele->type="Undulator";
  ele->z=zin;
  ele->zoff=0;
  ele->aw=0;
  ele->lambdau=0;
  ele->nwig=0;
  ele->kx=0;
  ele->ky=1;
  ele->ax=0;
  ele->ay=0;
  ele->gradx=0;
  ele->grady=0;
  ele->helical=false;
  // current unused
  ele->paw=0;
  ele->pkx=0;
  ele->pky=1;
  ele->pdadx=0;
  ele->pdady=0;
  ele->phase=0;


  vector <string> par;
  string fld,val;
  bool haskx = false;

  this->chop(argument[idx],&par);
  for (int i=0;i<par.size();i++){
      if (par[i].length() < 1){  // in the case of C: COR = {};  one is getting a zero length string.
          continue;
      }
    size_t pos=par[i].find_first_of("=");
    if (pos==string::npos){
      if (rank==0){cout << "*** Warning: Ignoring invalid format: " << par[i] << " for element " << label[idx]<< endl;}
      continue;  
    }
    fld=par[i].substr(0,pos);
    val=par[i].erase(0,pos+1);
    this->trim(fld);
    this->trim(val);
    bool found=false;
    if (fld.compare("l")==0)    { ele->l=atof(val.c_str());  found=true; };
    if (fld.compare("lambdau")==0){ ele->lambdau=atof(val.c_str()); found=true; };
#if 0
    if (fld.compare("aw")==0) {
      /* determine aw value of this undulator */
      double this_aw=-1.;
      string this_aw_refname;
      // if the value string begins with "@", we assume it is a reference to a sequence
      this->reference(val, &this_aw, &this_aw_refname);
      if(!this_aw_refname.empty()) {
        // we got a reference to a sequence, query it
        this_aw=sm->getElement(this_aw_refname);
        if(0==rank) {
          cout << "*** Ref to series " << this_aw_refname << ", result is: " << this_aw << " ***" << endl;
        }
      }
      ele->aw=this_aw; // atof(val.c_str());
      found=true;
    }
#endif
    if(extractParameterValue(fld,val,sm,rank, "aw", &ele->aw)) {found=true;}
    if (fld.compare("aw_perp")==0)  { ele->paw=atof(val.c_str()); found=true; };
    if (fld.compare("nwig")==0) { ele->nwig=atof(val.c_str()); found=true; };
    if (fld.compare("kx")==0)   { ele->kx=atof(val.c_str()); found=true; haskx = true;};
    if (fld.compare("ky")==0)   { ele->ky=atof(val.c_str()); found=true; haskx = true;};
    if (fld.compare("kx_perp")==0)  { ele->pkx=atof(val.c_str()); found=true; };
    if (fld.compare("ky_perp")==0)  { ele->pky=atof(val.c_str()); found=true; };
    if (fld.compare("ax")==0)   { ele->ax=atof(val.c_str()); found=true; };
    if (fld.compare("ay")==0)   { ele->ay=atof(val.c_str()); found=true; };
    if (fld.compare("gradx")==0) { ele->gradx=atof(val.c_str()); found=true; };
    if (fld.compare("grady")==0) { ele->grady=atof(val.c_str()); found=true; };
    if (fld.compare("gradx_perp")==0){ ele->pdadx=atof(val.c_str()); found=true; };
    if (fld.compare("grady_perp")==0){ ele->pdady=atof(val.c_str()); found=true; };
    if (fld.compare("phase_perp")==0){ ele->phase=atof(val.c_str()); found=true; };
    if (fld.compare("helical")==0)   { ele->helical=atob(val.c_str()); found=true; };
    if (found==false){
      if (rank==0){cout << "*** Warning: Ignoring unknown parameter: " << fld << " for element " << label[idx]<< endl;}
    }
  }
  ele->l=ele->lambdau*ele->nwig;
  if (!haskx){
      if (ele->helical){
          ele->kx = 0.5;
          ele->ky = 0.5;
      } else {
          ele->kx =0;
          ele->ky = 1;
      }
  }
  return ele;
}


Corrector *LatticeParser::parseCorrector(int idx,int rank, double zin)
{
  Corrector *ele=new Corrector;

  ele->type="Corrector";
  ele->z=zin;
  ele->zoff=0;
  ele->l=0;
  ele->cx=0;
  ele->cy=0;
  vector <string> par;
  string fld,val;

  this->chop(argument[idx],&par);
  for (int i=0;i<par.size();i++){
      if (par[i].length() < 1){  // in the case of C: COR = {};  one is getting a zero length string.
          continue;
      }
    size_t pos=par[i].find_first_of("=");
    if (par[i].length() < 1){  // in the case of C: COR = {};  one is getting a zero length string.
        continue;
    }

    fld=par[i].substr(0,pos);
    val=par[i].erase(0,pos+1);
    this->trim(fld);
    bool found=false;
    if (fld.compare("l")==0)    { ele->l=atof(val.c_str());  found=true; };
    if (fld.compare("cx")==0)   { ele->cx=atof(val.c_str()); found=true; };
    if (fld.compare("cy")==0)   { ele->cy=atof(val.c_str()); found=true; };
    if (found==false){
      if (rank==0){cout << "*** Warning: Ignoring unknown parameter: " << fld << " for element " << label[idx]<< endl;}
    }
  }
  return ele;
}


Chicane *LatticeParser::parseChicane(int idx,int rank, double zin, SeriesManager *sm)
{
  Chicane *ele=new Chicane;

  ele->type="Chicane";
  ele->z=zin;
  ele->zoff=0;
  ele->l=0;
  ele->delay=0;
  ele->lb=0;
  ele->ld=0;
  vector <string> par;
  string fld,val;

  this->chop(argument[idx],&par);
  for (int i=0;i<par.size();i++){
      if (par[i].length() < 1){  // in the case of C: COR = {};  one is getting a zero length string.
          continue;
      }
    size_t pos=par[i].find_first_of("=");
    if (pos==string::npos){
      if (rank==0){cout << "*** Warning: Ignoring invalid format: " << par[i] << " for element " << label[idx]<< endl;}
      continue;  
    }
    fld=par[i].substr(0,pos);
    val=par[i].erase(0,pos+1);
    this->trim(fld);
    bool found=false;
    if (fld.compare("l")==0)    { ele->l=atof(val.c_str());  found=true; };
    if (extractParameterValue(fld,val,sm,rank, "delay", &ele->delay)){ found=true; };
    if (fld.compare("lb")==0)   { ele->lb=atof(val.c_str()); found=true; };
    if (fld.compare("ld")==0)   { ele->ld=atof(val.c_str()); found=true; };
    if (found==false){
      if (rank==0){cout << "*** Warning: Ignoring unknown parameter: " << fld << " for element " << label[idx]<< endl;}
    }
  }
  return ele;
}

Marker *LatticeParser::parseMarker(int idx,int rank, double zin)
{
  Marker *ele=new Marker;

  ele->type="Marker";
  ele->z=zin;
  ele->zoff=0;
  ele->l=0;
  ele->action=0;
  vector <string> par;
  string fld,val;

  this->chop(argument[idx],&par);
  for (int i=0;i<par.size();i++){
      if (par[i].length() < 1){  // in the case of C: COR = {};  one is getting a zero length string.
          continue;
      }
    size_t pos=par[i].find_first_of("=");
    if (pos==string::npos){
      if (rank==0){cout << "*** Warning: Ignoring invalid format: " << par[i] << " for element " << label[idx]<< endl;}
      continue;  
    }
    fld=par[i].substr(0,pos);
    val=par[i].erase(0,pos+1);
    this->trim(fld);
    bool found=false;
    int tag=atoi(val.c_str());
    if (tag !=0) { tag=1;}
    if (fld.compare("dumpfield")==0){ele->action|=1*tag;found=true;}
    if (fld.compare("dumpbeam")==0){ele->action|=2*tag;found=true;}
    if (fld.compare("sort")==0){ele->action|=4*tag;found=true;}
    if (fld.compare("stop")==0){ele->action|=8*tag;found=true;}
    if (found==false){
      if (rank==0){cout << "*** Warning: Ignoring unknown parameter: " << fld << " for element " << label[idx]<< endl;}
    }
  }
  return ele;
}



Drift *LatticeParser::parseDrift(int idx,int rank, double zin)
{
  Drift *ele=new Drift;

  ele->type="Drift";
  ele->z=zin;
  ele->zoff=0;
  ele->l=0;
  vector <string> par;
  string fld,val;

  this->chop(argument[idx],&par);
  for (int i=0;i<par.size();i++){
      if (par[i].length() < 1){  // in the case of C: COR = {};  one is getting a zero length string.
          continue;
      }
    size_t pos=par[i].find_first_of("=");
    if (pos==string::npos){
      if (rank==0){cout << "*** Warning: Ignoring invalid format: " << par[i] << " for element " << label[idx]<< endl;}
      continue;  
    }
    fld=par[i].substr(0,pos);
    val=par[i].erase(0,pos+1);
    this->trim(fld);
    bool found=false;
    if (fld.compare("l")==0) { ele->l=atof(val.c_str());  found=true; };
    if (found==false){
      if (rank==0){cout << "*** Warning: Ignoring unknown parameter: " << fld << " for element " << label[idx]<< endl;}
    }
  }
  return ele;
}


Quadrupole *LatticeParser::parseQuad(int idx,int rank, double zin)
{
  Quadrupole *ele=new Quadrupole;

  ele->type="Quadrupole";
  ele->z=zin;
  ele->zoff=0;
  ele->l=0;
  ele->k1=0;
  ele->dx=0;
  ele->dy=0;
  vector <string> par;
  string fld,val;

  this->chop(argument[idx],&par);
  for (int i=0;i<par.size();i++){
      if (par[i].length() < 1){  // in the case of C: COR = {};  one is getting a zero length string.
          continue;
      }
    size_t pos=par[i].find_first_of("=");
    if (pos==string::npos){
      if (rank==0){cout << "*** Warning: Ignoring invalid format: " << par[i] << " for element " << label[idx]<< endl;}
      continue;  
    }
    fld=par[i].substr(0,pos);
    val=par[i].erase(0,pos+1);
    this->trim(fld);
    bool found=false;
    if (fld.compare("l")==0) { ele->l=atof(val.c_str());  found=true; };
    if (fld.compare("k1")==0){ ele->k1=atof(val.c_str()); found=true; };
    if (fld.compare("dx")==0){ ele->dx=atof(val.c_str()); found=true; };
    if (fld.compare("dy")==0){ ele->dy=atof(val.c_str()); found=true; };
    if (found==false){
      if (rank==0){cout << "*** Warning: Ignoring unknown parameter: " << fld << " for element " << label[idx]<< endl;}
    }
  }
  return ele;
}

Phaseshifter *LatticeParser::parsePhaseshifter(int idx,int rank, double zin, SeriesManager *sm)
{
  Phaseshifter *ele=new Phaseshifter;

  ele->type="Phaseshifter";
  ele->z=zin;
  ele->zoff=0;
  ele->l=0;
  ele->phi=0;

  vector <string> par;
  string fld,val;

  this->chop(argument[idx],&par);
  for (int i=0;i<par.size();i++){
      if (par[i].length() < 1){  // in the case of C: COR = {};  one is getting a zero length string.
          continue;
      }
    size_t pos=par[i].find_first_of("=");
    if (pos==string::npos){
      if (rank==0){cout << "*** Warning: Ignoring invalid format: " << par[i] << " for element " << label[idx]<< endl;}
      continue;
    }
    fld=par[i].substr(0,pos);
    val=par[i].erase(0,pos+1);
    this->trim(fld);
    bool found=false;
    if (fld.compare("l")==0) { ele->l=atof(val.c_str());  found=true; };
    // if (fld.compare("phi")==0){ele->phi=atof(val.c_str());found=true; };
    if (extractParameterValue(fld,val,sm,rank, "phi", &ele->phi)) {found=true;}
    if (found==false){
      if (rank==0){cout << "*** Warning: Ignoring unknown parameter: " << fld << " for element " << label[idx]<< endl;}
    }
  }
  return ele;
}

bool LatticeParser::resolve(int idx, int recursion, int rank)
{
  if (recursion < 0){
    if (rank==0){  cout << "*** Error: Too many nested references" << endl;}
   return false;
  }
  
  size_t pos=argument[idx].find("ref");
  if (pos==string::npos){ return true;}
  vector<string> arg1;
  this->chop(argument[idx], &arg1);

  string ref;
  for (int i=0; i< arg1.size(); i++){
    pos=arg1[i].find("ref");
    if (pos !=string::npos){
      pos=arg1[i].find_first_of("=");
      ref=arg1[1].erase(0,pos+1);
      arg1.erase(arg1.begin()+i);
      break;
    }
  } 
  this->trim(ref); 
  int iref=findIndex(&label,ref);
  if (iref<0){
      if (rank==0) {cout << "*** Error: Unresolved reference: " << ref <<" for element: " << label[idx] << endl;}
      return false;
  }  

  if (this->resolve(iref,recursion-1,rank)==false) {return false;}

  argument[idx]=argument[iref];  
  for (int i=0; i< arg1.size(); i++){
    argument[idx].append(",");
    argument[idx].append(arg1[i]);
  }
 
  return true;



}

bool LatticeParser::unroll(int idx, int recursion,int rank){

 
  if (recursion < 0){
    if (rank==0){  cout << "*** Error: Too many nested elements in selected beamline" << endl;}
   return false;
  }
  vector<string> line;
  this->chop(argument[idx],&line);

  int refelesave=refele;  // safe the index to the last element. Is initialized to -1 at start

  for (int i=0; i<line.size();i++){
    int count = checkMultiplier(&line[i]);  // check for multiplier
    double off;
    bool resetpos = checkResetPosition(&line[i], &off);
    int ix=this->findIndex(&label,line[i]);

    if (ix <0){
      if (rank==0) {cout << "*** Error: Undefined element in beamline: " << line[i] << endl;}
      return false;
    }
    if (type[ix].compare("line")==0){
      for (int j=0; j < count ; j++){
	bool error=this->unroll(ix,recursion-1,rank);
	if (error==false) { return error; }
      }
    } else {
      for (int j=0; j < count ; j++){
        sequence.push_back(line[i]); 
	if (resetpos) {
	  zref.push_back(refelesave);
	  zoff.push_back(off);
	} else {
	  zref.push_back(refele);
	  zoff.push_back(0.);
	}
        refele++;	
      }
    }
  }

  line.clear();
  return true;

  /*

  int refelesave=refele;

  for (int i=0; i<line.size();i++){

    int count = checkMultiplier(&line[i]);  // check if a line has a multiplier
    bool resetpos = checkResetPosition(&line[i]);
    int ix=this->findIndex(&label,line[i]);
    if (ix <0){
      if (rank==0) {cout << "*** Error: Undefined element in beamline: " << line[i] << endl;}
      return false;
    }
    for (int j=0; j < count ; j++){
      if (type[ix].compare("line")==0){
        bool error=this->unroll(ix,recursion-1,rank);
        if (error==false) { return error; }
        if (resetpos) { refele=refelesave; }
      } else {        
        sequence.push_back(line[i]); 
        zref.push_back(refele);
        refele++;
        if (resetpos) { refele--; }
      }
    }
  }

  line.clear();
  return true;
*/

}


int LatticeParser::checkMultiplier(string *element){

  size_t pos=element->find_first_of("*");
  if (pos==string::npos){
    return 1;
  } else {
    string num=element->substr(0,pos);
    element->erase(0,pos+1);
    this->trim(num);
    this->trim(*element);
    return atoi(num.c_str());
  }
}

bool LatticeParser::checkResetPosition(string *element,double *off){
   
  *off=0;
  size_t pos=element->find_first_of("@");
  if (pos==string::npos){
    return false;
  }
  string num=element->substr(pos+1);
  element->erase(pos);
  this->trim(num);
  this->trim(*element);
  *off=atof(num.c_str());
  return true; 
} 

int LatticeParser::findIndex(vector<string> *list, string element){

  for (int i=list->size()-1;i>-1;i--){
    if (list->at(i).compare(element)==0) { return i; }
  }
    return -1;
}

// The data value 'pv' points to is only changed in case of success
bool LatticeParser::extractParameterValue(string fld, string val, SeriesManager *sm, int rank, string parameterName, double *pv)
{
   if (fld.compare(parameterName)!=0)
     return(false); // this is not the parameter we are looking for

   /* determine parameter value */
   double q=0.;
   string refname;
   // if the value string begins with "@", we assume it is a reference to a sequence
   this->reference(val, &q, &refname);
   if(!refname.empty()) {
     // we got a reference to a sequence, query it
     q=sm->getElement(refname);
     if(0==rank) {
       cout << "*** Ref to series " << refname << ", result is: " << parameterName<<"="<< q << " ***" << endl;
     }
   }
   *pv = q;
   // ele->aw=this_aw; // atof(val.c_str());
   return(true);
}
