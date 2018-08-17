#include "Parser.h"

Parser::Parser()
{
}

Parser::~Parser()
{
}

bool Parser::open(string file, int inrank, bool streaming)
{

  string instring;
  ostringstream os;
  fstream fin;

  rank=inrank;
  if (!streaming){
      fin.open(file.c_str(),ios_base::in);
      if (!fin){
         if (rank==0) {cout << "*** Error: Cannot open main input file: " << file << endl;}
         return false;
      }
      while(getline(fin,instring,'\n')){
	os << instring <<endl;
      }
      fin.close();
      input.str(os.str());
  } else {
    input.str(file);   // in the case that the input file is streamed
  }
  return true;

}



bool Parser::parse(string *element, map<string,string> *argument)
{

  string instring, list;
  bool accumulate=false;
  
  
  while(getline(input,instring)){    // read line
    this->trim(instring);

    // check for empty lines or comments
    if ((!instring.compare(0,1,"#") || instring.length() < 1)){ continue; } // skip comment and empty rows

    // check for terminating element and then return
    if ((instring.compare("&end")==0)|| (instring.compare("&END")==0) || (instring.compare("&End")==0)){
      if (!accumulate){
        if (rank==0){cout << "*** Error: Termination string outside element definition in input file" << endl;}    
	return false;
      }
      return this->fillMap(&list,argument);
    }

    // check whether element starts
    if (!instring.compare(0,1,"&")){
      if (accumulate){
        if (rank==0) {cout << "*** Error: Nested elements in main input file" << endl;}
	return false;
      }
      accumulate=true;
      for (int i=0; i<instring.size();i++){
	instring[i]=tolower(instring[i]);
      }
      *element=instring;   
      continue;

    }

    // add content if in element
    if (accumulate){
       list.append(instring);  // add all content into one string
       list.append(";");
    }
  }
  
  return false;

}

bool Parser::fillMap(string *list,map<string,string> *map){

  // splits the string into a map with pairs : key and value.
  // keys are converted to lower case
  
  map->clear();
  string key,val;
 
  size_t pos,tpos;
  while ((pos=list->find_first_of(";")) !=string::npos){  // split into individual lines
     val=list->substr(0,pos);
     list->erase(0,pos+1);
     this->trim(key);
     if (val.length()>0){
       tpos=val.find_first_of("=");
       if (tpos ==string::npos){
         if (rank==0){ cout << "*** Error: Invalid format " << val << " in input file" << endl;}
         return false;
       } else{
         key=val.substr(0,tpos);
         val.erase(0,tpos+1);
         this->trim(key);
         this->trim(val);
         map->insert(pair<string,string>(key,val));
       }
     }
  }
  
  return true;
}
