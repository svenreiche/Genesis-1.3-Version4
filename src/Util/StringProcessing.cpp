/*
 *  StringProcessing.cpp
 *  Genesis
 *
 *  Created by Sven Reiche on 9/10/10.
 *  Copyright 2010 Paul Scherrer Institut. All rights reserved.
 *
 */

#include "StringProcessing.h"

StringProcessing::StringProcessing(){}
StringProcessing::~StringProcessing(){}

//---------------------------------------------------------------------------
// some general string utility functions

void StringProcessing::chop(string str, vector<string> *list){
	
	size_t pos;
	
	list->clear();
	
	
	while((pos=str.find_first_of(","))!=string::npos){
	    list->push_back(str.substr(0,pos));
	    str.erase(0,pos+1);
	}
        list->push_back(str);
	
	for (int i=0; i<list->size();i++){
	  this->trim(list->at(i));
        }

	return;
}

void StringProcessing::trim(string &str){
    size_t startpos = str.find_first_not_of(" \t");
	size_t endpos = str.find_last_not_of(" \t");
	if(( string::npos == startpos ) || ( string::npos == endpos)){
		str="";
	} else {
		str = str.substr( startpos, endpos-startpos+1 );
	}
	return;
}


bool StringProcessing::atob(string in){
    
	bool ret=false;
	if ((in.compare("1")==0)||(in.compare("true")==0)||(in.compare("t")==0)) { ret=true; }
	return ret;
}


void StringProcessing::reference(string in , double *val, string *ref)
{
  size_t pos=in.find_first_of("@");
  if (pos!=string::npos){
    *ref=in.erase(0,pos+1);
    return;
  } else {
    *ref="";
    *val=atof(in.c_str());
    return;
  }
}
