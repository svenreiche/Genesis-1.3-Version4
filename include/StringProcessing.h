/*
 *  StringProcessing.h
 *  Genesis
 *
 *  Created by Sven Reiche on 9/10/10.
 *  Copyright 2010 Paul Scherrer Institut. All rights reserved.
 *
 */

#ifndef __GENESIS_STRINGPROCESSING__
#define __GENESIS_STRINGPROCESSING__

#include <string>
#include <vector>
#include <iostream>
#include <stdlib.h>

using namespace std;

class StringProcessing {
public:
	StringProcessing();
	~StringProcessing();
protected:
	void trim(string &); 
	void chop(string, vector<string> *);  
	bool atob(string);
	void reference(string, double *, string *);
	vector<string> slist;   // temporary working array

};

#endif
