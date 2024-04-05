#include <fstream>
#include <iostream>
#include <cstdio>
#include <string>
#include <map>

#include "g4_tb_util.h"

using namespace std;

TB_Cfg::TB_Cfg() { }

void TB_Cfg::eat_whitespaces(string& s)
{
	const char *ws = " \t";
	size_t startpos = s.find_first_not_of(ws);
	size_t endpos   = s.find_last_not_of(ws);
	if ((startpos==string::npos) || (endpos==string::npos)) {
		s = "";
		return;
	}
	s = s.substr(startpos, endpos-startpos+1);
}

bool TB_Cfg::update_param(const string key, const string value)
{
	const map<string, double *> param_double {
		{"dgrid",  &dgrid},
		{"lambda", &lambda},
		{"power",  &power},
		{"w0",     &w0}
	};
	const map<string, int *> param_int {
		{"nslice", &nslice},
		{"nz",     &nz},
		{"ngrid",  &ngrid}
	};
	const map<string, string *> param_string {
		{"libfile",   &libfile},
		{"parameter", &parameter}
	};

	/* PARAMETERS WITH EXPECTED VALUE OF TYPE 'double' */
	auto end_d = param_double.end();
	auto ele_d = param_double.find(key);
	if(ele_d!=end_d) {
		double *pd = ele_d->second;
		int r = sscanf(value.c_str(), "%lg", pd);
		if(r==0) {
			cout << "unable to parse value assigned to parameter " << key << endl;
			return(false);
		}
		return(true);
	}

	/* ... 'int' */
	auto end_i = param_int.end();
	auto ele_i = param_int.find(key);
	if(ele_i!=end_i) {
		int *pi = ele_i->second;
		int r = sscanf(value.c_str(), "%d", pi);
		if(r==0) {
			cout << "unable to parse value assigned to parameter " << key << endl;
			return(false);
		}
		return(true);
	}

	/* ... 'string' */
	auto end_s = param_string.end();
	auto ele_s = param_string.find(key);
	if(ele_s!=end_s) {
		string *ps = ele_s->second;
		*ps = value;
		return(true);
	}

	/* unknown parameter */
	cout << "unknown parameter " << key << endl;
	return(false);
}

bool TB_Cfg::update_from_stream(ifstream& ifs)
{
	int lcntr=1;
	string linebuf;
	while(getline(ifs, linebuf)) {
		size_t pos = linebuf.find_first_of("=");
		if (pos==string::npos) {
			cout << "line " << lcntr << ": error parsing line \"" << linebuf << "\"" << endl;
			continue;
		}
		string left = linebuf.substr(0,pos);
		string right = linebuf.substr(pos+1,string::npos);
		eat_whitespaces(left);
		eat_whitespaces(right);
		// cout << "L:" << left << ", R:" << right << endl;
		if(!update_param(left,right)) {
			cout << "line " << lcntr << ": error processing parameter " << left << endl;
		}

		lcntr++;
	}

	return(true);
}
