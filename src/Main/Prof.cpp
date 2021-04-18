/*
 * Profiling namelist element for GENESIS4.
 *
 * Christoph Lechner, European XFEL GmbH, 13-Apr-2021
 *
 * For standalone testing, compile with -DSTANDALONE, for example:
 *    g++ -Wall -DSTANDALONE -Iinclude/ -o proftest src/Main/Prof.cpp
 */
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <map>
#include <string>
#include <vector>
#include <time.h>
#ifndef STANDALONE
 #include <mpi.h>
#endif
#ifdef STANDALONE
 #include <unistd.h> // for 'sleep'
#endif
#include "Prof.h"

using namespace std;

Prof::Prof()
{
	mytime_t t;

	getnano(&t);
	t0_          = t;
	times_["t0"] = t;
}
Prof::~Prof() {}

void Prof::usage(void)
{
	cout << "List of keywords for PROFILE" << endl;
	cout << "&profile" << endl;
	cout << " barrier = false" << endl;
	cout << " label = <empty>" << endl;
	cout << " writetofile = <empty>" << endl;
	cout << "&end" << endl << endl;
}


/* copy of StringProcessing::atob, it is the only function from this class that is needed here */
bool Prof::atob(string in){
	bool ret=false;
	if ((in.compare("1")==0)||(in.compare("true")==0)||(in.compare("t")==0)) { ret=true; }
	return ret;
}

int Prof::getnano(mytime_t *outnano)
{
	const long long nano2sec = 1000000000;
	mytime_t nano;
	struct timespec ts;

	if(clock_gettime(CLOCK_REALTIME, &ts)<0) {
		cout << "call to clock_gettime failed!" << endl;
		return(-1);
	}

	nano  = ts.tv_sec;
	nano *= nano2sec;
	nano += ts.tv_nsec;

	*outnano = nano;
	return(0);
}

/*
 * Conversion of nanosecond values to string based on integer arithmetics
 * => avoid rounding errors as the numbers are very large
 */
int Prof::nano2str(char *buf, const int buflen, mytime_t nano)
{
	const long long nano2sec = 1000000000;
	long long sec;
	long usec;

	sec  = nano / nano2sec;
	usec = (nano % nano2sec) / 1000;	// microseconds is sufficient precision
	usec = labs(usec);			// if given negative time (for instance deltat), show negative sign only before decimal point
	return(snprintf(buf, buflen, "%lld.%06ld", sec, usec));
}

bool Prof::report_cmp(const std::pair<std::string, mytime_t> &p1, const std::pair<std::string, mytime_t> &p2)
{
	if(p1.second < p2.second)
		return(true);
	else if(p1.second == p2.second)
		if(p1.first < p2.first)
			return(true);

	return(false);
}
void Prof::report_core(FILE *fout, bool pretty)
{
	vector< pair<string,mytime_t> > v;
	vector<string> col_labels;
	size_t max_len_label;
	vector<string> col_t1;
	size_t max_len_col1;
	vector<string> col_t2;
	size_t max_len_col2;
	vector<string> col_t3;
	size_t max_len_col3;

	/* sort data in map according to time (and add current time) */
	pair<string,mytime_t> tmp_p;
	tmp_p.first = string("Current time");
	if(getnano(&tmp_p.second)==0)
		v.push_back(tmp_p);
#if 0
	for(map<string,mytime_t>::const_iterator it = times_.begin();
	    it != times_.end();
	    ++it)
	{
		v.push_back(*it);
	}
#endif
	copy(times_.begin(), times_.end(), back_inserter(v));
	sort(v.begin(), v.end(), report_cmp);

	if(pretty==false)
	{
		for(int i=0; i<v.size(); i++)
		{
			const int buflen=1000;
			char buf_t[buflen], buf_since_t0[buflen], buf_dt[buflen];
			mytime_t this_time, deltat;

			this_time = v[i].second;
			deltat=0;
			if(i>0) {
				deltat = this_time-v[i-1].second;
			}

			nano2str(buf_t,        buflen, this_time);
			nano2str(buf_since_t0, buflen, this_time-t0_);
			nano2str(buf_dt,       buflen, deltat);

			fprintf(fout, "\"%s\",%s,%s,%s\n",
				v[i].first.c_str(), buf_t, buf_since_t0, buf_dt);
		}
		return;
	}

	/* first round prepares pretty-printing => determine needed column widths */
	string col_hdr_label("label");
	string col_hdr1("t [s]");
	string col_hdr2("t-t0 [s]");
	string col_hdr3("deltat [s]");
	max_len_label = col_hdr_label.length();
	max_len_col1 = col_hdr1.length();
	max_len_col2 = col_hdr2.length();
	max_len_col3 = col_hdr3.length();
	for(int i=0; i<v.size(); i++)
	{
		const int buflen=1000;
		char buf[buflen];
		size_t len_this_label;
		mytime_t this_time;
		string tmpstr;

		this_time = v[i].second;
		len_this_label = v[i].first.length();
		max_len_label = (len_this_label > max_len_label) ? len_this_label : max_len_label;
		col_labels.push_back(v[i].first);

		nano2str(buf, buflen, this_time);
		tmpstr = buf;
		max_len_col1 = (tmpstr.length() > max_len_col1) ? tmpstr.length() : max_len_col1;
		col_t1.push_back(tmpstr);

		nano2str(buf, buflen, this_time-t0_);
		tmpstr = buf;
		max_len_col2 = (tmpstr.length() > max_len_col2) ? tmpstr.length() : max_len_col2;
		col_t2.push_back(tmpstr);

		if(i>0) {
			nano2str(buf, buflen, this_time-v[i-1].second);
			tmpstr = buf;
		} else {
			/* no deltat in first row */
			tmpstr = "";
		}
		max_len_col3 = (tmpstr.length() > max_len_col3) ? tmpstr.length() : max_len_col3;
		col_t3.push_back(tmpstr);
	}

	string sep1(max_len_label, '=');
	string sep2(max_len_col1, '=');
	string sep3(max_len_col2, '=');
	string sep4(max_len_col3, '=');
	fprintf(fout, "%-*s   %-*s   %-*s   %-*s\n",
		static_cast<int>(max_len_label), col_hdr_label.c_str(),
	        static_cast<int>(max_len_col1),  col_hdr1.c_str(),
	        static_cast<int>(max_len_col2),  col_hdr2.c_str(),
	        static_cast<int>(max_len_col3),  col_hdr3.c_str());
	fprintf(fout, "%s   %s   %s   %s\n", sep1.c_str(), sep2.c_str(), sep3.c_str(), sep4.c_str());
	for(int i=0; i<v.size(); i++)
	{
		fprintf(fout, "%-*s   %*s   %*s   %*s\n",
			static_cast<int>(max_len_label), col_labels[i].c_str(),
		        static_cast<int>(max_len_col1),  col_t1[i].c_str(),
		        static_cast<int>(max_len_col2),  col_t2[i].c_str(),
		        static_cast<int>(max_len_col3),  col_t3[i].c_str());
	}
	fprintf(fout, "%s   %s   %s   %s\n", sep1.c_str(), sep2.c_str(), sep3.c_str(), sep4.c_str());
}
void Prof::report(void)
{
	report_core(stdout, true);
}

bool Prof::init(int mpirank, int mpisize, map<string,string> *arg)
{
	bool do_barrier;
	bool do_record;
	string id_label;
	bool do_write;
	string fn;
	map<string,string>::iterator end=arg->end();

	do_barrier = do_record = do_write = false;
	if(arg->find("barrier") != end) {
		do_barrier = atob(arg->at("barrier"));
		arg->erase(arg->find("barrier"));
	}
	if(arg->find("label") != end) {
		id_label = arg->at("label");
		arg->erase(arg->find("label"));
		do_record = true;
	}
	if(arg->find("writetofile") != end) {
		fn = arg->at("writetofile");
		arg->erase(arg->find("writetofile"));
		do_write = true;
	}
	
	/* if only known/valid arguments were supplied, then 'arg' should be empty now... */
	if(arg->size())
	{
		if(mpirank==0)
		{
			cout << "*** Error unknown elements in &profile" << endl;
			usage();
		}
		return false;
	}


	if(do_barrier)
	{
#ifndef STANDALONE
		MPI_Barrier(MPI_COMM_WORLD);
#else
		cout << "binary for standalone testing: MPI_Barrier disabled" << endl;
#endif
	}

	if(do_record)
	{
		mytime_t t;

		getnano(&t);
		times_[id_label] = t;
	}

	if((do_write) && (mpirank==0))
	{
		FILE *fout;

		if(NULL==(fout=fopen(fn.c_str(), "w")))
		{
			cout << "Profiling: cannot open file " << fn << " for writing (ignoring)" << endl;
			return true;
		}
		report_core(fout, false);
		fclose(fout);
	}

	return true;
}

#ifdef STANDALONE
int main(void)
{
	Prof p;
	map<string,string> arg;

	/* emulate calls from GenMain loop, records current time */
	arg["label"] = "Hallo";
	p.init(0,0,&arg);

	sleep(1);

	arg.clear();
	arg["label"] = "after sleep(1)";
	p.init(0,0,&arg);



	p.report();

	/* test file report functionality */
	arg.clear();
	arg["writetofile"] = "profdemoout.txt";
	p.init(0,0,&arg);

	return(0);
}
#endif
