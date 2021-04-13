#ifndef __GENESIS_PROF_H
#define __GENESIS_PROF_H

#include <map>
#include <string>
#include <stdio.h>

typedef long long mytime_t;

class Prof
{
public:
	Prof();
	~Prof();
	bool init(int, int, std::map<std::string,std::string> *);
	void report(void);

private:
	void usage(void);

	bool atob(std::string);
	int getnano(mytime_t *);

	static bool report_cmp(const std::pair<std::string, mytime_t> &p1, const std::pair<std::string, mytime_t> &p2);
	void report_core(FILE *);

	mytime_t t0_; // point of reference, i.e. when this class was constructed
	std::map<std::string, mytime_t> times_;
};

#endif /* __GENESIS_PROF_H */
