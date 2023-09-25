#ifndef __GENESIS_REGPLUGIN_H
#define __GENESIS_REGPLUGIN_H

#include <iostream>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <map>
#include <fstream>
#include <complex>

#include "StringProcessing.h"
#include "Setup.h"

class AddPluginFieldDiag: public StringProcessing {
public:
	AddPluginFieldDiag();
	~AddPluginFieldDiag();

	bool init(int, int, std::map<std::string,std::string> *, Setup *);

private:
	void usage(void);
};

#endif // __GENESIS_REGPLUGIN_H
