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

enum PluginType {PLUGIN_FIELD, PLUGIN_BEAM};

class AddPluginCommon: public StringProcessing {
public:
	AddPluginCommon() = default;
	virtual ~AddPluginCommon() = default;

	bool init(int, int, std::map<std::string,std::string> *, Setup *);

private:
	void usage(void);
	
	virtual PluginType id(void) = 0;
	virtual std::string nameliststr(void) = 0;
};

class AddPluginFieldDiag: public AddPluginCommon {
public:
	AddPluginFieldDiag() = default;
	~AddPluginFieldDiag() = default;

	// bool init(int, int, std::map<std::string,std::string> *, Setup *);

private:
	// void usage(void);
	PluginType id(void);
	std::string nameliststr(void);
};

class AddPluginBeamDiag: public AddPluginCommon {
public:
	AddPluginBeamDiag() = default;
	~AddPluginBeamDiag() = default;

	// bool init(int, int, std::map<std::string,std::string> *, Setup *);

private:
	// void usage(void);
	PluginType id(void);
	std::string nameliststr(void);
};

#endif // __GENESIS_REGPLUGIN_H
