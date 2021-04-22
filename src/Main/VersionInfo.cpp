#include "VersionInfo.h"

VersionInfo::VersionInfo()  {};
VersionInfo::~VersionInfo() {};

int VersionInfo::Major(void) {
	return(4);
}

int VersionInfo::Minor(void) {
	return(4);
}

int VersionInfo::Rev(void) {
	return(0);
}

bool VersionInfo::isBeta(void) {
	return(true);
}
