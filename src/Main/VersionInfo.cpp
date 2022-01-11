#include "VersionInfo.h"
#include "build_info.h"

VersionInfo::VersionInfo()  {};
VersionInfo::~VersionInfo() {};

int VersionInfo::Major(void) {
	return(4);
}

int VersionInfo::Minor(void) {
	return(5);
}

int VersionInfo::Rev(void) {
	return(1);
}

bool VersionInfo::isBeta(void) {
	return(false);
}

const char * VersionInfo::BuildInfo(void) {
	return(build_info());
}
