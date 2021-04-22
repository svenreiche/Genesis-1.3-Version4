#ifndef __GENESIS_VERSIONINFO_H
#define __GENESIS_VERSIONINFO_H

class VersionInfo {
public:
	VersionInfo();
	~VersionInfo();

	int Major(void);
	int Minor(void);
	int Rev(void);
	bool isBeta(void);
};

#endif
