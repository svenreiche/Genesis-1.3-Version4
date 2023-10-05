#ifndef DIAGFIELDHOOKEDDEMO_H
#define DIAGFIELDHOOKEDDEMO_H

#include "DiagnosticPlugin.h"
#include "DiagnosticHookS.h"
#include "DiagUtil.h"

class DiagBeamHookedDemo: public DiagBeamHookedBase, public DiagUtil {
public:
	DiagBeamHookedDemo();
	~DiagBeamHookedDemo();

	void doit(DiagBeamHookData *);
	bool get_infos(DiagBeamHookInfos *);

private:
	// remember the names of the resources here
	std::vector<const char *> my_obj_names;
};

#endif
