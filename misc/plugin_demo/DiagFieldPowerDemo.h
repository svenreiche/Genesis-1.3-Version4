#ifndef DIAGFIELDHOOKEDDEMO_H
#define DIAGFIELDHOOKEDDEMO_H

#include "DiagnosticPlugin.h"
#include "DiagnosticHookS.h"
#include "DiagUtil.h"

class DiagFieldHookedDemo: public DiagFieldHookedBase, public DiagUtil {
public:
	DiagFieldHookedDemo();
	~DiagFieldHookedDemo();

	void doit(DiagFieldHookData *);
	bool get_infos(DiagFieldHookInfos *);

private:
	// remember the names of the resources here
	std::vector<const char *> my_obj_names;
};

#endif
