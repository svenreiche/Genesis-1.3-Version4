#ifndef G4__DIAGPLUGIN_H
#define G4__DIAGPLUGIN_H

// header file with definitions of data structures for information exchange
#include "DiagnosticHookS.h"

class DiagFieldHookedBase {
public:
//	DiagFieldHookedBase();
	virtual ~DiagFieldHookedBase() {};

	virtual bool get_infos(DiagFieldHookInfos *) = 0;
	virtual void doit(DiagFieldHookData *) = 0;
};

#endif // G4__DIAGPLUGIN_H
