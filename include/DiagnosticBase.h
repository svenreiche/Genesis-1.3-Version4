//
// Fundamental diagnostic classes from Diagnostic.h
//

#ifndef GENESIS_1_3_VERSION4_DIAGNOSTICBASE_H
#define GENESIS_1_3_VERSION4_DIAGNOSTICBASE_H

#include <string>
#include <array>
#include <vector>
#include <tuple>
#include <functional>
#include <memory>

#include "Field.h"
#include "Beam.h"
#include "Undulator.h"


#ifdef FFTW
#include <fftw3.h>
#include "FFTObj.h"
#endif


extern bool MPISingle;

// supplemental info for any field to be calculated. It helps to allocate the right size and automize
// the writing to the hdf5 file

struct OutputInfo{
    bool global {false};
    bool once {false};
    std::string units {" "};
};

// filters to control output. A true value means that the output is included in the output file
// these values are initialized and changed in the namelist setup, track and the upcoming output

struct FilterBeam{
    bool global {false};
    bool spatial {true};
    bool energy {true};
    bool current {false};
    bool twiss {false};
    bool auxiliar {true};
    int harm {1};
    bool exclharm {false};
};

struct FilterField{
    bool global {false};
    bool spatial {true};
    bool fft {true};
    bool intensity {true};
};

struct FilterDiagnostics{
    FilterBeam beam;
    FilterField field;
};

//---------------------------------------
// base class for diagnostics

class DiagBase{
protected:
    void storeValue(std::map<std::string,std::vector<double> >&val, std::string, unsigned long, double);
};

inline void DiagBase::storeValue(std::map<std::string,std::vector<double> >&container,std::string field, unsigned long index, double value)
{
    if (container.find(field) != container.end()){container[field].at(index)=value;}
}


class DiagBeamBase : public DiagBase{
protected:
    std::map<std::string, OutputInfo>  tags;  // holds a map with tags and units
    std::map<std::string,bool> filter;        // general map to store the selected flags
    bool global {false};

public:
    virtual ~DiagBeamBase() = default;
    DiagBeamBase() = default; // this is needed since the harmonics can be changed
    virtual std::map<std::string,OutputInfo> getTags(FilterDiagnostics &filter) =0;
    virtual void getValues(Beam *, std::map<std::string,std::vector<double> > &, int) =0;
};

class DiagFieldBase : public DiagBase{
protected:
    std::map<std::string, OutputInfo>  tags;  // holds a map with tags and units
    std::map<std::string,bool> filter;        // general map to store the selected flags
    bool global {false};
#ifdef FFTW
    map_fftobj fftobj;
#endif

public:
    virtual ~DiagFieldBase() = default;
    DiagFieldBase() = default;
    virtual std::map<std::string,OutputInfo> getTags(FilterDiagnostics &filter) =0;
    virtual void getValues(Field *, std::map<std::string,std::vector<double> > &, int) =0;
};





#endif //GENESIS_1_3_VERSION4_DIAGNOSTICBASE_H
