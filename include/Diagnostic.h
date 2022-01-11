//
// Created by reiche on 06.01.22.
//

#ifndef GENESIS_1_3_VERSION4_DIAGNOSTIC_H
#define GENESIS_1_3_VERSION4_DIAGNOSTIC_H

#include <string>
#include <array>
#include <vector>
#include <tuple>
#include <functional>
#include <memory>

#include "Field.h"
#include "Beam.h"
#include "Undulator.h"

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
    bool global {true};
    bool spatial {true};
    bool energy {true};
    bool current {false};
    bool auxiliar {true};
    int harm {1};
};

struct FilterField{
    bool global {true};
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

class DiagBeamBase{
protected:
    std::map<std::string, OutputInfo>  tags;  // holds a map with tags and units
    std::map<std::string,bool> filter;        // general map to store the selected flags
    bool global {false};
public:
    ~DiagBeamBase() = default;
    DiagBeamBase() = default; // this is needed since the harmonics can be changed
    virtual std::map<std::string,OutputInfo> getTags(FilterDiagnostics &filter) =0;
    virtual void getValues(Beam *, std::map<std::string,std::vector<double> > &, int) =0;
};

class DiagFieldBase{
protected:
    std::map<std::string, OutputInfo>  tags;  // holds a map with tags and units
    std::map<std::string,bool> filter;        // general map to store the selected flags
    bool global {false};
public:
    ~DiagFieldBase() = default;
    DiagFieldBase() = default; // this is needed since the harmonics can be changed
    virtual std::map<std::string,OutputInfo> getTags(FilterDiagnostics &filter) =0;
    virtual void getValues(Field *, std::map<std::string,std::vector<double> > &, int) =0;
};


//------------------------------------
// genesis official class for beam diagnostics

class DiagBeam: public DiagBeamBase{
private:
    unsigned int nharm {1};    // beam specific for calculating harmonics in the bunching
public:
    ~DiagBeam() = default;
    DiagBeam() = default; // this is needed since the harmonics can be changed
    std::map<std::string,OutputInfo> getTags(FilterDiagnostics &filter);
    void getValues(Beam *, std::map<std::string,std::vector<double> > &, int) ;
};

class DiagField: public DiagFieldBase{
public:
    ~DiagField() = default;
    DiagField() = default; // this is needed since the harmonics can be changed
    std::map<std::string,OutputInfo> getTags(FilterDiagnostics &filter);
    void getValues(Field *, std::map<std::string,std::vector<double> > &, int) ;
};

//----------------------------------------
// template for user defined beam diagnostics
class DiagBeamUser: public DiagBeamBase{
private:
    unsigned int nharm {1};
public:
    DiagBeamUser() = default;
    ~DiagBeamUser() = default;
    std::map<std::string,OutputInfo> getTags(FilterDiagnostics &filter);
    void getValues(Beam *, std::map<std::string,std::vector<double> > &, int);
};


// ----------------------------------------------
// class which manages all the diagnostic classes, acting as a wrapper for them, standardizing the interface.
class Diagnostic{
    std::array<DiagBeamBase*,2> dbeam = {new DiagBeam(), new DiagBeamUser()};
    std::array<DiagFieldBase*,1> dfield = {new DiagField()};
    int nz = 1;
    int ns = 1;
    int iz = 0;

public:
    Diagnostic() = default;
    virtual ~Diagnostic() = default;
    void init(int,int,int);
    void calc(Beam *, std::vector<Field*> *,Undulator *);
    std::vector<std::map<std::string,std::vector<double> > > val;
    std::vector<std::map<std::string,std::string > >units;
    std::vector<std::map<std::string,double> > svaldouble;
    std::vector<std::map<std::string,int> > svalint;
    std::vector<double> zout;
};


#endif //GENESIS_1_3_VERSION4_DIAGNOSTIC_H
