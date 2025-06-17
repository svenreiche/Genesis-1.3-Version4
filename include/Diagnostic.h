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

#ifdef FFTW
#include <fftw3.h>
#include "FFTObj.h"
#endif

//#include "Field.h"
//#include "Beam.h"
//#include "Undulator.h"
class Beam;
class Field;
class Undulator;
class Setup;

#include "DiagnosticBase.h"
#ifdef USE_DPI
  #include "DiagnosticHook.h"
#endif


//------------------------------------
// genesis official class for beam diagnostics


class DiagBeam: public DiagBeamBase{
private:
    unsigned int nharm {1};    // beam specific for calculating harmonics in the bunching
    bool exclharm {false};
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

#ifdef FFTW
    void cleanup_FFT_resources(void);
    int obtain_FFT_resources(int ngrid, complex<double> **in, complex<double> **out, fftw_plan *pp);
#endif
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

//----------------------------------------
// template for user defined field diagnostics
class DiagFieldUser: public DiagFieldBase{
public:
    DiagFieldUser() = default;
    ~DiagFieldUser() = default;
    std::map<std::string,OutputInfo> getTags(FilterDiagnostics &filter);
    void getValues(Field *, std::map<std::string,std::vector<double> > &, int);
};

// ----------------------------------------------
// class which manages all the diagnostic classes, acting as a wrapper for them, standardizing the interface.
class Diagnostic{
    std::vector<DiagBeamBase*>  dbeam  {new DiagBeam(), new DiagBeamUser()};
    std::vector<DiagFieldBase*> dfield {new DiagField(), new DiagFieldUser()};


    int nz = 1;
    int ns = 1;
    int iz = 0;
    int noff = 0;
    int ntotal=1;


public:
    Diagnostic();
    // virtual ~Diagnostic() = default;
    virtual ~Diagnostic();
    void init(int,int, int, int,int,bool,bool, FilterDiagnostics &);
    void calc(Beam *, std::vector<Field*> *,double);
    bool writeToOutputFile(Beam *, vector<Field*> *, Setup *, Undulator *);

    std::vector<std::map<std::string,std::vector<double> > > val;
    std::vector<std::map<std::string,std::string > >units;
    std::vector<std::map<std::string, bool> > single;
    std::vector<double> zout;

    bool add_beam_diag(DiagBeamBase *);
    bool add_field_diag(DiagFieldBase *);

private:
    void addOutput(int groupID, std::string key, std::string unit, std::vector<double> &data);

    bool time,scan;

    bool diag_can_add {true};

    int my_rank_;
};


#endif //GENESIS_1_3_VERSION4_DIAGNOSTIC_H
