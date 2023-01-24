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


#ifdef FFTW
#include <fftw3.h>
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
    bool auxiliar {true};
    int harm {1};
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



// For field diagnostics: helper class for management of FFT-related resources
#ifdef FFTW
class FFTObj {
public:
	FFTObj(int);
	~FFTObj();

	// The fftw_plan is a pointer, delete copy constructor and copy
	// assignment operator.
	// This avoids dangling pointers in object copies (fftw_plan contains
	// the data locations that were used to prepare the plan).
	FFTObj(const FFTObj&) = delete;
	FFTObj& operator= (const FFTObj&) = delete;

	std::complex<double> *in_ {nullptr};
	std::complex<double> *out_ {nullptr};

	// according to FFTW documentation, fftw_plan is an "opaque pointer type" (https://www.fftw.org/fftw3_doc/Using-Plans.html , 19.01.2023)
	fftw_plan p_;

private:
	int rank_;
	int ngrid_;
};
// map holding FFT-related resources for the values of ngrid
typedef std::map<int,FFTObj *> map_fftobj;
#endif

class DiagFieldBase{
protected:
    std::map<std::string, OutputInfo>  tags;  // holds a map with tags and units
    std::map<std::string,bool> filter;        // general map to store the selected flags
    bool global {false};
#ifdef FFTW
    map_fftobj fftobj;
#endif

public:
    ~DiagFieldBase() = default;
    DiagFieldBase() = default;
    virtual std::map<std::string,OutputInfo> getTags(FilterDiagnostics &filter) =0;
    virtual void getValues(Field *, std::map<std::string,std::vector<double> > &, int) =0;
//#ifdef FFTW
//    virtual void destroyFFTPlan()=0;
//#endif
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
// template for user defined beam diagnostics
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
    std::array<DiagBeamBase*,2> dbeam = {new DiagBeam(), new DiagBeamUser()};
    std::array<DiagFieldBase*,2> dfield = {new DiagField(), new DiagFieldUser()};

    int nz = 1;
    int ns = 1;
    int iz = 0;
    int noff = 0;
    int ntotal=1;


public:
    Diagnostic() = default;
    virtual ~Diagnostic() = default;
    void init(int,int, int, int,int,bool,bool, FilterDiagnostics &);
    void calc(Beam *, std::vector<Field*> *,double);
    void writeToOutputFile(std:: string, Beam *, vector<Field*> *, Undulator *);
    std::vector<std::map<std::string,std::vector<double> > > val;
    std::vector<std::map<std::string,std::string > >units;
    std::vector<std::map<std::string, bool> > single;
    std::vector<double> zout;
private:
    void addOutput(int groupID, std::string key, std::string unit, std::vector<double> &data);
    bool time,scan;
};


#endif //GENESIS_1_3_VERSION4_DIAGNOSTIC_H
