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


// base class for calculating values for the output. It must be inherited by all classes
//class DiagBase{
//public:
//    virtual unsigned int getCount() = 0;  // function to return number of calculated values per diagnostic step
//    virtual std::vector< std::string> getTags() = 0;
//    virtual std::vector<double> getValues() = 0
//};


// base class for beam diangostics
class DiagBeam{
protected:
    unsigned int npar {4};
    std::map<std::string,std::string>  tags;

private:
    unsigned int nharm {1};
    void defineTags();
public:
    ~DiagBeam() = default;
    DiagBeam() noexcept {this->defineTags();}; // this is needed since the harmonics can be changed
    auto getTags() { return tags;};
    void getValues(Beam *, std::map<std::string,std::vector<double> > &, int);
};


//class DiagField : public DiagBase{
//    unsigned int npar {1};
//public:
//    unsigned int getCount() {return npar;};
//    std::vector< std::string> getTags() { return {"power"};};
//};



// the class to manage the various subclasses
class Diagnostic{
    std::array<DiagBeam*,1> dbeam = {new DiagBeam()};
//    std::array<std::unique_ptr<DiagBase>,1> dfield = {std::unique_ptr<DiagBase>(new DiagField())};
    int nz = 1;
    int ns = 1;
    int iz = 0;

public:
    Diagnostic();
    virtual ~Diagnostic();
    void init(int,int);
    void calc(Beam *, std::vector<Field*> *,Undulator *);
    std::map<std::string,std::vector<double> > val;
    std::map<std::string,std::string > units;
};


#endif //GENESIS_1_3_VERSION4_DIAGNOSTIC_H
