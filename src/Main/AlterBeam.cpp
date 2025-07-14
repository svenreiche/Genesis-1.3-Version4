//
// Created by reiche on 4/17/25.
//

#include "AlterBeam.h"

AlterBeam::AlterBeam()= default;
AlterBeam::~AlterBeam()= default;
void AlterBeam::usage(){
    std::cout << "List of keywords for ALTER_BEAM" << std::endl;
    std::cout << "&alter_beam" << std::endl;
    std::cout << " double dgamma = 0 / reference" << std::endl;
    std::cout << " double phase = 0 / reference" << std::endl;
    std::cout << " double R56    = 0" << std::endl;
    std::cout << " double lambda = 800e-9 " << std::endl;
    std::cout << "&end" << std::endl << std::endl;
}

bool AlterBeam::init(int rank, int size, std::map<std::string,std::string> *arg, Beam *beam, Setup *setup, Time *time, Profile *prof) {
    double lambda = setup->getReferenceLength();   // reference length for theta
    double gammaref = setup->getReferenceEnergy();
    bool one4one = setup->getOne4One();            // check for one4one simulations
    bool dotime = time->isTime();


    // default value
    double dgamma = 0;
    double phase = 0;
    std::string dgammaref;
    std::string phaseref;
    double r56 = 0;
    double lambdaM = 800e-9;

    auto end = arg->end();
    if (arg->find("dgamma") != end) {
        this->reference(arg->at("dgamma"), &dgamma, &dgammaref);
        arg->erase(arg->find("dgamma"));
    }
    if (arg->find("phase") != end) {
        this->reference(arg->at("phase"), &phase, &phaseref);
        arg->erase(arg->find("phase"));
    }
    if (arg->find("r56") != end) {
        r56 = strtod(arg->at("r56").c_str(), nullptr);
        arg->erase(arg->find("r56"));
    }
    if (arg->find("lambda") != end) {
        lambdaM = strtod(arg->at("lambda").c_str(), nullptr);
        arg->erase(arg->find("lambda"));
    }

    if (!arg->empty()) {
        if (rank == 0) {
            std::cout << "*** Error: Unknown elements in &alter_beam" << std::endl;
            AlterBeam::usage();
        }
        return false;
    }

    // checking for wrong profiles
    std::string wrongProf;
    if (!prof->check(dgammaref)) { wrongProf = dgammaref; }
    if (!prof->check(phaseref)) { wrongProf = phaseref; }
    if (!wrongProf.empty()) {
        if (rank == 0) {
            std::cout << "*** Error: Unknown profile reference in &alter_beam: " << wrongProf << std::endl;
        }
        return false;
    }

    bool doAverage = !(dotime & one4one); // non-average calculation only possible with one4one simulation with an explicit time-profile.

    if (doAverage) {
        double rat = round(lambdaM / lambda);
        lambdaM = rat * lambda;
        if (rank == 0) {
            std::cout << "Adjusting modulation wavelength to " << lambdaM * 1e9 << " nm to yield integer harmonic of h = " << static_cast<int>(rat)
                      << std::endl;
        }
    }
    // get s position along the bunches
    std::vector<double> s;
    auto nslice = time->getPosition(&s);
    auto nsliceNode = beam->beam.size();
    auto ioffset = nsliceNode * rank;
    if (nsliceNode == 0) {
        std::cout << "*** Error: Beam must be defined before applying alteration " << std::endl;
        return false;
    }
    auto kr = 4. * asin(1.) / lambdaM;
    auto invk0 = lambda/(4. * asin(1.));
    for (size_t j = 0; j < nsliceNode; j++) {
        auto i = j + ioffset;
        auto delgam = prof->value(s.at(i), dgamma, dgammaref);
        auto rphase = prof->value(s.at(i), phase, phaseref);
        auto npart = beam->beam.at(j).size();
        if (doAverage) {
            for (int ip = 0; ip < npart; ip++) {
                auto theta = beam->beam.at(j).at(ip).theta;
                auto gamma = beam->beam.at(j).at(ip).gamma;
                // note that theta is treated as theta at longer wavelength!
                gamma += delgam * sin(theta + rphase);
                theta += kr * r56 * (gamma - gammaref) / gammaref;
                beam->beam.at(j).at(ip).gamma = gamma;
                beam->beam.at(j).at(ip).theta = theta * lambdaM / lambda;  // apply harmonic conversion
            }
        } else {
            for (int ip = 0; ip < npart; ip++) {
                auto theta = beam->beam.at(j).at(ip).theta;
                auto gamma = beam->beam.at(j).at(ip).gamma;
                auto sloc = s.at(i) + theta *invk0;
                auto delgamloc = prof->value(sloc, dgamma, dgammaref);
                auto rphaseloc = prof->value(sloc, phase, phaseref);
                gamma += delgamloc * sin(kr * sloc + rphaseloc);
                theta += lambdaM/lambda*kr * r56 * (gamma - gammaref) / gammaref;
                beam->beam.at(j).at(ip).gamma = gamma;
                beam->beam.at(j).at(ip).theta = theta;
            }
        }

    }
    return true;

}