//
// Created by reiche on 10.01.22.
//

#include "Diagnostic.h"

//  This source file works as a template  to allow users to add additional output.
// the definition is already given in the header file "Diagnostic.h"

//--------------------------------------------
// 1st mandatory definition: the registration of the output to be calculated.
// The core elements, which have to been defined are:
//      tags -> a map of output parameters. The key is the name as it appears as the dataset name in the output file name
//              Each tag holds a struct with two bools and a string. The bool "global" indicates that a collective calculation
//              over all MPI ranks are done. It will appear in the output file in the sub griup  "global" .
//              The bool "once" indicates that the parameter is calculated only at the beginning of the calculation for
//              the first integration step. The string holds the description for the attribute in the output file, defining
//              the unit of the output.
//     filter-> a map which can hold user defined flags (bools) for a possible filtering in the routing getValues. If some filters
//              are defined, depends on the user. The function is called with the argument filter_in, which holds some general flags
//              to control the output, defined in setup or track namelist.

std::map<std::string,OutputInfo> DiagBeamUser::getTags(FilterDiagnostics & filter_in) {

    tags.clear();
    filter.clear();

    // two variables inherited from the base class.
    nharm = filter_in.beam.harm;   // the limit to up which harmonics in the particle phase is analysed
    global = filter_in.beam.global; // user flag if global parameter are used or not

    // 2 example as template for user defined parameters:
    //         A: energy modulation (for each slice)
    //         B: total loss by external fields e.g. wakes (global variable)

    filter["template"] = false;     // <- set to true to enable output for the template parameters
    if (filter["template"]){
        tags["modulation"] = {false, false, "mc^2"};
        if (global) {
            tags["Global/loss"] = {true, false, "eV/m"};  // to avoid duplicated tags with same name,
                                                    // global should have a preceeding 'global/' in tag
        }
    }

    return tags;  // must be returned
}

void DiagBeamUser::getValues(Beam *beam, std::map<std::string,std::vector<double> >&val, int iz) {

    int ns = beam->beam.size();   // the number of slices per node
    int is = 0;    // counter for the current slice
    double g_norm = 0;   // global variable for normalizing the other global variables (here the current)
    double g_loss = 0;   // global variable for beam loss

    // loop over the slices in the electron bunch
    for (auto const &slice: beam->beam) {
        complex<double> emod = 0; // parameter for the energy modulation
        // loop over the particle in each slice
        for (auto const &par: slice) {
            emod + complex<double>(par.gamma * cos(par.theta), par.gamma * sin(par.theta));
        }

        double norm = 1.;
        if (!slice.empty()) {
            norm = 1. / static_cast<double>(slice.size()); // normalize with the number of particles per slice
        }
        emod *= norm;
        // gather info for the global variables
        g_norm += beam->current[is];
        g_loss += beam->current[is] * beam->eloss[is];

        // saving the calculated values into the records defined with getTags.
        // it is good practice to check whether a given tag exists, since they can vary with filter flags, such as for global parameters
        int idx = iz * ns + is;         // index for saving the data in the individual vectors of 'val'
        // for 'once' parameters it would be just idx=is
        // for 'global' parameters it would just iz
        if (val.find("modulation") != val.end()) { val["modulation"][idx] = abs(emod); }
    }

    //-------------------------------------------------
    // complete gathering from all nodes if MPI size is larger than one
    if (global) {
        int size = 1;
        if (!MPISingle) {
            MPI_Comm_size(MPI_COMM_WORLD, &size);
        }  // for future functionality to do scan functionality of time-dependent runs!
        // this will be better controlled in the future with a dedicated MPI class (RAII)
        if (size > 1) {  // more than one mode is running, we need to gather from all
            double temp;
            MPI_Allreduce(&g_norm, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            g_norm = temp;
            MPI_Allreduce(&g_loss, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            g_loss = temp;
        }
        double norm = 1;
        if (g_norm > 0) { norm = 1. / g_norm; } // normalize with the weighting function (for beam it is the current)
        g_loss *= norm;


        if (val.find("global/loss") != val.end()) { val["global/loss"][iz] = g_loss; }
    }
}