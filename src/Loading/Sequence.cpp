//
// Created by reiche on 4/21/23.
//

#include "Sequence.h"
#include "RandomU.h"

void SequenceRandom::init(double c_in, double dc_in, int seed_in, bool gauss_in){
    c0 = c_in;
    dc= dc_in;
    seed = seed_in;
    gauss = gauss_in;
    seq = new RandomU (seed);
}

double SequenceRandom::getElement(){
    if (gauss){
        return c0+dc*(erf.value(2*seq->getElement()));
    }
    return c0+dc*(2*seq->getElement()-1);
}

