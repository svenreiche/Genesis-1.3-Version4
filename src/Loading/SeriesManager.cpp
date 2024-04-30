//
// Created by reiche on 4/21/23.
// Additional code by clechner, Nov-2023
//

#include <iostream>
#include <map>
#include <string>
#include <algorithm>
#include <cctype>

#include <mpi.h>

#include "SeriesManager.h"
#include "Sequence.h"

using namespace std;

SeriesManager::SeriesManager() {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
}

void SeriesManager::add(string tag, Sequence *seq)
{
    // Inform user if there are upper-case characters in the sequence name
    // (even if the lattice file contains undulator elements with the
    // undulator parameter "aw = @AW", the function
    // LatticeParser::extractParameterValue is called with val="@aw",
    // this applies for instance to git commit id 3aa7ef8, date: 2024-04-18)
    bool has_uppercase = any_of(
        tag.begin(), tag.end(),
        [](unsigned char c) { return(isupper(c)); });
    if(has_uppercase && (0==rank_)) {
        cout << "Info: sequence name " << tag << " contains upper-case characters, this may cause issues" << endl;
    }

    bool exists = this->check(tag);
    if((exists) && (0==rank_)) {
        cout << "Info: sequence with name " << tag << " already defined, replacing it" << endl;
    }
    catalogue[tag]=seq;
}

bool SeriesManager::check(string tag)
{
    auto end = catalogue.end();
    auto ele = catalogue.find(tag);
    bool exists = (ele!=end);
    return exists;
}

double SeriesManager::getElement(string tag)
{
    // check if tag is in the list
    bool exists = this->check(tag);
    if (!exists) {
        if(0==rank_) {
            cout << "Cannot find requested series " << tag << endl;
        }
        abort(); // implement good error handling (currently it's developer code)
    }
    auto ele = catalogue.find(tag);
    // evaluate the Series
    Sequence *p_seq = ele->second;
    double v = p_seq->getElement();
    return(v);
}

