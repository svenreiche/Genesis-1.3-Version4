//
// Created by reiche on 4/21/23.
// Additional code by clechner, Nov-2023
//

#include <iostream>
#include <map>
#include <string>

#include "SeriesManager.h"
#include "Sequence.h"

using namespace std;

void SeriesManager::add(std::string tag, Sequence *seq){
    catalogue[tag]=seq;
}

double SeriesManager::getElement(std::string tag){
    // check if tag is in the list
    auto end = catalogue.end();
    auto ele = catalogue.find(tag);
    bool exists = (ele!=end);
    if (!exists) {
        cout << "Cannot find requested series " << tag << endl;
        abort(); // implement good error handling (currently it's developer code)
    }

    // evaluate the Series
    Sequence *p_seq = ele->second;
    double v = p_seq->getElement();
    return(v);
}

