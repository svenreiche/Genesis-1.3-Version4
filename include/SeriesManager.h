//
// Created by reiche on 4/21/23.
//

#include <string>
#include <map>

#ifndef GENESIS_1_3_VERSION4_SERIESMANAGER_H
#define GENESIS_1_3_VERSION4_SERIESMANAGER_H

#include "Sequence.h"

class SeriesManager {
    public:
        void add(std::string, Sequence *);
        double getElement(std::string);
    private:
        std::map<std::string,Sequence *> catalogue;
};

inline void SeriesManager::add(std::string tag, Sequence *seq){
    catalogue[tag]=seq;
}

inline double SeriesManager::getElement(std::string tag){
    // check if tag is in the list !!
    return catalogue[tag]->getElement();
}

#endif //GENESIS_1_3_VERSION4_SERIESMANAGER_H
