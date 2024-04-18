//
// Created by reiche on 4/21/23.
//

#include <string>
#include <map>

#ifndef GENESIS_1_3_VERSION4_SERIESMANAGER_H
#define GENESIS_1_3_VERSION4_SERIESMANAGER_H

class Sequence;

class SeriesManager {
    public:
        void add(std::string, Sequence *);
        double getElement(std::string);
        bool check(std::string);
    private:

        std::map<std::string,Sequence *> catalogue;
};

#endif //GENESIS_1_3_VERSION4_SERIESMANAGER_H
