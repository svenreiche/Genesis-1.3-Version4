//
// Created by reiche on 4/21/23.
//

#ifndef GENESIS_1_3_VERSION4_SERIESMANAGER_H
#define GENESIS_1_3_VERSION4_SERIESMANAGER_H

#include <string>
#include <map>

class Sequence;

class SeriesManager {
public:
	SeriesManager();
        void add(std::string, Sequence *);
        double getElement(std::string);
        bool check(std::string);
private:
        int rank_;
        std::map<std::string,Sequence *> catalogue;
};

#endif //GENESIS_1_3_VERSION4_SERIESMANAGER_H
