//
// Created by reiche on 4/21/23.
//

#ifndef GENESIS_1_3_VERSION4_SERIESPARSER_H
#define GENESIS_1_3_VERSION4_SERIESPARSER_H

#include "SeriesManager.h"
#include "StringProcessing.h"
#include <map>
#include <string>

class SeriesParser : public StringProcessing{
public:
    bool init(int, std::map<std::string,std::string> *, std::string, SeriesManager *);
    bool initConst(int, std::map<std::string,std::string> *, SeriesManager *);
    bool initPolynom(int, std::map<std::string,std::string> *, SeriesManager *);
    bool initPower(int, std::map<std::string,std::string> *, SeriesManager *);
    bool initRandom(int, std::map<std::string,std::string> *, SeriesManager *);
    bool initList(int, std::map<std::string,std::string> *, SeriesManager *);
    bool initFileList(int, std::map<std::string,std::string> *, SeriesManager *);

private:
	bool readfile(const std::string, std::vector<double> &, int);

    void usageConst();
    void usagePolynom();
    void usagePower();
    void usageRandom();
    void usageList();
    void usageFileList();
};


#endif //GENESIS_1_3_VERSION4_SERIESPARSER_H
