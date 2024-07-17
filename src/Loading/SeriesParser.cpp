//
// Created by reiche on 4/21/23.
//
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <vector>
#include "SeriesParser.h"
#include "Sequence.h" // contains implementation of several classes

bool SeriesParser::init(int rank, std::map<std::string,std::string> *arg, std::string element, SeriesManager *sm){
    // slit for the different cases
    if (element.compare("&sequence_const")==0){
        return this->initConst(rank, arg, sm);
    }
    if (element.compare("&sequence_polynom")==0){
        return this->initPolynom(rank, arg, sm);
    }
    if (element.compare("&sequence_power")==0){
        return this->initPower(rank, arg, sm);
    }
    if (element.compare("&sequence_random")==0){
        return this->initRandom(rank, arg, sm);
    }
    if (element.compare("&sequence_list")==0){
        return this->initList(rank, arg, sm);
    }
    if (element.compare("&sequence_filelist")==0){
        return this->initFileList(rank, arg, sm);
    }
    return false;
}

bool SeriesParser::initConst(int rank, std::map<std::string,std::string> *arg, SeriesManager *sm){
    std::string label="";
    double c0=0;
    std::map<std::string,std::string>::iterator end=arg->end();

    if (arg->find("label")!=end){label = arg->at("label");  arg->erase(arg->find("label"));}
    if (arg->find("c0")!=end)   {c0    = atof(arg->at("c0").c_str());  arg->erase(arg->find("c0"));}

    if (arg->size()!=0) {
        if (rank == 0) {
            std::cout << "*** Error: Unknown elements in &sequence_const" << std::endl;
            this->usageConst();
        }
        return false;
    }
    if (label.size()<1) {
        if (rank==0){
            std::cout << "*** Error: Label not defined in &sequence_const" << std::endl; this->usageConst();
        }
        return false;
    }
    if (rank==0){
        std::cout << "Adding sequence with label: " <<label << std::endl;
    }
    SequenceConst *seq = new SequenceConst;
    seq->init(c0);
    sm -> add(label,seq);
    return true;
}

bool SeriesParser::initPolynom(int rank, std::map<std::string,std::string> *arg, SeriesManager *sm){
    std::string label="";
    double c0=0,c1=0,c2=0,c3=0,c4=0;
    std::map<std::string,std::string>::iterator end=arg->end();

    if (arg->find("label")!=end){label = arg->at("label");  arg->erase(arg->find("label"));}
    if (arg->find("c0")!=end)   {c0    = atof(arg->at("c0").c_str());  arg->erase(arg->find("c0"));}
    if (arg->find("c1")!=end)   {c1    = atof(arg->at("c1").c_str());  arg->erase(arg->find("c1"));}
    if (arg->find("c2")!=end)   {c2    = atof(arg->at("c2").c_str());  arg->erase(arg->find("c2"));}
    if (arg->find("c3")!=end)   {c3    = atof(arg->at("c3").c_str());  arg->erase(arg->find("c3"));}
    if (arg->find("c4")!=end)   {c4    = atof(arg->at("c4").c_str());  arg->erase(arg->find("c4"));}

    if (arg->size()!=0) {
        if (rank == 0) {
            std::cout << "*** Error: Unknown elements in &sequence_polynom" << std::endl;
            this->usagePolynom();
        }
        return false;
    }
    if (label.size()<1) {
        if (rank==0){
            std::cout << "*** Error: Label not defined in &sequence_polynom" << std::endl;
            this->usagePolynom();
        }
        return false;
    }
    if (rank==0){
        std::cout << "Adding sequence with label: " <<label << std::endl;
    }
    SequencePolynom *seq = new SequencePolynom;
    seq->init(c0,c1,c2,c3,c4);
    sm -> add(label,seq);
    return true;
}

bool SeriesParser::initPower(int rank, std::map<std::string,std::string> *arg, SeriesManager *sm){
    std::string label="";
    double c0=0;
    double dc=0;
    double alpha=0;
    int n0=1;
    std::map<std::string,std::string>::iterator end=arg->end();

    if (arg->find("label")!=end){label = arg->at("label");  arg->erase(arg->find("label"));}
    if (arg->find("c0")!=end)   {c0    = atof(arg->at("c0").c_str());  arg->erase(arg->find("c0"));}
    if (arg->find("dc")!=end)   {dc    = atof(arg->at("dc").c_str());  arg->erase(arg->find("dc"));}
    if (arg->find("alpha")!=end){alpha = atof(arg->at("alpha").c_str());  arg->erase(arg->find("alpha"));}
    if (arg->find("n0")!=end)   {n0    = atoi(arg->at("n0").c_str());  arg->erase(arg->find("n0"));}

    if (arg->size()!=0) {
        if (rank == 0) {
            std::cout << "*** Error: Unknown elements in &sequence_power" << std::endl;
            this->usagePower();
        }
        return false;
    }
    if (label.size()<1) {
        if (rank==0){
            std::cout << "*** Error: Label not defined in &sequence_power" << std::endl; this->usageConst();
        }
        return false;
    }
    if (rank==0){
        std::cout << "Adding sequence with label: " <<label << std::endl;
    }

    SequencePower *seq = new SequencePower;
    seq->init(c0,dc,alpha,n0);
    sm -> add(label,seq);
    return true;
}

bool SeriesParser::initRandom(int rank, std::map<std::string,std::string> *arg, SeriesManager *sm){
    std::string label="";
    double c0=0;
    double dc=0;
    int seed;
    bool gauss;
    std::map<std::string,std::string>::iterator end=arg->end();

    if (arg->find("label")!=end){label = arg->at("label");  arg->erase(arg->find("label"));}
    if (arg->find("c0")!=end)   {c0    = atof(arg->at("c0").c_str());  arg->erase(arg->find("c0"));}
    if (arg->find("dc")!=end)   {dc    = atof(arg->at("dc").c_str());  arg->erase(arg->find("dc"));}
    if (arg->find("seed")!=end) {seed  = atoi(arg->at("seed").c_str());  arg->erase(arg->find("seed"));}
    if (arg->find("normal")!=end){gauss= atob(arg->at("normal").c_str());  arg->erase(arg->find("normal"));}

    if (arg->size()!=0) {
        if (rank == 0) {
            std::cout << "*** Error: Unknown elements in &sequence_random" << std::endl;
            this->usageRandom();
        }
        return false;
    }

    if (label.size()<1) {
        if (rank==0){
            std::cout << "*** Error: Label not defined in &sequence_random" << std::endl; this->usageConst();
        }
        return false;
    }
    if (rank==0){
        std::cout << "Adding sequence with label: " <<label << std::endl;
    }

    SequenceRandom *seq = new SequenceRandom;
    seq->init(c0,dc,seed,gauss);
    sm -> add(label,seq);
    return true;
}

bool SeriesParser::initList(int rank, std::map<std::string,std::string> *arg, SeriesManager *sm){
    std::string label="";
    double def =0;
    std::string valstr ="";
    std::map<std::string,std::string>::iterator end=arg->end();

    if (arg->find("label")!=end){label = arg->at("label");  arg->erase(arg->find("label"));}
    if (arg->find("default")!=end)   {def = atof(arg->at("default").c_str());  arg->erase(arg->find("default"));}
    if (arg->find("val")!=end){valstr = arg->at("val");  arg->erase(arg->find("val"));}

    std::vector<std::string> vals;
    this->chop(valstr,&vals);
    vector<double> fval;
    for (auto & ele : vals){
        if (ele.size() < 1){
            if (rank == 0) {
                std::cout << "*** Error: empty list element for val in &sequence_list" << std::endl;
                this->usageList();
            }
            return false;
        }
        fval.push_back(atof(ele.c_str()));
    }

    if (fval.size() < 1) {
        if (rank == 0) {
            std::cout << "*** Error: empty list for val in &sequence_list" << std::endl;
            this->usageList();
        }
        return false;
    }

    if (arg->size()!=0) {
        if (rank == 0) {
            std::cout << "*** Error: Unknown elements in &sequence_list" << std::endl;
            this->usageList();
        }
        return false;
    }
    if (label.size()<1) {
        if (rank==0){
            std::cout << "*** Error: Label not defined in &sequence_list" << std::endl; this->usageConst();
        }
        return false;
    }
    if (rank==0){
        std::cout << "Adding sequence with label: " <<label << std::endl;
    }
    SequenceList *seq = new SequenceList;
    seq->init(fval,def);
    sm -> add(label,seq);
    return true;
}



bool SeriesParser::initFileList(int rank, std::map<std::string,std::string> *arg, SeriesManager *sm)
{
    std::string label="";
    std::string filename="";
    vector<double> data;
    std::map<std::string,std::string>::iterator end=arg->end();

    if (arg->find("label")!=end){label   = arg->at("label"); arg->erase(arg->find("label"));}
    if (arg->find("file")!=end){filename = arg->at("file");  arg->erase(arg->find("file"));}

    if (arg->size()!=0) {
        if (rank == 0) {
            std::cout << "*** Error: Unknown elements in &sequence_filelist" << std::endl;
            this->usageFileList();
        }
        return false;
    }
    if (label.empty()) {
        if (rank==0){
            std::cout << "*** Error: Label not defined in &sequence_filelist" << std::endl; this->usageFileList();
        }
        return false;
    }
    if (filename.empty()) {
        if (rank==0){
            std::cout << "*** Error: Filename not defined in &sequence_filelist" << std::endl; this->usageFileList();
        }
        return false;
    }
    
    if(!readfile(filename, data, rank)) {
		if(rank==0) {
			cout << "*** Error while loading data from file in &sequence_filelist" << endl;
		}
		return(false);
	}
    if (rank==0){
        std::cout << "Adding sequence with label: " << label << ", data file: " << filename << ", length of data set is: " << data.size() << endl;
    }

    double default_value=0.; // Default value is currently 0 (cannot be adjusted)
    SequenceList *seq = new SequenceList;
    seq->init(data,default_value);
    sm -> add(label,seq);
    return true;
}
bool SeriesParser::readfile(const string fn, vector<double> &v, int rank)
{
	ifstream ifs;

	ifs.open(fn);
	if(!ifs) {
		if(rank==0) {
			cout << "Error: cannot open file " << fn << " for read access" << endl;
		}
		return(false);
	}
	if(rank==0) {
		cout << "Reading data from file " << fn << endl;
	}

	unsigned int linecntr=0;
	string linebuf;
	while(getline(ifs,linebuf))
	{
		linecntr++;

		stringstream ss(linebuf);
		double data=-1.;
		if (!(ss >> data)) {
			if(rank==0) {
				cout << "Issue parsing line " << linecntr << endl;
			}
			return(false);
		}
		v.push_back(data);
		// cout << "buf=" << linebuf << ", data=" << data << endl;
	}

	ifs.close();
	return(true);
}



// the individual usages calls

void SeriesParser::usageConst(){
    cout << "List of keywords for SEQUENCE_CONST" << endl;
    cout << "&sequence_const" << endl;
    cout << " string label = <empty>" << endl;
    cout << " double c0 = 0" << endl;
    cout << "&end" << endl << endl;
    return;
}
void SeriesParser::usagePolynom(){
    cout << "List of keywords for SEQUENCE_POLYNOM" << endl;
    cout << "&sequence_polynom" << endl;
    cout << " string label = <empty>" << endl;
    cout << " double c0 = 0" << endl;
    cout << " double c1 = 0" << endl;
    cout << " double c2 = 0" << endl;
    cout << " double c3 = 0" << endl;
    cout << " double c4 = 0" << endl;
    cout << "&end" << endl << endl;
    return;
}
void SeriesParser::usagePower(){
    cout << "List of keywords for SEQUENCE_POWER" << endl;
    cout << "&sequence_power" << endl;
    cout << " string label = <empty>" << endl;
    cout << " double c0 = 0" << endl;
    cout << " double dc = 0" << endl;
    cout << " double alpha = 0" << endl;
    cout << " int n0 = 1" << endl;
    cout << "&end" << endl << endl;
    return;
}
void SeriesParser::usageRandom(){
    cout << "List of keywords for SEQUENCE_RANDOM" << endl;
    cout << "&sequence_random" << endl;
    cout << " string label = <empty>" << endl;
    cout << " double c0 = 0" << endl;
    cout << " double dc = 0" << endl;
    cout << " int seed = 100" << endl;
    cout << " bool normal = true" << endl;
    cout << "&end" << endl << endl;
    return;
}
void SeriesParser::usageList(){
    cout << "List of keywords for SEQUENCE_LIST" << endl;
    cout << "&sequence_list" << endl;
    cout << " string label = <empty>" << endl;
    cout << " double val = [<empty>]" << endl;
    cout << " double default = 0" << endl;
    cout << "&end" << endl << endl;
    return;
}
void SeriesParser::usageFileList(){
    cout << "List of keywords for SEQUENCE_FILELIST" << endl;
    cout << "&sequence_filelist" << endl;
    cout << " string label = <empty>" << endl;
    cout << " string file = <empty>" << endl;
    cout << "&end" << endl << endl;
    return;
}
