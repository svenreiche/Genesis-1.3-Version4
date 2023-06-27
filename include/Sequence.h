// sequence.h: interface for the base clase of sequences such as Hammersley
//             or randomnumber
//
//////////////////////////////////////////////////////////////////////


#ifndef __GENESIS_SEQUENCE__
#define __GENESIS_SEQUENCE__
#include <cmath>

//#include "RandomU.h"
#include "Inverfc.h"

class Sequence
{
public:
	virtual ~Sequence() {};
	virtual double getElement() = 0;
	virtual void set(unsigned int) = 0;
};

class RandomU;

// derived sequences

class SequenceConst : public Sequence{
public:
    ~SequenceConst() {};
    double getElement();
    void set(unsigned int) {};
    void init(double );
private:
    double c;
    int i;
};
inline void SequenceConst::init(double c_in){ c = c_in; i = 0;}
inline double SequenceConst::getElement(){ return c;}

class SequencePower : public Sequence{
public:
    ~SequencePower() {};
    double getElement();
    void set(unsigned int) {};
    void init(double, double, double, int );
private:
    double c,dc,alpha;
    int  i,n;
};

inline void SequencePower::init(double c_in, double dc_in, double a_in, int n_in)
{
    c = c_in;dc = dc_in;alpha = a_in;i = 0;n = n_in;
}

inline double SequencePower::getElement() {
    i++;
    if (i < n) {
        return c;
    }
    return c + dc * pow(static_cast<double>(i - n), alpha);
}

class SequenceRandom : public Sequence
{
public:
    ~SequenceRandom(){};
    double getElement();
    void set(unsigned int) {}
    void init(double, double, int, bool);

private:
    double c0,dc;
    int seed;
    bool gauss;
    RandomU *seq;
    Inverfc erf;
};


#endif 
