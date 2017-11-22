// sequence.h: interface for the base clase of sequences such as Hammersley
//             or randomnumber
//
//////////////////////////////////////////////////////////////////////


#ifndef __GENESIS_SEQUENCE__
#define __GENESIS_SEQUENCE__


class Sequence
{
public:
	virtual double getElement() = 0;
	virtual void set(unsigned int) = 0;
};


#endif 
