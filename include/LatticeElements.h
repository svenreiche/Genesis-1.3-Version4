#ifndef __GENESIS__LATTICEELEMENTS__
#define __GENESIS__LATTICEELEMENTS__

#include <string>
#include <iostream>

//enum LatticeElement {Quadrupole, Drift, Dipole, Undulator, Marker};

using namespace std;

class Element
{
 public:
  Element(){};
  ~Element(){};
  double z,l;
  string type;
     //enum LatticeElement type;
};

class ID : public Element
{
 public:
  ID(){};
  ~ID(){};
  double aw,lambdau,kx,ky,nwig,ax,ay;
  double gradx,grady,pdadx,pdady; // gradients
  double paw,pkx,pky,phase;   // perpendicular phase
  bool helical;
};

class Phaseshifter : public Element
{
 public:
  Phaseshifter(){};
  ~Phaseshifter(){};
  double phi;
};
class Quadrupole : public Element
{
 public:
  Quadrupole(){};
  ~Quadrupole(){};
  double k1,dx,dy;
};

class Drift : public Element
{
 public:
  Drift(){};
  ~Drift(){};
};

class Marker : public Element
{
 public:
  Marker(){};
  ~Marker(){};
  int action;  // this is a bitwise coded element
};

class Chicane : public Element
{
 public:
  Chicane(){};
  ~Chicane(){};
  double lb,ld,delay;
}; 

class Corrector : public Element
{
 public:
  Corrector(){};
  ~Corrector(){};
  double cx,cy;
}; 

#endif
