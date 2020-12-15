/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Compute reciprocal space contribution to pressure profile using
   Ewald sums.
*/

#ifndef COMPUTEEWALD_H
#define COMPUTEEWALD_H

#include "ComputeHomePatches.h"
#include "NamdTypes.h"
#include "ComputeMgr.decl.h"
 
class floatcomplex {
public:
  float r, i;
  floatcomplex() {}
  floatcomplex(float c) : r(c), i(0) {}
  floatcomplex(float a, float b) : r(a), i(b) {}
  ~floatcomplex() {}

  // complex conjuate
  floatcomplex star() const 
  {
    return floatcomplex(r, -i);
  }
  // square modulus
  float mod2() const 
  {
    return r*r + i*i;
  }
  // scalar assignment
  floatcomplex &operator=(const float &c) 
  {
    r = c;
    i = 0;
    return *this;
  }
  // complex add
  floatcomplex &operator+=(const floatcomplex &d) 
  {
    r += d.r;
    i += d.i;
    return *this;
  }
  // scalar multiply 
  floatcomplex &operator*=(const float &c) 
  {
    r *= c;
    i *= c;
    return *this;
  }
  // complex multiply
  floatcomplex &operator*=(const floatcomplex &d) 
  {
    float re = r*d.r - i*d.i;
    float im = r*d.i + i*d.r;
    r = re;
    i = im;
    return *this;
  }
};

class ComputeEwaldMsg : public CMessage_ComputeEwaldMsg {
public:
  float *eik;
};

class EwaldParticle;
class ComputeMgr;
class SubmitReduction;

class ComputeEwald : public ComputeHomePatches {
public:
  ComputeEwald(ComputeID, ComputeMgr*);
  virtual ~ComputeEwald();
  void doWork();
  void recvData(ComputeEwaldMsg *);
  void recvResults(ComputeEwaldMsg *);

  int getMasterNode() const { return masterNode; }

private:
  ComputeMgr *comm;
  SubmitReduction *reduction;
  int masterNode;
  int numWorkingPes;
  int recvCount;

  // Cached copy of coordinates and charges
  EwaldParticle *localAtoms;
  int *localPartitions;
  int numLocalAtoms;

  // pressure profile arrays
  int pressureProfileSlabs;
  int pressureProfileAtomTypes;
  float pressureProfileMin, pressureProfileThickness;
  float *pressureProfileData;

  // maximum k in x, y, z directions
  int kxmax, kymax, kzmax;

  // number of k vectors is (kxmax+1) + (2*kymax+1) + (2*kzmax+1)
  int ktot;

  // Ewald coefficient
  float kappa;

  // summed eik
  float *eiktotal;

  // space for temporary arrays
  floatcomplex *eiky, *eikz;
  float *expx, *expy, *expz;
  float *Qk;

  // support for multiple atom types
  // table of mappings from atom type to grid number.  Table for atom
  // of type i starts at offset n*i, where n=#of atom types.
  int *gridsForAtomType;
  int numAtomTypes;


  // cache the Lattice in doWork.
  Lattice lattice;

  // Store structure factor from localAtoms in given array
  void compute_structurefactor(float *);

  // compute reciprocal space contribute to pressure using
  // summation over local particles.
  void computePprofile(const float *eik) const; 
};

#endif

