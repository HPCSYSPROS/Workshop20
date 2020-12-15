/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef GLOBALMASTERSYMMETRY_H
#define GLOBALMASTERSYMMETRY_H

#include "GlobalMaster.h"
#include "PDBData.h"
#include <map>
#include <vector>
#include "Matrix4Symmetry.h"

class GlobalMasterSymmetry : public GlobalMaster {
public:
  GlobalMasterSymmetry();
  ~GlobalMasterSymmetry();

private:
//  map <int, vector <Matrix4Symmetry> > matrices;
  std::vector <Matrix4Symmetry> matrices;
  std::map < int, Matrix4Symmetry > backmatrices;

  bool gluInvertMatrix(const BigReal [16], BigReal invOut[16]);
  void initialTransform();
  void backTransform();
 // void alignMonomers();
  void determineAverage();  
  void parseMatrix(int, char fileName []);
  void calculate();
  void parseAtoms(const char *file, int, int);

  std::map <int, std::map <int, int> > kdmap; //domain/atomid/k
  std::map <int, Position> positions; //atomid/position
  std::map <int, std::vector < BigReal * > > averagePos;
  std::map <int, BigReal *> backavg;
  std::map <int, std::vector<int>  > dmap;
  std::map <int, BigReal * > posmap;
  std::map <int, BigReal * > startmap;
  std::map <int, BigReal > kmap;
  std::map <int, std::vector <int> > simmap;
  std::map <int, int> bmap;


  BigReal K;
  const char *symmetrykfile;
  int currentStep, firstStep, lastStep, firstFullStep, lastFullStep;
  bool scaleForces;
};
#endif
