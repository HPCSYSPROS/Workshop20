/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Forwards atoms to master node for force evaluation.
*/

#ifndef GLOBALMASTEREASY_H
#define GLOBALMASTEREASY_H

#include "NamdTypes.h"

class Molecule;
class SubmitReduction;

class GlobalMasterEasy : public GlobalMaster {
protected:
  GlobalMasterEasy(const char *the_config_name);
  virtual ~GlobalMasterEasy();

  int getAtomID(const char *segid, int resid, const char *aname);
  int getNumAtoms(const char* segid, int resid); // 0 on error
  int getAtomID(const char *segid, int resid, int index);
  double getMass(int atomid);
  int requestAtom(int atomid);
  int getPosition(int atomid, Position &position);
  int addForce(int atomid, Force force);
  void addEnergy(BigReal);

  virtual void easy_init(const char *);
  virtual void easy_calc(void);

  char *config;
private:

  void initialize();
  virtual void calculate();

  Molecule *molecule;
  SubmitReduction *reduction;

  char *configName;
  BigReal energy;

};

#endif
