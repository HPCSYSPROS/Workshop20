/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Forwards atoms to master node for force evaluation.
*/

#if !defined(GLOBALMASTERFREEENERGY_H)
#define GLOBALMASTERFREEENERGY_H

class Molecule;
class SimParameters;

class GlobalMasterFreeEnergy : public GlobalMaster {
 public:
  GlobalMasterFreeEnergy();
  ~GlobalMasterFreeEnergy();
 private:
  virtual void calculate();
  void user_calculate();

  void initialize();
  void user_initialize();
  void update();
  Molecule *molecule;
  SimParameters *simParams;
  char *config;
 ARestraintManager  m_RestraintManager;
  ALambdaManager     m_LambdaManager;
public:
  // These all return -1 on error.
  int getAtomID(const char *segid, int resid, const char *aname);
  int getNumAtoms(const char* segid, int resid); // 0 on error
  int getAtomID(const char *segid, int resid, int index);
  double getMass(int atomid);
  int requestAtom(int atomid);
  int getPosition(int atomid, Position &position);
  int addForce(int atomid, Force force);
};

#endif
