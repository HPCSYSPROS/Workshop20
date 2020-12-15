/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Forwards atoms to master node for force evaluation.
*/

#include "InfoStream.h"
#include "NamdTypes.h"
#include "GlobalMaster.h"
#include "GlobalMasterMisc.h"

//#define DEBUGM
#define MIN_DEBUG_LEVEL 1
#include "Debug.h"

GlobalMasterMisc::GlobalMasterMisc() : GlobalMasterEasy("miscForcesScript") {
  // Initialize subclass
  easy_init(config);
  delete[] config;
}

GlobalMasterMisc::~GlobalMasterMisc() {
  ;
}

/* you may use the following in easy_init or easy_calc
  int getAtomID(const char *segid, int resid, const char *aname); // zero-based
  int getNumAtoms(const char* segid, int resid); // 0 on error
  int getAtomID(const char *segid, int resid, int index);
  double getMass(int atomid); // -1.0 on error
*/

/* you may use the following only in easy_init()
  int requestAtom(int atomid); // 0 on success, -1 on error
*/

void GlobalMasterMisc::easy_init(const char *config) {
  iout << iINFO << "  MISC FORCES CONFIG\n";
  iout << iINFO << "**********************\n";
  iout << config;
  iout << iINFO << "**********************\n" << endi;

  requestAtom(0);

  firstTime = true;
}

/* you may use the following only in easy_calc()
  int getPosition(int atomid, Position &position); // 0 on success, -1 on error
  int addForce(int atomid, Force force); // 0 on success, -1 on error
  void addEnergy(BigReal);
*/

/* EXAMPEL: harmonically restrain atom 0 to its original position */

void GlobalMasterMisc::easy_calc() {
  if(firstTime) {
    if(getPosition(0,originalPosition)) {
      NAMD_die("Couldn't get initial position.");
    }
    iout << iINFO << "Initial pos: " << originalPosition << "\n" << endi;
    firstTime = false;
  }

  Vector myp;
  BigReal k = 100.0;
  if(getPosition(0,myp)) {
    NAMD_die("Couldn't get position.");
  }
  iout << iINFO << "Atom 0 is at " << myp << "\n" << endi;
  addForce(0,-k * (myp-originalPosition));
  iout<<iINFO << "Adding force " << -k*(myp-originalPosition) << "\n" << endi;
  addEnergy(0.5 * k * (myp-originalPosition).length2());
}

