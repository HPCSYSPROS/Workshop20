/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "NamdTypes.h"
#include "GlobalMaster.h"
#include "GlobalMasterTest.h"

//#define DEBUGM
#define MIN_DEBUG_LEVEL 1
#include "Debug.h"

void GlobalMasterTest::calculate() {
  DebugM(3,"Test calculate called\n");
  if(requestedAtoms().size() != 2) {
    DebugM(3,"putting two atoms into request list");
    modifyRequestedAtoms().resize(0);
    modifyRequestedAtoms().add(0);
    modifyRequestedAtoms().add(1);
    modifyForcedAtoms().add(0);
    
    int s = appliedForces().size();
    modifyAppliedForces().resize(s+1);
    modifyAppliedForces().item(s).x = 2;
    modifyAppliedForces().item(s).y = 0;
    modifyAppliedForces().item(s).z = 0;
    return; // we can't expect to have the rights atoms yet, so go on
  }

  const AtomID *atomID_ptr = getAtomIdBegin();
  const AtomID *atomID_end = getAtomIdEnd();
  const Position *position_ptr = getAtomPositionBegin();
  
  if(atomID_end - atomID_ptr != 2) {
    DebugM(3,"found " << atomID_end - atomID_ptr << " atoms\n");
    NAMD_die("Wrong number of atoms.");
  }

  while(atomID_ptr != atomID_end) {
    DebugM(2,"Atom " << *atomID_ptr << " is at "
	   << position_ptr->x << " " << position_ptr->y << " "
	   << position_ptr->z << "\n");
    position_ptr ++;
    atomID_ptr ++;
  }
}

