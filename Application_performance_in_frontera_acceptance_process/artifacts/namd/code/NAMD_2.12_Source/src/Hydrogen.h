/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef HYDROGEN_H
#define HYDROGEN_H

#include "NamdTypes.h"
#include "SortableResizeArray.h"

// List maintaining the global atom indicies sorted by helix groups.
class HydrogenGroupID {
  public:
    AtomID atomID;      // global atom ID
    // isGP and atomsInGroup are determined when hydrogen bonds are found.
    int isGP;		// flag determining whether this atom is a group parent
    int atomsInGroup;   // positive number means parent of group.
                        // 0 means not in group.
    int waterVal;	// number of H bonded to O parent, 2 if water
    // although the Molecule object contains get_mother_atom(), we cannot
    // use it since Molecule.h, Node.h, and structure.h would have cyclical
    // include statements.
    int GPID;	// group parent ID, should be atomID if isGP is true
    // extension for migration groups
    int isMP;  // is this atom a migration group parent
    int MPID;  // migration group parent ID
    int atomsInMigrationGroup;
#ifdef MEM_OPT_VERSION
    int atomHoldIdxPar;
#endif

    HydrogenGroupID() {};
    ~HydrogenGroupID() {};

    int operator < (const HydrogenGroupID &a) const {
      int rval;
      int mp1 = ( isMP ? atomID : MPID );
      int mp2 = ( a.isMP ? a.atomID : a.MPID );
      if ( mp1 != mp2 ) {
        rval = ( mp1 < mp2 );
      } else { // same migration group, compare hydrogen groups
        int gp1 = ( isGP ? atomID : GPID );
        int gp2 = ( a.isGP ? a.atomID : a.GPID );
        if ( gp1 != gp2 ) {
          rval = ( gp1 < gp2 );
        // in the same group, check for group parent
        } else if ( a.isGP ) {  // compare to self should return 0
          rval = 0;
        } else if ( isGP ) {
          rval = 1;
        // neither is group parent, check for drude
        } else if ( a.atomID == gp1 + 1 ) {  // compare to self return 0
          rval = 0;
        } else if ( atomID == gp1 + 1 ) {
          rval = 1;
        // neither is drude
        } else {
	  rval = (atomID < a.atomID);
        }
      }
      return rval;
    }
};

typedef SortableResizeArray<HydrogenGroupID> HydrogenGroup ;

#endif

