/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "InfoStream.h"
#include "Node.h"
#include "Molecule.h"
#include "NamdTypes.h"

#include "GlobalMaster.h"
#include "GlobalMasterSMD.h"
#include "PDB.h"
#include "PDBData.h"

/* XXX necessary?
#include <iostream.h>
#if !defined(WIN32) || defined(__CYGWIN__)
#include <strstream.h>
#else
#include <strstrea.h>
#endif
#include <string.h>
*/

//#define DEBUGM
#define MIN_DEBUG_LEVEL 1
#include "Debug.h"

GlobalMasterSMD::GlobalMasterSMD(BigReal spring_constant,
	       BigReal transverse_spring_constant, BigReal velocity,
               const Vector direction, int output_frequency,
	       int first_timestep, const char *filename,
	       int numAtoms) {
  DebugM(3,"initialize called\n");
  moveVel = velocity;
  moveDir = direction;
  outputFreq = output_frequency;
  k = spring_constant;
  k2 = transverse_spring_constant;
  currentTime = first_timestep;

  parseAtoms(filename,numAtoms);
  iout << iINFO << requestedGroups()[0].size() << " SMD ATOMS\n" << endi;
  DebugM(1,"done with initialize\n");
}

void GlobalMasterSMD::parseAtoms(const char *file, int numTotalAtoms) {
  DebugM(3,"parseAtoms called\n");
  PDB smdpdb(file);
  int numatoms = smdpdb.num_atoms();
  if (numatoms < 1) 
    NAMD_die("No atoms found in SMDFile\n");
  if (numatoms != numTotalAtoms)
    NAMD_die("The number of atoms in SMDFile must be equal to the total number of atoms in the structure!");

  // add a single group
  modifyRequestedGroups().resize(1);

  // add all the required atoms to the group

  // also compute the center of mass of these atoms
  cm.x = cm.y = cm.z = 0;
  Molecule *mol = Node::Object()->molecule; // to get masses
  BigReal imass = 0; 

  for (int i=0; i<numatoms; i++) {
#ifdef MEM_OPT_VERSION
    PDBCoreData *atom = smdpdb.atom(i);
#else
    PDBAtom *atom = smdpdb.atom(i); // get an atom from the file
#endif    
    if (atom->occupancy()) { // if occupancy is not 0, then add it!
      // add the atom to the list
      modifyRequestedGroups()[0].add(i);

      // compute the center of mass
      BigReal mass = mol->atommass(i); 
      cm.x += atom->xcoor()*mass;
      cm.y += atom->ycoor()*mass;
      cm.z += atom->zcoor()*mass;
      imass += mass; 
    }
  }
  if (imass == 0) // we didn't find any!
    NAMD_die("SMDFile contained no SMD atoms (atoms w/ nonzero occupancy)\n");

  // now compute the center of mass
  cm /= imass;
  DebugM(1,"done with parseAtoms\n");
}

GlobalMasterSMD::~GlobalMasterSMD() { }

void GlobalMasterSMD::calculate() {
  DebugM(3,"calculate called\n");

  // make extra-sure we have been set up correctly
  if(!requestedGroups().size() == 1)
    NAMD_die("Internal error: uninitialized!");

  // set the correct forces
  Position curcm = *getGroupPositionBegin(); // get the center of mass
  DebugM(1,"Current CM "<<cm<<"\n");
  BigReal diff = (curcm - cm)*moveDir;
  // second term below is along transverse direction: -diff*moveDir + (diff*moveDir - (curcm-cm)) = -(curcm-cm)
  // so if k = k2 and moveVel = 0, we see that f = -k * (curcm - cm), the desired result
  Force f = k*(moveVel*currentTime - diff)*moveDir + k2*(diff*moveDir - (curcm - cm));
  modifyGroupForces().resize(1);
  modifyGroupForces()[0] = f;

  // print some output sometimes
  if (currentTime % outputFreq == 0)
    output(currentTime, curcm, f);

  // keep track of the number of timesteps
  currentTime++;

  DebugM(1,"done with calculate: force: "<<f<<"\n");
}

void GlobalMasterSMD::output(int t, Position p, Force f) {
  DebugM(3,"output called\n");

  if (t % (100*outputFreq) == 0) 
    iout << "SMDTITLE: TS   CURRENT_POSITION         FORCE\n" << endi; 
  iout << "SMD  " << t << ' ' << p << ' ' << f*PNPERKCALMOL << '\n' << endi;

  DebugM(1,"done with output\n");
}
