/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef COMPUTEGLOBALMSGS_H
#define COMPUTEGLOBALMSGS_H

#include "charm++.h"

#include "NamdTypes.h"
#include "Lattice.h"
#include "ComputeMgr.decl.h"

#if 0
class ComputeGlobalConfigMsg : public CMessage_ComputeGlobalConfigMsg {
public:
  // data members
  AtomIDList aid;
  AtomIDList gdef;  // group definitions

  // constructor and destructor
  ComputeGlobalConfigMsg(void);
  ~ComputeGlobalConfigMsg(void);

  // pack and unpack functions
  static void* pack(ComputeGlobalConfigMsg *msg);
  static ComputeGlobalConfigMsg* unpack(void *ptr);
};
#endif


class ComputeGlobalDataMsg : public CMessage_ComputeGlobalDataMsg {
public:
  // data members
  int step;
  int count;
  AtomIDList aid;
  PositionList p;
  PositionList gcom;  // group center of mass
  BigRealList gmass;  // group total mass
  AtomIDList fid;
  ForceList tf;
  ForceList gtf;  // group total force
  ResizeArray<Lattice> lat;

  // constructor and destructor
  ComputeGlobalDataMsg(void);
  ~ComputeGlobalDataMsg(void);

  // pack and unpack functions
  static void* pack(ComputeGlobalDataMsg *msg);
  static ComputeGlobalDataMsg* unpack(void *ptr);
};


class ComputeGlobalResultsMsg : public CMessage_ComputeGlobalResultsMsg {
public:
  // data members
  AtomIDList aid;
  ForceList f;
  ForceList gforce;  // force on group

  int seq;
  int totalforces;  // send total forces?
  int reconfig;

  /* If <resendCoordinates> is 1, this message indicates a request for
     another set of coordinates (a ComputeGlobalDataMessage) during
     this timestep.  It may be 1 even if reconfig was not set,
     though there is no particular reason to do that.  A 1 here also
     indicates that the ComputeGlobal should ignore any forces
     included in this message, and wait instead for the next Results
     Message to come in. */
  int resendCoordinates;
  
  AtomIDList newaid;
  AtomIDList newgdef;

  // constructor and destructor
  ComputeGlobalResultsMsg(void);
  ~ComputeGlobalResultsMsg(void);

  // pack and unpack functions
  static void* pack(ComputeGlobalResultsMsg *msg);
  static ComputeGlobalResultsMsg* unpack(void *ptr);

};


#endif

