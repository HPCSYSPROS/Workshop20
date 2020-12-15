/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef COMPUTECONSFORCEMSG_H
#define COMPUTECONSFORCEMSG_H

#include "charm++.h"

#include "NamdTypes.h"
#include "ComputeMgr.decl.h"

class ComputeConsForceMsg : public CMessage_ComputeConsForceMsg {
public:
  // data members
  AtomIDList aid;
  ForceList f;

  // constructor and destructor
  ComputeConsForceMsg() {}
  ~ComputeConsForceMsg() {}

  // pack and unpack functions
  static void* pack(ComputeConsForceMsg *msg);
  static ComputeConsForceMsg* unpack(void *ptr);
};

#endif

