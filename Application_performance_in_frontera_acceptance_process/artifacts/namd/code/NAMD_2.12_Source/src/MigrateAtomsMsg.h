/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Migration messages are sent from HomePatch to HomePatch
   with lists of atoms and atom information (if any) that
   need to be migrated.  A message must be sent from a
   neighbor even if null so that the HomePatch knows
   what atoms it will have before commencing a positionsReady()
   to its Computes.
*/

#ifndef MIGRATEATOMSMSG_H
#define MIGRATEATOMSMSG_H

#include "charm++.h"

#include "NamdTypes.h"
#include "SortedArray.h"
#include "Migration.h"
#include "PatchMgr.decl.h"

// Message which stores list of atoms and their data
// which are to be migrated from one patch to another.
// This message does not contain information that will change asynchronously
// It does not need to be prepacked
class MigrateAtomsMsg : public CMessage_MigrateAtomsMsg {
public:
  NodeID  fromNodeID;
  PatchID srcPatchID;
  PatchID destPatchID;
  MigrationList migrationList;

  MigrateAtomsMsg(void) { ; }

  // pack and unpack functions
  static void* pack(MigrateAtomsMsg* msg);
  static MigrateAtomsMsg* unpack(void *ptr);
};

class MigrateAtomsCombinedMsg : public CMessage_MigrateAtomsCombinedMsg
{
public:
  NodeID fromNodeID;
  ResizeArray<PatchID> srcPatchID;
  ResizeArray<PatchID> destPatchID;
  ResizeArray<int> numAtoms;
  int totalAtoms;
  MigrationList migrationList;

  MigrateAtomsCombinedMsg(void);
  ~MigrateAtomsCombinedMsg(void) { };

  void add(PatchID source, PatchID destination, MigrationList &m);
  void distribute(void);

  // pack and unpack functions
  static void* pack(MigrateAtomsCombinedMsg *msg);
  static MigrateAtomsCombinedMsg* unpack(void *ptr);
};

#endif

