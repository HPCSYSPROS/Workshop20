/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Methods are primarily for pack(ing) and unpack(ing) messages for Charm.
*/

#include "InfoStream.h"
#include "Migration.h"
#include "MigrateAtomsMsg.h"
#include "NamdTypes.h"
// #define DEBUGM
#define MIN_DEBUG_LEVEL 3
#include "Debug.h"

#include "PatchMgr.decl.h"
#include "PatchMap.h"
#include "HomePatch.h"
#include "packmsg.h"


PACK_MSG(MigrateAtomsMsg,
  PACK(fromNodeID);
  PACK(srcPatchID);
  PACK(destPatchID);
  PACK_RESIZE(migrationList);
)


MigrateAtomsCombinedMsg::MigrateAtomsCombinedMsg(void)
{
  fromNodeID = CkMyPe();
  totalAtoms = 0;
}

void MigrateAtomsCombinedMsg::
	add(PatchID source, PatchID destination, MigrationList &m)
{
  srcPatchID.add(source);
  destPatchID.add(destination);
  int n = m.size();
  numAtoms.add(n);
  totalAtoms += n;
  for ( int i = 0; i < n; ++i )
  {
    migrationList.add(m[i]);
  }
}


void MigrateAtomsCombinedMsg::distribute(void)
{
  int n = srcPatchID.size();
  int m = 0;
  for ( int i = 0; i < n; ++i )
  {
    MigrateAtomsMsg *msg = new MigrateAtomsMsg;
    msg->fromNodeID = fromNodeID;
    msg->srcPatchID = srcPatchID[i];
    msg->destPatchID = destPatchID[i];
    int l = numAtoms[i];
    {
      DebugM(3,"Distributing " << l << " atoms to patch " << msg->destPatchID << "\n");
      msg->migrationList.resize(l);
      for ( int j = 0; j < l; ++j ) msg->migrationList[j] = migrationList[m+j];
      m += l;
    }
    PatchMap::Object()->homePatch(msg->destPatchID)->depositMigration(msg);
  }
  if ( m != totalAtoms ) NAMD_bug("MigrateAtomsCombinedMsg::distribute bad atom count");
}


PACK_MSG(MigrateAtomsCombinedMsg,
  PACK(fromNodeID);
  PACK(totalAtoms);
  PACK_RESIZE(srcPatchID);
  PACK_RESIZE(destPatchID);
  PACK_RESIZE(numAtoms);
  PACK_RESIZE(migrationList);
)

