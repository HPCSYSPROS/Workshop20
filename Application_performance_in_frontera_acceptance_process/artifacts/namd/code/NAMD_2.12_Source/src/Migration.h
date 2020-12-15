/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef MIGRATION_H
#define MIGRATION_H

#include "InfoStream.h"
#include "NamdTypes.h"
#include "PatchTypes.h"

typedef FullAtom MigrationElem;

typedef ResizeArray<MigrationElem> MigrationList;
typedef ResizeArrayIter<MigrationElem> MigrationListIter;

struct MigrationInfo {
  PatchID destPatchID;
  NodeID  destNodeID;
  MigrationList mList;
};

#endif // MIGRATION_H

