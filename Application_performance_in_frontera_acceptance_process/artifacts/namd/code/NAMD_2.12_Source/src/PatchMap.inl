/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef PATCHMAP_INL
#define PATCHMAP_INL

#include "PatchMap.h"
#include "AtomMap.h"

//----------------------------------------------------------------------
inline PatchID PatchMap::assignToPatch(Position p, const Lattice &l)
{
  int ai, bi, ci;
  ScaledPosition s = l.scale(p);
  ai = (int)floor(((BigReal)aDim)*((s.x-aOrigin)/aLength));
  bi = (int)floor(((BigReal)bDim)*((s.y-bOrigin)/bLength));
  ci = (int)floor(((BigReal)cDim)*((s.z-cOrigin)/cLength));
  return pid(ai,bi,ci);
}

//----------------------------------------------------------------------
#define MODULO(I,J) ( (I)<0 ? ((J)-(-1*(I))%(J))%(J) : (I)%(J) )

inline int PatchMap::pid(int aIndex, int bIndex, int cIndex)
{
  if ( aPeriodic ) aIndex = MODULO(aIndex,aDim);
  else
  {
    if ( aIndex < 0 ) aIndex = 0;
    if ( aIndex >= aDim ) aIndex = aDim - 1;
  }
  if ( bPeriodic ) bIndex = MODULO(bIndex,bDim);
  else
  {
    if ( bIndex < 0 ) bIndex = 0;
    if ( bIndex >= bDim ) bIndex = bDim - 1;
  }
  if ( cPeriodic ) cIndex = MODULO(cIndex,cDim);
  else
  {
    if ( cIndex < 0 ) cIndex = 0;
    if ( cIndex >= cDim ) cIndex = cDim - 1;
  }
  return ((cIndex*bDim)+bIndex)*aDim + aIndex;
}

//----------------------------------------------------------------------
inline int PatchMap::downstream(int pid1, int pid2)
{
  register int ds;

  if ( pid1 == pid2 ) { ds = pid1; }

  else if ( pid1 == notUsed || pid2 == notUsed ) { ds =  notUsed; }

  else {
    register PatchData *pdat1 = &(patchData[pid1]);
    register PatchData *pdat2 = &(patchData[pid2]);

    // c
    register int k = pdat1->cIndex;
    register int k2 = pdat2->cIndex;
    if ( ( k ? k : cMaxIndex ) == k2 + 1 ) k = k2;

    // b
    register int j = pdat1->bIndex;
    register int j2 = pdat2->bIndex;
    if ( ( j ? j : bMaxIndex ) == j2 + 1 ) j = j2;

    // a
    register int i = pdat1->aIndex;
    register int i2 = pdat2->aIndex;
    if ( ( i ? i : aMaxIndex ) == i2 + 1 ) i = i2;

    ds = ((k*bDim)+j)*aDim + i;
  }

  return ds;
}

#endif

