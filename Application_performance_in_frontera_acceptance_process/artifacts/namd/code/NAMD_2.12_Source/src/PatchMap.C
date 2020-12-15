/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include <stddef.h>
#if !defined(WIN32) || defined(__CYGWIN__)
#include <unistd.h>
#endif
#include <stdio.h>

#include "InfoStream.h"
#include "ObjectArena.h"
#include "PatchMgr.h"
#include "PatchMap.inl"
#include "Patch.h"
#include "Lattice.h"
#include "HomePatchList.h"
#include "AtomMap.h"
#include "DataExchanger.h"
#include "memusage.h"
//#define DEBUGM
#define MIN_DEBUG_LEVEL 5
#include "Debug.h"

int *PatchMap::nPatchesOnNode = 0;
PatchMap::PatchData *PatchMap::patchData = 0;
ObjectArena<ComputeID> *PatchMap::computeIdArena = 0;

// Safe singleton creation
PatchMap *PatchMap::Instance() {
  if (CkpvAccess(PatchMap_instance) == 0) {
     CkpvAccess(PatchMap_instance) = new PatchMap;
  }
  return(CkpvAccess(PatchMap_instance));
}

PatchMap::PatchMap(void)
{
  nPatches = 0;
  nNodesWithPatches = 0;
  int npes = CkNumPes();
  if ( ! CkMyRank() ) {
    nPatchesOnNode = new int[npes];
    memset(nPatchesOnNode,0,npes*sizeof(int));
    patchData = NULL;
    computeIdArena = NULL;
  }
  patchBounds_a = 0;
  patchBounds_b = 0;
  patchBounds_c = 0;
  myPatch = 0;
  myHomePatch = 0;

  aDim = bDim = cDim = 0;
  aAway = bAway = cAway = 1;
  aPeriodic = bPeriodic = cPeriodic = 0;
  aMaxIndex = bMaxIndex = cMaxIndex = 0;
}

int PatchMap::sizeGrid(ScaledPosition xmin, ScaledPosition xmax,
				const Lattice &lattice, BigReal patchSize,
				double maxNumPatches, int staticAtomAssignment,
				int asplit, int bsplit, int csplit)
{
  aPeriodic = lattice.a_p();
  bPeriodic = lattice.b_p();
  cPeriodic = lattice.c_p();

  aAway = asplit;
  bAway = bsplit;
  cAway = csplit;

  int minNumPatches = 1;
  if ( aPeriodic ) minNumPatches *= aAway;
  if ( bPeriodic ) minNumPatches *= bAway;
  if ( cPeriodic ) minNumPatches *= cAway;
  if ( maxNumPatches < minNumPatches ) maxNumPatches = minNumPatches;

  if ( aPeriodic ) {
    BigReal sysDim = lattice.a_r().unit() * lattice.a();
    aDim = (int)(sysDim * aAway / patchSize);
  } else {
    BigReal sysDim = xmax.x - xmin.x;
    aDim = (int)(sysDim * aAway / patchSize);
    if ((aDim * patchSize) < (sysDim * aAway)) aDim++;
    if ( aDim < aAway + 1 ) aDim = aAway + 1;
  }

  if ( bPeriodic ) {
    BigReal sysDim = lattice.b_r().unit() * lattice.b();
    bDim = (int)(sysDim * bAway / patchSize);
  } else {
    BigReal sysDim = xmax.y - xmin.y;
    bDim = (int)(sysDim * bAway / patchSize);
    if ((bDim * patchSize) < (sysDim * bAway)) bDim++;
    if ( bDim < bAway + 1 ) bDim = bAway + 1;
  }

  if ( cPeriodic ) {
    BigReal sysDim = lattice.c_r().unit() * lattice.c();
    cDim = (int)(sysDim * cAway / patchSize);
  } else {
    BigReal sysDim = xmax.z - xmin.z;
    cDim = (int)(sysDim * cAway / patchSize);
    if ((cDim * patchSize) < (sysDim * cAway)) cDim++;
    if ( cDim < cAway + 1 ) cDim = cAway + 1;
  }

  if ( aDim < 0 || bDim < 0 || cDim < 0 ) {
    NAMD_die("Bug in PatchMap::sizeGrid - negative grid dimension.");
  }

  if ( staticAtomAssignment ) {
    if ( aPeriodic || bPeriodic || cPeriodic )
      NAMD_die("Static atom assignment is incompatible with periodic boundary conditions.");
    aDim = aAway + 1;
    bDim = bAway + 1;
    cDim = cAway + 1;
  }

  const int amin = (aPeriodic ? aAway : 1);
  const int bmin = (bPeriodic ? bAway : 1);
  const int cmin = (cPeriodic ? cAway : 1);

  // CkPrintf("searching %d-away %d-away %d-away max %d\n",aAway,bAway,cAway,(int)maxNumPatches);

  if ( aDim < amin ) aDim = amin;
  if ( bDim < bmin ) bDim = bmin;
  if ( cDim < cmin ) cDim = cmin;

  if ( maxNumPatches > aDim*bDim*cDim ) {
    maxNumPatches = aDim*bDim*cDim;
  }

  int abest = amin;
  int bbest = bmin;
  int cbest = cmin;
  int cdim = maxNumPatches;
  cdim /= aDim;  cdim /= bDim;
  if ( cdim < cmin ) cdim = cmin;
  for ( ; cdim <= cDim; ++cdim ) {
    int bdim = maxNumPatches;
    bdim /= aDim;  bdim /= cdim;
    if ( bdim < bmin ) bdim = bmin;
    for ( ; bdim <= bDim; ++bdim ) {
      int adim = maxNumPatches;
      adim /= bdim;  adim /= cdim;
      if ( adim < amin ) adim = amin;
      for ( ; adim <= aDim; ++adim ) {
        if ( adim*bdim*cdim > maxNumPatches ) break;
        // CkPrintf("testing %d * %d * %d == %d\n",adim,bdim,cdim,adim*bdim*cdim);
        if ( adim*bdim*cdim > abest*bbest*cbest ) {
          abest = adim;  bbest = bdim;  cbest = cdim;
        }
        if ( abest*bbest*cbest == maxNumPatches ) break;
      }
      if ( abest*bbest*cbest == maxNumPatches ) break;
    }
    if ( abest*bbest*cbest == maxNumPatches ) break;
  }
  aDim = abest;
  bDim = bbest;
  cDim = cbest;

  // CkPrintf("found %d * %d * %d == %d\n",aDim,bDim,cDim,aDim*bDim*cDim);
  return aDim*bDim*cDim;
}

void PatchMap::makePatches(ScaledPosition xmin, ScaledPosition xmax,
				const Lattice &lattice, BigReal patchSize,
				double maxNumPatches, int staticAtomAssignment,
				int replicaUniformPatchGrids, int lcpo,
				int asplit, int bsplit, int csplit)
{
  sizeGrid(xmin,xmax,lattice,patchSize,maxNumPatches,staticAtomAssignment,asplit,bsplit,csplit);

  if ( replicaUniformPatchGrids ) {
    int oldpcount = aDim * bDim * cDim;
    double dims[3];
    dims[0] = aDim;
    dims[1] = bDim;
    dims[2] = cDim;
    replica_min_double(dims,3);
    aDim = dims[0];
    bDim = dims[1];
    cDim = dims[2];
    int newpcount = aDim * bDim * cDim;
    if ( newpcount > oldpcount ) {
      NAMD_bug("replicaUniformPatchGrids increased patch count");
    }
    if ( newpcount < oldpcount ) {
      iout << iINFO << "PATCH GRID REDUCED TO BE UNIFORM ACROSS REPLICAS\n";
    }
  }

  iout << iINFO << "PATCH GRID IS ";
  iout << aDim;
  if ( aPeriodic ) iout << " (PERIODIC)";
  iout << " BY ";
  iout << bDim;
  if ( bPeriodic ) iout << " (PERIODIC)";
  iout << " BY ";
  iout << cDim;
  if ( cPeriodic ) iout << " (PERIODIC)";
  iout << "\n";
  iout << iINFO << "PATCH GRID IS ";
  iout << aAway << "-AWAY BY ";
  iout << bAway << "-AWAY BY ";
  iout << cAway << "-AWAY\n";
  iout << endi;

  aMaxIndex = ( ! aPeriodic || aDim == 2 ) ? 10000 : aDim;
  bMaxIndex = ( ! bPeriodic || bDim == 2 ) ? 10000 : bDim;
  cMaxIndex = ( ! cPeriodic || cDim == 2 ) ? 10000 : cDim;

  aLength = aPeriodic ? 1.0 :
      ( lcpo || aDim > aAway + 1 ? aDim * (patchSize / aAway) : xmax.x - xmin.x );
  bLength = bPeriodic ? 1.0 :
      ( lcpo || bDim > bAway + 1 ? bDim * (patchSize / bAway) : xmax.y - xmin.y );
  cLength = cPeriodic ? 1.0 :
      ( lcpo || cDim > cAway + 1 ? cDim * (patchSize / cAway) : xmax.z - xmin.z );

  aOrigin = aPeriodic ? -0.5 : 0.5 * (xmin.x + xmax.x - aLength);
  bOrigin = bPeriodic ? -0.5 : 0.5 * (xmin.y + xmax.y - bLength);
  cOrigin = cPeriodic ? -0.5 : 0.5 * (xmin.z + xmax.z - cLength);

  nPatches=aDim*bDim*cDim;
  patchData = new PatchData[nPatches];

  patchBounds_a = new BigReal[2*aDim+1];
  patchBounds_b = new BigReal[2*bDim+1];
  patchBounds_c = new BigReal[2*cDim+1];
  for ( int i=0; i<(2*aDim+1); ++i ) {
    patchBounds_a[i] = ((0.5*(double)i)/(double)aDim) * aLength + aOrigin;
  }
  for ( int i=0; i<(2*bDim+1); ++i ) {
    patchBounds_b[i] = ((0.5*(double)i)/(double)bDim) * bLength + bOrigin;
  }
  for ( int i=0; i<(2*cDim+1); ++i ) {
    patchBounds_c[i] = ((0.5*(double)i)/(double)cDim) * cLength + cOrigin;
  }

  for(int i=0; i<nPatches; ++i)
  {
    PatchData &p = patchData[i];
    p.basenode = -1;
    p.numCids = 0;
    p.aIndex = index_a(i);
    p.bIndex = index_b(i);
    p.cIndex = index_c(i);
#ifdef MEM_OPT_VERSION
    p.numAtoms = 0;
    p.numFixedAtoms = 0;
#endif
    p.numCids = 0;
    int max_computes = 30;
    p.cids = new int[max_computes];
    for ( int j = 0; j < max_computes; ++j ) p.cids[j] = -1;
    p.numCidsAllocated = max_computes;
  }

  if ( ! myPatch ) {
    myPatch = new Patch*[nPatches];
  }
  memset(myPatch,0,nPatches*sizeof(Patch*));
  if ( ! myHomePatch ) {
    myHomePatch = new HomePatch*[nPatches];
  }
  memset(myHomePatch,0,nPatches*sizeof(HomePatch*));
}

void PatchMap::checkMap(void)
{
  int patchCount=0;
  for (int i=0; i<nPatches; i++) {
    if (myPatch[i]) {
      patchCount++;
      if ( myPatch[i]->getPatchID() != i) {
	DebugM(4, "patchID("<<myPatch[i]->getPatchID()
	  <<") != patchID(" 
	  <<i<<")\n");
      }
    }
  }
  DebugM(4, "Patch Count = " <<patchCount<<"\n");
}
  

PatchMap::~PatchMap(void)
{
  if ( ! CkMyRank() ) {
    if (patchData && ! computeIdArena ) {
      for (int i=0; i<nPatches; i++) {
        delete [] patchData[i].cids;
      }
    }
    delete [] patchData;
    delete [] nPatchesOnNode;
    delete computeIdArena;
  }
  delete [] patchBounds_a;
  delete [] patchBounds_b;
  delete [] patchBounds_c;
  delete [] myPatch;
  delete [] myHomePatch;
}

#undef PACK
#define PACK(type,data) { memcpy(b, &data,sizeof(type)); b += sizeof(type); }
#define PACKN(type,data,cnt) { memcpy(b, data,(cnt)*sizeof(type)); b += (cnt)*sizeof(type); }

int PatchMap::packSize(void)
{
  int i, size = 0;
  size += 14 * sizeof(int) + 6 * sizeof(BigReal);
  size += (2*(aDim+bDim+cDim)+3) * sizeof(BigReal);
  size += CkNumPes() * sizeof(int);
  for(i=0;i<nPatches;++i)
  {
    size += sizeof(PatchData);
    size += patchData[i].numCids * sizeof(ComputeID);
  }
  return size;
}

void PatchMap::pack (char *buffer, int size)
{
  DebugM(4,"Packing PatchMap on node " << CkMyPe() << std::endl);
  int i,j;

  // fill in the data
  char *b = buffer;
  PACK(int,nPatches);
  DebugM(3,"nPatches = " << nPatches << std::endl);
  PACK(int,aDim); PACK(int,bDim); PACK(int,cDim);
  PACK(int,aAway); PACK(int,bAway); PACK(int,cAway);
  PACK(int,aPeriodic); PACK(int,bPeriodic); PACK(int,cPeriodic);
  PACK(int,aMaxIndex); PACK(int,bMaxIndex); PACK(int,cMaxIndex);
  PACK(BigReal,aOrigin); PACK(BigReal,bOrigin); PACK(BigReal,cOrigin);
  PACK(BigReal,aLength); PACK(BigReal,bLength); PACK(BigReal,cLength);
  PACK(int,nNodesWithPatches);
  PACKN(BigReal,patchBounds_a,2*aDim+1);
  PACKN(BigReal,patchBounds_b,2*bDim+1);
  PACKN(BigReal,patchBounds_c,2*cDim+1);
  PACKN(int,nPatchesOnNode,CkNumPes());
  PACKN(PatchData,patchData,nPatches);
  for(i=0;i<nPatches;++i) {
    DebugM(3,"Packing Patch " << i << " is on node " << patchData[i].node << 
	" with " << patchData[i].numCids << " cids.\n");
    PACKN(ComputeID,patchData[i].cids,patchData[i].numCids);
  }
  if ( buffer + size != b ) {
    NAMD_bug("PatchMap::pack does not match PatchMap::packSize");
  }
  //DebugM(3,buffer + size - b << " == 0 ?" << std::endl);
}

#undef UNPACK
#define UNPACK(type,data) { memcpy(&data, b, sizeof(type)); b += sizeof(type); }
#define UNPACKN(type,data,cnt) { memcpy(data, b, (cnt)*sizeof(type)); b += (cnt)*sizeof(type); }
#define SKIPN(type,cnt) { b += (cnt)*sizeof(type); }

void PatchMap::unpack (char *ptr)
{
  DebugM(4,"Unpacking PatchMap on node " << CkMyPe() << std::endl);

  int i,j;
  char *b = (char*)ptr;
  {
    // defeat some over-zealous compilers
    int nPatches_tmp;
    UNPACK(int,nPatches_tmp);
    nPatches = nPatches_tmp;
  }
  DebugM(3,"nPatches = " << nPatches << std::endl);

  if ( ! myPatch ) {
    myPatch = new Patch*[nPatches];
  }
  memset(myPatch,0,nPatches*sizeof(Patch*));
  if ( ! myHomePatch ) {
    myHomePatch = new HomePatch*[nPatches];
  }
  memset(myHomePatch,0,nPatches*sizeof(HomePatch*));

  UNPACK(int,aDim); UNPACK(int,bDim); UNPACK(int,cDim);
  UNPACK(int,aAway); UNPACK(int,bAway); UNPACK(int,cAway);
  UNPACK(int,aPeriodic); UNPACK(int,bPeriodic); UNPACK(int,cPeriodic);
  UNPACK(int,aMaxIndex); UNPACK(int,bMaxIndex); UNPACK(int,cMaxIndex);
  UNPACK(BigReal,aOrigin); UNPACK(BigReal,bOrigin); UNPACK(BigReal,cOrigin);
  UNPACK(BigReal,aLength); UNPACK(BigReal,bLength); UNPACK(BigReal,cLength);
  UNPACK(int,nNodesWithPatches);


//  CkPrintf("[%d] has bounds a %d b %d c %d npatches %d mem %d\n",CkMyPe(),aDim, bDim, cDim, nPatches, memusage_MB() );

  if ( ! patchBounds_a ) patchBounds_a = new BigReal[2*aDim+1];
  if ( ! patchBounds_b ) patchBounds_b = new BigReal[2*bDim+1];
  if ( ! patchBounds_c ) patchBounds_c = new BigReal[2*cDim+1];
  UNPACKN(BigReal,patchBounds_a,2*aDim+1);
  UNPACKN(BigReal,patchBounds_b,2*bDim+1);
  UNPACKN(BigReal,patchBounds_c,2*cDim+1);
 
  if ( CkMyRank() ) return;

  UNPACKN(int,nPatchesOnNode,CkNumPes());

  if ( ! patchData ) patchData = new PatchData[nPatches];
  else if ( ! computeIdArena ) {
    for(i=0;i<nPatches;++i) {
      delete [] patchData[i].cids;
    }
  }
  UNPACKN(PatchData,patchData,nPatches);

  delete computeIdArena;
  computeIdArena = new ObjectArena<ComputeID>;
  computeIdArena->setBlockSize(1024);

  for(i=0;i<nPatches;++i) {
    DebugM(3,"Unpacking Patch " << i << " is on node " << patchData[i].node << 
	" with " << patchData[i].numCids << " cids.\n");
    patchData[i].cids = computeIdArena->getNewArray(patchData[i].numCids);
    patchData[i].numCidsAllocated = patchData[i].numCids;
    UNPACKN(ComputeID,patchData[i].cids,patchData[i].numCids);
  }
}

//----------------------------------------------------------------------
int PatchMap::numHomePatches(void)
{
  return CkpvAccess(PatchMap_patchMgr)->homePatches.size();
}

//----------------------------------------------------------------------
HomePatchList *PatchMap::homePatchList() {
  return &(CkpvAccess(PatchMap_patchMgr)->homePatches);
}

//----------------------------------------------------------------------
void PatchMap::homePatchIDList(PatchIDList &pids) {
  pids.resize(0);
  int i;
  for ( i=0; i<nPatches; ++i ) {
    if ( patchData[i].node == CkMyPe() ) {
      pids.add(i);
    }
  }
}

//----------------------------------------------------------------------
void PatchMap::basePatchIDList(int pe, PatchIDList &pids) {
  pids.resize(0);
  int i;
  for ( i=0; i<nPatches; ++i ) {
    if ( patchData[i].basenode == pe ) {
      pids.add(i);
    }
  }
}

//----------------------------------------------------------------------
void PatchMap::assignNode(PatchID pid, NodeID node) {
  patchData[pid].node=node;
  if ( nPatchesOnNode[node] == 0 ) nNodesWithPatches += 1;
  nPatchesOnNode[node] += 1;
}

//----------------------------------------------------------------------
void PatchMap::assignBaseNode(PatchID pid, NodeID node) {
  patchData[pid].basenode=node;
}

void PatchMap::assignBaseNode(PatchID pid) {
  
  int i = 1;

  NodeID node = patchData[pid].node;

  if ( CkNumPes() > 2*nPatches+1 ) {

    int newnode =  ( CkNumPes() + node - 1 ) % CkNumPes();    
    bool success = 0;

    while ( i < CkNumPes() && !success) {
      if ( nPatchesOnNode[newnode] == 0 )
	success = 1;

      //we know till pid, we have assigned all base nodes
      for (int count = 0; count < pid; count ++)
	if (patchData[count].basenode > 0 && patchData[count].basenode == newnode) {
	  success = 0;
	  break;
	}
	  
      //no patch or a patche's base node on this newnode. this is a good node
      if (success) break;

      newnode = ( CkNumPes() + node - i - 1 ) % CkNumPes();
      i ++;
    }
    patchData[pid].basenode = newnode;

  } else {
    patchData[pid].basenode=node;
  }
}

//----------------------------------------------------------------------
void PatchMap::newCid(int pid, int cid)
{
  if (patchData[pid].numCids >= patchData[pid].numCidsAllocated)
  { // allocate more
//    NAMD_die("PatchMap::newCid - not enough compute ID's allocated.");
    ComputeID *old = patchData[pid].cids;
    patchData[pid].numCidsAllocated += 10;
    patchData[pid].cids = new int[patchData[pid].numCidsAllocated];
    int i;
    for (i=0; i<patchData[pid].numCids; i++) 
    	patchData[pid].cids[i] = old[i];
    for (i=patchData[pid].numCids; i<patchData[pid].numCidsAllocated; i++) 
	patchData[pid].cids[i] = -1;
    delete [] old;
  }
  patchData[pid].cids[patchData[pid].numCids]=cid;
  patchData[pid].numCids++;
}

//----------------------------------------------------------------------
int PatchMap::oneAwayNeighbors(int pid, PatchID *neighbor_ids)
{
  int xi, yi, zi;
  int xinc, yinc, zinc;
  int n=0;

  for(zinc=-1;zinc<=1;zinc++)
  {
    zi = patchData[pid].cIndex + zinc;
    if ((zi < 0) || (zi >= cDim))
      if ( ! cPeriodic ) continue;
    for(yinc=-1;yinc<=1;yinc++)
    {
      yi = patchData[pid].bIndex + yinc;
      if ((yi < 0) || (yi >= bDim))
	if ( ! bPeriodic ) continue;
      for(xinc=-1;xinc<=1;xinc++)
      {
	if ((xinc==0) && (yinc==0) && (zinc==0))
	  continue;

	xi = patchData[pid].aIndex + xinc;
	if ((xi < 0) || (xi >= aDim))
	  if ( ! aPeriodic ) continue;

	if (neighbor_ids)
	  neighbor_ids[n]=this->pid(xi,yi,zi);
#if 0
	if ( transform_ids )
	{
	  int xt = 0; if ( xi < 0 ) xt = -1; if ( xi >= aDim ) xt = 1;
	  int yt = 0; if ( yi < 0 ) yt = -1; if ( yi >= bDim ) yt = 1;
	  int zt = 0; if ( zi < 0 ) zt = -1; if ( zi >= cDim ) zt = 1;
	  transform_ids[n] = Lattice::index(xt,yt,zt);
	}
#endif
	n++;
      }
    }
  }
  DebugM(3,"Patch " << pid << " has " << n << " first neighbors.\n");
  return n;
}


//----------------------------------------------------------------------
// Only returns half of neighbors!
int PatchMap::oneOrTwoAwayNeighbors(int pid, PatchID *neighbor_ids, PatchID *downstream_ids, int *transform_ids)
{
  int xi, yi, zi;
  int xinc, yinc, zinc;
  int n=0;
  const int xs = patchData[pid].aIndex;
  const int ys = patchData[pid].bIndex;
  const int zs = patchData[pid].cIndex;

  for(zinc=0;zinc<=cAway;zinc++)
  {
    zi = zs + zinc;
    if ((zi < 0) || (zi >= cDim))
      if ( ! cPeriodic ) continue;
    for(yinc=(zinc>0 ? -bAway : 0);yinc<=bAway;yinc++)
    {
      yi = ys + yinc;
      if ((yi < 0) || (yi >= bDim))
	if ( ! bPeriodic ) continue;
      for(xinc=((zinc>0 || yinc>0) ? -aAway : 0);xinc<=aAway;xinc++)
      {
	if ((xinc==0) && (yinc==0) && (zinc==0))
	  continue;

	xi = xs + xinc;
	if ((xi < 0) || (xi >= aDim))
	  if ( ! aPeriodic ) continue;

	neighbor_ids[n] = this->pid(xi,yi,zi);
	if ( transform_ids )
	{
	  int xt = 0; if ( xi < 0 ) xt = -1; if ( xi >= aDim ) xt = 1;
	  int yt = 0; if ( yi < 0 ) yt = -1; if ( yi >= bDim ) yt = 1;
	  int zt = 0; if ( zi < 0 ) zt = -1; if ( zi >= cDim ) zt = 1;
	  transform_ids[n] = Lattice::index(xt,yt,zt);
	}
	if ( downstream_ids )
	{
	  int xd = ( xi < xs ? xi : xs );
	  int yd = ( yi < ys ? yi : ys );
	  int zd = ( zi < zs ? zi : zs );
	  downstream_ids[n] = this->pid(xd,yd,zd);
	}
	n++;
      }
    }
  }
  DebugM(3,"Patch " << pid << " has " << n << " second neighbors.\n");
  return n;
}

//----------------------------------------------------------------------
// Return all patches in corresponding octet (2x2x2 block of patches)
// regardless of periodic boundary conditions
// Used for LCPO Computes
int PatchMap::getPatchesInOctet(int pid, PatchID *pids, int *transform_ids)
{
  int xi, yi, zi;
  int xinc, yinc, zinc;
  int n=0;
  const int xs = patchData[pid].aIndex;
  const int ys = patchData[pid].bIndex;
  const int zs = patchData[pid].cIndex;

  for(zinc=0; zinc<2; zinc++) {
    zi = zs + zinc;
    for(yinc=0; yinc<2; yinc++) {
      yi = ys + yinc;
      for(xinc=0; xinc<2; xinc++) {
	      xi = xs + xinc;
        int aIndex = MODULO(xi,aDim);
        int bIndex = MODULO(yi,bDim);
        int cIndex = MODULO(zi,cDim);
        pids[n] = ((cIndex*bDim)+bIndex)*aDim + aIndex;
      	if ( transform_ids ) {
	        int xt = 0; if ( xi < 0 ) xt = -1; if ( xi >= aDim ) xt = 1;
	        int yt = 0; if ( yi < 0 ) yt = -1; if ( yi >= bDim ) yt = 1;
	        int zt = 0; if ( zi < 0 ) zt = -1; if ( zi >= cDim ) zt = 1;
	        transform_ids[n] = Lattice::index(xt,yt,zt);
	      }
	      n++;
      } // for x
    } // for y
  } // for z
  DebugM(3,"Patch " << pid << " has " << n << " second neighbors.\n");
  return n;
}


//----------------------------------------------------------------------
int PatchMap::upstreamNeighbors(int pid, PatchID *neighbor_ids)
{
  int xi, yi, zi;
  int xinc, yinc, zinc;
  int n=0;

  for(zinc=0;zinc<=1;zinc++)
  {
    zi = patchData[pid].cIndex + zinc;
    if ((zi < 0) || (zi >= cDim))
      if ( ! cPeriodic ) continue;
    for(yinc=0;yinc<=1;yinc++)
    {
      yi = patchData[pid].bIndex + yinc;
      if ((yi < 0) || (yi >= bDim))
	if ( ! bPeriodic ) continue;
      for(xinc=0;xinc<=1;xinc++)
      {
	if ((xinc==0) && (yinc==0) && (zinc==0))
	  continue;

	xi = patchData[pid].aIndex + xinc;
	if ((xi < 0) || (xi >= aDim))
	  if ( ! aPeriodic ) continue;

	if (neighbor_ids)
	  neighbor_ids[n]=this->pid(xi,yi,zi);
#if 0
	if ( transform_ids )
	{
	  int xt = 0; if ( xi < 0 ) xt = -1; if ( xi >= aDim ) xt = 1;
	  int yt = 0; if ( yi < 0 ) yt = -1; if ( yi >= bDim ) yt = 1;
	  int zt = 0; if ( zi < 0 ) zt = -1; if ( zi >= cDim ) zt = 1;
	  transform_ids[n] = Lattice::index(xt,yt,zt);
	}
#endif
	n++;
      }
    }
  }
  DebugM(3,"Patch " << pid << " has " << n << " upstream neighbors.\n");
  return n;
}

//----------------------------------------------------------------------
int PatchMap::downstreamNeighbors(int pid, PatchID *neighbor_ids)
{
  int xi, yi, zi;
  int xinc, yinc, zinc;
  int n=0;

  for(zinc=-1;zinc<=0;zinc++)
  {
    zi = patchData[pid].cIndex + zinc;
    if ((zi < 0) || (zi >= cDim))
      if ( ! cPeriodic ) continue;
    for(yinc=-1;yinc<=0;yinc++)
    {
      yi = patchData[pid].bIndex + yinc;
      if ((yi < 0) || (yi >= bDim))
	if ( ! bPeriodic ) continue;
      for(xinc=-1;xinc<=0;xinc++)
      {
	if ((xinc==0) && (yinc==0) && (zinc==0))
	  continue;

	xi = patchData[pid].aIndex + xinc;
	if ((xi < 0) || (xi >= aDim))
	  if ( ! aPeriodic ) continue;

	if (neighbor_ids)
	  neighbor_ids[n]=this->pid(xi,yi,zi);
#if 0
	if ( transform_ids )
	{
	  int xt = 0; if ( xi < 0 ) xt = -1; if ( xi >= aDim ) xt = 1;
	  int yt = 0; if ( yi < 0 ) yt = -1; if ( yi >= bDim ) yt = 1;
	  int zt = 0; if ( zi < 0 ) zt = -1; if ( zi >= cDim ) zt = 1;
	  transform_ids[n] = Lattice::index(xt,yt,zt);
	}
#endif
	n++;
      }
    }
  }
  DebugM(3,"Patch " << pid << " has " << n << " upstream neighbors.\n");
  return n;
}

//----------------------------------------------------------------------
void PatchMap::printPatchMap(void)
{
  CkPrintf("---------------------------------------");
  CkPrintf("---------------------------------------\n");

  CkPrintf("nPatches = %d\n",nPatches);
  for(int i=0;i<nPatches;i++)
  {
    CkPrintf("Patch %d:\n",i);
    CkPrintf("  node = %d\n",patchData[i].node);
    CkPrintf("  xi,yi,zi = %d, %d, %d\n",
	    patchData[i].aIndex,patchData[i].bIndex,patchData[i].cIndex);
    CkPrintf("  numCids = %d\n",patchData[i].numCids);
    CkPrintf("  numCidsAllocated = %d\n",patchData[i].numCidsAllocated);
    for(int j=0; j < patchData[i].numCids; j++)
    {
      CkPrintf(" %10d ",patchData[i].cids[j]);
      if (!((j+1) % 6))
	CkPrintf("\n");
    }
    CkPrintf("\n---------------------------------------");
    CkPrintf("---------------------------------------\n");
  }

}

//----------------------------------------------------------------------
void PatchMap::registerPatch(PatchID pid, HomePatch *pptr) {
  registerPatch(pid,(Patch*)pptr);
  if (myHomePatch[pid] != 0) {
    iout << iPE << iERRORF 
      << "homePatchID("<<pid<<") is being re-registered!\n" << endi;
  }
  myHomePatch[pid] = pptr;
}

//----------------------------------------------------------------------
void PatchMap::unregisterPatch(PatchID pid, HomePatch *pptr) {
  unregisterPatch(pid,(Patch*)pptr);
  if (pptr == myHomePatch[pid]) {
      DebugM(4, "UnregisterHomePatch("<<pid<<") at " << pptr << "\n");
      myHomePatch[pid] = NULL;
  }
}

//----------------------------------------------------------------------
void PatchMap::registerPatch(PatchID pid, Patch *pptr)
{
  if (myPatch[pid] != 0) {
    iout << iPE << iERRORF 
      << "patchID("<<pid<<") is being re-registered!\n" << endi;
  }
  myPatch[pid] = pptr;
}

//----------------------------------------------------------------------
void PatchMap::unregisterPatch(PatchID pid, Patch *pptr)
{
  if (pptr == myPatch[pid]) {
      DebugM(4, "UnregisterPatch("<<pid<<") at " << pptr << "\n");
      myPatch[pid] = NULL;
  }
}

