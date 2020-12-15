/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/
#include "ComputeLCPO.h"
#include "ReductionMgr.h"
#include "PatchMap.h"
#include "ComputeMgr.h"
#include "Molecule.h"
#include "Node.h"
#include "SimParameters.h"
#include "Debug.h"
#include "WorkDistrib.decl.h"
#include "Node.h"
#include "ComputeLCPO.h"
#include "Priorities.h"
#include "PatchMap.inl"
#include "Patch.h"
#include "ComputeMap.h"
#include "LdbCoordinator.h"
#include "common.h"
#include "time.h"

//#define COUNT_FLOPS

#ifdef COUNT_FLOPS
#define FLOPS(X) flops += X;
#else
#define FLOPS(X)
#endif

//#define MIN_DEBUG_LEVEL 4
// #define DEBUGM

ComputeLCPO::ComputeLCPO(ComputeID c, PatchID p[], int t[], 
		ComputeNonbondedWorkArrays* _workArrays,
		int minPartition, int maxPartition, int numPartitions, int numPatches)
  : Compute(c), workArrays(_workArrays),
    minPart(minPartition), maxPart(maxPartition),
    strideIg(numPartitions), numParts(numPartitions),
    maxAtomRadius(1.9+1.4)
  {

  reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);

  setNumPatches(8);
  SimParameters *simParams = Node::Object()->simParameters;
  surfTen = simParams->surface_tension;

  for (int i=0; i<getNumPatches(); i++) {
    patchID[i] = p[i];
    trans[i] = t[i];
    patch[i] = NULL;
    positionBox[i] = NULL;
    forceBox[i] = NULL;
    lcpoTypeBox[i] = NULL;
  } // for all patches
} // constructor

ComputeLCPO::~ComputeLCPO() {
  DebugM(4, "~ComputeLCPO("<<cid<<") numAtoms("<<patchID[0]<<") = " 
    << numAtoms[0] 
    << " numAtoms("<<patchID[1]<<") = " << numAtoms[1] << "\n" );
  DebugM(4, "~ComputeLCPO("<<cid<<") addr("<<patchID[0]<<") = " 
    << PatchMap::Object()->patch(patchID[0]) << " addr("<<patchID[1]<<") = "
    << PatchMap::Object()->patch(patchID[1]) << "\n");

  for (int i=0; i<getNumPatches(); i++) {
    if (positionBox[i] != NULL) {
      PatchMap::Object()->patch(patchID[i])->unregisterPositionPickup(this,
	 &positionBox[i]);
    }
    if (forceBox[i] != NULL) {
      PatchMap::Object()->patch(patchID[i])->unregisterForceDeposit(this,
		&forceBox[i]);
    }
    if (lcpoTypeBox[i] != NULL) {
      PatchMap::Object()->patch(patchID[i])->unregisterLcpoTypePickup(this,
		&lcpoTypeBox[i]);
    }
  }
  delete reduction;
} // destructor

void ComputeLCPO::initialize() {
  Compute::initialize();
    // How can we tell if BoxOwner has packed up and left?  Need a mechanism
    // to handle this or do we assume the Boxes have been dumped?
    PatchMap *patchMap = PatchMap::Object();
    //Check out Boxes
    for (int i=0; i<8; i++) {
      //invalid patch so don't even checkout boxes
	    if (positionBox[i] == NULL) { // We have yet to get boxes
        patch[i] = PatchMap::Object()->patch(patchID[i]);
	      if (!(patch[i] = PatchMap::Object()->patch(patchID[i]))) {
	        DebugM(5,"invalid patch(" << patchID[i] 
		      << ")  pointer!\n");
	      }
        positionBox[i] = patch[i]->registerPositionPickup(this);
        forceBox[i] = patch[i]->registerForceDeposit(this);
        lcpoTypeBox[i] = patch[i]->registerLcpoTypePickup(this);
        // will need to open a box full of lcpo parameters
	    }
      numAtoms[i] = patch[i]->getNumAtoms();
    } // for all patches

  // set priority
  basePriority = PATCH_PRIORITY(patchID[0]) + PROXY_RESULTS_PRIORITY;

  //get bounds of inner rectangular prism in octet
  bounds[0][0] = 0.5*(patchMap->min_a(patchID[0])+patchMap->max_a(patchID[0]));
  bounds[1][0] = 0.5*(patchMap->min_b(patchID[0])+patchMap->max_b(patchID[0]));
  bounds[2][0] = 0.5*(patchMap->min_c(patchID[0])+patchMap->max_c(patchID[0]));
  bounds[0][1] = 0.5*(patchMap->min_a(patchID[7])+patchMap->max_a(patchID[7]));
  bounds[1][1] = 0.5*(patchMap->min_b(patchID[7])+patchMap->max_b(patchID[7]));
  bounds[2][1] = 0.5*(patchMap->min_c(patchID[7])+patchMap->max_c(patchID[7]));

  //if only 1 patch in a dimenion, invalidate those patches
  int gsa = patchMap->gridsize_a();
  int gsb = patchMap->gridsize_b();
  int gsc = patchMap->gridsize_c();
  invalidPatch[0] = 0;
  invalidPatch[1] = 0;
  invalidPatch[2] = 0;
  invalidPatch[3] = 0;
  invalidPatch[4] = 0;
  invalidPatch[5] = 0;
  invalidPatch[6] = 0;
  invalidPatch[7] = 0;

  if (gsa==1) {
    //CkPrintf("ONLY 1 PATCH in A DIMENSION!\n");
    invalidPatch[1] = 1;
    invalidPatch[3] = 1;
    invalidPatch[5] = 1;
    invalidPatch[7] = 1;
  }
  if (gsb==1) {
    //CkPrintf("ONLY 1 PATCH in B DIMENSION!\n");
    invalidPatch[2] = 1;
    invalidPatch[3] = 1;
    invalidPatch[6] = 1;
    invalidPatch[7] = 1;
  }
  if (gsc==1) {
    //CkPrintf("ONLY 1 PATCH in C DIMENSION!\n");
    invalidPatch[4] = 1;
    invalidPatch[5] = 1;
    invalidPatch[6] = 1;
    invalidPatch[7] = 1;
  }
  //relative a,b,c index for 8 patches in ComputeLCPO
  int idx[8][3] = {
    { 0, 0, 0},
    { 1, 0, 0},
    { 0, 1, 0},
    { 1, 1, 0},
    { 0, 0, 1},
    { 1, 0, 1},
    { 0, 1, 1},
    { 1, 1, 1}    };
/*
  int i_a = patchMap->index_a(patchID[0]);
  int i_b = patchMap->index_b(patchID[0]);
  int i_c = patchMap->index_c(patchID[0]);
  CkPrintf("VALID[%d,%d,%d]=\n",i_a,i_b,i_c);
*/
  for (int pI = 0; pI < 8; pI++) {
    int iia = patchMap->index_a(patchID[pI]);
    int iib = patchMap->index_b(patchID[pI]);
    int iic = patchMap->index_c(patchID[pI]);
    for (int pJ = 0; pJ < 8; pJ++) {
      int jia = patchMap->index_a(patchID[pJ]);
      int jib = patchMap->index_b(patchID[pJ]);
      int jic = patchMap->index_c(patchID[pJ]);
      if (  ( gsa==1 && (jia>iia) != (idx[pJ][0]>idx[pI][0]) ) ||
            ( gsb==1 && (jib>iib) != (idx[pJ][1]>idx[pI][1]) ) ||
            ( gsc==1 && (jic>iic) != (idx[pJ][2]>idx[pI][2]) ) ||
            ( invalidPatch[pI] ) ||
            ( invalidPatch[pJ] )   )
        valid[pI][pJ] = 0;
      else
        valid[pI][pJ] = 1;
      //CkPrintf("%d ",valid[pI][pJ]);
    }
    //CkPrintf("\n");
  }
  //CkPrintf("\n");

} // initialize

void ComputeLCPO::atomUpdate() {
  for (int i=0; i<8; i++) {
	  numAtoms[i] = patch[i]->getNumAtoms();
  }
}

//---------------------------------------------------------------------
// doWork
//---------------------------------------------------------------------
void ComputeLCPO::doWork() {
  LdbCoordinator::Object()->startWork(ldObjHandle);
  for (int i=0; i<8; i++) {
    pos[i] = positionBox[i]->open();
    force[i] = forceBox[i]->open();
    posExt[i] = patch[i]->getCompAtomExtInfo();
    lcpoType[i] = lcpoTypeBox[i]->open();
  }

  doForce();

 // Inform load balancer
  LdbCoordinator::Object()->endWork(ldObjHandle);

  // Close up boxes
  for (int i=0; i<getNumPatches(); i++) {
    positionBox[i]->close(&pos[i]);
    forceBox[i]->close(&force[i]);
    lcpoTypeBox[i]->close(&lcpoType[i]);
  }
} // doWork


int ComputeLCPO::noWork() {

  if ( patch[0]->flags.doNonbonded) {
    return 0;  // work to do, enqueue as usual
  } else {

    // skip all boxes
    for (int i=0; i<8; i++) {
      positionBox[i]->skip();
      forceBox[i]->skip();
      lcpoTypeBox[i]->skip();
    }

    reduction->item(REDUCTION_COMPUTE_CHECKSUM) += 1.;
    reduction->submit();
    LdbCoordinator::Object()->skipWork(ldObjHandle);

    return 1;  // no work to do, do not enqueue
  }
  return 0;
} // noWork

// 1 - yes in bounds, 0 - not in bounds
// this does uniquely assign atoms to ComputeLCPO octets
int ComputeLCPO::isInBounds(Real x, Real y, Real z ) {

  //check x dimension
  if ( bounds[0][0] < bounds[0][1] ) { // internal
    if (x < bounds[0][0] || x >= bounds[0][1] )
      return 0;
  } else { // edge
    if (x < bounds[0][0] && x >= bounds[0][1] )
      return 0;
  }

  //check y dimension
  if ( bounds[1][0] < bounds[1][1] ) { // internal 
    if (y < bounds[1][0] || y >= bounds[1][1] )
      return 0;
  } else { // edge
    if (y < bounds[1][0] && y >= bounds[1][1] )
      return 0;
  }

  //check z dimension
  if ( bounds[2][0] < bounds[2][1] ) { // internal
    if (z < bounds[2][0] || z >= bounds[2][1] )
      return 0;
  } else { // edge
    if (z < bounds[2][0] && z >= bounds[2][1] )
      return 0;
  }

  return 1;
} // isInBounds

inline BigReal calcOverlap( BigReal r, Real ri, Real rj ) {
  return PI*ri*(2*ri-r-(ri*ri-rj*rj)/r);
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//// doForce
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void ComputeLCPO::doForce() {
  //CkPrintf("ComputeLCPO::doForce\n");
  step = patch[0]->flags.sequence;

  Real probeRadius = 1.4f;
  Real cutMargin = 0.f;//regenerating pairlists every step

  Position ngir, ngjr, ngkr;
  Real ri, rj, rk;
  BigReal dxij, dyij, dzij, r2ij;
  BigReal dxik, dyik, dzik, r2ik;
  BigReal dxjk, dyjk, dzjk, r2jk;

#ifdef COUNT_FLOPS
 int flops = 0;
#endif

//////////////////////////////////////////////////
// Build Pairlists
//////////////////////////////////////////////////
//generate pairlists every step since contain coordinates
if ( true ) {
  double t_start = 1.0*clock()/CLOCKS_PER_SEC;

  inAtomsPl.reset();
  lcpoNeighborList.reset();
  FLOPS(8);
  cut2 = 2*maxAtomRadius+cutMargin; cut2 *= cut2;
  maxAtomRadius = 0;
  //find in-bounds atoms in each patch
  for (int pI = 0; pI < 8; pI++) {
    if (invalidPatch[pI]) continue;
    if (numAtoms[pI] == 0) continue;

    int minIg = 0;
    for (int s = 0; s < minPart; s++) {
      minIg += pos[pI][minIg].hydrogenGroupSize;
      FLOPS(1)
    }
    strideIg = numParts;//stride through partitions
    plint *inAtoms = inAtomsPl.newlist(numAtoms[pI]);
    int numAtomsInBounds = 0;

    //iterate over heavy atoms only
    for ( int ngi = minIg; ngi < numAtoms[pI]; /* ngi */) {
      ngir = pos[pI][ngi].position;
      if ( isInBounds(ngir.x, ngir.y, ngir.z) && lcpoType[pI][ngi] > 0 ) {
        inAtoms[numAtomsInBounds++] = ngi;
        ri = probeRadius+lcpoParams[ lcpoType[pI][ngi] ][0];
        maxAtomRadius = (ri > maxAtomRadius) ? ri : maxAtomRadius;
        FLOPS(1);

        int maxAtoms = 0;
        for (int pJ = 0; pJ < 8; pJ++) {
          if (numAtoms[pJ] > 0) {
            maxAtoms += numAtoms[pJ];
          }
        }
        LCPOAtom *lcpoNeighbors = lcpoNeighborList.newlist(maxAtoms);
        int numLcpoNeighbors = 0;

        //find pairs of this inAtom from all 8 patches
        for (int pJ = 0; pJ < 8; pJ++) {
          if (invalidPatch[pJ]) continue;
          if (!valid[pI][pJ]) continue;

          // j atom pairs
          for ( int ngj = 0; ngj < numAtoms[pJ]; /* ngj */) {
            FLOPS(1)
            ngjr = pos[pJ][ngj].position;
            dxij = ngir.x - ngjr.x;
            dyij = ngir.y - ngjr.y;
            dzij = ngir.z - ngjr.z;

            // i-j coarse check if too far apart
            r2ij = dxij*dxij + dyij*dyij + dzij*dzij;
            FLOPS(8)
            if (r2ij < cut2 && r2ij > 0.01) {

              // i-j precise check if too far apart
              rj = probeRadius+lcpoParams[ lcpoType[pJ][ngj] ][0];
              FLOPS(5)
              BigReal rirjcutMargin2 = ri+rj+cutMargin;
              rirjcutMargin2 *= rirjcutMargin2;
              if (r2ij < rirjcutMargin2 && r2ij > 0.0001 &&
                  lcpoType[pJ][ngj] > 0) {
                lcpoNeighbors[numLcpoNeighbors].x = ngjr.x;
                lcpoNeighbors[numLcpoNeighbors].y = ngjr.y;
                lcpoNeighbors[numLcpoNeighbors].z = ngjr.z;
                lcpoNeighbors[numLcpoNeighbors].r = rj;
                lcpoNeighbors[numLcpoNeighbors].f =
                  &force[pJ]->f[Results::nbond][ngj];
                numLcpoNeighbors++;
                FLOPS(2)
                maxAtomRadius = (rj > maxAtomRadius) ? rj : maxAtomRadius;
              } // precise cutoff
            } // coarse cutoff
            //jump to next nonbonded group
            ngj += pos[pJ][ngj].hydrogenGroupSize;
            FLOPS(1)
          } // for j atoms
        } // for patches J
        lcpoNeighborList.newsize(numLcpoNeighbors);
      } // in bounds
      //jump to next nonbonded group for round-robin
      for (int s = 0; s < strideIg; s++) {
        ngi += pos[pI][ngi].hydrogenGroupSize;
        FLOPS(1)
      }
    } // for i atoms
    inAtomsPl.newsize(numAtomsInBounds);
  } // for patches I
#ifdef COUNT_FLOPS
  double t_stop = 1.0*clock()/CLOCKS_PER_SEC;
  CkPrintf("LCPO_TIME_P %7.3f Gflops %9d @ %f\n", flops*1e-9/(t_stop-t_start),flops,(t_stop-t_start));
#endif
}
#ifdef COUNT_FLOPS
  double t_start = 1.0*clock()/CLOCKS_PER_SEC;
  flops = 0;
#endif

  //reset pairlists
  inAtomsPl.reset();
  lcpoNeighborList.reset();
  cut2 = maxAtomRadius*2; cut2 *= cut2;

  //init values
  BigReal totalSurfaceArea = 0;

//////////////////////////////////////////////////
//////////////////////////////////////////////////
////
////   Perform LCPO Calculation
////
//////////////////////////////////////////////////
//////////////////////////////////////////////////
  //for each patch in octet
  for (int pI = 0; pI < 8; pI++) {
    if (invalidPatch[pI]) continue;
    if (numAtoms[pI] == 0) continue;
    plint *inAtoms;
    int numInAtoms;
    inAtomsPl.nextlist( &inAtoms, &numInAtoms );
    //for each inAtom in each patch
    for (int i = 0; i < numInAtoms; i++) {
      int iIndex = inAtoms[i];
      int idi = posExt[pI][iIndex].id;
      Real xi = pos[pI][iIndex].position.x;
      Real yi = pos[pI][iIndex].position.y;
      Real zi = pos[pI][iIndex].position.z;
      const Real *lcpoParamI = lcpoParams[ lcpoType[pI][iIndex] ];
      ri = probeRadius+lcpoParamI[0];
      FLOPS(1)

      Real P1 = lcpoParamI[1];
      Real P2 = lcpoParamI[2];
      Real P3 = lcpoParamI[3];
      Real P4 = lcpoParamI[4];

//////////////////////////////////////////////////
// S1
//////////////////////////////////////////////////
      BigReal S1 = 4.0*PI*ri*ri; // a
      FLOPS(3)

      //for surface area calculation
      BigReal AijSum = 0; // b
      BigReal AjkSum = 0; // c
      BigReal AjkjSum = 0; // d'
      BigReal AijAjkSum = 0; // d

      //for force calculation
      BigReal dAijdrijdxiSum    = 0.0;
      BigReal dAijdrijdyiSum    = 0.0;
      BigReal dAijdrijdziSum    = 0.0;
      BigReal dAijdrijdxiAjkSum = 0.0;
      BigReal dAijdrijdyiAjkSum = 0.0;
      BigReal dAijdrijdziAjkSum = 0.0;

//////////////////////////////////////////////////
// for J Atoms
//////////////////////////////////////////////////
      LCPOAtom *lcpoNeighbors;
      int numLcpoNeighbors;
      lcpoNeighborList.nextlist( &lcpoNeighbors, &numLcpoNeighbors );

      for (int j = 0; j < numLcpoNeighbors; j++) {
        Real xj = lcpoNeighbors[j].x;
        Real yj = lcpoNeighbors[j].y;
        Real zj = lcpoNeighbors[j].z;
        Real rj = lcpoNeighbors[j].r;

        // i-j coarse check if too far away
        dxij = xj-xi;
        dyij = yj-yi;
        dzij = zj-zi;
        r2ij = dxij*dxij + dyij*dyij + dzij*dzij;
        FLOPS(7);
        if (r2ij >= cut2 || r2ij < 0.01) { continue; }

        // i-j precise check if too far away
        FLOPS(5)
        BigReal rirj2 = ri+rj;
        rirj2 *= rirj2;
        if ( r2ij >= rirj2 ) { continue; }

        BigReal rij = sqrt(r2ij);
        BigReal rij_1 = 1.f / rij;
          
//////////////////////////////////////////////////
// S2
//////////////////////////////////////////////////
        BigReal Aij = calcOverlap(rij, ri, rj);
        AijSum += Aij;
        FLOPS(12)

        //for dAi_drj force calculation
        BigReal dAijdrij = PI*ri*(rij_1*rij_1*(ri*ri-rj*rj)-1);
        BigReal dAijdrijdxj = dAijdrij*dxij*rij_1; // g k' i' l'
        BigReal dAijdrijdyj = dAijdrij*dyij*rij_1;
        BigReal dAijdrijdzj = dAijdrij*dzij*rij_1;
        FLOPS(14)

        BigReal AjkjSum = 0; // i' l'
        BigReal dAjkdrjkdxjSum = 0.0;
        BigReal dAjkdrjkdyjSum = 0.0;
        BigReal dAjkdrjkdzjSum = 0.0;

//////////////////////////////////////////////////
// for K Atoms
//////////////////////////////////////////////////
        for (int k = 0; k < numLcpoNeighbors; k++) {
          Real xk = lcpoNeighbors[k].x;
          Real yk = lcpoNeighbors[k].y;
          Real zk = lcpoNeighbors[k].z;
          Real rk = lcpoNeighbors[k].r;

          // i-k coarse check if too far away
          dxik = xk-xi;
          dyik = yk-yi;
          dzik = zk-zi;
          r2ik = dxik*dxik + dyik*dyik + dzik*dzik;
          FLOPS(8)
          if (r2ik >= cut2 || r2ik < 0.01) { continue; }

          // j-k coarse check if too far away
          dxjk = xk-xj;
          dyjk = yk-yj;
          dzjk = zk-zj;
          r2jk = dxjk*dxjk + dyjk*dyjk + dzjk*dzjk;
          FLOPS(8)
          if (r2jk >= cut2 || r2jk < 0.01) { continue; }

          // i-k precise check if too far away
          FLOPS(3)
          BigReal rirk2 = ri+rk;
          rirk2 *= rirk2;
          if ( r2ik >= rirk2 ) { continue; }

          // j-k precise check if too far away
          FLOPS(2)
          BigReal rjrk2 = rj+rk;
          rjrk2 *= rjrk2;
          if ( r2jk >= rjrk2 ) { continue; }
          BigReal rjk  = sqrt(r2jk);

//////////////////////////////////////////////////
// S3
//////////////////////////////////////////////////
          BigReal rjk_1 = 1.0/rjk;
          BigReal Ajk = calcOverlap(rjk, rj, rk);
          FLOPS(12)
          AjkSum  += Ajk;
          AjkjSum += Ajk; // i' l'
          FLOPS(5)

//////////////////////////////////////////////////
// Force dAi_drk
//////////////////////////////////////////////////
          BigReal dAjkdrjk = PI*rj*rjk_1*(rjk_1*rjk_1*(rj*rj-rk*rk) - 1.f);//ef'
          BigReal dAjkdrjkdxj = -dAjkdrjk*dxjk; // e f h'
          BigReal dAjkdrjkdyj = -dAjkdrjk*dyjk;
          BigReal dAjkdrjkdzj = -dAjkdrjk*dzjk;
          lcpoNeighbors[k].f->x -= -dAjkdrjkdxj*(P3+P4*Aij)*surfTen; // e f
          lcpoNeighbors[k].f->y -= -dAjkdrjkdyj*(P3+P4*Aij)*surfTen;
          lcpoNeighbors[k].f->z -= -dAjkdrjkdzj*(P3+P4*Aij)*surfTen;

          dAjkdrjkdxjSum += dAjkdrjkdxj; // h j'
          dAjkdrjkdyjSum += dAjkdrjkdyj;
          dAjkdrjkdzjSum += dAjkdrjkdzj;
          FLOPS(34)

        } // k atoms
//////////////////////////////////////////////////
// S4
//////////////////////////////////////////////////
        AijAjkSum += Aij*AjkjSum;

//////////////////////////////////////////////////
// Force dAi_drj
//////////////////////////////////////////////////
        BigReal lastxj = dAijdrijdxj*AjkjSum + Aij*dAjkdrjkdxjSum; // i j
        BigReal lastyj = dAijdrijdyj*AjkjSum + Aij*dAjkdrjkdyjSum;
        BigReal lastzj = dAijdrijdzj*AjkjSum + Aij*dAjkdrjkdzjSum;
        BigReal dAidxj = (P2*dAijdrijdxj + P3*dAjkdrjkdxjSum + P4*lastxj);//ghij
        BigReal dAidyj = (P2*dAijdrijdyj + P3*dAjkdrjkdyjSum + P4*lastyj);
        BigReal dAidzj = (P2*dAijdrijdzj + P3*dAjkdrjkdzjSum + P4*lastzj);
        lcpoNeighbors[j].f->x -= dAidxj*surfTen;
        lcpoNeighbors[j].f->y -= dAidyj*surfTen;
        lcpoNeighbors[j].f->z -= dAidzj*surfTen;

        //for dAi_dri force calculation
        dAijdrijdxiSum -= dAijdrijdxj; // k
        dAijdrijdyiSum -= dAijdrijdyj;
        dAijdrijdziSum -= dAijdrijdzj;
        dAijdrijdxiAjkSum -= dAijdrijdxj*AjkjSum; // l
        dAijdrijdyiAjkSum -= dAijdrijdyj*AjkjSum;
        dAijdrijdziAjkSum -= dAijdrijdzj*AjkjSum;
        FLOPS(41)
      } // j atoms

//////////////////////////////////////////////////
// Force dAi_dri
//////////////////////////////////////////////////
      BigReal dAidxi = (P2*dAijdrijdxiSum + P4*dAijdrijdxiAjkSum); // k l
      BigReal dAidyi = (P2*dAijdrijdyiSum + P4*dAijdrijdyiAjkSum);
      BigReal dAidzi = (P2*dAijdrijdziSum + P4*dAijdrijdziAjkSum);
      force[pI]->f[Results::nbond][iIndex].x -= dAidxi*surfTen;
      force[pI]->f[Results::nbond][iIndex].y -= dAidyi*surfTen;
      force[pI]->f[Results::nbond][iIndex].z -= dAidzi*surfTen;

//////////////////////////////////////////////////
// Atom I Surface Area
//////////////////////////////////////////////////
      BigReal SAi = P1*S1 + P2*AijSum + P3*AjkSum + P4*AijAjkSum;
      //CkPrintf("SurfArea[%05d] = % 7.3f\n",idi,SAi);
      //SAi = (SAi > 0) ? SAi : 0;
      totalSurfaceArea += SAi;
      FLOPS(22)
    } // for inAtoms
  } // for patches I
#ifdef COUNT_FLOPS
  double t_stop = 1.0*clock()/CLOCKS_PER_SEC;
  CkPrintf("LCPO_TIME_F %7.3f Gflops %9d @ %f\n", 1e-9*flops/(t_stop-t_start),flops, (t_stop-t_start));
#endif

//////////////////////////////////////////////////
//  end calculation by submitting reduction
//////////////////////////////////////////////////

  reduction->item(REDUCTION_COMPUTE_CHECKSUM) += 1.;
  reduction->item(REDUCTION_ELECT_ENERGY) += totalSurfaceArea * surfTen;
  reduction->submit();

}//end do Force


// Lookup table for lcpo paramters
//indices 0 -> 22 are determined in Molecule.C
const Real ComputeLCPO::lcpoParams[23][5] = { //                      neigh
    { 0.00, 0.0000e+00,  0.0000e+00,  0.0000e+00, 0.0000e+00 }, //  0 H
    { 1.70, 7.7887e-01, -2.8063e-01, -1.2968e-03, 3.9328e-04 }, //  1 C sp3 1
    { 1.70, 5.6482e-01, -1.9608e-01, -1.0219e-03, 2.6580e-04 }, //  2 C sp3 2
    { 1.70, 2.3348e-01, -7.2627e-02, -2.0079e-04, 7.9670e-05 }, //  3 C sp3 3
    { 1.70, 0.0000e+00,  0.0000e+00,  0.0000e+00, 0.0000e+00 }, //  4 C sp3 4
    { 1.70, 5.1245e-01, -1.5966e-01, -1.9781e-04, 1.6392e-04 }, //  5 C sp2 2
    { 1.70, 7.0344e-02, -1.9015e-02, -2.2009e-05, 1.6875e-05 }, //  6 C sp2 3
    { 1.60, 7.7914e-01, -2.5262e-01, -1.6056e-03, 3.5071e-04 }, //  7 O sp3 1
    { 1.60, 4.9392e-01, -1.6038e-01, -1.5512e-04, 1.6453e-04 }, //  8 O sp3 2
    { 1.60, 6.8563e-01, -1.8680e-01, -1.3557e-03, 2.3743e-04 }, //  9 O sp2 1
    { 1.60, 8.8857e-01, -3.3421e-01, -1.8683e-03, 4.9372e-04 }, // 10 O=C-O
    { 1.65, 7.8602e-02, -2.9198e-01, -6.5370e-04, 3.6247e-04 }, // 11 N sp3 1
    { 1.65, 2.2599e-01, -3.6648e-02, -1.2297e-03, 8.0038e-05 }, // 12 N sp3 2
    { 1.65, 5.1481e-02, -1.2603e-02, -3.2006e-04, 2.4774e-05 }, // 13 N sp3 3
    { 1.65, 7.3511e-01, -2.2116e-01, -8.9148e-04, 2.5230e-04 }, // 14 N sp2 1
    { 1.65, 4.1102e-01, -1.2254e-01, -7.5448e-05, 1.1804e-04 }, // 15 N sp2 2
    { 1.65, 6.2577e-02, -1.7874e-02, -8.3120e-05, 1.9849e-05 }, // 16 N sp2 3
    { 1.90, 7.7220e-01, -2.6393e-01,  1.0629e-03, 2.1790e-04 }, // 17 S     1
    { 1.90, 5.4581e-01, -1.9477e-01, -1.2873e-03, 2.9247e-04 }, // 18 S     2
    { 1.90, 3.8650e-01, -1.8249e-01, -3.6598e-03, 4.2640e-04 }, // 19 P     3
    { 1.90, 3.8730e-02, -8.9339e-03,  8.3582e-06, 3.0381e-06 }, // 20 P     4
    { 1.80, 9.8318e-01, -4.0437e-01,  1.1249e-04, 4.9901e-04 }, // 21 Cl
    { 1.18, 4.9392e-01, -1.6038e-01, -1.5512e-04, 1.6453e-04 }  // 22 Mg
};
