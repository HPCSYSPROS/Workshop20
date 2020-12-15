/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "InfoStream.h"
#include "Patch.h"
#include "PatchMap.h"
#include "Compute.h"

#include "AtomMap.h"
#include "ComputeMap.h"
#include "Node.h"
#include "Molecule.h"
#include "SimParameters.h"
#include "ResizeArrayPrimIter.h"
#include "ComputeNonbondedMICKernel.h"
#include "Priorities.h"
#include "ReductionMgr.h"

#include "Sync.h"

#include <algorithm>

typedef ResizeArrayPrimIter<Compute*> ComputePtrListIter;

struct cptr_sortop_priority {
  bool operator() (Compute *i, Compute *j) const {
    return ( i->priority() < j->priority() );
  }
};

//#define  DEBUGM
#define MIN_DEBUG_LEVEL 4
#include "Debug.h"

Patch::~Patch() {
  delete atomMapper;
#ifndef NAMD_CUDA
  delete reduction;
#endif
}

Patch::Patch(PatchID pd) :
   lattice(flags.lattice),
   patchID(pd), numAtoms(0), numFixedAtoms(0),
   avgPositionPtrBegin(0), avgPositionPtrEnd(0),
   velocityPtrBegin(0), velocityPtrEnd(0),	// BEGIN LA, END LA
   positionBox(this,&Patch::positionBoxClosed,pd,0),
   avgPositionBox(this,&Patch::avgPositionBoxClosed,pd,3),
   velocityBox(this,&Patch::velocityBoxClosed,pd,4), // BEGIN LA, END LA
   psiSumBox(this,&Patch::psiSumBoxClosed,pd,5), // begin gbis
   intRadBox(this,&Patch::intRadBoxClosed,pd,6),
   bornRadBox(this,&Patch::bornRadBoxClosed,pd,7),
   dEdaSumBox(this,&Patch::dEdaSumBoxClosed,pd,8),
   dHdrPrefixBox(this,&Patch::dHdrPrefixBoxClosed,pd,9), //end gbis
   lcpoTypeBox(this,&Patch::lcpoTypeBoxClosed,pd,10),
   forceBox(this,&Patch::forceBoxClosed,pd,1),
   boxesOpen(0), _hasNewAtoms(0),
   computesSortedByPriority(0), firstHoldableCompute(0)

   // DMK - Atom Separation (water vs. non-water)
   #if NAMD_SeparateWaters != 0
     ,numWaterAtoms(-1)
   #endif
{
  //CkPrintf("GBIS: PatchCreated\n");
#ifdef REMOVE_PROXYDATAMSG_EXTRACOPY
    positionPtrBegin = 0;
    positionPtrEnd = 0;
#endif

	nChild = 0;
	child = NULL;
#ifdef NODEAWARE_PROXY_SPANNINGTREE
  #ifdef USE_NODEPATCHMGR
  nodeChildren = NULL;
  numNodeChild = 0;
  #endif
#endif

  lattice = Node::Object()->simParameters->lattice;
  atomMapper = new AtomMapper(pd);
#ifndef NAMD_CUDA
  reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);
#endif

  // DMK
  #if defined(NAMD_MIC)
    cudaAtomPtr = NULL;
    #if MIC_SUBMIT_ATOMS_ON_ARRIVAL != 0
      pthread_mutex_init(&mic_atomData_mutex, NULL);
      mic_atomData = NULL;
      mic_atomData_seq = -1;
      mic_atomData_allocSize_host = 0;
      for (int i = 0; i < MIC_MAX_DEVICES_PER_NODE; i++) {
        mic_atomData_prev[i] = NULL;
        mic_atomData_deviceSeq[i] = -1;
        mic_atomData_devicePtr[i] = 0;
        mic_atomData_allocSize_device[i] = 0;
      }
    #endif
  #endif
}

Box<Patch,CompAtom>* Patch::registerPositionPickup(Compute *cid)
{
   if ( computesSortedByPriority ) {
     positionComputeList.sort();
     computesSortedByPriority = 0;
   }
   //DebugM(4, "registerPositionPickupa("<<patchID<<") from " << cid->cid << "\n");
   if (positionComputeList.add(cid) < 0)
   {
     DebugM(7, "registerPositionPickup() failed for cid " << cid->cid << std::endl);
     return NULL;
   }
   return positionBox.checkOut(cid->cid);
}

void Patch::unregisterPositionPickup(Compute *cid, Box<Patch,CompAtom> **const box)
{
   if ( computesSortedByPriority ) {
     positionComputeList.sort();
     computesSortedByPriority = 0;
   }
   DebugM(4, "UnregisterPositionPickup from " << cid->cid << "\n");
   positionComputeList.del(cid);
   positionBox.checkIn(*box);
   *box = 0;
}

Box<Patch,CompAtom>* Patch::registerAvgPositionPickup(Compute *cid)
{
   //DebugM(4, "registerAvgPositionPickup("<<patchID<<") from " << cid->cid << "\n");
   return avgPositionBox.checkOut(cid->cid);
}

void Patch::unregisterAvgPositionPickup(Compute *cid, Box<Patch,CompAtom> **const box)
{
   DebugM(4, "UnregisterAvgPositionPickup from " << cid->cid << "\n");
   avgPositionBox.checkIn(*box);
   *box = 0;
}

// BEGIN LA
Box<Patch,CompAtom>* Patch::registerVelocityPickup(Compute *cid)
{
   //DebugM(4, "registerVelocityPickup("<<patchID<<") from " << cid->cid << "\n");
   return velocityBox.checkOut(cid->cid);
}

void Patch::unregisterVelocityPickup(Compute *cid, Box<Patch,CompAtom> **const box)
{
   DebugM(4, "UnregisterVelocityPickup from " << cid->cid << "\n");
   velocityBox.checkIn(*box);
   *box = 0;
}
// END LA

//begin gbis
//deposit, not pickup
Box<Patch,GBReal>* Patch::registerPsiSumDeposit(Compute *cid) {

  if (psiSumComputeList.add(cid) < 0) {
    DebugM(7, "registerPsiSumDeposit() failed for cid " << cid->cid << std::endl);
    DebugM(7, "  size of psiSumCompueList " << psiSumComputeList.size() << std::endl);
     return NULL;
  }
  return psiSumBox.checkOut(cid->cid);
}

void Patch::unregisterPsiSumDeposit(Compute *cid,Box<Patch,GBReal> **const box) {
  psiSumComputeList.del(cid);
  psiSumBox.checkIn(*box);
  *box = 0;
}
Box<Patch,Real>* Patch::registerIntRadPickup(Compute *cid) {
  return intRadBox.checkOut(cid->cid);
}
void Patch::unregisterIntRadPickup(Compute *cid,Box<Patch,Real> **const box) {
  intRadBox.checkIn(*box);
  *box = 0;
}

//LCPO
Box<Patch,int>* Patch::registerLcpoTypePickup(Compute *cid) {
  return lcpoTypeBox.checkOut(cid->cid);
}
void Patch::unregisterLcpoTypePickup(Compute *cid,Box<Patch,int> **const box) {
  lcpoTypeBox.checkIn(*box);
  *box = 0;
}

Box<Patch,Real>* Patch::registerBornRadPickup(Compute *cid) {
  return bornRadBox.checkOut(cid->cid);
}
void Patch::unregisterBornRadPickup(Compute *cid,Box<Patch,Real> **const box) {
  bornRadBox.checkIn(*box);
  *box = 0;
}

Box<Patch,GBReal>* Patch::registerDEdaSumDeposit(Compute *cid) {
  if (dEdaSumComputeList.add(cid) < 0) {
    DebugM(7, "registerDEdaSumDeposit() failed for cid " << cid->cid << std::endl);
    DebugM(7, "  size of dEdaSumCompueList " << dEdaSumComputeList.size() << std::endl);
     return NULL;
  }
  return dEdaSumBox.checkOut(cid->cid);
}
void Patch::unregisterDEdaSumDeposit(Compute *cid,Box<Patch,GBReal> **const box){
  dEdaSumComputeList.del(cid);
  dEdaSumBox.checkIn(*box);
  *box = 0;
}

Box<Patch,Real>* Patch::registerDHdrPrefixPickup(Compute *cid)
{
  return dHdrPrefixBox.checkOut(cid->cid);
}
void Patch::unregisterDHdrPrefixPickup(Compute *cid,Box<Patch,Real> **const box) {
  dHdrPrefixBox.checkIn(*box);
  *box = 0;
}
//end gbis

Box<Patch,Results>* Patch::registerForceDeposit(Compute *cid)
{
   if (forceComputeList.add(cid) < 0)
   {
     DebugM(7, "registerForceDeposit() failed for cid " << cid->cid << std::endl);
     DebugM(7, "  size of forceCompueList " << forceComputeList.size() << std::endl);
     return NULL;
   }
   return forceBox.checkOut(cid->cid);
}

void Patch::unregisterForceDeposit(Compute *cid, Box<Patch,Results> **const box)
{
   DebugM(4, "unregisterForceDeposit() computeID("<<cid<<")"<<std::endl);
   forceComputeList.del(cid);
   forceBox.checkIn(*box);
   *box = 0;
}

void Patch::positionBoxClosed(void)
{
   //positionPtrBegin = 0;
   this->boxClosed(0);
}

void Patch::forceBoxClosed(void)
{
   DebugM(4, "patchID("<<patchID<<") forceBoxClosed! call\n");
#ifndef NAMD_CUDA
   // calculate direct nonbonded virial from segregated forces and aggregate forces
   const Vector center = lattice.unscale( PatchMap::Object()->center(patchID) );
#ifdef REMOVE_PROXYDATAMSG_EXTRACOPY
   const CompAtom * const __restrict pd = positionPtrBegin;
#else
   const CompAtom * const __restrict pd = p.begin();
#endif
   const int n = numAtoms;
   int roff = 0;
   for ( int j1 = Results::nbond; j1 <= Results::slow; ++j1, roff += REDUCTION_VIRIAL_SLOW_XX - REDUCTION_VIRIAL_NBOND_XX ) {
      int j2 = j1 + ( Results::nbond_virial - Results::nbond );
      Force * __restrict f1 = results.f[j1];
      Force * __restrict f2 = results.f[j2];
      BigReal virial_xx = 0.;
      BigReal virial_xy = 0.;
      BigReal virial_xz = 0.;
      BigReal virial_yy = 0.;
      BigReal virial_yz = 0.;
      BigReal virial_zz = 0.;
#pragma simd reduction(+:virial_xx,virial_xy,virial_xz,virial_yy,virial_yz,virial_zz)
#pragma ivdep
      for ( int i=0; i<n; ++i ) {
        BigReal p_x = pd[i].position.x - center.x;
        BigReal p_y = pd[i].position.y - center.y;
        BigReal p_z = pd[i].position.z - center.z;
        BigReal f_x = f2[i].x;
        BigReal f_y = f2[i].y;
        BigReal f_z = f2[i].z;
        virial_xx += f_x * p_x;
        virial_xy += f_x * p_y;
        virial_xz += f_x * p_z;
        virial_yy += f_y * p_y;
        virial_yz += f_y * p_z;
        virial_zz += f_z * p_z;
        f1[i].x += f_x;
        f1[i].y += f_y;
        f1[i].z += f_z;
      }
      f[j2].resize(0);
      reduction->item(roff + REDUCTION_VIRIAL_NBOND_XX) += virial_xx;
      reduction->item(roff + REDUCTION_VIRIAL_NBOND_XY) += virial_xy;
      reduction->item(roff + REDUCTION_VIRIAL_NBOND_XZ) += virial_xz;
      reduction->item(roff + REDUCTION_VIRIAL_NBOND_YX) += virial_xy;
      reduction->item(roff + REDUCTION_VIRIAL_NBOND_YY) += virial_yy;
      reduction->item(roff + REDUCTION_VIRIAL_NBOND_YZ) += virial_yz;
      reduction->item(roff + REDUCTION_VIRIAL_NBOND_ZX) += virial_xz;
      reduction->item(roff + REDUCTION_VIRIAL_NBOND_ZY) += virial_yz;
      reduction->item(roff + REDUCTION_VIRIAL_NBOND_ZZ) += virial_zz;
   }
   reduction->submit();
#endif

   for (int j = 0; j < Results::maxNumForces; ++j )
   {
     results.f[j] = 0;
   }
   this->boxClosed(1);
}

void Patch::avgPositionBoxClosed(void)
{
   avgPositionPtrBegin = 0;
   this->boxClosed(3);
}

// BEGIN LA
void Patch::velocityBoxClosed(void)
{
   DebugM(4, "patchID("<<patchID<<") velocityBoxClosed! call\n");
   velocityPtrBegin = 0;
   this->boxClosed(4);	// ?? Don't know about number
}
// END LA

// void Patch::boxClosed(int box) is virtual

// begin gbis
void Patch::psiSumBoxClosed(void) {
  this->boxClosed(5);
}
void Patch::intRadBoxClosed(void) {
   //dHdrPrefixPtr = 0;
   this->boxClosed(6);
}
void Patch::bornRadBoxClosed(void) {
   //bornRadPtr = 0;
   this->boxClosed(7);
}
void Patch::dEdaSumBoxClosed(void) {
   //dEdaSumPtr = 0;
   this->boxClosed(8);
}
void Patch::dHdrPrefixBoxClosed(void) {
   //dHdrPrefixPtr = 0;
   this->boxClosed(9);
}
// end gbis

//LCPO
void Patch::lcpoTypeBoxClosed(void) {
   this->boxClosed(10);
}

void Patch::positionsReady(int doneMigration)
{
   DebugM(4,"Patch::positionsReady() - patchID(" << patchID <<")"<<std::endl );

   if ( doneMigration ){
// #ifdef REMOVE_PROXYDATAMSG_EXTRACOPY
//       AtomMap::Object()->registerIDs(patchID,positionPtrBegin,positionPtrEnd);       
// #else
       atomMapper->registerIDsCompAtomExt(pExt.begin(),pExt.end());
// #endif

   }

#ifdef NAMD_KNL
   {
     const Vector center = lattice.unscale( PatchMap::Object()->center(patchID) );
     const int n = numAtoms;
     pFlt.resize(n);
     CompAtomFlt * const pf = pFlt.begin();
#ifdef REMOVE_PROXYDATAMSG_EXTRACOPY
     const CompAtom * const pd = positionPtrBegin;
#else
     const CompAtom * const pd = p.begin();
#endif
     for ( int i=0; i<n; ++i ) {
       // need to subtract center in double precision, then assign to float
       pf[i].position.x = pd[i].position.x - center.x;
       pf[i].position.y = pd[i].position.y - center.y;
       pf[i].position.z = pd[i].position.z - center.z;
       pf[i].vdwType = pd[i].vdwType;
     }
   }
#endif

   boxesOpen = 2;
   if ( flags.doMolly ) boxesOpen++;
   // BEGIN LA
   if (flags.doLoweAndersen) {
       DebugM(4, "Patch::positionsReady, flags.doMolly = " << flags.doMolly << "\n");
       boxesOpen++;
   }
   // END LA
   _hasNewAtoms = (doneMigration != 0);

#if CMK_BLUEGENEL
   CmiNetworkProgressAfter (0);
#endif

   // Give all position pickup boxes access to positions
   //positionPtrBegin = p.begin();
#ifdef REMOVE_PROXYDATAMSG_EXTRACOPY
   positionBox.open(positionPtrBegin);
#else
   positionBox.open(p.begin());
#endif
   if ( flags.doMolly ) {
     //avgPositionPtrBegin = p_avg.begin();
     avgPositionBox.open(avgPositionPtrBegin);
   }
   
   // BEGIN LA
   if (flags.doLoweAndersen) {
       velocityBox.open(velocityPtrBegin);
   }
   // END LA
   // begin gbis
    if (flags.doGBIS) {
      boxesOpen += 5;
      //intRad should already be taken care of
      intRadBox.open(intRad.begin());
      psiSum.resize(numAtoms);//resize array
      psiSum.setall(0);
      psiSumBox.open(psiSum.begin());
      psiFin.resize(numAtoms);//has no box
      psiFin.setall(0);
      bornRad.resize(numAtoms);
      bornRad.setall(0);
      bornRadBox.open(bornRad.begin());
      dEdaSum.resize(numAtoms);//resize array
      dEdaSum.setall(0);
      dEdaSumBox.open(dEdaSum.begin());
      dHdrPrefix.resize(numAtoms);
      dHdrPrefix.setall(0);
      dHdrPrefixBox.open(dHdrPrefix.begin());
    }
   // end gbis

  //LCPO
  if (flags.doLCPO) {
    boxesOpen++;
    lcpoTypeBox.open(lcpoType.begin());
  }

#if CMK_BLUEGENEL
   CmiNetworkProgressAfter (0);
#endif
   
   // Give all force deposit boxes access to forces
   Force *forcePtr;
   for ( int j = 0; j < Results::maxNumForces; ++j )
    {
      f[j].resize(numAtoms);
      forcePtr = f[j].begin();
      memset (forcePtr, 0, sizeof (Force) * numAtoms);
      results.f[j] = forcePtr;
    }
   forceBox.open(&results);

   if ( ! computesSortedByPriority ) {
     if (positionComputeList.size() == 0 && PatchMap::Object()->node(patchID) != CkMyPe()) {
       iout << iINFO << "PATCH_COUNT: Patch " << patchID 
	    << " on PE " << CkMyPe() <<" home patch " 
	    << PatchMap::Object()->node(patchID)
	    << " does not have any computes\n" 
	    << endi;
     }

     cptr_sortop_priority so;
     std::sort(positionComputeList.begin(), positionComputeList.end(), so);
     computesSortedByPriority = 1;
     int i;
     for ( i=0; i<positionComputeList.size(); ++i ) {
       if ( positionComputeList[i]->priority() > PME_PRIORITY ) break;
     }
     firstHoldableCompute = i;
   }

   int seq = flags.sequence;

   // Iterate over compute objects that need to be informed we are ready
   Compute **cid = positionComputeList.begin();
   for ( int i=0; i < firstHoldableCompute; ++i, ++cid ) {
     (*cid)->patchReady(patchID,doneMigration,seq);
   }
   Compute **cend = positionComputeList.end();

   // gzheng
   if (Sync::Object()->holdComputes(patchID, cid, cend, doneMigration, seq)) {
     return;
   }

   for( ; cid != cend; cid++ ) {
     (*cid)->patchReady(patchID,doneMigration,seq);
   }
}

// begin gbis

void Patch::gbisP2Ready() {
 ComputePtrListIter cid(positionComputeList);

  int seq = flags.sequence;
  for(cid = cid.begin(); cid != cid.end(); cid++) {
    if ( (*cid)->type() == computeNonbondedSelfType ||
         (*cid)->type() == computeNonbondedPairType ||
         (*cid)->type() == computeNonbondedCUDAType
#ifdef NAMD_CUDA
         || (*cid)->type() == computeNonbondedCUDA2Type
#endif
         ) {
      (*cid)->gbisP2PatchReady(patchID,seq);
    }
  }
}

void Patch::gbisP3Ready() {

  ComputePtrListIter cid(positionComputeList);

  int seq = flags.sequence;
  for(cid = cid.begin(); cid != cid.end(); cid++) {
    if ( (*cid)->type() == computeNonbondedSelfType ||
         (*cid)->type() == computeNonbondedPairType ||
         (*cid)->type() == computeNonbondedCUDAType
#ifdef NAMD_CUDA
         || (*cid)->type() == computeNonbondedCUDA2Type
#endif
         ) {
      (*cid)->gbisP3PatchReady(patchID,seq);
    }
  }
}

//end gbis
