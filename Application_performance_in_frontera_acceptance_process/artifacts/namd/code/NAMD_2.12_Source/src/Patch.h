/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef PATCH_H
#define PATCH_H

#include "NamdTypes.h"
#include "OwnerBox.h"
#include "Box.h"
#include "UniqueSortedArray.h"
#include "Lattice.h"
#include "PatchTypes.h"

#ifdef NAMD_MIC
// defined here to avoid including ComputeNonbondedMICKernel.h
#define MIC_MAX_DEVICES_PER_NODE          ( 16 )
#endif

typedef SortedArray<Compute*> ComputePtrList;

class Compute;
class Sequencer;
class PatchMap;
class AtomMapper;
class SubmitReduction;

// This the base class of homepatches and proxy patches. It maintains
// common functions of these patches. These include managing dependences
// between compute (force) objects and the patch and updating atom map.

class Patch
{
  public:

     Patch(PatchID pd);
     int hasNewAtoms() { return _hasNewAtoms; }
     virtual ~Patch();

     // methods for use by Compute objects
     Box<Patch,CompAtom>* registerPositionPickup(Compute *cid);
     void unregisterPositionPickup(Compute *cid,
				   Box<Patch,CompAtom>**const box);
     Box<Patch,CompAtom>* registerAvgPositionPickup(Compute *cid);
     void unregisterAvgPositionPickup(Compute *cid,
				   Box<Patch,CompAtom>**const box);
     // BEGIN LA
     Box<Patch,CompAtom>* registerVelocityPickup(Compute *cid);
     void unregisterVelocityPickup(Compute *cid,
                                  Box<Patch,CompAtom>**const box);
     // END LA

    //begin gbis
    Box<Patch,Real>* registerIntRadPickup(Compute *cid);
    void unregisterIntRadPickup(Compute *cid, Box<Patch,Real>**const box);

    Box<Patch,GBReal>* registerPsiSumDeposit(Compute *cid);
    void unregisterPsiSumDeposit(Compute *cid, Box<Patch,GBReal>**const box);

    Box<Patch,Real>* registerBornRadPickup(Compute *cid);
    void unregisterBornRadPickup(Compute *cid, Box<Patch,Real>**const box);

    Box<Patch,GBReal>* registerDEdaSumDeposit(Compute *cid);
    void unregisterDEdaSumDeposit(Compute *cid,Box<Patch,GBReal> **const box);

    Box<Patch,Real>* registerDHdrPrefixPickup(Compute *cid);
    void unregisterDHdrPrefixPickup(Compute *cid, Box<Patch,Real>**const box);
     //end gbis

    //LCPO
    Box<Patch,int>* registerLcpoTypePickup(Compute *cid);
    void unregisterLcpoTypePickup(Compute *cid, Box<Patch,int>**const box);

     Box<Patch,Results>* registerForceDeposit(Compute *cid);
     void unregisterForceDeposit(Compute *cid, Box<Patch,Results> **const box);

     // methods for use by Sequencer or ProxyManager
     // void positionsReady(void) { positionsReady(0); }
     void positionsReady(int n=0);

     // methods for Box callbacks
     void positionBoxClosed(void);
     void forceBoxClosed(void);
     void avgPositionBoxClosed(void);
     // BEGIN LA
     void velocityBoxClosed(void);
     // END LA

     //begin gbis
     void intRadBoxClosed(void);// intrinsic radii
     void psiSumBoxClosed(void);// sum screening 
     void bornRadBoxClosed(void);// born radius
     void dEdaSumBoxClosed(void);// sum dEda contributions
     void dHdrPrefixBoxClosed(void);//dHdr prefix
     void gbisP2Ready();
     void gbisP3Ready();
     //end gbis

     //LCPO
     void lcpoTypeBoxClosed(void);

     int getNumAtoms() { return numAtoms; }

     // DMK - Atom Separation (water vs. non-water)
     #if NAMD_SeparateWaters != 0
       int getNumWaterAtoms() { return numWaterAtoms; }
     #endif

     int getNumFixedAtoms() { return numFixedAtoms; }  // not updated
     void setNumFixedAtoms(int numFixed) { numFixedAtoms=numFixed; }  // not updated
     PatchID getPatchID() { return patchID; }
     int getNumComputes() { return positionComputeList.size(); }

     CompAtomExt* getCompAtomExtInfo() { return pExt.begin(); }
#ifdef NAMD_KNL
     CompAtomFlt* getCompAtomFlt() { return pFlt.begin(); }
#endif
     CudaAtom* getCudaAtomList() { return cudaAtomPtr; }

     Lattice &lattice;
     Flags flags;

     // DMK - NOTE : Just placing the variables in public for now so only one location, move to protected if this actually helps performance
     #if defined(NAMD_MIC) // NOTE: Used for submit atoms on arrival
       pthread_mutex_t mic_atomData_mutex;
       void* mic_atomData;
       void* mic_atomData_prev[MIC_MAX_DEVICES_PER_NODE];
       int mic_atomData_seq;
       int mic_atomData_deviceSeq[MIC_MAX_DEVICES_PER_NODE];
       uint64_t mic_atomData_devicePtr[MIC_MAX_DEVICES_PER_NODE];
       int mic_atomData_allocSize_host;
       int mic_atomData_allocSize_device[MIC_MAX_DEVICES_PER_NODE];
     #endif

  protected:

     const PatchID patchID;
     int           numAtoms;
     int           numFixedAtoms;
     CompAtomList  p;
     CompAtomList  p_avg;
     // BEGIN LA
     CompAtomList  v;
     // END LA

     AtomMapper *atomMapper;

     // begin gbis
     RealList intRad;
     GBRealList psiSum;
     GBRealList psiFin;
     RealList bornRad;
     RealList dHdrPrefix;
     GBRealList dEdaSum;
     // end gbis

    //LCPO
    IntList lcpoType;

     // DMK - Atom Separation (water vs. non-water)
     #if NAMD_SeparateWaters != 0
       int numWaterAtoms;  // Set numWaters to the number of water atoms at
                           //   the lead of the atoms list.  If numWaters is
                           //   set to -1, this should indicate that
                           //   atoms has not been separated yet.
     #endif

     CompAtomExtList pExt;
#ifdef NAMD_KNL
     CompAtomFltList pFlt;
#endif

#ifdef REMOVE_PROXYDATAMSG_EXTRACOPY
     //1. Those fields are declared for reusing position info
     //inside the ProxyDataMsg msg at every step so that the
     //extra copy is avoided.
     //Regarding the CompAtomExt list inside the msg of ProxyAllMsg type
     //we cannot avoid the copy in the current scheme because this information
     //will be lost as the msg will be deleted at the next timestep. But the
     //overhead is amortized among the steps that atoms don't migrate
     //2. positionPtrBegin is better to be made 32-byte aligned so we could
     // have better cache performance in the force calculation part. This
     // is especially needed for BG/L machine.
     // --Chao Mei
     CompAtom      *positionPtrBegin;
     CompAtom      *positionPtrEnd;     
#endif
     CompAtom      *avgPositionPtrBegin;
     CompAtom      *avgPositionPtrEnd;

     // BEGIN LA
     CompAtom      *velocityPtrBegin;
     CompAtom      *velocityPtrEnd;
     // END LA

     CudaAtom      *cudaAtomPtr;

     ForceList     f[Results::maxNumForces];
     Results	   results;

     int computesSortedByPriority;
     int firstHoldableCompute;

     OwnerBox<Patch,CompAtom> positionBox;
     ComputePtrList              positionComputeList;
     OwnerBox<Patch,CompAtom> avgPositionBox;
     ComputePtrList              avgPositionComputeList;
     // BEGIN LA
     OwnerBox<Patch,CompAtom> velocityBox;
     ComputePtrList              velocityComputeList;
     // END LA

     //begin gbis
     OwnerBox<Patch,Real>    intRadBox;
     ComputePtrList           intRadComputeList;
     OwnerBox<Patch,GBReal>  psiSumBox;
     ComputePtrList           psiSumComputeList;
     OwnerBox<Patch,Real>    bornRadBox;
     ComputePtrList           bornRadComputeList;
     OwnerBox<Patch,GBReal>  dEdaSumBox;
     ComputePtrList           dEdaSumComputeList;
     OwnerBox<Patch,Real>    dHdrPrefixBox;
     ComputePtrList           dHdrPrefixComputeList;
     //end gbis

    //LCPO
     OwnerBox<Patch,int>    lcpoTypeBox;
     ComputePtrList          lcpoTypeComputeList;

     OwnerBox<Patch,Results>    forceBox;
     ComputePtrList              forceComputeList;

     virtual void boxClosed(int /* box */) = 0;
     int boxesOpen;

     int _hasNewAtoms;

#ifdef NODEAWARE_PROXY_SPANNINGTREE    
    //its own children in proxy tree
    #ifdef USE_NODEPATCHMGR
    //the immediate children (in terms of node id) also cotains two parts
    //as the above variable shows
    //If this patch has proxies residing on the same node, then the last entry
    //of "nodeChildren" stores this node id. It is same with that variable
    //in ProxyPatch      
    //If this proxy resides on node A, then the last entry
    //of "nodeChildren" has to be A
    //It is same with that variable in HomePatch
    int *nodeChildren;
    int numNodeChild;
    #endif
#endif
	int *child;
	int nChild;

  private:

    SubmitReduction *reduction;

};


#endif

