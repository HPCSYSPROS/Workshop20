/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/
 
/*****************************************************************************
 * $Source: /home/cvs/namd/cvsroot/namd2/src/LdbCoordinator.C,v $
 * $Author: jim $
 * $Date: 2016/03/02 21:33:06 $
 * $Revision: 1.127 $
 *****************************************************************************/

#include <stdlib.h>

#include "InfoStream.h"
#include "NamdCentLB.h"
#include "NamdHybridLB.h"
#include "NamdDummyLB.h"
#include "NamdNborLB.h"

#include "HomePatch.h"
#include "LdbCoordinator.decl.h"
#include "LdbCoordinator.h"
#include "NamdTypes.h"
#include "Node.h"
#include "SimParameters.h"
#include "PatchMap.inl"
#include "ComputeMap.h"
#include "ComputeNonbondedMICKernel.h"
//#define DEBUGM
#define MIN_DEBUG_LEVEL 3
#include "Debug.h"
#include "Controller.h"
#include "Sequencer.h"
#include "RefineOnly.h"
#include "ComputeMgr.h"
#include "Compute.h"
#include "packmsg.h"
#include "Sync.h"

#include "elements.h"
#include "ComputeMgr.decl.h"

#define DEBUG_LEVEL 4

#if CONVERSE_VERSION_ELAN
extern "C" void enableBlockingReceives();
extern "C" void disableBlockingReceives();
#endif

void LdbCoordinator_initproc() {
  // Set the load balancing period (in seconds).  Without this the
  // load balancing framework will hang until 1 second has passed
  // since the last load balancing, causing hiccups in very fast runs.
  // This is duplicated below for older versions, but putting it here
  // also fixes the first load balance.
  LBSetPeriod(1.0e-5);
}

void LdbCoordinator::staticMigrateFn(LDObjHandle handle, int dest)
{
   LdbCoordinator *ldbCoordinator = (LdbCoordinator *)LDOMUserData(handle.omhandle);
   ldbCoordinator->Migrate(handle,dest);
}

void LdbCoordinator::Migrate(LDObjHandle handle, int dest)
{
  LdbMigrateMsg* msg = new LdbMigrateMsg;
  msg->handle = handle;
  msg->from = CkMyPe();
  msg->to = dest;
  if ( msg->to != CkMyPe() ) {
    CProxy_LdbCoordinator ldbProxy(thisgroup);
    ldbProxy[CkMyPe()].RecvMigrate(msg);
  } else {
    ExpectMigrate(msg);
  }
}

void LdbCoordinator::staticStatsFn(LDOMHandle h, int state)
{
  CkPrintf("I'm supposed to set stats\n");
}

void LdbCoordinator::staticQueryEstLoadFn(LDOMHandle h)
{
  CkPrintf("I'm supposed to query load\n");
}

void LdbCoordinator::staticReceiveAtSync(void* data)
{

#if CONVERSE_VERSION_ELAN
    //disableBlockingReceives();
#endif

  ((LdbCoordinator*)data)->ReceiveAtSync();
}

void LdbCoordinator::ReceiveAtSync()
{
  theLbdb->RegisteringObjects(myHandle);
}

void LdbCoordinator::staticResumeFromSync(void* data)
{
  ((LdbCoordinator*)data)->ResumeFromSync();
}

void LdbCoordinator::ResumeFromSync()
{
  theLbdb->DoneRegisteringObjects(myHandle);
  CkCallback cb(CkIndex_LdbCoordinator::nodeDone(NULL), 0, thisgroup);
  contribute(0, NULL, CkReduction::random, cb);
}

LdbCoordinator::LdbCoordinator()
{
  if (CkpvAccess(LdbCoordinator_instance) == NULL) {
    CkpvAccess(LdbCoordinator_instance) = this;
  } else {
    NAMD_bug("LdbCoordinator instanced twice on same node!");
  }
  
#if 0
  // Create a load balancer
  if (CkMyPe() == 0) {
    //   CreateCentralLB();
    CreateNamdCentLB();
    //   CreateNamdNborLB();
  }
#endif

  collPes = 0;
  ldbCycleNum = 1;
  takingLdbData = 1;
  totalStepsDone = 0;
  nLocalComputes = nLocalPatches = 0;
  patchNAtoms = (int *) NULL;
  sequencerThreads = (Sequencer **) NULL;
  ldbStatsFP = NULL;
  computeArray = NULL;
  patchArray = NULL;
  processorArray = NULL;

  // Register self as an object manager for new charm++ balancer framework
  theLbdb = LBDatabase::Object(); 

  // Set the load balancing period (in seconds).  Without this the
  // load balancing framework will hang until 1 second has passed
  // since the last load balancing, causing hiccups in very fast runs.
  // Unfortunately, the clock is already set for the first load
  // balancing, but only +LBPeriod 1.0e-5 can fix that in older charm.
  // For newer versions this is handled in initproc above.

  theLbdb->SetLBPeriod(1.0e-5);

  myOMid.id.idx = 1;
  LDCallbacks cb = { (LDMigrateFn)staticMigrateFn,
		     (LDStatsFn)staticStatsFn,
		     (LDQueryEstLoadFn)staticQueryEstLoadFn
                   };
  myHandle = theLbdb->RegisterOM(myOMid,(void*)this,cb);

  // Add myself as a local barrier receiver, so I know when I might
  // be registering objects.
  theLbdb->AddLocalBarrierReceiver((LDBarrierFn)staticReceiveAtSync,
				   (void*)this);;

  // Also, add a local barrier client, to trigger load balancing
  ldBarrierHandle = theLbdb->
    AddLocalBarrierClient((LDResumeFn)staticResumeFromSync,
			  (void*)this);
  migrateMsgs = 0; // linked list
  numComputes = 0;
  reg_all_objs = 1;
}

LdbCoordinator::~LdbCoordinator(void)
{
  delete [] patchNAtoms;
  delete [] sequencerThreads;
  if (CkMyPe() == 0)
  {
    delete [] computeArray;
    delete [] patchArray;
    delete [] processorArray;
  }
  if (ldbStatsFP)
    fclose(ldbStatsFP);

}

void LdbCoordinator::createLoadBalancer()
{
  const SimParameters *simParams = Node::Object()->simParameters;

  // Create hierarchical or centralized load balancers
  // Currently centralized is the default
  if (simParams->ldBalancer == LDBAL_CENTRALIZED) {
    CkPrintf("LDB: Central LB being created...\n");
    CreateNamdCentLB();
  } else if (simParams->ldBalancer == LDBAL_HYBRID) {
    CkPrintf("LDB: Hybrid LB being created...\n");
    CreateNamdHybridLB();
  }
}

void LdbCoordinator::initialize(PatchMap *pMap, ComputeMap *cMap, int reinit)
{
  const SimParameters *simParams = Node::Object()->simParameters;

#if 0
  static int lbcreated = 0; // XXX static variables are unsafe for SMP
  // PE0 first time Create a load balancer
  if (CkMyPe() == 0 && !lbcreated) {
    if (simParams->ldbStrategy == LDBSTRAT_ALGNBOR) 
      CreateNamdNborLB();
    else {
      //   CreateCentralLB();
      CreateNamdCentLB();
    }
    lbcreated = 1;
  }
#endif

  //  DebugM(10,"stepsPerLdbCycle initialized\n");
  stepsPerLdbCycle = simParams->ldbPeriod;
  firstLdbStep = simParams->firstLdbStep;
  int lastLdbStep = simParams->lastLdbStep;
  int stepsPerCycle = simParams->stepsPerCycle;

  computeMap = cMap;
  patchMap = pMap;

  // Set the number of received messages correctly for node 0

  nStatsMessagesExpected = Node::Object()->numNodes();
  nStatsMessagesReceived = 0;

  if (patchNAtoms) 
    delete [] patchNAtoms;  // Depends on delete NULL to do nothing
  nPatches = patchMap->numPatches();
  patchNAtoms = new int[nPatches];

  typedef Sequencer *seqPtr;

  if ( ! reinit ) {
    delete [] sequencerThreads;  // Depends on delete NULL to do nothing
    sequencerThreads = new seqPtr[nPatches];
  }

  nLocalPatches=0;

  int i;
  for(i=0;i<nPatches;i++)
  {
    if (patchMap->node(i) == Node::Object()->myid())
    {
      nLocalPatches++;
      patchNAtoms[i]=0;
    } else {
      patchNAtoms[i]=-1;
    }
    if ( ! reinit ) sequencerThreads[i]=NULL;
  }
  if ( ! reinit ) controllerThread = NULL;
  if (nLocalPatches != patchMap->numHomePatches())
    NAMD_die("Disaggreement in patchMap data.\n");
 
  const int oldNumComputes = numComputes;
  nLocalComputes = 0;
  numComputes = computeMap->numComputes();

  for(i=0;i<numComputes;i++)  {
    if ( (computeMap->node(i) == Node::Object()->myid())
	 && ( 0
              #if (defined(NAMD_CUDA) || defined(NAMD_MIC))
                #if defined(NAMD_MIC)
                  || ((computeMap->type(i) == computeNonbondedSelfType) && (computeMap->directToDevice(i) == 0))
                  || ((computeMap->type(i) == computeNonbondedPairType) && (computeMap->directToDevice(i) == 0))
                #endif
              #else
	      || (computeMap->type(i) == computeNonbondedSelfType)
	      || (computeMap->type(i) == computeNonbondedPairType)
#endif
	      || (computeMap->type(i) == computeLCPOType)
	      || (computeMap->type(i) == computeSelfExclsType)
	      || (computeMap->type(i) == computeSelfBondsType)
	      || (computeMap->type(i) == computeSelfAnglesType)
	      || (computeMap->type(i) == computeSelfDihedralsType)
	      || (computeMap->type(i) == computeSelfImpropersType)
	      || (computeMap->type(i) == computeSelfTholeType)
	      || (computeMap->type(i) == computeSelfAnisoType)
	      || (computeMap->type(i) == computeSelfCrosstermsType)

                 || (computeMap->type(i) == computeBondsType)
                 || (computeMap->type(i) == computeExclsType)
                 || (computeMap->type(i) == computeAnglesType)
                 || (computeMap->type(i) == computeDihedralsType)
                 || (computeMap->type(i) == computeImpropersType)
                 || (computeMap->type(i) == computeTholeType)
                 || (computeMap->type(i) == computeAnisoType)
                 || (computeMap->type(i) == computeCrosstermsType)
	      // JLai
	         || (computeMap->type(i) == computeGromacsPairType)
	         || (computeMap->type(i) == computeSelfGromacsPairType)
	) ) {
      nLocalComputes++;
    }
  }
  
  // New LB frameworks registration

  // Allocate data structure to save incoming migrations.  Processor
  // zero will get all migrations

  // If this is the first time through, we need it register patches
  if (ldbCycleNum == reg_all_objs) {
    if ( 1 ) { // ( Node::Object()->simParameters->ldBalancer == LDBAL_CENTRALIZED ) {
      reg_all_objs = 3;
    }
    // Tell the lbdb that I'm registering objects, until I'm done
    // registering them.
    theLbdb->RegisteringObjects(myHandle);
    
   if ( ldbCycleNum == 1 ) {
    patchHandles = new LDObjHandle[nLocalPatches];
    int patch_count=0;
    int i;
    for(i=0;i<nPatches;i++)
      if (patchMap->node(i) == Node::Object()->myid()) {
	LDObjid elemID;
	elemID.id[0] = i;
	elemID.id[1] = elemID.id[2] = elemID.id[3] = -2;

	if (patch_count >= nLocalPatches) {
    NAMD_bug("LdbCoordinator found too many local patches!");
	}
        HomePatch *p = patchMap->homePatch(i);
        p->ldObjHandle = 
	patchHandles[patch_count] 
	  = theLbdb->RegisterObj(myHandle,elemID,0,0);
	patch_count++;

      }
   }
  
    if ( numComputes > oldNumComputes ) {
      // Register computes
      for(i=oldNumComputes; i<numComputes; i++)  {
	if ( computeMap->node(i) == Node::Object()->myid())
        {
	  if ( 0
               #if (defined(NAMD_CUDA) || defined(NAMD_MIC))
                 #if defined(NAMD_MIC)
                   || ((computeMap->type(i) == computeNonbondedSelfType) && (computeMap->directToDevice(i) == 0))
                   || ((computeMap->type(i) == computeNonbondedPairType) && (computeMap->directToDevice(i) == 0))
                 #endif
               #else
	          || (computeMap->type(i) == computeNonbondedSelfType)
	          || (computeMap->type(i) == computeNonbondedPairType)
               #endif
	          || (computeMap->type(i) == computeLCPOType)
	          || (computeMap->type(i) == computeSelfExclsType)
	          || (computeMap->type(i) == computeSelfBondsType)
	          || (computeMap->type(i) == computeSelfAnglesType)
	          || (computeMap->type(i) == computeSelfDihedralsType)
	          || (computeMap->type(i) == computeSelfImpropersType)
	          || (computeMap->type(i) == computeSelfTholeType)
	          || (computeMap->type(i) == computeSelfAnisoType)
	          || (computeMap->type(i) == computeSelfCrosstermsType)
	       // JLai
	          || (computeMap->type(i) == computeSelfGromacsPairType)
	       // End of JLai
		)  {
	  // Register the object with the load balancer
	  // Store the depended patch IDs in the rest of the element ID
	  LDObjid elemID;
	  elemID.id[0] = i;
	
	  if (computeMap->numPids(i) > 2)
	    elemID.id[3] = computeMap->pid(i,2);
	  else elemID.id[3] = -1;

	  if (computeMap->numPids(i) > 1)
	    elemID.id[2] =  computeMap->pid(i,1);
	  else elemID.id[2] = -1;

	  if (computeMap->numPids(i) > 0)
	    elemID.id[1] =  computeMap->pid(i,0);
	  else elemID.id[1] = -1;

          Compute *c = computeMap->compute(i);
          if ( ! c ) NAMD_bug("LdbCoordinator::initialize() null compute pointer");

          c->ldObjHandle = theLbdb->RegisterObj(myHandle,elemID,0,1);
          }
          else if ( (computeMap->type(i) == computeBondsType)
                 || (computeMap->type(i) == computeExclsType)
                 || (computeMap->type(i) == computeAnglesType)
                 || (computeMap->type(i) == computeDihedralsType)
                 || (computeMap->type(i) == computeImpropersType)
                 || (computeMap->type(i) == computeTholeType)
                 || (computeMap->type(i) == computeAnisoType)
                 || (computeMap->type(i) == computeCrosstermsType)
		 // JLai
		 || (computeMap->type(i) == computeGromacsPairType)
                 // End of JLai
               ) {
	  // Register the object with the load balancer
	  // Store the depended patch IDs in the rest of the element ID
	  LDObjid elemID;
	  elemID.id[0] = i;
	
	  elemID.id[1] = elemID.id[2] = elemID.id[3] = -3;

          Compute *c = computeMap->compute(i);
          if ( ! c ) NAMD_bug("LdbCoordinator::initialize() null compute pointer");

          c->ldObjHandle = theLbdb->RegisterObj(myHandle,elemID,0,0);
          }
	}
      }
    }
    theLbdb->DoneRegisteringObjects(myHandle);
  }

  // process saved migration messages, if any
  while ( migrateMsgs ) {
    LdbMigrateMsg *m = migrateMsgs;
    migrateMsgs = m->next;
    Compute *c = computeMap->compute(m->handle.id.id[0]);
    if ( ! c ) NAMD_bug("LdbCoordinator::initialize() null compute pointer 2");
    c->ldObjHandle = m->handle;
    delete m;
  }

  // Fixup to take care of the extra timestep at startup
  // This is pretty ugly here, but it makes the count correct
  
  // iout << "LDB Cycle Num: " << ldbCycleNum << "\n";

 if ( 1 ) { // ( simParams->ldBalancer == LDBAL_CENTRALIZED ) {
  if (ldbCycleNum == 1 || ldbCycleNum == 3) {
    numStepsToRun = stepsPerCycle;
    totalStepsDone += numStepsToRun;
    takingLdbData = 0;
    theLbdb->CollectStatsOff();
  } else if (ldbCycleNum == 2 || ldbCycleNum == 4) {
    numStepsToRun = firstLdbStep - stepsPerCycle;
    while ( numStepsToRun <= 0 ) numStepsToRun += stepsPerCycle;
    totalStepsDone += numStepsToRun;
    takingLdbData = 1;
    theLbdb->CollectStatsOn();
  } else if ( (ldbCycleNum <= 6) || !takingLdbData )
  {
    totalStepsDone += firstLdbStep;
    if(lastLdbStep != -1 && totalStepsDone > lastLdbStep) {
      numStepsToRun = -1;
      takingLdbData = 0;
      theLbdb->CollectStatsOff();
    } else {
      numStepsToRun = firstLdbStep;
      takingLdbData = 1;
      theLbdb->CollectStatsOn();
    }
  }
  else 
  {
    totalStepsDone += stepsPerLdbCycle - firstLdbStep;
    if(lastLdbStep != -1 && totalStepsDone > lastLdbStep) {
      numStepsToRun = -1;
      takingLdbData = 0;
      theLbdb->CollectStatsOff();
    } else {
      numStepsToRun = stepsPerLdbCycle - firstLdbStep;
      takingLdbData = 0;
      theLbdb->CollectStatsOff();
    }
  }
 } else {
  if (ldbCycleNum==1)
  {
    totalStepsDone += firstLdbStep;
    numStepsToRun = firstLdbStep;
    takingLdbData = 0;
    theLbdb->CollectStatsOff();
  }
  else if ( (ldbCycleNum <= 4) || !takingLdbData )
  {
    totalStepsDone += firstLdbStep;
    if(lastLdbStep != -1 && totalStepsDone > lastLdbStep) {
      numStepsToRun = -1;
      takingLdbData = 0;
      theLbdb->CollectStatsOff();
    } else {
      numStepsToRun = firstLdbStep;
      takingLdbData = 1;
      theLbdb->CollectStatsOn();
    }
  }
  else 
  {
    totalStepsDone += stepsPerLdbCycle - firstLdbStep;
    if(lastLdbStep != -1 && totalStepsDone > lastLdbStep) {
      numStepsToRun = -1;
      takingLdbData = 0;
      theLbdb->CollectStatsOff();
    } else {
      numStepsToRun = stepsPerLdbCycle - firstLdbStep;
      takingLdbData = 0;
      theLbdb->CollectStatsOff();
    }
  }
 }

/*-----------------------------------------------------------------------------*
 * --------------------------------------------------------------------------- *
 * Comments inserted by Abhinav to clarify relation between ldbCycleNum,       *
 * load balancing step numbers (printed by the step() function) and            *
 * tracing of the steps                                                        *
 * --------------------------------------------------------------------------- *
 * If trace is turned off in the beginning, then tracing is turned on          *
 * at ldbCycleNum = 4 and turned off at ldbCycleNum = 8. ldbCycleNum can       *
 * be adjusted by specifying firstLdbStep and ldbPeriod which are set by       *
 * default to 5*stepspercycle and 200*stepspercycle if not specified.          *
 *                                                                             *
 * If we choose firstLdbStep = 20 and ldbPeriod = 100, we have the             *
 * following timeline (for these particular numbers):                          *
 *                                                                             *
 * Tracing         :  <------ off ------><------------- on -----------><-- off *
 * Ldb Step() No   :              1     2     3        4      5       6      7 *
 * Iteration Steps : 00====20====40====60====80======160====180=====260====280 *
 * ldbCycleNum     :  1     2     3     4     5        6      7       8      9 *
 * Instrumention   :          Inst  Inst  Inst           Inst            Inst  *
 * LDB Strategy    :              TLB  RLB   RLB            RLB            RLB *
 *                                                                             *
 * TLB = TorusLB                                                               *
 * RLB = RefineTorusLB                                                         *
 * Inst = Instrumentation Phase (no real load balancing)                       *
 * --------------------------------------------------------------------------- *
 *-----------------------------------------------------------------------------*
 */
#if 0 //replaced by traceBarrier at Controller and Sequencer
  if (traceAvailable()) {
    static int specialTracing = 0; // XXX static variables are unsafe for SMP
    if (ldbCycleNum == 1 && traceIsOn() == 0)  specialTracing = 1;
    if (specialTracing) {
      if (ldbCycleNum == 4) traceBegin();
      if (ldbCycleNum == 8) traceEnd();
    }
  }
#endif

  nPatchesReported = 0;
  nPatchesExpected = nLocalPatches;
  nComputesReported = 0;
  nComputesExpected = nLocalComputes * numStepsToRun;
  controllerReported = 0;
  controllerExpected = ! CkMyPe();

  if (simParams->multigratorOn) {
    // Add the number of pressure cycles into nComputesExpected:
    // Pressure cycle is done when !(step % simParams->multigratorPressureFreq) = true
    // step = Current step
    int step = totalStepsDone - numStepsToRun;
    int freq = simParams->multigratorPressureFreq;
    // dstep = Number of steps we have to take until next pressure cycle
    int dstep = 0;
    if ((step % freq) != 0) dstep = freq - (step % freq);
    step += dstep;
    if (step < totalStepsDone) {
      int numPressureCycles = 1 + ((totalStepsDone-step-1)/freq);
      if (step==0) numPressureCycles--;
      // if (CkMyPe()==2) fprintf(stderr, "step %d totalStepsDone %d numPressureCycles %d\n",
      //   step, totalStepsDone, numPressureCycles);
      nComputesExpected += 2*nLocalComputes*numPressureCycles;
    }
  }

  if (CkMyPe() == 0)
  {
    if (computeArray == NULL)
      computeArray = new computeInfo[numComputes];
    if (patchArray == NULL)
      patchArray = new patchInfo[nPatches];
    if (processorArray == NULL)
      processorArray = new processorInfo[CkNumPes()];
  }
    
  theLbdb->ClearLoads();
}

void LdbCoordinator::patchLoad(PatchID id, int nAtoms, int /* timestep */)
{
  CmiAssert( id >=0 && id < nPatches);
  if (patchNAtoms[id] != -1) {
    patchNAtoms[id] = nAtoms;
    nPatchesReported++;
  } else {
    DebugM(10, "::patchLoad() Unexpected patch reporting in\n");
  }
}

void LdbCoordinator::rebalance(Sequencer *seq, PatchID pid)
{
  if (Node::Object()->simParameters->ldBalancer == LDBAL_NONE)
    return;

  sequencerThreads[pid] = seq;
  seq->suspend();
}

void LdbCoordinator::rebalance(Controller *c)
{
  if (Node::Object()->simParameters->ldBalancer == LDBAL_NONE)
    return;

  iout << "LDB: ============= START OF LOAD BALANCING ============== " << CmiWallTimer() << "\n" << endi;
  DebugM(3, "Controller reached load balance barrier.\n");
  controllerReported = 1;
  controllerThread = c;

  CProxy_LdbCoordinator(thisgroup).barrier();

  CthSuspend();
}

void LdbCoordinator::barrier(void)
{
  if ( (nPatchesReported != nPatchesExpected) 
       || (nComputesReported != nComputesExpected)
       || (controllerReported != controllerExpected) )
  {
    NAMD_bug("Load balancer received wrong number of events.\n");
  }

  theLbdb->AtLocalBarrier(ldBarrierHandle);
}

void LdbCoordinator::nodeDone(CkReductionMsg *msg)
{
  delete msg;

  iout << "LDB: ============== END OF LOAD BALANCING =============== " << CmiWallTimer() << "\n" << endi;
  if ( takingLdbData ) {
      ExecuteMigrations();
  } else {
      updateComputesReady();
  }
}

void LdbCoordinator::ExecuteMigrations(void)
{
 // computeMgr->updateComputes() call only on Node(0) i.e. right here
  // This will barrier for all Nodes - (i.e. Computes must be
  // here and with proxies before anyone can start up

  CProxy_ComputeMgr cm(CkpvAccess(BOCclass_group).computeMgr);
  ComputeMgr *computeMgr = cm.ckLocalBranch();
  computeMgr->updateComputes(CkIndex_LdbCoordinator::
                             updateComputesReady(),thisgroup);
}

void LdbCoordinator::RecvMigrate(LdbMigrateMsg* m)
{
  // This method receives the migration from the framework,
  // unregisters it, and sends it to the destination PE

  if ( m->to != CkMyPe() ) {
    theLbdb->UnregisterObj(m->handle);

    CProxy_LdbCoordinator  ldbProxy(thisgroup);
    ldbProxy[m->to].ExpectMigrate(m);
  } else {
    ExpectMigrate(m);
  }
}

void LdbCoordinator::ExpectMigrate(LdbMigrateMsg* m)
{
  if ( m->from != CkMyPe() ) {
    m->handle = theLbdb->RegisterObj(myHandle,m->handle.id,0,1);
    theLbdb->Migrated(m->handle);
  }

  m->next = migrateMsgs;
  migrateMsgs = m;
}

void LdbCoordinator::updateComputesReady() {
  DebugM(3,"updateComputesReady()\n");

  CProxy_LdbCoordinator(thisgroup).resume();
  CkStartQD(CkIndex_LdbCoordinator::resumeReady((CkQdMsg*)0),&thishandle);
}

void LdbCoordinator::resume(void)
{
  DebugM(3,"resume()\n");
  //  printLocalLdbReport();

  ldbCycleNum++;
  initialize(PatchMap::Object(),ComputeMap::Object(),1);

  Sync::Object()->openSync();
}

void LdbCoordinator::resumeReady(CkQdMsg *msg) {

  iout << "LDB: =============== DONE WITH MIGRATION ================ " << CmiWallTimer() << "\n" << endi;
  DebugM(3,"resumeReady()\n");
  delete msg;

  CProxy_LdbCoordinator(thisgroup).resume2();
}

void LdbCoordinator::resume2(void)
{
  DebugM(3,"resume2()\n");

#if CONVERSE_VERSION_ELAN
  //  enableBlockingReceives();
#endif

  awakenSequencers();
}

void LdbCoordinator::awakenSequencers()
{
  if (controllerThread)
  {
    controllerThread->awaken();
    controllerThread = NULL;
  }
  for(int i=0; i < patchMap->numPatches(); i++)
  {
    if (sequencerThreads[i])
    {
      sequencerThreads[i]->awaken();
    }
    sequencerThreads[i]= NULL;
  }
}

// Figure out which proxies we will definitely create on other
// nodes, without regard for non-bonded computes.  This code is swiped
// from ProxyMgr, and changes there probable need to be propagated here.

int LdbCoordinator::requiredProxies(PatchID id, int neighborNodes[])
{
  PatchID neighbors[1 + PatchMap::MaxOneAway];
  neighbors[0] = id;
  int numNeighbors = 1 + patchMap->downstreamNeighbors(id,neighbors+1);

  int nProxyNodes = 0;
  int myNode = patchMap->node(id);
  for ( int i = 0; i < numNeighbors; ++i ) {
    const int proxyNode = patchMap->basenode(neighbors[i]);
    if ( proxyNode != myNode ) {
      int j;
      for ( j = 0; j < nProxyNodes; ++j ) {
        if ( neighborNodes[j] == proxyNode ) break;
      }
      if ( j == nProxyNodes ) {
        neighborNodes[nProxyNodes] = proxyNode;
        nProxyNodes++;
      }
    }
  }
  return nProxyNodes;
}

void LdbCoordinator::printLocalLdbReport(void)
{
  char outputBuf[255];
  char *curLoc;

  CkPrintf("%d:Patch report:\n",CkMyPe());
  
  curLoc = outputBuf;
  int i,j=0;
  for(i=0; i<patchMap->numPatches(); i++)
  {
    if (patchNAtoms[i] != -1)
    {
      curLoc += sprintf(curLoc,"%5d: %5d ",i,patchNAtoms[i]);
      j++;
    } 
    if (((j % 4) == 0) && j)
    {
      curLoc = outputBuf;
      CkPrintf("[%d]%s\n",CkMyPe(),outputBuf);
      j=0;
    }
  }

  CkPrintf("%d:Compute report:\n",CkMyPe());
  
  curLoc = outputBuf;
  j=0;
}

void LdbCoordinator::printRequiredProxies(PatchID id, FILE *fp)
{
  // Check all two-away neighbors.
  // This is really just one-away neighbors, since 
  // two-away always returns zero: RKB
  int neighborNodes[PatchMap::MaxOneAway + PatchMap::MaxTwoAway];
  const int nProxyNodes = requiredProxies(id,neighborNodes);

  fprintf(fp,"%4d ",nProxyNodes);

  for(int i=0;i<nProxyNodes;i++)
    fprintf(fp,"%4d ",neighborNodes[i]);
}

void LdbCoordinator::sendCollectLoads(CollectLoadsMsg *msg) {
  CProxy_LdbCoordinator(thisgroup)[0].collectLoads(msg);
}

void LdbCoordinator::collectLoads(CollectLoadsMsg *msg) {
  // CkPrintf("LdbCoordinator::collectLoads recv %d-%d\n", msg->firstPe, msg->lastPe);
  if ( collPes == 0 ) {
    reverted = 0;
    initTotalProxies = 0;
    finalTotalProxies = 0;
    initMaxPeProxies = 0;
    finalMaxPeProxies = 0;
    initMaxPatchProxies = 0;
    finalMaxPatchProxies = 0;
    initTime = 0;
    finalTime = 0;
    initMemory = 0;
    finalMemory = 0;
    initAvgPeLoad = 0;
    finalAvgPeLoad = 0;
    initMaxPeLoad = 0;
    finalMaxPeLoad = 0;
  }
  int numPes = msg->lastPe - msg->firstPe + 1;
  collPes += numPes;
#define COLL_MAX(F) if ( msg->F > F ) F = msg->F;
#define COLL_AVG(F) F += msg->F * (double) numPes / (double) CkNumPes();
#define COLL_SUM(F) F += msg->F;
  COLL_SUM(reverted)
  COLL_SUM(initTotalProxies)
  COLL_SUM(finalTotalProxies)
  COLL_MAX(initMaxPeProxies)
  COLL_MAX(finalMaxPeProxies)
  COLL_MAX(initMaxPatchProxies)
  COLL_MAX(finalMaxPatchProxies)
  if ( (msg->finalTime - msg->initTime) > (finalTime - initTime) ) {
    initTime = msg->initTime;
    finalTime = msg->finalTime;
  }
  COLL_MAX(initMemory)
  COLL_MAX(finalMemory)
  COLL_AVG(initAvgPeLoad)
  COLL_AVG(finalAvgPeLoad)
  COLL_MAX(initMaxPeLoad)
  COLL_MAX(finalMaxPeLoad)

  if ( collPes == CkNumPes() ) {
    collPes = 0;
    iout << "LDB: TIME " << initTime << " LOAD: AVG " << initAvgPeLoad
      << " MAX " << initMaxPeLoad << "  PROXIES: TOTAL " << initTotalProxies << " MAXPE " <<
      initMaxPeProxies << " MAXPATCH " << initMaxPatchProxies << " " << "None"
      << " MEM: " << initMemory << " MB\n";
    if ( reverted ) iout << "LDB: Reverting to original mapping on " << reverted << " balancers\n";
    iout << "LDB: TIME " << finalTime << " LOAD: AVG " << finalAvgPeLoad
      << " MAX " << finalMaxPeLoad << "  PROXIES: TOTAL " << finalTotalProxies << " MAXPE " <<
      finalMaxPeProxies << " MAXPATCH " << finalMaxPatchProxies << " " << msg->strategyName
      << " MEM: " << finalMemory << " MB\n";
    iout << endi;
    fflush(stdout);
  }

  delete msg;
}

#include "LdbCoordinator.def.h"
