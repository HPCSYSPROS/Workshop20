/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Toplevel routines for initializing a Node for a simulation
   one Node per Pe (processor element).
*/

#if !defined(WIN32) || defined(__CYGWIN__)
#include <unistd.h>
#endif
#include "InfoStream.h"
#include "Node.decl.h"
#include "Node.h"
#ifdef DPMTA
#include <pvm3.h>
#endif

#include "ProcessorPrivate.h"

#define MIN_DEBUG_LEVEL 3
//#define DEBUGM
#include "Debug.h"

#include <stdio.h>
#include <converse.h>
#include "memusage.h"
#include "IMDOutput.h"
#include "Lattice.h"
#include "ComputeMsmMsa.h"  // needed for MsmMsaData definition
#include "ComputeMsm.h"     // needed for MsmInitMsg definition
#include "main.decl.h"
#include "main.h"
#include "WorkDistrib.h"
#include "PatchMgr.h"
#include "Patch.h"
#include "Compute.h"
#include "ComputeMap.h"
#include "ComputeMgr.h"
#include "Molecule.h"
#include "HomePatchList.h"
#include "AtomMap.h"
#include "Sequencer.h"
#include "Controller.h"
#include "NamdState.h"
#include "Output.h"
#include "ProxyMgr.h"
#include "PatchMap.h"
#include "PatchMap.inl"
#include "Parameters.h"
#include "SimParameters.h"
#include "Communicate.h"
#include "LdbCoordinator.h"
#include "ScriptTcl.h"
#include "ComputeMgr.decl.h"
#include "ComputePmeMgr.decl.h"
// #ifdef NAMD_CUDA
#include "ComputeCUDAMgr.decl.h"
#include "ComputePmeCUDAMgr.decl.h"
// #endif
#include "ComputeGridForceMgr.decl.h"
#include "OptPmeMgr.decl.h"
#include "Sync.h"
#include "BackEnd.h"
#include "PDB.h"
#include "packmsg.h"
#include "CollectionMgr.decl.h"
#include "ParallelIOMgr.decl.h"
#include "Vector.h"
// BEGIN LA
#include "Random.h"
// END LA

#if(CMK_CCS_AVAILABLE && CMK_WEB_MODE)
extern "C" void CApplicationInit();
#endif

#include "DumpBench.h"

class CheckpointMsg : public CMessage_CheckpointMsg {
public:
  int task;
  int replica;
  Controller::checkpoint checkpoint;
  char *key;
};

extern "C" {
  void recvCheckpointCReq_handler(envelope*);
  void recvCheckpointCAck_handler(envelope*);
}

#ifdef CMK_BALANCED_INJECTION_API
#include "ckBIconfig.h"
#endif

#include "CollectionMgr.h"
#include "CollectionMaster.h"
#include "CollectionMgr.decl.h"
#include "CollectionMaster.decl.h"

#if USE_HPM
extern "C" void HPM_Init(int);
extern "C" void HPM_Start(char *label, int);
extern "C" void HPM_Stop(char *label, int);
extern "C" void HPM_Print(int, int);
#endif

#if defined(NAMD_MIC)
  extern void mic_dumpHostDeviceComputeMap();
  extern void mic_initHostDeviceLDB();
#endif

#ifdef MEASURE_NAMD_WITH_PAPI
#include "papi.h"
#if CMK_SMP
#include <pthread.h>
#endif
#define NUM_PAPI_EVENTS 6
CkpvDeclare(int *, papiEvents);

#define MEASURE_PAPI_SPP 1
#define MEASURE_PAPI_CACHE 0
#define MEASURE_PAPI_FLOPS 0

static void namdInitPapiCounters(){
	if(CkMyRank()==0){
		//only initialize per OS process (i.e. a charm node)
		int retval = PAPI_library_init(PAPI_VER_CURRENT);
		if(retval != PAPI_VER_CURRENT) {
			if(CkMyPe()==0){
				NAMD_die("PAPI library is not compatitible!");
			}
		}
	#if CMK_SMP
		//now only consider systems that are compatible with POSIX
		if(PAPI_thread_init(pthread_self)!=PAPI_OK) {
			if(CkMyPe()==0){
				NAMD_die("Multi-thread mode in PAPI could not be initialized!");
			}
		}
	#endif
	}
	CkpvInitialize(int *, papiEvents);
	CkpvAccess(papiEvents) = new int[NUM_PAPI_EVENTS+1];

#if MEASURE_PAPI_CACHE
	if(PAPI_query_event(PAPI_L1_DCM)==PAPI_OK) {
		CkpvAccess(papiEvents)[0] = PAPI_L1_DCM;
	}else{
		if(CkMyPe()==0){
			CkPrintf("WARNING: PAPI_L1_DCM doesn't exsit on this platform!\n");			
		}
		//if not default to PAPI_TOT_INS
		CkpvAccess(papiEvents)[0] = PAPI_TOT_INS;
	}

	if(PAPI_query_event(PAPI_L2_DCM)==PAPI_OK) {
		CkpvAccess(papiEvents)[1] = PAPI_L2_DCM;
	}else{
		//if not default to PAPI_TOT_CYC
		CkpvAccess(papiEvents)[1] = PAPI_TOT_CYC;
	}	
#elif MEASURE_PAPI_FLOPS
	if(PAPI_query_event(PAPI_FP_INS)==PAPI_OK) {
		CkpvAccess(papiEvents)[0] = PAPI_FP_INS;
	}else{
		if(CkMyPe()==0){
			CkPrintf("WARNING: PAPI_FP_INS doesn't exsit on this platform!\n");
		}
		//if not default to PAPI_TOT_INS
		CkpvAccess(papiEvents)[0] = PAPI_TOT_INS;
	}

	if(PAPI_query_event(PAPI_FMA_INS)==PAPI_OK) {
		CkpvAccess(papiEvents)[1] = PAPI_FMA_INS;
	}else{
		//if not default to PAPI_TOT_CYC
		CkpvAccess(papiEvents)[1] = PAPI_TOT_CYC;
	}
#elif MEASURE_PAPI_SPP
/* for SPP we record these
1) PAPI_FP_OPS
2) PAPI_TOT_INS
3) perf::PERF_COUNT_HW_CACHE_LL:MISS
4) DATA_PREFETCHER:ALL
5) PAPI_L1_DCA
6) INSTRUCTION_FETCH_STALL
7) PAPI_TOT_CYC, and 
8) real (wall) time
*/
	int papiEventSet = PAPI_NULL; 
	if (PAPI_create_eventset(&papiEventSet) != PAPI_OK) {
	  CmiAbort("PAPI failed to create event set!\n");
	}

	if(PAPI_query_event(PAPI_FP_OPS)==PAPI_OK) {
		CkpvAccess(papiEvents)[0] = PAPI_FP_OPS;
	}else{
		if(CkMyPe()==0){
			CkAbort("WARNING: PAPI_FP_OPS doesn't exist on this platform!");
		}
	}
	if(PAPI_query_event(PAPI_TOT_INS)==PAPI_OK) {
		CkpvAccess(papiEvents)[1] = PAPI_TOT_INS;
	}else{
		if(CkMyPe()==0){
			CkAbort("WARNING: PAPI_TOT_INS doesn't exist on this platform!");
		}
	}
	int EventCode;
	int ret;
	ret=PAPI_event_name_to_code("perf::PERF_COUNT_HW_CACHE_LL:MISS",&EventCode);
	if(ret==PAPI_OK && PAPI_query_event(EventCode)==PAPI_OK) {
	  CkpvAccess(papiEvents)[2] = EventCode;
	}else{
		if(CkMyPe()==0){
			CkAbort("WARNING: perf::PERF_COUNT_HW_CACHE_LL:MISS doesn't exist on this platform!");
		}
	}
	ret=PAPI_event_name_to_code("DATA_PREFETCHER:ALL",&EventCode);
	if(ret==PAPI_OK && PAPI_query_event(EventCode)==PAPI_OK) {
	  CkpvAccess(papiEvents)[3] = EventCode;
	}else{
		if(CkMyPe()==0){
			CkAbort("WARNING: DATA_PREFETCHER:ALL doesn't exist on this platform!");
		}
	}
	if(PAPI_query_event(PAPI_L1_DCA)==PAPI_OK) {
		CkpvAccess(papiEvents)[4] = PAPI_L1_DCA;
	}else{
		if(CkMyPe()==0){
			CkAbort("WARNING: PAPI_L1_DCA doesn't exist on this platform!");
		}
	}
	/*	ret=PAPI_event_name_to_code("INSTRUCTION_FETCH_STALL",&EventCode);
	if(ret==PAPI_OK && PAPI_query_event(EventCode)==PAPI_OK) {
	  CkpvAccess(papiEvents)[5] = EventCode;
	}else{
		if(CkMyPe()==0){
			CkAbort("WARNING: INSTRUCTION_FETCH_STALL doesn't exist on this platform!");
		}
	}
	*/
	if(PAPI_query_event(PAPI_TOT_CYC)==PAPI_OK) {
		CkpvAccess(papiEvents)[5] = PAPI_TOT_CYC;
	}else{
		if(CkMyPe()==0){
			CkAbort("WARNING: PAPI_TOT_CYC doesn't exist on this platform!");
		}
	}
	for(int i=0;i<NUM_PAPI_EVENTS;i++)
	  {
	    int papiRetValue=PAPI_add_events(papiEventSet, &CkpvAccess(papiEvents)[i],1);
	    if (papiRetValue != PAPI_OK) {
	      CkPrintf("failure for event %d\n",i);
	      if (papiRetValue == PAPI_ECNFLCT) {
		CmiAbort("PAPI events conflict! Please re-assign event types!\n");
	      } else {
		CmiAbort("PAPI failed to add designated events!\n");
	      }
	    }
	    
	  }
#endif
}
#endif

#ifdef OPENATOM_VERSION
static void startOA(){(char inDriverFile[1024], char inPhysicsFile[1024], CkCallback doneCB)
{
  CProxy_oaSetup moaInstance = CProxy_oaSetup::ckNew(inDriverFile, inPhysicsFile, doneCB);
}
#endif //OPENATOM_VERSION

//======================================================================
// Public Functions

//----------------------------------------------------------------------

int eventEndOfTimeStep;
double startupTime;

//----------------------------------------------------------------------
// BOC constructor
Node::Node(GroupInitMsg *msg)
{    
  DebugM(4,"Creating Node\n");
#if(CMK_CCS_AVAILABLE && CMK_WEB_MODE)
  CApplicationInit();
#endif
  if (CkpvAccess(Node_instance) == 0) {
    CkpvAccess(Node_instance) = this;
    eventEndOfTimeStep = traceRegisterUserEvent("EndOfTimeStep");
  } else {
    NAMD_bug("Node::Node() - another instance of Node exists!");
  }

  CkpvAccess(BOCclass_group) = msg->group;
  delete msg;

  CkpvAccess(BOCclass_group).node = thisgroup;

  recvCheckpointCReq_index = CmiRegisterHandler((CmiHandler)recvCheckpointCReq_handler);
  recvCheckpointCAck_index = CmiRegisterHandler((CmiHandler)recvCheckpointCAck_handler);

  startupPhase = 0;

  molecule = NULL;
  parameters = NULL;
  simParameters = NULL;
  configList = NULL;
  pdb = NULL;
  state = NULL;
  output = NULL;
  imd = new IMDOutput;
  colvars = 0;

#if USE_HPM
  // assumes that this will be done only on BG/P
  TopoManager *tmgr = new TopoManager();
  int x, y, z;
  tmgr->rankToCoordinates(CkMyPe(), x, y, z, localRankOnNode);
  delete tmgr;
#endif

  specialTracing = traceAvailable() && (traceIsOn()==0);

  DebugM(4,"Creating PatchMap, AtomMap, ComputeMap\n");
  patchMap = PatchMap::Instance();
  atomMap = AtomMap::Instance();
  if ( CkMyRank() == 0 ) ComputeMap::Instance();

  //Note: Binding BOC vars such as workDistrib has been moved
  //to the 1st phase of startup because the in-order message delivery
  //is not always guaranteed --Chao Mei
#ifdef CMK_BALANCED_INJECTION_API
  if(CkMyRank() == 0){
    balancedInjectionLevel=ck_get_GNI_BIConfig();
    // CkPrintf("[%d] get retrieved BI=%d\n",CkMyPe(),balancedInjectionLevel);
    ck_set_GNI_BIConfig(20);
    // CkPrintf("[%d] set retrieved BI=%d\n",CkMyPe(),ck_get_GNI_BIConfig());
  }
#endif

}

//----------------------------------------------------------------------
// ~Node(void) needs to clean up everything.

Node::~Node(void)
{
  delete output;
  delete computeMap;
  delete atomMap;
  delete patchMap;
  delete CkpvAccess(comm);
  // BEGIN LA
  delete rand;
  // END LA
#ifdef MEASURE_NAMD_WITH_PAPI
  delete CkpvAccess(papiEvents);
#endif
}

void Node::bindBocVars(){
    DebugM(4,"Binding to BOC's\n");
    CProxy_PatchMgr pm(CkpvAccess(BOCclass_group).patchMgr);
    patchMgr = pm.ckLocalBranch();
    CProxy_ProxyMgr prm(CkpvAccess(BOCclass_group).proxyMgr);
    proxyMgr = prm.ckLocalBranch();
    CProxy_WorkDistrib wd(CkpvAccess(BOCclass_group).workDistrib);
    workDistrib = wd.ckLocalBranch();
    CProxy_ComputeMgr cm(CkpvAccess(BOCclass_group).computeMgr);
    computeMgr = cm.ckLocalBranch();
    CProxy_LdbCoordinator lc(CkpvAccess(BOCclass_group).ldbCoordinator);
    ldbCoordinator = lc.ckLocalBranch();
  #ifdef MEM_OPT_VERSION      
    CProxy_ParallelIOMgr io(CkpvAccess(BOCclass_group).ioMgr);
    ioMgr = io.ckLocalBranch();
  #endif

}

//----------------------------------------------------------------------
// Malloc Test Sequence
void Node::mallocTest(int step) {
  int MB = 1024*1024;
  int size = 100;
  char* foo = (char*) malloc(size*MB);
  if ( ! foo ) {
    char buf[256];
    sprintf(buf,"Malloc fails on Pe %d at %d MB.\n",CkMyPe(),step*size);
    NAMD_die(buf);
  }
  memset(foo,0,size*MB*sizeof(char));
}

void Node::mallocTestQd(CkQdMsg *qmsg) {
  delete qmsg;
  if ( mallocTest_size ) {
    CkPrintf("All PEs successfully allocated %d MB.\n", 100*mallocTest_size);
  } else {
    CkPrintf("Starting malloc test on all PEs.\n");
  }
  fflush(stdout);
  ++mallocTest_size;
  CkStartQD(CkIndex_Node::mallocTestQd((CkQdMsg*)0),&thishandle);
  (CProxy_Node(CkpvAccess(BOCclass_group).node)).mallocTest(mallocTest_size);
}

//----------------------------------------------------------------------
// Startup Sequence

void Node::messageStartUp() {
  (CProxy_Node(CkpvAccess(BOCclass_group).node)).startup();
}

void Node::startUp(CkQdMsg *qmsg) {
  delete qmsg;
  (CProxy_Node(CkpvAccess(BOCclass_group).node)).startup();
}

SimParameters *node_simParameters;
Parameters *node_parameters;
Molecule *node_molecule;

extern void registerUserEventsForAllComputeObjs(void);

void Node::startup() {
  int gotoRun = false;
  double newTime;

  if (!CkMyPe()) {
    if (!startupPhase) {
      iout << iINFO << "\n";
      startupTime = CmiWallTimer();
      iout << iINFO << "Entering startup at " << startupTime << " s, ";
    } else {
      newTime = CmiWallTimer();
      iout << iINFO << "Startup phase " << startupPhase-1 << " took "
	   << newTime - startupTime << " s, ";
      startupTime = newTime;
    }
    iout << memusage_MB() << " MB of memory in use\n" << endi;
    fflush(stdout);
  }
  switch (startupPhase) {

  case 0:
    computeMap = ComputeMap::Object();
    namdOneCommInit(); // Namd1.X style
  break;

  case 1:
      bindBocVars();

    // send & receive molecule, simparameters... (Namd1.X style)
    if (CkMyPe()) {
      namdOneRecv();
    } else {
      namdOneSend();
    }
  break;

  case 2:
    // fix up one-per-node objects (for SMP version)
    simParameters = node_simParameters;
    parameters = node_parameters;
    molecule = node_molecule;

    SimParameters::nonbonded_select();
    
    #if !CMK_SMP || ! USE_CKLOOP
    //the CkLoop library should be only used in SMP mode
    simParameters->useCkLoop = 0;
    #else
    if ( CkNumPes() < 2 * CkNumNodes() ) simParameters->useCkLoop = 0;
    #endif


    if ( simParameters->mallocTest ) {
      if (!CkMyPe()) {
        mallocTest_size = 0;
        CkStartQD(CkIndex_Node::mallocTestQd((CkQdMsg*)0),&thishandle);
      }
      return;
    }

      
	#ifdef MEASURE_NAMD_WITH_PAPI
	if(simParameters->papiMeasure) namdInitPapiCounters();	
	#endif
    
    #ifdef MEM_OPT_VERSION
    //At this point, each Node object has received the simParameters,
    //parameters and the atom signatures info from the master Node
    //(proc 0). It's time to initialize the parallel IO manager and
    //read the binary per-atom file --Chao Mei

    //Step 1: initialize the parallel IO manager per Node
    ioMgr->initialize(this);
    #endif

  break;

  case 3:

    #ifdef MEM_OPT_VERSION
    //Step 2: read the binary per-atom files (signater index, coordinates etc.)
    ioMgr->readPerAtomInfo();
    #endif

  break;

  case 4:

    #ifdef MEM_OPT_VERSION
    //Step 3: update counters of tuples and exclusions inside Molecule object
    ioMgr->updateMolInfo();

    //Step 4: prepare distributing the atoms to neighboring procs if necessary
    ioMgr->migrateAtomsMGrp();

    //step 5: initialize patchMap and send it to every other processors
    //to decide atoms to patch distribution on every input processor
    if(!CkMyPe()) {
        workDistrib->patchMapInit(); // create space division
        workDistrib->sendPatchMap();
    }
    #endif

    #if USE_HPM
    HPM_Init(localRankOnNode);
    #endif    

    // take care of inital thread setting
    threadInit();

    // create blank AtomMap
    AtomMap::Object()->allocateMap(molecule->numAtoms);

    if (!CkMyPe()) {
      if (simParameters->useOptPME)
	CkpvAccess(BOCclass_group).computePmeMgr = CProxy_OptPmeMgr::ckNew();
      else 
#ifdef NAMD_CUDA
      if (simParameters->usePMECUDA) {
        // computePmeCUDAMgr was created in BackEnd.C
        // This empty branch is to avoid initializing ComputePmeMgr
      } else
#endif
      if (simParameters->PMEOn) {
        CkpvAccess(BOCclass_group).computePmeMgr = CProxy_ComputePmeMgr::ckNew();
      }
        #ifdef OPENATOM_VERSION
        if ( simParameters->openatomOn ) { 
          CkpvAccess(BOCclass_group).computeMoaMgr = CProxy_ComputeMoaMgr::ckNew();
        }
        #endif // OPENATOM_VERSION

    }
    
    #ifdef OPENATOM_VERSION
    if ( simParameters->openatomOn ) {
      // if ( ! CkMyPe() ) { 
        CkCallback doneMoaStart(CkIndexmain::doneMoaSetup(), thishandle); 
        startOA(simParameters->moaDriverFile, simParameters->moaPhysicsFile, doneMoaStart);
      // }
    }
    #endif // OPENATOM_VERSION
  
    // BEGIN LA
    rand = new Random(simParameters->randomSeed);
    rand->split(CkMyPe(), CkNumPes());
    // END LA

  break;

  case 5:
    #ifdef MEM_OPT_VERSION
    //Now, every input proc has received all the atoms necessary
    //to decide the patches those atoms belong to
    
    //step 1: integrate the migrated atoms into the atom list that
    //contains the initally distributed atoms, and sort the atoms
    //based on hydrogenList value
    ioMgr->integrateMigratedAtoms();

    //step 2: integrate the cluster size of each atom on each output proc
    ioMgr->integrateClusterSize();

    //step 3: calculate the number of atoms in each patch on every
    //input procs (atoms belonging to a patch may lie on different
    //procs), and reduce such info on proc 0. Such info is required
    //for determing which node a particular patch is assigned to.
    ioMgr->calcAtomsInEachPatch();

    //set to false to re-send PatchMap later
    workDistrib->setPatchMapArrived(false);
    #endif
    break;
  case 6:     
    if(simParameters->isSendSpanningTreeOn()) {				
			ProxyMgr::Object()->setSendSpanning();
    }
    if(simParameters->isRecvSpanningTreeOn()) {				
			ProxyMgr::Object()->setRecvSpanning();
    }
    if(simParameters->proxyTreeBranchFactor) {
			ProxyMgr::Object()->setProxyTreeBranchFactor(simParameters->proxyTreeBranchFactor);
    }
    #ifdef PROCTRACE_DEBUG
    DebugFileTrace::Instance("procTrace");
    #endif

    if (!CkMyPe()) {
      output = new Output; // create output object just on PE(0)

      #ifndef MEM_OPT_VERSION
      workDistrib->patchMapInit(); // create space division
      workDistrib->createHomePatches(); // load atoms into HomePatch(es)
      #endif
      
      workDistrib->assignNodeToPatch();
      workDistrib->mapComputes();	  
      //ComputeMap::Object()->printComputeMap();

      // For MIC runs, take the additional step after the compute map has been created to
      //   assign the various computes to either the host or the device.  This info will
      //   be distributed across the PEs.
      #if defined(NAMD_MIC)
        mic_initHostDeviceLDB();
      #endif
	  
	  if(simParameters->simulateInitialMapping) {
    	  iout << iINFO << "Simulating initial mapping with " << simParameters->simulatedPEs
			  << " PEs with " << simParameters->simulatedNodeSize << " PEs per node\n" << endi;
		  outputPatchComputeMaps("init_mapping", 0);
		  iout << iINFO << "Simulating initial mapping is done, now NAMD exits\n" << endi;
		  BackEnd::exit();
	  }

      registerUserEventsForAllComputeObjs();

      //in MEM_OPT_VERSION, patchMap is resent
      //because they have been updated since creation including
      //#atoms per patch, the proc a patch should stay etc. --Chao Mei
      workDistrib->sendPatchMap();
      #if defined(NODEAWARE_PROXY_SPANNINGTREE) && defined(USE_NODEPATCHMGR)
      CProxy_NodeProxyMgr npm(CkpvAccess(BOCclass_group).nodeProxyMgr);
      //a node broadcast
      npm.createProxyInfo(PatchMap::Object()->numPatches());
      #endif
    }
    {
        #if defined(NODEAWARE_PROXY_SPANNINGTREE) && defined(USE_NODEPATCHMGR)
        CProxy_NodeProxyMgr npm(CkpvAccess(BOCclass_group).nodeProxyMgr);
        if(CkMyRank()==0) {
            //just need to register once
            npm[CkMyNode()].ckLocalBranch()->registerLocalProxyMgr(CkpvAccess(BOCclass_group).proxyMgr);
        }
        npm[CkMyNode()].ckLocalBranch()->registerLocalPatchMap(CkMyRank(), PatchMap::Object());
        #endif
    }
  break;

  case 7:
    if ( simParameters->PMEOn ) {
      if ( simParameters->useOptPME ) {
	CProxy_OptPmeMgr pme(CkpvAccess(BOCclass_group).computePmeMgr);
	pme[CkMyPe()].initialize(new CkQdMsg);
      }
      else {
        #ifdef OPENATOM_VERSION
        if ( simParameters->openatomOn ) { 
          CProxy_ComputeMoaMgr moa(CkpvAccess(BOCclass_group).computeMoaMgr); 
          moa[CkMyPe()].initialize(new CkQdMsg);
        }
        #endif // OPENATOM_VERSION
#ifdef NAMD_CUDA
        if ( simParameters->usePMECUDA ) {
          if(CkMyRank()==0) {
            CProxy_ComputePmeCUDAMgr pme(CkpvAccess(BOCclass_group).computePmeCUDAMgr);
            pme[CkMyNode()].initialize(new CkQdMsg);
          }
        } else 
#endif
        {
          CProxy_ComputePmeMgr pme(CkpvAccess(BOCclass_group).computePmeMgr);
          pme[CkMyPe()].initialize(new CkQdMsg);          
        }
      }
    }
#ifdef NAMD_CUDA
    if ( simParameters->useCUDA2 && CkMyRank()==0 ) {
      CProxy_ComputeCUDAMgr nb(CkpvAccess(BOCclass_group).computeCUDAMgr);
      nb[CkMyNode()].initialize(new CkQdMsg);
    }
#endif

#ifdef CHARM_HAS_MSA
    else if ( simParameters->MSMOn && ! simParameters->MsmSerialOn ) {
      CProxy_ComputeMsmMsaMgr msm(CkpvAccess(BOCclass_group).computeMsmMsaMgr);
      msm[CkMyPe()].initialize(new CkQdMsg);
    }
#else
    else if ( simParameters->MSMOn && ! simParameters->MsmSerialOn ) {
      CProxy_ComputeMsmMgr msm(CkpvAccess(BOCclass_group).computeMsmMgr);
      MsmInitMsg *msg = new MsmInitMsg;
      Lattice lattice = simParameters->lattice;  // system lattice vectors
      ScaledPosition smin=0, smax=0;
      if (lattice.a_p() && lattice.b_p() && lattice.c_p()) {
        msg->smin = smin;
        msg->smax = smax;
        msm[CkMyPe()].initialize(msg);  // call from my own PE
      }
      else if ( ! CkMyPe() ) {
        pdb->get_extremes(smin, smax);  // only available on PE 0
        msg->smin = smin;
        msg->smax = smax;
        msm.initialize(msg);  // broadcast to chare group
      }

      /*
      CProxy_Node nd(CkpvAccess(BOCclass_group).node);
      Node *node = nd.ckLocalBranch();
      ScaledPosition smin, smax;
      node->pdb->get_extremes(smin, smax);
      msg->smin = smin;                       // extreme positions in system
      msg->smax = smax;
      msm[CkMyPe()].initialize(msg);
      */
    }
#endif

    if (!CkMyPe()) {
      workDistrib->sendComputeMap();
    }

    #ifdef MEM_OPT_VERSION
    //migrate atoms to HomePatch processors
    ioMgr->sendAtomsToHomePatchProcs();
    #endif
    break;
    
  case 8:
    if ( simParameters->PMEOn ) {
      if ( simParameters->useOptPME ) {
	CProxy_OptPmeMgr pme(CkpvAccess(BOCclass_group).computePmeMgr);
	pme[CkMyPe()].initialize_pencils(new CkQdMsg);
      }
      else {
        #ifdef OPENATOM_VERSION
        if ( simParameters->openatomOn ) { 
          CProxy_ComputeMoaMgr moa(CkpvAccess(BOCclass_group).computeMoaMgr); 
          moa[CkMyPe()].initWorkers(new CkQdMsg);
        }
        #endif // OPENATOM_VERSION
#ifdef NAMD_CUDA
        if ( simParameters->usePMECUDA ) {
          if(CkMyRank()==0) {
            CProxy_ComputePmeCUDAMgr pme(CkpvAccess(BOCclass_group).computePmeCUDAMgr);
            pme[CkMyNode()].initialize_pencils(new CkQdMsg);
          }
        } else
#endif
        {
          CProxy_ComputePmeMgr pme(CkpvAccess(BOCclass_group).computePmeMgr);
          pme[CkMyPe()].initialize_pencils(new CkQdMsg);          
        }
      }
    }
#ifdef CHARM_HAS_MSA
    else if ( simParameters->MSMOn && ! simParameters->MsmSerialOn ) {
      CProxy_ComputeMsmMsaMgr msm(CkpvAccess(BOCclass_group).computeMsmMsaMgr);
      msm[CkMyPe()].initWorkers(new CkQdMsg);
    }
#else
    else if ( simParameters->MSMOn && ! simParameters->MsmSerialOn ) {
      CProxy_ComputeMsmMgr msm(CkpvAccess(BOCclass_group).computeMsmMgr);
      msm[CkMyPe()].update(new CkQdMsg);
    }
#endif

    #ifdef MEM_OPT_VERSION
    //Now every processor has all the atoms it needs to create the HomePatches.
    //The HomePatches are created in parallel on every home patch procs.
    ioMgr->createHomePatches();
    #else
    if (!CkMyPe()) {
      workDistrib->distributeHomePatches();          
    }
    #endif
  break;

  case 9:
    if ( simParameters->PMEOn ) {
      if ( simParameters->useOptPME ) {
	CProxy_OptPmeMgr pme(CkpvAccess(BOCclass_group).computePmeMgr);
	pme[CkMyPe()].activate_pencils(new CkQdMsg);
      }
      else {
        #ifdef OPENATOM_VERSION
        if ( simParameters->openatomOn ) { 
          CProxy_ComputeMoaMgr moa(CkpvAccess(BOCclass_group).computeMoaMgr); 
          moa[CkMyPe()].startWorkers(new CkQdMsg);
        }
        #endif // OPENATOM_VERSION
#ifdef NAMD_CUDA
        if ( simParameters->usePMECUDA ) {
          if(CkMyRank()==0) {
            CProxy_ComputePmeCUDAMgr pme(CkpvAccess(BOCclass_group).computePmeCUDAMgr);
            pme[CkMyNode()].activate_pencils(new CkQdMsg);
          }
        } else
#endif
        {
          CProxy_ComputePmeMgr pme(CkpvAccess(BOCclass_group).computePmeMgr);
          pme[CkMyPe()].activate_pencils(new CkQdMsg);          
        }
      }
    }
#ifdef CHARM_HAS_MSA
    else if ( simParameters->MSMOn && ! simParameters->MsmSerialOn ) {
      CProxy_ComputeMsmMsaMgr msm(CkpvAccess(BOCclass_group).computeMsmMsaMgr);
      msm[CkMyPe()].startWorkers(new CkQdMsg);
    }
#else
    /*
    else if ( simParameters->MSMOn && ! simParameters->MsmSerialOn ) {
      CProxy_ComputeMsmMgr msm(CkpvAccess(BOCclass_group).computeMsmMgr);
      //msm[CkMyPe()].startWorkers(new CkQdMsg);
    }
    */
#endif

    proxyMgr->createProxies();  // need Home patches before this
    if (!CkMyPe()) LdbCoordinator::Object()->createLoadBalancer();

#ifdef NAMD_TCL
    // TclInitSubsystems() has a race condition so we create one interp per node here
    if (CkMyPe() && CkMyNodeSize() > 1 && ! CkMyRank()) Tcl_DeleteInterp(Tcl_CreateInterp());
#endif

#ifdef USE_NODEPATCHMGR
	//at this point, PatchMap info has been recved on PEs. It is time to create
	//the home patch spanning tree for receiving proxy list info
	if(proxyMgr->getSendSpanning() || proxyMgr->getRecvSpanning()) {
		if(CkMyRank()==0) {
			CProxy_NodeProxyMgr npm(CkpvAccess(BOCclass_group).nodeProxyMgr);
			npm[CkMyNode()].ckLocalBranch()->createSTForHomePatches(PatchMap::Object());
		}
	}
#endif

  break;

  case 10:

    // DMK - DEBUG - If, in MIC runs, the debug option to dump all the compute maps to files
    //   for debugging/verification purposes has been enabled, have each PE do so now.
    #if defined(NAMD_MIC)
      mic_dumpHostDeviceComputeMap();
    #endif

    if (!CkMyPe()) {
      iout << iINFO << "CREATING " << ComputeMap::Object()->numComputes()
           << " COMPUTE OBJECTS\n" << endi;
    }
    DebugM(4,"Creating Computes\n");
    computeMgr->createComputes(ComputeMap::Object());
    DebugM(4,"Building Sequencers\n");
    buildSequencers();
    DebugM(4,"Initializing LDB\n");
    LdbCoordinator::Object()->initialize(PatchMap::Object(),ComputeMap::Object());
  break;

  case 11:
    // computes may create proxies on the fly so put these in separate phase
    Sync::Object()->openSync();  // decide if to open local Sync 
    if (proxySendSpanning || proxyRecvSpanning ) proxyMgr->buildProxySpanningTree();
#ifdef CMK_BALANCED_INJECTION_API
    if(CkMyRank() == 0){
      // CkPrintf("[%d] get retrieved BI=%d\n",CkMyPe(),balancedInjectionLevel);
      ck_set_GNI_BIConfig(balancedInjectionLevel);
      // CkPrintf("[%d] set retrieved BI=%d\n",CkMyPe(),ck_get_GNI_BIConfig());
    }
#endif

  break;

  case 12:
    {
	//For debugging
	/*if(!CkMyPe()){
	FILE *dumpFile = fopen("/tmp/NAMD_Bench.dump", "w");
	dumpbench(dumpFile);
	NAMD_die("Normal execution\n");
	}*/
    }
    #ifdef MEM_OPT_VERSION
    //free space in the Molecule object that are not used anymore
    ioMgr->freeMolSpace();
    #endif
    gotoRun = true;
  break;

  default:
    NAMD_bug("Startup Phase has a bug - check case statement");
  break;

  }

  startupPhase++;
  if (!CkMyPe()) {
    if (!gotoRun) {
      CkStartQD(CkIndex_Node::startUp((CkQdMsg*)0),&thishandle);
    } else {
      Node::messageRun();
    }
  }
}

#ifdef OPENATOM_VERSION
void Node::doneMoaStart()
{
#ifdef OPENATOM_VERSION_DEBUG
  CkPrintf("doneMoaStart executed on processor %d.\n", CkMyPe() );
#endif //OPENATOM_VERSION_DEBUG
}
#endif //OPENATOM_VERSION

void Node::namdOneCommInit()
{
  if (CkpvAccess(comm) == NULL) {
    CkpvAccess(comm) = new Communicate();
#ifdef DPMTA
    pvmc_init();
#endif
  }
}

// Namd 1.X style Send/Recv of simulation information

void Node::namdOneRecv() {
  if ( CmiMyRank() ) return;

  MIStream *conv_msg;

  // Receive molecule and simulation parameter information
  simParameters = node_simParameters = new SimParameters;
  //****** BEGIN CHARMM/XPLOR type changes
  parameters = node_parameters = new Parameters();
  //****** END CHARMM/XPLOR type changes
  molecule = node_molecule = new Molecule(simParameters,parameters);

  DebugM(4, "Getting SimParameters\n");
  conv_msg = CkpvAccess(comm)->newInputStream(0, SIMPARAMSTAG);
  simParameters->receive_SimParameters(conv_msg);

  DebugM(4, "Getting Parameters\n");
  conv_msg = CkpvAccess(comm)->newInputStream(0, STATICPARAMSTAG);
  parameters->receive_Parameters(conv_msg);

  DebugM(4, "Getting Molecule\n");
  conv_msg = CkpvAccess(comm)->newInputStream(0, MOLECULETAG);
  // Modified by JLai -- 10.21.11
  molecule->receive_Molecule(conv_msg);
  if(simParameters->goForcesOn) {
    iout << iINFO << "Compute Nodes receiving GoMolecule Information" << "\n" << endi;
    conv_msg = CkpvAccess(comm)->newInputStream(0, MOLECULETAG);
    molecule->receive_GoMolecule(conv_msg);
  } 
  // End of modification
  DebugM(4, "Done Receiving\n");
}

void Node::namdOneSend() {
  node_simParameters = simParameters;
  node_parameters = parameters;
  node_molecule = molecule;

  MOStream *conv_msg;
  // I'm Pe(0) so I send what I know
  DebugM(4, "Sending SimParameters\n");  
  conv_msg = CkpvAccess(comm)->newOutputStream(ALLBUTME, SIMPARAMSTAG, BUFSIZE);
  simParameters->send_SimParameters(conv_msg);

  DebugM(4, "Sending Parameters\n");
  conv_msg = CkpvAccess(comm)->newOutputStream(ALLBUTME, STATICPARAMSTAG, BUFSIZE);
  parameters->send_Parameters(conv_msg);

  DebugM(4, "Sending Molecule\n");
  int bufSize = BUFSIZE;
  if(molecule->numAtoms>=1000000) bufSize = 16*BUFSIZE;
  conv_msg = CkpvAccess(comm)->newOutputStream(ALLBUTME, MOLECULETAG, bufSize);
  // Modified by JLai -- 10.21.11
  molecule->send_Molecule(conv_msg);
  
  if(simParameters->goForcesOn) {
    iout << iINFO <<  "Master Node sending GoMolecule Information" << "\n" << endi;
    conv_msg = CkpvAccess(comm)->newOutputStream(ALLBUTME, MOLECULETAG, bufSize);
    molecule->send_GoMolecule(conv_msg);
  } // End of modification
}


void Node::reloadStructure(const char *fname, const char *pdbname) {
  delete molecule;
  molecule = state->molecule = 0;
  delete pdb;
  pdb = state->pdb = 0;
  state->loadStructure(fname,pdbname,1);
  this->molecule = state->molecule;
  this->pdb = state->pdb;
  CProxy_Node nodeProxy(thisgroup);
  nodeProxy.resendMolecule();
}


void Node::resendMolecule() {
  if ( CmiMyRank() ) {
    return;
  }
  if ( CmiMyPe() == 0 ) {
    int bufSize = BUFSIZE;
    MOStream *conv_msg;
    conv_msg = CkpvAccess(comm)->newOutputStream(ALLBUTME, STATICPARAMSTAG, bufSize);
    parameters->send_Parameters(conv_msg);
    if(molecule->numAtoms>=1000000) bufSize = 16*BUFSIZE;
    conv_msg = CkpvAccess(comm)->newOutputStream(ALLBUTME, MOLECULETAG, bufSize);
    molecule->send_Molecule(conv_msg);
  } else {
    MIStream *conv_msg;
    delete parameters;
    parameters = new Parameters;
    conv_msg = CkpvAccess(comm)->newInputStream(0, STATICPARAMSTAG);
    parameters->receive_Parameters(conv_msg);
    delete molecule;
    molecule = new Molecule(simParameters,parameters);
    conv_msg = CkpvAccess(comm)->newInputStream(0, MOLECULETAG);
    molecule->receive_Molecule(conv_msg);
  }
  node_parameters = parameters;
  node_molecule = molecule;
  SimParameters::nonbonded_select();
  computeMgr->sendBuildCudaExclusions();
  CProxy_Node nodeProxy(thisgroup);
  for ( int i=0; i<CmiMyNodeSize(); ++i ) {
    nodeProxy[CmiMyPe()+i].resendMolecule2();
  }
}

void Node::resendMolecule2() {
  parameters = node_parameters;
  molecule = node_molecule;
  AtomMap::Object()->allocateMap(molecule->numAtoms);
}


// Initial thread setup

void Node::threadInit() {
  // Thread initialization
  if (CthImplemented()) {
    CthSetStrategyDefault(CthSelf());
  } else {
    NAMD_bug("Node::startup() Oh no, tiny elvis, threads not implemented");
  }
}

//
void Node::buildSequencers() {
  HomePatchList *hpl = PatchMap::Object()->homePatchList();
  ResizeArrayIter<HomePatchElem> ai(*hpl);

  // Controller object is only on Pe(0)
  if ( ! CkMyPe() ) {
    Controller *controller = new Controller(state);
    state->useController(controller);
  }

  // Assign Sequencer to all HomePatch(es)
  for (ai=ai.begin(); ai != ai.end(); ai++) {
    HomePatch *patch = (*ai).patch;
    Sequencer *sequencer = new Sequencer(patch);
    patch->useSequencer(sequencer);
  }
}



//-----------------------------------------------------------------------
// Node run() - broadcast to all nodes
//-----------------------------------------------------------------------
void Node::messageRun() {
  (CProxy_Node(CkpvAccess(BOCclass_group).node)).run();
}


//-----------------------------------------------------------------------
// run(void) runs the specified simulation for the specified number of
// steps, overriding the contents of the configuration file
//-----------------------------------------------------------------------
void Node::run()
{
  // Start Controller (aka scalar Sequencer) on Pe(0)
//  printf("\n\n I am in Node.C in run method about to call  state->runController\n\n");
  if ( ! CkMyPe() ) {
    state->runController();
  }

  DebugM(4, "Starting Sequencers\n");
  // Run Sequencer on each HomePatch - i.e. start simulation
  HomePatchList *hpl = PatchMap::Object()->homePatchList();
  ResizeArrayIter<HomePatchElem> ai(*hpl);
  for (ai=ai.begin(); ai != ai.end(); ai++) {
    HomePatch *patch = (*ai).patch;
//CkPrintf("Proc#%d in Node calling Sequencer ",CkMyPe());
    patch->runSequencer();
  }

  if (!CkMyPe()) {
    double newTime = CmiWallTimer();
    iout << iINFO << "Startup phase " << startupPhase-1 << " took "
	 << newTime - startupTime << " s, "
	 << memusage_MB() << " MB of memory in use\n";
    iout << iINFO << "Finished startup at " << newTime << " s, "
	 << memusage_MB() << " MB of memory in use\n\n" << endi;
    fflush(stdout);
  }
  
}


//-----------------------------------------------------------------------
// Node scriptBarrier() - twiddle parameters with simulation halted
//-----------------------------------------------------------------------

void Node::enableScriptBarrier() {
  CkStartQD(CkIndex_Node::scriptBarrier((CkQdMsg*)0),&thishandle);
}

void Node::scriptBarrier(CkQdMsg *qmsg) {
  delete qmsg;
  //script->awaken();
}

void Node::scriptParam(ScriptParamMsg *msg) {
  simParameters->scriptSet(msg->param,msg->value);
  delete msg;
}

void Node::reloadCharges(const char *filename) {
  FILE *file = fopen(filename,"r");
  if ( ! file ) NAMD_die("node::reloadCharges():Error opening charge file.");

  int n = molecule->numAtoms;
  float *charge = new float[n];

  for ( int i = 0; i < n; ++i ) {
    if ( ! fscanf(file,"%f",&charge[i]) )
      NAMD_die("Node::reloadCharges():Not enough numbers in charge file.");
  }

  fclose(file);
  CProxy_Node(thisgroup).reloadCharges(charge,n);
  delete [] charge;
}

void Node::reloadCharges(float charge[], int n) {
  molecule->reloadCharges(charge,n);
}


// BEGIN gf
void Node::reloadGridforceGrid(const char * key) {
    DebugM(4, "reloadGridforceGrid(const char*) called on node " << CkMyPe() << "\n" << endi);
    
    int gridnum;
    MGridforceParams *mgridParams;
    if (key == NULL) {
	gridnum = simParameters->mgridforcelist.index_for_key(MGRIDFORCEPARAMS_DEFAULTKEY);
	mgridParams = simParameters->mgridforcelist.find_key(MGRIDFORCEPARAMS_DEFAULTKEY);
    } else {
	gridnum = simParameters->mgridforcelist.index_for_key(key);
	mgridParams = simParameters->mgridforcelist.find_key(key);
    }
    
    if (gridnum < 0 || mgridParams == NULL) {
	NAMD_die("Node::reloadGridforceGrid(const char*):Could not find grid.");
    }
    
    GridforceGrid *grid = molecule->get_gridfrc_grid(gridnum);
    if (grid == NULL) {
	NAMD_bug("Node::reloadGridforceGrid(const char*):grid not found");
    }
    grid->reinitialize(simParameters, mgridParams);
    
    CProxy_Node(thisgroup).reloadGridforceGrid(gridnum);
    
    DebugM(4, "reloadGridforceGrid(const char*) finished\n" << endi);
}

void Node::updateGridScale(char* key, Vector scale) {
    DebugM(4, "updateGridScale(char*, Vector) called on node " << CkMyPe() << "\n" << endi);
    
    int gridnum;
    MGridforceParams* mgridParams;
    if (key == NULL) {
	gridnum = simParameters->mgridforcelist.index_for_key(MGRIDFORCEPARAMS_DEFAULTKEY);
	mgridParams = simParameters->mgridforcelist.find_key(MGRIDFORCEPARAMS_DEFAULTKEY);
    } else {
	gridnum = simParameters->mgridforcelist.index_for_key(key);
	mgridParams = simParameters->mgridforcelist.find_key(key);
    }

    if (gridnum < 0 || mgridParams == NULL) {
	NAMD_die("Node::updateGridScale(char*, Vector): Could not find grid.");
    }
    
    GridforceGrid* grid = molecule->get_gridfrc_grid(gridnum);
    if (grid == NULL) {
	NAMD_bug("Node::updateGridScale(char*, Vector): grid not found");
    }
    CProxy_Node(thisgroup).updateGridScale(gridnum, scale.x, scale.y, scale.z);
    
    DebugM(4, "updateGridScale(char*, Vector) finished\n" << endi);
}
void Node::updateGridScale(int gridnum, float sx, float sy, float sz) {
    if (CmiMyRank()) return;
    DebugM(4, "updateGridScale(char*, int, float, float, float) called on node " << CkMyPe() << "\n" << endi);
       
    GridforceGrid *grid = molecule->get_gridfrc_grid(gridnum);
    if (grid == NULL) {
	NAMD_bug("Node::updateGridScale(char*, int, float, float, float):grid not found");
    }
    
    Vector scale(sx,sy,sz);
    simParameters->mgridforcelist.at_index(gridnum)->gridforceScale = scale;
    grid->set_scale( scale );

    DebugM(4, "updateGridScale(char*, int, float, float, float) finished\n" << endi);
}

void Node::reloadGridforceGrid(int gridnum) {
    if (CmiMyRank()) return;
    DebugM(4, "reloadGridforceGrid(int) called on node " << CkMyPe() << "\n" << endi);
    
    GridforceGrid *grid = molecule->get_gridfrc_grid(gridnum);
    if (grid == NULL) {
	NAMD_bug("Node::reloadGridforceGrid(int):grid not found");
    }
    
    if (CkMyPe()) {
	// not node 0 -> receive grid
	DebugM(4, "Receiving grid\n");
	
	delete grid;
	
	MIStream *msg = CkpvAccess(comm)->newInputStream(0, GRIDFORCEGRIDTAG);
	grid = GridforceGrid::unpack_grid(gridnum, msg);
	molecule->set_gridfrc_grid(gridnum, grid);
	delete msg;
    } else {
	// node 0 -> send grid
	DebugM(4, "Sending grid\n");
	
	MOStream *msg = CkpvAccess(comm)->newOutputStream(ALLBUTME, GRIDFORCEGRIDTAG, BUFSIZE);
	GridforceGrid::pack_grid(grid, msg);
	msg->end();
	delete msg;
    }
    
    DebugM(4, "reloadGridforceGrid(int) finished\n" << endi);
}
// END gf


// initiating replica
void Node::sendCheckpointReq(int remote, const char *key, int task, Lattice &lat, ControllerState &cs) {
  CheckpointMsg *msg = new (1+strlen(key),0) CheckpointMsg;
  msg->replica = CmiMyPartition();
  msg->task = task;
  msg->checkpoint.lattice = lat;
  msg->checkpoint.state = cs;
  strcpy(msg->key,key);
  envelope *env = UsrToEnv(CheckpointMsg::pack(msg));
  CmiSetHandler(env,recvCheckpointCReq_index);
#if CMK_HAS_PARTITION
  CmiInterSyncSendAndFree(CkMyPe(),remote,env->getTotalsize(),(char*)env);
#else
  CmiSyncSendAndFree(CkMyPe(),env->getTotalsize(),(char*)env);
#endif
}

// responding replica
extern "C" {
  void recvCheckpointCReq_handler(envelope *env) {
    Node::Object()->recvCheckpointReq(CheckpointMsg::unpack(EnvToUsr(env)));
  }
}

// responding replica
void Node::recvCheckpointReq(CheckpointMsg *msg) {
  state->controller->recvCheckpointReq(msg->key,msg->task,msg->checkpoint);

  int remote = msg->replica;
  msg->replica = CmiMyPartition();
  envelope *env = UsrToEnv(CheckpointMsg::pack(msg));
  CmiSetHandler(env,recvCheckpointCAck_index);
#if CMK_HAS_PARTITION
  CmiInterSyncSendAndFree(CkMyPe(),remote,env->getTotalsize(),(char*)env);
#else
  CmiSyncSendAndFree(CkMyPe(),env->getTotalsize(),(char*)env);
#endif
}

// initiating replica
extern "C" {
  void recvCheckpointCAck_handler(envelope *env) {
    Node::Object()->recvCheckpointAck(CheckpointMsg::unpack(EnvToUsr(env)));
  }
}

// initiating replica
void Node::recvCheckpointAck(CheckpointMsg *msg) {
  state->controller->recvCheckpointAck(msg->checkpoint);
  delete msg;
}


void Node::sendEnableExitScheduler(void) {
  //CmiPrintf("sendEnableExitScheduler\n");
  CkQdMsg *msg = new CkQdMsg;
  CProxy_Node nodeProxy(thisgroup);
  nodeProxy[0].recvEnableExitScheduler(msg);
}

void Node::recvEnableExitScheduler(CkQdMsg *msg) {
  //CmiPrintf("recvEnableExitScheduler\n");
  delete msg;
  enableExitScheduler();
}

void Node::enableExitScheduler(void) {
  if ( CkMyPe() ) {
    sendEnableExitScheduler();
  } else {
    CkStartQD(CkIndex_Node::exitScheduler((CkQdMsg*)0),&thishandle);
  }
}

void Node::exitScheduler(CkQdMsg *msg) {
  //CmiPrintf("exitScheduler %d\n",CkMyPe());
  CsdExitScheduler();
  delete msg;
}

void Node::sendEnableEarlyExit(void) {
  CkQdMsg *msg = new CkQdMsg;
  CProxy_Node nodeProxy(thisgroup);
  nodeProxy[0].recvEnableEarlyExit(msg);
}

void Node::recvEnableEarlyExit(CkQdMsg *msg) {
  delete msg;
  enableEarlyExit();
}

void Node::enableEarlyExit(void) {
  if ( CkMyPe() ) {
    sendEnableEarlyExit();
  } else {
    CkStartQD(CkIndex_Node::earlyExit((CkQdMsg*)0),&thishandle);
  }
}

void Node::earlyExit(CkQdMsg *msg) {
  iout << iERROR << "Exiting prematurely; see error messages above.\n" << endi;
  if ( CmiNumPartitions() > 1 ) NAMD_quit("Exiting prematurely; see error messages above.");
  BackEnd::exit();
  delete msg;
}


//------------------------------------------------------------------------
// Some odd utilities
//------------------------------------------------------------------------
void Node::saveMolDataPointers(NamdState *state)
{
  this->molecule = state->molecule;
  this->parameters = state->parameters;
  this->simParameters = state->simParameters;
  this->configList = state->configList;
  this->pdb = state->pdb;
  this->state = state;
}

// entry methods for BG/P HPM (performance counters) library
void Node::startHPM() {
#if USE_HPM
  HPM_Start("500 steps", localRankOnNode);
#endif
}

void Node::stopHPM() {
#if USE_HPM
  HPM_Stop("500 steps", localRankOnNode);
  HPM_Print(CkMyPe(), localRankOnNode);
#endif
}

void Node::traceBarrier(int turnOnTrace, int step){
	curTimeStep = step;
	if(turnOnTrace) traceBegin();
	else traceEnd();

    if(turnOnTrace) CmiTurnOnStats();
    else CmiTurnOffStats();

	//CkPrintf("traceBarrier (%d) at step %d called on proc %d\n", turnOnTrace, step, CkMyPe());	
	CProxy_Node nd(CkpvAccess(BOCclass_group).node);
	CkCallback cb(CkIndex_Node::resumeAfterTraceBarrier(NULL), nd[0]);
	contribute(0, NULL, CkReduction::sum_int, cb);
	
}

void Node::resumeAfterTraceBarrier(CkReductionMsg *msg){
	CmiAssert(CmiMyPe()==0);
	delete msg;	
	state->controller->resumeAfterTraceBarrier(curTimeStep);
}

void Node::papiMeasureBarrier(int turnOnMeasure, int step){
#ifdef MEASURE_NAMD_WITH_PAPI
	curMFlopStep = step;
	double results[NUM_PAPI_EVENTS+1];

	if(turnOnMeasure){		
	  CkpvAccess(papiEvents)[NUM_PAPI_EVENTS]=CmiWallTimer();

	  long long counters[NUM_PAPI_EVENTS+1];
	  int ret=PAPI_start_counters(CkpvAccess(papiEvents), NUM_PAPI_EVENTS);
	  if(ret==PAPI_OK)
	    {
	      //	      CkPrintf("traceBarrier start counters (%d) at step %d called on proc %d\n", turnOnMeasure, step, CkMyPe());
	    }
	  else
	    {
	      CkPrintf("error PAPI_start_counters (%d) at step %d called on proc %d\n",ret , step, CkMyPe());
	    }
	  if(PAPI_read_counters(counters, NUM_PAPI_EVENTS)!=PAPI_OK)
	    {
	      CkPrintf("error PAPI_read_counters %d\n",PAPI_read_counters(counters, NUM_PAPI_EVENTS));
	    };
	}else{
	  long long counters[NUM_PAPI_EVENTS+1];
	  for(int i=0;i<NUM_PAPI_EVENTS;i++)  counters[i]=0LL;
	  if(PAPI_read_counters(counters, NUM_PAPI_EVENTS)==PAPI_OK)
	    {
#if !MEASURE_PAPI_SPP
	      results[0] = (double)counters[0]/1e6;
	      results[1] = (double)counters[1]/1e6;
#else
	      for(int i=0;i<NUM_PAPI_EVENTS;i++)  results[i] = counters[i]/1e6;
#endif
	      //	      for(int i=0;i<NUM_PAPI_EVENTS;i++) CkPrintf("[%d] counter %d is %ld\n",CkMyPe(),i,counters[i]);
	    }
	  else
	    {
	      //	      CkPrintf("error PAPI_read_counters %d\n",PAPI_read_counters(counters, NUM_PAPI_EVENTS));
	    }
	  //	  CkPrintf("traceBarrier stop counters (%d) at step %d called on proc %d\n", turnOnMeasure, step, CkMyPe());
		
	  PAPI_stop_counters(counters, NUM_PAPI_EVENTS);	
	}
	if(CkMyPe()==0)
	  //	    CkPrintf("traceBarrier (%d) at step %d called on proc %d\n", turnOnMeasure, step, CkMyPe());
	results[NUM_PAPI_EVENTS]=CkpvAccess(papiEvents)[NUM_PAPI_EVENTS]; //starttime
	CProxy_Node nd(CkpvAccess(BOCclass_group).node);
	CkCallback cb(CkIndex_Node::resumeAfterPapiMeasureBarrier(NULL), nd[0]);
	contribute(sizeof(double)*(NUM_PAPI_EVENTS+1), &results, CkReduction::sum_double, cb);	
#endif
}

void Node::resumeAfterPapiMeasureBarrier(CkReductionMsg *msg){
#ifdef MEASURE_NAMD_WITH_PAPI
  
	if(simParameters->papiMeasureStartStep != curMFlopStep) {
		double *results = (double *)msg->getData();
		double endtime=CmiWallTimer();
		int bstep = simParameters->papiMeasureStartStep;
		int estep = bstep + simParameters->numPapiMeasureSteps;
#if MEASURE_PAPI_SPP
		CkPrintf("SPP INFO: PAPI_FP_OPS timestep %d to %d is %lf(1e6)\n", bstep,estep,results[0]);
		CkPrintf("SPP INFO: PAPI_TOT_INS timestep %d to %d is %lf(1e6)\n", bstep,estep,results[1]);
		CkPrintf("SPP INFO: perf::PERF_COUNT_HW_CACHE_LL:MISS timestep %d to %d is %lf(1e6)\n", bstep,estep,results[2]);
		CkPrintf("SPP INFO: DATA_PREFETCHER:ALL timestep %d to %d is %lf(1e6)\n", bstep,estep,results[3]);
		CkPrintf("SPP INFO: PAPI_L1_DCA timestep %d to %d is %lf(1e6)\n", bstep,estep,results[4]);
		CkPrintf("SPP INFO: PAPI_TOT_CYC timestep %d to % is %lf(1e6)\n", bstep,estep,results[5]);
		//		CkPrintf("SPP INFO: INSTRUCTION_FETCH_STALL timestep %d to %d is %lf(1e6)\n", bstep,estep,results[6]);
		//		CkPrintf("SPP INFO: WALLtime timestep %d to %d is %lf\n", bstep,estep,endtime-results[NUM_PAPI_EVENTS]/CkNumPes());
		CkPrintf("SPP INFO: WALLtime timestep %d to %d is %lf\n", bstep,estep,endtime-results[NUM_PAPI_EVENTS]);
		CkPrintf("SPP INFO: endtime %lf avgtime %lf tottime %lf\n", endtime,results[NUM_PAPI_EVENTS]/CkNumPes(),results[NUM_PAPI_EVENTS] );
#else
		if(CkpvAccess(papiEvents)[0] == PAPI_FP_INS){
			double totalFPIns = results[0];
			if(CkpvAccess(papiEvents)[1] == PAPI_FMA_INS) totalFPIns += (results[1]*2);
			CkPrintf("FLOPS INFO: from timestep %d to %d, the total FP instruction of NAMD is %lf(x1e6) per processor\n", 
					 bstep, estep, totalFPIns/CkNumPes());
		}else{
			char nameBuf[PAPI_MAX_STR_LEN];
			CkPrintf("PAPI COUNTERS INFO: from timestep %d to %d, ", 
					 bstep, estep);
			for(int i=0; i<NUM_PAPI_EVENTS; i++) {
				PAPI_event_code_to_name(CkpvAccess(papiEvents)[i], nameBuf);
				CkPrintf("%s is %lf(x1e6), ", nameBuf, results[i]/CkNumPes());
			}
			CkPrintf("per processor\n");
		}		
#endif
	}
	delete msg;	
	state->controller->resumeAfterPapiMeasureBarrier(curMFlopStep);
#endif
}

extern char *gNAMDBinaryName;
void Node::outputPatchComputeMaps(const char *filename, int tag){
	if(!simParameters->outputMaps && !simParameters->simulateInitialMapping) return;

	int numpes = CkNumPes();
	int nodesize = CkMyNodeSize();
	if(simParameters->simulateInitialMapping) {
		numpes = simParameters->simulatedPEs;
		nodesize = simParameters->simulatedNodeSize;
	}

	char fname[128];
	sprintf(fname, "mapdump_%s.%d_%d_%d_%s", filename, numpes, nodesize, tag, gNAMDBinaryName);

	FILE *fp = fopen(fname, "w");
	if(fp == NULL) {
		NAMD_die("Error in outputing PatchMap and ComputeMap info!\n");
		return;
	}
	PatchMap *pMap = PatchMap::Object();
	ComputeMap *cMap = ComputeMap::Object();
	int numPatches = pMap->numPatches();
	int numComputes = cMap->numComputes();
	fprintf(fp, "%d %d %d %d %d %d %d\n", numpes, nodesize, numPatches, numComputes, 
			pMap->gridsize_a(), pMap->gridsize_b(), pMap->gridsize_c());
	//output PatchMap info
	for(int i=0; i<numPatches; i++) {
	#ifdef MEM_OPT_VERSION
		fprintf(fp, "%d %d\n", pMap->numAtoms(i), pMap->node(i));
	#else
		fprintf(fp, "%d %d\n", pMap->patch(i)->getNumAtoms(), pMap->node(i));
	#endif
	}

	//output ComputeMap info
	for(int i=0; i<numComputes; i++) {		
		fprintf(fp, "%d %d %d %d\n", cMap->node(i), cMap->type(i), cMap->pid(i,0), cMap->pid(i,1));		
	}
}


//======================================================================
// Private functions

#include "Node.def.h"

