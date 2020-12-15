/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "BackEnd.h"
#include "ProcessorPrivate.h"
#include "common.h"
#include "Node.h"
#include "memusage.h"

#include <new>
#if defined(WIN32) && !defined(__CYGWIN__)
#include <new.h>
#endif

#include "Lattice.h"
#include "ComputeMoa.h" 
#include "ComputeMsmMsa.h"  // needed for MsmMsaData definition
#include "main.decl.h"
#include "main.h"
#include "BOCgroup.h"
#include "WorkDistrib.decl.h"
#include "ProxyMgr.decl.h"
#include "PatchMgr.decl.h"
#include "DataExchanger.decl.h"
#ifdef CHARM_HAS_MSA
#include "ComputeMgr.decl.h"
#endif
#include "ReductionMgr.decl.h"
#include "CollectionMgr.decl.h"
#include "CollectionMaster.decl.h"
#include "CollectionMgr.h"
#include "CollectionMaster.h"
#include "BroadcastMgr.decl.h"
#include "LdbCoordinator.decl.h"
#include "Sync.decl.h"

#ifdef MEM_OPT_VERSION
#include "ParallelIOMgr.decl.h"
#endif

#ifdef NAMD_TCL
#include <tcl.h>
#endif

extern void _initCharm(int, char**);

float cpuTime_start;
float wallTime_start;

CkpvStaticDeclare(int,exitSchedHndlr);

extern "C" void exit_sched(void* msg)
{
  //  CmiPrintf("Exiting scheduler on %d\n",CmiMyPe());
  CsdExitScheduler();
}

static void register_exit_sched(void)
{
  CkpvInitialize(int,exitSchedHndlr);
  CkpvAccess(exitSchedHndlr) = CmiRegisterHandler((CmiHandler)exit_sched);
}

void BackEnd::ExitSchedOn(int pe)
{
  void* msg = CmiAlloc(CmiMsgHeaderSizeBytes);
  CmiSetHandler(msg,CkpvAccess(exitSchedHndlr));
  CmiSyncSendAndFree(pe,CmiMsgHeaderSizeBytes,(char *)msg);
}

#if defined(WIN32) && !defined(__CYGWIN__) && !defined(__MINGW_H)
int NAMD_new_handler(size_t) {
#else
void NAMD_new_handler() {
#endif
  char tmp[100];
  sprintf(tmp,"Memory allocation failed on processor %d.",CmiMyPe());
  NAMD_die(tmp);
#if defined(WIN32) && !defined(__CYGWIN__) && !defined(__MINGW_H)
  return 0;
#endif
}

void topo_getargs(char**);
void cuda_getargs(char**);
// void cuda_initialize();
void mic_getargs(char**);
// void mic_initialize();

// called on all procs
void all_init(int argc, char **argv)
{
#if defined(WIN32) && !defined(__CYGWIN__) && !defined(__MINGW_H)
  _set_new_handler(NAMD_new_handler);
#else
  std::set_new_handler(NAMD_new_handler);
#endif
  ProcessorPrivateInit();
  register_exit_sched();
  CmiGetArgFlag(argv, "+idlepoll");  // remove +idlepoll if it's still there
  topo_getargs(argv);
#ifdef NAMD_CUDA
  cuda_getargs(argv);
  argc = CmiGetArgc(argv);
#endif
#ifdef NAMD_MIC
  CmiGetArgFlag(argv, "+idlepoll");  // remove +idlepoll if it's still there
  mic_getargs(argv);
  argc = CmiGetArgc(argv);
#endif
  
  _initCharm(argc, argv);  // message main Chare

//#if 0  // moved to WorkDistrib
//#ifdef NAMD_CUDA
//  if ( CkMyPe() < CkNumPes() ) cuda_initialize();
//#endif
//#ifdef NAMD_MIC
//  if ( CkMyPe() < CkNumPes() ) mic_initialize();
//#endif
//#endif
}

extern void after_backend_init(int argc, char **argv);
void master_init(int argc, char **argv);

// called on slave procs
void slave_init(int argc, char **argv)
{
#if CMK_SMP
  //the original main thread could now be a comm thread
  //and a slave thread could now be the main thread,
  //so we have to do the master initialization here
  if(CmiMyRank()==0){
    master_init(argc, argv);
    if(CmiMyPe()==0)
      after_backend_init(argc, argv);
    return;
  }
#endif

  all_init(argc, argv);

  if (CkMyRank() < CkMyNodeSize()) 	// skip the communication thread
    CsdScheduler(-1);
}

void master_init(int argc, char **argv){
  cpuTime_start = CmiCpuTimer();
  wallTime_start = CmiWallTimer();
  if ( CmiMyPe() ) {
    all_init(argc, argv);
    CsdScheduler(-1);
    ConverseExit();  // should never return
  }

  all_init(argc, argv);

  // Create branch-office chares
  BOCgroup group;
  group.workDistrib = CProxy_WorkDistrib::ckNew();
  group.proxyMgr = CProxy_ProxyMgr::ckNew();
  group.patchMgr = CProxy_PatchMgr::ckNew();
  group.computeMgr = CProxy_ComputeMgr::ckNew();
  group.reductionMgr = CProxy_ReductionMgr::ckNew();
  // group.computePmeMgr set in constructor during startup
  group.nodePmeMgr = CProxy_NodePmeMgr::ckNew();
#ifdef NAMD_CUDA
  group.computePmeCUDAMgr = CProxy_ComputePmeCUDAMgr::ckNew();
  group.computeCUDAMgr = CProxy_ComputeCUDAMgr::ckNew();
#endif
#ifdef OPENATOM_VERSION
  group.computeMoaMgr = CProxy_ComputeMoaMgr::ckNew();
#endif // OPENATOM_VERSION
  group.computeExtMgr = CProxy_ComputeExtMgr::ckNew();
  group.computeQMMgr = CProxy_ComputeQMMgr::ckNew();
  group.computeGBISserMgr = CProxy_ComputeGBISserMgr::ckNew();
  group.computeFmmSerialMgr = CProxy_ComputeFmmSerialMgr::ckNew();
  group.computeMsmSerialMgr = CProxy_ComputeMsmSerialMgr::ckNew();
#ifdef CHARM_HAS_MSA
  group.computeMsmMsaMgr = CProxy_ComputeMsmMsaMgr::ckNew();
#endif
  group.computeMsmMgr = CProxy_ComputeMsmMgr::ckNew();
  // Charm CkMulticast library module
  group.multicastMgr = CProxy_CkMulticastMgr::ckNew();
#ifdef MEM_OPT_VERSION
  group.ioMgr=CProxy_ParallelIOMgr::ckNew();
#endif

  group.sync = CProxy_Sync::ckNew();

  #ifdef USE_NODEPATCHMGR
  group.nodeProxyMgr = CProxy_NodeProxyMgr::ckNew();
  #endif
  
#if     CMK_SMP && USE_CKLOOP
  group.ckLoop = CkLoop_Init();
#endif

  CkChareID collectionMaster = CProxy_CollectionMaster::ckNew(0);  
  SlaveInitMsg *initmsg7 = new SlaveInitMsg;
  initmsg7->master = collectionMaster;
  group.collectionMgr = CProxy_CollectionMgr::ckNew(initmsg7);

  group.broadcastMgr = CProxy_BroadcastMgr::ckNew();
  group.ldbCoordinator = CProxy_LdbCoordinator::ckNew();

  group.dataExchanger = CProxy_DataExchanger::ckNew();

  GroupInitMsg *msg = new GroupInitMsg;
  msg->group = group;
  CProxy_Node::ckNew(msg);
 
}

char *gNAMDBinaryName = NULL;
// called by main on one or all procs
void BackEnd::init(int argc, char **argv) {

  gNAMDBinaryName = argv[0]+strlen(argv[0])-1;
  while(gNAMDBinaryName != argv[0]){
    if(*gNAMDBinaryName=='/' || *gNAMDBinaryName=='\\'){
      gNAMDBinaryName++;
      break;
    }
    gNAMDBinaryName--;
  }

#if defined(NAMD_CUDA) || defined(NAMD_MIC)
  // look for but don't remove +idlepoll on command line
  int idlepoll = 0;
  for ( int i = 0; i < argc; ++i ) {
    if ( 0==strcmp(argv[i],"+idlepoll") ) {
      idlepoll = 1;
      break;
    }
  }
#endif

  ConverseInit(argc, argv, slave_init, 1, 1);  // calls slave_init on others

// idlepoll only matters for non-smp UDP layer
#if (defined(NAMD_CUDA) || defined(NAMD_MIC)) && CMK_NET_VERSION && CMK_SHARED_VARS_UNAVAILABLE && CMK_WHEN_PROCESSOR_IDLE_USLEEP && ! CMK_USE_IBVERBS && ! CMK_USE_TCP
  if ( ! idlepoll ) {
    NAMD_die("Please add +idlepoll to command line for proper performance.");
  }
#endif

  master_init(argc, argv);
}

void cuda_finalize();

// called on proc 0 by front end
void BackEnd::exit(void) {
  float cpuTime = CmiCpuTimer() - cpuTime_start;
  float wallTime = CmiWallTimer() - wallTime_start;
  CmiPrintf("====================================================\n\n"
	    "WallClock: %f  CPUTime: %f  Memory: %f MB\n",
	    wallTime, cpuTime, memusage_MB());
#ifdef NAMD_TCL
  Tcl_Finalize();
#endif
#ifdef NAMD_CUDA
  cuda_finalize();
#ifdef __APPLE__
#if 0 && CMK_MULTICORE
  CmiPrintf("EXITING ABNORMALLY TO AVOID HANGING CUDA RUNTIME THREADS\n");
  ::exit(0);
#endif
#endif
#endif
#ifdef NAMD_MIC
#if 0 && CMK_MULTICORE
  CmiPrintf("EXITING ABNORMALLY TO AVOID HANGING MIC OFFLOAD THREADS\n");
#pragma offload target(mic)
  {
    ::exit(0);
  }
#endif
#endif
  CkExit();
}

// start scheduler
void BackEnd::suspend(void) {
  CsdScheduler(-1);
}

// start quiescence detection to return to front end
void BackEnd::awaken(void) {
  Node::Object()->enableExitScheduler();
}

// start QD and scheduler
void BackEnd::barrier(void) {
  awaken();
  suspend();
}

