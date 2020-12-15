/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifdef NAMD_FFTW
//#define MANUAL_DEBUG_FFTW3 1
#ifdef NAMD_FFTW_3
#include <fftw3.h>
#else
// fftw2 doesn't have these defined
#define fftwf_malloc fftw_malloc
#define fftwf_free fftw_free
#ifdef NAMD_FFTW_NO_TYPE_PREFIX
#include <fftw.h>
#include <rfftw.h>
#else
#include <sfftw.h>
#include <srfftw.h>
#endif
#endif
#endif

#include <vector>
#include <algorithm>
#include <deque>
using namespace std;

#include "InfoStream.h"
#include "Node.h"
#include "PatchMap.h"
#include "PatchMap.inl"
#include "AtomMap.h"
#include "ComputePme.h"
#include "ComputePmeMgr.decl.h"
#include "PmeBase.inl"
#include "PmeRealSpace.h"
#include "PmeKSpace.h"
#include "ComputeNonbondedUtil.h"
#include "PatchMgr.h"
#include "Molecule.h"
#include "ReductionMgr.h"
#include "ComputeMgr.h"
#include "ComputeMgr.decl.h"
// #define DEBUGM
#define MIN_DEBUG_LEVEL 3
#include "Debug.h"
#include "SimParameters.h"
#include "WorkDistrib.h"
#include "varsizemsg.h"
#include "Random.h"
#include "ckhashtable.h"
#include "Priorities.h"

#include "ComputeMoa.h"
#include "ComputeMoaMgr.decl.h" 

//#define     USE_RANDOM_TOPO         1

//#define USE_TOPO_SFC                    1
//#define     USE_CKLOOP                1
//#include "TopoManager.h"

#include "DeviceCUDA.h"
#ifdef NAMD_CUDA
#include <cuda_runtime.h>
#include <cuda.h>
void cuda_errcheck(const char *msg);
#ifdef WIN32
#define __thread __declspec(thread)
#endif
extern __thread DeviceCUDA *deviceCUDA;
#endif

#include "ComputePmeCUDAKernel.h"

#ifndef SQRT_PI
#define SQRT_PI 1.7724538509055160273 /* mathematica 15 digits*/
#endif

#if CMK_PERSISTENT_COMM 
#define USE_PERSISTENT      1
#endif

#if USE_PERSISTENT
#define Z_PERSIST 1
#define Y_PERSIST 1
#define X_PERSIST 1
#endif

#if defined(NAMD_CUDA) && defined(MEM_OPT_VERSION)
#define USE_NODE_PAR_RECEIVE    1
#endif

char *pencilPMEProcessors;

class PmeAckMsg : public CMessage_PmeAckMsg {
};

class PmeGridMsg : public CMessage_PmeGridMsg {
public:

  int sourceNode;
  int sequence;
  int hasData;
  Lattice lattice;
  int start;
  int len;
  int zlistlen;
  int *zlist;
  char *fgrid;
  float *qgrid;
  CkArrayIndex3D destElem;
};

class PmeTransMsg : public CMessage_PmeTransMsg {
public:

  int sourceNode;
  int sequence;
  int hasData;
  Lattice lattice;
  int x_start;
  int nx;
  float *qgrid;
  CkArrayIndex3D destElem;
};

class PmeSharedTransMsg : public CMessage_PmeSharedTransMsg {
public:
  PmeTransMsg *msg;
  int *count;
  CmiNodeLock lock;
};

class PmeUntransMsg : public CMessage_PmeUntransMsg {
public:

  int sourceNode;
  int y_start;
  int ny;
  float *qgrid;
  CkArrayIndex3D destElem;
};

class PmeSharedUntransMsg : public CMessage_PmeSharedUntransMsg {
public:
  PmeUntransMsg *msg;
  int *count;
  CmiNodeLock lock;
};

class PmeEvirMsg : public CMessage_PmeEvirMsg {
public:
  PmeReduction *evir;
};

class PmePencilMap : public CBase_PmePencilMap {
public:
  PmePencilMap(int i_a, int i_b, int n_b, int n, int *d)
    : ia(i_a), ib(i_b), nb(n_b),
      size(n), data(newcopyint(n,d)) {
  }
  virtual int registerArray(CkArrayIndexMax&, CkArrayID) {
    //Return an ``arrayHdl'', given some information about the array
    return 0;
  }
  virtual int procNum(int, const CkArrayIndex &i) {
    //Return the home processor number for this element of this array
    return data[ i.data()[ia] * nb + i.data()[ib] ];
  }
  virtual void populateInitial(int, CkArrayIndexMax &, void *msg, CkArrMgr *mgr) {
    int mype = CkMyPe();
    for ( int i=0; i < size; ++i ) {
      if ( data[i] == mype ) {
        CkArrayIndex3D ai(0,0,0);
        ai.data()[ia] = i / nb;
        ai.data()[ib] = i % nb;
        if ( procNum(0,ai) != mype ) NAMD_bug("PmePencilMap is inconsistent");
        if ( ! msg ) NAMD_bug("PmePencilMap multiple pencils on a pe?");
        mgr->insertInitial(ai,msg);
        msg = 0;
      }
    }
    mgr->doneInserting();
    if ( msg ) CkFreeMsg(msg);
  }
private:
  const int ia, ib, nb, size;
  const int* const data;
  static int* newcopyint(int n, int *d) {
    int *newd = new int[n];
    memcpy(newd, d, n*sizeof(int));
    return newd;
  }
};

// use this idiom since messages don't have copy constructors
struct PmePencilInitMsgData {
  PmeGrid grid;
  int xBlocks, yBlocks, zBlocks;
  CProxy_PmeXPencil xPencil;
  CProxy_PmeYPencil yPencil;
  CProxy_PmeZPencil zPencil;
  CProxy_ComputePmeMgr pmeProxy;
  CProxy_NodePmeMgr pmeNodeProxy;
  CProxy_PmePencilMap xm;
  CProxy_PmePencilMap ym;
  CProxy_PmePencilMap zm;
};

class PmePencilInitMsg : public CMessage_PmePencilInitMsg {
public:
   PmePencilInitMsg(PmePencilInitMsgData &d) { data = d; }
   PmePencilInitMsgData data;
};


struct LocalPmeInfo {
  int nx, x_start;
  int ny_after_transpose, y_start_after_transpose;
};

struct NodePmeInfo {
  int npe, pe_start, real_node;
};


static int findRecipEvirPe() {
    PatchMap *patchMap = PatchMap::Object();
    {
      int mype = CkMyPe();
      if ( patchMap->numPatchesOnNode(mype) ) {
        return mype; 
      }
    }
    {
      int node = CmiMyNode();
      int firstpe = CmiNodeFirst(node);
      int nodeSize = CmiNodeSize(node);
      int myrank = CkMyRank();
      for ( int i=0; i<nodeSize; ++i ) {
        int pe = firstpe + (myrank+i)%nodeSize;
        if ( patchMap->numPatchesOnNode(pe) ) {
          return pe;
        }
      }
    }
    {
      int *pelist;
      int nodeSize;
      CmiGetPesOnPhysicalNode(CmiPhysicalNodeID(CkMyPe()), &pelist, &nodeSize);
      int myrank;
      for ( int i=0; i<nodeSize; ++i ) {
        if ( pelist[i] == CkMyPe() ) myrank = i;
      }
      for ( int i=0; i<nodeSize; ++i ) {
        int pe = pelist[(myrank+i)%nodeSize];
        if ( patchMap->numPatchesOnNode(pe) ) {
          return pe;
        }
      }
    }
    {
      int mype = CkMyPe();
      int npes = CkNumPes();
      for ( int i=0; i<npes; ++i ) {
        int pe = (mype+i)%npes;
        if ( patchMap->numPatchesOnNode(pe) ) {
          return pe;
        }
      }
    }
    NAMD_bug("findRecipEvirPe() failed!");
    return -999;  // should never happen
}


//Assigns gridPeMap and transPeMap to different set of processors.
void generatePmePeList2(int *gridPeMap, int numGridPes, int *transPeMap, int numTransPes){
  int ncpus = CkNumPes();
  
  for ( int i=0; i<numGridPes; ++i ) {
    gridPeMap[i] = WorkDistrib::peDiffuseOrdering[ncpus - numGridPes + i];
  }
  std::sort(gridPeMap,gridPeMap+numGridPes);
  int firstTransPe = ncpus - numGridPes - numTransPes;
  if ( firstTransPe < 0 ) {
    firstTransPe = 0;
    // 0 should be first in list, skip if possible
    if ( ncpus > numTransPes ) firstTransPe = 1;
  }
  for ( int i=0; i<numTransPes; ++i ) {
    transPeMap[i] = WorkDistrib::peDiffuseOrdering[firstTransPe + i];
  }
  std::sort(transPeMap,transPeMap+numTransPes);
}

#if USE_TOPOMAP 
//Topology aware PME allocation
bool generateBGLORBPmePeList(int *pemap, int numPes, int *block_pes=0, 
			     int nbpes=0);
#endif


int compare_bit_reversed(int a, int b) {
  int d = a ^ b;
  int c = 1;
  if ( d ) while ( ! (d & c) ) {
    c = c << 1;
  }
  return (a & c) - (b & c);
}

inline bool less_than_bit_reversed(int a, int b) {
  int d = a ^ b;
  int c = 1;
  if ( d ) while ( ! (d & c) ) {
    c = c << 1;
  }
  return d && (b & c);
}

struct sortop_bit_reversed {
  inline bool operator() (int a, int b) const {
    return less_than_bit_reversed(a,b);
  }
};

struct ijpair {
  int i,j;
  ijpair() {;}
  ijpair(int I, int J) : i(I), j(J) {;}
};

struct ijpair_sortop_bit_reversed {
  inline bool operator() (const ijpair &a, const ijpair &b) const {
    return ( less_than_bit_reversed(a.i,b.i)
             || ( (a.i == b.i) && less_than_bit_reversed(a.j,b.j) ) );
  }
};

class ComputePmeMgr : public CBase_ComputePmeMgr {
public:
  friend class ComputePme;
  friend class NodePmeMgr;
  ComputePmeMgr();
  ~ComputePmeMgr();

  void initialize(CkQdMsg*);
  void initialize_pencils(CkQdMsg*);
  void activate_pencils(CkQdMsg*);
  void recvArrays(CProxy_PmeXPencil, CProxy_PmeYPencil, CProxy_PmeZPencil);
  void initialize_computes();

  void sendData(Lattice &, int sequence);
  void sendDataPart(int first, int last, Lattice &, int sequence, int sourcepe, int errors);
  Lattice *sendDataHelper_lattice;
  int sendDataHelper_sequence;
  int sendDataHelper_sourcepe;
  int sendDataHelper_errors;
  void sendPencils(Lattice &, int sequence);
  void sendPencilsPart(int first, int last, Lattice &, int sequence, int sourcepe);
  void recvGrid(PmeGridMsg *);
  void gridCalc1(void);
  void sendTransBarrier(void);
  void sendTrans(void);
  void fwdSharedTrans(PmeTransMsg *);
  void recvSharedTrans(PmeSharedTransMsg *);
  void sendDataHelper(int);
  void sendPencilsHelper(int);
  void recvTrans(PmeTransMsg *);
  void procTrans(PmeTransMsg *);
  void gridCalc2(void);
  #ifdef OPENATOM_VERSION
  void gridCalc2Moa(void);
  #endif // OPENATOM_VERSION
  void gridCalc2R(void);
  void fwdSharedUntrans(PmeUntransMsg *);
  void recvSharedUntrans(PmeSharedUntransMsg *);
  void sendUntrans(void);
  void recvUntrans(PmeUntransMsg *);
  void procUntrans(PmeUntransMsg *);
  void gridCalc3(void);
  void sendUngrid(void);
  void recvUngrid(PmeGridMsg *);
  void recvAck(PmeAckMsg *);
  void copyResults(PmeGridMsg *);
  void copyPencils(PmeGridMsg *);
  void ungridCalc(void);
  void recvRecipEvir(PmeEvirMsg *);
  void addRecipEvirClient(void);
  void submitReductions();

#if 0 && USE_PERSISTENT
  void setup_recvgrid_persistent();
#endif

  static CmiNodeLock fftw_plan_lock;
  CmiNodeLock pmemgr_lock;  // for accessing this object from other threads

#ifdef NAMD_CUDA
  float *a_data_host;
  float *a_data_dev;
  float *f_data_host;
  float *f_data_dev;
  int cuda_atoms_count;
  int cuda_atoms_alloc;
  static CmiNodeLock cuda_lock;
  void chargeGridSubmitted(Lattice &lattice, int sequence);
  cudaEvent_t end_charges;
  cudaEvent_t *end_forces;
  int forces_count;
  int forces_done_count;
  double charges_time;
  double forces_time;
  int check_charges_count;
  int check_forces_count;
  int master_pe;
  int this_pe;

  void cuda_submit_charges(Lattice &lattice, int sequence);
  struct cuda_submit_charges_args {
    ComputePmeMgr *mgr; Lattice *lattice; int sequence;
  };
  static std::deque<cuda_submit_charges_args> cuda_submit_charges_deque;
  static bool cuda_busy;

  int chargeGridSubmittedCount;
  void sendChargeGridReady();
#endif
  Lattice *saved_lattice;  // saved by chargeGridSubmitted
  int saved_sequence;      // saved by chargeGridSubmitted
  void pollChargeGridReady();
  void pollForcesReady();
  void recvChargeGridReady();
  void chargeGridReady(Lattice &lattice, int sequence);

  ResizeArray<ComputePme*> pmeComputes;

private:

#if 0 && USE_PERSISTENT
  PersistentHandle   *recvGrid_handle;
#endif

  CProxy_ComputePmeMgr pmeProxy;
  CProxy_ComputePmeMgr pmeProxyDir;
  CProxy_NodePmeMgr pmeNodeProxy;
  NodePmeMgr *nodePmeMgr;
  ComputePmeMgr *masterPmeMgr;
  
  void addCompute(ComputePme *c) {
    if ( ! pmeComputes.size() ) initialize_computes();
    pmeComputes.add(c);
    c->setMgr(this);
  }

  ResizeArray<ComputePme*> heldComputes;
  PmeGrid myGrid;
  Lattice lattice;
  PmeKSpace *myKSpace;
  float *qgrid;
  float *kgrid;

#ifdef NAMD_FFTW
#ifdef NAMD_FFTW_3
  fftwf_plan *forward_plan_x, *backward_plan_x;
  fftwf_plan *forward_plan_yz, *backward_plan_yz;
  fftwf_complex *work;
#else
  fftw_plan forward_plan_x, backward_plan_x;
  rfftwnd_plan forward_plan_yz, backward_plan_yz;
  fftw_complex *work;
#endif
#else
  float *work;
#endif

  int qsize, fsize, bsize;
  int alchOn, alchFepOn, alchThermIntOn, lesOn, lesFactor, pairOn, selfOn, numGrids;
  int alchDecouple;
  int offload;
  BigReal alchElecLambdaStart;
  BigReal alchLambda;  // set on each step in ComputePme::ungridForces()

  float **q_arr;
  // q_list and q_count not used for offload
  float **q_list;
  int q_count;
  char *f_arr;
  char *fz_arr;
  PmeReduction evir[PME_MAX_EVALS];
  SubmitReduction *reduction;

  int noWorkCount;
  int doWorkCount;
  int ungridForcesCount;

#ifdef NAMD_CUDA
#define NUM_STREAMS 1
  cudaStream_t streams[NUM_STREAMS];
  int stream;

  float **q_arr_dev;
  float **v_arr_dev;
  float *q_data_host;
  float *q_data_dev;
  float *v_data_dev;
  int *ffz_host;
  int *ffz_dev;
  int q_data_size;
  int ffz_size;

  int f_data_mgr_alloc;
  float *f_data_mgr_host;
  float *f_data_mgr_dev;
  float **afn_host;
  float **afn_dev;

  float *bspline_coeffs_dev;
  float *bspline_dcoeffs_dev;
#endif
  int recipEvirCount;   // used in compute only
  int recipEvirClients; // used in compute only
  int recipEvirPe;      // used in trans only
  
  LocalPmeInfo *localInfo;
  NodePmeInfo *gridNodeInfo;
  NodePmeInfo *transNodeInfo;
  int qgrid_size;
  int qgrid_start;
  int qgrid_len;
  int fgrid_start;
  int fgrid_len;

  int numSources;
  int numGridPes;
  int numTransPes;
  int numGridNodes;
  int numTransNodes;
  int numDestRecipPes;
  int myGridPe, myGridNode;
  int myTransPe, myTransNode;
  int *gridPeMap;
  int *transPeMap;
  int *recipPeDest;
  int *gridPeOrder;
  int *gridNodeOrder;
  int *transNodeOrder;
  int grid_count;
  int trans_count;
  int untrans_count;
  int ungrid_count;
  PmeGridMsg **gridmsg_reuse;
  PmeReduction recip_evir2[PME_MAX_EVALS];

  int compute_sequence;  // set from patch computes, used for priorities
  int grid_sequence;  // set from grid messages, used for priorities
  int useBarrier;
  int sendTransBarrier_received;

  int usePencils;
  int xBlocks, yBlocks, zBlocks;
  CProxy_PmeXPencil xPencil;
  CProxy_PmeYPencil yPencil;
  CProxy_PmeZPencil zPencil;
  char *pencilActive;
  ijpair *activePencils;
  int numPencilsActive;
  int strayChargeErrors;
};

ResizeArray<ComputePme*>& getComputes(ComputePmeMgr *mgr) {
    return mgr->pmeComputes ;
}

  CmiNodeLock ComputePmeMgr::fftw_plan_lock;
#ifdef NAMD_CUDA
  CmiNodeLock ComputePmeMgr::cuda_lock;
  std::deque<ComputePmeMgr::cuda_submit_charges_args> ComputePmeMgr::cuda_submit_charges_deque;
  bool ComputePmeMgr::cuda_busy;
#endif

int isPmeProcessor(int p){ 
  SimParameters *simParams = Node::Object()->simParameters;
  if (simParams->usePMECUDA) {
    return 0;
  } else {
    return pencilPMEProcessors[p];
  }
}

class NodePmeMgr : public CBase_NodePmeMgr {
public:
  friend class ComputePmeMgr;
  friend class ComputePme;
  NodePmeMgr();
  ~NodePmeMgr();
  void initialize();
  void sendDataHelper(int);
  void sendPencilsHelper(int);
  void recvTrans(PmeTransMsg *);
  void recvUntrans(PmeUntransMsg *);
  void registerXPencil(CkArrayIndex3D, PmeXPencil *);
  void registerYPencil(CkArrayIndex3D, PmeYPencil *);
  void registerZPencil(CkArrayIndex3D, PmeZPencil *);
  void recvXTrans(PmeTransMsg *);
  void recvYTrans(PmeTransMsg *);
  void recvYUntrans(PmeUntransMsg *);
  void recvZGrid(PmeGridMsg *);
  void recvZUntrans(PmeUntransMsg *);

  void recvUngrid(PmeGridMsg *);

  void recvPencilMapProxies(CProxy_PmePencilMap _xm, CProxy_PmePencilMap _ym, CProxy_PmePencilMap _zm){
      xm=_xm; ym=_ym; zm=_zm;
  }
  CProxy_PmePencilMap xm;
  CProxy_PmePencilMap ym;
  CProxy_PmePencilMap zm;

private:
  CProxy_ComputePmeMgr mgrProxy;
  ComputePmeMgr *mgrObject;
  ComputePmeMgr **mgrObjects;
#ifdef NAMD_CUDA
  ComputePmeMgr *masterPmeMgr;
  int master_pe;
#endif
  CProxy_PmeXPencil xPencil;
  CProxy_PmeYPencil yPencil;
  CProxy_PmeZPencil zPencil;
  CkHashtableT<CkArrayIndex3D,PmeXPencil*> xPencilObj;
  CkHashtableT<CkArrayIndex3D,PmeYPencil*> yPencilObj;
  CkHashtableT<CkArrayIndex3D,PmeZPencil*> zPencilObj;  

#ifdef NAMD_CUDA
  cudaEvent_t end_charge_memset;
  cudaEvent_t end_all_pme_kernels;
  cudaEvent_t end_potential_memcpy;
#endif
};

NodePmeMgr::NodePmeMgr() {
  mgrObjects = new ComputePmeMgr*[CkMyNodeSize()];
}

NodePmeMgr::~NodePmeMgr() {
  delete [] mgrObjects;
}

void NodePmeMgr::initialize() {
  CProxy_ComputePmeMgr proxy = CkpvAccess(BOCclass_group).computePmeMgr;
  mgrObjects[CkMyRank()] = proxy.ckLocalBranch();
  if ( CkMyRank() == 0 ) {
    mgrProxy = proxy;
    mgrObject = proxy.ckLocalBranch();
  }
}

void NodePmeMgr::recvTrans(PmeTransMsg *msg) {
  mgrObject->fwdSharedTrans(msg);
}

void NodePmeMgr::recvUntrans(PmeUntransMsg *msg) {
  mgrObject->fwdSharedUntrans(msg);
}

void NodePmeMgr::recvUngrid(PmeGridMsg *msg) {
#ifdef NAMD_CUDA
  masterPmeMgr->recvUngrid(msg);
#else
  NAMD_bug("NodePmeMgr::recvUngrid called in non-CUDA build.");
#endif
}

void NodePmeMgr::registerXPencil(CkArrayIndex3D idx, PmeXPencil *obj)
{
  CmiLock(ComputePmeMgr::fftw_plan_lock);
  xPencilObj.put(idx)=obj;
  CmiUnlock(ComputePmeMgr::fftw_plan_lock);
}
void NodePmeMgr::registerYPencil(CkArrayIndex3D idx, PmeYPencil *obj)
{
  CmiLock(ComputePmeMgr::fftw_plan_lock);
  yPencilObj.put(idx)=obj;
  CmiUnlock(ComputePmeMgr::fftw_plan_lock);
}
void NodePmeMgr::registerZPencil(CkArrayIndex3D idx, PmeZPencil *obj)
{
  CmiLock(ComputePmeMgr::fftw_plan_lock);
  zPencilObj.put(idx)=obj;
  CmiUnlock(ComputePmeMgr::fftw_plan_lock);
}

ComputePmeMgr::ComputePmeMgr() : pmeProxy(thisgroup), 
				 pmeProxyDir(thisgroup) {

  CkpvAccess(BOCclass_group).computePmeMgr = thisgroup;
  pmeNodeProxy = CkpvAccess(BOCclass_group).nodePmeMgr;
  nodePmeMgr = pmeNodeProxy[CkMyNode()].ckLocalBranch();

  pmeNodeProxy.ckLocalBranch()->initialize();

  if ( CmiMyRank() == 0 ) {
    fftw_plan_lock = CmiCreateLock();
  }
  pmemgr_lock = CmiCreateLock();

  myKSpace = 0;
  kgrid = 0;
  work = 0;
  grid_count = 0;
  trans_count = 0;
  untrans_count = 0;
  ungrid_count = 0;
  gridmsg_reuse= new PmeGridMsg*[CkNumPes()];
  useBarrier = 0;
  sendTransBarrier_received = 0;
  usePencils = 0;

#ifdef NAMD_CUDA
 // offload has not been set so this happens on every run
  if ( CmiMyRank() == 0 ) {
    cuda_lock = CmiCreateLock();
  }

#if CUDA_VERSION >= 5050
  int leastPriority, greatestPriority;
  cudaDeviceGetStreamPriorityRange(&leastPriority, &greatestPriority);
  cuda_errcheck("in cudaDeviceGetStreamPriorityRange");
  //if ( CkMyNode() == 0 ) {
  //  CkPrintf("Pe %d PME CUDA stream priority range %d %d\n", CkMyPe(), leastPriority, greatestPriority);
  //}
#define CUDA_STREAM_CREATE(X) cudaStreamCreateWithPriority(X,cudaStreamDefault,greatestPriority)
#else
#define CUDA_STREAM_CREATE(X) cudaStreamCreate(X)
#endif

  stream = 0;
  for ( int i=0; i<NUM_STREAMS; ++i ) {
#if 1
    CUDA_STREAM_CREATE(&streams[i]);
    cuda_errcheck("cudaStreamCreate");
#else
  streams[i] = 0;  // XXXX Testing!!!
#endif
  }

  this_pe = CkMyPe();
 
  cudaEventCreateWithFlags(&end_charges,cudaEventDisableTiming);
  end_forces = 0;
  check_charges_count = 0;
  check_forces_count = 0;
  chargeGridSubmittedCount = 0;

  cuda_atoms_count = 0;
  cuda_atoms_alloc = 0;

  f_data_mgr_alloc = 0;
  f_data_mgr_host = 0;
  f_data_mgr_dev = 0;
  afn_host = 0;
  afn_dev = 0;

#define CUDA_EVENT_ID_PME_CHARGES 80
#define CUDA_EVENT_ID_PME_FORCES 81
#define CUDA_EVENT_ID_PME_TICK 82
#define CUDA_EVENT_ID_PME_COPY 83
#define CUDA_EVENT_ID_PME_KERNEL 84
  if ( 0 == CkMyPe() ) {
    traceRegisterUserEvent("CUDA PME charges", CUDA_EVENT_ID_PME_CHARGES);
    traceRegisterUserEvent("CUDA PME forces", CUDA_EVENT_ID_PME_FORCES);
    traceRegisterUserEvent("CUDA PME tick", CUDA_EVENT_ID_PME_TICK);
    traceRegisterUserEvent("CUDA PME memcpy", CUDA_EVENT_ID_PME_COPY);
    traceRegisterUserEvent("CUDA PME kernel", CUDA_EVENT_ID_PME_KERNEL);
  }
#endif
  recipEvirCount = 0;
  recipEvirClients = 0;
  recipEvirPe = -999;
}


void ComputePmeMgr::recvArrays(
	CProxy_PmeXPencil x, CProxy_PmeYPencil y, CProxy_PmeZPencil z) {
  xPencil = x;  yPencil = y;  zPencil = z;
  
    if(CmiMyRank()==0)
    {
      pmeNodeProxy.ckLocalBranch()->xPencil=x;
      pmeNodeProxy.ckLocalBranch()->yPencil=y;
      pmeNodeProxy.ckLocalBranch()->zPencil=z;
    }
}

#if USE_TOPO_SFC
 struct Coord
  {
    int x, y, z;
    Coord(): x(0), y(0), z(0) {}
    Coord(int a, int b, int c): x(a), y(b), z(c) {}
  };
  extern void SFC_grid(int xdim, int ydim, int zdim, int xdim1, int ydim1, int zdim1, vector<Coord> &result);

  void sort_sfc(SortableResizeArray<int> &procs, TopoManager &tmgr, vector<Coord> &result)
  {
     SortableResizeArray<int> newprocs(procs.size());
     int num = 0;
     for (int i=0; i<result.size(); i++) {
       Coord &c = result[i];
       for (int j=0; j<procs.size(); j++) {
         int pe = procs[j];
         int x,y,z,t;
         tmgr.rankToCoordinates(pe, x, y, z, t);    
         if (x==c.x && y==c.y && z==c.z)
           newprocs[num++] = pe;
       }
     } 
     CmiAssert(newprocs.size() == procs.size());
     procs = newprocs;
  }

  int find_level_grid(int x) 
  {
     int a = sqrt(x);
     int b;
     for (; a>0; a--) {
       if (x%a == 0) break;
     }
     if (a==1) a = x;
     b = x/a;
     //return a>b?a:b;
     return b;
  }
  CmiNodeLock tmgr_lock;
#endif

void Pme_init()
{
#if USE_TOPO_SFC
  if (CkMyRank() == 0) 
    tmgr_lock = CmiCreateLock();
#endif
}

void ComputePmeMgr::initialize(CkQdMsg *msg) {
  delete msg;

  localInfo = new LocalPmeInfo[CkNumPes()];
  gridNodeInfo = new NodePmeInfo[CkNumNodes()];
  transNodeInfo = new NodePmeInfo[CkNumNodes()];
  gridPeMap = new int[CkNumPes()];
  transPeMap = new int[CkNumPes()];
  recipPeDest = new int[CkNumPes()];
  gridPeOrder = new int[CkNumPes()];
  gridNodeOrder = new int[CkNumNodes()];
  transNodeOrder = new int[CkNumNodes()];

  if (CkMyRank() == 0) {
    pencilPMEProcessors = new char [CkNumPes()];
    memset (pencilPMEProcessors, 0, sizeof(char) * CkNumPes());
  }

  SimParameters *simParams = Node::Object()->simParameters;
  PatchMap *patchMap = PatchMap::Object();

  offload = simParams->PMEOffload;
#ifdef NAMD_CUDA
  if ( offload && ! deviceCUDA->one_device_per_node() ) {
    NAMD_die("PME offload requires exactly one CUDA device per process.  Use \"PMEOffload no\".");
  }
  if ( offload ) {
    int dev;
    cudaGetDevice(&dev);
    cuda_errcheck("in cudaGetDevice");
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, dev);
    cuda_errcheck("in cudaGetDeviceProperties");
    if ( deviceProp.major < 2 )
      NAMD_die("PME offload requires CUDA device of compute capability 2.0 or higher.  Use \"PMEOffload no\".");
  }
#endif

  alchLambda = -1.;  // illegal value to catch if not updated

  alchOn = simParams->alchOn;
  alchFepOn = simParams->alchFepOn;
  alchThermIntOn = simParams->alchThermIntOn;
  alchDecouple = alchOn && simParams->alchDecouple;
  alchElecLambdaStart = alchOn ? simParams->alchElecLambdaStart : 0;
  if (alchOn) {
    numGrids = 2;
    if (alchDecouple) numGrids += 2;
    if (alchElecLambdaStart || alchThermIntOn) numGrids ++;
  }
  else numGrids = 1;
  lesOn = simParams->lesOn;
  useBarrier = simParams->PMEBarrier;
  if ( lesOn ) {
    lesFactor = simParams->lesFactor;
    numGrids = lesFactor;
  }
  selfOn = 0;
  pairOn = simParams->pairInteractionOn;
  if ( pairOn ) {
    selfOn = simParams->pairInteractionSelf;
    if ( selfOn ) pairOn = 0;  // make pairOn and selfOn exclusive
    numGrids = selfOn ? 1 : 3;
  }

  if ( numGrids != 1 || simParams->PMEPencils == 0 ) usePencils = 0;
  else if ( simParams->PMEPencils > 0 ) usePencils = 1;
  else {
    int nrps = simParams->PMEProcessors;
    if ( nrps <= 0 ) nrps = CkNumPes();
    if ( nrps > CkNumPes() ) nrps = CkNumPes();
    int dimx = simParams->PMEGridSizeX;
    int dimy = simParams->PMEGridSizeY;
    int maxslabs = 1 + (dimx - 1) / simParams->PMEMinSlices;
    if ( maxslabs > nrps ) maxslabs = nrps;
    int maxpencils = ( simParams->PMEGridSizeX * simParams->PMEGridSizeY
		* simParams->PMEGridSizeZ ) / simParams->PMEMinPoints;
    if ( maxpencils > nrps ) maxpencils = nrps;
    if ( maxpencils > 3 * maxslabs ) usePencils = 1;
    else usePencils = 0;
  }

  if ( usePencils ) {
    int nrps = simParams->PMEProcessors;
    if ( nrps <= 0 ) nrps = CkNumPes();
    if ( nrps > CkNumPes() ) nrps = CkNumPes();
    if ( simParams->PMEPencils > 1 &&
         simParams->PMEPencils * simParams->PMEPencils <= nrps ) {
      xBlocks = yBlocks = zBlocks = simParams->PMEPencils;
    } else {
      int nb2 = ( simParams->PMEGridSizeX * simParams->PMEGridSizeY
		* simParams->PMEGridSizeZ ) / simParams->PMEMinPoints;
      if ( nb2 > nrps ) nb2 = nrps;
      if ( nb2 < 1 ) nb2 = 1;
      int nb = (int) sqrt((float)nb2);
      if ( nb < 1 ) nb = 1;
      xBlocks = zBlocks = nb;
      yBlocks = nb2 / nb;
    }

    if ( simParams->PMEPencilsX > 0 ) xBlocks = simParams->PMEPencilsX;
    if ( simParams->PMEPencilsY > 0 ) yBlocks = simParams->PMEPencilsY;
    if ( simParams->PMEPencilsZ > 0 ) zBlocks = simParams->PMEPencilsZ;

    int dimx = simParams->PMEGridSizeX;
    int bx = 1 + ( dimx - 1 ) / xBlocks;
    xBlocks = 1 + ( dimx - 1 ) / bx;

    int dimy = simParams->PMEGridSizeY;
    int by = 1 + ( dimy - 1 ) / yBlocks;
    yBlocks = 1 + ( dimy - 1 ) / by;

    int dimz = simParams->PMEGridSizeZ / 2 + 1;  // complex
    int bz = 1 + ( dimz - 1 ) / zBlocks;
    zBlocks = 1 + ( dimz - 1 ) / bz;

    if ( xBlocks * yBlocks > CkNumPes() ) {
      NAMD_die("PME pencils xBlocks * yBlocks > numPes");
    }
    if ( xBlocks * zBlocks > CkNumPes() ) {
      NAMD_die("PME pencils xBlocks * zBlocks > numPes");
    }
    if ( yBlocks * zBlocks > CkNumPes() ) {
      NAMD_die("PME pencils yBlocks * zBlocks > numPes");
    }

    if ( ! CkMyPe() ) {
      iout << iINFO << "PME using " << xBlocks << " x " <<
        yBlocks << " x " << zBlocks <<
        " pencil grid for FFT and reciprocal sum.\n" << endi;
    }
  } else { // usePencils

  {  // decide how many pes to use for reciprocal sum

    // rules based on work available
    int minslices = simParams->PMEMinSlices;
    int dimx = simParams->PMEGridSizeX;
    int nrpx = ( dimx + minslices - 1 ) / minslices;
    int dimy = simParams->PMEGridSizeY;
    int nrpy = ( dimy + minslices - 1 ) / minslices;

    // rules based on processors available
    int nrpp = CkNumPes();
    // if ( nrpp > 32 ) nrpp = 32;  // cap to limit messages
    if ( nrpp < nrpx ) nrpx = nrpp;
    if ( nrpp < nrpy ) nrpy = nrpp;

    // user override
    int nrps = simParams->PMEProcessors;
    if ( nrps > CkNumPes() ) nrps = CkNumPes();
    if ( nrps > 0 ) nrpx = nrps;
    if ( nrps > 0 ) nrpy = nrps;

    // make sure there aren't any totally empty processors
    int bx = ( dimx + nrpx - 1 ) / nrpx;
    nrpx = ( dimx + bx - 1 ) / bx;
    int by = ( dimy + nrpy - 1 ) / nrpy;
    nrpy = ( dimy + by - 1 ) / by;
    if ( bx != ( dimx + nrpx - 1 ) / nrpx )
      NAMD_bug("Error in selecting number of PME processors.");
    if ( by != ( dimy + nrpy - 1 ) / nrpy )
      NAMD_bug("Error in selecting number of PME processors.");

    numGridPes = nrpx;
    numTransPes = nrpy;
  }
  if ( ! CkMyPe() ) {
    iout << iINFO << "PME using " << numGridPes << " and " << numTransPes <<
      " processors for FFT and reciprocal sum.\n" << endi;
  }

  int sum_npes = numTransPes + numGridPes;
  int max_npes = (numTransPes > numGridPes)?numTransPes:numGridPes;

#if 0 // USE_TOPOMAP
  /* This code is being disabled permanently for slab PME on Blue Gene machines */
  PatchMap * pmap = PatchMap::Object();
  
  int patch_pes = pmap->numNodesWithPatches();
  TopoManager tmgr;
  if(tmgr.hasMultipleProcsPerNode())
    patch_pes *= 2;

  bool done = false;
  if(CkNumPes() > 2*sum_npes + patch_pes) {    
    done = generateBGLORBPmePeList(transPeMap, numTransPes);
    done &= generateBGLORBPmePeList(gridPeMap, numGridPes, transPeMap, numTransPes);    
  }
  else 
    if(CkNumPes() > 2 *max_npes + patch_pes) {
      done = generateBGLORBPmePeList(transPeMap, max_npes);
      gridPeMap = transPeMap;
    }

  if (!done)
#endif
    {
      //generatePmePeList(transPeMap, max_npes);
      //gridPeMap = transPeMap;
      generatePmePeList2(gridPeMap, numGridPes, transPeMap, numTransPes);
    }
  
  if ( ! CkMyPe() ) {
    iout << iINFO << "PME GRID LOCATIONS:";
    int i;
    for ( i=0; i<numGridPes && i<10; ++i ) {
      iout << " " << gridPeMap[i];
    }
    if ( i < numGridPes ) iout << " ...";
    iout << "\n" << endi;
    iout << iINFO << "PME TRANS LOCATIONS:";
    for ( i=0; i<numTransPes && i<10; ++i ) {
      iout << " " << transPeMap[i];
    }
    if ( i < numTransPes ) iout << " ...";
    iout << "\n" << endi;
  }

  // sort based on nodes and physical nodes
  std::sort(gridPeMap,gridPeMap+numGridPes,WorkDistrib::pe_sortop_compact());

  myGridPe = -1;
  myGridNode = -1;
  int i = 0;
  int node = -1;
  int real_node = -1;
  for ( i=0; i<numGridPes; ++i ) {
    if ( gridPeMap[i] == CkMyPe() ) myGridPe = i;
    if (CkMyRank() == 0) pencilPMEProcessors[gridPeMap[i]] |= 1;
    int real_node_i = CkNodeOf(gridPeMap[i]);
    if ( real_node_i == real_node ) {
      gridNodeInfo[node].npe += 1;
    } else {
      real_node = real_node_i;
      ++node;
      gridNodeInfo[node].real_node = real_node;
      gridNodeInfo[node].pe_start = i;
      gridNodeInfo[node].npe = 1;
    }
    if ( CkMyNode() == real_node_i ) myGridNode = node;
  }
  numGridNodes = node + 1;
  myTransPe = -1;
  myTransNode = -1;
  node = -1;
  real_node = -1;
  for ( i=0; i<numTransPes; ++i ) {
    if ( transPeMap[i] == CkMyPe() ) myTransPe = i;
    if (CkMyRank() == 0) pencilPMEProcessors[transPeMap[i]] |= 2;
    int real_node_i = CkNodeOf(transPeMap[i]);
    if ( real_node_i == real_node ) {
      transNodeInfo[node].npe += 1;
    } else {
      real_node = real_node_i;
      ++node;
      transNodeInfo[node].real_node = real_node;
      transNodeInfo[node].pe_start = i;
      transNodeInfo[node].npe = 1;
    }
    if ( CkMyNode() == real_node_i ) myTransNode = node;
  }
  numTransNodes = node + 1;

  if ( ! CkMyPe() ) {
    iout << iINFO << "PME USING " << numGridNodes << " GRID NODES AND "
         << numTransNodes << " TRANS NODES\n" << endi;
  }

  { // generate random orderings for grid and trans messages
    int i;
    for ( i = 0; i < numGridPes; ++i ) {
      gridPeOrder[i] = i;
    }
    Random rand(CkMyPe());
    if ( myGridPe < 0 ) {
      rand.reorder(gridPeOrder,numGridPes);
    } else {  // self last
      gridPeOrder[myGridPe] = numGridPes-1;
      gridPeOrder[numGridPes-1] = myGridPe;
      rand.reorder(gridPeOrder,numGridPes-1);
    } 
    for ( i = 0; i < numGridNodes; ++i ) {
      gridNodeOrder[i] = i;
    }
    if ( myGridNode < 0 ) {
      rand.reorder(gridNodeOrder,numGridNodes);
    } else {  // self last
      gridNodeOrder[myGridNode] = numGridNodes-1;
      gridNodeOrder[numGridNodes-1] = myGridNode;
      rand.reorder(gridNodeOrder,numGridNodes-1);
    }
    for ( i = 0; i < numTransNodes; ++i ) {
      transNodeOrder[i] = i;
    }
    if ( myTransNode < 0 ) {
      rand.reorder(transNodeOrder,numTransNodes);
    } else {  // self last
      transNodeOrder[myTransNode] = numTransNodes-1;
      transNodeOrder[numTransNodes-1] = myTransNode;
      rand.reorder(transNodeOrder,numTransNodes-1);
    }
  }
  
  } // ! usePencils

  myGrid.K1 = simParams->PMEGridSizeX;
  myGrid.K2 = simParams->PMEGridSizeY;
  myGrid.K3 = simParams->PMEGridSizeZ;
  myGrid.order = simParams->PMEInterpOrder;
  myGrid.dim2 = myGrid.K2;
  myGrid.dim3 = 2 * (myGrid.K3/2 + 1);

  if ( ! usePencils ) {
    myGrid.block1 = ( myGrid.K1 + numGridPes - 1 ) / numGridPes;
    myGrid.block2 = ( myGrid.K2 + numTransPes - 1 ) / numTransPes;
    myGrid.block3 = myGrid.dim3 / 2;  // complex
  }

  if ( usePencils ) {
    myGrid.block1 = ( myGrid.K1 + xBlocks - 1 ) / xBlocks;
    myGrid.block2 = ( myGrid.K2 + yBlocks - 1 ) / yBlocks;
    myGrid.block3 = ( myGrid.K3/2 + 1 + zBlocks - 1 ) / zBlocks;  // complex


      int pe = 0;
      int x,y,z;

      		SortableResizeArray<int> zprocs(xBlocks*yBlocks);
      		SortableResizeArray<int> yprocs(xBlocks*zBlocks);
      		SortableResizeArray<int> xprocs(yBlocks*zBlocks);
      
      		// decide which pes to use by bit reversal and patch use
      		int i;
      		int ncpus = CkNumPes();
      		SortableResizeArray<int> patches, nopatches, pmeprocs;
      		PatchMap *pmap = PatchMap::Object();
      		for ( int icpu=0; icpu<ncpus; ++icpu ) {
        		int ri = WorkDistrib::peDiffuseOrdering[icpu];
        		if ( ri ) { // keep 0 for special case
          			if ( pmap->numPatchesOnNode(ri) ) patches.add(ri);
          			else nopatches.add(ri);
        		}
      		}

#if USE_RANDOM_TOPO
            Random rand(CkMyPe());
            int *tmp = new int[patches.size()];
            int nn = patches.size();
            for (i=0;i<nn;i++)  tmp[i] = patches[i];
            rand.reorder(tmp, nn);
            patches.resize(0);
            for (i=0;i<nn;i++)  patches.add(tmp[i]);
            delete [] tmp;
            tmp = new int[nopatches.size()];
            nn = nopatches.size();
            for (i=0;i<nn;i++)  tmp[i] = nopatches[i];
            rand.reorder(tmp, nn);
            nopatches.resize(0);
            for (i=0;i<nn;i++)  nopatches.add(tmp[i]);
            delete [] tmp;
#endif

      		// only use zero if it eliminates overloading or has patches
      		int useZero = 0;
      		int npens = xBlocks*yBlocks;
      		if ( npens % ncpus == 0 ) useZero = 1;
      		if ( npens == nopatches.size() + 1 ) useZero = 1;
      		npens += xBlocks*zBlocks;
      		if ( npens % ncpus == 0 ) useZero = 1;
      		if ( npens == nopatches.size() + 1 ) useZero = 1;
      		npens += yBlocks*zBlocks;
      		if ( npens % ncpus == 0 ) useZero = 1;
      		if ( npens == nopatches.size() + 1 ) useZero = 1;

      		// add nopatches then patches in reversed order
      		for ( i=nopatches.size()-1; i>=0; --i ) pmeprocs.add(nopatches[i]);
      		if ( useZero && ! pmap->numPatchesOnNode(0) ) pmeprocs.add(0);
      		for ( i=patches.size()-1; i>=0; --i ) pmeprocs.add(patches[i]);
      		if ( pmap->numPatchesOnNode(0) ) pmeprocs.add(0);
  
      		int npes = pmeprocs.size();
      		for ( i=0; i<xBlocks*yBlocks; ++i, ++pe ) zprocs[i] = pmeprocs[pe%npes];
		if ( i>1 && zprocs[0] == zprocs[i-1] ) zprocs[0] = 0;
#if !USE_RANDOM_TOPO
      		zprocs.sort();
#endif
      		for ( i=0; i<xBlocks*zBlocks; ++i, ++pe ) yprocs[i] = pmeprocs[pe%npes];
		if ( i>1 && yprocs[0] == yprocs[i-1] ) yprocs[0] = 0;
#if !USE_RANDOM_TOPO
      		yprocs.sort();
#endif
      for ( i=0; i<yBlocks*zBlocks; ++i, ++pe ) xprocs[i] = pmeprocs[pe%npes];
      if ( i>1 && xprocs[0] == xprocs[i-1] ) xprocs[0] = 0;
#if !USE_RANDOM_TOPO
      xprocs.sort();
#endif

#if USE_TOPO_SFC
  CmiLock(tmgr_lock);
  //{
  TopoManager tmgr;
  int xdim = tmgr.getDimNX();
  int ydim = tmgr.getDimNY();
  int zdim = tmgr.getDimNZ();
  int xdim1 = find_level_grid(xdim);
  int ydim1 = find_level_grid(ydim);
  int zdim1 = find_level_grid(zdim);
  if(CkMyPe() == 0)
      printf("xdim: %d %d %d, %d %d %d\n", xdim, ydim, zdim, xdim1, ydim1, zdim1);

  vector<Coord> result;
  SFC_grid(xdim, ydim, zdim, xdim1, ydim1, zdim1, result);
  sort_sfc(xprocs, tmgr, result);
  sort_sfc(yprocs, tmgr, result);
  sort_sfc(zprocs, tmgr, result);
  //}
  CmiUnlock(tmgr_lock);
#endif


		if(CkMyPe() == 0){  
	      iout << iINFO << "PME Z PENCIL LOCATIONS:";
          for ( i=0; i<zprocs.size() && i<10; ++i ) {
#if USE_TOPO_SFC
              int x,y,z,t;
              tmgr.rankToCoordinates(zprocs[i], x,y, z, t);
              iout << " " << zprocs[i] << "(" << x << " " << y << " " << z << ")";
#else
              iout << " " << zprocs[i];
#endif
          }
    	  if ( i < zprocs.size() ) iout << " ...";
	      iout << "\n" << endi;
		}

    if (CkMyRank() == 0) {
      for (pe=0, x = 0; x < xBlocks; ++x)
	for (y = 0; y < yBlocks; ++y, ++pe ) {
	  pencilPMEProcessors[zprocs[pe]] = 1;
	}
    }
     
		if(CkMyPe() == 0){  
	      iout << iINFO << "PME Y PENCIL LOCATIONS:";
    	  for ( i=0; i<yprocs.size() && i<10; ++i ) {
#if USE_TOPO_SFC
              int x,y,z,t;
              tmgr.rankToCoordinates(yprocs[i], x,y, z, t);
              iout << " " << yprocs[i] << "(" << x << " " << y << " " << z << ")";
#else
              iout << " " << yprocs[i];
#endif
          }
    	  if ( i < yprocs.size() ) iout << " ...";
	      iout << "\n" << endi;
		}

    if (CkMyRank() == 0) {
      for (pe=0, z = 0; z < zBlocks; ++z )
	for (x = 0; x < xBlocks; ++x, ++pe ) {
	  pencilPMEProcessors[yprocs[pe]] = 1;
	}
    }
    
		if(CkMyPe() == 0){  
      		iout << iINFO << "PME X PENCIL LOCATIONS:";
		    for ( i=0; i<xprocs.size() && i<10; ++i ) {
#if USE_TOPO_SFC
                int x,y,z,t;
                tmgr.rankToCoordinates(xprocs[i], x,y, z, t);
                iout << " " << xprocs[i] << "(" << x << "  " << y << " " << z << ")";
#else
                iout << " " << xprocs[i];
#endif
            }
      		if ( i < xprocs.size() ) iout << " ...";
      		iout << "\n" << endi;
		}

    if (CkMyRank() == 0) {
      for (pe=0, y = 0; y < yBlocks; ++y )	
	for (z = 0; z < zBlocks; ++z, ++pe ) {
	  pencilPMEProcessors[xprocs[pe]] = 1;
	}
    }
	

	// creating the pencil arrays
	if ( CkMyPe() == 0 ){
#if !USE_RANDOM_TOPO
	// std::sort(zprocs.begin(),zprocs.end(),WorkDistrib::pe_sortop_compact());
	WorkDistrib::sortPmePes(zprocs.begin(),xBlocks,yBlocks);
	std::sort(yprocs.begin(),yprocs.end(),WorkDistrib::pe_sortop_compact());
	std::sort(xprocs.begin(),xprocs.end(),WorkDistrib::pe_sortop_compact());
#endif
#if 1
        CProxy_PmePencilMap zm = CProxy_PmePencilMap::ckNew(0,1,yBlocks,xBlocks*yBlocks,zprocs.begin());
        CProxy_PmePencilMap ym;
        if ( simParams->PMEPencilsYLayout )
          ym = CProxy_PmePencilMap::ckNew(0,2,zBlocks,zBlocks*xBlocks,yprocs.begin()); // new
        else
          ym = CProxy_PmePencilMap::ckNew(2,0,xBlocks,zBlocks*xBlocks,yprocs.begin()); // old
        CProxy_PmePencilMap xm;
        if ( simParams->PMEPencilsXLayout )
          xm = CProxy_PmePencilMap::ckNew(2,1,yBlocks,yBlocks*zBlocks,xprocs.begin()); // new
        else
          xm = CProxy_PmePencilMap::ckNew(1,2,zBlocks,yBlocks*zBlocks,xprocs.begin()); // old
        pmeNodeProxy.recvPencilMapProxies(xm,ym,zm);
        CkArrayOptions zo(xBlocks,yBlocks,1);  zo.setMap(zm);
        CkArrayOptions yo(xBlocks,1,zBlocks);  yo.setMap(ym);
        CkArrayOptions xo(1,yBlocks,zBlocks);  xo.setMap(xm);
        zo.setAnytimeMigration(false);  zo.setStaticInsertion(true);
        yo.setAnytimeMigration(false);  yo.setStaticInsertion(true);
        xo.setAnytimeMigration(false);  xo.setStaticInsertion(true);
        zPencil = CProxy_PmeZPencil::ckNew(zo);  // (xBlocks,yBlocks,1);
        yPencil = CProxy_PmeYPencil::ckNew(yo);  // (xBlocks,1,zBlocks);
        xPencil = CProxy_PmeXPencil::ckNew(xo);  // (1,yBlocks,zBlocks);
#else
	zPencil = CProxy_PmeZPencil::ckNew();  // (xBlocks,yBlocks,1);
      	yPencil = CProxy_PmeYPencil::ckNew();  // (xBlocks,1,zBlocks);
      	xPencil = CProxy_PmeXPencil::ckNew();  // (1,yBlocks,zBlocks);

		for (pe=0, x = 0; x < xBlocks; ++x)
			for (y = 0; y < yBlocks; ++y, ++pe ) {
	  			zPencil(x,y,0).insert(zprocs[pe]);
			}
      	zPencil.doneInserting();

		for (pe=0, x = 0; x < xBlocks; ++x)
			for (z = 0; z < zBlocks; ++z, ++pe ) {
	  			yPencil(x,0,z).insert(yprocs[pe]);
			}
      	yPencil.doneInserting();


		for (pe=0, y = 0; y < yBlocks; ++y )	
			for (z = 0; z < zBlocks; ++z, ++pe ) {
	  			xPencil(0,y,z).insert(xprocs[pe]);
			}
		xPencil.doneInserting();     
#endif

		pmeProxy.recvArrays(xPencil,yPencil,zPencil);
		PmePencilInitMsgData msgdata;
		msgdata.grid = myGrid;
		msgdata.xBlocks = xBlocks;
		msgdata.yBlocks = yBlocks;
		msgdata.zBlocks = zBlocks;
		msgdata.xPencil = xPencil;
		msgdata.yPencil = yPencil;
		msgdata.zPencil = zPencil;
		msgdata.pmeProxy = pmeProxyDir;
        msgdata.pmeNodeProxy = pmeNodeProxy;
        msgdata.xm = xm;
        msgdata.ym = ym;
        msgdata.zm = zm;
		xPencil.init(new PmePencilInitMsg(msgdata));
		yPencil.init(new PmePencilInitMsg(msgdata));
		zPencil.init(new PmePencilInitMsg(msgdata));
	}

    return;  // continue in initialize_pencils() at next startup stage
  }


  int pe;
  int nx = 0;
  for ( pe = 0; pe < numGridPes; ++pe ) {
    localInfo[pe].x_start = nx;
    nx += myGrid.block1;
    if ( nx > myGrid.K1 ) nx = myGrid.K1;
    localInfo[pe].nx = nx - localInfo[pe].x_start;
  }
  int ny = 0;
  for ( pe = 0; pe < numTransPes; ++pe ) {
    localInfo[pe].y_start_after_transpose = ny;
    ny += myGrid.block2;
    if ( ny > myGrid.K2 ) ny = myGrid.K2;
    localInfo[pe].ny_after_transpose =
			ny - localInfo[pe].y_start_after_transpose;
  }

  {  // decide how many pes this node exchanges charges with

  PatchMap *patchMap = PatchMap::Object();
  Lattice lattice = simParams->lattice;
  BigReal sysdima = lattice.a_r().unit() * lattice.a();
  BigReal cutoff = simParams->cutoff;
  BigReal patchdim = simParams->patchDimension;
  int numPatches = patchMap->numPatches();
  int numNodes = CkNumPes();
  int *source_flags = new int[numNodes];
  int node;
  for ( node=0; node<numNodes; ++node ) {
    source_flags[node] = 0;
    recipPeDest[node] = 0;
  }

  // // make sure that we don't get ahead of ourselves on this node
  // if ( CkMyPe() < numPatches && myRecipPe >= 0 ) {
  //   source_flags[CkMyPe()] = 1;
  //   recipPeDest[myRecipPe] = 1;
  // }

  for ( int pid=0; pid < numPatches; ++pid ) {
    int pnode = patchMap->node(pid);
#ifdef NAMD_CUDA
    if ( offload ) pnode = CkNodeFirst(CkNodeOf(pnode));
#endif
    int shift1 = (myGrid.K1 + myGrid.order - 1)/2;
    BigReal minx = patchMap->min_a(pid);
    BigReal maxx = patchMap->max_a(pid);
    BigReal margina = 0.5 * ( patchdim - cutoff ) / sysdima;
    // min1 (max1) is smallest (largest) grid line for this patch
    int min1 = ((int) floor(myGrid.K1 * (minx - margina))) + shift1 - myGrid.order + 1;
    int max1 = ((int) floor(myGrid.K1 * (maxx + margina))) + shift1;
    for ( int i=min1; i<=max1; ++i ) {
      int ix = i;
      while ( ix >= myGrid.K1 ) ix -= myGrid.K1;
      while ( ix < 0 ) ix += myGrid.K1;
      // set source_flags[pnode] if this patch sends to our node
      if ( myGridPe >= 0 && ix >= localInfo[myGridPe].x_start &&
           ix < localInfo[myGridPe].x_start + localInfo[myGridPe].nx ) {
        source_flags[pnode] = 1;
      }
      // set dest_flags[] for node that our patch sends to
#ifdef NAMD_CUDA
      if ( offload ) {
        if ( pnode == CkNodeFirst(CkMyNode()) ) {
          recipPeDest[ix / myGrid.block1] = 1;
        }
      } else
#endif
      if ( pnode == CkMyPe() ) {
        recipPeDest[ix / myGrid.block1] = 1;
      }
    }
  }

  int numSourcesSamePhysicalNode = 0;
  numSources = 0;
  numDestRecipPes = 0;
  for ( node=0; node<numNodes; ++node ) {
    if ( source_flags[node] ) ++numSources;
    if ( recipPeDest[node] ) ++numDestRecipPes;
    if ( source_flags[node] && CmiPeOnSamePhysicalNode(node,CkMyPe()) ) ++numSourcesSamePhysicalNode;
  }

#if 0
  if ( numSources ) {
    CkPrintf("pe %5d pme %5d of %5d on same physical node\n",
            CkMyPe(), numSourcesSamePhysicalNode, numSources);
    iout << iINFO << "PME " << CkMyPe() << " sources:";
    for ( node=0; node<numNodes; ++node ) {
      if ( source_flags[node] ) iout << " " << node;
    }
    iout << "\n" << endi;
  }
#endif

  delete [] source_flags;

  // CkPrintf("PME on node %d has %d sources and %d destinations\n",
  //           CkMyPe(), numSources, numDestRecipPes);

  }  // decide how many pes this node exchanges charges with (end)

  ungrid_count = numDestRecipPes;

  sendTransBarrier_received = 0;

  if ( myGridPe < 0 && myTransPe < 0 ) return;
  // the following only for nodes doing reciprocal sum

  if ( myTransPe >= 0 ) {
    recipEvirPe = findRecipEvirPe();
    pmeProxy[recipEvirPe].addRecipEvirClient();
  }

  if ( myTransPe >= 0 ) {
      int k2_start = localInfo[myTransPe].y_start_after_transpose;
      int k2_end = k2_start + localInfo[myTransPe].ny_after_transpose;
      #ifdef OPENATOM_VERSION
      if ( simParams->openatomOn ) { 
        CProxy_ComputeMoaMgr moaProxy(CkpvAccess(BOCclass_group).computeMoaMgr);
        myKSpace = new PmeKSpace(myGrid, k2_start, k2_end, 0, myGrid.dim3/2, moaProxy);
      } else {
        myKSpace = new PmeKSpace(myGrid, k2_start, k2_end, 0, myGrid.dim3/2);
      }
      #else  // OPENATOM_VERSION
      myKSpace = new PmeKSpace(myGrid, k2_start, k2_end, 0, myGrid.dim3/2);
      #endif // OPENATOM_VERSION
  }

  int local_size = myGrid.block1 * myGrid.K2 * myGrid.dim3;
  int local_size_2 = myGrid.block2 * myGrid.K1 * myGrid.dim3;
  if ( local_size < local_size_2 ) local_size = local_size_2;
  qgrid = new float[local_size*numGrids];
  if ( numGridPes > 1 || numTransPes > 1 ) {
    kgrid = new float[local_size*numGrids];
  } else {
    kgrid = qgrid;
  }
  qgrid_size = local_size;

  if ( myGridPe >= 0 ) {
  qgrid_start = localInfo[myGridPe].x_start * myGrid.K2 * myGrid.dim3;
  qgrid_len = localInfo[myGridPe].nx * myGrid.K2 * myGrid.dim3;
  fgrid_start = localInfo[myGridPe].x_start * myGrid.K2;
  fgrid_len = localInfo[myGridPe].nx * myGrid.K2;
  }

  int n[3]; n[0] = myGrid.K1; n[1] = myGrid.K2; n[2] = myGrid.K3;
#ifdef NAMD_FFTW
  CmiLock(fftw_plan_lock);
#ifdef NAMD_FFTW_3
  work = new fftwf_complex[n[0]];
  int fftwFlags = simParams->FFTWPatient ? FFTW_PATIENT  : simParams->FFTWEstimate ? FFTW_ESTIMATE  : FFTW_MEASURE ;
  if ( myGridPe >= 0 ) {
    forward_plan_yz=new fftwf_plan[numGrids];
    backward_plan_yz=new fftwf_plan[numGrids];
  }
  if ( myTransPe >= 0 ) {
    forward_plan_x=new fftwf_plan[numGrids];
    backward_plan_x=new fftwf_plan[numGrids];
  }
  /* need one plan per grid */
  if ( ! CkMyPe() ) iout << iINFO << "Optimizing 4 FFT steps.  1..." << endi;
  if ( myGridPe >= 0 ) {
    for( int g=0; g<numGrids; g++)
      {
	forward_plan_yz[g] = fftwf_plan_many_dft_r2c(2, n+1, 
						     localInfo[myGridPe].nx,
						     qgrid + qgrid_size * g,
						     NULL,
						     1,
						     myGrid.dim2 * myGrid.dim3,
						     (fftwf_complex *) 
						     (qgrid + qgrid_size * g),
						     NULL,
						     1,
						     myGrid.dim2 * (myGrid.dim3/2),
						     fftwFlags);
      }
  }
  int zdim = myGrid.dim3;
  int xStride=localInfo[myTransPe].ny_after_transpose *( myGrid.dim3 / 2);
  if ( ! CkMyPe() ) iout << " 2..." << endi;
  if ( myTransPe >= 0 ) {
    for( int g=0; g<numGrids; g++)
      {

	forward_plan_x[g] = fftwf_plan_many_dft(1, n, xStride,
						(fftwf_complex *)
						(kgrid+qgrid_size*g),
						NULL,
						xStride,
						1,
						(fftwf_complex *)
						(kgrid+qgrid_size*g),
						NULL,
						xStride,
						1,
						FFTW_FORWARD,fftwFlags);
	
      }
  }
  if ( ! CkMyPe() ) iout << " 3..." << endi;
  if ( myTransPe >= 0 ) {
    for( int g=0; g<numGrids; g++)
      {
	backward_plan_x[g] = fftwf_plan_many_dft(1, n, xStride,
						 (fftwf_complex *)
						 (kgrid+qgrid_size*g),
						 NULL,
						 xStride,
						 1,
						 (fftwf_complex *)
						 (kgrid+qgrid_size*g),
						 NULL,
						 xStride,
						 1,
						 FFTW_BACKWARD, fftwFlags);

      }
  }
  if ( ! CkMyPe() ) iout << " 4..." << endi;
  if ( myGridPe >= 0 ) {
    for( int g=0; g<numGrids; g++)
      {
	backward_plan_yz[g] = fftwf_plan_many_dft_c2r(2, n+1, 
						      localInfo[myGridPe].nx,
						      (fftwf_complex *)
						      (qgrid + qgrid_size * g),
						      NULL,
						      1,
						      myGrid.dim2*(myGrid.dim3/2),
						      qgrid + qgrid_size * g,
						      NULL,
						      1,
						      myGrid.dim2 * myGrid.dim3,
						      fftwFlags);
      }
  }
  if ( ! CkMyPe() ) iout << "   Done.\n" << endi;

#else
  work = new fftw_complex[n[0]];

  if ( ! CkMyPe() ) iout << iINFO << "Optimizing 4 FFT steps.  1..." << endi;
  if ( myGridPe >= 0 ) {
  forward_plan_yz = rfftwnd_create_plan_specific(2, n+1, FFTW_REAL_TO_COMPLEX,
	( simParams->FFTWEstimate ? FFTW_ESTIMATE : FFTW_MEASURE )
	| FFTW_IN_PLACE | FFTW_USE_WISDOM, qgrid, 1, 0, 0);
  }
  if ( ! CkMyPe() ) iout << " 2..." << endi;
  if ( myTransPe >= 0 ) {
      forward_plan_x = fftw_create_plan_specific(n[0], FFTW_FORWARD,
	( simParams->FFTWEstimate ? FFTW_ESTIMATE : FFTW_MEASURE )
	| FFTW_IN_PLACE | FFTW_USE_WISDOM, (fftw_complex *) kgrid,
	localInfo[myTransPe].ny_after_transpose * myGrid.dim3 / 2, work, 1);
  }
  if ( ! CkMyPe() ) iout << " 3..." << endi;
  if ( myTransPe >= 0 ) {
  backward_plan_x = fftw_create_plan_specific(n[0], FFTW_BACKWARD,
	( simParams->FFTWEstimate ? FFTW_ESTIMATE : FFTW_MEASURE )
	| FFTW_IN_PLACE | FFTW_USE_WISDOM, (fftw_complex *) kgrid,
	localInfo[myTransPe].ny_after_transpose * myGrid.dim3 / 2, work, 1);
  }
  if ( ! CkMyPe() ) iout << " 4..." << endi;
  if ( myGridPe >= 0 ) {
  backward_plan_yz = rfftwnd_create_plan_specific(2, n+1, FFTW_COMPLEX_TO_REAL,
	( simParams->FFTWEstimate ? FFTW_ESTIMATE : FFTW_MEASURE )
	| FFTW_IN_PLACE | FFTW_USE_WISDOM, qgrid, 1, 0, 0);
  }
  if ( ! CkMyPe() ) iout << "   Done.\n" << endi;
#endif
  CmiUnlock(fftw_plan_lock);
#else
  NAMD_die("Sorry, FFTW must be compiled in to use PME.");
#endif

  if ( myGridPe >= 0 && numSources == 0 )
		NAMD_bug("PME grid elements exist without sources.");
  grid_count = numSources;
  memset( (void*) qgrid, 0, qgrid_size * numGrids * sizeof(float) );
  trans_count = numGridPes;
}



void ComputePmeMgr::initialize_pencils(CkQdMsg *msg) {
  delete msg;
  if ( ! usePencils ) return;

  SimParameters *simParams = Node::Object()->simParameters;

  PatchMap *patchMap = PatchMap::Object();
  Lattice lattice = simParams->lattice;
  BigReal sysdima = lattice.a_r().unit() * lattice.a();
  BigReal sysdimb = lattice.b_r().unit() * lattice.b();
  BigReal cutoff = simParams->cutoff;
  BigReal patchdim = simParams->patchDimension;
  int numPatches = patchMap->numPatches();

  pencilActive = new char[xBlocks*yBlocks];
  for ( int i=0; i<xBlocks; ++i ) {
    for ( int j=0; j<yBlocks; ++j ) {
      pencilActive[i*yBlocks+j] = 0;
    }
  }

  for ( int pid=0; pid < numPatches; ++pid ) {
    int pnode = patchMap->node(pid);
#ifdef NAMD_CUDA
    if ( offload ) {
      if ( CkNodeOf(pnode) != CkMyNode() ) continue;
    } else
#endif
    if ( pnode != CkMyPe() ) continue;

    int shift1 = (myGrid.K1 + myGrid.order - 1)/2;
    int shift2 = (myGrid.K2 + myGrid.order - 1)/2;

    BigReal minx = patchMap->min_a(pid);
    BigReal maxx = patchMap->max_a(pid);
    BigReal margina = 0.5 * ( patchdim - cutoff ) / sysdima;
    // min1 (max1) is smallest (largest) grid line for this patch
    int min1 = ((int) floor(myGrid.K1 * (minx - margina))) + shift1 - myGrid.order + 1;
    int max1 = ((int) floor(myGrid.K1 * (maxx + margina))) + shift1;

    BigReal miny = patchMap->min_b(pid);
    BigReal maxy = patchMap->max_b(pid);
    BigReal marginb = 0.5 * ( patchdim - cutoff ) / sysdimb;
    // min2 (max2) is smallest (largest) grid line for this patch
    int min2 = ((int) floor(myGrid.K2 * (miny - marginb))) + shift2 - myGrid.order + 1;
    int max2 = ((int) floor(myGrid.K2 * (maxy + marginb))) + shift2;

    for ( int i=min1; i<=max1; ++i ) {
      int ix = i;
      while ( ix >= myGrid.K1 ) ix -= myGrid.K1;
      while ( ix < 0 ) ix += myGrid.K1;
      for ( int j=min2; j<=max2; ++j ) {
        int jy = j;
        while ( jy >= myGrid.K2 ) jy -= myGrid.K2;
        while ( jy < 0 ) jy += myGrid.K2;
        pencilActive[(ix / myGrid.block1)*yBlocks + (jy / myGrid.block2)] = 1;
      }
    }
  }

  numPencilsActive = 0;
  for ( int i=0; i<xBlocks; ++i ) {
    for ( int j=0; j<yBlocks; ++j ) {
      if ( pencilActive[i*yBlocks+j] ) {
        ++numPencilsActive;
#ifdef NAMD_CUDA
        if ( CkMyPe() == deviceCUDA->getMasterPe() || ! offload )
#endif
        zPencil(i,j,0).dummyRecvGrid(CkMyPe(),0);
      }
    }
  }
  activePencils = new ijpair[numPencilsActive];
  numPencilsActive = 0;
  for ( int i=0; i<xBlocks; ++i ) {
    for ( int j=0; j<yBlocks; ++j ) {
      if ( pencilActive[i*yBlocks+j] ) {
        activePencils[numPencilsActive++] = ijpair(i,j);
      }
    }
  }
  if ( simParams->PMESendOrder ) {
    std::sort(activePencils,activePencils+numPencilsActive,ijpair_sortop_bit_reversed());
  } else {
    Random rand(CkMyPe());
    rand.reorder(activePencils,numPencilsActive);
  }
  //if ( numPencilsActive ) {
  //  CkPrintf("node %d sending to %d pencils\n", CkMyPe(), numPencilsActive);
  //}

  ungrid_count = numPencilsActive;
}


void ComputePmeMgr::activate_pencils(CkQdMsg *msg) {
  if ( ! usePencils ) return;
  if ( CkMyPe() == 0 ) zPencil.dummyRecvGrid(CkMyPe(),1);
}


ComputePmeMgr::~ComputePmeMgr() {

  if ( CmiMyRank() == 0 ) {
    CmiDestroyLock(fftw_plan_lock);
  }
  CmiDestroyLock(pmemgr_lock);

  delete myKSpace;
  delete [] localInfo;
  delete [] gridNodeInfo;
  delete [] transNodeInfo;
  delete [] gridPeMap;
  delete [] transPeMap;
  delete [] recipPeDest;
  delete [] gridPeOrder;
  delete [] gridNodeOrder;
  delete [] transNodeOrder;
  delete [] qgrid;
  if ( kgrid != qgrid ) delete [] kgrid;
  delete [] work;
  delete [] gridmsg_reuse;

 if ( ! offload ) {
  for (int i=0; i<q_count; ++i) {
    delete [] q_list[i];
  }
  delete [] q_list;
  delete [] fz_arr;
 }
  delete [] f_arr;
  delete [] q_arr;
}

void ComputePmeMgr::recvGrid(PmeGridMsg *msg) {
  // CkPrintf("recvGrid from %d on Pe(%d)\n",msg->sourceNode,CkMyPe());
  if ( grid_count == 0 ) {
    NAMD_bug("Message order failure in ComputePmeMgr::recvGrid\n");
  }
  if ( grid_count == numSources ) {
    lattice = msg->lattice;
    grid_sequence = msg->sequence;
  }

  int zdim = myGrid.dim3;
  int zlistlen = msg->zlistlen;
  int *zlist = msg->zlist;
  float *qmsg = msg->qgrid;
  for ( int g=0; g<numGrids; ++g ) {
    char *f = msg->fgrid + fgrid_len * g;
    float *q = qgrid + qgrid_size * g;
    for ( int i=0; i<fgrid_len; ++i ) {
      if ( f[i] ) {
        for ( int k=0; k<zlistlen; ++k ) {
          q[zlist[k]] += *(qmsg++);
        }
      }
      q += zdim;
    }
  }

  gridmsg_reuse[numSources-grid_count] = msg;
  --grid_count;

  if ( grid_count == 0 ) {
    pmeProxyDir[CkMyPe()].gridCalc1();
    if ( useBarrier ) pmeProxyDir[0].sendTransBarrier();
  }
}
#ifdef MANUAL_DEBUG_FFTW3

/* utility functions for manual debugging */
void dumpMatrixFloat(const char *infilename, float *matrix, int xdim, int ydim, int zdim,int pe)
{

  char fmt[1000];
  char filename[1000];
  strncpy(fmt,infilename,999);
  strncat(fmt,"_%d.out",999);
  sprintf(filename,fmt, pe);
  FILE *loutfile = fopen(filename, "w");
#ifdef PAIRCALC_TEST_DUMP
  fprintf(loutfile,"%d\n",ydim);
#endif
  fprintf(loutfile,"%d %d %d\n",xdim,ydim, zdim);
  for(int i=0;i<xdim;i++)
    for(int j=0;j<ydim;j++)
      for(int k=0;k<zdim;k++)
	fprintf(loutfile,"%d %d %d %.8f\n",i,j,k,matrix[i*zdim*ydim+j*zdim +k]);
  fclose(loutfile);

}

void dumpMatrixFloat3(const char *infilename, float *matrix, int xdim, int ydim, int zdim,int x, int y, int z)
{
  char fmt[1000];
  char filename[1000];
  strncpy(fmt,infilename,999);
  strncat(fmt,"_%d_%d_%d.out",999);
  sprintf(filename,fmt, x,y,z);
  FILE *loutfile = fopen(filename, "w");
  CkAssert(loutfile!=NULL);
  CkPrintf("opened %s for dump\n",filename);
  fprintf(loutfile,"%d %d %d\n",xdim,ydim, zdim);
  for(int i=0;i<xdim;i++)
    for(int j=0;j<ydim;j++)
      for(int k=0;k<zdim;k++)
	fprintf(loutfile,"%d %d %d %.8f\n",i,j,k,matrix[i*zdim*ydim+j*zdim +k]);
  fclose(loutfile);
}

#endif

void ComputePmeMgr::gridCalc1(void) {
  // CkPrintf("gridCalc1 on Pe(%d)\n",CkMyPe());

#ifdef NAMD_FFTW
  for ( int g=0; g<numGrids; ++g ) {
#ifdef NAMD_FFTW_3
    fftwf_execute(forward_plan_yz[g]);
#else
    rfftwnd_real_to_complex(forward_plan_yz, localInfo[myGridPe].nx,
	qgrid + qgrid_size * g, 1, myGrid.dim2 * myGrid.dim3, 0, 0, 0);
#endif

  }
#endif

  if ( ! useBarrier ) pmeProxyDir[CkMyPe()].sendTrans();
}

void ComputePmeMgr::sendTransBarrier(void) {
  sendTransBarrier_received += 1;
  // CkPrintf("sendTransBarrier on %d %d\n",myGridPe,numGridPes-sendTransBarrier_received);
  if ( sendTransBarrier_received < numGridPes ) return;
  sendTransBarrier_received = 0;
  for ( int i=0; i<numGridPes; ++i ) {
    pmeProxyDir[gridPeMap[i]].sendTrans();
  }
}

void ComputePmeMgr::sendTrans(void) {
  // CkPrintf("sendTrans on Pe(%d)\n",CkMyPe());

  // send data for transpose
  int zdim = myGrid.dim3;
  int nx = localInfo[myGridPe].nx;
  int x_start = localInfo[myGridPe].x_start;
  int slicelen = myGrid.K2 * zdim;

  ComputePmeMgr **mgrObjects = pmeNodeProxy.ckLocalBranch()->mgrObjects;

#if CMK_BLUEGENEL
  CmiNetworkProgressAfter (0);
#endif

  for (int j=0; j<numTransNodes; j++) {
    int node = transNodeOrder[j];  // different order on each node
    int pe = transNodeInfo[node].pe_start;
    int npe = transNodeInfo[node].npe;
    int totlen = 0;
    if ( node != myTransNode ) for (int i=0; i<npe; ++i, ++pe) {
      LocalPmeInfo &li = localInfo[pe];
      int cpylen = li.ny_after_transpose * zdim;
      totlen += cpylen;
    }
    PmeTransMsg *newmsg = new (nx * totlen * numGrids,
				PRIORITY_SIZE) PmeTransMsg;
    newmsg->sourceNode = myGridPe;
    newmsg->lattice = lattice;
    newmsg->x_start = x_start;
    newmsg->nx = nx;
    for ( int g=0; g<numGrids; ++g ) {
      float *qmsg = newmsg->qgrid + nx * totlen * g;
      pe = transNodeInfo[node].pe_start;
      for (int i=0; i<npe; ++i, ++pe) {
        LocalPmeInfo &li = localInfo[pe];
        int cpylen = li.ny_after_transpose * zdim;
        if ( node == myTransNode ) {
          ComputePmeMgr *m = mgrObjects[CkRankOf(transPeMap[pe])];
          qmsg = m->kgrid + m->qgrid_size * g + x_start*cpylen;
        }
        float *q = qgrid + qgrid_size * g + li.y_start_after_transpose * zdim;
        for ( int x = 0; x < nx; ++x ) {
          CmiMemcpy((void*)qmsg, (void*)q, cpylen*sizeof(float));
          q += slicelen;
          qmsg += cpylen;
        }
      }
    }
    newmsg->sequence = grid_sequence;
    SET_PRIORITY(newmsg,grid_sequence,PME_TRANS_PRIORITY)
    if ( node == myTransNode ) newmsg->nx = 0;
    if ( npe > 1 ) {
      if ( node == myTransNode ) fwdSharedTrans(newmsg);
      else pmeNodeProxy[transNodeInfo[node].real_node].recvTrans(newmsg);
    } else pmeProxy[transPeMap[transNodeInfo[node].pe_start]].recvTrans(newmsg);
  }
 
  untrans_count = numTransPes;

}

void ComputePmeMgr::fwdSharedTrans(PmeTransMsg *msg) {
  // CkPrintf("fwdSharedTrans on Pe(%d)\n",CkMyPe());
  int pe = transNodeInfo[myTransNode].pe_start;
  int npe = transNodeInfo[myTransNode].npe;
  CmiNodeLock lock = CmiCreateLock();
  int *count = new int; *count = npe;
  for (int i=0; i<npe; ++i, ++pe) {
    PmeSharedTransMsg *shmsg = new (PRIORITY_SIZE) PmeSharedTransMsg;
    SET_PRIORITY(shmsg,msg->sequence,PME_TRANS_PRIORITY)
    shmsg->msg = msg;
    shmsg->count = count;
    shmsg->lock = lock;
    pmeProxy[transPeMap[pe]].recvSharedTrans(shmsg);
  }
}

void ComputePmeMgr::recvSharedTrans(PmeSharedTransMsg *msg) {
  procTrans(msg->msg);
  CmiLock(msg->lock);
  int count = --(*msg->count);
  CmiUnlock(msg->lock);
  if ( count == 0 ) {
    CmiDestroyLock(msg->lock);
    delete msg->count;
    delete msg->msg;
  }
  delete msg;
}

void ComputePmeMgr::recvTrans(PmeTransMsg *msg) {
  procTrans(msg);
  delete msg;
}

void ComputePmeMgr::procTrans(PmeTransMsg *msg) {
  // CkPrintf("procTrans on Pe(%d)\n",CkMyPe());
  if ( trans_count == numGridPes ) {
    lattice = msg->lattice;
    grid_sequence = msg->sequence;
  }

 if ( msg->nx ) {
  int zdim = myGrid.dim3;
  NodePmeInfo &nodeInfo(transNodeInfo[myTransNode]);
  int first_pe = nodeInfo.pe_start;
  int last_pe = first_pe+nodeInfo.npe-1;
  int y_skip = localInfo[myTransPe].y_start_after_transpose
             - localInfo[first_pe].y_start_after_transpose;
  int ny_msg = localInfo[last_pe].y_start_after_transpose
             + localInfo[last_pe].ny_after_transpose
             - localInfo[first_pe].y_start_after_transpose;
  int ny = localInfo[myTransPe].ny_after_transpose;
  int x_start = msg->x_start;
  int nx = msg->nx;
  for ( int g=0; g<numGrids; ++g ) {
    CmiMemcpy((void*)(kgrid + qgrid_size * g + x_start*ny*zdim),
	(void*)(msg->qgrid + nx*(ny_msg*g+y_skip)*zdim),
	nx*ny*zdim*sizeof(float));
  }
 }

  --trans_count;

  if ( trans_count == 0 ) {
    pmeProxyDir[CkMyPe()].gridCalc2();
  }
}

void ComputePmeMgr::gridCalc2(void) {
  // CkPrintf("gridCalc2 on Pe(%d)\n",CkMyPe());

#if CMK_BLUEGENEL
  CmiNetworkProgressAfter (0);
#endif

  int zdim = myGrid.dim3;
  // int y_start = localInfo[myTransPe].y_start_after_transpose;
  int ny = localInfo[myTransPe].ny_after_transpose;

  for ( int g=0; g<numGrids; ++g ) {
    // finish forward FFT (x dimension)
#ifdef NAMD_FFTW
#ifdef NAMD_FFTW_3
    fftwf_execute(forward_plan_x[g]);
#else
    fftw(forward_plan_x, ny * zdim / 2, (fftw_complex *)(kgrid+qgrid_size*g),
	ny * zdim / 2, 1, work, 1, 0);
#endif
#endif
  }

#ifdef OPENATOM_VERSION
    if ( ! simParams -> openatomOn ) { 
#endif // OPENATOM_VERSION
      gridCalc2R();
#ifdef OPENATOM_VERSION
    } else {
      gridCalc2Moa();
    }
#endif // OPENATOM_VERSION
}

#ifdef OPENATOM_VERSION
void ComputePmeMgr::gridCalc2Moa(void) {

  int zdim = myGrid.dim3;
  // int y_start = localInfo[myTransPe].y_start_after_transpose;
  int ny = localInfo[myTransPe].ny_after_transpose;

  SimParameters *simParams = Node::Object()->simParameters;

  CProxy_ComputeMoaMgr moaProxy(CkpvAccess(BOCclass_group).computeMoaMgr);

  for ( int g=0; g<numGrids; ++g ) {
    #ifdef OPENATOM_VERSION_DEBUG 
    CkPrintf("Sending recQ on processor %d \n", CkMyPe());
    for ( int i=0; i<=(ny * zdim / 2); ++i) 
    {
      CkPrintf("PE, g,fftw_q,k*q*g, kgrid, qgrid_size value %d pre-send = %d, %d, %f %f, %d, \n", i, CkMyPe(), g, (kgrid+qgrid_size*g)[i], kgrid[i], qgrid_size);
    }
    #endif // OPENATOM_VERSION_DEBUG
//     mqcpProxy[CkMyPe()].recvQ((ny * zdim / 2),((fftw_complex *)(kgrid+qgrid_size*g)));
    CkCallback resumePme(CkIndex_ComputePmeMgr::gridCalc2R(), thishandle);
    moaProxy[CkMyPe()].recvQ(g,numGrids,(ny * zdim / 2),(kgrid+qgrid_size*g), resumePme);
  }
}
#endif // OPENATOM_VERSION

void ComputePmeMgr::gridCalc2R(void) {

  int useCkLoop = 0;
#if CMK_SMP && USE_CKLOOP
  if ( Node::Object()->simParameters->useCkLoop >= CKLOOP_CTRL_PME_KSPACE
       && CkNumPes() >= 2 * numTransPes ) {
    useCkLoop = 1;
  }
#endif

  int zdim = myGrid.dim3;
  // int y_start = localInfo[myTransPe].y_start_after_transpose;
  int ny = localInfo[myTransPe].ny_after_transpose;

  for ( int g=0; g<numGrids; ++g ) {
    // reciprocal space portion of PME
    BigReal ewaldcof = ComputeNonbondedUtil::ewaldcof;
    recip_evir2[g][0] = myKSpace->compute_energy(kgrid+qgrid_size*g,
			lattice, ewaldcof, &(recip_evir2[g][1]), useCkLoop);
    // CkPrintf("Ewald reciprocal energy = %f\n", recip_evir2[g][0]);

    // start backward FFT (x dimension)

#ifdef NAMD_FFTW
#ifdef NAMD_FFTW_3
    fftwf_execute(backward_plan_x[g]);
#else
    fftw(backward_plan_x, ny * zdim / 2, (fftw_complex *)(kgrid+qgrid_size*g),
	ny * zdim / 2, 1, work, 1, 0);
#endif
#endif
  }
  
  pmeProxyDir[CkMyPe()].sendUntrans();
}

void ComputePmeMgr::sendUntrans(void) {

  int zdim = myGrid.dim3;
  int y_start = localInfo[myTransPe].y_start_after_transpose;
  int ny = localInfo[myTransPe].ny_after_transpose;
  int slicelen = myGrid.K2 * zdim;

  ComputePmeMgr **mgrObjects = pmeNodeProxy.ckLocalBranch()->mgrObjects;

  { // send energy and virial
    PmeEvirMsg *newmsg = new (numGrids, PRIORITY_SIZE) PmeEvirMsg;
    for ( int g=0; g<numGrids; ++g ) {
      newmsg->evir[g] = recip_evir2[g];
    }
    SET_PRIORITY(newmsg,grid_sequence,PME_UNGRID_PRIORITY)
    CmiEnableUrgentSend(1);
    pmeProxy[recipEvirPe].recvRecipEvir(newmsg);
    CmiEnableUrgentSend(0);
  }

#if CMK_BLUEGENEL
  CmiNetworkProgressAfter (0);
#endif

  // send data for reverse transpose
  for (int j=0; j<numGridNodes; j++) {
    int node = gridNodeOrder[j];  // different order on each node
    int pe = gridNodeInfo[node].pe_start;
    int npe = gridNodeInfo[node].npe;
    int totlen = 0;
    if ( node != myGridNode ) for (int i=0; i<npe; ++i, ++pe) {
      LocalPmeInfo &li = localInfo[pe];
      int cpylen = li.nx * zdim;
      totlen += cpylen;
    }
    PmeUntransMsg *newmsg = new (ny * totlen * numGrids, PRIORITY_SIZE) PmeUntransMsg;
    newmsg->sourceNode = myTransPe;
    newmsg->y_start = y_start;
    newmsg->ny = ny;
    for ( int g=0; g<numGrids; ++g ) {
      float *qmsg = newmsg->qgrid + ny * totlen * g;
      pe = gridNodeInfo[node].pe_start;
      for (int i=0; i<npe; ++i, ++pe) {
        LocalPmeInfo &li = localInfo[pe];
        if ( node == myGridNode ) {
          ComputePmeMgr *m = mgrObjects[CkRankOf(gridPeMap[pe])];
          qmsg = m->qgrid + m->qgrid_size * g + y_start * zdim;
          float *q = kgrid + qgrid_size*g + li.x_start*ny*zdim;
          int cpylen = ny * zdim;
          for ( int x = 0; x < li.nx; ++x ) {
            CmiMemcpy((void*)qmsg, (void*)q, cpylen*sizeof(float));
            q += cpylen;
            qmsg += slicelen;
          }
        } else {
          CmiMemcpy((void*)qmsg,
		(void*)(kgrid + qgrid_size*g + li.x_start*ny*zdim),
		li.nx*ny*zdim*sizeof(float));
          qmsg += li.nx*ny*zdim;
        }
      }
    }
    SET_PRIORITY(newmsg,grid_sequence,PME_UNTRANS_PRIORITY)
    if ( node == myGridNode ) newmsg->ny = 0;
    if ( npe > 1 ) {
      if ( node == myGridNode ) fwdSharedUntrans(newmsg);
      else pmeNodeProxy[gridNodeInfo[node].real_node].recvUntrans(newmsg);
    } else pmeProxy[gridPeMap[gridNodeInfo[node].pe_start]].recvUntrans(newmsg);
  }

  trans_count = numGridPes;
}

void ComputePmeMgr::fwdSharedUntrans(PmeUntransMsg *msg) {
  int pe = gridNodeInfo[myGridNode].pe_start;
  int npe = gridNodeInfo[myGridNode].npe;
  CmiNodeLock lock = CmiCreateLock();
  int *count = new int; *count = npe;
  for (int i=0; i<npe; ++i, ++pe) {
    PmeSharedUntransMsg *shmsg = new PmeSharedUntransMsg;
    shmsg->msg = msg;
    shmsg->count = count;
    shmsg->lock = lock;
    pmeProxy[gridPeMap[pe]].recvSharedUntrans(shmsg);
  }
}

void ComputePmeMgr::recvSharedUntrans(PmeSharedUntransMsg *msg) {
  procUntrans(msg->msg);
  CmiLock(msg->lock);
  int count = --(*msg->count);
  CmiUnlock(msg->lock);
  if ( count == 0 ) {
    CmiDestroyLock(msg->lock);
    delete msg->count;
    delete msg->msg;
  }
  delete msg;
}

void ComputePmeMgr::recvUntrans(PmeUntransMsg *msg) {
  procUntrans(msg);
  delete msg;
}

void ComputePmeMgr::procUntrans(PmeUntransMsg *msg) {
  // CkPrintf("recvUntrans on Pe(%d)\n",CkMyPe());

#if CMK_BLUEGENEL
  CmiNetworkProgressAfter (0);
#endif

  NodePmeInfo &nodeInfo(gridNodeInfo[myGridNode]);
  int first_pe = nodeInfo.pe_start;
  int g;

 if ( msg->ny ) {
  int zdim = myGrid.dim3;
  int last_pe = first_pe+nodeInfo.npe-1;
  int x_skip = localInfo[myGridPe].x_start
             - localInfo[first_pe].x_start;
  int nx_msg = localInfo[last_pe].x_start
             + localInfo[last_pe].nx
             - localInfo[first_pe].x_start;
  int nx = localInfo[myGridPe].nx;
  int y_start = msg->y_start;
  int ny = msg->ny;
  int slicelen = myGrid.K2 * zdim;
  int cpylen = ny * zdim;
  for ( g=0; g<numGrids; ++g ) {
    float *q = qgrid + qgrid_size * g + y_start * zdim;
    float *qmsg = msg->qgrid + (nx_msg*g+x_skip) * cpylen;
    for ( int x = 0; x < nx; ++x ) {
      CmiMemcpy((void*)q, (void*)qmsg, cpylen*sizeof(float));
      q += slicelen;
      qmsg += cpylen;
    }
  }
 }

  --untrans_count;

  if ( untrans_count == 0 ) {
    pmeProxyDir[CkMyPe()].gridCalc3();
  }
}

void ComputePmeMgr::gridCalc3(void) {
  // CkPrintf("gridCalc3 on Pe(%d)\n",CkMyPe());

  // finish backward FFT
#ifdef NAMD_FFTW

  for ( int g=0; g<numGrids; ++g ) {
#ifdef NAMD_FFTW_3
    fftwf_execute(backward_plan_yz[g]);
#else
    rfftwnd_complex_to_real(backward_plan_yz, localInfo[myGridPe].nx,
	(fftw_complex *) (qgrid + qgrid_size * g),
	1, myGrid.dim2 * myGrid.dim3 / 2, 0, 0, 0);
#endif
  }

#endif

  pmeProxyDir[CkMyPe()].sendUngrid();
}

void ComputePmeMgr::sendUngrid(void) {

#ifdef NAMD_CUDA
  const int UNGRID_PRIORITY = ( offload ? PME_OFFLOAD_UNGRID_PRIORITY : PME_UNGRID_PRIORITY );
#else
  const int UNGRID_PRIORITY = PME_UNGRID_PRIORITY ;
#endif

  for ( int j=0; j<numSources; ++j ) {
    // int msglen = qgrid_len;
    PmeGridMsg *newmsg = gridmsg_reuse[j];
    int pe = newmsg->sourceNode;
    int zdim = myGrid.dim3;
    int flen = newmsg->len;
    int fstart = newmsg->start;
    int zlistlen = newmsg->zlistlen;
    int *zlist = newmsg->zlist;
    float *qmsg = newmsg->qgrid;
    for ( int g=0; g<numGrids; ++g ) {
      char *f = newmsg->fgrid + fgrid_len * g;
      float *q = qgrid + qgrid_size * g + (fstart-fgrid_start) * zdim;
      for ( int i=0; i<flen; ++i ) {
        if ( f[i] ) {
          for ( int k=0; k<zlistlen; ++k ) {
            *(qmsg++) = q[zlist[k]];
          }
        }
        q += zdim;
      }
    }
    newmsg->sourceNode = myGridPe;

    SET_PRIORITY(newmsg,grid_sequence,UNGRID_PRIORITY)
    CmiEnableUrgentSend(1);
#ifdef NAMD_CUDA
    if ( offload ) {
      pmeNodeProxy[CkNodeOf(pe)].recvUngrid(newmsg);
    } else
#endif
    pmeProxyDir[pe].recvUngrid(newmsg);
    CmiEnableUrgentSend(0);
  }
  grid_count = numSources;
  memset( (void*) qgrid, 0, qgrid_size * numGrids * sizeof(float) );
}

void ComputePmeMgr::recvUngrid(PmeGridMsg *msg) {
  // CkPrintf("recvUngrid on Pe(%d)\n",CkMyPe());
#ifdef NAMD_CUDA
  if ( ! offload )  // would need lock
#endif
  if ( ungrid_count == 0 ) {
    NAMD_bug("Message order failure in ComputePmeMgr::recvUngrid\n");
  }

  if ( usePencils ) copyPencils(msg);
  else copyResults(msg);
  delete msg;
  recvAck(0);
}

void ComputePmeMgr::recvAck(PmeAckMsg *msg) {
  if ( msg ) delete msg;
#ifdef NAMD_CUDA
  if ( offload ) {
    CmiLock(cuda_lock);
    if ( ungrid_count == 0 ) {
      NAMD_bug("Message order failure in ComputePmeMgr::recvUngrid\n");
    }
    int uc = --ungrid_count;
    CmiUnlock(cuda_lock);

    if ( uc == 0 ) {
      pmeProxyDir[master_pe].ungridCalc();
    }
    return;
  }
#endif
  --ungrid_count;

  if ( ungrid_count == 0 ) {
    pmeProxyDir[CkMyPe()].ungridCalc();
  }
}

#ifdef NAMD_CUDA
#define count_limit 1000000
#define CUDA_POLL(FN,ARG) CcdCallFnAfter(FN,ARG,0.1)
#define EVENT_STRIDE 10

extern "C" void CcdCallBacksReset(void *ignored,double curWallTime);  // fix Charm++

void cuda_check_pme_forces(void *arg, double walltime) {
  ComputePmeMgr *argp = (ComputePmeMgr *) arg;

 while ( 1 ) { // process multiple events per call
  cudaError_t err = cudaEventQuery(argp->end_forces[argp->forces_done_count/EVENT_STRIDE]);
  if ( err == cudaSuccess ) {
    argp->check_forces_count = 0;
    for ( int i=0; i<EVENT_STRIDE; ++i ) {
      WorkDistrib::messageEnqueueWork(argp->pmeComputes[argp->forces_done_count]);
      if ( ++(argp->forces_done_count) == argp->forces_count ) break;
    }
    if ( argp->forces_done_count == argp->forces_count ) { // last event
      traceUserBracketEvent(CUDA_EVENT_ID_PME_FORCES,argp->forces_time,walltime);
      argp->forces_time = walltime - argp->forces_time;
      //CkPrintf("cuda_check_pme_forces forces_time == %f\n", argp->forces_time);
      return;
    } else { // more events
      continue; // check next event
    }
  } else if ( err != cudaErrorNotReady ) {
    cuda_errcheck("in cuda_check_pme_forces");
    NAMD_bug("cuda_errcheck missed error in cuda_check_pme_forces");
  } else if ( ++(argp->check_forces_count) >= count_limit ) {
    char errmsg[256];
    sprintf(errmsg,"cuda_check_pme_forces polled %d times over %f s on seq %d",
            argp->check_forces_count, walltime - argp->forces_time,
            argp->saved_sequence);
    cuda_errcheck(errmsg);
    NAMD_die(errmsg);
  } else {
    break; // call again
  }
 } // while ( 1 )
 CcdCallBacksReset(0,walltime);  // fix Charm++
 CUDA_POLL(cuda_check_pme_forces, arg);
}
#endif // NAMD_CUDA

void ComputePmeMgr::ungridCalc(void) {
  // CkPrintf("ungridCalc on Pe(%d)\n",CkMyPe());

  ungridForcesCount = pmeComputes.size();

#ifdef NAMD_CUDA
 if ( offload ) {
  //CmiLock(cuda_lock);

  if ( this == masterPmeMgr ) {
    double before = CmiWallTimer();
    cudaMemcpyAsync(v_data_dev, q_data_host, q_data_size, cudaMemcpyHostToDevice, 0 /*streams[stream]*/);
    cudaEventRecord(nodePmeMgr->end_potential_memcpy, 0 /*streams[stream]*/);
    traceUserBracketEvent(CUDA_EVENT_ID_PME_COPY,before,CmiWallTimer());

    const int myrank = CkMyRank();
    for ( int i=0; i<CkMyNodeSize(); ++i ) {
      if ( myrank != i && nodePmeMgr->mgrObjects[i]->pmeComputes.size() ) {
        nodePmeMgr->mgrObjects[i]->ungridCalc();
      }
    }
    if ( ! pmeComputes.size() ) return;
  }

  if ( ! end_forces ) {
    int n=(pmeComputes.size()-1)/EVENT_STRIDE+1;
    end_forces = new cudaEvent_t[n];
    for ( int i=0; i<n; ++i ) {
      cudaEventCreateWithFlags(&end_forces[i],cudaEventDisableTiming);
    }
  }

  const int pcsz = pmeComputes.size();
  if ( ! afn_host ) {
    cudaMallocHost((void**) &afn_host, 3*pcsz*sizeof(float*));
    cudaMalloc((void**) &afn_dev, 3*pcsz*sizeof(float*));
    cuda_errcheck("malloc params for pme");
  }
  int totn = 0;
  for ( int i=0; i<pcsz; ++i ) {
    int n = pmeComputes[i]->numGridAtoms[0];
    totn += n;
  }
  if ( totn > f_data_mgr_alloc ) {
    if ( f_data_mgr_alloc ) {
      CkPrintf("Expanding CUDA forces allocation because %d > %d\n", totn, f_data_mgr_alloc);
      cudaFree(f_data_mgr_dev);
      cudaFreeHost(f_data_mgr_host);
    }
    f_data_mgr_alloc = 1.2 * (totn + 100);
    cudaMalloc((void**) &f_data_mgr_dev, 3*f_data_mgr_alloc*sizeof(float));
    cudaMallocHost((void**) &f_data_mgr_host, 3*f_data_mgr_alloc*sizeof(float));
    cuda_errcheck("malloc forces for pme");
  }
  // CkPrintf("pe %d pcsz %d totn %d alloc %d\n", CkMyPe(), pcsz, totn, f_data_mgr_alloc);
  float *f_dev = f_data_mgr_dev;
  float *f_host = f_data_mgr_host;
  for ( int i=0; i<pcsz; ++i ) {
    int n = pmeComputes[i]->numGridAtoms[0];
    pmeComputes[i]->f_data_dev = f_dev;
    pmeComputes[i]->f_data_host = f_host;
    afn_host[3*i  ] = a_data_dev + 7 * pmeComputes[i]->cuda_atoms_offset;
    afn_host[3*i+1] = f_dev;
    afn_host[3*i+2] = f_dev + n;  // avoid type conversion issues
    f_dev += 3*n;
    f_host += 3*n;
  }
  //CmiLock(cuda_lock);
  double before = CmiWallTimer();
  cudaMemcpyAsync(afn_dev, afn_host, 3*pcsz*sizeof(float*), cudaMemcpyHostToDevice, streams[stream]);
  traceUserBracketEvent(CUDA_EVENT_ID_PME_COPY,before,CmiWallTimer());
  cudaStreamWaitEvent(streams[stream], nodePmeMgr->end_potential_memcpy, 0);
  traceUserEvent(CUDA_EVENT_ID_PME_TICK);

  for ( int i=0; i<pcsz; ++i ) {
    // cudaMemsetAsync(pmeComputes[i]->f_data_dev, 0, 3*n*sizeof(float), streams[stream]);
    if ( i%EVENT_STRIDE == 0 ) {
      int dimy = pcsz - i;
      if ( dimy > EVENT_STRIDE ) dimy = EVENT_STRIDE;
      int maxn = 0;
      int subtotn = 0;
      for ( int j=0; j<dimy; ++j ) {
        int n = pmeComputes[i+j]->numGridAtoms[0];
        subtotn += n;
        if ( n > maxn ) maxn = n;
      }
      // CkPrintf("pe %d dimy %d maxn %d subtotn %d\n", CkMyPe(), dimy, maxn, subtotn);
      before = CmiWallTimer();
      cuda_pme_forces(
        bspline_coeffs_dev,
        v_arr_dev, afn_dev+3*i, dimy, maxn, /*
        pmeComputes[i]->a_data_dev,
        pmeComputes[i]->f_data_dev,
        n, */ myGrid.K1, myGrid.K2, myGrid.K3, myGrid.order,
        streams[stream]);
      traceUserBracketEvent(CUDA_EVENT_ID_PME_KERNEL,before,CmiWallTimer());
      before = CmiWallTimer();
      cudaMemcpyAsync(pmeComputes[i]->f_data_host, pmeComputes[i]->f_data_dev, 3*subtotn*sizeof(float),
        cudaMemcpyDeviceToHost, streams[stream]);
      traceUserBracketEvent(CUDA_EVENT_ID_PME_COPY,before,CmiWallTimer());
      cudaEventRecord(end_forces[i/EVENT_STRIDE], streams[stream]);
      traceUserEvent(CUDA_EVENT_ID_PME_TICK);
    }
    // CkPrintf("pe %d c %d natoms %d fdev %lld fhost %lld\n", CkMyPe(), i, (int64)afn_host[3*i+2], pmeComputes[i]->f_data_dev, pmeComputes[i]->f_data_host);
  }
  //CmiUnlock(cuda_lock);
 } else
#endif // NAMD_CUDA
 {
  for ( int i=0; i<pmeComputes.size(); ++i ) {
    WorkDistrib::messageEnqueueWork(pmeComputes[i]);
    // pmeComputes[i]->ungridForces();
  }
 }
  // submitReductions();  // must follow all ungridForces()

#ifdef NAMD_CUDA
 if ( offload ) {
  forces_time = CmiWallTimer();
  forces_count = ungridForcesCount;
  forces_done_count = 0;
  pmeProxy[this_pe].pollForcesReady();
 }
#endif

  ungrid_count = (usePencils ? numPencilsActive : numDestRecipPes );
}

void ComputePmeMgr::pollForcesReady() {
#ifdef NAMD_CUDA
  CcdCallBacksReset(0,CmiWallTimer());  // fix Charm++
  CUDA_POLL(cuda_check_pme_forces,this);
#else
  NAMD_bug("ComputePmeMgr::pollForcesReady() called in non-CUDA build.");
#endif
}

void ComputePme::atomUpdate() { atomsChanged = 1; }

ComputePme::ComputePme(ComputeID c, PatchID pid) : Compute(c), patchID(pid)
{
  DebugM(4,"ComputePme created.\n");
  basePriority = PME_PRIORITY;
  setNumPatches(1);

  CProxy_ComputePmeMgr::ckLocalBranch(
	CkpvAccess(BOCclass_group).computePmeMgr)->addCompute(this);

  SimParameters *simParams = Node::Object()->simParameters;

  qmForcesOn =  simParams->qmForcesOn;
  offload = simParams->PMEOffload;

  alchOn = simParams->alchOn;
  alchFepOn = simParams->alchFepOn;
  alchThermIntOn = simParams->alchThermIntOn;
  alchDecouple = alchOn && simParams->alchDecouple;
  alchElecLambdaStart = alchOn ? simParams->alchElecLambdaStart : 0;
            
  if (alchOn) {
    numGrids = 2;
    if (alchDecouple) numGrids += 2;
    if (alchElecLambdaStart || alchThermIntOn) numGrids ++;
  }
  else numGrids = 1;
  lesOn = simParams->lesOn;
  if ( lesOn ) {
    lesFactor = simParams->lesFactor;
    numGrids = lesFactor;
  }
  selfOn = 0;
  pairOn = simParams->pairInteractionOn;
  if ( pairOn ) {
    selfOn = simParams->pairInteractionSelf;
    if ( selfOn ) pairOn = 0;  // make pairOn and selfOn exclusive
    numGrids = selfOn ? 1 : 3;
  }

  myGrid.K1 = simParams->PMEGridSizeX;
  myGrid.K2 = simParams->PMEGridSizeY;
  myGrid.K3 = simParams->PMEGridSizeZ;
  myGrid.order = simParams->PMEInterpOrder;
  myGrid.dim2 = myGrid.K2;
  myGrid.dim3 = 2 * (myGrid.K3/2 + 1);

#ifdef NAMD_CUDA
  cuda_atoms_offset = 0;
  f_data_host = 0;
  f_data_dev = 0;
 if ( ! offload )
#endif
 {
  for ( int g=0; g<numGrids; ++g ) myRealSpace[g] = new PmeRealSpace(myGrid);
 }

  atomsChanged = 0;
  
  qmLoclIndx = 0;
  qmLocalCharges = 0;
}

void ComputePme::initialize() {
  if (!(patch = PatchMap::Object()->patch(patchID))) {
    NAMD_bug("ComputePme used with unknown patch.");
  }
  positionBox = patch->registerPositionPickup(this);
  avgPositionBox = patch->registerAvgPositionPickup(this);
  forceBox = patch->registerForceDeposit(this);
#ifdef NAMD_CUDA
 if ( offload ) {
  myMgr->cuda_atoms_count += patch->getNumAtoms();
 }
#endif
}

void ComputePmeMgr::initialize_computes() {

  noWorkCount = 0;
  doWorkCount = 0;
  ungridForcesCount = 0;

  reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);

  SimParameters *simParams = Node::Object()->simParameters;

  strayChargeErrors = 0;

#ifdef NAMD_CUDA
 PatchMap *patchMap = PatchMap::Object();
 int pe = master_pe = CkNodeFirst(CkMyNode());
 for ( int i=0; i<CkMyNodeSize(); ++i, ++pe ) {
    if ( ! patchMap->numPatchesOnNode(master_pe) ) master_pe = pe;
    if ( ! patchMap->numPatchesOnNode(pe) ) continue;
    if ( master_pe < 1 && pe != deviceCUDA->getMasterPe() ) master_pe = pe;
    if ( master_pe == deviceCUDA->getMasterPe() ) master_pe = pe;
    if ( WorkDistrib::pe_sortop_diffuse()(pe,master_pe)
        && pe != deviceCUDA->getMasterPe() ) {
      master_pe = pe;
    }
 }
 if ( ! patchMap->numPatchesOnNode(master_pe) ) {
   NAMD_bug("ComputePmeMgr::initialize_computes() master_pe has no patches.");
 }

 masterPmeMgr = nodePmeMgr->mgrObjects[master_pe - CkNodeFirst(CkMyNode())];
 bool cudaFirst = 1;
 if ( offload ) {
  CmiLock(cuda_lock);
  cudaFirst = ! masterPmeMgr->chargeGridSubmittedCount++;
 }

 if ( cudaFirst ) {
  nodePmeMgr->master_pe = master_pe;
  nodePmeMgr->masterPmeMgr = masterPmeMgr;
 }
#endif

  qsize = myGrid.K1 * myGrid.dim2 * myGrid.dim3;
  fsize = myGrid.K1 * myGrid.dim2;
  if ( myGrid.K2 != myGrid.dim2 ) NAMD_bug("PME myGrid.K2 != myGrid.dim2");
#ifdef NAMD_CUDA
 if ( ! offload )
#endif
 {
  q_arr = new float*[fsize*numGrids];
  memset( (void*) q_arr, 0, fsize*numGrids * sizeof(float*) );
  q_list = new float*[fsize*numGrids];
  memset( (void*) q_list, 0, fsize*numGrids * sizeof(float*) );
  q_count = 0;
 }

#ifdef NAMD_CUDA
 if ( cudaFirst || ! offload ) {
#endif
  f_arr = new char[fsize*numGrids];
  // memset to non-zero value has race condition on BlueGene/Q
  // memset( (void*) f_arr, 2, fsize*numGrids * sizeof(char) );
  for ( int n=fsize*numGrids, i=0; i<n; ++i ) f_arr[i] = 2;

  for ( int g=0; g<numGrids; ++g ) {
    char *f = f_arr + g*fsize;
    if ( usePencils ) {
      int K1 = myGrid.K1;
      int K2 = myGrid.K2;
      int block1 = ( K1 + xBlocks - 1 ) / xBlocks;
      int block2 = ( K2 + yBlocks - 1 ) / yBlocks;
      int dim2 = myGrid.dim2;
      for (int ap=0; ap<numPencilsActive; ++ap) {
        int ib = activePencils[ap].i;
        int jb = activePencils[ap].j;
        int ibegin = ib*block1;
        int iend = ibegin + block1;  if ( iend > K1 ) iend = K1;
        int jbegin = jb*block2;
        int jend = jbegin + block2;  if ( jend > K2 ) jend = K2;
        int flen = numGrids * (iend - ibegin) * (jend - jbegin);
        for ( int i=ibegin; i<iend; ++i ) {
          for ( int j=jbegin; j<jend; ++j ) {
            f[i*dim2+j] = 0;
          }
        }
      }
    } else {
      int block1 = ( myGrid.K1 + numGridPes - 1 ) / numGridPes;
      bsize = block1 * myGrid.dim2 * myGrid.dim3;
      for (int pe=0; pe<numGridPes; pe++) {
        if ( ! recipPeDest[pe] ) continue;
        int start = pe * bsize;
        int len = bsize;
        if ( start >= qsize ) { start = 0; len = 0; }
        if ( start + len > qsize ) { len = qsize - start; }
        int zdim = myGrid.dim3;
        int fstart = start / zdim;
        int flen = len / zdim;
        memset(f + fstart, 0, flen*sizeof(char));
        // CkPrintf("pe %d enabled slabs %d to %d\n", CkMyPe(), fstart/myGrid.dim2, (fstart+flen)/myGrid.dim2-1);
      }
    }
  }
#ifdef NAMD_CUDA
 }
 if ( offload ) {
 if ( cudaFirst ) {

  int f_alloc_count = 0;
  for ( int n=fsize, i=0; i<n; ++i ) {
    if ( f_arr[i] == 0 ) {
      ++f_alloc_count;
    }
  }
  // CkPrintf("pe %d f_alloc_count == %d (%d slabs)\n", CkMyPe(), f_alloc_count, f_alloc_count/myGrid.dim2);

  q_arr = new float*[fsize*numGrids];
  memset( (void*) q_arr, 0, fsize*numGrids * sizeof(float*) );

  float **q_arr_dev_host = new float*[fsize];
  cudaMalloc((void**) &q_arr_dev, fsize * sizeof(float*));

  float **v_arr_dev_host = new float*[fsize];
  cudaMalloc((void**) &v_arr_dev, fsize * sizeof(float*));

  int q_stride = myGrid.K3+myGrid.order-1;
  q_data_size = f_alloc_count * q_stride * sizeof(float);
  ffz_size = (fsize + q_stride) * sizeof(int);

  // tack ffz onto end of q_data to allow merged transfer
  cudaMallocHost((void**) &q_data_host, q_data_size+ffz_size);
  ffz_host = (int*)(((char*)q_data_host) + q_data_size);
  cudaMalloc((void**) &q_data_dev, q_data_size+ffz_size);
  ffz_dev = (int*)(((char*)q_data_dev) + q_data_size);
  cudaMalloc((void**) &v_data_dev, q_data_size);
  cuda_errcheck("malloc grid data for pme");
  cudaMemset(q_data_dev, 0, q_data_size + ffz_size);  // for first time
  cudaEventCreateWithFlags(&(nodePmeMgr->end_charge_memset),cudaEventDisableTiming);
  cudaEventRecord(nodePmeMgr->end_charge_memset, 0);
  cudaEventCreateWithFlags(&(nodePmeMgr->end_all_pme_kernels),cudaEventDisableTiming);
  cudaEventCreateWithFlags(&(nodePmeMgr->end_potential_memcpy),cudaEventDisableTiming);

  f_alloc_count = 0;
  for ( int n=fsize, i=0; i<n; ++i ) {
    if ( f_arr[i] == 0 ) {
      q_arr[i] = q_data_host + f_alloc_count * q_stride;
      q_arr_dev_host[i] = q_data_dev + f_alloc_count * q_stride;
      v_arr_dev_host[i] = v_data_dev + f_alloc_count * q_stride;
      ++f_alloc_count;
    } else {
      q_arr[i] = 0;
      q_arr_dev_host[i] = 0;
      v_arr_dev_host[i] = 0;
    }
  }

  cudaMemcpy(q_arr_dev, q_arr_dev_host, fsize * sizeof(float*), cudaMemcpyHostToDevice);
  cudaMemcpy(v_arr_dev, v_arr_dev_host, fsize * sizeof(float*), cudaMemcpyHostToDevice);
  delete [] q_arr_dev_host;
  delete [] v_arr_dev_host;
  delete [] f_arr;
  f_arr = new char[fsize + q_stride];
  fz_arr = f_arr + fsize;
  memset(f_arr, 0, fsize + q_stride);
  memset(ffz_host, 0, (fsize + q_stride)*sizeof(int));

  cuda_errcheck("initialize grid data for pme");

  cuda_init_bspline_coeffs(&bspline_coeffs_dev, &bspline_dcoeffs_dev, myGrid.order);
  cuda_errcheck("initialize bspline coefficients for pme");

#define XCOPY(X) masterPmeMgr->X = X;
  XCOPY(bspline_coeffs_dev)
  XCOPY(bspline_dcoeffs_dev)
  XCOPY(q_arr)
  XCOPY(q_arr_dev)
  XCOPY(v_arr_dev)
  XCOPY(q_data_size)
  XCOPY(q_data_host)
  XCOPY(q_data_dev)
  XCOPY(v_data_dev)
  XCOPY(ffz_size)
  XCOPY(ffz_host)
  XCOPY(ffz_dev)
  XCOPY(f_arr)
  XCOPY(fz_arr)
#undef XCOPY
  //CkPrintf("pe %d init first\n", CkMyPe());
 } else { // cudaFirst
  //CkPrintf("pe %d init later\n", CkMyPe());
#define XCOPY(X) X = masterPmeMgr->X;
  XCOPY(bspline_coeffs_dev)
  XCOPY(bspline_dcoeffs_dev)
  XCOPY(q_arr)
  XCOPY(q_arr_dev)
  XCOPY(v_arr_dev)
  XCOPY(q_data_size)
  XCOPY(q_data_host)
  XCOPY(q_data_dev)
  XCOPY(v_data_dev)
  XCOPY(ffz_size)
  XCOPY(ffz_host)
  XCOPY(ffz_dev)
  XCOPY(f_arr)
  XCOPY(fz_arr)
#undef XCOPY
 } // cudaFirst
  CmiUnlock(cuda_lock);
 } else // offload
#endif // NAMD_CUDA
 {
  fz_arr = new char[myGrid.K3+myGrid.order-1];
 }

#if 0 && USE_PERSISTENT
  recvGrid_handle = NULL;
#endif
}

ComputePme::~ComputePme()
{
#ifdef NAMD_CUDA
  if ( ! offload )
#endif
  {
    for ( int g=0; g<numGrids; ++g ) delete myRealSpace[g];
  }
}

#if 0 && USE_PERSISTENT 
void ComputePmeMgr::setup_recvgrid_persistent() 
{
    int K1 = myGrid.K1;
    int K2 = myGrid.K2;
    int dim2 = myGrid.dim2;
    int dim3 = myGrid.dim3;
    int block1 = myGrid.block1;
    int block2 = myGrid.block2;

    CkArray *zPencil_local = zPencil.ckLocalBranch();
    recvGrid_handle = (PersistentHandle*) malloc( sizeof(PersistentHandle) * numPencilsActive);
    for (int ap=0; ap<numPencilsActive; ++ap) {
        int ib = activePencils[ap].i;
        int jb = activePencils[ap].j;
        int ibegin = ib*block1;
        int iend = ibegin + block1;  if ( iend > K1 ) iend = K1;
        int jbegin = jb*block2;
        int jend = jbegin + block2;  if ( jend > K2 ) jend = K2;
        int flen = numGrids * (iend - ibegin) * (jend - jbegin);
        // f is changing
        int fcount = 0;
        for ( int g=0; g<numGrids; ++g ) {
            char *f = f_arr + g*fsize;
            for ( int i=ibegin; i<iend; ++i ) {
                for ( int j=jbegin; j<jend; ++j ) {
                    fcount += f[i*dim2+j];
                }
            }
        }
        int zlistlen = 0;
        for ( int i=0; i<myGrid.K3; ++i ) {
            if ( fz_arr[i] ) ++zlistlen;
        }
        int hd = ( fcount? 1 : 0 );  // has data?
        int peer = zPencil_local->homePe(CkArrayIndex3D(ib, jb, 0));
        int compress_start = sizeof(PmeGridMsg ) + sizeof(envelope) + sizeof(int)*hd*zlistlen + sizeof(char)*hd*flen +sizeof(PmeReduction)*hd*numGrids ;
        int compress_size = sizeof(float)*hd*fcount*zlistlen;
        int size = compress_start +  compress_size  + PRIORITY_SIZE/8+6;
        recvGrid_handle[ap] =  CmiCreateCompressPersistentSize(peer, size, compress_start, compress_size, CMI_FLOATING);
    }
}
#endif

int ComputePme::noWork() {

  if ( patch->flags.doFullElectrostatics ) {
    // In QM/MM simulations, atom charges form QM regions need special treatment.
    if ( qmForcesOn ) {
        return 1;
    }
    if ( ! myMgr->ungridForcesCount && ! myMgr->recipEvirCount ) return 0;  // work to do, enqueue as usual
    myMgr->heldComputes.add(this);
    return 1;  // don't enqueue yet
  }

  positionBox->skip();
  forceBox->skip();

  if ( ++(myMgr->noWorkCount) == myMgr->pmeComputes.size() ) {
    myMgr->noWorkCount = 0;
    myMgr->reduction->submit();
  }

  atomsChanged = 0;

  return 1;  // no work for this step
}

void ComputePmeMgr::addRecipEvirClient() {
  ++recipEvirClients;
}

void ComputePmeMgr::recvRecipEvir(PmeEvirMsg *msg) {
  if ( ! pmeComputes.size() ) NAMD_bug("ComputePmeMgr::recvRecipEvir() called on pe without patches");
  for ( int g=0; g<numGrids; ++g ) {
    evir[g] += msg->evir[g];
  }
  delete msg;
  // CkPrintf("recvRecipEvir pe %d %d %d\n", CkMyPe(), ungridForcesCount, recipEvirCount);
  if ( ! --recipEvirCount && ! ungridForcesCount ) submitReductions();
}

void ComputePme::doQMWork() {
    
//     iout << CkMyPe() << ") ----> PME doQMWork.\n" << endi ;
    
    
    int numQMAtms = Node::Object()->molecule->get_numQMAtoms();
    const Real *qmAtmChrg = Node::Object()->molecule->get_qmAtmChrg() ;
    const int *qmAtmIndx = Node::Object()->molecule->get_qmAtmIndx() ;
    const Real *qmAtomGroup = Node::Object()->molecule->get_qmAtomGroup() ;
    
    const CompAtomExt *xExt = patch->getCompAtomExtInfo();
    
    // Determine number of qm atoms in this patch for the current step.
    numLocalQMAtoms = 0;
    for (int paIter=0; paIter<patch->getNumAtoms(); paIter++) {
        if ( qmAtomGroup[xExt[paIter].id] != 0 ) {
            numLocalQMAtoms++;
        }
    }
    
    // We prepare a charge vector with QM charges for use in the PME calculation.
    
    // Clears data from last step, if there is any.
    if (qmLoclIndx != 0)
        delete [] qmLoclIndx;
    if (qmLocalCharges != 0)
        delete [] qmLocalCharges;
    
    qmLoclIndx = new int[numLocalQMAtoms] ;
    qmLocalCharges = new Real[numLocalQMAtoms] ;
    
    // I am assuming there will be (in general) more QM atoms among all QM groups
    // than MM atoms in a patch.
    int procAtms = 0;
    
    for (int paIter=0; paIter<patch->getNumAtoms(); paIter++) {
        
        for (int i=0; i<numQMAtms; i++) {
            
            if (qmAtmIndx[i] == xExt[paIter].id) {
                
                qmLoclIndx[procAtms] = paIter ;
                qmLocalCharges[procAtms] = qmAtmChrg[i];
                
                procAtms++;
                break;
            }
            
        }
        
        if (procAtms == numLocalQMAtoms)
            break;
    }
    
    doWork();
    return ;
}

void ComputePme::doWork()
{
  DebugM(4,"Entering ComputePme::doWork().\n");

  if ( basePriority >= COMPUTE_HOME_PRIORITY ) {
#ifdef NAMD_CUDA
    basePriority = ( offload ? PME_OFFLOAD_PRIORITY : PME_PRIORITY );
#else
    basePriority = PME_PRIORITY;
#endif
    ungridForces();
    // CkPrintf("doWork 2 pe %d %d %d\n", CkMyPe(), myMgr->ungridForcesCount, myMgr->recipEvirCount);
    if ( ! --(myMgr->ungridForcesCount) && ! myMgr->recipEvirCount ) myMgr->submitReductions();
    return;
  }
  basePriority = COMPUTE_HOME_PRIORITY + PATCH_PRIORITY(patchID);
  // CkPrintf("doWork 1 pe %d %d %d\n", CkMyPe(), myMgr->ungridForcesCount, myMgr->recipEvirCount);

#ifdef TRACE_COMPUTE_OBJECTS
    double traceObjStartTime = CmiWallTimer();
#endif

  // allocate storage
  numLocalAtoms = patch->getNumAtoms();

  Lattice &lattice = patch->flags.lattice;

  localData_alloc.resize(numLocalAtoms*(numGrids+ ((numGrids>1 || selfOn)?1:0)));
  localData = localData_alloc.begin();
  localPartition_alloc.resize(numLocalAtoms);
  localPartition = localPartition_alloc.begin();

  int g;
  for ( g=0; g<numGrids; ++g ) {
    localGridData[g] = localData + numLocalAtoms*(g+1);
  }

  // get positions and charges
  PmeParticle * data_ptr = localData;
  unsigned char * part_ptr = localPartition;
  const BigReal coulomb_sqrt = sqrt( COULOMB * ComputeNonbondedUtil::scaling
				* ComputeNonbondedUtil::dielectric_1 );

  {
    CompAtom *x = positionBox->open();
    // CompAtomExt *xExt = patch->getCompAtomExtInfo();
    if ( patch->flags.doMolly ) {
      positionBox->close(&x);
      x = avgPositionBox->open();
    }
    int numAtoms = patch->getNumAtoms();

    for(int i=0; i<numAtoms; ++i)
    {
      data_ptr->x = x[i].position.x;
      data_ptr->y = x[i].position.y;
      data_ptr->z = x[i].position.z;
      data_ptr->cg = coulomb_sqrt * x[i].charge;
      ++data_ptr;
      *part_ptr = x[i].partition;
      ++part_ptr;
    }

    // QM loop to overwrite charges of QM atoms.
    // They are zero for NAMD, but are updated in ComputeQM.
    if ( qmForcesOn ) {
        
        for(int i=0; i<numLocalQMAtoms; ++i)
        {
          localData[qmLoclIndx[i]].cg = coulomb_sqrt * qmLocalCharges[i];
        }
        
    }
    
    if ( patch->flags.doMolly ) { avgPositionBox->close(&x); }
    else { positionBox->close(&x); }
  }

  // copy to other grids if needed
  if ( (alchOn && (!alchDecouple)) || lesOn ) {
    for ( g=0; g<numGrids; ++g ) {
      PmeParticle *lgd = localGridData[g];
      int nga = 0;
      for(int i=0; i<numLocalAtoms; ++i) {
        if ( localPartition[i] == 0 || localPartition[i] == (g+1) ) {
          // for FEP/TI: grid 0 gets non-alch + partition 1;
          // grid 1 gets non-alch + partition 2;
          // grid 2 (only if called for with numGrids=3) gets only non-alch
          lgd[nga++] = localData[i];
        }
      }
      numGridAtoms[g] = nga;
    }
  } else if ( alchOn && alchDecouple) {
    // alchemical decoupling: four grids
    // g=0: partition 0 and partition 1
    // g=1: partition 0 and partition 2
    // g=2: only partition 1 atoms
    // g=3: only partition 2 atoms
    // plus one grid g=4, only partition 0, if numGrids=5
    for ( g=0; g<2; ++g ) {  // same as before for first 2
      PmeParticle *lgd = localGridData[g];
      int nga = 0;
      for(int i=0; i<numLocalAtoms; ++i) {
        if ( localPartition[i] == 0 || localPartition[i] == (g+1) ) {
          lgd[nga++] = localData[i];
        }
      }
      numGridAtoms[g] = nga;
    }
    for (g=2 ; g<4 ; ++g ) {  // only alchemical atoms for these 2
      PmeParticle *lgd = localGridData[g];
      int nga = 0;
      for(int i=0; i<numLocalAtoms; ++i) {
        if ( localPartition[i] == (g-1) ) {
          lgd[nga++] = localData[i];
        }
      }
      numGridAtoms[g] = nga;
    }
    for (g=4 ; g<numGrids ; ++g ) {  // only non-alchemical atoms 
      // numGrids=5 only if alchElecLambdaStart > 0
      PmeParticle *lgd = localGridData[g];
      int nga = 0;
      for(int i=0; i<numLocalAtoms; ++i) {
        if ( localPartition[i] == 0 ) {
          lgd[nga++] = localData[i];
        }
      }
      numGridAtoms[g] = nga;
    }
  } else if ( selfOn ) {
    if ( numGrids != 1 ) NAMD_bug("ComputePme::doWork assertion 1 failed");
    g = 0;
    PmeParticle *lgd = localGridData[g];
    int nga = 0;
    for(int i=0; i<numLocalAtoms; ++i) {
      if ( localPartition[i] == 1 ) {
        lgd[nga++] = localData[i];
      }
    }
    numGridAtoms[g] = nga;
  } else if ( pairOn ) {
    if ( numGrids != 3 ) NAMD_bug("ComputePme::doWork assertion 2 failed");
    g = 0;
    PmeParticle *lgd = localGridData[g];
    int nga = 0;
    for(int i=0; i<numLocalAtoms; ++i) {
      if ( localPartition[i] == 1 || localPartition[i] == 2 ) {
        lgd[nga++] = localData[i];
      }
    }
    numGridAtoms[g] = nga;
    for ( g=1; g<3; ++g ) {
      PmeParticle *lgd = localGridData[g];
      int nga = 0;
      for(int i=0; i<numLocalAtoms; ++i) {
        if ( localPartition[i] == g ) {
          lgd[nga++] = localData[i];
        }
      }
      numGridAtoms[g] = nga;
    }
  } else {
    if ( numGrids != 1 ) NAMD_bug("ComputePme::doWork assertion 3 failed");
    localGridData[0] = localData;
    numGridAtoms[0] = numLocalAtoms;
  }

 if ( ! myMgr->doWorkCount ) {
  myMgr->doWorkCount = myMgr->pmeComputes.size();

#ifdef NAMD_CUDA
 if ( !  offload )
#endif // NAMD_CUDA
 {
  memset( (void*) myMgr->fz_arr, 0, (myGrid.K3+myGrid.order-1) * sizeof(char) );

  for (int i=0; i<myMgr->q_count; ++i) {
    memset( (void*) (myMgr->q_list[i]), 0, (myGrid.K3+myGrid.order-1) * sizeof(float) );
  }
 }

  for ( g=0; g<numGrids; ++g ) {
    myMgr->evir[g] = 0;
  }

  myMgr->strayChargeErrors = 0;

  myMgr->compute_sequence = sequence();
 }

  if ( sequence() != myMgr->compute_sequence ) NAMD_bug("ComputePme sequence mismatch in doWork()");

  int strayChargeErrors = 0;

  // calculate self energy
  BigReal ewaldcof = ComputeNonbondedUtil::ewaldcof;
  for ( g=0; g<numGrids; ++g ) {
    BigReal selfEnergy = 0;
    data_ptr = localGridData[g];
    int i;
    for(i=0; i<numGridAtoms[g]; ++i)
    {
      selfEnergy += data_ptr->cg * data_ptr->cg;
      ++data_ptr;
    }
    selfEnergy *= -1. * ewaldcof / SQRT_PI;
    myMgr->evir[g][0] += selfEnergy;

    float **q = myMgr->q_arr + g*myMgr->fsize;
    char *f = myMgr->f_arr + g*myMgr->fsize;

    scale_coordinates(localGridData[g], numGridAtoms[g], lattice, myGrid);
#ifdef NAMD_CUDA
   if ( offload ) {
    if ( myMgr->cuda_atoms_alloc == 0 ) {  // first call
      int na = myMgr->cuda_atoms_alloc = 1.2 * (myMgr->cuda_atoms_count + 1000);
      cuda_errcheck("before malloc atom data for pme");
      cudaMallocHost((void**) &(myMgr->a_data_host), 7*na*sizeof(float));
      cudaMalloc((void**) &(myMgr->a_data_dev), 7*na*sizeof(float));
      cuda_errcheck("malloc atom data for pme");
      myMgr->cuda_atoms_count = 0;
    }
    cuda_atoms_offset = myMgr->cuda_atoms_count;
    int n = numGridAtoms[g];
    myMgr->cuda_atoms_count += n;
    if ( myMgr->cuda_atoms_count > myMgr->cuda_atoms_alloc ) {
      CkPrintf("Pe %d expanding CUDA PME atoms allocation because %d > %d\n",
			CkMyPe(), myMgr->cuda_atoms_count, myMgr->cuda_atoms_alloc);
      cuda_errcheck("before malloc expanded atom data for pme");
      int na = myMgr->cuda_atoms_alloc = 1.2 * (myMgr->cuda_atoms_count + 1000);
      const float *a_data_host_old = myMgr->a_data_host;
      cudaMallocHost((void**) &(myMgr->a_data_host), 7*na*sizeof(float));
      cuda_errcheck("malloc expanded host atom data for pme");
      memcpy(myMgr->a_data_host, a_data_host_old, 7*cuda_atoms_offset*sizeof(float));
      cudaFreeHost((void*) a_data_host_old);
      cuda_errcheck("free expanded host atom data for pme");
      cudaFree(myMgr->a_data_dev);
      cuda_errcheck("free expanded dev atom data for pme");
      cudaMalloc((void**) &(myMgr->a_data_dev), 7*na*sizeof(float));
      cuda_errcheck("malloc expanded dev atom data for pme");
    }
    float *a_data_host = myMgr->a_data_host + 7 * cuda_atoms_offset;
    data_ptr = localGridData[g];
    double order_1 = myGrid.order - 1;
    double K1 = myGrid.K1;
    double K2 = myGrid.K2;
    double K3 = myGrid.K3;
    int found_negative = 0;
    for ( int i=0; i<n; ++i ) {
      if ( data_ptr[i].x < 0 || data_ptr[i].y < 0 || data_ptr[i].z < 0 ) {
        found_negative = 1;
        // CkPrintf("low coord: %f %f %f\n", data_ptr[i].x, data_ptr[i].y, data_ptr[i].z);
      }
      double x_int = (int) data_ptr[i].x;
      double y_int = (int) data_ptr[i].y;
      double z_int = (int) data_ptr[i].z;
      a_data_host[7*i  ] = data_ptr[i].x - x_int;  // subtract in double precision
      a_data_host[7*i+1] = data_ptr[i].y - y_int;
      a_data_host[7*i+2] = data_ptr[i].z - z_int;
      a_data_host[7*i+3] = data_ptr[i].cg;
      x_int -= order_1;  if ( x_int < 0 ) x_int += K1;
      y_int -= order_1;  if ( y_int < 0 ) y_int += K2;
      z_int -= order_1;  if ( z_int < 0 ) z_int += K3;
      a_data_host[7*i+4] = x_int;
      a_data_host[7*i+5] = y_int;
      a_data_host[7*i+6] = z_int;
    }
    if ( found_negative ) NAMD_bug("found negative atom coordinate in ComputePme::doWork");
   } else
#endif // NAMD_CUDA
   {
    myRealSpace[g]->set_num_atoms(numGridAtoms[g]);
    myRealSpace[g]->fill_charges(q, myMgr->q_list, myMgr->q_count, strayChargeErrors, f, myMgr->fz_arr, localGridData[g]);
   }
  }
  myMgr->strayChargeErrors += strayChargeErrors;

#ifdef TRACE_COMPUTE_OBJECTS
    traceUserBracketEvent(TRACE_COMPOBJ_IDOFFSET+this->cid, traceObjStartTime, CmiWallTimer());
#endif

 if ( --(myMgr->doWorkCount) == 0 ) {
// cudaDeviceSynchronize();  // XXXX
#ifdef NAMD_CUDA
  if ( offload ) {
    ComputePmeMgr::cuda_submit_charges_args args;
    args.mgr = myMgr;
    args.lattice = &lattice;
    args.sequence = sequence();
    CmiLock(ComputePmeMgr::cuda_lock);
    if ( ComputePmeMgr::cuda_busy ) {
      ComputePmeMgr::cuda_submit_charges_deque.push_back(args);
    } else if ( CkMyPe() == deviceCUDA->getMasterPe() ) {
      // avoid adding work to nonbonded data preparation pe
      args.mgr->cuda_submit_charges(*args.lattice, args.sequence);
    } else {
      ComputePmeMgr::cuda_busy = true;
      while ( 1 ) {
        CmiUnlock(ComputePmeMgr::cuda_lock);
        args.mgr->cuda_submit_charges(*args.lattice, args.sequence);
        CmiLock(ComputePmeMgr::cuda_lock);
        if ( ComputePmeMgr::cuda_submit_charges_deque.size() ) {
          args = ComputePmeMgr::cuda_submit_charges_deque.front();
          ComputePmeMgr::cuda_submit_charges_deque.pop_front();
        } else {
          ComputePmeMgr::cuda_busy = false;
          break;
        }
      }
    }
    CmiUnlock(ComputePmeMgr::cuda_lock);
  } else
#endif // NAMD_CUDA
  {
    myMgr->chargeGridReady(lattice,sequence());
  }
 }
 atomsChanged = 0;
}

#ifdef NAMD_CUDA

void ComputePmeMgr::cuda_submit_charges(Lattice &lattice, int sequence) {

    int n = cuda_atoms_count;
    //CkPrintf("pe %d cuda_atoms_count %d\n", CkMyPe(), cuda_atoms_count);
    cuda_atoms_count = 0;

    const double before = CmiWallTimer();
    cudaMemcpyAsync(a_data_dev, a_data_host, 7*n*sizeof(float),
                          cudaMemcpyHostToDevice, streams[stream]);
    const double after = CmiWallTimer();

    cudaStreamWaitEvent(streams[stream], nodePmeMgr->end_charge_memset, 0);

    cuda_pme_charges(
      bspline_coeffs_dev,
      q_arr_dev, ffz_dev, ffz_dev + fsize,
      a_data_dev, n,
      myGrid.K1, myGrid.K2, myGrid.K3, myGrid.order,
      streams[stream]);
    const double after2 = CmiWallTimer();

    chargeGridSubmitted(lattice,sequence);  // must be inside lock

    masterPmeMgr->charges_time = before;
    traceUserBracketEvent(CUDA_EVENT_ID_PME_COPY,before,after);
    traceUserBracketEvent(CUDA_EVENT_ID_PME_KERNEL,after,after2);
}

void cuda_check_pme_charges(void *arg, double walltime) {
  ComputePmeMgr *argp = (ComputePmeMgr *) arg;

  cudaError_t err = cudaEventQuery(argp->end_charges);
  if ( err == cudaSuccess ) {
    traceUserBracketEvent(CUDA_EVENT_ID_PME_CHARGES,argp->charges_time,walltime);
    argp->charges_time = walltime - argp->charges_time;
    argp->sendChargeGridReady();
    argp->check_charges_count = 0;
  } else if ( err != cudaErrorNotReady ) {
    cuda_errcheck("in cuda_check_pme_charges");
    NAMD_bug("cuda_errcheck missed error in cuda_check_pme_charges");
  } else if ( ++(argp->check_charges_count) >= count_limit ) {
    char errmsg[256];
    sprintf(errmsg,"cuda_check_pme_charges polled %d times over %f s on seq %d",
            argp->check_charges_count, walltime - argp->charges_time,
            argp->saved_sequence);
    cuda_errcheck(errmsg);
    NAMD_die(errmsg);
  } else {
    CcdCallBacksReset(0,walltime);  // fix Charm++
    CUDA_POLL(cuda_check_pme_charges, arg);
  }
}

void ComputePmeMgr::chargeGridSubmitted(Lattice &lattice, int sequence) {
  saved_lattice = &lattice;
  saved_sequence = sequence;

  // cudaDeviceSynchronize();  //  XXXX TESTING
  //int q_stride = myGrid.K3+myGrid.order-1;
  //for (int n=fsize+q_stride, j=0; j<n; ++j) {
  //  if ( ffz_host[j] != 0 && ffz_host[j] != 1 ) {
  //    CkPrintf("pre-memcpy flag %d/%d == %d on pe %d in ComputePmeMgr::chargeGridReady\n", j, n, ffz_host[j], CkMyPe());
  //  }
  //}
  //CmiLock(cuda_lock);

 if ( --(masterPmeMgr->chargeGridSubmittedCount) == 0 ) {
  double before = CmiWallTimer();
  cudaEventRecord(nodePmeMgr->end_all_pme_kernels, 0);  // when all streams complete
  cudaStreamWaitEvent(streams[stream], nodePmeMgr->end_all_pme_kernels, 0);
  cudaMemcpyAsync(q_data_host, q_data_dev, q_data_size+ffz_size,
                        cudaMemcpyDeviceToHost, streams[stream]);
  traceUserBracketEvent(CUDA_EVENT_ID_PME_COPY,before,CmiWallTimer());
  cudaEventRecord(masterPmeMgr->end_charges, streams[stream]);
  cudaMemsetAsync(q_data_dev, 0, q_data_size + ffz_size, streams[stream]);  // for next time
  cudaEventRecord(nodePmeMgr->end_charge_memset, streams[stream]);
  //CmiUnlock(cuda_lock);
  // cudaDeviceSynchronize();  //  XXXX TESTING
  // cuda_errcheck("after memcpy grid to host");

  SimParameters *simParams = Node::Object()->simParameters;
  if ( ! simParams->useCUDA2 ) {
    CProxy_ComputeMgr cm(CkpvAccess(BOCclass_group).computeMgr);
    cm[deviceCUDA->getMasterPe()].recvYieldDevice(-1);
  }

  pmeProxy[master_pe].pollChargeGridReady();
 }
}

void ComputePmeMgr::sendChargeGridReady() {
  for ( int i=0; i<CkMyNodeSize(); ++i ) {
    ComputePmeMgr *mgr = nodePmeMgr->mgrObjects[i];
    int cs = mgr->pmeComputes.size();
    if ( cs ) {
      mgr->ungridForcesCount = cs;
      mgr->recipEvirCount = mgr->recipEvirClients;
      masterPmeMgr->chargeGridSubmittedCount++;
    }
  }
  pmeProxy[master_pe].recvChargeGridReady();
}
#endif // NAMD_CUDA

void ComputePmeMgr::pollChargeGridReady() {
#ifdef NAMD_CUDA
  CcdCallBacksReset(0,CmiWallTimer());  // fix Charm++
  CUDA_POLL(cuda_check_pme_charges,this);
#else
  NAMD_bug("ComputePmeMgr::pollChargeGridReady() called in non-CUDA build.");
#endif
}

void ComputePmeMgr::recvChargeGridReady() {
  chargeGridReady(*saved_lattice,saved_sequence);
}

void ComputePmeMgr::chargeGridReady(Lattice &lattice, int sequence) {

#ifdef NAMD_CUDA
 if ( offload ) {
  int errcount = 0;
  int q_stride = myGrid.K3+myGrid.order-1;
  for (int n=fsize+q_stride, j=fsize; j<n; ++j) {
    f_arr[j] = ffz_host[j];
    if ( ffz_host[j] & ~1 ) ++errcount;
  }
  if ( errcount ) NAMD_bug("bad flag in ComputePmeMgr::chargeGridReady");
 }
#endif
  recipEvirCount = recipEvirClients;
  ungridForcesCount = pmeComputes.size();

  for (int j=0; j<myGrid.order-1; ++j) {
    fz_arr[j] |= fz_arr[myGrid.K3+j];
  }

  if ( usePencils ) {
    sendPencils(lattice,sequence);
  } else {
    sendData(lattice,sequence);
  }
}


void ComputePmeMgr::sendPencilsPart(int first, int last, Lattice &lattice, int sequence, int sourcepe) {

  // iout << "Sending charge grid for " << numLocalAtoms << " atoms to FFT on " << iPE << ".\n" << endi;

#if 0 && USE_PERSISTENT
    if (recvGrid_handle== NULL) setup_recvgrid_persistent();
#endif
  int K1 = myGrid.K1;
  int K2 = myGrid.K2;
  int dim2 = myGrid.dim2;
  int dim3 = myGrid.dim3;
  int block1 = myGrid.block1;
  int block2 = myGrid.block2;

  // int savedMessages = 0;
  NodePmeMgr *npMgr = pmeNodeProxy[CkMyNode()].ckLocalBranch();

  for (int ap=first; ap<=last; ++ap) {
    int ib = activePencils[ap].i;
    int jb = activePencils[ap].j;
    int ibegin = ib*block1;
    int iend = ibegin + block1;  if ( iend > K1 ) iend = K1;
    int jbegin = jb*block2;
    int jend = jbegin + block2;  if ( jend > K2 ) jend = K2;
    int flen = numGrids * (iend - ibegin) * (jend - jbegin);

    int fcount = 0;
    for ( int g=0; g<numGrids; ++g ) {
      char *f = f_arr + g*fsize;
#ifdef NAMD_CUDA
     if ( offload ) {
      int errcount = 0;
      for ( int i=ibegin; i<iend; ++i ) {
       for ( int j=jbegin; j<jend; ++j ) {
        int k = i*dim2+j;
        f[k] = ffz_host[k];
        fcount += f[k];
        if ( ffz_host[k] & ~1 ) ++errcount;
       }
      }
      if ( errcount ) NAMD_bug("bad flag in ComputePmeMgr::sendPencilsPart");
     } else
#endif
      for ( int i=ibegin; i<iend; ++i ) {
       for ( int j=jbegin; j<jend; ++j ) {
        fcount += f[i*dim2+j];
       }
      }
    }

#ifdef NETWORK_PROGRESS
    CmiNetworkProgress();
#endif

    if ( ! pencilActive[ib*yBlocks+jb] )
      NAMD_bug("PME activePencils list inconsistent");

    int zlistlen = 0;
    for ( int i=0; i<myGrid.K3; ++i ) {
      if ( fz_arr[i] ) ++zlistlen;
    }

    int hd = ( fcount? 1 : 0 );  // has data?
    // if ( ! hd ) ++savedMessages;

    
    PmeGridMsg *msg = new ( hd*zlistlen, hd*flen,
	hd*fcount*zlistlen, PRIORITY_SIZE) PmeGridMsg;
    msg->sourceNode = sourcepe;
    msg->hasData = hd;
    msg->lattice = lattice;
   if ( hd ) {
#if 0
    msg->start = fstart;
    msg->len = flen;
#else
    msg->start = -1;   // obsolete?
    msg->len = -1;   // obsolete?
#endif
    msg->zlistlen = zlistlen;
    int *zlist = msg->zlist;
    zlistlen = 0;
    for ( int i=0; i<myGrid.K3; ++i ) {
      if ( fz_arr[i] ) zlist[zlistlen++] = i;
    }
    char *fmsg = msg->fgrid;
    float *qmsg = msg->qgrid;
    for ( int g=0; g<numGrids; ++g ) {
      char *f = f_arr + g*fsize;
      float **q = q_arr + g*fsize;
      for ( int i=ibegin; i<iend; ++i ) {
       for ( int j=jbegin; j<jend; ++j ) {
        *(fmsg++) = f[i*dim2+j];
        if( f[i*dim2+j] ) {
          for (int h=0; h<myGrid.order-1; ++h) {
            q[i*dim2+j][h] += q[i*dim2+j][myGrid.K3+h];
          }
          for ( int k=0; k<zlistlen; ++k ) {
            *(qmsg++) = q[i*dim2+j][zlist[k]];
          }
        }
       }
      }
    }
   }

    msg->sequence = compute_sequence;
    SET_PRIORITY(msg,compute_sequence,PME_GRID_PRIORITY)
    CmiEnableUrgentSend(1);
#if USE_NODE_PAR_RECEIVE
    msg->destElem=CkArrayIndex3D(ib,jb,0);
    CProxy_PmePencilMap lzm = npMgr->zm;
    int destproc = lzm.ckLocalBranch()->procNum(0, msg->destElem);
    int destnode = CmiNodeOf(destproc);
    
#if  0 
    CmiUsePersistentHandle(&recvGrid_handle[ap], 1);
#endif
    pmeNodeProxy[destnode].recvZGrid(msg);
#if 0 
    CmiUsePersistentHandle(NULL, 0);
#endif
#else
#if 0 
    CmiUsePersistentHandle(&recvGrid_handle[ap], 1);
#endif
    zPencil(ib,jb,0).recvGrid(msg);
#if 0 
    CmiUsePersistentHandle(NULL, 0);
#endif
#endif
    CmiEnableUrgentSend(0);
  }


  // if ( savedMessages ) {
  //   CkPrintf("Pe %d eliminated %d PME messages\n",CkMyPe(),savedMessages);
  // }

}


void ComputePmeMgr::sendPencilsHelper(int iter) {
  nodePmeMgr->sendPencilsHelper(iter);
}

void NodePmeMgr::sendPencilsHelper(int iter) {
#ifdef NAMD_CUDA
  ComputePmeMgr *obj = masterPmeMgr;
  obj->sendPencilsPart(iter, iter, *obj->sendDataHelper_lattice, obj->sendDataHelper_sequence, obj->sendDataHelper_sourcepe);
#else
  NAMD_bug("NodePmeMgr::sendPencilsHelper called in non-CUDA build");
#endif
}

void ComputePmeMgr::sendPencils(Lattice &lattice, int sequence) {

  sendDataHelper_lattice = &lattice;
  sendDataHelper_sequence = sequence;
  sendDataHelper_sourcepe = CkMyPe();

#ifdef NAMD_CUDA
  if ( offload ) {
    for ( int ap=0; ap < numPencilsActive; ++ap ) {
#if CMK_MULTICORE
      // nodegroup messages on multicore are delivered to sending pe, or pe 0 if expedited
      int ib = activePencils[ap].i;
      int jb = activePencils[ap].j;
      int destproc = nodePmeMgr->zm.ckLocalBranch()->procNum(0, CkArrayIndex3D(ib,jb,0));
      pmeProxy[destproc].sendPencilsHelper(ap);
#else
      pmeNodeProxy[CkMyNode()].sendPencilsHelper(ap);
#endif
    }
  } else
#endif
  {
    sendPencilsPart(0,numPencilsActive-1,lattice,sequence,CkMyPe());
  }

  if ( strayChargeErrors ) {
   strayChargeErrors = 0;
   iout << iERROR << "Stray PME grid charges detected: "
 	<< CkMyPe() << " sending to (x,y)";
   int K1 = myGrid.K1;
   int K2 = myGrid.K2;
   int dim2 = myGrid.dim2;
   int block1 = myGrid.block1;
   int block2 = myGrid.block2;
   for (int ib=0; ib<xBlocks; ++ib) {
    for (int jb=0; jb<yBlocks; ++jb) {
     int ibegin = ib*block1;
     int iend = ibegin + block1;  if ( iend > K1 ) iend = K1;
     int jbegin = jb*block2;
     int jend = jbegin + block2;  if ( jend > K2 ) jend = K2;
     int flen = numGrids * (iend - ibegin) * (jend - jbegin);

     for ( int g=0; g<numGrids; ++g ) {
       char *f = f_arr + g*fsize;
       if ( ! pencilActive[ib*yBlocks+jb] ) {
           for ( int i=ibegin; i<iend; ++i ) {
            for ( int j=jbegin; j<jend; ++j ) {
             if ( f[i*dim2+j] == 3 ) {
               f[i*dim2+j] = 2;
               iout << " (" << i << "," << j << ")";
             }
            }
           }
       }
     }
    }
   }
   iout << "\n" << endi;
  }
 
}


void ComputePmeMgr::copyPencils(PmeGridMsg *msg) {

  int K1 = myGrid.K1;
  int K2 = myGrid.K2;
  int dim2 = myGrid.dim2;
  int dim3 = myGrid.dim3;
  int block1 = myGrid.block1;
  int block2 = myGrid.block2;

  // msg->sourceNode = thisIndex.x * initdata.yBlocks + thisIndex.y;
  int ib = msg->sourceNode / yBlocks;
  int jb = msg->sourceNode % yBlocks;

  int ibegin = ib*block1;
  int iend = ibegin + block1;  if ( iend > K1 ) iend = K1;
  int jbegin = jb*block2;
  int jend = jbegin + block2;  if ( jend > K2 ) jend = K2;

  int zlistlen = msg->zlistlen;
  int *zlist = msg->zlist;
  float *qmsg = msg->qgrid;
  int g;
  for ( g=0; g<numGrids; ++g ) {
    char *f = f_arr + g*fsize;
    float **q = q_arr + g*fsize;
    for ( int i=ibegin; i<iend; ++i ) {
     for ( int j=jbegin; j<jend; ++j ) {
      if( f[i*dim2+j] ) {
        f[i*dim2+j] = 0;
        for ( int k=0; k<zlistlen; ++k ) {
          q[i*dim2+j][zlist[k]] = *(qmsg++);
        }
        for (int h=0; h<myGrid.order-1; ++h) {
          q[i*dim2+j][myGrid.K3+h] = q[i*dim2+j][h];
        }
      }
     }
    }
  }
}


void ComputePmeMgr::sendDataPart(int first, int last, Lattice &lattice, int sequence, int sourcepe, int errors) {

  // iout << "Sending charge grid for " << numLocalAtoms << " atoms to FFT on " << iPE << ".\n" << endi;

  bsize = myGrid.block1 * myGrid.dim2 * myGrid.dim3;

  CProxy_ComputePmeMgr pmeProxy(CkpvAccess(BOCclass_group).computePmeMgr);
  for (int j=first; j<=last; j++) {
    int pe = gridPeOrder[j];  // different order
    if ( ! recipPeDest[pe] && ! errors ) continue;
    int start = pe * bsize;
    int len = bsize;
    if ( start >= qsize ) { start = 0; len = 0; }
    if ( start + len > qsize ) { len = qsize - start; }
    int zdim = myGrid.dim3;
    int fstart = start / zdim;
    int flen = len / zdim;
    int fcount = 0;
    int i;

    int g;
    for ( g=0; g<numGrids; ++g ) {
      char *f = f_arr + fstart + g*fsize;
#ifdef NAMD_CUDA
     if ( offload ) {
      int errcount = 0;
      for ( i=0; i<flen; ++i ) {
        f[i] = ffz_host[fstart+i];
        fcount += f[i];
        if ( ffz_host[fstart+i] & ~1 ) ++errcount;
      }
      if ( errcount ) NAMD_bug("bad flag in ComputePmeMgr::sendDataPart");
     } else
#endif
      for ( i=0; i<flen; ++i ) {
        fcount += f[i];
      }
      if ( ! recipPeDest[pe] ) {
        int errfound = 0;
        for ( i=0; i<flen; ++i ) {
          if ( f[i] == 3 ) {
            errfound = 1;
            break;
          }
        }
        if ( errfound ) {
          iout << iERROR << "Stray PME grid charges detected: "
		<< sourcepe << " sending to " << gridPeMap[pe] << " for planes";
          int iz = -1;
          for ( i=0; i<flen; ++i ) {
            if ( f[i] == 3 ) {
              f[i] = 2;
              int jz = (i+fstart)/myGrid.K2;
              if ( iz != jz ) { iout << " " << jz;  iz = jz; }
            }
          }
          iout << "\n" << endi;
        }
      }
    }

#ifdef NETWORK_PROGRESS
    CmiNetworkProgress();
#endif

    if ( ! recipPeDest[pe] ) continue;

    int zlistlen = 0;
    for ( i=0; i<myGrid.K3; ++i ) {
      if ( fz_arr[i] ) ++zlistlen;
    }

    PmeGridMsg *msg = new (zlistlen, flen*numGrids,
				fcount*zlistlen, PRIORITY_SIZE) PmeGridMsg;

    msg->sourceNode = sourcepe;
    msg->lattice = lattice;
    msg->start = fstart;
    msg->len = flen;
    msg->zlistlen = zlistlen;
    int *zlist = msg->zlist;
    zlistlen = 0;
    for ( i=0; i<myGrid.K3; ++i ) {
      if ( fz_arr[i] ) zlist[zlistlen++] = i;
    }
    float *qmsg = msg->qgrid;
    for ( g=0; g<numGrids; ++g ) {
      char *f = f_arr + fstart + g*fsize;
      CmiMemcpy((void*)(msg->fgrid+g*flen),(void*)f,flen*sizeof(char));
      float **q = q_arr + fstart + g*fsize;
      for ( i=0; i<flen; ++i ) {
        if ( f[i] ) {
          for (int h=0; h<myGrid.order-1; ++h) {
            q[i][h] += q[i][myGrid.K3+h];
          }
          for ( int k=0; k<zlistlen; ++k ) {
            *(qmsg++) = q[i][zlist[k]];
          }
        }
      }
    }

    msg->sequence = compute_sequence;
    SET_PRIORITY(msg,compute_sequence,PME_GRID_PRIORITY)
    pmeProxy[gridPeMap[pe]].recvGrid(msg);
  }

}

void ComputePmeMgr::sendDataHelper(int iter) {
  nodePmeMgr->sendDataHelper(iter);
}

void NodePmeMgr::sendDataHelper(int iter) {
#ifdef NAMD_CUDA
  ComputePmeMgr *obj = masterPmeMgr;
  obj->sendDataPart(iter, iter, *obj->sendDataHelper_lattice, obj->sendDataHelper_sequence, obj->sendDataHelper_sourcepe, obj->sendDataHelper_errors);
#else
  NAMD_bug("NodePmeMgr::sendDataHelper called in non-CUDA build");
#endif
}

void ComputePmeMgr::sendData(Lattice &lattice, int sequence) {

  sendDataHelper_lattice = &lattice;
  sendDataHelper_sequence = sequence;
  sendDataHelper_sourcepe = CkMyPe();
  sendDataHelper_errors = strayChargeErrors;
  strayChargeErrors = 0;

#ifdef NAMD_CUDA
  if ( offload ) {
    for ( int i=0; i < numGridPes; ++i ) {
      int pe = gridPeOrder[i];  // different order
      if ( ! recipPeDest[pe] && ! sendDataHelper_errors ) continue;
#if CMK_MULTICORE
      // nodegroup messages on multicore are delivered to sending pe, or pe 0 if expedited
      pmeProxy[gridPeMap[pe]].sendDataHelper(i);
#else
      pmeNodeProxy[CkMyNode()].sendDataHelper(i);
#endif
    }
  } else
#endif
  {
    sendDataPart(0,numGridPes-1,lattice,sequence,CkMyPe(),sendDataHelper_errors);
  }
 
}

void ComputePmeMgr::copyResults(PmeGridMsg *msg) {

  int zdim = myGrid.dim3;
  int flen = msg->len;
  int fstart = msg->start;
  int zlistlen = msg->zlistlen;
  int *zlist = msg->zlist;
  float *qmsg = msg->qgrid;
  int g;
  for ( g=0; g<numGrids; ++g ) {
    char *f = msg->fgrid + g*flen;
    float **q = q_arr + fstart + g*fsize;
    for ( int i=0; i<flen; ++i ) {
      if ( f[i] ) {
        f[i] = 0;
        for ( int k=0; k<zlistlen; ++k ) {
          q[i][zlist[k]] = *(qmsg++);
        }
        for (int h=0; h<myGrid.order-1; ++h) {
          q[i][myGrid.K3+h] = q[i][h];
        }
      }
    }
  }
}

void ComputePme::ungridForces() {

    if ( sequence() != myMgr->compute_sequence ) NAMD_bug("ComputePme sequence mismatch in ungridForces()");
 
    SimParameters *simParams = Node::Object()->simParameters;

    localResults_alloc.resize(numLocalAtoms* ((numGrids>1 || selfOn)?2:1));
    Vector *localResults = localResults_alloc.begin();
    Vector *gridResults;
    if ( alchOn || lesOn || selfOn || pairOn ) {
      for(int i=0; i<numLocalAtoms; ++i) { localResults[i] = 0.; }
      gridResults = localResults + numLocalAtoms;
    } else {
      gridResults = localResults;
    }

    Vector pairForce = 0.;
    Lattice &lattice = patch->flags.lattice;
    int g = 0;
    if(!simParams->commOnly) {
    for ( g=0; g<numGrids; ++g ) {
#ifdef NETWORK_PROGRESS
      CmiNetworkProgress();
#endif

#ifdef NAMD_CUDA
      if ( offload ) {
	for ( int n=numGridAtoms[g], i=0; i<n; ++i ) {
	  gridResults[i].x = f_data_host[3*i];
	  gridResults[i].y = f_data_host[3*i+1];
	  gridResults[i].z = f_data_host[3*i+2];
	}
      } else
#endif // NAMD_CUDA
	{
	  myRealSpace[g]->compute_forces(myMgr->q_arr+g*myMgr->fsize, localGridData[g], gridResults);
	}
      scale_forces(gridResults, numGridAtoms[g], lattice);
      
      if (alchOn) {
        float scale = 1.;
        BigReal elecLambdaUp, elecLambdaDown;
        if ( simParams->alchFepWhamOn ) {
	  if ( simParams->alchFepElecOn ) {
            elecLambdaUp = simParams->alchElecLambda;
            elecLambdaDown = 1.0 - simParams->alchElecLambda;
	  }
	  else {
            elecLambdaUp = 0.0;
            elecLambdaDown = 1.0;
	  }
        }
        else {
          BigReal alchLambda = simParams->getCurrentLambda(patch->flags.step);
          myMgr->alchLambda = alchLambda;
	  elecLambdaUp = simParams->getElecLambda(alchLambda);
	  elecLambdaDown = simParams->getElecLambda(1. - alchLambda);
        }
	
        if ( g == 0 ) scale = elecLambdaUp;
        else if ( g == 1 ) scale = elecLambdaDown;
        else if ( g == 2 ) scale = (elecLambdaUp + elecLambdaDown - 1)*(-1);

        if (alchDecouple) {
          if ( g == 2 ) scale = 1 - elecLambdaUp;
          else if ( g == 3 ) scale = 1 - elecLambdaDown;
          else if ( g == 4 ) scale = (elecLambdaUp + elecLambdaDown - 1)*(-1);
        }
        int nga = 0;
        if (!alchDecouple) {
          for(int i=0; i<numLocalAtoms; ++i) {
            if ( localPartition[i] == 0 || localPartition[i] == (g+1) ) {
              // (g=2: only partition 0)
              localResults[i] += gridResults[nga++] * scale;
            }
          }
        }
        else {  // alchDecouple
          if ( g < 2 ) {
            for(int i=0; i<numLocalAtoms; ++i) {
              if ( localPartition[i] == 0 || localPartition[i] == (g+1) ) {
                // g = 0: partition 0 or partition 1
                // g = 1: partition 0 or partition 2
                localResults[i] += gridResults[nga++] * scale;
              }
            }
          }
          else {
            for(int i=0; i<numLocalAtoms; ++i) {
              if ( localPartition[i] == (g-1) || localPartition[i] == (g-4)) {
                // g = 2: partition 1 only
                // g = 3: partition 2 only
                // g = 4: partition 0 only
                localResults[i] += gridResults[nga++] * scale;
              }
            }
          }
        }
      } else if ( lesOn ) {
        float scale = 1.;
        if ( alchFepOn ) {
	  if(simParams->alchFepWhamOn) {
	    if(simParams->alchFepElecOn) {
	      if ( g == 0 ) scale = simParams->alchElecLambda;
	      else if ( g == 1 ) scale = 1. - simParams->alchElecLambda;
	    }
	    else {
	      if ( g == 0 ) scale = 0.0;
	      else if ( g == 1 ) scale = 1.0;
	    }
	  }
	  else {
            BigReal alchLambda = simParams->getCurrentLambda(patch->flags.step);
            myMgr->alchLambda = alchLambda;
            if ( g == 0 ) scale = alchLambda;
            else if ( g == 1 ) scale = 1. - alchLambda;
	  }
        } else if ( lesOn ) {
          scale = 1.0 / (float)lesFactor;
        }
        int nga = 0;
        for(int i=0; i<numLocalAtoms; ++i) {
          if ( localPartition[i] == 0 || localPartition[i] == (g+1) ) {
            localResults[i] += gridResults[nga++] * scale;
          }
        }
      } else if ( selfOn ) {
        PmeParticle *lgd = localGridData[g];
        int nga = 0;
        for(int i=0; i<numLocalAtoms; ++i) {
          if ( localPartition[i] == 1 ) {
            pairForce += gridResults[nga];  // should add up to almost zero
            localResults[i] += gridResults[nga++];
          }
        }
      } else if ( pairOn ) {
        if ( g == 0 ) {
          int nga = 0;
          for(int i=0; i<numLocalAtoms; ++i) {
            if ( localPartition[i] == 1 ) {
              pairForce += gridResults[nga];
            }
            if ( localPartition[i] == 1 || localPartition[i] == 2 ) {
              localResults[i] += gridResults[nga++];
            }
          }
        } else if ( g == 1 ) {
          int nga = 0;
          for(int i=0; i<numLocalAtoms; ++i) {
            if ( localPartition[i] == g ) {
              pairForce -= gridResults[nga];  // should add up to almost zero
              localResults[i] -= gridResults[nga++];
            }
          }
        } else {
          int nga = 0;
          for(int i=0; i<numLocalAtoms; ++i) {
            if ( localPartition[i] == g ) {
              localResults[i] -= gridResults[nga++];
            }
         }
        }
      }
    }
    }

    Vector *results_ptr = localResults;
    
    // add in forces
    {
      Results *r = forceBox->open();
      Force *f = r->f[Results::slow];
      int numAtoms = patch->getNumAtoms();

      if ( ! myMgr->strayChargeErrors && ! simParams->commOnly ) {
        for(int i=0; i<numAtoms; ++i) {
          f[i].x += results_ptr->x;
          f[i].y += results_ptr->y;
          f[i].z += results_ptr->z;
          ++results_ptr;
        }
      }
      forceBox->close(&r);
    }

    if ( pairOn || selfOn ) {
        ADD_VECTOR_OBJECT(myMgr->reduction,REDUCTION_PAIR_ELECT_FORCE,pairForce);
    }

}

void ComputePmeMgr::submitReductions() {

    SimParameters *simParams = Node::Object()->simParameters;

    for ( int g=0; g<numGrids; ++g ) {
      float scale = 1.;
      if (alchOn) {
        BigReal elecLambdaUp, elecLambdaDown;
        if( simParams->alchFepWhamOn ) {
          if( simParams->alchFepElecOn ) {
            elecLambdaUp = simParams->alchElecLambda;
            elecLambdaDown = 1.0 - simParams->alchElecLambda;
          }
          else {
            elecLambdaUp = 0.0;
            elecLambdaDown = 1.0;
          }
        }
        else {
          // alchLambda set on each step in ComputePme::ungridForces()
          if ( alchLambda < 0 || alchLambda > 1 ) {
            NAMD_bug("ComputePmeMgr::submitReductions alchLambda out of range");
          }
          elecLambdaUp = simParams->getElecLambda(alchLambda);
          elecLambdaDown = simParams->getElecLambda(1-alchLambda);
        }
        if ( g == 0 ) scale = elecLambdaUp;
        else if ( g == 1 ) scale = elecLambdaDown;
        else if ( g == 2 ) scale = (elecLambdaUp + elecLambdaDown - 1)*(-1);
        if (alchDecouple) {
          if ( g == 2 ) scale = 1-elecLambdaUp;
          else if ( g == 3 ) scale = 1-elecLambdaDown;
          else if ( g == 4 ) scale = (elecLambdaUp + elecLambdaDown - 1)*(-1);
        }
      } else if ( lesOn ) {
        scale = 1.0 / lesFactor;
      } else if ( pairOn ) {
        scale = ( g == 0 ? 1. : -1. );
      }
      reduction->item(REDUCTION_ELECT_ENERGY_SLOW) += evir[g][0] * scale;
      reduction->item(REDUCTION_VIRIAL_SLOW_XX) += evir[g][1] * scale;
      reduction->item(REDUCTION_VIRIAL_SLOW_XY) += evir[g][2] * scale;
      reduction->item(REDUCTION_VIRIAL_SLOW_XZ) += evir[g][3] * scale;
      reduction->item(REDUCTION_VIRIAL_SLOW_YX) += evir[g][2] * scale;
      reduction->item(REDUCTION_VIRIAL_SLOW_YY) += evir[g][4] * scale;
      reduction->item(REDUCTION_VIRIAL_SLOW_YZ) += evir[g][5] * scale;
      reduction->item(REDUCTION_VIRIAL_SLOW_ZX) += evir[g][3] * scale;
      reduction->item(REDUCTION_VIRIAL_SLOW_ZY) += evir[g][5] * scale;
      reduction->item(REDUCTION_VIRIAL_SLOW_ZZ) += evir[g][6] * scale;

      float scale2 = 0.;

      // why is this declared/defined again here?
      SimParameters *simParams = Node::Object()->simParameters;

      if (alchFepOn) {
      	BigReal elecLambda2Up=0.0, elecLambda2Down=0.0;
        if(simParams->alchFepWhamOn) {
          if(simParams->alchFepElecOn) {
            elecLambda2Up = simParams->alchElecLambda;
            elecLambda2Down =  1.0 - simParams->alchElecLambda;
          }
          else {
            elecLambda2Up = 0.0;
            elecLambda2Down =  1.0;
          }
        }
        else {
          elecLambda2Up = simParams->getElecLambda(simParams->alchLambda2);
          elecLambda2Down = simParams->getElecLambda(1.-simParams->alchLambda2);
        }
        
        if ( g == 0 ) scale2 = elecLambda2Up;
        else if ( g == 1 ) scale2 = elecLambda2Down;
        else if ( g == 2 ) scale2 = (elecLambda2Up + elecLambda2Down - 1)*(-1);
        if (alchDecouple && g == 2 ) scale2 = 1 - elecLambda2Up;
        else if (alchDecouple && g == 3 ) scale2 = 1 - elecLambda2Down;
        else if (alchDecouple && g == 4 ) scale2 = (elecLambda2Up + elecLambda2Down - 1)*(-1);
      }
      if(simParams->alchFepWhamOn && simParams->alchFepElecOn)	{	// FEP with wham post-process
      	if( g==0 )	scale2 = scale + 1.0;
      	else if( g==1 )	scale2 = scale - 1.0;
      	else if( g==2 )	scale2 = scale - 1.0;
      	else if( g==3 )	scale2 = scale + 1.0;
      }
      reduction->item(REDUCTION_ELECT_ENERGY_SLOW_F) += evir[g][0] * scale2;
      
      if (alchThermIntOn) {
        
        // no decoupling:
        // part. 1 <-> all of system except partition 2: g[0] - g[2] 
        // (interactions between all atoms [partition 0 OR partition 1], 
        // minus all [within partition 0])
        // U = elecLambdaUp * (U[0] - U[2])
        // dU/dl = U[0] - U[2];
        
        // part. 2 <-> all of system except partition 1: g[1] - g[2] 
        // (interactions between all atoms [partition 0 OR partition 2], 
        // minus all [within partition 0])
        // U = elecLambdaDown * (U[1] - U[2])
        // dU/dl = U[1] - U[2];

        // alchDecouple:
        // part. 1 <-> part. 0: g[0] - g[2] - g[4] 
        // (interactions between all atoms [partition 0 OR partition 1]
        // minus all [within partition 1] minus all [within partition 0]
        // U = elecLambdaUp * (U[0] - U[4]) + (1-elecLambdaUp)* U[2]
        // dU/dl = U[0] - U[2] - U[4];

        // part. 2 <-> part. 0: g[1] - g[3] - g[4] 
        // (interactions between all atoms [partition 0 OR partition 2]
        // minus all [within partition 2] minus all [within partition 0]
        // U = elecLambdaDown * (U[1] - U[4]) + (1-elecLambdaDown)* U[3]
        // dU/dl = U[1] - U[3] - U[4];
        
        
        if ( g == 0 ) reduction->item(REDUCTION_ELECT_ENERGY_PME_TI_1) += evir[g][0];
        if ( g == 1 ) reduction->item(REDUCTION_ELECT_ENERGY_PME_TI_2) += evir[g][0];
        if (!alchDecouple) {
          if ( g == 2 ) reduction->item(REDUCTION_ELECT_ENERGY_PME_TI_1) -= evir[g][0];
          if ( g == 2 ) reduction->item(REDUCTION_ELECT_ENERGY_PME_TI_2) -= evir[g][0];
        }
        else {  // alchDecouple
          if ( g == 2 ) reduction->item(REDUCTION_ELECT_ENERGY_PME_TI_1) -= evir[g][0];
          if ( g == 3 ) reduction->item(REDUCTION_ELECT_ENERGY_PME_TI_2) -= evir[g][0];
          if ( g == 4 ) reduction->item(REDUCTION_ELECT_ENERGY_PME_TI_1) -= evir[g][0];
          if ( g == 4 ) reduction->item(REDUCTION_ELECT_ENERGY_PME_TI_2) -= evir[g][0];
        }
      }
    }

    alchLambda = -1.;  // illegal value to catch if not updated

    reduction->item(REDUCTION_STRAY_CHARGE_ERRORS) += strayChargeErrors;
    reduction->submit();

  for ( int i=0; i<heldComputes.size(); ++i ) {
    WorkDistrib::messageEnqueueWork(heldComputes[i]);
  }
  heldComputes.resize(0);
}

#if USE_TOPOMAP 

#define NPRIMES 8
const static unsigned int NAMDPrimes[] = {
  3,
  5,
  7,
  11,
  13,
  17,
  19,
  23,  
  29,
  31,
  37,
  59,
  73,
  93,
  113,
  157,
  307,
  617,
  1217                  //This should b enough for 64K nodes of BGL. 
};

#include "RecBisection.h"

/***-----------------------------------------------------**********
    The Orthogonal Recursive Bisection strategy, which allocates PME
    objects close to the patches they communicate, and at the same
    time spreads them around the grid 
****----------------------------------------------------------****/

bool generateBGLORBPmePeList(int *pemap, int numPes, 
			     int *block_pes, int nbpes) {

  PatchMap *pmap = PatchMap::Object();
  int *pmemap = new int [CkNumPes()];

  if (pemap == NULL)
    return false;

  TopoManager tmgr;

  memset(pmemap, 0, sizeof(int) * CkNumPes());

  for(int count = 0; count < CkNumPes(); count++) {
    if(count < nbpes)
      pmemap[block_pes[count]] = 1;
    
    if(pmap->numPatchesOnNode(count)) {
      pmemap[count] = 1;
      
      //Assumes an XYZT mapping !!
      if(tmgr.hasMultipleProcsPerNode()) {
	pmemap[(count + CkNumPes()/2)% CkNumPes()] = 1;
      }
    }
  }

  if(numPes + nbpes + pmap->numNodesWithPatches() > CkNumPes())
    //NAMD_bug("PME ORB Allocator: Processors Unavailable\n");
    return false;

  CProxy_Node nd(CkpvAccess(BOCclass_group).node);
  Node *node = nd.ckLocalBranch();
  SimParameters *simParams = node->simParameters;

  //first split PME processors into patch groups

  int xsize = 0, ysize = 0, zsize = 0;

  xsize = tmgr.getDimNX();
  ysize = tmgr.getDimNY();
  zsize = tmgr.getDimNZ();
  
  int nx = xsize, ny = ysize, nz = zsize;
  DimensionMap dm;
  
  dm.x = 0;
  dm.y = 1;
  dm.z = 2;
  
  findOptimalDimensions(xsize, ysize, zsize, nx, ny, nz, dm);

  //group size processors have to be allocated to each YZ plane
  int group_size = numPes/nx;
  if(numPes % nx)
    group_size ++;

  int my_prime = NAMDPrimes[0];
  int density = (ny * nz)/group_size + 1;
  int count = 0;
  
  //Choose a suitable prime Number
  for(count = 0; count < NPRIMES; count ++) {
    //Find a prime just greater than the density
    if(density < NAMDPrimes[count]) {
      my_prime = NAMDPrimes[count];
      break;
    }      
  }
  
  if(count == NPRIMES)
    my_prime = NAMDPrimes[NPRIMES-1];

  //int gcount = numPes/2;
  int gcount = 0;
  int npme_pes = 0;
  
  int coord[3];

  for(int x = 0; x < nx; x++) {
    coord[0] = (x + nx/2)%nx;
    
    for(count=0; count < group_size && npme_pes < numPes; count++) {
      int dest = (count + 1) * my_prime;      
      dest = dest % (ny * nz);      
      
      coord[2] = dest / ny;
      coord[1] = dest - coord[2] * ny;
      
      //Locate where in the actual grid the processor is
      int destPe = coord[dm.x] + coord[dm.y] * xsize + 
	coord[dm.z] * xsize* ysize;
      
      if(pmemap[destPe] == 0) {
        pemap[gcount++] = destPe;
        pmemap[destPe] = 1;
	
	if(tmgr.hasMultipleProcsPerNode())
	  pmemap[(destPe + CkNumPes()/2) % CkNumPes()] = 1;	

        npme_pes ++;
      }
      else {
        for(int pos = 1; pos < ny * nz; pos++) {
          
          coord[2] += pos / ny;
          coord[1] += pos % ny;
          
          coord[2] = coord[2] % nz;
          coord[1] = coord[1] % ny;       
          
          int newdest = coord[dm.x] + coord[dm.y] * xsize + 
	    coord[dm.z] * xsize * ysize;
          
          if(pmemap[newdest] == 0) {
            pemap[gcount++] = newdest;
            pmemap[newdest] = 1;
	    
	    if(tmgr.hasMultipleProcsPerNode())
	      pmemap[(newdest + CkNumPes()/2) % CkNumPes()] = 1;	
	    
            npme_pes ++;
            break;
          }
        }
      }      
    }   
    
    if(gcount == numPes)
      gcount = 0;    
    
    if(npme_pes >= numPes)
      break;
  }
  
  delete [] pmemap;
  
  if(npme_pes != numPes)
    //NAMD_bug("ORB PME allocator failed\n");
    return false;

  return true;
}

#endif

template <class T> class PmePencil : public T {
public:
  PmePencil() {
    data = 0;
    work = 0;
    send_order = 0;
    needs_reply = 0;
#if USE_PERSISTENT
    trans_handle = untrans_handle = ungrid_handle = NULL;
#endif
  }
  ~PmePencil() {
#ifdef NAMD_FFTW
    fftwf_free(data);
#endif
    delete [] work;
    delete [] send_order;
    delete [] needs_reply;
  }
  void base_init(PmePencilInitMsg *msg) {
    imsg=0;
    imsgb=0;
    hasData=0;
    initdata = msg->data;
  }
  void order_init(int nBlocks) {
    send_order = new int[nBlocks];
    for ( int i=0; i<nBlocks; ++i ) send_order[i] = i;
    if ( Node::Object()->simParameters->PMESendOrder ) {
      std::sort(send_order,send_order+nBlocks,sortop_bit_reversed());
    } else {
      Random rand(CkMyPe());
      rand.reorder(send_order,nBlocks);
    }
    needs_reply = new int[nBlocks];
    offload = Node::Object()->simParameters->PMEOffload;
  }
  PmePencilInitMsgData initdata;
  Lattice lattice;
  PmeReduction evir;
  int sequence;  // used for priorities
  int imsg;  // used in sdag code
  int imsgb;  // Node par uses distinct counter for back path
  int hasData;  // used in message elimination
  int offload;
  float *data;
  float *work;
  int *send_order;
  int *needs_reply;
#if USE_PERSISTENT
  PersistentHandle *trans_handle;
  PersistentHandle *untrans_handle;
  PersistentHandle *ungrid_handle;
#endif
};

class PmeZPencil : public PmePencil<CBase_PmeZPencil> {
public:
    PmeZPencil_SDAG_CODE
    PmeZPencil() { __sdag_init(); setMigratable(false); }
    PmeZPencil(CkMigrateMessage *) { __sdag_init();  setMigratable (false); imsg=imsgb=0;}
	~PmeZPencil() {
	#ifdef NAMD_FFTW
	#ifdef NAMD_FFTW_3
		delete [] forward_plans;
		delete [] backward_plans;
	#endif
	#endif
	}
    void fft_init();
    void recv_grid(const PmeGridMsg *);
    void forward_fft();
    void send_trans();
	void send_subset_trans(int fromIdx, int toIdx);
    void recv_untrans(const PmeUntransMsg *);
    void node_process_untrans(PmeUntransMsg *);
    void node_process_grid(PmeGridMsg *);
    void backward_fft();
	void send_ungrid(PmeGridMsg *);
	void send_all_ungrid();
	void send_subset_ungrid(int fromIdx, int toIdx, int specialIdx);
private:
    ResizeArray<PmeGridMsg *> grid_msgs;
    ResizeArray<int> work_zlist;
#ifdef NAMD_FFTW
#ifdef NAMD_FFTW_3
    fftwf_plan forward_plan, backward_plan;

	//for ckloop usage
	int numPlans;
	fftwf_plan *forward_plans, *backward_plans;
#else
    rfftwnd_plan forward_plan, backward_plan;
#endif
#endif

    int nx, ny;
#if USE_PERSISTENT
    void setup_persistent() {
      int hd = 1;// ( hasData ? 1 : 0 );
      int zBlocks = initdata.zBlocks;
      int block3 = initdata.grid.block3;
      int dim3 = initdata.grid.dim3;
      CkArray *yPencil_local = initdata.yPencil.ckLocalBranch();
      CmiAssert(yPencil_local);
      trans_handle = (PersistentHandle*) malloc( sizeof(PersistentHandle) * zBlocks);
      for ( int isend=0; isend<zBlocks; ++isend ) {
          int kb = send_order[isend];
          int nz1 = block3;
          if ( (kb+1)*block3 > dim3/2 ) nz1 = dim3/2 - kb*block3;
          int peer = yPencil_local->homePe(CkArrayIndex3D(thisIndex.x, 0, kb));
          int size = sizeof(PmeTransMsg) + sizeof(float)*hd*nx*ny*nz1*2 +sizeof( envelope)+PRIORITY_SIZE/8+24;
          int compress_start = sizeof(PmeTransMsg)+sizeof(envelope);
          int compress_size = sizeof(float)*hd*nx*ny*nz1*2;
          trans_handle[isend] = CmiCreateCompressPersistentSize(peer, size, compress_start, compress_size, CMI_FLOATING);
      }
    }
    
    void setup_ungrid_persistent() 
    {
       ungrid_handle = (PersistentHandle*) malloc( sizeof(PersistentHandle) * grid_msgs.size());
       for ( imsg=0; imsg < grid_msgs.size(); ++imsg ) {
           int peer = grid_msgs[imsg]->sourceNode;
           //ungrid_handle[imsg] = CmiCreatePersistent(peer, 0); 
       }
    }
#endif
};

class PmeYPencil : public PmePencil<CBase_PmeYPencil> {
public:
    PmeYPencil_SDAG_CODE
    PmeYPencil() { __sdag_init(); setMigratable(false); imsg=imsgb=0;}
    PmeYPencil(CkMigrateMessage *) { __sdag_init(); }
    void fft_init();
    void recv_trans(const PmeTransMsg *);
    void forward_fft();
	void forward_subset_fft(int fromIdx, int toIdx);
    void send_trans();
	void send_subset_trans(int fromIdx, int toIdx);
    void recv_untrans(const PmeUntransMsg *);    
    void node_process_trans(PmeTransMsg *);
    void node_process_untrans(PmeUntransMsg *);
    void backward_fft();
	void backward_subset_fft(int fromIdx, int toIdx);
    void send_untrans();
    void send_subset_untrans(int fromIdx, int toIdx, int evirIdx);
private:
#ifdef NAMD_FFTW
#ifdef NAMD_FFTW_3
    fftwf_plan forward_plan, backward_plan;
#else
    fftw_plan forward_plan, backward_plan;
#endif
#endif

    int nx, nz;
#if USE_PERSISTENT
    void setup_persistent() {
      int yBlocks = initdata.yBlocks;
      int block2 = initdata.grid.block2;
      int K2 = initdata.grid.K2;
      int hd = 1;
      CkArray *xPencil_local = initdata.xPencil.ckLocalBranch();
      trans_handle = (PersistentHandle*) malloc( sizeof(PersistentHandle) * yBlocks);
      for ( int isend=0; isend<yBlocks; ++isend ) {
          int jb = send_order[isend];
          int ny1 = block2;
          if ( (jb+1)*block2 > K2 ) ny1 = K2 - jb*block2;
          int peer = xPencil_local->homePe(CkArrayIndex3D(0, jb, thisIndex.z));
          int size = sizeof(PmeTransMsg) + sizeof(float)*hd*nx*ny1*nz*2 +sizeof( envelope) + PRIORITY_SIZE/8+24;
          int compress_start = sizeof(PmeTransMsg)+sizeof( envelope);
          int compress_size = sizeof(float)*hd*nx*ny1*nz*2; 
          trans_handle[isend] = CmiCreateCompressPersistentSize(peer, size, compress_start, compress_size, CMI_FLOATING);
      }

      CkArray *zPencil_local = initdata.zPencil.ckLocalBranch();
      untrans_handle = (PersistentHandle*) malloc( sizeof(PersistentHandle) * yBlocks);
      for ( int isend=0; isend<yBlocks; ++isend ) {
          int jb = send_order[isend];
          int ny1 = block2;
          if ( (jb+1)*block2 > K2 ) ny1 = K2 - jb*block2;
          int peer = zPencil_local->homePe(CkArrayIndex3D(thisIndex.x, jb, 0));
          int size= sizeof(PmeUntransMsg) + sizeof(float)*nx*ny1*nz*2 + sizeof( envelope) + PRIORITY_SIZE/8+24;
          int compress_start = sizeof(PmeUntransMsg) + sizeof( envelope); 
          int compress_size = sizeof(float)*nx*ny1*nz*2;
          untrans_handle[isend] = CmiCreateCompressPersistentSize(peer, size,  compress_start, compress_size, CMI_FLOATING);
      }
    }
#endif
};

class PmeXPencil : public PmePencil<CBase_PmeXPencil> {
public:
    PmeXPencil_SDAG_CODE
    PmeXPencil() { __sdag_init();  myKSpace = 0; setMigratable(false); imsg=imsgb=0; recipEvirPe = -999; }
    PmeXPencil(CkMigrateMessage *) { __sdag_init(); }
	~PmeXPencil() {
	#ifdef NAMD_FFTW
	#ifdef NAMD_FFTW_3
		delete [] forward_plans;
		delete [] backward_plans;
	#endif
	#endif
	}
    void fft_init();
    void recv_trans(const PmeTransMsg *);
    void forward_fft();
    void pme_kspace();
    void backward_fft();
    void send_untrans();
	void send_subset_untrans(int fromIdx, int toIdx, int evirIdx);
    void node_process_trans(PmeTransMsg *);
#ifdef NAMD_FFTW
#ifdef NAMD_FFTW_3
    fftwf_plan forward_plan, backward_plan;

	int numPlans;
	fftwf_plan *forward_plans, *backward_plans;
#else
    fftw_plan forward_plan, backward_plan;
#endif
#endif
    int ny, nz;
    int recipEvirPe;
    void evir_init();
    PmeKSpace *myKSpace;
#if USE_PERSISTENT
    void  setup_persistent() {
      int xBlocks = initdata.xBlocks;
      int block1 = initdata.grid.block1;
      int K1 = initdata.grid.K1;
      CkArray *yPencil_local = initdata.yPencil.ckLocalBranch();
      untrans_handle = (PersistentHandle*) malloc( sizeof(PersistentHandle) * xBlocks);
      for ( int isend=0; isend<xBlocks; ++isend ) {
          int ib = send_order[isend];
          int nx1 = block1;
          if ( (ib+1)*block1 > K1 ) nx1 = K1 - ib*block1;
          int peer = yPencil_local->procNum(CkArrayIndex3D(ib, 0, thisIndex.z));
          int size = sizeof(PmeUntransMsg) +
              sizeof(float)*nx1*ny*nz*2 +sizeof( envelope) + PRIORITY_SIZE/8+24; 
          int compress_start = sizeof(PmeUntransMsg) + sizeof( envelope); 
          int compress_size = sizeof(float)*nx1*ny*nz*2;
          untrans_handle[isend] = CmiCreateCompressPersistentSize(peer, size, compress_start, compress_size, CMI_FLOATING);
      }
    }
#endif

};

void PmeXPencil::evir_init() {
  recipEvirPe = findRecipEvirPe();
  initdata.pmeProxy[recipEvirPe].addRecipEvirClient();
}

void PmeZPencil::fft_init() {
  CProxy_Node nd(CkpvAccess(BOCclass_group).node);
  Node *node = nd.ckLocalBranch();
  SimParameters *simParams = node->simParameters;

#if USE_NODE_PAR_RECEIVE
  ((NodePmeMgr *)CkLocalNodeBranch(initdata.pmeNodeProxy))->registerZPencil(thisIndex,this);
#endif

  int K1 = initdata.grid.K1;
  int K2 = initdata.grid.K2;
  int K3 = initdata.grid.K3;
  int dim3 = initdata.grid.dim3;
  int block1 = initdata.grid.block1;
  int block2 = initdata.grid.block2;

  nx = block1;
  if ( (thisIndex.x + 1) * block1 > K1 ) nx = K1 - thisIndex.x * block1;
  ny = block2;
  if ( (thisIndex.y + 1) * block2 > K2 ) ny = K2 - thisIndex.y * block2;

#ifdef NAMD_FFTW
  CmiLock(ComputePmeMgr::fftw_plan_lock);

  data = (float *) fftwf_malloc( sizeof(float) *nx*ny*dim3);
  work = new float[dim3];

  order_init(initdata.zBlocks);

#ifdef NAMD_FFTW_3
  /* need array of sizes for the how many */

  int fftwFlags = simParams->FFTWPatient ? FFTW_PATIENT  : simParams->FFTWEstimate ? FFTW_ESTIMATE  : FFTW_MEASURE ;
  int sizeLines=nx*ny;
  int planLineSizes[1];
  planLineSizes[0]=K3;
  int ndim=initdata.grid.dim3; // storage space is initdata.grid.dim3
  int ndimHalf=ndim/2;
  forward_plan = fftwf_plan_many_dft_r2c(1, planLineSizes, sizeLines,
					 (float *) data, NULL, 1, 
					 ndim,
					 (fftwf_complex *) data, NULL, 1,
					 ndimHalf,
					 fftwFlags);

  backward_plan = fftwf_plan_many_dft_c2r(1, planLineSizes, sizeLines,
					  (fftwf_complex *) data, NULL, 1, 
					  ndimHalf,
					  (float *) data, NULL, 1, 
					  ndim,
					  fftwFlags);
#if     CMK_SMP && USE_CKLOOP
  if(simParams->useCkLoop) {
	  //How many FFT plans to be created? The grain-size issue!!.
	  //Currently, I am choosing the min(nx, ny) to be coarse-grain
	  numPlans = (nx<=ny?nx:ny);
          if ( numPlans < CkMyNodeSize() ) numPlans = (nx>=ny?nx:ny);
          if ( numPlans < CkMyNodeSize() ) numPlans = sizeLines;
	  int howmany = sizeLines/numPlans;
	  forward_plans = new fftwf_plan[numPlans];
	  backward_plans = new fftwf_plan[numPlans];
	  for(int i=0; i<numPlans; i++) {
		  int dimStride = i*ndim*howmany;
		  int dimHalfStride = i*ndimHalf*howmany;
		  forward_plans[i] = fftwf_plan_many_dft_r2c(1, planLineSizes, howmany,
													 ((float *)data)+dimStride, NULL, 1,
													 ndim,
													 ((fftwf_complex *)data)+dimHalfStride, NULL, 1,
													 ndimHalf,
													 fftwFlags);

		  backward_plans[i] = fftwf_plan_many_dft_c2r(1, planLineSizes, howmany,
													 ((fftwf_complex *)data)+dimHalfStride, NULL, 1,
													 ndimHalf,
													 ((float *)data)+dimStride, NULL, 1,
													 ndim,
													 fftwFlags);
	  }
  }else 
#endif 
  {
	  forward_plans = NULL;
	  backward_plans = NULL;
  }
#else
  forward_plan = rfftwnd_create_plan_specific(1, &K3, FFTW_REAL_TO_COMPLEX,
	( simParams->FFTWEstimate ? FFTW_ESTIMATE : FFTW_MEASURE )
	| FFTW_IN_PLACE | FFTW_USE_WISDOM, data, 1, work, 1);
  backward_plan = rfftwnd_create_plan_specific(1, &K3, FFTW_COMPLEX_TO_REAL,
	( simParams->FFTWEstimate ? FFTW_ESTIMATE : FFTW_MEASURE )
	| FFTW_IN_PLACE | FFTW_USE_WISDOM, data, 1, work, 1);
#endif
  CmiUnlock(ComputePmeMgr::fftw_plan_lock);
#else
  NAMD_die("Sorry, FFTW must be compiled in to use PME.");
#endif

#if USE_NODE_PAR_RECEIVE
    evir = 0.;
    memset(data, 0, sizeof(float) * nx*ny*dim3);
#endif
}

void PmeYPencil::fft_init() {
  CProxy_Node nd(CkpvAccess(BOCclass_group).node);
  Node *node = nd.ckLocalBranch();
  SimParameters *simParams = node->simParameters;

#if USE_NODE_PAR_RECEIVE
  ((NodePmeMgr *)CkLocalNodeBranch(initdata.pmeNodeProxy))->registerYPencil(thisIndex,this);
#endif

  int K1 = initdata.grid.K1;
  int K2 = initdata.grid.K2;
  int dim2 = initdata.grid.dim2;
  int dim3 = initdata.grid.dim3;
  int block1 = initdata.grid.block1;
  int block3 = initdata.grid.block3;

  nx = block1;
  if ( (thisIndex.x + 1) * block1 > K1 ) nx = K1 - thisIndex.x * block1;
  nz = block3;
  if ( (thisIndex.z+1)*block3 > dim3/2 ) nz = dim3/2 - thisIndex.z*block3;

#ifdef NAMD_FFTW
  CmiLock(ComputePmeMgr::fftw_plan_lock);

  data = (float *) fftwf_malloc( sizeof(float) * nx*dim2*nz*2);
  work = new float[2*K2];

  order_init(initdata.yBlocks);

#ifdef NAMD_FFTW_3
  /* need array of sizes for the dimensions */
  /* ideally this should be implementable as a single multidimensional
   *  plan, but that has proven tricky to implement, so we maintain the
   *  loop of 1d plan executions. */
  int sizeLines=nz;
  int planLineSizes[1];
  planLineSizes[0]=K2;
  int fftwFlags = simParams->FFTWPatient ? FFTW_PATIENT  : simParams->FFTWEstimate ? FFTW_ESTIMATE  : FFTW_MEASURE ;
  forward_plan = fftwf_plan_many_dft(1, planLineSizes, sizeLines, 
				     (fftwf_complex *) data, NULL, sizeLines, 1,
				     (fftwf_complex *) data, NULL, sizeLines, 1,
				     FFTW_FORWARD, 
				     fftwFlags);
  backward_plan = fftwf_plan_many_dft(1, planLineSizes, sizeLines, 
				     (fftwf_complex *) data, NULL, sizeLines, 1,
				     (fftwf_complex *) data, NULL, sizeLines, 1,
				     FFTW_BACKWARD, 
				      fftwFlags);
  CkAssert(forward_plan != NULL);
  CkAssert(backward_plan != NULL);
#else
  forward_plan = fftw_create_plan_specific(K2, FFTW_FORWARD,
	( simParams->FFTWEstimate ? FFTW_ESTIMATE : FFTW_MEASURE )
	| FFTW_IN_PLACE | FFTW_USE_WISDOM, (fftw_complex *) data,
	nz, (fftw_complex *) work, 1);
  backward_plan = fftw_create_plan_specific(K2, FFTW_BACKWARD,
	( simParams->FFTWEstimate ? FFTW_ESTIMATE : FFTW_MEASURE )
	| FFTW_IN_PLACE | FFTW_USE_WISDOM, (fftw_complex *) data,
	nz, (fftw_complex *) work, 1);
#endif
  CmiUnlock(ComputePmeMgr::fftw_plan_lock);
#else
  NAMD_die("Sorry, FFTW must be compiled in to use PME.");
#endif

#if USE_NODE_PAR_RECEIVE
  evir = 0;
  CmiMemoryWriteFence();
#endif
}

void PmeYPencil::node_process_trans(PmeTransMsg *msg)
{
  if ( msg->hasData ) hasData = 1;
  needs_reply[msg->sourceNode] = msg->hasData;
  recv_trans(msg);
  int limsg;
  CmiMemoryAtomicFetchAndInc(imsg,limsg);
  if(limsg+1 == initdata.yBlocks)
    {
      if ( hasData ) {
        forward_fft();
      }
      send_trans();
      if( ! hasData)
        {
          send_untrans(); //todo, what is up with the recvAck in SDAG version?
        }
      imsg=0;
      CmiMemoryWriteFence();
    }
}

void PmeYPencil::node_process_untrans(PmeUntransMsg *msg)
{
  recv_untrans(msg);
  int limsg;
  CmiMemoryAtomicFetchAndInc(imsgb,limsg);
  if(limsg+1 == initdata.yBlocks)
    {
      backward_fft();
      send_untrans();
      imsgb=0;
      CmiMemoryWriteFence();
    }
}

#define DEBUG_NODE_PAR_RECV 0

void NodePmeMgr::recvXTrans(PmeTransMsg *msg) {
  //  CkPrintf("[%d] NodePmeMgr recvXTrans for %d %d %d\n",CkMyPe(),msg->destElem.index[0],msg->destElem.index[1],msg->destElem.index[2]);
  PmeXPencil *target=xPencilObj.get(msg->destElem);
#if DEBUG_NODE_PAR_RECV
  if(target == NULL)
    CkAbort("xpencil in recvXTrans not found, debug registeration");
#endif  
    target->node_process_trans(msg);
  delete msg;
}


void NodePmeMgr::recvYTrans(PmeTransMsg *msg) {
  //  CkPrintf("[%d] NodePmeMgr recvYTrans for %d %d %d\n",CkMyPe(),msg->destElem.index[0],msg->destElem.index[1],msg->destElem.index[2]);
  PmeYPencil *target=yPencilObj.get(msg->destElem);
#if DEBUG_NODE_PAR_RECV
  if(target == NULL)
    CkAbort("ypencil in recvYTrans not found, debug registeration");
#endif  
    target->node_process_trans(msg);
  delete msg;
 }
void NodePmeMgr::recvYUntrans(PmeUntransMsg *msg) {
  //  CkPrintf("[%d] NodePmeMgr recvYUntrans for %d %d %d\n",CkMyPe(),msg->destElem.index[0],msg->destElem.index[1],msg->destElem.index[2]);
  PmeYPencil *target=yPencilObj.get(msg->destElem);
#if DEBUG_NODE_PAR_RECV  
  if(target == NULL)
    CkAbort("ypencil in recvYUntrans not found, debug registeration");
#endif  
    target->node_process_untrans(msg);
  delete msg;
 }
void NodePmeMgr::recvZUntrans(PmeUntransMsg *msg) {
  //CkPrintf("[%d] NodePmeMgr recvZUntrans for %d %d %d\n",CkMyPe(),msg->destElem.index[0],msg->destElem.index[1],msg->destElem.index[2]);
  PmeZPencil *target=zPencilObj.get(msg->destElem);
#if DEBUG_NODE_PAR_RECV
  if(target == NULL)
    CkAbort("zpencil in recvZUntrans not found, debug registeration");
#endif
  target->node_process_untrans(msg);
  delete msg;
}

void NodePmeMgr::recvZGrid(PmeGridMsg *msg) {
  //CkPrintf("[%d] NodePmeMgr %p recvGrid for %d %d %d\n",CkMyPe(),this,msg->destElem.index[0],msg->destElem.index[1],msg->destElem.index[2]);
  PmeZPencil *target=zPencilObj.get(msg->destElem);
#if DEBUG_NODE_PAR_RECV
  if(target == NULL){
    CkAbort("zpencil in recvZGrid not found, debug registeration");
  }
#endif
  target->node_process_grid(msg); //msg is stored inside node_proces_grid
}

void PmeXPencil::fft_init() {
  CProxy_Node nd(CkpvAccess(BOCclass_group).node);
  Node *node = nd.ckLocalBranch();
  SimParameters *simParams = node->simParameters;
#if USE_NODE_PAR_RECEIVE
  ((NodePmeMgr *)CkLocalNodeBranch(initdata.pmeNodeProxy))->registerXPencil(thisIndex,this);
#endif

  int K1 = initdata.grid.K1;
  int K2 = initdata.grid.K2;
  int dim3 = initdata.grid.dim3;
  int block2 = initdata.grid.block2;
  int block3 = initdata.grid.block3;

  ny = block2;
  if ( (thisIndex.y + 1) * block2 > K2 ) ny = K2 - thisIndex.y * block2;
  nz = block3;
  if ( (thisIndex.z+1)*block3 > dim3/2 ) nz = dim3/2 - thisIndex.z*block3;

#ifdef NAMD_FFTW
  CmiLock(ComputePmeMgr::fftw_plan_lock);

  data = (float *) fftwf_malloc( sizeof(float) * K1*ny*nz*2);
  work = new float[2*K1];

  order_init(initdata.xBlocks);

#ifdef NAMD_FFTW_3
  /* need array of sizes for the how many */
  int fftwFlags = simParams->FFTWPatient ? FFTW_PATIENT  : simParams->FFTWEstimate ? FFTW_ESTIMATE  : FFTW_MEASURE ;
  int sizeLines=ny*nz;
  int planLineSizes[1];
  planLineSizes[0]=K1;
  forward_plan = fftwf_plan_many_dft(1, planLineSizes, sizeLines,
				     (fftwf_complex *) data, NULL, sizeLines, 1,
				     (fftwf_complex *) data, NULL, sizeLines, 1,
				   FFTW_FORWARD,
				     fftwFlags);
  backward_plan = fftwf_plan_many_dft(1, planLineSizes, sizeLines,
				     (fftwf_complex *) data, NULL, sizeLines, 1,
				     (fftwf_complex *) data, NULL, sizeLines, 1,
					  FFTW_BACKWARD,
				      fftwFlags);

#if     CMK_SMP && USE_CKLOOP
  if(simParams->useCkLoop) {
	  //How many FFT plans to be created? The grain-size issue!!.
	  //Currently, I am choosing the min(nx, ny) to be coarse-grain
	  numPlans = (ny<=nz?ny:nz);
          // limit attempted parallelism due to false sharing
          //if ( numPlans < CkMyNodeSize() ) numPlans = (ny>=nz?ny:nz);
          //if ( numPlans < CkMyNodeSize() ) numPlans = sizeLines;
          if ( sizeLines/numPlans < 4 ) numPlans = 1;
	  int howmany = sizeLines/numPlans;
	  forward_plans = new fftwf_plan[numPlans];
	  backward_plans = new fftwf_plan[numPlans];
	  for(int i=0; i<numPlans; i++) {
		  int curStride = i*howmany;		  
		  forward_plans[i] = fftwf_plan_many_dft(1, planLineSizes, howmany,
													 ((fftwf_complex *)data)+curStride, NULL, sizeLines, 1,
													 ((fftwf_complex *)data)+curStride, NULL, sizeLines, 1,
													FFTW_FORWARD,
													 fftwFlags);

		  backward_plans[i] = fftwf_plan_many_dft(1, planLineSizes, howmany,
													 ((fftwf_complex *)data)+curStride, NULL, sizeLines, 1,
													 ((fftwf_complex *)data)+curStride, NULL, sizeLines, 1,
													  FFTW_BACKWARD,
													 fftwFlags);
	  }
  }else
#endif
  {
	  forward_plans = NULL;
	  backward_plans = NULL;
  }
#else
  forward_plan = fftw_create_plan_specific(K1, FFTW_FORWARD,
	( simParams->FFTWEstimate ? FFTW_ESTIMATE : FFTW_MEASURE )
	| FFTW_IN_PLACE | FFTW_USE_WISDOM, (fftw_complex *) data,
	ny*nz, (fftw_complex *) work, 1);
  backward_plan = fftw_create_plan_specific(K1, FFTW_BACKWARD,
	( simParams->FFTWEstimate ? FFTW_ESTIMATE : FFTW_MEASURE )
	| FFTW_IN_PLACE | FFTW_USE_WISDOM, (fftw_complex *) data,
	ny*nz, (fftw_complex *) work, 1);
#endif
  CmiUnlock(ComputePmeMgr::fftw_plan_lock);
#else
  NAMD_die("Sorry, FFTW must be compiled in to use PME.");
#endif

  myKSpace = new PmeKSpace(initdata.grid,
		thisIndex.y*block2, thisIndex.y*block2 + ny,
		thisIndex.z*block3, thisIndex.z*block3 + nz);

}

// #define FFTCHECK   // run a grid of integers through the fft
// #define ZEROCHECK  // check for suspicious zeros in fft

void PmeZPencil::recv_grid(const PmeGridMsg *msg) {

  int dim3 = initdata.grid.dim3;
  if ( imsg == 0 ) {
    lattice = msg->lattice;
    sequence = msg->sequence;
#if ! USE_NODE_PAR_RECEIVE
    memset(data, 0, sizeof(float)*nx*ny*dim3);
#endif
  }

  if ( ! msg->hasData ) return;

  int zlistlen = msg->zlistlen;
#ifdef NAMD_KNL
  int * __restrict msg_zlist = msg->zlist;
  int * __restrict zlist = work_zlist.begin();
  __assume_aligned(zlist,64);
  for ( int k=0; k<zlistlen; ++k ) {
    zlist[k] = msg_zlist[k];
  }
#else
  int * __restrict zlist = msg->zlist;
#endif
  char * __restrict fmsg = msg->fgrid;
  float * __restrict qmsg = msg->qgrid;
  float * __restrict d = data;
  int numGrids = 1;  // pencil FFT doesn't support multiple grids
  for ( int g=0; g<numGrids; ++g ) {
    for ( int i=0; i<nx; ++i ) {
     for ( int j=0; j<ny; ++j, d += dim3 ) {
      if( *(fmsg++) ) {
        #pragma ivdep
        for ( int k=0; k<zlistlen; ++k ) {
          d[zlist[k]] += *(qmsg++);
        }
      }
     }
    }
  }
}

static inline void PmeXZPencilFFT(int first, int last, void *result, int paraNum, void *param){
#ifdef NAMD_FFTW
#ifdef NAMD_FFTW_3    
    fftwf_plan *plans = (fftwf_plan *)param;
    for(int i=first; i<=last; i++) fftwf_execute(plans[i]);
#endif
#endif        
}

void PmeZPencil::forward_fft() {
  evir = 0.;
#ifdef FFTCHECK
  int dim3 = initdata.grid.dim3;
  int K3 = initdata.grid.K3;
  float std_base = 100. * (thisIndex.x+1.) + 10. * (thisIndex.y+1.);
  float *d = data;
  for ( int i=0; i<nx; ++i ) {
   for ( int j=0; j<ny; ++j, d += dim3 ) {
    for ( int k=0; k<dim3; ++k ) {
      d[k] = 10. * (10. * (10. * std_base + i) + j) + k;
    }
   }
  }
#endif
#ifdef NAMD_FFTW
#ifdef MANUAL_DEBUG_FFTW3
  dumpMatrixFloat3("fw_z_b", data, nx, ny, initdata.grid.dim3, thisIndex.x, thisIndex.y, thisIndex.z);
#endif
#ifdef NAMD_FFTW_3
#if     CMK_SMP && USE_CKLOOP
  int useCkLoop = Node::Object()->simParameters->useCkLoop;
  if(useCkLoop>=CKLOOP_CTRL_PME_FORWARDFFT
     && CkNumPes() >= 2 * initdata.xBlocks * initdata.yBlocks) {
          //for(int i=0; i<numPlans; i++) fftwf_execute(forward_plans[i]);
          //transform the above loop
          CkLoop_Parallelize(PmeXZPencilFFT, 1, (void *)forward_plans, CkMyNodeSize(), 0, numPlans-1); //sync
          return;
  }
#endif
  fftwf_execute(forward_plan);
#else
  rfftwnd_real_to_complex(forward_plan, nx*ny,
	data, 1, initdata.grid.dim3, (fftw_complex *) work, 1, 0);
#endif
#ifdef MANUAL_DEBUG_FFTW3
  dumpMatrixFloat3("fw_z_a", data, nx, ny, initdata.grid.dim3, thisIndex.x, thisIndex.y, thisIndex.z);
#endif

#endif
#ifdef ZEROCHECK
  int dim3 = initdata.grid.dim3;
  int K3 = initdata.grid.K3;
  float *d = data;
  for ( int i=0; i<nx; ++i ) {
   for ( int j=0; j<ny; ++j, d += dim3 ) {
    for ( int k=0; k<dim3; ++k ) {
      if ( d[k] == 0. ) CkPrintf("0 in Z at %d %d %d %d %d %d %d %d %d\n",
	thisIndex.x, thisIndex.y, i, j, k, nx, ny, dim3);
    }
   }
  }
#endif
}

/* A single task for partitioned PmeZPencil::send_trans work */
static inline void PmeZPencilSendTrans(int first, int last, void *result, int paraNum, void *param){
	PmeZPencil *zpencil = (PmeZPencil *)param;
	zpencil->send_subset_trans(first, last);	
}

void PmeZPencil::send_subset_trans(int fromIdx, int toIdx){
	int zBlocks = initdata.zBlocks;
	int block3 = initdata.grid.block3;
	int dim3 = initdata.grid.dim3;
	for ( int isend=fromIdx; isend<=toIdx; ++isend ) {
	  int kb = send_order[isend];
	  int nz = block3;
	  if ( (kb+1)*block3 > dim3/2 ) nz = dim3/2 - kb*block3;
	  int hd = ( hasData ? 1 : 0 );
	  PmeTransMsg *msg = new (hd*nx*ny*nz*2,PRIORITY_SIZE) PmeTransMsg;
	  msg->lattice = lattice;
	  msg->sourceNode = thisIndex.y;
	  msg->hasData = hasData;
	  msg->nx = ny;
	 if ( hasData ) {
	  float *md = msg->qgrid;
	  const float *d = data;
	  for ( int i=0; i<nx; ++i ) {
	   for ( int j=0; j<ny; ++j, d += dim3 ) {
		for ( int k=kb*block3; k<(kb*block3+nz); ++k ) {
		  *(md++) = d[2*k];
		  *(md++) = d[2*k+1];
		}
	   }
	  }
	 }
	  msg->sequence = sequence;
	  SET_PRIORITY(msg,sequence,PME_TRANS_PRIORITY)

    CmiEnableUrgentSend(1);
#if USE_NODE_PAR_RECEIVE
      msg->destElem=CkArrayIndex3D(thisIndex.x,0,kb);
#if Y_PERSIST 
      CmiUsePersistentHandle(&trans_handle[isend], 1);
#endif
      initdata.pmeNodeProxy[CmiNodeOf(initdata.ym.ckLocalBranch()->procNum(0,msg->destElem))].recvYTrans(msg);
#if Y_PERSIST 
      CmiUsePersistentHandle(NULL, 0);
#endif    
#else
#if Y_PERSIST 
      CmiUsePersistentHandle(&trans_handle[isend], 1);
#endif
      initdata.yPencil(thisIndex.x,0,kb).recvTrans(msg);
#if Y_PERSIST 
      CmiUsePersistentHandle(NULL, 0);
#endif    
#endif
    CmiEnableUrgentSend(0);
    }
}

void PmeZPencil::send_trans() {
#if USE_PERSISTENT
    if (trans_handle == NULL) setup_persistent();
#endif
#if     CMK_SMP && USE_CKLOOP
	Bool useCkLoop = Node::Object()->simParameters->useCkLoop;
	if(useCkLoop>=CKLOOP_CTRL_PME_SENDTRANS
           && CkNumPes() >= 2 * initdata.xBlocks * initdata.yBlocks) {
		/**
		 * Basically, this function call could be converted into 
		 * a for-loop of: 
		 * for(int i=0; i<=initdata.zBlocks-1; i++) 
		 * send_subset_trans(i,i); 
		 */
		//send_subset_trans(0, initdata.zBlocks-1);
		CkLoop_Parallelize(PmeZPencilSendTrans, 1, (void *)this, CkMyNodeSize(), 0, initdata.zBlocks-1, 1); //not sync
		return;
	}
#endif
  int zBlocks = initdata.zBlocks;
  int block3 = initdata.grid.block3;
  int dim3 = initdata.grid.dim3;
  for ( int isend=0; isend<zBlocks; ++isend ) {
    int kb = send_order[isend];
    int nz = block3;
    if ( (kb+1)*block3 > dim3/2 ) nz = dim3/2 - kb*block3;
    int hd = ( hasData ? 1 : 0 );
    PmeTransMsg *msg = new (hd*nx*ny*nz*2,PRIORITY_SIZE) PmeTransMsg;
    msg->lattice = lattice;
    msg->sourceNode = thisIndex.y;
    msg->hasData = hasData;
    msg->nx = ny;
   if ( hasData ) {
    float *md = msg->qgrid;
    const float *d = data;
    for ( int i=0; i<nx; ++i ) {
     for ( int j=0; j<ny; ++j, d += dim3 ) {
      for ( int k=kb*block3; k<(kb*block3+nz); ++k ) {
        *(md++) = d[2*k];
        *(md++) = d[2*k+1];
      }
     }
    }
   }
    msg->sequence = sequence;
    SET_PRIORITY(msg,sequence,PME_TRANS_PRIORITY)

    CmiEnableUrgentSend(1);
#if USE_NODE_PAR_RECEIVE
    msg->destElem=CkArrayIndex3D(thisIndex.x,0,kb);
#if Y_PERSIST 
    CmiUsePersistentHandle(&trans_handle[isend], 1);
#endif
    initdata.pmeNodeProxy[CmiNodeOf(initdata.ym.ckLocalBranch()->procNum(0,msg->destElem))].recvYTrans(msg);
#if Y_PERSIST 
    CmiUsePersistentHandle(NULL, 0);
#endif    
#else
#if Y_PERSIST 
    CmiUsePersistentHandle(&trans_handle[isend], 1);
#endif
    initdata.yPencil(thisIndex.x,0,kb).recvTrans(msg);
#if Y_PERSIST 
    CmiUsePersistentHandle(NULL, 0);
#endif    
#endif
    CmiEnableUrgentSend(0);
  }
}

void PmeYPencil::recv_trans(const PmeTransMsg *msg) {
  if ( imsg == 0 ) {
    lattice = msg->lattice;
    sequence = msg->sequence;
  }
  int block2 = initdata.grid.block2;
  int K2 = initdata.grid.K2;
  int jb = msg->sourceNode;
  int ny = msg->nx;
 if ( msg->hasData ) {
  const float *md = msg->qgrid;
  float *d = data;
  for ( int i=0; i<nx; ++i, d += K2*nz*2 ) {
   for ( int j=jb*block2; j<(jb*block2+ny); ++j ) {
    for ( int k=0; k<nz; ++k ) {
#ifdef ZEROCHECK
      if ( (*md) == 0. ) CkPrintf("0 in ZY at %d %d %d %d %d %d %d %d %d\n",
	thisIndex.x, jb, thisIndex.z, i, j, k, nx, ny, nz);
#endif
      d[2*(j*nz+k)] = *(md++);
      d[2*(j*nz+k)+1] = *(md++);
    }
   }
  }
 } else {
  float *d = data;
  for ( int i=0; i<nx; ++i, d += K2*nz*2 ) {
   for ( int j=jb*block2; j<(jb*block2+ny); ++j ) {
    for ( int k=0; k<nz; ++k ) {
      d[2*(j*nz+k)] = 0;
      d[2*(j*nz+k)+1] = 0;
    }
   }
  }
 }
}

static inline void PmeYPencilForwardFFT(int first, int last, void *result, int paraNum, void *param){
        PmeYPencil *ypencil = (PmeYPencil *)param;
        ypencil->forward_subset_fft(first, last);
}
void PmeYPencil::forward_subset_fft(int fromIdx, int toIdx) {
#ifdef NAMD_FFTW
#ifdef NAMD_FFTW_3
	for(int i=fromIdx; i<=toIdx; i++){
		fftwf_execute_dft(forward_plan, ((fftwf_complex *) data) + i 
		      * nz * initdata.grid.K2, 	
		      ((fftwf_complex *) data) + i * nz * initdata.grid.K2);
	}
#endif
#endif
}

void PmeYPencil::forward_fft() {
    evir = 0.;
#ifdef NAMD_FFTW
#ifdef MANUAL_DEBUG_FFTW3
  dumpMatrixFloat3("fw_y_b", data, nx, initdata.grid.K2, nz, thisIndex.x, thisIndex.y, thisIndex.z);
#endif
  
#ifdef NAMD_FFTW_3
#if     CMK_SMP && USE_CKLOOP
  int useCkLoop = Node::Object()->simParameters->useCkLoop;
  if(useCkLoop>=CKLOOP_CTRL_PME_FORWARDFFT
     && CkNumPes() >= 2 * initdata.xBlocks * initdata.zBlocks) {
	  CkLoop_Parallelize(PmeYPencilForwardFFT, 1, (void *)this, CkMyNodeSize(), 0, nx-1); //sync
	  return;
  }
#endif
  //the above is a transformation of the following loop using CkLoop
  for ( int i=0; i<nx; ++i ) {
    fftwf_execute_dft(forward_plan, ((fftwf_complex *) data) + i 
		      * nz * initdata.grid.K2, 	
		      ((fftwf_complex *) data) + i * nz * initdata.grid.K2);
  }
#else
  for ( int i=0; i<nx; ++i ) {
    fftw(forward_plan, nz,
	((fftw_complex *) data) + i * nz * initdata.grid.K2,
	nz, 1, (fftw_complex *) work, 1, 0);
  }
#endif
#ifdef MANUAL_DEBUG_FFTW3
  dumpMatrixFloat3("fw_y_a", data, nx, initdata.grid.dim2, nz, thisIndex.x, thisIndex.y, thisIndex.z);
#endif

#endif
}

static inline void PmeYPencilSendTrans(int first, int last, void *result, int paraNum, void *param){
	PmeYPencil *ypencil = (PmeYPencil *)param;
	ypencil->send_subset_trans(first, last);
}

void PmeYPencil::send_subset_trans(int fromIdx, int toIdx){
	int yBlocks = initdata.yBlocks;
	int block2 = initdata.grid.block2;
	int K2 = initdata.grid.K2;
    for ( int isend=fromIdx; isend<=toIdx; ++isend ) {
	  int jb = send_order[isend];
	  int ny = block2;
	  if ( (jb+1)*block2 > K2 ) ny = K2 - jb*block2;
	  int hd = ( hasData ? 1 : 0 );
	  PmeTransMsg *msg = new (hd*nx*ny*nz*2,PRIORITY_SIZE) PmeTransMsg;
	  msg->lattice = lattice;
	  msg->sourceNode = thisIndex.x;
	  msg->hasData = hasData;
	  msg->nx = nx;
	 if ( hasData ) {
	  float *md = msg->qgrid;
	  const float *d = data;
	  for ( int i=0; i<nx; ++i, d += K2*nz*2 ) {
	   for ( int j=jb*block2; j<(jb*block2+ny); ++j ) {
		for ( int k=0; k<nz; ++k ) {
		  *(md++) = d[2*(j*nz+k)];
		  *(md++) = d[2*(j*nz+k)+1];
  #ifdef ZEROCHECK
		  if ( *(md-2) == 0. ) CkPrintf("send 0 in YX at %d %d %d %d %d %d %d %d %d\n",
	  thisIndex.x, jb, thisIndex.z, i, j, k, nx, ny, nz);
  #endif
		}
	   }
	  }
	  if ( md != msg->qgrid + nx*ny*nz*2 ) CkPrintf("error in YX at %d %d %d\n",
	  thisIndex.x, jb, thisIndex.z);
	 }
	  msg->sequence = sequence;
	  SET_PRIORITY(msg,sequence,PME_TRANS2_PRIORITY)
      CmiEnableUrgentSend(1);
#if USE_NODE_PAR_RECEIVE
      msg->destElem=CkArrayIndex3D(0,jb,thisIndex.z);
#if X_PERSIST 
      CmiUsePersistentHandle(&trans_handle[isend], 1);
#endif
      initdata.pmeNodeProxy[CmiNodeOf(initdata.xm.ckLocalBranch()->procNum(0,msg->destElem))].recvXTrans(msg);   
#if X_PERSIST 
      CmiUsePersistentHandle(NULL, 0);
#endif
#else      
#if X_PERSIST 
      CmiUsePersistentHandle(&trans_handle[isend], 1);
#endif
      initdata.xPencil(0,jb,thisIndex.z).recvTrans(msg);
#if X_PERSIST 
      CmiUsePersistentHandle(NULL, 0);
#endif
#endif
      CmiEnableUrgentSend(0);
	}
}

void PmeYPencil::send_trans() {
#if USE_PERSISTENT
    if (trans_handle == NULL) setup_persistent();
#endif
#if     CMK_SMP && USE_CKLOOP
	Bool useCkLoop = Node::Object()->simParameters->useCkLoop;
	if(useCkLoop>=CKLOOP_CTRL_PME_SENDTRANS
           && CkNumPes() >= 2 * initdata.xBlocks * initdata.zBlocks) {
		/**
		 * Basically, this function call could be converted into 
		 * a for-loop of: 
		 * for(int i=0; i<=initdata.yBlocks; i++) 
		 * send_subset_trans(i,i); 
		 */
		//send_subset_trans(0, initdata.yBlocks-1);
		CkLoop_Parallelize(PmeYPencilSendTrans, 1, (void *)this, CkMyNodeSize(), 0, initdata.yBlocks-1, 1); //not sync
		return;
	}
#endif
  int yBlocks = initdata.yBlocks;
  int block2 = initdata.grid.block2;
  int K2 = initdata.grid.K2;
  for ( int isend=0; isend<yBlocks; ++isend ) {
    int jb = send_order[isend];
    int ny = block2;
    if ( (jb+1)*block2 > K2 ) ny = K2 - jb*block2;
    int hd = ( hasData ? 1 : 0 );
    PmeTransMsg *msg = new (hd*nx*ny*nz*2,PRIORITY_SIZE) PmeTransMsg;
    msg->lattice = lattice;
    msg->sourceNode = thisIndex.x;
    msg->hasData = hasData;
    msg->nx = nx;
   if ( hasData ) {
    float *md = msg->qgrid;
    const float *d = data;
    for ( int i=0; i<nx; ++i, d += K2*nz*2 ) {
     for ( int j=jb*block2; j<(jb*block2+ny); ++j ) {
      for ( int k=0; k<nz; ++k ) {
        *(md++) = d[2*(j*nz+k)];
        *(md++) = d[2*(j*nz+k)+1];
#ifdef ZEROCHECK
        if ( *(md-2) == 0. ) CkPrintf("send 0 in YX at %d %d %d %d %d %d %d %d %d\n",
	thisIndex.x, jb, thisIndex.z, i, j, k, nx, ny, nz);
#endif
      }
     }
    }
    if ( md != msg->qgrid + nx*ny*nz*2 ) CkPrintf("error in YX at %d %d %d\n",
	thisIndex.x, jb, thisIndex.z);
   }
    msg->sequence = sequence;
    SET_PRIORITY(msg,sequence,PME_TRANS2_PRIORITY)
    CmiEnableUrgentSend(1);
#if USE_NODE_PAR_RECEIVE
    msg->destElem=CkArrayIndex3D(0,jb,thisIndex.z);
#if X_PERSIST 
        CmiUsePersistentHandle(&trans_handle[isend], 1);
#endif
    initdata.pmeNodeProxy[CmiNodeOf(initdata.xm.ckLocalBranch()->procNum(0,msg->destElem))].recvXTrans(msg);   
#if X_PERSIST 
        CmiUsePersistentHandle(NULL, 0);
#endif
#else
#if X_PERSIST 
        CmiUsePersistentHandle(&trans_handle[isend], 1);
#endif
    initdata.xPencil(0,jb,thisIndex.z).recvTrans(msg);
#if X_PERSIST 
        CmiUsePersistentHandle(NULL, 0);
#endif
    
#endif
    CmiEnableUrgentSend(0);
  }
}

void PmeXPencil::node_process_trans(PmeTransMsg *msg)
{
  if(msg->hasData) hasData=1;
  needs_reply[msg->sourceNode] = msg->hasData;
  recv_trans(msg);
  int limsg;
  CmiMemoryAtomicFetchAndInc(imsg,limsg);
  if(limsg+1 == initdata.xBlocks)
    {
      if(hasData){
        forward_fft();
        pme_kspace();
        backward_fft();
      }
      send_untrans();
      imsg=0;
      CmiMemoryWriteFence();
    }
}

void PmeXPencil::recv_trans(const PmeTransMsg *msg) {
  if ( imsg == 0 ) {
    lattice = msg->lattice;
    sequence = msg->sequence;
  }
  int block1 = initdata.grid.block1;
  int K1 = initdata.grid.K1;
  int ib = msg->sourceNode;
  int nx = msg->nx;
 if ( msg->hasData ) {
  const float *md = msg->qgrid;
  for ( int i=ib*block1; i<(ib*block1+nx); ++i ) {
   float *d = data + i*ny*nz*2;
   for ( int j=0; j<ny; ++j, d += nz*2 ) {
    for ( int k=0; k<nz; ++k ) {
#ifdef ZEROCHECK
      if ( (*md) == 0. ) CkPrintf("0 in YX at %d %d %d %d %d %d %d %d %d\n",
	ib, thisIndex.y, thisIndex.z, i, j, k, nx, ny, nz);
#endif
      d[2*k] = *(md++);
      d[2*k+1] = *(md++);
    }
   }
  }
 } else {
  for ( int i=ib*block1; i<(ib*block1+nx); ++i ) {
   float *d = data + i*ny*nz*2;
   for ( int j=0; j<ny; ++j, d += nz*2 ) {
    for ( int k=0; k<nz; ++k ) {
      d[2*k] = 0;
      d[2*k+1] = 0;
    }
   }
  }
 }
}

void PmeXPencil::forward_fft() {
#ifdef NAMD_FFTW

#ifdef MANUAL_DEBUG_FFTW3
  dumpMatrixFloat3("fw_x_b", data, initdata.grid.K1, ny, nz, thisIndex.x, thisIndex.y, thisIndex.z);
#endif

#ifdef NAMD_FFTW_3
#if     CMK_SMP && USE_CKLOOP
  int useCkLoop = Node::Object()->simParameters->useCkLoop;
  if(useCkLoop>=CKLOOP_CTRL_PME_FORWARDFFT
     && CkNumPes() >= 2 * initdata.yBlocks * initdata.zBlocks) {
	  //for(int i=0; i<numPlans; i++) fftwf_execute(forward_plans[i]);
	  //transform the above loop
	  CkLoop_Parallelize(PmeXZPencilFFT, 1, (void *)forward_plans, CkMyNodeSize(), 0, numPlans-1); //sync
	  return;
  }
#endif
  fftwf_execute(forward_plan);
#else
  fftw(forward_plan, ny*nz,
	((fftw_complex *) data), ny*nz, 1, (fftw_complex *) work, 1, 0);
#endif
#ifdef MANUAL_DEBUG_FFTW3
  dumpMatrixFloat3("fw_x_a", data, initdata.grid.K1, ny, nz, thisIndex.x, thisIndex.y, thisIndex.z);
#endif

#endif
}

void PmeXPencil::pme_kspace() {

  evir = 0.;

#ifdef FFTCHECK
  return;
#endif

  BigReal ewaldcof = ComputeNonbondedUtil::ewaldcof;

  int useCkLoop = 0;
#if CMK_SMP && USE_CKLOOP
  if ( Node::Object()->simParameters->useCkLoop >= CKLOOP_CTRL_PME_KSPACE
       && CkNumPes() >= 2 * initdata.yBlocks * initdata.zBlocks ) {
    useCkLoop = 1;
  }
#endif

  int numGrids = 1;
  for ( int g=0; g<numGrids; ++g ) {
    evir[0] = myKSpace->compute_energy(data+0*g,
		lattice, ewaldcof, &(evir[1]), useCkLoop);
  }
  
#if USE_NODE_PAR_RECEIVE
    CmiMemoryWriteFence();
#endif
}

void PmeXPencil::backward_fft() {
#ifdef NAMD_FFTW
#ifdef MANUAL_DEBUG_FFTW3
  dumpMatrixFloat3("bw_x_b", data, initdata.grid.K1, ny, nz, thisIndex.x, thisIndex.y, thisIndex.z);
#endif

#ifdef NAMD_FFTW_3
#if     CMK_SMP && USE_CKLOOP
  int useCkLoop = Node::Object()->simParameters->useCkLoop;
  if(useCkLoop>=CKLOOP_CTRL_PME_BACKWARDFFT
     && CkNumPes() >= 2 * initdata.yBlocks * initdata.zBlocks) {
          //for(int i=0; i<numPlans; i++) fftwf_execute(backward_plans[i]);
          //transform the above loop
          CkLoop_Parallelize(PmeXZPencilFFT, 1, (void *)backward_plans, CkMyNodeSize(), 0, numPlans-1); //sync
          return;
  }
#endif
  fftwf_execute(backward_plan);
#else
  fftw(backward_plan, ny*nz,
	((fftw_complex *) data), ny*nz, 1, (fftw_complex *) work, 1, 0);
#endif
#ifdef MANUAL_DEBUG_FFTW3
  dumpMatrixFloat3("bw_x_a", data, initdata.grid.K1, ny, nz, thisIndex.x, thisIndex.y, thisIndex.z);
#endif
#endif
}

static inline void PmeXPencilSendUntrans(int first, int last, void *result, int paraNum, void *param){
	int evirIdx = paraNum;
	PmeXPencil *xpencil = (PmeXPencil *)param;
	xpencil->send_subset_untrans(first, last, evirIdx);
}

void PmeXPencil::send_subset_untrans(int fromIdx, int toIdx, int evirIdx){
	int xBlocks = initdata.xBlocks;
	int block1 = initdata.grid.block1;	
	int K1 = initdata.grid.K1;

	int ackL=0, ackH=-1;
	int unL=0, unH=-1;
	int send_evir=0;
	if(fromIdx >= evirIdx+1) {
		//send PmeUntransMsg with has_evir=0
		unL = fromIdx;
		unH = toIdx;		
	} else if(toIdx <= evirIdx-1) {
		//send PmeAckMsg
		ackL=fromIdx;
		ackH=toIdx;		
	} else {
		//partially send PmeAckMsg and partially send PmeUntransMsg
		ackL=fromIdx;
		ackH=evirIdx-1;
		send_evir=1;
		unL=evirIdx+1;
		unH=toIdx;
	}

	for(int isend=ackL; isend<=ackH; isend++) {
		//send PmeAckMsg
        CmiEnableUrgentSend(1);
		int ib = send_order[isend];
		PmeAckMsg *msg = new (PRIORITY_SIZE) PmeAckMsg;
		SET_PRIORITY(msg,sequence,PME_UNTRANS_PRIORITY)
		initdata.yPencil(ib,0,thisIndex.z).recvAck(msg);
        CmiEnableUrgentSend(0);
    }

    CmiEnableUrgentSend(1);
	//send PmeUntransMsg with has_evir=1
	if(send_evir) {
		int ib = send_order[evirIdx];
		int nx = block1;
		if ( (ib+1)*block1 > K1 ) nx = K1 - ib*block1;
		PmeUntransMsg *msg = new (nx*ny*nz*2,PRIORITY_SIZE) PmeUntransMsg;		
		msg->sourceNode = thisIndex.y;
		msg->ny = ny;
		float *md = msg->qgrid;
		for ( int i=ib*block1; i<(ib*block1+nx); ++i ) {
			float *d = data + i*ny*nz*2;
			for ( int j=0; j<ny; ++j, d += nz*2 ) {
				for ( int k=0; k<nz; ++k ) {
					*(md++) = d[2*k];
					*(md++) = d[2*k+1];
				}
			}
		}
		SET_PRIORITY(msg,sequence,PME_UNTRANS_PRIORITY)
#if USE_NODE_PAR_RECEIVE
        msg->destElem=CkArrayIndex3D(ib,0, thisIndex.z);
        initdata.pmeNodeProxy[CmiNodeOf(initdata.ym.ckLocalBranch()->procNum(0,msg->destElem))].recvYUntrans(msg);
#else
        initdata.yPencil(ib,0,thisIndex.z).recvUntrans(msg);
#endif
	 }
    CmiEnableUrgentSend(0);
	
	//send PmeUntransMsg with has_evir=0
	for(int isend=unL; isend<=unH; isend++) {
		int ib = send_order[isend];
		int nx = block1;
		if ( (ib+1)*block1 > K1 ) nx = K1 - ib*block1;
		PmeUntransMsg *msg = new (nx*ny*nz*2,PRIORITY_SIZE) PmeUntransMsg;
		msg->sourceNode = thisIndex.y;
		msg->ny = ny;
		float *md = msg->qgrid;
		for ( int i=ib*block1; i<(ib*block1+nx); ++i ) {
			float *d = data + i*ny*nz*2;
			for ( int j=0; j<ny; ++j, d += nz*2 ) {
				for ( int k=0; k<nz; ++k ) {
					*(md++) = d[2*k];
					*(md++) = d[2*k+1];
				}
			}
		}
		SET_PRIORITY(msg,sequence,PME_UNTRANS_PRIORITY)
        CmiEnableUrgentSend(1);
#if USE_NODE_PAR_RECEIVE
        msg->destElem=CkArrayIndex3D(ib,0, thisIndex.z);
#if Y_PERSIST 
        CmiUsePersistentHandle(&untrans_handle[isend], 1);
#endif
        initdata.pmeNodeProxy[CmiNodeOf(initdata.ym.ckLocalBranch()->procNum(0,msg->destElem))].recvYUntrans(msg);
#if Y_PERSIST 
        CmiUsePersistentHandle(NULL, 0);
#endif
#else
#if Y_PERSIST 
  //      CmiUsePersistentHandle(&untrans_handle[isend], 1);
#endif
        initdata.yPencil(ib,0,thisIndex.z).recvUntrans(msg);
#if Y_PERSIST 
   //     CmiUsePersistentHandle(NULL, 0);
#endif
#endif
        CmiEnableUrgentSend(0);
	}
}

void PmeXPencil::send_untrans() {

  { // send energy and virial
    int numGrids = 1;
    PmeEvirMsg *newmsg = new (numGrids, PRIORITY_SIZE) PmeEvirMsg;
    newmsg->evir[0] = evir;
    SET_PRIORITY(newmsg,sequence,PME_UNGRID_PRIORITY)
    CmiEnableUrgentSend(1);
    initdata.pmeProxy[recipEvirPe].recvRecipEvir(newmsg);
    CmiEnableUrgentSend(0);
  }

#if USE_PERSISTENT
  if (untrans_handle == NULL) setup_persistent();
#endif
#if     CMK_SMP && USE_CKLOOP
  Bool useCkLoop = Node::Object()->simParameters->useCkLoop;
  if(useCkLoop>=CKLOOP_CTRL_PME_SENDUNTRANS
     && CkNumPes() >= 2 * initdata.yBlocks * initdata.zBlocks) {
	  	int xBlocks = initdata.xBlocks;
		int evirIdx = 0;
		for ( int isend=0; isend<xBlocks; ++isend ) {
			int ib = send_order[isend];
			if (needs_reply[ib]) {
				evirIdx = isend;
				break;
			}
		}

		//basically: 
		//[0,evirIdx-1]->send PmeAckMsg
		//evirIdx->send PmeUntransMsg with has_evir=1
		//[evirIdx+1, xBlocks-1]->send PmeUntransMsg with has_evir=0
		//send_subset_untrans(0, xBlocks-1, evirIdx);
#if USE_NODE_PAR_RECEIVE
		//CkLoop_Parallelize(PmeXPencilSendUntrans, evirIdx, (void *)this, CkMyNodeSize(), 0, xBlocks-1, 1); //has to sync
		CkLoop_Parallelize(PmeXPencilSendUntrans, evirIdx, (void *)this, xBlocks, 0, xBlocks-1, 1); //has to sync
#else
        //CkLoop_Parallelize(PmeXPencilSendUntrans, evirIdx, (void *)this, CkMyNodeSize(), 0, xBlocks-1, 0); //not sync
		CkLoop_Parallelize(PmeXPencilSendUntrans, evirIdx, (void *)this, xBlocks, 0, xBlocks-1, 0); //not sync
#endif        
		return;
  }
#endif
  int xBlocks = initdata.xBlocks;
  int block1 = initdata.grid.block1;
  int K1 = initdata.grid.K1;
  int send_evir = 1;
  for ( int isend=0; isend<xBlocks; ++isend ) {
    int ib = send_order[isend];
    if ( ! needs_reply[ib] ) {
      PmeAckMsg *msg = new (PRIORITY_SIZE) PmeAckMsg;
      CmiEnableUrgentSend(1);
      SET_PRIORITY(msg,sequence,PME_UNTRANS_PRIORITY)
      initdata.yPencil(ib,0,thisIndex.z).recvAck(msg);
      CmiEnableUrgentSend(0);
      continue;
    }
    int nx = block1;
    if ( (ib+1)*block1 > K1 ) nx = K1 - ib*block1;
    PmeUntransMsg *msg = new (nx*ny*nz*2,PRIORITY_SIZE) PmeUntransMsg;
    if ( send_evir ) {
      send_evir = 0;
    }
    msg->sourceNode = thisIndex.y;
    msg->ny = ny;
    float *md = msg->qgrid;
    for ( int i=ib*block1; i<(ib*block1+nx); ++i ) {
     float *d = data + i*ny*nz*2;
     for ( int j=0; j<ny; ++j, d += nz*2 ) {
      for ( int k=0; k<nz; ++k ) {
        *(md++) = d[2*k];
        *(md++) = d[2*k+1];
      }
     }
    }
    SET_PRIORITY(msg,sequence,PME_UNTRANS_PRIORITY)

    CmiEnableUrgentSend(1);
#if USE_NODE_PAR_RECEIVE
    msg->destElem=CkArrayIndex3D(ib,0, thisIndex.z);
#if Y_PERSIST 
    CmiUsePersistentHandle(&untrans_handle[isend], 1);
#endif
    initdata.pmeNodeProxy[CmiNodeOf(initdata.ym.ckLocalBranch()->procNum(0,msg->destElem))].recvYUntrans(msg);
#if Y_PERSIST 
    CmiUsePersistentHandle(NULL, 0);
#endif
#else
#if Y_PERSIST 
    CmiUsePersistentHandle(&untrans_handle[isend], 1);
#endif
    initdata.yPencil(ib,0,thisIndex.z).recvUntrans(msg);
#if Y_PERSIST 
    CmiUsePersistentHandle(NULL, 0);
#endif
#endif
    CmiEnableUrgentSend(0);
  }
}

void PmeYPencil::recv_untrans(const PmeUntransMsg *msg) {
  int block2 = initdata.grid.block2;
  int K2 = initdata.grid.K2;
  int jb = msg->sourceNode;
  int ny = msg->ny;
  const float *md = msg->qgrid;
  float *d = data;
  for ( int i=0; i<nx; ++i, d += K2*nz*2 ) {
#if CMK_BLUEGENEL
    CmiNetworkProgress();
#endif   
    for ( int j=jb*block2; j<(jb*block2+ny); ++j ) {
      for ( int k=0; k<nz; ++k ) {
#ifdef ZEROCHECK
	if ( (*md) == 0. ) CkPrintf("0 in XY at %d %d %d %d %d %d %d %d %d\n",
				    thisIndex.x, jb, thisIndex.z, i, j, k, nx, ny, nz);
#endif
	d[2*(j*nz+k)] = *(md++);
	d[2*(j*nz+k)+1] = *(md++);
      }
    }
  }
}

static inline void PmeYPencilBackwardFFT(int first, int last, void *result, int paraNum, void *param){
	PmeYPencil *ypencil = (PmeYPencil *)param;
	ypencil->backward_subset_fft(first, last);
}

void PmeYPencil::backward_subset_fft(int fromIdx, int toIdx) {
#ifdef NAMD_FFTW
#ifdef NAMD_FFTW_3
	for(int i=fromIdx; i<=toIdx; i++){
		fftwf_execute_dft(backward_plan, 	
						  ((fftwf_complex *) data) + i * nz * initdata.grid.K2,    	
						  ((fftwf_complex *) data) + i * nz * initdata.grid.K2);
	}
#endif
#endif
}

void PmeYPencil::backward_fft() {
#ifdef NAMD_FFTW
#ifdef MANUAL_DEBUG_FFTW3
  dumpMatrixFloat3("bw_y_b", data, nx, initdata.grid.K2, nz, thisIndex.x, thisIndex.y, thisIndex.z);
#endif

#ifdef NAMD_FFTW_3
#if     CMK_SMP && USE_CKLOOP
  int useCkLoop = Node::Object()->simParameters->useCkLoop;
  if(useCkLoop>=CKLOOP_CTRL_PME_BACKWARDFFT
     && CkNumPes() >= 2 * initdata.xBlocks * initdata.zBlocks) {
	  CkLoop_Parallelize(PmeYPencilBackwardFFT, 1, (void *)this, CkMyNodeSize(), 0, nx-1); //sync
	  return;
  }
#endif
  //the above is a transformation of the following loop using CkLoop
  for ( int i=0; i<nx; ++i ) {
#if CMK_BLUEGENEL
	CmiNetworkProgress();
#endif
    fftwf_execute_dft(backward_plan, 	
					  ((fftwf_complex *) data) + i * nz * initdata.grid.K2,
					  ((fftwf_complex *) data) + i * nz * initdata.grid.K2);
  }
#else
	for ( int i=0; i<nx; ++i ) {
#if CMK_BLUEGENEL
	  CmiNetworkProgress();
#endif
		fftw(backward_plan, nz,
		((fftw_complex *) data) + i * nz * initdata.grid.K2,
		nz, 1, (fftw_complex *) work, 1, 0);
	}
#endif

#ifdef MANUAL_DEBUG_FFTW3
  dumpMatrixFloat3("bw_y_a", data, nx, initdata.grid.K2, nz, thisIndex.x, thisIndex.y, thisIndex.z);
#endif

#endif
}

static inline void PmeYPencilSendUntrans(int first, int last, void *result, int paraNum, void *param){
        int evirIdx = paraNum;
        PmeYPencil *ypencil = (PmeYPencil *)param;
        ypencil->send_subset_untrans(first, last, evirIdx);
}

void PmeYPencil::send_subset_untrans(int fromIdx, int toIdx, int evirIdx){
	int yBlocks = initdata.yBlocks;
	int block2 = initdata.grid.block2;	
	int K2 = initdata.grid.K2;

	int ackL=0, ackH=-1;
	int unL=0, unH=-1;
	int send_evir=0;
	if(fromIdx >= evirIdx+1) {
		//send PmeUntransMsg with has_evir=0
		unL = fromIdx;
		unH = toIdx;		
	} else if(toIdx <= evirIdx-1) {
		//send PmeAckMsg
		ackL=fromIdx;
		ackH=toIdx;		
	} else {
		//partially send PmeAckMsg and partially send PmeUntransMsg
		ackL=fromIdx;
		ackH=evirIdx-1;
		send_evir=1;
		unL=evirIdx+1;
		unH=toIdx;
	}

	for(int isend=ackL; isend<=ackH; isend++) {
		//send PmeAckMsg
        CmiEnableUrgentSend(1);
		int jb = send_order[isend];
		PmeAckMsg *msg = new (PRIORITY_SIZE) PmeAckMsg;
		SET_PRIORITY(msg,sequence,PME_UNTRANS2_PRIORITY)
		initdata.zPencil(thisIndex.x,jb,0).recvAck(msg);
        CmiEnableUrgentSend(0);
	}

    CmiEnableUrgentSend(1);
	//send PmeUntransMsg with has_evir=1
	if(send_evir) {
		int jb = send_order[evirIdx];
		int ny = block2;
		if ( (jb+1)*block2 > K2 ) ny = K2 - jb*block2;
		PmeUntransMsg *msg = new (nx*ny*nz*2,PRIORITY_SIZE) PmeUntransMsg;		
		msg->sourceNode = thisIndex.z;
		msg->ny = nz;
		float *md = msg->qgrid;
		const float *d = data;
		for ( int i=0; i<nx; ++i, d += K2*nz*2 ) {
			for ( int j=jb*block2; j<(jb*block2+ny); ++j ) {
				for ( int k=0; k<nz; ++k ) {
					*(md++) = d[2*(j*nz+k)];
					*(md++) = d[2*(j*nz+k)+1];
				}
			}
		}
		SET_PRIORITY(msg,sequence,PME_UNTRANS2_PRIORITY)
#if USE_NODE_PAR_RECEIVE
        msg->destElem=CkArrayIndex3D( thisIndex.x, jb, 0);
    //    CkPrintf("[%d] sending to %d %d %d recvZUntrans on node %d\n", CkMyPe(), thisIndex.x, jb, 0, CmiNodeOf(initdata.zm.ckLocalBranch()->procNum(0,msg->destElem)));
        initdata.pmeNodeProxy[CmiNodeOf(initdata.zm.ckLocalBranch()->procNum(0,msg->destElem))].recvZUntrans(msg);
#else
        initdata.zPencil(thisIndex.x,jb,0).recvUntrans(msg);
#endif
	}

    CmiEnableUrgentSend(0);
	//send PmeUntransMsg with has_evir=0
	for(int isend=unL; isend<=unH; isend++) {
		int jb = send_order[isend];
		int ny = block2;
		if ( (jb+1)*block2 > K2 ) ny = K2 - jb*block2;
		PmeUntransMsg *msg = new (nx*ny*nz*2,PRIORITY_SIZE) PmeUntransMsg;
		msg->sourceNode = thisIndex.z;
		msg->ny = nz;
		float *md = msg->qgrid;
		const float *d = data;
		for ( int i=0; i<nx; ++i, d += K2*nz*2 ) {
			for ( int j=jb*block2; j<(jb*block2+ny); ++j ) {
				for ( int k=0; k<nz; ++k ) {
					*(md++) = d[2*(j*nz+k)];
					*(md++) = d[2*(j*nz+k)+1];
				}
			}
		}
		SET_PRIORITY(msg,sequence,PME_UNTRANS2_PRIORITY)
            CmiEnableUrgentSend(1);
#if USE_NODE_PAR_RECEIVE
        msg->destElem=CkArrayIndex3D( thisIndex.x, jb, 0);
        //    CkPrintf("[%d] sending to %d %d %d recvZUntrans on node %d\n", CkMyPe(), thisIndex.x, jb, 0, CmiNodeOf(initdata.zm.ckLocalBranch()->procNum(0,msg->destElem)));
#if Z_PERSIST 
        CmiUsePersistentHandle(&untrans_handle[isend], 1);
#endif
        initdata.pmeNodeProxy[CmiNodeOf(initdata.zm.ckLocalBranch()->procNum(0,msg->destElem))].recvZUntrans(msg);
#if Z_PERSIST 
        CmiUsePersistentHandle(NULL, 0);
#endif
#else
#if Z_PERSIST 
        CmiUsePersistentHandle(&untrans_handle[isend], 1);
#endif
        initdata.zPencil(thisIndex.x,jb,0).recvUntrans(msg);
#if Z_PERSIST 
        CmiUsePersistentHandle(NULL, 0);
#endif
#endif
    CmiEnableUrgentSend(0);
	}
}

void PmeYPencil::send_untrans() {
#if USE_PERSISTENT
  if (untrans_handle == NULL) setup_persistent();
#endif
#if     CMK_SMP && USE_CKLOOP
  Bool useCkLoop = Node::Object()->simParameters->useCkLoop;
  if(useCkLoop>=CKLOOP_CTRL_PME_SENDUNTRANS
     && CkNumPes() >= 2 * initdata.xBlocks * initdata.zBlocks) {
	  int yBlocks = initdata.yBlocks;
	  int evirIdx = 0;
	  for ( int isend=0; isend<yBlocks; ++isend ) {
		  int jb = send_order[isend];
		  if (needs_reply[jb]) {
			  evirIdx = isend;
			  break;
		  }
	  }

	  //basically: 
	  //[0,evirIdx-1]->send PmeAckMsg
	  //evirIdx->send PmeUntransMsg with has_evir=1
	  //[evirIdx+1, yBlocks-1]->send PmeUntransMsg with has_evir=0
	  //send_subset_untrans(0, yBlocks-1, evirIdx);
#if USE_NODE_PAR_RECEIVE      
	  //CkLoop_Parallelize(PmeYPencilSendUntrans, evirIdx, (void *)this, CkMyNodeSize(), 0, yBlocks-1, 1); //sync
	  CkLoop_Parallelize(PmeYPencilSendUntrans, evirIdx, (void *)this, yBlocks, 0, yBlocks-1, 1);
      evir = 0.;
      CmiMemoryWriteFence();
#else
      //CkLoop_Parallelize(PmeYPencilSendUntrans, evirIdx, (void *)this, CkMyNodeSize(), 0, yBlocks-1, 0); //not sync
	  CkLoop_Parallelize(PmeYPencilSendUntrans, evirIdx, (void *)this, yBlocks, 0, yBlocks-1, 0); //not sync
#endif
	  return;
  }
#endif
  int yBlocks = initdata.yBlocks;
  int block2 = initdata.grid.block2;
  int K2 = initdata.grid.K2;
  int send_evir = 1;
  for ( int isend=0; isend<yBlocks; ++isend ) {
    int jb = send_order[isend];
    if ( ! needs_reply[jb] ) {
      PmeAckMsg *msg = new (PRIORITY_SIZE) PmeAckMsg;
      CmiEnableUrgentSend(1);
      SET_PRIORITY(msg,sequence,PME_UNTRANS2_PRIORITY)
      initdata.zPencil(thisIndex.x,jb,0).recvAck(msg);
      CmiEnableUrgentSend(0);
      continue;
    }
    int ny = block2;
    if ( (jb+1)*block2 > K2 ) ny = K2 - jb*block2;
    PmeUntransMsg *msg = new (nx*ny*nz*2,PRIORITY_SIZE) PmeUntransMsg;
    if ( send_evir ) {
      send_evir = 0;
    }
    msg->sourceNode = thisIndex.z;
    msg->ny = nz;
    float *md = msg->qgrid;
    const float *d = data;
    for ( int i=0; i<nx; ++i, d += K2*nz*2 ) {
     for ( int j=jb*block2; j<(jb*block2+ny); ++j ) {
      for ( int k=0; k<nz; ++k ) {
        *(md++) = d[2*(j*nz+k)];
        *(md++) = d[2*(j*nz+k)+1];
      }
     }
    }
    SET_PRIORITY(msg,sequence,PME_UNTRANS2_PRIORITY)

    CmiEnableUrgentSend(1);
#if USE_NODE_PAR_RECEIVE
    msg->destElem=CkArrayIndex3D( thisIndex.x, jb, 0);
    //    CkPrintf("[%d] sending to %d %d %d recvZUntrans on node %d\n", CkMyPe(), thisIndex.x, jb, 0, CmiNodeOf(initdata.zm.ckLocalBranch()->procNum(0,msg->destElem)));
#if Z_PERSIST 
    CmiUsePersistentHandle(&untrans_handle[isend], 1);
#endif
    initdata.pmeNodeProxy[CmiNodeOf(initdata.zm.ckLocalBranch()->procNum(0,msg->destElem))].recvZUntrans(msg);
#if Z_PERSIST
    CmiUsePersistentHandle(NULL, 0);
#endif
#else
#if Z_PERSIST 
    CmiUsePersistentHandle(&untrans_handle[isend], 1);
#endif
    initdata.zPencil(thisIndex.x,jb,0).recvUntrans(msg);
#if Z_PERSIST 
    CmiUsePersistentHandle(NULL, 0);
#endif
#endif    
    CmiEnableUrgentSend(0);
  }
  
#if USE_NODE_PAR_RECEIVE
  evir = 0.;
  CmiMemoryWriteFence();
#endif
}

void PmeZPencil::recv_untrans(const PmeUntransMsg *msg) {
#if ! USE_NODE_PAR_RECEIVE
    if(imsg==0) evir=0.;
#endif

  int block3 = initdata.grid.block3;
  int dim3 = initdata.grid.dim3;
  int kb = msg->sourceNode;
  int nz = msg->ny;
  const float *md = msg->qgrid;
  float *d = data;
  for ( int i=0; i<nx; ++i ) {
#if CMK_BLUEGENEL
    CmiNetworkProgress();
#endif   
    for ( int j=0; j<ny; ++j, d += dim3 ) {
      for ( int k=kb*block3; k<(kb*block3+nz); ++k ) {
#ifdef ZEROCHECK
	if ( (*md) == 0. ) CkPrintf("0 in YZ at %d %d %d %d %d %d %d %d %d\n",
				    thisIndex.x, thisIndex.y, kb, i, j, k, nx, ny, nz);
#endif
	d[2*k] = *(md++);
	d[2*k+1] = *(md++);
      }
    }
  }
}

void PmeZPencil::backward_fft() {
#ifdef NAMD_FFTW
#ifdef MANUAL_DEBUG_FFTW3
  dumpMatrixFloat3("bw_z_b", data, nx, ny, initdata.grid.dim3, thisIndex.x, thisIndex.y, thisIndex.z);
#endif
#ifdef NAMD_FFTW_3
#if     CMK_SMP && USE_CKLOOP
  int useCkLoop = Node::Object()->simParameters->useCkLoop;
  if(useCkLoop>=CKLOOP_CTRL_PME_BACKWARDFFT
     && CkNumPes() >= 2 * initdata.xBlocks * initdata.yBlocks) {
	  //for(int i=0; i<numPlans; i++) fftwf_execute(backward_plans[i]);
	  //transform the above loop
	  CkLoop_Parallelize(PmeXZPencilFFT, 1, (void *)backward_plans, CkMyNodeSize(), 0, numPlans-1); //sync
	  return;
  }
#endif
  fftwf_execute(backward_plan);
#else
  rfftwnd_complex_to_real(backward_plan, nx*ny,
	    (fftw_complex *) data, 1, initdata.grid.dim3/2, work, 1, 0);
#endif
#ifdef MANUAL_DEBUG_FFTW3
  dumpMatrixFloat3("bw_z_a", data, nx, ny, initdata.grid.dim3, thisIndex.x, thisIndex.y, thisIndex.z);
#endif

#endif
  
#if CMK_BLUEGENEL
  CmiNetworkProgress();
#endif

#ifdef FFTCHECK
  int dim3 = initdata.grid.dim3;
  int K1 = initdata.grid.K1;
  int K2 = initdata.grid.K2;
  int K3 = initdata.grid.K3;
  float scale = 1. / (1. * K1 * K2 * K3);
  float maxerr = 0.;
  float maxstd = 0.;
  int mi, mj, mk;  mi = mj = mk = -1;
  float std_base = 100. * (thisIndex.x+1.) + 10. * (thisIndex.y+1.);
  const float *d = data;
  for ( int i=0; i<nx; ++i ) {
   for ( int j=0; j<ny; ++j, d += dim3 ) {
    for ( int k=0; k<K3; ++k ) {
      float std = 10. * (10. * (10. * std_base + i) + j) + k;
      float err = scale * d[k] - std;
      if ( fabsf(err) > fabsf(maxerr) ) {
        maxerr = err;
        maxstd = std;
        mi = i;  mj = j;  mk = k;
      }
    }
   }
  }
  CkPrintf("pencil %d %d max error %f at %d %d %d (should be %f)\n",
		thisIndex.x, thisIndex.y, maxerr, mi, mj, mk, maxstd);
#endif

}

static inline void PmeZPencilSendUngrid(int first, int last, void *result, int paraNum, void *param){
	//to take advantage of the interface which allows 3 user params at most.
	//under such situtation, no new parameter list needs to be created!! -Chao Mei
	int specialIdx = paraNum;
	PmeZPencil *zpencil = (PmeZPencil *)param;
	zpencil->send_subset_ungrid(first, last, specialIdx);
}

void PmeZPencil::send_all_ungrid() {
/* 
//Original code: the transformation is to first extract the msg 
//idx that will has evir value set. -Chao Mei  
	int send_evir = 1;
	for (int imsg=0; imsg < grid_msgs.size(); ++imsg ) {
		PmeGridMsg *msg = grid_msgs[imsg];
		if ( msg->hasData ) {
			if ( send_evir ) {
				msg->evir[0] = evir;
				send_evir = 0;
			} else {
				msg->evir[0] = 0.;
			}
		}
		send_ungrid(msg);
	}
*/
	int evirIdx = 0;
	for(int imsg=0; imsg<grid_msgs.size(); imsg++) {
		if(grid_msgs[imsg]->hasData) {
			evirIdx = imsg;
			break;
		}
	}

#if     CMK_SMP && USE_CKLOOP
	Bool useCkLoop = Node::Object()->simParameters->useCkLoop;
	if(useCkLoop>=CKLOOP_CTRL_PME_SENDUNTRANS
           && CkNumPes() >= 2 * initdata.xBlocks * initdata.yBlocks) {
		//????What's the best value for numChunks?????
#if USE_NODE_PAR_RECEIVE        
		//CkLoop_Parallelize(PmeZPencilSendUngrid, evirIdx, (void *)this, CkMyNodeSize(), 0, grid_msgs.size()-1, 1); //has to sync
		CkLoop_Parallelize(PmeZPencilSendUngrid, evirIdx, (void *)this, grid_msgs.size(), 0, grid_msgs.size()-1, 1); //has to sync
#else
        //CkLoop_Parallelize(PmeZPencilSendUngrid, evirIdx, (void *)this, CkMyNodeSize(), 0, grid_msgs.size()-1, 0); //not sync
		CkLoop_Parallelize(PmeZPencilSendUngrid, evirIdx, (void *)this, grid_msgs.size(), 0, grid_msgs.size()-1, 0); //not sync
#endif        
		return;
	}
#endif
	send_subset_ungrid(0, grid_msgs.size()-1, evirIdx);
}

void PmeZPencil::send_subset_ungrid(int fromIdx, int toIdx, int specialIdx){
	for (int imsg=fromIdx; imsg <=toIdx; ++imsg ) {
		PmeGridMsg *msg = grid_msgs[imsg];
		send_ungrid(msg);
	}
}

void PmeZPencil::send_ungrid(PmeGridMsg *msg) {

#ifdef NAMD_CUDA
  const int UNGRID_PRIORITY = ( offload ? PME_OFFLOAD_UNGRID_PRIORITY : PME_UNGRID_PRIORITY );
#else
  const int UNGRID_PRIORITY = PME_UNGRID_PRIORITY ;
#endif

  int pe = msg->sourceNode;
  if ( ! msg->hasData ) {
    delete msg;
    PmeAckMsg *ackmsg = new (PRIORITY_SIZE) PmeAckMsg;
    SET_PRIORITY(ackmsg,sequence,UNGRID_PRIORITY)
    CmiEnableUrgentSend(1);
    initdata.pmeProxy[pe].recvAck(ackmsg);
    CmiEnableUrgentSend(0);
    return;
  }
  msg->sourceNode = thisIndex.x * initdata.yBlocks + thisIndex.y;
  int dim3 = initdata.grid.dim3;
  int zlistlen = msg->zlistlen;
  int *zlist = msg->zlist;
  char *fmsg = msg->fgrid;
  float *qmsg = msg->qgrid;
  float *d = data;
  int numGrids = 1;  // pencil FFT doesn't support multiple grids
  for ( int g=0; g<numGrids; ++g ) {
#if CMK_BLUEGENEL
    CmiNetworkProgress();
#endif    
    for ( int i=0; i<nx; ++i ) {
      for ( int j=0; j<ny; ++j, d += dim3 ) {
	if( *(fmsg++) ) {
	  for ( int k=0; k<zlistlen; ++k ) {
	    *(qmsg++) = d[zlist[k]];
	  }
	}
      }
    }
  }
  SET_PRIORITY(msg,sequence,UNGRID_PRIORITY)
    CmiEnableUrgentSend(1);
#ifdef NAMD_CUDA
    if ( offload ) {
      initdata.pmeNodeProxy[CkNodeOf(pe)].recvUngrid(msg);
    } else
#endif
  initdata.pmeProxy[pe].recvUngrid(msg);
    CmiEnableUrgentSend(0);
}

void PmeZPencil::node_process_grid(PmeGridMsg *msg)
{
#if USE_NODE_PAR_RECEIVE
  CmiLock(ComputePmeMgr::fftw_plan_lock);
  CmiMemoryReadFence();
#endif
  recv_grid(msg);
  if(msg->hasData) hasData=msg->hasData;
  int limsg;
  CmiMemoryAtomicFetchAndInc(imsg,limsg);
  grid_msgs[limsg] = msg;
  //  CkPrintf("[%d] PmeZPencil node_process_grid for %d %d %d has %d of %d imsg %d\n",CkMyPe(),thisIndex.x,thisIndex.y,thisIndex.z, limsg, grid_msgs.size(), imsg);      
  if(limsg+1 == grid_msgs.size())
    {

      if (hasData)
        {
          forward_fft();
        }
      send_trans();
      imsg=0;
      CmiMemoryWriteFence();
      //      CkPrintf("[%d] PmeZPencil grid node_zero imsg for %d %d %d\n",CkMyPe(),thisIndex.x,thisIndex.y,thisIndex.z);
    }
#if USE_NODE_PAR_RECEIVE
  CmiUnlock(ComputePmeMgr::fftw_plan_lock);
  CmiMemoryWriteFence();
#endif
}

void PmeZPencil::node_process_untrans(PmeUntransMsg *msg)
{
  recv_untrans(msg);
#if USE_NODE_PAR_RECEIVE
  CmiMemoryWriteFence();
  CmiLock(ComputePmeMgr::fftw_plan_lock);
#endif    
  int limsg;
  CmiMemoryAtomicFetchAndInc(imsgb,limsg);
  if(limsg+1 == initdata.zBlocks)
    {
#if USE_NODE_PAR_RECEIVE
      CmiMemoryReadFence();
#endif    
      if(hasData) // maybe this should be an assert
        {
          backward_fft();
        }
        
        send_all_ungrid();
    /*  int send_evir = 1;
      // TODO: this part should use Chao's output parallelization
      for ( limsg=0; limsg < grid_msgs.size(); ++limsg ) {
        PmeGridMsg *omsg = grid_msgs[limsg];
        if ( omsg->hasData ) {
          if ( send_evir ) {
            omsg->evir[0] = evir;
            send_evir = 0;
          } else {
            omsg->evir[0] = 0.;
          }
        }
        send_ungrid(omsg);
      } */
      imsgb=0;
      evir = 0;
      memset(data, 0, sizeof(float) * nx*ny* initdata.grid.dim3); 
      CmiMemoryWriteFence();
      //      CkPrintf("[%d] PmeZPencil untrans node_zero imsg for %d %d %d\n",CkMyPe(),thisIndex.x,thisIndex.y,thisIndex.z);
    }
#if USE_NODE_PAR_RECEIVE
  CmiUnlock(ComputePmeMgr::fftw_plan_lock);
#endif
}


#include "ComputePmeMgr.def.h"

