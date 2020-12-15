
#include "common.h"
#include "charm++.h"

#ifdef NAMD_CUDA
#include <cuda_runtime.h>
#include <cuda.h>
#endif

#include "WorkDistrib.h"
#include "ComputeMgr.h"
#include "ProxyMgr.h"
#include "ComputeNonbondedCUDAKernel.h"
#include "ComputeNonbondedCUDA.h"
#include "LJTable.h"
#include "ObjectArena.h"
#include "SortAtoms.h"
#include "Priorities.h"
#include <algorithm>

#include "NamdTypes.h"
#include "DeviceCUDA.h"
#include "CudaUtils.h"

//#define PRINT_GBIS
#undef PRINT_GBIS

#ifdef NAMD_CUDA

#ifdef WIN32
#define __thread __declspec(thread)
#endif

extern __thread int max_grid_size;

extern __thread cudaStream_t stream;
extern __thread cudaStream_t stream2;

extern __thread DeviceCUDA *deviceCUDA;

void cuda_errcheck(const char *msg) {
  cudaError_t err;
  if ((err = cudaGetLastError()) != cudaSuccess) {
    char host[128];
#ifdef NOHOSTNAME
    sprintf(host,"physical node %d", CmiPhysicalNodeID(CkMyPe()));
#else
    gethostname(host, 128);  host[127] = 0;
#endif
    char devstr[128] = "";
    int devnum;
    if ( cudaGetDevice(&devnum) == cudaSuccess ) {
      sprintf(devstr, " device %d", devnum);
    }
    char errmsg[1024];
    sprintf(errmsg,"CUDA error %s on Pe %d (%s%s): %s", msg, CkMyPe(), host, devstr, cudaGetErrorString(err));
    NAMD_die(errmsg);
  }
}

/*
void cuda_die(const char *msg) {
    char host[128];
#ifdef NOHOSTNAME
    sprintf(host,"physical node %d", CmiPhysicalNodeID(CkMyPe()));
#else
    gethostname(host, 128);  host[127] = 0;
#endif
    char devstr[128] = "";
    int devnum;
    if ( cudaGetDevice(&devnum) == cudaSuccess ) {
      sprintf(devstr, " device %d", devnum);
    }
    char errmsg[1024];
    sprintf(errmsg,"CUDA error on Pe %d (%s%s): %s", CkMyPe(), host, devstr, msg);
    NAMD_die(errmsg);
}
*/

/*
char *devicelist;
static __thread int usedevicelist;
static __thread int ignoresharing;
static __thread int mergegrids;
static __thread int nomergegrids;
static __thread int nostreaming;

void cuda_getargs(char **argv) {
  devicelist = 0;
  usedevicelist = CmiGetArgStringDesc(argv, "+devices", &devicelist,
	"comma-delimited list of CUDA device numbers such as 0,2,1,2");
  ignoresharing = CmiGetArgFlag(argv, "+ignoresharing");
  mergegrids = CmiGetArgFlag(argv, "+mergegrids");
  nomergegrids = CmiGetArgFlag(argv, "+nomergegrids");
  if ( mergegrids && nomergegrids ) NAMD_die("Do not specify both +mergegrids and +nomergegrids");
  nostreaming = CmiGetArgFlag(argv, "+nostreaming");
}

static __thread int shared_gpu;
static __thread int first_pe_sharing_gpu;
static __thread int next_pe_sharing_gpu;
static __thread int devicePe;
static __thread int numPesSharingDevice;
static __thread int *pesSharingDevice;

static __thread int gpu_is_mine;

int cuda_device_pe() { return devicePe; }

bool cuda_device_shared_with_pe(int pe) {
  for ( int i=0; i<numPesSharingDevice; ++i ) {
    if ( pesSharingDevice[i] == pe ) return true;
  }
  return false;
}

bool one_cuda_device_per_node() {
  if ( numPesSharingDevice != CkMyNodeSize() ) return false;
  int numPesOnNodeSharingDevice = 0;
  for ( int i=0; i<numPesSharingDevice; ++i ) {
    if ( CkNodeOf(pesSharingDevice[i]) == CkMyNode() ) {
      ++numPesOnNodeSharingDevice;
    }
  }
  return ( numPesOnNodeSharingDevice == CkMyNodeSize() );
}
*/

static inline bool sortop_bitreverse(int a, int b) {
  if ( a == b ) return 0; 
  for ( int bit = 1; bit; bit *= 2 ) {
    if ( (a&bit) != (b&bit) ) return ((a&bit) < (b&bit));
  }
  return 0;
}

/*
BASE
2 types (remote & local)
16 pes per node
3 phases (1, 2, 3)
*/

/*
#define CUDA_EVENT_ID_POLL_REMOTE 98
#define CUDA_TRACE_POLL_REMOTE \
  traceUserEvent(CUDA_EVENT_ID_POLL_REMOTE)
#define CUDA_EVENT_ID_POLL_LOCAL 99
#define CUDA_TRACE_POLL_LOCAL \
  traceUserEvent(CUDA_EVENT_ID_POLL_LOCAL)
#define CUDA_EVENT_ID_BASE 100
#define CUDA_TRACE_REMOTE(START,END) \
  do { int dev; cudaGetDevice(&dev); traceUserBracketEvent( \
       CUDA_EVENT_ID_BASE + 2 * dev, START, END); } while (0)
#define CUDA_TRACE_LOCAL(START,END) \
  do { int dev; cudaGetDevice(&dev); traceUserBracketEvent( \
       CUDA_EVENT_ID_BASE + 2 * dev + 1, START, END); } while (0)

void cuda_register_user_events() {

  traceRegisterUserEvent("CUDA poll remote", CUDA_EVENT_ID_POLL_REMOTE);
  traceRegisterUserEvent("CUDA poll local", CUDA_EVENT_ID_POLL_LOCAL);

#define REGISTER_DEVICE_EVENTS(DEV) \
  traceRegisterUserEvent("CUDA device " #DEV " remote", CUDA_EVENT_ID_BASE + 2 * DEV); \
  traceRegisterUserEvent("CUDA device " #DEV " local", CUDA_EVENT_ID_BASE + 2 * DEV + 1);

  REGISTER_DEVICE_EVENTS(0)
  REGISTER_DEVICE_EVENTS(1)
  REGISTER_DEVICE_EVENTS(2)
  REGISTER_DEVICE_EVENTS(3)
  REGISTER_DEVICE_EVENTS(4)
  REGISTER_DEVICE_EVENTS(5)
  REGISTER_DEVICE_EVENTS(6)
  REGISTER_DEVICE_EVENTS(7)
  REGISTER_DEVICE_EVENTS(8)
  REGISTER_DEVICE_EVENTS(9)
  REGISTER_DEVICE_EVENTS(10)
  REGISTER_DEVICE_EVENTS(11)
  REGISTER_DEVICE_EVENTS(12)
  REGISTER_DEVICE_EVENTS(13)
  REGISTER_DEVICE_EVENTS(14)
  REGISTER_DEVICE_EVENTS(15)

}

void cuda_initialize() {

  if ( 0 == CkMyPe() ) cuda_register_user_events();

  if ( 0 == CkMyPe() ) CkPrintf("Info: Built with CUDA version %d\n", CUDA_VERSION);

  char host[128];
#ifdef NOHOSTNAME
  sprintf(host,"physical node %d", CmiPhysicalNodeID(CkMyPe()));
#else
  gethostname(host, 128);  host[127] = 0;
#endif

  int myPhysicalNodeID = CmiPhysicalNodeID(CkMyPe());
  int myRankInPhysicalNode;
  int numPesOnPhysicalNode;
  int *pesOnPhysicalNode;
  CmiGetPesOnPhysicalNode(myPhysicalNodeID,
                           &pesOnPhysicalNode,&numPesOnPhysicalNode);

  {
    int i;
    for ( i=0; i < numPesOnPhysicalNode; ++i ) {
      if ( i && (pesOnPhysicalNode[i] <= pesOnPhysicalNode[i-1]) ) {
        i = numPesOnPhysicalNode;
        break;
      }
      if ( pesOnPhysicalNode[i] == CkMyPe() ) break;
    }
    if ( i == numPesOnPhysicalNode || i != CmiPhysicalRank(CkMyPe()) ) {
      CkPrintf("Bad result from CmiGetPesOnPhysicalNode!\n");
      for ( i=0; i < numPesOnPhysicalNode; ++i ) {
        CkPrintf("pe %d physnode rank %d of %d is %d\n", CkMyPe(),
          i, numPesOnPhysicalNode, pesOnPhysicalNode[i]);
      }
      myRankInPhysicalNode = 0;
      numPesOnPhysicalNode = 1;
      pesOnPhysicalNode = new int[1];
      pesOnPhysicalNode[0] = CkMyPe();
    } else {
      myRankInPhysicalNode = i;
    }
  }
  // CkPrintf("Pe %d ranks %d in physical node\n",CkMyPe(),myRankInPhysicalNode);

  int deviceCount = 0;
  cudaGetDeviceCount(&deviceCount);
  cuda_errcheck("in cudaGetDeviceCount");
  if ( deviceCount <= 0 ) {
    cuda_die("No CUDA devices found.");
  }

  int *devices;
  int ndevices = 0;
  int nexclusive = 0;
  if ( usedevicelist ) {
    devices = new int[strlen(devicelist)];
    int i = 0;
    while ( devicelist[i] ) {
      ndevices += sscanf(devicelist+i,"%d",devices+ndevices);
      while ( devicelist[i] && isdigit(devicelist[i]) ) ++i;
      while ( devicelist[i] && ! isdigit(devicelist[i]) ) ++i;
    }
  } else {
    if ( ! CkMyPe() ) {
      CkPrintf("Did not find +devices i,j,k,... argument, using all\n");
    }
    devices = new int[deviceCount];
    for ( int i=0; i<deviceCount; ++i ) {
      int dev = i % deviceCount;
#if CUDA_VERSION >= 2020
      cudaDeviceProp deviceProp;
      cudaGetDeviceProperties(&deviceProp, dev);
      cuda_errcheck("in cudaGetDeviceProperties");
      if ( deviceProp.computeMode != cudaComputeModeProhibited
           && (deviceProp.major >= 3)
           && deviceProp.canMapHostMemory
           && ( (deviceProp.multiProcessorCount > 2) ||
                ((ndevices==0)&&(CkNumNodes()==1)) ) // exclude weak cards
         ) {
        devices[ndevices++] = dev;
      }
      if ( deviceProp.computeMode == cudaComputeModeExclusive ) {
        ++nexclusive;
      }
#else
      devices[ndevices++] = dev;
#endif
    }
  }

  if ( ! ndevices ) {
    cuda_die("all devices are in prohibited mode, of compute capability < 3.0, unable to map host memory, too small, or otherwise unusable");
  }

  shared_gpu = 0;
  gpu_is_mine = 1;
  first_pe_sharing_gpu = CkMyPe();
  next_pe_sharing_gpu = CkMyPe();

 // if ( (ndevices >= numPesOnPhysicalNode) || (nexclusive == 0) ) 
  {

  int dev;
  if ( numPesOnPhysicalNode > 1 ) {
    int myDeviceRank = myRankInPhysicalNode * ndevices / numPesOnPhysicalNode;
    dev = devices[myDeviceRank];
    devicePe = CkMyPe();
    if ( ignoresharing ) {
      pesSharingDevice = new int[1];
      pesSharingDevice[0] = CkMyPe();
      numPesSharingDevice = 1;
    } else {
      pesSharingDevice = new int[numPesOnPhysicalNode];
      devicePe = -1;
      numPesSharingDevice = 0;
      for ( int i = 0; i < numPesOnPhysicalNode; ++i ) {
        if ( i * ndevices / numPesOnPhysicalNode == myDeviceRank ) {
          int thisPe = pesOnPhysicalNode[i];
          pesSharingDevice[numPesSharingDevice++] = thisPe;
          if ( devicePe < 1 ) devicePe = thisPe;
          if ( WorkDistrib::pe_sortop_diffuse()(thisPe,devicePe) ) devicePe = thisPe;
        }
      }
      for ( int j = 0; j < ndevices; ++j ) {
        if ( devices[j] == dev && j != myDeviceRank ) shared_gpu = 1;
      }
    }
    if ( shared_gpu && devicePe == CkMyPe() ) {
      if ( CmiPhysicalNodeID(devicePe) < 2 )
      CkPrintf("Pe %d sharing CUDA device %d\n", CkMyPe(), dev);
    }
  } else {  // in case phys node code is lying
    dev = devices[CkMyPe() % ndevices];
    devicePe = CkMyPe();
    pesSharingDevice = new int[1];
    pesSharingDevice[0] = CkMyPe();
    numPesSharingDevice = 1;
  }

  if ( devicePe != CkMyPe() ) {
    if ( CmiPhysicalNodeID(devicePe) < 2 )
    CkPrintf("Pe %d physical rank %d will use CUDA device of pe %d\n",
             CkMyPe(), myRankInPhysicalNode, devicePe);

    // for PME only
    cudaError_t err = cudaSetDevice(dev);
    if ( err != cudaSuccess ) {
      char errmsg[1024];
      sprintf(errmsg,"CUDA error binding to device %d on pe %d: %s",
  			dev, CkMyPe(), cudaGetErrorString(err));
      NAMD_die(errmsg);
    }

    return;
  }

  // disable token-passing but don't submit local until remote finished
  // if shared_gpu is true, otherwise submit all work immediately
  first_pe_sharing_gpu = CkMyPe();
  next_pe_sharing_gpu = CkMyPe();

  gpu_is_mine = ( first_pe_sharing_gpu == CkMyPe() ); 

  if ( dev >= deviceCount ) {
    char buf[256];
    sprintf(buf,"Pe %d unable to bind to CUDA device %d on %s because only %d devices are present",
		CkMyPe(), dev, host, deviceCount);
    NAMD_die(buf);
  }

  cudaError_t err;
  cudaDeviceProp deviceProp;
  err = cudaGetDeviceProperties(&deviceProp, dev);
  if (err == cudaSuccess) {
    if ( CmiPhysicalNodeID(devicePe) < 2 )
    CkPrintf("Pe %d physical rank %d binding to CUDA device %d on %s: '%s'  Mem: %dMB  Rev: %d.%d\n",
             CkMyPe(), myRankInPhysicalNode, dev, host,
             deviceProp.name, deviceProp.totalGlobalMem / (1024*1024),
             deviceProp.major, deviceProp.minor);

    err = cudaSetDevice(dev);
  }
  if ( err != cudaSuccess) {
    char errmsg[1024];
    sprintf(errmsg,"CUDA error binding to device %d on pe %d: %s",
			dev, CkMyPe(), cudaGetErrorString(err));
    NAMD_die(errmsg);
  }

 }  // just let CUDA pick a device for us

  cudaSetDeviceFlags(cudaDeviceMapHost);
  cuda_errcheck("in cudaSetDeviceFlags");

  int dev;
  if ( cudaGetDevice(&dev) == cudaSuccess ) {
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, dev);
    cuda_errcheck("in cudaGetDeviceProperties");
    if ( deviceProp.computeMode == cudaComputeModeProhibited )
      cuda_die("device in prohibited mode");
    if ( deviceProp.major < 3 )
      cuda_die("device not of compute capability 3.0 or higher");
    if ( ! deviceProp.canMapHostMemory )
      cuda_die("device cannot map host memory");
  }

}
*/

static __thread ComputeNonbondedCUDA* cudaCompute = 0;
static __thread ComputeMgr *computeMgr = 0;

void send_build_cuda_force_table() {
  computeMgr->sendBuildCudaForceTable();
}

void build_cuda_force_table() {
  if ( deviceCUDA->getMasterPe() != CkMyPe() ) return;
  ComputeNonbondedCUDA::build_lj_table();
  ComputeNonbondedCUDA::build_force_table();
}

void ComputeNonbondedCUDA::build_lj_table() {  // static

  const LJTable* const ljTable = ComputeNonbondedUtil:: ljTable;
  const int dim = ljTable->get_table_dim();

  // round dim up to odd multiple of 16
  int tsize = (((dim+16+31)/32)*32)-16;
  if ( tsize < dim ) NAMD_bug("ComputeNonbondedCUDA::build_lj_table bad tsize");

  float2 *t = new float2[tsize*tsize];
  float2 *row = t;
  for ( int i=0; i<dim; ++i, row += tsize ) {
    for ( int j=0; j<dim; ++j ) {
      const LJTable::TableEntry *e = ljTable->table_val(i,j);
      row[j].x = e->A * scaling;
      row[j].y = e->B * scaling;
    }
  }

  cuda_bind_lj_table(t,tsize);
  delete [] t;

  if ( ! CmiPhysicalNodeID(CkMyPe()) ) {
    CkPrintf("Info: Updated CUDA LJ table with %d x %d elements.\n", dim, dim);
  }
}

void ComputeNonbondedCUDA::build_force_table() {  // static

  float4 t[FORCE_TABLE_SIZE];
  float4 et[FORCE_TABLE_SIZE];  // energy table

  const BigReal r2_delta = ComputeNonbondedUtil:: r2_delta;
  const int r2_delta_exp = ComputeNonbondedUtil:: r2_delta_exp;
  // const int r2_delta_expc = 64 * (r2_delta_exp - 127);
  const int r2_delta_expc = 64 * (r2_delta_exp - 1023);

  double r2list[FORCE_TABLE_SIZE];  // double to match cpu code
  for ( int i=1; i<FORCE_TABLE_SIZE; ++i ) {
    double r = ((double) FORCE_TABLE_SIZE) / ( (double) i + 0.5 );
    r2list[i] = r*r + r2_delta;
  }

  union { double f; int32 i[2]; } byte_order_test;
  byte_order_test.f = 1.0;  // should occupy high-order bits only
  int32 *r2iilist = (int32*)r2list + ( byte_order_test.i[0] ? 0 : 1 );

  for ( int i=1; i<FORCE_TABLE_SIZE; ++i ) {
    double r = ((double) FORCE_TABLE_SIZE) / ( (double) i + 0.5 );
    int table_i = (r2iilist[2*i] >> 14) + r2_delta_expc;  // table_i >= 0

    if ( r > cutoff ) {
      t[i].x = 0.;
      t[i].y = 0.;
      t[i].z = 0.;
      t[i].w = 0.;
      et[i].x = 0.;
      et[i].y = 0.;
      et[i].z = 0.;
      et[i].w = 0.;
      continue;
    }

    BigReal diffa = r2list[i] - r2_table[table_i];

    // coulomb 1/r or fast force
    // t[i].x = 1. / (r2 * r);  // -1/r * d/dr r^-1
    {
      BigReal table_a = fast_table[4*table_i];
      BigReal table_b = fast_table[4*table_i+1];
      BigReal table_c = fast_table[4*table_i+2];
      BigReal table_d = fast_table[4*table_i+3];
      BigReal grad =
		( 3. * table_d * diffa + 2. * table_c ) * diffa + table_b;
      t[i].x = 2. * grad;
      BigReal ener = table_a + diffa *
		( ( table_d * diffa + table_c ) * diffa + table_b);
      et[i].x = ener;
    }


    // pme correction for slow force
    // t[i].w = 0.;
    {
      BigReal table_a = scor_table[4*table_i];
      BigReal table_b = scor_table[4*table_i+1];
      BigReal table_c = scor_table[4*table_i+2];
      BigReal table_d = scor_table[4*table_i+3];
      BigReal grad =
		( 3. * table_d * diffa + 2. * table_c ) * diffa + table_b;
      t[i].w = 2. * grad;
      BigReal ener = table_a + diffa *
		( ( table_d * diffa + table_c ) * diffa + table_b);
      et[i].w = ener;
    }


    // vdw 1/r^6
    // t[i].y = 6. / (r8);  // -1/r * d/dr r^-6
    {
      BigReal table_a = vdwb_table[4*table_i];
      BigReal table_b = vdwb_table[4*table_i+1];
      BigReal table_c = vdwb_table[4*table_i+2];
      BigReal table_d = vdwb_table[4*table_i+3];
      BigReal grad =
		( 3. * table_d * diffa + 2. * table_c ) * diffa + table_b;
      t[i].y = 2. * -1. * grad;
      BigReal ener = table_a + diffa *
		( ( table_d * diffa + table_c ) * diffa + table_b);
      et[i].y = -1. * ener;
    }


    // vdw 1/r^12
    // t[i].z = 12e / (r8 * r4 * r2);  // -1/r * d/dr r^-12
    {
      BigReal table_a = vdwa_table[4*table_i];
      BigReal table_b = vdwa_table[4*table_i+1];
      BigReal table_c = vdwa_table[4*table_i+2];
      BigReal table_d = vdwa_table[4*table_i+3];
      BigReal grad =
		( 3. * table_d * diffa + 2. * table_c ) * diffa + table_b;
      t[i].z = 2. * grad;
      BigReal ener = table_a + diffa *
		( ( table_d * diffa + table_c ) * diffa + table_b);
      et[i].z = ener;
    }

    // CkPrintf("%d %g %g %g %g %g %g\n", i, r, diffa,
    //   t[i].x, t[i].y, t[i].z, t[i].w);

/*
    double r2 = r * r;
    double r4 = r2 * r2;
    double r8 = r4 * r4;

    t[i].x = 1. / (r2 * r);  // -1/r * d/dr r^-1
    t[i].y = 6. / (r8);  // -1/r * d/dr r^-6
    t[i].z = 12. / (r8 * r4 * r2);  // -1/r * d/dr r^-12
    t[i].w = 0.;
*/
  }

  t[0].x = 0.f;
  t[0].y = 0.f;
  t[0].z = 0.f;
  t[0].w = 0.f;
  et[0].x = et[1].x;
  et[0].y = et[1].y;
  et[0].z = et[1].z;
  et[0].w = et[1].w;

  cuda_bind_force_table(t,et);

  if ( ! CmiPhysicalNodeID(CkMyPe()) ) {
    CkPrintf("Info: Updated CUDA force table with %d elements.\n", FORCE_TABLE_SIZE);
  }
}

struct exlist_sortop {
  bool operator() (int32 *li, int32 *lj) {
    return ( li[1] < lj[1] );
  }
};

void build_cuda_exclusions() {
  if ( deviceCUDA->getMasterPe() != CkMyPe() ) return;
  ComputeNonbondedCUDA::build_exclusions();
}

static __thread int2 *exclusionsByAtom;

void ComputeNonbondedCUDA::build_exclusions() {
  Molecule *mol = Node::Object()->molecule;

#ifdef MEM_OPT_VERSION
  int natoms = mol->exclSigPoolSize;
#else
  int natoms = mol->numAtoms; 
#endif

  delete [] exclusionsByAtom;
  exclusionsByAtom = new int2[natoms];

  // create unique sorted lists

  ObjectArena<int32> listArena;
  ResizeArray<int32*> unique_lists;
  int32 **listsByAtom = new int32*[natoms];
  SortableResizeArray<int32> curList;
  for ( int i=0; i<natoms; ++i ) {
    curList.resize(0);
    curList.add(0);  // always excluded from self
#ifdef MEM_OPT_VERSION
    const ExclusionSignature *sig = mol->exclSigPool + i;
    int n = sig->fullExclCnt;
    for ( int j=0; j<n; ++j ) { curList.add(sig->fullOffset[j]); }
    n += 1;
#else
    const int32 *mol_list = mol->get_full_exclusions_for_atom(i);
    int n = mol_list[0] + 1;
    for ( int j=1; j<n; ++j ) {
      curList.add(mol_list[j] - i);
    }
#endif
    curList.sort();

    int j;
    for ( j=0; j<unique_lists.size(); ++j ) {
      if ( n != unique_lists[j][0] ) continue;  // no match
      int k;
      for ( k=0; k<n; ++k ) {
        if ( unique_lists[j][k+3] != curList[k] ) break;
      }
      if ( k == n ) break;  // found match
    }
    if ( j == unique_lists.size() ) {  // no match
      int32 *list = listArena.getNewArray(n+3);
      list[0] = n;
      int maxdiff = 0;
      maxdiff = -1 * curList[0];
      if ( curList[n-1] > maxdiff ) maxdiff = curList[n-1];
      list[1] = maxdiff;
      for ( int k=0; k<n; ++k ) {
        list[k+3] = curList[k];
      }
      unique_lists.add(list);
    }
    listsByAtom[i] = unique_lists[j];
  }
  // sort lists by maxdiff
  std::stable_sort(unique_lists.begin(), unique_lists.end(), exlist_sortop());
  long int totalbits = 0;
  int nlists = unique_lists.size();
  for ( int j=0; j<nlists; ++j ) {
    int32 *list = unique_lists[j];
    int maxdiff = list[1];
    list[2] = totalbits + maxdiff;
    totalbits += 2*maxdiff + 1;
  }
  for ( int i=0; i<natoms; ++i ) {
    exclusionsByAtom[i].x = listsByAtom[i][1];  // maxdiff
    exclusionsByAtom[i].y = listsByAtom[i][2];  // start
  }
  delete [] listsByAtom;

  if ( totalbits & 31 ) totalbits += ( 32 - ( totalbits & 31 ) );

  {
    long int bytesneeded = totalbits / 8;
    if ( ! CmiPhysicalNodeID(CkMyPe()) ) {
    CkPrintf("Info: Found %d unique exclusion lists needing %ld bytes\n",
		unique_lists.size(), bytesneeded);
    }

    long int bytesavail = MAX_EXCLUSIONS * sizeof(unsigned int);
    if ( bytesneeded > bytesavail ) {
      char errmsg[512];
      sprintf(errmsg,"Found %d unique exclusion lists needing %ld bytes "
                     "but only %ld bytes can be addressed with 32-bit int.",
                     unique_lists.size(), bytesneeded, bytesavail);
      NAMD_die(errmsg);
    }
  }

#define SET_EXCL(EXCL,BASE,DIFF) \
         (EXCL)[((BASE)+(DIFF))>>5] |= (1<<(((BASE)+(DIFF))&31))

  unsigned int *exclusion_bits = new unsigned int[totalbits/32];
  memset(exclusion_bits, 0, totalbits/8);

  long int base = 0;
  for ( int i=0; i<unique_lists.size(); ++i ) {
    base += unique_lists[i][1];
    if ( unique_lists[i][2] != (int32)base ) {
      NAMD_bug("ComputeNonbondedCUDA::build_exclusions base != stored");
    }
    int n = unique_lists[i][0];
    for ( int j=0; j<n; ++j ) {
      SET_EXCL(exclusion_bits,base,unique_lists[i][j+3]);
    }
    base += unique_lists[i][1] + 1;
  }

  cuda_bind_exclusions(exclusion_bits, totalbits/32);

  delete [] exclusion_bits;
}


void register_cuda_compute_self(ComputeID c, PatchID pid) {

  if ( ! cudaCompute ) NAMD_bug("register_self called early");

  cudaCompute->requirePatch(pid);

  ComputeNonbondedCUDA::compute_record cr;
  cr.c = c;
  cr.pid[0] = pid;  cr.pid[1] = pid;
  cr.offset = 0.;
  if ( cudaCompute->patchRecords[pid].isLocal ) {
    cudaCompute->localComputeRecords.add(cr);
  } else {
    cudaCompute->remoteComputeRecords.add(cr);
  }
}

void register_cuda_compute_pair(ComputeID c, PatchID pid[], int t[]) {

  if ( ! cudaCompute ) NAMD_bug("register_pair called early");
 
  cudaCompute->requirePatch(pid[0]);
  cudaCompute->requirePatch(pid[1]);

  ComputeNonbondedCUDA::compute_record cr;
  cr.c = c; 
  cr.pid[0] = pid[0];  cr.pid[1] = pid[1];

  int t1 = t[0];
  int t2 = t[1];
  Vector offset = cudaCompute->patchMap->center(pid[0])
                - cudaCompute->patchMap->center(pid[1]);
  offset.x += (t1%3-1) - (t2%3-1);
  offset.y += ((t1/3)%3-1) - ((t2/3)%3-1);
  offset.z += (t1/9-1) - (t2/9-1);
  cr.offset = offset;

  if ( cudaCompute->patchRecords[pid[0]].isLocal ) {
    cudaCompute->localComputeRecords.add(cr);
  } else {
    cudaCompute->remoteComputeRecords.add(cr);
  }
}

void unregister_cuda_compute(ComputeID c) {  // static

  NAMD_bug("unregister_compute unimplemented");

}

// static __thread cudaEvent_t start_upload;
static __thread cudaEvent_t start_calc;
static __thread cudaEvent_t end_remote_download;
static __thread cudaEvent_t end_local_download;

static __thread ResizeArray<patch_pair> *patch_pairs_ptr;
static __thread ResizeArray<int> *patch_pair_num_ptr;

void init_arrays();

ComputeNonbondedCUDA::ComputeNonbondedCUDA(ComputeID c, ComputeMgr *mgr,
					   ComputeNonbondedCUDA *m, int idx) : Compute(c), slaveIndex(idx) {
#ifdef PRINT_GBIS
   CkPrintf("C.N.CUDA[%d]::constructor cid=%d\n", CkMyPe(), c);
#endif

  if ( sizeof(patch_pair) & 15 ) NAMD_bug("sizeof(patch_pair) % 16 != 0");
  if ( sizeof(atom) & 15 ) NAMD_bug("sizeof(atom) % 16 != 0");
  if ( sizeof(atom_param) & 15 ) NAMD_bug("sizeof(atom_param) % 16 != 0");

  // CkPrintf("create ComputeNonbondedCUDA\n");
  master = m ? m : this;
  cudaCompute = this;
  computeMgr = mgr;
  patchMap = PatchMap::Object();
  atomMap = AtomMap::Object();
  reduction = 0;

  SimParameters *params = Node::Object()->simParameters;
  if (params->pressureProfileOn) {
    NAMD_die("pressure profile not supported in CUDA");
  }

  atomsChanged = 1;
  computesChanged = 1;
  patchPairsReordered = 1;

  pairlistsValid = 0;
  pairlistTolerance = 0.;
  usePairlists = 0;
  savePairlists = 0;
  plcutoff2 = 0.;

  workStarted = 0;
  basePriority = PROXY_DATA_PRIORITY;
  localWorkMsg2 = new (PRIORITY_SIZE) LocalWorkMsg;

  // Zero array sizes and pointers
  init_arrays();

  atoms_size = 0;
  atoms = NULL;

  forces_size = 0;
  forces = NULL;
  
  slow_forces_size = 0;
  slow_forces = NULL;

  psiSumH_size = 0;
  psiSumH = NULL;

  dEdaSumH_size = 0;
  dEdaSumH = NULL;

  if ( master != this ) { // I am slave
    masterPe = master->masterPe;
    master->slaves[slaveIndex] = this;
    if ( master->slavePes[slaveIndex] != CkMyPe() ) {
      NAMD_bug("ComputeNonbondedCUDA slavePes[slaveIndex] != CkMyPe");
    }
    deviceID = master->deviceID;
    registerPatches();
    return;
  }
  masterPe = CkMyPe();

  const bool streaming = ! (deviceCUDA->getNoStreaming() || params->GBISOn);
  if ( streaming && ! deviceCUDA->getSharedGpu() && ! deviceCUDA->getNoMergeGrids() )
    deviceCUDA->setMergeGrids(1);

  // Sanity checks for New CUDA kernel
  if (params->GBISOn) {
    // GBIS
    if (deviceCUDA->getNoMergeGrids()) {
      NAMD_die("CUDA kernel cannot use +nomergegrids with GBIS simulations");
    }
    // Set mergegrids ON as long as user hasn't defined +nomergegrids
    if (!deviceCUDA->getNoMergeGrids()) deviceCUDA->setMergeGrids(1);
    // Final sanity check
    if (!deviceCUDA->getMergeGrids()) NAMD_die("CUDA GBIS kernel final sanity check failed");
  } else {
    // non-GBIS
    if (deviceCUDA->getNoStreaming() && !deviceCUDA->getMergeGrids()) {
      NAMD_die("CUDA kernel requires +mergegrids if +nostreaming is used");
    }
  }

#if CUDA_VERSION >= 5050
  int leastPriority, greatestPriority;
  cudaDeviceGetStreamPriorityRange(&leastPriority, &greatestPriority);
  cuda_errcheck("in cudaDeviceGetStreamPriorityRange");
  if ( leastPriority != greatestPriority ) {
    if ( CkMyNode() == 0 ) {
      int dev = deviceCUDA->getDeviceID();
      CkPrintf("CUDA device %d stream priority range %d %d\n", dev, leastPriority, greatestPriority);
    }
    if ( deviceCUDA->getMergeGrids() && params->PMEOn && params->PMEOffload && !params->usePMECUDA) {
      greatestPriority = leastPriority;
    }
    if (params->usePMECUDA) greatestPriority = leastPriority;
    cudaStreamCreateWithPriority(&stream,cudaStreamDefault,greatestPriority);
    cudaStreamCreateWithPriority(&stream2,cudaStreamDefault,leastPriority);
  } else
#endif
  {
    cudaStreamCreate(&stream);
    cuda_errcheck("cudaStreamCreate");
    int dev = deviceCUDA->getDeviceID();
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, dev);
    cuda_errcheck("cudaGetDeviceProperties");
    if ( deviceProp.concurrentKernels && deviceProp.major > 2 ) {
      if ( CkMyNode() == 0 ) CkPrintf("CUDA device %d supports concurrent kernels.\n", dev);
      cudaStreamCreate(&stream2);
    } else {
      stream2 = stream;
    }
  }
  cuda_errcheck("cudaStreamCreate");

  // Get GPU device ID
  deviceID = deviceCUDA->getDeviceID();

  cuda_init();
  if ( max_grid_size < 65535 ) NAMD_bug("bad CUDA max_grid_size");
  build_exclusions();
  // cudaEventCreate(&start_upload);
  cudaEventCreateWithFlags(&start_calc,cudaEventDisableTiming);
  cudaEventCreateWithFlags(&end_remote_download,cudaEventDisableTiming);
  cudaEventCreateWithFlags(&end_local_download,cudaEventDisableTiming);

  patchRecords.resize(patchMap->numPatches());
  patch_pairs_ptr = new ResizeArray<patch_pair>;
  patch_pair_num_ptr = new ResizeArray<int>;

  if ( params->PMEOn && params->PMEOffload && !params->usePMECUDA) deviceCUDA->setGpuIsMine(0);
}


ComputeNonbondedCUDA::~ComputeNonbondedCUDA() { ; }

void ComputeNonbondedCUDA::requirePatch(int pid) {

  computesChanged = 1;
  patch_record &pr = patchRecords[pid];
  if ( pr.refCount == 0 ) {
    pr.isSamePhysicalNode = ( CmiPhysicalNodeID(patchMap->node(pid)) == CmiPhysicalNodeID(CkMyPe()) );
    pr.isSameNode = ( CkNodeOf(patchMap->node(pid)) == CkMyNode() );
    if ( deviceCUDA->getMergeGrids() ) {
      pr.isLocal = 0;
    } else if ( CkNumNodes() < 2 ) {
      pr.isLocal = 1 & ( 1 ^ patchMap->index_a(pid) ^
         patchMap->index_b(pid) ^ patchMap->index_c(pid) );
    } else {
      pr.isLocal = pr.isSameNode;
    }
    if ( pr.isLocal ) {
      localActivePatches.add(pid);
    } else {
      remoteActivePatches.add(pid);
    }
    activePatches.add(pid);
    pr.patchID = pid;
    pr.hostPe = -1;
    pr.x = NULL;
    pr.xExt = NULL;
    pr.r = NULL;
    pr.f = NULL;
    pr.intRad      = NULL;
    pr.psiSum      = NULL;
    pr.bornRad     = NULL;
    pr.dEdaSum     = NULL;
    pr.dHdrPrefix  = NULL;
  }
  pr.refCount += 1;
}

void ComputeNonbondedCUDA::registerPatches() {

  SimParameters *simParams = Node::Object()->simParameters;
  int npatches = master->activePatches.size();
  int *pids = master->activePatches.begin();
  patch_record *recs = master->patchRecords.begin();
  for ( int i=0; i<npatches; ++i ) {
    int pid = pids[i];
    patch_record &pr = recs[pid];
    if ( pr.hostPe == CkMyPe() ) {
      pr.slave = this;
      pr.msg = new (PRIORITY_SIZE) FinishWorkMsg;
      hostedPatches.add(pid);
      if ( pr.isLocal ) {
        localHostedPatches.add(pid);
      } else {
        remoteHostedPatches.add(pid);
      }
      ProxyMgr::Object()->createProxy(pid);
      pr.p = patchMap->patch(pid);
      pr.positionBox = pr.p->registerPositionPickup(this);
      pr.forceBox = pr.p->registerForceDeposit(this);
      if (simParams->GBISOn) {
        pr.intRadBox      = pr.p->registerIntRadPickup(this);
        pr.psiSumBox      = pr.p->registerPsiSumDeposit(this);
        pr.bornRadBox     = pr.p->registerBornRadPickup(this);
        pr.dEdaSumBox     = pr.p->registerDEdaSumDeposit(this);
        pr.dHdrPrefixBox  = pr.p->registerDHdrPrefixPickup(this);
      }
    }
  }
  if ( master == this ) setNumPatches(activePatches.size());
  else setNumPatches(hostedPatches.size());
  if ( CmiPhysicalNodeID(CkMyPe()) < 2 )
  CkPrintf("Pe %d hosts %d local and %d remote patches for pe %d\n", CkMyPe(), localHostedPatches.size(), remoteHostedPatches.size(), masterPe);
}

struct pid_sortop_reverse_priority {
  bool operator() (int pidj, int pidi) {  // i and j reversed
    int ppi = PATCH_PRIORITY(pidi);
    int ppj = PATCH_PRIORITY(pidj);
    if ( ppi != ppj ) return ppi < ppj;
    return pidi < pidj;
  }
};

void ComputeNonbondedCUDA::assignPatches() {

  int *pesOnNodeSharingDevice = new int[CkMyNodeSize()];
  int numPesOnNodeSharingDevice = 0;
  int masterIndex = -1;
  for ( int i=0; i<deviceCUDA->getNumPesSharingDevice(); ++i ) {
    int pe = deviceCUDA->getPesSharingDevice(i);
    if ( pe == CkMyPe() ) masterIndex = numPesOnNodeSharingDevice;
    if ( CkNodeOf(pe) == CkMyNode() ) {
      pesOnNodeSharingDevice[numPesOnNodeSharingDevice++] = pe;
    }
  }

  int npatches = activePatches.size();

  if ( npatches ) {
    reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);
  }

  // calculate priority rank of local home patch within pe
  {
    ResizeArray< ResizeArray<int> > homePatchByRank(CkMyNodeSize());
    for ( int i=0; i<npatches; ++i ) {
      int pid = activePatches[i];
      int homePe = patchMap->node(pid);
      if ( CkNodeOf(homePe) == CkMyNode() ) {
        homePatchByRank[CkRankOf(homePe)].add(pid);
      }
    }
    for ( int i=0; i<CkMyNodeSize(); ++i ) {
      pid_sortop_reverse_priority so;
      std::sort(homePatchByRank[i].begin(),homePatchByRank[i].end(),so);
      int masterBoost = ( CkMyRank() == i ? 2 : 0 );
      for ( int j=0; j<homePatchByRank[i].size(); ++j ) {
        int pid = homePatchByRank[i][j];
        patchRecords[pid].reversePriorityRankInPe = j + masterBoost;
      }
    }
  }

  int *count = new int[npatches];
  memset(count, 0, sizeof(int)*npatches);
  int *pcount = new int[numPesOnNodeSharingDevice];
  memset(pcount, 0, sizeof(int)*numPesOnNodeSharingDevice);
  int *rankpcount = new int[CkMyNodeSize()];
  memset(rankpcount, 0, sizeof(int)*CkMyNodeSize());
  char *table = new char[npatches*numPesOnNodeSharingDevice];
  memset(table, 0, npatches*numPesOnNodeSharingDevice);

  int unassignedpatches = npatches;

  if ( 0 ) { // assign all to device pe
    for ( int i=0; i<npatches; ++i ) {
      int pid = activePatches[i];
      patch_record &pr = patchRecords[pid];
      pr.hostPe = CkMyPe();
    }
    unassignedpatches = 0;
    pcount[masterIndex] = npatches;
  } else 

  // assign if home pe and build table of natural proxies
  for ( int i=0; i<npatches; ++i ) {
    int pid = activePatches[i];
    patch_record &pr = patchRecords[pid];
    int homePe = patchMap->node(pid);
    for ( int j=0; j<numPesOnNodeSharingDevice; ++j ) {
      int pe = pesOnNodeSharingDevice[j];
      if ( pe == homePe ) {
        pr.hostPe = pe;  --unassignedpatches;
        pcount[j] += 1;
      }
      if ( PatchMap::ObjectOnPe(pe)->patch(pid) ) {
        table[i*numPesOnNodeSharingDevice+j] = 1;
      }
    }
    if ( pr.hostPe == -1 && CkNodeOf(homePe) == CkMyNode() ) {
      pr.hostPe = homePe;  --unassignedpatches;
      rankpcount[CkRankOf(homePe)] += 1;
    }
  }
  // assign if only one pe has a required proxy
  int assignj = 0;
  for ( int i=0; i<npatches; ++i ) {
    int pid = activePatches[i];
    patch_record &pr = patchRecords[pid];
    if ( pr.hostPe != -1 ) continue;
    int c = 0;
    int lastj;
    for ( int j=0; j<numPesOnNodeSharingDevice; ++j ) {
      if ( table[i*numPesOnNodeSharingDevice+j] ) { ++c; lastj=j; }
    }
    count[i] = c;
    if ( c == 1 ) {
      pr.hostPe = pesOnNodeSharingDevice[lastj];
      --unassignedpatches;
      pcount[lastj] += 1;
    }
  }
  while ( unassignedpatches ) {
    int i;
    for ( i=0; i<npatches; ++i ) {
      if ( ! table[i*numPesOnNodeSharingDevice+assignj] ) continue;
      int pid = activePatches[i];
      patch_record &pr = patchRecords[pid];
      if ( pr.hostPe != -1 ) continue;
      pr.hostPe = pesOnNodeSharingDevice[assignj];
      --unassignedpatches;
      pcount[assignj] += 1;
      if ( ++assignj == numPesOnNodeSharingDevice ) assignj = 0;
      break;
    }
    if ( i<npatches ) continue;  // start search again
    for ( i=0; i<npatches; ++i ) {
      int pid = activePatches[i];
      patch_record &pr = patchRecords[pid];
      if ( pr.hostPe != -1 ) continue;
      if ( count[i] ) continue;
      pr.hostPe = pesOnNodeSharingDevice[assignj];
      --unassignedpatches;
      pcount[assignj] += 1;
      if ( ++assignj == numPesOnNodeSharingDevice ) assignj = 0;
      break;
    }
    if ( i<npatches ) continue;  // start search again
    if ( ++assignj == numPesOnNodeSharingDevice ) assignj = 0;
  }

  for ( int i=0; i<npatches; ++i ) {
    int pid = activePatches[i];
    patch_record &pr = patchRecords[pid];
    // CkPrintf("Pe %d patch %d hostPe %d\n", CkMyPe(), pid, pr.hostPe);
  }

  slavePes = new int[CkMyNodeSize()];
  slaves = new ComputeNonbondedCUDA*[CkMyNodeSize()];
  numSlaves = 0;
  for ( int j=0; j<numPesOnNodeSharingDevice; ++j ) {
    int pe = pesOnNodeSharingDevice[j];
    int rank = pe - CkNodeFirst(CkMyNode());
    // CkPrintf("host %d sharing %d pe %d rank %d pcount %d rankpcount %d\n",
    //          CkMyPe(),j,pe,rank,pcount[j],rankpcount[rank]);
    if ( pe == CkMyPe() ) continue;
    if ( ! pcount[j] && ! rankpcount[rank] ) continue;
    rankpcount[rank] = 0;  // skip in rank loop below
    slavePes[numSlaves] = pe;
    computeMgr->sendCreateNonbondedCUDASlave(pe,numSlaves);
    ++numSlaves;
  }
  for ( int j=0; j<CkMyNodeSize(); ++j ) {
    int pe = CkNodeFirst(CkMyNode()) + j;
    // CkPrintf("host %d rank %d pe %d rankpcount %d\n",
    //          CkMyPe(),j,pe,rankpcount[j]);
    if ( ! rankpcount[j] ) continue;
    if ( pe == CkMyPe() ) continue;
    slavePes[numSlaves] = pe;
    computeMgr->sendCreateNonbondedCUDASlave(pe,numSlaves);
    ++numSlaves;
  }
  registerPatches();

  delete [] pesOnNodeSharingDevice;
  delete [] count;
  delete [] pcount;
  delete [] rankpcount;
  delete [] table;
}

static __thread int atom_params_size;
static __thread atom_param* atom_params;

static __thread int vdw_types_size;
static __thread int* vdw_types;

static __thread int dummy_size;
static __thread float* dummy_dev;

static __thread int force_ready_queue_size;
static __thread int *force_ready_queue;
static __thread int force_ready_queue_len;
static __thread int force_ready_queue_next;

static __thread int block_order_size;
static __thread int *block_order;

static __thread int num_atoms;
static __thread int num_local_atoms;
static __thread int num_remote_atoms;

static __thread int virials_size;
static __thread float *virials;
static __thread int num_virials;
// NOTE: slow_virials is a computed pointer from "virials" -do not deallocate
static __thread float *slow_virials;

static __thread int energy_gbis_size;
static __thread float *energy_gbis;

//GBIS host pointers
static __thread int intRad0H_size;
static __thread float *intRad0H;
static __thread int intRadSH_size;
static __thread float *intRadSH;
static __thread int bornRadH_size;
static __thread float *bornRadH;
static __thread int dHdrPrefixH_size;
static __thread float *dHdrPrefixH;

static __thread int cuda_timer_count;
static __thread double cuda_timer_total;
static __thread double kernel_time;
static __thread double remote_submit_time;
static __thread double local_submit_time;

#define CUDA_POLL(FN,ARG) CcdCallFnAfter(FN,ARG,0.1)

extern "C" void CcdCallBacksReset(void *ignored,double curWallTime);  // fix Charm++

#ifdef PRINT_GBIS
#define GBISP(...) CkPrintf(__VA_ARGS__);
#else
#define GBISP(...)
#endif

#define count_limit 1000000
static __thread int check_count;
static __thread int check_remote_count;
static __thread int check_local_count;

void init_arrays() {

  atom_params_size = 0;
  atom_params = NULL;

  vdw_types_size = 0;
  vdw_types = NULL;
  
  dummy_size = 0;
  dummy_dev = NULL;

  force_ready_queue_size = 0;
  force_ready_queue = NULL;
  force_ready_queue_len = 0;
  force_ready_queue_next = 0;
  
  block_order_size = 0;
  block_order = NULL;
  
  num_atoms = 0;
  num_local_atoms = 0;
  num_remote_atoms = 0;

  virials_size = 0;
  virials = NULL;
  num_virials = 0;

  energy_gbis_size = 0;
  energy_gbis = NULL;

  intRad0H_size = 0;
  intRad0H = NULL;
  intRadSH_size = 0;
  intRadSH = NULL;
  bornRadH_size = 0;
  bornRadH = NULL;
  dHdrPrefixH_size = 0;
  dHdrPrefixH = NULL;

}

void cuda_check_progress(void *arg, double walltime) {
  CUDA_TRACE_POLL_REMOTE;

  int flindex;
  int poll_again = 1;
  while ( -1 != (flindex = force_ready_queue[force_ready_queue_next]) ) {
    //    CkPrintf("Pe %d forces ready %d is index %d at %lf\n",
    //	     CkMyPe(), force_ready_queue_next, flindex, walltime);
    force_ready_queue[force_ready_queue_next] = -1;
    ++force_ready_queue_next;
    check_count = 0;
    if ( force_ready_queue_next == force_ready_queue_len ) {
      poll_again = 0;
      CUDA_TRACE_LOCAL(kernel_time,walltime);
      kernel_time = walltime - kernel_time;
      // need to guarantee this finishes before the last patch message!
      ((ComputeNonbondedCUDA *) arg)->workStarted = 0;
      ((ComputeNonbondedCUDA *) arg)->finishReductions();
    }
    ((ComputeNonbondedCUDA *) arg)->messageFinishPatch(flindex);
    if ( force_ready_queue_next == force_ready_queue_len ) break;
  }
  if ( ++check_count >= count_limit ) {
    char errmsg[256];
    sprintf(errmsg,"cuda_check_progress polled %d times over %f s on step %d",
            check_count, walltime - remote_submit_time,
            ((ComputeNonbondedCUDA *) arg)->step);
    cuda_errcheck(errmsg);
    NAMD_die(errmsg);
  }
  if ( poll_again ) {
    CcdCallBacksReset(0,walltime);  // fix Charm++
    CUDA_POLL(cuda_check_progress, arg);
  }
}

void cuda_check_remote_progress(void *arg, double walltime) {

  CUDA_TRACE_POLL_REMOTE;
  cudaError_t err = cudaEventQuery(end_remote_download);
  if ( err == cudaSuccess ) {
    local_submit_time = walltime;
    CUDA_TRACE_REMOTE(remote_submit_time,local_submit_time);
    if ( deviceCUDA->getMergeGrids() ) {  // no local
      kernel_time = local_submit_time - kernel_time;
    }
    check_remote_count = 0;
    cuda_errcheck("at cuda remote stream completed");
    WorkDistrib::messageFinishCUDA((ComputeNonbondedCUDA *) arg);
  } else if ( err != cudaErrorNotReady ) {
    cuda_errcheck("in cuda_check_remote_progress");
    NAMD_bug("cuda_errcheck missed error in cuda_check_remote_progress");
  } else if ( ++check_remote_count >= count_limit ) {
    char errmsg[256];
    sprintf(errmsg,"cuda_check_remote_progress polled %d times over %f s on step %d",
            check_remote_count, walltime - remote_submit_time,
            ((ComputeNonbondedCUDA *) arg)->step);
    cuda_errcheck(errmsg);
    NAMD_die(errmsg);
  } else if ( check_local_count ) {
    NAMD_bug("nonzero check_local_count in cuda_check_remote_progress");
  } else {
    CcdCallBacksReset(0,walltime);  // fix Charm++
    CUDA_POLL(cuda_check_remote_progress, arg);
  }
}

void cuda_check_local_progress(void *arg, double walltime) {

  CUDA_TRACE_POLL_LOCAL;
  cudaError_t err = cudaEventQuery(end_local_download);
  if ( err == cudaSuccess ) {
    CUDA_TRACE_LOCAL(local_submit_time,walltime);
    kernel_time = walltime - kernel_time;
    check_local_count = 0;
    cuda_errcheck("at cuda local stream completed");
    WorkDistrib::messageFinishCUDA((ComputeNonbondedCUDA *) arg);
  } else if ( err != cudaErrorNotReady ) {
    cuda_errcheck("in cuda_check_local_progress");
    NAMD_bug("cuda_errcheck missed error in cuda_check_local_progress");
  } else if ( ++check_local_count >= count_limit ) {
    char errmsg[256];
    sprintf(errmsg,"cuda_check_local_progress polled %d times over %f s on step %d",
            check_local_count, walltime - local_submit_time,
            ((ComputeNonbondedCUDA *) arg)->step);
    cuda_errcheck(errmsg);
    NAMD_die(errmsg);
  } else if ( check_remote_count ) {
    NAMD_bug("nonzero check_remote_count in cuda_check_local_progress");
  } else {
    CcdCallBacksReset(0,walltime);  // fix Charm++
    CUDA_POLL(cuda_check_local_progress, arg);
  }
}

#if 0
// don't use this one unless timer is part of stream, above is better
void cuda_check_progress(void *arg, double walltime) {
  if ( cuda_stream_finished() ) {
    kernel_time = walltime - kernel_time;
    CcdCallBacksReset(0,walltime);  // fix Charm++
    CUDA_POLL(ccd_index);
    // ((ComputeNonbondedCUDA *) arg)->finishWork();
    WorkDistrib::messageEnqueueWork((ComputeNonbondedCUDA *) arg);
  }
}
#endif

void ComputeNonbondedCUDA::atomUpdate() {
  //fprintf(stderr, "%d ComputeNonbondedCUDA::atomUpdate\n",CkMyPe());
  atomsChanged = 1;
}

static __thread int kernel_launch_state = 0;

struct cr_sortop_distance {
  const Lattice &l;
  cr_sortop_distance(const Lattice &lattice) : l(lattice) { }
  bool operator() (ComputeNonbondedCUDA::compute_record i,
			ComputeNonbondedCUDA::compute_record j) {
    Vector a = l.a();
    Vector b = l.b();
    Vector c = l.c();
    BigReal ri = (i.offset.x * a + i.offset.y * b + i.offset.z * c).length2();
    BigReal rj = (j.offset.x * a + j.offset.y * b + j.offset.z * c).length2();
    return ( ri < rj );
  }
};

struct cr_sortop_reverse_priority {
  cr_sortop_distance &distop;
  const ComputeNonbondedCUDA::patch_record *pr;
  cr_sortop_reverse_priority(cr_sortop_distance &sod,
       const ComputeNonbondedCUDA::patch_record *patchrecs) : distop(sod), pr(patchrecs) { }
  bool pid_compare_priority(int pidi, int pidj) {
    const ComputeNonbondedCUDA::patch_record &pri = pr[pidi];
    const ComputeNonbondedCUDA::patch_record &prj = pr[pidj];
    if ( pri.isSamePhysicalNode && ! prj.isSamePhysicalNode ) return 0;
    if ( prj.isSamePhysicalNode && ! pri.isSamePhysicalNode ) return 1;
    if ( pri.isSameNode && ! prj.isSameNode ) return 0;
    if ( prj.isSameNode && ! pri.isSameNode ) return 1;
    if ( pri.isSameNode ) {  // and prj.isSameNode
      int rpri = pri.reversePriorityRankInPe;
      int rprj = prj.reversePriorityRankInPe;
      if ( rpri != rprj ) return rpri > rprj;
      return sortop_bitreverse(CkRankOf(pri.hostPe),CkRankOf(prj.hostPe));
    }
    int ppi = PATCH_PRIORITY(pidi);
    int ppj = PATCH_PRIORITY(pidj);
    if ( ppi != ppj ) return ppi < ppj;
    return pidi < pidj;
  }
  bool operator() (ComputeNonbondedCUDA::compute_record j,
			ComputeNonbondedCUDA::compute_record i) {  // i and j reversed
    int pidi = pid_compare_priority(i.pid[0],i.pid[1]) ? i.pid[0] : i.pid[1];
    int pidj = pid_compare_priority(j.pid[0],j.pid[1]) ? j.pid[0] : j.pid[1];
    if ( pidi != pidj ) return pid_compare_priority(pidi, pidj);
    return distop(i,j);
  }
};

void ComputeNonbondedCUDA::skip() {
  //fprintf(stderr, "ComputeNonbondedCUDA::skip()\n");
  SimParameters *simParams = Node::Object()->simParameters;
  for ( int i=0; i<hostedPatches.size(); ++i ) {
    patch_record &pr = master->patchRecords[hostedPatches[i]];
    pr.positionBox->skip();
    pr.forceBox->skip();
    if (simParams->GBISOn) {
      pr.intRadBox->skip();
      pr.psiSumBox->skip();
      pr.bornRadBox->skip();
      pr.dEdaSumBox->skip();
      pr.dHdrPrefixBox->skip();
    }
  }
}

int ComputeNonbondedCUDA::noWork() {

  SimParameters *simParams = Node::Object()->simParameters;
  Flags &flags = master->patchRecords[hostedPatches[0]].p->flags;
  lattice = flags.lattice;
  doSlow = flags.doFullElectrostatics;
  doEnergy = flags.doEnergy;
  step = flags.step;

  if ( ! flags.doNonbonded ) {
    GBISP("GBIS[%d] noWork() don't do nonbonded\n",CkMyPe());
    if ( master != this ) {
      computeMgr->sendNonbondedCUDASlaveReady(masterPe,
			hostedPatches.size(),atomsChanged,sequence());
    } else {
      for ( int i = 0; i < numSlaves; ++i ) {
        computeMgr->sendNonbondedCUDASlaveSkip(slaves[i],slavePes[i]);
      }
      skip();
    }
    if ( reduction ) {
       reduction->submit();
    }

    return 1;
  }

  for ( int i=0; i<hostedPatches.size(); ++i ) {
    patch_record &pr = master->patchRecords[hostedPatches[i]];
    if (!simParams->GBISOn || gbisPhase == 1) {
      GBISP("GBIS[%d] noWork() P0[%d] open()\n",CkMyPe(), pr.patchID);
      pr.x = pr.positionBox->open();
      pr.xExt = pr.p->getCompAtomExtInfo();
    }

    if (simParams->GBISOn) {
      if (gbisPhase == 1) {
        GBISP("GBIS[%d] noWork() P1[%d] open()\n",CkMyPe(),pr.patchID);
        pr.intRad     = pr.intRadBox->open();
        pr.psiSum     = pr.psiSumBox->open();
      } else if (gbisPhase == 2) {
        GBISP("GBIS[%d] noWork() P2[%d] open()\n",CkMyPe(),pr.patchID);
        pr.bornRad    = pr.bornRadBox->open();
        pr.dEdaSum    = pr.dEdaSumBox->open();
      } else if (gbisPhase == 3) {
        GBISP("GBIS[%d] noWork() P3[%d] open()\n",CkMyPe(),pr.patchID);
        pr.dHdrPrefix = pr.dHdrPrefixBox->open();
      }
      GBISP("opened GBIS boxes");
    }
  }

  if ( master == this ) return 0; //work to do, enqueue as usual

  // message masterPe
  computeMgr->sendNonbondedCUDASlaveReady(masterPe,
			hostedPatches.size(),atomsChanged,sequence());
  atomsChanged = 0;

  workStarted = 1;

  return 1;
}

void ComputeNonbondedCUDA::doWork() {
GBISP("C.N.CUDA[%d]::doWork: seq %d, phase %d, workStarted %d, atomsChanged %d\n", \
CkMyPe(), sequence(), gbisPhase, workStarted, atomsChanged);

  // Set GPU device ID
  cudaCheck(cudaSetDevice(deviceID));

  ResizeArray<patch_pair> &patch_pairs(*patch_pairs_ptr);
  ResizeArray<int> &patch_pair_num(*patch_pair_num_ptr);

  if ( workStarted ) { //if work already started, check if finished
    if ( finishWork() ) {  // finished
      workStarted = 0;
      basePriority = PROXY_DATA_PRIORITY;  // higher to aid overlap
    } else {  // need to call again
      workStarted = 2;
      basePriority = PROXY_RESULTS_PRIORITY;  // lower for local
      if ( master == this && kernel_launch_state > 2 ) {
        cuda_check_local_progress(this,0.);  // launches polling
      }
    }
    return;
  }

  workStarted = 1;
  basePriority = COMPUTE_PROXY_PRIORITY;

  Molecule *mol = Node::Object()->molecule;
  Parameters *params = Node::Object()->parameters;
  SimParameters *simParams = Node::Object()->simParameters;

  //execute only during GBIS phase 1, or if not using GBIS
  if (!simParams->GBISOn || gbisPhase == 1) {

    //bind new patches to GPU
    if ( atomsChanged || computesChanged ) {
      int npatches = activePatches.size();

      pairlistsValid = 0;
      pairlistTolerance = 0.;

      if ( computesChanged ) {
        computesChanged = 0;

        if ( ! dummy_size ) {
          dummy_size = 1024*1024;
          cudaMalloc((void**)&dummy_dev,dummy_size);
          cuda_errcheck("in cudaMalloc dummy_dev");
        }

        bool did_realloc = reallocate_host<int>(&force_ready_queue, &force_ready_queue_size, npatches,
          1.2f, cudaHostAllocMapped);
        if (did_realloc) {
          for (int k=0; k < force_ready_queue_size; ++k)
            force_ready_queue[k] = -1;
        }
        force_ready_queue_len = npatches;
        reallocate_host<int>(&block_order, &block_order_size,
          2*(localComputeRecords.size()+remoteComputeRecords.size()),
          1.2f, cudaHostAllocMapped);
    
        num_virials = npatches;
        reallocate_host<float>(&virials, &virials_size, 2*16*num_virials,
          1.2f, cudaHostAllocMapped);
        slow_virials = virials + 16*num_virials;

        reallocate_host<float>(&energy_gbis, &energy_gbis_size, npatches,
          1.2f, cudaHostAllocMapped);
        for (int i  = 0; i < energy_gbis_size; i++) energy_gbis[i] = 0.f;

        int *ap = activePatches.begin();
        for ( int i=0; i<localActivePatches.size(); ++i ) {
          *(ap++) = localActivePatches[i];
        }
        for ( int i=0; i<remoteActivePatches.size(); ++i ) {
          *(ap++) = remoteActivePatches[i];
        }

        // sort computes by distance between patches
        cr_sortop_distance so(lattice);
        std::stable_sort(localComputeRecords.begin(),localComputeRecords.end(),so);
        std::stable_sort(remoteComputeRecords.begin(),remoteComputeRecords.end(),so);

        const bool streaming = ! (deviceCUDA->getNoStreaming() || simParams->GBISOn);

        if ( streaming ) {
          // iout << "Reverse-priority sorting...\n" << endi;
          cr_sortop_reverse_priority sorp(so,patchRecords.begin());
          std::stable_sort(localComputeRecords.begin(),localComputeRecords.end(),sorp);
          std::stable_sort(remoteComputeRecords.begin(),remoteComputeRecords.end(),sorp);
          patchPairsReordered = 0;
          //patchPairsReordered = 1;
      // int len = remoteComputeRecords.size();
      // for ( int i=0; i<len; ++i ) {
      //   iout << "reverse_order " << i << " " << remoteComputeRecords[i].pid[0] << "\n";
      // }
      // int len2 = localComputeRecords.size();
      // for ( int i=0; i<len2; ++i ) {
      //   iout << "reverse_order " << (i+len) << " " << localComputeRecords[i].pid[0] << "\n";
      // }
      // iout << endi;
        } else {
          patchPairsReordered = 1;
        }
 
        int nlc = localComputeRecords.size();
        int nrc = remoteComputeRecords.size();
        computeRecords.resize(nlc+nrc);
        compute_record *cr = computeRecords.begin();
        for ( int i=0; i<nrc; ++i ) {
          *(cr++) = remoteComputeRecords[i];
        }
        for ( int i=0; i<nlc; ++i ) {
          *(cr++) = localComputeRecords[i];
        }

        // patch_pair_num[i] = number of patch pairs that involve patch i
        patch_pair_num.resize(npatches);
        for ( int i=0; i<npatches; ++i ) {
          patchRecords[activePatches[i]].localIndex = i;
          patch_pair_num[i] = 0;
        }

        int ncomputes = computeRecords.size();
        patch_pairs.resize(ncomputes);
        for ( int i=0; i<ncomputes; ++i ) {
          ComputeNonbondedCUDA::compute_record &cr = computeRecords[i];
          int lp1 = patchRecords[cr.pid[0]].localIndex;
          int lp2 = patchRecords[cr.pid[1]].localIndex;
          patch_pair_num[lp1]++;
          if (lp1 != lp2) patch_pair_num[lp2]++;
          patch_pair &pp = patch_pairs[i];
          pp.offset.x = cr.offset.x;
          pp.offset.y = cr.offset.y;
          pp.offset.z = cr.offset.z;
        }

        for ( int i=0; i<ncomputes; ++i ) {
          ComputeNonbondedCUDA::compute_record &cr = computeRecords[i];
          int lp1 = patchRecords[cr.pid[0]].localIndex;
          int lp2 = patchRecords[cr.pid[1]].localIndex;
          patch_pair &pp = patch_pairs[i];
          pp.patch1_ind = lp1;
          pp.patch2_ind = lp2;
          pp.patch1_num_pairs = patch_pair_num[lp1];
          pp.patch2_num_pairs = patch_pair_num[lp2];
        }

      if ( CmiPhysicalNodeID(CkMyPe()) < 2 ) {
        CkPrintf("Pe %d has %d local and %d remote patches and %d local and %d remote computes.\n",
	         CkMyPe(), localActivePatches.size(), remoteActivePatches.size(),
	         localComputeRecords.size(), remoteComputeRecords.size());
      }
    }  // computesChanged
    else if ( ! patchPairsReordered ) {
      patchPairsReordered = 1;
      int len = patch_pairs.size();
      int nlc = localComputeRecords.size();
      int nrc = remoteComputeRecords.size();
      if ( len != nlc + nrc ) NAMD_bug("array size mismatch in ComputeNonbondedCUDA reordering");
      ResizeArray<ComputeNonbondedCUDA::compute_record> new_computeRecords(len);
      ResizeArray<patch_pair> new_patch_pairs(len);
      int irc=nrc;
      int ilc=len;
      for ( int i=0; i<len; ++i ) {
        int boi = block_order[i];
        int dest;
        if ( boi < nrc ) { dest = --irc; } else { dest = --ilc; }
        new_computeRecords[dest] = computeRecords[boi];
        new_patch_pairs[dest] = patch_pairs[boi];
      }
      if ( irc != 0 || ilc != nrc ) NAMD_bug("block index mismatch in ComputeNonbondedCUDA reordering");
      computeRecords.swap(new_computeRecords);
      patch_pairs.swap(new_patch_pairs);
    }

    int istart = 0;
    int i;
    for ( i=0; i<npatches; ++i ) {
      if ( i == localActivePatches.size() ) {
        num_local_atoms = istart;
      }
      patch_record &pr = patchRecords[activePatches[i]];
      pr.localStart = istart;
      int natoms = pr.p->getNumAtoms();
      int nfreeatoms = natoms;
      if ( fixedAtomsOn ) {
        const CompAtomExt *aExt = pr.xExt;
        for ( int j=0; j<natoms; ++j ) {
          if ( aExt[j].atomFixed ) --nfreeatoms;
        }
      }
      pr.numAtoms = natoms;
      pr.numFreeAtoms = nfreeatoms;
      istart += natoms;
      istart += 16 - (natoms & 15);
    }
    if ( i == localActivePatches.size() ) {
      num_local_atoms = istart;
    }
    num_atoms = istart;
    num_remote_atoms = num_atoms - num_local_atoms;
    reallocate_host<atom_param>(&atom_params, &atom_params_size, num_atoms, 1.2f);
    reallocate_host<int>(&vdw_types, &vdw_types_size, num_atoms, 1.2f);
    reallocate_host<CudaAtom>(&atoms, &atoms_size, num_atoms, 1.2f);
    reallocate_host<float4>(&forces, &forces_size, num_atoms, 1.2f, cudaHostAllocMapped);
    reallocate_host<float4>(&slow_forces, &slow_forces_size, num_atoms, 1.2f, cudaHostAllocMapped);
    if (simParams->GBISOn) {
      reallocate_host<float>(&intRad0H, &intRad0H_size, num_atoms, 1.2f);
      reallocate_host<float>(&intRadSH, &intRadSH_size, num_atoms, 1.2f);
      reallocate_host<GBReal>(&psiSumH, &psiSumH_size, num_atoms, 1.2f, cudaHostAllocMapped);
      reallocate_host<float>(&bornRadH, &bornRadH_size, num_atoms, 1.2f);
      reallocate_host<GBReal>(&dEdaSumH, &dEdaSumH_size, num_atoms, 1.2f, cudaHostAllocMapped);
      reallocate_host<float>(&dHdrPrefixH, &dHdrPrefixH_size, num_atoms, 1.2f);
    }

    int bfstart = 0;
    int exclmask_start = 0;
    int ncomputes = computeRecords.size();
    for ( int i=0; i<ncomputes; ++i ) {
      ComputeNonbondedCUDA::compute_record &cr = computeRecords[i];
      int p1 = cr.pid[0];
      int p2 = cr.pid[1];
      patch_pair &pp = patch_pairs[i];
      pp.patch1_start = patchRecords[p1].localStart;
      pp.patch1_size  = patchRecords[p1].numAtoms;
      pp.patch1_free_size = patchRecords[p1].numFreeAtoms;
      pp.patch2_start = patchRecords[p2].localStart;
      pp.patch2_size  = patchRecords[p2].numAtoms;
      pp.patch2_free_size = patchRecords[p2].numFreeAtoms;
      pp.plist_start = bfstart;
      // size1*size2 = number of patch pairs
      int size1 = (pp.patch1_size-1)/WARPSIZE+1;
      int size2 = (pp.patch2_size-1)/WARPSIZE+1;
      pp.plist_size = (size1*size2-1)/32+1;
      bfstart += pp.plist_size;
      pp.exclmask_start = exclmask_start;
      exclmask_start += size1*size2;
    } //for ncomputes

    cuda_bind_patch_pairs(patch_pairs.begin(), patch_pairs.size(),
      activePatches.size(), num_atoms, bfstart, 
      exclmask_start);

  }  // atomsChanged || computesChanged

  double charge_scaling = sqrt(COULOMB * scaling * dielectric_1);

  Flags &flags = patchRecords[hostedPatches[0]].p->flags;
  float maxAtomMovement = 0.;
  float maxPatchTolerance = 0.;

  for ( int i=0; i<activePatches.size(); ++i ) {
    patch_record &pr = patchRecords[activePatches[i]];

    float maxMove = pr.p->flags.maxAtomMovement;
    if ( maxMove > maxAtomMovement ) maxAtomMovement = maxMove;

    float maxTol = pr.p->flags.pairlistTolerance;
    if ( maxTol > maxPatchTolerance ) maxPatchTolerance = maxTol;

    int start = pr.localStart;
    int n = pr.numAtoms;
    const CompAtom *a = pr.x;
    const CompAtomExt *aExt = pr.xExt;
    if ( atomsChanged ) {

      atom_param *ap = atom_params + start;
      for ( int k=0; k<n; ++k ) {
        int j = aExt[k].sortOrder;
        ap[k].vdw_type = a[j].vdwType;
        vdw_types[start + k] = a[j].vdwType;
        ap[k].index = aExt[j].id;
#ifdef MEM_OPT_VERSION
        ap[k].excl_index = exclusionsByAtom[aExt[j].exclId].y;
        ap[k].excl_maxdiff = exclusionsByAtom[aExt[j].exclId].x;
#else // ! MEM_OPT_VERSION
        ap[k].excl_index = exclusionsByAtom[aExt[j].id].y;
        ap[k].excl_maxdiff = exclusionsByAtom[aExt[j].id].x;
#endif // MEM_OPT_VERSION
      }
    }
    {
#if 1
      const CudaAtom *ac = pr.p->getCudaAtomList();
      CudaAtom *ap = atoms + start;
      memcpy(ap, ac, sizeof(atom)*n);
#else
      Vector center =
        pr.p->flags.lattice.unscale(cudaCompute->patchMap->center(pr.patchID));
      atom *ap = atoms + start;
      for ( int k=0; k<n; ++k ) {
        int j = aExt[k].sortOrder;
        ap[k].position.x = a[j].position.x - center.x;
        ap[k].position.y = a[j].position.y - center.y;
        ap[k].position.z = a[j].position.z - center.z;
        ap[k].charge = charge_scaling * a[j].charge;
      }
#endif
    }
  }

  savePairlists = 0;
  usePairlists = 0;
  if ( flags.savePairlists ) {
    savePairlists = 1;
    usePairlists = 1;
  } else if ( flags.usePairlists ) {
    if ( ! pairlistsValid ||
         ( 2. * maxAtomMovement > pairlistTolerance ) ) {
      reduction->item(REDUCTION_PAIRLIST_WARNINGS) += 1;
    } else {
      usePairlists = 1;
    }
  }
  if ( ! usePairlists ) {
    pairlistsValid = 0;
  }
  float plcutoff = cutoff;
  if ( savePairlists ) {
    pairlistsValid = 1;
    pairlistTolerance = 2. * maxPatchTolerance;
    plcutoff += pairlistTolerance;
  }
  plcutoff2 = plcutoff * plcutoff;

  //CkPrintf("plcutoff = %f  listTolerance = %f  save = %d  use = %d\n",
  //  plcutoff, pairlistTolerance, savePairlists, usePairlists);

#if 0
  // calculate warp divergence
  if ( 1 ) {
    Flags &flags = patchRecords[hostedPatches[0]].p->flags;
    Lattice &lattice = flags.lattice;
    float3 lata, latb, latc;
    lata.x = lattice.a().x;
    lata.y = lattice.a().y;
    lata.z = lattice.a().z;
    latb.x = lattice.b().x;
    latb.y = lattice.b().y;
    latb.z = lattice.b().z;
    latc.x = lattice.c().x;
    latc.y = lattice.c().y;
    latc.z = lattice.c().z;

    int ncomputes = computeRecords.size();
    for ( int ic=0; ic<ncomputes; ++ic ) {
      patch_pair &pp = patch_pairs[ic];
      atom *a1 = atoms + pp.patch1_atom_start;
      int n1 = pp.patch1_size;
      atom *a2 = atoms + pp.patch2_atom_start;
      int n2 = pp.patch2_size;
      float offx = pp.offset.x * lata.x
               + pp.offset.y * latb.x
               + pp.offset.z * latc.x;
      float offy = pp.offset.x * lata.y
               + pp.offset.y * latb.y
               + pp.offset.z * latc.y;
      float offz = pp.offset.x * lata.z
               + pp.offset.y * latb.z
               + pp.offset.z * latc.z;
      // CkPrintf("%f %f %f\n", offx, offy, offz);
      int atoms_tried = 0;
      int blocks_tried = 0;
      int atoms_used = 0;
      int blocks_used = 0;
      for ( int ii=0; ii<n1; ii+=32 ) {  // warps
        for ( int jj=0; jj<n2; jj+=16 ) {  // shared atom loads
          int block_used = 0;
          for ( int j=jj; j<jj+16 && j<n2; ++j ) {  // shared atoms
            int atom_used = 0;
            for ( int i=ii; i<ii+32 && i<n1; ++i ) {  // threads
              float dx = offx + a1[i].position.x - a2[j].position.x;
              float dy = offy + a1[i].position.y - a2[j].position.y;
              float dz = offz + a1[i].position.z - a2[j].position.z;
              float r2 = dx*dx + dy*dy + dz*dz;
              if ( r2 < cutoff2 ) atom_used = 1;
            }
            ++atoms_tried;
            if ( atom_used ) { block_used = 1; ++atoms_used; }
          }
          ++blocks_tried;
          if ( block_used ) { ++blocks_used; }
        }
      }
      CkPrintf("blocks = %d/%d (%f)  atoms = %d/%d (%f)\n",
                blocks_used, blocks_tried, blocks_used/(float)blocks_tried,
                atoms_used, atoms_tried, atoms_used/(float)atoms_tried);
    }
  }
#endif

  } // !GBISOn || gbisPhase == 1

  //Do GBIS
  if (simParams->GBISOn) {
    //open GBIS boxes depending on phase
    for ( int i=0; i<activePatches.size(); ++i ) {
      patch_record &pr = master->patchRecords[activePatches[i]];
      GBISP("doWork[%d] accessing arrays for P%d\n",CkMyPe(),gbisPhase);
      if (gbisPhase == 1) {
        //Copy GBIS intRadius to Host
        if (atomsChanged) {
          float *intRad0 = intRad0H + pr.localStart;
          float *intRadS = intRadSH + pr.localStart;
          for ( int k=0; k<pr.numAtoms; ++k ) {
            int j = pr.xExt[k].sortOrder;
            intRad0[k] = pr.intRad[2*j+0];
            intRadS[k] = pr.intRad[2*j+1];
          }
        }
      } else if (gbisPhase == 2) {
        float *bornRad = bornRadH + pr.localStart;
        for ( int k=0; k<pr.numAtoms; ++k ) {
          int j = pr.xExt[k].sortOrder;
          bornRad[k] = pr.bornRad[j];
        }
      } else if (gbisPhase == 3) {
        float *dHdrPrefix = dHdrPrefixH + pr.localStart;
        for ( int k=0; k<pr.numAtoms; ++k ) {
          int j = pr.xExt[k].sortOrder;
          dHdrPrefix[k] = pr.dHdrPrefix[j];
        }
      } // end phases
    } // end for patches
  } // if GBISOn

  kernel_time = CkWallTimer();
  kernel_launch_state = 1;
  //if ( gpu_is_mine || ! doSlow ) recvYieldDevice(-1);
  if ( deviceCUDA->getGpuIsMine() || ! doSlow ) recvYieldDevice(-1);
}

void cuda_check_remote_calc(void *arg, double walltime) {
  // in theory we only need end_remote_calc, but overlap isn't reliable
  // if ( cudaEventQuery(end_remote_calc) == cudaSuccess ) {
  if ( cudaEventQuery(end_remote_download) == cudaSuccess ) {
// CkPrintf("Pe %d yielding to %d after remote calc\n", CkMyPe(), next_pe_sharing_gpu);
    computeMgr->sendYieldDevice(deviceCUDA->getNextPeSharingGpu());
// CkPrintf("Pe %d yielded to %d after remote calc\n", CkMyPe(), next_pe_sharing_gpu);
  } else {
    CcdCallBacksReset(0,walltime);  // fix Charm++
    CUDA_POLL(cuda_check_remote_calc, arg);
  }
}

void cuda_check_local_calc(void *arg, double walltime) {
  // in theory we only need end_local_calc, but overlap isn't reliable
  // if ( cudaEventQuery(end_local_calc) == cudaSuccess ) {
  if ( cudaEventQuery(end_local_download) == cudaSuccess ) {
// CkPrintf("Pe %d yielding to %d after local calc\n", CkMyPe(), next_pe_sharing_gpu);
    computeMgr->sendYieldDevice(deviceCUDA->getNextPeSharingGpu());
// CkPrintf("Pe %d yielded to %d after local calc\n", CkMyPe(), next_pe_sharing_gpu);
  } else {
    CcdCallBacksReset(0,walltime);  // fix Charm++
    CUDA_POLL(cuda_check_local_calc, arg);
  }
}

// computeMgr->sendYieldDevice(next_pe_sharing_gpu);

void ComputeNonbondedCUDA::recvYieldDevice(int pe) {
GBISP("C.N.CUDA[%d]::recvYieldDevice: seq %d, workStarted %d, \
gbisPhase %d, kls %d, from pe %d\n", CkMyPe(), sequence(), \
workStarted, gbisPhase, kernel_launch_state, pe)

  float3 lata, latb, latc;
  lata.x = lattice.a().x;
  lata.y = lattice.a().y;
  lata.z = lattice.a().z;
  latb.x = lattice.b().x;
  latb.y = lattice.b().y;
  latb.z = lattice.b().z;
  latc.x = lattice.c().x;
  latc.y = lattice.c().y;
  latc.z = lattice.c().z;
  SimParameters *simParams = Node::Object()->simParameters;

  // Set GPU device ID
  cudaSetDevice(deviceID);

  const bool streaming = ! (deviceCUDA->getNoStreaming() || simParams->GBISOn);

  double walltime;
  if ( kernel_launch_state == 1 || kernel_launch_state == 2 ) {
    walltime = CkWallTimer();
    CcdCallBacksReset(0,walltime);  // fix Charm++
  }

  switch ( kernel_launch_state ) {
////////////////////////////////////////////////////////////
// Remote
  case 1:
GBISP("C.N.CUDA[%d]::recvYieldDeviceR: case 1\n", CkMyPe())
    ++kernel_launch_state;
    //gpu_is_mine = 0;
    deviceCUDA->setGpuIsMine(0);
    remote_submit_time = walltime;

    if (!simParams->GBISOn || gbisPhase == 1) {
      // cudaEventRecord(start_upload, stream);
      if ( atomsChanged ) {
        cuda_bind_atom_params(atom_params);
        cuda_bind_vdw_types(vdw_types);
      }
      if ( simParams->GBISOn) {
        cuda_bind_GBIS_psiSum(psiSumH);
	       if ( atomsChanged ) {
	         cuda_bind_GBIS_intRad(intRad0H, intRadSH);
	       }
      }
      atomsChanged = 0;
      cuda_bind_atoms((const atom *)atoms);
      cuda_bind_forces(forces, slow_forces);
      cuda_bind_virials(virials, force_ready_queue, block_order);
      if ( simParams->GBISOn) {
        cuda_bind_GBIS_energy(energy_gbis);
      }
      if ( stream2 != stream ) cudaEventRecord(start_calc, stream);
      //call CUDA Kernels

      cuda_nonbonded_forces(lata, latb, latc, cutoff2, plcutoff2,
			    0,remoteComputeRecords.size(),
			    remoteComputeRecords.size()+localComputeRecords.size(),
			    doSlow, doEnergy, usePairlists, savePairlists, streaming, ! patchPairsReordered, stream);
      if (simParams->GBISOn) {
        cuda_GBIS_P1(
          0,remoteComputeRecords.size(),
		      localActivePatches.size(),remoteActivePatches.size(),
		      simParams->alpha_cutoff-simParams->fsMax,
		      simParams->coulomb_radius_offset,
		      lata, latb, latc, stream);
      }
      //cuda_load_forces(forces, (doSlow ? slow_forces : 0 ),
      //    num_local_atom_records,num_remote_atom_records);
      if ( ( ! streaming ) || ( deviceCUDA->getSharedGpu() && ! deviceCUDA->getMergeGrids() ) ) {
        cudaEventRecord(end_remote_download, stream);
      }
      if ( streaming ) {
        force_ready_queue_next = 0;
        CUDA_POLL(cuda_check_progress,this);
      } else {
        CUDA_POLL(cuda_check_remote_progress,this);
      }
      if ( deviceCUDA->getSharedGpu() && ! deviceCUDA->getMergeGrids() ) {
        CUDA_POLL(cuda_check_remote_calc,this);
        break;
      }
    } // !GBIS or gbisPhase==1
    if (simParams->GBISOn) {
      if (gbisPhase == 1) {
        //GBIS P1 Kernel launched in previous code block
      } else if (gbisPhase == 2) {
GBISP("C.N.CUDA[%d]::recvYieldDeviceR: <<<P2>>>\n", CkMyPe())
        // cudaEventRecord(start_upload, stream);
        cuda_bind_GBIS_bornRad(bornRadH);
        cuda_bind_GBIS_dEdaSum(dEdaSumH);
        if ( stream2 != stream ) cudaEventRecord(start_calc, stream);
        cuda_GBIS_P2(
          0,remoteComputeRecords.size(),
          localActivePatches.size(),remoteActivePatches.size(),
          (simParams->alpha_cutoff-simParams->fsMax), simParams->cutoff,
          simParams->nonbondedScaling, simParams->kappa,
          (simParams->switchingActive ? simParams->switchingDist : -1.0),
          simParams->dielectric, simParams->solvent_dielectric,
          lata, latb, latc,
          doEnergy, doSlow, stream
          );
        cudaEventRecord(end_remote_download, stream);
        CUDA_POLL(cuda_check_remote_progress,this);
      } else if (gbisPhase == 3) {
GBISP("C.N.CUDA[%d]::recvYieldDeviceR: <<<P3>>>\n", CkMyPe())
        // cudaEventRecord(start_upload, stream);
        cuda_bind_GBIS_dHdrPrefix(dHdrPrefixH);
        if ( stream2 != stream ) cudaEventRecord(start_calc, stream);
        if (doSlow)
        cuda_GBIS_P3(
          0,remoteComputeRecords.size(),
          localActivePatches.size(),remoteActivePatches.size(),
          (simParams->alpha_cutoff-simParams->fsMax),
          simParams->coulomb_radius_offset,
          simParams->nonbondedScaling,
          lata, latb, latc, stream
          );
        cudaEventRecord(end_remote_download, stream);
        CUDA_POLL(cuda_check_remote_progress,this);
      }
    }

////////////////////////////////////////////////////////////
// Local
  case 2:
GBISP("C.N.CUDA[%d]::recvYieldDeviceL: case 2\n", CkMyPe())
    ++kernel_launch_state;
    //gpu_is_mine = 0;
    deviceCUDA->setGpuIsMine(0);

    if ( stream2 != stream ) {
      // needed to ensure that upload is finished
      cudaStreamWaitEvent(stream2, start_calc, 0);
      // priorities do not prevent local from launching ahead
      // of remote, so delay local stream with a small memset
      cudaMemsetAsync(dummy_dev, 0, dummy_size, stream2);
    }

    if (!simParams->GBISOn || gbisPhase == 1) {

      cuda_nonbonded_forces(lata, latb, latc, cutoff2, plcutoff2,
			    remoteComputeRecords.size(),localComputeRecords.size(),
			    remoteComputeRecords.size()+localComputeRecords.size(),
			    doSlow, doEnergy, usePairlists, savePairlists, streaming, ! patchPairsReordered, stream2);
      if (simParams->GBISOn) {
        cuda_GBIS_P1(
		     remoteComputeRecords.size(),localComputeRecords.size(),
		     0,localActivePatches.size(),
		     simParams->alpha_cutoff-simParams->fsMax,
		     simParams->coulomb_radius_offset,
		     lata, latb, latc, stream2 );
      }
      //cuda_load_forces(forces, (doSlow ? slow_forces : 0 ),
      //    0,num_local_atom_records);
      //cuda_load_virials(virials, doSlow);  // slow_virials follows virials
      if ( ( ! streaming ) || ( deviceCUDA->getSharedGpu() && ! deviceCUDA->getMergeGrids() ) ) {
        cudaEventRecord(end_local_download, stream2);
      }
      if ( ! streaming && workStarted == 2 ) {
        GBISP("C.N.CUDA[%d]::recvYieldDeviceL: adding POLL \
cuda_check_local_progress\n", CkMyPe())
        CUDA_POLL(cuda_check_local_progress,this);
      }
      if ( deviceCUDA->getSharedGpu() && ! deviceCUDA->getMergeGrids() ) {
	GBISP("C.N.CUDA[%d]::recvYieldDeviceL: adding POLL \
cuda_check_local_calc\n", CkMyPe())
	  CUDA_POLL(cuda_check_local_calc,this);
	break;
      }

    } // !GBIS or gbisPhase==1
    if (simParams->GBISOn) {
      if (gbisPhase == 1) {
        //GBIS P1 Kernel launched in previous code block
      } else if (gbisPhase == 2) {
GBISP("C.N.CUDA[%d]::recvYieldDeviceL: calling <<<P2>>>\n", CkMyPe())
        cuda_GBIS_P2(
          remoteComputeRecords.size(),localComputeRecords.size(),
          0,localActivePatches.size(),
          (simParams->alpha_cutoff-simParams->fsMax), simParams->cutoff,
          simParams->nonbondedScaling, simParams->kappa,
          (simParams->switchingActive ? simParams->switchingDist : -1.0),
          simParams->dielectric, simParams->solvent_dielectric,
          lata, latb, latc,
          doEnergy, doSlow, stream2
          );
        cudaEventRecord(end_local_download, stream2);
        if ( workStarted == 2 ) {
          CUDA_POLL(cuda_check_local_progress,this);
        }
      } else if (gbisPhase == 3) {
GBISP("C.N.CUDA[%d]::recvYieldDeviceL: calling <<<P3>>>\n", CkMyPe())
        if (doSlow)
        cuda_GBIS_P3(
          remoteComputeRecords.size(),localComputeRecords.size(),
          0,localActivePatches.size(),
          (simParams->alpha_cutoff-simParams->fsMax),
          simParams->coulomb_radius_offset,
          simParams->nonbondedScaling,
          lata, latb, latc, stream2
          );
        cudaEventRecord(end_local_download, stream2);
        if ( workStarted == 2 ) {
          CUDA_POLL(cuda_check_local_progress,this);
        }
      } // phases
    } // GBISOn
    if ( simParams->PMEOn && simParams->PMEOffload  && !simParams->usePMECUDA) break;

  default:
GBISP("C.N.CUDA[%d]::recvYieldDevice: case default\n", CkMyPe())
    //gpu_is_mine = 1;
    deviceCUDA->setGpuIsMine(1);
    break;
  } // switch
GBISP("C.N.CUDA[%d]::recvYieldDevice: DONE\n", CkMyPe())
}

void ComputeNonbondedCUDA::messageFinishPatch(int flindex) {
  int pid = activePatches[flindex];
  patch_record &pr = patchRecords[pid];
  //fprintf(stderr, "%d ComputeNonbondedCUDA::messageFinishPatch %d\n",CkMyPe(),pr.hostPe);
  computeMgr->sendNonbondedCUDASlaveEnqueuePatch(pr.slave,pr.hostPe,sequence(),PROXY_DATA_PRIORITY,flindex,pr.msg);
}

void ComputeNonbondedCUDA::finishPatch(int flindex) {
  //fprintf(stderr, "%d ComputeNonbondedCUDA::finishPatch \n",CkMyPe());
  patch_record &pr = master->patchRecords[master->activePatches[flindex]];
  pr.r = pr.forceBox->open();
  finishPatch(pr);
  pr.positionBox->close(&(pr.x));
  pr.forceBox->close(&(pr.r));
}

void ComputeNonbondedCUDA::finishPatch(patch_record &pr) {
  int start = pr.localStart;
  const CompAtomExt *aExt = pr.xExt;
  int nfree = pr.numAtoms;
  pr.f = pr.r->f[Results::nbond];
  Force *f = pr.f;
  Force *f_slow = pr.r->f[Results::slow];
  const CompAtom *a = pr.x;
  // int *ao = atom_order + start;
  float4 *af = master->forces + start;
  float4 *af_slow = master->slow_forces + start;
  // only free atoms return forces
  for ( int k=0; k<nfree; ++k ) {
    int j = aExt[k].sortOrder;
    f[j].x += af[k].x;
    f[j].y += af[k].y;
    f[j].z += af[k].z;
    // wcount += af[k].w;
    // virial += af[k].w;
    if ( doSlow ) {
      f_slow[j].x += af_slow[k].x;
      f_slow[j].y += af_slow[k].y;
      f_slow[j].z += af_slow[k].z;
      // virial_slow += af_slow[k].w;
    }
  }
}

//dtanner
int ComputeNonbondedCUDA::finishWork() {
  //fprintf(stderr, "%d ComputeNonbondedCUDA::finishWork() \n",CkMyPe());

  if ( master == this ) {
    for ( int i = 0; i < numSlaves; ++i ) {
      computeMgr->sendNonbondedCUDASlaveEnqueue(slaves[i],slavePes[i],sequence(),priority(),workStarted);
    }
  }

GBISP("C.N.CUDA[%d]::fnWork: workStarted %d, phase %d\n", \
CkMyPe(), workStarted, gbisPhase)

  Molecule *mol = Node::Object()->molecule;
  SimParameters *simParams = Node::Object()->simParameters;

  ResizeArray<int> &patches( workStarted == 1 ?
				remoteHostedPatches : localHostedPatches );

  // long long int wcount = 0;
  // double virial = 0.;
  // double virial_slow = 0.;
  for ( int i=0; i<patches.size(); ++i ) {
    // CkPrintf("Pe %d patch %d of %d pid %d\n",CkMyPe(),i,patches.size(),patches[i]);
    patch_record &pr = master->patchRecords[patches[i]];

    if ( !simParams->GBISOn || gbisPhase == 1 ) {
      patch_record &pr = master->patchRecords[patches[i]];
GBISP("GBIS[%d] fnWork() P0[%d] force.open()\n",CkMyPe(), pr.patchID);
      pr.r = pr.forceBox->open();
    } // !GBISOn || gbisPhase==1

    int start = pr.localStart;
    const CompAtomExt *aExt = pr.xExt;
    if ( !simParams->GBISOn || gbisPhase == 3 ) {
      finishPatch(pr);
    } // !GBISOn || gbisPhase == 3

#if 0
    if ( i % 31 == 0 ) for ( int j=0; j<3; ++j ) {
      CkPrintf("Pe %d patch %d atom %d (%f %f %f) force %f\n", CkMyPe(), i,
	j, pr.x[j].position.x, pr.x[j].position.y, pr.x[j].position.z,
	af[j].w);
    }
#endif

    //Close Boxes depending on Phase
    if (simParams->GBISOn) {
      if (gbisPhase == 1) {
        //Copy dEdaSum from Host to Patch Box
        GBReal *psiSumMaster = master->psiSumH + start;
        for ( int k=0; k<pr.numAtoms; ++k ) {
          int j = aExt[k].sortOrder;
          pr.psiSum[j] += psiSumMaster[k];
        }
GBISP("C.N.CUDA[%d]::fnWork: P1 psiSum.close()\n", CkMyPe());
        pr.psiSumBox->close(&(pr.psiSum));

      } else if (gbisPhase == 2) {
        //Copy dEdaSum from Host to Patch Box
        GBReal *dEdaSumMaster = master->dEdaSumH + start;
        for ( int k=0; k<pr.numAtoms; ++k ) {
          int j = aExt[k].sortOrder;
          pr.dEdaSum[j] += dEdaSumMaster[k];
        }
GBISP("C.N.CUDA[%d]::fnWork: P2 dEdaSum.close()\n", CkMyPe());
        pr.dEdaSumBox->close(&(pr.dEdaSum));

      } else if (gbisPhase == 3) {
GBISP("C.N.CUDA[%d]::fnWork: P3 all.close()\n", CkMyPe());
        pr.intRadBox->close(&(pr.intRad)); //box 6
        pr.bornRadBox->close(&(pr.bornRad)); //box 7
        pr.dHdrPrefixBox->close(&(pr.dHdrPrefix)); //box 9
        pr.positionBox->close(&(pr.x)); //box 0
        pr.forceBox->close(&(pr.r));
      } //end phases
    } else { //not GBIS
GBISP("C.N.CUDA[%d]::fnWork: pos/force.close()\n", CkMyPe());
      pr.positionBox->close(&(pr.x));
      pr.forceBox->close(&(pr.r));
    }
  }// end for

#if 0
  virial *= (-1./6.);
  reduction->item(REDUCTION_VIRIAL_NBOND_XX) += virial;
  reduction->item(REDUCTION_VIRIAL_NBOND_YY) += virial;
  reduction->item(REDUCTION_VIRIAL_NBOND_ZZ) += virial;
  if ( doSlow ) {
    virial_slow *= (-1./6.);
    reduction->item(REDUCTION_VIRIAL_SLOW_XX) += virial_slow;
    reduction->item(REDUCTION_VIRIAL_SLOW_YY) += virial_slow;
    reduction->item(REDUCTION_VIRIAL_SLOW_ZZ) += virial_slow;
  }
#endif

  if ( workStarted == 1 && ! deviceCUDA->getMergeGrids() &&
       ( localHostedPatches.size() || master == this ) ) {
    GBISP("not finished, call again\n");
    return 0;  // not finished, call again
  }

  if ( master != this ) {  // finished
    GBISP("finished\n");
    if (simParams->GBISOn) gbisPhase = 1 + (gbisPhase % 3);//1->2->3->1...
    atomsChanged = 0;
    return 1;
  }

  cuda_timer_total += kernel_time;

  if ( !simParams->GBISOn || gbisPhase == 3 ) {

    atomsChanged = 0;
    finishReductions();

  } // !GBISOn || gbisPhase==3  

  // Next GBIS Phase
GBISP("C.N.CUDA[%d]::fnWork: incrementing phase\n", CkMyPe())
    if (simParams->GBISOn) gbisPhase = 1 + (gbisPhase % 3);//1->2->3->1...

  GBISP("C.N.CUDA[%d] finished ready for next step\n",CkMyPe());
  return 1;  // finished and ready for next step
}


void ComputeNonbondedCUDA::finishReductions() {
  //fprintf(stderr, "%d ComputeNonbondedCUDA::finishReductions \n",CkMyPe());

  basePriority = PROXY_DATA_PRIORITY;  // higher to aid overlap

  SimParameters *simParams = Node::Object()->simParameters;

  if ( 0 && patchPairsReordered && patchPairsReordered < 3 ) {
    if ( patchPairsReordered++ == 2) {
      int patch_len = patchRecords.size();
      ResizeArray<int> plast(patch_len);
      for ( int i=0; i<patch_len; ++i ) {
        plast[i] = -1;
      }
      int order_len = localComputeRecords.size()+remoteComputeRecords.size();
      for ( int i=0; i<order_len; ++i ) {
        plast[computeRecords[block_order[i]].pid[0]] = i;
        iout << "block_order " << i << " " << block_order[i] << " " << computeRecords[block_order[i]].pid[0] << "\n";
      }
      iout << endi;
      for ( int i=0; i<patch_len; ++i ) {
        iout << "patch_last " << i << " " << plast[i] << "\n";
      }
      iout << endi;
    }
  }

  {
    Tensor virial_tensor;
    BigReal energyv = 0.;
    BigReal energye = 0.;
    BigReal energys = 0.;
    int nexcluded = 0;
    for ( int i = 0; i < num_virials; ++i ) {
      virial_tensor.xx += virials[16*i];
      virial_tensor.xy += virials[16*i+1];
      virial_tensor.xz += virials[16*i+2];
      virial_tensor.yx += virials[16*i+3];
      virial_tensor.yy += virials[16*i+4];
      virial_tensor.yz += virials[16*i+5];
      virial_tensor.zx += virials[16*i+6];
      virial_tensor.zy += virials[16*i+7];
      virial_tensor.zz += virials[16*i+8];
      energyv += virials[16*i+9];
      energye += virials[16*i+10];
      energys += virials[16*i+11];
      nexcluded += ((int *)virials)[16*i+12];
      if (simParams->GBISOn) {
        energye += energy_gbis[i];
      }
    }
    reduction->item(REDUCTION_EXCLUSION_CHECKSUM_CUDA) += nexcluded;
    ADD_TENSOR_OBJECT(reduction,REDUCTION_VIRIAL_NBOND,virial_tensor);
    if ( doEnergy ) {
      reduction->item(REDUCTION_LJ_ENERGY) += energyv;
      reduction->item(REDUCTION_ELECT_ENERGY) += energye;
      reduction->item(REDUCTION_ELECT_ENERGY_SLOW) += energys;
    }
  }
  if ( doSlow ) {
    Tensor virial_slow_tensor;
    for ( int i = 0; i < num_virials; ++i ) {
      virial_slow_tensor.xx += slow_virials[16*i];
      virial_slow_tensor.xy += slow_virials[16*i+1];
      virial_slow_tensor.xz += slow_virials[16*i+2];
      virial_slow_tensor.yx += slow_virials[16*i+3];
      virial_slow_tensor.yy += slow_virials[16*i+4];
      virial_slow_tensor.yz += slow_virials[16*i+5];
      virial_slow_tensor.zx += slow_virials[16*i+6];
      virial_slow_tensor.zy += slow_virials[16*i+7];
      virial_slow_tensor.zz += slow_virials[16*i+8];
    }
    ADD_TENSOR_OBJECT(reduction,REDUCTION_VIRIAL_SLOW,virial_slow_tensor);
  }

  reduction->submit();

  cuda_timer_count++;
  if ( simParams->outputCudaTiming &&
	step % simParams->outputCudaTiming == 0 ) {

    // int natoms = mol->numAtoms; 
    // double wpa = wcount;  wpa /= natoms;

    // CkPrintf("Pe %d CUDA kernel %f ms, total %f ms, wpa %f\n", CkMyPe(),
    // 	kernel_time * 1.e3, time * 1.e3, wpa);

#if 0
    float upload_ms, remote_calc_ms;
    float local_calc_ms, total_ms;
    cuda_errcheck("before event timers");
    cudaEventElapsedTime(&upload_ms, start_upload, start_calc);
    cuda_errcheck("in event timer 1");
    cudaEventElapsedTime(&remote_calc_ms, start_calc, end_remote_download);
    cuda_errcheck("in event timer 2");
    cudaEventElapsedTime(&local_calc_ms, end_remote_download, end_local_download);
    cuda_errcheck("in event timer 3");
    cudaEventElapsedTime(&total_ms, start_upload, end_local_download);
    cuda_errcheck("in event timer 4");
    cuda_errcheck("in event timers");

    CkPrintf("CUDA EVENT TIMING: %d %f %f %f %f\n",
	     CkMyPe(), upload_ms, remote_calc_ms,
	     local_calc_ms, total_ms);
#endif

    if ( cuda_timer_count >= simParams->outputCudaTiming ) {
      cuda_timer_total /= cuda_timer_count;
      CkPrintf("CUDA TIMING: %d  %f ms/step on node %d\n",
	       step, cuda_timer_total * 1.e3, CkMyPe());
    }
    cuda_timer_count = 0;
    cuda_timer_total = 0;
  }

}


#endif  // NAMD_CUDA

