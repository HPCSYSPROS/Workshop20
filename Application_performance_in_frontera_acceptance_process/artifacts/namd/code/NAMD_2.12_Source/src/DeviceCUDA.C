
#include "common.h"
#include "charm++.h"
#include "DeviceCUDA.h"
#include "WorkDistrib.h"
#include "CudaUtils.h"

#ifdef NAMD_CUDA

#include <cuda_runtime.h>
#include <cuda.h>

#ifdef WIN32
#define __thread __declspec(thread)
#endif

// Global storage for CUDA devices
__thread DeviceCUDA *deviceCUDA;

void cuda_initialize() {
	deviceCUDA = new DeviceCUDA();
	deviceCUDA->initialize();
}

// kill all service threads
void cuda_finalize() {
    int ndevs = 0;
    cudaGetDeviceCount(&ndevs);
    for ( int dev=0; dev < ndevs; ++dev ) {
        cudaSetDevice(dev);
        cudaDeviceReset();
    }
}

// -------------------------------------------------------------------------------------------------
// Called from BackEnd.C by all processes to read command line arguments
// These argument settings are used by DeviceCUDA -class
// -------------------------------------------------------------------------------------------------
struct cuda_args_t {
	char *devicelist;
	int usedevicelist;
  int devicesperreplica;
	int ignoresharing;
	int mergegrids;
	int nomergegrids;
	int nostreaming;
};

static __thread cuda_args_t cuda_args;

void cuda_getargs(char **argv) {
  cuda_args.devicelist = 0;
  cuda_args.usedevicelist = CmiGetArgStringDesc(argv, "+devices", &cuda_args.devicelist,
		"comma-delimited list of CUDA device numbers such as 0,2,1,2");
  cuda_args.devicesperreplica = 0;
  CmiGetArgInt(argv, "+devicesperreplica", &cuda_args.devicesperreplica);
  if ( cuda_args.devicesperreplica < 0 ) NAMD_die("Devices per replica must be positive\n");
  cuda_args.ignoresharing = CmiGetArgFlag(argv, "+ignoresharing");
  cuda_args.mergegrids = CmiGetArgFlag(argv, "+mergegrids");
  cuda_args.nomergegrids = CmiGetArgFlag(argv, "+nomergegrids");
  if ( cuda_args.mergegrids && cuda_args.nomergegrids ) NAMD_die("Do not specify both +mergegrids and +nomergegrids");
  cuda_args.nostreaming = CmiGetArgFlag(argv, "+nostreaming");
}
// -------------------------------------------------------------------------------------------------

// Node-wide list of device IDs for every rank
#define MAX_NUM_RANKS 2048
int deviceIDList[MAX_NUM_RANKS];
// Node-wide of master PEs for every device ID
#define MAX_NUM_DEVICES 256
int masterPeList[MAX_NUM_DEVICES];

// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------

//
// Class creator
//
DeviceCUDA::DeviceCUDA() : deviceProps(NULL), devices(NULL) {}

//
// Initalize device
//
void DeviceCUDA::initialize() {
	// Copy command-line arguments into class
	this->devicelist = cuda_args.devicelist;
	this->usedevicelist = cuda_args.usedevicelist;
  this->devicesperreplica = cuda_args.devicesperreplica;
	this->ignoresharing = cuda_args.ignoresharing;
	this->mergegrids = cuda_args.mergegrids;
	this->nomergegrids = cuda_args.nomergegrids;
	this->nostreaming = cuda_args.nostreaming;

  if (CkMyPe() == 0) register_user_events();

  if (CkMyPe() == 0) CkPrintf("Info: Built with CUDA version %d\n", CUDA_VERSION);

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

  deviceCount = 0;
  cudaCheck(cudaGetDeviceCount(&deviceCount));
  if ( deviceCount <= 0 ) {
    cudaDie("No CUDA devices found.");
  }

  // Store all device props
  deviceProps = new cudaDeviceProp[deviceCount];
  for ( int i=0; i<deviceCount; ++i ) {
    cudaCheck(cudaGetDeviceProperties(&deviceProps[i], i));
  }

  ndevices = 0;
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
      cudaCheck(cudaGetDeviceProperties(&deviceProp, dev));
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
    cudaDie("all devices are in prohibited mode, of compute capability < 3.0, unable to map host memory, too small, or otherwise unusable");
  }

  if ( devicesperreplica > 0 ) {
    if ( devicesperreplica > ndevices ) {
      NAMD_die("More devices per partition requested than devices are available");
    }
    int *olddevices = devices;
    devices = new int[devicesperreplica];
    for ( int i=0; i<devicesperreplica; ++i ) {
      int mypart = CmiMyPartition();
      devices[i] = olddevices[(i+devicesperreplica*mypart)%ndevices];
    }
    ndevices = devicesperreplica;
    delete [] olddevices;
  }

  sharedGpu = 0;
  gpuIsMine = 1;
  int firstPeSharingGpu = CkMyPe();
  nextPeSharingGpu = CkMyPe();

 {

    int dev;
    if ( numPesOnPhysicalNode > 1 ) {
      int myDeviceRank = myRankInPhysicalNode * ndevices / numPesOnPhysicalNode;
      dev = devices[myDeviceRank];
      masterPe = CkMyPe();
      if ( ignoresharing ) {
        pesSharingDevice = new int[1];
        pesSharingDevice[0] = CkMyPe();
        numPesSharingDevice = 1;
      } else {
        pesSharingDevice = new int[numPesOnPhysicalNode];
        masterPe = -1;
        numPesSharingDevice = 0;
        for ( int i = 0; i < numPesOnPhysicalNode; ++i ) {
          if ( i * ndevices / numPesOnPhysicalNode == myDeviceRank ) {
            int thisPe = pesOnPhysicalNode[i];
            pesSharingDevice[numPesSharingDevice++] = thisPe;
            if ( masterPe < 1 ) masterPe = thisPe;
            if ( WorkDistrib::pe_sortop_diffuse()(thisPe,masterPe) ) masterPe = thisPe;
          }
        }
        for ( int j = 0; j < ndevices; ++j ) {
          if ( devices[j] == dev && j != myDeviceRank ) sharedGpu = 1;
        }
      }
      if ( sharedGpu && masterPe == CkMyPe() ) {
        if ( CmiPhysicalNodeID(masterPe) < 2 )
        CkPrintf("Pe %d sharing CUDA device %d\n", CkMyPe(), dev);
      }
    } else {  // in case phys node code is lying
      dev = devices[CkMyPe() % ndevices];
      masterPe = CkMyPe();
      pesSharingDevice = new int[1];
      pesSharingDevice[0] = CkMyPe();
      numPesSharingDevice = 1;
    }

    deviceID = dev;

    // Store device IDs to node-wide list
    if (CkMyRank() >= MAX_NUM_RANKS)
      NAMD_die("Maximum number of ranks (2048) per node exceeded");
    deviceIDList[CkMyRank()] = deviceID;

    if ( masterPe != CkMyPe() ) {
      if ( CmiPhysicalNodeID(masterPe) < 2 )
      CkPrintf("Pe %d physical rank %d will use CUDA device of pe %d\n",
               CkMyPe(), myRankInPhysicalNode, masterPe);
      // for PME only
      cudaCheck(cudaSetDevice(dev));
      return;
    }

    // Store master PEs for every device ID to node-wide list
    if (CkMyRank() >= MAX_NUM_DEVICES)
      NAMD_die("Maximum number of CUDA devices (256) per node exceeded");
    masterPeList[deviceID] = masterPe;
    // Set masterPe values to -1 for devices that do not exist.
    // Only master Pe with deviceID == devices[0] does the writing
    if (deviceID == devices[0]) {
      // For device IDs 0...deviceCount-1, check if it is in the devices[0...deviceCount-1]
      for (int i=0;i < deviceCount;i++) {
        bool deviceOK = false;
        for (int j=0;j < deviceCount;j++) {
          if (devices[j] == i) deviceOK = true;
        }
        if (!deviceOK) masterPeList[i] = -1;
      }
      // Device IDs deviceCount ... MAX_NUM_DEVICES are not possible, just set them to -1
      for (int i=deviceCount;i < MAX_NUM_DEVICES;i++) {
        masterPeList[i] = -1;
      }
    }

    // disable token-passing but don't submit local until remote finished
    // if shared_gpu is true, otherwise submit all work immediately
    firstPeSharingGpu = CkMyPe();
    nextPeSharingGpu = CkMyPe();

    gpuIsMine = ( firstPeSharingGpu == CkMyPe() ); 

    if ( dev >= deviceCount ) {
      char buf[256];
      sprintf(buf,"Pe %d unable to bind to CUDA device %d on %s because only %d devices are present",
  		CkMyPe(), dev, host, deviceCount);
      NAMD_die(buf);
    }

    cudaDeviceProp deviceProp;
    cudaCheck(cudaGetDeviceProperties(&deviceProp, dev));
    if ( CmiPhysicalNodeID(masterPe) < 2 )
    	CkPrintf("Pe %d physical rank %d binding to CUDA device %d on %s: '%s'  Mem: %dMB  Rev: %d.%d\n",
               CkMyPe(), myRankInPhysicalNode, dev, host,
               deviceProp.name, deviceProp.totalGlobalMem / (1024*1024),
               deviceProp.major, deviceProp.minor);

    cudaCheck(cudaSetDevice(dev));

  }  // just let CUDA pick a device for us

  {
    cudaCheck(cudaSetDeviceFlags(cudaDeviceMapHost));

    int dev;
    cudaCheck(cudaGetDevice(&dev));
    deviceID = dev;
    cudaDeviceProp deviceProp;
    cudaCheck(cudaGetDeviceProperties(&deviceProp, dev));
    if ( deviceProp.computeMode == cudaComputeModeProhibited )
      cudaDie("device in prohibited mode");
    if ( deviceProp.major < 3 )
      cudaDie("device not of compute capability 3.0 or higher");
    if ( ! deviceProp.canMapHostMemory )
      cudaDie("device cannot map host memory");
  #ifndef DISABLE_CUDA_TEXTURE_OBJECTS
    if (deviceProp.major < 3)
      cudaDie("CUDA texture objects require compute capability 3.0 or higher.\nUse DISABLE_CUDA_TEXTURE_OBJECTS to disable texture objects.");
  #endif
    extern int read_CUDA_ARCH();
    cuda_arch = read_CUDA_ARCH();
  }
}

//
// Class destructor
//
DeviceCUDA::~DeviceCUDA() {
  if (deviceProps != NULL) delete [] deviceProps;
  if (devices != NULL) delete [] devices;
	delete [] pesSharingDevice;
}

//
// Return device ID for pe. Assumes all nodes are the same
//
int DeviceCUDA::getDeviceIDforPe(int pe) {
  return deviceIDList[CkRankOf(pe) % CkMyNodeSize()];
}

//
// Returns master PE for the device ID, or -1 if device not found
//
int DeviceCUDA::getMasterPeForDeviceID(int deviceID) {
  return masterPeList[deviceID % deviceCount];
}

//
// Returns true if process "pe" shares this device
//
bool DeviceCUDA::device_shared_with_pe(int pe) {
  for ( int i=0; i<numPesSharingDevice; ++i ) {
    if ( pesSharingDevice[i] == pe ) return true;
  }
  return false;
}

//
// Returns true if there is single device per node
//
bool DeviceCUDA::one_device_per_node() {
  if ( numPesSharingDevice != CkMyNodeSize() ) return false;
  int numPesOnNodeSharingDevice = 0;
  for ( int i=0; i<numPesSharingDevice; ++i ) {
    if ( CkNodeOf(pesSharingDevice[i]) == CkMyNode() ) {
      ++numPesOnNodeSharingDevice;
    }
  }
  return ( numPesOnNodeSharingDevice == CkMyNodeSize() );
}

int DeviceCUDA::getMaxNumThreads() {
  int dev;
  cudaCheck(cudaGetDevice(&dev));
  return deviceProps[dev].maxThreadsPerBlock;
}

int DeviceCUDA::getMaxNumBlocks() {
  int dev;
  cudaCheck(cudaGetDevice(&dev));
  return deviceProps[dev].maxGridSize[0];
}

/*
BASE
2 types (remote & local)
16 pes per node
3 phases (1, 2, 3)
*/

void DeviceCUDA::register_user_events() {

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

#endif  // NAMD_CUDA

