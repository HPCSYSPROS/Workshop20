#include "Node.h"
#include "Priorities.h"
#include "ComputeNonbondedUtil.h"
#include "CudaPmeSolverUtil.h"
#include "ComputePmeCUDAMgr.h"
#include "ComputePmeCUDAMgr.decl.h"
#include "CudaPmeSolver.h"
#include "DeviceCUDA.h"

#ifdef NAMD_CUDA
#ifdef WIN32
#define __thread __declspec(thread)
#endif
extern __thread DeviceCUDA *deviceCUDA;
//#define DISABLE_P2P

void CudaPmePencilXYZ::initialize(CudaPmeXYZInitMsg *msg) {
  pmeGrid = msg->pmeGrid;
  delete msg;
}

//
// CUDA specific initialization
//
void CudaPmePencilXYZ::initializeDevice(InitDeviceMsg *msg) {
  // Store device proxy
  deviceProxy = msg->deviceProxy;
  delete msg;
  int deviceID = deviceProxy.ckLocalBranch()->getDeviceID();
  cudaStream_t stream = deviceProxy.ckLocalBranch()->getStream();
  CProxy_ComputePmeCUDAMgr mgrProxy = deviceProxy.ckLocalBranch()->getMgrProxy();
  // Setup fftCompute and pmeKSpaceCompute
  fftCompute = new CudaFFTCompute(deviceID, stream);
  pmeKSpaceCompute = new CudaPmeKSpaceCompute(pmeGrid, Perm_cX_Y_Z, 0, 0, 
    ComputeNonbondedUtil::ewaldcof,
    deviceCUDA->get_cuda_arch(),
    deviceID, stream);
}

void CudaPmePencilXYZ::backwardDone() {
  deviceProxy[CkMyNode()].gatherForce();
  ((CudaPmeKSpaceCompute *)pmeKSpaceCompute)->energyAndVirialSetCallback(this);

  // ((CudaPmeKSpaceCompute *)pmeKSpaceCompute)->waitEnergyAndVirial();
  // submitReductions();
  // deviceProxy[CkMyNode()].gatherForce();
}

void CudaPmePencilXYZ::energyAndVirialDone() {
  submitReductions();
  // deviceProxy[CkMyNode()].gatherForce();
}

//###########################################################################
//###########################################################################
//###########################################################################

void CudaPmePencilXY::initialize(CudaPmeXYInitMsg *msg) {
  pmeGrid = msg->pmeGrid;
  pmePencilZ = msg->pmePencilZ;
  zMap = msg->zMap;

  delete msg;

  initBlockSizes();
}

CudaPmePencilXY::~CudaPmePencilXY() {
  if (eventCreated) cudaCheck(cudaEventDestroy(event));
}

//
// CUDA specific initialization
//
void CudaPmePencilXY::initializeDevice(InitDeviceMsg *msg) {
  // Store device proxy
  deviceProxy = msg->deviceProxy;
  delete msg;
  deviceID = deviceProxy.ckLocalBranch()->getDeviceID();
  stream = deviceProxy.ckLocalBranch()->getStream();
  CProxy_ComputePmeCUDAMgr mgrProxy = deviceProxy.ckLocalBranch()->getMgrProxy();
  // Setup fftCompute and pmeKSpaceCompute
  fftCompute = new CudaFFTCompute(deviceID, stream);
  pmeTranspose = new CudaPmeTranspose(pmeGrid, Perm_cX_Y_Z, 0, thisIndex.z, deviceID, stream);  

  deviceBuffers.resize(pmeGrid.xBlocks, DeviceBuffer(-1, false, NULL));
  numDeviceBuffers = 0;

  // Create event. NOTE: Events are tied to devices, hence the cudaSetDevice() here
  cudaCheck(cudaSetDevice(deviceID));
  cudaCheck(cudaEventCreateWithFlags(&event, cudaEventDisableTiming));
  eventCreated = true;

/*
  bool useMultiGPUfft = true;
  bool allDeviceOnSameNode = true;
  for (int x=0;x < pmeGrid.xBlocks;x++) {
    int pe = zMap.ckLocalBranch()->procNum(0, CkArrayIndex3D(x,0,0));
    allDeviceOnSameNode &= (CkNodeOf(pe) == CkMyNode());
  }

  if (useMultiGPUfft && allDeviceOnSameNode && pmeGrid.xBlocks > 1) {



  } else {
*/

  for (int x=0;x < pmeGrid.xBlocks;x++) {
    int pe = zMap.ckLocalBranch()->procNum(0, CkArrayIndex3D(x,0,0));
    if (CkNodeOf(pe) == CkMyNode()) {
      // Get device ID on a device on this node
      int deviceID0 = mgrProxy.ckLocalBranch()->getDeviceIDPencilZ(x, 0);
      // Check for Peer-to-Peer access
      int canAccessPeer = 0;
      if (deviceID != deviceID0) {
        cudaCheck(cudaSetDevice(deviceID));
        cudaCheck(cudaDeviceCanAccessPeer(&canAccessPeer, deviceID, deviceID0));
#ifdef DISABLE_P2P
        canAccessPeer = 0;
#endif
        if (canAccessPeer) {
          unsigned int flags = 0;
          cudaCheck(cudaDeviceEnablePeerAccess(deviceID0, flags));
          // fprintf(stderr, "device %d can access device %d\n", deviceID, deviceID0);
        }
      }
      numDeviceBuffers++;
      deviceBuffers[x] = DeviceBuffer(deviceID0, canAccessPeer, NULL);
      pmePencilZ(x,0,0).getDeviceBuffer(thisIndex.z, (deviceID0 == deviceID) || canAccessPeer, thisProxy);
    }
  }

  // }

}

//
// CUDA specific start
//
void CudaPmePencilXY::start() {
  recvDeviceBuffers();
}

void CudaPmePencilXY::setDeviceBuffers() {
  std::vector<float2*> dataPtrs(pmeGrid.xBlocks, (float2*)0);
  for (int x=0;x < pmeGrid.xBlocks;x++) {
    if (deviceBuffers[x].data != NULL) {
      if (deviceBuffers[x].deviceID == deviceID || deviceBuffers[x].isPeerDevice) {
        // Device buffer on same device => directly transpose into destination pencil
        dataPtrs[x] = deviceBuffers[x].data;
        // Otherwise, when device buffer on different device on same node => transpose locally and then 
        // use cudaMemcpy3DPeerAsync to perform the copying
      }
    }
  }
  ((CudaPmeTranspose *)pmeTranspose)->setDataPtrsZXY(dataPtrs, (float2 *)fftCompute->getDataDst());
}

float2* CudaPmePencilXY::getData(const int i, const bool sameDevice) {
  float2* data;
#ifndef P2P_ENABLE_3D
  if (sameDevice) {
    int i0, i1, j0, j1, k0, k1;
    getBlockDim(pmeGrid, Perm_cX_Y_Z, i, 0, 0, i0, i1, j0, j1, k0, k1);
    data = (float2 *)fftCompute->getDataDst() + i0;
  } else {
    data = ((CudaPmeTranspose *)pmeTranspose)->getBuffer(i);
  }
#else
  int i0, i1, j0, j1, k0, k1;
  getBlockDim(pmeGrid, Perm_cX_Y_Z, i, 0, 0, i0, i1, j0, j1, k0, k1);
  data = (float2 *)fftCompute->getDataDst() + i0;
#endif
  return data;
}

void CudaPmePencilXY::backwardDone() {
  deviceProxy[CkMyNode()].gatherForce();
}

void CudaPmePencilXY::forwardDone() {
  // Transpose locally
  pmeTranspose->transposeXYZtoZXY((float2 *)fftCompute->getDataDst());

  // Direct Device-To-Device communication within node
  if (numDeviceBuffers > 0) {
    // Copy data
    for (int x=0;x < pmeGrid.xBlocks;x++) {
      if (deviceBuffers[x].data != NULL) {
        if (deviceBuffers[x].deviceID != deviceID && !deviceBuffers[x].isPeerDevice) {
          ((CudaPmeTranspose *)pmeTranspose)->copyDataToPeerDeviceZXY(x, deviceBuffers[x].deviceID,
            Perm_Z_cX_Y, deviceBuffers[x].data);
        }
      }
    }
    // Record event for this pencil
    cudaCheck(cudaEventRecord(event, stream));
    // Send empty message
    for (int x=0;x < pmeGrid.xBlocks;x++) {
      if (deviceBuffers[x].data != NULL) {
        PmeBlockMsg* msg = new (0, PRIORITY_SIZE) PmeBlockMsg();
        msg->dataSize = 0;
        msg->x = x;
        msg->y = thisIndex.y;
        msg->z = thisIndex.z;
        msg->doEnergy = doEnergy;
        msg->doVirial = doVirial;
        msg->lattice  = lattice;
        msg->numStrayAtoms = numStrayAtoms;
        pmePencilZ(x,0,0).recvBlock(msg);
      }
    }
  }

  // Copy-Via-Host communication
  for (int x=0;x < pmeGrid.xBlocks;x++) {
    if (deviceBuffers[x].data == NULL) {
      PmeBlockMsg* msg = new (blockSizes[x], PRIORITY_SIZE) PmeBlockMsg();
      msg->dataSize = blockSizes[x];
      msg->x = x;
      msg->y = thisIndex.y;
      msg->z = thisIndex.z;
      msg->doEnergy = doEnergy;
      msg->doVirial = doVirial;
      msg->lattice  = lattice;
      msg->numStrayAtoms = numStrayAtoms;
      ((CudaPmeTranspose *)pmeTranspose)->copyDataDeviceToHost(x, msg->data, msg->dataSize);
      ((CudaPmeTranspose *)pmeTranspose)->waitStreamSynchronize();
      pmePencilZ(x,0,0).recvBlock(msg);
    }
  }
}

void CudaPmePencilXY::recvDataFromZ(PmeBlockMsg *msg) {
  if (msg->dataSize != 0) {
    // Buffer is coming from a different node
    ((CudaPmeTranspose *)pmeTranspose)->copyDataHostToDevice(msg->x, msg->data, (float2 *)fftCompute->getDataDst());
  } else {
    // Buffer is coming from the same node
    // Wait for event that was recorded on the sending pencil
    // device ID = deviceBuffers[msg->x].deviceID
    // event     = deviceBuffers[msg->x].event
    cudaCheck(cudaStreamWaitEvent(stream, deviceBuffers[msg->x].event, 0));
#ifndef P2P_ENABLE_3D
    if (deviceBuffers[msg->x].data != NULL && deviceBuffers[msg->x].deviceID != deviceID && !deviceBuffers[msg->x].isPeerDevice) {
      // Data is in temporary device buffer, copy it into final fft-buffer
      ((CudaPmeTranspose *)pmeTranspose)->copyDataDeviceToDevice(msg->x, (float2 *)fftCompute->getDataDst());
    }
#endif
  }
  delete msg;
}

//###########################################################################
//###########################################################################
//###########################################################################

void CudaPmePencilX::initialize(CudaPmeXInitMsg *msg) {
  pmeGrid = msg->pmeGrid;
  pmePencilY = msg->pmePencilY;
  yMap = msg->yMap;

  delete msg;

  initBlockSizes();

}

CudaPmePencilX::~CudaPmePencilX() {
  if (eventCreated) cudaCheck(cudaEventDestroy(event));
}

//
// CUDA specific initialization
//
void CudaPmePencilX::initializeDevice(InitDeviceMsg *msg) {
  // Store device proxy
  deviceProxy = msg->deviceProxy;
  delete msg;
  deviceID = deviceProxy.ckLocalBranch()->getDeviceID();
  stream = deviceProxy.ckLocalBranch()->getStream();
  CProxy_ComputePmeCUDAMgr mgrProxy = deviceProxy.ckLocalBranch()->getMgrProxy();
  // Setup fftCompute and pmeKSpaceCompute
  fftCompute = new CudaFFTCompute(deviceID, stream);
  pmeTranspose = new CudaPmeTranspose(pmeGrid, Perm_cX_Y_Z, thisIndex.y, thisIndex.z, deviceID, stream);  

  // Create event. NOTE: Events are tied to devices, hence the cudaSetDevice() here
  cudaCheck(cudaSetDevice(deviceID));
  cudaCheck(cudaEventCreateWithFlags(&event, cudaEventDisableTiming));
  eventCreated = true;

  deviceBuffers.resize(pmeGrid.xBlocks, DeviceBuffer(-1, false, NULL));
  numDeviceBuffers = 0;

  for (int x=0;x < pmeGrid.xBlocks;x++) {
    int pe = yMap.ckLocalBranch()->procNum(0, CkArrayIndex3D(x,0,thisIndex.z));
    if (CkNodeOf(pe) == CkMyNode()) {
      int deviceID0 = mgrProxy.ckLocalBranch()->getDeviceIDPencilY(x, thisIndex.z);
      numDeviceBuffers++;
      deviceBuffers[x] = DeviceBuffer(deviceID0, false, NULL);
      pmePencilY(x,0,thisIndex.z).getDeviceBuffer(thisIndex.y, (deviceID0 == deviceID), thisProxy);
    }
  }

}

//
// CUDA specific start
//
void CudaPmePencilX::start() {
  recvDeviceBuffers();
}

//
// Setup direct device buffers
//
void CudaPmePencilX::setDeviceBuffers() {
  std::vector<float2*> dataPtrs(pmeGrid.xBlocks, (float2*)0);
  for (int x=0;x < pmeGrid.xBlocks;x++) {
    if (deviceBuffers[x].data != NULL) {
      if (deviceBuffers[x].deviceID == deviceID) {
        // Device buffer on same device => directly transpose into destination pencil
        dataPtrs[x] = deviceBuffers[x].data;
        // Otherwise, when device buffer on different device on same node => transpose locally and then 
        // use cudaMemcpy3DPeerAsync to perform the copying
      }
    }
  }
  ((CudaPmeTranspose *)pmeTranspose)->setDataPtrsYZX(dataPtrs, (float2 *)fftCompute->getDataDst());
}

float2* CudaPmePencilX::getData(const int i, const bool sameDevice) {
  float2* data;
#ifndef P2P_ENABLE_3D
  if (sameDevice) {
    int i0, i1, j0, j1, k0, k1;
    getBlockDim(pmeGrid, Perm_cX_Y_Z, i, 0, 0, i0, i1, j0, j1, k0, k1);
    data = (float2 *)fftCompute->getDataDst() + i0;
  } else {
    data = ((CudaPmeTranspose *)pmeTranspose)->getBuffer(i);
  }
#else
  int i0, i1, j0, j1, k0, k1;
  getBlockDim(pmeGrid, Perm_cX_Y_Z, i, 0, 0, i0, i1, j0, j1, k0, k1);
  data = (float2 *)fftCompute->getDataDst() + i0;
#endif
  return data;
}

void CudaPmePencilX::backwardDone() {
  deviceProxy[CkMyNode()].gatherForce();
}

void CudaPmePencilX::forwardDone() {
  if (pmeTranspose == NULL)
    NAMD_bug("CudaPmePencilX::forwardDone, pmeTranspose not initialized");
  if (blockSizes.size() == 0)
    NAMD_bug("CudaPmePencilX::forwardDone, blockSizes not initialized");
  // Transpose locally
  pmeTranspose->transposeXYZtoYZX((float2 *)fftCompute->getDataDst());

  // Send data to y-pencils that share the same z-coordinate. There are pmeGrid.xBlocks of them
  // Direct-Device-To-Device communication
  if (numDeviceBuffers > 0) {
    // Copy data
    for (int x=0;x < pmeGrid.xBlocks;x++) {
      if (deviceBuffers[x].data != NULL) {
        if (deviceBuffers[x].deviceID != deviceID) {
          ((CudaPmeTranspose *)pmeTranspose)->copyDataToPeerDeviceYZX(x, deviceBuffers[x].deviceID,
            Perm_Y_Z_cX, deviceBuffers[x].data);
        }
      }
    }
    // Record event for this pencil
    cudaCheck(cudaEventRecord(event, stream));
    // Send empty messages
    for (int x=0;x < pmeGrid.xBlocks;x++) {
      if (deviceBuffers[x].data != NULL) {
        PmeBlockMsg* msg = new (0, PRIORITY_SIZE) PmeBlockMsg();
        msg->dataSize = 0;
        msg->x = x;
        msg->y = thisIndex.y;
        msg->z = thisIndex.z;
        msg->doEnergy = doEnergy;
        msg->doVirial = doVirial;
        msg->lattice  = lattice;
        msg->numStrayAtoms = numStrayAtoms;
        pmePencilY(x,0,thisIndex.z).recvBlock(msg);     
      }
    }
  }

  // Copy-To-Host communication
  for (int x=0;x < pmeGrid.xBlocks;x++) {
    if (deviceBuffers[x].data == NULL) {
      PmeBlockMsg* msg = new (blockSizes[x], PRIORITY_SIZE) PmeBlockMsg();
      msg->dataSize = blockSizes[x];
      msg->x = x;
      msg->y = thisIndex.y;
      msg->z = thisIndex.z;
      msg->doEnergy = doEnergy;
      msg->doVirial = doVirial;
      msg->lattice  = lattice;
      msg->numStrayAtoms = numStrayAtoms;
      ((CudaPmeTranspose *)pmeTranspose)->copyDataDeviceToHost(x, msg->data, msg->dataSize);
      ((CudaPmeTranspose *)pmeTranspose)->waitStreamSynchronize();
      pmePencilY(x,0,thisIndex.z).recvBlock(msg);
    }
  }
}

void CudaPmePencilX::recvDataFromY(PmeBlockMsg *msg) {
  if (msg->dataSize != 0) {
    // Buffer is coming from a different node
    ((CudaPmeTranspose *)pmeTranspose)->copyDataHostToDevice(msg->x, msg->data, (float2 *)fftCompute->getDataDst());
  } else {
    // Buffer is coming from the same node
    // Wait for event that was recorded on the sending pencil
    // device ID = deviceBuffers[msg->x].deviceID
    // event     = deviceBuffers[msg->x].event
    cudaCheck(cudaStreamWaitEvent(stream, deviceBuffers[msg->x].event, 0));
#ifndef P2P_ENABLE_3D
    if (deviceBuffers[msg->x].data != NULL && deviceBuffers[msg->x].deviceID != deviceID) {
      // Data is in temporary device buffer, copy it into final fft-buffer
      ((CudaPmeTranspose *)pmeTranspose)->copyDataDeviceToDevice(msg->x, (float2 *)fftCompute->getDataDst());
    }
#endif
  }
  delete msg;
}

//###########################################################################
//###########################################################################
//###########################################################################

void CudaPmePencilY::initialize(CudaPmeXInitMsg *msg) {
  pmeGrid = msg->pmeGrid;
  pmePencilX = msg->pmePencilX;
  pmePencilZ = msg->pmePencilZ;
  xMap = msg->xMap;
  zMap = msg->zMap;

  delete msg;

  initBlockSizes();
}

CudaPmePencilY::~CudaPmePencilY() {
  if (eventCreated) cudaCheck(cudaEventDestroy(event));
}

//
// CUDA specific initialization
//
void CudaPmePencilY::initializeDevice(InitDeviceMsg2 *msg) {
  // Get device proxy
  // CProxy_ComputePmeCUDADevice deviceProxy = msg->deviceProxy;
  deviceID = msg->deviceID;
  stream = msg->stream;
  CProxy_ComputePmeCUDAMgr mgrProxy = msg->mgrProxy;
  delete msg;
  // deviceID = deviceProxy.ckLocalBranch()->getDeviceID();
  // cudaStream_t stream = deviceProxy.ckLocalBranch()->getStream();
  // CProxy_ComputePmeCUDAMgr mgrProxy = deviceProxy.ckLocalBranch()->getMgrProxy();
  // Setup fftCompute and pmeKSpaceCompute
  fftCompute = new CudaFFTCompute(deviceID, stream);
  pmeTranspose = new CudaPmeTranspose(pmeGrid, Perm_Y_Z_cX, thisIndex.z, thisIndex.x, deviceID, stream);

  // Create event. NOTE: Events are tied to devices, hence the cudaSetDevice() here
  cudaCheck(cudaSetDevice(deviceID));
  cudaCheck(cudaEventCreateWithFlags(&event, cudaEventDisableTiming));
  eventCreated = true;

  deviceBuffersZ.resize(pmeGrid.yBlocks, DeviceBuffer(-1, false, NULL));
  deviceBuffersX.resize(pmeGrid.yBlocks, DeviceBuffer(-1, false, NULL));
  numDeviceBuffersZ = 0;
  numDeviceBuffersX = 0;

  for (int y=0;y < pmeGrid.yBlocks;y++) {
    int pe;
    pe = zMap.ckLocalBranch()->procNum(0, CkArrayIndex3D(thisIndex.x, y, 0));
    if (CkNodeOf(pe) == CkMyNode()) {
      int deviceID0 = mgrProxy.ckLocalBranch()->getDeviceIDPencilZ(thisIndex.x, y);
      numDeviceBuffersZ++;
      deviceBuffersZ[y] = DeviceBuffer(deviceID0, false, NULL);
      pmePencilZ(thisIndex.x, y, 0).getDeviceBuffer(thisIndex.z, (deviceID0 == deviceID), thisProxy);
    }
    pe = xMap.ckLocalBranch()->procNum(0, CkArrayIndex3D(0, y, thisIndex.z));
    if (CkNodeOf(pe) == CkMyNode()) {
      int deviceID0 = mgrProxy.ckLocalBranch()->getDeviceIDPencilX(y, thisIndex.z);
      numDeviceBuffersX++;
      deviceBuffersX[y] = DeviceBuffer(deviceID0, false, NULL);
      pmePencilX(0, y, thisIndex.z).getDeviceBuffer(thisIndex.x, (deviceID0 == deviceID), thisProxy);
    }
  }

}

//
// CUDA specific start
//
void CudaPmePencilY::start() {
  recvDeviceBuffers();
}

//
// Setup direct device buffers
//
void CudaPmePencilY::setDeviceBuffers() {
  std::vector<float2*> dataPtrsYZX(pmeGrid.yBlocks, (float2*)0);
  std::vector<float2*> dataPtrsZXY(pmeGrid.yBlocks, (float2*)0);
  for (int y=0;y < pmeGrid.yBlocks;y++) {
    if (deviceBuffersZ[y].data != NULL) {
      if (deviceBuffersZ[y].deviceID == deviceID) {
        dataPtrsYZX[y] = deviceBuffersZ[y].data;
      }
    }
    if (deviceBuffersX[y].data != NULL) {
      if (deviceBuffersX[y].deviceID == deviceID) {
        dataPtrsZXY[y] = deviceBuffersX[y].data;
      }
    }
  }
  ((CudaPmeTranspose *)pmeTranspose)->setDataPtrsYZX(dataPtrsYZX, (float2 *)fftCompute->getDataDst());
  ((CudaPmeTranspose *)pmeTranspose)->setDataPtrsZXY(dataPtrsZXY, (float2 *)fftCompute->getDataSrc());
}

float2* CudaPmePencilY::getDataForX(const int i, const bool sameDevice) {
  float2* data;
#ifndef P2P_ENABLE_3D
  if (sameDevice) {
    int i0, i1, j0, j1, k0, k1;
    getBlockDim(pmeGrid, Perm_Y_Z_cX, i, 0, 0, i0, i1, j0, j1, k0, k1);
    data = (float2 *)fftCompute->getDataSrc() + i0;
  } else {
    data = ((CudaPmeTranspose *)pmeTranspose)->getBuffer(i);
  }
#else
  int i0, i1, j0, j1, k0, k1;
  getBlockDim(pmeGrid, Perm_Y_Z_cX, i, 0, 0, i0, i1, j0, j1, k0, k1);
  data = (float2 *)fftCompute->getDataSrc() + i0;
#endif
  return data;
}

float2* CudaPmePencilY::getDataForZ(const int i, const bool sameDevice) {
  float2* data;
#ifndef P2P_ENABLE_3D
  if (sameDevice) {
    int i0, i1, j0, j1, k0, k1;
    getBlockDim(pmeGrid, Perm_Y_Z_cX, i, 0, 0, i0, i1, j0, j1, k0, k1);
    data = (float2 *)fftCompute->getDataDst() + i0;
  } else {
    data = ((CudaPmeTranspose *)pmeTranspose)->getBuffer(i);
  }
#else
  int i0, i1, j0, j1, k0, k1;
  getBlockDim(pmeGrid, Perm_Y_Z_cX, i, 0, 0, i0, i1, j0, j1, k0, k1);
  data = (float2 *)fftCompute->getDataDst() + i0;
#endif
  return data;
}

void CudaPmePencilY::backwardDone() {
  // Transpose locally
  pmeTranspose->transposeXYZtoZXY((float2 *)fftCompute->getDataSrc());

  // Send data to x-pencils that share the same x-coordinate. There are pmeGrid.yBlocks of them
  // Direct-Device-To-Device communication
  if (numDeviceBuffersX > 0) {
    for (int y=0;y < pmeGrid.yBlocks;y++) {
      if (deviceBuffersX[y].data != NULL) {
        if (deviceBuffersX[y].deviceID != deviceID) {
          ((CudaPmeTranspose *)pmeTranspose)->copyDataToPeerDeviceZXY(y, deviceBuffersX[y].deviceID,
            Perm_cX_Y_Z, deviceBuffersX[y].data);
        }
      }
    }
    // Record event for this pencil
    cudaCheck(cudaEventRecord(event, stream));
    // Send empty message
    for (int y=0;y < pmeGrid.yBlocks;y++) {
      if (deviceBuffersX[y].data != NULL) {
        PmeBlockMsg* msg = new (0, PRIORITY_SIZE) PmeBlockMsg();
        msg->dataSize = 0;
        msg->x = thisIndex.x;
        msg->y = y;
        msg->z = thisIndex.z;
        pmePencilX(0,y,thisIndex.z).recvBlock(msg);
      }
    }
  }

  // Copy via host
  for (int y=0;y < pmeGrid.yBlocks;y++) {
    if (deviceBuffersX[y].data == NULL) {
      PmeBlockMsg* msg = new (blockSizes[y], PRIORITY_SIZE) PmeBlockMsg();
      msg->dataSize = blockSizes[y];
      msg->x = thisIndex.x;
      msg->y = y;
      msg->z = thisIndex.z;
      ((CudaPmeTranspose *)pmeTranspose)->copyDataDeviceToHost(y, msg->data, msg->dataSize);
      ((CudaPmeTranspose *)pmeTranspose)->waitStreamSynchronize();
      pmePencilX(0,y,thisIndex.z).recvBlock(msg);
    }
  }
}

void CudaPmePencilY::forwardDone() {
  if (pmeTranspose == NULL)
    NAMD_bug("CudaPmePencilY::forwardDone, pmeTranspose not initialized");
  if (blockSizes.size() == 0)
    NAMD_bug("CudaPmePencilY::forwardDone, blockSizes not initialized");

  // Transpose locally
  pmeTranspose->transposeXYZtoYZX((float2 *)fftCompute->getDataDst());

  // Send data to z-pencils that share the same x-coordinate. There are pmeGrid.yBlocks of them
  // Direct-Device-To-Device communication
  if (numDeviceBuffersZ > 0) {
    for (int y=0;y < pmeGrid.yBlocks;y++) {
      if (deviceBuffersZ[y].data != NULL) {
        if (deviceBuffersZ[y].deviceID != deviceID) {
          ((CudaPmeTranspose *)pmeTranspose)->copyDataToPeerDeviceYZX(y, deviceBuffersZ[y].deviceID,
            Perm_Z_cX_Y, deviceBuffersZ[y].data);
        }
      }
    }
    // Record event for this pencil
    cudaCheck(cudaEventRecord(event, stream));
    // Send empty message
    for (int y=0;y < pmeGrid.yBlocks;y++) {
      if (deviceBuffersZ[y].data != NULL) {
        PmeBlockMsg* msg = new (0, PRIORITY_SIZE) PmeBlockMsg();
        msg->dataSize = 0;
        msg->x = thisIndex.x;
        msg->y = y;
        msg->z = thisIndex.z;
        msg->doEnergy = doEnergy;
        msg->doVirial = doVirial;
        msg->lattice  = lattice;
        msg->numStrayAtoms = numStrayAtoms;
        pmePencilZ(thisIndex.x,y,0).recvBlock(msg);
      }
    }
  }

  // Copy-To-Host communication
  for (int y=0;y < pmeGrid.yBlocks;y++) {
    if (deviceBuffersZ[y].data == NULL) {
      PmeBlockMsg* msg = new (blockSizes[y], PRIORITY_SIZE) PmeBlockMsg();
      msg->dataSize = blockSizes[y];
      msg->x = thisIndex.x;
      msg->y = y;
      msg->z = thisIndex.z;
      msg->doEnergy = doEnergy;
      msg->doVirial = doVirial;
      msg->lattice  = lattice;
      msg->numStrayAtoms = numStrayAtoms;
      ((CudaPmeTranspose *)pmeTranspose)->copyDataDeviceToHost(y, msg->data, msg->dataSize);
      ((CudaPmeTranspose *)pmeTranspose)->waitStreamSynchronize();
      pmePencilZ(thisIndex.x,y,0).recvBlock(msg);
    }
  }
}

void CudaPmePencilY::recvDataFromX(PmeBlockMsg *msg) {
  if (msg->dataSize != 0) {
    // Buffer is coming from a different node
    ((CudaPmeTranspose *)pmeTranspose)->copyDataHostToDevice(msg->y, msg->data, (float2 *)fftCompute->getDataSrc());
  } else {
    // Buffer is coming from the same node
    // Wait for event that was recorded on the sending pencil
    // device ID = deviceBuffersX[msg->y].deviceID
    // event     = deviceBuffersX[msg->y].event
    cudaCheck(cudaStreamWaitEvent(stream, deviceBuffersX[msg->y].event, 0));
#ifndef P2P_ENABLE_3D
    if (deviceBuffersX[msg->y].data != NULL && deviceBuffersX[msg->y].deviceID != deviceID) {
      // Data is in temporary device buffer, copy it into final fft-buffer
      ((CudaPmeTranspose *)pmeTranspose)->copyDataDeviceToDevice(msg->y, (float2 *)fftCompute->getDataSrc());
    }
#endif
  }
  delete msg;
}

void CudaPmePencilY::recvDataFromZ(PmeBlockMsg *msg) {
  if (msg->dataSize != 0) {
    // Buffer is coming from a different node
    ((CudaPmeTranspose *)pmeTranspose)->copyDataHostToDevice(msg->y, msg->data, (float2 *)fftCompute->getDataDst());
  } else {
    // Buffer is coming from the same node
    // Wait for event that was recorded on the sending pencil
    // device ID = deviceBuffersZ[msg->y].deviceID
    // event     = deviceBuffersZ[msg->y].event
    cudaCheck(cudaStreamWaitEvent(stream, deviceBuffersZ[msg->y].event, 0));
#ifndef P2P_ENABLE_3D
    if (deviceBuffersZ[msg->y].data != NULL && deviceBuffersZ[msg->y].deviceID != deviceID) {
      // Data is in temporary device buffer, copy it into final fft-buffer
      ((CudaPmeTranspose *)pmeTranspose)->copyDataDeviceToDevice(msg->y, (float2 *)fftCompute->getDataDst());
    }
#endif
  }
  delete msg;
}

//###########################################################################
//###########################################################################
//###########################################################################

void CudaPmePencilZ::initialize(CudaPmeXInitMsg *msg) {
  useXYslab = false;
  pmeGrid = msg->pmeGrid;
  pmePencilY = msg->pmePencilY;
  yMap = msg->yMap;

  delete msg;

  initBlockSizes();
}

void CudaPmePencilZ::initialize(CudaPmeXYInitMsg *msg) {
  useXYslab = true;
  pmeGrid = msg->pmeGrid;
  pmePencilXY = msg->pmePencilXY;
  xyMap = msg->xyMap;

  delete msg;

  initBlockSizes();
}

CudaPmePencilZ::~CudaPmePencilZ() {
  if (eventCreated) cudaCheck(cudaEventDestroy(event));
}

//
// CUDA specific initialization
//
void CudaPmePencilZ::initializeDevice(InitDeviceMsg2 *msg) {
  // Get device proxy
  // CProxy_ComputePmeCUDADevice deviceProxy = msg->deviceProxy;
  deviceID = msg->deviceID;
  stream = msg->stream;
  CProxy_ComputePmeCUDAMgr mgrProxy = msg->mgrProxy;
  delete msg;
  // deviceID = deviceProxy.ckLocalBranch()->getDeviceID();
  // cudaStream_t stream = deviceProxy.ckLocalBranch()->getStream();
  // CProxy_ComputePmeCUDAMgr mgrProxy = deviceProxy.ckLocalBranch()->getMgrProxy();
  // Setup fftCompute and pmeKSpaceCompute
  fftCompute = new CudaFFTCompute(deviceID, stream);
  pmeTranspose = new CudaPmeTranspose(pmeGrid, Perm_Z_cX_Y, thisIndex.x, thisIndex.y, deviceID, stream);
  pmeKSpaceCompute = new CudaPmeKSpaceCompute(pmeGrid, Perm_Z_cX_Y, thisIndex.x, thisIndex.y,
    ComputeNonbondedUtil::ewaldcof,
    deviceCUDA->get_cuda_arch(),
    deviceID, stream);

  // Create event. NOTE: Events are tied to devices, hence the cudaSetDevice() here
  cudaCheck(cudaSetDevice(deviceID));
  cudaCheck(cudaEventCreateWithFlags(&event, cudaEventDisableTiming));
  eventCreated = true;

  deviceBuffers.resize(pmeGrid.zBlocks, DeviceBuffer(-1, false, NULL));
  numDeviceBuffers = 0;

  if (useXYslab) {
    for (int z=0;z < pmeGrid.zBlocks;z++) {
      int pe = xyMap.ckLocalBranch()->procNum(0, CkArrayIndex3D(0,0,z));
      if (CkNodeOf(pe) == CkMyNode()) {
        int deviceID0 = mgrProxy.ckLocalBranch()->getDeviceIDPencilX(0, z);
        // Check for Peer-to-Peer access
        int canAccessPeer = 0;
        if (deviceID != deviceID0) {
          cudaCheck(cudaSetDevice(deviceID));
          cudaCheck(cudaDeviceCanAccessPeer(&canAccessPeer, deviceID, deviceID0));
        }
#ifdef DISABLE_P2P
        canAccessPeer = 0;
#endif
        numDeviceBuffers++;
        deviceBuffers[z] = DeviceBuffer(deviceID0, canAccessPeer, NULL);
        pmePencilXY(0,0,z).getDeviceBuffer(thisIndex.x, (deviceID0 == deviceID) || canAccessPeer, thisProxy);
      }
    }
  } else {
    for (int z=0;z < pmeGrid.zBlocks;z++) {
      int pe = yMap.ckLocalBranch()->procNum(0, CkArrayIndex3D(thisIndex.x,0,z));
      if (CkNodeOf(pe) == CkMyNode()) {
        int deviceID0 = mgrProxy.ckLocalBranch()->getDeviceIDPencilY(thisIndex.x, z);
        numDeviceBuffers++;
        deviceBuffers[z] = DeviceBuffer(deviceID0, false, NULL);
        pmePencilY(thisIndex.x,0,z).getDeviceBuffer(thisIndex.y, (deviceID0 == deviceID), thisProxy);
      }
    }
  }

}

//
// CUDA specific start
//
void CudaPmePencilZ::start() {
  recvDeviceBuffers();
}

void CudaPmePencilZ::setDeviceBuffers() {
  std::vector<float2*> dataPtrs(pmeGrid.zBlocks, (float2*)0);
  for (int z=0;z < pmeGrid.zBlocks;z++) {
    if (deviceBuffers[z].data != NULL) {
      if (deviceBuffers[z].deviceID == deviceID || deviceBuffers[z].isPeerDevice) {
        dataPtrs[z] = deviceBuffers[z].data;
      }
    }
  }
  if (useXYslab) {
    ((CudaPmeTranspose *)pmeTranspose)->setDataPtrsYZX(dataPtrs, (float2 *)fftCompute->getDataSrc());
  } else {
    ((CudaPmeTranspose *)pmeTranspose)->setDataPtrsZXY(dataPtrs, (float2 *)fftCompute->getDataSrc());
  }
}

float2* CudaPmePencilZ::getData(const int i, const bool sameDevice) {
  float2* data;
#ifndef P2P_ENABLE_3D
  if (sameDevice) {
    int i0, i1, j0, j1, k0, k1;
    getBlockDim(pmeGrid, Perm_Z_cX_Y, i, 0, 0, i0, i1, j0, j1, k0, k1);
    data = (float2 *)fftCompute->getDataSrc() + i0;
  } else {
    data = ((CudaPmeTranspose *)pmeTranspose)->getBuffer(i);
  }
#else
  int i0, i1, j0, j1, k0, k1;
  getBlockDim(pmeGrid, Perm_Z_cX_Y, i, 0, 0, i0, i1, j0, j1, k0, k1);
  data = (float2 *)fftCompute->getDataSrc() + i0;
#endif
  return data;
}

void CudaPmePencilZ::backwardDone() {
  // Transpose locally
  if (useXYslab) {
    pmeTranspose->transposeXYZtoYZX((float2 *)fftCompute->getDataSrc());
  } else {
    pmeTranspose->transposeXYZtoZXY((float2 *)fftCompute->getDataSrc());   
  }

  if (useXYslab) {
    // Direct-Device-To-Device communication
    if (numDeviceBuffers > 0) {
      for (int z=0;z < pmeGrid.zBlocks;z++) {
        if (deviceBuffers[z].data != NULL) {
          if (deviceBuffers[z].deviceID != deviceID && !deviceBuffers[z].isPeerDevice) {
            ((CudaPmeTranspose *)pmeTranspose)->copyDataToPeerDeviceYZX(z, deviceBuffers[z].deviceID,
              Perm_cX_Y_Z, deviceBuffers[z].data);
          }
        }
      }
      // Record event for this pencil
      cudaCheck(cudaEventRecord(event, stream));
      // Send empty message
      for (int z=0;z < pmeGrid.zBlocks;z++) {
        if (deviceBuffers[z].data != NULL) {
          PmeBlockMsg* msg = new (0, PRIORITY_SIZE) PmeBlockMsg();
          msg->dataSize = 0;
          msg->x = thisIndex.x;
          msg->y = thisIndex.y;
          msg->z = z;
          pmePencilXY(0,0,z).recvBlock(msg);
        }
      }
    }

    // Copy-To-Host communication
    for (int z=0;z < pmeGrid.zBlocks;z++) {
      if (deviceBuffers[z].data == NULL) {
        PmeBlockMsg* msg = new (blockSizes[z], PRIORITY_SIZE) PmeBlockMsg();
        msg->dataSize = blockSizes[z];
        msg->x = thisIndex.x;
        msg->y = thisIndex.y;
        msg->z = z;
        ((CudaPmeTranspose *)pmeTranspose)->copyDataDeviceToHost(z, msg->data, msg->dataSize);
        ((CudaPmeTranspose *)pmeTranspose)->waitStreamSynchronize();
        pmePencilXY(0,0,z).recvBlock(msg);
      }
    }
  } else {
    // Send data to y-pencils that share the same x-coordinate. There are pmeGrid.zBlocks of them
    // Direct-Device-To-Device communication
    if (numDeviceBuffers > 0) {
      for (int z=0;z < pmeGrid.zBlocks;z++) {
        if (deviceBuffers[z].data != NULL) {
          if (deviceBuffers[z].deviceID != deviceID) {
            ((CudaPmeTranspose *)pmeTranspose)->copyDataToPeerDeviceZXY(z, deviceBuffers[z].deviceID,
              Perm_Y_Z_cX, deviceBuffers[z].data);
          }
        }
      }
      // Record event for this pencil
      cudaCheck(cudaEventRecord(event, stream));
      // Send empty message
      for (int z=0;z < pmeGrid.zBlocks;z++) {
        if (deviceBuffers[z].data != NULL) {
          PmeBlockMsg* msg = new (0, PRIORITY_SIZE) PmeBlockMsg();
          msg->dataSize = 0;
          msg->x = thisIndex.x;
          msg->y = thisIndex.y;
          msg->z = z;
          pmePencilY(thisIndex.x,0,z).recvBlock(msg);
        }
      }
    }

    // Copy-To-Host communication
    for (int z=0;z < pmeGrid.zBlocks;z++) {
      if (deviceBuffers[z].data == NULL) {
        PmeBlockMsg* msg = new (blockSizes[z], PRIORITY_SIZE) PmeBlockMsg();
        msg->dataSize = blockSizes[z];
        msg->x = thisIndex.x;
        msg->y = thisIndex.y;
        msg->z = z;
        ((CudaPmeTranspose *)pmeTranspose)->copyDataDeviceToHost(z, msg->data, msg->dataSize);
        ((CudaPmeTranspose *)pmeTranspose)->waitStreamSynchronize();
        pmePencilY(thisIndex.x,0,z).recvBlock(msg);
      }
    }
  }

  // Submit reductions
  ((CudaPmeKSpaceCompute *)pmeKSpaceCompute)->energyAndVirialSetCallback(this);
  // ((CudaPmeKSpaceCompute *)pmeKSpaceCompute)->waitEnergyAndVirial();
  // submitReductions();
}

void CudaPmePencilZ::energyAndVirialDone() {
  submitReductions();
}

void CudaPmePencilZ::recvDataFromY(PmeBlockMsg *msg) {
  // NOTE: No need to synchronize stream here since memory copies are in the stream
  if (msg->dataSize != 0) {
    // Buffer is coming from a different node
    ((CudaPmeTranspose *)pmeTranspose)->copyDataHostToDevice(msg->z, msg->data, (float2 *)fftCompute->getDataSrc());
  } else {
    // Buffer is coming from the same node
    // Wait for event that was recorded on the sending pencil
    // device ID = deviceBuffers[msg->z].deviceID
    // event     = deviceBuffers[msg->z].event
    cudaCheck(cudaStreamWaitEvent(stream, deviceBuffers[msg->z].event, 0));
#ifndef P2P_ENABLE_3D
    if (deviceBuffers[msg->z].data != NULL && deviceBuffers[msg->z].deviceID != deviceID && !deviceBuffers[msg->z].isPeerDevice) {
      // Data is in temporary device buffer, copy it into final fft-buffer
      ((CudaPmeTranspose *)pmeTranspose)->copyDataDeviceToDevice(msg->z, (float2 *)fftCompute->getDataSrc());
    }
#endif
  }
  delete msg;
}
#endif // NAMD_CUDA

#include "CudaPmeSolver.def.h"
