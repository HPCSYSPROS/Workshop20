#ifndef CUDAPMESOLVER_H
#define CUDAPMESOLVER_H
#include "PmeSolver.h"
#include "CudaPmeSolver.decl.h"

#ifdef NAMD_CUDA
class CudaPmeXYZInitMsg : public CMessage_CudaPmeXYZInitMsg {
public:
	CudaPmeXYZInitMsg(PmeGrid& pmeGrid) : pmeGrid(pmeGrid) {}
	PmeGrid pmeGrid;
};

class CudaPmeXYInitMsg : public CMessage_CudaPmeXYInitMsg {
public:
	CudaPmeXYInitMsg(PmeGrid& pmeGrid, CProxy_CudaPmePencilXY& pmePencilXY, CProxy_CudaPmePencilZ& pmePencilZ,
		CProxy_PmePencilXYMap& xyMap, CProxy_PmePencilXMap& zMap) : 
		pmeGrid(pmeGrid), pmePencilXY(pmePencilXY), pmePencilZ(pmePencilZ), xyMap(xyMap), zMap(zMap) {}
	PmeGrid pmeGrid;
  CProxy_CudaPmePencilXY pmePencilXY;
  CProxy_CudaPmePencilZ pmePencilZ;
  CProxy_PmePencilXMap zMap;
  CProxy_PmePencilXYMap xyMap;
};

class CudaPmeXInitMsg : public CMessage_CudaPmeXInitMsg {
public:
	CudaPmeXInitMsg(PmeGrid& pmeGrid,
		CProxy_CudaPmePencilX& pmePencilX, CProxy_CudaPmePencilY& pmePencilY, CProxy_CudaPmePencilZ& pmePencilZ,
		CProxy_PmePencilXMap& xMap, CProxy_PmePencilXMap& yMap, CProxy_PmePencilXMap& zMap) : 
		pmeGrid(pmeGrid), pmePencilX(pmePencilX), pmePencilY(pmePencilY), pmePencilZ(pmePencilZ),
		xMap(xMap), yMap(yMap), zMap(zMap) {}
	PmeGrid pmeGrid;
  CProxy_CudaPmePencilX pmePencilX;
  CProxy_CudaPmePencilY pmePencilY;
  CProxy_CudaPmePencilZ pmePencilZ;
  CProxy_PmePencilXMap xMap;
  CProxy_PmePencilXMap yMap;
  CProxy_PmePencilXMap zMap;
};

class InitDeviceMsg : public CMessage_InitDeviceMsg {
public:
	InitDeviceMsg(CProxy_ComputePmeCUDADevice deviceProxy) : deviceProxy(deviceProxy) {}
	CProxy_ComputePmeCUDADevice deviceProxy;
};

class InitDeviceMsg2 : public CMessage_InitDeviceMsg2 {
public:
	InitDeviceMsg2(int deviceID, cudaStream_t stream, CProxy_ComputePmeCUDAMgr mgrProxy) : 
	deviceID(deviceID), stream(stream), mgrProxy(mgrProxy) {}
	int deviceID;
	cudaStream_t stream;
	CProxy_ComputePmeCUDAMgr mgrProxy;
};

class CudaPmePencilXYZ : public CBase_CudaPmePencilXYZ {
public:
	CudaPmePencilXYZ() {}
	CudaPmePencilXYZ(CkMigrateMessage *m) {}
	void initialize(CudaPmeXYZInitMsg *msg);
	void initializeDevice(InitDeviceMsg *msg);
	void energyAndVirialDone();
private:
	void backwardDone();
  CProxy_ComputePmeCUDADevice deviceProxy;
};

struct DeviceBuffer {
	DeviceBuffer(int deviceID, bool isPeerDevice, float2* data) : deviceID(deviceID), isPeerDevice(isPeerDevice), data(data) {}
	bool isPeerDevice;
	int deviceID;
	cudaEvent_t event;
	float2 *data;
};

class DeviceDataMsg : public CMessage_DeviceDataMsg {
public:
	DeviceDataMsg(int i, cudaEvent_t event, float2 *data) : i(i), event(event), data(data) {}
	int i;
	cudaEvent_t event;
	float2 *data;
};

class CudaPmePencilXY : public CBase_CudaPmePencilXY {
public:
	CudaPmePencilXY_SDAG_CODE
	CudaPmePencilXY() : numGetDeviceBuffer(0), eventCreated(false) {}
	CudaPmePencilXY(CkMigrateMessage *m) : numGetDeviceBuffer(0), eventCreated(false) {}
	~CudaPmePencilXY();
	void initialize(CudaPmeXYInitMsg *msg);
	void initializeDevice(InitDeviceMsg *msg);
private:
	void forwardDone();
	void backwardDone();
	void recvDataFromZ(PmeBlockMsg *msg);
	void start();
	void setDeviceBuffers();
	float2* getData(const int i, const bool sameDevice);
	int deviceID;
	cudaStream_t stream;
	cudaEvent_t event;
	bool eventCreated;
	int imsgZ;
	int numDeviceBuffers;
	int numGetDeviceBuffer;
	std::vector<DeviceBuffer> deviceBuffers;
  CProxy_ComputePmeCUDADevice deviceProxy;
  CProxy_CudaPmePencilZ pmePencilZ;
  CProxy_PmePencilXMap zMap;
};

class CudaPmePencilX : public CBase_CudaPmePencilX {
public:
	CudaPmePencilX_SDAG_CODE
	CudaPmePencilX() : numGetDeviceBuffer(0), eventCreated(false) {}
	CudaPmePencilX(CkMigrateMessage *m) : numGetDeviceBuffer(0), eventCreated(false) {}
	~CudaPmePencilX();
	void initialize(CudaPmeXInitMsg *msg);
	void initializeDevice(InitDeviceMsg *msg);
private:
	void forwardDone();
	void backwardDone();
	void recvDataFromY(PmeBlockMsg *msg);
	void start();
	void setDeviceBuffers();
	float2* getData(const int i, const bool sameDevice);
	int deviceID;
	cudaStream_t stream;
	cudaEvent_t event;
	bool eventCreated;
	int imsgY;
	int numDeviceBuffers;
	int numGetDeviceBuffer;
	std::vector<DeviceBuffer> deviceBuffers;
  CProxy_ComputePmeCUDADevice deviceProxy;
  CProxy_CudaPmePencilY pmePencilY;
  CProxy_PmePencilXMap yMap;
};

class CudaPmePencilY : public CBase_CudaPmePencilY {
public:
	CudaPmePencilY_SDAG_CODE
	CudaPmePencilY() : numGetDeviceBufferZ(0), numGetDeviceBufferX(0), eventCreated(false) {}
	CudaPmePencilY(CkMigrateMessage *m) : numGetDeviceBufferZ(0), numGetDeviceBufferX(0), eventCreated(false) {}
	~CudaPmePencilY();
	void initialize(CudaPmeXInitMsg *msg);
	void initializeDevice(InitDeviceMsg2 *msg);
private:
	void forwardDone();
	void backwardDone();
	void recvDataFromX(PmeBlockMsg *msg);
	void recvDataFromZ(PmeBlockMsg *msg);
	void start();
	void setDeviceBuffers();
	float2* getDataForX(const int i, const bool sameDevice);
	float2* getDataForZ(const int i, const bool sameDevice);
	int deviceID;
	cudaStream_t stream;
	cudaEvent_t event;
	bool eventCreated;
	int imsgZ, imsgX;
	int imsgZZ, imsgXX;
	int numGetDeviceBufferZ;
	int numGetDeviceBufferX;
	int numDeviceBuffersZ;
	int numDeviceBuffersX;
	std::vector<DeviceBuffer> deviceBuffersZ;
	std::vector<DeviceBuffer> deviceBuffersX;
  CProxy_CudaPmePencilX pmePencilX;
  CProxy_CudaPmePencilZ pmePencilZ;
  CProxy_PmePencilXMap xMap;
  CProxy_PmePencilXMap zMap;
};

class CudaPmePencilZ : public CBase_CudaPmePencilZ {
public:
	CudaPmePencilZ_SDAG_CODE
	CudaPmePencilZ() : numGetDeviceBufferY(0), numGetDeviceBufferXY(0), eventCreated(false) {}
	CudaPmePencilZ(CkMigrateMessage *m) : numGetDeviceBufferY(0), numGetDeviceBufferXY(0), eventCreated(false) {}
	~CudaPmePencilZ();
	void initialize(CudaPmeXInitMsg *msg);
	void initialize(CudaPmeXYInitMsg *msg);
	void initializeDevice(InitDeviceMsg2 *msg);
	void energyAndVirialDone();
private:
	void backwardDone();
	void recvDataFromY(PmeBlockMsg *msg);
	void start();
	void setDeviceBuffers();
	float2* getData(const int i, const bool sameDevice);
	int deviceID;
	cudaStream_t stream;
	cudaEvent_t event;
	bool eventCreated;
	int imsgY;
	int numDeviceBuffers;
	int numGetDeviceBufferY;
	std::vector<DeviceBuffer> deviceBuffers;
  CProxy_CudaPmePencilY pmePencilY;
  CProxy_PmePencilXMap yMap;

	bool useXYslab;
	int numGetDeviceBufferXY;
  CProxy_CudaPmePencilXY pmePencilXY;
  CProxy_PmePencilXYMap xyMap;
};

#endif // NAMD_CUDA
#endif //CUDAPMESOLVER_H