#ifndef CUDACOMPUTEGBISKERNEL_H
#define CUDACOMPUTEGBISKERNEL_H
#ifdef NAMD_CUDA

class CudaComputeGBISKernel {
private:

  int deviceID;

  float* intRad0;
  int intRad0Size;

  float* intRadS;
  int intRadSSize;

  float* psiSum;
  int psiSumSize;

  float* bornRad;
  int bornRadSize;

  float* dEdaSum;
  int dEdaSumSize;

  float* dHdrPrefix;
  int dHdrPrefixSize;

public:
	CudaComputeGBISKernel(int deviceID);
	~CudaComputeGBISKernel();

	void updateIntRad(const int atomStorageSize, float* intRad0H, float* intRadSH,
		cudaStream_t stream);

  void updateBornRad(const int atomStorageSize, float* bornRadH, cudaStream_t stream);

  void update_dHdrPrefix(const int atomStorageSize, float* dHdrPrefixH, cudaStream_t stream);

	void GBISphase1(CudaTileListKernel& tlKernel, const int atomStorageSize,
		const float latticeX, const float latticeY, const float latticeZ, const float a_cut, float* h_psiSum,
  	cudaStream_t stream);

  void GBISphase2(CudaTileListKernel& tlKernel, const int atomStorageSize,
    const bool doEnergy, const bool doSlow,
    const float latticeX, const float latticeY, const float latticeZ,
    const float r_cut, const float scaling, const float kappa, const float smoothDist,
    const float epsilon_p, const float epsilon_s,
    float4* d_forces,
    float* h_dEdaSum, cudaStream_t stream);

  void GBISphase3(CudaTileListKernel& tlKernel, const int atomStorageSize,
    const float latticeX, const float latticeY, const float latticeZ, const float a_cut,
    float4* d_forces,
    cudaStream_t stream);

};

#endif // NAMD_CUDA
#endif //CUDACOMPUTEGBISKERNEL_H
