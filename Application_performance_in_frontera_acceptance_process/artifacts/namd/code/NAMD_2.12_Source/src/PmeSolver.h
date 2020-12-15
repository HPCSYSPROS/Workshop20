#ifndef PMESOLVER_H
#define PMESOLVER_H
#include <vector>
#include "ReductionMgr.h"
#include "PatchMap.h"
#include "PmeSolverUtil.h"
#include "PmeSolver.decl.h"

class PmePencilXYZMap : public CBase_PmePencilXYZMap {
public:
  PmePencilXYZMap(int pe) : pe(pe) {
  }
  //PmePencilXYZMap(CkMigrateMessage *m) {}
  int registerArray(CkArrayIndex& numElements, CkArrayID aid) {
    return 0;
  }
  virtual int procNum(int, const CkArrayIndex& idx) {
    return pe;
  }
  virtual void populateInitial(int, CkArrayOptions &, void *msg, CkArrMgr *mgr) {
    if (pe == CkMyPe()) {
      if ( ! msg ) NAMD_bug("PmePencilXYZMap::populateInitial, multiple pencils on a pe?");
      mgr->insertInitial(CkArrayIndex1D(0), msg);
      msg = NULL;
    }
    mgr->doneInserting();
    if (msg != NULL) CkFreeMsg(msg);
  }
private:
  const int pe;
};

class PmePencilXMap : public CBase_PmePencilXMap {
public:
  PmePencilXMap(int ia, int ib, int width, std::vector<int>& pes) : ia(ia), ib(ib), width(width), pes(pes) {}
  int registerArray(CkArrayIndex& numElements, CkArrayID aid) {
    return 0;
  }
  virtual int procNum(int, const CkArrayIndex& idx) {
    int ind = idx.data()[ia] + idx.data()[ib] * width;
    if (ind < 0 || ind >= pes.size())
      NAMD_bug("PmePencilXMap::procNum, index out of bounds");
    return pes[ind];
  }
  virtual void populateInitial(int, CkArrayOptions &, void *msg, CkArrMgr *mgr) {
    for (int i=0;i < pes.size();i++) {
      if (pes[i] == CkMyPe()) {
        if ( msg == NULL ) NAMD_bug("PmePencilXMap::populateInitial, multiple pencils on a pe?");
        CkArrayIndex3D ai(0,0,0);
        ai.data()[ib] = i / width;
        ai.data()[ia] = i % width;
        //fprintf(stderr, "Pe %d i %d at %d %d\n", pes[i], i, ai.data()[ia], ai.data()[ib]);
        if ( procNum(0,ai) != CkMyPe() ) NAMD_bug("PmePencilXMap::populateInitial, map is inconsistent");
        mgr->insertInitial(ai,msg);
        msg = NULL;
      }
    }
    mgr->doneInserting();
    if (msg != NULL) CkFreeMsg(msg);
  }
private:
  // Index of CkArrayIndex data()
  const int ia, ib;
  // Width of the 2D array in pes
  const int width;
  // List of Pes. Index is given by pes[i + j*width]
  const std::vector<int> pes;
};

class PmePencilXYMap : public CBase_PmePencilXYMap {
public:
  PmePencilXYMap(std::vector<int>& pes) : pes(pes) {}
  int registerArray(CkArrayIndex& numElements, CkArrayID aid) {
    return 0;
  }
  virtual int procNum(int, const CkArrayIndex& idx) {
    int ind = idx.data()[2];
    if (ind < 0 || ind >= pes.size())
      NAMD_bug("PmePencilXYMap::procNum, index out of bounds");
    return pes[ind];
  }
  virtual void populateInitial(int, CkArrayOptions &, void *msg, CkArrMgr *mgr) {
    for (int i=0;i < pes.size();i++) {
      if (pes[i] == CkMyPe()) {
        if ( msg == NULL ) NAMD_bug("PmePencilXYMap::populateInitial, multiple pencils on a pe?");
        CkArrayIndex3D ai(0,0,0);
        ai.data()[2] = i;
        if ( procNum(0,ai) != CkMyPe() ) NAMD_bug("PmePencilXYMap::populateInitial, map is inconsistent");
        mgr->insertInitial(ai, msg);
        msg = NULL;
      }
    }
    mgr->doneInserting();
    if (msg != NULL) CkFreeMsg(msg);
  }
private:
  // List of Pes.
  const std::vector<int> pes;
};

class PmeStartMsg : public CMessage_PmeStartMsg {
public:
  float *data;
  int dataSize;
  int device;
};

class PmeRunMsg : public CMessage_PmeRunMsg {
public:
  bool doEnergy, doVirial;
  int numStrayAtoms;
  Lattice lattice;
};

class PmeSelfEnergyMsg : public CMessage_PmeSelfEnergyMsg {
public:
  double energy;
};

class PmeDoneMsg : public CMessage_PmeDoneMsg {
public:
  PmeDoneMsg(int i, int j) : i(i), j(j) {}
  int i, j;
};

class PmeBlockMsg : public CMessage_PmeBlockMsg {
public:
  float2 *data;
  int dataSize;
  int x, y, z;
  bool doEnergy, doVirial;
  int numStrayAtoms;
  Lattice lattice;
};

class PmePencilXYZ : public CBase_PmePencilXYZ {
public:
  PmePencilXYZ_SDAG_CODE
  PmePencilXYZ();
  PmePencilXYZ(CkMigrateMessage *m);
  virtual ~PmePencilXYZ();
  void skip();
protected:
  PmeGrid pmeGrid;
  bool doEnergy, doVirial;
  FFTCompute* fftCompute;
  PmeKSpaceCompute* pmeKSpaceCompute;
  Lattice lattice;
  int numStrayAtoms;
  virtual void backwardDone();
  void submitReductions();
private:
  void forwardFFT();
  void backwardFFT();
  void forwardDone();
  void initFFT(PmeStartMsg *msg);

  int numSelfEnergyRecv;
  bool firstIter;
  double selfEnergy;
  SubmitReduction* reduction;

};

class PmePencilXY : public CBase_PmePencilXY {
public:
  PmePencilXY_SDAG_CODE
  PmePencilXY();
  PmePencilXY(CkMigrateMessage *m);
  virtual ~PmePencilXY();
protected:
  PmeGrid pmeGrid;
  bool doEnergy, doVirial;
  FFTCompute* fftCompute;
  PmeTranspose* pmeTranspose;
  std::vector<int> blockSizes;
  Lattice lattice;
  int numStrayAtoms;
  void initBlockSizes();
  int imsg;
private:
  void forwardFFT();
  void backwardFFT();
  void initFFT(PmeStartMsg *msg);
  virtual void forwardDone();
  virtual void backwardDone();
  virtual void recvDataFromZ(PmeBlockMsg *msg);
  virtual void start();

};

class PmePencilX : public CBase_PmePencilX {
public:
  PmePencilX_SDAG_CODE
  PmePencilX();
  PmePencilX(CkMigrateMessage *m);
  virtual ~PmePencilX();
protected:
  PmeGrid pmeGrid;
  bool doEnergy, doVirial;
  FFTCompute* fftCompute;
  PmeTranspose* pmeTranspose;
  std::vector<int> blockSizes;
  Lattice lattice;
  int numStrayAtoms;
  void initBlockSizes();
  int imsg;
private:
  void forwardFFT();
  void backwardFFT();
  void initFFT(PmeStartMsg *msg);
  virtual void forwardDone();
  virtual void backwardDone();
  virtual void recvDataFromY(PmeBlockMsg *msg);
  virtual void start();

};

class PmePencilY : public CBase_PmePencilY {
public:
  PmePencilY_SDAG_CODE
  PmePencilY();
  PmePencilY(CkMigrateMessage *m);
  virtual ~PmePencilY();
protected:
  PmeGrid pmeGrid;
  bool doEnergy, doVirial;
  FFTCompute* fftCompute;
  PmeTranspose* pmeTranspose;
  std::vector<int> blockSizes;
  Lattice lattice;
  int numStrayAtoms;
  void initBlockSizes();
  int imsg;
private:
  void forwardFFT();
  void backwardFFT();
  void initFFT(PmeStartMsg *msg);
  virtual void forwardDone();
  virtual void backwardDone();
  virtual void recvDataFromX(PmeBlockMsg *msg);
  virtual void recvDataFromZ(PmeBlockMsg *msg);
  virtual void start();

};

class PmePencilZ : public CBase_PmePencilZ {
public:
  PmePencilZ_SDAG_CODE
  PmePencilZ();
  PmePencilZ(CkMigrateMessage *m);
  virtual ~PmePencilZ();
  void skip();
protected:
  PmeGrid pmeGrid;
  bool doEnergy, doVirial;
  FFTCompute* fftCompute;
  PmeTranspose* pmeTranspose;
  PmeKSpaceCompute* pmeKSpaceCompute;
  std::vector<int> blockSizes;
  Lattice lattice;
  int numStrayAtoms;
  void initBlockSizes();
  void submitReductions();
  int imsg;
private:
  void forwardFFT();
  void backwardFFT();
  void forwardDone();
  void initFFT(PmeStartMsg *msg);
  virtual void backwardDone();
  virtual void recvDataFromY(PmeBlockMsg *msg);
  virtual void start();

  int numSelfEnergyRecv;
  bool firstIter;
  double selfEnergy;
  SubmitReduction* reduction;

};

// #define CK_TEMPLATES_ONLY
// #include "PmeSolver.def.h"
// #undef CK_TEMPLATES_ONLY

#endif // PMESOLVER_H