#include <stdio.h>
#include "Priorities.h"
struct float2 {float x,y;};
#include "PmeSolver.h"

//
// Data flow for PmePencilXYZ
//
// dataSrc [xyz]   dataDst
//
// dataDst [solve] dataDst
//
// dataDst [xyz]   dataSrc
//
// dataSrc [force]
//

PmePencilXYZ::PmePencilXYZ() {
  __sdag_init();
  setMigratable(false);
  fftCompute = NULL;
  pmeKSpaceCompute = NULL;
  reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);
  numSelfEnergyRecv = 0;
  selfEnergy = 0.0;
  doEnergy = false;
  doVirial = false;
}

PmePencilXYZ::PmePencilXYZ(CkMigrateMessage *m) {
  NAMD_bug("PmePencilXYZ cannot be migrated");
  //__sdag_init();
  // setMigratable(false);
  // fftCompute = NULL;
  // pmeKSpaceCompute = NULL;
  // reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);
}

PmePencilXYZ::~PmePencilXYZ() {
  if (fftCompute != NULL) delete fftCompute;
  if (pmeKSpaceCompute != NULL) delete pmeKSpaceCompute;
  delete reduction;
}

void PmePencilXYZ::initFFT(PmeStartMsg *msg) {
  if (fftCompute == NULL)
    NAMD_bug("PmePencilXYZ::initFFT, fftCompute not initialized");
  fftCompute->init(msg->data, msg->dataSize, NULL, 0, Perm_X_Y_Z, pmeGrid, 3, 0, 0, 0);
}

void PmePencilXYZ::forwardFFT() {
  if (fftCompute == NULL)
    NAMD_bug("PmePencilXYZ::forwardFFT, fftCompute not initialized");
  fftCompute->forward();
}

void PmePencilXYZ::backwardFFT() {
  if (fftCompute == NULL)
    NAMD_bug("PmePencilXYZ::backwardFFT, fftCompute not initialized");
  fftCompute->backward();
}

void PmePencilXYZ::forwardDone() {
  if (pmeKSpaceCompute == NULL)
    NAMD_bug("PmePencilXYZ::forwardDone, pmeKSpaceCompute not initialized");
  pmeKSpaceCompute->solve(lattice, doEnergy, doVirial, fftCompute->getDataDst());
}

void PmePencilXYZ::backwardDone() {
  NAMD_bug("PmePencilXYZ::backwardDone(), base class method called");
}

void PmePencilXYZ::submitReductions() {
  if (pmeKSpaceCompute == NULL)
    NAMD_bug("PmePencilXYZ::submitReductions, pmeKSpaceCompute not initialized");
  double virial[9];
  double energy = pmeKSpaceCompute->getEnergy();
  pmeKSpaceCompute->getVirial(virial);
  reduction->item(REDUCTION_ELECT_ENERGY_SLOW) += energy + selfEnergy;
  reduction->item(REDUCTION_VIRIAL_SLOW_XX) += virial[0];
  reduction->item(REDUCTION_VIRIAL_SLOW_XY) += virial[1];
  reduction->item(REDUCTION_VIRIAL_SLOW_XZ) += virial[2];
  reduction->item(REDUCTION_VIRIAL_SLOW_YX) += virial[3];
  reduction->item(REDUCTION_VIRIAL_SLOW_YY) += virial[4];
  reduction->item(REDUCTION_VIRIAL_SLOW_YZ) += virial[5];
  reduction->item(REDUCTION_VIRIAL_SLOW_ZX) += virial[6];
  reduction->item(REDUCTION_VIRIAL_SLOW_ZY) += virial[7];
  reduction->item(REDUCTION_VIRIAL_SLOW_ZZ) += virial[8];
  reduction->item(REDUCTION_STRAY_CHARGE_ERRORS) += numStrayAtoms;
  reduction->submit();
}

void PmePencilXYZ::skip() {
  reduction->submit();  
}

//###########################################################################
//###########################################################################
//###########################################################################

//
// Data flow for PmePencilXY & PmePencilZ
//
// dataSrc(XY) [xy]     dataDst(XY)
//
// dataDst(XY) [transp] dataSrc(Z)
//---------------------------------
//
// dataSrc(Z)  [z]      dataDst(Z)
//
// dataDst(Z)  [solve]  dataDst(Z)
//
// dataDst(Z)  [z]      dataSrc(Z)
//
// dataSrc(Z)  [transp] dataDst(XY)
//---------------------------------
//
// dataDst(XY) [xy]     dataSrc(XY)
//
// dataSrc(XY) [force]
//

PmePencilXY::PmePencilXY() {
  __sdag_init();
  setMigratable(false);
  fftCompute = NULL;
  pmeTranspose = NULL;
}

PmePencilXY::PmePencilXY(CkMigrateMessage *m) {
  NAMD_bug("PmePencilXY cannot be migrated");
//__sdag_init();
  // setMigratable(false);
  // fftCompute = NULL;
  // pmeTranspose = NULL;
}

PmePencilXY::~PmePencilXY() {
  if (fftCompute != NULL) delete fftCompute;
  if (pmeTranspose != NULL) delete pmeTranspose;
}

void PmePencilXY::initFFT(PmeStartMsg *msg) {
  if (fftCompute == NULL)
    NAMD_bug("PmePencilXY::initFFT, fftCompute not initialized");
  fftCompute->init(msg->data, msg->dataSize,  NULL, 0, Perm_X_Y_Z, pmeGrid, 2, 0, thisIndex.z, 0);
}

void PmePencilXY::forwardFFT() {
  if (fftCompute == NULL)
    NAMD_bug("PmePencilXY::forwardFFT, fftCompute not initialized");
  fftCompute->forward();
}

void PmePencilXY::backwardFFT() {
  if (fftCompute == NULL)
    NAMD_bug("PmePencilXY::backwardFFT, fftCompute not initialized");
  fftCompute->backward();
}

void PmePencilXY::initBlockSizes() {
  blockSizes.resize(pmeGrid.xBlocks);
  for (int x=0;x < pmeGrid.xBlocks;x++) {
    int i0, i1, j0, j1, k0, k1;
    getBlockDim(pmeGrid, Perm_cX_Y_Z, x, 0, thisIndex.z,
      i0, i1, j0, j1, k0, k1);
    int size = (i1-i0+1)*(j1-j0+1)*(k1-k0+1);
    blockSizes[x] = size;
  }
}

void PmePencilXY::forwardDone() {
  NAMD_bug("PmePencilXY::forwardDone(), base class method called");
}

void PmePencilXY::backwardDone() {
  NAMD_bug("PmePencilXY::backwardDone(), base class method called");
}

void PmePencilXY::recvDataFromZ(PmeBlockMsg *msg) {
  NAMD_bug("PmePencilXY::recvDataFromZ(), base class method called");
}

void PmePencilXY::start() {
  NAMD_bug("PmePencilXY::start(), base class method called");
}

//###########################################################################
//###########################################################################
//###########################################################################

//
// Data flow for PmePencilX & PmePencilX & PmePencilZ
//
// dataSrc(X) [x]     dataDst(X)
//
// dataDst(X) [transp] dataSrc(Y)
//---------------------------------
//
// dataSrc(Y) [y]      dataDst(Y)
//
// dataDst(Y) [transp] dataSrc(Z)
//---------------------------------
//
// dataSrc(Z) [z]      dataDst(Z)
//
// dataDst(Z) [solve]  dataDst(Z)
//
// dataDst(Z) [z]      dataSrc(Z)
//
// dataSrc(Z) [transp] dataDst(Y)
//---------------------------------
//
// dataDst(Y) [y]      dataSrc(Y)
//
// dataSrc(Y) [transp] dataDst(X)
//---------------------------------
//
// dataDst(X) [x]      dataSrc(X)
//
// dataSrc(X) [force]
//

PmePencilX::PmePencilX() {
  __sdag_init();
  setMigratable(false);
  fftCompute = NULL;
  pmeTranspose = NULL;
  numStrayAtoms = 0;
}

PmePencilX::PmePencilX(CkMigrateMessage *m) {
  NAMD_bug("PmePencilX cannot be migrated");
//__sdag_init();
  // setMigratable(false);
  // fftCompute = NULL;
  // pmeTranspose = NULL;
}

PmePencilX::~PmePencilX() {
  if (fftCompute != NULL) delete fftCompute;
  if (pmeTranspose != NULL) delete pmeTranspose;
}

void PmePencilX::initFFT(PmeStartMsg *msg) {
  if (fftCompute == NULL)
    NAMD_bug("PmePencilX::initFFT, fftCompute not initialized");
  fftCompute->init(msg->data, msg->dataSize,  NULL, 0, Perm_X_Y_Z, pmeGrid, 1, thisIndex.y, thisIndex.z, 0);
}

void PmePencilX::forwardFFT() {
  if (fftCompute == NULL)
    NAMD_bug("PmePencilX::forwardFFT, fftCompute not initialized");
  fftCompute->forward();
}

void PmePencilX::backwardFFT() {
  if (fftCompute == NULL)
    NAMD_bug("PmePencilX::backwardFFT, fftCompute not initialized");
  fftCompute->backward();
}

void PmePencilX::initBlockSizes() {
  blockSizes.resize(pmeGrid.xBlocks);
  for (int x=0;x < pmeGrid.xBlocks;x++) {
    int i0, i1, j0, j1, k0, k1;
    getBlockDim(pmeGrid, Perm_cX_Y_Z, x, thisIndex.y, thisIndex.z,
      i0, i1, j0, j1, k0, k1);
    int size = (i1-i0+1)*(j1-j0+1)*(k1-k0+1);
    blockSizes[x] = size;
  }
}

void PmePencilX::forwardDone() {
  NAMD_bug("PmePencilX::forwardDone(), base class method called");
}

void PmePencilX::backwardDone() {
  NAMD_bug("PmePencilX::backwardDone(), base class method called");
}

void PmePencilX::recvDataFromY(PmeBlockMsg *msg) {
  NAMD_bug("PmePencilX::recvDataFromY(), base class method called");
}

void PmePencilX::start() {
  NAMD_bug("PmePencilX::start(), base class method called");
}

//###########################################################################
//###########################################################################
//###########################################################################

PmePencilY::PmePencilY() {
  __sdag_init();
  setMigratable(false);
  fftCompute = NULL;
  pmeTranspose = NULL;
  numStrayAtoms = 0;
}

PmePencilY::PmePencilY(CkMigrateMessage *m) {
  NAMD_bug("PmePencilY cannot be migrated");
  // __sdag_init();
  // setMigratable(false);
  // fftCompute = NULL;
  // pmeTranspose = NULL;
}

PmePencilY::~PmePencilY() {
  if (fftCompute != NULL) delete fftCompute;
  if (pmeTranspose != NULL) delete pmeTranspose;
}

void PmePencilY::initFFT(PmeStartMsg *msg) {
  if (fftCompute == NULL)
    NAMD_bug("PmePencilY::initFFT, fftCompute not initialized");
  fftCompute->init(msg->data, msg->dataSize,  NULL, 0, Perm_Y_Z_cX, pmeGrid, 1, thisIndex.z, thisIndex.x, 0);
}

void PmePencilY::forwardFFT() {
  if (fftCompute == NULL)
    NAMD_bug("PmePencilY::forwardFFT, fftCompute not initialized");
  fftCompute->forward();
}

void PmePencilY::backwardFFT() {
  if (fftCompute == NULL)
    NAMD_bug("PmePencilY::backwardFFT, fftCompute not initialized");
  fftCompute->backward();
}

void PmePencilY::initBlockSizes() {
  blockSizes.resize(pmeGrid.yBlocks);
  for (int y=0;y < pmeGrid.yBlocks;y++) {
    int i0, i1, j0, j1, k0, k1;
    getBlockDim(pmeGrid, Perm_Y_Z_cX, y, thisIndex.z, thisIndex.x,
      i0, i1, j0, j1, k0, k1);
    int size = (i1-i0+1)*(j1-j0+1)*(k1-k0+1);
    blockSizes[y] = size;
  }
}

void PmePencilY::forwardDone() {
  NAMD_bug("PmePencilY::forwardDone(), base class method called");
}

void PmePencilY::backwardDone() {
  NAMD_bug("PmePencilY::backwardDone(), base class method called");
}

void PmePencilY::recvDataFromX(PmeBlockMsg *msg) {
  NAMD_bug("PmePencilY::recvDataFromX(), base class method called");
}

void PmePencilY::recvDataFromZ(PmeBlockMsg *msg) {
  NAMD_bug("PmePencilY::recvDataFromZ(), base class method called");
}

void PmePencilY::start() {
  NAMD_bug("PmePencilY::start(), base class method called");
}

//###########################################################################
//###########################################################################
//###########################################################################

PmePencilZ::PmePencilZ() {
  __sdag_init();
  setMigratable(false);
  fftCompute = NULL;
  pmeTranspose = NULL;
  pmeKSpaceCompute = NULL;
  reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);
  numSelfEnergyRecv = 0;
  selfEnergy = 0.0;
  doEnergy = false;
  doVirial = false;
  numStrayAtoms = 0;
}

PmePencilZ::PmePencilZ(CkMigrateMessage *m) {
  NAMD_bug("PmePencilZ cannot be migrated");
  //__sdag_init();
  // setMigratable(false);
  // fftCompute = NULL;
  // pmeTranspose = NULL;
}

PmePencilZ::~PmePencilZ() {
  if (fftCompute != NULL) delete fftCompute;
  if (pmeTranspose != NULL) delete pmeTranspose;
  if (pmeKSpaceCompute != NULL) delete pmeKSpaceCompute;
  delete reduction;
}

void PmePencilZ::initFFT(PmeStartMsg *msg) {
  if (fftCompute == NULL)
    NAMD_bug("PmePencilZ::initFFT, fftCompute not initialized");
  fftCompute->init(msg->data, msg->dataSize,  NULL, 0, Perm_Z_cX_Y, pmeGrid, 1, thisIndex.x, thisIndex.y, 0); 
}

void PmePencilZ::forwardFFT() {
  if (fftCompute == NULL)
    NAMD_bug("PmePencilZ::forwardFFT, fftCompute not initialized");
  fftCompute->forward();
}

void PmePencilZ::backwardFFT() {
  if (fftCompute == NULL)
    NAMD_bug("PmePencilZ::backwardFFT, fftCompute not initialized");
  fftCompute->backward();
}

void PmePencilZ::initBlockSizes() {
  blockSizes.resize(pmeGrid.zBlocks);
  for (int z=0;z < pmeGrid.zBlocks;z++) {
    int i0, i1, j0, j1, k0, k1;
    getBlockDim(pmeGrid, Perm_Z_cX_Y, z, thisIndex.x, thisIndex.y,
      i0, i1, j0, j1, k0, k1);
    int size = (i1-i0+1)*(j1-j0+1)*(k1-k0+1);
    blockSizes[z] = size;
  }
}

void PmePencilZ::forwardDone() {
  if (pmeKSpaceCompute == NULL)
    NAMD_bug("PmePencilZ::forwardDone, pmeKSpaceCompute not initialized");
  pmeKSpaceCompute->solve(lattice, doEnergy, doVirial, fftCompute->getDataDst());
}

void PmePencilZ::submitReductions() {
  if (pmeKSpaceCompute == NULL)
    NAMD_bug("PmePencilZ::submitReductions, pmeKSpaceCompute not initialized");
  double virial[9];
  double energy = pmeKSpaceCompute->getEnergy();
  // fprintf(stderr, "PmePencilZ::submitReductions(), numStrayAtoms %d\n", numStrayAtoms);
  pmeKSpaceCompute->getVirial(virial);
  reduction->item(REDUCTION_ELECT_ENERGY_SLOW) += energy + selfEnergy;
  reduction->item(REDUCTION_VIRIAL_SLOW_XX) += virial[0];
  reduction->item(REDUCTION_VIRIAL_SLOW_XY) += virial[1];
  reduction->item(REDUCTION_VIRIAL_SLOW_XZ) += virial[2];
  reduction->item(REDUCTION_VIRIAL_SLOW_YX) += virial[3];
  reduction->item(REDUCTION_VIRIAL_SLOW_YY) += virial[4];
  reduction->item(REDUCTION_VIRIAL_SLOW_YZ) += virial[5];
  reduction->item(REDUCTION_VIRIAL_SLOW_ZX) += virial[6];
  reduction->item(REDUCTION_VIRIAL_SLOW_ZY) += virial[7];
  reduction->item(REDUCTION_VIRIAL_SLOW_ZZ) += virial[8];
  reduction->item(REDUCTION_STRAY_CHARGE_ERRORS) += numStrayAtoms;
  reduction->submit();
  numStrayAtoms = 0;
}

void PmePencilZ::backwardDone() {
  NAMD_bug("PmePencilZ::backwardDone(), base class method called");
}

void PmePencilZ::recvDataFromY(PmeBlockMsg *msg) {
  NAMD_bug("PmePencilY::recvDataFromY(), base class method called");
}

void PmePencilZ::start() {
  NAMD_bug("PmePencilZ::start(), base class method called");
}

void PmePencilZ::skip() {
  reduction->submit();  
}

#include "PmeSolver.def.h"
