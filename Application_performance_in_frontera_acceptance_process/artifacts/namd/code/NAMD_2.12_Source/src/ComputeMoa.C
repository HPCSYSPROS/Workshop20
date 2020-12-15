/* *
 * ***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
 * ***  The Board of Trustees of the University of Illinois.
 * ***  All rights reserved.
 * **
 * =====================================================================================
 *
 *       Filename:  ComputeMoa.C
 *
 *    Description:  
 *
 *        Version:  Stub File
 *        Created:  03/08/2012 03:55:27 PM
 *       Revision:  
 *       Compiler:  charm++ 
 *
 *         Author:  Christopher B Harrison, Ph.D.
 *         Email:   charris5@gmail.com, char1@illinois.edu, char@ks.uiuc.edu
 *        Company:  
 *
 * =====================================================================================
 */

#include "InfoStream.h"
#include "Node.h"
#include "PDB.h"
#include "PatchMap.h"
#include "PatchMap.inl"
#include "ComputeMoa.h"
#include "ComputeMoaMgr.decl.h"
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
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include "pup_stl.h"

#ifdef CHARM_HAS_MSA

#ifdef OPENATOM_VERSION


void MoaData::pup(PUP::er &p)
{
  p|K1, p|K2, p|K3;
  p|order;
  p|orig_x, p|orig_y, p|orig_z;

  p|num_clients;
  p|num_clients_qh, p|num_clients_bh;

  p|qh; 
  p|sh;
  p|bh;

  p|k1r, p|k2r, p|k3r;
  p|gbqs;
  p|hasLQ;
  p|hasLB; 
  p|hasLS;
}

void MoaData::print()
{
#if 0
  printf("MoaData:\n");
  printf("K1,K2,K3 = %d, %d, %d\n", K1, K2, K3);
  printf("order = %d\n", order);
  printf("origin = %g %g %g\n", orig_x, orig_y, orig_z);
#endif
}

class ComputeMoaMgr : public CBase_ComputeMoaMgr {
public:
  ComputeMoaMgr();                    // entry
  ~ComputeMoaMgr();

  CkCallback doneMoa;

  void initialize(CkQdMsg *);         // entry with message

  void setCompute(ComputeMoa *c) { moaCompute = c; }  // local

  void recvMoaData(const MoaData &);  // entry with parameter marshalling

  void initWorkers(CkQdMsg *);        // entry with message
  void startWorkers(CkQdMsg *);       // entry with message
  
  void recvQ(int g, int numgrids, int nq, float qgrid[nq], CkCallback doneMoa);
  
  void recvB(int K2_start, int K2_end, int K3_start, int K3_end, int K1_len, int K2_len, int K3_len, double bm1[K1_len], double bm2[K2_len], double bm3[K3_len], int order);     // entry with message
  
  const MoaData &getMoaData() const { return moaData; }  // local

private:
  CProxy_ComputeMoaMgr moaProxy;
  ComputeMoa *moaCompute;

  MoaData moaData;

  CProxy_MoaPatch moaPatchProxy;   // 1D chare array 
  CProxy_MoaS moaSProxy;   // 1D chare array 
  
  Moa3Grid::Write qhwrit // write charge grid
  Moa3Grid::Accum qhacc;  // accumulate charge grid 
  Moa3Grid::Read qhread;  // read charge grid 
  Moa3Grid::Write shwrit;  // writes S grid 
  Moa3Grid::Accum shacc;  // accumulate S grid 
  Moa3Grid::Read shread;  // read S grid 
  Moa3Grid::Write bhwrit; // writes b1, b2, b3 values
  Moa1Grid::Read bhread; // read b1, b2, b3 values
};

// class MoaPatch : public CBase_MoaPatch {
// public:
//   MoaPatch(Moa3Grid &qh, Moa3Grid &eh);
//   MoaPatch(CkMigrateMessage *m) { }
//   void compute();
// private:
// };

class MoaS : public CBase_MoaS {
public:
  MoaS(Moa3Grid &q, Moa3Grid &b);
  MoaS(CkMigrateMessage *m) { }
  ~MoaS();
  void compute();
private:
  Moa3Grid qh, bh;
  Moa3Grid::Write shwrite;
  Moa3Grid::Read qhread, bhread;
  int moaS_setup;
  int lia, lja, lka; // local indices of sub-lattice
  int lib, ljb, lkb; // local indices of sub-lattice
  SubmitReduction *reduction;
};


ComputeMoa::ComputeMoa(ComputeID c) :
  ComputeHomePatches(c)
{
#ifdef OPENATOM_VERSION_DEBUG
  CkPrintf("Instantiated ComputeMoa %d on PE: %d, at %f.\n", thisIndex, CkMyPe(), CkWallTimer() );
#endif //OPENATOM_VERSION_DEBUG

  CProxy_ComputeMoaMgr::ckLocalBranch(
      CkpvAccess(BOCclass_group).computeMoaMgr)->setCompute(this);
  SimParameters *simParams = Node::Object()->simParameters;
  // reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);
}

ComputeMoa::~ComputeMoa()
{
}

void ComputeMoa::doWork()
{
#ifdef OPENATOM_VERSION_DEBUG
  CkPrintf("Starting ComputeMoa::doWork() %d on PE: %d, at %f.\n", thisIndex, CkMyPe(), CkWallTimer() );
#endif //OPENATOM_VERSION_DEBUG

  // allocate message

  CProxy_ComputeMoaMgr moaProxy(CkpvAccess(BOCclass_group).computeMoaMgr);
  moaS_Proxy[CkMyPe()].compute(new CkQdMsg);

}

ComputeMoaMgr::ComputeMoaMgr() :
  moaProxy(thisgroup), moaCompute(0)
{
#ifdef OPENATOM_VERSION_DEBUG
  CkPrintf("Instantiating ComputeMoaMgr %d on PE: %d, at %f.\n", thisIndex, CkMyPe(), CkWallTimer() );
#endif //OPENATOM_VERSION_DEBUG
  CkpvAccess(BOCclass_group).computeMoaMgr = thisgroup;
}

ComputeMoaMgr::~ComputeMoaMgr()
{
#ifdef OPENATOM_VERSION_DEBUG
  CkPrintf("Destructing ComputeMoaMgr %d on PE: %d, at %f.\n", thisIndex, CkMyPe(), CkWallTimer() );
#endif //OPENATOM_VERSION_DEBUG

}


void ComputeMoaMgr::initialize(CkQdMsg *msg)
{
#ifdef OPENATOM_VERSION_DEBUG
  CkPrintf("Initializing ComputeMoaMgr %d on PE: %d, at %f.\n", thisIndex, CkMyPe(), CkWallTimer() );
#endif //OPENATOM_VERSION_DEBUG

  delete msg;
  if (CkMyPe() != 0) return;  // initialize only on PE 0, broadcast MoaData

  // initialize MOA here
  SimParameters *simParams = Node::Object()->simParameters;
  Molecule *mol = Node::Object()->molecule;

  // setup grids and vectors
  MoaData &p = moaData;
  // p.hasLQ=0;
  // p.hasLB=0;

  extern int isPmeProcessor(int);
  int i, num_clients=0;
  for (i=0;i<CkNumPes(), i++)
  {
    if (isPmeProcessor(i)) p.num_clients++;
  }

  // RESTORE VDW PARAMETERS 

  // RESTORE MOLECULE OBJECT PARAMETERS fepAtomFlag atoms vdw 0 and charge bCol 

  int done, alldone;


#ifdef MOA_DEBUG
  CkPrintf("global num clients = %d\n", p.num_clients);
  for (i = 0;  i < p.num_clients;  i++) {
    CkPrintf("num clients bh[%d] = %d\n", i, p.num_clients_bh[i]);
    CkPrintf("num clients qh[%d] = %d\n", i, p.num_clients_qh[i]);
    CkPrintf("num clients sh[%d] = %d\n", i, p.num_clients_bh[i]);
  }
#endif

  p.hasLB.resize(p.num_clients);
  p.hasLQ.resize(p.num_clients);
  p.bh.resize(p.num_clients);
  p.qh.resize(p.num_clients);
  p.sh.resize(p.num_clients);

  for (i = 0;  i < p.num_clients;  i++) {
    ia = p.k1r[i].nx;
    ib = p.k1r[i].ny; 
    ja = p.k2r[i].nx;
    jb = p.k2r[i].ny;
    ka = p.k3r[i].nx;
    kb = p.k3r[i].ny;

    p.bh[i] = Moa3Grid(ia, ib, ja, jb, ka, kb, p.num_clients_qh[i]);
    p.qh[i] = Moa3Grid(ia, ib, ja, jb, ka, kb, p.num_clients_qh[i]);
    p.sh[i] = Moa3Grid(ia, ib, ja, jb, ka, kb, p.num_clients_sh[i]);

  }
  moaProxy.recvMoaData(moaData);  // broadcast MoaData to chare group

#ifdef OPENATOM_VERSION_DEBUG
  CkPrintf("Initialized ComputeMoaMgr %d on PE: %d, at %f.\n", thisIndex, CkMyPe(), CkWallTimer() );
#endif //OPENATOM_VERSION_DEBUG

}

void ComputeMoaMgr::recvB(int K2_start, int K2_end, int K3_start, int K3_end, int K1_len, int K2_len, int K3_len, double bm1[K1_len], double bm2[K2_len], double bm3[K3_len], int order)
{
#ifdef OPENATOM_VERSION_DEBUG
  CkPrintf("Started recvB() %d on PE: %d, at %f.\n", thisIndex, CkMyPe(), CkWallTimer() );
#endif //OPENATOM_VERSION_DEBUG

//  const MoaData &p = CProxy_ComputeMoaMgr::ckLocalBranch( CkpvAccess(BOCclass_group).computeMoaMgr)->getMoaData();

  MoaData &p = moaData;
  int j,m,n;

  p.k1r[thisIndex].nx = 0;
  p.k1r[thisIndex].ny = K1_len;
  p.k2r[thisIndex].nx = K2_start;
  p.k2r[thisIndex].ny = K2_end;
  p.k3r[thisIndex].nx = K3_start;
  p.k3r[thisIndex].ny = k3_end;

  p.bh[thisIndex].enroll();
  bhwrit = p.bh[thisIndex].syncToWrite();
  
  for ( j=0; j < K1_len; j++ ) {
    bhwrit.set(j,K2_start,K3_start) = bm1[j];
    for ( m=K2_start; m < K2_end; m++ ) {
      bhwrit.set(j,m,K3_start) = bm2[m];
      for ( n=K3_start; n < K3_end; n++ ) {
        bhwrit.set(j,m,n) = bm3[n];
      }
    }
  }

  bhread = bhwrit.syncToRead();
  
  p.hasLB[thisIndex]=1;
  p.gbqs.nx++;
  
}

void ComputeMoaMgr::recvQ(int g, int numgrids, int nq, float qgrid[nq], CkCallback doneMoa)
{
#ifdef OPENATOM_VERSION_DEBUG
  CkPrintf("Started recvQ() %d on PE: %d, at %f.\n", thisIndex, CkMyPe(), CkWallTimer() );
#endif //OPENATOM_VERSION_DEBUG

//  const MoaData &p = CProxy_ComputeMoaMgr::ckLocalBranch( CkpvAccess(BOCclass_group).computeMoaMgr)->getMoaData();
  MoaData &p = moaData;

  p.qh[g].enroll(p.num_clients_qh[thisIndex]);
  Eqhacc = p.qh[g].getInitialAccum();

  int i,j,m,n;
  
  int K1_len = p.k1r[thisIndex].ny;
  int K2s = p.k2r[thisIndex].nx;
  int K2e = p.k2r[thisIndex].ny;
  int K3s = p.k3r[thisIndex].nx;
  int K3e = p.k3r[thisIndex].ny;
  
  // qhwrit = p.qh[thisIndex].syncToWrite();
  bhread = p.bh[thisIndex].syncToRead();

//  for ( i=0, i < nq, i+2) {
//    qhrwrit.set(i)=qgrid[i];
//    qhiwrit.set(i)=qgrid[i+1];
//  }

  
  for ( j=CkIndex(); j < K1_len; j++ ) {
    qhwrit.set(j,K2s,K3s) = qgrid[j];
    for ( m=K2s; m < K2e; m++ ) {
      bhwrit.set(j,m,K3s) = qgrid[m];
      for ( n=K3s; n < K3e; n++ ) {
        bhwrit.set(j,m,n) = qgrid[n];
      }
    }
  }

 
  // for ( j = nq*g; j <= (nq * g + nq); j++) 
  for ( i = 0; i <= nq ; i++) {
    p.qhwrite[g][i] = qgrid[i];
  }

  if (g = numgrids) {
    p.gbqs.ny++; 
    p.hasLQ[thisIndex]=1;  
  }

}

void ComputeMoaMgr::sendS2OA()
{
  // syncToRead BGrid, pack it
  // syncToRead SGrid, pack it
  // send msg to openatom containing BGrid and SGrid
  if (CkMyPe() != 0) return;

#ifdef OPENATOM_VERSION_DEBUG
  CkPrintf("Started sendS2OA() %d on PE: %d, at %f.\n", thisIndex, CkMyPe(), CkWallTimer() );
#endif //OPENATOM_VERSION_DEBUG

  CProxy_CP_LargeSP_RhoGSpacePlane oaRGSP_Proxy;

  CkCallback resumeMoa(CkIndex_ComputeMoaMgr::recvSGrid(), thishandle);
  oaRGSP_Proxy[CkMyPe()].recvMDSg(g, ns, sg, resumeMoa);

}

void ComputeMoaMgr::recvSGrid(int g, int nq, double oaSg[nq])
{
  if (CkMyPe() != 0) return;
#ifdef OPENATOM_VERSION_DEBUG
  CkPrintf("Started recvSGrid() %d on PE: %d, at %f.\n", thisIndex, CkMyPe(), CkWallTimer() );
#endif //OPENATOM_VERSION_DEBUG

  MoaData &p = moaData;
  numgrids = p.numgrids;


}

void ComputeMoaMgr::initWorkers(CkQdMsg *msg)
{
  delete msg;
  //printf("MOA initWorkers PE=%d\n", CkMyPe());
  if (CkMyPe() != 0) return;
  // PE 0 creates the compute chare arrays

  MoaData &p = moaData;
  int n;
  moaPatchProxy = CProxy_MoaPatch::ckNew();  // create empty chare array
#ifdef OPENATOM_VERSION_DEBUG
  CkPrintf("Instantiated MoaPatch chare %d on PE: %d, at %f.\n", thisIndex, CkMyPe(), CkWallTimer() );
#endif // OPENATOM_VERSION_DEBUG

  moaSProxy = CProxy_MoaS::ckNew(p.qh[0], p.eh[0], p.num_energy_chares);

#ifdef OPENATOM_VERSION_DEBUG
  CkPrintf("Created MoaS chare %d at PE: %d, at %f.\n", thisIndex, CkMyPe(), CkWallTimer() );
#endif // OPENATOM_VERSION_DEBUG
}


void ComputeMoaMgr::startWorkers(CkQdMsg *msg)
{
  delete msg;
  if (CkMyPe() != 0) return;
  // only PE 0 can start workers;
#ifdef OPENATOM_VERSION_DEBUG
  CkPrintf("Starting MoaMgr Workers on chare %d at PE: %d, at %f.\n", thisIndex, CkMyPe(), CkWallTimer() );
#endif // OPENATOM_VERSION_DEBUG

  // moaPatchProxy.compute();
  moaSProxy.compute();
}


//  /* Archaic soln */
//  MoaPatch::MoaPatch(Moa3Grid &qh, Moa3Grid &sh, Moa3Grid &bh)
//  {
//  #ifdef OPENATOM_VERSION_DEBUG
//    CkPrintf("Instantiated MoaPatch %d on PE: %d, at %f.\n", thisIndex, CkMyPe(), CkWallTimer() );
//  #endif // OPENATOM_VERSION_DEBUG
//  
//    const MoaData &p = CProxy_ComputeMoaMgr::ckLocalBranch(CkpvAccess(BOCclass_group).computeMoaMgr)->getMoaData();
//    int pIndex = thisIndex;
//    const Int3 &dg = p.dim_gridcutoff_chares[pIndex];
//    calcS = CProxy_MoaSGrid::ckNew(pIndex, qh, eh, dg.nx, dg.ny, dg.nz);
//  }


MoaS::MoaS(Moa3Grid &qh, Moa3Grid &bh)
{
#ifdef OPENATOM_VERSION_DEBUG
  CkPrintf("Instantiating MoaS %d on PE: %d, at %f.\n", thisIndex, CkMyPe(), CkWallTimer() );
#endif //OPENATOM_VERSION_DEBUG
  
  const MoaData &p = CProxy_ComputeMoaMgr::ckLocalBranch( CkpvAccess(BOCclass_group).computeMoaMgr)->getMoaData(); 

  
  lia = p.k1r.nx;
  lib = p.k1r.ny;
  lja = p.k2r.nx;
  ljb = p.k2r.ny;
  lka = p.k3r.nx;
  lkb = p.k3r.ny;
  buff_len = (lib - lia) + (ljb - lja) + (lkb - lka);
  qhbuff = new float[buff_len];
  bhbuff = new int[buff_len];

  moaS_setup = 0;

  reduction = ReductionMgr:Object()->willSubmit(REDUCTIONS_BASIC);
}

MoaS::~MoaS()
{
#ifdef OPENATOM_VERSION_DEBUG
  CkPrintf("Destructing MoaS %d on PE: %d, at %f.\n", thisIndex, CkMyPe(), CkWallTimer() );
#endif //OPENATOM_VERSION_DEBUG

}

MoaS::compute()
{
#ifdef OPENATOM_VERSION_DEBUG
  CkPrintf("Starting MoaS compute %d on PE: %d, at %f.\n", thisIndex, CkMyPe(), CkWallTimer() );
#endif //OPENATOM_VERSION_DEBUG

  const MoaData &p = CProxy_ComputeMoaMgr::ckLocalBranch( CkpvAccess(BOCclass_group).computeMoaMgr)->getMoaData(); 
   
  if ( ! moaS_setup ) {
    bh.enroll(p.num_clients);
    qh[0].enroll(p.num_clients);
    sh.enroll(p.num_clients);
   
    bhread = bh.syncToRead();
    qhacc = qh.getInitialAccum();
    qhread = qhacc.syncToRead();
    shwrit = sh.getInitialWrite();
    moaS_setup = 1;
  }

  bhread = bh.syncToRead();

  int j,m,n,index;

  for (index = 0, j = lia; j <= lib; j++) {
    for (m = lja; m <= ljb; m++) {
      for (n = lka; n <= lkb; n++, index++) {
        bhbuff[index] = bhread.get(j,m,n);
      }
    }
  }
  
  qhacc = qhread.syncToEAccum();
  qhread = qhacc.syncToRead();
  for (index = 0, j = lia; j <= lib; j++) {
    for (m = lja; m <= ljb; m++) {
      for (n = lka; n <= lkb; n++, index++) {
        qhbuff[index] = qhread.get(j,m,n);
      }
    }
  }

  shwrit = sh.getInitialWrite();
  for (index = 0, j = lia; j <= lib; j++) {
    for (m = lja; m <= ljb; m++) {
      for (n = lka; n <= lkb; n++, index++) {
        shwrite.set(j,m,n) = qhbuff[index] * bhbuff[index] * bhbuff[index+1] * bhbuff[index+2]
      }
    }
  }

}


#include "ComputeMoaMgr.def.h"

#endif // OPENTATOM_VERSION

#endif // CHARM_HAS_MSA


