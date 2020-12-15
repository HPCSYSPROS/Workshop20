/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   HomePatch is the key distributed source/sink of Atom data
   including positions, velocities and forces applied
*/

#ifndef HOMEPATCH_H
#define HOMEPATCH_H

#include "charm++.h"

#include "NamdTypes.h"
#include "Patch.h"
#include "PatchMap.h"

#include "MigrateAtomsMsg.h"
#include "main.h"
#include "common.h"
#include "Migration.h"
#include "Settle.h"

#include <string>
#include <map>

class RegisterProxyMsg;
class UnregisterProxyMsg;
class ProxyResultVarsizeMsg;
class ProxyResultMsg;
class ProxyCombinedResultRawMsg;
class Sequencer;
class SubmitReduction;
class ProxyGBISP1ResultMsg;
class ProxyGBISP2ResultMsg;
class CheckpointAtomsMsg;
class ExchangeAtomsMsg;

class ProxyNodeAwareSpanningTreeMsg;

class ComputeQMMgr;

class HomePatch : public Patch {
  friend class PatchMgr;
  friend class Sequencer;
  friend class ComputeGlobal;

private: 

  HomePatch(PatchID, int atomCnt);
  // for PatchMgr to use only
  HomePatch(PatchID, FullAtomList&);

//  HomePatch(PatchID, int atomCnt);

  void reinitAtoms(FullAtomList&);
  ScaledPosition min, max, center;
  BigReal aAwayDist, bAwayDist, cAwayDist;

  Bool doAtomUpdate;  // atom changes other than migration

  //Note: If new proxies are added to this HomePatch
  // after load balancing, and it is not the immediate step
  // after atom migration (where ProxyAllMsg will be sent), 
  // then the CompAtomExt list has to be resent with the 
  // ProxyDataMsg (the normal proxy msg when atoms don't 
  // migrate), otherwise, program will crash without such 
  // information when doing force calculations --Chao Mei
  Bool isNewProxyAdded;
  int numGBISP1Arrived, numGBISP2Arrived, numGBISP3Arrived;
  bool phase1BoxClosedCalled;
  bool phase2BoxClosedCalled;
  bool phase3BoxClosedCalled;

public:
  ~HomePatch();

  // Message from ProxyPatch (via ProxyMgr) which registers its existence
  void registerProxy(RegisterProxyMsg *);
  // opposite of above
  void unregisterProxy(UnregisterProxyMsg *);

  // ProxyPatch sends Forces back to here (via ProxyMgr)  
  void receiveResults(ProxyResultVarsizeMsg *msg);
  void receiveResults(ProxyResultMsg *msg);     
  //gbis receiving results from intermediate phases
  void receiveResult(ProxyGBISP1ResultMsg *msg);//after P1
  void receiveResult(ProxyGBISP2ResultMsg *msg);//after P2
  
  //direct function calls, not as entry methods
  void receiveResults(ProxyCombinedResultRawMsg *msg);

  // AtomMigration messages passes from neighbor HomePatches to here.
  void depositMigration(MigrateAtomsMsg *);

  // Bind a Sequencer to this HomePatch
  void useSequencer(Sequencer *sequencerPtr);
  // start simulation over this Patch of atoms
  void runSequencer(void);
  
  //--------------------------------------------------------------------
  // methods for Sequencer to use
  //

  // Signal HomePatch that positions stored are to be now to be used
  void positionsReady(int doMigration=0);
  int marginViolations;

  // methods to implement integration
  void saveForce(const int ftag = Results::normal);
  void addForceToMomentum(const BigReal, const int ftag = Results::normal,
				const int useSaved = 0);
  void addForceToMomentum3(const BigReal timestep1, const int ftag1, const int useSaved1,
    const BigReal timestep2, const int ftag2, const int useSaved2,
    const BigReal timestep3, const int ftag3, const int useSaved3);
  void addVelocityToPosition(const BigReal);

  // impose hard wall constraint on Drude bond length
  int hardWallDrude(const BigReal, Tensor *virial, SubmitReduction *);

  // methods for rigidBonds
  struct RattleList {
    int ig;
    int icnt;
  };

  std::vector<int> settleList;
  std::vector<RattleList> rattleList;
  std::vector<RattleParam> rattleParam;
  std::vector<int> noconstList;

  bool rattleListValid;

  // Array to store new positions and velocities. Allocated in "buildRattleList" to size numAtoms
  std::vector<Vector> velNew;
  std::vector<Vector> posNew;

  void addRattleForce(const BigReal invdt, Tensor& wc);

  void buildRattleList();
  int rattle1old(const BigReal, Tensor *virial, SubmitReduction *);
  int rattle1(const BigReal, Tensor *virial, SubmitReduction *);
  void rattle2(const BigReal, Tensor *virial);
  void minimize_rattle2(const BigReal, Tensor *virial, bool forces=false);

  // methods for mollified impluse (MOLLY)
  void mollyAverage();
  void mollyMollify(Tensor *virial);
//  Bool average(Vector qtilde[],const Vector q[],BigReal lambda[],const int n,const int m, const BigReal imass[], const BigReal length2[], const int ial[], const int ilb[], const Vector qji[], const BigReal tolf, const int ntrial);
//  void mollify(Vector qtilde[],const Vector q0[],const BigReal lambda[], Vector force[],const int n, const int m, const BigReal imass[],const int ial[],const int ibl[],const Vector refab[]); 
  
  // BEGIN LA
  void loweAndersenVelocities();
  void loweAndersenFinish();
  // END LA

  void setGBISIntrinsicRadii();
  void gbisComputeAfterP1();//calculate bornRad
  void gbisComputeAfterP2();//calculate dHdrPrefix or self energies
  void gbisP2Ready();
  void gbisP3Ready();

  //LCPO
  void setLcpoType();

  // methods for CONTRA, etc
  void checkpoint(void);
  void revert(void);

  void exchangeCheckpoint(int scriptTask, int &bpc);
  void recvCheckpointReq(int task, const char *key, int replica, int pe);
  void recvCheckpointLoad(CheckpointAtomsMsg *msg);
  void recvCheckpointStore(CheckpointAtomsMsg *msg);
  void recvCheckpointAck();
  int checkpoint_task;
  struct checkpoint_t {
    Lattice lattice;
    int berendsenPressure_count;
    int numAtoms;
    ResizeArray<FullAtom> atoms;
  };
  std::map<std::string,checkpoint_t*> checkpoints;

  // replica exchange
  void exchangeAtoms(int scriptTask);
  void recvExchangeReq(int req);
  void recvExchangeMsg(ExchangeAtomsMsg *msg);
  int exchange_dst;
  int exchange_src;
  int exchange_req;
  ExchangeAtomsMsg *exchange_msg;

  // methods for QM (ExtForces replacement)
  void replaceForces(ExtForce *f);

  void qmSwapAtoms();
  
  // load-balancing trigger
  void submitLoadStats(int timestep);

  // for ComputeHomePatches
  FullAtomList &getAtomList() { return (atom); }

#ifdef NODEAWARE_PROXY_SPANNINGTREE
  // build spanning tree for proxy nodes
  void buildNodeAwareSpanningTree(void);
  void setupChildrenFromProxySpanningTree();
#else
    // build spanning tree for proxy nodes
  void buildSpanningTree(void);
#endif

  void sendNodeAwareSpanningTree();
  void recvNodeAwareSpanningTree(ProxyNodeAwareSpanningTreeMsg *msg);

  void sendSpanningTree();
  void recvSpanningTree(int *t, int n);


  void sendProxies();

#if USE_TOPOMAP 
  int findSubroots(int dim, int* subroots, int psize, int* pidscopy);
#endif

  LDObjHandle ldObjHandle;
protected:
  virtual void boxClosed(int);

  // Internal Atom Migration methods and data
  void doPairlistCheck();
  void doGroupSizeCheck();
  void doMarginCheck();
  void doAtomMigration();
  int inMigration;
  int numMlBuf;
  MigrateAtomsMsg *msgbuf[PatchMap::MaxOneAway];
  
private:
  // Store of Atom-wise variables
  FullAtomList  atom;
  ForceList f_saved[Results::maxNumForces];
  ExtForce *replacementForces;

  CudaAtomList cudaAtomList;

  // DMK - Atom Separation (water vs. non-water)
  #if NAMD_SeparateWaters != 0
    FullAtomList tempAtom;  // A temporary array used to sort waters
                            //   from non-waters in the atom array
    void separateAtoms();   // Function to separate the atoms currently in atoms.
    void mergeAtomList(FullAtomList &al);  // Function to combine and separate
                                           //   the atoms in al with atoms.
  #endif


  // checkpointed state
  FullAtomList  checkpoint_atom;
  Lattice  checkpoint_lattice;

  // DMK - Atom Separation (water vs. non-water)
  #if NAMD_SeparateWaters != 0
    int checkpoint_numWaterAtoms;
  #endif


  // checkPairlist data
  CompAtomList doPairlistCheck_positions;
  Lattice doPairlistCheck_lattice;
  BigReal doPairlistCheck_newTolerance;

  // MOLLY data
  ResizeArray<BigReal> molly_lambda;
  
  // List of Proxies
  NodeIDList proxy;
  
  Sequencer  *sequencer;

  // Needed for initialization
  int patchMapRead;
  void readPatchMap();

  // Atom Migration internals
  int allMigrationIn;
  int migrationSuspended;
  int patchMigrationCounter;
  int numNeighbors;
  MigrationInfo realInfo[PatchMap::MaxOneAway];
  MigrationInfo *mInfo[3][3][3];

#ifdef NODEAWARE_PROXY_SPANNINGTREE
  //the whole spanning tree for all the proxies this home patch has
  proxyTreeNodeList ptnTree;
  //the immediate children (recording pe ids) containing two parts: 
  //one part of them all belong to the physical node this home patch
  // resides on; the other part of pes belong to all external nodes.
  /* Moved to Patch.h */ 
  //int *children;
  //int numChild;
#else
  NodeIDList tree;              // the whole tree
  int *child;	// spanning tree of proxies - immediate children
  int nChild;
#endif

  // Cached settle1 parameters
  int settle_initialized;
  BigReal settle_mOrmT; BigReal settle_mHrmT; BigReal settle_ra;
  BigReal settle_rb; BigReal settle_rc; BigReal settle_rra;

  // Drude lone pairs
  void redistrib_lonepair_forces(const int, Tensor *);

  // PLF -- for TIP4P
  //void redistrib_tip4p_force(Vector&, Vector&, Vector&, Vector&, int, Tensor*);
  void redistrib_tip4p_forces(const int, Tensor*);
  void tip4_omrepos(Vector*, Vector*, Vector*, BigReal);
  void init_tip4();

  // Drude SWM4
  void redistrib_swm4_forces(const int, Tensor*);
  void swm4_omrepos(Vector*, Vector*, Vector*, BigReal);
  void init_swm4();

  // reposition a lone pair using its host atoms and additional parameters
  void reposition_lonepair(
      Vector& ri, const Vector& rj, const Vector& rk, const Vector& rl,
      Real distance, Real angle, Real dihedral);

  void reposition_all_lonepairs(void);

  // general redistribution of lone pair forces to host atoms
  void redistrib_lp_force(
      Vector& fi, Vector& fj, Vector& fk, Vector& fl,
      const Vector& ri, const Vector& rj, const Vector& rk, const Vector& rl,
      Tensor *virial, int midpt);

  // use for both TIP4P and SWM4 water
  void redistrib_lp_water_force(
      Vector& f_ox, Vector& f_h1, Vector& f_h2, Vector& f_lp,
      const Vector& p_ox, const Vector& p_h1, const Vector& p_h2,
      const Vector& p_lp, Tensor *virial);

  BigReal r_om, r_ohc;
  void write_tip4_props(void);

  int isProxyChanged;

#if CMK_PERSISTENT_COMM
  PersistentHandle *localphs;
  int nphs;
#endif
};

#endif

