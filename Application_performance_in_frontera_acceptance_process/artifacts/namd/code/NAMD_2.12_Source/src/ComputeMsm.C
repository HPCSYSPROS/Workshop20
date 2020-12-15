/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "InfoStream.h"
#include "Node.h"
#include "PDB.h"
#include "PatchMap.h"
#include "PatchMap.inl"
#include "AtomMap.h"
#include "ComputeMsm.h"
#include "PatchMgr.h"
#include "Molecule.h"
#include "ReductionMgr.h"
#include "ComputeMgr.h"
#include "ComputeMgr.decl.h"
#include "Debug.h"
#include "SimParameters.h"
#include "WorkDistrib.h"
#include "Priorities.h"
#include "varsizemsg.h"
//#include "ckmulticast.h"
#include <stdio.h>
#include "MsmMap.h"

// use multicast reduction of grids from sections of MsmGridCutoff
#define MSM_REDUCE_GRID
//#undef MSM_REDUCE_GRID

// use the decomposition of grid cutoff to create more work units
#define MSM_GRID_CUTOFF_DECOMP
//#undef MSM_GRID_CUTOFF_DECOMP

// skip over pairs of blocks that do not actually interact
#define MSM_SKIP_TOO_DISTANT_BLOCKS
//#undef MSM_SKIP_TOO_DISTANT_BLOCKS

// skip over pairs of blocks whose overlap is beyond nonzero gc sphere
// this search is more expensive than MSM_SKIP_TOO_DISTANT_BLOCKS
// and does not eliminate many block pairs
#define MSM_SKIP_BEYOND_SPHERE
//#undef MSM_SKIP_BEYOND_SPHERE

// node aware mapping of chare arrays
#define MSM_NODE_MAPPING
//#undef MSM_NODE_MAPPING

#define MSM_NODE_MAPPING_STATS
#undef MSM_NODE_MAPPING_STATS

// top of hierarchy calculates smaller blocks of charge to 
// unfolded image blocks of potential, up to the desired block size,
// then sums the unfolded images of potential back into the 
// actual potential block, thereby greatly reducing the number of 
// block pairs that would otherwise be scheduled
#define MSM_FOLD_FACTOR
//#undef MSM_FOLD_FACTOR

// report timings for compute routines
// for diagnostic purposes only
#define MSM_TIMING
#undef MSM_TIMING

// report profiling for compute routines
// for diagnostic purposes only
#define MSM_PROFILING
#undef MSM_PROFILING

// use fixed size grid message
// XXX probably does not work anymore
#define MSM_FIXED_SIZE_GRID_MSG
#undef MSM_FIXED_SIZE_GRID_MSG

// turn off computation
// for diagnostic purposes only
//#define MSM_COMM_ONLY

// print diagnostics for memory alignment (for grid cutoff calculation)
// for diagnostic purposes only
#define DEBUG_MEMORY_ALIGNMENT
#undef DEBUG_MEMORY_ALIGNMENT


//
// This is the main message that gets passed between compute chares.
// It is used to bundle blocks of charge (sendUp and send to MsmGridCutoff) 
// and blocks of potential (sendAcross, sendDown, and sendPatch).  
//
// Higher priority has a numerically lower value.  
//
// The priorities are set as follows:
//
//   sendUp priority = level+1
//
//   (send to MsmGridCutoff) and sendAcross priority
//       = nlevels + 2*(nlevels - level) - 1
//
//   sendDown and sendPatch priority
//       = nlevels + 2*(nlevels - level)
//
// This puts the priority on going up the hierarchy before going across 
// and puts the priority on finishing the top levels and down before 
// finishing the lower levels.
//

class GridMsg : public CkMcastBaseMsg, public CMessage_GridMsg {
  public:
    char *gdata;
    int idnum;
    int nlower_i;
    int nlower_j;
    int nlower_k;
    int nextent_i;
    int nextent_j;
    int nextent_k;
    int nbytes;
    int seqnum;  // sequence number is used for message priority

    // put a grid into an allocated message to be sent
    template <class T>
    void put(const msm::Grid<T>& g, int id, int seq) {
      idnum = id;
      nlower_i = g.lower().i;
      nlower_j = g.lower().j;
      nlower_k = g.lower().k;
      nextent_i = g.extent().i;
      nextent_j = g.extent().j;
      nextent_k = g.extent().k;
      nbytes = g.data().len()*sizeof(T);
      seqnum = seq;
      memcpy(gdata, g.data().buffer(), nbytes);
    }

    // get the grid from a received message
    template <class T>
    void get(msm::Grid<T>& g, int& id, int& seq) {
      id = idnum;
      g.set(nlower_i, nextent_i, nlower_j, nextent_j,
          nlower_k, nextent_k);
      seq = seqnum;
      ASSERT(g.data().len()*sizeof(T) == nbytes);
      memcpy(g.data().buffer(), gdata, nbytes);
    }
};


class MsmBlockProxyMsg : public CMessage_MsmBlockProxyMsg {
  public:
    enum { maxlevels = 32 };
    char msmBlockProxyData[maxlevels*sizeof(CProxy_MsmBlock)];
    int nlevels;

    // put an array into an allocated message to be sent
    void put(const msm::Array<CProxy_MsmBlock>& a) {
      nlevels = a.len();
      if (nlevels > maxlevels) {
        NAMD_die("Exceeded maximum number of MSM levels\n");
      }
      memcpy(msmBlockProxyData, a.buffer(), nlevels*sizeof(CProxy_MsmBlock));
    }

    // get the array from a received message
    void get(msm::Array<CProxy_MsmBlock>& a) {
      a.resize(nlevels);
      memcpy(a.buffer(), msmBlockProxyData, nlevels*sizeof(CProxy_MsmBlock));
    }
};


class MsmC1HermiteBlockProxyMsg : public CMessage_MsmC1HermiteBlockProxyMsg {
  public:
    enum { maxlevels = 32 };
    char msmBlockProxyData[maxlevels*sizeof(CProxy_MsmC1HermiteBlock)];
    int nlevels;

    // put an array into an allocated message to be sent
    void put(const msm::Array<CProxy_MsmC1HermiteBlock>& a) {
      nlevels = a.len();
      if (nlevels > maxlevels) {
        NAMD_die("Exceeded maximum number of MSM levels\n");
      }
      memcpy(msmBlockProxyData, a.buffer(),
          nlevels*sizeof(CProxy_MsmC1HermiteBlock));
    }

    // get the array from a received message
    void get(msm::Array<CProxy_MsmC1HermiteBlock>& a) {
      a.resize(nlevels);
      memcpy(a.buffer(), msmBlockProxyData,
          nlevels*sizeof(CProxy_MsmC1HermiteBlock));
    }
};


class MsmGridCutoffProxyMsg : public CMessage_MsmGridCutoffProxyMsg {
  public:
    char msmGridCutoffProxyData[sizeof(CProxy_MsmGridCutoff)];

    // put proxy into an allocated message to be sent
    void put(const CProxy_MsmGridCutoff *p) {
      memcpy(msmGridCutoffProxyData, p, sizeof(CProxy_MsmGridCutoff));
    }

    // get the proxy from a received message
    void get(CProxy_MsmGridCutoff *p) {
      memcpy(p, msmGridCutoffProxyData, sizeof(CProxy_MsmGridCutoff));
    }
};


class MsmC1HermiteGridCutoffProxyMsg :
  public CMessage_MsmC1HermiteGridCutoffProxyMsg
{
  public:
    char msmGridCutoffProxyData[sizeof(CProxy_MsmC1HermiteGridCutoff)];

    // put proxy into an allocated message to be sent
    void put(const CProxy_MsmC1HermiteGridCutoff *p) {
      memcpy(msmGridCutoffProxyData, p,
          sizeof(CProxy_MsmC1HermiteGridCutoff));
    }

    // get the proxy from a received message
    void get(CProxy_MsmC1HermiteGridCutoff *p) {
      memcpy(p, msmGridCutoffProxyData,
          sizeof(CProxy_MsmC1HermiteGridCutoff));
    }
};


class MsmGridCutoffInitMsg : public CMessage_MsmGridCutoffInitMsg {
  public:
    msm::BlockIndex qhBlockIndex;  // charge block index
    msm::BlockSend ehBlockSend;    // potential block sending address
    MsmGridCutoffInitMsg(const msm::BlockIndex& i, const msm::BlockSend& b)
      : qhBlockIndex(i), ehBlockSend(b) { }
};


class MsmGridCutoffSetupMsg :
  public CkMcastBaseMsg, public CMessage_MsmGridCutoffSetupMsg
{
  public:
    char msmBlockElementProxyData[sizeof(CProxyElement_MsmBlock)];

    // put proxy into an allocated message to be sent
    void put(
        const CProxyElement_MsmBlock *q //,
        ) {
      memcpy(msmBlockElementProxyData, q, sizeof(CProxyElement_MsmBlock));
    }

    // get the proxy from a received message
    void get(
        CProxyElement_MsmBlock *q //,
        ) {
      memcpy(q, msmBlockElementProxyData, sizeof(CProxyElement_MsmBlock));
    }
};


class MsmC1HermiteGridCutoffSetupMsg :
  public CkMcastBaseMsg, public CMessage_MsmC1HermiteGridCutoffSetupMsg
{
  public:
    char msmBlockElementProxyData[sizeof(CProxyElement_MsmC1HermiteBlock)];

    // put proxy into an allocated message to be sent
    void put(
        const CProxyElement_MsmC1HermiteBlock *q //,
        ) {
      memcpy(msmBlockElementProxyData, q,
          sizeof(CProxyElement_MsmC1HermiteBlock));
    }

    // get the proxy from a received message
    void get(
        CProxyElement_MsmC1HermiteBlock *q //,
        ) {
      memcpy(q, msmBlockElementProxyData,
          sizeof(CProxyElement_MsmC1HermiteBlock));
    }
};


// Used only when MSM_TIMING is defined
class MsmTimer : public CBase_MsmTimer {
  public:
    enum { ANTERP=0, INTERP, RESTRICT, PROLONGATE, GRIDCUTOFF, COMM, MAX };

    MsmTimer() {
      for (int i = 0;  i < MAX;  i++)  timing[i] = 0;
    }
    void done(double tm[], int n) {
      for (int i = 0;  i < MAX;  i++)  timing[i] = tm[i];
      print();
    }
    void print() {
      CkPrintf("MSM timings:\n");
      CkPrintf("   anterpolation   %8.6f sec\n", timing[ANTERP]);
      CkPrintf("   interpolation   %8.6f sec\n", timing[INTERP]);
      CkPrintf("   restriction     %8.6f sec\n", timing[RESTRICT]);
      CkPrintf("   prolongation    %8.6f sec\n", timing[PROLONGATE]);
      CkPrintf("   grid cutoff     %8.6f sec\n", timing[GRIDCUTOFF]);
      CkPrintf("   communication   %8.6f sec\n", timing[COMM]);
    }

    double timing[MAX];
};


// Used only when MSM_PROFILING is defined
class MsmProfiler : public CBase_MsmProfiler {
  public:
    enum { MAX = MSM_MAX_BLOCK_SIZE+1 };

    MsmProfiler() {
      for (int i = 0;  i < MAX;  i++)  xloopcnt[i] = 0;
    }
    void done(int lc[], int n) {
      for (int i = 0;  i < MAX;  i++)  xloopcnt[i] = lc[i];
      print();
    }
    void print() {
      int sum = 0;
      for (int i = 0;  i < MAX;  i++)  sum += xloopcnt[i];
      CkPrintf("MSM profiling:\n");
      CkPrintf("   total executions of inner loop:   %d\n", sum);
      for (int i = 0;  i < MAX;  i++) {
        CkPrintf("   executing %d times:   %d  (%5.2f%%)\n",
            i, xloopcnt[i], 100*double(xloopcnt[i])/sum);
      }
    }

    int xloopcnt[MAX];
};


// used with PriorityQueue
// when determining work mapped to node or PE
struct WorkIndex {
  float work;
  int index;
  WorkIndex() : work(0), index(0) { }
  WorkIndex(float w, int i) : work(w), index(i) { }
  int operator<=(const WorkIndex& wn) {
    return (work <= wn.work);
  }
};


//////////////////////////////////////////////////////////////////////////////
//
//  ComputeMsmMgr
//  chare group containing MSM parameters and constants;
//  one chare object per PE
//

class ComputeMsmMgr : public CBase_ComputeMsmMgr {
  friend struct msm::PatchData;
  friend class MsmBlock;
  //friend class MsmGridCutoff;
  friend class MsmBlockMap;
  friend class MsmGridCutoffMap;

public:
  ComputeMsmMgr();                    // entry
  ~ComputeMsmMgr();

  void initialize(MsmInitMsg *);      // entry with message
private:
  void initialize2();                 // split in two
public:

  void recvMsmBlockProxy(MsmBlockProxyMsg *);  // entry with message
  void recvMsmGridCutoffProxy(MsmGridCutoffProxyMsg *);  // entry with message

  void recvMsmC1HermiteBlockProxy(MsmC1HermiteBlockProxyMsg *);
    // entry with message
  void recvMsmC1HermiteGridCutoffProxy(MsmC1HermiteGridCutoffProxyMsg *);
    // entry with message

  void update(CkQdMsg *);             // entry with message

  void compute(msm::Array<int>& patchIDList);
                                      // called by local ComputeMsm object

  void addPotential(GridMsg *);  // entry with message
  void doneCompute();  // called by each local patch

#ifdef MSM_TIMING
  void initTiming() {
    for (int i = 0;  i < MsmTimer::MAX;  i++)  msmTiming[i] = 0;
    cntTiming = 0;
  }
  // every local object being timed should call this during initialization
  void addTiming() {
    numTiming++;
  }
  // object calls before being migrated
  void subtractTiming() {
    numTiming--;
  }
  void doneTiming() {
    if (++cntTiming >= numTiming) {
      CkCallback cb(CkReductionTarget(MsmTimer, done), msmTimer);
      contribute(MsmTimer::MAX*sizeof(double), msmTiming,
          CkReduction::sum_double, cb);
      initTiming();
    }
  }
#endif

#ifdef MSM_PROFILING
  void initProfiling() {
    for (int i = 0;  i < MsmProfiler::MAX;  i++)  xLoopCnt[i] = 0;
    cntProfiling = 0;
  }
  // every local object being profiled should call this during initialization
  void addProfiling() {
    numProfiling++;
  }
  // object calls before being migrated
  void subtractProfiling() {
    numProfiling--;
  }
  void doneProfiling() {
    if (++cntProfiling >= numProfiling) {
      CkCallback cb(CkReductionTarget(MsmProfiler, done), msmProfiler);
      contribute(MsmProfiler::MAX*sizeof(int), xLoopCnt,
          CkReduction::sum_int, cb);
      initProfiling();  // reset accumulators for next visit
    }
  }
#endif

  void setCompute(ComputeMsm *c) { msmCompute = c;  c->setMgr(this); } // local

  msm::PatchPtrArray& patchPtrArray() { return patchPtr; }

  msm::Map& mapData() { return map; }

  int numLevels() const { return nlevels; }

  // sign(n) = -1 if n < 0,  0 if n == 0,  or  1 if n > 0
  static inline int sign(int n) {
    return (n < 0 ? -1 : (n > 0 ? 1 : 0));
  }

//private:
  void setup_hgrid_1d(BigReal len, BigReal& hh, int& nn,
      int& ia, int& ib, int isperiodic);
  void setup_periodic_blocksize(int& bsize, int n);

  CProxy_ComputeMsmMgr msmProxy;
  ComputeMsm *msmCompute;

  msm::Array<CProxy_MsmBlock> msmBlock;
  msm::Array<CProxy_MsmC1HermiteBlock> msmC1HermiteBlock;

  CProxy_MsmGridCutoff msmGridCutoff;
  CProxy_MsmC1HermiteGridCutoff msmC1HermiteGridCutoff;
  int numGridCutoff;  // length of msmGridCutoff chare array

  msm::Map map;

  // find patch by patchID
  // array is length number of patches, initialized to NULL
  // allocate PatchData for only those patches on this PE
  msm::PatchPtrArray patchPtr;

  // allocate subgrid used for receiving message data in addPotential()
  // and sending on to PatchData::addPotential()
  msm::Grid<Float> subgrid;
  msm::Grid<C1Vector> subgrid_c1hermite;

#ifdef MSM_NODE_MAPPING
  msm::Array<int> blockAssign;
  msm::Array<int> gcutAssign;
  //msm::Array<int> nodecnt;
  int blockFlatIndex(int level, int i, int j, int k) {
    int n = 0;
    for (int l = 0;  l < level;  l++) {
      n += map.blockLevel[l].nn();
    }
    return (n + map.blockLevel[level].flatindex(i,j,k));
  }
  float calcBlockWork(const msm::BlockDiagram& b) {
    // XXX ratio of work for MsmBlock to MsmGridCutoff?
    const float scalingFactor = 3;
    const int volumeFullBlock = map.bsx[0] * map.bsy[0] * map.bsz[0];
    msm::Ivec gn;
    if (approx == C1HERMITE) {
      gn = map.gc_c1hermite[0].extent();
    }
    else {
      gn = map.gc[0].extent();
    }
    const int volumeFullCutoff = (map.bsx[0] + gn.i - 1) *
      (map.bsy[0] + gn.j - 1) * (map.bsz[0] + gn.k - 1);
    msm::Ivec n = b.nrange.extent();
    int volumeBlock = n.i * n.j * n.k;
    msm::Ivec nc = b.nrangeCutoff.extent();
    int volumeCutoff = nc.i * nc.j * nc.k;
    return( scalingFactor * (float(volumeBlock) / volumeFullBlock) *
        (float(volumeCutoff) / volumeFullCutoff) );
  }
  float calcGcutWork(const msm::BlockSend& bs) {
    const int volumeFullBlock = map.bsx[0] * map.bsy[0] * map.bsz[0];
    msm::Ivec n = bs.nrange_wrap.extent();;
    int volumeBlock = n.i * n.j * n.k;
    return( float(volumeBlock) / volumeFullBlock );
  }
#endif

  // sum local virial factors
  msm::Grid<Float> gvsum;
  int numVirialContrib;
  int cntVirialContrib;
  enum { VXX=0, VXY, VXZ, VYY, VYZ, VZZ, VMAX };
  Float virial[VMAX];

  void initVirialContrib() {
    gvsum.reset(0);
    cntVirialContrib = 0;
  }
  void addVirialContrib() {
    numVirialContrib++;
  }
  void subtractVirialContrib() {
    numVirialContrib--;
  }
  void doneVirialContrib() {
    if (++cntVirialContrib >= numVirialContrib) {
      // reduce all gvsum contributions into virial tensor
      for (int n = 0;  n < VMAX;  n++) { virial[n] = 0; }
      int ia = gvsum.ia();
      int ib = gvsum.ib();
      int ja = gvsum.ja();
      int jb = gvsum.jb();
      int ka = gvsum.ka();
      int kb = gvsum.kb();
      for (int k = ka;  k <= kb;  k++) {
        for (int j = ja;  j <= jb;  j++) {
          for (int i = ia;  i <= ib;  i++) {
            Float cu = Float(i);
            Float cv = Float(j);
            Float cw = Float(k);
            Float c = gvsum(i,j,k);
            Float vx = cu*hufx + cv*hvfx + cw*hwfx;
            Float vy = cu*hufy + cv*hvfy + cw*hwfy;
            Float vz = cu*hufz + cv*hvfz + cw*hwfz;
            virial[VXX] -= c * vx * vx;
            virial[VXY] -= c * vx * vy;
            virial[VXZ] -= c * vx * vz;
            virial[VYY] -= c * vy * vy;
            virial[VYZ] -= c * vy * vz;
            virial[VZZ] -= c * vz * vz;
          }
        }
      }
      initVirialContrib();
    }
  }

#ifdef MSM_TIMING
  CProxy_MsmTimer msmTimer;
  double msmTiming[MsmTimer::MAX];
  int numTiming;  // total number of objects being timed
  int cntTiming;  // count the objects as they provide timing results
  CkCallback *cbTiming;
#endif

#ifdef MSM_PROFILING
  CProxy_MsmProfiler msmProfiler;
  int xLoopCnt[MsmProfiler::MAX];
  int numProfiling;  // total number of objects being profiled
  int cntProfiling;  // count the objects as they provide profiling results
  CkCallback *cbProfiling;
#endif

  Vector c, u, v, w;    // rescaled center and lattice vectors
  Vector ru, rv, rw;    // row vectors to transform to unit space
  int ispu, ispv, ispw; // is periodic along u, v, w?

  Lattice lattice;      // keep local copy of lattice
  ScaledPosition smin;  // keep min values for non-periodic dimensions
  ScaledPosition smax;  // keep max values for non-periodic dimensions
  BigReal gridspacing;  // preferred grid spacing
  BigReal padding;      // padding for non-periodic boundaries
  BigReal gridScalingFactor;  // scaling for Hermite interpolation
  BigReal a;            // cutoff distance
  BigReal hxlen, hylen, hzlen;  // first level grid spacings along basis vectors
  BigReal hxlen_1, hylen_1, hzlen_1;  // inverses of grid spacings
  Vector hu, hv, hw;    // first level grid spacing vectors
  Float hufx, hufy, hufz, hvfx, hvfy, hvfz, hwfx, hwfy, hwfz;
  int nhx, nhy, nhz;    // number of h spacings that cover cell
  int approx;           // ID for approximation
  int split;            // ID for splitting
  int nlevels;          // number of grid levels
  int dispersion;       // calculating dispersion forces?
  BigReal gzero;        // self energy factor from splitting

  Vector sglower;       // lower corner of grid in scaled space
                        // corresponds to index (0,0,0)

  BigReal shx, shy, shz;  // grid spacings in scaled space
  BigReal shx_1, shy_1, shz_1;
  Vector sx_shx;          // row vector to transform interpolated force x
  Vector sy_shy;          // row vector to transform interpolated force y
  Vector sz_shz;          // row vector to transform interpolated force z
  Float srx_x, srx_y, srx_z;  // float version of sx_shx
  Float sry_x, sry_y, sry_z;  // float version of sy_shy
  Float srz_x, srz_y, srz_z;  // float version of sz_shz

  int s_edge;
  int omega;

  enum Approx { CUBIC=0, QUINTIC, QUINTIC2,
    SEPTIC, SEPTIC3, NONIC, NONIC4, C1HERMITE, NUM_APPROX };

  enum Split { TAYLOR2=0, TAYLOR3, TAYLOR4,
    TAYLOR5, TAYLOR6, TAYLOR7, TAYLOR8,
    TAYLOR2_DISP, TAYLOR3_DISP, TAYLOR4_DISP, TAYLOR5_DISP,
    TAYLOR6_DISP, TAYLOR7_DISP, TAYLOR8_DISP, NUM_SPLIT };

  enum {
    // Approximation formulas with up to degree 9 polynomials.
    MAX_POLY_DEGREE = 9,

    // Max stencil length for polynomial approximation.
    MAX_NSTENCIL_SIZE = (2*MAX_POLY_DEGREE + 1),

    // Max stencil length when skipping zeros
    // (almost half entries are zero for interpolating polynomials).
    MAX_NSTENCIL_SKIP_ZERO = (MAX_POLY_DEGREE + 2),

    // Number of scalar approximation formulaes
    NUM_APPROX_FORMS = (NONIC4 - CUBIC) + 1
  };

  // Degree of polynomial basis function Phi.
  static const int PolyDegree[NUM_APPROX];

  // The stencil array lengths below.
  static const int Nstencil[NUM_APPROX];

  // Index offsets from the stencil-centered grid element, to get
  // to the correct contributing grid element.
  static const int IndexOffset[NUM_APPROX][MAX_NSTENCIL_SKIP_ZERO];

  // The grid transfer stencils for the non-factored restriction and
  // prolongation procedures.
  static const Float PhiStencil[NUM_APPROX_FORMS][MAX_NSTENCIL_SKIP_ZERO];

  // Calculate the smoothing function and its derivative:
  // g(R) and (d/dR)g(R), where R=r/a.
  // Use double precision for calculating the MSM constant weights 
  // and coefficients.  The resulting coefficents to be used in 
  // the repeatedly called algorithm are stored in single precision.
  static void splitting(BigReal& g, BigReal& dg, BigReal r_a, int _split) {
    BigReal s = r_a * r_a;  // s = (r/a)^2, assuming 0 <= s <= 1
    switch (_split) {
      case TAYLOR2:
        g = 1 + (s-1)*(-1./2 + (s-1)*(3./8));
        dg = (2*r_a)*(-1./2 + (s-1)*(3./4));
        break;
      case TAYLOR3:
        g = 1 + (s-1)*(-1./2 + (s-1)*(3./8 + (s-1)*(-5./16)));
        dg = (2*r_a)*(-1./2 + (s-1)*(3./4 + (s-1)*(-15./16)));
        break;
      case TAYLOR4:
        g = 1 + (s-1)*(-1./2 + (s-1)*(3./8 + (s-1)*(-5./16
                + (s-1)*(35./128))));
        dg = (2*r_a)*(-1./2 + (s-1)*(3./4 + (s-1)*(-15./16
                + (s-1)*(35./32))));
        break;
      case TAYLOR5:
        g = 1 + (s-1)*(-1./2 + (s-1)*(3./8 + (s-1)*(-5./16
                + (s-1)*(35./128 + (s-1)*(-63./256)))));
        dg = (2*r_a)*(-1./2 + (s-1)*(3./4 + (s-1)*(-15./16
                + (s-1)*(35./32 + (s-1)*(-315./256)))));
        break;
      case TAYLOR6:
        g = 1 + (s-1)*(-1./2 + (s-1)*(3./8 + (s-1)*(-5./16
                + (s-1)*(35./128 + (s-1)*(-63./256
                    + (s-1)*(231./1024))))));
        dg = (2*r_a)*(-1./2 + (s-1)*(3./4 + (s-1)*(-15./16
                + (s-1)*(35./32 + (s-1)*(-315./256
                    + (s-1)*(693./512))))));
        break;
      case TAYLOR7:
        g = 1 + (s-1)*(-1./2 + (s-1)*(3./8 + (s-1)*(-5./16
            + (s-1)*(35./128 + (s-1)*(-63./256
                + (s-1)*(231./1024 + (s-1)*(-429./2048)))))));
        dg = (2*r_a)*(-1./2 + (s-1)*(3./4 + (s-1)*(-15./16
                + (s-1)*(35./32 + (s-1)*(-315./256
                    + (s-1)*(693./512 + (s-1)*(-3003./2048)))))));
        break;
      case TAYLOR8:
        g = 1 + (s-1)*(-1./2 + (s-1)*(3./8 + (s-1)*(-5./16
                + (s-1)*(35./128 + (s-1)*(-63./256
                    + (s-1)*(231./1024 + (s-1)*(-429./2048
                        + (s-1)*(6435./32768))))))));
        dg = (2*r_a)*(-1./2 + (s-1)*(3./4 + (s-1)*(-15./16
                + (s-1)*(35./32 + (s-1)*(-315./256
                    + (s-1)*(693./512 + (s-1)*(-3003./2048
                        + (s-1)*(6435./4096))))))));
        break;
      case TAYLOR2_DISP:
        g = 1 + (s-1)*(-3 + (s-1)*(6));
        dg = (2*r_a)*(-3 + (s-1)*(12));
        break;
      case TAYLOR3_DISP:
        g = 1 + (s-1)*(-3 + (s-1)*(6 + (s-1)*(-10)));
        dg = (2*r_a)*(-3 + (s-1)*(12 + (s-1)*(-30)));
        break;
      case TAYLOR4_DISP:
        g = 1 + (s-1)*(-3 + (s-1)*(6 + (s-1)*(-10 + (s-1)*(15))));
        dg = (2*r_a)*(-3 + (s-1)*(12 + (s-1)*(-30 + (s-1)*(60))));
        break;
      case TAYLOR5_DISP:
        g = 1 + (s-1)*(-3 + (s-1)*(6 + (s-1)*(-10
                + (s-1)*(15 + (s-1)*(-21)))));
        dg = (2*r_a)*(-3 + (s-1)*(12 + (s-1)*(-30
                + (s-1)*(60 + (s-1)*(-105)))));
        break;
      case TAYLOR6_DISP:
        g = 1 + (s-1)*(-3 + (s-1)*(6 + (s-1)*(-10
                + (s-1)*(15 + (s-1)*(-21 + (s-1)*(28))))));
        dg = (2*r_a)*(-3 + (s-1)*(12 + (s-1)*(-30
                + (s-1)*(60 + (s-1)*(-105 + (s-1)*(168))))));
        break;
      case TAYLOR7_DISP:
        g = 1 + (s-1)*(-3 + (s-1)*(6 + (s-1)*(-10
                + (s-1)*(15 + (s-1)*(-21 + (s-1)*(28
                      + (s-1)*(-36)))))));
        dg = (2*r_a)*(-3 + (s-1)*(12 + (s-1)*(-30
                + (s-1)*(60 + (s-1)*(-105 + (s-1)*(168
                      + (s-1)*(-252)))))));
        break;
      case TAYLOR8_DISP:
        g = 1 + (s-1)*(-3 + (s-1)*(6 + (s-1)*(-10
                + (s-1)*(15 + (s-1)*(-21 + (s-1)*(28
                      + (s-1)*(-36 + (s-1)*(45))))))));
        dg = (2*r_a)*(-3 + (s-1)*(12 + (s-1)*(-30
                + (s-1)*(60 + (s-1)*(-105 + (s-1)*(168
                      + (s-1)*(-252 + (s-1)*(360))))))));
        break;
      default:
        NAMD_die("Unknown MSM splitting.");
    } // switch
  } // splitting()

  void stencil_1d(Float phi[], Float t) {
    switch (approx) {
      case CUBIC:
        phi[0] = 0.5f * (1 - t) * (2 - t) * (2 - t);
        t--;
        phi[1] = (1 - t) * (1 + t - 1.5f * t * t);
        t--;
        phi[2] = (1 + t) * (1 - t - 1.5f * t * t);
        t--;
        phi[3] = 0.5f * (1 + t) * (2 + t) * (2 + t);
        break;
      case QUINTIC:
        phi[0] = (1.f/24) * (1-t) * (2-t) * (3-t) * (3-t) * (4-t);
        t--;
        phi[1] = (1-t) * (2-t) * (3-t) * ((1.f/6)
            + t * (0.375f - (5.f/24)*t));
        t--;
        phi[2] = (1-t*t) * (2-t) * (0.5f + t * (0.25f - (5.f/12)*t));
        t--;
        phi[3] = (1-t*t) * (2+t) * (0.5f - t * (0.25f + (5.f/12)*t));
        t--;
        phi[4] = (1+t) * (2+t) * (3+t) * ((1.f/6)
            - t * (0.375f + (5.f/24)*t));
        t--;
        phi[5] = (1.f/24) * (1+t) * (2+t) * (3+t) * (3+t) * (4+t);
        break;
      case QUINTIC2:
        phi[0] = (1.f/24) * (3-t) * (3-t) * (3-t) * (t-2) * (5*t-8);
        t--;
        phi[1] = (-1.f/24) * (2-t) * (t-1) * (-48+t*(153+t*(-114+t*25)));
        t--;
        phi[2] = (1.f/12) * (1-t) * (12+t*(12+t*(-3+t*(-38+t*25))));
        t--;
        phi[3] = (1.f/12) * (1+t) * (12+t*(-12+t*(-3+t*(38+t*25))));
        t--;
        phi[4] = (-1.f/24) * (2+t) * (t+1) * (48+t*(153+t*(114+t*25)));
        t--;
        phi[5] = (1.f/24) * (3+t) * (3+t) * (3+t) * (t+2) * (5*t+8);
        break;
      case SEPTIC:
        phi[0] = (-1.f/720)*(t-1)*(t-2)*(t-3)*(t-4)*(t-4)*(t-5)*(t-6);
        t--;
        phi[1] = (1.f/720)*(t-1)*(t-2)*(t-3)*(t-4)*(t-5)*(-6+t*(-20+7*t));
        t--;
        phi[2] = (-1.f/240)*(t*t-1)*(t-2)*(t-3)*(t-4)*(-10+t*(-12+7*t));
        t--;
        phi[3] = (1.f/144)*(t*t-1)*(t*t-4)*(t-3)*(-12+t*(-4+7*t));
        t--;
        phi[4] = (-1.f/144)*(t*t-1)*(t*t-4)*(t+3)*(-12+t*(4+7*t));
        t--;
        phi[5] = (1.f/240)*(t*t-1)*(t+2)*(t+3)*(t+4)*(-10+t*(12+7*t));
        t--;
        phi[6] = (-1.f/720)*(t+1)*(t+2)*(t+3)*(t+4)*(t+5)*(-6+t*(20+7*t));
        t--;
        phi[7] = (1.f/720)*(t+1)*(t+2)*(t+3)*(t+4)*(t+4)*(t+5)*(t+6);
        break;
      case SEPTIC3:
        phi[0] = (3632.f/5) + t*((-7456.f/5) + t*((58786.f/45) + t*(-633
                + t*((26383.f/144) + t*((-22807.f/720) + t*((727.f/240)
                      + t*(-89.f/720)))))));
        t--;
        phi[1] = -440 + t*((25949.f/20) + t*((-117131.f/72) + t*((2247.f/2)
                + t*((-66437.f/144) + t*((81109.f/720) + t*((-727.f/48)
                      + t*(623.f/720)))))));
        t--;
        phi[2] = (138.f/5) + t*((-8617.f/60) + t*((12873.f/40) + t*((-791.f/2)
                + t*((4557.f/16) + t*((-9583.f/80) + t*((2181.f/80)
                      + t*(-623.f/240)))))));
        t--;
        phi[3] = 1 + t*t*((-49.f/36) + t*t*((-959.f/144) + t*((2569.f/144)
                + t*((-727.f/48) + t*(623.f/144)))));
        t--;
        phi[4] = 1 + t*t*((-49.f/36) + t*t*((-959.f/144) + t*((-2569.f/144)
                + t*((-727.f/48) + t*(-623.f/144)))));
        t--;
        phi[5] = (138.f/5) + t*((8617.f/60) + t*((12873.f/40) + t*((791.f/2)
                + t*((4557.f/16) + t*((9583.f/80) + t*((2181.f/80)
                      + t*(623.f/240)))))));
        t--;
        phi[6] = -440 + t*((-25949.f/20) + t*((-117131.f/72) + t*((-2247.f/2)
                + t*((-66437.f/144) + t*((-81109.f/720) + t*((-727.f/48)
                      + t*(-623.f/720)))))));
        t--;
        phi[7] = (3632.f/5) + t*((7456.f/5) + t*((58786.f/45) + t*(633
                + t*((26383.f/144) + t*((22807.f/720) + t*((727.f/240)
                      + t*(89.f/720)))))));
        break;
      case NONIC:
        phi[0] = (-1.f/40320)*(t-8)*(t-7)*(t-6)*(t-5)*(t-5)*(t-4)*(t-3)*
          (t-2)*(t-1);
        t--;
        phi[1] = (1.f/40320)*(t-7)*(t-6)*(t-5)*(t-4)*(t-3)*(t-2)*(t-1)*
          (-8+t*(-35+9*t));
        t--;
        phi[2] = (-1.f/10080)*(t-6)*(t-5)*(t-4)*(t-3)*(t-2)*(t-1)*(t+1)*
          (-14+t*(-25+9*t));
        t--;
        phi[3] = (1.f/1440)*(t-5)*(t-4)*(t-3)*(t-2)*(t-1)*(t+1)*(t+2)*
          (-6+t*(-5+3*t));
        t--;
        phi[4] = (-1.f/2880)*(t-4)*(t-3)*(t-2)*(t-1)*(t+1)*(t+2)*(t+3)*
          (-20+t*(-5+9*t));
        t--;
        phi[5] = (1.f/2880)*(t-3)*(t-2)*(t-1)*(t+1)*(t+2)*(t+3)*(t+4)*
          (-20+t*(5+9*t));
        t--;
        phi[6] = (-1.f/1440)*(t-2)*(t-1)*(t+1)*(t+2)*(t+3)*(t+4)*(t+5)*
          (-6+t*(5+3*t));
        t--;
        phi[7] = (1.f/10080)*(t-1)*(t+1)*(t+2)*(t+3)*(t+4)*(t+5)*(t+6)*
          (-14+t*(25+9*t));
        t--;
        phi[8] = (-1.f/40320)*(t+1)*(t+2)*(t+3)*(t+4)*(t+5)*(t+6)*(t+7)*
          (-8+t*(35+9*t));
        t--;
        phi[9] = (1.f/40320)*(t+1)*(t+2)*(t+3)*(t+4)*(t+5)*(t+5)*(t+6)*
          (t+7)*(t+8);
        break;
      case NONIC4:
      { // begin grouping to define local variables
        double Tphi[10];
        double T=t;
        Tphi[0] = 439375./7+T*(-64188125./504+T*(231125375./2016
              +T*(-17306975./288+T*(7761805./384+T*(-2895587./640
                    +T*(129391./192+T*(-259715./4032+T*(28909./8064
                          +T*(-3569./40320)))))))));
        T--;
        Tphi[1] = -56375+T*(8314091./56+T*(-49901303./288+T*(3763529./32
                +T*(-19648027./384+T*(9469163./640+T*(-545977./192
                      +T*(156927./448+T*(-28909./1152
                          +T*(3569./4480)))))))));
        T--;
        Tphi[2] = 68776./7+T*(-1038011./28+T*(31157515./504+T*(-956669./16
                +T*(3548009./96+T*(-2422263./160+T*(197255./48
                      +T*(-19959./28+T*(144545./2016
                          +T*(-3569./1120)))))))));
        T--;
        Tphi[3] = -154+T*(12757./12+T*(-230123./72+T*(264481./48
                +T*(-576499./96+T*(686147./160+T*(-96277./48
                      +T*(14221./24+T*(-28909./288+T*(3569./480)))))))));
        T--;
        Tphi[4] = 1+T*T*(-205./144+T*T*(91./192+T*(-6181./320
                +T*(6337./96+T*(-2745./32+T*(28909./576
                      +T*(-3569./320)))))));
        T--;
        Tphi[5] = 1+T*T*(-205./144+T*T*(91./192+T*(6181./320
                +T*(6337./96+T*(2745./32+T*(28909./576
                      +T*(3569./320)))))));
        T--;
        Tphi[6] = -154+T*(-12757./12+T*(-230123./72+T*(-264481./48
                +T*(-576499./96+T*(-686147./160+T*(-96277./48
                      +T*(-14221./24+T*(-28909./288+T*(-3569./480)))))))));
        T--;
        Tphi[7] = 68776./7+T*(1038011./28+T*(31157515./504+T*(956669./16
                +T*(3548009./96+T*(2422263./160+T*(197255./48
                      +T*(19959./28+T*(144545./2016+T*(3569./1120)))))))));
        T--;
        Tphi[8] = -56375+T*(-8314091./56+T*(-49901303./288+T*(-3763529./32
                +T*(-19648027./384+T*(-9469163./640+T*(-545977./192
                      +T*(-156927./448+T*(-28909./1152
                          +T*(-3569./4480)))))))));
        T--;
        Tphi[9] = 439375./7+T*(64188125./504+T*(231125375./2016
              +T*(17306975./288+T*(7761805./384+T*(2895587./640
                    +T*(129391./192+T*(259715./4032+T*(28909./8064
                          +T*(3569./40320)))))))));
        for (int i=0;  i < 10;  i++) {
          phi[i] = Float(Tphi[i]);
        }
      } // end grouping to define local variables
        break;
      default:
        NAMD_die("Unknown MSM approximation.");
    } // switch
  } // stencil_1d()

  void d_stencil_1d(Float dphi[], Float phi[], Float t, Float h_1) {
    switch (approx) {
      case CUBIC:
        phi[0] = 0.5f * (1 - t) * (2 - t) * (2 - t);
        dphi[0] = (1.5f * t - 2) * (2 - t) * h_1;
        t--;
        phi[1] = (1 - t) * (1 + t - 1.5f * t * t);
        dphi[1] = (-5 + 4.5f * t) * t * h_1;
        t--;
        phi[2] = (1 + t) * (1 - t - 1.5f * t * t);
        dphi[2] = (-5 - 4.5f * t) * t * h_1;
        t--;
        phi[3] = 0.5f * (1 + t) * (2 + t) * (2 + t);
        dphi[3] = (1.5f * t + 2) * (2 + t) * h_1;
        break;
      case QUINTIC:
        phi[0] = (1.f/24) * (1-t) * (2-t) * (3-t) * (3-t) * (4-t);
        dphi[0] = ((-1.f/24) * ((3-t) * (3-t) * (14 + t * (-14 + 3*t))
              + 2 * (1-t) * (2-t) * (3-t) * (4-t))) * h_1;
        t--;
        phi[1] = (1-t) * (2-t) * (3-t) * ((1.f/6)
            + t * (0.375f - (5.f/24)*t));
        dphi[1] = (-((1.f/6) + t * (0.375f - (5.f/24)*t)) *
            (11 + t * (-12 + 3*t)) + (1-t) * (2-t) * (3-t) *
            (0.375f - (5.f/12)*t)) * h_1;
        t--;
        phi[2] = (1-t*t) * (2-t) * (0.5f + t * (0.25f - (5.f/12)*t));
        dphi[2] = (-(0.5f + t * (0.25f - (5.f/12)*t)) * (1 + t * (4 - 3*t))
            + (1-t*t) * (2-t) * (0.25f - (5.f/6)*t)) * h_1;
        t--;
        phi[3] = (1-t*t) * (2+t) * (0.5f - t * (0.25f + (5.f/12)*t));
        dphi[3] = ((0.5f + t * (-0.25f - (5.f/12)*t)) * (1 + t * (-4 - 3*t))
            - (1-t*t) * (2+t) * (0.25f + (5.f/6)*t)) * h_1;
        t--;
        phi[4] = (1+t) * (2+t) * (3+t) * ((1.f/6)
            - t * (0.375f + (5.f/24)*t));
        dphi[4] = (((1.f/6) + t * (-0.375f - (5.f/24)*t)) *
            (11 + t * (12 + 3*t)) - (1+t) * (2+t) * (3+t) *
            (0.375f + (5.f/12)*t)) * h_1;
        t--;
        phi[5] = (1.f/24) * (1+t) * (2+t) * (3+t) * (3+t) * (4+t);
        dphi[5] = ((1.f/24) * ((3+t) * (3+t) * (14 + t * (14 + 3*t))
              + 2 * (1+t) * (2+t) * (3+t) * (4+t))) * h_1;
        break;
      case QUINTIC2:
        phi[0] = (1.f/24) * (3-t) * (3-t) * (3-t) * (t-2) * (5*t-8);
        dphi[0] = ((1.f/24) * (3-t) * (3-t) * ((3-t)*(5*t-8)
              - 3*(t-2)*(5*t-8) + 5*(t-2)*(3-t))) * h_1;
        t--;
        phi[1] = (-1.f/24) * (2-t) * (t-1) * (-48+t*(153+t*(-114+t*25)));
        dphi[1] = ((-1.f/24) * ((2-t)*(-48+t*(153+t*(-114+t*25)))
              - (t-1)* (-48+t*(153+t*(-114+t*25)))
              + (2-t)*(t-1)*(153+t*(-228+t*75)))) * h_1;
        t--;
        phi[2] = (1.f/12) * (1-t) * (12+t*(12+t*(-3+t*(-38+t*25))));
        dphi[2] = ((1.f/12) * (-(12+t*(12+t*(-3+t*(-38+t*25))))
              + (1-t)*(12+t*(-6+t*(-114+t*100))))) * h_1;
        t--;
        phi[3] = (1.f/12) * (1+t) * (12+t*(-12+t*(-3+t*(38+t*25))));
        dphi[3] = ((1.f/12) * ((12+t*(-12+t*(-3+t*(38+t*25))))
              + (1+t)*(-12+t*(-6+t*(114+t*100))))) * h_1;
        t--;
        phi[4] = (-1.f/24) * (2+t) * (t+1) * (48+t*(153+t*(114+t*25)));
        dphi[4] = ((-1.f/24) * ((2+t)*(48+t*(153+t*(114+t*25)))
              + (t+1)* (48+t*(153+t*(114+t*25)))
              + (2+t)*(t+1)*(153+t*(228+t*75)))) * h_1;
        t--;
        phi[5] = (1.f/24) * (3+t) * (3+t) * (3+t) * (t+2) * (5*t+8);
        dphi[5] = ((1.f/24) * (3+t) * (3+t) * ((3+t)*(5*t+8)
              + 3*(t+2)*(5*t+8) + 5*(t+2)*(3+t))) * h_1;
        break;
      case SEPTIC:
        phi[0] = (-1.f/720)*(t-1)*(t-2)*(t-3)*(t-4)*(t-4)*(t-5)*(t-6);
        dphi[0] = (-1.f/720)*(t-4)*(-1944+t*(3644+t*(-2512+t*(807
                  +t*(-122+t*7))))) * h_1;
        t--;
        phi[1] = (1.f/720)*(t-1)*(t-2)*(t-3)*(t-4)*(t-5)*(-6+t*(-20+7*t));
        dphi[1] = (1.f/720)*(756+t*(-9940+t*(17724+t*(-12740+t*(4445
                    +t*(-750+t*49)))))) * h_1;
        t--;
        phi[2] = (-1.f/240)*(t*t-1)*(t-2)*(t-3)*(t-4)*(-10+t*(-12+7*t));
        dphi[2] = (-1.f/240)*(-28+t*(1260+t*(-756+t*(-1260+t*(1365
                    +t*(-450+t*49)))))) * h_1;
        t--;
        phi[3] = (1.f/144)*(t*t-1)*(t*t-4)*(t-3)*(-12+t*(-4+7*t));
        dphi[3] = (1.f/144)*t*(-560+t*(84+t*(644+t*(-175
                  +t*(-150+t*49))))) * h_1;
        t--;
        phi[4] = (-1.f/144)*(t*t-1)*(t*t-4)*(t+3)*(-12+t*(4+7*t));
        dphi[4] = (-1.f/144)*t*(560+t*(84+t*(-644+t*(-175
                  +t*(150+t*49))))) * h_1;
        t--;
        phi[5] = (1.f/240)*(t*t-1)*(t+2)*(t+3)*(t+4)*(-10+t*(12+7*t));
        dphi[5] = (1.f/240)*(-28+t*(-1260+t*(-756+t*(1260+t*(1365
                    +t*(450+t*49)))))) * h_1;
        t--;
        phi[6] = (-1.f/720)*(t+1)*(t+2)*(t+3)*(t+4)*(t+5)*(-6+t*(20+7*t));
        dphi[6] = (-1.f/720)*(756+t*(9940+t*(17724+t*(12740+t*(4445
                    +t*(750+t*49)))))) * h_1;
        t--;
        phi[7] = (1.f/720)*(t+1)*(t+2)*(t+3)*(t+4)*(t+4)*(t+5)*(t+6);
        dphi[7] = (1.f/720)*(t+4)*(1944+t*(3644+t*(2512+t*(807
                  +t*(122+t*7))))) * h_1;
        break;
      case SEPTIC3:
        phi[0] = (3632.f/5) + t*((-7456.f/5) + t*((58786.f/45) + t*(-633
                + t*((26383.f/144) + t*((-22807.f/720) + t*((727.f/240)
                      + t*(-89.f/720)))))));
        dphi[0] = ((-7456.f/5) + t*((117572.f/45) + t*(-1899
                + t*((26383.f/36) + t*((-22807.f/144) + t*((727.f/40)
                      + t*(-623.f/720))))))) * h_1;
        t--;
        phi[1] = -440 + t*((25949.f/20) + t*((-117131.f/72) + t*((2247.f/2)
                + t*((-66437.f/144) + t*((81109.f/720) + t*((-727.f/48)
                      + t*(623.f/720)))))));
        dphi[1] = ((25949.f/20) + t*((-117131.f/36) + t*((6741.f/2)
                + t*((-66437.f/36) + t*((81109.f/144) + t*((-727.f/8)
                      + t*(4361.f/720))))))) * h_1;
        t--;
        phi[2] = (138.f/5) + t*((-8617.f/60) + t*((12873.f/40) + t*((-791.f/2)
                + t*((4557.f/16) + t*((-9583.f/80) + t*((2181.f/80)
                      + t*(-623.f/240)))))));
        dphi[2] = ((-8617.f/60) + t*((12873.f/20) + t*((-2373.f/2)
                + t*((4557.f/4) + t*((-9583.f/16) + t*((6543.f/40)
                      + t*(-4361.f/240))))))) * h_1;
        t--;
        phi[3] = 1 + t*t*((-49.f/36) + t*t*((-959.f/144) + t*((2569.f/144)
                + t*((-727.f/48) + t*(623.f/144)))));
        dphi[3] = (t*((-49.f/18) + t*t*((-959.f/36) + t*((12845.f/144)
                  + t*((-727.f/8) + t*(4361.f/144)))))) * h_1;
        t--;
        phi[4] = 1 + t*t*((-49.f/36) + t*t*((-959.f/144) + t*((-2569.f/144)
                + t*((-727.f/48) + t*(-623.f/144)))));
        dphi[4] = (t*((-49.f/18) + t*t*((-959.f/36) + t*((-12845.f/144)
                  + t*((-727.f/8) + t*(-4361.f/144)))))) * h_1;
        t--;
        phi[5] = (138.f/5) + t*((8617.f/60) + t*((12873.f/40) + t*((791.f/2)
                + t*((4557.f/16) + t*((9583.f/80) + t*((2181.f/80)
                      + t*(623.f/240)))))));
        dphi[5] = ((8617.f/60) + t*((12873.f/20) + t*((2373.f/2)
                + t*((4557.f/4) + t*((9583.f/16) + t*((6543.f/40)
                      + t*(4361.f/240))))))) * h_1;
        t--;
        phi[6] = -440 + t*((-25949.f/20) + t*((-117131.f/72) + t*((-2247.f/2)
                + t*((-66437.f/144) + t*((-81109.f/720) + t*((-727.f/48)
                      + t*(-623.f/720)))))));
        dphi[6] = ((-25949.f/20) + t*((-117131.f/36) + t*((-6741.f/2)
                + t*((-66437.f/36) + t*((-81109.f/144) + t*((-727.f/8)
                      + t*(-4361.f/720))))))) * h_1;
        t--;
        phi[7] = (3632.f/5) + t*((7456.f/5) + t*((58786.f/45) + t*(633
                + t*((26383.f/144) + t*((22807.f/720) + t*((727.f/240)
                      + t*(89.f/720)))))));
        dphi[7] = ((7456.f/5) + t*((117572.f/45) + t*(1899
                + t*((26383.f/36) + t*((22807.f/144) + t*((727.f/40)
                      + t*(623.f/720))))))) * h_1;
        break;
      case NONIC:
        phi[0] = (-1.f/40320)*(t-8)*(t-7)*(t-6)*(t-5)*(t-5)*(t-4)*(t-3)*
          (t-2)*(t-1);
        dphi[0] = (-1.f/40320)*(t-5)*(-117648+t*(256552+t*(-221416
                +t*(99340+t*(-25261+t*(3667+t*(-283+t*9)))))))*h_1;
        t--;
        phi[1] = (1.f/40320)*(t-7)*(t-6)*(t-5)*(t-4)*(t-3)*(t-2)*(t-1)*
          (-8+t*(-35+9*t));
        dphi[1] = (1.f/40320)*(71856+t*(-795368+t*(1569240+t*(-1357692
                  +t*(634725+t*(-172116+t*(27090+t*(-2296+t*81))))))))*h_1;
        t--;
        phi[2] = (-1.f/10080)*(t-6)*(t-5)*(t-4)*(t-3)*(t-2)*(t-1)*(t+1)*
          (-14+t*(-25+9*t));
        dphi[2] = (1.f/10080)*(3384+t*(-69080+t*(55026
                +t*(62580+t*(-99225+t*(51660+t*(-13104+t*(1640
                          +t*(-81)))))))))*h_1;
        t--;
        phi[3] = (1.f/1440)*(t-5)*(t-4)*(t-3)*(t-2)*(t-1)*(t+1)*(t+2)*
          (-6+t*(-5+3*t));
        dphi[3] = (1.f/1440)*(72+t*(-6344+t*(2070
                +t*(7644+t*(-4725+t*(-828+t*(1260+t*(-328+t*27))))))))*h_1;
        t--;
        phi[4] = (-1.f/2880)*(t-4)*(t-3)*(t-2)*(t-1)*(t+1)*(t+2)*(t+3)*
          (-20+t*(-5+9*t));
        dphi[4] = (-1.f/2880)*t*(10792+t*(-972+t*(-12516
                +t*(2205+t*(3924+t*(-882+t*(-328+t*81)))))))*h_1;
        t--;
        phi[5] = (1.f/2880)*(t-3)*(t-2)*(t-1)*(t+1)*(t+2)*(t+3)*(t+4)*
          (-20+t*(5+9*t));
        dphi[5] = (1.f/2880)*t*(-10792+t*(-972+t*(12516
                +t*(2205+t*(-3924+t*(-882+t*(328+t*81)))))))*h_1;
        t--;
        phi[6] = (-1.f/1440)*(t-2)*(t-1)*(t+1)*(t+2)*(t+3)*(t+4)*(t+5)*
          (-6+t*(5+3*t));
        dphi[6] = (1.f/1440)*(-72+t*(-6344+t*(-2070
                +t*(7644+t*(4725+t*(-828+t*(-1260+t*(-328+t*(-27)))))))))*h_1;
        t--;
        phi[7] = (1.f/10080)*(t-1)*(t+1)*(t+2)*(t+3)*(t+4)*(t+5)*(t+6)*
          (-14+t*(25+9*t));
        dphi[7] = (1.f/10080)*(-3384+t*(-69080+t*(-55026
                +t*(62580+t*(99225+t*(51660+t*(13104+t*(1640+t*81))))))))*h_1;
        t--;
        phi[8] = (-1.f/40320)*(t+1)*(t+2)*(t+3)*(t+4)*(t+5)*(t+6)*(t+7)*
          (-8+t*(35+9*t));
        dphi[8] = (-1.f/40320)*(71856+t*(795368+t*(1569240
                +t*(1357692+t*(634725+t*(172116+t*(27090+t*(2296
                          +t*81))))))))*h_1;
        t--;
        phi[9] = (1.f/40320)*(t+1)*(t+2)*(t+3)*(t+4)*(t+5)*(t+5)*(t+6)*
          (t+7)*(t+8);
        dphi[9] = (1.f/40320)*(t+5)*(117648+t*(256552+t*(221416
                +t*(99340+t*(25261+t*(3667+t*(283+t*9)))))))*h_1;
        break;
      case NONIC4:
      { // begin grouping to define local variables
        double Tphi[10], Tdphi[10];
        double T=t;
        Tphi[0] = 439375./7+T*(-64188125./504+T*(231125375./2016
              +T*(-17306975./288+T*(7761805./384+T*(-2895587./640
                    +T*(129391./192+T*(-259715./4032+T*(28909./8064
                          +T*(-3569./40320)))))))));
        Tdphi[0] = (-64188125./504+T*(231125375./1008
              +T*(-17306975./96+T*(7761805./96+T*(-2895587./128
                    +T*(129391./32+T*(-259715./576+T*(28909./1008
                          +T*(-3569./4480))))))))) * h_1;
        T--;
        Tphi[1] = -56375+T*(8314091./56+T*(-49901303./288+T*(3763529./32
                +T*(-19648027./384+T*(9469163./640+T*(-545977./192
                      +T*(156927./448+T*(-28909./1152
                          +T*(3569./4480)))))))));
        Tdphi[1] = (8314091./56+T*(-49901303./144+T*(11290587./32
                +T*(-19648027./96+T*(9469163./128+T*(-545977./32
                      +T*(156927./64+T*(-28909./144
                          +T*(32121./4480))))))))) * h_1;
        T--;
        Tphi[2] = 68776./7+T*(-1038011./28+T*(31157515./504+T*(-956669./16
                +T*(3548009./96+T*(-2422263./160+T*(197255./48
                      +T*(-19959./28+T*(144545./2016
                          +T*(-3569./1120)))))))));
        Tdphi[2] = (-1038011./28+T*(31157515./252+T*(-2870007./16
                +T*(3548009./24+T*(-2422263./32+T*(197255./8
                      +T*(-19959./4+T*(144545./252
                          +T*(-32121./1120))))))))) * h_1;
        T--;
        Tphi[3] = -154+T*(12757./12+T*(-230123./72+T*(264481./48
                +T*(-576499./96+T*(686147./160+T*(-96277./48
                      +T*(14221./24+T*(-28909./288+T*(3569./480)))))))));
        Tdphi[3] = (12757./12+T*(-230123./36+T*(264481./16
                +T*(-576499./24+T*(686147./32+T*(-96277./8
                      +T*(99547./24+T*(-28909./36
                          +T*(10707./160))))))))) * h_1;
        T--;
        Tphi[4] = 1+T*T*(-205./144+T*T*(91./192+T*(-6181./320
                +T*(6337./96+T*(-2745./32+T*(28909./576
                      +T*(-3569./320)))))));
        Tdphi[4] = T*(-205./72+T*T*(91./48+T*(-6181./64
                +T*(6337./16+T*(-19215./32+T*(28909./72
                      +T*(-32121./320))))))) * h_1;
        T--;
        Tphi[5] = 1+T*T*(-205./144+T*T*(91./192+T*(6181./320
                +T*(6337./96+T*(2745./32+T*(28909./576
                      +T*(3569./320)))))));
        Tdphi[5] = T*(-205./72+T*T*(91./48+T*(6181./64
                +T*(6337./16+T*(19215./32+T*(28909./72
                      +T*(32121./320))))))) * h_1;
        T--;
        Tphi[6] = -154+T*(-12757./12+T*(-230123./72+T*(-264481./48
                +T*(-576499./96+T*(-686147./160+T*(-96277./48
                      +T*(-14221./24+T*(-28909./288+T*(-3569./480)))))))));
        Tdphi[6] = (-12757./12+T*(-230123./36+T*(-264481./16
                +T*(-576499./24+T*(-686147./32+T*(-96277./8
                      +T*(-99547./24+T*(-28909./36
                          +T*(-10707./160))))))))) * h_1;
        T--;
        Tphi[7] = 68776./7+T*(1038011./28+T*(31157515./504+T*(956669./16
                +T*(3548009./96+T*(2422263./160+T*(197255./48
                      +T*(19959./28+T*(144545./2016+T*(3569./1120)))))))));
        Tdphi[7] = (1038011./28+T*(31157515./252+T*(2870007./16
                +T*(3548009./24+T*(2422263./32+T*(197255./8
                      +T*(19959./4+T*(144545./252
                          +T*(32121./1120))))))))) * h_1;
        T--;
        Tphi[8] = -56375+T*(-8314091./56+T*(-49901303./288+T*(-3763529./32
                +T*(-19648027./384+T*(-9469163./640+T*(-545977./192
                      +T*(-156927./448+T*(-28909./1152
                          +T*(-3569./4480)))))))));
        Tdphi[8] = (-8314091./56+T*(-49901303./144+T*(-11290587./32
                +T*(-19648027./96+T*(-9469163./128+T*(-545977./32
                      +T*(-156927./64+T*(-28909./144
                          +T*(-32121./4480))))))))) * h_1;
        T--;
        Tphi[9] = 439375./7+T*(64188125./504+T*(231125375./2016
              +T*(17306975./288+T*(7761805./384+T*(2895587./640
                    +T*(129391./192+T*(259715./4032+T*(28909./8064
                          +T*(3569./40320)))))))));
        Tdphi[9] = (64188125./504+T*(231125375./1008
              +T*(17306975./96+T*(7761805./96+T*(2895587./128
                    +T*(129391./32+T*(259715./576+T*(28909./1008
                          +T*(3569./4480))))))))) * h_1;
        for (int i=0;  i < 10;  i++) {
          phi[i] = Float(Tphi[i]);
          dphi[i] = Float(Tdphi[i]);
        }
      } // end grouping to define local variables
        break;
      default:
        NAMD_die("Unknown MSM approximation.");
    } // switch
  } // d_stencil_1d()

  void stencil_1d_c1hermite(Float phi[], Float psi[], Float t, Float h) {
    phi[0] = (1 - t) * (1 - t) * (1 + 2*t);
    psi[0] = h * t * (1 - t) * (1 - t);
    t--;
    phi[1] = (1 + t) * (1 + t) * (1 - 2*t);
    psi[1] = h * t * (1 + t) * (1 + t);
  }

  void d_stencil_1d_c1hermite(
      Float dphi[], Float phi[], Float dpsi[], Float psi[],
      Float t, Float h, Float h_1) {
    phi[0] = (1 - t) * (1 - t) * (1 + 2*t);
    dphi[0] = -6 * t * (1 - t) * h_1;
    psi[0] = h * t * (1 - t) * (1 - t);
    dpsi[0] = (1 - t) * (1 - 3*t);
    t--;
    phi[1] = (1 + t) * (1 + t) * (1 - 2*t);
    dphi[1] = -6 * t * (1 + t) * h_1;
    psi[1] = h * t * (1 + t) * (1 + t);
    dpsi[1] = (1 + t) * (1 + 3*t);
  }

  static void ndsplitting(BigReal pg[], BigReal s, int n, int _split) {
    int k = 0;
    if (k == n) return;
    if (s <= 1) {
      // compute derivatives of smoothed part
      switch (_split) {
        case TAYLOR2:
          pg[k++] = 1 + (s-1)*(-1./2 + (s-1)*(3./8));
          if (k == n) break;
          pg[k++] = -1./2 + (s-1)*(3./4);
          if (k == n) break;
          pg[k++] = 3./4;
          break;
        case TAYLOR3:
          pg[k++] = 1 + (s-1)*(-1./2 + (s-1)*(3./8 + (s-1)*(-5./16)));
          if (k == n) break;
          pg[k++] = -1./2 + (s-1)*(3./4 + (s-1)*(-15./16));
          if (k == n) break;
          pg[k++] = 3./4 + (s-1)*(-15./8);
          if (k == n) break;
          pg[k++] = -15./8;
          break;
        case TAYLOR4:
          pg[k++] = 1 + (s-1)*(-1./2 + (s-1)*(3./8 + (s-1)*(-5./16
                  + (s-1)*(35./128))));
          if (k == n) break;
          pg[k++] = -1./2 + (s-1)*(3./4 + (s-1)*(-15./16 + (s-1)*(35./32)));
          if (k == n) break;
          pg[k++] = 3./4 + (s-1)*(-15./8 + (s-1)*(105./32));
          if (k == n) break;
          pg[k++] = -15./8 + (s-1)*(105./16);
          if (k == n) break;
          pg[k++] = 105./16;
          break;
        case TAYLOR5:
          pg[k++] = 1 + (s-1)*(-1./2 + (s-1)*(3./8 + (s-1)*(-5./16
                  + (s-1)*(35./128 + (s-1)*(-63./256)))));
          if (k == n) break;
          pg[k++] = -1./2 + (s-1)*(3./4 + (s-1)*(-15./16 + (s-1)*(35./32
                  + (s-1)*(-315./256))));
          if (k == n) break;
          pg[k++] = 3./4 + (s-1)*(-15./8 + (s-1)*(105./32 + (s-1)*(-315./64)));
          if (k == n) break;
          pg[k++] = -15./8 + (s-1)*(105./16 + (s-1)*(-945./64));
          if (k == n) break;
          pg[k++] = 105./16 + (s-1)*(-945./32);
          if (k == n) break;
          pg[k++] = -945./32;
          break;
        case TAYLOR6:
          pg[k++] = 1 + (s-1)*(-1./2 + (s-1)*(3./8 + (s-1)*(-5./16
                  + (s-1)*(35./128 + (s-1)*(-63./256 + (s-1)*(231./1024))))));
          if (k == n) break;
          pg[k++] = -1./2 + (s-1)*(3./4 + (s-1)*(-15./16 + (s-1)*(35./32
                  + (s-1)*(-315./256 + (s-1)*(693./512)))));
          if (k == n) break;
          pg[k++] = 3./4 + (s-1)*(-15./8 + (s-1)*(105./32 + (s-1)*(-315./64
                  + (s-1)*(3465./512))));
          if (k == n) break;
          pg[k++] = -15./8 + (s-1)*(105./16 + (s-1)*(-945./64
                + (s-1)*(3465./128)));
          if (k == n) break;
          pg[k++] = 105./16 + (s-1)*(-945./32 + (s-1)*(10395./128));
          if (k == n) break;
          pg[k++] = -945./32 + (s-1)*(10395./64);
          if (k == n) break;
          pg[k++] = 10395./64;
          break;
        case TAYLOR7:
          pg[k++] = 1 + (s-1)*(-1./2 + (s-1)*(3./8 + (s-1)*(-5./16
                  + (s-1)*(35./128 + (s-1)*(-63./256
                      + (s-1)*(231./1024 + (s-1)*(-429./2048)))))));
          if (k == n) break;
          pg[k++] = -1./2 + (s-1)*(3./4 + (s-1)*(-15./16 + (s-1)*(35./32
                  + (s-1)*(-315./256 + (s-1)*(693./512
                      + (s-1)*(-3003./2048))))));
          if (k == n) break;
          pg[k++] = 3./4 + (s-1)*(-15./8 + (s-1)*(105./32 + (s-1)*(-315./64
                  + (s-1)*(3465./512 + (s-1)*(-9009./1024)))));
          if (k == n) break;
          pg[k++] = -15./8 + (s-1)*(105./16 + (s-1)*(-945./64 + (s-1)*(3465./128
                  + (s-1)*(-45045./1024))));
          if (k == n) break;
          pg[k++] = 105./16 + (s-1)*(-945./32 + (s-1)*(10395./128
                + (s-1)*(-45045./256)));
          if (k == n) break;
          pg[k++] = -945./32 + (s-1)*(10395./64 + (s-1)*(-135135./256));
          if (k == n) break;
          pg[k++] = 10395./64 + (s-1)*(-135135./128);
          if (k == n) break;
          pg[k++] = -135135./128;
          break;
        case TAYLOR8:
          pg[k++] = 1 + (s-1)*(-1./2 + (s-1)*(3./8 + (s-1)*(-5./16
                  + (s-1)*(35./128 + (s-1)*(-63./256
                      + (s-1)*(231./1024 + (s-1)*(-429./2048
                          + (s-1)*(6435./32768))))))));
          if (k == n) break;
          pg[k++] = -1./2 + (s-1)*(3./4 + (s-1)*(-15./16 + (s-1)*(35./32
                  + (s-1)*(-315./256 + (s-1)*(693./512
                      + (s-1)*(-3003./2048 + (s-1)*(6435./4096)))))));
          if (k == n) break;
          pg[k++] = 3./4 + (s-1)*(-15./8 + (s-1)*(105./32 + (s-1)*(-315./64
                  + (s-1)*(3465./512 + (s-1)*(-9009./1024
                      + (s-1)*(45045./4096))))));
          if (k == n) break;
          pg[k++] = -15./8 + (s-1)*(105./16 + (s-1)*(-945./64 + (s-1)*(3465./128
                  + (s-1)*(-45045./1024 + (s-1)*(135135./2048)))));
          if (k == n) break;
          pg[k++] = 105./16 + (s-1)*(-945./32 + (s-1)*(10395./128
                + (s-1)*(-45045./256 + (s-1)*(675675./2048))));
          if (k == n) break;
          pg[k++] = -945./32 + (s-1)*(10395./64 + (s-1)*(-135135./256
                + (s-1)*(675675./512)));
          if (k == n) break;
          pg[k++] = 10395./64 + (s-1)*(-135135./128 + (s-1)*(2027025./512));
          if (k == n) break;
          pg[k++] = -135135./128 + (s-1)*(2027025./256);
          if (k == n) break;
          pg[k++] = 2027025./256;
          break;
        default:
          NAMD_die("Unknown MSM splitting.");
      }
    } // if (s <= 1)
    else { // (s > 1)
      // compute derivatives of s^(-1/2)
      const BigReal s_1 = 1./s;
      BigReal s_p = sqrt(s_1);
      BigReal p = -0.5;
      BigReal _c = 1;
      pg[k++] = _c * s_p;
      while (k < n) {
        s_p *= s_1;
        _c *= p;
        p -= 1;
        pg[k++] = _c * s_p;
      }
    } // else (s > 1)
    // higher derivatives are zero
    while (k < n) pg[k++] = 0;
  } // ndsplitting()


  static void gc_c1hermite_elem_accum(C1Matrix& matrix, BigReal _c,
      Vector rv, BigReal _a, int _split) {
    const BigReal a_1 = 1./_a;
    const BigReal a_2 = a_1 * a_1;
    const BigReal s = (rv * rv) * a_2;
    const BigReal dx = -2 * rv.x * a_2;  // ds/dx
    const BigReal dy = -2 * rv.y * a_2;  // ds/dy
    const BigReal dz = -2 * rv.z * a_2;  // ds/dz
    const BigReal dd = 2 * a_2;  // d^2s/dx^2 = d^2s/dy^2 = d^2s/dz^2
    BigReal tmp;
    enum { nderiv = C1_VECTOR_SIZE-1 };
    BigReal p[nderiv];
    Float *g = matrix.melem;

    // multiply entire matrix by this coefficient
    _c = _c * a_1;

    // compute derivatives (d/ds)^k of splitting g(s), s=r^2
    ndsplitting(p, s, nderiv, _split);

    // weight 0
    tmp = _c * p[0];
    g[C1INDEX(D000,D000)] += tmp;

    // weight 1
    tmp = _c * p[1] * dx;
    g[C1INDEX(D100,D000)] += tmp;
    g[C1INDEX(D000,D100)] -= tmp;

    tmp = _c * p[1] * dy;
    g[C1INDEX(D010,D000)] += tmp;
    g[C1INDEX(D000,D010)] -= tmp;

    tmp = _c * p[1] * dz;
    g[C1INDEX(D001,D000)] += tmp;
    g[C1INDEX(D000,D001)] -= tmp;

    // C1 splitting returns here

    // weight 2
    tmp = _c * p[2] * dx * dy;
    g[C1INDEX(D110,D000)] += tmp;
    g[C1INDEX(D000,D110)] += tmp;
    g[C1INDEX(D100,D010)] -= tmp;
    g[C1INDEX(D010,D100)] -= tmp;

    tmp = _c * p[2] * dx * dz;
    g[C1INDEX(D101,D000)] += tmp;
    g[C1INDEX(D000,D101)] += tmp;
    g[C1INDEX(D100,D001)] -= tmp;
    g[C1INDEX(D001,D100)] -= tmp;

    tmp = _c * p[2] * dy * dz;
    g[C1INDEX(D011,D000)] += tmp;
    g[C1INDEX(D000,D011)] += tmp;
    g[C1INDEX(D010,D001)] -= tmp;
    g[C1INDEX(D001,D010)] -= tmp;

    tmp = _c * (p[2] * dx*dx + p[1] * dd);
    g[C1INDEX(D100,D100)] -= tmp;
    tmp = _c * (p[2] * dy*dy + p[1] * dd);
    g[C1INDEX(D010,D010)] -= tmp;
    tmp = _c * (p[2] * dz*dz + p[1] * dd);
    g[C1INDEX(D001,D001)] -= tmp;

    // C2 splitting returns here
    if (_split == TAYLOR2) return;

    // weight 3
    tmp = _c * p[3] * dx * dy * dz;
    g[C1INDEX(D111,D000)] += tmp;
    g[C1INDEX(D110,D001)] -= tmp;
    g[C1INDEX(D101,D010)] -= tmp;
    g[C1INDEX(D011,D100)] -= tmp;
    g[C1INDEX(D100,D011)] += tmp;
    g[C1INDEX(D010,D101)] += tmp;
    g[C1INDEX(D001,D110)] += tmp;
    g[C1INDEX(D000,D111)] -= tmp;

    tmp = _c * (p[3] * dx*dx * dy + p[2] * dd * dy);
    g[C1INDEX(D110,D100)] -= tmp;
    g[C1INDEX(D100,D110)] += tmp;

    tmp = _c * (p[3] * dx*dx * dz + p[2] * dd * dz);
    g[C1INDEX(D101,D100)] -= tmp;
    g[C1INDEX(D100,D101)] += tmp;

    tmp = _c * (p[3] * dy*dy * dx + p[2] * dd * dx);
    g[C1INDEX(D110,D010)] -= tmp;
    g[C1INDEX(D010,D110)] += tmp;

    tmp = _c * (p[3] * dy*dy * dz + p[2] * dd * dz);
    g[C1INDEX(D011,D010)] -= tmp;
    g[C1INDEX(D010,D011)] += tmp;

    tmp = _c * (p[3] * dz*dz * dx + p[2] * dd * dx);
    g[C1INDEX(D101,D001)] -= tmp;
    g[C1INDEX(D001,D101)] += tmp;

    tmp = _c * (p[3] * dz*dz * dy + p[2] * dd * dy);
    g[C1INDEX(D011,D001)] -= tmp;
    g[C1INDEX(D001,D011)] += tmp;

    // C3 splitting returns here
    if (_split == TAYLOR3) return;

    // weight 4
    tmp = _c * (p[4] * dx*dx * dy * dz + p[3] * dd * dy * dz);
    g[C1INDEX(D111,D100)] -= tmp;
    g[C1INDEX(D100,D111)] -= tmp;
    g[C1INDEX(D110,D101)] += tmp;
    g[C1INDEX(D101,D110)] += tmp;

    tmp = _c * (p[4] * dy*dy * dx * dz + p[3] * dd * dx * dz);
    g[C1INDEX(D111,D010)] -= tmp;
    g[C1INDEX(D010,D111)] -= tmp;
    g[C1INDEX(D110,D011)] += tmp;
    g[C1INDEX(D011,D110)] += tmp;

    tmp = _c * (p[4] * dz*dz * dx * dy + p[3] * dd * dx * dy);
    g[C1INDEX(D111,D001)] -= tmp;
    g[C1INDEX(D001,D111)] -= tmp;
    g[C1INDEX(D101,D011)] += tmp;
    g[C1INDEX(D011,D101)] += tmp;

    tmp = _c * (p[4] * dx*dx * dy*dy + p[3] * dx*dx * dd
        + p[3] * dd * dy*dy + p[2] * dd * dd);
    g[C1INDEX(D110,D110)] += tmp;
    tmp = _c * (p[4] * dx*dx * dz*dz + p[3] * dx*dx * dd
        + p[3] * dd * dz*dz + p[2] * dd * dd);
    g[C1INDEX(D101,D101)] += tmp;
    tmp = _c * (p[4] * dy*dy * dz*dz + p[3] * dy*dy * dd
        + p[3] * dd * dz*dz + p[2] * dd * dd);
    g[C1INDEX(D011,D011)] += tmp;

    // C4 splitting returns here
    if (_split == TAYLOR4) return;

    // weight 5
    tmp = _c * (p[5] * dx*dx * dy*dy * dz + p[4] * dx*dx * dd * dz
        + p[4] * dd * dy*dy * dz + p[3] * dd * dd * dz);
    g[C1INDEX(D111,D110)] += tmp;
    g[C1INDEX(D110,D111)] -= tmp;

    tmp = _c * (p[5] * dx*dx * dz*dz * dy + p[4] * dx*dx * dd * dy
        + p[4] * dd * dz*dz * dy + p[3] * dd * dd * dy);
    g[C1INDEX(D111,D101)] += tmp;
    g[C1INDEX(D101,D111)] -= tmp;

    tmp = _c * (p[5] * dy*dy * dz*dz * dx + p[4] * dy*dy * dd * dx
        + p[4] * dd * dz*dz * dx + p[3] * dd * dd * dx);
    g[C1INDEX(D111,D011)] += tmp;
    g[C1INDEX(D011,D111)] -= tmp;

    // C5 splitting returns here
    if (_split == TAYLOR5) return;

    // weight 6
    tmp = _c * (p[6] * dx*dx * dy*dy * dz*dz + p[5] * dx*dx * dy*dy * dd
        + p[5] * dx*dx * dd * dz*dz + p[5] * dd * dy*dy * dz*dz
        + p[4] * dx*dx * dd * dd + p[4] * dd * dy*dy * dd
        + p[4] * dd * dd * dz*dz + p[3] * dd * dd * dd);
    g[C1INDEX(D111,D111)] -= tmp;

    // calculate full matrix for C6 or higher splitting

  } // gc_c1hermite_elem_accum()


}; // ComputeMsmMgr


// Degree of polynomial basis function Phi.
// For the purpose of finding the stencil width, Hermite interpolation 
// sets this value to 1.
const int ComputeMsmMgr::PolyDegree[NUM_APPROX] = {
  3, 5, 5, 7, 7, 9, 9, 1,
};

// The stencil array lengths below.
const int ComputeMsmMgr::Nstencil[NUM_APPROX] = {
  5, 7, 7, 9, 9, 11, 11, 3,
};

// Index offsets from the stencil-centered grid element, to get
// to the correct contributing grid element.
const int
ComputeMsmMgr::IndexOffset[NUM_APPROX][MAX_NSTENCIL_SKIP_ZERO] = {
  // cubic
  {-3, -1, 0, 1, 3},

  // quintic C1
  {-5, -3, -1, 0, 1, 3, 5},

  // quintic C2  (same as quintic C1)
  {-5, -3, -1, 0, 1, 3, 5},

  // septic C1
  {-7, -5, -3, -1, 0, 1, 3, 5, 7},

  // septic C3  (same as septic C1)
  {-7, -5, -3, -1, 0, 1, 3, 5, 7},

  // nonic C1
  {-9, -7, -5, -3, -1, 0, 1, 3, 5, 7, 9},

  // nonic C4  (same as nonic C1)
  {-9, -7, -5, -3, -1, 0, 1, 3, 5, 7, 9},

  // C1 Hermite
  {-1, 0, 1},
};

// The grid transfer stencils for the non-factored restriction and
// prolongation procedures.
const Float
ComputeMsmMgr::PhiStencil[NUM_APPROX_FORMS][MAX_NSTENCIL_SKIP_ZERO] = {
  // cubic
  {-1.f/16, 9.f/16, 1, 9.f/16, -1.f/16},

  // quintic C1
  {3.f/256, -25.f/256, 75.f/128, 1, 75.f/128, -25.f/256, 3.f/256},

  // quintic C2  (same as quintic C1)
  {3.f/256, -25.f/256, 75.f/128, 1, 75.f/128, -25.f/256, 3.f/256},

  // septic C1
  { -5.f/2048, 49.f/2048, -245.f/2048, 1225.f/2048, 1, 1225.f/2048,
    -245.f/2048, 49.f/2048, -5.f/2048 },

  // septic C3  (same as septic C3)
  { -5.f/2048, 49.f/2048, -245.f/2048, 1225.f/2048, 1, 1225.f/2048,
    -245.f/2048, 49.f/2048, -5.f/2048 },

  // nonic C1
  { 35.f/65536, -405.f/65536, 567.f/16384, -2205.f/16384, 
    19845.f/32768, 1, 19845.f/32768, -2205.f/16384, 567.f/16384, 
    -405.f/65536, 35.f/65536 },

  // nonic C4  (same as nonic C1)
  { 35.f/65536, -405.f/65536, 567.f/16384, -2205.f/16384, 
    19845.f/32768, 1, 19845.f/32768, -2205.f/16384, 567.f/16384, 
    -405.f/65536, 35.f/65536 },
};


// Designates PE assignment for static load balancing of 
// MsmBlock-related arrays
class MsmBlockMap : public CkArrayMap {
  private:
    ComputeMsmMgr *mgrLocal;
    int *penum;
    int level;
  public:
    MsmBlockMap(int lvl) {
      mgrLocal = CProxy_ComputeMsmMgr::ckLocalBranch(
          CkpvAccess(BOCclass_group).computeMsmMgr);
#ifdef MSM_NODE_MAPPING
      penum = mgrLocal->blockAssign.buffer();
#else
      penum = 0;
#endif
      level = lvl;
    }
    MsmBlockMap(CkMigrateMessage *m) { }
    int registerArray(CkArrayIndex& numElements, CkArrayID aid) {
      return 0;
    }
    int procNum(int /*arrayHdl*/, const CkArrayIndex &idx) {
      int *pn = (int *)idx.data();
#ifdef MSM_NODE_MAPPING
      int n = mgrLocal->blockFlatIndex(level, pn[0], pn[1], pn[2]);
      return penum[n];
#else
      return 0;
#endif
    }
};


// Designates PE assignment for static load balancing of 
// MsmGridCutoff-related arrays
class MsmGridCutoffMap : public CkArrayMap {
  private:
    int *penum;
  public:
    MsmGridCutoffMap() {
      ComputeMsmMgr *mgrLocal = CProxy_ComputeMsmMgr::ckLocalBranch(
          CkpvAccess(BOCclass_group).computeMsmMgr);
#ifdef MSM_NODE_MAPPING
      penum = mgrLocal->gcutAssign.buffer();
#else
      penum = 0;
#endif
    }
    int registerArray(CkArrayIndex& numElements, CkArrayID aid) {
      return 0;
    }
    int procNum(int /*arrayHdl*/, const CkArrayIndex &idx) {
#if 1
      int n = *((int *)idx.data());
#ifdef MSM_NODE_MAPPING
      return penum[n];
#else
      return 0;
#endif
#else
      return 0;  // XXX to test load balancing
#endif
    }
};


namespace msm {

  //
  // PatchData
  //
  // Performs anterpolation and interpolation algorithms.
  //
  // Surround each NAMD patch with enough grid points to perform 
  // anterpolation and interpolation without having to do any 
  // grid wrapping.  This does not give a partitioning of the 
  // MSM finest level grid --- rather, the edges of adjacent 
  // PatchData grids will overlap or contain image points along 
  // the periodic boundaries.  
  //
 
  struct PatchData {
    ComputeMsmMgr *mgr;
    Map *map;
    PatchDiagram *pd;
    AtomCoordArray coord;
    ForceArray force;
    Grid<Float> qh;
    Grid<Float> eh;
    Grid<Float> subgrid;
    Grid<C1Vector> qh_c1hermite;
    Grid<C1Vector> eh_c1hermite;
    Grid<C1Vector> subgrid_c1hermite;
    BigReal energy;
    //BigReal virial[3][3];
    int cntRecvs;
    int patchID;
    int sequence;  // from Compute object for message priority

    AtomCoordArray& coordArray() { return coord; }
    ForceArray& forceArray() { return force; }

    PatchData(ComputeMsmMgr *pmgr, int pid);
    void init(int natoms);

    void anterpolation();
    void sendCharge();
    void addPotential(const Grid<Float>& epart);
    void interpolation();

    void anterpolationC1Hermite();
    void sendChargeC1Hermite();
    void addPotentialC1Hermite(const Grid<C1Vector>& epart);
    void interpolationC1Hermite();
  };

} // namespace msm


/////////////////
//
// MsmGridCutoff
//
// Performs grid cutoff part of the computation.
//
// The grid cutoff part is the most computationally intensive part 
// of MSM.  The templated MsmGridCutoffKernel class takes Vtype 
// for charge and potential data (generalizes to vector for Hermite
// interpolation) and takes Mtype for the pre-computed grid coefficient 
// weights (generalizes to matrix for Hermite interpolation).
//

template <class Vtype, class Mtype>
class MsmGridCutoffKernel {
  public:
    ComputeMsmMgr *mgrLocal;     // for quick access to data
    msm::Map *map;
    msm::BlockIndex qhblockIndex;  // source of charges
    msm::BlockSend ehblockSend;    // destination for potentials
    int eia, eib, eja, ejb, eka, ekb, eni, enj, enk;  // for "fold factor"
    int isfold;  // for "fold factor"
    msm::Grid<Vtype> qh;
    msm::Grid<Vtype> eh;
    msm::Grid<Vtype> ehfold;  // for "fold factor"
    const msm::Grid<Mtype> *pgc;
    const msm::Grid<Mtype> *pgvc;
    int priority;
    int sequence;

    MsmGridCutoffKernel() { init(); }

    void init() {
      isfold = 0;
      mgrLocal = CProxy_ComputeMsmMgr::ckLocalBranch(
          CkpvAccess(BOCclass_group).computeMsmMgr);
      map = &(mgrLocal->mapData());
      mgrLocal->addVirialContrib();
#ifdef MSM_TIMING
      mgrLocal->addTiming();
#endif
#ifdef MSM_PROFILING
      mgrLocal->addProfiling();
#endif
    }

#ifdef MSM_MIGRATION
    void pup(PUP::er& p) {
#ifdef MSM_TIMING
      mgrLocal->subtractTiming();
#endif
#ifdef MSM_PROFILING
      mgrLocal->subtractProfiling();
#endif
      p | qhblockIndex;
      p | ehblockSend;
      p | eia, p | eib, p | eja, p | ejb, p | eka, p | ekb;
      p | eni, p | enj, p | enk;
      p | isfold;
    }
#endif // MSM_MIGRATION

    void setup(MsmGridCutoffInitMsg *bmsg) {
      qhblockIndex = bmsg->qhBlockIndex;
      ehblockSend = bmsg->ehBlockSend;
      delete bmsg;

      // set message priority
      priority = mgrLocal->nlevels
        + 2*(mgrLocal->nlevels - ehblockSend.nblock_wrap.level) - 1;
      // allocate qh buffer
      qh.init(map->blockLevel[qhblockIndex.level](qhblockIndex.n).nrange);
      // allocate eh buffer
      eh.init(ehblockSend.nrange);
      // preprocess "fold factor" if active for this level
      if (map->foldfactor[qhblockIndex.level].active) {
        // allocate ehfold buffer
        ehfold = eh;
        // set index range of potentials
        eia = eh.ia();
        eib = eh.ib();
        eja = eh.ja();
        ejb = eh.jb();
        eka = eh.ka();
        ekb = eh.kb();
        eni = eh.ni();
        enj = eh.nj();
        enk = eh.nk();
        if (map->blockLevel[qhblockIndex.level].nn() == 1) {
          if (map->ispx) { eia = qh.ia();  eib = qh.ib();  eni = qh.ni(); }
          if (map->ispy) { eja = qh.ja();  ejb = qh.jb();  enj = qh.nj(); }
          if (map->ispz) { eka = qh.ka();  ekb = qh.kb();  enk = qh.nk(); }
        }
        else {
          // find destination block index
          int level = qhblockIndex.level;
          msm::BlockIndex bn = map->blockOfGridIndex(
              ehblockSend.nrange_wrap.lower(), level);
          map->wrapBlockIndex(bn);
          if (map->ispx) {
            eia = bn.n.i * map->bsx[level];
            eib = eia + qh.ni() - 1;
            eni = qh.ni();
          }
          if (map->ispy) {
            eja = bn.n.j * map->bsy[level];
            ejb = eja + qh.nj() - 1;
            enj = qh.nj();
          }
          if (map->ispz) {
            eka = bn.n.k * map->bsz[level];
            ekb = eka + qh.nk() - 1;
            enk = qh.nk();
          }
        }
        isfold = 1;
      } // if fold factor
    } // setup()

    void setupWeights(
        const msm::Grid<Mtype> *ptrgc,
        const msm::Grid<Mtype> *ptrgvc
        ) {
      pgc = ptrgc;
      pgvc = ptrgvc;
    } // setupWeights()


    void compute(GridMsg *gmsg) {
#ifdef MSM_TIMING
      double startTime, stopTime;
      startTime = CkWallTimer();
#endif
      //
      // receive block of charges
      //
      int pid;
      // qh is resized only the first time, memory allocation persists
      gmsg->get(qh, pid, sequence);
      delete gmsg;
#ifdef MSM_TIMING
      stopTime = CkWallTimer();
      mgrLocal->msmTiming[MsmTimer::COMM] += stopTime - startTime;
#endif

      //
      // grid cutoff calculation
      // this charge block -> this potential block
      //

#ifdef MSM_TIMING
      startTime = stopTime;
#endif
      // resets indexing on block
      eh.init(ehblockSend.nrange);  // (always have to re-init nrange for eh)
      eh.reset(0);
      // index range of weights
      int gia = pgc->ia();
      int gib = pgc->ib();
      int gja = pgc->ja();
      int gjb = pgc->jb();
      int gka = pgc->ka();
      int gkb = pgc->kb();
      int gni = pgc->ni();
      int gnj = pgc->nj();
      // index range of charge grid
      int qia = qh.ia();
      int qib = qh.ib();
      int qja = qh.ja();
      int qjb = qh.jb();
      int qka = qh.ka();
      int qkb = qh.kb();
      int qni = qh.ni();
      int qnj = qh.nj();
      // index range of potentials
      int ia = eh.ia();
      int ib = eh.ib();
      int ja = eh.ja();
      int jb = eh.jb();
      int ka = eh.ka();
      int kb = eh.kb();

      int index = 0;

      // access buffers directly
      const Mtype *gcbuffer = pgc->data().buffer();
      //const Mtype *gvcbuffer = pgvc->data().buffer();
      const Vtype *qhbuffer = qh.data().buffer();
      Vtype *ehbuffer = eh.data().buffer();
      //Vtype *gvsumbuffer = mgrLocal->gvsum.data().buffer();

#ifndef MSM_COMM_ONLY
      // loop over potentials
      for (int k = ka;  k <= kb;  k++) {
        // clip charges to weights along k
        int mka = ( qka >= gka + k ? qka : gka + k );
        int mkb = ( qkb <= gkb + k ? qkb : gkb + k );

        for (int j = ja;  j <= jb;  j++) {
          // clip charges to weights along j
          int mja = ( qja >= gja + j ? qja : gja + j );
          int mjb = ( qjb <= gjb + j ? qjb : gjb + j );

          for (int i = ia;  i <= ib;  i++) {
            // clip charges to weights along i
            int mia = ( qia >= gia + i ? qia : gia + i );
            int mib = ( qib <= gib + i ? qib : gib + i );

            // accumulate sum to this eh point
            Vtype ehsum = 0;

#if 0
            // loop over charge grid
            for (int qk = mka;  qk <= mkb;  qk++) {
              int qkoff = (qk - qka) * qnj;
              int gkoff = ((qk-k) - gka) * gnj;

              for (int qj = mja;  qj <= mjb;  qj++) {
                int qjkoff = (qkoff + qj - qja) * qni;
                int gjkoff = (gkoff + (qj-j) - gja) * gni;

// help the vectorizer make reasonable decisions
#if defined(__INTEL_COMPILER)
#pragma vector always 
#endif
                for (int qi = mia;  qi <= mib;  qi++) {
                  int qijkoff = qjkoff + qi - qia;
                  int gijkoff = gjkoff + (qi-i) - gia;

                  ehsum += gcbuffer[gijkoff] * qhbuffer[qijkoff];
                }
              }
            } // end loop over charge grid
#else

#if 0
            // loop over charge grid
            int nn = mib - mia + 1;
            for (int qk = mka;  qk <= mkb;  qk++) {
              int qkoff = (qk - qka) * qnj;
              int gkoff = ((qk-k) - gka) * gnj;

              for (int qj = mja;  qj <= mjb;  qj++) {
                int qjkoff = (qkoff + qj - qja) * qni;
                int gjkoff = (gkoff + (qj-j) - gja) * gni;

                const Float *qbuf = qhbuffer + (qjkoff - qia + mia);
                const Float *gbuf = gcbuffer + (gjkoff - i - gia + mia);
#ifdef MSM_PROFILING
                mgrLocal->xLoopCnt[nn]++;
#endif
// help the vectorizer make reasonable decisions
#if defined(__INTEL_COMPILER)
#pragma vector always 
#endif
                for (int ii = 0;  ii < nn;  ii++) {
                  ehsum += gbuf[ii] * qbuf[ii];
                }
              }
            } // end loop over charge grid
#else
            // loop over charge grid
            int nn = mib - mia + 1;
            if (nn == 8) {  // hard coded inner loop = 8
              int qnji = qnj * qni;
              int qkoff = -qka*qnji - qja*qni - qia + mia;
              int gnji = gnj * gni;
              int gkoff = (-k-gka)*gnji + (-j-gja)*gni - i - gia + mia;

              for (int qk = mka;  qk <= mkb;  qk++) {
                int qjkoff = qkoff + qk*qnji;
                int gjkoff = gkoff + qk*gnji;

                for (int qj = mja;  qj <= mjb;  qj++) {
                  const Vtype *qbuf = qhbuffer + (qjkoff + qj*qni);
                  const Mtype *gbuf = gcbuffer + (gjkoff + qj*gni);
                  //const Mtype *gvcbuf = gvcbuffer + (gjkoff + qj*gni);
                  //Vtype *gvsumbuf = gvsumbuffer + (gjkoff + qj*gni);
#ifdef MSM_PROFILING
                  mgrLocal->xLoopCnt[nn]++;
#endif
// help the vectorizer make reasonable decisions
#if defined(__INTEL_COMPILER)
#pragma vector always 
#endif
                  for (int ii = 0;  ii < 8;  ii++) {
                    ehsum += gbuf[ii] * qbuf[ii];
                    //gvsumbuf[ii] += qbuf[ii] * qbuf[ii] * gvcbuf[ii];
                  }
                }
              } // end loop over charge grid
            }
            else {  // variable length inner loop < 8
              int qnji = qnj * qni;
              int qkoff = -qka*qnji - qja*qni - qia + mia;
              int gnji = gnj * gni;
              int gkoff = (-k-gka)*gnji + (-j-gja)*gni - i - gia + mia;

              for (int qk = mka;  qk <= mkb;  qk++) {
                int qjkoff = qkoff + qk*qnji;
                int gjkoff = gkoff + qk*gnji;

                for (int qj = mja;  qj <= mjb;  qj++) {
                  const Vtype *qbuf = qhbuffer + (qjkoff + qj*qni);
                  const Mtype *gbuf = gcbuffer + (gjkoff + qj*gni);
                  //const Mtype *gvcbuf = gvcbuffer + (gjkoff + qj*gni);
                  //Vtype *gvsumbuf = gvsumbuffer + (gjkoff + qj*gni);
#ifdef MSM_PROFILING
                  mgrLocal->xLoopCnt[nn]++;
#endif
// help the vectorizer make reasonable decisions
#if defined(__INTEL_COMPILER)
#pragma vector always 
#endif
                  for (int ii = 0;  ii < nn;  ii++) {
                    ehsum += gbuf[ii] * qbuf[ii];
                    //gvsumbuf[ii] += qbuf[ii] * qbuf[ii] * gvcbuf[ii];
                  }
                }
              } // end loop over charge grid
            }
#endif // 0

#endif // 0

            ehbuffer[index] = ehsum;
            index++;
          }
        }
      } // end loop over potentials
#endif // !MSM_COMM_ONLY

#ifdef MSM_PROFILING
      mgrLocal->doneProfiling();
#endif

      //
      // send block of potentials
      //

#ifdef MSM_FOLD_FACTOR
      // if "fold factor" is active for this level,
      // need to sum unfolded potential grid back into periodic grid
      if (isfold) {
        // copy unfolded grid
        ehfold = eh;
        // reset eh indexing to correctly folded size
        eh.set(eia, eni, eja, enj, eka, enk);
        eh.reset(0);
#ifdef DEBUG_MSM_GRID
        printf("level=%d   ehfold:  [%d..%d] x [%d..%d] x [%d..%d]  "
            "(%d x %d x %d)\n"
                "              eh:  [%d..%d] x [%d..%d] x [%d..%d]  "
            "(%d x %d x %d)\n"
               "         eh lower:  %d %d %d\n",
            qhblockIndex.level,
            ehfold.ia(), ehfold.ib(), 
            ehfold.ja(), ehfold.jb(),
            ehfold.ka(), ehfold.kb(),
            ehfold.ni(), ehfold.nj(), ehfold.nk(),
            eh.ia(), eh.ib(), 
            eh.ja(), eh.jb(),
            eh.ka(), eh.kb(),
            eh.ni(), eh.nj(), eh.nk(),
            ehblockSend.nrange_wrap.lower().i,
            ehblockSend.nrange_wrap.lower().j,
            ehblockSend.nrange_wrap.lower().k
            );
#endif
        const Vtype *ehfoldbuf = ehfold.data().buffer();
        Vtype *ehbuf = eh.data().buffer();
        // now we "fold" eh by calculating the
        // wrap around sum of ehfold into correctly sized eh
        int index = 0;
        for (int k = ka;  k <= kb;  k++) {
          int kk = k;
          if      (kk < eka)  do { kk += enk; } while (kk < eka);
          else if (kk > ekb)  do { kk -= enk; } while (kk > ekb);
          int koff = (kk - eka) * enj;
          for (int j = ja;  j <= jb;  j++) {
            int jj = j;
            if      (jj < eja)  do { jj += enj; } while (jj < eja);
            else if (jj > ejb)  do { jj -= enj; } while (jj > ejb);
            int jkoff = (koff + (jj - eja)) * eni;
            for (int i = ia;  i <= ib;  i++, index++) {
              int ii = i;
              if      (ii < eia)  do { ii += eni; } while (ii < eia);
              else if (ii > eib)  do { ii -= eni; } while (ii > eib);
              int ijkoff = jkoff + (ii - eia);
              ehbuf[ijkoff] += ehfoldbuf[index];
            }
          }
        }
      }
      else {
        // shift grid index range to its true (wrapped) values
        eh.updateLower( ehblockSend.nrange_wrap.lower() );
      }
#else    // !MSM_FOLD_FACTOR
      // shift grid index range to its true (wrapped) values
      eh.updateLower( ehblockSend.nrange_wrap.lower() );
#endif   // MSM_FOLD_FACTOR

#ifdef MSM_TIMING
      stopTime = CkWallTimer();
      mgrLocal->msmTiming[MsmTimer::GRIDCUTOFF] += stopTime - startTime;
#endif
    } // compute()

};


//
// MsmGridCutoff wraps kernel template for approximations 
// that involve only function values (e.g., CUBIC, QUINTIC).
// Elements of 1D chare array.
//
class MsmGridCutoff :
  public CBase_MsmGridCutoff,
  public MsmGridCutoffKernel<Float,Float>
{
  public:
    CProxyElement_MsmBlock msmBlockElementProxy;  // root of reduction
    CkSectionInfo cookie;  // need to save cookie for section reduction
#ifdef MSM_REDUCE_GRID
    msm::Grid<Float> ehfull;
#endif // MSM_REDUCE_GRID

    MsmGridCutoff() { }

    MsmGridCutoff(CkMigrateMessage *m)
#if  ! defined(MSM_MIGRATION)
    { }
#else // MSM_MIGRATION
      : CBase_MsmGridCutoff(m) {
#ifdef DEBUG_MSM_MIGRATE
      printf("MsmGridCutoff element %d migrated to processor %d\n",
          thisIndex, CkMyPe());
#endif
      init();
      // access type dependent constants from map
      MsmGridCutoffKernel<Float,Float>::setupWeights(
          &(map->gc[ehblockSend.nblock_wrap.level]),
          &(map->gvc[ehblockSend.nblock_wrap.level])
          );
    }

    virtual void pup(PUP::er& p) {
#ifdef DEBUG_MSM_MIGRATE
      printf("MsmGridCutoff element %d pupped on processor %d\n",
          thisIndex, CkMyPe());
#endif
      CBase_MsmGridCutoff::pup(p);  // pack our superclass
      MsmGridCutoffKernel<Float,Float>::pup(p);
    }
#endif // MSM_MIGRATION

    void init() {
      MsmGridCutoffKernel<Float,Float>::init();
    }

    void setup(MsmGridCutoffInitMsg *bmsg) {
      // base class consumes this init proxy  message
      MsmGridCutoffKernel<Float,Float>::setup(bmsg);
      // access type dependent constants from map
      MsmGridCutoffKernel<Float,Float>::setupWeights(
          &(map->gc[ehblockSend.nblock_wrap.level]),
          &(map->gvc[ehblockSend.nblock_wrap.level])
          );
#ifdef MSM_REDUCE_GRID
      // allocate full buffer space needed for section reduction
      int level = ehblockSend.nblock_wrap.level;
      int i = ehblockSend.nblock_wrap.n.i;
      int j = ehblockSend.nblock_wrap.n.j;
      int k = ehblockSend.nblock_wrap.n.k;
      ehfull.init( map->blockLevel[level](i,j,k).nrange );
#endif // MSM_REDUCE_GRID
#ifdef DEBUG_MSM_GRID
      printf("MsmGridCutoff[%d]:  setup()"
          " send to level=%d block=(%d,%d,%d)\n",
          thisIndex, ehblockSend.nblock_wrap.level,
          ehblockSend.nblock_wrap.n.i,
          ehblockSend.nblock_wrap.n.j,
          ehblockSend.nblock_wrap.n.k);
#endif
    }

    void setupSections(MsmGridCutoffSetupMsg *msg) {
#ifdef DEBUG_MSM_GRID
      CkPrintf("MSM GRID CUTOFF %d setup section on PE %d\n",
          thisIndex, CkMyPe());
#endif
      CkGetSectionInfo(cookie, msg);  // init the cookie
      msg->get(&msmBlockElementProxy);  // get proxy to MsmBlock
      delete msg;
    }

    void compute(GridMsg *gmsg) {
#ifdef DEBUG_MSM_GRID
      printf("MsmGridCutoff %d:  compute()\n", thisIndex);
#endif
      // base class consumes this grid message
      MsmGridCutoffKernel<Float,Float>::compute(gmsg);

#ifdef MSM_TIMING
      double startTime, stopTime;
      startTime = CkWallTimer();
#endif
#ifdef MSM_REDUCE_GRID

      // perform section reduction over potential grids
      CProxy_CkMulticastMgr mcastProxy =
        CkpvAccess(BOCclass_group).multicastMgr;
      CkMulticastMgr *mcastPtr =
        CProxy_CkMulticastMgr(mcastProxy).ckLocalBranch();
      CkCallback cb(CkIndex_MsmBlock::sumReducedPotential(NULL),
          msmBlockElementProxy);
      // sum into "full" sized buffer needed for contribute
      ehfull.reset(0);
      ehfull += eh;
      mcastPtr->contribute(
          ehfull.nn() * sizeof(Float), ehfull.data().buffer(), 
          CkReduction::sum_float, cookie, cb);

#else
      // place eh into message
      const msm::BlockIndex& bindex = ehblockSend.nblock_wrap;
      int msgsz = eh.data().len() * sizeof(Float);
      GridMsg *gm = new(msgsz, sizeof(int)) GridMsg;
      SET_PRIORITY(gm, sequence, priority);
      gm->put(eh, bindex.level, sequence);
      // lookup in ComputeMsmMgr proxy array by level
      mgrLocal->msmBlock[bindex.level](
          bindex.n.i, bindex.n.j, bindex.n.k).addPotential(gm);

#endif // MSM_REDUCE_GRID

#ifdef MSM_TIMING
      stopTime = CkWallTimer();
      mgrLocal->msmTiming[MsmTimer::COMM] += stopTime - startTime;
      mgrLocal->doneTiming();
#endif
    } // compute()

}; // MsmGridCutoff


//
// MsmC1HermiteGridCutoff wraps kernel template for
// C1 Hermite approximation.  Elements of 1D chare array.
//
class MsmC1HermiteGridCutoff :
  public CBase_MsmC1HermiteGridCutoff,
  public MsmGridCutoffKernel<C1Vector,C1Matrix>
{
  public:
    CProxyElement_MsmC1HermiteBlock msmBlockElementProxy;  // root of reduction
    CkSectionInfo cookie;  // need to save cookie for section reduction
#ifdef MSM_REDUCE_GRID
    msm::Grid<C1Vector> ehfull;
#endif // MSM_REDUCE_GRID

    MsmC1HermiteGridCutoff() { }

    MsmC1HermiteGridCutoff(CkMigrateMessage *m)
#if  ! defined(MSM_MIGRATION)
    { }
#else // MSM_MIGRATION
      : CBase_MsmC1HermiteGridCutoff(m) {
#ifdef DEBUG_MSM_MIGRATE
      printf("MsmC1HermiteGridCutoff element %d migrated to processor %d\n",
          thisIndex, CkMyPe());
#endif
      init();
      // access type dependent constants from map
      MsmGridCutoffKernel<C1Vector,C1Matrix>::setupWeights(
          &(map->gc_c1hermite[ehblockSend.nblock_wrap.level]),
          NULL
          );
    }

    virtual void pup(PUP::er& p) {
#ifdef DEBUG_MSM_MIGRATE
      printf("MsmC1HermiteGridCutoff element %d pupped on processor %d\n",
          thisIndex, CkMyPe());
#endif
      CBase_MsmC1HermiteGridCutoff::pup(p);  // pack our superclass
      MsmGridCutoffKernel<C1Vector,C1Matrix>::pup(p);
    }
#endif // MSM_MIGRATION

    void init() {
      MsmGridCutoffKernel<C1Vector,C1Matrix>::init();
    }

    void setup(MsmGridCutoffInitMsg *bmsg) {
      // base class consumes this init proxy  message
      MsmGridCutoffKernel<C1Vector,C1Matrix>::setup(bmsg);
      // access type dependent constants from map
      MsmGridCutoffKernel<C1Vector,C1Matrix>::setupWeights(
          &(map->gc_c1hermite[ehblockSend.nblock_wrap.level]),
          NULL
          );
#ifdef DEBUG_MSM_GRID
      printf("MsmC1HermiteGridCutoff[%d]:  setup()"
          " send to level=%d block=(%d,%d,%d)\n",
          thisIndex, ehblockSend.nblock_wrap.level,
          ehblockSend.nblock_wrap.n.i,
          ehblockSend.nblock_wrap.n.j,
          ehblockSend.nblock_wrap.n.k);
#endif
#ifdef MSM_REDUCE_GRID
      // allocate full buffer space needed for section reduction
      int level = ehblockSend.nblock_wrap.level;
      int i = ehblockSend.nblock_wrap.n.i;
      int j = ehblockSend.nblock_wrap.n.j;
      int k = ehblockSend.nblock_wrap.n.k;
      ehfull.init( map->blockLevel[level](i,j,k).nrange );
#endif // MSM_REDUCE_GRID
    }

    void setupSections(MsmC1HermiteGridCutoffSetupMsg *msg) {
#ifdef DEBUG_MSM_GRID
      CkPrintf("MSM C1 HERMITE GRID CUTOFF %d setup section on PE %d\n",
          thisIndex, CkMyPe());
#endif
      CkGetSectionInfo(cookie, msg);  // init the cookie
      msg->get(&msmBlockElementProxy);  // get proxy to MsmC1HermiteBlock
      delete msg;
    }

    void compute(GridMsg *gmsg) {
#ifdef DEBUG_MSM_GRID
      printf("MsmC1HermiteGridCutoff %d:  compute()\n", thisIndex);
#endif
#if 0
      // base class consumes this grid message
      MsmGridCutoffKernel<C1Vector,C1Matrix>::compute(gmsg);
#else
      compute_specialized(gmsg);
#endif

#ifdef MSM_TIMING
      double startTime, stopTime;
      startTime = CkWallTimer();
#endif
#ifdef MSM_REDUCE_GRID

      // perform section reduction over potential grids
      CProxy_CkMulticastMgr mcastProxy =
        CkpvAccess(BOCclass_group).multicastMgr;
      CkMulticastMgr *mcastPtr =
        CProxy_CkMulticastMgr(mcastProxy).ckLocalBranch();
      CkCallback cb(CkIndex_MsmC1HermiteBlock::sumReducedPotential(NULL),
          msmBlockElementProxy);
      // sum into "full" sized buffer needed for contribute
      ehfull.reset(0);
      ehfull += eh;
      mcastPtr->contribute(
          ehfull.nn() * sizeof(C1Vector), ehfull.data().buffer(), 
          CkReduction::sum_float, cookie, cb);

#else
      // place eh into message
      const msm::BlockIndex& bindex = ehblockSend.nblock_wrap;
      int msgsz = eh.data().len() * sizeof(C1Vector);
      GridMsg *gm = new(msgsz, sizeof(int)) GridMsg;
      SET_PRIORITY(gm, sequence, priority);
      gm->put(eh, bindex.level, sequence);
      // lookup in ComputeMsmMgr proxy array by level
      mgrLocal->msmC1HermiteBlock[bindex.level](
          bindex.n.i, bindex.n.j, bindex.n.k).addPotential(gm);

#endif // MSM_REDUCE_GRID

#ifdef MSM_TIMING
      stopTime = CkWallTimer();
      mgrLocal->msmTiming[MsmTimer::COMM] += stopTime - startTime;
      mgrLocal->doneTiming();
#endif
    } // compute()

    // try to improve performance of the major computational part
    void compute_specialized(GridMsg *gmsg);

}; // MsmC1HermiteGridCutoff

void MsmC1HermiteGridCutoff::compute_specialized(GridMsg *gmsg) {
#ifdef MSM_TIMING
      double startTime, stopTime;
      startTime = CkWallTimer();
#endif
      //
      // receive block of charges
      //
      int pid;
      // qh is resized only the first time, memory allocation persists
      gmsg->get(qh, pid, sequence);
      delete gmsg;
#ifdef MSM_TIMING
      stopTime = CkWallTimer();
      mgrLocal->msmTiming[MsmTimer::COMM] += stopTime - startTime;
#endif

      //
      // grid cutoff calculation
      // this charge block -> this potential block
      //

#ifdef MSM_TIMING
      startTime = stopTime;
#endif
      // resets indexing on block
      eh.init(ehblockSend.nrange);  // (always have to re-init nrange for eh)
      eh.reset(0);
      // index range of weights
      int gia = pgc->ia();
      int gib = pgc->ib();
      int gja = pgc->ja();
      int gjb = pgc->jb();
      int gka = pgc->ka();
      int gkb = pgc->kb();
      int gni = pgc->ni();
      int gnj = pgc->nj();
      // index range of charge grid
      int qia = qh.ia();
      int qib = qh.ib();
      int qja = qh.ja();
      int qjb = qh.jb();
      int qka = qh.ka();
      int qkb = qh.kb();
      int qni = qh.ni();
      int qnj = qh.nj();
      // index range of potentials
      int ia = eh.ia();
      int ib = eh.ib();
      int ja = eh.ja();
      int jb = eh.jb();
      int ka = eh.ka();
      int kb = eh.kb();

      int index = 0;

      // access buffers directly
      const C1Matrix *gcbuffer = pgc->data().buffer();
      const C1Vector *qhbuffer = qh.data().buffer();
      C1Vector *ehbuffer = eh.data().buffer();
#ifdef DEBUG_MEMORY_ALIGNMENT
      printf("gcbuffer mem:  addr=%p  div32=%lu  mod32=%lu\n",
          gcbuffer,
          (unsigned long)(gcbuffer)/32,
          (unsigned long)(gcbuffer)%32);
      printf("qhbuffer mem:  addr=%p  div32=%lu  mod32=%lu\n",
          qhbuffer,
          (unsigned long)(qhbuffer)/32,
          (unsigned long)(qhbuffer)%32);
      printf("ehbuffer mem:  addr=%p  div32=%lu  mod32=%lu\n",
          ehbuffer,
          (unsigned long)(ehbuffer)/32,
          (unsigned long)(ehbuffer)%32);
#endif

#ifndef MSM_COMM_ONLY
      // loop over potentials
      for (int k = ka;  k <= kb;  k++) {
        // clip charges to weights along k
        int mka = ( qka >= gka + k ? qka : gka + k );
        int mkb = ( qkb <= gkb + k ? qkb : gkb + k );

        for (int j = ja;  j <= jb;  j++) {
          // clip charges to weights along j
          int mja = ( qja >= gja + j ? qja : gja + j );
          int mjb = ( qjb <= gjb + j ? qjb : gjb + j );

          for (int i = ia;  i <= ib;  i++) {
            // clip charges to weights along i
            int mia = ( qia >= gia + i ? qia : gia + i );
            int mib = ( qib <= gib + i ? qib : gib + i );

            // accumulate sum to this eh point
            C1Vector ehsum = 0;

            // loop over charge grid
            int nn = mib - mia + 1;

            {
              int qnji = qnj * qni;
              int qkoff = -qka*qnji - qja*qni - qia + mia;
              int gnji = gnj * gni;
              int gkoff = (-k-gka)*gnji + (-j-gja)*gni - i - gia + mia;

              for (int qk = mka;  qk <= mkb;  qk++) {
                int qjkoff = qkoff + qk*qnji;
                int gjkoff = gkoff + qk*gnji;

                for (int qj = mja;  qj <= mjb;  qj++) {
                  const C1Vector *qbuf = qhbuffer + (qjkoff + qj*qni);
                  const C1Matrix *gbuf = gcbuffer + (gjkoff + qj*gni);
#ifdef MSM_PROFILING
                  mgrLocal->xLoopCnt[nn]++;
#endif
// help the vectorizer make reasonable decisions
#if defined(__INTEL_COMPILER)
#pragma vector always 
#endif
                  for (int ii = 0;  ii < nn;  ii++) {

#if 0
                    ehsum += gbuf[ii] * qbuf[ii];
#else
                    // skip matvec when matrix is 0
                    // first matrix element tells us if this is the case
                    if ( *((int *)(gbuf)) != 0) {

                      // expand matrix-vector multiply
#if defined(__INTEL_COMPILER)
#pragma vector always
#endif
                      for (int km=0, jm=0;  jm < C1_VECTOR_SIZE;  jm++) {
                        for (int im=0;  im < C1_VECTOR_SIZE;  im++, km++) {
                          ehsum.velem[jm] += gbuf->melem[km] * qbuf->velem[im];
                        }
                      }
                    } // if
                    gbuf++;
                    qbuf++;
#endif
                  }
                }
              } // end loop over charge grid

            }

            ehbuffer[index] = ehsum;
            index++;
          }
        }
      } // end loop over potentials
#endif // !MSM_COMM_ONLY

#ifdef MSM_PROFILING
      mgrLocal->doneProfiling();
#endif

      //
      // send block of potentials
      //

#ifdef MSM_FOLD_FACTOR
      // if "fold factor" is active for this level,
      // need to sum unfolded potential grid back into periodic grid
      if (isfold) {
        // copy unfolded grid
        ehfold = eh;
        // reset eh indexing to correctly folded size
        eh.set(eia, eni, eja, enj, eka, enk);
        eh.reset(0);
#ifdef DEBUG_MSM_GRID
        printf("level=%d   ehfold:  [%d..%d] x [%d..%d] x [%d..%d]  "
            "(%d x %d x %d)\n"
                "              eh:  [%d..%d] x [%d..%d] x [%d..%d]  "
            "(%d x %d x %d)\n"
               "         eh lower:  %d %d %d\n",
            qhblockIndex.level,
            ehfold.ia(), ehfold.ib(), 
            ehfold.ja(), ehfold.jb(),
            ehfold.ka(), ehfold.kb(),
            ehfold.ni(), ehfold.nj(), ehfold.nk(),
            eh.ia(), eh.ib(), 
            eh.ja(), eh.jb(),
            eh.ka(), eh.kb(),
            eh.ni(), eh.nj(), eh.nk(),
            ehblockSend.nrange_wrap.lower().i,
            ehblockSend.nrange_wrap.lower().j,
            ehblockSend.nrange_wrap.lower().k
            );
#endif
        const C1Vector *ehfoldbuf = ehfold.data().buffer();
        C1Vector *ehbuf = eh.data().buffer();
        // now we "fold" eh by calculating the
        // wrap around sum of ehfold into correctly sized eh
        int index = 0;
        for (int k = ka;  k <= kb;  k++) {
          int kk = k;
          if      (kk < eka)  do { kk += enk; } while (kk < eka);
          else if (kk > ekb)  do { kk -= enk; } while (kk > ekb);
          int koff = (kk - eka) * enj;
          for (int j = ja;  j <= jb;  j++) {
            int jj = j;
            if      (jj < eja)  do { jj += enj; } while (jj < eja);
            else if (jj > ejb)  do { jj -= enj; } while (jj > ejb);
            int jkoff = (koff + (jj - eja)) * eni;
            for (int i = ia;  i <= ib;  i++, index++) {
              int ii = i;
              if      (ii < eia)  do { ii += eni; } while (ii < eia);
              else if (ii > eib)  do { ii -= eni; } while (ii > eib);
              int ijkoff = jkoff + (ii - eia);
              ehbuf[ijkoff] += ehfoldbuf[index];
            }
          }
        }
      }
      else {
        // shift grid index range to its true (wrapped) values
        eh.updateLower( ehblockSend.nrange_wrap.lower() );
      }
#else    // !MSM_FOLD_FACTOR
      // shift grid index range to its true (wrapped) values
      eh.updateLower( ehblockSend.nrange_wrap.lower() );
#endif   // MSM_FOLD_FACTOR

#ifdef MSM_TIMING
      stopTime = CkWallTimer();
      mgrLocal->msmTiming[MsmTimer::GRIDCUTOFF] += stopTime - startTime;
#endif
} // MsmC1HermiteGridCutoff::compute_specialized()

// MsmGridCutoff
//
/////////////////


/////////////////
//
// MsmBlock
//
// Performs restriction and prolongation.
//
// Each level of the MSM grid hierarchy is partitioned into MsmBlocks,
// holding both charge and potential grid blocks.
//
// The MsmBlockKernel provides templated routines for the MSM 
// restriction and prolongation algorithms.  Overall is very small 
// part of computational work (less than 2% total for C1 Hermite, 
// less than 4% total for cubic).
// XXX Could be made faster with factored restriction and prolongation 
// algorithms --- especially important for higher order or for 
// generalizing to coarser grid spacing that is not 2h.
// XXX Haven't yet determined factorization for C1 Hermite.
//
// The classes that inherit from MsmBlockKernel provide 
// 3D chare array elements for each level with significant management:
// - receive and sum charges from below
//   (either PatchData or lower level MsmBlock)
// - calculate restriction to 2h grid
// - send up (if not on highest level)
// - section broadcast to MsmGridCutoff
// - receive and sum potentials from above and from 
//   section reduction of MsmGridCutoff
// - calculate prolongation to (1/2)h grid and send down,
//   OR send to PatchData
//
// XXX Grid cutoff calculation below is now replaced with 
// MsmGridCutoff to provide enough parallel work units.
// 

template <class Vtype, class Mtype>
class MsmBlockKernel {
  public:
    CProxy_ComputeMsmMgr mgrProxy;
    ComputeMsmMgr *mgrLocal;  // for quick access to data
    msm::Map *map;
    msm::BlockDiagram *bd;
    msm::Grid<Vtype> qh;
    msm::Grid<Vtype> eh;
#ifndef MSM_GRID_CUTOFF_DECOMP
    const msm::Grid<Mtype> *gcWeights;
    msm::Grid<Vtype> ehCutoff;
#endif
    const msm::Grid<Mtype> *resStencil;
    const msm::Grid<Mtype> *proStencil;
    msm::Grid<Vtype> qhRestricted;
    msm::Grid<Vtype> ehProlongated;
    int cntRecvsCharge;
    int cntRecvsPotential;
    msm::BlockIndex blockIndex;
 
    msm::Grid<Vtype> subgrid;

    int sequence;  // from incoming message for message priority

    MsmBlockKernel(const msm::BlockIndex&);
    MsmBlockKernel(CkMigrateMessage *m) { }

    void init();

#ifndef MSM_GRID_CUTOFF_DECOMP
    void setupStencils(
        const msm::Grid<Mtype> *res,
        const msm::Grid<Mtype> *pro,
        const msm::Grid<Mtype> *gc
        )
    {
      resStencil = res;
      proStencil = pro;
      gcWeights = gc;
    }
#else
    void setupStencils(
        const msm::Grid<Mtype> *res,
        const msm::Grid<Mtype> *pro
        )
    {
      resStencil = res;
      proStencil = pro;
    }
#endif

    void restrictionKernel();
#ifndef MSM_GRID_CUTOFF_DECOMP
    void gridCutoffKernel();
#endif
    void prolongationKernel();

}; // class MsmBlockKernel<Vtype,Mtype>

template <class Vtype, class Mtype>
MsmBlockKernel<Vtype,Mtype>::MsmBlockKernel(const msm::BlockIndex& bindex) {
  blockIndex = bindex;
  mgrProxy = CProxy_ComputeMsmMgr(CkpvAccess(BOCclass_group).computeMsmMgr);
  mgrLocal = CProxy_ComputeMsmMgr::ckLocalBranch(
      CkpvAccess(BOCclass_group).computeMsmMgr);
  map = &(mgrLocal->mapData());
  bd = &(map->blockLevel[blockIndex.level](blockIndex.n));
  qh.init( bd->nrange );
  eh.init( bd->nrange );
#ifndef MSM_GRID_CUTOFF_DECOMP
  ehCutoff.init( bd->nrangeCutoff );
#endif
  qhRestricted.init( bd->nrangeRestricted );
  ehProlongated.init( bd->nrangeProlongated );
#ifdef DEBUG_MSM_GRID
  printf("MsmBlockKernel level=%d, n=%d %d %d:  constructor\n",
      blockIndex.level, blockIndex.n.i, blockIndex.n.j, blockIndex.n.k);
#endif
#ifdef MSM_TIMING
  mgrLocal->addTiming();
#endif
  init();
} // MsmBlockKernel<Vtype,Mtype>::MsmBlockKernel()


template <class Vtype, class Mtype>
void MsmBlockKernel<Vtype,Mtype>::init() {
  qh.reset(0);
  eh.reset(0);
#ifndef MSM_GRID_CUTOFF_DECOMP
  ehCutoff.reset(0);
#endif
  qhRestricted.reset(0);
  ehProlongated.reset(0);
  cntRecvsCharge = 0;
  cntRecvsPotential = 0;
} // MsmBlockKernel<Vtype,Mtype>::init()


template <class Vtype, class Mtype>
void MsmBlockKernel<Vtype,Mtype>::restrictionKernel()
{
#ifdef DEBUG_MSM_GRID
  printf("MsmBlockKernel level=%d, id=%d %d %d:  restriction\n",
      blockIndex.level, blockIndex.n.i, blockIndex.n.j, blockIndex.n.k);
#endif

#ifdef MSM_TIMING
  double startTime, stopTime;
  startTime = CkWallTimer();
#endif

#ifndef MSM_COMM_ONLY
  // stencil data for approximating charge on restricted grid
  const int approx = mgrLocal->approx;
  const int nstencil = ComputeMsmMgr::Nstencil[approx];
  const int *offset = ComputeMsmMgr::IndexOffset[approx];
  const msm::Grid<Mtype>& res = *resStencil;

  // index range for h grid charges
  int ia1 = qh.ia();
  int ib1 = qh.ib();
  int ja1 = qh.ja();
  int jb1 = qh.jb();
  int ka1 = qh.ka();
  int kb1 = qh.kb();

  // index range for restricted (2h) grid charges
  int ia2 = qhRestricted.ia();
  int ib2 = qhRestricted.ib();
  int ja2 = qhRestricted.ja();
  int jb2 = qhRestricted.jb();
  int ka2 = qhRestricted.ka();
  int kb2 = qhRestricted.kb();

  // reset grid
  qhRestricted.reset(0);

  // loop over restricted (2h) grid
  for (int k2 = ka2;  k2 <= kb2;  k2++) {
    int k1 = 2 * k2;
    for (int j2 = ja2;  j2 <= jb2;  j2++) {
      int j1 = 2 * j2;
      for (int i2 = ia2;  i2 <= ib2;  i2++) {
        int i1 = 2 * i2;

        // loop over stencils on h grid
        Vtype& q2hsum = qhRestricted(i2,j2,k2);

        for (int k = 0;  k < nstencil;  k++) {
          int kn = k1 + offset[k];
          if      (kn < ka1) continue;
          else if (kn > kb1) break;

          for (int j = 0;  j < nstencil;  j++) {
            int jn = j1 + offset[j];
            if      (jn < ja1) continue;
            else if (jn > jb1) break;

            for (int i = 0;  i < nstencil;  i++) {
              int in = i1 + offset[i];
              if      (in < ia1) continue;
              else if (in > ib1) break;

              q2hsum += res(i,j,k) * qh(in,jn,kn);
            }
          }
        } // end loop over stencils on h grid

      }
    }
  } // end loop over restricted (2h) grid
#else
  qhRestricted.reset(0);
#endif // !MSM_COMM_ONLY

#ifdef MSM_TIMING
  stopTime = CkWallTimer();
  mgrLocal->msmTiming[MsmTimer::RESTRICT] += stopTime - startTime;
#endif
} // MsmBlockKernel<Vtype,Mtype>::restrictionKernel()


#ifndef MSM_GRID_CUTOFF_DECOMP
template <class Vtype, class Mtype>
void MsmBlockKernel<Vtype,Mtype>::gridCutoffKernel()
{
#ifdef DEBUG_MSM_GRID
  printf("MsmBlockKernel level=%d, id=%d %d %d:  grid cutoff\n",
      blockIndex.level, blockIndex.n.i, blockIndex.n.j, blockIndex.n.k);
#endif
#ifdef MSM_TIMING
  double startTime, stopTime;
  startTime = CkWallTimer();
#endif
#ifndef MSM_COMM_ONLY
  // need grid of weights for this level
  msm::Grid<Mtype>& gc = *gcWeights;
  // index range of weights
  int gia = gc.ia();
  int gib = gc.ib();
  int gja = gc.ja();
  int gjb = gc.jb();
  int gka = gc.ka();
  int gkb = gc.kb();
  // index range of charge grid
  int qia = qh.ia();
  int qib = qh.ib();
  int qja = qh.ja();
  int qjb = qh.jb();
  int qka = qh.ka();
  int qkb = qh.kb();
  // index range of potentials
  int ia = ehCutoff.ia();
  int ib = ehCutoff.ib();
  int ja = ehCutoff.ja();
  int jb = ehCutoff.jb();
  int ka = ehCutoff.ka();
  int kb = ehCutoff.kb();
  // reset grid
  ehCutoff.reset(0);
  // loop over potentials
  for (int k = ka;  k <= kb;  k++) {
    for (int j = ja;  j <= jb;  j++) {
      for (int i = ia;  i <= ib;  i++) {
        // clip charges to weights
        int mia = ( qia >= gia + i ? qia : gia + i );
        int mib = ( qib <= gib + i ? qib : gib + i );
        int mja = ( qja >= gja + j ? qja : gja + j );
        int mjb = ( qjb <= gjb + j ? qjb : gjb + j );
        int mka = ( qka >= gka + k ? qka : gka + k );
        int mkb = ( qkb <= gkb + k ? qkb : gkb + k );
        // accumulate sum to this eh point
        Vtype& ehsum = ehCutoff(i,j,k);
        // loop over smaller charge grid
        for (int qk = mka;  qk <= mkb;  qk++) {
          for (int qj = mja;  qj <= mjb;  qj++) {
            for (int qi = mia;  qi <= mib;  qi++) {
              ehsum += gc(qi-i, qj-j, qk-k) * qh(qi,qj,qk);
            }
          }
        } // end loop over smaller charge grid

      }
    }
  } // end loop over potentials
#else
  ehCutoff.reset(0);
#endif // !MSM_COMM_ONLY
#ifdef MSM_TIMING
  stopTime = CkWallTimer();
  mgrLocal->msmTiming[MsmTimer::GRIDCUTOFF] += stopTime - startTime;
#endif
} // MsmBlockKernel<Vtype,Mtype>::gridCutoffKernel()
#endif // MSM_GRID_CUTOFF_DECOMP


template <class Vtype, class Mtype>
void MsmBlockKernel<Vtype,Mtype>::prolongationKernel()
{
#ifdef DEBUG_MSM_GRID
  printf("MsmBlockKernel level=%d, id=%d %d %d:  prolongation\n",
      blockIndex.level, blockIndex.n.i, blockIndex.n.j, blockIndex.n.k);
#endif

#ifdef MSM_TIMING
  double startTime, stopTime;
  startTime = CkWallTimer();
#endif
#ifndef MSM_COMM_ONLY
  // stencil data for approximating potential on prolongated grid
  const int approx = mgrLocal->approx;
  const int nstencil = ComputeMsmMgr::Nstencil[approx];
  const int *offset = ComputeMsmMgr::IndexOffset[approx];
  const msm::Grid<Mtype>& pro = *proStencil;

  // index range for prolongated h grid potentials
  int ia1 = ehProlongated.ia();
  int ib1 = ehProlongated.ib();
  int ja1 = ehProlongated.ja();
  int jb1 = ehProlongated.jb();
  int ka1 = ehProlongated.ka();
  int kb1 = ehProlongated.kb();

  // index range for 2h grid potentials
  int ia2 = eh.ia();
  int ib2 = eh.ib();
  int ja2 = eh.ja();
  int jb2 = eh.jb();
  int ka2 = eh.ka();
  int kb2 = eh.kb();

  // loop over 2h grid
  for (int k2 = ka2;  k2 <= kb2;  k2++) {
    int k1 = 2 * k2;
    for (int j2 = ja2;  j2 <= jb2;  j2++) {
      int j1 = 2 * j2;
      for (int i2 = ia2;  i2 <= ib2;  i2++) {
        int i1 = 2 * i2;

        // loop over stencils on prolongated h grid
        for (int k = 0;  k < nstencil;  k++) {
          int kn = k1 + offset[k];
          if      (kn < ka1) continue;
          else if (kn > kb1) break;

          for (int j = 0;  j < nstencil;  j++) {
            int jn = j1 + offset[j];
            if      (jn < ja1) continue;
            else if (jn > jb1) break;

            for (int i = 0;  i < nstencil;  i++) {
              int in = i1 + offset[i];
              if      (in < ia1) continue;
              else if (in > ib1) break;

              ehProlongated(in,jn,kn) += pro(i,j,k) * eh(i2,j2,k2);
            }
          }
        } // end loop over stencils on prolongated h grid

      }
    }
  } // end loop over 2h grid
#else
  ehProlongated.reset(0);
#endif // !MSM_COMM_ONLY
#ifdef MSM_TIMING
  stopTime = CkWallTimer();
  mgrLocal->msmTiming[MsmTimer::PROLONGATE] += stopTime - startTime;
#endif
} // MsmBlockKernel<Vtype,Mtype>::prolongationKernel()


//
// MsmBlock handles grids of function values only
// (for cubic, quintic, etc., approximation)
//
class MsmBlock :
  public CBase_MsmBlock,
  public MsmBlockKernel<Float,Float>
{
  public:
    CProxySection_MsmGridCutoff msmGridCutoffBroadcast;
    CProxySection_MsmGridCutoff msmGridCutoffReduction;
 
    MsmBlock(int level) :
      MsmBlockKernel<Float,Float>(
          msm::BlockIndex(level,
            msm::Ivec(thisIndex.x, thisIndex.y, thisIndex.z))
          )
    {
#ifndef MSM_GRID_CUTOFF_DECOMP
      setupStencils(&(map->grespro), &(map->grespro), &(map->gc[level]));
#else
      setupStencils(&(map->grespro), &(map->grespro));
#endif
    }
    MsmBlock(CkMigrateMessage *m) : MsmBlockKernel<Float,Float>(m) { }

    void setupSections();

    void sumReducedPotential(CkReductionMsg *msg) {
#ifdef MSM_TIMING
      double startTime, stopTime;
      startTime = CkWallTimer();
#endif
      msm::Grid<Float> ehfull;
      ehfull.init( msm::IndexRange(eh) );
      memcpy(ehfull.data().buffer(), msg->getData(), msg->getSize());
      delete msg;
      int priority = mgrLocal->nlevels
        + 2*(mgrLocal->nlevels - blockIndex.level)-1;
      int msgsz = ehfull.data().len() * sizeof(Float);
      GridMsg *gm = new(msgsz, sizeof(int)) GridMsg;
      SET_PRIORITY(gm, sequence, MSM_PRIORITY + priority);
      gm->put(ehfull, blockIndex.level, sequence);  // send my level
#ifdef MSM_TIMING
      stopTime = CkWallTimer();
      mgrLocal->msmTiming[MsmTimer::COMM] += stopTime - startTime;
#endif
      addPotential(gm);
    }

    void addCharge(GridMsg *);  // entry

    void restriction() {
      restrictionKernel();
      sendUpCharge();
    }
    void sendUpCharge();
    void gridCutoff();
#ifndef MSM_GRID_CUTOFF_DECOMP
    void sendAcrossPotential();
#endif

    void addPotential(GridMsg *);  // entry

    void prolongation() {
      prolongationKernel();
      sendDownPotential();
    }
    void sendDownPotential();
    void sendPatch();
}; // class MsmBlock


void MsmBlock::setupSections()
{
#ifdef DEBUG_MSM_GRID
  CkPrintf("LEVEL %d MSM BLOCK (%d,%d,%d):  "
      "creating broadcast section on PE %d\n",
      blockIndex.level, thisIndex.x, thisIndex.y, thisIndex.z, CkMyPe());
#endif
  CkVec<CkArrayIndex1D> elems;
  for (int n = 0;  n < bd->indexGridCutoff.len();  n++) {
    elems.push_back(CkArrayIndex1D( bd->indexGridCutoff[n] ));
  }
  msmGridCutoffBroadcast = CProxySection_MsmGridCutoff::ckNew(
      mgrLocal->msmGridCutoff, elems.getVec(), elems.size()
      );
  CProxy_CkMulticastMgr mcastProxy = CkpvAccess(BOCclass_group).multicastMgr;
  CkMulticastMgr *mcastPtr = CProxy_CkMulticastMgr(mcastProxy).ckLocalBranch();
  msmGridCutoffBroadcast.ckSectionDelegate(mcastPtr);
  
#ifdef DEBUG_MSM_GRID
  char s[1024];
  sprintf(s, "LEVEL %d MSM BLOCK (%d,%d,%d):  "
      "creating reduction section on PE %d\n",
      blockIndex.level, thisIndex.x, thisIndex.y, thisIndex.z, CkMyPe());
#endif
  CkVec<CkArrayIndex1D> elems2;
#ifdef DEBUG_MSM_GRID
  strcat(s, "receiving from MsmGridCutoff ID:");
#endif
  for (int n = 0;  n < bd->recvGridCutoff.len();  n++) {
#ifdef DEBUG_MSM_GRID
    char t[20];
    sprintf(t, "  %d", bd->recvGridCutoff[n]);
    strcat(s, t);
#endif
    elems2.push_back(CkArrayIndex1D( bd->recvGridCutoff[n] ));
  }
#ifdef DEBUG_MSM_GRID
  strcat(s, "\n");
  CkPrintf(s);
#endif
  msmGridCutoffReduction = CProxySection_MsmGridCutoff::ckNew(
      mgrLocal->msmGridCutoff, elems2.getVec(), elems2.size()
      );
  msmGridCutoffReduction.ckSectionDelegate(mcastPtr);
  MsmGridCutoffSetupMsg *msg = new MsmGridCutoffSetupMsg;
  CProxyElement_MsmBlock thisElementProxy = thisProxy(thisIndex);
  msg->put(&thisElementProxy);

  msmGridCutoffReduction.setupSections(msg);  // broadcast to entire section

  /* XXX alternatively, setup default reduction client
   *
  mcastPtr->setReductionClient(msmGridCutoffReduction,
      new CkCallback(CkIndex_MsmBlock::myReductionEntry(NULL),
        thisElementProxy));
   *
   */
}


void MsmBlock::addCharge(GridMsg *gm)
{
#ifdef MSM_TIMING
  double startTime, stopTime;
  startTime = CkWallTimer();
#endif
  int pid;
  gm->get(subgrid, pid, sequence);
  delete gm;
  qh += subgrid;
#ifdef MSM_TIMING
  stopTime = CkWallTimer();
  mgrLocal->msmTiming[MsmTimer::COMM] += stopTime - startTime;
#endif
  if (++cntRecvsCharge == bd->numRecvsCharge) {
    int nlevels = mgrLocal->numLevels();
    if (blockIndex.level < nlevels-1) {
      restriction();
    }
    gridCutoff();
  }
} // MsmBlock::addCharge()


void MsmBlock::sendUpCharge()
{
#ifdef MSM_TIMING
  double startTime, stopTime;
  startTime = CkWallTimer();
#endif
  int lnext = blockIndex.level + 1;
  // buffer portions of grid to send to Blocks on next level
  for (int n = 0;  n < bd->sendUp.len();  n++) {
    // initialize the proper subgrid indexing range
    subgrid.init( bd->sendUp[n].nrange );
    // extract the values from the larger grid into the subgrid
    qhRestricted.extract(subgrid);
    // translate the subgrid indexing range to match the MSM block
    subgrid.updateLower( bd->sendUp[n].nrange_wrap.lower() );
    // add the subgrid charges into the block
    msm::BlockIndex& bindex = bd->sendUp[n].nblock_wrap;
    ASSERT(bindex.level == lnext);
    // place subgrid into message
    // SET MESSAGE PRIORITY
    int msgsz = subgrid.nn() * sizeof(Float);
    GridMsg *gm = new(msgsz, sizeof(int)) GridMsg;
    SET_PRIORITY(gm, sequence, MSM_PRIORITY + lnext);
    gm->put(subgrid, blockIndex.level, sequence);  // send my level
    // lookup in ComputeMsmMgr proxy array by level
    mgrLocal->msmBlock[lnext](
        bindex.n.i, bindex.n.j, bindex.n.k).addCharge(gm);
  } // for
#ifdef MSM_TIMING
  stopTime = CkWallTimer();
  mgrLocal->msmTiming[MsmTimer::COMM] += stopTime - startTime;
#endif
} // MsmBlock::sendUpCharge()


void MsmBlock::gridCutoff()
{
#ifdef DEBUG_MSM_GRID
  printf("MsmBlock level=%d, id=%d %d %d:  grid cutoff\n",
      blockIndex.level, blockIndex.n.i, blockIndex.n.j, blockIndex.n.k);
#endif
#ifndef MSM_GRID_CUTOFF_DECOMP
  gridCutoffKernel();
  sendAcrossPotential();
#else // MSM_GRID_CUTOFF_DECOMP

  // send charge block to MsmGridCutoff compute objects
#ifdef MSM_TIMING
  double startTime, stopTime;
  startTime = CkWallTimer();
#endif
  int priority = mgrLocal->nlevels + 2*(mgrLocal->nlevels - blockIndex.level)-1;
  int msgsz = qh.data().len() * sizeof(Float);
  int len = bd->indexGridCutoff.len();

#if 0
  // send charge message to each MsmGridCutoff compute element in list
  for (int n = 0;  n < len;  n++) {
#ifdef MSM_TIMING
    startTime = CkWallTimer();
#endif
    int index = bd->indexGridCutoff[n];
    GridMsg *gm = new(msgsz, sizeof(int)) GridMsg;
    SET_PRIORITY(gm, sequence, MSM_PRIORITY + priority);
    gm->put(qh, blockIndex.level, sequence);  // send my level
#ifdef MSM_TIMING
    stopTime = CkWallTimer();
    mgrLocal->msmTiming[MsmTimer::COMM] += stopTime - startTime;
#endif
    mgrLocal->msmGridCutoff[index].compute(gm);
  }
#else

  // broadcast charge message to section
  GridMsg *gm = new(msgsz, sizeof(int)) GridMsg;
  SET_PRIORITY(gm, sequence, MSM_PRIORITY + priority);
  gm->put(qh, blockIndex.level, sequence);  // send my level
  msmGridCutoffBroadcast.compute(gm);
#ifdef MSM_TIMING
  stopTime = CkWallTimer();
  mgrLocal->msmTiming[MsmTimer::COMM] += stopTime - startTime;
#endif

#endif // 0

#endif // MSM_GRID_CUTOFF_DECOMP

} // MsmBlock::gridCutoff()


#ifndef MSM_GRID_CUTOFF_DECOMP
void MsmBlock::sendAcrossPotential()
{
#ifdef MSM_TIMING
  double startTime, stopTime;
  startTime = CkWallTimer();
#endif
  int lnext = blockIndex.level;
  int priority = mgrLocal->nlevels + 2*(mgrLocal->nlevels - blockIndex.level)-1;
  // buffer portions of grid to send to Blocks on this level
  for (int n = 0;  n < bd->sendAcross.len();  n++) {
    // initialize the proper subgrid indexing range
    subgrid.init( bd->sendAcross[n].nrange );
    // extract the values from the larger grid into the subgrid
    ehCutoff.extract(subgrid);
    // translate the subgrid indexing range to match the MSM block
    subgrid.updateLower( bd->sendAcross[n].nrange_wrap.lower() );
    // add the subgrid charges into the block
    msm::BlockIndex& bindex = bd->sendAcross[n].nblock_wrap;
    ASSERT(bindex.level == lnext);
    // place subgrid into message
    int msgsz = subgrid.nn() * sizeof(Float);
    GridMsg *gm = new(msgsz, sizeof(int)) GridMsg;
    SET_PRIORITY(gm, sequence, MSM_PRIORITY + priority);
    gm->put(subgrid, blockIndex.level, sequence);  // send my level
    // lookup in ComputeMsmMgr proxy array by level
    mgrLocal->msmBlock[lnext](
        bindex.n.i, bindex.n.j, bindex.n.k).addPotential(gm);
  } // for
#ifdef MSM_TIMING
  stopTime = CkWallTimer();
  mgrLocal->msmTiming[MsmTimer::COMM] += stopTime - startTime;
#endif
} // MsmBlock::sendAcrossPotential()
#endif


void MsmBlock::addPotential(GridMsg *gm)
{
#ifdef MSM_TIMING
  double startTime, stopTime;
  startTime = CkWallTimer();
#endif
  int pid;
  int pseq;
  gm->get(subgrid, pid, pseq);  // receive sender's level
  delete gm;
  eh += subgrid;
#ifdef MSM_TIMING
  stopTime = CkWallTimer();
  mgrLocal->msmTiming[MsmTimer::COMM] += stopTime - startTime;
#endif
  if (++cntRecvsPotential == bd->numRecvsPotential) {
    if (blockIndex.level > 0) {
      prolongation();
    }
    else {
      sendPatch();
    }
  }
} // MsmBlock::addPotential()


void MsmBlock::sendDownPotential()
{
#ifdef MSM_TIMING
  double startTime, stopTime;
  startTime = CkWallTimer();
#endif
  int lnext = blockIndex.level - 1;
  int priority = mgrLocal->nlevels + 2*(mgrLocal->nlevels - blockIndex.level);
  // buffer portions of grid to send to Blocks on next level
  for (int n = 0;  n < bd->sendDown.len();  n++) {
    // initialize the proper subgrid indexing range
    subgrid.init( bd->sendDown[n].nrange );
    // extract the values from the larger grid into the subgrid
    ehProlongated.extract(subgrid);
    // translate the subgrid indexing range to match the MSM block
    subgrid.updateLower( bd->sendDown[n].nrange_wrap.lower() );
    // add the subgrid charges into the block
    msm::BlockIndex& bindex = bd->sendDown[n].nblock_wrap;
    ASSERT(bindex.level == lnext);
    // place subgrid into message
    int msgsz = subgrid.nn() * sizeof(Float);
    GridMsg *gm = new(msgsz, sizeof(int)) GridMsg;
    SET_PRIORITY(gm, sequence, MSM_PRIORITY + priority);
    gm->put(subgrid, blockIndex.level, sequence);  // send my level
    // lookup in ComputeMsmMgr proxy array by level
    mgrLocal->msmBlock[lnext](
        bindex.n.i, bindex.n.j, bindex.n.k).addPotential(gm);
  } // for
#ifdef MSM_TIMING
  stopTime = CkWallTimer();
  mgrLocal->msmTiming[MsmTimer::COMM] += stopTime - startTime;
  mgrLocal->doneTiming();
#endif
  init();  // reinitialize for next computation
} // MsmBlock::sendDownPotential()


void MsmBlock::sendPatch()
{
#ifdef MSM_TIMING
  double startTime, stopTime;
  startTime = CkWallTimer();
#endif
  int lnext = blockIndex.level;
  int priority = mgrLocal->nlevels + 2*(mgrLocal->nlevels - blockIndex.level);
  ASSERT(lnext == 0);
  // buffer portions of grid to send to Blocks on next level
  for (int n = 0;  n < bd->sendPatch.len();  n++) {
    // initialize the proper subgrid indexing range
    subgrid.init( bd->sendPatch[n].nrange );
    // extract the values from the larger grid into the subgrid
    eh.extract(subgrid);
    // translate the subgrid indexing range to match the MSM block
    subgrid.updateLower( bd->sendPatch[n].nrange_unwrap.lower() );
    // add the subgrid charges into the block, need its patch ID
    int pid = bd->sendPatch[n].patchID;
    // place subgrid into message
    int msgsz = subgrid.nn() * sizeof(Float);
    GridMsg *gm = new(msgsz, sizeof(int)) GridMsg;
    SET_PRIORITY(gm, sequence, MSM_PRIORITY + priority);
    gm->put(subgrid, pid, sequence);  // send patch ID
    // lookup which PE has this patch
    PatchMap *pm = PatchMap::Object();
    int pe = pm->node(pid);
    mgrProxy[pe].addPotential(gm);
  }
#ifdef MSM_TIMING
  stopTime = CkWallTimer();
  mgrLocal->msmTiming[MsmTimer::COMM] += stopTime - startTime;
  mgrLocal->doneTiming();
#endif
  init();  // reinitialize for next computation
} // MsmBlock::sendPatch()


//
// MsmC1HermiteBlock handles grids of vector elements
// for C1 Hermite approximation
//
class MsmC1HermiteBlock :
  public CBase_MsmC1HermiteBlock,
  public MsmBlockKernel<C1Vector,C1Matrix>
{
  public:
    CProxySection_MsmC1HermiteGridCutoff msmGridCutoffBroadcast;
    CProxySection_MsmC1HermiteGridCutoff msmGridCutoffReduction;
 
    MsmC1HermiteBlock(int level) :
      MsmBlockKernel<C1Vector,C1Matrix>(
          msm::BlockIndex(level,
            msm::Ivec(thisIndex.x, thisIndex.y, thisIndex.z))
          )
    {
      int isfirstlevel = (level == 0);
      int istoplevel = (level == map->gridrange.len()-1);
      const msm::Grid<C1Matrix> *res =
        (istoplevel ? NULL : &(map->gres_c1hermite[level]));
      const msm::Grid<C1Matrix> *pro =
        (isfirstlevel ? NULL : &(map->gpro_c1hermite[level-1]));
#ifndef MSM_GRID_CUTOFF_DECOMP
      const msm::Grid<C1Matrix> *gc = &(map->gc_c1hermite[level]);
      setupStencils(res, pro, gc);
#else
      setupStencils(res, pro);
#endif
    }
    MsmC1HermiteBlock(CkMigrateMessage *m) :
      MsmBlockKernel<C1Vector,C1Matrix>(m) { }

    void setupSections();

    void sumReducedPotential(CkReductionMsg *msg) {
#ifdef MSM_TIMING
      double startTime, stopTime;
      startTime = CkWallTimer();
#endif
      msm::Grid<C1Vector> ehfull;
      ehfull.init( msm::IndexRange(eh) );
      memcpy(ehfull.data().buffer(), msg->getData(), msg->getSize());
      delete msg;
      int priority = mgrLocal->nlevels
        + 2*(mgrLocal->nlevels - blockIndex.level)-1;
      int msgsz = ehfull.data().len() * sizeof(C1Vector);
      GridMsg *gm = new(msgsz, sizeof(int)) GridMsg;
      SET_PRIORITY(gm, sequence, MSM_PRIORITY + priority);
      gm->put(ehfull, blockIndex.level, sequence);  // send my level
#ifdef MSM_TIMING
      stopTime = CkWallTimer();
      mgrLocal->msmTiming[MsmTimer::COMM] += stopTime - startTime;
#endif
      addPotential(gm);
    }

    void addCharge(GridMsg *);  // entry

    void restriction() {
      restrictionKernel();
      sendUpCharge();
    }
    void sendUpCharge();
    void gridCutoff();
#ifndef MSM_GRID_CUTOFF_DECOMP
    void sendAcrossPotential();
#endif

    void addPotential(GridMsg *);  // entry

    void prolongation() {
      prolongationKernel();
      sendDownPotential();
    }
    void sendDownPotential();
    void sendPatch();
}; // class MsmC1HermiteBlock


void MsmC1HermiteBlock::setupSections()
{
#ifdef DEBUG_MSM_GRID
  CkPrintf("LEVEL %d MSM C1 HERMITE BLOCK (%d,%d,%d):  "
      "creating broadcast section on PE %d\n",
      blockIndex.level, thisIndex.x, thisIndex.y, thisIndex.z, CkMyPe());
#endif
  CkVec<CkArrayIndex1D> elems;
  for (int n = 0;  n < bd->indexGridCutoff.len();  n++) {
    elems.push_back(CkArrayIndex1D( bd->indexGridCutoff[n] ));
  }
  msmGridCutoffBroadcast = CProxySection_MsmC1HermiteGridCutoff::ckNew(
      mgrLocal->msmC1HermiteGridCutoff, elems.getVec(), elems.size()
      );
  CProxy_CkMulticastMgr mcastProxy = CkpvAccess(BOCclass_group).multicastMgr;
  CkMulticastMgr *mcastPtr = CProxy_CkMulticastMgr(mcastProxy).ckLocalBranch();
  msmGridCutoffBroadcast.ckSectionDelegate(mcastPtr);
  
#ifdef DEBUG_MSM_GRID
  char s[1024];
  sprintf(s, "LEVEL %d MSM C1 HERMITE BLOCK (%d,%d,%d):  "
      "creating reduction section on PE %d\n",
      blockIndex.level, thisIndex.x, thisIndex.y, thisIndex.z, CkMyPe());
#endif
  CkVec<CkArrayIndex1D> elems2;
#ifdef DEBUG_MSM_GRID
  strcat(s, "receiving from MsmC1HermiteGridCutoff ID:");
#endif
  for (int n = 0;  n < bd->recvGridCutoff.len();  n++) {
#ifdef DEBUG_MSM_GRID
    char t[20];
    sprintf(t, "  %d", bd->recvGridCutoff[n]);
    strcat(s, t);
#endif
    elems2.push_back(CkArrayIndex1D( bd->recvGridCutoff[n] ));
  }
#ifdef DEBUG_MSM_GRID
  strcat(s, "\n");
  CkPrintf(s);
#endif
  msmGridCutoffReduction = CProxySection_MsmC1HermiteGridCutoff::ckNew(
      mgrLocal->msmC1HermiteGridCutoff, elems2.getVec(), elems2.size()
      );
  msmGridCutoffReduction.ckSectionDelegate(mcastPtr);
  MsmC1HermiteGridCutoffSetupMsg *msg = new MsmC1HermiteGridCutoffSetupMsg;
  CProxyElement_MsmC1HermiteBlock thisElementProxy = thisProxy(thisIndex);
  msg->put(&thisElementProxy);

  msmGridCutoffReduction.setupSections(msg);  // broadcast to entire section

  /* XXX alternatively, setup default reduction client
   *
  mcastPtr->setReductionClient(msmGridCutoffReduction,
      new CkCallback(CkIndex_MsmC1HermiteBlock::myReductionEntry(NULL),
        thisElementProxy));
   *
   */
} // MsmC1HermiteBlock::setupSections()


void MsmC1HermiteBlock::addCharge(GridMsg *gm)
{
#ifdef MSM_TIMING
  double startTime, stopTime;
  startTime = CkWallTimer();
#endif
  int pid;
  gm->get(subgrid, pid, sequence);
  delete gm;
  qh += subgrid;
#ifdef MSM_TIMING
  stopTime = CkWallTimer();
  mgrLocal->msmTiming[MsmTimer::COMM] += stopTime - startTime;
#endif
  if (++cntRecvsCharge == bd->numRecvsCharge) {
    int nlevels = mgrLocal->numLevels();
    if (blockIndex.level < nlevels-1) {
      restriction();
    }
    gridCutoff();
  }
} // MsmC1HermiteBlock::addCharge()


void MsmC1HermiteBlock::sendUpCharge()
{
#ifdef MSM_TIMING
  double startTime, stopTime;
  startTime = CkWallTimer();
#endif
  int lnext = blockIndex.level + 1;
  // buffer portions of grid to send to Blocks on next level
  for (int n = 0;  n < bd->sendUp.len();  n++) {
    // initialize the proper subgrid indexing range
    subgrid.init( bd->sendUp[n].nrange );
    // extract the values from the larger grid into the subgrid
    qhRestricted.extract(subgrid);
    // translate the subgrid indexing range to match the MSM block
    subgrid.updateLower( bd->sendUp[n].nrange_wrap.lower() );
    // add the subgrid charges into the block
    msm::BlockIndex& bindex = bd->sendUp[n].nblock_wrap;
    ASSERT(bindex.level == lnext);
    // place subgrid into message
    // SET MESSAGE PRIORITY
    int msgsz = subgrid.nn() * sizeof(C1Vector);
    GridMsg *gm = new(msgsz, sizeof(int)) GridMsg;
    SET_PRIORITY(gm, sequence, MSM_PRIORITY + lnext);
    gm->put(subgrid, blockIndex.level, sequence);  // send my level
    // lookup in ComputeMsmMgr proxy array by level
    mgrLocal->msmC1HermiteBlock[lnext](
        bindex.n.i, bindex.n.j, bindex.n.k).addCharge(gm);
  } // for
#ifdef MSM_TIMING
  stopTime = CkWallTimer();
  mgrLocal->msmTiming[MsmTimer::COMM] += stopTime - startTime;
#endif
} // MsmC1HermiteBlock::sendUpCharge()


void MsmC1HermiteBlock::gridCutoff()
{
#ifdef DEBUG_MSM_GRID
  printf("MsmC1HermiteBlock level=%d, id=%d %d %d:  grid cutoff\n",
      blockIndex.level, blockIndex.n.i, blockIndex.n.j, blockIndex.n.k);
#endif
#ifndef MSM_GRID_CUTOFF_DECOMP
  gridCutoffKernel();
  sendAcrossPotential();
#else // MSM_GRID_CUTOFF_DECOMP

  // send charge block to MsmGridCutoff compute objects
#ifdef MSM_TIMING
  double startTime, stopTime;
  startTime = CkWallTimer();
#endif
  int priority = mgrLocal->nlevels + 2*(mgrLocal->nlevels - blockIndex.level)-1;
  int msgsz = qh.data().len() * sizeof(C1Vector);
  int len = bd->indexGridCutoff.len();

#if 0
  // send charge message to each MsmGridCutoff compute element in list
  for (int n = 0;  n < len;  n++) {
#ifdef MSM_TIMING
    startTime = CkWallTimer();
#endif
    int index = bd->indexGridCutoff[n];
    GridMsg *gm = new(msgsz, sizeof(int)) GridMsg;
    SET_PRIORITY(gm, sequence, MSM_PRIORITY + priority);
    gm->put(qh, blockIndex.level, sequence);  // send my level
#ifdef MSM_TIMING
    stopTime = CkWallTimer();
    mgrLocal->msmTiming[MsmTimer::COMM] += stopTime - startTime;
#endif
    mgrLocal->msmGridCutoff[index].compute(gm);
  }
#else

  // broadcast charge message to section
  GridMsg *gm = new(msgsz, sizeof(int)) GridMsg;
  SET_PRIORITY(gm, sequence, MSM_PRIORITY + priority);
  gm->put(qh, blockIndex.level, sequence);  // send my level
  msmGridCutoffBroadcast.compute(gm);
#ifdef MSM_TIMING
  stopTime = CkWallTimer();
  mgrLocal->msmTiming[MsmTimer::COMM] += stopTime - startTime;
#endif

#endif // 0

#endif // MSM_GRID_CUTOFF_DECOMP

} // MsmC1HermiteBlock::gridCutoff()


#ifndef MSM_GRID_CUTOFF_DECOMP
void MsmC1HermiteBlock::sendAcrossPotential()
{
#ifdef MSM_TIMING
  double startTime, stopTime;
  startTime = CkWallTimer();
#endif
  int lnext = blockIndex.level;
  int priority = mgrLocal->nlevels + 2*(mgrLocal->nlevels - blockIndex.level)-1;
  // buffer portions of grid to send to Blocks on this level
  for (int n = 0;  n < bd->sendAcross.len();  n++) {
    // initialize the proper subgrid indexing range
    subgrid.init( bd->sendAcross[n].nrange );
    // extract the values from the larger grid into the subgrid
    ehCutoff.extract(subgrid);
    // translate the subgrid indexing range to match the MSM block
    subgrid.updateLower( bd->sendAcross[n].nrange_wrap.lower() );
    // add the subgrid charges into the block
    msm::BlockIndex& bindex = bd->sendAcross[n].nblock_wrap;
    ASSERT(bindex.level == lnext);
    // place subgrid into message
    int msgsz = subgrid.nn() * sizeof(C1Vector);
    GridMsg *gm = new(msgsz, sizeof(int)) GridMsg;
    SET_PRIORITY(gm, sequence, MSM_PRIORITY + priority);
    gm->put(subgrid, blockIndex.level, sequence);  // send my level
    // lookup in ComputeMsmMgr proxy array by level
    mgrLocal->msmC1HermiteBlock[lnext](
        bindex.n.i, bindex.n.j, bindex.n.k).addPotential(gm);
  } // for
#ifdef MSM_TIMING
  stopTime = CkWallTimer();
  mgrLocal->msmTiming[MsmTimer::COMM] += stopTime - startTime;
#endif
} // MsmC1HermiteBlock::sendAcrossPotential()
#endif


void MsmC1HermiteBlock::addPotential(GridMsg *gm)
{
#ifdef MSM_TIMING
  double startTime, stopTime;
  startTime = CkWallTimer();
#endif
  int pid;
  int pseq;
  gm->get(subgrid, pid, pseq);  // receive sender's level
  delete gm;
  eh += subgrid;
#ifdef MSM_TIMING
  stopTime = CkWallTimer();
  mgrLocal->msmTiming[MsmTimer::COMM] += stopTime - startTime;
#endif
  if (++cntRecvsPotential == bd->numRecvsPotential) {
    if (blockIndex.level > 0) {
      prolongation();
    }
    else {
      sendPatch();
    }
  }
} // MsmC1HermiteBlock::addPotential()


void MsmC1HermiteBlock::sendDownPotential()
{
#ifdef MSM_TIMING
  double startTime, stopTime;
  startTime = CkWallTimer();
#endif
  int lnext = blockIndex.level - 1;
  int priority = mgrLocal->nlevels + 2*(mgrLocal->nlevels - blockIndex.level);
  // buffer portions of grid to send to Blocks on next level
  for (int n = 0;  n < bd->sendDown.len();  n++) {
    // initialize the proper subgrid indexing range
    subgrid.init( bd->sendDown[n].nrange );
    // extract the values from the larger grid into the subgrid
    ehProlongated.extract(subgrid);
    // translate the subgrid indexing range to match the MSM block
    subgrid.updateLower( bd->sendDown[n].nrange_wrap.lower() );
    // add the subgrid charges into the block
    msm::BlockIndex& bindex = bd->sendDown[n].nblock_wrap;
    ASSERT(bindex.level == lnext);
    // place subgrid into message
    int msgsz = subgrid.nn() * sizeof(C1Vector);
    GridMsg *gm = new(msgsz, sizeof(int)) GridMsg;
    SET_PRIORITY(gm, sequence, MSM_PRIORITY + priority);
    gm->put(subgrid, blockIndex.level, sequence);  // send my level
    // lookup in ComputeMsmMgr proxy array by level
    mgrLocal->msmC1HermiteBlock[lnext](
        bindex.n.i, bindex.n.j, bindex.n.k).addPotential(gm);
  } // for
#ifdef MSM_TIMING
  stopTime = CkWallTimer();
  mgrLocal->msmTiming[MsmTimer::COMM] += stopTime - startTime;
  mgrLocal->doneTiming();
#endif
  init();  // reinitialize for next computation
} // MsmC1HermiteBlock::sendDownPotential()


void MsmC1HermiteBlock::sendPatch()
{
#ifdef MSM_TIMING
  double startTime, stopTime;
  startTime = CkWallTimer();
#endif
  int lnext = blockIndex.level;
  int priority = mgrLocal->nlevels + 2*(mgrLocal->nlevels - blockIndex.level);
  ASSERT(lnext == 0);
  // buffer portions of grid to send to Blocks on next level
  for (int n = 0;  n < bd->sendPatch.len();  n++) {
    // initialize the proper subgrid indexing range
    subgrid.init( bd->sendPatch[n].nrange );
    // extract the values from the larger grid into the subgrid
    eh.extract(subgrid);
    // translate the subgrid indexing range to match the MSM block
    subgrid.updateLower( bd->sendPatch[n].nrange_unwrap.lower() );
    // add the subgrid charges into the block, need its patch ID
    int pid = bd->sendPatch[n].patchID;
    // place subgrid into message
    int msgsz = subgrid.nn() * sizeof(C1Vector);
    GridMsg *gm = new(msgsz, sizeof(int)) GridMsg;
    SET_PRIORITY(gm, sequence, MSM_PRIORITY + priority);
    gm->put(subgrid, pid, sequence);  // send patch ID
    // lookup which PE has this patch
    PatchMap *pm = PatchMap::Object();
    int pe = pm->node(pid);
    mgrProxy[pe].addPotential(gm);
  }
#ifdef MSM_TIMING
  stopTime = CkWallTimer();
  mgrLocal->msmTiming[MsmTimer::COMM] += stopTime - startTime;
  mgrLocal->doneTiming();
#endif
  init();  // reinitialize for next computation
} // MsmC1HermiteBlock::sendPatch()


// MsmBlock
//
//////////////////


ComputeMsmMgr::ComputeMsmMgr() :
  msmProxy(thisgroup), msmCompute(0)
{
#ifdef DEBUG_MSM_VERBOSE
  printf("ComputeMsmMgr:  (constructor) PE %d\n", CkMyPe());
#endif
  CkpvAccess(BOCclass_group).computeMsmMgr = thisgroup;

#ifdef MSM_TIMING
  if (CkMyPe() == 0) {
    msmTimer = CProxy_MsmTimer::ckNew();
  }
  initTiming();
#endif
#ifdef MSM_PROFILING
  if (CkMyPe() == 0) {
    msmProfiler = CProxy_MsmProfiler::ckNew();
  }
  initProfiling();
#endif
}

ComputeMsmMgr::~ComputeMsmMgr()
{
#ifdef DEBUG_MSM_VERBOSE
  printf("ComputeMsmMgr:  (destructor) PE %d\n", CkMyPe());
#endif
  // free memory?
}


//
// Given basis vector length "len" (and with user grid spacing)
// If using periodic boundary conditions along this basis vector,
// h is calculated to be close to desired grid spacing such that 
// nn = 2^k or 3*2^k.  For non-periodic boundaries, we can set h
// to the desired grid spacing, and set ia and ib to pad 1/2 the 
// interpolating stencil width.  
//
void ComputeMsmMgr::setup_hgrid_1d(BigReal len, BigReal& hh, int& nn,
    int& ia, int& ib, int isperiodic)
{
  ASSERT(gridspacing > 0);
  if (isperiodic) {
    const BigReal hmin = (4./5) * gridspacing;
    const BigReal hmax = 1.5 * hmin;
    hh = len;
    nn = 1;  // start with one grid point across length
    while (hh >= hmax) {
      hh *= 0.5;  // halve spacing and double grid points
      nn <<= 1;
    }
    if (hh < hmin) {
      if (nn < 4) {
        NAMD_die("Basis vector is too short or MSM grid spacing is too large");
      }
      hh *= (4./3);  // scale hh by 4/3 and nn by 3/4
      nn >>= 2;
      nn *= 3;
    }
    // now we have:  hmin <= h < hmax,
    // where nn is a power of 2 times no more than one power of 3
    ia = 0;
    ib = nn-1;
  }
  else {
    hh = gridspacing;
    // Instead of "nn = (int) ceil(len / hh);"
    // len is divisible by hh, up to roundoff error, so round to closest nn
    nn = (int) floor(len/hh + 0.5);
    ia = -s_edge;
    ib = nn + s_edge;
  }
} // ComputeMsmMgr::setup_hgrid_1d()


// make sure that block sizes divide evenly into periodic dimensions
// call only for periodic dimensions
void ComputeMsmMgr::setup_periodic_blocksize(int& bsize, int n)
{
  if (n % bsize != 0) {
    // n is either 2^k or 3*2^k
    int newbsize = 1;
    if (n % 3 == 0) newbsize = 3;
    while (newbsize < bsize && newbsize < n) newbsize *= 2;
    if (bsize < newbsize) newbsize /= 2;
    if (n % newbsize != 0) {
      NAMD_die("MSM grid size for periodic dimensions must be "
          "a power of 2 times at most one power of 3");
    }
    bsize = newbsize;
  }
  return;
}


//
// This is the major routine that sets everything up for MSM based on 
// 1. cell basis vectors and/or max and min coordinates plus padding
// 2. cutoff and MSM-related parameters from SimParameter
// Includes determining grid spacings along periodic dimensions, 
// determining grid dimensions and number of levels for system,
// then calculating all needed coefficients for grid cutoff part
// and grid transfer parts (restriction and prolongation).
//
// Then sets up Map for parallel decomposition based on 
// MSM block size parameters from SimParameter.
//
// Then determines chare array element placement of MsmBlock and 
// MsmGridCutoff arrays based on number of PEs and number of nodes.
//
// Then allocates (on PE 0) MsmBlock (3D chare arrays, one per level) 
// and MsmGridCutoff (one 1D chare array for all block-block interactions)
// and then broadcasts array proxies across group.
//
void ComputeMsmMgr::initialize(MsmInitMsg *msg)
{
#ifdef DEBUG_MSM_VERBOSE
  printf("ComputeMsmMgr:  initialize() PE %d\n", CkMyPe());
#endif

  smin = msg->smin;
  smax = msg->smax;
  delete msg;

#if 0
  printf("PE%d: initializing MSM\n", CkMyPe());
#endif

  SimParameters *simParams = Node::Object()->simParameters;

  // get required sim params, check validity
  lattice = simParams->lattice;

  // set user-defined extent of system
  Vector rmin(simParams->MSMxmin, simParams->MSMymin, simParams->MSMzmin);
  Vector rmax(simParams->MSMxmax, simParams->MSMymax, simParams->MSMzmax);
  Vector sdmin = lattice.scale(rmin);
  Vector sdmax = lattice.scale(rmax);
  // swap coordinates between min and max to correct for possible rotation
  if (sdmin.x > sdmax.x) { double t=sdmin.x; sdmin.x=sdmax.x; sdmax.x=t; }
  if (sdmin.y > sdmax.y) { double t=sdmin.y; sdmin.y=sdmax.y; sdmax.y=t; }
  if (sdmin.z > sdmax.z) { double t=sdmin.z; sdmin.z=sdmax.z; sdmax.z=t; }
  // extend smin, smax by user-defined extent, where appropriate
  if ( ! lattice.a_p() && (sdmin.x != 0 || sdmax.x != 0)) {
    if (sdmin.x < smin.x) {
      smin.x = sdmin.x;
      if (CkMyPe() == 0) {
        iout << iINFO << "MSM extending minimum X to "
          << simParams->MSMxmin << " A\n" << endi;
      }
    }
    if (sdmax.x > smax.x) {
      smax.x = sdmax.x;
      if (CkMyPe() == 0) {
        iout << iINFO << "MSM extending maximum X to "
          << simParams->MSMxmax << " A\n" << endi;
      }
    }
  }
  if ( ! lattice.b_p() && (sdmin.y != 0 || sdmax.y != 0)) {
    if (sdmin.y < smin.y) {
      smin.y = sdmin.y;
      if (CkMyPe() == 0) {
        iout << iINFO << "MSM extending minimum Y to "
          << simParams->MSMymin << " A\n" << endi;
      }
    }
    if (sdmax.y > smax.y) {
      smax.y = sdmax.y;
      if (CkMyPe() == 0) {
        iout << iINFO << "MSM extending maximum Y to "
          << simParams->MSMymax << " A\n" << endi;
      }
    }
  }
  if ( ! lattice.c_p() && (sdmin.z != 0 || sdmax.z != 0)) {
    if (sdmin.z < smin.z) {
      smin.z = sdmin.z;
      if (CkMyPe() == 0) {
        iout << iINFO << "MSM extending minimum Z to "
          << simParams->MSMzmin << " A\n" << endi;
      }
    }
    if (sdmax.z > smax.z) {
      smax.z = sdmax.z;
      if (CkMyPe() == 0) {
        iout << iINFO << "MSM extending maximum Z to "
          << simParams->MSMzmax << " A\n" << endi;
      }
    }
  }

#ifdef DEBUG_MSM_VERBOSE
  printf("smin = %g %g %g  smax = %g %g %g\n",
      smin.x, smin.y, smin.z, smax.x, smax.y, smax.z);
#endif

  approx = simParams->MSMApprox;
  if (approx < 0 || approx >= NUM_APPROX) {
    NAMD_die("MSM: unknown approximation requested (MSMApprox)");
  }

  split = simParams->MSMSplit;
  if (split < 0 || split >= NUM_SPLIT) {
    NAMD_die("MSM: unknown splitting requested (MSMSplit)");
  }

  if (CkMyPe() == 0) {
    const char *approx_str, *split_str;
    switch (approx) {
      case CUBIC:      approx_str = "C1 cubic";    break;
      case QUINTIC:    approx_str = "C1 quintic";  break;
      case QUINTIC2:   approx_str = "C2 quintic";  break;
      case SEPTIC:     approx_str = "C1 septic";   break;
      case SEPTIC3:    approx_str = "C3 septic";   break;
      case NONIC:      approx_str = "C1 nonic";    break;
      case NONIC4:     approx_str = "C4 nonic";    break;
      case C1HERMITE:  approx_str = "C1 Hermite";  break;
      default:         approx_str = "unknown";     break;
    }
    switch (split) {
      case TAYLOR2:  split_str = "C2 Taylor";   break;
      case TAYLOR3:  split_str = "C3 Taylor";   break;
      case TAYLOR4:  split_str = "C4 Taylor";   break;
      case TAYLOR5:  split_str = "C5 Taylor";   break;
      case TAYLOR6:  split_str = "C6 Taylor";   break;
      case TAYLOR7:  split_str = "C7 Taylor";   break;
      case TAYLOR8:  split_str = "C8 Taylor";   break;
      default:       split_str = "unknown";     break;
    }
    iout << iINFO << "MSM using "
                  << approx_str << " interpolation\n";
    iout << iINFO << "MSM using "
                  << split_str << " splitting function\n";
  }

  a = simParams->cutoff;

  if (approx == C1HERMITE) {
    gridScalingFactor = 2;
  }
  else {
    gridScalingFactor = 1;
  }

  gridspacing = gridScalingFactor * simParams->MSMGridSpacing;
  if (gridspacing <= 0) {
    NAMD_die("MSM: grid spacing must be greater than 0 (MSMGridSpacing)");
  }
  else if (gridspacing >= a) {
    NAMD_die("MSM: grid spacing must be less than cutoff (MSMGridSpacing)");
  }

  padding = gridScalingFactor * simParams->MSMPadding;
  if (padding < 0) {
    NAMD_die("MSM: padding must be non-negative (MSMPadding)");
  }

  // set maximum number of levels (default 0 adapts levels to system)
  nlevels = simParams->MSMLevels;

  // XXX dispersion unused for now
  dispersion = 0;
  if ( ! dispersion && split >= TAYLOR2_DISP) {
    NAMD_die("MSM: requested splitting for long-range dispersion "
        "(not implemented)");
  }

  // set block sizes for grid decomposition
  int bsx = simParams->MSMBlockSizeX / int(gridScalingFactor);
  int bsy = simParams->MSMBlockSizeY / int(gridScalingFactor);
  int bsz = simParams->MSMBlockSizeZ / int(gridScalingFactor);
  if (bsx <= 0 || bsy <= 0 || bsz <= 0) {
    NAMD_die("MSM: invalid block size requested (MSMBlockSize[XYZ])");
  }
#ifdef MSM_FIXED_SIZE_GRID_MSG
  else if (bsx * bsy * bsz > MSM_MAX_BLOCK_VOLUME) {
    NAMD_die("MSM: requested block size (MSMBlockSize[XYZ]) too big");
  }
#endif
  if (CkMyPe() == 0) {
    iout << iINFO << "MSM block size decomposition along X is "
                  << bsx << " grid points\n";
    iout << iINFO << "MSM block size decomposition along Y is "
                  << bsy << " grid points\n";
    iout << iINFO << "MSM block size decomposition along Z is "
                  << bsz << " grid points\n";
  }

  s_edge = (PolyDegree[approx] - 1) / 2;  // stencil edge size
  omega = 2 * PolyDegree[approx];         // smallest non-periodic grid length

  BigReal xlen, ylen, zlen;
  Vector sgmin, sgmax;  // grid min and max, in scaled coordinates
  int ispx = lattice.a_p();
  int ispy = lattice.b_p();
  int ispz = lattice.c_p();
  int ispany = (ispx || ispy || ispz);  // is there any periodicity?

  if (ispx) {  // periodic along basis vector
    xlen = lattice.a().length();
    sgmax.x = 0.5;
    sgmin.x = -0.5;
  }
  else {  // non-periodic
    sgmax.x = smax.x + padding;  // pad the edges
    sgmin.x = smin.x - padding;
    ASSERT(gridspacing > 0);
    // restrict center to be on a grid point
    BigReal mupper = ceil(sgmax.x / (2*gridspacing));
    BigReal mlower = floor(sgmin.x / (2*gridspacing));
    sgmax.x = 2*gridspacing*mupper;
    sgmin.x = 2*gridspacing*mlower;
    xlen = sgmax.x - sgmin.x;
  }
#ifdef DEBUG_MSM_VERBOSE
  printf("xlen = %g   sgmin.x = %g   sgmax.x = %g\n", xlen, sgmin.x, sgmax.x);
#endif

  if (ispy) {  // periodic along basis vector
    ylen = lattice.b().length();
    sgmax.y = 0.5;
    sgmin.y = -0.5;
  }
  else {  // non-periodic
    sgmax.y = smax.y + padding;  // pad the edges
    sgmin.y = smin.y - padding;
    ASSERT(gridspacing > 0);
    // restrict center to be on a grid point
    BigReal mupper = ceil(sgmax.y / (2*gridspacing));
    BigReal mlower = floor(sgmin.y / (2*gridspacing));
    sgmax.y = 2*gridspacing*mupper;
    sgmin.y = 2*gridspacing*mlower;
    ylen = sgmax.y - sgmin.y;
  }
#ifdef DEBUG_MSM_VERBOSE
  printf("ylen = %g   sgmin.y = %g   sgmax.y = %g\n", ylen, sgmin.y, sgmax.y);
#endif

  if (ispz) {  // periodic along basis vector
    zlen = lattice.c().length();
    sgmax.z = 0.5;
    sgmin.z = -0.5;
  }
  else {  // non-periodic
    sgmax.z = smax.z + padding;  // pad the edges
    sgmin.z = smin.z - padding;
    ASSERT(gridspacing > 0);
    // restrict center to be on a grid point
    BigReal mupper = ceil(sgmax.z / (2*gridspacing));
    BigReal mlower = floor(sgmin.z / (2*gridspacing));
    sgmax.z = 2*gridspacing*mupper;
    sgmin.z = 2*gridspacing*mlower;
    zlen = sgmax.z - sgmin.z;
  }
#ifdef DEBUG_MSM_VERBOSE
  printf("zlen = %g   sgmin.z = %g   sgmax.z = %g\n", zlen, sgmin.z, sgmax.z);
#endif
  sglower = sgmin;

  int ia, ib, ja, jb, ka, kb;
  setup_hgrid_1d(xlen, hxlen, nhx, ia, ib, ispx);
  setup_hgrid_1d(ylen, hylen, nhy, ja, jb, ispy);
  setup_hgrid_1d(zlen, hzlen, nhz, ka, kb, ispz);
  hxlen_1 = 1 / hxlen;
  hylen_1 = 1 / hylen;
  hzlen_1 = 1 / hzlen;
  if (CkMyPe() == 0) {
    if (ispx || ispy || ispz) {
      iout << iINFO << "MSM grid spacing along X is "<< hxlen << " A\n";
      iout << iINFO << "MSM grid spacing along Y is "<< hylen << " A\n";
      iout << iINFO << "MSM grid spacing along Z is "<< hzlen << " A\n";
    }
    else {
      iout << iINFO << "MSM grid spacing is " << gridspacing << " A\n";
    }
    if ( ! ispx || ! ispy || ! ispz ) {
      iout << iINFO<<"MSM non-periodic padding is "<< padding << " A\n";
    }
  }

  int ni = ib - ia + 1;
  int nj = jb - ja + 1;
  int nk = kb - ka + 1;
  int n;

#if 0
  // reserve temp space for factored grid transfer operation
  n = (nk > omega ? nk : omega);  // row along z-dimension
  lzd.resize(n);
  n *= (nj > omega ? nj : omega);  // plane along yz-dimensions
  lyzd.resize(n);
#endif

  int lastnelems = 1;
  int smallestnbox = 1;
  int isclamped = 0;
  int maxlevels = nlevels;  // user-defined number of levels

#ifdef DEBUG_MSM_VERBOSE
  printf("maxlevels = %d\n", maxlevels);
#endif
  if (nlevels <= 0) {  // instead we set number of levels
    n = ni;
    if (n < nj) n = nj;
    if (n < nk) n = nk;
    for (maxlevels = 1;  n > 0;  n >>= 1)  maxlevels++;
    if (ispany == 0) {  // no periodicity
      // use rule of thumb 3/4 diameter of grid cutoff sphere
      int ngci = (int) ceil(3*a / hxlen) - 1;
      int ngcj = (int) ceil(3*a / hylen) - 1;
      int ngck = (int) ceil(3*a / hzlen) - 1;
      int omega3 = omega * omega * omega;
      int nhalf = (int) sqrt((double)ni * nj * nk);
      lastnelems = (nhalf > omega3 ? nhalf : omega3);
      smallestnbox = ngci * ngcj * ngck;  // smaller grids don't reduce work
      isclamped = 1;
    }
  }
#ifdef DEBUG_MSM_VERBOSE
  printf("maxlevels = %d\n", maxlevels);
#endif

  // allocate space for storing grid dimensions for each level
  map.gridrange.resize(maxlevels);

  // set periodicity flags
  map.ispx = ispx;
  map.ispy = ispy;
  map.ispz = ispz;

  int level = 0;
  int done = 0;
  int alldone = 0;
  do {
    map.gridrange[level].setbounds(ia, ib, ja, jb, ka, kb);

    // Msm index?

    if (++level == nlevels)  done |= 0x07;  // user limit on levels

    if (isclamped) {
      int nelems = ni * nj * nk;
      if (nelems <= lastnelems)    done |= 0x07;
      if (nelems <= smallestnbox)  done |= 0x07;
    }

    alldone = (done == 0x07);  // make sure all dimensions are done

    if (ispx) {
      ni >>= 1;
      ib = ni-1;
      if (ni & 1)        done |= 0x07;  // == 3 or 1
      else if (ni == 2)  done |= 0x01;  // can do one more
    }
    else {
      ia = -((-ia+1)/2) - s_edge;
      ib = (ib+1)/2 + s_edge;
      ni = ib - ia + 1;
      if (ni <= omega)   done |= 0x01;  // can do more restrictions
    }

    if (ispy) {
      nj >>= 1;
      jb = nj-1;
      if (nj & 1)        done |= 0x07;  // == 3 or 1
      else if (nj == 2)  done |= 0x02;  // can do one more
    }
    else {
      ja = -((-ja+1)/2) - s_edge;
      jb = (jb+1)/2 + s_edge;
      nj = jb - ja + 1;
      if (nj <= omega)   done |= 0x02;  // can do more restrictions
    }

    if (ispz) {
      nk >>= 1;
      kb = nk-1;
      if (nk & 1)        done |= 0x07;  // == 3 or 1
      else if (nk == 2)  done |= 0x04;  // can do one more
    }
    else {
      ka = -((-ka+1)/2) - s_edge;
      kb = (kb+1)/2 + s_edge;
      nk = kb - ka + 1;
      if (nk <= omega)   done |= 0x04;  // can do more restrictions
    }
  } while ( ! alldone );
  nlevels = level;

  // for periodic boundaries, don't visit top level (all 0)
  // toplevel visited only for all nonperiodic boundaries
  int toplevel = (ispany ? nlevels : nlevels - 1);

  // resize down to the actual number of levels (does not change alloc)
  map.gridrange.resize(nlevels);

  // print out some information about MSM
  if (CkMyPe() == 0) {
    iout << iINFO << "MSM using " << nlevels << " levels\n";
    for (n = 0;  n < nlevels;  n++) {
      char s[100];
      snprintf(s, sizeof(s), "    level %d:  "
          "[%d..%d] x [%d..%d] x [%d..%d]\n", n,
          map.gridrange[n].ia(), map.gridrange[n].ib(),
          map.gridrange[n].ja(), map.gridrange[n].jb(),
          map.gridrange[n].ka(), map.gridrange[n].kb());
      iout << iINFO << s;
    }
    iout << endi;
  }

  // find grid spacing basis vectors
  hu = hxlen * lattice.a().unit();
  hv = hylen * lattice.b().unit();
  hw = hzlen * lattice.c().unit();
  hufx = Float(hu.x);
  hufy = Float(hu.y);
  hufz = Float(hu.z);
  hvfx = Float(hv.x);
  hvfy = Float(hv.y);
  hvfz = Float(hv.z);
  hwfx = Float(hw.x);
  hwfy = Float(hw.y);
  hwfz = Float(hw.z);

  ru = lattice.a_r();
  rv = lattice.b_r();
  rw = lattice.c_r();

  // determine grid spacings in scaled space
  shx = ru * hu;
  shy = rv * hv;
  shz = rw * hw;
  shx_1 = 1 / shx;
  shy_1 = 1 / shy;
  shz_1 = 1 / shz;

  // row vectors to transform interpolated force back to real space
  // XXX Is not needed.
  sx_shx = shx_1 * Vector(ru.x, rv.x, rw.x);
  sy_shy = shy_1 * Vector(ru.y, rv.y, rw.y);
  sz_shz = shz_1 * Vector(ru.z, rv.z, rw.z);
  srx_x = Float(sx_shx.x);
  srx_y = Float(sx_shx.y);
  srx_z = Float(sx_shx.z);
  sry_x = Float(sy_shy.x);
  sry_y = Float(sy_shy.y);
  sry_z = Float(sy_shy.z);
  srz_x = Float(sz_shz.x);
  srz_y = Float(sz_shz.y);
  srz_z = Float(sz_shz.z);

  Vector pu = cross(hv, hw);
  BigReal s = (hu * pu) / (pu * pu);
  pu *= s;  // pu is orthogonal projection of hu onto hv CROSS hw

  Vector pv = cross(hw, hu);
  s = (hv * pv) / (pv * pv);
  pv *= s;  // pv is orthogonal projection of hv onto hw CROSS hu

  Vector pw = cross(hu, hv);
  s = (hw * pw) / (pw * pw);
  pw *= s;  // pw is orthogonal projection of hw onto hu CROSS hv

  // radii for parallelepiped of weights enclosing grid cutoff sphere
  ni = (int) ceil(2*a / pu.length() ) - 1;
  nj = (int) ceil(2*a / pv.length() ) - 1;
  nk = (int) ceil(2*a / pw.length() ) - 1;

  Float scaling = 1;
  Float scaling_factor = 0.5f;
  BigReal a_1 = 1/a;
  BigReal a_p = a_1;
  if (dispersion) {
    a_p = a_p * a_p * a_p;   // = 1/a^3
    a_p = a_p * a_p;         // = 1/a^6
    scaling_factor = 1.f/64;  // = 1/2^6
  }
  int i, j, k;
  if (approx < C1HERMITE) {
    // resize gc and gvc constants for number of levels
    map.gc.resize(nlevels);
    map.gvc.resize(nlevels);

    for (level = 0;  level < toplevel;  level++) {
      map.gc[level].setbounds(-ni, ni, -nj, nj, -nk, nk);
      map.gvc[level].setbounds(-ni, ni, -nj, nj, -nk, nk);

      for (k = -nk;  k <= nk;  k++) {
        for (j = -nj;  j <= nj;  j++) {
          for (i = -ni;  i <= ni;  i++) {
            if (level == 0) {
              BigReal s, t, gs=0, gt=0, g=0, dgs=0, dgt=0, dg=0;
              BigReal vlen = (i*hu + j*hv + k*hw).length();
              s = vlen * a_1;
              t = 0.5 * s;
              if (t >= 1) {
                g = 0;
                dg = 0;
              }
              else {
                splitting(gt, dgt, t, split);
                if (s >= 1) {
                  BigReal s_1 = 1/s;
                  if (dispersion) {
                    gs = s_1 * s_1 * s_1;  // = 1/s^3
                    gs = gs * gs;  // = 1/s^6
                    dgs = -6 * gs * s_1;
                  }
                  else {
                    gs = s_1;
                    dgs = -gs * s_1;
                  }
                }
                else {
                  splitting(gs, dgs, s, split);
                }
                g = (gs - scaling_factor * gt) * a_p;
                BigReal c=0;
                if (i || j || k) {
                  c = a_p * a_1 / vlen;
                }
                dg = 0.5 * (dgs - 0.5*scaling_factor * dgt) * c;

                // Msm index?

              }
              map.gc[0](i,j,k) = Float(g);
              map.gvc[0](i,j,k) = Float(dg);
            } // if level 0
            else {
              map.gc[level](i,j,k) = scaling * map.gc[0](i,j,k);
              map.gvc[level](i,j,k) = scaling * map.gvc[0](i,j,k);
            }

          } // for i
        } // for j
      } // for k
      scaling *= scaling_factor;

    } // for level

    // for summed virial factors
    gvsum.setbounds(-ni, ni, -nj, nj, -nk, nk);
    // make sure final virial sum is initialized to 0
    for (i = 0;  i < VMAX;  i++) { virial[i] = 0; }

    if (toplevel < nlevels) {
      // nonperiodic along all basis vector directions
      // calculate top level weights where all grid points
      // interact with each other
      ni = map.gridrange[toplevel].ni();
      nj = map.gridrange[toplevel].nj();
      nk = map.gridrange[toplevel].nk();
      map.gc[toplevel].setbounds(-ni, ni, -nj, nj, -nk, nk);

      // Msm index?

      for (k = -nk;  k <= nk;  k++) {
        for (j = -nj;  j <= nj;  j++) {
          for (i = -ni;  i <= ni;  i++) {
            BigReal s, gs, d;
            BigReal vlen = (i*hu + j*hv + k*hw).length();
            s = vlen * a_1;
            if (s >= 1) {
              BigReal s_1 = 1/s;
              if (dispersion) {
                gs = s_1 * s_1 * s_1;  // = 1/s^3
                gs = gs * gs;  // = 1/s^6
              }
              else {
                gs = s_1;
              }
            }
            else {
              splitting(gs, d, s, split);
            }
            map.gc[toplevel](i,j,k) = scaling * Float(gs * a_p);
          } // for i
        } // for j
      } // for k
    } // if toplevel

    // generate grespro stencil
    const int nstencil = Nstencil[approx];
    const Float *phi = PhiStencil[approx];
    map.grespro.set(0, nstencil, 0, nstencil, 0, nstencil);
    for (k = 0;  k < nstencil;  k++) {
      for (j = 0;  j < nstencil;  j++) {
        for (i = 0;  i < nstencil;  i++) {
          map.grespro(i,j,k) = phi[i] * phi[j] * phi[k];
        }
      }
    }

  } // end if approx < C1HERMITE
  else {
    // C1HERMITE
    // resize gc_c1hermite constants for number of levels
    map.gc_c1hermite.resize(nlevels);
    scaling = 1;

    for (level = 0;  level < toplevel;  level++) {

      Vector hmu = scaling * hu;
      Vector hmv = scaling * hv;
      Vector hmw = scaling * hw;
      BigReal am = scaling * a;

      map.gc_c1hermite[level].setbounds(-ni, ni, -nj, nj, -nk, nk);

      for (k = -nk;  k <= nk;  k++) {
        for (j = -nj;  j <= nj;  j++) {
          for (i = -ni;  i <= ni;  i++) {
            C1Matrix& m = map.gc_c1hermite[level](i,j,k);
            Vector rv = i*hmu + j*hmv + k*hmw;
            BigReal r2 = rv * rv;
            m.set(0);
            if (r2 < 4*am*am) {
              // accumulate D( g_{a}(0,r) ) term for this level
              gc_c1hermite_elem_accum(m, 1, rv, am, split);
              // accumulate D( -g_{2a}(0,r) ) term for this level
              gc_c1hermite_elem_accum(m, -1, rv, 2*am, split);
            } // if within cutoff
          }
        }
      } // end loop over gc_c1hermite elements for this level
      scaling *= 2;  // double grid spacing and cutoff at each iteration

    } // end loop over levels

    if (toplevel < nlevels) {
      Vector hmu = scaling * hu;
      Vector hmv = scaling * hv;
      Vector hmw = scaling * hw;
      BigReal am = scaling * a;

      // nonperiodic along all basis vector directions
      // calculate top level weights where all grid points
      // interact with each other
      ni = map.gridrange[toplevel].ni();
      nj = map.gridrange[toplevel].nj();
      nk = map.gridrange[toplevel].nk();
      map.gc_c1hermite[toplevel].setbounds(-ni, ni, -nj, nj, -nk, nk);

      for (k = -nk;  k <= nk;  k++) {
        for (j = -nj;  j <= nj;  j++) {
          for (i = -ni;  i <= ni;  i++) {
            C1Matrix& m = map.gc_c1hermite[level](i,j,k);
            Vector rv = i*hmu + j*hmv + k*hmw;
            m.set(0);
            // accumulate D( g_{a}(0,r) ) term for this level
            gc_c1hermite_elem_accum(m, 1, rv, am, split);
          }
        }
      } // end loop over gc_c1hermite elements for top level

    } // end if top level

    // C1 Hermite restriction and prolongation stencils
    map.gres_c1hermite.resize(nlevels-1);
    map.gpro_c1hermite.resize(nlevels-1);

    enum {
      ND = 3,    // stencil diameter
      NR = ND/2  // stencil radius
    };

    // the master basis functions PHI0 and PHI1 for the 3-point stencil
    // and their derivatives DPHI0 and DPHI1
    const double  PHI0[ND] = { 0.5, 1, 0.5 };
    const double DPHI0[ND] = { 1.5, 0, -1.5 };
    const double  PHI1[ND] = { -0.125, 0, 0.125 };
    const double DPHI1[ND] = { -0.25, 1, -0.25 };

    // for intermediate calculations
    double  xphi0_base_array[ND];
    double dxphi0_base_array[ND];
    double  yphi0_base_array[ND];
    double dyphi0_base_array[ND];
    double  zphi0_base_array[ND];
    double dzphi0_base_array[ND];
    double  xphi1_base_array[ND];
    double dxphi1_base_array[ND];
    double  yphi1_base_array[ND];
    double dyphi1_base_array[ND];
    double  zphi1_base_array[ND];
    double dzphi1_base_array[ND];
    // will point to center of stencil arrays
    double *xphi0, *dxphi0, *xphi1, *dxphi1;
    double *yphi0, *dyphi0, *yphi1, *dyphi1;
    double *zphi0, *dzphi0, *zphi1, *dzphi1;

    for (n = 0;  n < ND;  n++) {
      xphi0_base_array[n]  = PHI0[n];
      dxphi0_base_array[n] = hxlen_1 * DPHI0[n];  // scale by grid spacing
      xphi1_base_array[n]  = hxlen * PHI1[n];     // scale by grid spacing
      dxphi1_base_array[n] = DPHI1[n];
      yphi0_base_array[n]  = PHI0[n];
      dyphi0_base_array[n] = hylen_1 * DPHI0[n];  // scale by grid spacing
      yphi1_base_array[n]  = hylen * PHI1[n];     // scale by grid spacing
      dyphi1_base_array[n] = DPHI1[n];
      zphi0_base_array[n]  = PHI0[n];
      dzphi0_base_array[n] = hzlen_1 * DPHI0[n];  // scale by grid spacing
      zphi1_base_array[n]  = hzlen * PHI1[n];     // scale by grid spacing
      dzphi1_base_array[n] = DPHI1[n];
    }
    xphi0  =  xphi0_base_array + NR;  // point into center of arrays
    dxphi0 = dxphi0_base_array + NR;
    xphi1  =  xphi1_base_array + NR;
    dxphi1 = dxphi1_base_array + NR;
    yphi0  =  yphi0_base_array + NR;
    dyphi0 = dyphi0_base_array + NR;
    yphi1  =  yphi1_base_array + NR;
    dyphi1 = dyphi1_base_array + NR;
    zphi0  =  zphi0_base_array + NR;
    dzphi0 = dzphi0_base_array + NR;
    zphi1  =  zphi1_base_array + NR;
    dzphi1 = dzphi1_base_array + NR;

    for (level = 0;  level < nlevels-1;  level++) {
      // allocate space for restriction and prolongation stencils
      map.gres_c1hermite[level].set(0, ND, 0, ND, 0, ND);
      map.gpro_c1hermite[level].set(0, ND, 0, ND, 0, ND);

      // scale up to next level grid spacing
      //
      // have to determine for each dimension whether or not 
      // a periodic grid spacing has increased 
      // (equivalent to if there are fewer grid points)
      for (n = -NR;  n <= NR;  n++) {
        if ( ! ispx ||
              map.gridrange[level+1].ni() < map.gridrange[level].ni() ) {
          dxphi0[n] *= 0.5;
          xphi1[n] *= 2;
        }
        if ( ! ispy ||
              map.gridrange[level+1].nj() < map.gridrange[level].nj() ) {
          dyphi0[n] *= 0.5;
          yphi1[n] *= 2;
        }
        if ( ! ispz ||
              map.gridrange[level+1].nk() < map.gridrange[level].nk() ) {
          dzphi0[n] *= 0.5;
          zphi1[n] *= 2;
        }
      }

      // loop over restriction stencil matrices
      // calculate from partial derivatives
      for (k = -NR;  k <= NR;  k++) {
        for (j = -NR;  j <= NR;  j++) {
          for (i = -NR;  i <= NR;  i++) {
            Float *t = map.gres_c1hermite[level](i+NR,j+NR,k+NR).melem;

            t[C1INDEX(D000,D000)] =  xphi0[i] *  yphi0[j]  *  zphi0[k];
            t[C1INDEX(D000,D100)] = dxphi0[i] *  yphi0[j]  *  zphi0[k];
            t[C1INDEX(D000,D010)] =  xphi0[i] * dyphi0[j]  *  zphi0[k];
            t[C1INDEX(D000,D001)] =  xphi0[i] *  yphi0[j]  * dzphi0[k];
            t[C1INDEX(D000,D110)] = dxphi0[i] * dyphi0[j]  *  zphi0[k];
            t[C1INDEX(D000,D101)] = dxphi0[i] *  yphi0[j]  * dzphi0[k];
            t[C1INDEX(D000,D011)] =  xphi0[i] * dyphi0[j]  * dzphi0[k];
            t[C1INDEX(D000,D111)] = dxphi0[i] * dyphi0[j]  * dzphi0[k];

            t[C1INDEX(D100,D000)] =  xphi1[i] *  yphi0[j]  *  zphi0[k];
            t[C1INDEX(D100,D100)] = dxphi1[i] *  yphi0[j]  *  zphi0[k];
            t[C1INDEX(D100,D010)] =  xphi1[i] * dyphi0[j]  *  zphi0[k];
            t[C1INDEX(D100,D001)] =  xphi1[i] *  yphi0[j]  * dzphi0[k];
            t[C1INDEX(D100,D110)] = dxphi1[i] * dyphi0[j]  *  zphi0[k];
            t[C1INDEX(D100,D101)] = dxphi1[i] *  yphi0[j]  * dzphi0[k];
            t[C1INDEX(D100,D011)] =  xphi1[i] * dyphi0[j]  * dzphi0[k];
            t[C1INDEX(D100,D111)] = dxphi1[i] * dyphi0[j]  * dzphi0[k];

            t[C1INDEX(D010,D000)] =  xphi0[i] *  yphi1[j]  *  zphi0[k];
            t[C1INDEX(D010,D100)] = dxphi0[i] *  yphi1[j]  *  zphi0[k];
            t[C1INDEX(D010,D010)] =  xphi0[i] * dyphi1[j]  *  zphi0[k];
            t[C1INDEX(D010,D001)] =  xphi0[i] *  yphi1[j]  * dzphi0[k];
            t[C1INDEX(D010,D110)] = dxphi0[i] * dyphi1[j]  *  zphi0[k];
            t[C1INDEX(D010,D101)] = dxphi0[i] *  yphi1[j]  * dzphi0[k];
            t[C1INDEX(D010,D011)] =  xphi0[i] * dyphi1[j]  * dzphi0[k];
            t[C1INDEX(D010,D111)] = dxphi0[i] * dyphi1[j]  * dzphi0[k];

            t[C1INDEX(D001,D000)] =  xphi0[i] *  yphi0[j]  *  zphi1[k];
            t[C1INDEX(D001,D100)] = dxphi0[i] *  yphi0[j]  *  zphi1[k];
            t[C1INDEX(D001,D010)] =  xphi0[i] * dyphi0[j]  *  zphi1[k];
            t[C1INDEX(D001,D001)] =  xphi0[i] *  yphi0[j]  * dzphi1[k];
            t[C1INDEX(D001,D110)] = dxphi0[i] * dyphi0[j]  *  zphi1[k];
            t[C1INDEX(D001,D101)] = dxphi0[i] *  yphi0[j]  * dzphi1[k];
            t[C1INDEX(D001,D011)] =  xphi0[i] * dyphi0[j]  * dzphi1[k];
            t[C1INDEX(D001,D111)] = dxphi0[i] * dyphi0[j]  * dzphi1[k];

            t[C1INDEX(D110,D000)] =  xphi1[i] *  yphi1[j]  *  zphi0[k];
            t[C1INDEX(D110,D100)] = dxphi1[i] *  yphi1[j]  *  zphi0[k];
            t[C1INDEX(D110,D010)] =  xphi1[i] * dyphi1[j]  *  zphi0[k];
            t[C1INDEX(D110,D001)] =  xphi1[i] *  yphi1[j]  * dzphi0[k];
            t[C1INDEX(D110,D110)] = dxphi1[i] * dyphi1[j]  *  zphi0[k];
            t[C1INDEX(D110,D101)] = dxphi1[i] *  yphi1[j]  * dzphi0[k];
            t[C1INDEX(D110,D011)] =  xphi1[i] * dyphi1[j]  * dzphi0[k];
            t[C1INDEX(D110,D111)] = dxphi1[i] * dyphi1[j]  * dzphi0[k];

            t[C1INDEX(D101,D000)] =  xphi1[i] *  yphi0[j]  *  zphi1[k];
            t[C1INDEX(D101,D100)] = dxphi1[i] *  yphi0[j]  *  zphi1[k];
            t[C1INDEX(D101,D010)] =  xphi1[i] * dyphi0[j]  *  zphi1[k];
            t[C1INDEX(D101,D001)] =  xphi1[i] *  yphi0[j]  * dzphi1[k];
            t[C1INDEX(D101,D110)] = dxphi1[i] * dyphi0[j]  *  zphi1[k];
            t[C1INDEX(D101,D101)] = dxphi1[i] *  yphi0[j]  * dzphi1[k];
            t[C1INDEX(D101,D011)] =  xphi1[i] * dyphi0[j]  * dzphi1[k];
            t[C1INDEX(D101,D111)] = dxphi1[i] * dyphi0[j]  * dzphi1[k];

            t[C1INDEX(D011,D000)] =  xphi0[i] *  yphi1[j]  *  zphi1[k];
            t[C1INDEX(D011,D100)] = dxphi0[i] *  yphi1[j]  *  zphi1[k];
            t[C1INDEX(D011,D010)] =  xphi0[i] * dyphi1[j]  *  zphi1[k];
            t[C1INDEX(D011,D001)] =  xphi0[i] *  yphi1[j]  * dzphi1[k];
            t[C1INDEX(D011,D110)] = dxphi0[i] * dyphi1[j]  *  zphi1[k];
            t[C1INDEX(D011,D101)] = dxphi0[i] *  yphi1[j]  * dzphi1[k];
            t[C1INDEX(D011,D011)] =  xphi0[i] * dyphi1[j]  * dzphi1[k];
            t[C1INDEX(D011,D111)] = dxphi0[i] * dyphi1[j]  * dzphi1[k];

            t[C1INDEX(D111,D000)] =  xphi1[i] *  yphi1[j]  *  zphi1[k];
            t[C1INDEX(D111,D100)] = dxphi1[i] *  yphi1[j]  *  zphi1[k];
            t[C1INDEX(D111,D010)] =  xphi1[i] * dyphi1[j]  *  zphi1[k];
            t[C1INDEX(D111,D001)] =  xphi1[i] *  yphi1[j]  * dzphi1[k];
            t[C1INDEX(D111,D110)] = dxphi1[i] * dyphi1[j]  *  zphi1[k];
            t[C1INDEX(D111,D101)] = dxphi1[i] *  yphi1[j]  * dzphi1[k];
            t[C1INDEX(D111,D011)] =  xphi1[i] * dyphi1[j]  * dzphi1[k];
            t[C1INDEX(D111,D111)] = dxphi1[i] * dyphi1[j]  * dzphi1[k];
          }
        }
      } // end loops over restriction stencil matrices

      // loop over prolongation stencil matrices
      // prolongation stencil matrices are the transpose of restriction
      for (k = -NR;  k <= NR;  k++) {
        for (j = -NR;  j <= NR;  j++) {
          for (i = -NR;  i <= NR;  i++) {
            Float *t = map.gres_c1hermite[level](i+NR,j+NR,k+NR).melem;
            Float *tt = map.gpro_c1hermite[level](i+NR,j+NR,k+NR).melem;

            tt[C1INDEX(D000,D000)] = t[C1INDEX(D000,D000)];
            tt[C1INDEX(D000,D100)] = t[C1INDEX(D100,D000)];
            tt[C1INDEX(D000,D010)] = t[C1INDEX(D010,D000)];
            tt[C1INDEX(D000,D001)] = t[C1INDEX(D001,D000)];
            tt[C1INDEX(D000,D110)] = t[C1INDEX(D110,D000)];
            tt[C1INDEX(D000,D101)] = t[C1INDEX(D101,D000)];
            tt[C1INDEX(D000,D011)] = t[C1INDEX(D011,D000)];
            tt[C1INDEX(D000,D111)] = t[C1INDEX(D111,D000)];

            tt[C1INDEX(D100,D000)] = t[C1INDEX(D000,D100)];
            tt[C1INDEX(D100,D100)] = t[C1INDEX(D100,D100)];
            tt[C1INDEX(D100,D010)] = t[C1INDEX(D010,D100)];
            tt[C1INDEX(D100,D001)] = t[C1INDEX(D001,D100)];
            tt[C1INDEX(D100,D110)] = t[C1INDEX(D110,D100)];
            tt[C1INDEX(D100,D101)] = t[C1INDEX(D101,D100)];
            tt[C1INDEX(D100,D011)] = t[C1INDEX(D011,D100)];
            tt[C1INDEX(D100,D111)] = t[C1INDEX(D111,D100)];

            tt[C1INDEX(D010,D000)] = t[C1INDEX(D000,D010)];
            tt[C1INDEX(D010,D100)] = t[C1INDEX(D100,D010)];
            tt[C1INDEX(D010,D010)] = t[C1INDEX(D010,D010)];
            tt[C1INDEX(D010,D001)] = t[C1INDEX(D001,D010)];
            tt[C1INDEX(D010,D110)] = t[C1INDEX(D110,D010)];
            tt[C1INDEX(D010,D101)] = t[C1INDEX(D101,D010)];
            tt[C1INDEX(D010,D011)] = t[C1INDEX(D011,D010)];
            tt[C1INDEX(D010,D111)] = t[C1INDEX(D111,D010)];

            tt[C1INDEX(D001,D000)] = t[C1INDEX(D000,D001)];
            tt[C1INDEX(D001,D100)] = t[C1INDEX(D100,D001)];
            tt[C1INDEX(D001,D010)] = t[C1INDEX(D010,D001)];
            tt[C1INDEX(D001,D001)] = t[C1INDEX(D001,D001)];
            tt[C1INDEX(D001,D110)] = t[C1INDEX(D110,D001)];
            tt[C1INDEX(D001,D101)] = t[C1INDEX(D101,D001)];
            tt[C1INDEX(D001,D011)] = t[C1INDEX(D011,D001)];
            tt[C1INDEX(D001,D111)] = t[C1INDEX(D111,D001)];

            tt[C1INDEX(D110,D000)] = t[C1INDEX(D000,D110)];
            tt[C1INDEX(D110,D100)] = t[C1INDEX(D100,D110)];
            tt[C1INDEX(D110,D010)] = t[C1INDEX(D010,D110)];
            tt[C1INDEX(D110,D001)] = t[C1INDEX(D001,D110)];
            tt[C1INDEX(D110,D110)] = t[C1INDEX(D110,D110)];
            tt[C1INDEX(D110,D101)] = t[C1INDEX(D101,D110)];
            tt[C1INDEX(D110,D011)] = t[C1INDEX(D011,D110)];
            tt[C1INDEX(D110,D111)] = t[C1INDEX(D111,D110)];

            tt[C1INDEX(D101,D000)] = t[C1INDEX(D000,D101)];
            tt[C1INDEX(D101,D100)] = t[C1INDEX(D100,D101)];
            tt[C1INDEX(D101,D010)] = t[C1INDEX(D010,D101)];
            tt[C1INDEX(D101,D001)] = t[C1INDEX(D001,D101)];
            tt[C1INDEX(D101,D110)] = t[C1INDEX(D110,D101)];
            tt[C1INDEX(D101,D101)] = t[C1INDEX(D101,D101)];
            tt[C1INDEX(D101,D011)] = t[C1INDEX(D011,D101)];
            tt[C1INDEX(D101,D111)] = t[C1INDEX(D111,D101)];

            tt[C1INDEX(D011,D000)] = t[C1INDEX(D000,D011)];
            tt[C1INDEX(D011,D100)] = t[C1INDEX(D100,D011)];
            tt[C1INDEX(D011,D010)] = t[C1INDEX(D010,D011)];
            tt[C1INDEX(D011,D001)] = t[C1INDEX(D001,D011)];
            tt[C1INDEX(D011,D110)] = t[C1INDEX(D110,D011)];
            tt[C1INDEX(D011,D101)] = t[C1INDEX(D101,D011)];
            tt[C1INDEX(D011,D011)] = t[C1INDEX(D011,D011)];
            tt[C1INDEX(D011,D111)] = t[C1INDEX(D111,D011)];

            tt[C1INDEX(D111,D000)] = t[C1INDEX(D000,D111)];
            tt[C1INDEX(D111,D100)] = t[C1INDEX(D100,D111)];
            tt[C1INDEX(D111,D010)] = t[C1INDEX(D010,D111)];
            tt[C1INDEX(D111,D001)] = t[C1INDEX(D001,D111)];
            tt[C1INDEX(D111,D110)] = t[C1INDEX(D110,D111)];
            tt[C1INDEX(D111,D101)] = t[C1INDEX(D101,D111)];
            tt[C1INDEX(D111,D011)] = t[C1INDEX(D011,D111)];
            tt[C1INDEX(D111,D111)] = t[C1INDEX(D111,D111)];
          }
        }
      } // end loops over prolongation stencil matrices

    } // end loop over levels

  } // end if C1HERMITE

  // calculate self energy factor for splitting
  BigReal gs=0, d=0;
  splitting(gs, d, 0, split);
  gzero = gs * a_p;

  if (CkMyPe() == 0) {
    iout << iINFO << "MSM finished calculating stencils\n" << endi;
  }

  // allocate map for patches
  PatchMap *pm = PatchMap::Object();
  int numpatches = pm->numPatches();
  map.patchList.resize(numpatches);
#ifdef DEBUG_MSM_VERBOSE
  printf("numPatches = %d\n", numpatches);
#endif

  // allocate map for blocks for each grid level
  map.blockLevel.resize(nlevels);
  map.bsx.resize(nlevels);
  map.bsy.resize(nlevels);
  map.bsz.resize(nlevels);
#ifdef MSM_FOLD_FACTOR
  map.foldfactor.resize(nlevels);
#endif
  for (level = 0;  level < nlevels;  level++) {
    msm::IndexRange& g = map.gridrange[level];
    msm::Grid<msm::BlockDiagram>& b = map.blockLevel[level];
    int gia = g.ia();
    int gni = g.ni();
    int gja = g.ja();
    int gnj = g.nj();
    int gka = g.ka();
    int gnk = g.nk();
    map.bsx[level] = bsx;
    map.bsy[level] = bsy;
    map.bsz[level] = bsz;
    if (/* map.bsx[level] < gni ||
        map.bsy[level] < gnj ||
        map.bsz[level] < gnk */ 1) {
      // make sure that block sizes divide evenly into periodic dimensions
      if (ispx) setup_periodic_blocksize(map.bsx[level], gni);
      if (ispy) setup_periodic_blocksize(map.bsy[level], gnj);
      if (ispz) setup_periodic_blocksize(map.bsz[level], gnk);
#ifdef MSM_DEBUG_VERBOSE
      if (CkMyPe() == 0) {
        printf("level = %d\n  map.bs* = %d %d %d  gn* = %d %d %d\n",
            level, map.bsx[level], map.bsy[level], map.bsz[level],gni,gnj,gnk);
      }
#endif
      // subdivide grid into multiple blocks
      //   == ceil(gni / bsx), etc.
      int bni = (gni / map.bsx[level]) + (gni % map.bsx[level] != 0);
      int bnj = (gnj / map.bsy[level]) + (gnj % map.bsy[level] != 0);
      int bnk = (gnk / map.bsz[level]) + (gnk % map.bsz[level] != 0);
#ifdef MSM_FOLD_FACTOR
      if (/* level > 2 && */ (bni == 1 || bnj == 1 || bnk == 1)) {
        map.foldfactor[level].set(bsx / gni, bsy / gnj, bsz / gnk);
#if 0
        if (CkMyPe() == 0) {
          printf("Setting MSM FoldFactor level %d:  %d %d %d\n",
              level, bsx / gni, bsy / gnj, bsz / gnk);
        }
#endif
      }
#endif
      b.set(0, bni, 0, bnj, 0, bnk);
      for (k = 0;  k < bnk;  k++) {
        for (j = 0;  j < bnj;  j++) {
          for (i = 0;  i < bni;  i++) {
            b(i,j,k).reset();
            int ia = gia + i*map.bsx[level];
            int ib = ia + map.bsx[level] - 1;
            int ja = gja + j*map.bsy[level];
            int jb = ja + map.bsy[level] - 1;
            int ka = gka + k*map.bsz[level];
            int kb = ka + map.bsz[level] - 1;
            if (ib >= gia + gni) ib = gia + gni - 1;
            if (jb >= gja + gnj) jb = gja + gnj - 1;
            if (kb >= gka + gnk) kb = gka + gnk - 1;
            b(i,j,k).nrange.setbounds(ia, ib, ja, jb, ka, kb);
          }
        }
      }
    }
    /*
    else {
      // make entire grid into single block
      b.set(0, 1, 0, 1, 0, 1);
      b(0,0,0).reset();
      b(0,0,0).nrange.set(gia, gni, gja, gnj, gka, gnk);
      // set correct block dimensions
      map.bsx[level] = gni;
      map.bsy[level] = gnj;
      map.bsz[level] = gnk;
    }
    */
  }
  //CkExit();
#ifdef DEBUG_MSM_VERBOSE
  printf("Done allocating map for grid levels\n");
  printf("Grid level decomposition:\n");
  for (level = 0;  level < nlevels;  level++) {
    msm::Grid<msm::BlockDiagram>& b = map.blockLevel[level];
    int bia = b.ia();
    int bib = b.ib();
    int bja = b.ja();
    int bjb = b.jb();
    int bka = b.ka();
    int bkb = b.kb();
    for (k = bka;  k <= bkb;  k++) {
      for (j = bja;  j <= bjb;  j++) {
        for (i = bia;  i <= bib;  i++) {
          int ia = b(i,j,k).nrange.ia();
          int ib = b(i,j,k).nrange.ib();
          int ja = b(i,j,k).nrange.ja();
          int jb = b(i,j,k).nrange.jb();
          int ka = b(i,j,k).nrange.ka();
          int kb = b(i,j,k).nrange.kb();
          printf("level=%d  id=%d %d %d  [%d..%d] x [%d..%d] x [%d..%d]"
              " --> %d\n",
              level, i, j, k, ia, ib, ja, jb, ka, kb,
              b(i,j,k).nrange.nn());
        }
      }
    }
  }
#endif
  if (CkMyPe() == 0) {
    iout << iINFO << "MSM finished creating map for grid levels\n" << endi;
  }

  initialize2();
}

void ComputeMsmMgr::initialize2()
{
  SimParameters *simParams = Node::Object()->simParameters;
  PatchMap *pm = PatchMap::Object();
  int numpatches = pm->numPatches();
  int i, j, k, n, level;

  // initialize grid of PatchDiagram
  // a = cutoff
  BigReal sysdima = lattice.a_r().unit() * lattice.a();
  BigReal sysdimb = lattice.b_r().unit() * lattice.b();
  BigReal sysdimc = lattice.c_r().unit() * lattice.c();
  BigReal patchdim = simParams->patchDimension;
  BigReal xmargin = 0.5 * (patchdim - a) / sysdima;
  BigReal ymargin = 0.5 * (patchdim - a) / sysdimb;
  BigReal zmargin = 0.5 * (patchdim - a) / sysdimc;
#if 0
  // set min and max grid indices for patch covering
  // for non-periodic boundaries they conform to grid
  // periodic permits wrapping, so set to min/max for int
  int ia_min = (lattice.a_p() ? INT_MIN : map.gridrange[0].ia());
  int ib_max = (lattice.a_p() ? INT_MAX : map.gridrange[0].ib());
  int ja_min = (lattice.b_p() ? INT_MIN : map.gridrange[0].ja());
  int jb_max = (lattice.b_p() ? INT_MAX : map.gridrange[0].jb());
  int ka_min = (lattice.c_p() ? INT_MIN : map.gridrange[0].ka());
  int kb_max = (lattice.c_p() ? INT_MAX : map.gridrange[0].kb());
#endif
  int pid;
  for (pid = 0;  pid < numpatches;  pid++) {
    // shortcut reference to this patch diagram
    msm::PatchDiagram& p = map.patchList[pid];
    p.reset();
    // find extent of patch including margin
    BigReal xmin = pm->min_a(pid) - xmargin;
    BigReal xmax = pm->max_a(pid) + xmargin;
    BigReal ymin = pm->min_b(pid) - ymargin;
    BigReal ymax = pm->max_b(pid) + ymargin;
    BigReal zmin = pm->min_c(pid) - zmargin;
    BigReal zmax = pm->max_c(pid) + zmargin;
    // find grid point covering of patch plus outer edge stencil
    int ia = int(floor((xmin - sglower.x) * shx_1)) - s_edge;
    int ib = int(floor((xmax - sglower.x) * shx_1)) + 1 + s_edge;
    int ja = int(floor((ymin - sglower.y) * shy_1)) - s_edge;
    int jb = int(floor((ymax - sglower.y) * shy_1)) + 1 + s_edge;
    int ka = int(floor((zmin - sglower.z) * shz_1)) - s_edge;
    int kb = int(floor((zmax - sglower.z) * shz_1)) + 1 + s_edge;
    // for edge patches along non-periodic boundaries
    // clamp subgrid to full grid boundaries
    if ( ! lattice.a_p() ) {  // non-periodic along lattice basis vector a
      int mi = pm->index_a(pid);
      if (mi == 0)                   ia = map.gridrange[0].ia();
      if (mi == pm->gridsize_a()-1)  ib = map.gridrange[0].ib();
    }
    if ( ! lattice.b_p() ) {  // non-periodic along lattice basis vector b
      int mj = pm->index_b(pid);
      if (mj == 0)                   ja = map.gridrange[0].ja();
      if (mj == pm->gridsize_b()-1)  jb = map.gridrange[0].jb();
    }
    if ( ! lattice.c_p() ) {  // non-periodic along lattice basis vector a
      int mk = pm->index_c(pid);
      if (mk == 0)                   ka = map.gridrange[0].ka();
      if (mk == pm->gridsize_c()-1)  kb = map.gridrange[0].kb();
    }
#if 0
    // truncate subgrid covering to grid min/max
    // so that subgrid does not extend beyond full grid
    // works for both periodic and non-periodic boundary conditions
    if (ia < ia_min)  ia = ia_min;
    if (ib > ib_max)  ib = ib_max;
    if (ja < ja_min)  ja = ja_min;
    if (jb > jb_max)  jb = jb_max;
    if (ka < ka_min)  ka = ka_min;
    if (kb > kb_max)  kb = kb_max;
    // check for edge patch and extend subgrid to grid min/max
    // so that subgrid fully covers up to the edge of full grid
    int mi = pm->index_a(pid);
    int mj = pm->index_b(pid);
    int mk = pm->index_c(pid);
    int npi = pm->gridsize_a();
    int npj = pm->gridsize_b();
    int npk = pm->gridsize_c();
    if (mi == 0)      ia = ia_min;
    if (mi == npi-1)  ib = ib_max;
    if (mj == 0)      ja = ja_min;
    if (mj == npj-1)  jb = jb_max;
    if (mk == 0)      ka = ka_min;
    if (mk == npk-1)  kb = kb_max;
#endif
#if 0
    printf("patch %d:  grid covering:  [%d..%d] x [%d..%d] x [%d..%d]\n",
        pid, ia, ib, ja, jb, ka, kb);
    fflush(stdout);
#endif
    // set the index range for this patch's surrounding grid points
    p.nrange.setbounds(ia,ib,ja,jb,ka,kb);
    // find lower and upper blocks of MSM h-grid
    msm::BlockIndex blower = map.blockOfGridIndex(msm::Ivec(ia,ja,ka),0);
    msm::BlockIndex bupper = map.blockOfGridIndex(msm::Ivec(ib,jb,kb),0);
    int maxarrlen = (bupper.n.i - blower.n.i + 1) *
      (bupper.n.j - blower.n.j + 1) * (bupper.n.k - blower.n.k + 1);
    p.send.setmax(maxarrlen);  // allocate space for send array
    // loop over the blocks
#if 0
    printf("blower: level=%d  n=%d %d %d   bupper: level=%d  n=%d %d %d\n",
        blower.level, blower.n.i, blower.n.j, blower.n.k,
        bupper.level, bupper.n.i, bupper.n.j, bupper.n.k);
    fflush(stdout);
#endif
    for (int kk = blower.n.k;  kk <= bupper.n.k;  kk++) {
      for (int jj = blower.n.j;  jj <= bupper.n.j;  jj++) {
        for (int ii = blower.n.i;  ii <= bupper.n.i;  ii++) {
#if 0
          printf("ii=%d  jj=%d  kk=%d\n", ii, jj, kk);
          fflush(stdout);
#endif
          // determine actual block and range to send to
          msm::BlockSend bs;
          bs.nblock.n = msm::Ivec(ii,jj,kk);
          bs.nblock.level = 0;
          bs.nrange = map.clipBlockToIndexRange(bs.nblock, p.nrange);
          map.wrapBlockSend(bs);  // determine wrapping to true block index
          p.send.append(bs);  // append this block to the send array
          // increment counter for receive block
          map.blockLevel[0](bs.nblock_wrap.n).numRecvsCharge++;
          // initialize patch send back from this block
          msm::PatchSend ps;
          ps.nrange = bs.nrange_wrap;
          ps.nrange_unwrap = bs.nrange;
          ps.patchID = pid;
          map.blockLevel[0](bs.nblock_wrap.n).sendPatch.append(ps);
          // increment number of receives back to this patch
          p.numRecvs++;
        }
      }
    }
    // number of receives should be same as number of sends
    ASSERT(p.numRecvs == p.send.len() );
  }
#ifdef DEBUG_MSM_VERBOSE
if (CkMyPe() == 0) {
  printf("Done allocating map for patches\n");
  printf("Patch level decomposition:\n");
  for (pid = 0;  pid < numpatches;  pid++) {
    msm::PatchDiagram& p = map.patchList[pid];
    int ia = p.nrange.ia();
    int ib = p.nrange.ib();
    int ja = p.nrange.ja();
    int jb = p.nrange.jb();
    int ka = p.nrange.ka();
    int kb = p.nrange.kb();
    printf("patch id=%d  [%d..%d] x [%d..%d] x [%d..%d]\n",
        pid, ia, ib, ja, jb, ka, kb);
  }
}
#endif
  if (CkMyPe() == 0) {
    iout << iINFO << "MSM finished creating map for patches\n" << endi;
  }

  // initialize grid of BlockDiagram for each level
  int polydeg = PolyDegree[approx];
  numGridCutoff = 0;
  for (level = 0;  level < nlevels;  level++) {
    msm::Grid<msm::BlockDiagram>& b = map.blockLevel[level];
    int bni = b.ni();
    int bnj = b.nj();
    int bnk = b.nk();
#ifdef MSM_SKIP_BEYOND_SPHERE
    int gia, gib, gja, gjb, gka, gkb;
    if (approx == C1HERMITE) {
      gia = map.gc_c1hermite[level].ia();
      gib = map.gc_c1hermite[level].ib();
      gja = map.gc_c1hermite[level].ja();
      gjb = map.gc_c1hermite[level].jb();
      gka = map.gc_c1hermite[level].ka();
      gkb = map.gc_c1hermite[level].kb();
    }
    else {
      gia = map.gc[level].ia();
      gib = map.gc[level].ib();
      gja = map.gc[level].ja();
      gjb = map.gc[level].jb();
      gka = map.gc[level].ka();
      gkb = map.gc[level].kb();
    }
#endif
#ifdef MSM_SKIP_TOO_DISTANT_BLOCKS
    int bsx = map.bsx[level];
    int bsy = map.bsy[level];
    int bsz = map.bsz[level];
#endif
#ifdef MSM_FOLD_FACTOR
    if (map.foldfactor[level].active) {
      bsx *= map.foldfactor[level].numrep.i;
      bsy *= map.foldfactor[level].numrep.j;
      bsz *= map.foldfactor[level].numrep.k;
    }
#endif
    for (k = 0;  k < bnk;  k++) {
      for (j = 0;  j < bnj;  j++) {
        for (i = 0;  i < bni;  i++) {

          // Grid cutoff calculation, sendAcross
          int ia = b(i,j,k).nrange.ia();
          int ib = b(i,j,k).nrange.ib();
          int ja = b(i,j,k).nrange.ja();
          int jb = b(i,j,k).nrange.jb();
          int ka = b(i,j,k).nrange.ka();
          int kb = b(i,j,k).nrange.kb();
          if (approx == C1HERMITE) {
            ia += map.gc_c1hermite[level].ia();
            ib += map.gc_c1hermite[level].ib();
            ja += map.gc_c1hermite[level].ja();
            jb += map.gc_c1hermite[level].jb();
            ka += map.gc_c1hermite[level].ka();
            kb += map.gc_c1hermite[level].kb();
          }
          else {
            ia += map.gc[level].ia();
            ib += map.gc[level].ib();
            ja += map.gc[level].ja();
            jb += map.gc[level].jb();
            ka += map.gc[level].ka();
            kb += map.gc[level].kb();
          }
          msm::Ivec na = map.clipIndexToLevel(msm::Ivec(ia,ja,ka), level);
          msm::Ivec nb = map.clipIndexToLevel(msm::Ivec(ib,jb,kb), level);
          b(i,j,k).nrangeCutoff.setbounds(na.i, nb.i, na.j, nb.j, na.k, nb.k);
          // determine sendAcross blocks
#ifdef MSM_FOLD_FACTOR
          msm::BlockIndex blower = map.blockOfGridIndexFold(na, level);
          msm::BlockIndex bupper = map.blockOfGridIndexFold(nb, level);
#else
          msm::BlockIndex blower = map.blockOfGridIndex(na, level);
          msm::BlockIndex bupper = map.blockOfGridIndex(nb, level);
#endif
          int maxarrlen = (bupper.n.i - blower.n.i + 1) *
            (bupper.n.j - blower.n.j + 1) * (bupper.n.k - blower.n.k + 1);
          b(i,j,k).sendAcross.setmax(maxarrlen);  // allocate send array
          b(i,j,k).indexGridCutoff.setmax(maxarrlen);  // alloc indexing
          b(i,j,k).recvGridCutoff.setmax(maxarrlen);  // alloc indexing
          // loop over sendAcross blocks
          int ii, jj, kk;
#if 0
          {
            msm::IndexRange& bn = b(i,j,k).nrange;
            printf("ME %4d   [%d..%d] x [%d..%d] x [%d..%d]\n",
                bn.nn(),
                bn.ia(), bn.ib(),
                bn.ja(), bn.jb(),
                bn.ka(), bn.kb());
          }
#endif
          for (kk = blower.n.k;  kk <= bupper.n.k;  kk++) {
            for (jj = blower.n.j;  jj <= bupper.n.j;  jj++) {
              for (ii = blower.n.i;  ii <= bupper.n.i;  ii++) {
#ifdef MSM_SKIP_TOO_DISTANT_BLOCKS
                // make sure that block (ii,jj,kk) interacts with (i,j,k)
                int si = sign(ii-i);
                int sj = sign(jj-j);
                int sk = sign(kk-k);
                int di = (ii-i)*bsx + si*(1-bsx);
                int dj = (jj-j)*bsy + sj*(1-bsy);
                int dk = (kk-k)*bsz + sk*(1-bsz);
                Vector d = di*hu + dj*hv + dk*hw;
                if (d.length2() >= 4*a*a) continue;
#endif
                // determine actual block and range to send to
                msm::BlockSend bs;
                bs.nblock.n = msm::Ivec(ii,jj,kk);
                bs.nblock.level = level;
#ifdef MSM_FOLD_FACTOR
                bs.nrange = map.clipBlockToIndexRangeFold(bs.nblock,
                    b(i,j,k).nrangeCutoff);
                map.wrapBlockSendFold(bs);  // wrap to true block index
#else
                bs.nrange = map.clipBlockToIndexRange(bs.nblock,
                    b(i,j,k).nrangeCutoff);
                map.wrapBlockSend(bs);  // wrap to true block index
#endif
#ifdef MSM_SKIP_BEYOND_SPHERE
#if 0
                printf("send to volume %4d   [%d..%d] x [%d..%d] x [%d..%d]\n",
                    bs.nrange.nn(),
                    bs.nrange.ia(), bs.nrange.ib(),
                    bs.nrange.ja(), bs.nrange.jb(),
                    bs.nrange.ka(), bs.nrange.kb());
#endif
                msm::IndexRange& bm = b(i,j,k).nrange;
                msm::IndexRange& bn = bs.nrange;
                int qia = bm.ia();
                int qib = bm.ib();
                int qja = bm.ja();
                int qjb = bm.jb();
                int qka = bm.ka();
                int qkb = bm.kb();
                int inc_in = (bn.ni() > 1 ? bn.ni()-1 : 1);
                int inc_jn = (bn.nj() > 1 ? bn.nj()-1 : 1);
                int inc_kn = (bn.nk() > 1 ? bn.nk()-1 : 1);
                // loop over corner points of potential grid
                int iscalc = 0;
                for (int kn = bn.ka();  kn <= bn.kb();  kn += inc_kn) {
                  for (int jn = bn.ja();  jn <= bn.jb();  jn += inc_jn) {
                    for (int in = bn.ia();  in <= bn.ib();  in += inc_in) {
                      // clip charges to weights
                      int mia = ( qia >= gia + in ? qia : gia + in );
                      int mib = ( qib <= gib + in ? qib : gib + in );
                      int mja = ( qja >= gja + jn ? qja : gja + jn );
                      int mjb = ( qjb <= gjb + jn ? qjb : gjb + jn );
                      int mka = ( qka >= gka + kn ? qka : gka + kn );
                      int mkb = ( qkb <= gkb + kn ? qkb : gkb + kn );
                      int inc_im = (mib-mia > 0 ? mib-mia : 1);
                      int inc_jm = (mjb-mja > 0 ? mjb-mja : 1);
                      int inc_km = (mkb-mka > 0 ? mkb-mka : 1);

                      // loop over corner points of charge grid
                      for (int km = mka;  km <= mkb;  km += inc_km) {
                        for (int jm = mja;  jm <= mjb;  jm += inc_jm) {
                          for (int im = mia;  im <= mib;  im += inc_im) {

                            Float g;
                            if (approx == C1HERMITE) {
                              g = map.gc_c1hermite[level](im-in,jm-jn,km-kn).melem[0];
                            }
                            else {
                              g = map.gc[level](im-in,jm-jn,km-kn);
                            }
                            iscalc |= (g != 0);
                          }
                        }
                      }

                    }
                  }
                }
                if ( ! iscalc) {
                  //printf("SKIPPING\n");  // XXX
                  continue;  // skip because overlap is beyond nonzero gc sphere
                }
#endif
                b(i,j,k).sendAcross.append(bs);
                b(i,j,k).indexGridCutoff.append(numGridCutoff);
                // receiving block records this grid cutoff ID
                b(bs.nblock_wrap.n).recvGridCutoff.append(numGridCutoff);
                // increment counter for receive block
                b(bs.nblock_wrap.n).numRecvsPotential++;

                numGridCutoff++;  // one MsmGridCutoff for each send across
              }
            }
          } // end loop over sendAcross blocks

          // Restriction, sendUp
          if (level < nlevels-1) {
            int ia2, ib2, ja2, jb2, ka2, kb2;
            ia = b(i,j,k).nrange.ia();
            ib = b(i,j,k).nrange.ib();
            ja = b(i,j,k).nrange.ja();
            jb = b(i,j,k).nrange.jb();
            ka = b(i,j,k).nrange.ka();
            kb = b(i,j,k).nrange.kb();
            // determine expansion of h-grid onto 2h-grid
            if ( ia==ib && ((ia & 1)==0) ) {
              ia2 = ib2 = ia / 2;
            }
            else {
              ia2 = (ia / 2) - ((polydeg+1) / 2) + 1;
              ib2 = ((ib+1) / 2) + ((polydeg+1) / 2) - 1;
            }
            if ( ja==jb && ((ja & 1)==0) ) {
              ja2 = jb2 = ja / 2;
            }
            else {
              ja2 = (ja / 2) - ((polydeg+1) / 2) + 1;
              jb2 = ((jb+1) / 2) + ((polydeg+1) / 2) - 1;
            }
            if ( ka==kb && ((ka & 1)==0) ) {
              ka2 = kb2 = ka / 2;
            }
            else {
              ka2 = (ka / 2) - ((polydeg+1) / 2) + 1;
              kb2 = ((kb+1) / 2) + ((polydeg+1) / 2) - 1;
            }
            // clip to boundaries of 2h-grid
            msm::Ivec na2, nb2;
            na2 = map.clipIndexToLevel(msm::Ivec(ia2,ja2,ka2), level+1);
            nb2 = map.clipIndexToLevel(msm::Ivec(ib2,jb2,kb2), level+1);
            b(i,j,k).nrangeRestricted.setbounds(na2.i, nb2.i, na2.j, nb2.j,
                na2.k, nb2.k);
            // determine sendUp blocks
            msm::BlockIndex blower = map.blockOfGridIndex(na2, level+1);
            msm::BlockIndex bupper = map.blockOfGridIndex(nb2, level+1);
            int maxarrlen = (bupper.n.i - blower.n.i + 1) *
              (bupper.n.j - blower.n.j + 1) * (bupper.n.k - blower.n.k + 1);
            b(i,j,k).sendUp.setmax(maxarrlen);  // allocate send array
            // loop over sendUp blocks
            int ii, jj, kk;
            for (kk = blower.n.k;  kk <= bupper.n.k;  kk++) {
              for (jj = blower.n.j;  jj <= bupper.n.j;  jj++) {
                for (ii = blower.n.i;  ii <= bupper.n.i;  ii++) {
                  // determine actual block and range to send to
                  msm::BlockSend bs;
                  bs.nblock.n = msm::Ivec(ii,jj,kk);
                  bs.nblock.level = level+1;
                  bs.nrange = map.clipBlockToIndexRange(bs.nblock,
                      b(i,j,k).nrangeRestricted);
                  map.wrapBlockSend(bs);  // wrap to true block index
                  b(i,j,k).sendUp.append(bs);
                  // increment counter for receive block
                  map.blockLevel[level+1](bs.nblock_wrap.n).numRecvsCharge++;
                }
              }
            } // end loop over sendUp blocks

          } // end if restriction

          // Prolongation, sendDown
          if (level > 0) {
            int ia2 = b(i,j,k).nrange.ia();
            int ib2 = b(i,j,k).nrange.ib();
            int ja2 = b(i,j,k).nrange.ja();
            int jb2 = b(i,j,k).nrange.jb();
            int ka2 = b(i,j,k).nrange.ka();
            int kb2 = b(i,j,k).nrange.kb();
            // determine expansion of 2h-grid onto h-grid
            ia = 2*ia2 - polydeg;
            ib = 2*ib2 + polydeg;
            ja = 2*ja2 - polydeg;
            jb = 2*jb2 + polydeg;
            ka = 2*ka2 - polydeg;
            kb = 2*kb2 + polydeg;
            // clip to boundaries of h-grid
            msm::Ivec na, nb;
            na = map.clipIndexToLevel(msm::Ivec(ia,ja,ka), level-1);
            nb = map.clipIndexToLevel(msm::Ivec(ib,jb,kb), level-1);
            b(i,j,k).nrangeProlongated.setbounds(na.i, nb.i, na.j, nb.j,
                na.k, nb.k);
            // determine sendDown blocks
            msm::BlockIndex blower = map.blockOfGridIndex(na, level-1);
            msm::BlockIndex bupper = map.blockOfGridIndex(nb, level-1);
            int maxarrlen = (bupper.n.i - blower.n.i + 1) *
              (bupper.n.j - blower.n.j + 1) * (bupper.n.k - blower.n.k + 1);
            b(i,j,k).sendDown.setmax(maxarrlen);  // allocate send array
            // loop over sendUp blocks
            int ii, jj, kk;
            for (kk = blower.n.k;  kk <= bupper.n.k;  kk++) {
              for (jj = blower.n.j;  jj <= bupper.n.j;  jj++) {
                for (ii = blower.n.i;  ii <= bupper.n.i;  ii++) {
                  // determine actual block and range to send to
                  msm::BlockSend bs;
                  bs.nblock.n = msm::Ivec(ii,jj,kk);
                  bs.nblock.level = level-1;
                  bs.nrange = map.clipBlockToIndexRange(bs.nblock,
                      b(i,j,k).nrangeProlongated);
                  map.wrapBlockSend(bs);  // wrap to true block index
                  b(i,j,k).sendDown.append(bs);
                  // increment counter for receive block
                  map.blockLevel[level-1](bs.nblock_wrap.n).numRecvsPotential++;
                }
              }
            } // end loop over sendDown blocks

          } // end if prolongation

#ifdef MSM_REDUCE_GRID
          // using a reduction decreases the number of messages
          // from MsmGridCutoff elements to just 1
          b(i,j,k).numRecvsPotential -= ( b(i,j,k).indexGridCutoff.len() - 1 );
#endif

        }
      }
    } // end loop over block diagram

  } // end loop over levels
  // end of Map setup

  // allocate chare arrays

  if (1) {
    PatchMap *pm = PatchMap::Object();
    patchPtr.resize( pm->numPatches() );
    for (int i = 0;  i < pm->numPatches();  i++) {
      patchPtr[i] = NULL;
    }
#ifdef DEBUG_MSM_VERBOSE
    printf("Allocating patchPtr array length %d\n", pm->numPatches());
#endif
    if (CkMyPe() == 0) {
      iout << iINFO << "MSM has " << pm->numPatches()
                    << " interpolation / anterpolation objects"
                    << " (one per patch)\n" << endi;
    }
  }

#ifdef MSM_NODE_MAPPING
  if (1) {
    // Node aware initial assignment of chares
    //
    // Create map object for each 3D chare array of MsmBlock and the
    // 1D chare array of MsmGridCutoff.  Design map to equally distribute
    // blocks across nodes, assigned to node PEs in round robin manner.  
    // Attempt to reduce internode communication bandwidth by assigning 
    // each MsmGridCutoff element to either its source node or its 
    // destination node, again assigned to node PEs in round robin manner.
#if 0
    // for testing
#if 0
    int numNodes = 16;
    int numPes = 512;
#else
    int numNodes = 32;
    int numPes = 1024;
#endif
#else
    int numNodes = CkNumNodes();
    int numPes = CkNumPes();
#endif
    int numPesPerNode = numPes / numNodes;
    int numBlocks = 0;  // find total number of blocks
    for (level = 0;  level < nlevels;  level++) {
      numBlocks += map.blockLevel[level].nn();
    }

    // final result is arrays for blocks and gcuts, each with pe number
    blockAssign.resize(numBlocks);
    gcutAssign.resize(numGridCutoff);
    //printf("XXX numBlocks = %d\n", numBlocks);
    //printf("XXX numGridCutoff = %d\n", numGridCutoff);

    msm::Array<float> blockWork(numBlocks);
    msm::Array<float> gcutWork(numGridCutoff);

    msm::Array<float> nodeWork(numNodes);
    nodeWork.reset(0);
#ifdef MSM_NODE_MAPPING_STATS
    msm::Array<float> peWork(numPes);
    peWork.reset(0);
#endif

    msm::PriorityQueue<WorkIndex> nodeQueue(numNodes);
    for (n = 0;  n < numNodes;  n++) {
      nodeQueue.insert(WorkIndex(0, n));
    }

    int bindex = 0;  // index for block array
    for (level = 0;  level < nlevels;  level++) {
      msm::Grid<msm::BlockDiagram>& b = map.blockLevel[level];
      int bni = b.ni();
      int bnj = b.nj();
      int bnk = b.nk();
      for (k = 0;  k < bnk;  k++) { // for all blocks
        for (j = 0;  j < bnj;  j++) {
          for (i = 0;  i < bni;  i++) {
            WorkIndex wn;
            nodeQueue.remove(wn);
            float bw = calcBlockWork(b(i,j,k));
            blockAssign[bindex] = wn.index;
            nodeWork[wn.index] += bw;
            wn.work += bw;
            blockWork[bindex] = bw;
            nodeQueue.insert(wn);
            bindex++;
          }
        }
      } // end for all blocks
    } // end for all levels

#if 0
    for (n = 0;  n < numBlocks;  n++) {
      WorkIndex wn;
      nodeQueue.remove(wn);
      float bw = calcBlockWork(n);
      blockAssign[n] = wn.index;
      nodeWork[wn.index] += bw;
      wn.work += bw;
      blockWork[n] = bw;
      nodeQueue.insert(wn);
    }
#endif

    // assign grid cutoff objects to nodes (gcutAssign)
    // choose whichever of source or destination node has less work
    int gindex = 0;  // index for grid cutoff array
    for (level = 0;  level < nlevels;  level++) { // for all levels
      msm::Grid<msm::BlockDiagram>& b = map.blockLevel[level];
      int bni = b.ni();
      int bnj = b.nj();
      int bnk = b.nk();
      for (k = 0;  k < bnk;  k++) { // for all blocks
        for (j = 0;  j < bnj;  j++) {
          for (i = 0;  i < bni;  i++) {
            int isrc = blockFlatIndex(level, i, j, k);
            int nsrc = blockAssign[isrc];  // node block isrc is assigned
            int numSendAcross = b(i,j,k).sendAcross.len();
            ASSERT( numSendAcross == b(i,j,k).indexGridCutoff.len() );
            for (n = 0;  n < numSendAcross;  n++) {
              msm::BlockSend& bs = b(i,j,k).sendAcross[n];
              msm::BlockIndex& bn = bs.nblock_wrap;
              int idest = blockFlatIndex(level, bn.n.i, bn.n.j, bn.n.k);
              int ndest = blockAssign[idest];  // node block idest is assigned
              gcutWork[gindex] = calcGcutWork(bs);
              if (nodeWork[nsrc] <= nodeWork[ndest]) {
                gcutAssign[gindex] = nsrc;
                nodeWork[nsrc] += gcutWork[gindex];
              }
              else {
                gcutAssign[gindex] = ndest;
                nodeWork[ndest] += gcutWork[gindex];
              }
              gindex++;
            } // end for numSendAcross
          }
        }
      } // end for all blocks
    } // end for all levels

    msm::Array< msm::PriorityQueue<WorkIndex> > peQueue(numNodes);
    for (n = 0;  n < numNodes;  n++) {
      peQueue[n].init(numPesPerNode);
      for (int poff = 0;  poff < numPesPerNode;  poff++) {
        peQueue[n].insert(WorkIndex(0, n*numPesPerNode + poff));
      }
    }

    for (n = 0;  n < numBlocks;  n++) {
      WorkIndex wn;
      int node = blockAssign[n];
      peQueue[node].remove(wn);
      blockAssign[n] = wn.index;
      wn.work += blockWork[n];
      peQueue[node].insert(wn);
#ifdef MSM_NODE_MAPPING_STATS
      peWork[wn.index] += blockWork[n];
#endif
    }

    for (n = 0;  n < numGridCutoff;  n++) {
      WorkIndex wn;
      int node = gcutAssign[n];
      peQueue[node].remove(wn);
      gcutAssign[n] = wn.index;
      wn.work += gcutWork[n];
      peQueue[node].insert(wn);
#ifdef MSM_NODE_MAPPING_STATS
      peWork[wn.index] += gcutWork[n];
#endif
    }

#ifdef MSM_NODE_MAPPING_STATS
    if (CkMyPe() == 0) {
      printf("Mapping of MSM work (showing scaled estimated work units):\n");
      for (n = 0;  n < numNodes;  n++) {
        printf("    node %d   work %8.3f\n", n, nodeWork[n]);
        for (int poff = 0;  poff < numPesPerNode;  poff++) {
          int p = n*numPesPerNode + poff;
          printf("        pe %d     work %8.3f\n", p, peWork[p]);
        }
      }
      //CkExit();
    }
#endif

#if 0
    int numBlocks = 0;  // find total number of blocks
    for (level = 0;  level < nlevels;  level++) {
      numBlocks += map.blockLevel[level].nn();
    }

    // final result is arrays for blocks and gcuts, each with pe number
    blockAssign.resize(numBlocks);
    gcutAssign.resize(numGridCutoff);

    nodecnt.resize(numNodes);

    // assign blocks to nodes
    // the following algorithm divides as evenly as possible the
    // blocks across the nodes
    int r = numBlocks % numNodes;
    int q = numBlocks / numNodes;
    int qp = q + 1;
    for (n = 0;  n < numNodes - r;  n++) {
      int moffset = n * q;
      for (int m = 0;  m < q;  m++) {
        blockAssign[moffset + m] = n;
      }
      nodecnt[n] = q;
    }
    for ( ;  n < numNodes;  n++) {
      int moffset = (numNodes - r)*q + (n - (numNodes - r))*qp;
      for (int m = 0;  m < qp;  m++) {
        blockAssign[moffset + m] = n;
      }
      nodecnt[n] = qp;
    }
#if 0
    if (CkMyPe() == 0) {
      CkPrintf("%d objects to %d nodes\n", q, numNodes-r);
      if (r != 0) {
        CkPrintf("%d objects to %d nodes\n", qp, r);
      }
      CkPrintf("%d  =?  %d\n", (numNodes-r)*q + r*qp, numBlocks);
    }
#endif

    // assign grid cutoff objects to nodes (gcutAssign)
    // choose whichever of source or destination node has less work
    int gindex = 0;  // index for grid cutoff array
    for (level = 0;  level < nlevels;  level++) { // for all levels
      msm::Grid<msm::BlockDiagram>& b = map.blockLevel[level];
      int bni = b.ni();
      int bnj = b.nj();
      int bnk = b.nk();
      for (k = 0;  k < bnk;  k++) { // for all blocks
        for (j = 0;  j < bnj;  j++) {
          for (i = 0;  i < bni;  i++) {
            int isrc = blockFlatIndex(level, i, j, k);
            int nsrc = blockAssign[isrc];  // node block isrc is assigned
            int numSendAcross = b(i,j,k).sendAcross.len();
            ASSERT( numSendAcross == b(i,j,k).indexGridCutoff.len() );
            for (n = 0;  n < numSendAcross;  n++) {
              msm::BlockIndex &bn = b(i,j,k).sendAcross[n].nblock_wrap;
              int idest = blockFlatIndex(level, bn.n.i, bn.n.j, bn.n.k);
              int ndest = blockAssign[idest];  // node block idest is assigned
              // assign this grid cutoff work to least subscribed node
              if (nodecnt[nsrc] <= nodecnt[ndest]) {
                gcutAssign[gindex] = nsrc;
                nodecnt[nsrc]++;
              }
              else {
                gcutAssign[gindex] = ndest;
                nodecnt[ndest]++;
              }
              gindex++;
            } // end for numSendAcross
          }
        }
      } // end for all blocks
    } // end for all levels

    // now change the node assignments into PE assignments
    // use round robin assignment to PEs within each node
    int ppn = numPes / numNodes;  // num PEs per node
    // reset nodecnt - this array will now store PE offset for that node
    for (n = 0;  n < numNodes;  n++)  nodecnt[n] = 0;
    for (n = 0;  n < numBlocks;  n++) {
      int node = blockAssign[n];
      blockAssign[n] = node * ppn + nodecnt[node];  // PE number within node
      nodecnt[node]++;  // increment to next PE
      if (nodecnt[node] >= ppn)  nodecnt[node] = 0;  // with wrap around
    }
    for (n = 0;  n < numGridCutoff;  n++) {
      int node = gcutAssign[n];
      gcutAssign[n] = node * ppn + nodecnt[node];  // PE number within node
      nodecnt[node]++;  // increment to next PE
      if (nodecnt[node] >= ppn)  nodecnt[node] = 0;  // with wrap around
    }

    // print mapping
#if 0
    if (CkMyPe() == 0) {
      for (n = 0;  n < numBlocks;  n++) {
        CkPrintf("block %d:   node=%d  pe=%d\n",
            n, blockAssign[n]/ppn, blockAssign[n]);
      }
#if 0
      for (n = 0;  n < numGridCutoff;  n++) {
        CkPrintf("grid cutoff %d:   node=%d  pe=%d\n",
            n, gcutAssign[n]/ppn, gcutAssign[n]);
      }
#endif
    }
#endif

#endif // 0

  } // end node aware initial assignment of chares
#endif // MSM_NODE_MAPPING

  if (CkMyPe() == 0) {

    // on PE 0, create 3D chare array of MsmBlock for each level;
    // broadcast this array of proxies to the rest of the group
    if (approx == C1HERMITE) {
      msmC1HermiteBlock.resize(nlevels);
    }
    else {
      msmBlock.resize(nlevels);
    }
    for (level = 0;  level < nlevels;  level++) {
      int ni = map.blockLevel[level].ni();
      int nj = map.blockLevel[level].nj();
      int nk = map.blockLevel[level].nk();
#ifdef MSM_NODE_MAPPING
      CkPrintf("Using MsmBlockMap for level %d\n", level);
      CProxy_MsmBlockMap blockMap = CProxy_MsmBlockMap::ckNew(level);
      CkArrayOptions opts(ni, nj, nk);
      opts.setMap(blockMap);
      if (approx == C1HERMITE) {
        msmC1HermiteBlock[level] =
          CProxy_MsmC1HermiteBlock::ckNew(level, opts);
      }
      else {
        msmBlock[level] = CProxy_MsmBlock::ckNew(level, opts);
      }
#else
      if (approx == C1HERMITE) {
        msmC1HermiteBlock[level] =
          CProxy_MsmC1HermiteBlock::ckNew(level, ni, nj, nk);
      }
      else {
        msmBlock[level] = CProxy_MsmBlock::ckNew(level, ni, nj, nk);
      }
#endif
#ifdef DEBUG_MSM_VERBOSE
      printf("Create MsmBlock[%d] 3D chare array ( %d x %d x %d )\n",
          level, ni, nj, nk);
#endif
      char msg[128];
      int nijk = ni * nj * nk;
      sprintf(msg, "MSM grid level %d decomposed into %d block%s"
          " ( %d x %d x %d )\n",
          level, nijk, (nijk==1 ? "" : "s"), ni, nj, nk);
      iout << iINFO << msg;
    }
    if (approx == C1HERMITE) {
      MsmC1HermiteBlockProxyMsg *msg = new MsmC1HermiteBlockProxyMsg;
      msg->put(msmC1HermiteBlock);
      msmProxy.recvMsmC1HermiteBlockProxy(msg);  // broadcast
    }
    else {
      MsmBlockProxyMsg *msg = new MsmBlockProxyMsg;
      msg->put(msmBlock);
      msmProxy.recvMsmBlockProxy(msg);  // broadcast
    }

#ifdef MSM_GRID_CUTOFF_DECOMP
    // on PE 0, create 1D chare array of MsmGridCutoff
    // broadcast this array proxy to the rest of the group
#ifdef MSM_NODE_MAPPING
    CkPrintf("Using MsmGridCutoffMap\n");
    CProxy_MsmGridCutoffMap gcutMap = CProxy_MsmGridCutoffMap::ckNew();
    CkArrayOptions optsgcut(numGridCutoff);
    optsgcut.setMap(gcutMap);
    if (approx == C1HERMITE) {
      msmC1HermiteGridCutoff = CProxy_MsmC1HermiteGridCutoff::ckNew(optsgcut);
    }
    else {
      msmGridCutoff = CProxy_MsmGridCutoff::ckNew(optsgcut);
    }
#else
    if (approx == C1HERMITE) {
      msmC1HermiteGridCutoff =
        CProxy_MsmC1HermiteGridCutoff::ckNew(numGridCutoff);
    }
    else {
      msmGridCutoff = CProxy_MsmGridCutoff::ckNew(numGridCutoff);
    }
#endif
    if (approx == C1HERMITE) {
      MsmC1HermiteGridCutoffProxyMsg *gcmsg =
        new MsmC1HermiteGridCutoffProxyMsg;
      gcmsg->put(&msmC1HermiteGridCutoff);
      msmProxy.recvMsmC1HermiteGridCutoffProxy(gcmsg);
    }
    else {
      MsmGridCutoffProxyMsg *gcmsg = new MsmGridCutoffProxyMsg;
      gcmsg->put(&msmGridCutoff);
      msmProxy.recvMsmGridCutoffProxy(gcmsg);
    }

    // XXX PE 0 initializes each MsmGridCutoff
    // one-to-many
    // for M length chare array, better for each PE to initialize M/P?
    for (level = 0;  level < nlevels;  level++) { // for all levels
      msm::Grid<msm::BlockDiagram>& b = map.blockLevel[level];
      int bni = b.ni();
      int bnj = b.nj();
      int bnk = b.nk();
      for (k = 0;  k < bnk;  k++) { // for all blocks
        for (j = 0;  j < bnj;  j++) {
          for (i = 0;  i < bni;  i++) {
            // source for charges
            msm::BlockIndex bi = msm::BlockIndex(level, msm::Ivec(i,j,k));
            int numSendAcross = b(i,j,k).sendAcross.len();
            ASSERT( numSendAcross == b(i,j,k).indexGridCutoff.len() );
            // for this source, loop over destinations for potentials
            for (n = 0;  n < numSendAcross;  n++) {
              msm::BlockSend &bs = b(i,j,k).sendAcross[n];
              int index = b(i,j,k).indexGridCutoff[n];
              MsmGridCutoffInitMsg *bsmsg = new MsmGridCutoffInitMsg(bi, bs);
              if (approx == C1HERMITE) {
                msmC1HermiteGridCutoff[index].setup(bsmsg);
              }
              else {
                msmGridCutoff[index].setup(bsmsg);
              }
            } // traverse sendAcross, indexGridCutoff arrays

          }
        }
      } // end for all blocks

    } // end for all levels

    iout << iINFO << "MSM grid cutoff calculation decomposed into "
      << numGridCutoff << " work objects\n";
#endif
    iout << endi;
  }

#ifdef DEBUG_MSM_VERBOSE
  printf("end of initialization\n");
#endif
} // ComputeMsmMgr::initialize()

void ComputeMsmMgr::recvMsmBlockProxy(MsmBlockProxyMsg *msg)
{
  msg->get(msmBlock);
  delete(msg);
}

void ComputeMsmMgr::recvMsmGridCutoffProxy(MsmGridCutoffProxyMsg *msg)
{
  msg->get(&msmGridCutoff);
  delete(msg);
}

void ComputeMsmMgr::recvMsmC1HermiteBlockProxy(
    MsmC1HermiteBlockProxyMsg *msg
    )
{
  msg->get(msmC1HermiteBlock);
  delete(msg);
}

void ComputeMsmMgr::recvMsmC1HermiteGridCutoffProxy(
    MsmC1HermiteGridCutoffProxyMsg *msg
    )
{
  msg->get(&msmC1HermiteGridCutoff);
  delete(msg);
}

void ComputeMsmMgr::update(CkQdMsg *msg)
{
#ifdef DEBUG_MSM_VERBOSE
  printf("ComputeMsmMgr:  update() PE %d\n", CkMyPe());
#endif
  delete msg;

  // have to setup sections AFTER initialization is finished
  if (CkMyPe() == 0) {
    for (int level = 0;  level < nlevels;  level++) {
      if (approx == C1HERMITE) {
        msmC1HermiteBlock[level].setupSections();
      }
      else {
        msmBlock[level].setupSections();
      }
    }
  }

  // XXX how do update for constant pressure simulation?
}


void ComputeMsmMgr::compute(msm::Array<int>& patchIDList)
{
#ifdef DEBUG_MSM_VERBOSE
  printf("ComputeMsmMgr:  compute() PE=%d\n", CkMyPe());
#endif

  int n; 
  for (n = 0;  n < patchIDList.len();  n++) {
    int patchID = patchIDList[n];
    if (patchPtr[patchID] == NULL) {
      char msg[100];
      snprintf(msg, sizeof(msg),
          "Expected MSM data for patch %d does not exist on PE %d",
          patchID, CkMyPe());
      NAMD_die(msg);
    }
    if (approx == C1HERMITE) {
      patchPtr[patchID]->anterpolationC1Hermite();
    }
    else {
      patchPtr[patchID]->anterpolation();
    }
    // all else should follow from here
  }
  return;
}


void ComputeMsmMgr::addPotential(GridMsg *gm)
{
  int pid;  // receive patch ID
  int pseq;
  if (approx == C1HERMITE) {
    gm->get(subgrid_c1hermite, pid, pseq);
  }
  else {
    gm->get(subgrid, pid, pseq);
  }
  delete gm;
  if (patchPtr[pid] == NULL) {
    char msg[100];
    snprintf(msg, sizeof(msg), "Expecting patch %d to exist on PE %d",
        pid, CkMyPe());
    NAMD_die(msg);
  }
  if (approx == C1HERMITE) {
    patchPtr[pid]->addPotentialC1Hermite(subgrid_c1hermite);
  }
  else {
    patchPtr[pid]->addPotential(subgrid);
  }
}


void ComputeMsmMgr::doneCompute()
{
  msmCompute->saveResults();
}


//////////////////////////////////////////////////////////////////////////////
//
//  ComputeMsm
//  MSM compute objects, starts and finishes calculation;
//  there is up to one compute object per PE
//

ComputeMsm::ComputeMsm(ComputeID c) : ComputeHomePatches(c)
{
  CProxy_ComputeMsmMgr::ckLocalBranch(
      CkpvAccess(BOCclass_group).computeMsmMgr)->setCompute(this);
  SimParameters *simParams = Node::Object()->simParameters;
  qscaling = sqrtf(COULOMB / simParams->dielectric);
  reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);
#ifdef DEBUG_MSM_VERBOSE
  printf("ComputeMsm:  (constructor) PE=%d\n", CkMyPe());
#endif
}

ComputeMsm::~ComputeMsm()
{
  // free memory
#ifdef DEBUG_MSM_VERBOSE
  printf("ComputeMsm:  (destructor) PE=%d\n", CkMyPe());
#endif
}

void ComputeMsm::doWork()
{
  // for each patch do stuff
#ifdef DEBUG_MSM_VERBOSE
  printf("ComputeMsm:  doWork() PE=%d\n", CkMyPe());
#endif

#if 0
#ifdef MSM_TIMING
  myMgr->initTiming();
#endif
#ifdef MSM_PROFILING
  myMgr->initProfiling();
#endif
#endif

  // patchList is inherited from ComputeHomePatches
  ResizeArrayIter<PatchElem> ap(patchList);
  numLocalPatches = patchList.size();
  cntLocalPatches = 0;
  ASSERT(cntLocalPatches < numLocalPatches);

#ifdef DEBUG_MSM_VERBOSE
  printf("patchList size = %d\n", patchList.size() );
#endif

  // Skip computations if nothing to do.
  if ( ! patchList[0].p->flags.doFullElectrostatics ) {
    for (ap = ap.begin();  ap != ap.end();  ap++) {
      CompAtom *x = (*ap).positionBox->open();
      Results *r = (*ap).forceBox->open();
      (*ap).positionBox->close(&x);
      (*ap).forceBox->close(&r);
    }
    reduction->submit();
    return;
  }
  msm::Map& map = myMgr->mapData();
  // This is the patchPtr array for MSM; any local patch will be set up
  // with a non-NULL pointer to its supporting data structure.
  msm::PatchPtrArray& patchPtr = myMgr->patchPtrArray();
  // also store just a list of IDs for the local patches
  msm::Array<int> patchIDList(numLocalPatches);
  patchIDList.resize(0);  // to use append on pre-allocated array buffer
  int cnt=0, n;
  for (ap = ap.begin();  ap != ap.end();  ap++) {
    CompAtom *x = (*ap).positionBox->open();
    CompAtomExt *xExt = (*ap).p->getCompAtomExtInfo();
    if ( patchList[0].p->flags.doMolly ) {
      (*ap).positionBox->close(&x);
      x = (*ap).avgPositionBox->open();
    }
    int numAtoms = (*ap).p->getNumAtoms();
    int patchID = (*ap).patchID;
    patchIDList.append(patchID);
    if (patchPtr[patchID] == NULL) {
      // create PatchData if it doesn't exist for this patchID
      patchPtr[patchID] = new msm::PatchData(myMgr, patchID);
#ifdef DEBUG_MSM_VERBOSE
      printf("Creating new PatchData:  patchID=%d  PE=%d\n",
          patchID, CkMyPe());
#endif
    }
    msm::PatchData& patch = *(patchPtr[patchID]);
    patch.init(numAtoms);
    msm::AtomCoordArray& coord = patch.coordArray();
    ASSERT(coord.len() == numAtoms);
    for (n = 0;  n < numAtoms;  n++) {
      coord[n].position = x[n].position;
      coord[n].charge = qscaling * x[n].charge;
      coord[n].id = xExt[n].id;
    }
    if ( patchList[0].p->flags.doMolly ) {
      (*ap).avgPositionBox->close(&x);
    }
    else {
      (*ap).positionBox->close(&x);
    }
    patch.sequence = sequence();
  }

  myMgr->compute(patchIDList);
}

void ComputeMsm::saveResults()
{
  if (++cntLocalPatches != numLocalPatches) return;

  // NAMD patches
  ResizeArrayIter<PatchElem> ap(patchList);
#ifdef DEBUG_MSM
  for (ap = ap.begin();  ap != ap.end();  ap++) {
    int patchID = (*ap).patchID;
    ASSERT(myMgr->patchPtrArray()[patchID]->cntRecvs ==
        myMgr->mapData().patchList[patchID].numRecvs);
  }
#endif

  // get results from ComputeMsmMgr
  msm::PatchPtrArray& patchPtr = myMgr->patchPtrArray();

#ifdef DEBUG_MSM_VERBOSE
  printf("ComputeMsm:  saveResults() PE=%d\n", CkMyPe());
#endif
  // store force updates
  // submit reductions

  // add in forces
  int cnt=0, n;
  for (ap = ap.begin(); ap != ap.end(); ap++) {
    Results *r = (*ap).forceBox->open();
    Force *f = r->f[Results::slow];
    int numAtoms = (*ap).p->getNumAtoms();
    int patchID = (*ap).patchID;
    if (patchPtr[patchID] == NULL) {
      char msg[100];
      snprintf(msg, sizeof(msg), "Expecting patch %d to exist on PE %d",
          patchID, CkMyPe());
      NAMD_die(msg);
    }
    msm::PatchData& patch = *(patchPtr[patchID]);
    ASSERT(numAtoms == patch.force.len() );
    for (n = 0;  n < numAtoms;  n++) {
      f[n] += patch.force[n];
    }
    (*ap).forceBox->close(&r);

    reduction->item(REDUCTION_ELECT_ENERGY_SLOW) += patch.energy;
//    reduction->item(REDUCTION_VIRIAL_SLOW_XX) += patch.virial[0][0];
//    reduction->item(REDUCTION_VIRIAL_SLOW_XY) += patch.virial[0][1];
//    reduction->item(REDUCTION_VIRIAL_SLOW_XZ) += patch.virial[0][2];
//    reduction->item(REDUCTION_VIRIAL_SLOW_YX) += patch.virial[1][0];
//    reduction->item(REDUCTION_VIRIAL_SLOW_YY) += patch.virial[1][1];
//    reduction->item(REDUCTION_VIRIAL_SLOW_YZ) += patch.virial[1][2];
//    reduction->item(REDUCTION_VIRIAL_SLOW_ZX) += patch.virial[2][0];
//    reduction->item(REDUCTION_VIRIAL_SLOW_ZY) += patch.virial[2][1];
//    reduction->item(REDUCTION_VIRIAL_SLOW_ZZ) += patch.virial[2][2];
    Float *virial = myMgr->virial;
    reduction->item(REDUCTION_VIRIAL_SLOW_XX) += virial[ComputeMsmMgr::VXX];
    reduction->item(REDUCTION_VIRIAL_SLOW_XY) += virial[ComputeMsmMgr::VXY];
    reduction->item(REDUCTION_VIRIAL_SLOW_XZ) += virial[ComputeMsmMgr::VXZ];
    reduction->item(REDUCTION_VIRIAL_SLOW_YX) += virial[ComputeMsmMgr::VXY];
    reduction->item(REDUCTION_VIRIAL_SLOW_YY) += virial[ComputeMsmMgr::VYY];
    reduction->item(REDUCTION_VIRIAL_SLOW_YZ) += virial[ComputeMsmMgr::VYZ];
    reduction->item(REDUCTION_VIRIAL_SLOW_ZX) += virial[ComputeMsmMgr::VXZ];
    reduction->item(REDUCTION_VIRIAL_SLOW_ZY) += virial[ComputeMsmMgr::VYZ];
    reduction->item(REDUCTION_VIRIAL_SLOW_ZZ) += virial[ComputeMsmMgr::VZZ];
  }
  reduction->submit();
}

// method definitions for PatchData
namespace msm {

  PatchData::PatchData(ComputeMsmMgr *pmgr, int pid) {
    mgr = pmgr;
    map = &(mgr->mapData());
    patchID = pid;
    //PatchMap *pm = PatchMap::Object();
    pd = &(map->patchList[pid]);
    if (mgr->approx == ComputeMsmMgr::C1HERMITE) {
      qh_c1hermite.init(pd->nrange);
      eh_c1hermite.init(pd->nrange);
      subgrid_c1hermite.resize(map->bsx[0] * map->bsy[0] * map->bsz[0]);
    }
    else {
      qh.init(pd->nrange);
      eh.init(pd->nrange);
      subgrid.resize(map->bsx[0] * map->bsy[0] * map->bsz[0]);
    }
#ifdef MSM_TIMING
    mgr->addTiming();
#endif
  }

  void PatchData::init(int natoms) {
    coord.resize(natoms);
    force.resize(natoms);
    cntRecvs = 0;
    energy = 0;
    //memset(virial, 0, 3*3*sizeof(BigReal));
    for (int i = 0;  i < natoms;  i++)  force[i] = 0;
    if (mgr->approx == ComputeMsmMgr::C1HERMITE) {
      qh_c1hermite.reset(0);
      eh_c1hermite.reset(0);
    }
    else {
      qh.reset(0);
      eh.reset(0);
    }
  }

  void PatchData::anterpolation() {
#ifdef DEBUG_MSM_GRID
    printf("patchID %d:  anterpolation\n", patchID);
#endif

#ifdef MSM_TIMING
    double startTime, stopTime;
    startTime = CkWallTimer();
#endif
#ifndef MSM_COMM_ONLY
    Float xphi[ComputeMsmMgr::MAX_POLY_DEGREE+1];
    Float yphi[ComputeMsmMgr::MAX_POLY_DEGREE+1];
    Float zphi[ComputeMsmMgr::MAX_POLY_DEGREE+1];

    const Double rs_edge = Double( mgr->s_edge );
    const int s_size = ComputeMsmMgr::PolyDegree[mgr->approx] + 1;

    const int ia = qh.ia();
    const int ib = qh.ib();
    const int ja = qh.ja();
    const int jb = qh.jb();
    const int ka = qh.ka();
    const int kb = qh.kb();
    const int ni = qh.ni();
    const int nj = qh.nj();
    Float *qhbuffer = qh.data().buffer();

    // loop over atoms
    for (int n = 0;  n < coord.len();  n++) {
      Float q = coord[n].charge;
      if (0==q) continue;

      ScaledPosition s = mgr->lattice.scale(coord[n].position);

      BigReal sx_hx = (s.x - mgr->sglower.x) * mgr->shx_1;
      BigReal sy_hy = (s.y - mgr->sglower.y) * mgr->shy_1;
      BigReal sz_hz = (s.z - mgr->sglower.z) * mgr->shz_1;

      BigReal xlo = floor(sx_hx) - rs_edge;
      BigReal ylo = floor(sy_hy) - rs_edge;
      BigReal zlo = floor(sz_hz) - rs_edge;

      // calculate Phi stencils along each dimension
      Float xdelta = Float(sx_hx - xlo);
      mgr->stencil_1d(xphi, xdelta);
      Float ydelta = Float(sy_hy - ylo);
      mgr->stencil_1d(yphi, ydelta);
      Float zdelta = Float(sz_hz - zlo);
      mgr->stencil_1d(zphi, zdelta);

      int ilo = int(xlo);
      int jlo = int(ylo);
      int klo = int(zlo);

      // test to see if stencil is within edges of grid
      int iswithin = ( ia <= ilo && (ilo+(s_size-1)) <= ib &&
                       ja <= jlo && (jlo+(s_size-1)) <= jb &&
                       ka <= klo && (klo+(s_size-1)) <= kb );

      if ( ! iswithin ) {
#if 0
        printf("PE %d:  atom %d:  pos= %g %g %g  patchID=%d\n",
            CkMyPe(), coord[n].id,
            coord[n].position.x, coord[n].position.y, coord[n].position.z,
            patchID);
        printf("PE %d:  atom subgrid [%d..%d] x [%d..%d] x [%d..%d]\n",
            CkMyPe(),
            ilo, ilo+s_size-1, jlo, jlo+s_size-1, klo, klo+s_size-1);
        printf("PE %d:  patch grid [%d..%d] x [%d..%d] x [%d..%d]\n",
            CkMyPe(),
            ia, ib, ja, jb, ka, kb);
#endif
        char msg[100];
        snprintf(msg, sizeof(msg), "Atom %d is outside of the MSM grid.",
            coord[n].id);
        NAMD_die(msg);
      }

      // determine charge on cube of grid points around atom
      for (int k = 0;  k < s_size;  k++) {
        int koff = ((k+klo) - ka) * nj;
        Float ck = zphi[k] * q;
        for (int j = 0;  j < s_size;  j++) {
          int jkoff = (koff + (j+jlo) - ja) * ni;
          Float cjk = yphi[j] * ck;
          for (int i = 0;  i < s_size;  i++) {
            int ijkoff = jkoff + (i+ilo) - ia;
            qhbuffer[ijkoff] += xphi[i] * cjk;
          }
        }
      }

    } // end loop over atoms
#endif // !MSM_COMM_ONLY
#ifdef MSM_TIMING
    stopTime = CkWallTimer();
    mgr->msmTiming[MsmTimer::ANTERP] += stopTime - startTime;
#endif

    sendCharge();
  }

  void PatchData::sendCharge() {
#ifdef MSM_TIMING
    double startTime, stopTime;
#endif
    int priority = 1;
    // buffer portions of grid to send to Blocks on level 0
    // allocate the largest buffer space we'll need
    //Grid<BigReal> subgrid;
    //subgrid.resize(map->bsx[0] * map->bsy[0] * map->bsz[0]);
    for (int n = 0;  n < pd->send.len();  n++) {
#ifdef MSM_TIMING
      startTime = CkWallTimer();
#endif
      // initialize the proper subgrid indexing range
      subgrid.init( pd->send[n].nrange );
      // extract the values from the larger grid into the subgrid
      qh.extract(subgrid);
      // translate the subgrid indexing range to match the MSM block
      subgrid.updateLower( pd->send[n].nrange_wrap.lower() );
      // add the subgrid charges into the block
      BlockIndex& bindex = pd->send[n].nblock_wrap;
      // place subgrid into message
      int msgsz = subgrid.data().len() * sizeof(Float);
      GridMsg *gm = new(msgsz, sizeof(int)) GridMsg;
      SET_PRIORITY(gm, sequence, MSM_PRIORITY + priority);
      gm->put(subgrid, bindex.level, sequence);
#ifdef MSM_TIMING
      stopTime = CkWallTimer();
      mgr->msmTiming[MsmTimer::COMM] += stopTime - startTime;
#endif
      mgr->msmBlock[bindex.level](
          bindex.n.i, bindex.n.j, bindex.n.k).addCharge(gm);
    }
  }

  void PatchData::addPotential(const Grid<Float>& epart) {
#ifdef MSM_TIMING
    double startTime, stopTime;
    startTime = CkWallTimer();
#endif
    eh += epart;
#ifdef MSM_TIMING
    stopTime = CkWallTimer();
    mgr->msmTiming[MsmTimer::COMM] += stopTime - startTime;
#endif
    if (++cntRecvs == pd->numRecvs) {
      interpolation();
    }
  }

  void PatchData::interpolation() {
#ifdef DEBUG_MSM_GRID
    printf("patchID %d:  interpolation\n", patchID);
#endif

#ifdef MSM_TIMING
    double startTime, stopTime;
    startTime = CkWallTimer();
#endif
#ifndef MSM_COMM_ONLY
    BigReal energy_self = 0;

    Float xphi[ComputeMsmMgr::MAX_POLY_DEGREE+1];
    Float yphi[ComputeMsmMgr::MAX_POLY_DEGREE+1];
    Float zphi[ComputeMsmMgr::MAX_POLY_DEGREE+1];
    Float dxphi[ComputeMsmMgr::MAX_POLY_DEGREE+1];
    Float dyphi[ComputeMsmMgr::MAX_POLY_DEGREE+1];
    Float dzphi[ComputeMsmMgr::MAX_POLY_DEGREE+1];

    const Double rs_edge = Double( mgr->s_edge );
    const int s_size = ComputeMsmMgr::PolyDegree[mgr->approx] + 1;

    const Float hx_1 = Float(mgr->hxlen_1);  // real space inverse grid spacing
    const Float hy_1 = Float(mgr->hylen_1);
    const Float hz_1 = Float(mgr->hzlen_1);

    const int ia = eh.ia();
    const int ib = eh.ib();
    const int ja = eh.ja();
    const int jb = eh.jb();
    const int ka = eh.ka();
    const int kb = eh.kb();
    const int ni = eh.ni();
    const int nj = eh.nj();
    Float *ehbuffer = eh.data().buffer();

    // loop over atoms
    for (int n = 0;  n < coord.len();  n++) {
      Float q = coord[n].charge;
      if (0==q) continue;

      ScaledPosition s = mgr->lattice.scale(coord[n].position);

      BigReal sx_hx = (s.x - mgr->sglower.x) * mgr->shx_1;
      BigReal sy_hy = (s.y - mgr->sglower.y) * mgr->shy_1;
      BigReal sz_hz = (s.z - mgr->sglower.z) * mgr->shz_1;

      BigReal xlo = floor(sx_hx) - rs_edge;
      BigReal ylo = floor(sy_hy) - rs_edge;
      BigReal zlo = floor(sz_hz) - rs_edge;

      // calculate Phi stencils along each dimension
      Float xdelta = Float(sx_hx - xlo);
      mgr->d_stencil_1d(dxphi, xphi, xdelta, hx_1);
      Float ydelta = Float(sy_hy - ylo);
      mgr->d_stencil_1d(dyphi, yphi, ydelta, hy_1);
      Float zdelta = Float(sz_hz - zlo);
      mgr->d_stencil_1d(dzphi, zphi, zdelta, hz_1);

      int ilo = int(xlo);
      int jlo = int(ylo);
      int klo = int(zlo);

#if 0
      // XXX don't need to test twice!

      // test to see if stencil is within edges of grid
      int iswithin = ( ia <= ilo && (ilo+(s_size-1)) <= ib &&
                       ja <= jlo && (jlo+(s_size-1)) <= jb &&
                       ka <= klo && (klo+(s_size-1)) <= kb );

      if ( ! iswithin ) {
        char msg[100];
        snprintf(msg, sizeof(msg), "Atom %d is outside of the MSM grid.",
            coord[n].id);
        NAMD_die(msg);
      }
#endif

      // determine force on atom from surrounding potential grid points
      //Force f = 0;
      //BigReal e = 0;
      Float fx=0, fy=0, fz=0, e=0;
      for (int k = 0;  k < s_size;  k++) {
        int koff = ((k+klo) - ka) * nj;
        for (int j = 0;  j < s_size;  j++) {
          int jkoff = (koff + (j+jlo) - ja) * ni;
          Float cx = yphi[j] * zphi[k];
          Float cy = dyphi[j] * zphi[k];
          Float cz = yphi[j] * dzphi[k];
          for (int i = 0;  i < s_size;  i++) {
            int ijkoff = jkoff + (i+ilo) - ia;
            Float ec = ehbuffer[ijkoff];
            fx += ec * dxphi[i] * cx;
            fy += ec * xphi[i] * cy;
            fz += ec * xphi[i] * cz;
            e += ec * xphi[i] * cx;
          }
        }
      }

#if 0
      force[n].x -= q * (mgr->srx_x * fx + mgr->srx_y * fy + mgr->srx_z * fz);
      force[n].y -= q * (mgr->sry_x * fx + mgr->sry_y * fy + mgr->sry_z * fz);
      force[n].z -= q * (mgr->srz_x * fx + mgr->srz_y * fy + mgr->srz_z * fz);
#endif
      force[n].x -= q * fx;
      force[n].y -= q * fy;
      force[n].z -= q * fz;
      energy += q * e;
      energy_self += q * q;

    } // end loop over atoms

    energy_self *= mgr->gzero;
    energy -= energy_self;
    energy *= 0.5;
#endif // !MSM_COMM_ONLY
#ifdef MSM_TIMING
    stopTime = CkWallTimer();
    mgr->msmTiming[MsmTimer::INTERP] += stopTime - startTime;
    mgr->doneTiming();
#endif
    mgr->doneCompute();
  }

  void PatchData::anterpolationC1Hermite() {
#ifdef DEBUG_MSM_GRID
    printf("patchID %d:  anterpolationC1Hermite\n", patchID);
#endif

#ifdef MSM_TIMING
    double startTime, stopTime;
    startTime = CkWallTimer();
#endif
#ifndef MSM_COMM_ONLY
    Float xphi[2], xpsi[2];
    Float yphi[2], ypsi[2];
    Float zphi[2], zpsi[2];

    const Float hx = Float(mgr->hxlen);  // real space grid spacing
    const Float hy = Float(mgr->hylen);
    const Float hz = Float(mgr->hzlen);

    const int ia = qh_c1hermite.ia();
    const int ib = qh_c1hermite.ib();
    const int ja = qh_c1hermite.ja();
    const int jb = qh_c1hermite.jb();
    const int ka = qh_c1hermite.ka();
    const int kb = qh_c1hermite.kb();
    const int ni = qh_c1hermite.ni();
    const int nj = qh_c1hermite.nj();
    C1Vector *qhbuffer = qh_c1hermite.data().buffer();

    // loop over atoms
    for (int n = 0;  n < coord.len();  n++) {
      Float q = coord[n].charge;
      if (0==q) continue;

      ScaledPosition s = mgr->lattice.scale(coord[n].position);

      BigReal sx_hx = (s.x - mgr->sglower.x) * mgr->shx_1;
      BigReal sy_hy = (s.y - mgr->sglower.y) * mgr->shy_1;
      BigReal sz_hz = (s.z - mgr->sglower.z) * mgr->shz_1;

      BigReal xlo = floor(sx_hx);
      BigReal ylo = floor(sy_hy);
      BigReal zlo = floor(sz_hz);

      // calculate Phi stencils along each dimension
      Float xdelta = Float(sx_hx - xlo);
      mgr->stencil_1d_c1hermite(xphi, xpsi, xdelta, hx);
      Float ydelta = Float(sy_hy - ylo);
      mgr->stencil_1d_c1hermite(yphi, ypsi, ydelta, hy);
      Float zdelta = Float(sz_hz - zlo);
      mgr->stencil_1d_c1hermite(zphi, zpsi, zdelta, hz);

      int ilo = int(xlo);
      int jlo = int(ylo);
      int klo = int(zlo);

      // test to see if stencil is within edges of grid
      int iswithin = ( ia <= ilo && ilo < ib &&
                       ja <= jlo && jlo < jb &&
                       ka <= klo && klo < kb );

      if ( ! iswithin ) {
        char msg[100];
        snprintf(msg, sizeof(msg), "Atom %d is outside of the MSM grid.",
            coord[n].id);
        NAMD_die(msg);
      }

      // determine charge on cube of grid points around atom
      for (int k = 0;  k < 2;  k++) {
        int koff = ((k+klo) - ka) * nj;
        Float c_zphi = zphi[k] * q;
        Float c_zpsi = zpsi[k] * q;
        for (int j = 0;  j < 2;  j++) {
          int jkoff = (koff + (j+jlo) - ja) * ni;
          Float c_yphi_zphi = yphi[j] * c_zphi;
          Float c_ypsi_zphi = ypsi[j] * c_zphi;
          Float c_yphi_zpsi = yphi[j] * c_zpsi;
          Float c_ypsi_zpsi = ypsi[j] * c_zpsi;
          for (int i = 0;  i < 2;  i++) {
            int ijkoff = jkoff + (i+ilo) - ia;
            qhbuffer[ijkoff].velem[D000] += xphi[i] * c_yphi_zphi;
            qhbuffer[ijkoff].velem[D100] += xpsi[i] * c_yphi_zphi;
            qhbuffer[ijkoff].velem[D010] += xphi[i] * c_ypsi_zphi;
            qhbuffer[ijkoff].velem[D001] += xphi[i] * c_yphi_zpsi;
            qhbuffer[ijkoff].velem[D110] += xpsi[i] * c_ypsi_zphi;
            qhbuffer[ijkoff].velem[D101] += xpsi[i] * c_yphi_zpsi;
            qhbuffer[ijkoff].velem[D011] += xphi[i] * c_ypsi_zpsi;
            qhbuffer[ijkoff].velem[D111] += xpsi[i] * c_ypsi_zpsi;
          }
        }
      }

    } // end loop over atoms

#endif // !MSM_COMM_ONLY
#ifdef MSM_TIMING
    stopTime = CkWallTimer();
    mgr->msmTiming[MsmTimer::ANTERP] += stopTime - startTime;
#endif

    sendChargeC1Hermite();
  }

  void PatchData::sendChargeC1Hermite() {
#ifdef MSM_TIMING
    double startTime, stopTime;
#endif
    int priority = 1;
    // buffer portions of grid to send to Blocks on level 0
    for (int n = 0;  n < pd->send.len();  n++) {
#ifdef MSM_TIMING
      startTime = CkWallTimer();
#endif
      // initialize the proper subgrid indexing range
      subgrid_c1hermite.init( pd->send[n].nrange );
      // extract the values from the larger grid into the subgrid
      qh_c1hermite.extract(subgrid_c1hermite);
      // translate the subgrid indexing range to match the MSM block
      subgrid_c1hermite.updateLower( pd->send[n].nrange_wrap.lower() );
      // add the subgrid charges into the block
      BlockIndex& bindex = pd->send[n].nblock_wrap;
      // place subgrid into message
      int msgsz = subgrid_c1hermite.data().len() * sizeof(C1Vector);
      GridMsg *gm = new(msgsz, sizeof(int)) GridMsg;
      SET_PRIORITY(gm, sequence, MSM_PRIORITY + priority);
      gm->put(subgrid_c1hermite, bindex.level, sequence);
#ifdef MSM_TIMING
      stopTime = CkWallTimer();
      mgr->msmTiming[MsmTimer::COMM] += stopTime - startTime;
#endif
      mgr->msmC1HermiteBlock[bindex.level](
          bindex.n.i, bindex.n.j, bindex.n.k).addCharge(gm);
    }
  }

  void PatchData::addPotentialC1Hermite(const Grid<C1Vector>& epart) {
#ifdef MSM_TIMING
    double startTime, stopTime;
    startTime = CkWallTimer();
#endif
    eh_c1hermite += epart;
#ifdef MSM_TIMING
    stopTime = CkWallTimer();
    mgr->msmTiming[MsmTimer::COMM] += stopTime - startTime;
#endif
    if (++cntRecvs == pd->numRecvs) {
      interpolationC1Hermite();
    }
  }

  void PatchData::interpolationC1Hermite() {
#ifdef DEBUG_MSM_GRID
    printf("patchID %d:  interpolation\n", patchID);
#endif

#ifdef MSM_TIMING
    double startTime, stopTime;
    startTime = CkWallTimer();
#endif
#ifndef MSM_COMM_ONLY
    BigReal energy_self = 0;

    Float xphi[2], dxphi[2], xpsi[2], dxpsi[2];
    Float yphi[2], dyphi[2], ypsi[2], dypsi[2];
    Float zphi[2], dzphi[2], zpsi[2], dzpsi[2];

    const Float hx = Float(mgr->hxlen);      // real space grid spacing
    const Float hy = Float(mgr->hylen);
    const Float hz = Float(mgr->hzlen);

    const Float hx_1 = Float(mgr->hxlen_1);  // real space inverse grid spacing
    const Float hy_1 = Float(mgr->hylen_1);
    const Float hz_1 = Float(mgr->hzlen_1);

    const int ia = eh_c1hermite.ia();
    const int ib = eh_c1hermite.ib();
    const int ja = eh_c1hermite.ja();
    const int jb = eh_c1hermite.jb();
    const int ka = eh_c1hermite.ka();
    const int kb = eh_c1hermite.kb();
    const int ni = eh_c1hermite.ni();
    const int nj = eh_c1hermite.nj();
    C1Vector *ehbuffer = eh_c1hermite.data().buffer();

    // loop over atoms
    for (int n = 0;  n < coord.len();  n++) {
      Float q = coord[n].charge;
      if (0==q) continue;

      ScaledPosition s = mgr->lattice.scale(coord[n].position);

      BigReal sx_hx = (s.x - mgr->sglower.x) * mgr->shx_1;
      BigReal sy_hy = (s.y - mgr->sglower.y) * mgr->shy_1;
      BigReal sz_hz = (s.z - mgr->sglower.z) * mgr->shz_1;

      BigReal xlo = floor(sx_hx);
      BigReal ylo = floor(sy_hy);
      BigReal zlo = floor(sz_hz);

      // calculate Phi stencils along each dimension
      Float xdelta = Float(sx_hx - xlo);
      mgr->d_stencil_1d_c1hermite(dxphi, xphi, dxpsi, xpsi,
          xdelta, hx, hx_1);
      Float ydelta = Float(sy_hy - ylo);
      mgr->d_stencil_1d_c1hermite(dyphi, yphi, dypsi, ypsi,
          ydelta, hy, hy_1);
      Float zdelta = Float(sz_hz - zlo);
      mgr->d_stencil_1d_c1hermite(dzphi, zphi, dzpsi, zpsi,
          zdelta, hz, hz_1);

      int ilo = int(xlo);
      int jlo = int(ylo);
      int klo = int(zlo);

#if 0
      // XXX don't need to test twice!

      // test to see if stencil is within edges of grid
      int iswithin = ( ia <= ilo && ilo < ib &&
                       ja <= jlo && jlo < jb &&
                       ka <= klo && klo < kb );

      if ( ! iswithin ) {
        char msg[100];
        snprintf(msg, sizeof(msg), "Atom %d is outside of the MSM grid.",
            coord[n].id);
        NAMD_die(msg);
      }
#endif

      // determine force on atom from surrounding potential grid points
      Float fx=0, fy=0, fz=0, e=0;
      for (int k = 0;  k < 2;  k++) {
        int koff = ((k+klo) - ka) * nj;
        for (int j = 0;  j < 2;  j++) {
          int jkoff = (koff + (j+jlo) - ja) * ni;
          Float c_yphi_zphi = yphi[j] * zphi[k];
          Float c_ypsi_zphi = ypsi[j] * zphi[k];
          Float c_yphi_zpsi = yphi[j] * zpsi[k];
          Float c_ypsi_zpsi = ypsi[j] * zpsi[k];
          Float c_yphi_dzphi = yphi[j] * dzphi[k];
          Float c_ypsi_dzphi = ypsi[j] * dzphi[k];
          Float c_yphi_dzpsi = yphi[j] * dzpsi[k];
          Float c_ypsi_dzpsi = ypsi[j] * dzpsi[k];
          Float c_dyphi_zphi = dyphi[j] * zphi[k];
          Float c_dypsi_zphi = dypsi[j] * zphi[k];
          Float c_dyphi_zpsi = dyphi[j] * zpsi[k];
          Float c_dypsi_zpsi = dypsi[j] * zpsi[k];
          for (int i = 0;  i < 2;  i++) {
            int ijkoff = jkoff + (i+ilo) - ia;
            fx += dxphi[i] * (c_yphi_zphi * ehbuffer[ijkoff].velem[D000]
                + c_ypsi_zphi * ehbuffer[ijkoff].velem[D010]
                + c_yphi_zpsi * ehbuffer[ijkoff].velem[D001]
                + c_ypsi_zpsi * ehbuffer[ijkoff].velem[D011])
              + dxpsi[i] * (c_yphi_zphi * ehbuffer[ijkoff].velem[D100]
                  + c_ypsi_zphi * ehbuffer[ijkoff].velem[D110]
                  + c_yphi_zpsi * ehbuffer[ijkoff].velem[D101]
                  + c_ypsi_zpsi * ehbuffer[ijkoff].velem[D111]);
            fy += xphi[i] * (c_dyphi_zphi * ehbuffer[ijkoff].velem[D000]
                + c_dypsi_zphi * ehbuffer[ijkoff].velem[D010]
                + c_dyphi_zpsi * ehbuffer[ijkoff].velem[D001]
                + c_dypsi_zpsi * ehbuffer[ijkoff].velem[D011])
              + xpsi[i] * (c_dyphi_zphi * ehbuffer[ijkoff].velem[D100]
                  + c_dypsi_zphi * ehbuffer[ijkoff].velem[D110]
                  + c_dyphi_zpsi * ehbuffer[ijkoff].velem[D101]
                  + c_dypsi_zpsi * ehbuffer[ijkoff].velem[D111]);
            fz += xphi[i] * (c_yphi_dzphi * ehbuffer[ijkoff].velem[D000]
                + c_ypsi_dzphi * ehbuffer[ijkoff].velem[D010]
                + c_yphi_dzpsi * ehbuffer[ijkoff].velem[D001]
                + c_ypsi_dzpsi * ehbuffer[ijkoff].velem[D011])
              + xpsi[i] * (c_yphi_dzphi * ehbuffer[ijkoff].velem[D100]
                  + c_ypsi_dzphi * ehbuffer[ijkoff].velem[D110]
                  + c_yphi_dzpsi * ehbuffer[ijkoff].velem[D101]
                  + c_ypsi_dzpsi * ehbuffer[ijkoff].velem[D111]);
            e += xphi[i] * (c_yphi_zphi * ehbuffer[ijkoff].velem[D000]
                + c_ypsi_zphi * ehbuffer[ijkoff].velem[D010]
                + c_yphi_zpsi * ehbuffer[ijkoff].velem[D001]
                + c_ypsi_zpsi * ehbuffer[ijkoff].velem[D011])
              + xpsi[i] * (c_yphi_zphi * ehbuffer[ijkoff].velem[D100]
                  + c_ypsi_zphi * ehbuffer[ijkoff].velem[D110]
                  + c_yphi_zpsi * ehbuffer[ijkoff].velem[D101]
                  + c_ypsi_zpsi * ehbuffer[ijkoff].velem[D111]);
          }
        }
      }

#if 0
      force[n].x -= q * (mgr->srx_x * fx + mgr->srx_y * fy + mgr->srx_z * fz);
      force[n].y -= q * (mgr->sry_x * fx + mgr->sry_y * fy + mgr->sry_z * fz);
      force[n].z -= q * (mgr->srz_x * fx + mgr->srz_y * fy + mgr->srz_z * fz);
#endif
      force[n].x -= q * fx;
      force[n].y -= q * fy;
      force[n].z -= q * fz;
      energy += q * e;
      energy_self += q * q;

    } // end loop over atoms

    energy_self *= mgr->gzero;
    energy -= energy_self;
    energy *= 0.5;
#endif // !MSM_COMM_ONLY
#ifdef MSM_TIMING
    stopTime = CkWallTimer();
    mgr->msmTiming[MsmTimer::INTERP] += stopTime - startTime;
    mgr->doneTiming();
#endif
    mgr->doneCompute();
  }

} // namespace msm


#include "ComputeMsmMgr.def.h"
