/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifdef NAMD_FFTW
#ifdef NAMD_FFTW_3
#include <fftw3.h>
#else
#ifdef NAMD_FFTW_NO_TYPE_PREFIX
#include <fftw.h>
#include <rfftw.h>
#else
#include <sfftw.h>
#include <srfftw.h>
#endif
#endif
#endif

#include <assert.h>

#include "InfoStream.h"
#include "Node.h"
#include "PatchMap.h"
#include "PatchMap.inl"
#include "AtomMap.h"
#include "OptPme.h"
#include "OptPmeMgr.decl.h"
#include "OptPmeRealSpace.h"
#include "PmeKSpace.h"
#include "ComputeNonbondedUtil.h"
#include "PatchMgr.h"
#include "Molecule.h"
#include "ReductionMgr.h"
//#include "ComputeMgr.h"
#include "ComputeMgr.decl.h"
// #define DEBUGM
#define MIN_DEBUG_LEVEL 3
#include "Debug.h"
#include "SimParameters.h"
#include "WorkDistrib.h"
#include "varsizemsg.h"
#include "Random.h"
#include "Priorities.h"
#include "PmeBase.inl"

extern char *pencilPMEProcessors;

#include "fftlib.h"
#include "fftmap.h"

//Very large integer
int many_to_many_start = 0x7fffffff;

class OptPmeMgr : public CBase_OptPmeMgr {
public:
  friend class OptPmeCompute;
  OptPmeMgr();
  ~OptPmeMgr();

  void initialize(CkQdMsg*);
  void initialize_pencils(CkQdMsg*);
  void activate_pencils(CkQdMsg*);
  void recvArrays(CProxy_OptPmeXPencil, CProxy_OptPmeYPencil, CProxy_OptPmeZPencil);

  void recvUngrid(OptPmeGridMsg *);
  void ungridCalc(OptPmeDummyMsg *);
  void ungridCalc_subcompute(OptPmeSubComputeMsg *);
  void ungridCalc_subcompute_done(OptPmeSubComputeMsg *);
  void doWorkOnPeer(OptPmeSubComputeMsg *);
  void recvEvir (CkReductionMsg *msg);

  void setCompute(OptPmeCompute *c) { pmeCompute = c; c->setMgr(this); }

private:
  CProxy_OptPmeMgr pmeProxy;
  CProxy_OptPmeMgr pmeProxyDir;
  OptPmeCompute *pmeCompute;
  PmeGrid myGrid;
  PmeKSpace *myKSpace;

  CProxy_OptPmeXPencil xPencil;
  CProxy_OptPmeYPencil yPencil;
  CProxy_OptPmeZPencil zPencil;
  int    numPencilsActive;
  int    ungrid_count;
  SubmitReduction *reduction;
  int    _iter;
  void   *handle;
  bool   constant_pressure;    //Does the simulation need constant pressure
  int    subcompute_count;

  int peersAllocated;
  int peers [SUBCOMPUTE_NPAR];
  OptPmeSubComputeMsg *subcompute_msgs[SUBCOMPUTE_NPAR];
};


void pme_f2d (double *dst, float *src, int N);
void pme_d2f (float *dst, double *src, int N);

static inline void initializePmeGrid (SimParameters *simParams, PmeGrid &grid) {
    int xBlocks = 0, yBlocks = 0, zBlocks= 0;

    if ( simParams->PMEPencils > 1 ) {
      xBlocks = yBlocks = zBlocks = simParams->PMEPencils;
    } else {
      int nb2 = ( simParams->PMEGridSizeX * simParams->PMEGridSizeY
		  * simParams->PMEGridSizeZ ) / simParams->PMEMinPoints;
      if ( nb2 > CkNumPes() ) nb2 = CkNumPes();
      if ( nb2 < 1 ) nb2 = 1;
      int nb = (int) sqrt((float)nb2);
      if ( nb < 1 ) nb = 1;
      xBlocks = zBlocks = nb;
      yBlocks = nb2 / nb;
    }
    
    int dimx = simParams->PMEGridSizeX;
    int bx = 1 + ( dimx - 1 ) / xBlocks;
    xBlocks = 1 + ( dimx - 1 ) / bx;
    
    int dimy = simParams->PMEGridSizeY;
    int by = 1 + ( dimy - 1 ) / yBlocks;
    yBlocks = 1 + ( dimy - 1 ) / by;
    
    int dimz = simParams->PMEGridSizeZ / 2 + 1;  // complex
    int bz = 1 + ( dimz - 1 ) / zBlocks;
    zBlocks = 1 + ( dimz - 1 ) / bz;

    grid.xBlocks = xBlocks;
    grid.yBlocks = yBlocks;
    grid.zBlocks = zBlocks;

    grid.K1 = simParams->PMEGridSizeX;
    grid.K2 = simParams->PMEGridSizeY;
    grid.K3 = simParams->PMEGridSizeZ;
    grid.order = simParams->PMEInterpOrder;
    grid.dim2 = grid.K2;
    grid.dim3 = 2 * (grid.K3/2 + 1);

    grid.block1 = ( grid.K1 + xBlocks - 1 ) / xBlocks;
    grid.block2 = ( grid.K2 + yBlocks - 1 ) / yBlocks;
    grid.block3 = ( grid.K3/2 + 1 + zBlocks - 1 ) / zBlocks;  // complex
}

#ifdef NAMD_FFTW
static CmiNodeLock fftw_plan_lock;
#endif


static inline void scale_n_copy_coordinates(CompAtom *x, PmeParticle p[], 
					    int &N, 
					    Lattice &lattice, PmeGrid grid,
					    //double **qline,
					    double xmin, double xlen,
					    double ymin, double ylen,
					    double zmin, double zlen,
					    int &scg) {
  Vector origin = lattice.origin();
  Vector recip1 = lattice.a_r();
  Vector recip2 = lattice.b_r();
  Vector recip3 = lattice.c_r();
  double ox = origin.x;
  double oy = origin.y;
  double oz = origin.z;
  double r1x = recip1.x;
  double r1y = recip1.y;
  double r1z = recip1.z;
  double r2x = recip2.x;
  double r2y = recip2.y;
  double r2z = recip2.z;
  double r3x = recip3.x;
  double r3y = recip3.y;
  double r3z = recip3.z;
  int K1 = grid.K1;
  int K2 = grid.K2;
  int K3 = grid.K3;

  const BigReal coulomb_sqrt = sqrt( COULOMB * ComputeNonbondedUtil::scaling
				     * ComputeNonbondedUtil::dielectric_1 );


  int natoms = 0;
  for (int i=0; i<N; i++) {
    double px = x[i].position.x - ox;
    double py = x[i].position.y - oy;
    double pz = x[i].position.z - oz;
    double sx = px*r1x + py*r1y + pz*r1z;
    double sy = px*r2x + py*r2y + pz*r2z;
    double sz = px*r3x + py*r3y + pz*r3z;
    p[natoms].x = K1 * ( sx - floor(sx) );
    p[natoms].y = K2 * ( sy - floor(sy) );
    p[natoms].z = K3 * ( sz - floor(sz) );
#ifndef ARCH_POWERPC
    //  Check for rare rounding condition where K * ( 1 - epsilon ) == K      
    //  which was observed with g++ on Intel x86 architecture. 
    if ( p[natoms].x == K1 ) p[natoms].x = 0;
    if ( p[natoms].y == K2 ) p[natoms].y = 0;
    if ( p[natoms].z == K3 ) p[natoms].z = 0;
#endif

#if 1 //stray charge detection
    BigReal u1,u2,u3;
    u1 = (int) (p[natoms].x - xmin);
    if (u1 >= grid.K1) u1 -= grid.K1;    
    u2 = (int) (p[natoms].y - ymin);
    if (u2 >= grid.K2) u2 -= grid.K2;    
    u3 = (int) (p[natoms].z - zmin);
    if (u3 >= grid.K3) u3 -= grid.K3;
    
    if ( (u1 < 0.0) || (u1 >= xlen) || 
	 (u2 < 0.0) || (u2 >= ylen) || 
	 (u3 < 0.0) || (u3 >= zlen) ) {
      scg ++;
      continue;
    }
#endif
    
    p[natoms].cg = coulomb_sqrt * x[i].charge;
    natoms ++;
  }
  N = natoms;
}


OptPmeMgr::OptPmeMgr() : pmeProxy(thisgroup), 
				 pmeProxyDir(thisgroup), pmeCompute(0) {

  CkpvAccess(BOCclass_group).computePmeMgr = thisgroup;

  myKSpace = 0;
  ungrid_count = 0;
  peersAllocated = 0;

#ifdef NAMD_FFTW
  if ( CmiMyRank() == 0 ) {
    fftw_plan_lock = CmiCreateLock();
  }
#endif    
}


void OptPmeMgr::recvArrays(CProxy_OptPmeXPencil x, CProxy_OptPmeYPencil y, CProxy_OptPmeZPencil z) {
  xPencil = x;  yPencil = y;  zPencil = z;
}

void OptPmeMgr::initialize(CkQdMsg *msg) {
    delete msg;

    _iter = 0;

    handle = CmiDirect_manytomany_allocate_handle ();

    SimParameters *simParams = Node::Object()->simParameters;
    PatchMap *patchMap = PatchMap::Object();
    
    initializePmeGrid (simParams, myGrid);    

    if (simParams->langevinPistonOn || simParams->berendsenPressureOn)	
      constant_pressure = true;
    else
      constant_pressure = false;      

    bool useManyToMany = simParams->useManyToMany;
    //Many-to-many requires that patches and pmepencils are all on different processors
    //int npes = patchMap->numPatches() + 
    //         myGrid.xBlocks *  myGrid.yBlocks + 
    //         myGrid.zBlocks *  myGrid.xBlocks +
    //         myGrid.yBlocks *  myGrid.zBlocks;
    
    int npes = patchMap->numPatches();
    if (npes < myGrid.xBlocks *  myGrid.yBlocks)
      npes = myGrid.xBlocks *  myGrid.yBlocks;
    if (npes <  myGrid.zBlocks *  myGrid.xBlocks)
      npes = myGrid.zBlocks *  myGrid.xBlocks;
    if (npes < myGrid.yBlocks *  myGrid.zBlocks)
      npes = myGrid.yBlocks *  myGrid.zBlocks;
    
   if (npes >= CkNumPes()) {
      if (CkMyPe() == 0)
	printf ("Warning : Not enough processors for the many-to-many optimization \n");      
      useManyToMany = false;
    }
    
    if (useManyToMany)  {
      if (CkMyPe() == 0)
	printf ("Enabling the Many-to-many optimization\n");
      //defaults to max integer
      many_to_many_start = MANY_TO_MANY_START;
    }

    if (CkMyRank() == 0) { //create the pencil pme processor map
      pencilPMEProcessors = new char [CkNumPes()];
      memset (pencilPMEProcessors, 0, sizeof(char) * CkNumPes());
    }

    if ( CkMyPe() == 0) {
      iout << iINFO << "PME using " << myGrid.xBlocks << " x " <<
        myGrid.yBlocks << " x " << myGrid.zBlocks <<
        " pencil grid for FFT and reciprocal sum.\n" << endi;
      
      CProxy_OptPmePencilMapZ   mapz;      
      CProxy_OptPmePencilMapY   mapy;
      CProxy_OptPmePencilMapX   mapx;
      
      mapz = CProxy_OptPmePencilMapZ::ckNew(myGrid.xBlocks, myGrid.yBlocks, myGrid.zBlocks);      
      mapy = CProxy_OptPmePencilMapY::ckNew(myGrid.xBlocks, myGrid.yBlocks, myGrid.zBlocks);
      mapx = CProxy_OptPmePencilMapX::ckNew(myGrid.xBlocks, myGrid.yBlocks, myGrid.zBlocks);
      
      CkArrayOptions optsz;
      optsz.setMap (mapz);
      CkArrayOptions optsy;
      optsy.setMap (mapy);
      CkArrayOptions optsx;
      optsx.setMap (mapx);
      
      zPencil = CProxy_OptPmeZPencil::ckNew(optsz);  
      yPencil = CProxy_OptPmeYPencil::ckNew(optsy);  
      xPencil = CProxy_OptPmeXPencil::ckNew(optsx);  
      
      int x,y,z;
      for (x = 0; x < myGrid.xBlocks; ++x)
	for (y = 0; y < myGrid.yBlocks; ++y ) {
	  zPencil(x,y,0).insert();
	}
      zPencil.doneInserting();
      
      for (z = 0; z < myGrid.zBlocks; ++z )
	for (x = 0; x < myGrid.xBlocks; ++x ) {
	  yPencil(x,0,z).insert();
	}
      yPencil.doneInserting();
      
      for (y = 0; y < myGrid.yBlocks; ++y )	
	for (z = 0; z < myGrid.zBlocks; ++z ) {
	  xPencil(0,y,z).insert();
	}
      xPencil.doneInserting();      
      
      pmeProxy.recvArrays(xPencil,yPencil,zPencil);
      OptPmePencilInitMsgData msgdata;
      msgdata.grid = myGrid;
      msgdata.xBlocks = myGrid.xBlocks;
      msgdata.yBlocks = myGrid.yBlocks;
      msgdata.zBlocks = myGrid.zBlocks;
      msgdata.xPencil = xPencil;
      msgdata.yPencil = yPencil;
      msgdata.zPencil = zPencil;
      msgdata.constant_pressure = constant_pressure;

      CkCallback cb (CkIndex_OptPmeMgr::recvEvir(NULL), thisProxy[0]);
      msgdata.cb_energy = cb;

      msgdata.pmeProxy = pmeProxyDir;
      xPencil.init(new OptPmePencilInitMsg(msgdata));
      yPencil.init(new OptPmePencilInitMsg(msgdata));
      zPencil.init(new OptPmePencilInitMsg(msgdata));
     
#if 0 
      reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);
#endif

#ifndef NAMD_FFTW
      NAMD_die("Sorry, FFTW must be compiled in to use PME.");
#endif
    }

}

void OptPmeMgr::initialize_pencils(CkQdMsg *msg) {
  delete msg;

  SimParameters *simParams = Node::Object()->simParameters;

  PatchMap *patchMap = PatchMap::Object();
  Lattice lattice = simParams->lattice;
  BigReal sysdima = lattice.a_r().unit() * lattice.a();
  BigReal sysdimb = lattice.b_r().unit() * lattice.b();
  BigReal cutoff = simParams->cutoff;
  BigReal patchdim = simParams->patchDimension;
  int numPatches = patchMap->numPatches();

  //fprintf(stderr, "Node %d PE %d trying to allocate %d bytes\n", CmiMyNode(), CmiMyPe(), myGrid.xBlocks*myGrid.yBlocks);

  char *pencilActive = new char[myGrid.xBlocks*myGrid.yBlocks];
  for ( int i=0; i<myGrid.xBlocks; ++i ) {
    for ( int j=0; j<myGrid.yBlocks; ++j ) {
      pencilActive[i*myGrid.yBlocks+j] = 0;
    }
  }

  //Right now we only support one patch per processor
  assert (patchMap->numPatchesOnNode(CkMyPe()) <= 1);
  for ( int pid=0; pid < numPatches; ++pid ) {
    int pnode = patchMap->node(pid);
    if ( pnode != CkMyPe() ) continue;

    BigReal minx = patchMap->min_a(pid);
    BigReal maxx = patchMap->max_a(pid);
    BigReal margina = 0.5 * ( patchdim - cutoff ) / sysdima;
    // min1 (max1) is smallest (largest) grid line for this patch
    int min1 = ((int) floor(myGrid.K1 * (minx - margina))) - myGrid.order + 1;
    int max1 = ((int) floor(myGrid.K1 * (maxx + margina)));

    BigReal miny = patchMap->min_b(pid);
    BigReal maxy = patchMap->max_b(pid);
    BigReal marginb = 0.5 * ( patchdim - cutoff ) / sysdimb;
    // min2 (max2) is smallest (largest) grid line for this patch
    int min2 = ((int) floor(myGrid.K2 * (miny - marginb))) - myGrid.order + 1;
    int max2 = ((int) floor(myGrid.K2 * (maxy + marginb)));

    for ( int i=min1; i<=max1; ++i ) {
      int ix = i;
      while ( ix >= myGrid.K1 ) ix -= myGrid.K1;
      while ( ix < 0 ) ix += myGrid.K1;
      for ( int j=min2; j<=max2; ++j ) {
        int jy = j;
        while ( jy >= myGrid.K2 ) jy -= myGrid.K2;
        while ( jy < 0 ) jy += myGrid.K2;
        pencilActive[(ix / myGrid.block1)*myGrid.yBlocks + (jy / myGrid.block2)] = 1;
      }
    }
  }

  numPencilsActive = 0;
  for ( int i=0; i<myGrid.xBlocks; ++i ) {
    for ( int j=0; j<myGrid.yBlocks; ++j ) {
      if ( pencilActive[i*myGrid.yBlocks+j] ) {
        ++numPencilsActive;

        zPencil(i,j,0).dummyRecvGrid(CkMyPe(),0);
      }
    }
  }

  ungrid_count = numPencilsActive;
  delete [] pencilActive;  
}


void OptPmeMgr::activate_pencils(CkQdMsg *msg) {
  if ( CkMyPe() == 0 ) zPencil.dummyRecvGrid(CkMyPe(),1);
}


OptPmeMgr::~OptPmeMgr() {
  delete myKSpace;
}

void OptPmeMgr::recvUngrid(OptPmeGridMsg *msg) {
  if ( ungrid_count == 0 ) {
    NAMD_bug("Message order failure in OptPmeMgr::recvUngrid\n");
  }
    
  pmeCompute->copyPencils(msg);
  delete msg;
  --ungrid_count;

  if ( ungrid_count == 0 ) {
    //CkPrintf("recvUngrid on Pe(%d)\n",CkMyPe());
    ungridCalc(NULL);
  }
}

void OptPmeMgr::ungridCalc(OptPmeDummyMsg *dmsg) {    
  pmeCompute->ungridForces_init();
  if ( CmiMyNodeSize() >= SUBCOMPUTE_NPAR ) {
    int npar = SUBCOMPUTE_NPAR;
    OptPmeSubComputeMsg *smsg = NULL;

    if (!peersAllocated) {
      peersAllocated = 1;
      int next_rank = CmiMyRank();   
      PatchMap *patchMap = PatchMap::Object();          

      for (int i = 1; i < npar; ++i) {      
	smsg = new (PRIORITY_SIZE) OptPmeSubComputeMsg;
	subcompute_msgs[i] = smsg;
	smsg->src_pe  = CkMyPe();
	smsg->compute = pmeCompute;

	next_rank ++;
	if (next_rank >= CmiMyNodeSize())
	  next_rank = 0;
	int n = 0;
	int nr = next_rank;
	while(n < CmiMyNodeSize() &&
	      patchMap->numPatchesOnNode(CmiNodeFirst(CmiMyNode())+nr) > 0)
	{
	  nr ++;
	  if (nr >= CmiMyNodeSize())
	    nr = 0;
	  n++;
	}
	if (n < CmiMyNodeSize()) 
	  next_rank = nr;  //we are successful, so save this rank
	
	smsg->dest = next_rank;
      }

      //Local subcompute msg
      smsg = new (PRIORITY_SIZE) OptPmeSubComputeMsg;
      subcompute_msgs[0] = smsg;
      smsg->src_pe  = CkMyPe();
      smsg->compute = pmeCompute;
      smsg->dest    = CmiMyRank();
    }

    int start  = 0;
    int nlocal = pmeCompute->getNumLocalAtoms();
    //CmiAssert (npar <= nlocal);
    if (nlocal < npar)
      npar = nlocal;
    if (npar == 0)
      npar = 1;
    int n_per_iter = nlocal / npar;
    //We dont handle the case where there are very few atoms
    subcompute_count = npar;

    for (int i = 0; i < npar; ++i) {      
      smsg = subcompute_msgs[i];
      smsg->start   = start;
      smsg->end     = start + n_per_iter;
      start += n_per_iter;
      if (i == npar - 1)
	smsg->end = nlocal;
      pmeProxy[CmiNodeFirst(CmiMyNode())+smsg->dest].ungridCalc_subcompute(smsg);
    }    
  }
  else {
    pmeCompute->ungridForces_compute(0, 0);
    pmeCompute->ungridForces_finalize();
    ungrid_count = numPencilsActive; 
  }
}

void OptPmeMgr::ungridCalc_subcompute(OptPmeSubComputeMsg *msg){ 
  OptPmeCompute *compute = (OptPmeCompute *) msg->compute;
  compute->ungridForces_compute(msg->start, msg->end);
  pmeProxy[msg->src_pe].ungridCalc_subcompute_done(msg);
}

void OptPmeMgr::ungridCalc_subcompute_done(OptPmeSubComputeMsg *msg){  
  subcompute_count --;
  //delete msg; //message pointers saved
  if (subcompute_count == 0) {
    pmeCompute->ungridForces_finalize();
    ungrid_count = numPencilsActive; 
  }
}

void OptPmeMgr::doWorkOnPeer(OptPmeSubComputeMsg *msg) {
  OptPmeCompute *compute = (OptPmeCompute *) msg->compute;
  compute->doWorkOnPeer();
  //  delete msg; //saved in compute
}

void OptPmeMgr::recvEvir (CkReductionMsg *msg) {

  assert (CkMyPe() == 0);

  double *data = (double *) msg->getData();
  assert (msg->getSize() == 7 * sizeof(double));

  //printf ("[%d]: Received Evir\n", CkMyPe());

  double scale = 1.;
  reduction->item(REDUCTION_ELECT_ENERGY_SLOW) += data[0] * scale;
  reduction->item(REDUCTION_VIRIAL_SLOW_XX) += data[1] * scale;
  reduction->item(REDUCTION_VIRIAL_SLOW_XY) += data[2] * scale;
  reduction->item(REDUCTION_VIRIAL_SLOW_XZ) += data[3] * scale;
  reduction->item(REDUCTION_VIRIAL_SLOW_YX) += data[2] * scale;
  reduction->item(REDUCTION_VIRIAL_SLOW_YY) += data[4] * scale;
  reduction->item(REDUCTION_VIRIAL_SLOW_YZ) += data[5] * scale;
  reduction->item(REDUCTION_VIRIAL_SLOW_ZX) += data[3] * scale;
  reduction->item(REDUCTION_VIRIAL_SLOW_ZY) += data[5] * scale;
  reduction->item(REDUCTION_VIRIAL_SLOW_ZZ) += data[6] * scale;   

  delete msg;

  SimParameters *simParams = Node::Object()->simParameters;
  int fef = simParams->fullElectFrequency;
  for (int i = 0; i < fef; i++) {
    reduction->submit();
  }
}

OptPmeCompute::OptPmeCompute(ComputeID c) :
  ComputeHomePatches(c)
{
  DebugM(4,"OptPmeCompute created.\n");

  CProxy_OptPmeMgr::ckLocalBranch(
	CkpvAccess(BOCclass_group).computePmeMgr)->setCompute(this);

  _initialized = false;

  useAvgPositions = 1;
  localResults = NULL;

  reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);
  SimParameters *simParams = Node::Object()->simParameters;
}

void recv_ungrid_done (void *m) {
  OptPmeDummyMsg *msg =  (OptPmeDummyMsg *) m;  
  CProxy_OptPmeMgr pmeProxy (CkpvAccess(BOCclass_group).computePmeMgr);
  pmeProxy[msg->to_pe].ungridCalc (msg);
}  


void OptPmeCompute::resetPatchCoordinates (const Lattice &lattice) {
  PatchMap *patchMap = PatchMap::Object();    
  SimParameters *simParams = Node::Object()->simParameters;

  BigReal sysdima = lattice.a_r().unit() * lattice.a();
  BigReal sysdimb = lattice.b_r().unit() * lattice.b();
  BigReal sysdimc = lattice.c_r().unit() * lattice.c();
  BigReal cutoff = simParams->cutoff;
  BigReal patchdim = simParams->patchDimension;

  int pid = patchList[0].patchID; 
  assert (patchList.size() == 1);

  //printf ("Patch[%d]: zstart = %d, zlen = %d, maxz %f, minz %f, marginec %f\n", pid, zstart, zlen, maxz, minz, marginc);

  BigReal minx = patchMap->min_a(pid);
  BigReal maxx = patchMap->max_a(pid);
  BigReal margina = 0.5 * ( patchdim - cutoff ) / sysdima;
  // min1 (max1) is smallest (largest) grid line for this patch
  int min1 = ((int) floor(myGrid.K1 * (minx - margina))) - myGrid.order + 1;
  int max1 = ((int) floor(myGrid.K1 * (maxx + margina)));

  BigReal miny = patchMap->min_b(pid);
  BigReal maxy = patchMap->max_b(pid);
  BigReal marginb = 0.5 * ( patchdim - cutoff ) / sysdimb;
  // min2 (max2) is smallest (largest) grid line for this patch
  int min2 = ((int) floor(myGrid.K2 * (miny - marginb))) - myGrid.order + 1;
  int max2 = ((int) floor(myGrid.K2 * (maxy + marginb)));

  BigReal minz = patchMap->min_c(pid);
  BigReal maxz = patchMap->max_c(pid);
  BigReal marginc = 0.5 * ( patchdim - cutoff ) / sysdimc;
  // min3 (max3) is smallest (largest) grid line for this patch
  int min3 = ((int) floor(myGrid.K3 * (minz - marginc))) - myGrid.order + 1;
  int max3 = ((int) floor(myGrid.K3 * (maxz + marginc)));

  xlen = max1 - min1 + 1;
  xstart = min1;  
  
  ylen = max2 - min2 + 1;
  ystart = min2;  
  
  zlen = max3 - min3 + 1;
  zstart = min3;  
}

void OptPmeCompute::initializeOptPmeCompute () {
  
  _initialized = true;
  
  strayChargeErrors = 0;

  SimParameters *simParams = Node::Object()->simParameters;
  PatchMap *patchMap = PatchMap::Object();

  initializePmeGrid (simParams, myGrid);

  alchFepOn    = simParams->alchFepOn;
  alchThermIntOn = simParams->alchThermIntOn;
  lesOn    = simParams->lesOn;
  pairOn   = simParams->pairInteractionOn;

  assert (!alchFepOn);
  assert (!alchThermIntOn);
  assert (!lesOn);
  assert (!pairOn);
  
  qsize = myGrid.K1 * myGrid.dim2 * myGrid.dim3;
  fsize = myGrid.K1 * myGrid.dim2;
  q_arr = new double*[fsize];
  memset( (void*) q_arr, 0, fsize * sizeof(double*) );

  assert (myMgr != NULL);

  const Lattice lattice = simParams->lattice;
  resetPatchCoordinates (lattice);
  
  PencilElement *pencilActive = new PencilElement [myGrid.xBlocks * myGrid.yBlocks];
  memset (pencilActive, 0, myGrid.xBlocks * myGrid.yBlocks * sizeof(PencilElement));

  for ( int i=xstart; i<=xstart+xlen-1; ++i ) {
    int ix = i;
    while ( ix >= myGrid.K1 ) ix -= myGrid.K1;
    while ( ix < 0 ) ix += myGrid.K1;
    for ( int j=ystart; j<=ystart+ylen-1; ++j ) {
      int jy = j;
      while ( jy >= myGrid.K2 ) jy -= myGrid.K2;
      while ( jy < 0 ) jy += myGrid.K2;
      
      int pencil_idx = (ix / myGrid.block1)*myGrid.yBlocks + (jy / myGrid.block2);
      //If not initialized yet, initialize this pencil
      if (! pencilActive[pencil_idx].isActive) {
	pencilActive[pencil_idx].isActive = 1;
	pencilActive[pencil_idx].xmin = ix;
	pencilActive[pencil_idx].xmax = ix;
	pencilActive[pencil_idx].ymin = jy;
	pencilActive[pencil_idx].ymax = jy;
	pencilActive[pencil_idx].ib = ix / myGrid.block1;
	pencilActive[pencil_idx].jb = jy / myGrid.block2;
      }
      else { //pencil has been initialized, update the min and max
	if (pencilActive[pencil_idx].xmin > ix)
	  pencilActive[pencil_idx].xmin = ix;
	if (pencilActive[pencil_idx].xmax < ix)
	  pencilActive[pencil_idx].xmax = ix;	
	if (pencilActive[pencil_idx].ymin > jy)
	  pencilActive[pencil_idx].ymin = jy;
	if (pencilActive[pencil_idx].ymax < jy)
	  pencilActive[pencil_idx].ymax = jy;
      }
    }
  }
  
  int nactive = 0; 
  for (int ib=0; ib<myGrid.xBlocks; ++ib) 
    for (int jb=0; jb<myGrid.yBlocks; ++jb) 
      if (pencilActive[ib * myGrid.yBlocks + jb].isActive) 
	nactive ++;

  assert (nactive == myMgr->numPencilsActive);  
  
  nzlines = xlen * ylen;
  zline_storage = new double [nzlines * zlen];
  sp_zstorage = new float [nzlines * zlen];  

  //printf ("%d: Allocate %d bytes of storage for %d PME Pencils\n", CkMyPe(), 
  //	  nzlines * zlen * (int)sizeof(double), nactive);
  
  assert (zline_storage != NULL);
  double * zblock = zline_storage;
  
  for (int ib=0; ib<myGrid.xBlocks; ++ib) {
    for (int jb=0; jb<myGrid.yBlocks; ++jb) {
      int index = ib * myGrid.yBlocks + jb;
      if (pencilActive[index].isActive) {
	pencilActive[index].data = zblock;
	for (int i = pencilActive[index].xmin; i <= pencilActive[index].xmax; i++)
	  for (int j = pencilActive[index].ymin; j <= pencilActive[index].ymax; j++) {	    
	    q_arr[i*myGrid.dim2 + j] = zblock;
	    zblock += zlen;
	    assert ((char *) zblock <= (char *)(zline_storage + nzlines*zlen));
	  }
      }
    }
  }
  
  pencilVec.resize (nactive);
  nactive = 0;
  for (int ib=0; ib<myGrid.xBlocks; ++ib) 
    for (int jb=0; jb<myGrid.yBlocks; ++jb) 
      if (pencilActive[ib*myGrid.yBlocks + jb].isActive) 
	pencilVec[nactive ++] = pencilActive [ib * myGrid.yBlocks + jb];

  //We dont need the sparse array anymore
  delete [] pencilActive;
  
  /******************************* Initialize Many to Many ***********************************/
  OptPmeDummyMsg *m = new (PRIORITY_SIZE) OptPmeDummyMsg;
  m->to_pe = CkMyPe();
  CmiDirect_manytomany_initialize_recvbase (myMgr->handle, PHASE_UG, 
					    recv_ungrid_done, m, (char *)sp_zstorage, 
					    myGrid.xBlocks*myGrid.yBlocks, -1); 
  
  CkCallback cbi (CkCallback::ignore);
  CmiDirect_manytomany_initialize_sendbase (myMgr->handle, PHASE_GR, NULL, NULL, 
					    (char *)sp_zstorage, 
					    pencilVec.size(), patchList[0].patchID);

  for (int idx = 0; idx < pencilVec.size(); idx++) {
    int ib =  pencilVec[idx].ib;
    int jb =  pencilVec[idx].jb;
    double * data = pencilVec[idx].data;
    int offset = (data - zline_storage)*sizeof(float);
    int xlen   = pencilVec[idx].xmax - pencilVec[idx].xmin + 1;
    int ylen   = pencilVec[idx].ymax - pencilVec[idx].ymin + 1;
    int fcount = xlen * ylen * zlen;

    CkArrayIndex3D index (ib, jb, 0);
    CProxy_OptPmePencilMapZ zproxy (global_map_z);
    int pe = zproxy.ckLocalBranch()->procNum(0, index);
    CmiDirect_manytomany_initialize_send (myMgr->handle, PHASE_GR, idx, 
					   offset, fcount*sizeof(float), pe);
    
    int srcnode = ib * myGrid.yBlocks + jb;
    CmiDirect_manytomany_initialize_recv (myMgr->handle, PHASE_UG, 
					  srcnode, offset, fcount*sizeof(float), pe);
  }
  /********************************** End initialize many to many ****************************/

}

OptPmeCompute::~OptPmeCompute()
{
  delete [] zline_storage;
  delete [] sp_zstorage;
  delete [] q_arr;
}


void OptPmeCompute::doWork()
{
  DebugM(4,"Entering OptPmeCompute::doWork().\n");

#ifdef TRACE_COMPUTE_OBJECTS
  double traceObjStartTime = CmiWallTimer();
#endif

  if (!_initialized) initializeOptPmeCompute();

  ResizeArrayIter<PatchElem> ap(patchList);

  // Skip computations if nothing to do.
  if ( ! patchList[0].p->flags.doFullElectrostatics )
  {
    for (ap = ap.begin(); ap != ap.end(); ap++) {
      CompAtom *x = (*ap).positionBox->open();
      Results *r = (*ap).forceBox->open();
      (*ap).positionBox->close(&x);
      (*ap).forceBox->close(&r);
    }
    reduction->submit();
    return;
  }

  myMgr->_iter ++;  //this is a pme step

  // allocate storage
  numLocalAtoms = 0;
  for (ap = ap.begin(); ap != ap.end(); ap++) {
    numLocalAtoms += (*ap).p->getNumAtoms();
  }

  Lattice &lattice = patchList[0].p->flags.lattice;

  localData = new PmeParticle[numLocalAtoms];
  // get positions and charges
  PmeParticle * data_ptr = localData;

  int natoms = 0;
  //  if (myMgr->constant_pressure)
  //resetPatchCoordinates(lattice);  //Update patch coordinates with new lattice

  for (ap = ap.begin(); ap != ap.end(); ap++) {
#ifdef NETWORK_PROGRESS
    CmiNetworkProgress();
#endif
    
    CompAtom *x = (*ap).positionBox->open();
    if ( patchList[0].p->flags.doMolly ) {
      (*ap).positionBox->close(&x);
      x = (*ap).avgPositionBox->open();
    }
    int numAtoms = (*ap).p->getNumAtoms();        
    int order_1  = myGrid.order - 1;
    scale_n_copy_coordinates(x, localData, numAtoms, 
			     lattice, myGrid,
			     xstart + order_1, xlen - order_1,
			     ystart + order_1, ylen - order_1,
			     zstart + order_1, zlen - order_1,
			     strayChargeErrors);
    natoms += numAtoms;
    
    if ( patchList[0].p->flags.doMolly ) { (*ap).avgPositionBox->close(&x); }
    else { (*ap).positionBox->close(&x); }
  }

  numLocalAtoms = natoms;  //Exclude all atoms out of range

  // calculate self energy
  BigReal ewaldcof = ComputeNonbondedUtil::ewaldcof;
  evir = 0;
  BigReal selfEnergy = 0;
  data_ptr = localData;
  int i;
  for(i=0; i<numLocalAtoms; ++i)
  {
    selfEnergy += data_ptr->cg * data_ptr->cg;
    ++data_ptr;
  }
  selfEnergy *= -1. * ewaldcof / SQRT_PI;
  evir[0] += selfEnergy;

#if 0
  if (myMgr->_iter > many_to_many_start) {
    OptPmeSubComputeMsg *smsg = myMgr->subcompute_msgs[1]; //not self
    CProxy_OptPmeMgr pmeProxy (CkpvAccess(BOCclass_group).computePmeMgr);
    pmeProxy[CmiNodeFirst(CmiMyNode())+smsg->dest].doWorkOnPeer(smsg);
  }
  else
#endif
    doWorkOnPeer();

#ifdef TRACE_COMPUTE_OBJECTS
  traceUserBracketEvent(TRACE_COMPOBJ_IDOFFSET+this->cid, traceObjStartTime, CmiWallTimer());
#endif

}

void OptPmeCompute::doWorkOnPeer()
{
  Lattice &lattice = patchList[0].p->flags.lattice;
  double **q = q_arr;
  memset( (void*) zline_storage, 0, zlen * nzlines * sizeof(double) );

  myRealSpace = new OptPmeRealSpace(myGrid,numLocalAtoms);
  if (!strayChargeErrors)
    myRealSpace->fill_charges(q, localData, zstart, zlen);

  if (myMgr->constant_pressure && patchList[0].patchID == 0)
    myMgr->xPencil.recvLattice (lattice);
  
  if (myMgr->_iter <= many_to_many_start)
    sendPencils();
  else {
    pme_d2f (sp_zstorage, zline_storage, nzlines * zlen);
    CmiDirect_manytomany_start (myMgr->handle, PHASE_GR);
  }
}



void OptPmeCompute::sendPencils() {  

  //iout << iPE << " Sending charge grid for " << numLocalAtoms << " atoms to FFT with " << myMgr->numPencilsActive << " messages" <<".\n" << endi;

  int xBlocks = myGrid.xBlocks;
  int yBlocks = myGrid.yBlocks;
  int zBlocks = myGrid.zBlocks;

  int K1 = myGrid.K1;
  int K2 = myGrid.K2;
  int dim2 = myGrid.dim2;
  int dim3 = myGrid.dim3;
  int block1 = myGrid.block1;
  int block2 = myGrid.block2;

  //Lattice lattice = patchList[0].p->flags.lattice;

  int nactive = 0;  
  for (int idx = 0; idx < pencilVec.size(); idx++) {
    int xstart = pencilVec[idx].xmin;
    int ystart = pencilVec[idx].ymin;
    int xlen   = pencilVec[idx].xmax - pencilVec[idx].xmin + 1;
    int ylen   = pencilVec[idx].ymax - pencilVec[idx].ymin + 1;
    int ib     = pencilVec[idx].ib;
    int jb     = pencilVec[idx].jb;
    double *data = pencilVec[idx].data;
    
    int fcount = xlen * ylen * zlen;
    OptPmeGridMsg *msg = new (fcount, PRIORITY_SIZE) OptPmeGridMsg;
    msg->zstart = zstart;
    msg->zlen   = zlen;
    msg->xstart = xstart;
    msg->xlen   = xlen;
    msg->ystart = ystart;
    msg->ylen   = ylen;
    msg->sourceNode = CkMyPe();
    msg->patchID    = patchList[0].patchID;
    
    float *qmsg = msg->qgrid;	
#pragma disjoint (*data, *qmsg)
#pragma unroll(8)
    for ( int k=0; k< fcount; ++k ) 
      *(qmsg++) = data[k];    
    
    myMgr->zPencil(ib,jb,0).recvGrid(msg);
  }
}


void OptPmeCompute::copyPencils(OptPmeGridMsg *msg) {

  if (!_initialized) initializeOptPmeCompute();

  int ibegin = msg->xstart;
  int iend   = msg->xstart + msg->xlen;
  int jbegin = msg->ystart;
  int jend   = msg->ylen;
  int fcount = zlen * msg->xlen * msg->ylen;

  float *qmsg = msg->qgrid;
  double *data = q_arr[ibegin * myGrid.dim2 + jbegin];  

#pragma disjoint (*qmsg, *data)
#pragma unroll(8)
  for ( int k=0; k<fcount; ++k ) 
    data[k] = *(qmsg++);
}

void OptPmeCompute::ungridForces_init() {

    //printf ("%d: In OptPMECompute::ungridforces_init\n", CkMyPe());

    if (myMgr->_iter > many_to_many_start)
      pme_f2d (zline_storage, sp_zstorage, nzlines * zlen);

    localResults = new Vector[numLocalAtoms];
    //memset (localResults, 0, sizeof (Vector) * numLocalAtoms);
}

void OptPmeCompute::ungridForces_compute(int    istart,
					 int    iend) 
{
    Vector *gridResults;
    gridResults = localResults;

    if (iend == 0)
      iend = numLocalAtoms;
    
    SimParameters *simParams = Node::Object()->simParameters;
    Vector pairForce = 0.;
    Lattice &lattice = patchList[0].p->flags.lattice;
    if(!simParams->commOnly) {
#ifdef NETWORK_PROGRESS
      CmiNetworkProgress();
#endif      
      if (!strayChargeErrors) {
	myRealSpace->compute_forces(q_arr, localData, gridResults, 
				    zstart, zlen, istart, iend);
	scale_forces(gridResults + istart, iend - istart, lattice);
      }
    }
}

void OptPmeCompute::ungridForces_finalize() {
    SimParameters *simParams = Node::Object()->simParameters;
    delete myRealSpace;
    
    delete [] localData;
    //    delete [] localPartition;
    
    Vector *results_ptr = localResults;
    ResizeArrayIter<PatchElem> ap(patchList);
    
    // add in forces
    for (ap = ap.begin(); ap != ap.end(); ap++) {
      Results *r = (*ap).forceBox->open();
      Force *f = r->f[Results::slow];
      int numAtoms = (*ap).p->getNumAtoms();
      
      if ( ! strayChargeErrors && ! simParams->commOnly ) {
        for(int i=0; i<numAtoms; ++i) {
          f[i].x += results_ptr->x;
          f[i].y += results_ptr->y;
          f[i].z += results_ptr->z;
          ++results_ptr;
        }
      }
  
      (*ap).forceBox->close(&r);
    }

    delete [] localResults;
   
    double scale = 1.;

    reduction->item(REDUCTION_ELECT_ENERGY_SLOW) += evir[0] * scale;
    reduction->item(REDUCTION_STRAY_CHARGE_ERRORS) += strayChargeErrors;
    strayChargeErrors = 0;
    reduction->submit();
}



#include "fftlib.C"
#include "OptPmeMgr.def.h"

void pme_f2d (double *dst, float *src, int N) {
#pragma disjoint (*src, *dst)  
#pragma unroll (8)
  for (int i = 0; i < N; i++) 
    dst[i] = src[i];
}

void pme_d2f (float *dst, double *src, int N) {
#pragma disjoint (*src, *dst)  
#pragma unroll (8)
  for (int i = 0; i < N; i++) 
    dst[i] = src[i];
}
