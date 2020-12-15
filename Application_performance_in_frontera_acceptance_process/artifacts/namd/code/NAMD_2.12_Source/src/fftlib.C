
#include "fftlib.h"

#include <assert.h>

void call_ck_cb (void *arg) {
  CkCallbackWrapper *cw = (CkCallbackWrapper *)arg;
  cw->cb.send (cw->msg);
}

void call_ck_cb_recv_grid (void *arg) {
  CkCallbackWrapper *cw = (CkCallbackWrapper *)arg;
  OptPmeDummyMsg *msg = (OptPmeDummyMsg *)cw->msg;
  OptPmeZPencil *zp = (OptPmeZPencil *)cw->array;
  zp->many_to_manyRecvGrid(msg);
}

void call_ck_cb_recv_trans_y (void *arg) {
  CkCallbackWrapper *cw = (CkCallbackWrapper *)arg;
  OptPmeDummyMsg *msg = (OptPmeDummyMsg *)cw->msg;
  OptPmeYPencil *yp = (OptPmeYPencil *)cw->array;
  yp->many_to_manyRecvTrans(msg);
}

void call_ck_cb_recv_trans_x (void *arg) {
  CkCallbackWrapper *cw = (CkCallbackWrapper *)arg;
  OptPmeDummyMsg *msg = (OptPmeDummyMsg *)cw->msg;
  OptPmeXPencil *xp = (OptPmeXPencil *)cw->array;
  xp->many_to_manyRecvTrans(msg);
}

void call_ck_cb_recv_untrans_y (void *arg) {
  CkCallbackWrapper *cw = (CkCallbackWrapper *)arg;
  OptPmeDummyMsg *msg = (OptPmeDummyMsg *)cw->msg;
  OptPmeYPencil *yp = (OptPmeYPencil *)cw->array;
  yp->many_to_manyRecvUntrans(msg);
}

void call_ck_cb_recv_untrans_z (void *arg) {
  CkCallbackWrapper *cw = (CkCallbackWrapper *)arg;
  OptPmeDummyMsg *msg = (OptPmeDummyMsg *)cw->msg;
  OptPmeZPencil *zp = (OptPmeZPencil *)cw->array;
  zp->many_to_manyRecvUntrans(msg);
}

void OptPmeZPencil::fft_init() {
  
  //printf ("Initialize zpencil [%d,%d], on pd %d\n", thisIndex.x, thisIndex.y, CkMyPe());

  CProxy_Node nd(CkpvAccess(BOCclass_group).node);
  Node *node = nd.ckLocalBranch();
  SimParameters *simParams = node->simParameters;

  int K1 = initdata.grid.K1;
  int K2 = initdata.grid.K2;
  int K3 = initdata.grid.K3;
  int dim3 = initdata.grid.dim3;
  int block1 = initdata.grid.block1;
  int block2 = initdata.grid.block2;

  nx = block1;
  if ( (thisIndex.x + 1) * block1 > K1 ) nx = K1 - thisIndex.x * block1;
  ny = block2;
  if ( (thisIndex.y + 1) * block2 > K2 ) ny = K2 - thisIndex.y * block2;

  data = new float[nx*ny*dim3];
  many_to_many_data = new float[nx*ny*dim3];
  many_to_many_nb = new int [initdata.zBlocks];
  work = new float[dim3];

  memset(data, 0, sizeof(float) * nx*ny*dim3);
  memset(many_to_many_data, 0, sizeof(float) * nx*ny*dim3);

  order_init(initdata.zBlocks);

#ifdef NAMD_FFTW
  CmiLock(fftw_plan_lock);
#ifdef NAMD_FFTW_3
  /* need array of sizes for the how many */
  int numLines=nx*ny;
  int planLineSizes[1];
  planLineSizes[0]=K3;
  CkAbort("what are we doing in here?");
  forward_plan = fftwf_plan_many_dft_r2c(1, planLineSizes, numLines,
				     (float *) data, NULL, 1, 
					 initdata.grid.dim3,
					 (fftwf_complex *) data, NULL, 1, 0,
				   ( simParams->FFTWEstimate ? FFTW_ESTIMATE 
				     : FFTW_MEASURE ));
  backward_plan = fftwf_plan_many_dft_c2r(1, planLineSizes, numLines,
				     (fftwf_complex *) data, NULL, 1, 
					 initdata.grid.dim3/2,
				     (float *) data, NULL, 1, 0,
				   ( simParams->FFTWEstimate ? FFTW_ESTIMATE 
				     : FFTW_MEASURE ));
#else

  forward_plan = rfftwnd_create_plan_specific(1, &K3, FFTW_REAL_TO_COMPLEX,
	( simParams->FFTWEstimate ? FFTW_ESTIMATE : FFTW_MEASURE )
	| FFTW_IN_PLACE | FFTW_USE_WISDOM, data, 1, work, 1);
  backward_plan = rfftwnd_create_plan_specific(1, &K3, FFTW_COMPLEX_TO_REAL,
	( simParams->FFTWEstimate ? FFTW_ESTIMATE : FFTW_MEASURE )
	| FFTW_IN_PLACE | FFTW_USE_WISDOM, data, 1, work, 1);
#endif
  CmiUnlock(fftw_plan_lock);
#else
  NAMD_die("Sorry, FFTW must be compiled in to use PME.");
#endif

  handle = CmiDirect_manytomany_allocate_handle();
}

void OptPmeYPencil::fft_init() {
  CProxy_Node nd(CkpvAccess(BOCclass_group).node);
  Node *node = nd.ckLocalBranch();
  SimParameters *simParams = node->simParameters;

  //  printf ("Initialize ypencil [%d,%d], on pd %d\n", thisIndex.x, thisIndex.y, CkMyPe());

  int K1 = initdata.grid.K1;
  int K2 = initdata.grid.K2;
  int dim2 = initdata.grid.dim2;
  int dim3 = initdata.grid.dim3;
  int block1 = initdata.grid.block1;
  int block3 = initdata.grid.block3;

  nx = block1;
  if ( (thisIndex.x + 1) * block1 > K1 ) nx = K1 - thisIndex.x * block1;
  nz = block3;
  if ( (thisIndex.z+1)*block3 > dim3/2 ) nz = dim3/2 - thisIndex.z*block3;

  data = new float[nx*dim2*nz*2];
  many_to_many_data = new float[nx*dim2*nz*2];
  many_to_many_nb = new int [initdata.yBlocks];
  work = new float[2*K2];

  memset(many_to_many_data, 0, sizeof(float) * nx*dim2*nz*2);

  order_init(initdata.yBlocks);

#ifdef NAMD_FFTW
  CmiLock(fftw_plan_lock);
#ifdef NAMD_FFTW_3
  /* need array of sizes for the dimensions */
  int numLines=nz;
  int planLineSizes[2];
  planLineSizes[0]=initdata.grid.K2;
  planLineSizes[1]=nz;
  forward_plan = fftwf_plan_many_dft(2, planLineSizes, numLines, 
				     (fftwf_complex *) data, NULL, nz, 1,
				     (fftwf_complex *) data, NULL, 1, 0,
				     FFTW_FORWARD, 
				     ( simParams->FFTWEstimate 
				       ? FFTW_ESTIMATE : FFTW_MEASURE ));
  backward_plan = fftwf_plan_many_dft(2, planLineSizes, numLines, 
				     (fftwf_complex *) data, NULL, nz, 1,
				     (fftwf_complex *) data, NULL, 1, 0,
				     FFTW_FORWARD, 
				     ( simParams->FFTWEstimate 
				       ? FFTW_ESTIMATE : FFTW_MEASURE ));
#else

  forward_plan = fftw_create_plan_specific(K2, FFTW_FORWARD,
	( simParams->FFTWEstimate ? FFTW_ESTIMATE : FFTW_MEASURE )
	| FFTW_IN_PLACE | FFTW_USE_WISDOM, (fftw_complex *) data,
	nz, (fftw_complex *) work, 1);
  backward_plan = fftw_create_plan_specific(K2, FFTW_BACKWARD,
	( simParams->FFTWEstimate ? FFTW_ESTIMATE : FFTW_MEASURE )
	| FFTW_IN_PLACE | FFTW_USE_WISDOM, (fftw_complex *) data,
	nz, (fftw_complex *) work, 1);
#endif
  CmiUnlock(fftw_plan_lock);
#else
  NAMD_die("Sorry, FFTW must be compiled in to use PME.");
#endif

  handle = CmiDirect_manytomany_allocate_handle();
  initialize_manytomany();
}


void OptPmeXPencil::fft_init() {
  CProxy_Node nd(CkpvAccess(BOCclass_group).node);
  Node *node = nd.ckLocalBranch();
  SimParameters *simParams = node->simParameters;

  //  printf ("Initialize xpencil [%d,%d], on pd %d\n", thisIndex.x, thisIndex.y, CkMyPe());

  lattice = simParams->lattice;

  int K1 = initdata.grid.K1;
  int K2 = initdata.grid.K2;
  int dim3 = initdata.grid.dim3;
  int block2 = initdata.grid.block2;
  int block3 = initdata.grid.block3;

  ny = block2;
  if ( (thisIndex.y + 1) * block2 > K2 ) ny = K2 - thisIndex.y * block2;
  nz = block3;
  if ( (thisIndex.z+1)*block3 > dim3/2 ) nz = dim3/2 - thisIndex.z*block3;

  data = new float[K1*block2*block3*2];
  many_to_many_data = new float[K1*block2*block3*2];
  many_to_many_nb = new int [initdata.xBlocks];
  work = new float[2*K1];

  memset(many_to_many_data, 0, sizeof(float) * K1*block2*block3*2);

  order_init(initdata.xBlocks);

#ifdef NAMD_FFTW
  CmiLock(fftw_plan_lock);
#ifdef NAMD_FFTW_3
  /* need array of sizes for the how many */
  int numLines=ny*nz;
  int planLineSizes[1];
  planLineSizes[0]=K1;
  forward_plan = fftwf_plan_many_dft(1, planLineSizes, numLines,
				     (fftwf_complex *) data, NULL, K1, 1,
				     (fftwf_complex *) data, NULL, 1, 0,
				   FFTW_FORWARD,
				   ( simParams->FFTWEstimate ? FFTW_ESTIMATE 
				     : FFTW_MEASURE ));
  backward_plan = fftwf_plan_many_dft(1, planLineSizes, numLines,
				     (fftwf_complex *) data, NULL, K1, 1,
				     (fftwf_complex *) data, NULL, 1, 0,
				   FFTW_BACKWARD,
				   ( simParams->FFTWEstimate ? FFTW_ESTIMATE 
				     : FFTW_MEASURE ));
#else

  forward_plan = fftw_create_plan_specific(K1, FFTW_FORWARD,
	( simParams->FFTWEstimate ? FFTW_ESTIMATE : FFTW_MEASURE )
	| FFTW_IN_PLACE | FFTW_USE_WISDOM, (fftw_complex *) data,
	ny*nz, (fftw_complex *) work, 1);
  backward_plan = fftw_create_plan_specific(K1, FFTW_BACKWARD,
	( simParams->FFTWEstimate ? FFTW_ESTIMATE : FFTW_MEASURE )
	| FFTW_IN_PLACE | FFTW_USE_WISDOM, (fftw_complex *) data,
	ny*nz, (fftw_complex *) work, 1);

  CmiUnlock(fftw_plan_lock);
#endif
#else
  NAMD_die("Sorry, FFTW must be compiled in to use PME.");
#endif

  myKSpace = new PmeKSpace(initdata.grid,
		thisIndex.y*block2, thisIndex.y*block2 + ny,
		thisIndex.z*block3, thisIndex.z*block3 + nz);

  handle = CmiDirect_manytomany_allocate_handle();
  initialize_manytomany();

  constant_pressure = initdata.constant_pressure;

  reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);
}

//#define FFTCHECK   // run a grid of integers through the fft
// #define ZEROCHECK  // check for suspicious zeros in fft

void OptPmeZPencil::recv_grid(const OptPmeGridMsg *msg) {

  int dim3 = initdata.grid.dim3;
  if ( imsg == 0 ) {
    memset(data, 0, sizeof(float) * nx*ny*dim3);
  }

  int xstart = msg->xstart - thisIndex.x * initdata.grid.block1;
  int ystart = msg->ystart - thisIndex.y * initdata.grid.block2;
  assert (xstart >= 0);
  assert (ystart >= 0);
  int xlen   = msg->xlen;
  int ylen   = msg->ylen;    
  assert (xstart + xlen <= nx);
  assert (ystart + ylen <= ny);
  
  float *qmsg = msg->qgrid;
  float *d = data;
  int zlen = msg->zlen;
  int zstart = msg->zstart;
  int zmax = zstart + zlen - 1;
  int k = 0;

  int K3 = initdata.grid.K3; 
  int K3_1 = initdata.grid.K3 - 1; 
  while (zstart < 0) {
    zstart += K3;
    zmax += K3;
  }

  for ( int i=xstart; i<xstart+xlen; ++i ) {
    for ( int j=ystart; j<ystart+ylen; ++j) {
      float *d = data + (i * ny + j) * dim3; 
      for ( k=zstart; k<=zmax; ++k ) {
	int kz = k;
	//if (kz >= K3) kz -= K3;
	kz = kz - ((unsigned)(K3_1 - kz)>>31)*K3;
	//assert (kz >= 0);
	//assert (kz <  K3);
	d[kz] += *(qmsg++);
      }
    }
  }

  if (_iter == MANY_TO_MANY_SETUP) {
    if (imsg == 0)
      m2m_recv_grid = new PatchGridElem [grid_msgs.size()];    
    
    m2m_recv_grid[imsg].xstart = xstart;
    m2m_recv_grid[imsg].xlen   = xlen;
    m2m_recv_grid[imsg].ystart = ystart;
    m2m_recv_grid[imsg].ylen   = ylen;
    m2m_recv_grid[imsg].zstart = zstart;
    m2m_recv_grid[imsg].zlen   = zlen;
    m2m_recv_grid[imsg].patchID = msg->patchID;
    
    if ( imsg == grid_msgs.size() - 1) 
        initialize_manytomany ();
  }
}


void OptPmeZPencil::forward_fft() {
#ifdef FFTCHECK
  int dim3 = initdata.grid.dim3;
  int K3 = initdata.grid.K3;
  float std_base = 100. * (thisIndex.x+1.) + 10. * (thisIndex.y+1.);
  float *d = data;
  for ( int i=0; i<nx; ++i ) {
   for ( int j=0; j<ny; ++j, d += dim3 ) {
    for ( int k=0; k<dim3; ++k ) {
      d[k] = 10. * (10. * (10. * std_base + i) + j) + k;
    }
   }
  }
#endif
#ifdef NAMD_FFTW
#ifdef NAMD_FFTW_3
  fftwf_execute(forward_plan);
#else
  rfftwnd_real_to_complex(forward_plan, nx*ny,
			  data, 1, initdata.grid.dim3, (fftw_complex *) work, 1, 0);
#endif
#endif
}

void OptPmeZPencil::send_trans() {
  int zBlocks = initdata.zBlocks;
  int block3 = initdata.grid.block3;
  int dim3 = initdata.grid.dim3;

  //int offset = 0;
  for ( int isend=0; isend<zBlocks; ++isend ) {
    int kb = send_order[isend];
    int nz = block3;
    if ( (kb+1)*block3 > dim3/2 ) nz = dim3/2 - kb*block3;
    //assert (nz > 0);
    OptPmeFFTMsg *msg = new (nx*ny*nz*2,PRIORITY_SIZE) OptPmeFFTMsg;
    msg->sourceNode = thisIndex.y;
    msg->nx = ny;
    float *md = msg->qgrid;
    const float *d = data;
    for ( int i=0; i<nx; ++i ) {
      for ( int j=0; j<ny; ++j, d += dim3 ) {
	for ( int k=kb*block3; k<(kb*block3+nz); ++k ) {
	  *(md++) = d[2*k];
	  *(md++) = d[2*k+1];
	}
      }
    }          
    //    printf ("%d, %d: Zpencil Sending trans to %d, %d\n", thisIndex.x, thisIndex.y, thisIndex.x, kb);
    initdata.yPencil(thisIndex.x,0,kb).recvTrans(msg);
  }
}


void OptPmeYPencil::recv_trans(const OptPmeFFTMsg *msg) {

  int block2 = initdata.grid.block2;
  int K2 = initdata.grid.K2;
  int jb = msg->sourceNode;
  int ny = msg->nx;
  const float *md = msg->qgrid;
  float *d = data;
  for ( int i=0; i<nx; ++i, d += K2*nz*2 ) {
    for ( int j=jb*block2; j<(jb*block2+ny); ++j ) {
      for ( int k=0; k<nz; ++k ) {
	d[2*(j*nz+k)] = *(md++);
	d[2*(j*nz+k)+1] = *(md++);
      }
    }
  }
} 

void OptPmeYPencil::forward_fft() {
#ifdef NAMD_FFTW
#ifdef NAMD_FFTW_3
  fftwf_execute(forward_plan);
#else
  for ( int i=0; i<nx; ++i ) {
    fftw(forward_plan, nz,
	 ((fftw_complex *) data) + i * nz * initdata.grid.K2,
	 nz, 1, (fftw_complex *) work, 1, 0);
  }
#endif
#endif
}

void OptPmeYPencil::send_trans() {
  int yBlocks = initdata.yBlocks;
  int block2 = initdata.grid.block2;
  int K2 = initdata.grid.K2;
  
  for ( int isend=0; isend<yBlocks; ++isend ) {
    int jb = send_order[isend];
    int ny = block2;
    if ( (jb+1)*block2 > K2 ) ny = K2 - jb*block2;
    OptPmeFFTMsg *msg = new (nx*ny*nz*2,PRIORITY_SIZE) OptPmeFFTMsg;
    msg->sourceNode = thisIndex.x;
    msg->nx = nx;
    float *md = msg->qgrid;
    const float *d = data;
    for ( int i=0; i<nx; ++i, d += K2*nz*2 ) {
      for ( int j=jb*block2; j<(jb*block2+ny); ++j ) {
	for ( int k=0; k<nz; ++k ) {
	  *(md++) = d[2*(j*nz+k)];
	  *(md++) = d[2*(j*nz+k)+1];
	}
      }
    }
    if ( md != msg->qgrid + nx*ny*nz*2 ) CkPrintf("error in YX at %d %d %d\n",
						  thisIndex.x, jb, thisIndex.z);
    
    //printf ("%d, %d: Ypencil Sending trans to %d, %d\n", thisIndex.z, thisIndex.x, jb, thisIndex.z);
    initdata.xPencil(0,jb,thisIndex.z).recvTrans(msg);
  }
}



void OptPmeXPencil::recv_trans(const OptPmeFFTMsg *msg) {

  int block1 = initdata.grid.block1;
  int K1 = initdata.grid.K1;
  int ib = msg->sourceNode;
  int nx = msg->nx;
  const float *md = msg->qgrid;
  for ( int i=ib*block1; i<(ib*block1+nx); ++i ) {
    float *d = data + i*ny*nz*2;
    for ( int j=0; j<ny; ++j, d += nz*2 ) {
      for ( int k=0; k<nz; ++k ) {
	d[2*k] = *(md++);
	d[2*k+1] = *(md++);
      }
    }
  }
}



void OptPmeXPencil::forward_fft() {
#ifdef NAMD_FFTW
#ifdef NAMD_FFTW_3
  fftwf_execute(forward_plan);
#else

  fftw(forward_plan, ny*nz,
       ((fftw_complex *) data), ny*nz, 1, (fftw_complex *) work, 1, 0);
#endif
#endif
}

void OptPmeXPencil::pme_kspace() {
  evir = 0.;  //set evir to 0

#ifdef FFTCHECK
  return;
#endif

  BigReal ewaldcof = ComputeNonbondedUtil::ewaldcof;
  evir[0] = myKSpace->compute_energy(data,
				     lattice, ewaldcof, &(evir[1]), 0);

  //contribute (7*sizeof(double), evir.begin(), CkReduction::sum_double, initdata.cb_energy);
}

void OptPmeXPencil::submit_evir() {
  double * cdata = (double *) evir.begin();
  reduction->item(REDUCTION_ELECT_ENERGY_SLOW) += cdata[0];
  reduction->item(REDUCTION_VIRIAL_SLOW_XX) += cdata[1];
  reduction->item(REDUCTION_VIRIAL_SLOW_XY) += cdata[2];
  reduction->item(REDUCTION_VIRIAL_SLOW_XZ) += cdata[3];
  reduction->item(REDUCTION_VIRIAL_SLOW_YX) += cdata[2];
  reduction->item(REDUCTION_VIRIAL_SLOW_YY) += cdata[4];
  reduction->item(REDUCTION_VIRIAL_SLOW_YZ) += cdata[5];
  reduction->item(REDUCTION_VIRIAL_SLOW_ZX) += cdata[3];
  reduction->item(REDUCTION_VIRIAL_SLOW_ZY) += cdata[5];
  reduction->item(REDUCTION_VIRIAL_SLOW_ZZ) += cdata[6];   

  SimParameters *simParams = Node::Object()->simParameters;
  int fef = simParams->fullElectFrequency;
  for (int i = 0; i < fef; i++) {
    reduction->submit();
  }
}

void OptPmeXPencil::backward_fft() {
#ifdef NAMD_FFTW
#ifdef NAMD_FFTW_3
  fftwf_execute(backward_plan);
#else

  fftw(backward_plan, ny*nz,
       ((fftw_complex *) data), ny*nz, 1, (fftw_complex *) work, 1, 0);
#endif
#endif
}

void OptPmeXPencil::send_untrans() {
  int xBlocks = initdata.xBlocks;
  int block1 = initdata.grid.block1;
  int K1 = initdata.grid.K1;

  for ( int isend=0; isend<xBlocks; ++isend ) {
    int ib = send_order[isend];
    int nx = block1;
    if ( (ib+1)*block1 > K1 ) nx = K1 - ib*block1;
    OptPmeFFTMsg *msg = new (nx*ny*nz*2,PRIORITY_SIZE) OptPmeFFTMsg;    
    msg->sourceNode = thisIndex.y;
    msg->nx = ny;
    float *md = msg->qgrid;
    for ( int i=ib*block1; i<(ib*block1+nx); ++i ) {
      float *d = data + i*ny*nz*2;
      for ( int j=0; j<ny; ++j, d += nz*2 ) {
	for ( int k=0; k<nz; ++k ) {
	  *(md++) = d[2*k];
	  *(md++) = d[2*k+1];
	}
      }
    }
    
    //printf ("%d, %d: Xpencil Sending untrans to %d, %d\n", thisIndex.y, thisIndex.z, ib, thisIndex.z);
    initdata.yPencil(ib,0,thisIndex.z).recvUntrans(msg);
  }
}


void OptPmeYPencil::recv_untrans(const OptPmeFFTMsg *msg) {

  int block2 = initdata.grid.block2;
  int K2 = initdata.grid.K2;
  int jb = msg->sourceNode;
  int ny = msg->nx;
  const float *md = msg->qgrid;
  float *d = data;
  for ( int i=0; i<nx; ++i, d += K2*nz*2 ) {
    for ( int j=jb*block2; j<(jb*block2+ny); ++j ) {
      for ( int k=0; k<nz; ++k ) {
	d[2*(j*nz+k)] = *(md++);
	d[2*(j*nz+k)+1] = *(md++);
      }
    }
  }
} 



void OptPmeYPencil::backward_fft() {
#ifdef NAMD_FFTW
#ifdef NAMD_FFTW_3
  fftwf_execute(backward_plan);
#else
  for ( int i=0; i<nx; ++i ) {
#if CMK_BLUEGENEL
    CmiNetworkProgress();
#endif

    fftw(backward_plan, nz,
	 ((fftw_complex *) data) + i * nz * initdata.grid.K2,
	 nz, 1, (fftw_complex *) work, 1, 0);
  }
#endif
#endif
}

void OptPmeYPencil::send_untrans() {
  //printf ("%d, %d: In ypencil send_untrans called once \n", thisIndex.x, thisIndex.z);

  int yBlocks = initdata.yBlocks;
  int block2 = initdata.grid.block2;
  int K2 = initdata.grid.K2;

  for ( int isend=0; isend<yBlocks; ++isend ) {
    int jb = send_order[isend];
    //if ( ! needs_reply[jb] ) continue;
    int ny = block2;
    if ( (jb+1)*block2 > K2 ) ny = K2 - jb*block2;
    OptPmeFFTMsg *msg = new (nx*ny*nz*2,PRIORITY_SIZE) OptPmeFFTMsg;
    msg->sourceNode = thisIndex.z;
    msg->nx = nz;
    float *md = msg->qgrid;
    const float *d = data;
    for ( int i=0; i<nx; ++i, d += K2*nz*2 ) {
     for ( int j=jb*block2; j<(jb*block2+ny); ++j ) {
      for ( int k=0; k<nz; ++k ) {
        *(md++) = d[2*(j*nz+k)];
        *(md++) = d[2*(j*nz+k)+1];
      }
     }
    }

    //printf ("%d, %d: Sending untrans to %d, %d\n", thisIndex.z, thisIndex.x, thisIndex.x, jb);
    initdata.zPencil(thisIndex.x,jb,0).recvUntrans(msg);
  }
}

void OptPmeZPencil::recv_untrans(const OptPmeFFTMsg *msg) {

  //printf ("%d, %d In recv untrans\n", thisIndex.x, thisIndex.y);

  int block3 = initdata.grid.block3;
  int dim3 = initdata.grid.dim3;
  int kb = msg->sourceNode;
  int nz = msg->nx;
  const float *md = msg->qgrid;
  float *d = data;
  for ( int i=0; i<nx; ++i ) {
    for ( int j=0; j<ny; ++j, d += dim3 ) {
      for ( int k=kb*block3; k<(kb*block3+nz); ++k ) {
	d[2*k] = *(md++);
	d[2*k+1] = *(md++);
      }
    }
  }
}


void OptPmeZPencil::backward_fft() {
#ifdef NAMD_FFTW
#ifdef NAMD_FFTW_3
  fftwf_execute(backward_plan);
#else
  rfftwnd_complex_to_real(backward_plan, nx*ny,
			  (fftw_complex *) data, 1, initdata.grid.dim3/2, work, 1, 0);
#endif
#endif
  
#if CMK_BLUEGENEL
  CmiNetworkProgress();
#endif

#ifdef FFTCHECK
  int dim3 = initdata.grid.dim3;
  int K1 = initdata.grid.K1;
  int K2 = initdata.grid.K2;
  int K3 = initdata.grid.K3;
  float scale = 1. / (1. * K1 * K2 * K3);
  float maxerr = 0.;
  float maxstd = 0.;
  int mi, mj, mk;  mi = mj = mk = -1;
  float std_base = 100. * (thisIndex.x+1.) + 10. * (thisIndex.y+1.);
  const float *d = data;
  for ( int i=0; i<nx; ++i ) {
   for ( int j=0; j<ny; ++j, d += dim3 ) {
    for ( int k=0; k<K3; ++k ) {
      float std = 10. * (10. * (10. * std_base + i) + j) + k;
      float err = scale * d[k] - std;
      if ( fabsf(err) > fabsf(maxerr) ) {
        maxerr = err;
        maxstd = std;
        mi = i;  mj = j;  mk = k;
      }
    }
   }
  }
  CkPrintf("pencil %d %d max error %f at %d %d %d (should be %f)\n",
	   thisIndex.x, thisIndex.y, maxerr, mi, mj, mk, maxstd);
#endif
}

void OptPmeZPencil::send_ungrid(OptPmeGridMsg *msg) {
  int pe = msg->sourceNode;
  msg->sourceNode = thisIndex.x * initdata.yBlocks + thisIndex.y;

  int dim3 = initdata.grid.dim3;
  int xstart = msg->xstart - thisIndex.x * initdata.grid.block1;
  int ystart = msg->ystart - thisIndex.y * initdata.grid.block2;
  assert (xstart >= 0);
  assert (ystart >= 0);
  int xlen   = msg->xlen;
  int ylen   = msg->ylen;
  assert (xstart + xlen <= nx);
  assert (ystart + ylen <= ny);

  float *qmsg = msg->qgrid;
  float *d = data;
  int zstart = msg->zstart;
  int zlen = msg->zlen;
  int zmax = zstart + zlen - 1;
  
  int K3_1 = initdata.grid.K3 - 1; 
  int K3 = initdata.grid.K3;
  if (zstart < 0) {
    zstart += K3;
    zmax += K3;
  }
  
  int k = 0;
  for ( int i=xstart; i<xstart+xlen; ++i ) {
    for ( int j=ystart; j<ystart+ylen; ++j) {
      float *d = data + (i * ny + j) * dim3;
      for ( k=zstart; k<=zmax; ++k ) {
	int kz = k;
	//if (kz >= K3) kz -= K3;
	kz = kz - ((unsigned)(K3_1 - kz)>>31)*K3;
	*(qmsg++) = d[kz];
      }
    }
  }
  
  initdata.pmeProxy[pe].recvUngrid(msg);
}


//////////////////////////////////////////////////////////////////////
//////////////////////  Many to Many Implementation //////////////////
//////////////////////////////////////////////////////////////////////

void OptPmeZPencil::many_to_many_recv_grid () {
  int dim3 = initdata.grid.dim3;  
  memset(data, 0, sizeof(float) * nx*ny*dim3);
  int K3 = initdata.grid.K3;
  int K3_1 = initdata.grid.K3 - 1;
  int k = 0;
  for (int idx = 0; idx < grid_msgs.size(); idx++) {
    int xstart = m2m_recv_grid[idx].xstart; 
    int xlen   = m2m_recv_grid[idx].xlen; 
    int ystart = m2m_recv_grid[idx].ystart; 
    int ylen   = m2m_recv_grid[idx].ylen; 
    int zstart = m2m_recv_grid[idx].zstart; 
    int zlen   = m2m_recv_grid[idx].zlen;   
    
    float *qmsg = m2m_recv_grid[idx].data;
    for ( int i=xstart; i<xstart+xlen; ++i ) {
      for ( int j=ystart; j<ystart+ylen; ++j) {
	float *d = data + (i * ny + j) * dim3; 	

#pragma disjoint (*qmsg, *d)
#pragma unroll (4)
	for ( k=zstart; k<zstart+zlen; ++k ) {
	  int kz = k;
	  //if (kz >= K3) kz -= K3;
	  kz = kz - ((unsigned)(K3_1 - kz)>>31)*K3;
	  //assert (kz >= 0);
	  //assert (kz <  K3);
	  d[kz] += *(qmsg++);
	}
      }
    } 
  }
}

void OptPmeZPencil::many_to_many_send_trans() {
  int zBlocks = initdata.zBlocks;
  int block3 = initdata.grid.block3;
  int dim3 = initdata.grid.dim3;

  float *md = many_to_many_data;
  if (single_pencil) {
    const float *d = data;
    for ( int kb=0; kb<zBlocks; ++kb ) {
      *(md++) = d[2*kb];
      *(md++) = d[2*kb+1];
    }
  }
  else {
    int nz = block3;
    for ( int kb=0; kb<zBlocks; ++kb ) {
      nz = block3;
      if ( (kb+1)*block3 > dim3/2 ) nz = dim3/2 - kb*block3;
      
      const float *d = data;
      for ( int i=0; i<nx; ++i ) {
	for ( int j=0; j<ny; ++j, d += dim3 ) {
	  for ( int k=kb*block3; k<(kb*block3+nz); ++k ) {
	    *(md++) = d[2*k];
	    *(md++) = d[2*k+1];
	  }
	}
      }
    }
  }

  CmiDirect_manytomany_start (handle, PHASE_YF);
}

void OptPmeYPencil::many_to_many_recv_trans () {  
  int block2 = initdata.grid.block2;
  int K2 = initdata.grid.K2;
  
  const float *md = many_to_many_data;
  if (single_pencil) {
    float *d = data;
    for (int jb = 0; jb < initdata.yBlocks; jb++ ) {
      d[2*jb]    = *(md++);
      d[2*jb +1] = *(md++);
    }
  }
  else {
    for (int jb = 0; jb < initdata.yBlocks; jb++ ) {
      int ny = many_to_many_nb[jb];  
      float *d = data;
      for ( int i=0; i<nx; ++i, d += K2*nz*2 ) {
	for ( int j=jb*block2; j<(jb*block2+ny); ++j ) {
	  for ( int k=0; k<nz; ++k ) {
	    d[2*(j*nz+k)] = *(md++);
	    d[2*(j*nz+k)+1] = *(md++);
	  }
	}
      }
    }
  }
}

void OptPmeYPencil::many_to_many_send(int phase) {
  int yBlocks = initdata.yBlocks;
  int block2 = initdata.grid.block2;
  int K2 = initdata.grid.K2;
  
  float *md = many_to_many_data;
  int ny = block2;
  if (single_pencil) {
    const float *d = data;
    for ( int jb=0; jb<yBlocks; ++jb ) {
      *(md++) = d[2*jb];
      *(md++) = d[2*jb+1];
    }
  }
  else {
    for ( int jb=0; jb<yBlocks; ++jb ) {
      ny = block2;
      if ( (jb+1)*block2 > K2 ) ny = K2 - jb*block2;
      
      const float *d = data;
      for ( int i=0; i<nx; ++i, d += K2*nz*2 ) {
	for ( int j=jb*block2; j<(jb*block2+ny); ++j ) {
	  for ( int k=0; k<nz; ++k ) {
	    *(md++) = d[2*(j*nz+k)];
	    *(md++) = d[2*(j*nz+k)+1];
	  }
	}
      }
    }
  }

  CmiDirect_manytomany_start (handle, phase);  
}

void OptPmeXPencil::many_to_many_recv_trans () {
  
  int block1 = initdata.grid.block1;
  int K1 = initdata.grid.K1;
  
  const float *md = many_to_many_data;
  if (single_pencil) {
    for (int ib =0; ib < initdata.xBlocks; ib++ ) {
      data[2*ib]   = *(md++);
      data[2*ib+1] = *(md++);
    }
  }
  else {
    for (int ib =0; ib < initdata.xBlocks; ib++ ) {
      int nx = many_to_many_nb[ib];    
      for ( int i=ib*block1; i<(ib*block1+nx); ++i ) {
	float *d = data + i*ny*nz*2;
	for ( int j=0; j<ny; ++j, d += nz*2 ) {
	  for ( int k=0; k<nz; ++k ) {
	    d[2*k] = *(md++);
	    d[2*k+1] = *(md++);
	  }
	}
      }
    }
  }
}

void OptPmeXPencil::many_to_many_send_untrans() {
  int xBlocks = initdata.xBlocks;
  int block1 = initdata.grid.block1;
  int K1 = initdata.grid.K1;
  
  int     nx = block1;
  float * md = many_to_many_data;
  if (single_pencil) {
    float *d = data;
    for ( int ib=0; ib<xBlocks; ++ib ) {
      *(md++) = d[2*ib];
      *(md++) = d[2*ib+1];
    }
  }
  else {
    for ( int ib=0; ib<xBlocks; ++ib ) {
      nx = block1;
      if ( (ib+1)*block1 > K1 ) nx = K1 - ib*block1;
      
      for ( int i=ib*block1; i<(ib*block1+nx); ++i ) {
	float *d = data + i*ny*nz*2;
	for ( int j=0; j<ny; ++j, d += nz*2 ) {
	  for ( int k=0; k<nz; ++k ) {
	    *(md++) = d[2*k];
	    *(md++) = d[2*k+1];
	  }
	}
      }
    }
  }
  CmiDirect_manytomany_start (handle, PHASE_YB);
}

void OptPmeYPencil::many_to_many_recv_untrans () {  
  int block2 = initdata.grid.block2;
  int K2 = initdata.grid.K2;
  
  const float *md = many_to_many_data;
  if (single_pencil) {
    float *d = data;
    for (int jb = 0; jb < initdata.yBlocks; jb++ ) {        
      d[2*jb]   = *(md++);
      d[2*jb+1] = *(md++);
    }
  }
  else {
    for (int jb = 0; jb < initdata.yBlocks; jb++ ) {
      int ny = many_to_many_nb[jb];  
      float *d = data;
      for ( int i=0; i<nx; ++i, d += K2*nz*2 ) {
	for ( int j=jb*block2; j<(jb*block2+ny); ++j ) {
	  for ( int k=0; k<nz; ++k ) {
	    d[2*(j*nz+k)] = *(md++);
	    d[2*(j*nz+k)+1] = *(md++);
	  }
	}
      }
    }
  }
}

void OptPmeZPencil::many_to_many_recv_untrans() {  
  //printf ("%d, %d In recv untrans\n", thisIndex.x, thisIndex.y);                                                                             
  int block3 = initdata.grid.block3;
  int dim3 = initdata.grid.dim3;
  
  const float *md = many_to_many_data;
  if (single_pencil) {
    float *d = data;
    for (int kb = 0; kb < initdata.zBlocks; kb++ ) {
      d[2*kb]   = *(md++);
      d[2*kb+1] = *(md++);
    }
  }
  else {
    for (int kb = 0; kb < initdata.zBlocks; kb++ ) {
      int nz = many_to_many_nb[kb];
      float *d = data;
      for ( int i=0; i<nx; ++i ) {
	for ( int j=0; j<ny; ++j, d += dim3 ) {
	  for ( int k=kb*block3; k<(kb*block3+nz); ++k ) {
	    d[2*k] = *(md++);
	    d[2*k+1] = *(md++);
	  }
	}
      }
    }
  }
}

void OptPmeZPencil::many_to_many_send_ungrid () {
  int dim3 = initdata.grid.dim3;
  int K3 = initdata.grid.K3;
  int K3_1 = initdata.grid.K3 - 1;
  int k = 0;
  for (int idx = 0; idx < grid_msgs.size(); idx++) {
    int xstart = m2m_recv_grid[idx].xstart; 
    int xlen   = m2m_recv_grid[idx].xlen; 
    int ystart = m2m_recv_grid[idx].ystart; 
    int ylen   = m2m_recv_grid[idx].ylen; 
    int zstart = m2m_recv_grid[idx].zstart; 
    int zlen   = m2m_recv_grid[idx].zlen;   
    
    float *qmsg = m2m_recv_grid[idx].data;
    for ( int i=xstart; i<xstart+xlen; ++i ) {
      for ( int j=ystart; j<ystart+ylen; ++j) {
	float *d = data + (i * ny + j) * dim3; 	

#pragma disjoint (*d, *qmsg)
#pragma unroll (4)
	for ( k=zstart; k<zstart+zlen; ++k ) {
	  int kz = k;
	  //if (kz >= K3) kz -= K3;
	  kz = kz - ((unsigned)(K3_1 - kz)>>31)*K3;
	  //assert (kz >= 0);
	  //assert (kz <  K3);
	  *(qmsg++) = d[kz];
	}
      }
    } 
  }
  
  CmiDirect_manytomany_start (handle, PHASE_UG);
}



void  OptPmeZPencil::initialize_manytomany () {  
  int idx = 0;
  int totalcount = 0;
  for (idx = 0; idx < grid_msgs.size(); idx ++) 
    totalcount += m2m_recv_grid[idx].xlen * m2m_recv_grid[idx].ylen * m2m_recv_grid[idx].zlen;
  many_to_many_gr_data = new float [totalcount];
  
  CkArrayIndex3D aidx (thisIndex.x, thisIndex.y, thisIndex.z);
  cbw_recvgrid.cb = CkCallback (CkIndex_OptPmeZPencil::many_to_manyRecvGrid(NULL), aidx, thisProxy.ckGetArrayID());
  cbw_recvgrid.msg = new (PRIORITY_SIZE) OptPmeDummyMsg;
  cbw_recvgrid.array = this;
  PatchMap *patchMap = PatchMap::Object();
  CmiDirect_manytomany_initialize_recvbase (handle, PHASE_GR, 
					    call_ck_cb_recv_grid, 
					    &cbw_recvgrid, 
					    (char *)many_to_many_gr_data, 
					    patchMap->numPatches(), -1);
  
  CmiDirect_manytomany_initialize_sendbase (handle, PHASE_UG, NULL, NULL, (char *)many_to_many_gr_data, 
					    grid_msgs.size(), thisIndex.x *initdata.yBlocks + thisIndex.y);
  
  int offset = 0;
  for (idx = 0; idx < grid_msgs.size(); idx ++) {
    m2m_recv_grid[idx].data = (float *) ((char *)many_to_many_gr_data + offset);
    int fcount = m2m_recv_grid[idx].xlen * m2m_recv_grid[idx].ylen * m2m_recv_grid[idx].zlen;
    CmiDirect_manytomany_initialize_recv (handle, PHASE_GR, m2m_recv_grid[idx].patchID, 
					  offset, fcount *sizeof(float), patchMap->node(m2m_recv_grid[idx].patchID));
    CmiDirect_manytomany_initialize_send (handle, PHASE_UG, idx, offset, fcount *sizeof(float), 
					  patchMap->node(m2m_recv_grid[idx].patchID));
    offset += fcount * sizeof(float);
  }        

  int zBlocks = initdata.zBlocks;
  int block3 = initdata.grid.block3;
  int dim3 = initdata.grid.dim3;

  //Initialize send trans
  CmiDirect_manytomany_initialize_sendbase (handle, PHASE_YF, NULL, NULL, (char *)many_to_many_data, zBlocks, thisIndex.y);

  //Initialize recv untrans
  cbw_recvuntrans.cb = CkCallback(CkIndex_OptPmeZPencil::many_to_manyRecvUntrans(NULL), aidx, thisProxy.ckGetArrayID());
  cbw_recvuntrans.msg    = new (PRIORITY_SIZE) OptPmeDummyMsg;
  cbw_recvuntrans.array  = this;
  CmiDirect_manytomany_initialize_recvbase (handle, PHASE_ZB, 
					    call_ck_cb_recv_untrans_z, 
					    &cbw_recvuntrans, 
					    (char *)many_to_many_data, 
					    initdata.zBlocks, -1);
  single_pencil = false;
  if (nx == 1 && ny == 1 && zBlocks <= dim3/2 && block3==1)
    single_pencil = true;
  
  for ( int kb=0; kb<zBlocks; ++kb ) {
    int nz = block3;
    if ( (kb+1)*block3 > dim3/2 ) {
      single_pencil = false;
      nz = dim3/2 - kb*block3;
    }
    
    //Initialize send trans
    CkArrayIndex3D index (thisIndex.x,0,kb);
    CProxy_OptPmePencilMapY yproxy (global_map_y);
    int pe = yproxy.ckLocalBranch()->procNum(0, index);
    CmiDirect_manytomany_initialize_send (handle, PHASE_YF, kb, kb*2*nx*ny*block3*sizeof(float), 2*nx*ny*nz*sizeof(float), pe);
    
    //Initialize recv untrans
    CmiDirect_manytomany_initialize_recv (handle, PHASE_ZB, kb, kb*nx*ny*block3*2*sizeof(float), nx*ny*nz*2*sizeof(float), pe);
    many_to_many_nb [kb] = nz;
  }

}

void  OptPmeYPencil::initialize_manytomany () {
  int yBlocks = initdata.yBlocks;
  int block2 = initdata.grid.block2;
  int K2 = initdata.grid.K2;

  //send trans
  CmiDirect_manytomany_initialize_sendbase (handle, PHASE_XF, NULL, NULL, (char *)many_to_many_data, yBlocks, thisIndex.x);
  //send untrans
  CmiDirect_manytomany_initialize_sendbase (handle, PHASE_ZB, NULL, NULL, (char *)many_to_many_data, yBlocks, thisIndex.z);

  //recv trans
  CkArrayIndex3D idx (thisIndex.x, thisIndex.y, thisIndex.z);
  cbw_recvtrans.cb = CkCallback (CkIndex_OptPmeYPencil::many_to_manyRecvTrans(NULL), idx, thisProxy.ckGetArrayID());
  cbw_recvtrans.msg    = new (PRIORITY_SIZE) OptPmeDummyMsg;
  cbw_recvtrans.array  = this;
  CmiDirect_manytomany_initialize_recvbase (handle, PHASE_YF, 
					    call_ck_cb_recv_trans_y, 
					    &cbw_recvtrans, 
					    (char *)many_to_many_data, 
					    initdata.yBlocks, -1);

  //recv untrans
  cbw_recvuntrans.cb = CkCallback(CkIndex_OptPmeYPencil::many_to_manyRecvUntrans(NULL), idx, thisProxy.ckGetArrayID());
  cbw_recvuntrans.msg    = new (PRIORITY_SIZE) OptPmeDummyMsg;
  cbw_recvuntrans.array  = this;
  CmiDirect_manytomany_initialize_recvbase (handle, PHASE_YB, 
					    call_ck_cb_recv_untrans_y, 
					    &cbw_recvuntrans, 
					    (char *)many_to_many_data, 
					    initdata.yBlocks, -1);

  single_pencil = false;
  if (nz == 1 && nx == 1 && yBlocks <= K2 && block2==1)
    single_pencil = true;
  
  for ( int jb=0; jb<yBlocks; ++jb ) {
    int ny = block2;
    if ( (jb+1)*block2 > K2 ) {
      single_pencil = false;
      ny = K2 - jb*block2;
    }

    //send untrans
    CkArrayIndex3D index (thisIndex.x,jb,0);
    CProxy_OptPmePencilMapZ zproxy (global_map_z);
    int pe = zproxy.ckLocalBranch()->procNum(0, index);
    CmiDirect_manytomany_initialize_send (handle, PHASE_ZB, jb, 2*nx*block2*nz*jb*sizeof(float), 2*nx*ny*nz*sizeof(float), pe);

    //recv trans
    CmiDirect_manytomany_initialize_recv (handle, PHASE_YF, jb, jb*nx*block2*nz*2*sizeof(float), nx*ny*nz*2*sizeof(float), pe);
    
    //send trans
    index = CkArrayIndex3D (0,jb,thisIndex.z);
    CProxy_OptPmePencilMapX xproxy (global_map_x);
    pe = xproxy.ckLocalBranch()->procNum(0, index);
    CmiDirect_manytomany_initialize_send (handle, PHASE_XF, jb, 2*nx*block2*nz*jb*sizeof(float), 2*nx*ny*nz*sizeof(float), pe);
    
    //Recv untrans
    CmiDirect_manytomany_initialize_recv (handle, PHASE_YB, jb, jb*nx*block2*nz*2*sizeof(float), nx*ny*nz*2*sizeof(float), pe);

    many_to_many_nb [jb] = ny;
  }
}

void  OptPmeXPencil::initialize_manytomany () {  
  int xBlocks = initdata.xBlocks;
  int block1 = initdata.grid.block1;
  int K1 = initdata.grid.K1;
  
  CmiDirect_manytomany_initialize_sendbase (handle, PHASE_YB, NULL, NULL, (char *)many_to_many_data, xBlocks, thisIndex.y);
  CkArrayIndex3D idx (thisIndex.x, thisIndex.y, thisIndex.z);
  cbw_recvtrans.cb = CkCallback(CkIndex_OptPmeXPencil::many_to_manyRecvTrans(NULL), idx, thisProxy.ckGetArrayID());
  cbw_recvtrans.msg    = new (PRIORITY_SIZE) OptPmeDummyMsg;
  cbw_recvtrans.array  = this;
  CmiDirect_manytomany_initialize_recvbase (handle, PHASE_XF, 
					    call_ck_cb_recv_trans_x, 
					    &cbw_recvtrans,  
					    (char *)many_to_many_data, 
					    initdata.xBlocks, -1);

  single_pencil = false;
  if (ny == 1 && nz == 1 && xBlocks <= K1 && block1==1)
    single_pencil = true;

  for ( int ib=0; ib<xBlocks; ++ib ) {
    int nx = block1;
    if ( (ib+1)*block1 > K1 ) {
      single_pencil = false;
      nx = K1 - ib*block1;
    }
    
    CkArrayIndex3D index (ib,0,thisIndex.z);
    CProxy_OptPmePencilMapY yproxy (global_map_y);
    int pe = yproxy.ckLocalBranch()->procNum(0, index);
    CmiDirect_manytomany_initialize_send (handle, PHASE_YB, ib, 2*block1*ny*nz*ib*sizeof(float), 2*nx*ny*nz*sizeof(float), pe);
    
    CmiDirect_manytomany_initialize_recv (handle, PHASE_XF, ib, ib*block1*ny*nz*2*sizeof(float), nx*ny*nz*2*sizeof(float), pe);
    many_to_many_nb [ib] = nx;
  }
}


////////////////////////////////////////////////////////////////////////////
///////////////////////// End Many To Many Implementation //////////////////
////////////////////////////////////////////////////////////////////////////


#include "PmeFFTLib.def.h"
