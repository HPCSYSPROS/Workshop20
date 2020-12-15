/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include <string.h>
#include "PmeBase.inl"
#include "PmeRealSpace.h"
#include "Node.h"
#include "SimParameters.h"

PmeRealSpace::PmeRealSpace(PmeGrid grid)
  : myGrid(grid) {
}

PmeRealSpace::~PmeRealSpace() {
}

void PmeRealSpace::set_num_atoms(int natoms) {
  N = natoms;
  int order = myGrid.order;
  M_alloc.resize(3*N*order);
  M = M_alloc.begin();
  dM_alloc.resize(3*N*order);
  dM = dM_alloc.begin();
}

template <int order>
void PmeRealSpace::fill_b_spline(PmeParticle p[]) {
  float fr[3]; 
  float *Mi, *dMi;
  int i, stride;

  stride = 3*order;
  Mi = M; dMi = dM;
  for (i=0; i<N; i++) {
    fr[0] = (float)(p[i].x - (double)(int)(p[i].x));  // subtract in double precision
    fr[1] = (float)(p[i].y - (double)(int)(p[i].y));
    fr[2] = (float)(p[i].z - (double)(int)(p[i].z));
    compute_b_spline(fr, Mi, dMi, order);
    Mi += stride;
    dMi += stride;
  }
}

void PmeRealSpace::fill_charges(float **q_arr, float **q_arr_list, int &q_arr_count, 
                       int &stray_count, char *f_arr, char *fz_arr, PmeParticle p[]) {

  switch (myGrid.order) {
  case 4:
    fill_charges_order4(q_arr, q_arr_list, q_arr_count, stray_count, f_arr, fz_arr, p);
    break;
  case 6:
    fill_charges_order<6>(q_arr, q_arr_list, q_arr_count, stray_count, f_arr, fz_arr, p);
    break;
  case 8:
    fill_charges_order<8>(q_arr, q_arr_list, q_arr_count, stray_count, f_arr, fz_arr, p);
    break;
  case 10:
    fill_charges_order<10>(q_arr, q_arr_list, q_arr_count, stray_count, f_arr, fz_arr, p);
    break;
  default: NAMD_die("unsupported PMEInterpOrder");
  }

}
  
template <int order>
void PmeRealSpace::fill_charges_order(float **q_arr, float **q_arr_list, int &q_arr_count, 
                       int &stray_count, char *f_arr, char *fz_arr, PmeParticle p[]) {
  
  int i, j, k, l;
  int stride;
  int K1, K2, K3, dim2, dim3;

  if ( order != myGrid.order ) NAMD_bug("fill_charges_order template mismatch");

  if ( order == 4 ) {
    fill_charges_order4(q_arr, q_arr_list, q_arr_count, stray_count, f_arr, fz_arr, p);
    return;
  }

  float * __restrict Mi;
  Mi = M;
  K1=myGrid.K1; K2=myGrid.K2; K3=myGrid.K3; dim2=myGrid.dim2; dim3=myGrid.dim3;
  stride = 3*order;

  fill_b_spline<order>(p);

  for (i=0; i<N; i++) {
    float q;
    int u1, u2, u2i, u3i;
    q = p[i].cg;
    u1 = (int)(p[i].x);
    u2i = (int)(p[i].y);
    u3i = (int)(p[i].z);
    u1 -= order;
    u2i -= order;
    u3i -= order;
    u3i += 1;
    if ( u3i < 0 ) u3i += K3;
    for (j=0; j<order; j++) {
      float m1;
      int ind1;
      m1 = Mi[j]*q;
      u1++;
      ind1 = (u1 + (u1 < 0 ? K1 : 0))*dim2;
      u2 = u2i;
      for (k=0; k<order; k++) {
        float m1m2;
	int ind2;
        m1m2 = m1*Mi[order+k];
	u2++;
	ind2 = ind1 + (u2 + (u2 < 0 ? K2 : 0));
	float * __restrict qline = q_arr[ind2];
	if ( ! qline ) {
          if ( f_arr[ind2] ) {
	    f_arr[ind2] = 3;
            ++stray_count;
            continue;
          }
	  qline = q_arr[ind2] = q_arr_list[q_arr_count++]
					= new float[K3+order-1];
	  memset( (void*) qline, 0, (K3+order-1) * sizeof(float) );
	}
	f_arr[ind2] = 1;
        for (l=0; l<order; l++) {
          qline[u3i+l] += m1m2 * Mi[2*order + l];
        }
      }
    }
    Mi += stride;
    for (l=0; l<order; l++) {
      int u3 = u3i + l;
      int ind = u3 + (u3 < 0 ? K3 : 0);
      fz_arr[ind] = 1;
    }
  }
}

void PmeRealSpace::compute_forces(const float * const *q_arr,
				const PmeParticle p[], Vector f[]) {

  switch (myGrid.order) {
  case 4:
    compute_forces_order4(q_arr, p, f);
    break;
  case 6:
    compute_forces_order<6>(q_arr, p, f);
    break;
  case 8:
    compute_forces_order<8>(q_arr, p, f);
    break;
  case 10:
    compute_forces_order<10>(q_arr, p, f);
    break;
  default: NAMD_die("unsupported PMEInterpOrder");
  }

}
  
template <int order>
void PmeRealSpace::compute_forces_order(const float * const *q_arr,
				const PmeParticle p[], Vector f[]) {
  
  int i, j, k, l, stride;
  float f1, f2, f3;
  float *Mi, *dMi;
  int K1, K2, K3, dim2;

  if ( order != myGrid.order ) NAMD_bug("compute_forces_order template mismatch");

  if ( order == 4 ) {
    compute_forces_order4(q_arr, p, f);
    return;
  }

  K1=myGrid.K1; K2=myGrid.K2; K3=myGrid.K3; dim2=myGrid.dim2;
  stride=3*order;
  Mi = M; dMi = dM;
 
  for (i=0; i<N; i++) {
    float q;
    int u1, u2, u2i, u3i;
    q = p[i].cg;
    f1=f2=f3=0.0;
    u1 = (int)(p[i].x);
    u2i = (int)(p[i].y);
    u3i = (int)(p[i].z);
    u1 -= order;
    u2i -= order;
    u3i -= order;
    u3i += 1;
    if ( u3i < 0 ) u3i += K3;
    for (j=0; j<order; j++) {
      float m1, d1;
      int ind1;
      m1=Mi[j]*q;
      d1=K1*dMi[j]*q;
      u1++;
      ind1 = (u1 + (u1 < 0 ? K1 : 0))*dim2; 
      u2 = u2i;
      for (k=0; k<order; k++) {
        float m2, d2, m1m2, m1d2, d1m2;
	int ind2;
        m2=Mi[order+k];
	d2=K2*dMi[order+k];
	m1m2=m1*m2;
	m1d2=m1*d2;
	d1m2=d1*m2;
	u2++;
	ind2 = ind1 + (u2 + (u2 < 0 ? K2 : 0));
	const float *qline = q_arr[ind2];
	if ( ! qline ) continue;
        for (l=0; l<order; l++) {
	  float term, m3, d3;
	  m3=Mi[2*order+l];
	  d3=K3*dMi[2*order+l];
	  term = qline[u3i+l];
	  f1 -= d1m2 * m3 * term;
	  f2 -= m1d2 * m3 * term;
	  f3 -= m1m2 * d3 * term;
        }
      }
    }
    Mi += stride;
    dMi += stride;
    f[i].x = f1;
    f[i].y = f2;
    f[i].z = f3;
  }
}

// this should definitely help the compiler
#define order 4

void PmeRealSpace::fill_charges_order4(float **q_arr, float **q_arr_list, int &q_arr_count, 
                       int &stray_count, char *f_arr, char *fz_arr, PmeParticle p[]) {
  
  int i, j, k, l;
  int stride;
  int K1, K2, K3, dim2, dim3;
  float * __restrict Mi, * __restrict dMi;
  Mi = M; dMi = dM;
  K1=myGrid.K1; K2=myGrid.K2; K3=myGrid.K3; dim2=myGrid.dim2; dim3=myGrid.dim3;
  // order = myGrid.order;
  stride = 3*order;

  // fill_b_spline(p);   leave this off!  replaced with code below

#ifdef NAMD_CUDA
  for ( int istart = 0; istart < N; istart += 1000 ) {
   int iend = istart + 1000;
   if ( iend > N ) iend = N;
   CmiNetworkProgress();
   for (i=istart; i<iend; i++) {
#else
  {
   for (i=0; i<N; i++) {
#endif
    float q;
    int u1, u2, u2i, u3i;
    u1 = (int)p[i].x;
    u2i = (int)p[i].y;
    u3i = (int)p[i].z;
    float fr1 = (float)(p[i].x - (double)u1);  // subtract in double precision
    float fr2 = (float)(p[i].y - (double)u2i);
    float fr3 = (float)(p[i].z - (double)u3i);
    q = p[i].cg;

    // calculate b_spline for order = 4
    Mi[0] = ( ( (-1./6.) * fr1 + 0.5 ) * fr1 - 0.5 ) * fr1 + (1./6.);
    Mi[1] = ( ( 0.5 * fr1 - 1.0 ) * fr1 ) * fr1 + (2./3.);
    Mi[2] = ( ( -0.5 * fr1 + 0.5 ) * fr1 + 0.5 ) * fr1 + (1./6.);
    Mi[3] = (1./6.) * fr1 * fr1 * fr1;
    dMi[0] = ( -0.5 * fr1 + 1.0 )* fr1 - 0.5;
    dMi[1] = ( 1.5 * fr1 - 2.0 ) * fr1;
    dMi[2] = ( -1.5 * fr1 + 1.0 ) * fr1 + 0.5;
    dMi[3] = 0.5 * fr1 * fr1;
    Mi[4] = ( ( (-1./6.) * fr2 + 0.5 ) * fr2 - 0.5 ) * fr2 + (1./6.);
    Mi[5] = ( ( 0.5 * fr2 - 1.0 ) * fr2 ) * fr2 + (2./3.);
    Mi[6] = ( ( -0.5 * fr2 + 0.5 ) * fr2 + 0.5 ) * fr2 + (1./6.);
    Mi[7] = (1./6.) * fr2 * fr2 * fr2;
    dMi[4] = ( -0.5 * fr2 + 1.0 )* fr2 - 0.5;
    dMi[5] = ( 1.5 * fr2 - 2.0 ) * fr2;
    dMi[6] = ( -1.5 * fr2 + 1.0 ) * fr2 + 0.5;
    dMi[7] = 0.5 * fr2 * fr2;
    Mi[8] = ( ( (-1./6.) * fr3 + 0.5 ) * fr3 - 0.5 ) * fr3 + (1./6.);
    Mi[9] = ( ( 0.5 * fr3 - 1.0 ) * fr3 ) * fr3 + (2./3.);
    Mi[10] = ( ( -0.5 * fr3 + 0.5 ) * fr3 + 0.5 ) * fr3 + (1./6.);
    Mi[11] = (1./6.) * fr3 * fr3 * fr3;
    dMi[8] = ( -0.5 * fr3 + 1.0 )* fr3 - 0.5;
    dMi[9] = ( 1.5 * fr3 - 2.0 ) * fr3;
    dMi[10] = ( -1.5 * fr3 + 1.0 ) * fr3 + 0.5;
    dMi[11] = 0.5 * fr3 * fr3;

    u1 -= order;
    u2i -= order;
    u3i -= order;
    u3i += 1;
    if ( u3i < 0 ) u3i += K3;
    for (j=0; j<order; j++) {
      float m1;
      int ind1;
      m1 = Mi[j]*q;
      u1++;
      ind1 = (u1 + (u1 < 0 ? K1 : 0))*dim2;
      u2 = u2i;
      for (k=0; k<order; k++) {
        float m1m2;
	int ind2;
        m1m2 = m1*Mi[order+k];
	u2++;
	ind2 = ind1 + (u2 + (u2 < 0 ? K2 : 0));
	float * __restrict qline = q_arr[ind2];
	if ( ! qline ) {
          if ( f_arr[ind2] ) {
	    f_arr[ind2] = 3;
            ++stray_count;
            continue;
          }
	  qline = q_arr[ind2] = q_arr_list[q_arr_count++]
					= new float[K3+order-1];
	  memset( (void*) qline, 0, (K3+order-1) * sizeof(float) );
	}
	f_arr[ind2] = 1;
        for (l=0; l<order; l++) {
          qline[u3i+l] += m1m2 * Mi[2*order + l];
        }
      }
    }
    Mi += stride;
    dMi += stride;
    for (l=0; l<order; l++) {
      int u3 = u3i + l;
      int ind = u3 + (u3 < 0 ? K3 : 0);
      fz_arr[ind] = 1;
    }
   }
  }
}

static inline void compute_forces_order4_helper(int first, int last, void *result, int paraNum, void *param){
    void **params = (void **)param;
    PmeRealSpace *rs = (PmeRealSpace *)params[0];
    const float * const *q_arr = (const float * const *)params[1];
    const PmeParticle *p = (const PmeParticle *)params[2];
    Vector *f = (Vector *)params[3];
    rs->compute_forces_order4_partial(first, last, q_arr, p, f);
}

void PmeRealSpace::compute_forces_order4(const float * const *q_arr,
				const PmeParticle p[], Vector f[]) {
  
  int i, j, k, l, stride;
  float f1, f2, f3;
  float *Mi, *dMi;
  int K1, K2, K3, dim2;

#if     CMK_SMP && USE_CKLOOP
  int useCkLoop = Node::Object()->simParameters->useCkLoop;
  if(useCkLoop>=CKLOOP_CTRL_PME_UNGRIDCALC){          
      //compute_forces_order4_partial(0, N-1, q_arr, p, f);
      void *params[] = {(void *)this, (void *)q_arr, (void *)p, (void *)f};
      CkLoop_Parallelize(compute_forces_order4_helper, 4, (void *)params, CkMyNodeSize(), 0, N-1);
      return;
  }
#endif
  K1=myGrid.K1; K2=myGrid.K2; K3=myGrid.K3; dim2=myGrid.dim2;
  // order = myGrid.order;
  stride=3*order;
  Mi = M; dMi = dM;
 
  for (i=0; i<N; i++) {
    float q;
    int u1, u2, u2i, u3i;
    q = p[i].cg;
    f1=f2=f3=0.0;
    u1 = (int)(p[i].x);
    u2i = (int)(p[i].y);
    u3i = (int)(p[i].z);
    u1 -= order;
    u2i -= order;
    u3i -= order;
    u3i += 1;
    if ( u3i < 0 ) u3i += K3;
    for (j=0; j<order; j++) {
      float m1, d1;
      int ind1;
      m1=Mi[j]*q;
      d1=K1*dMi[j]*q;
      u1++;
      ind1 = (u1 + (u1 < 0 ? K1 : 0))*dim2; 
      u2 = u2i;
      for (k=0; k<order; k++) {
        float m2, d2, m1m2, m1d2, d1m2;
	int ind2;
        m2=Mi[order+k];
	d2=K2*dMi[order+k];
	m1m2=m1*m2;
	m1d2=m1*d2;
	d1m2=d1*m2;
	u2++;
	ind2 = ind1 + (u2 + (u2 < 0 ? K2 : 0));
	const float *qline = q_arr[ind2];
	if ( ! qline ) continue;
        for (l=0; l<order; l++) {
	  float term, m3, d3;
	  m3=Mi[2*order+l];
	  d3=K3*dMi[2*order+l];
	  term = qline[u3i+l];
	  f1 -= d1m2 * m3 * term;
	  f2 -= m1d2 * m3 * term;
	  f3 -= m1m2 * d3 * term;
        }
      }
    }
    Mi += stride;
    dMi += stride;
    f[i].x = f1;
    f[i].y = f2;
    f[i].z = f3;
  }
}

void PmeRealSpace::compute_forces_order4_partial(int first, int last, 
                const float * const *q_arr,
				const PmeParticle p[], Vector f[]) {
  
  int i, j, k, l, stride;
  float f1, f2, f3;
  float *Mi, *dMi;
  int K1, K2, K3, dim2;

  K1=myGrid.K1; K2=myGrid.K2; K3=myGrid.K3; dim2=myGrid.dim2;
  // order = myGrid.order;
  stride=3*order;
  Mi = M; dMi = dM;
 
  for (i=first; i<=last; i++) {
    Mi = M + i*stride;
    dMi = dM + i*stride;
    float q;
    int u1, u2, u2i, u3i;
    q = p[i].cg;
    f1=f2=f3=0.0;
    u1 = (int)(p[i].x);
    u2i = (int)(p[i].y);
    u3i = (int)(p[i].z);
    u1 -= order;
    u2i -= order;
    u3i -= order;
    u3i += 1;
    if ( u3i < 0 ) u3i += K3;
    for (j=0; j<order; j++) {
      float m1, d1;
      int ind1;
      m1=Mi[j]*q;
      d1=K1*dMi[j]*q;
      u1++;
      ind1 = (u1 + (u1 < 0 ? K1 : 0))*dim2; 
      u2 = u2i;
      for (k=0; k<order; k++) {
        float m2, d2, m1m2, m1d2, d1m2;
	int ind2;
        m2=Mi[order+k];
	d2=K2*dMi[order+k];
	m1m2=m1*m2;
	m1d2=m1*d2;
	d1m2=d1*m2;
	u2++;
	ind2 = ind1 + (u2 + (u2 < 0 ? K2 : 0));
	const float *qline = q_arr[ind2];
	if ( ! qline ) continue;
        for (l=0; l<order; l++) {
	  float term, m3, d3;
	  m3=Mi[2*order+l];
	  d3=K3*dMi[2*order+l];
	  term = qline[u3i+l];
	  f1 -= d1m2 * m3 * term;
	  f2 -= m1d2 * m3 * term;
	  f3 -= m1m2 * d3 * term;
        }
      }
    }
    f[i].x = f1;
    f[i].y = f2;
    f[i].z = f3;
  }
}
