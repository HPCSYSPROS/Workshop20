/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include <string.h>
#include "OptPmeRealSpace.h"
#include <assert.h>

OptPmeRealSpace::OptPmeRealSpace(PmeGrid grid,int natoms)
  : myGrid(grid), N(natoms) {
  int order = myGrid.order;
  M = new double[3*N*order];
  dM = new double[3*N*order];
}

OptPmeRealSpace::~OptPmeRealSpace() {
  delete [] M;
  delete [] dM;
}


// this should definitely help the compiler
#define order 4

void OptPmeRealSpace::fill_charges(double **q_arr, PmeParticle p[], int zstart, int zlen) {
  
  int i, j, k, l;
  int stride;
  int K1, K2, K3, dim2, dim3;
  int32  K3_1;
  double *Mi, *dMi;
  Mi = M; dMi = dM;
  K1=myGrid.K1; K2=myGrid.K2; K3=myGrid.K3; dim2=myGrid.dim2; dim3=myGrid.dim3;
  K3_1 = K3 - 1;
  // order = myGrid.order;
  stride = 3*order;

  for (i=0; i<N; i++) {
    double q;
    int32 u1, u2, u2i, u3i;
    double fr1 = p[i].x;
    double fr2 = p[i].y;
    double fr3 = p[i].z;
    q = p[i].cg;
    u1  = (int)fr1;
    u2i = (int)fr2;
    u3i = (int)fr3;
    fr1 -= u1;
    fr2 -= u2i;
    fr3 -= u3i;

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

    u1  -= order;
    u2i -= order;
    u3i -= order;
    u3i += 1;
    for (j=0; j<order; j++) {
      double m1;
      int32 ind1;
      m1 = Mi[j]*q;
      u1++;
      //ind1 = (u1 + (u1 < 0 ? K1 : 0))*dim2;
      ind1 = (u1 + ((unsigned)u1 >>31)*K1)*dim2;
      u2 = u2i;
      for (k=0; k<order; k++) {
        double m1m2;
	int32 ind2;
        m1m2 = m1*Mi[order+k];
	u2++;

	//ind2 = ind1 + (u2 + (u2 < 0 ? K2 : 0));	
	ind2 = ind1 + (u2 + ((unsigned)u2 >>31)*K2);	
	//assert (ind2 < K1 * K2);	
	double *qline = q_arr[ind2];
	
	//if (qline == NULL)
	//printf ("qline [%d, %d] = NULL\n", (u1 + (u1 < 0 ? K1 : 0)), (u2 + (u2 < 0 ? K2 : 0)));
	
	//assert (qline != NULL);
        for (l=0; l<order; l++) {
	  double m3;
	  int32 ind;
	  m3 = Mi[2*order + l];
	  int32 u3 = u3i + l;
          ind = u3 - zstart;
	  //if (ind >= K3) ind -= K3;
	  ind = ind - ((unsigned)(K3_1 - ind) >>31)*K3;

	  //if (ind < 0 || ind >= zlen) printf ("ind >= zlen  (%d,%d,%d)\n", ind, zlen, zstart); 
	  
	  //assert (ind < zlen && ind >= 0);
          qline[ind] += m1m2*m3; 
        }
      }
    }
    Mi += stride;
    dMi += stride;
  }
}

void OptPmeRealSpace::compute_forces(const double * const  * q_arr,
				     const PmeParticle       p[], 
				     Vector                  f[],
				     int                     zstart, 
				     int                     zlen,
				     int                     istart,
				     int                     iend) {
  
  int i, j, k, l, stride;
  double f1, f2, f3;
  double *Mi, *dMi;
  int K1, K2, K3, dim2;
  int32 K3_1;

  if (iend == 0)
    iend = N;

  K1=myGrid.K1; K2=myGrid.K2; K3=myGrid.K3; dim2=myGrid.dim2;
  K3_1 = K3 - 1;
  // order = myGrid.order;
  stride=3*order;
  Mi = M + stride*istart;
  dMi = dM + stride*istart;
 
  for (i=istart; i<iend; i++) {
    double q;
    int32 u1, u2, u2i, u3i;
    q = p[i].cg;
    f1=f2=f3=0.0;
    u1 = (int)(p[i].x);
    u2i = (int)(p[i].y);
    u3i = (int)(p[i].z);
    u1 -= order;
    u2i -= order;
    u3i -= order;
    u3i += 1;
    for (j=0; j<order; j++) {
      double m1, d1;
      int32 ind1;
      m1=Mi[j]*q;
      d1=K1*dMi[j]*q;
      u1++;
      //ind1 = (u1 + (u1 < 0 ? K1 : 0))*dim2; 
      ind1 = (u1 + ((unsigned)u1 >>31)*K1)*dim2;
      u2 = u2i;
      for (k=0; k<order; k++) {
        double m2, d2, m1m2, m1d2, d1m2;
	int32 ind2;
        m2=Mi[order+k];
	d2=K2*dMi[order+k];
	m1m2=m1*m2;
	m1d2=m1*d2;
	d1m2=d1*m2;
	u2++;
	//ind2 = ind1 + (u2 + (u2 < 0 ? K2 : 0));
	ind2 = ind1 + (u2 + ((unsigned)u2 >>31)*K2);	
	const double *qline = q_arr[ind2];
        for (l=0; l<order; l++) {
	  double term, m3, d3;
	  int32 ind;
	  m3=Mi[2*order+l];
	  d3=K3*dMi[2*order+l];
	  int32 u3 = u3i + l;
	  ind = u3 - zstart;	  
	  //if (ind >= K3) ind -= K3;
	  ind = ind - ((unsigned)(K3_1 -ind) >>31)*K3;
	  
	  //if (ind < 0 || ind >= zlen) printf ("ind >= zlen  (%d,%d,%d)\n", ind, zlen, zstart); 	  
	  //assert (ind < zlen && ind >= 0);

	  term = qline[ind];
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
