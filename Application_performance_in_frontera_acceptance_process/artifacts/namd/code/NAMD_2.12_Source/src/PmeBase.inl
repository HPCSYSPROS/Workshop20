/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef PME_BASE_INL__
#define PME_BASE_INL__

#include "PmeBase.h"
#include <math.h>
#include <charm++.h>
#include "Vector.h"
#include "Lattice.h"


static inline void scale_coordinates(PmeParticle p[], int N, Lattice lattice, PmeGrid grid) {
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
  double shift1 = ((K1 + grid.order - 1)/2)/(double)K1;
  double shift2 = ((K2 + grid.order - 1)/2)/(double)K2;
  double shift3 = ((K3 + grid.order - 1)/2)/(double)K3;

  for (int i=0; i<N; i++) {
    double px = p[i].x - ox;
    double py = p[i].y - oy;
    double pz = p[i].z - oz;
    double sx = shift1 + px*r1x + py*r1y + pz*r1z;
    double sy = shift2 + px*r2x + py*r2y + pz*r2z;
    double sz = shift3 + px*r3x + py*r3y + pz*r3z;
    p[i].x = K1 * ( sx - floor(sx) );
    p[i].y = K2 * ( sy - floor(sy) );
    p[i].z = K3 * ( sz - floor(sz) );
    //  Check for rare rounding condition where K * ( 1 - epsilon ) == K
    //  which was observed with g++ on Intel x86 architecture.
    if ( p[i].x == K1 ) p[i].x = 0;
    if ( p[i].y == K2 ) p[i].y = 0;
    if ( p[i].z == K3 ) p[i].z = 0;
  }
}


static inline void scale_forces(Vector f[], int N, Lattice &lattice) {
  Vector recip1 = lattice.a_r();
  Vector recip2 = lattice.b_r();
  Vector recip3 = lattice.c_r();
  double r1x = recip1.x;
  double r1y = recip1.y;
  double r1z = recip1.z;
  double r2x = recip2.x;
  double r2y = recip2.y;
  double r2z = recip2.z;
  double r3x = recip3.x;
  double r3y = recip3.y;
  double r3z = recip3.z;

  for (int i=0; i<N; i++) {
    double f1 = f[i].x;
    double f2 = f[i].y;
    double f3 = f[i].z;
    f[i].x = f1*r1x + f2*r2x + f3*r3x;
    f[i].y = f1*r1y + f2*r2y + f3*r3y;
    f[i].z = f1*r1z + f2*r2z + f3*r3z;
  }
}


template <class REAL>
static inline void compute_b_spline(REAL * __restrict frac, REAL *M, REAL *dM, int order) {
  int j, n;
  REAL x,y,z,x1,y1,z1, div;
  REAL * __restrict Mx, * __restrict My, * __restrict Mz, * __restrict dMx, * __restrict dMy, * __restrict dMz;
  Mx=M-1; My=M+order-1; Mz=M+2*order-1;
  dMx=dM-1; dMy =dM+order-1; dMz=dM+2*order-1;
  x=frac[0];
  y=frac[1];
  z=frac[2];
  x1=1.0-x; y1=1.0-y; z1=1.0-z;
  /* Do n=3 case first */
  Mx[1]=.5*x1*x1;
  Mx[2]=x1*x + .5;
  Mx[3]=0.5*x*x;
  Mx[order]=0.0;
  My[1]=.5*y1*y1;
  My[2]=y1*y + .5;
  My[3]=0.5*y*y;
  My[order]=0.0;
  Mz[1]=.5*z1*z1;
  Mz[2]=z1*z + .5;
  Mz[3]=0.5*z*z;
  Mz[order]=0.0;

  /* Recursively fill in the rest.  */
  for (n=4; n<=order-1; n++) {
    REAL div=1.0/(n-1);
    int j;
    Mx[n] = x*div*Mx[n-1];
    My[n] = y*div*My[n-1];
    Mz[n] = z*div*Mz[n-1];
    for (j=1; j<=n-2; j++) {
      Mx[n-j] = ((x+j)*Mx[n-j-1] + (n-x-j)*Mx[n-j])*div;
      My[n-j] = ((y+j)*My[n-j-1] + (n-y-j)*My[n-j])*div;
      Mz[n-j] = ((z+j)*Mz[n-j-1] + (n-z-j)*Mz[n-j])*div;
    }
    Mx[1] *= (1.0-x)*div;
    My[1] *= (1.0-y)*div;
    Mz[1] *= (1.0-z)*div;
  }
  /* Now do the derivatives  */
  dMx[1]=-Mx[1]; dMy[1]=-My[1]; dMz[1]=-Mz[1];
  for (j=2; j <= order; j++) {
    dMx[j] = Mx[j-1] - Mx[j];
    dMy[j] = My[j-1] - My[j];
    dMz[j] = Mz[j-1] - Mz[j];
  }
  /* Now finish the job!    */
  div=1.0/(order-1);
  Mx[order] = x*div*Mx[order-1];
  My[order] = y*div*My[order-1];
  Mz[order] = z*div*Mz[order-1];
  for (j=1; j<=order-2; j++) {
    Mx[order-j] = ((x+j)*Mx[order-j-1] + (order-x-j)*Mx[order-j])*div;
    My[order-j] = ((y+j)*My[order-j-1] + (order-y-j)*My[order-j])*div;
    Mz[order-j] = ((z+j)*Mz[order-j-1] + (order-z-j)*Mz[order-j])*div;
  }
  Mx[1] *= (1.0-x)*div;
  My[1] *= (1.0-y)*div;
  Mz[1] *= (1.0-z)*div;
}


#endif

