/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

//#include "InfoStream.h"
#include "Settle.h"
#include <string.h>
#include <math.h>
//#include <charm++.h> // for CkPrintf

#if defined(__SSE2__) && ! defined(NAMD_DISABLE_SSE)
#include <emmintrin.h>  // SSE2
#if defined(__INTEL_COMPILER)
#define __align(X) __declspec(align(X) )
#elif defined(__PGI)
#define __align(X)  __attribute__((aligned(X) ))
#define MISSING_mm_cvtsd_f64
#elif defined(__GNUC__)
#define __align(X)  __attribute__((aligned(X) ))
#if (__GNUC__ < 4)
#define MISSING_mm_cvtsd_f64
#endif
#else
#define __align(X) __declspec(align(X) )
#endif
#endif

//
// XXX static and global variables are unsafe for shared memory builds.
// The global and static vars should be eliminated.
// Unfortunately, the routines that use these below are actually
// in use in NAMD.
//

// Initialize various properties of the waters
// settle1() assumes all waters are identical, 
// and will generate bad results if they are not.
void settle1init(BigReal pmO, BigReal pmH, BigReal hhdist, BigReal ohdist,
                 BigReal &mOrmT, BigReal &mHrmT, BigReal &ra,
                 BigReal &rb, BigReal &rc, BigReal &rra) {

    BigReal rmT = 1.0 / (pmO+pmH+pmH);
    mOrmT = pmO * rmT;
    mHrmT = pmH * rmT;
    BigReal t1 = 0.5*pmO/pmH;
    rc = 0.5*hhdist;
    ra = sqrt(ohdist*ohdist-rc*rc)/(1.0+t1);
    rb = t1*ra;
    rra = 1.0 / ra;
}


int settle1(const Vector *ref, Vector *pos, Vector *vel, BigReal invdt,
                 BigReal mOrmT, BigReal mHrmT, BigReal ra,
                 BigReal rb, BigReal rc, BigReal rra) {
#if defined(__SSE2__) && ! defined(NAMD_DISABLE_SSE)
  // SSE acceleration of some of the costly parts of settle using
  // the Intel C/C++ classes.  This implementation uses the SSE units
  // less efficiency than is potentially possible, but in order to do
  // better, the settle algorithm will have to be vectorized and operate
  // on multiple waters at a time.  Doing so could give us the ability to
  // do two (double precison) or four (single precision) waters at a time.
  // This code achieves a modest speedup without the need to reorganize
  // the NAMD structure.  Once we have water molecules sorted in a single
  // block we can do far better.

  // vectors in the plane of the original positions
  Vector b0, c0;

  __m128d REF0xy = _mm_loadu_pd((double *) &ref[0].x);  // ref0.y and ref0.x
  __m128d REF1xy = _mm_loadu_pd((double *) &ref[1].x);  // ref1.y and ref1.x

  __m128d B0xy = _mm_sub_pd(REF1xy, REF0xy);
  _mm_storeu_pd((double *) &b0.x, B0xy);
  b0.z = ref[1].z - ref[0].z;

  __m128d REF2xy = _mm_loadu_pd((double *) &ref[2].x);  // ref2.y and ref2.x

  __m128d C0xy = _mm_sub_pd(REF2xy, REF0xy);
  _mm_storeu_pd((double *) &c0.x, C0xy);
  c0.z = ref[2].z - ref[0].z;

  // new center of mass
  // Vector d0 = pos[0] * mOrmT + ((pos[1] + pos[2]) * mHrmT);
  __align(16) Vector a1;
  __align(16) Vector b1;
  __align(16) Vector c1;
  __align(16) Vector d0;

  __m128d POS1xy = _mm_loadu_pd((double *) &pos[1].x);
  __m128d POS2xy = _mm_loadu_pd((double *) &pos[2].x);
  __m128d PMHrmTxy = _mm_mul_pd(_mm_add_pd(POS1xy, POS2xy), _mm_set1_pd(mHrmT));

  __m128d POS0xy = _mm_loadu_pd((double *) &pos[0].x);
  __m128d PMOrmTxy = _mm_mul_pd(POS0xy, _mm_set1_pd(mOrmT));
  __m128d D0xy = _mm_add_pd(PMOrmTxy, PMHrmTxy);

  d0.z = pos[0].z * mOrmT + ((pos[1].z + pos[2].z) * mHrmT);
  a1.z = pos[0].z - d0.z;
  b1.z = pos[1].z - d0.z;
  c1.z = pos[2].z - d0.z;

  __m128d A1xy = _mm_sub_pd(POS0xy, D0xy);
  _mm_store_pd((double *) &a1.x, A1xy); // must be aligned

  __m128d B1xy = _mm_sub_pd(POS1xy, D0xy);
  _mm_store_pd((double *) &b1.x, B1xy); // must be aligned

  __m128d C1xy = _mm_sub_pd(POS2xy, D0xy);
  _mm_store_pd((double *) &c1.x, C1xy); // must be aligned

  _mm_store_pd((double *) &d0.x, D0xy); // must be aligned
  
  // Vectors describing transformation from original coordinate system to
  // the 'primed' coordinate system as in the diagram.  
  Vector n0 = cross(b0, c0);
  Vector n1 = cross(a1, n0); 
  Vector n2 = cross(n0, n1); 
#else
  // vectors in the plane of the original positions
  Vector b0 = ref[1]-ref[0];
  Vector c0 = ref[2]-ref[0];
  
  // new center of mass
  Vector d0 = pos[0]*mOrmT + ((pos[1] + pos[2])*mHrmT);
 
  Vector a1 = pos[0] - d0;
  Vector b1 = pos[1] - d0;
  Vector c1 = pos[2] - d0;
  
  // Vectors describing transformation from original coordinate system to
  // the 'primed' coordinate system as in the diagram.  
  Vector n0 = cross(b0, c0);
  Vector n1 = cross(a1, n0); 
  Vector n2 = cross(n0, n1); 
#endif

#if defined(__SSE2__) && ! defined(NAMD_DISABLE_SSE) && ! defined(MISSING_mm_cvtsd_f64)
  __m128d l1 = _mm_set_pd(n0.x, n0.y);
  l1 = _mm_mul_pd(l1, l1);
  // n0.x^2 + n0.y^2
  double l1xy0 = _mm_cvtsd_f64(_mm_add_sd(l1, _mm_shuffle_pd(l1, l1, 1)));

  __m128d l3 = _mm_set_pd(n1.y, n1.z);
  l3 = _mm_mul_pd(l3, l3);
  // n1.y^2 + n1.z^2
  double l3yz1 = _mm_cvtsd_f64(_mm_add_sd(l3, _mm_shuffle_pd(l3, l3, 1)));

  __m128d l2 = _mm_set_pd(n1.x, n0.z);
  // len(n1)^2 and len(n0)^2 
  __m128d ts01 = _mm_add_pd(_mm_set_pd(l3yz1, l1xy0), _mm_mul_pd(l2, l2));

  __m128d l4 = _mm_set_pd(n2.x, n2.y);
  l4 = _mm_mul_pd(l4, l4);
  // n2.x^2 + n2.y^2
  double l4xy2 = _mm_cvtsd_f64(_mm_add_sd(l4, _mm_shuffle_pd(l4, l4, 1)));
  double ts2 = l4xy2 + (n2.z * n2.z);              // len(n2)^2

  double  invlens[4];
  // since rsqrt_nr() doesn't work with current compiler
  // this is the next best option 
  static const __m128d fvecd1p0 = _mm_set1_pd(1.0);

  // 1/len(n1) and 1/len(n0)
  __m128d invlen12 = _mm_div_pd(fvecd1p0, _mm_sqrt_pd(ts01));

  // invlens[0]=1/len(n0), invlens[1]=1/len(n1)
  _mm_storeu_pd(invlens, invlen12);

  n0 = n0 * invlens[0];

  // shuffle the order of operations around from the normal algorithm so
  // that we can double pump sqrt() with n2 and cosphi at the same time
  // these components are usually computed down in the canonical water block
  BigReal A1Z = n0 * a1;
  BigReal sinphi = A1Z * rra;
  BigReal tmp = 1.0-sinphi*sinphi;

  __m128d n2cosphi = _mm_sqrt_pd(_mm_set_pd(tmp, ts2));
  // invlens[2] = 1/len(n2), invlens[3] = cosphi
  _mm_storeu_pd(invlens+2, n2cosphi);

  n1 = n1 * invlens[1];
  n2 = n2 * (1.0 / invlens[2]);
  BigReal cosphi = invlens[3];

  b0 = Vector(n1*b0, n2*b0, n0*b0); // note: b0.z is never referenced again
  c0 = Vector(n1*c0, n2*c0, n0*c0); // note: c0.z is never referenced again
 
  b1 = Vector(n1*b1, n2*b1, n0*b1);
  c1 = Vector(n1*c1, n2*c1, n0*c1);

  // now we can compute positions of canonical water 
  BigReal sinpsi = (b1.z - c1.z)/(2.0*rc*cosphi);
  tmp = 1.0-sinpsi*sinpsi;
  BigReal cospsi = sqrt(tmp);
#else
  n0 = n0.unit();
  n1 = n1.unit();
  n2 = n2.unit();

  b0 = Vector(n1*b0, n2*b0, n0*b0); // note: b0.z is never referenced again
  c0 = Vector(n1*c0, n2*c0, n0*c0); // note: c0.z is never referenced again
 
  BigReal A1Z = n0 * a1;
  b1 = Vector(n1*b1, n2*b1, n0*b1);
  c1 = Vector(n1*c1, n2*c1, n0*c1);

  // now we can compute positions of canonical water 
  BigReal sinphi = A1Z * rra;
  BigReal tmp = 1.0-sinphi*sinphi;
  BigReal cosphi = sqrt(tmp);
  BigReal sinpsi = (b1.z - c1.z)/(2.0*rc*cosphi);
  tmp = 1.0-sinpsi*sinpsi;
  BigReal cospsi = sqrt(tmp);
#endif

  BigReal rbphi = -rb*cosphi;
  BigReal tmp1 = rc*sinpsi*sinphi;
  BigReal tmp2 = rc*sinpsi*cosphi;
 
  Vector a2(0, ra*cosphi, ra*sinphi);
  Vector b2(-rc*cospsi, rbphi - tmp1, -rb*sinphi + tmp2);
  Vector c2( rc*cosphi, rbphi + tmp1, -rb*sinphi - tmp2);

  // there are no a0 terms because we've already subtracted the term off 
  // when we first defined b0 and c0.
  BigReal alpha = b2.x*(b0.x - c0.x) + b0.y*b2.y + c0.y*c2.y;
  BigReal beta  = b2.x*(c0.y - b0.y) + b0.x*b2.y + c0.x*c2.y;
  BigReal gama  = b0.x*b1.y - b1.x*b0.y + c0.x*c1.y - c1.x*c0.y;
 
  BigReal a2b2 = alpha*alpha + beta*beta;
  BigReal sintheta = (alpha*gama - beta*sqrt(a2b2 - gama*gama))/a2b2;
  BigReal costheta = sqrt(1.0 - sintheta*sintheta);
  
#if 0
  Vector a3( -a2.y*sintheta, 
              a2.y*costheta,
              a2.z);
  Vector b3(b2.x*costheta - b2.y*sintheta,
              b2.x*sintheta + b2.y*costheta,
              b2.z);
  Vector c3(c2.x*costheta - c2.y*sintheta,
              c2.x*sintheta + c2.y*costheta,
              c2.z);
  
#else
  Vector a3( -a2.y*sintheta, 
              a2.y*costheta,
              A1Z);
  Vector b3(b2.x*costheta - b2.y*sintheta,
              b2.x*sintheta + b2.y*costheta,
              b1.z);
  Vector c3(-b2.x*costheta - c2.y*sintheta,
            -b2.x*sintheta + c2.y*costheta,
              c1.z);

#endif

  // undo the transformation; generate new normal vectors from the transpose.
  Vector m1(n1.x, n2.x, n0.x);
  Vector m2(n1.y, n2.y, n0.y);
  Vector m0(n1.z, n2.z, n0.z);

  pos[0] = Vector(a3*m1, a3*m2, a3*m0) + d0;
  pos[1] = Vector(b3*m1, b3*m2, b3*m0) + d0;
  pos[2] = Vector(c3*m1, c3*m2, c3*m0) + d0;

  // dt can be negative during startup!
  if (invdt != 0) {
    vel[0] = (pos[0]-ref[0])*invdt;
    vel[1] = (pos[1]-ref[1])*invdt;
    vel[2] = (pos[2]-ref[2])*invdt;
  }

  return 0;
}

//
// Settle multiple waters using SIMD
//
template <int veclen>
void settle1_SIMD(const Vector *ref, Vector *pos,
  BigReal mOrmT, BigReal mHrmT, BigReal ra,
  BigReal rb, BigReal rc, BigReal rra) {

  BigReal ref0xt[veclen];
  BigReal ref0yt[veclen];
  BigReal ref0zt[veclen];
  BigReal ref1xt[veclen];
  BigReal ref1yt[veclen];
  BigReal ref1zt[veclen];
  BigReal ref2xt[veclen];
  BigReal ref2yt[veclen];
  BigReal ref2zt[veclen];

  BigReal pos0xt[veclen];
  BigReal pos0yt[veclen];
  BigReal pos0zt[veclen];
  BigReal pos1xt[veclen];
  BigReal pos1yt[veclen];
  BigReal pos1zt[veclen];
  BigReal pos2xt[veclen];
  BigReal pos2yt[veclen];
  BigReal pos2zt[veclen];

  for (int i=0;i < veclen;i++) {
    ref0xt[i] = ref[i*3+0].x;
    ref0yt[i] = ref[i*3+0].y;
    ref0zt[i] = ref[i*3+0].z;
    ref1xt[i] = ref[i*3+1].x;
    ref1yt[i] = ref[i*3+1].y;
    ref1zt[i] = ref[i*3+1].z;
    ref2xt[i] = ref[i*3+2].x;
    ref2yt[i] = ref[i*3+2].y;
    ref2zt[i] = ref[i*3+2].z;

    pos0xt[i] = pos[i*3+0].x;
    pos0yt[i] = pos[i*3+0].y;
    pos0zt[i] = pos[i*3+0].z;
    pos1xt[i] = pos[i*3+1].x;
    pos1yt[i] = pos[i*3+1].y;
    pos1zt[i] = pos[i*3+1].z;
    pos2xt[i] = pos[i*3+2].x;
    pos2yt[i] = pos[i*3+2].y;
    pos2zt[i] = pos[i*3+2].z;
  }

#pragma simd assert
  for (int i=0;i < veclen;i++) {

    BigReal ref0x = ref0xt[i];
    BigReal ref0y = ref0yt[i];
    BigReal ref0z = ref0zt[i];
    BigReal ref1x = ref1xt[i];
    BigReal ref1y = ref1yt[i];
    BigReal ref1z = ref1zt[i];
    BigReal ref2x = ref2xt[i];
    BigReal ref2y = ref2yt[i];
    BigReal ref2z = ref2zt[i];

    BigReal pos0x = pos0xt[i];
    BigReal pos0y = pos0yt[i];
    BigReal pos0z = pos0zt[i];
    BigReal pos1x = pos1xt[i];
    BigReal pos1y = pos1yt[i];
    BigReal pos1z = pos1zt[i];
    BigReal pos2x = pos2xt[i];
    BigReal pos2y = pos2yt[i];
    BigReal pos2z = pos2zt[i];

    // vectors in the plane of the original positions
    BigReal b0x = ref1x - ref0x;
    BigReal b0y = ref1y - ref0y;
    BigReal b0z = ref1z - ref0z;

    BigReal c0x = ref2x - ref0x;
    BigReal c0y = ref2y - ref0y;
    BigReal c0z = ref2z - ref0z;
    
    // new center of mass
    BigReal d0x = pos0x*mOrmT + ((pos1x + pos2x)*mHrmT);
    BigReal d0y = pos0y*mOrmT + ((pos1y + pos2y)*mHrmT);
    BigReal d0z = pos0z*mOrmT + ((pos1z + pos2z)*mHrmT);
   
    BigReal a1x = pos0x - d0x;
    BigReal a1y = pos0y - d0y;
    BigReal a1z = pos0z - d0z;

    BigReal b1x = pos1x - d0x;
    BigReal b1y = pos1y - d0y;
    BigReal b1z = pos1z - d0z;

    BigReal c1x = pos2x - d0x;
    BigReal c1y = pos2y - d0y;
    BigReal c1z = pos2z - d0z;
    
    // Vectors describing transformation from original coordinate system to
    // the 'primed' coordinate system as in the diagram.
    // n0 = b0 x c0
    BigReal n0x = b0y*c0z-c0y*b0z;
    BigReal n0y = c0x*b0z-b0x*c0z;
    BigReal n0z = b0x*c0y-c0x*b0y;

    // n1 = a1 x n0
    BigReal n1x = a1y*n0z-n0y*a1z;
    BigReal n1y = n0x*a1z-a1x*n0z;
    BigReal n1z = a1x*n0y-n0x*a1y;

    // n2 = n0 x n1
    BigReal n2x = n0y*n1z-n1y*n0z;
    BigReal n2y = n1x*n0z-n0x*n1z;
    BigReal n2z = n0x*n1y-n1x*n0y;

    // Normalize n0
    BigReal n0inv = 1.0/sqrt(n0x*n0x + n0y*n0y + n0z*n0z);
    n0x *= n0inv;
    n0y *= n0inv;
    n0z *= n0inv;

    BigReal n1inv = 1.0/sqrt(n1x*n1x + n1y*n1y + n1z*n1z);
    n1x *= n1inv;
    n1y *= n1inv;
    n1z *= n1inv;

    BigReal n2inv = 1.0/sqrt(n2x*n2x + n2y*n2y + n2z*n2z);
    n2x *= n2inv;
    n2y *= n2inv;
    n2z *= n2inv;

    //b0 = Vector(n1*b0, n2*b0, n0*b0); // note: b0.z is never referenced again
    BigReal n1b0 = n1x*b0x + n1y*b0y + n1z*b0z;
    BigReal n2b0 = n2x*b0x + n2y*b0y + n2z*b0z;

    //c0 = Vector(n1*c0, n2*c0, n0*c0); // note: c0.z is never referenced again
    BigReal n1c0 = n1x*c0x + n1y*c0y + n1z*c0z;
    BigReal n2c0 = n2x*c0x + n2y*c0y + n2z*c0z;
   
    BigReal A1Z = n0x*a1x + n0y*a1y + n0z*a1z;
    
    //b1 = Vector(n1*b1, n2*b1, n0*b1);
    BigReal n1b1 = n1x*b1x + n1y*b1y + n1z*b1z;
    BigReal n2b1 = n2x*b1x + n2y*b1y + n2z*b1z;
    BigReal n0b1 = n0x*b1x + n0y*b1y + n0z*b1z;

    //c1 = Vector(n1*c1, n2*c1, n0*c1);
    BigReal n1c1 = n1x*c1x + n1y*c1y + n1z*c1z;
    BigReal n2c1 = n2x*c1x + n2y*c1y + n2z*c1z;
    BigReal n0c1 = n0x*c1x + n0y*c1y + n0z*c1z;

    // now we can compute positions of canonical water 
    BigReal sinphi = A1Z * rra;
    BigReal tmp = 1.0-sinphi*sinphi;
    BigReal cosphi = sqrt(tmp);
    BigReal sinpsi = (n0b1 - n0c1)/(2.0*rc*cosphi);
    tmp = 1.0-sinpsi*sinpsi;
    BigReal cospsi = sqrt(tmp);

    BigReal rbphi = -rb*cosphi;
    BigReal tmp1 = rc*sinpsi*sinphi;
    BigReal tmp2 = rc*sinpsi*cosphi;
   
    //Vector a2(0, ra*cosphi, ra*sinphi);
    BigReal a2y = ra*cosphi;

    //Vector b2(-rc*cospsi, rbphi - tmp1, -rb*sinphi + tmp2);
    BigReal b2x = -rc*cospsi;
    BigReal b2y = rbphi - tmp1;

    //Vector c2( rc*cosphi, rbphi + tmp1, -rb*sinphi - tmp2);
    BigReal c2y = rbphi + tmp1;

    // there are no a0 terms because we've already subtracted the term off 
    // when we first defined b0 and c0.
    BigReal alpha = b2x*(n1b0 - n1c0) + n2b0*b2y + n2c0*c2y;
    BigReal beta  = b2x*(n2c0 - n2b0) + n1b0*b2y + n1c0*c2y;
    BigReal gama  = n1b0*n2b1 - n1b1*n2b0 + n1c0*n2c1 - n1c1*n2c0;
   
    BigReal a2b2 = alpha*alpha + beta*beta;
    BigReal sintheta = (alpha*gama - beta*sqrt(a2b2 - gama*gama))/a2b2;
    BigReal costheta = sqrt(1.0 - sintheta*sintheta);
    
    //Vector a3( -a2y*sintheta, 
    //            a2y*costheta,
    //            A1Z);
    BigReal a3x = -a2y*sintheta;
    BigReal a3y = a2y*costheta;
    BigReal a3z = A1Z;

    // Vector b3(b2x*costheta - b2y*sintheta,
    //             b2x*sintheta + b2y*costheta,
    //             n0b1);
    BigReal b3x = b2x*costheta - b2y*sintheta;
    BigReal b3y = b2x*sintheta + b2y*costheta;
    BigReal b3z = n0b1;

    // Vector c3(-b2x*costheta - c2y*sintheta,
    //           -b2x*sintheta + c2y*costheta,
    //             n0c1);
    BigReal c3x = -b2x*costheta - c2y*sintheta;
    BigReal c3y = -b2x*sintheta + c2y*costheta;
    BigReal c3z = n0c1;

    // undo the transformation; generate new normal vectors from the transpose.
    // Vector m1(n1.x, n2.x, n0.x);
    BigReal m1x = n1x;
    BigReal m1y = n2x;
    BigReal m1z = n0x;

    // Vector m2(n1.y, n2.y, n0.y);
    BigReal m2x = n1y;
    BigReal m2y = n2y;
    BigReal m2z = n0y;

    // Vector m0(n1.z, n2.z, n0.z);
    BigReal m0x = n1z;
    BigReal m0y = n2z;
    BigReal m0z = n0z;

    //pos[i*3+0] = Vector(a3*m1, a3*m2, a3*m0) + d0;
    pos0x = a3x*m1x + a3y*m1y + a3z*m1z + d0x;
    pos0y = a3x*m2x + a3y*m2y + a3z*m2z + d0y;
    pos0z = a3x*m0x + a3y*m0y + a3z*m0z + d0z;

    // pos[i*3+1] = Vector(b3*m1, b3*m2, b3*m0) + d0;
    pos1x = b3x*m1x + b3y*m1y + b3z*m1z + d0x;
    pos1y = b3x*m2x + b3y*m2y + b3z*m2z + d0y;
    pos1z = b3x*m0x + b3y*m0y + b3z*m0z + d0z;

    // pos[i*3+2] = Vector(c3*m1, c3*m2, c3*m0) + d0;
    pos2x = c3x*m1x + c3y*m1y + c3z*m1z + d0x;
    pos2y = c3x*m2x + c3y*m2y + c3z*m2z + d0y;
    pos2z = c3x*m0x + c3y*m0y + c3z*m0z + d0z;

    pos0xt[i] = pos0x;
    pos0yt[i] = pos0y;
    pos0zt[i] = pos0z;
    pos1xt[i] = pos1x;
    pos1yt[i] = pos1y;
    pos1zt[i] = pos1z;
    pos2xt[i] = pos2x;
    pos2yt[i] = pos2y;
    pos2zt[i] = pos2z;
  }

  for (int i=0;i < veclen;i++) {
    pos[i*3+0].x = pos0xt[i];
    pos[i*3+0].y = pos0yt[i];
    pos[i*3+0].z = pos0zt[i];
    pos[i*3+1].x = pos1xt[i];
    pos[i*3+1].y = pos1yt[i];
    pos[i*3+1].z = pos1zt[i];
    pos[i*3+2].x = pos2xt[i];
    pos[i*3+2].y = pos2yt[i];
    pos[i*3+2].z = pos2zt[i];
  }

}

//
// Rattle pair of atoms
//
template <int veclen>
void rattlePair(const RattleParam* rattleParam,
  const BigReal *refx, const BigReal *refy, const BigReal *refz,
  BigReal *posx, BigReal *posy, BigReal *posz) {

  int a = rattleParam[0].ia;
  int b = rattleParam[0].ib;
  BigReal pabx = posx[a] - posx[b];
  BigReal paby = posy[a] - posy[b];
  BigReal pabz = posz[a] - posz[b];
  BigReal pabsq = pabx*pabx + paby*paby + pabz*pabz;
  BigReal rabsq = rattleParam[0].dsq;
  BigReal diffsq = rabsq - pabsq;
  BigReal rabx = refx[a] - refx[b];
  BigReal raby = refy[a] - refy[b];
  BigReal rabz = refz[a] - refz[b];

  BigReal refsq = rabx*rabx + raby*raby + rabz*rabz;

  BigReal rpab = rabx*pabx + raby*paby + rabz*pabz;

  BigReal rma = rattleParam[0].rma;
  BigReal rmb = rattleParam[0].rmb;

  BigReal gab = (-rpab + sqrt(rpab*rpab + refsq*diffsq))/(refsq*(rma + rmb));

  BigReal dpx = rabx * gab;
  BigReal dpy = raby * gab;
  BigReal dpz = rabz * gab;
  posx[a] += rma * dpx;
  posy[a] += rma * dpy;
  posz[a] += rma * dpz;
  posx[b] -= rmb * dpx;
  posy[b] -= rmb * dpy;
  posz[b] -= rmb * dpz;

}

void rattleN(const int icnt, const RattleParam* rattleParam,
  const BigReal *refx, const BigReal *refy, const BigReal *refz,
  BigReal *posx, BigReal *posy, BigReal *posz,
  const BigReal tol2, const int maxiter,
  bool& done, bool& consFailure) {

  for (int iter = 0; iter < maxiter; ++iter ) {
    done = true;
    consFailure = false;
    for (int i = 0; i < icnt; ++i ) {
      int a = rattleParam[i].ia;
      int b = rattleParam[i].ib;
      BigReal pabx = posx[a] - posx[b];
      BigReal paby = posy[a] - posy[b];
      BigReal pabz = posz[a] - posz[b];
      BigReal pabsq = pabx*pabx + paby*paby + pabz*pabz;
      BigReal rabsq = rattleParam[i].dsq;
      BigReal diffsq = rabsq - pabsq;
      if ( fabs(diffsq) > (rabsq * tol2) ) {
        BigReal rabx = refx[a] - refx[b];
        BigReal raby = refy[a] - refy[b];
        BigReal rabz = refz[a] - refz[b];
        BigReal rpab = rabx*pabx + raby*paby + rabz*pabz;
        if ( rpab < ( rabsq * 1.0e-6 ) ) {
          done = false;
          consFailure = true;
          continue;
        }
        BigReal rma = rattleParam[i].rma;
        BigReal rmb = rattleParam[i].rmb;
        BigReal gab = diffsq / ( 2.0 * ( rma + rmb ) * rpab );
        BigReal dpx = rabx * gab;
        BigReal dpy = raby * gab;
        BigReal dpz = rabz * gab;
        posx[a] += rma * dpx;
        posy[a] += rma * dpy;
        posz[a] += rma * dpz;
        posx[b] -= rmb * dpx;
        posy[b] -= rmb * dpy;
        posz[b] -= rmb * dpz;
        done = false;
      }
    }
    if ( done ) break;
  }

}

//
// Explicit instances of templated methods
//
template void rattlePair<1>(const RattleParam* rattleParam,
  const BigReal *refx, const BigReal *refy, const BigReal *refz,
  BigReal *posx, BigReal *posy, BigReal *posz);
template void settle1_SIMD<2>(const Vector *ref, Vector *pos,
  BigReal mOrmT, BigReal mHrmT, BigReal ra,
  BigReal rb, BigReal rc, BigReal rra);
template void settle1_SIMD<1>(const Vector *ref, Vector *pos,
  BigReal mOrmT, BigReal mHrmT, BigReal ra,
  BigReal rb, BigReal rc, BigReal rra);

static int settlev(const Vector *pos, BigReal ma, BigReal mb, Vector *vel,
				   BigReal dt, Tensor *virial) {
  
  Vector rAB = pos[1]-pos[0];
  Vector rBC = pos[2]-pos[1];
  Vector rCA = pos[0]-pos[2];
 
  Vector AB = rAB.unit();
  Vector BC = rBC.unit();
  Vector CA = rCA.unit();
  
  BigReal cosA = -AB * CA;
  BigReal cosB = -BC * AB;
  BigReal cosC = -CA * BC;

  BigReal vab = (vel[1]-vel[0])*AB;
  BigReal vbc = (vel[2]-vel[1])*BC;
  BigReal vca = (vel[0]-vel[2])*CA;

  BigReal mab = ma+mb;
  
  BigReal d = (2*mab*mab + 2*ma*mb*cosA*cosB*cosC - 2*mb*mb*cosA*cosA
               - ma*mab*(cosB*cosB + cosC*cosC))*0.5/mb;

  BigReal tab = (vab*(2*mab - ma*cosC*cosC) +
                vbc*(mb*cosC*cosA - mab*cosB) +
                vca*(ma*cosB*cosC - 2*mb*cosA))*ma/d;
            
  BigReal tbc = (vbc*(mab*mab - mb*mb*cosA*cosA) +
                vca*ma*(mb*cosA*cosB - mab*cosC) +
                vab*ma*(mb*cosC*cosA - mab*cosB))/d;
  
  BigReal tca = (vca*(2*mab - ma*cosB*cosB) +
                vab*(ma*cosB*cosC - 2*mb*cosA) +
                vbc*(mb*cosA*cosB - mab*cosC))*ma/d;
 
  Vector ga = tab*AB - tca*CA;
  Vector gb = tbc*BC - tab*AB;
  Vector gc = tca*CA - tbc*BC;
#if 0
  if (virial) {
    *virial += 0.5*outer(tab, rAB)/dt;
    *virial += 0.5*outer(tbc, rBC)/dt;
    *virial += 0.5*outer(tca, rCA)/dt;
  }
#endif
  vel[0] += (0.5/ma)*ga;
  vel[1] += (0.5/mb)*gb;
  vel[2] += (0.5/mb)*gc;

  return 0;
}

int settle2(BigReal mO, BigReal mH, const Vector *pos,
                  Vector *vel, BigReal dt, Tensor *virial) {

  settlev(pos, mO, mH, vel, dt, virial);
  return 0;
}

