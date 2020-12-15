#include <cmath>
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

#include "GRHydro_Reconstruct_drv_impl.hh"
#include "GRHydro_MP5Reconstruct.hh"

using namespace std;

template <typename T> static inline T SQR (T const & x) { return x*x; }

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

#define MIN3(a,b,c) (MIN(a, MIN(b,c)))
#define MIN4(a,b,c,d) (MIN(a, MIN(b,MIN(c,d))))
#define MAX3(a,b,c) (MAX(a, MAX(b,c)))

#define MINMOD(x,y)  \
     (0.5*(copysign(1.0,x) + copysign(1.0,y)) * MIN(fabs(x), fabs(y)))

#define MINMOD4(w,x,y,z) \
     (0.125*( copysign(1.0,w)+copysign(1.0,x) )*fabs( (copysign(1.0,w)+copysign(1.0,y)) * (copysign(1.0,w)+copysign(1.0,z)) )*MIN4(fabs(w), fabs(x), fabs(y), fabs(z)))

     
static inline CCTK_REAL MP5(const CCTK_REAL am2, 
                            const CCTK_REAL am1, 
                            const CCTK_REAL a, 
                            const CCTK_REAL ap1, 
                            const CCTK_REAL ap2,
                            const CCTK_REAL anorm,
                            const CCTK_REAL mp5_eps,
                            const CCTK_REAL mp5_alpha
                            )
{
   const CCTK_REAL vl = (2.0*am2 - 13.0*am1 + 47.0*a + 27.0*ap1 - 3.0*ap2) * 1.0/60.0;
   const CCTK_REAL vmp = a + MINMOD( ap1-a, mp5_alpha*(a-am1) );
   if ((vl-a)*(vl-vmp) <= mp5_eps*anorm)
      return vl;
   else {
      const CCTK_REAL djm1 = am2 -2.0*am1 + a;
      const CCTK_REAL dj   = am1 -2.0*a + ap1;
      const CCTK_REAL djp1 = a -2.0*ap1 + ap2;
      const CCTK_REAL dm4jph = MINMOD4(4.0*dj-djp1, 4.0*djp1-dj, dj, djp1);
      const CCTK_REAL dm4jmh = MINMOD4(4.0*dj-djm1, 4.0*djm1-dj, dj, djm1);
      const CCTK_REAL vul = a + mp5_alpha*(a-am1);
      const CCTK_REAL vav = 0.5*(a+ap1);
      const CCTK_REAL vmd = vav - 0.5*dm4jph;
      const CCTK_REAL vlc = a + 0.5*(a-am1) + 4.0/3.0*dm4jmh;
      const CCTK_REAL vmin = MAX(MIN3(a,ap1,vmd), MIN3(a,vul,vlc));
      const CCTK_REAL vmax = MIN(MAX3(a,ap1,vmd), MAX3(a,vul,vlc));
      return vl + MINMOD(vmin-vl, vmax-vl);
   }
   return 0;
}


/**
   MP5 reconstruction operator.
*/
template <bool do_MP5_adaptive_epsilon>
template <int dir>
inline void
GRHydro_MP5Reconstruct1d_cxx<do_MP5_adaptive_epsilon>::
                   apply(const int nx,
                         const CCTK_REAL* const restrict a,
                         CCTK_REAL* const restrict aminus,
                         CCTK_REAL* const restrict aplus,
                         const cGH* const cctkGH,
                         const int j, const int k
                        )
{
#define A(i_) (a[ijk[i_]])
#define Aplus(i_) (aplus[ijk[i_]])
#define Aminus(i_) (aminus[ijk[i_]])

   DECLARE_CCTK_PARAMETERS;
   
   for (int i=GRHydro_stencil-1; i < nx-GRHydro_stencil+1; ++i)
   {
      const int ijk[5] = {
                            dir ==0 ? CCTK_GFINDEX3D(cctkGH, i-2, j, k) : dir ==1 ? CCTK_GFINDEX3D(cctkGH, j, i-2, k) : CCTK_GFINDEX3D(cctkGH, j, k, i-2), 
                            dir ==0 ? CCTK_GFINDEX3D(cctkGH, i-1, j, k) : dir ==1 ? CCTK_GFINDEX3D(cctkGH, j, i-1, k) : CCTK_GFINDEX3D(cctkGH, j, k, i-1),
                            dir ==0 ? CCTK_GFINDEX3D(cctkGH, i  , j, k) : dir ==1 ? CCTK_GFINDEX3D(cctkGH, j, i  , k) : CCTK_GFINDEX3D(cctkGH, j, k, i  ),
                            dir ==0 ? CCTK_GFINDEX3D(cctkGH, i+1, j, k) : dir ==1 ? CCTK_GFINDEX3D(cctkGH, j, i+1, k) : CCTK_GFINDEX3D(cctkGH, j, k, i+1),
                            dir ==0 ? CCTK_GFINDEX3D(cctkGH, i+2, j, k) : dir ==1 ? CCTK_GFINDEX3D(cctkGH, j, i+2, k) : CCTK_GFINDEX3D(cctkGH, j, k, i+2)
                         };
                       
      if (!do_MP5_adaptive_epsilon) {
                         
         Aplus(2)  = MP5(A(0), A(1), A(2), A(3), A(4), 1.0, mp5_eps, mp5_alpha);
         Aminus(2) = MP5(A(4), A(3), A(2), A(1), A(0), 1.0, mp5_eps, mp5_alpha);
   
      } else {
         
         const CCTK_REAL anorm = sqrt(SQR(A(0)) + SQR(A(1)) + SQR(A(2)) + SQR(A(3)) + SQR(A(4)));
         
         Aplus(2)  = MP5(A(0), A(1), A(2), A(3), A(4), anorm, mp5_eps, mp5_alpha);
         Aminus(2) = MP5(A(4), A(3), A(2), A(1), A(0), anorm, mp5_eps, mp5_alpha);
         
      }

   }
}

// instantiate all copies we need, this way different operators can be compiled
// in parallel. This must match the select routine in GRHydro_Reconstruct.cc
template class GRHydro_MP5Reconstruct1d_cxx<false>;
template class GRHydro_MP5Reconstruct1d_cxx<true>;

INSTANTIATE_RECONSTRUCTION_OPERATOR(GRHydro_MP5Reconstruct1d_cxx<false>)
INSTANTIATE_RECONSTRUCTION_OPERATOR(GRHydro_MP5Reconstruct1d_cxx<true >)
