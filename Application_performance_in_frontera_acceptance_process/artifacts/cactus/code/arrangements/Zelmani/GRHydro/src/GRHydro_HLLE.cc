#include <cmath>
#include <vector>
#include <iostream>
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "SpaceMask.h"

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
#define SQR(x) ((x)*(x))

// build up multiple element MIN/MAX from pairwise recursive comparison
// these are ~10 times faster than the MIN(MIN(MIN(... approach presumably
// since they allow for out-of-order evalutation of intermediate results
#define MIN3(a,b,c) MIN(MIN(a,b),c)
#define MIN4(a,b,c,d) MIN(MIN(a,b),MIN(c,d))
#define MIN5(a,b,c,d,e) MIN(MIN4(a,b,c,d),e)
#define MIN6(a,b,c,d,e,f) MIN(MIN4(a,b,c,d),MIN(e,f))
#define MIN7(a,b,c,d,e,f,g) MIN(MIN4(a,b,c,d),MIN3(e,f,g))
#define MIN8(a,b,c,d,e,f,g,h) MIN(MIN4(a,b,c,d),MIN4(e,f,g,h))
#define MIN9(a,b,c,d,e,f,g,h,i) MIN(MIN8(a,b,c,d,e,f,g,h),i)
#define MIN10(a,b,c,d,e,f,g,h,i,j) MIN(MIN8(a,b,c,d,e,f,g,h),MIN(i,j))
#define MIN11(a,b,c,d,e,f,g,h,i,j,k) MIN(MIN8(a,b,c,d,e,f,g,h),MIN3(i,j,k))

#define MAX3(a,b,c) MAX(MAX(a,b),c)
#define MAX4(a,b,c,d) MAX(MAX(a,b),MAX(c,d))
#define MAX5(a,b,c,d,e) MAX(MAX4(a,b,c,d),e)
#define MAX6(a,b,c,d,e,f) MAX(MAX4(a,b,c,d),MAX(e,f))
#define MAX7(a,b,c,d,e,f,g) MAX(MAX4(a,b,c,d),MAX3(e,f,g))
#define MAX8(a,b,c,d,e,f,g,h) MAX(MAX4(a,b,c,d),MAX4(e,f,g,h))
#define MAX9(a,b,c,d,e,f,g,h,i) MAX(MAX8(a,b,c,d,e,f,g,h),i)
#define MAX10(a,b,c,d,e,f,g,h,i,j) MAX(MAX8(a,b,c,d,e,f,g,h),MAX(i,j))
#define MAX11(a,b,c,d,e,f,g,h,i,j,k) MAX(MAX8(a,b,c,d,e,f,g,h),MAX3(i,j,k))

using namespace std;

// some prototypes
extern "C"
CCTK_INT GRHydro_UseGeneralCoordinates(const cGH * cctkGH);

static inline CCTK_REAL SpatialDeterminantC(const CCTK_REAL gxx, const CCTK_REAL gxy,
                                            const CCTK_REAL gxz, const CCTK_REAL gyy,
                                            const CCTK_REAL gyz, const CCTK_REAL gzz);
static inline void UpperMetricC(CCTK_REAL& restrict uxx,
                                CCTK_REAL& restrict uxy,
                                CCTK_REAL& restrict uxz,
                                CCTK_REAL& restrict uyy,
                                CCTK_REAL& restrict uyz,
                                CCTK_REAL& restrict uzz,
                                const CCTK_REAL det,
                                const CCTK_REAL gxx,
                                const CCTK_REAL gxy,
                                const CCTK_REAL gxz,
                                const CCTK_REAL gyy,
                                const CCTK_REAL gyz,
                                const CCTK_REAL gzz);

static inline CCTK_REAL max10(const CCTK_REAL * const restrict a)
{
  return MAX11(0., a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8], a[9]);
}

static inline CCTK_REAL max9(const CCTK_REAL a,
                          const CCTK_REAL b,
                          const CCTK_REAL c,
                          const CCTK_REAL d,
                          const CCTK_REAL e,
                          const CCTK_REAL f,
                          const CCTK_REAL g,
                          const CCTK_REAL h,
                          const CCTK_REAL i)

{
  return MAX9(a, b, c, d, e, f, g, h, i);
}

static inline CCTK_REAL min10(const CCTK_REAL * const restrict a)
{
  return MIN11(0., a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8], a[9]);
}

static inline void num_x_flux_cc(const CCTK_REAL dens,
                                 const CCTK_REAL sx,
                                 const CCTK_REAL sy,
                                 const CCTK_REAL sz,
                                 const CCTK_REAL tau,
                                 CCTK_REAL& densflux,
                                 CCTK_REAL& sxflux,
                                 CCTK_REAL& syflux,
                                 CCTK_REAL& szflux,
                                 CCTK_REAL& tauflux,
                                 const CCTK_REAL vel,
                                 const CCTK_REAL press,
                                 const CCTK_REAL sdet,
                                 const CCTK_REAL alp,
                                 const CCTK_REAL beta)
{

  const CCTK_REAL velmbetainvalp = vel - beta / alp;
  const CCTK_REAL sqrtdetpress = sdet*press;

  densflux = dens * velmbetainvalp;
  sxflux   = sx * velmbetainvalp + sqrtdetpress;
  syflux   = sy * velmbetainvalp;
  szflux   = sz * velmbetainvalp;
  tauflux  = tau * velmbetainvalp + sqrtdetpress*vel;

  return;
}

static inline void num_x_fluxM_cc(const CCTK_REAL dens,
                                 const CCTK_REAL sx,
                                 const CCTK_REAL sy,
                                 const CCTK_REAL sz,
                                 const CCTK_REAL tau,
                                 const CCTK_REAL Bx,
                                 const CCTK_REAL By,
                                 const CCTK_REAL Bz,
                                 CCTK_REAL& densflux,
                                 CCTK_REAL& sxflux,
                                 CCTK_REAL& syflux,
                                 CCTK_REAL& szflux,
                                 CCTK_REAL& tauflux,
                                 CCTK_REAL& Bxflux,
                                 CCTK_REAL& Byflux,
                                 CCTK_REAL& Bzflux,
                                 const CCTK_REAL velx,
                                 const CCTK_REAL vely,
                                 const CCTK_REAL velz,
                                 const CCTK_REAL press,
                                 const CCTK_REAL sdet,
                                 const CCTK_REAL alp,
                                 const CCTK_REAL betax,
                                 const CCTK_REAL betay,
                                 const CCTK_REAL betaz,
                                 const CCTK_REAL bsubx,
                                 const CCTK_REAL bsuby,
                                 const CCTK_REAL bsubz,
                                 const CCTK_REAL b2,
                                 const CCTK_REAL ab0,
                                 const CCTK_REAL w)
{

  const CCTK_REAL velmbetainvalpx = velx - betax / alp;
  const CCTK_REAL velmbetainvalpy = vely - betay / alp;
  const CCTK_REAL velmbetainvalpz = velz - betaz / alp;
  const CCTK_REAL pressstar = press + 0.5*b2;
  const CCTK_REAL sqrtdetpressstar = sdet*pressstar;

  densflux = dens * velmbetainvalpx;
  sxflux   = sx * velmbetainvalpx + sqrtdetpressstar - bsubx*Bx/w;
  syflux   = sy * velmbetainvalpx - bsuby*Bx/w;
  szflux   = sz * velmbetainvalpx - bsubz*Bx/w;
  tauflux  = tau * velmbetainvalpx + sqrtdetpressstar*velx - ab0*Bx/w;
  Bxflux = 0.0;
  Byflux = By * velmbetainvalpx - Bx*velmbetainvalpy;
  Bzflux = Bz * velmbetainvalpx - Bx*velmbetainvalpz;
  return;
}

static inline void eigenvalues_cc(const CCTK_REAL rho,
                                  const CCTK_REAL velx,
                                  const CCTK_REAL vely,
                                  const CCTK_REAL velz,
                                  const CCTK_REAL eps,
                                  const CCTK_REAL press,
                                  const CCTK_REAL cs2,
                                  const CCTK_REAL w,
                                  const CCTK_REAL gxx,
                                  const CCTK_REAL gxy,
                                  const CCTK_REAL gxz,
                                  const CCTK_REAL gyy,
                                  const CCTK_REAL gyz,
                                  const CCTK_REAL gzz,
                                  const CCTK_REAL alp,
                                  const CCTK_REAL beta,
                                  const CCTK_REAL u,
                                  CCTK_REAL * const restrict lam)
{

  const CCTK_REAL vlowx = gxx*velx + gxy*vely + gxz*velz;
  const CCTK_REAL vlowy = gxy*velx + gyy*vely + gyz*velz;
  const CCTK_REAL vlowz = gxz*velx + gyz*vely + gzz*velz;
  const CCTK_REAL v2 = vlowx*velx + vlowy*vely + vlowz*velz;

  const CCTK_REAL boa = beta/alp;
  lam[1] = velx - boa;
  lam[2] = lam[1];
  lam[3] = lam[1];
  const CCTK_REAL lam_tmp1 = 1.0e0/(1.0-v2*cs2);
  const CCTK_REAL lam_tmp2 = sqrt(cs2*(1.0-v2)*
                               (u*(1.0-v2*cs2) - velx*velx*(1.0-cs2)));
  const CCTK_REAL lam_tmp3 = velx*(1.0-cs2);
  lam[0] = (lam_tmp3 - lam_tmp2)*lam_tmp1 - boa;
  lam[4] = (lam_tmp3 + lam_tmp2)*lam_tmp1 - boa;

  return;
}

static inline void eigenvaluesM_cc(const CCTK_REAL rho,
                                  const CCTK_REAL velx,
                                  const CCTK_REAL vely,
                                  const CCTK_REAL velz,
                                  const CCTK_REAL eps,
                                  const CCTK_REAL press,
                                  const CCTK_REAL Bx,
                                  const CCTK_REAL By,
                                  const CCTK_REAL Bz,
                                  const CCTK_REAL cs2,
                                  const CCTK_REAL w,
                                  const CCTK_REAL gxx,
                                  const CCTK_REAL gxy,
                                  const CCTK_REAL gxz,
                                  const CCTK_REAL gyy,
                                  const CCTK_REAL gyz,
                                  const CCTK_REAL gzz,
                                  const CCTK_REAL alp,
                                  const CCTK_REAL beta,
                                  const CCTK_REAL u,
                                  CCTK_REAL * const restrict lam)
{

  const CCTK_REAL vlowx = gxx*velx + gxy*vely + gxz*velz;
  const CCTK_REAL vlowy = gxy*velx + gyy*vely + gyz*velz;
  const CCTK_REAL vlowz = gxz*velx + gyz*vely + gzz*velz;
  const CCTK_REAL v2 = vlowx*velx + vlowy*vely + vlowz*velz;

  const CCTK_REAL boa = beta/alp;
  const CCTK_REAL Bxlow = gxx*Bx + gxy*By + gxz*Bz;
  const CCTK_REAL Bylow = gxy*Bx + gyy*By + gyz*Bz;
  const CCTK_REAL Bzlow = gxz*Bx + gyz*By + gzz*Bz;
  const CCTK_REAL B2 =Bxlow*Bx+Bylow*By+Bzlow*Bz;
  const CCTK_REAL Bdotv =Bxlow*velx+Bylow*vely+Bzlow*velz;
  const CCTK_REAL Bdotv2 =Bdotv*Bdotv;
  const CCTK_REAL w2 = w*w;
  const CCTK_REAL b2 =B2/w2+Bdotv2;
  const CCTK_REAL rhos = rho*(1.0+eps)+press+b2;
  const CCTK_REAL va2 = b2/rhos;
  const CCTK_REAL u2 = va2+cs2*(1.0-va2);

  lam[1] = velx - boa;
  lam[2] = lam[1];
  lam[3] = lam[1];

  const CCTK_REAL lam_tmp1 = (velx*(1.0-u2) - sqrt(u2*(1.0-v2)*
                               (u*(1.0-v2*u2) - velx*velx*(1.0-u2))))/(1.0-v2*u2);
  const CCTK_REAL lam_tmp2 = (velx*(1.0-u2) + sqrt(u2*(1.0-v2)*
                               (u*(1.0-v2*u2) - velx*velx*(1.0-u2))))/(1.0-v2*u2);

  lam[0] = lam_tmp1 - boa;
  lam[4] = lam_tmp2 - boa;

  return;
}

// helpers for H viscosity
static inline CCTK_REAL etaX(const cGH *cctkGH,
                          const int i, const int j, const int k,
                          const CCTK_REAL * const restrict vel,
                          const CCTK_REAL * const restrict cs) {

  const int vidxr = CCTK_VECTGFINDEX3D(cctkGH,i+1,j,k,0);
  const int vidx = CCTK_VECTGFINDEX3D(cctkGH,i,j,k,0);
  const int idxr = CCTK_GFINDEX3D(cctkGH,i+1,j,k);
  const int idx = CCTK_GFINDEX3D(cctkGH,i,j,k);

  return 0.5 * (fabs(vel[vidxr]-vel[vidx]) + fabs(cs[idxr]-cs[idx]));
}

static inline CCTK_REAL etaY(const cGH *cctkGH,
                          const int i, const int j, const int k,
                          const CCTK_REAL * const restrict vel,
                          const CCTK_REAL * const restrict cs) {

  const int vidxr = CCTK_VECTGFINDEX3D(cctkGH,i,j+1,k,1);
  const int vidx = CCTK_VECTGFINDEX3D(cctkGH,i,j,k,1);
  const int idxr = CCTK_GFINDEX3D(cctkGH,i,j+1,k);
  const int idx = CCTK_GFINDEX3D(cctkGH,i,j,k);

  return 0.5 * (fabs(vel[vidxr]-vel[vidx]) + fabs(cs[idxr]-cs[idx]));
}

static inline CCTK_REAL etaZ(const cGH *cctkGH,
                          const int i, const int j, const int k,
                          const CCTK_REAL * const restrict vel,
                          const CCTK_REAL * const restrict cs) {

  const int vidxr = CCTK_VECTGFINDEX3D(cctkGH,i,j,k+1,2);
  const int vidx = CCTK_VECTGFINDEX3D(cctkGH,i,j,k,2);
  const int idxr = CCTK_GFINDEX3D(cctkGH,i,j,k+1);
  const int idx = CCTK_GFINDEX3D(cctkGH,i,j,k);

  return 0.5 * (fabs(vel[vidxr]-vel[vidx]) + fabs(cs[idxr]-cs[idx]));
}

static inline void vlow_blow(const CCTK_REAL gxx,const CCTK_REAL gxy,const CCTK_REAL gxz,const CCTK_REAL gyy,const CCTK_REAL gyz,const CCTK_REAL gzz,
     const CCTK_REAL velx,const CCTK_REAL vely,const CCTK_REAL velz,const CCTK_REAL Bvecx,const CCTK_REAL Bvecy,const CCTK_REAL Bvecz,
     CCTK_REAL& ab0,CCTK_REAL& b2,CCTK_REAL& w,CCTK_REAL& bxlow,CCTK_REAL& bylow,CCTK_REAL& bzlow)
{

  // Calculates v_i (see Anton Eq. 5) and B_i (Bvecxlow)- undensitized!
  // calculates B^i v_i [Anton eq. 44] and b^2 [LHS of Anton eq. 46]
  // Calculates w (Lorentz factor) as (1-v^i v_i)^{-1/2}
  // Calculates b_i (bxlow)

  // vel_i  = g_ij v^j
  // B_i = g_ij B^i
  const CCTK_REAL velxlow = gxx*velx + gxy*vely + gxz*velz;
  const CCTK_REAL velylow = gxy*velx + gyy*vely + gyz*velz;
  const CCTK_REAL velzlow = gxz*velx + gyz*vely + gzz*velz;
  const CCTK_REAL Bvecxlow = gxx*Bvecx + gxy*Bvecy + gxz*Bvecz;
  const CCTK_REAL Bvecylow = gxy*Bvecx + gyy*Bvecy + gyz*Bvecz;
  const CCTK_REAL Bveczlow = gxz*Bvecx + gyz*Bvecy + gzz*Bvecz;

  // B^i v_i (= b^0/u^0)
  const CCTK_REAL Bdotv = velxlow*Bvecx+velylow*Bvecy+velzlow*Bvecz;

  // v^2 = v_i v^i; w=(1-v^2)^{-1/2}
  const CCTK_REAL v2 = velxlow*velx + velylow*vely + velzlow*velz;
  w = 1./sqrt(1.-v2);

  // b^2 = B^i B_i / w^2 + (b^0/u^0)^2
  b2=(Bvecx*Bvecxlow+Bvecy*Bvecylow+Bvecz*Bveczlow)/SQR(w)+SQR(Bdotv);

  // b_i = B_i/w +w*(B dot v)*v_i
  bxlow = Bvecxlow/w+w*Bdotv*velxlow;
  bylow = Bvecylow/w+w*Bdotv*velylow;
  bzlow = Bveczlow/w+w*Bdotv*velzlow;

  // alpha b_0 = w_lorentz B^i vel_i
  ab0 = w * Bdotv;
}

template<const int fdir, const bool do_ye, const bool do_hvisc, const bool do_tracers, const bool do_mhd, const bool do_clean_divergence>
void GRHydro_HLLE_CC_LL(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // save memory when multipatch is not used
  const CCTK_REAL * restrict g11, * restrict g12, * restrict g13;
  const CCTK_REAL * restrict g22, * restrict g23, * restrict g33;
  const CCTK_REAL * restrict beta1, * restrict beta2, * restrict beta3;
  const CCTK_REAL * restrict vup;
  const int offsetx = fdir==1, offsety = fdir==2, offsetz = fdir==3;
  assert(offsetx+offsety+offsetz == 1);

  if(GRHydro_disable_hydro_update) return;

  if(GRHydro_UseGeneralCoordinates(cctkGH)) {
    g11   = gaa;
    g12   = gab;
    g13   = gac;
    g22   = gbb;
    g23   = gbc;
    g33   = gcc;
    beta1 = betaa;
    beta2 = betab;
    beta3 = betac;
    vup   = lvel;
   } else {
    g11   = gxx;
    g12   = gxy;
    g13   = gxz;
    g22   = gyy;
    g23   = gyz;
    g33   = gzz;
    beta1 = betax;
    beta2 = betay;
    beta3 = betaz;
    vup   = vel;
  }

  int type_bits, trivial;
  if(fdir==1) {
    type_bits = SpaceMask_GetTypeBits("Hydro_RiemannProblemX");
    trivial = SpaceMask_GetStateBits("Hydro_RiemannProblemX",
                                            "trivial");
  } else if (fdir==2) {
    type_bits = SpaceMask_GetTypeBits("Hydro_RiemannProblemY");
    trivial = SpaceMask_GetStateBits("Hydro_RiemannProblemY",
                                         "trivial");
  } else if (fdir==3) {
    type_bits = SpaceMask_GetTypeBits("Hydro_RiemannProblemZ");
    trivial = SpaceMask_GetStateBits("Hydro_RiemannProblemZ",
                                         "trivial");
  } else {
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                "Flux direction %d not 1,2,3.", fdir);
  }

#pragma omp parallel for
  for(int k = GRHydro_stencil-1; k < cctk_lsh[2]-GRHydro_stencil+ transport_constraints*(1-offsetz); ++k)
    for(int j = GRHydro_stencil-1; j < cctk_lsh[1]-GRHydro_stencil+ transport_constraints*(1-offsety); ++j)
      for(int i = GRHydro_stencil-1; i < cctk_lsh[0]-GRHydro_stencil+ transport_constraints*(1-offsetx); ++i) {

        const int idx = CCTK_GFINDEX3D(cctkGH,i,j,k);
        const int idxr =
          CCTK_GFINDEX3D(cctkGH,i+offsetx,j+offsety,k+offsetz);

        // C++98 does not allow variable sized arrays (C99 does)
        // C++98 does however allow the use of a const-variable as the aray
        // size as long as its value is known at compile time
        // Unfortunately the Intel compiler (v11 at least) produces lots of
        // warnings about out of bounds array accesses unless the arrays are
        // declared with the maximum possible size.
        const int activesize = do_mhd ? 8 : 5;

        CCTK_REAL consp[8];
        CCTK_REAL consm_i[8];
        CCTK_REAL fplus[8];
        CCTK_REAL fminus[8];
        CCTK_REAL f1[8];
        CCTK_REAL qdiff[8];

        CCTK_REAL psidcp, psidcm_i;
        CCTK_REAL psidcfplus, psidcfminus;
        CCTK_REAL psidcf1;
        
        CCTK_REAL lamplusminus[10];
        CCTK_REAL * const lamminus = lamplusminus;
        CCTK_REAL * const lamplus = lamplusminus+5;
        CCTK_REAL avg_beta;

        consp[0] = densplus[idx];
        consp[1] = sxplus[idx];
        consp[2] = syplus[idx];
        consp[3] = szplus[idx];
        consp[4] = tauplus[idx];
        if (do_mhd) {
          consp[5] = Bconsxplus[idx];
          consp[6] = Bconsyplus[idx];
          consp[7] = Bconszplus[idx];
        }
        if (do_clean_divergence)
          psidcp = psidcplus[idx];
        consm_i[0] = densminus[idxr];
        consm_i[1] = sxminus[idxr];
        consm_i[2] = syminus[idxr];
        consm_i[3] = szminus[idxr];
        consm_i[4] = tauminus[idxr];
        if (do_mhd) {
          consm_i[5] = Bconsxminus[idxr];
          consm_i[6] = Bconsyminus[idxr];
          consm_i[7] = Bconszminus[idxr];
        }
        if (do_clean_divergence)
          psidcm_i = psidcminus[idxr];

        if(fdir==1) 
          avg_beta = 0.5 * (beta1[idxr]+beta1[idx]);
        else if(fdir==2)
          avg_beta = 0.5 * (beta2[idxr]+beta2[idx]);
        else 
          avg_beta = 0.5 * (beta3[idxr]+beta3[idx]);

        const CCTK_REAL avg_betax = 0.5 * (beta1[idxr]+beta1[idx]);
        const CCTK_REAL avg_betay = 0.5 * (beta2[idxr]+beta2[idx]);
        const CCTK_REAL avg_betaz = 0.5 * (beta3[idxr]+beta3[idx]);
        const CCTK_REAL avg_alp = 0.5 * (alp[idxr]+alp[idx]);

        const CCTK_REAL gxxh = 0.5 * (g11[idx] + g11[idxr]);
        const CCTK_REAL gxyh = 0.5 * (g12[idx] + g12[idxr]);
        const CCTK_REAL gxzh = 0.5 * (g13[idx] + g13[idxr]);
        const CCTK_REAL gyyh = 0.5 * (g22[idx] + g22[idxr]);
        const CCTK_REAL gyzh = 0.5 * (g23[idx] + g23[idxr]);
        const CCTK_REAL gzzh = 0.5 * (g33[idx] + g33[idxr]);
  
        const CCTK_REAL avg_detr =
          SpatialDeterminantC(gxxh,gxyh,gxzh,gyyh,gyzh,gzzh);
        const CCTK_REAL savg_detr = sqrt(avg_detr);

	CCTK_REAL uxxh,uxyh,uxzh,uyyh,uyzh,uzzh;
	UpperMetricC(uxxh,uxyh,uxzh,uyyh,uyzh,uzzh,
                        avg_detr,gxxh,gxyh,gxzh,
	                gyyh, gyzh, gzzh);

        // these need to be initialized otherwise the compiler complains even
        // if they are only ever used when do_mhd is true
        CCTK_REAL ab0m = -42.0;
        CCTK_REAL ab0p = -42.0;
        CCTK_REAL b2m = -42.0;
        CCTK_REAL b2p = -42.0;
        CCTK_REAL wm = -42.0;
        CCTK_REAL wp = -42.0;
        CCTK_REAL blowxm = -42.0;
        CCTK_REAL blowym = -42.0;
        CCTK_REAL blowzm = -42.0;
        CCTK_REAL blowxp = -42.0;
        CCTK_REAL blowyp = -42.0;
        CCTK_REAL blowzp = -42.0;

        if (do_mhd) {
          // TODO: we do actually store w_lorentz on the grid and could use it
          vlow_blow(gxxh,gxyh,gxzh,gyyh,gyzh,gzzh,
               velxplus[idx],velyplus[idx],velzplus[idx],Bvecxplus[idx],Bvecyplus[idx],Bveczplus[idx],
               ab0p,b2p,wp,blowxp,blowyp,blowzp);

          vlow_blow(gxxh,gxyh,gxzh,gyyh,gyzh,gzzh,
               velxminus[idxr],velyminus[idxr],velzminus[idxr],Bvecxminus[idxr],Bvecyminus[idxr],Bveczminus[idxr],
               ab0m,b2m,wm,blowxm,blowym,blowzm);
        }

        // If the Riemann problem is trivial, just calculate the fluxes
        // from the left state and skip to the next cell
        if ( SpaceMask_CheckStateBits(space_mask, idx, type_bits, trivial) ) {

          if (do_mhd) {

            if (fdir==1) {
              num_x_fluxM_cc(consp[0],consp[1],consp[2],consp[3],consp[4],consp[5],consp[6],consp[7],
                          f1[0],f1[1],f1[2],f1[3],f1[4],f1[5],f1[6],f1[7],
                          velxplus[idx],velyplus[idx],velzplus[idx],
                          pressplus[idx],savg_detr,
                          avg_alp,avg_betax,avg_betay,avg_betaz,blowxp,blowyp,blowzp,b2p,ab0p,wp);
              if (do_clean_divergence) {
                f1[5] += savg_detr*uxxh*psidcp - consp[5]*avg_betax/avg_alp;
                f1[6] += savg_detr*uxyh*psidcp - consp[5]*avg_betay/avg_alp;
                f1[7] += savg_detr*uxzh*psidcp - consp[5]*avg_betaz/avg_alp;
                psidcf1 = consp[5]/savg_detr - psidcp*avg_betax/avg_alp;
              }
            } else if(fdir==2) {
              num_x_fluxM_cc(consp[0],consp[2],consp[3],consp[1],consp[4],consp[6],consp[7],consp[5],
                          f1[0],f1[2],f1[3],f1[1],f1[4],f1[6],f1[7],f1[5],
                          velyplus[idx],velzplus[idx],velxplus[idx],
                          pressplus[idx],savg_detr,
                          avg_alp,avg_betay,avg_betaz,avg_betax,blowyp,blowzp,blowxp,b2p,ab0p,wp);
              if (do_clean_divergence) {
                f1[5] += savg_detr*uxyh*psidcp - consp[6]*avg_betax/avg_alp;
                f1[6] += savg_detr*uyyh*psidcp - consp[6]*avg_betay/avg_alp;
                f1[7] += savg_detr*uyzh*psidcp - consp[6]*avg_betaz/avg_alp;
                psidcf1 = consp[6]/savg_detr - psidcp*avg_betay/avg_alp;
              }
            } else {
              num_x_fluxM_cc(consp[0],consp[3],consp[1],consp[2],consp[4],consp[7],consp[5],consp[6],
                          f1[0],f1[3],f1[1],f1[2],f1[4],f1[7],f1[5],f1[6],
                          velzplus[idx],velxplus[idx],velyplus[idx],
                          pressplus[idx],savg_detr,
                          avg_alp,avg_betaz,avg_betax,avg_betay,blowzp,blowxp,blowyp,b2p,ab0p,wp);
              if (do_clean_divergence) {
                f1[5] += savg_detr*uxzh*psidcp - consp[7]*avg_betax/avg_alp;
                f1[6] += savg_detr*uyzh*psidcp - consp[7]*avg_betay/avg_alp;
                f1[7] += savg_detr*uzzh*psidcp - consp[7]*avg_betaz/avg_alp;
                psidcf1 = consp[7]/savg_detr - psidcp*avg_betaz/avg_alp;
              }
            }
          } else {
            if(fdir==1)
              num_x_flux_cc(consp[0],consp[1],consp[2],consp[3],consp[4],
                            f1[0],f1[1],f1[2],f1[3],f1[4],
                            velxplus[idx],pressplus[idx],savg_detr,
                            avg_alp,avg_beta);
            else if (fdir==2)
              num_x_flux_cc(consp[0],consp[2],consp[3],consp[1],consp[4],
                            f1[0],f1[2],f1[3],f1[1],f1[4],
                            velyplus[idx],pressplus[idx],savg_detr,
                            avg_alp,avg_beta);
            else
              num_x_flux_cc(consp[0],consp[3],consp[1],consp[2],consp[4],
                            f1[0],f1[3],f1[1],f1[2],f1[4],
                            velzplus[idx],pressplus[idx],savg_detr,
                            avg_alp,avg_beta);
          }

        } else {

          // normal case -- full HLLE Riemann solution

          CCTK_REAL usendh;
          if(fdir==1)
            usendh = uxxh;
          else if(fdir==2)
            usendh = uyyh;
          else
            usendh = uzzh;

          // calculate the jumps in the conserved variables
          for(int m=0;m<activesize;++m)
            qdiff[m] = consm_i[m] - consp[m];

          if (do_mhd) {

            if (fdir==1) {
              eigenvaluesM_cc(rhominus[idxr],
                           velxminus[idxr],
                           velyminus[idxr],
                           velzminus[idxr],
                           epsminus[idxr],
                           pressminus[idxr],
                           Bvecxminus[idxr],
                           Bvecyminus[idxr],
                           Bveczminus[idxr],
                           cs2minus[idxr],
                           w_lorentzminus[idxr],
                           gxxh,gxyh,gxzh,gyyh,gyzh,gzzh,
                           avg_alp,avg_beta,usendh,
                           lamminus);
              eigenvaluesM_cc(rhoplus[idx],
                           velxplus[idx],
                           velyplus[idx],
                           velzplus[idx],
                           epsplus[idx],
                           pressplus[idx],
                           Bvecxplus[idx],
                           Bvecyplus[idx],
                           Bveczplus[idx],
                           cs2plus[idx],
                           w_lorentzplus[idx],
                           gxxh,gxyh,gxzh,gyyh,gyzh,gzzh,
                           avg_alp,avg_beta,usendh,
                           lamplus);

              num_x_fluxM_cc(consp[0],consp[1],consp[2],consp[3],consp[4],consp[5],consp[6],consp[7],
                          fplus[0],fplus[1],fplus[2],fplus[3],fplus[4],fplus[5],fplus[6],fplus[7],
                          velxplus[idx],velyplus[idx],velzplus[idx],
                          pressplus[idx],savg_detr,
                          avg_alp,avg_betax,avg_betay,avg_betaz,blowxp,blowyp,blowzp,b2p,ab0p,wp);
              num_x_fluxM_cc(consm_i[0],consm_i[1],consm_i[2],consm_i[3],consm_i[4],consm_i[5],consm_i[6],consm_i[7],
                          fminus[0],fminus[1],fminus[2],fminus[3],fminus[4],fminus[5],fminus[6],fminus[7],
                          velxminus[idxr],velyminus[idxr],velzminus[idxr],
                          pressminus[idxr],savg_detr,
                          avg_alp,avg_betax,avg_betay,avg_betaz,blowxm,blowym,blowzm,b2m,ab0m,wm);
              if (do_clean_divergence) {
                fminus[5] += savg_detr*uxxh*psidcm_i - consm_i[5]*avg_betax/avg_alp;
                fminus[6] += savg_detr*uxyh*psidcm_i - consm_i[5]*avg_betay/avg_alp;
                fminus[7] += savg_detr*uxzh*psidcm_i - consm_i[5]*avg_betaz/avg_alp;
                fplus[5] += savg_detr*uxxh*psidcp - consp[5]*avg_betax/avg_alp;
                fplus[6] += savg_detr*uxyh*psidcp - consp[5]*avg_betay/avg_alp;
                fplus[7] += savg_detr*uxzh*psidcp - consp[5]*avg_betaz/avg_alp;
                psidcfminus = consm_i[5]/savg_detr - psidcm_i*avg_betax/avg_alp;
                psidcfplus = consp[5]/savg_detr - psidcp*avg_betax/avg_alp;
              }
            } else if(fdir==2) {
              eigenvaluesM_cc(rhominus[idxr],
                           velyminus[idxr],
                           velzminus[idxr],
                           velxminus[idxr],
                           epsminus[idxr],
                           pressminus[idxr],
                           Bvecyminus[idxr],
                           Bveczminus[idxr],
                           Bvecxminus[idxr],
                           cs2minus[idxr],
                           w_lorentzminus[idxr],
                           gyyh,gyzh,gxyh,gzzh,gxzh,gxxh,
                           avg_alp,avg_beta,usendh,
                           lamminus);
              eigenvaluesM_cc(rhoplus[idx],
                           velyplus[idx],
                           velzplus[idx],
                           velxplus[idx],
                           epsplus[idx],
                           pressplus[idx],
                           Bvecyplus[idx],
                           Bveczplus[idx],
                           Bvecxplus[idx],
                           cs2plus[idx],
                           w_lorentzplus[idx],
                           gyyh,gyzh,gxyh,gzzh,gxzh,gxxh,
                           avg_alp,avg_beta,usendh,
                           lamplus);
              num_x_fluxM_cc(consp[0],consp[2],consp[3],consp[1],consp[4],consp[6],consp[7],consp[5],
                          fplus[0],fplus[2],fplus[3],fplus[1],fplus[4],fplus[6],fplus[7],fplus[5],
                          velyplus[idx],velzplus[idx],velxplus[idx],
                          pressplus[idx],savg_detr,
                          avg_alp,avg_betay,avg_betaz,avg_betax,blowyp,blowzp,blowxp,b2p,ab0p,wp);
              num_x_fluxM_cc(consm_i[0],consm_i[2],consm_i[3],consm_i[1],consm_i[4],consm_i[6],consm_i[7],consm_i[5],
                          fminus[0],fminus[2],fminus[3],fminus[1],fminus[4],fminus[6],fminus[7],fminus[5],
                          velyminus[idxr],velzminus[idxr],velxminus[idxr],
                          pressminus[idxr],savg_detr,
                          avg_alp,avg_betay,avg_betaz,avg_betax,blowym,blowzm,blowxm,b2m,ab0m,wm);
              if (do_clean_divergence) {
                fminus[5] += savg_detr*uxyh*psidcm_i - consm_i[6]*avg_betax/avg_alp;
                fminus[6] += savg_detr*uyyh*psidcm_i - consm_i[6]*avg_betay/avg_alp;
                fminus[7] += savg_detr*uyzh*psidcm_i - consm_i[6]*avg_betaz/avg_alp;
                fplus[5] += savg_detr*uxyh*psidcp - consp[6]*avg_betax/avg_alp;
                fplus[6] += savg_detr*uyyh*psidcp - consp[6]*avg_betay/avg_alp;
                fplus[7] += savg_detr*uyzh*psidcp - consp[6]*avg_betaz/avg_alp;
                psidcfminus = consm_i[6]/savg_detr - psidcm_i*avg_betay/avg_alp;
                psidcfplus = consp[6]/savg_detr - psidcp*avg_betay/avg_alp;
              }
            } else {
              eigenvaluesM_cc(rhominus[idxr],
                           velzminus[idxr],
                           velxminus[idxr],
                           velyminus[idxr],
                           epsminus[idxr],
                           pressminus[idxr],
                           Bveczminus[idxr],
                           Bvecxminus[idxr],
                           Bvecyminus[idxr],
                           cs2minus[idxr],
                           w_lorentzminus[idxr],
                           gzzh,gxzh,gyzh,gxxh,gxyh,gyyh,
                           avg_alp,avg_beta,usendh,
                           lamminus);
              eigenvaluesM_cc(rhoplus[idx],
                           velzplus[idx],
                           velxplus[idx],
                           velyplus[idx],
                           epsplus[idx],
                           pressplus[idx],
                           Bveczplus[idx],
                           Bvecxplus[idx],
                           Bvecyplus[idx],
                           cs2plus[idx],
                           w_lorentzplus[idx],
                           gzzh,gxzh,gyzh,gxxh,gxyh,gyyh,
                           avg_alp,avg_beta,usendh,
                           lamplus);
              num_x_fluxM_cc(consp[0],consp[3],consp[1],consp[2],consp[4],consp[7],consp[5],consp[6],
                          fplus[0],fplus[3],fplus[1],fplus[2],fplus[4],fplus[7],fplus[5],fplus[6],
                          velzplus[idx],velxplus[idx],velyplus[idx],
                          pressplus[idx],savg_detr,
                          avg_alp,avg_betaz,avg_betax,avg_betay,blowzp,blowxp,blowyp,b2p,ab0p,wp);
              num_x_fluxM_cc(consm_i[0],consm_i[3],consm_i[1],consm_i[2],consm_i[4],consm_i[7],consm_i[5],consm_i[6],
                          fminus[0],fminus[3],fminus[1],fminus[2],fminus[4],fminus[7],fminus[5],fminus[6],
                          velzminus[idxr],velxminus[idxr],velyminus[idxr],
                          pressminus[idxr],savg_detr,
                          avg_alp,avg_betaz,avg_betax,avg_betay,blowzm,blowxm,blowym,b2m,ab0m,wm);
              if (do_clean_divergence) {
                fminus[5] += savg_detr*uxzh*psidcm_i - consm_i[7]*avg_betax/avg_alp;
                fminus[6] += savg_detr*uyzh*psidcm_i - consm_i[7]*avg_betay/avg_alp;
                fminus[7] += savg_detr*uzzh*psidcm_i - consm_i[7]*avg_betaz/avg_alp;
                fplus[5] += savg_detr*uxzh*psidcp - consp[7]*avg_betax/avg_alp;
                fplus[6] += savg_detr*uyzh*psidcp - consp[7]*avg_betay/avg_alp;
                fplus[7] += savg_detr*uzzh*psidcp - consp[7]*avg_betaz/avg_alp;
                psidcfminus = consm_i[7]/savg_detr - psidcm_i*avg_betaz/avg_alp;
                psidcfplus = consp[7]/savg_detr - psidcp*avg_betaz/avg_alp;
              }
            } 
          } else {

            if (fdir==1) {
              eigenvalues_cc(rhominus[idxr],
                           velxminus[idxr],
                           velyminus[idxr],
                           velzminus[idxr],
                           epsminus[idxr],
                           pressminus[idxr],
                           cs2minus[idxr],
                           w_lorentzminus[idxr],
                           gxxh,gxyh,gxzh,gyyh,gyzh,gzzh,
                           avg_alp,avg_beta,usendh,
                           lamminus);
              eigenvalues_cc(rhoplus[idx],
                           velxplus[idx],
                           velyplus[idx],
                           velzplus[idx],
                           epsplus[idx],
                           pressplus[idx],
                           cs2plus[idx],
                           w_lorentzplus[idx],
                           gxxh,gxyh,gxzh,gyyh,gyzh,gzzh,
                           avg_alp,avg_beta,usendh,
                           lamplus);
              num_x_flux_cc(consp[0],consp[1],consp[2],consp[3],consp[4],
                          fplus[0],fplus[1],fplus[2],fplus[3],fplus[4],
                          velxplus[idx],pressplus[idx],savg_detr,
                          avg_alp,avg_beta);
              num_x_flux_cc(consm_i[0],consm_i[1],consm_i[2],consm_i[3],consm_i[4],
                          fminus[0],fminus[1],fminus[2],fminus[3],fminus[4],
                          velxminus[idxr],pressminus[idxr],savg_detr,
                          avg_alp,avg_beta);
            } else if(fdir==2) {
              eigenvalues_cc(rhominus[idxr],
                           velyminus[idxr],
                           velzminus[idxr],
                           velxminus[idxr],
                           epsminus[idxr],
                           pressminus[idxr],
                           cs2minus[idxr],
                           w_lorentzminus[idxr],
                           gyyh,gyzh,gxyh,gzzh,gxzh,gxxh,
                           avg_alp,avg_beta,usendh,
                           lamminus);
              eigenvalues_cc(rhoplus[idx],
                           velyplus[idx],
                           velzplus[idx],
                           velxplus[idx],
                           epsplus[idx],
                           pressplus[idx],
                           cs2plus[idx],
                           w_lorentzplus[idx],
                           gyyh,gyzh,gxyh,gzzh,gxzh,gxxh,
                           avg_alp,avg_beta,usendh,
                           lamplus);
              num_x_flux_cc(consp[0],consp[2],consp[3],consp[1],consp[4],
                          fplus[0],fplus[2],fplus[3],fplus[1],fplus[4],
                          velyplus[idx],pressplus[idx],savg_detr,
                          avg_alp,avg_beta);
              num_x_flux_cc(consm_i[0],consm_i[2],consm_i[3],consm_i[1],consm_i[4],
                          fminus[0],fminus[2],fminus[3],fminus[1],fminus[4],
                          velyminus[idxr],pressminus[idxr],savg_detr,
                          avg_alp,avg_beta);
            } else {
              eigenvalues_cc(rhominus[idxr],
                           velzminus[idxr],
                           velxminus[idxr],
                           velyminus[idxr],
                           epsminus[idxr],
                           pressminus[idxr],
                           cs2minus[idxr],
                           w_lorentzminus[idxr],
                           gzzh,gxzh,gyzh,gxxh,gxyh,gyyh,
                           avg_alp,avg_beta,usendh,
                           lamminus);
              eigenvalues_cc(rhoplus[idx],
                           velzplus[idx],
                           velxplus[idx],
                           velyplus[idx],
                           epsplus[idx],
                           pressplus[idx],
                           cs2plus[idx],
                           w_lorentzplus[idx],
                           gzzh,gxzh,gyzh,gxxh,gxyh,gyyh,
                           avg_alp,avg_beta,usendh,
                           lamplus);
              num_x_flux_cc(consp[0],consp[3],consp[1],consp[2],consp[4],
                          fplus[0],fplus[3],fplus[1],fplus[2],fplus[4],
                          velzplus[idx],pressplus[idx],savg_detr,
                          avg_alp,avg_beta);
              num_x_flux_cc(consm_i[0],consm_i[3],consm_i[1],consm_i[2],consm_i[4],
                          fminus[0],fminus[3],fminus[1],fminus[2],fminus[4],
                          velzminus[idxr],pressminus[idxr],savg_detr,
                          avg_alp,avg_beta);
            }
          }
          // H-viscosity
          if(do_hvisc) {
            CCTK_REAL etabar = 0.0;

            if(fdir==1) {
              etabar = max9(etaX(cctkGH,i,j,k,vup,eos_c),
                            etaY(cctkGH,i,j,k,vup,eos_c),
                            etaY(cctkGH,i+1,j,k,vup,eos_c),
                            etaY(cctkGH,i,j-1,k,vup,eos_c),
                            etaY(cctkGH,i+1,j-1,k,vup,eos_c),
                            etaZ(cctkGH,i,j,k,vup,eos_c),
                            etaZ(cctkGH,i+1,j,k,vup,eos_c),
                            etaZ(cctkGH,i,j,k-1,vup,eos_c),
                            etaZ(cctkGH,i+1,j,k-1,vup,eos_c));
            } else if(fdir==2) {
              etabar = max9(etaY(cctkGH,i,j,k,vup,eos_c),
                            etaX(cctkGH,i,j,k,vup,eos_c),
                            etaX(cctkGH,i,j+1,k,vup,eos_c),
                            etaX(cctkGH,i-1,j,k,vup,eos_c),
                            etaX(cctkGH,i-1,j+1,k,vup,eos_c),
                            etaZ(cctkGH,i,j,k,vup,eos_c),
                            etaZ(cctkGH,i,j+1,k,vup,eos_c),
                            etaZ(cctkGH,i,j,k-1,vup,eos_c),
                            etaZ(cctkGH,i,j+1,k-1,vup,eos_c));
            } else if(fdir==3) {
              etabar = max9(etaZ(cctkGH,i,j,k,vup,eos_c),
                            etaX(cctkGH,i,j,k,vup,eos_c),
                            etaX(cctkGH,i,j,k+1,vup,eos_c),
                            etaX(cctkGH,i-1,j,k,vup,eos_c),
                            etaX(cctkGH,i-1,j,k+1,vup,eos_c),
                            etaY(cctkGH,i,j,k,vup,eos_c),
                            etaY(cctkGH,i,j,k+1,vup,eos_c),
                            etaY(cctkGH,i,j-1,k,vup,eos_c),
                            etaY(cctkGH,i,j-1,k+1,vup,eos_c));
            } else {
              CCTK_ERROR("Flux direction not x,y,z");
            }
            // modify eigenvalues of Roe's matrix by computed H viscosity
            for(int m=0;m<10;m++) {
              lamplusminus[m]  = copysign(1.0,lamplusminus[m])*MAX(fabs(lamplusminus[m]),etabar);
              // copysign is C99 and the equivalent of FORTRAN sign(a,b)
            }

          }

          // Find minimum and maximum wavespeeds
          CCTK_REAL charmax = max10(lamplusminus);
          CCTK_REAL charmin = min10(lamplusminus);

          // Introduce some optional dissipation
          // by decreasing/increasing min/max wavespeeds
          // min10/max10 guarantee: charmax >= 0, charmin <= 0
          // Anton Eq. 60: light-like eigenvalues are
          // lambda = +/- sqrt(gxx) - beta^x
          const CCTK_REAL susendh = sqrt(usendh);
          charmin = MAX(-avg_beta/avg_alp-susendh,charmin *
                        (1. + riemann_wavespeed_factor));
          charmax = MIN(-avg_beta/avg_alp+susendh,charmax *
                        (1. + riemann_wavespeed_factor));

          CCTK_REAL icharpm = 1.0 / (charmax-charmin);
          for(int m=0;m<activesize; ++m) {
              f1[m] = (charmax * fplus[m] - charmin * fminus[m] +
                       charmax*charmin * qdiff[m]) * icharpm;
          }

          if (do_clean_divergence) {
            const CCTK_REAL psidcqdiff = psidcm_i - psidcp;
            // The fastest speed for psidc comes from the condition
            // that the normal vector to the characteristic hypersurface
            // be spacelike (Eq. 60 of Anton et al.)

            const CCTK_REAL charmax_dc = +sqrt(usendh) - avg_beta/avg_alp;
            const CCTK_REAL charmin_dc = -sqrt(usendh) - avg_beta/avg_alp;
            const CCTK_REAL charpm_dc = charmax_dc - charmin_dc;

            psidcf1 = (charmax_dc * psidcfplus - charmin_dc * psidcfminus +
                       charmax_dc * charmin_dc * psidcqdiff) / charpm_dc;
            f1[4+fdir] = (charmax_dc * fplus[4+fdir]
                        - charmin_dc * fminus[4+fdir] +
                          charmax_dc * charmin_dc * qdiff[4+fdir]) / charpm_dc;
          }
        } // else trivial

        densflux[idx] = f1[0];
        sxflux[idx]   = f1[1];
        syflux[idx]   = f1[2];
        szflux[idx]   = f1[3];
        tauflux[idx]  = f1[4];
        if (do_mhd) {
          Bconsxflux[idx]   = f1[5];
          Bconsyflux[idx]   = f1[6];
          Bconszflux[idx]   = f1[7];
        }
        if (do_clean_divergence)
          psidcflux[idx] = psidcf1;

        if(do_ye) {
          if(densflux[idx] > 0.0)
            Y_e_con_flux[idx] = Y_e_plus[idx] * densflux[idx];
          else
            Y_e_con_flux[idx] = Y_e_minus[idxr] * densflux[idx];
        }

        if(do_tracers) {
          if(densflux[idx] > 0.0)
            for(int m=0;m<number_of_tracers; ++m) {
              const int idx4 = CCTK_VECTGFINDEX3D(cctkGH,i,j,k,m);
              cons_tracerflux[idx4] = cons_tracerplus[idx4]
                * densflux[idx];
            }
          else
            for(int m=0;m<number_of_tracers; ++m) {
              const int idx4 = CCTK_VECTGFINDEX3D(cctkGH,i,j,k,m);
              const int idx4r =
                CCTK_VECTGFINDEX3D(cctkGH,i+offsetx,j+offsety,k+offsetz,m);
              cons_tracerflux[idx4] = cons_tracerminus[idx4r]
                * densflux[idx];
            }
        } // do tracers


      } // end big loop
}


/**
   Select a particular flavor of HLLE Riemann solver.
*/
template <int dir>
struct riemann {

  static inline void select(const bool apply_Hvisc,
                            const bool evolve_Y_e,
                            const bool evolve_tracer,
                            const bool evolve_MHD,
                            const bool clean_divergence,
                             CCTK_ARGUMENTS)
  {
    if(apply_Hvisc)
      if(evolve_Y_e)
        if(evolve_tracer) 
          if(evolve_MHD) 
            if(clean_divergence)
              GRHydro_HLLE_CC_LL<dir,true,true,true,true,true>(CCTK_PASS_CTOC);
            else
              GRHydro_HLLE_CC_LL<dir,true,true,true,true,false>(CCTK_PASS_CTOC);
          else
            GRHydro_HLLE_CC_LL<dir,true,true,true,false,false>(CCTK_PASS_CTOC);
        else
          if(evolve_MHD) 
            if(clean_divergence)
              GRHydro_HLLE_CC_LL<dir,true,true,false,true,true>(CCTK_PASS_CTOC);
            else
              GRHydro_HLLE_CC_LL<dir,true,true,false,true,false>(CCTK_PASS_CTOC);
          else
            GRHydro_HLLE_CC_LL<dir,true,true,false,false,false>(CCTK_PASS_CTOC);
      else
        if(evolve_tracer)
          if(evolve_MHD) 
            if(clean_divergence)
              GRHydro_HLLE_CC_LL<dir,false,true,true,true,true>(CCTK_PASS_CTOC);
            else
              GRHydro_HLLE_CC_LL<dir,false,true,true,true,false>(CCTK_PASS_CTOC);
          else
            GRHydro_HLLE_CC_LL<dir,false,true,true,false,false>(CCTK_PASS_CTOC);
        else
          if(evolve_MHD) 
            if(clean_divergence)
              GRHydro_HLLE_CC_LL<dir,false,true,false,true,true>(CCTK_PASS_CTOC);
            else
              GRHydro_HLLE_CC_LL<dir,false,true,false,true,false>(CCTK_PASS_CTOC);
          else
            GRHydro_HLLE_CC_LL<dir,false,true,false,false,false>(CCTK_PASS_CTOC);
    else
      if(evolve_Y_e)
        if(evolve_tracer)
          if(evolve_MHD) 
            if(clean_divergence)
              GRHydro_HLLE_CC_LL<dir,true,false,true,true,true>(CCTK_PASS_CTOC);
            else
              GRHydro_HLLE_CC_LL<dir,true,false,true,true,false>(CCTK_PASS_CTOC);
          else
            GRHydro_HLLE_CC_LL<dir,true,false,true,false,false>(CCTK_PASS_CTOC);
        else
          if(evolve_MHD) 
            if(clean_divergence)
              GRHydro_HLLE_CC_LL<dir,true,false,false,true,true>(CCTK_PASS_CTOC);
            else
              GRHydro_HLLE_CC_LL<dir,true,false,false,true,false>(CCTK_PASS_CTOC);
          else
            GRHydro_HLLE_CC_LL<dir,true,false,false,false,false>(CCTK_PASS_CTOC);
      else
        if(evolve_tracer) 
          if(evolve_MHD) 
            if(clean_divergence)
              GRHydro_HLLE_CC_LL<dir,false,false,true,true,true>(CCTK_PASS_CTOC);
            else
              GRHydro_HLLE_CC_LL<dir,false,false,true,true,false>(CCTK_PASS_CTOC);
          else
            GRHydro_HLLE_CC_LL<dir,false,false,true,false,false>(CCTK_PASS_CTOC);
        else
          if(evolve_MHD) 
            if(clean_divergence)
              GRHydro_HLLE_CC_LL<dir,false,false,false,true,true>(CCTK_PASS_CTOC);
            else
              GRHydro_HLLE_CC_LL<dir,false,false,false,true,false>(CCTK_PASS_CTOC);
          else
            GRHydro_HLLE_CC_LL<dir,false,false,false,false,false>(CCTK_PASS_CTOC);
  }

};


void GRHydro_HLLE_CC(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // flux direction 1
  if(*flux_direction==1)

    riemann<1>::select(apply_H_viscosity, *evolve_Y_e, evolve_tracer, *evolve_MHD, clean_divergence, CCTK_PASS_CTOC);

  // flux direction 2
  else if(*flux_direction==2)

    riemann<2>::select(apply_H_viscosity, *evolve_Y_e, evolve_tracer, *evolve_MHD, clean_divergence, CCTK_PASS_CTOC);

    // flux direction 3
  else if(*flux_direction==3)

    riemann<3>::select(apply_H_viscosity, *evolve_Y_e, evolve_tracer, *evolve_MHD, clean_divergence, CCTK_PASS_CTOC);

  else
    CCTK_ERROR("Illegal flux direction!");

  return;
}

extern "C"
CCTK_FCALL void CCTK_FNAME(GRHydro_HLLE_CC_F2C)(cGH ** p_cctkGH) {
  GRHydro_HLLE_CC(*p_cctkGH);
}


extern "C"
void H_viscosity_calc_cs_cc(CCTK_ARGUMENTS) {

  // This is a preperatory routine for H viscosity.
  // All it does is fill a speed of sound helper array

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const size_t n = size_t(cctk_ash[0]*cctk_ash[1]*cctk_ash[2]);
  std::vector<CCTK_INT> keyerr(n);
  CCTK_INT anyerr = 0;
  int keytemp = 0;

  if(!*evolve_temper) {
#pragma omp parallel for
    for(int k=0;k<cctk_lsh[2];k++)
      for(int j=0;j<cctk_lsh[1];j++) {
        const int i = CCTK_GFINDEX3D(cctkGH,0,j,k);
        CCTK_REAL cs2[cctk_ash[0]];
        EOS_Omni_cs2(*GRHydro_eos_handle,keytemp,GRHydro_eos_rf_prec,cctk_lsh[0],
                     &(rho[i]),&(eps[i]),NULL,NULL,&cs2[i],
                     &(keyerr[i]),&anyerr);

        for(int ii=0;ii<cctk_lsh[0];++ii) {
          eos_c[ii+i] = sqrt(cs2[ii]);
        }
      }
  } else {
#pragma omp parallel for
    for(int k=0;k<cctk_lsh[2];k++)
      for(int j=0;j<cctk_lsh[1];j++) {
        const int i = CCTK_GFINDEX3D(cctkGH,0,j,k);
        CCTK_REAL cs2[cctk_ash[0]];
        EOS_Omni_cs2(*GRHydro_eos_handle,keytemp,GRHydro_eos_rf_prec,cctk_lsh[0],
                     &rho[i],&eps[i],&temperature[i],&Y_e[i],
                     &cs2[i],&keyerr[i],&anyerr);
        for(int ii=0;ii<cctk_lsh[0];++ii) {
          eos_c[ii+i] = sqrt(cs2[ii]);
        }
      } // for
    if(anyerr) {
      for(size_t i=0;i<n;i++) {
        if(keyerr[i] != 0) {
#pragma omp critical
          {
            CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                       "rl: %d i,x,y,z: %d %15.6E %15.6E %15.6E, keyerr: %d",
                       int(*GRHydro_reflevel), int(i), x[i], y[i], z[i], int(keyerr[i]));
            CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                       "rl: %d r,t,ye: %15.6E %15.6E %15.6E, keyerr: %d",
                       int(*GRHydro_reflevel), rho[i], temperature[i], Y_e[i], int(keyerr[i]));
          }
        }
      } // for i, i<n
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Aborting!");
    } // anyerr
  } // if evolve_temper

}

static inline void UpperMetricC(CCTK_REAL& restrict uxx,
                                CCTK_REAL& restrict uxy,
                                CCTK_REAL& restrict uxz,
                                CCTK_REAL& restrict uyy,
                                CCTK_REAL& restrict uyz,
                                CCTK_REAL& restrict uzz,
                                const CCTK_REAL det,
                                const CCTK_REAL gxx,
                                const CCTK_REAL gxy,
                                const CCTK_REAL gxz,
                                const CCTK_REAL gyy,
                                const CCTK_REAL gyz,
                                const CCTK_REAL gzz)
{
  const CCTK_REAL invdet = 1.0 / det;
  uxx = (-gyz*gyz + gyy*gzz)*invdet;
  uxy = (gxz*gyz - gxy*gzz)*invdet;
  uxz = (-gxz*gyy + gxy*gyz)*invdet;
  uyy = (-gxz*gxz + gxx*gzz)*invdet;
  uyz = (gxy*gxz - gxx*gyz)*invdet;
  uzz = (-gxy*gxy + gxx*gyy)*invdet;

  return;
}

static inline CCTK_REAL SpatialDeterminantC(const CCTK_REAL gxx, const CCTK_REAL gxy,
                                            const CCTK_REAL gxz, const CCTK_REAL gyy,
                                            const CCTK_REAL gyz, const CCTK_REAL gzz)
{
  return -gxz*gxz*gyy + 2.0*gxy*gxz*gyz - gxx*gyz*gyz
    - gxy*gxy*gzz + gxx*gyy*gzz;
}

