 /*@@
   @file      GRHydro_Source.cc
   @date      Nov 29, 2013
   @author    Christian Reisswig, Joshua Faber, Scott Noble, Bruno Mundim, Ian Hawke
   @desc 
   The geometric source terms for the matter evolution
   @enddesc 
 @@*/




#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

#include <cstring>
#include <cassert>

#define velx vup
#define vely (&vup[N])
#define velz (&vup[2*N])
#define Bvecx Bprim
#define Bvecy (&Bprim[N])
#define Bvecz (&Bprim[2*N])
#define Avecx Avec
#define Avecy (&Avec[N])
#define Avecz (&Avec[2*N])
#define Avecrhsx Avecrhs
#define Avecrhsy (&Avecrhs[N])
#define Avecrhsz (&Avecrhs[2*N])
#define Bconsx Bcons
#define Bconsy (&Bcons[N])
#define Bconsz (&Bcons[2*N])
#define Bconsrhsx Bconsrhs
#define Bconsrhsy (&Bconsrhs[N])
#define Bconsrhsz (&Bconsrhs[2*N])


extern "C" CCTK_INT GRHydro_UseGeneralCoordinates(const cGH * cctkGH);

template <typename T> static inline T SQR (T const & x) { return x*x; }


struct alldiff2 {
  static void apply(const cGH* const restrict cctkGH,
                    CCTK_REAL dvars[][3],
                    const CCTK_REAL* const restrict * const restrict vars, 
                    const int i, const int j, const int k,
                    const CCTK_REAL* restrict const ih,
                    const int nvars)
  {
     for (int n=0; n < nvars; ++n) {
        dvars[n][0] = (vars[n][CCTK_GFINDEX3D(cctkGH, i+1,j,k)] - vars[n][CCTK_GFINDEX3D(cctkGH, i-1,j,k)]) * 0.5 * ih[0];
        dvars[n][1] = (vars[n][CCTK_GFINDEX3D(cctkGH, i,j+1,k)] - vars[n][CCTK_GFINDEX3D(cctkGH, i,j-1,k)]) * 0.5 * ih[1];
        dvars[n][2] = (vars[n][CCTK_GFINDEX3D(cctkGH, i,j,k+1)] - vars[n][CCTK_GFINDEX3D(cctkGH, i,j,k-1)]) * 0.5 * ih[2];
     }
  }
};

struct alldiff4 {
  static void apply(const cGH* const restrict cctkGH,
                    CCTK_REAL dvars[][3],
                    const CCTK_REAL* const restrict * const restrict vars, 
                    const int i, const int j, const int k,
                    const CCTK_REAL* restrict const ih,
                    const int nvars)
  {
     for (int n=0; n < nvars; ++n) {
        dvars[n][0] = (        vars[n][CCTK_GFINDEX3D(cctkGH, i-2,j,k)] - vars[n][CCTK_GFINDEX3D(cctkGH, i+2,j,k)]
                      + 8.0 * (vars[n][CCTK_GFINDEX3D(cctkGH, i+1,j,k)] - vars[n][CCTK_GFINDEX3D(cctkGH, i-1,j,k)])) * (1.0 / 12.0) * ih[0];
        dvars[n][1] = (        vars[n][CCTK_GFINDEX3D(cctkGH, i,j-2,k)] - vars[n][CCTK_GFINDEX3D(cctkGH, i,j+2,k)]
                      + 8.0 * (vars[n][CCTK_GFINDEX3D(cctkGH, i,j+1,k)] - vars[n][CCTK_GFINDEX3D(cctkGH, i,j-1,k)])) * (1.0 / 12.0) * ih[1];
        dvars[n][2] = (        vars[n][CCTK_GFINDEX3D(cctkGH, i,j,k-2)] - vars[n][CCTK_GFINDEX3D(cctkGH, i,j,k+2)]
                      + 8.0 * (vars[n][CCTK_GFINDEX3D(cctkGH, i,j,k+1)] - vars[n][CCTK_GFINDEX3D(cctkGH, i,j,k-1)])) * (1.0 / 12.0) * ih[2];
     }
  }
};


static inline void UpperMetric(CCTK_REAL& restrict uxx,
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


template <class alldiff,
          bool do_mhd,
          bool do_clean_divergence,
          bool do_divergence_flux,
          bool do_Avec>
static void SourceTerms_LL(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  /*
     Set up multipatch stuff
  */

  const CCTK_REAL * restrict vup;
  const CCTK_REAL * restrict Bprim;
  const CCTK_REAL * restrict g11;
  const CCTK_REAL * restrict g12;
  const CCTK_REAL * restrict g13;
  const CCTK_REAL * restrict g22;
  const CCTK_REAL * restrict g23;
  const CCTK_REAL * restrict g33;
  const CCTK_REAL * restrict k11;
  const CCTK_REAL * restrict k12;
  const CCTK_REAL * restrict k13;
  const CCTK_REAL * restrict k22;
  const CCTK_REAL * restrict k23;
  const CCTK_REAL * restrict k33;
  const CCTK_REAL * restrict beta1;
  const CCTK_REAL * restrict beta2;
  const CCTK_REAL * restrict beta3;

  //Multipatch related pointers
  if(GRHydro_UseGeneralCoordinates(cctkGH)) {
    vup=lvel;
    g11=gaa; g12=gab; g13=gac;
    g22=gbb; g23=gbc;
    g33=gcc;
    k11=kaa; k12=kab; k13=kac;
    k22=kbb; k23=kbc;
    k33=kcc;
    Bprim=lBvec;
    beta1 = betaa;
    beta2 = betab;
    beta3 = betac;
  } else {
    vup=vel;
    g11=gxx; g12=gxy; g13=gxz;
    g22=gyy; g23=gyz;
    g33=gzz;
    k11=kxx; k12=kxy; k13=kxz;
    k22=kyy; k23=kyz;
    k33=kzz;
    Bprim=Bvec;
    beta1 = betax;
    beta2 = betay;
    beta3 = betaz;
  }

  // TODO: for vectorization we might want to use cctk_ash (have to make sure
  // sensible values exist everywhere)
  const int nx=cctk_lsh[0];
  const int ny=cctk_lsh[1];
  const int nz=cctk_lsh[2];

  const int N = cctk_ash[0]*cctk_ash[1]*cctk_ash[2];

  const CCTK_REAL one = 1.00;
  const CCTK_REAL two = 2.00;
  const CCTK_REAL half = 0.50;
  const CCTK_REAL dx = CCTK_DELTA_SPACE(0);
  const CCTK_REAL dy = CCTK_DELTA_SPACE(1);
  const CCTK_REAL dz = CCTK_DELTA_SPACE(2);
  
  const CCTK_REAL ih[3] = { 1.0/dx, 1.0/dy, 1.0/dz };

  memset(densrhs, 0, sizeof(densrhs[0])*N);
  assert(densrhs[0] == 0.);
  memset(srhs, 0, sizeof(srhs[0])*3*N);
  assert(srhs[0] == 0.);
  memset(taurhs, 0, sizeof(taurhs[0])*N);
  assert(taurhs[0] == 0.);
 
#pragma omp parallel for
  for (int k=GRHydro_stencil; k < nz-GRHydro_stencil; ++k) {
    for (int j=GRHydro_stencil; j < ny-GRHydro_stencil; ++j) {
      for (int i=GRHydro_stencil; i < nx-GRHydro_stencil; ++i) {
        
        const int ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
        
        // TODO: these are not needed.
        const CCTK_REAL localgxx = g11[ijk];
        const CCTK_REAL localgxy = g12[ijk];
        const CCTK_REAL localgxz = g13[ijk];
        const CCTK_REAL localgyy = g22[ijk];
        const CCTK_REAL localgyz = g23[ijk];
        const CCTK_REAL localgzz = g33[ijk];

        const CCTK_REAL sqrtdet = sdetg[ijk];
        const CCTK_REAL invsqrtdet = 1./sqrtdet;
        
        CCTK_REAL uxx, uxy, uxz, uyy, uyz, uzz;
        UpperMetric(uxx, uxy, uxz, uyy, uyz, uzz, sqrtdet*sqrtdet, localgxx,
                    localgxy, localgxz, localgyy, localgyz, localgzz);
        
        const CCTK_REAL shiftx = beta1[ijk];
        const CCTK_REAL shifty = beta2[ijk];
        const CCTK_REAL shiftz = beta3[ijk];

        //  Derivatives of the lapse, metric and shift

        const int nvars = 10;
        const CCTK_REAL* const restrict vars[nvars] = { beta1, beta2, beta3, alp,
                                                     g11, g12, g13,
                                                     g22, g23,
                                                     g33 };
        CCTK_REAL dvars[nvars][3];
        
        alldiff::apply(cctkGH, dvars, vars, i, j, k, ih, nvars);
        
        const CCTK_REAL dx_betax = dvars[0][0];
        const CCTK_REAL dx_betay = dvars[1][0];
        const CCTK_REAL dx_betaz = dvars[2][0];
           
        const CCTK_REAL dy_betax = dvars[0][1];
        const CCTK_REAL dy_betay = dvars[1][1];
        const CCTK_REAL dy_betaz = dvars[2][1];
           
        const CCTK_REAL dz_betax = dvars[0][2];
        const CCTK_REAL dz_betay = dvars[1][2];
        const CCTK_REAL dz_betaz = dvars[2][2];
        
        const CCTK_REAL dx_alp = dvars[3][0];
        const CCTK_REAL dy_alp = dvars[3][1];
        const CCTK_REAL dz_alp = dvars[3][2];

        const CCTK_REAL dx_gxx = dvars[4][0];
        const CCTK_REAL dx_gxy = dvars[5][0];
        const CCTK_REAL dx_gxz = dvars[6][0];
        const CCTK_REAL dx_gyy = dvars[7][0];
        const CCTK_REAL dx_gyz = dvars[8][0];
        const CCTK_REAL dx_gzz = dvars[9][0];
        const CCTK_REAL dy_gxx = dvars[4][1];
        const CCTK_REAL dy_gxy = dvars[5][1];
        const CCTK_REAL dy_gxz = dvars[6][1];
        const CCTK_REAL dy_gyy = dvars[7][1];
        const CCTK_REAL dy_gyz = dvars[8][1];
        const CCTK_REAL dy_gzz = dvars[9][1];
        const CCTK_REAL dz_gxx = dvars[4][2];
        const CCTK_REAL dz_gxy = dvars[5][2];
        const CCTK_REAL dz_gxz = dvars[6][2];
        const CCTK_REAL dz_gyy = dvars[7][2];
        const CCTK_REAL dz_gyz = dvars[8][2];
        const CCTK_REAL dz_gzz = dvars[9][2];


        const CCTK_REAL invalp = 1.0 / alp[ijk];
        const CCTK_REAL invalp2 = SQR(invalp);
        const CCTK_REAL velxshift = velx[ijk] - shiftx*invalp;
        const CCTK_REAL velyshift = vely[ijk] - shifty*invalp;
        const CCTK_REAL velzshift = velz[ijk] - shiftz*invalp;

        // vel_i  = g_ij v^j
        // B_i = g_ij B^i

        const CCTK_REAL vlowx = g11[ijk]*velx[ijk] + g12[ijk]*vely[ijk] + g13[ijk]*velz[ijk];
        const CCTK_REAL vlowy = g12[ijk]*velx[ijk] + g22[ijk]*vely[ijk] + g23[ijk]*velz[ijk];
        const CCTK_REAL vlowz = g13[ijk]*velx[ijk] + g23[ijk]*vely[ijk] + g33[ijk]*velz[ijk];
        const CCTK_REAL Bvecxlow = do_mhd ? g11[ijk]*Bvecx[ijk] + g12[ijk]*Bvecy[ijk] + g13[ijk]*Bvecz[ijk] : 0;
        const CCTK_REAL Bvecylow = do_mhd ? g12[ijk]*Bvecx[ijk] + g22[ijk]*Bvecy[ijk] + g23[ijk]*Bvecz[ijk] : 0;
        const CCTK_REAL Bveczlow = do_mhd ? g13[ijk]*Bvecx[ijk] + g23[ijk]*Bvecy[ijk] + g33[ijk]*Bvecz[ijk] : 0;

        //  B^i v_i (= b^0/u^0)
        const CCTK_REAL Bdotv = do_mhd ? vlowx*Bvecx[ijk]+vlowy*Bvecy[ijk]+vlowz*Bvecz[ijk] : 0;

        // v^2 = v_i v^i; w=(1-v^2)^{-1/2}

        const CCTK_REAL v2 = vlowx*velx[ijk] + vlowy*vely[ijk] + vlowz*velz[ijk];
        const CCTK_REAL invw = sqrt(1.0-v2);
        const CCTK_REAL w = 1./invw;

        // b^2 = B^i B_i / w^2 + (b^0/u^0)^2

        const CCTK_REAL b2 = do_mhd ? (Bvecx[ijk]*Bvecxlow+Bvecy[ijk]*Bvecylow+Bvecz[ijk]*Bveczlow)*SQR(invw)+SQR(Bdotv) : 0;

        // b_i = B_i/w +w*(B dot v)*v_i
        const CCTK_REAL bxlow = do_mhd ? Bvecxlow*invw+w*Bdotv*vlowx : 0;
        const CCTK_REAL bylow = do_mhd ? Bvecylow*invw+w*Bdotv*vlowy : 0;
        const CCTK_REAL bzlow = do_mhd ? Bveczlow*invw+w*Bdotv*vlowz : 0;


        // These are the contravariant components
        const CCTK_REAL bt = do_mhd ? w*invalp*Bdotv : 0;
        const CCTK_REAL bx = do_mhd ? Bvecx[ijk]*invw+w*Bdotv*velxshift : 0;
        const CCTK_REAL by = do_mhd ? Bvecy[ijk]*invw+w*Bdotv*velyshift : 0;
        const CCTK_REAL bz = do_mhd ? Bvecz[ijk]*invw+w*Bdotv*velzshift : 0;

        // TODO: all of these can be expressed much more easily in terms of the
        // conservatives
        const CCTK_REAL rhohstarW2 = (rho[ijk]*(one + eps[ijk]) + (press[ijk] + b2)) *
                                       SQR(w);
        const CCTK_REAL pstar = press[ijk]+0.50*b2;

        //  For a change, these are T^{ij}

        // TODO: strict IEEE compliance does not let the compiler remove "+ 0"
        // terms, so we have to do something else here
        const CCTK_REAL t00 = (rhohstarW2 - pstar)*invalp2-SQR(bt);
        const CCTK_REAL t0x = rhohstarW2*velxshift*invalp +
             pstar*shiftx*invalp2-bt*bx;
        const CCTK_REAL t0y = rhohstarW2*velyshift*invalp +
             pstar*shifty*invalp2-bt*by;
        const CCTK_REAL t0z = rhohstarW2*velzshift*invalp +
             pstar*shiftz*invalp2-bt*bz;
        const CCTK_REAL txx = rhohstarW2*velxshift*velxshift +
             pstar*(uxx - shiftx*shiftx*invalp2)-SQR(bx);
        const CCTK_REAL txy = rhohstarW2*velxshift*velyshift +
             pstar*(uxy - shiftx*shifty*invalp2)-bx*by;
        const CCTK_REAL txz = rhohstarW2*velxshift*velzshift +
             pstar*(uxz - shiftx*shiftz*invalp2)-bx*bz;
        const CCTK_REAL tyy = rhohstarW2*velyshift*velyshift +
             pstar*(uyy - shifty*shifty*invalp2)-SQR(by);
        const CCTK_REAL tyz = rhohstarW2*velyshift*velzshift +
             pstar*(uyz - shifty*shiftz*invalp2)-by*bz;
        const CCTK_REAL tzz = rhohstarW2*velzshift*velzshift +
             pstar*(uzz - shiftz*shiftz*invalp2)-SQR(bz);

//        Contract the shift with the extrinsic curvature

        const CCTK_REAL shiftshiftk = shiftx*shiftx*k11[ijk] +
                                      shifty*shifty*k22[ijk] +
                                      shiftz*shiftz*k33[ijk] +
             two*(shiftx*shifty*k12[ijk] +
                  shiftx*shiftz*k13[ijk] +
                  shifty*shiftz*k23[ijk]);

        const CCTK_REAL shiftkx = shiftx*k11[ijk] + shifty*k12[ijk] + shiftz*k13[ijk];
        const CCTK_REAL shiftky = shiftx*k12[ijk] + shifty*k22[ijk] + shiftz*k23[ijk];
        const CCTK_REAL shiftkz = shiftx*k13[ijk] + shifty*k23[ijk] + shiftz*k33[ijk];

//        Contract the matter terms with the extrinsic curvature

        const CCTK_REAL sumTK = txx*k11[ijk] + tyy*k22[ijk] + tzz*k33[ijk]
                         + two*(txy*k12[ijk] + txz*k13[ijk] + tyz*k23[ijk]);

//        Update term for tau
        
        const CCTK_REAL tau_source = t00*
             (shiftshiftk - (shiftx*dx_alp + shifty*dy_alp + shiftz*dz_alp) )
             + t0x*(-dx_alp + two*shiftkx)
             + t0y*(-dy_alp + two*shiftky)
             + t0z*(-dz_alp + two*shiftkz)
             + sumTK;

//        The following looks very little like the terms in the
//        standard papers. Take a look in the ThornGuide to see why
//        it is really the same thing.

//        Contract the shift with derivatives of the metric

        const CCTK_REAL halfshiftdgx = half*(shiftx*shiftx*dx_gxx +
             shifty*shifty*dx_gyy + shiftz*shiftz*dx_gzz) +
             shiftx*shifty*dx_gxy + shiftx*shiftz*dx_gxz +
             shifty*shiftz*dx_gyz;
        const CCTK_REAL halfshiftdgy = half*(shiftx*shiftx*dy_gxx +
             shifty*shifty*dy_gyy + shiftz*shiftz*dy_gzz) +
             shiftx*shifty*dy_gxy + shiftx*shiftz*dy_gxz +
             shifty*shiftz*dy_gyz;
        const CCTK_REAL halfshiftdgz = half*(shiftx*shiftx*dz_gxx +
             shifty*shifty*dz_gyy + shiftz*shiftz*dz_gzz) +
             shiftx*shifty*dz_gxy + shiftx*shiftz*dz_gxz +
             shifty*shiftz*dz_gyz;

//        Contract the matter with derivatives of the metric

        const CCTK_REAL halfTdgx = half*(txx*dx_gxx + tyy*dx_gyy + tzz*dx_gzz) +
             txy*dx_gxy + txz*dx_gxz + tyz*dx_gyz;
        const CCTK_REAL halfTdgy = half*(txx*dy_gxx + tyy*dy_gyy + tzz*dy_gzz) +
             txy*dy_gxy + txz*dy_gxz + tyz*dy_gyz;
        const CCTK_REAL halfTdgz = half*(txx*dz_gxx + tyy*dz_gyy + tzz*dz_gzz) +
             txy*dz_gxy + txz*dz_gxz + tyz*dz_gyz;

     
       const CCTK_REAL sx_source = t00*
             (halfshiftdgx - alp[ijk]*dx_alp) + halfTdgx +
             t0x*(shiftx*dx_gxx + shifty*dx_gxy + shiftz*dx_gxz) +
             t0y*(shiftx*dx_gxy + shifty*dx_gyy + shiftz*dx_gyz) +
             t0z*(shiftx*dx_gxz + shifty*dx_gyz + shiftz*dx_gzz) +
             rhohstarW2*invalp*(vlowx*dx_betax + vlowy*dx_betay + vlowz*dx_betaz) -
             bt*(bxlow*dx_betax + bylow*dx_betay + bzlow*dx_betaz);
        
       const CCTK_REAL sy_source = t00*
             (halfshiftdgy - alp[ijk]*dy_alp) + halfTdgy +
             t0x*(shiftx*dy_gxx + shifty*dy_gxy + shiftz*dy_gxz) +
             t0y*(shiftx*dy_gxy + shifty*dy_gyy + shiftz*dy_gyz) +
             t0z*(shiftx*dy_gxz + shifty*dy_gyz + shiftz*dy_gzz) +
             rhohstarW2*invalp*(vlowx*dy_betax + vlowy*dy_betay + vlowz*dy_betaz) -
             bt*(bxlow*dy_betax + bylow*dy_betay + bzlow*dy_betaz);

       const CCTK_REAL sz_source = t00*
             (halfshiftdgz - alp[ijk]*dz_alp) + halfTdgz +
             t0x*(shiftx*dz_gxx + shifty*dz_gxy + shiftz*dz_gxz) +
             t0y*(shiftx*dz_gxy + shifty*dz_gyy + shiftz*dz_gyz) +
             t0z*(shiftx*dz_gxz + shifty*dz_gyz + shiftz*dz_gzz) +
             rhohstarW2*invalp*(vlowx*dz_betax + vlowy*dz_betay + vlowz*dz_betaz) -
             bt*(bxlow*dz_betax + bylow*dz_betay + bzlow*dz_betaz);

        densrhs[ijk] = 0.0;
        srhs[ijk]        = alp[ijk]*sqrtdet*sx_source;
        srhs[ijk + N]    = alp[ijk]*sqrtdet*sy_source;
        srhs[ijk + 2*N]  = alp[ijk]*sqrtdet*sz_source;
        taurhs[ijk]      = alp[ijk]*sqrtdet*tau_source;

        if (do_Avec) {

          // B^i and A^i both live in cell centers currently
          const CCTK_REAL Avcx_source = alp[ijk]*sqrtdet*(velyshift*Bvecz[ijk] - velzshift*Bvecy[ijk]);
          const CCTK_REAL Avcy_source = alp[ijk]*sqrtdet*(velzshift*Bvecx[ijk] - velxshift*Bvecz[ijk]);
          const CCTK_REAL Avcz_source = alp[ijk]*sqrtdet*(velxshift*Bvecy[ijk] - velyshift*Bvecx[ijk]);

          Avecrhsx[ijk] = Avcx_source;
          Avecrhsy[ijk] = Avcy_source;
          Avecrhsz[ijk] = Avcz_source;

        }

        if(do_clean_divergence) {
   
           // g^{jk} d_i g_{kj} = d_i (g) / det
           const CCTK_REAL dx_det_bydet = uxx*dx_gxx + uyy*dx_gyy + uzz*dx_gzz +
                two*(uxy*dx_gxy+uxz*dx_gxz+uyz*dx_gyz);
           const CCTK_REAL dy_det_bydet = uxx*dy_gxx + uyy*dy_gyy + uzz*dy_gzz +
                two*(uxy*dy_gxy+uxz*dy_gxz+uyz*dy_gyz);
           const CCTK_REAL dz_det_bydet = uxx*dz_gxx + uyy*dz_gyy + uzz*dz_gzz +
                two*(uxy*dz_gxy+uxz*dz_gxz+uyz*dz_gyz);

           // g^{ik} d_k g_{li}
           const CCTK_REAL gdg_x = uxx*dx_gxx + uxy*dy_gxx + uxz*dz_gxx +
                   uxy*dx_gxy + uyy*dy_gxy + uyz*dz_gxy +
                   uxz*dx_gxz + uyz*dy_gxz + uzz*dz_gxz;

           const CCTK_REAL gdg_y = uxx*dx_gxy + uxy*dy_gxy + uxz*dz_gxy +
                   uxy*dx_gyy + uyy*dy_gyy + uyz*dz_gyy +
                   uxz*dx_gyz + uyz*dy_gyz + uzz*dz_gyz;

           const CCTK_REAL gdg_z = uxx*dx_gxz + uxy*dy_gxz + uxz*dz_gxz +
                   uxy*dx_gyz + uyy*dy_gyz + uyz*dz_gyz +
                   uxz*dx_gzz + uyz*dy_gzz + uzz*dz_gzz;

           CCTK_REAL bvcx_source, bvcy_source, bvcz_source;
           if(do_divergence_flux) {
             psidcrhs[ijk] = -one * (kap_dc*alp[ijk] + 
                  dx_betax + dy_betay + dz_betaz ) * psidc[ijk] + 
                  Bconsx[ijk] * (dx_alp - half*alp[ijk] * 
                   ( uxx*dx_gxx + uyy*dx_gyy + uzz*dx_gzz + two*uxy*dx_gxy + 
                     two*uxz*dx_gxz + two*uyz*dx_gyz ) )*invsqrtdet + 
                  Bconsy[ijk] * (dy_alp - half*alp[ijk] * 
                   ( uxx*dy_gxx + uyy*dy_gyy + uzz*dy_gzz + two*uxy*dy_gxy + 
                     two*uxz*dy_gxz + two*uyz*dy_gyz ) )*invsqrtdet + 
                  Bconsz[ijk] * (dz_alp - half*alp[ijk] * 
                   ( uxx*dz_gxx + uyy*dz_gyy + uzz*dz_gzz + two*uxy*dz_gxy + 
                     two*uxz*dz_gxz + two*uyz*dz_gyz ) )*invsqrtdet;

             bvcx_source = -one * ( Bconsx[ijk]*dx_betax + 
                  Bconsy[ijk]*dy_betax + Bconsz[ijk]*dz_betax ) + 
                  psidc[ijk]*sqrtdet*(( uxx*dx_alp+uxy*dy_alp+uxz*dz_alp ) + 
                  alp[ijk]*(half*( uxx*dx_det_bydet + 
                    uxy*dy_det_bydet + uxz*dz_det_bydet) - 
                  ( uxx*gdg_x + uxy*gdg_y + uxz*gdg_z )));

             bvcy_source = -one * ( Bconsx[ijk]*dx_betay + 
                  Bconsy[ijk]*dy_betay + Bconsz[ijk]*dz_betay ) + 
                  psidc[ijk]*sqrtdet*(( uxy*dx_alp+uyy*dy_alp+uyz*dz_alp ) + 
                  alp[ijk]*(half*( uxy*dx_det_bydet + 
                    uyy*dy_det_bydet + uyz*dz_det_bydet ) - 
                  ( uxy*gdg_x + uyy*gdg_y + uyz*gdg_z )));

             bvcz_source = -one * ( Bconsx[ijk]*dx_betaz + 
                  Bconsy[ijk]*dy_betaz + Bconsz[ijk]*dz_betaz ) + 
                  psidc[ijk]*sqrtdet*(( uxz*dx_alp+uyz*dy_alp+uzz*dz_alp ) + 
                  alp[ijk]*(half*( uxz*dx_det_bydet + 
                    uyz*dy_det_bydet + uzz*dz_det_bydet ) - 
                  ( uxz*gdg_x + uyz*gdg_y + uzz*gdg_z )));
           } else {
            const int dcnvars = 4;
            const CCTK_REAL* const restrict dcvars[dcnvars] =
              { Bconsx, Bconsy, Bconsz, psidc };
            CCTK_REAL dcdvars[dcnvars][3];
            
            // we really only need divB and d_i psidc
            alldiff::apply(cctkGH, dcdvars, dcvars, i, j, k, ih, dcnvars);
            
            const CCTK_REAL dx_bconsx = dcdvars[0][0];
            const CCTK_REAL dy_bconsy = dcdvars[1][1];
            const CCTK_REAL dz_bconsz = dcdvars[2][2];
            
            const CCTK_REAL dx_psi = dcdvars[3][0];
            const CCTK_REAL dy_psi = dcdvars[3][1];
            const CCTK_REAL dz_psi = dcdvars[3][2];

             psidcrhs[ijk] = -one*kap_dc*alp[ijk]*psidc[ijk] + beta1[ijk]*dx_psi + beta2[ijk]*dy_psi + beta3[ijk]*dz_psi -
                               alp[ijk]*invsqrtdet*(dx_bconsx + dy_bconsy+dz_bconsz);

             bvcx_source = beta1[ijk]*dx_bconsx + beta1[ijk]*dy_bconsy + beta1[ijk]*dz_bconsz - 
                           alp[ijk]*sqrtdet*(uxx*dx_psi+uxy*dy_psi+uxz*dz_psi);
             bvcy_source = beta2[ijk]*dx_bconsx + beta2[ijk]*dy_bconsy + beta2[ijk]*dz_bconsz - 
                           alp[ijk]*sqrtdet*(uxy*dx_psi+uyy*dy_psi+uyz*dz_psi);
             bvcz_source = beta3[ijk]*dx_bconsx + beta3[ijk]*dy_bconsy + beta3[ijk]*dz_bconsz - 
                           alp[ijk]*sqrtdet*(uxz*dx_psi+uyz*dy_psi+uzz*dz_psi);
           }

           Bconsrhsx[ijk] = bvcx_source;
           Bconsrhsy[ijk] = bvcy_source;
           Bconsrhsz[ijk] = bvcz_source;
        } // if(do_clean_divergence)

      } // for(i
    } // for(j
  } // for(k

#if(0) // poison edges of domain
  static int last_iteration_seen = -1;
  static int reflevel = -1;
  static int mol_substep;
  if(last_iteration_seen != cctk_iteration || reflevel != *GRHydro_reflevel) {
    last_iteration_seen = cctk_iteration;
    reflevel = *GRHydro_reflevel;
    mol_substep = 0;
  } else {
    mol_substep += 1;
  }
  for(k = 0 ; k < GRHydro_stencil*mol_substep ; ++k) {
    for(j = 0 ; j < cctk_ash[1] ; ++j) {
      for(i = 0 ; i < cctk_ash[0] ; ++i) {
        dens[CCTK_GFINDEX3D(i,j,k)] = -1e100
        Scon[CCTK_VECTGFINDEX3D(i,j,k,0)] = -1e100
        Scon[CCTK_VECTGFINDEX3D(i,j,k,1)] = -1e100
        Scon[CCTK_VECTGFINDEX3D(i,j,k,2)] = -1e100
        tau[CCTK_GFINDEX3D(i,j,k)] = -1e100
        if(do_mhd) {
          Bcons[CCTK_VECTGFINDEX3D(i,j,k,0)] = -1e100
          Bcons[CCTK_VECTGFINDEX3D(i,j,k,1)] = -1e100
          Bcons[CCTK_VECTGFINDEX3D(i,j,k,2)] = -1e100
        }
        if(do_clean_divergene)
          psidc[CCTK_GFINDEX3D(i,j,k)] = -1e100
      }
    }
  }
  for(k = cctk_ash[2]-GRHydro_stencil*mol_substep ; k < cctk_ash[2] ; ++k) {
    for(j = 0 ; j < cctk_ash[1] ; ++j) {
      for(i = 0 ; i < cctk_ash[0] ; ++i) {
        dens[CCTK_GFINDEX3D(i,j,k)] = -1e100
        Scon[CCTK_VECTGFINDEX3D(i,j,k,0)] = -1e100
        Scon[CCTK_VECTGFINDEX3D(i,j,k,1)] = -1e100
        Scon[CCTK_VECTGFINDEX3D(i,j,k,2)] = -1e100
        tau[CCTK_GFINDEX3D(i,j,k)] = -1e100
        if(do_mhd) {
          Bcons[CCTK_VECTGFINDEX3D(i,j,k,0)] = -1e100
          Bcons[CCTK_VECTGFINDEX3D(i,j,k,1)] = -1e100
          Bcons[CCTK_VECTGFINDEX3D(i,j,k,2)] = -1e100
        }
        if(do_clean_divergene)
          psidc[CCTK_GFINDEX3D(i,j,k)] = -1e100
      }
    }
  }
  for(i = 0 ; i < GRHydro_stencil*mol_substep ; ++i) {
    for(k = 0 ; k < cctk_ash[2] ; ++k) {
      for(j = 0 ; j < cctk_ash[1] ; ++j) {
        dens[CCTK_GFINDEX3D(i,j,k)] = -1e100
        Scon[CCTK_VECTGFINDEX3D(i,j,k,0)] = -1e100
        Scon[CCTK_VECTGFINDEX3D(i,j,k,1)] = -1e100
        Scon[CCTK_VECTGFINDEX3D(i,j,k,2)] = -1e100
        tau[CCTK_GFINDEX3D(i,j,k)] = -1e100
        if(do_mhd) {
          Bcons[CCTK_VECTGFINDEX3D(i,j,k,0)] = -1e100
          Bcons[CCTK_VECTGFINDEX3D(i,j,k,1)] = -1e100
          Bcons[CCTK_VECTGFINDEX3D(i,j,k,2)] = -1e100
        }
        if(do_clean_divergene)
          psidc[CCTK_GFINDEX3D(i,j,k)] = -1e100
      }
    }
  }
  for(i = cctk_ash[0]-GRHydro_stencil*mol_substep ; i < cctk_ash[0] ; ++i) {
    for(k = 0 ; k < cctk_ash[2] ; ++k) {
      for(j = 0 ; j < cctk_ash[1] ; ++j) {
        dens[CCTK_GFINDEX3D(i,j,k)] = -1e100
        Scon[CCTK_VECTGFINDEX3D(i,j,k,0)] = -1e100
        Scon[CCTK_VECTGFINDEX3D(i,j,k,1)] = -1e100
        Scon[CCTK_VECTGFINDEX3D(i,j,k,2)] = -1e100
        tau[CCTK_GFINDEX3D(i,j,k)] = -1e100
        if(do_mhd) {
          Bcons[CCTK_VECTGFINDEX3D(i,j,k,0)] = -1e100
          Bcons[CCTK_VECTGFINDEX3D(i,j,k,1)] = -1e100
          Bcons[CCTK_VECTGFINDEX3D(i,j,k,2)] = -1e100
        }
        if(do_clean_divergene)
          psidc[CCTK_GFINDEX3D(i,j,k)] = -1e100
      }
    }
  }
  for(j = 0 ; j < GRHydro_stencil*mol_substep ; ++j) {
    for(i = 0 ; i < cctk_ash[0] ; ++i) {
      for(k = 0 ; k < cctk_ash[2] ; ++k) {
        dens[CCTK_GFINDEX3D(i,j,k)] = -1e100
        Scon[CCTK_VECTGFINDEX3D(i,j,k,0)] = -1e100
        Scon[CCTK_VECTGFINDEX3D(i,j,k,1)] = -1e100
        Scon[CCTK_VECTGFINDEX3D(i,j,k,2)] = -1e100
        tau[CCTK_GFINDEX3D(i,j,k)] = -1e100
        if(do_mhd) {
          Bcons[CCTK_VECTGFINDEX3D(i,j,k,0)] = -1e100
          Bcons[CCTK_VECTGFINDEX3D(i,j,k,1)] = -1e100
          Bcons[CCTK_VECTGFINDEX3D(i,j,k,2)] = -1e100
        }
        if(do_clean_divergene)
          psidc[CCTK_GFINDEX3D(i,j,k)] = -1e100
      }
    }
  }
  for(j = cctk_ash[1]-GRHydro_stencil*mol_substep ; j < cctk_ash[1] ; ++j) {
    for(i = 0 ; i < cctk_ash[0] ; ++i) {
      for(k = 0 ; k < cctk_ash[2] ; ++k) {
        dens[CCTK_GFINDEX3D(i,j,k)] = -1e100
        Scon[CCTK_VECTGFINDEX3D(i,j,k,0)] = -1e100
        Scon[CCTK_VECTGFINDEX3D(i,j,k,1)] = -1e100
        Scon[CCTK_VECTGFINDEX3D(i,j,k,2)] = -1e100
        tau[CCTK_GFINDEX3D(i,j,k)] = -1e100
        if(do_mhd) {
          Bcons[CCTK_VECTGFINDEX3D(i,j,k,0)] = -1e100
          Bcons[CCTK_VECTGFINDEX3D(i,j,k,1)] = -1e100
          Bcons[CCTK_VECTGFINDEX3D(i,j,k,2)] = -1e100
        }
        if(do_clean_divergene)
          psidc[CCTK_GFINDEX3D(i,j,k)] = -1e100
      }
    }
  }
#endif
}

extern "C" void SourceTerms(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  void (*sources)(CCTK_ARGUMENTS) = NULL;
  const bool do_Avec = CCTK_EQUALS(Bvec_evolution_method, "GRHydro_Avec");
#define T(e,t) (!!(e)==(t))
#define T4(t1,t2,t3,t4)             \
        (T(*evolve_MHD,t1)      &&  \
         T(clean_divergence,t2) &&  \
         T(clean_divergence && divergence_flux,t3)  &&  \
         T(do_Avec,t4))
#define SET_SOURCES(df,t1,t2,t3,t4) \
        if(T4(t1,t2,t3,t4))         \
          sources=SourceTerms_LL<df,t1,t2,t3,t4>;
  enum do_t {noMHD=0,MHD=1,noDC=0,DC=1,noDCflux=0,DCflux=1,noAVec=0,AVec=1};

  // order of arguments:
  // do_MHD do_clean_divergence do_divergence_flux do_Avec
  if(spatial_order == 4) {
    // no MHD
    SET_SOURCES(alldiff4,noMHD,noDC,noDCflux,noAVec)

    // not all combinations make sense for MHD, we list the ones that do
    SET_SOURCES(alldiff4,  MHD,noDC,noDCflux,noAVec)
    SET_SOURCES(alldiff4,  MHD,  DC,noDCflux,noAVec)
    SET_SOURCES(alldiff4,  MHD,  DC,  DCflux,noAVec)
    SET_SOURCES(alldiff4,  MHD,noDC,noDCflux,  AVec)
  } else if(spatial_order == 2) {
    // no MHD
    SET_SOURCES(alldiff2,noMHD,noDC,noDCflux,noAVec)

    // not all combinations make sense for MHD, we list the ones that do
    SET_SOURCES(alldiff2,  MHD,noDC,noDCflux,noAVec)
    SET_SOURCES(alldiff2,  MHD,  DC,noDCflux,noAVec)
    SET_SOURCES(alldiff2,  MHD,  DC,  DCflux,noAVec)
    SET_SOURCES(alldiff2,  MHD,noDC,noDCflux,  AVec)
  }

  if(sources == NULL) {
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING, 
                "No source terms can be computed for combination "
                "do_mhd = %d "
                "do_clean_divergence = %d "
                "do_divergence_flux = %d "
                "do_Avec = %d spatial_order = %d",
                int(*evolve_MHD),
                int(clean_divergence),
                int(divergence_flux),
                int(do_Avec), int(spatial_order));
  }

  // do some non-time critial operations here to reduce the number of templates 
  if(transport_constraints) {
    const int N = cctk_ash[0]*cctk_ash[1]*cctk_ash[2];
    // not much need to OpenMP parallelize this since it is pure memory access
    // and memset is highly optimized
    memset(Evec, 0, sizeof(Evec[0])*3*N);
    // I think IEEE requires an all zero byte double to be 0.
    assert(Evec[0] == 0.);
  }
  if(*evolve_Y_e) {
    const int N = cctk_ash[0]*cctk_ash[1]*cctk_ash[2];
    // not much need to OpenMP parallelize this since it is pure memory access
    // and memset is highly optimized
    memset(Y_e_con_rhs, 0, sizeof(Y_e_con_rhs[0])*N);
    // I think IEEE requires an all zero byte double to be 0.
    assert(Y_e_con_rhs[0] == 0.);
  }
  if(number_of_tracers>0) {
    const int N = cctk_ash[0]*cctk_ash[1]*cctk_ash[2];
    // not much need to OpenMP parallelize this since it is pure memory access
    // and memset is highly optimized
    memset(cons_tracerrhs, 0, sizeof(cons_tracerrhs[0])*number_of_tracers*N);
    // I think IEEE requires an all zero byte double to be 0.
    assert(cons_tracerrhs[0] == 0.);
  }
  if(*evolve_entropy) {
    const int N = cctk_ash[0]*cctk_ash[1]*cctk_ash[2];
    // not much need to OpenMP parallelize this since it is pure memory access
    // and memset is highly optimized
    memset(entropyrhs, 0, sizeof(entropyrhs[0])*N);
    // I think IEEE requires an all zero byte double to be 0.
    assert(entropyrhs[0] == 0.);
  }
  if(track_divB) {
    const int N = cctk_ash[0]*cctk_ash[1]*cctk_ash[2];
    // not much need to OpenMP parallelize this since it is pure memory access
    // and memset is highly optimized
    memset(divB, 0, sizeof(divB[0])*N);
    // I think IEEE requires an all zero byte double to be 0.
    assert(divB[0] == 0.);
  }

  sources(CCTK_PASS_CTOC);
}


