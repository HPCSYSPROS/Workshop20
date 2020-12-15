#ifndef _GRHYDRO_RECONSTRUCT_DRV_IMPL_H
#define _GRHYDRO_RECONSTRUCT_DRV_IMPL_H

#include <cmath>
#include <vector>
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

#include "SpaceMask.h"

#include "GRHydro_Reconstruct_drv_cxx.hh"

#define velx (&vup[0])
#define vely (&vup[N])
#define velz (&vup[2*N])

#define Bvecx (&Bprim[0])
#define Bvecy (&Bprim[N])
#define Bvecz (&Bprim[2*N])


using namespace std;

extern "C" CCTK_INT GRHydro_UseGeneralCoordinates(const cGH * cctkGH);

/// transform back to vel if Wv was used
template <int dir>
static inline void undo_Wv(const int nx,
                   CCTK_REAL* const velxminus,
                   CCTK_REAL* const velyminus,
                   CCTK_REAL* const velzminus,
                   CCTK_REAL* const velxplus,
                   CCTK_REAL* const velyplus,
                   CCTK_REAL* const velzplus,
                   const CCTK_REAL* const gxx,
                   const CCTK_REAL* const gxy,
                   const CCTK_REAL* const gxz,
                   const CCTK_REAL* const gyy,
                   const CCTK_REAL* const gyz,
                   const CCTK_REAL* const gzz,
                   const cGH* const cctkGH,
                   const int j, const int k)
{
   DECLARE_CCTK_PARAMETERS;

   // TODO: check index bounds!!
   for (int i = GRHydro_stencil-1; i < nx - GRHydro_stencil+1; ++i) {
   
      const int ijk = dir == 0 ? CCTK_GFINDEX3D(cctkGH, i, j, k) : dir == 1 ? CCTK_GFINDEX3D(cctkGH, j, i, k) : CCTK_GFINDEX3D(cctkGH, j, k, i);
   
      // divide out the Loretnz factor obtained from w_lorentz =
      // sqrt(1+g_{ij} w^i w^j) for both the
      // plus and minus quantities this should by construction ensure
      // that any Lorentz factor calculated from them later on is
      // physical (ie. > 1.d0)
      {
      const int ijk_m = dir == 0 ? CCTK_GFINDEX3D(cctkGH, i-1, j, k) : dir == 1 ? CCTK_GFINDEX3D(cctkGH, j, i-1, k) : CCTK_GFINDEX3D(cctkGH, j, k, i-1);
      
      const CCTK_REAL agxx = 0.5*( gxx[ijk] + gxx[ijk_m] );
      const CCTK_REAL agxy = 0.5*( gxy[ijk] + gxy[ijk_m] );
      const CCTK_REAL agxz = 0.5*( gxz[ijk] + gxz[ijk_m] );
      const CCTK_REAL agyy = 0.5*( gyy[ijk] + gyy[ijk_m] );
      const CCTK_REAL agyz = 0.5*( gyz[ijk] + gyz[ijk_m] );
      const CCTK_REAL agzz = 0.5*( gzz[ijk] + gzz[ijk_m] );
      const CCTK_REAL w = sqrt( 1.0 + agxx*velxminus[ijk]*velxminus[ijk] +
                                      agyy*velyminus[ijk]*velyminus[ijk] +
                                      agzz*velzminus[ijk]*velzminus[ijk] +
                                  2.0*agxy*velxminus[ijk]*velyminus[ijk] +
                                  2.0*agxz*velxminus[ijk]*velzminus[ijk] +
                                  2.0*agyz*velyminus[ijk]*velzminus[ijk] );
      velxminus[ijk] = velxminus[ijk]/w;
      velyminus[ijk] = velyminus[ijk]/w;
      velzminus[ijk] = velzminus[ijk]/w;
      }

      {
      const int ijk_p = dir == 0 ? CCTK_GFINDEX3D(cctkGH, i+1, j, k) : dir == 1 ? CCTK_GFINDEX3D(cctkGH, j, i+1, k) : CCTK_GFINDEX3D(cctkGH, j, k, i+1);
      
      const CCTK_REAL agxx = 0.5*( gxx[ijk] + gxx[ijk_p] );
      const CCTK_REAL agxy = 0.5*( gxy[ijk] + gxy[ijk_p] );
      const CCTK_REAL agxz = 0.5*( gxz[ijk] + gxz[ijk_p] );
      const CCTK_REAL agyy = 0.5*( gyy[ijk] + gyy[ijk_p] );
      const CCTK_REAL agyz = 0.5*( gyz[ijk] + gyz[ijk_p] );
      const CCTK_REAL agzz = 0.5*( gzz[ijk] + gzz[ijk_p] );
      const CCTK_REAL w = sqrt( 1.0 + agxx*velxplus[ijk]*velxplus[ijk] +
                                      agyy*velyplus[ijk]*velyplus[ijk] +
                                      agzz*velzplus[ijk]*velzplus[ijk] +
                                  2.0*agxy*velxplus[ijk]*velyplus[ijk] +
                                  2.0*agxz*velxplus[ijk]*velzplus[ijk] +
                                  2.0*agyz*velyplus[ijk]*velzplus[ijk] );
     velxplus[ijk] = velxplus[ijk]/w;
     velyplus[ijk] = velyplus[ijk]/w;
     velzplus[ijk] = velzplus[ijk]/w;
     }
   }
}


/*
  Cases that must be considered:
  * basic hydro
  * hydro + temperature + ye
  * hydro + ye
  * basic mhd
  * mhd + temperature + ye 
  * mhd + ye 
 */

// TODO: remove all templatization from this reoutine since they are only used
// in the 2d loops and bloat the source code by a factor of 12
// TODO: measure if 1d slicing slows things down since this would buy a factor
// of 3 in source code size again
template <bool do_mhd,
          bool do_Ye,
          bool do_temp,
          bool do_reconstruct_Wv,
          bool do_clean_divergence,
          class RECONSTRUCT>                   // the reconstruction operator
void GRHydro_Reconstruct_drv_cxx(cGH const * const restrict cctkGH)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL * restrict vup;
  const CCTK_REAL * restrict Bprim;
  const CCTK_REAL * restrict g11;
  const CCTK_REAL * restrict g12;
  const CCTK_REAL * restrict g13;
  const CCTK_REAL * restrict g22;
  const CCTK_REAL * restrict g23;
  const CCTK_REAL * restrict g33;

  //Multipatch related pointers
  if(GRHydro_UseGeneralCoordinates(cctkGH)) {
    vup=lvel;
    g11=gaa; g12=gab; g13=gac;
    g22=gbb; g23=gbc;
    g33=gcc;
    Bprim=lBvec;
  } else {
    vup=vel;
    g11=gxx; g12=gxy; g13=gxz;
    g22=gyy; g23=gyz;
    g33=gzz;
    Bprim=Bvec;
  }

  const int nx=cctk_lsh[0];
  const int ny=cctk_lsh[1];
  const int nz=cctk_lsh[2];

  const int N = nx*ny*nz;

  int type_bitsx,type_bitsy,type_bitsz;
  // UNUSED:   int trivialx,trivialy,trivialz;
  int not_trivialx,not_trivialy,not_trivialz;

  type_bitsx = SpaceMask_GetTypeBits("Hydro_RiemannProblemX");
  //  trivialx = SpaceMask_GetStateBits("Hydro_RiemannProblemX","trivial");
  not_trivialx = SpaceMask_GetStateBits("Hydro_RiemannProblemX","not_trivial");

  type_bitsy = SpaceMask_GetTypeBits("Hydro_RiemannProblemY");
  //  trivialy = SpaceMask_GetStateBits("Hydro_RiemannProblemY","trivial");
  not_trivialy = SpaceMask_GetStateBits("Hydro_RiemannProblemY","not_trivial");

  type_bitsz = SpaceMask_GetTypeBits("Hydro_RiemannProblemZ");
  //  trivialz = SpaceMask_GetStateBits("Hydro_RiemannProblemZ","trivial");
  not_trivialz = SpaceMask_GetStateBits("Hydro_RiemannProblemZ","not_trivial");
  
  
  vector<CCTK_REAL> Wvx(0);
  vector<CCTK_REAL> Wvy(0);
  vector<CCTK_REAL> Wvz(0);
  // allocate temp memory for Wv reconstruction
  // TODO: measure impact of C++ intializing the vector to zero and write
  // custom wrapper around new/delete if required
  // template<class T> class pod_vector {
  //   pod_vector() {data = NULL;};
  //   ~pod_vecotr() {delete[] data;};
  //   T& operator[](size_t i) {return data[i];};
  //   void allocate(size_t sz) {data = new T[sz];};
  // private:
  //   pod_vector(const pod_vector&);
  //   pod_vector& operator=(const pod_vector&);
  // };
  if (do_reconstruct_Wv) {
     Wvx.resize(N, 0);
     Wvy.resize(N, 0);
     Wvz.resize(N, 0);
     #pragma omp parallel for
     for (int i=0; i < N; ++i) {
        Wvx[i] = w_lorentz[i]*velx[i];
        Wvy[i] = w_lorentz[i]*vely[i];
        Wvz[i] = w_lorentz[i]*velz[i];
     }
  }

// Reconstruction starts:
  if (flux_direction[0] == 1) {
    #pragma omp parallel for
    for(int k = GRHydro_stencil-1; k < nz - GRHydro_stencil+1+transport_constraints; ++k) {
      for(int j = GRHydro_stencil-1; j < ny - GRHydro_stencil+1+transport_constraints; ++j) {

        RECONSTRUCT::template apply<0>(nx, rho, rhominus, rhoplus, cctkGH, j, k);
        RECONSTRUCT::template apply<0>(nx, eps, epsminus, epsplus, cctkGH, j, k);

        if (do_reconstruct_Wv) {
           RECONSTRUCT::template apply<0>(nx, &Wvx[0], velxminus, velxplus, cctkGH, j, k);
           RECONSTRUCT::template apply<0>(nx, &Wvy[0], velyminus, velyplus, cctkGH, j, k);
           RECONSTRUCT::template apply<0>(nx, &Wvz[0], velzminus, velzplus, cctkGH, j, k);
           undo_Wv<0>(nx, velxminus, velyminus, velzminus,
                          velxplus, velyplus, velzplus,
                          g11, g12, g13, g22, g23, g33,
                          cctkGH, j, k);
        } else {
           RECONSTRUCT::template apply<0>(nx, velx, velxminus, velxplus, cctkGH, j, k);
           RECONSTRUCT::template apply<0>(nx, vely, velyminus, velyplus, cctkGH, j, k);
           RECONSTRUCT::template apply<0>(nx, velz, velzminus, velzplus, cctkGH, j, k);
        }
        
        if (do_Ye)
           RECONSTRUCT::template apply<0>(nx, Y_e, Y_e_minus, Y_e_plus, cctkGH, j, k);
        if (do_temp)
           RECONSTRUCT::template apply<0>(nx, temperature, tempminus, tempplus, cctkGH, j, k);
        if (do_mhd) {
           RECONSTRUCT::template apply<0>(nx, Bvecx, Bvecxminus, Bvecxplus, cctkGH, j, k);
           RECONSTRUCT::template apply<0>(nx, Bvecy, Bvecyminus, Bvecyplus, cctkGH, j, k);
           RECONSTRUCT::template apply<0>(nx, Bvecz, Bveczminus, Bveczplus, cctkGH, j, k);
           if (do_clean_divergence)
              RECONSTRUCT::template apply<0>(nx, psidc, psidcminus, psidcplus, cctkGH, j, k);
        }

        for (int i=0; i<nx; ++i) {
          SpaceMask_SetStateBits(space_mask, i+j*nx+k*nx*ny, type_bitsx, not_trivialx);
        }
      }
    }  

  } else if (flux_direction[0]==2) {
    //Make sure to doublecheck the bounds here!
    #pragma omp parallel for
    for(int k = GRHydro_stencil-1; k < nz - GRHydro_stencil+1+transport_constraints; ++k) {
      for(int j = GRHydro_stencil-1; j < nx - GRHydro_stencil+1+transport_constraints; ++j) {
        
        RECONSTRUCT::template apply<1>(ny, rho, rhominus, rhoplus, cctkGH, j, k);
        RECONSTRUCT::template apply<1>(ny, eps, epsminus, epsplus, cctkGH, j, k);

        if (do_reconstruct_Wv) {
           RECONSTRUCT::template apply<1>(ny, &Wvx[0], velxminus, velxplus, cctkGH, j, k);
           RECONSTRUCT::template apply<1>(ny, &Wvy[0], velyminus, velyplus, cctkGH, j, k);
           RECONSTRUCT::template apply<1>(ny, &Wvz[0], velzminus, velzplus, cctkGH, j, k);
           undo_Wv<1>(ny, velyminus, velzminus, velxminus,
                          velyplus, velzplus, velxplus,
                          g22, g23, g12, g33, g13, g11,
                          cctkGH, j, k);
        } else {
           RECONSTRUCT::template apply<1>(ny, velx, velxminus, velxplus, cctkGH, j, k);
           RECONSTRUCT::template apply<1>(ny, vely, velyminus, velyplus, cctkGH, j, k);
           RECONSTRUCT::template apply<1>(ny, velz, velzminus, velzplus, cctkGH, j, k);
        }
        
        if (do_Ye)
           RECONSTRUCT::template apply<1>(ny, Y_e, Y_e_minus, Y_e_plus, cctkGH, j, k);
        if (do_temp)
           RECONSTRUCT::template apply<1>(ny, temperature, tempminus, tempplus, cctkGH, j, k);
        if (do_mhd) {
           RECONSTRUCT::template apply<1>(ny, Bvecx, Bvecxminus, Bvecxplus, cctkGH, j, k);
           RECONSTRUCT::template apply<1>(ny, Bvecy, Bvecyminus, Bvecyplus, cctkGH, j, k);
           RECONSTRUCT::template apply<1>(ny, Bvecz, Bveczminus, Bveczplus, cctkGH, j, k);
           if (do_clean_divergence)
              RECONSTRUCT::template apply<1>(ny, psidc, psidcminus, psidcplus, cctkGH, j, k);
        }

        for (int i=0; i<ny; ++i) {
          SpaceMask_SetStateBits(space_mask, j+i*nx+k*nx*ny, type_bitsy, not_trivialy);
        }
        
      }
    }
    
  } else if (flux_direction[0]==3) {
    //Make sure to doublecheck the bounds here!
    #pragma omp parallel for
    for(int k = GRHydro_stencil-1; k < ny - GRHydro_stencil+1+transport_constraints; ++k) {
      for(int j = GRHydro_stencil-1; j < nx - GRHydro_stencil+1+transport_constraints; ++j) {
        
        RECONSTRUCT::template apply<2>(nz, rho, rhominus, rhoplus, cctkGH, j, k);
        RECONSTRUCT::template apply<2>(nz, eps, epsminus, epsplus, cctkGH, j, k);

        if (do_reconstruct_Wv) {
           RECONSTRUCT::template apply<2>(nz, &Wvx[0], velxminus, velxplus, cctkGH, j, k);
           RECONSTRUCT::template apply<2>(nz, &Wvy[0], velyminus, velyplus, cctkGH, j, k);
           RECONSTRUCT::template apply<2>(nz, &Wvz[0], velzminus, velzplus, cctkGH, j, k);
           undo_Wv<2>(nz, velzminus, velxminus, velyminus,
                          velzplus, velxplus, velyplus,
                          g33, g13, g23, g11, g12, g22,
                          cctkGH, j, k);
        } else {
           RECONSTRUCT::template apply<2>(nz, velx, velxminus, velxplus, cctkGH, j, k);
           RECONSTRUCT::template apply<2>(nz, vely, velyminus, velyplus, cctkGH, j, k);
           RECONSTRUCT::template apply<2>(nz, velz, velzminus, velzplus, cctkGH, j, k);
        }
        
        if (do_Ye)
           RECONSTRUCT::template apply<2>(nz, Y_e, Y_e_minus, Y_e_plus, cctkGH, j, k);
        if (do_temp)
           RECONSTRUCT::template apply<2>(nz, temperature, tempminus, tempplus, cctkGH, j, k);
        if (do_mhd) {
           RECONSTRUCT::template apply<2>(nz, Bvecx, Bvecxminus, Bvecxplus, cctkGH, j, k);
           RECONSTRUCT::template apply<2>(nz, Bvecy, Bvecyminus, Bvecyplus, cctkGH, j, k);
           RECONSTRUCT::template apply<2>(nz, Bvecz, Bveczminus, Bveczplus, cctkGH, j, k);
           if (do_clean_divergence)
              RECONSTRUCT::template apply<2>(nz, psidc, psidcminus, psidcplus, cctkGH, j, k);
        }
        
        for (int i=0; i<nz; ++i) {
          SpaceMask_SetStateBits(space_mask, j+k*nx+i*nx*ny, type_bitsz, not_trivialz);
        }
        
      }
    }
    
  }
}

#endif // _GRHYDRO_RECONSTRUCT_DRV_IMPL_H
