 /*@@
   @file      GRHydro_Reconstruct.cc
   @date      Thu Dec 2013
   @author    Christian Reisswig
   @desc 
   Wrapper routine to perform the reconstruction.using tmeplate metaprogramming
   @enddesc 
 @@*/

#include <iostream>
#include <cassert>
#include <cstring>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "SpaceMask.h"

#include "GRHydro_Functions.h"

#include "GRHydro_Reconstruct_drv_cxx.hh"
#include "GRHydro_WENOReconstruct.hh"
#include "GRHydro_TrivialReconstruct.hh"
#include "GRHydro_TVDReconstruct.hh"
#include "GRHydro_MP5Reconstruct.hh"

#define velx (&vup[0])
#define vely (&vup[N])
#define velz (&vup[2*N])

#define Bvecx (&Bprim[0])
#define Bvecy (&Bprim[N])
#define Bvecz (&Bprim[2*N])

/*
  External routines
*/
extern "C"
CCTK_INT GRHydro_UseGeneralCoordinates(const cGH * cctkGH);

/*
  This is needed right now because we want to also call this from Fotran
*/
extern "C"
void CCTK_FCALL CCTK_FNAME(GRHydro_PPMReconstruct_drv_opt)(cGH const * const restrict & restrict cctkGH);

extern "C"
void CCTK_FCALL CCTK_FNAME(primitive2conservativeM)(cGH const * const restrict & restrict cctkGH);


/*
   Undoing Wv pointwise.
   This is necessary for atmo and excision treatment.
*/
inline void undo_Wv_ptw(const CCTK_REAL* const restrict gxx,
                   const CCTK_REAL* const restrict gxy,
                   const CCTK_REAL* const restrict gxz,
                   const CCTK_REAL* const restrict gyy,
                   const CCTK_REAL* const restrict gyz,
                   const CCTK_REAL* const restrict gzz,
                   CCTK_REAL* const restrict vx,
                   CCTK_REAL* const restrict vy,
                   CCTK_REAL* const restrict vz,
                   const int ijk, const int ijk_pm)
{
  const CCTK_REAL agxx = 0.5*( gxx[ijk] + gxx[ijk_pm] );
  const CCTK_REAL agxy = 0.5*( gxy[ijk] + gxy[ijk_pm] );
  const CCTK_REAL agxz = 0.5*( gxz[ijk] + gxz[ijk_pm] );
  const CCTK_REAL agyy = 0.5*( gyy[ijk] + gyy[ijk_pm] );
  const CCTK_REAL agyz = 0.5*( gyz[ijk] + gyz[ijk_pm] );
  const CCTK_REAL agzz = 0.5*( gzz[ijk] + gzz[ijk_pm] );
  const CCTK_REAL w = sqrt( 1.0 + agxx*vx[ijk]*vx[ijk] +
                                  agyy*vy[ijk]*vy[ijk] +
                                  agzz*vz[ijk]*vz[ijk] +
                              2.0*agxy*vx[ijk]*vy[ijk] +
                              2.0*agxz*vx[ijk]*vz[ijk] +
                              2.0*agyz*vy[ijk]*vz[ijk] );
                 
  vx[ijk] = vx[ijk]/w;
  vy[ijk] = vy[ijk]/w;
  vz[ijk] = vz[ijk]/w;
}


/**
   This class selects a reconstruction variant.
*/
template <class RECONSTRUCT>
struct reconstruct
{
   static inline void select(const bool do_mhd,
                             const bool do_Ye,
                             const bool do_temp,
                             const bool reconstruct_Wv,
                             const bool reconstruct_pressure,
                             const bool clean_divergence,
                             cGH const * const restrict & restrict cctkGH)
   {
       GRHydro_Reconstruct_drv_cxx <RECONSTRUCT> (cctkGH, do_mhd, do_Ye,
                                                  do_temp, reconstruct_Wv,
                                                  reconstruct_pressure,
                                                  clean_divergence);
   }
};


/**
   This is a new template metaprogramming way of calling various reconstruction methods.
   It tries to avoid duplicated code as much as possible
*/

extern "C" void Reconstruction_cxx(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const bool do_mhd = CCTK_EQUALS(Bvec_evolution_method, "GRHydro");
  const bool do_Ye = CCTK_EQUALS(Y_e_evolution_method, "GRHydro");
  const bool do_temp = CCTK_EQUALS(temperature_evolution_method, "GRHydro");
  const bool do_clean_divergence = clean_divergence;

  /*
     Set up multipatch stuff
  */

  const CCTK_REAL * restrict vup;
  const CCTK_REAL * restrict Bprim;

  //Multipatch related pointers
  if(GRHydro_UseGeneralCoordinates(cctkGH)) {
    vup=lvel;
    Bprim=lBvec;
  } else {
    vup=vel;
    Bprim=Bvec;
  }
  
  const int nx=cctk_lsh[0];
  const int ny=cctk_lsh[1];
  const int nz=cctk_lsh[2];

  const int N = nx*ny*nz;

  int type_bits;
  int not_trivial;
  int trivial;

  if (*flux_direction == 1) {
     type_bits = SpaceMask_GetTypeBits("Hydro_RiemannProblemX");
     not_trivial = SpaceMask_GetStateBits("Hydro_RiemannProblemX","not_trivial");
     trivial     = SpaceMask_GetStateBits("Hydro_RiemannProblemX","trivial");
  } else if (*flux_direction == 2) {
     type_bits = SpaceMask_GetTypeBits("Hydro_RiemannProblemY");
     not_trivial = SpaceMask_GetStateBits("Hydro_RiemannProblemY","not_trivial");
     trivial     = SpaceMask_GetStateBits("Hydro_RiemannProblemY","trivial");
  } else if (*flux_direction == 3) {
     type_bits = SpaceMask_GetTypeBits("Hydro_RiemannProblemZ");
     not_trivial = SpaceMask_GetStateBits("Hydro_RiemannProblemZ","not_trivial");
     trivial     = SpaceMask_GetStateBits("Hydro_RiemannProblemZ","trivial");
  }

  CCTK_REAL local_min_tracer = 1e42;
  if (evolve_tracer) {
    if (use_min_tracer)
      local_min_tracer = min_tracer;
    else
      local_min_tracer = 0.0;
  }

  /*
     First initialize all variables by zeroth order trivial reconstruction!
  */


  memcpy(rhominus,rho,sizeof(CCTK_REAL)*N);
  memcpy(velxminus,vel+0*N,sizeof(CCTK_REAL)*N);
  memcpy(velyminus,vel+1*N,sizeof(CCTK_REAL)*N);
  memcpy(velzminus,vel+2*N,sizeof(CCTK_REAL)*N);
  if(reconstruct_pressure)
    memcpy(pressminus,press,sizeof(CCTK_REAL)*N);
  else
    memcpy(epsminus,eps,sizeof(CCTK_REAL)*N);
  if(do_Ye)
    memcpy(Y_e_minus,Y_e,sizeof(CCTK_REAL)*N);
  if(do_temp)
    memcpy(tempminus,temperature,sizeof(CCTK_REAL)*N);

  memcpy(rhoplus,rho,sizeof(CCTK_REAL)*N);
  memcpy(velxplus,vel+0*N,sizeof(CCTK_REAL)*N);
  memcpy(velyplus,vel+1*N,sizeof(CCTK_REAL)*N);
  memcpy(velzplus,vel+2*N,sizeof(CCTK_REAL)*N);
  if(reconstruct_pressure)
    memcpy(pressplus,press,sizeof(CCTK_REAL)*N);
  else
    memcpy(epsplus,eps,sizeof(CCTK_REAL)*N);
  if(do_Ye)
    memcpy(Y_e_plus,Y_e,sizeof(CCTK_REAL)*N);
  if(do_temp)
    memcpy(tempplus,temperature,sizeof(CCTK_REAL)*N);

  /*
     Now do the real reconstruction!
  */

  if (CCTK_EQUALS(recon_method,"tvd")) {
    // this handles MHD and non-MHD

    if (CCTK_EQUALS(tvd_limiter, "minmod")) {
       
       reconstruct<GRHydro_TVDReconstruct1d<tvd::minmod> >::select(do_mhd, do_Ye, do_temp, reconstruct_Wv, reconstruct_pressure, clean_divergence, cctkGH);
       
    } else if (CCTK_EQUALS(tvd_limiter, "vanleerMC2")) {
       
       reconstruct<GRHydro_TVDReconstruct1d<tvd::mc2> >::select(do_mhd, do_Ye, do_temp, reconstruct_Wv, reconstruct_pressure, clean_divergence, cctkGH);
       
    } else if (CCTK_EQUALS(tvd_limiter, "superbee")) {
       
       reconstruct<GRHydro_TVDReconstruct1d<tvd::superbee> >::select(do_mhd, do_Ye, do_temp, reconstruct_Wv, reconstruct_pressure, clean_divergence, cctkGH);
    
    } else {
       CCTK_ERROR("TVD limiter not recognized!");
    }

  } else if (CCTK_EQUALS(recon_method,"ppm")) {
    // PPM is special. It cannot reconstruct individual functions.
    // Rather, all function are reconstructed at once!
    // Hence we have a specialized calling function.
    
    CCTK_FNAME(GRHydro_PPMReconstruct_drv_opt)(cctkGH);

  } else if (CCTK_EQUALS(recon_method,"eno")) {
    // this handles MHD and non-MHD

    // TODO!!
    CCTK_ERROR("Sorry. Eno reconstruction not yet implemented!");

  } else if (CCTK_EQUALS(recon_method,"weno")) {
    // this handles MHD and non-MHD

    if (weno_adaptive_epsilon) {
       
       reconstruct<GRHydro_WENOReconstruct1d_cxx<false, true> >::select(do_mhd, do_Ye, do_temp, reconstruct_Wv, reconstruct_pressure, clean_divergence, cctkGH);
          
    } else {
       
       reconstruct<GRHydro_WENOReconstruct1d_cxx<false, false> >::select(do_mhd, do_Ye, do_temp, reconstruct_Wv, reconstruct_pressure, clean_divergence, cctkGH);
    
    }

  } else if (CCTK_EQUALS(recon_method, "weno-z")) {

       reconstruct<GRHydro_WENOReconstruct1d_cxx<true, false> >::select(do_mhd, do_Ye, do_temp, reconstruct_Wv, reconstruct_pressure, clean_divergence, cctkGH);

  } else if (CCTK_EQUALS(recon_method,"mp5")) {
    // this handles MHD and non-MHD

    if (mp5_adaptive_eps)
       reconstruct<GRHydro_MP5Reconstruct1d_cxx<true> >::select(do_mhd, do_Ye, do_temp, reconstruct_Wv, reconstruct_pressure, clean_divergence, cctkGH);
    else
       reconstruct<GRHydro_MP5Reconstruct1d_cxx<false> >::select(do_mhd, do_Ye, do_temp, reconstruct_Wv, reconstruct_pressure, clean_divergence, cctkGH);

  } else
    CCTK_ERROR("Reconstruction method not recognized!");


  


  /*
     Here goes the atmosphere treatment!
  */
  #pragma omp parallel for
  for (int k=GRHydro_stencil-1; k < cctk_lsh[2] - GRHydro_stencil+1+transport_constraints*(1-*zoffset); ++k)
     for (int j=GRHydro_stencil-1; j < cctk_lsh[1] - GRHydro_stencil+1+transport_constraints*(1-*yoffset); ++j)
        for (int i=GRHydro_stencil-1; i < cctk_lsh[0] - GRHydro_stencil+1+transport_constraints*(1-*xoffset); ++i) {
           const int ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
           const int ijk_m = CCTK_GFINDEX3D(cctkGH, i-*xoffset, j-*yoffset, k-*zoffset);
           const int ijk_p = CCTK_GFINDEX3D(cctkGH, i+*xoffset, j+*yoffset, k+*zoffset);
           
           // Whenever we modify a previously reconstructed value,
           // we need to change the trivial flag!
           bool change_riemann_trivial_flag = false;
           
           if (rhominus[ijk] < *GRHydro_rho_min || rhoplus[ijk] < *GRHydro_rho_min || (!do_temp && (epsplus[ijk] < 0.0 || epsminus[ijk] < 0.0)) )
           {
              rhominus[ijk] = rho[ijk];
              rhoplus[ijk]  = rho[ijk];
              
              epsminus[ijk] = eps[ijk];
              epsplus[ijk]  = eps[ijk];
              
              if (reconstruct_Wv) {
                 velxminus[ijk] = w_lorentz[ijk]*velx[ijk];
                 velxplus[ijk]  = w_lorentz[ijk]*velx[ijk];
                 velyminus[ijk] = w_lorentz[ijk]*vely[ijk];
                 velyplus[ijk]  = w_lorentz[ijk]*vely[ijk];
                 velzminus[ijk] = w_lorentz[ijk]*velz[ijk];
                 velzplus[ijk]  = w_lorentz[ijk]*velz[ijk];
                 
                 undo_Wv_ptw(gxx, gxy, gxz, gyy, gyz, gzz,
                             velxplus, velyplus, velzplus, ijk, ijk_p);
                 undo_Wv_ptw(gxx, gxy, gxz, gyy, gyz, gzz,
                             velxminus, velyminus, velzminus, ijk, ijk_m);
                 
              } else {
                 velxminus[ijk] = velx[ijk];
                 velxplus[ijk]  = velx[ijk];
                 velyminus[ijk] = vely[ijk];
                 velyplus[ijk]  = vely[ijk];
                 velzminus[ijk] = velz[ijk];
                 velzplus[ijk]  = velz[ijk];
              }
              
              if (do_Ye) {
                 Y_e_plus[ijk]  = Y_e[ijk];
                 Y_e_minus[ijk] = Y_e[ijk];
              }
              if (do_temp) {
                 tempminus[ijk] = temperature[ijk];
                 tempplus[ijk] = temperature[ijk];
              }
              if (do_mhd) {
                 Bvecxminus[ijk] = Bvecx[ijk];
                 Bvecxplus[ijk]  = Bvecx[ijk];
                 Bvecyminus[ijk] = Bvecy[ijk];
                 Bvecyplus[ijk]  = Bvecy[ijk];
                 Bveczminus[ijk] = Bvecz[ijk];
                 Bveczplus[ijk]  = Bvecz[ijk];
                 if (do_clean_divergence) {
                    psidcminus[ijk] = psidc[ijk];
                    psidcplus[ijk] = psidc[ijk];
                 }
              }
              if (evolve_tracer) {
                 for (int l=0; l < number_of_tracers; ++l) {
                    if (tracerminus[ijk + l*N] < local_min_tracer || tracerplus[ijk + l*N] < local_min_tracer) {
                       tracerminus[ijk + l*N] = tracer[ijk + l*N];
                       tracerplus[ijk + l*N]  = tracer[ijk + l*N];
                    }
                 }
              }
              
              change_riemann_trivial_flag = true;
           }
           
//           if (!do_temp) {
//               if (epsplus[ijk] < 0.0) {
//                  epsplus[ijk] = eps[ijk];
//                  change_riemann_trivial_flag = true;
//               }
//               if (epsminus[ijk] < 0.0) {
//                  epsminus[ijk] = eps[ijk];
//                  change_riemann_trivial_flag = true;
//               }
//               if (epsplus[ijk] < 0.0 || epsminus[ijk] < 0.0) {
//                  epsplus[ijk] = eps[ijk];
//                  epsminus[ijk] = eps[ijk];
//                  change_riemann_trivial_flag = true;
//               }
//               
//            }
           
           if (change_riemann_trivial_flag) {
              // Riemann problem might not be trivial anymore!
              SpaceMask_SetStateBits(space_mask, ijk_m, type_bits, not_trivial);
              SpaceMask_SetStateBits(space_mask, ijk,   type_bits, not_trivial);
           }
        }

        
  /*
     Here goes the excision treatment.
     If any of the points in the reconstruction stencil is excised,
     we switch to trivial reconstruction.
     If the current point itself is excised, we mark the Riemann problem as "trivial".
     If the current point is just one point away from an excision boundary, we use the current cell values
     to provide an "extrapolated" guess for the reconstructed state at the excision boundary
  */
  if (GRHydro_enable_internal_excision) {
     
      #pragma omp parallel for
      for (int k=GRHydro_stencil-1; k < cctk_lsh[2] - GRHydro_stencil+1+transport_constraints*(1-*zoffset); ++k)
         for (int j=GRHydro_stencil-1; j < cctk_lsh[1] - GRHydro_stencil+1+transport_constraints*(1-*yoffset); ++j)
            for (int i=GRHydro_stencil-1; i < cctk_lsh[0] - GRHydro_stencil+1+transport_constraints*(1-*xoffset); ++i) {
               
               const int ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
               const int ijk_m = CCTK_GFINDEX3D(cctkGH, i-*xoffset, j-*yoffset, k-*zoffset);
               const int ijk_p = CCTK_GFINDEX3D(cctkGH, i+*xoffset, j+*yoffset, k+*zoffset);
               
               bool excise = false;
               bool first_order = false;
               bool extrapolate = false;
               int extrap_ijk = 0;
               
               // Whenever we modify the previously reconstructed value,
               // we need to change the trivial flag!
               bool change_riemann_trivial_flag = false;
               
               // for an excised point in the reconstruction stencil radius of the current point, we need to switch down to first order
               //for (int l=-(GRHydro_stencil-1); l < GRHydro_stencil-1; ++l) {
               for (int l=-(GRHydro_stencil-1); l < GRHydro_stencil-1; ++l) {
                  const int ijk_l = CCTK_GFINDEX3D(cctkGH, i+*xoffset*l, j+*yoffset*l, k+*zoffset*l);
                  if (hydro_excision_mask[ijk_l]) {
                     first_order = true;
                     break;
                  }
               }
               
               // if the current point is excised, then we will mark the Riemann problem as trivial.
               if (hydro_excision_mask[ijk]) {
                  excise = true;
               
                  if (!hydro_excision_mask[ijk_m]) {
                     extrapolate = true;
                     extrap_ijk = ijk_m;
                  }
               
                  if (!hydro_excision_mask[ijk_p]) {
                     extrapolate = true;
                     extrap_ijk = ijk_p;
                  }
               }
               
               if (first_order && !extrapolate)
               {
                  rhominus[ijk] = rho[ijk];
                  rhoplus[ijk]  = rho[ijk];
                  
                  epsminus[ijk] = eps[ijk];
                  epsplus[ijk]  = eps[ijk];
                  
                  if (reconstruct_Wv) {
                     velxminus[ijk] = w_lorentz[ijk]*velx[ijk];
                     velxplus[ijk]  = w_lorentz[ijk]*velx[ijk];
                     velyminus[ijk] = w_lorentz[ijk]*vely[ijk];
                     velyplus[ijk]  = w_lorentz[ijk]*vely[ijk];
                     velzminus[ijk] = w_lorentz[ijk]*velz[ijk];
                     velzplus[ijk]  = w_lorentz[ijk]*velz[ijk];
                 
                     undo_Wv_ptw(gxx, gxy, gxz, gyy, gyz, gzz,
                             velxplus, velyplus, velzplus, ijk, ijk_p);
                     undo_Wv_ptw(gxx, gxy, gxz, gyy, gyz, gzz,
                             velxminus, velyminus, velzminus, ijk, ijk_m);
                  } else {
                     velxminus[ijk] = velx[ijk];
                     velxplus[ijk]  = velx[ijk];
                     velyminus[ijk] = vely[ijk];
                     velyplus[ijk]  = vely[ijk];
                     velzminus[ijk] = velz[ijk];
                     velzplus[ijk]  = velz[ijk];
                  }
                  
                  if (do_Ye) {
                     Y_e_plus[ijk]  = Y_e[ijk];
                     Y_e_minus[ijk] = Y_e[ijk];
                  }
                  if (do_temp) {
                    tempminus[ijk] = temperature[ijk];
                    tempplus[ijk] = temperature[ijk];
                  }
                  if (do_mhd) {
                    Bvecxminus[ijk] = Bvecx[ijk];
                    Bvecxplus[ijk]  = Bvecx[ijk];
                    Bvecyminus[ijk] = Bvecy[ijk];
                    Bvecyplus[ijk]  = Bvecy[ijk];
                    Bveczminus[ijk] = Bvecz[ijk];
                    Bveczplus[ijk]  = Bvecz[ijk];
                    if (do_clean_divergence) {
                       psidcminus[ijk] = psidc[ijk];
                       psidcplus[ijk] = psidc[ijk];
                    }
                  }
                  if (evolve_tracer) {
                     for (int l=0; l < number_of_tracers; ++l) {
                        tracerminus[ijk + l*N] = tracer[ijk + l*N];
                        tracerplus[ijk + l*N]  = tracer[ijk + l*N];
                     }
                  }

                  // Riemann problem might not be trivial anymore!
                  change_riemann_trivial_flag = true;
               } 
               else if (extrapolate)
               {
                  rhominus[ijk] = rho[extrap_ijk];
                  rhoplus[ijk]  = rho[extrap_ijk];
                  
                  epsminus[ijk] = eps[extrap_ijk];
                  epsplus[ijk]  = eps[extrap_ijk];
                  
                  if (reconstruct_Wv) {
                     velxminus[ijk] = w_lorentz[extrap_ijk]*velx[extrap_ijk];
                     velxplus[ijk]  = w_lorentz[extrap_ijk]*velx[extrap_ijk];
                     velyminus[ijk] = w_lorentz[extrap_ijk]*vely[extrap_ijk];
                     velyplus[ijk]  = w_lorentz[extrap_ijk]*vely[extrap_ijk];
                     velzminus[ijk] = w_lorentz[extrap_ijk]*velz[extrap_ijk];
                     velzplus[ijk]  = w_lorentz[extrap_ijk]*velz[extrap_ijk];
                 
                     undo_Wv_ptw(gxx, gxy, gxz, gyy, gyz, gzz,
                             velxplus, velyplus, velzplus, ijk, ijk_p);
                     undo_Wv_ptw(gxx, gxy, gxz, gyy, gyz, gzz,
                             velxminus, velyminus, velzminus, ijk, ijk_m);
                  } else {
                     velxminus[ijk] = velx[extrap_ijk];
                     velxplus[ijk]  = velx[extrap_ijk];
                     velyminus[ijk] = vely[extrap_ijk];
                     velyplus[ijk]  = vely[extrap_ijk];
                     velzminus[ijk] = velz[extrap_ijk];
                     velzplus[ijk]  = velz[extrap_ijk];
                  }
                  
                  if (do_Ye) {
                     Y_e_plus[ijk]  = Y_e[extrap_ijk];
                     Y_e_minus[ijk] = Y_e[extrap_ijk];
                  }
                  if (do_temp) {
                    tempminus[ijk] = temperature[extrap_ijk];
                    tempplus[ijk] = temperature[extrap_ijk];
                  }
                  if (do_mhd) {
                    Bvecxminus[ijk] = Bvecx[extrap_ijk];
                    Bvecxplus[ijk]  = Bvecx[extrap_ijk];
                    Bvecyminus[ijk] = Bvecy[extrap_ijk];
                    Bvecyplus[ijk]  = Bvecy[extrap_ijk];
                    Bveczminus[ijk] = Bvecz[extrap_ijk];
                    Bveczplus[ijk]  = Bvecz[extrap_ijk];
                    if (do_clean_divergence) {
                       psidcminus[ijk] = psidc[extrap_ijk];
                       psidcplus[ijk] = psidc[extrap_ijk];
                    }
                  }
                  if (evolve_tracer) {
                     for (int l=0; l < number_of_tracers; ++l) {
                        tracerminus[ijk + l*N] = tracer[extrap_ijk + l*N];
                        tracerplus[ijk + l*N]  = tracer[extrap_ijk + l*N];
                     }
                  }

                  // Riemann problem might not be trivial anymore!
                  change_riemann_trivial_flag = true;
               }

               if (change_riemann_trivial_flag && !excise) {
                  SpaceMask_SetStateBits(space_mask, ijk_m, type_bits, not_trivial);
                  SpaceMask_SetStateBits(space_mask, ijk,   type_bits, not_trivial);
               } else if (excise) {
                  SpaceMask_SetStateBits(space_mask, ijk_m, type_bits, trivial);
                  SpaceMask_SetStateBits(space_mask, ijk,   type_bits, trivial);
               }
            }

      
   }

  /*
     Here goes the prim2con call for reconstructed variables!
  */

  if (CCTK_EQUALS(recon_vars, "primitive") || CCTK_EQUALS(recon_method, "ppm")) {
    GRHydro_Primitive2Conservative_CC(CCTK_PASS_CTOC);
  } else if (CCTK_EQUALS(recon_vars,"conservative")) {
    CCTK_ERROR("Conservative reconstruction currently not supported!");
  } else
    CCTK_ERROR("Variable type to reconstruct not recognized.");


}

