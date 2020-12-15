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
void CCTK_FCALL CCTK_FNAME(GRHydro_PPMReconstruct_drv_opt)(cGH /*const*/ * const restrict & /*restrict*/ cctkGH);

extern "C"
void CCTK_FCALL CCTK_FNAME(primitive2conservativeM)(cGH /*const*/ * const restrict & /*restrict*/ cctkGH);


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
                             const bool clean_divergence,
                             cGH const * const restrict cctkGH)
   {
       if (!do_mhd && !do_Ye && !do_temp && reconstruct_Wv)
          GRHydro_Reconstruct_drv_cxx <false, false, false, true, false, RECONSTRUCT > (cctkGH);
       else if (do_mhd && !do_Ye && !do_temp && reconstruct_Wv && !clean_divergence)
          GRHydro_Reconstruct_drv_cxx <true, false, false, true, false, RECONSTRUCT > (cctkGH);
       else if (do_mhd && do_Ye && do_temp && reconstruct_Wv && !clean_divergence)
          GRHydro_Reconstruct_drv_cxx <true, true, true, true, false, RECONSTRUCT  > (cctkGH);
       else if (do_mhd && !do_Ye && !do_temp && reconstruct_Wv && clean_divergence)
          GRHydro_Reconstruct_drv_cxx <true, false, false, true, true, RECONSTRUCT  > (cctkGH);
       else if (do_mhd && do_Ye && do_temp && reconstruct_Wv && clean_divergence)
          GRHydro_Reconstruct_drv_cxx <true, true, true, true, true, RECONSTRUCT  > (cctkGH);
       else if (!do_mhd && do_Ye && do_temp && reconstruct_Wv)
          GRHydro_Reconstruct_drv_cxx <false, true, true, true, false, RECONSTRUCT  > (cctkGH);
       else if (!do_mhd && !do_Ye && !do_temp && !reconstruct_Wv)
          GRHydro_Reconstruct_drv_cxx <false, false, false, false, false, RECONSTRUCT  > (cctkGH);
       else if (do_mhd && !do_Ye && !do_temp && !reconstruct_Wv && !clean_divergence)
          GRHydro_Reconstruct_drv_cxx <true, false, false, false, false, RECONSTRUCT  > (cctkGH);
       else if (do_mhd && do_Ye && do_temp && !reconstruct_Wv && !clean_divergence)
          GRHydro_Reconstruct_drv_cxx <true, true, true, false, false, RECONSTRUCT  > (cctkGH);
       else if (do_mhd && !do_Ye && !do_temp && !reconstruct_Wv && clean_divergence)
          GRHydro_Reconstruct_drv_cxx <true, false, false, false, true, RECONSTRUCT  > (cctkGH);
       else if (do_mhd && do_Ye && do_temp && !reconstruct_Wv && clean_divergence)
          GRHydro_Reconstruct_drv_cxx <true, true, true, false, true, RECONSTRUCT  > (cctkGH);
       else if (!do_mhd && do_Ye && do_temp && !reconstruct_Wv)
          GRHydro_Reconstruct_drv_cxx <false, true, true, false, false, RECONSTRUCT  > (cctkGH);
       else
          assert(0&&"Unsupported template arguments");
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
     First initialize all variables by first order trivial reconstruction!
  */


  if (!do_mhd && !do_Ye && !do_temp)
      GRHydro_Reconstruct_drv_cxx <false, false, false, false, false, GRHydro_TrivialReconstruct1d > (cctkGH);
  else if (do_mhd && !do_Ye && !do_temp && !clean_divergence)
      GRHydro_Reconstruct_drv_cxx <true, false, false, false, false, GRHydro_TrivialReconstruct1d > (cctkGH);
  else if (do_mhd && !do_Ye && !do_temp && clean_divergence)
      GRHydro_Reconstruct_drv_cxx <true, false, false, false, true, GRHydro_TrivialReconstruct1d > (cctkGH);
  else if (do_mhd && do_Ye && do_temp && !clean_divergence)
      GRHydro_Reconstruct_drv_cxx <true, true, true, false, false, GRHydro_TrivialReconstruct1d > (cctkGH);
  else if (do_mhd && do_Ye && do_temp && clean_divergence)
      GRHydro_Reconstruct_drv_cxx <true, true, true, false, true, GRHydro_TrivialReconstruct1d > (cctkGH);
  else if (!do_mhd && do_Ye && do_temp)
      GRHydro_Reconstruct_drv_cxx <false, true, true, false, false, GRHydro_TrivialReconstruct1d > (cctkGH);
  else
      CCTK_ERROR("Don't know how to initialize reconstructed variables");


  /*
     Now do the real reconstruction!
  */

  if (CCTK_EQUALS(recon_method,"tvd")) {
    // this handles MHD and non-MHD

    if (CCTK_EQUALS(tvd_limiter, "minmod")) {
       
       reconstruct<GRHydro_TVDReconstruct1d<tvd::minmod> >::select(do_mhd, do_Ye, do_temp, reconstruct_Wv, clean_divergence, cctkGH);
       
    } else if (CCTK_EQUALS(tvd_limiter, "vanleerMC2")) {
       
       reconstruct<GRHydro_TVDReconstruct1d<tvd::mc2> >::select(do_mhd, do_Ye, do_temp, reconstruct_Wv, clean_divergence, cctkGH);
       
    } else if (CCTK_EQUALS(tvd_limiter, "superbee")) {
       
       reconstruct<GRHydro_TVDReconstruct1d<tvd::superbee> >::select(do_mhd, do_Ye, do_temp, reconstruct_Wv, clean_divergence, cctkGH);
    
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
       
       reconstruct<GRHydro_WENOReconstruct1d_cxx<false, true> >::select(do_mhd, do_Ye, do_temp, reconstruct_Wv, clean_divergence, cctkGH);
          
    } else {
       
       reconstruct<GRHydro_WENOReconstruct1d_cxx<false, false> >::select(do_mhd, do_Ye, do_temp, reconstruct_Wv, clean_divergence, cctkGH);
    
    }

  } else if (CCTK_EQUALS(recon_method, "weno-z")) {

       reconstruct<GRHydro_WENOReconstruct1d_cxx<true, false> >::select(do_mhd, do_Ye, do_temp, reconstruct_Wv, clean_divergence, cctkGH);

  } else if (CCTK_EQUALS(recon_method,"mp5")) {
    // this handles MHD and non-MHD

    if (mp5_adaptive_eps)
       reconstruct<GRHydro_MP5Reconstruct1d_cxx<true> >::select(do_mhd, do_Ye, do_temp, reconstruct_Wv, clean_divergence, cctkGH);
    else
       reconstruct<GRHydro_MP5Reconstruct1d_cxx<false> >::select(do_mhd, do_Ye, do_temp, reconstruct_Wv, clean_divergence, cctkGH);

  } else
    CCTK_ERROR("Reconstruction method not recognized!");


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
               for (int l=-(GRHydro_stencil-1); l < GRHydro_stencil-1; ++l) {
                  const int ijk_l = CCTK_GFINDEX3D(cctkGH, i+*xoffset*l, j+*yoffset*l, k+*zoffset*l);
                  if (hydro_excision_mask[ijk_l]) {
                     first_order = true;
                     break;
                  }
               }
               
               // if the current point is excised, then we will mark the Riemann problem as trivial.
               if (hydro_excision_mask[ijk])
                  excise = true;
               
               if (!excise && hydro_excision_mask[ijk_m]) {
                  extrapolate = true;
                  extrap_ijk = ijk_m;
               }
               
               if (!excise && hydro_excision_mask[ijk_p]) {
                  extrapolate = true;
                  extrap_ijk = ijk_p;
               }
               
               if (first_order)
               {
                  rhominus[ijk] = rho[ijk];
                  rhoplus[ijk]  = rho[ijk];
                  
                  epsminus[ijk] = eps[ijk];
                  epsplus[ijk]  = eps[ijk];
                  
                  velxminus[ijk] = velx[ijk];
                  velxplus[ijk]  = velx[ijk];
                  velyminus[ijk] = vely[ijk];
                  velyplus[ijk]  = vely[ijk];
                  velzminus[ijk] = velz[ijk];
                  velzplus[ijk]  = velz[ijk];
                  
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
               
               if (extrapolate)
               {
                  rhominus[extrap_ijk] = rho[ijk];
                  rhoplus[extrap_ijk]  = rho[ijk];
                  
                  epsminus[extrap_ijk] = eps[ijk];
                  epsplus[extrap_ijk]  = eps[ijk];
                  
                  velxminus[extrap_ijk] = velx[ijk];
                  velxplus[extrap_ijk]  = velx[ijk];
                  velyminus[extrap_ijk] = vely[ijk];
                  velyplus[extrap_ijk]  = vely[ijk];
                  velzminus[extrap_ijk] = velz[ijk];
                  velzplus[extrap_ijk]  = velz[ijk];
                  
                  if (do_Ye) {
                     Y_e_plus[extrap_ijk]  = Y_e[ijk];
                     Y_e_minus[extrap_ijk] = Y_e[ijk];
                  }
                  if (do_temp) {
                    tempminus[extrap_ijk] = temperature[ijk];
                    tempplus[extrap_ijk] = temperature[ijk];
                  }
                  if (do_mhd) {
                    Bvecxminus[extrap_ijk] = Bvecx[ijk];
                    Bvecxplus[extrap_ijk]  = Bvecx[ijk];
                    Bvecyminus[extrap_ijk] = Bvecy[ijk];
                    Bvecyplus[extrap_ijk]  = Bvecy[ijk];
                    Bveczminus[extrap_ijk] = Bvecz[ijk];
                    Bveczplus[extrap_ijk]  = Bvecz[ijk];
                    if (do_clean_divergence) {
                       psidcminus[extrap_ijk] = psidc[ijk];
                       psidcplus[extrap_ijk] = psidc[ijk];
                    }
                  }
                  if (evolve_tracer) {
                     for (int l=0; l < number_of_tracers; ++l) {
                        tracerminus[extrap_ijk + l*N] = tracer[ijk + l*N];
                        tracerplus[extrap_ijk + l*N]  = tracer[ijk + l*N];
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
     Here goes the atmosphere treatment!
  */
  #pragma omp parallel for
  for (int k=GRHydro_stencil-1; k < cctk_lsh[2] - GRHydro_stencil+1+transport_constraints*(1-*zoffset); ++k)
     for (int j=GRHydro_stencil-1; j < cctk_lsh[1] - GRHydro_stencil+1+transport_constraints*(1-*yoffset); ++j)
        for (int i=GRHydro_stencil-1; i < cctk_lsh[0] - GRHydro_stencil+1+transport_constraints*(1-*xoffset); ++i) {
           const int ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
           const int ijk_m = CCTK_GFINDEX3D(cctkGH, i-*xoffset, j-*yoffset, k-*zoffset);
           
           // Whenever we modify a previously reconstructed value,
           // we need to change the trivial flag!
           bool change_riemann_trivial_flag = false;
           
           if (rhominus[ijk] < *GRHydro_rho_min || rhoplus[ijk] < *GRHydro_rho_min)
           {
              rhominus[ijk] = rho[ijk];
              rhoplus[ijk]  = rho[ijk];
              
              epsminus[ijk] = eps[ijk];
              epsplus[ijk]  = eps[ijk];
              
              velxminus[ijk] = velx[ijk];
              velxplus[ijk]  = velx[ijk];
              velyminus[ijk] = vely[ijk];
              velyplus[ijk]  = vely[ijk];
              velzminus[ijk] = velz[ijk];
              velzplus[ijk]  = velz[ijk];
              // TODO: add reconstruct_Wv treatment
              
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
           
           if (!do_temp) {
              if (epsplus[ijk] < 0.0) {
                 epsplus[ijk] = eps[ijk];
                 change_riemann_trivial_flag = true;
              }
              if (epsminus[ijk] < 0.0) {
                 epsminus[ijk] = eps[ijk];
                 change_riemann_trivial_flag = true;
              }
           }
           
           if (change_riemann_trivial_flag) {
              // Riemann problem might not be trivial anymore!
              SpaceMask_SetStateBits(space_mask, ijk_m, type_bits, not_trivial);
              SpaceMask_SetStateBits(space_mask, ijk,   type_bits, not_trivial);
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

