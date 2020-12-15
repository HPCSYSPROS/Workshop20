#ifndef _GRHYDRO_RECONSTRUCT_DRV_CXX_H
#define _GRHYDRO_RECONSTRUCT_DRV_CXX_H

#include "cctk.h"

/*
  Cases that must be considered:
  * basic hydro
  * hydro + temperature + ye
  * hydro + ye
  * basic mhd
  * mhd + temperature + ye 
  * mhd + ye 
 */

template <bool do_mhd,
          bool do_Ye,
          bool do_temp,
          bool do_reconstruct_Wv,
          bool do_clean_divergence,
          class RECONSTRUCT>                   // the reconstruction operator
void GRHydro_Reconstruct_drv_cxx(cGH const * const restrict cctkGH);

// helper macro to instantiate all required permuations of the template options
// this must match GRHydro_Reconstruct's reconstruct::select routine
#define INSTANTIATE_RECONSTRUCTION_OPERATOR(RECONSTRUCT)                       \
template void                                                                  \
GRHydro_Reconstruct_drv_cxx <false, false, false, true, false, RECONSTRUCT> (  \
  cGH const * const restrict cctkGH);                                          \
template void                                                                  \
GRHydro_Reconstruct_drv_cxx <true, false, false, true, false, RECONSTRUCT> (   \
  cGH const * const restrict cctkGH);                                          \
template void                                                                  \
GRHydro_Reconstruct_drv_cxx <true, true, true, true, false, RECONSTRUCT> (     \
  cGH const * const restrict cctkGH);                                          \
template void                                                                  \
GRHydro_Reconstruct_drv_cxx <true, false, false, true, true, RECONSTRUCT> (    \
  cGH const * const restrict cctkGH);                                          \
template void                                                                  \
GRHydro_Reconstruct_drv_cxx <true, true, true, true, true, RECONSTRUCT> (      \
  cGH const * const restrict cctkGH);                                          \
template void                                                                  \
GRHydro_Reconstruct_drv_cxx <false, true, true, true, false, RECONSTRUCT> (    \
  cGH const * const restrict cctkGH);                                          \
template void                                                                  \
GRHydro_Reconstruct_drv_cxx <false, false, false, false, false, RECONSTRUCT> ( \
  cGH const * const restrict cctkGH);                                          \
template void                                                                  \
GRHydro_Reconstruct_drv_cxx <true, false, false, false, false, RECONSTRUCT> (  \
  cGH const * const restrict cctkGH);                                          \
template void                                                                  \
GRHydro_Reconstruct_drv_cxx <true, true, true, false, false, RECONSTRUCT> (    \
  cGH const * const restrict cctkGH);                                          \
template void                                                                  \
GRHydro_Reconstruct_drv_cxx <true, false, false, false, true, RECONSTRUCT> (   \
  cGH const * const restrict cctkGH);                                          \
template void                                                                  \
GRHydro_Reconstruct_drv_cxx <true, true, true, false, true, RECONSTRUCT> (     \
  cGH const * const restrict cctkGH);                                          \
template void                                                                  \
GRHydro_Reconstruct_drv_cxx <false, true, true, false, false, RECONSTRUCT> (   \
  cGH const * const restrict cctkGH);

#endif // _GRHYDRO_RECONSTRUCT_DRV_CXX_H
