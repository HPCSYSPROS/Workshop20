#ifndef _GRHYDRO_RECONSTRUCT_DRV_CXX_H
#define _GRHYDRO_RECONSTRUCT_DRV_CXX_H

#include "cctk.h"

template <class RECONSTRUCT>                   // the reconstruction operator
void GRHydro_Reconstruct_drv_cxx(cGH const * const restrict cctkGH,
                                 bool do_mhd,
                                 bool do_Ye,
                                 bool do_temp,
                                 bool do_reconstruct_Wv,
                                 bool do_reconstruct_pressure,
                                 bool do_clean_divergence);

// helper macro to instantiate all required permuations of the template options
// this must match GRHydro_Reconstruct's reconstruct::select routine
#define INSTANTIATE_RECONSTRUCTION_OPERATOR(RECONSTRUCT)                       \
template void                                                                  \
GRHydro_Reconstruct_drv_cxx <RECONSTRUCT> (                                    \
  cGH const * const restrict cctkGH,                                           \
  bool do_mhd,                                                                 \
  bool do_Ye,                                                                  \
  bool do_temp,                                                                \
  bool do_reconstruct_Wv,                                                      \
  bool do_reconstruct_pressure,                                                \
  bool do_clean_divergence);

#endif // _GRHYDRO_RECONSTRUCT_DRV_CXX_H
