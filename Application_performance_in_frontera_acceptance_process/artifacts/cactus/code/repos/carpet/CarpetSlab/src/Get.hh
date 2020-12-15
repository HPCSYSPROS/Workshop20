#ifndef CARPETSLAB_GET_HH
#define CARPETSLAB_GET_HH

#include "cctk.h"

namespace CarpetSlab {

extern "C" CCTK_INT CarpetSlab_Get(CCTK_POINTER_TO_CONST const cctkGH,
                                   CCTK_INT const mapping_handle,
                                   CCTK_INT const proc, CCTK_INT const vindex,
                                   CCTK_INT const timelevel,
                                   CCTK_INT const hdatatype,
                                   CCTK_POINTER const hdata);

extern "C" CCTK_INT
CarpetSlab_GetList(CCTK_POINTER_TO_CONST const cctkGH,
                   CCTK_INT const mapping_handle, CCTK_INT const num_arrays,
                   CCTK_INT const *const procs, CCTK_INT const *const vindices,
                   CCTK_INT const *const timelevels,
                   CCTK_INT const *const hdatatypes,
                   CCTK_POINTER const *const hdata, CCTK_INT *const retvals);

} // namespace CarpetSlab

#endif // !defined(CARPETSLAB_GET_HH)
