#include <cassert>

#include "cctk.h"

#include "carpet.hh"

#include "mapping.hh"
#include "slab.hh"
#include "Get.hh"

namespace CarpetSlab {

using namespace Carpet;

CCTK_INT
CarpetSlab_Get(CCTK_POINTER_TO_CONST const cctkGH_,
               CCTK_INT const mapping_handle, CCTK_INT const proc,
               CCTK_INT const vindex, CCTK_INT const timelevel,
               CCTK_INT const hdatatype, CCTK_POINTER const hdata) {
  cGH const *const cctkGH = (cGH const *)cctkGH_;

  // Check arguments
  assert(cctkGH);
  assert(mapping_handle >= 0);
  assert(proc == -1 || (proc >= 0 && proc < CCTK_nProcs(cctkGH)));
  assert(vindex >= 0 && vindex < CCTK_NumVars());
  assert(timelevel >= 0);
  assert(hdatatype >= 0);
  assert(hdata);

  // Get mapping
  const mapping *const mp = RetrieveMapping(mapping_handle);
  assert(mp);

  // Calculate total size
  size_t size = 1;
  for (size_t d = 0; d < (size_t)mp->hdim; ++d) {
    size *= mp->length[d];
  }

  // Get type size
  size_t const sz = CCTK_VarTypeSize(hdatatype);
  assert(sz > 0);

  // Forward call
  FillSlab(cctkGH, proc, vindex, timelevel, mp->hdim, &mp->origin[0],
           &mp->dirs[0], &mp->stride[0], &mp->length[0], hdata);

  return 0;
}

CCTK_INT
CarpetSlab_GetList(CCTK_POINTER_TO_CONST const cctkGH_,
                   CCTK_INT const mapping_handle, CCTK_INT const num_arrays,
                   CCTK_INT const *const procs, CCTK_INT const *const vindices,
                   CCTK_INT const *const timelevels,
                   CCTK_INT const *const hdatatypes,
                   CCTK_POINTER const *const hdata, CCTK_INT *const retvals) {
  cGH const *const cctkGH = (cGH const *)cctkGH_;

  // Check arguments
  assert(cctkGH);
  assert(mapping_handle >= 0);
  assert(num_arrays >= 0);
  assert(procs);
  assert(vindices);
  assert(timelevels);
  assert(hdatatypes);
  assert(hdata);
  assert(retvals);

  // Remember whether there were errors
  bool everyting_okay = true;

  // Loop over all slabs
  for (int n = 0; n < num_arrays; ++n) {
    // Forward call
    retvals[n] = CarpetSlab::CarpetSlab_Get(cctkGH, mapping_handle, procs[n],
                                            vindices[n], timelevels[n],
                                            hdatatypes[n], hdata[n]);
    everyting_okay = everyting_okay && retvals[n] == 0;
  }

  return everyting_okay ? num_arrays : -1;
}

} // namespace CarpetSlab
