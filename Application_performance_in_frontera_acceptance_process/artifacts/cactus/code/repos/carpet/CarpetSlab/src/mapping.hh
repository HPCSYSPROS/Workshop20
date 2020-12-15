#ifndef CARPETSLAB_MAPPING_HH
#define CARPETSLAB_MAPPING_HH

#include <vector>

#include "cctk.h"

namespace CarpetSlab {

// Mapping object
// (just store the mapping)
struct mapping {
  int vindex;
  int hdim;
  vector<int> origin; // [vdim]
  vector<int> dirs;   // [hdim]
  vector<int> stride; // [hdim]
  vector<int> length; // [hdim]
};

int StoreMapping(mapping *const mp);

mapping *RetrieveMapping(int const table);

void DeleteMapping(int const table);

typedef CCTK_INT (*conversion_fn_ptr)(
    CCTK_INT const nelems, CCTK_INT const src_stride, CCTK_INT const dst_stride,
    CCTK_INT const src_type, CCTK_INT const dst_type,
    CCTK_POINTER_TO_CONST const from, CCTK_POINTER const to);

extern "C" CCTK_INT CarpetSlab_LocalMappingByIndex(
    CCTK_POINTER_TO_CONST const cctkGH, CCTK_INT const vindex,
    CCTK_INT const hdim, CCTK_INT const *const direction,
    CCTK_INT const *const origin, CCTK_INT const *const extent,
    CCTK_INT const *const downsample, CCTK_INT const table_handle,
    conversion_fn_ptr const conversion_fn, CCTK_INT *const hsize_local,
    CCTK_INT *const hsize_global, CCTK_INT *const hoffset_global);

extern "C" CCTK_INT CarpetSlab_GlobalMappingByIndex(
    CCTK_POINTER_TO_CONST const cctkGH, CCTK_INT const vindex,
    CCTK_INT const hdim, CCTK_INT const *const direction,
    CCTK_INT const *const origin, CCTK_INT const *const extent,
    CCTK_INT const *const downsample, CCTK_INT const table_handle,
    conversion_fn_ptr const conversion_fn, CCTK_INT *const hsize);

extern "C" CCTK_INT CarpetSlab_FreeMapping(CCTK_INT const mapping_handle);

} // namespace CarpetSlab

#endif // !defined(CARPETSLAB_MAPPING_HH)
