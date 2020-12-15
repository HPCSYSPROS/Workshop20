#ifndef CARPETSLAB_H
#define CARPETSLAB_H

#include "cctk.h"

#ifdef __cplusplus
namespace CarpetSlab {
extern "C" {
#endif

CCTK_INT CarpetSlab_Get(CCTK_POINTER_TO_CONST const cctkGH,
                        CCTK_INT const mapping_handle, CCTK_INT const proc,
                        CCTK_INT const vindex, CCTK_INT const timelevel,
                        CCTK_INT const hdatatype, CCTK_POINTER const hdata);

CCTK_INT
CarpetSlab_GetList(CCTK_POINTER_TO_CONST const cctkGH,
                   CCTK_INT const mapping_handle, CCTK_INT const num_arrays,
                   CCTK_INT const *const procs, CCTK_INT const *const vindices,
                   CCTK_INT const *const timelevels,
                   CCTK_INT const *const hdatatypes,
                   CCTK_POINTER const *const hdata, CCTK_INT *const retvals);

CCTK_INT CarpetSlab_LocalMappingByIndex(
    CCTK_POINTER_TO_CONST const cctkGH, CCTK_INT const vindex,
    CCTK_INT const hdim, CCTK_INT const *const direction,
    CCTK_INT const *const origin, CCTK_INT const *const extent,
    CCTK_INT const *const downsample, CCTK_INT const table_handle,
    CCTK_INT (*const conversion_fn)(
        CCTK_INT const nelems, CCTK_INT const src_stride,
        CCTK_INT const dst_stride, CCTK_INT const src_type,
        CCTK_INT const dst_type, void const *const from, CCTK_POINTER const to),
    CCTK_INT *const hsize_local, CCTK_INT *const hsize_global,
    CCTK_INT *const hoffset_global);

CCTK_INT CarpetSlab_GlobalMappingByIndex(
    CCTK_POINTER_TO_CONST const cctkGH, CCTK_INT const vindex,
    CCTK_INT const hdim, CCTK_INT const *const direction,
    CCTK_INT const *const origin, CCTK_INT const *const extent,
    CCTK_INT const *const downsample, CCTK_INT const table_handle,
    CCTK_INT (*const conversion_fn)(
        CCTK_INT const nelems, CCTK_INT const src_stride,
        CCTK_INT const dst_stride, CCTK_INT const src_type,
        CCTK_INT const dst_type, void const *const from, CCTK_POINTER const to),
    CCTK_INT *const hsize);

CCTK_INT CarpetSlab_FreeMapping(CCTK_INT const mapping_handle);

/* Old interface -- don't use */
int Hyperslab_GetHyperslab(const cGH *const GH, const int target_proc,
                           const int vindex, const int vtimelvl, const int hdim,
                           const int global_startpoint[/*vdim*/],
                           const int directions[/*vdim*/],
                           const int lengths[/*hdim*/],
                           const int downsample[/*hdim*/], void **const hdata,
                           int hsize[/*hdim*/]);

#ifdef __cplusplus
} /* extern "C" */
} /* namespace CarpetSlab */
#endif

#endif /* !defined(CARPETSLAB_H) */
