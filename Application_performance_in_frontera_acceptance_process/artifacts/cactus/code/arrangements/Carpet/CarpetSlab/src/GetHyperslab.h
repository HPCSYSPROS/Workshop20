#ifndef CARPETSLAB_GETHYPERSLAB_H
#define CARPETSLAB_GETHYPERSLAB_H

#include "cctk.h"

#ifdef __cplusplus
namespace CarpetSlab {
extern "C" {
#endif

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

#endif /* !defined(CARPETSLAB_GETHYPERSLAB_H) */
