#ifndef CARPETSLAB_SLAB_HH
#define CARPETSLAB_SLAB_HH

#include "cctk.h"

namespace CarpetSlab {

void FillSlab(const cGH *const cgh, const int dest_proc, const int n,
              const int tl, const int hdim, const int origin[/*vdim*/],
              const int dirs[/*hdim*/], const int stride[/*hdim*/],
              const int length[/*hdim*/], void *const hdata);

} // namespace CarpetSlab

#endif // !defined(CARPETSLAB_SLAB_HH)
