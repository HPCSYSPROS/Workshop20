#ifndef CARPETSLAB_GETHYPERSLAB_HH
#define CARPETSLAB_GETHYPERSLAB_HH

#include "GetHyperslab.h"

namespace CarpetSlab {

void *GetSlab(const cGH *const cgh, const int dest_proc, const int n,
              const int tl, const int hdim, const int origin[/*vdim*/],
              const int dirs[/*hdim*/], const int stride[/*hdim*/],
              const int length[/*hdim*/]);

} // namespace CarpetSlab

#endif // !defined(CARPETSLAB_GETHYPERSLAB_HH)
