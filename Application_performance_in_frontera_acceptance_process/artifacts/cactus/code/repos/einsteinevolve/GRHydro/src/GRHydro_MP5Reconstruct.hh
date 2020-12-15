#ifndef _GRHYDRO_MP5RECONSTRUCT_HH
#define _GRHYDRO_MP5RECONSTRUCT_HH

#include "cctk.h"

using namespace std;


/**
   MP5 reconstruction operator.
*/
template <bool do_MP5_adaptive_epsilon>
struct GRHydro_MP5Reconstruct1d_cxx {

template <int dir>
static inline void apply(const int nx,
                         const CCTK_REAL* const restrict a,
                         CCTK_REAL* const restrict aminus,
                         CCTK_REAL* const restrict aplus,
                         const cGH* const cctkGH,
                         const int j, const int k
                        );
}; // end struct

#endif // _GRHYDRO_MP5RECONSTRUCT_HH
