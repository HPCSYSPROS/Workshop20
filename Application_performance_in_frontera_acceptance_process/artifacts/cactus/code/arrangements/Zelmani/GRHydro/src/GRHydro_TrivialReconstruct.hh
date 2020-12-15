#ifndef _GRHYDRO_TRIVIALRECONSTRUCT_HH
#define _GRHYDRO_TRIVIALRECONSTRUCT_HH

#include "cctk.h"

/**
   Trivial first-order reconstruction operator.
*/
struct GRHydro_TrivialReconstruct1d {
template <int dir>
static inline void apply(const int nx,
                         const CCTK_REAL* const restrict a,
                         CCTK_REAL* const restrict aminus,
                         CCTK_REAL* const restrict aplus,
                         const cGH* const cctkGH,
                         const int j, const int k
                        );
}; // end struct

#endif // _GRHYDRO_TRIVIALRECONSTRUCT_HH
