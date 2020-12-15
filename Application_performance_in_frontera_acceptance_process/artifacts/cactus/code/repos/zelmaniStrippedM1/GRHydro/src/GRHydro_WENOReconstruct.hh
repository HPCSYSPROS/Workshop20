#ifndef _GRHYDRO_WENORECONSTRUCT_HH
#define _GRHYDRO_WENORECONSTRUCT_HH

#include "cctk.h"

using namespace std;

/**
   WENO5 reconstruction operator.
   Supports standard WENO5 (with and without adaptive epsilon), and WENO-z.
*/
template <bool do_wenoz, bool do_adaptive_epsilon>
struct GRHydro_WENOReconstruct1d_cxx {
template <int dir>
static inline void apply(const int nx,
                         const CCTK_REAL* const restrict a,
                         CCTK_REAL* const restrict aminus,
                         CCTK_REAL* const restrict aplus,
                         const cGH* const cctkGH,
                         const int j, const int k
                        );

}; // end struct

#endif // _GRHYDRO_WENORECONSTRUCT_HH
