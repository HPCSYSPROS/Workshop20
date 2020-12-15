#ifndef _GRHYDRO_TVDRECONSTRUCT_HH
#define _GRHYDRO_TVDRECONSTRUCT_HH

#include "cctk.h"

/**
   TVD reconstruction operator.
   limiter_type: 0 = minmod
                 1 = mc2
                 2 = superbee
*/
namespace tvd {
  enum limiter_type {minmod = 0, mc2 = 1, superbee = 2, numlimiters};
}
template <tvd::limiter_type limiter>
struct GRHydro_TVDReconstruct1d {

  static inline CCTK_REAL minmod_func(const CCTK_REAL dupw, 
				    const CCTK_REAL dloc);
  static inline CCTK_REAL min3(const CCTK_REAL a, 
			    const CCTK_REAL b,
			    const CCTK_REAL c);
  template <int dir>
  static inline void apply(const int nx,
                           const CCTK_REAL* const restrict a,
                           CCTK_REAL* const restrict aminus,
                           CCTK_REAL* const restrict aplus,
                           const cGH* const cctkGH,
                           const int j, const int k
                          );
}; // end struct

#endif // _GRHYDRO_TVDRECONSTRUCT_HH
