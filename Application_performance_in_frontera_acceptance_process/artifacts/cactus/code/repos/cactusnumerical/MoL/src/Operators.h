#ifndef OPERATORS_H
#define OPERATORS_H

#include <cctk.h>

CCTK_INT
MoL_LinearCombination(cGH const *const cctkGH,
                      CCTK_INT   const var,
                      CCTK_INT   const rl,
                      CCTK_INT   const tl,
                      CCTK_REAL  const scale,
                      CCTK_INT   const srcs[],
                      CCTK_INT   const tls[],
                      CCTK_REAL  const facts[],
                      CCTK_INT   const nsrcs);

#endif  // #ifndef OPERATORS_H
