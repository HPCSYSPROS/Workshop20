#ifndef _GRHYDRO_FUNCTIONS_H_
#define _GRHYDRO_FUNCTIONS_H_

#include "cctk.h"

#ifdef __cplusplus
extern "C" {
#endif

void GRHydro_Primitive2Conservative_CC(CCTK_ARGUMENTS);

#ifdef __cplusplus
}
#endif

#endif // _GRHYDRO_FUNCTIONS_H_
