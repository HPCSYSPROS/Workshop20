/* $Header$ */

#ifndef TENSORTYPES_IMP_H
#define TENSORTYPES_IMP_H

#include "cctk_Arguments.h"

#include "tensortypes.h"



extern struct tensor const TT_scalar;
extern struct tensor const TT_vector;
extern struct tensor const TT_tensor;
extern struct tensor const TT_symmtensor;
extern struct tensor const TT_symmtensor3a;
extern struct tensor const TT_symmtensor3b;
extern struct tensor const TT_symmtensor3c;
extern struct tensor const TT_symmtensor4;

extern struct tensor const TT_4scalar;
extern struct tensor const TT_4vector;
extern struct tensor const TT_4symmtensor;

extern struct tensor const TT_weylscalars;
extern struct tensor const TT_weylscalars_real;

extern int const TT_numtensors;
extern struct tensor const * const TT_alltensors[];



void
CheckTensorType (struct tensor const * restrict const atensor);

void
CheckTensorTypes (CCTK_ARGUMENTS);



#endif /* TENSORTYPES_IMP_H */
