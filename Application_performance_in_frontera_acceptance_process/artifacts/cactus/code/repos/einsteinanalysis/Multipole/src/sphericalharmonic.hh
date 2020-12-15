
#ifndef __sphericalharmonic_h
#define __sphericalharmonic_h
#include "cctk.h"

void Multipole_HarmonicSetup(int s, int l, int m,
                             int array_size, CCTK_REAL const th[], CCTK_REAL const ph[],
                             CCTK_REAL reY[], CCTK_REAL imY[]);


#endif
