#ifndef STENCIL_H
#define STENCIL_H

#include "cctk.h"

void SBP_determine_onesided_stencil (const cGH * cctkGH, int * onesided);

#endif /* #ifndef STENCIL_H */
