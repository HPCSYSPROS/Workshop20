#ifndef VECTORS_HH
#define VECTORS_HH

#include <cassert>
#include <cmath>
#include <cstdlib>

  // Default vector implementation, does not vectorise
#include "vectors-default.hh"
  
#if 0
  // Intel SSE vector instructions
#include "vectors-intel.hh"
  
  // Power (Altivec) vector instructions
#include "vectors-power.hh"
#endif

#endif  // #ifndef VECTORS_HH
