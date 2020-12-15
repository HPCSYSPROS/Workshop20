#include "vectors.h"

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

extern "C"
int Vectors_Startup(void)
{
  CCTK_VInfo(CCTK_THORNSTRING,
             "Using vector size %d for architecture %s",
             CCTK_REAL_VEC_SIZE, vec_architecture);
  return 0;
}
