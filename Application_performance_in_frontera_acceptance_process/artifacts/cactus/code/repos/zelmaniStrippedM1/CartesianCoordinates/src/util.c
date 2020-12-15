#include "cctk.h"

CCTK_INT
CartesianCoordinates_GetBbox(CCTK_POINTER_TO_CONST const GH,
                             CCTK_INT const size,
                             CCTK_INT * const bbox)
{
  const cGH* const cctkGH = (const cGH* const) GH;  

  for (int i=0; i<6; ++i)
    bbox[i] = cctkGH->cctk_bbox[i];

  return 0;
}


CCTK_INT
CartesianCoordinates_GetMap(CCTK_POINTER_TO_CONST const GH)
{
  const cGH* const cctkGH = (const cGH* const) GH;  
  return 0;
}
