/*@@
  @file      Courant.c
  @date      September 4 1999
  @author    Gabrielle Allen
  @desc
             Specification of timestep using Courant condition
  @enddesc
@@*/

#include <math.h>
#include <stdlib.h>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

void Time_Courant(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  /* Calculate the minimum grid spacing */
  CCTK_REAL min_spacing = cctk_delta_space[0];

  if (cctk_dim >= 2) {
    min_spacing =
        (min_spacing < cctk_delta_space[1] ? min_spacing : cctk_delta_space[1]);
  }

  if (cctk_dim >= 3) {
    min_spacing =
        (min_spacing < cctk_delta_space[2] ? min_spacing : cctk_delta_space[2]);
  }

  if (cctk_dim >= 4) {
    CCTK_WARN(0, "Time Step not defined for greater than 4 dimensions");
  }

  /* Calculate the courant timestep */
  if (CCTK_Equals(timestep_method, "courant_time")) {
    *courant_dt = courant_fac *(*courant_min_time) / sqrt((double)cctk_dim);
  } else if (CCTK_Equals(timestep_method, "courant_speed")) {
    *courant_dt = courant_fac *min_spacing / (*courant_wave_speed) /
                  sqrt((double)cctk_dim);
  }

  if (CCTK_Equals(terminate, "time") || CCTK_Equals(terminate, "both")) {
    if (cctkGH->cctk_time + *courant_dt > cctk_final_time) {
      *courant_dt = (1 + 1.e-10) * (cctk_final_time - cctkGH->cctk_time);
    }
  }

  /* Set the Cactus timestep */

  if (!timestep_outonly) {
    cctkGH->cctk_delta_time = *courant_dt;
    if (verbose) {
      CCTK_VInfo(CCTK_THORNSTRING, "Time step set to %g", CCTK_DELTA_TIME);
    }
  } else {
    cctkGH->cctk_delta_time = dtfac * min_spacing;
    if (cctkGH->cctk_iteration % timestep_outevery == 0) {
      CCTK_VInfo(CCTK_THORNSTRING, "Courant timestep would be %g", *courant_dt);
    }
  }
}
