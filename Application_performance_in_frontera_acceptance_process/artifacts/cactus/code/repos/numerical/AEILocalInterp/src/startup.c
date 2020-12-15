 /*@@
   @file      startup.c
   @date      Thu Feb 21 16:25:35 CET 2002
   @author    Jonathan Thornburg <jthorn@aei.mpg.de>
   @desc
              Startup routines for thorn AEILocalInterp
   @enddesc
   @version   $Id$
 @@*/

#include <stdlib.h>
#include <limits.h>

#include "cctk.h"
#include "cctk_Interp.h"
#include "InterpLocalUniform.h"

/* the rcs ID and its dummy function to use it */
static const char *rcsid = "$Header$";
CCTK_FILEVERSION(AEIThorns_AEILocalInterp_src_startup_c)


/*@@
  @routine   LocalInterp_GPU_Startup
  @date      Thu Feb 21 16:27:41 CET 2002
  @author    Jonathan Thornburg <jthorn@aei.mpg.de>
  @desc      This is the startup routine for thorn AEILocalInterp.
             It registers the interpolation operators.
  @enddesc
  @@*/
int AEILocalInterp_U_Startup(void)
{
CCTK_InterpRegisterOpLocalUniform(AEILocalInterp_U_Lagrange_TP,
			  "Lagrange polynomial interpolation (tensor product)",
				  CCTK_THORNSTRING);
CCTK_InterpRegisterOpLocalUniform(AEILocalInterp_U_Lagrange_MD,
			  "Lagrange polynomial interpolation (maximum degree)",
				  CCTK_THORNSTRING);

CCTK_InterpRegisterOpLocalUniform(AEILocalInterp_U_Hermite,
				  "Hermite polynomial interpolation",
				  CCTK_THORNSTRING);

/* synonym operator names for backwards compatability */
CCTK_InterpRegisterOpLocalUniform(AEILocalInterp_U_Lagrange_TP,
				  "Lagrange polynomial interpolation",
				  CCTK_THORNSTRING);
CCTK_InterpRegisterOpLocalUniform(AEILocalInterp_U_Lagrange_TP,
                                  "generalized polynomial interpolation",
				  CCTK_THORNSTRING);

return 0;
}
