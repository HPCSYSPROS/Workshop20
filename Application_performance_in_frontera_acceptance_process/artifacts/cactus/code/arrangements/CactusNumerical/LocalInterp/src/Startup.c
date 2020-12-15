/*@@
  @file      Startup.c
  @date      Sun Jul 04 1999
  @author    Thomas Radke
  @desc
             Startup routines for LocalInterp/UniformCartesian
  @enddesc

  @history
  @date      22 Jan 2002
  @author    Jonathan Thornburg
  @hdesc     Move all local-interpolation code from LocalInterp to here
  @endhistory

  @version   $Header$
  @@*/

#include "cctk.h"
#include "Interpolate.h"

/* the rcs ID and its dummy function to use it */
static const char *rcsid = "$Header$";
CCTK_FILEVERSION(CactusBase_LocalInterp_Startup_c)


/********************************************************************
 ***************    External Routine Prototypes   *******************
 ********************************************************************/
int LocalInterp_Startup(void);


/*@@
  @routine   LocalInterp_Startup
  @date      Sun Jul 04 1999
  @author    Thomas Radke
  @desc
             The startup registration routine for LocalInterp.
             Registers the "uniform cartesian" interpolation operator
             with the flesh.
  @enddesc
  @calls     CCTK_InterpRegisterOpLocalUniform
  @@*/
int LocalInterp_Startup (void)
{
  CCTK_InterpRegisterOpLocalUniform (LocalInterp_InterpLocalUniform,
                                     "uniform cartesian", CCTK_THORNSTRING);
  return 0;
}
