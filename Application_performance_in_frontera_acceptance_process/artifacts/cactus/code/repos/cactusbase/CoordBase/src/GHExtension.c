/*@@
  @file      GHExtension.c
  @date      Sun Sept 22 2002
  @author    Gabrielle Allen and David Rideout
  @desc
             CoordBase GHExtension setup
  @enddesc
  @version   $Id$
@@*/

#include <stdlib.h>

#include "cctk.h"
#include "CoordBase.h"
#include "coordbaseGH.h"
#include "util_Hash.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusBase_CoordBase_GHExtension_c);

/********************************************************************
 *********************  Scheduled Routine Prototypes  ***************
 ********************************************************************/

int CoordBase_Startup(void);

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

static void *CoordBase_SetupGH(tFleshConfig *config, int conv_level, cGH *GH);

/********************************************************************
 *********************     External Routines   **********************
 ********************************************************************/

/*@@
  @routine   CoordBase_Startup
  @date      Sunday Sept 22 2002
  @author    Gabrielle Allen
  @desc
             The startup registration routine for CoordBase.
             Registers the GH extension needed for CoordBase.
  @enddesc
  @calls     CCTK_RegisterGHExtension
             CCTK_RegisterGHExtensionSetupGH
@@*/
int CoordBase_Startup(void) {
  int GHex_handle;

  GHex_handle = CCTK_RegisterGHExtension("CoordBase");
  CCTK_RegisterGHExtensionSetupGH(GHex_handle, CoordBase_SetupGH);

  return 0;
}

/********************************************************************
 ********************    Internal Routines   ************************
 ********************************************************************/

/*@@
  @routine   CoordBase_SetupGH
  @date      Sun Sept 22 2002
  @author    Gabrielle Allen
  @desc
             Allocates CoordBase's GH extension structure
  @enddesc

  @calls     Util_HashCreate

  @var       config
  @vdesc     the CCTK configuration as provided by the flesh
  @vtype     tFleshConfig *
  @vio       unused
  @endvar
  @var       conv_level
  @vdesc     the convergence level
  @vtype     int
  @vio       unused
  @endvar
  @var       GH
  @vdesc     Pointer to CCTK grid hierarchy
  @vtype     cGH *
  @vio       in
  @endvar

  @returntype void *
  @returndesc
              pointer to the allocated GH extension structure
  @endreturndesc
@@*/

static void *CoordBase_SetupGH(tFleshConfig *config, int conv_level, cGH *GH) {
  int maxdim;
  int *default_coord_systems;
  coordbaseGH *myGH;

  /* suppress compiler warnings about unused variables */
  (void)(config + 0);
  (void)(conv_level + 0);
  (void)(GH + 0);

  maxdim = CCTK_MaxDim();

  /* allocate the GH extension and its components */
  myGH = malloc(sizeof(coordbaseGH));
  default_coord_systems = malloc(maxdim * sizeof(int));
  if (myGH && default_coord_systems) {
    myGH->default_coord_systems = default_coord_systems;
    while (--maxdim >= 0) {
      myGH->default_coord_systems[maxdim] = COORDERROR_NOSYSTEM;
    }

    myGH->coordsystems = Util_HashCreate(8);
  } else {
    CCTK_WARN(0, "CoordBase_SetupGH: Unable to allocate memory for GH "
                 "extension");
    free(default_coord_systems);
    myGH = NULL;
  }

  return (myGH);
}
