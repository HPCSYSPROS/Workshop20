/*@@
  @file      Startup.c
  @date      Mon Mar 15 15:48:42 1999
  @author    Gerd Lanfermann
  @desc
             Startup file to register the GHextension and coordinates
  @enddesc
  @version   $Id$
@@*/

#include <stdlib.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "util_Table.h"

#include "CoordBase.h"
#include "Symmetry.h"

/* the rcs ID and its dummy function to use it */
static const char *rcsid = "$Header$";
CCTK_FILEVERSION(CactusBase_CartGrid3D_Startup_c);

/********************************************************************
 ********************    External Prototypes   **********************
 ********************************************************************/
int SymmetryStartup(void);
void RegisterCartGrid3DCoords(CCTK_ARGUMENTS);

/********************************************************************
 ********************    Internal Prototypes   **********************
 ********************************************************************/
static void *SetupGH(tFleshConfig *config, int convlevel, cGH *GH);

/********************************************************************
 ********************    External Routines   ************************
 ********************************************************************/
/*@@
  @routine    SymmetryStartup
  @date       Mon Mar 15 15:49:16 1999
  @author     Gerd Lanfermann
  @desc
              Routine registers the GH extension for CartGrid3D
              along with its setup routine.
  @enddesc
  @calls      CCTK_RegisterGHExtension
              CCTK_RegisterGHExtensionSetupGH
  @history
  @endhistory

  @returntype void
@@*/
int SymmetryStartup(void) {
  CCTK_RegisterGHExtensionSetupGH(CCTK_RegisterGHExtension("Symmetry"),
                                  SetupGH);
  return 0;
}

/*@@
  @routine    RegisterCartGrid3DCoords
  @date
  @author     Gabrielle Allen
  @desc
              Routine registers the coordinates provided by
              CartGrid3D, using both the new and old APIs.  The old
              CCTK_ API is deprecated.
  @enddesc
  @calls      CCTK_CoordRegisterSystem
              CCTK_CoordRegisterData

  @returntype int
  @returndesc
              0 for success, or negative in case of errors
  @endreturndesc
@@*/
void RegisterCartGrid3DCoords(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int ierr, coord_system_handle;

  /* Register coordinate systems */
  ierr = Coord_SystemRegister(cctkGH, 3, "cart3d");
  ierr += Coord_SystemRegister(cctkGH, 2, "cart2d");
  ierr += Coord_SystemRegister(cctkGH, 1, "cart1d");
  if (ierr < 0) {
    CCTK_WARN(0, "Error registering cartnd coordinate systems");
  } else {

    /* Register coordinates for cart3d */
    coord_system_handle = Coord_SystemHandle(cctkGH, "cart3d");
    if (coord_system_handle < 0) {
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Error obtaining system handle for cart3d");
    }
    ierr = Coord_CoordRegister(cctkGH, coord_system_handle, 1, "x");
    ierr += Coord_CoordRegister(cctkGH, coord_system_handle, 2, "y");
    ierr += Coord_CoordRegister(cctkGH, coord_system_handle, 3, "z");
    if (ierr < 0) {
      CCTK_WARN(0, "Error registering cart3d coordinates");
    }

    /* Fill out rest of coordinate system table for cart3d */
    ierr = Util_TableSetString(coord_system_handle, "uniform", "TYPE");
    if (ierr < 0) {
      CCTK_WARN(1, "Error registering cart3d type");
    }

    /* Register coordinates for cart2d */
    coord_system_handle = Coord_SystemHandle(cctkGH, "cart2d");
    if (coord_system_handle < 0) {
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Error obtaining system handle for cart2d");
    }
    ierr = Coord_CoordRegister(cctkGH, coord_system_handle, 1, "x");
    ierr += Coord_CoordRegister(cctkGH, coord_system_handle, 2, "y");
    if (ierr < 0) {
      CCTK_WARN(0, "Error registering cart2d coordinates");
    }

    /* Fill out rest of coordinate system table for cart2d */
    ierr = Util_TableSetString(coord_system_handle, "uniform", "TYPE");
    if (ierr < 0) {
      CCTK_WARN(1, "Error registering cart2d type");
    }

    /* Register coordinate for cart1d */
    coord_system_handle = Coord_SystemHandle(cctkGH, "cart1d");
    if (coord_system_handle < 0) {
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Error obtaining system handle for cart1d");
    }
    ierr = Coord_CoordRegister(cctkGH, coord_system_handle, 1, "x");
    if (ierr < 0) {
      CCTK_WARN(0, "Error registering cart1d coordinate");
    }

    /* Fill out rest of coordinate system table for cart1d */
    ierr = Util_TableSetString(coord_system_handle, "uniform", "TYPE");
    if (ierr < 0) {
      CCTK_WARN(1, "Error registering cart1d type");
    }

    /* Register cartnd as the default coordinate systems */
    if (register_default_coordinate_systems) {
      ierr = Coord_SetDefaultSystem(cctkGH, "cart3d");
      ierr += Coord_SetDefaultSystem(cctkGH, "cart2d");
      ierr += Coord_SetDefaultSystem(cctkGH, "cart1d");
      if (ierr < 0) {
        CCTK_WARN(1, "Error registering cartnd as default coordinate systems");
      }
    }
  }

  /* Register coordinates under the old API */
  CCTK_CoordRegisterSystem(3, "cart3d");
  CCTK_CoordRegisterSystem(3, "spher3d");

  if (CCTK_CoordRegisterData(1, "grid::x", "x", "cart3d") < 0) {
    CCTK_WARN(1, "Problem with registering coordinate x");
  }
  if (CCTK_CoordRegisterData(2, "grid::y", "y", "cart3d") < 0) {
    CCTK_WARN(1, "Problem with registering coordinate y");
  }
  if (CCTK_CoordRegisterData(3, "grid::z", "z", "cart3d") < 0) {
    CCTK_WARN(1, "Problem with registering coordinate z");
  }
  if (CCTK_CoordRegisterData(1, "grid::r", "r", "spher3d") < 0) {
    CCTK_WARN(1, "Problem with registering coordinate r");
  }
}

/********************************************************************
 ********************    Internal Routines   ************************
 ********************************************************************/
static void *SetupGH(tFleshConfig *config, int convlevel, cGH *GH) {
  int i, j, maxdim, numvars;
  SymmetryGHex *myGH;

  /* avoid compiler warnings about unused arguments */
  (void)(config + 0);
  (void)(convlevel + 0);
  (void)(GH + 0);

  maxdim = CCTK_MaxDim();
  numvars = CCTK_NumVars();

  /* allocate the GH extension */
  myGH = (SymmetryGHex *)malloc(sizeof(SymmetryGHex));
  if (myGH) {
    /* allocation for the number of grid functions */
    myGH->GFSym = (int **)malloc(numvars * sizeof(int *));

    /* allocation for the number of dimensions*/
    for (i = 0; i < numvars; i++) {
      myGH->GFSym[i] = (int *)malloc(2 * maxdim * sizeof(int));

      for (j = 0; j < 2 * maxdim; j++) {
        myGH->GFSym[i][j] = GFSYM_UNSET; /* not set */
      }
    }
  }

  return (myGH);
}
