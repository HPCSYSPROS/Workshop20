/*@@
  @file      CartGrid3D.c
  @date      Thu Oct  7 13:20:06 1999
  @author    Tom Goodale
  @desc
             Set up coordinates for a 3D Cartesian grid.
             C Conversion of Fortran routine written by Gab.
  @version   $Id$
  @enddesc
@@*/

#include <stdio.h>
#include <math.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "util_Table.h"

#include "Symmetry.h"
#include "CoordBase.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusBase_CartGrid3D_CartGrid3D_c);

/********************************************************************
 *********************  Macro Definitions  **************************
 ********************************************************************/
#define max(a, b) ((a) > (b) ? (a) : (b))
#define SQR(a) ((a) * (a))

/********************************************************************
 *********************  Scheduled Routine Prototypes  ***************
 ********************************************************************/
void CartGrid3D_SetRanges(CCTK_ARGUMENTS);
void CartGrid3D_SetCoordinates(CCTK_ARGUMENTS);

/********************************************************************
 *********************  External Routine Prototypes  ****************
 ********************************************************************/
void DecodeSymParameters3D(int sym[6]);

/********************************************************************
 *********************  Local Routine Prototypes  *******************
 ********************************************************************/

/********************************************************************
 *********************  Scheduled Routines  *************************
 ********************************************************************/
/*@@
  @routine    CartGrid3D_SetRanges
  @date       Oct 1999?
  @author     Tom Goodale?  Gabrielle Allen? Thomas Radke
  @date       Mon 3 Jan 2005
  @author     Thomas Radke
  @desc
              Sets up ranges for Cartesian coordinates.
  @enddesc
  @calls      DecodeSymParameters3D, CCTK_Equals, CCTK_WARN,
              CCTK_CoordRegisterRange, CCTK_CoordRegisterRangePhysIndex,
              CCTK_INFO, Coord_CoordHandle, Util_TableSetReal,
              Util_TableSetString, Util_TableSetInt

  @var        CCTK_ARGUMENTS
  @vdesc      Cactus argument list
  @vtype
  @vio        in/out
  @endvar
@@*/
void CartGrid3D_SetRanges(CCTK_ARGUMENTS) {
  int i;
  int coord_handle, ierr;
  CCTK_REAL this_delta[3], origin[3], min1[3], max1[3];
  CCTK_REAL *coarse_delta[3];
  double lower[3], upper[3];
  int domainsym[6], cntstag[3], loweri[3], upperi[3], do_periodic[3];
  char coord_name[16];
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_EQUALS(set_coordinate_ranges_on, "first level")) {
    /* Ranges must be set up only once, and this must happen on the
       coarse grid.  However, the coarse grid itself may not actually
       exist; in this case, use the coarsest existing grid.  We assume
       that this is the first grid for which this routine is
       called.  */
    static int is_coarsest_refinement_level = 1;
    if (!is_coarsest_refinement_level) {
      return;
    }
    is_coarsest_refinement_level = 0;
  } else {
    /* Ranges need to be set up only once (or once per map), on the
       coarsest refinement level */
    int const is_coarsest_refinement_level =
        cctk_levfac[0] == 1 && cctk_levfac[1] == 1 && cctk_levfac[2] == 1;
    if (!is_coarsest_refinement_level) {
      return;
    }
  }

  coarse_delta[0] = coarse_dx;
  coarse_delta[1] = coarse_dy;
  coarse_delta[2] = coarse_dz;

  /* Calculate the coordinate ranges only for the coarsest level */
  /* Avoid origin?  Default is yes */
  cntstag[0] = no_origin && no_originx && avoid_origin && avoid_originx;
  cntstag[1] = no_origin && no_originy && avoid_origin && avoid_originy;
  cntstag[2] = no_origin && no_originz && avoid_origin && avoid_originz;

  /* Determine symmetries of domain */
  DecodeSymParameters3D(domainsym);

  do_periodic[0] = periodic && periodic_x;
  do_periodic[1] = periodic && periodic_y;
  do_periodic[2] = periodic && periodic_z;

  /* Calculate physical indices, using symmetries and periodicity */
  for (i = 0; i < 3; i++) {
    loweri[i] = 0;
    upperi[i] = cctk_gsh[i] - 1;

    if (domainsym[2 * i + 0] || do_periodic[i]) {
      loweri[i] += cctk_nghostzones[i];
    }
    if (domainsym[2 * i + 1] || do_periodic[i]) {
      upperi[i] -= cctk_nghostzones[i];
    }
  }

  /****************************************************************
   *
   *  BYRANGE
   *
   *  User gives: minimum and maximum values of coordinates and
   *              the number of gridpoints on the coarse grid
   *
   ***************************************************************/
  /**************************************************************
   *
   *   BOX (-0.5 to 0.5)
   *
   *   User gives: number of gridpoints on the coarse grid
   *
   **************************************************************/

  if (CCTK_Equals(type, "byrange") || CCTK_Equals(type, "box")) {
    if (CCTK_Equals(type, "box")) {
      /*  Coordinates are all -0.5 to 0.5 */
      min1[0] = min1[1] = min1[2] = -0.5;
      max1[0] = max1[1] = max1[2] = 0.5;
    } else {
      if (xyzmin != -424242) {
        min1[0] = min1[1] = min1[2] = xyzmin;
      } else {
        min1[0] = xmin;
        min1[1] = ymin;
        min1[2] = zmin;
      }

      if (xyzmax != -424242) {
        max1[0] = max1[1] = max1[2] = xyzmax;
      } else {
        max1[0] = xmax;
        max1[1] = ymax;
        max1[2] = zmax;
      }
    }

    /* Grid spacing on coarsest grid */
    for (i = 0; i < 3; i++) {
      if (domainsym[2 * i + 0]) {
        if (cntstag[i]) {
          *coarse_delta[i] =
              max1[i] / (cctk_gsh[i] - cctk_nghostzones[i] - 0.5);
          origin[i] = -(cctk_nghostzones[i] - 0.5) * *coarse_delta[i];
        } else {
          *coarse_delta[i] = max1[i] / (cctk_gsh[i] - cctk_nghostzones[i] - 1);
          origin[i] = -cctk_nghostzones[i] * *coarse_delta[i];
        }
      } else if (domainsym[2 * i + 1]) {
        if (cntstag[i]) {
          *coarse_delta[i] =
              fabs(min1[i]) / (cctk_gsh[i] - cctk_nghostzones[i] - 0.5);
        } else {
          *coarse_delta[i] =
              fabs(min1[i]) / (cctk_gsh[i] - cctk_nghostzones[i] - 1);
        }
        origin[i] = min1[i];
      } else {
        if (cntstag[i]) {
          CCTK_VWarn(4, __LINE__, __FILE__, CCTK_THORNSTRING,
                     "Ignoring request to avoid origin in %c-direction, "
                     "it is not relevant for this grid type",
                     'x' + i);
        }
        *coarse_delta[i] = (max1[i] - min1[i]) / max(cctk_gsh[i] - 1, 1);
        origin[i] = min1[i];
      }

      this_delta[i] = *coarse_delta[i];
    }
  }

  /**************************************************************
   *  BYSPACING
   *
   *  User gives: grid spacing on the coarsest GH and
   *              the number of gridpoints on the coarsest GH
   *
   **************************************************************/
  else if (CCTK_Equals(type, "byspacing")) {
    /* Dx, Dy, Dx on the coarsest grid */
    if (dxyz > 0) {
      *coarse_delta[0] = *coarse_delta[1] = *coarse_delta[2] = dxyz;
    } else {
      *coarse_delta[0] = dx;
      *coarse_delta[1] = dy;
      *coarse_delta[2] = dz;
    }

    for (i = 0; i < 3; i++) {
      this_delta[i] = *coarse_delta[i];

      /* Set minimum values of coordinates */
      if (domainsym[2 * i + 0]) {
        origin[i] = -(cctk_nghostzones[i] - cntstag[i] * 0.5);
      } else if (domainsym[2 * i + 1]) {
        origin[i] = -(cctk_gsh[i] - 1 - cctk_nghostzones[i] + cntstag[i] * 0.5);
      } else {
        origin[i] = -0.5 * (cctk_gsh[i] - 1 - cntstag[i] * cctk_gsh[i] % 2);
      }
      origin[i] *= this_delta[i];
    }
  }

  /**************************************************************
   *  COORDBASE
   *
   *  CoordBase gives: grid spacing on the coarsest GH and
   *                   minimum and maximum values of coordinates and
   *                   the number of gridpoints on the coarsest GH
   *
   **************************************************************/
  else if (CCTK_Equals(type, "coordbase")) {
    CCTK_REAL physical_min[3];
    CCTK_REAL physical_max[3];
    CCTK_REAL interior_min[3];
    CCTK_REAL interior_max[3];
    CCTK_REAL exterior_min[3];
    CCTK_REAL exterior_max[3];
    CCTK_REAL spacing[3];
    int d;

    ierr = GetDomainSpecification(3, physical_min, physical_max, interior_min,
                                  interior_max, exterior_min, exterior_max,
                                  spacing);
    if (ierr)
      CCTK_WARN(0, "error returned from function GetDomainSpecification");

    /* Adapt to convergence level */
    for (d = 0; d < 3; ++d) {
      spacing[d] *= pow(cctkGH->cctk_convfac, cctkGH->cctk_convlevel);
    }

    ierr = ConvertFromPhysicalBoundary(3, physical_min, physical_max,
                                       interior_min, interior_max, exterior_min,
                                       exterior_max, spacing);
    if (ierr)
      CCTK_WARN(0, "error returned from function ConvertFromPhysicalBoundary");

    for (d = 0; d < 3; ++d) {
      origin[d] = exterior_min[d];
      this_delta[d] = spacing[d];
      *coarse_delta[d] = this_delta[d];
    }
  }
  /**************************************************************
   *  MULTIPATCH
   *
   *  MultiPatch gives: grid spacing on the coarsest GH and
   *                    minimum and maximum values of coordinates and
   *                    the number of gridpoints on the coarsest GH
   *
   **************************************************************/
  else if (CCTK_Equals(type, "multipatch")) {
    CCTK_REAL physical_min[3];
    CCTK_REAL physical_max[3];
    CCTK_REAL interior_min[3];
    CCTK_REAL interior_max[3];
    CCTK_REAL exterior_min[3];
    CCTK_REAL exterior_max[3];
    CCTK_REAL spacing[3];
    CCTK_INT map;
    int d;

    map = MultiPatch_GetMap(cctkGH);
    if (map < 0)
      CCTK_WARN(0, "error returned from function MultiPatch_GetMap");

    ierr = MultiPatch_GetDomainSpecification(
        map, 3, physical_min, physical_max, interior_min, interior_max,
        exterior_min, exterior_max, spacing);
    if (ierr)
      CCTK_WARN(
          0, "error returned from function MultiPatch_GetDomainSpecification");

    /* Adapt to convergence level */
    for (d = 0; d < 3; ++d) {
      spacing[d] *= pow(cctkGH->cctk_convfac, cctkGH->cctk_convlevel);
    }

    if (CCTK_IsFunctionAliased("MultiPatch_ConvertFromPhysicalBoundary")) {
      ierr = MultiPatch_ConvertFromPhysicalBoundary(
          map, 3, physical_min, physical_max, interior_min, interior_max,
          exterior_min, exterior_max, spacing);
      if (ierr)
        CCTK_WARN(0, "error returned from function "
                     "MultiPatch_ConvertFromPhysicalBoundary");
    } else {
      ierr = ConvertFromPhysicalBoundary(3, physical_min, physical_max,
                                         interior_min, interior_max,
                                         exterior_min, exterior_max, spacing);
      if (ierr)
        CCTK_WARN(0,
                  "error returned from function ConvertFromPhysicalBoundary");
    }

    for (d = 0; d < 3; ++d) {
      origin[d] = exterior_min[d];
      this_delta[d] = spacing[d];
      *coarse_delta[d] = this_delta[d];
    }
  }

  else {
    if (0)
      CCTK_WARN(0, "type is out of bounds");
  }

  /* Register the coordinate ranges */
  for (i = 0; i < 3; i++) {
    cctkGH->cctk_delta_space[i] = this_delta[i];
    cctkGH->cctk_origin_space[i] = origin[i];

    lower[i] = origin[i];
    upper[i] = origin[i] + this_delta[i] * (cctk_gsh[i] - 1);
    if (CCTK_CoordRegisterRange(cctkGH, lower[i], upper[i], i + 1, NULL,
                                "cart3d") < 0) {
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Failed to register %c-coordinate computational range",
                 'x' + i);
    }
    if (CCTK_CoordRegisterRangePhysIndex(cctkGH, loweri[i], upperi[i], i + 1,
                                         NULL, "cart3d") < 0) {
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Failed to register %c-coordinate physical range", 'x' + i);
    }
  }

  CCTK_INFO("Grid Spacings:");
  CCTK_VInfo(CCTK_THORNSTRING, "dx=>%12.7e  dy=>%12.7e  dz=>%12.7e",
             (double)cctk_delta_space[0], (double)cctk_delta_space[1],
             (double)cctk_delta_space[2]);
  CCTK_INFO("Computational Coordinates:");
  CCTK_VInfo(CCTK_THORNSTRING,
             "x=>[%6.3f,%6.3f]  y=>[%6.3f,%6.3f]  z=>[%6.3f,%6.3f]", lower[0],
             upper[0], lower[1], upper[1], lower[2], upper[2]);
  CCTK_INFO("Indices of Physical Coordinates:");
  CCTK_VInfo(CCTK_THORNSTRING, "x=>[%d,%d]  y=>[%d,%d]  z=>[%d,%d]", loweri[0],
             upperi[0], loweri[1], upperi[1], loweri[2], upperi[2]);

  if ((domainsym[0] == GFSYM_ROTATION_Y || domainsym[2] == GFSYM_ROTATION_X) &&
      (lower[2] + upper[2] > 1e-12)) {
    CCTK_WARN(0, "minimum z must equal maximum z for rotation symmetry");
  }

  if ((domainsym[0] == GFSYM_ROTATION_Z || domainsym[4] == GFSYM_ROTATION_X) &&
      (lower[1] + upper[1] > 1e-12)) {
    CCTK_WARN(0, "minimum y must equal maximum y for rotation symmetry");
  }

  if ((domainsym[2] == GFSYM_ROTATION_Z || domainsym[4] == GFSYM_ROTATION_Y) &&
      (lower[0] - upper[0] > 1e-12)) {
    CCTK_WARN(0, "minimum x must equal maximum x for rotation symmetry");
  }

  /* cart3d */
  for (i = 0; i < 3; i++) {
    sprintf(coord_name, "%c", 'x' + i);

    coord_handle = Coord_CoordHandle(cctkGH, coord_name, "cart3d");
    if (coord_handle < 0) {
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Error retreiving coordinate handle for '%s' of cart3d",
                 coord_name);
    }
    sprintf(coord_name, "grid::%c", 'x' + i);
    ierr = Util_TableSetInt(coord_handle, loweri[i], "PHYSICALMIN");
    ierr += Util_TableSetReal(coord_handle, lower[i], "COMPMIN");
    ierr += Util_TableSetInt(coord_handle, upperi[i], "PHYSICALMAX");
    ierr += Util_TableSetReal(coord_handle, upper[i], "COMPMAX");
    ierr += Util_TableSetString(coord_handle, "uniform", "TYPE");
    ierr += Util_TableSetString(coord_handle, "no", "TIMEDEPENDENT");
    ierr += Util_TableSetString(coord_handle, "CCTK_REAL", "DATATYPE");
    ierr +=
        Util_TableSetInt(coord_handle, CCTK_VarIndex(coord_name), "GAINDEX");
    ierr += Util_TableSetReal(coord_handle, cctk_delta_space[i], "DELTA");
  }

  /* cart2d */
  for (i = 0; i < 2; i++) {
    sprintf(coord_name, "%c", 'x' + i);
    coord_handle = Coord_CoordHandle(cctkGH, coord_name, "cart2d");
    if (coord_handle < 0) {
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Error retreiving coordinate handle for '%s' of cart2d",
                 coord_name);
    }
    sprintf(coord_name, "grid::%c", 'x' + i);
    ierr = Util_TableSetReal(coord_handle, lower[i], "PHYSICALMIN"); /*??*/
    ierr += Util_TableSetReal(coord_handle, lower[i], "COMPMIN");
    ierr += Util_TableSetReal(coord_handle, upper[i], "PHYSICALMAX"); /*??*/
    ierr += Util_TableSetReal(coord_handle, upper[i], "COMPMAX");
    ierr += Util_TableSetString(coord_handle, "uniform", "TYPE");
    ierr += Util_TableSetString(coord_handle, "no", "TIMEDEPENDENT");
    ierr += Util_TableSetString(coord_handle, "CCTK_REAL", "DATATYPE");
    ierr +=
        Util_TableSetInt(coord_handle, CCTK_VarIndex(coord_name), "GAINDEX");
    ierr += Util_TableSetReal(coord_handle, cctk_delta_space[i], "DELTA");
  }

  /* cart1d */
  coord_handle = Coord_CoordHandle(cctkGH, "x", "cart1d");
  if (coord_handle < 0) {
    CCTK_WARN(0, "Error retreiving coordinate handle for x of cart1d");
  }
  ierr = Util_TableSetReal(coord_handle, lower[0], "PHYSICALMIN"); /*??*/
  ierr += Util_TableSetReal(coord_handle, lower[0], "COMPMIN");
  ierr += Util_TableSetReal(coord_handle, upper[0], "PHYSICALMAX"); /*??*/
  ierr += Util_TableSetReal(coord_handle, upper[0], "COMPMAX");
  ierr += Util_TableSetString(coord_handle, "uniform", "TYPE");
  ierr += Util_TableSetString(coord_handle, "no", "TIMEDEPENDENT");
  ierr += Util_TableSetString(coord_handle, "CCTK_REAL", "DATATYPE");
  ierr += Util_TableSetInt(coord_handle, CCTK_VarIndex("grid::x"), "GAINDEX");
  ierr += Util_TableSetReal(coord_handle, cctk_delta_space[0], "DELTA");

  /* Set up coordinate tables */
  /* Should this be done in a function?
     WriteCoordinateTable(cctkGH, "cart3d"); */
}

/*@@
  @routine    CartGrid3D_SetCoordinates
  @date       2004-06-17
  @author     Christian Ott
  @desc
              Sets up Cartesian coordinates.
  @enddesc
  @var        CCTK_ARGUMENTS
  @vdesc      Cactus argument list
  @vtype
  @vio        in/out
  @endvar
@@*/
void CartGrid3D_SetCoordinates(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;

  /* CCTK_VInfo(CCTK_THORNSTRING,"Resetting coordinates after regridding."); */

  for (int k = 0; k < cctk_lsh[2]; k++) {
    for (int j = 0; j < cctk_lsh[1]; j++) {
      for (int i = 0; i < cctk_lsh[0]; i++) {
        int idx = CCTK_GFINDEX3D(cctkGH, i, j, k);
        x[idx] =
            CCTK_DELTA_SPACE(0) * (i + cctk_lbnd[0]) + CCTK_ORIGIN_SPACE(0);
        y[idx] =
            CCTK_DELTA_SPACE(1) * (j + cctk_lbnd[1]) + CCTK_ORIGIN_SPACE(1);
        z[idx] =
            CCTK_DELTA_SPACE(2) * (k + cctk_lbnd[2]) + CCTK_ORIGIN_SPACE(2);
        r[idx] = sqrt(SQR(x[idx]) + SQR(y[idx]) + SQR(z[idx]));
      }
    }
  }
}
