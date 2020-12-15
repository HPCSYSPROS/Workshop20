/*@@
  @file      ParamCheck.c
  @date      Thu Oct  7 17:11:44 1999
  @author    Tom Goodale
  @desc
  C version of Gab's paramcheck stuff
  @enddesc
@@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusBase_CartGrid3D_ParamCheck_c)

void ParamCheck_CartGrid3D(CCTK_ARGUMENTS);

/*@@
  @routine    ParamCheckCartGrid3D
  @date      Tue Feb 23 1999
  @author    Gabrielle Allen
  @desc
  Check parameters for CartGrid3D
  @enddesc
  @calls
  @calledby
  @history
  @hdate Thu Oct  7 17:23:15 1999 @hauthor Tom Goodale
  @hdesc Converted to C
  @endhistory

@@*/
void ParamCheck_CartGrid3D(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  int iflag;

  iflag = 0;

  if (CCTK_Equals(type, "byrange")) {
    if (CCTK_Equals(domain, "octant")) {
      iflag++;
    } else if (CCTK_Equals(domain, "quadrant")) {
      iflag++;
    } else if (CCTK_Equals(domain, "quadrant_reflect_rotate")) {
      iflag++;
    } else if (CCTK_Equals(domain, "bitant")) {
      iflag++;
    } else if (CCTK_Equals(domain, "bitant_rotate")) {
      iflag++;
    } else if (CCTK_Equals(domain, "full")) {
      iflag++;
    }

  } else if (CCTK_Equals(type, "byspacing")) {
    if (CCTK_Equals(domain, "bitant")) {
      iflag++;
    }
    if (CCTK_Equals(domain, "bitant_rotate")) {
      iflag++;
    } else if (CCTK_Equals(domain, "quadrant")) {
      iflag++;
    } else if (CCTK_Equals(domain, "quadrant_reflect_rotate")) {
      iflag++;
    } else if (CCTK_Equals(domain, "octant")) {
      iflag++;
    } else if (CCTK_Equals(domain, "full")) {
      iflag++;
    }

  } else if (CCTK_Equals(type, "coordbase") ||
             CCTK_Equals(type, "multipatch")) {
    if (CCTK_IsFunctionAliased("GetDomainSpecification")) {
      if (CCTK_Equals(domain, "bitant")) {
        iflag++;
      }
      if (CCTK_Equals(domain, "bitant_rotate")) {
        iflag++;
      } else if (CCTK_Equals(domain, "quadrant")) {
        iflag++;
      } else if (CCTK_Equals(domain, "quadrant_reflect_rotate")) {
        iflag++;
      } else if (CCTK_Equals(domain, "octant")) {
        iflag++;
      } else if (CCTK_Equals(domain, "full")) {
        iflag++;
      }
    }
  } else if (CCTK_Equals(type, "box")) {
    iflag++;

    if (!CCTK_Equals(domain, "full"))
      CCTK_PARAMWARN("No symmetries can be used with box grid");
  }

  /* No grid was set up */

  if (iflag != 1) {
    CCTK_PARAMWARN("No grid set up in CartGrid3D");
  }

  if (CCTK_Equals(domain, "bitant_rotate")) {
    if (CCTK_nProcs(cctkGH) != 1)
      CCTK_PARAMWARN("domain 'bitant_rotate' only works on a single processor");

    if (CCTK_Equals(bitant_plane, "xy") && CCTK_Equals(rotation_axis, "z"))
      CCTK_PARAMWARN(
          "rotation_axis=\"z\" is incompatible with bitant_plane=\"xy\"");

    if (CCTK_Equals(bitant_plane, "xz") && CCTK_Equals(rotation_axis, "y"))
      CCTK_PARAMWARN(
          "rotation_axis=\"y\" is incompatible with bitant_plane=\"xz\"");

    if (CCTK_Equals(bitant_plane, "yz") && CCTK_Equals(rotation_axis, "x"))
      CCTK_PARAMWARN(
          "rotation_axis=\"x\" is incompatible with bitant_plane=\"yz\"");
  }

  if (CCTK_Equals(domain, "quadrant_reflect_rotate")) {
    if (CCTK_nProcs(cctkGH) != 1)
      CCTK_PARAMWARN(
          "domain 'quadrant_reflect_rotate' only works on a single processor");

    if (CCTK_Equals(quadrant_direction, "x") && CCTK_Equals(rotation_axis, "x"))
      CCTK_PARAMWARN(
          "rotation_axis=\"x\" is incompatible with quadrant_direction=\"x\"");

    if (CCTK_Equals(quadrant_direction, "y") && CCTK_Equals(rotation_axis, "y"))
      CCTK_PARAMWARN(
          "rotation_axis=\"y\" is incompatible with quadrant_direction=\"y\"");

    if (CCTK_Equals(quadrant_direction, "z") && CCTK_Equals(rotation_axis, "z"))
      CCTK_PARAMWARN(
          "rotation_axis=\"z\" is incompatible with quadrant_direction=\"z\"");
  }

  return;
}
