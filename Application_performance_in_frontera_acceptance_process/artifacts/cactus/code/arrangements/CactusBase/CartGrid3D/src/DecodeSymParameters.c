
/*@@
  @file    DecodeSymParameters.c
  @date    Wed May 10 18:58:00 EST 2000
  @author  Erik Schnetter
  @desc
           Decode the symmetry parameters.
  @enddesc
  @version $Id$
@@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Symmetry.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusBase_CartGrid3D_DecodeSymParameters_c)

void DecodeSymParameters3D(int sym[6]);
void CCTK_FCALL CCTK_FNAME(DecodeSymParameters3D)(int sym[6]);

/*@@
   @routine    DecodeSymParameters3D
   @date       Thu May 11 11:49:08 2000
   @author     Erik Schnetter
   @desc
      Decode the Symmetry parameters.
      returns the symmetry flags (yes/no=1/0)
      in the array sym
   @enddesc
@@*/

void DecodeSymParameters3D(int sym[6]) {
  DECLARE_CCTK_PARAMETERS

  /* The default is as set by the explicit symmetry parameters */
  /* lower faces */
  sym[0] = symmetry_xmin;
  sym[2] = symmetry_ymin;
  sym[4] = symmetry_zmin;

  /* upper faces */
  sym[1] = symmetry_xmax;
  sym[3] = symmetry_ymax;
  sym[5] = symmetry_zmax;

  /* The default can be overridden by bitant, quadrant, and octant mode */
  if (CCTK_Equals(domain, "bitant")) {
    if (CCTK_Equals(bitant_plane, "xy")) {
      sym[4] = GFSYM_REFLECTION;
    } else if (CCTK_Equals(bitant_plane, "xz")) {
      sym[2] = GFSYM_REFLECTION;
    } else if (CCTK_Equals(bitant_plane, "yz")) {
      sym[0] = GFSYM_REFLECTION;
    }
  } else if (CCTK_Equals(domain, "bitant_rotate")) {
    if (CCTK_Equals(bitant_plane, "xy")) {
      if (CCTK_Equals(rotation_axis, "y"))
        sym[4] = GFSYM_ROTATION_Y;
      else if (CCTK_Equals(rotation_axis, "x"))
        sym[4] = GFSYM_ROTATION_X;
    } else if (CCTK_Equals(bitant_plane, "xz")) {
      if (CCTK_Equals(rotation_axis, "x"))
        sym[2] = GFSYM_ROTATION_X;
      else if (CCTK_Equals(rotation_axis, "z"))
        sym[2] = GFSYM_ROTATION_Z;
    } else if (CCTK_Equals(bitant_plane, "yz")) {
      if (CCTK_Equals(rotation_axis, "y"))
        sym[0] = GFSYM_ROTATION_Y;
      else if (CCTK_Equals(rotation_axis, "z"))
        sym[0] = GFSYM_ROTATION_Z;
    }
  } else if (CCTK_Equals(domain, "quadrant")) {
    if (CCTK_Equals(quadrant_direction, "x")) {
      sym[2] = GFSYM_REFLECTION;
      sym[4] = GFSYM_REFLECTION;
    } else if (CCTK_Equals(quadrant_direction, "y")) {
      sym[0] = GFSYM_REFLECTION;
      sym[4] = GFSYM_REFLECTION;
    } else if (CCTK_Equals(quadrant_direction, "z")) {
      sym[0] = GFSYM_REFLECTION;
      sym[2] = GFSYM_REFLECTION;
    }
  } else if (CCTK_Equals(domain, "quadrant_reflect_rotate")) {
    if (CCTK_Equals(quadrant_direction, "x")) {
      if (CCTK_Equals(rotation_axis, "y")) {
        sym[2] = GFSYM_REFLECTION;
        sym[4] = GFSYM_ROTATION_Y;
      } else if (CCTK_Equals(rotation_axis, "z")) {
        sym[2] = GFSYM_ROTATION_Z;
        sym[4] = GFSYM_REFLECTION;
      }
    } else if (CCTK_Equals(quadrant_direction, "y")) {
      if (CCTK_Equals(rotation_axis, "x")) {
        sym[0] = GFSYM_REFLECTION;
        sym[4] = GFSYM_ROTATION_X;
      }
      if (CCTK_Equals(rotation_axis, "z")) {
        sym[0] = GFSYM_ROTATION_Z;
        sym[4] = GFSYM_REFLECTION;
      }
    } else if (CCTK_Equals(quadrant_direction, "z")) {
      if (CCTK_Equals(rotation_axis, "x")) {
        sym[0] = GFSYM_REFLECTION;
        sym[2] = GFSYM_ROTATION_X;
      }
      if (CCTK_Equals(rotation_axis, "y")) {
        sym[0] = GFSYM_ROTATION_Y;
        sym[2] = GFSYM_REFLECTION;
      }
    }
  } else if (CCTK_Equals(domain, "octant")) {
    sym[0] = GFSYM_REFLECTION;
    sym[2] = GFSYM_REFLECTION;
    sym[4] = GFSYM_REFLECTION;
  }
}

void CCTK_FCALL CCTK_FNAME(DecodeSymParameters3D)(int sym[6]) {
  DecodeSymParameters3D(sym);
}
