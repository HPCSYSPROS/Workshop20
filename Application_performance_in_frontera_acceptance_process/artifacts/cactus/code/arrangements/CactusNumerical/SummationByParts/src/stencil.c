#include "cctk.h"
#include "cctk_Parameters.h"

#include "util_ErrorCodes.h"
#include "util_Table.h"

#include "stencil.h"



/* Determine whether a boundary with the symmetry handle symbnd is a
   "regular" symmetry boundary (where the SBP stencils should not be
   modified), or an outer boundary (or a multi-patch boundary), where
   the SBP stencils need to be modified.  */

void SBP_determine_onesided_stencil (const cGH * cctkGH, int * onesided)
{
  DECLARE_CCTK_PARAMETERS;
  
  int symtable;
  int mp_sym_handle;
  CCTK_INT symbnd[6];
  CCTK_INT bbox[6];
  int ierr;
  int d;
  
  symtable = SymmetryTableHandleForGrid (cctkGH);
  if (symtable<0) {
    CCTK_WARN(0,"Cannot get symmetry table handle -- maybe thorn SymBase is not active?");
  }
  
  ierr = Util_TableGetIntArray (symtable, 6, symbnd, "symmetry_handle");
  if (ierr!=6) {
    CCTK_WARN(0,"Cannot get symmetry handles");
  }

  mp_sym_handle = SymmetryHandleOfName ( "multipatch" );
  
  if (mp_sym_handle >= 0) {
    /* We are using a multi-patch system */
    ierr = MultiPatch_GetBbox (cctkGH, 6, bbox);
    if (ierr) {
      CCTK_WARN(0,"Cannot get multi-patch bbox information");
    }
  }
  
  for (d=0; d<6; ++d) {
    if (! cctkGH->cctk_bbox[d]) {
      /* This is an inter-processor boundary */
      onesided[d] = 0;
    } else {
      /* On an outer boundary (which is not a symmetry boundary), it
         is either symbnd < 0, or it is a boundary claimed by the MP
         infrastructure and bbox != 0 */
      if (symbnd[d] < 0 ||
          (mp_sym_handle >= 0 && symbnd[d] == mp_sym_handle && bbox[d]))
      {
        /* Use one-sided stencils near outer boundaries if the user
           wants it so */
        onesided[d] = onesided_outer_boundaries;
      } else {
        /* If the symmetry boundary is a multi-patch boundary, then
           symbnd = mp_sym_handle.  However, we can only check this if
           a multi-patch thorn is active, i.e., when mp_sym_handle >=
           0 */
        if (mp_sym_handle >= 0 && symbnd[d] == mp_sym_handle) {
          /* Use one-sided stencils near inter-patch boundaries if the
             user wants it so */
          onesided[d] = onesided_interpatch_boundaries;
        } else {
          /* Always use centred stencils near regular symmetry
             boundaries (e.g. a reflection symmetry) */
          onesided[d] = 0;
        }
      }
    }
  }
}

CCTK_FCALL void CCTK_FNAME (SBP_determine_onesided_stencil) (CCTK_POINTER_TO_CONST * cctkGH, int * onesided)
{
  SBP_determine_onesided_stencil (* cctkGH, onesided);
}
