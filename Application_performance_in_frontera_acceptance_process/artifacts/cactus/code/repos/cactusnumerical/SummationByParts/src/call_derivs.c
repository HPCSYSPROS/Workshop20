#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#include "util_ErrorCodes.h"
#include "util_Table.h"

#include "stencil.h"


void DiffGv ( const CCTK_POINTER_TO_CONST cctkGH_, const CCTK_INT dir,
              const CCTK_REAL *var, CCTK_REAL *dvar,
              const CCTK_INT table_handle )
{
  cGH const * restrict const cctkGH = cctkGH_;
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS

  CCTK_INT i, ni, nj, nk, gsize, loc_order;
  CCTK_REAL delta;
  CCTK_INT bb[2], loc_offset[2];
  int onesided[6];
  CCTK_INT enforce_centered[6];
  CCTK_INT offset[6];
  CCTK_INT dummy[6];
  int nelements;

  /*
  CCTK_INT nboundaryzones[6];
  CCTK_INT is_internal[6];
  CCTK_INT is_staggered[6];
  */

  void CCTK_FCALL CCTK_FNAME(deriv_gf_2_1)(const CCTK_REAL *var,
                                       const CCTK_INT *ni,
                                       const CCTK_INT *nj,
                                       const CCTK_INT *nk,
                                       const CCTK_INT *dir,
                                       const CCTK_INT *bb,
                                       const CCTK_INT *gsize,
                                       const CCTK_INT *offset,
                                       const CCTK_REAL *delta,
                                       CCTK_REAL *dvar);
  void CCTK_FCALL CCTK_FNAME(deriv_gf_4_2)(const CCTK_REAL *var,
                                       const CCTK_INT *ni,
                                       const CCTK_INT *nj,
                                       const CCTK_INT *nk,
                                       const CCTK_INT *dir,
                                       const CCTK_INT *bb,
                                       const CCTK_INT *gsize,
                                       const CCTK_INT *offset,
                                       const CCTK_REAL *delta,
                                       CCTK_REAL *dvar);
  void CCTK_FCALL CCTK_FNAME(deriv_gf_4_3)(const CCTK_REAL *var,
                                       const CCTK_INT *ni,
                                       const CCTK_INT *nj,
                                       const CCTK_INT *nk,
                                       const CCTK_INT *dir,
                                       const CCTK_INT *bb,
                                       const CCTK_INT *gsize,
                                       const CCTK_INT *offset,
                                       const CCTK_REAL *delta,
                                       CCTK_REAL *dvar);
  void CCTK_FCALL CCTK_FNAME(deriv_gf_4_3_opt)(const CCTK_REAL *var,
                                       const CCTK_INT *ni,
                                       const CCTK_INT *nj,
                                       const CCTK_INT *nk,
                                       const CCTK_INT *dir,
                                       const CCTK_INT *bb,
                                       const CCTK_INT *gsize,
                                       const CCTK_INT *offset,
                                       const CCTK_REAL *delta,
                                       CCTK_REAL *dvar);
  void CCTK_FCALL CCTK_FNAME(deriv_gf_6_3)(const CCTK_REAL *var,
                                       const CCTK_INT *ni,
                                       const CCTK_INT *nj,
                                       const CCTK_INT *nk,
                                       const CCTK_INT *dir,
                                       const CCTK_INT *bb,
                                       const CCTK_INT *gsize,
                                       const CCTK_INT *offset,
                                       const CCTK_REAL *delta,
                                       CCTK_REAL *dvar);
  void CCTK_FCALL CCTK_FNAME(deriv_gf_6_3_opt)(const CCTK_REAL *var,
                                       const CCTK_INT *ni,
                                       const CCTK_INT *nj,
                                       const CCTK_INT *nk,
                                       const CCTK_INT *dir,
                                       const CCTK_INT *bb,
                                       const CCTK_INT *gsize,
                                       const CCTK_INT *offset,
                                       const CCTK_REAL *delta,
                                       CCTK_REAL *dvar);
  void CCTK_FCALL CCTK_FNAME(deriv_gf_6_5)(const CCTK_REAL *var,
                                       const CCTK_INT *ni,
                                       const CCTK_INT *nj,
                                       const CCTK_INT *nk,
                                       const CCTK_INT *dir,
                                       const CCTK_INT *bb,
                                       const CCTK_INT *gsize,
                                       const CCTK_INT *offset,
                                       const CCTK_REAL *delta,
                                       CCTK_REAL *dvar);
  void CCTK_FCALL CCTK_FNAME(deriv_gf_6_5_opt)(const CCTK_REAL *var,
                                       const CCTK_INT *ni,
                                       const CCTK_INT *nj,
                                       const CCTK_INT *nk,
                                       const CCTK_INT *dir,
                                       const CCTK_INT *bb,
                                       const CCTK_INT *gsize,
                                       const CCTK_INT *offset,
                                       const CCTK_REAL *delta,
                                       CCTK_REAL *dvar);
  void CCTK_FCALL CCTK_FNAME(deriv_gf_8_4)(const CCTK_REAL *var,
                                       const CCTK_INT *ni,
                                       const CCTK_INT *nj,
                                       const CCTK_INT *nk,
                                       const CCTK_INT *dir,
                                       const CCTK_INT *bb,
                                       const CCTK_INT *gsize,
                                       const CCTK_INT *offset,
                                       const CCTK_REAL *delta,
                                       CCTK_REAL *dvar);
  void CCTK_FCALL CCTK_FNAME(deriv_gf_8_4_opt)(const CCTK_REAL *var,
                                       const CCTK_INT *ni,
                                       const CCTK_INT *nj,
                                       const CCTK_INT *nk,
                                       const CCTK_INT *dir,
                                       const CCTK_INT *bb,
                                       const CCTK_INT *gsize,
                                       const CCTK_INT *offset,
                                       const CCTK_REAL *delta,
                                       CCTK_REAL *dvar);
  
  void get_shiftout ( const CCTK_POINTER_TO_CONST cctkGH_, CCTK_INT *offset);

  void CCTK_FCALL CCTK_FNAME(SBP_Poisoning) (
                                       const CCTK_INT *ni, 
                                       const CCTK_INT *nj, 
                                       const CCTK_INT *nk, 
                                       const CCTK_INT *bb, 
                                       const CCTK_INT *offset, 
                                       CCTK_REAL *dvar);

  ni = cctk_lsh[0]; nj = cctk_lsh[1]; nk = cctk_lsh[2];
 
  if ( table_handle >=0 ) {
    nelements = Util_TableGetInt ( table_handle, &loc_order, "order" );
    if ( nelements == UTIL_ERROR_TABLE_NO_SUCH_KEY ) {
      loc_order = order;
    } else if ( nelements != 1) {
      CCTK_WARN (0, "The options table has an entry \"order\", but it does not have the right properties");
    }
  } else {
    loc_order = order;
  }

  for (i=0; i<6; i++) {
    enforce_centered[i] = 0;
  }

  if ( table_handle >=0 ) {
    nelements = Util_TableGetIntArray ( table_handle, 6, dummy, "enforce_centered" );
    if ( nelements == UTIL_ERROR_TABLE_NO_SUCH_KEY ) {
    } else if ( nelements != 6) {
      CCTK_WARN (0, "The options table has an entry \"enforce_centered\", but it does not have 6 elements");
    } else {
      for (i=0; i<6; i++) {
        enforce_centered[i] = dummy[i];
      }
    }
  }

  /* get values of boundary_shiftout_* from CoordBase  */
  for (i=0; i<6; i++) {
    offset[i] = 0;
  }
  get_shiftout (cctkGH_, offset);
  /*GetBoundarySpecification (6, nboundaryzones, is_internal, is_staggered, offset);*/

  SBP_determine_onesided_stencil (cctkGH, onesided);

  switch(dir) {
  case 0: {
    delta = CCTK_DELTA_SPACE(0);
    if ( enforce_centered[0] ) {
      bb[0] = 0;
    } else {
      bb[0] = onesided[0];
    }
    loc_offset[0] = offset[0];
    if ( enforce_centered[1] ) {
      bb[1] = 0;
    } else {
      bb[1] = onesided[1];
    }
    loc_offset[1] = offset[1];
    gsize = cctk_nghostzones[0];
    break;
  }
  case 1: {
    delta = CCTK_DELTA_SPACE(1);
    if ( enforce_centered[2] ) {
      bb[0] = 0;
    } else {
      bb[0] = onesided[2];
    }
    loc_offset[0] = offset[2];
    if ( enforce_centered[3] ) {
      bb[0] = 0;
    } else {
      bb[1] = onesided[3];
    }
    loc_offset[1] = offset[3];
    gsize = cctk_nghostzones[1];
    break;
  }
  case 2: {
    delta = CCTK_DELTA_SPACE(2);
    if ( enforce_centered[4] ) {
      bb[0] = 0;
    } else {
      bb[0] = onesided[4];
    }
    loc_offset[0] = offset[4];
    if ( enforce_centered[5] ) {
      bb[0] = 0;
    } else {
      bb[1] = onesided[5];
    }
    loc_offset[1] = offset[5];
    gsize = cctk_nghostzones[2];
    break;
  }
  default:
    CCTK_WARN (0, "Wrong direction specified");
  }

  if ( CCTK_Equals(norm_type,"Diagonal") ) {
    if ( CCTK_Equals(operator_type,"Minimal Bandwidth") ) {
      switch(loc_order) {
      case 2: {
        CCTK_FNAME(deriv_gf_2_1)(var,&ni,&nj,&nk,&dir,bb,&gsize,loc_offset,&delta,dvar);
        break;
      }
      case 4: {
        CCTK_FNAME(deriv_gf_4_2)(var,&ni,&nj,&nk,&dir,bb,&gsize,loc_offset,&delta,dvar);
        break;
      }
      case 6: {
        CCTK_FNAME(deriv_gf_6_3)(var,&ni,&nj,&nk,&dir,bb,&gsize,loc_offset,&delta,dvar);
        break;
      }
      case 8: {
        CCTK_FNAME(deriv_gf_8_4)(var,&ni,&nj,&nk,&dir,bb,&gsize,loc_offset,&delta,dvar);
        break;
      }
      default:
        CCTK_WARN (0, "Unknown stencil specified");
      }
    } else {
      switch(loc_order) {
      case 2: {
        CCTK_FNAME(deriv_gf_2_1)(var,&ni,&nj,&nk,&dir,bb,&gsize,loc_offset,&delta,dvar);
        break;
      }
      case 4: {
        CCTK_FNAME(deriv_gf_4_2)(var,&ni,&nj,&nk,&dir,bb,&gsize,loc_offset,&delta,dvar);
        break;
      }
      case 6: {
        CCTK_FNAME(deriv_gf_6_3_opt)(var,&ni,&nj,&nk,&dir,bb,&gsize,loc_offset,&delta,dvar);
        break;
      }
      case 8: {
        CCTK_FNAME(deriv_gf_8_4_opt)(var,&ni,&nj,&nk,&dir,bb,&gsize,loc_offset,&delta,dvar);
        break;
      }
      default:
        CCTK_WARN (0, "Unknown stencil specified");
      }
    }
  } else {
    if ( CCTK_Equals(operator_type,"Minimal Bandwidth") ) {
      switch(loc_order) {
      case 4: {
        CCTK_FNAME(deriv_gf_4_3)(var,&ni,&nj,&nk,&dir,bb,&gsize,loc_offset,&delta,dvar);
        break;
      }
      case 6: {
        CCTK_FNAME(deriv_gf_6_5)(var,&ni,&nj,&nk,&dir,bb,&gsize,loc_offset,&delta,dvar);
        break;
      }
      default:
        CCTK_WARN (0, "Unknown stencil specified");
      }
    } else {
      switch(loc_order) {
      case 4: {
        CCTK_FNAME(deriv_gf_4_3_opt)(var,&ni,&nj,&nk,&dir,bb,&gsize,loc_offset,&delta,dvar);
        break;
      }
      case 6: {
        CCTK_FNAME(deriv_gf_6_5_opt)(var,&ni,&nj,&nk,&dir,bb,&gsize,loc_offset,&delta,dvar);
        break;
      }
      default:
        CCTK_WARN (0, "Unknown stencil specified");
      }
    }
  }

  if (poison_derivatives) {
      CCTK_FNAME(SBP_Poisoning)(&ni,&nj,&nk,bb,offset,dvar);
  }
} 
