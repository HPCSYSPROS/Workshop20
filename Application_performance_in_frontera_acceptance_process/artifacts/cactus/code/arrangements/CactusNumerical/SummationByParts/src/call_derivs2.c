#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#include "util_ErrorCodes.h"
#include "util_Table.h"

#include "stencil.h"



void DiffGv2 ( const CCTK_POINTER_TO_CONST cctkGH_, const CCTK_INT dir1,
               const CCTK_INT dir2, const CCTK_REAL *var, CCTK_REAL *dvar2,
               const CCTK_INT table_handle )
{
  cGH const * restrict const cctkGH = cctkGH_;
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS

  CCTK_INT ni, nj, nk, gsize, loc_order, dir;
  CCTK_REAL delta;
  CCTK_INT bb[2];
  int onesided[6];
  int nelements;

  void CCTK_FCALL CCTK_FNAME(deriv2_gf_2_1)(const CCTK_REAL *var,
                                       const CCTK_INT *ni,
                                       const CCTK_INT *nj,
                                       const CCTK_INT *nk,
                                       const CCTK_INT *dir,
                                       const CCTK_INT *bb,
                                       const CCTK_INT *gsize,
                                       const CCTK_REAL *delta,
                                       CCTK_REAL *dvar2);
  void CCTK_FCALL CCTK_FNAME(deriv2_gf_4_2)(const CCTK_REAL *var,
                                       const CCTK_INT *ni,
                                       const CCTK_INT *nj,
                                       const CCTK_INT *nk,
                                       const CCTK_INT *dir,
                                       const CCTK_INT *bb,
                                       const CCTK_INT *gsize,
                                       const CCTK_REAL *delta,
                                       CCTK_REAL *dvar2);
  void CCTK_FCALL CCTK_FNAME(deriv2_gf_4_2_opt)(const CCTK_REAL *var,
                                       const CCTK_INT *ni,
                                       const CCTK_INT *nj,
                                       const CCTK_INT *nk,
                                       const CCTK_INT *dir,
                                       const CCTK_INT *bb,
                                       const CCTK_INT *gsize,
                                       const CCTK_REAL *delta,
                                       CCTK_REAL *dvar2);
  void CCTK_FCALL CCTK_FNAME(deriv2_gf_6_3)(const CCTK_REAL *var,
                                       const CCTK_INT *ni,
                                       const CCTK_INT *nj,
                                       const CCTK_INT *nk,
                                       const CCTK_INT *dir,
                                       const CCTK_INT *bb,
                                       const CCTK_INT *gsize,
                                       const CCTK_REAL *delta,
                                       CCTK_REAL *dvar2);
  void CCTK_FCALL CCTK_FNAME(deriv2_gf_8_4)(const CCTK_REAL *var,
                                       const CCTK_INT *ni,
                                       const CCTK_INT *nj,
                                       const CCTK_INT *nk,
                                       const CCTK_INT *dir,
                                       const CCTK_INT *bb,
                                       const CCTK_INT *gsize,
                                       const CCTK_REAL *delta,
                                       CCTK_REAL *dvar2);
  void CCTK_FCALL CCTK_FNAME(deriv2_mixed)(const CCTK_POINTER_TO_CONST *cctkGH,
                                       const CCTK_INT *dir1,
                                       const CCTK_INT *dir2,
                                       const CCTK_REAL *var,
                                       const CCTK_INT *ni,
                                       const CCTK_INT *nj,
                                       const CCTK_INT *nk,
                                       CCTK_REAL *dvar2,
                                       const CCTK_INT *table_handle);

  
  ni = cctk_lsh[0]; nj = cctk_lsh[1]; nk = cctk_lsh[2];
 
  if ( dir1 != dir2 ) {
    CCTK_FNAME(deriv2_mixed)(&cctkGH_,&dir1,&dir2,var,&ni,&nj,&nk,dvar2,&table_handle);
    return;
  } else {
    dir = dir1;
  }

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

  SBP_determine_onesided_stencil (cctkGH, onesided);

  switch(dir) {
  case 0: {
    delta = CCTK_DELTA_SPACE(0);
    bb[0] = onesided[0];
    bb[1] = onesided[1];
    gsize = cctk_nghostzones[0];
    break;
  }
  case 1: {
    delta = CCTK_DELTA_SPACE(1);
    bb[0] = onesided[2];
    bb[1] = onesided[3];
    gsize = cctk_nghostzones[1];
    break;
  }
  case 2: {
    delta = CCTK_DELTA_SPACE(2);
    bb[0] = onesided[4];
    bb[1] = onesided[5];
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
        CCTK_FNAME(deriv2_gf_2_1)(var,&ni,&nj,&nk,&dir,bb,&gsize,&delta,dvar2);
        break;
      }
      case 4: {
        CCTK_FNAME(deriv2_gf_4_2)(var,&ni,&nj,&nk,&dir,bb,&gsize,&delta,dvar2);
        break;
      }
      case 6: {
        CCTK_FNAME(deriv2_gf_6_3)(var,&ni,&nj,&nk,&dir,bb,&gsize,&delta,dvar2);
        break;
      }
      case 8: {
        CCTK_FNAME(deriv2_gf_8_4)(var,&ni,&nj,&nk,&dir,bb,&gsize,&delta,dvar2);
        break;
      }
      default:
      CCTK_WARN (0, "Unknown 2nd derivative stencil specified");
      }
    } else {
      switch(loc_order) {
      case 2: {
        CCTK_FNAME(deriv2_gf_2_1)(var,&ni,&nj,&nk,&dir,bb,&gsize,&delta,dvar2);
        break;
      }
      case 4: {
        CCTK_FNAME(deriv2_gf_4_2_opt)(var,&ni,&nj,&nk,&dir,bb,&gsize,&delta,dvar2);
        break;
      }
      case 6: {
        CCTK_FNAME(deriv2_gf_6_3)(var,&ni,&nj,&nk,&dir,bb,&gsize,&delta,dvar2);
        break;
      }
      case 8: {
        CCTK_FNAME(deriv2_gf_8_4)(var,&ni,&nj,&nk,&dir,bb,&gsize,&delta,dvar2);
        break;
      }
      default:
        CCTK_WARN (0, "Unknown stencil specified");
      }
    }
  } else {
    switch(loc_order) {
    default:
      CCTK_WARN (0, "Unknown stencil specified");
    }
  }
} 
