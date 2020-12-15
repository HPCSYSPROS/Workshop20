#include <assert.h>
#include <stdio.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#include "util_ErrorCodes.h"
#include "util_Table.h"

#include "stencil.h"



void DiffCoeff2 ( const CCTK_POINTER_TO_CONST cctkGH_, const CCTK_INT dir,
                 const CCTK_INT nsize, CCTK_INT *imin, CCTK_INT *imax,
                 CCTK_REAL *q, const CCTK_INT table_handle  )
{
  cGH const * restrict const cctkGH = cctkGH_;
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS

  CCTK_INT ni, nj, nk, gsize, loc_order;
  CCTK_INT bb[2];
  int onesided[6];
  int nelements;

  void CCTK_FCALL CCTK_FNAME(set_coeff2_2_1)(const CCTK_INT *nsize,
                                             const CCTK_INT *loc_order,
                                             const CCTK_INT *bb,
                                             const CCTK_INT *gsize,
                                             CCTK_INT *imin,
                                             CCTK_INT *imax,
                                             CCTK_REAL *q);
  void CCTK_FCALL CCTK_FNAME(set_coeff2_2)(const CCTK_INT *nsize,
                                           const CCTK_INT *loc_order,
                                           const CCTK_INT *bb,
                                           const CCTK_INT *gsize,
                                           CCTK_INT *imin,
                                           CCTK_INT *imax,
                                           CCTK_REAL *q);
  void CCTK_FCALL CCTK_FNAME(set_coeff2_4_2)(const CCTK_INT *nsize,
                                             const CCTK_INT *loc_order,
                                             const CCTK_INT *bb,
                                             const CCTK_INT *gsize,
                                             CCTK_INT *imin,
                                             CCTK_INT *imax,
                                             CCTK_REAL *q);
  void CCTK_FCALL CCTK_FNAME(set_coeff2_4_2_opt)(const CCTK_INT *nsize,
                                             const CCTK_INT *loc_order,
                                             const CCTK_INT *bb,
                                             const CCTK_INT *gsize,
                                             CCTK_INT *imin,
                                             CCTK_INT *imax,
                                             CCTK_REAL *q);
  void CCTK_FCALL CCTK_FNAME(set_coeff2_4)(const CCTK_INT *nsize,
                                           const CCTK_INT *loc_order,
                                           const CCTK_INT *bb,
                                           const CCTK_INT *gsize,
                                           CCTK_INT *imin,
                                           CCTK_INT *imax,
                                           CCTK_REAL *q);
  void CCTK_FCALL CCTK_FNAME(set_coeff2_6_3)(const CCTK_INT *nsize,
                                             const CCTK_INT *loc_order,
                                             const CCTK_INT *bb,
                                             const CCTK_INT *gsize,
                                             CCTK_INT *imin,
                                             CCTK_INT *imax,
                                             CCTK_REAL *q);
  void CCTK_FCALL CCTK_FNAME(set_coeff2_6)(const CCTK_INT *nsize,
                                           const CCTK_INT *loc_order,
                                           const CCTK_INT *bb,
                                           const CCTK_INT *gsize,
                                           CCTK_INT *imin,
                                           CCTK_INT *imax,
                                           CCTK_REAL *q);
  void CCTK_FCALL CCTK_FNAME(set_coeff2_8_4)(const CCTK_INT *nsize,
                                             const CCTK_INT *loc_order,
                                             const CCTK_INT *bb,
                                             const CCTK_INT *gsize,
                                             CCTK_INT *imin,
                                             CCTK_INT *imax,
                                             CCTK_REAL *q);
  void CCTK_FCALL CCTK_FNAME(set_coeff2_8)(const CCTK_INT *nsize,
                                           const CCTK_INT *loc_order,
                                           const CCTK_INT *bb,
                                           const CCTK_INT *gsize,
                                           CCTK_INT *imin,
                                           CCTK_INT *imax,
                                           CCTK_REAL *q);

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

  SBP_determine_onesided_stencil (cctkGH, onesided);

  switch(dir) {
  case 0: {
    assert(nsize==ni);
    bb[0] = onesided[0]; bb[1] = onesided[1];
    gsize = cctk_nghostzones[0];
    break;
  }
  case 1: {
    assert(nsize==nj);
    bb[0] = onesided[2]; bb[1] = onesided[3];
    gsize = cctk_nghostzones[1];
    break;
  }
  case 2: {
    assert(nsize==nk);
    bb[0] = onesided[4]; bb[1] = onesided[5];
    gsize = cctk_nghostzones[2];
    break;
  }
  default:
    CCTK_WARN (0, "Wrong direction specified");
  }

  if ( CCTK_Equals(norm_type,"Diagonal") ) {
    if ( sbp_2nd_deriv ) {
      if ( CCTK_Equals(operator_type,"Minimal Bandwidth") ) {
        switch(loc_order) {
        case 2: {
          CCTK_FNAME(set_coeff2_2_1)(&nsize,&loc_order,bb,&gsize,imin,imax,q);
          break;
        }
        case 4: {
          CCTK_FNAME(set_coeff2_4_2)(&nsize,&loc_order,bb,&gsize,imin,imax,q);
          break;
        }
        case 6: {
          CCTK_FNAME(set_coeff2_6_3)(&nsize,&loc_order,bb,&gsize,imin,imax,q);
          break;
        }
        case 8: {
          CCTK_FNAME(set_coeff2_8_4)(&nsize,&loc_order,bb,&gsize,imin,imax,q);
          break;
        }
        default:
          CCTK_WARN (0, "Unknown stencil specified");
        }
      } else {
        switch(loc_order) {
        case 2: {
          CCTK_FNAME(set_coeff2_2_1)(&nsize,&loc_order,bb,&gsize,imin,imax,q);
          break;
        }
        case 4: {
          CCTK_FNAME(set_coeff2_4_2_opt)(&nsize,&loc_order,bb,&gsize,imin,imax,q);
          break;
        }
        case 6: {
          CCTK_FNAME(set_coeff2_6_3)(&nsize,&loc_order,bb,&gsize,imin,imax,q);
          break;
        }
        case 8: {
          CCTK_FNAME(set_coeff2_8_4)(&nsize,&loc_order,bb,&gsize,imin,imax,q);
          break;
        }
        default:
          CCTK_WARN (0, "Unknown stencil specified");
        }
      }
    } else {
      switch(loc_order) {
      case 2: {
        CCTK_FNAME(set_coeff2_2)(&nsize,&loc_order,bb,&gsize,imin,imax,q);
        break;
      }
      case 4: {
        CCTK_FNAME(set_coeff2_4)(&nsize,&loc_order,bb,&gsize,imin,imax,q);
        break;
      }
      case 6: {
        CCTK_FNAME(set_coeff2_6)(&nsize,&loc_order,bb,&gsize,imin,imax,q);
        break;
      }
      case 8: {
        CCTK_FNAME(set_coeff2_8)(&nsize,&loc_order,bb,&gsize,imin,imax,q);
        break;
      }
      default:
        CCTK_WARN (0, "Unknown 2nd derivative stencil specified");
      }
    }
  } else {
    CCTK_WARN (0, "Second derivatives not implemented for restricted full norm");
  }
}
