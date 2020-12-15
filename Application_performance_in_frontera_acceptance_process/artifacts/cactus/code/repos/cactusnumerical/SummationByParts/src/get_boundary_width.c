#include <assert.h>
#include <stdio.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#include "util_ErrorCodes.h"
#include "util_Table.h"

#include "stencil.h"


void Get_Bound_Width ( const CCTK_POINTER_TO_CONST cctkGH_, CCTK_INT *bsize,
                       const CCTK_INT table_handle )
{
  cGH const * restrict const cctkGH = cctkGH_;
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS

  int onesided[6];
  CCTK_INT loc_order; 
  CCTK_INT gsize[6];
  int nelements;

  SBP_determine_onesided_stencil (cctkGH, onesided);

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

  gsize[0] = cctk_nghostzones[0];
  gsize[1] = cctk_nghostzones[0];
  gsize[2] = cctk_nghostzones[1];
  gsize[3] = cctk_nghostzones[1];
  gsize[4] = cctk_nghostzones[2];
  gsize[5] = cctk_nghostzones[2];

  if ( CCTK_Equals(norm_type,"Diagonal") ) {
    if ( sbp_1st_deriv ) {
      switch(loc_order) {
      case 2: {
        for (int d = 0; d<6; d++) bsize[d] = (onesided[d])?1:gsize[d];
        break;
      }
      case 4: {
        for (int d = 0; d<6; d++) bsize[d] = (onesided[d])?4:gsize[d];
        break;
      }
      case 6: {
        for (int d = 0; d<6; d++) bsize[d] = (onesided[d])?6:gsize[d];
        break;
      }
      case 8: {
        for (int d = 0; d<6; d++) bsize[d] = (onesided[d])?8:gsize[d];
        break;
      }
      default:
        CCTK_WARN (0, "Unknown stencil specified");
      }
    } else {
     switch(loc_order) {
      case 2: {
        for (int d = 0; d<6; d++) bsize[d] = (onesided[d])?1:gsize[d];
        break;
      }
      case 4: {
        for (int d = 0; d<6; d++) bsize[d] = (onesided[d])?2:gsize[d];
        break;
      }
      case 6: {
        for (int d = 0; d<6; d++) bsize[d] = (onesided[d])?3:gsize[d];
        break;
      }
      case 8: {
        for (int d = 0; d<6; d++) bsize[d] = (onesided[d])?4:gsize[d];
        break;
      }
      default:
        CCTK_WARN (0, "Unknown 1st derivative stencil specified");
      }
    }
  } else {
    switch(loc_order) {
    case 4: {
      for (int d = 0; d<6; d++) bsize[d] = (onesided[d])?5:gsize[d];
      break;
    }
    case 6: {
      for (int d = 0; d<6; d++) bsize[d] = (onesided[d])?7:gsize[d];
      break;
    }
    default:
      CCTK_WARN (0, "Unknown stencil specified");
    }
  }
} 
