#include <assert.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#include "stencil.h"



void DiffGf ( const CCTK_POINTER_TO_CONST cctkGH_, const CCTK_INT dir,
              const char *var_name, const char *dvar_name )
{
  cGH const * restrict const cctkGH = cctkGH_;
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS

  CCTK_REAL *var, *dvar;
  CCTK_INT ni, nj, nk, gsize;
  CCTK_REAL delta;
  CCTK_INT bb[2];
  int onesided[6];
  void CCTK_FCALL CCTK_FNAME(deriv_gf_2_1)(const CCTK_REAL *var,
                                       const CCTK_INT *ni,
                                       const CCTK_INT *nj,
                                       const CCTK_INT *nk,
                                       const CCTK_INT *dir,
                                       const CCTK_INT *bb,
                                       const CCTK_INT *gsize,
                                       const CCTK_REAL *delta,
                                       CCTK_REAL *dvar);
  void CCTK_FCALL CCTK_FNAME(deriv_gf_4_2)(const CCTK_REAL *var,
                                       const CCTK_INT *ni,
                                       const CCTK_INT *nj,
                                       const CCTK_INT *nk,
                                       const CCTK_INT *dir,
                                       const CCTK_INT *bb,
                                       const CCTK_INT *gsize,
                                       const CCTK_REAL *delta,
                                       CCTK_REAL *dvar);
  void CCTK_FCALL CCTK_FNAME(deriv_gf_4_3)(const CCTK_REAL *var,
                                       const CCTK_INT *ni,
                                       const CCTK_INT *nj,
                                       const CCTK_INT *nk,
                                       const CCTK_INT *dir,
                                       const CCTK_INT *bb,
                                       const CCTK_INT *gsize,
                                       const CCTK_REAL *delta,
                                       CCTK_REAL *dvar);
  void CCTK_FCALL CCTK_FNAME(deriv_gf_6_3)(const CCTK_REAL *var,
                                       const CCTK_INT *ni,
                                       const CCTK_INT *nj,
                                       const CCTK_INT *nk,
                                       const CCTK_INT *dir,
                                       const CCTK_INT *bb,
                                       const CCTK_INT *gsize,
                                       const CCTK_REAL *delta,
                                       CCTK_REAL *dvar);
  void CCTK_FCALL CCTK_FNAME(deriv_gf_8_4)(const CCTK_REAL *var,
                                       const CCTK_INT *ni,
                                       const CCTK_INT *nj,
                                       const CCTK_INT *nk,
                                       const CCTK_INT *dir,
                                       const CCTK_INT *bb,
                                       const CCTK_INT *gsize,
                                       const CCTK_REAL *delta,
                                       CCTK_REAL *dvar);

  
  ni = cctk_lsh[0]; nj = cctk_lsh[1]; nk = cctk_lsh[2];
  
  SBP_determine_onesided_stencil (cctkGH, onesided);

  switch(dir) {
  case 0: {
    delta = CCTK_DELTA_SPACE(0);
    bb[0] = onesided[0]; bb[1] = onesided[1];
    gsize = cctk_nghostzones[0];
    break;
  }
  case 1: {
    delta = CCTK_DELTA_SPACE(1);
    bb[0] = onesided[2]; bb[1] = onesided[3];
    gsize = cctk_nghostzones[1];
    break;
  }
  case 2: {
    delta = CCTK_DELTA_SPACE(2);
    bb[0] = onesided[4]; bb[1] = onesided[5];
    gsize = cctk_nghostzones[2];
    break;
  }
  default:
      assert(0);
  }

  var = (CCTK_REAL *)(CCTK_VarDataPtr(cctkGH,0,var_name));
  dvar = (CCTK_REAL *)(CCTK_VarDataPtr(cctkGH,0,dvar_name));

  if ( CCTK_Equals(norm_type,"Diagonal") ) {
    switch(order) {
    case 2: {
      CCTK_FNAME(deriv_gf_2_1)(var,&ni,&nj,&nk,&dir,bb,&gsize,&delta,dvar);
      break;
    }
    case 4: {
      CCTK_FNAME(deriv_gf_4_2)(var,&ni,&nj,&nk,&dir,bb,&gsize,&delta,dvar);
      break;
    }
    case 6: {
      CCTK_FNAME(deriv_gf_6_3)(var,&ni,&nj,&nk,&dir,bb,&gsize,&delta,dvar);
      break;
    }
    case 8: {
      CCTK_FNAME(deriv_gf_8_4)(var,&ni,&nj,&nk,&dir,bb,&gsize,&delta,dvar);
      break;
    }
    default:
      assert(0);
    }
  } else {
    switch(order) {
    case 4: {
      CCTK_FNAME(deriv_gf_4_3)(var,&ni,&nj,&nk,&dir,bb,&gsize,&delta,dvar);
      break;
    }
    default:
      assert(0);
    }
  }
} 
