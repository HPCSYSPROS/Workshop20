/* $Header$ */

#include <assert.h>
#include <stdlib.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

#include "util_ErrorCodes.h"
#include "util_Table.h"

#include "stencil.h"

void CCTK_FCALL CCTK_FNAME(dissipation_4_3_opt) (const CCTK_REAL *var,
                                             const CCTK_INT *lsh,
                                             const CCTK_INT *gsh,
                                             const CCTK_INT *lbnd,
                                             const CCTK_INT *bbox,
                                             const CCTK_INT *gsize,
                                             const CCTK_INT *offset,
                                             const CCTK_REAL *dx,
                                             const CCTK_REAL *epsdis,
                                             const CCTK_REAL *diss_fraction,
                                             const CCTK_INT *npatches,
                                             const CCTK_INT *patch,
                                             CCTK_REAL *rhs);
void CCTK_FCALL CCTK_FNAME(dissipation_4_3_alt) (const CCTK_REAL *var,
                                             const CCTK_INT *lsh,
                                             const CCTK_INT *gsh,
                                             const CCTK_INT *lbnd,
                                             const CCTK_INT *bbox,
                                             const CCTK_INT *gsize,
                                             const CCTK_INT *offset,
                                             const CCTK_REAL *dx,
                                             const CCTK_REAL *epsdis,
                                             const CCTK_REAL *diss_fraction,
                                             const CCTK_INT *npatches,
                                             const CCTK_INT *patch,
                                             CCTK_REAL *rhs);
void CCTK_FCALL CCTK_FNAME(dissipation_6_5_opt) (const CCTK_REAL *var,
                                             const CCTK_INT *ni,
                                             const CCTK_INT *nj,
                                             const CCTK_INT *nk,
                                             const CCTK_INT *bbox,
                                             const CCTK_INT *gsize,
                                             const CCTK_INT *offset,
                                             const CCTK_REAL *dx,
                                             const CCTK_REAL *epsdis,
                                             const CCTK_REAL *diss_fraction,
                                             const CCTK_INT *npatches,
                                             const CCTK_INT *patch,
                                             CCTK_REAL *rhs);
void CCTK_FCALL CCTK_FNAME(dissipation_6_5_alt) (const CCTK_REAL *var,
                                             const CCTK_INT *ni,
                                             const CCTK_INT *nj,
                                             const CCTK_INT *nk,
                                             const CCTK_INT *bbox,
                                             const CCTK_INT *gsize,
                                             const CCTK_INT *offset,
                                             const CCTK_REAL *dx,
                                             const CCTK_REAL *epsdis,
                                             const CCTK_REAL *diss_fraction,
                                             const CCTK_INT *npatches,
                                             const CCTK_INT *patch,
                                             CCTK_REAL *rhs);
void CCTK_FCALL CCTK_FNAME(dissipation_2_1) (const CCTK_REAL *var,
                                             const CCTK_INT *ni,
                                             const CCTK_INT *nj,
                                             const CCTK_INT *nk,
                                             const CCTK_INT *bbox,
                                             const CCTK_INT *gsize,
                                             const CCTK_INT *offset,
                                             const CCTK_REAL *dx,
                                             const CCTK_REAL *epsdis,
                                             CCTK_REAL *rhs);
void CCTK_FCALL CCTK_FNAME(dissipation_2_1_alt) (const CCTK_REAL *var,
                                             const CCTK_INT *ni,
                                             const CCTK_INT *nj,
                                             const CCTK_INT *nk,
                                             const CCTK_INT *bbox,
                                             const CCTK_INT *gsize,
                                             const CCTK_INT *offset,
                                             const CCTK_REAL *dx,
                                             const CCTK_REAL *epsdis,
                                             CCTK_REAL *rhs);
void CCTK_FCALL CCTK_FNAME(dissipation_2_1_delta) (const CCTK_REAL *var,
                                             const CCTK_INT *ni,
                                             const CCTK_INT *nj,
                                             const CCTK_INT *nk,
                                             const CCTK_INT *bbox,
                                             const CCTK_INT *gsize,
                                             const CCTK_INT *offset,
                                             const CCTK_REAL *dx,
                                             const CCTK_REAL *dy,
                                             const CCTK_REAL *dz,
                                             const CCTK_REAL *epsdis,
                                             CCTK_REAL *rhs);
void CCTK_FCALL CCTK_FNAME(dissipation_4_2) (const CCTK_REAL *var,
                                             const CCTK_INT *ni,
                                             const CCTK_INT *nj,
                                             const CCTK_INT *nk,
                                             const CCTK_INT *bbox,
                                             const CCTK_INT *gsize,
                                             const CCTK_INT *offset,
                                             const CCTK_REAL *dx,
                                             const CCTK_REAL *epsdis,
                                             CCTK_REAL *rhs);
void CCTK_FCALL CCTK_FNAME(dissipation_4_2_alt) (const CCTK_REAL *var,
                                             const CCTK_INT *ni,
                                             const CCTK_INT *nj,
                                             const CCTK_INT *nk,
                                             const CCTK_INT *bbox,
                                             const CCTK_INT *gsize,
                                             const CCTK_INT *offset,
                                             const CCTK_REAL *dx,
                                             const CCTK_REAL *epsdis,
                                             CCTK_REAL *rhs);
void CCTK_FCALL CCTK_FNAME(dissipation_4_2_delta) (const CCTK_REAL *var,
                                             const CCTK_INT *ni,
                                             const CCTK_INT *nj,
                                             const CCTK_INT *nk,
                                             const CCTK_INT *bbox,
                                             const CCTK_INT *gsize,
                                             const CCTK_INT *offset,
                                             const CCTK_REAL *dx,
                                             const CCTK_REAL *dy,
                                             const CCTK_REAL *dz,
                                             const CCTK_REAL *epsdis,
                                             CCTK_REAL *rhs);
void CCTK_FCALL CCTK_FNAME(dissipation_6_3) (const CCTK_REAL *var,
                                             const CCTK_INT *ni,
                                             const CCTK_INT *nj,
                                             const CCTK_INT *nk,
                                             const CCTK_INT *bbox,
                                             const CCTK_INT *gsize,
                                             const CCTK_INT *offset,
                                             const CCTK_REAL *dx,
                                             const CCTK_REAL *epsdis,
                                             CCTK_REAL *rhs);
void CCTK_FCALL CCTK_FNAME(dissipation_6_3_alt) (const CCTK_REAL *var,
                                             const CCTK_INT *ni,
                                             const CCTK_INT *nj,
                                             const CCTK_INT *nk,
                                             const CCTK_INT *bbox,
                                             const CCTK_INT *gsize,
                                             const CCTK_INT *offset,
                                             const CCTK_REAL *dx,
                                             const CCTK_REAL *epsdis,
                                             CCTK_REAL *rhs);
void CCTK_FCALL CCTK_FNAME(dissipation_6_3_delta) (const CCTK_REAL *var,
                                             const CCTK_INT *ni,
                                             const CCTK_INT *nj,
                                             const CCTK_INT *nk,
                                             const CCTK_INT *bbox,
                                             const CCTK_INT *gsize,
                                             const CCTK_INT *offset,
                                             const CCTK_REAL *dx,
                                             const CCTK_REAL *dy,
                                             const CCTK_REAL *dz,
                                             const CCTK_REAL *epsdis,
                                             CCTK_REAL *rhs);
void CCTK_FCALL CCTK_FNAME(dissipation_8_4) (const CCTK_REAL *var,
                                             const CCTK_INT *ni,
                                             const CCTK_INT *nj,
                                             const CCTK_INT *nk,
                                             const CCTK_INT *bbox,
                                             const CCTK_INT *gsize,
                                             const CCTK_INT *offset,
                                             const CCTK_REAL *dx,
                                             const CCTK_REAL *epsdis,
                                             CCTK_REAL *rhs);
void CCTK_FCALL CCTK_FNAME(dissipation_8_4_alt) (const CCTK_REAL *var,
                                             const CCTK_INT *ni,
                                             const CCTK_INT *nj,
                                             const CCTK_INT *nk,
                                             const CCTK_INT *bbox,
                                             const CCTK_INT *gsize,
                                             const CCTK_INT *offset,
                                             const CCTK_REAL *dx,
                                             const CCTK_REAL *epsdis,
                                             CCTK_REAL *rhs);
void CCTK_FCALL CCTK_FNAME(dissipation_8_4_delta) (const CCTK_REAL *var,
                                             const CCTK_INT *ni,
                                             const CCTK_INT *nj,
                                             const CCTK_INT *nk,
                                             const CCTK_INT *bbox,
                                             const CCTK_INT *gsize,
                                             const CCTK_INT *offset,
                                             const CCTK_REAL *dx,
                                             const CCTK_REAL *dy,
                                             const CCTK_REAL *dz,
                                             const CCTK_REAL *epsdis,
                                             CCTK_REAL *rhs);
void get_shiftout ( const CCTK_POINTER_TO_CONST cctkGH_, CCTK_INT *offset);

void CCTK_FCALL CCTK_FNAME(SBP_Poisoning) (  const CCTK_INT *ni, 
                                             const CCTK_INT *nj, 
                                             const CCTK_INT *nk, 
                                             const CCTK_INT *bb, 
                                             const CCTK_INT *offset, 
                                             CCTK_REAL *dvar);


static void
apply (int const varindex, char const * const optstring, void * const arg);

void
SBP_DissipationAdd (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_TraverseString (vars, apply, cctkGH, CCTK_GROUP_OR_VAR);
}

void
apply (int const varindex, char const * const optstring, void * const arg)
{
  cGH const * const cctkGH = (cGH const *) arg;
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  int rhsindex;
  int vargroup, rhsgroup;
  cGroup vardata, rhsdata;
  CCTK_REAL const * varptr;
  CCTK_REAL       * rhsptr;
  CCTK_INT gsize[3];
  CCTK_REAL dx[3];
  int d;
  int ierr;
  CCTK_INT bbox[6];
  int onesided[6];
  CCTK_INT npatches, patch;
  CCTK_INT offset[6];
  int i;
  
  assert (varindex >= 0);
  
  for (d=0; d<3; ++d) {
    dx[d] = CCTK_DELTA_SPACE(d)*h_scaling[d];
    gsize[d] = cctk_nghostzones[d];
  }

  /* get values of boundary_shiftout_* from CoordBase  */
  for (i=0; i<6; i++) {
    offset[i] = 0;
  }
  get_shiftout (cctkGH, offset);
  /*GetBoundarySpecification (6, nboundaryzones, is_internal, is_staggered, offset);*/

  SBP_determine_onesided_stencil (cctkGH, onesided);
  for (d=0; d<6; ++d) {
    bbox[d] = onesided[d];
  }
  rhsindex = MoLQueryEvolvedRHS (varindex);
  if (rhsindex < 0) {
    char * const fullvarname = CCTK_FullName (varindex);
    assert (fullvarname);
    CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                "There is no RHS variable registered with MoL for the evolved variable \"%s\"",
                fullvarname);
    free (fullvarname);
  }
  assert (rhsindex >= 0);
  
/*  if (verbose) {
    char * const fullvarname = CCTK_FullName (varindex);
    char * const fullrhsname = CCTK_FullName (rhsindex);
    assert (fullvarname);
    assert (fullrhsname);
    CCTK_VInfo (CCTK_THORNSTRING,
                "Applying dissipation to \"%s\" (RHS \"%s\")",
                fullvarname, fullrhsname);
    free (fullvarname);
    free (fullrhsname);
  } */
  
  vargroup = CCTK_GroupIndexFromVarI (varindex);
  assert (vargroup >= 0);
  rhsgroup = CCTK_GroupIndexFromVarI (rhsindex);
  assert (rhsgroup >= 0);
  
  ierr = CCTK_GroupData (vargroup, &vardata);
  assert (!ierr);
  ierr = CCTK_GroupData (rhsgroup, &rhsdata);
  assert (!ierr);
  
  assert (vardata.grouptype == CCTK_GF);
  assert (vardata.vartype == CCTK_VARIABLE_REAL);
  assert (vardata.dim == cctk_dim);
  assert (rhsdata.grouptype == CCTK_GF);
  assert (rhsdata.vartype == CCTK_VARIABLE_REAL);
  assert (rhsdata.dim == cctk_dim);
  
  varptr = CCTK_VarDataPtrI (cctkGH, 0, varindex);
  assert (varptr);
  rhsptr = CCTK_VarDataPtrI (cctkGH, 0, rhsindex);
  assert (rhsptr);

  CCTK_INT lsh[3], gsh[3], lbnd[3];
  for (int d=0; d<3; ++d) {
    lsh[d] = cctk_lsh[d];
    gsh[d] = cctk_gsh[d];
    lbnd[d] = cctk_lbnd[d];
  }
  
  if ( CCTK_Equals(norm_type,"Diagonal") ) {
    switch(order) {
    case 2: {
      if ( CCTK_Equals(dissipation_type,"Mattson-Svard-Nordstrom") ) {
        CCTK_FNAME(dissipation_2_1)
          (varptr, &lsh[0], &lsh[1], &lsh[2],
                            bbox, gsize, offset, dx, &epsdis, rhsptr);
      } else {
        if (!use_variable_deltas) {
          CCTK_FNAME(dissipation_2_1_alt)
            (varptr, &lsh[0], &lsh[1], &lsh[2],
                              bbox, gsize, offset, dx, &epsdis, rhsptr);
        } else {
          CCTK_FNAME(dissipation_2_1_delta)
            (varptr, &lsh[0], &lsh[1], &lsh[2],
                              bbox, gsize, offset,
                              sbp_dx, sbp_dy, sbp_dz, &epsdis, rhsptr);
        }
      }
      break;
    }
    case 4: {
      if ( CCTK_Equals(dissipation_type,"Mattson-Svard-Nordstrom") ) {
        CCTK_FNAME(dissipation_4_2)
          (varptr, &lsh[0], &lsh[1], &lsh[2],
                            bbox, gsize, offset, dx, &epsdis, rhsptr);
      } else {
        if (!use_variable_deltas) {
          CCTK_FNAME(dissipation_4_2_alt)
            (varptr, &lsh[0], &lsh[1], &lsh[2],
                              bbox, gsize, offset, dx, &epsdis, rhsptr);
        } else {
          CCTK_FNAME(dissipation_4_2_delta)
            (varptr, &lsh[0], &lsh[1], &lsh[2],
                              bbox, gsize, offset,
                              sbp_dx, sbp_dy, sbp_dz, &epsdis, rhsptr);
        }
      }
      break;
    }
    case 6: {
      if ( CCTK_Equals(dissipation_type,"Mattson-Svard-Nordstrom") ) {
        CCTK_FNAME(dissipation_6_3)
          (varptr, &lsh[0], &lsh[1], &lsh[2],
                            bbox, gsize, offset, dx, &epsdis, rhsptr);
      } else {
        if (!use_variable_deltas) {
          CCTK_FNAME(dissipation_6_3_alt)
            (varptr, &lsh[0], &lsh[1], &lsh[2],
                              bbox, gsize, offset, dx, &epsdis, rhsptr);
        } else {
          CCTK_FNAME(dissipation_6_3_delta)
            (varptr, &lsh[0], &lsh[1], &lsh[2],
                              bbox, gsize, offset,
                              sbp_dx, sbp_dy, sbp_dz, &epsdis, rhsptr);
        }
      }
      break;
    }
    case 8: {
      if ( CCTK_Equals(dissipation_type,"Mattson-Svard-Nordstrom") ) {
        CCTK_FNAME(dissipation_8_4)
          (varptr, &lsh[0], &lsh[1], &lsh[2],
                            bbox, gsize, offset, dx, &epsdis, rhsptr);
      } else {
        if (!use_variable_deltas) {
          CCTK_FNAME(dissipation_8_4_alt)
            (varptr, &lsh[0], &lsh[1], &lsh[2],
                              bbox, gsize, offset, dx, &epsdis, rhsptr);
        } else {
          CCTK_FNAME(dissipation_8_4_delta)
            (varptr, &lsh[0], &lsh[1], &lsh[2],
                              bbox, gsize, offset,
                              sbp_dx, sbp_dy, sbp_dz, &epsdis, rhsptr);
        }
      }
      break;
    }
    default:
      assert(0);
    }
  } else {
    npatches = MultiPatch_GetMaps ( cctkGH );
    patch = MultiPatch_GetMap ( cctkGH ) + 1;
    switch(order) {
    case 4: {
      if ( CCTK_Equals(dissipation_type,"Mattson-Svard-Nordstrom") ) {
      CCTK_FNAME(dissipation_4_3_opt)
        (varptr, lsh, gsh, lbnd, bbox, gsize, offset,
                      dx, &epsdis, diss_fraction, &npatches, &patch, rhsptr);
      } else {
      CCTK_FNAME(dissipation_4_3_alt)
        (varptr, lsh, gsh, lbnd, bbox, gsize, offset,
                      dx, &epsdis, diss_fraction, &npatches, &patch, rhsptr);
      }
      break;
    }
    case 6: {
      if ( CCTK_Equals(dissipation_type,"Mattson-Svard-Nordstrom") ) {
        CCTK_FNAME(dissipation_6_5_opt)
          (varptr, lsh, gsh, lbnd, bbox, gsize, offset,
                        dx, &epsdis, diss_fraction, &npatches, &patch, rhsptr);
      } else {
        CCTK_FNAME(dissipation_6_5_alt)
          (varptr, lsh, gsh, lbnd, bbox, gsize, offset,
                        dx, &epsdis, diss_fraction, &npatches, &patch, rhsptr);
      }
      break;
    }
    default:
      assert(0);
    }
  }

  if (poison_dissipation) {
      CCTK_FNAME(SBP_Poisoning)(&lsh[0], &lsh[1], &lsh[2],
                                bbox, offset, rhsptr);
  }
}
