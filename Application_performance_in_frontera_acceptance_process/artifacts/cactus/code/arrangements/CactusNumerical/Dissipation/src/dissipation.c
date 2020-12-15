/* $Header$ */

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void CCTK_FCALL
CCTK_FNAME(apply_dissipation) (CCTK_REAL const * const var,
                               CCTK_REAL       * const rhs,
                               int       const * const ni,
                               int       const * const nj,
                               int       const * const nk,
                               CCTK_REAL const * const dx,
                               CCTK_INT  const * const order,
                               CCTK_REAL const * const epsdis);

static void
apply (int const varindex, char const * const optstring, void * const arg);

void
dissipation_add (CCTK_ARGUMENTS)
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
  CCTK_REAL dx[3];
  int d;
  int ierr;
  
  assert (varindex >= 0);
  
  for (d=0; d<cctk_dim; ++d) {
    if (cctk_nghostzones[d] < (order+1)/2) {
      CCTK_WARN (0, "This thorn requires at least (order+1)/2 ghost zones");
    }
  }
  
  for (d=0; d<3; ++d) {
    dx[d] = CCTK_DELTA_SPACE(d);
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
  
  if (verbose) {
    char * const fullvarname = CCTK_FullName (varindex);
    char * const fullrhsname = CCTK_FullName (rhsindex);
    assert (fullvarname);
    assert (fullrhsname);
    CCTK_VInfo (CCTK_THORNSTRING,
                "Applying dissipation to \"%s\" (RHS \"%s\")",
                fullvarname, fullrhsname);
    free (fullvarname);
    free (fullrhsname);
  }
  
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

  CCTK_FNAME(apply_dissipation)
    (varptr, rhsptr, &cctk_lsh[0], &cctk_lsh[1], &cctk_lsh[2],
     dx, &order, epsdisA);
}
