/*@@
  @file    InitSymBound.F
  @date    March 1999
  @author  Gerd Lanfermann
  @desc
           Sets the symmetries for the Einstein grid functions
  @enddesc
@@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

#include "Symmetry.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusEinstein_Einstein_InitSymBound_c)

void Einstein_InitSymBound(CCTK_ARGUMENTS);

/*@@
  @routine    Einstein_InitSymBound
  @date       March 1999
  @author     Gerd Lanfermann
  @desc
              Sets the symmetries for the Einstein grid functions
  @enddesc
  @calls
  @calledby
  @history

  @endhistory

@@*/

void Einstein_InitSymBound(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS

  int one;
  int sym[3];

  one = 1;

  sym[0] = one;
  sym[1] = one;
  sym[2] = one;
  SetCartSymVN(cctkGH, sym, "admbase::gxx");
  SetCartSymVN(cctkGH, sym, "admbase::gyy");
  SetCartSymVN(cctkGH, sym, "admbase::gzz");
  sym[0] = -one;
  sym[1] = -one;
  sym[2] = one;
  SetCartSymVN(cctkGH, sym, "admbase::gxy");
  sym[0] = -one;
  sym[1] = one;
  sym[2] = -one;
  SetCartSymVN(cctkGH, sym, "admbase::gxz");
  sym[0] = one;
  sym[1] = -one;
  sym[2] = -one;
  SetCartSymVN(cctkGH, sym, "admbase::gyz");

  sym[0] = one;
  sym[1] = one;
  sym[2] = one;
  SetCartSymVN(cctkGH, sym, "admbase::kxx");
  SetCartSymVN(cctkGH, sym, "admbase::kyy");
  SetCartSymVN(cctkGH, sym, "admbase::kzz");
  sym[0] = -one;
  sym[1] = -one;
  sym[2] = one;
  SetCartSymVN(cctkGH, sym, "admbase::kxy");
  sym[0] = -one;
  sym[1] = one;
  sym[2] = -one;
  SetCartSymVN(cctkGH, sym, "admbase::kxz");
  sym[0] = one;
  sym[1] = -one;
  sym[2] = -one;
  SetCartSymVN(cctkGH, sym, "admbase::kyz");

  sym[0] = one;
  sym[1] = one;
  sym[2] = one;
  SetCartSymVN(cctkGH, sym, "admbase::alp");

  sym[0] = -one;
  sym[1] = one;
  sym[2] = one;
  SetCartSymVN(cctkGH, sym, "admbase::betax");
  sym[0] = one;
  sym[1] = -one;
  sym[2] = one;
  SetCartSymVN(cctkGH, sym, "admbase::betay");
  sym[0] = one;
  sym[1] = one;
  sym[2] = -one;
  SetCartSymVN(cctkGH, sym, "admbase::betaz");

  sym[0] = one;
  sym[1] = one;
  sym[2] = one;
  SetCartSymVN(cctkGH, sym, "admbase::dtalp");

  sym[0] = -one;
  sym[1] = one;
  sym[2] = one;
  SetCartSymVN(cctkGH, sym, "admbase::dtbetax");
  sym[0] = one;
  sym[1] = -one;
  sym[2] = one;
  SetCartSymVN(cctkGH, sym, "admbase::dtbetay");
  sym[0] = one;
  sym[1] = one;
  sym[2] = -one;
  SetCartSymVN(cctkGH, sym, "admbase::dtbetaz");

  return;
}

/* A macro for selecting boundary conditions and checking for errors */
static void select_bc(cGH const *const cctkGH, char const *const groupname,
                      int const stencil) {
  DECLARE_CCTK_PARAMETERS;

  int const ierr =
      Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, stencil, -1, groupname,
                                admbase_boundary_condition);
  if (ierr < 0) {
    CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Failed to select boundary conditions for group %s", groupname);
  }
}

/* Select boundary conditions on ADMBase variables */
void ADMBase_Boundaries(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT width[6];
  CCTK_INT ierr;
  CCTK_INT is_internal[6];
  CCTK_INT is_staggered[6];
  CCTK_INT shiftout[6];
  CCTK_INT stencil = 0, i;

  ierr =
      GetBoundarySpecification(6, width, is_internal, is_staggered, shiftout);
  if (ierr < 0) {
    CCTK_WARN(0, "Could not get the boundary specification");
  }
  for (i = 0; i < 6; i++) {
    stencil = stencil > width[i] ? stencil : width[i];
  }

  if (CCTK_EQUALS(evolution_method, "none") ||
      CCTK_EQUALS(evolution_method, "static")) {
    select_bc(cctkGH, "ADMBase::metric", stencil);
    select_bc(cctkGH, "ADMBase::curv", stencil);
  }

  if (CCTK_EQUALS(lapse_evolution_method, "static"))
    select_bc(cctkGH, "ADMBase::lapse", stencil);

  if (!CCTK_EQUALS(initial_dtlapse, "none") &&
      CCTK_EQUALS(dtlapse_evolution_method, "static"))
    select_bc(cctkGH, "ADMBase::dtlapse", stencil);

  if (!CCTK_EQUALS(initial_shift, "none") &&
      CCTK_EQUALS(shift_evolution_method, "static"))
    select_bc(cctkGH, "ADMBase::shift", stencil);

  if (!CCTK_EQUALS(initial_dtshift, "none") &&
      CCTK_EQUALS(dtshift_evolution_method, "static"))
    select_bc(cctkGH, "ADMBase::dtshift", stencil);
}
