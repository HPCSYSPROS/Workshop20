/*@@
  @file      Check.c
  @date      19 May 2003
  @author    David Rideout
  @desc
  Check that the dimension of any grid variables is not greater than
  what the faces specification can handle (currently 15).
  @enddesc
  @history
  @hdate
  @hauthor
  @hdesc
  @endhistory
  @version   $Header$
@@*/

#include "cctk.h"
#include "cctk_Arguments.h"

#include <assert.h>
#include <stdlib.h>

#include "Boundary.h"

/* the rcs ID and its dummy function to use it */
static const char *rcsid = "$Header$";
CCTK_FILEVERSION(CactusBase_Boundary_Check_c);

/* Maximum dimension for grid variables to be selected for BC.
   Restriction is representation of set of faces in a 32 bit signed integer. */
#define MAX_DIM 15

/********************************************************************
 ********************    Scheduled Routines   ***********************
 ********************************************************************/

/*@@
   @routine    Boundary_Check
   @date       19 May 2003
   @author     David Rideout
   @desc
   Check that the dimension of any grid variables is not greater than
   what the faces specification can handle (currently 15).
   @enddesc
   @calls
   @history
   Note that in future versions of Cactus this will require the GH
   passed through CCTK_ARGUMENTS.
   @endhistory
   @var        CCTK_ARGUMENTS
   @vdesc      standard Cactus argument list
   @vtype      various
   @vio        in
   @endvar
   @returntype int
   @returndesc
               0 success
              -1 cactus executable contains variable group with too many
                dimensions
   @endreturndesc
@@*/

void Boundary_Check(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;

  if (CCTK_MaxDim() > 15) {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Executable contains variable group with more than 15 "
               "dimensions.  Thorn Boundary will not function properly for "
               "these variable groups.");
  }
}

/********************************************************************
 ********************* Externally visible helpers *******************
 ********************************************************************/

/*@@
  @routine    BndSanityCheckWidths
  @date       Wed Dec 11 17:54:25 PST 2013
  @author     Roland Haas
  @desc
              Checks that boundary width is less then 100 grid points. Aborts
              if this is not the case since this likely points to a user
              error.
  @enddesc
  @calls
  @history
  @endhistory
  @var        GH
  @vdesc      cctkGH *
  @vtype      CCTK_POINTER_TO_CONST
  @vio        in
  @endvar
  @var        varindex
  @vdesc      Cactus variable to which symmetry is applid
  @vtype      CCTK_INT
  @vio        in
  @endvar
  @var        dim
  @vdesc      Dimensionality of variable, there are 2*dim faces.
  @vtype      CCTK_INT
  @vio        int
  @endvar
  @var        boundary_widths
  @vdesc      width of boundary at each face
  @vtype      const CCTK_INT *
  @vio        in
  @endvar
  @var        bcname
  @vdesc      user readable name of Boundary condition
  @vtype      const char *
  @vio        in
  @endvar
  @returntype none
  @returndesc
              aborts on error
  @endreturndesc
@@*/
void BndSanityCheckWidths(const cGH *GH, CCTK_INT varindex, CCTK_INT dim,
                          const CCTK_INT *boundary_widths, const char *bcname) {
  for (int i = 0; i < 2 * dim; i++) {
    // due to the special meaning of width<0 in Boundary this does not catch
    // width<0 in the SelectVarForBCs call, but only negative BOUNDARY_WIDTH
    // table entries
    if (boundary_widths[i] > 100 || boundary_widths[i] < 0) {
      char *fullname = CCTK_FullName(varindex);
      const char dims[3] = "xyz";
      const char *cond = boundary_widths[i] < 0 ? "<0" : ">100";
      assert(dim < (int)sizeof(dims));
      CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                  "Tried to register a physical boundary condition \"%s\" for "
                  "variable \"%s\" with a boundary zone width %s in the %c "
                  "direction -- this looks like an error.",
                  bcname, fullname, cond, dims[dim]);
      free(fullname);
    }
  }
}
