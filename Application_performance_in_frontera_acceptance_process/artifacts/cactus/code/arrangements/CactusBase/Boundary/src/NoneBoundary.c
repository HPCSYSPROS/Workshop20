/*@@
  @file      NoneBoundary.c
  @date      Sat Jan  4 05:43:35 CET 2003
  @author    David Rideout
  @desc
             Implements 'no boundary condition'.  This can be used to
             inform the Boundary thorn that this variable should have
             boundary conditions applied to it, but then not actually
             apply any local physical boundary condition.  This is
             useful e.g. if the physical boundary is being updated in
             some other manner, but the symmetry boundaries should be
             updated normally.

             BndNone is also used as a dummy local physical boundary
             condition, when the true physical bc is non-local.
  @enddesc
  @history
  @hdate
  @hauthor
  @hdesc
  @endhistory
  @version   $Header$
@@*/

#include "cctk.h"
#include "Boundary.h"

/* the rcs ID and its dummy function to use it */
static const char *rcsid = "$Header$";
CCTK_FILEVERSION(CactusBase_Boundary_NoneBoundary_c);

/********************************************************************
 ********************    External Routines   ************************
 ********************************************************************/

/*@@
   @routine    BndNone
   @date       4 Jan 2003
   @author     David Rideout
   @desc
               Function which handles 'None' boundary condition
   @enddesc
   @calls
   @history
   @endhistory
   @var        GH
   @vdesc      Pointer to CCTK grid hierarchy
   @vtype      const cGH *
   @vio        in
   @endvar
   @var        num_vars
   @vdesc      number of variables passed in through var_indices[]
   @vtype      CCTK_INT
   @vio        in
   @endvar
   @var        var_indices
   @vdesc      array of variable indicies to which to apply this boundary
               condition
   @vtype      CCTK_INT *
   @vio        in
   @endvar
   @var        faces
   @vdesc      array of set of faces to which to apply the bc
   @vtype      CCTK_INT
   @vio        in
   @endvar
   @var        widths
   @vdesc      array of boundary widths for each variable
   @vtype      CCTK_INT
   @vio        in
   @endvar
   @var        table_handles
   @vdesc      array of table handles which hold extra arguments
   @vtype      CCTK_INT
   @vio        in
   @endvar
   @returntype CCTK_INT
   @returndesc
               0 success
   @endreturndesc
@@*/

CCTK_INT BndNone(const cGH *GH, CCTK_INT num_vars, CCTK_INT *var_indices,
                 CCTK_INT *faces, CCTK_INT *widths, CCTK_INT *table_handles) {
#ifdef DEBUG
  printf("BndNone(): got passed GH=%p, num_vars=%d, var_indices[0]=%d, "
         "table_handles[0]=%d\n",
         (const void *)GH, num_vars, var_indices[0], table_handles[0]);
#endif

  /* Ignore all input arguments */
  GH = GH;
  num_vars = num_vars;
  var_indices = var_indices;
  faces = faces;
  widths = widths;
  table_handles = table_handles;

  return 0;
}
