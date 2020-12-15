/*@@
  @file      RadiationBoundary.c
  @date      Mon Mar 15 15:09:00 1999
  @author    Miguel Alcubierre, Gabrielle Allen, Gerd Lanfermann
  @desc
             Routines for applying radiation boundary conditions

             The radiative boundary condition that is implemented is:

               f  =  f0  +  u(r - v*t) / r  +  h(r + v*t) / r

             That is, I assume outgoing radial waves with a 1/r
             fall off, and the correct asymptotic value f0, plus
             I include the possibility of incoming waves as well
             (these incoming waves should be modeled somehow).

             The condition above leads to the differential equation:

               (x / r) d f  +  v d f  + v x (f - f0) / r^2  =  v x H / r^2
                 i      t         i        i                      i

             where x_i is the normal direction to the given boundaries,
             and H = 2 dh(s)/ds.

             So at a given boundary I only worry about derivatives in
             the normal direction.  Notice that u(r-v*t) has dissapeared,
             but we still do not know the value of H.

             To get H what I do is the following:  I evaluate the
             expression one point in from the boundary and solve for H
             there.  We now need a way of extrapolation H to the boundary.
             For this I assume that H falls off as a power law:

               H = k/r**n   =>  d H  =  - n H/r
                                 i

             The value of n is is defined by the parameter "radpower".
             If this parameter is negative, H is forced to be zero (this
             corresponds to pure outgoing waves and is the default).

             The behaviour I have observed is the following:  Using H=0
             is very stable, but has a very bad initial transient. Taking
             n to be 0 or positive improves the initial behaviour considerably,
             but introduces a drift that can kill the evolution at very late
             times.  Empirically, the best value I have found is n=2, for
             which the initial behaviour is very nice, and the late time drift
             is quite small.

             Another problem with this condition is that it does not
             use the physical characteristic speed, but rather it assumes
             a wave speed of v, so the boundaries should be out in
             the region where the characteristic speed is constant.
             Notice that this speed does not have to be 1.  For gauge
             quantities {alpha,phi,trK} we can have a different asymptotic
             speed, which is why the value of v is passed as a parameter.
  @enddesc
  @history
  @hdate     unknown
  @hauthor   Gerd Lanfermann
  @hdesc     Ported to Cactus 4.0
  @hdate     Fri 6 Apr 2001
  @hauthor   Thomas Radke
  @hdesc     BC routines generalized for applying to arbitrary CCTK data types
  @endhistory
  @version   $Id$
@@*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cctk.h"
#include "util_Table.h"
#include "util_ErrorCodes.h"
#include "cctk_FortranString.h"
#include "cctk_Parameters.h"

#include "Boundary.h"

/* #define DEBUG */

/* the rcs ID and its dummy function to use it */
static const char *rcsid = "$Header$";
CCTK_FILEVERSION(CactusBase_Boundary_RadiationBoundary_c);

static int ApplyBndRadiative(const cGH *GH, int stencil_dir,
                             const CCTK_INT *stencil_alldirs, int dir,
                             CCTK_REAL var0, CCTK_REAL speed,
                             CCTK_INT first_var_to, CCTK_INT first_var_from,
                             int num_vars);
static int OldApplyBndRadiative(const cGH *GH, int stencil_dir,
                                const CCTK_INT *stencil_alldirs, int dir,
                                CCTK_REAL var0, CCTK_REAL speed,
                                int first_var_to, int first_var_from,
                                int num_vars);

/********************************************************************
 ********************    External Routines   ************************
 ********************************************************************/

/*@@
   @routine    BndRadiative
   @date       6 Nov 2002
   @author     David Rideout
   @desc
               Top level function which is registered as handling
               the Radiative boundary condition
   @enddesc
   @calls      ApplyBndRadiative

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
               return code of @seeroutine ApplyBndRadiative
               -21 error reading boundary width array from table
               -22 wrong size boundary width array in table
   @endreturndesc
@@*/

CCTK_INT BndRadiative(const cGH *GH, CCTK_INT num_vars, CCTK_INT *vars,
                      CCTK_INT *faces, CCTK_INT *widths, CCTK_INT *tables) {
  int i, j, k, gi, gdim, max_gdim, err, retval;

  /* variables to pass to ApplyBndRadiative */
  CCTK_INT *width_alldirs; /* width of boundary in all directions */
  int dir;                 /* direction in which to apply bc */
  CCTK_REAL limit, speed;
  CCTK_INT
      prev_time_level; /* variable index which holds the previous time level */

#ifdef DEBUG
  printf("BndRadiative() got passed: GH=%p, num_vars=%d:\n", (const void *)GH,
         num_vars);
  printf("var index  var name  table handle\n");
  for (i = 0; i < num_vars; ++i) {
    printf("%d  %12s  %d\n", vars[i], CCTK_VarName(vars[i]), tables[i]);
  }
  printf("end of table\n");

/*  CCTK_WARN(0, "stopping code"); */
#endif

  /* Initialize variables */
  retval = 0;
  width_alldirs = NULL;
  max_gdim = 0;

  /* loop through variables, j at a time */
  for (i = 0; i < num_vars; i += j) {
    /* find other adjacent vars which are selected for identical bcs */
    j = 1;
    /* Since GFs are allowed to have different staggering, the best we
       can do is find variables of the same group which are selected
       for identical bcs.  If all GFs had the same staggering then we
       could groups many GFs together. */
    gi = CCTK_GroupIndexFromVarI(vars[i]);
#ifdef DEBUG
    printf("starting increment computation with group %d:\n", gi);
    printf("this group had %d members\n", CCTK_NumVarsInGroupI(gi));
#endif
    while (i + j < num_vars && vars[i + j] == vars[i] + j &&
           CCTK_GroupIndexFromVarI(vars[i + j]) == gi &&
           tables[i + j] == tables[i] && faces[i + j] == faces[i] &&
           widths[i + j] == widths[i]) {
      ++j;
    }

    /* Check to see if faces specification is valid */
    if (faces[i] != CCTK_ALL_FACES) {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Faces specification %d for Radiative boundary conditions on "
                 "%s is not implemented yet.  "
                 "Applying Radiative bcs to all (external) faces.",
                 (int)faces[i], CCTK_VarName(vars[i]));
    }
    dir = 0; /* apply bc to all faces */

    /* Set up default arguments for ApplyBndRadiative */
    /* Defaults for remainder of arguments */
    limit = 0.;
    prev_time_level = vars[i];
    speed = 1.;

    /* Look on table for possible non-default arguments
     * (If any of these table look-ups fail, the value will be unchanged
     * from its default value)
     */
    /* Asymptotic value of function at infinity */
    err = Util_TableGetReal(tables[i], &limit, "LIMIT");
    if (err == UTIL_ERROR_BAD_HANDLE) {
      CCTK_VWarn(5, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Invalid table handle passed for Radiative boundary "
                 "conditions for %s.  Using all default values.",
                 CCTK_VarName(vars[i]));
    } else {
      /* Wave speed */
      Util_TableGetReal(tables[i], &speed, "SPEED");
    }

    /* Determine boundary width on all faces */
    /* allocate memory for buffer */
    gdim = CCTK_GroupDimI(gi);
    if (gdim > max_gdim) {
      width_alldirs = realloc(width_alldirs, 2 * gdim * sizeof(CCTK_INT));
      max_gdim = gdim;
    }

    /* fill it with values, either from table or the boundary_width
       parameter */
    if (widths[i] < 0) {
      err = Util_TableGetIntArray(tables[i], gdim, width_alldirs,
                                  "BOUNDARY_WIDTH");
      if (err < 0) {
        CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Error %d when reading boundary width array from table "
                   "for %s",
                   err, CCTK_VarName(vars[i]));
        return -21;
      } else if (err != 2 * gdim) {
        CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Boundary width array for %s has %d elements, but %d "
                   "expected",
                   CCTK_VarName(vars[i]), err, 2 * gdim);
        return -22;
      }
    } else {
      for (k = 0; k < 2 * gdim; ++k) {
        width_alldirs[k] = widths[i];
      }
    }

    /* Apply the boundary condition */
    if ((retval = ApplyBndRadiative(GH, 0, width_alldirs, dir, limit, speed,
                                    vars[i], prev_time_level, j)) < 0) {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "ApplyBndRadiative() returned %d", retval);
    }
  }
#ifdef DEBUG
  printf("BndRadiative(): returning %d\n", retval);
#endif
  free(width_alldirs);

  return retval;
}

/* prototypes for external C routines are declared in header Boundary.h
   here only follow the fortran wrapper prototypes */
void CCTK_FCALL
    CCTK_FNAME(BndRadiativeDirGI)(int *ierr, const cGH **GH,
                                  const int *stencil_size, const int *dir,
                                  const CCTK_REAL *var0, const CCTK_REAL *speed,
                                  const int *gi_to, const int *gi_from);
void CCTK_FCALL
    CCTK_FNAME(BndRadiativeGI)(int *ierr, const cGH **GH, const int *stencil,
                               const CCTK_REAL *var0, const CCTK_REAL *speed,
                               const int *gi_to, const int *gi_from);
void CCTK_FCALL
    CCTK_FNAME(BndRadiativeDirGN)(int *ierr, const cGH **GH,
                                  const int *stencil_size, const int *dir,
                                  const CCTK_REAL *var0, const CCTK_REAL *speed,
                                  TWO_FORTSTRING_ARG);
void CCTK_FCALL
    CCTK_FNAME(BndRadiativeGN)(int *ierr, const cGH **GH, const int *stencil,
                               const CCTK_REAL *var0, const CCTK_REAL *speed,
                               TWO_FORTSTRING_ARG);
void CCTK_FCALL
    CCTK_FNAME(BndRadiativeDirVI)(int *ierr, const cGH **GH,
                                  const int *stencil_size, const int *dir,
                                  const CCTK_REAL *var0, const CCTK_REAL *speed,
                                  const int *vi_to, const int *vi_from);
void CCTK_FCALL
    CCTK_FNAME(BndRadiativeVI)(int *ierr, const cGH **GH, const int *stencil,
                               const CCTK_REAL *var0, const CCTK_REAL *speed,
                               const int *vi_to, const int *vi_from);
void CCTK_FCALL
    CCTK_FNAME(BndRadiativeDirVN)(int *ierr, const cGH **GH,
                                  const int *stencil_size, const int *dir,
                                  const CCTK_REAL *var0, const CCTK_REAL *speed,
                                  TWO_FORTSTRING_ARG);
void CCTK_FCALL
    CCTK_FNAME(BndRadiativeVN)(int *ierr, const cGH **GH, const int *stencil,
                               const CCTK_REAL *var0, const CCTK_REAL *speed,
                               TWO_FORTSTRING_ARG);

/********************************************************************
 ********************    Internal Routines   ************************
 ********************************************************************/
/*@@
   @routine    BndRadiativeDirGI
   @date       Sat Jan 20 2001
   @author     Gabrielle Allen
   @desc
               Aply radiative BC's by group index in given direction
   @enddesc
   @calls      ApplyBndRadiative

   @var        GH
   @vdesc      Pointer to CCTK grid hierarchy
   @vtype      const cGH *
   @vio        in
   @endvar
   @var        stencil_size
   @vdesc      stencil size in this direction
   @vtype      int
   @vio        in
   @endvar
   @var        dir
   @vdesc      direction to apply BC
   @vtype      int
   @vio        in
   @endvar
   @var        var0
   @vdesc      asymptotic value of function at infinity
   @vtype      CCTK_REAL
   @vio        in
   @endvar
   @var        speed
   @vdesc      wave speed
   @vtype      CCTK_REAL
   @vio        in
   @endvar
   @var        gi_to
   @vdesc      index of group to apply BC to
   @vtype      int
   @vio        in
   @endvar
   @var        gi_from
   @vdesc      index of group to apply BC from
   @vtype      int
   @vio        in
   @endvar

   @returntype int
   @returndesc
               return code of @seeroutine ApplyBndRadiative <BR>
               -1 if invalid group indices are given
   @endreturndesc
@@*/
int BndRadiativeDirGI(const cGH *GH, int stencil_size, int dir, CCTK_REAL var0,
                      CCTK_REAL speed, int gi_to, int gi_from) {
  int first_vi_to, first_vi_from, retval;

  first_vi_to = CCTK_FirstVarIndexI(gi_to);
  first_vi_from = CCTK_FirstVarIndexI(gi_from);
  if (first_vi_to >= 0 && first_vi_from >= 0) {
    retval =
        ApplyBndRadiative(GH, stencil_size, NULL, dir, var0, speed, first_vi_to,
                          first_vi_from, CCTK_NumVarsInGroupI(gi_to));
  } else {
    CCTK_VWarn(2, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Invalid group indices %d and/or %d in BndRadiativeDirGI", gi_to,
               gi_from);
    retval = -1;
  }

  return (retval);
}

void CCTK_FCALL
    CCTK_FNAME(BndRadiativeDirGI)(int *ierr, const cGH **GH,
                                  const int *stencil_size, const int *dir,
                                  const CCTK_REAL *var0, const CCTK_REAL *speed,
                                  const int *gi_to, const int *gi_from) {
  *ierr = BndRadiativeDirGI(*GH, *stencil_size, *dir, *var0, *speed, *gi_to,
                            *gi_from);
}

/*@@
   @routine    BndRadiativeGI
   @date       Tue Jul 18 18:04:07 2000
   @author     Gerd Lanfermann
   @desc
               Aply radiative BC's by group index
   @enddesc
   @calls      ApplyBndRadiative

   @var        GH
   @vdesc      Pointer to CCTK grid hierarchy
   @vtype      const cGH *
   @vio        in
   @endvar
   @var        stencil
   @vdesc      stencil width
   @vtype      int [ dimension of group ]
   @vio        in
   @endvar
   @var        var0
   @vdesc      asymptotic value of function at infinity
   @vtype      CCTK_REAL
   @vio        in
   @endvar
   @var        speed
   @vdesc      wave speed
   @vtype      CCTK_REAL
   @vio        in
   @endvar
   @var        gi_to
   @vdesc      index of group to apply BC to
   @vtype      int
   @vio        in
   @endvar
   @var        gi_from
   @vdesc      index of group to apply BC from
   @vtype      int
   @vio        in
   @endvar

   @returntype int
   @returndesc
               return code of @seeroutine ApplyBndRadiative
               -1 if invalid group indices are given
   @endreturndesc
@@*/
int BndRadiativeGI(const cGH *GH, const int *stencil, CCTK_REAL var0,
                   CCTK_REAL speed, int gi_to, int gi_from) {
  int first_vi_to, first_vi_from, retval;

  first_vi_to = CCTK_FirstVarIndexI(gi_to);
  first_vi_from = CCTK_FirstVarIndexI(gi_from);
  if (first_vi_to >= 0 && first_vi_from >= 0) {
    retval = OldApplyBndRadiative(GH, -1, (const CCTK_INT *)stencil, 0, var0,
                                  speed, first_vi_to, first_vi_from,
                                  CCTK_NumVarsInGroupI(gi_to));
  } else {
    CCTK_VWarn(2, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Invalid group indices %d and/or %d in BndRadiativeGI", gi_to,
               gi_from);
    retval = -1;
  }

  return (retval);
}

void CCTK_FCALL
    CCTK_FNAME(BndRadiativeGI)(int *ierr, const cGH **GH, const int *stencil,
                               const CCTK_REAL *var0, const CCTK_REAL *speed,
                               const int *gi_to, const int *gi_from) {
  *ierr = BndRadiativeGI(*GH, stencil, *var0, *speed, *gi_to, *gi_from);
}

/* ===================================================================== */

/*@@
   @routine    BndRadiativeDirGN
   @date       Sat Jan 20 2001
   @author     Gabrielle Allen
   @desc
               Apply radiative BC's by group name in given direction
   @enddesc
   @calls      BndRadiativeDirGI

   @var        GH
   @vdesc      Pointer to CCTK grid hierarchy
   @vtype      const cGH *
   @vio        in
   @endvar
   @var        stencil_size
   @vdesc      stencil size in this direction
   @vtype      int
   @vio        in
   @endvar
   @var        dir
   @vdesc      direction to apply BC
   @vtype      int
   @vio        in
   @endvar
   @var        var0
   @vdesc      asymptotic value of function at infinity
   @vtype      CCTK_REAL
   @vio        in
   @endvar
   @var        speed
   @vdesc      wave speed
   @vtype      CCTK_REAL
   @vio        in
   @endvar
   @var        gname_to
   @vdesc      name of group to apply BC to
   @vtype      const char *
   @vio        in
   @endvar
   @var        gname_from
   @vdesc      name of group to apply BC from
   @vtype      const char *
   @vio        in
   @endvar

   @returntype int
   @returndesc
               return code of @seeroutine BndRadiativeDirGI
               -1 if invalid group names are given
   @endreturndesc
@@*/
int BndRadiativeDirGN(const cGH *GH, int stencil_size, int dir, CCTK_REAL var0,
                      CCTK_REAL speed, const char *gname_to,
                      const char *gname_from) {
  int gi_to, gi_from, retval;

  gi_to = CCTK_GroupIndex(gname_to);
  gi_from = CCTK_GroupIndex(gname_from);
  if (gi_to >= 0 && gi_from >= 0) {
    retval =
        BndRadiativeDirGI(GH, stencil_size, dir, var0, speed, gi_to, gi_from);
  } else {
    CCTK_VWarn(2, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Invalid group names '%s' and/or '%s' in BndRadiativeDirGN",
               gname_to, gname_from);
    retval = -1;
  }

  return (retval);
}

void CCTK_FCALL
    CCTK_FNAME(BndRadiativeDirGN)(int *ierr, const cGH **GH,
                                  const int *stencil_size, const int *dir,
                                  const CCTK_REAL *var0, const CCTK_REAL *speed,
                                  TWO_FORTSTRING_ARG) {
  TWO_FORTSTRINGS_CREATE(gname_to, gname_from)
  *ierr = BndRadiativeDirGN(*GH, *stencil_size, *dir, *var0, *speed, gname_to,
                            gname_from);
  free(gname_to);
  free(gname_from);
}

/*@@
   @routine    BndRadiativeGN
   @date       Tue Jul 18 18:04:07 2000
   @author     Gerd Lanfermann
   @desc
               Aply radiative BC's by group name
   @enddesc
   @calls      BndRadiativeGI

   @var        GH
   @vdesc      Pointer to CCTK grid hierarchy
   @vtype      const cGH *
   @vio        in
   @endvar
   @var        stencil
   @vdesc      stencil width
   @vtype      int [ dimension of group ]
   @vio        in
   @endvar
   @var        var0
   @vdesc      asymptotic value of function at infinity
   @vtype      CCTK_REAL
   @vio        in
   @endvar
   @var        speed
   @vdesc      wave speed
   @vtype      CCTK_REAL
   @vio        in
   @endvar
   @var        gname_to
   @vdesc      name of group to apply BC to
   @vtype      const char *
   @vio        in
   @endvar
   @var        gname_from
   @vdesc      name of group to apply BC from
   @vtype      const char *
   @vio        in
   @endvar

   @returntype int
   @returndesc
               return code of @seeroutine BndRadiativeGI <BR>
               -1 if invalid group names are given
   @endreturndesc
@@*/
int BndRadiativeGN(const cGH *GH, const int *stencil, CCTK_REAL var0,
                   CCTK_REAL speed, const char *gname_to,
                   const char *gname_from) {
  int gi_to, gi_from, retval;

  gi_to = CCTK_GroupIndex(gname_to);
  gi_from = CCTK_GroupIndex(gname_from);
  if (gi_to >= 0 && gi_from >= 0) {
    retval = BndRadiativeGI(GH, stencil, var0, speed, gi_to, gi_from);
  } else {
    CCTK_VWarn(2, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Invalid group names '%s' and/or '%s' in BndRadiativeGN",
               gname_to, gname_from);
    retval = -1;
  }

  return (retval);
}

void CCTK_FCALL
    CCTK_FNAME(BndRadiativeGN)(int *ierr, const cGH **GH, const int *stencil,
                               const CCTK_REAL *var0, const CCTK_REAL *speed,
                               TWO_FORTSTRING_ARG) {
  TWO_FORTSTRINGS_CREATE(gname_to, gname_from)
  *ierr = BndRadiativeGN(*GH, stencil, *var0, *speed, gname_to, gname_from);
  free(gname_to);
  free(gname_from);
}

/* ===================================================================== */

/*@@
   @routine    BndRadiativeDirVI
   @date       Sat Jan 20 2001
   @author     Gabrielle Allen
   @desc
               Apply radiative BC's by variable index in given direction
   @enddesc
   @calls      ApplyBndRadiative

   @var        GH
   @vdesc      Pointer to CCTK grid hierarchy
   @vtype      const cGH *
   @vio        in
   @endvar
   @var        stencil_size
   @vdesc      stencil size in this direction
   @vtype      int
   @vio        in
   @endvar
   @var        dir
   @vdesc      direction to apply BC
   @vtype      int
   @vio        in
   @endvar
   @var        var0
   @vdesc      asymptotic value of function at infinity
   @vtype      CCTK_REAL
   @vio        in
   @endvar
   @var        speed
   @vdesc      wave speed
   @vtype      CCTK_REAL
   @vio        in
   @endvar
   @var        vi_to
   @vdesc      index of variable to apply BC to
   @vtype      int
   @vio        in
   @endvar
   @var        vi_from
   @vdesc      index of variable to apply BC from
   @vtype      int
   @vio        in
   @endvar

   @returntype int
   @returndesc
               return code of @seeroutine ApplyBndRadiative <BR>
               -1 if invalid variable indices are given
   @endreturndesc
@@*/
int BndRadiativeDirVI(const cGH *GH, int stencil_size, int dir, CCTK_REAL var0,
                      CCTK_REAL speed, int vi_to, int vi_from) {
  int retval, num_vars;

  num_vars = CCTK_NumVars();
  if (vi_to >= 0 && vi_to < num_vars && vi_from >= 0 && vi_from < num_vars) {
    retval = ApplyBndRadiative(GH, stencil_size, NULL, dir, var0, speed, vi_to,
                               vi_from, 1);
  } else {
    CCTK_VWarn(2, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Invalid variable indices %d and/or %d in BndRadiativeDirVI",
               vi_to, vi_from);
    retval = -1;
  }

  return (retval);
}

void CCTK_FCALL
    CCTK_FNAME(BndRadiativeDirVI)(int *ierr, const cGH **GH,
                                  const int *stencil_size, const int *dir,
                                  const CCTK_REAL *var0, const CCTK_REAL *speed,
                                  const int *vi_to, const int *vi_from) {
  *ierr = BndRadiativeDirVI(*GH, *stencil_size, *dir, *var0, *speed, *vi_to,
                            *vi_from);
}

/*@@
   @routine    BndRadiativeVI
   @date       Tue Jul 18 18:04:07 2000
   @author     Gerd Lanfermann
   @desc
               Apply radiative BC's by variable index
   @enddesc
   @calls      ApplyBndRadiative

   @var        GH
   @vdesc      Pointer to CCTK grid hierarchy
   @vtype      const cGH *
   @vio        in
   @endvar
   @var        stencil
   @vdesc      stencil width
   @vtype      int [ dimension of variable ]
   @vio        in
   @endvar
   @var        var0
   @vdesc      asymptotic value of function at infinity
   @vtype      CCTK_REAL
   @vio        in
   @endvar
   @var        speed
   @vdesc      wave speed
   @vtype      CCTK_REAL
   @vio        in
   @endvar
   @var        vi_to
   @vdesc      index of variable to apply BC to
   @vtype      int
   @vio        in
   @endvar
   @var        vi_from
   @vdesc      index of variable to apply BC from
   @vtype      int
   @vio        in
   @endvar

   @returntype int
   @returndesc
               return code of @seeroutine ApplyBndRadiative <BR>
               -1 if invalid variable indices are given
   @endreturndesc
@@*/
int BndRadiativeVI(const cGH *GH, const int *stencil, CCTK_REAL var0,
                   CCTK_REAL speed, int vi_to, int vi_from) {
  int retval, num_vars;

  num_vars = CCTK_NumVars();
  if (vi_to >= 0 && vi_to < num_vars && vi_from >= 0 && vi_from < num_vars) {
    retval = OldApplyBndRadiative(GH, -1, (const CCTK_INT *)stencil, 0, var0,
                                  speed, vi_to, vi_from, 1);
  } else {
    CCTK_VWarn(2, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Invalid variable indices %d and/or %d in BndRadiativeVI", vi_to,
               vi_from);
    retval = -1;
  }

  return (retval);
}

void CCTK_FCALL
    CCTK_FNAME(BndRadiativeVI)(int *ierr, const cGH **GH, const int *stencil,
                               const CCTK_REAL *var0, const CCTK_REAL *speed,
                               const int *vi_to, const int *vi_from) {
  *ierr = BndRadiativeVI(*GH, stencil, *var0, *speed, *vi_to, *vi_from);
}

/* ======================================================================= */

/*@@
   @routine    BndRadiativeDirVN
   @date       Sat Jan 20 2001
   @author     Gabrielle Allen
   @desc
               Apply radiative BC's by variable name in given direction
   @enddesc
   @calls      BndRadiativeDirVI

   @var        GH
   @vdesc      Pointer to CCTK grid hierarchy
   @vtype      const cGH *
   @vio        in
   @endvar
   @var        stencil_size
   @vdesc      stencil size in this direction
   @vtype      int
   @vio        in
   @endvar
   @var        dir
   @vdesc      direction to apply BC
   @vtype      int
   @vio        in
   @endvar
   @var        var0
   @vdesc      asymptotic value of function at infinity
   @vtype      CCTK_REAL
   @vio        in
   @endvar
   @var        speed
   @vdesc      wave speed
   @vtype      CCTK_REAL
   @vio        in
   @endvar
   @var        vname_to
   @vdesc      name of variable to apply BC to
   @vtype      const char *
   @vio        in
   @endvar
   @var        vname_from
   @vdesc      name of variable to apply BC from
   @vtype      const char *
   @vio        in
   @endvar

   @returntype int
   @returndesc
               return code of @seeroutine BndRadiativeDirVI <BR>
               -1 if invalid variable names are given
   @endreturndesc
@@*/
int BndRadiativeDirVN(const cGH *GH, int stencil_size, int dir, CCTK_REAL var0,
                      CCTK_REAL speed, const char *vname_to,
                      const char *vname_from) {
  int vi_to, vi_from, num_vars, retval;

  vi_to = CCTK_VarIndex(vname_to);
  vi_from = CCTK_VarIndex(vname_from);
  num_vars = CCTK_NumVars();

  if (vi_to >= 0 && vi_to < num_vars && vi_from >= 0 && vi_from < num_vars) {
    retval =
        BndRadiativeDirVI(GH, stencil_size, dir, var0, speed, vi_to, vi_from);
  } else {
    CCTK_VWarn(2, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Invalid variable names '%s' and/or '%s' in BndRadiativeDirVN",
               vname_to, vname_from);
    retval = -1;
  }

  return (retval);
}

void CCTK_FCALL
    CCTK_FNAME(BndRadiativeDirVN)(int *ierr, const cGH **GH,
                                  const int *stencil_size, const int *dir,
                                  const CCTK_REAL *var0, const CCTK_REAL *speed,
                                  TWO_FORTSTRING_ARG) {
  TWO_FORTSTRINGS_CREATE(vname_to, vname_from)
  *ierr = BndRadiativeDirVN(*GH, *stencil_size, *dir, *var0, *speed, vname_to,
                            vname_from);
  free(vname_to);
  free(vname_from);
}

/*@@
   @routine    BndRadiativeVN
   @date       Tue Jul 18 18:04:07 2000
   @author     Gerd Lanfermann
   @desc
               Apply radiative BC's by variable name
   @enddesc
   @calls      BndRadiativeVI

   @var        GH
   @vdesc      Pointer to CCTK grid hierarchy
   @vtype      const cGH *
   @vio        in
   @endvar
   @var        stencil
   @vdesc      stencil width
   @vtype      int [ dimension of variable ]
   @vio        in
   @endvar
   @var        var0
   @vdesc      asymptotic value of function at infinity
   @vtype      CCTK_REAL
   @vio        in
   @endvar
   @var        speed
   @vdesc      wave speed
   @vtype      CCTK_REAL
   @vio        in
   @endvar
   @var        vname_to
   @vdesc      name of variable to apply BC to
   @vtype      const char *
   @vio        in
   @endvar
   @var        vname_from
   @vdesc      name of variable to apply BC from
   @vtype      const char *
   @vio        in
   @endvar

   @returntype int
   @returndesc
               return code of @seeroutine BndRadiativeVI <BR>
               -1 if invalid variable names are given
   @endreturndesc
@@*/
int BndRadiativeVN(const cGH *GH, const int *stencil, CCTK_REAL var0,
                   CCTK_REAL speed, const char *vname_to,
                   const char *vname_from) {
  int vi_to, vi_from, num_vars, retval;

  vi_to = CCTK_VarIndex(vname_to);
  vi_from = CCTK_VarIndex(vname_from);
  num_vars = CCTK_NumVars();

  if (vi_to >= 0 && vi_to < num_vars && vi_from >= 0 && vi_from < num_vars) {
    retval = BndRadiativeVI(GH, stencil, var0, speed, vi_to, vi_from);
  } else {
    CCTK_VWarn(2, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Invalid variable names '%s' and/or '%s' in BndRadiativeVN",
               vname_to, vname_from);
    retval = -1;
  }

  return (retval);
}

void CCTK_FCALL
    CCTK_FNAME(BndRadiativeVN)(int *ierr, const cGH **GH, const int *stencil,
                               const CCTK_REAL *var0, const CCTK_REAL *speed,
                               TWO_FORTSTRING_ARG) {
  TWO_FORTSTRINGS_CREATE(vname_to, vname_from)
  *ierr = BndRadiativeVN(*GH, stencil, *var0, *speed, vname_to, vname_from);
  free(vname_to);
  free(vname_from);
}

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/

/* shortcut for multiplying a variable with itself */
#define SQR(a) ((a) * (a))

/* the maximum dimension we can deal with */
#define MAXDIM 3

/*@@
   @routine    LOWER_RADIATIVE_BOUNDARY_3D
   @date       Mon 9 Apr 2001
   @author     Thomas Radke
   @desc
               Macro to apply radiative BC to a lower bound of a 3D variable
   @enddesc

   @var        istart, jstart, kstart
   @vdesc      start index for the x,y,z direction
   @vtype      int
   @vio        in
   @endvar
   @var        dim
   @vdesc      dimension to apply BC
   @vtype      int
   @vio        in
   @endvar
   @var        cctk_type
   @vdesc      CCTK datatypes of the source and target variable
   @vtype      <cctk_type>
   @vio        in
   @endvar
@@*/
#define LOWER_RADIATIVE_BOUNDARY_3D(istart, jstart, kstart, dim, cctk_type)    \
  {                                                                            \
    int _i, _j, _k;                                                            \
    int _0 = 0 * offset[dim], _1 = 1 * offset[dim], _2 = 2 * offset[dim];      \
                                                                               \
    for (_k = kstart - 1; _k >= 0; _k--) {                                     \
      for (_j = jstart - 1; _j >= 0; _j--) {                                   \
        int _idx = CCTK_GFINDEX3D(GH, istart - 1, _j, _k);                     \
        const CCTK_REAL *_r = xyzr[MAXDIM] + _idx, *_xyz = xyzr[dim] + _idx;   \
        cctk_type *_to = (cctk_type *)to_ptr + _idx;                           \
        const cctk_type *_from = (const cctk_type *)from_ptr + _idx;           \
                                                                               \
        for (_i = istart - 1; _i >= 0; _i--) {                                 \
          CCTK_REAL _r0_inv = 1 / _r[_0], _r1_inv = 1 / _r[_1];                \
                                                                               \
          if (radpower > 0) {                                                  \
            CCTK_REAL H;                                                       \
                                                                               \
            H = 0.25 * radpower * dxyz[dim] *                                  \
                (_xyz[_0] * SQR(_r0_inv) + _xyz[_1] * SQR(_r1_inv));           \
            H = (1 + H) / (1 - H);                                             \
            H *= dtv * (0.25 * (_to[_1] + _to[_2] + _from[_1] + _from[_2]) -   \
                        var0) +                                                \
                 0.5 * (_r[_1] * (_to[_1] - _from[_1]) +                       \
                        _r[_2] * (_to[_2] - _from[_2])) +                      \
                 0.25 * (_to[_2] - _to[_1] + _from[_2] - _from[_1]) *          \
                     rho[dim] *                                                \
                     (SQR(_r[_1]) / _xyz[_1] + SQR(_r[_2]) / _xyz[_2]);        \
            dtvvar0H = dtvvar0 + H;                                            \
          }                                                                    \
                                                                               \
          _to[_0] = (cctk_type)(                                               \
              (dtvvar0H *                                                      \
                   (_xyz[_0] * SQR(_r0_inv) + _xyz[_1] * SQR(_r1_inv)) -       \
               _to[_1] *                                                       \
                   (rho[dim] + _xyz[_1] * _r1_inv * (1 + dtvh * _r1_inv)) +    \
               _from[_0] *                                                     \
                   (rho[dim] + _xyz[_0] * _r0_inv * (1 - dtvh * _r0_inv)) -    \
               _from[_1] *                                                     \
                   (rho[dim] - _xyz[_1] * _r1_inv * (1 - dtvh * _r1_inv))) /   \
              (-rho[dim] + _xyz[_0] * _r0_inv * (1 + dtvh * _r0_inv)));        \
          _r--;                                                                \
          _xyz--;                                                              \
          _to--;                                                               \
          _from--;                                                             \
        }                                                                      \
      }                                                                        \
    }                                                                          \
  }

/*@@
   @routine    UPPER_RADIATIVE_BOUNDARY_3D
   @date       Mon 9 Apr 2001
   @author     Thomas Radke
   @desc
               Macro to apply radiative BC to an upper bound of a 3D variable
   @enddesc

   @var        istart, jstart, kstart
   @vdesc      start index for the x,y,z direction
   @vtype      int
   @vio        in
   @endvar
   @var        dim
   @vdesc      dimension to apply BC
   @vtype      int
   @vio        in
   @endvar
   @var        cctk_type
   @vdesc      CCTK datatypes of the source and target variable
   @vtype      <cctk_type>
   @vio        in
   @endvar
@@*/
#define UPPER_RADIATIVE_BOUNDARY_3D(istart, jstart, kstart, dim, cctk_type)    \
  {                                                                            \
    int _i, _j, _k;                                                            \
    int _0 = -0 * offset[dim], _1 = -1 * offset[dim], _2 = -2 * offset[dim];   \
                                                                               \
    for (_k = kstart; _k < GH->cctk_lsh[2]; _k++) {                            \
      for (_j = jstart; _j < GH->cctk_lsh[1]; _j++) {                          \
        int _idx = CCTK_GFINDEX3D(GH, istart, _j, _k);                         \
        const CCTK_REAL *_r = xyzr[MAXDIM] + _idx, *_xyz = xyzr[dim] + _idx;   \
        cctk_type *_to = (cctk_type *)to_ptr + _idx;                           \
        const cctk_type *_from = (const cctk_type *)from_ptr + _idx;           \
                                                                               \
        for (_i = istart; _i < GH->cctk_lsh[0]; _i++) {                        \
          CCTK_REAL _r0_inv = 1 / _r[_0], _r1_inv = 1 / _r[_1];                \
                                                                               \
          if (radpower > 0) {                                                  \
            CCTK_REAL H;                                                       \
                                                                               \
            H = 0.25 * radpower * dxyz[dim] *                                  \
                (_xyz[_0] * SQR(_r0_inv) + _xyz[_1] * SQR(_r1_inv));           \
            H = (1 - H) / (1 + H);                                             \
            H *= dtv * (0.25 * (_to[_1] + _to[_2] + _from[_1] + _from[_2]) -   \
                        var0) +                                                \
                 0.5 * (_r[_1] * (_to[_1] - _from[_1]) +                       \
                        _r[_2] * (_to[_2] - _from[_2])) +                      \
                 0.25 * (_to[_1] - _to[_2] + _from[_1] - _from[_2]) *          \
                     rho[dim] *                                                \
                     (SQR(_r[_1]) / _xyz[_1] + SQR(_r[_2]) / _xyz[_2]);        \
            dtvvar0H = dtvvar0 + H;                                            \
          }                                                                    \
                                                                               \
          _to[_0] = (cctk_type)(                                               \
              (dtvvar0H *                                                      \
                   (_xyz[_0] * (SQR(_r0_inv)) + _xyz[_1] * (SQR(_r1_inv))) +   \
               _to[_1] *                                                       \
                   (rho[dim] - _xyz[_1] * _r1_inv * (1 + dtvh * _r1_inv)) +    \
               _from[_0] *                                                     \
                   (-rho[dim] + _xyz[_0] * _r0_inv * (1 - dtvh * _r0_inv)) +   \
               _from[_1] *                                                     \
                   (rho[dim] + _xyz[_1] * _r1_inv * (1 - dtvh * _r1_inv))) /   \
              (rho[dim] + _xyz[_0] * _r0_inv * (1 + dtvh * _r0_inv)));         \
          _r++;                                                                \
          _xyz++;                                                              \
          _to++;                                                               \
          _from++;                                                             \
        }                                                                      \
      }                                                                        \
    }                                                                          \
  }

/*@@
   @routine    RADIATIVE_BOUNDARY
   @date       Mon 9 Apr 2001
   @author     Thomas Radke
   @desc
               Macro to apply radiative BC to a variable
               Currently it is limited to 3D variables only.
   @enddesc
   @calls      LOWER_RADIATIVE_BOUNDARY_3D
               UPPER_RADIATIVE_BOUNDARY_3D

   @var        lsh
   @vdesc      local shape of the variable
   @vtype      int [ dim ]
   @vio        in
   @endvar
   @var        stencil
   @vdesc      stencils in every direction
   @vtype      int [ 2*dim ]
   @vio        in
   @endvar
   @var        cctk_type
   @vdesc      CCTK datatypes of the source and target variable
   @vtype      <cctk_type>
   @vio        in
   @endvar
@@*/
#define RADIATIVE_BOUNDARY(lsh, stencil, cctk_type)                            \
  {                                                                            \
    /* check the dimensionality */                                             \
    if (gdim != 3) {                                                           \
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,                      \
                 "ApplyBndRadiative: variable dimension of %d not supported",  \
                 gdim);                                                        \
      return (-5);                                                             \
    }                                                                          \
                                                                               \
    /* Lower x-bound */                                                        \
    if (doBC[0]) {                                                             \
      LOWER_RADIATIVE_BOUNDARY_3D(stencil[0], lsh[1], lsh[2], 0, cctk_type);   \
    }                                                                          \
                                                                               \
    /* Upper x-bound */                                                        \
    if (doBC[1]) {                                                             \
      UPPER_RADIATIVE_BOUNDARY_3D(lsh[0] - stencil[1], 0, 0, 0, cctk_type);    \
    }                                                                          \
                                                                               \
    /* Lower y-bound */                                                        \
    if (doBC[2]) {                                                             \
      LOWER_RADIATIVE_BOUNDARY_3D(lsh[0], stencil[2], lsh[2], 1, cctk_type);   \
    }                                                                          \
                                                                               \
    /* Upper y-bound */                                                        \
    if (doBC[3]) {                                                             \
      UPPER_RADIATIVE_BOUNDARY_3D(0, lsh[1] - stencil[3], 0, 1, cctk_type);    \
    }                                                                          \
                                                                               \
    /* Lower z-bound */                                                        \
    if (doBC[4]) {                                                             \
      LOWER_RADIATIVE_BOUNDARY_3D(lsh[0], lsh[1], stencil[4], 2, cctk_type);   \
    }                                                                          \
                                                                               \
    /* Upper z-bound */                                                        \
    if (doBC[5]) {                                                             \
      UPPER_RADIATIVE_BOUNDARY_3D(0, 0, lsh[2] - stencil[5], 2, cctk_type);    \
    }                                                                          \
  }

/*@@
   @routine    ApplyBndRadiative
   @date       Tue Jul 18 18:04:07 2000
   @author     Gerd Lanfermann
   @desc
               Apply radiation boundary conditions to a group of grid functions
               given by their indices
               This routine is called by the various BndRadiativeXXX wrappers.

               Although it is currently limited to handle 3D variables only
               it can easily be extended for other dimensions
               by adapting the appropriate macros.
   @enddesc

   @var        GH
   @vdesc      Pointer to CCTK grid hierarchy
   @vtype      const cGH *
   @vio        in
   @endvar
   @var        width_dir
   @vdesc      boundary width in direction dir
   @vtype      int
   @vio        in
   @endvar
   @var        in_widths
   @vdesc      boundary widths for all directions
   @vtype      const CCTK_INT [ dimension of variable(s) ]
   @vio        in
   @endvar
   @var        dir
   @vdesc      direction to copy boundaries (0 for copying all directions)
   @vtype      int
   @vio        in
   @endvar
   @var        var0
   @vdesc      asymptotic value of function at infinity
   @vtype      CCTK_REAL
   @vio        in
   @endvar
   @var        speed
   @vdesc      wave speed
   @vtype      CCTK_REAL
   @vio        in
   @endvar
   @var        first_var_to
   @vdesc      index of first variable to copy boundaries to
   @vtype      CCTK_INT
   @vio        in
   @endvar
   @var        first_var_from
   @vdesc      index of first variable to copy boundaries from
   @vtype      CCTK_INT
   @vio        in
   @endvar
   @var        num_vars
   @vdesc      number of variables
   @vtype      int
   @vio        in
   @endvar
   @calls      CCTK_VarTypeI
               CCTK_GroupDimFromVarI
               RADIATIVE_BOUNDARY
   @history
   @hdate      Mon 9 Apr 2001
   @hauthor    Thomas Radke
   @hdesc      Merged separate routines for 1D, 2D, and 3D
               into a single generic routine
   @endhistory

   @returntype int
   @returndesc
                0 for success
               -1 if abs(direction) is greater than variables' dimension
               -2 if variable dimension is not supported
               -3 if NULL pointer passed as boundary width array
               -4 if variable type is not supported
               -5 if variable dimension is other than 3D
               -6 if a coordinate is not found
   @endreturndesc
@@*/
static int ApplyBndRadiative(const cGH *GH, int width_dir,
                             const CCTK_INT *in_widths, int dir, CCTK_REAL var0,
                             CCTK_REAL speed, CCTK_INT first_var_to,
                             CCTK_INT first_var_from, int num_vars) {
  int i, gdim, indx;
  int var_to, var_from;
  int timelvl_from;
  char coord_system_name[10];
  CCTK_REAL dxyz[MAXDIM], rho[MAXDIM];
  const CCTK_REAL *xyzr[MAXDIM + 1];
  CCTK_INT doBC[2 * MAXDIM], widths[2 * MAXDIM], offset[MAXDIM];
  CCTK_INT symtable;
  CCTK_INT symbnd[2 * MAXDIM];
  CCTK_INT is_physical[2 * MAXDIM];
  CCTK_INT ierr;
  CCTK_REAL dtv, dtvh, dtvvar0, dtvvar0H;
  void *to_ptr;
  const void *from_ptr;
  DECLARE_CCTK_PARAMETERS

  /* check the direction parameter */
  if (abs(dir) > MAXDIM) {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "ApplyBndRadiative: direction %d is greater than maximum "
               "dimension %d",
               dir, MAXDIM);
    return (-1);
  }

  /* get the dimensionality */
  gdim = CCTK_GroupDimFromVarI(first_var_to);

  /* check the dimensionality */
  if (gdim > MAXDIM) {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "ApplyBndRadiative: variable dimension of %d not supported",
               gdim);
    return (-2);
  }

  /* set up boundary width array */
  if (dir) {
    widths[2 * (abs(dir) - 1)] = width_dir;
    widths[2 * (abs(dir) - 1) + 1] = width_dir;
  } else if (in_widths) {
    memcpy(widths, in_widths, 2 * gdim * sizeof *widths);
  } else {
    CCTK_WARN(1, "ApplyBndRadiative: NULL pointer passed "
                 "for boundary width array");
    return (-3);
  }

  /* sanity check on width of boundary,  */
  BndSanityCheckWidths(GH, first_var_to, gdim, widths, "Radiative");

  /* Use next time level, if available */
  timelvl_from = 0;
  if (CCTK_DeclaredTimeLevelsVI(first_var_from) > 1) {
    timelvl_from = 1;
  }

  /* Find Courant parameters. */
  dtv = speed * GH->cctk_delta_time;
  dtvh = 0.5 * dtv;
  dtvvar0 = dtv * var0;
  dtvvar0H = dtvvar0;

  sprintf(coord_system_name, "cart%dd", gdim);
  for (i = 0; i < gdim; i++) {
    /* Radiative boundaries need the underlying Cartesian coordinates */
    indx = CCTK_CoordIndex(i + 1, NULL, coord_system_name);
    if (indx < 0) {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Coordinate for system %s not found", coord_system_name);
      return (-6);
    }
    xyzr[i] = GH->data[indx][0];

    /* According to the Cactus spec, the true delta_space values for a
       grid are calculated as follows: */
    dxyz[i] = GH->cctk_delta_space[i] / GH->cctk_levfac[i];

    rho[i] = dtv / dxyz[i];

    offset[i] = i == 0 ? 1 : offset[i - 1] * GH->cctk_lsh[i - 1];
  }

  /* Append r grid variable to end of xyzr[] array */
  sprintf(coord_system_name, "spher%dd", gdim);
  indx = CCTK_CoordIndex(-1, "r", coord_system_name);
  if (indx < 0) {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Coordinate for system %s not found", coord_system_name);
    return (-6);
  }
  xyzr[MAXDIM] = GH->data[indx][0];

  /* see if we have a physical boundary */
  symtable = SymmetryTableHandleForGrid(GH);
  if (symtable < 0)
    CCTK_WARN(0, "internal error");
  ierr = Util_TableGetIntArray(symtable, 2 * gdim, symbnd, "symmetry_handle");
  if (ierr != 2 * gdim)
    CCTK_WARN(0, "internal error");
  for (i = 0; i < 2 * gdim; i++) {
    is_physical[i] = symbnd[i] < 0;
  }

  /* now loop over all variables */
  for (var_to = first_var_to, var_from = first_var_from;
       var_to < first_var_to + num_vars; var_to++, var_from++) {
    to_ptr = GH->data[var_to][0];
    from_ptr = GH->data[var_from][timelvl_from];

    /* Apply condition if:
       + boundary is a physical boundary
       + boundary is an outer boundary
       + have enough grid points
    */
    for (i = 0; i < 2 * MAXDIM; i++) {
      doBC[i] = is_physical[i];
    }
    for (i = 0; i < MAXDIM; i++) {
      doBC[i * 2] &= GH->cctk_lsh[i] > widths[i * 2] && GH->cctk_bbox[i * 2];
      doBC[i * 2 + 1] &=
          GH->cctk_lsh[i] > widths[i * 2 + 1] && GH->cctk_bbox[i * 2 + 1];
      if (dir != 0) {
        doBC[i * 2] &= (dir < 0 && (i + 1 == abs(dir)));
        doBC[i * 2 + 1] &= (dir > 0 && (i + 1 == abs(dir)));
      }
    }

    switch (CCTK_VarTypeI(var_to)) {
    case CCTK_VARIABLE_REAL:
      RADIATIVE_BOUNDARY(GH->cctk_lsh, widths, CCTK_REAL);
      break;

#ifdef HAVE_CCTK_REAL4
    case CCTK_VARIABLE_REAL4:
      RADIATIVE_BOUNDARY(GH->cctk_lsh, widths, CCTK_REAL4);
      break;
#endif

#ifdef HAVE_CCTK_REAL8
    case CCTK_VARIABLE_REAL8:
      RADIATIVE_BOUNDARY(GH->cctk_lsh, widths, CCTK_REAL8);
      break;
#endif

#ifdef HAVE_CCTK_REAL16
    case CCTK_VARIABLE_REAL16:
      RADIATIVE_BOUNDARY(GH->cctk_lsh, widths, CCTK_REAL16);
      break;
#endif

    default:
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Unsupported variable type %d for variable '%s'",
                 CCTK_VarTypeI(var_to), CCTK_VarName(var_to));
      return (-4);
    }
  }

  return (0);
}

/*@@
   @routine    OldApplyBndRadiative
   @date       5 May 2003
   @author     David Rideout
   @desc
               The new boundary API expects a 2d-element array for the
               boundary_widths (d=dimension of grid variable), while
               the old API expects a d-element array.  This function
               converts the old array to the new format.
   @enddesc

   @var        GH
   @vdesc      Pointer to CCTK grid hierarchy
   @vtype      const cGH *
   @vio        in
   @endvar
   @var        width_dir
   @vdesc      boundary width in direction dir
   @vtype      int
   @vio        in
   @endvar
   @var        stencil_alldirs
   @vdesc      boundary widths for all directions
   @vtype      const CCTK_INT [ dimension of variable(s) ]
   @vio        in
   @endvar
   @var        dir
   @vdesc      direction to copy boundaries (0 for copying all directions)
   @vtype      int
   @vio        in
   @endvar
   @var        var0
   @vdesc      asymptotic value of function at infinity
   @vtype      CCTK_REAL
   @vio        in
   @endvar
   @var        speed
   @vdesc      wave speed
   @vtype      CCTK_REAL
   @vio        in
   @endvar
   @var        first_var_to
   @vdesc      index of first variable to copy boundaries to
   @vtype      int
   @vio        in
   @endvar
   @var        first_var_from
   @vdesc      index of first variable to copy boundaries from
   @vtype      int
   @vio        in
   @endvar
   @var        num_vars
   @vdesc      number of variables
   @vtype      int
   @vio        in
   @endvar

   @calls      CCTK_GroupIndexFromVarI
               ApplyBndRadiative
   @returntype int
   @returndesc
               returncode from @seeroutine ApplyBndRadiative
   @endreturndesc
@@*/

static int OldApplyBndRadiative(const cGH *GH, int width_dir,
                                const CCTK_INT *stencil_alldirs, int dir,
                                CCTK_REAL var0, CCTK_REAL speed,
                                int first_var_to, int first_var_from,
                                int num_vars) {
  int i, dim, retval;
  CCTK_INT *boundary_widths;
  static int warned;

  /* Convert stencil_alldirs to new format */
  dim = CCTK_GroupDimFromVarI(first_var_to);
  boundary_widths = malloc(2 * dim * sizeof(CCTK_INT));
  for (i = 0; i < 2 * dim; ++i) {
    boundary_widths[i] = stencil_alldirs[i / 2];
  }

  /* Bug people for using the old interface */
  if (!warned) {
    CCTK_VWarn(2, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Copied older d-element array of boundary widths into the "
               "newer 2d-element format.  Please use the new boundary "
               "interface to avoid this.");
    warned = 1;
  }

  /* Call ApplyBnd... with new boundary width array */
  retval = ApplyBndRadiative(GH, width_dir, boundary_widths, dir, var0, speed,
                             first_var_to, first_var_from, num_vars);

  free(boundary_widths);
  return retval;
}
