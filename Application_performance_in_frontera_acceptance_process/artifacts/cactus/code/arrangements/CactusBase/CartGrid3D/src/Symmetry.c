/*@@
  @file      SymmetryWrappers.c
  @date      April 2000
  @author    Gerd Lanfermann
  @desc
             Routines to apply the 1/2/3D Symmetries for
             all symmetry domains (octant/bitant/quadrant).
  @enddesc
  @history
  @hdate     Sat 02 Nov 2002
  @hauthor   Thomas Radke
  @hdesc     routines generalized for applying to arbitrary CCTK data types
  @endhistory
  @version   $Id$
@@*/

#include <stdlib.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "Symmetry.h"

/* the rcs ID and its dummy function to use it */
static const char *rcsid = "$Header$";
CCTK_FILEVERSION(CactusBase_CartGrid3D_Symmetry_c)

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/
static int ApplySymmetry(const cGH *GH, int gindex, int first_vindex,
                         int numvars);

/********************************************************************
 ******************* Fortran Wrapper Prototypes *********************
 ********************************************************************/
void CCTK_FCALL CCTK_FNAME(CartSymGI)(int *ierr, const cGH **GH, int *gindex);
void CCTK_FCALL
    CCTK_FNAME(CartSymGN)(int *ierr, const cGH **GH, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(CartSymVI)(int *ierr, const cGH **GH, int *vindex);
void CCTK_FCALL
    CCTK_FNAME(CartSymVN)(int *ierr, const cGH **GH, ONE_FORTSTRING_ARG);

/********************************************************************
 ****************** External Routine Prototypes *********************
 ********************************************************************/
void DecodeSymParameters3D(int sym[6]);
void CartGrid3D_ApplyBC(CCTK_ARGUMENTS);

/*@@
   @routine    CartSymGI
   @date       April 2000
   @author     Gerd Lanfermann
   @desc
               Apply symmetry boundary routines by group index
   @enddesc
   @calls      ApplySymmetry

   @var        GH
   @vdesc      Pointer to CCTK grid hierarchy
   @vtype      const cGH *
   @vio        in
   @endvar
   @var        gindex
   @vdesc      index of group to apply symmetry BC
   @vtype      int
   @vio        in
   @endvar

   @returntype int
   @returndesc
               return code of @seeroutine ApplySymmetry <BR>
               -1 if invalid group index was given
   @endreturndesc
@@*/
int CartSymGI(const cGH *GH, int gindex) {
  int numvars, first_vindex, retval;

  numvars = CCTK_NumVarsInGroupI(gindex);
  first_vindex = CCTK_FirstVarIndexI(gindex);
  if (numvars > 0 && first_vindex >= 0) {
    char *groupname = CCTK_GroupName(gindex);
    if (!groupname)
      CCTK_VWarn(0, __LINE__, __FILE__, "CartGrid3D",
                 "error returned from function CCTK_GroupName");
    CCTK_VWarn(3, __LINE__, __FILE__, CCTK_THORNSTRING,
               "You should not call the symmetry boundary condition routines "
               "for the group \"%s\" through the CartSym* routines any more.  "
               "The symmetry boundary conditions are now applied automatically "
               "when a physical boundary condition is applied.",
               groupname);
    free(groupname);
    retval = ApplySymmetry(GH, gindex, first_vindex, numvars);
  } else {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Invalid group index %d in CartSymGI", gindex);
    retval = -1;
  }

  return (retval);
}

void CCTK_FCALL CCTK_FNAME(CartSymGI)(int *ierr, const cGH **GH, int *gindex) {
  *ierr = CartSymGI(*GH, *gindex);
}

/*@@
   @routine    CartSymGN
   @date       April 2000
   @author     Gerd Lanfermann
   @desc
               Apply symmetry boundary routines by group name
   @enddesc
   @calls      ApplySymmetry

   @var        GH
   @vdesc      Pointer to CCTK grid hierarchy
   @vtype      const cGH *
   @vio        in
   @endvar
   @var        gname
   @vdesc      name of group to apply symmetry BC
   @vtype      const char *
   @vio        in
   @endvar

   @returntype int
   @returndesc
               return code of @seeroutine ApplySymmetry <BR>
               -1 if invalid group name was given
   @endreturndesc
@@*/
int CartSymGN(const cGH *GH, const char *gname) {
  int gindex, retval;

  gindex = CCTK_GroupIndex(gname);
  if (gindex >= 0) {
    retval = CartSymGI(GH, gindex);
  } else {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Invalid group name '%s' in CartSymGN", gname);
    retval = -1;
  }

  return (retval);
}

void CCTK_FCALL
    CCTK_FNAME(CartSymGN)(int *ierr, const cGH **GH, ONE_FORTSTRING_ARG) {
  ONE_FORTSTRING_CREATE(gname)
  *ierr = CartSymGN(*GH, gname);
  free(gname);
}

/*@@
   @routine    CartSymVI
   @date       April 2000
   @author     Gerd Lanfermann
   @desc
               Apply symmetry boundary routines by variable index
   @enddesc
   @calls      ApplySymmetry

   @var        GH
   @vdesc      Pointer to CCTK grid hierarchy
   @vtype      const cGH *
   @vio        in
   @endvar
   @var        gindex
   @vdesc      index of variable to apply symmetry BC
   @vtype      int
   @vio        in
   @endvar

   @returntype int
   @returndesc
               return code of @seeroutine ApplySymmetry <BR>
               -1 if invalid variable index was given
   @endreturndesc
@@*/
int CartSymVI(const cGH *GH, int vindex) {
  int retval, gindex;

  gindex = CCTK_GroupIndexFromVarI(vindex);
  if (gindex >= 0) {
    char *fullname = CCTK_FullName(vindex);
    if (!fullname)
      CCTK_VWarn(0, __LINE__, __FILE__, "CartGrid3D",
                 "error returned from function CCTK_FullName");
    CCTK_VWarn(3, __LINE__, __FILE__, CCTK_THORNSTRING,
               "You should not call the symmetry boundary condition routines "
               "for the variable \"%s\" through the CartSym* routines any "
               "more.  The symmetry boundary conditions are now applied "
               "automatically when a physical boundary condition is applied.",
               fullname);
    free(fullname);
    retval = ApplySymmetry(GH, gindex, vindex, 1);
  } else {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Invalid variable index %d in CartSymVI", vindex);
    retval = -1;
  }

  return (retval);
}

void CCTK_FCALL CCTK_FNAME(CartSymVI)(int *ierr, const cGH **GH, int *vindex) {
  *ierr = CartSymVI(*GH, *vindex);
}

/*@@
   @routine    CartSymVN
   @date       April 2000
   @author     Gerd Lanfermann
   @desc
               Apply symmetry boundary routines by variable name
   @enddesc
   @calls      ApplySymmetry

   @var        GH
   @vdesc      Pointer to CCTK grid hierarchy
   @vtype      const cGH *
   @vio        in
   @endvar
   @var        gname
   @vdesc      name of variable to apply symmetry BC
   @vtype      const char *
   @vio        in
   @endvar

   @returntype int
   @returndesc
               return code of @seeroutine ApplySymmetry <BR>
               -1 if invalid variable name was given
   @endreturndesc
@@*/
int CartSymVN(const cGH *GH, const char *vname) {
  int vindex, retval;

  vindex = CCTK_VarIndex(vname);
  if (vindex >= 0) {
    retval = CartSymVI(GH, vindex);
  } else {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Invalid variable name '%s' in CartSymVN", vname);
    retval = -1;
  }

  return (retval);
}

void CCTK_FCALL
    CCTK_FNAME(CartSymVN)(int *ierr, const cGH **GH, ONE_FORTSTRING_ARG) {
  ONE_FORTSTRING_CREATE(vname)
  *ierr = CartSymVN(*GH, vname);
  free(vname);
}

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/
/* macro to compute the linear index of a 3D point */
#define INDEX_3D(ash, i, j, k) ((i) + (ash)[0] * ((j) + (ash)[1] * (k)))

/*@@
   @routine    SYMMETRY_BOUNDARY
   @date       Sat 02 Nov 2002
   @author     Thomas Radke
   @desc
               Macro to apply symmetry boundary conditions to a variable
               of given datatype
               Currently it is limited up to 3D variables only.
   @enddesc

   @var        cctk_type
   @vdesc      CCTK datatype of the variable
   @vtype      <cctk_type>
   @vio        in
   @endvar
   @var        itype
   @vdesc      integral CCTK datatype of the variable (used for typecasting)
   @vtype      <cctk_type>
   @vio        in
   @endvar
@@*/
#define APPLY_LOWER(dir, ii, jj, kk, jjend, kkend, iii, jjj, kkk, itype)       \
  {                                                                            \
    for (kk = 0; kk < kkend; kk++) {                                           \
      for (jj = 0; jj < jjend; jj++) {                                         \
        for (ii = 0; ii < gdata.nghostzones[dir / 2]; ii++) {                  \
          _var[INDEX_3D(ash, i, j, k)] = (itype)(                              \
              GFSym[vindex][dir] * _var[INDEX_3D(ash, iii, jjj, kkk)]);        \
        }                                                                      \
      }                                                                        \
    }                                                                          \
  }

#define APPLY_UPPER(dir, ii, jj, kk, jjend, kkend, iii, jjj, kkk, itype)       \
  {                                                                            \
    for (kk = 0; kk < kkend; kk++) {                                           \
      for (jj = 0; jj < jjend; jj++) {                                         \
        for (ii = lsh[dir / 2] - gdata.nghostzones[dir / 2];                   \
             ii < lsh[dir / 2]; ii++) {                                        \
          _var[INDEX_3D(ash, i, j, k)] = (itype)(                              \
              GFSym[vindex][dir] * _var[INDEX_3D(ash, iii, jjj, kkk)]);        \
        }                                                                      \
      }                                                                        \
    }                                                                          \
  }

#define SYMMETRY_BOUNDARY(cctk_type, itype)                                    \
  {                                                                            \
    cctk_type *_var = GH->data[vindex][0];                                     \
                                                                               \
    switch (group_static_data.dim) {                                           \
    case 3:                                                                    \
      /* apply symmetry to the z faces */                                      \
      if (doSym[4] == GFSYM_REFLECTION) {                                      \
        APPLY_LOWER(4, k, j, i, lsh[1], lsh[0], i, j, offset[4] - k, itype);   \
      }                                                                        \
      if (doSym[5] == GFSYM_REFLECTION) {                                      \
        APPLY_UPPER(5, k, j, i, lsh[1], lsh[0], i, j, offset[5] - k, itype);   \
      }                                                                        \
      if (doSym[4] == GFSYM_ROTATION_X) {                                      \
        APPLY_LOWER(4, k, j, i, lsh[1], lsh[0], i, lsh[1] - j - 1,             \
                    offset[4] - k, itype);                                     \
      }                                                                        \
      if (doSym[5] == GFSYM_ROTATION_X) {                                      \
        APPLY_UPPER(5, k, j, i, lsh[1], lsh[0], i, lsh[1] - j - 1,             \
                    offset[5] - k, itype);                                     \
      }                                                                        \
      if (doSym[4] == GFSYM_ROTATION_Y) {                                      \
        APPLY_LOWER(4, k, j, i, lsh[1], lsh[0], lsh[0] - i - 1, j,             \
                    offset[4] - k, itype);                                     \
      }                                                                        \
      if (doSym[5] == GFSYM_ROTATION_Y) {                                      \
        APPLY_UPPER(5, k, j, i, lsh[1], lsh[0], lsh[0] - i - 1, j,             \
                    offset[5] - k, itype);                                     \
      }                                                                        \
    /* FALL THROUGH */                                                         \
    case 2:                                                                    \
      /* apply symmetry to the y faces */                                      \
      if (doSym[2] == GFSYM_REFLECTION) {                                      \
        APPLY_LOWER(2, j, k, i, lsh[2], lsh[0], i, offset[2] - j, k, itype);   \
      }                                                                        \
      if (doSym[3] == GFSYM_REFLECTION) {                                      \
        APPLY_UPPER(3, j, k, i, lsh[2], lsh[0], i, offset[3] - j, k, itype);   \
      }                                                                        \
      if (doSym[2] == GFSYM_ROTATION_Z) {                                      \
        APPLY_LOWER(2, j, k, i, lsh[2], lsh[0], lsh[0] - i - 1, offset[2] - j, \
                    k, itype);                                                 \
      }                                                                        \
      if (doSym[3] == GFSYM_ROTATION_Z) {                                      \
        APPLY_UPPER(3, j, k, i, lsh[2], lsh[0], lsh[0] - i - 1, offset[3] - j, \
                    k, itype);                                                 \
      }                                                                        \
      if (group_static_data.dim > 2) {                                         \
        if (doSym[2] == GFSYM_ROTATION_X) {                                    \
          APPLY_LOWER(2, j, k, i, lsh[2], lsh[0], i, offset[2] - j,            \
                      lsh[2] - k - 1, itype);                                  \
        }                                                                      \
        if (doSym[3] == GFSYM_ROTATION_X) {                                    \
          APPLY_UPPER(3, j, k, i, lsh[2], lsh[0], i, offset[3] - j,            \
                      lsh[2] - k - 1, itype);                                  \
        }                                                                      \
      }                                                                        \
    /* FALL THROUGH */                                                         \
    case 1:                                                                    \
      /* apply symmetry to the x faces */                                      \
      if (doSym[0] == GFSYM_REFLECTION) {                                      \
        APPLY_LOWER(0, i, j, k, lsh[1], lsh[2], offset[0] - i, j, k, itype);   \
      }                                                                        \
      if (doSym[1] == GFSYM_REFLECTION) {                                      \
        APPLY_UPPER(1, i, j, k, lsh[1], lsh[2], offset[1] - i, j, k, itype);   \
      }                                                                        \
      if (group_static_data.dim > 1) {                                         \
        if (doSym[0] == GFSYM_ROTATION_Z) {                                    \
          APPLY_LOWER(0, i, j, k, lsh[1], lsh[2], offset[0] - i,               \
                      lsh[1] - j - 1, k, itype);                               \
        }                                                                      \
        if (doSym[1] == GFSYM_ROTATION_Z) {                                    \
          APPLY_UPPER(1, i, j, k, lsh[1], lsh[2], offset[1] - i,               \
                      lsh[1] - j - 1, k, itype);                               \
        }                                                                      \
      }                                                                        \
      if (group_static_data.dim > 2) {                                         \
        if (doSym[0] == GFSYM_ROTATION_Y) {                                    \
          APPLY_LOWER(0, i, j, k, lsh[1], lsh[2], offset[0] - i, j,            \
                      lsh[2] - k - 1, itype);                                  \
        }                                                                      \
        if (doSym[1] == GFSYM_ROTATION_Y) {                                    \
          APPLY_UPPER(1, i, j, k, lsh[1], lsh[2], offset[1] - i, j,            \
                      lsh[2] - k - 1, itype);                                  \
        }                                                                      \
      }                                                                        \
    /* FALL THROUGH */                                                         \
    default:                                                                   \
      ;                                                                        \
    }                                                                          \
  }

/* Function to apply above macros. */
#define SYMMETRY_FUNCTION(cctk_type, integral_type)                            \
  static void cctk_type##_SymBC(                                               \
      const cGH *GH, const int vindex, const int *doSym, const int *offset,    \
      const int *lsh, const int *ash, const cGroup group_static_data,          \
      const cGroupDynamicData gdata, int **GFSym) {                            \
    int i, j, k;                                                               \
    SYMMETRY_BOUNDARY(cctk_type, integral_type);                               \
  }

/* Create functions to apply macros.
 * This is much easier for the optiser to deal with.
 * E.g. on some machines we can't compile this file if we use the macros
 * directly in the switch statement in ApplySymmetry.
 */

SYMMETRY_FUNCTION(CCTK_COMPLEX, CCTK_REAL)
#ifdef HAVE_CCTK_COMPLEX8
SYMMETRY_FUNCTION(CCTK_COMPLEX8, CCTK_REAL4)
#endif
#ifdef HAVE_CCTK_COMPLEX16
SYMMETRY_FUNCTION(CCTK_COMPLEX16, CCTK_REAL8)
#endif
#ifdef HAVE_CCTK_COMPLEX32
SYMMETRY_FUNCTION(CCTK_COMPLEX32, CCTK_REAL16)
#endif

SYMMETRY_FUNCTION(CCTK_BYTE, CCTK_BYTE)
SYMMETRY_FUNCTION(CCTK_INT, CCTK_INT)
#ifdef HAVE_CCTK_INT1
SYMMETRY_FUNCTION(CCTK_INT1, CCTK_INT1)
#endif
#ifdef HAVE_CCTK_INT2
SYMMETRY_FUNCTION(CCTK_INT2, CCTK_INT2)
#endif
#ifdef HAVE_CCTK_INT4
SYMMETRY_FUNCTION(CCTK_INT4, CCTK_INT4)
#endif
#ifdef HAVE_CCTK_INT8
SYMMETRY_FUNCTION(CCTK_INT8, CCTK_INT8)
#endif
#ifdef HAVE_CCTK_INT16
SYMMETRY_FUNCTION(CCTK_INT16, CCTK_INT16)
#endif
SYMMETRY_FUNCTION(CCTK_REAL, CCTK_REAL)
#ifdef HAVE_CCTK_REAL4
SYMMETRY_FUNCTION(CCTK_REAL4, CCTK_REAL4)
#endif
#ifdef HAVE_CCTK_REAL8
SYMMETRY_FUNCTION(CCTK_REAL8, CCTK_REAL8)
#endif
#ifdef HAVE_CCTK_REAL16
SYMMETRY_FUNCTION(CCTK_REAL16, CCTK_REAL16)
#endif

#define CALL_SYMBC(cctk_type)                                                  \
  cctk_type##_SymBC(GH, vindex, doSym, offset, lsh, ash, group_static_data,    \
                    gdata, GFSym)

/*@@
   @routine    ApplySymmetry
   @date       Thu Mar  2 11:02:10 2000
   @author     Gerd Lanfermann
   @desc
               Apply symmetry boundary conditions to a group of grid variables
               This routine is called by the various CartSymXXX wrappers.
   @enddesc

   @var        GH
   @vdesc      Pointer to CCTK grid hierarchy
   @vtype      const cGH *
   @vio        in
   @endvar
   @var        gindex
   @vdesc      group index of the variables to apply symmetry BCs
   @vtype      int
   @vio        in
   @endvar
   @var        first_var
   @vdesc      index of first variable to apply symmetry BCs
   @vtype      int
   @vio        in
   @endvar
   @var        num_vars
   @vdesc      number of variables
   @vtype      int
   @vio        in
   @endvar

   @calls      CCTK_GroupData
               CCTK_GroupDynamicData
               SYMMETRY_BOUNDARY
   @history
   @hdate      Sat 02 Nov 2002
   @hauthor    Thomas Radke
   @hdesc      Merged separate routines for 1D, 2D, and 3D
               into a single generic routine
   @endhistory

   @returntype int
   @returndesc
                0 for success, or<BR>
               -1 if group dimension is not supported<BR>
               -2 if group datatype is not supported
   @endreturndesc
@@*/
static int ApplySymmetry(const cGH *GH, int gindex, int first_vindex,
                         int numvars) {
  int i, dim, vindex, retval;
  int **GFSym;
  int domainsym[2 * MAX_DIM], doSym[2 * MAX_DIM], offset[2 * MAX_DIM];
  int lsh[MAX_DIM], ash[MAX_DIM], cntstag[MAX_DIM];
  cGroup group_static_data;
  cGroupDynamicData gdata;
  DECLARE_CCTK_PARAMETERS

  DecodeSymParameters3D(domainsym);

  /* check if any symmetries are to be applied */
  for (i = 0; i < 2 * MAX_DIM; i++) {
    if (domainsym[i]) {
      break;
    }
  }
  if (i == 2 * MAX_DIM) {
    return (0);
  }

  /* get the group's static and dynamic data structure */
  CCTK_GroupData(gindex, &group_static_data);
  CCTK_GroupDynamicData(GH, gindex, &gdata);
  if (group_static_data.dim <= 0 || group_static_data.dim > MAX_DIM) {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "ApplySymmetry: group dimension must 1, 2, or 3");
    return (-1);
  }

  /* Avoid origin?  Default is yes */
  cntstag[0] = no_origin && no_originx && avoid_origin && avoid_originx;
  cntstag[1] = no_origin && no_originy && avoid_origin && avoid_originy;
  cntstag[2] = no_origin && no_originz && avoid_origin && avoid_originz;

  /* initialize array for variables with less dimensions than MAX_DIM
     so that we can use the INDEX_3D macro later on */
  for (i = 0; i < MAX_DIM; i++) {
    if (i < group_static_data.dim) {
      lsh[i] = gdata.lsh[i];
      ash[i] = gdata.ash[i];
    } else {
      lsh[i] = 1;
      ash[i] = 1;
    }
    offset[2 * i + 0] = 2 * gdata.nghostzones[i] - cntstag[i];
    offset[2 * i + 1] = 2 * (lsh[i] - 1) - offset[2 * i + 0];
  }

  GFSym = ((SymmetryGHex *)CCTK_GHExtension(GH, "Symmetry"))->GFSym;

  /* Apply Symmetries if:
     + the Symmetry is activated (== NOT NOSYM)
     + the Symmetry is set (== NOT UNSET)
     + the length in the direction is more than 1 grid point
     + the processor has a lower/upper physical boundary.
     Whether a grid allows a symmetry along a direction (e.g. octant=all)
     is part if the Symmetry Setup process.
  */
  retval = 0;
  for (vindex = first_vindex; vindex < first_vindex + numvars; vindex++) {
    for (dim = 0; dim < 2 * group_static_data.dim; dim++) {
      if (GFSym[vindex][dim] == GFSYM_UNSET) {
        CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Symmetries unspecified for '%s'", CCTK_VarName(vindex));
      }

      doSym[dim] = (GFSym[vindex][dim] != GFSYM_NOSYM &&
                    GFSym[vindex][dim] != GFSYM_UNSET && lsh[dim / 2] > 1 &&
                    GH->cctk_bbox[dim])
                       ? domainsym[dim]
                       : 0;
    }

    switch (group_static_data.vartype) {
    case CCTK_VARIABLE_BYTE:
      CALL_SYMBC(CCTK_BYTE);
      break;

    case CCTK_VARIABLE_INT:
      CALL_SYMBC(CCTK_INT);
      break;

    case CCTK_VARIABLE_REAL:
      CALL_SYMBC(CCTK_REAL);
      break;

    case CCTK_VARIABLE_COMPLEX:
      CALL_SYMBC(CCTK_COMPLEX);
      break;

#ifdef HAVE_CCTK_INT1
    case CCTK_VARIABLE_INT1:
      CALL_SYMBC(CCTK_INT1);
      break;
#endif

#ifdef HAVE_CCTK_INT2
    case CCTK_VARIABLE_INT2:
      CALL_SYMBC(CCTK_INT2);
      break;
#endif

#ifdef HAVE_CCTK_INT4
    case CCTK_VARIABLE_INT4:
      CALL_SYMBC(CCTK_INT4);
      break;
#endif

#ifdef HAVE_CCTK_INT8
    case CCTK_VARIABLE_INT8:
      CALL_SYMBC(CCTK_INT8);
      break;
#endif

#ifdef HAVE_CCTK_INT16
    case CCTK_VARIABLE_INT16:
      CALL_SYMBC(CCTK_INT16);
      break;
#endif

#ifdef HAVE_CCTK_REAL4
    case CCTK_VARIABLE_REAL4:
      CALL_SYMBC(CCTK_REAL4);
      break;

    case CCTK_VARIABLE_COMPLEX8:
      CALL_SYMBC(CCTK_COMPLEX8);
      break;
#endif

#ifdef HAVE_CCTK_REAL8
    case CCTK_VARIABLE_REAL8:
      CALL_SYMBC(CCTK_REAL8);
      break;

    case CCTK_VARIABLE_COMPLEX16:
      CALL_SYMBC(CCTK_COMPLEX16);
      break;
#endif

#ifdef HAVE_CCTK_REAL16
    case CCTK_VARIABLE_REAL16:
      CALL_SYMBC(CCTK_REAL16);
      break;

    case CCTK_VARIABLE_COMPLEX32:
      CALL_SYMBC(CCTK_COMPLEX32);
      break;
#endif

    default:
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Unsupported variable type %d for variable '%s'",
                 CCTK_VarTypeI(vindex), CCTK_VarName(vindex));
      retval = -2;
    }
  }

  return (retval);
}

/*@@
   @routine    CartGrid3D_ApplyBC
   @date       Sat Feb 07
   @author     Erik Schnetter
   @desc       Apply the symmetry boundary conditions
   @enddesc
@@*/

void CartGrid3D_ApplyBC(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;

  int nvars;
  CCTK_INT *restrict indices;
  CCTK_INT *restrict faces;
  CCTK_INT *restrict widths;
  CCTK_INT *restrict tables;
  int vi;
  int gi;
  int i;
  int ierr;

  if (!cctkGH)
    CCTK_VWarn(0, __LINE__, __FILE__, "CartGrid3D", "cctkGH undeclared");

  nvars = Boundary_SelectedGVs(cctkGH, 0, 0, 0, 0, 0, 0);
  if (nvars < 0)
    CCTK_VWarn(0, __LINE__, __FILE__, "CartGrid3D",
               "error returned from function Boundary_selectedGVs");

  indices = malloc(nvars * sizeof *indices);
  if (!(nvars == 0 || indices))
    CCTK_VWarn(0, __LINE__, __FILE__, "CartGrid3D",
               "error in function CartGrid3D_ApplyBC");
  faces = malloc(nvars * sizeof *faces);
  if (!(nvars == 0 || faces))
    CCTK_VWarn(0, __LINE__, __FILE__, "CartGrid3D",
               "error in function CartGrid3D_ApplyBC");
  widths = malloc(nvars * sizeof *widths);
  if (!(nvars == 0 || widths))
    CCTK_VWarn(0, __LINE__, __FILE__, "CartGrid3D",
               "error in function CartGrid3D_ApplyBC");
  tables = malloc(nvars * sizeof *tables);
  if (!(nvars == 0 || tables))
    CCTK_VWarn(0, __LINE__, __FILE__, "CartGrid3D",
               "error in function CartGrid3D_ApplyBC");

  ierr = Boundary_SelectedGVs(cctkGH, nvars, indices, faces, widths, tables, 0);
  if (!(ierr == nvars))
    CCTK_VWarn(0, __LINE__, __FILE__, "CartGrid3D",
               "error in function CartGrid3D_ApplyBC");

  for (i = 0; i < nvars; ++i) {
    vi = indices[i];
    if (!(vi >= 0 && vi < CCTK_NumVars()))
      CCTK_VWarn(0, __LINE__, __FILE__, "CartGrid3D",
                 "error in function CartGrid3D_ApplyBC");

    gi = CCTK_GroupIndexFromVarI(vi);
    if ((gi < 0))
      CCTK_VWarn(0, __LINE__, __FILE__, "CartGrid3D",
                 "error in function CartGrid3D_ApplyBC");

    ierr = ApplySymmetry(cctkGH, gi, vi, 1);
    if (ierr)
      CCTK_VWarn(0, __LINE__, __FILE__, "CartGrid3D",
                 "error in function CartGrid3D_ApplyBC");
  }

  free(indices);
  free(faces);
  free(widths);
  free(tables);
}
