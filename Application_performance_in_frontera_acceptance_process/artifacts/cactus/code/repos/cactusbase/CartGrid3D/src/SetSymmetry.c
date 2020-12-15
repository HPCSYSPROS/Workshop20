
/*@@
   @file      Symmetry.c
   @date      Mon Mar 15 15:09:00 1999
   @author    Gerd Lanfermann
   @desc
     This file contains the routines for registering and applying symmetry
     boundary conditions
   @enddesc
 @@*/

#include <stdlib.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_FortranString.h"
#include "Symmetry.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusBase_CartGrid3D_SetSymmetry_c)

/********************************************************************
 *********************     External Routines   **********************
 ********************************************************************/
void DecodeSymParameters3D(int sym[6]);
void CCTK_FCALL CCTK_FNAME(SetCartSymVI)(int *ierr, const cGH **GH,
                                         const int *sym, const int *vi);
void CCTK_FCALL CCTK_FNAME(SetCartSymVN)(int *ierr, const cGH **GH,
                                         const int *sym, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME(SetCartSymGI)(int *ierr, const cGH **GH,
                                         const int *sym, const int *gi);
void CCTK_FCALL CCTK_FNAME(SetCartSymGN)(int *ierr, const cGH **GH,
                                         const int *sym, ONE_FORTSTRING_ARG);

/*@@
  @routine    SetCartSymmetry
  @date       Mon Mar 15 15:10:58 1999
  @author     Gerd Lanfermann
  @desc
              This routine sets the GH extension (EinsteinBoundGHex *bGHex),
              which describes the symmetry  boundary type of each GF. Takes
              the name of the GF ("implementation::gfname") and the
              symmetry operators sx,sy,sz and inserts them in the array bGHex.
              These values will looked up by the application routines
              SymmetryWrappers
  @enddesc
  @history
              enhanced by E.Schnetter
  @endhistory
@@*/
int SetCartSymVI(const cGH *GH, const int *sym, int vi) {
  int domainsym[MAX_FACE];
  SymmetryGHex *sGHex;
  int dir;
  DECLARE_CCTK_PARAMETERS

  if (vi < 0 || vi >= CCTK_NumVars()) {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Invalid variable index %d in SetCartSymVI", vi);
    return (-1);
  }

  /* Pointer to the SymmetryGHextension */
  sGHex = (SymmetryGHex *)CCTK_GHExtension(GH, "Symmetry");

/* Reference the hash table in the GHex and tell it what kind of
   symmetry is being applied
   (depending on sym and the grid layout)
   If there is no symmetry necessary,set ESYM_NOSYM
   When we apply a symmetry and find ESYM_UNSET, something went wrong!
 */

#ifdef SYM_DEBUG
  printf("SetSymmetry: %s   [%d,%d,%d]\n", CCTK_VarName(vi), sym[0], sym[1],
         sym[2]);
#endif

  DecodeSymParameters3D(domainsym);
  for (dir = 0; dir < MAX_FACE; ++dir) {
    if (domainsym[dir] == GFSYM_REFLECTION) {
      sGHex->GFSym[vi][dir] = sym[dir / 2];
    } else if (domainsym[dir] == GFSYM_ROTATION_X) {
      sGHex->GFSym[vi][dir] = sym[1] * sym[2];
    } else if (domainsym[dir] == GFSYM_ROTATION_Y) {
      sGHex->GFSym[vi][dir] = sym[0] * sym[2];
    } else if (domainsym[dir] == GFSYM_ROTATION_Z) {
      sGHex->GFSym[vi][dir] = sym[0] * sym[1];
    } else {
      sGHex->GFSym[vi][dir] = GFSYM_NOSYM;
    }
  }

#ifdef SYM_DEBUG
  printf("SetSymmetry: %s   [%d,%d,%d]\n\n", CCTK_VarName(vi),
         sGHex->GFSym[vi][0], sGHex->GFSym[vi][2], sGHex->GFSym[vi][4]);
#endif

  return (0);
}

void CCTK_FCALL CCTK_FNAME(SetCartSymVI)(int *ierr, const cGH **GH,
                                         const int *sym, const int *vi) {
  *ierr = SetCartSymVI(*GH, sym, *vi);
}

/*@@
  @routine    SetCartSymVN
  @date       Thu May 11 13:32:55 2000
  @author     Gerd Lanfermann
  @desc
     Applies symmetry boundary conditions from
     variable index
  @enddesc
@@*/
int SetCartSymVN(const cGH *GH, const int *sym, const char *vn) {
  int vi, retval;

  vi = CCTK_VarIndex(vn);
  if (vi >= 0) {
    retval = SetCartSymVI(GH, sym, vi);
  } else {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Unknown variable '%s' in SetCartSymVN", vn);
    retval = -1;
  }

  return (retval);
}

void CCTK_FCALL CCTK_FNAME(SetCartSymVN)(int *ierr, const cGH **GH,
                                         const int *sym, ONE_FORTSTRING_ARG) {
  ONE_FORTSTRING_CREATE(vn)
  *ierr = SetCartSymVN(*GH, sym, vn);
  free(vn);
}

/*@@
  @routine SetCartSymGI
  @date
  @author  Gerd Lanfermann
  @desc
     Applies symmetry boundary conditions from
     Group index
  @enddesc
@@*/
int SetCartSymGI(const cGH *GH, const int *sym, int gi) {
  int domainsym[MAX_FACE];
  SymmetryGHex *sGHex;
  int first_vari, numvars, vi;
  int dir;
  DECLARE_CCTK_PARAMETERS

  sGHex = (SymmetryGHex *)CCTK_GHExtension(GH, "Symmetry");

  first_vari = CCTK_FirstVarIndexI(gi);
  numvars = CCTK_NumVarsInGroupI(gi);

  if (first_vari < 0) {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Cannot find group %s (grp.index: %d) in SetCartSymGI",
               CCTK_GroupName(gi), first_vari);
    return (-1);
  }

  /* Reference the hash table in the GHex and tell it what kind of
     symmetry is being applied
     (depending on sym and the grid layout)
     If there is no symmetry necessary,set ESYM_NOSYM
     When we apply a symmetry and find ESYM_UNSET, something went wrong!
   */
  for (vi = first_vari; vi < first_vari + numvars; vi++) {

#ifdef SYM_DEBUG
    printf("SetSymmetry: %s   [%d,%d,%d]\n", CCTK_VarName(vi), sym[0], sym[1],
           sym[2]);
#endif

    DecodeSymParameters3D(domainsym);
    for (dir = 0; dir < MAX_FACE; dir++) {
      if (domainsym[dir] == GFSYM_REFLECTION) {
        sGHex->GFSym[vi][dir] = sym[dir / 2];
      } else if (domainsym[dir] == GFSYM_ROTATION_X) {
        sGHex->GFSym[vi][dir] = sym[1] * sym[2];
      } else if (domainsym[dir] == GFSYM_ROTATION_Y) {
        sGHex->GFSym[vi][dir] = sym[0] * sym[2];
      } else if (domainsym[dir] == GFSYM_ROTATION_Z) {
        sGHex->GFSym[vi][dir] = sym[0] * sym[1];
      } else {
        sGHex->GFSym[vi][dir] = GFSYM_NOSYM;
      }
    }

#ifdef SYM_DEBUG
    printf("SetSymmetry: %s   [%d,%d,%d]\n\n", CCTK_VarName(vi),
           sGHex->GFSym[vi][0], sGHex->GFSym[vi][2], sGHex->GFSym[vi][4]);
#endif
  }
  return (0);
}

void CCTK_FCALL CCTK_FNAME(SetCartSymGI)(int *ierr, const cGH **GH,
                                         const int *sym, const int *gi) {
  *ierr = SetCartSymGI(*GH, sym, *gi);
}

/*@@
  @routine
  @date
  @author
  @desc
     Applies symmetry boundary conditions from
     "Implementation::Groupname"
  @enddesc
@@*/
int SetCartSymGN(const cGH *GH, const int *sym, const char *gn) {
  int gi, retval;

  gi = CCTK_GroupIndex(gn);
  if (gi >= 0) {
    retval = SetCartSymGI(GH, sym, gi);
  } else {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Unknown group '%s' in SetCartSymGN", gn);
    retval = -1;
  }

  return (retval);
}

void CCTK_FCALL CCTK_FNAME(SetCartSymGN)(int *ierr, const cGH **GH,
                                         const int *sym, ONE_FORTSTRING_ARG) {
  ONE_FORTSTRING_CREATE(gn)
  *ierr = SetCartSymGN(*GH, sym, gn);
  free(gn);
}
