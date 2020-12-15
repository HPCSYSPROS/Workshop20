/*@@
   @header    Symmetry.h
   @date      Sun 7th Mar 1999
   @author    Gerd Lanfermann
   @desc
              The extensions to the GH structure for 3D grid symmetry treatment.

              We'll have six int array for every GF, which holds a flag
              for which symmetry or (physical) bnd-condition to apply
              at the grid faces.

                * These tables are set by SetSymmetry(GF,int,int,int)
                  during initialization.
                * Default values ?
                * The information is used during evolution
                  by Einstein_DoBound(GF), Einstein_DoSym(GF).
   @enddesc
   @version   $Header$
 @@*/

#ifndef _SYMMETRY_H_
#define _SYMMETRY_H_ 1

#include "cctk.h"

#define GFSYM_UNKNOWN 0
#define GFSYM_UNSET -41
#define GFSYM_NOSYM -42

#define GFSYM_REFLECTION 1
#define GFSYM_ROTATION_X 2
#define GFSYM_ROTATION_Y 3
#define GFSYM_ROTATION_Z 4

#define MAX_DIM 3
#define MAX_FACE 6

typedef struct Symmetry {
  /* Symmetry[0..GF-1][0..dim-1] */
  /* in each direction [0,..dim-1], this will hold the symmetry
     operation across that plane, iff the grid layout requires this.
     this compares to the {sx,sy,sz} of Cactus3.2  */
  int **GFSym;
} SymmetryGHex;

#ifdef __cplusplus
extern "C" {
#endif

int GetCartSymVI(const cGH *GH, int *sym, int varindexi);
int GetCartSymVN(const cGH *GH, int *sym, const char *varname);

int SetCartSymVI(const cGH *GH, const int *sym, int varindex);
int SetCartSymGI(const cGH *GH, const int *sym, int groupindex);
int SetCartSymVN(const cGH *GH, const int *sym, const char *varname);
int SetCartSymGN(const cGH *GH, const int *sym, const char *groupname);

int CartSymVI(const cGH *GH, int varindex);
int CartSymGI(const cGH *GH, int groupindexi);
int CartSymVN(const cGH *GH, const char *varname);
int CartSymGN(const cGH *GH, const char *groupname);

#ifdef __cplusplus
}
#endif

#endif /* _SYMMETRY_H_ */
