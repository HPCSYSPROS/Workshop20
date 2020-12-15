/*@@
  @file      Boundary.h
  @date      Tue Sep 26 11:50:46 2000
  @author    Gerd Lanfermann
  @desc
             Prototypes for boundary routines
  @enddesc
  @version   $Header$
@@*/

#ifndef _BOUNDARY_H_
#define _BOUNDARY_H_

#ifdef __cplusplus
extern "C" {
#endif

/* data type for pointer to function which implements a physical boundary
   condition: */
typedef CCTK_INT (*phys_bc_fn_ptr)(CCTK_POINTER_TO_CONST, const CCTK_INT,
                                   const CCTK_INT *, const CCTK_INT *,
                                   const CCTK_INT *, const CCTK_INT *);

/* check boundary withd and abort if unlikely large (>100 points) */
void BndSanityCheckWidths(const cGH *GH, CCTK_INT varindex, CCTK_INT dim,
                          const CCTK_INT *boundary_widths, const char *bcname);

/* prototype for routine registed as providing 'None' boundary condition */
CCTK_INT BndNone(const cGH *GH, CCTK_INT num_vars, CCTK_INT *var_indicies,
                 CCTK_INT *faces, CCTK_INT *boundary_widths,
                 CCTK_INT *table_handle);

/* Scalar boundaries */

/* prototype for routine registed as providing 'Scalar' boundary condition */
CCTK_INT BndScalar(const cGH *GH, CCTK_INT num_vars, CCTK_INT *var_indicies,
                   CCTK_INT *faces, CCTK_INT *boundary_widths,
                   CCTK_INT *table_handle);

int BndScalarDirGI(const cGH *GH, int stencil_size, int dir, CCTK_REAL var0,
                   int gi);
int BndScalarDirGN(const cGH *GH, int stencil_size, int dir, CCTK_REAL var0,
                   const char *gname);
int BndScalarDirVI(const cGH *GH, int stencil_size, int dir, CCTK_REAL var0,
                   int vi);
int BndScalarDirVN(const cGH *GH, int stencil_size, int dir, CCTK_REAL var0,
                   const char *vname);

int BndScalarGI(const cGH *GH, const int *stencil, CCTK_REAL var0, int gi);
int BndScalarGN(const cGH *GH, const int *stencil, CCTK_REAL var0,
                const char *gname);
int BndScalarVI(const cGH *GH, const int *stencil, CCTK_REAL var0, int vi);
int BndScalarVN(const cGH *GH, const int *stencil, CCTK_REAL var0,
                const char *vname);

/* Copying boundaries */

/* prototype for routine registed as providing 'Copy' boundary condition */
CCTK_INT BndCopy(const cGH *GH, CCTK_INT num_vars, CCTK_INT *var_indicies,
                 CCTK_INT *faces, CCTK_INT *boundary_widths,
                 CCTK_INT *table_handle);

int BndCopyDirGI(const cGH *GH, int stencil_size, int dir, int gi_to,
                 int gi_from);
int BndCopyDirGN(const cGH *GH, int stencil_size, int dir, const char *gname_to,
                 const char *gname_from);
int BndCopyDirVI(const cGH *GH, int stencil_size, int dir, int vi_to,
                 int vi_from);
int BndCopyDirVN(const cGH *GH, int stencil_size, int dir, const char *vname_to,
                 const char *vname_from);

int BndCopyGI(const cGH *GH, const int *stencil, int gi_to, int gi_from);
int BndCopyGN(const cGH *GH, const int *stencil, const char *gname_to,
              const char *gname_from);
int BndCopyVI(const cGH *GH, const int *stencil, int vi_to, int vi_from);
int BndCopyVN(const cGH *GH, const int *stencil, const char *vname_to,
              const char *vname_from);

/* Static boundaries  */

/* prototype for routine registed as providing 'Static' boundary condition */
CCTK_INT BndStatic(const cGH *GH, CCTK_INT num_vars, CCTK_INT *var_indicies,
                   CCTK_INT *faces, CCTK_INT *boundary_widths,
                   CCTK_INT *table_handle);

int BndStaticDirGI(const cGH *GH, int stencil_size, int dir, int gi);
int BndStaticDirGN(const cGH *GH, int stencil_size, int dir, const char *gname);
int BndStaticDirVI(const cGH *GH, int stencil_size, int dir, int vi);
int BndStaticDirVN(const cGH *GH, int stencil_size, int dir, const char *vname);

int BndStaticGI(const cGH *GH, const int *stencil, int gi);
int BndStaticGN(const cGH *GH, const int *stencil, const char *gname);
int BndStaticVI(const cGH *GH, const int *stencil, int vi);
int BndStaticVN(const cGH *GH, const int *stencil, const char *vname);

/* Radiative boundaries */

/* prototype for routine registed as providing 'Radiative' boundary conditions
 */
CCTK_INT BndRadiative(const cGH *GH, CCTK_INT num_vars, CCTK_INT *var_indicies,
                      CCTK_INT *faces, CCTK_INT *boundary_widths,
                      CCTK_INT *table_handle);

int BndRadiativeDirGI(const cGH *GH, int stencil_size, int dir, CCTK_REAL var0,
                      CCTK_REAL v0, int gi, int gi_p);
int BndRadiativeDirGN(const cGH *GH, int stencil_size, int dir, CCTK_REAL var0,
                      CCTK_REAL v0, const char *gname_to,
                      const char *gname_from);
int BndRadiativeDirVI(const cGH *GH, int stencil_size, int dir, CCTK_REAL var0,
                      CCTK_REAL v0, int vi, int vi_p);
int BndRadiativeDirVN(const cGH *GH, int stencil_size, int dir, CCTK_REAL var0,
                      CCTK_REAL v0, const char *vname_to,
                      const char *vname_from);

int BndRadiativeGI(const cGH *GH, const int *stencil, CCTK_REAL var0,
                   CCTK_REAL v0, int gi, int gi_p);
int BndRadiativeGN(const cGH *GH, const int *stencil, CCTK_REAL var0,
                   CCTK_REAL v0, const char *gname_to, const char *gname_from);
int BndRadiativeVI(const cGH *GH, const int *stencil, CCTK_REAL var0,
                   CCTK_REAL v0, int vi, int vi_p);
int BndRadiativeVN(const cGH *GH, const int *stencil, CCTK_REAL var0,
                   CCTK_REAL v0, const char *vname_to, const char *vname_from);

/* Robin boundaries */

/* prototype for routine registed as providing 'Robin' boundary condition */
CCTK_INT BndRobin(const cGH *GH, CCTK_INT num_vars, CCTK_INT *var_indicies,
                  CCTK_INT *faces, CCTK_INT *boundary_widths,
                  CCTK_INT *table_handle);

int BndRobinGI(const cGH *GH, const int *stencil, CCTK_REAL finf, int npow,
               int gi);
int BndRobinGN(const cGH *GH, const int *stencil, CCTK_REAL finf, int npow,
               const char *gname);
int BndRobinVI(const cGH *GH, const int *stencil, CCTK_REAL finf, int npow,
               int vi);
int BndRobinVN(const cGH *GH, const int *stencil, CCTK_REAL finf, int npow,
               const char *vname);

/* Flat boundaries */

/* prototype for routine registed as providing 'Flat' boundary condition */
CCTK_INT BndFlat(const cGH *GH, CCTK_INT num_vars, CCTK_INT *var_indicies,
                 CCTK_INT *faces, CCTK_INT *boundary_widths,
                 CCTK_INT *table_handle);

int BndFlatDirGI(const cGH *GH, int stencil_size, int dir, int gi);
int BndFlatDirGN(const cGH *GH, int stencil_size, int dir, const char *gname);
int BndFlatDirVI(const cGH *GH, int stencil_size, int dir, int vi);
int BndFlatDirVN(const cGH *GH, int stencil_size, int dir, const char *vname);

int BndFlatGI(const cGH *GH, const int *stencil, int gi);
int BndFlatGN(const cGH *GH, const int *stencil, const char *gname);
int BndFlatVI(const cGH *GH, const int *stencil, int vi);
int BndFlatVN(const cGH *GH, const int *stencil, const char *vname);

#ifdef __cplusplus
}
#endif

#endif /* _BOUNDARY_H_ */
