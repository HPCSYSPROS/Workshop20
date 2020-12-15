/*@@
  @file      SymBase.h
  @author    Erik Schnetter
  @date      2004/03/07 09:48:53
  @desc
             Declarations for thorn SymBase
  @version   $Header$
  @enddesc
@@*/

#ifndef SYMBASE_H
#define SYMBASE_H

#include "cctk.h"
#include "cctk_Arguments.h"

/* SymBase's GH extension */
struct SymBase {
  /* The table handles below identify tables that have the following
     entries:

     CCTK_INT symmetry_handle[2*dim]
     CCTK_INT symmetry_zone_width[2*dim]
     CCTK_FPOINTER symmetry_interpolate[2*dim]

     dim is here the dimension of the group containing the variable,
     i.e. there is one entry per face.  The true type of the
     symmetry_interpolate entries is SymmetryInterpolateFPointer.

     (Given that user code is not allowed to modify these tables, and
     that this structure is defined here in this header file, there is
     no real need to use tables for this.  It would be possible to put
     these data directly into this structure, thereby simplifying user
     code and not make anything less flexible.)
  */

  /* Grid symmetry table handle */
  int sym_table;

  /* Grid array symmetry table handles, one handle per grid
     variable */
  int *array_sym_tables;
};

typedef CCTK_INT (*SymmetryInterpolateFPointer)(
    CCTK_POINTER_TO_CONST const cctkGH_, CCTK_INT const N_dims,
    CCTK_INT const local_interp_handle, CCTK_INT const param_table_handle,
    CCTK_INT const coord_system_handle, CCTK_INT const N_interp_points,
    CCTK_INT const interp_coords_type,
    CCTK_POINTER_TO_CONST const interp_coords[], CCTK_INT const N_input_arrays,
    CCTK_INT const input_array_indices[], CCTK_INT const N_output_arrays,
    CCTK_INT const output_array_types[], CCTK_POINTER const output_arrays[],
    CCTK_INT const faces);

/* Number of registered symmetries */
extern int SymBase_num_symmetries;

/* The names of these symmetries */
extern const char **SymBase_symmetry_names;

/* Startup.c */
int SymBase_Startup(void);
void *SymBase_Setup(tFleshConfig *const config, int const convlev,
                    cGH *const cctkGH);

/* Handles.c */
CCTK_INT SymBase_SymmetryRegister(CCTK_STRING const sym_name);
CCTK_INT SymBase_SymmetryHandleOfName(CCTK_STRING const sym_name);
CCTK_POINTER_TO_CONST SymBase_SymmetryNameOfHandle(CCTK_INT const sym_handle);

/* Faces.c */
CCTK_INT
SymBase_SymmetryRegisterFaces(CCTK_INT const sym_table,
                              CCTK_INT const group_dim,
                              CCTK_INT const sym_handle,
                              CCTK_INT const *const which_faces,
                              CCTK_INT const *const new_symmetry_zone_width);
CCTK_INT
SymBase_SymmetryRegisterGrid(CCTK_POINTER const cctkGH_,
                             CCTK_INT const sym_handle,
                             CCTK_INT const *const which_faces,
                             CCTK_INT const *const new_symmetry_zone_width);
CCTK_INT
SymBase_SymmetryRegisterGI(CCTK_POINTER const cctkGH_,
                           CCTK_INT const sym_handle,
                           CCTK_INT const *const which_faces,
                           CCTK_INT const *const new_symmetry_zone_width,
                           CCTK_INT const group_index);
CCTK_INT
SymBase_SymmetryRegisterGN(CCTK_POINTER const cctkGH_,
                           CCTK_INT const sym_handle,
                           CCTK_INT const *const which_faces,
                           CCTK_INT const *const new_symmetry_zone_width,
                           CCTK_STRING const group_name);

CCTK_INT
SymBase_SymmetryRegisterInterpolatorFaces(
    CCTK_INT const sym_table, CCTK_INT const group_dim,
    CCTK_INT const sym_handle,
    SymmetryInterpolateFPointer const new_symmetry_interpolate);
CCTK_INT
SymBase_SymmetryRegisterGridInterpolator(
    CCTK_POINTER const cctkGH_, CCTK_INT const sym_handle,
    SymmetryInterpolateFPointer const new_symmetry_interpolate);

/* Table.c */
CCTK_INT
SymBase_SymmetryTableHandleForGrid(CCTK_POINTER_TO_CONST const cctkGH_);
CCTK_INT
SymBase_SymmetryTableHandleForGI(CCTK_POINTER_TO_CONST const cctkGH_,
                                 CCTK_INT const group_index);
CCTK_INT
SymBase_SymmetryTableHandleForGN(CCTK_POINTER_TO_CONST const cctkGH_,
                                 CCTK_STRING const group_name);

/* Interpolate.c */
CCTK_INT
SymBase_SymmetryInterpolate(
    CCTK_POINTER_TO_CONST const cctkGH_, CCTK_INT const N_dims,
    CCTK_INT const local_interp_handle, CCTK_INT const param_table_handle,
    CCTK_INT const coord_system_handle, CCTK_INT const N_interp_points,
    CCTK_INT const interp_coords_type,
    CCTK_POINTER_TO_CONST const interp_coords[], CCTK_INT const N_input_arrays,
    CCTK_INT const input_array_indices[], CCTK_INT const N_output_arrays,
    CCTK_INT const output_array_types[], CCTK_POINTER const output_arrays[]);
CCTK_INT
SymBase_SymmetryInterpolateFaces(
    CCTK_POINTER_TO_CONST const cctkGH_, CCTK_INT const N_dims,
    CCTK_INT const local_interp_handle, CCTK_INT const param_table_handle,
    CCTK_INT const coord_system_handle, CCTK_INT const N_interp_points,
    CCTK_INT const interp_coords_type,
    CCTK_POINTER_TO_CONST const interp_coords[], CCTK_INT const N_input_arrays,
    CCTK_INT const input_array_indices[], CCTK_INT const N_output_arrays,
    CCTK_INT const output_array_types[], CCTK_POINTER const output_arrays[],
    CCTK_INT const faces);

/* Statistics.c */
void SymBase_Statistics(CCTK_ARGUMENTS);

/* Check.c */
void SymBase_Check(CCTK_ARGUMENTS);

#endif /* ! defined SYMBASE_H */
