#ifndef CARPETINTERP_HH
#define CARPETINTERP_HH

#include "cctk.h"

#include "interp.h"

namespace CarpetInterp {

int InterpGridArrays(
    cGH const *const cGH, int const N_dims, int const local_interp_handle,
    int const param_table_handle, int const coord_system_handle,
    int const N_interp_points, int const interp_coords_type_code,
    void const *const interp_coords[], int const N_input_arrays,
    CCTK_INT const input_array_variable_indices[], int const N_output_arrays,
    CCTK_INT const output_array_type_codes[], void *const output_arrays[]);

} // namespace CarpetInterp

#endif // !defined(CARPETINTERP_HH)
