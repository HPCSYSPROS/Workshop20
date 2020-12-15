#include <cctk.h>

#include "fasterp.hh"

namespace CarpetInterp2 {

extern "C" CCTK_INT CarpetInterp2_InterpGridArrays(
    CCTK_POINTER_TO_CONST const cctkGH_, CCTK_INT const N_dims,
    CCTK_INT const order, CCTK_INT const N_interp_points,
    CCTK_POINTER_TO_CONST const interp_coords_, CCTK_INT const N_input_arrays,
    CCTK_INT const *const input_array_indices, CCTK_INT const N_output_arrays,
    CCTK_POINTER const output_arrays_) {
  // Check input values and convert types
  cGH const *const cctkGH = static_cast<cGH const *>(cctkGH_);
  assert(cctkGH);

  assert(N_dims == dim);

  assert(N_interp_points >= 0);

  CCTK_REAL const *const *const interp_coords =
      static_cast<CCTK_REAL const *const *>(interp_coords_);
  assert(interp_coords);
  for (int d = 0; d < dim; ++d) {
    assert(N_interp_points == 0 or interp_coords[d]);
  }

  assert(N_input_arrays >= 0);

  assert(input_array_indices);
  for (int n = 0; n < N_input_arrays; ++n) {
    assert(input_array_indices[n] >= 0 and
           input_array_indices[n] < CCTK_NumVars());
  }

  assert(N_output_arrays >= 0);
  assert(N_output_arrays == N_input_arrays);

  CCTK_REAL *const *const output_arrays =
      static_cast<CCTK_REAL *const *>(output_arrays_);
  assert(output_arrays);
  for (int n = 0; n < N_output_arrays; ++n) {
    assert(N_output_arrays == 0 or output_arrays[n]);
  }

  // Set up interpolation
  fasterp_glocs_t locations(N_interp_points);
  for (int d = 0; d < dim; ++d) {
    for (int i = 0; i < N_interp_points; ++i) {
      locations.coords[d].AT(i) = interp_coords[d][i];
    }
  }
  fasterp_setup_t const setup(cctkGH, locations, order);

  // Interpolate
  vector<int> varinds(N_input_arrays);
  for (int n = 0; n < N_input_arrays; ++n) {
    varinds.AT(n) = input_array_indices[n];
  }
  vector<CCTK_REAL *> values(N_input_arrays);
  for (int n = 0; n < N_input_arrays; ++n) {
    values.AT(n) = output_arrays[n];
  }
  setup.interpolate(cctkGH, varinds, values);

  // Done
  return 0;
}

extern "C" CCTK_INT CarpetInterp2_MultiPatchInterpGridArrays(
    CCTK_POINTER_TO_CONST const cctkGH_, CCTK_INT const N_dims,
    CCTK_INT const order, CCTK_INT const N_interp_points,
    CCTK_INT const *const interp_maps,
    CCTK_POINTER_TO_CONST const interp_coords_, CCTK_INT const N_input_arrays,
    CCTK_INT const *const input_array_indices, CCTK_INT const N_output_arrays,
    CCTK_POINTER const output_arrays_) {
  // Check input values and convert types
  cGH const *const cctkGH = static_cast<cGH const *>(cctkGH_);
  assert(cctkGH);

  assert(N_dims == dim);

  assert(N_interp_points >= 0);

  assert(interp_maps);
  CCTK_REAL const *const *const interp_coords =
      static_cast<CCTK_REAL const *const *>(interp_coords_);
  assert(interp_coords);
  for (int d = 0; d < dim; ++d) {
    assert(N_interp_points == 0 or interp_coords[d]);
  }

  assert(N_input_arrays >= 0);

  assert(input_array_indices);
  for (int n = 0; n < N_input_arrays; ++n) {
    assert(input_array_indices[n] >= 0 and
           input_array_indices[n] < CCTK_NumVars());
  }

  assert(N_output_arrays >= 0);
  assert(N_output_arrays == N_input_arrays);

  CCTK_REAL *const *const output_arrays =
      static_cast<CCTK_REAL *const *>(output_arrays_);
  assert(output_arrays);
  for (int n = 0; n < N_output_arrays; ++n) {
    assert(N_output_arrays == 0 or output_arrays[n]);
  }

  // Set up interpolation
  fasterp_llocs_t locations(N_interp_points);
  for (int i = 0; i < N_interp_points; ++i) {
    locations.maps.AT(i) = interp_maps[i];
  }
  for (int d = 0; d < dim; ++d) {
    for (int i = 0; i < N_interp_points; ++i) {
      locations.coords[d].AT(i) = interp_coords[d][i];
    }
  }
  fasterp_setup_t const setup(cctkGH, locations, order);

  // Interpolate
  vector<int> varinds(N_input_arrays);
  for (int n = 0; n < N_input_arrays; ++n) {
    varinds.AT(n) = input_array_indices[n];
  }
  vector<CCTK_REAL *> values(N_input_arrays);
  for (int n = 0; n < N_input_arrays; ++n) {
    values.AT(n) = output_arrays[n];
  }
  setup.interpolate(cctkGH, varinds, values);

  // Done
  return 0;
}

extern "C" CCTK_POINTER CarpetInterp2_Interp2GridArraysSetup(
    CCTK_POINTER_TO_CONST const cctkGH_, CCTK_INT const N_dims,
    CCTK_INT const order, CCTK_INT const N_interp_points,
    CCTK_POINTER_TO_CONST const interp_coords_) {
  // Check input values and convert types
  cGH const *const cctkGH = static_cast<cGH const *>(cctkGH_);
  assert(cctkGH);

  assert(N_dims == dim);

  assert(N_interp_points >= 0);

  CCTK_REAL const *const *const interp_coords =
      static_cast<CCTK_REAL const *const *>(interp_coords_);
  assert(interp_coords);
  for (int d = 0; d < dim; ++d) {
    assert(N_interp_points == 0 or interp_coords[d]);
  }

  // Set up interpolation
  fasterp_glocs_t locations(N_interp_points);
  for (int d = 0; d < dim; ++d) {
    for (int i = 0; i < N_interp_points; ++i) {
      locations.coords[d].AT(i) = interp_coords[d][i];
    }
  }
  fasterp_setup_t *const setup = new fasterp_setup_t(cctkGH, locations, order);

  // Done
  return setup;
}

extern "C" CCTK_INT CarpetInterp2_Interp2GridArrays(
    CCTK_POINTER_TO_CONST const cctkGH_, CCTK_POINTER_TO_CONST const setup_,
    CCTK_INT const N_input_arrays, CCTK_INT const *const input_array_indices,
    CCTK_INT const N_output_arrays, CCTK_POINTER const output_arrays_) {
  // Check input values and convert types
  cGH const *const cctkGH = static_cast<cGH const *>(cctkGH_);
  assert(cctkGH);

  fasterp_setup_t const *const setup =
      static_cast<fasterp_setup_t const *>(setup_);
  assert(setup);

  assert(N_input_arrays >= 0);

  assert(input_array_indices);
  for (int n = 0; n < N_input_arrays; ++n) {
    assert(input_array_indices[n] >= 0 and
           input_array_indices[n] < CCTK_NumVars());
  }

  assert(N_output_arrays >= 0);
  assert(N_output_arrays == N_input_arrays);

  CCTK_REAL *const *const output_arrays =
      static_cast<CCTK_REAL *const *>(output_arrays_);
  assert(output_arrays);
  for (int n = 0; n < N_output_arrays; ++n) {
    assert(N_output_arrays == 0 or output_arrays[n]);
  }

  // Interpolate
  vector<int> varinds(N_input_arrays);
  for (int n = 0; n < N_input_arrays; ++n) {
    varinds.AT(n) = input_array_indices[n];
  }
  vector<CCTK_REAL *> values(N_input_arrays);
  for (int n = 0; n < N_input_arrays; ++n) {
    values.AT(n) = output_arrays[n];
  }
  setup->interpolate(cctkGH, varinds, values);

  // Done
  return 0;
}

extern "C" CCTK_INT
CarpetInterp2_Interp2GridArraysFree(CCTK_POINTER_TO_CONST const cctkGH_,
                                    CCTK_POINTER const setup_) {
  // Check input values and convert types
  cGH const *const cctkGH = static_cast<cGH const *>(cctkGH_);
  assert(cctkGH);

  fasterp_setup_t const *const setup =
      static_cast<fasterp_setup_t const *>(setup_);
  assert(setup);

  delete setup;

  // Done
  return 0;
}

} // namespace CarpetInterp2
