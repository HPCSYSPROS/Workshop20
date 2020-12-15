#include <cctk.h>
#include <cctk_Parameters.h>
#include <util_ErrorCodes.h>
#include <util_Table.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <map>
#include <vector>

#ifdef CCTK_MPI
#include <mpi.h>
#else
#include "nompi.h"
#endif

#include "bbox.hh"
#include "data.hh"
#include "defs.hh"
#include "dist.hh"
#include "ggf.hh"
#include "timestat.hh"
#include "typeprops.hh"
#include "vect.hh"

#include "carpet.hh"

#include "interp.hh"

//////////////////////////////////////////////////////////////////////////
// For a general overview on the implementation
// please see CarpetInterp's thorn documentation.
//////////////////////////////////////////////////////////////////////////

namespace CarpetInterp {

using namespace std;
using namespace Carpet;

// Return a unique index for component c on reflevel rl, map m and processor p
#define component_idx(p, m, rl, c)                                             \
  component_idx_(p, m, rl, c, minrl, maxrl, maxncomps)
static inline int component_idx_(int const p, int const m, int const rl,
                                 int const c, int const minrl, int const maxrl,
                                 int const maxncomps) {
  assert(p >= 0 and p < dist::size());
  assert(m >= 0 and m < maps);
  assert(rl >= minrl and rl < maxrl);
  assert(c >= 0 and c < maxncomps);
  int const local_idx = ((rl - minrl) * maps + m) * maxncomps + c;
  int const global_idx = p * (maxrl - minrl) * maps * maxncomps + local_idx;
  return global_idx;
}

static int extract_parameter_table_options(
    cGH const *const cctkGH, int const param_table_handle,
    int const N_interp_points, int const N_input_arrays,
    int const N_output_arrays, bool &want_global_mode, bool &have_source_map,
    vector<int> &num_time_derivs, int &prolongation_order_time,
    CCTK_REAL &current_time, CCTK_REAL &delta_time,
    vector<CCTK_INT> &source_map, vector<CCTK_INT> &operand_indices,
    vector<CCTK_INT> &time_deriv_order);

static void map_points(cGH const *const cctkGH, int const coord_system_handle,
                       int const coord_group, int const ml, int const minrl,
                       int const maxrl, int const maxncomps, int const N_dims,
                       int const ndims, int const npoints,
                       vector<CCTK_INT> &source_map,
                       void const *const coords_list[],
                       CCTK_REAL const *const coords, vector<int> &procs,
                       vector<int> &sendcnt, vector<int> &rlev,
                       vector<int> &home, std::map<int, int> &homecntsmap,
                       vector<int> &homecnts);

static void interpolate_components(
    cGH const *const cctkGH, int const coord_system_handle,
    int const coord_group, int const minrl, int const maxrl,
    int const maxncomps, bool const want_global_mode,
    int const prolongation_order_time, int const N_dims,
    vector<int> const &homecnts, std::map<int, int> const &homecntsmap,
    vector<int> const &recvcnt, vector<CCTK_REAL *> const &coords,
    vector<char *> const &outputs, CCTK_INT *const per_proc_statuses,
    CCTK_INT *const per_proc_retvals, vector<CCTK_INT> const &operand_indices,
    vector<CCTK_INT> const &time_deriv_order,
    vector<int> const &num_time_derivs, CCTK_INT const local_interp_handle,
    CCTK_INT const param_table_handle, CCTK_REAL const current_time,
    CCTK_REAL const delta_time, int const N_input_arrays,
    int const N_output_arrays, CCTK_INT const output_array_type_codes[],
    CCTK_INT const input_array_variable_indices[]);

static void interpolate_single_component(
    cGH const *const cctkGH, int const coord_system_handle,
    int const coord_group, int const N_dims, int const npoints,
    CCTK_REAL const *const coords, char *const outputs,
    CCTK_INT &overall_status, CCTK_INT &overall_retval,
    vector<CCTK_INT> const &operand_indices,
    vector<CCTK_INT> const &time_deriv_order,
    vector<CCTK_INT> const &interp_num_time_levels,
    CCTK_INT const local_interp_handle, CCTK_INT const param_table_handle,
    int rl, int m, int lc, vector<int> const &num_tl,
    vector<bool> const &need_time_interp, CCTK_REAL const current_time,
    CCTK_REAL const delta_time, int const N_input_arrays,
    int const N_output_arrays, CCTK_INT const output_array_type_codes[],
    CCTK_INT const input_array_variable_indices[]);

int CarpetInterpStartup() {
  CCTK_OverloadInterpGridArrays(InterpGridArrays);
  return 0;
}

int InterpGridArrays(
    cGH const *const cctkGH, int const N_dims, int const local_interp_handle,
    int const param_table_handle, int const coord_system_handle,
    int const N_interp_points, int const interp_coords_type_code,
    void const *const coords[], int const N_input_arrays,
    CCTK_INT const input_array_variable_indices[], int const N_output_arrays,
    CCTK_INT const output_array_type_codes[], void *const output_arrays[]) {
  if (CCTK_IsFunctionAliased("SymmetryInterpolate")) {
    return SymmetryInterpolate(
        cctkGH, N_dims, local_interp_handle, param_table_handle,
        coord_system_handle, N_interp_points, interp_coords_type_code, coords,
        N_input_arrays, input_array_variable_indices, N_output_arrays,
        output_array_type_codes, output_arrays);
  } else {
    return Carpet_DriverInterpolate(
        cctkGH, N_dims, local_interp_handle, param_table_handle,
        coord_system_handle, N_interp_points, interp_coords_type_code, coords,
        N_input_arrays, input_array_variable_indices, N_output_arrays,
        output_array_type_codes, output_arrays);
  }
}

extern "C" CCTK_INT Carpet_DriverInterpolate(
    CCTK_POINTER_TO_CONST const cctkGH_, CCTK_INT const N_dims,
    CCTK_INT const local_interp_handle, CCTK_INT const param_table_handle,
    CCTK_INT const coord_system_handle, CCTK_INT const N_interp_points,
    CCTK_INT const interp_coords_type_code,
    CCTK_POINTER_TO_CONST const coords_list[], CCTK_INT const N_input_arrays,
    CCTK_INT const input_array_variable_indices[],
    CCTK_INT const N_output_arrays, CCTK_INT const output_array_type_codes[],
    CCTK_POINTER const output_arrays[]) {
  DECLARE_CCTK_PARAMETERS;

  static Timer *timer_CDI = NULL;
  if (not timer_CDI) {
    timer_CDI = new Timer("CarpetInterp::Carpet_DriverInterpolate");
  }
  timer_CDI->start();

  cGH const *const cctkGH = static_cast<cGH const *>(cctkGH_);
  assert(cctkGH);
  assert(0 <= N_dims and N_dims <= dim);

  // Check input arrays
  int coord_group = -1;
  cGroupDynamicData coord_group_data;
  int real_N_input_arrays = 0;
  for (int n = 0; n < N_input_arrays; n++) {

    // Negative indices are ignored
    const int vindex = input_array_variable_indices[n];
    if (vindex < 0) {
      continue;
    }
    ++real_N_input_arrays;

    const int group = CCTK_GroupIndexFromVarI(vindex);
    if (group < 0) {
      CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "input array variable %d = %d is not a valid "
                 "variable index",
                 n, vindex);
    }
    const int gtype = CCTK_GroupTypeI(group);
    if (gtype != CCTK_GF and gtype != CCTK_ARRAY) {
      CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "input array variable %d is not of type CCTK_GF or "
                 "CCTK_ARRY",
                 n);
    }
    if (coord_group < 0) {
      coord_group = group;
      CCTK_GroupDynamicData(cctkGH, coord_group, &coord_group_data);
    } else {
      // Check that group and coord_group have the same layout
      cGroupDynamicData gdata;
      CCTK_GroupDynamicData(cctkGH, group, &gdata);
      const int size = gdata.dim * sizeof(int);
      if (gdata.dim != coord_group_data.dim or
          memcmp(gdata.lsh, coord_group_data.lsh, size) or
          memcmp(gdata.lbnd, coord_group_data.lbnd, size) or
          memcmp(gdata.ubnd, coord_group_data.ubnd, size) or
          memcmp(gdata.bbox, coord_group_data.bbox, 2 * size) or
          memcmp(gdata.nghostzones, coord_group_data.nghostzones, size)) {
        CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "input array variable %d has different layout than "
                   "the underlying coordinate system",
                   n);
      }
    }
  }
  if (real_N_input_arrays == 0) {
    // When there are no interpolation variables, use a pseudo group
    coord_group = CCTK_GroupIndex("grid::coordinates");
  }
  assert(coord_group >= 0);

  // Check output arrays
  assert(N_output_arrays > 0);
  const int output_array_type = output_array_type_codes[0];
  for (int n = 1; n < N_output_arrays; n++) {
    if (output_array_type != output_array_type_codes[n]) {
      CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Currently all output arrays have to have the same datatype. "
                 "Array 0 has type '%s' but array %d has type '%s'",
                 CCTK_VarTypeName(output_array_type), n,
                 CCTK_VarTypeName(output_array_type_codes[n]));
    }
  }

  if (is_meta_mode()) {
    CCTK_WARN(CCTK_WARN_ABORT,
              "It is not possible to interpolate in meta mode");
  }

  // Multiple convergence levels are not supported
  assert(mglevels == 1);
  int const ml = 0;

  assert(N_interp_points >= 0);
  assert(coords_list);
  for (int d = 0; d < N_dims; ++d) {
    assert(N_interp_points == 0 or coords_list[d]);
  }

  if (interp_coords_type_code != CCTK_VARIABLE_REAL) {
    CCTK_WARN(CCTK_WARN_ABORT,
              "CarpetInterp does not support interpolation "
              "coordinates other than datatype CCTK_VARIABLE_REAL");
  }

  assert(N_output_arrays >= 0);
  if (N_interp_points > 0) {
    assert(output_arrays);
    for (int j = 0; j < N_output_arrays; ++j) {
      assert(output_arrays[j]);
      for (int jj = 0; jj < j; ++jj) {
        assert(output_arrays[j] != output_arrays[jj]);
      }
    }
  }

  if (barriers) {
    Carpet::NamedBarrier(cctkGH, 696681976,
                         "CarpetInterp::Carpet_DriverInterpolate");
  }

  //////////////////////////////////////////////////////////////////////
  // Extract parameter table options:
  //   - source map
  //   - output array operand indices
  //   - time interpolation order
  //////////////////////////////////////////////////////////////////////
  bool want_global_mode;
  vector<CCTK_INT> source_map(N_interp_points);
  vector<CCTK_INT> operand_indices(N_output_arrays);
  vector<CCTK_INT> time_deriv_order(N_output_arrays);
  bool have_source_map;
  vector<int> num_time_derivs;
  CCTK_REAL current_time, delta_time;
  int prolongation_order_time;

  {
    int const iret = extract_parameter_table_options(
        cctkGH, param_table_handle, N_interp_points, N_input_arrays,
        N_output_arrays, want_global_mode, have_source_map, num_time_derivs,
        prolongation_order_time, current_time, delta_time, source_map,
        operand_indices, time_deriv_order);
    if (iret < 0) {
      timer_CDI->stop(0);
      return iret;
    }
  }

  // Find range of refinement levels
  assert(maps > 0);
  for (int m = 1; m < maps; ++m) {
    assert(arrdata.AT(coord_group).AT(0).hh->reflevels() ==
           arrdata.AT(coord_group).AT(m).hh->reflevels());
  }
  int const minrl = want_global_mode ? 0 : reflevel;
  int const maxrl = want_global_mode
                        ? arrdata.AT(coord_group).AT(0).hh->reflevels()
                        : reflevel + 1;

  // Find maximum number of components over all levels and maps
  int maxncomps = 0;
  for (int rl = minrl; rl < maxrl; ++rl) {
    for (int m = 0; m < maps; ++m) {
      maxncomps =
          max(maxncomps, arrdata.AT(coord_group).AT(m).hh->components(rl));
    }
  }

  //////////////////////////////////////////////////////////////////////
  // Map interpolation points to processors
  //////////////////////////////////////////////////////////////////////
  vector<int> sendcnt(dist::size());
  vector<int> dstprocs(N_interp_points); // which processor owns point n
  vector<int> rlev(N_interp_points);     // refinement level of point n
  vector<int> home(N_interp_points);     // component of point n
  std::map<int, int> homecntsmap;        // components hash map
  vector<int> allhomecnts;               // number of points in component
                                         // homecntsmap.find(c)
  int const ndims = have_source_map ? N_dims + 1 : N_dims;

  // Each point from coord_list is mapped onto the processor
  // that owns it (dstprocs)
  // Also accumulate the number of points per processor (sendcnt)
  map_points(cctkGH, coord_system_handle, coord_group, ml, minrl, maxrl,
             maxncomps, N_dims, ndims, N_interp_points, source_map, coords_list,
             NULL, dstprocs, sendcnt, rlev, home, homecntsmap, allhomecnts);
  // cout << "CarpetInterp: sendcnt=" << sendcnt << endl;

  //////////////////////////////////////////////////////////////////////
  // Communicate the number of points each processor is going to communicate
  //////////////////////////////////////////////////////////////////////

  // recvcnt denotes the number of points
  // that this processor is to receive from others
  vector<int> recvcnt(dist::size());
  {
    static Timer *timer = NULL;
    if (not timer) {
      timer = new Timer("CarpetInterp::send_npoints");
    }
    timer->start();
    MPI_Alltoall(&sendcnt[0], 1, dist::mpi_datatype(sendcnt[0]), &recvcnt[0], 1,
                 dist::mpi_datatype(recvcnt[0]), dist::comm());
    timer->stop(dist::size() * sizeof(CCTK_INT));
  }
  // cout << "CarpetInterp: recvcnt=" << recvcnt << endl;

  //////////////////////////////////////////////////////////////////////
  // Communicate the interpolation coordinates
  //////////////////////////////////////////////////////////////////////

  // Set up counts and displacements for MPI_Alltoallv()
  // N_points_local is the total number of points to receive
  // and thus the total number of points to interpolate on this processor
  int N_points_local = recvcnt.AT(0);
  vector<int> senddispl(dist::size());
  vector<int> recvdispl(dist::size());
  for (int p = 1; p < dist::size(); p++) {
    N_points_local += recvcnt.AT(p);
    senddispl.AT(p) = senddispl.AT(p - 1) + sendcnt.AT(p - 1);
    recvdispl.AT(p) = recvdispl.AT(p - 1) + recvcnt.AT(p - 1);
  }

  // remember the position of each point in the original input arrays
  vector<int> indices(N_interp_points);
  {
    // totalhomecnts is the accumulated number of points over all components
    vector<int> totalhomecnts(allhomecnts.size());
    if (totalhomecnts.size() > 0) {
      totalhomecnts.AT(0) = 0;
      for (size_t c = 1; c < totalhomecnts.size(); c++) {
        totalhomecnts.AT(c) = totalhomecnts.AT(c - 1) + allhomecnts.AT(c - 1);
      }
    }

    vector<int> tmpcnts(allhomecnts.size());
#pragma omp parallel for
    for (int n = 0; n < N_interp_points; n++) {
      int const cidx = component_idx(dstprocs.AT(n), source_map.AT(n),
                                     rlev.AT(n), home.AT(n));
      std::map<int, int>::const_iterator it = homecntsmap.find(cidx);
      assert(it != homecntsmap.end());
      int const idx = it->second;
      assert(idx < (int)totalhomecnts.size());
      int mytmpcnt;
#pragma omp critical
      { mytmpcnt = tmpcnts.AT(idx)++; }
      indices.AT(n) = totalhomecnts.AT(idx) + mytmpcnt;
    }
    assert(tmpcnts == allhomecnts);
  }

  // Allocate the communication send buffer
  // and copy the input coordinates and the source map (if necessary) into it
  vector<CCTK_REAL> coords_buffer(ndims * N_interp_points);
#pragma omp parallel for
  for (int n = 0; n < N_interp_points; n++) {
    int const idx = indices.AT(n);
    assert((idx + 1) * ndims <= (int)coords_buffer.size());
    for (int d = 0; d < N_dims; d++) {
      coords_buffer[d + ndims * idx] =
          static_cast<CCTK_REAL const *>(coords_list[d])[n];
    }
    if (have_source_map) {
      coords_buffer[N_dims + ndims * idx] = source_map.AT(n);
    }
  }

  // Send this processor's points to owning processors,
  // receive other processors' points for interpolation on this processor
  {
    vector<CCTK_REAL> tmp(ndims * N_points_local);
    {
      int const sendbufsize = (int)coords_buffer.size();
      int const recvbufsize = (int)tmp.size();
      for (int n = 0; n < dist::size(); ++n) {
        assert(sendcnt.AT(n) >= 0);
        assert(senddispl.AT(n) >= 0);
        assert(senddispl.AT(n) + sendcnt.AT(n) <= sendbufsize);
        assert(recvcnt.AT(n) >= 0);
        assert(recvdispl.AT(n) >= 0);
        assert(recvdispl.AT(n) + recvcnt.AT(n) <= recvbufsize);
      }
    }
#ifndef _NDEBUG
#pragma omp parallel for
    for (int i = 0; i < (int)tmp.size(); ++i) {
      tmp.AT(i) = poison;
    }
#endif
    {
      MPI_Datatype const datatype = dist::mpi_datatype(tmp[0]);
      MPI_Datatype vdatatype;
      MPI_Type_vector(1, ndims, 0, datatype, &vdatatype);
      MPI_Type_commit(&vdatatype);

      static Timer *timer = NULL;
      if (not timer) {
        timer = new Timer("CarpetInterp::send_coordinates");
      }
      timer->start();
      MPI_Alltoallv(&coords_buffer[0], &sendcnt[0], &senddispl[0], vdatatype,
                    &tmp[0], &recvcnt[0], &recvdispl[0], vdatatype,
                    dist::comm());
      timer->stop(coords_buffer.size() * sizeof(CCTK_REAL));

      MPI_Type_free(&vdatatype);
    }
#ifndef _NDEBUG
    {
      vector<bool> filled(N_points_local, false);
      for (int n = 0; n < (int)dist::size(); ++n) {
        //#pragma omp parallel for
        for (int i = 0; i < recvcnt.AT(n); ++i) {
          assert(not filled.AT(recvdispl.AT(n) + i));
          filled.AT(recvdispl.AT(n) + i) = true;
        }
      }
      bool error = false;
      for (int i = 0; i < (int)filled.size(); ++i) {
        error = error or not filled.AT(i);
      }
      if (error) {
        cout << "error" << endl;
        cout << "recvdispl: " << recvdispl << endl;
        cout << "recvcnt: " << recvcnt << endl;
        cout << "filled: " << filled << endl;
      }
#pragma omp parallel for
      for (int i = 0; i < (int)filled.size(); ++i) {
        assert(filled.AT(i));
      }
    }
#pragma omp parallel for
    for (int i = 0; i < (int)tmp.size(); ++i) {
      assert(tmp.AT(i) != poison);
    }
#endif
    coords_buffer.swap(tmp);
  }

  //////////////////////////////////////////////////////////////////////
  // Set up (local) source map
  //////////////////////////////////////////////////////////////////////
  if (have_source_map) {
    source_map.resize(N_points_local);
#pragma omp parallel for
    for (int n = 0; n < N_points_local; n++) {
      source_map.AT(n) = static_cast<int>(coords_buffer[ndims * n + N_dims]);
    }
  } else {
    // No explicit source map specified
    if (Carpet::map != -1) {
      // Interpolate from the current map
      source_map.resize(N_points_local, Carpet::map);
    } else {
      // Interpolate from map 0 if this is the only one
      // (for backwards compatibility)
      assert(maps == 1);
      source_map.resize(N_points_local, 0);
    }
  }

  //////////////////////////////////////////////////////////////////////
  // Map (local) interpolation points to components
  //////////////////////////////////////////////////////////////////////
  // Remember from where point n came
  vector<int> srcprocs(N_points_local); // which processor requested point n
  {
    int offset = 0;
    for (int p = 0; p < (int)recvcnt.size(); p++) {
      for (int n = 0; n < recvcnt.AT(p); n++) {
        srcprocs.AT(offset++) = p;
      }
    }
    assert(offset == N_points_local);
  }

  // map each point in coords_buffer onto its component (srcproc, rlev, home)
  // also accumulate the number of points in each component (homecnts)
  //
  // In order to be somewhat more efficient, we could also map
  // all processors' (rlev, home) points onto the same component
  // and call the local interpolator on that (thus saving potentially
  // nprocs-1 calls).
  // The reason why this hasn't been implemented here is because
  // CCTK_InterpGridArrays() is supposed to return a value which is
  // the minimum over all local interpolator invocations' return values
  // for the points requested by this processor. This overall minimum
  // isn't defined anymore if we call the local interpolator on a component
  // now containing points from different processors.
  // (One could argue though that the per-point status codes as returned
  // by the local interpolator could be used to determine the global
  // interpolator return value instead.)
  rlev.resize(N_points_local); // reflevel of (local) point n
  home.resize(N_points_local); // component of (local) point n
  vector<int> homecnts;        // number of points in component
                               // homecntsmap.find(c)
  map_points(cctkGH, coord_system_handle, coord_group, ml, minrl, maxrl,
             maxncomps, N_dims, ndims, N_points_local, source_map, NULL,
             &coords_buffer.front(), srcprocs, sendcnt, rlev, home, homecntsmap,
             homecnts);

  // Free some memory
  source_map.clear();
  srcprocs.clear();

  // Reorder the coordinates from <N_dims>-tupels into <N_dims> vectors
  // as expected by CCTK_InterpLocalUniform()
  {
    int offset = 0;
    vector<CCTK_REAL> tmp(coords_buffer.size());
    for (size_t c = 0; c < homecnts.size(); c++) {
#pragma omp parallel for
      for (int n = 0; n < homecnts.AT(c); n++) {
        for (int d = 0; d < N_dims; d++) {
          tmp.AT(d * homecnts.AT(c) + n + N_dims * offset) =
              coords_buffer.AT((n + offset) * ndims + d);
        }
      }
      offset += homecnts.AT(c);
    }
    assert(offset == N_points_local);
    coords_buffer.swap(tmp);
  }

  //////////////////////////////////////////////////////////////////////
  // Do the local interpolation on individual components
  //////////////////////////////////////////////////////////////////////
  const int vtype = output_array_type_codes[0];
  const int vtypesize = CCTK_VarTypeSize(vtype);
  assert(vtypesize > 0);
  vector<char> outputs_buffer(N_points_local * N_output_arrays * vtypesize);
  vector<char *> outputs(homecnts.size(), &outputs_buffer.front());
  vector<CCTK_REAL *> coords(homecnts.size(), &coords_buffer.front());
  vector<CCTK_INT> status_and_retval_buffer(2 * dist::size(), 0);
  CCTK_INT *per_proc_statuses = &status_and_retval_buffer.front();
  CCTK_INT *per_proc_retvals = per_proc_statuses + dist::size();

  // Set up the per-component coordinates and output arrays as offsets
  // into the single communication buffers
  {
    int offset = 0;
    for (size_t c = 0; c < homecnts.size(); c++) {
      coords[c] += N_dims * offset;
      outputs[c] += N_output_arrays * offset * vtypesize;
      offset += homecnts[c];
    }
    assert(offset == N_points_local);
  }

  interpolate_components(
      cctkGH, coord_system_handle, coord_group, minrl, maxrl, maxncomps,
      want_global_mode, prolongation_order_time, N_dims, homecnts, homecntsmap,
      recvcnt, coords, outputs, per_proc_statuses, per_proc_retvals,
      operand_indices, time_deriv_order, num_time_derivs, local_interp_handle,
      param_table_handle, current_time, delta_time, N_input_arrays,
      N_output_arrays, output_array_type_codes, input_array_variable_indices);

  // Free some memory
  coords_buffer.clear();
  coords.clear();
  homecnts.clear();
  home.clear();
  rlev.clear();

  //////////////////////////////////////////////////////////////////////
  // Communicate interpolation results
  //////////////////////////////////////////////////////////////////////

  {
    vector<char> tmp(N_interp_points * N_output_arrays * vtypesize);

    MPI_Datatype datatype;
    switch (specific_cactus_type(vtype)) {
#define TYPECASE(N, T)                                                         \
  case N: {                                                                    \
    T dummy;                                                                   \
    datatype = dist::mpi_datatype(dummy);                                      \
    break;                                                                     \
  }
#include "typecase.hh"
#undef TYPECASE
    default: {
      CCTK_WARN(CCTK_WARN_ABORT, "invalid datatype");
      abort();
    }
    }
    if (datatype == MPI_DATATYPE_NULL) {
      CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "MPI datatype for Cactus datatype %d is not defined", vtype);
    }
    MPI_Datatype vdatatype;
    MPI_Type_vector(1, N_output_arrays, 0, datatype, &vdatatype);
    MPI_Type_commit(&vdatatype);

    static Timer *timer = NULL;
    if (not timer) {
      timer = new Timer("CarpetInterp::recv_points");
    }
    timer->start();
    // distribute the results the same way as the coordinates were gathered
    // simply by interchanging the send/recv counts/displacements
    MPI_Alltoallv(&outputs_buffer[0], &recvcnt[0], &recvdispl[0], vdatatype,
                  &tmp[0], &sendcnt[0], &senddispl[0], vdatatype, dist::comm());
    timer->stop(N_interp_points * N_output_arrays * vtypesize);

    MPI_Type_free(&vdatatype);

    outputs_buffer.swap(tmp);
  }

  //////////////////////////////////////////////////////////////////////
  // Communicate interpolation status codes and return values
  //////////////////////////////////////////////////////////////////////
  {
    // A processor's overall status and return code
    // is defined as the minimum over all local interpolator status and
    // return codes across all processors for that processor
    vector<CCTK_INT> tmp(status_and_retval_buffer.size());
    { assert(status_and_retval_buffer.size() == tmp.size()); }
    MPI_Allreduce(&status_and_retval_buffer[0], &tmp[0], tmp.size(),
                  dist::mpi_datatype(tmp[0]), MPI_MIN, dist::comm());
    status_and_retval_buffer.swap(tmp);
    per_proc_statuses = &status_and_retval_buffer.front();
    per_proc_retvals = per_proc_statuses + dist::size();
  }

  //////////////////////////////////////////////////////////////////////
  // Finally, sort the received outputs back into the caller's arrays
  //////////////////////////////////////////////////////////////////////

  // Sorting is done with the help of the inverse indices vector
  vector<int> reverse_indices(indices.size());
#pragma omp parallel for
  for (int i = 0; i < (int)indices.size(); i++) {
    reverse_indices[indices[i]] = i;
  }

  for (int d = 0; d < N_output_arrays; d++) {
    char *output_array = static_cast<char *>(output_arrays[d]);
    int offset = 0;
    for (int c = 0, i = 0; c < (int)allhomecnts.size(); c++) {
      assert((int)(allhomecnts.AT(c) * (d + 1) + offset) <=
             N_output_arrays * N_interp_points);
      assert(d >= 0);
#pragma omp parallel for
      for (int n = 0; n < allhomecnts.AT(c); n++) {
        {
          assert(reverse_indices.AT(i + n) >= 0 and
                 reverse_indices.AT(i + n) < N_interp_points);
          assert(allhomecnts.AT(c) >= 0);
          assert(allhomecnts.AT(c) * d + offset + n <
                 (int)outputs_buffer.size() / vtypesize);
        }
        memcpy(output_array + reverse_indices.AT(i + n) * vtypesize,
               &outputs_buffer.front() +
                   (allhomecnts.AT(c) * d + offset + n) * vtypesize,
               vtypesize);
      }
      i += allhomecnts.AT(c);
      offset += N_output_arrays * allhomecnts.AT(c);
    }
    assert(offset == N_output_arrays * N_interp_points);
  }

  // set this processor's overall local interpolator status code
  int ierr =
      Util_TableSetInt(param_table_handle, per_proc_statuses[dist::rank()],
                       "local_interpolator_status");
  assert(ierr >= 0);

  // Done.
  {
    timer_CDI->stop(0);
    int const iret = per_proc_retvals[dist::rank()];
    return iret;
  }
}

static int extract_parameter_table_options(
    cGH const *const cctkGH, int const param_table_handle,
    int const N_interp_points, int const N_input_arrays,
    int const N_output_arrays, bool &want_global_mode, bool &have_source_map,
    vector<int> &num_time_derivs, int &prolongation_order_time,
    CCTK_REAL &current_time, CCTK_REAL &delta_time,
    vector<CCTK_INT> &source_map, vector<CCTK_INT> &operand_indices,
    vector<CCTK_INT> &time_deriv_order) {
  DECLARE_CCTK_PARAMETERS;

  int iret;

  // Do we want to interpolate in global mode, i.e., from all
  // refinement levels?
  CCTK_INT want_global_mode1;
  iret = Util_TableGetInt(param_table_handle, &want_global_mode1,
                          "want_global_mode");
  if (iret == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
    want_global_mode = is_global_mode();
  } else if (iret < 0) {
    CCTK_WARN(CCTK_WARN_ALERT, "internal error");
    return -1;
  } else if (iret != 1) {
    CCTK_WARN(CCTK_WARN_ALERT, "internal error");
    return -1;
  } else {
    want_global_mode = want_global_mode1;
  }

  // Find source map
  assert((int)source_map.size() == N_interp_points);
  iret = Util_TableGetIntArray(param_table_handle, N_interp_points,
                               &source_map.front(), "source_map");
  have_source_map = not(iret == UTIL_ERROR_TABLE_NO_SUCH_KEY);
  if (not have_source_map) {
    // No explicit source map specified
    if (Carpet::map != -1) {
      // Interpolate from the current map
      source_map.assign(source_map.size(), Carpet::map);
    } else if (maps == 1) {
      // Interpolate from map 0 if this is the only one
      // (for backwards compatibility)
      source_map.assign(source_map.size(), 0);
    } else {
      CCTK_WARN(CCTK_WARN_ALERT, "No source map specified");
      return -1;
    }
  } else if (iret < 0) {
    CCTK_WARN(CCTK_WARN_ALERT, "internal error");
    return -1;
  } else if (iret != N_interp_points) {
    CCTK_WARN(CCTK_WARN_ALERT, "Source map array has wrong size");
    return -1;
  } else {
    iret = Util_TableGetIntArray(param_table_handle, source_map.size(),
                                 &source_map.front(), "source_map");
    assert(iret == (int)source_map.size());

#ifndef _NDEBUG
// Check source map
#pragma omp parallel for
    for (int n = 0; n < (int)source_map.size(); ++n) {
      if (not(source_map[n] >= 0 and source_map[n] < maps)) {
        cout << "CI: n=" << n << " map=" << source_map[n] << endl;
      }
      assert(source_map[n] >= 0 and source_map[n] < maps);
    }
#endif
  }

  // Find the time interpolation order
  int partype;
  void const *const parptr =
      CCTK_ParameterGet("prolongation_order_time", "Carpet", &partype);
  assert(parptr);
  assert(partype == PARAMETER_INTEGER);
  prolongation_order_time = *(CCTK_INT const *)parptr;

  current_time = cctkGH->cctk_time;
  delta_time = cctkGH->cctk_delta_time;

  iret = Util_TableGetIntArray(param_table_handle, N_output_arrays,
                               &time_deriv_order.front(), "time_deriv_order");
  if (iret == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
    time_deriv_order.assign(time_deriv_order.size(), 0);
    num_time_derivs.assign(N_output_arrays, 0);
  } else {
    assert(iret == N_output_arrays);
    num_time_derivs.resize(N_output_arrays);
    for (int m = 0; m < N_output_arrays; ++m) {
      num_time_derivs.AT(m) = time_deriv_order[m];
    }
  }

  // Find output variable indices
  iret = Util_TableGetIntArray(param_table_handle, N_output_arrays,
                               &operand_indices.front(), "operand_indices");
  if (iret == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
    assert(N_output_arrays == N_input_arrays);
    for (int m = 0; m < N_output_arrays; ++m) {
      operand_indices[m] = m;
    }
  } else {
    assert(iret == N_output_arrays);
  }

  return 0;
}

// Find the component and integer index to which a grid point
// belongs.  This uses a linear search over all components, which
// does NOT scale with the number of components.
static void
find_location_linear(gh const *restrict const hh, rvect const &restrict pos,
                     rvect const &restrict lower, rvect const &restrict upper,
                     rvect const &restrict delta, int const ml, int const minrl,
                     int const maxrl, int &restrict rl, int &restrict c) {
  // cout << "CarpetInterp: assign: m=" << m << " pos=" << pos << endl;

  assert(ml >= 0 and ml < mglevels);
  assert(minrl >= 0 and minrl < maxrl and maxrl <= reflevels);

  CCTK_REAL const rone = 1.0;
  CCTK_REAL const rhalf = rone / 2;

  if (all(pos >= lower and pos <= upper)) {
    // The point is within the domain

    // Try finer levels first
    for (rl = maxrl - 1; rl >= minrl; --rl) {

      ivect const &istride = hh->baseextent(ml, rl).stride();
      rvect const level_delta =
          delta / rvect(spacereffacts.AT(rl)) * rvect(ipow(mgfact, ml));
      ivect const ipos =
          ivect(floor((pos - lower) / level_delta + rhalf)) * istride;

      if (hh->refcent == cell_centered) {
        assert(all(istride % 2 == 0));
      }

      gh::cregs const &regs = hh->regions.AT(ml).AT(rl);

      // Search all components linearly
      for (c = 0; c < int(regs.size()); ++c) {
        region_t const &reg = regs.AT(c);
        if (reg.extent.contains(ipos)) {
          // We found the refinement level, component, and index to
          // which this grid point belongs
          return;
        }
      }
    }
  }

  // The point does not belong to any component.  This should happen
  // only rarely.
  rl = -1;
  c = -1;
}

// Find the component and integer index to which a grid point
// belongs.  This uses a tree search over the superregions in the
// grid struction, which should scale reasonably (O(n log n)) better
// with the number of componets components.
static void
find_location_tree(gh const *restrict const hh, rvect const &restrict pos,
                   rvect const &restrict lower, rvect const &restrict upper,
                   rvect const &restrict delta, int const ml, int const minrl,
                   int const maxrl, int &restrict rl, int &restrict c) {
  // cout << "CarpetInterp: assign: m=" << m << " pos=" << pos << endl;

  assert(ml >= 0 and ml < mglevels);
  assert(minrl >= 0 and minrl < maxrl and maxrl <= reflevels);

  CCTK_REAL const rone = 1.0;
  CCTK_REAL const rhalf = rone / 2;

  if (all(pos >= lower and pos <= upper)) {
    // The point is within the domain

    // Try finer levels first
    for (rl = maxrl - 1; rl >= minrl; --rl) {

      ivect const &istride = hh->baseextent(ml, rl).stride();
      if (hh->refcent == cell_centered) {
        assert(all(istride % 2 == 0));
      }

      rvect const level_delta =
          delta / rvect(spacereffacts.AT(rl)) * rvect(ipow(mgfact, ml));
      ivect const ipos =
          ivect(floor((pos - lower) / level_delta + rhalf)) * istride;

      gh::cregs const &regs = hh->superregions.AT(rl);

      // Search all superregions linearly.  Each superregion
      // corresponds to a "refined region", and the number of
      // superregions is thus presumably independent of the number
      // of processors.
      for (size_t r = 0; r < regs.size(); ++r) {
        region_t const &reg = regs.AT(r);
        if (reg.extent.contains(ipos)) {
          // We found the superregion to which this grid point
          // belongs

          // Search the superregion hierarchically
          pseudoregion_t const *const preg = reg.processors->search(ipos);
          assert(preg);

          // We now know the refinement level, component, and index
          // to which this grid point belongs
          c = preg->component;
          return;
        }
      }
    }
  }

  // The point does not belong to any component.  This should happen
  // only rarely.
  rl = -1;
  c = -1;
}

static void map_points(cGH const *const cctkGH, int const coord_system_handle,
                       int const coord_group, int const ml, int const minrl,
                       int const maxrl, int const maxncomps, int const N_dims,
                       int const ndims, int const npoints,
                       vector<CCTK_INT> &source_map,
                       void const *const coords_list[],
                       CCTK_REAL const *const coords, vector<int> &procs,
                       vector<int> &sendcnt, vector<int> &rlev,
                       vector<int> &home, std::map<int, int> &homecntsmap,
                       vector<int> &homecnts) {
  DECLARE_CCTK_PARAMETERS;

  static Timer *timer = NULL;
  if (not timer)
    timer = new Timer("CarpetInterp::map_points");
  timer->start();

  assert(npoints == 0 or coords or coords_list);
  bool const map_onto_processors = coords_list != NULL;

  assert((int)procs.size() == npoints);
  assert((int)sendcnt.size() == dist::size());
  assert((int)rlev.size() == npoints);
  assert((int)home.size() == npoints);
  assert((int)source_map.size() == npoints);

  // Find out about the coordinates: origin and delta for the Carpet
  // grid indices
  vector<rvect> lower(maps);
  vector<rvect> upper(maps);
  vector<rvect> delta(maps); // spacing on coarse grid

  int const grouptype = CCTK_GroupTypeI(coord_group);
  switch (grouptype) {

  case CCTK_GF: {
    for (int m = 0; m < Carpet::maps; ++m) {
      jvect gsh;
      GetCoordRange(cctkGH, m, mglevel, dim, &gsh[0], &lower.AT(m)[0],
                    &upper.AT(m)[0], &delta.AT(m)[0]);
    }
    break;
  }

  case CCTK_SCALAR:
  case CCTK_ARRAY: {

    rvect coord_lower, coord_upper;
    char const *const coord_system_name =
        CCTK_CoordSystemName(coord_system_handle);

    for (int d = 0; d < N_dims; ++d) {
      int const iret = CCTK_CoordRange(cctkGH, &coord_lower[d], &coord_upper[d],
                                       d + 1, NULL, coord_system_name);
      assert(iret == 0);
    }

    assert(arrdata.AT(coord_group).size() == 1);
    int const m = 0;
    ibbox const &baseextent =
        arrdata.AT(coord_group).AT(m).hh->baseextents.AT(mglevel).AT(0);
    lower.AT(m) = coord_lower;
    upper.AT(m) = coord_upper;
    delta.AT(m) = ((coord_upper - coord_lower) /
                   rvect(baseextent.upper() - baseextent.lower()));
    break;
  }

  default:
    assert(0);
  }

  // Clear the components hash map
  homecntsmap.clear();

// Assign interpolation points to processors/components
#pragma omp parallel for
  for (int n = 0; n < npoints; ++n) {

    CCTK_INT &m = source_map.AT(n);
    int rl, c = -1;
    rvect pos;
    gh const *hh = NULL;

    if (m >= 0) {

      hh = arrdata.AT(coord_group).AT(m).hh;

      // Find the grid point closest to the interpolation point
      for (int d = 0; d < N_dims; ++d) {
        if (map_onto_processors) {
          pos[d] = static_cast<CCTK_REAL const *>(coords_list[d])[n];
        } else {
          pos[d] = coords[ndims * n + d];
        }
      }

      // Find the component that this grid point belongs to

      // Calculate rl, c, and proc
      if (not tree_search) {
        find_location_linear(hh, pos, lower.AT(m), upper.AT(m), delta.AT(m), ml,
                             minrl, maxrl, rl, c);
      } else {
        find_location_tree(hh, pos, lower.AT(m), upper.AT(m), delta.AT(m), ml,
                           minrl, maxrl, rl, c);
        if (check_tree_search) {
          int rl2, c2;
          find_location_linear(hh, pos, lower.AT(m), upper.AT(m), delta.AT(m),
                               ml, minrl, maxrl, rl2, c2);
          if (rl2 != rl or c2 != c) {
#pragma omp critical
            CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                       "Inconsistent search result from find_location_tree for "
                       "interpolation point #%d at [%g,%g,%g] of patch #%d is "
                       "not on any component",
                       n, (double)pos[0], (double)pos[1], (double)pos[2],
                       (int)m);
          }
        }
      }

    } // if m >= 0

    if (c == -1) {
      // The point could not be mapped onto any component

      // Warn only once, namely when mapping points onto processors.
      // (This routine is called twice; first to map points onto
      // processors, then to map points onto components.)
      if (map_onto_processors) {
#pragma omp critical
        CCTK_VWarn(CCTK_WARN_PICKY, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Interpolation point #%d at [%g,%g,%g] of patch #%d is not "
                   "on any component",
                   n, (double)pos[0], (double)pos[1], (double)pos[2], (int)m);
      }

      // Map the point (arbitrarily) to the first component of the
      // coarsest grid
      // TODO: Handle these points explicitly later on
      rl = minrl;
      c = 0;

      // Find a patch which exists on this processor
      for (m = 0; m < maps; ++m) {
        hh = arrdata.AT(coord_group).AT(m).hh;
        if (hh->components(rl) > c)
          break;
      }
      assert(m < maps);
    }

#ifndef _NDEBUG
    if (not(rl >= minrl and rl < maxrl) or
        not(c >= 0 and c < hh->components(rl))) {
#pragma omp critical
      cout << "CI: m=" << m << " rl=" << rl << " c=" << c
           << " ext=" << hh->extent(ml, rl, c) << endl;
    }
    assert(rl >= minrl and rl < maxrl);
    assert(c >= 0 and c < hh->components(rl));
#endif

    if (map_onto_processors) {
      int const proc = hh->processor(rl, c);
      procs.AT(n) = proc;
      int &this_sendcnt = sendcnt.AT(proc);
#pragma omp atomic
      ++this_sendcnt;
    }
    rlev.AT(n) = rl;
    home.AT(n) = c;
    int const cidx = component_idx(procs.AT(n), m, rl, c);
#pragma omp critical
    {
      // Increase counter, creating a new hash element if necessary
      homecntsmap[cidx]++;
    }

  } // for n

  // allocate and fill the (dense) homecnts vector from the hash map
  homecnts.resize(homecntsmap.size());
  {
    int c = 0;
    for (std::map<int, int>::iterator it = homecntsmap.begin();
         it != homecntsmap.end(); it++) {
      // store the number of points of this component in homecnts
      assert(it->second > 0);
      homecnts.AT(c) = it->second;
      // save the homecnts index in the hash map
      it->second = c++;
    }
    assert(c == (int)homecnts.size());
  }

  timer->stop(npoints);
}

static void interpolate_components(
    cGH const *const cctkGH, int const coord_system_handle,
    int const coord_group, int const minrl, int const maxrl,
    int const maxncomps, bool const want_global_mode,
    int const prolongation_order_time, int const N_dims,
    vector<int> const &homecnts, std::map<int, int> const &homecntsmap,
    vector<int> const &recvcnt, vector<CCTK_REAL *> const &coords,
    vector<char *> const &outputs, CCTK_INT *const per_proc_statuses,
    CCTK_INT *const per_proc_retvals, vector<CCTK_INT> const &operand_indices,
    vector<CCTK_INT> const &time_deriv_order,
    vector<int> const &num_time_derivs, CCTK_INT const local_interp_handle,
    CCTK_INT const param_table_handle, CCTK_REAL const current_time,
    CCTK_REAL const delta_time, int const N_input_arrays,
    int const N_output_arrays, CCTK_INT const output_array_type_codes[],
    CCTK_INT const input_array_variable_indices[]) {
  static Timer timer("CarpetInterp::interpolate_components");
  timer.start();

  // Find out about the number of time levels for interpolation
  // for all input arrays
  vector<CCTK_INT> interp_num_time_levels(N_input_arrays, 0);
  for (int n = 0; n < N_input_arrays; n++) {
    int const vi = input_array_variable_indices[n];
    if (vi >= 0) {
      int const gi = CCTK_GroupIndexFromVarI(vi);
      assert(gi >= 0 and gi < CCTK_NumGroups());
      int const table = CCTK_GroupTagsTableI(gi);
      assert(table >= 0);

      Util_TableGetInt(table, &interp_num_time_levels[n],
                       "InterpNumTimelevels");
    }
  }

#ifndef _NDEBUG
  // Ensure that this processor is only supposed to interpolate
  // points from maps and components that are actually located on
  // this processor
  for (int rl = minrl; rl < maxrl; ++rl) {
    for (int m = 0; m < maps; ++m) {
      gh const *const hh = arrdata.AT(coord_group).AT(m).hh;
      for (int c = 0; c < hh->components(rl); ++c) {
        for (int p = 0; p < dist::size(); ++p) {
          int const cidx = component_idx(p, m, rl, c);
          if (homecntsmap.find(cidx) != homecntsmap.end()) {
            assert(hh->is_local(rl, c));
          }
        }
      }
    }
  }
#endif

  int const tl = 0;
  for (int rl = minrl; rl < maxrl; ++rl) {

    // Number of neccessary time levels
    // CCTK_REAL const level_time = cctkGH->cctk_time;
    CCTK_REAL const level_time = tt->get_time(mglevel, rl, tl);
    vector<int> num_tl(N_input_arrays, 0);
    vector<bool> need_time_interp(N_output_arrays);
    for (int m = 0; m < N_output_arrays; ++m) {
      need_time_interp.AT(m) = num_time_derivs.AT(m) > 0 or
                               (fabs(current_time - level_time) >
                                1e-12 * (fabs(level_time) + fabs(current_time) +
                                         fabs(cctkGH->cctk_delta_time)));
      assert(not(not want_global_mode and num_time_derivs.AT(m) == 0 and
                 need_time_interp.AT(m)));
      int const n = operand_indices.AT(m);
      num_tl.AT(n) =
          max(num_tl.AT(n),
              (need_time_interp.AT(m)
                   ? max(num_time_derivs.AT(m) + 1, prolongation_order_time + 1)
                   : 1));
    }

    // TODO: Loop only over those maps and components that exist for
    // this variable group
    for (int m = 0; m < maps; ++m) {
      for (int lc = 0; lc < vhh.AT(m)->local_components(rl); ++lc) {
        int const c = vhh.AT(m)->get_component(rl, lc);
        for (int p = 0; p < dist::size(); ++p) {
          // short cut if there's nothing to interpolate for processor p
          if (recvcnt[p] <= 0)
            continue;

          int const cidx = component_idx(p, m, rl, c);
          std::map<int, int>::const_iterator it = homecntsmap.find(cidx);
          if (it != homecntsmap.end()) {
            int const idx = it->second;
            assert(idx < (int)homecnts.size());
            interpolate_single_component(
                cctkGH, coord_system_handle, coord_group, N_dims,
                homecnts.AT(idx), coords.AT(idx), outputs.AT(idx),
                per_proc_statuses[p], per_proc_retvals[p], operand_indices,
                time_deriv_order, interp_num_time_levels, local_interp_handle,
                param_table_handle, rl, m, lc, num_tl, need_time_interp,
                current_time, delta_time, N_input_arrays, N_output_arrays,
                output_array_type_codes, input_array_variable_indices);
          }
        } // for p
      }   // for lc
    }     // for m
  }       // for rl

  timer.stop(0);
}

///////////////////////////////////////////////////////////////////////////

/** A list of time values corresoponding to the last few timelevels
 * on the given patch.
 *
 * Values are determined by secret magic.
 * */
class InterpolationTimes : private vector<CCTK_REAL> {
public:
  InterpolationTimes(int const rl, int const num_timelevels_)
      : vector<CCTK_REAL>(num_timelevels_) {
    for (int tl = 0; tl < num_timelevels_; ++tl) {
      at(tl) = tt->get_time(mglevel, rl, tl);
    }
  }

  CCTK_REAL const &operator[](size_type i) const { return at(i); }

  int num_timelevels() const { return (int)size(); }
};

/** Represents interpolation coefficients, or weights, for a given
 *  derivative order and set  of time levels.
 *
 * */
class InterpolationWeights : private vector<CCTK_REAL> {
public:
  InterpolationWeights(CCTK_INT deriv_order, const InterpolationTimes &times,
                       const CCTK_REAL current_time, const CCTK_REAL delta_time)
      : vector<CCTK_REAL>(times.num_timelevels()) {
    CCTK_INT num_timelevels = times.num_timelevels();
    // Initialise weights
    switch (deriv_order) {
    case 0:
      switch (num_timelevels) {
      case 1:
        for_no_interp(times, current_time, delta_time);
        break;
      case 2:
        for_linear_2_pt_interp(times, current_time, delta_time);
        break;
      case 3:
        for_quadratic_3_pt_interp(times, current_time, delta_time);
        break;
      default:
        assert(0);
      }
      break;

    case 1:
      switch (num_timelevels) {
      case 2:
        for_1st_deriv_linear_2_pt_interp_1st_order(times, current_time,
                                                   delta_time);
        break;
      case 3:
        for_1st_deriv_quadratic_3_pt_interp_2nd_order(times, current_time,
                                                      delta_time);
        break;
      default:
        assert(0);
      }
      break;

    case 2:
      switch (num_timelevels) {
      case 3:
        for_2nd_deriv_quadratic_3_pt_interp_2nd_order(times, current_time,
                                                      delta_time);
        break;
      default:
        assert(0);
      }
      break;

    default:
      assert(0);
    } // switch time_deriv_order
  }

  CCTK_REAL const &operator[](size_type i) const { return at(i); }

private:
  void for_no_interp(const InterpolationTimes &t, const CCTK_REAL t0,
                     const CCTK_REAL dt) {
    // We have to assume that any GF with one timelevel is constant in time
    at(0) = 1.0;
  }

  void for_linear_2_pt_interp(const InterpolationTimes &t, const CCTK_REAL t0,
                              const CCTK_REAL dt) {
    at(0) = (t0 - t[1]) / (t[0] - t[1]);
    at(1) = (t0 - t[0]) / (t[1] - t[0]);
  }

  void for_quadratic_3_pt_interp(const InterpolationTimes &t,
                                 const CCTK_REAL t0, const CCTK_REAL dt) {
    at(0) = (t0 - t[1]) * (t0 - t[2]) / ((t[0] - t[1]) * (t[0] - t[2]));
    at(1) = (t0 - t[0]) * (t0 - t[2]) / ((t[1] - t[0]) * (t[1] - t[2]));
    at(2) = (t0 - t[0]) * (t0 - t[1]) / ((t[2] - t[0]) * (t[2] - t[1]));
  }

  void for_1st_deriv_linear_2_pt_interp_1st_order(const InterpolationTimes &t,
                                                  const CCTK_REAL t0,
                                                  const CCTK_REAL dt) {
    at(0) = 1 / (dt * (t[0] - t[1]));
    at(1) = 1 / (dt * (t[1] - t[0]));
  }

  void for_1st_deriv_quadratic_3_pt_interp_2nd_order(
      const InterpolationTimes &t, const CCTK_REAL t0, const CCTK_REAL dt) {
    at(0) = ((t0 - t[2]) + (t0 - t[1])) / (dt * (t[0] - t[1]) * (t[0] - t[2]));
    at(1) = ((t0 - t[2]) + (t0 - t[0])) / (dt * (t[1] - t[0]) * (t[1] - t[2]));
    at(2) = ((t0 - t[1]) + (t0 - t[0])) / (dt * (t[2] - t[0]) * (t[2] - t[1]));
  }

  void for_2nd_deriv_quadratic_3_pt_interp_2nd_order(
      const InterpolationTimes &t, const CCTK_REAL t0, const CCTK_REAL dt) {
    at(0) = 2 / (dt * dt * (t[0] - t[1]) * (t[0] - t[2]));
    at(1) = 2 / (dt * dt * (t[1] - t[0]) * (t[1] - t[2]));
    at(2) = 2 / (dt * dt * (t[2] - t[0]) * (t[2] - t[1]));
  }
};

static void interpolate_single_component(
    cGH const *const cctkGH, int const coord_system_handle,
    int const coord_group, int const N_dims, int const npoints,
    CCTK_REAL const *const coords, char *const outputs,
    CCTK_INT &overall_status, CCTK_INT &overall_retval,
    vector<CCTK_INT> const &operand_indices,
    vector<CCTK_INT> const &time_deriv_order,
    vector<CCTK_INT> const &interp_num_time_levels,
    CCTK_INT const local_interp_handle, CCTK_INT const param_table_handle,
    int const rl, int const m, int const lc, vector<int> const &num_tl,
    vector<bool> const &need_time_interp, CCTK_REAL const current_time,
    CCTK_REAL const delta_time, int const N_input_arrays,
    int const N_output_arrays, CCTK_INT const output_array_type_codes[],
    CCTK_INT const input_array_variable_indices[]) {
  static Timer timer("CarpetInterp::interpolate_single_component");
  timer.start();

  // Find out the datatypes of all input grid arrays
  vector<CCTK_INT> input_array_type_codes(N_input_arrays, -1);
  vector<const void *> input_arrays(N_input_arrays);
  for (int n = 0; n < N_input_arrays; ++n) {
    if (input_array_variable_indices[n] >= 0) {
      input_array_type_codes[n] =
          CCTK_VarTypeI(input_array_variable_indices[n]);
    }
  }

  // Do the processor-local interpolation
  // Find out about the local geometry
  rvect lower, upper, delta;

  // Get global origin and spacing of the underlying coordinate system
  int const grouptype = CCTK_GroupTypeI(coord_group);
  switch (grouptype) {

  case CCTK_GF: {
    jvect gsh;
    GetCoordRange(cctkGH, m, mglevel, dim, &gsh[0], &lower[0], &upper[0],
                  &delta[0]);
    break;
  }

  case CCTK_SCALAR:
  case CCTK_ARRAY: {
#ifdef NEW_COORD_API
    int const iret1 = Util_TableGetRealArray(coord_system_handle, N_dims,
                                             &lower[0], "COMPMIN");
    int const iret2 =
        Util_TableGetRealArray(coord_system_handle, N_dims, &delta[0], "DELTA");
    assert(iret1 == N_dims and iret2 == N_dims);
#else
    char const *const coord_system_name =
        CCTK_CoordSystemName(coord_system_handle);
    assert(CCTK_CoordSystemDim(coord_system_name) >= N_dims);

    for (int d = 0; d < N_dims; ++d) {
      int const iret = CCTK_CoordRange(cctkGH, &lower[d], &upper[d], d + 1,
                                       NULL, coord_system_name);
      assert(iret == 0);
    }

    // int const m = 0;
    assert(m == 0); // We may be looping over too many maps
    // delta for the Carpet grid indices
    ibbox const &baseextent =
        arrdata.AT(coord_group).AT(m).hh->baseextents.AT(mglevel).AT(0);
    delta = (upper - lower) / rvect(baseextent.upper() - baseextent.lower());
#endif
    break;
  }

  default:
    assert(0);
  }

  // Get processor-local origin and spacing
  // cGroupDynamicData coord_group_data;
  // CCTK_GroupDynamicData (cctkGH, coord_group, &coord_group_data);
  // To do: do this via hh->bases instead
  const ibbox &coarseext = vhh.AT(m)->baseextents.AT(mglevel).AT(0);
  const ibbox &baseext = vhh.AT(m)->baseextents.AT(mglevel).AT(rl);
  int const c = vhh.AT(m)->get_component(rl, lc);
  const ibbox &ext = vdd.AT(m)->light_boxes.AT(mglevel).AT(rl).AT(c).exterior;
  // ivect const lsh = ext.shape() / ext.stride();
  for (int d = 0; d < N_dims; ++d) {
    // if (grouptype == CCTK_GF) {
    //   assert (maxspacereflevelfact[d] % cctkGH->cctk_levfac[d] == 0);
    //   delta[d] /= cctkGH->cctk_levfac[d];
    //   lower[d] += (delta[d] *
    //                (cctkGH->cctk_lbnd[d] +
    //                 1.0 * cctkGH->cctk_levoff[d] /
    //                 cctkGH->cctk_levoffdenom[d]));
    // } else {
    //   lower[d] += delta[d] * cctkGH->cctk_lbnd[d];
    // }
    if (grouptype == CCTK_GF) {
      assert(maxspacereflevelfact[d] % spacereffacts.AT(rl)[d] == 0);
      delta[d] /= spacereffacts.AT(rl)[d];
      ivect const lbnd = (ext.lower() - baseext.lower()) / ext.stride();
      ivect const levoff = baseext.lower() - coarseext.lower();
      ivect const levoffdenom = baseext.stride();
      lower[d] += delta[d] * (lbnd[d] + 1.0 * levoff[d] / levoffdenom[d]);
    } else {
      ivect const lbnd = (ext.lower() - baseext.lower()) / ext.stride();
      lower[d] += delta[d] * lbnd[d];
    }
  }

  void const *tmp_coords[dim];
  for (int d = 0; d < N_dims; ++d) {
    tmp_coords[d] = coords + d * npoints;
  }

  int max_num_tl = 0;
  for (int n = 0; n < N_input_arrays; ++n) {
    max_num_tl = max(max_num_tl, num_tl.AT(n));
  }
  vector<vector<void *> > tmp_output_arrays(max_num_tl);

  for (int tl = 0; tl < max_num_tl; ++tl) {

    for (int n = 0; n < N_input_arrays; ++n) {

      int const vi = input_array_variable_indices[n];
      assert(vi == -1 or (vi >= 0 and vi < CCTK_NumVars()));

      if (vi == -1) {
        input_arrays[n] = NULL;
      } else {
        int const interp_num_tl = interp_num_time_levels.AT(n) > 0
                                      ? interp_num_time_levels.AT(n)
                                      : num_tl.AT(n);

        // Do a dummy interpolation from a later timelevel
        // if the desired timelevel does not exist
        int const my_tl = tl >= interp_num_tl ? 0 : tl;
        assert(my_tl < num_tl.AT(n));

        // Are there enough time levels?
        int const gi = CCTK_GroupIndexFromVarI(vi);
        int const active_tl =
            groupdata.AT(gi).activetimelevels.AT(mglevel).AT(rl);
        if (active_tl <= my_tl) {
          char *const fullname = CCTK_FullName(vi);
          CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                     "Grid function \"%s\" has only %d active time levels on "
                     "refinement level %d; this is not enough for time "
                     "interpolation",
                     fullname, active_tl, rl);
          free(fullname);
        }

        int const vi0 = CCTK_FirstVarIndexI(gi);
        ggf const *const ff = arrdata.AT(gi).AT(m).data.AT(vi - vi0);
        void const *const ptr =
            ff->data_pointer(my_tl, rl, lc, mglevel)->storage();
        input_arrays[n] = ptr;
      }
    } // for input arrays

    tmp_output_arrays[tl].resize(N_output_arrays);
    for (int j = 0; j < N_output_arrays; ++j) {
      if (need_time_interp.AT(j)) {
        if (output_array_type_codes[j] != CCTK_VARIABLE_REAL) {
          CCTK_WARN(CCTK_WARN_ABORT,
                    "time interpolation into output arrays of datatype "
                    "other than CCTK_VARIABLE_REAL is not supported");
        }
        tmp_output_arrays[tl][j] = new CCTK_REAL[npoints];
      } else {
        const int vartypesize = CCTK_VarTypeSize(output_array_type_codes[j]);
        tmp_output_arrays[tl][j] = outputs + j * npoints * vartypesize;
      }
    }

    vector<CCTK_INT> per_point_status(npoints);
    int ierr = Util_TableSetPointer(
        param_table_handle, &per_point_status.front(), "per_point_status");
    assert(ierr >= 0);

    vector<CCTK_INT> lsh(N_dims);
    for (int d = 0; d < N_dims; ++d)
      lsh.AT(d) = (ext.shape() / ext.stride())[d];
    const int retval = CCTK_InterpLocalUniform(
        N_dims, local_interp_handle, param_table_handle, &lower[0], &delta[0],
        npoints, CCTK_VARIABLE_REAL, tmp_coords, N_input_arrays, &lsh[0],
        &input_array_type_codes.front(), &input_arrays.front(), N_output_arrays,
        output_array_type_codes, &tmp_output_arrays[tl].front());
    if (retval) {
      CCTK_VWarn(CCTK_WARN_DEBUG, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "The local interpolator returned the error code %d", retval);
    }

    overall_retval = min(overall_retval, (CCTK_INT)retval);
    for (int n = 0; n < (int)per_point_status.size(); n++) {
      overall_status = min(overall_status, per_point_status[n]);
    }
    ierr = Util_TableDeleteKey(param_table_handle, "per_point_status");
    assert(not ierr);

  } // for tl

  // Interpolate in time, if neccessary
  for (int j = 0; j < N_output_arrays; ++j) {
    if (need_time_interp.AT(j)) {

      // Find input array for this output array
      int const n = operand_indices.AT(j);
      CCTK_INT const deriv_order = time_deriv_order.AT(j);

      int const interp_num_tl = interp_num_time_levels.AT(n) > 0
                                    ? interp_num_time_levels.AT(n)
                                    : num_tl.AT(n);
      const InterpolationTimes times(rl, interp_num_tl);
      const InterpolationWeights tfacs(deriv_order, times, current_time,
                                       delta_time);

      // Interpolate
      assert(output_array_type_codes[j] == CCTK_VARIABLE_REAL);
      for (int k = 0; k < npoints; ++k) {
        CCTK_REAL &dest = static_cast<CCTK_REAL *>(
            static_cast<void *>(outputs))[k + j * npoints];
        dest = 0;
        for (int tl = 0; tl < interp_num_tl; ++tl) {
          CCTK_REAL const src =
              ((CCTK_REAL const *)tmp_output_arrays[tl][j])[k];
          dest += tfacs[tl] * src;
        }
      }

      for (int tl = 0; tl < max_num_tl; ++tl) {
        delete[](CCTK_REAL *)tmp_output_arrays[tl][j];
      }

    } // if need_time_interp
  }   // for j

  timer.stop(0);
}

} // namespace CarpetInterp
