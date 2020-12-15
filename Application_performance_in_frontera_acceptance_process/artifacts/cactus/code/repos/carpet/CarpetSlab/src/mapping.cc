#include <cassert>

#include "cctk.h"

#include "util_Table.h"

#include "carpet.hh"

#include "mapping.hh"

namespace CarpetSlab {

using namespace Carpet;

int StoreMapping(mapping *const mp) {
  int const table = Util_TableCreate(UTIL_TABLE_FLAGS_DEFAULT);
  assert(table >= 0);
  int const ierr = Util_TableSetPointer(table, mp, "mapping");
  assert(ierr >= 0);
  return table;
}

mapping *RetrieveMapping(int const table) {
  CCTK_POINTER mp;
  int const ierr = Util_TableGetPointer(table, &mp, "mapping");
  assert(ierr >= 0);
  return (mapping *)mp;
}

void DeleteMapping(int const table) {
  int const ierr = Util_TableDestroy(table);
  assert(ierr >= 0);
}

CCTK_INT
CarpetSlab_LocalMappingByIndex(
    CCTK_POINTER_TO_CONST const cctkGH_, CCTK_INT const vindex,
    CCTK_INT const hdim, CCTK_INT const *const direction,
    CCTK_INT const *const origin, CCTK_INT const *const extent,
    CCTK_INT const *const downsample_, CCTK_INT const table_handle,
    conversion_fn_ptr const conversion_fn, CCTK_INT *const hsize_local,
    CCTK_INT *const hsize_global, CCTK_INT *const hoffset_global) {
  CCTK_WARN(0, "not implemented");
  return 0;
}

CCTK_INT
CarpetSlab_GlobalMappingByIndex(
    CCTK_POINTER_TO_CONST const cctkGH_, CCTK_INT const vindex,
    CCTK_INT const hdim, CCTK_INT const *const direction,
    CCTK_INT const *const origin, CCTK_INT const *const extent,
    CCTK_INT const *const downsample_, CCTK_INT const table_handle,
    conversion_fn_ptr const conversion_fn, CCTK_INT *const hsize) {
  cGH const *const cctkGH = (cGH const *)cctkGH_;

  // Check arguments
  assert(cctkGH);
  assert(vindex >= 0 && vindex < CCTK_NumVars());
  assert(hdim >= 0 && hdim <= dim);
  assert(direction);
  assert(origin);
  assert(extent);
  // assert (downsample);
  // assert (table_handle>=0);
  assert(hsize);

  // Get more information
  int const vdim = CCTK_GroupDimFromVarI(vindex);
  assert(vdim >= 0 && vdim <= dim);
  assert(hdim <= vdim);

  // Not implemented
  assert(!conversion_fn);

  // Allocate memory
  mapping *mp = new mapping;

  // Calculate more convenient representation of the direction
  vector<int> dirs(hdim);
  for (int d = 0; d < hdim; ++d) {
    for (int dd = 0; dd < vdim; ++dd) {
      if (direction[d * vdim + dd] != 0) {
        dirs[d] = dd + 1;
        goto found;
      }
    }
    assert(0);
  found:;
    for (int dd = 0; dd < vdim; ++dd) {
      assert((direction[d * vdim + dd] != 0) == (dirs[d] == dd + 1));
    }
    for (int dd = 0; dd < d; ++dd) {
      assert(dirs[dd] != dirs[d]);
    }
  }

  // Calculate lengths
  vector<CCTK_INT> downsample(hdim);
  for (int dd = 0; dd < hdim; ++dd) {
    downsample[dd] = downsample_ ? downsample_[dd] : 1;
    if (extent[dd] < 0) {
      int const oldmap = Carpet::map;
      int const grouptype = CCTK_GroupTypeFromVarI(vindex);
      assert(grouptype >= 0);
      if (grouptype == CCTK_GF) {
        assert(reflevel >= 0);
        assert(oldmap >= 0 || maps == 1);
        if (oldmap == -1) {
          enter_singlemap_mode(const_cast<cGH *>(cctkGH), 0, grouptype);
        }
      }
      int gsh[dim];
      int ierr = CCTK_GroupgshVI(cctkGH, dim, gsh, vindex);
      assert(!ierr);
      if (grouptype == CCTK_GF) {
        if (oldmap == -1) {
          leave_singlemap_mode(const_cast<cGH *>(cctkGH));
        }
      }
      const int totlen = gsh[dirs[dd] - 1];
      assert(totlen >= 0);
      // Partial argument check
      assert(origin[dirs[dd] - 1] >= 0);
      assert(origin[dirs[dd] - 1] <= totlen);
      assert(downsample[dd] > 0);
      hsize[dd] = (totlen - origin[dirs[dd] - 1]) / downsample[dd];
    } else {
      hsize[dd] = extent[dd];
    }
    assert(hsize[dd] >= 0);
  }

  // Store information
  mp->vindex = vindex;
  mp->hdim = hdim;
  mp->origin.resize(vdim);
  mp->dirs.resize(hdim);
  mp->stride.resize(hdim);
  mp->length.resize(hdim);
  for (size_t d = 0; d < (size_t)vdim; ++d) {
    mp->origin[d] = origin[d];
  }
  for (size_t d = 0; d < (size_t)hdim; ++d) {
    mp->dirs[d] = dirs[d];
    mp->stride[d] = downsample[d];
    mp->length[d] = hsize[d];
  }

  return StoreMapping(mp);
}

CCTK_INT
CarpetSlab_FreeMapping(CCTK_INT const mapping_handle) {
  // Check arguments
  assert(mapping_handle >= 0);

  // Get mapping
  mapping *mp = RetrieveMapping(mapping_handle);
  assert(mp);

  // Delete storage
  DeleteMapping(mapping_handle);

  delete mp;

  return 0;
}

} // namespace CarpetSlab
