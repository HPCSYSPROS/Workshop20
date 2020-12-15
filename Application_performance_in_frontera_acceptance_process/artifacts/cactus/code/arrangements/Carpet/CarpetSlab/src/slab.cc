#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <algorithm>
#include <limits>
#include <vector>

#include "cctk.h"

#include "util_Table.h"

#include "bbox.hh"
#include "bboxset.hh"
#include "dh.hh"
#include "gdata.hh"
#include "gh.hh"
#include "ggf.hh"
#include "vect.hh"

#include "carpet.hh"

#include "mapping.hh"
#include "slab.hh"

namespace CarpetSlab {

using namespace std;
using namespace Carpet;

void FillSlab(const cGH *const cgh, const int dest_proc, const int n,
              const int ti, const int hdim, const int origin[/*vdim*/],
              const int dirs[/*hdim*/], const int stride[/*hdim*/],
              const int length[/*hdim*/], void *const hdata) {
  int ierr;

  // Check Cactus grid hierarchy
  assert(cgh);

  // Check destination processor
  assert(dest_proc >= -1 && dest_proc < CCTK_nProcs(cgh));

  // Check variable index
  assert(n >= 0 && n < CCTK_NumVars());

  // Get info about variable
  const int group = CCTK_GroupIndexFromVarI(n);
  assert(group >= 0);
  const int n0 = CCTK_FirstVarIndexI(group);
  assert(n0 >= 0);
  const int var = n - n0;
  assert(var >= 0);

  // Get info about group
  cGroup gp;
  ierr = CCTK_GroupData(group, &gp);
  assert(!ierr);
  assert(gp.dim <= dim);
  assert(CCTK_QueryGroupStorageI(cgh, group));
  const int typesize = CCTK_VarTypeSize(gp.vartype);
  assert(typesize > 0);

  if (gp.grouptype == CCTK_GF && reflevel == -1) {
    CCTK_WARN(0, "It is not possible to use hyperslabbing for a grid function "
                 "in meta mode or global mode (use singlemap mode instead)");
  }
  const int rl = gp.grouptype == CCTK_GF ? reflevel : 0;
  assert(rl >= 0);

  if (gp.grouptype == CCTK_GF && Carpet::map == -1 && maps > 1) {
    CCTK_WARN(0, "It is not possible to use hyperslabbing for a grid function "
                 "in level mode when there are multiple maps (use singlemap "
                 "mode instead, or make sure that there is only one map)");
  }
  const int m = gp.grouptype == CCTK_GF ? (maps > 1 ? Carpet::map : 0) : 0;
  assert(m >= 0);

  const int oldmap = Carpet::map;
  if (gp.grouptype == CCTK_GF && oldmap == -1) {
    enter_singlemap_mode(const_cast<cGH *>(cgh), m, gp.grouptype);
  }

  // Check dimension
  assert(hdim >= 0 && hdim <= gp.dim);

  // Get more info about group
  cGroupDynamicData gd;
  ierr = CCTK_GroupDynamicData(cgh, group, &gd);
  assert(!ierr);
  const vect<int, dim> sizes = vect<int, dim>::ref(gd.gsh);
  for (int d = 0; d < dim; ++d) {
    assert(sizes[d] >= 0);
  }

  // Check timelevel
  const int num_tl = gp.numtimelevels;
  assert(ti >= 0 && ti < num_tl);
  const int tl = -ti;

  // Check origin
  for (int d = 0; d < dim; ++d) {
    assert(origin[d] >= 0 && origin[d] <= sizes[d]);
  }

  // Check directions
  for (int dd = 0; dd < hdim; ++dd) {
    assert(dirs[dd] >= 1 && dirs[dd] <= dim);
  }

  // Check stride
  for (int dd = 0; dd < hdim; ++dd) {
    assert(stride[dd] > 0);
  }

  // Check length
  for (int dd = 0; dd < hdim; ++dd) {
    assert(length[dd] >= 0);
  }

  // Check extent
  for (int dd = 0; dd < hdim; ++dd) {
    assert(origin[dirs[dd] - 1] + length[dd] <= sizes[dirs[dd] - 1]);
  }

  // Get insider information about variable
  const gh *myhh;
  const dh *mydd;
  const ggf *myff;
  assert(group < (int)arrdata.size());
  myhh = arrdata.at(group).at(m).hh;
  assert(myhh);
  mydd = arrdata.at(group).at(m).dd;
  assert(mydd);
  assert(var < (int)arrdata.at(group).at(m).data.size());
  myff = arrdata.at(group).at(m).data.at(var);
  assert(myff);

  // Detemine collecting processor
  const int collect_proc = dest_proc < 0 ? 0 : dest_proc;

  // Determine own rank
  const int rank = CCTK_MyProc(cgh);

  // Sanity check
  // (if this fails, someone requested an insane number of grid points)
  {
    int max = numeric_limits<int>::max();
    for (int dd = 0; dd < hdim; ++dd) {
      assert(length[dd] >= 0 && length[dd] <= max);
      if (length[dd] > 0)
        max /= length[dd];
    }
    assert(typesize <= max);
  }

  // Calculate global size
  int totalsize = 1;
  for (int dd = 0; dd < hdim; ++dd) {
    totalsize *= length[dd];
  }

  // Allocate memory
  assert(hdata);
  if (dest_proc == -1 || rank == dest_proc) {
    memset(hdata, 0, totalsize * typesize);
  }

  // Get sample data
  const gdata *mydata;
  mydata = myff->data_pointer(tl, rl, 0, 0);

  // Stride of data in memory
  const vect<int, dim> str = mydata->extent().stride();

  // Stride of collected data
  vect<int, dim> hstr = str;
  for (int dd = 0; dd < hdim; ++dd) {
    hstr[dirs[dd] - 1] *= stride[dd];
  }

  // Lower bound of collected data
  vect<int, dim> hlb(0);
  for (int d = 0; d < gp.dim; ++d) {
    hlb[d] = origin[d] * str[d];
  }

  // Upper bound of collected data
  vect<int, dim> hub = hlb;
  for (int dd = 0; dd < hdim; ++dd) {
    hub[dirs[dd] - 1] += (length[dd] - 1) * hstr[dirs[dd] - 1];
  }

  // Calculate extent to collect
  const bbox<int, dim> hextent(hlb, hub, hstr);
  assert(hextent.size() == totalsize);

  // Create collector data object
  void *myhdata = rank == collect_proc ? hdata : 0;
  size_t const mymemsize = totalsize * typesize;
  gdata *const alldata = mydata->make_typed(-1, error_centered, op_sync);
  alldata->allocate(hextent, collect_proc, myhdata, mymemsize);

  // Done with the temporary stuff
  mydata = 0;

  for (comm_state state; !state.done(); state.step()) {

    // Loop over all components, copying data from them
    BEGIN_COMPONENT_LOOP(cgh, gp.grouptype) {

      // Get data object
      mydata = myff->data_pointer(tl, rl, component, mglevel);

      // Calculate overlapping extents
      const bboxset<int, dim> myextents =
          mydd->light_boxes.at(mglevel).at(rl).at(component).interior & hextent;

      // Loop over overlapping extents
      for (bboxset<int, dim>::const_iterator ext_iter = myextents.begin();
           ext_iter != myextents.end(); ++ext_iter) {

        // Copy data
        int const proc = myhh->processor(reflevel, component);
        alldata->copy_from(state, mydata, *ext_iter, *ext_iter, NULL,
                           collect_proc, proc);
      }
    }
    END_COMPONENT_LOOP;

  } // for step

  // Copy result to all processors
  if (dest_proc == -1) {
    vector<gdata *> tmpdata(CCTK_nProcs(cgh));

    for (int proc = 0; proc < CCTK_nProcs(cgh); ++proc) {
      if (proc != collect_proc) {
        void *myhdata = rank == proc ? hdata : 0;
        tmpdata.at(proc) = mydata->make_typed(-1, error_centered, op_sync);
        tmpdata.at(proc)->allocate(alldata->extent(), proc, myhdata);
      }
    }

    for (comm_state state; not state.done(); state.step()) {
      for (int proc = 0; proc < CCTK_nProcs(cgh); ++proc) {
        if (proc != collect_proc) {
          tmpdata.at(proc)->copy_from(state, alldata, alldata->extent(),
                                      alldata->extent(), NULL, proc,
                                      collect_proc);
        }
      }
    }

    for (int proc = 0; proc < CCTK_nProcs(cgh); ++proc) {
      if (proc != collect_proc) {
        delete tmpdata.at(proc);
      }
    }

  } // Copy result

  if (gp.grouptype == CCTK_GF && oldmap == -1) {
    leave_singlemap_mode(const_cast<cGH *>(cgh));
  }

  delete alldata;
}

} // namespace CarpetSlab
