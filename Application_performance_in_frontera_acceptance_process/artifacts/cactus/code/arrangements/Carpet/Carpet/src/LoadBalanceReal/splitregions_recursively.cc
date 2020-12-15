#include <cctk.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <list>
#include <sstream>
#include <string>
#include <vector>

#include <sys/stat.h>
#include <sys/types.h>

#ifdef CCTK_MPI
#include <mpi.h>
#else
#include "nompi.h"
#endif

#include <Timer.hh>

#include <loopcontrol.h>

#include <bbox.hh>
#include <bboxset.hh>
#include <defs.hh>
#include <dh.hh>
#include <gh.hh>
#include <region.hh>
#include <vect.hh>

#include <carpet.hh>
#include <modes.hh>
#include <variables.hh>

namespace Carpet {

using namespace std;

// Helper routines for spliting regions automatically

// The cost for a region, assuming a cost of 1 per interior point
static rvect cost(region_t const &reg) {
  DECLARE_CCTK_PARAMETERS;
  static rvect costfactor;
  static bool initialised = false;
  if (not initialised) {
    costfactor = rvect(1.0);
    if (dim > 0)
      costfactor[0] = 1.0 / aspect_ratio_x;
    if (dim > 1)
      costfactor[1] = 1.0 / aspect_ratio_y;
    if (dim > 2)
      costfactor[2] = 1.0 / aspect_ratio_z;
  }
  if (reg.extent.empty())
    return rvect(0);
  return rvect(reg.extent.shape() / reg.extent.stride()) * costfactor;
}

struct f_range {
  int lower, upper, stride;
};

struct f_bbox {
  f_range dim[3];

  f_bbox() {}
  f_bbox(ibbox const &box) {
    assert(::dim == 3);
    for (int d = 0; d < 3; ++d) {
      dim[d].lower = box.lower()[d];
      dim[d].upper = box.upper()[d];
      dim[d].stride = box.stride()[d];
    }
  }
  /*explicit*/ operator ibbox() const {
    ivect lower, upper, stride;
    assert(::dim == 3);
    for (int d = 0; d < 3; ++d) {
      lower[d] = dim[d].lower;
      upper[d] = dim[d].upper;
      stride[d] = dim[d].stride;
    }
    return ibbox(lower, upper, stride);
  }
};

struct f_boundary {
  int obound[2][3];

  f_boundary() {}
  f_boundary(b2vect const &ob) {
    assert(::dim == 3);
    for (int d = 0; d < 3; ++d) {
      for (int f = 0; f < 2; ++f) {
        obound[f][d] = ob[f][d];
      }
    }
  }
  /*explicit*/ operator b2vect() const {
    b2vect ob;
    assert(::dim == 3);
    for (int d = 0; d < 3; ++d) {
      for (int f = 0; f < 2; ++f) {
        ob[f][d] = obound[f][d];
      }
    }
    return ob;
  }
};

struct f_superregion2slim {
  f_bbox extent;
  f_boundary outer_boundaries;
  int map;
  int processor;

  /*explicit*/ operator region_t() const {
    region_t reg;
    reg.extent = ibbox(extent);
    reg.outer_boundaries = b2vect(outer_boundaries);
    reg.map = map;
    reg.processor = processor;
    reg.processors = NULL;
    return reg;
  }

  /*explicit*/ operator pseudoregion_t() const {
    pseudoregion_t preg;
    preg.extent = ibbox(extent);
    preg.component = processor;
    return preg;
  }
};

extern "C" CCTK_FCALL void CCTK_FNAME(splitregions_recursively)(
    CCTK_POINTER const &cxx_superregs, int const &nsuperregs,
    CCTK_POINTER const &cxx_regs, int const &nprocs, int const &ghostsize,
    CCTK_REAL const &alpha, int const &limit_size, int const &procid);

void SplitRegionsMaps_Recursively(cGH const *const cctkGH,
                                  vector<vector<region_t> > &superregss,
                                  vector<vector<region_t> > &regss) {
  DECLARE_CCTK_PARAMETERS;

  if (recompose_verbose)
    cout << "SRMR enter" << endl;

  int const nmaps = superregss.size();
  assert(int(regss.size()) == nmaps);
  int map_offset = 1000000000;
  int max_map = 0;
  vector<bool> have_map(maps, false);
  for (int m = 0; m < nmaps; ++m) {
    for (int r = 0; r < int(superregss.AT(m).size()); ++r) {
      map_offset = min(map_offset, superregss.AT(m).AT(r).map);
      max_map = max(max_map, superregss.AT(m).AT(r).map);
      have_map.AT(superregss.AT(m).AT(r).map) = true;
    }
  }
  // Apparently this is not always the case:
  // assert (max_map - map_offset == nmaps - 1);

  int nsuperregs = 0;
  for (int m = 0; m < nmaps; ++m) {
    nsuperregs += superregss.AT(m).size();
  }
  if (recompose_verbose)
    cout << "SRMR nsuperregs " << nsuperregs << endl;

  // Something to do?
  if (nsuperregs == 0)
    return;

  // Collect slices
  vector<region_t> superregs;
  {
    for (int m = 0; m < nmaps; ++m) {
      combine_regions(superregss.AT(m), superregs);
    }
    nsuperregs = superregs.size();

    // If the last region was removed, add a new empty region again.
    // A set of regions (corresponding to a refinement level or a
    // grid array) cannot be empty.
    if (nsuperregs == 0) {
      assert(nmaps == 1); // we should only be here for grid
                          // arrays
      region_t reg;
      reg.extent = ibbox(ivect(0), ivect(-1), ivect(1));
      reg.outer_boundaries = b2vect(bvect(true), bvect(true));
      reg.map = 0;
      superregs.push_back(reg);
      nsuperregs = superregs.size();
    }
  }

  // Create a mapping from this list of regions to maps, and set the
  // map to the superregion index instead
  vector<int> superreg_maps(nsuperregs);
  for (int r = 0; r < nsuperregs; ++r) {
    superreg_maps.AT(r) = superregs.AT(r).map;
    superregs.AT(r).map = r;
  }

  // Handle granularity: artificially reduce the size of the
  // superregions
  for (int r = 0; r < nsuperregs; ++r) {
    region_t &superreg = superregs.AT(r);
    ibbox &ext = superregs.AT(r).extent;
    ivect lo = ext.lower();
    ivect hi = ext.upper() + ext.stride();
    ivect str = ext.stride();
    // cut off outer boundaries
    lo +=
        (ivect(superreg.outer_boundaries[0]) * int(granularity_boundary) * str);
    hi -=
        (ivect(superreg.outer_boundaries[1]) * int(granularity_boundary) * str);
    // multiply stride by granularity
    str *= granularity;
    assert(all((hi - lo) % str == 0));
    ext = ibbox(lo, hi - str, str);
  }

  int const real_nprocs = CCTK_nProcs(cctkGH);
  if (recompose_verbose)
    cout << "SRMR real_nprocs " << real_nprocs << endl;

  // Deactivate some processors if there are too many
  int nprocs;
  if (min_points_per_proc == 0) {
    nprocs = real_nprocs;
  } else {
    CCTK_REAL mycost = 0;
    for (int r = 0; r < nsuperregs; ++r) {
      mycost += prod(cost(superregs.AT(r)));
    }
    int const goodnprocs = int(floor(mycost / min_points_per_proc));
    nprocs = max(1, min(real_nprocs, goodnprocs));
  }
  if (recompose_verbose)
    cout << "SRMR nprocs " << nprocs << endl;

  // Distribute load
  CCTK_POINTER const cxx_superregs = &superregs;
  assert((int)superregs.size() == nsuperregs);
  vector<region_t> regs;
  // regs.reserve (...);
  CCTK_POINTER const cxx_regs = &regs;
  int const ghostsize = vdd.AT(0)->ghost_widths.AT(0)[0][0];
  for (int m = 0; m < maps; ++m) {
    for (int rl = 0; rl < reflevels; ++rl) {
      for (int f = 0; f < 2; ++f) {
        for (int d = 0; d < dim; ++d) {
          assert(vdd.AT(m)->ghost_widths.AT(rl)[f][d] == ghostsize);
        }
      }
    }
  }
  CCTK_REAL const alpha = ghost_zone_cost;
  int const limit_size = true;
  int const procid = CCTK_MyProc(cctkGH);
  CCTK_FNAME(splitregions_recursively)
  (cxx_superregs, nsuperregs, cxx_regs, nprocs, ghostsize, alpha, limit_size,
   procid);
  int nregs = regs.size();

  // Handle granularity: increase the size of the superregions and
  // regions again
  for (int r = 0; r < nsuperregs; ++r) {
    region_t &superreg = superregs.AT(r);
    ibbox &ext = superregs.AT(r).extent;
    ivect lo = ext.lower();
    ivect hi = ext.upper() + ext.stride();
    ivect str = ext.stride();
    // divide stride by granularity
    assert(all(str % int(granularity) == 0));
    str /= granularity;
    // add outer boundaries again
    lo -=
        (ivect(superreg.outer_boundaries[0]) * int(granularity_boundary) * str);
    hi +=
        (ivect(superreg.outer_boundaries[1]) * int(granularity_boundary) * str);
    ext = ibbox(lo, hi - str, str);
  }
  for (int r = 0; r < nregs; ++r) {
    region_t &reg = regs.AT(r);
    ibbox &ext = regs.AT(r).extent;
    ivect lo = ext.lower();
    ivect hi = ext.upper() + ext.stride();
    ivect str = ext.stride();
    // divide stride by granularity
    assert(all(str % int(granularity) == 0));
    str /= int(granularity);
    // add outer boundaries again
    lo -= ivect(reg.outer_boundaries[0]) * int(granularity_boundary) * str;
    hi += ivect(reg.outer_boundaries[1]) * int(granularity_boundary) * str;
    ext = ibbox(lo, hi - str, str);
  }

  if (same_number_of_components_on_each_process) {
    // Ensure all processes have the same number of components
    vector<int> ncomps(nprocs, 0);
    for (int r = 0; r < nregs; ++r) {
      int const p = regs.AT(r).processor;
      assert(p >= 0);
      ++ncomps.AT(p);
    }
    int maxncomps = 0;
    int sumncomps = 0;
    for (int p = 0; p < nprocs; ++p) {
      maxncomps = max(maxncomps, ncomps.AT(p));
      sumncomps += ncomps.AT(p);
    }
    int const missingcomps = maxncomps * nprocs - sumncomps;
    if (missingcomps > 0) {
      // Invent a dummy component
      ibbox const &ext = superregss.AT(0).AT(0).extent;
      region_t dummy;
      dummy.extent =
          ibbox(ext.lower(), ext.lower() - ext.stride(), ext.stride());
      assert(dummy.extent.empty());
      dummy.outer_boundaries = b2vect(true);
      dummy.map = nmaps - 1; // arbitrary choice
      // Insert dummy regions at the end
      regs.resize(nregs + missingcomps, dummy);
      for (int p = 0; p < nprocs; ++p) {
        for (int i = ncomps.AT(p); i < maxncomps; ++i) {
          regs.AT(nregs++).processor = p;
        }
      }
      assert(nregs == int(regs.size()));
      // Insert a superregion
      // TODO: Do we need this? Should we skip this?
      pseudoregion_t sample;
      sample.extent = dummy.extent;
      // sample.component remains unset (component will be set later)
      if (missingcomps == 1) {
        dummy.processors = new ipfulltree(sample);
      } else {
        int const dir = 0; // arbitrary choice
        vector<int> bounds(missingcomps + 1);
        for (int n = 0; n < missingcomps + 1; ++n) {
          bounds.AT(n) = dummy.extent.lower()[dir];
        }
        vector<ipfulltree *> subtrees(missingcomps);
        for (int n = 0; n < missingcomps; ++n) {
          subtrees.AT(n) = new ipfulltree(sample);
        }
        dummy.processors = new ipfulltree(dir, bounds, subtrees);
      }
      superregs.push_back(dummy);
    }
  }

  // Allocate regions, saving the old regions for debugging or
  // self-checking
  vector<vector<region_t> > old_superregss;
  swap(superregss, old_superregss);
  superregss.resize(old_superregss.size());
  assert((int)regss.size() == nmaps);
  for (int m = 0; m < nmaps; ++m) {
    assert(regss.AT(m).empty());
    // regss.AT(m).reserve (...);
    // superregss.AT(m).clear();
    assert(superregss.AT(m).empty());
    // superregss.AT(m).reserve (...);
  }
  // Assign regions
  for (int r = 0; r < nsuperregs; ++r) {
    superregs.AT(r).map = superreg_maps.AT(r); // correct map
    int const m = superregs.AT(r).map - map_offset;
    assert(m >= 0 and m < nmaps);
    superregss.AT(m).push_back(superregs.AT(r));
  }
  // Renumber components
  for (int m = 0; m < nmaps; ++m) {
    int c = 0;
    for (size_t r = 0; r < superregss.AT(m).size(); ++r) {
      ipfulltree *const procs = superregss.AT(m).AT(r).processors;
      assert(procs);
      for (ipfulltree::iterator it(*procs); not it.done(); ++it) {
        (*it).payload().component = c++;
      }
    }
  }
  for (int r = 0; r < nregs; ++r) {
    int const s = regs.AT(r).map;
    assert(s >= 0 and s < nsuperregs);
    regs.AT(r).map = superreg_maps.AT(s); // correct map
    int const m = regs.AT(r).map - map_offset;
    assert(m >= 0 and m < nmaps);
    regss.AT(m).push_back(regs.AT(r));
  }
  // Output regions
  if (recompose_verbose) {
    region_t::full_output = true;
    cout << "SRMR superregss " << superregss << endl;
    cout << "SRMR regss " << regss << endl;
    region_t::full_output = false;
  }

  // Consistency check
  bool has_error = false;
  vector<ibset> all_old_superregss(maps);
  for (size_t m = 0; m < superregss.size(); ++m) {
    for (size_t r = 0; r < old_superregss.AT(m).size(); ++r) {
      region_t const &reg = old_superregss.AT(m).AT(r);
      if (all_old_superregss.AT(reg.map).intersects(reg.extent)) {
        has_error = true;
        cout << "SRMR: old_superregss:\n"
             << "m=" << m << " r=" << r << " reg=" << reg << "\n";
      }
      all_old_superregss.AT(reg.map) += reg.extent;
    }
  }
  vector<ibset> all_superregss(maps);
  for (size_t m = 0; m < superregss.size(); ++m) {
    for (size_t r = 0; r < superregss.AT(m).size(); ++r) {
      region_t const &reg = superregss.AT(m).AT(r);
      if (all_superregss.AT(reg.map).intersects(reg.extent)) {
        has_error = true;
        cout << "SRMR: all_superregss:\n"
             << "m=" << m << " r=" << r << " reg=" << reg << "\n";
      }
      all_superregss.AT(reg.map) += reg.extent;
    }
  }
  for (int m = 0; m < maps; ++m) {
    if (not(all_superregss.AT(m) == all_old_superregss.AT(m))) {
      has_error = true;
      cout << "SRMR: all_superregss m=" << m << "\n";
    }
  }
  vector<ibset> all_regss(maps);
  for (size_t m = 0; m < regss.size(); ++m) {
    for (size_t r = 0; r < regss.AT(m).size(); ++r) {
      region_t const &reg = regss.AT(m).AT(r);
      if (all_regss.AT(reg.map).intersects(reg.extent)) {
        has_error = true;
        cout << "SRMR: all_regss:\n"
             << "m=" << m << " r=" << r << " reg=" << reg << "\n";
      }
      all_regss.AT(reg.map) += reg.extent;
    }
  }
  for (int m = 0; m < maps; ++m) {
    if (not(all_regss.AT(m) == all_old_superregss.AT(m))) {
      has_error = true;
      cout << "SRMR: all_regss m=" << m << "\n";
    }
  }
  if (has_error) {
    region_t::full_output = true;
    cout << "SRMR: all_old_superregss=" << all_old_superregss << "\n"
         << "SRMR: all_superregss=" << all_superregss << "\n"
         << "SRMR: all_regss=" << all_regss << "\n";
    region_t::full_output = false;
    CCTK_WARN(CCTK_WARN_ABORT, "Internal error");
  }

  if (recompose_verbose)
    cout << "SRMR exit" << endl;
}

extern "C" CCTK_FCALL void
    CCTK_FNAME(carpet_get_region)(CCTK_POINTER &cxx_superregs, int const &i,
                                  CCTK_POINTER &cxx_superreg) {
  vector<region_t> &superregs = *static_cast<vector<region_t> *>(cxx_superregs);
  region_t &superreg = superregs.AT(i);
  cxx_superreg = &superreg;
}

extern "C" CCTK_FCALL void
    CCTK_FNAME(carpet_get_bbox)(CCTK_POINTER &cxx_superreg, f_bbox &box,
                                f_boundary &obound) {
  region_t &superreg = *static_cast<region_t *>(cxx_superreg);
  box = f_bbox(superreg.extent);
  obound = f_boundary(superreg.outer_boundaries);
}

extern "C" CCTK_FCALL void
    CCTK_FNAME(carpet_insert_region)(CCTK_POINTER &cxx_regs,
                                     f_superregion2slim const &reg) {
  vector<region_t> &regs = *static_cast<vector<region_t> *>(cxx_regs);
  regs.push_back(region_t(reg));
}

extern "C" CCTK_FCALL void CCTK_FNAME(carpet_create_tree_branch)(
    int const &nch, int const &dir, int const fbounds[],
    CCTK_POINTER cxx_subtrees[], CCTK_POINTER &cxx_tree) {
  vector<int> bounds(nch + 1);
  vector<ipfulltree *> subtrees(nch);
  for (int i = 0; i < nch + 1; ++i) {
    bounds.AT(i) = fbounds[i];
  }
  for (int i = 0; i < nch; ++i) {
    ipfulltree *const tree = static_cast<ipfulltree *>(cxx_subtrees[i]);
    assert(tree->invariant());
    subtrees.AT(i) = tree;
  }
  cxx_tree = new ipfulltree(dir, bounds, subtrees);
}

extern "C" CCTK_FCALL void
    CCTK_FNAME(carpet_create_tree_leaf)(f_superregion2slim const &sreg,
                                        CCTK_POINTER &cxx_tree) {
  cxx_tree = new ipfulltree(pseudoregion_t(sreg));
}

extern "C" CCTK_FCALL void
    CCTK_FNAME(carpet_set_tree)(CCTK_POINTER &cxx_superreg,
                                CCTK_POINTER &cxx_tree) {
  region_t &superreg = *static_cast<region_t *>(cxx_superreg);
  ipfulltree *tree = static_cast<ipfulltree *>(cxx_tree);
  assert(not superreg.processors);
  superreg.processors = tree;
}

} // namespace Carpet
