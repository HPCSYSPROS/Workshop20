#include <cassert>
#include <cstdlib>
#include <iostream>

#include "bboxset.hh"
#include "defs.hh"
#include "dist.hh"
#include "mpi_string.hh"
#include "region.hh"

using namespace std;

region_t::region_t() : processor(-1), processors(NULL) { assert(invariant()); }

region_t::region_t(region_t const &a) {
  assert(a.invariant());
  extent = a.extent;
  outer_boundaries = a.outer_boundaries;
  map = a.map;
  processor = a.processor;
  if (a.processors == NULL) {
    processors = NULL;
  } else {
    processors = new ipfulltree(*a.processors);
  }
  assert(invariant());
}

region_t &region_t::operator=(region_t const &a) {
  if (&a == this)
    return *this; // nothing to do
  assert(invariant());
  if (processors != NULL) {
    delete processors;
  }
  assert(a.invariant());
  extent = a.extent;
  outer_boundaries = a.outer_boundaries;
  map = a.map;
  processor = a.processor;
  if (a.processors == NULL) {
    processors = NULL;
  } else {
    processors = new ipfulltree(*a.processors);
  }
  assert(invariant());
  return *this;
}

region_t::~region_t() {
  assert(invariant());
  if (processors != NULL) {
    delete processors;
  }
}

bool region_t::invariant() const {
  if (processor >= 0 and processors != NULL)
    return false;
  return true;
}

// Compare two regions for equality.
bool operator==(region_t const &a, region_t const &b) {
  return a.extent == b.extent and
         all(all(a.outer_boundaries == b.outer_boundaries)) and
         a.map == b.map and a.processor == b.processor and
         ((a.processors == NULL and b.processors == NULL) or
          (a.processors != NULL and b.processors != NULL and
           *a.processors == *b.processors));
}

// Assign a load to a region
CCTK_REAL
region_t::load() const { return extent.size(); }

// Split a region into two
region_t region_t::split(CCTK_REAL const ratio_new_over_old) {
  assert(ratio_new_over_old >= 0 and ratio_new_over_old <= 1);
  if (extent.empty()) {
    // Don't do anything for empty regions
    return *this;
  }
  // Choose a direction (prefer the z direction)
  int const idir = maxloc1(extent.shape());
  int const np = extent.shape()[idir];
  // Keep the lower part, and split off the upper part
  int const new_np = floor(np * ratio_new_over_old + 0.5);
  int const keep_np = np - new_np;
  // Calculate new region extents
  ivect const lo = extent.lower();
  ivect const up = extent.upper();
  ivect const str = extent.stride();
  ivect const locut = lo + str * ivect::dir(idir) * keep_np;
  ivect const upcut = up - str * ivect::dir(idir) * new_np;

  // Copy the region
  region_t newreg = *this;
  // Set new extents
  extent = ibbox(lo, upcut, str);
  newreg.extent = ibbox(locut, up, str);

  // Mark cutting boundary as not outer boundary
  outer_boundaries[idir][1] = false;
  newreg.outer_boundaries[idir][0] = false;

  return newreg;
}

// Combine a collection of regions. Regions can be combined if they
// abutt on boundaries which are not outer boundaries, ignoring the
// process distribution. This should lead to a canonical
// representations of collections of regions.
//
// We use vectors to represent the collection, but we could also use
// other containers.  oldregs is read, newregs is added-to.  newregs
// is not cleared.
void combine_regions(vector<region_t> const &oldregs,
                     vector<region_t> &newregs) {
  // Find the union of all regions' bounding boxes, and the union of
  // all regions' outer boundaries.  Represent the boundaries as the
  // outermost layer of grid points of the corresponding bounding
  // boxes.
  int const m = oldregs.empty() ? -1 : oldregs.AT(0).map;
  ibset comps;
  ibset cobnds[2][dim];
  for (vector<region_t>::const_iterator ri = oldregs.begin();
       ri != oldregs.end(); ++ri) {
    region_t const &reg = *ri;
    assert(reg.map == m);
    comps += reg.extent;
    for (int f = 0; f < 2; ++f) {
      for (int d = 0; d < dim; ++d) {
        if (reg.outer_boundaries[f][d]) {
          ibbox bnd = reg.extent;
          ivect lo = bnd.lower();
          ivect up = bnd.upper();
          if (f == 0) {
            up[d] = lo[d];
          } else {
            lo[d] = up[d];
          }
          bnd = ibbox(lo, up, bnd.stride());
          cobnds[f][d] += bnd;
        }
      }
    }
  }

  // Reserve (generous) memory for the result
  size_t const needsize = newregs.size() + comps.setsize();
  if (newregs.capacity() < needsize) {
    newregs.reserve(1000 + 2 * needsize);
  }

  // Insert the regions
  for (ibset::const_iterator ci = comps.begin(); ci != comps.end(); ++ci) {
    ibbox const &c = *ci;
    b2vect obnds;
    for (int f = 0; f < 2; ++f) {
      for (int d = 0; d < dim; ++d) {
        obnds[f][d] = cobnds[f][d].intersects(c);
        if (obnds[f][d]) {
          ivect lo = c.lower();
          ivect up = c.upper();
          if (f)
            lo[d] = up[d];
          else
            up[d] = lo[d];
          ibbox const cbnds(lo, up, c.stride());
          if (not((cobnds[f][d] & ibset(c)) == ibset(cbnds))) {
            cout << "cobnds[f][d] = " << cobnds[f][d] << endl
                 << "ibset(c) = " << ibset(c) << endl
                 << "(cobnds[f][d] & ibset(c)) = " << (cobnds[f][d] & ibset(c))
                 << endl
                 << "ibset(cbnds) = " << ibset(cbnds) << endl;
          }
          assert((cobnds[f][d] & ibset(c)) == ibset(cbnds));
        }
      }
    }

    region_t reg;
    reg.extent = c;
    reg.outer_boundaries = obnds;
    reg.map = m;
    reg.processor = -1;
    reg.processors = NULL;
    newregs.push_back(reg);
  }
}

size_t memoryof(region_t const &reg) {
  return memoryof(reg.extent) + memoryof(reg.outer_boundaries) +
         memoryof(reg.map) + memoryof(reg.processor) +
         memoryof(reg.processors) +
         (reg.processors != NULL ? memoryof(*reg.processors) : 0);
}

istream &operator>>(istream &is, region_t &reg) {
  skipws(is);
  consume(is, "region_t");
  skipws(is);
  consume(is, '(');

  skipws(is);
  consume(is, "extent");
  skipws(is);
  consume(is, '=');
  is >> reg.extent;
  skipws(is);
  consume(is, ',');

  skipws(is);
  consume(is, "outer_boundaries");
  skipws(is);
  consume(is, '=');
  is >> reg.outer_boundaries;
  skipws(is);
  consume(is, ',');

  skipws(is);
  consume(is, "map");
  skipws(is);
  consume(is, '=');
  is >> reg.map;
  skipws(is);
  consume(is, ',');

  skipws(is);
  consume(is, "processor");
  skipws(is);
  consume(is, '=');
  is >> reg.processor;
  skipws(is);
  consume(is, ')');

  reg.processors = NULL;

  return is;
}

bool region_t::check_region(bool const is_superregion) const {
  assert(invariant());
  // TODO: empty regions need to be allowed (e.g. for empty grid
  // arrays, or when the domain decomposition assigns no points to a
  // particular process)
  // if (extent.empty()) {
  //   CCTK_WARN(CCTK_WARN_PICKY, "extent is empty");
  //   return false;
  // }
  assert(map >= 0);

  if (is_superregion) {
    // Checking a superregion: a tree without process assignments

    if (processor != -1) {
      CCTK_WARN(CCTK_WARN_PICKY, "process number is defined");
      return false;
    }
    if (not processors) {
      CCTK_WARN(CCTK_WARN_PICKY, "process tree not defined");
      return false;
    }

    // TODO: check outer_boundaries as well
    ibset child_extents;
    check_children(*processors, map, 0, child_extents);
    if (child_extents != extent) {
      CCTK_WARN(CCTK_WARN_PICKY, "child extents not equal to total extent");
      return false;
    }

  } else {
    // Checking a regular region: no tree structure, but has a process
    // assignment

    if (processor < 0 or processor >= dist::size()) {
      CCTK_WARN(CCTK_WARN_PICKY, "process number not defined");
      return false;
    }
    if (processors) {
      CCTK_WARN(CCTK_WARN_PICKY, "process tree is defined");
      return false;
    }
  }

  return true;
}

bool region_t::check_children(ipfulltree const &tree, int const parent_map,
                              int const level, ibset &child_extents) const {
  if (tree.empty()) {
    CCTK_VWarn(CCTK_WARN_PICKY, __LINE__, __FILE__, CCTK_THORNSTRING,
               "level %d: tree is empty", level);
    return false;
  } else if (tree.is_branch()) {
    for (ipfulltree::const_iterator i(tree); not i.done(); ++i) {
      bool const bres =
          check_children(*i, parent_map, level + 1, child_extents);
      if (not bres) {
        CCTK_VWarn(CCTK_WARN_PICKY, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "level %d: error in children", level);
        return false;
      }
    }
  } else if (tree.is_leaf()) {
    ibbox const &ext = tree.payload().extent;
    if (ext.empty()) {
      CCTK_VWarn(CCTK_WARN_PICKY, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "level %d: tree leaf extent is empty", level);
      return false;
    }
    if (map != parent_map) {
      CCTK_VWarn(CCTK_WARN_PICKY, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "level %d: tree leaf map=%d, differs from parent map=%d",
                 level, map, parent_map);
      return false;
    }
    if (child_extents.intersects(ext)) {
      CCTK_VWarn(CCTK_WARN_PICKY, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "level %d: tree leaf extent overlaps other tree extents",
                 level);
      return false;
    }
    child_extents += ext;
  } else {
    assert(0);
  }
  return true;
}

bool region_t::full_output = false;

ostream &operator<<(ostream &os, region_t const &reg) {
  os << "region_t("
     << "extent=" << reg.extent << ","
     << "outer_boundaries=" << reg.outer_boundaries << ","
     << "map=" << reg.map << ","
     << "processor=" << reg.processor;
  if (region_t::full_output) {
    os << ","
       << "processors=";
    if (reg.processors) {
      os << *reg.processors;
    } else {
      os << "NULL";
    }
  }
  os << ")";
  return os;
}

// Create an MPI datatype for a pseudoretion
MPI_Datatype mpi_datatype(pseudoregion_t const &) {
  static bool initialised = false;
  static MPI_Datatype newtype;
  if (not initialised) {
    static pseudoregion_t s;
#define ENTRY(type, name)                                                      \
  {                                                                            \
    sizeof s.name / sizeof(type),         /* count elements */                 \
        (char *) & s.name - (char *) & s, /* offsetof doesn't work (why?) */   \
        dist::mpi_datatype<type>(),       /* find MPI datatype */              \
        STRINGIFY(name),                  /* field name */                     \
        STRINGIFY(type),                  /* type name */                      \
  }
    dist::mpi_struct_descr_t const descr[] = {
        ENTRY(int, extent),
        ENTRY(int, component),
        {1, sizeof s, MPI_UB, "MPI_UB", "MPI_UB"}};
#undef ENTRY
    newtype = dist::create_mpi_datatype(sizeof descr / sizeof descr[0], descr,
                                        "pseudoregion_t", sizeof s);
    initialised = true;
  }
  return newtype;
}

MPI_Datatype mpi_datatype(sendrecv_pseudoregion_t const &) {
  static bool initialised = false;
  static MPI_Datatype newtype;
  if (not initialised) {
    static sendrecv_pseudoregion_t s;
#define ENTRY(type, name)                                                      \
  {                                                                            \
    sizeof s.name / sizeof(type),         /* count elements */                 \
        (char *) & s.name - (char *) & s, /* offsetof doesn't work (why?) */   \
        dist::mpi_datatype<type>(),       /* find MPI datatype */              \
        STRINGIFY(name),                  /* field name */                     \
        STRINGIFY(type),                  /* type name */                      \
  }
    dist::mpi_struct_descr_t const descr[] = {
        ENTRY(pseudoregion_t, send),
        ENTRY(pseudoregion_t, recv),
        {1, sizeof s, MPI_UB, "MPI_UB", "MPI_UB"}};
#undef ENTRY
    newtype = dist::create_mpi_datatype(sizeof descr / sizeof descr[0], descr,
                                        "sendrecv_pseudoregion_t", sizeof s);
    initialised = true;
  }
  return newtype;
}

// Compare two pseudoregions for equality.
bool operator==(pseudoregion_t const &a, pseudoregion_t const &b) {
  return a.extent == b.extent and a.component == b.component;
}

istream &operator>>(istream &is, pseudoregion_t &p) {
  try {
    skipws(is);
    consume(is, "(ext:");
    is >> p.extent;
    skipws(is);
    consume(is, ",c:");
    is >> p.component;
    skipws(is);
    consume(is, ")");
  } catch (input_error &err) {
    cout << "Input error while reading a pseudoregion_t" << endl;
    throw err;
  }
  return is;
}

istream &operator>>(istream &is, sendrecv_pseudoregion_t &srp) {
  try {
    skipws(is);
    consume(is, "(send:");
    is >> srp.send;
    consume(is, ",recv:");
    is >> srp.recv;
    consume(is, ")");
  } catch (input_error &err) {
    cout << "Input error while reading a sendrecv_pseudoregion_t" << endl;
    throw err;
  }
  return is;
}

ostream &operator<<(ostream &os, pseudoregion_t const &p) {
  return os << "(ext:" << p.extent << ",c:" << p.component << ")";
}

ostream &operator<<(ostream &os, sendrecv_pseudoregion_t const &srp) {
  return os << "(send:" << srp.send << ",recv:" << srp.recv << ")";
}

namespace CarpetLib {

template vector<sendrecv_pseudoregion_t>
alltoallv1(MPI_Comm comm, vector<vector<sendrecv_pseudoregion_t> > const &data);
}
