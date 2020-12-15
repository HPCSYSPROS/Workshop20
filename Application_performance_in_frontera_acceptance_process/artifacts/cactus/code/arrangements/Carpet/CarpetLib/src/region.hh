#ifndef REGION_HH
#define REGION_HH

#include <iostream>
#include <vector>

#include "defs.hh"
#include "dist.hh"
#include "bbox.hh"
#include "bboxset.hh"
#include "fulltree.hh"
#include "vect.hh"

// Region description
struct region_t {
  ibbox extent;            // extent
  b2vect outer_boundaries; // outer boundaries
  int map;                 // map to which this region belongs
  int processor;           // process number
  ipfulltree *processors;  // process decomposition

  region_t();
  region_t(region_t const &a);
  region_t &operator=(region_t const &a);
  ~region_t();

  bool invariant() const CCTK_MEMBER_ATTRIBUTE_PURE;

  // For regridding
  CCTK_REAL load() const CCTK_MEMBER_ATTRIBUTE_PURE;
  region_t split(CCTK_REAL ratio_new_over_old);

  // Check whether a region is defined consistently
  bool check_region(bool is_superregion) const;

private:
  bool check_children(ipfulltree const &tree, int parent_map, int level,
                      ibset &child_extents) const;

public:
  // Output process decomposition? (Off by default.)
  static bool full_output;
};

bool operator==(region_t const &a, region_t const &b) CCTK_ATTRIBUTE_PURE;
inline bool operator!=(region_t const &a, region_t const &b) {
  return not(a == b);
}

void combine_regions(vector<region_t> const &oldregs,
                     vector<region_t> &newregs);

size_t memoryof(region_t const &reg) CCTK_ATTRIBUTE_PURE;

istream &operator>>(istream &is, region_t &reg);
ostream &operator<<(ostream &os, region_t const &reg);
void fulloutput(ostream &os, region_t const &reg);

// A pseudoregion is almost a region; it is a bbox that belongs to a
// certain component.  Pseudoregions are a compact way to store
// information about what components needs to send data to what other
// components during synchronisation or regridding.
struct pseudoregion_t {
  ibbox extent;
  int component;
  pseudoregion_t() {}
  pseudoregion_t(pseudoregion_t const &p)
      : extent(p.extent), component(p.component) {}
  pseudoregion_t(ibbox const &extent_, int const component_)
      : extent(extent_), component(component_) {}
};

MPI_Datatype mpi_datatype(pseudoregion_t const &) CCTK_ATTRIBUTE_PURE;
namespace dist {
template <> inline MPI_Datatype mpi_datatype<pseudoregion_t>() {
  pseudoregion_t dummy;
  return mpi_datatype(dummy);
}
}

bool operator==(pseudoregion_t const &a,
                pseudoregion_t const &b) CCTK_ATTRIBUTE_PURE;
inline bool operator!=(pseudoregion_t const &a, pseudoregion_t const &b) {
  return not(a == b);
}

inline size_t memoryof(pseudoregion_t const &p) {
  return memoryof(p.extent) + memoryof(p.component);
}

istream &operator>>(istream &is, pseudoregion_t &p);
ostream &operator<<(ostream &os, pseudoregion_t const &p);

struct sendrecv_pseudoregion_t {
  pseudoregion_t send, recv;
  sendrecv_pseudoregion_t() {}
  sendrecv_pseudoregion_t(sendrecv_pseudoregion_t const &srp)
      : send(srp.send), recv(srp.recv) {}
  sendrecv_pseudoregion_t(ibbox const &send_extent, int const send_component,
                          ibbox const &recv_extent, int const recv_component)
      : send(pseudoregion_t(send_extent, send_component)),
        recv(pseudoregion_t(recv_extent, recv_component)) {}
};

MPI_Datatype mpi_datatype(sendrecv_pseudoregion_t const &) CCTK_ATTRIBUTE_PURE;
namespace dist {
template <> inline MPI_Datatype mpi_datatype<sendrecv_pseudoregion_t>() {
  sendrecv_pseudoregion_t dummy;
  return mpi_datatype(dummy);
}
}

inline size_t memoryof(sendrecv_pseudoregion_t const &srp) {
  return memoryof(srp.send) + memoryof(srp.recv);
}

istream &operator>>(istream &os, sendrecv_pseudoregion_t &srp);
ostream &operator<<(ostream &os, sendrecv_pseudoregion_t const &srp);

#endif // #ifndef REGION_HH
