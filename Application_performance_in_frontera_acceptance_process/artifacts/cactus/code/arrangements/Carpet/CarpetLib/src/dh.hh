#ifndef DH_HH
#define DH_HH

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <set>
#include <map>
#include <string>
#include <vector>

#include "bbox.hh"
#include "bboxset.hh"
#include "defs.hh"
#include "gh.hh"
#include "region.hh"
#include "vect.hh"

using namespace std;

#define CARPET_HAVE_BUFFER_WIDTHS
#define CARPET_HAVE_OVERLAP_WIDTHS

// Forward declaration
class ggf;
class dh;

// A data hierarchy (grid hierarchy plus ghost zones)
class dh {

  static set<dh *> alldh;

  // Types
public:
  typedef vector<pseudoregion_t> pvect;
  typedef vector<sendrecv_pseudoregion_t> srpvect;

  struct light_dboxes {

    // Region description:

    ibbox exterior; // whole region (including boundaries)
    ibbox owned;    // evolved in time
    ibbox interior; // owned, plus outer boundary
    // Region statistics:
    size_type exterior_size, owned_size, active_size;

    size_t memory() const CCTK_MEMBER_ATTRIBUTE_PURE;
    istream &input(istream &is);
    ostream &output(ostream &os) const;
  };

  struct local_dboxes {

    // Information about the processor-local region:

    ibset buffers;  // buffer zones
    ibset overlaps; // overlap zones
    ibset active;   // owned minus (buffers + overlaps)
#if 0
    vector<ibset> buffers_stepped; // buffer zones [substep]
#endif

// Mask
#if 0
    ibset fine_active;
#endif
    vector<ibset> prolongation_boundary;
    vector<ibset> restriction_boundary;

    ibset restricted_region; // filled by restriction
    // Does not influence anything but the restricted_region
    ibset unused_region;

    // Refluxing (these live on the respective faces)
    vect<vect<ibset, 2>, dim> coarse_boundary;
    vect<vect<ibset, 2>, dim> fine_boundary;
#if 0 // OFFSETS,SIZE
    vect<vect<vector<int>,2>,dim> coarse_boundary_offsets;
    vect<vect<vector<int>,2>,dim> fine_boundary_offsets;
    vect<vect<int,2>,dim> coarse_boundary_size;
    vect<vect<int,2>,dim> fine_boundary_size;
#endif

    size_t memory() const CCTK_MEMBER_ATTRIBUTE_PURE;
    istream &input(istream &is);
    ostream &output(ostream &os) const;
  };

  struct level_dboxes {

    // Level description:

    // ibset exterior;
    // ibset outer_boundaries;
    // ibset communicated;
    // ibset boundaries;
    // ibset owned;
    ibset buffers;
    // ibset overlaps;
    ibset active;
    // ibset bndref;

    size_t memory() const CCTK_MEMBER_ATTRIBUTE_PURE;
    istream &input(istream &is);
    ostream &output(ostream &os) const;
  };

  struct full_dboxes {

    // Complete region description:

    ibbox exterior; // whole region (including boundaries)

    b2vect is_outer_boundary;
    ibset outer_boundaries; // outer boundary
    ibbox communicated;     // exterior without outer boundary

    ibset boundaries; // ghost zones
    ibbox owned;      // evolved in time

    ibset buffers;  // buffer zones
    ibset overlaps; // overlap zones
    ibset active;   // owned minus (buffers plus overlaps)

    ibset sync;   // filled by synchronisation
    ibset bndref; // filled by boundary prolongation

    // For Cactus: (these are like boundary or owned, but include the
    // outer boundary)
    ibset ghosts;   // ghost zones, as seen from Cactus
    ibbox interior; // interior (without ghost zones)

    bool operator==(full_dboxes const &b) const CCTK_MEMBER_ATTRIBUTE_PURE;
    bool operator!=(full_dboxes const &b) const { return not operator==(b); }

    size_t memory() const CCTK_MEMBER_ATTRIBUTE_PURE;
    istream &input(istream &is);
    ostream &output(ostream &os) const;
  };

  struct fast_dboxes {

    // Communication schedule:

    srpvect fast_mg_rest_sendrecv;
    srpvect fast_mg_prol_sendrecv;
    srpvect fast_ref_prol_sendrecv;
    srpvect fast_ref_rest_sendrecv;
    srpvect fast_sync_sendrecv;
    srpvect fast_ref_bnd_prol_sendrecv;

    // refluxing
    // vect<vect<srpvect,2>,dim> fast_ref_refl_sendrecv;
    srpvect fast_ref_refl_sendrecv_0_0;
    srpvect fast_ref_refl_sendrecv_0_1;
    srpvect fast_ref_refl_sendrecv_1_0;
    srpvect fast_ref_refl_sendrecv_1_1;
    srpvect fast_ref_refl_sendrecv_2_0;
    srpvect fast_ref_refl_sendrecv_2_1;
    // Note: Unfortunately we can't use an array of srpvects since
    // this doesn't work with C++ member pointers.  We instead define
    // them explicitly above (bah!), and maintain a static array with
    // member pointers for easier access.
    static vect<vect<srpvect fast_dboxes::*, 2>, dim> fast_ref_refl_sendrecv;
    static void init_fast_ref_refl_sendrecv();

    // vect<vect<srpvect,2>,dim> fast_ref_refl_prol_sendrecv;
    srpvect fast_ref_refl_prol_sendrecv_0_0;
    srpvect fast_ref_refl_prol_sendrecv_0_1;
    srpvect fast_ref_refl_prol_sendrecv_1_0;
    srpvect fast_ref_refl_prol_sendrecv_1_1;
    srpvect fast_ref_refl_prol_sendrecv_2_0;
    srpvect fast_ref_refl_prol_sendrecv_2_1;
    static vect<vect<srpvect fast_dboxes::*, 2>, dim>
        fast_ref_refl_prol_sendrecv;
    static void init_fast_ref_refl_prol_sendrecv();

    // Regridding schedule:

    bool do_init; // the srpvects below are only defined
                  // if this is true
    srpvect fast_old2new_sync_sendrecv;
    srpvect fast_old2new_ref_prol_sendrecv;

    bool operator==(fast_dboxes const &b) const CCTK_MEMBER_ATTRIBUTE_PURE;
    bool operator!=(fast_dboxes const &b) const { return not operator==(b); }

    size_t memory() const CCTK_MEMBER_ATTRIBUTE_PURE;
    istream &input(istream &is);
    ostream &output(ostream &os) const;
  };

  typedef vector<light_dboxes> light_cboxes; // ... for each component
  typedef vector<light_cboxes> light_rboxes; // ... for each refinement level
  typedef vector<light_rboxes> light_mboxes; // ... for each multigrid level

  typedef vector<local_dboxes> local_cboxes; // ... for each component
  typedef vector<local_cboxes> local_rboxes; // ... for each refinement level
  typedef vector<local_rboxes> local_mboxes; // ... for each multigrid level

  typedef vector<level_dboxes> level_rboxes; // ... for each refinement level
  typedef vector<level_rboxes> level_mboxes; // ... for each multigrid level

  typedef vector<full_dboxes> full_cboxes; // ... for each component
  typedef vector<full_cboxes> full_rboxes; // ... for each refinement level
  typedef vector<full_rboxes> full_mboxes; // ... for each multigrid level

  typedef vector<fast_dboxes> fast_rboxes; // ... for each refinement level
  typedef vector<fast_rboxes> fast_mboxes; // ... for each multigrid level

private:
  void setup_bboxes();

public: // should be readonly
  // Fields
  gh &h; // hierarchy

#if 0
  i2vect ghost_width;           // number of ghost zones
  i2vect buffer_width;          // number of buffer zones
  int prolongation_order_space; // order of spatial prolongation operator
#endif
  vector<i2vect> ghost_widths;           // number of ghost zones [rl]
  vector<i2vect> buffer_widths;          // number of buffer zones [rl]
  vector<i2vect> overlap_widths;         // number of overlap zones [rl]
  vector<int> prolongation_orders_space; // order of spatial
                                         // prolongation operator [rl]

  light_mboxes light_boxes; // grid hierarchy [ml][rl][c]
  local_mboxes local_boxes; // grid hierarchy [ml][rl][lc]
  level_mboxes level_boxes; // grid hierarchy [ml][rl]
  fast_mboxes fast_boxes;   // grid hierarchy [ml][rl][p]

private:
  // this needs to be sorted by varindex so that when iterating through the
  // container in order with a forward iterator, vector leaders are processed
  // first
  map<int, ggf *> gfs; // all grid functions

public:
  // Constructors
  dh(gh &h, vector<i2vect> const &ghost_widths,
     vector<i2vect> const &buffer_widths, vector<i2vect> const &overlap_widths,
     vector<int> const &prolongation_orders_space);

  // Destructors
  ~dh();

  // Helpers
  int prolongation_stencil_size(int rl) const CCTK_MEMBER_ATTRIBUTE_PURE;

  // Modifiers
  void regrid(bool do_init);
  void regrid_free(bool do_init);
  void recompose(int rl, bool do_prolongate);

private:
  int this_proc(int rl, int c) const CCTK_MEMBER_ATTRIBUTE_PURE;
  bool on_this_proc(int rl, int c) const CCTK_MEMBER_ATTRIBUTE_PURE;
  bool on_this_proc(int rl, int c, int cc) const CCTK_MEMBER_ATTRIBUTE_PURE;
  int this_oldproc(int rl, int c) const CCTK_MEMBER_ATTRIBUTE_PURE;
  bool on_this_oldproc(int rl, int c) const CCTK_MEMBER_ATTRIBUTE_PURE;

  static void broadcast_schedule(vector<fast_dboxes> &fast_level_otherprocs,
                                 fast_dboxes &fast_level,
                                 srpvect fast_dboxes::*const schedule_item);

public:
  // Grid function management
  void insert(ggf *f);
  void erase(ggf *f);

  // Output
  size_t memory() const CCTK_MEMBER_ATTRIBUTE_PURE;
  static size_t allmemory() CCTK_MEMBER_ATTRIBUTE_PURE;
  ostream &output(ostream &os) const;
};

MPI_Datatype mpi_datatype(dh::light_dboxes const &) CCTK_ATTRIBUTE_CONST;
MPI_Datatype mpi_datatype(dh::fast_dboxes const &) CCTK_ATTRIBUTE_CONST;
namespace dist {
template <> inline MPI_Datatype mpi_datatype<dh::light_dboxes>() {
  dh::light_dboxes dummy;
  return mpi_datatype(dummy);
}
template <> inline MPI_Datatype mpi_datatype<dh::fast_dboxes>() {
  dh::fast_dboxes dummy;
  return mpi_datatype(dummy);
}
}

inline size_t memoryof(dh::light_dboxes const &b) { return b.memory(); }

inline size_t memoryof(dh::local_dboxes const &b) { return b.memory(); }

inline size_t memoryof(dh::level_dboxes const &b) { return b.memory(); }

inline size_t memoryof(dh::full_dboxes const &b) { return b.memory(); }

inline size_t memoryof(dh::fast_dboxes const &b) { return b.memory(); }

inline size_t memoryof(dh const &d) { return d.memory(); }

inline istream &operator>>(istream &is, dh::light_dboxes &b) {
  return b.input(is);
}

inline istream &operator>>(istream &is, dh::local_dboxes &b) {
  return b.input(is);
}

inline istream &operator>>(istream &is, dh::level_dboxes &b) {
  return b.input(is);
}

inline istream &operator>>(istream &is, dh::full_dboxes &b) {
  return b.input(is);
}

inline istream &operator>>(istream &is, dh::fast_dboxes &b) {
  return b.input(is);
}

inline ostream &operator<<(ostream &os, dh::light_dboxes const &b) {
  return b.output(os);
}

inline ostream &operator<<(ostream &os, dh::local_dboxes const &b) {
  return b.output(os);
}

inline ostream &operator<<(ostream &os, dh::level_dboxes const &b) {
  return b.output(os);
}

inline ostream &operator<<(ostream &os, dh::full_dboxes const &b) {
  return b.output(os);
}

inline ostream &operator<<(ostream &os, dh::fast_dboxes const &b) {
  return b.output(os);
}

inline ostream &operator<<(ostream &os, dh const &d) { return d.output(os); }

#endif // DH_HH
