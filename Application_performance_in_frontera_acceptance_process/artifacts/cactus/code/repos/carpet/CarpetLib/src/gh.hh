#ifndef GH_HH
#define GH_HH

#include <cassert>
#include <iostream>
#include <set>
#include <vector>

#include "bbox.hh"
#include "bboxset.hh"
#include "defs.hh"
#include "dist.hh"
#include "region.hh"
#include "vect.hh"

using namespace std;

// Forward declaration
class dh;
class th;
class gh;

// A refinement hierarchy, where higher levels are finer than the base
// level.  The extents do not include ghost zones.
class gh {

  static set<gh *> allgh;

public:
  // Types
  typedef vector<region_t> cregs; // ... for each component
  typedef vector<cregs> rregs;    // ... for each refinement level
  typedef vector<rregs> mregs;    // ... for each multigrid level

public: // should be readonly
  // Fields
  const vector<ivect> reffacts; // refinement factors
  const centering refcent;      // vertex or cell centered

  const int mgfact;       // default multigrid factor
  const centering mgcent; // default (vertex or cell centered)

  vector<vector<ibbox> > baseextents; // [ml][rl]
  const i2vect boundary_width;

private:
  vector<vector<int> > global_components_;     // [rl][lc]
  vector<vector<int> > local_components_;      // [rl][c]
  vector<vector<int> > old_global_components_; // [rl][lc]
  vector<vector<int> > old_local_components_;  // [rl][c]
public:
  // Extents of the regions before distributing them over the
  // processes
  rregs superregions;

  mregs regions;    // extents and properties of all grids
  mregs oldregions; // extents and properties of all grids

  set<dh *> dhs; // all data hierarchies
  set<th *> ths; // all time hierarchies

public:
  // Constructors
  gh(vector<ivect> const &reffacts, centering refcent, int mgfact,
     centering mgcent, vector<vector<ibbox> > const &baseextents,
     i2vect const &boundary_width);

  // Destructors
  ~gh();

  // Modifiers
  void regrid(rregs const &superregs, mregs const &regs, bool do_init);
  void regrid_free(bool do_init);
  bool recompose(int rl, bool do_prolongate);

private:
  bool level_did_change(int rl) const CCTK_MEMBER_ATTRIBUTE_PURE;

  // Accessors

public:
  ibbox const &extent(const int ml, const int rl, const int c) const {
    return regions.AT(ml).AT(rl).AT(c).extent;
  }

  ibbox const &baseextent(const int ml, const int rl) const {
    return baseextents.AT(ml).AT(rl);
  }

  b2vect const &outer_boundaries(const int rl, const int c) const {
    return regions.AT(0).AT(rl).AT(c).outer_boundaries;
  }

  int processor(const int rl, const int c) const CCTK_MEMBER_ATTRIBUTE_PURE {
    return regions.AT(0).AT(rl).AT(c).processor;
  }

  int old_processor(const int rl,
                    const int c) const CCTK_MEMBER_ATTRIBUTE_PURE {
    return oldregions.AT(0).AT(rl).AT(c).processor;
  }

  int mglevels() const CCTK_MEMBER_ATTRIBUTE_PURE {
    return (int)regions.size();
  }

  int reflevels() const CCTK_MEMBER_ATTRIBUTE_PURE {
    if (mglevels() == 0)
      return 0;
    return (int)regions.AT(0).size();
  }

  int components(const int rl) const CCTK_MEMBER_ATTRIBUTE_PURE {
    return (int)regions.AT(0).AT(rl).size();
  }

  bool is_local(const int rl, const int c) const CCTK_MEMBER_ATTRIBUTE_PURE {
    return processor(rl, c) == dist::rank();
  }

  int local_components(int rl) const CCTK_MEMBER_ATTRIBUTE_PURE;
  int get_component(int rl, int lc) const CCTK_MEMBER_ATTRIBUTE_PURE;
  int get_local_component(int rl, int c) const CCTK_MEMBER_ATTRIBUTE_PURE;

  bool old_is_local(const int rl,
                    const int c) const CCTK_MEMBER_ATTRIBUTE_PURE {
    return old_processor(rl, c) == dist::rank();
  }

  int old_local_components(int rl) const CCTK_MEMBER_ATTRIBUTE_PURE;
  int get_old_component(int rl, int lc) const CCTK_MEMBER_ATTRIBUTE_PURE;
  int get_old_local_component(int rl, int c) const CCTK_MEMBER_ATTRIBUTE_PURE;

#if 0
  // Convert between index positions and coordinate positions
  rvect ipos2rpos (ivect const & ipos,
                   rvect const & origin, rvect const & scale,
                   int const ml, int const rl) const CCTK_MEMBER_ATTRIBUTE_PURE;
  ivect rpos2ipos (rvect const & rpos,
                   rvect const & origin, rvect const & scale,
                   int const ml, int const rl) const CCTK_MEMBER_ATTRIBUTE_PURE;
  ivect rpos2ipos1 (rvect const & rpos,
                    rvect const & origin, rvect const & scale,
                    int const ml, int const rl) const
    CCTK_MEMBER_ATTRIBUTE_PURE;
#endif

  void locate_position(rvect const &rpos, int const ml, int const minrl,
                       int const maxrl, int &rl, int &c,
                       ivect &aligned_ipos) const CCTK_MEMBER_ATTRIBUTE_PURE;

  void locate_position(ivect const &ipos, int const ml, int const minrl,
                       int const maxrl, int &rl, int &c,
                       ivect &aligned_ipos) const CCTK_MEMBER_ATTRIBUTE_PURE;

  // Time hierarchy management
  void insert(th *t);
  void erase(th *t);

  // Data hierarchy management
  void insert(dh *d);
  void erase(dh *d);

  // Output
  size_t memory() const CCTK_MEMBER_ATTRIBUTE_PURE;
  static size_t allmemory() CCTK_MEMBER_ATTRIBUTE_PURE;
  ostream &output(ostream &os) const;

private:
  void do_output_bboxes(ostream &os) const;
  void do_output_bases(ostream &os) const;
};

inline size_t memoryof(gh const &g) { return g.memory(); }

inline ostream &operator<<(ostream &os, gh const &h) { return h.output(os); }

#endif // GH_HH
