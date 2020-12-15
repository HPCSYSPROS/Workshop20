#ifndef TH_HH
#define TH_HH

#include <cctk.h>

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "defs.hh"
#include "gh.hh"

using namespace std;

// Forward declaration
class th;

// Input
istream &operator>>(istream &is, th &t);
// Output
ostream &operator<<(ostream &os, const th &t);

// The time hierarchy (information about the current time)
class th {

  static set<th *> allth;

public: // should be readonly
  // Fields
  gh &h; // hierarchy

  bool const time_interpolation_during_regridding;

  int timelevels; // const

private:
  vector<int> reffacts; // const

  vector<vector<vector<CCTK_REAL> > > times; // current times [ml][rl][tl]
  vector<vector<CCTK_REAL> > deltas;         // time steps [ml][rl]

public:
  // Constructors
  th(gh &h, bool time_interpolation_during_regridding);

  th(gh &h, int timelevels, vector<int> const &reffacts,
     bool time_interpolation_during_regridding);

  // Destructors
  ~th();

  // Modifiers
  void regrid();
  void regrid_free();

  // Time management
  void set_time(int const ml, int const rl, int const tl, CCTK_REAL const &t) {
    assert(ml >= 0 and ml < h.mglevels());
    assert(rl >= 0 and rl < h.reflevels());
    assert(tl >= 0 and tl < timelevels);
    // assert (isfinite(t));
    times.AT(ml).AT(rl).AT(tl) = t;
  }

  CCTK_REAL get_time(int const ml, int const rl, int const tl) const {
    assert(ml >= 0 and ml < h.mglevels());
    assert(rl >= 0 and rl < h.reflevels());
    assert(tl >= 0 and tl < timelevels);
    CCTK_REAL const t = times.AT(ml).AT(rl).AT(tl);
    // assert (isfinite(t));
    return t;
  }

  void set_delta(int const ml, int const rl, CCTK_REAL const &dt) {
    assert(ml >= 0 and ml < h.mglevels());
    assert(rl >= 0 and rl < h.reflevels());
    // assert (isfinite(dt));
    deltas.AT(ml).AT(rl) = dt;
  }

  CCTK_REAL get_delta(int const ml, int const rl) const {
    assert(ml >= 0 and ml < h.mglevels());
    assert(rl >= 0 and rl < h.reflevels());
    CCTK_REAL const dt = deltas.AT(ml).AT(rl);
    // assert (isfinite(dt));
    return dt;
  }

  void advance_time(int const ml, int const rl);
  void retreat_time(int const ml, int const rl);
  void flip_timelevels(int const ml, int const rl);

  // Output
  size_t memory() const CCTK_MEMBER_ATTRIBUTE_PURE;
  static size_t allmemory() CCTK_MEMBER_ATTRIBUTE_PURE;
  istream &input(istream &is);
  ostream &output(ostream &os) const;
};

inline size_t memoryof(th const &t) { return t.memory(); }
inline istream &operator>>(istream &is, th &t) { return t.input(is); }
inline ostream &operator<<(ostream &os, const th &t) { return t.output(os); }

#endif // TH_HH
