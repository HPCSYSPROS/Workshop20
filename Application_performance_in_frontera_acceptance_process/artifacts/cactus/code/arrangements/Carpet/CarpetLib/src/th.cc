#include <cctk.h>

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "defs.hh"
#include "gh.hh"

#include "th.hh"

using namespace std;

set<th *> th::allth;

// Constructors
th::th(gh &h_, bool const time_interpolation_during_regridding_)
    : h(h_), time_interpolation_during_regridding(
                 time_interpolation_during_regridding_),
      timelevels(0) {
  reffacts.resize(1, 1);
  allth.insert(this);
  h.insert(this);
}

th::th(gh &h_, int const timelevels_, vector<int> const &reffacts_,
       bool const time_interpolation_during_regridding_)
    : h(h_), time_interpolation_during_regridding(
                 time_interpolation_during_regridding_),
      timelevels(timelevels_), reffacts(reffacts_) {
  assert(timelevels_ > 0);
  assert(reffacts.size() >= 1);
  assert(reffacts.front() == 1);
  for (size_t n = 1; n < reffacts.size(); ++n) {
    assert(reffacts.AT(n) >= reffacts.AT(n - 1));
    assert(reffacts.AT(n) % reffacts.AT(n - 1) == 0);
  }
  allth.insert(this);
  h.insert(this);
}

// Destructors
th::~th() {
  h.erase(this);
  allth.erase(this);
}

// Modifiers
void th::regrid() {
  CCTK_REAL const basetime = 0.0;
  CCTK_REAL const basedelta = 1.0;

  const int old_mglevels = times.size();
  times.resize(h.mglevels());
  deltas.resize(h.mglevels());
  for (int ml = 0; ml < h.mglevels(); ++ml) {
    const int old_reflevels = times.AT(ml).size();
    times.AT(ml).resize(h.reflevels());
    deltas.AT(ml).resize(h.reflevels());
    for (int rl = 0; rl < h.reflevels(); ++rl) {
      if (ml == 0) {
        deltas.AT(ml).AT(rl) = basedelta / reffacts.AT(rl);
      } else {
        deltas.AT(ml).AT(rl) = deltas.AT(ml - 1).AT(rl) * h.mgfact;
      }
      if (old_mglevels == 0) {
        times.AT(ml).AT(rl).resize(timelevels);
        for (int tl = 0; tl < timelevels; ++tl) {
          times.AT(ml).AT(rl).AT(tl) = basetime - tl * deltas.AT(ml).AT(rl);
        }
      } else if (rl < old_reflevels) {
        // do nothing
      } else if (rl == 0) {
        assert(0);
        times.AT(ml).AT(rl).resize(timelevels);
        for (int tl = 0; tl < timelevels; ++tl) {
          times.AT(ml).AT(rl).AT(tl) = basetime - tl * deltas.AT(ml).AT(rl);
        }
      } else {
        if (time_interpolation_during_regridding) {
          // We probably don't want to do this, but it is nice for
          // compatibility
          times.AT(ml).AT(rl).resize(timelevels);
          for (int tl = 0; tl < timelevels; ++tl) {
            // linear interpolation between the two surrounding coarse
            // grid times
            assert(reffacts.AT(rl) % reffacts.AT(rl - 1) == 0);
            int const rf = reffacts.AT(rl) / reffacts.AT(rl - 1);
            int const ctl = tl / rf;
            int const mtl = tl % rf;
            if (mtl == 0) {
              times.AT(ml).AT(rl).AT(tl) = times.AT(ml).AT(rl - 1).AT(ctl);
            } else {
              assert(ctl + 1 < timelevels);
              CCTK_REAL const alpha = (CCTK_REAL)mtl / rf;
              assert(alpha > 0 and alpha < 1);
              times.AT(ml).AT(rl).AT(tl) =
                  (1 - alpha) * times.AT(ml).AT(rl - 1).AT(ctl) +
                  (alpha)*times.AT(ml).AT(rl - 1).AT(ctl + 1);
            }
          }
        } else {
          times.AT(ml).AT(rl) = times.AT(ml).AT(rl - 1);
        }
      }
    }
  }
  for (int ml = 0; ml < h.mglevels(); ++ml) {
    for (int rl = 0; rl < h.reflevels(); ++rl) {
      for (int tl = 1; tl < timelevels; ++tl) {
        assert(times.AT(ml).AT(rl).AT(tl) < times.AT(ml).AT(rl).AT(tl - 1));
      }
    }
  }
  for (int ml = 0; ml < h.mglevels(); ++ml) {
    for (int rl = 0; rl < h.reflevels(); ++rl) {
      // assert (isfinite(deltas.AT(ml).AT(rl)));
      for (int tl = 0; tl < timelevels; ++tl) {
        // assert (isfinite(times.AT(ml).AT(rl).AT(tl)));
      }
    }
  }
}

void th::regrid_free() {}

void th::advance_time(int const ml, int const rl) {
  for (int tl = timelevels - 1; tl > 0; --tl) {
    set_time(ml, rl, tl, get_time(ml, rl, tl - 1));
  }
  set_time(ml, rl, 0, get_time(ml, rl, 0) + get_delta(ml, rl));
}

void th::retreat_time(int const ml, int const rl) {
  CCTK_REAL const t = get_time(ml, rl, 0);
  for (int tl = 0; tl < timelevels - 1; ++tl) {
    set_time(ml, rl, tl, get_time(ml, rl, tl + 1));
  }
  set_time(ml, rl, timelevels - 1, t);
}

void th::flip_timelevels(int const ml, int const rl) {
  for (int tl = 1; tl < (timelevels + 1) / 2; ++tl) {
    int const tl2 = timelevels - tl;
    CCTK_REAL const t = get_time(ml, rl, tl);
    CCTK_REAL const t2 = get_time(ml, rl, tl2);
    set_time(ml, rl, tl, t2);
    set_time(ml, rl, tl2, t);
  }
  set_delta(ml, rl, -get_delta(ml, rl));
}

// Memory usage
size_t th::memory() const {
  return memoryof(reffacts) + memoryof(times) + memoryof(deltas);
}

size_t th::allmemory() {
  size_t mem = memoryof(allth);
  for (set<th *>::const_iterator thi = allth.begin(); thi != allth.end();
       ++thi) {
    mem += memoryof(**thi);
  }
  return mem;
}

// Input
istream &th::input(istream &is) {
  skipws(is);
  consume(is, "th:{");
  consume(is, "timelevels=");
  is >> timelevels;
  consume(is, ",");
  consume(is, "reffacts=");
  is >> reffacts;
  consume(is, ",");
  consume(is, "times=");
  is >> times;
  consume(is, ",");
  consume(is, "deltas=");
  is >> deltas;
  consume(is, "}");
  return is;
}

// Output
ostream &th::output(ostream &os) const {
  os << "th:{"
     << "timelevels=" << timelevels << ","
     << "reffacts=" << reffacts << ","
     << "times=" << times << ","
     << "deltas=" << deltas << "}";
  return os;
}
