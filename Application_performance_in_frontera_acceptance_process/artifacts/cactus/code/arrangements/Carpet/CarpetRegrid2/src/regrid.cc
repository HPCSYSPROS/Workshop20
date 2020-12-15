#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <typeinfo>
#include <vector>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <Timer.hh>

#include <bbox.hh>
#include <bboxset.hh>
#include <defs.hh>
#include <dh.hh>
#include <gh.hh>
#include <region.hh>
#include <vect.hh>

#include <carpet.hh>

#include "amr.hh"
#include "boundary.hh"
#include "indexing.hh"
#include "property.hh"

namespace CarpetRegrid2 {

using namespace std;
using namespace Carpet;

struct centre_description {
  int _num_levels;
  int _active;
  rvect _position;
  vector<rvect> _radius;

  centre_description(cGH const *cctkGH, int n);
};

centre_description::centre_description(cGH const *const cctkGH, int const n) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  assert(n >= 0 and n < num_centres);

  bool found_error = false;

  int lsh[2];
  getvectorindex2(cctkGH, "CarpetRegrid2::radii", lsh);

  this->_num_levels = num_levels[n];
  this->_active = active[n];
  this->_position = rvect(position_x[n], position_y[n], position_z[n]);
  if (any(not isfinite(this->_position))) {
    CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
               "The position of region %d is [%g,%g,%g], which is not finite",
               n + 1, double(this->_position[0]), double(this->_position[1]),
               double(this->_position[2]));
    found_error = true;
  }
  this->_radius.resize(this->_num_levels);
  this->_radius.at(0) = rvect(-1.0, -1.0, -1.0); // unused
  for (int rl = 1; rl < this->_num_levels; ++rl) {
    int const ind = index2(lsh, rl, n);
    CCTK_REAL const rx = radius_x[ind] < 0.0 ? radius[ind] : radius_x[ind];
    CCTK_REAL const ry = radius_y[ind] < 0.0 ? radius[ind] : radius_y[ind];
    CCTK_REAL const rz = radius_z[ind] < 0.0 ? radius[ind] : radius_z[ind];
    rvect const rad(rx, ry, rz);
    if (any(not isfinite(rad))) {
      CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "The radius of refinement level %d of region %d is "
                 "[%g,%g,%g], which is not finite",
                 rl, n + 1, double(rad[0]), double(rad[1]), double(rad[2]));
      found_error = true;
    }
    if (any(rad < CCTK_REAL(0))) {
      CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "The radius of refinement level %d of region %d is "
                 "[%g,%g,%g], which is not non-negative",
                 rl, n + 1, double(rad[0]), double(rad[1]), double(rad[2]));
      found_error = true;
    }
    this->_radius.at(rl) = rad;
  }

  if (this->_num_levels > maxreflevels) {
    CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Region %d has %d levels active, which is larger than the "
               "maximum number of refinement levels %d",
               n + 1, this->_num_levels, maxreflevels);
    found_error = true;
  }

  if (found_error) {
    CCTK_WARN(CCTK_WARN_ABORT, "Errors found in grid structure specification");
  }
}

extern "C" {
CCTK_INT
CarpetRegrid2_Regrid(CCTK_POINTER_TO_CONST const cctkGH_,
                     CCTK_POINTER const superregss_, CCTK_POINTER const regsss_,
                     CCTK_INT const force);

CCTK_INT
CarpetRegrid2_RegridMaps(CCTK_POINTER_TO_CONST const cctkGH_,
                         CCTK_POINTER const superregsss_,
                         CCTK_POINTER const regssss_, CCTK_INT const force);
}

void Regrid(cGH const *const cctkGH, gh::rregs &regss) {
  DECLARE_CCTK_PARAMETERS;

  if (verbose or veryverbose) {
    CCTK_VInfo(CCTK_THORNSTRING,
               "Regridding level %d map %d at iteration %d time %g", reflevel,
               Carpet::map, cctkGH->cctk_iteration, cctkGH->cctk_time);
  }

  assert(is_singlemap_mode());
  gh const &hh = *vhh.at(Carpet::map);
  dh const &dd = *vdd.at(Carpet::map);

  //
  // Find extent of domain
  //

  // This requires that CoordBase is used
  // TODO: check this (for Carpet, and maybe also for CartGrid3D)
  // TODO: (the check that these two are consistent should be in
  // Carpet)

  jjvect nboundaryzones, is_internal;
  jjvect is_staggered, shiftout;
  get_boundary_specification(nboundaryzones, is_internal, is_staggered,
                             shiftout);

  rvect physical_lower, physical_upper;
  rvect spacing;
  get_physical_boundary(physical_lower, physical_upper, spacing);

  // Adapt spacing for convergence level
  spacing *= ipow((CCTK_REAL)mgfact, basemglevel);

  rvect exterior_lower, exterior_upper;
  calculate_exterior_boundary(physical_lower, physical_upper, exterior_lower,
                              exterior_upper, spacing);

  //
  // Determine which refinement levels may be changed
  //
  int min_rl = 1; // we cannot change the coarsest level
  if (freeze_unaligned_levels or freeze_unaligned_parent_levels) {
    while (min_rl < int(regss.size())) {
// Increase min_rl until we find a level that can be changed
#if 0
        // TODO: think about this a bit more
        // TODO: use this taper-checking also in Comm.cc
        bool in_sync = true;
        if (freeze_unaligned_parent_levels) {
          int const parent_do_every =
            ipow(mgfact, mglevel) *
            (maxtimereflevelfact / timereffacts.at(min_rl-1));
          in_sync =
            cctkGH->cctk_iteration == 0 or
            (cctkGH->cctk_iteration-1) % parent_do_every == 0;
        }
#else
      bool in_sync = true;
      if (freeze_unaligned_parent_levels) {
        // Assume that non-existing levels are always aligned
        if (min_rl < reflevels) {
          CCTK_REAL const mytime = tt->get_time(mglevel, min_rl, 0);
          CCTK_REAL const parenttime = tt->get_time(mglevel, min_rl - 1, 0);
          CCTK_REAL const eps = 1.0e-12;
          in_sync = abs(mytime - parenttime) <= eps * abs(delta_time);
        }
      }
#endif
      int const do_every = ipow(mgfact, mglevel) *
                           (maxtimereflevelfact / timereffacts.at(min_rl));
      bool const is_active = cctkGH->cctk_iteration == 0 or
                             (cctkGH->cctk_iteration - 1) % do_every == 0;
      if (is_active and in_sync)
        break;
      ++min_rl;
    }
    if (verbose or veryverbose) {
      CCTK_VInfo(CCTK_THORNSTRING, "Regridding levels %d and up", min_rl);
    }
  }

  //
  // Calculate the union of the bounding boxes for all levels
  //

  // The set of refined regions
  vector<ibset> regions(min_rl);

  // Set up coarse levels
  for (int rl = 0; rl < min_rl; ++rl) {
    if (veryverbose) {
      cout << "Refinement level " << rl << ": will not be changed" << endl;
    }
    for (size_t c = 0; c < regss.at(rl).size(); ++c) {
      regions.at(rl) += regss.at(rl).at(c).extent;
    }
    if (veryverbose) {
      cout << "Refinement level " << rl << ": regions are " << regions.at(rl)
           << endl;
    }
  }

  // Refine only patch 0
  if (Carpet::map == 0) {

    // Loop over all centres
    for (int n = 0; n < num_centres; ++n) {
      centre_description centre(cctkGH, n);
      if (centre._active) {

        // Loop over all levels for this centre
        for (int rl = min_rl; rl < centre._num_levels; ++rl) {

          if (veryverbose) {
            cout << "Refinement level " << rl
                 << ": importing refined region settings..." << endl;
          }

          level_boundary const bnd(hh, dd, rl);

          // Calculate a bbox for this region
          rvect const rmin = centre._position - centre._radius.at(rl);
          rvect const rmax = centre._position + centre._radius.at(rl);

          if (veryverbose) {
            cout << "Centre " << n + 1 << " refinement level " << rl
                 << ": coordinate region is (" << rmin << ":" << rmax << ")\n";
          }

          // Convert to an integer bbox
          ivect const istride = hh.baseextent(0, rl).stride();
          ivect const imin = rpos2ipos(rmin, bnd.origin, bnd.scale, hh, rl) -
                             int(boundary_shiftout) * istride;
          ivect const imax = rpos2ipos1(rmax, bnd.origin, bnd.scale, hh, rl) +
                             int(boundary_shiftout) * istride;

          if (veryverbose) {
            cout << "Centre " << n + 1 << " refinement level " << rl
                 << ": integer region is (" << imin << ":" << imax << ")\n";
          }

          ibbox const region(imin, imax, istride);

          // Add this region to the list of regions
          if (static_cast<int>(regions.size()) < rl + 1) {
            regions.resize(rl + 1);
          }
          regions.at(rl) |= region;

          if (veryverbose) {
            cout << "Refinement level " << rl << ": preliminary regions are "
                 << regions.at(rl) << endl;
          }

        } // for rl

      } // if centre is active
    }   // for n

    if (adaptive_refinement) {
      // Loop over all levels
      for (int rl = min_rl; rl < min(maxreflevels, reflevels + 1); ++rl) {
        if (static_cast<int>(regions.size()) < rl + 1) {
          regions.resize(rl + 1);
        }
        evaluate_level_mask(cctkGH, regions, rl);
        if (regions.at(rl).empty()) {
          // If there are no refined regions on this level, truncate
          // the refinement hierarchy, and stop
          assert(static_cast<int>(regions.size()) == rl + 1);
          regions.resize(rl);
          break;
        }
      }
    }

  } // if map==0

  // We need to check and/or enforce certain properties of the grid
  // structure, such as e.g. proper nesting. Unfortunately,
  // enforcing one of the properties may invalidate another. We
  // therefore abstract the concept of a property into a class
  // "property", and enforce them round-robin until the grid
  // structure does not change any more. The order in which the
  // properties are enforced is relevant (i.e. it influence the
  // final grid structure).

  // Properties to be applied (without testing), only once, in the
  // beginning
  vector<property *> once_properties;
  once_properties.push_back(new proper_nesting());
  once_properties.push_back(new add_buffers());

  // Properties to be enforced "until all is well"
  vector<property *> properties;
  properties.push_back(new combine_regions());
  properties.push_back(new snap_coarse());
  properties.push_back(new rotsym90());
  properties.push_back(new rotsym180());
  properties.push_back(new parsym());
  properties.push_back(new periodic<0>());
  properties.push_back(new periodic<1>());
  properties.push_back(new periodic<2>());
  properties.push_back(new boundary_clip());

  // Properties to be tested (and not enforced) in the end
  vector<property *> final_properties;
  final_properties.push_back(new in_domain());
  final_properties.push_back(new is_symmetric());

  // TODO: Check that the coarse grid contains all finer grids. To
  // do this, apply proper_nesting to the coarse grid, and then
  // check whether it grew -- if so, this is an error.

  regss.resize(regions.size());

  // Loop over all levels from finest to coarsest
  assert(regions.size() > 0);
  for (int rl = regions.size() - 1; rl >= min_rl; --rl) {

    // Sanity check
    assert(not regions.at(rl).empty());

    // Determine boundary locations of this level
    level_boundary const bnd(hh, dd, rl);

    if (veryverbose) {
      cout << "Refinement level " << rl << ":\n"
           << "   Original regions are " << regions.at(rl) << "\n";
    }

    // Apply "once" properties unconditionally
    for (vector<property *>::iterator pi = once_properties.begin();
         pi != once_properties.end(); ++pi) {
      (*pi)->enforce(hh, dd, bnd, regions, rl);
    }

    // Enforce properties on this level
    for (int count = 0;; ++count) {
      CCTK_VInfo(CCTK_THORNSTRING,
                 "Enforcing grid structure properties, iteration %d", count);
      bool done_enforcing = true;
      ibset const old_regions = regions.at(rl);
      for (vector<property *>::iterator pi = properties.begin();
           pi != properties.end(); ++pi) {
        bool const test_satisfied = (*pi)->test(hh, dd, bnd, regions, rl);
        if (not test_satisfied) {
          (*pi)->enforce(hh, dd, bnd, regions, rl);
        }
        done_enforcing = done_enforcing and test_satisfied;
      }

      if (done_enforcing)
        break;

      if (regions.at(rl) == old_regions) {
        CCTK_WARN(CCTK_WARN_ABORT, "Could not enforce grid structure "
                                   "properties (not making any progress); "
                                   "giving up");
      }
      if (count == 10) {
        CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Could not enforce grid structure properties after %d "
                   "iterations; giving up",
                   count);
      }
      if (count != 0) {
        // This may not be true. However, the previous version of
        // the code assumed this, so we want to know whether this
        // fails.
        CCTK_WARN(CCTK_WARN_ALERT, "Could not enforce grid structure "
                                   "properties in this round, starting another "
                                   "round");
      }
    }

    if (veryverbose) {
      cout << "Refinement level " << rl << ":\n"
           << "   Final regions are " << regions.at(rl) << "\n";
    }

    // Test final properties
    bool all_is_well = true;
    for (vector<property *>::const_iterator pi = final_properties.begin();
         pi != final_properties.end(); ++pi) {
      bool const test_satisfied = (*pi)->test(hh, dd, bnd, regions, rl);
      if (not test_satisfied) {
        cout << "Necessary property " << typeid(*pi).name()
             << " does not hold\n";
      }
      all_is_well = all_is_well and test_satisfied;
    }
    if (not all_is_well) {
      CCTK_WARN(CCTK_WARN_ABORT,
                "Not all necessary grid structure properties are holding");
    }

    // Create a vector of bboxes for this level
    gh::cregs regs;
    for (ibset::const_iterator ibb = regions.AT(rl).begin();
         ibb != regions.AT(rl).end(); ++ibb) {
      ibbox const &bb = *ibb;
      assert(bb.is_contained_in(hh.baseextent(0, rl)));

      bvect const lower_is_outer = bb.lower() <= bnd.level_physical_ilower;
      bvect const upper_is_outer = bb.upper() >= bnd.level_physical_iupper;

      b2vect const ob(lower_is_outer, upper_is_outer);

      region_t reg;
      reg.extent = bb;
      reg.map = Carpet::map;
      reg.outer_boundaries = ob;
      regs.push_back(reg);
    }

    regss.AT(rl) = regs;

  } // for rl

  // Delete properties
  for (vector<property *>::iterator pi = once_properties.begin();
       pi != once_properties.end(); ++pi) {
    delete *pi;
  }
  for (vector<property *>::iterator pi = properties.begin();
       pi != properties.end(); ++pi) {
    delete *pi;
  }
  for (vector<property *>::iterator pi = final_properties.begin();
       pi != final_properties.end(); ++pi) {
    delete *pi;
  }
}

CCTK_INT
CarpetRegrid2_Regrid(CCTK_POINTER_TO_CONST const cctkGH_,
                     CCTK_POINTER const superregss_, CCTK_POINTER const regsss_,
                     CCTK_INT const force) {
  cGH const *const cctkGH = static_cast<cGH const *>(cctkGH_);

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  char const *const where = "CarpetRegrid2::Regrid";
  static Timers::Timer timer(where);
  timer.start();

  assert(is_singlemap_mode());

  // Decide whether to change the grid hierarchy
  bool do_recompose;
  if (force) {
    do_recompose = true;
  } else {
    if (regrid_every == -1) {
      do_recompose = false;
    } else if (regrid_every == 0) {
      do_recompose = cctk_iteration == 0;
    } else {
      // Regrid at most once per iteration
      do_recompose =
          (cctk_iteration == 0 or
           (cctk_iteration > 0 and (cctk_iteration - 1) % regrid_every == 0 and
            (cctk_iteration > *last_iteration or
             (cctk_iteration == *last_iteration and Carpet::map > *last_map))));
    }
  }

  if (verbose or veryverbose) {
    if (do_recompose) {
      for (int n = 0; n < num_centres; ++n) {
        if (active[n]) {
          CCTK_VInfo(CCTK_THORNSTRING,
                     "Centre %d is at position [%g,%g,%g] with %d levels",
                     n + 1, static_cast<double>(position_x[n]),
                     static_cast<double>(position_y[n]),
                     static_cast<double>(position_z[n]),
                     static_cast<int>(num_levels[n]));
        } else {
          CCTK_VInfo(CCTK_THORNSTRING, "Centre %d is not active", n + 1);
        }
      }
    }
  }

  if (not force and do_recompose and *last_iteration != -1) {

    int lsh[2];
    getvectorindex2(cctkGH, "CarpetRegrid2::radii", lsh);

    // Regrid only if the regions have changed sufficiently
    do_recompose = false;
    if (adaptive_refinement)
      do_recompose = true;
    if (not do_recompose) {
      for (int n = 0; n < num_centres; ++n) {

        // Regrid if a region became active or inactive
        do_recompose = active[n] != old_active[n];
        if (do_recompose)
          break;

        // Check only active regions
        if (not active[n])
          continue;

        // Regrid if the number of levels changed
        do_recompose = num_levels[n] != old_num_levels[n];
        if (do_recompose)
          break;

        // Regrid if the positions have changed sufficiently
        CCTK_REAL const dist2 = pow(position_x[n] - old_position_x[n], 2) +
                                pow(position_y[n] - old_position_y[n], 2) +
                                pow(position_z[n] - old_position_z[n], 2);
        CCTK_REAL mindist;
        switch (n) {
        case 0:
          mindist = movement_threshold_1;
          break;
        case 1:
          mindist = movement_threshold_2;
          break;
        case 2:
          mindist = movement_threshold_3;
          break;
        case 3:
          mindist = movement_threshold_4;
          break;
        case 4:
          mindist = movement_threshold_5;
          break;
        case 5:
          mindist = movement_threshold_6;
          break;
        case 6:
          mindist = movement_threshold_7;
          break;
        case 7:
          mindist = movement_threshold_8;
          break;
        case 8:
          mindist = movement_threshold_9;
          break;
        case 9:
          mindist = movement_threshold_10;
          break;
        default:
          assert(0);
        }
        do_recompose = dist2 > pow(mindist, 2);
        if (do_recompose)
          break;

        // Regrid if the radii have changed sufficiently
        for (int rl = 1; rl < num_levels[n]; ++rl) {
          int const ind = index2(lsh, rl, n);
          CCTK_REAL const rx = radius_x[ind] < 0 ? radius[ind] : radius_x[ind];
          CCTK_REAL const ry = radius_y[ind] < 0 ? radius[ind] : radius_y[ind];
          CCTK_REAL const rz = radius_z[ind] < 0 ? radius[ind] : radius_z[ind];
          rvect const rad(rx, ry, rz);
          rvect const oldrad(old_radius_x[ind], old_radius_y[ind],
                             old_radius_z[ind]);
          CCTK_REAL const drfac =
              sqrt(sum(ipow(rad - oldrad, 2))) / sqrt(sum(ipow(oldrad, 2)));
          CCTK_REAL mindrfac;
          switch (n) {
          case 0:
            mindrfac = radius_rel_change_threshold_1;
            break;
          case 1:
            mindrfac = radius_rel_change_threshold_2;
            break;
          case 2:
            mindrfac = radius_rel_change_threshold_3;
            break;
          case 3:
            mindrfac = radius_rel_change_threshold_4;
            break;
          case 4:
            mindrfac = radius_rel_change_threshold_5;
            break;
          case 5:
            mindrfac = radius_rel_change_threshold_6;
            break;
          case 6:
            mindrfac = radius_rel_change_threshold_7;
            break;
          case 7:
            mindrfac = radius_rel_change_threshold_8;
            break;
          case 8:
            mindrfac = radius_rel_change_threshold_9;
            break;
          case 9:
            mindrfac = radius_rel_change_threshold_10;
            break;
          default:
            assert(0);
          }
          do_recompose = drfac > mindrfac;
          if (do_recompose)
            break;
        } // for rl
        if (do_recompose)
          break;

      } // for n
      if (verbose or veryverbose) {
        if (not do_recompose) {
          CCTK_INFO("Refined regions have not changed sufficiently; skipping "
                    "regridding");
        }
      }
    }
  }

  if (do_recompose) {
    *last_iteration = cctk_iteration;
    *last_map = Carpet::map;
  }

  if (do_recompose) {

    gh::rregs &superregss = *static_cast<gh::rregs *>(superregss_);
    gh::mregs &regsss = *static_cast<gh::mregs *>(regsss_);

    Regrid(cctkGH, superregss);

    // Make multiprocessor aware
    vector<vector<region_t> > regss(superregss.size());
    for (size_t rl = 0; rl < regss.size(); ++rl) {
      SplitRegions(cctkGH, superregss.at(rl), regss.at(rl));
    } // for rl

    // Make multigrid aware
    MakeMultigridBoxes(cctkGH, Carpet::map, regss, regsss);

    // Remember current positions
    for (int n = 0; n < num_centres; ++n) {
      old_active[n] = active[n];

      old_num_levels[n] = num_levels[n];

      old_position_x[n] = position_x[n];
      old_position_y[n] = position_y[n];
      old_position_z[n] = position_z[n];

      old_radius_x[n] = radius_x[n] < 0 ? radius[n] : radius_x[n];
      old_radius_y[n] = radius_y[n] < 0 ? radius[n] : radius_y[n];
      old_radius_z[n] = radius_z[n] < 0 ? radius[n] : radius_z[n];
    }

  } // if do_recompose

  timer.stop();

  return do_recompose;
}

CCTK_INT
CarpetRegrid2_RegridMaps(CCTK_POINTER_TO_CONST const cctkGH_,
                         CCTK_POINTER const superregsss_,
                         CCTK_POINTER const regssss_, CCTK_INT const force) {
  cGH const *const cctkGH = static_cast<cGH const *>(cctkGH_);

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  char const *const where = "CarpetRegrid2::RegridMaps";
  static Timers::Timer timer(where);
  timer.start();

  assert(is_level_mode());

  // Decide whether to change the grid hierarchy
  assert(*last_map == -1); // ensure this remains unused
  bool do_recompose;
  if (force) {
    do_recompose = true;
  } else {
    if (regrid_every == -1) {
      do_recompose = false;
    } else if (regrid_every == 0) {
      do_recompose = cctk_iteration == 0;
    } else {
      // Regrid at most once per iteration
      do_recompose =
          (cctk_iteration == 0 or
           (cctk_iteration > 0 and (cctk_iteration - 1) % regrid_every == 0 and
            cctk_iteration > *last_iteration));
    }
  }

  if (verbose or veryverbose) {
    if (do_recompose) {
      for (int n = 0; n < num_centres; ++n) {
        if (active[n]) {
          CCTK_VInfo(CCTK_THORNSTRING,
                     "Centre %d is at position [%g,%g,%g] with %d levels",
                     n + 1, static_cast<double>(position_x[n]),
                     static_cast<double>(position_y[n]),
                     static_cast<double>(position_z[n]),
                     static_cast<int>(num_levels[n]));
        } else {
          CCTK_VInfo(CCTK_THORNSTRING, "Centre %d is not active", n + 1);
        }
      }
    }
  }

  if (not force and do_recompose and *last_iteration != -1) {

    int lsh[2];
    getvectorindex2(cctkGH, "CarpetRegrid2::radii", lsh);

    // Regrid only if the regions have changed sufficiently
    do_recompose = false;
    if (adaptive_refinement)
      do_recompose = true;
    if (not do_recompose) {
      for (int n = 0; n < num_centres; ++n) {

        // When debugging, sneakily add a new level, but skip the
        // initial regrid, and the regrid before the first time step
        if (add_levels_automatically and cctk_iteration > 1) {
          num_levels[n] = min(int(num_levels[n] + 1), maxreflevels);
          CCTK_VInfo(CCTK_THORNSTRING,
                     "Increasing number of levels of centre %d to %d (it=%d)",
                     n + 1, static_cast<int>(num_levels[n]), cctk_iteration);
        }

        // Regrid if a region became active or inactive
        do_recompose = active[n] != old_active[n];
        if (do_recompose)
          break;

        // Check only active regions
        if (not active[n])
          continue;

        // Regrid if the number of levels changed
        do_recompose = num_levels[n] != old_num_levels[n];
        if (do_recompose)
          break;

        // Regrid if the positions have changed sufficiently
        CCTK_REAL const dist2 = pow(position_x[n] - old_position_x[n], 2) +
                                pow(position_y[n] - old_position_y[n], 2) +
                                pow(position_z[n] - old_position_z[n], 2);
        CCTK_REAL mindist;
        switch (n) {
        case 0:
          mindist = movement_threshold_1;
          break;
        case 1:
          mindist = movement_threshold_2;
          break;
        case 2:
          mindist = movement_threshold_3;
          break;
        case 3:
          mindist = movement_threshold_4;
          break;
        case 4:
          mindist = movement_threshold_5;
          break;
        case 5:
          mindist = movement_threshold_6;
          break;
        case 6:
          mindist = movement_threshold_7;
          break;
        case 7:
          mindist = movement_threshold_8;
          break;
        case 8:
          mindist = movement_threshold_9;
          break;
        case 9:
          mindist = movement_threshold_10;
          break;
        default:
          assert(0);
        }
        do_recompose = dist2 > pow(mindist, 2);
        if (do_recompose)
          break;

        // Regrid if the radii have changed sufficiently
        for (int rl = 1; rl < num_levels[n]; ++rl) {
          int const ind = index2(lsh, rl, n);
          CCTK_REAL const rx = radius_x[ind] < 0 ? radius[ind] : radius_x[ind];
          CCTK_REAL const ry = radius_y[ind] < 0 ? radius[ind] : radius_y[ind];
          CCTK_REAL const rz = radius_z[ind] < 0 ? radius[ind] : radius_z[ind];
          rvect const rad(rx, ry, rz);
          rvect const oldrad(old_radius_x[ind], old_radius_y[ind],
                             old_radius_z[ind]);
          CCTK_REAL const drfac =
              sqrt(sum(ipow(rad - oldrad, 2))) / sqrt(sum(ipow(oldrad, 2)));
          CCTK_REAL mindrfac;
          switch (n) {
          case 0:
            mindrfac = radius_rel_change_threshold_1;
            break;
          case 1:
            mindrfac = radius_rel_change_threshold_2;
            break;
          case 2:
            mindrfac = radius_rel_change_threshold_3;
            break;
          case 3:
            mindrfac = radius_rel_change_threshold_4;
            break;
          case 4:
            mindrfac = radius_rel_change_threshold_5;
            break;
          case 5:
            mindrfac = radius_rel_change_threshold_6;
            break;
          case 6:
            mindrfac = radius_rel_change_threshold_7;
            break;
          case 7:
            mindrfac = radius_rel_change_threshold_8;
            break;
          case 8:
            mindrfac = radius_rel_change_threshold_9;
            break;
          case 9:
            mindrfac = radius_rel_change_threshold_10;
            break;
          default:
            assert(0);
          }
          do_recompose = drfac > mindrfac;
          if (do_recompose)
            break;
        } // for rl
        if (do_recompose)
          break;

      } // for n
      if (verbose or veryverbose) {
        if (not do_recompose) {
          CCTK_INFO("Refined regions have not changed sufficiently; skipping "
                    "regridding");
        }
      }
    }
  }

  if (do_recompose) {
    *last_iteration = cctk_iteration;
  }

  if (do_recompose) {

    vector<gh::rregs> &superregsss =
        *static_cast<vector<gh::rregs> *>(superregsss_);
    vector<gh::mregs> &regssss = *static_cast<vector<gh::mregs> *>(regssss_);

    vector<gh::rregs> old_superregsss;
    vector<gh::mregs> old_regssss;
    old_superregsss.swap(superregsss);
    old_regssss.swap(regssss);

    superregsss = old_superregsss; // only needed for output
    BEGIN_MAP_LOOP(cctkGH, CCTK_GF) {
      Regrid(cctkGH, superregsss.at(Carpet::map));
    }
    END_MAP_LOOP;

    vector<vector<vector<region_t> > > regsss(maps);

    // Count levels
    vector<int> rls(maps);
    for (int m = 0; m < maps; ++m) {
      rls.at(m) = superregsss.at(m).size();
    }
    int maxrl = 0;
    for (int m = 0; m < maps; ++m) {
      maxrl = max(maxrl, rls.at(m));
    }
    // All maps must have the same number of levels
    for (int m = 0; m < maps; ++m) {
      superregsss.at(m).resize(maxrl);
      regsss.at(m).resize(maxrl);
    }

    // Make multiprocessor aware
    bool any_level_did_change = false;
    for (int rl = 0; rl < maxrl; ++rl) {
      bool level_did_change = false;
      for (int m = 0; m < maps; ++m) {
        level_did_change =
            level_did_change or int(old_superregsss.at(m).size()) <= rl or
            superregsss.at(m).at(rl) != old_superregsss.at(m).at(rl);
      }
      any_level_did_change = any_level_did_change or level_did_change;

      if (level_did_change) {
        // The level changed: perform domain decomposition

        vector<vector<region_t> > superregss(maps);
        for (int m = 0; m < maps; ++m) {
          superregss.at(m) = superregsss.at(m).at(rl);
        }
        vector<vector<region_t> > regss(maps);
        SplitRegionsMaps(cctkGH, superregss, regss);

        for (int m = 0; m < maps; ++m) {
          superregsss.at(m).at(rl) = superregss.at(m);
          regsss.at(m).at(rl) = regss.at(m);
        }

      } else {
        // The level did not actually change: re-use the old domain
        // decomposition

        for (int m = 0; m < maps; ++m) {
          superregsss.at(m).at(rl).swap(old_superregsss.at(m).at(rl));
          int const ml = 0;
          regsss.at(m).at(rl).swap(old_regssss.at(m).at(ml).at(rl));
        }

      } // if level did change
    }   // for rl

    // Make multigrid aware
    MakeMultigridBoxesMaps(cctkGH, regsss, regssss);

    int lsh[2];
    getvectorindex2(cctkGH, "CarpetRegrid2::radii", lsh);

    // Remember current positions
    for (int n = 0; n < num_centres; ++n) {
      old_active[n] = active[n];

      old_num_levels[n] = num_levels[n];

      old_position_x[n] = position_x[n];
      old_position_y[n] = position_y[n];
      old_position_z[n] = position_z[n];

      for (int rl = 1; rl < num_levels[n]; ++rl) {
        int const ind = index2(lsh, rl, n);
        old_radius_x[ind] = radius_x[ind] < 0 ? radius[ind] : radius_x[ind];
        old_radius_y[ind] = radius_y[ind] < 0 ? radius[ind] : radius_y[ind];
        old_radius_z[ind] = radius_z[ind] < 0 ? radius[ind] : radius_z[ind];
      }
    }

    do_recompose = do_recompose and any_level_did_change;

  } // if do_recompose

  timer.stop();

  return do_recompose;
}

} // namespace CarpetRegrid2
