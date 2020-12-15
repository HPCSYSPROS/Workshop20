#include <algorithm>
#include <cassert>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>

#include <cctk.h>
#include <cctk_Parameters.h>

#include <Timer.hh>

#include <cacheinfo.hh>
#include <defs.hh>
#include <gdata.hh>
#include <ggf.hh>

#include <carpet.hh>

namespace Carpet {

using namespace std;

// Local routine:
// Access a const int* as ivect reference (without the const)
static inline ivect &ivect_ref(const int *const iptr) {
  return ivect::ref(const_cast<int *>(iptr));
}

//
// Mode indicators
//

bool is_meta_mode() {
  if (mglevel == -1)
    assert(reflevel == -1 and map == -1 and component == -1);
  return mglevel == -1;
}

bool is_global_mode() {
  if (mglevel == -1 and reflevel != -1)
    assert(map == -1 and component == -1);
  return mglevel != -1 and reflevel == -1;
}

bool is_level_mode() {
  if (mglevel != -1 and reflevel != -1 and map == -1)
    assert(component == -1);
  return mglevel != -1 and reflevel != -1 and map == -1;
}

bool is_singlemap_mode() {
  return mglevel != -1 and reflevel != -1 and map != -1 and component == -1;
}

bool is_local_mode() {
  return mglevel != -1 and reflevel != -1 and map != -1 and component != -1;
  //       assert (mglevel>=0 and mglevel<mglevels);
  //       assert (reflevel>=0 and reflevel<reflevels);
  //       assert (map>=0 and map<maps);
  //       assert (vhh.AT(map)->local_components(reflevel)==1 or component==-1);
}

//
// Mode setting
//

// Set global mode

void enter_global_mode(cGH *const cctkGH, int const ml) {
  DECLARE_CCTK_PARAMETERS;

  assert(is_meta_mode());
  assert(ml >= 0 and ml < mglevels);
  Checkpoint("Entering global mode");

  static Timers::Timer ftimer("enter_global_mode", 1);
  ftimer.start();

  mglevel = ml;
  mglevelfact = ipow(mgfact, mglevel);
  // TODO: this could also just be "mglevel" instead
  cctkGH->cctk_convlevel = basemglevel + mglevel;

// Set mode
#ifdef HAVE_CGH_CCTK_MODE
  cctkGH->cctk_mode = CCTK_MODE_GLOBAL;
#endif

  // Set time delta
  if (not adaptive_stepsize) {
    cctkGH->cctk_delta_time = delta_time * mglevelfact;
  } else {
    delta_time = cctkGH->cctk_delta_time;
  }
  if (maps == 1) {
    // Set space delta
    for (int d = 0; d < dim; ++d) {
      cctkGH->cctk_origin_space[d] = origin_space.AT(0).AT(mglevel)[d];
      cctkGH->cctk_delta_space[d] = delta_space.AT(0)[d] * mglevelfact;
    }
  }

  // Set array information
  for (int group = 0; group < CCTK_NumGroups(); ++group) {
    cGroup gp;
    check(not CCTK_GroupData(group, &gp));
    if (gp.grouptype != CCTK_GF) {

      const int rl = 0;
      const int m = 0;
      const int c = dist::rank();

      const ibbox &baseext =
          arrdata.AT(group).AT(m).hh->baseextents.AT(ml).AT(rl);
      const ibbox &ext =
          arrdata.AT(group).AT(m).dd->light_boxes.AT(ml).AT(rl).AT(c).exterior;
      const b2vect &obnds = arrdata.AT(group).AT(m).hh->outer_boundaries(rl, c);

      cGroupDynamicData &info = groupdata.AT(group).info;

      ivect_ref(info.nghostzones) =
          arrdata.AT(group).AT(m).dd->ghost_widths.AT(0)[0];
      assert(all(ivect_ref(info.nghostzones) ==
                 arrdata.AT(group).AT(m).dd->ghost_widths.AT(0)[1]));
      ivect_ref(info.gsh) = baseext.shape() / baseext.stride();
      ivect_ref(info.lbnd) = (ext.lower() - baseext.lower()) / ext.stride();
      ivect_ref(info.ubnd) = (ext.upper() - baseext.lower()) / ext.stride();
      ivect_ref(info.lsh) = ext.shape() / ext.stride();
      ivect_ref(info.ash) = pad_shape(ext);
      // CCTK_DISTRIB_CONSTANT groups are handled by extending an extra
      // dimension to nprocs elements. Here we undo this so that the user
      // code does not see this trick.
      if (gp.disttype == CCTK_DISTRIB_CONSTANT) {
        int const dir = min(dim - 1, gp.dim == 0 ? 1 : gp.dim);
        ivect &gsh = ivect_ref(info.gsh);
        ivect &lbnd = ivect_ref(info.lbnd);
        ivect &ubnd = ivect_ref(info.ubnd);
        ivect const &lsh = ivect_ref(info.lsh);
        gsh[dir] = lsh[dir];
        lbnd[dir] = 0;
        ubnd[dir] = lsh[dir] - 1;
      }
      for (int d = 0; d < dim; ++d) {
        const_cast<int *>(info.bbox)[2 * d] = obnds[0][d];
        const_cast<int *>(info.bbox)[2 * d + 1] = obnds[1][d];
      }
      info.activetimelevels =
          groupdata.AT(group).activetimelevels.AT(mglevel).AT(0);

      for (int d = 0; d < dim; ++d) {
        assert(info.lsh[d] >= 0);
        assert(info.lsh[d] <= info.gsh[d]);
        assert(info.lbnd[d] >= 0);
        assert(info.lbnd[d] <= info.ubnd[d] + 1);
        assert(info.ubnd[d] < info.gsh[d]);
        assert(info.lbnd[d] + info.lsh[d] - 1 == info.ubnd[d]);
        assert(info.lbnd[d] <= info.ubnd[d] + 1);
        assert(info.lsh[d] <= info.ash[d]);
      }

      const int numvars = CCTK_NumVarsInGroupI(group);
      if (numvars > 0) {
        const int firstvar = CCTK_FirstVarIndexI(group);
        assert(firstvar >= 0);
        // TODO: modify flesh to allow changing max_tl
        const int active_tl = info.activetimelevels;
        assert(active_tl >= 0);

        assert(arrdata.AT(group).AT(m).hh->is_local(rl, c));

        for (int var = 0; var < numvars; ++var) {
          assert(firstvar + var < CCTK_NumVars());
          ggf *const ff = arrdata.AT(group).AT(m).data.AT(var);
          for (int tl = 0; tl < info.maxtimelevels; ++tl) {
            if (ff and tl < active_tl) {
              int const lc = 0;
              gdata *const data = ff->data_pointer(tl, rl, lc, ml);
              assert(data);
              cctkGH->data[firstvar + var][tl] = data->storage();
            } else {
              cctkGH->data[firstvar + var][tl] = NULL;
            }
          }
        }
      }

    } // if grouptype
  }   // for group

  ftimer.stop();

  Timers::Timer timer("global mode", 1);
  timer.start();

  assert(is_global_mode());
}

void leave_global_mode(cGH *const cctkGH) {
  DECLARE_CCTK_PARAMETERS;

  assert(is_global_mode() or is_meta_mode());

  if (mglevel == -1)
    return; // early return

  Checkpoint("Leaving global mode");

  Timers::Timer timer("global mode", 1);
  timer.stop();

  Timers::Timer ftimer("leave_global_mode", 1);
  ftimer.start();

  // Unset time delta
  if (not adaptive_stepsize) {
    cctkGH->cctk_delta_time = 0.0;
  }
  if (maps == 1) {
    // Save and unset space delta
    for (int d = 0; d < dim; ++d) {
      origin_space.AT(0).AT(mglevel)[d] = cctkGH->cctk_origin_space[d];
      delta_space.AT(mglevel)[d] = cctkGH->cctk_delta_space[d] / mglevelfact;
      cctkGH->cctk_origin_space[d] = -424242.0;
      cctkGH->cctk_delta_space[d] = -424242.0;
    }
  }

  CCTK_INT const deadbeef = get_deadbeef();

  // Set array information
  for (int group = 0; group < CCTK_NumGroups(); ++group) {
    if (CCTK_GroupTypeI(group) != CCTK_GF) {

      const int m = 0;

      cGroupDynamicData &info = groupdata.AT(group).info;

      // ivect_ref(info.nghostzones) = deadbeef;
      ivect_ref(info.nghostzones) =
          arrdata.AT(group).AT(m).dd->ghost_widths.AT(0)[0];
      ivect_ref(info.gsh) = deadbeef;
      ivect_ref(info.lbnd) = -deadbeef;
      ivect_ref(info.ubnd) = deadbeef;
      ivect_ref(info.lsh) = deadbeef;
      for (int d = 0; d < dim; ++d) {
        const_cast<int *>(info.bbox)[2 * d] = deadbeef;
        const_cast<int *>(info.bbox)[2 * d + 1] = deadbeef;
      }
      ivect_ref(info.ash) = deadbeef;
      info.activetimelevels = deadbeef;

      const int numvars = CCTK_NumVarsInGroupI(group);
      if (numvars > 0) {
        const int firstvar = CCTK_FirstVarIndexI(group);
        assert(firstvar >= 0);
        const int max_tl = groupdata.AT(group).info.maxtimelevels;
        assert(max_tl >= 0);

        assert(group < (int)arrdata.size());
        for (int var = 0; var < numvars; ++var) {
          assert(firstvar + var < CCTK_NumVars());
          for (int tl = 0; tl < max_tl; ++tl) {
            cctkGH->data[firstvar + var][tl] = NULL;
          }
        }
      }

    } // if grouptype
  }   // for group

// Set mode
#ifdef HAVE_CGH_CCTK_MODE
  cctkGH->cctk_mode = CCTK_MODE_META;
#endif

  mglevel = -1;
  mglevelfact = -deadbeef;
  cctkGH->cctk_convlevel = -deadbeef;

  ftimer.stop();

  assert(is_meta_mode());
}

// Set level mode

void enter_level_mode(cGH *const cctkGH, int const rl) {
  DECLARE_CCTK_PARAMETERS;

  assert(is_global_mode());
  assert(rl >= 0 and rl < reflevels);
  Checkpoint("Entering level mode");

  static Timers::Timer ftimer("enter_level_mode", 1);
  ftimer.start();

// Set mode
#ifdef HAVE_CGH_CCTK_MODE
  cctkGH->cctk_mode = CCTK_MODE_LEVEL;
#endif

  reflevel = rl;
  timereflevelfact = timereffacts.AT(reflevel);
  spacereflevelfact = spacereffacts.AT(reflevel);
  ivect_ref(cctkGH->cctk_levfac) = spacereflevelfact;
  cctkGH->cctk_timefac = timereflevelfact;

  // Set number of time levels
  for (int group = 0; group < CCTK_NumGroups(); ++group) {
    if (CCTK_GroupTypeI(group) == CCTK_GF) {
      cGroupDynamicData &info = groupdata.AT(group).info;
      info.activetimelevels =
          groupdata.AT(group).activetimelevels.AT(mglevel).AT(reflevel);
    }
  }

  // Set current time
  // assert (mglevel>=0 and mglevel<(int)leveltimes.size());
  // assert (reflevel>=0 and reflevel<(int)leveltimes.AT(mglevel).size());
  // if (not adaptive_stepsize) {
  //   cctkGH->cctk_time = leveltimes.AT(mglevel).AT(reflevel);
  // } else {
  //   leveltimes.AT(mglevel).AT(reflevel) = cctkGH->cctk_time;
  // }
  if (not adaptive_stepsize) {
    cctkGH->cctk_time = tt->get_time(mglevel, reflevel, timelevel);
  } else {
    tt->set_time(mglevel, reflevel, timelevel, cctkGH->cctk_time);
  }

  ftimer.stop();

  ostringstream mode_s;
  mode_s << "level(" << reflevel << ")";
  Timers::Timer timer(mode_s.str(), 1);
  timer.start();

  assert(is_level_mode());
}

void leave_level_mode(cGH *const cctkGH) {
  DECLARE_CCTK_PARAMETERS;

  assert(is_level_mode() or is_global_mode());

  if (reflevel == -1)
    return; // early return

  Checkpoint("Leaving level mode");

  ostringstream mode_s;
  mode_s << "level(" << reflevel << ")";
  Timers::Timer timer(mode_s.str(), 1);
  timer.stop();

  static Timers::Timer ftimer("leave_level_mode", 1);
  ftimer.start();

  CCTK_INT const deadbeef = get_deadbeef();

  // Save and unset current time
  // assert (mglevel>=0 and mglevel<(int)leveltimes.size());
  // assert (reflevel>=0 and reflevel<(int)leveltimes.AT(mglevel).size());
  // leveltimes.AT(mglevel).AT(reflevel) = cctkGH->cctk_time;
  // if (not adaptive_stepsize) {
  //   cctkGH->cctk_time = global_time;
  // } else {
  //   global_time = cctkGH->cctk_time;
  // }
  tt->set_time(mglevel, reflevel, timelevel, cctkGH->cctk_time);
  if (not adaptive_stepsize) {
    cctkGH->cctk_time = global_time;
  } else {
    global_time = cctkGH->cctk_time;
  }

  // Unset number of time levels
  for (int group = 0; group < CCTK_NumGroups(); ++group) {
    if (CCTK_GroupTypeI(group) == CCTK_GF) {
      cGroupDynamicData &info = groupdata.AT(group).info;
      info.activetimelevels = deadbeef;
    }
  }

// Set mode
#ifdef HAVE_CGH_CCTK_MODE
  cctkGH->cctk_mode = CCTK_MODE_GLOBAL;
#endif

  reflevel = -1;
  timereflevelfact = timereffacts.AT(reflevels - 1);
  // TODO: use spacereffacts.AT (reflevel - 1) instead?
  spacereflevelfact = ivect(-deadbeef);
  ivect_ref(cctkGH->cctk_levfac) = spacereflevelfact;
  cctkGH->cctk_timefac = timereflevelfact;

  ftimer.stop();

  assert(is_global_mode());
}

// Set singlemap mode

void enter_singlemap_mode(cGH *const cctkGH, int const m, int const grouptype) {
  DECLARE_CCTK_PARAMETERS

  assert(is_level_mode());
  assert(mc_grouptype == -1);
  assert(m >= 0 and m < maps);
  assert(grouptype == CCTK_SCALAR or grouptype == CCTK_ARRAY or
         grouptype == CCTK_GF);
  Checkpoint("Entering singlemap mode");

  static Timers::Timer ftimer("enter_singlemap_mode", 1);
  ftimer.start();

// Set mode
#ifdef HAVE_CGH_CCTK_MODE
  cctkGH->cctk_mode = CCTK_MODE_SINGLEMAP;
#endif

  assert(mc_grouptype == -1);
  mc_grouptype = grouptype;
  carpetGH.map = map = m;

  if (mc_grouptype == CCTK_GF) {

    if (maps > 1) {
      // Set space delta
      for (int d = 0; d < dim; ++d) {
        cctkGH->cctk_origin_space[d] = origin_space.AT(map).AT(mglevel)[d];
        cctkGH->cctk_delta_space[d] = delta_space.AT(map)[d] * mglevelfact;
      }
    }

    // Set grid shape
    const ibbox &coarseext = vhh.AT(map)->baseextents.AT(mglevel).AT(0);
    const ibbox &baseext = vhh.AT(map)->baseextents.AT(mglevel).AT(reflevel);
    //       assert (all (baseext.lower() % baseext.stride() == 0));

    ivect_ref(cctkGH->cctk_levoff) = baseext.lower() - coarseext.lower();
    ivect_ref(cctkGH->cctk_levoffdenom) = baseext.stride();
    ivect_ref(cctkGH->cctk_gsh) = baseext.shape() / baseext.stride();
    assert(all(vdd.AT(map)->ghost_widths.AT(reflevel)[0] ==
               vdd.AT(map)->ghost_widths.AT(reflevel)[1]));
    ivect_ref(cctkGH->cctk_nghostzones) =
        vdd.AT(map)->ghost_widths.AT(reflevel)[0];

    for (int group = 0; group < CCTK_NumGroups(); ++group) {
      if (CCTK_GroupTypeI(group) == CCTK_GF) {
        cGroupDynamicData &info = groupdata.AT(group).info;
        ivect_ref(info.gsh) = ivect_ref(cctkGH->cctk_gsh);
        ivect_ref(info.nghostzones) = ivect_ref(cctkGH->cctk_nghostzones);
      }
    }

  } // if mc_grouptype

  ftimer.stop();

  if (include_maps_in_mode_timer_tree) {
    ostringstream mode_s;
    mode_s << "map(" << map << ")";
    Timers::Timer timer(mode_s.str(), 1);
    timer.start();
  }

  assert(is_singlemap_mode());
}

void leave_singlemap_mode(cGH *const cctkGH) {
  DECLARE_CCTK_PARAMETERS;

  assert(is_singlemap_mode() or is_level_mode());

  if (map == -1)
    return; // early return

  Checkpoint("Leaving singlemap mode");

  if (include_maps_in_mode_timer_tree) {
    ostringstream mode_s;
    mode_s << "map(" << map << ")";
    Timers::Timer timer(mode_s.str(), 1);
    timer.stop();
  }

  static Timers::Timer ftimer("leave_singlemap_mode", 1);
  ftimer.start();

  assert(mc_grouptype == CCTK_SCALAR or mc_grouptype == CCTK_ARRAY or
         mc_grouptype == CCTK_GF);

  if (mc_grouptype == CCTK_GF) {

    CCTK_INT const deadbeef = get_deadbeef();

    // Save space delta
    // (Do this early and often, so that interpolation has access to
    // the correct values right away.)
    for (int d = 0; d < dim; ++d) {
      origin_space.AT(map).AT(mglevel)[d] = cctkGH->cctk_origin_space[d];
      delta_space.AT(map)[d] = cctkGH->cctk_delta_space[d] / mglevelfact;
    }
    if (maps > 1) {
      // Unset space delta
      for (int d = 0; d < dim; ++d) {
        cctkGH->cctk_origin_space[d] = -424242.0;
        cctkGH->cctk_delta_space[d] = -424242.0;
      }
    }

    // Unset grid shape
    ivect_ref(cctkGH->cctk_levoff) = deadbeef;
    ivect_ref(cctkGH->cctk_levoffdenom) = 0;
    ivect_ref(cctkGH->cctk_gsh) = deadbeef;
    ivect_ref(cctkGH->cctk_nghostzones) = deadbeef;

    for (int group = 0; group < CCTK_NumGroups(); ++group) {
      if (CCTK_GroupTypeI(group) == CCTK_GF) {
        cGroupDynamicData &info = groupdata.AT(group).info;
        ivect_ref(info.gsh) = ivect_ref(cctkGH->cctk_gsh);
        ivect_ref(info.nghostzones) = ivect_ref(cctkGH->cctk_nghostzones);
      }
    }

  } // if mc_grouptype

// Set mode
#ifdef HAVE_CGH_CCTK_MODE
  cctkGH->cctk_mode = CCTK_MODE_LEVEL;
#endif

  mc_grouptype = -1;
  carpetGH.map = map = -1;

  ftimer.stop();

  assert(is_level_mode());
}

// Set local mode

void enter_local_mode(cGH *const cctkGH, int const c, int const lc,
                      int const grouptype) {
  DECLARE_CCTK_PARAMETERS

  assert(is_singlemap_mode());
  if (mc_grouptype == CCTK_GF) {
    assert(c >= 0 and c < vhh.AT(map)->components(reflevel));
    assert(lc == -1 or
           (lc >= 0 and lc < vhh.AT(map)->local_components(reflevel)));
  } else {
    assert(c >= 0 and c < dist::size());
    assert(lc == -1 or lc == 0);
  }
  Checkpoint("Entering local mode");

  static Timers::Timer ftimer("enter_local_mode", 1);
  ftimer.start();

// Set mode
#ifdef HAVE_CGH_CCTK_MODE
  cctkGH->cctk_mode = CCTK_MODE_LOCAL;
#endif

  assert(grouptype == mc_grouptype);
  component = c;
  local_component = lc;

  if (mc_grouptype == CCTK_GF) {

    // Set cGH fields
    const ibbox &baseext = vhh.AT(map)->baseextents.AT(mglevel).AT(reflevel);
    const ibbox &ext = vdd.AT(map)
                           ->light_boxes.AT(mglevel)
                           .AT(reflevel)
                           .AT(component)
                           .exterior;
    const b2vect &obnds = vhh.AT(map)->outer_boundaries(reflevel, component);

    ivect_ref(cctkGH->cctk_lbnd) =
        (ext.lower() - baseext.lower()) / ext.stride();
    ivect_ref(cctkGH->cctk_ubnd) =
        (ext.upper() - baseext.lower()) / ext.stride();
    ivect_ref(cctkGH->cctk_lsh) = ext.shape() / ext.stride();
    ivect_ref(cctkGH->cctk_ash) = pad_shape(ext);

    for (int d = 0; d < dim; ++d) {
      cctkGH->cctk_bbox[2 * d] = obnds[0][d];
      cctkGH->cctk_bbox[2 * d + 1] = obnds[1][d];
    }

    for (int d = 0; d < dim; ++d) {
      cctkGH->cctk_from[d] = 0;
      cctkGH->cctk_to[d] = cctkGH->cctk_lsh[d];
    }

    for (int d = 0; d < dim; ++d) {
      assert(cctkGH->cctk_lsh[d] >= 0);
      assert(cctkGH->cctk_lsh[d] <= cctkGH->cctk_gsh[d]);
      assert(cctkGH->cctk_lbnd[d] >= 0);
      assert(cctkGH->cctk_lbnd[d] <= cctkGH->cctk_ubnd[d] + 1);
      assert(cctkGH->cctk_ubnd[d] < cctkGH->cctk_gsh[d]);
      assert(cctkGH->cctk_lbnd[d] + cctkGH->cctk_lsh[d] - 1 ==
             cctkGH->cctk_ubnd[d]);
      assert(cctkGH->cctk_lbnd[d] <= cctkGH->cctk_ubnd[d] + 1);
      assert(cctkGH->cctk_lsh[d] <= cctkGH->cctk_ash[d]);
      assert(cctkGH->cctk_from[d] >= 0);
      assert(cctkGH->cctk_from[d] <= cctkGH->cctk_to[d]);
      assert(cctkGH->cctk_to[d] <= cctkGH->cctk_lsh[d]);
    }

    for (int group = 0; group < CCTK_NumGroups(); ++group) {
      if (CCTK_GroupTypeI(group) == CCTK_GF) {

        cGroupDynamicData &info = groupdata.AT(group).info;

        ivect_ref(info.lbnd) = ivect_ref(cctkGH->cctk_lbnd);
        ivect_ref(info.ubnd) = ivect_ref(cctkGH->cctk_ubnd);
        ivect_ref(info.lsh) = ivect_ref(cctkGH->cctk_lsh);
        ivect_ref(info.ash) = ivect_ref(cctkGH->cctk_ash);

        for (int d = 0; d < dim; ++d) {
          const_cast<int *>(info.bbox)[2 * d] = cctkGH->cctk_bbox[2 * d];
          const_cast<int *>(info.bbox)[2 * d + 1] =
              cctkGH->cctk_bbox[2 * d + 1];
        }

        if (local_component != -1) {
          const int numvars = CCTK_NumVarsInGroupI(group);
          if (numvars > 0) {
            const int firstvar = CCTK_FirstVarIndexI(group);
            assert(firstvar >= 0);
            // TODO: Modify flesh so the max_tl can be updated
            const int active_tl =
                groupdata.AT(group).activetimelevels.AT(mglevel).AT(reflevel);
            assert(active_tl >= 0);
            int available_tl;
            int tl_offset;
            if (active_tl == 0) {
              // group has no storage
              available_tl = active_tl;
              tl_offset = 0;
            } else if (do_allow_past_timelevels) {
              // regular case; all timelevels are accessible
              available_tl = active_tl;
              tl_offset = 0;
            } else {
              // only one timelevel is accessible
              available_tl = 1;
              if (timelevel < active_tl) {
                // timelevel "timelevel" exists
                tl_offset = timelevel;
              } else {
                // timelevel "timelevel" does not exist
                tl_offset = active_tl - 1;
              }
            }

            // assert (vhh.AT(map)->is_local(reflevel,component));

            assert(group < (int)arrdata.size());
            for (int var = 0; var < numvars; ++var) {
              assert(firstvar + var < CCTK_NumVars());
              ggf *const ff = arrdata.AT(group).AT(map).data.AT(var);
              for (int tl = 0; tl < info.maxtimelevels; ++tl) {
                if (ff and tl < available_tl) {
                  gdata *const data = ff->data_pointer(
                      tl_offset + tl, reflevel, local_component, mglevel);
                  assert(data);
                  cctkGH->data[firstvar + var][tl] = data->storage();
                } else {
                  cctkGH->data[firstvar + var][tl] = NULL;
                }
              }
            }
          }
        }

      } // if grouptype
    }   // for group

  } // if mc_grouptype

  ftimer.stop();

  if (include_local_mode_in_mode_timer_tree) {
    Timers::Timer timer("local", 1);
    timer.start();
  }

  assert(is_local_mode());
}

void leave_local_mode(cGH *const cctkGH) {
  DECLARE_CCTK_PARAMETERS;

  assert(is_local_mode() or is_singlemap_mode());

  if (component == -1)
    return; // early return

  Checkpoint("Leaving local mode");

  if (include_local_mode_in_mode_timer_tree) {
    Timers::Timer timer("local", 1);
    timer.stop();
  }

  static Timers::Timer ftimer("leave_local_mode", 1);
  ftimer.start();

  if (mc_grouptype == CCTK_GF) {

    CCTK_INT const deadbeef = get_deadbeef();

    // Unset cGH fields
    ivect_ref(cctkGH->cctk_lbnd) = -deadbeef;
    ivect_ref(cctkGH->cctk_ubnd) = deadbeef;
    ivect_ref(cctkGH->cctk_from) = -deadbeef;
    ivect_ref(cctkGH->cctk_to) = deadbeef;
    ivect_ref(cctkGH->cctk_lsh) = deadbeef;
    ivect_ref(cctkGH->cctk_ash) = deadbeef;

    for (int d = 0; d < dim; ++d) {
      cctkGH->cctk_bbox[2 * d] = deadbeef;
      cctkGH->cctk_bbox[2 * d + 1] = deadbeef;
    }

    for (int group = 0; group < CCTK_NumGroups(); ++group) {
      if (CCTK_GroupTypeI(group) == CCTK_GF) {

        cGroupDynamicData &info = groupdata.AT(group).info;

        ivect_ref(info.lbnd) = ivect_ref(cctkGH->cctk_lbnd);
        ivect_ref(info.ubnd) = ivect_ref(cctkGH->cctk_ubnd);
        ivect_ref(info.lsh) = ivect_ref(cctkGH->cctk_lsh);
        ivect_ref(info.ash) = ivect_ref(cctkGH->cctk_ash);

        for (int d = 0; d < dim; ++d) {
          const_cast<int *>(info.bbox)[2 * d] = cctkGH->cctk_bbox[2 * d];
          const_cast<int *>(info.bbox)[2 * d + 1] =
              cctkGH->cctk_bbox[2 * d + 1];
        }

        if (local_component != -1) {
          const int numvars = CCTK_NumVarsInGroupI(group);
          if (numvars > 0) {
            const int firstvar = CCTK_FirstVarIndexI(group);
            assert(firstvar >= 0);
            const int active_tl =
                groupdata.AT(group).activetimelevels.AT(mglevel).AT(reflevel);
            assert(active_tl >= 0);

            assert(group < (int)arrdata.size());
            for (int var = 0; var < numvars; ++var) {
              assert(firstvar + var < CCTK_NumVars());
              for (int tl = 0; tl < active_tl; ++tl) {
                cctkGH->data[firstvar + var][tl] = NULL;
              }
            }
          }
        }

      } // if grouptype
    }   // for group

  } // if mc_grouptype

// Set mode
#ifdef HAVE_CGH_CCTK_MODE
  cctkGH->cctk_mode = CCTK_MODE_SINGLEMAP;
#endif

  component = -1;
  local_component = -1;

  ftimer.stop();

  assert(is_singlemap_mode());
}

//
// Mode iterators
//

// mglevel iterator

mglevel_iterator::mglevel_iterator(cGH const *const cctkGH_)
    : cctkGH(const_cast<cGH *>(cctkGH_)), ml(mglevels - 1) {
  enter_global_mode(cctkGH, ml);
}

mglevel_iterator::~mglevel_iterator() { leave_global_mode(cctkGH); }

bool mglevel_iterator::done() const { return ml < 0; }

void mglevel_iterator::step() {
  --ml;
  if (not done()) {
    leave_global_mode(cctkGH);
    enter_global_mode(cctkGH, ml);
  }
}

// reverse mglevel iterator

reverse_mglevel_iterator::reverse_mglevel_iterator(cGH const *const cctkGH_)
    : cctkGH(const_cast<cGH *>(cctkGH_)), ml(0) {
  enter_global_mode(cctkGH, ml);
}

reverse_mglevel_iterator::~reverse_mglevel_iterator() {
  leave_global_mode(cctkGH);
}

bool reverse_mglevel_iterator::done() const { return ml >= mglevels; }

void reverse_mglevel_iterator::step() {
  ++ml;
  if (not done()) {
    leave_global_mode(cctkGH);
    enter_global_mode(cctkGH, ml);
  }
}

// reflevel iterator

reflevel_iterator::reflevel_iterator(cGH const *const cctkGH_)
    : cctkGH(const_cast<cGH *>(cctkGH_)), rl(0) {
  enter_level_mode(cctkGH, rl);
}

reflevel_iterator::~reflevel_iterator() { leave_level_mode(cctkGH); }

bool reflevel_iterator::done() const { return rl >= reflevels; }

void reflevel_iterator::step() {
  ++rl;
  if (not done()) {
    leave_level_mode(cctkGH);
    enter_level_mode(cctkGH, rl);
  }
}

// reverse reflevel iterator

reverse_reflevel_iterator::reverse_reflevel_iterator(cGH const *const cctkGH_)
    : cctkGH(const_cast<cGH *>(cctkGH_)), rl(reflevels - 1) {
  enter_level_mode(cctkGH, rl);
}

reverse_reflevel_iterator::~reverse_reflevel_iterator() {
  leave_level_mode(cctkGH);
}

bool reverse_reflevel_iterator::done() const { return rl < 0; }

void reverse_reflevel_iterator::step() {
  --rl;
  if (not done()) {
    leave_level_mode(cctkGH);
    enter_level_mode(cctkGH, rl);
  }
}

// map iterator

map_iterator::map_iterator(cGH const *const cctkGH_, int const grouptype_)
    : cctkGH(const_cast<cGH *>(cctkGH_)), grouptype(grouptype_), m(0) {
  assert(grouptype == CCTK_GF or grouptype == CCTK_ARRAY or
         grouptype == CCTK_SCALAR);
  enter_singlemap_mode(cctkGH, m, grouptype);
}

map_iterator::~map_iterator() { leave_singlemap_mode(cctkGH); }

bool map_iterator::done() const {
  int const maxm = grouptype == CCTK_GF ? maps : 1;
  return m >= maxm;
}

void map_iterator::step() {
  ++m;
  if (not done()) {
    leave_singlemap_mode(cctkGH);
    enter_singlemap_mode(cctkGH, m, grouptype);
  }
}

// processor local map iterator

local_map_iterator::local_map_iterator(cGH const *const cctkGH_,
                                       int const grouptype_)
    : cctkGH(const_cast<cGH *>(cctkGH_)), grouptype(grouptype_), m(0) {
  assert(grouptype == CCTK_GF or grouptype == CCTK_ARRAY or
         grouptype == CCTK_SCALAR);
  enter_singlemap_mode(cctkGH, m, grouptype);
}

local_map_iterator::~local_map_iterator() { leave_singlemap_mode(cctkGH); }

bool local_map_iterator::done() const {
  int const maxm = grouptype == CCTK_GF ? maps : 1;
  return m >= maxm;
}

void local_map_iterator::step() {
  while (true) {
    ++m;
    if (done())
      break;
    int const maxlc =
        grouptype == CCTK_GF ? vhh.AT(m)->local_components(reflevel) : 1;
    if (maxlc > 0)
      break;
  }
  if (not done()) {
    leave_singlemap_mode(cctkGH);
    enter_singlemap_mode(cctkGH, m, grouptype);
  }
}

// local mode iterator

component_iterator::component_iterator(cGH const *const cctkGH_,
                                       int const grouptype_)
    : cctkGH(const_cast<cGH *>(cctkGH_)), grouptype(grouptype_), c(-1) {
  assert(grouptype == CCTK_GF or grouptype == CCTK_ARRAY or
         grouptype == CCTK_SCALAR);
  step();
}

component_iterator::~component_iterator() { leave_local_mode(cctkGH); }

bool component_iterator::done() const {
  int const maxc =
      (grouptype == CCTK_GF ? vhh.AT(map)->components(reflevel) : dist::size());
  return c >= maxc;
}

void component_iterator::step() {
  ++c;
  if (not done()) {
    leave_local_mode(cctkGH);
    int const lc =
        (grouptype == CCTK_GF ? vhh.AT(map)->get_local_component(reflevel, c)
                              : (c == dist::rank() ? 0 : -1));
    enter_local_mode(cctkGH, c, lc, grouptype);
  }
}

// processor local mode iterator

local_component_iterator::local_component_iterator(cGH const *const cctkGH_,
                                                   int const grouptype_)
    : cctkGH(const_cast<cGH *>(cctkGH_)), grouptype(grouptype_), lc(-1) {
  assert(grouptype == CCTK_GF or grouptype == CCTK_ARRAY or
         grouptype == CCTK_SCALAR);
  assert(is_singlemap_mode());
  step();
}

local_component_iterator::~local_component_iterator() {
  leave_local_mode(cctkGH);
}

bool local_component_iterator::done() const {
  int const maxlc =
      (grouptype == CCTK_GF ? vhh.AT(map)->local_components(reflevel) : 1);
  return lc >= maxlc;
}

void local_component_iterator::step() {
  ++lc;
  if (not done()) {
    int const c =
        (grouptype == CCTK_GF ? vhh.AT(map)->get_component(reflevel, lc)
                              : dist::rank());
    leave_local_mode(cctkGH);
    enter_local_mode(cctkGH, c, lc, grouptype);
  }
}

//
// Mode escapes
//

// Singlemap escape

singlemap_escape::singlemap_escape(cGH const *const cctkGH_)
    : cctkGH(const_cast<cGH *>(cctkGH_)), c(component), lc(local_component) {
  assert(not is_meta_mode());
  assert(not is_global_mode());
  assert(not is_level_mode());
  if (not is_singlemap_mode()) {
    leave_local_mode(cctkGH);
  }
}

singlemap_escape::~singlemap_escape() {
  assert(is_singlemap_mode());
  if (c != -1) {
    enter_local_mode(cctkGH, c, lc, mc_grouptype);
  }
}

// Level escape

level_escape::level_escape(cGH const *const cctkGH_)
    : cctkGH(const_cast<cGH *>(cctkGH_)), grouptype(mc_grouptype), m(map),
      c(component), lc(local_component) {
  assert(not is_meta_mode());
  assert(not is_global_mode());
  if (not is_level_mode()) {
    if (not is_singlemap_mode()) {
      leave_local_mode(cctkGH);
    }
    leave_singlemap_mode(cctkGH);
  }
}

level_escape::~level_escape() {
  assert(is_level_mode());
  if (m != -1) {
    enter_singlemap_mode(cctkGH, m, grouptype);
    if (c != -1) {
      enter_local_mode(cctkGH, c, lc, grouptype);
    }
  }
}

// Global escape

global_escape::global_escape(cGH const *const cctkGH_)
    : cctkGH(const_cast<cGH *>(cctkGH_)), rl(reflevel), grouptype(mc_grouptype),
      m(map), c(component), lc(local_component) {
  assert(not is_meta_mode());
  if (not is_global_mode()) {
    if (not is_level_mode()) {
      if (not is_singlemap_mode()) {
        leave_local_mode(cctkGH);
      }
      leave_singlemap_mode(cctkGH);
    }
    leave_level_mode(cctkGH);
  }
}

global_escape::~global_escape() {
  assert(is_global_mode());
  if (rl != -1) {
    enter_level_mode(cctkGH, rl);
    if (m != -1) {
      enter_singlemap_mode(cctkGH, m, grouptype);
      if (c != -1) {
        enter_local_mode(cctkGH, c, lc, grouptype);
      }
    }
  }
}

// Meta escape

meta_escape::meta_escape(cGH const *const cctkGH_)
    : cctkGH(const_cast<cGH *>(cctkGH_)), ml(mglevel), rl(reflevel),
      grouptype(mc_grouptype), m(map), c(component), lc(local_component) {
  if (not is_meta_mode()) {
    if (not is_global_mode()) {
      if (not is_level_mode()) {
        if (not is_singlemap_mode()) {
          leave_local_mode(cctkGH);
        }
        leave_singlemap_mode(cctkGH);
      }
      leave_level_mode(cctkGH);
    }
    leave_global_mode(cctkGH);
  }
}

meta_escape::~meta_escape() {
  assert(is_meta_mode());
  if (ml != -1) {
    enter_global_mode(cctkGH, ml);
    if (rl != -1) {
      enter_level_mode(cctkGH, rl);
      if (m != -1) {
        enter_singlemap_mode(cctkGH, m, grouptype);
        if (c != -1) {
          enter_local_mode(cctkGH, c, lc, grouptype);
        }
      }
    }
  }
}

//
// Call functions in specific modes
//

int CallLocalFunction(cGH *const cctkGH,
                      void (*const function)(cGH *const cctkGH)) {
  if (is_meta_mode()) {
    BEGIN_MGLEVEL_LOOP(cctkGH) {
      BEGIN_REFLEVEL_LOOP(cctkGH) {
        BEGIN_LOCAL_MAP_LOOP(cctkGH, CCTK_GF) {
          BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, CCTK_GF) { function(cctkGH); }
          END_LOCAL_COMPONENT_LOOP;
        }
        END_LOCAL_MAP_LOOP;
      }
      END_REFLEVEL_LOOP;
    }
    END_MGLEVEL_LOOP;
  } else if (is_global_mode()) {
    BEGIN_REFLEVEL_LOOP(cctkGH) {
      BEGIN_LOCAL_MAP_LOOP(cctkGH, CCTK_GF) {
        BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, CCTK_GF) { function(cctkGH); }
        END_LOCAL_COMPONENT_LOOP;
      }
      END_LOCAL_MAP_LOOP;
    }
    END_REFLEVEL_LOOP;
  } else if (is_level_mode()) {
    BEGIN_LOCAL_MAP_LOOP(cctkGH, CCTK_GF) {
      BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, CCTK_GF) { function(cctkGH); }
      END_LOCAL_COMPONENT_LOOP;
    }
    END_LOCAL_MAP_LOOP;
  } else if (is_singlemap_mode()) {
    BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, CCTK_GF) { function(cctkGH); }
    END_LOCAL_COMPONENT_LOOP;
  } else {
    function(cctkGH);
  }
  return 0;
}

int CallSinglemapFunction(cGH *const cctkGH,
                          void (*const function)(cGH *const cctkGH)) {
  if (is_meta_mode()) {
    BEGIN_MGLEVEL_LOOP(cctkGH) {
      BEGIN_REFLEVEL_LOOP(cctkGH) {
        BEGIN_MAP_LOOP(cctkGH, CCTK_GF) { function(cctkGH); }
        END_MAP_LOOP;
      }
      END_REFLEVEL_LOOP;
    }
    END_MGLEVEL_LOOP;
  } else if (is_global_mode()) {
    BEGIN_REFLEVEL_LOOP(cctkGH) {
      BEGIN_MAP_LOOP(cctkGH, CCTK_GF) { function(cctkGH); }
      END_MAP_LOOP;
    }
    END_REFLEVEL_LOOP;
  } else if (is_level_mode()) {
    BEGIN_MAP_LOOP(cctkGH, CCTK_GF) { function(cctkGH); }
    END_MAP_LOOP;
  } else {
    BEGIN_SINGLEMAP_MODE(cctkGH) { function(cctkGH); }
    END_SINGLEMAP_MODE;
  }
  return 0;
}

int CallLevelFunction(cGH *const cctkGH,
                      void (*const function)(cGH *const cctkGH)) {
  if (is_meta_mode()) {
    BEGIN_MGLEVEL_LOOP(cctkGH) {
      BEGIN_REFLEVEL_LOOP(cctkGH) { function(cctkGH); }
      END_REFLEVEL_LOOP;
    }
    END_MGLEVEL_LOOP;
  } else if (is_global_mode()) {
    BEGIN_REFLEVEL_LOOP(cctkGH) { function(cctkGH); }
    END_REFLEVEL_LOOP;
  } else {
    BEGIN_LEVEL_MODE(cctkGH) { function(cctkGH); }
    END_LEVEL_MODE;
  }
  return 0;
}

int CallGlobalFunction(cGH *const cctkGH,
                       void (*const function)(cGH *const cctkGH)) {
  if (is_meta_mode()) {
    BEGIN_MGLEVEL_LOOP(cctkGH) { function(cctkGH); }
    END_MGLEVEL_LOOP;
  } else {
    BEGIN_GLOBAL_MODE(cctkGH) { function(cctkGH); }
    END_GLOBAL_MODE;
  }
  return 0;
}

int CallMetaFunction(cGH *const cctkGH,
                     void (*const function)(cGH *const cctkGH)) {
  BEGIN_META_MODE(cctkGH) { function(cctkGH); }
  END_META_MODE;
  return 0;
}

//
// Mode setters
//

// mglevel setter

mglevel_setter::mglevel_setter(cGH const *const cctkGH_, int const ml)
    : cctkGH(const_cast<cGH *>(cctkGH_)) {
  assert(is_meta_mode());
  enter_global_mode(cctkGH, ml);
}

mglevel_setter::~mglevel_setter() { leave_global_mode(cctkGH); }

// reflevel setter

reflevel_setter::reflevel_setter(cGH const *const cctkGH_, int const rl)
    : cctkGH(const_cast<cGH *>(cctkGH_)) {
  assert(is_global_mode());
  enter_level_mode(cctkGH, rl);
}

reflevel_setter::~reflevel_setter() { leave_level_mode(cctkGH); }

// map setter

map_setter::map_setter(cGH const *const cctkGH_, int const m,
                       int const grouptype)
    : cctkGH(const_cast<cGH *>(cctkGH_)) {
  assert(is_level_mode());
  enter_singlemap_mode(cctkGH, m, grouptype);
}

map_setter::~map_setter() { leave_singlemap_mode(cctkGH); }

// component setter

component_setter::component_setter(cGH const *const cctkGH_, int const c,
                                   int const grouptype)
    : cctkGH(const_cast<cGH *>(cctkGH_)) {
  assert(is_singlemap_mode());
  int const lc = (c != -1 ? (grouptype == CCTK_GF
                                 ? vhh.AT(map)->get_local_component(reflevel, c)
                                 : (c == dist::rank() ? 0 : -1))
                          : -1);
  enter_local_mode(cctkGH, c, lc, grouptype);
}

component_setter::~component_setter() { leave_local_mode(cctkGH); }

//
// Loop over a bounding box
//

void ibbox2iminimax(ibbox const &ext, // component extent
                    ibbox const &box, // this bbox
                    ivect &imin, ivect &imax) {
  assert(all((box.lower() - ext.lower()) >= 0));
  assert(all((box.upper() - ext.lower() + ext.stride()) >= 0));
  assert(all((box.lower() - ext.lower()) % ext.stride() == 0));
  assert(all((box.upper() - ext.lower() + ext.stride()) % ext.stride() == 0));

  imin = (box.lower() - ext.lower()) / ext.stride();
  imax = (box.upper() - ext.lower() + ext.stride()) / ext.stride();

  assert(all(0 <= imin));
  assert(box.empty() xor all(imin <= imax));
}

//
// Call a scheduling group
//

int CallScheduleGroup(cGH *const cctkGH, const char *const group) {
  CCTK_ScheduleTraverse(group, cctkGH, CallFunction);
  return 0;
}

extern "C" void CCTK_FCALL CCTK_FNAME(CallScheduleGroup)(
    int *const ierr, cGH *const *const cctkGH, ONE_FORTSTRING_ARG) {
  ONE_FORTSTRING_CREATE(group);
  *ierr = CallScheduleGroup(*cctkGH, group);
  free(group);
}

} // namespace Carpet
