#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <map>
#include <string>
#include <sstream>

#include <cctk.h>
#include <cctk_Parameters.h>
#include <cctki_GHExtensions.h>

#include <Requirements.hh>

#include <Timer.hh>

#include <gh.hh>

#include <carpet.hh>

#include "adler32.hh"

namespace Carpet {

using namespace std;

static void CallScheduledFunction(char const *restrict time_and_mode,
                                  void *function, cFunctionData *attribute,
                                  void *data, Timers::Timer &user_timer);

static void SyncGroupsInScheduleBlock(cFunctionData *attribute, cGH *cctkGH,
                                      vector<int> const &sync_groups,
                                      Timers::Timer &sync_timer);

// check if scheduled function overwrote any of the poison surrounding
// allocated memory
static void CheckFence(cGH const *const cctkGH, cFunctionData *attribute);

/// Traverse one function on all components of one refinement level
/// of one multigrid level.
int CallFunction(void *function,           ///< the function to call
                 cFunctionData *attribute, ///< attributes of the function
                 void *data) ///< private data for CCTK_CallFunction
{
  DECLARE_CCTK_PARAMETERS;

  static Timers::Timer total_timer("CallFunction");
  static Timers::Timer user_timer("thorns");
  static Timers::Timer sync_timer("syncs");

  total_timer.start();

  cGH *cctkGH = static_cast<cGH *>(data);

  assert(
      int(not not attribute->meta) + int(not not attribute->meta_early) +
          int(not not attribute->meta_late) + int(not not attribute->global) +
          int(not not attribute->global_early) +
          int(not not attribute->global_late) + int(not not attribute->level) +
          int(not not attribute->singlemap) + int(not not attribute->local) <=
      1);

  assert(not not attribute->loop_global + not not attribute->loop_level +
             not not attribute->loop_singlemap +
             not not attribute->loop_local <=
         1);

  // Create list of all groups that need to be synchronised
  vector<int> sync_groups;
  sync_groups.reserve(attribute->n_SyncGroups);
  for (int g = 0; g < attribute->n_SyncGroups; g++) {
    const int group = attribute->SyncGroups[g];
    if (CCTK_NumVarsInGroupI(group) > 0) {
      // don't add empty groups from the list
      sync_groups.push_back(group);
    }
  }

  if (attribute->meta or attribute->meta_early or attribute->meta_late or
      is_meta_mode()) {
    // Convtest operation

    if ((attribute->meta and do_meta_mode) or
        (attribute->meta_early and do_early_meta_mode) or
        (attribute->meta_late and do_late_meta_mode) or is_meta_mode()) {
      if (attribute->loop_local) {
        BEGIN_META_MODE(cctkGH) {
          BEGIN_MGLEVEL_LOOP(cctkGH) {
            BEGIN_REFLEVEL_LOOP(cctkGH) {
              BEGIN_LOCAL_MAP_LOOP(cctkGH, CCTK_GF) {
                BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, CCTK_GF) {
                  CallScheduledFunction("Meta time local mode", function,
                                        attribute, data, user_timer);
                }
                END_LOCAL_COMPONENT_LOOP;
              }
              END_LOCAL_MAP_LOOP;
              if (not sync_groups.empty()) {
                SyncGroupsInScheduleBlock(attribute, cctkGH, sync_groups,
                                          sync_timer);
              }
            }
            END_REFLEVEL_LOOP;
          }
          END_MGLEVEL_LOOP;
        }
        END_META_MODE;
      } else if (attribute->loop_singlemap) {
        BEGIN_META_MODE(cctkGH) {
          BEGIN_MGLEVEL_LOOP(cctkGH) {
            BEGIN_REFLEVEL_LOOP(cctkGH) {
              BEGIN_MAP_LOOP(cctkGH, CCTK_GF) {
                CallScheduledFunction("Meta time singlemap mode", function,
                                      attribute, data, user_timer);
              }
              END_MAP_LOOP;
              if (not sync_groups.empty()) {
                SyncGroupsInScheduleBlock(attribute, cctkGH, sync_groups,
                                          sync_timer);
              }
            }
            END_REFLEVEL_LOOP;
          }
          END_MGLEVEL_LOOP;
        }
        END_META_MODE;
      } else if (attribute->loop_level) {
        BEGIN_META_MODE(cctkGH) {
          BEGIN_MGLEVEL_LOOP(cctkGH) {
            BEGIN_REFLEVEL_LOOP(cctkGH) {
              CallScheduledFunction("Meta time level mode", function, attribute,
                                    data, user_timer);
              if (not sync_groups.empty()) {
                SyncGroupsInScheduleBlock(attribute, cctkGH, sync_groups,
                                          sync_timer);
              }
            }
            END_REFLEVEL_LOOP;
          }
          END_MGLEVEL_LOOP;
        }
        END_META_MODE;
      } else if (attribute->loop_global) {
        BEGIN_META_MODE(cctkGH) {
          BEGIN_MGLEVEL_LOOP(cctkGH) {
            CallScheduledFunction("Meta time global mode", function, attribute,
                                  data, user_timer);
            if (not sync_groups.empty()) {
              BEGIN_REFLEVEL_LOOP(cctkGH) {
                SyncGroupsInScheduleBlock(attribute, cctkGH, sync_groups,
                                          sync_timer);
              }
              END_REFLEVEL_LOOP;
            }
          }
          END_MGLEVEL_LOOP;
        }
        END_META_MODE;
      } else {
        BEGIN_META_MODE(cctkGH) {
          CallScheduledFunction("Meta mode", function, attribute, data,
                                user_timer);
          if (not sync_groups.empty()) {
            BEGIN_MGLEVEL_LOOP(cctkGH) {
              BEGIN_REFLEVEL_LOOP(cctkGH) {
                SyncGroupsInScheduleBlock(attribute, cctkGH, sync_groups,
                                          sync_timer);
              }
              END_REFLEVEL_LOOP;
            }
            END_MGLEVEL_LOOP;
          }
        }
        END_META_MODE;
      }
    }

  } else if (attribute->global or attribute->global_early or
             attribute->global_late or is_global_mode()) {
    // Global operation: call once

    if ((attribute->global and do_global_mode) or
        (attribute->global_early and do_early_global_mode) or
        (attribute->global_late and do_late_global_mode) or is_global_mode()) {
      if (attribute->loop_local) {
        BEGIN_GLOBAL_MODE(cctkGH) {
          BEGIN_REFLEVEL_LOOP(cctkGH) {
            BEGIN_LOCAL_MAP_LOOP(cctkGH, CCTK_GF) {
              BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, CCTK_GF) {
                CallScheduledFunction("Global time local mode", function,
                                      attribute, data, user_timer);
              }
              END_LOCAL_COMPONENT_LOOP;
            }
            END_LOCAL_MAP_LOOP;
            if (not sync_groups.empty()) {
              SyncGroupsInScheduleBlock(attribute, cctkGH, sync_groups,
                                        sync_timer);
            }
          }
          END_REFLEVEL_LOOP;
        }
        END_GLOBAL_MODE;
      } else if (attribute->loop_singlemap) {
        BEGIN_GLOBAL_MODE(cctkGH) {
          BEGIN_REFLEVEL_LOOP(cctkGH) {
            BEGIN_MAP_LOOP(cctkGH, CCTK_GF) {
              CallScheduledFunction("Global time singlemap mode", function,
                                    attribute, data, user_timer);
            }
            END_MAP_LOOP;
            if (not sync_groups.empty()) {
              SyncGroupsInScheduleBlock(attribute, cctkGH, sync_groups,
                                        sync_timer);
            }
          }
          END_REFLEVEL_LOOP;
        }
        END_GLOBAL_MODE;
      } else if (attribute->loop_level) {
        BEGIN_GLOBAL_MODE(cctkGH) {
          BEGIN_REFLEVEL_LOOP(cctkGH) {
            CallScheduledFunction("Global time level mode", function, attribute,
                                  data, user_timer);
            if (not sync_groups.empty()) {
              SyncGroupsInScheduleBlock(attribute, cctkGH, sync_groups,
                                        sync_timer);
            }
          }
          END_REFLEVEL_LOOP;
        }
        END_GLOBAL_MODE;
      } else {
        BEGIN_GLOBAL_MODE(cctkGH) {
          CallScheduledFunction("Global mode", function, attribute, data,
                                user_timer);
          if (not sync_groups.empty()) {
            BEGIN_REFLEVEL_LOOP(cctkGH) {
              SyncGroupsInScheduleBlock(attribute, cctkGH, sync_groups,
                                        sync_timer);
            }
            END_REFLEVEL_LOOP;
          }
        }
        END_GLOBAL_MODE;
      }
    }

  } else if (attribute->level) {
    // Level operation: call once per refinement level

    if (attribute->loop_local) {
      BEGIN_LOCAL_MAP_LOOP(cctkGH, CCTK_GF) {
        BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, CCTK_GF) {
          CallScheduledFunction("Level time local mode", function, attribute,
                                data, user_timer);
        }
        END_LOCAL_COMPONENT_LOOP;
      }
      END_LOCAL_MAP_LOOP;
    } else if (attribute->loop_singlemap) {
      BEGIN_MAP_LOOP(cctkGH, CCTK_GF) {
        CallScheduledFunction("Level time singlemap mode", function, attribute,
                              data, user_timer);
      }
      END_MAP_LOOP;
    } else {
      CallScheduledFunction("Level mode", function, attribute, data,
                            user_timer);
    }
    if (not sync_groups.empty()) {
      SyncGroupsInScheduleBlock(attribute, cctkGH, sync_groups, sync_timer);
    }

  } else if (attribute->singlemap) {
    // Single map operation: call once per refinement level and map

    if (attribute->loop_local) {
      BEGIN_LOCAL_MAP_LOOP(cctkGH, CCTK_GF) {
        BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, CCTK_GF) {
          CallScheduledFunction("Singlemap time local mode", function,
                                attribute, data, user_timer);
        }
        END_LOCAL_COMPONENT_LOOP;
      }
      END_LOCAL_MAP_LOOP;
    } else {
      BEGIN_MAP_LOOP(cctkGH, CCTK_GF) {
        CallScheduledFunction("Singlemap mode", function, attribute, data,
                              user_timer);
      }
      END_MAP_LOOP;
    }
    if (not sync_groups.empty()) {
      SyncGroupsInScheduleBlock(attribute, cctkGH, sync_groups, sync_timer);
    }

  } else {
    // Local operation: call once per component

    BEGIN_LOCAL_MAP_LOOP(cctkGH, CCTK_GF) {
      BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, CCTK_GF) {
        CallScheduledFunction("Local mode", function, attribute, data,
                              user_timer);
      }
      END_LOCAL_COMPONENT_LOOP;
    }
    END_LOCAL_MAP_LOOP;
    if (not sync_groups.empty()) {
      SyncGroupsInScheduleBlock(attribute, cctkGH, sync_groups, sync_timer);
    }
  }

  if (schedule_barriers) {
    // Create an ID that is almost unique for this scheduled
    // function call
    stringstream buf;
    buf << cctkGH->cctk_iteration << "\n";
    buf << attribute->meta << attribute->meta_early << attribute->meta_late
        << attribute->global << attribute->global_early
        << attribute->global_late << attribute->level << attribute->singlemap
        << attribute->local << "\n";
    buf << attribute->where << "\n";
    buf << attribute->thorn << "\n";
    buf << attribute->routine << "\n";
    string const str = buf.str();
    int const id = adler32(str.c_str(), str.length());
    static Timers::Timer barrier_timer("barrier");
    barrier_timer.start();
    Carpet::NamedBarrier(NULL, id, "Carpet::CallFunction");
    barrier_timer.stop();
  }

  total_timer.stop();

  // The return value indicates whether the grid functions have been
  // synchronised.
  // 0: let the flesh do the synchronisation
  // 1: we did the synchronisation
  return 1;
}

void CallScheduledFunction(char const *restrict const time_and_mode,
                           void *const function, cFunctionData *const attribute,
                           void *const data, Timers::Timer &user_timer) {
  cGH const *const cctkGH = static_cast<cGH const *>(data);
  Checkpoint("%s call at %s to %s::%s", time_and_mode, attribute->where,
             attribute->thorn, attribute->routine);
  int const skip = CallBeforeRoutines(cctkGH, function, attribute, data);
  if (not skip) {
    Timers::Timer timer(attribute->routine);

    // Save the time step size
    CCTK_REAL const saved_cctk_delta_time = cctkGH->cctk_delta_time;

    user_timer.start();
#ifdef REQUIREMENTS_HH
    Requirements::BeforeRoutine(attribute, cctkGH->cctk_iteration, reflevel,
                                map, timelevel, timelevel_offset);
#endif
    timer.start();
    if (CCTK_IsFunctionAliased("Accelerator_PreCallFunction")) {
      Timers::Timer pre_timer("PreCall");
      pre_timer.start();
      Accelerator_PreCallFunction(cctkGH, attribute);
      pre_timer.stop();
    }
    int const res = CCTK_CallFunction(function, attribute, data);
    assert(res == 0);
    if (CCTK_IsFunctionAliased("Accelerator_PostCallFunction")) {
      Timers::Timer post_timer("PostCall");
      post_timer.start();
      Accelerator_PostCallFunction(cctkGH, attribute);
      post_timer.stop();
    }
    timer.stop();
    CheckFence(cctkGH, attribute);
#ifdef REQUIREMENTS_HH
    Requirements::AfterRoutine(attribute, cctkGH->cctk_iteration, reflevel, map,
                               timelevel, timelevel_offset);
#endif
    user_timer.stop();

    // Manage the time step size. If the time step size changes
    // during initialisation, assume it is thorn Time, and update
    // the time hierarchy. If it changes during evolution, assume it
    // is MoL, and do nothing.
    if (cctkGH->cctk_iteration == 0 and
        cctkGH->cctk_delta_time != saved_cctk_delta_time) {
      // The user changed cctk_delta_time during initialisation --
      // update our internals and the time hierarchy
      bool const is_global = attribute->meta or attribute->meta_early or
                             attribute->meta_late or attribute->global or
                             attribute->global_early or attribute->global_late;
      delta_time = cctkGH->cctk_delta_time / mglevelfact *
                   (is_global ? 1.0 : timereflevelfact);
      for (int ml = 0; ml < mglevels; ++ml) {
        for (int rl = 0; rl < reflevels; ++rl) {
          // Update the time delta
          CCTK_REAL const dt =
              delta_time / timereffacts.AT(rl) * ipow(mgfact, ml);
          tt->set_delta(ml, rl, dt);
          CCTK_REAL const t0 = tt->get_time(ml, rl, 0);
          // Update the times of the past timelevels
          for (int tl = 1; tl < timelevels; ++tl) {
            CCTK_REAL const t = t0 - tl * dt;
            tt->set_time(ml, rl, tl, t);
          }
        }
      }
    }
  }
  CallAfterRoutines(cctkGH, function, attribute, data);
}

void SyncGroupsInScheduleBlock(cFunctionData *attribute, cGH *cctkGH,
                               vector<int> const &sync_groups,
                               Timers::Timer &sync_timer) {
  DECLARE_CCTK_PARAMETERS;

  if (sync_barriers) {
    // Create an ID that is almost unique for this scheduled
    // function call
    // stringstream buf;
    // buf << cctkGH->cctk_iteration << "\n";
    // buf << attribute->meta
    //     << attribute->meta_early
    //     << attribute->meta_late
    //     << attribute->global
    //     << attribute->global_early
    //     << attribute->global_late
    //     << attribute->level
    //     << attribute->singlemap
    //     << attribute->local << "\n";
    // buf << attribute->where << "\n";
    // buf << attribute->thorn << "\n";
    // buf << attribute->routine << " sync\n";
    // string const str = buf.str();
    // int const id = adler32(str.c_str(), str.length());
    int const id = 101;
    static Timers::Timer barrier_timer("pre_sync_barrier");
    barrier_timer.start();
    Carpet::NamedBarrier(NULL, id, "Carpet::Sync");
    barrier_timer.stop();
  }

  sync_timer.start();
  SyncProlongateGroups(cctkGH, sync_groups, attribute);
  sync_timer.stop();

  if (sync_barriers) {
    // Create an ID that is almost unique for this scheduled
    // function call
    // stringstream buf;
    // buf << cctkGH->cctk_iteration << "\n";
    // buf << attribute->meta
    //     << attribute->meta_early
    //     << attribute->meta_late
    //     << attribute->global
    //     << attribute->global_early
    //     << attribute->global_late
    //     << attribute->level
    //     << attribute->singlemap
    //     << attribute->local << "\n";
    // buf << attribute->where << "\n";
    // buf << attribute->thorn << "\n";
    // buf << attribute->routine << " sync\n";
    // string const str = buf.str();
    // int const id = adler32(str.c_str(), str.length());
    int const id = 102;
    static Timers::Timer barrier_timer("post_sync_barrier");
    barrier_timer.start();
    Carpet::NamedBarrier(NULL, id, "Carpet::Sync");
    barrier_timer.stop();
  }
}

void CheckFence(cGH const *const cctkGH, cFunctionData *attribute) {
  DECLARE_CCTK_PARAMETERS;

  if (is_meta_mode())
    return; // meta mode has no accessible variables
  // enumerating the grid variables is expensive
  if (not gdata::fence_is_energized())
    return;

  Timers::Timer timer("FenceCheck");
  timer.start();

  for (int group = 0; group < CCTK_NumGroups(); ++group) {
    if (CCTK_QueryGroupStorageI(cctkGH, group)) {
      int const nvar = CCTK_NumVarsInGroupI(group);
      if (nvar == 0)
        continue;
      int const n0 = CCTK_FirstVarIndexI(group);
      assert(n0 >= 0);

      int const num_tl = CCTK_ActiveTimeLevelsVI(cctkGH, n0);
      if (num_tl == 0)
        continue;
      int const min_tl = 0;
      int const max_tl = num_tl - 1;

      int const grouptype = CCTK_GroupTypeI(group);
      if (grouptype == CCTK_GF && not is_local_mode())
        continue;
      assert(not is_meta_mode());

      // FIXME: query is_XXX() functions
      const bool is_array = grouptype == CCTK_ARRAY or grouptype == CCTK_SCALAR;
      const int m = is_array ? 0 : map;
      const int rl = is_array ? 0 : reflevel;
      const int ml = is_array ? 0 : mglevel;
      const int lc = is_array ? 0 : local_component;

      {
        char *const groupname = CCTK_GroupName(group);
        Checkpoint("FenceCheck \"%s\"", groupname);
        free(groupname);
      }

      for (int var = 0; var < nvar; ++var) {
        assert(n0 + var < CCTK_NumVars());
        for (int tl = min_tl; tl <= max_tl; ++tl) {
          // TODO: turns this into a call of ggf taken a callback routine
          ggf *const ff = arrdata.AT(group).AT(m).data.AT(var);
          assert(ff);
          gdata *const data = ff->data_pointer(tl, rl, lc, ml);
          for (int f = 0; f < 2; ++f) {
            if (not data->check_fence(f)) {
              char *fullname = CCTK_FullName(n0 + var);
              CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                          "At iteration %d: timelevel %d, component %d, map "
                          "%d, refinement level %d of the variable \"%s\" "
                          "contains was written to beyond to %s bound of the "
                          "allocated memory by %s::%s at %s.",
                          cctkGH->cctk_iteration, tl, component, m, rl,
                          fullname, f ? "lower" : "upper", attribute->thorn,
                          attribute->routine, attribute->where);
              free(fullname);
            }
          } // for f
        }   // for tl
      }     // for var
    }
  }
  timer.stop();
}

} // namespace Carpet
