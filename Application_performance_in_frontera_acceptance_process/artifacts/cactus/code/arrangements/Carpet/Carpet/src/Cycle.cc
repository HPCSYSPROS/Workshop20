#include <cassert>
#include <cstdlib>

#include <cctk.h>
#include <cctk_Parameters.h>

#include <Requirements.hh>

#include <Timer.hh>

#include <ggf.hh>
#include <gh.hh>

#include <carpet.hh>

namespace Carpet {

using namespace std;

void CycleTimeLevels(cGH *const cctkGH) {
  DECLARE_CCTK_PARAMETERS;

  Timers::Timer timer("CycleTimeLevels");
  timer.start();

  Checkpoint("CycleTimeLevels");
  assert(is_level_mode());

  assert(timelevel == 0);
  tt->advance_time(mglevel, reflevel);
  if (not adaptive_stepsize) {
    cctkGH->cctk_time = tt->get_time(mglevel, reflevel, timelevel);
  }

#ifdef REQUIREMENTS_HH
  Requirements::Cycle(reflevel);
#endif

  int errors = 0;

  for (int group = 0; group < CCTK_NumGroups(); ++group) {
    if (CCTK_QueryGroupStorageI(cctkGH, group)) {

      int const activetimelevels = CCTK_ActiveTimeLevelsGI(cctkGH, group);

      if (activetimelevels > 1) {
        switch (groupdata.AT(group).transport_operator) {
        case op_Lagrange:
        case op_ENO:
        case op_WENO:
        case op_Lagrange_monotone:
          if (activetimelevels > 1) {
            if (activetimelevels < prolongation_order_time + 1) {
              char *const groupname = CCTK_GroupName(group);
              CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,
                         "Group \"%s\" has %d only active time levels.  Groups "
                         "with more than one active time level need at least "
                         "%d active time levels for prolongation_order_time=%d",
                         groupname, activetimelevels,
                         int(prolongation_order_time + 1),
                         int(prolongation_order_time));
              free(groupname);
              ++errors;
            }
          }
          break;
        default:; // do nothing
        }

        switch (CCTK_GroupTypeI(group)) {

        case CCTK_GF:
          assert(reflevel >= 0 and reflevel < reflevels);
          for (int m = 0; m < (int)arrdata.AT(group).size(); ++m) {
            for (int var = 0; var < CCTK_NumVarsInGroupI(group); ++var) {
              arrdata.AT(group).AT(m).data.AT(var)->cycle_all(reflevel,
                                                              mglevel);
            }
          }
          break;

        case CCTK_SCALAR:
        case CCTK_ARRAY:
          if (do_global_mode) {
            int const numtimelevels = CCTK_MaxActiveTimeLevelsVI(cctkGH, group);
            int const firstvarindex = CCTK_FirstVarIndexI(group);
            for (int var = 0; var < CCTK_NumVarsInGroupI(group); ++var) {
              arrdata.AT(group).AT(0).data.AT(var)->cycle_all(0, mglevel);
              {
                int const varindex = firstvarindex + var;
                for (int tl = 0; tl < numtimelevels; ++tl) {
                  if (tl < groupdata.AT(group).info.activetimelevels) {
                    ggf *const ff = arrdata.AT(group).AT(0).data.AT(var);
                    void *const ptr = ff->data_pointer(tl, 0, 0, 0)->storage();
                    cctkGH->data[varindex][tl] = ptr;
                  } else {
                    cctkGH->data[varindex][tl] = NULL;
                  }
                }
              }
            }
          }
          break;

        default:
          assert(0);
        } // switch grouptype
      }   // if activetimelevels > 1
    }     // if storage
  }       // for group

  if (CCTK_IsFunctionAliased("Accelerator_Cycle")) {
    Accelerator_Cycle(cctkGH);
  }

  if (errors > 0) {
    CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Errors in %d groups detected; aborting", errors);
  }
  timer.stop();
}

void UncycleTimeLevels(cGH *const cctkGH) {
  DECLARE_CCTK_PARAMETERS;

  Checkpoint("UncycleTimeLevels");
  assert(is_level_mode());

  assert(timelevel == 0);
  tt->retreat_time(mglevel, reflevel);
  assert(not adaptive_stepsize);
  cctkGH->cctk_time = tt->get_time(mglevel, reflevel, timelevel);

  for (int group = 0; group < CCTK_NumGroups(); ++group) {
    if (CCTK_QueryGroupStorageI(cctkGH, group)) {

      switch (CCTK_GroupTypeI(group)) {

      case CCTK_GF:
        assert(reflevel >= 0 and reflevel < reflevels);
        for (int m = 0; m < (int)arrdata.AT(group).size(); ++m) {
          for (int var = 0; var < CCTK_NumVarsInGroupI(group); ++var) {
            arrdata.AT(group).AT(m).data.AT(var)->uncycle_all(reflevel,
                                                              mglevel);
          }
        }
        break;

      case CCTK_SCALAR:
      case CCTK_ARRAY:
        if (do_global_mode) {
          int const numtimelevels = CCTK_MaxActiveTimeLevelsVI(cctkGH, group);
          int const firstvarindex = CCTK_FirstVarIndexI(group);
          for (int var = 0; var < CCTK_NumVarsInGroupI(group); ++var) {
            arrdata.AT(group).AT(0).data.AT(var)->uncycle_all(0, mglevel);
            {
              int const varindex = firstvarindex + var;
              for (int tl = 0; tl < numtimelevels; ++tl) {
                if (tl < groupdata.AT(group).info.activetimelevels) {
                  ggf *const ff = arrdata.AT(group).AT(0).data.AT(var);
                  void *const ptr = ff->data_pointer(tl, 0, 0, 0)->storage();
                  cctkGH->data[varindex][tl] = ptr;
                } else {
                  cctkGH->data[varindex][tl] = NULL;
                }
              }
            }
          }
        }
        break;

      default:
        assert(0);
      } // switch grouptype
    }   // if storage
  }     // for group
}

void FlipTimeLevels(cGH *const cctkGH) {
  DECLARE_CCTK_PARAMETERS;

  Checkpoint("FlipTimeLevels");
  assert(is_level_mode());

  assert(timelevel == 0);
  tt->flip_timelevels(mglevel, reflevel);
  assert(not adaptive_stepsize);
  cctkGH->cctk_time = tt->get_time(mglevel, reflevel, timelevel);
  cctkGH->cctk_delta_time *= -1;

  for (int group = 0; group < CCTK_NumGroups(); ++group) {
    if (CCTK_QueryGroupStorageI(cctkGH, group)) {

      switch (CCTK_GroupTypeI(group)) {

      case CCTK_GF:
        assert(reflevel >= 0 and reflevel < reflevels);
        for (int m = 0; m < (int)arrdata.AT(group).size(); ++m) {
          for (int var = 0; var < CCTK_NumVarsInGroupI(group); ++var) {
            arrdata.AT(group).AT(m).data.AT(var)->flip_all(reflevel, mglevel);
          }
        }
        break;

      case CCTK_SCALAR:
      case CCTK_ARRAY:
        if (do_global_mode) {
          int const numtimelevels = CCTK_MaxActiveTimeLevelsVI(cctkGH, group);
          int const firstvarindex = CCTK_FirstVarIndexI(group);
          for (int var = 0; var < CCTK_NumVarsInGroupI(group); ++var) {
            arrdata.AT(group).AT(0).data.AT(var)->flip_all(0, mglevel);
            {
              int const varindex = firstvarindex + var;
              for (int tl = 0; tl < numtimelevels; ++tl) {
                if (tl < groupdata.AT(group).info.activetimelevels) {
                  ggf *const ff = arrdata.AT(group).AT(0).data.AT(var);
                  void *const ptr = ff->data_pointer(tl, 0, 0, 0)->storage();
                  cctkGH->data[varindex][tl] = ptr;
                } else {
                  cctkGH->data[varindex][tl] = NULL;
                }
              }
            }
          }
        }
        break;

      default:
        assert(0);
      } // switch grouptype
    }   // if storage
  }     // for group
}

void FillTimeLevels(const cGH *const cctkGH) {
  Checkpoint("FillTimeLevels");
  assert(is_level_mode());

  for (int group = 0; group < CCTK_NumGroups(); ++group) {
    if (CCTK_QueryGroupStorageI(cctkGH, group)) {
      switch (CCTK_GroupTypeI(group)) {

      case CCTK_GF:
        assert(reflevel >= 0 and reflevel < reflevels);
        for (int m = 0; m < (int)arrdata.AT(group).size(); ++m) {
          const bool have_accel =
              CCTK_IsFunctionAliased("Accelerator_NotifyDataModified");
          vector<CCTK_INT> vis, rls, tls;
          const int varn = CCTK_NumVarsInGroupI(group);
          for (int var = 0; var < varn; ++var) {
            arrdata.AT(group).AT(m).data.AT(var)->fill_all(reflevel, mglevel);
            if (have_accel) {
              const int var0 = CCTK_FirstVarIndexI(group);
              const int num_tl =
                  arrdata.AT(group).AT(m).data.AT(var)->timelevels(mglevel,
                                                                   reflevel);
              for (int tl = 1; tl < num_tl; tl++) {
                vis.push_back(var0 + var);
                rls.push_back(reflevel);
                tls.push_back(tl);
              }
            }
          }

          if (have_accel) {
            const CCTK_INT on_device = 0;
            Accelerator_NotifyDataModified(cctkGH, &vis.front(), &rls.front(),
                                           &tls.front(), vis.size(), on_device);
          }
        }
        break;

      case CCTK_SCALAR:
      case CCTK_ARRAY:
        if (do_global_mode) {
          for (int var = 0; var < CCTK_NumVarsInGroupI(group); ++var) {
            arrdata.AT(group).AT(0).data.AT(var)->fill_all(0, mglevel);
          }
        }
        break;

      default:
        assert(0);
      } // switch grouptype
    }   // if storage
  }     // for group
}

} // namespace Carpet
