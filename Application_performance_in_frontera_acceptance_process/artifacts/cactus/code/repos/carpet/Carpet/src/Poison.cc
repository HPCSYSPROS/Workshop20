#include <cassert>
#include <cstdlib>
#include <cstring>

#include <cctk.h>
#include <cctk_Parameters.h>
#include <util_ErrorCodes.h>
#include <util_Table.h>

#include <Timer.hh>

#include <defs.hh>
#include <typeprops.hh>

#include <carpet.hh>

namespace Carpet {

using namespace std;

// The parameter where specifies which time levels should be
// poisoned.  what specifies what kind of grid variables should be
// poisoned.
void Poison(cGH const *const cctkGH, checktimes const where, int const what) {
  DECLARE_CCTK_PARAMETERS;

  assert(what == 0 or what == CCTK_GF or what == CCTK_ARRAY);

  if (not poison_new_timelevels)
    return;

  Timers::Timer timer("Poison");
  timer.start();

  for (int group = 0; group < CCTK_NumGroups(); ++group) {
    if (CCTK_QueryGroupStorageI(cctkGH, group)) {
      int const grouptype = CCTK_GroupTypeI(group);
      if (what == 0 or (what == CCTK_GF and grouptype == CCTK_GF) or
          (what == CCTK_ARRAY and
           (grouptype == CCTK_ARRAY or grouptype == CCTK_SCALAR))) {
        PoisonGroup(cctkGH, group, where);
      }
    } // if has storage
  }   // for group
  timer.stop();
}

void PoisonGroup(cGH const *const cctkGH, int const group,
                 checktimes const where) {
  DECLARE_CCTK_PARAMETERS;

  assert(group >= 0 and group < CCTK_NumGroups());

  if (not poison_new_timelevels)
    return;

  if (not CCTK_QueryGroupStorageI(cctkGH, group)) {
    char *const groupname = CCTK_GroupName(group);
    CCTK_VWarn(2, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Cannot poison group \"%s\" because it has no storage",
               groupname);
    free(groupname);
    return;
  }

  int const nvar = CCTK_NumVarsInGroupI(group);
  if (nvar == 0)
    return;
  int const n0 = CCTK_FirstVarIndexI(group);
  assert(n0 >= 0);
  int const sz = CCTK_VarTypeSize(CCTK_VarTypeI(n0));
  assert(sz > 0);

  int const table = CCTK_GroupTagsTableI(group);
  assert(table >= 0);
  bool persistent;
  char buf[100];
  int const ilen = Util_TableGetString(table, sizeof buf, buf, "Persistent");
  if (ilen > 0) {
    if (CCTK_EQUALS(buf, "yes")) {
      persistent = true;
    } else if (CCTK_EQUALS(buf, "no")) {
      persistent = false;
    } else {
      assert(0);
    }
  } else if (ilen == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
    // default
    persistent = true;
  } else {
    assert(0);
  }

  int const num_tl = CCTK_ActiveTimeLevelsVI(cctkGH, n0);
  assert(num_tl > 0);
  int const min_tl = min_timelevel(where, num_tl, persistent);
  int const max_tl = max_timelevel(where, num_tl, persistent);

  if (min_tl <= max_tl) {

    {
      char *const groupname = CCTK_GroupName(group);
      Checkpoint("PoisonGroup \"%s\"", groupname);
      free(groupname);
    }

    CCTK_INT const poison_value = get_poison_value();

    int const grouptype = CCTK_GroupTypeI(group);

    BEGIN_LOCAL_MAP_LOOP(cctkGH, grouptype) {
      BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, grouptype) {

        ivect size(1);
        int const gpdim = groupdata.AT(group).info.dim;
        for (int d = 0; d < gpdim; ++d) {
          size[d] = groupdata.AT(group).info.ash[d];
        }
        int const np = prod(size);

        for (int var = 0; var < nvar; ++var) {
          int const n = n0 + var;
          for (int tl = min_tl; tl <= max_tl; ++tl) {
            memset(cctkGH->data[n][tl], poison_value, size_t(np) * sz);
          } // for tl
        }   // for var
      }
      END_LOCAL_COMPONENT_LOOP;
    }
    END_LOCAL_MAP_LOOP;

  } // if tl
}

void PoisonCheck(cGH const *const cctkGH, checktimes const where) {
  DECLARE_CCTK_PARAMETERS;

  if (not check_for_poison)
    return;

  Checkpoint("PoisonCheck");
  Timers::Timer timer("PoisonCheck");
  timer.start();

  for (int group = 0; group < CCTK_NumGroups(); ++group) {
    int const nvar = CCTK_NumVarsInGroupI(group);
    if (nvar > 0 and CCTK_QueryGroupStorageI(cctkGH, group)) {

      int const grouptype = CCTK_GroupTypeI(group);
      int const n0 = CCTK_FirstVarIndexI(group);
      assert(n0 >= 0);
      int const tp = CCTK_VarTypeI(n0);
      int const gpdim = groupdata.AT(group).info.dim;

      int const table = CCTK_GroupTagsTableI(group);
      assert(table >= 0);
      bool persistent;
      char buf[100];
      int const ilen =
          Util_TableGetString(table, sizeof buf, buf, "Persistent");
      if (ilen > 0) {
        if (CCTK_EQUALS(buf, "yes")) {
          persistent = true;
        } else if (CCTK_EQUALS(buf, "no")) {
          persistent = false;
        } else {
          assert(0);
        }
      } else if (ilen == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
        // default
        persistent = true;
      } else {
        assert(0);
      }

      CCTK_INT const poison_value = get_poison_value();

      int const num_tl = CCTK_ActiveTimeLevelsVI(cctkGH, n0);
      assert(num_tl > 0);
      int const min_tl = min_timelevel(where, num_tl, persistent);
      int const max_tl = max_timelevel(where, num_tl, persistent);

      BEGIN_LOCAL_MAP_LOOP(cctkGH, grouptype) {
        BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, grouptype) {

          ivect size(1), asize(1);
          for (int d = 0; d < gpdim; ++d) {
            size[d] = groupdata.AT(group).info.lsh[d];
            asize[d] = groupdata.AT(group).info.ash[d];
          }
          int const np = prod(size);

          for (int var = 0; var < nvar; ++var) {
            int const n = n0 + var;

            for (int tl = min_tl; tl <= max_tl; ++tl) {

              const void *const data = cctkGH->data[n][tl];
              int numpoison = 0;
              for (int k = 0; k < size[2]; ++k) {
                for (int j = 0; j < size[1]; ++j) {
                  for (int i = 0; i < size[0]; ++i) {
                    int const idx = i + asize[0] * (j + asize[1] * k);
                    bool poisoned = false;
                    switch (specific_cactus_type(tp)) {
#define TYPECASE(N, T)                                                         \
  case N: {                                                                    \
    T worm;                                                                    \
    memset(&worm, poison_value, sizeof worm);                                  \
    const T &val = ((const T *)data)[idx];                                     \
    poisoned = memcmp(&worm, &val, sizeof worm) == 0;                          \
    break;                                                                     \
  }
#include "typecase.hh"
#undef TYPECASE
                    default:
                      UnsupportedVarType(n);
                    }
                    if (poisoned) {
                      ++numpoison;
                      if (max_poison_locations == -1 or
                          numpoison <= max_poison_locations) {
                        char *fullname = CCTK_FullName(n);
                        CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                                   "At iteration %d: timelevel %d, component "
                                   "%d, map %d, refinement level %d of the "
                                   "variable \"%s\" contains poison at "
                                   "[%d,%d,%d]",
                                   cctkGH->cctk_iteration, tl, component, map,
                                   reflevel, fullname, i, j, k);
                        free(fullname);
                      }
                    } // if poisoned
                  }   // for i
                }     // for j
              }       // for k
              if (max_poison_locations != -1 and
                  numpoison > max_poison_locations) {
                char *fullname = CCTK_FullName(n);
                CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                           "At iteration %d: timelevel %d, component %d, map "
                           "%d, refinement level %d of the variable \"%s\" "
                           "contains poison at %d of %d locations; not all "
                           "locations were printed",
                           cctkGH->cctk_iteration, tl, component, map, reflevel,
                           fullname, numpoison, np);
                free(fullname);
              } else if (numpoison > 0) {
                CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                           "Found poison at %d of %d locations", numpoison, np);
              }

            } // for tl
          }   // for var
        }
        END_LOCAL_COMPONENT_LOOP;
      }
      END_LOCAL_MAP_LOOP;

    } // if has storage
  }   // for group
  timer.stop();
}

} // namespace Carpet
