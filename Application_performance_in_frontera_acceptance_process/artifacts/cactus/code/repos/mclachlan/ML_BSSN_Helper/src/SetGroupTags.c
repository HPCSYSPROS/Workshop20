#include <cctk.h>
#include <cctk_Parameters.h>
#include <util_Table.h>

#include <assert.h>
#include <stdbool.h>

#define THORN_ML_BSSN // "ML_BSSN" will be replaced
#ifdef THORN_ML_CCZ4
#define HAVE_THETA
#endif

static void set_group_tags(bool checkpoint, bool persistent, bool prolongate,
                           const char *gn);

int ML_BSSN_SetGroupTags(void) {
  DECLARE_CCTK_PARAMETERS;

  // global describes whether these variables are valid on the whole
  // grid hierarchy, or only on the current level

  bool global = timelevels > 1;
  set_group_tags(global, global, global, "ADMBase::metric");
  set_group_tags(global, global, global, "ADMBase::curv");
  set_group_tags(global, global, global, "ADMBase::lapse");
  set_group_tags(global, global, global, "ADMBase::shift");
  set_group_tags(global, global, global, "ADMBase::dtlapse");
  set_group_tags(global, global, global, "ADMBase::dtshift");

  bool other_global = other_timelevels > 1;
  set_group_tags(other_global, other_global, other_global,
                 "ML_BSSN::ML_cons_detg");
  set_group_tags(other_global, other_global, other_global,
                 "ML_BSSN::ML_cons_Gamma");
  set_group_tags(other_global, other_global, other_global,
                 "ML_BSSN::ML_cons_traceA");
  set_group_tags(other_global, other_global, other_global, "ML_BSSN::ML_Ham");
  set_group_tags(other_global, other_global, other_global, "ML_BSSN::ML_mom");

  // RHS variables do not have valid boundaries
  bool rhs_global = rhs_timelevels > 1;
  set_group_tags(rhs_global, rhs_global, false, "ML_BSSN::ML_log_confacrhs");
  set_group_tags(rhs_global, rhs_global, false, "ML_BSSN::ML_metricrhs");
  set_group_tags(rhs_global, rhs_global, false, "ML_BSSN::ML_Gammarhs");

#ifdef HAVE_THETA
  set_group_tags(rhs_global, rhs_global, false, "ML_BSSN::ML_Thetarhs");
#endif
  set_group_tags(rhs_global, rhs_global, false, "ML_BSSN::ML_trace_curvrhs");
  set_group_tags(rhs_global, rhs_global, false, "ML_BSSN::ML_curvrhs");
  set_group_tags(rhs_global, rhs_global, false, "ML_BSSN::ML_lapserhs");
  set_group_tags(rhs_global, rhs_global, false, "ML_BSSN::ML_dtlapserhs");
  set_group_tags(rhs_global, rhs_global, false, "ML_BSSN::ML_shiftrhs");
  set_group_tags(rhs_global, rhs_global, false, "ML_BSSN::ML_dtshiftrhs");

  return 0;
}

static void set_group_tags(bool checkpoint, bool persistent, bool prolongate,
                           const char *gn) {
  assert(gn);

  int gi = CCTK_GroupIndex(gn);
  assert(gi >= 0);

  int table = CCTK_GroupTagsTableI(gi);
  assert(table >= 0);

  if (!checkpoint) {
    int ierr = Util_TableSetString(table, "no", "Checkpoint");
    assert(!ierr);
  }

  if (!persistent) {
    int ierr = Util_TableSetString(table, "no", "Persistent");
    assert(!ierr);
  }

/* Ian Hinder says that Kranc handles this now automatically. After
   confirming this, remove this code (and this note). */
#if 0
  if (!prolongate) {
    int iret;
    iret = Util_TableDeleteKey(table, "ProlongationParameter");
    assert(iret == 0 || iret == UTIL_ERROR_TABLE_NO_SUCH_KEY);
    iret = Util_TableDeleteKey(table, "Prolongation");
    assert(iret == 0 || iret == UTIL_ERROR_TABLE_NO_SUCH_KEY);
    int ierr = Util_TableSetString(table, "none", "Prolongation");
    /* 0 means key did not exist before; 1 means key existed before
       and has now been reset */
    assert(ierr == 0 || ierr == 1);
  }
#endif
}
