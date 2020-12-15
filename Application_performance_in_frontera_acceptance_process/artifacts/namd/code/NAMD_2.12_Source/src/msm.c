/* msm.c */

#include "msm_defn.h"


NL_Msm *NL_msm_create(void) {
  NL_Msm *pm = (NL_Msm *) calloc(1, sizeof(NL_Msm));
  if (NULL == pm) return NULL;

  pm->timer = wkf_timer_create();
  if (NULL == pm->timer) return NULL;
  pm->timer_longrng = wkf_timer_create();
  if (NULL == pm->timer) return NULL;

  /* method parameters are modified with msm_configure() */
  pm->gridspacing = DEFAULT_GRIDSPACING;
  pm->approx = DEFAULT_APPROX;
  pm->split = DEFAULT_SPLIT;
  pm->nlevels = DEFAULT_NLEVELS;

  pm->density = DEFAULT_DENSITY;
  pm->binfill = DEFAULT_BINFILL;
  pm->nbinslots = DEFAULT_NBINSLOTS;
  return pm;
}


void NL_msm_destroy(NL_Msm *pm) {
  NL_msm_cleanup(pm);
  wkf_timer_destroy(pm->timer);
  wkf_timer_destroy(pm->timer_longrng);
  free(pm);
}


int NL_msm_configure(
    NL_Msm *pm,          /**< the MSM solver object */
    double gridspacing,  /**< grid spacing for first grid level */
    int approx,          /**< which approximation method */
    int split,           /**< which splitting */
    int nlevels          /**< number of grid levels to use */
    ) {
  if (gridspacing > 0) pm->gridspacing = gridspacing;
  else if (gridspacing < 0) return NL_MSM_ERROR_PARAM;
  if (approx >= 0 && approx < NL_MSM_APPROX_END) pm->approx = approx;
  else return NL_MSM_ERROR_PARAM;
  if (split >= 0 && split < NL_MSM_SPLIT_END) pm->split = split;
  else return NL_MSM_ERROR_PARAM;
  if (nlevels >= 0) pm->nlevels = nlevels;
  else return NL_MSM_ERROR_PARAM;

  return NL_MSM_SUCCESS;
}


int NL_msm_compute_force(
    NL_Msm *pm,          /**< the MSM solver object */
    double *felec,       /**< electrostatic forces x/y/z for each atom */
    double *uelec,       /**< electrostatic potential energy */
    const double *atom,  /**< positions and charge x/y/z/q for each atom */
    int natoms           /**< number of atoms */
    ) {
  int rc;
  double time_delta;

  if (pm->msmflags & NL_MSM_COMPUTE_SPREC) {
    return NL_MSM_ERROR_PARAM;
  }

  /* store buffer pointers from caller */
  pm->felec = felec;
  pm->atom = atom;
  pm->numatoms = natoms;
  pm->uelec = 0;

  if (pm->msmflags & NL_MSM_COMPUTE_SHORT_RANGE) {
    wkf_timer_start(pm->timer);
    rc = NL_msm_compute_short_range(pm);
    if (rc) return rc;
    wkf_timer_stop(pm->timer);
    time_delta = wkf_timer_time(pm->timer);
    if (pm->report_timings) {
      printf("MSM short-range part:  %6.3f sec\n", time_delta);
    }
  }
  if (pm->msmflags & NL_MSM_COMPUTE_LONG_RANGE) {
    wkf_timer_start(pm->timer);
    rc = NL_msm_compute_long_range(pm);
    if (rc) return rc;
    wkf_timer_stop(pm->timer);
    time_delta = wkf_timer_time(pm->timer);
    if (pm->report_timings) {
      printf("MSM long-range part:   %6.3f sec\n", time_delta);
    }
  }
  *uelec += pm->uelec;  /* add to caller's potential */

  return 0;
}


int NL_msm_compute_force_sprec(
    NL_Msm *pm,          /**< the MSM solver object */
    float *felec_f,      /**< electrostatic forces x/y/z for each atom */
    float *uelec_f,      /**< electrostatic potential energy */
    const float *atom_f, /**< positions and charge x/y/z/q for each atom */
    int natoms           /**< number of atoms */
    ) {
  int rc;
  double time_delta;

  if ((pm->msmflags & NL_MSM_COMPUTE_SPREC) == 0) {
    return NL_MSM_ERROR_PARAM;
  }

  /* store buffer pointers from caller */
  pm->felec_f = felec_f;
  pm->atom_f = atom_f;
  pm->numatoms = natoms;
  pm->uelec = 0;

  if (pm->msmflags & NL_MSM_COMPUTE_SHORT_RANGE) {
    wkf_timer_start(pm->timer);
    rc = NL_msm_compute_short_range_sprec(pm);
    if (rc) return rc;
    wkf_timer_stop(pm->timer);
    time_delta = wkf_timer_time(pm->timer);
    if (pm->report_timings) {
      printf("MSM short-range part:  %6.3f sec\n", time_delta);
    }
  }
  if (pm->msmflags & NL_MSM_COMPUTE_LONG_RANGE) {
    wkf_timer_start(pm->timer);
    rc = NL_msm_compute_long_range_sprec(pm);
    if (rc) return rc;
    wkf_timer_stop(pm->timer);
    time_delta = wkf_timer_time(pm->timer);
    if (pm->report_timings) {
      printf("MSM long-range part:   %6.3f sec\n", time_delta);
    }
  }
  *uelec_f += (float) pm->uelec;  /* add to caller's potential */

  return 0;
}


/** Order must be the same as APPROX enum in msm.h */
static const char *ApproxName[] = {
  "cubic",
  "quintic",
  "quintic2",
  "septic",
  "septic3",
  "nonic",
  "nonic4",
  "bspline",
};

int NL_msm_approx(const char *name) {
  int i;
  if (name) {
    for (i = 0;  i < NELEMS(ApproxName);  i++) {
      if (strcasecmp(name, ApproxName[i]) == 0) return i;
    }
  }
  return -1;
}

const char *NL_msm_approx_name(int approx) {
  if (approx >= 0 && approx < NELEMS(ApproxName)) {
    return ApproxName[approx];
  }
  else return NULL;
}


/** Order must be the same as SPLIT enum in msm.h */
static const char *SplitName[] = {
  "Taylor2",
  "Taylor3",
  "Taylor4",
  "Taylor5",
  "Taylor6",
  "Taylor7",
  "Taylor8",
  "Taylor1",
  "sigma2_3",
  "sigma3_5",
  "sigma4_6",
  "sigma4_7",
  "sigma5_8",
  "sigma5_9",
  "sigma6_9",
  "sigma6_10",
  "sigma6_11",
  "sigma7_11",
  "sigma7_12",
  "sigma7_13",
  "sigma8_12",
  "sigma8_13",
  "sigma8_14",
  "sigma8_15",
  "sigma2_6",
  "switch1_2",
  "switch3_4",
  "switch7_8",
};

int NL_msm_split(const char *name) {
  int i;
  if (name) {
    for (i = 0;  i < NELEMS(SplitName);  i++) {
      if (strcasecmp(name, SplitName[i]) == 0) return i;
    }
  }
  return -1;
}

const char *NL_msm_split_name(int split) {
  if (split >= 0 && split < NELEMS(SplitName)) {
    return SplitName[split];
  }
  else return NULL;
}
