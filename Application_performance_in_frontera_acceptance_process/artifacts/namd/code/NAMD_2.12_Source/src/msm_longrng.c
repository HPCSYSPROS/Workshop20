/** msm_longrng.c
 *
 * Compute the long-range part of MSM forces.
 *
 * The entire procedure can be depicted as an inverse V-cycle with bars
 * across each level to represent the short-range gridded calculations.  
 * The algorithm first performs anterpolation of the charge to the finest
 * level grid.  Each intermediate grid level, below the coarsest level at
 * the top, calculates a more slowly varying long-range part of the
 * interaction by restricting the charges to a grid with twice the spacing.
 * Each grid level also has a short-range cutoff part between charges that
 * contributes to the potentials.  Also summed to the potentials at each
 * intermediate grid level are the long-range contributions prolongated from
 * a grid with twice the spacing.  Finally, the forces are interpolated
 * (or approximated if the basis functions do not interpolate) from the
 * finest level grid of potentials.  The steps of the algorithm are 
 * detailed in section 2.3 of the thesis.
 *
 * The restriction and prolongation operations use the factored approach
 * having a work constant which is linear in the degree of the polynomial,
 * i.e. O(pM), where p is the degree and M is the number of grid points.  
 *
 * The algorithm below uses stencil widths expressed as a function of
 * the degree of the basis function polynomial Phi.  The constants for
 * restriction and prolongation stencils are stored in an array indexed
 * by the APPROX enum defined in msm.h header.  Calculations of the
 * grid stencils of the 1D Phi and its derivative, needed for interpolation
 * and anterpolation, are expressed as macros that use the APPROX enum
 * to select the polynomial pieces from a switch statement.  Periodic
 * boundaries are handled by wrapping around the edges of the grid,
 * as discussed in section 6.2 of the thesis.
 *
 * XXX The short-range gridded calculations contribute to the virial,
 * but it is not calculated.
 *
 * XXX The factored restriction and prolongation procedures are not
 * suitable for parallel and GPU computation due to the sequential 
 * dependence along each dimension of the grid.
 */

#include "msm_defn.h"

static int anterpolation(NL_Msm *);
static int interpolation(NL_Msm *);
static int restriction(NL_Msm *, int level);
static int prolongation(NL_Msm *, int level);
static int restriction_factored(NL_Msm *, int level);
static int prolongation_factored(NL_Msm *, int level);
static int gridcutoff(NL_Msm *, int level);


int NL_msm_compute_long_range(NL_Msm *msm) {
  double time_delta;
  int rc = 0;
  int level;
  int nlevels = msm->nlevels;
  int use_cuda_gridcutoff = (msm->msmflags & NL_MSM_COMPUTE_CUDA_GRID_CUTOFF);
  int fallback_cpu = (msm->msmflags & NL_MSM_COMPUTE_CUDA_FALL_BACK);
  int use_nonfactored = (msm->msmflags & NL_MSM_COMPUTE_NONFACTORED);

  wkf_timer_start(msm->timer_longrng);
  rc = anterpolation(msm);
  if (rc) return rc;
  wkf_timer_stop(msm->timer_longrng);
  time_delta = wkf_timer_time(msm->timer_longrng);
  if (msm->report_timings) {
    printf("MSM anterpolation:  %6.3f sec\n", time_delta);
  }

  wkf_timer_start(msm->timer_longrng);
  if (use_nonfactored) {
    for (level = 0;  level < nlevels - 1;  level++) {
      rc = restriction(msm, level);
      if (rc) return rc;
    }
  }
  else {
    for (level = 0;  level < nlevels - 1;  level++) {
      rc = restriction_factored(msm, level);
      if (rc) return rc;
    }
  }
  wkf_timer_stop(msm->timer_longrng);
  time_delta = wkf_timer_time(msm->timer_longrng);
  if (msm->report_timings) {
    printf("MSM restriction:    %6.3f sec\n", time_delta);
  }

  wkf_timer_start(msm->timer_longrng);
  if (use_cuda_gridcutoff && nlevels > 1) {
#ifdef NL_USE_CUDA
    if ((rc = NL_msm_cuda_condense_qgrids(msm)) != NL_MSM_SUCCESS ||
       	(rc = NL_msm_cuda_compute_gridcutoff(msm)) != NL_MSM_SUCCESS ||
       	(rc = NL_msm_cuda_expand_egrids(msm)) != NL_MSM_SUCCESS) {
      if (fallback_cpu) {
        printf("Falling back on CPU for grid cutoff computation\n");
        use_cuda_gridcutoff = 0;
      }
      else return rc;
    }
#else
    if (fallback_cpu) {
      printf("Falling back on CPU for grid cutoff computation\n");
      use_cuda_gridcutoff = 0;
    }
    else return NL_MSM_ERROR_SUPPORT;
#endif
  }

  if ( ! use_cuda_gridcutoff ) {
    for (level = 0;  level < nlevels - 1;  level++) {
      rc = gridcutoff(msm, level);
      if (rc) return rc;
    }
  }

  rc = gridcutoff(msm, level);  /* top level */
  if (rc) return rc;

  wkf_timer_stop(msm->timer_longrng);
  time_delta = wkf_timer_time(msm->timer_longrng);
  if (msm->report_timings) {
    printf("MSM grid cutoff:    %6.3f sec\n", time_delta);
  }

  wkf_timer_start(msm->timer_longrng);
  if (use_nonfactored) {
    for (level--;  level >= 0;  level--) {
      rc = prolongation(msm, level);
      if (rc) return rc;
    }
  }
  else {
    for (level--;  level >= 0;  level--) {
      rc = prolongation_factored(msm, level);
      if (rc) return rc;
    }
  }
  wkf_timer_stop(msm->timer_longrng);
  time_delta = wkf_timer_time(msm->timer_longrng);
  if (msm->report_timings) {
    printf("MSM prolongation:   %6.3f sec\n", time_delta);
  }

  wkf_timer_start(msm->timer_longrng);
  rc = interpolation(msm);
  if (rc) return rc;
  wkf_timer_stop(msm->timer_longrng);
  time_delta = wkf_timer_time(msm->timer_longrng);
  if (msm->report_timings) {
    printf("MSM interpolation:  %6.3f sec\n", time_delta);
  }

  return 0;
}


/** Approximation formulaes are up to degree 9 polynomials. */
enum { MAX_POLY_DEGREE = 9 };

/** Degree of polynomial basis function Phi.
 * Must be listed in same order as APPROX enum from msm.h */
static const int PolyDegree[NL_MSM_APPROX_END] = {
  3, 5, 5, 7, 7, 9, 9, 3,
};

/** The grid transfer stencils for factored restriction and prolongation.
 * Must be listed in same order as APPROX enum from msm.h */
static const double
PhiStencilFactored[NL_MSM_APPROX_END][2*MAX_POLY_DEGREE+1] = {
  /* cubic */
  {-1./16, 0, 9./16, 1, 9./16, 0, -1./16},

  /* quintic C1 */
  {3./256, 0, -25./256, 0, 75./128, 1, 75./128, 0, -25./256, 0, 3./256},

  /* quintic C2  (same as quintic C1) */
  {3./256, 0, -25./256, 0, 75./128, 1, 75./128, 0, -25./256, 0, 3./256},

  /* septic C1 */
  { -5./2048, 0, 49./2048, 0, -245./2048, 0, 1225./2048, 1, 1225./2048,
    0, -245./2048, 0, 49./2048, 0, -5./2048 },

  /* septic C3  (same as septic C3) */
  { -5./2048, 0, 49./2048, 0, -245./2048, 0, 1225./2048, 1, 1225./2048,
    0, -245./2048, 0, 49./2048, 0, -5./2048 },

  /* nonic C1 */
  { 35./65536, 0, -405./65536, 0, 567./16384, 0, -2205./16384, 0,
    19845./32768, 1, 19845./32768, 0, -2205./16384, 0, 567./16384, 0,
    -405./65536, 0, 35./65536 },

  /* nonic C4  (same as nonic C1) */
  { 35./65536, 0, -405./65536, 0, 567./16384, 0, -2205./16384, 0,
    19845./32768, 1, 19845./32768, 0, -2205./16384, 0, 567./16384, 0,
    -405./65536, 0, 35./65536 },

  /* bspline */
  { 1./48, 1./6, 23./48, 2./3, 23./48, 1./6, 1./48 },
};


/** Max stencil length is basically PolyDegree+2 for those approximations
 * that interpolate.  (We skip over zero in the complete stencils above.) */
enum { MAX_NSTENCIL = 11 };

/** The stencil array lengths below. */
static const int Nstencil[NL_MSM_APPROX_END] = {
  5, 7, 7, 9, 9, 11, 11, 7
};

/** Index offsets from the stencil-centered grid element, to get
 * to the correct contributing grid element. */
static const int IndexOffset[NL_MSM_APPROX_END][MAX_NSTENCIL] = {
  /* cubic */
  {-3, -1, 0, 1, 3},

  /* quintic C1 */
  {-5, -3, -1, 0, 1, 3, 5},

  /* quintic C2  (same as quintic C1) */
  {-5, -3, -1, 0, 1, 3, 5},

  /* septic C1 */
  {-7, -5, -3, -1, 0, 1, 3, 5, 7},

  /* septic C3  (same as septic C3) */
  {-7, -5, -3, -1, 0, 1, 3, 5, 7},

  /* nonic C1 */
  {-9, -7, -5, -3, -1, 0, 1, 3, 5, 7, 9},

  /* nonic C4  (same as nonic C1) */
  {-9, -7, -5, -3, -1, 0, 1, 3, 5, 7, 9},

  /* bspline */
  {-3, -2, -1, 0, 1, 2, 3},
};

/** The grid transfer stencils for the non-factored restriction and
 * prolongation procedures. */
static const double PhiStencil[NL_MSM_APPROX_END][MAX_NSTENCIL] = {
  /* cubic */
  {-1./16, 9./16, 1, 9./16, -1./16},

  /* quintic C1 */
  {3./256, -25./256, 75./128, 1, 75./128, -25./256, 3./256},

  /* quintic C2  (same as quintic C1) */
  {3./256, -25./256, 75./128, 1, 75./128, -25./256, 3./256},

  /* septic C1 */
  { -5./2048, 49./2048, -245./2048, 1225./2048, 1, 1225./2048,
    -245./2048, 49./2048, -5./2048 },

  /* septic C3  (same as septic C3) */
  { -5./2048, 49./2048, -245./2048, 1225./2048, 1, 1225./2048,
    -245./2048, 49./2048, -5./2048 },

  /* nonic C1 */
  { 35./65536, -405./65536, 567./16384, -2205./16384, 
    19845./32768, 1, 19845./32768, -2205./16384, 567./16384, 
    -405./65536, 35./65536 },

  /* nonic C4  (same as nonic C1) */
  { 35./65536, -405./65536, 567./16384, -2205./16384, 
    19845./32768, 1, 19845./32768, -2205./16384, 567./16384, 
    -405./65536, 35./65536 },

  /* bspline */
  { 1./48, 1./6, 23./48, 2./3, 23./48, 1./6, 1./48 },
};


/** Calculate the stencil of basis function values of Phi.
 * phi - stencil array (up to size MAX_POLY_DEGREE+1)
 * delta - normalized distance of atom from lowest grid point in stencil
 * approx - APPROX enum value from msm.h
 */
#define STENCIL_1D(_phi, _delta, _approx) \
  do { \
    double *phi = _phi; \
    double t = _delta; \
    switch (_approx) { \
      case NL_MSM_APPROX_CUBIC: \
        phi[0] = 0.5 * (1 - t) * (2 - t) * (2 - t); \
        t--; \
        phi[1] = (1 - t) * (1 + t - 1.5 * t * t); \
        t--; \
        phi[2] = (1 + t) * (1 - t - 1.5 * t * t); \
        t--; \
        phi[3] = 0.5 * (1 + t) * (2 + t) * (2 + t); \
        break; \
      case NL_MSM_APPROX_QUINTIC: \
        phi[0] = (1./24) * (1-t) * (2-t) * (3-t) * (3-t) * (4-t); \
        t--; \
        phi[1] = (1-t) * (2-t) * (3-t) * ((1./6) + t * (0.375 - (5./24)*t)); \
        t--; \
        phi[2] = (1-t*t) * (2-t) * (0.5 + t * (0.25 - (5./12)*t)); \
        t--; \
        phi[3] = (1-t*t) * (2+t) * (0.5 - t * (0.25 + (5./12)*t)); \
        t--; \
        phi[4] = (1+t) * (2+t) * (3+t) * ((1./6) - t * (0.375 + (5./24)*t)); \
        t--; \
        phi[5] = (1./24) * (1+t) * (2+t) * (3+t) * (3+t) * (4+t); \
        break; \
      case NL_MSM_APPROX_QUINTIC2: \
        phi[0] = (1./24) * (3-t) * (3-t) * (3-t) * (t-2) * (5*t-8); \
        t--; \
        phi[1] = (-1./24) * (2-t) * (t-1) * (-48+t*(153+t*(-114+t*25))); \
        t--; \
        phi[2] = (1./12) * (1-t) * (12+t*(12+t*(-3+t*(-38+t*25)))); \
        t--; \
        phi[3] = (1./12) * (1+t) * (12+t*(-12+t*(-3+t*(38+t*25)))); \
        t--; \
        phi[4] = (-1./24) * (2+t) * (t+1) * (48+t*(153+t*(114+t*25))); \
        t--; \
        phi[5] = (1./24) * (3+t) * (3+t) * (3+t) * (t+2) * (5*t+8); \
        break; \
      case NL_MSM_APPROX_SEPTIC: \
        phi[0] = (-1./720)*(t-1)*(t-2)*(t-3)*(t-4)*(t-4)*(t-5)*(t-6); \
        t--; \
        phi[1] = (1./720)*(t-1)*(t-2)*(t-3)*(t-4)*(t-5)*(-6+t*(-20+7*t)); \
        t--; \
        phi[2] = (-1./240)*(t*t-1)*(t-2)*(t-3)*(t-4)*(-10+t*(-12+7*t)); \
        t--; \
        phi[3] = (1./144)*(t*t-1)*(t*t-4)*(t-3)*(-12+t*(-4+7*t)); \
        t--; \
        phi[4] = (-1./144)*(t*t-1)*(t*t-4)*(t+3)*(-12+t*(4+7*t)); \
        t--; \
        phi[5] = (1./240)*(t*t-1)*(t+2)*(t+3)*(t+4)*(-10+t*(12+7*t)); \
        t--; \
        phi[6] = (-1./720)*(t+1)*(t+2)*(t+3)*(t+4)*(t+5)*(-6+t*(20+7*t)); \
        t--; \
        phi[7] = (1./720)*(t+1)*(t+2)*(t+3)*(t+4)*(t+4)*(t+5)*(t+6); \
        break; \
      case NL_MSM_APPROX_SEPTIC3: \
        phi[0] = (3632./5) + t*((-7456./5) + t*((58786./45) + t*(-633 \
                + t*((26383./144) + t*((-22807./720) + t*((727./240) \
                      + t*(-89./720))))))); \
        t--; \
        phi[1] = -440 + t*((25949./20) + t*((-117131./72) + t*((2247./2) \
                + t*((-66437./144) + t*((81109./720) + t*((-727./48) \
                      + t*(623./720))))))); \
        t--; \
        phi[2] = (138./5) + t*((-8617./60) + t*((12873./40) + t*((-791./2) \
                + t*((4557./16) + t*((-9583./80) + t*((2181./80) \
                      + t*(-623./240))))))); \
        t--; \
        phi[3] = 1 + t*t*((-49./36) + t*t*((-959./144) + t*((2569./144) \
                + t*((-727./48) + t*(623./144))))); \
        t--; \
        phi[4] = 1 + t*t*((-49./36) + t*t*((-959./144) + t*((-2569./144) \
                + t*((-727./48) + t*(-623./144))))); \
        t--; \
        phi[5] = (138./5) + t*((8617./60) + t*((12873./40) + t*((791./2) \
                + t*((4557./16) + t*((9583./80) + t*((2181./80) \
                      + t*(623./240))))))); \
        t--; \
        phi[6] = -440 + t*((-25949./20) + t*((-117131./72) + t*((-2247./2) \
                + t*((-66437./144) + t*((-81109./720) + t*((-727./48) \
                      + t*(-623./720))))))); \
        t--; \
        phi[7] = (3632./5) + t*((7456./5) + t*((58786./45) + t*(633 \
                + t*((26383./144) + t*((22807./720) + t*((727./240) \
                      + t*(89./720))))))); \
        break; \
      case NL_MSM_APPROX_NONIC: \
        phi[0] = (-1./40320)*(t-8)*(t-7)*(t-6)*(t-5)*(t-5)*(t-4)*(t-3)* \
          (t-2)*(t-1); \
        t--; \
        phi[1] = (1./40320)*(t-7)*(t-6)*(t-5)*(t-4)*(t-3)*(t-2)*(t-1)* \
          (-8+t*(-35+9*t)); \
        t--; \
        phi[2] = (-1./10080)*(t-6)*(t-5)*(t-4)*(t-3)*(t-2)*(t-1)*(t+1)* \
          (-14+t*(-25+9*t)); \
        t--; \
        phi[3] = (1./1440)*(t-5)*(t-4)*(t-3)*(t-2)*(t-1)*(t+1)*(t+2)* \
          (-6+t*(-5+3*t)); \
        t--; \
        phi[4] = (-1./2880)*(t-4)*(t-3)*(t-2)*(t-1)*(t+1)*(t+2)*(t+3)* \
          (-20+t*(-5+9*t)); \
        t--; \
        phi[5] = (1./2880)*(t-3)*(t-2)*(t-1)*(t+1)*(t+2)*(t+3)*(t+4)* \
          (-20+t*(5+9*t)); \
        t--; \
        phi[6] = (-1./1440)*(t-2)*(t-1)*(t+1)*(t+2)*(t+3)*(t+4)*(t+5)* \
          (-6+t*(5+3*t)); \
        t--; \
        phi[7] = (1./10080)*(t-1)*(t+1)*(t+2)*(t+3)*(t+4)*(t+5)*(t+6)* \
          (-14+t*(25+9*t)); \
        t--; \
        phi[8] = (-1./40320)*(t+1)*(t+2)*(t+3)*(t+4)*(t+5)*(t+6)*(t+7)* \
          (-8+t*(35+9*t)); \
        t--; \
        phi[9] = (1./40320)*(t+1)*(t+2)*(t+3)*(t+4)*(t+5)*(t+5)*(t+6)* \
          (t+7)*(t+8); \
        break; \
      case NL_MSM_APPROX_NONIC4: \
        phi[0] = 439375./7+t*(-64188125./504+t*(231125375./2016 \
              +t*(-17306975./288+t*(7761805./384+t*(-2895587./640 \
                    +t*(129391./192+t*(-259715./4032+t*(28909./8064 \
                          +t*(-3569./40320))))))))); \
        t--; \
        phi[1] = -56375+t*(8314091./56+t*(-49901303./288+t*(3763529./32 \
                +t*(-19648027./384+t*(9469163./640+t*(-545977./192 \
                      +t*(156927./448+t*(-28909./1152 \
                          +t*(3569./4480))))))))); \
        t--; \
        phi[2] = 68776./7+t*(-1038011./28+t*(31157515./504+t*(-956669./16 \
                +t*(3548009./96+t*(-2422263./160+t*(197255./48 \
                      +t*(-19959./28+t*(144545./2016 \
                          +t*(-3569./1120))))))))); \
        t--; \
        phi[3] = -154+t*(12757./12+t*(-230123./72+t*(264481./48 \
                +t*(-576499./96+t*(686147./160+t*(-96277./48 \
                      +t*(14221./24+t*(-28909./288+t*(3569./480))))))))); \
        t--; \
        phi[4] = 1+t*t*(-205./144+t*t*(91./192+t*(-6181./320 \
                +t*(6337./96+t*(-2745./32+t*(28909./576 \
                      +t*(-3569./320))))))); \
        t--; \
        phi[5] = 1+t*t*(-205./144+t*t*(91./192+t*(6181./320 \
                +t*(6337./96+t*(2745./32+t*(28909./576 \
                      +t*(3569./320))))))); \
        t--; \
        phi[6] = -154+t*(-12757./12+t*(-230123./72+t*(-264481./48 \
                +t*(-576499./96+t*(-686147./160+t*(-96277./48 \
                      +t*(-14221./24+t*(-28909./288+t*(-3569./480))))))))); \
        t--; \
        phi[7] = 68776./7+t*(1038011./28+t*(31157515./504+t*(956669./16 \
                +t*(3548009./96+t*(2422263./160+t*(197255./48 \
                      +t*(19959./28+t*(144545./2016+t*(3569./1120))))))))); \
        t--; \
        phi[8] = -56375+t*(-8314091./56+t*(-49901303./288+t*(-3763529./32 \
                +t*(-19648027./384+t*(-9469163./640+t*(-545977./192 \
                      +t*(-156927./448+t*(-28909./1152 \
                          +t*(-3569./4480))))))))); \
        t--; \
        phi[9] = 439375./7+t*(64188125./504+t*(231125375./2016 \
              +t*(17306975./288+t*(7761805./384+t*(2895587./640 \
                    +t*(129391./192+t*(259715./4032+t*(28909./8064 \
                          +t*(3569./40320))))))))); \
        break; \
      case NL_MSM_APPROX_BSPLINE: \
        phi[0] = (1./6) * (2-t) * (2-t) * (2-t); \
        t--; \
        phi[1] = (2./3) + t*t*(-1 + 0.5*t); \
        t--; \
        phi[2] = (2./3) + t*t*(-1 - 0.5*t); \
        t--; \
        phi[3] = (1./6) * (2+t) * (2+t) * (2+t); \
        break; \
      default: \
        return NL_MSM_ERROR_SUPPORT; \
    } \
  } while (0)
  /* closing ';' from use as function call */


/** Calculate the stencil of basis function and derivatives of (1/h)Phi.
 * dphi - stencil array (up to size MAX_POLY_DEGREE+1)
 * phi - stencil array (up to size MAX_POLY_DEGREE+1)
 * h_1 - 1/h, h is the grid spacing
 * delta - normalized distance of atom from lowest grid point in stencil
 * approx - APPROX enum value from msm.h
 */
#define D_STENCIL_1D(_dphi, _phi, _h_1, _delta, _approx) \
  do { \
    double *dphi = _dphi; \
    double *phi = _phi; \
    double h_1 = _h_1; \
    double t = _delta; \
    switch (_approx) { \
      case NL_MSM_APPROX_CUBIC: \
        phi[0] = 0.5 * (1 - t) * (2 - t) * (2 - t); \
        dphi[0] = (1.5 * t - 2) * (2 - t) * h_1; \
        t--; \
        phi[1] = (1 - t) * (1 + t - 1.5 * t * t); \
        dphi[1] = (-5 + 4.5 * t) * t * h_1; \
        t--; \
        phi[2] = (1 + t) * (1 - t - 1.5 * t * t); \
        dphi[2] = (-5 - 4.5 * t) * t * h_1; \
        t--; \
        phi[3] = 0.5 * (1 + t) * (2 + t) * (2 + t); \
        dphi[3] = (1.5 * t + 2) * (2 + t) * h_1; \
        break; \
      case NL_MSM_APPROX_QUINTIC: \
        phi[0] = (1./24) * (1-t) * (2-t) * (3-t) * (3-t) * (4-t); \
        dphi[0] = ((-1./24) * ((3-t) * (3-t) * (14 + t * (-14 + 3*t)) \
              + 2 * (1-t) * (2-t) * (3-t) * (4-t))) * h_1; \
        t--; \
        phi[1] = (1-t) * (2-t) * (3-t) * ((1./6) + t * (0.375 - (5./24)*t)); \
        dphi[1] = (-((1./6) + t * (0.375 - (5./24)*t)) * (11 + t * (-12 + 3*t))\
            + (1-t) * (2-t) * (3-t) * (0.375 - (5./12)*t)) * h_1; \
        t--; \
        phi[2] = (1-t*t) * (2-t) * (0.5 + t * (0.25 - (5./12)*t)); \
        dphi[2] = (-(0.5 + t * (0.25 - (5./12)*t)) * (1 + t * (4 - 3*t)) \
            + (1-t*t) * (2-t) * (0.25 - (5./6)*t)) * h_1; \
        t--; \
        phi[3] = (1-t*t) * (2+t) * (0.5 - t * (0.25 + (5./12)*t)); \
        dphi[3] = ((0.5 + t * (-0.25 - (5./12)*t)) * (1 + t * (-4 - 3*t)) \
            - (1-t*t) * (2+t) * (0.25 + (5./6)*t)) * h_1; \
        t--; \
        phi[4] = (1+t) * (2+t) * (3+t) * ((1./6) - t * (0.375 + (5./24)*t)); \
        dphi[4] = (((1./6) + t * (-0.375 - (5./24)*t)) * (11 + t * (12 + 3*t)) \
            - (1+t) * (2+t) * (3+t) * (0.375 + (5./12)*t)) * h_1; \
        t--; \
        phi[5] = (1./24) * (1+t) * (2+t) * (3+t) * (3+t) * (4+t); \
        dphi[5] = ((1./24) * ((3+t) * (3+t) * (14 + t * (14 + 3*t)) \
              + 2 * (1+t) * (2+t) * (3+t) * (4+t))) * h_1; \
        break; \
      case NL_MSM_APPROX_QUINTIC2: \
        phi[0] = (1./24) * (3-t) * (3-t) * (3-t) * (t-2) * (5*t-8); \
        dphi[0] = ((1./24) * (3-t) * (3-t) * ((3-t)*(5*t-8) - 3*(t-2)*(5*t-8) \
              + 5*(t-2)*(3-t))) * h_1; \
        t--; \
        phi[1] = (-1./24) * (2-t) * (t-1) * (-48+t*(153+t*(-114+t*25))); \
        dphi[1] = ((-1./24) * ((2-t)*(-48+t*(153+t*(-114+t*25))) - (t-1)* \
              (-48+t*(153+t*(-114+t*25))) \
              + (2-t)*(t-1)*(153+t*(-228+t*75)))) * h_1; \
        t--; \
        phi[2] = (1./12) * (1-t) * (12+t*(12+t*(-3+t*(-38+t*25)))); \
        dphi[2] = ((1./12) * (-(12+t*(12+t*(-3+t*(-38+t*25)))) \
              + (1-t)*(12+t*(-6+t*(-114+t*100))))) * h_1; \
        t--; \
        phi[3] = (1./12) * (1+t) * (12+t*(-12+t*(-3+t*(38+t*25)))); \
        dphi[3] = ((1./12) * ((12+t*(-12+t*(-3+t*(38+t*25)))) \
              + (1+t)*(-12+t*(-6+t*(114+t*100))))) * h_1; \
        t--; \
        phi[4] = (-1./24) * (2+t) * (t+1) * (48+t*(153+t*(114+t*25))); \
        dphi[4] = ((-1./24) * ((2+t)*(48+t*(153+t*(114+t*25))) + (t+1)* \
              (48+t*(153+t*(114+t*25))) \
              + (2+t)*(t+1)*(153+t*(228+t*75)))) * h_1; \
        t--; \
        phi[5] = (1./24) * (3+t) * (3+t) * (3+t) * (t+2) * (5*t+8); \
        dphi[5] = ((1./24) * (3+t) * (3+t) * ((3+t)*(5*t+8) + 3*(t+2)*(5*t+8) \
              + 5*(t+2)*(3+t))) * h_1; \
        break; \
      case NL_MSM_APPROX_SEPTIC: \
        phi[0] = (-1./720)*(t-1)*(t-2)*(t-3)*(t-4)*(t-4)*(t-5)*(t-6); \
        dphi[0] = (-1./720)*(t-4)*(-1944+t*(3644+t*(-2512+t*(807 \
                  +t*(-122+t*7))))) * h_1; \
        t--; \
        phi[1] = (1./720)*(t-1)*(t-2)*(t-3)*(t-4)*(t-5)*(-6+t*(-20+7*t)); \
        dphi[1] = (1./720)*(756+t*(-9940+t*(17724+t*(-12740+t*(4445 \
                    +t*(-750+t*49)))))) * h_1; \
        t--; \
        phi[2] = (-1./240)*(t*t-1)*(t-2)*(t-3)*(t-4)*(-10+t*(-12+7*t)); \
        dphi[2] = (-1./240)*(-28+t*(1260+t*(-756+t*(-1260+t*(1365 \
                    +t*(-450+t*49)))))) * h_1; \
        t--; \
        phi[3] = (1./144)*(t*t-1)*(t*t-4)*(t-3)*(-12+t*(-4+7*t)); \
        dphi[3] = (1./144)*t*(-560+t*(84+t*(644+t*(-175+t*(-150+t*49))))) * \
          h_1; \
        t--; \
        phi[4] = (-1./144)*(t*t-1)*(t*t-4)*(t+3)*(-12+t*(4+7*t)); \
        dphi[4] = (-1./144)*t*(560+t*(84+t*(-644+t*(-175+t*(150+t*49))))) * \
          h_1; \
        t--; \
        phi[5] = (1./240)*(t*t-1)*(t+2)*(t+3)*(t+4)*(-10+t*(12+7*t)); \
        dphi[5] = (1./240)*(-28+t*(-1260+t*(-756+t*(1260+t*(1365 \
                    +t*(450+t*49)))))) * h_1; \
        t--; \
        phi[6] = (-1./720)*(t+1)*(t+2)*(t+3)*(t+4)*(t+5)*(-6+t*(20+7*t)); \
        dphi[6] = (-1./720)*(756+t*(9940+t*(17724+t*(12740+t*(4445 \
                    +t*(750+t*49)))))) * h_1; \
        t--; \
        phi[7] = (1./720)*(t+1)*(t+2)*(t+3)*(t+4)*(t+4)*(t+5)*(t+6); \
        dphi[7] = (1./720)*(t+4)*(1944+t*(3644+t*(2512+t*(807 \
                  +t*(122+t*7))))) * h_1; \
        break; \
      case NL_MSM_APPROX_SEPTIC3: \
        phi[0] = (3632./5) + t*((-7456./5) + t*((58786./45) + t*(-633 \
                + t*((26383./144) + t*((-22807./720) + t*((727./240) \
                      + t*(-89./720))))))); \
        dphi[0] = ((-7456./5) + t*((117572./45) + t*(-1899 + t*((26383./36) \
                  + t*((-22807./144) + t*((727./40) + t*(-623./720)))))))*h_1; \
        t--; \
        phi[1] = -440 + t*((25949./20) + t*((-117131./72) + t*((2247./2) \
                + t*((-66437./144) + t*((81109./720) + t*((-727./48) \
                      + t*(623./720))))))); \
        dphi[1] = ((25949./20) + t*((-117131./36) + t*((6741./2) \
                + t*((-66437./36) + t*((81109./144) + t*((-727./8) \
                      + t*(4361./720))))))) * h_1; \
        t--; \
        phi[2] = (138./5) + t*((-8617./60) + t*((12873./40) + t*((-791./2) \
                + t*((4557./16) + t*((-9583./80) + t*((2181./80) \
                      + t*(-623./240))))))); \
        dphi[2] = ((-8617./60) + t*((12873./20) + t*((-2373./2) + t*((4557./4) \
                  + t*((-9583./16) + t*((6543./40) + t*(-4361./240)))))))*h_1; \
        t--; \
        phi[3] = 1 + t*t*((-49./36) + t*t*((-959./144) + t*((2569./144) \
                + t*((-727./48) + t*(623./144))))); \
        dphi[3] = (t*((-49./18) + t*t*((-959./36) + t*((12845./144) \
                  + t*((-727./8) + t*(4361./144)))))) * h_1; \
        t--; \
        phi[4] = 1 + t*t*((-49./36) + t*t*((-959./144) + t*((-2569./144) \
                + t*((-727./48) + t*(-623./144))))); \
        dphi[4] = (t*((-49./18) + t*t*((-959./36) + t*((-12845./144) \
                  + t*((-727./8) + t*(-4361./144)))))) * h_1; \
        t--; \
        phi[5] = (138./5) + t*((8617./60) + t*((12873./40) + t*((791./2) \
                + t*((4557./16) + t*((9583./80) + t*((2181./80) \
                      + t*(623./240))))))); \
        dphi[5] = ((8617./60) + t*((12873./20) + t*((2373./2) + t*((4557./4) \
                  + t*((9583./16) + t*((6543./40) + t*(4361./240)))))))*h_1; \
        t--; \
        phi[6] = -440 + t*((-25949./20) + t*((-117131./72) + t*((-2247./2) \
                + t*((-66437./144) + t*((-81109./720) + t*((-727./48) \
                      + t*(-623./720))))))); \
        dphi[6] = ((-25949./20) + t*((-117131./36) + t*((-6741./2) \
                + t*((-66437./36) + t*((-81109./144) + t*((-727./8) \
                      + t*(-4361./720))))))) * h_1; \
        t--; \
        phi[7] = (3632./5) + t*((7456./5) + t*((58786./45) + t*(633 \
                + t*((26383./144) + t*((22807./720) + t*((727./240) \
                      + t*(89./720))))))); \
        dphi[7] = ((7456./5) + t*((117572./45) + t*(1899 + t*((26383./36) \
                  + t*((22807./144) + t*((727./40) + t*(623./720)))))))*h_1; \
        break; \
      case NL_MSM_APPROX_NONIC: \
        phi[0] = (-1./40320)*(t-8)*(t-7)*(t-6)*(t-5)*(t-5)*(t-4)*(t-3)* \
          (t-2)*(t-1); \
        dphi[0] = (-1./40320)*(t-5)*(-117648+t*(256552+t*(-221416 \
                +t*(99340+t*(-25261+t*(3667+t*(-283+t*9)))))))*h_1; \
        t--; \
        phi[1] = (1./40320)*(t-7)*(t-6)*(t-5)*(t-4)*(t-3)*(t-2)*(t-1)* \
          (-8+t*(-35+9*t)); \
        dphi[1] = (1./40320)*(71856+t*(-795368+t*(1569240+t*(-1357692 \
                  +t*(634725+t*(-172116+t*(27090+t*(-2296+t*81))))))))*h_1; \
        t--; \
        phi[2] = (-1./10080)*(t-6)*(t-5)*(t-4)*(t-3)*(t-2)*(t-1)*(t+1)* \
          (-14+t*(-25+9*t)); \
        dphi[2] = (1./10080)*(3384+t*(-69080+t*(55026 \
                +t*(62580+t*(-99225+t*(51660+t*(-13104+t*(1640 \
                          +t*(-81)))))))))*h_1; \
        t--; \
        phi[3] = (1./1440)*(t-5)*(t-4)*(t-3)*(t-2)*(t-1)*(t+1)*(t+2)* \
          (-6+t*(-5+3*t)); \
        dphi[3] = (1./1440)*(72+t*(-6344+t*(2070 \
                +t*(7644+t*(-4725+t*(-828+t*(1260+t*(-328+t*27))))))))*h_1; \
        t--; \
        phi[4] = (-1./2880)*(t-4)*(t-3)*(t-2)*(t-1)*(t+1)*(t+2)*(t+3)* \
          (-20+t*(-5+9*t)); \
        dphi[4] = (-1./2880)*t*(10792+t*(-972+t*(-12516 \
                +t*(2205+t*(3924+t*(-882+t*(-328+t*81)))))))*h_1; \
        t--; \
        phi[5] = (1./2880)*(t-3)*(t-2)*(t-1)*(t+1)*(t+2)*(t+3)*(t+4)* \
          (-20+t*(5+9*t)); \
        dphi[5] = (1./2880)*t*(-10792+t*(-972+t*(12516 \
                +t*(2205+t*(-3924+t*(-882+t*(328+t*81)))))))*h_1; \
        t--; \
        phi[6] = (-1./1440)*(t-2)*(t-1)*(t+1)*(t+2)*(t+3)*(t+4)*(t+5)* \
          (-6+t*(5+3*t)); \
        dphi[6] = (1./1440)*(-72+t*(-6344+t*(-2070 \
                +t*(7644+t*(4725+t*(-828+t*(-1260+t*(-328+t*(-27)))))))))*h_1; \
        t--; \
        phi[7] = (1./10080)*(t-1)*(t+1)*(t+2)*(t+3)*(t+4)*(t+5)*(t+6)* \
          (-14+t*(25+9*t)); \
        dphi[7] = (1./10080)*(-3384+t*(-69080+t*(-55026 \
                +t*(62580+t*(99225+t*(51660+t*(13104+t*(1640+t*81))))))))*h_1; \
        t--; \
        phi[8] = (-1./40320)*(t+1)*(t+2)*(t+3)*(t+4)*(t+5)*(t+6)*(t+7)* \
          (-8+t*(35+9*t)); \
        dphi[8] = (-1./40320)*(71856+t*(795368+t*(1569240 \
                +t*(1357692+t*(634725+t*(172116+t*(27090+t*(2296 \
                          +t*81))))))))*h_1; \
        t--; \
        phi[9] = (1./40320)*(t+1)*(t+2)*(t+3)*(t+4)*(t+5)*(t+5)*(t+6)* \
          (t+7)*(t+8); \
        dphi[9] = (1./40320)*(t+5)*(117648+t*(256552+t*(221416 \
                +t*(99340+t*(25261+t*(3667+t*(283+t*9)))))))*h_1; \
        break; \
      case NL_MSM_APPROX_NONIC4: \
        phi[0] = 439375./7+t*(-64188125./504+t*(231125375./2016 \
              +t*(-17306975./288+t*(7761805./384+t*(-2895587./640 \
                    +t*(129391./192+t*(-259715./4032+t*(28909./8064 \
                          +t*(-3569./40320))))))))); \
        dphi[0] = (-64188125./504+t*(231125375./1008 \
              +t*(-17306975./96+t*(7761805./96+t*(-2895587./128 \
                    +t*(129391./32+t*(-259715./576+t*(28909./1008 \
                          +t*(-3569./4480)))))))))*h_1; \
        t--; \
        phi[1] = -56375+t*(8314091./56+t*(-49901303./288+t*(3763529./32 \
                +t*(-19648027./384+t*(9469163./640+t*(-545977./192 \
                      +t*(156927./448+t*(-28909./1152 \
                          +t*(3569./4480))))))))); \
        dphi[1] = (8314091./56+t*(-49901303./144+t*(11290587./32 \
                +t*(-19648027./96+t*(9469163./128+t*(-545977./32 \
                      +t*(156927./64+t*(-28909./144 \
                          +t*(32121./4480)))))))))*h_1; \
        t--; \
        phi[2] = 68776./7+t*(-1038011./28+t*(31157515./504+t*(-956669./16 \
                +t*(3548009./96+t*(-2422263./160+t*(197255./48 \
                      +t*(-19959./28+t*(144545./2016 \
                          +t*(-3569./1120))))))))); \
        dphi[2] = (-1038011./28+t*(31157515./252+t*(-2870007./16 \
                +t*(3548009./24+t*(-2422263./32+t*(197255./8 \
                      +t*(-19959./4+t*(144545./252 \
                          +t*(-32121./1120)))))))))*h_1; \
        t--; \
        phi[3] = -154+t*(12757./12+t*(-230123./72+t*(264481./48 \
                +t*(-576499./96+t*(686147./160+t*(-96277./48 \
                      +t*(14221./24+t*(-28909./288+t*(3569./480))))))))); \
        dphi[3] = (12757./12+t*(-230123./36+t*(264481./16 \
                +t*(-576499./24+t*(686147./32+t*(-96277./8 \
                      +t*(99547./24+t*(-28909./36 \
                          +t*(10707./160)))))))))*h_1; \
        t--; \
        phi[4] = 1+t*t*(-205./144+t*t*(91./192+t*(-6181./320 \
                +t*(6337./96+t*(-2745./32+t*(28909./576 \
                      +t*(-3569./320))))))); \
        dphi[4] = t*(-205./72+t*t*(91./48+t*(-6181./64 \
                +t*(6337./16+t*(-19215./32+t*(28909./72 \
                      +t*(-32121./320)))))))*h_1; \
        t--; \
        phi[5] = 1+t*t*(-205./144+t*t*(91./192+t*(6181./320 \
                +t*(6337./96+t*(2745./32+t*(28909./576 \
                      +t*(3569./320))))))); \
        dphi[5] = t*(-205./72+t*t*(91./48+t*(6181./64 \
                +t*(6337./16+t*(19215./32+t*(28909./72 \
                      +t*(32121./320)))))))*h_1; \
        t--; \
        phi[6] = -154+t*(-12757./12+t*(-230123./72+t*(-264481./48 \
                +t*(-576499./96+t*(-686147./160+t*(-96277./48 \
                      +t*(-14221./24+t*(-28909./288+t*(-3569./480))))))))); \
        dphi[6] = (-12757./12+t*(-230123./36+t*(-264481./16 \
                +t*(-576499./24+t*(-686147./32+t*(-96277./8 \
                      +t*(-99547./24+t*(-28909./36 \
                          +t*(-10707./160)))))))))*h_1; \
        t--; \
        phi[7] = 68776./7+t*(1038011./28+t*(31157515./504+t*(956669./16 \
                +t*(3548009./96+t*(2422263./160+t*(197255./48 \
                      +t*(19959./28+t*(144545./2016+t*(3569./1120))))))))); \
        dphi[7] = (1038011./28+t*(31157515./252+t*(2870007./16 \
                +t*(3548009./24+t*(2422263./32+t*(197255./8 \
                      +t*(19959./4+t*(144545./252 \
                          +t*(32121./1120)))))))))*h_1; \
        t--; \
        phi[8] = -56375+t*(-8314091./56+t*(-49901303./288+t*(-3763529./32 \
                +t*(-19648027./384+t*(-9469163./640+t*(-545977./192 \
                      +t*(-156927./448+t*(-28909./1152 \
                          +t*(-3569./4480))))))))); \
        dphi[8] = (-8314091./56+t*(-49901303./144+t*(-11290587./32 \
                +t*(-19648027./96+t*(-9469163./128+t*(-545977./32 \
                      +t*(-156927./64+t*(-28909./144 \
                          +t*(-32121./4480)))))))))*h_1; \
        t--; \
        phi[9] = 439375./7+t*(64188125./504+t*(231125375./2016 \
              +t*(17306975./288+t*(7761805./384+t*(2895587./640 \
                    +t*(129391./192+t*(259715./4032+t*(28909./8064 \
                          +t*(3569./40320))))))))); \
        dphi[9] = (64188125./504+t*(231125375./1008 \
              +t*(17306975./96+t*(7761805./96+t*(2895587./128 \
                    +t*(129391./32+t*(259715./576+t*(28909./1008 \
                          +t*(3569./4480)))))))))*h_1; \
        break; \
      case NL_MSM_APPROX_BSPLINE: \
        phi[0] = (1./6) * (2-t) * (2-t) * (2-t); \
        dphi[0] = -0.5 * (2-t) * (2-t) * h_1; \
        t--; \
        phi[1] = (2./3) + t*t*(-1 + 0.5*t); \
        dphi[1] = t*(-2 + 1.5*t) * h_1; \
        t--; \
        phi[2] = (2./3) + t*t*(-1 - 0.5*t); \
        dphi[2] = t*(-2 - 1.5*t) * h_1; \
        t--; \
        phi[3] = (1./6) * (2+t) * (2+t) * (2+t); \
        dphi[3] = 0.5 * (2+t) * (2+t) * h_1; \
        break; \
      default: \
        return NL_MSM_ERROR_SUPPORT; \
    } \
  } while (0)
  /* closing ';' from use as function call */


int anterpolation(NL_Msm *msm)
{
  const double *atom = msm->atom;
  const int natoms = msm->numatoms;

  double xphi[MAX_POLY_DEGREE+1];  /* Phi stencil along x-dimension */
  double yphi[MAX_POLY_DEGREE+1];  /* Phi stencil along y-dimension */
  double zphi[MAX_POLY_DEGREE+1];  /* Phi stencil along z-dimension */
  double rx_hx, ry_hy, rz_hz;      /* normalized distance from atom to origin */
#ifndef TEST_INLINING
  double delta;                    /* normalized distance to stencil point */
#else
  double t;                    /* normalized distance to stencil point */
#endif
  double ck, cjk;
  const double hx_1 = 1/msm->hx;
  const double hy_1 = 1/msm->hy;
  const double hz_1 = 1/msm->hz;
  const double xm0 = msm->gx;      /* grid origin x */
  const double ym0 = msm->gy;      /* grid origin y */
  const double zm0 = msm->gz;      /* grid origin z */
  double q;

  NL_Msmgrid_double *qhgrid = &(msm->qh[0]);
  double *qh = qhgrid->data;
  const int ni = qhgrid->ni;
  const int nj = qhgrid->nj;
  const int nk = qhgrid->nk;
  const int ia = qhgrid->i0;
  const int ib = ia + ni - 1;
  const int ja = qhgrid->j0;
  const int jb = ja + nj - 1;
  const int ka = qhgrid->k0;
  const int kb = ka + nk - 1;

  const int ispx = (msm->msmflags & NL_MSM_PERIODIC_VEC1);
  const int ispy = (msm->msmflags & NL_MSM_PERIODIC_VEC2);
  const int ispz = (msm->msmflags & NL_MSM_PERIODIC_VEC3);
  int iswithin;

  int n, i, j, k, ilo, jlo, klo, koff;
  int jkoff, index;

  const int approx = msm->approx;
  const int s_size = PolyDegree[approx] + 1;         /* stencil size */
  const int s_edge = (PolyDegree[approx] - 1) >> 1;  /* stencil "edge" size */

  GRID_ZERO(qhgrid);

  for (n = 0;  n < natoms;  n++) {

    /* atomic charge */
    q = atom[4*n + 3];
    if (0==q) continue;

    /* distance between atom and origin measured in grid points */
    rx_hx = (atom[4*n    ] - xm0) * hx_1;
    ry_hy = (atom[4*n + 1] - ym0) * hy_1;
    rz_hz = (atom[4*n + 2] - zm0) * hz_1;

    /* find smallest numbered grid point in stencil */
    ilo = (int) floor(rx_hx) - s_edge;
    jlo = (int) floor(ry_hy) - s_edge;
    klo = (int) floor(rz_hz) - s_edge;

    /* calculate Phi stencils along each dimension */
#ifndef TEST_INLINING
    delta = rx_hx - (double) ilo;
    STENCIL_1D(xphi, delta, approx);
    delta = ry_hy - (double) jlo;
    STENCIL_1D(yphi, delta, approx);
    delta = rz_hz - (double) klo;
    STENCIL_1D(zphi, delta, approx);
#else
    t = rx_hx - (double) ilo;
        xphi[0] = 0.5 * (1 - t) * (2 - t) * (2 - t); \
        t--; \
        xphi[1] = (1 - t) * (1 + t - 1.5 * t * t); \
        t--; \
        xphi[2] = (1 + t) * (1 - t - 1.5 * t * t); \
        t--; \
        xphi[3] = 0.5 * (1 + t) * (2 + t) * (2 + t); \

    t = ry_hy - (double) jlo;
        yphi[0] = 0.5 * (1 - t) * (2 - t) * (2 - t); \
        t--; \
        yphi[1] = (1 - t) * (1 + t - 1.5 * t * t); \
        t--; \
        yphi[2] = (1 + t) * (1 - t - 1.5 * t * t); \
        t--; \
        yphi[3] = 0.5 * (1 + t) * (2 + t) * (2 + t); \

    t = rz_hz - (double) klo;
        zphi[0] = 0.5 * (1 - t) * (2 - t) * (2 - t); \
        t--; \
        zphi[1] = (1 - t) * (1 + t - 1.5 * t * t); \
        t--; \
        zphi[2] = (1 + t) * (1 - t - 1.5 * t * t); \
        t--; \
        zphi[3] = 0.5 * (1 + t) * (2 + t) * (2 + t); \

#endif

    /* test to see if stencil is within edges of grid */
    iswithin = ( ia <= ilo && (ilo+(s_size-1)) <= ib &&
                 ja <= jlo && (jlo+(s_size-1)) <= jb &&
                 ka <= klo && (klo+(s_size-1)) <= kb );

    if ( iswithin ) {  /* no wrapping needed */

      /* determine charge on cube of grid points around atom */
      for (k = 0;  k < s_size;  k++) {
        koff = (k + klo) * nj;
        ck = zphi[k] * q;
        for (j = 0;  j < s_size;  j++) {
          jkoff = (koff + (j + jlo)) * ni;
          cjk = yphi[j] * ck;
          for (i = 0;  i < s_size;  i++) {
            index = jkoff + (i + ilo);
            GRID_INDEX_CHECK(qhgrid, i+ilo, j+jlo, k+klo);
            ASSERT(GRID_INDEX(qhgrid, i+ilo, j+jlo, k+klo) == index);
            qh[index] += xphi[i] * cjk;
          }
        }
      }
    } /* if */

    else {  /* requires wrapping around grid */
      int ip, jp, kp;

      /* adjust ilo, jlo, klo so they are within grid indexing */
      if (ispx) {
        if      (ilo < ia) do { ilo += ni; } while (ilo < ia);
        else if (ilo > ib) do { ilo -= ni; } while (ilo > ib);
      }
      else if (ilo < ia || (ilo+(s_size-1)) > ib) {
        return NL_MSM_ERROR_RANGE;
      }

      if (ispy) {
        if      (jlo < ja) do { jlo += nj; } while (jlo < ja);
        else if (jlo > jb) do { jlo -= nj; } while (jlo > jb);
      }
      else if (jlo < ja || (jlo+(s_size-1)) > jb) {
        return NL_MSM_ERROR_RANGE;
      }

      if (ispz) {
        if      (klo < ka) do { klo += nk; } while (klo < ka);
        else if (klo > kb) do { klo -= nk; } while (klo > kb);
      }
      else if (klo < ka || (klo+(s_size-1)) > kb) {
        return NL_MSM_ERROR_RANGE;
      }

      /* determine charge on cube of grid points around atom, with wrapping */
      for (k = 0, kp = klo;  k < s_size;  k++, kp++) {
        if (kp > kb) kp = ka;  /* wrap stencil around grid */
        koff = kp * nj;
        ck = zphi[k] * q;
        for (j = 0, jp = jlo;  j < s_size;  j++, jp++) {
          if (jp > jb) jp = ja;  /* wrap stencil around grid */
          jkoff = (koff + jp) * ni;
          cjk = yphi[j] * ck;
          for (i = 0, ip = ilo;  i < s_size;  i++, ip++) {
            if (ip > ib) ip = ia;  /* wrap stencil around grid */
            index = jkoff + ip;
            GRID_INDEX_CHECK(qhgrid, ip, jp, kp);
            ASSERT(GRID_INDEX(qhgrid, ip, jp, kp) == index);
            qh[index] += xphi[i] * cjk;
          }
        }
      }
    } /* else */

  } /* end loop over atoms */

  return NL_MSM_SUCCESS;
} /* anterpolation */


int interpolation(NL_Msm *msm)
{
  double *f = msm->felec;
  const double *atom = msm->atom;
  const int natoms = msm->numatoms;

  double xphi[MAX_POLY_DEGREE+1];  /* Phi stencil along x-dimension */
  double yphi[MAX_POLY_DEGREE+1];  /* Phi stencil along y-dimension */
  double zphi[MAX_POLY_DEGREE+1];  /* Phi stencil along z-dimension */
  double dxphi[MAX_POLY_DEGREE+1]; /* derivative of Phi along x-dimension */
  double dyphi[MAX_POLY_DEGREE+1]; /* derivative of Phi along y-dimension */
  double dzphi[MAX_POLY_DEGREE+1]; /* derivative of Phi along z-dimension */
  double rx_hx, ry_hy, rz_hz;      /* normalized distance from atom to origin */
#ifndef TEST_INLINING
  double delta;                    /* normalized distance to stencil point */
#else
  double t;                    /* normalized distance to stencil point */
#endif
  const double hx_1 = 1/msm->hx;
  const double hy_1 = 1/msm->hy;
  const double hz_1 = 1/msm->hz;
  const double xm0 = msm->gx;      /* grid origin x */
  const double ym0 = msm->gy;      /* grid origin y */
  const double zm0 = msm->gz;      /* grid origin z */
  double q;
  double u = 0;

  const NL_Msmgrid_double *ehgrid = &(msm->eh[0]);
  const double *eh = ehgrid->data;
  const double *ebuf = ehgrid->buffer;
  const NL_Msmgrid_double *qhgrid = &(msm->qh[0]);
  const double *qbuf = qhgrid->buffer;
  const int ni = ehgrid->ni;
  const int nj = ehgrid->nj;
  const int nk = ehgrid->nk;
  const int ia = ehgrid->i0;
  const int ib = ia + ni - 1;
  const int ja = ehgrid->j0;
  const int jb = ja + nj - 1;
  const int ka = ehgrid->k0;
  const int kb = ka + nk - 1;

  const int ispx = (msm->msmflags & NL_MSM_PERIODIC_VEC1);
  const int ispy = (msm->msmflags & NL_MSM_PERIODIC_VEC2);
  const int ispz = (msm->msmflags & NL_MSM_PERIODIC_VEC3);
  int iswithin;

  double fx, fy, fz, cx, cy, cz;

  int n, i, j, k, ilo, jlo, klo, koff;
  long jkoff, index;
  const int nn = (ni*nj) * nk;

  const int approx = msm->approx;
  const int s_size = PolyDegree[approx] + 1;         /* stencil size */
  const int s_edge = (PolyDegree[approx] - 1) >> 1;  /* stencil "edge" size */

  for (n = 0;  n < natoms;  n++) {

    /* atomic charge */
    q = atom[4*n + 3];
    if (0==q) continue;

    /* distance between atom and origin measured in grid points */
    rx_hx = (atom[4*n    ] - xm0) * hx_1;
    ry_hy = (atom[4*n + 1] - ym0) * hy_1;
    rz_hz = (atom[4*n + 2] - zm0) * hz_1;

    /* find smallest numbered grid point in stencil */
    ilo = (int) floor(rx_hx) - s_edge;
    jlo = (int) floor(ry_hy) - s_edge;
    klo = (int) floor(rz_hz) - s_edge;

    /* calculate Phi stencils and derivatives along each dimension */
#ifndef TEST_INLINING
    delta = rx_hx - (double) ilo;
    //STENCIL_1D(xphi, delta, approx);
    D_STENCIL_1D(dxphi, xphi, hx_1, delta, approx);
    delta = ry_hy - (double) jlo;
    //STENCIL_1D(yphi, delta, approx);
    D_STENCIL_1D(dyphi, yphi, hy_1, delta, approx);
    delta = rz_hz - (double) klo;
    //STENCIL_1D(zphi, delta, approx);
    D_STENCIL_1D(dzphi, zphi, hz_1, delta, approx);
#else
    t = rx_hx - (double) ilo;
        xphi[0] = 0.5 * (1 - t) * (2 - t) * (2 - t); \
        dxphi[0] = (1.5 * t - 2) * (2 - t) * hx_1; \
        t--; \
        xphi[1] = (1 - t) * (1 + t - 1.5 * t * t); \
        dxphi[1] = (-5 + 4.5 * t) * t * hx_1; \
        t--; \
        xphi[2] = (1 + t) * (1 - t - 1.5 * t * t); \
        dxphi[2] = (-5 - 4.5 * t) * t * hx_1; \
        t--; \
        xphi[3] = 0.5 * (1 + t) * (2 + t) * (2 + t); \
        dxphi[3] = (1.5 * t + 2) * (2 + t) * hx_1; \

    t = ry_hy - (double) jlo;
        yphi[0] = 0.5 * (1 - t) * (2 - t) * (2 - t); \
        dyphi[0] = (1.5 * t - 2) * (2 - t) * hy_1; \
        t--; \
        yphi[1] = (1 - t) * (1 + t - 1.5 * t * t); \
        dyphi[1] = (-5 + 4.5 * t) * t * hy_1; \
        t--; \
        yphi[2] = (1 + t) * (1 - t - 1.5 * t * t); \
        dyphi[2] = (-5 - 4.5 * t) * t * hy_1; \
        t--; \
        yphi[3] = 0.5 * (1 + t) * (2 + t) * (2 + t); \
        dyphi[3] = (1.5 * t + 2) * (2 + t) * hy_1; \

    t = rz_hz - (double) klo;
        zphi[0] = 0.5 * (1 - t) * (2 - t) * (2 - t); \
        dzphi[0] = (1.5 * t - 2) * (2 - t) * hz_1; \
        t--; \
        zphi[1] = (1 - t) * (1 + t - 1.5 * t * t); \
        dzphi[1] = (-5 + 4.5 * t) * t * hz_1; \
        t--; \
        zphi[2] = (1 + t) * (1 - t - 1.5 * t * t); \
        dzphi[2] = (-5 - 4.5 * t) * t * hz_1; \
        t--; \
        zphi[3] = 0.5 * (1 + t) * (2 + t) * (2 + t); \
        dzphi[3] = (1.5 * t + 2) * (2 + t) * hz_1; \

#endif

    /* test to see if stencil is within edges of grid */
    iswithin = ( ia <= ilo && (ilo+(s_size-1)) <= ib &&
                 ja <= jlo && (jlo+(s_size-1)) <= jb &&
                 ka <= klo && (klo+(s_size-1)) <= kb );

    if ( iswithin ) {  /* no wrapping needed */

      /* determine force on atom from cube of grid point potentials */
      fx = fy = fz = 0;
      for (k = 0;  k < s_size;  k++) {
        koff = (k + klo) * nj;
        for (j = 0;  j < s_size;  j++) {
          jkoff = (koff + (j + jlo)) * ni;
          cx = yphi[j] * zphi[k];
          cy = dyphi[j] * zphi[k];
          cz = yphi[j] * dzphi[k];
          for (i = 0;  i < s_size;  i++) {
            index = jkoff + (i + ilo);
            GRID_INDEX_CHECK(ehgrid, i+ilo, j+jlo, k+klo);
            ASSERT(GRID_INDEX(ehgrid, i+ilo, j+jlo, k+klo) == index);
            fx += eh[index] * dxphi[i] * cx;
            fy += eh[index] * xphi[i] * cy;
            fz += eh[index] * xphi[i] * cz;
          }
        }
      }
    } /* if */

    else {  /* requires wrapping around grid */
      int ip, jp, kp;

      /* adjust ilo, jlo, klo so they are within grid indexing */
      if (ispx) {
        if      (ilo < ia) do { ilo += ni; } while (ilo < ia);
        else if (ilo > ib) do { ilo -= ni; } while (ilo > ib);
      }
      else if (ilo < ia || (ilo+(s_size-1)) > ib) {
        return NL_MSM_ERROR_RANGE;
      }

      if (ispy) {
        if      (jlo < ja) do { jlo += nj; } while (jlo < ja);
        else if (jlo > jb) do { jlo -= nj; } while (jlo > jb);
      }
      else if (jlo < ja || (jlo+(s_size-1)) > jb) {
        return NL_MSM_ERROR_RANGE;
      }

      if (ispz) {
        if      (klo < ka) do { klo += nk; } while (klo < ka);
        else if (klo > kb) do { klo -= nk; } while (klo > kb);
      }
      else if (klo < ka || (klo+(s_size-1)) > kb) {
        return NL_MSM_ERROR_RANGE;
      }

      /* determine force on atom from cube of grid point potentials, wrapping */
      fx = fy = fz = 0;
      for (k = 0, kp = klo;  k < s_size;  k++, kp++) {
        if (kp > kb) kp = ka;  /* wrap stencil around grid */
        koff = kp * nj;
        for (j = 0, jp = jlo;  j < s_size;  j++, jp++) {
          if (jp > jb) jp = ja;  /* wrap stencil around grid */
          jkoff = (koff + jp) * ni;
          cx = yphi[j] * zphi[k];
          cy = dyphi[j] * zphi[k];
          cz = yphi[j] * dzphi[k];
          for (i = 0, ip = ilo;  i < s_size;  i++, ip++) {
            if (ip > ib) ip = ia;  /* wrap stencil around grid */
            index = jkoff + ip;
            GRID_INDEX_CHECK(ehgrid, ip, jp, kp);
            ASSERT(GRID_INDEX(ehgrid, ip, jp, kp) == index);
            fx += eh[index] * dxphi[i] * cx;
            fy += eh[index] * xphi[i] * cy;
            fz += eh[index] * xphi[i] * cz;
          }
        }
      }
    } /* else */

    /* update force */
    f[3*n  ] -= q * fx;
    f[3*n+1] -= q * fy;
    f[3*n+2] -= q * fz;
//    printf("force[%d] = %g %g %g\n", n, -q*fx, -q*fy, -q*fz);

  } /* end loop over atoms */

  /* compute potential, subtract self potential */
  u = 0;
  if (1) {
    for (n = 0;  n < natoms;  n++) {
      double q = atom[4*n + 3];
      u += q * q;
    }
    u *= -msm->gzero;
  }
  for (index = 0;  index < nn;  index++) {
    u += qbuf[index] * ebuf[index];
  }
  msm->uelec += 0.5 * u;
//  printf("MSM long-range potential:  %g\n", 0.5 * u);

  return NL_MSM_SUCCESS;
} /* interpolation() */


int restriction_factored(NL_Msm *msm, int level) {
  /* grids of charge, finer grid and coarser grid */
  const NL_Msmgrid_double *qhgrid = &(msm->qh[level]);
  const double *qh = qhgrid->data;
  NL_Msmgrid_double *q2hgrid = &(msm->qh[level+1]);
  double *q2h = q2hgrid->data;

  /* finer grid index ranges and dimensions */
  const int ia1 = qhgrid->i0;             /* lowest x-index */
  const int ib1 = ia1 + qhgrid->ni - 1;   /* highest x-index */
  const int ja1 = qhgrid->j0;             /* lowest y-index */
  const int jb1 = ja1 + qhgrid->nj - 1;   /* highest y-index */
  const int ka1 = qhgrid->k0;             /* lowest z-index */
  const int kb1 = ka1 + qhgrid->nk - 1;   /* highest z-index */
  const int ni1 = qhgrid->ni;             /* length along x-dim */
  const int nj1 = qhgrid->nj;             /* length along y-dim */
  const int nk1 = qhgrid->nk;             /* length along z-dim */

  /* coarser grid index ranges and dimensions */
  const int ia2 = q2hgrid->i0;            /* lowest x-index */
  const int ib2 = ia2 + q2hgrid->ni - 1;  /* highest x-index */
  const int ja2 = q2hgrid->j0;            /* lowest y-index */
  const int jb2 = ja2 + q2hgrid->nj - 1;  /* highest y-index */
  const int ka2 = q2hgrid->k0;            /* lowest z-index */
  const int kb2 = ka2 + q2hgrid->nk - 1;  /* highest z-index */
  const int nrow_q2 = q2hgrid->ni;
  const int nstride_q2 = nrow_q2 * q2hgrid->nj;

  const int ispx = (msm->msmflags & NL_MSM_PERIODIC_VEC1);
  const int ispy = (msm->msmflags & NL_MSM_PERIODIC_VEC2);
  const int ispz = (msm->msmflags & NL_MSM_PERIODIC_VEC3);

  /* set buffer using indexing offset, so that indexing matches qh grid */
  double *qzd = msm->lzd + (-ka1);
  double *qyzd = msm->lyzd + (-ka1*nj1 + -ja1);
  double qsum;

  const double *phi = NULL;

  int i2, j2, k2;
  int im, jm, km;
  int i, j, k;
  int index_plane_q2, index_q2;
  int index_jk, offset_k;
  int offset;

  const double *phi_factored = PhiStencilFactored[msm->approx];
  const int r_stencil = PolyDegree[msm->approx];  /* "radius" of stencil */
  const int diam_stencil = 2*r_stencil + 1;       /* "diameter" of stencil */

  for (i2 = ia2;  i2 <= ib2;  i2++) {

    for (k = ka1;  k <= kb1;  k++) {
      offset_k = k * nj1;

      for (j = ja1;  j <= jb1;  j++) {
        index_jk = offset_k + j;
        offset = index_jk * ni1;
        im = (i2 << 1);  /* = 2*i2 */
        qsum = 0;
        if ( ! ispx ) {  /* nonperiodic */
          int lower = im - r_stencil;
          int upper = im + r_stencil;
          if (lower < ia1) lower = ia1;
          if (upper > ib1) upper = ib1;  /* clip edges */
          phi = phi_factored + r_stencil;  /* center of stencil */
          for (i = lower;  i <= upper;  i++) {
            qsum += phi[i-im] * qh[offset + i];
          }
        }
        else {  /* periodic */
          int ip = im - r_stencil;  /* index at left end of stencil */
          if (ip < ia1) do { ip += ni1; } while (ip < ia1);  /* start inside */
          phi = phi_factored;  /* left end of stencil */
          for (i = 0;  i < diam_stencil;  i++, ip++) {
            if (ip > ib1) ip = ia1;  /* wrap around edge of grid */
            qsum += phi[i] * qh[offset + ip];
          }
        }
        qyzd[index_jk] = qsum;
      } /* for j */

    } /* for k */

    for (j2 = ja2;  j2 <= jb2;  j2++) {
      index_plane_q2 = j2 * nrow_q2 + i2;

      for (k = ka1;  k <= kb1;  k++) {
        offset = k * nj1;
        jm = (j2 << 1);  /* = 2*j2 */
        qsum = 0;
        if ( ! ispy ) {  /* nonperiodic */
          int lower = jm - r_stencil;
          int upper = jm + r_stencil;
          if (lower < ja1) lower = ja1;
          if (upper > jb1) upper = jb1;  /* clip edges */
          phi = phi_factored + r_stencil;  /* center of stencil */
          for (j = lower;  j <= upper;  j++) {
            qsum += phi[j-jm] * qyzd[offset + j];
          }
        }
        else {  /* periodic */
          int jp = jm - r_stencil;  /* index at left end of stencil */
          if (jp < ja1) do { jp += nj1; } while (jp < ja1);  /* start inside */
          phi = phi_factored;  /* left end of stencil */
          for (j = 0;  j < diam_stencil;  j++, jp++) {
            if (jp > jb1) jp = ja1;  /* wrap around edge of grid */
            qsum += phi[j] * qyzd[offset + jp];
          }
        }
        qzd[k] = qsum;
      } /* for k */

      for (k2 = ka2;  k2 <= kb2;  k2++) {
        index_q2 = k2 * nstride_q2 + index_plane_q2;
        km = (k2 << 1);  /* = 2*k2 */
        qsum = 0;
        if ( ! ispz ) {  /* nonperiodic */
          int lower = km - r_stencil;
          int upper = km + r_stencil;
          if (lower < ka1) lower = ka1;
          if (upper > kb1) upper = kb1;  /* clip edges */
          phi = phi_factored + r_stencil;  /* center of stencil */
          for (k = lower;  k <= upper;  k++) {
            qsum += phi[k-km] * qzd[k];
          }
        }
        else {  /* periodic */
          int kp = km - r_stencil;  /* index at left end of stencil */
          if (kp < ka1) do { kp += nk1; } while (kp < ka1);  /* start inside */
          phi = phi_factored;  /* left end of stencil */
          for (k = 0;  k < diam_stencil;  k++, kp++) {
            if (kp > kb1) kp = ka1;  /* wrap around edge of grid */
            qsum += phi[k] * qzd[kp];
          }
        }
        q2h[index_q2] = qsum;
      } /* for k2 */

    } /* for j2 */

  } /* for i2 */

  return NL_MSM_SUCCESS;
} /* restriction_factored */


int prolongation_factored(NL_Msm *msm, int level) {
  /* grids of potential, finer grid and coarser grid */
  NL_Msmgrid_double *ehgrid = &(msm->eh[level]);
  double *eh = ehgrid->data;
  const NL_Msmgrid_double *e2hgrid = &(msm->eh[level+1]);
  const double *e2h = e2hgrid->data;

  /* finer grid index ranges and dimensions */
  const int ia1 = ehgrid->i0;             /* lowest x-index */
  const int ib1 = ia1 + ehgrid->ni - 1;   /* highest x-index */
  const int ja1 = ehgrid->j0;             /* lowest y-index */
  const int jb1 = ja1 + ehgrid->nj - 1;   /* highest y-index */
  const int ka1 = ehgrid->k0;             /* lowest z-index */
  const int kb1 = ka1 + ehgrid->nk - 1;   /* highest z-index */
  const int ni1 = ehgrid->ni;             /* length along x-dim */
  const int nj1 = ehgrid->nj;             /* length along y-dim */
  const int nk1 = ehgrid->nk;             /* length along z-dim */

  /* coarser grid index ranges and dimensions */
  const int ia2 = e2hgrid->i0;            /* lowest x-index */
  const int ib2 = ia2 + e2hgrid->ni - 1;  /* highest x-index */
  const int ja2 = e2hgrid->j0;            /* lowest y-index */
  const int jb2 = ja2 + e2hgrid->nj - 1;  /* highest y-index */
  const int ka2 = e2hgrid->k0;            /* lowest z-index */
  const int kb2 = ka2 + e2hgrid->nk - 1;  /* highest z-index */
  const int nrow_e2 = e2hgrid->ni;
  const int nstride_e2 = nrow_e2 * e2hgrid->nj;

  const int ispx = (msm->msmflags & NL_MSM_PERIODIC_VEC1);
  const int ispy = (msm->msmflags & NL_MSM_PERIODIC_VEC2);
  const int ispz = (msm->msmflags & NL_MSM_PERIODIC_VEC3);

  /* set buffer using indexing offset, so that indexing matches eh grid */
  double *ezd = msm->lzd + (-ka1);
  double *eyzd = msm->lyzd + (-ka1*nj1 + -ja1);

  const size_t size_lzd = nk1 * sizeof(double);
  const size_t size_lyzd = nj1 * nk1 * sizeof(double);

  const double *phi = NULL;

  int i2, j2, k2;
  int im, jm, km;
  int i, j, k;
  int index_plane_e2, index_e2;
  int index_jk, offset_k;
  int offset;

  const double *phi_factored = PhiStencilFactored[msm->approx];
  const int r_stencil = PolyDegree[msm->approx];  /* "radius" of stencil */
  const int diam_stencil = 2*r_stencil + 1;       /* "diameter" of stencil */

  for (i2 = ia2;  i2 <= ib2;  i2++) {
    memset(msm->lyzd, 0, size_lyzd);

    for (j2 = ja2;  j2 <= jb2;  j2++) {
      memset(msm->lzd, 0, size_lzd);
      index_plane_e2 = j2 * nrow_e2 + i2;

      for (k2 = ka2;  k2 <= kb2;  k2++) {
        index_e2 = k2 * nstride_e2 + index_plane_e2;
        km = (k2 << 1);  /* = 2*k2 */
        if ( ! ispz ) {  /* nonperiodic */
          int lower = km - r_stencil;
          int upper = km + r_stencil;
          if (lower < ka1) lower = ka1;
          if (upper > kb1) upper = kb1;  /* clip edges */
          phi = phi_factored + r_stencil;  /* center of stencil */
          for (k = lower;  k <= upper;  k++) {
            ezd[k] += phi[k-km] * e2h[index_e2];
          }
        }
        else {  /* periodic */
          int kp = km - r_stencil;  /* index at left end of stencil */
          if (kp < ka1) do { kp += nk1; } while (kp < ka1);  /* start inside */
          phi = phi_factored;  /* left end of stencil */
          for (k = 0;  k < diam_stencil;  k++, kp++) {
            if (kp > kb1) kp = ka1;  /* wrap around edge of grid */
            ezd[kp] += phi[k] * e2h[index_e2];
          }
        }
      } /* for k2 */

      for (k = ka1;  k <= kb1;  k++) {
        offset = k * nj1;
        jm = (j2 << 1);  /* = 2*j2 */
        if ( ! ispy ) {  /* nonperiodic */
          int lower = jm - r_stencil;
          int upper = jm + r_stencil;
          if (lower < ja1) lower = ja1;
          if (upper > jb1) upper = jb1;  /* clip edges */
          phi = phi_factored + r_stencil;  /* center of stencil */
          for (j = lower;  j <= upper;  j++) {
            eyzd[offset + j] += phi[j-jm] * ezd[k];
          }
        }
        else {  /* periodic */
          int jp = jm - r_stencil;  /* index at left end of stencil */
          if (jp < ja1) do { jp += nj1; } while (jp < ja1);  /* start inside */
          phi = phi_factored;  /* left end of stencil */
          for (j = 0;  j < diam_stencil;  j++, jp++) {
            if (jp > jb1) jp = ja1;  /* wrap around edge of grid */
            eyzd[offset + jp] += phi[j] * ezd[k];
          }
        }
      } /* for k */

    } /* for j2 */

    for (k = ka1;  k <= kb1;  k++) {
      offset_k = k * nj1;

      for (j = ja1;  j <= jb1;  j++) {
        index_jk = offset_k + j;
        offset = index_jk * ni1;
        im = (i2 << 1);  /* = 2*i2 */
        if ( ! ispx ) {  /* nonperiodic */
          int lower = im - r_stencil;
          int upper = im + r_stencil;
          if (lower < ia1) lower = ia1;
          if (upper > ib1) upper = ib1;  /* clip edges */
          phi = phi_factored + r_stencil;  /* center of stencil */
          for (i = lower;  i <= upper;  i++) {
            eh[offset + i] += phi[i-im] * eyzd[index_jk];
          }
        }
        else {  /* periodic */
          int ip = im - r_stencil;  /* index at left end of stencil */
          if (ip < ia1) do { ip += ni1; } while (ip < ia1);  /* start inside */
          phi = phi_factored;  /* left end of stencil */
          for (i = 0;  i < diam_stencil;  i++, ip++) {
            if (ip > ib1) ip = ia1;  /* wrap around edge of grid */
            eh[offset + ip] += phi[i] * eyzd[index_jk];
          }
        }
      } /* for j */

    } /* for k */

  } /* for i2 */

  return NL_MSM_SUCCESS;
} /* prolongation_factored */


int restriction(NL_Msm *msm, int level) {
  /* grids of charge, finer grid and coarser grid */
  const NL_Msmgrid_double *qhgrid = &(msm->qh[level]);
  const double *qh = qhgrid->data;        /* index the offset data buffer */
  NL_Msmgrid_double *q2hgrid = &(msm->qh[level+1]);
  double *q2h_buffer = q2hgrid->buffer;   /* index the raw buffer */

  /* finer grid index ranges and dimensions */
  const int ia1 = qhgrid->i0;             /* lowest x-index */
  const int ib1 = ia1 + qhgrid->ni - 1;   /* highest x-index */
  const int ja1 = qhgrid->j0;             /* lowest y-index */
  const int jb1 = ja1 + qhgrid->nj - 1;   /* highest y-index */
  const int ka1 = qhgrid->k0;             /* lowest z-index */
  const int kb1 = ka1 + qhgrid->nk - 1;   /* highest z-index */
  const int ni1 = qhgrid->ni;             /* length along x-dim */
  const int nj1 = qhgrid->nj;             /* length along y-dim */
  const int nk1 = qhgrid->nk;             /* length along z-dim */

  /* coarser grid index ranges and dimensions */
  const int ia2 = q2hgrid->i0;            /* lowest x-index */
  const int ib2 = ia2 + q2hgrid->ni - 1;  /* highest x-index */
  const int ja2 = q2hgrid->j0;            /* lowest y-index */
  const int jb2 = ja2 + q2hgrid->nj - 1;  /* highest y-index */
  const int ka2 = q2hgrid->k0;            /* lowest z-index */
  const int kb2 = ka2 + q2hgrid->nk - 1;  /* highest z-index */

  const int ispx = (msm->msmflags & NL_MSM_PERIODIC_VEC1);
  const int ispy = (msm->msmflags & NL_MSM_PERIODIC_VEC2);
  const int ispz = (msm->msmflags & NL_MSM_PERIODIC_VEC3);

  const int nstencil = Nstencil[msm->approx];
  const int *offset = IndexOffset[msm->approx];
  const double *phi = PhiStencil[msm->approx];

  double q2h_sum, cjk;

  int i1, j1, k1, index1, jk1off, k1off;
  int i2, j2, k2;
  int index2;
  int i, j, k;

  for (index2 = 0, k2 = ka2;  k2 <= kb2;  k2++) {
    k1 = k2 * 2;
    for (j2 = ja2;  j2 <= jb2;  j2++) {
      j1 = j2 * 2;
      for (i2 = ia2;  i2 <= ib2;  i2++, index2++) {
        i1 = i2 * 2;

        q2h_sum = 0;
        for (k = 0;  k < nstencil;  k++) {
          k1off = k1 + offset[k];
          if (k1off < ka1) {
            if (ispz) do { k1off += nk1; } while (k1off < ka1);
            else continue;
          }
          else if (k1off > kb1) {
            if (ispz) do { k1off -= nk1; } while (k1off > kb1);
            else break;
          }
          k1off *= nj1;
          for (j = 0;  j < nstencil;  j++) {
            jk1off = j1 + offset[j];
            if (jk1off < ja1) {
              if (ispy) do { jk1off += nj1; } while (jk1off < ja1);
              else continue;
            }
            else if (jk1off > jb1) {
              if (ispy) do { jk1off -= nj1; } while (jk1off > jb1);
              else break;
            }
            jk1off = (k1off + jk1off) * ni1;
            cjk = phi[j] * phi[k];
            for (i = 0;  i < nstencil;  i++) {
              index1 = i1 + offset[i];
              if (index1 < ia1) {
                if (ispx) do { index1 += ni1; } while (index1 < ia1);
                else continue;
              }
              else if (index1 > ib1) {
                if (ispx) do { index1 -= ni1; } while (index1 > ib1);
                else break;
              }
              index1 += jk1off;
              q2h_sum += qh[index1] * phi[i] * cjk;
            }
          }
        } /* end loop over finer grid stencil */

        q2h_buffer[index2] = q2h_sum;  /* store charge to coarser grid */
      }
    }
  } /* end loop over each coarser grid point */

  return NL_MSM_SUCCESS;
}


int prolongation(NL_Msm *msm, int level) {
  /* grids of charge, finer grid and coarser grid */
  NL_Msmgrid_double *ehgrid = &(msm->eh[level]);
  double *eh = ehgrid->data;        /* index the offset data buffer */
  const NL_Msmgrid_double *e2hgrid = &(msm->eh[level+1]);
  const double *e2h_buffer = e2hgrid->buffer;   /* index the raw buffer */

  /* finer grid index ranges and dimensions */
  const int ia1 = ehgrid->i0;             /* lowest x-index */
  const int ib1 = ia1 + ehgrid->ni - 1;   /* highest x-index */
  const int ja1 = ehgrid->j0;             /* lowest y-index */
  const int jb1 = ja1 + ehgrid->nj - 1;   /* highest y-index */
  const int ka1 = ehgrid->k0;             /* lowest z-index */
  const int kb1 = ka1 + ehgrid->nk - 1;   /* highest z-index */
  const int ni1 = ehgrid->ni;             /* length along x-dim */
  const int nj1 = ehgrid->nj;             /* length along y-dim */
  const int nk1 = ehgrid->nk;             /* length along z-dim */

  /* coarser grid index ranges and dimensions */
  const int ia2 = e2hgrid->i0;            /* lowest x-index */
  const int ib2 = ia2 + e2hgrid->ni - 1;  /* highest x-index */
  const int ja2 = e2hgrid->j0;            /* lowest y-index */
  const int jb2 = ja2 + e2hgrid->nj - 1;  /* highest y-index */
  const int ka2 = e2hgrid->k0;            /* lowest z-index */
  const int kb2 = ka2 + e2hgrid->nk - 1;  /* highest z-index */

  const int ispx = (msm->msmflags & NL_MSM_PERIODIC_VEC1);
  const int ispy = (msm->msmflags & NL_MSM_PERIODIC_VEC2);
  const int ispz = (msm->msmflags & NL_MSM_PERIODIC_VEC3);

  const int nstencil = Nstencil[msm->approx];
  const int *offset = IndexOffset[msm->approx];
  const double *phi = PhiStencil[msm->approx];

  double cjk;

  int i1, j1, k1, index1, jk1off, k1off;
  int i2, j2, k2;
  int index2;
  int i, j, k;

  for (index2 = 0, k2 = ka2;  k2 <= kb2;  k2++) {
    k1 = k2 * 2;
    for (j2 = ja2;  j2 <= jb2;  j2++) {
      j1 = j2 * 2;
      for (i2 = ia2;  i2 <= ib2;  i2++, index2++) {
        i1 = i2 * 2;

        for (k = 0;  k < nstencil;  k++) {
          k1off = k1 + offset[k];
          if (k1off < ka1) {
            if (ispz) do { k1off += nk1; } while (k1off < ka1);
            else continue;
          }
          else if (k1off > kb1) {
            if (ispz) do { k1off -= nk1; } while (k1off > kb1);
            else break;
          }
          k1off *= nj1;
          for (j = 0;  j < nstencil;  j++) {
            jk1off = j1 + offset[j];
            if (jk1off < ja1) {
              if (ispy) do { jk1off += nj1; } while (jk1off < ja1);
              else continue;
            }
            else if (jk1off > jb1) {
              if (ispy) do { jk1off -= nj1; } while (jk1off > jb1);
              else break;
            }
            jk1off = (k1off + jk1off) * ni1;
            cjk = phi[j] * phi[k];
            for (i = 0;  i < nstencil;  i++) {
              index1 = i1 + offset[i];
              if (index1 < ia1) {
                if (ispx) do { index1 += ni1; } while (index1 < ia1);
                else continue;
              }
              else if (index1 > ib1) {
                if (ispx) do { index1 -= ni1; } while (index1 > ib1);
                else break;
              }
              index1 += jk1off;
              eh[index1] += e2h_buffer[index2] * phi[i] * cjk;
            }
          }
        } /* end loop over finer grid stencil */

      }
    }
  } /* end loop over each coarser grid point */

  return NL_MSM_SUCCESS;
}


int gridcutoff(NL_Msm *msm, int level)
{
  double eh_sum;

  /* grids of charge and potential */
  const NL_Msmgrid_double *qhgrid = &(msm->qh[level]);
  const double *qh = qhgrid->data;
  NL_Msmgrid_double *ehgrid = &(msm->eh[level]);
  double *eh = ehgrid->data;
  const int ia = qhgrid->i0;            /* lowest x-index */
  const int ib = ia + qhgrid->ni - 1;   /* highest x-index */
  const int ja = qhgrid->j0;            /* lowest y-index */
  const int jb = ja + qhgrid->nj - 1;   /* highest y-index */
  const int ka = qhgrid->k0;            /* lowest z-index */
  const int kb = ka + qhgrid->nk - 1;   /* highest z-index */
  const int ni = qhgrid->ni;            /* length along x-dim */
  const int nj = qhgrid->nj;            /* length along y-dim */
  const int nk = qhgrid->nk;            /* length along z-dim */

  /* grid of weights for pairwise grid point interactions within cutoff */
  const NL_Msmgrid_double *gcgrid = &(msm->gc[level]);
  const double *gc = gcgrid->data;
  const int gia = gcgrid->i0;            /* lowest x-index */
  const int gib = gia + gcgrid->ni - 1;  /* highest x-index */
  const int gja = gcgrid->j0;            /* lowest y-index */
  const int gjb = gja + gcgrid->nj - 1;  /* highest y-index */
  const int gka = gcgrid->k0;            /* lowest z-index */
  const int gkb = gka + gcgrid->nk - 1;  /* highest z-index */
  const int gni = gcgrid->ni;            /* length along x-dim */
  const int gnj = gcgrid->nj;            /* length along y-dim */

  const int ispx = (msm->msmflags & NL_MSM_PERIODIC_VEC1);
  const int ispy = (msm->msmflags & NL_MSM_PERIODIC_VEC2);
  const int ispz = (msm->msmflags & NL_MSM_PERIODIC_VEC3);

  const int ispnone = !(ispx || ispy || ispz);

  int i, j, k;
  int gia_clip, gib_clip;
  int gja_clip, gjb_clip;
  int gka_clip, gkb_clip;
  int koff;
  long jkoff, index;
  int id, jd, kd;
  int knoff;
  long jknoff, nindex;
  int kgoff, jkgoff, ngindex;

  if ( ispnone ) {  /* non-periodic boundaries */

    /* loop over all grid points */
    for (k = ka;  k <= kb;  k++) {

      /* clip gc ranges to keep offset for k index within grid */
      gka_clip = (k + gka < ka ? ka - k : gka);
      gkb_clip = (k + gkb > kb ? kb - k : gkb);

      koff = k * nj;  /* find eh flat index */

      for (j = ja;  j <= jb;  j++) {

        /* clip gc ranges to keep offset for j index within grid */
        gja_clip = (j + gja < ja ? ja - j : gja);
        gjb_clip = (j + gjb > jb ? jb - j : gjb);

        jkoff = (koff + j) * ni;  /* find eh flat index */

        for (i = ia;  i <= ib;  i++) {

          /* clip gc ranges to keep offset for i index within grid */
          gia_clip = (i + gia < ia ? ia - i : gia);
          gib_clip = (i + gib > ib ? ib - i : gib);

          index = jkoff + i;  /* eh flat index */

          /* sum over "sphere" of weighted charge */
          eh_sum = 0;
          for (kd = gka_clip;  kd <= gkb_clip;  kd++) {
            knoff = (k + kd) * nj;  /* find qh flat index */
            kgoff = kd * gnj;       /* find gc flat index */

            for (jd = gja_clip;  jd <= gjb_clip;  jd++) {
              jknoff = (knoff + (j + jd)) * ni;  /* find qh flat index */
              jkgoff = (kgoff + jd) * gni;       /* find gc flat index */

              for (id = gia_clip;  id <= gib_clip;  id++) {
                nindex = jknoff + (i + id);  /* qh flat index */
                ngindex = jkgoff + id;       /* gc flat index */

                GRID_INDEX_CHECK(qhgrid, i+id, j+jd, k+kd);
                ASSERT(GRID_INDEX(qhgrid, i+id, j+jd, k+kd) == nindex);

                GRID_INDEX_CHECK(gcgrid, id, jd, kd);
                ASSERT(GRID_INDEX(gcgrid, id, jd, kd) == ngindex);

                eh_sum += qh[nindex] * gc[ngindex];  /* sum weighted charge */
              }
            }
          } /* end loop over "sphere" of charge */

          GRID_INDEX_CHECK(ehgrid, i, j, k);
          ASSERT(GRID_INDEX(ehgrid, i, j, k) == index);
          eh[index] = eh_sum;  /* store potential */
        }
      }
    } /* end loop over all grid points */

  } /* if nonperiodic boundaries */
  else {
    /* some boundary is periodic */
    int ilo, jlo, klo;
    int ip, jp, kp;

    /* loop over all grid points */
    for (k = ka;  k <= kb;  k++) {
      klo = k + gka;
      if ( ! ispz ) {  /* nonperiodic z */
        /* clip gc ranges to keep offset for k index within grid */
        gka_clip = (k + gka < ka ? ka - k : gka);
        gkb_clip = (k + gkb > kb ? kb - k : gkb);
        if (klo < ka) klo = ka;  /* keep lowest qh index within grid */
      }
      else {  /* periodic z */
        gka_clip = gka;
        gkb_clip = gkb;
        if (klo < ka) do { klo += nk; } while (klo < ka);
      }
      ASSERT(klo <= kb);

      koff = k * nj;  /* find eh flat index */

      for (j = ja;  j <= jb;  j++) {
        jlo = j + gja;
        if ( ! ispy ) {  /* nonperiodic y */
          /* clip gc ranges to keep offset for j index within grid */
          gja_clip = (j + gja < ja ? ja - j : gja);
          gjb_clip = (j + gjb > jb ? jb - j : gjb);
          if (jlo < ja) jlo = ja;  /* keep lowest qh index within grid */
        }
        else {  /* periodic y */
          gja_clip = gja;
          gjb_clip = gjb;
          if (jlo < ja) do { jlo += nj; } while (jlo < ja);
        }
        ASSERT(jlo <= jb);

        jkoff = (koff + j) * ni;  /* find eh flat index */

        for (i = ia;  i <= ib;  i++) {
          ilo = i + gia;
          if ( ! ispx ) {  /* nonperiodic x */
            /* clip gc ranges to keep offset for i index within grid */
            gia_clip = (i + gia < ia ? ia - i : gia);
            gib_clip = (i + gib > ib ? ib - i : gib);
            if (ilo < ia) ilo = ia;  /* keep lowest qh index within grid */
          }
          else {  /* periodic x */
            gia_clip = gia;
            gib_clip = gib;
            if (ilo < ia) do { ilo += ni; } while (ilo < ia);
          }
          ASSERT(ilo <= ib);

          index = jkoff + i;  /* eh flat index */

          /* sum over "sphere" of weighted charge */
          eh_sum = 0;
          for (kd = gka_clip, kp = klo;  kd <= gkb_clip;  kd++, kp++) {
            /* clipping makes conditional always fail for nonperiodic */
            if (kp > kb) kp = ka;  /* wrap z direction */
            knoff = kp * nj;       /* find qh flat index */
            kgoff = kd * gnj;      /* find gc flat index */

            for (jd = gja_clip, jp = jlo;  jd <= gjb_clip;  jd++, jp++) {
              /* clipping makes conditional always fail for nonperiodic */
              if (jp > jb) jp = ja;              /* wrap y direction */
              jknoff = (knoff + jp) * ni;  /* find qh flat index */
              jkgoff = (kgoff + jd) * gni;       /* find gc flat index */

              for (id = gia_clip, ip = ilo;  id <= gib_clip;  id++, ip++) {
                /* clipping makes conditional always fail for nonperiodic */
                if (ip > ib) ip = ia;   /* wrap x direction */
                nindex = jknoff +  ip;  /* qh flat index */
                ngindex = jkgoff + id;  /* gc flat index */

                GRID_INDEX_CHECK(qhgrid, ip, jp, kp);
                ASSERT(GRID_INDEX(qhgrid, ip, jp, kp) == nindex);

                GRID_INDEX_CHECK(gcgrid, id, jd, kd);
                ASSERT(GRID_INDEX(gcgrid, id, jd, kd) == ngindex);

                eh_sum += qh[nindex] * gc[ngindex];  /* sum weighted charge */

              }
            }
          } /* end loop over "sphere" of charge */

          GRID_INDEX_CHECK(ehgrid, i, j, k);
          ASSERT(GRID_INDEX(ehgrid, i, j, k) == index);
          eh[index] = eh_sum;  /* store potential */
        }
      }
    } /* end loop over all grid points */

  } /* else some boundary is periodic */

  return NL_MSM_SUCCESS;
}
