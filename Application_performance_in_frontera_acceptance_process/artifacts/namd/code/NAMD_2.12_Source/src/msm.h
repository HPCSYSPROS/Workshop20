/* msm.h */

#ifndef NL_MSM_H
#define NL_MSM_H

#ifdef __cplusplus
extern "C" {
#endif


  /** Private MSM data structure */
  struct NL_Msm_t;
  typedef struct NL_Msm_t NL_Msm;

  /** Allocate MSM solver. */
  NL_Msm *NL_msm_create(void);

  /** Free MSM solver. */
  void NL_msm_destroy(NL_Msm *);

  /** Setup the MSM solver for the molecular system.  These parameters come
   * directly from the simulation.  Here the cutoff provides some control
   * over the accuracy.  The region of space containing the atoms is defined 
   * as a parallelepiped, and the "isperiodic" flag determines whether or
   * not each individual cell dimension has periodic boundary conditions.  
   * The defined cell is expected to contain the atoms regardless of the 
   * choice of boundary conditions.  The cell volume must contain the
   * cutoff-sphere.  
   *
   * XXX For now, the MSM solver permits only cell vectors that align with 
   * the x-, y-, and z-axis, respectively. */
  int NL_msm_setup(
      NL_Msm *msm,         /**< the MSM solver object */
      double cutoff,       /**< cutoff distance for short-range part */
      double cellvec1[3],  /**< cell basis vector 1 */
      double cellvec2[3],  /**< cell basis vector 2 */
      double cellvec3[3],  /**< cell basis vector 3 */
      double cellcenter[3],/**< center of cell */
      int msmflags         /**< flags for periodicity and calculation */
      );

  /** Use in setup for msmflags, bitwise OR the options. */
  enum {
    NL_MSM_PERIODIC_NONE     = 0x000,  /**< aperiodic cell */
    NL_MSM_PERIODIC_VEC1     = 0x001,  /**< periodic along basis vector 1 */
    NL_MSM_PERIODIC_VEC2     = 0x002,  /**< periodic along basis vector 2 */
    NL_MSM_PERIODIC_VEC3     = 0x004,  /**< periodic along basis vector 3 */
    NL_MSM_PERIODIC_ALL      = 0x007,  /**< periodic along all basis vectors */

    NL_MSM_COMPUTE_SHORT_RANGE
                             = 0x008,  /**< compute the short-range part */
    NL_MSM_COMPUTE_LONG_RANGE
                             = 0x010,  /**< compute the long-range part */
    NL_MSM_COMPUTE_ALL       = 0x018,  /**< compute both parts */

    NL_MSM_COMPUTE_SPREC     = 0x020,  /**< use single precision (faster) */
    NL_MSM_COMPUTE_1AWAY     = 0x040,  /**< 1-away neighborhood (slower) */
    NL_MSM_COMPUTE_NONFACTORED
                             = 0x080,  /**< use non-factored restriction and
                                         prolongation procedures (slower) */

    NL_MSM_REPORT_TIMINGS    = 0x100,  /**< report timings */

    NL_MSM_COMPUTE_CUDA_FALL_BACK
                             = 0x200,  /**< fall back on CPU
                                         if configuration unsupported by CUDA */
    NL_MSM_COMPUTE_CUDA_GRID_CUTOFF
                             = 0x400,  /**< use CUDA for grid cutoff */
    NL_MSM_COMPUTE_CUDA_SHORT_RANGE
                             = 0x800,  /**< use CUDA for short range */
    NL_MSM_COMPUTE_CUDA      = 0xE00,  /**< all CUDA flags set on */

    NL_MSM_ALL_FLAGS         = 0xFFF   /**< all flags set on */
  };

  /** Return codes. */
  enum {
    NL_MSM_SUCCESS = 0,    /**< success */
    NL_MSM_ERROR_MALLOC,   /**< can't allocate memory */
    NL_MSM_ERROR_PARAM,    /**< incorrect parameters */
    NL_MSM_ERROR_RANGE,    /**< atom is outside of cell */
    NL_MSM_ERROR_SUPPORT,  /**< unsupported feature */
    NL_MSM_ERROR_CUDA,     /**< CUDA failure */
    NL_MSM_ERROR_END       /**< (for internal use) */
  };

  /** Advanced configuration of MSM.  The ratio of cutoff / gridspacing 
   * is generally kept between about 2.5 and 6.  For atomic lengths in
   * Angstroms, good results are demonstrated for a cutoff distance
   * between 8 and 12 A and the default grid spacing of 2.5 A.
   *
   * The optimal pairings of approx and split have been demonstrated to be: 
   * - CUBIC with TAYLOR2 (the default),
   * - QUINTIC with TAYLOR3,
   * - SEPTIC with TAYLOR4,
   * - NONIC with TAYLOR5.
   *
   * Use QUINTIC2, SEPTIC3, NONIC4 paired with the same splitting functions
   * above if greater continuity is desired.
   *
   * Set "nlevels" to 0 (the default) to fully adapt the number of
   * grid levels to the system cell size. */
  int NL_msm_configure(
      NL_Msm *msm,         /**< the MSM solver object */
      double gridspacing,  /**< grid spacing for first grid level */
      int approx,          /**< which approximation method */
      int split,           /**< which splitting */
      int nlevels          /**< number of grid levels to use */
      );

  /** Compute the electrostatic forces and potential energy for the 
   * array of atoms.  The felec array and the location pointed to by uelec
   * is expected to be initialized before the call.  As stated, the atoms
   * must be within the defined cell. */
  int NL_msm_compute_force(
      NL_Msm *msm,         /**< the MSM solver object */
      double *felec,       /**< electrostatic forces x/y/z for each atom */
      double *uelec,       /**< electrostatic potential energy */
      const double *atom,  /**< positions and charge x/y/z/q for each atom */
      int natoms           /**< number of atoms */
      );

  /** Same as NL_msm_compute_force() except for single precision calculation.
   * Call this only if NL_MSM_COMPUTE_SPREC flag was selected. */
  int NL_msm_compute_force_sprec(
      NL_Msm *msm,         /**< the MSM solver object */
      float *felec,        /**< electrostatic forces x/y/z for each atom */
      float *uelec,        /**< electrostatic potential energy */
      const float *atom,   /**< positions and charge x/y/z/q for each atom */
      int natoms           /**< number of atoms */
      );

#if 0
  /** Compute the electrostatic potential map for the array of atoms. 
   * Each maplen component is (0,1], where all 1's mean the map extends 
   * over the entire cell.  The maplen with mapcenter arguments allow the 
   * map to be defined as any fractional sub-cell whose sides have same
   * orientation as the system cell.  The map must fit within the cell.  */
  int NL_msm_compute_epotmap(
      NL_Msm *msm,         /**< the MSM solver object */
      double *epotmap,     /**< electrostatic potential map */
      int mx,              /**< size of map along cell basis vector 1 */
      int my,              /**< size of map along cell basis vector 2 */
      int mz,              /**< size of map along cell basis vector 3 */
      double maplen[3],    /**< fraction of map along cell basis vectors */
      double mapcenter[3], /**< center of map */
      const double *atom,  /**< positions and charge x/y/z/q for each atom */
      int natoms           /**< number of atoms */
      );
#endif

  /** MSM approximation/interpolation methods.  (Default is CUBIC.) */
  enum {
    NL_MSM_APPROX_CUBIC = 0, /**< C1 degree 3 poly (numerical Hermite) */
    NL_MSM_APPROX_QUINTIC,   /**< C1 degree 5 poly (linear blend of degree 4) */
    NL_MSM_APPROX_QUINTIC2,  /**< C2 degree 5 poly */
    NL_MSM_APPROX_SEPTIC,    /**< C1 degree 7 poly (linear blend of degree 6) */
    NL_MSM_APPROX_SEPTIC3,   /**< C3 degree 7 poly */
    NL_MSM_APPROX_NONIC,     /**< C1 degree 9 poly (linear blend of degree 8) */
    NL_MSM_APPROX_NONIC4,    /**< C4 degree 9 poly */
    NL_MSM_APPROX_BSPLINE,   /**< C2 degree 3 B-spline */
    NL_MSM_APPROX_END        /**< (for internal use) */
  };

  /** MSM splitting functions to smooth the 1/r potential.
   * (Default is TAYLOR2.) 
   *
   * The Taylor splittings are the family of C^k Taylor polynomials of
   * s^(-1/2) about s=1.  Discussed in section 5.1 of the thesis.
   * TAYLOR1 is listed last so that reasonable default parameters
   * (CUBIC approximation with TAYLOR2 splitting) are automatically
   * selected from zeroing the values.  
   *
   * The sigma(k,d) splittings generalize the Taylor splittings above
   * for C^k continuity at r=1 and degree d polynomial in r.
   * (The Taylor splittings are sigma(d/2,d).)  The ones here are
   * uniquely determined.  Discussed in section 5.2 of the thesis.  
   *
   * The special sigma(2,6) splitting is from section 5.3 of the thesis.  
   * This is a proof-of-concept where the polynomial in even powers of r
   * has an extra degree of freedom chosen to minimize the error bound.
   *
   * Elimiate the self-force artifact through the use of a switching
   * function between 1/r and a polynomial reproduced exactly by CUBIC
   * interpolation.  Discussed in section 5.4 of thesis.
   */
  enum {
    NL_MSM_SPLIT_TAYLOR2 = 0,  /**< C2 Taylor */
    NL_MSM_SPLIT_TAYLOR3,      /**< C3 Taylor */
    NL_MSM_SPLIT_TAYLOR4,      /**< C4 Taylor */
    NL_MSM_SPLIT_TAYLOR5,      /**< C5 Taylor */
    NL_MSM_SPLIT_TAYLOR6,      /**< C6 Taylor */
    NL_MSM_SPLIT_TAYLOR7,      /**< C7 Taylor */
    NL_MSM_SPLIT_TAYLOR8,      /**< C8 Taylor */
    NL_MSM_SPLIT_TAYLOR1,      /**< C1 Taylor */

    NL_MSM_SPLIT_SIGMA2_3,     /**< C2, degree 3 (the "perfect" smoothing) */
    NL_MSM_SPLIT_SIGMA3_5,     /**< C3, degree 5 */
    NL_MSM_SPLIT_SIGMA4_6,     /**< C4, degree 6 */
    NL_MSM_SPLIT_SIGMA4_7,     /**< C4, degree 7 */
    NL_MSM_SPLIT_SIGMA5_8,     /**< C5, degree 8 */
    NL_MSM_SPLIT_SIGMA5_9,     /**< C5, degree 9 */
    NL_MSM_SPLIT_SIGMA6_9,     /**< C6, degree 9 */
    NL_MSM_SPLIT_SIGMA6_10,    /**< C6, degree 10 */
    NL_MSM_SPLIT_SIGMA6_11,    /**< C6, degree 11 */
    NL_MSM_SPLIT_SIGMA7_11,    /**< C7, degree 11 */
    NL_MSM_SPLIT_SIGMA7_12,    /**< C7, degree 12 */
    NL_MSM_SPLIT_SIGMA7_13,    /**< C7, degree 13 */
    NL_MSM_SPLIT_SIGMA8_12,    /**< C8, degree 12 */
    NL_MSM_SPLIT_SIGMA8_13,    /**< C8, degree 13 */
    NL_MSM_SPLIT_SIGMA8_14,    /**< C8, degree 14 */
    NL_MSM_SPLIT_SIGMA8_15,    /**< C8, degree 15 */

    NL_MSM_SPLIT_SIGMA2_6,     /**< C2, degree 6, even powers of r,
                                 chosen to minimize error bound */

    NL_MSM_SPLIT_SWITCH1_2,    /**< C2, switching at r/a=1/2 */
    NL_MSM_SPLIT_SWITCH3_4,    /**< C2, switching at r/a=3/4 */
    NL_MSM_SPLIT_SWITCH7_8,    /**< C2, switching at r/a=7/8 */

    NL_MSM_SPLIT_END           /**< (for internal use) */
  };

  /** Helper function to determine APPROX enum constant from string name. */
  int NL_msm_approx(const char *name);

  /** Helper function to determine SPLIT enum constant from string name. */
  int NL_msm_split(const char *name);

  /** Helper function returning string name for APPROX enum constant. */
  const char *NL_msm_approx_name(int approx);

  /** Helper function returning string name for SPLIT enum constant. */
  const char *NL_msm_split_name(int split);


#ifdef __cplusplus
}
#endif

#endif /* MSM_H */
