#ifndef _CCTK_LOOP_H_
#define _CCTK_LOOP_H_

/* WARNING: This file is auto-generated. Do not edit. */
/* Edit cctk_Loop.h.pl instead, and then re-generate this file via */
/*    perl cctk_Loop.h.pl */
/* Documentation can also be found in "cctk_Loop.h.pl". */

#ifdef CCODE
#  include <cctk_Config.h>
#  include <cctk_WarnLevel.h>
#  include <cGH.h>
#  include <assert.h>

#  ifndef CCTK_DISABLE_OMP_COLLAPSE
#    define CCTK_PRAGMA_OMP_FOR_COLLAPSE_1 _Pragma("omp for collapse(1)")
#    define CCTK_PRAGMA_OMP_FOR_COLLAPSE_2 _Pragma("omp for collapse(2)")
#    define CCTK_PRAGMA_OMP_FOR_COLLAPSE_3 _Pragma("omp for collapse(3)")
#    define CCTK_PRAGMA_OMP_FOR_COLLAPSE_4 _Pragma("omp for collapse(4)")
#  else
#    define CCTK_PRAGMA_OMP_FOR_COLLAPSE_1 _Pragma("omp for")
#    define CCTK_PRAGMA_OMP_FOR_COLLAPSE_2 _Pragma("omp for")
#    define CCTK_PRAGMA_OMP_FOR_COLLAPSE_3 _Pragma("omp for")
#    define CCTK_PRAGMA_OMP_FOR_COLLAPSE_4 _Pragma("omp for")
#  endif
#endif /* #ifdef CCODE */



/* 1D */

#ifdef CCODE

/* LOOP */

#define CCTK_LOOP1_NORMAL(name, \
                          i, \
                          ni, \
                          cctki0_idir_, \
                          cctki0_imin_, \
                          cctki0_imax_, \
                          cctki0_iash_) \
  CCTK_LOOP1STR_NORMAL(name, \
                       i, \
                       ni, \
                       (cctki0_idir_), \
                       (cctki0_imin_), \
                       (cctki0_imax_), \
                       (cctki0_iash_), \
                       cctki0_dummy_imin,cctki0_dummy_imax, 1) \

#define CCTK_ENDLOOP1_NORMAL(name) \
  CCTK_ENDLOOP1STR_NORMAL(name) \

#define CCTK_LOOP1STR_NORMAL(name, \
                             i, \
                             ni, \
                             cctki0_idir_, \
                             cctki0_imin_, \
                             cctki0_imax_, \
                             cctki0_iash_, \
                             imin,imax, cctki0_istr_) \
  do { \
    typedef int cctki0_loop1_normal_##name; \
    const int cctki0_idir = (cctki0_idir_); \
    const int cctki0_imin = (cctki0_imin_); \
    const int cctki0_imax = (cctki0_imax_); \
    const int cctki0_iash CCTK_ATTRIBUTE_UNUSED = (cctki0_iash_); \
    const int cctki0_istr = (cctki0_istr_); \
    assert(cctki0_istr>0 && (cctki0_istr & (cctki0_istr-1)) == 0); \
    const int imin CCTK_ATTRIBUTE_UNUSED = cctki0_imin; \
    const int imax CCTK_ATTRIBUTE_UNUSED = cctki0_imax; \
     \
     \
    const int cctki0_ioff = (cctki0_imin) & (cctki0_istr-1); \
    _Pragma("omp for") \
    for (int i=cctki0_imin-cctki0_ioff; i<cctki0_imax; i+=cctki0_istr) { \
      const int ni CCTK_ATTRIBUTE_UNUSED = cctki0_idir<0 ? i+1 : cctki0_idir==0 ? 0 : cctki0_imax-i; \
      { \

#define CCTK_ENDLOOP1STR_NORMAL(name) \
      } \
    } \
    typedef cctki0_loop1_normal_##name cctki0_ensure_proper_nesting CCTK_ATTRIBUTE_UNUSED; \
  } while (0) \



#define CCTK_LOOP1(name, \
                   i, \
                   cctki1_imin_, \
                   cctki1_imax_, \
                   cctki1_iash_) \
  CCTK_LOOP1STR(name, \
                i, \
                (cctki1_imin_), \
                (cctki1_imax_), \
                (cctki1_iash_), \
                cctki1_dummy_imin,cctki1_dummy_imax, 1) \

#define CCTK_ENDLOOP1(name) \
  CCTK_ENDLOOP1STR(name) \

#define CCTK_LOOP1STR(name, \
                      i, \
                      cctki1_imin_, \
                      cctki1_imax_, \
                      cctki1_iash_, \
                      imin,imax, cctki1_istr_) \
  CCTK_LOOP1STR_NORMAL(name, \
                       i, \
                       cctki1_ni, \
                       0, \
                       (cctki1_imin_), \
                       (cctki1_imax_), \
                       (cctki1_iash_), \
                       imin,imax, (cctki1_istr_)) \

#define CCTK_ENDLOOP1STR(name) \
  CCTK_ENDLOOP1STR_NORMAL(name) \



/* LOOP_INTERIOR */

#define CCTK_LOOP1_INTERIOR(name, cctki2_cctkGH_, \
                            i, \
                            cctki2_iblo_, \
                            cctki2_ibhi_) \
  CCTK_LOOP1STR_INTERIOR(name, (cctki2_cctkGH_), \
                         i, \
                         (cctki2_iblo_), \
                         (cctki2_ibhi_), \
                         cctki2_dummy_imin,cctki2_dummy_imax, 1) \

#define CCTK_ENDLOOP1_INTERIOR(name) \
  CCTK_ENDLOOP1STR_INTERIOR(name) \

#define CCTK_LOOP1STR_INTERIOR(name, cctki2_cctkGH_, \
                               i, \
                               cctki2_iblo_, \
                               cctki2_ibhi_, \
                               imin,imax, cctki2_istr_) \
  do { \
    typedef int cctki2_loop1_interior_##name; \
    cGH const *restrict const cctki2_cctkGH = (cctki2_cctkGH_); \
    if (cctki2_cctkGH->cctk_dim != 1) { \
      _Pragma("omp critical") \
      CCTK_WARN(CCTK_WARN_ABORT, \
                "The macro CCTK_LOOP1_INTERIOR can only be used in 1 dimensions"); \
    } \
    CCTK_LOOP1STR(name##_interior, \
                  i, \
                  (cctki2_iblo_), \
                  cctki2_cctkGH->cctk_lsh[0]-(cctki2_ibhi_), \
                  cctki2_cctkGH->cctk_ash[0], \
                  imin,imax, (cctki2_istr_)) { \

#define CCTK_ENDLOOP1STR_INTERIOR(name) \
    } CCTK_ENDLOOP1STR(name##_interior); \
    typedef cctki2_loop1_interior_##name cctki2_ensure_proper_nesting CCTK_ATTRIBUTE_UNUSED; \
  } while(0) \



/* LOOP_BOUNDARIES */

#define CCTK_LOOP1_BOUNDARIES(name, cctki2_cctkGH_, \
                              i, \
                              ni, \
                              cctki2_iblo_, \
                              cctki2_ibhi_, \
                              cctki2_ibboxlo_, \
                              cctki2_ibboxhi_) \
  CCTK_LOOP1STR_BOUNDARIES(name, (cctki2_cctkGH_), \
                           i, \
                           ni, \
                           (cctki2_iblo_), \
                           (cctki2_ibhi_), \
                           (cctki2_ibboxlo_), \
                           (cctki2_ibboxhi_), \
                           cctki2_dummy_imin,cctki2_dummy_imax, 1) \

#define CCTK_ENDLOOP1_BOUNDARIES(name) \
  CCTK_ENDLOOP1STR_BOUNDARIES(name) \

#define CCTK_LOOP1STR_BOUNDARIES(name, cctki2_cctkGH_, \
                                 i, \
                                 ni, \
                                 cctki2_iblo_, \
                                 cctki2_ibhi_, \
                                 cctki2_ibboxlo_, \
                                 cctki2_ibboxhi_, \
                                 imin,imax, cctki2_istr_) \
  do { \
    typedef int cctki2_loop1_boundaries_##name; \
    cGH const *restrict const cctki2_cctkGH = (cctki2_cctkGH_); \
    if (cctki2_cctkGH->cctk_dim != 1) { \
      _Pragma("omp critical") \
      CCTK_WARN(CCTK_WARN_ABORT, \
                "The macro CCTK_LOOP1_BOUNDARIES can only be used in 1 dimensions"); \
    } \
    const int cctki2_blo[] = { (cctki2_iblo_) }; \
    const int cctki2_bhi[] = { (cctki2_ibhi_) }; \
    const int cctki2_bbox[] = { (cctki2_ibboxlo_), (cctki2_ibboxhi_) }; \
    const int cctki2_lsh[] = { cctki2_cctkGH->cctk_lsh[0] }; \
    const int cctki2_istr CCTK_ATTRIBUTE_UNUSED = (cctki2_istr_); \
    for (int cctki2_idir=-1; cctki2_idir<=+1; ++cctki2_idir) { \
      const int cctki2_any_bbox = \
        (cctki2_idir<0 ? cctki2_bbox[0] : 0) || (cctki2_idir>0 ? cctki2_bbox[1] : 0); \
      if (cctki2_any_bbox) { \
        const int cctki2_bmin[] = { \
          cctki2_idir<0 ? 0 : cctki2_idir==0 ? cctki2_blo[0] : cctki2_lsh[0] - cctki2_bhi[0], \
        }; \
        const int cctki2_bmax[] = { \
          cctki2_idir<0 ? cctki2_blo[0] : cctki2_idir==0 ? cctki2_lsh[0] - cctki2_bhi[0] : cctki2_lsh[0], \
        }; \
        CCTK_LOOP1STR_NORMAL(name##_boundaries, \
                             i, \
                             ni, \
                             cctki2_idir, \
                             cctki2_bmin[0], \
                             cctki2_bmax[0], \
                             cctki2_cctkGH->cctk_ash[0], \
                             imin,imax, cctki2_istr) { \

#define CCTK_ENDLOOP1STR_BOUNDARIES(name) \
        } CCTK_ENDLOOP1STR_NORMAL(name##_boundaries); \
      } /* if bbox */ \
    } /* for dir */ \
    typedef cctki2_loop1_boundaries_##name cctki2_ensure_proper_nesting CCTK_ATTRIBUTE_UNUSED; \
  } while (0) \



/* LOOP_INTBOUNDARIES */

#define CCTK_LOOP1_INTBOUNDARIES(name, cctki2_cctkGH_, \
                                 i, \
                                 ni, \
                                 cctki2_iblo_, \
                                 cctki2_ibhi_, \
                                 cctki2_ibboxlo_, \
                                 cctki2_ibboxhi_) \
  CCTK_LOOP1STR_INTBOUNDARIES(name, (cctki2_cctkGH_), \
                              i, \
                              ni, \
                              (cctki2_iblo_), \
                              (cctki2_ibhi_), \
                              (cctki2_ibboxlo_), \
                              (cctki2_ibboxhi_), \
                              cctki2_dummy_imin,cctki2_dummy_imax, 1) \

#define CCTK_ENDLOOP1_INTBOUNDARIES(name) \
  CCTK_ENDLOOP1STR_INTBOUNDARIES(name) \

#define CCTK_LOOP1STR_INTBOUNDARIES(name, cctki2_cctkGH_, \
                                    i, \
                                    ni, \
                                    cctki2_iblo_, \
                                    cctki2_ibhi_, \
                                    cctki2_ibboxlo_, \
                                    cctki2_ibboxhi_, \
                                    imin,imax, cctki2_istr_) \
  do { \
    typedef int cctki2_loop1_intboundaries_##name; \
    cGH const *restrict const cctki2_cctkGH = (cctki2_cctkGH_); \
    if (cctki2_cctkGH->cctk_dim != 1) { \
      _Pragma("omp critical") \
      CCTK_WARN(CCTK_WARN_ABORT, \
                "The macro CCTK_LOOP1_INTBOUNDARIES can only be used in 1 dimensions"); \
    } \
    const int cctki2_blo[] = { (cctki2_iblo_) }; \
    const int cctki2_bhi[] = { (cctki2_ibhi_) }; \
    const int cctki2_bbox[] = { (cctki2_ibboxlo_), (cctki2_ibboxhi_) }; \
    const int cctki2_lsh[] = { cctki2_cctkGH->cctk_lsh[0] }; \
    const int cctki2_istr CCTK_ATTRIBUTE_UNUSED = (cctki2_istr_); \
    for (int cctki2_idir=-1; cctki2_idir<=+1; ++cctki2_idir) { \
      const int cctki2_any_bbox = \
        (cctki2_idir<0 ? cctki2_bbox[0] : 0) || (cctki2_idir>0 ? cctki2_bbox[1] : 0); \
      const int cctki2_all_bbox = \
        (cctki2_idir<0 ? cctki2_bbox[0] : 1) && (cctki2_idir>0 ? cctki2_bbox[1] : 1); \
      if (cctki2_all_bbox && cctki2_any_bbox) { \
        const int cctki2_bmin[] = { \
          cctki2_idir<0 ? 0 : cctki2_idir==0 ? cctki2_blo[0] : cctki2_lsh[0] - cctki2_bhi[0], \
        }; \
        const int cctki2_bmax[] = { \
          cctki2_idir<0 ? cctki2_blo[0] : cctki2_idir==0 ? cctki2_lsh[0] - cctki2_bhi[0] : cctki2_lsh[0], \
        }; \
        CCTK_LOOP1STR_NORMAL(name##_intboundaries, \
                             i, \
                             ni, \
                             cctki2_idir, \
                             cctki2_bmin[0], \
                             cctki2_bmax[0], \
                             cctki2_cctkGH->cctk_ash[0], \
                             imin,imax, cctki2_istr) { \

#define CCTK_ENDLOOP1STR_INTBOUNDARIES(name) \
        } CCTK_ENDLOOP1STR_NORMAL(name##_intboundaries); \
      } /* if bbox */ \
    } /* for dir */ \
    typedef cctki2_loop1_intboundaries_##name cctki2_ensure_proper_nesting CCTK_ATTRIBUTE_UNUSED; \
  } while (0) \



/* LOOP_ALL */

#define CCTK_LOOP1_ALL(name, cctki3_cctkGH_, \
                       i) \
  CCTK_LOOP1STR_ALL(name, (cctki3_cctkGH_), \
                    i, \
                    cctki3_dummy_imin,cctki3_dummy_imax, 1) \

#define CCTK_ENDLOOP1_ALL(name) \
  CCTK_ENDLOOP1STR_ALL(name) \

#define CCTK_LOOP1STR_ALL(name, cctki3_cctkGH_, \
                          i, \
                          imin,imax, cctki3_istr_) \
  do { \
    typedef int cctki3_loop1_all_##name; \
    cGH const *restrict const cctki3_cctkGH = (cctki3_cctkGH_); \
    if (cctki3_cctkGH->cctk_dim != 1) { \
      _Pragma("omp critical") \
      CCTK_WARN(CCTK_WARN_ABORT, \
                "The macro CCTK_LOOP1_ALL can only be used in 1 dimensions"); \
    } \
    CCTK_LOOP1STR(name##_all, \
                  i, \
                  0, \
                  cctki3_cctkGH->cctk_lsh[0], \
                  cctki3_cctkGH->cctk_ash[0], \
                  imin,imax, (cctki3_istr_)) { \

#define CCTK_ENDLOOP1STR_ALL(name) \
    } CCTK_ENDLOOP1STR(name##_all); \
    typedef cctki3_loop1_all_##name cctki3_ensure_proper_nesting CCTK_ATTRIBUTE_UNUSED; \
  } while (0) \



/* LOOP_INT */

#define CCTK_LOOP1_INT(name, cctki3_cctkGH_, \
                       i) \
  CCTK_LOOP1STR_INT(name, (cctki3_cctkGH_), \
                    i, \
                    cctki3_dummy_imin,cctki3_dummy_imax, 1) \

#define CCTK_ENDLOOP1_INT(name) \
  CCTK_ENDLOOP1STR_INT(name) \

#define CCTK_LOOP1STR_INT(name, cctki3_cctkGH_, \
                          i, \
                          imin,imax, cctki3_istr_) \
  do { \
    typedef int cctki3_loop1_int_##name; \
    cGH const *restrict const cctki3_cctkGH = (cctki3_cctkGH_); \
    if (cctki3_cctkGH->cctk_dim != 1) { \
      _Pragma("omp critical") \
      CCTK_WARN(CCTK_WARN_ABORT, \
                "The macro CCTK_LOOP1_INT can only be used in 1 dimensions"); \
    } \
    CCTK_INT cctki3_bndsize    [2]; \
    CCTK_INT cctki3_is_ghostbnd[2]; \
    CCTK_INT cctki3_is_symbnd  [2]; \
    CCTK_INT cctki3_is_physbnd [2]; \
    _Pragma("omp single copyprivate(cctki3_bndsize)") \
    GetBoundarySizesAndTypes \
      (cctki3_cctkGH, 2, cctki3_bndsize, cctki3_is_ghostbnd, cctki3_is_symbnd, cctki3_is_physbnd); \
    CCTK_LOOP1STR_INTERIOR(name##_int, \
                           cctki3_cctkGH, \
                           i, \
                           cctki3_bndsize[0], \
                           cctki3_bndsize[1], \
                           imin,imax, (cctki3_istr_)) { \

#define CCTK_ENDLOOP1STR_INT(name) \
    } CCTK_ENDLOOP1STR_INTERIOR(name##_int); \
    typedef cctki3_loop1_int_##name cctki3_ensure_proper_nesting CCTK_ATTRIBUTE_UNUSED; \
  } while (0) \



/* LOOP_BND */

#define CCTK_LOOP1_BND(name, cctki3_cctkGH_, \
                       i, \
                       ni) \
  CCTK_LOOP1STR_BND(name, (cctki3_cctkGH_), \
                    i, \
                    ni, \
                    cctki3_dummy_imin,cctki3_dummy_imax, 1) \

#define CCTK_ENDLOOP1_BND(name) \
  CCTK_ENDLOOP1STR_BND(name) \

#define CCTK_LOOP1STR_BND(name, cctki3_cctkGH_, \
                          i, \
                          ni, \
                          imin,imax, cctki3_istr_) \
  do { \
    typedef int cctki3_loop1_bnd_##name; \
    cGH const *restrict const cctki3_cctkGH = (cctki3_cctkGH_); \
    if (cctki3_cctkGH->cctk_dim != 1) { \
      _Pragma("omp critical") \
      CCTK_WARN(CCTK_WARN_ABORT, \
                "The macro CCTK_LOOP1_BND can only be used in 1 dimensions"); \
    } \
    CCTK_INT cctki3_bndsize    [2]; \
    CCTK_INT cctki3_is_ghostbnd[2]; \
    CCTK_INT cctki3_is_symbnd  [2]; \
    CCTK_INT cctki3_is_physbnd [2]; \
    _Pragma("omp single copyprivate(cctki3_bndsize, cctki3_is_physbnd)") \
    GetBoundarySizesAndTypes \
      (cctki3_cctkGH, 2, cctki3_bndsize, cctki3_is_ghostbnd, cctki3_is_symbnd, cctki3_is_physbnd); \
    CCTK_LOOP1STR_BOUNDARIES(name##_bnd, \
                             cctki3_cctkGH, \
                             i, \
                             ni, \
                             cctki3_bndsize[0], \
                             cctki3_bndsize[1], \
                             cctki3_is_physbnd[0], \
                             cctki3_is_physbnd[1], \
                             imin,imax, (cctki3_istr_)) { \

#define CCTK_ENDLOOP1STR_BND(name) \
    } CCTK_ENDLOOP1STR_BOUNDARIES(name##_bnd); \
    typedef cctki3_loop1_bnd_##name cctki3_ensure_proper_nesting CCTK_ATTRIBUTE_UNUSED; \
  } while (0) \



/* LOOP_INTBND */

#define CCTK_LOOP1_INTBND(name, cctki3_cctkGH_, \
                           i, \
                           ni) \
  CCTK_LOOP1STR_INTBND(name, (cctki3_cctkGH_), \
                        i, \
                        ni, \
                        cctki3_dummy_imin,cctki3_dummy_imax, 1) \

#define CCTK_ENDLOOP1_INTBND(name) \
  CCTK_ENDLOOP1STR_INTBND(name) \

#define CCTK_LOOP1STR_INTBND(name, cctki3_cctkGH_, \
                              i, \
                              ni, \
                              imin,imax, cctki3_istr_) \
  do { \
    typedef int cctki3_loop1_intbnd_##name; \
    cGH const *restrict const cctki3_cctkGH = (cctki3_cctkGH_); \
    if (cctki3_cctkGH->cctk_dim != 1) { \
      _Pragma("omp critical") \
      CCTK_WARN(CCTK_WARN_ABORT, \
                "The macro CCTK_LOOP1_INTBND can only be used in 1 dimensions"); \
    } \
    CCTK_INT cctki3_bndsize    [2]; \
    CCTK_INT cctki3_is_ghostbnd[2]; \
    CCTK_INT cctki3_is_symbnd  [2]; \
    CCTK_INT cctki3_is_physbnd [2]; \
    _Pragma("omp single copyprivate(cctki3_bndsize, cctki3_is_physbnd)") \
    GetBoundarySizesAndTypes \
      (cctki3_cctkGH, 2, cctki3_bndsize, cctki3_is_ghostbnd, cctki3_is_symbnd, cctki3_is_physbnd); \
    CCTK_LOOP1STR_INTBOUNDARIES(name##_intbnd, \
                                cctki3_cctkGH, \
                                i, \
                                ni, \
                                cctki3_bndsize[0], \
                                cctki3_bndsize[1], \
                                cctki3_is_physbnd[0], \
                                cctki3_is_physbnd[1], \
                                imin,imax, (cctki3_istr_)) { \

#define CCTK_ENDLOOP1STR_INTBND(name) \
    } CCTK_ENDLOOP1STR_INTBOUNDARIES(name##_intbnd); \
    typedef cctki3_loop1_intbnd_##name cctki3_ensure_proper_nesting CCTK_ATTRIBUTE_UNUSED; \
  } while (0) \

#endif /* #ifdef CCODE */



#ifdef FCODE

/* LOOP */

#define CCTK_LOOP1_NORMAL_DECLARE(name) \
   CCTK_LOOP1STR_NORMAL_DECLARE(name) \
   && integer :: name/**/0_dummy_imin, name/**/0_dummy_imax \

#define CCTK_LOOP1_NORMAL_OMP_PRIVATE(name) \
   CCTK_LOOP1STR_NORMAL_OMP_PRIVATE(name) \

#define CCTK_LOOP1_NORMAL(name, \
                          i, \
                          ni, \
                          cctki0_idir, \
                          cctki0_imin, \
                          cctki0_imax, \
                          cctki0_iash) \
   CCTK_LOOP1STR_NORMAL(name, \
                        i, \
                        ni, \
                        cctki0_idir, \
                        cctki0_imin, \
                        cctki0_imax, \
                        cctki0_iash, \
                        name/**/0_dummy_imin,name/**/0_dummy_imax, 1) \

#define CCTK_ENDLOOP1_NORMAL(name) \
   CCTK_ENDLOOP1STR_NORMAL(name) \

#define CCTK_LOOP1STR_NORMAL_DECLARE(name) \
   && integer :: name/**/0_idir \
   && integer :: name/**/0_imin \
   && integer :: name/**/0_imax \
   && integer :: name/**/0_iash \
   && integer :: name/**/0_istr \

#define CCTK_LOOP1STR_NORMAL_OMP_PRIVATE(name) \
   && !$omp private (i) \
   && !$omp private (ni) \

#define CCTK_LOOP1STR_NORMAL(name, \
                             i, \
                             ni, \
                             cctki0_idir, \
                             cctki0_imin, \
                             cctki0_imax, \
                             cctki0_iash, \
                             imin,imax, cctki0_istr) \
   && name/**/0_idir = cctki0_idir \
   && name/**/0_imin = cctki0_imin \
   && name/**/0_imax = cctki0_imax \
   && name/**/0_iash = cctki0_iash \
   && name/**/0_istr = cctki0_istr \
   && imin = name/**/0_imin \
   && imax = name/**/0_imax \
   && !$omp do \
   && do i = name/**/0_imin - modulo((imin), name/**/0_istr), name/**/0_imax, name/**/0_istr \
   &&    ni = 0 \
   &&    if (name/**/0_idir < 0) ni = i \
   &&    if (name/**/0_idir > 0) ni = name/**/0_imax+1-i \

#define CCTK_ENDLOOP1STR_NORMAL(name) \
   && end do \



#define CCTK_LOOP1_DECLARE(name) \
   CCTK_LOOP1STR_DECLARE(name) \
   && integer :: name/**/1_dummy_imin, name/**/1_dummy_imax \

#define CCTK_LOOP1_OMP_PRIVATE(name) \
   CCTK_LOOP1STR_OMP_PRIVATE(name) \

#define CCTK_LOOP1(name, \
                   i, \
                   cctki0_imin, \
                   cctki0_imax, \
                   cctki0_iash) \
   CCTK_LOOP1STR(name, \
                 i, \
                 cctki0_imin, \
                 cctki0_imax, \
                 cctki0_iash, \
                 name/**/1_dummy_imin,name/**/1_dummy_imax, 1) \

#define CCTK_ENDLOOP1(name) \
   CCTK_ENDLOOP1STR(name) \

#define CCTK_LOOP1STR_DECLARE(name) \
   CCTK_LOOP1STR_NORMAL_DECLARE(name) \
   && integer :: name/**/1_ni \

#define CCTK_LOOP1STR_OMP_PRIVATE(name) \
   CCTK_LOOP1STR_NORMAL_OMP_PRIVATE(name) \

#define CCTK_LOOP1STR(name, \
                      i, \
                      cctki1_imin, \
                      cctki1_imax, \
                      cctki1_iash, \
                      imin,imax, cctki1_istr) \
   CCTK_LOOP1STR_NORMAL(name, \
                        i, \
                        name/**/1_ni, \
                        0, \
                        cctki1_imin, \
                        cctki1_imax, \
                        cctki1_iash, \
                        imin,imax, cctki1_istr) \

#define CCTK_ENDLOOP1STR(name) \
   CCTK_ENDLOOP1STR_NORMAL(name) \



/* LOOP_INTERIOR */

#define CCTK_LOOP1_INTERIOR_DECLARE(name) \
   CCTK_LOOP1STR_INTERIOR_DECLARE(name) \
   && integer :: name/**/2_dummy_imin, name/**/2_dummy_imax \

#define CCTK_LOOP1_INTERIOR_OMP_PRIVATE(name) \
   CCTK_LOOP1STR_INTERIOR_OMP_PRIVATE(name) \

#define CCTK_LOOP1_INTERIOR(name, \
                            i, \
                            cctki2_iblo, \
                            cctki2_ibhi) \
   CCTK_LOOP1STR_INTERIOR(name, \
                          i, \
                          cctki2_iblo, \
                          cctki2_ibhi, \
                          name/**/2_dummy_imin,name/**/2_dummy_imax, 1) \

#define CCTK_ENDLOOP1_INTERIOR(name) \
   CCTK_ENDLOOP1STR_INTERIOR(name) \

#define CCTK_LOOP1STR_INTERIOR_DECLARE(name) \
   CCTK_LOOP1STR_DECLARE(name/**/_interior) \

#define CCTK_LOOP1STR_INTERIOR_OMP_PRIVATE(name) \
   CCTK_LOOP1STR_OMP_PRIVATE(name/**/_interior) \

#define CCTK_LOOP1STR_INTERIOR(name, \
                               i, \
                               cctki2_iblo, \
                               cctki2_ibhi, \
                               imin,imax, cctki2_istr) \
   CCTK_LOOP1STR(name/**/_interior, \
                 i, \
                 (cctki2_iblo)+1, \
                 cctk_lsh(1)-(cctki2_ibhi), \
                 cctk_ash(1), \
                 imin,imax, cctki2_istr) \

#define CCTK_ENDLOOP1STR_INTERIOR(name) \
   CCTK_ENDLOOP1STR(name/**/_interior) \



/* LOOP_BOUNDARIES */

#define CCTK_LOOP1_BOUNDARIES_DECLARE(name) \
   CCTK_LOOP1STR_BOUNDARIES_DECLARE(name) \
   && integer :: name/**/2_dummy_imin, name/**/2_dummy_imax \

#define CCTK_LOOP1_BOUNDARIES_OMP_PRIVATE(name) \
   CCTK_LOOP1STR_BOUNDARIES_OMP_PRIVATE(name) \

#define CCTK_LOOP1_BOUNDARIES(name, \
                              i, \
                              ni, \
                              cctki2_iblo, \
                              cctki2_ibhi, \
                              cctki2_ibboxlo, \
                              cctki2_ibboxhi) \
   CCTK_LOOP1STR_BOUNDARIES(name, \
                            i, \
                            ni, \
                            cctki2_iblo, \
                            cctki2_ibhi, \
                            cctki2_ibboxlo, \
                            cctki2_ibboxhi, \
                            name/**/2_dummy_imin,name/**/2_dummy_imax, 1) \

#define CCTK_ENDLOOP1_BOUNDARIES(name) \
   CCTK_ENDLOOP1STR_BOUNDARIES(name) \

#define CCTK_LOOP1STR_BOUNDARIES_DECLARE(name) \
   CCTK_LOOP1STR_NORMAL_DECLARE(name/**/_boundaries) \
   && integer :: name/**/2_blo(1), name/**/2_bhi(1) \
   && integer :: name/**/2_bboxlo(1), name/**/2_bboxhi(1) \
   && integer :: name/**/2_istr \
   && integer :: name/**/2_idir \
   && logical :: name/**/2_any_bbox \
   && integer :: name/**/2_bmin(1), name/**/2_bmax(1) \

#define CCTK_LOOP1STR_BOUNDARIES_OMP_PRIVATE(name) \
   CCTK_LOOP1STR_NORMAL_OMP_PRIVATE(name/**/_boundaries) \
   && !$omp private (name/**/2_bmin, name/**/2_bmax) \

#define CCTK_LOOP1STR_BOUNDARIES(name, \
                                 i, \
                                 ni, \
                                 cctki2_iblo, \
                                 cctki2_ibhi, \
                                 cctki2_ibboxlo, \
                                 cctki2_ibboxhi, \
                                 imin,imax, cctki2_istr) \
   && name/**/2_blo = (/ cctki2_iblo /) \
   && name/**/2_bhi = (/ cctki2_ibhi /) \
   && name/**/2_bboxlo = (/ cctki2_ibboxlo /) \
   && name/**/2_bboxhi = (/ cctki2_ibboxhi /) \
   && name/**/2_istr = (cctki2_istr) \
   && do name/**/2_idir=-1, +1 \
   &&     name/**/2_any_bbox = .false. \
   &&     if (name/**/2_idir<0) name/**/2_any_bbox = name/**/2_any_bbox .or. name/**/2_bboxlo(1) /= 0 \
   &&     if (name/**/2_idir>0) name/**/2_any_bbox = name/**/2_any_bbox .or. name/**/2_bboxhi(1) /= 0 \
   &&    if (name/**/2_any_bbox) then \
   &&       name/**/2_bmin(1) = name/**/2_blo(1)+1 \
   &&       if (name/**/2_idir<0) name/**/2_bmin(1) = 1 \
   &&       if (name/**/2_idir>0) name/**/2_bmin(1) = cctk_lsh(1) - name/**/2_bhi(1) \
   &&       name/**/2_bmax(1) = cctk_lsh(1) - name/**/2_bhi(1) \
   &&       if (name/**/2_idir<0) name/**/2_bmax(1) = name/**/2_blo(1) \
   &&       if (name/**/2_idir>0) name/**/2_bmax(1) = cctk_lsh(1) \
   &&       CCTK_LOOP1STR_NORMAL(name/**/_boundaries, \
                                 i, \
                                 ni, \
                                 name/**/2_idir, \
                                 name/**/2_bmin(1), \
                                 name/**/2_bmax(1), \
                                 cctk_ash(1), \
                                 imin,imax, name/**/2_istr) \

#define CCTK_ENDLOOP1STR_BOUNDARIES(name) \
            CCTK_ENDLOOP1STR_NORMAL(name/**/_boundaries) \
   &&    end if /* bbox */ \
   && end do /* dir */ \



/* LOOP_INTBOUNDARIES */

#define CCTK_LOOP1_INTBOUNDARIES_DECLARE(name) \
   CCTK_LOOP1STR_INTBOUNDARIES_DECLARE(name) \
   && integer :: name/**/2_dummy_imin, name/**/2_dummy_imax \

#define CCTK_LOOP1_INTBOUNDARIES_OMP_PRIVATE(name) \
   CCTK_LOOP1STR_INTBOUNDARIES_OMP_PRIVATE(name) \

#define CCTK_LOOP1_INTBOUNDARIES(name, \
                                 i, \
                                 ni, \
                                 cctki2_iblo, \
                                 cctki2_ibhi, \
                                 cctki2_ibboxlo, \
                                 cctki2_ibboxhi) \
   CCTK_LOOP1STR_INTBOUNDARIES(name, \
                               i, \
                               ni, \
                               cctki2_iblo, \
                               cctki2_ibhi, \
                               cctki2_ibboxlo, \
                               cctki2_ibboxhi, \
                               name/**/2_dummy_imin,name/**/2_dummy_max, 1) \

#define CCTK_ENDLOOP1_INTBOUNDARIES(name) \
   CCTK_ENDLOOP1STR_INTBOUNDARIES(name) \

#define CCTK_LOOP1STR_INTBOUNDARIES_DECLARE(name) \
   CCTK_LOOP1STR_NORMAL_DECLARE(name/**/_intboundaries) \
   && integer :: name/**/2_blo(1), name/**/2_bhi(1) \
   && integer :: name/**/2_bboxlo(1), name/**/2_bboxhi(1) \
   && integer :: name/**/2_istr \
   && integer :: name/**/2_idir \
   && logical :: name/**/2_any_bbox, name/**/2_all_bbox \
   && integer :: name/**/2_bmin(1), name/**/2_bmax(1) \

#define CCTK_LOOP1STR_INTBOUNDARIES_OMP_PRIVATE(name) \
   CCTK_LOOP1STR_NORMAL_OMP_PRIVATE(name/**/_intboundaries) \
   && !$omp private (name/**/2_any_bbox, name/**/2_all_bbox) \
   && !$omp private (name/**/2_bmin, name/**/2_bmax \

#define CCTK_LOOP1STR_INTBOUNDARIES(name, \
                                    i, \
                                    ni, \
                                    cctki2_iblo, \
                                    cctki2_ibhi, \
                                    cctki2_ibboxlo, \
                                    cctki2_ibboxhi, \
                                    imin,imax, cctki2_istr) \
   && name/**/2_blo = (/ cctki2_iblo /) \
   && name/**/2_bhi = (/ cctki2_ibhi /) \
   && name/**/2_bboxlo = (/ cctki2_ibboxlo /) \
   && name/**/2_bboxhi = (/ cctki2_ibboxhi /) \
   && name/**/2_istr = (cctki2_istr) \
   && do name/**/2_idir=-1, +1 \
   &&     name/**/2_any_bbox = .false. \
   &&     if (name/**/2_idir<0) name/**/2_any_bbox = name/**/2_any_bbox .or. name/**/2_bboxlo(1) /= 0 \
   &&     if (name/**/2_idir>0) name/**/2_any_bbox = name/**/2_any_bbox .or. name/**/2_bboxhi(1) /= 0 \
   &&     name/**/2_all_bbox = .true. \
   &&     if (name/**/2_idir<0) name/**/2_all_bbox = name/**/2_all_bbox .and. name/**/2_bboxlo(1) /= 0 \
   &&     if (name/**/2_idir>0) name/**/2_all_bbox = name/**/2_all_bbox .and. name/**/2_bboxhi(1) /= 0 \
   &&    if (name/**/2_all_bbox .and. name/**/2_any_bbox) then \
   &&       name/**/2_bmin(1) = name/**/2_blo(1)+1 \
   &&       if (name/**/2_idir<0) name/**/2_bmin(1) = 1 \
   &&       if (name/**/2_idir>0) name/**/2_bmin(1) = cctk_lsh(1) - name/**/2_bhi(1) \
   &&       name/**/2_bmax(1) = cctk_lsh(1) - name/**/2_bhi(1) \
   &&       if (name/**/2_idir<0) name/**/2_bmax(1) = name/**/2_blo(1) \
   &&       if (name/**/2_idir>0) name/**/2_bmax(1) = cctk_lsh(1) \
   &&       CCTK_LOOP1STR_NORMAL(name/**/_intboundaries, \
                                 i, \
                                 ni, \
                                 name/**/2_idir, \
                                 name/**/2_bmin(1), \
                                 name/**/2_bmax(1), \
                                 cctk_ash(1), \
                                 imin,imax, name/**/2_istr) \

#define CCTK_ENDLOOP1STR_INTBOUNDARIES(name) \
            CCTK_ENDLOOP1STR_NORMAL(name/**/_intboundaries) \
   &&    end if /* bbox */ \
   && end do /* dir */ \



/* LOOP_ALL */

#define CCTK_LOOP1_ALL_DECLARE(name) \
   CCTK_LOOP1STR_ALL_DECLARE(name) \
   && integer :: name/**/3_dummy_imin, name/**/3_dummy_imax \

#define CCTK_LOOP1_ALL_OMP_PRIVATE(name) \
   CCTK_LOOP1STR_ALL_OMP_PRIVATE(name) \

#define CCTK_LOOP1_ALL(name, \
                       i) \
   CCTK_LOOP1STR_ALL(name, \
                     i, \
                     name/**/3_dummy_imin,name/**/3_dummy_imax, 1) \

#define CCTK_ENDLOOP1_ALL(name) \
   CCTK_ENDLOOP1STR_ALL(name) \

#define CCTK_LOOP1STR_ALL_DECLARE(name) \
   CCTK_LOOP1STR_DECLARE(name/**/_all) \

#define CCTK_LOOP1STR_ALL_OMP_PRIVATE(name) \
   CCTK_LOOP1STR_OMP_PRIVATE(name/**/_all) \

#define CCTK_LOOP1STR_ALL(name, \
                          i, \
                          imin,imax, cctki3_istr) \
   CCTK_LOOP1STR(name/**/_all, \
                 i, \
                 1, \
                 cctk_lsh(1), \
                 cctk_ash(1), \
                 imin,imax, cctki3_istr) \

#define CCTK_ENDLOOP1STR_ALL(name) \
   CCTK_ENDLOOP1STR(name/**/_all) \



/* LOOP_INT */

#define CCTK_LOOP1_INT_DECLARE(name) \
   CCTK_LOOP1STR_INT_DECLARE(name) \
   && integer :: name/**/3_dummy_imin, name/**/3_dummy_imax \

#define CCTK_LOOP1_INT_OMP_PRIVATE(name) \
   CCTK_LOOP1STR_INT_OMP_PRIVATE(name) \

#define CCTK_LOOP1_INT(name, \
                        i) \
   CCTK_LOOP1STR_INT(name, \
                      i, \
                      name/**/3_dummy_imin,name/**/3_dummy_imax, 1) \

#define CCTK_ENDLOOP1_INT(name) \
   CCTK_ENDLOOP1STR_INT(name) \

#define CCTK_LOOP1STR_INT_DECLARE(name) \
   CCTK_LOOP1STR_INTERIOR_DECLARE(name/**/_int) \
   && CCTK_INT, parameter :: name/**/3_isize = 2 \
   && CCTK_INT :: name/**/3_bndsize    (2) \
   && CCTK_INT :: name/**/3_is_ghostbnd(2) \
   && CCTK_INT :: name/**/3_is_symbnd  (2) \
   && CCTK_INT :: name/**/3_is_physbnd (2) \
   && CCTK_INT :: name/**/3_ierr \

#define CCTK_LOOP1STR_INT_OMP_PRIVATE(name) \
   CCTK_LOOP1STR_INTERIOR_OMP_PRIVATE(name/**/_int) \

#define CCTK_LOOP1STR_INT(name, \
                          i, \
                          imin,imax, cctki3_istr) \
   && !$omp single \
   && name/**/3_ierr = GetBoundarySizesAndTypes \
         (cctkGH, name/**/3_isize, name/**/3_bndsize, name/**/3_is_ghostbnd, name/**/3_is_symbnd, name/**/3_is_physbnd) \
   && !$omp end single copyprivate(name/**/3_bndsize) \
   && CCTK_LOOP1STR_INTERIOR(name/**/_int, \
                             i, \
                             int(name/**/3_bndsize(1+1)), \
                             int(name/**/3_bndsize(2)), \
                             imin,imax, (cctki3_istr)) \

#define CCTK_ENDLOOP1STR_INT(name) \
      CCTK_ENDLOOP1STR_INTERIOR(name/**/_int) \



/* LOOP_BND */

#define CCTK_LOOP1_BND_DECLARE(name) \
   CCTK_LOOP1STR_BND_DECLARE(name) \
   && integer :: name/**/3_dummy_imin, name/**/3_dummy_imax \

#define CCTK_LOOP1_BND_OMP_PRIVATE(name) \
   CCTK_LOOP1STR_BND_OMP_PRIVATE(name) \

#define CCTK_LOOP1_BND(name, \
                       i, \
                       ni) \
   CCTK_LOOP1STR_BND(name, \
                     i, \
                     ni, \
                     name/**/3_dummy_imin,name/**/3_dummy_imax, 1) \

#define CCTK_ENDLOOP1_BND(name) \
   CCTK_ENDLOOP1STR_BND(name) \

#define CCTK_LOOP1STR_BND_DECLARE(name) \
   CCTK_LOOP1STR_BOUNDARIES_DECLARE(name/**/_bnd) \
   && CCTK_INT, parameter :: name/**/3_isize = 2 \
   && CCTK_INT :: name/**/3_bndsize    (2) \
   && CCTK_INT :: name/**/3_is_ghostbnd(2) \
   && CCTK_INT :: name/**/3_is_symbnd  (2) \
   && CCTK_INT :: name/**/3_is_physbnd (2) \
   && CCTK_INT :: name/**/3_ierr \

#define CCTK_LOOP1STR_BND_OMP_PRIVATE(name) \
   CCTK_LOOP1STR_BOUNDARIES_OMP_PRIVATE(name/**/_bnd) \

#define CCTK_LOOP1STR_BND(name, \
                          i, \
                          ni, \
                          imin,imax, cctki3_istr) \
   && !$omp single \
   && name/**/3_ierr = GetBoundarySizesAndTypes \
         (cctkGH, name/**/3_isize, name/**/3_bndsize, name/**/3_is_ghostbnd, name/**/3_is_symbnd, name/**/3_is_physbnd) \
   && !$omp end single copyprivate(name/**/3_bndsize, name/**/3_is_physbnd) \
   && CCTK_LOOP1STR_BOUNDARIES(name/**/_bnd, \
                               i, \
                               ni, \
                               int(name/**/3_bndsize(1))+1, \
                               int(name/**/3_bndsize(2)), \
                               int(name/**/3_is_physbnd(1)), \
                               int(name/**/3_is_physbnd(2)), \
                               imin,imax, (cctki3_istr)) \

#define CCTK_ENDLOOP1STR_BND(name) \
      CCTK_ENDLOOP1STR_BOUNDARIES(name/**/_bnd) \






/* LOOP_INTBND */

#define CCTK_LOOP1_INTBND_DECLARE(name) \
   CCTK_LOOP1STR_INTBND_DECLARE(name) \
   && integer :: name/**/3_dummy_imin, name/**/3_dummy_imax \

#define CCTK_LOOP1_INTBND_OMP_PRIVATE(name) \
   CCTK_LOOP1STR_INTBND_OMP_PRIVATE(name) \

#define CCTK_LOOP1_INTBND(name, \
                          i, \
                          ni) \
   CCTK_LOOP1STR_INTBND(name, \
                        i, \
                        ni, \
                        name/**/3_dummy_imin,name/**/3_dummy_imax, 1) \

#define CCTK_ENDLOOP1_INTBND(name) \
   CCTK_ENDLOOP1STR_INTBND(name) \

#define CCTK_LOOP1STR_INTBND_DECLARE(name) \
   CCTK_LOOP1STR_INTBOUNDARIES_DECLARE(name/**/_bnd) \
   && CCTK_INT, parameter :: name/**/3_isize = 2 \
   && CCTK_INT :: name/**/3_bndsize    (2) \
   && CCTK_INT :: name/**/3_is_ghostbnd(2) \
   && CCTK_INT :: name/**/3_is_symbnd  (2) \
   && CCTK_INT :: name/**/3_is_physbnd (2) \
   && CCTK_INT :: name/**/3_ierr \

#define CCTK_LOOP1STR_INTBND_OMP_PRIVATE(name) \
   CCTK_LOOP1STR_INTBOUNDARIES_OMP_PRIVATE(name/**/_bnd) \

#define CCTK_LOOP1STR_INTBND(name, \
                             i, \
                             ni, \
                             imin,imax, cctki3_istr) \
   && !$omp single \
   && name/**/3_ierr = GetBoundarySizesAndTypes \
         (cctkGH, name/**/3_isize, name/**/3_bndsize, name/**/3_is_ghostbnd, name/**/3_is_symbnd, name/**/3_is_physbnd) \
   && !$omp end single copyprivate(name/**/3_bndsize, name/**/3_is_physbnd) \
   && CCTK_LOOP1STR_INTBOUNDARIES(name/**/_bnd, \
                                  i, \
                                  ni, \
                                  int(name/**/3_bndsize(1+1)), \
                                  int(name/**/3_bndsize(2)), \
                                  int(name/**/3_is_physbnd(1)), \
                                  int(name/**/3_is_physbnd(2)), \
                                  imin,imax, (cctki3_istr)) \

#define CCTK_ENDLOOP1STR_INTBND(name) \
      CCTK_ENDLOOP1STR_INTBOUNDARIES(name/**/_bnd) \

#endif /* #ifdef FCODE */



/* 2D */

#ifdef CCODE

/* LOOP */

#define CCTK_LOOP2_NORMAL(name, \
                          i,j, \
                          ni,nj, \
                          cctki0_idir_,cctki0_jdir_, \
                          cctki0_imin_,cctki0_jmin_, \
                          cctki0_imax_,cctki0_jmax_, \
                          cctki0_iash_,cctki0_jash_) \
  CCTK_LOOP2STR_NORMAL(name, \
                       i,j, \
                       ni,nj, \
                       (cctki0_idir_),(cctki0_jdir_), \
                       (cctki0_imin_),(cctki0_jmin_), \
                       (cctki0_imax_),(cctki0_jmax_), \
                       (cctki0_iash_),(cctki0_jash_), \
                       cctki0_dummy_imin,cctki0_dummy_imax, 1) \

#define CCTK_ENDLOOP2_NORMAL(name) \
  CCTK_ENDLOOP2STR_NORMAL(name) \

#define CCTK_LOOP2STR_NORMAL(name, \
                             i,j, \
                             ni,nj, \
                             cctki0_idir_,cctki0_jdir_, \
                             cctki0_imin_,cctki0_jmin_, \
                             cctki0_imax_,cctki0_jmax_, \
                             cctki0_iash_,cctki0_jash_, \
                             imin,imax, cctki0_istr_) \
  do { \
    typedef int cctki0_loop2_normal_##name; \
    const int cctki0_idir = (cctki0_idir_); \
    const int cctki0_jdir = (cctki0_jdir_); \
    const int cctki0_imin = (cctki0_imin_); \
    const int cctki0_jmin = (cctki0_jmin_); \
    const int cctki0_imax = (cctki0_imax_); \
    const int cctki0_jmax = (cctki0_jmax_); \
    const int cctki0_iash CCTK_ATTRIBUTE_UNUSED = (cctki0_iash_); \
    const int cctki0_jash CCTK_ATTRIBUTE_UNUSED = (cctki0_jash_); \
    const int cctki0_istr = (cctki0_istr_); \
    assert(cctki0_istr>0 && (cctki0_istr & (cctki0_istr-1)) == 0); \
    const int imin CCTK_ATTRIBUTE_UNUSED = cctki0_imin; \
    const int imax CCTK_ATTRIBUTE_UNUSED = cctki0_imax; \
    CCTK_PRAGMA_OMP_FOR_COLLAPSE_1 \
    for (int j=cctki0_jmin; j<cctki0_jmax; ++j) { \
     \
    const int cctki0_ioff = (cctki0_imin+cctki0_iash*(j)) & (cctki0_istr-1); \
     \
    for (int i=cctki0_imin-cctki0_ioff; i<cctki0_imax; i+=cctki0_istr) { \
      const int ni CCTK_ATTRIBUTE_UNUSED = cctki0_idir<0 ? i+1 : cctki0_idir==0 ? 0 : cctki0_imax-i; \
      const int nj CCTK_ATTRIBUTE_UNUSED = cctki0_jdir<0 ? j+1 : cctki0_jdir==0 ? 0 : cctki0_jmax-j; \
      { \

#define CCTK_ENDLOOP2STR_NORMAL(name) \
      } \
    } \
    } \
    typedef cctki0_loop2_normal_##name cctki0_ensure_proper_nesting CCTK_ATTRIBUTE_UNUSED; \
  } while (0) \



#define CCTK_LOOP2(name, \
                   i,j, \
                   cctki1_imin_,cctki1_jmin_, \
                   cctki1_imax_,cctki1_jmax_, \
                   cctki1_iash_,cctki1_jash_) \
  CCTK_LOOP2STR(name, \
                i,j, \
                (cctki1_imin_),(cctki1_jmin_), \
                (cctki1_imax_),(cctki1_jmax_), \
                (cctki1_iash_),(cctki1_jash_), \
                cctki1_dummy_imin,cctki1_dummy_imax, 1) \

#define CCTK_ENDLOOP2(name) \
  CCTK_ENDLOOP2STR(name) \

#define CCTK_LOOP2STR(name, \
                      i,j, \
                      cctki1_imin_,cctki1_jmin_, \
                      cctki1_imax_,cctki1_jmax_, \
                      cctki1_iash_,cctki1_jash_, \
                      imin,imax, cctki1_istr_) \
  CCTK_LOOP2STR_NORMAL(name, \
                       i,j, \
                       cctki1_ni,cctki1_nj, \
                       0,0, \
                       (cctki1_imin_),(cctki1_jmin_), \
                       (cctki1_imax_),(cctki1_jmax_), \
                       (cctki1_iash_),(cctki1_jash_), \
                       imin,imax, (cctki1_istr_)) \

#define CCTK_ENDLOOP2STR(name) \
  CCTK_ENDLOOP2STR_NORMAL(name) \



/* LOOP_INTERIOR */

#define CCTK_LOOP2_INTERIOR(name, cctki2_cctkGH_, \
                            i,j, \
                            cctki2_iblo_,cctki2_jblo_, \
                            cctki2_ibhi_,cctki2_jbhi_) \
  CCTK_LOOP2STR_INTERIOR(name, (cctki2_cctkGH_), \
                         i,j, \
                         (cctki2_iblo_),(cctki2_jblo_), \
                         (cctki2_ibhi_),(cctki2_jbhi_), \
                         cctki2_dummy_imin,cctki2_dummy_imax, 1) \

#define CCTK_ENDLOOP2_INTERIOR(name) \
  CCTK_ENDLOOP2STR_INTERIOR(name) \

#define CCTK_LOOP2STR_INTERIOR(name, cctki2_cctkGH_, \
                               i,j, \
                               cctki2_iblo_,cctki2_jblo_, \
                               cctki2_ibhi_,cctki2_jbhi_, \
                               imin,imax, cctki2_istr_) \
  do { \
    typedef int cctki2_loop2_interior_##name; \
    cGH const *restrict const cctki2_cctkGH = (cctki2_cctkGH_); \
    if (cctki2_cctkGH->cctk_dim != 2) { \
      _Pragma("omp critical") \
      CCTK_WARN(CCTK_WARN_ABORT, \
                "The macro CCTK_LOOP2_INTERIOR can only be used in 2 dimensions"); \
    } \
    CCTK_LOOP2STR(name##_interior, \
                  i,j, \
                  (cctki2_iblo_),(cctki2_jblo_), \
                  cctki2_cctkGH->cctk_lsh[0]-(cctki2_ibhi_), \
                  cctki2_cctkGH->cctk_lsh[1]-(cctki2_jbhi_), \
                  cctki2_cctkGH->cctk_ash[0], \
                  cctki2_cctkGH->cctk_ash[1], \
                  imin,imax, (cctki2_istr_)) { \

#define CCTK_ENDLOOP2STR_INTERIOR(name) \
    } CCTK_ENDLOOP2STR(name##_interior); \
    typedef cctki2_loop2_interior_##name cctki2_ensure_proper_nesting CCTK_ATTRIBUTE_UNUSED; \
  } while(0) \



/* LOOP_BOUNDARIES */

#define CCTK_LOOP2_BOUNDARIES(name, cctki2_cctkGH_, \
                              i,j, \
                              ni,nj, \
                              cctki2_iblo_,cctki2_jblo_, \
                              cctki2_ibhi_,cctki2_jbhi_, \
                              cctki2_ibboxlo_,cctki2_jbboxlo_, \
                              cctki2_ibboxhi_,cctki2_jbboxhi_) \
  CCTK_LOOP2STR_BOUNDARIES(name, (cctki2_cctkGH_), \
                           i,j, \
                           ni,nj, \
                           (cctki2_iblo_),(cctki2_jblo_), \
                           (cctki2_ibhi_),(cctki2_jbhi_), \
                           (cctki2_ibboxlo_),(cctki2_jbboxlo_), \
                           (cctki2_ibboxhi_),(cctki2_jbboxhi_), \
                           cctki2_dummy_imin,cctki2_dummy_imax, 1) \

#define CCTK_ENDLOOP2_BOUNDARIES(name) \
  CCTK_ENDLOOP2STR_BOUNDARIES(name) \

#define CCTK_LOOP2STR_BOUNDARIES(name, cctki2_cctkGH_, \
                                 i,j, \
                                 ni,nj, \
                                 cctki2_iblo_,cctki2_jblo_, \
                                 cctki2_ibhi_,cctki2_jbhi_, \
                                 cctki2_ibboxlo_,cctki2_jbboxlo_, \
                                 cctki2_ibboxhi_,cctki2_jbboxhi_, \
                                 imin,imax, cctki2_istr_) \
  do { \
    typedef int cctki2_loop2_boundaries_##name; \
    cGH const *restrict const cctki2_cctkGH = (cctki2_cctkGH_); \
    if (cctki2_cctkGH->cctk_dim != 2) { \
      _Pragma("omp critical") \
      CCTK_WARN(CCTK_WARN_ABORT, \
                "The macro CCTK_LOOP2_BOUNDARIES can only be used in 2 dimensions"); \
    } \
    const int cctki2_blo[] = { (cctki2_iblo_),(cctki2_jblo_) }; \
    const int cctki2_bhi[] = { (cctki2_ibhi_),(cctki2_jbhi_) }; \
    const int cctki2_bbox[] = { (cctki2_ibboxlo_), (cctki2_ibboxhi_),(cctki2_jbboxlo_), (cctki2_jbboxhi_) }; \
    const int cctki2_lsh[] = { cctki2_cctkGH->cctk_lsh[0],cctki2_cctkGH->cctk_lsh[1] }; \
    const int cctki2_istr CCTK_ATTRIBUTE_UNUSED = (cctki2_istr_); \
    for (int cctki2_jdir=-1; cctki2_jdir<=+1; ++cctki2_jdir) { \
    for (int cctki2_idir=-1; cctki2_idir<=+1; ++cctki2_idir) { \
      const int cctki2_any_bbox = \
        (cctki2_idir<0 ? cctki2_bbox[0] : 0) || (cctki2_idir>0 ? cctki2_bbox[1] : 0) || \
        (cctki2_jdir<0 ? cctki2_bbox[2] : 0) || (cctki2_jdir>0 ? cctki2_bbox[3] : 0); \
      if (cctki2_any_bbox) { \
        const int cctki2_bmin[] = { \
          cctki2_idir<0 ? 0 : cctki2_idir==0 ? cctki2_blo[0] : cctki2_lsh[0] - cctki2_bhi[0], \
          cctki2_jdir<0 ? 0 : cctki2_jdir==0 ? cctki2_blo[1] : cctki2_lsh[1] - cctki2_bhi[1], \
        }; \
        const int cctki2_bmax[] = { \
          cctki2_idir<0 ? cctki2_blo[0] : cctki2_idir==0 ? cctki2_lsh[0] - cctki2_bhi[0] : cctki2_lsh[0], \
          cctki2_jdir<0 ? cctki2_blo[1] : cctki2_jdir==0 ? cctki2_lsh[1] - cctki2_bhi[1] : cctki2_lsh[1], \
        }; \
        CCTK_LOOP2STR_NORMAL(name##_boundaries, \
                             i,j, \
                             ni,nj, \
                             cctki2_idir,cctki2_jdir, \
                             cctki2_bmin[0],cctki2_bmin[1], \
                             cctki2_bmax[0],cctki2_bmax[1], \
                             cctki2_cctkGH->cctk_ash[0], \
                             cctki2_cctkGH->cctk_ash[1], \
                             imin,imax, cctki2_istr) { \

#define CCTK_ENDLOOP2STR_BOUNDARIES(name) \
        } CCTK_ENDLOOP2STR_NORMAL(name##_boundaries); \
      } /* if bbox */ \
    } /* for dir */ \
    } /* for dir */ \
    typedef cctki2_loop2_boundaries_##name cctki2_ensure_proper_nesting CCTK_ATTRIBUTE_UNUSED; \
  } while (0) \



/* LOOP_INTBOUNDARIES */

#define CCTK_LOOP2_INTBOUNDARIES(name, cctki2_cctkGH_, \
                                 i,j, \
                                 ni,nj, \
                                 cctki2_iblo_,cctki2_jblo_, \
                                 cctki2_ibhi_,cctki2_jbhi_, \
                                 cctki2_ibboxlo_,cctki2_jbboxlo_, \
                                 cctki2_ibboxhi_,cctki2_jbboxhi_) \
  CCTK_LOOP2STR_INTBOUNDARIES(name, (cctki2_cctkGH_), \
                              i,j, \
                              ni,nj, \
                              (cctki2_iblo_),(cctki2_jblo_), \
                              (cctki2_ibhi_),(cctki2_jbhi_), \
                              (cctki2_ibboxlo_),(cctki2_jbboxlo_), \
                              (cctki2_ibboxhi_),(cctki2_jbboxhi_), \
                              cctki2_dummy_imin,cctki2_dummy_imax, 1) \

#define CCTK_ENDLOOP2_INTBOUNDARIES(name) \
  CCTK_ENDLOOP2STR_INTBOUNDARIES(name) \

#define CCTK_LOOP2STR_INTBOUNDARIES(name, cctki2_cctkGH_, \
                                    i,j, \
                                    ni,nj, \
                                    cctki2_iblo_,cctki2_jblo_, \
                                    cctki2_ibhi_,cctki2_jbhi_, \
                                    cctki2_ibboxlo_,cctki2_jbboxlo_, \
                                    cctki2_ibboxhi_,cctki2_jbboxhi_, \
                                    imin,imax, cctki2_istr_) \
  do { \
    typedef int cctki2_loop2_intboundaries_##name; \
    cGH const *restrict const cctki2_cctkGH = (cctki2_cctkGH_); \
    if (cctki2_cctkGH->cctk_dim != 2) { \
      _Pragma("omp critical") \
      CCTK_WARN(CCTK_WARN_ABORT, \
                "The macro CCTK_LOOP2_INTBOUNDARIES can only be used in 2 dimensions"); \
    } \
    const int cctki2_blo[] = { (cctki2_iblo_),(cctki2_jblo_) }; \
    const int cctki2_bhi[] = { (cctki2_ibhi_),(cctki2_jbhi_) }; \
    const int cctki2_bbox[] = { (cctki2_ibboxlo_), (cctki2_ibboxhi_),(cctki2_jbboxlo_), (cctki2_jbboxhi_) }; \
    const int cctki2_lsh[] = { cctki2_cctkGH->cctk_lsh[0],cctki2_cctkGH->cctk_lsh[1] }; \
    const int cctki2_istr CCTK_ATTRIBUTE_UNUSED = (cctki2_istr_); \
    for (int cctki2_jdir=-1; cctki2_jdir<=+1; ++cctki2_jdir) { \
    for (int cctki2_idir=-1; cctki2_idir<=+1; ++cctki2_idir) { \
      const int cctki2_any_bbox = \
        (cctki2_idir<0 ? cctki2_bbox[0] : 0) || (cctki2_idir>0 ? cctki2_bbox[1] : 0) || \
        (cctki2_jdir<0 ? cctki2_bbox[2] : 0) || (cctki2_jdir>0 ? cctki2_bbox[3] : 0); \
      const int cctki2_all_bbox = \
        (cctki2_idir<0 ? cctki2_bbox[0] : 1) && (cctki2_idir>0 ? cctki2_bbox[1] : 1) && \
        (cctki2_jdir<0 ? cctki2_bbox[2] : 1) && (cctki2_jdir>0 ? cctki2_bbox[3] : 1); \
      if (cctki2_all_bbox && cctki2_any_bbox) { \
        const int cctki2_bmin[] = { \
          cctki2_idir<0 ? 0 : cctki2_idir==0 ? cctki2_blo[0] : cctki2_lsh[0] - cctki2_bhi[0], \
          cctki2_jdir<0 ? 0 : cctki2_jdir==0 ? cctki2_blo[1] : cctki2_lsh[1] - cctki2_bhi[1], \
        }; \
        const int cctki2_bmax[] = { \
          cctki2_idir<0 ? cctki2_blo[0] : cctki2_idir==0 ? cctki2_lsh[0] - cctki2_bhi[0] : cctki2_lsh[0], \
          cctki2_jdir<0 ? cctki2_blo[1] : cctki2_jdir==0 ? cctki2_lsh[1] - cctki2_bhi[1] : cctki2_lsh[1], \
        }; \
        CCTK_LOOP2STR_NORMAL(name##_intboundaries, \
                             i,j, \
                             ni,nj, \
                             cctki2_idir,cctki2_jdir, \
                             cctki2_bmin[0],cctki2_bmin[1], \
                             cctki2_bmax[0],cctki2_bmax[1], \
                             cctki2_cctkGH->cctk_ash[0], \
                             cctki2_cctkGH->cctk_ash[1], \
                             imin,imax, cctki2_istr) { \

#define CCTK_ENDLOOP2STR_INTBOUNDARIES(name) \
        } CCTK_ENDLOOP2STR_NORMAL(name##_intboundaries); \
      } /* if bbox */ \
    } /* for dir */ \
    } /* for dir */ \
    typedef cctki2_loop2_intboundaries_##name cctki2_ensure_proper_nesting CCTK_ATTRIBUTE_UNUSED; \
  } while (0) \



/* LOOP_ALL */

#define CCTK_LOOP2_ALL(name, cctki3_cctkGH_, \
                       i,j) \
  CCTK_LOOP2STR_ALL(name, (cctki3_cctkGH_), \
                    i,j, \
                    cctki3_dummy_imin,cctki3_dummy_imax, 1) \

#define CCTK_ENDLOOP2_ALL(name) \
  CCTK_ENDLOOP2STR_ALL(name) \

#define CCTK_LOOP2STR_ALL(name, cctki3_cctkGH_, \
                          i,j, \
                          imin,imax, cctki3_istr_) \
  do { \
    typedef int cctki3_loop2_all_##name; \
    cGH const *restrict const cctki3_cctkGH = (cctki3_cctkGH_); \
    if (cctki3_cctkGH->cctk_dim != 2) { \
      _Pragma("omp critical") \
      CCTK_WARN(CCTK_WARN_ABORT, \
                "The macro CCTK_LOOP2_ALL can only be used in 2 dimensions"); \
    } \
    CCTK_LOOP2STR(name##_all, \
                  i,j, \
                  0,0, \
                  cctki3_cctkGH->cctk_lsh[0], \
                  cctki3_cctkGH->cctk_lsh[1], \
                  cctki3_cctkGH->cctk_ash[0], \
                  cctki3_cctkGH->cctk_ash[1], \
                  imin,imax, (cctki3_istr_)) { \

#define CCTK_ENDLOOP2STR_ALL(name) \
    } CCTK_ENDLOOP2STR(name##_all); \
    typedef cctki3_loop2_all_##name cctki3_ensure_proper_nesting CCTK_ATTRIBUTE_UNUSED; \
  } while (0) \



/* LOOP_INT */

#define CCTK_LOOP2_INT(name, cctki3_cctkGH_, \
                       i,j) \
  CCTK_LOOP2STR_INT(name, (cctki3_cctkGH_), \
                    i,j, \
                    cctki3_dummy_imin,cctki3_dummy_imax, 1) \

#define CCTK_ENDLOOP2_INT(name) \
  CCTK_ENDLOOP2STR_INT(name) \

#define CCTK_LOOP2STR_INT(name, cctki3_cctkGH_, \
                          i,j, \
                          imin,imax, cctki3_istr_) \
  do { \
    typedef int cctki3_loop2_int_##name; \
    cGH const *restrict const cctki3_cctkGH = (cctki3_cctkGH_); \
    if (cctki3_cctkGH->cctk_dim != 2) { \
      _Pragma("omp critical") \
      CCTK_WARN(CCTK_WARN_ABORT, \
                "The macro CCTK_LOOP2_INT can only be used in 2 dimensions"); \
    } \
    CCTK_INT cctki3_bndsize    [4]; \
    CCTK_INT cctki3_is_ghostbnd[4]; \
    CCTK_INT cctki3_is_symbnd  [4]; \
    CCTK_INT cctki3_is_physbnd [4]; \
    _Pragma("omp single copyprivate(cctki3_bndsize)") \
    GetBoundarySizesAndTypes \
      (cctki3_cctkGH, 4, cctki3_bndsize, cctki3_is_ghostbnd, cctki3_is_symbnd, cctki3_is_physbnd); \
    CCTK_LOOP2STR_INTERIOR(name##_int, \
                           cctki3_cctkGH, \
                           i,j, \
                           cctki3_bndsize[0],cctki3_bndsize[2], \
                           cctki3_bndsize[1],cctki3_bndsize[3], \
                           imin,imax, (cctki3_istr_)) { \

#define CCTK_ENDLOOP2STR_INT(name) \
    } CCTK_ENDLOOP2STR_INTERIOR(name##_int); \
    typedef cctki3_loop2_int_##name cctki3_ensure_proper_nesting CCTK_ATTRIBUTE_UNUSED; \
  } while (0) \



/* LOOP_BND */

#define CCTK_LOOP2_BND(name, cctki3_cctkGH_, \
                       i,j, \
                       ni,nj) \
  CCTK_LOOP2STR_BND(name, (cctki3_cctkGH_), \
                    i,j, \
                    ni,nj, \
                    cctki3_dummy_imin,cctki3_dummy_imax, 1) \

#define CCTK_ENDLOOP2_BND(name) \
  CCTK_ENDLOOP2STR_BND(name) \

#define CCTK_LOOP2STR_BND(name, cctki3_cctkGH_, \
                          i,j, \
                          ni,nj, \
                          imin,imax, cctki3_istr_) \
  do { \
    typedef int cctki3_loop2_bnd_##name; \
    cGH const *restrict const cctki3_cctkGH = (cctki3_cctkGH_); \
    if (cctki3_cctkGH->cctk_dim != 2) { \
      _Pragma("omp critical") \
      CCTK_WARN(CCTK_WARN_ABORT, \
                "The macro CCTK_LOOP2_BND can only be used in 2 dimensions"); \
    } \
    CCTK_INT cctki3_bndsize    [4]; \
    CCTK_INT cctki3_is_ghostbnd[4]; \
    CCTK_INT cctki3_is_symbnd  [4]; \
    CCTK_INT cctki3_is_physbnd [4]; \
    _Pragma("omp single copyprivate(cctki3_bndsize, cctki3_is_physbnd)") \
    GetBoundarySizesAndTypes \
      (cctki3_cctkGH, 4, cctki3_bndsize, cctki3_is_ghostbnd, cctki3_is_symbnd, cctki3_is_physbnd); \
    CCTK_LOOP2STR_BOUNDARIES(name##_bnd, \
                             cctki3_cctkGH, \
                             i,j, \
                             ni,nj, \
                             cctki3_bndsize[0],cctki3_bndsize[2], \
                             cctki3_bndsize[1],cctki3_bndsize[3], \
                             cctki3_is_physbnd[0],cctki3_is_physbnd[2], \
                             cctki3_is_physbnd[1],cctki3_is_physbnd[3], \
                             imin,imax, (cctki3_istr_)) { \

#define CCTK_ENDLOOP2STR_BND(name) \
    } CCTK_ENDLOOP2STR_BOUNDARIES(name##_bnd); \
    typedef cctki3_loop2_bnd_##name cctki3_ensure_proper_nesting CCTK_ATTRIBUTE_UNUSED; \
  } while (0) \



/* LOOP_INTBND */

#define CCTK_LOOP2_INTBND(name, cctki3_cctkGH_, \
                           i,j, \
                           ni,nj) \
  CCTK_LOOP2STR_INTBND(name, (cctki3_cctkGH_), \
                        i,j, \
                        ni,nj, \
                        cctki3_dummy_imin,cctki3_dummy_imax, 1) \

#define CCTK_ENDLOOP2_INTBND(name) \
  CCTK_ENDLOOP2STR_INTBND(name) \

#define CCTK_LOOP2STR_INTBND(name, cctki3_cctkGH_, \
                              i,j, \
                              ni,nj, \
                              imin,imax, cctki3_istr_) \
  do { \
    typedef int cctki3_loop2_intbnd_##name; \
    cGH const *restrict const cctki3_cctkGH = (cctki3_cctkGH_); \
    if (cctki3_cctkGH->cctk_dim != 2) { \
      _Pragma("omp critical") \
      CCTK_WARN(CCTK_WARN_ABORT, \
                "The macro CCTK_LOOP2_INTBND can only be used in 2 dimensions"); \
    } \
    CCTK_INT cctki3_bndsize    [4]; \
    CCTK_INT cctki3_is_ghostbnd[4]; \
    CCTK_INT cctki3_is_symbnd  [4]; \
    CCTK_INT cctki3_is_physbnd [4]; \
    _Pragma("omp single copyprivate(cctki3_bndsize, cctki3_is_physbnd)") \
    GetBoundarySizesAndTypes \
      (cctki3_cctkGH, 4, cctki3_bndsize, cctki3_is_ghostbnd, cctki3_is_symbnd, cctki3_is_physbnd); \
    CCTK_LOOP2STR_INTBOUNDARIES(name##_intbnd, \
                                cctki3_cctkGH, \
                                i,j, \
                                ni,nj, \
                                cctki3_bndsize[0],cctki3_bndsize[2], \
                                cctki3_bndsize[1],cctki3_bndsize[3], \
                                cctki3_is_physbnd[0],cctki3_is_physbnd[2], \
                                cctki3_is_physbnd[1],cctki3_is_physbnd[3], \
                                imin,imax, (cctki3_istr_)) { \

#define CCTK_ENDLOOP2STR_INTBND(name) \
    } CCTK_ENDLOOP2STR_INTBOUNDARIES(name##_intbnd); \
    typedef cctki3_loop2_intbnd_##name cctki3_ensure_proper_nesting CCTK_ATTRIBUTE_UNUSED; \
  } while (0) \

#endif /* #ifdef CCODE */



#ifdef FCODE

/* LOOP */

#define CCTK_LOOP2_NORMAL_DECLARE(name) \
   CCTK_LOOP2STR_NORMAL_DECLARE(name) \
   && integer :: name/**/0_dummy_imin, name/**/0_dummy_imax \

#define CCTK_LOOP2_NORMAL_OMP_PRIVATE(name) \
   CCTK_LOOP2STR_NORMAL_OMP_PRIVATE(name) \

#define CCTK_LOOP2_NORMAL(name, \
                          i,j, \
                          ni,nj, \
                          cctki0_idir,cctki0_jdir, \
                          cctki0_imin,cctki0_jmin, \
                          cctki0_imax,cctki0_jmax, \
                          cctki0_iash,cctki0_jash) \
   CCTK_LOOP2STR_NORMAL(name, \
                        i,j, \
                        ni,nj, \
                        cctki0_idir,cctki0_jdir, \
                        cctki0_imin,cctki0_jmin, \
                        cctki0_imax,cctki0_jmax, \
                        cctki0_iash,cctki0_jash, \
                        name/**/0_dummy_imin,name/**/0_dummy_imax, 1) \

#define CCTK_ENDLOOP2_NORMAL(name) \
   CCTK_ENDLOOP2STR_NORMAL(name) \

#define CCTK_LOOP2STR_NORMAL_DECLARE(name) \
   && integer :: name/**/0_idir,name/**/0_jdir \
   && integer :: name/**/0_imin,name/**/0_jmin \
   && integer :: name/**/0_imax,name/**/0_jmax \
   && integer :: name/**/0_iash,name/**/0_jash \
   && integer :: name/**/0_istr \

#define CCTK_LOOP2STR_NORMAL_OMP_PRIVATE(name) \
   && !$omp private (i,j) \
   && !$omp private (ni,nj) \

#define CCTK_LOOP2STR_NORMAL(name, \
                             i,j, \
                             ni,nj, \
                             cctki0_idir,cctki0_jdir, \
                             cctki0_imin,cctki0_jmin, \
                             cctki0_imax,cctki0_jmax, \
                             cctki0_iash,cctki0_jash, \
                             imin,imax, cctki0_istr) \
   && name/**/0_idir = cctki0_idir \
   && name/**/0_jdir = cctki0_jdir \
   && name/**/0_imin = cctki0_imin \
   && name/**/0_jmin = cctki0_jmin \
   && name/**/0_imax = cctki0_imax \
   && name/**/0_jmax = cctki0_jmax \
   && name/**/0_iash = cctki0_iash \
   && name/**/0_jash = cctki0_jash \
   && name/**/0_istr = cctki0_istr \
   && imin = name/**/0_imin \
   && imax = name/**/0_imax \
   && !$omp do collapse(1) \
   && do j = name/**/0_jmin, name/**/0_jmax \
   && do i = name/**/0_imin - modulo((imin+name/**/0_iash*(j)), name/**/0_istr), name/**/0_imax, name/**/0_istr \
   &&    ni = 0 \
   &&    nj = 0 \
   &&    if (name/**/0_idir < 0) ni = i \
   &&    if (name/**/0_jdir < 0) nj = j \
   &&    if (name/**/0_idir > 0) ni = name/**/0_imax+1-i \
   &&    if (name/**/0_jdir > 0) nj = name/**/0_jmax+1-j \

#define CCTK_ENDLOOP2STR_NORMAL(name) \
   && end do \
   && end do \



#define CCTK_LOOP2_DECLARE(name) \
   CCTK_LOOP2STR_DECLARE(name) \
   && integer :: name/**/1_dummy_imin, name/**/1_dummy_imax \

#define CCTK_LOOP2_OMP_PRIVATE(name) \
   CCTK_LOOP2STR_OMP_PRIVATE(name) \

#define CCTK_LOOP2(name, \
                   i,j, \
                   cctki0_imin,cctki0_jmin, \
                   cctki0_imax,cctki0_jmax, \
                   cctki0_iash,cctki0_jash) \
   CCTK_LOOP2STR(name, \
                 i,j, \
                 cctki0_imin,cctki0_jmin, \
                 cctki0_imax,cctki0_jmax, \
                 cctki0_iash,cctki0_jash, \
                 name/**/1_dummy_imin,name/**/1_dummy_imax, 1) \

#define CCTK_ENDLOOP2(name) \
   CCTK_ENDLOOP2STR(name) \

#define CCTK_LOOP2STR_DECLARE(name) \
   CCTK_LOOP2STR_NORMAL_DECLARE(name) \
   && integer :: name/**/1_ni,name/**/1_nj \

#define CCTK_LOOP2STR_OMP_PRIVATE(name) \
   CCTK_LOOP2STR_NORMAL_OMP_PRIVATE(name) \

#define CCTK_LOOP2STR(name, \
                      i,j, \
                      cctki1_imin,cctki1_jmin, \
                      cctki1_imax,cctki1_jmax, \
                      cctki1_iash,cctki1_jash, \
                      imin,imax, cctki1_istr) \
   CCTK_LOOP2STR_NORMAL(name, \
                        i,j, \
                        name/**/1_ni,name/**/1_nj, \
                        0,0, \
                        cctki1_imin,cctki1_jmin, \
                        cctki1_imax,cctki1_jmax, \
                        cctki1_iash,cctki1_jash, \
                        imin,imax, cctki1_istr) \

#define CCTK_ENDLOOP2STR(name) \
   CCTK_ENDLOOP2STR_NORMAL(name) \



/* LOOP_INTERIOR */

#define CCTK_LOOP2_INTERIOR_DECLARE(name) \
   CCTK_LOOP2STR_INTERIOR_DECLARE(name) \
   && integer :: name/**/2_dummy_imin, name/**/2_dummy_imax \

#define CCTK_LOOP2_INTERIOR_OMP_PRIVATE(name) \
   CCTK_LOOP2STR_INTERIOR_OMP_PRIVATE(name) \

#define CCTK_LOOP2_INTERIOR(name, \
                            i,j, \
                            cctki2_iblo,cctki2_jblo, \
                            cctki2_ibhi,cctki2_jbhi) \
   CCTK_LOOP2STR_INTERIOR(name, \
                          i,j, \
                          cctki2_iblo,cctki2_jblo, \
                          cctki2_ibhi,cctki2_jbhi, \
                          name/**/2_dummy_imin,name/**/2_dummy_imax, 1) \

#define CCTK_ENDLOOP2_INTERIOR(name) \
   CCTK_ENDLOOP2STR_INTERIOR(name) \

#define CCTK_LOOP2STR_INTERIOR_DECLARE(name) \
   CCTK_LOOP2STR_DECLARE(name/**/_interior) \

#define CCTK_LOOP2STR_INTERIOR_OMP_PRIVATE(name) \
   CCTK_LOOP2STR_OMP_PRIVATE(name/**/_interior) \

#define CCTK_LOOP2STR_INTERIOR(name, \
                               i,j, \
                               cctki2_iblo,cctki2_jblo, \
                               cctki2_ibhi,cctki2_jbhi, \
                               imin,imax, cctki2_istr) \
   CCTK_LOOP2STR(name/**/_interior, \
                 i,j, \
                 (cctki2_iblo)+1, \
                 (cctki2_jblo)+1, \
                 cctk_lsh(1)-(cctki2_ibhi), \
                 cctk_lsh(2)-(cctki2_jbhi), \
                 cctk_ash(1),cctk_ash(2), \
                 imin,imax, cctki2_istr) \

#define CCTK_ENDLOOP2STR_INTERIOR(name) \
   CCTK_ENDLOOP2STR(name/**/_interior) \



/* LOOP_BOUNDARIES */

#define CCTK_LOOP2_BOUNDARIES_DECLARE(name) \
   CCTK_LOOP2STR_BOUNDARIES_DECLARE(name) \
   && integer :: name/**/2_dummy_imin, name/**/2_dummy_imax \

#define CCTK_LOOP2_BOUNDARIES_OMP_PRIVATE(name) \
   CCTK_LOOP2STR_BOUNDARIES_OMP_PRIVATE(name) \

#define CCTK_LOOP2_BOUNDARIES(name, \
                              i,j, \
                              ni,nj, \
                              cctki2_iblo,cctki2_jblo, \
                              cctki2_ibhi,cctki2_jbhi, \
                              cctki2_ibboxlo,cctki2_jbboxlo, \
                              cctki2_ibboxhi,cctki2_jbboxhi) \
   CCTK_LOOP2STR_BOUNDARIES(name, \
                            i,j, \
                            ni,nj, \
                            cctki2_iblo,cctki2_jblo, \
                            cctki2_ibhi,cctki2_jbhi, \
                            cctki2_ibboxlo,cctki2_jbboxlo, \
                            cctki2_ibboxhi,cctki2_jbboxhi, \
                            name/**/2_dummy_imin,name/**/2_dummy_imax, 1) \

#define CCTK_ENDLOOP2_BOUNDARIES(name) \
   CCTK_ENDLOOP2STR_BOUNDARIES(name) \

#define CCTK_LOOP2STR_BOUNDARIES_DECLARE(name) \
   CCTK_LOOP2STR_NORMAL_DECLARE(name/**/_boundaries) \
   && integer :: name/**/2_blo(2), name/**/2_bhi(2) \
   && integer :: name/**/2_bboxlo(2), name/**/2_bboxhi(2) \
   && integer :: name/**/2_istr \
   && integer :: name/**/2_idir \
   && integer :: name/**/2_jdir \
   && logical :: name/**/2_any_bbox \
   && integer :: name/**/2_bmin(2), name/**/2_bmax(2) \

#define CCTK_LOOP2STR_BOUNDARIES_OMP_PRIVATE(name) \
   CCTK_LOOP2STR_NORMAL_OMP_PRIVATE(name/**/_boundaries) \
   && !$omp private (name/**/2_bmin, name/**/2_bmax) \

#define CCTK_LOOP2STR_BOUNDARIES(name, \
                                 i,j, \
                                 ni,nj, \
                                 cctki2_iblo,cctki2_jblo, \
                                 cctki2_ibhi,cctki2_jbhi, \
                                 cctki2_ibboxlo,cctki2_jbboxlo, \
                                 cctki2_ibboxhi,cctki2_jbboxhi, \
                                 imin,imax, cctki2_istr) \
   && name/**/2_blo = (/ cctki2_iblo,cctki2_jblo /) \
   && name/**/2_bhi = (/ cctki2_ibhi,cctki2_jbhi /) \
   && name/**/2_bboxlo = (/ cctki2_ibboxlo,cctki2_jbboxlo /) \
   && name/**/2_bboxhi = (/ cctki2_ibboxhi,cctki2_jbboxhi /) \
   && name/**/2_istr = (cctki2_istr) \
   && do name/**/2_jdir=-1, +1 \
   && do name/**/2_idir=-1, +1 \
   &&     name/**/2_any_bbox = .false. \
   &&     if (name/**/2_idir<0) name/**/2_any_bbox = name/**/2_any_bbox .or. name/**/2_bboxlo(1) /= 0 \
   &&     if (name/**/2_jdir<0) name/**/2_any_bbox = name/**/2_any_bbox .or. name/**/2_bboxlo(2) /= 0 \
   &&     if (name/**/2_idir>0) name/**/2_any_bbox = name/**/2_any_bbox .or. name/**/2_bboxhi(1) /= 0 \
   &&     if (name/**/2_jdir>0) name/**/2_any_bbox = name/**/2_any_bbox .or. name/**/2_bboxhi(2) /= 0 \
   &&    if (name/**/2_any_bbox) then \
   &&       name/**/2_bmin(1) = name/**/2_blo(1)+1 \
   &&       name/**/2_bmin(2) = name/**/2_blo(2)+1 \
   &&       if (name/**/2_idir<0) name/**/2_bmin(1) = 1 \
   &&       if (name/**/2_jdir<0) name/**/2_bmin(2) = 1 \
   &&       if (name/**/2_idir>0) name/**/2_bmin(1) = cctk_lsh(1) - name/**/2_bhi(1) \
   &&       if (name/**/2_jdir>0) name/**/2_bmin(2) = cctk_lsh(2) - name/**/2_bhi(2) \
   &&       name/**/2_bmax(1) = cctk_lsh(1) - name/**/2_bhi(1) \
   &&       name/**/2_bmax(2) = cctk_lsh(2) - name/**/2_bhi(2) \
   &&       if (name/**/2_idir<0) name/**/2_bmax(1) = name/**/2_blo(1) \
   &&       if (name/**/2_jdir<0) name/**/2_bmax(2) = name/**/2_blo(2) \
   &&       if (name/**/2_idir>0) name/**/2_bmax(1) = cctk_lsh(1) \
   &&       if (name/**/2_jdir>0) name/**/2_bmax(2) = cctk_lsh(2) \
   &&       CCTK_LOOP2STR_NORMAL(name/**/_boundaries, \
                                 i,j, \
                                 ni,nj, \
                                 name/**/2_idir,name/**/2_jdir, \
                                 name/**/2_bmin(1),name/**/2_bmin(2), \
                                 name/**/2_bmax(1),name/**/2_bmax(2), \
                                 cctk_ash(1), \
                                 cctk_ash(2), \
                                 imin,imax, name/**/2_istr) \

#define CCTK_ENDLOOP2STR_BOUNDARIES(name) \
            CCTK_ENDLOOP2STR_NORMAL(name/**/_boundaries) \
   &&    end if /* bbox */ \
   && end do /* dir */ \
   && end do /* dir */ \



/* LOOP_INTBOUNDARIES */

#define CCTK_LOOP2_INTBOUNDARIES_DECLARE(name) \
   CCTK_LOOP2STR_INTBOUNDARIES_DECLARE(name) \
   && integer :: name/**/2_dummy_imin, name/**/2_dummy_imax \

#define CCTK_LOOP2_INTBOUNDARIES_OMP_PRIVATE(name) \
   CCTK_LOOP2STR_INTBOUNDARIES_OMP_PRIVATE(name) \

#define CCTK_LOOP2_INTBOUNDARIES(name, \
                                 i,j, \
                                 ni,nj, \
                                 cctki2_iblo,cctki2_jblo, \
                                 cctki2_ibhi,cctki2_jbhi, \
                                 cctki2_ibboxlo,cctki2_jbboxlo, \
                                 cctki2_ibboxhi,cctki2_jbboxhi) \
   CCTK_LOOP2STR_INTBOUNDARIES(name, \
                               i,j, \
                               ni,nj, \
                               cctki2_iblo,cctki2_jblo, \
                               cctki2_ibhi,cctki2_jbhi, \
                               cctki2_ibboxlo,cctki2_jbboxlo, \
                               cctki2_ibboxhi,cctki2_jbboxhi, \
                               name/**/2_dummy_imin,name/**/2_dummy_max, 1) \

#define CCTK_ENDLOOP2_INTBOUNDARIES(name) \
   CCTK_ENDLOOP2STR_INTBOUNDARIES(name) \

#define CCTK_LOOP2STR_INTBOUNDARIES_DECLARE(name) \
   CCTK_LOOP2STR_NORMAL_DECLARE(name/**/_intboundaries) \
   && integer :: name/**/2_blo(2), name/**/2_bhi(2) \
   && integer :: name/**/2_bboxlo(2), name/**/2_bboxhi(2) \
   && integer :: name/**/2_istr \
   && integer :: name/**/2_idir \
   && integer :: name/**/2_jdir \
   && logical :: name/**/2_any_bbox, name/**/2_all_bbox \
   && integer :: name/**/2_bmin(2), name/**/2_bmax(2) \

#define CCTK_LOOP2STR_INTBOUNDARIES_OMP_PRIVATE(name) \
   CCTK_LOOP2STR_NORMAL_OMP_PRIVATE(name/**/_intboundaries) \
   && !$omp private (name/**/2_any_bbox, name/**/2_all_bbox) \
   && !$omp private (name/**/2_bmin, name/**/2_bmax \

#define CCTK_LOOP2STR_INTBOUNDARIES(name, \
                                    i,j, \
                                    ni,nj, \
                                    cctki2_iblo,cctki2_jblo, \
                                    cctki2_ibhi,cctki2_jbhi, \
                                    cctki2_ibboxlo,cctki2_jbboxlo, \
                                    cctki2_ibboxhi,cctki2_jbboxhi, \
                                    imin,imax, cctki2_istr) \
   && name/**/2_blo = (/ cctki2_iblo,cctki2_jblo /) \
   && name/**/2_bhi = (/ cctki2_ibhi,cctki2_jbhi /) \
   && name/**/2_bboxlo = (/ cctki2_ibboxlo,cctki2_jbboxlo /) \
   && name/**/2_bboxhi = (/ cctki2_ibboxhi,cctki2_jbboxhi /) \
   && name/**/2_istr = (cctki2_istr) \
   && do name/**/2_jdir=-1, +1 \
   && do name/**/2_idir=-1, +1 \
   &&     name/**/2_any_bbox = .false. \
   &&     if (name/**/2_idir<0) name/**/2_any_bbox = name/**/2_any_bbox .or. name/**/2_bboxlo(1) /= 0 \
   &&     if (name/**/2_jdir<0) name/**/2_any_bbox = name/**/2_any_bbox .or. name/**/2_bboxlo(2) /= 0 \
   &&     if (name/**/2_idir>0) name/**/2_any_bbox = name/**/2_any_bbox .or. name/**/2_bboxhi(1) /= 0 \
   &&     if (name/**/2_jdir>0) name/**/2_any_bbox = name/**/2_any_bbox .or. name/**/2_bboxhi(2) /= 0 \
   &&     name/**/2_all_bbox = .true. \
   &&     if (name/**/2_idir<0) name/**/2_all_bbox = name/**/2_all_bbox .and. name/**/2_bboxlo(1) /= 0 \
   &&     if (name/**/2_jdir<0) name/**/2_all_bbox = name/**/2_all_bbox .and. name/**/2_bboxlo(2) /= 0 \
   &&     if (name/**/2_idir>0) name/**/2_all_bbox = name/**/2_all_bbox .and. name/**/2_bboxhi(1) /= 0 \
   &&     if (name/**/2_jdir>0) name/**/2_all_bbox = name/**/2_all_bbox .and. name/**/2_bboxhi(2) /= 0 \
   &&    if (name/**/2_all_bbox .and. name/**/2_any_bbox) then \
   &&       name/**/2_bmin(1) = name/**/2_blo(1)+1 \
   &&       name/**/2_bmin(2) = name/**/2_blo(2)+1 \
   &&       if (name/**/2_idir<0) name/**/2_bmin(1) = 1 \
   &&       if (name/**/2_jdir<0) name/**/2_bmin(2) = 1 \
   &&       if (name/**/2_idir>0) name/**/2_bmin(1) = cctk_lsh(1) - name/**/2_bhi(1) \
   &&       if (name/**/2_jdir>0) name/**/2_bmin(2) = cctk_lsh(2) - name/**/2_bhi(2) \
   &&       name/**/2_bmax(1) = cctk_lsh(1) - name/**/2_bhi(1) \
   &&       name/**/2_bmax(2) = cctk_lsh(2) - name/**/2_bhi(2) \
   &&       if (name/**/2_idir<0) name/**/2_bmax(1) = name/**/2_blo(1) \
   &&       if (name/**/2_jdir<0) name/**/2_bmax(2) = name/**/2_blo(2) \
   &&       if (name/**/2_idir>0) name/**/2_bmax(1) = cctk_lsh(1) \
   &&       if (name/**/2_jdir>0) name/**/2_bmax(2) = cctk_lsh(2) \
   &&       CCTK_LOOP2STR_NORMAL(name/**/_intboundaries, \
                                 i,j, \
                                 ni,nj, \
                                 name/**/2_idir,name/**/2_jdir, \
                                 name/**/2_bmin(1),name/**/2_bmin(2), \
                                 name/**/2_bmax(1),name/**/2_bmax(2), \
                                 cctk_ash(1), \
                                 cctk_ash(2), \
                                 imin,imax, name/**/2_istr) \

#define CCTK_ENDLOOP2STR_INTBOUNDARIES(name) \
            CCTK_ENDLOOP2STR_NORMAL(name/**/_intboundaries) \
   &&    end if /* bbox */ \
   && end do /* dir */ \
   && end do /* dir */ \



/* LOOP_ALL */

#define CCTK_LOOP2_ALL_DECLARE(name) \
   CCTK_LOOP2STR_ALL_DECLARE(name) \
   && integer :: name/**/3_dummy_imin, name/**/3_dummy_imax \

#define CCTK_LOOP2_ALL_OMP_PRIVATE(name) \
   CCTK_LOOP2STR_ALL_OMP_PRIVATE(name) \

#define CCTK_LOOP2_ALL(name, \
                       i,j) \
   CCTK_LOOP2STR_ALL(name, \
                     i,j, \
                     name/**/3_dummy_imin,name/**/3_dummy_imax, 1) \

#define CCTK_ENDLOOP2_ALL(name) \
   CCTK_ENDLOOP2STR_ALL(name) \

#define CCTK_LOOP2STR_ALL_DECLARE(name) \
   CCTK_LOOP2STR_DECLARE(name/**/_all) \

#define CCTK_LOOP2STR_ALL_OMP_PRIVATE(name) \
   CCTK_LOOP2STR_OMP_PRIVATE(name/**/_all) \

#define CCTK_LOOP2STR_ALL(name, \
                          i,j, \
                          imin,imax, cctki3_istr) \
   CCTK_LOOP2STR(name/**/_all, \
                 i,j, \
                 1,1, \
                 cctk_lsh(1),cctk_lsh(2), \
                 cctk_ash(1),cctk_ash(2), \
                 imin,imax, cctki3_istr) \

#define CCTK_ENDLOOP2STR_ALL(name) \
   CCTK_ENDLOOP2STR(name/**/_all) \



/* LOOP_INT */

#define CCTK_LOOP2_INT_DECLARE(name) \
   CCTK_LOOP2STR_INT_DECLARE(name) \
   && integer :: name/**/3_dummy_imin, name/**/3_dummy_imax \

#define CCTK_LOOP2_INT_OMP_PRIVATE(name) \
   CCTK_LOOP2STR_INT_OMP_PRIVATE(name) \

#define CCTK_LOOP2_INT(name, \
                        i,j) \
   CCTK_LOOP2STR_INT(name, \
                      i,j, \
                      name/**/3_dummy_imin,name/**/3_dummy_imax, 1) \

#define CCTK_ENDLOOP2_INT(name) \
   CCTK_ENDLOOP2STR_INT(name) \

#define CCTK_LOOP2STR_INT_DECLARE(name) \
   CCTK_LOOP2STR_INTERIOR_DECLARE(name/**/_int) \
   && CCTK_INT, parameter :: name/**/3_isize = 4 \
   && CCTK_INT :: name/**/3_bndsize    (4) \
   && CCTK_INT :: name/**/3_is_ghostbnd(4) \
   && CCTK_INT :: name/**/3_is_symbnd  (4) \
   && CCTK_INT :: name/**/3_is_physbnd (4) \
   && CCTK_INT :: name/**/3_ierr \

#define CCTK_LOOP2STR_INT_OMP_PRIVATE(name) \
   CCTK_LOOP2STR_INTERIOR_OMP_PRIVATE(name/**/_int) \

#define CCTK_LOOP2STR_INT(name, \
                          i,j, \
                          imin,imax, cctki3_istr) \
   && !$omp single \
   && name/**/3_ierr = GetBoundarySizesAndTypes \
         (cctkGH, name/**/3_isize, name/**/3_bndsize, name/**/3_is_ghostbnd, name/**/3_is_symbnd, name/**/3_is_physbnd) \
   && !$omp end single copyprivate(name/**/3_bndsize) \
   && CCTK_LOOP2STR_INTERIOR(name/**/_int, \
                             i,j, \
                             int(name/**/3_bndsize(1+1)),int(name/**/3_bndsize(3+1)), \
                             int(name/**/3_bndsize(2)),int(name/**/3_bndsize(4)), \
                             imin,imax, (cctki3_istr)) \

#define CCTK_ENDLOOP2STR_INT(name) \
      CCTK_ENDLOOP2STR_INTERIOR(name/**/_int) \



/* LOOP_BND */

#define CCTK_LOOP2_BND_DECLARE(name) \
   CCTK_LOOP2STR_BND_DECLARE(name) \
   && integer :: name/**/3_dummy_imin, name/**/3_dummy_imax \

#define CCTK_LOOP2_BND_OMP_PRIVATE(name) \
   CCTK_LOOP2STR_BND_OMP_PRIVATE(name) \

#define CCTK_LOOP2_BND(name, \
                       i,j, \
                       ni,nj) \
   CCTK_LOOP2STR_BND(name, \
                     i,j, \
                     ni,nj, \
                     name/**/3_dummy_imin,name/**/3_dummy_imax, 1) \

#define CCTK_ENDLOOP2_BND(name) \
   CCTK_ENDLOOP2STR_BND(name) \

#define CCTK_LOOP2STR_BND_DECLARE(name) \
   CCTK_LOOP2STR_BOUNDARIES_DECLARE(name/**/_bnd) \
   && CCTK_INT, parameter :: name/**/3_isize = 4 \
   && CCTK_INT :: name/**/3_bndsize    (4) \
   && CCTK_INT :: name/**/3_is_ghostbnd(4) \
   && CCTK_INT :: name/**/3_is_symbnd  (4) \
   && CCTK_INT :: name/**/3_is_physbnd (4) \
   && CCTK_INT :: name/**/3_ierr \

#define CCTK_LOOP2STR_BND_OMP_PRIVATE(name) \
   CCTK_LOOP2STR_BOUNDARIES_OMP_PRIVATE(name/**/_bnd) \

#define CCTK_LOOP2STR_BND(name, \
                          i,j, \
                          ni,nj, \
                          imin,imax, cctki3_istr) \
   && !$omp single \
   && name/**/3_ierr = GetBoundarySizesAndTypes \
         (cctkGH, name/**/3_isize, name/**/3_bndsize, name/**/3_is_ghostbnd, name/**/3_is_symbnd, name/**/3_is_physbnd) \
   && !$omp end single copyprivate(name/**/3_bndsize, name/**/3_is_physbnd) \
   && CCTK_LOOP2STR_BOUNDARIES(name/**/_bnd, \
                               i,j, \
                               ni,nj, \
                               int(name/**/3_bndsize(1))+1,int(name/**/3_bndsize(3))+1, \
                               int(name/**/3_bndsize(2)),int(name/**/3_bndsize(4)), \
                               int(name/**/3_is_physbnd(1)),int(name/**/3_is_physbnd(3)), \
                               int(name/**/3_is_physbnd(2)),int(name/**/3_is_physbnd(4)), \
                               imin,imax, (cctki3_istr)) \

#define CCTK_ENDLOOP2STR_BND(name) \
      CCTK_ENDLOOP2STR_BOUNDARIES(name/**/_bnd) \






/* LOOP_INTBND */

#define CCTK_LOOP2_INTBND_DECLARE(name) \
   CCTK_LOOP2STR_INTBND_DECLARE(name) \
   && integer :: name/**/3_dummy_imin, name/**/3_dummy_imax \

#define CCTK_LOOP2_INTBND_OMP_PRIVATE(name) \
   CCTK_LOOP2STR_INTBND_OMP_PRIVATE(name) \

#define CCTK_LOOP2_INTBND(name, \
                          i,j, \
                          ni,nj) \
   CCTK_LOOP2STR_INTBND(name, \
                        i,j, \
                        ni,nj, \
                        name/**/3_dummy_imin,name/**/3_dummy_imax, 1) \

#define CCTK_ENDLOOP2_INTBND(name) \
   CCTK_ENDLOOP2STR_INTBND(name) \

#define CCTK_LOOP2STR_INTBND_DECLARE(name) \
   CCTK_LOOP2STR_INTBOUNDARIES_DECLARE(name/**/_bnd) \
   && CCTK_INT, parameter :: name/**/3_isize = 4 \
   && CCTK_INT :: name/**/3_bndsize    (4) \
   && CCTK_INT :: name/**/3_is_ghostbnd(4) \
   && CCTK_INT :: name/**/3_is_symbnd  (4) \
   && CCTK_INT :: name/**/3_is_physbnd (4) \
   && CCTK_INT :: name/**/3_ierr \

#define CCTK_LOOP2STR_INTBND_OMP_PRIVATE(name) \
   CCTK_LOOP2STR_INTBOUNDARIES_OMP_PRIVATE(name/**/_bnd) \

#define CCTK_LOOP2STR_INTBND(name, \
                             i,j, \
                             ni,nj, \
                             imin,imax, cctki3_istr) \
   && !$omp single \
   && name/**/3_ierr = GetBoundarySizesAndTypes \
         (cctkGH, name/**/3_isize, name/**/3_bndsize, name/**/3_is_ghostbnd, name/**/3_is_symbnd, name/**/3_is_physbnd) \
   && !$omp end single copyprivate(name/**/3_bndsize, name/**/3_is_physbnd) \
   && CCTK_LOOP2STR_INTBOUNDARIES(name/**/_bnd, \
                                  i,j, \
                                  ni,nj, \
                                  int(name/**/3_bndsize(1+1)),int(name/**/3_bndsize(3+1)), \
                                  int(name/**/3_bndsize(2)),int(name/**/3_bndsize(4)), \
                                  int(name/**/3_is_physbnd(1)),int(name/**/3_is_physbnd(3)), \
                                  int(name/**/3_is_physbnd(2)),int(name/**/3_is_physbnd(4)), \
                                  imin,imax, (cctki3_istr)) \

#define CCTK_ENDLOOP2STR_INTBND(name) \
      CCTK_ENDLOOP2STR_INTBOUNDARIES(name/**/_bnd) \

#endif /* #ifdef FCODE */



/* 3D */

#ifdef CCODE

/* LOOP */

#define CCTK_LOOP3_NORMAL(name, \
                          i,j,k, \
                          ni,nj,nk, \
                          cctki0_idir_,cctki0_jdir_,cctki0_kdir_, \
                          cctki0_imin_,cctki0_jmin_,cctki0_kmin_, \
                          cctki0_imax_,cctki0_jmax_,cctki0_kmax_, \
                          cctki0_iash_,cctki0_jash_,cctki0_kash_) \
  CCTK_LOOP3STR_NORMAL(name, \
                       i,j,k, \
                       ni,nj,nk, \
                       (cctki0_idir_),(cctki0_jdir_),(cctki0_kdir_), \
                       (cctki0_imin_),(cctki0_jmin_),(cctki0_kmin_), \
                       (cctki0_imax_),(cctki0_jmax_),(cctki0_kmax_), \
                       (cctki0_iash_),(cctki0_jash_),(cctki0_kash_), \
                       cctki0_dummy_imin,cctki0_dummy_imax, 1) \

#define CCTK_ENDLOOP3_NORMAL(name) \
  CCTK_ENDLOOP3STR_NORMAL(name) \

#define CCTK_LOOP3STR_NORMAL(name, \
                             i,j,k, \
                             ni,nj,nk, \
                             cctki0_idir_,cctki0_jdir_,cctki0_kdir_, \
                             cctki0_imin_,cctki0_jmin_,cctki0_kmin_, \
                             cctki0_imax_,cctki0_jmax_,cctki0_kmax_, \
                             cctki0_iash_,cctki0_jash_,cctki0_kash_, \
                             imin,imax, cctki0_istr_) \
  do { \
    typedef int cctki0_loop3_normal_##name; \
    const int cctki0_idir = (cctki0_idir_); \
    const int cctki0_jdir = (cctki0_jdir_); \
    const int cctki0_kdir = (cctki0_kdir_); \
    const int cctki0_imin = (cctki0_imin_); \
    const int cctki0_jmin = (cctki0_jmin_); \
    const int cctki0_kmin = (cctki0_kmin_); \
    const int cctki0_imax = (cctki0_imax_); \
    const int cctki0_jmax = (cctki0_jmax_); \
    const int cctki0_kmax = (cctki0_kmax_); \
    const int cctki0_iash CCTK_ATTRIBUTE_UNUSED = (cctki0_iash_); \
    const int cctki0_jash CCTK_ATTRIBUTE_UNUSED = (cctki0_jash_); \
    const int cctki0_kash CCTK_ATTRIBUTE_UNUSED = (cctki0_kash_); \
    const int cctki0_istr = (cctki0_istr_); \
    assert(cctki0_istr>0 && (cctki0_istr & (cctki0_istr-1)) == 0); \
    const int imin CCTK_ATTRIBUTE_UNUSED = cctki0_imin; \
    const int imax CCTK_ATTRIBUTE_UNUSED = cctki0_imax; \
    CCTK_PRAGMA_OMP_FOR_COLLAPSE_2 \
    for (int k=cctki0_kmin; k<cctki0_kmax; ++k) { \
    for (int j=cctki0_jmin; j<cctki0_jmax; ++j) { \
     \
    const int cctki0_ioff = (cctki0_imin+cctki0_iash*(j+cctki0_jash*(k))) & (cctki0_istr-1); \
     \
    for (int i=cctki0_imin-cctki0_ioff; i<cctki0_imax; i+=cctki0_istr) { \
      const int ni CCTK_ATTRIBUTE_UNUSED = cctki0_idir<0 ? i+1 : cctki0_idir==0 ? 0 : cctki0_imax-i; \
      const int nj CCTK_ATTRIBUTE_UNUSED = cctki0_jdir<0 ? j+1 : cctki0_jdir==0 ? 0 : cctki0_jmax-j; \
      const int nk CCTK_ATTRIBUTE_UNUSED = cctki0_kdir<0 ? k+1 : cctki0_kdir==0 ? 0 : cctki0_kmax-k; \
      { \

#define CCTK_ENDLOOP3STR_NORMAL(name) \
      } \
    } \
    } \
    } \
    typedef cctki0_loop3_normal_##name cctki0_ensure_proper_nesting CCTK_ATTRIBUTE_UNUSED; \
  } while (0) \



#define CCTK_LOOP3(name, \
                   i,j,k, \
                   cctki1_imin_,cctki1_jmin_,cctki1_kmin_, \
                   cctki1_imax_,cctki1_jmax_,cctki1_kmax_, \
                   cctki1_iash_,cctki1_jash_,cctki1_kash_) \
  CCTK_LOOP3STR(name, \
                i,j,k, \
                (cctki1_imin_),(cctki1_jmin_),(cctki1_kmin_), \
                (cctki1_imax_),(cctki1_jmax_),(cctki1_kmax_), \
                (cctki1_iash_),(cctki1_jash_),(cctki1_kash_), \
                cctki1_dummy_imin,cctki1_dummy_imax, 1) \

#define CCTK_ENDLOOP3(name) \
  CCTK_ENDLOOP3STR(name) \

#define CCTK_LOOP3STR(name, \
                      i,j,k, \
                      cctki1_imin_,cctki1_jmin_,cctki1_kmin_, \
                      cctki1_imax_,cctki1_jmax_,cctki1_kmax_, \
                      cctki1_iash_,cctki1_jash_,cctki1_kash_, \
                      imin,imax, cctki1_istr_) \
  CCTK_LOOP3STR_NORMAL(name, \
                       i,j,k, \
                       cctki1_ni,cctki1_nj,cctki1_nk, \
                       0,0,0, \
                       (cctki1_imin_),(cctki1_jmin_),(cctki1_kmin_), \
                       (cctki1_imax_),(cctki1_jmax_),(cctki1_kmax_), \
                       (cctki1_iash_),(cctki1_jash_),(cctki1_kash_), \
                       imin,imax, (cctki1_istr_)) \

#define CCTK_ENDLOOP3STR(name) \
  CCTK_ENDLOOP3STR_NORMAL(name) \



/* LOOP_INTERIOR */

#define CCTK_LOOP3_INTERIOR(name, cctki2_cctkGH_, \
                            i,j,k, \
                            cctki2_iblo_,cctki2_jblo_,cctki2_kblo_, \
                            cctki2_ibhi_,cctki2_jbhi_,cctki2_kbhi_) \
  CCTK_LOOP3STR_INTERIOR(name, (cctki2_cctkGH_), \
                         i,j,k, \
                         (cctki2_iblo_),(cctki2_jblo_),(cctki2_kblo_), \
                         (cctki2_ibhi_),(cctki2_jbhi_),(cctki2_kbhi_), \
                         cctki2_dummy_imin,cctki2_dummy_imax, 1) \

#define CCTK_ENDLOOP3_INTERIOR(name) \
  CCTK_ENDLOOP3STR_INTERIOR(name) \

#define CCTK_LOOP3STR_INTERIOR(name, cctki2_cctkGH_, \
                               i,j,k, \
                               cctki2_iblo_,cctki2_jblo_,cctki2_kblo_, \
                               cctki2_ibhi_,cctki2_jbhi_,cctki2_kbhi_, \
                               imin,imax, cctki2_istr_) \
  do { \
    typedef int cctki2_loop3_interior_##name; \
    cGH const *restrict const cctki2_cctkGH = (cctki2_cctkGH_); \
    if (cctki2_cctkGH->cctk_dim != 3) { \
      _Pragma("omp critical") \
      CCTK_WARN(CCTK_WARN_ABORT, \
                "The macro CCTK_LOOP3_INTERIOR can only be used in 3 dimensions"); \
    } \
    CCTK_LOOP3STR(name##_interior, \
                  i,j,k, \
                  (cctki2_iblo_),(cctki2_jblo_),(cctki2_kblo_), \
                  cctki2_cctkGH->cctk_lsh[0]-(cctki2_ibhi_), \
                  cctki2_cctkGH->cctk_lsh[1]-(cctki2_jbhi_), \
                  cctki2_cctkGH->cctk_lsh[2]-(cctki2_kbhi_), \
                  cctki2_cctkGH->cctk_ash[0], \
                  cctki2_cctkGH->cctk_ash[1], \
                  cctki2_cctkGH->cctk_ash[2], \
                  imin,imax, (cctki2_istr_)) { \

#define CCTK_ENDLOOP3STR_INTERIOR(name) \
    } CCTK_ENDLOOP3STR(name##_interior); \
    typedef cctki2_loop3_interior_##name cctki2_ensure_proper_nesting CCTK_ATTRIBUTE_UNUSED; \
  } while(0) \



/* LOOP_BOUNDARIES */

#define CCTK_LOOP3_BOUNDARIES(name, cctki2_cctkGH_, \
                              i,j,k, \
                              ni,nj,nk, \
                              cctki2_iblo_,cctki2_jblo_,cctki2_kblo_, \
                              cctki2_ibhi_,cctki2_jbhi_,cctki2_kbhi_, \
                              cctki2_ibboxlo_,cctki2_jbboxlo_,cctki2_kbboxlo_, \
                              cctki2_ibboxhi_,cctki2_jbboxhi_,cctki2_kbboxhi_) \
  CCTK_LOOP3STR_BOUNDARIES(name, (cctki2_cctkGH_), \
                           i,j,k, \
                           ni,nj,nk, \
                           (cctki2_iblo_),(cctki2_jblo_),(cctki2_kblo_), \
                           (cctki2_ibhi_),(cctki2_jbhi_),(cctki2_kbhi_), \
                           (cctki2_ibboxlo_),(cctki2_jbboxlo_),(cctki2_kbboxlo_), \
                           (cctki2_ibboxhi_),(cctki2_jbboxhi_),(cctki2_kbboxhi_), \
                           cctki2_dummy_imin,cctki2_dummy_imax, 1) \

#define CCTK_ENDLOOP3_BOUNDARIES(name) \
  CCTK_ENDLOOP3STR_BOUNDARIES(name) \

#define CCTK_LOOP3STR_BOUNDARIES(name, cctki2_cctkGH_, \
                                 i,j,k, \
                                 ni,nj,nk, \
                                 cctki2_iblo_,cctki2_jblo_,cctki2_kblo_, \
                                 cctki2_ibhi_,cctki2_jbhi_,cctki2_kbhi_, \
                                 cctki2_ibboxlo_,cctki2_jbboxlo_,cctki2_kbboxlo_, \
                                 cctki2_ibboxhi_,cctki2_jbboxhi_,cctki2_kbboxhi_, \
                                 imin,imax, cctki2_istr_) \
  do { \
    typedef int cctki2_loop3_boundaries_##name; \
    cGH const *restrict const cctki2_cctkGH = (cctki2_cctkGH_); \
    if (cctki2_cctkGH->cctk_dim != 3) { \
      _Pragma("omp critical") \
      CCTK_WARN(CCTK_WARN_ABORT, \
                "The macro CCTK_LOOP3_BOUNDARIES can only be used in 3 dimensions"); \
    } \
    const int cctki2_blo[] = { (cctki2_iblo_),(cctki2_jblo_),(cctki2_kblo_) }; \
    const int cctki2_bhi[] = { (cctki2_ibhi_),(cctki2_jbhi_),(cctki2_kbhi_) }; \
    const int cctki2_bbox[] = { (cctki2_ibboxlo_), (cctki2_ibboxhi_),(cctki2_jbboxlo_), (cctki2_jbboxhi_),(cctki2_kbboxlo_), (cctki2_kbboxhi_) }; \
    const int cctki2_lsh[] = { cctki2_cctkGH->cctk_lsh[0],cctki2_cctkGH->cctk_lsh[1],cctki2_cctkGH->cctk_lsh[2] }; \
    const int cctki2_istr CCTK_ATTRIBUTE_UNUSED = (cctki2_istr_); \
    for (int cctki2_kdir=-1; cctki2_kdir<=+1; ++cctki2_kdir) { \
    for (int cctki2_jdir=-1; cctki2_jdir<=+1; ++cctki2_jdir) { \
    for (int cctki2_idir=-1; cctki2_idir<=+1; ++cctki2_idir) { \
      const int cctki2_any_bbox = \
        (cctki2_idir<0 ? cctki2_bbox[0] : 0) || (cctki2_idir>0 ? cctki2_bbox[1] : 0) || \
        (cctki2_jdir<0 ? cctki2_bbox[2] : 0) || (cctki2_jdir>0 ? cctki2_bbox[3] : 0) || \
        (cctki2_kdir<0 ? cctki2_bbox[4] : 0) || (cctki2_kdir>0 ? cctki2_bbox[5] : 0); \
      if (cctki2_any_bbox) { \
        const int cctki2_bmin[] = { \
          cctki2_idir<0 ? 0 : cctki2_idir==0 ? cctki2_blo[0] : cctki2_lsh[0] - cctki2_bhi[0], \
          cctki2_jdir<0 ? 0 : cctki2_jdir==0 ? cctki2_blo[1] : cctki2_lsh[1] - cctki2_bhi[1], \
          cctki2_kdir<0 ? 0 : cctki2_kdir==0 ? cctki2_blo[2] : cctki2_lsh[2] - cctki2_bhi[2], \
        }; \
        const int cctki2_bmax[] = { \
          cctki2_idir<0 ? cctki2_blo[0] : cctki2_idir==0 ? cctki2_lsh[0] - cctki2_bhi[0] : cctki2_lsh[0], \
          cctki2_jdir<0 ? cctki2_blo[1] : cctki2_jdir==0 ? cctki2_lsh[1] - cctki2_bhi[1] : cctki2_lsh[1], \
          cctki2_kdir<0 ? cctki2_blo[2] : cctki2_kdir==0 ? cctki2_lsh[2] - cctki2_bhi[2] : cctki2_lsh[2], \
        }; \
        CCTK_LOOP3STR_NORMAL(name##_boundaries, \
                             i,j,k, \
                             ni,nj,nk, \
                             cctki2_idir,cctki2_jdir,cctki2_kdir, \
                             cctki2_bmin[0],cctki2_bmin[1],cctki2_bmin[2], \
                             cctki2_bmax[0],cctki2_bmax[1],cctki2_bmax[2], \
                             cctki2_cctkGH->cctk_ash[0], \
                             cctki2_cctkGH->cctk_ash[1], \
                             cctki2_cctkGH->cctk_ash[2], \
                             imin,imax, cctki2_istr) { \

#define CCTK_ENDLOOP3STR_BOUNDARIES(name) \
        } CCTK_ENDLOOP3STR_NORMAL(name##_boundaries); \
      } /* if bbox */ \
    } /* for dir */ \
    } /* for dir */ \
    } /* for dir */ \
    typedef cctki2_loop3_boundaries_##name cctki2_ensure_proper_nesting CCTK_ATTRIBUTE_UNUSED; \
  } while (0) \



/* LOOP_INTBOUNDARIES */

#define CCTK_LOOP3_INTBOUNDARIES(name, cctki2_cctkGH_, \
                                 i,j,k, \
                                 ni,nj,nk, \
                                 cctki2_iblo_,cctki2_jblo_,cctki2_kblo_, \
                                 cctki2_ibhi_,cctki2_jbhi_,cctki2_kbhi_, \
                                 cctki2_ibboxlo_,cctki2_jbboxlo_,cctki2_kbboxlo_, \
                                 cctki2_ibboxhi_,cctki2_jbboxhi_,cctki2_kbboxhi_) \
  CCTK_LOOP3STR_INTBOUNDARIES(name, (cctki2_cctkGH_), \
                              i,j,k, \
                              ni,nj,nk, \
                              (cctki2_iblo_),(cctki2_jblo_),(cctki2_kblo_), \
                              (cctki2_ibhi_),(cctki2_jbhi_),(cctki2_kbhi_), \
                              (cctki2_ibboxlo_),(cctki2_jbboxlo_),(cctki2_kbboxlo_), \
                              (cctki2_ibboxhi_),(cctki2_jbboxhi_),(cctki2_kbboxhi_), \
                              cctki2_dummy_imin,cctki2_dummy_imax, 1) \

#define CCTK_ENDLOOP3_INTBOUNDARIES(name) \
  CCTK_ENDLOOP3STR_INTBOUNDARIES(name) \

#define CCTK_LOOP3STR_INTBOUNDARIES(name, cctki2_cctkGH_, \
                                    i,j,k, \
                                    ni,nj,nk, \
                                    cctki2_iblo_,cctki2_jblo_,cctki2_kblo_, \
                                    cctki2_ibhi_,cctki2_jbhi_,cctki2_kbhi_, \
                                    cctki2_ibboxlo_,cctki2_jbboxlo_,cctki2_kbboxlo_, \
                                    cctki2_ibboxhi_,cctki2_jbboxhi_,cctki2_kbboxhi_, \
                                    imin,imax, cctki2_istr_) \
  do { \
    typedef int cctki2_loop3_intboundaries_##name; \
    cGH const *restrict const cctki2_cctkGH = (cctki2_cctkGH_); \
    if (cctki2_cctkGH->cctk_dim != 3) { \
      _Pragma("omp critical") \
      CCTK_WARN(CCTK_WARN_ABORT, \
                "The macro CCTK_LOOP3_INTBOUNDARIES can only be used in 3 dimensions"); \
    } \
    const int cctki2_blo[] = { (cctki2_iblo_),(cctki2_jblo_),(cctki2_kblo_) }; \
    const int cctki2_bhi[] = { (cctki2_ibhi_),(cctki2_jbhi_),(cctki2_kbhi_) }; \
    const int cctki2_bbox[] = { (cctki2_ibboxlo_), (cctki2_ibboxhi_),(cctki2_jbboxlo_), (cctki2_jbboxhi_),(cctki2_kbboxlo_), (cctki2_kbboxhi_) }; \
    const int cctki2_lsh[] = { cctki2_cctkGH->cctk_lsh[0],cctki2_cctkGH->cctk_lsh[1],cctki2_cctkGH->cctk_lsh[2] }; \
    const int cctki2_istr CCTK_ATTRIBUTE_UNUSED = (cctki2_istr_); \
    for (int cctki2_kdir=-1; cctki2_kdir<=+1; ++cctki2_kdir) { \
    for (int cctki2_jdir=-1; cctki2_jdir<=+1; ++cctki2_jdir) { \
    for (int cctki2_idir=-1; cctki2_idir<=+1; ++cctki2_idir) { \
      const int cctki2_any_bbox = \
        (cctki2_idir<0 ? cctki2_bbox[0] : 0) || (cctki2_idir>0 ? cctki2_bbox[1] : 0) || \
        (cctki2_jdir<0 ? cctki2_bbox[2] : 0) || (cctki2_jdir>0 ? cctki2_bbox[3] : 0) || \
        (cctki2_kdir<0 ? cctki2_bbox[4] : 0) || (cctki2_kdir>0 ? cctki2_bbox[5] : 0); \
      const int cctki2_all_bbox = \
        (cctki2_idir<0 ? cctki2_bbox[0] : 1) && (cctki2_idir>0 ? cctki2_bbox[1] : 1) && \
        (cctki2_jdir<0 ? cctki2_bbox[2] : 1) && (cctki2_jdir>0 ? cctki2_bbox[3] : 1) && \
        (cctki2_kdir<0 ? cctki2_bbox[4] : 1) && (cctki2_kdir>0 ? cctki2_bbox[5] : 1); \
      if (cctki2_all_bbox && cctki2_any_bbox) { \
        const int cctki2_bmin[] = { \
          cctki2_idir<0 ? 0 : cctki2_idir==0 ? cctki2_blo[0] : cctki2_lsh[0] - cctki2_bhi[0], \
          cctki2_jdir<0 ? 0 : cctki2_jdir==0 ? cctki2_blo[1] : cctki2_lsh[1] - cctki2_bhi[1], \
          cctki2_kdir<0 ? 0 : cctki2_kdir==0 ? cctki2_blo[2] : cctki2_lsh[2] - cctki2_bhi[2], \
        }; \
        const int cctki2_bmax[] = { \
          cctki2_idir<0 ? cctki2_blo[0] : cctki2_idir==0 ? cctki2_lsh[0] - cctki2_bhi[0] : cctki2_lsh[0], \
          cctki2_jdir<0 ? cctki2_blo[1] : cctki2_jdir==0 ? cctki2_lsh[1] - cctki2_bhi[1] : cctki2_lsh[1], \
          cctki2_kdir<0 ? cctki2_blo[2] : cctki2_kdir==0 ? cctki2_lsh[2] - cctki2_bhi[2] : cctki2_lsh[2], \
        }; \
        CCTK_LOOP3STR_NORMAL(name##_intboundaries, \
                             i,j,k, \
                             ni,nj,nk, \
                             cctki2_idir,cctki2_jdir,cctki2_kdir, \
                             cctki2_bmin[0],cctki2_bmin[1],cctki2_bmin[2], \
                             cctki2_bmax[0],cctki2_bmax[1],cctki2_bmax[2], \
                             cctki2_cctkGH->cctk_ash[0], \
                             cctki2_cctkGH->cctk_ash[1], \
                             cctki2_cctkGH->cctk_ash[2], \
                             imin,imax, cctki2_istr) { \

#define CCTK_ENDLOOP3STR_INTBOUNDARIES(name) \
        } CCTK_ENDLOOP3STR_NORMAL(name##_intboundaries); \
      } /* if bbox */ \
    } /* for dir */ \
    } /* for dir */ \
    } /* for dir */ \
    typedef cctki2_loop3_intboundaries_##name cctki2_ensure_proper_nesting CCTK_ATTRIBUTE_UNUSED; \
  } while (0) \



/* LOOP_ALL */

#define CCTK_LOOP3_ALL(name, cctki3_cctkGH_, \
                       i,j,k) \
  CCTK_LOOP3STR_ALL(name, (cctki3_cctkGH_), \
                    i,j,k, \
                    cctki3_dummy_imin,cctki3_dummy_imax, 1) \

#define CCTK_ENDLOOP3_ALL(name) \
  CCTK_ENDLOOP3STR_ALL(name) \

#define CCTK_LOOP3STR_ALL(name, cctki3_cctkGH_, \
                          i,j,k, \
                          imin,imax, cctki3_istr_) \
  do { \
    typedef int cctki3_loop3_all_##name; \
    cGH const *restrict const cctki3_cctkGH = (cctki3_cctkGH_); \
    if (cctki3_cctkGH->cctk_dim != 3) { \
      _Pragma("omp critical") \
      CCTK_WARN(CCTK_WARN_ABORT, \
                "The macro CCTK_LOOP3_ALL can only be used in 3 dimensions"); \
    } \
    CCTK_LOOP3STR(name##_all, \
                  i,j,k, \
                  0,0,0, \
                  cctki3_cctkGH->cctk_lsh[0], \
                  cctki3_cctkGH->cctk_lsh[1], \
                  cctki3_cctkGH->cctk_lsh[2], \
                  cctki3_cctkGH->cctk_ash[0], \
                  cctki3_cctkGH->cctk_ash[1], \
                  cctki3_cctkGH->cctk_ash[2], \
                  imin,imax, (cctki3_istr_)) { \

#define CCTK_ENDLOOP3STR_ALL(name) \
    } CCTK_ENDLOOP3STR(name##_all); \
    typedef cctki3_loop3_all_##name cctki3_ensure_proper_nesting CCTK_ATTRIBUTE_UNUSED; \
  } while (0) \



/* LOOP_INT */

#define CCTK_LOOP3_INT(name, cctki3_cctkGH_, \
                       i,j,k) \
  CCTK_LOOP3STR_INT(name, (cctki3_cctkGH_), \
                    i,j,k, \
                    cctki3_dummy_imin,cctki3_dummy_imax, 1) \

#define CCTK_ENDLOOP3_INT(name) \
  CCTK_ENDLOOP3STR_INT(name) \

#define CCTK_LOOP3STR_INT(name, cctki3_cctkGH_, \
                          i,j,k, \
                          imin,imax, cctki3_istr_) \
  do { \
    typedef int cctki3_loop3_int_##name; \
    cGH const *restrict const cctki3_cctkGH = (cctki3_cctkGH_); \
    if (cctki3_cctkGH->cctk_dim != 3) { \
      _Pragma("omp critical") \
      CCTK_WARN(CCTK_WARN_ABORT, \
                "The macro CCTK_LOOP3_INT can only be used in 3 dimensions"); \
    } \
    CCTK_INT cctki3_bndsize    [6]; \
    CCTK_INT cctki3_is_ghostbnd[6]; \
    CCTK_INT cctki3_is_symbnd  [6]; \
    CCTK_INT cctki3_is_physbnd [6]; \
    _Pragma("omp single copyprivate(cctki3_bndsize)") \
    GetBoundarySizesAndTypes \
      (cctki3_cctkGH, 6, cctki3_bndsize, cctki3_is_ghostbnd, cctki3_is_symbnd, cctki3_is_physbnd); \
    CCTK_LOOP3STR_INTERIOR(name##_int, \
                           cctki3_cctkGH, \
                           i,j,k, \
                           cctki3_bndsize[0],cctki3_bndsize[2],cctki3_bndsize[4], \
                           cctki3_bndsize[1],cctki3_bndsize[3],cctki3_bndsize[5], \
                           imin,imax, (cctki3_istr_)) { \

#define CCTK_ENDLOOP3STR_INT(name) \
    } CCTK_ENDLOOP3STR_INTERIOR(name##_int); \
    typedef cctki3_loop3_int_##name cctki3_ensure_proper_nesting CCTK_ATTRIBUTE_UNUSED; \
  } while (0) \



/* LOOP_BND */

#define CCTK_LOOP3_BND(name, cctki3_cctkGH_, \
                       i,j,k, \
                       ni,nj,nk) \
  CCTK_LOOP3STR_BND(name, (cctki3_cctkGH_), \
                    i,j,k, \
                    ni,nj,nk, \
                    cctki3_dummy_imin,cctki3_dummy_imax, 1) \

#define CCTK_ENDLOOP3_BND(name) \
  CCTK_ENDLOOP3STR_BND(name) \

#define CCTK_LOOP3STR_BND(name, cctki3_cctkGH_, \
                          i,j,k, \
                          ni,nj,nk, \
                          imin,imax, cctki3_istr_) \
  do { \
    typedef int cctki3_loop3_bnd_##name; \
    cGH const *restrict const cctki3_cctkGH = (cctki3_cctkGH_); \
    if (cctki3_cctkGH->cctk_dim != 3) { \
      _Pragma("omp critical") \
      CCTK_WARN(CCTK_WARN_ABORT, \
                "The macro CCTK_LOOP3_BND can only be used in 3 dimensions"); \
    } \
    CCTK_INT cctki3_bndsize    [6]; \
    CCTK_INT cctki3_is_ghostbnd[6]; \
    CCTK_INT cctki3_is_symbnd  [6]; \
    CCTK_INT cctki3_is_physbnd [6]; \
    _Pragma("omp single copyprivate(cctki3_bndsize, cctki3_is_physbnd)") \
    GetBoundarySizesAndTypes \
      (cctki3_cctkGH, 6, cctki3_bndsize, cctki3_is_ghostbnd, cctki3_is_symbnd, cctki3_is_physbnd); \
    CCTK_LOOP3STR_BOUNDARIES(name##_bnd, \
                             cctki3_cctkGH, \
                             i,j,k, \
                             ni,nj,nk, \
                             cctki3_bndsize[0],cctki3_bndsize[2],cctki3_bndsize[4], \
                             cctki3_bndsize[1],cctki3_bndsize[3],cctki3_bndsize[5], \
                             cctki3_is_physbnd[0],cctki3_is_physbnd[2],cctki3_is_physbnd[4], \
                             cctki3_is_physbnd[1],cctki3_is_physbnd[3],cctki3_is_physbnd[5], \
                             imin,imax, (cctki3_istr_)) { \

#define CCTK_ENDLOOP3STR_BND(name) \
    } CCTK_ENDLOOP3STR_BOUNDARIES(name##_bnd); \
    typedef cctki3_loop3_bnd_##name cctki3_ensure_proper_nesting CCTK_ATTRIBUTE_UNUSED; \
  } while (0) \



/* LOOP_INTBND */

#define CCTK_LOOP3_INTBND(name, cctki3_cctkGH_, \
                           i,j,k, \
                           ni,nj,nk) \
  CCTK_LOOP3STR_INTBND(name, (cctki3_cctkGH_), \
                        i,j,k, \
                        ni,nj,nk, \
                        cctki3_dummy_imin,cctki3_dummy_imax, 1) \

#define CCTK_ENDLOOP3_INTBND(name) \
  CCTK_ENDLOOP3STR_INTBND(name) \

#define CCTK_LOOP3STR_INTBND(name, cctki3_cctkGH_, \
                              i,j,k, \
                              ni,nj,nk, \
                              imin,imax, cctki3_istr_) \
  do { \
    typedef int cctki3_loop3_intbnd_##name; \
    cGH const *restrict const cctki3_cctkGH = (cctki3_cctkGH_); \
    if (cctki3_cctkGH->cctk_dim != 3) { \
      _Pragma("omp critical") \
      CCTK_WARN(CCTK_WARN_ABORT, \
                "The macro CCTK_LOOP3_INTBND can only be used in 3 dimensions"); \
    } \
    CCTK_INT cctki3_bndsize    [6]; \
    CCTK_INT cctki3_is_ghostbnd[6]; \
    CCTK_INT cctki3_is_symbnd  [6]; \
    CCTK_INT cctki3_is_physbnd [6]; \
    _Pragma("omp single copyprivate(cctki3_bndsize, cctki3_is_physbnd)") \
    GetBoundarySizesAndTypes \
      (cctki3_cctkGH, 6, cctki3_bndsize, cctki3_is_ghostbnd, cctki3_is_symbnd, cctki3_is_physbnd); \
    CCTK_LOOP3STR_INTBOUNDARIES(name##_intbnd, \
                                cctki3_cctkGH, \
                                i,j,k, \
                                ni,nj,nk, \
                                cctki3_bndsize[0],cctki3_bndsize[2],cctki3_bndsize[4], \
                                cctki3_bndsize[1],cctki3_bndsize[3],cctki3_bndsize[5], \
                                cctki3_is_physbnd[0],cctki3_is_physbnd[2],cctki3_is_physbnd[4], \
                                cctki3_is_physbnd[1],cctki3_is_physbnd[3],cctki3_is_physbnd[5], \
                                imin,imax, (cctki3_istr_)) { \

#define CCTK_ENDLOOP3STR_INTBND(name) \
    } CCTK_ENDLOOP3STR_INTBOUNDARIES(name##_intbnd); \
    typedef cctki3_loop3_intbnd_##name cctki3_ensure_proper_nesting CCTK_ATTRIBUTE_UNUSED; \
  } while (0) \

#endif /* #ifdef CCODE */



#ifdef FCODE

/* LOOP */

#define CCTK_LOOP3_NORMAL_DECLARE(name) \
   CCTK_LOOP3STR_NORMAL_DECLARE(name) \
   && integer :: name/**/0_dummy_imin, name/**/0_dummy_imax \

#define CCTK_LOOP3_NORMAL_OMP_PRIVATE(name) \
   CCTK_LOOP3STR_NORMAL_OMP_PRIVATE(name) \

#define CCTK_LOOP3_NORMAL(name, \
                          i,j,k, \
                          ni,nj,nk, \
                          cctki0_idir,cctki0_jdir,cctki0_kdir, \
                          cctki0_imin,cctki0_jmin,cctki0_kmin, \
                          cctki0_imax,cctki0_jmax,cctki0_kmax, \
                          cctki0_iash,cctki0_jash,cctki0_kash) \
   CCTK_LOOP3STR_NORMAL(name, \
                        i,j,k, \
                        ni,nj,nk, \
                        cctki0_idir,cctki0_jdir,cctki0_kdir, \
                        cctki0_imin,cctki0_jmin,cctki0_kmin, \
                        cctki0_imax,cctki0_jmax,cctki0_kmax, \
                        cctki0_iash,cctki0_jash,cctki0_kash, \
                        name/**/0_dummy_imin,name/**/0_dummy_imax, 1) \

#define CCTK_ENDLOOP3_NORMAL(name) \
   CCTK_ENDLOOP3STR_NORMAL(name) \

#define CCTK_LOOP3STR_NORMAL_DECLARE(name) \
   && integer :: name/**/0_idir,name/**/0_jdir,name/**/0_kdir \
   && integer :: name/**/0_imin,name/**/0_jmin,name/**/0_kmin \
   && integer :: name/**/0_imax,name/**/0_jmax,name/**/0_kmax \
   && integer :: name/**/0_iash,name/**/0_jash,name/**/0_kash \
   && integer :: name/**/0_istr \

#define CCTK_LOOP3STR_NORMAL_OMP_PRIVATE(name) \
   && !$omp private (i,j,k) \
   && !$omp private (ni,nj,nk) \

#define CCTK_LOOP3STR_NORMAL(name, \
                             i,j,k, \
                             ni,nj,nk, \
                             cctki0_idir,cctki0_jdir,cctki0_kdir, \
                             cctki0_imin,cctki0_jmin,cctki0_kmin, \
                             cctki0_imax,cctki0_jmax,cctki0_kmax, \
                             cctki0_iash,cctki0_jash,cctki0_kash, \
                             imin,imax, cctki0_istr) \
   && name/**/0_idir = cctki0_idir \
   && name/**/0_jdir = cctki0_jdir \
   && name/**/0_kdir = cctki0_kdir \
   && name/**/0_imin = cctki0_imin \
   && name/**/0_jmin = cctki0_jmin \
   && name/**/0_kmin = cctki0_kmin \
   && name/**/0_imax = cctki0_imax \
   && name/**/0_jmax = cctki0_jmax \
   && name/**/0_kmax = cctki0_kmax \
   && name/**/0_iash = cctki0_iash \
   && name/**/0_jash = cctki0_jash \
   && name/**/0_kash = cctki0_kash \
   && name/**/0_istr = cctki0_istr \
   && imin = name/**/0_imin \
   && imax = name/**/0_imax \
   && !$omp do collapse(2) \
   && do k = name/**/0_kmin, name/**/0_kmax \
   && do j = name/**/0_jmin, name/**/0_jmax \
   && do i = name/**/0_imin - modulo((imin+name/**/0_iash*(j+name/**/0_jash*(k))), name/**/0_istr), name/**/0_imax, name/**/0_istr \
   &&    ni = 0 \
   &&    nj = 0 \
   &&    nk = 0 \
   &&    if (name/**/0_idir < 0) ni = i \
   &&    if (name/**/0_jdir < 0) nj = j \
   &&    if (name/**/0_kdir < 0) nk = k \
   &&    if (name/**/0_idir > 0) ni = name/**/0_imax+1-i \
   &&    if (name/**/0_jdir > 0) nj = name/**/0_jmax+1-j \
   &&    if (name/**/0_kdir > 0) nk = name/**/0_kmax+1-k \

#define CCTK_ENDLOOP3STR_NORMAL(name) \
   && end do \
   && end do \
   && end do \



#define CCTK_LOOP3_DECLARE(name) \
   CCTK_LOOP3STR_DECLARE(name) \
   && integer :: name/**/1_dummy_imin, name/**/1_dummy_imax \

#define CCTK_LOOP3_OMP_PRIVATE(name) \
   CCTK_LOOP3STR_OMP_PRIVATE(name) \

#define CCTK_LOOP3(name, \
                   i,j,k, \
                   cctki0_imin,cctki0_jmin,cctki0_kmin, \
                   cctki0_imax,cctki0_jmax,cctki0_kmax, \
                   cctki0_iash,cctki0_jash,cctki0_kash) \
   CCTK_LOOP3STR(name, \
                 i,j,k, \
                 cctki0_imin,cctki0_jmin,cctki0_kmin, \
                 cctki0_imax,cctki0_jmax,cctki0_kmax, \
                 cctki0_iash,cctki0_jash,cctki0_kash, \
                 name/**/1_dummy_imin,name/**/1_dummy_imax, 1) \

#define CCTK_ENDLOOP3(name) \
   CCTK_ENDLOOP3STR(name) \

#define CCTK_LOOP3STR_DECLARE(name) \
   CCTK_LOOP3STR_NORMAL_DECLARE(name) \
   && integer :: name/**/1_ni,name/**/1_nj,name/**/1_nk \

#define CCTK_LOOP3STR_OMP_PRIVATE(name) \
   CCTK_LOOP3STR_NORMAL_OMP_PRIVATE(name) \

#define CCTK_LOOP3STR(name, \
                      i,j,k, \
                      cctki1_imin,cctki1_jmin,cctki1_kmin, \
                      cctki1_imax,cctki1_jmax,cctki1_kmax, \
                      cctki1_iash,cctki1_jash,cctki1_kash, \
                      imin,imax, cctki1_istr) \
   CCTK_LOOP3STR_NORMAL(name, \
                        i,j,k, \
                        name/**/1_ni,name/**/1_nj,name/**/1_nk, \
                        0,0,0, \
                        cctki1_imin,cctki1_jmin,cctki1_kmin, \
                        cctki1_imax,cctki1_jmax,cctki1_kmax, \
                        cctki1_iash,cctki1_jash,cctki1_kash, \
                        imin,imax, cctki1_istr) \

#define CCTK_ENDLOOP3STR(name) \
   CCTK_ENDLOOP3STR_NORMAL(name) \



/* LOOP_INTERIOR */

#define CCTK_LOOP3_INTERIOR_DECLARE(name) \
   CCTK_LOOP3STR_INTERIOR_DECLARE(name) \
   && integer :: name/**/2_dummy_imin, name/**/2_dummy_imax \

#define CCTK_LOOP3_INTERIOR_OMP_PRIVATE(name) \
   CCTK_LOOP3STR_INTERIOR_OMP_PRIVATE(name) \

#define CCTK_LOOP3_INTERIOR(name, \
                            i,j,k, \
                            cctki2_iblo,cctki2_jblo,cctki2_kblo, \
                            cctki2_ibhi,cctki2_jbhi,cctki2_kbhi) \
   CCTK_LOOP3STR_INTERIOR(name, \
                          i,j,k, \
                          cctki2_iblo,cctki2_jblo,cctki2_kblo, \
                          cctki2_ibhi,cctki2_jbhi,cctki2_kbhi, \
                          name/**/2_dummy_imin,name/**/2_dummy_imax, 1) \

#define CCTK_ENDLOOP3_INTERIOR(name) \
   CCTK_ENDLOOP3STR_INTERIOR(name) \

#define CCTK_LOOP3STR_INTERIOR_DECLARE(name) \
   CCTK_LOOP3STR_DECLARE(name/**/_interior) \

#define CCTK_LOOP3STR_INTERIOR_OMP_PRIVATE(name) \
   CCTK_LOOP3STR_OMP_PRIVATE(name/**/_interior) \

#define CCTK_LOOP3STR_INTERIOR(name, \
                               i,j,k, \
                               cctki2_iblo,cctki2_jblo,cctki2_kblo, \
                               cctki2_ibhi,cctki2_jbhi,cctki2_kbhi, \
                               imin,imax, cctki2_istr) \
   CCTK_LOOP3STR(name/**/_interior, \
                 i,j,k, \
                 (cctki2_iblo)+1, \
                 (cctki2_jblo)+1, \
                 (cctki2_kblo)+1, \
                 cctk_lsh(1)-(cctki2_ibhi), \
                 cctk_lsh(2)-(cctki2_jbhi), \
                 cctk_lsh(3)-(cctki2_kbhi), \
                 cctk_ash(1),cctk_ash(2),cctk_ash(3), \
                 imin,imax, cctki2_istr) \

#define CCTK_ENDLOOP3STR_INTERIOR(name) \
   CCTK_ENDLOOP3STR(name/**/_interior) \



/* LOOP_BOUNDARIES */

#define CCTK_LOOP3_BOUNDARIES_DECLARE(name) \
   CCTK_LOOP3STR_BOUNDARIES_DECLARE(name) \
   && integer :: name/**/2_dummy_imin, name/**/2_dummy_imax \

#define CCTK_LOOP3_BOUNDARIES_OMP_PRIVATE(name) \
   CCTK_LOOP3STR_BOUNDARIES_OMP_PRIVATE(name) \

#define CCTK_LOOP3_BOUNDARIES(name, \
                              i,j,k, \
                              ni,nj,nk, \
                              cctki2_iblo,cctki2_jblo,cctki2_kblo, \
                              cctki2_ibhi,cctki2_jbhi,cctki2_kbhi, \
                              cctki2_ibboxlo,cctki2_jbboxlo,cctki2_kbboxlo, \
                              cctki2_ibboxhi,cctki2_jbboxhi,cctki2_kbboxhi) \
   CCTK_LOOP3STR_BOUNDARIES(name, \
                            i,j,k, \
                            ni,nj,nk, \
                            cctki2_iblo,cctki2_jblo,cctki2_kblo, \
                            cctki2_ibhi,cctki2_jbhi,cctki2_kbhi, \
                            cctki2_ibboxlo,cctki2_jbboxlo,cctki2_kbboxlo, \
                            cctki2_ibboxhi,cctki2_jbboxhi,cctki2_kbboxhi, \
                            name/**/2_dummy_imin,name/**/2_dummy_imax, 1) \

#define CCTK_ENDLOOP3_BOUNDARIES(name) \
   CCTK_ENDLOOP3STR_BOUNDARIES(name) \

#define CCTK_LOOP3STR_BOUNDARIES_DECLARE(name) \
   CCTK_LOOP3STR_NORMAL_DECLARE(name/**/_boundaries) \
   && integer :: name/**/2_blo(3), name/**/2_bhi(3) \
   && integer :: name/**/2_bboxlo(3), name/**/2_bboxhi(3) \
   && integer :: name/**/2_istr \
   && integer :: name/**/2_idir \
   && integer :: name/**/2_jdir \
   && integer :: name/**/2_kdir \
   && logical :: name/**/2_any_bbox \
   && integer :: name/**/2_bmin(3), name/**/2_bmax(3) \

#define CCTK_LOOP3STR_BOUNDARIES_OMP_PRIVATE(name) \
   CCTK_LOOP3STR_NORMAL_OMP_PRIVATE(name/**/_boundaries) \
   && !$omp private (name/**/2_bmin, name/**/2_bmax) \

#define CCTK_LOOP3STR_BOUNDARIES(name, \
                                 i,j,k, \
                                 ni,nj,nk, \
                                 cctki2_iblo,cctki2_jblo,cctki2_kblo, \
                                 cctki2_ibhi,cctki2_jbhi,cctki2_kbhi, \
                                 cctki2_ibboxlo,cctki2_jbboxlo,cctki2_kbboxlo, \
                                 cctki2_ibboxhi,cctki2_jbboxhi,cctki2_kbboxhi, \
                                 imin,imax, cctki2_istr) \
   && name/**/2_blo = (/ cctki2_iblo,cctki2_jblo,cctki2_kblo /) \
   && name/**/2_bhi = (/ cctki2_ibhi,cctki2_jbhi,cctki2_kbhi /) \
   && name/**/2_bboxlo = (/ cctki2_ibboxlo,cctki2_jbboxlo,cctki2_kbboxlo /) \
   && name/**/2_bboxhi = (/ cctki2_ibboxhi,cctki2_jbboxhi,cctki2_kbboxhi /) \
   && name/**/2_istr = (cctki2_istr) \
   && do name/**/2_kdir=-1, +1 \
   && do name/**/2_jdir=-1, +1 \
   && do name/**/2_idir=-1, +1 \
   &&     name/**/2_any_bbox = .false. \
   &&     if (name/**/2_idir<0) name/**/2_any_bbox = name/**/2_any_bbox .or. name/**/2_bboxlo(1) /= 0 \
   &&     if (name/**/2_jdir<0) name/**/2_any_bbox = name/**/2_any_bbox .or. name/**/2_bboxlo(2) /= 0 \
   &&     if (name/**/2_kdir<0) name/**/2_any_bbox = name/**/2_any_bbox .or. name/**/2_bboxlo(3) /= 0 \
   &&     if (name/**/2_idir>0) name/**/2_any_bbox = name/**/2_any_bbox .or. name/**/2_bboxhi(1) /= 0 \
   &&     if (name/**/2_jdir>0) name/**/2_any_bbox = name/**/2_any_bbox .or. name/**/2_bboxhi(2) /= 0 \
   &&     if (name/**/2_kdir>0) name/**/2_any_bbox = name/**/2_any_bbox .or. name/**/2_bboxhi(3) /= 0 \
   &&    if (name/**/2_any_bbox) then \
   &&       name/**/2_bmin(1) = name/**/2_blo(1)+1 \
   &&       name/**/2_bmin(2) = name/**/2_blo(2)+1 \
   &&       name/**/2_bmin(3) = name/**/2_blo(3)+1 \
   &&       if (name/**/2_idir<0) name/**/2_bmin(1) = 1 \
   &&       if (name/**/2_jdir<0) name/**/2_bmin(2) = 1 \
   &&       if (name/**/2_kdir<0) name/**/2_bmin(3) = 1 \
   &&       if (name/**/2_idir>0) name/**/2_bmin(1) = cctk_lsh(1) - name/**/2_bhi(1) \
   &&       if (name/**/2_jdir>0) name/**/2_bmin(2) = cctk_lsh(2) - name/**/2_bhi(2) \
   &&       if (name/**/2_kdir>0) name/**/2_bmin(3) = cctk_lsh(3) - name/**/2_bhi(3) \
   &&       name/**/2_bmax(1) = cctk_lsh(1) - name/**/2_bhi(1) \
   &&       name/**/2_bmax(2) = cctk_lsh(2) - name/**/2_bhi(2) \
   &&       name/**/2_bmax(3) = cctk_lsh(3) - name/**/2_bhi(3) \
   &&       if (name/**/2_idir<0) name/**/2_bmax(1) = name/**/2_blo(1) \
   &&       if (name/**/2_jdir<0) name/**/2_bmax(2) = name/**/2_blo(2) \
   &&       if (name/**/2_kdir<0) name/**/2_bmax(3) = name/**/2_blo(3) \
   &&       if (name/**/2_idir>0) name/**/2_bmax(1) = cctk_lsh(1) \
   &&       if (name/**/2_jdir>0) name/**/2_bmax(2) = cctk_lsh(2) \
   &&       if (name/**/2_kdir>0) name/**/2_bmax(3) = cctk_lsh(3) \
   &&       CCTK_LOOP3STR_NORMAL(name/**/_boundaries, \
                                 i,j,k, \
                                 ni,nj,nk, \
                                 name/**/2_idir,name/**/2_jdir,name/**/2_kdir, \
                                 name/**/2_bmin(1),name/**/2_bmin(2),name/**/2_bmin(3), \
                                 name/**/2_bmax(1),name/**/2_bmax(2),name/**/2_bmax(3), \
                                 cctk_ash(1), \
                                 cctk_ash(2), \
                                 cctk_ash(3), \
                                 imin,imax, name/**/2_istr) \

#define CCTK_ENDLOOP3STR_BOUNDARIES(name) \
            CCTK_ENDLOOP3STR_NORMAL(name/**/_boundaries) \
   &&    end if /* bbox */ \
   && end do /* dir */ \
   && end do /* dir */ \
   && end do /* dir */ \



/* LOOP_INTBOUNDARIES */

#define CCTK_LOOP3_INTBOUNDARIES_DECLARE(name) \
   CCTK_LOOP3STR_INTBOUNDARIES_DECLARE(name) \
   && integer :: name/**/2_dummy_imin, name/**/2_dummy_imax \

#define CCTK_LOOP3_INTBOUNDARIES_OMP_PRIVATE(name) \
   CCTK_LOOP3STR_INTBOUNDARIES_OMP_PRIVATE(name) \

#define CCTK_LOOP3_INTBOUNDARIES(name, \
                                 i,j,k, \
                                 ni,nj,nk, \
                                 cctki2_iblo,cctki2_jblo,cctki2_kblo, \
                                 cctki2_ibhi,cctki2_jbhi,cctki2_kbhi, \
                                 cctki2_ibboxlo,cctki2_jbboxlo,cctki2_kbboxlo, \
                                 cctki2_ibboxhi,cctki2_jbboxhi,cctki2_kbboxhi) \
   CCTK_LOOP3STR_INTBOUNDARIES(name, \
                               i,j,k, \
                               ni,nj,nk, \
                               cctki2_iblo,cctki2_jblo,cctki2_kblo, \
                               cctki2_ibhi,cctki2_jbhi,cctki2_kbhi, \
                               cctki2_ibboxlo,cctki2_jbboxlo,cctki2_kbboxlo, \
                               cctki2_ibboxhi,cctki2_jbboxhi,cctki2_kbboxhi, \
                               name/**/2_dummy_imin,name/**/2_dummy_max, 1) \

#define CCTK_ENDLOOP3_INTBOUNDARIES(name) \
   CCTK_ENDLOOP3STR_INTBOUNDARIES(name) \

#define CCTK_LOOP3STR_INTBOUNDARIES_DECLARE(name) \
   CCTK_LOOP3STR_NORMAL_DECLARE(name/**/_intboundaries) \
   && integer :: name/**/2_blo(3), name/**/2_bhi(3) \
   && integer :: name/**/2_bboxlo(3), name/**/2_bboxhi(3) \
   && integer :: name/**/2_istr \
   && integer :: name/**/2_idir \
   && integer :: name/**/2_jdir \
   && integer :: name/**/2_kdir \
   && logical :: name/**/2_any_bbox, name/**/2_all_bbox \
   && integer :: name/**/2_bmin(3), name/**/2_bmax(3) \

#define CCTK_LOOP3STR_INTBOUNDARIES_OMP_PRIVATE(name) \
   CCTK_LOOP3STR_NORMAL_OMP_PRIVATE(name/**/_intboundaries) \
   && !$omp private (name/**/2_any_bbox, name/**/2_all_bbox) \
   && !$omp private (name/**/2_bmin, name/**/2_bmax \

#define CCTK_LOOP3STR_INTBOUNDARIES(name, \
                                    i,j,k, \
                                    ni,nj,nk, \
                                    cctki2_iblo,cctki2_jblo,cctki2_kblo, \
                                    cctki2_ibhi,cctki2_jbhi,cctki2_kbhi, \
                                    cctki2_ibboxlo,cctki2_jbboxlo,cctki2_kbboxlo, \
                                    cctki2_ibboxhi,cctki2_jbboxhi,cctki2_kbboxhi, \
                                    imin,imax, cctki2_istr) \
   && name/**/2_blo = (/ cctki2_iblo,cctki2_jblo,cctki2_kblo /) \
   && name/**/2_bhi = (/ cctki2_ibhi,cctki2_jbhi,cctki2_kbhi /) \
   && name/**/2_bboxlo = (/ cctki2_ibboxlo,cctki2_jbboxlo,cctki2_kbboxlo /) \
   && name/**/2_bboxhi = (/ cctki2_ibboxhi,cctki2_jbboxhi,cctki2_kbboxhi /) \
   && name/**/2_istr = (cctki2_istr) \
   && do name/**/2_kdir=-1, +1 \
   && do name/**/2_jdir=-1, +1 \
   && do name/**/2_idir=-1, +1 \
   &&     name/**/2_any_bbox = .false. \
   &&     if (name/**/2_idir<0) name/**/2_any_bbox = name/**/2_any_bbox .or. name/**/2_bboxlo(1) /= 0 \
   &&     if (name/**/2_jdir<0) name/**/2_any_bbox = name/**/2_any_bbox .or. name/**/2_bboxlo(2) /= 0 \
   &&     if (name/**/2_kdir<0) name/**/2_any_bbox = name/**/2_any_bbox .or. name/**/2_bboxlo(3) /= 0 \
   &&     if (name/**/2_idir>0) name/**/2_any_bbox = name/**/2_any_bbox .or. name/**/2_bboxhi(1) /= 0 \
   &&     if (name/**/2_jdir>0) name/**/2_any_bbox = name/**/2_any_bbox .or. name/**/2_bboxhi(2) /= 0 \
   &&     if (name/**/2_kdir>0) name/**/2_any_bbox = name/**/2_any_bbox .or. name/**/2_bboxhi(3) /= 0 \
   &&     name/**/2_all_bbox = .true. \
   &&     if (name/**/2_idir<0) name/**/2_all_bbox = name/**/2_all_bbox .and. name/**/2_bboxlo(1) /= 0 \
   &&     if (name/**/2_jdir<0) name/**/2_all_bbox = name/**/2_all_bbox .and. name/**/2_bboxlo(2) /= 0 \
   &&     if (name/**/2_kdir<0) name/**/2_all_bbox = name/**/2_all_bbox .and. name/**/2_bboxlo(3) /= 0 \
   &&     if (name/**/2_idir>0) name/**/2_all_bbox = name/**/2_all_bbox .and. name/**/2_bboxhi(1) /= 0 \
   &&     if (name/**/2_jdir>0) name/**/2_all_bbox = name/**/2_all_bbox .and. name/**/2_bboxhi(2) /= 0 \
   &&     if (name/**/2_kdir>0) name/**/2_all_bbox = name/**/2_all_bbox .and. name/**/2_bboxhi(3) /= 0 \
   &&    if (name/**/2_all_bbox .and. name/**/2_any_bbox) then \
   &&       name/**/2_bmin(1) = name/**/2_blo(1)+1 \
   &&       name/**/2_bmin(2) = name/**/2_blo(2)+1 \
   &&       name/**/2_bmin(3) = name/**/2_blo(3)+1 \
   &&       if (name/**/2_idir<0) name/**/2_bmin(1) = 1 \
   &&       if (name/**/2_jdir<0) name/**/2_bmin(2) = 1 \
   &&       if (name/**/2_kdir<0) name/**/2_bmin(3) = 1 \
   &&       if (name/**/2_idir>0) name/**/2_bmin(1) = cctk_lsh(1) - name/**/2_bhi(1) \
   &&       if (name/**/2_jdir>0) name/**/2_bmin(2) = cctk_lsh(2) - name/**/2_bhi(2) \
   &&       if (name/**/2_kdir>0) name/**/2_bmin(3) = cctk_lsh(3) - name/**/2_bhi(3) \
   &&       name/**/2_bmax(1) = cctk_lsh(1) - name/**/2_bhi(1) \
   &&       name/**/2_bmax(2) = cctk_lsh(2) - name/**/2_bhi(2) \
   &&       name/**/2_bmax(3) = cctk_lsh(3) - name/**/2_bhi(3) \
   &&       if (name/**/2_idir<0) name/**/2_bmax(1) = name/**/2_blo(1) \
   &&       if (name/**/2_jdir<0) name/**/2_bmax(2) = name/**/2_blo(2) \
   &&       if (name/**/2_kdir<0) name/**/2_bmax(3) = name/**/2_blo(3) \
   &&       if (name/**/2_idir>0) name/**/2_bmax(1) = cctk_lsh(1) \
   &&       if (name/**/2_jdir>0) name/**/2_bmax(2) = cctk_lsh(2) \
   &&       if (name/**/2_kdir>0) name/**/2_bmax(3) = cctk_lsh(3) \
   &&       CCTK_LOOP3STR_NORMAL(name/**/_intboundaries, \
                                 i,j,k, \
                                 ni,nj,nk, \
                                 name/**/2_idir,name/**/2_jdir,name/**/2_kdir, \
                                 name/**/2_bmin(1),name/**/2_bmin(2),name/**/2_bmin(3), \
                                 name/**/2_bmax(1),name/**/2_bmax(2),name/**/2_bmax(3), \
                                 cctk_ash(1), \
                                 cctk_ash(2), \
                                 cctk_ash(3), \
                                 imin,imax, name/**/2_istr) \

#define CCTK_ENDLOOP3STR_INTBOUNDARIES(name) \
            CCTK_ENDLOOP3STR_NORMAL(name/**/_intboundaries) \
   &&    end if /* bbox */ \
   && end do /* dir */ \
   && end do /* dir */ \
   && end do /* dir */ \



/* LOOP_ALL */

#define CCTK_LOOP3_ALL_DECLARE(name) \
   CCTK_LOOP3STR_ALL_DECLARE(name) \
   && integer :: name/**/3_dummy_imin, name/**/3_dummy_imax \

#define CCTK_LOOP3_ALL_OMP_PRIVATE(name) \
   CCTK_LOOP3STR_ALL_OMP_PRIVATE(name) \

#define CCTK_LOOP3_ALL(name, \
                       i,j,k) \
   CCTK_LOOP3STR_ALL(name, \
                     i,j,k, \
                     name/**/3_dummy_imin,name/**/3_dummy_imax, 1) \

#define CCTK_ENDLOOP3_ALL(name) \
   CCTK_ENDLOOP3STR_ALL(name) \

#define CCTK_LOOP3STR_ALL_DECLARE(name) \
   CCTK_LOOP3STR_DECLARE(name/**/_all) \

#define CCTK_LOOP3STR_ALL_OMP_PRIVATE(name) \
   CCTK_LOOP3STR_OMP_PRIVATE(name/**/_all) \

#define CCTK_LOOP3STR_ALL(name, \
                          i,j,k, \
                          imin,imax, cctki3_istr) \
   CCTK_LOOP3STR(name/**/_all, \
                 i,j,k, \
                 1,1,1, \
                 cctk_lsh(1),cctk_lsh(2),cctk_lsh(3), \
                 cctk_ash(1),cctk_ash(2),cctk_ash(3), \
                 imin,imax, cctki3_istr) \

#define CCTK_ENDLOOP3STR_ALL(name) \
   CCTK_ENDLOOP3STR(name/**/_all) \



/* LOOP_INT */

#define CCTK_LOOP3_INT_DECLARE(name) \
   CCTK_LOOP3STR_INT_DECLARE(name) \
   && integer :: name/**/3_dummy_imin, name/**/3_dummy_imax \

#define CCTK_LOOP3_INT_OMP_PRIVATE(name) \
   CCTK_LOOP3STR_INT_OMP_PRIVATE(name) \

#define CCTK_LOOP3_INT(name, \
                        i,j,k) \
   CCTK_LOOP3STR_INT(name, \
                      i,j,k, \
                      name/**/3_dummy_imin,name/**/3_dummy_imax, 1) \

#define CCTK_ENDLOOP3_INT(name) \
   CCTK_ENDLOOP3STR_INT(name) \

#define CCTK_LOOP3STR_INT_DECLARE(name) \
   CCTK_LOOP3STR_INTERIOR_DECLARE(name/**/_int) \
   && CCTK_INT, parameter :: name/**/3_isize = 6 \
   && CCTK_INT :: name/**/3_bndsize    (6) \
   && CCTK_INT :: name/**/3_is_ghostbnd(6) \
   && CCTK_INT :: name/**/3_is_symbnd  (6) \
   && CCTK_INT :: name/**/3_is_physbnd (6) \
   && CCTK_INT :: name/**/3_ierr \

#define CCTK_LOOP3STR_INT_OMP_PRIVATE(name) \
   CCTK_LOOP3STR_INTERIOR_OMP_PRIVATE(name/**/_int) \

#define CCTK_LOOP3STR_INT(name, \
                          i,j,k, \
                          imin,imax, cctki3_istr) \
   && !$omp single \
   && name/**/3_ierr = GetBoundarySizesAndTypes \
         (cctkGH, name/**/3_isize, name/**/3_bndsize, name/**/3_is_ghostbnd, name/**/3_is_symbnd, name/**/3_is_physbnd) \
   && !$omp end single copyprivate(name/**/3_bndsize) \
   && CCTK_LOOP3STR_INTERIOR(name/**/_int, \
                             i,j,k, \
                             int(name/**/3_bndsize(1+1)),int(name/**/3_bndsize(3+1)),int(name/**/3_bndsize(5+1)), \
                             int(name/**/3_bndsize(2)),int(name/**/3_bndsize(4)),int(name/**/3_bndsize(6)), \
                             imin,imax, (cctki3_istr)) \

#define CCTK_ENDLOOP3STR_INT(name) \
      CCTK_ENDLOOP3STR_INTERIOR(name/**/_int) \



/* LOOP_BND */

#define CCTK_LOOP3_BND_DECLARE(name) \
   CCTK_LOOP3STR_BND_DECLARE(name) \
   && integer :: name/**/3_dummy_imin, name/**/3_dummy_imax \

#define CCTK_LOOP3_BND_OMP_PRIVATE(name) \
   CCTK_LOOP3STR_BND_OMP_PRIVATE(name) \

#define CCTK_LOOP3_BND(name, \
                       i,j,k, \
                       ni,nj,nk) \
   CCTK_LOOP3STR_BND(name, \
                     i,j,k, \
                     ni,nj,nk, \
                     name/**/3_dummy_imin,name/**/3_dummy_imax, 1) \

#define CCTK_ENDLOOP3_BND(name) \
   CCTK_ENDLOOP3STR_BND(name) \

#define CCTK_LOOP3STR_BND_DECLARE(name) \
   CCTK_LOOP3STR_BOUNDARIES_DECLARE(name/**/_bnd) \
   && CCTK_INT, parameter :: name/**/3_isize = 6 \
   && CCTK_INT :: name/**/3_bndsize    (6) \
   && CCTK_INT :: name/**/3_is_ghostbnd(6) \
   && CCTK_INT :: name/**/3_is_symbnd  (6) \
   && CCTK_INT :: name/**/3_is_physbnd (6) \
   && CCTK_INT :: name/**/3_ierr \

#define CCTK_LOOP3STR_BND_OMP_PRIVATE(name) \
   CCTK_LOOP3STR_BOUNDARIES_OMP_PRIVATE(name/**/_bnd) \

#define CCTK_LOOP3STR_BND(name, \
                          i,j,k, \
                          ni,nj,nk, \
                          imin,imax, cctki3_istr) \
   && !$omp single \
   && name/**/3_ierr = GetBoundarySizesAndTypes \
         (cctkGH, name/**/3_isize, name/**/3_bndsize, name/**/3_is_ghostbnd, name/**/3_is_symbnd, name/**/3_is_physbnd) \
   && !$omp end single copyprivate(name/**/3_bndsize, name/**/3_is_physbnd) \
   && CCTK_LOOP3STR_BOUNDARIES(name/**/_bnd, \
                               i,j,k, \
                               ni,nj,nk, \
                               int(name/**/3_bndsize(1))+1,int(name/**/3_bndsize(3))+1,int(name/**/3_bndsize(5))+1, \
                               int(name/**/3_bndsize(2)),int(name/**/3_bndsize(4)),int(name/**/3_bndsize(6)), \
                               int(name/**/3_is_physbnd(1)),int(name/**/3_is_physbnd(3)),int(name/**/3_is_physbnd(5)), \
                               int(name/**/3_is_physbnd(2)),int(name/**/3_is_physbnd(4)),int(name/**/3_is_physbnd(6)), \
                               imin,imax, (cctki3_istr)) \

#define CCTK_ENDLOOP3STR_BND(name) \
      CCTK_ENDLOOP3STR_BOUNDARIES(name/**/_bnd) \






/* LOOP_INTBND */

#define CCTK_LOOP3_INTBND_DECLARE(name) \
   CCTK_LOOP3STR_INTBND_DECLARE(name) \
   && integer :: name/**/3_dummy_imin, name/**/3_dummy_imax \

#define CCTK_LOOP3_INTBND_OMP_PRIVATE(name) \
   CCTK_LOOP3STR_INTBND_OMP_PRIVATE(name) \

#define CCTK_LOOP3_INTBND(name, \
                          i,j,k, \
                          ni,nj,nk) \
   CCTK_LOOP3STR_INTBND(name, \
                        i,j,k, \
                        ni,nj,nk, \
                        name/**/3_dummy_imin,name/**/3_dummy_imax, 1) \

#define CCTK_ENDLOOP3_INTBND(name) \
   CCTK_ENDLOOP3STR_INTBND(name) \

#define CCTK_LOOP3STR_INTBND_DECLARE(name) \
   CCTK_LOOP3STR_INTBOUNDARIES_DECLARE(name/**/_bnd) \
   && CCTK_INT, parameter :: name/**/3_isize = 6 \
   && CCTK_INT :: name/**/3_bndsize    (6) \
   && CCTK_INT :: name/**/3_is_ghostbnd(6) \
   && CCTK_INT :: name/**/3_is_symbnd  (6) \
   && CCTK_INT :: name/**/3_is_physbnd (6) \
   && CCTK_INT :: name/**/3_ierr \

#define CCTK_LOOP3STR_INTBND_OMP_PRIVATE(name) \
   CCTK_LOOP3STR_INTBOUNDARIES_OMP_PRIVATE(name/**/_bnd) \

#define CCTK_LOOP3STR_INTBND(name, \
                             i,j,k, \
                             ni,nj,nk, \
                             imin,imax, cctki3_istr) \
   && !$omp single \
   && name/**/3_ierr = GetBoundarySizesAndTypes \
         (cctkGH, name/**/3_isize, name/**/3_bndsize, name/**/3_is_ghostbnd, name/**/3_is_symbnd, name/**/3_is_physbnd) \
   && !$omp end single copyprivate(name/**/3_bndsize, name/**/3_is_physbnd) \
   && CCTK_LOOP3STR_INTBOUNDARIES(name/**/_bnd, \
                                  i,j,k, \
                                  ni,nj,nk, \
                                  int(name/**/3_bndsize(1+1)),int(name/**/3_bndsize(3+1)),int(name/**/3_bndsize(5+1)), \
                                  int(name/**/3_bndsize(2)),int(name/**/3_bndsize(4)),int(name/**/3_bndsize(6)), \
                                  int(name/**/3_is_physbnd(1)),int(name/**/3_is_physbnd(3)),int(name/**/3_is_physbnd(5)), \
                                  int(name/**/3_is_physbnd(2)),int(name/**/3_is_physbnd(4)),int(name/**/3_is_physbnd(6)), \
                                  imin,imax, (cctki3_istr)) \

#define CCTK_ENDLOOP3STR_INTBND(name) \
      CCTK_ENDLOOP3STR_INTBOUNDARIES(name/**/_bnd) \

#endif /* #ifdef FCODE */



/* 4D */

#ifdef CCODE

/* LOOP */

#define CCTK_LOOP4_NORMAL(name, \
                          i,j,k,l, \
                          ni,nj,nk,nl, \
                          cctki0_idir_,cctki0_jdir_,cctki0_kdir_,cctki0_ldir_, \
                          cctki0_imin_,cctki0_jmin_,cctki0_kmin_,cctki0_lmin_, \
                          cctki0_imax_,cctki0_jmax_,cctki0_kmax_,cctki0_lmax_, \
                          cctki0_iash_,cctki0_jash_,cctki0_kash_,cctki0_lash_) \
  CCTK_LOOP4STR_NORMAL(name, \
                       i,j,k,l, \
                       ni,nj,nk,nl, \
                       (cctki0_idir_),(cctki0_jdir_),(cctki0_kdir_),(cctki0_ldir_), \
                       (cctki0_imin_),(cctki0_jmin_),(cctki0_kmin_),(cctki0_lmin_), \
                       (cctki0_imax_),(cctki0_jmax_),(cctki0_kmax_),(cctki0_lmax_), \
                       (cctki0_iash_),(cctki0_jash_),(cctki0_kash_),(cctki0_lash_), \
                       cctki0_dummy_imin,cctki0_dummy_imax, 1) \

#define CCTK_ENDLOOP4_NORMAL(name) \
  CCTK_ENDLOOP4STR_NORMAL(name) \

#define CCTK_LOOP4STR_NORMAL(name, \
                             i,j,k,l, \
                             ni,nj,nk,nl, \
                             cctki0_idir_,cctki0_jdir_,cctki0_kdir_,cctki0_ldir_, \
                             cctki0_imin_,cctki0_jmin_,cctki0_kmin_,cctki0_lmin_, \
                             cctki0_imax_,cctki0_jmax_,cctki0_kmax_,cctki0_lmax_, \
                             cctki0_iash_,cctki0_jash_,cctki0_kash_,cctki0_lash_, \
                             imin,imax, cctki0_istr_) \
  do { \
    typedef int cctki0_loop4_normal_##name; \
    const int cctki0_idir = (cctki0_idir_); \
    const int cctki0_jdir = (cctki0_jdir_); \
    const int cctki0_kdir = (cctki0_kdir_); \
    const int cctki0_ldir = (cctki0_ldir_); \
    const int cctki0_imin = (cctki0_imin_); \
    const int cctki0_jmin = (cctki0_jmin_); \
    const int cctki0_kmin = (cctki0_kmin_); \
    const int cctki0_lmin = (cctki0_lmin_); \
    const int cctki0_imax = (cctki0_imax_); \
    const int cctki0_jmax = (cctki0_jmax_); \
    const int cctki0_kmax = (cctki0_kmax_); \
    const int cctki0_lmax = (cctki0_lmax_); \
    const int cctki0_iash CCTK_ATTRIBUTE_UNUSED = (cctki0_iash_); \
    const int cctki0_jash CCTK_ATTRIBUTE_UNUSED = (cctki0_jash_); \
    const int cctki0_kash CCTK_ATTRIBUTE_UNUSED = (cctki0_kash_); \
    const int cctki0_lash CCTK_ATTRIBUTE_UNUSED = (cctki0_lash_); \
    const int cctki0_istr = (cctki0_istr_); \
    assert(cctki0_istr>0 && (cctki0_istr & (cctki0_istr-1)) == 0); \
    const int imin CCTK_ATTRIBUTE_UNUSED = cctki0_imin; \
    const int imax CCTK_ATTRIBUTE_UNUSED = cctki0_imax; \
    CCTK_PRAGMA_OMP_FOR_COLLAPSE_3 \
    for (int l=cctki0_lmin; l<cctki0_lmax; ++l) { \
    for (int k=cctki0_kmin; k<cctki0_kmax; ++k) { \
    for (int j=cctki0_jmin; j<cctki0_jmax; ++j) { \
     \
    const int cctki0_ioff = (cctki0_imin+cctki0_iash*(j+cctki0_jash*(k+cctki0_kash*(l)))) & (cctki0_istr-1); \
     \
    for (int i=cctki0_imin-cctki0_ioff; i<cctki0_imax; i+=cctki0_istr) { \
      const int ni CCTK_ATTRIBUTE_UNUSED = cctki0_idir<0 ? i+1 : cctki0_idir==0 ? 0 : cctki0_imax-i; \
      const int nj CCTK_ATTRIBUTE_UNUSED = cctki0_jdir<0 ? j+1 : cctki0_jdir==0 ? 0 : cctki0_jmax-j; \
      const int nk CCTK_ATTRIBUTE_UNUSED = cctki0_kdir<0 ? k+1 : cctki0_kdir==0 ? 0 : cctki0_kmax-k; \
      const int nl CCTK_ATTRIBUTE_UNUSED = cctki0_ldir<0 ? l+1 : cctki0_ldir==0 ? 0 : cctki0_lmax-l; \
      { \

#define CCTK_ENDLOOP4STR_NORMAL(name) \
      } \
    } \
    } \
    } \
    } \
    typedef cctki0_loop4_normal_##name cctki0_ensure_proper_nesting CCTK_ATTRIBUTE_UNUSED; \
  } while (0) \



#define CCTK_LOOP4(name, \
                   i,j,k,l, \
                   cctki1_imin_,cctki1_jmin_,cctki1_kmin_,cctki1_lmin_, \
                   cctki1_imax_,cctki1_jmax_,cctki1_kmax_,cctki1_lmax_, \
                   cctki1_iash_,cctki1_jash_,cctki1_kash_,cctki1_lash_) \
  CCTK_LOOP4STR(name, \
                i,j,k,l, \
                (cctki1_imin_),(cctki1_jmin_),(cctki1_kmin_),(cctki1_lmin_), \
                (cctki1_imax_),(cctki1_jmax_),(cctki1_kmax_),(cctki1_lmax_), \
                (cctki1_iash_),(cctki1_jash_),(cctki1_kash_),(cctki1_lash_), \
                cctki1_dummy_imin,cctki1_dummy_imax, 1) \

#define CCTK_ENDLOOP4(name) \
  CCTK_ENDLOOP4STR(name) \

#define CCTK_LOOP4STR(name, \
                      i,j,k,l, \
                      cctki1_imin_,cctki1_jmin_,cctki1_kmin_,cctki1_lmin_, \
                      cctki1_imax_,cctki1_jmax_,cctki1_kmax_,cctki1_lmax_, \
                      cctki1_iash_,cctki1_jash_,cctki1_kash_,cctki1_lash_, \
                      imin,imax, cctki1_istr_) \
  CCTK_LOOP4STR_NORMAL(name, \
                       i,j,k,l, \
                       cctki1_ni,cctki1_nj,cctki1_nk,cctki1_nl, \
                       0,0,0,0, \
                       (cctki1_imin_),(cctki1_jmin_),(cctki1_kmin_),(cctki1_lmin_), \
                       (cctki1_imax_),(cctki1_jmax_),(cctki1_kmax_),(cctki1_lmax_), \
                       (cctki1_iash_),(cctki1_jash_),(cctki1_kash_),(cctki1_lash_), \
                       imin,imax, (cctki1_istr_)) \

#define CCTK_ENDLOOP4STR(name) \
  CCTK_ENDLOOP4STR_NORMAL(name) \



/* LOOP_INTERIOR */

#define CCTK_LOOP4_INTERIOR(name, cctki2_cctkGH_, \
                            i,j,k,l, \
                            cctki2_iblo_,cctki2_jblo_,cctki2_kblo_,cctki2_lblo_, \
                            cctki2_ibhi_,cctki2_jbhi_,cctki2_kbhi_,cctki2_lbhi_) \
  CCTK_LOOP4STR_INTERIOR(name, (cctki2_cctkGH_), \
                         i,j,k,l, \
                         (cctki2_iblo_),(cctki2_jblo_),(cctki2_kblo_),(cctki2_lblo_), \
                         (cctki2_ibhi_),(cctki2_jbhi_),(cctki2_kbhi_),(cctki2_lbhi_), \
                         cctki2_dummy_imin,cctki2_dummy_imax, 1) \

#define CCTK_ENDLOOP4_INTERIOR(name) \
  CCTK_ENDLOOP4STR_INTERIOR(name) \

#define CCTK_LOOP4STR_INTERIOR(name, cctki2_cctkGH_, \
                               i,j,k,l, \
                               cctki2_iblo_,cctki2_jblo_,cctki2_kblo_,cctki2_lblo_, \
                               cctki2_ibhi_,cctki2_jbhi_,cctki2_kbhi_,cctki2_lbhi_, \
                               imin,imax, cctki2_istr_) \
  do { \
    typedef int cctki2_loop4_interior_##name; \
    cGH const *restrict const cctki2_cctkGH = (cctki2_cctkGH_); \
    if (cctki2_cctkGH->cctk_dim != 4) { \
      _Pragma("omp critical") \
      CCTK_WARN(CCTK_WARN_ABORT, \
                "The macro CCTK_LOOP4_INTERIOR can only be used in 4 dimensions"); \
    } \
    CCTK_LOOP4STR(name##_interior, \
                  i,j,k,l, \
                  (cctki2_iblo_),(cctki2_jblo_),(cctki2_kblo_),(cctki2_lblo_), \
                  cctki2_cctkGH->cctk_lsh[0]-(cctki2_ibhi_), \
                  cctki2_cctkGH->cctk_lsh[1]-(cctki2_jbhi_), \
                  cctki2_cctkGH->cctk_lsh[2]-(cctki2_kbhi_), \
                  cctki2_cctkGH->cctk_lsh[3]-(cctki2_lbhi_), \
                  cctki2_cctkGH->cctk_ash[0], \
                  cctki2_cctkGH->cctk_ash[1], \
                  cctki2_cctkGH->cctk_ash[2], \
                  cctki2_cctkGH->cctk_ash[3], \
                  imin,imax, (cctki2_istr_)) { \

#define CCTK_ENDLOOP4STR_INTERIOR(name) \
    } CCTK_ENDLOOP4STR(name##_interior); \
    typedef cctki2_loop4_interior_##name cctki2_ensure_proper_nesting CCTK_ATTRIBUTE_UNUSED; \
  } while(0) \



/* LOOP_BOUNDARIES */

#define CCTK_LOOP4_BOUNDARIES(name, cctki2_cctkGH_, \
                              i,j,k,l, \
                              ni,nj,nk,nl, \
                              cctki2_iblo_,cctki2_jblo_,cctki2_kblo_,cctki2_lblo_, \
                              cctki2_ibhi_,cctki2_jbhi_,cctki2_kbhi_,cctki2_lbhi_, \
                              cctki2_ibboxlo_,cctki2_jbboxlo_,cctki2_kbboxlo_,cctki2_lbboxlo_, \
                              cctki2_ibboxhi_,cctki2_jbboxhi_,cctki2_kbboxhi_,cctki2_lbboxhi_) \
  CCTK_LOOP4STR_BOUNDARIES(name, (cctki2_cctkGH_), \
                           i,j,k,l, \
                           ni,nj,nk,nl, \
                           (cctki2_iblo_),(cctki2_jblo_),(cctki2_kblo_),(cctki2_lblo_), \
                           (cctki2_ibhi_),(cctki2_jbhi_),(cctki2_kbhi_),(cctki2_lbhi_), \
                           (cctki2_ibboxlo_),(cctki2_jbboxlo_),(cctki2_kbboxlo_),(cctki2_lbboxlo_), \
                           (cctki2_ibboxhi_),(cctki2_jbboxhi_),(cctki2_kbboxhi_),(cctki2_lbboxhi_), \
                           cctki2_dummy_imin,cctki2_dummy_imax, 1) \

#define CCTK_ENDLOOP4_BOUNDARIES(name) \
  CCTK_ENDLOOP4STR_BOUNDARIES(name) \

#define CCTK_LOOP4STR_BOUNDARIES(name, cctki2_cctkGH_, \
                                 i,j,k,l, \
                                 ni,nj,nk,nl, \
                                 cctki2_iblo_,cctki2_jblo_,cctki2_kblo_,cctki2_lblo_, \
                                 cctki2_ibhi_,cctki2_jbhi_,cctki2_kbhi_,cctki2_lbhi_, \
                                 cctki2_ibboxlo_,cctki2_jbboxlo_,cctki2_kbboxlo_,cctki2_lbboxlo_, \
                                 cctki2_ibboxhi_,cctki2_jbboxhi_,cctki2_kbboxhi_,cctki2_lbboxhi_, \
                                 imin,imax, cctki2_istr_) \
  do { \
    typedef int cctki2_loop4_boundaries_##name; \
    cGH const *restrict const cctki2_cctkGH = (cctki2_cctkGH_); \
    if (cctki2_cctkGH->cctk_dim != 4) { \
      _Pragma("omp critical") \
      CCTK_WARN(CCTK_WARN_ABORT, \
                "The macro CCTK_LOOP4_BOUNDARIES can only be used in 4 dimensions"); \
    } \
    const int cctki2_blo[] = { (cctki2_iblo_),(cctki2_jblo_),(cctki2_kblo_),(cctki2_lblo_) }; \
    const int cctki2_bhi[] = { (cctki2_ibhi_),(cctki2_jbhi_),(cctki2_kbhi_),(cctki2_lbhi_) }; \
    const int cctki2_bbox[] = { (cctki2_ibboxlo_), (cctki2_ibboxhi_),(cctki2_jbboxlo_), (cctki2_jbboxhi_),(cctki2_kbboxlo_), (cctki2_kbboxhi_),(cctki2_lbboxlo_), (cctki2_lbboxhi_) }; \
    const int cctki2_lsh[] = { cctki2_cctkGH->cctk_lsh[0],cctki2_cctkGH->cctk_lsh[1],cctki2_cctkGH->cctk_lsh[2],cctki2_cctkGH->cctk_lsh[3] }; \
    const int cctki2_istr CCTK_ATTRIBUTE_UNUSED = (cctki2_istr_); \
    for (int cctki2_ldir=-1; cctki2_ldir<=+1; ++cctki2_ldir) { \
    for (int cctki2_kdir=-1; cctki2_kdir<=+1; ++cctki2_kdir) { \
    for (int cctki2_jdir=-1; cctki2_jdir<=+1; ++cctki2_jdir) { \
    for (int cctki2_idir=-1; cctki2_idir<=+1; ++cctki2_idir) { \
      const int cctki2_any_bbox = \
        (cctki2_idir<0 ? cctki2_bbox[0] : 0) || (cctki2_idir>0 ? cctki2_bbox[1] : 0) || \
        (cctki2_jdir<0 ? cctki2_bbox[2] : 0) || (cctki2_jdir>0 ? cctki2_bbox[3] : 0) || \
        (cctki2_kdir<0 ? cctki2_bbox[4] : 0) || (cctki2_kdir>0 ? cctki2_bbox[5] : 0) || \
        (cctki2_ldir<0 ? cctki2_bbox[6] : 0) || (cctki2_ldir>0 ? cctki2_bbox[7] : 0); \
      if (cctki2_any_bbox) { \
        const int cctki2_bmin[] = { \
          cctki2_idir<0 ? 0 : cctki2_idir==0 ? cctki2_blo[0] : cctki2_lsh[0] - cctki2_bhi[0], \
          cctki2_jdir<0 ? 0 : cctki2_jdir==0 ? cctki2_blo[1] : cctki2_lsh[1] - cctki2_bhi[1], \
          cctki2_kdir<0 ? 0 : cctki2_kdir==0 ? cctki2_blo[2] : cctki2_lsh[2] - cctki2_bhi[2], \
          cctki2_ldir<0 ? 0 : cctki2_ldir==0 ? cctki2_blo[3] : cctki2_lsh[3] - cctki2_bhi[3], \
        }; \
        const int cctki2_bmax[] = { \
          cctki2_idir<0 ? cctki2_blo[0] : cctki2_idir==0 ? cctki2_lsh[0] - cctki2_bhi[0] : cctki2_lsh[0], \
          cctki2_jdir<0 ? cctki2_blo[1] : cctki2_jdir==0 ? cctki2_lsh[1] - cctki2_bhi[1] : cctki2_lsh[1], \
          cctki2_kdir<0 ? cctki2_blo[2] : cctki2_kdir==0 ? cctki2_lsh[2] - cctki2_bhi[2] : cctki2_lsh[2], \
          cctki2_ldir<0 ? cctki2_blo[3] : cctki2_ldir==0 ? cctki2_lsh[3] - cctki2_bhi[3] : cctki2_lsh[3], \
        }; \
        CCTK_LOOP4STR_NORMAL(name##_boundaries, \
                             i,j,k,l, \
                             ni,nj,nk,nl, \
                             cctki2_idir,cctki2_jdir,cctki2_kdir,cctki2_ldir, \
                             cctki2_bmin[0],cctki2_bmin[1],cctki2_bmin[2],cctki2_bmin[3], \
                             cctki2_bmax[0],cctki2_bmax[1],cctki2_bmax[2],cctki2_bmax[3], \
                             cctki2_cctkGH->cctk_ash[0], \
                             cctki2_cctkGH->cctk_ash[1], \
                             cctki2_cctkGH->cctk_ash[2], \
                             cctki2_cctkGH->cctk_ash[3], \
                             imin,imax, cctki2_istr) { \

#define CCTK_ENDLOOP4STR_BOUNDARIES(name) \
        } CCTK_ENDLOOP4STR_NORMAL(name##_boundaries); \
      } /* if bbox */ \
    } /* for dir */ \
    } /* for dir */ \
    } /* for dir */ \
    } /* for dir */ \
    typedef cctki2_loop4_boundaries_##name cctki2_ensure_proper_nesting CCTK_ATTRIBUTE_UNUSED; \
  } while (0) \



/* LOOP_INTBOUNDARIES */

#define CCTK_LOOP4_INTBOUNDARIES(name, cctki2_cctkGH_, \
                                 i,j,k,l, \
                                 ni,nj,nk,nl, \
                                 cctki2_iblo_,cctki2_jblo_,cctki2_kblo_,cctki2_lblo_, \
                                 cctki2_ibhi_,cctki2_jbhi_,cctki2_kbhi_,cctki2_lbhi_, \
                                 cctki2_ibboxlo_,cctki2_jbboxlo_,cctki2_kbboxlo_,cctki2_lbboxlo_, \
                                 cctki2_ibboxhi_,cctki2_jbboxhi_,cctki2_kbboxhi_,cctki2_lbboxhi_) \
  CCTK_LOOP4STR_INTBOUNDARIES(name, (cctki2_cctkGH_), \
                              i,j,k,l, \
                              ni,nj,nk,nl, \
                              (cctki2_iblo_),(cctki2_jblo_),(cctki2_kblo_),(cctki2_lblo_), \
                              (cctki2_ibhi_),(cctki2_jbhi_),(cctki2_kbhi_),(cctki2_lbhi_), \
                              (cctki2_ibboxlo_),(cctki2_jbboxlo_),(cctki2_kbboxlo_),(cctki2_lbboxlo_), \
                              (cctki2_ibboxhi_),(cctki2_jbboxhi_),(cctki2_kbboxhi_),(cctki2_lbboxhi_), \
                              cctki2_dummy_imin,cctki2_dummy_imax, 1) \

#define CCTK_ENDLOOP4_INTBOUNDARIES(name) \
  CCTK_ENDLOOP4STR_INTBOUNDARIES(name) \

#define CCTK_LOOP4STR_INTBOUNDARIES(name, cctki2_cctkGH_, \
                                    i,j,k,l, \
                                    ni,nj,nk,nl, \
                                    cctki2_iblo_,cctki2_jblo_,cctki2_kblo_,cctki2_lblo_, \
                                    cctki2_ibhi_,cctki2_jbhi_,cctki2_kbhi_,cctki2_lbhi_, \
                                    cctki2_ibboxlo_,cctki2_jbboxlo_,cctki2_kbboxlo_,cctki2_lbboxlo_, \
                                    cctki2_ibboxhi_,cctki2_jbboxhi_,cctki2_kbboxhi_,cctki2_lbboxhi_, \
                                    imin,imax, cctki2_istr_) \
  do { \
    typedef int cctki2_loop4_intboundaries_##name; \
    cGH const *restrict const cctki2_cctkGH = (cctki2_cctkGH_); \
    if (cctki2_cctkGH->cctk_dim != 4) { \
      _Pragma("omp critical") \
      CCTK_WARN(CCTK_WARN_ABORT, \
                "The macro CCTK_LOOP4_INTBOUNDARIES can only be used in 4 dimensions"); \
    } \
    const int cctki2_blo[] = { (cctki2_iblo_),(cctki2_jblo_),(cctki2_kblo_),(cctki2_lblo_) }; \
    const int cctki2_bhi[] = { (cctki2_ibhi_),(cctki2_jbhi_),(cctki2_kbhi_),(cctki2_lbhi_) }; \
    const int cctki2_bbox[] = { (cctki2_ibboxlo_), (cctki2_ibboxhi_),(cctki2_jbboxlo_), (cctki2_jbboxhi_),(cctki2_kbboxlo_), (cctki2_kbboxhi_),(cctki2_lbboxlo_), (cctki2_lbboxhi_) }; \
    const int cctki2_lsh[] = { cctki2_cctkGH->cctk_lsh[0],cctki2_cctkGH->cctk_lsh[1],cctki2_cctkGH->cctk_lsh[2],cctki2_cctkGH->cctk_lsh[3] }; \
    const int cctki2_istr CCTK_ATTRIBUTE_UNUSED = (cctki2_istr_); \
    for (int cctki2_ldir=-1; cctki2_ldir<=+1; ++cctki2_ldir) { \
    for (int cctki2_kdir=-1; cctki2_kdir<=+1; ++cctki2_kdir) { \
    for (int cctki2_jdir=-1; cctki2_jdir<=+1; ++cctki2_jdir) { \
    for (int cctki2_idir=-1; cctki2_idir<=+1; ++cctki2_idir) { \
      const int cctki2_any_bbox = \
        (cctki2_idir<0 ? cctki2_bbox[0] : 0) || (cctki2_idir>0 ? cctki2_bbox[1] : 0) || \
        (cctki2_jdir<0 ? cctki2_bbox[2] : 0) || (cctki2_jdir>0 ? cctki2_bbox[3] : 0) || \
        (cctki2_kdir<0 ? cctki2_bbox[4] : 0) || (cctki2_kdir>0 ? cctki2_bbox[5] : 0) || \
        (cctki2_ldir<0 ? cctki2_bbox[6] : 0) || (cctki2_ldir>0 ? cctki2_bbox[7] : 0); \
      const int cctki2_all_bbox = \
        (cctki2_idir<0 ? cctki2_bbox[0] : 1) && (cctki2_idir>0 ? cctki2_bbox[1] : 1) && \
        (cctki2_jdir<0 ? cctki2_bbox[2] : 1) && (cctki2_jdir>0 ? cctki2_bbox[3] : 1) && \
        (cctki2_kdir<0 ? cctki2_bbox[4] : 1) && (cctki2_kdir>0 ? cctki2_bbox[5] : 1) && \
        (cctki2_ldir<0 ? cctki2_bbox[6] : 1) && (cctki2_ldir>0 ? cctki2_bbox[7] : 1); \
      if (cctki2_all_bbox && cctki2_any_bbox) { \
        const int cctki2_bmin[] = { \
          cctki2_idir<0 ? 0 : cctki2_idir==0 ? cctki2_blo[0] : cctki2_lsh[0] - cctki2_bhi[0], \
          cctki2_jdir<0 ? 0 : cctki2_jdir==0 ? cctki2_blo[1] : cctki2_lsh[1] - cctki2_bhi[1], \
          cctki2_kdir<0 ? 0 : cctki2_kdir==0 ? cctki2_blo[2] : cctki2_lsh[2] - cctki2_bhi[2], \
          cctki2_ldir<0 ? 0 : cctki2_ldir==0 ? cctki2_blo[3] : cctki2_lsh[3] - cctki2_bhi[3], \
        }; \
        const int cctki2_bmax[] = { \
          cctki2_idir<0 ? cctki2_blo[0] : cctki2_idir==0 ? cctki2_lsh[0] - cctki2_bhi[0] : cctki2_lsh[0], \
          cctki2_jdir<0 ? cctki2_blo[1] : cctki2_jdir==0 ? cctki2_lsh[1] - cctki2_bhi[1] : cctki2_lsh[1], \
          cctki2_kdir<0 ? cctki2_blo[2] : cctki2_kdir==0 ? cctki2_lsh[2] - cctki2_bhi[2] : cctki2_lsh[2], \
          cctki2_ldir<0 ? cctki2_blo[3] : cctki2_ldir==0 ? cctki2_lsh[3] - cctki2_bhi[3] : cctki2_lsh[3], \
        }; \
        CCTK_LOOP4STR_NORMAL(name##_intboundaries, \
                             i,j,k,l, \
                             ni,nj,nk,nl, \
                             cctki2_idir,cctki2_jdir,cctki2_kdir,cctki2_ldir, \
                             cctki2_bmin[0],cctki2_bmin[1],cctki2_bmin[2],cctki2_bmin[3], \
                             cctki2_bmax[0],cctki2_bmax[1],cctki2_bmax[2],cctki2_bmax[3], \
                             cctki2_cctkGH->cctk_ash[0], \
                             cctki2_cctkGH->cctk_ash[1], \
                             cctki2_cctkGH->cctk_ash[2], \
                             cctki2_cctkGH->cctk_ash[3], \
                             imin,imax, cctki2_istr) { \

#define CCTK_ENDLOOP4STR_INTBOUNDARIES(name) \
        } CCTK_ENDLOOP4STR_NORMAL(name##_intboundaries); \
      } /* if bbox */ \
    } /* for dir */ \
    } /* for dir */ \
    } /* for dir */ \
    } /* for dir */ \
    typedef cctki2_loop4_intboundaries_##name cctki2_ensure_proper_nesting CCTK_ATTRIBUTE_UNUSED; \
  } while (0) \



/* LOOP_ALL */

#define CCTK_LOOP4_ALL(name, cctki3_cctkGH_, \
                       i,j,k,l) \
  CCTK_LOOP4STR_ALL(name, (cctki3_cctkGH_), \
                    i,j,k,l, \
                    cctki3_dummy_imin,cctki3_dummy_imax, 1) \

#define CCTK_ENDLOOP4_ALL(name) \
  CCTK_ENDLOOP4STR_ALL(name) \

#define CCTK_LOOP4STR_ALL(name, cctki3_cctkGH_, \
                          i,j,k,l, \
                          imin,imax, cctki3_istr_) \
  do { \
    typedef int cctki3_loop4_all_##name; \
    cGH const *restrict const cctki3_cctkGH = (cctki3_cctkGH_); \
    if (cctki3_cctkGH->cctk_dim != 4) { \
      _Pragma("omp critical") \
      CCTK_WARN(CCTK_WARN_ABORT, \
                "The macro CCTK_LOOP4_ALL can only be used in 4 dimensions"); \
    } \
    CCTK_LOOP4STR(name##_all, \
                  i,j,k,l, \
                  0,0,0,0, \
                  cctki3_cctkGH->cctk_lsh[0], \
                  cctki3_cctkGH->cctk_lsh[1], \
                  cctki3_cctkGH->cctk_lsh[2], \
                  cctki3_cctkGH->cctk_lsh[3], \
                  cctki3_cctkGH->cctk_ash[0], \
                  cctki3_cctkGH->cctk_ash[1], \
                  cctki3_cctkGH->cctk_ash[2], \
                  cctki3_cctkGH->cctk_ash[3], \
                  imin,imax, (cctki3_istr_)) { \

#define CCTK_ENDLOOP4STR_ALL(name) \
    } CCTK_ENDLOOP4STR(name##_all); \
    typedef cctki3_loop4_all_##name cctki3_ensure_proper_nesting CCTK_ATTRIBUTE_UNUSED; \
  } while (0) \



/* LOOP_INT */

#define CCTK_LOOP4_INT(name, cctki3_cctkGH_, \
                       i,j,k,l) \
  CCTK_LOOP4STR_INT(name, (cctki3_cctkGH_), \
                    i,j,k,l, \
                    cctki3_dummy_imin,cctki3_dummy_imax, 1) \

#define CCTK_ENDLOOP4_INT(name) \
  CCTK_ENDLOOP4STR_INT(name) \

#define CCTK_LOOP4STR_INT(name, cctki3_cctkGH_, \
                          i,j,k,l, \
                          imin,imax, cctki3_istr_) \
  do { \
    typedef int cctki3_loop4_int_##name; \
    cGH const *restrict const cctki3_cctkGH = (cctki3_cctkGH_); \
    if (cctki3_cctkGH->cctk_dim != 4) { \
      _Pragma("omp critical") \
      CCTK_WARN(CCTK_WARN_ABORT, \
                "The macro CCTK_LOOP4_INT can only be used in 4 dimensions"); \
    } \
    CCTK_INT cctki3_bndsize    [8]; \
    CCTK_INT cctki3_is_ghostbnd[8]; \
    CCTK_INT cctki3_is_symbnd  [8]; \
    CCTK_INT cctki3_is_physbnd [8]; \
    _Pragma("omp single copyprivate(cctki3_bndsize)") \
    GetBoundarySizesAndTypes \
      (cctki3_cctkGH, 8, cctki3_bndsize, cctki3_is_ghostbnd, cctki3_is_symbnd, cctki3_is_physbnd); \
    CCTK_LOOP4STR_INTERIOR(name##_int, \
                           cctki3_cctkGH, \
                           i,j,k,l, \
                           cctki3_bndsize[0],cctki3_bndsize[2],cctki3_bndsize[4],cctki3_bndsize[6], \
                           cctki3_bndsize[1],cctki3_bndsize[3],cctki3_bndsize[5],cctki3_bndsize[7], \
                           imin,imax, (cctki3_istr_)) { \

#define CCTK_ENDLOOP4STR_INT(name) \
    } CCTK_ENDLOOP4STR_INTERIOR(name##_int); \
    typedef cctki3_loop4_int_##name cctki3_ensure_proper_nesting CCTK_ATTRIBUTE_UNUSED; \
  } while (0) \



/* LOOP_BND */

#define CCTK_LOOP4_BND(name, cctki3_cctkGH_, \
                       i,j,k,l, \
                       ni,nj,nk,nl) \
  CCTK_LOOP4STR_BND(name, (cctki3_cctkGH_), \
                    i,j,k,l, \
                    ni,nj,nk,nl, \
                    cctki3_dummy_imin,cctki3_dummy_imax, 1) \

#define CCTK_ENDLOOP4_BND(name) \
  CCTK_ENDLOOP4STR_BND(name) \

#define CCTK_LOOP4STR_BND(name, cctki3_cctkGH_, \
                          i,j,k,l, \
                          ni,nj,nk,nl, \
                          imin,imax, cctki3_istr_) \
  do { \
    typedef int cctki3_loop4_bnd_##name; \
    cGH const *restrict const cctki3_cctkGH = (cctki3_cctkGH_); \
    if (cctki3_cctkGH->cctk_dim != 4) { \
      _Pragma("omp critical") \
      CCTK_WARN(CCTK_WARN_ABORT, \
                "The macro CCTK_LOOP4_BND can only be used in 4 dimensions"); \
    } \
    CCTK_INT cctki3_bndsize    [8]; \
    CCTK_INT cctki3_is_ghostbnd[8]; \
    CCTK_INT cctki3_is_symbnd  [8]; \
    CCTK_INT cctki3_is_physbnd [8]; \
    _Pragma("omp single copyprivate(cctki3_bndsize, cctki3_is_physbnd)") \
    GetBoundarySizesAndTypes \
      (cctki3_cctkGH, 8, cctki3_bndsize, cctki3_is_ghostbnd, cctki3_is_symbnd, cctki3_is_physbnd); \
    CCTK_LOOP4STR_BOUNDARIES(name##_bnd, \
                             cctki3_cctkGH, \
                             i,j,k,l, \
                             ni,nj,nk,nl, \
                             cctki3_bndsize[0],cctki3_bndsize[2],cctki3_bndsize[4],cctki3_bndsize[6], \
                             cctki3_bndsize[1],cctki3_bndsize[3],cctki3_bndsize[5],cctki3_bndsize[7], \
                             cctki3_is_physbnd[0],cctki3_is_physbnd[2],cctki3_is_physbnd[4],cctki3_is_physbnd[6], \
                             cctki3_is_physbnd[1],cctki3_is_physbnd[3],cctki3_is_physbnd[5],cctki3_is_physbnd[7], \
                             imin,imax, (cctki3_istr_)) { \

#define CCTK_ENDLOOP4STR_BND(name) \
    } CCTK_ENDLOOP4STR_BOUNDARIES(name##_bnd); \
    typedef cctki3_loop4_bnd_##name cctki3_ensure_proper_nesting CCTK_ATTRIBUTE_UNUSED; \
  } while (0) \



/* LOOP_INTBND */

#define CCTK_LOOP4_INTBND(name, cctki3_cctkGH_, \
                           i,j,k,l, \
                           ni,nj,nk,nl) \
  CCTK_LOOP4STR_INTBND(name, (cctki3_cctkGH_), \
                        i,j,k,l, \
                        ni,nj,nk,nl, \
                        cctki3_dummy_imin,cctki3_dummy_imax, 1) \

#define CCTK_ENDLOOP4_INTBND(name) \
  CCTK_ENDLOOP4STR_INTBND(name) \

#define CCTK_LOOP4STR_INTBND(name, cctki3_cctkGH_, \
                              i,j,k,l, \
                              ni,nj,nk,nl, \
                              imin,imax, cctki3_istr_) \
  do { \
    typedef int cctki3_loop4_intbnd_##name; \
    cGH const *restrict const cctki3_cctkGH = (cctki3_cctkGH_); \
    if (cctki3_cctkGH->cctk_dim != 4) { \
      _Pragma("omp critical") \
      CCTK_WARN(CCTK_WARN_ABORT, \
                "The macro CCTK_LOOP4_INTBND can only be used in 4 dimensions"); \
    } \
    CCTK_INT cctki3_bndsize    [8]; \
    CCTK_INT cctki3_is_ghostbnd[8]; \
    CCTK_INT cctki3_is_symbnd  [8]; \
    CCTK_INT cctki3_is_physbnd [8]; \
    _Pragma("omp single copyprivate(cctki3_bndsize, cctki3_is_physbnd)") \
    GetBoundarySizesAndTypes \
      (cctki3_cctkGH, 8, cctki3_bndsize, cctki3_is_ghostbnd, cctki3_is_symbnd, cctki3_is_physbnd); \
    CCTK_LOOP4STR_INTBOUNDARIES(name##_intbnd, \
                                cctki3_cctkGH, \
                                i,j,k,l, \
                                ni,nj,nk,nl, \
                                cctki3_bndsize[0],cctki3_bndsize[2],cctki3_bndsize[4],cctki3_bndsize[6], \
                                cctki3_bndsize[1],cctki3_bndsize[3],cctki3_bndsize[5],cctki3_bndsize[7], \
                                cctki3_is_physbnd[0],cctki3_is_physbnd[2],cctki3_is_physbnd[4],cctki3_is_physbnd[6], \
                                cctki3_is_physbnd[1],cctki3_is_physbnd[3],cctki3_is_physbnd[5],cctki3_is_physbnd[7], \
                                imin,imax, (cctki3_istr_)) { \

#define CCTK_ENDLOOP4STR_INTBND(name) \
    } CCTK_ENDLOOP4STR_INTBOUNDARIES(name##_intbnd); \
    typedef cctki3_loop4_intbnd_##name cctki3_ensure_proper_nesting CCTK_ATTRIBUTE_UNUSED; \
  } while (0) \

#endif /* #ifdef CCODE */



#ifdef FCODE

/* LOOP */

#define CCTK_LOOP4_NORMAL_DECLARE(name) \
   CCTK_LOOP4STR_NORMAL_DECLARE(name) \
   && integer :: name/**/0_dummy_imin, name/**/0_dummy_imax \

#define CCTK_LOOP4_NORMAL_OMP_PRIVATE(name) \
   CCTK_LOOP4STR_NORMAL_OMP_PRIVATE(name) \

#define CCTK_LOOP4_NORMAL(name, \
                          i,j,k,l, \
                          ni,nj,nk,nl, \
                          cctki0_idir,cctki0_jdir,cctki0_kdir,cctki0_ldir, \
                          cctki0_imin,cctki0_jmin,cctki0_kmin,cctki0_lmin, \
                          cctki0_imax,cctki0_jmax,cctki0_kmax,cctki0_lmax, \
                          cctki0_iash,cctki0_jash,cctki0_kash,cctki0_lash) \
   CCTK_LOOP4STR_NORMAL(name, \
                        i,j,k,l, \
                        ni,nj,nk,nl, \
                        cctki0_idir,cctki0_jdir,cctki0_kdir,cctki0_ldir, \
                        cctki0_imin,cctki0_jmin,cctki0_kmin,cctki0_lmin, \
                        cctki0_imax,cctki0_jmax,cctki0_kmax,cctki0_lmax, \
                        cctki0_iash,cctki0_jash,cctki0_kash,cctki0_lash, \
                        name/**/0_dummy_imin,name/**/0_dummy_imax, 1) \

#define CCTK_ENDLOOP4_NORMAL(name) \
   CCTK_ENDLOOP4STR_NORMAL(name) \

#define CCTK_LOOP4STR_NORMAL_DECLARE(name) \
   && integer :: name/**/0_idir,name/**/0_jdir,name/**/0_kdir,name/**/0_ldir \
   && integer :: name/**/0_imin,name/**/0_jmin,name/**/0_kmin,name/**/0_lmin \
   && integer :: name/**/0_imax,name/**/0_jmax,name/**/0_kmax,name/**/0_lmax \
   && integer :: name/**/0_iash,name/**/0_jash,name/**/0_kash,name/**/0_lash \
   && integer :: name/**/0_istr \

#define CCTK_LOOP4STR_NORMAL_OMP_PRIVATE(name) \
   && !$omp private (i,j,k,l) \
   && !$omp private (ni,nj,nk,nl) \

#define CCTK_LOOP4STR_NORMAL(name, \
                             i,j,k,l, \
                             ni,nj,nk,nl, \
                             cctki0_idir,cctki0_jdir,cctki0_kdir,cctki0_ldir, \
                             cctki0_imin,cctki0_jmin,cctki0_kmin,cctki0_lmin, \
                             cctki0_imax,cctki0_jmax,cctki0_kmax,cctki0_lmax, \
                             cctki0_iash,cctki0_jash,cctki0_kash,cctki0_lash, \
                             imin,imax, cctki0_istr) \
   && name/**/0_idir = cctki0_idir \
   && name/**/0_jdir = cctki0_jdir \
   && name/**/0_kdir = cctki0_kdir \
   && name/**/0_ldir = cctki0_ldir \
   && name/**/0_imin = cctki0_imin \
   && name/**/0_jmin = cctki0_jmin \
   && name/**/0_kmin = cctki0_kmin \
   && name/**/0_lmin = cctki0_lmin \
   && name/**/0_imax = cctki0_imax \
   && name/**/0_jmax = cctki0_jmax \
   && name/**/0_kmax = cctki0_kmax \
   && name/**/0_lmax = cctki0_lmax \
   && name/**/0_iash = cctki0_iash \
   && name/**/0_jash = cctki0_jash \
   && name/**/0_kash = cctki0_kash \
   && name/**/0_lash = cctki0_lash \
   && name/**/0_istr = cctki0_istr \
   && imin = name/**/0_imin \
   && imax = name/**/0_imax \
   && !$omp do collapse(3) \
   && do l = name/**/0_lmin, name/**/0_lmax \
   && do k = name/**/0_kmin, name/**/0_kmax \
   && do j = name/**/0_jmin, name/**/0_jmax \
   && do i = name/**/0_imin - modulo((imin+name/**/0_iash*(j+name/**/0_jash*(k+name/**/0_kash*(l)))), name/**/0_istr), name/**/0_imax, name/**/0_istr \
   &&    ni = 0 \
   &&    nj = 0 \
   &&    nk = 0 \
   &&    nl = 0 \
   &&    if (name/**/0_idir < 0) ni = i \
   &&    if (name/**/0_jdir < 0) nj = j \
   &&    if (name/**/0_kdir < 0) nk = k \
   &&    if (name/**/0_ldir < 0) nl = l \
   &&    if (name/**/0_idir > 0) ni = name/**/0_imax+1-i \
   &&    if (name/**/0_jdir > 0) nj = name/**/0_jmax+1-j \
   &&    if (name/**/0_kdir > 0) nk = name/**/0_kmax+1-k \
   &&    if (name/**/0_ldir > 0) nl = name/**/0_lmax+1-l \

#define CCTK_ENDLOOP4STR_NORMAL(name) \
   && end do \
   && end do \
   && end do \
   && end do \



#define CCTK_LOOP4_DECLARE(name) \
   CCTK_LOOP4STR_DECLARE(name) \
   && integer :: name/**/1_dummy_imin, name/**/1_dummy_imax \

#define CCTK_LOOP4_OMP_PRIVATE(name) \
   CCTK_LOOP4STR_OMP_PRIVATE(name) \

#define CCTK_LOOP4(name, \
                   i,j,k,l, \
                   cctki0_imin,cctki0_jmin,cctki0_kmin,cctki0_lmin, \
                   cctki0_imax,cctki0_jmax,cctki0_kmax,cctki0_lmax, \
                   cctki0_iash,cctki0_jash,cctki0_kash,cctki0_lash) \
   CCTK_LOOP4STR(name, \
                 i,j,k,l, \
                 cctki0_imin,cctki0_jmin,cctki0_kmin,cctki0_lmin, \
                 cctki0_imax,cctki0_jmax,cctki0_kmax,cctki0_lmax, \
                 cctki0_iash,cctki0_jash,cctki0_kash,cctki0_lash, \
                 name/**/1_dummy_imin,name/**/1_dummy_imax, 1) \

#define CCTK_ENDLOOP4(name) \
   CCTK_ENDLOOP4STR(name) \

#define CCTK_LOOP4STR_DECLARE(name) \
   CCTK_LOOP4STR_NORMAL_DECLARE(name) \
   && integer :: name/**/1_ni,name/**/1_nj,name/**/1_nk,name/**/1_nl \

#define CCTK_LOOP4STR_OMP_PRIVATE(name) \
   CCTK_LOOP4STR_NORMAL_OMP_PRIVATE(name) \

#define CCTK_LOOP4STR(name, \
                      i,j,k,l, \
                      cctki1_imin,cctki1_jmin,cctki1_kmin,cctki1_lmin, \
                      cctki1_imax,cctki1_jmax,cctki1_kmax,cctki1_lmax, \
                      cctki1_iash,cctki1_jash,cctki1_kash,cctki1_lash, \
                      imin,imax, cctki1_istr) \
   CCTK_LOOP4STR_NORMAL(name, \
                        i,j,k,l, \
                        name/**/1_ni,name/**/1_nj,name/**/1_nk,name/**/1_nl, \
                        0,0,0,0, \
                        cctki1_imin,cctki1_jmin,cctki1_kmin,cctki1_lmin, \
                        cctki1_imax,cctki1_jmax,cctki1_kmax,cctki1_lmax, \
                        cctki1_iash,cctki1_jash,cctki1_kash,cctki1_lash, \
                        imin,imax, cctki1_istr) \

#define CCTK_ENDLOOP4STR(name) \
   CCTK_ENDLOOP4STR_NORMAL(name) \



/* LOOP_INTERIOR */

#define CCTK_LOOP4_INTERIOR_DECLARE(name) \
   CCTK_LOOP4STR_INTERIOR_DECLARE(name) \
   && integer :: name/**/2_dummy_imin, name/**/2_dummy_imax \

#define CCTK_LOOP4_INTERIOR_OMP_PRIVATE(name) \
   CCTK_LOOP4STR_INTERIOR_OMP_PRIVATE(name) \

#define CCTK_LOOP4_INTERIOR(name, \
                            i,j,k,l, \
                            cctki2_iblo,cctki2_jblo,cctki2_kblo,cctki2_lblo, \
                            cctki2_ibhi,cctki2_jbhi,cctki2_kbhi,cctki2_lbhi) \
   CCTK_LOOP4STR_INTERIOR(name, \
                          i,j,k,l, \
                          cctki2_iblo,cctki2_jblo,cctki2_kblo,cctki2_lblo, \
                          cctki2_ibhi,cctki2_jbhi,cctki2_kbhi,cctki2_lbhi, \
                          name/**/2_dummy_imin,name/**/2_dummy_imax, 1) \

#define CCTK_ENDLOOP4_INTERIOR(name) \
   CCTK_ENDLOOP4STR_INTERIOR(name) \

#define CCTK_LOOP4STR_INTERIOR_DECLARE(name) \
   CCTK_LOOP4STR_DECLARE(name/**/_interior) \

#define CCTK_LOOP4STR_INTERIOR_OMP_PRIVATE(name) \
   CCTK_LOOP4STR_OMP_PRIVATE(name/**/_interior) \

#define CCTK_LOOP4STR_INTERIOR(name, \
                               i,j,k,l, \
                               cctki2_iblo,cctki2_jblo,cctki2_kblo,cctki2_lblo, \
                               cctki2_ibhi,cctki2_jbhi,cctki2_kbhi,cctki2_lbhi, \
                               imin,imax, cctki2_istr) \
   CCTK_LOOP4STR(name/**/_interior, \
                 i,j,k,l, \
                 (cctki2_iblo)+1, \
                 (cctki2_jblo)+1, \
                 (cctki2_kblo)+1, \
                 (cctki2_lblo)+1, \
                 cctk_lsh(1)-(cctki2_ibhi), \
                 cctk_lsh(2)-(cctki2_jbhi), \
                 cctk_lsh(3)-(cctki2_kbhi), \
                 cctk_lsh(4)-(cctki2_lbhi), \
                 cctk_ash(1),cctk_ash(2),cctk_ash(3),cctk_ash(4), \
                 imin,imax, cctki2_istr) \

#define CCTK_ENDLOOP4STR_INTERIOR(name) \
   CCTK_ENDLOOP4STR(name/**/_interior) \



/* LOOP_BOUNDARIES */

#define CCTK_LOOP4_BOUNDARIES_DECLARE(name) \
   CCTK_LOOP4STR_BOUNDARIES_DECLARE(name) \
   && integer :: name/**/2_dummy_imin, name/**/2_dummy_imax \

#define CCTK_LOOP4_BOUNDARIES_OMP_PRIVATE(name) \
   CCTK_LOOP4STR_BOUNDARIES_OMP_PRIVATE(name) \

#define CCTK_LOOP4_BOUNDARIES(name, \
                              i,j,k,l, \
                              ni,nj,nk,nl, \
                              cctki2_iblo,cctki2_jblo,cctki2_kblo,cctki2_lblo, \
                              cctki2_ibhi,cctki2_jbhi,cctki2_kbhi,cctki2_lbhi, \
                              cctki2_ibboxlo,cctki2_jbboxlo,cctki2_kbboxlo,cctki2_lbboxlo, \
                              cctki2_ibboxhi,cctki2_jbboxhi,cctki2_kbboxhi,cctki2_lbboxhi) \
   CCTK_LOOP4STR_BOUNDARIES(name, \
                            i,j,k,l, \
                            ni,nj,nk,nl, \
                            cctki2_iblo,cctki2_jblo,cctki2_kblo,cctki2_lblo, \
                            cctki2_ibhi,cctki2_jbhi,cctki2_kbhi,cctki2_lbhi, \
                            cctki2_ibboxlo,cctki2_jbboxlo,cctki2_kbboxlo,cctki2_lbboxlo, \
                            cctki2_ibboxhi,cctki2_jbboxhi,cctki2_kbboxhi,cctki2_lbboxhi, \
                            name/**/2_dummy_imin,name/**/2_dummy_imax, 1) \

#define CCTK_ENDLOOP4_BOUNDARIES(name) \
   CCTK_ENDLOOP4STR_BOUNDARIES(name) \

#define CCTK_LOOP4STR_BOUNDARIES_DECLARE(name) \
   CCTK_LOOP4STR_NORMAL_DECLARE(name/**/_boundaries) \
   && integer :: name/**/2_blo(4), name/**/2_bhi(4) \
   && integer :: name/**/2_bboxlo(4), name/**/2_bboxhi(4) \
   && integer :: name/**/2_istr \
   && integer :: name/**/2_idir \
   && integer :: name/**/2_jdir \
   && integer :: name/**/2_kdir \
   && integer :: name/**/2_ldir \
   && logical :: name/**/2_any_bbox \
   && integer :: name/**/2_bmin(4), name/**/2_bmax(4) \

#define CCTK_LOOP4STR_BOUNDARIES_OMP_PRIVATE(name) \
   CCTK_LOOP4STR_NORMAL_OMP_PRIVATE(name/**/_boundaries) \
   && !$omp private (name/**/2_bmin, name/**/2_bmax) \

#define CCTK_LOOP4STR_BOUNDARIES(name, \
                                 i,j,k,l, \
                                 ni,nj,nk,nl, \
                                 cctki2_iblo,cctki2_jblo,cctki2_kblo,cctki2_lblo, \
                                 cctki2_ibhi,cctki2_jbhi,cctki2_kbhi,cctki2_lbhi, \
                                 cctki2_ibboxlo,cctki2_jbboxlo,cctki2_kbboxlo,cctki2_lbboxlo, \
                                 cctki2_ibboxhi,cctki2_jbboxhi,cctki2_kbboxhi,cctki2_lbboxhi, \
                                 imin,imax, cctki2_istr) \
   && name/**/2_blo = (/ cctki2_iblo,cctki2_jblo,cctki2_kblo,cctki2_lblo /) \
   && name/**/2_bhi = (/ cctki2_ibhi,cctki2_jbhi,cctki2_kbhi,cctki2_lbhi /) \
   && name/**/2_bboxlo = (/ cctki2_ibboxlo,cctki2_jbboxlo,cctki2_kbboxlo,cctki2_lbboxlo /) \
   && name/**/2_bboxhi = (/ cctki2_ibboxhi,cctki2_jbboxhi,cctki2_kbboxhi,cctki2_lbboxhi /) \
   && name/**/2_istr = (cctki2_istr) \
   && do name/**/2_ldir=-1, +1 \
   && do name/**/2_kdir=-1, +1 \
   && do name/**/2_jdir=-1, +1 \
   && do name/**/2_idir=-1, +1 \
   &&     name/**/2_any_bbox = .false. \
   &&     if (name/**/2_idir<0) name/**/2_any_bbox = name/**/2_any_bbox .or. name/**/2_bboxlo(1) /= 0 \
   &&     if (name/**/2_jdir<0) name/**/2_any_bbox = name/**/2_any_bbox .or. name/**/2_bboxlo(2) /= 0 \
   &&     if (name/**/2_kdir<0) name/**/2_any_bbox = name/**/2_any_bbox .or. name/**/2_bboxlo(3) /= 0 \
   &&     if (name/**/2_ldir<0) name/**/2_any_bbox = name/**/2_any_bbox .or. name/**/2_bboxlo(4) /= 0 \
   &&     if (name/**/2_idir>0) name/**/2_any_bbox = name/**/2_any_bbox .or. name/**/2_bboxhi(1) /= 0 \
   &&     if (name/**/2_jdir>0) name/**/2_any_bbox = name/**/2_any_bbox .or. name/**/2_bboxhi(2) /= 0 \
   &&     if (name/**/2_kdir>0) name/**/2_any_bbox = name/**/2_any_bbox .or. name/**/2_bboxhi(3) /= 0 \
   &&     if (name/**/2_ldir>0) name/**/2_any_bbox = name/**/2_any_bbox .or. name/**/2_bboxhi(4) /= 0 \
   &&    if (name/**/2_any_bbox) then \
   &&       name/**/2_bmin(1) = name/**/2_blo(1)+1 \
   &&       name/**/2_bmin(2) = name/**/2_blo(2)+1 \
   &&       name/**/2_bmin(3) = name/**/2_blo(3)+1 \
   &&       name/**/2_bmin(4) = name/**/2_blo(4)+1 \
   &&       if (name/**/2_idir<0) name/**/2_bmin(1) = 1 \
   &&       if (name/**/2_jdir<0) name/**/2_bmin(2) = 1 \
   &&       if (name/**/2_kdir<0) name/**/2_bmin(3) = 1 \
   &&       if (name/**/2_ldir<0) name/**/2_bmin(4) = 1 \
   &&       if (name/**/2_idir>0) name/**/2_bmin(1) = cctk_lsh(1) - name/**/2_bhi(1) \
   &&       if (name/**/2_jdir>0) name/**/2_bmin(2) = cctk_lsh(2) - name/**/2_bhi(2) \
   &&       if (name/**/2_kdir>0) name/**/2_bmin(3) = cctk_lsh(3) - name/**/2_bhi(3) \
   &&       if (name/**/2_ldir>0) name/**/2_bmin(4) = cctk_lsh(4) - name/**/2_bhi(4) \
   &&       name/**/2_bmax(1) = cctk_lsh(1) - name/**/2_bhi(1) \
   &&       name/**/2_bmax(2) = cctk_lsh(2) - name/**/2_bhi(2) \
   &&       name/**/2_bmax(3) = cctk_lsh(3) - name/**/2_bhi(3) \
   &&       name/**/2_bmax(4) = cctk_lsh(4) - name/**/2_bhi(4) \
   &&       if (name/**/2_idir<0) name/**/2_bmax(1) = name/**/2_blo(1) \
   &&       if (name/**/2_jdir<0) name/**/2_bmax(2) = name/**/2_blo(2) \
   &&       if (name/**/2_kdir<0) name/**/2_bmax(3) = name/**/2_blo(3) \
   &&       if (name/**/2_ldir<0) name/**/2_bmax(4) = name/**/2_blo(4) \
   &&       if (name/**/2_idir>0) name/**/2_bmax(1) = cctk_lsh(1) \
   &&       if (name/**/2_jdir>0) name/**/2_bmax(2) = cctk_lsh(2) \
   &&       if (name/**/2_kdir>0) name/**/2_bmax(3) = cctk_lsh(3) \
   &&       if (name/**/2_ldir>0) name/**/2_bmax(4) = cctk_lsh(4) \
   &&       CCTK_LOOP4STR_NORMAL(name/**/_boundaries, \
                                 i,j,k,l, \
                                 ni,nj,nk,nl, \
                                 name/**/2_idir,name/**/2_jdir,name/**/2_kdir,name/**/2_ldir, \
                                 name/**/2_bmin(1),name/**/2_bmin(2),name/**/2_bmin(3),name/**/2_bmin(4), \
                                 name/**/2_bmax(1),name/**/2_bmax(2),name/**/2_bmax(3),name/**/2_bmax(4), \
                                 cctk_ash(1), \
                                 cctk_ash(2), \
                                 cctk_ash(3), \
                                 cctk_ash(4), \
                                 imin,imax, name/**/2_istr) \

#define CCTK_ENDLOOP4STR_BOUNDARIES(name) \
            CCTK_ENDLOOP4STR_NORMAL(name/**/_boundaries) \
   &&    end if /* bbox */ \
   && end do /* dir */ \
   && end do /* dir */ \
   && end do /* dir */ \
   && end do /* dir */ \



/* LOOP_INTBOUNDARIES */

#define CCTK_LOOP4_INTBOUNDARIES_DECLARE(name) \
   CCTK_LOOP4STR_INTBOUNDARIES_DECLARE(name) \
   && integer :: name/**/2_dummy_imin, name/**/2_dummy_imax \

#define CCTK_LOOP4_INTBOUNDARIES_OMP_PRIVATE(name) \
   CCTK_LOOP4STR_INTBOUNDARIES_OMP_PRIVATE(name) \

#define CCTK_LOOP4_INTBOUNDARIES(name, \
                                 i,j,k,l, \
                                 ni,nj,nk,nl, \
                                 cctki2_iblo,cctki2_jblo,cctki2_kblo,cctki2_lblo, \
                                 cctki2_ibhi,cctki2_jbhi,cctki2_kbhi,cctki2_lbhi, \
                                 cctki2_ibboxlo,cctki2_jbboxlo,cctki2_kbboxlo,cctki2_lbboxlo, \
                                 cctki2_ibboxhi,cctki2_jbboxhi,cctki2_kbboxhi,cctki2_lbboxhi) \
   CCTK_LOOP4STR_INTBOUNDARIES(name, \
                               i,j,k,l, \
                               ni,nj,nk,nl, \
                               cctki2_iblo,cctki2_jblo,cctki2_kblo,cctki2_lblo, \
                               cctki2_ibhi,cctki2_jbhi,cctki2_kbhi,cctki2_lbhi, \
                               cctki2_ibboxlo,cctki2_jbboxlo,cctki2_kbboxlo,cctki2_lbboxlo, \
                               cctki2_ibboxhi,cctki2_jbboxhi,cctki2_kbboxhi,cctki2_lbboxhi, \
                               name/**/2_dummy_imin,name/**/2_dummy_max, 1) \

#define CCTK_ENDLOOP4_INTBOUNDARIES(name) \
   CCTK_ENDLOOP4STR_INTBOUNDARIES(name) \

#define CCTK_LOOP4STR_INTBOUNDARIES_DECLARE(name) \
   CCTK_LOOP4STR_NORMAL_DECLARE(name/**/_intboundaries) \
   && integer :: name/**/2_blo(4), name/**/2_bhi(4) \
   && integer :: name/**/2_bboxlo(4), name/**/2_bboxhi(4) \
   && integer :: name/**/2_istr \
   && integer :: name/**/2_idir \
   && integer :: name/**/2_jdir \
   && integer :: name/**/2_kdir \
   && integer :: name/**/2_ldir \
   && logical :: name/**/2_any_bbox, name/**/2_all_bbox \
   && integer :: name/**/2_bmin(4), name/**/2_bmax(4) \

#define CCTK_LOOP4STR_INTBOUNDARIES_OMP_PRIVATE(name) \
   CCTK_LOOP4STR_NORMAL_OMP_PRIVATE(name/**/_intboundaries) \
   && !$omp private (name/**/2_any_bbox, name/**/2_all_bbox) \
   && !$omp private (name/**/2_bmin, name/**/2_bmax \

#define CCTK_LOOP4STR_INTBOUNDARIES(name, \
                                    i,j,k,l, \
                                    ni,nj,nk,nl, \
                                    cctki2_iblo,cctki2_jblo,cctki2_kblo,cctki2_lblo, \
                                    cctki2_ibhi,cctki2_jbhi,cctki2_kbhi,cctki2_lbhi, \
                                    cctki2_ibboxlo,cctki2_jbboxlo,cctki2_kbboxlo,cctki2_lbboxlo, \
                                    cctki2_ibboxhi,cctki2_jbboxhi,cctki2_kbboxhi,cctki2_lbboxhi, \
                                    imin,imax, cctki2_istr) \
   && name/**/2_blo = (/ cctki2_iblo,cctki2_jblo,cctki2_kblo,cctki2_lblo /) \
   && name/**/2_bhi = (/ cctki2_ibhi,cctki2_jbhi,cctki2_kbhi,cctki2_lbhi /) \
   && name/**/2_bboxlo = (/ cctki2_ibboxlo,cctki2_jbboxlo,cctki2_kbboxlo,cctki2_lbboxlo /) \
   && name/**/2_bboxhi = (/ cctki2_ibboxhi,cctki2_jbboxhi,cctki2_kbboxhi,cctki2_lbboxhi /) \
   && name/**/2_istr = (cctki2_istr) \
   && do name/**/2_ldir=-1, +1 \
   && do name/**/2_kdir=-1, +1 \
   && do name/**/2_jdir=-1, +1 \
   && do name/**/2_idir=-1, +1 \
   &&     name/**/2_any_bbox = .false. \
   &&     if (name/**/2_idir<0) name/**/2_any_bbox = name/**/2_any_bbox .or. name/**/2_bboxlo(1) /= 0 \
   &&     if (name/**/2_jdir<0) name/**/2_any_bbox = name/**/2_any_bbox .or. name/**/2_bboxlo(2) /= 0 \
   &&     if (name/**/2_kdir<0) name/**/2_any_bbox = name/**/2_any_bbox .or. name/**/2_bboxlo(3) /= 0 \
   &&     if (name/**/2_ldir<0) name/**/2_any_bbox = name/**/2_any_bbox .or. name/**/2_bboxlo(4) /= 0 \
   &&     if (name/**/2_idir>0) name/**/2_any_bbox = name/**/2_any_bbox .or. name/**/2_bboxhi(1) /= 0 \
   &&     if (name/**/2_jdir>0) name/**/2_any_bbox = name/**/2_any_bbox .or. name/**/2_bboxhi(2) /= 0 \
   &&     if (name/**/2_kdir>0) name/**/2_any_bbox = name/**/2_any_bbox .or. name/**/2_bboxhi(3) /= 0 \
   &&     if (name/**/2_ldir>0) name/**/2_any_bbox = name/**/2_any_bbox .or. name/**/2_bboxhi(4) /= 0 \
   &&     name/**/2_all_bbox = .true. \
   &&     if (name/**/2_idir<0) name/**/2_all_bbox = name/**/2_all_bbox .and. name/**/2_bboxlo(1) /= 0 \
   &&     if (name/**/2_jdir<0) name/**/2_all_bbox = name/**/2_all_bbox .and. name/**/2_bboxlo(2) /= 0 \
   &&     if (name/**/2_kdir<0) name/**/2_all_bbox = name/**/2_all_bbox .and. name/**/2_bboxlo(3) /= 0 \
   &&     if (name/**/2_ldir<0) name/**/2_all_bbox = name/**/2_all_bbox .and. name/**/2_bboxlo(4) /= 0 \
   &&     if (name/**/2_idir>0) name/**/2_all_bbox = name/**/2_all_bbox .and. name/**/2_bboxhi(1) /= 0 \
   &&     if (name/**/2_jdir>0) name/**/2_all_bbox = name/**/2_all_bbox .and. name/**/2_bboxhi(2) /= 0 \
   &&     if (name/**/2_kdir>0) name/**/2_all_bbox = name/**/2_all_bbox .and. name/**/2_bboxhi(3) /= 0 \
   &&     if (name/**/2_ldir>0) name/**/2_all_bbox = name/**/2_all_bbox .and. name/**/2_bboxhi(4) /= 0 \
   &&    if (name/**/2_all_bbox .and. name/**/2_any_bbox) then \
   &&       name/**/2_bmin(1) = name/**/2_blo(1)+1 \
   &&       name/**/2_bmin(2) = name/**/2_blo(2)+1 \
   &&       name/**/2_bmin(3) = name/**/2_blo(3)+1 \
   &&       name/**/2_bmin(4) = name/**/2_blo(4)+1 \
   &&       if (name/**/2_idir<0) name/**/2_bmin(1) = 1 \
   &&       if (name/**/2_jdir<0) name/**/2_bmin(2) = 1 \
   &&       if (name/**/2_kdir<0) name/**/2_bmin(3) = 1 \
   &&       if (name/**/2_ldir<0) name/**/2_bmin(4) = 1 \
   &&       if (name/**/2_idir>0) name/**/2_bmin(1) = cctk_lsh(1) - name/**/2_bhi(1) \
   &&       if (name/**/2_jdir>0) name/**/2_bmin(2) = cctk_lsh(2) - name/**/2_bhi(2) \
   &&       if (name/**/2_kdir>0) name/**/2_bmin(3) = cctk_lsh(3) - name/**/2_bhi(3) \
   &&       if (name/**/2_ldir>0) name/**/2_bmin(4) = cctk_lsh(4) - name/**/2_bhi(4) \
   &&       name/**/2_bmax(1) = cctk_lsh(1) - name/**/2_bhi(1) \
   &&       name/**/2_bmax(2) = cctk_lsh(2) - name/**/2_bhi(2) \
   &&       name/**/2_bmax(3) = cctk_lsh(3) - name/**/2_bhi(3) \
   &&       name/**/2_bmax(4) = cctk_lsh(4) - name/**/2_bhi(4) \
   &&       if (name/**/2_idir<0) name/**/2_bmax(1) = name/**/2_blo(1) \
   &&       if (name/**/2_jdir<0) name/**/2_bmax(2) = name/**/2_blo(2) \
   &&       if (name/**/2_kdir<0) name/**/2_bmax(3) = name/**/2_blo(3) \
   &&       if (name/**/2_ldir<0) name/**/2_bmax(4) = name/**/2_blo(4) \
   &&       if (name/**/2_idir>0) name/**/2_bmax(1) = cctk_lsh(1) \
   &&       if (name/**/2_jdir>0) name/**/2_bmax(2) = cctk_lsh(2) \
   &&       if (name/**/2_kdir>0) name/**/2_bmax(3) = cctk_lsh(3) \
   &&       if (name/**/2_ldir>0) name/**/2_bmax(4) = cctk_lsh(4) \
   &&       CCTK_LOOP4STR_NORMAL(name/**/_intboundaries, \
                                 i,j,k,l, \
                                 ni,nj,nk,nl, \
                                 name/**/2_idir,name/**/2_jdir,name/**/2_kdir,name/**/2_ldir, \
                                 name/**/2_bmin(1),name/**/2_bmin(2),name/**/2_bmin(3),name/**/2_bmin(4), \
                                 name/**/2_bmax(1),name/**/2_bmax(2),name/**/2_bmax(3),name/**/2_bmax(4), \
                                 cctk_ash(1), \
                                 cctk_ash(2), \
                                 cctk_ash(3), \
                                 cctk_ash(4), \
                                 imin,imax, name/**/2_istr) \

#define CCTK_ENDLOOP4STR_INTBOUNDARIES(name) \
            CCTK_ENDLOOP4STR_NORMAL(name/**/_intboundaries) \
   &&    end if /* bbox */ \
   && end do /* dir */ \
   && end do /* dir */ \
   && end do /* dir */ \
   && end do /* dir */ \



/* LOOP_ALL */

#define CCTK_LOOP4_ALL_DECLARE(name) \
   CCTK_LOOP4STR_ALL_DECLARE(name) \
   && integer :: name/**/3_dummy_imin, name/**/3_dummy_imax \

#define CCTK_LOOP4_ALL_OMP_PRIVATE(name) \
   CCTK_LOOP4STR_ALL_OMP_PRIVATE(name) \

#define CCTK_LOOP4_ALL(name, \
                       i,j,k,l) \
   CCTK_LOOP4STR_ALL(name, \
                     i,j,k,l, \
                     name/**/3_dummy_imin,name/**/3_dummy_imax, 1) \

#define CCTK_ENDLOOP4_ALL(name) \
   CCTK_ENDLOOP4STR_ALL(name) \

#define CCTK_LOOP4STR_ALL_DECLARE(name) \
   CCTK_LOOP4STR_DECLARE(name/**/_all) \

#define CCTK_LOOP4STR_ALL_OMP_PRIVATE(name) \
   CCTK_LOOP4STR_OMP_PRIVATE(name/**/_all) \

#define CCTK_LOOP4STR_ALL(name, \
                          i,j,k,l, \
                          imin,imax, cctki3_istr) \
   CCTK_LOOP4STR(name/**/_all, \
                 i,j,k,l, \
                 1,1,1,1, \
                 cctk_lsh(1),cctk_lsh(2),cctk_lsh(3),cctk_lsh(4), \
                 cctk_ash(1),cctk_ash(2),cctk_ash(3),cctk_ash(4), \
                 imin,imax, cctki3_istr) \

#define CCTK_ENDLOOP4STR_ALL(name) \
   CCTK_ENDLOOP4STR(name/**/_all) \



/* LOOP_INT */

#define CCTK_LOOP4_INT_DECLARE(name) \
   CCTK_LOOP4STR_INT_DECLARE(name) \
   && integer :: name/**/3_dummy_imin, name/**/3_dummy_imax \

#define CCTK_LOOP4_INT_OMP_PRIVATE(name) \
   CCTK_LOOP4STR_INT_OMP_PRIVATE(name) \

#define CCTK_LOOP4_INT(name, \
                        i,j,k,l) \
   CCTK_LOOP4STR_INT(name, \
                      i,j,k,l, \
                      name/**/3_dummy_imin,name/**/3_dummy_imax, 1) \

#define CCTK_ENDLOOP4_INT(name) \
   CCTK_ENDLOOP4STR_INT(name) \

#define CCTK_LOOP4STR_INT_DECLARE(name) \
   CCTK_LOOP4STR_INTERIOR_DECLARE(name/**/_int) \
   && CCTK_INT, parameter :: name/**/3_isize = 8 \
   && CCTK_INT :: name/**/3_bndsize    (8) \
   && CCTK_INT :: name/**/3_is_ghostbnd(8) \
   && CCTK_INT :: name/**/3_is_symbnd  (8) \
   && CCTK_INT :: name/**/3_is_physbnd (8) \
   && CCTK_INT :: name/**/3_ierr \

#define CCTK_LOOP4STR_INT_OMP_PRIVATE(name) \
   CCTK_LOOP4STR_INTERIOR_OMP_PRIVATE(name/**/_int) \

#define CCTK_LOOP4STR_INT(name, \
                          i,j,k,l, \
                          imin,imax, cctki3_istr) \
   && !$omp single \
   && name/**/3_ierr = GetBoundarySizesAndTypes \
         (cctkGH, name/**/3_isize, name/**/3_bndsize, name/**/3_is_ghostbnd, name/**/3_is_symbnd, name/**/3_is_physbnd) \
   && !$omp end single copyprivate(name/**/3_bndsize) \
   && CCTK_LOOP4STR_INTERIOR(name/**/_int, \
                             i,j,k,l, \
                             int(name/**/3_bndsize(1+1)),int(name/**/3_bndsize(3+1)),int(name/**/3_bndsize(5+1)),int(name/**/3_bndsize(7+1)), \
                             int(name/**/3_bndsize(2)),int(name/**/3_bndsize(4)),int(name/**/3_bndsize(6)),int(name/**/3_bndsize(8)), \
                             imin,imax, (cctki3_istr)) \

#define CCTK_ENDLOOP4STR_INT(name) \
      CCTK_ENDLOOP4STR_INTERIOR(name/**/_int) \



/* LOOP_BND */

#define CCTK_LOOP4_BND_DECLARE(name) \
   CCTK_LOOP4STR_BND_DECLARE(name) \
   && integer :: name/**/3_dummy_imin, name/**/3_dummy_imax \

#define CCTK_LOOP4_BND_OMP_PRIVATE(name) \
   CCTK_LOOP4STR_BND_OMP_PRIVATE(name) \

#define CCTK_LOOP4_BND(name, \
                       i,j,k,l, \
                       ni,nj,nk,nl) \
   CCTK_LOOP4STR_BND(name, \
                     i,j,k,l, \
                     ni,nj,nk,nl, \
                     name/**/3_dummy_imin,name/**/3_dummy_imax, 1) \

#define CCTK_ENDLOOP4_BND(name) \
   CCTK_ENDLOOP4STR_BND(name) \

#define CCTK_LOOP4STR_BND_DECLARE(name) \
   CCTK_LOOP4STR_BOUNDARIES_DECLARE(name/**/_bnd) \
   && CCTK_INT, parameter :: name/**/3_isize = 8 \
   && CCTK_INT :: name/**/3_bndsize    (8) \
   && CCTK_INT :: name/**/3_is_ghostbnd(8) \
   && CCTK_INT :: name/**/3_is_symbnd  (8) \
   && CCTK_INT :: name/**/3_is_physbnd (8) \
   && CCTK_INT :: name/**/3_ierr \

#define CCTK_LOOP4STR_BND_OMP_PRIVATE(name) \
   CCTK_LOOP4STR_BOUNDARIES_OMP_PRIVATE(name/**/_bnd) \

#define CCTK_LOOP4STR_BND(name, \
                          i,j,k,l, \
                          ni,nj,nk,nl, \
                          imin,imax, cctki3_istr) \
   && !$omp single \
   && name/**/3_ierr = GetBoundarySizesAndTypes \
         (cctkGH, name/**/3_isize, name/**/3_bndsize, name/**/3_is_ghostbnd, name/**/3_is_symbnd, name/**/3_is_physbnd) \
   && !$omp end single copyprivate(name/**/3_bndsize, name/**/3_is_physbnd) \
   && CCTK_LOOP4STR_BOUNDARIES(name/**/_bnd, \
                               i,j,k,l, \
                               ni,nj,nk,nl, \
                               int(name/**/3_bndsize(1))+1,int(name/**/3_bndsize(3))+1,int(name/**/3_bndsize(5))+1,int(name/**/3_bndsize(7))+1, \
                               int(name/**/3_bndsize(2)),int(name/**/3_bndsize(4)),int(name/**/3_bndsize(6)),int(name/**/3_bndsize(8)), \
                               int(name/**/3_is_physbnd(1)),int(name/**/3_is_physbnd(3)),int(name/**/3_is_physbnd(5)),int(name/**/3_is_physbnd(7)), \
                               int(name/**/3_is_physbnd(2)),int(name/**/3_is_physbnd(4)),int(name/**/3_is_physbnd(6)),int(name/**/3_is_physbnd(8)), \
                               imin,imax, (cctki3_istr)) \

#define CCTK_ENDLOOP4STR_BND(name) \
      CCTK_ENDLOOP4STR_BOUNDARIES(name/**/_bnd) \






/* LOOP_INTBND */

#define CCTK_LOOP4_INTBND_DECLARE(name) \
   CCTK_LOOP4STR_INTBND_DECLARE(name) \
   && integer :: name/**/3_dummy_imin, name/**/3_dummy_imax \

#define CCTK_LOOP4_INTBND_OMP_PRIVATE(name) \
   CCTK_LOOP4STR_INTBND_OMP_PRIVATE(name) \

#define CCTK_LOOP4_INTBND(name, \
                          i,j,k,l, \
                          ni,nj,nk,nl) \
   CCTK_LOOP4STR_INTBND(name, \
                        i,j,k,l, \
                        ni,nj,nk,nl, \
                        name/**/3_dummy_imin,name/**/3_dummy_imax, 1) \

#define CCTK_ENDLOOP4_INTBND(name) \
   CCTK_ENDLOOP4STR_INTBND(name) \

#define CCTK_LOOP4STR_INTBND_DECLARE(name) \
   CCTK_LOOP4STR_INTBOUNDARIES_DECLARE(name/**/_bnd) \
   && CCTK_INT, parameter :: name/**/3_isize = 8 \
   && CCTK_INT :: name/**/3_bndsize    (8) \
   && CCTK_INT :: name/**/3_is_ghostbnd(8) \
   && CCTK_INT :: name/**/3_is_symbnd  (8) \
   && CCTK_INT :: name/**/3_is_physbnd (8) \
   && CCTK_INT :: name/**/3_ierr \

#define CCTK_LOOP4STR_INTBND_OMP_PRIVATE(name) \
   CCTK_LOOP4STR_INTBOUNDARIES_OMP_PRIVATE(name/**/_bnd) \

#define CCTK_LOOP4STR_INTBND(name, \
                             i,j,k,l, \
                             ni,nj,nk,nl, \
                             imin,imax, cctki3_istr) \
   && !$omp single \
   && name/**/3_ierr = GetBoundarySizesAndTypes \
         (cctkGH, name/**/3_isize, name/**/3_bndsize, name/**/3_is_ghostbnd, name/**/3_is_symbnd, name/**/3_is_physbnd) \
   && !$omp end single copyprivate(name/**/3_bndsize, name/**/3_is_physbnd) \
   && CCTK_LOOP4STR_INTBOUNDARIES(name/**/_bnd, \
                                  i,j,k,l, \
                                  ni,nj,nk,nl, \
                                  int(name/**/3_bndsize(1+1)),int(name/**/3_bndsize(3+1)),int(name/**/3_bndsize(5+1)),int(name/**/3_bndsize(7+1)), \
                                  int(name/**/3_bndsize(2)),int(name/**/3_bndsize(4)),int(name/**/3_bndsize(6)),int(name/**/3_bndsize(8)), \
                                  int(name/**/3_is_physbnd(1)),int(name/**/3_is_physbnd(3)),int(name/**/3_is_physbnd(5)),int(name/**/3_is_physbnd(7)), \
                                  int(name/**/3_is_physbnd(2)),int(name/**/3_is_physbnd(4)),int(name/**/3_is_physbnd(6)),int(name/**/3_is_physbnd(8)), \
                                  imin,imax, (cctki3_istr)) \

#define CCTK_ENDLOOP4STR_INTBND(name) \
      CCTK_ENDLOOP4STR_INTBOUNDARIES(name/**/_bnd) \

#endif /* #ifdef FCODE */



#endif /* #ifndef _CCTK_LOOP_H_ */
