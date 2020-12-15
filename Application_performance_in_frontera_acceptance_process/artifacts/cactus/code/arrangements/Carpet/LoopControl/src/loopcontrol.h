#ifndef LOOPCONTROL_H
#define LOOPCONTROL_H

/* This file uses the namespace LC_* for macros and lc_* for C
   identifiers. */

#define LC_DIM 3

#ifdef CCODE

#include <assert.h>
#include <stddef.h>
#include <stdlib.h>

#include <cctk.h>

/* #define lc_assert(x) ((void)0) */
#define lc_assert(x) assert(x)

#ifdef __cplusplus
extern "C" {
#endif

static inline ptrdiff_t lc_min(const ptrdiff_t i, const ptrdiff_t j) {
  return i < j ? i : j;
}

static inline ptrdiff_t lc_max(const ptrdiff_t i, const ptrdiff_t j) {
  return i > j ? i : j;
}

struct lc_thread_info_t;
struct lc_fine_thread_comm_t;

struct lc_descr_t;

typedef struct { ptrdiff_t v[LC_DIM]; } lc_vec_t;

typedef struct {
  /* Traverse pos from min (inclusive) to max (exclusive) with a
     stride of step. Equivalently, traverse idx from 0 (inclusive)
     to count (exclusive). */
  lc_vec_t min, max, step, pos;
  lc_vec_t count, idx;
} lc_space_t;

typedef struct {
  /* memory layout */
  lc_vec_t ash;

  /* overall loop bounds */
  lc_space_t overall;

  /* coarse threads */
  lc_space_t coarse_thread;
  /* shared between all threads */
  struct lc_thread_info_t *coarse_thread_info_ptr;
  int coarse_thread_done;

  /* coarse loop; count, idx, pos are undefined */
  lc_space_t coarse_loop;

  /* fine loop; count, idx, pos are undefined */
  lc_space_t fine_loop;

  /* fine threads; min, max, pos are undefined */
  lc_space_t fine_thread;

  /* selftest: shared between all threads */
  unsigned char *selftest_array;
} lc_control_t;

void lc_descr_init(struct lc_descr_t **descr, const char *name,
                   const char *file, int line);
void lc_control_init(lc_control_t *restrict control, struct lc_descr_t *descr,
                     ptrdiff_t imin, ptrdiff_t jmin, ptrdiff_t kmin,
                     ptrdiff_t imax, ptrdiff_t jmax, ptrdiff_t kmax,
                     ptrdiff_t iash, ptrdiff_t jash, ptrdiff_t kash,
                     ptrdiff_t istr);
void lc_control_finish(lc_control_t *restrict control,
                       struct lc_descr_t *descr);

void lc_thread_init(lc_control_t *restrict control);
int lc_thread_done(const lc_control_t *restrict control);
void lc_thread_step(lc_control_t *restrict control);

void lc_selftest_set(const lc_control_t *restrict control, ptrdiff_t imin,
                     ptrdiff_t imax, ptrdiff_t istr, ptrdiff_t i, ptrdiff_t j,
                     ptrdiff_t k);

#define LC_COARSE_SETUP(D)                                                     \
  lc_control.coarse_loop.min.v[D] = lc_control.coarse_thread.pos.v[D];         \
  lc_control.coarse_loop.max.v[D] = lc_min(                                    \
      lc_control.coarse_thread.max.v[D],                                       \
      lc_control.coarse_loop.min.v[D] + lc_control.coarse_thread.step.v[D]);   \
  const ptrdiff_t lc_cmin##D = lc_control.coarse_loop.min.v[D];                \
  const ptrdiff_t lc_cmax##D = lc_control.coarse_loop.max.v[D];                \
  const ptrdiff_t lc_cstep##D = lc_control.coarse_loop.step.v[D];
#define LC_COARSE_LOOP(D)                                                      \
  for (ptrdiff_t lc_cpos##D = lc_cmin##D; lc_cpos##D < lc_cmax##D;             \
       lc_cpos##D += lc_cstep##D)

#define LC_FINE_SETUP(D)                                                       \
  lc_control.fine_loop.min.v[D] = lc_cpos##D;                                  \
  lc_control.fine_loop.max.v[D] = lc_min(                                      \
      lc_control.coarse_loop.max.v[D],                                         \
      lc_control.fine_loop.min.v[D] + lc_control.coarse_loop.step.v[D]);       \
  /*const*/ ptrdiff_t lc_fmin##D = lc_control.fine_loop.min.v[D];              \
  /*const*/ ptrdiff_t lc_fmax##D = lc_control.fine_loop.max.v[D];              \
  const ptrdiff_t lc_fstep##D = lc_control.fine_loop.step.v[D];                \
  const ptrdiff_t lc_ftoff##D =                                                \
      lc_control.fine_thread.idx.v[D] * lc_control.fine_thread.step.v[D];
#define LC_FINE_LOOP(I, NI, D)                                                 \
  for (ptrdiff_t I = lc_fmin##D + lc_ftoff##D; I < lc_fmax##D;                 \
       I += lc_fstep##D) {                                                     \
    const ptrdiff_t NI CCTK_ATTRIBUTE_UNUSED =                                 \
        CCTK_BUILTIN_EXPECT(lc_dir##D == 0, 1)                                 \
            ? 0                                                                \
            : lc_dir##D < 0 ? I + 1 : lc_control.overall.max.v[D] - I;

#if VECTORISE && VECTORISE_ALIGNED_ARRAYS
/* Arrays are aligned: fmin0 is the aligned loop boundary; keep it,
   and set up imin to be the intended loop boundary */
#define LC_ALIGN(i, j, k, vec_imin, vec_imax)                                  \
  const ptrdiff_t vec_imin = lc_max(lc_control.overall.min.v[0], lc_fmin0);    \
  const ptrdiff_t vec_imax = lc_fmax0;                                         \
  lc_assert(lc_fmin0 >= 0);                                                    \
  lc_assert(lc_fmin0 < lc_fmax0);                                              \
  lc_assert(lc_fmax0 <= lc_control.overall.max.v[0]);                          \
  const ptrdiff_t lc_iminpos = lc_fmin0 + lc_ash0 * (j + lc_ash1 * k);         \
  const ptrdiff_t lc_iminoffset CCTK_ATTRIBUTE_UNUSED = lc_iminpos % lc_str0;  \
  const int lc_fmax0_is_outer = lc_fmax0 == lc_control.overall.max.v[0];       \
  const ptrdiff_t lc_imaxpos = lc_fmax0 + lc_ash0 * (j + lc_ash1 * k);         \
  const ptrdiff_t lc_imaxoffset CCTK_ATTRIBUTE_UNUSED = lc_imaxpos % lc_str0;  \
  lc_assert(lc_iminoffset == 0);                                               \
  if (!lc_fmax0_is_outer)                                                      \
    lc_assert(lc_imaxoffset == 0);                                             \
  lc_assert(vec_imin >= lc_control.overall.min.v[0]);                          \
  lc_assert(vec_imax <= lc_control.overall.max.v[0]);                          \
  lc_assert(vec_imin >= lc_fmin0);                                             \
  lc_assert(vec_imax <= lc_fmax0);                                             \
  lc_assert(vec_imin < vec_imax);
#else
/* Arrays are not aligned: fine.min[0] and fine.max[0] are the
   intended loop boundaries; override fmin0 and fmax0 to be aligned,
   and set imin and imax; this may move the fine loop boundaries
   except at outer boundaries to avoid partial stores */
#define LC_ALIGN(i, j, k, vec_imin, vec_imax)                                  \
  lc_fmin0 = lc_control.fine_loop.min.v[0];                                    \
  lc_fmax0 = lc_control.fine_loop.max.v[0];                                    \
  ptrdiff_t vec_imin = lc_fmin0;                                               \
  ptrdiff_t vec_imax = lc_fmax0;                                               \
  lc_assert(lc_fmin0 >= 0);                                                    \
  lc_assert(lc_fmin0 < lc_fmax0);                                              \
  lc_assert(lc_fmin0 >= lc_control.overall.min.v[0]);                          \
  lc_assert(lc_fmax0 <= lc_control.overall.max.v[0]);                          \
  const int lc_fmin0_is_outer = lc_fmin0 == lc_control.overall.min.v[0];       \
  const int lc_fmax0_is_outer = lc_fmax0 == lc_control.overall.max.v[0];       \
  const ptrdiff_t lc_iminpos = lc_fmin0 + lc_ash0 * (j + lc_ash1 * k);         \
  const ptrdiff_t lc_iminoffset = lc_iminpos % lc_str0;                        \
  const ptrdiff_t lc_imaxpos = lc_fmax0 + lc_ash0 * (j + lc_ash1 * k);         \
  const ptrdiff_t lc_imaxoffset = lc_imaxpos % lc_str0;                        \
  lc_fmin0 -= lc_iminoffset;                                                   \
  if (!lc_fmax0_is_outer)                                                      \
    lc_fmax0 -= lc_imaxoffset;                                                 \
  lc_assert(lc_fmin0 < lc_fmax0);                                              \
  if (!lc_fmin0_is_outer)                                                      \
    vec_imin = lc_fmin0;                                                       \
  if (!lc_fmax0_is_outer)                                                      \
    vec_imax = lc_fmax0;                                                       \
  lc_assert(vec_imin >= lc_control.overall.min.v[0]);                          \
  lc_assert(vec_imax <= lc_control.overall.max.v[0]);                          \
  lc_assert(vec_imin >= lc_fmin0);                                             \
  lc_assert(vec_imax <= lc_fmax0);                                             \
  lc_assert(vec_imin < vec_imax);
#endif

#define LC_SELFTEST(i, j, k, vec_imin, vec_imax)                               \
  if (CCTK_BUILTIN_EXPECT(lc_control.selftest_array != NULL, 0)) {             \
    lc_selftest_set(&lc_control, vec_imin, vec_imax, lc_str0, i, j, k);        \
  }

#define LC_LOOP3STR_NORMAL(name, i, j, k, ni, nj, nk, idir_, jdir_, kdir_,     \
                           imin_, jmin_, kmin_, imax_, jmax_, kmax_, iash_,    \
                           jash_, kash_, vec_imin, vec_imax, istr_)            \
  do {                                                                         \
    typedef int lc_loop3vec_##name;                                            \
                                                                               \
    const ptrdiff_t lc_dir0 CCTK_ATTRIBUTE_UNUSED = (idir_);                   \
    const ptrdiff_t lc_dir1 CCTK_ATTRIBUTE_UNUSED = (jdir_);                   \
    const ptrdiff_t lc_dir2 CCTK_ATTRIBUTE_UNUSED = (kdir_);                   \
                                                                               \
    const ptrdiff_t lc_ash0 CCTK_ATTRIBUTE_UNUSED = (iash_);                   \
    const ptrdiff_t lc_ash1 CCTK_ATTRIBUTE_UNUSED = (jash_);                   \
    const ptrdiff_t lc_ash2 CCTK_ATTRIBUTE_UNUSED = (kash_);                   \
                                                                               \
    const ptrdiff_t lc_str0 CCTK_ATTRIBUTE_UNUSED = (istr_);                   \
                                                                               \
    static struct lc_descr_t *lc_descr = NULL;                                 \
    lc_descr_init(&lc_descr, #name, __FILE__, __LINE__);                       \
                                                                               \
    lc_control_t lc_control;                                                   \
    lc_control_init(&lc_control, lc_descr, (imin_), (jmin_), (kmin_), (imax_), \
                    (jmax_), (kmax_), lc_ash0, lc_ash1, lc_ash2, lc_str0);     \
                                                                               \
    /* Multithreading */                                                       \
    for (lc_thread_init(&lc_control); !lc_thread_done(&lc_control);            \
         lc_thread_step(&lc_control)) {                                        \
                                                                               \
      /* Coarse loops */                                                       \
      LC_COARSE_SETUP(2)                                                       \
      LC_COARSE_SETUP(1)                                                       \
      LC_COARSE_SETUP(0)                                                       \
      LC_COARSE_LOOP(2) {                                                      \
        LC_COARSE_LOOP(1) {                                                    \
          LC_COARSE_LOOP(0) {                                                  \
                                                                               \
            /* Fine loops */                                                   \
            LC_FINE_SETUP(2)                                                   \
            LC_FINE_SETUP(1)                                                   \
            LC_FINE_SETUP(0)                                                   \
            LC_FINE_LOOP(k, nk, 2) {                                           \
              LC_FINE_LOOP(j, nj, 1) {                                         \
                LC_ALIGN(i, j, k, vec_imin, vec_imax)                          \
                LC_FINE_LOOP(i, ni, 0) {                                       \
                  LC_SELFTEST(i, j, k, vec_imin, vec_imax) {

#define LC_ENDLOOP3STR_NORMAL(name)                                            \
  } /* body */                                                                 \
  }                                                                            \
  }                                                                            \
  }                                                                            \
  }                                                                            \
  }                                                                            \
  } /* fine */                                                                 \
  }                                                                            \
  }                                                                            \
  } /* coarse */                                                               \
  } /* multithreading */                                                       \
  lc_control_finish(&lc_control, lc_descr);                                    \
  typedef lc_loop3vec_##name lc_ensure_proper_nesting CCTK_ATTRIBUTE_UNUSED;   \
  }                                                                            \
  while (0)

/* Definitions to ensure compatibility with earlier versions of
   LoopControl */
#define LC_LOOP3VEC(name, i, j, k, imin, jmin, kmin, imax, jmax, kmax, iash,   \
                    jash, kash, vec_imin, vec_imax, istr)                      \
  LC_LOOP3STR_NORMAL(name, i, j, k, lc_ni, lc_nj, lc_nk, 0, 0, 0, imin, jmin,  \
                     kmin, imax, jmax, kmax, iash, jash, kash, vec_imin,       \
                     vec_imax, istr)
#define LC_ENDLOOP3VEC(name) LC_ENDLOOP3STR_NORMAL(name)

#define LC_LOOP3(name, i, j, k, imin, jmin, kmin, imax, jmax, kmax, iash,      \
                 jash, kash)                                                   \
  LC_LOOP3VEC(name, i, j, k, imin, jmin, kmin, imax, jmax, kmax, iash, jash,   \
              kash, lc_vec_imin, lc_vec_imax, 1)
#define LC_ENDLOOP3(name) LC_ENDLOOP3VEC(name)

/* Replace CCTK_LOOP macros */
#if !defined CCTK_LOOP3STR_NORMAL || !defined CCTK_ENDLOOP3STR_NORMAL
#error "internal error"
#endif
#undef CCTK_LOOP3STR_NORMAL
#undef CCTK_ENDLOOP3STR_NORMAL
#define CCTK_LOOP3STR_NORMAL(name, i, j, k, ni, nj, nk, idir, jdir, kdir,      \
                             imin, jmin, kmin, imax, jmax, kmax, iash, jash,   \
                             kash, vec_imin, vec_imax, istr)                   \
  LC_LOOP3STR_NORMAL(name, i, j, k, ni, nj, nk, idir, jdir, kdir, imin, jmin,  \
                     kmin, imax, jmax, kmax, iash, jash, kash, vec_imin,       \
                     vec_imax, istr)
#define CCTK_ENDLOOP3STR_NORMAL(name) LC_ENDLOOP3STR_NORMAL(name)

#ifdef __cplusplus
}
#endif

#endif /* #ifdef CCODE */

#ifdef FCODE
#include "loopcontrol_fortran.h"
#endif

#endif /* #ifndef LOOPCONTROL_H */
