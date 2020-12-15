#ifndef OPERATOR_PROTOTYPES_3D
#define OPERATOR_PROTOTYPES_3D

#include <cctk.h>

#include <cstdlib>

#include "defs.hh"
#include "bbox.hh"
#include "vect.hh"

#include "operator_prototypes.hh"

namespace CarpetLib {

using namespace std;

static inline size_t index3(size_t const i, size_t const j, size_t const k,
                            size_t const padexti, size_t const padextj,
                            size_t const padextk, size_t const exti,
                            size_t const extj, size_t const extk) {
#ifdef CARPET_DEBUG
  assert(static_cast<ptrdiff_t>(i) >= 0 and i < exti);
  assert(static_cast<ptrdiff_t>(j) >= 0 and j < extj);
  assert(static_cast<ptrdiff_t>(k) >= 0 and k < extk);
#endif

  return i + padexti * (j + padextj * k);
}

static inline size_t offset3(size_t const i, size_t const j, size_t const k,
                             size_t const exti, size_t const extj,
                             size_t const extk) {
  return i + exti * (j + extj * k);
}

static int const dim3 = 3;

typedef vect<bool, dim3> bvect3;
typedef vect<int, dim3> ivect3;
typedef bbox<int, dim3> ibbox3;

template <typename T>
void copy_3d(T const *restrict const src, ivect3 const &restrict srcpadext,
             ivect3 const &restrict srcext, T *restrict const dst,
             ivect3 const &restrict dstpadext, ivect3 const &restrict dstext,
             ibbox3 const &restrict srcbbox, ibbox3 const &restrict dstbbox,
             ibbox3 const &restrict srcregbbox,
             ibbox3 const &restrict dstregbbox, void *extraargs);

template <typename T, int ORDER>
void prolongate_3d_rf2(T const *restrict const src,
                       ivect3 const &restrict srcpadext,
                       ivect3 const &restrict srcext, T *restrict const dst,
                       ivect3 const &restrict dstpadext,
                       ivect3 const &restrict dstext,
                       ibbox3 const &restrict srcbbox,
                       ibbox3 const &restrict dstbbox,
                       ibbox3 const &restrict srcregbbox,
                       ibbox3 const &restrict dstregbbox, void *extraargs);

template <typename T, int ORDER>
void prolongate_3d_stagger011(
    T const *restrict const src, ivect3 const &restrict srcpadext,
    ivect3 const &restrict srcext, T *restrict const dst,
    ivect3 const &restrict dstpadext, ivect3 const &restrict dstext,
    ibbox3 const &restrict srcbbox, ibbox3 const &restrict dstbbox,
    ibbox3 const &restrict srcregbbox, ibbox3 const &restrict dstregbbox,
    void *extraargs);

template <typename T, int ORDER>
void prolongate_3d_stagger101(
    T const *restrict const src, ivect3 const &restrict srcpadext,
    ivect3 const &restrict srcext, T *restrict const dst,
    ivect3 const &restrict dstpadext, ivect3 const &restrict dstext,
    ibbox3 const &restrict srcbbox, ibbox3 const &restrict dstbbox,
    ibbox3 const &restrict srcregbbox, ibbox3 const &restrict dstregbbox,
    void *extraargs);

template <typename T, int ORDER>
void prolongate_3d_stagger110(
    T const *restrict const src, ivect3 const &restrict srcpadext,
    ivect3 const &restrict srcext, T *restrict const dst,
    ivect3 const &restrict dstpadext, ivect3 const &restrict dstext,
    ibbox3 const &restrict srcbbox, ibbox3 const &restrict dstbbox,
    ibbox3 const &restrict srcregbbox, ibbox3 const &restrict dstregbbox,
    void *extraargs);

template <typename T, int ORDER>
void prolongate_3d_stagger111(
    T const *restrict const src, ivect3 const &restrict srcpadext,
    ivect3 const &restrict srcext, T *restrict const dst,
    ivect3 const &restrict dstpadext, ivect3 const &restrict dstext,
    ibbox3 const &restrict srcbbox, ibbox3 const &restrict dstbbox,
    ibbox3 const &restrict srcregbbox, ibbox3 const &restrict dstregbbox,
    void *extraargs);

template <typename T>
void prolongate_3d_o5_monotone_rf2(
    T const *restrict const src, ivect3 const &restrict srcpadext,
    ivect3 const &restrict srcext, T *restrict const dst,
    ivect3 const &restrict dstpadext, ivect3 const &restrict dstext,
    ibbox3 const &restrict srcbbox, ibbox3 const &restrict dstbbox,
    ibbox3 const &restrict srcregbbox, ibbox3 const &restrict dstregbbox,
    void *extraargs);

template <typename T, int ORDER>
void prolongate_3d_cc_rf2(T const *restrict const src,
                          ivect3 const &restrict srcpadext,
                          ivect3 const &restrict srcext, T *restrict const dst,
                          ivect3 const &restrict dstpadext,
                          ivect3 const &restrict dstext,
                          ibbox3 const &restrict srcbbox,
                          ibbox3 const &restrict dstbbox,
                          ibbox3 const &restrict srcregbbox,
                          ibbox3 const &restrict dstregbbox, void *extraargs);

template <typename T, int ORDER>
void prolongate_3d_cc_eno_rf2(
    T const *restrict const src, ivect3 const &restrict srcpadext,
    ivect3 const &restrict srcext, T *restrict const dst,
    ivect3 const &restrict dstpadext, ivect3 const &restrict dstext,
    ibbox3 const &restrict srcbbox, ibbox3 const &restrict dstbbox,
    ibbox3 const &restrict srcregbbox, ibbox3 const &restrict dstregbbox,
    void *extraargs);

template <typename T, int ORDER>
void prolongate_3d_cc_eno_rf2(
    T const *restrict const src, ivect3 const &restrict srcpadext,
    ivect3 const &restrict srcext, T *restrict const dst,
    ivect3 const &restrict dstpadext, ivect3 const &restrict dstext,
    ibbox3 const &restrict srcbbox, ibbox3 const &restrict dstbbox,
    ibbox3 const &restrict regbbox);

template <typename T, int ORDER>
void prolongate_3d_dgfe_rf2(
    T const *restrict const src, ivect3 const &restrict srcpadext,
    ivect3 const &restrict srcext, T *restrict const dst,
    ivect3 const &restrict dstpadext, ivect3 const &restrict dstext,
    ibbox3 const &restrict srcbbox, ibbox3 const &restrict dstbbox,
    ibbox3 const &restrict srcregbbox, ibbox3 const &restrict dstregbbox,
    void *extraargs);

template <typename T>
void restrict_3d_rf2(T const *restrict const src,
                     ivect3 const &restrict srcpadext,
                     ivect3 const &restrict srcext, T *restrict const dst,
                     ivect3 const &restrict dstpadext,
                     ivect3 const &restrict dstext,
                     ibbox3 const &restrict srcbbox,
                     ibbox3 const &restrict dstbbox,
                     ibbox3 const &restrict srcregbbox,
                     ibbox3 const &restrict dstregbbox, void *extraargs);

template <typename T>
void interpolate_3d_2tl(T const *restrict const src1, CCTK_REAL const t1,
                        T const *restrict const src2, CCTK_REAL const t2,
                        ivect3 const &restrict srcpadext,
                        ivect3 const &restrict srcext, T *restrict const dst,
                        CCTK_REAL const t, ivect3 const &restrict dstpadext,
                        ivect3 const &restrict dstext,
                        ibbox3 const &restrict srcbbox,
                        ibbox3 const &restrict dstbbox,
                        ibbox3 const &restrict srcregbbox,
                        ibbox3 const &restrict dstregbbox, void *extraargs);

template <typename T>
void interpolate_3d_3tl(T const *restrict const src1, CCTK_REAL const t1,
                        T const *restrict const src2, CCTK_REAL const t2,
                        T const *restrict const src3, CCTK_REAL const t3,
                        ivect3 const &restrict srcpadext,
                        ivect3 const &restrict srcext, T *restrict const dst,
                        CCTK_REAL const t, ivect3 const &restrict dstpadext,
                        ivect3 const &restrict dstext,
                        ibbox3 const &restrict srcbbox,
                        ibbox3 const &restrict dstbbox,
                        ibbox3 const &restrict srcregbbox,
                        ibbox3 const &restrict dstregbbox, void *extraargs);

template <typename T>
void interpolate_3d_4tl(T const *restrict const src1, CCTK_REAL const t1,
                        T const *restrict const src2, CCTK_REAL const t2,
                        T const *restrict const src3, CCTK_REAL const t3,
                        T const *restrict const src4, CCTK_REAL const t4,
                        ivect3 const &restrict srcpadext,
                        ivect3 const &restrict srcext, T *restrict const dst,
                        CCTK_REAL const t, ivect3 const &restrict dstpadext,
                        ivect3 const &restrict dstext,
                        ibbox3 const &restrict srcbbox,
                        ibbox3 const &restrict dstbbox,
                        ibbox3 const &restrict srcregbbox,
                        ibbox3 const &restrict dstregbbox, void *extraargs);

template <typename T>
void interpolate_3d_5tl(T const *restrict const src1, CCTK_REAL const t1,
                        T const *restrict const src2, CCTK_REAL const t2,
                        T const *restrict const src3, CCTK_REAL const t3,
                        T const *restrict const src4, CCTK_REAL const t4,
                        T const *restrict const src5, CCTK_REAL const t5,
                        ivect3 const &restrict srcpadext,
                        ivect3 const &restrict srcext, T *restrict const dst,
                        CCTK_REAL const t, ivect3 const &restrict dstpadext,
                        ivect3 const &restrict dstext,
                        ibbox3 const &restrict srcbbox,
                        ibbox3 const &restrict dstbbox,
                        ibbox3 const &restrict srcregbbox,
                        ibbox3 const &restrict dstregbbox, void *extraargs);

template <typename T>
void interpolate_eno_3d_3tl(
    T const *restrict const src1, CCTK_REAL const t1,
    T const *restrict const src2, CCTK_REAL const t2,
    T const *restrict const src3, CCTK_REAL const t3,
    ivect3 const &restrict srcpadext, ivect3 const &restrict srcext,
    T *restrict const dst, CCTK_REAL const t, ivect3 const &restrict dstpadext,
    ivect3 const &restrict dstext, ibbox3 const &restrict srcbbox,
    ibbox3 const &restrict dstbbox, ibbox3 const &restrict srcregbbox,
    ibbox3 const &restrict dstregbbox, void *extraargs);

template <typename T>
void restrict_3d_cc_rf2(T const *restrict const src,
                        ivect3 const &restrict srcpadext,
                        ivect3 const &restrict srcext, T *restrict const dst,
                        ivect3 const &restrict dstpadext,
                        ivect3 const &restrict dstext,
                        ibbox3 const &restrict srcbbox,
                        ibbox3 const &restrict dstbbox,
                        ibbox3 const &restrict srcregbbox,
                        ibbox3 const &restrict dstregbbox, void *extraargs);

template <typename T>
void restrict_3d_cc_o3_rf2(T const *restrict const src,
                           ivect3 const &restrict srcpadext,
                           ivect3 const &restrict srcext, T *restrict const dst,
                           ivect3 const &restrict dstpadext,
                           ivect3 const &restrict dstext,
                           ibbox3 const &restrict srcbbox,
                           ibbox3 const &restrict dstbbox,
                           ibbox3 const &restrict srcregbbox,
                           ibbox3 const &restrict dstregbbox, void *extraargs);

template <typename T>
void restrict_3d_cc_o5_rf2(T const *restrict const src,
                           ivect3 const &restrict srcpadext,
                           ivect3 const &restrict srcext, T *restrict const dst,
                           ivect3 const &restrict dstpadext,
                           ivect3 const &restrict dstext,
                           ibbox3 const &restrict srcbbox,
                           ibbox3 const &restrict dstbbox,
                           ibbox3 const &restrict srcregbbox,
                           ibbox3 const &restrict dstregbbox, void *extraargs);

template <typename T, int centi, int centj, int centk>
void restrict_3d_vc_rf2(T const *restrict const src,
                        ivect3 const &restrict srcpadext,
                        ivect3 const &restrict srcext, T *restrict const dst,
                        ivect3 const &restrict dstpadext,
                        ivect3 const &restrict dstext,
                        ibbox3 const &restrict srcbbox,
                        ibbox3 const &restrict dstbbox,
                        ibbox3 const &restrict srcregbbox,
                        ibbox3 const &restrict dstregbbox, void *extraargs);

template <typename T, int ORDER>
void restrict_3d_dgfe_rf2(T const *restrict const src,
                          ivect3 const &restrict srcpadext,
                          ivect3 const &restrict srcext, T *restrict const dst,
                          ivect3 const &restrict dstpadext,
                          ivect3 const &restrict dstext,
                          ibbox3 const &restrict srcbbox,
                          ibbox3 const &restrict dstbbox,
                          ibbox3 const &restrict srcregbbox,
                          ibbox3 const &restrict dstregbbox, void *extraargs);

template <typename t>
void restrict_3d_stagger011(
    t const *restrict const src, ivect3 const &restrict srcpadext,
    ivect3 const &restrict srcext, t *restrict const dst,
    ivect3 const &restrict dstpadext, ivect3 const &restrict dstext,
    ibbox3 const &restrict srcbbox, ibbox3 const &restrict dstbbox,
    ibbox3 const &restrict srcregbbox, ibbox3 const &restrict dstregbbox,
    void *extraargs);

template <typename t>
void restrict_3d_stagger101(
    t const *restrict const src, ivect3 const &restrict srcpadext,
    ivect3 const &restrict srcext, t *restrict const dst,
    ivect3 const &restrict dstpadext, ivect3 const &restrict dstext,
    ibbox3 const &restrict srcbbox, ibbox3 const &restrict dstbbox,
    ibbox3 const &restrict srcregbbox, ibbox3 const &restrict dstregbbox,
    void *extraargs);

template <typename t>
void restrict_3d_stagger110(
    t const *restrict const src, ivect3 const &restrict srcpadext,
    ivect3 const &restrict srcext, t *restrict const dst,
    ivect3 const &restrict dstpadext, ivect3 const &restrict dstext,
    ibbox3 const &restrict srcbbox, ibbox3 const &restrict dstbbox,
    ibbox3 const &restrict srcregbbox, ibbox3 const &restrict dstregbbox,
    void *extraargs);

template <typename t>
void restrict_3d_stagger111(
    t const *restrict const src, ivect3 const &restrict srcpadext,
    ivect3 const &restrict srcext, t *restrict const dst,
    ivect3 const &restrict dstpadext, ivect3 const &restrict dstext,
    ibbox3 const &restrict srcbbox, ibbox3 const &restrict dstbbox,
    ibbox3 const &restrict srcregbbox, ibbox3 const &restrict dstregbbox,
    void *extraargs);

} // namespace CarpetLib

#endif // #ifndef OPERATOR_PROTOTYPES_3D
