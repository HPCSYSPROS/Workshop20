#include <cctk.h>
#include <cctk_Parameters.h>

#include <loopcontrol.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>

#include "operator_prototypes_3d.hh"
#include "typeprops.hh"

using namespace std;

namespace CarpetLib {

#define SRCIND3(i, j, k)                                                       \
  index3(srcioff + (i), srcjoff + (j), srckoff + (k), srcipadext, srcjpadext,  \
         srckpadext, srciext, srcjext, srckext)
#define DSTIND3(i, j, k)                                                       \
  index3(dstioff + (i), dstjoff + (j), dstkoff + (k), dstipadext, dstjpadext,  \
         dstkpadext, dstiext, dstjext, dstkext)

template <typename T> inline T min3(T const &x, T const &y, T const &z) {
  return min(x, min(y, z));
}

template <typename T> inline T max3(T const &x, T const &y, T const &z) {
  return max(x, max(y, z));
}

template <typename T>
void interpolate_eno_3d_3tl(
    T const *restrict const src1, CCTK_REAL const t1,
    T const *restrict const src2, CCTK_REAL const t2,
    T const *restrict const src3, CCTK_REAL const t3,
    ivect3 const &restrict srcpadext, ivect3 const &restrict srcext,
    T *restrict const dst, CCTK_REAL const t, ivect3 const &restrict dstpadext,
    ivect3 const &restrict dstext, ibbox3 const &restrict srcbbox,
    ibbox3 const &restrict dstbbox, ibbox3 const &restrict,
    ibbox3 const &restrict regbbox, void *extraargs) {
  assert(not extraargs);

  typedef typename typeprops<T>::real RT;

  if (any(srcbbox.stride() != regbbox.stride() or
          dstbbox.stride() != regbbox.stride())) {
    CCTK_WARN(0, "Internal error: strides disagree");
  }

  if (any(srcbbox.stride() != dstbbox.stride())) {
    CCTK_WARN(0, "Internal error: strides disagree");
  }

  // This could be handled, but is likely to point to an error
  // elsewhere
  if (regbbox.empty()) {
    CCTK_WARN(0, "Internal error: region extent is empty");
  }

  if (not regbbox.is_contained_in(srcbbox) or
      not regbbox.is_contained_in(dstbbox)) {
    CCTK_WARN(0,
              "Internal error: region extent is not contained in array extent");
  }

  ivect3 const regext = regbbox.shape() / regbbox.stride();
  assert(all((regbbox.lower() - srcbbox.lower()) % srcbbox.stride() == 0));
  ivect3 const srcoff = (regbbox.lower() - srcbbox.lower()) / srcbbox.stride();
  assert(all((regbbox.lower() - dstbbox.lower()) % dstbbox.stride() == 0));
  ivect3 const dstoff = (regbbox.lower() - dstbbox.lower()) / dstbbox.stride();

  ptrdiff_t const srcipadext = srcpadext[0];
  ptrdiff_t const srcjpadext = srcpadext[1];
  ptrdiff_t const srckpadext = srcpadext[2];

  ptrdiff_t const dstipadext = dstpadext[0];
  ptrdiff_t const dstjpadext = dstpadext[1];
  ptrdiff_t const dstkpadext = dstpadext[2];

  ptrdiff_t const srciext = srcext[0];
  ptrdiff_t const srcjext = srcext[1];
  ptrdiff_t const srckext = srcext[2];

  ptrdiff_t const dstiext = dstext[0];
  ptrdiff_t const dstjext = dstext[1];
  ptrdiff_t const dstkext = dstext[2];

  ptrdiff_t const regiext = regext[0];
  ptrdiff_t const regjext = regext[1];
  ptrdiff_t const regkext = regext[2];

  ptrdiff_t const srcioff = srcoff[0];
  ptrdiff_t const srcjoff = srcoff[1];
  ptrdiff_t const srckoff = srcoff[2];

  ptrdiff_t const dstioff = dstoff[0];
  ptrdiff_t const dstjoff = dstoff[1];
  ptrdiff_t const dstkoff = dstoff[2];

  // Quadratic (second order) interpolation

  RT const tmin = min3(t1, t2, t3);
  RT const tmax = max3(t1, t2, t3);
  RT const eps = 1.0e-10 * (tmax - tmin);

  if (fabs(t1 - t2) < eps or fabs(t1 - t3) < eps or fabs(t2 - t3) < eps) {
    CCTK_WARN(0, "Internal error: arrays have same time");
  }
  if (t < min3(t1, t2, t3) - eps or t > max3(t1, t2, t3) + eps) {
    CCTK_WARN(0, "Internal error: extrapolation in time");
  }

  // Calculate stencil coefficients for 3-point and 2-point
  // interpolations
  RT const s1fac3 = (t - t2) * (t - t3) / ((t1 - t2) * (t1 - t3));
  RT const s2fac3 = (t - t1) * (t - t3) / ((t2 - t1) * (t2 - t3));
  RT const s3fac3 = (t - t1) * (t - t2) / ((t3 - t1) * (t3 - t2));

  RT const s1fac2_12 = (t - t2) / (t1 - t2);
  RT const s2fac2_12 = (t - t1) / (t2 - t1);

  RT const s2fac2_23 = (t - t3) / (t2 - t3);
  RT const s3fac2_23 = (t - t2) / (t3 - t2);

  // Choose which two time levels should be used for linear
  // interpolation
  bool const use_12 = t >= min(t1, t2) - eps and t <= max(t1, t2) + eps;
  bool const use_23 = t >= min(t2, t3) - eps and t <= max(t2, t3) + eps;
  assert(use_12 or use_23);
// TODO: Instead of use_12, calculate 3 coefficents that perform
// the desired 2-point interpolation, which would avoid the if
// statement in the loop, simplifying the code.

// Loop over region
#pragma omp parallel
  CCTK_LOOP3(interpolate_end_3d_3tl, i, j, k, 0, 0, 0, regiext, regjext,
             regkext, srcipadext, srcjpadext, srckpadext) {

    T const s1 = src1[SRCIND3(i, j, k)];
    T const s2 = src2[SRCIND3(i, j, k)];
    T const s3 = src3[SRCIND3(i, j, k)];

    // 3-point interpolation
    T d = s1fac3 * s1 + s2fac3 * s2 + s3fac3 * s3;

    // If the 3-point interpolation leads to a new extremum,
    // fall back to 2-point interpolation instead
    if (d > max3(s1, s2, s3) or d < min3(s1, s2, s3)) {
      if (use_12) {
        d = s1fac2_12 * s1 + s2fac2_12 * s2;
      } else {
        d = s2fac2_23 * s2 + s3fac2_23 * s3;
      }
    }

    dst[DSTIND3(i, j, k)] = d;
  }
  CCTK_ENDLOOP3(interpolate_end_3d_3tl);
}

#ifdef HAVE_CCTK_COMPLEX8
template <>
void interpolate_eno_3d_3tl(
    CCTK_COMPLEX8 const *restrict const src1, CCTK_REAL const t1,
    CCTK_COMPLEX8 const *restrict const src2, CCTK_REAL const t2,
    CCTK_COMPLEX8 const *restrict const src3, CCTK_REAL const t3,
    ivect3 const &restrict srcpadext, ivect3 const &restrict srcext,
    CCTK_COMPLEX8 *restrict const dst, CCTK_REAL const t,
    ivect3 const &restrict dstpadext, ivect3 const &restrict dstext,
    ibbox3 const &restrict srcbbox, ibbox3 const &restrict dstbbox,
    ibbox3 const &restrict, ibbox3 const &restrict regbbox, void *extraargs) {
  assert(not extraargs);

  CCTK_WARN(CCTK_WARN_ABORT, "ENO for complex numbers is not supported");
}
#endif

#ifdef HAVE_CCTK_COMPLEX16
template <>
void interpolate_eno_3d_3tl(
    CCTK_COMPLEX16 const *restrict const src1, CCTK_REAL const t1,
    CCTK_COMPLEX16 const *restrict const src2, CCTK_REAL const t2,
    CCTK_COMPLEX16 const *restrict const src3, CCTK_REAL const t3,
    ivect3 const &restrict srcpadext, ivect3 const &restrict srcext,
    CCTK_COMPLEX16 *restrict const dst, CCTK_REAL const t,
    ivect3 const &restrict dstpadext, ivect3 const &restrict dstext,
    ibbox3 const &restrict srcbbox, ibbox3 const &restrict dstbbox,
    ibbox3 const &restrict, ibbox3 const &restrict regbbox, void *extraargs) {
  CCTK_WARN(CCTK_WARN_ABORT, "ENO for complex numbers is not supported");
}
#endif

#ifdef HAVE_CCTK_COMPLEX32
template <>
void interpolate_eno_3d_3tl(
    CCTK_COMPLEX32 const *restrict const src1, CCTK_REAL const t1,
    CCTK_COMPLEX32 const *restrict const src2, CCTK_REAL const t2,
    CCTK_COMPLEX32 const *restrict const src3, CCTK_REAL const t3,
    ivect3 const &restrict srcpadext, ivect3 const &restrict srcext,
    CCTK_COMPLEX32 *restrict const dst, CCTK_REAL const t,
    ivect3 const &restrict dstpadext, ivect3 const &restrict dstext,
    ibbox3 const &restrict srcbbox, ibbox3 const &restrict dstbbox,
    ibbox3 const &restrict, ibbox3 const &restrict regbbox, void *extraargs) {
  CCTK_WARN(CCTK_WARN_ABORT, "ENO for complex numbers is not supported");
}
#endif

#define TYPECASE(N, T)                                                         \
  template void interpolate_eno_3d_3tl(                                        \
      T const *restrict const src1, CCTK_REAL const t1,                        \
      T const *restrict const src2, CCTK_REAL const t2,                        \
      T const *restrict const src3, CCTK_REAL const t3,                        \
      ivect3 const &restrict srcpadext, ivect3 const &restrict srcext,         \
      T *restrict const dst, CCTK_REAL const t,                                \
      ivect3 const &restrict dstpadext, ivect3 const &restrict dstext,         \
      ibbox3 const &restrict srcbbox, ibbox3 const &restrict dstbbox,          \
      ibbox3 const &restrict, ibbox3 const &restrict regbbox,                  \
      void *extraargs);
#define CARPET_NO_COMPLEX
#include "typecase.hh"
#undef TYPECASE

} // namespace CarpetLib
