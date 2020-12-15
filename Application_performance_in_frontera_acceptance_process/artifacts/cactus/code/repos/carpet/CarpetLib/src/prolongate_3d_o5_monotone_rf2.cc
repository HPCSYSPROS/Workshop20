// This is meant to reproduce the prolongation algorithm used in the
// SACRA code (based on IH's interpretation of their papers and
// comments in talks, so it might be an idea for someone to talk to
// them! Of course, given that this is "general purpose" and SACRA is
// very specific in the variables converted, it probably won't be
// possible to get a perfect reproduction).
//
// The idea is to use fifth order Lagrange interpolation based on the
// nearest 6 points (in any one dimension). However, we must also
// ensure monotonicity. To do this we check that the result of the
// fifth order result (which is just copied from prolongate_3d_o5_rf2)
// is monotonic with respect to the relevant neighbours), and if not
// we impose linear interpolation instead (from prolongate_3d_o1_rf2).
//
// Note that this code does not work for complex GFs (due to the use
// of the max and min intrinsics).

#include <cctk.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>

#include "operator_prototypes_3d.hh"
#include "typeprops.hh"

using namespace std;

namespace CarpetLib {

#define SRCIND3(i, j, k)                                                       \
  index3(i, j, k, srcipadext, srcjpadext, srckpadext, srciext, srcjext, srckext)
#define DSTIND3(i, j, k)                                                       \
  index3(i, j, k, dstipadext, dstjpadext, dstkpadext, dstiext, dstjext, dstkext)

template <typename T>
inline T min4(T const &x1, T const &x2, T const &x3, T const &x4) {
  return min(min(x1, x2), min(x3, x4));
}

template <typename T>
inline T max4(T const &x1, T const &x2, T const &x3, T const &x4) {
  return max(max(x1, x2), max(x3, x4));
}

template <typename T>
inline T min8(T const &x1, T const &x2, T const &x3, T const &x4, T const &x5,
              T const &x6, T const &x7, T const &x8) {
  return min(min(min(x1, x2), min(x3, x4)), min(min(x5, x6), min(x7, x8)));
}

template <typename T>
inline T max8(T const &x1, T const &x2, T const &x3, T const &x4, T const &x5,
              T const &x6, T const &x7, T const &x8) {
  return max(max(max(x1, x2), max(x3, x4)), max(max(x5, x6), max(x7, x8)));
}

template <typename T>
void prolongate_3d_o5_monotone_rf2(
    T const *restrict const src, ivect3 const &restrict srcpadext,
    ivect3 const &restrict srcext, T *restrict const dst,
    ivect3 const &restrict dstpadext, ivect3 const &restrict dstext,
    ibbox3 const &restrict srcbbox, ibbox3 const &restrict dstbbox,
    ibbox3 const &restrict, ibbox3 const &restrict regbbox, void *extraargs) {
  assert(not extraargs);

  typedef typename typeprops<T>::real RT;

  if (any(srcbbox.stride() <= regbbox.stride() or
          dstbbox.stride() != regbbox.stride())) {
    CCTK_WARN(0, "Internal error: strides disagree");
  }

  if (any(srcbbox.stride() != reffact2 * dstbbox.stride())) {
    CCTK_WARN(
        0,
        "Internal error: source strides are not twice the destination strides");
  }

  // This could be handled, but is likely to point to an error
  // elsewhere
  if (regbbox.empty()) {
    CCTK_WARN(0, "Internal error: region extent is empty");
  }

  ivect3 const regext = regbbox.shape() / regbbox.stride();
  assert(all((regbbox.lower() - srcbbox.lower()) % regbbox.stride() == 0));
  ivect3 const srcoff = (regbbox.lower() - srcbbox.lower()) / regbbox.stride();
  assert(all((regbbox.lower() - dstbbox.lower()) % regbbox.stride() == 0));
  ivect3 const dstoff = (regbbox.lower() - dstbbox.lower()) / regbbox.stride();

  bvect3 const needoffsetlo = srcoff % reffact2 != 0 or regext > 1;
  bvect3 const needoffsethi =
      (srcoff + regext - 1) % reffact2 != 0 or regext > 1;
  ivect3 const offsetlo = either(needoffsetlo, 3, 0);
  ivect3 const offsethi = either(needoffsethi, 3, 0);

  if (not regbbox.expand(offsetlo, offsethi).is_contained_in(srcbbox) or
      not regbbox.is_contained_in(dstbbox)) {
    CCTK_WARN(0,
              "Internal error: region extent is not contained in array extent");
  }

  size_t const srcipadext = srcpadext[0];
  size_t const srcjpadext = srcpadext[1];
  size_t const srckpadext = srcpadext[2];

  size_t const dstipadext = dstpadext[0];
  size_t const dstjpadext = dstpadext[1];
  size_t const dstkpadext = dstpadext[2];

  size_t const srciext = srcext[0];
  size_t const srcjext = srcext[1];
  size_t const srckext = srcext[2];

  size_t const dstiext = dstext[0];
  size_t const dstjext = dstext[1];
  size_t const dstkext = dstext[2];

  size_t const regiext = regext[0];
  size_t const regjext = regext[1];
  size_t const regkext = regext[2];

  size_t const srcioff = srcoff[0];
  size_t const srcjoff = srcoff[1];
  size_t const srckoff = srcoff[2];

  size_t const dstioff = dstoff[0];
  size_t const dstjoff = dstoff[1];
  size_t const dstkoff = dstoff[2];

  size_t const fi = srcioff % 2;
  size_t const fj = srcjoff % 2;
  size_t const fk = srckoff % 2;

  size_t const i0 = srcioff / 2;
  size_t const j0 = srcjoff / 2;
  size_t const k0 = srckoff / 2;

  RT const one = 1;

  RT const f1 = 3 * one / 256;
  RT const f2 = -25 * one / 256;
  RT const f3 = 150 * one / 256;
  RT const f4 = 150 * one / 256;
  RT const f5 = -25 * one / 256;
  RT const f6 = 3 * one / 256;

  RT const o1_f1 = one / 2;
  RT const o1_f2 = one / 2;

  // Loop over fine region
  // Label scheme: l 8 fk fj fi

  size_t is, js, ks;
  size_t id, jd, kd;
  size_t i, j, k;

  // begin k loop
  k = 0;
  ks = k0;
  kd = dstkoff;
  if (fk == 0)
    goto l80;
  goto l81;

// begin j loop
l80:
  j = 0;
  js = j0;
  jd = dstjoff;
  if (fj == 0)
    goto l800;
  goto l801;

// begin i loop
l800:
  i = 0;
  is = i0;
  id = dstioff;
  if (fi == 0)
    goto l8000;
  goto l8001;

// kernel
l8000:
  dst[DSTIND3(id, jd, kd)] = src[SRCIND3(is, js, ks)];
  i = i + 1;
  id = id + 1;
  if (i < regiext)
    goto l8001;
  goto l900;

// kernel
l8001:
  dst[DSTIND3(id, jd, kd)] =
      +f1 * src[SRCIND3(is - 2, js, ks)] + f2 * src[SRCIND3(is - 1, js, ks)] +
      f3 * src[SRCIND3(is, js, ks)] + f4 * src[SRCIND3(is + 1, js, ks)] +
      f5 * src[SRCIND3(is + 2, js, ks)] + f6 * src[SRCIND3(is + 3, js, ks)];
  // Monotonicity enforcement
  if ((dst[DSTIND3(id, jd, kd)] >
       max(src[SRCIND3(is, js, ks)], src[SRCIND3(is + 1, js, ks)])) ||
      (dst[DSTIND3(id, jd, kd)] <
       min(src[SRCIND3(is, js, ks)], src[SRCIND3(is + 1, js, ks)]))) {
    dst[DSTIND3(id, jd, kd)] = +o1_f1 * src[SRCIND3(is, js, ks)] +
                               o1_f2 * src[SRCIND3(is + 1, js, ks)];
  }

  i = i + 1;
  id = id + 1;
  is = is + 1;
  if (i < regiext)
    goto l8000;
  goto l900;

// end i loop
l900:
  j = j + 1;
  jd = jd + 1;
  if (j < regjext)
    goto l801;
  goto l90;

// begin i loop
l801:
  i = 0;
  is = i0;
  id = dstioff;
  if (fi == 0)
    goto l8010;
  goto l8011;

// kernel
l8010:
  dst[DSTIND3(id, jd, kd)] =
      +f1 * src[SRCIND3(is, js - 2, ks)] + f2 * src[SRCIND3(is, js - 1, ks)] +
      f3 * src[SRCIND3(is, js, ks)] + f4 * src[SRCIND3(is, js + 1, ks)] +
      f5 * src[SRCIND3(is, js + 2, ks)] + f6 * src[SRCIND3(is, js + 3, ks)];
  // Monotonicity enforcement
  if ((dst[DSTIND3(id, jd, kd)] >
       max(src[SRCIND3(is, js, ks)], src[SRCIND3(is, js + 1, ks)])) ||
      (dst[DSTIND3(id, jd, kd)] <
       min(src[SRCIND3(is, js, ks)], src[SRCIND3(is, js + 1, ks)]))) {
    dst[DSTIND3(id, jd, kd)] = +o1_f1 * src[SRCIND3(is, js, ks)] +
                               o1_f2 * src[SRCIND3(is, js + 1, ks)];
  }
  i = i + 1;
  id = id + 1;
  if (i < regiext)
    goto l8011;
  goto l901;

// kernel
l8011:
  dst[DSTIND3(id, jd, kd)] = +f1 * f1 * src[SRCIND3(is - 2, js - 2, ks)] +
                             f2 * f1 * src[SRCIND3(is - 1, js - 2, ks)] +
                             f3 * f1 * src[SRCIND3(is, js - 2, ks)] +
                             f4 * f1 * src[SRCIND3(is + 1, js - 2, ks)] +
                             f5 * f1 * src[SRCIND3(is + 2, js - 2, ks)] +
                             f6 * f1 * src[SRCIND3(is + 3, js - 2, ks)] +
                             f1 * f2 * src[SRCIND3(is - 2, js - 1, ks)] +
                             f2 * f2 * src[SRCIND3(is - 1, js - 1, ks)] +
                             f3 * f2 * src[SRCIND3(is, js - 1, ks)] +
                             f4 * f2 * src[SRCIND3(is + 1, js - 1, ks)] +
                             f5 * f2 * src[SRCIND3(is + 2, js - 1, ks)] +
                             f6 * f2 * src[SRCIND3(is + 3, js - 1, ks)] +
                             f1 * f3 * src[SRCIND3(is - 2, js, ks)] +
                             f2 * f3 * src[SRCIND3(is - 1, js, ks)] +
                             f3 * f3 * src[SRCIND3(is, js, ks)] +
                             f4 * f3 * src[SRCIND3(is + 1, js, ks)] +
                             f5 * f3 * src[SRCIND3(is + 2, js, ks)] +
                             f6 * f3 * src[SRCIND3(is + 3, js, ks)] +
                             f1 * f4 * src[SRCIND3(is - 2, js + 1, ks)] +
                             f2 * f4 * src[SRCIND3(is - 1, js + 1, ks)] +
                             f3 * f4 * src[SRCIND3(is, js + 1, ks)] +
                             f4 * f4 * src[SRCIND3(is + 1, js + 1, ks)] +
                             f5 * f4 * src[SRCIND3(is + 2, js + 1, ks)] +
                             f6 * f4 * src[SRCIND3(is + 3, js + 1, ks)] +
                             f1 * f5 * src[SRCIND3(is - 2, js + 2, ks)] +
                             f2 * f5 * src[SRCIND3(is - 1, js + 2, ks)] +
                             f3 * f5 * src[SRCIND3(is, js + 2, ks)] +
                             f4 * f5 * src[SRCIND3(is + 1, js + 2, ks)] +
                             f5 * f5 * src[SRCIND3(is + 2, js + 2, ks)] +
                             f6 * f5 * src[SRCIND3(is + 3, js + 2, ks)] +
                             f1 * f6 * src[SRCIND3(is - 2, js + 3, ks)] +
                             f2 * f6 * src[SRCIND3(is - 1, js + 3, ks)] +
                             f3 * f6 * src[SRCIND3(is, js + 3, ks)] +
                             f4 * f6 * src[SRCIND3(is + 1, js + 3, ks)] +
                             f5 * f6 * src[SRCIND3(is + 2, js + 3, ks)] +
                             f6 * f6 * src[SRCIND3(is + 3, js + 3, ks)];
  // Monotonicity enforcement
  if ((dst[DSTIND3(id, jd, kd)] >
       max4(src[SRCIND3(is, js, ks)], src[SRCIND3(is + 1, js, ks)],
            src[SRCIND3(is, js + 1, ks)], src[SRCIND3(is + 1, js + 1, ks)])) ||
      (dst[DSTIND3(id, jd, kd)] <
       min4(src[SRCIND3(is, js, ks)], src[SRCIND3(is + 1, js, ks)],
            src[SRCIND3(is, js + 1, ks)], src[SRCIND3(is + 1, js + 1, ks)]))) {
    dst[DSTIND3(id, jd, kd)] = +o1_f1 * o1_f1 * src[SRCIND3(is, js, ks)] +
                               o1_f2 * o1_f1 * src[SRCIND3(is + 1, js, ks)] +
                               o1_f1 * o1_f2 * src[SRCIND3(is, js + 1, ks)] +
                               o1_f2 * o1_f2 * src[SRCIND3(is + 1, js + 1, ks)];
  }
  i = i + 1;
  id = id + 1;
  is = is + 1;
  if (i < regiext)
    goto l8010;
  goto l901;

// end i loop
l901:
  j = j + 1;
  jd = jd + 1;
  js = js + 1;
  if (j < regjext)
    goto l800;
  goto l90;

// end j loop
l90:
  k = k + 1;
  kd = kd + 1;
  if (k < regkext)
    goto l81;
  goto l9;

// begin j loop
l81:
  j = 0;
  js = j0;
  jd = dstjoff;
  if (fj == 0)
    goto l810;
  goto l811;

// begin i loop
l810:
  i = 0;
  is = i0;
  id = dstioff;
  if (fi == 0)
    goto l8100;
  goto l8101;

// kernel
l8100:
  dst[DSTIND3(id, jd, kd)] =
      +f1 * src[SRCIND3(is, js, ks - 2)] + f2 * src[SRCIND3(is, js, ks - 1)] +
      f3 * src[SRCIND3(is, js, ks)] + f4 * src[SRCIND3(is, js, ks + 1)] +
      f5 * src[SRCIND3(is, js, ks + 2)] + f6 * src[SRCIND3(is, js, ks + 3)];
  // Monotonicity enforcement
  if ((dst[DSTIND3(id, jd, kd)] >
       max(src[SRCIND3(is, js, ks)], src[SRCIND3(is, js, ks + 1)])) ||
      (dst[DSTIND3(id, jd, kd)] <
       min(src[SRCIND3(is, js, ks)], src[SRCIND3(is, js, ks + 1)]))) {
    dst[DSTIND3(id, jd, kd)] = +o1_f1 * src[SRCIND3(is, js, ks)] +
                               o1_f2 * src[SRCIND3(is, js, ks + 1)];
  }
  i = i + 1;
  id = id + 1;
  if (i < regiext)
    goto l8101;
  goto l910;

// kernel
l8101:
  dst[DSTIND3(id, jd, kd)] = +f1 * f1 * src[SRCIND3(is - 2, js, ks - 2)] +
                             f2 * f1 * src[SRCIND3(is - 1, js, ks - 2)] +
                             f3 * f1 * src[SRCIND3(is, js, ks - 2)] +
                             f4 * f1 * src[SRCIND3(is + 1, js, ks - 2)] +
                             f5 * f1 * src[SRCIND3(is + 2, js, ks - 2)] +
                             f6 * f1 * src[SRCIND3(is + 3, js, ks - 2)] +
                             f1 * f2 * src[SRCIND3(is - 2, js, ks - 1)] +
                             f2 * f2 * src[SRCIND3(is - 1, js, ks - 1)] +
                             f3 * f2 * src[SRCIND3(is, js, ks - 1)] +
                             f4 * f2 * src[SRCIND3(is + 1, js, ks - 1)] +
                             f5 * f2 * src[SRCIND3(is + 2, js, ks - 1)] +
                             f6 * f2 * src[SRCIND3(is + 3, js, ks - 1)] +
                             f1 * f3 * src[SRCIND3(is - 2, js, ks)] +
                             f2 * f3 * src[SRCIND3(is - 1, js, ks)] +
                             f3 * f3 * src[SRCIND3(is, js, ks)] +
                             f4 * f3 * src[SRCIND3(is + 1, js, ks)] +
                             f5 * f3 * src[SRCIND3(is + 2, js, ks)] +
                             f6 * f3 * src[SRCIND3(is + 3, js, ks)] +
                             f1 * f4 * src[SRCIND3(is - 2, js, ks + 1)] +
                             f2 * f4 * src[SRCIND3(is - 1, js, ks + 1)] +
                             f3 * f4 * src[SRCIND3(is, js, ks + 1)] +
                             f4 * f4 * src[SRCIND3(is + 1, js, ks + 1)] +
                             f5 * f4 * src[SRCIND3(is + 2, js, ks + 1)] +
                             f6 * f4 * src[SRCIND3(is + 3, js, ks + 1)] +
                             f1 * f5 * src[SRCIND3(is - 2, js, ks + 2)] +
                             f2 * f5 * src[SRCIND3(is - 1, js, ks + 2)] +
                             f3 * f5 * src[SRCIND3(is, js, ks + 2)] +
                             f4 * f5 * src[SRCIND3(is + 1, js, ks + 2)] +
                             f5 * f5 * src[SRCIND3(is + 2, js, ks + 2)] +
                             f6 * f5 * src[SRCIND3(is + 3, js, ks + 2)] +
                             f1 * f6 * src[SRCIND3(is - 2, js, ks + 3)] +
                             f2 * f6 * src[SRCIND3(is - 1, js, ks + 3)] +
                             f3 * f6 * src[SRCIND3(is, js, ks + 3)] +
                             f4 * f6 * src[SRCIND3(is + 1, js, ks + 3)] +
                             f5 * f6 * src[SRCIND3(is + 2, js, ks + 3)] +
                             f6 * f6 * src[SRCIND3(is + 3, js, ks + 3)];
  // Monotonicity enforcement
  if ((dst[DSTIND3(id, jd, kd)] >
       max4(src[SRCIND3(is, js, ks)], src[SRCIND3(is + 1, js, ks)],
            src[SRCIND3(is, js, ks + 1)], src[SRCIND3(is + 1, js, ks + 1)])) ||
      (dst[DSTIND3(id, jd, kd)] <
       min4(src[SRCIND3(is, js, ks)], src[SRCIND3(is + 1, js, ks)],
            src[SRCIND3(is, js, ks + 1)], src[SRCIND3(is + 1, js, ks + 1)]))) {
    dst[DSTIND3(id, jd, kd)] = +o1_f1 * o1_f1 * src[SRCIND3(is, js, ks)] +
                               o1_f2 * o1_f1 * src[SRCIND3(is + 1, js, ks)] +
                               o1_f1 * o1_f2 * src[SRCIND3(is, js, ks + 1)] +
                               o1_f2 * o1_f2 * src[SRCIND3(is + 1, js, ks + 1)];
  }
  i = i + 1;
  id = id + 1;
  is = is + 1;
  if (i < regiext)
    goto l8100;
  goto l910;

// end i loop
l910:
  j = j + 1;
  jd = jd + 1;
  if (j < regjext)
    goto l811;
  goto l91;

// begin i loop
l811:
  i = 0;
  is = i0;
  id = dstioff;
  if (fi == 0)
    goto l8110;
  goto l8111;

// kernel
l8110:
  dst[DSTIND3(id, jd, kd)] = +f1 * f1 * src[SRCIND3(is, js - 2, ks - 2)] +
                             f2 * f1 * src[SRCIND3(is, js - 1, ks - 2)] +
                             f3 * f1 * src[SRCIND3(is, js, ks - 2)] +
                             f4 * f1 * src[SRCIND3(is, js + 1, ks - 2)] +
                             f5 * f1 * src[SRCIND3(is, js + 2, ks - 2)] +
                             f6 * f1 * src[SRCIND3(is, js + 3, ks - 2)] +
                             f1 * f2 * src[SRCIND3(is, js - 2, ks - 1)] +
                             f2 * f2 * src[SRCIND3(is, js - 1, ks - 1)] +
                             f3 * f2 * src[SRCIND3(is, js, ks - 1)] +
                             f4 * f2 * src[SRCIND3(is, js + 1, ks - 1)] +
                             f5 * f2 * src[SRCIND3(is, js + 2, ks - 1)] +
                             f6 * f2 * src[SRCIND3(is, js + 3, ks - 1)] +
                             f1 * f3 * src[SRCIND3(is, js - 2, ks)] +
                             f2 * f3 * src[SRCIND3(is, js - 1, ks)] +
                             f3 * f3 * src[SRCIND3(is, js, ks)] +
                             f4 * f3 * src[SRCIND3(is, js + 1, ks)] +
                             f5 * f3 * src[SRCIND3(is, js + 2, ks)] +
                             f6 * f3 * src[SRCIND3(is, js + 3, ks)] +
                             f1 * f4 * src[SRCIND3(is, js - 2, ks + 1)] +
                             f2 * f4 * src[SRCIND3(is, js - 1, ks + 1)] +
                             f3 * f4 * src[SRCIND3(is, js, ks + 1)] +
                             f4 * f4 * src[SRCIND3(is, js + 1, ks + 1)] +
                             f5 * f4 * src[SRCIND3(is, js + 2, ks + 1)] +
                             f6 * f4 * src[SRCIND3(is, js + 3, ks + 1)] +
                             f1 * f5 * src[SRCIND3(is, js - 2, ks + 2)] +
                             f2 * f5 * src[SRCIND3(is, js - 1, ks + 2)] +
                             f3 * f5 * src[SRCIND3(is, js, ks + 2)] +
                             f4 * f5 * src[SRCIND3(is, js + 1, ks + 2)] +
                             f5 * f5 * src[SRCIND3(is, js + 2, ks + 2)] +
                             f6 * f5 * src[SRCIND3(is, js + 3, ks + 2)] +
                             f1 * f6 * src[SRCIND3(is, js - 2, ks + 3)] +
                             f2 * f6 * src[SRCIND3(is, js - 1, ks + 3)] +
                             f3 * f6 * src[SRCIND3(is, js, ks + 3)] +
                             f4 * f6 * src[SRCIND3(is, js + 1, ks + 3)] +
                             f5 * f6 * src[SRCIND3(is, js + 2, ks + 3)] +
                             f6 * f6 * src[SRCIND3(is, js + 3, ks + 3)];
  // Monotonicity enforcement
  if ((dst[DSTIND3(id, jd, kd)] >
       max4(src[SRCIND3(is, js, ks)], src[SRCIND3(is, js + 1, ks)],
            src[SRCIND3(is, js, ks + 1)], src[SRCIND3(is, js + 1, ks + 1)])) ||
      (dst[DSTIND3(id, jd, kd)] <
       min4(src[SRCIND3(is, js, ks)], src[SRCIND3(is, js + 1, ks)],
            src[SRCIND3(is, js, ks + 1)], src[SRCIND3(is, js + 1, ks + 1)]))) {
    dst[DSTIND3(id, jd, kd)] = +o1_f1 * o1_f1 * src[SRCIND3(is, js, ks)] +
                               o1_f2 * o1_f1 * src[SRCIND3(is, js + 1, ks)] +
                               o1_f1 * o1_f2 * src[SRCIND3(is, js, ks + 1)] +
                               o1_f2 * o1_f2 * src[SRCIND3(is, js + 1, ks + 1)];
  }
  i = i + 1;
  id = id + 1;
  if (i < regiext)
    goto l8111;
  goto l911;

// kernel
l8111 : {
  T const res1 = +f1 * f1 * f1 * src[SRCIND3(is - 2, js - 2, ks - 2)] +
                 f2 * f1 * f1 * src[SRCIND3(is - 1, js - 2, ks - 2)] +
                 f3 * f1 * f1 * src[SRCIND3(is, js - 2, ks - 2)] +
                 f4 * f1 * f1 * src[SRCIND3(is + 1, js - 2, ks - 2)] +
                 f5 * f1 * f1 * src[SRCIND3(is + 2, js - 2, ks - 2)] +
                 f6 * f1 * f1 * src[SRCIND3(is + 3, js - 2, ks - 2)] +
                 f1 * f2 * f1 * src[SRCIND3(is - 2, js - 1, ks - 2)] +
                 f2 * f2 * f1 * src[SRCIND3(is - 1, js - 1, ks - 2)] +
                 f3 * f2 * f1 * src[SRCIND3(is, js - 1, ks - 2)] +
                 f4 * f2 * f1 * src[SRCIND3(is + 1, js - 1, ks - 2)] +
                 f5 * f2 * f1 * src[SRCIND3(is + 2, js - 1, ks - 2)] +
                 f6 * f2 * f1 * src[SRCIND3(is + 3, js - 1, ks - 2)] +
                 f1 * f3 * f1 * src[SRCIND3(is - 2, js, ks - 2)] +
                 f2 * f3 * f1 * src[SRCIND3(is - 1, js, ks - 2)] +
                 f3 * f3 * f1 * src[SRCIND3(is, js, ks - 2)] +
                 f4 * f3 * f1 * src[SRCIND3(is + 1, js, ks - 2)] +
                 f5 * f3 * f1 * src[SRCIND3(is + 2, js, ks - 2)] +
                 f6 * f3 * f1 * src[SRCIND3(is + 3, js, ks - 2)] +
                 f1 * f4 * f1 * src[SRCIND3(is - 2, js + 1, ks - 2)] +
                 f2 * f4 * f1 * src[SRCIND3(is - 1, js + 1, ks - 2)] +
                 f3 * f4 * f1 * src[SRCIND3(is, js + 1, ks - 2)] +
                 f4 * f4 * f1 * src[SRCIND3(is + 1, js + 1, ks - 2)] +
                 f5 * f4 * f1 * src[SRCIND3(is + 2, js + 1, ks - 2)] +
                 f6 * f4 * f1 * src[SRCIND3(is + 3, js + 1, ks - 2)] +
                 f1 * f5 * f1 * src[SRCIND3(is - 2, js + 2, ks - 2)] +
                 f2 * f5 * f1 * src[SRCIND3(is - 1, js + 2, ks - 2)] +
                 f3 * f5 * f1 * src[SRCIND3(is, js + 2, ks - 2)] +
                 f4 * f5 * f1 * src[SRCIND3(is + 1, js + 2, ks - 2)] +
                 f5 * f5 * f1 * src[SRCIND3(is + 2, js + 2, ks - 2)] +
                 f6 * f5 * f1 * src[SRCIND3(is + 3, js + 2, ks - 2)] +
                 f1 * f6 * f1 * src[SRCIND3(is - 2, js + 3, ks - 2)] +
                 f2 * f6 * f1 * src[SRCIND3(is - 1, js + 3, ks - 2)] +
                 f3 * f6 * f1 * src[SRCIND3(is, js + 3, ks - 2)] +
                 f4 * f6 * f1 * src[SRCIND3(is + 1, js + 3, ks - 2)] +
                 f5 * f6 * f1 * src[SRCIND3(is + 2, js + 3, ks - 2)] +
                 f6 * f6 * f1 * src[SRCIND3(is + 3, js + 3, ks - 2)];
  T const res2 = +f1 * f1 * f2 * src[SRCIND3(is - 2, js - 2, ks - 1)] +
                 f2 * f1 * f2 * src[SRCIND3(is - 1, js - 2, ks - 1)] +
                 f3 * f1 * f2 * src[SRCIND3(is, js - 2, ks - 1)] +
                 f4 * f1 * f2 * src[SRCIND3(is + 1, js - 2, ks - 1)] +
                 f5 * f1 * f2 * src[SRCIND3(is + 2, js - 2, ks - 1)] +
                 f6 * f1 * f2 * src[SRCIND3(is + 3, js - 2, ks - 1)] +
                 f1 * f2 * f2 * src[SRCIND3(is - 2, js - 1, ks - 1)] +
                 f2 * f2 * f2 * src[SRCIND3(is - 1, js - 1, ks - 1)] +
                 f3 * f2 * f2 * src[SRCIND3(is, js - 1, ks - 1)] +
                 f4 * f2 * f2 * src[SRCIND3(is + 1, js - 1, ks - 1)] +
                 f5 * f2 * f2 * src[SRCIND3(is + 2, js - 1, ks - 1)] +
                 f6 * f2 * f2 * src[SRCIND3(is + 3, js - 1, ks - 1)] +
                 f1 * f3 * f2 * src[SRCIND3(is - 2, js, ks - 1)] +
                 f2 * f3 * f2 * src[SRCIND3(is - 1, js, ks - 1)] +
                 f3 * f3 * f2 * src[SRCIND3(is, js, ks - 1)] +
                 f4 * f3 * f2 * src[SRCIND3(is + 1, js, ks - 1)] +
                 f5 * f3 * f2 * src[SRCIND3(is + 2, js, ks - 1)] +
                 f6 * f3 * f2 * src[SRCIND3(is + 3, js, ks - 1)] +
                 f1 * f4 * f2 * src[SRCIND3(is - 2, js + 1, ks - 1)] +
                 f2 * f4 * f2 * src[SRCIND3(is - 1, js + 1, ks - 1)] +
                 f3 * f4 * f2 * src[SRCIND3(is, js + 1, ks - 1)] +
                 f4 * f4 * f2 * src[SRCIND3(is + 1, js + 1, ks - 1)] +
                 f5 * f4 * f2 * src[SRCIND3(is + 2, js + 1, ks - 1)] +
                 f6 * f4 * f2 * src[SRCIND3(is + 3, js + 1, ks - 1)] +
                 f1 * f5 * f2 * src[SRCIND3(is - 2, js + 2, ks - 1)] +
                 f2 * f5 * f2 * src[SRCIND3(is - 1, js + 2, ks - 1)] +
                 f3 * f5 * f2 * src[SRCIND3(is, js + 2, ks - 1)] +
                 f4 * f5 * f2 * src[SRCIND3(is + 1, js + 2, ks - 1)] +
                 f5 * f5 * f2 * src[SRCIND3(is + 2, js + 2, ks - 1)] +
                 f6 * f5 * f2 * src[SRCIND3(is + 3, js + 2, ks - 1)] +
                 f1 * f6 * f2 * src[SRCIND3(is - 2, js + 3, ks - 1)] +
                 f2 * f6 * f2 * src[SRCIND3(is - 1, js + 3, ks - 1)] +
                 f3 * f6 * f2 * src[SRCIND3(is, js + 3, ks - 1)] +
                 f4 * f6 * f2 * src[SRCIND3(is + 1, js + 3, ks - 1)] +
                 f5 * f6 * f2 * src[SRCIND3(is + 2, js + 3, ks - 1)] +
                 f6 * f6 * f2 * src[SRCIND3(is + 3, js + 3, ks - 1)];
  T const res3 = +f1 * f1 * f3 * src[SRCIND3(is - 2, js - 2, ks)] +
                 f2 * f1 * f3 * src[SRCIND3(is - 1, js - 2, ks)] +
                 f3 * f1 * f3 * src[SRCIND3(is, js - 2, ks)] +
                 f4 * f1 * f3 * src[SRCIND3(is + 1, js - 2, ks)] +
                 f5 * f1 * f3 * src[SRCIND3(is + 2, js - 2, ks)] +
                 f6 * f1 * f3 * src[SRCIND3(is + 3, js - 2, ks)] +
                 f1 * f2 * f3 * src[SRCIND3(is - 2, js - 1, ks)] +
                 f2 * f2 * f3 * src[SRCIND3(is - 1, js - 1, ks)] +
                 f3 * f2 * f3 * src[SRCIND3(is, js - 1, ks)] +
                 f4 * f2 * f3 * src[SRCIND3(is + 1, js - 1, ks)] +
                 f5 * f2 * f3 * src[SRCIND3(is + 2, js - 1, ks)] +
                 f6 * f2 * f3 * src[SRCIND3(is + 3, js - 1, ks)] +
                 f1 * f3 * f3 * src[SRCIND3(is - 2, js, ks)] +
                 f2 * f3 * f3 * src[SRCIND3(is - 1, js, ks)] +
                 f3 * f3 * f3 * src[SRCIND3(is, js, ks)] +
                 f4 * f3 * f3 * src[SRCIND3(is + 1, js, ks)] +
                 f5 * f3 * f3 * src[SRCIND3(is + 2, js, ks)] +
                 f6 * f3 * f3 * src[SRCIND3(is + 3, js, ks)] +
                 f1 * f4 * f3 * src[SRCIND3(is - 2, js + 1, ks)] +
                 f2 * f4 * f3 * src[SRCIND3(is - 1, js + 1, ks)] +
                 f3 * f4 * f3 * src[SRCIND3(is, js + 1, ks)] +
                 f4 * f4 * f3 * src[SRCIND3(is + 1, js + 1, ks)] +
                 f5 * f4 * f3 * src[SRCIND3(is + 2, js + 1, ks)] +
                 f6 * f4 * f3 * src[SRCIND3(is + 3, js + 1, ks)] +
                 f1 * f5 * f3 * src[SRCIND3(is - 2, js + 2, ks)] +
                 f2 * f5 * f3 * src[SRCIND3(is - 1, js + 2, ks)] +
                 f3 * f5 * f3 * src[SRCIND3(is, js + 2, ks)] +
                 f4 * f5 * f3 * src[SRCIND3(is + 1, js + 2, ks)] +
                 f5 * f5 * f3 * src[SRCIND3(is + 2, js + 2, ks)] +
                 f6 * f5 * f3 * src[SRCIND3(is + 3, js + 2, ks)] +
                 f1 * f6 * f3 * src[SRCIND3(is - 2, js + 3, ks)] +
                 f2 * f6 * f3 * src[SRCIND3(is - 1, js + 3, ks)] +
                 f3 * f6 * f3 * src[SRCIND3(is, js + 3, ks)] +
                 f4 * f6 * f3 * src[SRCIND3(is + 1, js + 3, ks)] +
                 f5 * f6 * f3 * src[SRCIND3(is + 2, js + 3, ks)] +
                 f6 * f6 * f3 * src[SRCIND3(is + 3, js + 3, ks)];
  T const res4 = +f1 * f1 * f4 * src[SRCIND3(is - 2, js - 2, ks + 1)] +
                 f2 * f1 * f4 * src[SRCIND3(is - 1, js - 2, ks + 1)] +
                 f3 * f1 * f4 * src[SRCIND3(is, js - 2, ks + 1)] +
                 f4 * f1 * f4 * src[SRCIND3(is + 1, js - 2, ks + 1)] +
                 f5 * f1 * f4 * src[SRCIND3(is + 2, js - 2, ks + 1)] +
                 f6 * f1 * f4 * src[SRCIND3(is + 3, js - 2, ks + 1)] +
                 f1 * f2 * f4 * src[SRCIND3(is - 2, js - 1, ks + 1)] +
                 f2 * f2 * f4 * src[SRCIND3(is - 1, js - 1, ks + 1)] +
                 f3 * f2 * f4 * src[SRCIND3(is, js - 1, ks + 1)] +
                 f4 * f2 * f4 * src[SRCIND3(is + 1, js - 1, ks + 1)] +
                 f5 * f2 * f4 * src[SRCIND3(is + 2, js - 1, ks + 1)] +
                 f6 * f2 * f4 * src[SRCIND3(is + 3, js - 1, ks + 1)] +
                 f1 * f3 * f4 * src[SRCIND3(is - 2, js, ks + 1)] +
                 f2 * f3 * f4 * src[SRCIND3(is - 1, js, ks + 1)] +
                 f3 * f3 * f4 * src[SRCIND3(is, js, ks + 1)] +
                 f4 * f3 * f4 * src[SRCIND3(is + 1, js, ks + 1)] +
                 f5 * f3 * f4 * src[SRCIND3(is + 2, js, ks + 1)] +
                 f6 * f3 * f4 * src[SRCIND3(is + 3, js, ks + 1)] +
                 f1 * f4 * f4 * src[SRCIND3(is - 2, js + 1, ks + 1)] +
                 f2 * f4 * f4 * src[SRCIND3(is - 1, js + 1, ks + 1)] +
                 f3 * f4 * f4 * src[SRCIND3(is, js + 1, ks + 1)] +
                 f4 * f4 * f4 * src[SRCIND3(is + 1, js + 1, ks + 1)] +
                 f5 * f4 * f4 * src[SRCIND3(is + 2, js + 1, ks + 1)] +
                 f6 * f4 * f4 * src[SRCIND3(is + 3, js + 1, ks + 1)] +
                 f1 * f5 * f4 * src[SRCIND3(is - 2, js + 2, ks + 1)] +
                 f2 * f5 * f4 * src[SRCIND3(is - 1, js + 2, ks + 1)] +
                 f3 * f5 * f4 * src[SRCIND3(is, js + 2, ks + 1)] +
                 f4 * f5 * f4 * src[SRCIND3(is + 1, js + 2, ks + 1)] +
                 f5 * f5 * f4 * src[SRCIND3(is + 2, js + 2, ks + 1)] +
                 f6 * f5 * f4 * src[SRCIND3(is + 3, js + 2, ks + 1)] +
                 f1 * f6 * f4 * src[SRCIND3(is - 2, js + 3, ks + 1)] +
                 f2 * f6 * f4 * src[SRCIND3(is - 1, js + 3, ks + 1)] +
                 f3 * f6 * f4 * src[SRCIND3(is, js + 3, ks + 1)] +
                 f4 * f6 * f4 * src[SRCIND3(is + 1, js + 3, ks + 1)] +
                 f5 * f6 * f4 * src[SRCIND3(is + 2, js + 3, ks + 1)] +
                 f6 * f6 * f4 * src[SRCIND3(is + 3, js + 3, ks + 1)];
  T const res5 = +f1 * f1 * f5 * src[SRCIND3(is - 2, js - 2, ks + 2)] +
                 f2 * f1 * f5 * src[SRCIND3(is - 1, js - 2, ks + 2)] +
                 f3 * f1 * f5 * src[SRCIND3(is, js - 2, ks + 2)] +
                 f4 * f1 * f5 * src[SRCIND3(is + 1, js - 2, ks + 2)] +
                 f5 * f1 * f5 * src[SRCIND3(is + 2, js - 2, ks + 2)] +
                 f6 * f1 * f5 * src[SRCIND3(is + 3, js - 2, ks + 2)] +
                 f1 * f2 * f5 * src[SRCIND3(is - 2, js - 1, ks + 2)] +
                 f2 * f2 * f5 * src[SRCIND3(is - 1, js - 1, ks + 2)] +
                 f3 * f2 * f5 * src[SRCIND3(is, js - 1, ks + 2)] +
                 f4 * f2 * f5 * src[SRCIND3(is + 1, js - 1, ks + 2)] +
                 f5 * f2 * f5 * src[SRCIND3(is + 2, js - 1, ks + 2)] +
                 f6 * f2 * f5 * src[SRCIND3(is + 3, js - 1, ks + 2)] +
                 f1 * f3 * f5 * src[SRCIND3(is - 2, js, ks + 2)] +
                 f2 * f3 * f5 * src[SRCIND3(is - 1, js, ks + 2)] +
                 f3 * f3 * f5 * src[SRCIND3(is, js, ks + 2)] +
                 f4 * f3 * f5 * src[SRCIND3(is + 1, js, ks + 2)] +
                 f5 * f3 * f5 * src[SRCIND3(is + 2, js, ks + 2)] +
                 f6 * f3 * f5 * src[SRCIND3(is + 3, js, ks + 2)] +
                 f1 * f4 * f5 * src[SRCIND3(is - 2, js + 1, ks + 2)] +
                 f2 * f4 * f5 * src[SRCIND3(is - 1, js + 1, ks + 2)] +
                 f3 * f4 * f5 * src[SRCIND3(is, js + 1, ks + 2)] +
                 f4 * f4 * f5 * src[SRCIND3(is + 1, js + 1, ks + 2)] +
                 f5 * f4 * f5 * src[SRCIND3(is + 2, js + 1, ks + 2)] +
                 f6 * f4 * f5 * src[SRCIND3(is + 3, js + 1, ks + 2)] +
                 f1 * f5 * f5 * src[SRCIND3(is - 2, js + 2, ks + 2)] +
                 f2 * f5 * f5 * src[SRCIND3(is - 1, js + 2, ks + 2)] +
                 f3 * f5 * f5 * src[SRCIND3(is, js + 2, ks + 2)] +
                 f4 * f5 * f5 * src[SRCIND3(is + 1, js + 2, ks + 2)] +
                 f5 * f5 * f5 * src[SRCIND3(is + 2, js + 2, ks + 2)] +
                 f6 * f5 * f5 * src[SRCIND3(is + 3, js + 2, ks + 2)] +
                 f1 * f6 * f5 * src[SRCIND3(is - 2, js + 3, ks + 2)] +
                 f2 * f6 * f5 * src[SRCIND3(is - 1, js + 3, ks + 2)] +
                 f3 * f6 * f5 * src[SRCIND3(is, js + 3, ks + 2)] +
                 f4 * f6 * f5 * src[SRCIND3(is + 1, js + 3, ks + 2)] +
                 f5 * f6 * f5 * src[SRCIND3(is + 2, js + 3, ks + 2)] +
                 f6 * f6 * f5 * src[SRCIND3(is + 3, js + 3, ks + 2)];
  T const res6 = +f1 * f1 * f6 * src[SRCIND3(is - 2, js - 2, ks + 3)] +
                 f2 * f1 * f6 * src[SRCIND3(is - 1, js - 2, ks + 3)] +
                 f3 * f1 * f6 * src[SRCIND3(is, js - 2, ks + 3)] +
                 f4 * f1 * f6 * src[SRCIND3(is + 1, js - 2, ks + 3)] +
                 f5 * f1 * f6 * src[SRCIND3(is + 2, js - 2, ks + 3)] +
                 f6 * f1 * f6 * src[SRCIND3(is + 3, js - 2, ks + 3)] +
                 f1 * f2 * f6 * src[SRCIND3(is - 2, js - 1, ks + 3)] +
                 f2 * f2 * f6 * src[SRCIND3(is - 1, js - 1, ks + 3)] +
                 f3 * f2 * f6 * src[SRCIND3(is, js - 1, ks + 3)] +
                 f4 * f2 * f6 * src[SRCIND3(is + 1, js - 1, ks + 3)] +
                 f5 * f2 * f6 * src[SRCIND3(is + 2, js - 1, ks + 3)] +
                 f6 * f2 * f6 * src[SRCIND3(is + 3, js - 1, ks + 3)] +
                 f1 * f3 * f6 * src[SRCIND3(is - 2, js, ks + 3)] +
                 f2 * f3 * f6 * src[SRCIND3(is - 1, js, ks + 3)] +
                 f3 * f3 * f6 * src[SRCIND3(is, js, ks + 3)] +
                 f4 * f3 * f6 * src[SRCIND3(is + 1, js, ks + 3)] +
                 f5 * f3 * f6 * src[SRCIND3(is + 2, js, ks + 3)] +
                 f6 * f3 * f6 * src[SRCIND3(is + 3, js, ks + 3)] +
                 f1 * f4 * f6 * src[SRCIND3(is - 2, js + 1, ks + 3)] +
                 f2 * f4 * f6 * src[SRCIND3(is - 1, js + 1, ks + 3)] +
                 f3 * f4 * f6 * src[SRCIND3(is, js + 1, ks + 3)] +
                 f4 * f4 * f6 * src[SRCIND3(is + 1, js + 1, ks + 3)] +
                 f5 * f4 * f6 * src[SRCIND3(is + 2, js + 1, ks + 3)] +
                 f6 * f4 * f6 * src[SRCIND3(is + 3, js + 1, ks + 3)] +
                 f1 * f5 * f6 * src[SRCIND3(is - 2, js + 2, ks + 3)] +
                 f2 * f5 * f6 * src[SRCIND3(is - 1, js + 2, ks + 3)] +
                 f3 * f5 * f6 * src[SRCIND3(is, js + 2, ks + 3)] +
                 f4 * f5 * f6 * src[SRCIND3(is + 1, js + 2, ks + 3)] +
                 f5 * f5 * f6 * src[SRCIND3(is + 2, js + 2, ks + 3)] +
                 f6 * f5 * f6 * src[SRCIND3(is + 3, js + 2, ks + 3)] +
                 f1 * f6 * f6 * src[SRCIND3(is - 2, js + 3, ks + 3)] +
                 f2 * f6 * f6 * src[SRCIND3(is - 1, js + 3, ks + 3)] +
                 f3 * f6 * f6 * src[SRCIND3(is, js + 3, ks + 3)] +
                 f4 * f6 * f6 * src[SRCIND3(is + 1, js + 3, ks + 3)] +
                 f5 * f6 * f6 * src[SRCIND3(is + 2, js + 3, ks + 3)] +
                 f6 * f6 * f6 * src[SRCIND3(is + 3, js + 3, ks + 3)];
  dst[DSTIND3(id, jd, kd)] = res1 + res2 + res3 + res4 + res5 + res6;
  // Monotonicity enforcement
  if ((dst[DSTIND3(id, jd, kd)] >
       max8(src[SRCIND3(is, js, ks)], src[SRCIND3(is + 1, js, ks)],
            src[SRCIND3(is, js + 1, ks)], src[SRCIND3(is, js, ks + 1)],
            src[SRCIND3(is + 1, js + 1, ks)], src[SRCIND3(is + 1, js, ks + 1)],
            src[SRCIND3(is, js + 1, ks + 1)],
            src[SRCIND3(is + 1, js + 1, ks + 1)])) ||
      (dst[DSTIND3(id, jd, kd)] <
       min8(src[SRCIND3(is, js, ks)], src[SRCIND3(is + 1, js, ks)],
            src[SRCIND3(is, js + 1, ks)], src[SRCIND3(is, js, ks + 1)],
            src[SRCIND3(is + 1, js + 1, ks)], src[SRCIND3(is + 1, js, ks + 1)],
            src[SRCIND3(is, js + 1, ks + 1)],
            src[SRCIND3(is + 1, js + 1, ks + 1)]))) {
    T const res1 = +o1_f1 * o1_f1 * o1_f1 * src[SRCIND3(is, js, ks)] +
                   o1_f2 * o1_f1 * o1_f1 * src[SRCIND3(is + 1, js, ks)] +
                   o1_f1 * o1_f2 * o1_f1 * src[SRCIND3(is, js + 1, ks)] +
                   o1_f2 * o1_f2 * o1_f1 * src[SRCIND3(is + 1, js + 1, ks)];
    T const res2 = +o1_f1 * o1_f1 * o1_f2 * src[SRCIND3(is, js, ks + 1)] +
                   o1_f2 * o1_f1 * o1_f2 * src[SRCIND3(is + 1, js, ks + 1)] +
                   o1_f1 * o1_f2 * o1_f2 * src[SRCIND3(is, js + 1, ks + 1)] +
                   o1_f2 * o1_f2 * o1_f2 * src[SRCIND3(is + 1, js + 1, ks + 1)];
    dst[DSTIND3(id, jd, kd)] = res1 + res2;
  }
}
  i = i + 1;
  id = id + 1;
  is = is + 1;
  if (i < regiext)
    goto l8110;
  goto l911;

// end i loop
l911:
  j = j + 1;
  jd = jd + 1;
  js = js + 1;
  if (j < regjext)
    goto l810;
  goto l91;

// end j loop
l91:
  k = k + 1;
  kd = kd + 1;
  ks = ks + 1;
  if (k < regkext)
    goto l80;
  goto l9;

// end k loop
l9:;
}

template <>
void prolongate_3d_o5_monotone_rf2(
    CCTK_COMPLEX const *restrict const src, ivect3 const &restrict srcpadext,
    ivect3 const &restrict srcext, CCTK_COMPLEX *restrict const dst,
    ivect3 const &restrict dstpadext, ivect3 const &restrict dstext,
    ibbox3 const &restrict srcbbox, ibbox3 const &restrict dstbbox,
    ibbox3 const &restrict, ibbox3 const &restrict regbbox, void *extraargs) {
  CCTK_WARN(0, "This should never be called!");
}

#define TYPECASE(N, T)                                                         \
  template void prolongate_3d_o5_monotone_rf2(                                 \
      T const *restrict const src, ivect3 const &restrict srcpadext,           \
      ivect3 const &restrict srcext, T *restrict const dst,                    \
      ivect3 const &restrict dstpadext, ivect3 const &restrict dstext,         \
      ibbox3 const &restrict srcbbox, ibbox3 const &restrict dstbbox,          \
      ibbox3 const &restrict, ibbox3 const &restrict regbbox,                  \
      void *extraargs);
#define CARPET_NO_COMPLEX
#define CARPET_NO_INT
#include "typecase.hh"
#undef TYPECASE

} // namespace CarpetLib
