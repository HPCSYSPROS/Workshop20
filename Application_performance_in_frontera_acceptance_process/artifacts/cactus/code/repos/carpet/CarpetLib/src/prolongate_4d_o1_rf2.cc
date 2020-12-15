#include <cctk.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>

#include "operator_prototypes_4d.hh"
#include "typeprops.hh"

using namespace std;

namespace CarpetLib {

#define SRCIND4(i, j, k, l)                                                    \
  index4(i, j, k, l, srcipadext, srcjpadext, srckpadext, srclpadext, srciext,  \
         srcjext, srckext, srclext)
#define DSTIND4(i, j, k, l)                                                    \
  index4(i, j, k, l, dstipadext, dstjpadext, dstkpadext, dstlpadext, dstiext,  \
         dstjext, dstkext, dstlext)

template <typename T>
void prolongate_4d_o1_rf2(
    T const *restrict const src, ivect4 const &restrict srcpadext,
    ivect4 const &restrict srcext, T *restrict const dst,
    ivect4 const &restrict dstpadext, ivect4 const &restrict dstext,
    ibbox4 const &restrict srcbbox, ibbox4 const &restrict dstbbox,
    ibbox4 const &restrict, ibbox4 const &restrict regbbox, void *extraargs) {
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

  ivect4 const regext = regbbox.shape() / regbbox.stride();
  assert(all((regbbox.lower() - srcbbox.lower()) % regbbox.stride() == 0));
  ivect4 const srcoff = (regbbox.lower() - srcbbox.lower()) / regbbox.stride();
  assert(all((regbbox.lower() - dstbbox.lower()) % regbbox.stride() == 0));
  ivect4 const dstoff = (regbbox.lower() - dstbbox.lower()) / regbbox.stride();

  bvect4 const needoffsetlo = srcoff % reffact2 != 0 or regext > 1;
  bvect4 const needoffsethi =
      (srcoff + regext - 1) % reffact2 != 0 or regext > 1;
  ivect4 const offsetlo = either(needoffsetlo, 1, 0);
  ivect4 const offsethi = either(needoffsethi, 1, 0);

  if (not regbbox.expand(offsetlo, offsethi).is_contained_in(srcbbox) or
      not regbbox.is_contained_in(dstbbox)) {
    CCTK_WARN(0,
              "Internal error: region extent is not contained in array extent");
  }

  size_t const srcipadext = srcpadext[0];
  size_t const srcjpadext = srcpadext[1];
  size_t const srckpadext = srcpadext[2];
  size_t const srclpadext = srcpadext[3];

  size_t const dstipadext = dstpadext[0];
  size_t const dstjpadext = dstpadext[1];
  size_t const dstkpadext = dstpadext[2];
  size_t const dstlpadext = dstpadext[3];

  size_t const srciext = srcext[0];
  size_t const srcjext = srcext[1];
  size_t const srckext = srcext[2];
  size_t const srclext = srcext[3];

  size_t const dstiext = dstext[0];
  size_t const dstjext = dstext[1];
  size_t const dstkext = dstext[2];
  size_t const dstlext = dstext[3];

  size_t const regiext = regext[0];
  size_t const regjext = regext[1];
  size_t const regkext = regext[2];
  size_t const reglext = regext[3];

  size_t const srcioff = srcoff[0];
  size_t const srcjoff = srcoff[1];
  size_t const srckoff = srcoff[2];
  size_t const srcloff = srcoff[3];

  size_t const dstioff = dstoff[0];
  size_t const dstjoff = dstoff[1];
  size_t const dstkoff = dstoff[2];
  size_t const dstloff = dstoff[3];

  size_t const fi = srcioff % 2;
  size_t const fj = srcjoff % 2;
  size_t const fk = srckoff % 2;
  size_t const fl = srcloff % 2;

  size_t const i0 = srcioff / 2;
  size_t const j0 = srcjoff / 2;
  size_t const k0 = srckoff / 2;
  size_t const l0 = srcloff / 2;

  RT const one = 1;

  RT const f1 = one / 2;
  RT const f2 = one / 2;

  // Loop over fine region
  // Label scheme: l 8 fl fk fj fi

  size_t is, js, ks, ls;
  size_t id, jd, kd, ld;
  size_t i, j, k, l;

  // begin l loop
  l = 0;
  ls = l0;
  ld = dstloff;
  if (fl == 0)
    goto l80;
  goto l81;

// begin k loop
l80:
  k = 0;
  ks = k0;
  kd = dstkoff;
  if (fk == 0)
    goto l800;
  goto l801;

// begin j loop
l800:
  j = 0;
  js = j0;
  jd = dstjoff;
  if (fj == 0)
    goto l8000;
  goto l8001;

// begin i loop
l8000:
  i = 0;
  is = i0;
  id = dstioff;
  if (fi == 0)
    goto l80000;
  goto l80001;

// kernel
l80000:
  dst[DSTIND4(id, jd, kd, ld)] = +src[SRCIND4(is, js, ks, ls)];
  i = i + 1;
  id = id + 1;
  if (i < regiext)
    goto l80001;
  goto l9000;

// kernel
l80001:
  dst[DSTIND4(id, jd, kd, ld)] = +f1 * src[SRCIND4(is, js, ks, ls)] +
                                 f2 * src[SRCIND4(is + 1, js, ks, ls)];
  i = i + 1;
  id = id + 1;
  is = is + 1;
  if (i < regiext)
    goto l80000;
  goto l9000;

// end i loop
l9000:
  j = j + 1;
  jd = jd + 1;
  if (j < regjext)
    goto l8001;
  goto l900;

// begin i loop
l8001:
  i = 0;
  is = i0;
  id = dstioff;
  if (fi == 0)
    goto l80010;
  goto l80011;

// kernel
l80010:
  dst[DSTIND4(id, jd, kd, ld)] = +f1 * src[SRCIND4(is, js, ks, ls)] +
                                 f2 * src[SRCIND4(is, js + 1, ks, ls)];
  i = i + 1;
  id = id + 1;
  if (i < regiext)
    goto l80011;
  goto l9001;

// kernel
l80011:
  dst[DSTIND4(id, jd, kd, ld)] = +f1 * f1 * src[SRCIND4(is, js, ks, ls)] +
                                 f2 * f1 * src[SRCIND4(is + 1, js, ks, ls)] +
                                 f1 * f2 * src[SRCIND4(is, js + 1, ks, ls)] +
                                 f2 * f2 * src[SRCIND4(is + 1, js + 1, ks, ls)];
  i = i + 1;
  id = id + 1;
  is = is + 1;
  if (i < regiext)
    goto l80010;
  goto l9001;

// end i loop
l9001:
  j = j + 1;
  jd = jd + 1;
  js = js + 1;
  if (j < regjext)
    goto l8000;
  goto l900;

// end j loop
l900:
  k = k + 1;
  kd = kd + 1;
  if (k < regkext)
    goto l800;
  goto l90;

// begin j loop
l801:
  j = 0;
  js = j0;
  jd = dstjoff;
  if (fj == 0)
    goto l8010;
  goto l8011;

// begin i loop
l8010:
  i = 0;
  is = i0;
  id = dstioff;
  if (fi == 0)
    goto l80100;
  goto l80101;

// kernel
l80100:
  dst[DSTIND4(id, jd, kd, ld)] = +f1 * src[SRCIND4(is, js, ks, ls)] +
                                 f2 * src[SRCIND4(is, js, ks + 1, ls)];
  i = i + 1;
  id = id + 1;
  if (i < regiext)
    goto l80101;
  goto l9010;

// kernel
l80101:
  dst[DSTIND4(id, jd, kd, ld)] = +f1 * f1 * src[SRCIND4(is, js, ks, ls)] +
                                 f2 * f1 * src[SRCIND4(is + 1, js, ks, ls)] +
                                 f1 * f2 * src[SRCIND4(is, js, ks + 1, ls)] +
                                 f2 * f2 * src[SRCIND4(is + 1, js, ks + 1, ls)];
  i = i + 1;
  id = id + 1;
  is = is + 1;
  if (i < regiext)
    goto l80100;
  goto l9010;

// end i loop
l9010:
  j = j + 1;
  jd = jd + 1;
  if (j < regjext)
    goto l8011;
  goto l901;

// begin i loop
l8011:
  i = 0;
  is = i0;
  id = dstioff;
  if (fi == 0)
    goto l80110;
  goto l80111;

// kernel
l80110:
  dst[DSTIND4(id, jd, kd, ld)] = +f1 * f1 * src[SRCIND4(is, js, ks, ls)] +
                                 f2 * f1 * src[SRCIND4(is, js + 1, ks, ls)] +
                                 f1 * f2 * src[SRCIND4(is, js, ks + 1, ls)] +
                                 f2 * f2 * src[SRCIND4(is, js + 1, ks + 1, ls)];
  i = i + 1;
  id = id + 1;
  if (i < regiext)
    goto l80111;
  goto l9011;

// kernel
l80111:
  dst[DSTIND4(id, jd, kd, ld)] =
      +f1 * f1 * f1 * src[SRCIND4(is, js, ks, ls)] +
      f2 * f1 * f1 * src[SRCIND4(is + 1, js, ks, ls)] +
      f1 * f2 * f1 * src[SRCIND4(is, js + 1, ks, ls)] +
      f2 * f2 * f1 * src[SRCIND4(is + 1, js + 1, ks, ls)] +
      f1 * f1 * f2 * src[SRCIND4(is, js, ks + 1, ls)] +
      f2 * f1 * f2 * src[SRCIND4(is + 1, js, ks + 1, ls)] +
      f1 * f2 * f2 * src[SRCIND4(is, js + 1, ks + 1, ls)] +
      f2 * f2 * f2 * src[SRCIND4(is + 1, js + 1, ks + 1, ls)];
  i = i + 1;
  id = id + 1;
  is = is + 1;
  if (i < regiext)
    goto l80110;
  goto l9011;

// end i loop
l9011:
  j = j + 1;
  jd = jd + 1;
  js = js + 1;
  if (j < regjext)
    goto l8010;
  goto l901;

// end j loop
l901:
  k = k + 1;
  kd = kd + 1;
  ks = ks + 1;
  if (k < regkext)
    goto l800;
  goto l90;

// end k loop
l90:
  l = l + 1;
  ld = ld + 1;
  ls = ls + 1;
  if (l < reglext)
    goto l81;
  goto l80;

// begin k loop
l81:
  k = 0;
  ks = k0;
  kd = dstkoff;
  if (fk == 0)
    goto l810;
  goto l811;

// begin j loop
l810:
  j = 0;
  js = j0;
  jd = dstjoff;
  if (fj == 0)
    goto l8100;
  goto l8101;

// begin i loop
l8100:
  i = 0;
  is = i0;
  id = dstioff;
  if (fi == 0)
    goto l81000;
  goto l81001;

// kernel
l81000:
  dst[DSTIND4(id, jd, kd, ld)] = +f1 * src[SRCIND4(is, js, ks, ls)] +
                                 f2 * src[SRCIND4(is, js, ks, ls + 1)];
  i = i + 1;
  id = id + 1;
  if (i < regiext)
    goto l81001;
  goto l9100;

// kernel
l81001:
  dst[DSTIND4(id, jd, kd, ld)] = +f1 * f1 * src[SRCIND4(is, js, ks, ls)] +
                                 f2 * f1 * src[SRCIND4(is + 1, js, ks, ls)] +
                                 f1 * f2 * src[SRCIND4(is, js, ks, ls + 1)] +
                                 f2 * f2 * src[SRCIND4(is + 1, js, ks, ls + 1)];
  i = i + 1;
  id = id + 1;
  is = is + 1;
  if (i < regiext)
    goto l81000;
  goto l9100;

// end i loop
l9100:
  j = j + 1;
  jd = jd + 1;
  if (j < regjext)
    goto l8101;
  goto l910;

// begin i loop
l8101:
  i = 0;
  is = i0;
  id = dstioff;
  if (fi == 0)
    goto l81010;
  goto l81011;

// kernel
l81010:
  dst[DSTIND4(id, jd, kd, ld)] = +f1 * f1 * src[SRCIND4(is, js, ks, ls)] +
                                 f2 * f1 * src[SRCIND4(is, js + 1, ks, ls)] +
                                 f1 * f2 * src[SRCIND4(is, js, ks, ls + 1)] +
                                 f2 * f2 * src[SRCIND4(is, js + 1, ks, ls + 1)];
  i = i + 1;
  id = id + 1;
  if (i < regiext)
    goto l81011;
  goto l9101;

// kernel
l81011:
  dst[DSTIND4(id, jd, kd, ld)] =
      +f1 * f1 * f1 * src[SRCIND4(is, js, ks, ls)] +
      f2 * f1 * f1 * src[SRCIND4(is + 1, js, ks, ls)] +
      f1 * f2 * f1 * src[SRCIND4(is, js + 1, ks, ls)] +
      f2 * f2 * f1 * src[SRCIND4(is + 1, js + 1, ks, ls)] +
      f1 * f1 * f2 * src[SRCIND4(is, js, ks, ls + 1)] +
      f2 * f1 * f2 * src[SRCIND4(is + 1, js, ks, ls + 1)] +
      f1 * f2 * f2 * src[SRCIND4(is, js + 1, ks, ls + 1)] +
      f2 * f2 * f2 * src[SRCIND4(is + 1, js + 1, ks, ls + 1)];
  i = i + 1;
  id = id + 1;
  is = is + 1;
  if (i < regiext)
    goto l81010;
  goto l9101;

// end i loop
l9101:
  j = j + 1;
  jd = jd + 1;
  js = js + 1;
  if (j < regjext)
    goto l8100;
  goto l910;

// end j loop
l910:
  k = k + 1;
  kd = kd + 1;
  if (k < regkext)
    goto l810;
  goto l91;

// begin j loop
l811:
  j = 0;
  js = j0;
  jd = dstjoff;
  if (fj == 0)
    goto l8110;
  goto l8111;

// begin i loop
l8110:
  i = 0;
  is = i0;
  id = dstioff;
  if (fi == 0)
    goto l81100;
  goto l81101;

// kernel
l81100:
  dst[DSTIND4(id, jd, kd, ld)] = +f1 * f1 * src[SRCIND4(is, js, ks, ls)] +
                                 f2 * f1 * src[SRCIND4(is, js, ks + 1, ls)] +
                                 f1 * f2 * src[SRCIND4(is, js, ks, ls + 1)] +
                                 f2 * f2 * src[SRCIND4(is, js, ks + 1, ls + 1)];
  i = i + 1;
  id = id + 1;
  if (i < regiext)
    goto l81101;
  goto l9110;

// kernel
l81101:
  dst[DSTIND4(id, jd, kd, ld)] =
      +f1 * f1 * f1 * src[SRCIND4(is, js, ks, ls)] +
      f2 * f1 * f1 * src[SRCIND4(is + 1, js, ks, ls)] +
      f1 * f2 * f1 * src[SRCIND4(is, js, ks + 1, ls)] +
      f2 * f2 * f1 * src[SRCIND4(is + 1, js, ks + 1, ls)] +
      f1 * f1 * f2 * src[SRCIND4(is, js, ks, ls + 1)] +
      f2 * f1 * f2 * src[SRCIND4(is + 1, js, ks, ls + 1)] +
      f1 * f2 * f2 * src[SRCIND4(is, js, ks + 1, ls + 1)] +
      f2 * f2 * f2 * src[SRCIND4(is + 1, js, ks + 1, ls + 1)];
  i = i + 1;
  id = id + 1;
  is = is + 1;
  if (i < regiext)
    goto l81100;
  goto l9110;

// end i loop
l9110:
  j = j + 1;
  jd = jd + 1;
  if (j < regjext)
    goto l8111;
  goto l911;

// begin i loop
l8111:
  i = 0;
  is = i0;
  id = dstioff;
  if (fi == 0)
    goto l81110;
  goto l81111;

// kernel
l81110:
  dst[DSTIND4(id, jd, kd, ld)] =
      +f1 * f1 * f1 * f1 * src[SRCIND4(is, js, ks, ls)] +
      f2 * f1 * f1 * f1 * src[SRCIND4(is, js + 1, ks, ls)] +
      f1 * f2 * f2 * f1 * src[SRCIND4(is, js, ks + 1, ls)] +
      f2 * f2 * f2 * f1 * src[SRCIND4(is, js + 1, ks + 1, ls)] +
      f1 * f1 * f1 * f2 * src[SRCIND4(is, js, ks, ls + 1)] +
      f2 * f1 * f1 * f2 * src[SRCIND4(is, js + 1, ks, ls + 1)] +
      f1 * f2 * f2 * f2 * src[SRCIND4(is, js, ks + 1, ls + 1)] +
      f2 * f2 * f2 * f2 * src[SRCIND4(is, js + 1, ks + 1, ls + 1)];
  i = i + 1;
  id = id + 1;
  if (i < regiext)
    goto l81111;
  goto l9111;

// kernel
l81111:
  dst[DSTIND4(id, jd, kd, ld)] =
      +f1 * f1 * f1 * f1 * src[SRCIND4(is, js, ks, ls)] +
      f2 * f1 * f1 * f1 * src[SRCIND4(is + 1, js, ks, ls)] +
      f1 * f2 * f1 * f1 * src[SRCIND4(is, js + 1, ks, ls)] +
      f2 * f2 * f1 * f1 * src[SRCIND4(is + 1, js + 1, ks, ls)] +
      f1 * f1 * f2 * f1 * src[SRCIND4(is, js, ks + 1, ls)] +
      f2 * f1 * f2 * f1 * src[SRCIND4(is + 1, js, ks + 1, ls)] +
      f1 * f2 * f2 * f1 * src[SRCIND4(is, js + 1, ks + 1, ls)] +
      f2 * f2 * f2 * f1 * src[SRCIND4(is + 1, js + 1, ks + 1, ls)] +
      f1 * f1 * f1 * f2 * src[SRCIND4(is, js, ks, ls + 1)] +
      f2 * f1 * f1 * f2 * src[SRCIND4(is + 1, js, ks, ls + 1)] +
      f1 * f2 * f1 * f2 * src[SRCIND4(is, js + 1, ks, ls + 1)] +
      f2 * f2 * f1 * f2 * src[SRCIND4(is + 1, js + 1, ks, ls + 1)] +
      f1 * f1 * f2 * f2 * src[SRCIND4(is, js, ks + 1, ls + 1)] +
      f2 * f1 * f2 * f2 * src[SRCIND4(is + 1, js, ks + 1, ls + 1)] +
      f1 * f2 * f2 * f2 * src[SRCIND4(is, js + 1, ks + 1, ls + 1)] +
      f2 * f2 * f2 * f2 * src[SRCIND4(is + 1, js + 1, ks + 1, ls + 1)];
  i = i + 1;
  id = id + 1;
  is = is + 1;
  if (i < regiext)
    goto l81110;
  goto l9111;

// end i loop
l9111:
  j = j + 1;
  jd = jd + 1;
  js = js + 1;
  if (j < regjext)
    goto l8110;
  goto l911;

// end j loop
l911:
  k = k + 1;
  kd = kd + 1;
  ks = ks + 1;
  if (k < regkext)
    goto l810;
  goto l91;

// end k loop
l91:
  l = l + 1;
  ld = ld + 1;
  ls = ls + 1;
  if (l < reglext)
    goto l81;
  goto l81;
}

#define TYPECASE(N, T)                                                         \
  template void prolongate_4d_o1_rf2(                                          \
      T const *restrict const src, ivect4 const &restrict srcpadext,           \
      ivect4 const &restrict srcext, T *restrict const dst,                    \
      ivect4 const &restrict dstpadext, ivect4 const &restrict dstext,         \
      ibbox4 const &restrict srcbbox, ibbox4 const &restrict dstbbox,          \
      ibbox4 const &restrict, ibbox4 const &restrict regbbox,                  \
      void *extraargs);
#include "typecase.hh"
#undef TYPECASE

} // CarpetLib
