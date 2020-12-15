#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>

#include "operator_prototypes_3d.hh"
#include "typeprops.hh"

using namespace std;

namespace CarpetLib {

#define SRCIND3(i, j, k)                                                       \
  index3(i, j, k, srcipadext, srcjpadext, srckpadext, srciext, srcjext, srckext)
#define DSTIND3(i, j, k)                                                       \
  index3(i, j, k, dstipadext, dstjpadext, dstkpadext, dstiext, dstjext, dstkext)
#define SRCOFF3(i, j, k) offset3(i, j, k, srciext, srcjext, srckext)
#define DSTOFF3(i, j, k) offset3(i, j, k, dstiext, dstjext, dstkext)

template <typename T, int ORDER>
void prolongate_3d_stagger011(
    T const *restrict const src, ivect3 const &restrict srcpadext,
    ivect3 const &restrict srcext, T *restrict const dst,
    ivect3 const &restrict dstpadext, ivect3 const &restrict dstext,
    ibbox3 const &restrict srcbbox, ibbox3 const &restrict dstbbox,
    ibbox3 const &restrict, ibbox3 const &restrict regbbox, void *extraargs) {
  DECLARE_CCTK_PARAMETERS;

  assert(not extraargs);

  if (any(srcbbox.stride() <= regbbox.stride() or
          dstbbox.stride() != regbbox.stride())) {
    CCTK_ERROR("Internal error: strides disagree");
  }

  if (any(srcbbox.stride() != reffact2 * dstbbox.stride())) {
    CCTK_ERROR(
        "Internal error: source strides are not twice the destination strides");
  }

  if (any(srcbbox.lower() % srcbbox.stride() != 0)) {
    CCTK_ERROR("Internal error: source bbox is not aligned with vertices");
  }
  if (any(dstbbox.lower() % dstbbox.stride() != 0)) {
    CCTK_ERROR("Internal error: destination bbox is not aligned with vertices");
  }
  if (any(regbbox.lower() % regbbox.stride() != 0)) {
    CCTK_ERROR("Internal error: prolongation region bbox is not aligned with "
               "vertices");
  }

  // This could be handled, but is likely to point to an error
  // elsewhere
  if (regbbox.empty()) {
    CCTK_ERROR("Internal error: region extent is empty");
  }

  if (not support_staggered_operators) {
    CCTK_ERROR("When using staggered restriction and prolongation operators "
               "you need to set CarpetLib::support_staggered_operators");
  }

  ivect3 const regext = regbbox.shape() / regbbox.stride();
  assert(all((regbbox.lower() - srcbbox.lower()) % regbbox.stride() == 0));
  ivect3 const srcoff = (regbbox.lower() - srcbbox.lower()) / regbbox.stride();
  assert(all((regbbox.lower() - dstbbox.lower()) % regbbox.stride() == 0));
  ivect3 const dstoff = (regbbox.lower() - dstbbox.lower()) / regbbox.stride();

  bvect3 const needoffsetlo = srcoff % reffact2 != 0;
  bvect3 const needoffsethi = (srcoff + regext - 1) % reffact2 != 0;
  ivect3 const offsetlo =
      either(needoffsetlo, ORDER / 2 + 1, either(regext > 1, ORDER / 2, 0));
  ivect3 const offsethi =
      either(needoffsethi, ORDER / 2 + 1, either(regext > 1, ORDER / 2, 0));

  if (not regbbox.expand(offsetlo, offsethi)
              .expanded_for(srcbbox)
              .expand(1, 1)
              .is_contained_in(srcbbox) or
      not regbbox.is_contained_in(dstbbox)) {
    cerr << "ORDER=" << ORDER << "\n"
         << "offsetlo=" << offsetlo << "\n"
         << "offsethi=" << offsethi << "\n"
         << "regbbox=" << regbbox << "\n"
         << "dstbbox=" << dstbbox << "\n"
         << "regbbox.expand=" << regbbox.expand(offsetlo, offsethi) << "\n"
         << "srcbbox=" << srcbbox << "\n";
    CCTK_ERROR(
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

#include "coeffs_prolongate_3d_stagger.hh"

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
  // Staggered code, where source and destination are off by 1/4 of a gridpoint,
  // in y and z directions:
  // Interpolate to (is, js-1/4, ks-1/4), where is,js,ks is a point on the
  // COARSE grid, since y and z are staggered
  dst[DSTIND3(id, jd, kd)] = typeprops<T>::fromreal(0);
  if (ORDER <= 2) {
    for (int ii = 1; ii <= ORDER + 1; ii++)
      for (int jj = 1; jj <= ORDER + 1; jj++) {
        dst[DSTIND3(id, jd, kd)] +=
            coeff[1][ii] * coeff[1][jj] *
            src[SRCIND3(is, js - 1 + (ii - 1), ks - 1 + (jj - 1))];
      }
  }
  if (ORDER == 3) {
    for (int ii = 1; ii <= ORDER; ii++)
      for (int jj = 1; jj <= ORDER + 1; jj++)
        for (int kk = 1; kk <= ORDER + 1; kk++) {
          dst[DSTIND3(id, jd, kd)] +=
              coeff[4][ii] * coeff[1][jj] * coeff[1][kk] *
              src[SRCIND3(is - 1 + (ii - 1), js - 2 + (jj - 1),
                          ks - 2 + (kk - 1))];
        }
  }
  if (ORDER == 4) {
    for (int ii = 1; ii <= ORDER + 1; ii++)
      for (int jj = 1; jj <= ORDER + 1; jj++)
        for (int kk = 1; kk <= ORDER + 1; kk++) {
          dst[DSTIND3(id, jd, kd)] +=
              coeff[4][ii] * coeff[1][jj] * coeff[1][kk] *
              src[SRCIND3(is - 2 + (ii - 1), js - 2 + (jj - 1),
                          ks - 2 + (kk - 1))];
        }
  }

  i = i + 1;
  id = id + 1;
  if (i < regiext)
    goto l8001;
  goto l900;

// kernel
l8001:
  // Staggered code, where source and destination are off by 1/4 of a gridpoint,
  // in y and z directions:
  // Interpolate to (is+1/2,js-1/4,ks-1/4), where is,js,ks is a point on the
  // COARSE grid, since y and z are staggered
  dst[DSTIND3(id, jd, kd)] = typeprops<T>::fromreal(0);
  if (ORDER <= 2) {
    for (int ii = 1; ii <= ORDER + 1; ii++)
      for (int jj = 1; jj <= ORDER + 1; jj++)
        for (int kk = 1; kk <= ORDER + 1; kk++) {
          dst[DSTIND3(id, jd, kd)] +=
              coeff[3][ii] * coeff[1][jj] * coeff[1][kk] *
              src[SRCIND3(is - 1 + (ii - 1), js - 1 + (jj - 1),
                          ks - 1 + (kk - 1))];
        }
  }
  if (ORDER == 3 || ORDER == 4) {
    for (int ii = 1; ii <= ORDER + 1; ii++)
      for (int jj = 1; jj <= ORDER + 1; jj++)
        for (int kk = 1; kk <= ORDER + 1; kk++) {
          dst[DSTIND3(id, jd, kd)] +=
              coeff[3][ii] * coeff[1][jj] * coeff[1][kk] *
              src[SRCIND3(is - 2 + (ii - 1) + (4 - ORDER), js - 2 + (jj - 1),
                          ks - 2 + (kk - 1))];
        }
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
  // Staggered code, where source and destination are off by 1/4 of a gridpoint,
  // in y and z directions:
  // Interpolate to (is,js+1/4,ks-1/4)
  dst[DSTIND3(id, jd, kd)] = typeprops<T>::fromreal(0);
  if (ORDER <= 2) {
    for (int ii = 1; ii <= ORDER + 1; ii++)
      for (int jj = 1; jj <= ORDER + 1; jj++) {
        dst[DSTIND3(id, jd, kd)] +=
            coeff[2][ii] * coeff[1][jj] *
            src[SRCIND3(is, js - 1 + (ii - 1), ks - 1 + (jj - 1))];
      }
  }
  if (ORDER == 3) {
    for (int ii = 1; ii <= ORDER; ii++)
      for (int jj = 1; jj <= ORDER + 1; jj++)
        for (int kk = 1; kk <= ORDER + 1; kk++) {
          dst[DSTIND3(id, jd, kd)] +=
              coeff[4][ii] * coeff[2][jj] * coeff[1][kk] *
              src[SRCIND3(is - 1 + (ii - 1), js - 2 + (jj - 1) + (4 - ORDER),
                          ks - 2 + (kk - 1))];
        }
  }
  if (ORDER == 4) {
    for (int ii = 1; ii <= ORDER + 1; ii++)
      for (int jj = 1; jj <= ORDER + 1; jj++)
        for (int kk = 1; kk <= ORDER + 1; kk++) {
          dst[DSTIND3(id, jd, kd)] +=
              coeff[4][ii] * coeff[2][jj] * coeff[1][kk] *
              src[SRCIND3(is - 2 + (ii - 1), js - 2 + (jj - 1) + (4 - ORDER),
                          ks - 2 + (kk - 1))];
        }
  }

  i = i + 1;
  id = id + 1;
  if (i < regiext)
    goto l8011;
  goto l901;

// kernel
l8011:
  // Staggered code, where source and destination are off by 1/4 of a gridpoint,
  // in y and z directions:
  // Interpolate to (is+1/2,js+1/4,ks-1/4)
  dst[DSTIND3(id, jd, kd)] = typeprops<T>::fromreal(0);
  if (ORDER <= 2) {
    for (int ii = 1; ii <= ORDER + 1; ii++)
      for (int jj = 1; jj <= ORDER + 1; jj++)
        for (int kk = 1; kk <= ORDER + 1; kk++) {
          dst[DSTIND3(id, jd, kd)] +=
              coeff[3][ii] * coeff[2][jj] * coeff[1][kk] *
              src[SRCIND3(is - 1 + (ii - 1), js - 1 + (jj - 1),
                          ks - 1 + (kk - 1))];
        }
  }
  if (ORDER == 3 || ORDER == 4) {
    for (int ii = 1; ii <= ORDER + 1; ii++)
      for (int jj = 1; jj <= ORDER + 1; jj++)
        for (int kk = 1; kk <= ORDER + 1; kk++) {
          dst[DSTIND3(id, jd, kd)] +=
              coeff[3][ii] * coeff[2][jj] * coeff[1][kk] *
              src[SRCIND3(is - 2 + (ii - 1) + (4 - ORDER),
                          js - 2 + (jj - 1) + (4 - ORDER), ks - 2 + (kk - 1))];
        }
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
  // Staggered code, where source and destination are off by 1/4 of a gridpoint,
  // in y and z directions:
  // Interpolate to (is,js-1/4,ks+1/4)
  dst[DSTIND3(id, jd, kd)] = typeprops<T>::fromreal(0);
  if (ORDER <= 2) {
    for (int ii = 1; ii <= ORDER + 1; ii++)
      for (int jj = 1; jj <= ORDER + 1; jj++) {
        dst[DSTIND3(id, jd, kd)] +=
            coeff[1][ii] * coeff[2][jj] *
            src[SRCIND3(is, js - 1 + (ii - 1), ks - 1 + (jj - 1))];
      }
  }
  if (ORDER == 3) {
    for (int ii = 1; ii <= ORDER; ii++)
      for (int jj = 1; jj <= ORDER + 1; jj++)
        for (int kk = 1; kk <= ORDER + 1; kk++) {
          dst[DSTIND3(id, jd, kd)] +=
              coeff[4][ii] * coeff[1][jj] * coeff[2][kk] *
              src[SRCIND3(is - 1 + (ii - 1), js - 2 + (jj - 1),
                          ks - 2 + (kk - 1) + (4 - ORDER))];
        }
  }
  if (ORDER == 4) {
    for (int ii = 1; ii <= ORDER + 1; ii++)
      for (int jj = 1; jj <= ORDER + 1; jj++)
        for (int kk = 1; kk <= ORDER + 1; kk++) {
          dst[DSTIND3(id, jd, kd)] +=
              coeff[4][ii] * coeff[1][jj] * coeff[2][kk] *
              src[SRCIND3(is - 2 + (ii - 1), js - 2 + (jj - 1),
                          ks - 2 + (kk - 1) + (4 - ORDER))];
        }
  }

  i = i + 1;
  id = id + 1;
  if (i < regiext)
    goto l8101;
  goto l910;

// kernel
l8101:
  // Staggered code, where source and destination are off by 1/4 of a gridpoint,
  // in y and z directions:
  // Interpolate to (is+1/2,js-1/4,ks+1/4)
  dst[DSTIND3(id, jd, kd)] = typeprops<T>::fromreal(0);
  if (ORDER <= 2) {
    for (int ii = 1; ii <= ORDER + 1; ii++)
      for (int jj = 1; jj <= ORDER + 1; jj++)
        for (int kk = 1; kk <= ORDER + 1; kk++) {
          dst[DSTIND3(id, jd, kd)] +=
              coeff[3][ii] * coeff[1][jj] * coeff[2][kk] *
              src[SRCIND3(is - 1 + (ii - 1), js - 1 + (jj - 1),
                          ks - 1 + (kk - 1))];
        }
  }
  if (ORDER == 3 || ORDER == 4) {
    for (int ii = 1; ii <= ORDER + 1; ii++)
      for (int jj = 1; jj <= ORDER + 1; jj++)
        for (int kk = 1; kk <= ORDER + 1; kk++) {
          dst[DSTIND3(id, jd, kd)] +=
              coeff[3][ii] * coeff[1][jj] * coeff[2][kk] *
              src[SRCIND3(is - 2 + (ii - 1) + (4 - ORDER), js - 2 + (jj - 1),
                          ks - 2 + (kk - 1) + (4 - ORDER))];
        }
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
  // Staggered code, where source and destination are off by 1/4 of a gridpoint,
  // in y and z directions:
  // Interpolate to (is,js+1/4,ks+1/4)
  dst[DSTIND3(id, jd, kd)] = typeprops<T>::fromreal(0);
  if (ORDER <= 2) {
    for (int ii = 1; ii <= ORDER + 1; ii++)
      for (int jj = 1; jj <= ORDER + 1; jj++) {
        dst[DSTIND3(id, jd, kd)] +=
            coeff[2][ii] * coeff[2][jj] *
            src[SRCIND3(is, js - 1 + (ii - 1), ks - 1 + (jj - 1))];
      }
  }
  if (ORDER == 3) {
    for (int ii = 1; ii <= ORDER; ii++)
      for (int jj = 1; jj <= ORDER + 1; jj++)
        for (int kk = 1; kk <= ORDER + 1; kk++) {
          dst[DSTIND3(id, jd, kd)] +=
              coeff[4][ii] * coeff[2][jj] * coeff[2][kk] *
              src[SRCIND3(is - 1 + (ii - 1), js - 2 + (jj - 1) + (4 - ORDER),
                          ks - 2 + (kk - 1) + (4 - ORDER))];
        }
  }
  if (ORDER == 4) {
    for (int ii = 1; ii <= ORDER + 1; ii++)
      for (int jj = 1; jj <= ORDER + 1; jj++)
        for (int kk = 1; kk <= ORDER + 1; kk++) {
          dst[DSTIND3(id, jd, kd)] +=
              coeff[4][ii] * coeff[2][jj] * coeff[2][kk] *
              src[SRCIND3(is - 2 + (ii - 1), js - 2 + (jj - 1) + (4 - ORDER),
                          ks - 2 + (kk - 1) + (4 - ORDER))];
        }
  }

  i = i + 1;
  id = id + 1;
  if (i < regiext)
    goto l8111;
  goto l911;

// kernel
l8111:
  // Staggered code, where source and destination are off by 1/4 of a gridpoint,
  // in y and z directions:
  // Interpolate to (is+1/2,js+1/4,ks+1/4)
  dst[DSTIND3(id, jd, kd)] = typeprops<T>::fromreal(0);
  if (ORDER <= 2) {
    for (int ii = 1; ii <= ORDER + 1; ii++)
      for (int jj = 1; jj <= ORDER + 1; jj++)
        for (int kk = 1; kk <= ORDER + 1; kk++) {
          dst[DSTIND3(id, jd, kd)] +=
              coeff[3][ii] * coeff[2][jj] * coeff[2][kk] *
              src[SRCIND3(is - 1 + (ii - 1), js - 1 + (jj - 1),
                          ks - 1 + (kk - 1))];
        }
  }
  if (ORDER == 3 || ORDER == 4) {
    for (int ii = 1; ii <= ORDER + 1; ii++)
      for (int jj = 1; jj <= ORDER + 1; jj++)
        for (int kk = 1; kk <= ORDER + 1; kk++) {
          dst[DSTIND3(id, jd, kd)] +=
              coeff[3][ii] * coeff[2][jj] * coeff[2][kk] *
              src[SRCIND3(is - 2 + (ii - 1) + (4 - ORDER),
                          js - 2 + (jj - 1) + (4 - ORDER),
                          ks - 2 + (kk - 1) + (4 - ORDER))];
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

#define TYPECASE(N, T)                                                         \
                                                                               \
  template void prolongate_3d_stagger011<T, 2>(                                \
      T const *restrict const src, ivect3 const &restrict srcpadext,           \
      ivect3 const &restrict srcext, T *restrict const dst,                    \
      ivect3 const &restrict dstpadext, ivect3 const &restrict dstext,         \
      ibbox3 const &restrict srcbbox, ibbox3 const &restrict dstbbox,          \
      ibbox3 const &restrict, ibbox3 const &restrict regbbox,                  \
      void *extraargs);                                                        \
                                                                               \
  template void prolongate_3d_stagger011<T, 3>(                                \
      T const *restrict const src, ivect3 const &restrict srcpadext,           \
      ivect3 const &restrict srcext, T *restrict const dst,                    \
      ivect3 const &restrict dstpadext, ivect3 const &restrict dstext,         \
      ibbox3 const &restrict srcbbox, ibbox3 const &restrict dstbbox,          \
      ibbox3 const &restrict, ibbox3 const &restrict regbbox,                  \
      void *extraargs);                                                        \
                                                                               \
  template void prolongate_3d_stagger011<T, 4>(                                \
      T const *restrict const src, ivect3 const &restrict srcpadext,           \
      ivect3 const &restrict srcext, T *restrict const dst,                    \
      ivect3 const &restrict dstpadext, ivect3 const &restrict dstext,         \
      ibbox3 const &restrict srcbbox, ibbox3 const &restrict dstbbox,          \
      ibbox3 const &restrict, ibbox3 const &restrict regbbox,                  \
      void *extraargs);

#define CARPET_NO_INT
#include "typecase.hh"
#undef TYPECASE

} // namespace CarpetLib
