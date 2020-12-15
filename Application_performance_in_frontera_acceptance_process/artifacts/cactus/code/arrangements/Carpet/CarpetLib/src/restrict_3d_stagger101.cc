#include <cctk.h>
#include <cctk_Parameters.h>

#include <loopcontrol.h>

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
  index3(srcioff + (i), srcjoff + (j), srckoff + (k), srcipadext, srcjpadext,  \
         srckpadext, srciext, srcjext, srckext)
#define DSTIND3(i, j, k)                                                       \
  index3(dstioff + (i), dstjoff + (j), dstkoff + (k), dstipadext, dstjpadext,  \
         dstkpadext, dstiext, dstjext, dstkext)

template <typename T>
void restrict_3d_stagger101(
    T const *restrict const src, ivect3 const &restrict srcpadext,
    ivect3 const &restrict srcext, T *restrict const dst,
    ivect3 const &restrict dstpadext, ivect3 const &restrict dstext,
    ibbox3 const &restrict srcbbox, ibbox3 const &restrict dstbbox,
    ibbox3 const &restrict, ibbox3 const &restrict regbbox, void *extraargs) {
  DECLARE_CCTK_PARAMETERS;

  assert(not extraargs);

  if (any(srcbbox.stride() >= regbbox.stride() or
          dstbbox.stride() != regbbox.stride())) {
    CCTK_ERROR("Internal error: strides disagree");
  }

  if (any(reffact2 * srcbbox.stride() != dstbbox.stride())) {
    CCTK_ERROR(
        "Internal error: destination strides are not twice the source strides");
  }

  // This could be handled, but is likely to point to an error
  // elsewhere
  if (regbbox.empty()) {
    CCTK_ERROR("Internal error: region extent is empty");
  }

  if (not regbbox.expanded_for(srcbbox).expand(1, 1).is_contained_in(srcbbox) or
      not regbbox.is_contained_in(dstbbox)) {
    cerr << "srcbbox: " << srcbbox << endl
         << "dstbbox: " << dstbbox << endl
         << "regbbox: " << regbbox << endl;
    CCTK_ERROR(
        "Internal error: region extent is not contained in array extent");
  }

  if (not support_staggered_operators) {
    CCTK_ERROR("When using staggered restriction and prolongation operators "
               "you need to set CarpetLib::support_staggered_operators");
  }

  ivect3 const regext = regbbox.shape() / regbbox.stride();
  assert(all((regbbox.lower() - srcbbox.lower()) % srcbbox.stride() == 0));
  ivect3 const srcoff = (regbbox.lower() - srcbbox.lower()) / srcbbox.stride();
  assert(all((regbbox.lower() - dstbbox.lower()) % dstbbox.stride() == 0));
  ivect3 const dstoff = (regbbox.lower() - dstbbox.lower()) / dstbbox.stride();

  int const srcipadext = srcpadext[0];
  int const srcjpadext = srcpadext[1];
  int const srckpadext = srcpadext[2];

  int const dstipadext = dstpadext[0];
  int const dstjpadext = dstpadext[1];
  int const dstkpadext = dstpadext[2];

  int const srciext = srcext[0];
  int const srcjext = srcext[1];
  int const srckext = srcext[2];

  int const dstiext = dstext[0];
  int const dstjext = dstext[1];
  int const dstkext = dstext[2];

  int const regiext = regext[0];
  int const regjext = regext[1];
  int const regkext = regext[2];

  int const srcioff = srcoff[0];
  int const srcjoff = srcoff[1];
  int const srckoff = srcoff[2];

  int const dstioff = dstoff[0];
  int const dstjoff = dstoff[1];
  int const dstkoff = dstoff[2];

#include "coeffs_restrict_3d_stagger.hh"

  if (not use_loopcontrol_in_operators) {

    // Loop over coarse region
    for (int k = 0; k < regkext; ++k) {
      for (int j = 0; j < regjext; ++j) {
        for (int i = 0; i < regiext; ++i) {

          dst[DSTIND3(i, j, k)] = typeprops<T>::fromreal(0);

          if (ORDER_STAG == 2) {
            for (int ii = 1; ii <= ORDER_STAG + 1; ii++)
              for (int jj = 1; jj <= ORDER_STAG + 1; jj++) {
                dst[DSTIND3(i, j, k)] +=
                    coeff[ii] * coeff[jj] *
                    src[SRCIND3(2 * i - 1 + (ii - 1), 2 * j,
                                2 * k - 1 + (jj - 1))];
              }
          }

          if (ORDER_STAG == 3) {
            for (int ii = 1; ii <= ORDER_STAG + 1; ii++)
              for (int jj = 1; jj <= ORDER_STAG; jj++)
                for (int kk = 1; kk <= ORDER_STAG + 1; kk++) {
                  dst[DSTIND3(i, j, k)] +=
                      coeff[ii] * coeff_i[jj] * coeff[kk] *
                      src[SRCIND3(2 * i - 1 + (ii - 1), 2 * j - 1 + (jj - 1),
                                  2 * k - 1 + (kk - 1))];
                }
          }

          if (ORDER_STAG == 4) {
            for (int ii = 1; ii <= ORDER_STAG + 1; ii++)
              for (int jj = 1; jj <= ORDER_STAG + 1; jj++)
                for (int kk = 1; kk <= ORDER_STAG + 1; kk++) {
                  dst[DSTIND3(i, j, k)] +=
                      coeff[ii] * coeff_i[jj] * coeff[kk] *
                      src[SRCIND3(2 * i - 2 + (ii - 1), 2 * j - 2 + (jj - 1),
                                  2 * k - 2 + (kk - 1))];
                }
          }
        }
      }
    }

  } else {

// Loop over coarse region
#pragma omp parallel
    CCTK_LOOP3(restrict_3d_stagger101, i, j, k, 0, 0, 0, regiext, regjext,
               regkext, dstipadext, dstjpadext, dstkpadext) {

      dst[DSTIND3(i, j, k)] = typeprops<T>::fromreal(0);

      if (ORDER_STAG == 2) {
        for (int ii = 1; ii <= ORDER_STAG + 1; ii++)
          for (int jj = 1; jj <= ORDER_STAG + 1; jj++) {
            dst[DSTIND3(i, j, k)] +=
                coeff[ii] * coeff[jj] *
                src[SRCIND3(2 * i - 1 + (ii - 1), 2 * j, 2 * k - 1 + (jj - 1))];
          }
      }

      if (ORDER_STAG == 3) {
        for (int ii = 1; ii <= ORDER_STAG + 1; ii++)
          for (int jj = 1; jj <= ORDER_STAG; jj++)
            for (int kk = 1; kk <= ORDER_STAG + 1; kk++) {
              dst[DSTIND3(i, j, k)] +=
                  coeff[ii] * coeff_i[jj] * coeff[kk] *
                  src[SRCIND3(2 * i - 1 + (ii - 1), 2 * j - 1 + (jj - 1),
                              2 * k - 1 + (kk - 1))];
            }
      }

      if (ORDER_STAG == 4) {
        for (int ii = 1; ii <= ORDER_STAG + 1; ii++)
          for (int jj = 1; jj <= ORDER_STAG + 1; jj++)
            for (int kk = 1; kk <= ORDER_STAG + 1; kk++) {
              dst[DSTIND3(i, j, k)] +=
                  coeff[ii] * coeff_i[jj] * coeff[kk] *
                  src[SRCIND3(2 * i - 2 + (ii - 1), 2 * j - 2 + (jj - 1),
                              2 * k - 2 + (kk - 1))];
            }
      }
    }
    CCTK_ENDLOOP3(restrict_3d_stagger101);
  }
}

#define TYPECASE(N, T)                                                         \
  template void restrict_3d_stagger101(                                        \
      T const *restrict const src, ivect3 const &restrict srcpadext,           \
      ivect3 const &restrict srcext, T *restrict const dst,                    \
      ivect3 const &restrict dstpadext, ivect3 const &restrict dstext,         \
      ibbox3 const &restrict srcbbox, ibbox3 const &restrict dstbbox,          \
      ibbox3 const &restrict, ibbox3 const &restrict regbbox,                  \
      void *extraargs);
#include "typecase.hh"
#undef TYPECASE

} // namespace CarpetLib
