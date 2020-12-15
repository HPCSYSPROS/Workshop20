#include <cctk.h>
#include <cctk_Parameters.h>

#include <loopcontrol.h>

#include <algorithm>
#include <cassert>
#include <cmath>

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

// This operator offers third-order accurate restriction operators for cell
// centered grids when use_higher_order_restriction is set. This interpolation
// is done for samples at the cell centres (so this is not a ppm scheme or
// anything like that).  It really only fits a polynomial of the form
//    f(x,y,z) = \sum_{i,j,k=0}^3 a_{i,j,k} x^i y^j z^k
// to the fine cells and evaluates at x=y=z=0. So it is good for the metric,
// but bad for matter (since it will destroy the conservation).  Because of
// this it should not be used for grid functions whose transport operator is
// not WENO or ENO which hopefully excludes all matter variables.

template <typename T>
void restrict_3d_cc_o3_rf2(
    T const *restrict const src, ivect3 const &restrict srcpadext,
    ivect3 const &restrict srcext, T *restrict const dst,
    ivect3 const &restrict dstpadext, ivect3 const &restrict dstext,
    ibbox3 const &restrict srcbbox, ibbox3 const &restrict dstbbox,
    ibbox3 const &restrict, ibbox3 const &restrict regbbox, void *extraargs) {
  assert(not extraargs);

  DECLARE_CCTK_PARAMETERS;

  typedef typename typeprops<T>::real RT;

  if (any(srcbbox.stride() >= regbbox.stride() or
          dstbbox.stride() != regbbox.stride())) {
    CCTK_WARN(0, "Internal error: strides disagree");
  }

  if (any(reffact2 * srcbbox.stride() != dstbbox.stride())) {
    CCTK_WARN(
        0,
        "Internal error: destination strides are not twice the source strides");
  }

  // This could be handled, but is likely to point to an error
  // elsewhere
  if (regbbox.empty()) {
    CCTK_WARN(0, "Internal error: region extent is empty");
  }

  if (not regbbox.expanded_for(srcbbox).expand(1, 1).is_contained_in(srcbbox) or
      not regbbox.is_contained_in(dstbbox)) {
    cerr << "srcbbox: " << srcbbox << endl
         << "dstbbox: " << dstbbox << endl
         << "regbbox: " << regbbox << endl;
    CCTK_WARN(0,
              "Internal error: region extent is not contained in array extent");
  }

  ivect3 const regext = regbbox.shape() / regbbox.stride();
  assert(all(srcbbox.stride() % 2 == 0));
  assert(all((regbbox.lower() - srcbbox.lower() - srcbbox.stride() / 2) %
                 srcbbox.stride() ==
             0));
  ivect3 const srcoff =
      (regbbox.lower() - srcbbox.lower() - srcbbox.stride() / 2) /
      srcbbox.stride();
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

  RT const den = 4096;

  RT const f1 = 1 / den;
  RT const f2 = 9 / den;
  RT const f3 = 81 / den;
  RT const f4 = 729 / den;

  if (not use_loopcontrol_in_operators) {

// Loop over coarse region
#pragma omp parallel for collapse(3)
    for (int k = 0; k < regkext; ++k) {
      for (int j = 0; j < regjext; ++j) {
        for (int i = 0; i < regiext; ++i) {

#ifdef CARPET_DEBUG
          if (not(2 * k + 2 + srckoff < srckext and
                  2 * j + 2 + srcjoff < srcjext and
                  2 * i + 2 + srcioff < srciext and 2 * k - 1 + srckoff >= 0 and
                  2 * j - 1 + srcjoff >= 0 and 2 * i - 1 + srcioff >= 0)) {
            cout << "restrict_3d_cc_o3_rf2.cc\n";
            cout << "regext " << regext << "\n";
            cout << "srcext " << srcext << "\n";
            cout << "srcbbox=" << srcbbox << "\n";
            cout << "dstbbox=" << dstbbox << "\n";
            cout << "regbbox=" << regbbox << "\n";
            cout << "i,j,k=" << i << " " << j << " " << k << "\n";
            assert(2 * k + 2 + srckoff < srckext);
            assert(2 * j + 2 + srcjoff < srcjext);
            assert(2 * i + 2 + srcioff < srciext);
            assert(2 * k - 1 + srckoff >= 0);
            assert(2 * j - 1 + srcjoff >= 0);
            assert(2 * i - 1 + srcioff >= 0);
          }
#endif
          dst[DSTIND3(i, j, k)] =
              -f1 * src[SRCIND3(2 * i - 1, 2 * j - 1, 2 * k - 1)] +
              f2 * src[SRCIND3(2 * i, 2 * j - 1, 2 * k - 1)] +
              f2 * src[SRCIND3(2 * i + 1, 2 * j - 1, 2 * k - 1)] -
              f1 * src[SRCIND3(2 * i + 2, 2 * j - 1, 2 * k - 1)] +
              f2 * src[SRCIND3(2 * i - 1, 2 * j, 2 * k - 1)] -
              f3 * src[SRCIND3(2 * i, 2 * j, 2 * k - 1)] -
              f3 * src[SRCIND3(2 * i + 1, 2 * j, 2 * k - 1)] +
              f2 * src[SRCIND3(2 * i + 2, 2 * j, 2 * k - 1)] +
              f2 * src[SRCIND3(2 * i - 1, 2 * j + 1, 2 * k - 1)] -
              f3 * src[SRCIND3(2 * i, 2 * j + 1, 2 * k - 1)] -
              f3 * src[SRCIND3(2 * i + 1, 2 * j + 1, 2 * k - 1)] +
              f2 * src[SRCIND3(2 * i + 2, 2 * j + 1, 2 * k - 1)] -
              f1 * src[SRCIND3(2 * i - 1, 2 * j + 2, 2 * k - 1)] +
              f2 * src[SRCIND3(2 * i, 2 * j + 2, 2 * k - 1)] +
              f2 * src[SRCIND3(2 * i + 1, 2 * j + 2, 2 * k - 1)] -
              f1 * src[SRCIND3(2 * i + 2, 2 * j + 2, 2 * k - 1)] +
              f2 * src[SRCIND3(2 * i - 1, 2 * j - 1, 2 * k)] -
              f3 * src[SRCIND3(2 * i, 2 * j - 1, 2 * k)] -
              f3 * src[SRCIND3(2 * i + 1, 2 * j - 1, 2 * k)] +
              f2 * src[SRCIND3(2 * i + 2, 2 * j - 1, 2 * k)] -
              f3 * src[SRCIND3(2 * i - 1, 2 * j, 2 * k)] +
              f4 * src[SRCIND3(2 * i, 2 * j, 2 * k)] +
              f4 * src[SRCIND3(2 * i + 1, 2 * j, 2 * k)] -
              f3 * src[SRCIND3(2 * i + 2, 2 * j, 2 * k)] -
              f3 * src[SRCIND3(2 * i - 1, 2 * j + 1, 2 * k)] +
              f4 * src[SRCIND3(2 * i, 2 * j + 1, 2 * k)] +
              f4 * src[SRCIND3(2 * i + 1, 2 * j + 1, 2 * k)] -
              f3 * src[SRCIND3(2 * i + 2, 2 * j + 1, 2 * k)] +
              f2 * src[SRCIND3(2 * i - 1, 2 * j + 2, 2 * k)] -
              f3 * src[SRCIND3(2 * i, 2 * j + 2, 2 * k)] -
              f3 * src[SRCIND3(2 * i + 1, 2 * j + 2, 2 * k)] +
              f2 * src[SRCIND3(2 * i + 2, 2 * j + 2, 2 * k)] +
              f2 * src[SRCIND3(2 * i - 1, 2 * j - 1, 2 * k + 1)] -
              f3 * src[SRCIND3(2 * i, 2 * j - 1, 2 * k + 1)] -
              f3 * src[SRCIND3(2 * i + 1, 2 * j - 1, 2 * k + 1)] +
              f2 * src[SRCIND3(2 * i + 2, 2 * j - 1, 2 * k + 1)] -
              f3 * src[SRCIND3(2 * i - 1, 2 * j, 2 * k + 1)] +
              f4 * src[SRCIND3(2 * i, 2 * j, 2 * k + 1)] +
              f4 * src[SRCIND3(2 * i + 1, 2 * j, 2 * k + 1)] -
              f3 * src[SRCIND3(2 * i + 2, 2 * j, 2 * k + 1)] -
              f3 * src[SRCIND3(2 * i - 1, 2 * j + 1, 2 * k + 1)] +
              f4 * src[SRCIND3(2 * i, 2 * j + 1, 2 * k + 1)] +
              f4 * src[SRCIND3(2 * i + 1, 2 * j + 1, 2 * k + 1)] -
              f3 * src[SRCIND3(2 * i + 2, 2 * j + 1, 2 * k + 1)] +
              f2 * src[SRCIND3(2 * i - 1, 2 * j + 2, 2 * k + 1)] -
              f3 * src[SRCIND3(2 * i, 2 * j + 2, 2 * k + 1)] -
              f3 * src[SRCIND3(2 * i + 1, 2 * j + 2, 2 * k + 1)] +
              f2 * src[SRCIND3(2 * i + 2, 2 * j + 2, 2 * k + 1)] -
              f1 * src[SRCIND3(2 * i - 1, 2 * j - 1, 2 * k + 2)] +
              f2 * src[SRCIND3(2 * i, 2 * j - 1, 2 * k + 2)] +
              f2 * src[SRCIND3(2 * i + 1, 2 * j - 1, 2 * k + 2)] -
              f1 * src[SRCIND3(2 * i + 2, 2 * j - 1, 2 * k + 2)] +
              f2 * src[SRCIND3(2 * i - 1, 2 * j, 2 * k + 2)] -
              f3 * src[SRCIND3(2 * i, 2 * j, 2 * k + 2)] -
              f3 * src[SRCIND3(2 * i + 1, 2 * j, 2 * k + 2)] +
              f2 * src[SRCIND3(2 * i + 2, 2 * j, 2 * k + 2)] +
              f2 * src[SRCIND3(2 * i - 1, 2 * j + 1, 2 * k + 2)] -
              f3 * src[SRCIND3(2 * i, 2 * j + 1, 2 * k + 2)] -
              f3 * src[SRCIND3(2 * i + 1, 2 * j + 1, 2 * k + 2)] +
              f2 * src[SRCIND3(2 * i + 2, 2 * j + 1, 2 * k + 2)] -
              f1 * src[SRCIND3(2 * i - 1, 2 * j + 2, 2 * k + 2)] +
              f2 * src[SRCIND3(2 * i, 2 * j + 2, 2 * k + 2)] +
              f2 * src[SRCIND3(2 * i + 1, 2 * j + 2, 2 * k + 2)] -
              f1 * src[SRCIND3(2 * i + 2, 2 * j + 2, 2 * k + 2)];
        }
      }
    }

  } else {

// Loop over coarse region
#pragma omp parallel
    CCTK_LOOP3(restrict_3d_cc_o3_rf2, i, j, k, 0, 0, 0, regiext, regjext,
               regkext, dstipadext, dstjpadext, dstkpadext) {

#ifdef CARPET_DEBUG
      if (not(2 * k + 2 + srckoff < srckext and
              2 * j + 2 + srcjoff < srcjext and
              2 * i + 2 + srcioff < srciext and 2 * k - 1 + srckoff >= 0 and
              2 * j - 1 + srcjoff >= 0 and 2 * i - 1 + srcioff >= 0)) {
        cout << "restrict_3d_cc_o3_rf2.cc\n";
        cout << "regext " << regext << "\n";
        cout << "srcext " << srcext << "\n";
        cout << "srcbbox=" << srcbbox << "\n";
        cout << "dstbbox=" << dstbbox << "\n";
        cout << "regbbox=" << regbbox << "\n";
        cout << "i,j,k=" << i << " " << j << " " << k << "\n";
        assert(2 * k + 2 + srckoff < srckext);
        assert(2 * j + 2 + srcjoff < srcjext);
        assert(2 * i + 2 + srcioff < srciext);
        assert(2 * k - 1 + srckoff >= 0);
        assert(2 * j - 1 + srcjoff >= 0);
        assert(2 * i - 1 + srcioff >= 0);
      }
#endif
      dst[DSTIND3(i, j, k)] =
          -f1 * src[SRCIND3(2 * i - 1, 2 * j - 1, 2 * k - 1)] +
          f2 * src[SRCIND3(2 * i, 2 * j - 1, 2 * k - 1)] +
          f2 * src[SRCIND3(2 * i + 1, 2 * j - 1, 2 * k - 1)] -
          f1 * src[SRCIND3(2 * i + 2, 2 * j - 1, 2 * k - 1)] +
          f2 * src[SRCIND3(2 * i - 1, 2 * j, 2 * k - 1)] -
          f3 * src[SRCIND3(2 * i, 2 * j, 2 * k - 1)] -
          f3 * src[SRCIND3(2 * i + 1, 2 * j, 2 * k - 1)] +
          f2 * src[SRCIND3(2 * i + 2, 2 * j, 2 * k - 1)] +
          f2 * src[SRCIND3(2 * i - 1, 2 * j + 1, 2 * k - 1)] -
          f3 * src[SRCIND3(2 * i, 2 * j + 1, 2 * k - 1)] -
          f3 * src[SRCIND3(2 * i + 1, 2 * j + 1, 2 * k - 1)] +
          f2 * src[SRCIND3(2 * i + 2, 2 * j + 1, 2 * k - 1)] -
          f1 * src[SRCIND3(2 * i - 1, 2 * j + 2, 2 * k - 1)] +
          f2 * src[SRCIND3(2 * i, 2 * j + 2, 2 * k - 1)] +
          f2 * src[SRCIND3(2 * i + 1, 2 * j + 2, 2 * k - 1)] -
          f1 * src[SRCIND3(2 * i + 2, 2 * j + 2, 2 * k - 1)] +
          f2 * src[SRCIND3(2 * i - 1, 2 * j - 1, 2 * k)] -
          f3 * src[SRCIND3(2 * i, 2 * j - 1, 2 * k)] -
          f3 * src[SRCIND3(2 * i + 1, 2 * j - 1, 2 * k)] +
          f2 * src[SRCIND3(2 * i + 2, 2 * j - 1, 2 * k)] -
          f3 * src[SRCIND3(2 * i - 1, 2 * j, 2 * k)] +
          f4 * src[SRCIND3(2 * i, 2 * j, 2 * k)] +
          f4 * src[SRCIND3(2 * i + 1, 2 * j, 2 * k)] -
          f3 * src[SRCIND3(2 * i + 2, 2 * j, 2 * k)] -
          f3 * src[SRCIND3(2 * i - 1, 2 * j + 1, 2 * k)] +
          f4 * src[SRCIND3(2 * i, 2 * j + 1, 2 * k)] +
          f4 * src[SRCIND3(2 * i + 1, 2 * j + 1, 2 * k)] -
          f3 * src[SRCIND3(2 * i + 2, 2 * j + 1, 2 * k)] +
          f2 * src[SRCIND3(2 * i - 1, 2 * j + 2, 2 * k)] -
          f3 * src[SRCIND3(2 * i, 2 * j + 2, 2 * k)] -
          f3 * src[SRCIND3(2 * i + 1, 2 * j + 2, 2 * k)] +
          f2 * src[SRCIND3(2 * i + 2, 2 * j + 2, 2 * k)] +
          f2 * src[SRCIND3(2 * i - 1, 2 * j - 1, 2 * k + 1)] -
          f3 * src[SRCIND3(2 * i, 2 * j - 1, 2 * k + 1)] -
          f3 * src[SRCIND3(2 * i + 1, 2 * j - 1, 2 * k + 1)] +
          f2 * src[SRCIND3(2 * i + 2, 2 * j - 1, 2 * k + 1)] -
          f3 * src[SRCIND3(2 * i - 1, 2 * j, 2 * k + 1)] +
          f4 * src[SRCIND3(2 * i, 2 * j, 2 * k + 1)] +
          f4 * src[SRCIND3(2 * i + 1, 2 * j, 2 * k + 1)] -
          f3 * src[SRCIND3(2 * i + 2, 2 * j, 2 * k + 1)] -
          f3 * src[SRCIND3(2 * i - 1, 2 * j + 1, 2 * k + 1)] +
          f4 * src[SRCIND3(2 * i, 2 * j + 1, 2 * k + 1)] +
          f4 * src[SRCIND3(2 * i + 1, 2 * j + 1, 2 * k + 1)] -
          f3 * src[SRCIND3(2 * i + 2, 2 * j + 1, 2 * k + 1)] +
          f2 * src[SRCIND3(2 * i - 1, 2 * j + 2, 2 * k + 1)] -
          f3 * src[SRCIND3(2 * i, 2 * j + 2, 2 * k + 1)] -
          f3 * src[SRCIND3(2 * i + 1, 2 * j + 2, 2 * k + 1)] +
          f2 * src[SRCIND3(2 * i + 2, 2 * j + 2, 2 * k + 1)] -
          f1 * src[SRCIND3(2 * i - 1, 2 * j - 1, 2 * k + 2)] +
          f2 * src[SRCIND3(2 * i, 2 * j - 1, 2 * k + 2)] +
          f2 * src[SRCIND3(2 * i + 1, 2 * j - 1, 2 * k + 2)] -
          f1 * src[SRCIND3(2 * i + 2, 2 * j - 1, 2 * k + 2)] +
          f2 * src[SRCIND3(2 * i - 1, 2 * j, 2 * k + 2)] -
          f3 * src[SRCIND3(2 * i, 2 * j, 2 * k + 2)] -
          f3 * src[SRCIND3(2 * i + 1, 2 * j, 2 * k + 2)] +
          f2 * src[SRCIND3(2 * i + 2, 2 * j, 2 * k + 2)] +
          f2 * src[SRCIND3(2 * i - 1, 2 * j + 1, 2 * k + 2)] -
          f3 * src[SRCIND3(2 * i, 2 * j + 1, 2 * k + 2)] -
          f3 * src[SRCIND3(2 * i + 1, 2 * j + 1, 2 * k + 2)] +
          f2 * src[SRCIND3(2 * i + 2, 2 * j + 1, 2 * k + 2)] -
          f1 * src[SRCIND3(2 * i - 1, 2 * j + 2, 2 * k + 2)] +
          f2 * src[SRCIND3(2 * i, 2 * j + 2, 2 * k + 2)] +
          f2 * src[SRCIND3(2 * i + 1, 2 * j + 2, 2 * k + 2)] -
          f1 * src[SRCIND3(2 * i + 2, 2 * j + 2, 2 * k + 2)];
    }
    CCTK_ENDLOOP3(restrict_3d_cc_o3_rf2);
  }
}

#define TYPECASE(N, T)                                                         \
  template void restrict_3d_cc_o3_rf2(                                         \
      T const *restrict const src, ivect3 const &restrict srcpadext,           \
      ivect3 const &restrict srcext, T *restrict const dst,                    \
      ivect3 const &restrict dstpadext, ivect3 const &restrict dstext,         \
      ibbox3 const &restrict srcbbox, ibbox3 const &restrict dstbbox,          \
      ibbox3 const &restrict, ibbox3 const &restrict regbbox,                  \
      void *extraargs);
#include "typecase.hh"
#undef TYPECASE

} // namespace CarpetLib
