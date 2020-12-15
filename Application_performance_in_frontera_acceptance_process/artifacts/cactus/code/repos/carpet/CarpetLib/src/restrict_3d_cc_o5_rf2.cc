#include <cctk.h>
#include <cctk_Parameters.h>

#include <loopcontrol.h>

#include <algorithm>
#include <cassert>
#include <cmath>

#include "gdata.hh"
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
#define SRCOFF3(i, j, k)                                                       \
  offset3(srcioff + (i), srcjoff + (j), srckoff + (k), srciext, srcjext,       \
          srckext)
#define DSTOFF3(i, j, k)                                                       \
  offset3(dstioff + (i), dstjoff + (j), dstkoff + (k), dstiext, dstjext,       \
          dstkext)

// This operator offers fifth-order accurate restriction operators for cell
// centered grids when use_higher_order_restriction is set. This interpolation
// is done for samples at the cell centres (so this is not a ppm scheme or
// anything like that).  It really only fits a polynomial of the form
//    f(x,y,z) = \sum_{i,j,k=0}^5 a_{i,j,k} x^i y^j z^k
// to the fine cells and evaluates at x=y=z=0. So it is good for the metric,
// but bad for matter (since it will destroy the conservation).  Because of
// this it should not be used for grid functions whose transport operator is
// not WENO or ENO which hopefully excludes all matter variables.
// TODO: use prolongate_3d_rf2 instead of writing new operators, its the same
// interpolation.

namespace {

// 1D restriction
template <typename T>
inline T restrict1(T const *restrict const p, size_t const d1) {
  typedef typename typeprops<T>::real RT;
  RT const den = 256;
  RT const f2 = 3 / den, f1 = 25 / den, f0 = 150 / den;
  T const res = +f2 * p[-2] - f1 * p[-1] + f0 * p[-0] + f0 * p[+1] -
                f1 * p[+2] + f2 * p[+3];
  return res;
}

// 2D restriction
template <typename T>
inline T restrict2(T const *restrict const p, size_t const d1,
                   size_t const d2) {
  typedef typename typeprops<T>::real RT;
  RT const den = 256;
  RT const f2 = 3 / den, f1 = 25 / den, f0 = 150 / den;
  T const res =
      +f2 * restrict1(p - 2 * d2, d1) - f1 * restrict1(p - 1 * d2, d1) +
      f0 * restrict1(p - 0 * d2, d1) + f0 * restrict1(p + 1 * d2, d1) -
      f1 * restrict1(p + 2 * d2, d1) + f2 * restrict1(p + 3 * d2, d1);
  return res;
}

// 3D restriction
template <typename T>
inline T restrict3(T const *restrict const p, size_t const d1, size_t const d2,
                   size_t const d3) {
  typedef typename typeprops<T>::real RT;
  RT const den = 256;
  RT const f2 = 3 / den, f1 = 25 / den, f0 = 150 / den;
  T const res =
      +f2 * restrict2(p - 2 * d3, d1, d2) - f1 * restrict2(p - 1 * d3, d1, d2) +
      f0 * restrict2(p - 0 * d3, d1, d2) + f0 * restrict2(p + 1 * d3, d1, d2) -
      f1 * restrict2(p + 2 * d3, d1, d2) + f2 * restrict2(p + 3 * d3, d1, d2);
  return res;
}

} // namespace

template <typename T>
void restrict_3d_cc_o5_rf2(
    T const *restrict const src, ivect3 const &restrict srcpadext,
    ivect3 const &restrict srcext, T *restrict const dst,
    ivect3 const &restrict dstpadext, ivect3 const &restrict dstext,
    ibbox3 const &restrict srcbbox, ibbox3 const &restrict dstbbox,
    ibbox3 const &restrict, ibbox3 const &restrict regbbox, void *extraargs) {
  assert(not extraargs);

  DECLARE_CCTK_PARAMETERS;

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

  // size_t const srcdi == SRCOFF3(1,0,0) - SRCOFF3(0,0,0);
  size_t const srcdi = 1;
  assert(srcdi == SRCOFF3(1, 0, 0) - SRCOFF3(0, 0, 0));
  size_t const srcdj = SRCOFF3(0, 1, 0) - SRCOFF3(0, 0, 0);
  size_t const srcdk = SRCOFF3(0, 0, 1) - SRCOFF3(0, 0, 0);

  if (not use_loopcontrol_in_operators) {

// Loop over coarse region
#pragma omp parallel for collapse(3)
    for (int k = 0; k < regkext; ++k) {
      for (int j = 0; j < regjext; ++j) {
        for (int i = 0; i < regiext; ++i) {

#ifdef CARPET_DEBUG
          if (not(2 * k + 3 + srckoff < srckext and
                  2 * j + 3 + srcjoff < srcjext and
                  2 * i + 3 + srcioff < srciext and 2 * k - 2 + srckoff >= 0 and
                  2 * j - 2 + srcjoff >= 0 and 2 * i - 2 + srcioff >= 0)) {
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
          dst[DSTIND3(i, j, k)] = CarpetLib::restrict3(
              &src[SRCIND3(2 * i, 2 * j, 2 * k)], srcdi, srcdj, srcdk);
        }
      }
    }

  } else {

// Loop over coarse region
#pragma omp parallel
    CCTK_LOOP3(restrict_3d_cc_o5_rf2, i, j, k, 0, 0, 0, regiext, regjext,
               regkext, dstipadext, dstjpadext, dstkpadext) {

#ifdef CARPET_DEBUG
      if (not(2 * k + 3 + srckoff < srckext and
              2 * j + 3 + srcjoff < srcjext and
              2 * i + 3 + srcioff < srciext and 2 * k - 2 + srckoff >= 0 and
              2 * j - 2 + srcjoff >= 0 and 2 * i - 2 + srcioff >= 0)) {
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
      dst[DSTIND3(i, j, k)] = CarpetLib::restrict3(
          &src[SRCIND3(2 * i, 2 * j, 2 * k)], srcdi, srcdj, srcdk);
    }
    CCTK_ENDLOOP3(restrict_3d_cc_o5_rf2);
  }
}

#define TYPECASE(N, T)                                                         \
  template void restrict_3d_cc_o5_rf2(                                         \
      T const *restrict const src, ivect3 const &restrict srcpadext,           \
      ivect3 const &restrict srcext, T *restrict const dst,                    \
      ivect3 const &restrict dstpadext, ivect3 const &restrict dstext,         \
      ibbox3 const &restrict srcbbox, ibbox3 const &restrict dstbbox,          \
      ibbox3 const &restrict, ibbox3 const &restrict regbbox,                  \
      void *extraargs);
#include "typecase.hh"
#undef TYPECASE

} // namespace CarpetLib
