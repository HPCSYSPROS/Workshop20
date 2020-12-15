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
#define SRCOFF3(i, j, k)                                                       \
  offset3(srcioff + (i), srcjoff + (j), srckoff + (k), srciext, srcjext,       \
          srckext)
#define DSTOFF3(i, j, k)                                                       \
  offset3(dstioff + (i), dstjoff + (j), dstkoff + (k), dstiext, dstjext,       \
          dstkext)

// 0D "restriction"
template <typename T> struct restrict0 {
  static inline T call(T const *restrict const p) { return *p; }
};

// 1D restriction
template <typename T, int centi> struct restrict1 {
  static inline T call(T const *restrict const p, size_t const d1);
};
template <typename T> struct restrict1<T, 0> {
  static inline T call(T const *restrict const p, size_t const d1) {
    T const res = restrict0<T>::call(p);
    return res;
  }
};
template <typename T> struct restrict1<T, 1> {
  static inline T call(T const *restrict const p, size_t const d1) {
    typedef typename typeprops<T>::real RT;
    RT const one = 1;
    RT const half = one / 2;
    T const res =
        half * restrict0<T>::call(p) + half * restrict0<T>::call(p + d1);
    return res;
  }
};

// 2D restriction
template <typename T, int centi, int centj> struct restrict2 {
  static inline T call(T const *restrict const p, size_t const d1,
                       size_t const d2);
};
template <typename T, int centi> struct restrict2<T, centi, 0> {
  static inline T call(T const *restrict const p, size_t const d1,
                       size_t const d2) {
    T const res = restrict1<T, centi>::call(p, d1);
    return res;
  }
};
template <typename T, int centi> struct restrict2<T, centi, 1> {
  static inline T call(T const *restrict const p, size_t const d1,
                       size_t const d2) {
    typedef typename typeprops<T>::real RT;
    RT const one = 1;
    RT const half = one / 2;
    T const res = half * restrict1<T, centi>::call(p, d1) +
                  half * restrict1<T, centi>::call(p + d2, d1);
    return res;
  }
};

// 3D restriction
template <typename T, int centi, int centj, int centk> struct restrict3 {
  static inline T call(T const *restrict const p, size_t const d1,
                       size_t const d2, size_t const d3);
};
template <typename T, int centi, int centj>
struct restrict3<T, centi, centj, 0> {
  static inline T call(T const *restrict const p, size_t const d1,
                       size_t const d2, size_t const d3) {
    T const res = restrict2<T, centi, centj>::call(p, d1, d2);
    return res;
  }
};
template <typename T, int centi, int centj>
struct restrict3<T, centi, centj, 1> {
  static inline T call(T const *restrict const p, size_t const d1,
                       size_t const d2, size_t const d3) {
    typedef typename typeprops<T>::real RT;
    RT const one = 1;
    RT const half = one / 2;
    T const res = half * restrict2<T, centi, centj>::call(p, d1, d2) +
                  half * restrict2<T, centi, centj>::call(p + d3, d1, d2);
    return res;
  }
};

template <typename T, int centi, int centj, int centk>
void restrict_3d_vc_rf2(T const *restrict const src,
                        ivect3 const &restrict srcpadext,
                        ivect3 const &restrict srcext, T *restrict const dst,
                        ivect3 const &restrict dstpadext,
                        ivect3 const &restrict dstext,
                        ibbox3 const &restrict srcbbox,
                        ibbox3 const &restrict dstbbox,
                        ibbox3 const &restrict srcregbbox,
                        ibbox3 const &restrict regbbox, void *extraargs) {
  assert(not extraargs);

  DECLARE_CCTK_PARAMETERS;

  // false: vertex centered, true: cell centered
  ivect const icent(centi, centj, centk);
  assert(all(icent == 0 or icent == 1));

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

  // shift vertex centered directions to lower interface (see
  // Refluxing for conventions)
  ivect const ivert(icent == 0);
  ibbox3 const unshifted_srcbbox = srcbbox.shift(-ivert, 2);
  ibbox3 const unshifted_dstbbox = dstbbox.shift(-ivert, 2);
  // ibbox3 const unshifted_srcregbbox = srcregbbox.shift(-ivert,2);
  ibbox3 const unshifted_regbbox = regbbox.shift(-ivert, 2);

  if (not unshifted_regbbox.expanded_for(unshifted_srcbbox)
              .is_contained_in(unshifted_srcbbox) or
      not unshifted_regbbox.is_contained_in(unshifted_dstbbox)) {
    cerr << "unshifted_srcbbox: " << unshifted_srcbbox << endl
         << "unshifted_dstbbox: " << unshifted_dstbbox << endl
         << "unshifted_regbbox: " << unshifted_regbbox << endl;
    CCTK_WARN(0,
              "Internal error: region extent is not contained in array extent");
  }

  ivect3 const regext = regbbox.shape() / regbbox.stride();
  assert(all(srcbbox.stride() % 2 == 0));
  if (not(all((regbbox.lower() - srcbbox.lower() - srcbbox.stride() / 2) %
                  srcbbox.stride() ==
              0))) {
    cout << "restrict_3d_vc_rf2.cc\n";
    cout << "regbbox=" << regbbox << "\n";
    cout << "srcbbox=" << srcbbox << "\n";
    cout << "icent=" << icent << "\n";
  }
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
          if (not(2 * k + centk < srckext and 2 * j + centj < srcjext and
                  2 * i + centi < srciext)) {
            cout << "restrict_3d_vc_rf2.cc\n";
            cout << "regext " << regext << "\n";
            cout << "srcext " << srcext << "\n";
            cout << "srcbbox=" << srcbbox << "\n";
            cout << "dstbbox=" << dstbbox << "\n";
            cout << "regbbox=" << regbbox << "\n";
            cout << "srcregbbox=" << srcregbbox << "\n";
            cout << "icent=" << icent << "\n";
          }
          assert(2 * k + centk < srckext and 2 * j + centj < srcjext and
                 2 * i + centi < srciext);
#endif

          dst[DSTIND3(i, j, k)] = restrict3<T, centi, centj, centk>::call(
              &src[SRCIND3(2 * i, 2 * j, 2 * k)], srcdi, srcdj, srcdk);
        }
      }
    }

  } else {

// Loop over coarse region
#pragma omp parallel
    CCTK_LOOP3(restrict_3d_vc_rf2, i, j, k, 0, 0, 0, regiext, regjext, regkext,
               dstipadext, dstjpadext, dstkpadext) {
#ifdef CARPET_DEBUG
      if (not(2 * k + centk < srckext and 2 * j + centj < srcjext and
              2 * i + centi < srciext)) {
        cout << "restrict_3d_vc_rf2.cc\n";
        cout << "regext " << regext << "\n";
        cout << "srcext " << srcext << "\n";
        cout << "srcbbox=" << srcbbox << "\n";
        cout << "dstbbox=" << dstbbox << "\n";
        cout << "regbbox=" << regbbox << "\n";
        cout << "srcregbbox=" << srcregbbox << "\n";
        cout << "icent=" << icent << "\n";
      }
      assert(2 * k + centk < srckext and 2 * j + centj < srcjext and
             2 * i + centi < srciext);
#endif

      dst[DSTIND3(i, j, k)] = restrict3<T, centi, centj, centk>::call(
          &src[SRCIND3(2 * i, 2 * j, 2 * k)], srcdi, srcdj, srcdk);
    }
    CCTK_ENDLOOP3(restrict_3d_vc_rf2);
  }
}

#define TYPECASE(N, T)                                                         \
  template void restrict_3d_vc_rf2<T, 0, 0, 0>(                                \
      T const *restrict const src, ivect3 const &restrict srcpadext,           \
      ivect3 const &restrict srcext, T *restrict const dst,                    \
      ivect3 const &restrict dstpadext, ivect3 const &restrict dstext,         \
      ibbox3 const &restrict srcbbox, ibbox3 const &restrict dstbbox,          \
      ibbox3 const &restrict, ibbox3 const &restrict regbbox,                  \
      void *extraargs);                                                        \
  template void restrict_3d_vc_rf2<T, 0, 0, 1>(                                \
      T const *restrict const src, ivect3 const &restrict srcpadext,           \
      ivect3 const &restrict srcext, T *restrict const dst,                    \
      ivect3 const &restrict dstpadext, ivect3 const &restrict dstext,         \
      ibbox3 const &restrict srcbbox, ibbox3 const &restrict dstbbox,          \
      ibbox3 const &restrict, ibbox3 const &restrict regbbox,                  \
      void *extraargs);                                                        \
  template void restrict_3d_vc_rf2<T, 0, 1, 0>(                                \
      T const *restrict const src, ivect3 const &restrict srcpadext,           \
      ivect3 const &restrict srcext, T *restrict const dst,                    \
      ivect3 const &restrict dstpadext, ivect3 const &restrict dstext,         \
      ibbox3 const &restrict srcbbox, ibbox3 const &restrict dstbbox,          \
      ibbox3 const &restrict, ibbox3 const &restrict regbbox,                  \
      void *extraargs);                                                        \
  template void restrict_3d_vc_rf2<T, 0, 1, 1>(                                \
      T const *restrict const src, ivect3 const &restrict srcpadext,           \
      ivect3 const &restrict srcext, T *restrict const dst,                    \
      ivect3 const &restrict dstpadext, ivect3 const &restrict dstext,         \
      ibbox3 const &restrict srcbbox, ibbox3 const &restrict dstbbox,          \
      ibbox3 const &restrict, ibbox3 const &restrict regbbox,                  \
      void *extraargs);                                                        \
  template void restrict_3d_vc_rf2<T, 1, 0, 0>(                                \
      T const *restrict const src, ivect3 const &restrict srcpadext,           \
      ivect3 const &restrict srcext, T *restrict const dst,                    \
      ivect3 const &restrict dstpadext, ivect3 const &restrict dstext,         \
      ibbox3 const &restrict srcbbox, ibbox3 const &restrict dstbbox,          \
      ibbox3 const &restrict, ibbox3 const &restrict regbbox,                  \
      void *extraargs);                                                        \
  template void restrict_3d_vc_rf2<T, 1, 0, 1>(                                \
      T const *restrict const src, ivect3 const &restrict srcpadext,           \
      ivect3 const &restrict srcext, T *restrict const dst,                    \
      ivect3 const &restrict dstpadext, ivect3 const &restrict dstext,         \
      ibbox3 const &restrict srcbbox, ibbox3 const &restrict dstbbox,          \
      ibbox3 const &restrict, ibbox3 const &restrict regbbox,                  \
      void *extraargs);                                                        \
  template void restrict_3d_vc_rf2<T, 1, 1, 0>(                                \
      T const *restrict const src, ivect3 const &restrict srcpadext,           \
      ivect3 const &restrict srcext, T *restrict const dst,                    \
      ivect3 const &restrict dstpadext, ivect3 const &restrict dstext,         \
      ibbox3 const &restrict srcbbox, ibbox3 const &restrict dstbbox,          \
      ibbox3 const &restrict, ibbox3 const &restrict regbbox,                  \
      void *extraargs);                                                        \
  template void restrict_3d_vc_rf2<T, 1, 1, 1>(                                \
      T const *restrict const src, ivect3 const &restrict srcpadext,           \
      ivect3 const &restrict srcext, T *restrict const dst,                    \
      ivect3 const &restrict dstpadext, ivect3 const &restrict dstext,         \
      ibbox3 const &restrict srcbbox, ibbox3 const &restrict dstbbox,          \
      ibbox3 const &restrict, ibbox3 const &restrict regbbox,                  \
      void *extraargs);
#include "typecase.hh"
#undef TYPECASE

} // namespace CarpetLib
