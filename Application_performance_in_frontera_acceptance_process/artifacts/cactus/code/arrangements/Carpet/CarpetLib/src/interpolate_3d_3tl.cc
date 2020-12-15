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

template <typename T>
void interpolate_3d_3tl(T const *restrict const src1, CCTK_REAL const t1,
                        T const *restrict const src2, CCTK_REAL const t2,
                        T const *restrict const src3, CCTK_REAL const t3,
                        ivect3 const &restrict srcpadext,
                        ivect3 const &restrict srcext, T *restrict const dst,
                        CCTK_REAL const t, ivect3 const &restrict dstpadext,
                        ivect3 const &restrict dstext,
                        ibbox3 const &restrict srcbbox,
                        ibbox3 const &restrict dstbbox, ibbox3 const &restrict,
                        ibbox3 const &restrict regbbox, void *extraargs) {
  assert(not extraargs);

  typedef typename typeprops<T>::real RT;

  if (any(srcbbox.stride() != regbbox.stride() or
          dstbbox.stride() != regbbox.stride())) {
    CCTK_ERROR("Internal error: strides disagree");
  }

  if (any(srcbbox.stride() != dstbbox.stride())) {
    CCTK_ERROR("Internal error: strides disagree");
  }

  // This could be handled, but is likely to point to an error
  // elsewhere
  if (regbbox.empty()) {
    CCTK_ERROR("Internal error: region extent is empty");
  }

  if (not regbbox.is_contained_in(srcbbox) or
      not regbbox.is_contained_in(dstbbox)) {
    CCTK_ERROR(
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

  RT const eps = 1.0e-10;

  if (fabs(t1 - t2) < eps or fabs(t1 - t3) < eps or fabs(t2 - t3) < eps) {
    CCTK_ERROR("Internal error: arrays have same time");
  }
  if (t < fmin(fmin(t1, t2), t3) - eps or t > fmax(fmax(t1, t2), t3) + eps) {
    CCTK_ERROR("Internal error: extrapolation in time");
  }

  RT const s1fac = (t - t2) * (t - t3) / ((t1 - t2) * (t1 - t3));
  RT const s2fac = (t - t1) * (t - t3) / ((t2 - t1) * (t2 - t3));
  RT const s3fac = (t - t1) * (t - t2) / ((t3 - t1) * (t3 - t2));

// Loop over region
#pragma omp parallel
  CCTK_LOOP3(interpolate_3d_3tl, i, j, k, 0, 0, 0, regiext, regjext, regkext,
             dstipadext, dstjpadext, dstkpadext) {

    dst[DSTIND3(i, j, k)] = +s1fac * src1[SRCIND3(i, j, k)] +
                            s2fac * src2[SRCIND3(i, j, k)] +
                            s3fac * src3[SRCIND3(i, j, k)];
  }
  CCTK_ENDLOOP3(interpolate_3d_3tl);
}

#define TYPECASE(N, T)                                                         \
  template void interpolate_3d_3tl(                                            \
      T const *restrict const src1, CCTK_REAL const t1,                        \
      T const *restrict const src2, CCTK_REAL const t2,                        \
      T const *restrict const src3, CCTK_REAL const t3,                        \
      ivect3 const &restrict srcpadext, ivect3 const &restrict srcext,         \
      T *restrict const dst, CCTK_REAL const t,                                \
      ivect3 const &restrict dstpadext, ivect3 const &restrict dstext,         \
      ibbox3 const &restrict srcbbox, ibbox3 const &restrict dstbbox,          \
      ibbox3 const &restrict, ibbox3 const &restrict regbbox,                  \
      void *extraargs);
#include "typecase.hh"
#undef TYPECASE

} // namespace CarpetLib
