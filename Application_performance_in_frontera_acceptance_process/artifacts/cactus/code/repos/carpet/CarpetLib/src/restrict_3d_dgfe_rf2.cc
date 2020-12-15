#include <cctk.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>

#include <hrscc.hh>

#include "operator_prototypes_3d.hh"
#include "typeprops.hh"

using namespace std;
#ifdef HRSCC_GLL_ELEMENT_HH
using namespace hrscc;
#endif

namespace CarpetLib {

#define SRCIND3(i, j, k)                                                       \
  index3(i, j, k, srcipadext, srcjpadext, srckpadext, srciext, srcjext, srckext)
#define DSTIND3(i, j, k)                                                       \
  index3(i, j, k, dstipadext, dstjpadext, dstkpadext, dstiext, dstjext, dstkext)
#define SRCOFF3(i, j, k) offset3(i, j, k, srciext, srcjext, srckext)
#define DSTOFF3(i, j, k) offset3(i, j, k, dstiext, dstjext, dstkext)

template <typename T, int ORDER>
void restrict_3d_dgfe_rf2(
    T const *restrict const src, ivect3 const &restrict srcpadext,
    ivect3 const &restrict srcext, T *restrict const dst,
    ivect3 const &restrict dstpadext, ivect3 const &restrict dstext,
    ibbox3 const &restrict srcbbox, ibbox3 const &restrict dstbbox,
    ibbox3 const &restrict, ibbox3 const &restrict regbbox,
    void *const extraargs) {
#ifdef HRSCC_GLL_ELEMENT_HH
  assert(not extraargs);

  static_assert(ORDER >= 0, "ORDER must be non-negative");

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

  if (not regbbox.expanded_for(srcbbox).is_contained_in(srcbbox) or
      not regbbox.is_contained_in(dstbbox)) {
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

  // int const srcdi = SRCOFF3(1,0,0) - SRCOFF3(0,0,0);
  int const srcdi = 1;
  assert(srcdi == SRCOFF3(1, 0, 0) - SRCOFF3(0, 0, 0));
  int const srcdj = SRCOFF3(0, 1, 0) - SRCOFF3(0, 0, 0);
  int const srcdk = SRCOFF3(0, 0, 1) - SRCOFF3(0, 0, 0);

  // int const dstdi = DSTOFF3(1,0,0) - DSTOFF3(0,0,0);
  int const dstdi = 1;
  assert(dstdi == DSTOFF3(1, 0, 0) - DSTOFF3(0, 0, 0));
  int const dstdj = DSTOFF3(0, 1, 0) - DSTOFF3(0, 0, 0);
  int const dstdk = DSTOFF3(0, 0, 1) - DSTOFF3(0, 0, 0);

  int const srcstr[3] = {srcdi, srcdj, srcdk};
  int const dststr[3] = {dstdi, dstdj, dstdk};

  // Ensure we traverse an integer number of elements
  assert(all(regext % (ORDER + 1) == 0));

// Loop over coarse region
#pragma omp parallel for collapse(3)
  // Zwicky's Intel compiler 11.1 ices on ptrdiff_t
  for (/*ptrdiff_t*/ int k = 0; k < regkext; k += ORDER + 1) {
    for (/*ptrdiff_t*/ int j = 0; j < regjext; j += ORDER + 1) {
      for (/*ptrdiff_t*/ int i = 0; i < regiext; i += ORDER + 1) {
        GLLElement<ORDER>::restrict_full(
            &src[SRCIND3(srcioff + 2 * i, srcjoff + 2 * j, srckoff + 2 * k)],
            srcstr, &dst[DSTIND3(dstioff + i, dstjoff + j, dstkoff + k)],
            dststr);
      }
    }
  }

#else
  // HRSCCore is not available
  assert(0);
#endif
}

#define TYPECASE(N, T)                                                         \
  template void restrict_3d_dgfe_rf2<T, 5>(                                    \
      T const *restrict const src, ivect3 const &restrict srcpadext,           \
      ivect3 const &restrict srcext, T *restrict const dst,                    \
      ivect3 const &restrict dstpadext, ivect3 const &restrict dstext,         \
      ibbox3 const &restrict srcbbox, ibbox3 const &restrict dstbbox,          \
      ibbox3 const &restrict, ibbox3 const &restrict regbbox,                  \
      void *extraargs);

#define CARPET_NO_INT
#define CARPET_NO_COMPLEX
#include "typecase.hh"
#undef TYPECASE

template <>
void restrict_3d_dgfe_rf2<CCTK_COMPLEX, 5>(
    CCTK_COMPLEX const *restrict const src, ivect3 const &restrict srcpadext,
    ivect3 const &restrict srcext, CCTK_COMPLEX *restrict const dst,
    ivect3 const &restrict dstpadext, ivect3 const &restrict dstext,
    ibbox3 const &restrict srcbbox, ibbox3 const &restrict dstbbox,
    ibbox3 const &restrict, ibbox3 const &restrict regbbox, void *extraargs) {
  assert(0);
}

} // namespace CarpetLib
