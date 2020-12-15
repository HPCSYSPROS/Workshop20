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
void prolongate_3d_dgfe_rf2(
    T const *restrict const src, ivect3 const &restrict srcpadext,
    ivect3 const &restrict srcext, T *restrict const dst,
    ivect3 const &restrict dstpadext, ivect3 const &restrict dstext,
    ibbox3 const &restrict srcbbox, ibbox3 const &restrict dstbbox,
    ibbox3 const &restrict, ibbox3 const &restrict regbbox,
    void *const extraargs) {
#ifdef HRSCC_GLL_ELEMENT_HH
  assert(not extraargs);

  static_assert(ORDER >= 0, "ORDER must be non-negative");

  if (any(srcbbox.stride() <= regbbox.stride() or
          dstbbox.stride() != regbbox.stride())) {
    CCTK_WARN(0, "Internal error: strides disagree");
  }

  if (any(srcbbox.stride() != reffact2 * dstbbox.stride())) {
    CCTK_WARN(
        0,
        "Internal error: source strides are not twice the destination strides");
  }

  if (any(dstbbox.stride() % 2 != 0)) {
    CCTK_WARN(0, "Internal error: destination strides are not even");
  }

  // This could be handled, but is likely to point to an error
  // elsewhere
  if (regbbox.empty()) {
    CCTK_WARN(0, "Internal error: region extent is empty");
  }

  bvect3 const is_upper_face =
      (regbbox.lower() - srcbbox.lower() + regbbox.stride() / 2) %
          srcbbox.stride() !=
      0;

  ivect3 const regext = regbbox.shape() / regbbox.stride();
  assert(all((regbbox.lower() - srcbbox.lower() +
              either(is_upper_face, -1, +1) * regbbox.stride() / 2) %
                 srcbbox.stride() ==
             0));
  ivect3 const srcoff = (regbbox.lower() - srcbbox.lower() +
                         either(is_upper_face, -1, +1) * regbbox.stride() / 2) /
                        srcbbox.stride();
  assert(all((regbbox.lower() - dstbbox.lower()) % regbbox.stride() == 0));
  ivect3 const dstoff = (regbbox.lower() - dstbbox.lower()) / regbbox.stride();

  if (not regbbox.is_contained_in(srcbbox) or
      not regbbox.is_contained_in(dstbbox)) {
    cerr << "ORDER=" << ORDER << "\n"
         << "regbbox=" << regbbox << "\n"
         << "dstbbox=" << dstbbox << "\n"
         << "regbbox=" << regbbox << "\n"
         << "srcbbox=" << srcbbox << "\n";
    CCTK_WARN(0,
              "Internal error: region extent is not contained in array extent");
  }

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

  if (regext[0] == 1) {
    // 2D prolongation on x face

    assert(not is_upper_face[1]);
    assert(not is_upper_face[2]);

    // Ensure we traverse an even integer number of elements
    assert(regext[1] % (2 * (ORDER + 1)) == 0);
    assert(regext[2] % (2 * (ORDER + 1)) == 0);

    int const srcdi = 1; // 2d face
    int const srcdj = SRCOFF3(0, 1, 0) - SRCOFF3(0, 0, 0);
    int const srcdk = SRCOFF3(0, 0, 1) - SRCOFF3(0, 0, 0);
    int const dstdi = 1; // 2d face
    int const dstdj = DSTOFF3(0, 1, 0) - DSTOFF3(0, 0, 0);
    int const dstdk = DSTOFF3(0, 0, 1) - DSTOFF3(0, 0, 0);
    int const srcstr2d[2] = {srcdj, srcdk};
    int const dststr2d[2] = {dstdj, dstdk};

    // Loop over fine region
    ptrdiff_t const i = 0;
#pragma omp parallel for collapse(2)
    // Zwicky's Intel compiler 11.1 ices on ptrdiff_t
    for (/*ptrdiff_t*/ int k = 0; k < regkext; k += 2 * (ORDER + 1)) {
      for (/*ptrdiff_t*/ int j = 0; j < regjext; j += 2 * (ORDER + 1)) {
        GLLElement<ORDER>::prolongate_2D(
            &src[SRCIND3(srcioff + i, srcjoff + j, srckoff + k)], srcstr2d,
            &dst[DSTIND3(dstioff + 2 * i, dstjoff + 2 * j, dstkoff + 2 * k)],
            dststr2d);
      }
    }

  } else if (regext[1] == 1) {
    // 2D prolongation on y face

    assert(not is_upper_face[0]);
    assert(not is_upper_face[2]);

    // Ensure we traverse an even integer number of elements
    assert(regext[0] % (2 * (ORDER + 1)) == 0);
    assert(regext[2] % (2 * (ORDER + 1)) == 0);

    // int const srcdi = SRCOFF3(1,0,0) - SRCOFF3(0,0,0);
    int const srcdi = 1;
    assert(srcdi == SRCOFF3(1, 0, 0) - SRCOFF3(0, 0, 0));
    int const srcdj = 1; // 2d face
    int const srcdk = SRCOFF3(0, 0, 1) - SRCOFF3(0, 0, 0);
    // int const dstdi = DSTOFF3(1,0,0) - DSTOFF3(0,0,0);
    int const dstdi = 1;
    assert(dstdi == DSTOFF3(1, 0, 0) - DSTOFF3(0, 0, 0));
    int const dstdj = 1; // 2d face
    int const dstdk = DSTOFF3(0, 0, 1) - DSTOFF3(0, 0, 0);
    int const srcstr2d[2] = {srcdi, srcdk};
    int const dststr2d[2] = {dstdi, dstdk};

    // Loop over fine region
    ptrdiff_t const j = 0;
#pragma omp parallel for collapse(2)
    // Zwicky's Intel compiler 11.1 ices on ptrdiff_t
    for (/*ptrdiff_t*/ int k = 0; k < regkext; k += 2 * (ORDER + 1)) {
      for (/*ptrdiff_t*/ int i = 0; i < regiext; i += 2 * (ORDER + 1)) {
        GLLElement<ORDER>::prolongate_2D(
            &src[SRCIND3(srcioff + i, srcjoff + j, srckoff + k)], srcstr2d,
            &dst[DSTIND3(dstioff + 2 * i, dstjoff + 2 * j, dstkoff + 2 * k)],
            dststr2d);
      }
    }

  } else if (regext[2] == 1) {
    // 2D prolongation on z face

    assert(not is_upper_face[0]);
    assert(not is_upper_face[1]);

    // Ensure we traverse an even integer number of elements
    assert(regext[0] % (2 * (ORDER + 1)) == 0);
    assert(regext[1] % (2 * (ORDER + 1)) == 0);

    // int const srcdi = SRCOFF3(1,0,0) - SRCOFF3(0,0,0);
    int const srcdi = 1;
    assert(srcdi == SRCOFF3(1, 0, 0) - SRCOFF3(0, 0, 0));
    int const srcdj = SRCOFF3(0, 1, 0) - SRCOFF3(0, 0, 0);
    int const srcdk = 1; // 2d face
    // int const dstdi = DSTOFF3(1,0,0) - DSTOFF3(0,0,0);
    int const dstdi = 1;
    assert(dstdi == DSTOFF3(1, 0, 0) - DSTOFF3(0, 0, 0));
    int const dstdj = DSTOFF3(0, 1, 0) - DSTOFF3(0, 0, 0);
    int const dstdk = 1; // 2d face
    int const srcstr2d[2] = {srcdi, srcdj};
    int const dststr2d[2] = {dstdi, dstdj};

    // Loop over fine region
    ptrdiff_t const k = 0;
#pragma omp parallel for collapse(2)
    // Zwicky's Intel compiler 11.1 ices on ptrdiff_t
    for (/*ptrdiff_t*/ int j = 0; j < regjext; j += 2 * (ORDER + 1)) {
      for (/*ptrdiff_t*/ int i = 0; i < regiext; i += 2 * (ORDER + 1)) {
        GLLElement<ORDER>::prolongate_2D(
            &src[SRCIND3(srcioff + i, srcjoff + j, srckoff + k)], srcstr2d,
            &dst[DSTIND3(dstioff + 2 * i, dstjoff + 2 * j, dstkoff + 2 * k)],
            dststr2d);
      }
    }

  } else {
    // 3D prolongation

    assert(not any(is_upper_face));

    // Ensure we traverse an even integer number of elements
    assert(all(regext % (2 * (ORDER + 1)) == 0));

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

// Loop over fine region
#pragma omp parallel for collapse(3)
    // Zwicky's Intel compiler 11.1 ices on ptrdiff_t
    for (/*ptrdiff_t*/ int k = 0; k < regkext; k += 2 * (ORDER + 1)) {
      for (/*ptrdiff_t*/ int j = 0; j < regjext; j += 2 * (ORDER + 1)) {
        for (/*ptrdiff_t*/ int i = 0; i < regiext; i += 2 * (ORDER + 1)) {
          GLLElement<ORDER>::prolongate_full(
              &src[SRCIND3(srcioff + i, srcjoff + j, srckoff + k)], srcstr,
              &dst[DSTIND3(dstioff + 2 * i, dstjoff + 2 * j, dstkoff + 2 * k)],
              dststr);
        }
      }
    }
  }

#else
  // HRSCCore is not available
  assert(0);
#endif
}

#define TYPECASE(N, T)                                                         \
                                                                               \
  template void prolongate_3d_dgfe_rf2<T, 5>(                                  \
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
void prolongate_3d_dgfe_rf2<CCTK_COMPLEX, 5>(
    CCTK_COMPLEX const *restrict const src, ivect3 const &restrict srcpadext,
    ivect3 const &restrict srcext, CCTK_COMPLEX *restrict const dst,
    ivect3 const &restrict dstpadext, ivect3 const &restrict dstext,
    ibbox3 const &restrict srcbbox, ibbox3 const &restrict dstbbox,
    ibbox3 const &restrict, ibbox3 const &restrict regbbox, void *extraargs) {
  assert(0);
}

} // namespace CarpetLib
