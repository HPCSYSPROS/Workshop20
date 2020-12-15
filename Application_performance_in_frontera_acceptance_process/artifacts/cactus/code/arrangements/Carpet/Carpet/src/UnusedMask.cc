#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <carpet.hh>

#include <loopcontrol.h>

namespace Carpet {

using namespace std;

void CarpetUnusedMask(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  dh const &dd = *vdd.AT(map);
  ibbox const &ext =
      dd.light_boxes.AT(mglevel).AT(reflevel).AT(component).exterior;
  ibset const &unused_region =
      dd.local_boxes.AT(mglevel).AT(reflevel).AT(local_component).unused_region;

  assert(dim == 3);
  // Zero out
  LOOP_OVER_BSET(cctkGH, ext, box, imin, imax) {
#pragma omp parallel
    CCTK_LOOP3(unused_mask_zero, i, j, k, imin[0], imin[1], imin[2], imax[0],
               imax[1], imax[2], cctk_ash[0], cctk_ash[1], cctk_ash[2]) {
      CCTK_INT i3D = CCTK_GFINDEX3D(cctkGH, i, j, k);
      carpet_unusedpoints_mask[i3D] = 0;
    }
    CCTK_ENDLOOP3(unused_mask_zero);
  }
  END_LOOP_OVER_BSET;

  // Set it where unused
  LOOP_OVER_BSET(cctkGH, unused_region, box, imin, imax) {
#pragma omp parallel
    CCTK_LOOP3(unused_mask_set, i, j, k, imin[0], imin[1], imin[2], imax[0],
               imax[1], imax[2], cctk_ash[0], cctk_ash[1], cctk_ash[2]) {
      CCTK_INT i3D = CCTK_GFINDEX3D(cctkGH, i, j, k);
      carpet_unusedpoints_mask[i3D] = 1;
    }
    CCTK_ENDLOOP3(unused_mask_set);
  }
  END_LOOP_OVER_BSET;

} // CarpetUnusedMask
}
