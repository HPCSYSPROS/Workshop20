#include <cassert>
#include <istream>
#include <sstream>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <dh.hh>

#include <carpet.hh>

#include <loopcontrol.h>

#include "bits.h"

namespace CarpetMask {

using namespace std;
using namespace Carpet;

/**
 * Reduce the weight on the current and the next coarser level to
 * make things consistent. Set the weight to 0 inside the active
 * region of the next coarser level, maybe to 1/2 near the boundary
 * of that region, and also to 1/2 near the prolongation boundary of
 * this level.
 */

extern "C" {
void CarpetMaskSetup(CCTK_ARGUMENTS);
}

void CarpetMaskSetup(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;

  if (not is_singlemap_mode()) {
    CCTK_ERROR("This routine may only be called in singlemap mode");
  }

  dh const &dd = *vdd.AT(Carpet::map);

  if (reflevel > 0) {
    ivect const reffact =
        spacereffacts.AT(reflevel) / spacereffacts.AT(reflevel - 1);
    assert(all(reffact == 2));
  }

  // Set prolongation boundaries and restricted region of this level
  BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, CCTK_GF) {
    DECLARE_CCTK_ARGUMENTS;

    // Set the weight in the interior of the not_active and the
    // fine_active regions to zero, and set the weight on the
    // boundary of the not_active and fine_active regions to 1/2.

    // All regions and boundaries are masking out convex shapes, i.e.
    // the remaining grid points from a concave shape. That is, each
    // masking step removes 1/8, 1/4, or 1/2 of a cell (in three
    // dimensions).

    dh::local_dboxes const &local_box =
        dd.local_boxes.AT(mglevel).AT(reflevel).AT(local_component);

    for (int neighbour = 0; neighbour < ipow(3, dim); ++neighbour) {
      ivect shift;
      {
        int itmp = neighbour;
        for (int d = 0; d < dim; ++d) {
          shift[d] = itmp % 3 - 1; // [-1 ... +1]
          itmp /= 3;
        }
        assert(itmp == 0);
      }

      ibset const &boxes = local_box.prolongation_boundary.AT(neighbour);
      ibset const &cfboxes = local_box.restriction_boundary.AT(neighbour);

      // Set up a bit mask that keeps the lower (when dir[d]=-1) or
      // upper (when dir[d]=+1) half of the bits in each direction d
      unsigned const bits = BMSK(dim);
      unsigned bmask = 0;
      for (int d = 0; d < dim; ++d) {
        for (unsigned b = 0; b < bits; ++b) {
          if ((shift[d] == -1 and BGET(b, d) == 0) or
              (shift[d] == +1 and BGET(b, d) != 0)) {
            bmask = BSET(bmask, b);
          }
        }
      }

      // Handle prolongation region (region that is prolongated from
      // the next coarser level)

      if (verbose) {
        ostringstream buf;
        buf << "Setting boundary " << shift << ": level " << reflevel
            << " prolongation region " << boxes << " to bmask 0x" << hex
            << bmask << dec;
        CCTK_INFO(buf.str().c_str());
      }

      LOOP_OVER_BSET(cctkGH, boxes, box, imin, imax) {
        if (verbose) {
          ostringstream buf;
          buf << "  Setting prolongation box " << imin << ":"
              << imax - ivect(1);
          CCTK_INFO(buf.str().c_str());
        }
#pragma omp parallel
        CCTK_LOOP3(CarpetMaskSetup_prolongation, i, j, k, imin[0], imin[1],
                   imin[2], imax[0], imax[1], imax[2], cctk_ash[0], cctk_ash[1],
                   cctk_ash[2]) {
          int const ind = CCTK_GFINDEX3D(cctkGH, i, j, k);
          iweight[ind] &= bmask;
        }
        CCTK_ENDLOOP3(CarpetMaskSetup_prolongation);
      }
      END_LOOP_OVER_BSET;

      // Handle restriction region (region that is restricted from the
      // next finer level)

      if (verbose) {
        ostringstream buf;
        buf << "Setting boundary " << shift << ": level " << reflevel
            << " restriction region " << cfboxes << " to bmask 0x" << hex
            << bmask << dec;
        CCTK_INFO(buf.str().c_str());
      }

      LOOP_OVER_BSET(cctkGH, cfboxes, box, imin, imax) {
        if (verbose) {
          ostringstream buf;
          buf << "  Setting restriction box " << imin << ":" << imax - ivect(1);
          CCTK_INFO(buf.str().c_str());
        }
#pragma omp parallel
        CCTK_LOOP3(CarpetMaskSetup_restriction, i, j, k, imin[0], imin[1],
                   imin[2], imax[0], imax[1], imax[2], cctk_ash[0], cctk_ash[1],
                   cctk_ash[2]) {
          int const ind = CCTK_GFINDEX3D(cctkGH, i, j, k);
          iweight[ind] &= bmask;
        }
        CCTK_ENDLOOP3(CarpetMaskSetup_restriction);
      }
      END_LOOP_OVER_BSET;

    } // for neighbours
  }
  END_LOCAL_COMPONENT_LOOP;
}

} // namespace CarpetMask
