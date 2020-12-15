#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cmath>

#include <carpet.hh>
#include <mpi_string.hh>

#include "boundary.hh"

namespace CarpetRegrid2 {

using namespace std;
using namespace Carpet;

void evaluate_level_mask(cGH const *restrict const cctkGH,
                         vector<ibset> &regions, int const rl) {
  DECLARE_CCTK_PARAMETERS;

  assert(is_singlemap_mode());
  int const m = Carpet::map;
  assert(rl > 0);
  BEGIN_GLOBAL_MODE(cctkGH) {
    ENTER_LEVEL_MODE(cctkGH, rl - 1) {
      ENTER_SINGLEMAP_MODE(cctkGH, m, CCTK_GF) {
        if (verbose or veryverbose) {
          cout << "Refinement level " << reflevel << ":\n";
        }

        gh const &hh = *vhh.AT(Carpet::map);

        ivect const gsh = ivect::ref(cctkGH->cctk_gsh);
        ivect const nghostzones = ivect::ref(cctkGH->cctk_nghostzones);
        if (veryverbose) {
          cout << "   gsh: " << gsh << "\n"
               << "   nghostzones: " << nghostzones << "\n";
        }

        // Determine the block size
        ivect block_size = adaptive_block_size;
        if (adaptive_block_size_x > 0)
          block_size[0] = adaptive_block_size_x;
        if (adaptive_block_size_y > 0)
          block_size[1] = adaptive_block_size_y;
        if (adaptive_block_size_z > 0)
          block_size[2] = adaptive_block_size_z;
        assert(all(block_size > 0));

        // Determine the block offset, i.e. the grid point index of
        // the coordinate origin. We don't start blocks at index
        // zero, so that one can change the location of the outer
        // boundary without influencing the refinement mechanism.
        // Instead, we align the block such that one of the blocks
        // begins at the coordinate origin.
        ivect origin;
        rvect rorigin;
        {
          DECLARE_CCTK_ARGUMENTS;
          for (int d = 0; d < dim; ++d) {
            // We round up, because the origin is located between
            // grid points for cell-centred grids, and we want the
            // next upper grid point in this case. We subtract 0.25
            // just for proper rounding.
            origin[d] =
                ceil(-CCTK_ORIGIN_SPACE(d) / CCTK_DELTA_SPACE(d) - 0.25);
            rorigin[d] = CCTK_ORIGIN_SPACE(d) + origin[d] * CCTK_DELTA_SPACE(d);
          }
        }
        // block_index := ([0...gsh] + block_offset) / block_size
        // [0...gsh] = block_index * block_size - block_offset
        ivect const block_offset =
            (block_size - origin % block_size) % block_size;
        ivect const num_blocks =
            (gsh + block_offset + block_size - 1) / block_size;
        if (veryverbose) {
          cout << "   origin: " << origin << "\n"
               << "   coordinates of origin: " << rorigin << "\n"
               << "   offset: " << block_offset << "\n"
               << "   num blocks: " << num_blocks << "\n";
        }
        assert(all(num_blocks > 0));

        // Mask containing all refined block indices (in arbitrary
        // order, and potentially even duplicated). This is really a
        // set, but we need to communicate this data structure, and
        // therefore a vector is a better choice.
        vector<ivect> mask;

        // For vertex centred grids, the blocks need to overlap by
        // one grid point, so that e.g. an 8^3 block is turned into
        // a 9^3 block. This will ensure symmetry of refined regions
        // in the simulation domain.
        int const overlap = hh.refcent == vertex_centered ? 1 : 0;

        BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, CCTK_GF) {
          DECLARE_CCTK_ARGUMENTS;

          if (veryverbose) {
            cout << "   component " << component << ":\n";
          }

          ivect const lbnd = ivect::ref(cctk_lbnd);
          ivect const lsh = ivect::ref(cctk_lsh);
          ivect const bboxlo(cctk_bbox[0], cctk_bbox[2], cctk_bbox[4]);
          ivect const bboxhi(cctk_bbox[1], cctk_bbox[3], cctk_bbox[5]);

          ivect const imin = 0 + either(bboxlo, 0, nghostzones);
          ivect const imax = lsh - either(bboxhi, 0, nghostzones);

          ivect const bmin =
              max(0, ((lbnd + imin + block_offset + block_size - overlap) /
                      block_size) -
                         1);
          ivect const bmax = (lbnd + imax + block_offset) / block_size;

          // Loop over all blocks
          int nblocks = 0;
          int nrefined = 0;
          for (int bk = bmin[2]; bk < bmax[2]; ++bk) {
            for (int bj = bmin[1]; bj < bmax[1]; ++bj) {
              for (int bi = bmin[0]; bi < bmax[0]; ++bi) {
                ivect const bind(bi, bj, bk);

                ivect const bimin =
                    max(imin, (bind)*block_size - block_offset - lbnd);
                ivect const bimax =
                    min(imax, (bind + 1) * block_size - block_offset - lbnd +
                                  overlap);

                bool refine = false;
                bool have_nan = false;

                // Loop over all points in this block
                for (int k = bimin[2]; k < bimax[2]; ++k) {
                  for (int j = bimin[1]; j < bimax[1]; ++j) {
                    for (int i = bimin[0]; i < bimax[0]; ++i) {
                      int const ind3d = CCTK_GFINDEX3D(cctkGH, i, j, k);
                      if (isfinite(level_mask[ind3d])) {
                        refine = refine or level_mask[ind3d] >= rl;
                      } else {
                        have_nan = true;
                      }
                    }
                  }
                }

                // Refine this block if any point in this block requires
                // refinement
                if (have_nan) {
                  cout << "      *** found nan in block " << bind << " ***\n";
                }
                if (refine) {
                  if (veryverbose) {
                    cout << "      refining block " << bind << "\n";
                  }
                  mask.push_back(bind);
                  ++nrefined;
                }
                ++nblocks;
              }
            }
          }

          if (veryverbose) {
            cout << "   refining " << nrefined << " "
                 << "of " << nblocks << " blocks on this component\n";
          }
        }
        END_LOCAL_COMPONENT_LOOP;

        // Combine this mask from all processes, and to all processes
        vector<ivect> const fullmask = allgatherv1(dist::comm(), mask);

        // Convert vector into ibset
        const ibbox &cbaseext = hh.baseextent(mglevel, rl - 1);
        const ibbox &baseext = hh.baseextent(mglevel, rl);

        ibset &region = regions.at(rl);
        for (vector<ivect>::const_iterator it = fullmask.begin();
             it != fullmask.end(); ++it) {
          ivect const &bind = *it;
          ivect const cstr = cbaseext.stride();
          ivect const clo =
              cbaseext.lower() + (bind * block_size - block_offset) * cstr;
          ivect const cup = clo + (block_size + overlap - 1) * cstr;
          ibbox const cblock(clo, cup, cstr);
          ibbox const block = cblock.expanded_for(baseext);
          region |= block;
        }

        if (verbose or veryverbose) {
          if (veryverbose) {
            cout << "   refined region is " << region << "\n"
                 << "   domain is " << baseext << "\n"
                 << "   next coarser level is " << regions.at(rl - 1) << "\n";
          }
          size_type const nrefined = region.size();
          size_type const npoints = baseext.size();
          size_type const ncoarse =
              regions.at(rl - 1).expanded_for(baseext).size();
          cout << "   refining " << (100.0 * nrefined / ncoarse)
               << "% of the next coarser level\n"
               << "   refining " << (100.0 * nrefined / npoints)
               << "% of the domain\n";
        }
      }
      LEAVE_SINGLEMAP_MODE;
    }
    LEAVE_LEVEL_MODE;
  }
  END_GLOBAL_MODE;
}

} // namespace CarpetRegrid2
