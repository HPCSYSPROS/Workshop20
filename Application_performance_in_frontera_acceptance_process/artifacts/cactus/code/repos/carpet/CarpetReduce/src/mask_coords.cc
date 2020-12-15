#include <assert.h>
#include <stdlib.h>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Functions.h>
#include <cctk_Parameters.h>

#include <loopcontrol.h>

#include "bits.h"

extern "C" void CoordBase_SetupMask(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT nboundaryzones[2 * cctk_dim];
  CCTK_INT is_internal[2 * cctk_dim];
  CCTK_INT is_staggered[2 * cctk_dim];
  CCTK_INT shiftout[2 * cctk_dim];
  CCTK_INT is_periodic[2 * cctk_dim];

  int bnd_points[2 * cctk_dim]; /* points outside the domain */
  int int_points[cctk_dim];     /* global interior points */

  int gmin[cctk_dim], gmax[cctk_dim]; /* global domain extent */
  int imin[cctk_dim], imax[cctk_dim]; /* domain extent */
  int bmin[cctk_dim], bmax[cctk_dim]; /* boundary extent */

  int ierr;

  int const reflevel = GetRefinementLevel(cctkGH);

  if (CCTK_IsFunctionAliased("MultiPatch_GetBoundarySpecification")) {
    int const m = MultiPatch_GetMap(cctkGH);
    assert(m >= 0);
    ierr = MultiPatch_GetBoundarySpecification(
        m, 2 * cctk_dim, nboundaryzones, is_internal, is_staggered, shiftout);
  } else {
    ierr = GetBoundarySpecification(2 * cctk_dim, nboundaryzones, is_internal,
                                    is_staggered, shiftout);
  }
  if (ierr != 0) {
    CCTK_ERROR("Could not get boundary specification");
  }

  for (int d = 0; d < cctk_dim; ++d) {
    for (int f = 0; f < 2; ++f) {
      is_periodic[2 * d + f] = 0;
    }
  }
  CCTK_INT const *const periodic = static_cast<CCTK_INT const *>(
      CCTK_ParameterGet("periodic", "PeriodicCarpet", NULL));
  CCTK_INT const *const periodic_x = static_cast<CCTK_INT const *>(
      CCTK_ParameterGet("periodic_x", "PeriodicCarpet", NULL));
  CCTK_INT const *const periodic_y = static_cast<CCTK_INT const *>(
      CCTK_ParameterGet("periodic_x", "PeriodicCarpet", NULL));
  CCTK_INT const *const periodic_z = static_cast<CCTK_INT const *>(
      CCTK_ParameterGet("periodic_x", "PeriodicCarpet", NULL));
  is_periodic[0] = is_periodic[1] =
      (periodic && *periodic) || (periodic_x && *periodic_x);
  is_periodic[2] = is_periodic[3] =
      (periodic && *periodic) || (periodic_y && *periodic_y);
  is_periodic[4] = is_periodic[5] =
      (periodic && *periodic) || (periodic_z && *periodic_z);

  /* Calculate the number of boundary points. This excludes points
     that are directly on the boundary. */
  for (int d = 0; d < cctk_dim; ++d) {
    for (int f = 0; f < 2; ++f) {

      if (is_internal[2 * d + f]) {
        /* The boundary extends inwards */
        bnd_points[2 * d + f] = shiftout[2 * d + f];
      } else {
        /* The boundary extends outwards */
        if (is_periodic[2 * d + f]) {
          bnd_points[2 * d + f] = nboundaryzones[2 * d + f];
        } else {
          bnd_points[2 * d + f] = nboundaryzones[2 * d + f] +
                                  shiftout[2 * d + f] +
                                  is_staggered[2 * d + f] - 1;
        }
      }
    }
  }

  /* Calculate the global extent of the domain */
  for (int d = 0; d < cctk_dim; ++d) {
    gmin[d] = 0;
    gmax[d] = cctk_gsh[d];
  }

  /* Ensure that the domain specification is consistent */
  for (int d = 0; d < cctk_dim; ++d) {
    int_points[d] =
        (gmax[d] - bnd_points[2 * d + 1]) - (gmin[d] + bnd_points[2 * d]);
    if (int_points[d] < 0) {
      CCTK_ERROR("Number of internal grid points is negative");
    }
  }

  /* Loop over all dimensions and faces */
  for (int d = 0; d < cctk_dim; ++d) {
    for (int f = 0; f < 2; ++f) {

      /* If this processor has the outer boundary */
      if (cctk_bbox[2 * d + f]) {

        if (bnd_points[2 * d + f] < 0 || bnd_points[2 * d + f] > cctk_lsh[d]) {
          CCTK_ERROR("Illegal number of boundary points");
        }

        /* Calculate the extent of the local part of the domain */
        for (int dd = 0; dd < cctk_dim; ++dd) {
          imin[dd] = 0;
          imax[dd] = cctk_lsh[dd];
        }

        /* Calculate the extent of the boundary */
        for (int dd = 0; dd < cctk_dim; ++dd) {
          bmin[dd] = imin[dd];
          bmax[dd] = imax[dd];
        }
        if (f == 0) {
          /* lower face */
          bmax[d] = imin[d] + bnd_points[2 * d + f];
        } else {
          /* upper face */
          bmin[d] = imax[d] - bnd_points[2 * d + f];
        }

        /* Loop over the boundary */
        if (verbose) {
          CCTK_VInfo(CCTK_THORNSTRING, "Setting boundary points in direction "
                                       "%d face %d to weight 0 on level %d",
                     d, f, reflevel);
        }
#pragma omp parallel
        CCTK_LOOP3(CoordBase_SetupMask_boundary, i, j, k, bmin[0], bmin[1],
                   bmin[2], bmax[0], bmax[1], bmax[2], cctk_ash[0], cctk_ash[1],
                   cctk_ash[2]) {
          int const ind = CCTK_GFINDEX3D(cctkGH, i, j, k);
          iweight[ind] = 0;
        }
        CCTK_ENDLOOP3(CoordBase_SetupMask_boundary);

        /* When the boundary is neither staggered nor periodic, then
           give the points directly on the boundary the weight 1/2 */
        if (!is_staggered[2 * d + f] && !is_periodic[2 * d + f]) {

          /* Since we are going to cut off 1/2 at each boundary, the
             domain size must be at least 1 to begin with */
          assert(int_points[d] >= 1);

          /* Check whether the domain is empty */
          if (int_points[d] == 1) {

            /* The domain is empty. The correct thing to do would be
               to set the weights to 0. But this is boring, because
               then there are no interior points left, and all
               reductions become trivial. Instead, leave the
               non-staggered boundary points alone, and assume that
               the user wants to perform a reduction in one fewer
               dimensions. */
            /* do nothing */

          } else {

            /* Adapt the extent of the boundary */
            if (f == 0) {
              /* lower face */
              bmin[d] = bmax[d];
              bmax[d] = bmin[d] + 1;
            } else {
              /* upper face */
              bmax[d] = bmin[d];
              bmin[d] = bmax[d] - 1;
            }

            /* Loop over the points next to boundary */
            /* Set up a bit mask that keeps the upper (f==0) or lower (f==1)
               half of the bits in each direction d */
            unsigned const bits = BMSK(cctk_dim);
            unsigned bmask = 0;
            for (unsigned b = 0; b < bits; ++b) {
              if (BGET(b, d) == !f) {
                bmask = BSET(bmask, b);
              }
            }
            assert(BCNT(bmask) == bits / 2);
            if (verbose) {
              CCTK_VInfo(CCTK_THORNSTRING, "Setting non-staggered boundary "
                                           "points on level %d in direction %d "
                                           "face %d to bmask 0x%x",
                         reflevel, d, f, bmask);
            }
#pragma omp parallel
            CCTK_LOOP3(CoordBase_SetupMask_boundary2, i, j, k, bmin[0], bmin[1],
                       bmin[2], bmax[0], bmax[1], bmax[2], cctk_ash[0],
                       cctk_ash[1], cctk_ash[2]) {
              int const ind = CCTK_GFINDEX3D(cctkGH, i, j, k);
              iweight[ind] &= bmask;
            }
            CCTK_ENDLOOP3(CoordBase_SetupMask_boundary2);

          } /* if the domain is not empty */

        } /* if the boundary is not staggered */

      } else { /* if this is a ghost boundary */

        /* Calculate the extent of the local part of the domain */
        for (int dd = 0; dd < cctk_dim; ++dd) {
          imin[dd] = 0;
          imax[dd] = cctk_lsh[dd];
        }

        /* Calculate the extent of the boundary */
        for (int dd = 0; dd < cctk_dim; ++dd) {
          bmin[dd] = imin[dd];
          bmax[dd] = imax[dd];
        }
        if (f == 0) {
          /* lower face */
          bmax[d] = imin[d] + cctk_nghostzones[d];
        } else {
          /* upper face */
          bmin[d] = imax[d] - cctk_nghostzones[d];
        }

        /* Loop over the boundary */
        if (verbose) {
          CCTK_VInfo(CCTK_THORNSTRING, "Setting ghost points in direction %d "
                                       "face %d to weight 0 on level %d",
                     d, f, reflevel);
        }
#pragma omp parallel
        CCTK_LOOP3(CoordBase_SetupMask_ghost, i, j, k, bmin[0], bmin[1],
                   bmin[2], bmax[0], bmax[1], bmax[2], cctk_ash[0], cctk_ash[1],
                   cctk_ash[2]) {
          int const ind = CCTK_GFINDEX3D(cctkGH, i, j, k);
          iweight[ind] = 0;
        }
        CCTK_ENDLOOP3(CoordBase_SetupMask_ghost);

      } /* if this is a ghost boundary */

    } /* loop over faces */
  }   /* loop over directions */
}
