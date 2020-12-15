#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <carpet.h>

#include <assert.h>
#include <math.h>

extern "C" void MaskBase_TestMask(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (verbose) {
    CCTK_INFO("Testing weight");
  }

  int const sum = CCTK_ReductionHandle("sum");
  assert(sum >= 0);

  int const proc = 0;

  int const weight_var = CCTK_VarIndex("CarpetReduce::weight");
  int const one_var = CCTK_VarIndex("CarpetReduce::one");
  assert(weight_var >= 0);
  assert(one_var >= 0);

  CCTK_REAL sum_weight, all_excised_cells;

  {
    int const ierr = CCTK_Reduce(cctkGH, proc, sum, 1, CCTK_VARIABLE_REAL,
                                 &sum_weight, 1, one_var);
    assert(ierr >= 0);
  }
  {
    int const ierr =
        CCTK_ReduceLocalScalar(cctkGH, proc, sum, excised_cells,
                               &all_excised_cells, CCTK_VARIABLE_REAL);
    assert(ierr >= 0);
    *excised_cells = all_excised_cells;
  }

  if (proc == -1 || CCTK_MyProc(cctkGH) == proc) {

    if (verbose) {
      CCTK_VInfo(CCTK_THORNSTRING, "Reduction weight sum: %.17g",
                 (double)sum_weight);
    }

    CCTK_REAL domain_volume = 0.0;

    int maps = 1;
    if (CCTK_IsFunctionAliased("MultiPatch_GetMaps")) {
      maps = MultiPatch_GetMaps(cctkGH);
    }

    for (int m = 0; m < maps; ++m) {

      CCTK_REAL physical_min[cctk_dim];
      CCTK_REAL physical_max[cctk_dim];
      CCTK_REAL interior_min[cctk_dim];
      CCTK_REAL interior_max[cctk_dim];
      CCTK_REAL exterior_min[cctk_dim];
      CCTK_REAL exterior_max[cctk_dim];
      CCTK_REAL spacing[cctk_dim];

      /* TODO: This requires that the domain was set up via CoordBase, not via
       * grid. Check this, and don't output any warnings if so. */
      if (CCTK_IsFunctionAliased("MultiPatch_GetDomainSpecification")) {
        int const ierr = MultiPatch_GetDomainSpecification(
            m, cctk_dim, physical_min, physical_max, interior_min, interior_max,
            exterior_min, exterior_max, spacing);
        assert(!ierr);
      } else {
        int const ierr = GetDomainSpecification(
            cctk_dim, physical_min, physical_max, interior_min, interior_max,
            exterior_min, exterior_max, spacing);
        assert(!ierr);
      }

      CCTK_REAL map_volume = 1.0;
      for (int d = 0; d < cctk_dim; ++d) {
        map_volume *= (physical_max[d] - physical_min[d]) / spacing[d];
      }

      if (verbose) {
        CCTK_VInfo(CCTK_THORNSTRING, "Volume of map #%d: %.17g", m,
                   (double)map_volume);
      }

      domain_volume += map_volume;

    } /* for m */

    /* Don't count excised regions towards the domain volume */
    domain_volume -= *excised_cells;

    if (verbose) {
      CCTK_VInfo(CCTK_THORNSTRING, "Simulation domain volume:  %.17g",
                 (double)domain_volume);
      CCTK_VInfo(CCTK_THORNSTRING, "Additional excised volume: %.17g",
                 (double)*excised_cells);
    }

    int const there_is_a_problem = fabs(sum_weight - domain_volume) >
                                   1.0e-12 * (sum_weight + domain_volume);

    if (there_is_a_problem) {
      if (!verbose) {
        CCTK_VInfo(CCTK_THORNSTRING, "Simulation domain volume:  %.17g",
                   (double)domain_volume);
        CCTK_VInfo(CCTK_THORNSTRING, "Additional excised volume: %.17g",
                   (double)*excised_cells);
        CCTK_VInfo(CCTK_THORNSTRING, "Reduction weight sum:      %.17g",
                   (double)sum_weight);
      }
      CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Simulation domain volume and reduction weight sum differ");
    }
  }

  if (!debug_iweight) {
    CCTK_DisableGroupStorage(cctkGH, "CarpetReduce::iweight");
    CCTK_DisableGroupStorage(cctkGH, "CarpetReduce::one");
  }
}
