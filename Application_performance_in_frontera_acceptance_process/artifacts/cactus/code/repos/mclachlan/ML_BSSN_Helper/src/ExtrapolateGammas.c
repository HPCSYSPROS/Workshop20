#include <cctk.h>
#include <cctk_Arguments.h>

#include <stdbool.h>

static void extrap(const cGH *cctkGH, CCTK_REAL *var);

void ML_BSSN_ExtrapolateGammas(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;

  extrap(cctkGH, Xt1);
  extrap(cctkGH, Xt2);
  extrap(cctkGH, Xt3);

  static bool have_A_index = false;
  static int A_index;
  if (!have_A_index) {
    A_index = CCTK_VarIndex("ML_BSSN::A");
    have_A_index = true;
  }
  if (A_index >= 0) {
    CCTK_REAL *A_ptr = CCTK_VarDataPtrI(cctkGH, 0, A_index);
    extrap(cctkGH, A_ptr);
  }

  static bool have_B_index = false;
  static int B_index;
  if (!have_B_index) {
    B_index = CCTK_VarIndex("ML_BSSN::B1");
    have_B_index = true;
  }
  if (B_index >= 0) {
    CCTK_REAL *B1_ptr = CCTK_VarDataPtrI(cctkGH, 0, B_index);
    CCTK_REAL *B2_ptr = CCTK_VarDataPtrI(cctkGH, 0, B_index + 1);
    CCTK_REAL *B3_ptr = CCTK_VarDataPtrI(cctkGH, 0, B_index + 2);
    extrap(cctkGH, B1_ptr);
    extrap(cctkGH, B2_ptr);
    extrap(cctkGH, B3_ptr);
  }
}

static void extrap(const cGH *cctkGH, CCTK_REAL *var) {
  ExtrapolateGammas(cctkGH, var);
}
