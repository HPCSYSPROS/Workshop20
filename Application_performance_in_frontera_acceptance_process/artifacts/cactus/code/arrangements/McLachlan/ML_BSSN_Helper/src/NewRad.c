#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <math.h>
#include <stdbool.h>

static void newrad(const cGH *cctkGH, const CCTK_REAL *var, CCTK_REAL *rhs,
                   CCTK_REAL var0, CCTK_REAL v0);

void ML_BSSN_NewRad(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_REAL v0 = sqrt(harmonicF);

  newrad(cctkGH, phi, phirhs, conformalMethod ? 1.0 : 0.0, v0);

  newrad(cctkGH, gt11, gt11rhs, 1.0, 1.0);
  newrad(cctkGH, gt12, gt12rhs, 0.0, 1.0);
  newrad(cctkGH, gt13, gt13rhs, 0.0, 1.0);
  newrad(cctkGH, gt22, gt22rhs, 1.0, 1.0);
  newrad(cctkGH, gt23, gt23rhs, 0.0, 1.0);
  newrad(cctkGH, gt33, gt33rhs, 1.0, 1.0);

  newrad(cctkGH, Xt1, Xt1rhs, 0.0, 1.0);
  newrad(cctkGH, Xt2, Xt2rhs, 0.0, 1.0);
  newrad(cctkGH, Xt3, Xt3rhs, 0.0, 1.0);

  static bool have_Theta_index = false;
  static int Theta_index, Thetarhs_index;
  if (!have_Theta_index) {
    Theta_index = CCTK_VarIndex("ML_BSSN::Theta");
    Thetarhs_index = CCTK_VarIndex("ML_BSSN::Thetarhs");
    have_Theta_index = true;
  }
  if (Theta_index >= 0) {
    CCTK_REAL *Theta_ptr = CCTK_VarDataPtrI(cctkGH, 0, Theta_index);
    CCTK_REAL *Thetarhs_ptr = CCTK_VarDataPtrI(cctkGH, 0, Thetarhs_index);
    newrad(cctkGH, Theta_ptr, Thetarhs_ptr, 0.0, v0);
  }

  newrad(cctkGH, trK, trKrhs, 0.0, v0);

  newrad(cctkGH, At11, At11rhs, 0.0, 1.0);
  newrad(cctkGH, At12, At12rhs, 0.0, 1.0);
  newrad(cctkGH, At13, At13rhs, 0.0, 1.0);
  newrad(cctkGH, At22, At22rhs, 0.0, 1.0);
  newrad(cctkGH, At23, At23rhs, 0.0, 1.0);
  newrad(cctkGH, At33, At33rhs, 0.0, 1.0);

  newrad(cctkGH, alpha, alpharhs, 1.0, v0);

  static bool have_A_index = false;
  static int A_index, Arhs_index;
  if (!have_A_index) {
    A_index = CCTK_VarIndex("ML_BSSN::A");
    Arhs_index = CCTK_VarIndex("ML_BSSN::Arhs");
    have_A_index = true;
  }
  if (A_index >= 0) {
    CCTK_REAL *A_ptr = CCTK_VarDataPtrI(cctkGH, 0, A_index);
    CCTK_REAL *Arhs_ptr = CCTK_VarDataPtrI(cctkGH, 0, Arhs_index);
    newrad(cctkGH, A_ptr, Arhs_ptr, 0.0, v0);
  }

  newrad(cctkGH, beta1, beta1rhs, 0.0, 1.0);
  newrad(cctkGH, beta2, beta2rhs, 0.0, 1.0);
  newrad(cctkGH, beta3, beta3rhs, 0.0, 1.0);

  static bool have_B_index = false;
  static int B_index, Brhs_index;
  if (!have_B_index) {
    B_index = CCTK_VarIndex("ML_BSSN::B1");
    Brhs_index = CCTK_VarIndex("ML_BSSN::B1rhs");
    have_B_index = true;
  }
  if (B_index >= 0) {
    CCTK_REAL *B1_ptr = CCTK_VarDataPtrI(cctkGH, 0, B_index);
    CCTK_REAL *B2_ptr = CCTK_VarDataPtrI(cctkGH, 0, B_index + 1);
    CCTK_REAL *B3_ptr = CCTK_VarDataPtrI(cctkGH, 0, B_index + 2);
    CCTK_REAL *B1rhs_ptr = CCTK_VarDataPtrI(cctkGH, 0, Brhs_index);
    CCTK_REAL *B2rhs_ptr = CCTK_VarDataPtrI(cctkGH, 0, Brhs_index + 1);
    CCTK_REAL *B3rhs_ptr = CCTK_VarDataPtrI(cctkGH, 0, Brhs_index + 2);
    newrad(cctkGH, B1_ptr, B1rhs_ptr, 0.0, 1.0);
    newrad(cctkGH, B2_ptr, B2rhs_ptr, 0.0, 1.0);
    newrad(cctkGH, B3_ptr, B3rhs_ptr, 0.0, 1.0);
  }
}

static void newrad(const cGH *cctkGH, const CCTK_REAL *var, CCTK_REAL *rhs,
                   CCTK_REAL var0, CCTK_REAL v0) {
  DECLARE_CCTK_PARAMETERS;

  NewRad_Apply(cctkGH, var, rhs, var0, v0, radpower);
}
