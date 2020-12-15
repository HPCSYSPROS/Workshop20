/*  File produced by Kranc */

#define KRANC_C

#include <algorithm>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Kranc.hh"
#include "Differencing.h"
#include "loopcontrol.h"
#include "vectors.h"

namespace ML_BSSN {

extern "C" void ML_BSSN_InitialADMBase2Interior_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (cctk_iteration % ML_BSSN_InitialADMBase2Interior_calc_every != ML_BSSN_InitialADMBase2Interior_calc_offset)
    return;
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::ML_dtlapse","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for ML_BSSN::ML_dtlapse.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::ML_dtshift","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for ML_BSSN::ML_dtshift.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::ML_Gamma","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for ML_BSSN::ML_Gamma.");
  return;
}

static void ML_BSSN_InitialADMBase2Interior_Body(const cGH* restrict const cctkGH, const int dir, const int face, const CCTK_REAL normal[3], const CCTK_REAL tangentA[3], const CCTK_REAL tangentB[3], const int imin[3], const int imax[3], const int n_subblock_gfs, CCTK_REAL* restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  /* Include user-supplied include files */
  /* Initialise finite differencing variables */
  const ptrdiff_t di CCTK_ATTRIBUTE_UNUSED = 1;
  const ptrdiff_t dj CCTK_ATTRIBUTE_UNUSED = 
    CCTK_GFINDEX3D(cctkGH,0,1,0) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  const ptrdiff_t dk CCTK_ATTRIBUTE_UNUSED = 
    CCTK_GFINDEX3D(cctkGH,0,0,1) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL) * di;
  const ptrdiff_t cdj CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL) * dj;
  const ptrdiff_t cdk CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL) * dk;
  const ptrdiff_t cctkLbnd1 CCTK_ATTRIBUTE_UNUSED = cctk_lbnd[0];
  const ptrdiff_t cctkLbnd2 CCTK_ATTRIBUTE_UNUSED = cctk_lbnd[1];
  const ptrdiff_t cctkLbnd3 CCTK_ATTRIBUTE_UNUSED = cctk_lbnd[2];
  const CCTK_REAL_VEC t CCTK_ATTRIBUTE_UNUSED = ToReal(cctk_time);
  const CCTK_REAL_VEC cctkOriginSpace1 CCTK_ATTRIBUTE_UNUSED = 
    ToReal(CCTK_ORIGIN_SPACE(0));
  const CCTK_REAL_VEC cctkOriginSpace2 CCTK_ATTRIBUTE_UNUSED = 
    ToReal(CCTK_ORIGIN_SPACE(1));
  const CCTK_REAL_VEC cctkOriginSpace3 CCTK_ATTRIBUTE_UNUSED = 
    ToReal(CCTK_ORIGIN_SPACE(2));
  const CCTK_REAL_VEC dt CCTK_ATTRIBUTE_UNUSED = 
    ToReal(CCTK_DELTA_TIME);
  const CCTK_REAL_VEC dx CCTK_ATTRIBUTE_UNUSED = 
    ToReal(CCTK_DELTA_SPACE(0));
  const CCTK_REAL_VEC dy CCTK_ATTRIBUTE_UNUSED = 
    ToReal(CCTK_DELTA_SPACE(1));
  const CCTK_REAL_VEC dz CCTK_ATTRIBUTE_UNUSED = 
    ToReal(CCTK_DELTA_SPACE(2));
  const CCTK_REAL_VEC dxi CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1),dx);
  const CCTK_REAL_VEC dyi CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1),dy);
  const CCTK_REAL_VEC dzi CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1),dz);
  const CCTK_REAL_VEC khalf CCTK_ATTRIBUTE_UNUSED = ToReal(0.5);
  const CCTK_REAL_VEC kthird CCTK_ATTRIBUTE_UNUSED = 
    ToReal(0.333333333333333333333333333333);
  const CCTK_REAL_VEC ktwothird CCTK_ATTRIBUTE_UNUSED = 
    ToReal(0.666666666666666666666666666667);
  const CCTK_REAL_VEC kfourthird CCTK_ATTRIBUTE_UNUSED = 
    ToReal(1.33333333333333333333333333333);
  const CCTK_REAL_VEC hdxi CCTK_ATTRIBUTE_UNUSED = 
    kmul(dxi,ToReal(0.5));
  const CCTK_REAL_VEC hdyi CCTK_ATTRIBUTE_UNUSED = 
    kmul(dyi,ToReal(0.5));
  const CCTK_REAL_VEC hdzi CCTK_ATTRIBUTE_UNUSED = 
    kmul(dzi,ToReal(0.5));
  /* Initialize predefined quantities */
  const CCTK_REAL_VEC p1o1024dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0009765625),dx);
  const CCTK_REAL_VEC p1o1024dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0009765625),dy);
  const CCTK_REAL_VEC p1o1024dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0009765625),dz);
  const CCTK_REAL_VEC p1o120dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00833333333333333333333333333333),dx);
  const CCTK_REAL_VEC p1o120dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00833333333333333333333333333333),dy);
  const CCTK_REAL_VEC p1o120dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00833333333333333333333333333333),dz);
  const CCTK_REAL_VEC p1o12dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0833333333333333333333333333333),dx);
  const CCTK_REAL_VEC p1o12dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0833333333333333333333333333333),dy);
  const CCTK_REAL_VEC p1o12dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0833333333333333333333333333333),dz);
  const CCTK_REAL_VEC p1o144dxdy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00694444444444444444444444444444),kmul(dx,dy));
  const CCTK_REAL_VEC p1o144dxdz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00694444444444444444444444444444),kmul(dx,dz));
  const CCTK_REAL_VEC p1o144dydz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00694444444444444444444444444444),kmul(dy,dz));
  const CCTK_REAL_VEC p1o1680dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.000595238095238095238095238095238),dx);
  const CCTK_REAL_VEC p1o1680dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.000595238095238095238095238095238),dy);
  const CCTK_REAL_VEC p1o1680dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.000595238095238095238095238095238),dz);
  const CCTK_REAL_VEC p1o180dx2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00555555555555555555555555555556),kmul(dx,dx));
  const CCTK_REAL_VEC p1o180dy2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00555555555555555555555555555556),kmul(dy,dy));
  const CCTK_REAL_VEC p1o180dz2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00555555555555555555555555555556),kmul(dz,dz));
  const CCTK_REAL_VEC p1o24dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0416666666666666666666666666667),dx);
  const CCTK_REAL_VEC p1o24dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0416666666666666666666666666667),dy);
  const CCTK_REAL_VEC p1o24dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0416666666666666666666666666667),dz);
  const CCTK_REAL_VEC p1o2dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.5),dx);
  const CCTK_REAL_VEC p1o2dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.5),dy);
  const CCTK_REAL_VEC p1o2dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.5),dz);
  const CCTK_REAL_VEC p1o3600dxdy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.000277777777777777777777777777778),kmul(dx,dy));
  const CCTK_REAL_VEC p1o3600dxdz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.000277777777777777777777777777778),kmul(dx,dz));
  const CCTK_REAL_VEC p1o3600dydz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.000277777777777777777777777777778),kmul(dy,dz));
  const CCTK_REAL_VEC p1o4dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.25),dx);
  const CCTK_REAL_VEC p1o4dxdy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.25),kmul(dx,dy));
  const CCTK_REAL_VEC p1o4dxdz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.25),kmul(dx,dz));
  const CCTK_REAL_VEC p1o4dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.25),dy);
  const CCTK_REAL_VEC p1o4dydz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.25),kmul(dy,dz));
  const CCTK_REAL_VEC p1o4dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.25),dz);
  const CCTK_REAL_VEC p1o5040dx2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.000198412698412698412698412698413),kmul(dx,dx));
  const CCTK_REAL_VEC p1o5040dy2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.000198412698412698412698412698413),kmul(dy,dy));
  const CCTK_REAL_VEC p1o5040dz2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.000198412698412698412698412698413),kmul(dz,dz));
  const CCTK_REAL_VEC p1o560dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00178571428571428571428571428571),dx);
  const CCTK_REAL_VEC p1o560dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00178571428571428571428571428571),dy);
  const CCTK_REAL_VEC p1o560dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00178571428571428571428571428571),dz);
  const CCTK_REAL_VEC p1o60dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0166666666666666666666666666667),dx);
  const CCTK_REAL_VEC p1o60dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0166666666666666666666666666667),dy);
  const CCTK_REAL_VEC p1o60dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0166666666666666666666666666667),dz);
  const CCTK_REAL_VEC p1o64dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.015625),dx);
  const CCTK_REAL_VEC p1o64dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.015625),dy);
  const CCTK_REAL_VEC p1o64dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.015625),dz);
  const CCTK_REAL_VEC p1o6dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.166666666666666666666666666667),dx);
  const CCTK_REAL_VEC p1o6dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.166666666666666666666666666667),dy);
  const CCTK_REAL_VEC p1o6dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.166666666666666666666666666667),dz);
  const CCTK_REAL_VEC p1o705600dxdy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1.41723356009070294784580498866e-6),kmul(dx,dy));
  const CCTK_REAL_VEC p1o705600dxdz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1.41723356009070294784580498866e-6),kmul(dx,dz));
  const CCTK_REAL_VEC p1o705600dydz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1.41723356009070294784580498866e-6),kmul(dy,dz));
  const CCTK_REAL_VEC p1o840dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00119047619047619047619047619048),dx);
  const CCTK_REAL_VEC p1o840dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00119047619047619047619047619048),dy);
  const CCTK_REAL_VEC p1o840dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00119047619047619047619047619048),dz);
  const CCTK_REAL_VEC p1odx2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1),kmul(dx,dx));
  const CCTK_REAL_VEC p1ody2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1),kmul(dy,dy));
  const CCTK_REAL_VEC p1odz2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1),kmul(dz,dz));
  const CCTK_REAL_VEC pm1o120dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.00833333333333333333333333333333),dx);
  const CCTK_REAL_VEC pm1o120dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.00833333333333333333333333333333),dy);
  const CCTK_REAL_VEC pm1o120dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.00833333333333333333333333333333),dz);
  const CCTK_REAL_VEC pm1o12dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.0833333333333333333333333333333),dx);
  const CCTK_REAL_VEC pm1o12dx2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.0833333333333333333333333333333),kmul(dx,dx));
  const CCTK_REAL_VEC pm1o12dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.0833333333333333333333333333333),dy);
  const CCTK_REAL_VEC pm1o12dy2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.0833333333333333333333333333333),kmul(dy,dy));
  const CCTK_REAL_VEC pm1o12dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.0833333333333333333333333333333),dz);
  const CCTK_REAL_VEC pm1o12dz2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.0833333333333333333333333333333),kmul(dz,dz));
  const CCTK_REAL_VEC pm1o16dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.0625),dx);
  const CCTK_REAL_VEC pm1o16dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.0625),dy);
  const CCTK_REAL_VEC pm1o16dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.0625),dz);
  const CCTK_REAL_VEC pm1o256dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.00390625),dx);
  const CCTK_REAL_VEC pm1o256dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.00390625),dy);
  const CCTK_REAL_VEC pm1o256dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.00390625),dz);
  const CCTK_REAL_VEC pm1o2dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.5),dx);
  const CCTK_REAL_VEC pm1o2dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.5),dy);
  const CCTK_REAL_VEC pm1o2dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.5),dz);
  const CCTK_REAL_VEC pm1o4dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.25),dx);
  const CCTK_REAL_VEC pm1o4dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.25),dy);
  const CCTK_REAL_VEC pm1o4dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.25),dz);
  const CCTK_REAL_VEC pm1o60dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.0166666666666666666666666666667),dx);
  const CCTK_REAL_VEC pm1o60dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.0166666666666666666666666666667),dy);
  const CCTK_REAL_VEC pm1o60dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.0166666666666666666666666666667),dz);
  /* Jacobian variable pointers */
  const bool usejacobian1 = (!CCTK_IsFunctionAliased("MultiPatch_GetMap") || MultiPatch_GetMap(cctkGH) != jacobian_identity_map)
                        && strlen(jacobian_group) > 0;
  const bool usejacobian = assume_use_jacobian>=0 ? assume_use_jacobian : usejacobian1;
  if (usejacobian && (strlen(jacobian_derivative_group) == 0))
  {
    CCTK_WARN(CCTK_WARN_ALERT, "GenericFD::jacobian_group and GenericFD::jacobian_derivative_group must both be set to valid group names");
  }
  
  const CCTK_REAL* restrict jacobian_ptrs[9];
  if (usejacobian) GroupDataPointers(cctkGH, jacobian_group,
                                                9, jacobian_ptrs);
  
  const CCTK_REAL* restrict const J11 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_ptrs[0] : 0;
  const CCTK_REAL* restrict const J12 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_ptrs[1] : 0;
  const CCTK_REAL* restrict const J13 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_ptrs[2] : 0;
  const CCTK_REAL* restrict const J21 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_ptrs[3] : 0;
  const CCTK_REAL* restrict const J22 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_ptrs[4] : 0;
  const CCTK_REAL* restrict const J23 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_ptrs[5] : 0;
  const CCTK_REAL* restrict const J31 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_ptrs[6] : 0;
  const CCTK_REAL* restrict const J32 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_ptrs[7] : 0;
  const CCTK_REAL* restrict const J33 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_ptrs[8] : 0;
  
  const CCTK_REAL* restrict jacobian_determinant_ptrs[1] CCTK_ATTRIBUTE_UNUSED;
  if (usejacobian && strlen(jacobian_determinant_group) > 0) GroupDataPointers(cctkGH, jacobian_determinant_group,
                                                1, jacobian_determinant_ptrs);
  
  const CCTK_REAL* restrict const detJ CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_determinant_ptrs[0] : 0;
  
  const CCTK_REAL* restrict jacobian_inverse_ptrs[9] CCTK_ATTRIBUTE_UNUSED;
  if (usejacobian && strlen(jacobian_inverse_group) > 0) GroupDataPointers(cctkGH, jacobian_inverse_group,
                                                9, jacobian_inverse_ptrs);
  
  const CCTK_REAL* restrict const iJ11 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_inverse_ptrs[0] : 0;
  const CCTK_REAL* restrict const iJ12 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_inverse_ptrs[1] : 0;
  const CCTK_REAL* restrict const iJ13 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_inverse_ptrs[2] : 0;
  const CCTK_REAL* restrict const iJ21 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_inverse_ptrs[3] : 0;
  const CCTK_REAL* restrict const iJ22 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_inverse_ptrs[4] : 0;
  const CCTK_REAL* restrict const iJ23 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_inverse_ptrs[5] : 0;
  const CCTK_REAL* restrict const iJ31 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_inverse_ptrs[6] : 0;
  const CCTK_REAL* restrict const iJ32 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_inverse_ptrs[7] : 0;
  const CCTK_REAL* restrict const iJ33 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_inverse_ptrs[8] : 0;
  
  const CCTK_REAL* restrict jacobian_derivative_ptrs[18] CCTK_ATTRIBUTE_UNUSED;
  if (usejacobian) GroupDataPointers(cctkGH, jacobian_derivative_group,
                                      18, jacobian_derivative_ptrs);
  
  const CCTK_REAL* restrict const dJ111 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[0] : 0;
  const CCTK_REAL* restrict const dJ112 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[1] : 0;
  const CCTK_REAL* restrict const dJ113 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[2] : 0;
  const CCTK_REAL* restrict const dJ122 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[3] : 0;
  const CCTK_REAL* restrict const dJ123 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[4] : 0;
  const CCTK_REAL* restrict const dJ133 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[5] : 0;
  const CCTK_REAL* restrict const dJ211 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[6] : 0;
  const CCTK_REAL* restrict const dJ212 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[7] : 0;
  const CCTK_REAL* restrict const dJ213 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[8] : 0;
  const CCTK_REAL* restrict const dJ222 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[9] : 0;
  const CCTK_REAL* restrict const dJ223 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[10] : 0;
  const CCTK_REAL* restrict const dJ233 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[11] : 0;
  const CCTK_REAL* restrict const dJ311 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[12] : 0;
  const CCTK_REAL* restrict const dJ312 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[13] : 0;
  const CCTK_REAL* restrict const dJ313 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[14] : 0;
  const CCTK_REAL* restrict const dJ322 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[15] : 0;
  const CCTK_REAL* restrict const dJ323 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[16] : 0;
  const CCTK_REAL* restrict const dJ333 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[17] : 0;
  /* Assign local copies of arrays functions */
  
  
  /* Calculate temporaries and arrays functions */
  /* Copy local copies back to grid functions */
  /* Loop over the grid points */
  const int imin0=imin[0];
  const int imin1=imin[1];
  const int imin2=imin[2];
  const int imax0=imax[0];
  const int imax1=imax[1];
  const int imax2=imax[2];
  #pragma omp parallel
  CCTK_LOOP3STR(ML_BSSN_InitialADMBase2Interior,
    i,j,k, imin0,imin1,imin2, imax0,imax1,imax2,
    cctk_ash[0],cctk_ash[1],cctk_ash[2],
    vecimin,vecimax, CCTK_REAL_VEC_SIZE)
  {
    const ptrdiff_t index CCTK_ATTRIBUTE_UNUSED = di*i + dj*j + dk*k;
    /* Assign local copies of grid functions */
    
    CCTK_REAL_VEC alpL CCTK_ATTRIBUTE_UNUSED = vec_load(alp[index]);
    CCTK_REAL_VEC betaxL CCTK_ATTRIBUTE_UNUSED = vec_load(betax[index]);
    CCTK_REAL_VEC betayL CCTK_ATTRIBUTE_UNUSED = vec_load(betay[index]);
    CCTK_REAL_VEC betazL CCTK_ATTRIBUTE_UNUSED = vec_load(betaz[index]);
    CCTK_REAL_VEC dtalpL CCTK_ATTRIBUTE_UNUSED = vec_load(dtalp[index]);
    CCTK_REAL_VEC dtbetaxL CCTK_ATTRIBUTE_UNUSED = 
      vec_load(dtbetax[index]);
    CCTK_REAL_VEC dtbetayL CCTK_ATTRIBUTE_UNUSED = 
      vec_load(dtbetay[index]);
    CCTK_REAL_VEC dtbetazL CCTK_ATTRIBUTE_UNUSED = 
      vec_load(dtbetaz[index]);
    CCTK_REAL_VEC gt11L CCTK_ATTRIBUTE_UNUSED = vec_load(gt11[index]);
    CCTK_REAL_VEC gt12L CCTK_ATTRIBUTE_UNUSED = vec_load(gt12[index]);
    CCTK_REAL_VEC gt13L CCTK_ATTRIBUTE_UNUSED = vec_load(gt13[index]);
    CCTK_REAL_VEC gt22L CCTK_ATTRIBUTE_UNUSED = vec_load(gt22[index]);
    CCTK_REAL_VEC gt23L CCTK_ATTRIBUTE_UNUSED = vec_load(gt23[index]);
    CCTK_REAL_VEC gt33L CCTK_ATTRIBUTE_UNUSED = vec_load(gt33[index]);
    CCTK_REAL_VEC gxxL CCTK_ATTRIBUTE_UNUSED = vec_load(gxx[index]);
    CCTK_REAL_VEC gxyL CCTK_ATTRIBUTE_UNUSED = vec_load(gxy[index]);
    CCTK_REAL_VEC gxzL CCTK_ATTRIBUTE_UNUSED = vec_load(gxz[index]);
    CCTK_REAL_VEC gyyL CCTK_ATTRIBUTE_UNUSED = vec_load(gyy[index]);
    CCTK_REAL_VEC gyzL CCTK_ATTRIBUTE_UNUSED = vec_load(gyz[index]);
    CCTK_REAL_VEC gzzL CCTK_ATTRIBUTE_UNUSED = vec_load(gzz[index]);
    CCTK_REAL_VEC phiL CCTK_ATTRIBUTE_UNUSED = vec_load(phi[index]);
    CCTK_REAL_VEC rL CCTK_ATTRIBUTE_UNUSED = vec_load(r[index]);
    
    
    CCTK_REAL_VEC J11L, J12L, J13L, J21L, J22L, J23L, J31L, J32L, J33L CCTK_ATTRIBUTE_UNUSED ;
    
    if (usejacobian)
    {
      J11L = vec_load(J11[index]);
      J12L = vec_load(J12[index]);
      J13L = vec_load(J13[index]);
      J21L = vec_load(J21[index]);
      J22L = vec_load(J22[index]);
      J23L = vec_load(J23[index]);
      J31L = vec_load(J31[index]);
      J32L = vec_load(J32[index]);
      J33L = vec_load(J33[index]);
    }
    /* Include user supplied include files */
    /* Precompute derivatives */
    CCTK_REAL_VEC PDupwindNthSymm1alp CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm2alp CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm3alp CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti1alp CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti2alp CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti3alp CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm1betax CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm2betax CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm3betax CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti1betax CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti2betax CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti3betax CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm1betay CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm2betay CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm3betay CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti1betay CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti2betay CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti3betay CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm1betaz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm2betaz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm3betaz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti1betaz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti2betaz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti3betaz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3gt33 CCTK_ATTRIBUTE_UNUSED;
    
    switch (fdOrder)
    {
      case 2:
      {
        PDupwindNthSymm1alp = PDupwindNthSymmfdOrder21(&alp[index]);
        PDupwindNthSymm2alp = PDupwindNthSymmfdOrder22(&alp[index]);
        PDupwindNthSymm3alp = PDupwindNthSymmfdOrder23(&alp[index]);
        PDupwindNthAnti1alp = PDupwindNthAntifdOrder21(&alp[index]);
        PDupwindNthAnti2alp = PDupwindNthAntifdOrder22(&alp[index]);
        PDupwindNthAnti3alp = PDupwindNthAntifdOrder23(&alp[index]);
        PDupwindNthSymm1betax = PDupwindNthSymmfdOrder21(&betax[index]);
        PDupwindNthSymm2betax = PDupwindNthSymmfdOrder22(&betax[index]);
        PDupwindNthSymm3betax = PDupwindNthSymmfdOrder23(&betax[index]);
        PDupwindNthAnti1betax = PDupwindNthAntifdOrder21(&betax[index]);
        PDupwindNthAnti2betax = PDupwindNthAntifdOrder22(&betax[index]);
        PDupwindNthAnti3betax = PDupwindNthAntifdOrder23(&betax[index]);
        PDupwindNthSymm1betay = PDupwindNthSymmfdOrder21(&betay[index]);
        PDupwindNthSymm2betay = PDupwindNthSymmfdOrder22(&betay[index]);
        PDupwindNthSymm3betay = PDupwindNthSymmfdOrder23(&betay[index]);
        PDupwindNthAnti1betay = PDupwindNthAntifdOrder21(&betay[index]);
        PDupwindNthAnti2betay = PDupwindNthAntifdOrder22(&betay[index]);
        PDupwindNthAnti3betay = PDupwindNthAntifdOrder23(&betay[index]);
        PDupwindNthSymm1betaz = PDupwindNthSymmfdOrder21(&betaz[index]);
        PDupwindNthSymm2betaz = PDupwindNthSymmfdOrder22(&betaz[index]);
        PDupwindNthSymm3betaz = PDupwindNthSymmfdOrder23(&betaz[index]);
        PDupwindNthAnti1betaz = PDupwindNthAntifdOrder21(&betaz[index]);
        PDupwindNthAnti2betaz = PDupwindNthAntifdOrder22(&betaz[index]);
        PDupwindNthAnti3betaz = PDupwindNthAntifdOrder23(&betaz[index]);
        PDstandardNth1gt11 = PDstandardNthfdOrder21(&gt11[index]);
        PDstandardNth2gt11 = PDstandardNthfdOrder22(&gt11[index]);
        PDstandardNth3gt11 = PDstandardNthfdOrder23(&gt11[index]);
        PDstandardNth1gt12 = PDstandardNthfdOrder21(&gt12[index]);
        PDstandardNth2gt12 = PDstandardNthfdOrder22(&gt12[index]);
        PDstandardNth3gt12 = PDstandardNthfdOrder23(&gt12[index]);
        PDstandardNth1gt13 = PDstandardNthfdOrder21(&gt13[index]);
        PDstandardNth2gt13 = PDstandardNthfdOrder22(&gt13[index]);
        PDstandardNth3gt13 = PDstandardNthfdOrder23(&gt13[index]);
        PDstandardNth1gt22 = PDstandardNthfdOrder21(&gt22[index]);
        PDstandardNth2gt22 = PDstandardNthfdOrder22(&gt22[index]);
        PDstandardNth3gt22 = PDstandardNthfdOrder23(&gt22[index]);
        PDstandardNth1gt23 = PDstandardNthfdOrder21(&gt23[index]);
        PDstandardNth2gt23 = PDstandardNthfdOrder22(&gt23[index]);
        PDstandardNth3gt23 = PDstandardNthfdOrder23(&gt23[index]);
        PDstandardNth1gt33 = PDstandardNthfdOrder21(&gt33[index]);
        PDstandardNth2gt33 = PDstandardNthfdOrder22(&gt33[index]);
        PDstandardNth3gt33 = PDstandardNthfdOrder23(&gt33[index]);
        break;
      }
      
      case 4:
      {
        PDupwindNthSymm1alp = PDupwindNthSymmfdOrder41(&alp[index]);
        PDupwindNthSymm2alp = PDupwindNthSymmfdOrder42(&alp[index]);
        PDupwindNthSymm3alp = PDupwindNthSymmfdOrder43(&alp[index]);
        PDupwindNthAnti1alp = PDupwindNthAntifdOrder41(&alp[index]);
        PDupwindNthAnti2alp = PDupwindNthAntifdOrder42(&alp[index]);
        PDupwindNthAnti3alp = PDupwindNthAntifdOrder43(&alp[index]);
        PDupwindNthSymm1betax = PDupwindNthSymmfdOrder41(&betax[index]);
        PDupwindNthSymm2betax = PDupwindNthSymmfdOrder42(&betax[index]);
        PDupwindNthSymm3betax = PDupwindNthSymmfdOrder43(&betax[index]);
        PDupwindNthAnti1betax = PDupwindNthAntifdOrder41(&betax[index]);
        PDupwindNthAnti2betax = PDupwindNthAntifdOrder42(&betax[index]);
        PDupwindNthAnti3betax = PDupwindNthAntifdOrder43(&betax[index]);
        PDupwindNthSymm1betay = PDupwindNthSymmfdOrder41(&betay[index]);
        PDupwindNthSymm2betay = PDupwindNthSymmfdOrder42(&betay[index]);
        PDupwindNthSymm3betay = PDupwindNthSymmfdOrder43(&betay[index]);
        PDupwindNthAnti1betay = PDupwindNthAntifdOrder41(&betay[index]);
        PDupwindNthAnti2betay = PDupwindNthAntifdOrder42(&betay[index]);
        PDupwindNthAnti3betay = PDupwindNthAntifdOrder43(&betay[index]);
        PDupwindNthSymm1betaz = PDupwindNthSymmfdOrder41(&betaz[index]);
        PDupwindNthSymm2betaz = PDupwindNthSymmfdOrder42(&betaz[index]);
        PDupwindNthSymm3betaz = PDupwindNthSymmfdOrder43(&betaz[index]);
        PDupwindNthAnti1betaz = PDupwindNthAntifdOrder41(&betaz[index]);
        PDupwindNthAnti2betaz = PDupwindNthAntifdOrder42(&betaz[index]);
        PDupwindNthAnti3betaz = PDupwindNthAntifdOrder43(&betaz[index]);
        PDstandardNth1gt11 = PDstandardNthfdOrder41(&gt11[index]);
        PDstandardNth2gt11 = PDstandardNthfdOrder42(&gt11[index]);
        PDstandardNth3gt11 = PDstandardNthfdOrder43(&gt11[index]);
        PDstandardNth1gt12 = PDstandardNthfdOrder41(&gt12[index]);
        PDstandardNth2gt12 = PDstandardNthfdOrder42(&gt12[index]);
        PDstandardNth3gt12 = PDstandardNthfdOrder43(&gt12[index]);
        PDstandardNth1gt13 = PDstandardNthfdOrder41(&gt13[index]);
        PDstandardNth2gt13 = PDstandardNthfdOrder42(&gt13[index]);
        PDstandardNth3gt13 = PDstandardNthfdOrder43(&gt13[index]);
        PDstandardNth1gt22 = PDstandardNthfdOrder41(&gt22[index]);
        PDstandardNth2gt22 = PDstandardNthfdOrder42(&gt22[index]);
        PDstandardNth3gt22 = PDstandardNthfdOrder43(&gt22[index]);
        PDstandardNth1gt23 = PDstandardNthfdOrder41(&gt23[index]);
        PDstandardNth2gt23 = PDstandardNthfdOrder42(&gt23[index]);
        PDstandardNth3gt23 = PDstandardNthfdOrder43(&gt23[index]);
        PDstandardNth1gt33 = PDstandardNthfdOrder41(&gt33[index]);
        PDstandardNth2gt33 = PDstandardNthfdOrder42(&gt33[index]);
        PDstandardNth3gt33 = PDstandardNthfdOrder43(&gt33[index]);
        break;
      }
      
      case 6:
      {
        PDupwindNthSymm1alp = PDupwindNthSymmfdOrder61(&alp[index]);
        PDupwindNthSymm2alp = PDupwindNthSymmfdOrder62(&alp[index]);
        PDupwindNthSymm3alp = PDupwindNthSymmfdOrder63(&alp[index]);
        PDupwindNthAnti1alp = PDupwindNthAntifdOrder61(&alp[index]);
        PDupwindNthAnti2alp = PDupwindNthAntifdOrder62(&alp[index]);
        PDupwindNthAnti3alp = PDupwindNthAntifdOrder63(&alp[index]);
        PDupwindNthSymm1betax = PDupwindNthSymmfdOrder61(&betax[index]);
        PDupwindNthSymm2betax = PDupwindNthSymmfdOrder62(&betax[index]);
        PDupwindNthSymm3betax = PDupwindNthSymmfdOrder63(&betax[index]);
        PDupwindNthAnti1betax = PDupwindNthAntifdOrder61(&betax[index]);
        PDupwindNthAnti2betax = PDupwindNthAntifdOrder62(&betax[index]);
        PDupwindNthAnti3betax = PDupwindNthAntifdOrder63(&betax[index]);
        PDupwindNthSymm1betay = PDupwindNthSymmfdOrder61(&betay[index]);
        PDupwindNthSymm2betay = PDupwindNthSymmfdOrder62(&betay[index]);
        PDupwindNthSymm3betay = PDupwindNthSymmfdOrder63(&betay[index]);
        PDupwindNthAnti1betay = PDupwindNthAntifdOrder61(&betay[index]);
        PDupwindNthAnti2betay = PDupwindNthAntifdOrder62(&betay[index]);
        PDupwindNthAnti3betay = PDupwindNthAntifdOrder63(&betay[index]);
        PDupwindNthSymm1betaz = PDupwindNthSymmfdOrder61(&betaz[index]);
        PDupwindNthSymm2betaz = PDupwindNthSymmfdOrder62(&betaz[index]);
        PDupwindNthSymm3betaz = PDupwindNthSymmfdOrder63(&betaz[index]);
        PDupwindNthAnti1betaz = PDupwindNthAntifdOrder61(&betaz[index]);
        PDupwindNthAnti2betaz = PDupwindNthAntifdOrder62(&betaz[index]);
        PDupwindNthAnti3betaz = PDupwindNthAntifdOrder63(&betaz[index]);
        PDstandardNth1gt11 = PDstandardNthfdOrder61(&gt11[index]);
        PDstandardNth2gt11 = PDstandardNthfdOrder62(&gt11[index]);
        PDstandardNth3gt11 = PDstandardNthfdOrder63(&gt11[index]);
        PDstandardNth1gt12 = PDstandardNthfdOrder61(&gt12[index]);
        PDstandardNth2gt12 = PDstandardNthfdOrder62(&gt12[index]);
        PDstandardNth3gt12 = PDstandardNthfdOrder63(&gt12[index]);
        PDstandardNth1gt13 = PDstandardNthfdOrder61(&gt13[index]);
        PDstandardNth2gt13 = PDstandardNthfdOrder62(&gt13[index]);
        PDstandardNth3gt13 = PDstandardNthfdOrder63(&gt13[index]);
        PDstandardNth1gt22 = PDstandardNthfdOrder61(&gt22[index]);
        PDstandardNth2gt22 = PDstandardNthfdOrder62(&gt22[index]);
        PDstandardNth3gt22 = PDstandardNthfdOrder63(&gt22[index]);
        PDstandardNth1gt23 = PDstandardNthfdOrder61(&gt23[index]);
        PDstandardNth2gt23 = PDstandardNthfdOrder62(&gt23[index]);
        PDstandardNth3gt23 = PDstandardNthfdOrder63(&gt23[index]);
        PDstandardNth1gt33 = PDstandardNthfdOrder61(&gt33[index]);
        PDstandardNth2gt33 = PDstandardNthfdOrder62(&gt33[index]);
        PDstandardNth3gt33 = PDstandardNthfdOrder63(&gt33[index]);
        break;
      }
      
      case 8:
      {
        PDupwindNthSymm1alp = PDupwindNthSymmfdOrder81(&alp[index]);
        PDupwindNthSymm2alp = PDupwindNthSymmfdOrder82(&alp[index]);
        PDupwindNthSymm3alp = PDupwindNthSymmfdOrder83(&alp[index]);
        PDupwindNthAnti1alp = PDupwindNthAntifdOrder81(&alp[index]);
        PDupwindNthAnti2alp = PDupwindNthAntifdOrder82(&alp[index]);
        PDupwindNthAnti3alp = PDupwindNthAntifdOrder83(&alp[index]);
        PDupwindNthSymm1betax = PDupwindNthSymmfdOrder81(&betax[index]);
        PDupwindNthSymm2betax = PDupwindNthSymmfdOrder82(&betax[index]);
        PDupwindNthSymm3betax = PDupwindNthSymmfdOrder83(&betax[index]);
        PDupwindNthAnti1betax = PDupwindNthAntifdOrder81(&betax[index]);
        PDupwindNthAnti2betax = PDupwindNthAntifdOrder82(&betax[index]);
        PDupwindNthAnti3betax = PDupwindNthAntifdOrder83(&betax[index]);
        PDupwindNthSymm1betay = PDupwindNthSymmfdOrder81(&betay[index]);
        PDupwindNthSymm2betay = PDupwindNthSymmfdOrder82(&betay[index]);
        PDupwindNthSymm3betay = PDupwindNthSymmfdOrder83(&betay[index]);
        PDupwindNthAnti1betay = PDupwindNthAntifdOrder81(&betay[index]);
        PDupwindNthAnti2betay = PDupwindNthAntifdOrder82(&betay[index]);
        PDupwindNthAnti3betay = PDupwindNthAntifdOrder83(&betay[index]);
        PDupwindNthSymm1betaz = PDupwindNthSymmfdOrder81(&betaz[index]);
        PDupwindNthSymm2betaz = PDupwindNthSymmfdOrder82(&betaz[index]);
        PDupwindNthSymm3betaz = PDupwindNthSymmfdOrder83(&betaz[index]);
        PDupwindNthAnti1betaz = PDupwindNthAntifdOrder81(&betaz[index]);
        PDupwindNthAnti2betaz = PDupwindNthAntifdOrder82(&betaz[index]);
        PDupwindNthAnti3betaz = PDupwindNthAntifdOrder83(&betaz[index]);
        PDstandardNth1gt11 = PDstandardNthfdOrder81(&gt11[index]);
        PDstandardNth2gt11 = PDstandardNthfdOrder82(&gt11[index]);
        PDstandardNth3gt11 = PDstandardNthfdOrder83(&gt11[index]);
        PDstandardNth1gt12 = PDstandardNthfdOrder81(&gt12[index]);
        PDstandardNth2gt12 = PDstandardNthfdOrder82(&gt12[index]);
        PDstandardNth3gt12 = PDstandardNthfdOrder83(&gt12[index]);
        PDstandardNth1gt13 = PDstandardNthfdOrder81(&gt13[index]);
        PDstandardNth2gt13 = PDstandardNthfdOrder82(&gt13[index]);
        PDstandardNth3gt13 = PDstandardNthfdOrder83(&gt13[index]);
        PDstandardNth1gt22 = PDstandardNthfdOrder81(&gt22[index]);
        PDstandardNth2gt22 = PDstandardNthfdOrder82(&gt22[index]);
        PDstandardNth3gt22 = PDstandardNthfdOrder83(&gt22[index]);
        PDstandardNth1gt23 = PDstandardNthfdOrder81(&gt23[index]);
        PDstandardNth2gt23 = PDstandardNthfdOrder82(&gt23[index]);
        PDstandardNth3gt23 = PDstandardNthfdOrder83(&gt23[index]);
        PDstandardNth1gt33 = PDstandardNthfdOrder81(&gt33[index]);
        PDstandardNth2gt33 = PDstandardNthfdOrder82(&gt33[index]);
        PDstandardNth3gt33 = PDstandardNthfdOrder83(&gt33[index]);
        break;
      }
      default:
        CCTK_BUILTIN_UNREACHABLE();
    }
    /* Calculate temporaries and grid functions */
    CCTK_REAL_VEC JacPDstandardNth1gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti1alp CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti1betax CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti1betay CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti1betaz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti2alp CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti2betax CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti2betay CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti2betaz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti3alp CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti3betax CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti3betay CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti3betaz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm1alp CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm1betax CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm1betay CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm1betaz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm2alp CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm2betax CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm2betay CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm2betaz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm3alp CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm3betax CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm3betay CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm3betaz CCTK_ATTRIBUTE_UNUSED;
    
    if (usejacobian)
    {
      JacPDstandardNth1gt11 = 
        kmadd(J11L,PDstandardNth1gt11,kmadd(J21L,PDstandardNth2gt11,kmul(J31L,PDstandardNth3gt11)));
      
      JacPDstandardNth1gt12 = 
        kmadd(J11L,PDstandardNth1gt12,kmadd(J21L,PDstandardNth2gt12,kmul(J31L,PDstandardNth3gt12)));
      
      JacPDstandardNth1gt13 = 
        kmadd(J11L,PDstandardNth1gt13,kmadd(J21L,PDstandardNth2gt13,kmul(J31L,PDstandardNth3gt13)));
      
      JacPDstandardNth1gt22 = 
        kmadd(J11L,PDstandardNth1gt22,kmadd(J21L,PDstandardNth2gt22,kmul(J31L,PDstandardNth3gt22)));
      
      JacPDstandardNth1gt23 = 
        kmadd(J11L,PDstandardNth1gt23,kmadd(J21L,PDstandardNth2gt23,kmul(J31L,PDstandardNth3gt23)));
      
      JacPDstandardNth1gt33 = 
        kmadd(J11L,PDstandardNth1gt33,kmadd(J21L,PDstandardNth2gt33,kmul(J31L,PDstandardNth3gt33)));
      
      JacPDstandardNth2gt11 = 
        kmadd(J12L,PDstandardNth1gt11,kmadd(J22L,PDstandardNth2gt11,kmul(J32L,PDstandardNth3gt11)));
      
      JacPDstandardNth2gt12 = 
        kmadd(J12L,PDstandardNth1gt12,kmadd(J22L,PDstandardNth2gt12,kmul(J32L,PDstandardNth3gt12)));
      
      JacPDstandardNth2gt13 = 
        kmadd(J12L,PDstandardNth1gt13,kmadd(J22L,PDstandardNth2gt13,kmul(J32L,PDstandardNth3gt13)));
      
      JacPDstandardNth2gt22 = 
        kmadd(J12L,PDstandardNth1gt22,kmadd(J22L,PDstandardNth2gt22,kmul(J32L,PDstandardNth3gt22)));
      
      JacPDstandardNth2gt23 = 
        kmadd(J12L,PDstandardNth1gt23,kmadd(J22L,PDstandardNth2gt23,kmul(J32L,PDstandardNth3gt23)));
      
      JacPDstandardNth2gt33 = 
        kmadd(J12L,PDstandardNth1gt33,kmadd(J22L,PDstandardNth2gt33,kmul(J32L,PDstandardNth3gt33)));
      
      JacPDstandardNth3gt11 = 
        kmadd(J13L,PDstandardNth1gt11,kmadd(J23L,PDstandardNth2gt11,kmul(J33L,PDstandardNth3gt11)));
      
      JacPDstandardNth3gt12 = 
        kmadd(J13L,PDstandardNth1gt12,kmadd(J23L,PDstandardNth2gt12,kmul(J33L,PDstandardNth3gt12)));
      
      JacPDstandardNth3gt13 = 
        kmadd(J13L,PDstandardNth1gt13,kmadd(J23L,PDstandardNth2gt13,kmul(J33L,PDstandardNth3gt13)));
      
      JacPDstandardNth3gt22 = 
        kmadd(J13L,PDstandardNth1gt22,kmadd(J23L,PDstandardNth2gt22,kmul(J33L,PDstandardNth3gt22)));
      
      JacPDstandardNth3gt23 = 
        kmadd(J13L,PDstandardNth1gt23,kmadd(J23L,PDstandardNth2gt23,kmul(J33L,PDstandardNth3gt23)));
      
      JacPDstandardNth3gt33 = 
        kmadd(J13L,PDstandardNth1gt33,kmadd(J23L,PDstandardNth2gt33,kmul(J33L,PDstandardNth3gt33)));
      
      JacPDupwindNthSymm1alp = 
        kmadd(J11L,PDupwindNthSymm1alp,kmadd(J21L,PDupwindNthSymm2alp,kmul(J31L,PDupwindNthSymm3alp)));
      
      JacPDupwindNthSymm1betax = 
        kmadd(J11L,PDupwindNthSymm1betax,kmadd(J21L,PDupwindNthSymm2betax,kmul(J31L,PDupwindNthSymm3betax)));
      
      JacPDupwindNthSymm1betay = 
        kmadd(J11L,PDupwindNthSymm1betay,kmadd(J21L,PDupwindNthSymm2betay,kmul(J31L,PDupwindNthSymm3betay)));
      
      JacPDupwindNthSymm1betaz = 
        kmadd(J11L,PDupwindNthSymm1betaz,kmadd(J21L,PDupwindNthSymm2betaz,kmul(J31L,PDupwindNthSymm3betaz)));
      
      JacPDupwindNthSymm2alp = 
        kmadd(J12L,PDupwindNthSymm1alp,kmadd(J22L,PDupwindNthSymm2alp,kmul(J32L,PDupwindNthSymm3alp)));
      
      JacPDupwindNthSymm2betax = 
        kmadd(J12L,PDupwindNthSymm1betax,kmadd(J22L,PDupwindNthSymm2betax,kmul(J32L,PDupwindNthSymm3betax)));
      
      JacPDupwindNthSymm2betay = 
        kmadd(J12L,PDupwindNthSymm1betay,kmadd(J22L,PDupwindNthSymm2betay,kmul(J32L,PDupwindNthSymm3betay)));
      
      JacPDupwindNthSymm2betaz = 
        kmadd(J12L,PDupwindNthSymm1betaz,kmadd(J22L,PDupwindNthSymm2betaz,kmul(J32L,PDupwindNthSymm3betaz)));
      
      JacPDupwindNthSymm3alp = 
        kmadd(J13L,PDupwindNthSymm1alp,kmadd(J23L,PDupwindNthSymm2alp,kmul(J33L,PDupwindNthSymm3alp)));
      
      JacPDupwindNthSymm3betax = 
        kmadd(J13L,PDupwindNthSymm1betax,kmadd(J23L,PDupwindNthSymm2betax,kmul(J33L,PDupwindNthSymm3betax)));
      
      JacPDupwindNthSymm3betay = 
        kmadd(J13L,PDupwindNthSymm1betay,kmadd(J23L,PDupwindNthSymm2betay,kmul(J33L,PDupwindNthSymm3betay)));
      
      JacPDupwindNthSymm3betaz = 
        kmadd(J13L,PDupwindNthSymm1betaz,kmadd(J23L,PDupwindNthSymm2betaz,kmul(J33L,PDupwindNthSymm3betaz)));
      
      JacPDupwindNthAnti1alp = 
        kmadd(J11L,PDupwindNthAnti1alp,kmadd(J21L,PDupwindNthAnti2alp,kmul(J31L,PDupwindNthAnti3alp)));
      
      JacPDupwindNthAnti1betax = 
        kmadd(J11L,PDupwindNthAnti1betax,kmadd(J21L,PDupwindNthAnti2betax,kmul(J31L,PDupwindNthAnti3betax)));
      
      JacPDupwindNthAnti1betay = 
        kmadd(J11L,PDupwindNthAnti1betay,kmadd(J21L,PDupwindNthAnti2betay,kmul(J31L,PDupwindNthAnti3betay)));
      
      JacPDupwindNthAnti1betaz = 
        kmadd(J11L,PDupwindNthAnti1betaz,kmadd(J21L,PDupwindNthAnti2betaz,kmul(J31L,PDupwindNthAnti3betaz)));
      
      JacPDupwindNthAnti2alp = 
        kmadd(J12L,PDupwindNthAnti1alp,kmadd(J22L,PDupwindNthAnti2alp,kmul(J32L,PDupwindNthAnti3alp)));
      
      JacPDupwindNthAnti2betax = 
        kmadd(J12L,PDupwindNthAnti1betax,kmadd(J22L,PDupwindNthAnti2betax,kmul(J32L,PDupwindNthAnti3betax)));
      
      JacPDupwindNthAnti2betay = 
        kmadd(J12L,PDupwindNthAnti1betay,kmadd(J22L,PDupwindNthAnti2betay,kmul(J32L,PDupwindNthAnti3betay)));
      
      JacPDupwindNthAnti2betaz = 
        kmadd(J12L,PDupwindNthAnti1betaz,kmadd(J22L,PDupwindNthAnti2betaz,kmul(J32L,PDupwindNthAnti3betaz)));
      
      JacPDupwindNthAnti3alp = 
        kmadd(J13L,PDupwindNthAnti1alp,kmadd(J23L,PDupwindNthAnti2alp,kmul(J33L,PDupwindNthAnti3alp)));
      
      JacPDupwindNthAnti3betax = 
        kmadd(J13L,PDupwindNthAnti1betax,kmadd(J23L,PDupwindNthAnti2betax,kmul(J33L,PDupwindNthAnti3betax)));
      
      JacPDupwindNthAnti3betay = 
        kmadd(J13L,PDupwindNthAnti1betay,kmadd(J23L,PDupwindNthAnti2betay,kmul(J33L,PDupwindNthAnti3betay)));
      
      JacPDupwindNthAnti3betaz = 
        kmadd(J13L,PDupwindNthAnti1betaz,kmadd(J23L,PDupwindNthAnti2betaz,kmul(J33L,PDupwindNthAnti3betaz)));
    }
    else
    {
      JacPDstandardNth1gt11 = PDstandardNth1gt11;
      
      JacPDstandardNth1gt12 = PDstandardNth1gt12;
      
      JacPDstandardNth1gt13 = PDstandardNth1gt13;
      
      JacPDstandardNth1gt22 = PDstandardNth1gt22;
      
      JacPDstandardNth1gt23 = PDstandardNth1gt23;
      
      JacPDstandardNth1gt33 = PDstandardNth1gt33;
      
      JacPDstandardNth2gt11 = PDstandardNth2gt11;
      
      JacPDstandardNth2gt12 = PDstandardNth2gt12;
      
      JacPDstandardNth2gt13 = PDstandardNth2gt13;
      
      JacPDstandardNth2gt22 = PDstandardNth2gt22;
      
      JacPDstandardNth2gt23 = PDstandardNth2gt23;
      
      JacPDstandardNth2gt33 = PDstandardNth2gt33;
      
      JacPDstandardNth3gt11 = PDstandardNth3gt11;
      
      JacPDstandardNth3gt12 = PDstandardNth3gt12;
      
      JacPDstandardNth3gt13 = PDstandardNth3gt13;
      
      JacPDstandardNth3gt22 = PDstandardNth3gt22;
      
      JacPDstandardNth3gt23 = PDstandardNth3gt23;
      
      JacPDstandardNth3gt33 = PDstandardNth3gt33;
      
      JacPDupwindNthSymm1alp = PDupwindNthSymm1alp;
      
      JacPDupwindNthSymm1betax = PDupwindNthSymm1betax;
      
      JacPDupwindNthSymm1betay = PDupwindNthSymm1betay;
      
      JacPDupwindNthSymm1betaz = PDupwindNthSymm1betaz;
      
      JacPDupwindNthSymm2alp = PDupwindNthSymm2alp;
      
      JacPDupwindNthSymm2betax = PDupwindNthSymm2betax;
      
      JacPDupwindNthSymm2betay = PDupwindNthSymm2betay;
      
      JacPDupwindNthSymm2betaz = PDupwindNthSymm2betaz;
      
      JacPDupwindNthSymm3alp = PDupwindNthSymm3alp;
      
      JacPDupwindNthSymm3betax = PDupwindNthSymm3betax;
      
      JacPDupwindNthSymm3betay = PDupwindNthSymm3betay;
      
      JacPDupwindNthSymm3betaz = PDupwindNthSymm3betaz;
      
      JacPDupwindNthAnti1alp = PDupwindNthAnti1alp;
      
      JacPDupwindNthAnti1betax = PDupwindNthAnti1betax;
      
      JacPDupwindNthAnti1betay = PDupwindNthAnti1betay;
      
      JacPDupwindNthAnti1betaz = PDupwindNthAnti1betaz;
      
      JacPDupwindNthAnti2alp = PDupwindNthAnti2alp;
      
      JacPDupwindNthAnti2betax = PDupwindNthAnti2betax;
      
      JacPDupwindNthAnti2betay = PDupwindNthAnti2betay;
      
      JacPDupwindNthAnti2betaz = PDupwindNthAnti2betaz;
      
      JacPDupwindNthAnti3alp = PDupwindNthAnti3alp;
      
      JacPDupwindNthAnti3betax = PDupwindNthAnti3betax;
      
      JacPDupwindNthAnti3betay = PDupwindNthAnti3betay;
      
      JacPDupwindNthAnti3betaz = PDupwindNthAnti3betaz;
    }
    
    CCTK_REAL_VEC g11 CCTK_ATTRIBUTE_UNUSED = gxxL;
    
    CCTK_REAL_VEC g12 CCTK_ATTRIBUTE_UNUSED = gxyL;
    
    CCTK_REAL_VEC g13 CCTK_ATTRIBUTE_UNUSED = gxzL;
    
    CCTK_REAL_VEC g22 CCTK_ATTRIBUTE_UNUSED = gyyL;
    
    CCTK_REAL_VEC g23 CCTK_ATTRIBUTE_UNUSED = gyzL;
    
    CCTK_REAL_VEC g33 CCTK_ATTRIBUTE_UNUSED = gzzL;
    
    CCTK_REAL_VEC detg CCTK_ATTRIBUTE_UNUSED = 
      kmadd(ToReal(2),kmul(g12,kmul(g13,g23)),knmsub(g33,kmul(g12,g12),kmsub(g22,kmsub(g11,g33,kmul(g13,g13)),kmul(g11,kmul(g23,g23)))));
    
    CCTK_REAL_VEC gu11 CCTK_ATTRIBUTE_UNUSED = 
      kdiv(kmsub(g22,g33,kmul(g23,g23)),detg);
    
    CCTK_REAL_VEC gu12 CCTK_ATTRIBUTE_UNUSED = 
      kdiv(kmsub(g13,g23,kmul(g12,g33)),detg);
    
    CCTK_REAL_VEC gu13 CCTK_ATTRIBUTE_UNUSED = 
      kdiv(kmsub(g12,g23,kmul(g13,g22)),detg);
    
    CCTK_REAL_VEC gu22 CCTK_ATTRIBUTE_UNUSED = 
      kdiv(kmsub(g11,g33,kmul(g13,g13)),detg);
    
    CCTK_REAL_VEC gu23 CCTK_ATTRIBUTE_UNUSED = 
      kdiv(kmsub(g12,g13,kmul(g11,g23)),detg);
    
    CCTK_REAL_VEC gu33 CCTK_ATTRIBUTE_UNUSED = 
      kdiv(kmsub(g11,g22,kmul(g12,g12)),detg);
    
    CCTK_REAL_VEC em4phi CCTK_ATTRIBUTE_UNUSED = IfThen(conformalMethod != 
      0,kmul(phiL,phiL),kexp(kmul(phiL,ToReal(-4))));
    
    CCTK_REAL_VEC gtu11 CCTK_ATTRIBUTE_UNUSED = kdiv(gu11,em4phi);
    
    CCTK_REAL_VEC gtu12 CCTK_ATTRIBUTE_UNUSED = kdiv(gu12,em4phi);
    
    CCTK_REAL_VEC gtu13 CCTK_ATTRIBUTE_UNUSED = kdiv(gu13,em4phi);
    
    CCTK_REAL_VEC gtu22 CCTK_ATTRIBUTE_UNUSED = kdiv(gu22,em4phi);
    
    CCTK_REAL_VEC gtu23 CCTK_ATTRIBUTE_UNUSED = kdiv(gu23,em4phi);
    
    CCTK_REAL_VEC gtu33 CCTK_ATTRIBUTE_UNUSED = kdiv(gu33,em4phi);
    
    CCTK_REAL_VEC Gt111 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gtu11,JacPDstandardNth1gt11,kmsub(ToReal(2),kmadd(gtu12,JacPDstandardNth1gt12,kmul(gtu13,JacPDstandardNth1gt13)),kmadd(gtu13,JacPDstandardNth3gt11,kmul(gtu12,JacPDstandardNth2gt11)))),ToReal(0.5));
    
    CCTK_REAL_VEC Gt211 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gtu12,JacPDstandardNth1gt11,kmsub(ToReal(2),kmadd(gtu22,JacPDstandardNth1gt12,kmul(gtu23,JacPDstandardNth1gt13)),kmadd(gtu23,JacPDstandardNth3gt11,kmul(gtu22,JacPDstandardNth2gt11)))),ToReal(0.5));
    
    CCTK_REAL_VEC Gt311 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gtu13,JacPDstandardNth1gt11,kmsub(ToReal(2),kmadd(gtu23,JacPDstandardNth1gt12,kmul(gtu33,JacPDstandardNth1gt13)),kmadd(gtu33,JacPDstandardNth3gt11,kmul(gtu23,JacPDstandardNth2gt11)))),ToReal(0.5));
    
    CCTK_REAL_VEC Gt112 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gtu12,JacPDstandardNth1gt22,kmadd(gtu11,JacPDstandardNth2gt11,kmul(gtu13,kadd(JacPDstandardNth1gt23,ksub(JacPDstandardNth2gt13,JacPDstandardNth3gt12))))),ToReal(0.5));
    
    CCTK_REAL_VEC Gt212 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gtu22,JacPDstandardNth1gt22,kmadd(gtu12,JacPDstandardNth2gt11,kmul(gtu23,kadd(JacPDstandardNth1gt23,ksub(JacPDstandardNth2gt13,JacPDstandardNth3gt12))))),ToReal(0.5));
    
    CCTK_REAL_VEC Gt312 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gtu23,JacPDstandardNth1gt22,kmadd(gtu13,JacPDstandardNth2gt11,kmul(gtu33,kadd(JacPDstandardNth1gt23,ksub(JacPDstandardNth2gt13,JacPDstandardNth3gt12))))),ToReal(0.5));
    
    CCTK_REAL_VEC Gt113 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gtu13,JacPDstandardNth1gt33,kmadd(gtu11,JacPDstandardNth3gt11,kmul(gtu12,kadd(JacPDstandardNth1gt23,ksub(JacPDstandardNth3gt12,JacPDstandardNth2gt13))))),ToReal(0.5));
    
    CCTK_REAL_VEC Gt213 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gtu23,JacPDstandardNth1gt33,kmadd(gtu12,JacPDstandardNth3gt11,kmul(gtu22,kadd(JacPDstandardNth1gt23,ksub(JacPDstandardNth3gt12,JacPDstandardNth2gt13))))),ToReal(0.5));
    
    CCTK_REAL_VEC Gt313 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gtu33,JacPDstandardNth1gt33,kmadd(gtu13,JacPDstandardNth3gt11,kmul(gtu23,kadd(JacPDstandardNth1gt23,ksub(JacPDstandardNth3gt12,JacPDstandardNth2gt13))))),ToReal(0.5));
    
    CCTK_REAL_VEC Gt122 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gtu11,kmsub(ToReal(2),JacPDstandardNth2gt12,JacPDstandardNth1gt22),kmadd(gtu12,JacPDstandardNth2gt22,kmul(gtu13,kmsub(ToReal(2),JacPDstandardNth2gt23,JacPDstandardNth3gt22)))),ToReal(0.5));
    
    CCTK_REAL_VEC Gt222 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gtu12,kmsub(ToReal(2),JacPDstandardNth2gt12,JacPDstandardNth1gt22),kmadd(gtu22,JacPDstandardNth2gt22,kmul(gtu23,kmsub(ToReal(2),JacPDstandardNth2gt23,JacPDstandardNth3gt22)))),ToReal(0.5));
    
    CCTK_REAL_VEC Gt322 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gtu13,kmsub(ToReal(2),JacPDstandardNth2gt12,JacPDstandardNth1gt22),kmadd(gtu23,JacPDstandardNth2gt22,kmul(gtu33,kmsub(ToReal(2),JacPDstandardNth2gt23,JacPDstandardNth3gt22)))),ToReal(0.5));
    
    CCTK_REAL_VEC Gt123 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gtu13,JacPDstandardNth2gt33,kmadd(gtu11,ksub(kadd(JacPDstandardNth2gt13,JacPDstandardNth3gt12),JacPDstandardNth1gt23),kmul(gtu12,JacPDstandardNth3gt22))),ToReal(0.5));
    
    CCTK_REAL_VEC Gt223 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gtu23,JacPDstandardNth2gt33,kmadd(gtu12,ksub(kadd(JacPDstandardNth2gt13,JacPDstandardNth3gt12),JacPDstandardNth1gt23),kmul(gtu22,JacPDstandardNth3gt22))),ToReal(0.5));
    
    CCTK_REAL_VEC Gt323 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gtu33,JacPDstandardNth2gt33,kmadd(gtu13,ksub(kadd(JacPDstandardNth2gt13,JacPDstandardNth3gt12),JacPDstandardNth1gt23),kmul(gtu23,JacPDstandardNth3gt22))),ToReal(0.5));
    
    CCTK_REAL_VEC Gt133 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gtu11,kmsub(ToReal(2),JacPDstandardNth3gt13,JacPDstandardNth1gt33),kmadd(gtu12,kmsub(ToReal(2),JacPDstandardNth3gt23,JacPDstandardNth2gt33),kmul(gtu13,JacPDstandardNth3gt33))),ToReal(0.5));
    
    CCTK_REAL_VEC Gt233 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gtu12,kmsub(ToReal(2),JacPDstandardNth3gt13,JacPDstandardNth1gt33),kmadd(gtu22,kmsub(ToReal(2),JacPDstandardNth3gt23,JacPDstandardNth2gt33),kmul(gtu23,JacPDstandardNth3gt33))),ToReal(0.5));
    
    CCTK_REAL_VEC Gt333 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gtu13,kmsub(ToReal(2),JacPDstandardNth3gt13,JacPDstandardNth1gt33),kmadd(gtu23,kmsub(ToReal(2),JacPDstandardNth3gt23,JacPDstandardNth2gt33),kmul(gtu33,JacPDstandardNth3gt33))),ToReal(0.5));
    
    CCTK_REAL_VEC Xt1L CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gt111,gtu11,kmadd(Gt122,gtu22,kmadd(ToReal(2),kmadd(Gt112,gtu12,kmadd(Gt113,gtu13,kmul(Gt123,gtu23))),kmul(Gt133,gtu33))));
    
    CCTK_REAL_VEC Xt2L CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gt211,gtu11,kmadd(Gt222,gtu22,kmadd(ToReal(2),kmadd(Gt212,gtu12,kmadd(Gt213,gtu13,kmul(Gt223,gtu23))),kmul(Gt233,gtu33))));
    
    CCTK_REAL_VEC Xt3L CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gt311,gtu11,kmadd(Gt322,gtu22,kmadd(ToReal(2),kmadd(Gt312,gtu12,kmadd(Gt313,gtu13,kmul(Gt323,gtu23))),kmul(Gt333,gtu33))));
    
    CCTK_REAL_VEC AL CCTK_ATTRIBUTE_UNUSED = IfThen(evolveA != 
      0,kmul(kmul(ksub(IfThen(advectLapse != 
      0,kmadd(betaxL,JacPDupwindNthAnti1alp,kmadd(betayL,JacPDupwindNthAnti2alp,kmadd(betazL,JacPDupwindNthAnti3alp,kmadd(JacPDupwindNthSymm1alp,kfabs(betaxL),kmadd(JacPDupwindNthSymm2alp,kfabs(betayL),kmul(JacPDupwindNthSymm3alp,kfabs(betazL))))))),ToReal(0)),dtalpL),kpow(alpL,-harmonicN)),ToReal(pow(harmonicF,-1))),ToReal(0));
    
    CCTK_REAL_VEC shiftGammaCoeffValue CCTK_ATTRIBUTE_UNUSED = 
      IfThen(useSpatialShiftGammaCoeff != 
      0,kmul(kfmin(ToReal(1),kexp(knmsub(ToReal(pow(spatialShiftGammaCoeffRadius,-1)),rL,ToReal(1)))),ToReal(shiftGammaCoeff)),ToReal(shiftGammaCoeff));
    
    CCTK_REAL_VEC B1L CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC B2L CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC B3L CCTK_ATTRIBUTE_UNUSED;
    
    if (evolveB != 0)
    {
      B1L = 
        kdiv(kmul(kpow(alpL,-shiftAlphaPower),ksub(dtbetaxL,IfThen(advectShift 
        != 
        0,kmadd(betaxL,JacPDupwindNthAnti1betax,kmadd(betayL,JacPDupwindNthAnti2betax,kmadd(betazL,JacPDupwindNthAnti3betax,kmadd(JacPDupwindNthSymm1betax,kfabs(betaxL),kmadd(JacPDupwindNthSymm2betax,kfabs(betayL),kmul(JacPDupwindNthSymm3betax,kfabs(betazL))))))),ToReal(0)))),shiftGammaCoeffValue);
      
      B2L = 
        kdiv(kmul(kpow(alpL,-shiftAlphaPower),ksub(dtbetayL,IfThen(advectShift 
        != 
        0,kmadd(betaxL,JacPDupwindNthAnti1betay,kmadd(betayL,JacPDupwindNthAnti2betay,kmadd(betazL,JacPDupwindNthAnti3betay,kmadd(JacPDupwindNthSymm1betay,kfabs(betaxL),kmadd(JacPDupwindNthSymm2betay,kfabs(betayL),kmul(JacPDupwindNthSymm3betay,kfabs(betazL))))))),ToReal(0)))),shiftGammaCoeffValue);
      
      B3L = 
        kdiv(kmul(kpow(alpL,-shiftAlphaPower),ksub(dtbetazL,IfThen(advectShift 
        != 
        0,kmadd(betaxL,JacPDupwindNthAnti1betaz,kmadd(betayL,JacPDupwindNthAnti2betaz,kmadd(betazL,JacPDupwindNthAnti3betaz,kmadd(JacPDupwindNthSymm1betaz,kfabs(betaxL),kmadd(JacPDupwindNthSymm2betaz,kfabs(betayL),kmul(JacPDupwindNthSymm3betaz,kfabs(betazL))))))),ToReal(0)))),shiftGammaCoeffValue);
    }
    else
    {
      B1L = ToReal(0);
      
      B2L = ToReal(0);
      
      B3L = ToReal(0);
    }
    /* Copy local copies back to grid functions */
    vec_store_partial_prepare(i,vecimin,vecimax);
    vec_store_nta_partial(A[index],AL);
    vec_store_nta_partial(B1[index],B1L);
    vec_store_nta_partial(B2[index],B2L);
    vec_store_nta_partial(B3[index],B3L);
    vec_store_nta_partial(Xt1[index],Xt1L);
    vec_store_nta_partial(Xt2[index],Xt2L);
    vec_store_nta_partial(Xt3[index],Xt3L);
  }
  CCTK_ENDLOOP3STR(ML_BSSN_InitialADMBase2Interior);
}
extern "C" void ML_BSSN_InitialADMBase2Interior(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_InitialADMBase2Interior_Body");
  }
  if (cctk_iteration % ML_BSSN_InitialADMBase2Interior_calc_every != ML_BSSN_InitialADMBase2Interior_calc_offset)
  {
    return;
  }
  
  const char* const groups[] = {
    "ADMBase::dtlapse",
    "ADMBase::dtshift",
    "ADMBase::lapse",
    "ADMBase::metric",
    "ADMBase::shift",
    "grid::coordinates",
    "ML_BSSN::ML_dtlapse",
    "ML_BSSN::ML_dtshift",
    "ML_BSSN::ML_Gamma",
    "ML_BSSN::ML_log_confac",
    "ML_BSSN::ML_metric"};
  AssertGroupStorage(cctkGH, "ML_BSSN_InitialADMBase2Interior", 11, groups);
  
  switch (fdOrder)
  {
    case 2:
    {
      EnsureStencilFits(cctkGH, "ML_BSSN_InitialADMBase2Interior", 2, 2, 2);
      break;
    }
    
    case 4:
    {
      EnsureStencilFits(cctkGH, "ML_BSSN_InitialADMBase2Interior", 3, 3, 3);
      break;
    }
    
    case 6:
    {
      EnsureStencilFits(cctkGH, "ML_BSSN_InitialADMBase2Interior", 4, 4, 4);
      break;
    }
    
    case 8:
    {
      EnsureStencilFits(cctkGH, "ML_BSSN_InitialADMBase2Interior", 5, 5, 5);
      break;
    }
    default:
      CCTK_BUILTIN_UNREACHABLE();
  }
  
  LoopOverInterior(cctkGH, ML_BSSN_InitialADMBase2Interior_Body);
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving ML_BSSN_InitialADMBase2Interior_Body");
  }
}

} // namespace ML_BSSN
