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

extern "C" void ML_BSSN_ADMBaseInterior_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (cctk_iteration % ML_BSSN_ADMBaseInterior_calc_every != ML_BSSN_ADMBaseInterior_calc_offset)
    return;
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ADMBase::dtlapse","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for ADMBase::dtlapse.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ADMBase::dtshift","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for ADMBase::dtshift.");
  return;
}

static void ML_BSSN_ADMBaseInterior_Body(const cGH* restrict const cctkGH, const int dir, const int face, const CCTK_REAL normal[3], const CCTK_REAL tangentA[3], const CCTK_REAL tangentB[3], const int imin[3], const int imax[3], const int n_subblock_gfs, CCTK_REAL* restrict const subblock_gfs[])
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
  CCTK_LOOP3STR(ML_BSSN_ADMBaseInterior,
    i,j,k, imin0,imin1,imin2, imax0,imax1,imax2,
    cctk_ash[0],cctk_ash[1],cctk_ash[2],
    vecimin,vecimax, CCTK_REAL_VEC_SIZE)
  {
    const ptrdiff_t index CCTK_ATTRIBUTE_UNUSED = di*i + dj*j + dk*k;
    /* Assign local copies of grid functions */
    
    CCTK_REAL_VEC AL CCTK_ATTRIBUTE_UNUSED = vec_load(A[index]);
    CCTK_REAL_VEC alphaL CCTK_ATTRIBUTE_UNUSED = vec_load(alpha[index]);
    CCTK_REAL_VEC B1L CCTK_ATTRIBUTE_UNUSED = vec_load(B1[index]);
    CCTK_REAL_VEC B2L CCTK_ATTRIBUTE_UNUSED = vec_load(B2[index]);
    CCTK_REAL_VEC B3L CCTK_ATTRIBUTE_UNUSED = vec_load(B3[index]);
    CCTK_REAL_VEC beta1L CCTK_ATTRIBUTE_UNUSED = vec_load(beta1[index]);
    CCTK_REAL_VEC beta2L CCTK_ATTRIBUTE_UNUSED = vec_load(beta2[index]);
    CCTK_REAL_VEC beta3L CCTK_ATTRIBUTE_UNUSED = vec_load(beta3[index]);
    CCTK_REAL_VEC gt11L CCTK_ATTRIBUTE_UNUSED = vec_load(gt11[index]);
    CCTK_REAL_VEC gt12L CCTK_ATTRIBUTE_UNUSED = vec_load(gt12[index]);
    CCTK_REAL_VEC gt13L CCTK_ATTRIBUTE_UNUSED = vec_load(gt13[index]);
    CCTK_REAL_VEC gt22L CCTK_ATTRIBUTE_UNUSED = vec_load(gt22[index]);
    CCTK_REAL_VEC gt23L CCTK_ATTRIBUTE_UNUSED = vec_load(gt23[index]);
    CCTK_REAL_VEC gt33L CCTK_ATTRIBUTE_UNUSED = vec_load(gt33[index]);
    CCTK_REAL_VEC phiL CCTK_ATTRIBUTE_UNUSED = vec_load(phi[index]);
    CCTK_REAL_VEC rL CCTK_ATTRIBUTE_UNUSED = vec_load(r[index]);
    CCTK_REAL_VEC trKL CCTK_ATTRIBUTE_UNUSED = vec_load(trK[index]);
    CCTK_REAL_VEC Xt1L CCTK_ATTRIBUTE_UNUSED = vec_load(Xt1[index]);
    CCTK_REAL_VEC Xt2L CCTK_ATTRIBUTE_UNUSED = vec_load(Xt2[index]);
    CCTK_REAL_VEC Xt3L CCTK_ATTRIBUTE_UNUSED = vec_load(Xt3[index]);
    
    
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
    CCTK_REAL_VEC PDstandardNth1alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm1alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm2alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm3alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti1alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti2alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti3alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm1beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm2beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm3beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti1beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti2beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti3beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm1beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm2beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm3beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti1beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti2beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti3beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm1beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm2beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthSymm3beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti1beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti2beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDupwindNthAnti3beta3 CCTK_ATTRIBUTE_UNUSED;
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
    CCTK_REAL_VEC PDstandardNth1phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3phi CCTK_ATTRIBUTE_UNUSED;
    
    switch (fdOrder)
    {
      case 2:
      {
        PDstandardNth1alpha = PDstandardNthfdOrder21(&alpha[index]);
        PDstandardNth2alpha = PDstandardNthfdOrder22(&alpha[index]);
        PDstandardNth3alpha = PDstandardNthfdOrder23(&alpha[index]);
        PDupwindNthSymm1alpha = PDupwindNthSymmfdOrder21(&alpha[index]);
        PDupwindNthSymm2alpha = PDupwindNthSymmfdOrder22(&alpha[index]);
        PDupwindNthSymm3alpha = PDupwindNthSymmfdOrder23(&alpha[index]);
        PDupwindNthAnti1alpha = PDupwindNthAntifdOrder21(&alpha[index]);
        PDupwindNthAnti2alpha = PDupwindNthAntifdOrder22(&alpha[index]);
        PDupwindNthAnti3alpha = PDupwindNthAntifdOrder23(&alpha[index]);
        PDupwindNthSymm1beta1 = PDupwindNthSymmfdOrder21(&beta1[index]);
        PDupwindNthSymm2beta1 = PDupwindNthSymmfdOrder22(&beta1[index]);
        PDupwindNthSymm3beta1 = PDupwindNthSymmfdOrder23(&beta1[index]);
        PDupwindNthAnti1beta1 = PDupwindNthAntifdOrder21(&beta1[index]);
        PDupwindNthAnti2beta1 = PDupwindNthAntifdOrder22(&beta1[index]);
        PDupwindNthAnti3beta1 = PDupwindNthAntifdOrder23(&beta1[index]);
        PDupwindNthSymm1beta2 = PDupwindNthSymmfdOrder21(&beta2[index]);
        PDupwindNthSymm2beta2 = PDupwindNthSymmfdOrder22(&beta2[index]);
        PDupwindNthSymm3beta2 = PDupwindNthSymmfdOrder23(&beta2[index]);
        PDupwindNthAnti1beta2 = PDupwindNthAntifdOrder21(&beta2[index]);
        PDupwindNthAnti2beta2 = PDupwindNthAntifdOrder22(&beta2[index]);
        PDupwindNthAnti3beta2 = PDupwindNthAntifdOrder23(&beta2[index]);
        PDupwindNthSymm1beta3 = PDupwindNthSymmfdOrder21(&beta3[index]);
        PDupwindNthSymm2beta3 = PDupwindNthSymmfdOrder22(&beta3[index]);
        PDupwindNthSymm3beta3 = PDupwindNthSymmfdOrder23(&beta3[index]);
        PDupwindNthAnti1beta3 = PDupwindNthAntifdOrder21(&beta3[index]);
        PDupwindNthAnti2beta3 = PDupwindNthAntifdOrder22(&beta3[index]);
        PDupwindNthAnti3beta3 = PDupwindNthAntifdOrder23(&beta3[index]);
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
        PDstandardNth1phi = PDstandardNthfdOrder21(&phi[index]);
        PDstandardNth2phi = PDstandardNthfdOrder22(&phi[index]);
        PDstandardNth3phi = PDstandardNthfdOrder23(&phi[index]);
        break;
      }
      
      case 4:
      {
        PDstandardNth1alpha = PDstandardNthfdOrder41(&alpha[index]);
        PDstandardNth2alpha = PDstandardNthfdOrder42(&alpha[index]);
        PDstandardNth3alpha = PDstandardNthfdOrder43(&alpha[index]);
        PDupwindNthSymm1alpha = PDupwindNthSymmfdOrder41(&alpha[index]);
        PDupwindNthSymm2alpha = PDupwindNthSymmfdOrder42(&alpha[index]);
        PDupwindNthSymm3alpha = PDupwindNthSymmfdOrder43(&alpha[index]);
        PDupwindNthAnti1alpha = PDupwindNthAntifdOrder41(&alpha[index]);
        PDupwindNthAnti2alpha = PDupwindNthAntifdOrder42(&alpha[index]);
        PDupwindNthAnti3alpha = PDupwindNthAntifdOrder43(&alpha[index]);
        PDupwindNthSymm1beta1 = PDupwindNthSymmfdOrder41(&beta1[index]);
        PDupwindNthSymm2beta1 = PDupwindNthSymmfdOrder42(&beta1[index]);
        PDupwindNthSymm3beta1 = PDupwindNthSymmfdOrder43(&beta1[index]);
        PDupwindNthAnti1beta1 = PDupwindNthAntifdOrder41(&beta1[index]);
        PDupwindNthAnti2beta1 = PDupwindNthAntifdOrder42(&beta1[index]);
        PDupwindNthAnti3beta1 = PDupwindNthAntifdOrder43(&beta1[index]);
        PDupwindNthSymm1beta2 = PDupwindNthSymmfdOrder41(&beta2[index]);
        PDupwindNthSymm2beta2 = PDupwindNthSymmfdOrder42(&beta2[index]);
        PDupwindNthSymm3beta2 = PDupwindNthSymmfdOrder43(&beta2[index]);
        PDupwindNthAnti1beta2 = PDupwindNthAntifdOrder41(&beta2[index]);
        PDupwindNthAnti2beta2 = PDupwindNthAntifdOrder42(&beta2[index]);
        PDupwindNthAnti3beta2 = PDupwindNthAntifdOrder43(&beta2[index]);
        PDupwindNthSymm1beta3 = PDupwindNthSymmfdOrder41(&beta3[index]);
        PDupwindNthSymm2beta3 = PDupwindNthSymmfdOrder42(&beta3[index]);
        PDupwindNthSymm3beta3 = PDupwindNthSymmfdOrder43(&beta3[index]);
        PDupwindNthAnti1beta3 = PDupwindNthAntifdOrder41(&beta3[index]);
        PDupwindNthAnti2beta3 = PDupwindNthAntifdOrder42(&beta3[index]);
        PDupwindNthAnti3beta3 = PDupwindNthAntifdOrder43(&beta3[index]);
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
        PDstandardNth1phi = PDstandardNthfdOrder41(&phi[index]);
        PDstandardNth2phi = PDstandardNthfdOrder42(&phi[index]);
        PDstandardNth3phi = PDstandardNthfdOrder43(&phi[index]);
        break;
      }
      
      case 6:
      {
        PDstandardNth1alpha = PDstandardNthfdOrder61(&alpha[index]);
        PDstandardNth2alpha = PDstandardNthfdOrder62(&alpha[index]);
        PDstandardNth3alpha = PDstandardNthfdOrder63(&alpha[index]);
        PDupwindNthSymm1alpha = PDupwindNthSymmfdOrder61(&alpha[index]);
        PDupwindNthSymm2alpha = PDupwindNthSymmfdOrder62(&alpha[index]);
        PDupwindNthSymm3alpha = PDupwindNthSymmfdOrder63(&alpha[index]);
        PDupwindNthAnti1alpha = PDupwindNthAntifdOrder61(&alpha[index]);
        PDupwindNthAnti2alpha = PDupwindNthAntifdOrder62(&alpha[index]);
        PDupwindNthAnti3alpha = PDupwindNthAntifdOrder63(&alpha[index]);
        PDupwindNthSymm1beta1 = PDupwindNthSymmfdOrder61(&beta1[index]);
        PDupwindNthSymm2beta1 = PDupwindNthSymmfdOrder62(&beta1[index]);
        PDupwindNthSymm3beta1 = PDupwindNthSymmfdOrder63(&beta1[index]);
        PDupwindNthAnti1beta1 = PDupwindNthAntifdOrder61(&beta1[index]);
        PDupwindNthAnti2beta1 = PDupwindNthAntifdOrder62(&beta1[index]);
        PDupwindNthAnti3beta1 = PDupwindNthAntifdOrder63(&beta1[index]);
        PDupwindNthSymm1beta2 = PDupwindNthSymmfdOrder61(&beta2[index]);
        PDupwindNthSymm2beta2 = PDupwindNthSymmfdOrder62(&beta2[index]);
        PDupwindNthSymm3beta2 = PDupwindNthSymmfdOrder63(&beta2[index]);
        PDupwindNthAnti1beta2 = PDupwindNthAntifdOrder61(&beta2[index]);
        PDupwindNthAnti2beta2 = PDupwindNthAntifdOrder62(&beta2[index]);
        PDupwindNthAnti3beta2 = PDupwindNthAntifdOrder63(&beta2[index]);
        PDupwindNthSymm1beta3 = PDupwindNthSymmfdOrder61(&beta3[index]);
        PDupwindNthSymm2beta3 = PDupwindNthSymmfdOrder62(&beta3[index]);
        PDupwindNthSymm3beta3 = PDupwindNthSymmfdOrder63(&beta3[index]);
        PDupwindNthAnti1beta3 = PDupwindNthAntifdOrder61(&beta3[index]);
        PDupwindNthAnti2beta3 = PDupwindNthAntifdOrder62(&beta3[index]);
        PDupwindNthAnti3beta3 = PDupwindNthAntifdOrder63(&beta3[index]);
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
        PDstandardNth1phi = PDstandardNthfdOrder61(&phi[index]);
        PDstandardNth2phi = PDstandardNthfdOrder62(&phi[index]);
        PDstandardNth3phi = PDstandardNthfdOrder63(&phi[index]);
        break;
      }
      
      case 8:
      {
        PDstandardNth1alpha = PDstandardNthfdOrder81(&alpha[index]);
        PDstandardNth2alpha = PDstandardNthfdOrder82(&alpha[index]);
        PDstandardNth3alpha = PDstandardNthfdOrder83(&alpha[index]);
        PDupwindNthSymm1alpha = PDupwindNthSymmfdOrder81(&alpha[index]);
        PDupwindNthSymm2alpha = PDupwindNthSymmfdOrder82(&alpha[index]);
        PDupwindNthSymm3alpha = PDupwindNthSymmfdOrder83(&alpha[index]);
        PDupwindNthAnti1alpha = PDupwindNthAntifdOrder81(&alpha[index]);
        PDupwindNthAnti2alpha = PDupwindNthAntifdOrder82(&alpha[index]);
        PDupwindNthAnti3alpha = PDupwindNthAntifdOrder83(&alpha[index]);
        PDupwindNthSymm1beta1 = PDupwindNthSymmfdOrder81(&beta1[index]);
        PDupwindNthSymm2beta1 = PDupwindNthSymmfdOrder82(&beta1[index]);
        PDupwindNthSymm3beta1 = PDupwindNthSymmfdOrder83(&beta1[index]);
        PDupwindNthAnti1beta1 = PDupwindNthAntifdOrder81(&beta1[index]);
        PDupwindNthAnti2beta1 = PDupwindNthAntifdOrder82(&beta1[index]);
        PDupwindNthAnti3beta1 = PDupwindNthAntifdOrder83(&beta1[index]);
        PDupwindNthSymm1beta2 = PDupwindNthSymmfdOrder81(&beta2[index]);
        PDupwindNthSymm2beta2 = PDupwindNthSymmfdOrder82(&beta2[index]);
        PDupwindNthSymm3beta2 = PDupwindNthSymmfdOrder83(&beta2[index]);
        PDupwindNthAnti1beta2 = PDupwindNthAntifdOrder81(&beta2[index]);
        PDupwindNthAnti2beta2 = PDupwindNthAntifdOrder82(&beta2[index]);
        PDupwindNthAnti3beta2 = PDupwindNthAntifdOrder83(&beta2[index]);
        PDupwindNthSymm1beta3 = PDupwindNthSymmfdOrder81(&beta3[index]);
        PDupwindNthSymm2beta3 = PDupwindNthSymmfdOrder82(&beta3[index]);
        PDupwindNthSymm3beta3 = PDupwindNthSymmfdOrder83(&beta3[index]);
        PDupwindNthAnti1beta3 = PDupwindNthAntifdOrder81(&beta3[index]);
        PDupwindNthAnti2beta3 = PDupwindNthAntifdOrder82(&beta3[index]);
        PDupwindNthAnti3beta3 = PDupwindNthAntifdOrder83(&beta3[index]);
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
        PDstandardNth1phi = PDstandardNthfdOrder81(&phi[index]);
        PDstandardNth2phi = PDstandardNthfdOrder82(&phi[index]);
        PDstandardNth3phi = PDstandardNthfdOrder83(&phi[index]);
        break;
      }
      default:
        CCTK_BUILTIN_UNREACHABLE();
    }
    /* Calculate temporaries and grid functions */
    CCTK_REAL_VEC JacPDstandardNth1alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti1alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti1beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti1beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti1beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti2alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti2beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti2beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti2beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti3alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti3beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti3beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthAnti3beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm1alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm1beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm1beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm1beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm2alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm2beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm2beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm2beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm3alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm3beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm3beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDupwindNthSymm3beta3 CCTK_ATTRIBUTE_UNUSED;
    
    if (usejacobian)
    {
      JacPDstandardNth1alpha = 
        kmadd(J11L,PDstandardNth1alpha,kmadd(J21L,PDstandardNth2alpha,kmul(J31L,PDstandardNth3alpha)));
      
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
      
      JacPDstandardNth1phi = 
        kmadd(J11L,PDstandardNth1phi,kmadd(J21L,PDstandardNth2phi,kmul(J31L,PDstandardNth3phi)));
      
      JacPDstandardNth2alpha = 
        kmadd(J12L,PDstandardNth1alpha,kmadd(J22L,PDstandardNth2alpha,kmul(J32L,PDstandardNth3alpha)));
      
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
      
      JacPDstandardNth2phi = 
        kmadd(J12L,PDstandardNth1phi,kmadd(J22L,PDstandardNth2phi,kmul(J32L,PDstandardNth3phi)));
      
      JacPDstandardNth3alpha = 
        kmadd(J13L,PDstandardNth1alpha,kmadd(J23L,PDstandardNth2alpha,kmul(J33L,PDstandardNth3alpha)));
      
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
      
      JacPDstandardNth3phi = 
        kmadd(J13L,PDstandardNth1phi,kmadd(J23L,PDstandardNth2phi,kmul(J33L,PDstandardNth3phi)));
      
      JacPDupwindNthSymm1alpha = 
        kmadd(J11L,PDupwindNthSymm1alpha,kmadd(J21L,PDupwindNthSymm2alpha,kmul(J31L,PDupwindNthSymm3alpha)));
      
      JacPDupwindNthSymm1beta1 = 
        kmadd(J11L,PDupwindNthSymm1beta1,kmadd(J21L,PDupwindNthSymm2beta1,kmul(J31L,PDupwindNthSymm3beta1)));
      
      JacPDupwindNthSymm1beta2 = 
        kmadd(J11L,PDupwindNthSymm1beta2,kmadd(J21L,PDupwindNthSymm2beta2,kmul(J31L,PDupwindNthSymm3beta2)));
      
      JacPDupwindNthSymm1beta3 = 
        kmadd(J11L,PDupwindNthSymm1beta3,kmadd(J21L,PDupwindNthSymm2beta3,kmul(J31L,PDupwindNthSymm3beta3)));
      
      JacPDupwindNthSymm2alpha = 
        kmadd(J12L,PDupwindNthSymm1alpha,kmadd(J22L,PDupwindNthSymm2alpha,kmul(J32L,PDupwindNthSymm3alpha)));
      
      JacPDupwindNthSymm2beta1 = 
        kmadd(J12L,PDupwindNthSymm1beta1,kmadd(J22L,PDupwindNthSymm2beta1,kmul(J32L,PDupwindNthSymm3beta1)));
      
      JacPDupwindNthSymm2beta2 = 
        kmadd(J12L,PDupwindNthSymm1beta2,kmadd(J22L,PDupwindNthSymm2beta2,kmul(J32L,PDupwindNthSymm3beta2)));
      
      JacPDupwindNthSymm2beta3 = 
        kmadd(J12L,PDupwindNthSymm1beta3,kmadd(J22L,PDupwindNthSymm2beta3,kmul(J32L,PDupwindNthSymm3beta3)));
      
      JacPDupwindNthSymm3alpha = 
        kmadd(J13L,PDupwindNthSymm1alpha,kmadd(J23L,PDupwindNthSymm2alpha,kmul(J33L,PDupwindNthSymm3alpha)));
      
      JacPDupwindNthSymm3beta1 = 
        kmadd(J13L,PDupwindNthSymm1beta1,kmadd(J23L,PDupwindNthSymm2beta1,kmul(J33L,PDupwindNthSymm3beta1)));
      
      JacPDupwindNthSymm3beta2 = 
        kmadd(J13L,PDupwindNthSymm1beta2,kmadd(J23L,PDupwindNthSymm2beta2,kmul(J33L,PDupwindNthSymm3beta2)));
      
      JacPDupwindNthSymm3beta3 = 
        kmadd(J13L,PDupwindNthSymm1beta3,kmadd(J23L,PDupwindNthSymm2beta3,kmul(J33L,PDupwindNthSymm3beta3)));
      
      JacPDupwindNthAnti1alpha = 
        kmadd(J11L,PDupwindNthAnti1alpha,kmadd(J21L,PDupwindNthAnti2alpha,kmul(J31L,PDupwindNthAnti3alpha)));
      
      JacPDupwindNthAnti1beta1 = 
        kmadd(J11L,PDupwindNthAnti1beta1,kmadd(J21L,PDupwindNthAnti2beta1,kmul(J31L,PDupwindNthAnti3beta1)));
      
      JacPDupwindNthAnti1beta2 = 
        kmadd(J11L,PDupwindNthAnti1beta2,kmadd(J21L,PDupwindNthAnti2beta2,kmul(J31L,PDupwindNthAnti3beta2)));
      
      JacPDupwindNthAnti1beta3 = 
        kmadd(J11L,PDupwindNthAnti1beta3,kmadd(J21L,PDupwindNthAnti2beta3,kmul(J31L,PDupwindNthAnti3beta3)));
      
      JacPDupwindNthAnti2alpha = 
        kmadd(J12L,PDupwindNthAnti1alpha,kmadd(J22L,PDupwindNthAnti2alpha,kmul(J32L,PDupwindNthAnti3alpha)));
      
      JacPDupwindNthAnti2beta1 = 
        kmadd(J12L,PDupwindNthAnti1beta1,kmadd(J22L,PDupwindNthAnti2beta1,kmul(J32L,PDupwindNthAnti3beta1)));
      
      JacPDupwindNthAnti2beta2 = 
        kmadd(J12L,PDupwindNthAnti1beta2,kmadd(J22L,PDupwindNthAnti2beta2,kmul(J32L,PDupwindNthAnti3beta2)));
      
      JacPDupwindNthAnti2beta3 = 
        kmadd(J12L,PDupwindNthAnti1beta3,kmadd(J22L,PDupwindNthAnti2beta3,kmul(J32L,PDupwindNthAnti3beta3)));
      
      JacPDupwindNthAnti3alpha = 
        kmadd(J13L,PDupwindNthAnti1alpha,kmadd(J23L,PDupwindNthAnti2alpha,kmul(J33L,PDupwindNthAnti3alpha)));
      
      JacPDupwindNthAnti3beta1 = 
        kmadd(J13L,PDupwindNthAnti1beta1,kmadd(J23L,PDupwindNthAnti2beta1,kmul(J33L,PDupwindNthAnti3beta1)));
      
      JacPDupwindNthAnti3beta2 = 
        kmadd(J13L,PDupwindNthAnti1beta2,kmadd(J23L,PDupwindNthAnti2beta2,kmul(J33L,PDupwindNthAnti3beta2)));
      
      JacPDupwindNthAnti3beta3 = 
        kmadd(J13L,PDupwindNthAnti1beta3,kmadd(J23L,PDupwindNthAnti2beta3,kmul(J33L,PDupwindNthAnti3beta3)));
    }
    else
    {
      JacPDstandardNth1alpha = PDstandardNth1alpha;
      
      JacPDstandardNth1gt11 = PDstandardNth1gt11;
      
      JacPDstandardNth1gt12 = PDstandardNth1gt12;
      
      JacPDstandardNth1gt13 = PDstandardNth1gt13;
      
      JacPDstandardNth1gt22 = PDstandardNth1gt22;
      
      JacPDstandardNth1gt23 = PDstandardNth1gt23;
      
      JacPDstandardNth1gt33 = PDstandardNth1gt33;
      
      JacPDstandardNth1phi = PDstandardNth1phi;
      
      JacPDstandardNth2alpha = PDstandardNth2alpha;
      
      JacPDstandardNth2gt11 = PDstandardNth2gt11;
      
      JacPDstandardNth2gt12 = PDstandardNth2gt12;
      
      JacPDstandardNth2gt13 = PDstandardNth2gt13;
      
      JacPDstandardNth2gt22 = PDstandardNth2gt22;
      
      JacPDstandardNth2gt23 = PDstandardNth2gt23;
      
      JacPDstandardNth2gt33 = PDstandardNth2gt33;
      
      JacPDstandardNth2phi = PDstandardNth2phi;
      
      JacPDstandardNth3alpha = PDstandardNth3alpha;
      
      JacPDstandardNth3gt11 = PDstandardNth3gt11;
      
      JacPDstandardNth3gt12 = PDstandardNth3gt12;
      
      JacPDstandardNth3gt13 = PDstandardNth3gt13;
      
      JacPDstandardNth3gt22 = PDstandardNth3gt22;
      
      JacPDstandardNth3gt23 = PDstandardNth3gt23;
      
      JacPDstandardNth3gt33 = PDstandardNth3gt33;
      
      JacPDstandardNth3phi = PDstandardNth3phi;
      
      JacPDupwindNthSymm1alpha = PDupwindNthSymm1alpha;
      
      JacPDupwindNthSymm1beta1 = PDupwindNthSymm1beta1;
      
      JacPDupwindNthSymm1beta2 = PDupwindNthSymm1beta2;
      
      JacPDupwindNthSymm1beta3 = PDupwindNthSymm1beta3;
      
      JacPDupwindNthSymm2alpha = PDupwindNthSymm2alpha;
      
      JacPDupwindNthSymm2beta1 = PDupwindNthSymm2beta1;
      
      JacPDupwindNthSymm2beta2 = PDupwindNthSymm2beta2;
      
      JacPDupwindNthSymm2beta3 = PDupwindNthSymm2beta3;
      
      JacPDupwindNthSymm3alpha = PDupwindNthSymm3alpha;
      
      JacPDupwindNthSymm3beta1 = PDupwindNthSymm3beta1;
      
      JacPDupwindNthSymm3beta2 = PDupwindNthSymm3beta2;
      
      JacPDupwindNthSymm3beta3 = PDupwindNthSymm3beta3;
      
      JacPDupwindNthAnti1alpha = PDupwindNthAnti1alpha;
      
      JacPDupwindNthAnti1beta1 = PDupwindNthAnti1beta1;
      
      JacPDupwindNthAnti1beta2 = PDupwindNthAnti1beta2;
      
      JacPDupwindNthAnti1beta3 = PDupwindNthAnti1beta3;
      
      JacPDupwindNthAnti2alpha = PDupwindNthAnti2alpha;
      
      JacPDupwindNthAnti2beta1 = PDupwindNthAnti2beta1;
      
      JacPDupwindNthAnti2beta2 = PDupwindNthAnti2beta2;
      
      JacPDupwindNthAnti2beta3 = PDupwindNthAnti2beta3;
      
      JacPDupwindNthAnti3alpha = PDupwindNthAnti3alpha;
      
      JacPDupwindNthAnti3beta1 = PDupwindNthAnti3beta1;
      
      JacPDupwindNthAnti3beta2 = PDupwindNthAnti3beta2;
      
      JacPDupwindNthAnti3beta3 = PDupwindNthAnti3beta3;
    }
    
    CCTK_REAL_VEC detgt CCTK_ATTRIBUTE_UNUSED = ToReal(1);
    
    CCTK_REAL_VEC gtu11 CCTK_ATTRIBUTE_UNUSED = 
      kdiv(kmsub(gt22L,gt33L,kmul(gt23L,gt23L)),detgt);
    
    CCTK_REAL_VEC gtu12 CCTK_ATTRIBUTE_UNUSED = 
      kdiv(kmsub(gt13L,gt23L,kmul(gt12L,gt33L)),detgt);
    
    CCTK_REAL_VEC gtu13 CCTK_ATTRIBUTE_UNUSED = 
      kdiv(kmsub(gt12L,gt23L,kmul(gt13L,gt22L)),detgt);
    
    CCTK_REAL_VEC gtu22 CCTK_ATTRIBUTE_UNUSED = 
      kdiv(kmsub(gt11L,gt33L,kmul(gt13L,gt13L)),detgt);
    
    CCTK_REAL_VEC gtu23 CCTK_ATTRIBUTE_UNUSED = 
      kdiv(kmsub(gt12L,gt13L,kmul(gt11L,gt23L)),detgt);
    
    CCTK_REAL_VEC gtu33 CCTK_ATTRIBUTE_UNUSED = 
      kdiv(kmsub(gt11L,gt22L,kmul(gt12L,gt12L)),detgt);
    
    CCTK_REAL_VEC em4phi CCTK_ATTRIBUTE_UNUSED = IfThen(conformalMethod != 
      0,kmul(phiL,phiL),kexp(kmul(phiL,ToReal(-4))));
    
    CCTK_REAL_VEC gu11 CCTK_ATTRIBUTE_UNUSED = kmul(em4phi,gtu11);
    
    CCTK_REAL_VEC gu12 CCTK_ATTRIBUTE_UNUSED = kmul(em4phi,gtu12);
    
    CCTK_REAL_VEC gu13 CCTK_ATTRIBUTE_UNUSED = kmul(em4phi,gtu13);
    
    CCTK_REAL_VEC gu22 CCTK_ATTRIBUTE_UNUSED = kmul(em4phi,gtu22);
    
    CCTK_REAL_VEC gu23 CCTK_ATTRIBUTE_UNUSED = kmul(em4phi,gtu23);
    
    CCTK_REAL_VEC gu33 CCTK_ATTRIBUTE_UNUSED = kmul(em4phi,gtu33);
    
    CCTK_REAL_VEC fac1 CCTK_ATTRIBUTE_UNUSED = IfThen(conformalMethod != 
      0,kdiv(ToReal(-0.5),phiL),ToReal(1));
    
    CCTK_REAL_VEC cdphi1 CCTK_ATTRIBUTE_UNUSED = 
      kmul(fac1,JacPDstandardNth1phi);
    
    CCTK_REAL_VEC cdphi2 CCTK_ATTRIBUTE_UNUSED = 
      kmul(fac1,JacPDstandardNth2phi);
    
    CCTK_REAL_VEC cdphi3 CCTK_ATTRIBUTE_UNUSED = 
      kmul(fac1,JacPDstandardNth3phi);
    
    CCTK_REAL_VEC dotalpha CCTK_ATTRIBUTE_UNUSED = 
      kneg(kmul(kmul(IfThen(evolveA != 
      0,AL,kmadd(ToReal(alphaDriver),kadd(ToReal(-1),alphaL),trKL)),kpow(alphaL,harmonicN)),ToReal(harmonicF)));
    
    CCTK_REAL_VEC betaDriverValue CCTK_ATTRIBUTE_UNUSED = 
      IfThen(useSpatialBetaDriver != 
      0,kdiv(ToReal(betaDriver*spatialBetaDriverRadius),kfmax(rL,ToReal(spatialBetaDriverRadius))),ToReal(betaDriver));
    
    CCTK_REAL_VEC shiftGammaCoeffValue CCTK_ATTRIBUTE_UNUSED = 
      IfThen(useSpatialShiftGammaCoeff != 
      0,kmul(kfmin(ToReal(1),kexp(knmsub(ToReal(pow(spatialShiftGammaCoeffRadius,-1)),rL,ToReal(1)))),ToReal(shiftGammaCoeff)),ToReal(shiftGammaCoeff));
    
    CCTK_REAL_VEC ddetgt1 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(gtu11,JacPDstandardNth1gt11,kmadd(gtu22,JacPDstandardNth1gt22,kmadd(ToReal(2),kmadd(gtu12,JacPDstandardNth1gt12,kmadd(gtu13,JacPDstandardNth1gt13,kmul(gtu23,JacPDstandardNth1gt23))),kmul(gtu33,JacPDstandardNth1gt33))));
    
    CCTK_REAL_VEC ddetgt2 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(gtu11,JacPDstandardNth2gt11,kmadd(gtu22,JacPDstandardNth2gt22,kmadd(ToReal(2),kmadd(gtu12,JacPDstandardNth2gt12,kmadd(gtu13,JacPDstandardNth2gt13,kmul(gtu23,JacPDstandardNth2gt23))),kmul(gtu33,JacPDstandardNth2gt33))));
    
    CCTK_REAL_VEC ddetgt3 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(gtu11,JacPDstandardNth3gt11,kmadd(gtu22,JacPDstandardNth3gt22,kmadd(ToReal(2),kmadd(gtu12,JacPDstandardNth3gt12,kmadd(gtu13,JacPDstandardNth3gt13,kmul(gtu23,JacPDstandardNth3gt23))),kmul(gtu33,JacPDstandardNth3gt33))));
    
    CCTK_REAL_VEC dotbeta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC dotbeta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC dotbeta3 CCTK_ATTRIBUTE_UNUSED;
    
    if (shiftFormulation == 0)
    {
      dotbeta1 = kmul(shiftGammaCoeffValue,kmul(IfThen(evolveB != 
        0,B1L,knmsub(beta1L,betaDriverValue,Xt1L)),kpow(alphaL,shiftAlphaPower)));
      
      dotbeta2 = kmul(shiftGammaCoeffValue,kmul(IfThen(evolveB != 
        0,B2L,knmsub(beta2L,betaDriverValue,Xt2L)),kpow(alphaL,shiftAlphaPower)));
      
      dotbeta3 = kmul(shiftGammaCoeffValue,kmul(IfThen(evolveB != 
        0,B3L,knmsub(beta3L,betaDriverValue,Xt3L)),kpow(alphaL,shiftAlphaPower)));
    }
    else
    {
      dotbeta1 = 
        kneg(kmul(alphaL,kmadd(gu11,kmadd(alphaL,kmadd(ToReal(2),cdphi1,kmsub(ToReal(0.5),ddetgt1,kmadd(gtu11,JacPDstandardNth1gt11,kmadd(gtu12,kadd(JacPDstandardNth1gt12,JacPDstandardNth2gt11),kmadd(gtu22,JacPDstandardNth2gt12,kmadd(gtu13,kadd(JacPDstandardNth1gt13,JacPDstandardNth3gt11),kmadd(gtu23,kadd(JacPDstandardNth2gt13,JacPDstandardNth3gt12),kmul(gtu33,JacPDstandardNth3gt13)))))))),JacPDstandardNth1alpha),kmadd(gu12,kmadd(alphaL,kmadd(ToReal(2),cdphi2,kmsub(ToReal(0.5),ddetgt2,kmadd(gtu11,JacPDstandardNth1gt12,kmadd(gtu12,kadd(JacPDstandardNth1gt22,JacPDstandardNth2gt12),kmadd(gtu22,JacPDstandardNth2gt22,kmadd(gtu13,kadd(JacPDstandardNth1gt23,JacPDstandardNth3gt12),kmadd(gtu23,kadd(JacPDstandardNth2gt23,JacPDstandardNth3gt22),kmul(gtu33,JacPDstandardNth3gt23)))))))),JacPDstandardNth2alpha),kmul(gu13,kmadd(alphaL,kmadd(ToReal(2),cdphi3,kmsub(ToReal(0.5),ddetgt3,kmadd(gtu11,JacPDstandardNth1gt13,kmadd(gtu12,kadd(JacPDstandardNth1gt23,JacPDstandardNth2gt13),kmadd(gtu22,JacPDstandardNth2gt23,kmadd(gtu13,kadd(JacPDstandardNth1gt33,JacPDstandardNth3gt13),kmadd(gtu23,kadd(JacPDstandardNth2gt33,JacPDstandardNth3gt23),kmul(gtu33,JacPDstandardNth3gt33)))))))),JacPDstandardNth3alpha))))));
      
      dotbeta2 = 
        kneg(kmul(alphaL,kmadd(gu12,kmadd(alphaL,kmadd(ToReal(2),cdphi1,kmsub(ToReal(0.5),ddetgt1,kmadd(gtu11,JacPDstandardNth1gt11,kmadd(gtu12,kadd(JacPDstandardNth1gt12,JacPDstandardNth2gt11),kmadd(gtu22,JacPDstandardNth2gt12,kmadd(gtu13,kadd(JacPDstandardNth1gt13,JacPDstandardNth3gt11),kmadd(gtu23,kadd(JacPDstandardNth2gt13,JacPDstandardNth3gt12),kmul(gtu33,JacPDstandardNth3gt13)))))))),JacPDstandardNth1alpha),kmadd(gu22,kmadd(alphaL,kmadd(ToReal(2),cdphi2,kmsub(ToReal(0.5),ddetgt2,kmadd(gtu11,JacPDstandardNth1gt12,kmadd(gtu12,kadd(JacPDstandardNth1gt22,JacPDstandardNth2gt12),kmadd(gtu22,JacPDstandardNth2gt22,kmadd(gtu13,kadd(JacPDstandardNth1gt23,JacPDstandardNth3gt12),kmadd(gtu23,kadd(JacPDstandardNth2gt23,JacPDstandardNth3gt22),kmul(gtu33,JacPDstandardNth3gt23)))))))),JacPDstandardNth2alpha),kmul(gu23,kmadd(alphaL,kmadd(ToReal(2),cdphi3,kmsub(ToReal(0.5),ddetgt3,kmadd(gtu11,JacPDstandardNth1gt13,kmadd(gtu12,kadd(JacPDstandardNth1gt23,JacPDstandardNth2gt13),kmadd(gtu22,JacPDstandardNth2gt23,kmadd(gtu13,kadd(JacPDstandardNth1gt33,JacPDstandardNth3gt13),kmadd(gtu23,kadd(JacPDstandardNth2gt33,JacPDstandardNth3gt23),kmul(gtu33,JacPDstandardNth3gt33)))))))),JacPDstandardNth3alpha))))));
      
      dotbeta3 = 
        kneg(kmul(alphaL,kmadd(gu13,kmadd(alphaL,kmadd(ToReal(2),cdphi1,kmsub(ToReal(0.5),ddetgt1,kmadd(gtu11,JacPDstandardNth1gt11,kmadd(gtu12,kadd(JacPDstandardNth1gt12,JacPDstandardNth2gt11),kmadd(gtu22,JacPDstandardNth2gt12,kmadd(gtu13,kadd(JacPDstandardNth1gt13,JacPDstandardNth3gt11),kmadd(gtu23,kadd(JacPDstandardNth2gt13,JacPDstandardNth3gt12),kmul(gtu33,JacPDstandardNth3gt13)))))))),JacPDstandardNth1alpha),kmadd(gu23,kmadd(alphaL,kmadd(ToReal(2),cdphi2,kmsub(ToReal(0.5),ddetgt2,kmadd(gtu11,JacPDstandardNth1gt12,kmadd(gtu12,kadd(JacPDstandardNth1gt22,JacPDstandardNth2gt12),kmadd(gtu22,JacPDstandardNth2gt22,kmadd(gtu13,kadd(JacPDstandardNth1gt23,JacPDstandardNth3gt12),kmadd(gtu23,kadd(JacPDstandardNth2gt23,JacPDstandardNth3gt22),kmul(gtu33,JacPDstandardNth3gt23)))))))),JacPDstandardNth2alpha),kmul(gu33,kmadd(alphaL,kmadd(ToReal(2),cdphi3,kmsub(ToReal(0.5),ddetgt3,kmadd(gtu11,JacPDstandardNth1gt13,kmadd(gtu12,kadd(JacPDstandardNth1gt23,JacPDstandardNth2gt13),kmadd(gtu22,JacPDstandardNth2gt23,kmadd(gtu13,kadd(JacPDstandardNth1gt33,JacPDstandardNth3gt13),kmadd(gtu23,kadd(JacPDstandardNth2gt33,JacPDstandardNth3gt23),kmul(gtu33,JacPDstandardNth3gt33)))))))),JacPDstandardNth3alpha))))));
    }
    
    CCTK_REAL_VEC dtalpL CCTK_ATTRIBUTE_UNUSED = 
      kadd(dotalpha,IfThen(advectLapse != 
      0,kmadd(beta1L,JacPDupwindNthAnti1alpha,kmadd(beta2L,JacPDupwindNthAnti2alpha,kmadd(beta3L,JacPDupwindNthAnti3alpha,kmadd(JacPDupwindNthSymm1alpha,kfabs(beta1L),kmadd(JacPDupwindNthSymm2alpha,kfabs(beta2L),kmul(JacPDupwindNthSymm3alpha,kfabs(beta3L))))))),ToReal(0)));
    
    CCTK_REAL_VEC dtbetaxL CCTK_ATTRIBUTE_UNUSED = 
      kadd(dotbeta1,IfThen(advectShift != 
      0,kmadd(beta1L,JacPDupwindNthAnti1beta1,kmadd(beta2L,JacPDupwindNthAnti2beta1,kmadd(beta3L,JacPDupwindNthAnti3beta1,kmadd(JacPDupwindNthSymm1beta1,kfabs(beta1L),kmadd(JacPDupwindNthSymm2beta1,kfabs(beta2L),kmul(JacPDupwindNthSymm3beta1,kfabs(beta3L))))))),ToReal(0)));
    
    CCTK_REAL_VEC dtbetayL CCTK_ATTRIBUTE_UNUSED = 
      kadd(dotbeta2,IfThen(advectShift != 
      0,kmadd(beta1L,JacPDupwindNthAnti1beta2,kmadd(beta2L,JacPDupwindNthAnti2beta2,kmadd(beta3L,JacPDupwindNthAnti3beta2,kmadd(JacPDupwindNthSymm1beta2,kfabs(beta1L),kmadd(JacPDupwindNthSymm2beta2,kfabs(beta2L),kmul(JacPDupwindNthSymm3beta2,kfabs(beta3L))))))),ToReal(0)));
    
    CCTK_REAL_VEC dtbetazL CCTK_ATTRIBUTE_UNUSED = 
      kadd(dotbeta3,IfThen(advectShift != 
      0,kmadd(beta1L,JacPDupwindNthAnti1beta3,kmadd(beta2L,JacPDupwindNthAnti2beta3,kmadd(beta3L,JacPDupwindNthAnti3beta3,kmadd(JacPDupwindNthSymm1beta3,kfabs(beta1L),kmadd(JacPDupwindNthSymm2beta3,kfabs(beta2L),kmul(JacPDupwindNthSymm3beta3,kfabs(beta3L))))))),ToReal(0)));
    /* Copy local copies back to grid functions */
    vec_store_partial_prepare(i,vecimin,vecimax);
    vec_store_nta_partial(dtalp[index],dtalpL);
    vec_store_nta_partial(dtbetax[index],dtbetaxL);
    vec_store_nta_partial(dtbetay[index],dtbetayL);
    vec_store_nta_partial(dtbetaz[index],dtbetazL);
  }
  CCTK_ENDLOOP3STR(ML_BSSN_ADMBaseInterior);
}
extern "C" void ML_BSSN_ADMBaseInterior(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_ADMBaseInterior_Body");
  }
  if (cctk_iteration % ML_BSSN_ADMBaseInterior_calc_every != ML_BSSN_ADMBaseInterior_calc_offset)
  {
    return;
  }
  
  const char* const groups[] = {
    "ADMBase::dtlapse",
    "ADMBase::dtshift",
    "grid::coordinates",
    "ML_BSSN::ML_dtlapse",
    "ML_BSSN::ML_dtshift",
    "ML_BSSN::ML_Gamma",
    "ML_BSSN::ML_lapse",
    "ML_BSSN::ML_log_confac",
    "ML_BSSN::ML_metric",
    "ML_BSSN::ML_shift",
    "ML_BSSN::ML_trace_curv"};
  AssertGroupStorage(cctkGH, "ML_BSSN_ADMBaseInterior", 11, groups);
  
  switch (fdOrder)
  {
    case 2:
    {
      EnsureStencilFits(cctkGH, "ML_BSSN_ADMBaseInterior", 2, 2, 2);
      break;
    }
    
    case 4:
    {
      EnsureStencilFits(cctkGH, "ML_BSSN_ADMBaseInterior", 3, 3, 3);
      break;
    }
    
    case 6:
    {
      EnsureStencilFits(cctkGH, "ML_BSSN_ADMBaseInterior", 4, 4, 4);
      break;
    }
    
    case 8:
    {
      EnsureStencilFits(cctkGH, "ML_BSSN_ADMBaseInterior", 5, 5, 5);
      break;
    }
    default:
      CCTK_BUILTIN_UNREACHABLE();
  }
  
  LoopOverInterior(cctkGH, ML_BSSN_ADMBaseInterior_Body);
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving ML_BSSN_ADMBaseInterior_Body");
  }
}

} // namespace ML_BSSN
