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

namespace WeylScal4 {

extern "C" void WeylScal4_psis_calc_4th_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (cctk_iteration % WeylScal4_psis_calc_4th_calc_every != WeylScal4_psis_calc_4th_calc_offset)
    return;
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "WeylScal4::Psi0i_group","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for WeylScal4::Psi0i_group.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "WeylScal4::Psi0r_group","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for WeylScal4::Psi0r_group.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "WeylScal4::Psi1i_group","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for WeylScal4::Psi1i_group.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "WeylScal4::Psi1r_group","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for WeylScal4::Psi1r_group.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "WeylScal4::Psi2i_group","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for WeylScal4::Psi2i_group.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "WeylScal4::Psi2r_group","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for WeylScal4::Psi2r_group.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "WeylScal4::Psi3i_group","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for WeylScal4::Psi3i_group.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "WeylScal4::Psi3r_group","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for WeylScal4::Psi3r_group.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "WeylScal4::Psi4i_group","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for WeylScal4::Psi4i_group.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "WeylScal4::Psi4r_group","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for WeylScal4::Psi4r_group.");
  return;
}

static void WeylScal4_psis_calc_4th_Body(const cGH* restrict const cctkGH, const int dir, const int face, const CCTK_REAL normal[3], const CCTK_REAL tangentA[3], const CCTK_REAL tangentB[3], const int imin[3], const int imax[3], const int n_subblock_gfs, CCTK_REAL* restrict const subblock_gfs[])
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
  const CCTK_REAL_VEC p1o12dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0833333333333333333333333333333),dx);
  const CCTK_REAL_VEC p1o12dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0833333333333333333333333333333),dy);
  const CCTK_REAL_VEC p1o12dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0833333333333333333333333333333),dz);
  const CCTK_REAL_VEC p1o144dxdy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00694444444444444444444444444444),kmul(dx,dy));
  const CCTK_REAL_VEC p1o144dxdz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00694444444444444444444444444444),kmul(dx,dz));
  const CCTK_REAL_VEC p1o144dydz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00694444444444444444444444444444),kmul(dy,dz));
  const CCTK_REAL_VEC p1o180dx2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00555555555555555555555555555556),kmul(dx,dx));
  const CCTK_REAL_VEC p1o180dy2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00555555555555555555555555555556),kmul(dy,dy));
  const CCTK_REAL_VEC p1o180dz2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00555555555555555555555555555556),kmul(dz,dz));
  const CCTK_REAL_VEC p1o2dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.5),dx);
  const CCTK_REAL_VEC p1o2dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.5),dy);
  const CCTK_REAL_VEC p1o2dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.5),dz);
  const CCTK_REAL_VEC p1o3600dxdy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.000277777777777777777777777777778),kmul(dx,dy));
  const CCTK_REAL_VEC p1o3600dxdz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.000277777777777777777777777777778),kmul(dx,dz));
  const CCTK_REAL_VEC p1o3600dydz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.000277777777777777777777777777778),kmul(dy,dz));
  const CCTK_REAL_VEC p1o4dxdy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.25),kmul(dx,dy));
  const CCTK_REAL_VEC p1o4dxdz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.25),kmul(dx,dz));
  const CCTK_REAL_VEC p1o4dydz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.25),kmul(dy,dz));
  const CCTK_REAL_VEC p1o5040dx2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.000198412698412698412698412698413),kmul(dx,dx));
  const CCTK_REAL_VEC p1o5040dy2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.000198412698412698412698412698413),kmul(dy,dy));
  const CCTK_REAL_VEC p1o5040dz2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.000198412698412698412698412698413),kmul(dz,dz));
  const CCTK_REAL_VEC p1o60dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0166666666666666666666666666667),dx);
  const CCTK_REAL_VEC p1o60dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0166666666666666666666666666667),dy);
  const CCTK_REAL_VEC p1o60dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0166666666666666666666666666667),dz);
  const CCTK_REAL_VEC p1o705600dxdy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1.41723356009070294784580498866e-6),kmul(dx,dy));
  const CCTK_REAL_VEC p1o705600dxdz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1.41723356009070294784580498866e-6),kmul(dx,dz));
  const CCTK_REAL_VEC p1o705600dydz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1.41723356009070294784580498866e-6),kmul(dy,dz));
  const CCTK_REAL_VEC p1o840dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00119047619047619047619047619048),dx);
  const CCTK_REAL_VEC p1o840dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00119047619047619047619047619048),dy);
  const CCTK_REAL_VEC p1o840dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00119047619047619047619047619048),dz);
  const CCTK_REAL_VEC p1odx2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1),kmul(dx,dx));
  const CCTK_REAL_VEC p1ody2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1),kmul(dy,dy));
  const CCTK_REAL_VEC p1odz2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1),kmul(dz,dz));
  const CCTK_REAL_VEC pm1o12dx2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.0833333333333333333333333333333),kmul(dx,dx));
  const CCTK_REAL_VEC pm1o12dy2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.0833333333333333333333333333333),kmul(dy,dy));
  const CCTK_REAL_VEC pm1o12dz2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.0833333333333333333333333333333),kmul(dz,dz));
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
  CCTK_LOOP3STR(WeylScal4_psis_calc_4th,
    i,j,k, imin0,imin1,imin2, imax0,imax1,imax2,
    cctk_ash[0],cctk_ash[1],cctk_ash[2],
    vecimin,vecimax, CCTK_REAL_VEC_SIZE)
  {
    const ptrdiff_t index CCTK_ATTRIBUTE_UNUSED = di*i + dj*j + dk*k;
    /* Assign local copies of grid functions */
    
    CCTK_REAL_VEC gxxL CCTK_ATTRIBUTE_UNUSED = vec_load(gxx[index]);
    CCTK_REAL_VEC gxyL CCTK_ATTRIBUTE_UNUSED = vec_load(gxy[index]);
    CCTK_REAL_VEC gxzL CCTK_ATTRIBUTE_UNUSED = vec_load(gxz[index]);
    CCTK_REAL_VEC gyyL CCTK_ATTRIBUTE_UNUSED = vec_load(gyy[index]);
    CCTK_REAL_VEC gyzL CCTK_ATTRIBUTE_UNUSED = vec_load(gyz[index]);
    CCTK_REAL_VEC gzzL CCTK_ATTRIBUTE_UNUSED = vec_load(gzz[index]);
    CCTK_REAL_VEC kxxL CCTK_ATTRIBUTE_UNUSED = vec_load(kxx[index]);
    CCTK_REAL_VEC kxyL CCTK_ATTRIBUTE_UNUSED = vec_load(kxy[index]);
    CCTK_REAL_VEC kxzL CCTK_ATTRIBUTE_UNUSED = vec_load(kxz[index]);
    CCTK_REAL_VEC kyyL CCTK_ATTRIBUTE_UNUSED = vec_load(kyy[index]);
    CCTK_REAL_VEC kyzL CCTK_ATTRIBUTE_UNUSED = vec_load(kyz[index]);
    CCTK_REAL_VEC kzzL CCTK_ATTRIBUTE_UNUSED = vec_load(kzz[index]);
    CCTK_REAL_VEC xL CCTK_ATTRIBUTE_UNUSED = vec_load(x[index]);
    CCTK_REAL_VEC yL CCTK_ATTRIBUTE_UNUSED = vec_load(y[index]);
    CCTK_REAL_VEC zL CCTK_ATTRIBUTE_UNUSED = vec_load(z[index]);
    
    
    CCTK_REAL_VEC dJ111L, dJ112L, dJ113L, dJ122L, dJ123L, dJ133L, dJ211L, dJ212L, dJ213L, dJ222L, dJ223L, dJ233L, dJ311L, dJ312L, dJ313L, dJ322L, dJ323L, dJ333L, J11L, J12L, J13L, J21L, J22L, J23L, J31L, J32L, J33L CCTK_ATTRIBUTE_UNUSED ;
    
    if (usejacobian)
    {
      dJ111L = vec_load(dJ111[index]);
      dJ112L = vec_load(dJ112[index]);
      dJ113L = vec_load(dJ113[index]);
      dJ122L = vec_load(dJ122[index]);
      dJ123L = vec_load(dJ123[index]);
      dJ133L = vec_load(dJ133[index]);
      dJ211L = vec_load(dJ211[index]);
      dJ212L = vec_load(dJ212[index]);
      dJ213L = vec_load(dJ213[index]);
      dJ222L = vec_load(dJ222[index]);
      dJ223L = vec_load(dJ223[index]);
      dJ233L = vec_load(dJ233[index]);
      dJ311L = vec_load(dJ311[index]);
      dJ312L = vec_load(dJ312[index]);
      dJ313L = vec_load(dJ313[index]);
      dJ322L = vec_load(dJ322[index]);
      dJ323L = vec_load(dJ323[index]);
      dJ333L = vec_load(dJ333[index]);
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
    CCTK_REAL_VEC PDstandard4th1gxx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th2gxx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th3gxx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th11gxx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th22gxx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th33gxx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th12gxx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th13gxx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th23gxx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th1gxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th2gxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th3gxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th11gxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th22gxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th33gxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th12gxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th13gxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th23gxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th1gxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th2gxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th3gxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th11gxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th22gxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th33gxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th12gxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th13gxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th23gxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th1gyy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th2gyy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th3gyy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th11gyy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th22gyy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th33gyy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th12gyy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th13gyy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th23gyy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th1gyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th2gyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th3gyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th11gyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th22gyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th33gyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th12gyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th13gyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th23gyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th1gzz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th2gzz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th3gzz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th11gzz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th22gzz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th33gzz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th12gzz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th13gzz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th23gzz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th1kxx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th2kxx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th3kxx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th1kxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th2kxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th3kxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th1kxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th2kxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th3kxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th1kyy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th2kyy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th3kyy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th1kyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th2kyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th3kyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th1kzz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th2kzz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandard4th3kzz CCTK_ATTRIBUTE_UNUSED;
    
    switch (fdOrder)
    {
      case 2:
      {
        PDstandard4th1gxx = PDstandard4th1(&gxx[index]);
        PDstandard4th2gxx = PDstandard4th2(&gxx[index]);
        PDstandard4th3gxx = PDstandard4th3(&gxx[index]);
        PDstandard4th11gxx = PDstandard4th11(&gxx[index]);
        PDstandard4th22gxx = PDstandard4th22(&gxx[index]);
        PDstandard4th33gxx = PDstandard4th33(&gxx[index]);
        PDstandard4th12gxx = PDstandard4th12(&gxx[index]);
        PDstandard4th13gxx = PDstandard4th13(&gxx[index]);
        PDstandard4th23gxx = PDstandard4th23(&gxx[index]);
        PDstandard4th1gxy = PDstandard4th1(&gxy[index]);
        PDstandard4th2gxy = PDstandard4th2(&gxy[index]);
        PDstandard4th3gxy = PDstandard4th3(&gxy[index]);
        PDstandard4th11gxy = PDstandard4th11(&gxy[index]);
        PDstandard4th22gxy = PDstandard4th22(&gxy[index]);
        PDstandard4th33gxy = PDstandard4th33(&gxy[index]);
        PDstandard4th12gxy = PDstandard4th12(&gxy[index]);
        PDstandard4th13gxy = PDstandard4th13(&gxy[index]);
        PDstandard4th23gxy = PDstandard4th23(&gxy[index]);
        PDstandard4th1gxz = PDstandard4th1(&gxz[index]);
        PDstandard4th2gxz = PDstandard4th2(&gxz[index]);
        PDstandard4th3gxz = PDstandard4th3(&gxz[index]);
        PDstandard4th11gxz = PDstandard4th11(&gxz[index]);
        PDstandard4th22gxz = PDstandard4th22(&gxz[index]);
        PDstandard4th33gxz = PDstandard4th33(&gxz[index]);
        PDstandard4th12gxz = PDstandard4th12(&gxz[index]);
        PDstandard4th13gxz = PDstandard4th13(&gxz[index]);
        PDstandard4th23gxz = PDstandard4th23(&gxz[index]);
        PDstandard4th1gyy = PDstandard4th1(&gyy[index]);
        PDstandard4th2gyy = PDstandard4th2(&gyy[index]);
        PDstandard4th3gyy = PDstandard4th3(&gyy[index]);
        PDstandard4th11gyy = PDstandard4th11(&gyy[index]);
        PDstandard4th22gyy = PDstandard4th22(&gyy[index]);
        PDstandard4th33gyy = PDstandard4th33(&gyy[index]);
        PDstandard4th12gyy = PDstandard4th12(&gyy[index]);
        PDstandard4th13gyy = PDstandard4th13(&gyy[index]);
        PDstandard4th23gyy = PDstandard4th23(&gyy[index]);
        PDstandard4th1gyz = PDstandard4th1(&gyz[index]);
        PDstandard4th2gyz = PDstandard4th2(&gyz[index]);
        PDstandard4th3gyz = PDstandard4th3(&gyz[index]);
        PDstandard4th11gyz = PDstandard4th11(&gyz[index]);
        PDstandard4th22gyz = PDstandard4th22(&gyz[index]);
        PDstandard4th33gyz = PDstandard4th33(&gyz[index]);
        PDstandard4th12gyz = PDstandard4th12(&gyz[index]);
        PDstandard4th13gyz = PDstandard4th13(&gyz[index]);
        PDstandard4th23gyz = PDstandard4th23(&gyz[index]);
        PDstandard4th1gzz = PDstandard4th1(&gzz[index]);
        PDstandard4th2gzz = PDstandard4th2(&gzz[index]);
        PDstandard4th3gzz = PDstandard4th3(&gzz[index]);
        PDstandard4th11gzz = PDstandard4th11(&gzz[index]);
        PDstandard4th22gzz = PDstandard4th22(&gzz[index]);
        PDstandard4th33gzz = PDstandard4th33(&gzz[index]);
        PDstandard4th12gzz = PDstandard4th12(&gzz[index]);
        PDstandard4th13gzz = PDstandard4th13(&gzz[index]);
        PDstandard4th23gzz = PDstandard4th23(&gzz[index]);
        PDstandard4th1kxx = PDstandard4th1(&kxx[index]);
        PDstandard4th2kxx = PDstandard4th2(&kxx[index]);
        PDstandard4th3kxx = PDstandard4th3(&kxx[index]);
        PDstandard4th1kxy = PDstandard4th1(&kxy[index]);
        PDstandard4th2kxy = PDstandard4th2(&kxy[index]);
        PDstandard4th3kxy = PDstandard4th3(&kxy[index]);
        PDstandard4th1kxz = PDstandard4th1(&kxz[index]);
        PDstandard4th2kxz = PDstandard4th2(&kxz[index]);
        PDstandard4th3kxz = PDstandard4th3(&kxz[index]);
        PDstandard4th1kyy = PDstandard4th1(&kyy[index]);
        PDstandard4th2kyy = PDstandard4th2(&kyy[index]);
        PDstandard4th3kyy = PDstandard4th3(&kyy[index]);
        PDstandard4th1kyz = PDstandard4th1(&kyz[index]);
        PDstandard4th2kyz = PDstandard4th2(&kyz[index]);
        PDstandard4th3kyz = PDstandard4th3(&kyz[index]);
        PDstandard4th1kzz = PDstandard4th1(&kzz[index]);
        PDstandard4th2kzz = PDstandard4th2(&kzz[index]);
        PDstandard4th3kzz = PDstandard4th3(&kzz[index]);
        break;
      }
      
      case 4:
      {
        PDstandard4th1gxx = PDstandard4th1(&gxx[index]);
        PDstandard4th2gxx = PDstandard4th2(&gxx[index]);
        PDstandard4th3gxx = PDstandard4th3(&gxx[index]);
        PDstandard4th11gxx = PDstandard4th11(&gxx[index]);
        PDstandard4th22gxx = PDstandard4th22(&gxx[index]);
        PDstandard4th33gxx = PDstandard4th33(&gxx[index]);
        PDstandard4th12gxx = PDstandard4th12(&gxx[index]);
        PDstandard4th13gxx = PDstandard4th13(&gxx[index]);
        PDstandard4th23gxx = PDstandard4th23(&gxx[index]);
        PDstandard4th1gxy = PDstandard4th1(&gxy[index]);
        PDstandard4th2gxy = PDstandard4th2(&gxy[index]);
        PDstandard4th3gxy = PDstandard4th3(&gxy[index]);
        PDstandard4th11gxy = PDstandard4th11(&gxy[index]);
        PDstandard4th22gxy = PDstandard4th22(&gxy[index]);
        PDstandard4th33gxy = PDstandard4th33(&gxy[index]);
        PDstandard4th12gxy = PDstandard4th12(&gxy[index]);
        PDstandard4th13gxy = PDstandard4th13(&gxy[index]);
        PDstandard4th23gxy = PDstandard4th23(&gxy[index]);
        PDstandard4th1gxz = PDstandard4th1(&gxz[index]);
        PDstandard4th2gxz = PDstandard4th2(&gxz[index]);
        PDstandard4th3gxz = PDstandard4th3(&gxz[index]);
        PDstandard4th11gxz = PDstandard4th11(&gxz[index]);
        PDstandard4th22gxz = PDstandard4th22(&gxz[index]);
        PDstandard4th33gxz = PDstandard4th33(&gxz[index]);
        PDstandard4th12gxz = PDstandard4th12(&gxz[index]);
        PDstandard4th13gxz = PDstandard4th13(&gxz[index]);
        PDstandard4th23gxz = PDstandard4th23(&gxz[index]);
        PDstandard4th1gyy = PDstandard4th1(&gyy[index]);
        PDstandard4th2gyy = PDstandard4th2(&gyy[index]);
        PDstandard4th3gyy = PDstandard4th3(&gyy[index]);
        PDstandard4th11gyy = PDstandard4th11(&gyy[index]);
        PDstandard4th22gyy = PDstandard4th22(&gyy[index]);
        PDstandard4th33gyy = PDstandard4th33(&gyy[index]);
        PDstandard4th12gyy = PDstandard4th12(&gyy[index]);
        PDstandard4th13gyy = PDstandard4th13(&gyy[index]);
        PDstandard4th23gyy = PDstandard4th23(&gyy[index]);
        PDstandard4th1gyz = PDstandard4th1(&gyz[index]);
        PDstandard4th2gyz = PDstandard4th2(&gyz[index]);
        PDstandard4th3gyz = PDstandard4th3(&gyz[index]);
        PDstandard4th11gyz = PDstandard4th11(&gyz[index]);
        PDstandard4th22gyz = PDstandard4th22(&gyz[index]);
        PDstandard4th33gyz = PDstandard4th33(&gyz[index]);
        PDstandard4th12gyz = PDstandard4th12(&gyz[index]);
        PDstandard4th13gyz = PDstandard4th13(&gyz[index]);
        PDstandard4th23gyz = PDstandard4th23(&gyz[index]);
        PDstandard4th1gzz = PDstandard4th1(&gzz[index]);
        PDstandard4th2gzz = PDstandard4th2(&gzz[index]);
        PDstandard4th3gzz = PDstandard4th3(&gzz[index]);
        PDstandard4th11gzz = PDstandard4th11(&gzz[index]);
        PDstandard4th22gzz = PDstandard4th22(&gzz[index]);
        PDstandard4th33gzz = PDstandard4th33(&gzz[index]);
        PDstandard4th12gzz = PDstandard4th12(&gzz[index]);
        PDstandard4th13gzz = PDstandard4th13(&gzz[index]);
        PDstandard4th23gzz = PDstandard4th23(&gzz[index]);
        PDstandard4th1kxx = PDstandard4th1(&kxx[index]);
        PDstandard4th2kxx = PDstandard4th2(&kxx[index]);
        PDstandard4th3kxx = PDstandard4th3(&kxx[index]);
        PDstandard4th1kxy = PDstandard4th1(&kxy[index]);
        PDstandard4th2kxy = PDstandard4th2(&kxy[index]);
        PDstandard4th3kxy = PDstandard4th3(&kxy[index]);
        PDstandard4th1kxz = PDstandard4th1(&kxz[index]);
        PDstandard4th2kxz = PDstandard4th2(&kxz[index]);
        PDstandard4th3kxz = PDstandard4th3(&kxz[index]);
        PDstandard4th1kyy = PDstandard4th1(&kyy[index]);
        PDstandard4th2kyy = PDstandard4th2(&kyy[index]);
        PDstandard4th3kyy = PDstandard4th3(&kyy[index]);
        PDstandard4th1kyz = PDstandard4th1(&kyz[index]);
        PDstandard4th2kyz = PDstandard4th2(&kyz[index]);
        PDstandard4th3kyz = PDstandard4th3(&kyz[index]);
        PDstandard4th1kzz = PDstandard4th1(&kzz[index]);
        PDstandard4th2kzz = PDstandard4th2(&kzz[index]);
        PDstandard4th3kzz = PDstandard4th3(&kzz[index]);
        break;
      }
      
      case 6:
      {
        PDstandard4th1gxx = PDstandard4th1(&gxx[index]);
        PDstandard4th2gxx = PDstandard4th2(&gxx[index]);
        PDstandard4th3gxx = PDstandard4th3(&gxx[index]);
        PDstandard4th11gxx = PDstandard4th11(&gxx[index]);
        PDstandard4th22gxx = PDstandard4th22(&gxx[index]);
        PDstandard4th33gxx = PDstandard4th33(&gxx[index]);
        PDstandard4th12gxx = PDstandard4th12(&gxx[index]);
        PDstandard4th13gxx = PDstandard4th13(&gxx[index]);
        PDstandard4th23gxx = PDstandard4th23(&gxx[index]);
        PDstandard4th1gxy = PDstandard4th1(&gxy[index]);
        PDstandard4th2gxy = PDstandard4th2(&gxy[index]);
        PDstandard4th3gxy = PDstandard4th3(&gxy[index]);
        PDstandard4th11gxy = PDstandard4th11(&gxy[index]);
        PDstandard4th22gxy = PDstandard4th22(&gxy[index]);
        PDstandard4th33gxy = PDstandard4th33(&gxy[index]);
        PDstandard4th12gxy = PDstandard4th12(&gxy[index]);
        PDstandard4th13gxy = PDstandard4th13(&gxy[index]);
        PDstandard4th23gxy = PDstandard4th23(&gxy[index]);
        PDstandard4th1gxz = PDstandard4th1(&gxz[index]);
        PDstandard4th2gxz = PDstandard4th2(&gxz[index]);
        PDstandard4th3gxz = PDstandard4th3(&gxz[index]);
        PDstandard4th11gxz = PDstandard4th11(&gxz[index]);
        PDstandard4th22gxz = PDstandard4th22(&gxz[index]);
        PDstandard4th33gxz = PDstandard4th33(&gxz[index]);
        PDstandard4th12gxz = PDstandard4th12(&gxz[index]);
        PDstandard4th13gxz = PDstandard4th13(&gxz[index]);
        PDstandard4th23gxz = PDstandard4th23(&gxz[index]);
        PDstandard4th1gyy = PDstandard4th1(&gyy[index]);
        PDstandard4th2gyy = PDstandard4th2(&gyy[index]);
        PDstandard4th3gyy = PDstandard4th3(&gyy[index]);
        PDstandard4th11gyy = PDstandard4th11(&gyy[index]);
        PDstandard4th22gyy = PDstandard4th22(&gyy[index]);
        PDstandard4th33gyy = PDstandard4th33(&gyy[index]);
        PDstandard4th12gyy = PDstandard4th12(&gyy[index]);
        PDstandard4th13gyy = PDstandard4th13(&gyy[index]);
        PDstandard4th23gyy = PDstandard4th23(&gyy[index]);
        PDstandard4th1gyz = PDstandard4th1(&gyz[index]);
        PDstandard4th2gyz = PDstandard4th2(&gyz[index]);
        PDstandard4th3gyz = PDstandard4th3(&gyz[index]);
        PDstandard4th11gyz = PDstandard4th11(&gyz[index]);
        PDstandard4th22gyz = PDstandard4th22(&gyz[index]);
        PDstandard4th33gyz = PDstandard4th33(&gyz[index]);
        PDstandard4th12gyz = PDstandard4th12(&gyz[index]);
        PDstandard4th13gyz = PDstandard4th13(&gyz[index]);
        PDstandard4th23gyz = PDstandard4th23(&gyz[index]);
        PDstandard4th1gzz = PDstandard4th1(&gzz[index]);
        PDstandard4th2gzz = PDstandard4th2(&gzz[index]);
        PDstandard4th3gzz = PDstandard4th3(&gzz[index]);
        PDstandard4th11gzz = PDstandard4th11(&gzz[index]);
        PDstandard4th22gzz = PDstandard4th22(&gzz[index]);
        PDstandard4th33gzz = PDstandard4th33(&gzz[index]);
        PDstandard4th12gzz = PDstandard4th12(&gzz[index]);
        PDstandard4th13gzz = PDstandard4th13(&gzz[index]);
        PDstandard4th23gzz = PDstandard4th23(&gzz[index]);
        PDstandard4th1kxx = PDstandard4th1(&kxx[index]);
        PDstandard4th2kxx = PDstandard4th2(&kxx[index]);
        PDstandard4th3kxx = PDstandard4th3(&kxx[index]);
        PDstandard4th1kxy = PDstandard4th1(&kxy[index]);
        PDstandard4th2kxy = PDstandard4th2(&kxy[index]);
        PDstandard4th3kxy = PDstandard4th3(&kxy[index]);
        PDstandard4th1kxz = PDstandard4th1(&kxz[index]);
        PDstandard4th2kxz = PDstandard4th2(&kxz[index]);
        PDstandard4th3kxz = PDstandard4th3(&kxz[index]);
        PDstandard4th1kyy = PDstandard4th1(&kyy[index]);
        PDstandard4th2kyy = PDstandard4th2(&kyy[index]);
        PDstandard4th3kyy = PDstandard4th3(&kyy[index]);
        PDstandard4th1kyz = PDstandard4th1(&kyz[index]);
        PDstandard4th2kyz = PDstandard4th2(&kyz[index]);
        PDstandard4th3kyz = PDstandard4th3(&kyz[index]);
        PDstandard4th1kzz = PDstandard4th1(&kzz[index]);
        PDstandard4th2kzz = PDstandard4th2(&kzz[index]);
        PDstandard4th3kzz = PDstandard4th3(&kzz[index]);
        break;
      }
      
      case 8:
      {
        PDstandard4th1gxx = PDstandard4th1(&gxx[index]);
        PDstandard4th2gxx = PDstandard4th2(&gxx[index]);
        PDstandard4th3gxx = PDstandard4th3(&gxx[index]);
        PDstandard4th11gxx = PDstandard4th11(&gxx[index]);
        PDstandard4th22gxx = PDstandard4th22(&gxx[index]);
        PDstandard4th33gxx = PDstandard4th33(&gxx[index]);
        PDstandard4th12gxx = PDstandard4th12(&gxx[index]);
        PDstandard4th13gxx = PDstandard4th13(&gxx[index]);
        PDstandard4th23gxx = PDstandard4th23(&gxx[index]);
        PDstandard4th1gxy = PDstandard4th1(&gxy[index]);
        PDstandard4th2gxy = PDstandard4th2(&gxy[index]);
        PDstandard4th3gxy = PDstandard4th3(&gxy[index]);
        PDstandard4th11gxy = PDstandard4th11(&gxy[index]);
        PDstandard4th22gxy = PDstandard4th22(&gxy[index]);
        PDstandard4th33gxy = PDstandard4th33(&gxy[index]);
        PDstandard4th12gxy = PDstandard4th12(&gxy[index]);
        PDstandard4th13gxy = PDstandard4th13(&gxy[index]);
        PDstandard4th23gxy = PDstandard4th23(&gxy[index]);
        PDstandard4th1gxz = PDstandard4th1(&gxz[index]);
        PDstandard4th2gxz = PDstandard4th2(&gxz[index]);
        PDstandard4th3gxz = PDstandard4th3(&gxz[index]);
        PDstandard4th11gxz = PDstandard4th11(&gxz[index]);
        PDstandard4th22gxz = PDstandard4th22(&gxz[index]);
        PDstandard4th33gxz = PDstandard4th33(&gxz[index]);
        PDstandard4th12gxz = PDstandard4th12(&gxz[index]);
        PDstandard4th13gxz = PDstandard4th13(&gxz[index]);
        PDstandard4th23gxz = PDstandard4th23(&gxz[index]);
        PDstandard4th1gyy = PDstandard4th1(&gyy[index]);
        PDstandard4th2gyy = PDstandard4th2(&gyy[index]);
        PDstandard4th3gyy = PDstandard4th3(&gyy[index]);
        PDstandard4th11gyy = PDstandard4th11(&gyy[index]);
        PDstandard4th22gyy = PDstandard4th22(&gyy[index]);
        PDstandard4th33gyy = PDstandard4th33(&gyy[index]);
        PDstandard4th12gyy = PDstandard4th12(&gyy[index]);
        PDstandard4th13gyy = PDstandard4th13(&gyy[index]);
        PDstandard4th23gyy = PDstandard4th23(&gyy[index]);
        PDstandard4th1gyz = PDstandard4th1(&gyz[index]);
        PDstandard4th2gyz = PDstandard4th2(&gyz[index]);
        PDstandard4th3gyz = PDstandard4th3(&gyz[index]);
        PDstandard4th11gyz = PDstandard4th11(&gyz[index]);
        PDstandard4th22gyz = PDstandard4th22(&gyz[index]);
        PDstandard4th33gyz = PDstandard4th33(&gyz[index]);
        PDstandard4th12gyz = PDstandard4th12(&gyz[index]);
        PDstandard4th13gyz = PDstandard4th13(&gyz[index]);
        PDstandard4th23gyz = PDstandard4th23(&gyz[index]);
        PDstandard4th1gzz = PDstandard4th1(&gzz[index]);
        PDstandard4th2gzz = PDstandard4th2(&gzz[index]);
        PDstandard4th3gzz = PDstandard4th3(&gzz[index]);
        PDstandard4th11gzz = PDstandard4th11(&gzz[index]);
        PDstandard4th22gzz = PDstandard4th22(&gzz[index]);
        PDstandard4th33gzz = PDstandard4th33(&gzz[index]);
        PDstandard4th12gzz = PDstandard4th12(&gzz[index]);
        PDstandard4th13gzz = PDstandard4th13(&gzz[index]);
        PDstandard4th23gzz = PDstandard4th23(&gzz[index]);
        PDstandard4th1kxx = PDstandard4th1(&kxx[index]);
        PDstandard4th2kxx = PDstandard4th2(&kxx[index]);
        PDstandard4th3kxx = PDstandard4th3(&kxx[index]);
        PDstandard4th1kxy = PDstandard4th1(&kxy[index]);
        PDstandard4th2kxy = PDstandard4th2(&kxy[index]);
        PDstandard4th3kxy = PDstandard4th3(&kxy[index]);
        PDstandard4th1kxz = PDstandard4th1(&kxz[index]);
        PDstandard4th2kxz = PDstandard4th2(&kxz[index]);
        PDstandard4th3kxz = PDstandard4th3(&kxz[index]);
        PDstandard4th1kyy = PDstandard4th1(&kyy[index]);
        PDstandard4th2kyy = PDstandard4th2(&kyy[index]);
        PDstandard4th3kyy = PDstandard4th3(&kyy[index]);
        PDstandard4th1kyz = PDstandard4th1(&kyz[index]);
        PDstandard4th2kyz = PDstandard4th2(&kyz[index]);
        PDstandard4th3kyz = PDstandard4th3(&kyz[index]);
        PDstandard4th1kzz = PDstandard4th1(&kzz[index]);
        PDstandard4th2kzz = PDstandard4th2(&kzz[index]);
        PDstandard4th3kzz = PDstandard4th3(&kzz[index]);
        break;
      }
      default:
        CCTK_BUILTIN_UNREACHABLE();
    }
    /* Calculate temporaries and grid functions */
    CCTK_REAL_VEC JacPDstandard4th11gyy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandard4th11gyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandard4th11gzz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandard4th12gxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandard4th12gxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandard4th12gyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandard4th12gzz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandard4th13gxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandard4th1gxx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandard4th1gxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandard4th1gxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandard4th1gyy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandard4th1gyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandard4th1gzz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandard4th1kxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandard4th1kxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandard4th1kyy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandard4th1kyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandard4th1kzz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandard4th21gxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandard4th22gxx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandard4th22gxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandard4th22gzz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandard4th23gxx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandard4th23gxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandard4th23gxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandard4th23gyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandard4th2gxx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandard4th2gxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandard4th2gxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandard4th2gyy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandard4th2gyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandard4th2gzz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandard4th2kxx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandard4th2kxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandard4th2kxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandard4th2kyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandard4th2kzz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandard4th31gxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandard4th31gxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandard4th31gyy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandard4th31gyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandard4th32gyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandard4th33gxx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandard4th33gxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandard4th33gyy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandard4th3gxx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandard4th3gxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandard4th3gxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandard4th3gyy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandard4th3gyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandard4th3gzz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandard4th3kxx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandard4th3kxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandard4th3kxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandard4th3kyy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandard4th3kyz CCTK_ATTRIBUTE_UNUSED;
    
    if (usejacobian)
    {
      JacPDstandard4th1gxx = 
        kmadd(J11L,PDstandard4th1gxx,kmadd(J21L,PDstandard4th2gxx,kmul(J31L,PDstandard4th3gxx)));
      
      JacPDstandard4th1gxy = 
        kmadd(J11L,PDstandard4th1gxy,kmadd(J21L,PDstandard4th2gxy,kmul(J31L,PDstandard4th3gxy)));
      
      JacPDstandard4th1gxz = 
        kmadd(J11L,PDstandard4th1gxz,kmadd(J21L,PDstandard4th2gxz,kmul(J31L,PDstandard4th3gxz)));
      
      JacPDstandard4th1gyy = 
        kmadd(J11L,PDstandard4th1gyy,kmadd(J21L,PDstandard4th2gyy,kmul(J31L,PDstandard4th3gyy)));
      
      JacPDstandard4th1gyz = 
        kmadd(J11L,PDstandard4th1gyz,kmadd(J21L,PDstandard4th2gyz,kmul(J31L,PDstandard4th3gyz)));
      
      JacPDstandard4th1gzz = 
        kmadd(J11L,PDstandard4th1gzz,kmadd(J21L,PDstandard4th2gzz,kmul(J31L,PDstandard4th3gzz)));
      
      JacPDstandard4th1kxy = 
        kmadd(J11L,PDstandard4th1kxy,kmadd(J21L,PDstandard4th2kxy,kmul(J31L,PDstandard4th3kxy)));
      
      JacPDstandard4th1kxz = 
        kmadd(J11L,PDstandard4th1kxz,kmadd(J21L,PDstandard4th2kxz,kmul(J31L,PDstandard4th3kxz)));
      
      JacPDstandard4th1kyy = 
        kmadd(J11L,PDstandard4th1kyy,kmadd(J21L,PDstandard4th2kyy,kmul(J31L,PDstandard4th3kyy)));
      
      JacPDstandard4th1kyz = 
        kmadd(J11L,PDstandard4th1kyz,kmadd(J21L,PDstandard4th2kyz,kmul(J31L,PDstandard4th3kyz)));
      
      JacPDstandard4th1kzz = 
        kmadd(J11L,PDstandard4th1kzz,kmadd(J21L,PDstandard4th2kzz,kmul(J31L,PDstandard4th3kzz)));
      
      JacPDstandard4th2gxx = 
        kmadd(J12L,PDstandard4th1gxx,kmadd(J22L,PDstandard4th2gxx,kmul(J32L,PDstandard4th3gxx)));
      
      JacPDstandard4th2gxy = 
        kmadd(J12L,PDstandard4th1gxy,kmadd(J22L,PDstandard4th2gxy,kmul(J32L,PDstandard4th3gxy)));
      
      JacPDstandard4th2gxz = 
        kmadd(J12L,PDstandard4th1gxz,kmadd(J22L,PDstandard4th2gxz,kmul(J32L,PDstandard4th3gxz)));
      
      JacPDstandard4th2gyy = 
        kmadd(J12L,PDstandard4th1gyy,kmadd(J22L,PDstandard4th2gyy,kmul(J32L,PDstandard4th3gyy)));
      
      JacPDstandard4th2gyz = 
        kmadd(J12L,PDstandard4th1gyz,kmadd(J22L,PDstandard4th2gyz,kmul(J32L,PDstandard4th3gyz)));
      
      JacPDstandard4th2gzz = 
        kmadd(J12L,PDstandard4th1gzz,kmadd(J22L,PDstandard4th2gzz,kmul(J32L,PDstandard4th3gzz)));
      
      JacPDstandard4th2kxx = 
        kmadd(J12L,PDstandard4th1kxx,kmadd(J22L,PDstandard4th2kxx,kmul(J32L,PDstandard4th3kxx)));
      
      JacPDstandard4th2kxy = 
        kmadd(J12L,PDstandard4th1kxy,kmadd(J22L,PDstandard4th2kxy,kmul(J32L,PDstandard4th3kxy)));
      
      JacPDstandard4th2kxz = 
        kmadd(J12L,PDstandard4th1kxz,kmadd(J22L,PDstandard4th2kxz,kmul(J32L,PDstandard4th3kxz)));
      
      JacPDstandard4th2kyz = 
        kmadd(J12L,PDstandard4th1kyz,kmadd(J22L,PDstandard4th2kyz,kmul(J32L,PDstandard4th3kyz)));
      
      JacPDstandard4th2kzz = 
        kmadd(J12L,PDstandard4th1kzz,kmadd(J22L,PDstandard4th2kzz,kmul(J32L,PDstandard4th3kzz)));
      
      JacPDstandard4th3gxx = 
        kmadd(J13L,PDstandard4th1gxx,kmadd(J23L,PDstandard4th2gxx,kmul(J33L,PDstandard4th3gxx)));
      
      JacPDstandard4th3gxy = 
        kmadd(J13L,PDstandard4th1gxy,kmadd(J23L,PDstandard4th2gxy,kmul(J33L,PDstandard4th3gxy)));
      
      JacPDstandard4th3gxz = 
        kmadd(J13L,PDstandard4th1gxz,kmadd(J23L,PDstandard4th2gxz,kmul(J33L,PDstandard4th3gxz)));
      
      JacPDstandard4th3gyy = 
        kmadd(J13L,PDstandard4th1gyy,kmadd(J23L,PDstandard4th2gyy,kmul(J33L,PDstandard4th3gyy)));
      
      JacPDstandard4th3gyz = 
        kmadd(J13L,PDstandard4th1gyz,kmadd(J23L,PDstandard4th2gyz,kmul(J33L,PDstandard4th3gyz)));
      
      JacPDstandard4th3gzz = 
        kmadd(J13L,PDstandard4th1gzz,kmadd(J23L,PDstandard4th2gzz,kmul(J33L,PDstandard4th3gzz)));
      
      JacPDstandard4th3kxx = 
        kmadd(J13L,PDstandard4th1kxx,kmadd(J23L,PDstandard4th2kxx,kmul(J33L,PDstandard4th3kxx)));
      
      JacPDstandard4th3kxy = 
        kmadd(J13L,PDstandard4th1kxy,kmadd(J23L,PDstandard4th2kxy,kmul(J33L,PDstandard4th3kxy)));
      
      JacPDstandard4th3kxz = 
        kmadd(J13L,PDstandard4th1kxz,kmadd(J23L,PDstandard4th2kxz,kmul(J33L,PDstandard4th3kxz)));
      
      JacPDstandard4th3kyy = 
        kmadd(J13L,PDstandard4th1kyy,kmadd(J23L,PDstandard4th2kyy,kmul(J33L,PDstandard4th3kyy)));
      
      JacPDstandard4th3kyz = 
        kmadd(J13L,PDstandard4th1kyz,kmadd(J23L,PDstandard4th2kyz,kmul(J33L,PDstandard4th3kyz)));
      
      JacPDstandard4th11gyy = 
        kmadd(dJ111L,PDstandard4th1gyy,kmadd(ToReal(2),kmadd(J11L,kmadd(J21L,PDstandard4th12gyy,kmul(J31L,PDstandard4th13gyy)),kmul(J21L,kmul(J31L,PDstandard4th23gyy))),kmadd(dJ211L,PDstandard4th2gyy,kmadd(dJ311L,PDstandard4th3gyy,kmadd(PDstandard4th11gyy,kmul(J11L,J11L),kmadd(PDstandard4th22gyy,kmul(J21L,J21L),kmul(PDstandard4th33gyy,kmul(J31L,J31L))))))));
      
      JacPDstandard4th11gyz = 
        kmadd(dJ111L,PDstandard4th1gyz,kmadd(ToReal(2),kmadd(J11L,kmadd(J21L,PDstandard4th12gyz,kmul(J31L,PDstandard4th13gyz)),kmul(J21L,kmul(J31L,PDstandard4th23gyz))),kmadd(dJ211L,PDstandard4th2gyz,kmadd(dJ311L,PDstandard4th3gyz,kmadd(PDstandard4th11gyz,kmul(J11L,J11L),kmadd(PDstandard4th22gyz,kmul(J21L,J21L),kmul(PDstandard4th33gyz,kmul(J31L,J31L))))))));
      
      JacPDstandard4th11gzz = 
        kmadd(dJ111L,PDstandard4th1gzz,kmadd(ToReal(2),kmadd(J11L,kmadd(J21L,PDstandard4th12gzz,kmul(J31L,PDstandard4th13gzz)),kmul(J21L,kmul(J31L,PDstandard4th23gzz))),kmadd(dJ211L,PDstandard4th2gzz,kmadd(dJ311L,PDstandard4th3gzz,kmadd(PDstandard4th11gzz,kmul(J11L,J11L),kmadd(PDstandard4th22gzz,kmul(J21L,J21L),kmul(PDstandard4th33gzz,kmul(J31L,J31L))))))));
      
      JacPDstandard4th22gxx = 
        kmadd(dJ122L,PDstandard4th1gxx,kmadd(ToReal(2),kmadd(J12L,kmadd(J22L,PDstandard4th12gxx,kmul(J32L,PDstandard4th13gxx)),kmul(J22L,kmul(J32L,PDstandard4th23gxx))),kmadd(dJ222L,PDstandard4th2gxx,kmadd(dJ322L,PDstandard4th3gxx,kmadd(PDstandard4th11gxx,kmul(J12L,J12L),kmadd(PDstandard4th22gxx,kmul(J22L,J22L),kmul(PDstandard4th33gxx,kmul(J32L,J32L))))))));
      
      JacPDstandard4th22gxz = 
        kmadd(dJ122L,PDstandard4th1gxz,kmadd(ToReal(2),kmadd(J12L,kmadd(J22L,PDstandard4th12gxz,kmul(J32L,PDstandard4th13gxz)),kmul(J22L,kmul(J32L,PDstandard4th23gxz))),kmadd(dJ222L,PDstandard4th2gxz,kmadd(dJ322L,PDstandard4th3gxz,kmadd(PDstandard4th11gxz,kmul(J12L,J12L),kmadd(PDstandard4th22gxz,kmul(J22L,J22L),kmul(PDstandard4th33gxz,kmul(J32L,J32L))))))));
      
      JacPDstandard4th22gzz = 
        kmadd(dJ122L,PDstandard4th1gzz,kmadd(ToReal(2),kmadd(J12L,kmadd(J22L,PDstandard4th12gzz,kmul(J32L,PDstandard4th13gzz)),kmul(J22L,kmul(J32L,PDstandard4th23gzz))),kmadd(dJ222L,PDstandard4th2gzz,kmadd(dJ322L,PDstandard4th3gzz,kmadd(PDstandard4th11gzz,kmul(J12L,J12L),kmadd(PDstandard4th22gzz,kmul(J22L,J22L),kmul(PDstandard4th33gzz,kmul(J32L,J32L))))))));
      
      JacPDstandard4th33gxx = 
        kmadd(dJ133L,PDstandard4th1gxx,kmadd(ToReal(2),kmadd(J13L,kmadd(J23L,PDstandard4th12gxx,kmul(J33L,PDstandard4th13gxx)),kmul(J23L,kmul(J33L,PDstandard4th23gxx))),kmadd(dJ233L,PDstandard4th2gxx,kmadd(dJ333L,PDstandard4th3gxx,kmadd(PDstandard4th11gxx,kmul(J13L,J13L),kmadd(PDstandard4th22gxx,kmul(J23L,J23L),kmul(PDstandard4th33gxx,kmul(J33L,J33L))))))));
      
      JacPDstandard4th33gxy = 
        kmadd(dJ133L,PDstandard4th1gxy,kmadd(ToReal(2),kmadd(J13L,kmadd(J23L,PDstandard4th12gxy,kmul(J33L,PDstandard4th13gxy)),kmul(J23L,kmul(J33L,PDstandard4th23gxy))),kmadd(dJ233L,PDstandard4th2gxy,kmadd(dJ333L,PDstandard4th3gxy,kmadd(PDstandard4th11gxy,kmul(J13L,J13L),kmadd(PDstandard4th22gxy,kmul(J23L,J23L),kmul(PDstandard4th33gxy,kmul(J33L,J33L))))))));
      
      JacPDstandard4th33gyy = 
        kmadd(dJ133L,PDstandard4th1gyy,kmadd(ToReal(2),kmadd(J13L,kmadd(J23L,PDstandard4th12gyy,kmul(J33L,PDstandard4th13gyy)),kmul(J23L,kmul(J33L,PDstandard4th23gyy))),kmadd(dJ233L,PDstandard4th2gyy,kmadd(dJ333L,PDstandard4th3gyy,kmadd(PDstandard4th11gyy,kmul(J13L,J13L),kmadd(PDstandard4th22gyy,kmul(J23L,J23L),kmul(PDstandard4th33gyy,kmul(J33L,J33L))))))));
      
      JacPDstandard4th12gxy = 
        kmadd(J12L,kmadd(J11L,PDstandard4th11gxy,kmadd(J21L,PDstandard4th12gxy,kmul(J31L,PDstandard4th13gxy))),kmadd(J11L,kmadd(J22L,PDstandard4th12gxy,kmul(J32L,PDstandard4th13gxy)),kmadd(dJ112L,PDstandard4th1gxy,kmadd(J22L,kmadd(J21L,PDstandard4th22gxy,kmul(J31L,PDstandard4th23gxy)),kmadd(dJ212L,PDstandard4th2gxy,kmadd(J32L,kmadd(J21L,PDstandard4th23gxy,kmul(J31L,PDstandard4th33gxy)),kmul(dJ312L,PDstandard4th3gxy)))))));
      
      JacPDstandard4th12gxz = 
        kmadd(J12L,kmadd(J11L,PDstandard4th11gxz,kmadd(J21L,PDstandard4th12gxz,kmul(J31L,PDstandard4th13gxz))),kmadd(J11L,kmadd(J22L,PDstandard4th12gxz,kmul(J32L,PDstandard4th13gxz)),kmadd(dJ112L,PDstandard4th1gxz,kmadd(J22L,kmadd(J21L,PDstandard4th22gxz,kmul(J31L,PDstandard4th23gxz)),kmadd(dJ212L,PDstandard4th2gxz,kmadd(J32L,kmadd(J21L,PDstandard4th23gxz,kmul(J31L,PDstandard4th33gxz)),kmul(dJ312L,PDstandard4th3gxz)))))));
      
      JacPDstandard4th12gyz = 
        kmadd(J12L,kmadd(J11L,PDstandard4th11gyz,kmadd(J21L,PDstandard4th12gyz,kmul(J31L,PDstandard4th13gyz))),kmadd(J11L,kmadd(J22L,PDstandard4th12gyz,kmul(J32L,PDstandard4th13gyz)),kmadd(dJ112L,PDstandard4th1gyz,kmadd(J22L,kmadd(J21L,PDstandard4th22gyz,kmul(J31L,PDstandard4th23gyz)),kmadd(dJ212L,PDstandard4th2gyz,kmadd(J32L,kmadd(J21L,PDstandard4th23gyz,kmul(J31L,PDstandard4th33gyz)),kmul(dJ312L,PDstandard4th3gyz)))))));
      
      JacPDstandard4th12gzz = 
        kmadd(J12L,kmadd(J11L,PDstandard4th11gzz,kmadd(J21L,PDstandard4th12gzz,kmul(J31L,PDstandard4th13gzz))),kmadd(J11L,kmadd(J22L,PDstandard4th12gzz,kmul(J32L,PDstandard4th13gzz)),kmadd(dJ112L,PDstandard4th1gzz,kmadd(J22L,kmadd(J21L,PDstandard4th22gzz,kmul(J31L,PDstandard4th23gzz)),kmadd(dJ212L,PDstandard4th2gzz,kmadd(J32L,kmadd(J21L,PDstandard4th23gzz,kmul(J31L,PDstandard4th33gzz)),kmul(dJ312L,PDstandard4th3gzz)))))));
      
      JacPDstandard4th13gxz = 
        kmadd(J13L,kmadd(J11L,PDstandard4th11gxz,kmadd(J21L,PDstandard4th12gxz,kmul(J31L,PDstandard4th13gxz))),kmadd(J11L,kmadd(J23L,PDstandard4th12gxz,kmul(J33L,PDstandard4th13gxz)),kmadd(dJ113L,PDstandard4th1gxz,kmadd(J23L,kmadd(J21L,PDstandard4th22gxz,kmul(J31L,PDstandard4th23gxz)),kmadd(dJ213L,PDstandard4th2gxz,kmadd(J33L,kmadd(J21L,PDstandard4th23gxz,kmul(J31L,PDstandard4th33gxz)),kmul(dJ313L,PDstandard4th3gxz)))))));
      
      JacPDstandard4th21gxy = 
        kmadd(J12L,kmadd(J11L,PDstandard4th11gxy,kmadd(J21L,PDstandard4th12gxy,kmul(J31L,PDstandard4th13gxy))),kmadd(J11L,kmadd(J22L,PDstandard4th12gxy,kmul(J32L,PDstandard4th13gxy)),kmadd(dJ112L,PDstandard4th1gxy,kmadd(J22L,kmadd(J21L,PDstandard4th22gxy,kmul(J31L,PDstandard4th23gxy)),kmadd(dJ212L,PDstandard4th2gxy,kmadd(J32L,kmadd(J21L,PDstandard4th23gxy,kmul(J31L,PDstandard4th33gxy)),kmul(dJ312L,PDstandard4th3gxy)))))));
      
      JacPDstandard4th23gxx = 
        kmadd(J13L,kmadd(J12L,PDstandard4th11gxx,kmadd(J22L,PDstandard4th12gxx,kmul(J32L,PDstandard4th13gxx))),kmadd(J12L,kmadd(J23L,PDstandard4th12gxx,kmul(J33L,PDstandard4th13gxx)),kmadd(dJ123L,PDstandard4th1gxx,kmadd(J23L,kmadd(J22L,PDstandard4th22gxx,kmul(J32L,PDstandard4th23gxx)),kmadd(dJ223L,PDstandard4th2gxx,kmadd(J33L,kmadd(J22L,PDstandard4th23gxx,kmul(J32L,PDstandard4th33gxx)),kmul(dJ323L,PDstandard4th3gxx)))))));
      
      JacPDstandard4th23gxy = 
        kmadd(J13L,kmadd(J12L,PDstandard4th11gxy,kmadd(J22L,PDstandard4th12gxy,kmul(J32L,PDstandard4th13gxy))),kmadd(J12L,kmadd(J23L,PDstandard4th12gxy,kmul(J33L,PDstandard4th13gxy)),kmadd(dJ123L,PDstandard4th1gxy,kmadd(J23L,kmadd(J22L,PDstandard4th22gxy,kmul(J32L,PDstandard4th23gxy)),kmadd(dJ223L,PDstandard4th2gxy,kmadd(J33L,kmadd(J22L,PDstandard4th23gxy,kmul(J32L,PDstandard4th33gxy)),kmul(dJ323L,PDstandard4th3gxy)))))));
      
      JacPDstandard4th23gxz = 
        kmadd(J13L,kmadd(J12L,PDstandard4th11gxz,kmadd(J22L,PDstandard4th12gxz,kmul(J32L,PDstandard4th13gxz))),kmadd(J12L,kmadd(J23L,PDstandard4th12gxz,kmul(J33L,PDstandard4th13gxz)),kmadd(dJ123L,PDstandard4th1gxz,kmadd(J23L,kmadd(J22L,PDstandard4th22gxz,kmul(J32L,PDstandard4th23gxz)),kmadd(dJ223L,PDstandard4th2gxz,kmadd(J33L,kmadd(J22L,PDstandard4th23gxz,kmul(J32L,PDstandard4th33gxz)),kmul(dJ323L,PDstandard4th3gxz)))))));
      
      JacPDstandard4th23gyz = 
        kmadd(J13L,kmadd(J12L,PDstandard4th11gyz,kmadd(J22L,PDstandard4th12gyz,kmul(J32L,PDstandard4th13gyz))),kmadd(J12L,kmadd(J23L,PDstandard4th12gyz,kmul(J33L,PDstandard4th13gyz)),kmadd(dJ123L,PDstandard4th1gyz,kmadd(J23L,kmadd(J22L,PDstandard4th22gyz,kmul(J32L,PDstandard4th23gyz)),kmadd(dJ223L,PDstandard4th2gyz,kmadd(J33L,kmadd(J22L,PDstandard4th23gyz,kmul(J32L,PDstandard4th33gyz)),kmul(dJ323L,PDstandard4th3gyz)))))));
      
      JacPDstandard4th31gxy = 
        kmadd(J13L,kmadd(J11L,PDstandard4th11gxy,kmadd(J21L,PDstandard4th12gxy,kmul(J31L,PDstandard4th13gxy))),kmadd(J11L,kmadd(J23L,PDstandard4th12gxy,kmul(J33L,PDstandard4th13gxy)),kmadd(dJ113L,PDstandard4th1gxy,kmadd(J23L,kmadd(J21L,PDstandard4th22gxy,kmul(J31L,PDstandard4th23gxy)),kmadd(dJ213L,PDstandard4th2gxy,kmadd(J33L,kmadd(J21L,PDstandard4th23gxy,kmul(J31L,PDstandard4th33gxy)),kmul(dJ313L,PDstandard4th3gxy)))))));
      
      JacPDstandard4th31gxz = 
        kmadd(J13L,kmadd(J11L,PDstandard4th11gxz,kmadd(J21L,PDstandard4th12gxz,kmul(J31L,PDstandard4th13gxz))),kmadd(J11L,kmadd(J23L,PDstandard4th12gxz,kmul(J33L,PDstandard4th13gxz)),kmadd(dJ113L,PDstandard4th1gxz,kmadd(J23L,kmadd(J21L,PDstandard4th22gxz,kmul(J31L,PDstandard4th23gxz)),kmadd(dJ213L,PDstandard4th2gxz,kmadd(J33L,kmadd(J21L,PDstandard4th23gxz,kmul(J31L,PDstandard4th33gxz)),kmul(dJ313L,PDstandard4th3gxz)))))));
      
      JacPDstandard4th31gyy = 
        kmadd(J13L,kmadd(J11L,PDstandard4th11gyy,kmadd(J21L,PDstandard4th12gyy,kmul(J31L,PDstandard4th13gyy))),kmadd(J11L,kmadd(J23L,PDstandard4th12gyy,kmul(J33L,PDstandard4th13gyy)),kmadd(dJ113L,PDstandard4th1gyy,kmadd(J23L,kmadd(J21L,PDstandard4th22gyy,kmul(J31L,PDstandard4th23gyy)),kmadd(dJ213L,PDstandard4th2gyy,kmadd(J33L,kmadd(J21L,PDstandard4th23gyy,kmul(J31L,PDstandard4th33gyy)),kmul(dJ313L,PDstandard4th3gyy)))))));
      
      JacPDstandard4th31gyz = 
        kmadd(J13L,kmadd(J11L,PDstandard4th11gyz,kmadd(J21L,PDstandard4th12gyz,kmul(J31L,PDstandard4th13gyz))),kmadd(J11L,kmadd(J23L,PDstandard4th12gyz,kmul(J33L,PDstandard4th13gyz)),kmadd(dJ113L,PDstandard4th1gyz,kmadd(J23L,kmadd(J21L,PDstandard4th22gyz,kmul(J31L,PDstandard4th23gyz)),kmadd(dJ213L,PDstandard4th2gyz,kmadd(J33L,kmadd(J21L,PDstandard4th23gyz,kmul(J31L,PDstandard4th33gyz)),kmul(dJ313L,PDstandard4th3gyz)))))));
      
      JacPDstandard4th32gyz = 
        kmadd(J13L,kmadd(J12L,PDstandard4th11gyz,kmadd(J22L,PDstandard4th12gyz,kmul(J32L,PDstandard4th13gyz))),kmadd(J12L,kmadd(J23L,PDstandard4th12gyz,kmul(J33L,PDstandard4th13gyz)),kmadd(dJ123L,PDstandard4th1gyz,kmadd(J23L,kmadd(J22L,PDstandard4th22gyz,kmul(J32L,PDstandard4th23gyz)),kmadd(dJ223L,PDstandard4th2gyz,kmadd(J33L,kmadd(J22L,PDstandard4th23gyz,kmul(J32L,PDstandard4th33gyz)),kmul(dJ323L,PDstandard4th3gyz)))))));
    }
    else
    {
      JacPDstandard4th1gxx = PDstandard4th1gxx;
      
      JacPDstandard4th1gxy = PDstandard4th1gxy;
      
      JacPDstandard4th1gxz = PDstandard4th1gxz;
      
      JacPDstandard4th1gyy = PDstandard4th1gyy;
      
      JacPDstandard4th1gyz = PDstandard4th1gyz;
      
      JacPDstandard4th1gzz = PDstandard4th1gzz;
      
      JacPDstandard4th1kxy = PDstandard4th1kxy;
      
      JacPDstandard4th1kxz = PDstandard4th1kxz;
      
      JacPDstandard4th1kyy = PDstandard4th1kyy;
      
      JacPDstandard4th1kyz = PDstandard4th1kyz;
      
      JacPDstandard4th1kzz = PDstandard4th1kzz;
      
      JacPDstandard4th2gxx = PDstandard4th2gxx;
      
      JacPDstandard4th2gxy = PDstandard4th2gxy;
      
      JacPDstandard4th2gxz = PDstandard4th2gxz;
      
      JacPDstandard4th2gyy = PDstandard4th2gyy;
      
      JacPDstandard4th2gyz = PDstandard4th2gyz;
      
      JacPDstandard4th2gzz = PDstandard4th2gzz;
      
      JacPDstandard4th2kxx = PDstandard4th2kxx;
      
      JacPDstandard4th2kxy = PDstandard4th2kxy;
      
      JacPDstandard4th2kxz = PDstandard4th2kxz;
      
      JacPDstandard4th2kyz = PDstandard4th2kyz;
      
      JacPDstandard4th2kzz = PDstandard4th2kzz;
      
      JacPDstandard4th3gxx = PDstandard4th3gxx;
      
      JacPDstandard4th3gxy = PDstandard4th3gxy;
      
      JacPDstandard4th3gxz = PDstandard4th3gxz;
      
      JacPDstandard4th3gyy = PDstandard4th3gyy;
      
      JacPDstandard4th3gyz = PDstandard4th3gyz;
      
      JacPDstandard4th3gzz = PDstandard4th3gzz;
      
      JacPDstandard4th3kxx = PDstandard4th3kxx;
      
      JacPDstandard4th3kxy = PDstandard4th3kxy;
      
      JacPDstandard4th3kxz = PDstandard4th3kxz;
      
      JacPDstandard4th3kyy = PDstandard4th3kyy;
      
      JacPDstandard4th3kyz = PDstandard4th3kyz;
      
      JacPDstandard4th11gyy = PDstandard4th11gyy;
      
      JacPDstandard4th11gyz = PDstandard4th11gyz;
      
      JacPDstandard4th11gzz = PDstandard4th11gzz;
      
      JacPDstandard4th22gxx = PDstandard4th22gxx;
      
      JacPDstandard4th22gxz = PDstandard4th22gxz;
      
      JacPDstandard4th22gzz = PDstandard4th22gzz;
      
      JacPDstandard4th33gxx = PDstandard4th33gxx;
      
      JacPDstandard4th33gxy = PDstandard4th33gxy;
      
      JacPDstandard4th33gyy = PDstandard4th33gyy;
      
      JacPDstandard4th12gxy = PDstandard4th12gxy;
      
      JacPDstandard4th12gxz = PDstandard4th12gxz;
      
      JacPDstandard4th12gyz = PDstandard4th12gyz;
      
      JacPDstandard4th12gzz = PDstandard4th12gzz;
      
      JacPDstandard4th13gxz = PDstandard4th13gxz;
      
      JacPDstandard4th21gxy = PDstandard4th12gxy;
      
      JacPDstandard4th23gxx = PDstandard4th23gxx;
      
      JacPDstandard4th23gxy = PDstandard4th23gxy;
      
      JacPDstandard4th23gxz = PDstandard4th23gxz;
      
      JacPDstandard4th23gyz = PDstandard4th23gyz;
      
      JacPDstandard4th31gxy = PDstandard4th13gxy;
      
      JacPDstandard4th31gxz = PDstandard4th13gxz;
      
      JacPDstandard4th31gyy = PDstandard4th13gyy;
      
      JacPDstandard4th31gyz = PDstandard4th13gyz;
      
      JacPDstandard4th32gyz = PDstandard4th23gyz;
    }
    
    CCTK_REAL_VEC detg CCTK_ATTRIBUTE_UNUSED = 
      kmadd(ToReal(2),kmul(gxyL,kmul(gxzL,gyzL)),knmsub(gzzL,kmul(gxyL,gxyL),kmsub(gyyL,kmsub(gxxL,gzzL,kmul(gxzL,gxzL)),kmul(gxxL,kmul(gyzL,gyzL)))));
    
    CCTK_REAL_VEC invdetg CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1),detg);
    
    CCTK_REAL_VEC gInv11 CCTK_ATTRIBUTE_UNUSED = 
      kmul(invdetg,kmsub(gyyL,gzzL,kmul(gyzL,gyzL)));
    
    CCTK_REAL_VEC gInv12 CCTK_ATTRIBUTE_UNUSED = 
      kmul(invdetg,kmsub(gxzL,gyzL,kmul(gxyL,gzzL)));
    
    CCTK_REAL_VEC gInv13 CCTK_ATTRIBUTE_UNUSED = 
      kmul(invdetg,kmsub(gxyL,gyzL,kmul(gxzL,gyyL)));
    
    CCTK_REAL_VEC gInv21 CCTK_ATTRIBUTE_UNUSED = 
      kmul(invdetg,kmsub(gxzL,gyzL,kmul(gxyL,gzzL)));
    
    CCTK_REAL_VEC gInv22 CCTK_ATTRIBUTE_UNUSED = 
      kmul(invdetg,kmsub(gxxL,gzzL,kmul(gxzL,gxzL)));
    
    CCTK_REAL_VEC gInv23 CCTK_ATTRIBUTE_UNUSED = 
      kmul(invdetg,kmsub(gxyL,gxzL,kmul(gxxL,gyzL)));
    
    CCTK_REAL_VEC gInv31 CCTK_ATTRIBUTE_UNUSED = 
      kmul(invdetg,kmsub(gxyL,gyzL,kmul(gxzL,gyyL)));
    
    CCTK_REAL_VEC gInv32 CCTK_ATTRIBUTE_UNUSED = 
      kmul(invdetg,kmsub(gxyL,gxzL,kmul(gxxL,gyzL)));
    
    CCTK_REAL_VEC gInv33 CCTK_ATTRIBUTE_UNUSED = 
      kmul(invdetg,kmsub(gxxL,gyyL,kmul(gxyL,gxyL)));
    
    CCTK_REAL_VEC gamma111 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gInv11,JacPDstandard4th1gxx,kmsub(ToReal(2),kmadd(gInv12,JacPDstandard4th1gxy,kmul(gInv13,JacPDstandard4th1gxz)),kmadd(gInv13,JacPDstandard4th3gxx,kmul(gInv12,JacPDstandard4th2gxx)))),ToReal(0.5));
    
    CCTK_REAL_VEC gamma211 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gInv21,JacPDstandard4th1gxx,kmsub(ToReal(2),kmadd(gInv22,JacPDstandard4th1gxy,kmul(gInv23,JacPDstandard4th1gxz)),kmadd(gInv23,JacPDstandard4th3gxx,kmul(gInv22,JacPDstandard4th2gxx)))),ToReal(0.5));
    
    CCTK_REAL_VEC gamma311 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gInv31,JacPDstandard4th1gxx,kmsub(ToReal(2),kmadd(gInv32,JacPDstandard4th1gxy,kmul(gInv33,JacPDstandard4th1gxz)),kmadd(gInv33,JacPDstandard4th3gxx,kmul(gInv32,JacPDstandard4th2gxx)))),ToReal(0.5));
    
    CCTK_REAL_VEC gamma121 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gInv12,JacPDstandard4th1gyy,kmadd(gInv11,JacPDstandard4th2gxx,kmul(gInv13,kadd(JacPDstandard4th1gyz,ksub(JacPDstandard4th2gxz,JacPDstandard4th3gxy))))),ToReal(0.5));
    
    CCTK_REAL_VEC gamma221 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gInv22,JacPDstandard4th1gyy,kmadd(gInv21,JacPDstandard4th2gxx,kmul(gInv23,kadd(JacPDstandard4th1gyz,ksub(JacPDstandard4th2gxz,JacPDstandard4th3gxy))))),ToReal(0.5));
    
    CCTK_REAL_VEC gamma321 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gInv32,JacPDstandard4th1gyy,kmadd(gInv31,JacPDstandard4th2gxx,kmul(gInv33,kadd(JacPDstandard4th1gyz,ksub(JacPDstandard4th2gxz,JacPDstandard4th3gxy))))),ToReal(0.5));
    
    CCTK_REAL_VEC gamma131 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gInv13,JacPDstandard4th1gzz,kmadd(gInv11,JacPDstandard4th3gxx,kmul(gInv12,kadd(JacPDstandard4th1gyz,ksub(JacPDstandard4th3gxy,JacPDstandard4th2gxz))))),ToReal(0.5));
    
    CCTK_REAL_VEC gamma231 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gInv23,JacPDstandard4th1gzz,kmadd(gInv21,JacPDstandard4th3gxx,kmul(gInv22,kadd(JacPDstandard4th1gyz,ksub(JacPDstandard4th3gxy,JacPDstandard4th2gxz))))),ToReal(0.5));
    
    CCTK_REAL_VEC gamma331 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gInv33,JacPDstandard4th1gzz,kmadd(gInv31,JacPDstandard4th3gxx,kmul(gInv32,kadd(JacPDstandard4th1gyz,ksub(JacPDstandard4th3gxy,JacPDstandard4th2gxz))))),ToReal(0.5));
    
    CCTK_REAL_VEC gamma122 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gInv11,kmsub(ToReal(2),JacPDstandard4th2gxy,JacPDstandard4th1gyy),kmadd(gInv12,JacPDstandard4th2gyy,kmul(gInv13,kmsub(ToReal(2),JacPDstandard4th2gyz,JacPDstandard4th3gyy)))),ToReal(0.5));
    
    CCTK_REAL_VEC gamma222 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gInv21,kmsub(ToReal(2),JacPDstandard4th2gxy,JacPDstandard4th1gyy),kmadd(gInv22,JacPDstandard4th2gyy,kmul(gInv23,kmsub(ToReal(2),JacPDstandard4th2gyz,JacPDstandard4th3gyy)))),ToReal(0.5));
    
    CCTK_REAL_VEC gamma322 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gInv31,kmsub(ToReal(2),JacPDstandard4th2gxy,JacPDstandard4th1gyy),kmadd(gInv32,JacPDstandard4th2gyy,kmul(gInv33,kmsub(ToReal(2),JacPDstandard4th2gyz,JacPDstandard4th3gyy)))),ToReal(0.5));
    
    CCTK_REAL_VEC gamma132 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gInv13,JacPDstandard4th2gzz,kmadd(gInv11,ksub(kadd(JacPDstandard4th2gxz,JacPDstandard4th3gxy),JacPDstandard4th1gyz),kmul(gInv12,JacPDstandard4th3gyy))),ToReal(0.5));
    
    CCTK_REAL_VEC gamma232 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gInv23,JacPDstandard4th2gzz,kmadd(gInv21,ksub(kadd(JacPDstandard4th2gxz,JacPDstandard4th3gxy),JacPDstandard4th1gyz),kmul(gInv22,JacPDstandard4th3gyy))),ToReal(0.5));
    
    CCTK_REAL_VEC gamma332 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gInv33,JacPDstandard4th2gzz,kmadd(gInv31,ksub(kadd(JacPDstandard4th2gxz,JacPDstandard4th3gxy),JacPDstandard4th1gyz),kmul(gInv32,JacPDstandard4th3gyy))),ToReal(0.5));
    
    CCTK_REAL_VEC gamma133 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gInv11,kmsub(ToReal(2),JacPDstandard4th3gxz,JacPDstandard4th1gzz),kmadd(gInv12,kmsub(ToReal(2),JacPDstandard4th3gyz,JacPDstandard4th2gzz),kmul(gInv13,JacPDstandard4th3gzz))),ToReal(0.5));
    
    CCTK_REAL_VEC gamma233 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gInv21,kmsub(ToReal(2),JacPDstandard4th3gxz,JacPDstandard4th1gzz),kmadd(gInv22,kmsub(ToReal(2),JacPDstandard4th3gyz,JacPDstandard4th2gzz),kmul(gInv23,JacPDstandard4th3gzz))),ToReal(0.5));
    
    CCTK_REAL_VEC gamma333 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gInv31,kmsub(ToReal(2),JacPDstandard4th3gxz,JacPDstandard4th1gzz),kmadd(gInv32,kmsub(ToReal(2),JacPDstandard4th3gyz,JacPDstandard4th2gzz),kmul(gInv33,JacPDstandard4th3gzz))),ToReal(0.5));
    
    CCTK_REAL_VEC xmoved CCTK_ATTRIBUTE_UNUSED = ksub(xL,ToReal(xorig));
    
    CCTK_REAL_VEC ymoved CCTK_ATTRIBUTE_UNUSED = ksub(yL,ToReal(yorig));
    
    CCTK_REAL_VEC zmoved CCTK_ATTRIBUTE_UNUSED = ksub(zL,ToReal(zorig));
    
    CCTK_REAL_VEC va1 CCTK_ATTRIBUTE_UNUSED = kneg(ymoved);
    
    CCTK_REAL_VEC va2 CCTK_ATTRIBUTE_UNUSED = kadd(xmoved,ToReal(offset));
    
    CCTK_REAL_VEC va3 CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC vb1 CCTK_ATTRIBUTE_UNUSED = kadd(xmoved,ToReal(offset));
    
    CCTK_REAL_VEC vb2 CCTK_ATTRIBUTE_UNUSED = ymoved;
    
    CCTK_REAL_VEC vb3 CCTK_ATTRIBUTE_UNUSED = zmoved;
    
    CCTK_REAL_VEC vc1 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(kmsub(gInv12,va3,kmul(gInv13,va2)),vb1,kmadd(kmsub(gInv13,va1,kmul(gInv11,va3)),vb2,kmul(vb3,kmsub(gInv11,va2,kmul(gInv12,va1))))),kpow(detg,0.5));
    
    CCTK_REAL_VEC vc2 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(kmsub(gInv22,va3,kmul(gInv23,va2)),vb1,kmadd(kmsub(gInv23,va1,kmul(gInv21,va3)),vb2,kmul(vb3,kmsub(gInv21,va2,kmul(gInv22,va1))))),kpow(detg,0.5));
    
    CCTK_REAL_VEC vc3 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(kmsub(gInv32,va3,kmul(gInv33,va2)),vb1,kmadd(kmsub(gInv33,va1,kmul(gInv31,va3)),vb2,kmul(vb3,kmsub(gInv31,va2,kmul(gInv32,va1))))),kpow(detg,0.5));
    
    CCTK_REAL_VEC wa1 CCTK_ATTRIBUTE_UNUSED = va1;
    
    CCTK_REAL_VEC wa2 CCTK_ATTRIBUTE_UNUSED = va2;
    
    CCTK_REAL_VEC wa3 CCTK_ATTRIBUTE_UNUSED = va3;
    
    CCTK_REAL_VEC omega11 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(ToReal(2),kmadd(gyzL,kmul(wa2,wa3),kmul(wa1,kmadd(gxyL,wa2,kmul(gxzL,wa3)))),kmadd(gxxL,kmul(wa1,wa1),kmadd(gyyL,kmul(wa2,wa2),kmul(gzzL,kmul(wa3,wa3)))));
    
    CCTK_REAL_VEC ea1 CCTK_ATTRIBUTE_UNUSED = 
      kmul(wa1,kpow(omega11,-0.5));
    
    CCTK_REAL_VEC ea2 CCTK_ATTRIBUTE_UNUSED = 
      kmul(wa2,kpow(omega11,-0.5));
    
    CCTK_REAL_VEC ea3 CCTK_ATTRIBUTE_UNUSED = 
      kmul(wa3,kpow(omega11,-0.5));
    
    CCTK_REAL_VEC omega12 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(ea1,kmadd(gxxL,vb1,kmadd(gxyL,vb2,kmul(gxzL,vb3))),kmadd(ea2,kmadd(gxyL,vb1,kmadd(gyyL,vb2,kmul(gyzL,vb3))),kmul(ea3,kmadd(gxzL,vb1,kmadd(gyzL,vb2,kmul(gzzL,vb3))))));
    
    CCTK_REAL_VEC wb1 CCTK_ATTRIBUTE_UNUSED = knmsub(ea1,omega12,vb1);
    
    CCTK_REAL_VEC wb2 CCTK_ATTRIBUTE_UNUSED = knmsub(ea2,omega12,vb2);
    
    CCTK_REAL_VEC wb3 CCTK_ATTRIBUTE_UNUSED = knmsub(ea3,omega12,vb3);
    
    CCTK_REAL_VEC omega22 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(ToReal(2),kmadd(gyzL,kmul(wb2,wb3),kmul(wb1,kmadd(gxyL,wb2,kmul(gxzL,wb3)))),kmadd(gxxL,kmul(wb1,wb1),kmadd(gyyL,kmul(wb2,wb2),kmul(gzzL,kmul(wb3,wb3)))));
    
    CCTK_REAL_VEC eb1 CCTK_ATTRIBUTE_UNUSED = 
      kmul(wb1,kpow(omega22,-0.5));
    
    CCTK_REAL_VEC eb2 CCTK_ATTRIBUTE_UNUSED = 
      kmul(wb2,kpow(omega22,-0.5));
    
    CCTK_REAL_VEC eb3 CCTK_ATTRIBUTE_UNUSED = 
      kmul(wb3,kpow(omega22,-0.5));
    
    CCTK_REAL_VEC omega13 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(ea1,kmadd(gxxL,vc1,kmadd(gxyL,vc2,kmul(gxzL,vc3))),kmadd(ea2,kmadd(gxyL,vc1,kmadd(gyyL,vc2,kmul(gyzL,vc3))),kmul(ea3,kmadd(gxzL,vc1,kmadd(gyzL,vc2,kmul(gzzL,vc3))))));
    
    CCTK_REAL_VEC omega23 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(eb1,kmadd(gxxL,vc1,kmadd(gxyL,vc2,kmul(gxzL,vc3))),kmadd(eb2,kmadd(gxyL,vc1,kmadd(gyyL,vc2,kmul(gyzL,vc3))),kmul(eb3,kmadd(gxzL,vc1,kmadd(gyzL,vc2,kmul(gzzL,vc3))))));
    
    CCTK_REAL_VEC wc1 CCTK_ATTRIBUTE_UNUSED = 
      knmsub(ea1,omega13,knmsub(eb1,omega23,vc1));
    
    CCTK_REAL_VEC wc2 CCTK_ATTRIBUTE_UNUSED = 
      knmsub(ea2,omega13,knmsub(eb2,omega23,vc2));
    
    CCTK_REAL_VEC wc3 CCTK_ATTRIBUTE_UNUSED = 
      knmsub(ea3,omega13,knmsub(eb3,omega23,vc3));
    
    CCTK_REAL_VEC omega33 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(ToReal(2),kmadd(gyzL,kmul(wc2,wc3),kmul(wc1,kmadd(gxyL,wc2,kmul(gxzL,wc3)))),kmadd(gxxL,kmul(wc1,wc1),kmadd(gyyL,kmul(wc2,wc2),kmul(gzzL,kmul(wc3,wc3)))));
    
    CCTK_REAL_VEC ec1 CCTK_ATTRIBUTE_UNUSED = 
      kmul(wc1,kpow(omega33,-0.5));
    
    CCTK_REAL_VEC ec2 CCTK_ATTRIBUTE_UNUSED = 
      kmul(wc2,kpow(omega33,-0.5));
    
    CCTK_REAL_VEC ec3 CCTK_ATTRIBUTE_UNUSED = 
      kmul(wc3,kpow(omega33,-0.5));
    
    CCTK_REAL_VEC isqrt2 CCTK_ATTRIBUTE_UNUSED = 
      ToReal(0.707106781186547524);
    
    CCTK_REAL_VEC ltet1 CCTK_ATTRIBUTE_UNUSED = kmul(eb1,isqrt2);
    
    CCTK_REAL_VEC ltet2 CCTK_ATTRIBUTE_UNUSED = kmul(eb2,isqrt2);
    
    CCTK_REAL_VEC ltet3 CCTK_ATTRIBUTE_UNUSED = kmul(eb3,isqrt2);
    
    CCTK_REAL_VEC n1 CCTK_ATTRIBUTE_UNUSED = kneg(kmul(eb1,isqrt2));
    
    CCTK_REAL_VEC n2 CCTK_ATTRIBUTE_UNUSED = kneg(kmul(eb2,isqrt2));
    
    CCTK_REAL_VEC n3 CCTK_ATTRIBUTE_UNUSED = kneg(kmul(eb3,isqrt2));
    
    CCTK_REAL_VEC rm1 CCTK_ATTRIBUTE_UNUSED = kmul(ec1,isqrt2);
    
    CCTK_REAL_VEC rm2 CCTK_ATTRIBUTE_UNUSED = kmul(ec2,isqrt2);
    
    CCTK_REAL_VEC rm3 CCTK_ATTRIBUTE_UNUSED = kmul(ec3,isqrt2);
    
    CCTK_REAL_VEC im1 CCTK_ATTRIBUTE_UNUSED = kmul(ea1,isqrt2);
    
    CCTK_REAL_VEC im2 CCTK_ATTRIBUTE_UNUSED = kmul(ea2,isqrt2);
    
    CCTK_REAL_VEC im3 CCTK_ATTRIBUTE_UNUSED = kmul(ea3,isqrt2);
    
    CCTK_REAL_VEC rmbar1 CCTK_ATTRIBUTE_UNUSED = kmul(ec1,isqrt2);
    
    CCTK_REAL_VEC rmbar2 CCTK_ATTRIBUTE_UNUSED = kmul(ec2,isqrt2);
    
    CCTK_REAL_VEC rmbar3 CCTK_ATTRIBUTE_UNUSED = kmul(ec3,isqrt2);
    
    CCTK_REAL_VEC imbar1 CCTK_ATTRIBUTE_UNUSED = kneg(kmul(ea1,isqrt2));
    
    CCTK_REAL_VEC imbar2 CCTK_ATTRIBUTE_UNUSED = kneg(kmul(ea2,isqrt2));
    
    CCTK_REAL_VEC imbar3 CCTK_ATTRIBUTE_UNUSED = kneg(kmul(ea3,isqrt2));
    
    CCTK_REAL_VEC nn CCTK_ATTRIBUTE_UNUSED = isqrt2;
    
    CCTK_REAL_VEC R1212 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(ToReal(2),kmadd(gamma121,kmadd(gxxL,gamma121,kmadd(gxyL,gamma221,kmul(gxzL,gamma321))),kmadd(gamma221,kmadd(gxyL,gamma121,kmadd(gyyL,gamma221,kmul(gyzL,gamma321))),kmul(gamma321,kmadd(gxzL,gamma121,kmadd(gyzL,gamma221,kmul(gzzL,gamma321)))))),kmadd(ToReal(-2),kmadd(gamma122,kmadd(gxxL,gamma111,kmadd(gxyL,gamma211,kmul(gxzL,gamma311))),kmadd(gamma222,kmadd(gxyL,gamma111,kmadd(gyyL,gamma211,kmul(gyzL,gamma311))),kmul(kmadd(gxzL,gamma111,kmadd(gyzL,gamma211,kmul(gzzL,gamma311))),gamma322))),ksub(kadd(JacPDstandard4th12gxy,ksub(JacPDstandard4th21gxy,JacPDstandard4th22gxx)),JacPDstandard4th11gyy))),ToReal(0.5));
    
    CCTK_REAL_VEC R1213 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(ToReal(2),kmadd(gamma121,kmadd(gxxL,gamma131,kmadd(gxyL,gamma231,kmul(gxzL,gamma331))),kmadd(gamma221,kmadd(gxyL,gamma131,kmadd(gyyL,gamma231,kmul(gyzL,gamma331))),kmul(gamma321,kmadd(gxzL,gamma131,kmadd(gyzL,gamma231,kmul(gzzL,gamma331)))))),kmadd(ToReal(-2),kmadd(gamma132,kmadd(gxxL,gamma111,kmadd(gxyL,gamma211,kmul(gxzL,gamma311))),kmadd(gamma232,kmadd(gxyL,gamma111,kmadd(gyyL,gamma211,kmul(gyzL,gamma311))),kmul(kmadd(gxzL,gamma111,kmadd(gyzL,gamma211,kmul(gzzL,gamma311))),gamma332))),ksub(kadd(JacPDstandard4th12gxz,ksub(JacPDstandard4th31gxy,JacPDstandard4th23gxx)),JacPDstandard4th11gyz))),ToReal(0.5));
    
    CCTK_REAL_VEC R1223 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(ToReal(2),kmadd(gamma122,kmadd(gxxL,gamma131,kmadd(gxyL,gamma231,kmul(gxzL,gamma331))),kmadd(gamma222,kmadd(gxyL,gamma131,kmadd(gyyL,gamma231,kmul(gyzL,gamma331))),kmul(gamma322,kmadd(gxzL,gamma131,kmadd(gyzL,gamma231,kmul(gzzL,gamma331)))))),kmadd(ToReal(-2),kmadd(gamma132,kmadd(gxxL,gamma121,kmadd(gxyL,gamma221,kmul(gxzL,gamma321))),kmadd(gamma232,kmadd(gxyL,gamma121,kmadd(gyyL,gamma221,kmul(gyzL,gamma321))),kmul(kmadd(gxzL,gamma121,kmadd(gyzL,gamma221,kmul(gzzL,gamma321))),gamma332))),ksub(kadd(JacPDstandard4th22gxz,ksub(JacPDstandard4th31gyy,JacPDstandard4th23gxy)),JacPDstandard4th12gyz))),ToReal(0.5));
    
    CCTK_REAL_VEC R1313 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(ToReal(2),kmadd(gamma131,kmadd(gxxL,gamma131,kmadd(gxyL,gamma231,kmul(gxzL,gamma331))),kmadd(gamma231,kmadd(gxyL,gamma131,kmadd(gyyL,gamma231,kmul(gyzL,gamma331))),kmul(gamma331,kmadd(gxzL,gamma131,kmadd(gyzL,gamma231,kmul(gzzL,gamma331)))))),kmadd(ToReal(-2),kmadd(gamma133,kmadd(gxxL,gamma111,kmadd(gxyL,gamma211,kmul(gxzL,gamma311))),kmadd(gamma233,kmadd(gxyL,gamma111,kmadd(gyyL,gamma211,kmul(gyzL,gamma311))),kmul(kmadd(gxzL,gamma111,kmadd(gyzL,gamma211,kmul(gzzL,gamma311))),gamma333))),ksub(kadd(JacPDstandard4th13gxz,ksub(JacPDstandard4th31gxz,JacPDstandard4th33gxx)),JacPDstandard4th11gzz))),ToReal(0.5));
    
    CCTK_REAL_VEC R1323 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(ToReal(2),kmadd(gamma132,kmadd(gxxL,gamma131,kmadd(gxyL,gamma231,kmul(gxzL,gamma331))),kmadd(gamma232,kmadd(gxyL,gamma131,kmadd(gyyL,gamma231,kmul(gyzL,gamma331))),kmul(kmadd(gxzL,gamma131,kmadd(gyzL,gamma231,kmul(gzzL,gamma331))),gamma332))),kmadd(ToReal(-2),kmadd(gamma133,kmadd(gxxL,gamma121,kmadd(gxyL,gamma221,kmul(gxzL,gamma321))),kmadd(gamma233,kmadd(gxyL,gamma121,kmadd(gyyL,gamma221,kmul(gyzL,gamma321))),kmul(kmadd(gxzL,gamma121,kmadd(gyzL,gamma221,kmul(gzzL,gamma321))),gamma333))),ksub(kadd(JacPDstandard4th23gxz,ksub(JacPDstandard4th31gyz,JacPDstandard4th33gxy)),JacPDstandard4th12gzz))),ToReal(0.5));
    
    CCTK_REAL_VEC R2323 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(ToReal(2),kmadd(gamma132,kmadd(gxxL,gamma132,kmadd(gxyL,gamma232,kmul(gxzL,gamma332))),kmadd(gamma232,kmadd(gxyL,gamma132,kmadd(gyyL,gamma232,kmul(gyzL,gamma332))),kmul(gamma332,kmadd(gxzL,gamma132,kmadd(gyzL,gamma232,kmul(gzzL,gamma332)))))),kmadd(ToReal(-2),kmadd(gamma133,kmadd(gxxL,gamma122,kmadd(gxyL,gamma222,kmul(gxzL,gamma322))),kmadd(gamma233,kmadd(gxyL,gamma122,kmadd(gyyL,gamma222,kmul(gyzL,gamma322))),kmul(kmadd(gxzL,gamma122,kmadd(gyzL,gamma222,kmul(gzzL,gamma322))),gamma333))),ksub(kadd(JacPDstandard4th23gyz,ksub(JacPDstandard4th32gyz,JacPDstandard4th33gyy)),JacPDstandard4th22gzz))),ToReal(0.5));
    
    CCTK_REAL_VEC R4p1212 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(kxxL,kyyL,knmsub(kxyL,kxyL,R1212));
    
    CCTK_REAL_VEC R4p1213 CCTK_ATTRIBUTE_UNUSED = 
      knmsub(kxyL,kxzL,kmadd(kxxL,kyzL,R1213));
    
    CCTK_REAL_VEC R4p1223 CCTK_ATTRIBUTE_UNUSED = 
      knmsub(kxzL,kyyL,kmadd(kxyL,kyzL,R1223));
    
    CCTK_REAL_VEC R4p1313 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(kxxL,kzzL,knmsub(kxzL,kxzL,R1313));
    
    CCTK_REAL_VEC R4p1323 CCTK_ATTRIBUTE_UNUSED = 
      knmsub(kxzL,kyzL,kmadd(kxyL,kzzL,R1323));
    
    CCTK_REAL_VEC R4p2323 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(kyyL,kzzL,knmsub(kyzL,kyzL,R2323));
    
    CCTK_REAL_VEC Ro111 CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC Ro112 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(kxxL,gamma121,knmsub(kyyL,gamma211,kmadd(kxyL,ksub(gamma221,gamma111),knmsub(kyzL,gamma311,kmadd(kxzL,gamma321,ksub(JacPDstandard4th1kxy,JacPDstandard4th2kxx))))));
    
    CCTK_REAL_VEC Ro113 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(kxxL,gamma131,knmsub(kyzL,gamma211,kmadd(kxyL,gamma231,knmsub(kzzL,gamma311,kmadd(kxzL,ksub(gamma331,gamma111),ksub(JacPDstandard4th1kxz,JacPDstandard4th3kxx))))));
    
    CCTK_REAL_VEC Ro121 CCTK_ATTRIBUTE_UNUSED = 
      knmsub(kxxL,gamma121,kmadd(kyyL,gamma211,kmadd(kxyL,ksub(gamma111,gamma221),kmadd(kyzL,gamma311,knmsub(kxzL,gamma321,ksub(JacPDstandard4th2kxx,JacPDstandard4th1kxy))))));
    
    CCTK_REAL_VEC Ro122 CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC Ro123 CCTK_ATTRIBUTE_UNUSED = 
      knmsub(kxzL,gamma121,kmadd(kxyL,gamma131,kmadd(kyyL,gamma231,knmsub(kzzL,gamma321,kmadd(kyzL,ksub(gamma331,gamma221),ksub(JacPDstandard4th2kxz,JacPDstandard4th3kxy))))));
    
    CCTK_REAL_VEC Ro131 CCTK_ATTRIBUTE_UNUSED = 
      knmsub(kxxL,gamma131,kmadd(kyzL,gamma211,knmsub(kxyL,gamma231,kmadd(kzzL,gamma311,kmadd(kxzL,ksub(gamma111,gamma331),ksub(JacPDstandard4th3kxx,JacPDstandard4th1kxz))))));
    
    CCTK_REAL_VEC Ro132 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(kxzL,gamma121,knmsub(kxyL,gamma131,knmsub(kyyL,gamma231,kmadd(kzzL,gamma321,kmadd(kyzL,ksub(gamma221,gamma331),ksub(JacPDstandard4th3kxy,JacPDstandard4th2kxz))))));
    
    CCTK_REAL_VEC Ro133 CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC Ro211 CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC Ro212 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(kxxL,gamma122,knmsub(kyyL,gamma221,kmadd(kxyL,ksub(gamma222,gamma121),knmsub(kyzL,gamma321,kmadd(kxzL,gamma322,ksub(JacPDstandard4th1kyy,JacPDstandard4th2kxy))))));
    
    CCTK_REAL_VEC Ro213 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(kxxL,gamma132,knmsub(kyzL,gamma221,kmadd(kxyL,gamma232,knmsub(kzzL,gamma321,kmadd(kxzL,ksub(gamma332,gamma121),ksub(JacPDstandard4th1kyz,JacPDstandard4th3kxy))))));
    
    CCTK_REAL_VEC Ro221 CCTK_ATTRIBUTE_UNUSED = 
      knmsub(kxxL,gamma122,kmadd(kyyL,gamma221,kmadd(kxyL,ksub(gamma121,gamma222),kmadd(kyzL,gamma321,knmsub(kxzL,gamma322,ksub(JacPDstandard4th2kxy,JacPDstandard4th1kyy))))));
    
    CCTK_REAL_VEC Ro222 CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC Ro223 CCTK_ATTRIBUTE_UNUSED = 
      knmsub(kxzL,gamma122,kmadd(kxyL,gamma132,kmadd(kyyL,gamma232,knmsub(kzzL,gamma322,kmadd(kyzL,ksub(gamma332,gamma222),ksub(JacPDstandard4th2kyz,JacPDstandard4th3kyy))))));
    
    CCTK_REAL_VEC Ro231 CCTK_ATTRIBUTE_UNUSED = 
      knmsub(kxxL,gamma132,kmadd(kyzL,gamma221,knmsub(kxyL,gamma232,kmadd(kzzL,gamma321,kmadd(kxzL,ksub(gamma121,gamma332),ksub(JacPDstandard4th3kxy,JacPDstandard4th1kyz))))));
    
    CCTK_REAL_VEC Ro232 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(kxzL,gamma122,knmsub(kxyL,gamma132,knmsub(kyyL,gamma232,kmadd(kzzL,gamma322,kmadd(kyzL,ksub(gamma222,gamma332),ksub(JacPDstandard4th3kyy,JacPDstandard4th2kyz))))));
    
    CCTK_REAL_VEC Ro233 CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC Ro311 CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC Ro312 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(kxxL,gamma132,knmsub(kyyL,gamma231,kmadd(kxyL,ksub(gamma232,gamma131),knmsub(kyzL,gamma331,kmadd(kxzL,gamma332,ksub(JacPDstandard4th1kyz,JacPDstandard4th2kxz))))));
    
    CCTK_REAL_VEC Ro313 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(kxxL,gamma133,knmsub(kyzL,gamma231,kmadd(kxyL,gamma233,knmsub(kzzL,gamma331,kmadd(kxzL,ksub(gamma333,gamma131),ksub(JacPDstandard4th1kzz,JacPDstandard4th3kxz))))));
    
    CCTK_REAL_VEC Ro321 CCTK_ATTRIBUTE_UNUSED = 
      knmsub(kxxL,gamma132,kmadd(kyyL,gamma231,kmadd(kxyL,ksub(gamma131,gamma232),kmadd(kyzL,gamma331,knmsub(kxzL,gamma332,ksub(JacPDstandard4th2kxz,JacPDstandard4th1kyz))))));
    
    CCTK_REAL_VEC Ro322 CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC Ro323 CCTK_ATTRIBUTE_UNUSED = 
      knmsub(kxzL,gamma132,kmadd(kxyL,gamma133,kmadd(kyyL,gamma233,knmsub(kzzL,gamma332,kmadd(kyzL,ksub(gamma333,gamma232),ksub(JacPDstandard4th2kzz,JacPDstandard4th3kyz))))));
    
    CCTK_REAL_VEC Ro331 CCTK_ATTRIBUTE_UNUSED = 
      knmsub(kxxL,gamma133,kmadd(kyzL,gamma231,knmsub(kxyL,gamma233,kmadd(kzzL,gamma331,kmadd(kxzL,ksub(gamma131,gamma333),ksub(JacPDstandard4th3kxz,JacPDstandard4th1kzz))))));
    
    CCTK_REAL_VEC Ro332 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(kxzL,gamma132,knmsub(kxyL,gamma133,knmsub(kyyL,gamma233,kmadd(kzzL,gamma332,kmadd(kyzL,ksub(gamma232,gamma333),ksub(JacPDstandard4th3kyz,JacPDstandard4th2kzz))))));
    
    CCTK_REAL_VEC Ro333 CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC Rojo11 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(kadd(gInv23,gInv32),knmsub(kxyL,kxzL,kmadd(kxxL,kyzL,R1213)),kmadd(gInv22,kmadd(kxxL,kyyL,knmsub(kxyL,kxyL,R1212)),kmul(gInv33,kmadd(kxxL,kzzL,knmsub(kxzL,kxzL,R1313)))));
    
    CCTK_REAL_VEC Rojo12 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(kmsub(kxyL,kxzL,kmul(kxxL,kyzL)),gInv13,kmadd(kmsub(kxyL,kyzL,kmul(kxzL,kyyL)),gInv32,knmsub(gInv21,R1212,knmsub(gInv31,R1213,kmadd(gInv23,R1223,kmadd(gInv33,knmsub(kxzL,kyzL,kmadd(kxyL,kzzL,R1323)),kmul(gInv12,kmsub(kxyL,kxyL,kmul(kxxL,kyyL)))))))));
    
    CCTK_REAL_VEC Rojo13 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(kmsub(kxyL,kxzL,kmul(kxxL,kyzL)),gInv12,kmadd(kmsub(kxzL,kyzL,kmul(kxyL,kzzL)),gInv23,knmsub(gInv21,R1213,kmadd(gInv22,kmsub(kxzL,kyyL,kmadd(kxyL,kyzL,R1223)),knmsub(gInv31,R1313,kmsub(gInv13,kmsub(kxzL,kxzL,kmul(kxxL,kzzL)),kmul(gInv32,R1323)))))));
    
    CCTK_REAL_VEC Rojo21 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(kmsub(kxyL,kyzL,kmul(kxzL,kyyL)),gInv23,kmadd(kmsub(kxyL,kxzL,kmul(kxxL,kyzL)),gInv31,knmsub(gInv12,R1212,knmsub(gInv13,R1213,kmadd(gInv32,R1223,kmadd(gInv33,knmsub(kxzL,kyzL,kmadd(kxyL,kzzL,R1323)),kmul(gInv21,kmsub(kxyL,kxyL,kmul(kxxL,kyyL)))))))));
    
    CCTK_REAL_VEC Rojo22 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(kadd(gInv13,gInv31),kmsub(kxzL,kyyL,kmadd(kxyL,kyzL,R1223)),kmadd(gInv11,kmadd(kxxL,kyyL,knmsub(kxyL,kxyL,R1212)),kmul(gInv33,kmadd(kyyL,kzzL,knmsub(kyzL,kyzL,R2323)))));
    
    CCTK_REAL_VEC Rojo23 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(kmsub(kxzL,kyzL,kmul(kxyL,kzzL)),gInv13,kmadd(kmsub(kxyL,kyzL,kmul(kxzL,kyyL)),gInv21,kmadd(gInv11,knmsub(kxyL,kxzL,kmadd(kxxL,kyzL,R1213)),kmadd(gInv12,R1223,knmsub(gInv31,R1323,kmsub(gInv23,kmsub(kyzL,kyzL,kmul(kyyL,kzzL)),kmul(gInv32,R2323)))))));
    
    CCTK_REAL_VEC Rojo31 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(kmsub(kxyL,kxzL,kmul(kxxL,kyzL)),gInv21,kmadd(kmsub(kxzL,kyzL,kmul(kxyL,kzzL)),gInv32,knmsub(gInv12,R1213,kmadd(gInv22,kmsub(kxzL,kyyL,kmadd(kxyL,kyzL,R1223)),knmsub(gInv13,R1313,kmsub(gInv31,kmsub(kxzL,kxzL,kmul(kxxL,kzzL)),kmul(gInv23,R1323)))))));
    
    CCTK_REAL_VEC Rojo32 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(kmsub(kxyL,kyzL,kmul(kxzL,kyyL)),gInv12,kmadd(kmsub(kxzL,kyzL,kmul(kxyL,kzzL)),gInv31,kmadd(gInv11,knmsub(kxyL,kxzL,kmadd(kxxL,kyzL,R1213)),kmadd(gInv21,R1223,knmsub(gInv13,R1323,kmsub(gInv32,kmsub(kyzL,kyzL,kmul(kyyL,kzzL)),kmul(gInv23,R2323)))))));
    
    CCTK_REAL_VEC Rojo33 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(kadd(gInv12,gInv21),knmsub(kxzL,kyzL,kmadd(kxyL,kzzL,R1323)),kmadd(gInv11,kmadd(kxxL,kzzL,knmsub(kxzL,kxzL,R1313)),kmul(gInv22,kmadd(kyyL,kzzL,knmsub(kyzL,kyzL,R2323)))));
    
    CCTK_REAL_VEC Psi4rL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(ToReal(2),kmul(kmsub(imbar1,imbar2,kmul(rmbar1,rmbar2)),kmsub(n1,kmadd(n2,R4p1212,kmul(n3,R4p1213)),kmul(n3,kmadd(n2,R4p1223,kmul(n3,R4p1323))))),kmadd(ToReal(2),kmul(kmadd(n1,kmsub(n2,R4p1223,kmul(n3,R4p1323)),kmsub(R4p1213,kmul(n1,n1),kmul(n2,kmul(n3,R4p2323)))),kmsub(rmbar2,rmbar3,kmul(imbar2,imbar3))),kmadd(ToReal(2),kmul(kmadd(n1,kmul(n2,R4p1213),kmadd(n1,kmul(n3,R4p1313),kmadd(n2,kmul(n3,R4p1323),kmul(R4p1223,kmul(n2,n2))))),kmsub(imbar1,imbar3,kmul(rmbar1,rmbar3))),knmsub(kmadd(ToReal(2),kmul(n2,kmul(n3,R4p1213)),kmadd(R4p1212,kmul(n2,n2),kmul(R4p1313,kmul(n3,n3)))),kmsub(imbar1,imbar1,kmul(rmbar1,rmbar1)),knmsub(kmadd(ToReal(-2),kmul(n1,kmul(n3,R4p1223)),kmadd(R4p1212,kmul(n1,n1),kmul(R4p2323,kmul(n3,n3)))),kmsub(imbar2,imbar2,kmul(rmbar2,rmbar2)),knmsub(kmul(nn,nn),kmadd(kmsub(imbar2,imbar3,kmul(rmbar2,rmbar3)),Rojo23,kmadd(imbar1,kmadd(imbar2,kadd(Rojo12,Rojo21),kmul(imbar3,kadd(Rojo13,Rojo31))),knmsub(rmbar1,kmadd(rmbar2,kadd(Rojo12,Rojo21),kmul(rmbar3,kadd(Rojo13,Rojo31))),kmadd(kmsub(imbar2,imbar3,kmul(rmbar2,rmbar3)),Rojo32,kmadd(Rojo11,kmsub(imbar1,imbar1,kmul(rmbar1,rmbar1)),kmadd(Rojo22,kmsub(imbar2,imbar2,kmul(rmbar2,rmbar2)),kmul(Rojo33,kmsub(imbar3,imbar3,kmul(rmbar3,rmbar3))))))))),kmsub(ToReal(2),kmul(nn,kmadd(kmsub(rmbar1,rmbar2,kmul(imbar1,imbar2)),kmadd(n1,Ro112,kmadd(n2,Ro122,kmul(n3,Ro132))),kmadd(kmsub(rmbar1,rmbar3,kmul(imbar1,imbar3)),kmadd(n1,Ro113,kmadd(n2,Ro123,kmul(n3,Ro133))),kmadd(kmsub(rmbar1,rmbar2,kmul(imbar1,imbar2)),kmadd(n1,Ro211,kmadd(n2,Ro221,kmul(n3,Ro231))),kmadd(kmsub(rmbar2,rmbar3,kmul(imbar2,imbar3)),kmadd(n1,Ro213,kmadd(n2,Ro223,kmul(n3,Ro233))),kmadd(kmsub(rmbar1,rmbar3,kmul(imbar1,imbar3)),kmadd(n1,Ro311,kmadd(n2,Ro321,kmul(n3,Ro331))),kmadd(kmsub(rmbar2,rmbar3,kmul(imbar2,imbar3)),kmadd(n1,Ro312,kmadd(n2,Ro322,kmul(n3,Ro332))),kmadd(kmadd(n1,Ro111,kmadd(n2,Ro121,kmul(n3,Ro131))),kmsub(rmbar1,rmbar1,kmul(imbar1,imbar1)),kmadd(kmadd(n1,Ro212,kmadd(n2,Ro222,kmul(n3,Ro232))),kmsub(rmbar2,rmbar2,kmul(imbar2,imbar2)),kmul(kmadd(n1,Ro313,kmadd(n2,Ro323,kmul(n3,Ro333))),kmsub(rmbar3,rmbar3,kmul(imbar3,imbar3)))))))))))),kmul(kmadd(ToReal(2),kmul(n1,kmul(n2,R4p1323)),kmadd(R4p1313,kmul(n1,n1),kmul(R4p2323,kmul(n2,n2)))),kmsub(imbar3,imbar3,kmul(rmbar3,rmbar3))))))))));
    
    CCTK_REAL_VEC Psi4iL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(ToReal(2),kmadd(kmsub(n1,kmadd(n2,R4p1212,kmul(n3,R4p1213)),kmul(n3,kmadd(n2,R4p1223,kmul(n3,R4p1323)))),kmadd(im2,rm1,kmul(im1,rm2)),kmadd(nn,knmsub(kmadd(im2,rm1,kmul(im1,rm2)),kmadd(n1,kadd(Ro112,Ro211),kmadd(n3,kadd(Ro132,Ro231),kmul(n2,kadd(Ro122,Ro221)))),kmadd(ToReal(-2),kmadd(im1,kmul(rm1,kmadd(n1,Ro111,kmadd(n2,Ro121,kmul(n3,Ro131)))),kmul(im2,kmul(rm2,kmadd(n1,Ro212,kmadd(n2,Ro222,kmul(n3,Ro232)))))),knmsub(kmadd(im3,rm1,kmul(im1,rm3)),kmadd(n1,kadd(Ro113,Ro311),kmadd(n3,kadd(Ro133,Ro331),kmul(n2,kadd(Ro123,Ro321)))),kmsub(ToReal(-2),kmul(im3,kmul(rm3,kmadd(n1,Ro313,kmadd(n2,Ro323,kmul(n3,Ro333))))),kmul(kmadd(im3,rm2,kmul(im2,rm3)),kmadd(n1,kadd(Ro213,Ro312),kmadd(n3,kadd(Ro233,Ro332),kmul(n2,kadd(Ro223,Ro322))))))))),kmul(kmadd(im3,rm1,kmul(im1,rm3)),kmadd(n1,kmadd(n2,R4p1213,kmul(n3,R4p1313)),kmadd(n2,kmul(n3,R4p1323),kmul(R4p1223,kmul(n2,n2))))))),kmsub(ToReal(-2),kmadd(kmadd(im3,rm2,kmul(im2,rm3)),kmadd(n1,kmsub(n2,R4p1223,kmul(n3,R4p1323)),kmsub(R4p1213,kmul(n1,n1),kmul(n2,kmul(n3,R4p2323)))),kmadd(im3,kmul(rm3,kmadd(ToReal(2),kmul(n1,kmul(n2,R4p1323)),kmadd(R4p1313,kmul(n1,n1),kmul(R4p2323,kmul(n2,n2))))),kmadd(im1,kmul(rm1,kmadd(ToReal(2),kmul(n2,kmul(n3,R4p1213)),kmadd(R4p1212,kmul(n2,n2),kmul(R4p1313,kmul(n3,n3))))),kmul(im2,kmul(rm2,kmadd(ToReal(-2),kmul(n1,kmul(n3,R4p1223)),kmadd(R4p1212,kmul(n1,n1),kmul(R4p2323,kmul(n3,n3))))))))),kmul(kmadd(im1,kmadd(ToReal(2),kmul(rm1,Rojo11),kmadd(rm2,kadd(Rojo12,Rojo21),kmul(rm3,kadd(Rojo13,Rojo31)))),kmadd(im2,kmadd(rm1,kadd(Rojo12,Rojo21),kmadd(ToReal(2),kmul(rm2,Rojo22),kmul(rm3,kadd(Rojo23,Rojo32)))),kmul(im3,kmadd(rm1,kadd(Rojo13,Rojo31),kmadd(rm2,kadd(Rojo23,Rojo32),kmul(kmul(rm3,Rojo33),ToReal(2))))))),kmul(nn,nn))));
    
    CCTK_REAL_VEC Psi3rL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(n1,knmsub(ltet1,kmadd(n2,kmul(R4p1212,rm2),kmadd(n3,kmul(R4p1213,rm2),kmadd(n2,kmul(R4p1213,rm3),kmul(n3,kmul(R4p1313,rm3))))),kmadd(ltet3,kmadd(n1,kmul(R4p1213,rm2),kmadd(n2,kmul(R4p1223,rm2),kmadd(n1,kmul(R4p1313,rm3),kmul(n2,kmul(R4p1323,rm3))))),kmadd(ltet2,kmsub(n1,kmadd(R4p1212,rm2,kmul(R4p1213,rm3)),kmul(n3,kmadd(R4p1223,rm2,kmul(R4p1323,rm3)))),kmadd(kmadd(ToReal(-2),ltet2,n2),kmul(nn,kmul(rm2,Ro221)),kmadd(nn,kmsub(kmsub(ksub(n1,ltet1),rm2,kmul(ltet2,rm1)),Ro121,kmul(Ro231,kmadd(ltet3,rm2,kmul(ksub(ltet2,n2),rm3)))),kmadd(nn,kmsub(kmsub(ksub(n1,ltet1),rm3,kmul(ltet3,rm1)),Ro131,kmul(Ro311,kmadd(ksub(ltet3,n3),rm1,kmul(ltet1,rm3)))),kmadd(nn,kmsub(kmadd(ToReal(-2),ltet1,n1),kmul(rm1,Ro111),kmadd(kmadd(ksub(ltet3,n3),rm2,kmul(ltet2,rm3)),Ro321,kmul(Ro211,kmadd(ksub(ltet2,n2),rm1,kmul(ltet1,rm2))))),kmul(kmadd(ToReal(-2),ltet3,n3),kmul(nn,kmul(rm3,Ro331)))))))))),kmadd(n2,kmadd(ltet1,kmsub(kmadd(n2,R4p1212,kmul(n3,R4p1213)),rm1,kmul(rm3,kmadd(n3,R4p1323,kmul(n2,R4p1223)))),kmadd(ltet3,knmsub(kmadd(n2,R4p1223,kmul(n1,R4p1213)),rm1,kmadd(n1,kmul(R4p1323,rm3),kmul(n2,kmul(R4p2323,rm3)))),kmadd(ltet2,kmadd(kmsub(n3,R4p1223,kmul(n1,R4p1212)),rm1,kmul(rm3,kmsub(n1,R4p1223,kmul(n3,R4p2323)))),kmul(nn,kmadd(kmadd(ToReal(-2),ltet1,n1),kmul(rm1,Ro112),knmsub(kmadd(ltet2,rm1,kmul(rm2,ksub(ltet1,n1))),Ro122,knmsub(kmadd(ltet3,rm1,kmul(rm3,ksub(ltet1,n1))),Ro132,kmadd(kmsub(ksub(n2,ltet2),rm1,kmul(ltet1,rm2)),Ro212,kmadd(kmadd(ToReal(-2),ltet2,n2),kmul(rm2,Ro222),kmadd(kmsub(ksub(n2,ltet2),rm3,kmul(ltet3,rm2)),Ro232,knmsub(kmadd(ksub(ltet3,n3),rm1,kmul(ltet1,rm3)),Ro312,kmsub(kmadd(ToReal(-2),ltet3,n3),kmul(rm3,Ro332),kmul(Ro322,kmadd(ksub(ltet3,n3),rm2,kmul(ltet2,rm3))))))))))))))),kmsub(n3,kmadd(ltet1,kmadd(n2,kmul(R4p1213,rm1),kmadd(n3,kmul(R4p1313,rm1),kmadd(n2,kmul(R4p1223,rm2),kmul(n3,kmul(R4p1323,rm2))))),knmsub(ltet3,kmadd(n1,kmul(R4p1313,rm1),kmadd(n2,kmul(R4p1323,rm1),kmadd(n1,kmul(R4p1323,rm2),kmul(n2,kmul(R4p2323,rm2))))),kmadd(ltet2,kmsub(n3,kmadd(R4p1323,rm1,kmul(R4p2323,rm2)),kmul(n1,kmadd(R4p1213,rm1,kmul(R4p1223,rm2)))),kmadd(kmadd(ToReal(-2),ltet2,n2),kmul(nn,kmul(rm2,Ro223)),kmadd(nn,kmsub(kmadd(ToReal(-2),ltet1,n1),kmul(rm1,Ro113),kmadd(kmadd(ltet2,rm1,kmul(rm2,ksub(ltet1,n1))),Ro123,kmadd(kmadd(ltet3,rm1,kmul(rm3,ksub(ltet1,n1))),Ro133,kmadd(kmadd(ltet3,rm2,kmul(rm3,ksub(ltet2,n2))),Ro233,kmul(Ro213,kmadd(ksub(ltet2,n2),rm1,kmul(ltet1,rm2))))))),kmadd(nn,kmul(Ro313,kmsub(ksub(n3,ltet3),rm1,kmul(ltet1,rm3))),kmadd(nn,kmul(Ro323,kmsub(ksub(n3,ltet3),rm2,kmul(ltet2,rm3))),kmul(kmadd(ToReal(-2),ltet3,n3),kmul(nn,kmul(rm3,Ro333)))))))))),kmul(kmadd(ksub(n1,ltet1),kmadd(rm1,Rojo11,kmadd(rm2,Rojo12,kmul(rm3,Rojo13))),kmadd(ksub(n2,ltet2),kmadd(rm1,Rojo21,kmadd(rm2,Rojo22,kmul(rm3,Rojo23))),kmul(kmadd(rm1,Rojo31,kmadd(rm2,Rojo32,kmul(rm3,Rojo33))),ksub(n3,ltet3)))),kmul(nn,nn)))));
    
    CCTK_REAL_VEC Psi3iL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(n1,knmsub(im2,kmadd(kmsub(ltet2,n1,kmul(ltet1,n2)),R4p1212,kmadd(kmsub(ltet3,n1,kmul(ltet1,n3)),R4p1213,kmul(R4p1223,kmsub(ltet3,n2,kmul(ltet2,n3))))),kmsub(im1,kmul(kmsub(ToReal(2),ltet1,n1),kmul(nn,Ro111)),kmul(im3,kmadd(kmsub(ltet2,n1,kmul(ltet1,n2)),R4p1213,kmadd(kmsub(ltet3,n1,kmul(ltet1,n3)),R4p1313,kmul(R4p1323,kmsub(ltet3,n2,kmul(ltet2,n3)))))))),kmadd(n2,kmsub(im1,kmadd(kmsub(ltet2,n1,kmul(ltet1,n2)),R4p1212,kmadd(kmsub(ltet3,n1,kmul(ltet1,n3)),R4p1213,kmadd(kmsub(ltet3,n2,kmul(ltet2,n3)),R4p1223,kmul(kmsub(ToReal(2),ltet1,n1),kmul(nn,Ro112))))),kmul(im3,kmadd(kmsub(ltet2,n1,kmul(ltet1,n2)),R4p1223,kmadd(kmsub(ltet3,n1,kmul(ltet1,n3)),R4p1323,kmul(R4p2323,kmsub(ltet3,n2,kmul(ltet2,n3))))))),kmadd(n3,kmsub(im1,kmadd(kmsub(ltet2,n1,kmul(ltet1,n2)),R4p1213,kmadd(kmsub(ltet3,n1,kmul(ltet1,n3)),R4p1313,kmadd(kmsub(ltet3,n2,kmul(ltet2,n3)),R4p1323,kmul(kmsub(ToReal(2),ltet1,n1),kmul(nn,Ro113))))),kmul(im2,kmadd(kmsub(ltet1,n2,kmul(ltet2,n1)),R4p1223,kmadd(kmsub(ltet1,n3,kmul(ltet3,n1)),R4p1323,kmul(R4p2323,kmsub(ltet2,n3,kmul(ltet3,n2))))))),kmadd(kmadd(im1,ltet2,kmul(im2,ksub(ltet1,n1))),kmul(n1,kmul(nn,Ro121)),kmadd(kmadd(im1,ltet2,kmul(im2,ksub(ltet1,n1))),kmul(n2,kmul(nn,Ro122)),kmadd(kmadd(im1,ltet2,kmul(im2,ksub(ltet1,n1))),kmul(n3,kmul(nn,Ro123)),kmadd(kmadd(im1,ltet3,kmul(im3,ksub(ltet1,n1))),kmul(n1,kmul(nn,Ro131)),kmadd(kmadd(im1,ltet3,kmul(im3,ksub(ltet1,n1))),kmul(n2,kmul(nn,Ro132)),kmadd(kmadd(im1,ltet3,kmul(im3,ksub(ltet1,n1))),kmul(n3,kmul(nn,Ro133)),kmadd(n1,kmul(kmadd(im2,ltet1,kmul(im1,ksub(ltet2,n2))),kmul(nn,Ro211)),kmadd(kmadd(im2,ltet1,kmul(im1,ksub(ltet2,n2))),kmul(n2,kmul(nn,Ro212)),kmadd(kmadd(im2,ltet1,kmul(im1,ksub(ltet2,n2))),kmul(n3,kmul(nn,Ro213)),kmadd(im2,kmul(n1,kmul(kmsub(ToReal(2),ltet2,n2),kmul(nn,Ro221))),kmadd(im2,kmul(kmsub(ToReal(2),ltet2,n2),kmul(n2,kmul(nn,Ro222))),kmadd(im2,kmul(kmsub(ToReal(2),ltet2,n2),kmul(n3,kmul(nn,Ro223))),kmadd(n1,kmul(kmadd(im2,ltet3,kmul(im3,ksub(ltet2,n2))),kmul(nn,Ro231)),kmadd(kmadd(im2,ltet3,kmul(im3,ksub(ltet2,n2))),kmul(n2,kmul(nn,Ro232)),kmadd(kmadd(im2,ltet3,kmul(im3,ksub(ltet2,n2))),kmul(n3,kmul(nn,Ro233)),kmadd(n1,kmul(kmadd(im3,ltet1,kmul(im1,ksub(ltet3,n3))),kmul(nn,Ro311)),kmadd(n2,kmul(kmadd(im3,ltet1,kmul(im1,ksub(ltet3,n3))),kmul(nn,Ro312)),kmadd(kmadd(im3,ltet1,kmul(im1,ksub(ltet3,n3))),kmul(n3,kmul(nn,Ro313)),kmadd(n1,kmul(kmadd(im3,ltet2,kmul(im2,ksub(ltet3,n3))),kmul(nn,Ro321)),kmadd(n2,kmul(kmadd(im3,ltet2,kmul(im2,ksub(ltet3,n3))),kmul(nn,Ro322)),kmadd(kmadd(im3,ltet2,kmul(im2,ksub(ltet3,n3))),kmul(n3,kmul(nn,Ro323)),kmadd(im3,kmul(n1,kmul(kmsub(ToReal(2),ltet3,n3),kmul(nn,Ro331))),kmadd(im3,kmul(n2,kmul(kmsub(ToReal(2),ltet3,n3),kmul(nn,Ro332))),kmadd(im3,kmul(kmsub(ToReal(2),ltet3,n3),kmul(n3,kmul(nn,Ro333))),kmul(kmadd(ksub(n1,ltet1),kmadd(im1,Rojo11,kmadd(im2,Rojo12,kmul(im3,Rojo13))),kmadd(ksub(n2,ltet2),kmadd(im1,Rojo21,kmadd(im2,Rojo22,kmul(im3,Rojo23))),kmul(kmadd(im1,Rojo31,kmadd(im2,Rojo32,kmul(im3,Rojo33))),ksub(n3,ltet3)))),kmul(nn,nn)))))))))))))))))))))))))))));
    
    CCTK_REAL_VEC Psi2rL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(kmadd(ltet2,kmsub(n1,R4p1212,kmul(n3,R4p1223)),kmul(ltet3,kmsub(n1,R4p1213,kmul(n3,R4p1323)))),kmadd(im1,im2,kmul(rm1,rm2)),kmadd(kmsub(ltet1,kmadd(n2,R4p1212,kmul(n3,R4p1213)),kmul(ltet3,kmadd(n2,R4p1223,kmul(n3,R4p1323)))),kmadd(im1,im2,kmul(rm1,rm2)),kmadd(kmadd(ltet2,kmul(n1,R4p1213),kmadd(ltet2,kmul(n2,R4p1223),kmadd(ltet3,kmul(n1,R4p1313),kmul(ltet3,kmul(n2,R4p1323))))),kmadd(im1,im3,kmul(rm1,rm3)),kmadd(kmadd(ltet1,kmul(n2,R4p1213),kmadd(ltet2,kmul(n2,R4p1223),kmadd(ltet1,kmul(n3,R4p1313),kmul(ltet2,kmul(n3,R4p1323))))),kmadd(im1,im3,kmul(rm1,rm3)),kmadd(knmsub(n1,kmadd(ltet2,R4p1223,kmul(ltet1,R4p1213)),kmadd(ltet1,knmsub(n1,R4p1213,kmsub(n3,R4p1323,kmul(n2,R4p1223))),kmadd(ltet2,kmul(n3,R4p2323),kmul(ltet3,kmadd(n1,R4p1323,kmul(n2,R4p2323)))))),kmadd(im2,im3,kmul(rm2,rm3)),knmsub(kmadd(ltet2,kmul(n2,R4p1212),kmadd(ltet3,kmul(n2,R4p1213),kmadd(ltet2,kmul(n3,R4p1213),kmul(ltet3,kmul(n3,R4p1313))))),kmadd(im1,im1,kmul(rm1,rm1)),knmsub(kmadd(n1,kmsub(ltet1,R4p1212,kmul(ltet3,R4p1223)),kmul(n3,kmsub(ltet3,R4p2323,kmul(ltet1,R4p1223)))),kmadd(im2,im2,kmul(rm2,rm2)),knmsub(kmadd(ltet1,kmul(n1,R4p1313),kmadd(ltet2,kmul(n1,R4p1323),kmadd(ltet1,kmul(n2,R4p1323),kmul(ltet2,kmul(n2,R4p2323))))),kmadd(im3,im3,kmul(rm3,rm3)),kmsub(nn,kmadd(ksub(n2,ltet2),kmul(Ro122,kmadd(im1,im2,kmul(rm1,rm2))),kmadd(kmadd(im1,kmsub(im2,n3,kmul(im3,ltet2)),kmul(rm1,kmsub(n3,rm2,kmul(ltet2,rm3)))),Ro123,kmadd(kmadd(im1,kmsub(im3,n2,kmul(im2,ltet3)),kmul(rm1,kmsub(n2,rm3,kmul(ltet3,rm2)))),Ro132,kmadd(ksub(n3,ltet3),kmul(Ro133,kmadd(im1,im3,kmul(rm1,rm3))),kmadd(ksub(n1,ltet1),kmul(Ro211,kmadd(im1,im2,kmul(rm1,rm2))),kmadd(kmadd(im2,kmsub(im1,n3,kmul(im3,ltet1)),kmul(rm2,kmsub(n3,rm1,kmul(ltet1,rm3)))),Ro213,kmadd(kmadd(im2,kmsub(im3,n1,kmul(im1,ltet3)),kmul(rm2,kmsub(n1,rm3,kmul(ltet3,rm1)))),Ro231,kmadd(ksub(n3,ltet3),kmul(Ro233,kmadd(im2,im3,kmul(rm2,rm3))),kmadd(ksub(n1,ltet1),kmul(Ro311,kmadd(im1,im3,kmul(rm1,rm3))),kmadd(kmadd(im3,kmsub(im1,n2,kmul(im2,ltet1)),kmul(rm3,kmsub(n2,rm1,kmul(ltet1,rm2)))),Ro312,kmadd(kmadd(im3,kmsub(im2,n1,kmul(im1,ltet2)),kmul(rm3,kmsub(n1,rm2,kmul(ltet2,rm1)))),Ro321,kmadd(ksub(n2,ltet2),kmul(Ro322,kmadd(im2,im3,kmul(rm2,rm3))),kmadd(ksub(n1,ltet1),kmul(Ro111,kmadd(im1,im1,kmul(rm1,rm1))),kmadd(Ro121,kmsub(n1,kmadd(im1,im2,kmul(rm1,rm2)),kmul(ltet2,kmadd(im1,im1,kmul(rm1,rm1)))),kmadd(Ro131,kmsub(n1,kmadd(im1,im3,kmul(rm1,rm3)),kmul(ltet3,kmadd(im1,im1,kmul(rm1,rm1)))),kmadd(Ro112,kmsub(n2,kmadd(im1,im1,kmul(rm1,rm1)),kmul(ltet1,kmadd(im1,im2,kmul(rm1,rm2)))),kmadd(Ro113,kmsub(n3,kmadd(im1,im1,kmul(rm1,rm1)),kmul(ltet1,kmadd(im1,im3,kmul(rm1,rm3)))),kmadd(ksub(n2,ltet2),kmul(Ro222,kmadd(im2,im2,kmul(rm2,rm2))),kmadd(Ro212,kmsub(n2,kmadd(im1,im2,kmul(rm1,rm2)),kmul(ltet1,kmadd(im2,im2,kmul(rm2,rm2)))),kmadd(Ro232,kmsub(n2,kmadd(im2,im3,kmul(rm2,rm3)),kmul(ltet3,kmadd(im2,im2,kmul(rm2,rm2)))),kmadd(Ro221,kmsub(n1,kmadd(im2,im2,kmul(rm2,rm2)),kmul(ltet2,kmadd(im1,im2,kmul(rm1,rm2)))),kmadd(Ro223,kmsub(n3,kmadd(im2,im2,kmul(rm2,rm2)),kmul(ltet2,kmadd(im2,im3,kmul(rm2,rm3)))),kmadd(ksub(n3,ltet3),kmul(Ro333,kmadd(im3,im3,kmul(rm3,rm3))),kmadd(Ro313,kmsub(n3,kmadd(im1,im3,kmul(rm1,rm3)),kmul(ltet1,kmadd(im3,im3,kmul(rm3,rm3)))),kmadd(Ro323,kmsub(n3,kmadd(im2,im3,kmul(rm2,rm3)),kmul(ltet2,kmadd(im3,im3,kmul(rm3,rm3)))),kmadd(Ro331,kmsub(n1,kmadd(im3,im3,kmul(rm3,rm3)),kmul(ltet3,kmadd(im1,im3,kmul(rm1,rm3)))),kmul(Ro332,kmsub(n2,kmadd(im3,im3,kmul(rm3,rm3)),kmul(ltet3,kmadd(im2,im3,kmul(rm2,rm3))))))))))))))))))))))))))))))),kmul(kmadd(im2,kmul(im3,Rojo23),kmadd(rm2,kmul(rm3,Rojo23),kmadd(im1,kmadd(im2,kadd(Rojo12,Rojo21),kmul(im3,kadd(Rojo13,Rojo31))),kmadd(rm1,kmadd(rm2,kadd(Rojo12,Rojo21),kmul(rm3,kadd(Rojo13,Rojo31))),kmadd(im2,kmul(im3,Rojo32),kmadd(rm2,kmul(rm3,Rojo32),kmadd(Rojo11,kmul(im1,im1),kmadd(Rojo22,kmul(im2,im2),kmadd(Rojo33,kmul(im3,im3),kmadd(Rojo11,kmul(rm1,rm1),kmadd(Rojo22,kmul(rm2,rm2),kmul(Rojo33,kmul(rm3,rm3))))))))))))),kmul(nn,nn)))))))))));
    
    CCTK_REAL_VEC Psi2iL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(knmsub(n1,kmadd(ltet3,R4p1213,kmul(ltet2,R4p1212)),kmadd(ltet1,kmadd(n2,R4p1212,kmul(n3,R4p1213)),kmul(R4p1223,kmsub(ltet2,n3,kmul(ltet3,n2))))),kmsub(im2,rm1,kmul(im1,rm2)),kmadd(kmadd(kmsub(ltet1,n2,kmul(ltet2,n1)),R4p1213,kmadd(kmsub(ltet1,n3,kmul(ltet3,n1)),R4p1313,kmul(R4p1323,kmsub(ltet2,n3,kmul(ltet3,n2))))),kmsub(im3,rm1,kmul(im1,rm3)),kmadd(knmsub(n1,kmadd(ltet2,R4p1223,kmul(ltet1,R4p1213)),kmadd(ltet1,kmadd(n1,R4p1213,kmul(n2,R4p1223)),kmadd(ltet1,kmul(n3,R4p1323),kmsub(ltet2,kmul(n3,R4p2323),kmul(ltet3,kmadd(n1,R4p1323,kmul(n2,R4p2323))))))),kmsub(im3,rm2,kmul(im2,rm3)),kmsub(kmadd(im1,kmadd(rm2,ksub(Rojo21,Rojo12),kmul(rm3,ksub(Rojo31,Rojo13))),kmadd(im3,kmadd(rm1,ksub(Rojo13,Rojo31),kmul(rm2,ksub(Rojo23,Rojo32))),kmul(im2,kmadd(rm1,ksub(Rojo12,Rojo21),kmul(rm3,ksub(Rojo32,Rojo23)))))),kmul(nn,nn),kmul(nn,kmadd(im1,knmsub(rm2,kmadd(kadd(ltet2,n2),Ro122,kmadd(ltet3,Ro132,kmul(n3,Ro123))),knmsub(rm3,kmadd(ltet2,Ro123,kmadd(kadd(ltet3,n3),Ro133,kmul(n2,Ro132))),kmadd(n2,kmul(rm2,Ro212),kmadd(n3,kmul(rm2,Ro213),kmadd(ltet2,kmul(rm2,Ro221),kmadd(ltet3,kmul(rm2,Ro231),kmadd(ltet1,kmadd(rm2,ksub(Ro211,Ro112),kmul(rm3,ksub(Ro311,Ro113))),kmadd(n1,kmadd(rm2,ksub(Ro211,Ro121),kmul(rm3,ksub(Ro311,Ro131))),kmadd(n2,kmul(rm3,Ro312),kmadd(n3,kmul(rm3,Ro313),kmadd(ltet2,kmul(rm3,Ro321),kmul(ltet3,kmul(rm3,Ro331))))))))))))),kmadd(im2,kmadd(rm1,kmadd(n2,ksub(Ro122,Ro212),kmadd(n3,ksub(Ro123,Ro213),kmadd(ltet2,ksub(Ro122,Ro221),kmul(ltet3,ksub(Ro132,Ro231))))),knmsub(rm3,kmadd(ltet2,Ro223,kmadd(kadd(ltet3,n3),Ro233,kmul(n2,Ro232))),kmadd(ltet1,kmadd(rm1,ksub(Ro112,Ro211),kmul(rm3,ksub(Ro312,Ro213))),kmadd(n1,kmadd(rm1,ksub(Ro121,Ro211),kmul(rm3,ksub(Ro321,Ro231))),kmadd(ltet2,kmul(rm3,Ro322),kmadd(n2,kmul(rm3,Ro322),kmadd(n3,kmul(rm3,Ro323),kmul(ltet3,kmul(rm3,Ro332))))))))),kmul(im3,kmadd(ltet1,kmadd(rm1,ksub(Ro113,Ro311),kmul(rm2,ksub(Ro213,Ro312))),kmadd(ltet2,kmadd(rm1,ksub(Ro123,Ro321),kmul(rm2,ksub(Ro223,Ro322))),kmadd(rm1,kmadd(kadd(ltet3,n3),Ro133,kmadd(n1,ksub(Ro131,Ro311),kmsub(n2,ksub(Ro132,Ro312),kmadd(ltet3,Ro331,kmul(n3,Ro313))))),kmul(rm2,kmadd(kadd(ltet3,n3),Ro233,kmadd(n1,ksub(Ro231,Ro321),kmsub(n2,ksub(Ro232,Ro322),kmadd(ltet3,Ro332,kmul(n3,Ro323)))))))))))))))));
    
    CCTK_REAL_VEC Psi1rL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(ltet1,kmadd(ltet1,kmadd(n2,kmul(R4p1212,rm2),kmadd(n3,kmul(R4p1213,rm2),kmadd(n2,kmul(R4p1213,rm3),kmul(n3,kmul(R4p1313,rm3))))),knmsub(ltet3,kmadd(n1,kmul(R4p1213,rm2),kmadd(n2,kmul(R4p1223,rm2),kmadd(n1,kmul(R4p1313,rm3),kmul(n2,kmul(R4p1323,rm3))))),kmadd(ltet2,kmsub(n3,kmadd(R4p1223,rm2,kmul(R4p1323,rm3)),kmul(n1,kmadd(R4p1212,rm2,kmul(R4p1213,rm3)))),kmadd(nn,kmul(Ro121,kmsub(ksub(ltet1,n1),rm2,kmul(n2,rm1))),kmadd(nn,kmul(Ro131,kmsub(ksub(ltet1,n1),rm3,kmul(n3,rm1))),kmadd(nn,kmul(Ro211,kmsub(ksub(ltet2,n2),rm1,kmul(n1,rm2))),kmadd(kmadd(ToReal(-2),n2,ltet2),kmul(nn,kmul(rm2,Ro221)),kmadd(nn,kmsub(kmadd(ToReal(-2),n1,ltet1),kmul(rm1,Ro111),kmul(Ro231,kmadd(n3,rm2,kmul(ksub(n2,ltet2),rm3)))),kmadd(nn,kmul(Ro311,kmsub(ksub(ltet3,n3),rm1,kmul(n1,rm3))),kmadd(nn,kmul(Ro321,kmsub(ksub(ltet3,n3),rm2,kmul(n2,rm3))),kmul(kmadd(ToReal(-2),n3,ltet3),kmul(nn,kmul(rm3,Ro331))))))))))))),kmadd(ltet2,kmadd(ltet1,knmsub(kmadd(n3,R4p1213,kmul(n2,R4p1212)),rm1,kmadd(n2,kmul(R4p1223,rm3),kmul(n3,kmul(R4p1323,rm3)))),kmadd(ltet3,kmsub(kmadd(n1,R4p1213,kmul(n2,R4p1223)),rm1,kmul(rm3,kmadd(n2,R4p2323,kmul(n1,R4p1323)))),kmadd(ltet2,kmadd(kmsub(n1,R4p1212,kmul(n3,R4p1223)),rm1,kmul(rm3,kmsub(n3,R4p2323,kmul(n1,R4p1223)))),kmul(nn,kmadd(kmadd(ToReal(-2),n1,ltet1),kmul(rm1,Ro112),knmsub(kmadd(n2,rm1,kmul(rm2,ksub(n1,ltet1))),Ro122,knmsub(kmadd(n3,rm1,kmul(rm3,ksub(n1,ltet1))),Ro132,kmadd(kmsub(ksub(ltet2,n2),rm1,kmul(n1,rm2)),Ro212,kmadd(kmadd(ToReal(-2),n2,ltet2),kmul(rm2,Ro222),kmadd(kmsub(ksub(ltet2,n2),rm3,kmul(n3,rm2)),Ro232,kmadd(kmsub(ksub(ltet3,n3),rm1,kmul(n1,rm3)),Ro312,kmadd(kmsub(ksub(ltet3,n3),rm2,kmul(n2,rm3)),Ro322,kmul(kmadd(ToReal(-2),n3,ltet3),kmul(rm3,Ro332)))))))))))))),kmsub(ltet3,knmsub(ltet1,kmadd(n2,kmul(R4p1213,rm1),kmadd(n3,kmul(R4p1313,rm1),kmadd(n2,kmul(R4p1223,rm2),kmul(n3,kmul(R4p1323,rm2))))),kmadd(ltet3,kmadd(n1,kmul(R4p1313,rm1),kmadd(n2,kmul(R4p1323,rm1),kmadd(n1,kmul(R4p1323,rm2),kmul(n2,kmul(R4p2323,rm2))))),kmadd(ltet2,kmsub(n1,kmadd(R4p1213,rm1,kmul(R4p1223,rm2)),kmul(n3,kmadd(R4p1323,rm1,kmul(R4p2323,rm2)))),kmadd(nn,kmul(Ro213,kmsub(ksub(ltet2,n2),rm1,kmul(n1,rm2))),kmadd(kmadd(ToReal(-2),n2,ltet2),kmul(nn,kmul(rm2,Ro223)),kmadd(nn,kmsub(kmadd(ToReal(-2),n1,ltet1),kmul(rm1,Ro113),kmadd(kmadd(n2,rm1,kmul(rm2,ksub(n1,ltet1))),Ro123,kmadd(kmadd(n3,rm2,kmul(rm3,ksub(n2,ltet2))),Ro233,kmul(Ro133,kmadd(n3,rm1,kmul(ksub(n1,ltet1),rm3)))))),kmadd(nn,kmul(Ro313,kmsub(ksub(ltet3,n3),rm1,kmul(n1,rm3))),kmadd(nn,kmul(Ro323,kmsub(ksub(ltet3,n3),rm2,kmul(n2,rm3))),kmul(kmadd(ToReal(-2),n3,ltet3),kmul(nn,kmul(rm3,Ro333))))))))))),kmul(kmadd(ksub(ltet1,n1),kmadd(rm1,Rojo11,kmadd(rm2,Rojo12,kmul(rm3,Rojo13))),kmadd(rm1,kmadd(ksub(ltet2,n2),Rojo21,kmul(Rojo31,ksub(ltet3,n3))),kmadd(rm2,kmadd(ksub(ltet2,n2),Rojo22,kmul(Rojo32,ksub(ltet3,n3))),kmul(rm3,kmadd(ksub(ltet2,n2),Rojo23,kmul(Rojo33,ksub(ltet3,n3))))))),kmul(nn,nn)))));
    
    CCTK_REAL_VEC Psi1iL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(ltet1,kmadd(im2,kmadd(kmsub(ltet1,n2,kmul(ltet2,n1)),R4p1212,kmadd(kmsub(ltet1,n3,kmul(ltet3,n1)),R4p1213,kmadd(kmsub(ltet2,n3,kmul(ltet3,n2)),R4p1223,kmul(kmadd(ToReal(-2),n2,ltet2),kmul(nn,Ro221))))),kmadd(nn,kmadd(im1,kmul(Ro111,kmadd(ToReal(-2),n1,ltet1)),kmadd(kmsub(im2,ksub(ltet1,n1),kmul(im1,n2)),Ro121,kmadd(kmsub(im3,ksub(ltet1,n1),kmul(im1,n3)),Ro131,knmsub(kmadd(im2,n1,kmul(im1,ksub(n2,ltet2))),Ro211,kmsub(kmsub(im3,ksub(ltet2,n2),kmul(im2,n3)),Ro231,kmadd(kmadd(im3,n2,kmul(im2,ksub(n3,ltet3))),Ro321,kmul(Ro311,kmadd(im3,n1,kmul(im1,ksub(n3,ltet3)))))))))),kmul(im3,kmadd(kmsub(ltet1,n2,kmul(ltet2,n1)),R4p1213,kmadd(kmsub(ltet1,n3,kmul(ltet3,n1)),R4p1313,kmadd(kmsub(ltet2,n3,kmul(ltet3,n2)),R4p1323,kmul(kmadd(ToReal(-2),n3,ltet3),kmul(nn,Ro331)))))))),kmadd(ltet2,kmadd(im1,kmadd(kmsub(ltet2,n1,kmul(ltet1,n2)),R4p1212,kmadd(kmsub(ltet3,n1,kmul(ltet1,n3)),R4p1213,kmadd(kmsub(ltet3,n2,kmul(ltet2,n3)),R4p1223,kmul(kmadd(ToReal(-2),n1,ltet1),kmul(nn,Ro112))))),kmadd(nn,kmadd(kmsub(im2,ksub(ltet1,n1),kmul(im1,n2)),Ro122,kmadd(kmsub(im3,ksub(ltet1,n1),kmul(im1,n3)),Ro132,kmadd(kmsub(im1,ksub(ltet2,n2),kmul(im2,n1)),Ro212,kmadd(im2,kmul(Ro222,kmadd(ToReal(-2),n2,ltet2)),kmsub(kmsub(im3,ksub(ltet2,n2),kmul(im2,n3)),Ro232,kmadd(kmadd(im3,n2,kmul(im2,ksub(n3,ltet3))),Ro322,kmul(Ro312,kmadd(im3,n1,kmul(im1,ksub(n3,ltet3)))))))))),kmul(im3,kmadd(kmsub(ltet1,n2,kmul(ltet2,n1)),R4p1223,kmadd(kmsub(ltet1,n3,kmul(ltet3,n1)),R4p1323,kmadd(kmsub(ltet2,n3,kmul(ltet3,n2)),R4p2323,kmul(kmadd(ToReal(-2),n3,ltet3),kmul(nn,Ro332)))))))),kmsub(ltet3,kmadd(im1,kmadd(kmsub(ltet2,n1,kmul(ltet1,n2)),R4p1213,kmadd(kmsub(ltet3,n1,kmul(ltet1,n3)),R4p1313,kmadd(kmsub(ltet3,n2,kmul(ltet2,n3)),R4p1323,kmul(kmadd(ToReal(-2),n1,ltet1),kmul(nn,Ro113))))),kmadd(im2,kmadd(kmsub(ltet2,n1,kmul(ltet1,n2)),R4p1223,kmadd(kmsub(ltet3,n1,kmul(ltet1,n3)),R4p1323,kmadd(kmsub(ltet3,n2,kmul(ltet2,n3)),R4p2323,kmul(kmadd(ToReal(-2),n2,ltet2),kmul(nn,Ro223))))),kmul(nn,kmadd(kmsub(im2,ksub(ltet1,n1),kmul(im1,n2)),Ro123,kmadd(kmsub(im3,ksub(ltet1,n1),kmul(im1,n3)),Ro133,knmsub(kmadd(im2,n1,kmul(im1,ksub(n2,ltet2))),Ro213,kmadd(kmsub(im3,ksub(ltet2,n2),kmul(im2,n3)),Ro233,kmadd(kmsub(im1,ksub(ltet3,n3),kmul(im3,n1)),Ro313,kmadd(kmsub(im2,ksub(ltet3,n3),kmul(im3,n2)),Ro323,kmul(im3,kmul(Ro333,kmadd(ToReal(-2),n3,ltet3)))))))))))),kmul(kmadd(im1,kmadd(ksub(ltet1,n1),Rojo11,kmadd(ksub(ltet2,n2),Rojo21,kmul(Rojo31,ksub(ltet3,n3)))),kmadd(im2,kmadd(ksub(ltet1,n1),Rojo12,kmadd(ksub(ltet2,n2),Rojo22,kmul(Rojo32,ksub(ltet3,n3)))),kmul(im3,kmadd(ksub(ltet1,n1),Rojo13,kmadd(ksub(ltet2,n2),Rojo23,kmul(Rojo33,ksub(ltet3,n3))))))),kmul(nn,nn)))));
    
    CCTK_REAL_VEC Psi0rL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(ToReal(2),kmul(kmsub(im1,im2,kmul(rm1,rm2)),kmsub(ltet1,kmadd(ltet2,R4p1212,kmul(ltet3,R4p1213)),kmul(ltet3,kmadd(ltet2,R4p1223,kmul(ltet3,R4p1323))))),kmadd(ToReal(2),kmul(kmadd(ltet1,kmsub(ltet2,R4p1223,kmul(ltet3,R4p1323)),kmsub(R4p1213,kmul(ltet1,ltet1),kmul(ltet2,kmul(ltet3,R4p2323)))),kmsub(rm2,rm3,kmul(im2,im3))),kmadd(ToReal(2),kmul(kmadd(ltet1,kmul(ltet2,R4p1213),kmadd(ltet1,kmul(ltet3,R4p1313),kmadd(ltet2,kmul(ltet3,R4p1323),kmul(R4p1223,kmul(ltet2,ltet2))))),kmsub(im1,im3,kmul(rm1,rm3))),knmsub(kmadd(ToReal(2),kmul(ltet2,kmul(ltet3,R4p1213)),kmadd(R4p1212,kmul(ltet2,ltet2),kmul(R4p1313,kmul(ltet3,ltet3)))),kmsub(im1,im1,kmul(rm1,rm1)),knmsub(kmadd(ToReal(-2),kmul(ltet1,kmul(ltet3,R4p1223)),kmadd(R4p1212,kmul(ltet1,ltet1),kmul(R4p2323,kmul(ltet3,ltet3)))),kmsub(im2,im2,kmul(rm2,rm2)),knmsub(kmul(nn,nn),kmadd(kmsub(im2,im3,kmul(rm2,rm3)),Rojo23,kmadd(im1,kmadd(im2,kadd(Rojo12,Rojo21),kmul(im3,kadd(Rojo13,Rojo31))),knmsub(rm1,kmadd(rm2,kadd(Rojo12,Rojo21),kmul(rm3,kadd(Rojo13,Rojo31))),kmadd(kmsub(im2,im3,kmul(rm2,rm3)),Rojo32,kmadd(Rojo11,kmsub(im1,im1,kmul(rm1,rm1)),kmadd(Rojo22,kmsub(im2,im2,kmul(rm2,rm2)),kmul(Rojo33,kmsub(im3,im3,kmul(rm3,rm3))))))))),kmsub(ToReal(2),kmul(nn,kmadd(kmsub(rm1,rm2,kmul(im1,im2)),kmadd(ltet1,Ro112,kmadd(ltet2,Ro122,kmul(ltet3,Ro132))),kmadd(kmsub(rm1,rm3,kmul(im1,im3)),kmadd(ltet1,Ro113,kmadd(ltet2,Ro123,kmul(ltet3,Ro133))),kmadd(kmsub(rm1,rm2,kmul(im1,im2)),kmadd(ltet1,Ro211,kmadd(ltet2,Ro221,kmul(ltet3,Ro231))),kmadd(kmsub(rm2,rm3,kmul(im2,im3)),kmadd(ltet1,Ro213,kmadd(ltet2,Ro223,kmul(ltet3,Ro233))),kmadd(kmsub(rm1,rm3,kmul(im1,im3)),kmadd(ltet1,Ro311,kmadd(ltet2,Ro321,kmul(ltet3,Ro331))),kmadd(kmsub(rm2,rm3,kmul(im2,im3)),kmadd(ltet1,Ro312,kmadd(ltet2,Ro322,kmul(ltet3,Ro332))),kmadd(kmadd(ltet1,Ro111,kmadd(ltet2,Ro121,kmul(ltet3,Ro131))),kmsub(rm1,rm1,kmul(im1,im1)),kmadd(kmadd(ltet1,Ro212,kmadd(ltet2,Ro222,kmul(ltet3,Ro232))),kmsub(rm2,rm2,kmul(im2,im2)),kmul(kmadd(ltet1,Ro313,kmadd(ltet2,Ro323,kmul(ltet3,Ro333))),kmsub(rm3,rm3,kmul(im3,im3)))))))))))),kmul(kmadd(ToReal(2),kmul(ltet1,kmul(ltet2,R4p1323)),kmadd(R4p1313,kmul(ltet1,ltet1),kmul(R4p2323,kmul(ltet2,ltet2)))),kmsub(im3,im3,kmul(rm3,rm3))))))))));
    
    CCTK_REAL_VEC Psi0iL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(ToReal(-2),kmul(kmadd(im3,rm1,kmul(im1,rm3)),kmadd(ltet1,kmadd(ltet2,R4p1213,kmul(ltet3,R4p1313)),kmadd(ltet2,kmul(ltet3,R4p1323),kmul(R4p1223,kmul(ltet2,ltet2))))),kmadd(ToReal(2),kmadd(kmsub(ltet3,kmadd(ltet2,R4p1223,kmul(ltet3,R4p1323)),kmul(ltet1,kmadd(ltet2,R4p1212,kmul(ltet3,R4p1213)))),kmadd(im2,rm1,kmul(im1,rm2)),kmadd(nn,kmadd(kmadd(im2,rm1,kmul(im1,rm2)),kmadd(ltet1,kadd(Ro112,Ro211),kmadd(ltet2,kadd(Ro122,Ro221),kmul(ltet3,kadd(Ro132,Ro231)))),kmadd(kmadd(im3,rm1,kmul(im1,rm3)),kmadd(ltet1,kadd(Ro113,Ro311),kmadd(ltet2,kadd(Ro123,Ro321),kmul(ltet3,kadd(Ro133,Ro331)))),kmadd(kmadd(im3,rm2,kmul(im2,rm3)),kmadd(ltet1,kadd(Ro213,Ro312),kmadd(ltet2,kadd(Ro223,Ro322),kmul(ltet3,kadd(Ro233,Ro332)))),kmul(kmadd(im1,kmul(rm1,kmadd(ltet1,Ro111,kmadd(ltet2,Ro121,kmul(ltet3,Ro131)))),kmadd(im2,kmul(rm2,kmadd(ltet1,Ro212,kmadd(ltet2,Ro222,kmul(ltet3,Ro232)))),kmul(im3,kmul(rm3,kmadd(ltet1,Ro313,kmadd(ltet2,Ro323,kmul(ltet3,Ro333))))))),ToReal(2))))),kmadd(kmadd(im3,rm2,kmul(im2,rm3)),kmadd(ltet1,kmsub(ltet2,R4p1223,kmul(ltet3,R4p1323)),kmsub(R4p1213,kmul(ltet1,ltet1),kmul(ltet2,kmul(ltet3,R4p2323)))),kmadd(im3,kmul(rm3,kmadd(ToReal(2),kmul(ltet1,kmul(ltet2,R4p1323)),kmadd(R4p1313,kmul(ltet1,ltet1),kmul(R4p2323,kmul(ltet2,ltet2))))),kmadd(im1,kmul(rm1,kmadd(ToReal(2),kmul(ltet2,kmul(ltet3,R4p1213)),kmadd(R4p1212,kmul(ltet2,ltet2),kmul(R4p1313,kmul(ltet3,ltet3))))),kmul(im2,kmul(rm2,kmadd(ToReal(-2),kmul(ltet1,kmul(ltet3,R4p1223)),kmadd(R4p1212,kmul(ltet1,ltet1),kmul(R4p2323,kmul(ltet3,ltet3))))))))))),kmul(kmadd(im1,kmadd(ToReal(2),kmul(rm1,Rojo11),kmadd(rm2,kadd(Rojo12,Rojo21),kmul(rm3,kadd(Rojo13,Rojo31)))),kmadd(im2,kmadd(rm1,kadd(Rojo12,Rojo21),kmadd(ToReal(2),kmul(rm2,Rojo22),kmul(rm3,kadd(Rojo23,Rojo32)))),kmul(im3,kmadd(rm1,kadd(Rojo13,Rojo31),kmadd(rm2,kadd(Rojo23,Rojo32),kmul(kmul(rm3,Rojo33),ToReal(2))))))),kmul(nn,nn))));
    /* Copy local copies back to grid functions */
    vec_store_partial_prepare(i,vecimin,vecimax);
    vec_store_nta_partial(Psi0i[index],Psi0iL);
    vec_store_nta_partial(Psi0r[index],Psi0rL);
    vec_store_nta_partial(Psi1i[index],Psi1iL);
    vec_store_nta_partial(Psi1r[index],Psi1rL);
    vec_store_nta_partial(Psi2i[index],Psi2iL);
    vec_store_nta_partial(Psi2r[index],Psi2rL);
    vec_store_nta_partial(Psi3i[index],Psi3iL);
    vec_store_nta_partial(Psi3r[index],Psi3rL);
    vec_store_nta_partial(Psi4i[index],Psi4iL);
    vec_store_nta_partial(Psi4r[index],Psi4rL);
  }
  CCTK_ENDLOOP3STR(WeylScal4_psis_calc_4th);
}
extern "C" void WeylScal4_psis_calc_4th(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering WeylScal4_psis_calc_4th_Body");
  }
  if (cctk_iteration % WeylScal4_psis_calc_4th_calc_every != WeylScal4_psis_calc_4th_calc_offset)
  {
    return;
  }
  
  const char* const groups[] = {
    "admbase::curv",
    "admbase::metric",
    "grid::coordinates",
    "WeylScal4::Psi0i_group",
    "WeylScal4::Psi0r_group",
    "WeylScal4::Psi1i_group",
    "WeylScal4::Psi1r_group",
    "WeylScal4::Psi2i_group",
    "WeylScal4::Psi2r_group",
    "WeylScal4::Psi3i_group",
    "WeylScal4::Psi3r_group",
    "WeylScal4::Psi4i_group",
    "WeylScal4::Psi4r_group"};
  AssertGroupStorage(cctkGH, "WeylScal4_psis_calc_4th", 13, groups);
  
  switch (fdOrder)
  {
    case 2:
    {
      EnsureStencilFits(cctkGH, "WeylScal4_psis_calc_4th", 2, 2, 2);
      break;
    }
    
    case 4:
    {
      EnsureStencilFits(cctkGH, "WeylScal4_psis_calc_4th", 2, 2, 2);
      break;
    }
    
    case 6:
    {
      EnsureStencilFits(cctkGH, "WeylScal4_psis_calc_4th", 2, 2, 2);
      break;
    }
    
    case 8:
    {
      EnsureStencilFits(cctkGH, "WeylScal4_psis_calc_4th", 2, 2, 2);
      break;
    }
    default:
      CCTK_BUILTIN_UNREACHABLE();
  }
  
  LoopOverInterior(cctkGH, WeylScal4_psis_calc_4th_Body);
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving WeylScal4_psis_calc_4th_Body");
  }
}

} // namespace WeylScal4
