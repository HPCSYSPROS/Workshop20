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


static void WeylScal4_invars_calc_4th_Body(const cGH* restrict const cctkGH, const int dir, const int face, const CCTK_REAL normal[3], const CCTK_REAL tangentA[3], const CCTK_REAL tangentB[3], const int imin[3], const int imax[3], const int n_subblock_gfs, CCTK_REAL* restrict const subblock_gfs[])
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
  CCTK_LOOP3STR(WeylScal4_invars_calc_4th,
    i,j,k, imin0,imin1,imin2, imax0,imax1,imax2,
    cctk_ash[0],cctk_ash[1],cctk_ash[2],
    vecimin,vecimax, CCTK_REAL_VEC_SIZE)
  {
    const ptrdiff_t index CCTK_ATTRIBUTE_UNUSED = di*i + dj*j + dk*k;
    /* Assign local copies of grid functions */
    
    CCTK_REAL_VEC Psi0iL CCTK_ATTRIBUTE_UNUSED = vec_load(Psi0i[index]);
    CCTK_REAL_VEC Psi0rL CCTK_ATTRIBUTE_UNUSED = vec_load(Psi0r[index]);
    CCTK_REAL_VEC Psi1iL CCTK_ATTRIBUTE_UNUSED = vec_load(Psi1i[index]);
    CCTK_REAL_VEC Psi1rL CCTK_ATTRIBUTE_UNUSED = vec_load(Psi1r[index]);
    CCTK_REAL_VEC Psi2iL CCTK_ATTRIBUTE_UNUSED = vec_load(Psi2i[index]);
    CCTK_REAL_VEC Psi2rL CCTK_ATTRIBUTE_UNUSED = vec_load(Psi2r[index]);
    CCTK_REAL_VEC Psi3iL CCTK_ATTRIBUTE_UNUSED = vec_load(Psi3i[index]);
    CCTK_REAL_VEC Psi3rL CCTK_ATTRIBUTE_UNUSED = vec_load(Psi3r[index]);
    CCTK_REAL_VEC Psi4iL CCTK_ATTRIBUTE_UNUSED = vec_load(Psi4i[index]);
    CCTK_REAL_VEC Psi4rL CCTK_ATTRIBUTE_UNUSED = vec_load(Psi4r[index]);
    
    
    /* Include user supplied include files */
    /* Precompute derivatives */
    
    switch (fdOrder)
    {
      case 2:
      {
        break;
      }
      
      case 4:
      {
        break;
      }
      
      case 6:
      {
        break;
      }
      
      case 8:
      {
        break;
      }
      default:
        CCTK_BUILTIN_UNREACHABLE();
    }
    /* Calculate temporaries and grid functions */
    CCTK_REAL_VEC curvIrL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(ToReal(4),kmul(Psi1iL,Psi3iL),kmadd(ToReal(-4),kmul(Psi1rL,Psi3rL),knmsub(Psi0iL,Psi4iL,kmadd(Psi0rL,Psi4rL,kmadd(ToReal(-3),kmul(Psi2iL,Psi2iL),kmul(kmul(Psi2rL,Psi2rL),ToReal(3)))))));
    
    CCTK_REAL_VEC curvIiL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(ToReal(6),kmul(Psi2iL,Psi2rL),kmadd(ToReal(-4),kmadd(Psi1rL,Psi3iL,kmul(Psi1iL,Psi3rL)),kmadd(Psi0rL,Psi4iL,kmul(Psi0iL,Psi4rL))));
    
    CCTK_REAL_VEC curvJrL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(ToReal(2),kmul(Psi0iL,kmul(Psi3iL,Psi3rL)),kmadd(ToReal(2),kmul(Psi1iL,kmul(Psi1rL,Psi4iL)),knmsub(Psi2iL,kmadd(ToReal(2),kmul(Psi1rL,Psi3iL),kmadd(ToReal(2),kmul(Psi1iL,Psi3rL),kmadd(Psi0rL,Psi4iL,kmul(Psi0iL,Psi4rL)))),kmadd(Psi4rL,kmsub(Psi1iL,Psi1iL,kmul(Psi1rL,Psi1rL)),kmadd(Psi2rL,kmadd(ToReal(-2),kmul(Psi1iL,Psi3iL),kmadd(ToReal(2),kmul(Psi1rL,Psi3rL),knmsub(Psi0iL,Psi4iL,kmadd(Psi0rL,Psi4rL,kmul(kmul(Psi2iL,Psi2iL),ToReal(3)))))),kmsub(Psi0rL,kmsub(Psi3iL,Psi3iL,kmul(Psi3rL,Psi3rL)),kmul(Psi2rL,kmul(Psi2rL,Psi2rL))))))));
    
    CCTK_REAL_VEC curvJiL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(ToReal(2),kmul(Psi2rL,kmadd(Psi1rL,Psi3iL,kmul(Psi1iL,Psi3rL))),kmadd(Psi0rL,kmadd(ToReal(-2),kmul(Psi3iL,Psi3rL),kmul(Psi2rL,Psi4iL)),kmadd(ToReal(-2),kmul(Psi1iL,kmul(Psi1rL,Psi4rL)),kmadd(Psi4iL,kmsub(Psi1iL,Psi1iL,kmul(Psi1rL,Psi1rL)),kmadd(Psi2iL,kmul(Psi2iL,Psi2iL),kmadd(Psi2iL,kmadd(ToReal(-2),kmul(Psi1iL,Psi3iL),kmadd(ToReal(2),kmul(Psi1rL,Psi3rL),knmsub(Psi0iL,Psi4iL,kmadd(Psi0rL,Psi4rL,kmul(kmul(Psi2rL,Psi2rL),ToReal(-3)))))),kmul(Psi0iL,kmadd(Psi2rL,Psi4rL,kmsub(Psi3iL,Psi3iL,kmul(Psi3rL,Psi3rL))))))))));
    
    CCTK_REAL_VEC curvJ1L CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(ToReal(-4),kmul(Psi1iL,Psi3iL),kmadd(ToReal(4),kmul(Psi1rL,Psi3rL),kmadd(Psi0iL,Psi4iL,knmsub(Psi0rL,Psi4rL,kmadd(ToReal(3),kmul(Psi2iL,Psi2iL),kmul(ToReal(-3),kmul(Psi2rL,Psi2rL))))))),ToReal(-16));
    
    CCTK_REAL_VEC curvJ2L CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(ToReal(-2),kmadd(Psi0iL,kmul(Psi3iL,Psi3rL),kmul(Psi1iL,kmul(Psi1rL,Psi4iL))),kmadd(Psi2iL,kmadd(ToReal(2),kmadd(Psi1rL,Psi3iL,kmul(Psi1iL,Psi3rL)),kmadd(Psi0rL,Psi4iL,kmul(Psi0iL,Psi4rL))),kmadd(Psi4rL,kmsub(Psi1rL,Psi1rL,kmul(Psi1iL,Psi1iL)),kmadd(Psi2rL,kmadd(ToReal(2),kmul(Psi1iL,Psi3iL),kmadd(ToReal(-2),kmul(Psi1rL,Psi3rL),kmadd(Psi0iL,Psi4iL,kmsub(ToReal(-3),kmul(Psi2iL,Psi2iL),kmul(Psi0rL,Psi4rL))))),kmadd(Psi2rL,kmul(Psi2rL,Psi2rL),kmul(Psi0rL,kmsub(Psi3rL,Psi3rL,kmul(Psi3iL,Psi3iL)))))))),ToReal(96));
    
    CCTK_REAL_VEC curvJ3L CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(kmadd(kmadd(ToReal(8),kmul(Psi0iL,Psi1iL),kmul(ToReal(-8),kmul(Psi0rL,Psi1rL))),Psi3rL,kmul(ToReal(-4),kmul(Psi0iL,kmul(Psi0rL,Psi4iL)))),Psi4rL,kmadd(ToReal(12),kmul(Psi2iL,kmul(Psi2rL,kmsub(ToReal(4),kmadd(Psi1rL,Psi3iL,kmul(Psi1iL,Psi3rL)),kmadd(Psi0iL,Psi4rL,kmul(Psi0rL,Psi4iL))))),kmadd(Psi1iL,kmul(Psi3iL,kmadd(ToReal(-64),kmul(Psi1rL,Psi3rL),kmadd(ToReal(-8),kmul(Psi0iL,Psi4iL),kmul(ToReal(8),kmul(Psi0rL,Psi4rL))))),kmadd(ToReal(8),kmadd(Psi0rL,kmul(kmadd(Psi1rL,Psi3iL,kmul(Psi1iL,Psi3rL)),Psi4iL),kmul(Psi0iL,kmul(Psi1rL,kmadd(Psi3rL,Psi4iL,kmul(Psi3iL,Psi4rL))))),kmadd(ToReal(6),kmul(kmadd(ToReal(4),kmul(Psi1iL,Psi3iL),kmadd(ToReal(-4),kmul(Psi1rL,Psi3rL),kmsub(Psi0rL,Psi4rL,kmul(Psi0iL,Psi4iL)))),kmul(Psi2rL,Psi2rL)),kmadd(ToReal(-6),kmul(kmul(Psi2iL,Psi2iL),kmadd(ToReal(4),kmul(Psi1iL,Psi3iL),kmadd(ToReal(-4),kmul(Psi1rL,Psi3rL),knmsub(Psi0iL,Psi4iL,kmadd(Psi0rL,Psi4rL,kmul(ToReal(9),kmul(Psi2rL,Psi2rL))))))),kmadd(ToReal(9),kmadd(kmul(Psi2iL,Psi2iL),kmul(Psi2iL,Psi2iL),kmul(kmul(Psi2rL,Psi2rL),kmul(Psi2rL,Psi2rL))),kmadd(kmadd(ToReal(16),kmul(Psi1iL,Psi1iL),kmul(ToReal(-16),kmul(Psi1rL,Psi1rL))),kmul(Psi3iL,Psi3iL),kmadd(kmadd(ToReal(-16),kmul(Psi1iL,Psi1iL),kmul(ToReal(16),kmul(Psi1rL,Psi1rL))),kmul(Psi3rL,Psi3rL),kmadd(kmsub(Psi0iL,Psi0iL,kmul(Psi0rL,Psi0rL)),kmul(Psi4iL,Psi4iL),kmul(kmsub(Psi0rL,Psi0rL,kmul(Psi0iL,Psi0iL)),kmul(Psi4rL,Psi4rL)))))))))))),ToReal(64));
    
    CCTK_REAL_VEC curvJ4L CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(Psi3iL,kmadd(Psi4iL,kmadd(ToReal(-2),kmul(Psi3rL,kmul(Psi0iL,Psi0iL)),kmul(ToReal(12),kmul(Psi1rL,kmul(Psi1iL,Psi1iL)))),kmul(ToReal(4),kmul(Psi4rL,kmul(Psi1iL,kmul(Psi1iL,Psi1iL))))),kmadd(ToReal(-2),kmul(kmadd(ToReal(5),kmul(Psi1iL,Psi3iL),kmadd(ToReal(-5),kmul(Psi1rL,Psi3rL),kmsub(Psi0iL,Psi4iL,kmul(Psi0rL,Psi4rL)))),kmul(Psi2rL,kmul(Psi2rL,Psi2rL))),kmadd(ToReal(-3),kmul(Psi2rL,kmul(kmul(Psi2rL,Psi2rL),kmul(Psi2rL,Psi2rL))),kmadd(kmadd(kmadd(ToReal(12),kmul(Psi0iL,Psi1iL),kmul(ToReal(-12),kmul(Psi0rL,Psi1rL))),Psi3rL,kmsub(ToReal(-2),kmul(Psi0iL,kmul(Psi0rL,Psi4iL)),kmul(Psi4rL,kmul(Psi0iL,Psi0iL)))),kmul(Psi3iL,Psi3iL),kmadd(Psi4rL,kmadd(Psi4iL,kmadd(ToReal(4),kmul(Psi0rL,kmul(Psi1iL,Psi1rL)),kmul(ToReal(-2),kmul(Psi0iL,kmul(Psi1iL,Psi1iL)))),kmadd(Psi3rL,kmadd(ToReal(4),kmul(Psi0iL,kmul(Psi0rL,Psi3iL)),kmul(ToReal(-12),kmul(Psi1rL,kmul(Psi1iL,Psi1iL)))),kmadd(kmul(Psi0rL,Psi0rL),kmsub(Psi3iL,Psi3iL,kmul(Psi3rL,Psi3rL)),kmul(kmul(Psi0iL,Psi0iL),kmul(Psi3rL,Psi3rL))))),kmadd(ToReal(-12),kmadd(Psi1iL,kmul(Psi3rL,kmul(Psi4iL,kmul(Psi1rL,Psi1rL))),kmul(Psi3iL,kmadd(Psi1iL,kmul(Psi4rL,kmul(Psi1rL,Psi1rL)),kmul(kmadd(Psi0rL,Psi1iL,kmul(Psi0iL,Psi1rL)),kmul(Psi3rL,Psi3rL))))),kmadd(ToReal(3),kmadd(kmul(Psi2rL,Psi2rL),kmadd(ToReal(2),kmadd(Psi0iL,kmul(Psi3iL,Psi3rL),kmul(Psi1iL,kmul(Psi1rL,Psi4iL))),kmadd(Psi4rL,kmsub(Psi1iL,Psi1iL,kmul(Psi1rL,Psi1rL)),kmul(Psi0rL,kmsub(Psi3iL,Psi3iL,kmul(Psi3rL,Psi3rL))))),kmul(kmul(Psi2iL,Psi2iL),kmadd(ToReal(-2),kmadd(Psi0iL,kmul(Psi3iL,Psi3rL),kmul(Psi1iL,kmul(Psi1rL,Psi4iL))),kmadd(ToReal(2),kmul(Psi2rL,kmadd(ToReal(5),kmul(Psi1iL,Psi3iL),kmadd(ToReal(-5),kmul(Psi1rL,Psi3rL),kmsub(Psi0iL,Psi4iL,kmul(Psi0rL,Psi4rL))))),kmadd(Psi4rL,kmsub(Psi1rL,Psi1rL,kmul(Psi1iL,Psi1iL)),kmadd(ToReal(10),kmul(Psi2rL,kmul(Psi2rL,Psi2rL)),kmul(Psi0rL,kmsub(Psi3rL,Psi3rL,kmul(Psi3iL,Psi3iL))))))))),kmadd(ToReal(-4),kmadd(Psi3iL,kmul(Psi4iL,kmul(Psi1rL,kmul(Psi1rL,Psi1rL))),kmul(Psi0iL,kmul(Psi1iL,kmul(Psi3rL,kmul(Psi3rL,Psi3rL))))),kmadd(ToReal(4),kmadd(Psi3rL,kmadd(Psi4iL,kmul(Psi1iL,kmul(Psi1iL,Psi1iL)),kmul(Psi4rL,kmul(Psi1rL,kmul(Psi1rL,Psi1rL)))),kmadd(kmadd(Psi0rL,Psi1iL,kmul(Psi0iL,Psi1rL)),kmul(Psi3iL,kmul(Psi3iL,Psi3iL)),kmul(Psi0rL,kmul(Psi1rL,kmul(Psi3rL,kmul(Psi3rL,Psi3rL)))))),kmadd(kmadd(ToReal(-2),kmul(Psi0iL,kmul(Psi1iL,Psi1rL)),kmul(Psi0rL,kmsub(Psi1rL,Psi1rL,kmul(Psi1iL,Psi1iL)))),kmul(Psi4iL,Psi4iL),kmadd(kmadd(ToReal(2),kmul(Psi0iL,kmul(Psi1iL,Psi1rL)),kmul(Psi0rL,kmsub(Psi1iL,Psi1iL,kmul(Psi1rL,Psi1rL)))),kmul(Psi4rL,Psi4rL),kmadd(Psi2rL,kmadd(ToReal(-4),kmul(Psi0iL,kmul(Psi0rL,kmul(Psi4iL,Psi4rL))),kmadd(ToReal(2),kmadd(Psi1rL,kmadd(kmadd(Psi0rL,Psi3iL,kmul(Psi0iL,Psi3rL)),Psi4iL,kmul(kmsub(Psi0iL,Psi3iL,kmul(Psi0rL,Psi3rL)),Psi4rL)),kmul(Psi1iL,kmadd(Psi3iL,kmsub(ToReal(16),kmul(Psi1rL,Psi3rL),kmul(Psi0iL,Psi4iL)),kmadd(Psi0iL,kmul(Psi3rL,Psi4rL),kmul(Psi0rL,kmadd(Psi3rL,Psi4iL,kmul(Psi3iL,Psi4rL))))))),kmadd(ToReal(-15),kmul(kmul(Psi2iL,Psi2iL),kmul(Psi2iL,Psi2iL)),kmadd(kmadd(ToReal(-8),kmul(Psi1iL,Psi1iL),kmul(ToReal(8),kmul(Psi1rL,Psi1rL))),kmsub(Psi3iL,Psi3iL,kmul(Psi3rL,Psi3rL)),kmadd(kmsub(Psi0iL,Psi0iL,kmul(Psi0rL,Psi0rL)),kmul(Psi4iL,Psi4iL),kmul(kmsub(Psi0rL,Psi0rL,kmul(Psi0iL,Psi0iL)),kmul(Psi4rL,Psi4rL))))))),kmul(ToReal(2),kmadd(kmadd(ToReal(5),kmadd(Psi1rL,Psi3iL,kmul(Psi1iL,Psi3rL)),kmadd(Psi0rL,Psi4iL,kmul(Psi0iL,Psi4rL))),kmul(Psi2iL,kmul(Psi2iL,Psi2iL)),kmadd(Psi4iL,kmadd(Psi3iL,kmul(Psi3rL,kmul(Psi0rL,Psi0rL)),kmul(Psi0iL,kmadd(Psi4rL,kmul(Psi1rL,Psi1rL),kmul(Psi0rL,kmul(Psi3rL,Psi3rL))))),kmul(Psi2iL,kmadd(kmadd(ToReal(8),kmul(Psi3iL,Psi3rL),kmul(ToReal(3),kmul(Psi2rL,Psi4iL))),kmul(Psi1rL,Psi1rL),kmadd(Psi3rL,kmadd(Psi3iL,kmadd(ToReal(6),kmul(Psi0rL,Psi2rL),kmul(ToReal(-8),kmul(Psi1iL,Psi1iL))),kmul(Psi1iL,kmadd(Psi0rL,Psi4rL,kmul(ToReal(-15),kmul(Psi2rL,Psi2rL))))),kmadd(Psi4iL,knmsub(Psi0rL,kmul(Psi1iL,Psi3iL),kmadd(Psi4rL,kmsub(Psi0iL,Psi0iL,kmul(Psi0rL,Psi0rL)),kmul(ToReal(-3),kmadd(Psi2rL,kmul(Psi1iL,Psi1iL),kmul(Psi0rL,kmul(Psi2rL,Psi2rL)))))),kmsub(Psi1rL,kmadd(Psi0iL,kmul(Psi3rL,Psi4rL),kmadd(Psi0rL,kmadd(Psi3rL,Psi4iL,kmul(Psi3iL,Psi4rL)),kmadd(Psi3iL,kmsub(ToReal(-15),kmul(Psi2rL,Psi2rL),kmul(Psi0iL,Psi4iL)),kmul(Psi1iL,kmadd(ToReal(6),kmul(Psi2rL,Psi4rL),kmadd(ToReal(-8),kmul(Psi3iL,Psi3iL),kmul(ToReal(8),kmul(Psi3rL,Psi3rL)))))))),kmul(Psi0iL,kmadd(Psi1iL,kmadd(Psi3rL,Psi4iL,kmul(Psi3iL,Psi4rL)),kmadd(ToReal(3),kmul(Psi4rL,kmul(Psi2rL,Psi2rL)),kmadd(Psi2rL,kmadd(ToReal(3),kmul(Psi3iL,Psi3iL),kmul(ToReal(-3),kmul(Psi3rL,Psi3rL))),kmul(Psi0rL,kmsub(Psi4rL,Psi4rL,kmul(Psi4iL,Psi4iL))))))))))))))))))))))))))),ToReal(-640));
    /* Copy local copies back to grid functions */
    vec_store_partial_prepare(i,vecimin,vecimax);
    vec_store_nta_partial(curvIi[index],curvIiL);
    vec_store_nta_partial(curvIr[index],curvIrL);
    vec_store_nta_partial(curvJ1[index],curvJ1L);
    vec_store_nta_partial(curvJ2[index],curvJ2L);
    vec_store_nta_partial(curvJ3[index],curvJ3L);
    vec_store_nta_partial(curvJ4[index],curvJ4L);
    vec_store_nta_partial(curvJi[index],curvJiL);
    vec_store_nta_partial(curvJr[index],curvJrL);
  }
  CCTK_ENDLOOP3STR(WeylScal4_invars_calc_4th);
}
extern "C" void WeylScal4_invars_calc_4th(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering WeylScal4_invars_calc_4th_Body");
  }
  if (cctk_iteration % WeylScal4_invars_calc_4th_calc_every != WeylScal4_invars_calc_4th_calc_offset)
  {
    return;
  }
  
  const char* const groups[] = {
    "WeylScal4::curvIi_group",
    "WeylScal4::curvIr_group",
    "WeylScal4::curvJ1_group",
    "WeylScal4::curvJ2_group",
    "WeylScal4::curvJ3_group",
    "WeylScal4::curvJ4_group",
    "WeylScal4::curvJi_group",
    "WeylScal4::curvJr_group",
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
  AssertGroupStorage(cctkGH, "WeylScal4_invars_calc_4th", 18, groups);
  
  switch (fdOrder)
  {
    case 2:
    {
      break;
    }
    
    case 4:
    {
      break;
    }
    
    case 6:
    {
      break;
    }
    
    case 8:
    {
      break;
    }
    default:
      CCTK_BUILTIN_UNREACHABLE();
  }
  
  LoopOverEverything(cctkGH, WeylScal4_invars_calc_4th_Body);
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving WeylScal4_invars_calc_4th_Body");
  }
}

} // namespace WeylScal4
