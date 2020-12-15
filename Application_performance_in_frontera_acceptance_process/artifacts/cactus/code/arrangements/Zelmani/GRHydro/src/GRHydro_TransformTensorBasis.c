#include <assert.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

#define SQR(x) ((x)*(x))

/* helper routine to let Fortran code test for the presence of Multipatch.
 * Required since Fortran cannot access grid scalars via pointers returned by CCTK_VarDataPtrI
 * Returns 0 when no general coordinates are used 
 * Returns 1 when general coordinates are used 
 */ 
CCTK_FCALL CCTK_INT CCTK_FNAME(GRHydro_UseGeneralCoordinates)(CCTK_POINTER ptr_to_cctkGH);
CCTK_INT GRHydro_UseGeneralCoordinates(const cGH * cctkGH);

CCTK_FCALL CCTK_INT CCTK_FNAME(GRHydro_UseGeneralCoordinates)(CCTK_POINTER ptr_to_cctkGH)
{
   return GRHydro_UseGeneralCoordinates(*(cGH **)ptr_to_cctkGH);
}

CCTK_INT GRHydro_UseGeneralCoordinates(const cGH * cctkGH)
{
   static CCTK_INT idxGeneral_coordinates = -1;
   static CCTK_INT coordinates_active_set = 0;
   static CCTK_INT coordinates_active;

   if (!coordinates_active_set) {
      coordinates_active = CCTK_IsImplementationActive("Coordinates");
      coordinates_active_set = 1;
   }

   if (idxGeneral_coordinates == -1)
   {
      idxGeneral_coordinates = CCTK_VarIndex("Coordinates::general_coordinates");
      assert(!coordinates_active || idxGeneral_coordinates >= 0);
   }

   /* are coordinates relly not Cartesian? */
   return coordinates_active && *(CCTK_INT *)CCTK_VarDataPtrI(cctkGH, 0, idxGeneral_coordinates);
}  

void GRHydroTransformADMToLocalBasis(CCTK_ARGUMENTS)
{
   DECLARE_CCTK_ARGUMENTS;
   DECLARE_CCTK_PARAMETERS;

   int i,j,k;
   static CCTK_INT idxJacobian = -1, idxiJacobian = -1;
   static CCTK_INT idxJacobian_state = -1, idxiJacobian_state = -1;
   static CCTK_INT idxGeneral_coordinates = -1;
   CCTK_REAL *J11, *J12, *J13, *J21, *J22, *J23, *J31, *J32, *J33;
   CCTK_REAL *iJ11, *iJ12, *iJ13, *iJ21, *iJ22, *iJ23, *iJ31, *iJ32, *iJ33;

   if (idxGeneral_coordinates == -1)
   {
      idxGeneral_coordinates = CCTK_VarIndex("Coordinates::general_coordinates");
      assert(idxGeneral_coordinates >= 0);
   }

   if (!*(CCTK_INT *)CCTK_VarDataPtrI(cctkGH, 0, idxGeneral_coordinates))
      return; /* corodinates are Cartesian */

   if (idxJacobian == -1)
   {
      idxJacobian = CCTK_FirstVarIndex("Coordinates::jacobian");
      assert(idxJacobian >= 0);
   }

   if (idxiJacobian == -1)
   {
      idxiJacobian = CCTK_FirstVarIndex("Coordinates::inverse_jacobian");
      assert(idxiJacobian >= 0);
   }

   if (idxJacobian_state == -1)
   {
      idxJacobian_state = CCTK_VarIndex("Coordinates::jacobian_state");
      assert(idxJacobian_state >= 0);
   }

   if (idxiJacobian_state == -1)
   {
      idxiJacobian_state = CCTK_VarIndex("Coordinates::inverse_jacobian_state");
      assert(idxiJacobian_state >= 0);
   }

   if (!*(CCTK_INT *)CCTK_VarDataPtrI(cctkGH, 0, idxJacobian_state))
      CCTK_WARN(0,"No storage for Jacobians allocated! Tell your coordinates implemetation to store Jacobians!");

   if (!*(CCTK_INT *)CCTK_VarDataPtrI(cctkGH, 0, idxiJacobian_state))
      CCTK_WARN(0,"No storage for inverse Jacobians allocated! Tell your coordinates implemetation to store inverse Jacobians!");

   /* Cactus guarantess that variables in a group are consecutive, and
    * ReflectionSymmetry etc. require the order we use */
   J11 = CCTK_VarDataPtrI(cctkGH, 0, idxJacobian+0);
   J12 = CCTK_VarDataPtrI(cctkGH, 0, idxJacobian+1);
   J13 = CCTK_VarDataPtrI(cctkGH, 0, idxJacobian+2);
   J21 = CCTK_VarDataPtrI(cctkGH, 0, idxJacobian+3);
   J22 = CCTK_VarDataPtrI(cctkGH, 0, idxJacobian+4);
   J23 = CCTK_VarDataPtrI(cctkGH, 0, idxJacobian+5);
   J31 = CCTK_VarDataPtrI(cctkGH, 0, idxJacobian+6);
   J32 = CCTK_VarDataPtrI(cctkGH, 0, idxJacobian+7);
   J33 = CCTK_VarDataPtrI(cctkGH, 0, idxJacobian+8);
   assert(J11 && J12 && J13 && J21 && J22 && J23 && J31 && J32 && J33);
   iJ11 = CCTK_VarDataPtrI(cctkGH, 0, idxiJacobian+0);
   iJ12 = CCTK_VarDataPtrI(cctkGH, 0, idxiJacobian+1);
   iJ13 = CCTK_VarDataPtrI(cctkGH, 0, idxiJacobian+2);
   iJ21 = CCTK_VarDataPtrI(cctkGH, 0, idxiJacobian+3);
   iJ22 = CCTK_VarDataPtrI(cctkGH, 0, idxiJacobian+4);
   iJ23 = CCTK_VarDataPtrI(cctkGH, 0, idxiJacobian+5);
   iJ31 = CCTK_VarDataPtrI(cctkGH, 0, idxiJacobian+6);
   iJ32 = CCTK_VarDataPtrI(cctkGH, 0, idxiJacobian+7);
   iJ33 = CCTK_VarDataPtrI(cctkGH, 0, idxiJacobian+8);
   assert(iJ11 && iJ12 && iJ13 && iJ21 && iJ22 && iJ23 && iJ31 && iJ32 && iJ33);

   #pragma omp parallel for private (i,j,k)
   for (k=0 ; k<cctk_lsh[2] ; k++)
   {
      for (j=0 ;  j<cctk_lsh[1] ; j++)
      {
         for (i=0 ;  i<cctk_lsh[0] ; i++)
         {
            CCTK_INT idx = CCTK_GFINDEX3D(cctkGH, i,j,k);

            gaa[idx] = SQR(iJ11[idx]) * gxx[idx] +
                 2.0 * iJ11[idx] * iJ21[idx] * gxy[idx] +
                 2.0 * iJ11[idx] * iJ31[idx] * gxz[idx] + SQR(iJ21[idx]) * gyy[idx] +
                 2.0 * iJ21[idx] * iJ31[idx] * gyz[idx] + SQR(iJ31[idx]) * gzz[idx];

            gab[idx] = iJ11[idx] * iJ12[idx] * gxx[idx] +
                 (iJ11[idx] * iJ22[idx] + iJ21[idx] * iJ12[idx]) * gxy[idx] +
                 (iJ11[idx] * iJ32[idx] + iJ31[idx] * iJ12[idx]) * gxz[idx] +
                 iJ21[idx] * iJ22[idx] * gyy[idx] + (iJ21[idx] * iJ32[idx] +
                 iJ31[idx] * iJ22[idx]) * gyz[idx] + iJ31[idx] * iJ32[idx] * gzz[idx];

            gac[idx] = iJ11[idx] * iJ13[idx] * gxx[idx] + (iJ11[idx] * iJ23[idx] +
                 iJ21[idx] * iJ13[idx]) * gxy[idx] +
                 (iJ11[idx] * iJ33[idx] + iJ31[idx] * iJ13[idx]) * gxz[idx] +
                 iJ21[idx] * iJ23[idx] * gyy[idx] +
                 (iJ21[idx] * iJ33[idx] + iJ31[idx] * iJ23[idx]) * gyz[idx] +
                 iJ31[idx] * iJ33[idx] * gzz[idx];

            gbb[idx] = SQR(iJ12[idx]) * gxx[idx] +
                 2.0 * iJ12[idx] * iJ22[idx] * gxy[idx] +
                 2.0 * iJ12[idx] * iJ32[idx] * gxz[idx] + SQR(iJ22[idx]) * gyy[idx] +
                 2.0 * iJ22[idx] * iJ32[idx] * gyz[idx] + SQR(iJ32[idx]) * gzz[idx];

            gbc[idx] = iJ12[idx] * iJ13[idx] * gxx[idx] +
                 (iJ12[idx] * iJ23[idx] + iJ22[idx] * iJ13[idx]) * gxy[idx] +
                 (iJ12[idx] * iJ33[idx] + iJ32[idx] * iJ13[idx]) * gxz[idx] +
                 iJ22[idx] * iJ23[idx] * gyy[idx] +
                 (iJ22[idx]* iJ33[idx] + iJ32[idx] * iJ23[idx]) * gyz[idx] +
                 iJ32[idx] * iJ33[idx] * gzz[idx];

            gcc[idx] = SQR(iJ13[idx]) * gxx[idx] +
                 2.0 * iJ13[idx] * iJ23[idx] * gxy[idx] +
                 2.0 * iJ13[idx] * iJ33[idx] * gxz[idx] +
                 SQR(iJ23[idx]) * gyy[idx] +
                 2.0 * iJ23[idx] * iJ33[idx] * gyz[idx] +
                 SQR(iJ33[idx]) * gzz[idx];

            /* Transform extrinsic curvature from global to local basis.
             * Since extrinsic has covariant indices, use inverse Jacobian */
            kaa[idx] = SQR(iJ11[idx]) * kxx[idx] +
                 2.0 * iJ11[idx] * iJ21[idx] * kxy[idx] +
                 2.0 * iJ11[idx] * iJ31[idx] * kxz[idx] +
                 SQR(iJ21[idx]) * kyy[idx] +
                 2.0 * iJ21[idx] * iJ31[idx] * kyz[idx] +
                 SQR(iJ31[idx]) * kzz[idx];

            kab[idx] = iJ11[idx] * iJ12[idx] * kxx[idx] +
                 (iJ11[idx] * iJ22[idx] + iJ21[idx] * iJ12[idx]) * kxy[idx] +
                 (iJ11[idx] * iJ32[idx] + iJ31[idx] * iJ12[idx]) * kxz[idx] +
                 iJ21[idx] * iJ22[idx] * kyy[idx] +
                 (iJ21[idx] * iJ32[idx] + iJ31[idx] * iJ22[idx]) * kyz[idx] +
                 iJ31[idx] * iJ32[idx] * kzz[idx];


            kac[idx] = iJ11[idx] * iJ13[idx] * kxx[idx] +
                 (iJ11[idx] * iJ23[idx] + iJ21[idx] * iJ13[idx]) * kxy[idx] +
                 (iJ11[idx] * iJ33[idx] + iJ31[idx] * iJ13[idx]) * kxz[idx] +
                 iJ21[idx] * iJ23[idx] * kyy[idx] +
                 (iJ21[idx] * iJ33[idx] + iJ31[idx] * iJ23[idx]) * kyz[idx] +
                 iJ31[idx] * iJ33[idx] * kzz[idx];

            kbb[idx] = SQR(iJ12[idx]) * kxx[idx] +
                 2.0 * iJ12[idx] * iJ22[idx] * kxy[idx] +
                 2.0 * iJ12[idx] * iJ32[idx] * kxz[idx] +
                 SQR(iJ22[idx]) * kyy[idx] +
                 2.0 * iJ22[idx] * iJ32[idx] * kyz[idx] +
                 SQR(iJ32[idx]) * kzz[idx];

            kbc[idx] = iJ12[idx] * iJ13[idx] * kxx[idx] +
                 (iJ12[idx] * iJ23[idx] + iJ22[idx] * iJ13[idx]) * kxy[idx] +
                 (iJ12[idx] * iJ33[idx] + iJ32[idx] * iJ13[idx]) * kxz[idx] +
                 iJ22[idx] * iJ23[idx] * kyy[idx] +
                 (iJ22[idx] * iJ33[idx] + iJ32[idx] * iJ23[idx]) * kyz[idx] +
                 iJ32[idx] * iJ33[idx] * kzz[idx];

            kcc[idx] = SQR(iJ13[idx]) * kxx[idx] +
                 2.0 * iJ13[idx] * iJ23[idx] * kxy[idx] +
                 2.0 * iJ13[idx] * iJ33[idx] * kxz[idx] +
                 SQR(iJ23[idx]) * kyy[idx] +
                 2.0 * iJ23[idx] * iJ33[idx] * kyz[idx] +
                 SQR(iJ33[idx]) * kzz[idx];
            
            /* Transform shift from global to local basis.
             * Since shift has contravariant index, use Jacobian */
            betaa[idx] = betax[idx]*J11[idx] + betay[idx]*J12[idx] + betaz[idx]*J13[idx];
            betab[idx] = betax[idx]*J21[idx] + betay[idx]*J22[idx] + betaz[idx]*J23[idx];
            betac[idx] = betax[idx]*J31[idx] + betay[idx]*J32[idx] + betaz[idx]*J33[idx];
         }
      }
   }

}



void GRHydroTransformPrimToLocalBasis(CCTK_ARGUMENTS)
{
   
   DECLARE_CCTK_ARGUMENTS;
   DECLARE_CCTK_PARAMETERS;

   CCTK_INT i,j,k;
   static CCTK_INT idxJacobian = -1, idxiJacobian = -1;
   static CCTK_INT idxGeneral_coordinates = -1;
   static CCTK_INT idxJacobian_state = -1, idxiJacobian_state = -1;
   CCTK_REAL *J11, *J12, *J13, *J21, *J22, *J23, *J31, *J32, *J33;
   CCTK_REAL *iJ11, *iJ12, *iJ13, *iJ21, *iJ22, *iJ23, *iJ31, *iJ32, *iJ33;

   if (!CCTK_IsImplementationActive("Coordinates"))
      return; /* Multipatch infrastructure not active */

   if (idxGeneral_coordinates == -1)
   {
      idxGeneral_coordinates = CCTK_VarIndex("Coordinates::general_coordinates");
      assert(idxGeneral_coordinates >= 0);
   }

   if (!*(CCTK_INT *)CCTK_VarDataPtrI(cctkGH, 0, idxGeneral_coordinates))
      return; /* corodinates are Cartesian */

   if (idxJacobian == -1)
   {
      idxJacobian = CCTK_FirstVarIndex("Coordinates::jacobian");
      assert(idxJacobian >= 0);
   }

   if (idxiJacobian == -1)
   {
      idxiJacobian = CCTK_FirstVarIndex("Coordinates::inverse_jacobian");
      assert(idxiJacobian >= 0);
   }

   if (idxJacobian_state == -1)
   {
      idxJacobian_state = CCTK_VarIndex("Coordinates::jacobian_state");
      assert(idxJacobian_state >= 0);
   }

   if (idxiJacobian_state == -1)
   {
      idxiJacobian_state = CCTK_VarIndex("Coordinates::inverse_jacobian_state");
      assert(idxiJacobian_state >= 0);
   }

   if (!*(CCTK_INT *)CCTK_VarDataPtrI(cctkGH, 0, idxJacobian_state))
      CCTK_WARN(0,"No storage for Jacobians allocated! Tell your coordinates implemetation to store Jacobians!");

   if (!*(CCTK_INT *)CCTK_VarDataPtrI(cctkGH, 0, idxiJacobian_state))
      CCTK_WARN(0,"No storage for inverse Jacobians allocated! Tell your coordinates implemetation to store inverse Jacobians!");

   /* CCTK_INFO("Transforming primitives to local basis."); */

   /* Cactus guarantess that variables in a group are consecutive, and
    * ReflectionSymmetry etc. require the order we use */
   J11 = CCTK_VarDataPtrI(cctkGH, 0, idxJacobian+0);
   J12 = CCTK_VarDataPtrI(cctkGH, 0, idxJacobian+1);
   J13 = CCTK_VarDataPtrI(cctkGH, 0, idxJacobian+2);
   J21 = CCTK_VarDataPtrI(cctkGH, 0, idxJacobian+3);
   J22 = CCTK_VarDataPtrI(cctkGH, 0, idxJacobian+4);
   J23 = CCTK_VarDataPtrI(cctkGH, 0, idxJacobian+5);
   J31 = CCTK_VarDataPtrI(cctkGH, 0, idxJacobian+6);
   J32 = CCTK_VarDataPtrI(cctkGH, 0, idxJacobian+7);
   J33 = CCTK_VarDataPtrI(cctkGH, 0, idxJacobian+8);
   assert(J11 && J12 && J13 && J21 && J22 && J23 && J31 && J32 && J33);
   iJ11 = CCTK_VarDataPtrI(cctkGH, 0, idxiJacobian+0);
   iJ12 = CCTK_VarDataPtrI(cctkGH, 0, idxiJacobian+1);
   iJ13 = CCTK_VarDataPtrI(cctkGH, 0, idxiJacobian+2);
   iJ21 = CCTK_VarDataPtrI(cctkGH, 0, idxiJacobian+3);
   iJ22 = CCTK_VarDataPtrI(cctkGH, 0, idxiJacobian+4);
   iJ23 = CCTK_VarDataPtrI(cctkGH, 0, idxiJacobian+5);
   iJ31 = CCTK_VarDataPtrI(cctkGH, 0, idxiJacobian+6);
   iJ32 = CCTK_VarDataPtrI(cctkGH, 0, idxiJacobian+7);
   iJ33 = CCTK_VarDataPtrI(cctkGH, 0, idxiJacobian+8);
   assert(iJ11 && iJ12 && iJ13 && iJ21 && iJ22 && iJ23 && iJ31 && iJ32 && iJ33);

   #pragma omp parallel for private (i,j,k)
   for (k=0 ; k<cctk_lsh[2] ; k++)
   {
      for (j=0 ;  j<cctk_lsh[1] ; j++)
      {
         for (i=0 ;  i<cctk_lsh[0] ; i++)
         {
            CCTK_INT idx = CCTK_GFINDEX3D(cctkGH, i,j,k);
            CCTK_INT idx1 = CCTK_VECTGFINDEX3D(cctkGH, i,j,k,0);
            CCTK_INT idx2 = CCTK_VECTGFINDEX3D(cctkGH, i,j,k,1);
            CCTK_INT idx3 = CCTK_VECTGFINDEX3D(cctkGH, i,j,k,2);
   
            /* Transform primitive velocity from global to local basis.
             * Since velocity has contravariant index, use Jacobian */
            lvel[idx1] = vel[idx1]*J11[idx] + vel[idx2]*J12[idx] + vel[idx3]*J13[idx];
            lvel[idx2] = vel[idx1]*J21[idx] + vel[idx2]*J22[idx] + vel[idx3]*J23[idx];
            lvel[idx3] = vel[idx1]*J31[idx] + vel[idx2]*J32[idx] + vel[idx3]*J33[idx];
            
            if(*evolve_MHD)
            {
              /* Transform primitive B-field from global to local basis.
               * Since B-field has contravariant index, use Jacobian */
              lBvec[idx1] = Bvec[idx1]*J11[idx] + Bvec[idx2]*J12[idx] + Bvec[idx3]*J13[idx];
              lBvec[idx2] = Bvec[idx1]*J21[idx] + Bvec[idx2]*J22[idx] + Bvec[idx3]*J23[idx];
              lBvec[idx3] = Bvec[idx1]*J31[idx] + Bvec[idx2]*J32[idx] + Bvec[idx3]*J33[idx];
            }
         }
      }
   }

}



void GRHydroTransformPrimToGlobalBasis(CCTK_ARGUMENTS)
{
   
   DECLARE_CCTK_ARGUMENTS;
   DECLARE_CCTK_PARAMETERS;

   CCTK_INT i,j,k;
   
   static CCTK_INT idxJacobian = -1, idxiJacobian = -1;
   static CCTK_INT idxJacobian_state = -1, idxiJacobian_state = -1;
   static CCTK_INT idxGeneral_coordinates = -1;
   CCTK_REAL *J11, *J12, *J13, *J21, *J22, *J23, *J31, *J32, *J33;
   CCTK_REAL *iJ11, *iJ12, *iJ13, *iJ21, *iJ22, *iJ23, *iJ31, *iJ32, *iJ33;

   if (idxGeneral_coordinates == -1)
   {
      idxGeneral_coordinates = CCTK_VarIndex("Coordinates::general_coordinates");
      assert(idxGeneral_coordinates >= 0);
   }

   if (!*(CCTK_INT *)CCTK_VarDataPtrI(cctkGH, 0, idxGeneral_coordinates))
      return; /* corodinates are Cartesian */

   if (idxJacobian == -1)
   {
      idxJacobian = CCTK_FirstVarIndex("Coordinates::jacobian");
      assert(idxJacobian >= 0);
   }

   if (idxiJacobian == -1)
   {
      idxiJacobian = CCTK_FirstVarIndex("Coordinates::inverse_jacobian");
      assert(idxiJacobian >= 0);
   }

   if (idxJacobian_state == -1)
   {
      idxJacobian_state = CCTK_VarIndex("Coordinates::jacobian_state");
      assert(idxJacobian_state >= 0);
   }

   if (idxiJacobian_state == -1)
   {
      idxiJacobian_state = CCTK_VarIndex("Coordinates::inverse_jacobian_state");
      assert(idxiJacobian_state >= 0);
   }

   if (!*(CCTK_INT *)CCTK_VarDataPtrI(cctkGH, 0, idxJacobian_state))
      CCTK_WARN(0,"No storage for Jacobians allocated! Tell your coordinates implemetation to store Jacobians!");

   if (!*(CCTK_INT *)CCTK_VarDataPtrI(cctkGH, 0, idxiJacobian_state))
      CCTK_WARN(0,"No storage for inverse Jacobians allocated! Tell your coordinates implemetation to store inverse Jacobians!");

   /* CCTK_INFO("Transforming primitives to global basis."); */

   /* Cactus guarantess that variables in a group are consecutive, and
    * ReflectionSymmetry etc. require the order we use */
   J11 = CCTK_VarDataPtrI(cctkGH, 0, idxJacobian+0);
   J12 = CCTK_VarDataPtrI(cctkGH, 0, idxJacobian+1);
   J13 = CCTK_VarDataPtrI(cctkGH, 0, idxJacobian+2);
   J21 = CCTK_VarDataPtrI(cctkGH, 0, idxJacobian+3);
   J22 = CCTK_VarDataPtrI(cctkGH, 0, idxJacobian+4);
   J23 = CCTK_VarDataPtrI(cctkGH, 0, idxJacobian+5);
   J31 = CCTK_VarDataPtrI(cctkGH, 0, idxJacobian+6);
   J32 = CCTK_VarDataPtrI(cctkGH, 0, idxJacobian+7);
   J33 = CCTK_VarDataPtrI(cctkGH, 0, idxJacobian+8);
   assert(J11 && J12 && J13 && J21 && J22 && J23 && J31 && J32 && J33);
   iJ11 = CCTK_VarDataPtrI(cctkGH, 0, idxiJacobian+0);
   iJ12 = CCTK_VarDataPtrI(cctkGH, 0, idxiJacobian+1);
   iJ13 = CCTK_VarDataPtrI(cctkGH, 0, idxiJacobian+2);
   iJ21 = CCTK_VarDataPtrI(cctkGH, 0, idxiJacobian+3);
   iJ22 = CCTK_VarDataPtrI(cctkGH, 0, idxiJacobian+4);
   iJ23 = CCTK_VarDataPtrI(cctkGH, 0, idxiJacobian+5);
   iJ31 = CCTK_VarDataPtrI(cctkGH, 0, idxiJacobian+6);
   iJ32 = CCTK_VarDataPtrI(cctkGH, 0, idxiJacobian+7);
   iJ33 = CCTK_VarDataPtrI(cctkGH, 0, idxiJacobian+8);
   assert(iJ11 && iJ12 && iJ13 && iJ21 && iJ22 && iJ23 && iJ31 && iJ32 && iJ33);

   #pragma omp parallel for private (i,j,k)
   for (k=0 ; k<cctk_lsh[2] ; k++)
   {
      for (j=0 ;  j<cctk_lsh[1] ; j++)
      {
         for (i=0 ;  i<cctk_lsh[0] ; i++)
         {
            CCTK_INT idx = CCTK_GFINDEX3D(cctkGH, i,j,k);
            CCTK_INT idx1 = CCTK_VECTGFINDEX3D(cctkGH, i,j,k,0);
            CCTK_INT idx2 = CCTK_VECTGFINDEX3D(cctkGH, i,j,k,1);
            CCTK_INT idx3 = CCTK_VECTGFINDEX3D(cctkGH, i,j,k,2);
   
            /* Transform primitive velocity from local to global basis.
             * Since velocity has contravariant index, use inverse Jacobian */
            vel[idx1] = lvel[idx1]*iJ11[idx] + lvel[idx2]*iJ12[idx] + lvel[idx3]*iJ13[idx];
            vel[idx2] = lvel[idx1]*iJ21[idx] + lvel[idx2]*iJ22[idx] + lvel[idx3]*iJ23[idx];
            vel[idx3] = lvel[idx1]*iJ31[idx] + lvel[idx2]*iJ32[idx] + lvel[idx3]*iJ33[idx];

            if(*evolve_MHD)
            {
              /* Transform primitive B-field from local to global basis.
               * Since B-field has contravariant index, use inverse Jacobian */
              Bvec[idx1] = lBvec[idx1]*iJ11[idx] + lBvec[idx2]*iJ12[idx] + lBvec[idx3]*iJ13[idx];
              Bvec[idx2] = lBvec[idx1]*iJ21[idx] + lBvec[idx2]*iJ22[idx] + lBvec[idx3]*iJ23[idx];
              Bvec[idx3] = lBvec[idx1]*iJ31[idx] + lBvec[idx2]*iJ32[idx] + lBvec[idx3]*iJ33[idx];
            }
            
         }
      }
   }

}





