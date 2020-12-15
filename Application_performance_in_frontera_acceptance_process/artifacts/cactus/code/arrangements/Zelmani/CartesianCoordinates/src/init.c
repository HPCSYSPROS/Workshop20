#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void CartesianCoordinates_Init(CCTK_ARGUMENTS)
{
   DECLARE_CCTK_ARGUMENTS;
   DECLARE_CCTK_PARAMETERS;
   
   *general_coordinates = 0;
   
   
   if (store_jacobian)
   {
      *jacobian_state = 1;
     
      for (int k=0; k < cctk_lsh[2]; ++k) {
	for (int j=0; j < cctk_lsh[1]; ++j) {
	  for (int i=0; i < cctk_lsh[0]; ++i)
	  {
	    const int ijk = CCTK_GFINDEX3D(cctkGH, i,j,k);
	    J11[ijk] = 1.0;
	    J12[ijk] = 0.0;
	    J13[ijk] = 0.0;

	    J21[ijk] = 0.0;
	    J22[ijk] = 1.0;
	    J23[ijk] = 0.0;

	    J31[ijk] = 0.0;
	    J32[ijk] = 0.0;
	    J33[ijk] = 1.0;
	    
	  }
        }
      }
   }
   else
   {
      *jacobian_state = 0;
   }
   
   if (store_jacobian_derivative)
   {
      *jacobian_derivative_state = 1;
     
      for (int k=0; k < cctk_lsh[2]; ++k) {
	for (int j=0; j < cctk_lsh[1]; ++j) {
	  for (int i=0; i < cctk_lsh[0]; ++i)
	  {
	    const int ijk = CCTK_GFINDEX3D(cctkGH, i,j,k);
	    
	    dJ111[ijk] = 0.0; dJ112[ijk] = 0.0; dJ113[ijk] = 0.0; dJ122[ijk] = 0.0; dJ123[ijk] = 0.0; dJ133[ijk] = 0.0;
            dJ211[ijk] = 0.0; dJ212[ijk] = 0.0; dJ213[ijk] = 0.0; dJ222[ijk] = 0.0; dJ223[ijk] = 0.0; dJ233[ijk] = 0.0;
            dJ311[ijk] = 0.0; dJ312[ijk] = 0.0; dJ313[ijk] = 0.0; dJ322[ijk] = 0.0; dJ323[ijk] = 0.0; dJ333[ijk] = 0.0;
	  }
        }
      }
   }
   
   if (store_inverse_jacobian)
   {
      *inverse_jacobian_state = 1;
     
      for (int k=0; k < cctk_lsh[2]; ++k) {
	for (int j=0; j < cctk_lsh[1]; ++j) {
	  for (int i=0; i < cctk_lsh[0]; ++i)
	  {
	    const int ijk = CCTK_GFINDEX3D(cctkGH, i,j,k);
	    iJ11[ijk] = 1.0;
	    iJ12[ijk] = 0.0;
	    iJ13[ijk] = 0.0;

	    iJ21[ijk] = 0.0;
	    iJ22[ijk] = 1.0;
	    iJ23[ijk] = 0.0;

	    iJ31[ijk] = 0.0;
	    iJ32[ijk] = 0.0;
	    iJ33[ijk] = 1.0;
	    
	  }
        }
      }
   }
   else
   {
      *inverse_jacobian_state = 0;
   }
   
   if (store_volume_form)
   {
      *volume_form_state = 1;
     
      for (int k=0; k < cctk_lsh[2]; ++k) {
	for (int j=0; j < cctk_lsh[1]; ++j) {
	  for (int i=0; i < cctk_lsh[0]; ++i)
	  {
	    const int ijk = CCTK_GFINDEX3D(cctkGH, i,j,k);
	    // this is the product of the _coarse_ grid spacing: during redcution, Carpet will automatically rescale
	    // the fine grids!
	    volume_form[ijk] = cctk_delta_space[0]*cctk_delta_space[1]*cctk_delta_space[2];
	  }
        }
      }
   }
   else
   {
      *volume_form_state = 0;
   }
}


