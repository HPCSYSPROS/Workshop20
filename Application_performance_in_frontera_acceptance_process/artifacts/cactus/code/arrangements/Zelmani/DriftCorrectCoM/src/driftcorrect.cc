#include <cstdio>
#include <string>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include <assert.h>

using namespace std;


extern "C" void dcm_correct_drift(CCTK_ARGUMENTS)
{
   DECLARE_CCTK_ARGUMENTS;
   DECLARE_CCTK_PARAMETERS;

   const CCTK_REAL local_eps = 1e-12;

   const CCTK_REAL x0[3] = { CCTK_ORIGIN_SPACE(0), 
                             CCTK_ORIGIN_SPACE(1),
                             CCTK_ORIGIN_SPACE(2)  };
   
   const CCTK_REAL dx[3] = { CCTK_DELTA_SPACE(0), 
                             CCTK_DELTA_SPACE(1),
                             CCTK_DELTA_SPACE(2)  };

   const CCTK_REAL dt = CCTK_DELTA_TIME;
   
   
   if (cctk_time < first_driftcorrect_time) 
      return;

   // Check for valid horizon data
   if (*dc4_calc_error > 1.0e-10) {
      CCTK_WARN (1, "No valid drift correction data found! Not doing anything!");
      return;
   } else {
     if(verbose_level > 1) {
       CCTK_INFO("Applying corrections!");
     }
   }

   
   int shiftRhsIdx = -1;
   
   if (CCTK_IsImplementationActive("CTGGauge")) {
      shiftRhsIdx = CCTK_FirstVarIndex("ADMBase::dtshift");
   } else if (CCTK_IsImplementationActive("ML_BSSN")) {
      shiftRhsIdx = CCTK_FirstVarIndex("ML_BSSN::ML_shiftrhs");
   } else {
      CCTK_WARN(0, "No valid gauge evolution thorn found!");
   }
   
   if (shiftRhsIdx < 0)
   {
      CCTK_WARN(0, "Error getting shift rhs pointers!");
      return;
   }
   
   CCTK_REAL* const shiftxRHS = static_cast<CCTK_REAL*>(CCTK_VarDataPtrI(cctkGH, 0, shiftRhsIdx+0));
   CCTK_REAL* const shiftyRHS = static_cast<CCTK_REAL*>(CCTK_VarDataPtrI(cctkGH, 0, shiftRhsIdx+1));
   CCTK_REAL* const shiftzRHS = static_cast<CCTK_REAL*>(CCTK_VarDataPtrI(cctkGH, 0, shiftRhsIdx+2));
   assert(shiftxRHS && shiftyRHS && shiftzRHS);
   
#pragma omp parallel for
   for (int k=0; k < cctk_lsh[2]; ++k) {
      for (int j=0; j < cctk_lsh[1]; ++j) {
         for (int i=0; i < cctk_lsh[0]; ++i) {
         
            const CCTK_INT ijk = CCTK_GFINDEX3D(cctkGH, i,j,k);
         
            const CCTK_REAL xpos[3] = { x0[0] + dx[0] * (cctk_lbnd[0] + i),
                                        x0[1] + dx[1] * (cctk_lbnd[1] + j),
                                        x0[2] + dx[2] * (cctk_lbnd[2] + k)  };
            
            const CCTK_REAL radius = sqrt(xpos[0]*xpos[0] + xpos[1]*xpos[1] + xpos[2]*xpos[2]);
            
            const CCTK_REAL lpos[3] = { xpos[0] - position_x, 
                                        xpos[1] - position_y, 
                                        xpos[2] - position_z };
                                        
            const CCTK_REAL lradius = sqrt(lpos[0]*lpos[0] + lpos[1]*lpos[1] + lpos[2]*lpos[2]);
                                  
            const CCTK_REAL lnormal[3] = { lpos[0] / (lradius + local_eps),
                                           lpos[1] / (lradius + local_eps),
                                           lpos[2] / (lradius + local_eps) };
            

            // Position correction
            if (do_position_correction != 0)
            {
               shiftxRHS[ijk] = shiftxRHS[ijk]
                    - (*dc4_delta_posx_dot2)*
                    exp(-lradius*lradius*position_correction_falloff);
               shiftyRHS[ijk] = shiftyRHS[ijk]
                    - (*dc4_delta_posy_dot2)*
                    exp(-lradius*lradius*position_correction_falloff);
               shiftzRHS[ijk] = shiftzRHS[ijk]
                    - (*dc4_delta_posz_dot2)*
                    exp(-lradius*lradius*position_correction_falloff);
            }
         }
      }
   }
}





