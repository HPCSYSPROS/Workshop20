 /*@@
   @file      GRHydro_Minima.F90
   @date      Mon Feb 25 11:43:36 2002
   @author    
   @desc 
   Sets up the scalars used for the atmosphere, before initial data.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

#include "GRHydro_Macros.h"

 /*@@
   @routine    GRHydro_Minima_Setup
   @date       Mon Feb 25 11:25:27 2002
   @author     Ian Hawke
   @desc 
   Before initial data, set up the scalar GRHydro_rho_min used for the atmosphere.
   This is computed only from parameters.
   @enddesc 
   @calls     
   @calledby   
   @history 
   Modified 30 Aug 2006 by Luca Baiotti
   @endhistory 

@@*/

subroutine GRHydro_Rho_Minima_Setup(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  if (initial_rho_abs_min > 0.0) then
    GRHydro_rho_min = initial_rho_abs_min
  else if (initial_rho_rel_min > 0.0) then 
    GRHydro_rho_min = GRHydro_rho_central * initial_rho_rel_min
  else if (rho_abs_min > 0.d0) then
    GRHydro_rho_min = rho_abs_min
  else
    GRHydro_rho_min = GRHydro_rho_central * rho_rel_min
  end if

  if (initial_atmosphere_factor > 0.0) GRHydro_rho_min = GRHydro_rho_min * initial_atmosphere_factor

  GRHydro_tau_min = tau_rel_min

  return

end subroutine GRHydro_Rho_Minima_Setup


 /*@@
   @routine    GRHydro_Check_Rho_Minimum
   @date       Mon Jul 7 16:35:45 2008
   @author     Luca Baiotti 
   @desc 
   Check whether at some point rho < GRHydro_rho_min and print a warning in case.
   @enddesc 
   @calls     
   @calledby   
   @history 
   @endhistory 
@@*/

subroutine GRHydro_Check_Rho_Minimum(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT i,j,k
  CCTK_REAL dummy1
  character(len=100) warnline

  do i=1,cctk_lsh(1)
    do j=1,cctk_lsh(2)
      do k=1,cctk_lsh(3)

        SET_ATMO_MIN(dummy1, GRHydro_rho_min, r(i,j,k))
        if (rho(i,j,k) < dummy1) then
          call CCTK_WARN(2,"rho<GRHydro_rho_min!!!")
          write(warnline,'(a28,i2)') 'on carpet reflevel: ',GRHydro_reflevel
          call CCTK_WARN(2,warnline)
          write(warnline,'(a25,g15.6)') 'GRHydro_rho_min: ', dummy1
          call CCTK_WARN(2,warnline)
          write(warnline,'(a25,g15.6)') 'rho: ',rho(i,j,k)
          call CCTK_WARN(2,warnline)
          write(warnline,'(a25,4g15.6)') 'coordinates: x,y,z,r:',x(i,j,k),y(i,j,k),z(i,j,k),r(i,j,k)
          call CCTK_WARN(2,warnline)
        end if

      end do
    end do
  end do

  return

end subroutine GRHydro_Check_Rho_Minimum


 /*@@
   @routine    GRHydro__Change_Rho_Minimum_At_Recovery
   @date       Thu Aug 14 17:11:32 2008
   @author     Luca Baiotti 
   @desc 
   Change, via a parameter, the value of GRHydro_rho_min at recovery.
   @enddesc 
   @calls     
   @calledby   
   @history 
   @endhistory 
@@*/

subroutine GRHydro_Change_Rho_Minimum_At_Recovery(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  GRHydro_rho_min = rho_abs_min_after_recovery

  return

end subroutine GRHydro_Change_Rho_Minimum_At_Recovery
