 /*@@
   @file      GRHydro_PreLoop.F90
   @date      Mon Feb 25 11:43:36 2002
   @author    
   @desc 
   Sets up various scalars used for efficiency reasons.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

 /*@@
   @routine    GRHydro_Scalar_Setup
   @date       Mon Feb 25 11:25:27 2002
   @author     
   @desc 
   Sets up the logical scalars from the parameters. These are
   solely used for efficiency.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine GRHydro_Scalar_Setup(CCTK_ARGUMENTS)

  USE GRHydro_Scalars

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  MINMOD = .false.
  MC2 = .false.
  SUPERBEE = .false.
  HLLE = .false.
  LLF = .false.

  if (CCTK_EQUALS(tvd_limiter,"minmod")) then
    MINMOD = .true.    
  else if (CCTK_EQUALS(tvd_limiter,"vanleerMC2")) then
    MC2 = .true.
  else if (CCTK_EQUALS(tvd_limiter,"Superbee")) then
    SUPERBEE = .true.
  else
    call CCTK_ERROR("TVD Limiter not recognized!")
    STOP
  end if

  PPM3=.false.
  PPM4=.false.

  if (CCTK_EQUALS(ppm_flatten,"stencil_3")) then
    PPM3=.true.
  else if (CCTK_EQUALS(ppm_flatten,"stencil_4")) then
    PPM4=.true.
  else
    call CCTK_ERROR("PPM Flattening Procedure not recognized!")
    STOP
  end if

  ANALYTICAL = .false.

  if (CCTK_EQUALS(left_eigenvectors,"analytical")) then
    ANALYTICAL = .true.
  else if (CCTK_EQUALS(left_eigenvectors,"numerical")) then
    ANALYTICAL = .false.
  else
    call CCTK_ERROR("Left Eigenvector Mode not recognized!")
    STOP
  end if

  FAST = .false.

  if (CCTK_EQUALS(numerical_viscosity,"fast")) then
    FAST = .true.
  else if (CCTK_EQUALS(numerical_viscosity,"normal")) then
    FAST = .false.
  else
    call CCTK_ERROR("Numerical Viscosity Mode not recognized!")
    STOP
  end if

  if (CCTK_EQUALS(riemann_solver,"HLLE")) then
    HLLE = .true.
  else if (CCTK_EQUALS(riemann_solver,"LLF")) then
    LLF = .true.
  end if

end subroutine GRHydro_Scalar_Setup



 /*@@
   @routine    GRHydro_DivBInit
   @date       Apr 06, 2011
   @author     Bruno Mundim
   @desc 
   Set divB=0 initially.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine GRHydro_DivBInit(CCTK_ARGUMENTS)
  
  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS

  divB=0.0

end subroutine GRHydro_DivBInit


