 /*@@
   @file      GRHydro_RiemannSolve.F90
   @date      Sat Jan 26 02:20:25 2002
   @author    
   @desc 
   A wrapper routine to call the correct Riemann solver
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "GRHydro_Macros.h"

 /*@@
   @routine    RiemannSolve
   @date       Sat Jan 26 02:20:48 2002
   @author     Pedro Montero, Ian Hawke
   @desc 
   A wrapper routine to switch between the different Riemann solvers.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine RiemannSolve(CCTK_ARGUMENTS)

  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  if (CCTK_EQUALS(riemann_solver,"HLLE")) then

    if (use_cxx_code.eq.0) then
       call GRHydro_HLLE(CCTK_PASS_FTOF)
    else
       call GRHydro_HLLE_CC_F2C(cctkGH)
    endif

    if (evolve_tracer .ne. 0) then
    
      call GRHydro_HLLE_Tracer(CCTK_PASS_FTOF)

    end if
    
  else if (CCTK_EQUALS(riemann_solver,"HLLC")) then   

    call GRHydro_HLLC(CCTK_PASS_FTOF)

    if (evolve_tracer .ne. 0) then
    
      call GRHydro_HLLE_Tracer(CCTK_PASS_FTOF)

    end if
    
  else if (CCTK_EQUALS(riemann_solver,"Roe")) then   
    
    call GRHydro_RoeSolve(CCTK_PASS_FTOF)

    if (evolve_tracer .ne. 0) then
    
      call GRHydro_HLLE_Tracer(CCTK_PASS_FTOF)

    end if
    
  else if (CCTK_EQUALS(riemann_solver,"Marquina")) then   
    
    call GRHydro_Marquina(CCTK_PASS_FTOF)

!!$    Tracers are built directly in to the Marquina solver

  end if

end subroutine RiemannSolve

 /*@@
   @routine    RiemannSolvePolytype
   @date       Tue Mar 19 11:40:20 2002
   @author     Ian Hawke
   @desc 
   The same as above, just specializing to polytropic type EOS.
   Currently there is no point to this routine right now.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/


subroutine RiemannSolvePolytype(CCTK_ARGUMENTS)

  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  if (CCTK_EQUALS(riemann_solver,"HLLE")) then
    
    call GRHydro_HLLE(CCTK_PASS_FTOF)

    if (evolve_tracer .ne. 0) then
    
      call GRHydro_HLLE_Tracer(CCTK_PASS_FTOF)

    end if
    
  else if (CCTK_EQUALS(riemann_solver,"Roe")) then   
    
    call GRHydro_RoeSolve(CCTK_PASS_FTOF)

    if (evolve_tracer .ne. 0) then
    
      call GRHydro_HLLE_Tracer(CCTK_PASS_FTOF)

    end if
    
  else if (CCTK_EQUALS(riemann_solver,"Marquina")) then   
    
    call GRHydro_Marquina(CCTK_PASS_FTOF)

!!$    Tracers are built directly in to the Marquina solver

  end if

end subroutine RiemannSolvePolytype




