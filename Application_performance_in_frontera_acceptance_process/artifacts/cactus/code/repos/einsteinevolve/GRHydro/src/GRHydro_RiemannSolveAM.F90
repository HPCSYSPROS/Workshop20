 /*@@
   @file      GRHydro_RiemannSolveAM.F90
   @date      Sep 1, 2010
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
   @routine    RiemannSolveAM
   @date       Sep 1, 2010
   @author     Joshua Faber, Scott Noble, Bruno Mundim, Pedro Montero, Ian Hawke
   @desc 
   A wrapper routine to switch between the different Riemann solvers.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine RiemannSolveAM(CCTK_ARGUMENTS)

  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  if (CCTK_EQUALS(riemann_solver,"HLLE").or.CCTK_EQUALS(riemann_solver,"LLF")) then
     
     call GRHydro_HLLE_AM(CCTK_PASS_FTOF)
     
     if (evolve_tracer .ne. 0) then
        
!!$ There are no special calls for tracers, which care not one whit about B-fields!    
!!$ Just call the standard version...
        
        call GRHydro_HLLE_Tracer(CCTK_PASS_FTOF)
        
     end if
     
!!$  else if (CCTK_EQUALS(riemann_solver,"Roe")) then   
!!$    
!!$    call GRHydro_RoeSolveAM(CCTK_PASS_FTOF)
!!$
!!$    if (evolve_tracer .ne. 0) then
!!$    
!!$      call GRHydro_HLLE_Tracer(CCTK_PASS_FTOF)
!!$
!!$    end if
!!$    
!!$  else if (CCTK_EQUALS(riemann_solver,"Marquina")) then   
!!$   
!!$    call GRHydro_MarquinaM(CCTK_PASS_FTOF)
     
!!$    Tracers are built directly in to the Marquina solver
     
  else
     
     call CCTK_ERROR("Roe and Marquina not implemented in MHD yet!!!")
     STOP
     
  end if
  
end subroutine RiemannSolveAM

 /*@@
   @routine    RiemannSolvePolytypeAM
   @date       Sep 1, 2010
   @author     Joshua Faber, Scott Noble, Bruno Mundim, Ian Hawke
   @desc 
   The same as above, just specializing to polytropic type EOS.
   Currently there is no point to this routine right now.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/


subroutine RiemannSolvePolytypeAM(CCTK_ARGUMENTS)

  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  if (CCTK_EQUALS(riemann_solver,"HLLE")) then
    
     call GRHydro_HLLE_AM(CCTK_PASS_FTOF)
     
     if (evolve_tracer .ne. 0) then
        
!!$ Call the non-MHD version - see above
        
        call GRHydro_HLLE_Tracer(CCTK_PASS_FTOF)
        
     end if
     
!!$  else if (CCTK_EQUALS(riemann_solver,"Roe")) then   
!!$    
!!$    call GRHydro_RoeSolve(CCTK_PASS_FTOF)
!!$
!!$    if (evolve_tracer .ne. 0) then
!!$    
!!$      call GRHydro_HLLE_Tracer(CCTK_PASS_FTOF)
!!$
!!$    end if
!!$    
!!$  else if (CCTK_EQUALS(riemann_solver,"Marquina")) then   
!!$    
!!$    call GRHydro_Marquina(CCTK_PASS_FTOF)
     
!!$    Tracers are built directly in to the Marquina solver
     
  else
     
     call CCTK_ERROR("Roe and Marquina not implemented in MHD yet!!!")
     STOP

  end if
  
end subroutine RiemannSolvePolytypeAM





