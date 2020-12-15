 /*@@
   @file      GRHydro_Boundaries.F90
   @date      Sat Jan 26 01:01:14 2002
   @author    
   @desc 
   The two routines for dealing with boundary conditions.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "GRHydro_Macros.h"

#include "util_Table.h"

 /*@@
   @routine    GRHydro_InitSymBound
   @date       Sat Jan 26 01:03:04 2002
   @author     Ian Hawke
   @desc 
   Sets up the symmetries at the boundaries of the hydrodynamical variables.
   @enddesc 
   @calls     
   @calledby   
   @history 
   Direct translation of routines from GR3D, GRAstro_Hydro, 
   written by Mark Miller, or WaveToy routines, or...
   @endhistory 

@@*/

subroutine GRHydro_InitSymBound(CCTK_ARGUMENTS)

  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer, dimension(3) :: sym
  integer :: ierr
  CCTK_INT :: GRHydro_UseGeneralCoordinates, general_coordinates

  ierr = 0
  sym(1) = 1
  sym(2) = 1
  sym(3) = 1

  general_coordinates = GRHydro_UseGeneralCoordinates(cctkGH)
  
  !if (sync_conserved_only .eq. 0) then
   call SetCartSymVN(ierr, cctkGH, sym, "HydroBase::rho")
   call SetCartSymVN(ierr, cctkGH, sym, "HydroBase::press")
   call SetCartSymVN(ierr, cctkGH, sym, "HydroBase::w_lorentz")
   call SetCartSymVN(ierr, cctkGH, sym, "HydroBase::eps")
  !endif
  call SetCartSymVN(ierr, cctkGH, sym, "GRHydro::dens")
  call SetCartSymVN(ierr, cctkGH, sym, "GRHydro::tau")
  call SetCartSymVN(ierr, cctkGH, sym, "GRHydro::GRHydro_C2P_failed")
  if(evolve_mhd.ne.0.and.clean_divergence.ne.0) then
     call SetCartSymVN(ierr, cctkGH, sym, "GRHydro::psidc")
  endif

!!$ handle multiple tracer variables
  if(evolve_tracer.ne.0) then
     call SetCartSymGN(ierr, cctkGH, sym, "GRHydro::GRHydro_tracers")
     call SetCartSymGN(ierr, cctkGH, sym, "GRHydro::GRHydro_cons_tracers")
  endif 

  if(evolve_y_e.ne.0) then
   !if (sync_conserved_only .eq. 0) then
     call SetCartSymGN(ierr, cctkGH, sym, "HydroBase::Y_e")
   !endif
     call SetCartSymGN(ierr, cctkGH, sym, "GRHydro::Y_e_con")
  endif

  if(evolve_temper.ne.0) then
    !if (sync_conserved_only .eq. 0) then
      call SetCartSymGN(ierr, cctkGH, sym, "HydroBase::temperature")
      call SetCartSymGN(ierr, cctkGH, sym, "HydroBase::entropy")
    !endif
  endif
  
  call SetCartSymVN(ierr, cctkGH, sym, "GRHydro::atmosphere_mask_real")
  
  sym(1) = -1
  sym(2) = 1
  sym(3) = 1
  
  !if (sync_conserved_only .eq. 0) then
      if (general_coordinates .ne. 0) then
         call SetCartSymVN(ierr, cctkGH, sym, "GRHydro::lvel[0]")
      else
         call SetCartSymVN(ierr, cctkGH, sym, "HydroBase::vel[0]")
      endif
  !endif
  call SetCartSymVN(ierr, cctkGH, sym, "GRHydro::scon[0]")
  if(evolve_mhd.ne.0) then
     !if (sync_conserved_only .eq. 0) then
      call SetCartSymVN(ierr, cctkGH, sym, "HydroBase::Bvec[0]")
     !endif
     call SetCartSymVN(ierr, cctkGH, sym, "GRHydro::Bcons[0]")
  endif
  
  sym(1) = 1
  sym(2) = -1
  sym(3) = 1
  
  !if (sync_conserved_only .eq. 0) then
      if (general_coordinates .ne. 0) then
         call SetCartSymVN(ierr, cctkGH, sym, "GRHydro::lvel[1]")
      else
         call SetCartSymVN(ierr, cctkGH, sym, "HydroBase::vel[1]")
      endif
  !endif
  call SetCartSymVN(ierr, cctkGH, sym, "GRHydro::scon[1]")
  if(evolve_mhd.ne.0) then
     !if (sync_conserved_only .eq. 0) then
         call SetCartSymVN(ierr, cctkGH, sym, "HydroBase::Bvec[1]")
     !endif
     call SetCartSymVN(ierr, cctkGH, sym, "GRHydro::Bcons[1]")
  endif

  sym(1) = 1
  sym(2) = 1
  sym(3) = -1
  
  !if (sync_conserved_only .eq. 0) then
      if (general_coordinates .ne. 0) then
         call SetCartSymVN(ierr, cctkGH, sym, "GRHydro::lvel[2]")
      else
         call SetCartSymVN(ierr, cctkGH, sym, "HydroBase::vel[2]")
      endif
  !endif
  call SetCartSymVN(ierr, cctkGH, sym, "GRHydro::scon[2]")
  if(evolve_mhd.ne.0) then
    !if (sync_conserved_only .eq. 0) then
      call SetCartSymVN(ierr, cctkGH, sym, "HydroBase::Bvec[2]")
    !endif 
    call SetCartSymVN(ierr, cctkGH, sym, "GRHydro::Bcons[2]")
  endif

! check that storage for shift is active
 if(shift_state.eq.0) then
    call CCTK_WARN(0,"shift_state = 0 (no shift storage) no longer supported!");
  endif

  
end subroutine GRHydro_InitSymBound

 /*@@
   @routine    GRHydro_Boundaries
   @date       Sat Jan 26 01:04:04 2002
   @author     
   @desc 
   Calls the appropriate boundary routines
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine GRHydro_Boundaries(CCTK_ARGUMENTS)
  
  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  integer, dimension(3) :: sw
  integer :: ierr
  CCTK_INT :: GRHydro_UseGeneralCoordinates, general_coordinates

  CCTK_INT, parameter :: faces=CCTK_ALL_FACES
  CCTK_INT, parameter :: ione=1
  
  ierr = 0 

  sw = GRHydro_stencil

  general_coordinates = GRHydro_UseGeneralCoordinates(cctkGH)

  if (verbose.eq.1) call CCTK_INFO("Selecting conserved BC (and primitive BC if selected)")

!!$Flat boundaries if required  

  if (CCTK_EQUALS(bound,"flat")) then
    ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
         "GRHydro::dens", "Flat")
    ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
         "GRHydro::tau", "Flat")
    ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
         "GRHydro::scon", "Flat")
    if (sync_conserved_only .eq. 0) then
      ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
            "HydroBase::w_lorentz", "Flat")
      ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
            "HydroBase::rho", "Flat")
      ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
            "HydroBase::press", "Flat")
      ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
            "HydroBase::eps", "Flat")
      if (general_coordinates .ne. 0) then
         ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
               "GRHydro::lvel", "Flat")
      else
         ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
               "HydroBase::vel", "Flat")
      endif
    endif
    if(evolve_mhd.ne.0) then
       if (sync_conserved_only .eq. 0) then
         if (general_coordinates .ne. 0) then
            ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
                  "GRHydro::lBvec", "Flat")
         else
            ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
                  "HydroBase::Bvec", "Flat")
         endif
       endif
       ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
            "GRHydro::Bcons", "Flat")
       if(clean_divergence.ne.0) then
          ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
               "GRHydro::psidc", "Flat")
       endif
    endif

    if (CCTK_EQUALS(Bvec_evolution_method, "GRHydro_Avec")) then 
      ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
             "HydroBase::Avec", "Flat")
    endif 

    if(evolve_tracer.ne.0) then 
      if (sync_conserved_only .eq. 0) then
         ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
               "GRHydro::GRHydro_tracers", "Flat")
      endif
      ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
            "GRHydro::GRHydro_cons_tracers", "Flat")
    endif

    if(evolve_y_e.ne.0) then
      if (sync_conserved_only .eq. 0) then
         ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
               "HydroBase::Y_e", "Flat")
      endif
      ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
            "GRHydro::Y_e_con", "Flat")
    endif

    if(evolve_temper.ne.0) then
      if (sync_conserved_only .eq. 0) then
         ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
               "HydroBase::temperature", "Flat")
         ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
               "HydroBase::entropy", "Flat")
      endif
    endif


  endif

  if (CCTK_EQUALS(bound,"none")) then
    ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
         "GRHydro::dens", "None")
    ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
         "GRHydro::tau", "None")
    ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
         "GRHydro::scon", "None")
    
    if (sync_conserved_only .eq. 0) then
      ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
            "HydroBase::w_lorentz", "None")
      ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
            "HydroBase::rho", "None")
      ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
            "HydroBase::press", "None")
      ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
            "HydroBase::eps", "None")
      if (general_coordinates .ne. 0) then
         ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
               "GRHydro::lvel", "None")
      else
         ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
               "HydroBase::vel", "None")
      endif
    endif
    
    if(evolve_mhd.ne.0) then
       if (sync_conserved_only .eq. 0) then
         if (general_coordinates .ne. 0) then
            ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
                  "GRHydro::lBvec", "None")
         else
            ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
                  "HydroBase::Bvec", "None")
         endif
       endif
       ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
            "GRHydro::Bcons", "None")
       if(clean_divergence.ne.0) then
          ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
               "GRHydro::psidc", "None")
       endif
    endif

    if (CCTK_EQUALS(Bvec_evolution_method, "GRHydro_Avec")) then 
      ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
             "HydroBase::Avec", "None")
      if (general_coordinates .ne. 0) then
        ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
               "GRHydro::lBvec", "Flat")
      else
        ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
               "HydroBase::Bvec", "Flat")
      endif
    endif 

    if(evolve_tracer.ne.0) then 
      if (sync_conserved_only .eq. 0) then 
         ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
               "GRHydro::GRHydro_tracers", "None")
      endif
      ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
            "GRHydro::GRHydro_cons_tracers", "None")
    endif

    if(evolve_y_e.ne.0) then
      if (sync_conserved_only .eq. 0) then
         ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
               "HydroBase::Y_e", "None")
      endif
      ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
            "GRHydro::Y_e_con", "None")
    endif

    if(evolve_temper.ne.0) then
      if (sync_conserved_only .eq. 0) then
         ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
               "HydroBase::temperature", "None")
         ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
               "HydroBase::entropy", "None")
      endif
    endif

  endif

  if (CCTK_EQUALS(bound,"scalar")) then
    call CCTK_WARN(0, "Until somebody uses this I see no reason to support it")
  end if

  if (ierr < 0) call CCTK_WARN(0, "problems with applying the chosen boundary condition")

end subroutine GRHydro_Boundaries



subroutine GRHydro_SelectPrimitiveInitialGuessesBoundaries(CCTK_ARGUMENTS)
  
  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  integer, dimension(3) :: sw
  integer :: ierr
  CCTK_INT :: GRHydro_UseGeneralCoordinates, general_coordinates

  CCTK_INT, parameter :: faces=CCTK_ALL_FACES
  CCTK_INT, parameter :: ione=1
  
  ierr = 0
  sw = GRHydro_stencil

  general_coordinates = GRHydro_UseGeneralCoordinates(cctkGH)

!!$Flat boundaries if required  

! The commented out code are those primitives which do not require and explicit
! initial guess for Con2Prim. The guesses are either not needed or computed from other quantities!

  if (verbose.eq.1) call CCTK_INFO("Selecting primitive BC")


  if (CCTK_EQUALS(bound,"flat")) then
    
    !ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
    !      "HydroBase::w_lorentz", "Flat")
    ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
          "HydroBase::rho", "Flat")
    !ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
    !      "HydroBase::press", "Flat")
    ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
          "HydroBase::eps", "Flat")
    if (general_coordinates .ne. 0) then
       ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
             "GRHydro::lvel", "Flat")
    else
       ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
             "HydroBase::vel", "Flat")
    endif
    if(evolve_mhd.ne.0) then
         if (general_coordinates .ne. 0) then
            ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
                  "GRHydro::lBvec", "Flat")
         else
            ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
                  "HydroBase::Bvec", "Flat")
         endif
    endif

    !if(evolve_tracer.ne.0) then 
    !  ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
    !         "GRHydro::GRHydro_tracers", "Flat")
    !endif

    !if(evolve_y_e.ne.0) then
    !  ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
    !         "HydroBase::Y_e", "Flat")
    !endif

    if(evolve_temper.ne.0) then
       ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
             "HydroBase::temperature", "Flat")
       !ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
       !      "HydroBase::entropy", "Flat")
    endif


  endif

  if (CCTK_EQUALS(bound,"none")) then      
    !ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
    !      "HydroBase::w_lorentz", "None")
    ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
          "HydroBase::rho", "None")
    !ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
    !      "HydroBase::press", "None")
    ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
          "HydroBase::eps", "None")
    if (general_coordinates .ne. 0) then
       ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
             "GRHydro::lvel", "None")
    else
       ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
             "HydroBase::vel", "None")
    endif
    
    if(evolve_mhd.ne.0) then
         if (general_coordinates .ne. 0) then
            ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
                  "GRHydro::lBvec", "None")
         else
            ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
                  "HydroBase::Bvec", "None")
         endif
       
    endif

    !if(evolve_tracer.ne.0) then 
    !   ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
    !         "GRHydro::GRHydro_tracers", "None")
    !  
    !endif

    !if(evolve_y_e.ne.0) then
    !   ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
    !         "HydroBase::Y_e", "None")
    !  
    !endif

    if(evolve_temper.ne.0) then
       ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
             "HydroBase::temperature", "None")
       !ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
       !      "HydroBase::entropy", "None")
    endif

  end if

  if (CCTK_EQUALS(bound,"scalar")) then
    call CCTK_WARN(0, "Until somebody uses this I see no reason to support it")
  end if

  if (ierr < 0) call CCTK_WARN(0, "problems with applying the chosen boundary condition")

end subroutine GRHydro_SelectPrimitiveInitialGuessesBoundaries








subroutine GRHydro_SelectPrimitiveBoundaries(CCTK_ARGUMENTS)
  
  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  integer, dimension(3) :: sw
  integer :: ierr
  CCTK_INT :: GRHydro_UseGeneralCoordinates, general_coordinates

  CCTK_INT, parameter :: faces=CCTK_ALL_FACES
  CCTK_INT, parameter :: ione=1
  
  ierr = 0
  sw = GRHydro_stencil

  general_coordinates = GRHydro_UseGeneralCoordinates(cctkGH)

!!$Flat boundaries if required  

  if (verbose.eq.1) call CCTK_INFO("Selecting primitive BC")


  if (CCTK_EQUALS(bound,"flat")) then
    
    ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
          "HydroBase::w_lorentz", "Flat")
    ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
          "HydroBase::rho", "Flat")
    ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
          "HydroBase::press", "Flat")
    ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
          "HydroBase::eps", "Flat")
    if (general_coordinates .ne. 0) then
       ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
             "GRHydro::lvel", "Flat")
    else
       ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
             "HydroBase::vel", "Flat")
    endif
    if(evolve_mhd.ne.0) then
         if (general_coordinates .ne. 0) then
            ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
                  "GRHydro::lBvec", "Flat")
         else
            ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
                  "HydroBase::Bvec", "Flat")
         endif
    endif

    if(evolve_tracer.ne.0) then 
             ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
                   "GRHydro::GRHydro_tracers", "Flat")
    endif

    if(evolve_y_e.ne.0) then
             ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
                   "HydroBase::Y_e", "Flat")
    endif

    if(evolve_temper.ne.0) then
       ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
             "HydroBase::temperature", "Flat")
       ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
             "HydroBase::entropy", "Flat")
    endif


  endif

  if (CCTK_EQUALS(bound,"none")) then      
    ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
          "HydroBase::w_lorentz", "None")
    ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
          "HydroBase::rho", "None")
    ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
          "HydroBase::press", "None")
    ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
          "HydroBase::eps", "None")
    if (general_coordinates .ne. 0) then
       ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
             "GRHydro::lvel", "None")
    else
       ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
             "HydroBase::vel", "None")
    endif
    
    if(evolve_mhd.ne.0) then
         if (general_coordinates .ne. 0) then
            ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
                  "GRHydro::lBvec", "None")
         else
            ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
                  "HydroBase::Bvec", "None")
         endif
       
    endif

    if(evolve_tracer.ne.0) then 
       ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
             "GRHydro::GRHydro_tracers", "None")
      
    endif

    if(evolve_y_e.ne.0) then
       ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
             "HydroBase::Y_e", "None")
      
    endif

    if(evolve_temper.ne.0) then
       ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
             "HydroBase::temperature", "None")
       ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
             "HydroBase::entropy", "None")
    endif

  end if

  if (CCTK_EQUALS(bound,"scalar")) then
    call CCTK_WARN(0, "Until somebody uses this I see no reason to support it")
  end if

  if (ierr < 0) call CCTK_WARN(0, "problems with applying the chosen boundary condition")

end subroutine GRHydro_SelectPrimitiveBoundaries








subroutine GRHydro_SelectAtmosphereMaskBoundaries(CCTK_ARGUMENTS)
  
  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  integer, dimension(3) :: sw
  integer :: ierr

  CCTK_INT, parameter :: faces=CCTK_ALL_FACES
  CCTK_INT, parameter :: ione=1
  
  ierr = 0
  sw = GRHydro_stencil


!!$Flat boundaries if required  

  if (verbose.eq.1) call CCTK_INFO("Selecting atmosphere mask BC")


  if (CCTK_EQUALS(bound,"flat")) then
      ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
         "GRHydro::GRHydro_atmosphere_mask_real", "Flat")
  endif

  if (CCTK_EQUALS(bound,"none")) then      
      ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, GRHydro_stencil, -ione, &
         "GRHydro::GRHydro_atmosphere_mask_real", "None")
  end if

  if (CCTK_EQUALS(bound,"scalar")) then
    call CCTK_WARN(0, "Until somebody uses this I see no reason to support it")
  end if

  if (ierr < 0) call CCTK_WARN(0, "problems with applying the chosen boundary condition")

end subroutine GRHydro_SelectAtmosphereMaskBoundaries

