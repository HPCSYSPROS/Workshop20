/*@@
   @file      GRHydro_EvolutionMask.F90
   @date      Sat Jul 14 15:38:02 PDT 2012
   @author    Roland Haas
   @desc 
   User level module and Fortran glue code to get access to
   CarpetEvolutionMask::evolution_mask based on runtime parameters.
   @enddesc 
 @@*/
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

 /*@@
   @routine    GRHydro_EvolutionMask
   @date       Sat Jul 14 15:39:37 PDT 2012
   @author     Roland Haas
   @desc
   User level module to get access to CarpetEvolutionMask::evolution_mask based
   on runtime parameters.
   @history

   @endhistory

@@*/
module GRHydro_EvolutionMask

      implicit none

CONTAINS
       /*@@
         @routine    GRHydro_DeclareEvolutionMask
         @date       Sat Jul 14 15:39:37 PDT 2012
         @author     Roland Haas
         @desc
         Stores a pointer to CarpetEvolutionMask::evoltuion_mask in its
         argument. This function is not thread safe.
         @enddesc
         @calls     
         @calledby
         @history
      
         @endhistory

         @var     cctkGH
         @vdesc   Cactus grid hierarchy
         @vtype   cGH *
         @vio     in
         @vcomment
         @endvar

         @var     evolution_mask
         @vdesc   Cray pointer
         @vtype   CCTK_POINTER
         @vio     out
         @vcomment
         @endvar
      
         @var     evolution_mask_valid
         @vdesc   set to 1 if evolution_mask is valid, 0 otherwise
         @vtype   CCTK_INT
         @vio     out
         @vcomment
         @endvar
      
         @returntype none
         @returndesc
         @endreturndesc
      @@*/
      subroutine GRHydro_DeclareEvolutionMask(cctkGH, evolution_mask, &
                                              evolution_mask_valid)

        implicit none

        DECLARE_CCTK_PARAMETERS
        DECLARE_CCTK_FUNCTIONS

        CCTK_POINTER_TO_CONST :: cctkGH
        CCTK_REAL :: dummy
        pointer (evolution_mask, dummy)
        CCTK_INT :: evolution_mask_valid

        integer, save :: evolution_mask_idx = -1
        logical :: try_use_mask
        integer :: evolution_mask_active

        call CCTK_IsImplementationActive(evolution_mask_active, &
                                         "CarpetEvolutionMask")
        try_use_mask = CCTK_EQUALS(use_evolution_mask, "always") .or.   &
                       (CCTK_EQUALS(use_evolution_mask, "auto") .and. &
                        evolution_mask_active .ne. 0)

        if (try_use_mask) then
          if (evolution_mask_idx .eq. -1) then
            call CCTK_VarIndex(evolution_mask_idx,&
                               "CarpetEvolutionMask::evolution_mask")
          end if
          call CCTK_VarDataPtrI(evolution_mask, cctkGH, 0, evolution_mask_idx)
          if (evolution_mask .eq. CCTK_NullPointer()) then
            call CCTK_Warn(CCTK_WARN_ABORT, "Could not get pointer to evolution_mask. Is CarpetEvolutionMask active?")
          end if
          evolution_mask_valid = 1
        else
          evolution_mask = CCTK_NullPointer()
          evolution_mask_valid = 0
        end if

      end subroutine
end module
