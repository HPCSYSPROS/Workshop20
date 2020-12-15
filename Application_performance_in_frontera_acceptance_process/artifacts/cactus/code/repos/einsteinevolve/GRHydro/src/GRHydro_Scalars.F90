 /*@@
   @file      GRHydro_Scalars.F90
   @date      Mon Feb 25 11:11:23 2002
   @author    
   @desc 
   Module containing various scalars to avoid having to use
   CCTK_EQUALS on keywords at every cell (slow).
   @enddesc 
 @@*/

#include "cctk.h"

 module GRHydro_Scalars

   implicit none

   LOGICAL, SAVE :: MINMOD, MC2, SUPERBEE, PPM3, PPM4
   LOGICAL, SAVE :: ANALYTICAL
   LOGICAL, SAVE :: FAST
   LOGICAL, SAVE :: HLLE, LLF

   
 end module GRHydro_Scalars
