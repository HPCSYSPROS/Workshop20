 /*@@
   @file      GRHydro_ENOScalars.F90
   @date      Sat Apr  6 17:35:42 2002
   @author    Ian Hawke
   @desc 
   Module containing the coefficient array for ENO reconstruction.
   @enddesc 
 @@*/

#include "cctk.h"

 module GRHydro_ENOScalars

   implicit none

   CCTK_REAL, allocatable, save :: eno_coeffs(:, :)
   
 end module GRHydro_ENOScalars
