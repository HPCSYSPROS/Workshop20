 /*@@
   @file      GRHydro_WENOScalars.F90
   @date      Fri Jan  3 2013
   @author    Ian Hawke, Christian Reisswig
   @desc 
   Module containing the coefficient array for WENO reconstruction.
   @enddesc 
 @@*/

#include "cctk.h"

 module GRHydro_WENOScalars

   implicit none

   CCTK_REAL, allocatable, save :: weno_coeffs(:, :)
   CCTK_REAL, allocatable, save :: beta_shu(:, :)
   
 end module GRHydro_WENOScalars
