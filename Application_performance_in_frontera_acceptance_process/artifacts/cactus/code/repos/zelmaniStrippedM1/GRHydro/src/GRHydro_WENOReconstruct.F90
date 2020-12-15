 /*@@
   @file      GRHydro_WENOReconstruct.F90
   @date      Fri Jan 3 2013
   @author    Ian Hawke, Christian Reisswig
   @desc 
   Routines to set up the coefficient array and to perform one dimensional 
   ENO reconstruction of arbitrary order.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

 /*@@
   @routine    GRHydro_WENOSetup
   @date       Fri Jan 3 2013
   @author     Christian Reisswig
   @desc 
   Sets up the coefficient array for WENO reconstruction. 
   Uses the notation of Shu, equation (2.21), in 
   ''High Order ENO and WENO Schemes for CFD''
   (see ThornGuide for full reference).
   One exception: (Shu) r -> (Here) i: avoiding name clash.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine GRHydro_WENOSetup(CCTK_ARGUMENTS)

  USE GRHydro_WENOScalars

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: allocstat

  if(.not.allocated(weno_coeffs)) then ! uses weno_coeffs as a sentinel
     ! Right now we hardcode to 5th order
     allocate(weno_coeffs(3, 5), STAT=allocstat)
     if (allocstat .ne. 0) call CCTK_WARN(0, "Failed to allocate WENO coefficient arrays!")
     allocate(beta_shu(3, 6), STAT=allocstat)
     if (allocstat .ne. 0) call CCTK_WARN(0, "Failed to allocate smoothness indicator stencil coefficient arrays!")
  endif

  ! Set stencils
  if (CCTK_EQUALS(recon_method,"weno-z")) then
      weno_coeffs(1,1) =   2.0d0/6.0d0
      weno_coeffs(1,2) =  -7.0d0/6.0d0
      weno_coeffs(1,3) =  11.0d0/6.0d0
      weno_coeffs(1,4) =   0.0d0
      weno_coeffs(1,5) =   0.0d0
      
      weno_coeffs(2,1) =   0.0d0
      weno_coeffs(2,2) =  -1.0d0/6.0d0
      weno_coeffs(2,3) =   5.0d0/6.0d0
      weno_coeffs(2,4) =   2.0d0/6.0d0
      weno_coeffs(2,5) =   0.0d0
      
      weno_coeffs(3,1) =   0.0d0
      weno_coeffs(3,2) =   0.0d0
      weno_coeffs(3,3) =   2.0d0/6.0d0
      weno_coeffs(3,4) =   5.0d0/6.0d0
      weno_coeffs(3,5) =  -1.0d0/6.0d0
  else
      weno_coeffs(1,1) =   3.0d0/8.0d0
      weno_coeffs(1,2) =  -5.0d0/4.0d0
      weno_coeffs(1,3) =  15.0d0/8.0d0
      weno_coeffs(1,4) =   0.0d0
      weno_coeffs(1,5) =   0.0d0
      
      weno_coeffs(2,1) =   0.0d0
      weno_coeffs(2,2) =  -1.0d0/8.0d0
      weno_coeffs(2,3) =   3.0d0/4.0d0
      weno_coeffs(2,4) =   3.0d0/8.0d0
      weno_coeffs(2,5) =   0.0d0
      
      weno_coeffs(3,1) =   0.0d0
      weno_coeffs(3,2) =   0.0d0
      weno_coeffs(3,3) =   3.0d0/8.0d0
      weno_coeffs(3,4) =   3.0d0/4.0d0
      weno_coeffs(3,5) =  -1.0d0/8.0d0
  endif
   
  ! Shu smoothness indicator stencil coefficients
  beta_shu(1,1) =   4.0d0/3.0d0
  beta_shu(1,2) = -19.0d0/3.0d0
  beta_shu(1,3) =  25.0d0/3.0d0
  beta_shu(1,4) =  11.0d0/3.0d0
  beta_shu(1,5) = -31.0d0/3.0d0
  beta_shu(1,6) =  10.0d0/3.0d0
  
  beta_shu(2,1) =   4.0d0/3.0d0
  beta_shu(2,2) = -13.0d0/3.0d0
  beta_shu(2,3) =  13.0d0/3.0d0
  beta_shu(2,4) =   5.0d0/3.0d0
  beta_shu(2,5) = -13.0d0/3.0d0
  beta_shu(2,6) =   4.0d0/3.0d0
  
  beta_shu(3,1) =  10.0d0/3.0d0
  beta_shu(3,2) = -31.0d0/3.0d0
  beta_shu(3,3) =  25.0d0/3.0d0
  beta_shu(3,4) =  11.0d0/3.0d0
  beta_shu(3,5) = -19.0d0/3.0d0
  beta_shu(3,6) =   4.0d0/3.0d0
  
end subroutine GRHydro_WENOSetup

 /*@@
   @routine    GRHydro_ENOShutdown
   @date       Fri Jan  3 2013
   @author     Christian Reisswig
   @desc 
   Deallocates the coefficient arrays
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine GRHydro_WENOShutdown(CCTK_ARGUMENTS)

  USE GRHydro_WENOScalars

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: deallocstat

  if(allocated(weno_coeffs)) then
     deallocate(weno_coeffs, STAT = deallocstat)
     if (deallocstat .ne. 0) call CCTK_WARN(0, "Failed to deallocate WENO coefficients.")
     deallocate(beta_shu, STAT = deallocstat)
     if (deallocstat .ne. 0) call CCTK_WARN(0, "Failed to deallocate shu smoothness indicator coefficients.")
  endif

end subroutine GRHydro_WENOShutdown

 /*@@
   @routine    GRHydro_ENOReconstruct1d
   @date       Fri Jan 3 2013
   @author     Christian Reisswig
   @desc 
   Perform a one dimensional reconstruction of a given array using WENO.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

#define SpaceMask_CheckStateBitsF90_1D(mask,i,type_bits,state_bits) \
  (iand(mask((i)),(type_bits)).eq.(state_bits))

subroutine GRHydro_WENOReconstruct1d(order, nx, v, vminus, vplus, trivial_rp, &
     hydro_excision_mask)
  USE GRHydro_WENOScalars

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  CCTK_INT :: order, nx, i, j
  CCTK_REAL, dimension(nx) :: v, vplus, vminus

  CCTK_INT, dimension(nx) :: hydro_excision_mask
  logical, dimension(nx) :: trivial_rp
  logical, dimension(nx) :: excise
  logical :: normal_weno
  logical :: wenoZ

  CCTK_REAL :: beta1, beta2, beta3, vnorm, betanorm
  CCTK_REAL :: wplus1, wplus2, wplus3, wbarplus1, wbarplus2, wbarplus3
  CCTK_REAL :: wminus1, wminus2, wminus3, wbarminus1, wbarminus2, wbarminus3

  ! driver initialized these with v already
  !vminus = 0.d0
  !vplus = 0.d0

  excise = .false.
  trivial_rp = .false.

  if (CCTK_EQUALS(recon_method,"weno-z")) then
     wenoZ = .true.
  else
     wenoZ = .false.
  endif
  
!!$    Initialize excision
  do i = 1, nx
    if (GRHydro_enable_internal_excision /= 0 .and. (hydro_excision_mask(i) .ne. 0)) then
      trivial_rp(i) = .true.
      excise(i) = .true.
      if (i > 1) then
        trivial_rp(i-1) = .true.
      end if
    end if
  end do

  do i = 3, nx-2
!!$    Handle excision
    normal_weno = .true.
    if (i < nx) then
     if (excise(i+1)) then
      vminus(i) = v(i)
      vplus(i) = v(i)
      normal_weno = .false.
     end if
    end if
    if (i > 1) then
     if (excise(i-1)) then
      vminus(i) = v(i)
      vplus(i) = v(i)
      normal_weno = .false.
     end if
    end if

    if (normal_weno) then

      if (wenoZ) then

         beta1 = 13.0d0/12.d0*((v(i-2)-2.0d0*v(i-1))+v(i))**2 + 1.0d0/4.d0*((v(i-2)-4.0d0*v(i-1))+3.0d0*v(i))**2
         beta2 = 13.0d0/12.d0*((v(i-1)+v(i+1))-2.0d0*v(i))**2 + 1.0d0/4.d0*(v(i-1)-v(i+1))**2
         beta3 = 13.0d0/12.d0*((v(i+2)-2.0d0*v(i+1))+v(i))**2 + 1.0d0/4.d0*((v(i+2)-4.0d0*v(i+1))+3.0d0*v(i))**2
      
      
      !!$    Compute weights according to weno-z alorithm
         wbarplus1 = 1.0d0/10.0d0 * (1.0d0 + abs(beta1-beta3) / (weno_eps + beta1))
         wbarplus2 = 3.0d0/5.0d0 * (1.0d0 + abs(beta1-beta3) / (weno_eps + beta2))
         wbarplus3 = 3.0d0/10.0d0 * (1.0d0 + abs(beta1-beta3) / (weno_eps + beta3))

         wbarminus1 = 3.0d0/10.0d0 * (1.0d0 + abs(beta1-beta3) / (weno_eps + beta1))
         wbarminus2 = 3.0d0/5.0d0 * (1.0d0 + abs(beta1-beta3) / (weno_eps + beta2))
         wbarminus3 = 1.0d0/10.0d0 * (1.0d0 + abs(beta1-beta3) / (weno_eps + beta3))

      
      else
      
   !!$    Compute smoothness indicators
         beta1 = (13.0d0/12.d0)*((v(i-2)-2.0d0*v(i-1))+v(i))**2 + (1.0d0/4.d0)*((v(i-2)-4.0d0*v(i-1))+3.0d0*v(i))**2
         beta2 = (13.0d0/12.d0)*((v(i-1)+v(i+1))-2.0d0*v(i))**2 + (1.0d0/4.d0)*(v(i-1)-v(i+1))**2
         beta3 = (13.0d0/12.d0)*((v(i+2)-2.0d0*v(i+1))+v(i))**2 + (1.0d0/4.d0)*((v(i+2)-4.0d0*v(i+1))+3.0d0*v(i))**2
      
   !!$    This is modification is suggested by Tchekhovskoy et al 2007 (WHAM code paper).
         if (weno_adaptive_epsilon.ne.0) then
            vnorm = ((v(i-2)**2 + v(i+2)**2) + (v(i-1)**2 + v(i+1)**2) + v(i)**2)
            
            beta1 = beta1 + 100.0d0*weno_eps*(vnorm + 1.0d0)
            beta2 = beta2 + 100.0d0*weno_eps*(vnorm + 1.0d0)
            beta3 = beta3 + 100.0d0*weno_eps*(vnorm + 1.0d0)
            
            betanorm = beta1 + beta3 + beta2
            
            beta1 = beta1 / betanorm
            beta2 = beta2 / betanorm
            beta3 = beta3 / betanorm
         endif
         
         wbarplus1 = 1.0d0/16.0d0 / (weno_eps + beta1)**2
         wbarplus2 = 5.0d0/8.0d0 / (weno_eps + beta2)**2
         wbarplus3 = 5.0d0/16.0d0 / (weno_eps + beta3)**2
      
         wbarminus1 = 5.0d0/16.0d0 / (weno_eps + beta1)**2
         wbarminus2 = 5.0d0/8.0d0 / (weno_eps + beta2)**2
         wbarminus3 = 1.0d0/16.0d0 / (weno_eps + beta3)**2

      endif
      
      wplus1 = wbarplus1 / ((wbarplus1 + wbarplus3) + wbarplus2)
      wplus2 = wbarplus2 / ((wbarplus1 + wbarplus3) + wbarplus2)
      wplus3 = wbarplus3 / ((wbarplus1 + wbarplus3) + wbarplus2)
    
      wminus1 = wbarminus1 / ((wbarminus1 + wbarminus3) + wbarminus2)
      wminus2 = wbarminus2 / ((wbarminus1 + wbarminus3) + wbarminus2)
      wminus3 = wbarminus3 / ((wbarminus1 + wbarminus3) + wbarminus2)
            
!!$    Calculate the reconstruction
      vplus(i) = 0d0
      vminus(i) = 0d0
      do j = 1, 5
         vplus(i) = vplus(i) + ((wplus1 * weno_coeffs(1,j) &
                               + wplus2 * weno_coeffs(2,j)) &
                               + wplus3 * weno_coeffs(3,j))*v(i-3+j)
         vminus(i) = vminus(i) + ((wminus3 * weno_coeffs(1,j) &
                                 + wminus2 * weno_coeffs(2,j)) &
                                 + wminus1 * weno_coeffs(3,j))*v(i-3+6-j)
      end do
      !vminus(i) = v(i)
      !vplus(i) = v(i)
            
    end if
  end do

  do i = 1, nx
    if (excise(i)) then
      if (i > 1) then
        if (.not. excise(i-1)) then
          vminus(i) = vplus(i-1)
        end if
      end if
      vplus(i) = vminus(i)
    end if
  end do
  do i = nx, 1, -1
    if (excise(i)) then
      if (i < nx) then
        if (.not. excise(i+1)) then
          vplus(i) = vminus(i+1)
        end if
      end if
      vminus(i) = vplus(i)
    end if
  end do

end subroutine GRHydro_WENOReconstruct1d
