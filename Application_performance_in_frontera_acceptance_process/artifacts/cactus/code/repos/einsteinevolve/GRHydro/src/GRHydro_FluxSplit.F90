 /*@@
   @file      GRHydro_FluxSplit.F90
   @date      Wed Mar  3 22:16:00 2004
   @author    Ian Hawke
   @desc 
   Flux split reconstruction using WENO5.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"
#include "GRHydro_Macros.h"

 /*@@
   @routine    GRHydro_FSAlpha
   @date       Sun Jan  9 12:25:43 2005
   @author     Ian Hawke
   @desc 
   Compute the maximum characteristic speed over all space
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/


subroutine GRHydro_FSAlpha(CCTK_ARGUMENTS)

  use GRHydro_Eigenproblem

  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer :: nx, ny, nz, i, j, k, ierr, max_handle
  CCTK_REAL :: uxx, uxy, uxz, uyy, uyz, uzz, sdet, beta
  CCTK_REAL, dimension(5) :: lambda
  CCTK_REAL :: alpha1_local, alpha2_local, alpha3_local, alpha4_local, &
       alpha5_local

  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)

  fs_alpha1 = 0.d0
  fs_alpha2 = 0.d0
  fs_alpha3 = 0.d0
  fs_alpha4 = 0.d0
  fs_alpha5 = 0.d0

  alpha1_local = 0.d0
  alpha2_local = 0.d0
  alpha3_local = 0.d0
  alpha4_local = 0.d0
  alpha5_local = 0.d0

  if (flux_direction == 1) then

    do k = 1, nz
      do j = 1, ny
        do i = 1, nx

           sdet = sdetg(i,j,k)
           call UpperMetric(uxx,uxy,uxz,uyy,uyz,uzz,sdet*sdet,&
                gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),&
                gyz(i,j,k),gzz(i,j,k))        

           beta = betax(i,j,k)
          
          call eigenvalues(GRHydro_eos_handle, &
               rho      (i,j,k), &
               vel      (i,j,k,1), &
               vel      (i,j,k,2), &
               vel      (i,j,k,3), &
               eps      (i,j,k), &
               w_lorentz(i,j,k), &
               lambda          , &
               gxx      (i,j,k), &
               gxy      (i,j,k), &
               gxz      (i,j,k), &
               gyy      (i,j,k), &
               gyz      (i,j,k), &
               gzz      (i,j,k), &
               uxx             , &
               alp      (i,j,k), &
               beta               )  

          alpha1_local = max(abs(lambda(1)), alpha1_local)
          alpha2_local = max(abs(lambda(2)), alpha2_local)
          alpha3_local = max(abs(lambda(3)), alpha3_local)
          alpha4_local = max(abs(lambda(4)), alpha4_local)
          alpha5_local = max(abs(lambda(5)), alpha5_local)

        end do
      end do
    end do

  else if (flux_direction == 2) then

    do k = 1, nz
      do j = 1, ny
        do i = 1, nx

           sdet = sdetg(i,j,k)
           call UpperMetric(uxx,uxy,uxz,uyy,uyz,uzz,sdet*sdet,&
                gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),&
                gyz(i,j,k),gzz(i,j,k))        

           beta = betay(i,j,k)
          
          call eigenvalues(GRHydro_eos_handle, &
               rho      (i,j,k), &
               vel      (i,j,k,1), &
               vel      (i,j,k,2), &
               vel      (i,j,k,3), &
               eps      (i,j,k), &
               w_lorentz(i,j,k), &
               lambda          , &
               gyy      (i,j,k), &
               gyz      (i,j,k), &
               gxy      (i,j,k), &
               gzz      (i,j,k), &
               gxz      (i,j,k), &
               gxx      (i,j,k), &
               uyy             , &
               alp      (i,j,k), &
               beta               )  

          alpha1_local = max(abs(lambda(1)), alpha1_local)
          alpha2_local = max(abs(lambda(2)), alpha2_local)
          alpha3_local = max(abs(lambda(3)), alpha3_local)
          alpha4_local = max(abs(lambda(4)), alpha4_local)
          alpha5_local = max(abs(lambda(5)), alpha5_local)

        end do
      end do
    end do

  else if (flux_direction == 3) then

    do k = 1, nz
      do j = 1, ny
        do i = 1, nx

           sdet = sdetg(i,j,k)
           call UpperMetric(uxx,uxy,uxz,uyy,uyz,uzz,sdet*sdet,&
                gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),&
                gyz(i,j,k),gzz(i,j,k))        

           beta = betaz(i,j,k)
          
          call eigenvalues(GRHydro_eos_handle, &
               rho      (i,j,k), &
               vel      (i,j,k,1), &
               vel      (i,j,k,2), &
               vel      (i,j,k,3), &
               eps      (i,j,k), &
               w_lorentz(i,j,k), &
               lambda          , &
               gzz      (i,j,k), &
               gxz      (i,j,k), &
               gyz      (i,j,k), &
               gxx      (i,j,k), &
               gxy      (i,j,k), &
               gyy      (i,j,k), &
               uzz             , &
               alp      (i,j,k), &
               beta               )  

          alpha1_local = max(abs(lambda(1)), alpha1_local)
          alpha2_local = max(abs(lambda(2)), alpha2_local)
          alpha3_local = max(abs(lambda(3)), alpha3_local)
          alpha4_local = max(abs(lambda(4)), alpha4_local)
          alpha5_local = max(abs(lambda(5)), alpha5_local)

        end do
      end do
    end do

  else

    call CCTK_ERROR("Flux direction not x,y,z")
    STOP

  end if

!!$  write(*,*) 'fs_alpha: local',alpha1_local,alpha2_local,alpha3_local,alpha4_local,alpha5_local
  
  call CCTK_ReductionHandle(max_handle, "maximum")

  call CCTK_ReduceLocScalar ( ierr, cctkGH, -1, max_handle, &
       alpha1_local, fs_alpha1, CCTK_VARIABLE_REAL )
  if ( ierr .ne. 0 ) then
    call CCTK_ERROR("Reduction of alpha1 failed")
    STOP
  end if
  call CCTK_ReduceLocScalar ( ierr, cctkGH, -1, max_handle, &
       alpha2_local, fs_alpha2, CCTK_VARIABLE_REAL )
  if ( ierr .ne. 0 ) then
    call CCTK_ERROR("Reduction of alpha2 failed")
    STOP
  end if
  call CCTK_ReduceLocScalar ( ierr, cctkGH, -1, max_handle, &
       alpha3_local, fs_alpha3, CCTK_VARIABLE_REAL )
  if ( ierr .ne. 0 ) then
    call CCTK_ERROR("Reduction of alpha3 failed")
    STOP
  end if
  call CCTK_ReduceLocScalar ( ierr, cctkGH, -1, max_handle, &
       alpha4_local, fs_alpha4, CCTK_VARIABLE_REAL )
  if ( ierr .ne. 0 ) then
    call CCTK_ERROR("Reduction of alpha4 failed")
    STOP
  end if
  call CCTK_ReduceLocScalar ( ierr, cctkGH, -1, max_handle, &
       alpha5_local, fs_alpha5, CCTK_VARIABLE_REAL )
  if ( ierr .ne. 0 ) then
    call CCTK_ERROR("Reduction of alpha5 failed")
    STOP
  end if

!!$  write(*,*) 'fs_alpha: global',fs_alpha1,fs_alpha2,fs_alpha3,fs_alpha4,fs_alpha5

!!$  fs_alpha1 = max(abs(fs_alpha1), &
!!$                  abs(fs_alpha2), &
!!$                  abs(fs_alpha3), &
!!$                  abs(fs_alpha4), &
!!$                  abs(fs_alpha5) )
!!$
!!$  fs_alpha2 = max(abs(fs_alpha1), &
!!$                  abs(fs_alpha2), &
!!$                  abs(fs_alpha3), &
!!$                  abs(fs_alpha4), &
!!$                  abs(fs_alpha5) )
!!$
!!$  fs_alpha3 = max(abs(fs_alpha1), &
!!$                  abs(fs_alpha2), &
!!$                  abs(fs_alpha3), &
!!$                  abs(fs_alpha4), &
!!$                  abs(fs_alpha5) )
!!$
!!$  fs_alpha4 = max(abs(fs_alpha1), &
!!$                  abs(fs_alpha2), &
!!$                  abs(fs_alpha3), &
!!$                  abs(fs_alpha4), &
!!$                  abs(fs_alpha5) )
!!$
!!$  fs_alpha5 = max(abs(fs_alpha1), &
!!$                  abs(fs_alpha2), &
!!$                  abs(fs_alpha3), &
!!$                  abs(fs_alpha4), &
!!$                  abs(fs_alpha5) )

!!$  write(*,*) 'fs_alpha: final',fs_alpha1,fs_alpha2,fs_alpha3,fs_alpha4,fs_alpha5
  
end subroutine GRHydro_FSAlpha

 /*@@
   @routine    GRHydro_SplitFlux
   @date       Wed Mar  3 22:17:14 2004
   @author     Ian Hawke
   @desc 
   Wrapper routine to split the flux. Synchronization is 
   unfortunately required afterwards.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/


subroutine GRHydro_SplitFlux(CCTK_ARGUMENTS)

  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT, parameter :: izero = 0
  integer, parameter :: ik = kind(izero)

  integer :: nx, ny, nz, i, j, k

  CCTK_REAL, dimension(:), allocatable :: upper, det, dummy
  CCTK_REAL :: uxx, uxy, uxz, uyy, uyz, uzz

  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)

  if (flux_direction == 1) then
    allocate(upper(nx), det(nx), dummy(nx))
    do k = GRHydro_stencil, cctk_lsh(3) - GRHydro_stencil + 2
      do j = GRHydro_stencil, cctk_lsh(2) - GRHydro_stencil + 2

         dummy = betax(:,j,k)

        do i = 1, cctk_lsh(1)
           det(i) = sdetg(i,j,k)*sdetg(i,j,k)
           call UpperMetric(uxx,uxy,uxz,uyy,uyz,uzz,det(i),&
                gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),&
                gyz(i,j,k),gzz(i,j,k))        
           upper(i) = uxx
        end do

        call GRHydro_SplitFlux_1D(GRHydro_eos_handle, int(nx,ik),&
             fs_alpha1, fs_alpha2, fs_alpha3, fs_alpha4, fs_alpha5, &
             gxx(:,j,k),gxy(:,j,k),gxz(:,j,k),&
             gyy(:,j,k),gyz(:,j,k),gzz(:,j,k),&
             upper, det,&
             alp(:,j,k),dummy,&
             rho(:,j,k),vel(:,j,k,1),vel(:,j,k,2),vel(:,j,k,3),press(:,j,k),&
             w_lorentz(:,j,k),eps(:,j,k),&
             dens(:,j,k),scon(:,j,k,1),scon(:,j,k,2),scon(:,j,k,3),tau(:,j,k),&
             densfplus(:,j,k), sxfplus(:,j,k), syfplus(:,j,k), &
             szfplus(:,j,k), taufplus(:,j,k), &
             densfminus(:,j,k), sxfminus(:,j,k), syfminus(:,j,k), &
             szfminus(:,j,k), taufminus(:,j,k), &
             densflux(:,j,k), sxflux(:,j,k), syflux(:,j,k), &
             szflux(:,j,k), tauflux(:,j,k))        
      end do
    end do
    deallocate(upper, det, dummy)
  else if (flux_direction == 2) then
    allocate(upper(ny), det(ny), dummy(ny))
    do k = GRHydro_stencil, cctk_lsh(3) - GRHydro_stencil + 2
      do i = GRHydro_stencil, cctk_lsh(1) - GRHydro_stencil + 2

         dummy = betay(i,:,k)

        do j = 1, cctk_lsh(2)
           det(j) = sdetg(i,j,k)**2
           call UpperMetric(uxx,uxy,uxz,uyy,uyz,uzz,det(j),&
                gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),&
                gyz(i,j,k),gzz(i,j,k))        
           upper(j) = uyy
        end do

        call GRHydro_SplitFlux_1D(GRHydro_eos_handle, int(ny,ik),&
             fs_alpha1, fs_alpha2, fs_alpha3, fs_alpha4, fs_alpha5, &
             gyy(i,:,k),gyz(i,:,k),gxy(i,:,k),&
             gzz(i,:,k),gxz(i,:,k),gxx(i,:,k),&
             upper, det,&
             alp(i,:,k),dummy,&
             rho(i,:,k),vel(i,:,k,2),vel(i,:,k,3),vel(i,:,k,1),press(i,:,k),&
             w_lorentz(i,:,k),eps(i,:,k),&
             dens(i,:,k),scon(i,:,k,2),scon(i,:,k,3),scon(i,:,k,1),tau(i,:,k),&
             densfplus(i,:,k), syfplus(i,:,k), szfplus(i,:,k), &
             sxfplus(i,:,k), taufplus(i,:,k), &
             densfminus(i,:,k), syfminus(i,:,k), szfminus(i,:,k), &
             sxfminus(i,:,k), taufminus(i,:,k), &
             densflux(i,:,k), syflux(i,:,k), szflux(i,:,k), &
             sxflux(i,:,k), tauflux(i,:,k))
      end do
    end do
    deallocate(upper, det, dummy)
  else if (flux_direction == 3) then
    allocate(upper(nz), det(nz), dummy(nz))
    do j = GRHydro_stencil, cctk_lsh(2) - GRHydro_stencil + 2
      do i = GRHydro_stencil, cctk_lsh(1) - GRHydro_stencil + 2

         dummy = betaz(i,j,:)
         
        do k = 1, cctk_lsh(3)
           det(k) = sdetg(i,j,k)*sdetg(i,j,k)
           call UpperMetric(uxx,uxy,uxz,uyy,uyz,uzz,det(k),&
                gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),&
                gyz(i,j,k),gzz(i,j,k))        
           upper(k) = uzz
        end do

        call GRHydro_SplitFlux_1D(GRHydro_eos_handle, int(nz,ik),&
             fs_alpha1, fs_alpha2, fs_alpha3, fs_alpha4, fs_alpha5, &
             gzz(i,j,:),gxz(i,j,:),gyz(i,j,:),&
             gxx(i,j,:),gxy(i,j,:),gyy(i,j,:),&
             upper, det,&
             alp(i,j,:),dummy,&
             rho(i,j,:),vel(i,j,:,3),vel(i,j,:,1),vel(i,j,:,2),press(i,j,:),&
             w_lorentz(i,j,:),eps(i,j,:),&
             dens(i,j,:),scon(i,j,:,3),scon(i,j,:,1),scon(i,j,:,2),tau(i,j,:),&
             densfplus(i,j,:), szfplus(i,j,:), sxfplus(i,j,:), &
             syfplus(i,j,:), taufplus(i,j,:), &
             densfminus(i,j,:), szfminus(i,j,:), sxfminus(i,j,:), &
             syfminus(i,j,:), taufminus(i,j,:), &
             densflux(i,j,:), szflux(i,j,:), sxflux(i,j,:), &
             syflux(i,j,:), tauflux(i,j,:))
      end do
    end do
    deallocate(upper, det, dummy)
  else
    call CCTK_ERROR("Flux direction not x,y,z")
    STOP
  end if

  return

end subroutine GRHydro_SplitFlux

 /*@@
   @routine    GRHydro_SplitFlux_1D
   @date       Wed Mar  3 22:55:36 2004
   @author     Ian Hawke
   @desc 
   Actually performs the flux splitting. 
   We really should be doing it on a characteristic basis.

   Note the direction rotations on the calling routine.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine GRHydro_SplitFlux_1D(handle, nx, &
     alpha1, alpha2, alpha3, alpha4, alpha5, &
     gxx, gxy, gxz, gyy, gyz, gzz, u, det, alp, beta, &
     rho, velx1, vely1, velz1, press, w_lorentz, eps, &
     dens, sx, sy, sz, tau, &
     densfplus, sxfplus, syfplus, szfplus, taufplus, &
     densfminus, sxfminus, syfminus, szfminus, taufminus, &
     densflux, sxflux, syflux, szflux, tauflux)

  use GRHydro_Eigenproblem

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT, parameter :: izero = 0
  integer, parameter :: ik = kind(izero)

  CCTK_REAL, parameter :: half = 0.5d0

  CCTK_INT :: i, nx, handle, ll
  CCTK_REAL, dimension(nx) :: gxx, gxy, gxz, gyy, gyz, gzz, &
       u, det, alp, beta, &
       rho, velx1, vely1, velz1, press, w_lorentz, eps, &
       dens, sx, sy, sz, tau, &
       densfplus, sxfplus, syfplus, szfplus, taufplus, &
       densfminus, sxfminus, syfminus, szfminus, taufminus, &
       densflux, sxflux, syflux, szflux, tauflux

  CCTK_REAL :: alpha1, alpha2, alpha3, alpha4, alpha5

  CCTK_REAL, dimension(:), allocatable :: tmp

  CCTK_REAL, dimension(:,:), allocatable :: evals
  CCTK_REAL, dimension(:,:,:), allocatable :: levecs, revecs

  CCTK_REAL, dimension(:), allocatable :: &
     char_v1, char_v2, char_v3, char_v4 ,char_v5, &
     char_f1_plus, char_f2_plus, char_f3_plus, &
     char_f4_plus ,char_f5_plus, &
     char_f1_minus, char_f2_minus, char_f3_minus, &
     char_f4_minus ,char_f5_minus, &
     char_f1, char_f2, char_f3, &
     char_f4 ,char_f5

  CCTK_REAL, dimension(5) :: lambda, alpha
  CCTK_REAL, dimension(5, 5) :: levec, revec

  densfplus = 0.d0
  sxfplus = 0.d0
  syfplus = 0.d0
  szfplus = 0.d0
  taufplus = 0.d0
  densfminus = 0.d0
  sxfminus = 0.d0
  syfminus = 0.d0
  szfminus = 0.d0
  taufminus = 0.d0

  allocate(evals(nx, 5), levecs(nx,5,5), revecs(nx,5,5))
  allocate(tmp(nx))
  allocate(char_v1(nx), char_v2(nx), char_v3(nx), char_v4 (nx),char_v5(nx), &
     char_f1_plus(nx), char_f2_plus(nx), char_f3_plus(nx), &
     char_f4_plus (nx),char_f5_plus(nx), &
     char_f1_minus(nx), char_f2_minus(nx), char_f3_minus(nx), &
     char_f4_minus (nx),char_f5_minus(nx), &
     char_f1(nx), char_f2(nx), char_f3(nx), &
     char_f4 (nx),char_f5(nx))

  do i = 1, nx-1

!!$    Arithmetic mean

!!$    Calculate the maximum eigenvalue and put it here

    call eigenproblem_leftright(handle, &
         half * (rho      (i) + rho      (i+1)), &
         half * (velx1    (i) + velx1    (i+1)), &
         half * (vely1    (i) + vely1    (i+1)), &
         half * (velz1    (i) + velz1    (i+1)), &
         half * (eps      (i) + eps      (i+1)), &
         half * (w_lorentz(i) + w_lorentz(i+1)), &
         half * (gxx      (i) + gxx      (i+1)), &
         half * (gxy      (i) + gxy      (i+1)), &
         half * (gxz      (i) + gxz      (i+1)), &
         half * (gyy      (i) + gyy      (i+1)), &
         half * (gyz      (i) + gyz      (i+1)), &
         half * (gzz      (i) + gzz      (i+1)), &
         half * (u        (i) + u        (i+1)), &
         half * (alp      (i) + alp      (i+1)), &
         half * (beta     (i) + beta     (i+1)), &
         lambda,&
         levec,&
         revec)

    evals(i,:) = lambda
    levecs(i,:,:) = levec
    revecs(i,:,:) = revec
    
  end do
  
  do i = 3, nx - 3

    alpha(1) = alpha1
    alpha(2) = alpha2
    alpha(3) = alpha3
    alpha(4) = alpha4
    alpha(5) = alpha5
    
    do ll = i - 2, i + 3

!!$    Initialize the pointwise fluxes temporarily into the minus

      call num_x_flux(dens(ll), sx(ll), sy(ll), sz(ll), tau(ll), &
           densfminus(ll), sxfminus(ll), syfminus(ll), szfminus(ll), &
           taufminus(ll), &
           velx1(ll), press(ll), det(ll), alp(ll), beta(ll))

      densfminus(ll) = densfminus(ll) * alp(ll)
      sxfminus(ll)   = sxfminus(ll)   * alp(ll)
      syfminus(ll)   = syfminus(ll)   * alp(ll)
      szfminus(ll)   = szfminus(ll)   * alp(ll)
      taufminus(ll)  = taufminus(ll)  * alp(ll)

      char_v1(ll) = levecs(i,1,1) * dens(ll) + &
           levecs(i,1,2) * sx(ll) + &
           levecs(i,1,3) * sy(ll) + &
           levecs(i,1,4) * sz(ll) + &
           levecs(i,1,5) * tau(ll)
      char_v2(ll) = levecs(i,2,1) * dens(ll) + &
           levecs(i,2,2) * sx(ll) + &
           levecs(i,2,3) * sy(ll) + &
           levecs(i,2,4) * sz(ll) + &
           levecs(i,2,5) * tau(ll)
      char_v3(ll) = levecs(i,3,1) * dens(ll) + &
           levecs(i,3,2) * sx(ll) + &
           levecs(i,3,3) * sy(ll) + &
           levecs(i,3,4) * sz(ll) + &
           levecs(i,3,5) * tau(ll)
      char_v4(ll) = levecs(i,4,1) * dens(ll) + &
           levecs(i,4,2) * sx(ll) + &
           levecs(i,4,3) * sy(ll) + &
           levecs(i,4,4) * sz(ll) + &
           levecs(i,4,5) * tau(ll)
      char_v5(ll) = levecs(i,5,1) * dens(ll) + &
           levecs(i,5,2) * sx(ll) + &
           levecs(i,5,3) * sy(ll) + &
           levecs(i,5,4) * sz(ll) + &
           levecs(i,5,5) * tau(ll)
      
      char_f1(ll) = levecs(i,1,1) * densfminus(ll) + &
           levecs(i,1,2) * sxfminus(ll) + &
           levecs(i,1,3) * syfminus(ll) + &
           levecs(i,1,4) * szfminus(ll) + &
           levecs(i,1,5) * taufminus(ll)
      char_f2(ll) = levecs(i,2,1) * densfminus(ll) + &
           levecs(i,2,2) * sxfminus(ll) + &
           levecs(i,2,3) * syfminus(ll) + &
           levecs(i,2,4) * szfminus(ll) + &
           levecs(i,2,5) * taufminus(ll)
      char_f3(ll) = levecs(i,3,1) * densfminus(ll) + &
           levecs(i,3,2) * sxfminus(ll) + &
           levecs(i,3,3) * syfminus(ll) + &
           levecs(i,3,4) * szfminus(ll) + &
           levecs(i,3,5) * taufminus(ll)
      char_f4(ll) = levecs(i,4,1) * densfminus(ll) + &
           levecs(i,4,2) * sxfminus(ll) + &
           levecs(i,4,3) * syfminus(ll) + &
           levecs(i,4,4) * szfminus(ll) + &
           levecs(i,4,5) * taufminus(ll)
      char_f5(ll) = levecs(i,5,1) * densfminus(ll) + &
           levecs(i,5,2) * sxfminus(ll) + &
           levecs(i,5,3) * syfminus(ll) + &
           levecs(i,5,4) * szfminus(ll) + &
           levecs(i,5,5) * taufminus(ll)      

!!$    Calculate the split

      char_f1_plus(ll)  = 0.5d0 * (char_f1(ll) + alpha(1) * char_v1(ll))
      char_f1_minus(ll) = 0.5d0 * (char_f1(ll) - alpha(1) * char_v1(ll))
      char_f2_plus(ll)  = 0.5d0 * (char_f2(ll) + alpha(2) * char_v2(ll))
      char_f2_minus(ll) = 0.5d0 * (char_f2(ll) - alpha(2) * char_v2(ll))
      char_f3_plus(ll)  = 0.5d0 * (char_f3(ll) + alpha(3) * char_v3(ll))
      char_f3_minus(ll) = 0.5d0 * (char_f3(ll) - alpha(3) * char_v3(ll))
      char_f4_plus(ll)  = 0.5d0 * (char_f4(ll) + alpha(4) * char_v4(ll))
      char_f4_minus(ll) = 0.5d0 * (char_f4(ll) - alpha(4) * char_v4(ll))
      char_f5_plus(ll)  = 0.5d0 * (char_f5(ll) + alpha(5) * char_v5(ll))
      char_f5_minus(ll) = 0.5d0 * (char_f5(ll) - alpha(5) * char_v5(ll))

    end do

!!$    Reconstruct the characteristic split fluxes
!!$    After reconstruction, combine to get the characteristic flux

    call GRHydro_WENO5_Left(5_ik, char_f1_plus(i-2:i+2), tmp(i-2:i+2))
    char_f1(i) = tmp(i)
    call GRHydro_WENO5_Right(5_ik, char_f1_minus(i-1:i+3), &
         tmp(i-1:i+3))
    char_f1(i) = char_f1(i) + tmp(i+1)
    
    call GRHydro_WENO5_Left(5_ik, char_f2_plus(i-2:i+2), tmp(i-2:i+2))
    char_f2(i) = tmp(i)
    call GRHydro_WENO5_Right(5_ik, char_f2_minus(i-1:i+3), &
         tmp(i-1:i+3))
    char_f2(i) = char_f2(i) + tmp(i+1)
    
    call GRHydro_WENO5_Left(5_ik, char_f3_plus(i-2:i+2), tmp(i-2:i+2))
    char_f3(i) = tmp(i)
    call GRHydro_WENO5_Right(5_ik, char_f3_minus(i-1:i+3), &
         tmp(i-1:i+3))
    char_f3(i) = char_f3(i) + tmp(i+1)
    
    call GRHydro_WENO5_Left(5_ik, char_f4_plus(i-2:i+2), tmp(i-2:i+2))
    char_f4(i) = tmp(i)
    call GRHydro_WENO5_Right(5_ik, char_f4_minus(i-1:i+3), &
         tmp(i-1:i+3))
    char_f4(i) = char_f4(i) + tmp(i+1)
    
    call GRHydro_WENO5_Left(5_ik, char_f5_plus(i-2:i+2), tmp(i-2:i+2))
    char_f5(i) = tmp(i)
    call GRHydro_WENO5_Right(5_ik, char_f5_minus(i-1:i+3), &
         tmp(i-1:i+3))
    char_f5(i) = char_f5(i) + tmp(i+1)

!!$    Compute physical fluxes
    
    densflux(i) = &
         revecs(i,1,1) * char_f1(i) + revecs(i,2,1) * char_f2(i) + &
         revecs(i,3,1) * char_f3(i) + revecs(i,4,1) * char_f4(i) + &
         revecs(i,5,1) * char_f5(i)
    sxflux(i) = &
         revecs(i,1,2) * char_f1(i) + revecs(i,2,2) * char_f2(i) + &
         revecs(i,3,2) * char_f3(i) + revecs(i,4,2) * char_f4(i) + &
         revecs(i,5,2) * char_f5(i)
    syflux(i) = &
         revecs(i,1,3) * char_f1(i) + revecs(i,2,3) * char_f2(i) + &
         revecs(i,3,3) * char_f3(i) + revecs(i,4,3) * char_f4(i) + &
         revecs(i,5,3) * char_f5(i)
    szflux(i) = &
         revecs(i,1,4) * char_f1(i) + revecs(i,2,4) * char_f2(i) + &
         revecs(i,3,4) * char_f3(i) + revecs(i,4,4) * char_f4(i) + &
         revecs(i,5,4) * char_f5(i)
    tauflux(i) = &
         revecs(i,1,5) * char_f1(i) + revecs(i,2,5) * char_f2(i) + &
         revecs(i,3,5) * char_f3(i) + revecs(i,4,5) * char_f4(i) + &
         revecs(i,5,5) * char_f5(i)

!!$    if (abs(i-200) < 5) then
!!$      write(*,*) i, 'alpha',alpha1,alpha2,alpha3,alpha4,alpha5
!!$      write(*,*) i, 'var',dens(i), sx(i), tau(i)
!!$      write(*,*) i, 'f',densflux(i),sxflux(i),tauflux(i)
!!$    end if
    
  end do

  deallocate(evals, levecs, revecs, tmp)

end subroutine GRHydro_SplitFlux_1D

 /*@@
   @routine    GRHydro_WENO5_Left
   @date       Wed Mar  3 22:57:54 2004
   @author     Ian Hawke
   @desc 
   Upwind biased WENO5.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine GRHydro_WENO5_Left(nx, v, vminus)

  implicit none

  CCTK_INT :: i, nx
  CCTK_REAL, dimension(nx) :: v, vminus
  CCTK_REAL :: one, two, three, four, five, six, seven, &
       ten, eleven, twelve, thirteen, &
       ThirteenByTwelve, Quarter
  parameter (one = 1)
  parameter (two = 2)
  parameter (three = 3)
  parameter (four = 4)
  parameter (five = 5)
  parameter (six = 6)
  parameter (seven = 7)
  parameter (ten = 10)
  parameter (eleven = 11)
  parameter (twelve = 12)
  parameter (thirteen = 13)
  parameter (ThirteenByTwelve = thirteen / twelve)
  parameter (Quarter = one / four)
  CCTK_REAL :: d0, d1, d2
  parameter (d0 = three / ten)
  parameter (d1 = six / ten)
  parameter (d2 = one / ten)
  CCTK_REAL :: c00,c01,c02,c10,c11,c12,c20,c21,c22
  parameter (c00 = two / six)
  parameter (c01 = five / six)
  parameter (c02 = -one / six)
  parameter (c10 = -one / six)
  parameter (c11 = five / six)
  parameter (c12 = two / six)
  parameter (c20 = two / six)
  parameter (c21 = -seven / six)
  parameter (c22 = eleven / six)

  CCTK_REAL :: beta0, beta1, beta2
  CCTK_REAL :: epsilon
  CCTK_REAL :: alpha0, alpha1, alpha2, alphasum
  CCTK_REAL :: w0, w1, w2
  CCTK_REAL :: v0plushalf, v1plushalf, v2plushalf

  epsilon = 1.d-6

  do i = 3, nx-2

    beta0 = ThirteenByTwelve * (v(i  ) - two * v(i+1) + v(i+2))**2 + &
         Quarter * (three * v(i  ) - four * v(i+1) +         v(i+2))**2
    beta1 = ThirteenByTwelve * (v(i-1) - two * v(i  ) + v(i+1))**2 + &
         Quarter * (        v(i-1)                 -         v(i+1))**2
    beta2 = ThirteenByTwelve * (v(i-2) - two * v(i-1) + v(i  ))**2 + &
         Quarter * (        v(i-2) - four * v(i-1) + three * v(i  ))**2
    
    alpha0 = d0 / (epsilon + beta0)**2
    alpha1 = d1 / (epsilon + beta1)**2
    alpha2 = d2 / (epsilon + beta2)**2
    
    alphasum = alpha0 + alpha1 + alpha2
    
    w0 = alpha0 / alphasum
    w1 = alpha1 / alphasum
    w2 = alpha2 / alphasum
    
    v0plushalf = c00 * v(i  ) + c01 * v(i+1) + c02 * v(i+2)
    v1plushalf = c10 * v(i-1) + c11 * v(i  ) + c12 * v(i+1)
    v2plushalf = c20 * v(i-2) + c21 * v(i-1) + c22 * v(i  )
    
    vminus(i) = w0 * v0plushalf + &
                w1 * v1plushalf + &
                w2 * v2plushalf
  
  end do
  
end subroutine GRHydro_WENO5_Left

 /*@@
   @routine    GRHydro_WENO5_Right
   @date       Wed Mar  3 22:58:10 2004
   @author     Ian Hawke
   @desc 
   Downwind biased WENO5.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine GRHydro_WENO5_Right(nx, v, vplus)

  implicit none

  CCTK_INT :: i, nx
  CCTK_REAL, dimension(nx) :: v, vplus
  CCTK_REAL :: one, two, three, four, five, six, seven, &
       ten, eleven, twelve, thirteen, &
       ThirteenByTwelve, Quarter
  parameter (one = 1)
  parameter (two = 2)
  parameter (three = 3)
  parameter (four = 4)
  parameter (five = 5)
  parameter (six = 6)
  parameter (seven = 7)
  parameter (ten = 10)
  parameter (eleven = 11)
  parameter (twelve = 12)
  parameter (thirteen = 13)
  parameter (ThirteenByTwelve = thirteen / twelve)
  parameter (Quarter = one / four)
  CCTK_REAL :: dtilde0, dtilde1, dtilde2
  parameter (dtilde0 = one / ten)
  parameter (dtilde1 = six / ten)
  parameter (dtilde2 = three / ten)
  CCTK_REAL :: ctilde00,ctilde01,ctilde02,ctilde10,ctilde11,&
       ctilde12,ctilde20,ctilde21,ctilde22
  parameter (ctilde00 = eleven / six)
  parameter (ctilde01 = -seven / six)
  parameter (ctilde02 = two / six)
  parameter (ctilde10 = two / six)
  parameter (ctilde11 = five / six)
  parameter (ctilde12 = -one / six)
  parameter (ctilde20 = -one / six)
  parameter (ctilde21 = five / six)
  parameter (ctilde22 = two / six)

  CCTK_REAL :: betatilde0, betatilde1, betatilde2
  CCTK_REAL :: epsilon
  CCTK_REAL :: alphatilde0, alphatilde1, alphatilde2, alphatildesum
  CCTK_REAL :: wtilde0, wtilde1, wtilde2
  CCTK_REAL :: v0minushalf, v1minushalf, v2minushalf

  epsilon = 1.d-6

  do i = 3, nx-2

    betatilde0 = ThirteenByTwelve * (v(i  ) - two * v(i+1) + v(i+2))**2 + &
         Quarter * (three * v(i  ) - four * v(i+1) +         v(i+2))**2
    betatilde1 = ThirteenByTwelve * (v(i-1) - two * v(i  ) + v(i+1))**2 + &
         Quarter * (        v(i-1)                 -         v(i+1))**2
    betatilde2 = ThirteenByTwelve * (v(i-2) - two * v(i-1) + v(i  ))**2 + &
         Quarter * (        v(i-2) - four * v(i-1) + three * v(i  ))**2
    
    alphatilde0 = dtilde0 / (epsilon + betatilde0)**2
    alphatilde1 = dtilde1 / (epsilon + betatilde1)**2
    alphatilde2 = dtilde2 / (epsilon + betatilde2)**2
    
    alphatildesum = alphatilde0 + alphatilde1 + alphatilde2
    
    wtilde0 = alphatilde0 / alphatildesum
    wtilde1 = alphatilde1 / alphatildesum
    wtilde2 = alphatilde2 / alphatildesum
    
    v0minushalf = ctilde00 * v(i  ) + ctilde01 * v(i+1) + ctilde02 * v(i+2)
    v1minushalf = ctilde10 * v(i-1) + ctilde11 * v(i  ) + ctilde12 * v(i+1)
    v2minushalf = ctilde20 * v(i-2) + ctilde21 * v(i-1) + ctilde22 * v(i  )
    
    vplus(i) = wtilde0 * v0minushalf + &
               wtilde1 * v1minushalf + &
               wtilde2 * v2minushalf
  
  end do
  
end subroutine GRHydro_WENO5_Right
