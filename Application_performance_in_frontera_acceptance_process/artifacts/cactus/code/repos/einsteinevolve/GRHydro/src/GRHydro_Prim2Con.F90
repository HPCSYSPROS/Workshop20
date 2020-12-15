  /*@@
   @file      primitive2conservative
   @date      Thu Jan  11 11:03:32 2002
   @author    Pedro Montero, Ian Hawke
   @desc 
   Primitive to conservative routine
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"
#include "GRHydro_Macros.h"

#define velx(i,j,k) vup(i,j,k,1)
#define vely(i,j,k) vup(i,j,k,2)
#define velz(i,j,k) vup(i,j,k,3)
#define sx(i,j,k) scon(i,j,k,1)
#define sy(i,j,k) scon(i,j,k,2)
#define sz(i,j,k) scon(i,j,k,3)


 /*@@
   @routine   primitive2conservative.f90 
   @date       Thu Jan 11 11:03:32 2002
   @author     Pedro Montero, Ian Hawke
   @desc 
   Converts primitive to conserved variables for the boundary extended data.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine primitive2conservative(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  CCTK_INT, parameter :: izero = 0
  integer, parameter :: ik = kind(izero)

  integer :: i, j, k
  CCTK_INT :: keytemp
  CCTK_REAL :: g11l,g12l,g13l,g22l,g23l,g33l,avg_sdetl,&
       g11r,g12r,g13r,g22r,g23r,g33r,avg_sdetr
  CCTK_REAL :: xtemp(1)

  ! save memory when MP is not used
  CCTK_INT :: GRHydro_UseGeneralCoordinates
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3) :: g11, g12, g13, g22, g23, g33
  pointer (pg11,g11), (pg12,g12), (pg13,g13), (pg22,g22), (pg23,g23), (pg33,g33)
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3,3) :: vup
  pointer (pvup,vup)

  if (GRHydro_UseGeneralCoordinates(cctkGH).ne.0) then
    pg11 = loc(gaa)
    pg12 = loc(gab)
    pg13 = loc(gac)
    pg22 = loc(gbb)
    pg23 = loc(gbc)
    pg33 = loc(gcc)
    pvup = loc(lvel)
  else
    pg11 = loc(gxx)
    pg12 = loc(gxy)
    pg13 = loc(gxz)
    pg22 = loc(gyy)
    pg23 = loc(gyz)
    pg33 = loc(gzz)
    pvup = loc(vel)
  end if

  if(evolve_temper.ne.1) then
     !$OMP PARALLEL DO PRIVATE(i, j, k, avg_sdetl, avg_sdetr,&
     !$OMP                      g11l,g12l,g13l,g22l,g23l,g33l, &
     !$OMP                      g11r,g12r,g13r,g22r,g23r,g33r)
     do k = GRHydro_stencil,cctk_lsh(3)-GRHydro_stencil+1
        do j = GRHydro_stencil,cctk_lsh(2)-GRHydro_stencil+1
           do i = GRHydro_stencil,cctk_lsh(1)-GRHydro_stencil+1
              
              g11l = 0.5d0 * (g11(i,j,k) + g11(i-xoffset,j-yoffset,k-zoffset))
              g12l = 0.5d0 * (g12(i,j,k) + g12(i-xoffset,j-yoffset,k-zoffset))
              g13l = 0.5d0 * (g13(i,j,k) + g13(i-xoffset,j-yoffset,k-zoffset))
              g22l = 0.5d0 * (g22(i,j,k) + g22(i-xoffset,j-yoffset,k-zoffset))
              g23l = 0.5d0 * (g23(i,j,k) + g23(i-xoffset,j-yoffset,k-zoffset))
              g33l = 0.5d0 * (g33(i,j,k) + g33(i-xoffset,j-yoffset,k-zoffset))
              g11r = 0.5d0 * (g11(i,j,k) + g11(i+xoffset,j+yoffset,k+zoffset))
              g12r = 0.5d0 * (g12(i,j,k) + g12(i+xoffset,j+yoffset,k+zoffset))
              g13r = 0.5d0 * (g13(i,j,k) + g13(i+xoffset,j+yoffset,k+zoffset))
              g22r = 0.5d0 * (g22(i,j,k) + g22(i+xoffset,j+yoffset,k+zoffset))
              g23r = 0.5d0 * (g23(i,j,k) + g23(i+xoffset,j+yoffset,k+zoffset))
              g33r = 0.5d0 * (g33(i,j,k) + g33(i+xoffset,j+yoffset,k+zoffset))
              
              avg_sdetl = sqrt(SPATIAL_DETERMINANT(g11l,g12l,g13l,g22l, g23l,g33l))
              avg_sdetr = sqrt(SPATIAL_DETERMINANT(g11r,g12r,g13r,g22r, g23r,g33r))
              
              call prim2con(GRHydro_eos_handle, g11l,g12l,g13l,g22l,& 
                   g23l,g33l, &
                   avg_sdetl,densminus(i,j,k),sxminus(i,j,k),&
                   syminus(i,j,k),szminus(i,j,k),tauminus(i,j,k),&
                   rhominus(i,j,k), &
                   velxminus(i,j,k),velyminus(i,j,k),velzminus(i,j,k),&
                   epsminus(i,j,k),pressminus(i,j,k),w_lorentzminus(i, j, k))
              
              call prim2con(GRHydro_eos_handle, g11r,g12r,g13r,&
                   g22r,g23r,g33r, &
                   avg_sdetr, densplus(i,j,k),sxplus(i,j,k),&
                   syplus(i,j,k),szplus(i,j ,k),tauplus(i,j,k),&
                   rhoplus(i,j,k),velxplus(i,j,k),velyplus(i,j,k),&
                   velzplus(i,j,k),epsplus(i,j,k),pressplus(i,j,k),&
                   w_lorentzplus(i,j,k)) 
              
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  else
     if(reconstruct_temper.ne.0) then
        keytemp = 1
        !$OMP PARALLEL DO PRIVATE(i, j, k, avg_sdetl, avg_sdetr, xtemp,&
        !$OMP                      g11l,g12l,g13l,g22l,g23l,g33l, &
        !$OMP                      g11r,g12r,g13r,g22r,g23r,g33r)
        do k = GRHydro_stencil,cctk_lsh(3)-GRHydro_stencil+1
           do j = GRHydro_stencil,cctk_lsh(2)-GRHydro_stencil+1
              do i = GRHydro_stencil,cctk_lsh(1)-GRHydro_stencil+1
                 
                 g11l = 0.5d0 * (g11(i,j,k) + g11(i-xoffset,j-yoffset,k-zoffset))
                 g12l = 0.5d0 * (g12(i,j,k) + g12(i-xoffset,j-yoffset,k-zoffset))
                 g13l = 0.5d0 * (g13(i,j,k) + g13(i-xoffset,j-yoffset,k-zoffset))
                 g22l = 0.5d0 * (g22(i,j,k) + g22(i-xoffset,j-yoffset,k-zoffset))
                 g23l = 0.5d0 * (g23(i,j,k) + g23(i-xoffset,j-yoffset,k-zoffset))
                 g33l = 0.5d0 * (g33(i,j,k) + g33(i-xoffset,j-yoffset,k-zoffset))
                 g11r = 0.5d0 * (g11(i,j,k) + g11(i+xoffset,j+yoffset,k+zoffset))
                 g12r = 0.5d0 * (g12(i,j,k) + g12(i+xoffset,j+yoffset,k+zoffset))
                 g13r = 0.5d0 * (g13(i,j,k) + g13(i+xoffset,j+yoffset,k+zoffset))
                 g22r = 0.5d0 * (g22(i,j,k) + g22(i+xoffset,j+yoffset,k+zoffset))
                 g23r = 0.5d0 * (g23(i,j,k) + g23(i+xoffset,j+yoffset,k+zoffset))
                 g33r = 0.5d0 * (g33(i,j,k) + g33(i+xoffset,j+yoffset,k+zoffset))
              
                 avg_sdetl = sqrt(SPATIAL_DETERMINANT(g11l,g12l,g13l,g22l, g23l,g33l))
                 avg_sdetr = sqrt(SPATIAL_DETERMINANT(g11r,g12r,g13r,g22r, g23r,g33r))

                 call prim2con_hot(GRHydro_eos_handle, keytemp, GRHydro_reflevel,&
                      int(cctk_iteration,ik),int(i,ik),int(j,ik),int(k,ik),x(i,j,k),y(i,j,k),z(i,j,k),&
                      r(i,j,k),&
                      g11l,g12l,g13l,g22l,& 
                      g23l,g33l, &
                      avg_sdetl,densminus(i,j,k),sxminus(i,j,k),&
                      syminus(i,j,k),szminus(i,j,k),tauminus(i,j,k),&
                      rhominus(i,j,k), &
                      velxminus(i,j,k),velyminus(i,j,k),velzminus(i,j,k),&
                      epsminus(i,j,k),pressminus(i,j,k),w_lorentzminus(i, j, k),&
                      tempminus(i,j,k),y_e_minus(i,j,k))

                 call prim2con_hot(GRHydro_eos_handle, keytemp, GRHydro_reflevel, &
                      int(cctk_iteration,ik),int(i,ik),int(j,ik),int(k,ik),x(i,j,k),y(i,j,k),z(i,j,k),&
                      r(i,j,k),&
                      g11r,g12r,g13r,&
                      g22r,g23r,g33r,&
                      avg_sdetr, densplus(i,j,k),sxplus(i,j,k),&
                      syplus(i,j,k),szplus(i,j ,k),tauplus(i,j,k),&
                      rhoplus(i,j,k),velxplus(i,j,k),velyplus(i,j,k),&
                      velzplus(i,j,k),epsplus(i,j,k),pressplus(i,j,k),&
                      w_lorentzplus(i,j,k),tempplus(i,j,k), &
                      y_e_plus(i,j,k)) 
                 
              end do
           end do
        end do
        !$OMP END PARALLEL DO
     else
        keytemp = 0
        !$OMP PARALLEL DO PRIVATE(i, j, k, avg_sdetl, avg_sdetr, xtemp,&
        !$OMP                      g11l,g12l,g13l,g22l,g23l,g33l, &
        !$OMP                      g11r,g12r,g13r,g22r,g23r,g33r)
        do k = GRHydro_stencil,cctk_lsh(3)-GRHydro_stencil+1
           do j = GRHydro_stencil,cctk_lsh(2)-GRHydro_stencil+1
              do i = GRHydro_stencil,cctk_lsh(1)-GRHydro_stencil+1
                 
                 g11l = 0.5d0 * (g11(i,j,k) + g11(i-xoffset,j-yoffset,k-zoffset))
                 g12l = 0.5d0 * (g12(i,j,k) + g12(i-xoffset,j-yoffset,k-zoffset))
                 g13l = 0.5d0 * (g13(i,j,k) + g13(i-xoffset,j-yoffset,k-zoffset))
                 g22l = 0.5d0 * (g22(i,j,k) + g22(i-xoffset,j-yoffset,k-zoffset))
                 g23l = 0.5d0 * (g23(i,j,k) + g23(i-xoffset,j-yoffset,k-zoffset))
                 g33l = 0.5d0 * (g33(i,j,k) + g33(i-xoffset,j-yoffset,k-zoffset))
                 g11r = 0.5d0 * (g11(i,j,k) + g11(i+xoffset,j+yoffset,k+zoffset))
                 g12r = 0.5d0 * (g12(i,j,k) + g12(i+xoffset,j+yoffset,k+zoffset))
                 g13r = 0.5d0 * (g13(i,j,k) + g13(i+xoffset,j+yoffset,k+zoffset))
                 g22r = 0.5d0 * (g22(i,j,k) + g22(i+xoffset,j+yoffset,k+zoffset))
                 g23r = 0.5d0 * (g23(i,j,k) + g23(i+xoffset,j+yoffset,k+zoffset))
                 g33r = 0.5d0 * (g33(i,j,k) + g33(i+xoffset,j+yoffset,k+zoffset))
              
                 avg_sdetl = sqrt(SPATIAL_DETERMINANT(g11l,g12l,g13l,g22l, g23l,g33l))
                 avg_sdetr = sqrt(SPATIAL_DETERMINANT(g11r,g12r,g13r,g22r, g23r,g33r))

                 xtemp(1) = 0.5d0*(temperature(i,j,k) + &
                      temperature(i-xoffset,j-yoffset,k-zoffset))
                 call prim2con_hot(GRHydro_eos_handle, keytemp, GRHydro_reflevel,&
                      int(cctk_iteration,ik),int(i,ik),int(j,ik),int(k,ik),x(i,j,k),y(i,j,k),z(i,j,k),&
                      r(i,j,k),&
                      g11l,g12l,g13l,g22l,& 
                      g23l,g33l, &
                      avg_sdetl,densminus(i,j,k),sxminus(i,j,k),&
                      syminus(i,j,k),szminus(i,j,k),tauminus(i,j,k),&
                      rhominus(i,j,k), &
                      velxminus(i,j,k),velyminus(i,j,k),velzminus(i,j,k),&
                      epsminus(i,j,k),pressminus(i,j,k),w_lorentzminus(i, j, k),&
                      xtemp,y_e_minus(i,j,k))

                 xtemp(1) = 0.5d0*(temperature(i,j,k) + &
                      temperature(i-xoffset,j-yoffset,k-zoffset))
                 call prim2con_hot(GRHydro_eos_handle, keytemp, GRHydro_reflevel, &
                      int(cctk_iteration,ik),int(i,ik),int(j,ik),int(k,ik),x(i,j,k),y(i,j,k),z(i,j,k),&
                      r(i,j,k),&
                      g11r,g12r,g13r,&
                      g22r,g23r,g33r,&
                      avg_sdetr, densplus(i,j,k),sxplus(i,j,k),&
                      syplus(i,j,k),szplus(i,j ,k),tauplus(i,j,k),&
                      rhoplus(i,j,k),velxplus(i,j,k),velyplus(i,j,k),&
                      velzplus(i,j,k),epsplus(i,j,k),pressplus(i,j,k),&
                      w_lorentzplus(i,j,k),xtemp, &
                      y_e_plus(i,j,k)) 
                 tempminus(i,j,k) = xtemp(1)
              
              end do
           end do
        end do
        !$OMP END PARALLEL DO
     endif
  endif

end subroutine primitive2conservative

 /*@@
   @routine    prim2con
   @date       Sat Jan 26 01:52:18 2002
   @author     Pedro Montero, Ian Hawke
   @desc 
   Converts from primitive to conservative at a single point
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine prim2con(handle, gxx, gxy, gxz, gyy, gyz, gzz, sdetg, ddens, &
     dsx, dsy, dsz, dtau , drho, dvelx, dvely, dvelz, deps, dpress, w) 
  
  implicit none

  DECLARE_CCTK_PARAMETERS
  !!DECLARE_CCTK_FUNCTIONS

  CCTK_REAL :: gxx, gxy, gxz, gyy, gyz, gzz, sdetg
  CCTK_REAL :: ddens, dsx, dsy, dsz, dtau, drho, dvelx, dvely, dvelz,&
       deps, dpress, w, vlowx, vlowy, vlowz   
  CCTK_INT :: handle

! begin EOS Omni vars
  CCTK_INT  :: n,keytemp,anyerr,keyerr(1) 
  CCTK_REAL :: xye,xtemp
  n = 1;keytemp = 0;anyerr = 0;keyerr(1) = 0
  xtemp = 0.0d0; xye = 0.0d0
! end EOS Omni vars

  w = 1.d0 / sqrt(1.d0 - (gxx*dvelx*dvelx + gyy*dvely*dvely + gzz &
       *dvelz*dvelz + 2*gxy*dvelx*dvely + 2*gxz*dvelx *dvelz + 2*gyz&
       *dvely*dvelz))  

  call EOS_Omni_press(handle,keytemp,GRHydro_eos_rf_prec,n,&
       drho,deps,xtemp,xye,dpress,keyerr,anyerr)  

  vlowx = gxx*dvelx + gxy*dvely + gxz*dvelz
  vlowy = gxy*dvelx + gyy*dvely + gyz*dvelz
  vlowz = gxz*dvelx + gyz*dvely + gzz*dvelz

  ddens = sdetg * drho * w 
  dsx = sdetg * (drho*(1+deps)+dpress)*w*w * vlowx
  dsy = sdetg * (drho*(1+deps)+dpress)*w*w * vlowy
  dsz = sdetg * (drho*(1+deps)+dpress)*w*w * vlowz
  dtau = sdetg * ((drho*(1+deps)+dpress)*w*w - dpress) - ddens 

end subroutine prim2con


subroutine prim2con_hot(handle, keytemp, GRHydro_reflevel, cctk_iteration, ii, jj, kk, &
     x, y, z, r, gxx, gxy, gxz, gyy, gyz, gzz, sdetg, ddens, &
     dsx, dsy, dsz, dtau , drho, dvelx, dvely, dvelz, deps, dpress, w, &
     temp,ye) 
  
  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_REAL :: gxx, gxy, gxz, gyy, gyz, gzz, sdetg
  CCTK_REAL :: ddens, dsx, dsy, dsz, dtau, drho(1), dvelx, dvely, dvelz,&
       deps(1), dpress(1), w, vlowx, vlowy, vlowz   
  CCTK_REAL :: temp(1),ye(1), x, y, z, r
  CCTK_INT :: handle, GRHydro_reflevel, cctk_iteration, ii, jj, kk
  CCTK_REAL :: h
  character(len=512) warnline

! begin EOS Omni vars
  CCTK_INT  :: n,keytemp,lkeytemp,anyerr,keyerr(1)
  CCTK_REAL :: temp0(1)
  n = 1;lkeytemp=keytemp;anyerr = 0;keyerr(1) = 0
! end EOS Omni vars

  temp0 = temp
  w = 1.d0 / sqrt(1.d0 - (gxx*dvelx*dvelx + gyy*dvely*dvely + gzz &
       *dvelz*dvelz + 2*gxy*dvelx*dvely + 2*gxz*dvelx *dvelz + 2*gyz&
       *dvely*dvelz))  

  ! make sure Y_e and temperature within allowed range
  ye = max(min(ye,GRHydro_Y_e_max),GRHydro_Y_e_min)
  temp = max(GRHydro_hot_atmo_temp,min(temp,GRHydro_max_temp))

  call EOS_Omni_press(handle,lkeytemp,GRHydro_eos_rf_prec,n,&
       drho,deps,temp,ye,dpress,keyerr,anyerr)  
  ! error handling
  if(anyerr.ne.0) then
     if(reconstruct_temper.ne.0) then
        if(keytemp.eq.1) then
           !$OMP CRITICAL
           call CCTK_WARN(1,"EOS error in prim2con_hot:")
           write(warnline,"(i8,4i5,1P10E15.6)") cctk_iteration,GRHydro_Reflevel,ii,jj,kk,x,y,z,r
           call CCTK_WARN(1,warnline)
           write(warnline,"(1P10E15.6)") drho,deps,temp,ye
           call CCTK_WARN(1,warnline)
           write(warnline,"(A7,i8,A10,i8)") "keyerr: ",keyerr(1)," keytemp:",lkeytemp
           call CCTK_ERROR(warnline)
           STOP
           !$OMP END CRITICAL
        else
           if(keyerr(1).eq.668 .and. temp(1).lt.10.0d0) then
              !$OMP CRITICAL
              write(warnline,"(A18,i7,4i4,1P10E15.6)") "p2c resetting T: ",cctk_iteration, &
                   GRHydro_reflevel,ii,jj,kk,drho,deps,temp,ye
              call CCTK_WARN(1,warnline)
              !$OMP END CRITICAL
              temp = GRHydro_hot_atmo_temp
              lkeytemp = 1
              call EOS_Omni_press(handle,lkeytemp,GRHydro_eos_rf_prec,n,&
                   drho,deps,temp,ye,dpress,keyerr,anyerr)  
              lkeytemp = 0
              if(anyerr.ne.0) then
                 !$OMP CRITICAL
                 write(warnline,"(A7,i8,A10,i8)") "keyerr: ",keyerr(1)," keytemp:",lkeytemp
                 call CCTK_WARN(1,warnline)
                 call CCTK_ERROR("Fatal EOS error in p2c!")
                 STOP
                 !$OMP END CRITICAL
              endif
           else
              !$OMP CRITICAL
              call CCTK_WARN(1,"EOS error in prim2con_hot:")
              write(warnline,"(i8,4i5,1P10E15.6)") cctk_iteration,GRHydro_Reflevel,ii,jj,kk,x,y,z,r
              call CCTK_WARN(1,warnline)
              write(warnline,"(1P10E15.6)") drho,deps,temp,ye
              call CCTK_WARN(1,warnline)
              write(warnline,"(A7,i8,A10,i8)") "keyerr: ",keyerr(1)," keytemp:",lkeytemp
              call CCTK_ERROR(warnline)
              STOP
              !$OMP END CRITICAL
           endif
        endif
     else
        if(GRHydro_reflevel.lt.GRHydro_c2p_warn_from_reflevel) then
           ! in this case (coarse grid error that is hopefully restricted
           ! away), we use the average temperature between cells and call
           ! the EOS with keytemp=1
           lkeytemp=1
           temp = temp0
           call EOS_Omni_press(handle,lkeytemp,GRHydro_eos_rf_prec,n,&
                drho,deps,temp,ye,dpress,keyerr,anyerr)  
           lkeytemp=0
           if(anyerr.ne.0) then
              !$OMP CRITICAL
              call CCTK_WARN(1,"EOS error in prim2con_hot:")
              write(warnline,"(i8,3i5,1P10E15.6)") cctk_iteration,ii,jj,kk,x,y,z,r
              call CCTK_WARN(1,warnline)
              write(warnline,"(1P10E15.6)") drho,deps,temp,ye
              call CCTK_WARN(1,warnline)
              write(warnline,"(A7,i8)") "code: ",keyerr(1)
              call CCTK_WARN(1,warnline)
              write(warnline,"(A10,i5)") "reflevel: ", GRHydro_reflevel
              call CCTK_WARN(1,warnline)
              !$OMP END CRITICAL
           endif
        else
           ! This is a way of recovering even on finer refinement levels:
           ! Use the average temperature at the interface instead of the
           ! reconstructed specific internal energy.
           if(GRHydro_eos_hot_prim2con_warn.ne.0) then
              !$OMP CRITICAL
              call CCTK_WARN(1,"EOS error in prim2con_hot: NOW using averaged temp!")
              write(warnline,"(i8,3i5,1P10E15.6)") cctk_iteration,ii,jj,kk,x,y,z,r
              call CCTK_WARN(1,warnline)
              write(warnline,"(1P10E15.6)") drho,deps,temp0,ye
              call CCTK_WARN(1,warnline)
              write(warnline,"(A7,i8)") "code: ",keyerr(1)
              call CCTK_WARN(1,warnline)
              write(warnline,"(A10,i5)") "reflevel: ", GRHydro_reflevel
              call CCTK_WARN(1,warnline)
              !$OMP END CRITICAL
           endif
           lkeytemp=1
           temp = temp0
           call EOS_Omni_press(handle,lkeytemp,GRHydro_eos_rf_prec,n,&
                drho,deps,temp,ye,dpress,keyerr,anyerr)  
           lkeytemp=0
           if(anyerr.ne.0) then
              !$OMP CRITICAL
              call CCTK_WARN(1,"EOS error in prim2con_hot")
              write(warnline,"(i8,3i5,1P10E15.6)") cctk_iteration,ii,jj,kk,x,y,z,r
              call CCTK_WARN(1,warnline)
              write(warnline,"(1P10E15.6)") drho,deps,temp,ye
              call CCTK_WARN(1,warnline)
              write(warnline,"(A7,i8)") "code: ",keyerr(1)
              call CCTK_WARN(1,warnline)
              write(warnline,"(A10,i5)") "reflevel: ", GRHydro_reflevel
              call CCTK_WARN(1,warnline)
              call CCTK_ERROR("Aborting!!!")
              STOP
              !$OMP END CRITICAL
           endif
        endif
     endif
  endif
  
  vlowx = gxx*dvelx + gxy*dvely + gxz*dvelz
  vlowy = gxy*dvelx + gyy*dvely + gyz*dvelz
  vlowz = gxz*dvelx + gyz*dvely + gzz*dvelz

  h = drho(1)*(1.0d0+deps(1))+dpress(1)
  ddens = sdetg * drho(1) * w 
  dsx = sdetg * h*w*w * vlowx
  dsy = sdetg * h*w*w * vlowy
  dsz = sdetg * h*w*w * vlowz
  dtau = sdetg * (h*w*w - dpress(1)) - ddens 


end subroutine prim2con_hot


 /*@@
   @routine    Primitive2ConservativeCells
   @date       Sun Mar 10 21:16:20 2002
   @author     
   @desc 
   Wrapper function that converts primitive to conservative at the 
     cell centres. 
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/


subroutine Primitive2ConservativeCells(CCTK_ARGUMENTS)

  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  CCTK_INT, parameter :: izero = 0
  integer, parameter :: ik = kind(izero)

  CCTK_INT :: i, j, k
  CCTK_INT :: keytemp

  character(len=512) :: warnline 
  ! save memory when MP is not used
  CCTK_INT :: GRHydro_UseGeneralCoordinates
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3) :: g11, g12, g13, g22, g23, g33
  pointer (pg11,g11), (pg12,g12), (pg13,g13), (pg22,g22), (pg23,g23), (pg33,g33)
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3,3) :: vup
  pointer (pvup,vup)

  if (GRHydro_UseGeneralCoordinates(cctkGH).ne.0) then
    pg11 = loc(gaa)
    pg12 = loc(gab)
    pg13 = loc(gac)
    pg22 = loc(gbb)
    pg23 = loc(gbc)
    pg33 = loc(gcc)
    pvup = loc(lvel)
  else
    pg11 = loc(gxx)
    pg12 = loc(gxy)
    pg13 = loc(gxz)
    pg22 = loc(gyy)
    pg23 = loc(gyz)
    pg33 = loc(gzz)
    pvup = loc(vel)
  end if

  if(evolve_temper.ne.1) then
     !$OMP PARALLEL DO PRIVATE(i, j, k)
     do k = 1,cctk_lsh(3)
        do j = 1,cctk_lsh(2)
           do i = 1,cctk_lsh(1)
          
              call prim2con(GRHydro_eos_handle,g11(i,j,k),&
                   g12(i,j,k),g13(i,j,k),&
                   g22(i,j,k),g23(i,j,k),g33(i,j,k),&
                   sdetg(i,j,k), dens(i,j,k),sx(i,j,k),sy(i,j,k),sz(i,j,k),&
                   tau(i,j,k),rho(i,j,k),velx(i,j,k),vely(i,j,k),velz(i,j,k),&
                   eps(i,j,k),press(i,j,k),w_lorentz(i,j,k))

           end do
        end do
     end do
     !$OMP END PARALLEL DO
  else
     keytemp = 0
     !$OMP PARALLEL DO PRIVATE(i, j, k)
     do k = 1,cctk_lsh(3)
        do j = 1,cctk_lsh(2)
           do i = 1,cctk_lsh(1)
          
              call prim2con_hot(GRHydro_eos_handle,keytemp,&
                   GRHydro_reflevel,int(cctk_iteration,ik),&
                   int(i,ik),int(j,ik),int(k,ik),&
                   x(i,j,k),y(i,j,k),z(i,j,k),r(i,j,k),g11(i,j,k),&
                   g12(i,j,k),g13(i,j,k),&
                   g22(i,j,k),g23(i,j,k),g33(i,j,k),&
                   sdetg(i,j,k), dens(i,j,k),sx(i,j,k),sy(i,j,k),sz(i,j,k),&
                   tau(i,j,k),rho(i,j,k),velx(i,j,k),vely(i,j,k),velz(i,j,k),&
                   eps(i,j,k),press(i,j,k),w_lorentz(i,j,k), &
                   temperature(i,j,k),y_e(i,j,k))
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  endif

  if(evolve_Y_e.ne.0) then
     !$OMP PARALLEL DO PRIVATE(i, j, k)
     do k = 1,cctk_lsh(3)
        do j = 1,cctk_lsh(2)
           do i = 1,cctk_lsh(1)
              Y_e_con(i,j,k) = Y_e(i,j,k) * dens(i,j,k)
           enddo
        enddo
     enddo
     !$OMP END PARALLEL DO
  endif


end subroutine Primitive2ConservativeCells


 /*@@
   @routine    Prim2ConservativePolytype
   @date       Tue Mar 19 22:52:21 2002
   @author     Ian Hawke
   @desc 
   Same as first routine, only for polytropes.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/


subroutine Prim2ConservativePolytype(CCTK_ARGUMENTS)

  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  integer :: i, j, k
  CCTK_REAL :: g11l,g12l,g13l,g22l,g23l,g33l,avg_sdetl,&
       g11r,g12r,g13r,g22r,g23r,g33r,avg_sdetr

  ! save memory when MP is not used
  CCTK_INT :: GRHydro_UseGeneralCoordinates
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3) :: g11, g12, g13, g22, g23, g33
  pointer (pg11,g11), (pg12,g12), (pg13,g13), (pg22,g22), (pg23,g23), (pg33,g33)
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3,3) :: vup
  pointer (pvup,vup)

  if (GRHydro_UseGeneralCoordinates(cctkGH).ne.0) then
    pg11 = loc(gaa)
    pg12 = loc(gab)
    pg13 = loc(gac)
    pg22 = loc(gbb)
    pg23 = loc(gbc)
    pg33 = loc(gcc)
    pvup = loc(lvel)
  else
    pg11 = loc(gxx)
    pg12 = loc(gxy)
    pg13 = loc(gxz)
    pg22 = loc(gyy)
    pg23 = loc(gyz)
    pg33 = loc(gzz)
    pvup = loc(vel)
  end if
  
  !$OMP PARALLEL DO PRIVATE(i, j, k, g11l,g12l,g13l,g22l,g23l,g33l,avg_sdetl,&
  !$OMP g11r,g12r,g13r,g22r,g23r,g33r,avg_sdetr)
  do k = GRHydro_stencil,cctk_lsh(3)-GRHydro_stencil+1
    do j = GRHydro_stencil,cctk_lsh(2)-GRHydro_stencil+1
      do i = GRHydro_stencil,cctk_lsh(1)-GRHydro_stencil+1
        
        g11l = 0.5d0 * (g11(i,j,k) + g11(i-xoffset,j-yoffset,k-zoffset))
        g12l = 0.5d0 * (g12(i,j,k) + g12(i-xoffset,j-yoffset,k-zoffset))
        g13l = 0.5d0 * (g13(i,j,k) + g13(i-xoffset,j-yoffset,k-zoffset))
        g22l = 0.5d0 * (g22(i,j,k) + g22(i-xoffset,j-yoffset,k-zoffset))
        g23l = 0.5d0 * (g23(i,j,k) + g23(i-xoffset,j-yoffset,k-zoffset))
        g33l = 0.5d0 * (g33(i,j,k) + g33(i-xoffset,j-yoffset,k-zoffset))
        g11r = 0.5d0 * (g11(i,j,k) + g11(i+xoffset,j+yoffset,k+zoffset))
        g12r = 0.5d0 * (g12(i,j,k) + g12(i+xoffset,j+yoffset,k+zoffset))
        g13r = 0.5d0 * (g13(i,j,k) + g13(i+xoffset,j+yoffset,k+zoffset))
        g22r = 0.5d0 * (g22(i,j,k) + g22(i+xoffset,j+yoffset,k+zoffset))
        g23r = 0.5d0 * (g23(i,j,k) + g23(i+xoffset,j+yoffset,k+zoffset))
        g33r = 0.5d0 * (g33(i,j,k) + g33(i+xoffset,j+yoffset,k+zoffset))

        avg_sdetl = sqrt(SPATIAL_DETERMINANT(g11l,g12l,g13l,g22l, g23l,g33l))
        avg_sdetr = sqrt(SPATIAL_DETERMINANT(g11r,g12r,g13r,g22r, g23r,g33r))

        call prim2conpolytype(GRHydro_eos_handle, g11l,g12l,g13l,&
             g22l,g23l,g33l, &
             avg_sdetl,densminus(i,j,k),sxminus(i,j,k),&
             syminus(i,j,k),szminus(i,j,k),tauminus(i,j,k),rhominus(i,j,k), &
             velxminus(i,j,k),velyminus(i,j,k),velzminus(i,j,k),&
             epsminus(i,j,k),pressminus(i,j,k),w_lorentzminus(i, j, k))
 
        call prim2conpolytype(GRHydro_eos_handle, g11r,g12r,g13r,&
             g22r,g23r,g33r, &
             avg_sdetr, densplus(i,j,k),sxplus(i,j,k),&
             syplus(i,j,k),szplus(i,j ,k),tauplus(i,j,k),&
             rhoplus(i,j,k),velxplus(i,j,k),velyplus(i,j,k),&
             velzplus(i,j,k),epsplus(i,j,k),pressplus(i,j,k),&
             w_lorentzplus(i,j,k)) 
      end do
    end do
  end do
  !$OMP END PARALLEL DO


  ! Note on Y_e: We use Y_e_plus and Y_e_minus directly
  ! in the Riemann solver. That's why it is not necessary
  ! to do a prim2con for Y_e

end subroutine Prim2ConservativePolytype

 /*@@
   @routine    prim2conpolytype
   @date       Sat Jan 26 01:52:18 2002
   @author     Pedro Montero, Ian Hawke
   @desc 
   Converts from primitive to conservative at a single point
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine prim2conpolytype(handle, gxx, gxy, gxz, gyy, gyz, &
     gzz, sdetg, ddens, &
     dsx, dsy, dsz, dtau , drho, dvelx, dvely, dvelz, deps, dpress, w) 
  
  implicit none
  
  DECLARE_CCTK_PARAMETERS
  !!DECLARE_CCTK_FUNCTIONS

  CCTK_REAL :: gxx, gxy, gxz, gyy, gyz, gzz, sdetg
  CCTK_REAL :: ddens, dsx, dsy, dsz, dtau, drho, dvelx, dvely, dvelz,&
       deps, dpress, w_tmp, w, vlowx, vlowy, vlowz, sqrtdet   
  CCTK_INT :: handle
  character(len=256) NaN_WarnLine

! begin EOS Omni vars
  CCTK_INT  :: n, keytemp, anyerr, keyerr(1)
  CCTK_REAL :: xpress,xeps,xtemp,xye
  n = 1; keytemp = 0; anyerr = 0; keyerr(1) = 0
  xpress = 0.0d0; xeps = 0.0d0; xtemp = 0.0d0; xye = 0.0d0
! end EOS Omni vars
  
  w_tmp = gxx*dvelx*dvelx + gyy*dvely*dvely + gzz *dvelz*dvelz + &
          2*gxy*dvelx*dvely + 2*gxz*dvelx*dvelz + 2*gyz*dvely*dvelz
  if (w_tmp .ge. 1.d0) then
    ! In theory this should not happen, and even when accepting the fact
    ! that numerically it can, one might be tempted to set w to some large
    ! value in that case. However, this would lead to completely bogus
    ! and hard to trace wrong values below. There is no good value to
    ! choose in this case, but something small is probably the best of
    ! all bad choices.
    !$OMP CRITICAL
    write(NaN_WarnLine,'(a80,2g15.6)') 'Infinite Lorentz factor reset. rho, w_tmp: ', drho, w_tmp
    call CCTK_WARN(GRHydro_NaN_verbose, NaN_WarnLine)
    !$OMP END CRITICAL
    w = 1.d-20
  else
    w = 1.d0 / sqrt(1.d0 - w_tmp)
  endif

  call EOS_Omni_press(handle,keytemp,GRHydro_eos_rf_prec,n,&
       drho,xeps,xtemp,xye,dpress,keyerr,anyerr)

  call EOS_Omni_EpsFromPress(handle,keytemp,GRHydro_eos_rf_prec,n,&
       drho,xeps,xtemp,xye,dpress,deps,keyerr,anyerr)

  vlowx = gxx*dvelx + gxy*dvely + gxz*dvelz
  vlowy = gxy*dvelx + gyy*dvely + gyz*dvelz
  vlowz = gxz*dvelx + gyz*dvely + gzz*dvelz

  ddens = sdetg * drho * w 
  dsx = sdetg * (drho*(1+deps)+dpress)*w*w * vlowx
  dsy = sdetg * (drho*(1+deps)+dpress)*w*w * vlowy
  dsz = sdetg * (drho*(1+deps)+dpress)*w*w * vlowz
  dtau = sdetg * ((drho*(1+deps)+dpress)*w*w - dpress) - ddens 

end subroutine prim2conpolytype


 /*@@
   @routine    Primitive2ConservativePolyCells
   @date       Sun Mar 10 21:16:20 2002
   @author     
   @desc 
   Wrapper function that converts primitive to conservative at the 
     cell centres. 
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/


subroutine Primitive2ConservativePolyCells(CCTK_ARGUMENTS)

  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  CCTK_INT :: i, j, k

  ! save memory when MP is not used
  CCTK_INT :: GRHydro_UseGeneralCoordinates
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3) :: g11, g12, g13, g22, g23, g33
  pointer (pg11,g11), (pg12,g12), (pg13,g13), (pg22,g22), (pg23,g23), (pg33,g33)
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3,3) :: vup
  pointer (pvup,vup)

  if (GRHydro_UseGeneralCoordinates(cctkGH).ne.0) then
    pg11 = loc(gaa)
    pg12 = loc(gab)
    pg13 = loc(gac)
    pg22 = loc(gbb)
    pg23 = loc(gbc)
    pg33 = loc(gcc)
    pvup = loc(lvel)
  else
    pg11 = loc(gxx)
    pg12 = loc(gxy)
    pg13 = loc(gxz)
    pg22 = loc(gyy)
    pg23 = loc(gyz)
    pg33 = loc(gzz)
    pvup = loc(vel)
  end if
  
  !$OMP PARALLEL DO PRIVATE(i, j, k)
     do k = 1,cctk_lsh(3)
        do j = 1,cctk_lsh(2)
           do i = 1,cctk_lsh(1)
        
        call prim2conpolytype(GRHydro_eos_handle,g11(i,j,k),g12(i,j,k),&
             g13(i,j,k),g22(i,j,k),g23(i,j,k),g33(i,j,k),&
             sdetg(i,j,k), dens(i,j,k),sx(i,j,k),sy(i,j,k),sz(i,j,k),&
             tau(i,j,k),rho(i,j,k),velx(i,j,k),vely(i,j,k),velz(i,j,k),&
             eps(i,j,k),press(i,j,k),w_lorentz(i,j,k))

      end do
    end do
  end do
  !$OMP END PARALLEL DO

  if(evolve_Y_e.ne.0) then
     !$OMP PARALLEL DO PRIVATE(i, j, k)
     do k = 1,cctk_lsh(3)
        do j = 1,cctk_lsh(2)
           do i = 1,cctk_lsh(1)
              Y_e_con(i,j,k) = Y_e(i,j,k) * dens(i,j,k)
           enddo
        enddo
     enddo
     !$OMP END PARALLEL DO
  endif

end subroutine Primitive2ConservativePolyCells

 /*@@
   @routine    Prim2ConservativeTracer
   @date       Mon Mar  8 13:32:32 2004
   @author     Ian Hawke
   @desc 
   Gets the conserved tracer variable from the primitive.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine Prim2ConservativeTracer(CCTK_ARGUMENTS)

  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  integer :: i, j, k
  CCTK_REAL :: g11l,g12l,g13l,g22l,g23l,g33l,avg_sdetl,&
       g11r,g12r,g13r,g22r,g23r,g33r,avg_sdetr

  ! save memory when MP is not used
  CCTK_INT :: GRHydro_UseGeneralCoordinates
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3) :: g11, g12, g13, g22, g23, g33
  pointer (pg11,g11), (pg12,g12), (pg13,g13), (pg22,g22), (pg23,g23), (pg33,g33)

  if (GRHydro_UseGeneralCoordinates(cctkGH).ne.0) then
    pg11 = loc(gaa)
    pg12 = loc(gab)
    pg13 = loc(gac)
    pg22 = loc(gbb)
    pg23 = loc(gbc)
    pg33 = loc(gcc)
  else
    pg11 = loc(gxx)
    pg12 = loc(gxy)
    pg13 = loc(gxz)
    pg22 = loc(gyy)
    pg23 = loc(gyz)
    pg33 = loc(gzz)
  end if
  
  do k = GRHydro_stencil,cctk_lsh(3)-GRHydro_stencil+1
    do j = GRHydro_stencil,cctk_lsh(2)-GRHydro_stencil+1
      do i = GRHydro_stencil,cctk_lsh(1)-GRHydro_stencil+1
        
        g11l = 0.5d0 * (g11(i,j,k) + g11(i-xoffset,j-yoffset,k-zoffset))
        g12l = 0.5d0 * (g12(i,j,k) + g12(i-xoffset,j-yoffset,k-zoffset))
        g13l = 0.5d0 * (g13(i,j,k) + g13(i-xoffset,j-yoffset,k-zoffset))
        g22l = 0.5d0 * (g22(i,j,k) + g22(i-xoffset,j-yoffset,k-zoffset))
        g23l = 0.5d0 * (g23(i,j,k) + g23(i-xoffset,j-yoffset,k-zoffset))
        g33l = 0.5d0 * (g33(i,j,k) + g33(i-xoffset,j-yoffset,k-zoffset))
        g11r = 0.5d0 * (g11(i,j,k) + g11(i+xoffset,j+yoffset,k+zoffset))
        g12r = 0.5d0 * (g12(i,j,k) + g12(i+xoffset,j+yoffset,k+zoffset))
        g13r = 0.5d0 * (g13(i,j,k) + g13(i+xoffset,j+yoffset,k+zoffset))
        g22r = 0.5d0 * (g22(i,j,k) + g22(i+xoffset,j+yoffset,k+zoffset))
        g23r = 0.5d0 * (g23(i,j,k) + g23(i+xoffset,j+yoffset,k+zoffset))
        g33r = 0.5d0 * (g33(i,j,k) + g33(i+xoffset,j+yoffset,k+zoffset))

        avg_sdetl = sqrt(SPATIAL_DETERMINANT(g11l,g12l,g13l,g22l, g23l,g33l))
        avg_sdetr = sqrt(SPATIAL_DETERMINANT(g11r,g12r,g13r,g22r, g23r,g33r))
        cons_tracerplus(i,j,k,:) = tracerplus(i,j,k,:) * &
             avg_sdetr * rhoplus(i,j,k) / &
             sqrt(1.d0 - &
                   (g11r * velxplus(i,j,k)**2 + &
                    g22r * velyplus(i,j,k)**2 + &
                    g33r * velzplus(i,j,k)**2 + &
                    2.d0 * (g12r * velxplus(i,j,k) * velyplus(i,j,k) + &
                            g13r * velxplus(i,j,k) * velzplus(i,j,k) + &
                            g23r * velyplus(i,j,k) * velzplus(i,j,k) ) ) )
        cons_tracerminus(i,j,k,:) = tracerminus(i,j,k,:) * &
             avg_sdetl * rhominus(i,j,k) / &
             sqrt(1.d0 - &
                   (g11l * velxminus(i,j,k)**2 + &
                    g22l * velyminus(i,j,k)**2 + &
                    g33l * velzminus(i,j,k)**2 + &
                    2.d0 * (g12l * velxminus(i,j,k) * velyminus(i,j,k) + &
                            g13l * velxminus(i,j,k) * velzminus(i,j,k) + &
                            g23l * velyminus(i,j,k) * velzminus(i,j,k) ) ) )

      end do
    end do
  end do

end subroutine Prim2ConservativeTracer

