  /*@@
   @file      primitive2conservative
   @date      Aug 31, 2010
   @author    Joshua Faber, Scott Noble, Bruno Mundim, Pedro Montero, Ian Hawke
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
#define Bvecx(i,j,k) Bprim(i,j,k,1)
#define Bvecy(i,j,k) Bprim(i,j,k,2)
#define Bvecz(i,j,k) Bprim(i,j,k,3)
#define Bconsx(i,j,k) Bcons(i,j,k,1)
#define Bconsy(i,j,k) Bcons(i,j,k,2)
#define Bconsz(i,j,k) Bcons(i,j,k,3)

#define DOT(x1,y1,z1,x2,y2,z2)   ( DOTP(gxx,gxy,gxz,gyy,gyz,gzz,x1,y1,z1,x2,y2,z2) )
#define DOT2(x1,y1,z1)           ( DOTP2(gxx,gxy,gxz,gyy,gyz,gzz,x1,y1,z1) )
#define DOTPT(x1,y1,z1,x2,y2,z2) ( DOTP(gxxpt,gxypt,gxzpt,gyypt,gyzpt,gzzpt,x1,y1,z1,x2,y2,z2) )
#define DOTPT2(x1,y1,z1)         ( DOTP2(gxxpt,gxypt,gxzpt,gyypt,gyzpt,gzzpt,x1,y1,z1) )

 /*@@
   @routine   primitive2conservativeM
   @date      Aug 31
   @author    Joshua Faber, Scott Noble, Bruno Mundim, Pedro Montero, Ian Hawke
   @desc 
   Converts primitive to conserved variables for the boundary extended data.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine primitive2conservativeM(CCTK_ARGUMENTS)

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  CCTK_INT, parameter :: izero = 0
  integer, parameter :: ik = kind(izero)

  integer :: i, j, k
  CCTK_REAL :: g11l,g12l,g13l,g22l,g23l,g33l,avg_sdetl,&
       g11r,g12r,g13r,g22r,g23r,g33r,avg_sdetr
  CCTK_REAL :: xtemp(1)
  character(len=256) NaN_WarnLine
  ! save memory when MP is not used
  CCTK_INT :: GRHydro_UseGeneralCoordinates
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3) :: g11, g12, g13, g22, g23, g33
  pointer (pg11,g11), (pg12,g12), (pg13,g13), (pg22,g22), (pg23,g23), (pg33,g33)
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3,3) :: vup
  pointer (pvup,vup)
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3,3) :: Bprim
  pointer (pBprim,Bprim)

  if (GRHydro_UseGeneralCoordinates(cctkGH).ne.0) then
    pg11 = loc(gaa)
    pg12 = loc(gab)
    pg13 = loc(gac)
    pg22 = loc(gbb)
    pg23 = loc(gbc)
    pg33 = loc(gcc)
    pvup = loc(lvel)
    pBprim = loc(lBvec)
  else
    pg11 = loc(gxx)
    pg12 = loc(gxy)
    pg13 = loc(gxz)
    pg22 = loc(gyy)
    pg23 = loc(gyz)
    pg33 = loc(gzz)
    pvup = loc(vel)
    pBprim = loc(Bvec)
  end if
  
  ! constraint transport needs to be able to average fluxes in the directions
  ! other that flux_direction, which in turn need the primitives on interfaces

  if(evolve_temper.ne.1) then

  !$OMP PARALLEL DO PRIVATE(k,j,i,avg_sdetl,avg_sdetr,&
  !$OMP                      g11l,g12l,g13l,g22l,g23l,g33l,&
  !$OMP                      g11r,g12r,g13r,g22r,g23r,g33r)
  do k = GRHydro_stencil, cctk_lsh(3)-GRHydro_stencil+1 + transport_constraints*(1-zoffset)
    do j = GRHydro_stencil, cctk_lsh(2)-GRHydro_stencil+1 + transport_constraints*(1-yoffset)
      do i = GRHydro_stencil, cctk_lsh(1)-GRHydro_stencil+1 + transport_constraints*(1-xoffset)
        
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

        call prim2conM(GRHydro_eos_handle, g11l,g12l,g13l,g22l,g23l,g33l, &
             avg_sdetl,densminus(i,j,k),sxminus(i,j,k),&
             syminus(i,j,k),szminus(i,j,k),tauminus(i,j,k),&
             Bconsxminus(i,j,k),Bconsyminus(i,j,k),Bconszminus(i,j,k), rhominus(i,j,k), &
             velxminus(i,j,k),velyminus(i,j,k),velzminus(i,j,k),&
             epsminus(i,j,k),pressminus(i,j,k),Bvecxminus(i,j,k), &
             Bvecyminus(i,j,k), Bveczminus(i,j,k), w_lorentzminus(i, j, k))
 
        call prim2conM(GRHydro_eos_handle, g11r,g12r,g13r,g22r,g23r,g33r, &
             avg_sdetr, densplus(i,j,k),sxplus(i,j,k),&
             syplus(i,j,k),szplus(i,j ,k),tauplus(i,j,k),&
             Bconsxplus(i,j,k),Bconsyplus(i,j,k),Bconszplus(i,j,k), &
             rhoplus(i,j,k),velxplus(i,j,k),velyplus(i,j,k),&
             velzplus(i,j,k),epsplus(i,j,k),pressplus(i,j,k),&
             Bvecxplus(i,j,k), Bvecyplus(i,j,k), Bveczplus(i,j,k), &
             w_lorentzplus(i,j,k)) 

        if(evolve_entropy.ne.0) then
          entropyconsminus(i,j,k) = avg_sdetl*entropyminus(i,j,k)*w_lorentzminus(i,j,k)
          entropyconsplus(i,j,k) = avg_sdetr*entropyplus(i,j,k)*w_lorentzplus(i,j,k)
        end if

      end do
    end do
  end do
  !$OMP END PARALLEL DO

  else

     !$OMP PARALLEL DO PRIVATE(i, j, k, avg_sdetl, avg_sdetr, xtemp,&
     !$OMP                      g11l,g12l,g13l,g22l,g23l,g33l, &
     !$OMP                      g11r,g12r,g13r,g22r,g23r,g33r)
     do k = GRHydro_stencil, cctk_lsh(3)-GRHydro_stencil+1 + transport_constraints*(1-zoffset)
       do j = GRHydro_stencil, cctk_lsh(2)-GRHydro_stencil+1 + transport_constraints*(1-yoffset)
         do i = GRHydro_stencil, cctk_lsh(1)-GRHydro_stencil+1 + transport_constraints*(1-xoffset)
              
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
              
              avg_sdetl = sqrt(SPATIAL_DETERMINANT(g11l,g12l,g13l,g22l,g23l,g33l))
              avg_sdetr = sqrt(SPATIAL_DETERMINANT(g11r,g12r,g13r,g22r,g23r,g33r))

              if (y_e_minus(i,j,k) .le. 0.0d0 .or. y_e_plus(i,j,k) .le. 0.0d0) then
              !$OMP CRITICAL
                  write(NaN_WarnLine,'(a100,7g15.6)') '(y_e_minus,y_e_plus,x,y,z,rho)', y_e(i,j,k), y_e_minus(i,j,k), y_e_plus(i,j,k), x(i,j,k),y(i,j,k),z(i,j,k),rho(i,j,k)
                  call CCTK_WARN(1, NaN_WarnLine)
              !$OMP END CRITICAL
              endif 
              call prim2conM_hot(GRHydro_eos_handle, GRHydro_reflevel,&
                   int(i,ik),int(j,ik),int(k,ik),x(i,j,k),y(i,j,k),z(i,j,k), g11l,g12l,g13l,g22l,g23l,g33l, &
                   avg_sdetl,densminus(i,j,k),sxminus(i,j,k),&
                   syminus(i,j,k),szminus(i,j,k),tauminus(i,j,k),&
                   Bconsxminus(i,j,k),Bconsyminus(i,j,k),Bconszminus(i,j,k), rhominus(i,j,k), &
                   velxminus(i,j,k),velyminus(i,j,k),velzminus(i,j,k),&
                   epsminus(i,j,k),pressminus(i,j,k),Bvecxminus(i,j,k), &
                   Bvecyminus(i,j,k), Bveczminus(i,j,k), w_lorentzminus(i, j, k), tempminus(i,j,k), y_e_minus(i,j,k))
  
              call prim2conM_hot(GRHydro_eos_handle, GRHydro_reflevel,&
                   int(i,ik),int(j,ik),int(k,ik),x(i,j,k),y(i,j,k),z(i,j,k), g11r,g12r,g13r,g22r,g23r,g33r, &
                   avg_sdetr, densplus(i,j,k),sxplus(i,j,k),&
                   syplus(i,j,k),szplus(i,j,k),tauplus(i,j,k),&
                   Bconsxplus(i,j,k),Bconsyplus(i,j,k),Bconszplus(i,j,k),&
                   rhoplus(i,j,k),velxplus(i,j,k),velyplus(i,j,k),&
                   velzplus(i,j,k),epsplus(i,j,k),pressplus(i,j,k),&
                   Bvecxplus(i,j,k), Bvecyplus(i,j,k), Bveczplus(i,j,k), &
                   w_lorentzplus(i,j,k), tempplus(i,j,k), y_e_plus(i,j,k)) 
             
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  endif

end subroutine primitive2conservativeM

 /*@@
   @routine    prim2conM
   @date       Aug 31, 2010
   @author     Joshua Faber, Scott Noble, Bruno Mundim, Pedro Montero, Ian Hawke
   @desc 
   Converts from primitive to conservative at a single point
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine prim2conM_hot(handle, GRHydro_reflevel, ii, jj, kk, &
     x, y, z, gxx, gxy, gxz, gyy, gyz, gzz, sdet, ddens, &
     dsx, dsy, dsz, dtau , dBconsx, dBconsy, dBconsz, drho, &
     dvelx, dvely, dvelz, deps, dpress, dBvcx, dBvcy, dBvcz, w, temp, ye) 
  
  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_REAL :: gxx, gxy, gxz, gyy, gyz, gzz, sdet
  CCTK_REAL, dimension(1) :: ddens, dsx, dsy, dsz, dtau, &
       dBconsx, dBconsy, dBconsz, &
       drho, dvelx, dvely, dvelz,&
       deps, dpress, dBvcx, dBvcy, dBvcz, vlowx, vlowy, vlowz   
  CCTK_REAL :: temp(1),ye(1), yein(1), x, y, z
  CCTK_INT :: handle, GRHydro_reflevel, ii, jj, kk
  CCTK_REAL :: w
  CCTK_REAL, dimension(1) :: Bdotv,ab0,b2,blowx,blowy,blowz
  character(len=256) NaN_WarnLine
  character(len=512) warnline

! begin EOS Omni vars
  CCTK_INT  :: n, keytemp, anyerr, keyerr(1)
 ! CCTK_REAL :: xpress(1),xeps(1),xtemp(1),xye(1)
  CCTK_REAL :: temp0(1)
  n = 1; keytemp = reconstruct_temper; anyerr = 0; keyerr(1) = 0
  !xpress = 0.0d0; xeps = 0.0d0; xtemp = 0.0d0; xye = 0.0d0
! end EOS Omni vars
  
  temp0 = temp
  w = 1.d0 / sqrt(1.d0 - DOT2(dvelx(1),dvely(1),dvelz(1)))

!!$ BEGIN: Check for NaN value
  if (w .ne. w) then
    !$OMP CRITICAL
    write(NaN_WarnLine,'(a100,6g15.6)') 'NaN produced in sqrt(): (gxx,gxy,gxz,gyy,gyz,gzz)', gxx, gxy, gxz, gyy, gyz, gzz
    call CCTK_WARN(1, NaN_WarnLine)
    write(NaN_WarnLine,'(a100,3g15.6)') 'NaN produced in sqrt(): (dvelx,dvely,dvelz)', dvelx, dvely, dvelz
    call CCTK_WARN(1, NaN_WarnLine)
    write(NaN_WarnLine,'(a100,3g15.6)') 'NaN produced in sqrt(): (x,y,z)', x, y, z
    call CCTK_WARN(1, NaN_WarnLine)
    write(NaN_WarnLine,'(a100,g15.6)') 'NaN produced in sqrt(): v2', DOT2(dvelx(1),dvely(1),dvelz(1))
    call CCTK_WARN(GRHydro_NaN_verbose, NaN_WarnLine)
    !$OMP END CRITICAL
  endif
!!$ END: Check for NaN value

  ye = max(min(ye,GRHydro_Y_e_max),GRHydro_Y_e_min)
 
  yein = ye 
 
  call EOS_Omni_press(handle,keytemp,GRHydro_eos_rf_prec,n,&
       drho,deps,temp,ye,dpress,keyerr,anyerr)  
  ! error handling
  if(anyerr.ne.0) then
     if(GRHydro_reflevel.lt.GRHydro_c2p_warn_from_reflevel) then
        ! in this case (coarse grid error that is hopefully restricted
        ! away), we use the average temperature between cells and call
        ! the EOS with keytemp=1
        keytemp=1-reconstruct_temper
        temp = temp0
        call EOS_Omni_press(handle,keytemp,GRHydro_eos_rf_prec,n,&
             drho,deps,temp,ye,dpress,keyerr,anyerr)  
        keytemp=reconstruct_temper
     else
        ! This is a way of recovering even on finer refinement levels:
        ! Use the average temperature at the interface instead of the
        ! reconstructed specific internal energy.
        keytemp=1-reconstruct_temper
        temp = temp0
        call EOS_Omni_press(handle,keytemp,GRHydro_eos_rf_prec,n,&
             drho,deps,temp,ye,dpress,keyerr,anyerr)  
        keytemp=reconstruct_temper
        if(anyerr.ne.0) then
           !$OMP CRITICAL
           call CCTK_WARN(1,"EOS error in prim2con_hot")
           write(warnline,"(3i5,1P10E15.6)") ii,jj,kk,x,y,z
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

  vlowx = gxx*dvelx + gxy*dvely + gxz*dvelz
  vlowy = gxy*dvelx + gyy*dvely + gyz*dvelz
  vlowz = gxz*dvelx + gyz*dvely + gzz*dvelz

  Bdotv=DOT(dvelx,dvely,dvelz,dBvcx,dBvcy,dBvcz)
  ab0=w*Bdotv
  b2 = DOT2(dBvcx,dBvcy,dBvcz)/w**2+Bdotv**2
  blowx = (gxx*dBvcx + gxy*dBvcy + gxz*dBvcz)/w + &
       w*Bdotv*vlowx
  blowy = (gxy*dBvcx + gyy*dBvcy + gyz*dBvcz)/w + &
       w*Bdotv*vlowy
  blowz = (gxz*dBvcx + gyz*dBvcy + gzz*dBvcz)/w + &
       w*Bdotv*vlowz

  ddens = sdet * drho * w 
  dsx = sdet * ((drho*(1+deps)+dpress+b2)*w*w * vlowx - &
       ab0*blowx)
  dsy = sdet * ((drho*(1+deps)+dpress+b2)*w*w * vlowy - &
       ab0*blowy)
  dsz = sdet * ((drho*(1+deps)+dpress+b2)*w*w * vlowz - &
       ab0*blowz)
  dtau = sdet * ((drho*(1+deps)+dpress+b2)*w*w - dpress-b2/2.0-ab0**2) - ddens 

  dBconsx = sdet*dBvcx
  dBconsy = sdet*dBvcy
  dBconsz = sdet*dBvcz

  !!$OMP CRITICAL
  !write(NaN_WarnLine,'(a100,3g15.6)') '(dens out, sqrt(det), w:)', ddens,sqrt(det),w
  !call CCTK_WARN(1, NaN_WarnLine)
  !!$OMP END CRITICAL
end subroutine prim2conM_hot


subroutine prim2conM(handle, gxx, gxy, gxz, gyy, gyz, gzz, sdet, ddens, &
     dsx, dsy, dsz, dtau , dBconsx, dBconsy, dBconsz, drho, &
     dvelx, dvely, dvelz, deps, dpress, dBvcx, dBvcy, dBvcz, w) 
  
  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_REAL :: gxx, gxy, gxz, gyy, gyz, gzz, sdet
  CCTK_REAL, dimension(1) :: ddens, dsx, dsy, dsz, dtau, &
       dBconsx, dBconsy, dBconsz, &
       drho, dvelx, dvely, dvelz,&
       deps, dpress, dBvcx, dBvcy, dBvcz, vlowx, vlowy, vlowz   
  CCTK_REAL :: w
  CCTK_REAL, dimension(1) :: Bdotv,ab0,b2,blowx,blowy,blowz
  CCTK_INT :: handle
  character(len=256) NaN_WarnLine

! begin EOS Omni vars
  CCTK_INT  :: n, keytemp, anyerr, keyerr(1)
  CCTK_REAL :: xpress(1),xeps(1),xtemp(1),xye(1)
  n = 1; keytemp = 0; anyerr = 0; keyerr(1) = 0
  xpress = 0.0d0; xeps = 0.0d0; xtemp = 0.0d0; xye = 0.0d0
! end EOS Omni vars
  
  w = 1.d0 / sqrt(1.d0 - DOT2(dvelx(1),dvely(1),dvelz(1)))

!!$ BEGIN: Check for NaN value
  if (w .ne. w) then
    !$OMP CRITICAL
    write(NaN_WarnLine,'(a100,6g15.6)') 'NaN produced in sqrt(): (gxx,gxy,gxz,gyy,gyz,gzz)', gxx, gxy, gxz, gyy, gyz, gzz
    call CCTK_WARN(GRHydro_NaN_verbose, NaN_WarnLine)
    write(NaN_WarnLine,'(a100,3g15.6)') 'NaN produced in sqrt(): (dvelx,dvely,dvelz)', dvelx, dvely, dvelz
    call CCTK_WARN(GRHydro_NaN_verbose, NaN_WarnLine)
    !$OMP END CRITICAL
  endif
!!$ END: Check for NaN value

  call EOS_Omni_press(handle,keytemp,GRHydro_eos_rf_prec,n,&
       drho,deps,xtemp,xye,dpress,keyerr,anyerr)  

  vlowx = gxx*dvelx + gxy*dvely + gxz*dvelz
  vlowy = gxy*dvelx + gyy*dvely + gyz*dvelz
  vlowz = gxz*dvelx + gyz*dvely + gzz*dvelz

  Bdotv=DOT(dvelx,dvely,dvelz,dBvcx,dBvcy,dBvcz)
  ab0=w*Bdotv
  b2 = DOT2(dBvcx,dBvcy,dBvcz)/w**2+Bdotv**2
  blowx = (gxx*dBvcx + gxy*dBvcy + gxz*dBvcz)/w + &
       w*Bdotv*vlowx
  blowy = (gxy*dBvcx + gyy*dBvcy + gyz*dBvcz)/w + &
       w*Bdotv*vlowy
  blowz = (gxz*dBvcx + gyz*dBvcy + gzz*dBvcz)/w + &
       w*Bdotv*vlowz

  ddens = sdet * drho * w 
  dsx = sdet * ((drho*(1+deps)+dpress+b2)*w*w * vlowx - &
       ab0*blowx)
  dsy = sdet * ((drho*(1+deps)+dpress+b2)*w*w * vlowy - &
       ab0*blowy)
  dsz = sdet * ((drho*(1+deps)+dpress+b2)*w*w * vlowz - &
       ab0*blowz)
  dtau = sdet * ((drho*(1+deps)+dpress+b2)*w*w - dpress-b2/2.0-ab0**2) - ddens 

  dBconsx = sdet*dBvcx
  dBconsy = sdet*dBvcy
  dBconsz = sdet*dBvcz

end subroutine prim2conM


 /*@@
   @routine    Primitive2ConservativeCellsM
   @date       Aug 31, 2010
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

subroutine Primitive2ConservativeCellsM(CCTK_ARGUMENTS)

  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  CCTK_INT :: i, j, k
  CCTK_REAL :: sdet
  CCTK_REAL :: maxtau0
  ! save memory when MP is not used
  CCTK_INT :: GRHydro_UseGeneralCoordinates
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3) :: g11, g12, g13, g22, g23, g33
  pointer (pg11,g11), (pg12,g12), (pg13,g13), (pg22,g22), (pg23,g23), (pg33,g33)
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3,3) :: vup
  pointer (pvup,vup)
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3,3) :: Bprim
  pointer (pBprim,Bprim)

  if (GRHydro_UseGeneralCoordinates(cctkGH).ne.0) then
    pg11 = loc(gaa)
    pg12 = loc(gab)
    pg13 = loc(gac)
    pg22 = loc(gbb)
    pg23 = loc(gbc)
    pg33 = loc(gcc)
    pvup = loc(lvel)
    pBprim = loc(lBvec)
  else
    pg11 = loc(gxx)
    pg12 = loc(gxy)
    pg13 = loc(gxz)
    pg22 = loc(gyy)
    pg23 = loc(gyz)
    pg33 = loc(gzz)
    pvup = loc(vel)
    pBprim = loc(Bvec)
  end if
  
  
  maxtau0 = -1.0d60
  
  if(evolve_temper.ne.1) then

  !$OMP PARALLEL DO PRIVATE(k,j,i,sdet), REDUCTION(MAX:maxtau0)
  do k = 1,cctk_lsh(3)
    do j = 1,cctk_lsh(2)
      do i = 1,cctk_lsh(1)
        
         sdet = sdetg(i,j,k)

        call prim2conM(GRHydro_eos_handle,g11(i,j,k),&
             g12(i,j,k),g13(i,j,k),&
             g22(i,j,k),g23(i,j,k),g33(i,j,k),&
             sdet, dens(i,j,k),sx(i,j,k),sy(i,j,k),sz(i,j,k),&
             tau(i,j,k),Bconsx(i,j,k),Bconsy(i,j,k),Bconsz(i,j,k),&
             rho(i,j,k),velx(i,j,k),vely(i,j,k),velz(i,j,k),&
             eps(i,j,k),press(i,j,k),Bvecx(i,j,k), &
             Bvecy(i,j,k), Bvecz(i,j,k), w_lorentz(i,j,k))

        if(evolve_entropy.ne.0) then
          entropycons(i,j,k) = sdet*entropy(i,j,k)*w_lorentz(i,j,k)
        end if

        maxtau0 = max(maxtau0,tau(i,j,k))

      end do
    end do
  end do
  !$OMP END PARALLEL DO

  ! TODO: to actually reduce GRHydro_tau_min over all Carpet components
  ! we need to modify Carpet looping to allow this function to be called
  ! on all AMR levels before calling any other function. The best would be
  ! to set a special bin where functions would be called on all levels first
  ! instead of calling all functions per level. The workaround for this problem
  ! is to set GRHydro_tau_min to a user specified value as it was set in 
  ! GRHydro_Minima.F90. Once this issue is solved, then uncomment the line
  ! below and create two other routines to be run in global mode so that 
  ! GRHydro_tau_min can be properly initialized and reduced.

  !GRHydro_tau_min = GRHydro_tau_min * maxtau0 

  else 
  !$OMP PARALLEL DO PRIVATE(k,j,i,sdet)
  do k = 1,cctk_lsh(3)
    do j = 1,cctk_lsh(2)
      do i = 1,cctk_lsh(1)
        
         sdet = sdetg(i,j,k)

         call prim2conM_hot(GRHydro_eos_handle,GRHydro_reflevel,&
              i,j,k,x(i,j,k),y(i,j,k),z(i,j,k),&
              g11(i,j,k),g12(i,j,k),g13(i,j,k),&
              g22(i,j,k),g23(i,j,k),g33(i,j,k),&
              sdet, dens(i,j,k),sx(i,j,k),sy(i,j,k),sz(i,j,k),&
              tau(i,j,k),Bconsx(i,j,k),Bconsy(i,j,k),Bconsz(i,j,k),&
              rho(i,j,k),velx(i,j,k),vely(i,j,k),velz(i,j,k),&
              eps(i,j,k),press(i,j,k),Bvecx(i,j,k), &
              Bvecy(i,j,k), Bvecz(i,j,k), w_lorentz(i,j,k), temperature(i,j,k), Y_e(i,j,k))

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



end subroutine Primitive2ConservativeCellsM


 /*@@
   @routine    Prim2ConservativePolytypeM
   @date       Aug 31, 2010
   @author     Joshua Faber, Scott Noble, Bruno Mundim, Ian Hawke
   @desc 
   Same as first routine, only for polytropes.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/


subroutine Prim2ConservativePolytypeM(CCTK_ARGUMENTS)

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
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3,3) :: Bprim
  pointer (pBprim,Bprim)

  if (GRHydro_UseGeneralCoordinates(cctkGH).ne.0) then
    pg11 = loc(gaa)
    pg12 = loc(gab)
    pg13 = loc(gac)
    pg22 = loc(gbb)
    pg23 = loc(gbc)
    pg33 = loc(gcc)
    pvup = loc(lvel)
    pBprim = loc(lBvec)
  else
    pg11 = loc(gxx)
    pg12 = loc(gxy)
    pg13 = loc(gxz)
    pg22 = loc(gyy)
    pg23 = loc(gyz)
    pg33 = loc(gzz)
    pvup = loc(vel)
    pBprim = loc(Bvec)
  end if
  
  
  ! constraint transport needs to be able to average fluxes in the directions
  ! other that flux_direction, which in turn need the primitives on interfaces
  !$OMP PARALLEL DO PRIVATE(i, j, k, g11l,g12l,g13l,g22l,g23l,g33l,avg_sdetl,&
  !$OMP g11r,g12r,g13r,g22r,g23r,g33r,avg_sdetr)
  do k = GRHydro_stencil, cctk_lsh(3)-GRHydro_stencil+1 + transport_constraints*(1-zoffset)
    do j = GRHydro_stencil, cctk_lsh(2)-GRHydro_stencil+1 + transport_constraints*(1-yoffset)
      do i = GRHydro_stencil, cctk_lsh(1)-GRHydro_stencil+1 + transport_constraints*(1-xoffset)
        
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

        call prim2conpolytypeM(GRHydro_eos_handle, g11l,g12l,g13l,&
             g22l,g23l,g33l, &
             avg_sdetl,densminus(i,j,k),sxminus(i,j,k),&
             syminus(i,j,k),szminus(i,j,k),tauminus(i,j,k),&
             Bconsxminus(i,j,k),Bconsyminus(i,j,k),Bconszminus(i,j,k),&
             rhominus(i,j,k),&
             velxminus(i,j,k),velyminus(i,j,k),velzminus(i,j,k),&
             epsminus(i,j,k),pressminus(i,j,k),Bvecxminus(i,j,k), &
             Bvecyminus(i,j,k),Bveczminus(i,j,k),w_lorentzminus(i, j, k))
 
        call prim2conpolytypeM(GRHydro_eos_handle, g11r,g12r,g13r,&
             g22r,g23r,g33r, &
             avg_sdetr,densplus(i,j,k),sxplus(i,j,k),&
             syplus(i,j,k),szplus(i,j,k),tauplus(i,j,k),&
             Bconsxplus(i,j,k),Bconsyplus(i,j,k),Bconszplus(i,j,k),&
             rhoplus(i,j,k),&
             velxplus(i,j,k),velyplus(i,j,k),velzplus(i,j,k),&
             epsplus(i,j,k),pressplus(i,j,k),Bvecxplus(i,j,k), &
             Bvecyplus(i,j,k),Bveczplus(i,j,k),w_lorentzplus(i, j, k))

        if (densminus(i,j,k) .ne. densminus(i,j,k)) then
        !$OMP CRITICAL
        call CCTK_WARN(1, "NaN in densminus(i,j,k) (Prim2Con)")
        !$OMP END CRITICAL
        endif
      end do
    end do
  end do
  !$OMP END PARALLEL DO
end subroutine Prim2ConservativePolytypeM

 /*@@
   @routine    prim2conpolytypeM
   @date       Aug 31, 2010
   @author     Joshua Faber, Scott Noble, Bruno Mundim, Pedro Montero, Ian Hawke
   @desc 
   Converts from primitive to conservative at a single point
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine prim2conpolytypeM(handle, gxx, gxy, gxz, gyy, gyz, &
     gzz, sdet, ddens, dsx, dsy, dsz, dtau, dBconsx, dBconsy, dBconsz, &
     drho, dvelx, dvely, dvelz, deps, dpress, dBvcx, dBvcy, dBvcz, w) 
  
  implicit none
  
  DECLARE_CCTK_PARAMETERS
  !!DECLARE_CCTK_FUNCTIONS

  CCTK_REAL :: gxx, gxy, gxz, gyy, gyz, gzz, sdet
  CCTK_REAL :: ddens, dsx, dsy, dsz, dtau, dBconsx, dBconsy, dBconsz, &
       drho, dvelx, dvely, dvelz, deps, dpress, dBvcx, dBvcy, dBvcz, &
       w_tmp, w, vlowx, vlowy, vlowz
  CCTK_INT :: handle
  CCTK_REAL :: Bdotv,ab0,b2,blowx,blowy,blowz
  character(len=256) NaN_WarnLine

! begin EOS Omni vars
  CCTK_INT  :: n, keytemp, anyerr, keyerr(1)
  CCTK_REAL :: xpress,xeps,xtemp,xye
  n = 1; keytemp = 0; anyerr = 0; keyerr(1) = 0
  xpress = 0.0d0; xeps = 0.0d0; xtemp = 0.0d0; xye = 0.0d0
! end EOS Omni vars
  
  w_tmp = DOT2(dvelx,dvely,dvelz)

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

  Bdotv=DOT(dvelx,dvely,dvelz,dBvcx,dBvcy,dBvcz)
  ab0=w*Bdotv
  b2 = DOT2(dBvcx,dBvcy,dBvcz)/w**2+Bdotv**2
  blowx = (gxx*dBvcx + gxy*dBvcy + gxz*dBvcz)/w + &
       w*Bdotv*vlowx
  blowy = (gxy*dBvcx + gyy*dBvcy + gyz*dBvcz)/w + &
       w*Bdotv*vlowy
  blowz = (gxz*dBvcx + gyz*dBvcy + gzz*dBvcz)/w + &
       w*Bdotv*vlowz

  ddens = sdet * drho * w 
  dsx = sdet * ((drho*(1+deps)+dpress+b2)*w*w * vlowx - &
       ab0*blowx)
  dsy = sdet * ((drho*(1+deps)+dpress+b2)*w*w * vlowy - &
       ab0*blowy)
  dsz = sdet * ((drho*(1+deps)+dpress+b2)*w*w * vlowz - &
       ab0*blowz)
  dtau = sdet * ((drho*(1+deps)+dpress+b2)*w*w - dpress-b2/2.0-ab0**2) - ddens 
  
  dBconsx = sdet*dBvcx
  dBconsy = sdet*dBvcy
  dBconsz = sdet*dBvcz

end subroutine prim2conpolytypeM


 /*@@
   @routine    Primitive2ConservativePolyCellsM
   @date       Aug 31, 2010
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


subroutine Primitive2ConservativePolyCellsM(CCTK_ARGUMENTS)

  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  CCTK_INT :: i, j, k
  CCTK_REAL :: sdet
  ! save memory when MP is not used
  CCTK_INT :: GRHydro_UseGeneralCoordinates
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3) :: g11, g12, g13, g22, g23, g33
  pointer (pg11,g11), (pg12,g12), (pg13,g13), (pg22,g22), (pg23,g23), (pg33,g33)
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3,3) :: vup
  pointer (pvup,vup)
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3,3) :: Bprim
  pointer (pBprim,Bprim)

  if (GRHydro_UseGeneralCoordinates(cctkGH).ne.0) then
    pg11 = loc(gaa)
    pg12 = loc(gab)
    pg13 = loc(gac)
    pg22 = loc(gbb)
    pg23 = loc(gbc)
    pg33 = loc(gcc)
    pvup = loc(lvel)
    pBprim = loc(lBvec)
  else
    pg11 = loc(gxx)
    pg12 = loc(gxy)
    pg13 = loc(gxz)
    pg22 = loc(gyy)
    pg23 = loc(gyz)
    pg33 = loc(gzz)
    pvup = loc(vel)
    pBprim = loc(Bvec)
  end if
  
  !$OMP PARALLEL DO PRIVATE(k,j,i,sdet)
  do k = 1,cctk_lsh(3)
    do j = 1,cctk_lsh(2)
      do i = 1,cctk_lsh(1)
        
        sdet = sdetg(i,j,k)

        call prim2conpolytypeM(GRHydro_eos_handle,g11(i,j,k),g12(i,j,k),&
             g13(i,j,k),g22(i,j,k),g23(i,j,k),g33(i,j,k),&
             sdet, dens(i,j,k),sx(i,j,k),sy(i,j,k),sz(i,j,k),&
             tau(i,j,k),Bconsx(i,j,k),Bconsy(i,j,k),Bconsz(i,j,k),&
             rho(i,j,k),velx(i,j,k),vely(i,j,k),velz(i,j,k),&
             eps(i,j,k),press(i,j,k),Bvecx(i,j,k), Bvecy(i,j,k), &
             Bvecz(i,j,k), w_lorentz(i,j,k))

      end do
    end do
  end do
  !$OMP END PARALLEL DO

  if(evolve_Y_e.ne.0) then
     !$OMP PARALLEL DO PRIVATE(i, j, k)
     do k = GRHydro_stencil,cctk_lsh(3)-GRHydro_stencil+1
        do j = GRHydro_stencil,cctk_lsh(2)-GRHydro_stencil+1
           do i = GRHydro_stencil,cctk_lsh(1)-GRHydro_stencil+1
              Y_e_con(i,j,k) = Y_e(i,j,k) * dens(i,j,k)
           enddo
        enddo
     enddo
     !$OMP END PARALLEL DO
  endif


end subroutine Primitive2ConservativePolyCellsM

!!$
!!$ Prim2Con doesn't change for tracers with the addition of a B-field!
!!$

