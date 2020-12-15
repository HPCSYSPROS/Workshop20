#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "SpaceMask.h"
#include "GRHydro_Macros.h"


subroutine Conservative2PrimitiveHot(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  CCTK_INT, parameter :: izero = 0
  integer, parameter :: ik = kind(izero)

  integer :: i, j, k, itracer, nx, ny, nz, myproc
  CCTK_REAL :: uxx, uxy, uxz, uyy, uyz, uzz, sdet, pmin, epsmin, dummy1, dummy2
  CCTK_REAL :: vlowx, vlowy, vlowz
  logical :: epsnegative
  character*512 :: warnline
  
  CCTK_REAL :: local_min_tracer
  CCTK_REAL :: local_perc_ptol
  integer :: reflevel

  ! save memory when MP is not used
  CCTK_INT :: GRHydro_UseGeneralCoordinates
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3) :: g11, g12, g13, g22, g23, g33
  pointer (pg11,g11), (pg12,g12), (pg13,g13), (pg22,g22), (pg23,g23), (pg33,g33)
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3,3) :: vup
  pointer (pvup,vup)

! begin EOS Omni vars
  CCTK_INT  :: n,keytemp,anyerr,keyerr(1), reset_to_atmo
  CCTK_REAL :: xpress(1),xeps(1),xtemp(1),xye(1),xrho(1)
  n = 1;keytemp = 0;anyerr = 0;keyerr(1) = 0
  xpress = 0.0d0;xeps = 0.0d0;xtemp = 0.0d0;xye = 0.0d0
! end EOS Omni vars

  myproc = CCTK_MyProc(cctkGH)
  reflevel = int(log10(dble(cctk_levfac(1)))/log10(2.0d0))

  ! save memory when MP is not used
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

  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)
  
  if (use_min_tracer .ne. 0) then
    local_min_tracer = min_tracer
  else
    local_min_tracer = 0.0d0
  end if

  ! this is a poly call
  xrho(1) = GRHydro_rho_min
  call EOS_Omni_press(GRHydro_polytrope_handle,keytemp,GRHydro_eos_rf_prec,n,&
       xrho,xeps,xtemp,xye,xpress,keyerr,anyerr)
  call EOS_Omni_EpsFromPress(GRHydro_polytrope_handle,keytemp,GRHydro_eos_rf_prec,n,&
       xrho,xeps,xtemp,xye,xpress,xeps,keyerr,anyerr)
  pmin = xpress(1)
  epsmin = xeps(1)

  !$OMP PARALLEL DO PRIVATE(i,j,k,itracer,&
  !$OMP uxx, uxy, uxz, uyy, uyz, uzz, sdet, epsnegative, anyerr, keyerr, keytemp,&
  !$OMP warnline, dummy1, dummy2,reset_to_atmo)
  do k = 1, nz 
    do j = 1, ny 
      do i = 1, nx


        !do not compute if in atmosphere
        if (atmosphere_mask(i,j,k) .gt. 0) cycle
         
        epsnegative = .false.

        sdet = sdetg(i,j,k)

        call UpperMetric(uxx,uxy,uxz,uyy,uyz,uzz,sdet*sdet,&
             g11(i,j,k),g12(i,j,k),g13(i,j,k),g22(i,j,k),&
             g23(i,j,k),g33(i,j,k))        

        if (evolve_tracer .ne. 0) then
           do itracer=1,number_of_tracers
              call Con2Prim_ptTracer(cons_tracer(i,j,k,itracer), tracer(i,j,k,itracer), &
                   dens(i,j,k))

              if (use_min_tracer .ne. 0) then
                if (tracer(i,j,k,itracer) .le. local_min_tracer) then
                  tracer(i,j,k,itracer) = local_min_tracer
                end if
              end if

           enddo
           
        endif
        
        if(evolve_Y_e.ne.0) then
           Y_e(i,j,k) = max(min(Y_e_con(i,j,k) / dens(i,j,k),GRHydro_Y_e_max),&
                GRHydro_Y_e_min)
        endif

        reset_to_atmo = 0
        IF_BELOW_ATMO(dens(i,j,k), sdet*GRHydro_rho_min, GRHydro_atmo_tolerance, r(i,j,k)) then
           reset_to_atmo = 1
        endif

        if (reset_to_atmo .gt. 0 .or. (GRHydro_enable_internal_excision /= 0 .and. hydro_excision_mask(i,j,k) .gt. 0)) then
           SET_ATMO_MIN(dens(i,j,k), sdet*GRHydro_rho_min, r(i,j,k)) 
           SET_ATMO_MIN(rho(i,j,k), GRHydro_rho_min, r(i,j,k)) 
           scon(i,j,k,:) = 0.d0
           vup(i,j,k,:) = 0.d0
           w_lorentz(i,j,k) = 1.d0

           ! set hot atmosphere values
           temperature(i,j,k) = grhydro_hot_atmo_temp
           y_e(i,j,k) = grhydro_hot_atmo_Y_e
           y_e_con(i,j,k) = y_e(i,j,k) * dens(i,j,k)
           keytemp = 1
           call EOS_Omni_press(GRHydro_eos_handle,keytemp,GRHydro_eos_rf_prec,n,&
                rho(i,j,k),eps(i,j,k),temperature(i,j,k),y_e(i,j,k),&
                press(i,j,k),keyerr,anyerr)
           keytemp = 0
           ! w_lorentz=1, so the expression for tau reduces to:
           tau(i,j,k)  = sdet * (rho(i,j,k)+rho(i,j,k)*eps(i,j,k)) - dens(i,j,k)
           
           cycle

          end if

          call Con2Prim_pt_hot3(int(cctk_iteration,ik),myproc,int(i,ik),int(j,ik),int(k,ik),GRHydro_eos_handle,&
               dens(i,j,k),scon(i,j,k,1),&
               scon(i,j,k,2),scon(i,j,k,3),tau(i,j,k),Y_e_con(i,j,k),rho(i,j,k),vup(i,j,k,1),&
               vup(i,j,k,2),vup(i,j,k,3),eps(i,j,k),temperature(i,j,k),y_e(i,j,k),&
               press(i,j,k),w_lorentz(i,j,k), &
               uxx,uxy,uxz,uyy,uyz,uzz,sdet,x(i,j,k),y(i,j,k), &
               z(i,j,k),r(i,j,k),epsnegative,GRHydro_rho_min,pmin, epsmin, & 
               int(reflevel,ik), GRHydro_C2P_failed(i,j,k), GRHydro_perc_ptol)

          if(temperature(i,j,k).gt.GRHydro_max_temp) then
             !$OMP CRITICAL
             call CCTK_WARN(1,"C2P: Temperature too high")
             write(warnline,"(i8)") cctk_iteration
             call CCTK_WARN(1,warnline)
             write(warnline,"(4i5,1P10E15.6)") GRHydro_Reflevel,i,j,k,x(i,j,k),y(i,j,k),z(i,j,k)
             call CCTK_WARN(1,warnline)
             write(warnline,"(1P10E15.6)") dens(i,j,k),scon(i,j,k,1:3),&
                  tau(i,j,k),w_lorentz(i,j,k)
             call CCTK_WARN(1,warnline)
             write(warnline,"(1P10E15.6)") rho(i,j,k),eps(i,j,k),&
                  temperature(i,j,k),Y_e(i,j,k)
             call CCTK_WARN(1,warnline)
             write(warnline,"(A7,i8)") "code: ",keyerr(1)
             call CCTK_WARN(1,warnline)
             write(warnline,"(A10,i5)") "reflevel: ", reflevel
             call CCTK_WARN(1,warnline)
             call CCTK_ERROR("Aborting!!!")
             STOP
             !$OMP END CRITICAL
          endif


          if( abs(GRHydro_C2P_failed(i,j,k)-2.0d0) .lt. 1.0d-10) then
             ! this means c2p did not converge.
             ! In this case, we attempt to call c2p with a reduced
             ! accuracy requirement; if it fails again, we abort
             GRHydro_C2P_failed(i,j,k) = 0
             if(temperature(i,j,k) .gt. GRHydro_hot_atmo_temp) then
                local_perc_ptol = GRHydro_perc_ptol*100.0d0
             else
                ! If we are in the extrapolation regime for the EOS,
                ! we accept a larger c2p error, since we will be resetting
                ! the temperature anyway. The error scale here is chosen 
                ! somewhat arbitrarily to be 0.01% in pressnew/pressold.
                local_perc_ptol = 1.0d-4
             endif
             call Con2Prim_pt_hot3(int(cctk_iteration,ik),myproc,int(i,ik),int(j,ik),int(k,ik),GRHydro_eos_handle,&
                  dens(i,j,k),scon(i,j,k,1),&
                  scon(i,j,k,2),scon(i,j,k,3),tau(i,j,k),Y_e_con(i,j,k),rho(i,j,k),&
                  vup(i,j,k,1),vup(i,j,k,2), vup(i,j,k,3),eps(i,j,k),&
                  temperature(i,j,k),y_e(i,j,k),press(i,j,k),w_lorentz(i,j,k), &
                  uxx,uxy,uxz,uyy,uyz,uzz,sdet,x(i,j,k),y(i,j,k), &
                  z(i,j,k),r(i,j,k),epsnegative,GRHydro_rho_min,pmin, epsmin, & 
                  int(reflevel,ik), GRHydro_C2P_failed(i,j,k), local_perc_ptol)
             if( abs(GRHydro_C2P_failed(i,j,k)-2.0d0) .lt. 1.0d-10) then
                !$OMP CRITICAL
                if (reflevel.ge.GRHydro_c2p_warn_from_reflevel) then
                   call CCTK_WARN(1,"Convergence problem in c2p")
                   write(warnline,"(A10,i5)") "reflevel: ",reflevel
                   call CCTK_WARN(1,warnline)
                   write(warnline,"(A10,i10)") "iteration: ",cctk_iteration
                   call CCTK_WARN(1,warnline)
                   write(warnline,"(A10,F5.2)") "error: ", GRHydro_C2P_failed(i,j,k)
                   call CCTK_WARN(1,warnline)
                   write(warnline,"(3i5,1P10E15.6)") i,j,k,x(i,j,k),y(i,j,k),z(i,j,k),r(i,j,k)
                   call CCTK_WARN(1,warnline)
                   write(warnline,"(1P10E15.6)") rho(i,j,k),dens(i,j,k),eps(i,j,k),&
                        temperature(i,j,k),y_e(i,j,k)
                   call CCTK_WARN(1,warnline)
                   call CCTK_ERROR("Aborting!!!")
                   STOP
                endif
                !$OMP END CRITICAL
             endif
          endif

          if( abs(GRHydro_C2P_failed(i,j,k)-3.0d0) .lt. 1.0d-10 .or. &
               temperature(i,j,k).lt.GRHydro_hot_atmo_temp) then
             ! dropped out off the EOS table in temperature or below
             ! the temperature of the atmosphere.
             ! Now reset this point to minimum temperature.
             GRHydro_C2P_failed(i,j,k) = 0
             if (reflevel.ge.GRHydro_c2p_warn_from_reflevel) then
                !$OMP CRITICAL
                write(warnline,"(A10,i7,4i5,1P10E15.6)") "reset T:",&
                     cctk_iteration,GRHydro_Reflevel,&
                     i,j,k,rho(i,j,k),&
                     eps(i,j,k),temperature(i,j,k),y_e(i,j,k)
                call CCTK_WARN(1,warnline)
                write(warnline,"(A10,i7,4i5,1P10E15.6)") "reset T:",&
                     cctk_iteration,GRHydro_Reflevel,&
                     i,j,k,x(i,j,k),y(i,j,k),z(i,j,k),r(i,j,k)
                call CCTK_WARN(1,warnline)
                !$OMP END CRITICAL
             endif
             temperature(i,j,k) = GRHydro_hot_atmo_temp
             keytemp = 1
             call EOS_Omni_press(GRHydro_eos_handle,keytemp,GRHydro_eos_rf_prec,n,&
                  rho(i,j,k),eps(i,j,k),temperature(i,j,k),y_e(i,j,k),&
                  press(i,j,k),keyerr,anyerr)
             keytemp = 0
             if(anyerr.ne.0) then
                !$OMP CRITICAL
                call CCTK_ERROR("EOS Problem in C2P hot!!!")
                STOP
                !$OMP END CRITICAL
             endif
             ! prim2con
             w_lorentz(i,j,k) = &
                  1.d0 / sqrt(1.d0 - (gxx(i,j,k)*vel(i,j,k,1)**2 &
                  + gyy(i,j,k)*vel(i,j,k,2)**2 &
                  + gzz(i,j,k) *vel(i,j,k,3)**2 &
                  + 2.0d0*gxy(i,j,k)*vel(i,j,k,1)*vel(i,j,k,2) &
                  + 2.0d0*gxz(i,j,k)*vel(i,j,k,1)*vel(i,j,k,3)  &
                  + 2.0d0*gyz(i,j,k)*vel(i,j,k,2)*vel(i,j,k,3)))
             vlowx = gxx(i,j,k)*vel(i,j,k,1) &
                  + gxy(i,j,k)*vel(i,j,k,2)  &
                  + gxz(i,j,k)*vel(i,j,k,3)
             vlowy = gxy(i,j,k)*vel(i,j,k,1) &
                  + gyy(i,j,k)*vel(i,j,k,2)  &
                  + gyz(i,j,k)*vel(i,j,k,3)
             vlowz = gxz(i,j,k)*vel(i,j,k,1) &
                  + gyz(i,j,k)*vel(i,j,k,2)  &
                  + gzz(i,j,k)*vel(i,j,k,3)
             scon(i,j,k,1) = sdet * (rho(i,j,k)*(1.0d0+eps(i,j,k)) &
                  + press(i,j,k))*w_lorentz(i,j,k)**2 * vlowx
             scon(i,j,k,2) = sdet * (rho(i,j,k)*(1.0d0+eps(i,j,k)) &
                  +press(i,j,k))*w_lorentz(i,j,k)**2 * vlowy
             scon(i,j,k,3) = sdet * (rho(i,j,k)*(1.0d0+eps(i,j,k)) &
                  + press(i,j,k))*w_lorentz(i,j,k)**2 * vlowz
             tau(i,j,k) = sdet * ((rho(i,j,k)*(1.0d0+eps(i,j,k)) &
                  + press(i,j,k))*w_lorentz(i,j,k)**2 - press(i,j,k)) &
                  - dens(i,j,k)
          endif

      end do
    end do
  end do
  !$OMP END PARALLEL DO

  return
  
end subroutine Conservative2PrimitiveHot


subroutine Con2Prim_pt_hot3(cctk_iteration, myproc, ii,jj,kk,handle, dens, &
     sx, sy, sz, tau, ye_con, rho, velx, vely, &
     velz, epsilon, temp, ye, press, w_lorentz, uxx, uxy, uxz, uyy, &
     uyz, uzz, sdet, x, y, z, r, epsnegative, GRHydro_rho_min, pmin, epsmin, &
     GRHydro_reflevel, GRHydro_C2P_failed, local_perc_ptol)
  
  implicit none
  
  DECLARE_CCTK_PARAMETERS
  !!DECLARE_CCTK_FUNCTIONS

  CCTK_REAL dens, sx, sy, sz, tau, ye_con, rho, velx, vely, velz, epsilon, &
       press, uxx, uxy, uxz, uyy, uyz, uzz, invsdet, sdet, w_lorentz, x, &
       y, z, r, GRHydro_rho_min
  CCTK_REAL temp, ye
  CCTK_REAL s2, f, df, vlowx, vlowy, vlowz
  CCTK_INT cctk_iteration, ii,jj,kk,count, i, handle, GRHydro_reflevel
  CCTK_REAL GRHydro_C2P_failed
  CCTK_REAL udens, usx, usy, usz, utau, pold, pnew, &
            temp1, drhobydpress, depsbydpress, dpressbydeps, dpressbydrho, pmin, epsmin
  CCTK_REAL pminl,plow,tmp, dummy1, dummy2
  CCTK_REAL local_perc_ptol
  CCTK_REAL dpf

  integer :: myproc

  character(len=256) warnline
  logical epsnegative, mustbisect

  integer :: failwarnmode 
  integer :: failinfomode 

! begin EOS Omni vars
  CCTK_INT  :: n,keytemp,anyerr,keyerr(1)
  CCTK_REAL :: xpress,temp0
  CCTK_INT  :: nfudgemax,nf
  n=1;keytemp=0;anyerr=0;keyerr(1)=0
  temp0 = 0.0d0;xpress = 0.0d0
  nf=0;nfudgemax=30
! end EOS Omni vars

  mustbisect = .false.

  failinfomode = 1
  failwarnmode = 0

! set pmin to something sensible:
  pminl = 1.0d-28

  if(con2prim_oct_hack.ne.0.and.&
       (x .lt. -1.0d-12 .or.&
        y .lt. -1.0d-12 .or.&
        z .lt. -1.0d-12)) then
     failwarnmode = 2
     failinfomode = 2  
  endif

!!$  Undensitize the variables 
  invsdet = 1.0d0/sdet
  udens = dens * invsdet
  usx = sx * invsdet
  usy = sy * invsdet
  usz = sz * invsdet
  utau = tau * invsdet
  s2 = usx*usx*uxx + usy*usy*uyy + usz*usz*uzz + 2.*usx*usy*uxy + &
       2.*usx*usz*uxz + 2.*usy*usz*uyz

!!$  Set initial guess for pressure:
  temp0 = temp
  call EOS_Omni_press(handle,keytemp,GRHydro_eos_rf_prec,n,&
                      rho,epsilon,temp,ye,xpress,keyerr,anyerr)
  ! This is the lowest admissible pressure, otherwise we are off the physical regime
  ! Sometimes, we may end up with bad initial guesses. 
  ! In that case we need to start off with something
  ! reasonable at least
  plow = max(pminl, sqrt(s2) - utau - udens)
  ! Start out with some reasonable initial guess
  pold = max(plow+1.d-10,xpress)
  ! error handling
  if(anyerr.ne.0) then
     if(keyerr(1).eq.668) then
        continue
!        GRHydro_C2P_failed = 3.0d0
!        write(warnline,"(A3,3i3,1P10E15.6)") "*1", ii,jj,kk,rho,epsilon,temp,ye,xpress
!        call CCTK_WARN(1,warnline)
     else 
        !$OMP CRITICAL
        call CCTK_WARN(failinfomode,"EOS error in c2p 0")
        write(warnline,"(4i5,1P10E15.6)") GRHydro_Reflevel,ii,jj,kk,x,y,z
        call CCTK_WARN(failinfomode,warnline)
        write(warnline,"(1P10E15.6)") rho,dens,epsilon,temp,temp0,ye
        call CCTK_WARN(failinfomode,warnline)
        write(warnline,"(A7,i8)") "code: ",keyerr(1)
        call CCTK_WARN(failinfomode,warnline)
        write(warnline,"(A10,i5)") "reflevel: ", GRHydro_reflevel
        call CCTK_WARN(failinfomode,warnline)
        call CCTK_WARN(failwarnmode,"Aborting!!!")
        !$OMP END CRITICAL
     endif
  endif

!!$  Check that the variables have a chance of being physical

  if( (utau + pold + udens)**2 - s2 .le. 0.0d0) then

    if (c2p_reset_pressure .ne. 0) then
      pold = sqrt(s2 + c2p_reset_pressure_to_value) - utau - udens
    else 
      !$OMP CRITICAL
!!!!      call CCTK_WARN(GRHydro_NaN_verbose, "c2p failed and being told not to reset the pressure")
      GRHydro_C2P_failed = 1.0d0
      !$OMP END CRITICAL
    endif
    
  endif
  
!!$  Calculate rho and epsilon 

!define temporary variables to speed up

  rho = udens * sqrt( (utau + pold + udens)**2 - s2)/(utau + pold + udens)
  w_lorentz = (utau + pold + udens) / sqrt( (utau + pold + udens)**2 - s2)
  epsilon = (sqrt( (utau + pold + udens)**2 - s2) - pold * w_lorentz - &
       udens)/udens
  
!!$  Calculate the function

  call EOS_Omni_press(handle,keytemp,GRHydro_eos_rf_prec,n,&
                      rho,epsilon,temp,ye,xpress,keyerr,anyerr)
  ! error handling
  if(anyerr.ne.0) then
     if(keyerr(1).eq.668) then
        continue
!        GRHydro_C2P_failed = 3.0d0
!        write(warnline,"(A3,3i3,1P10E15.6)") "*2", ii,jj,kk,rho,epsilon,temp,ye,xpress
!        call CCTK_WARN(1,warnline)
     else 
        !$OMP CRITICAL
        call CCTK_WARN(failinfomode,"EOS error in c2p 1")
        write(warnline,"(4i5,1P10E15.6)") GRHydro_Reflevel,ii,jj,kk,x,y,z
        call CCTK_WARN(failinfomode,warnline)
        write(warnline,"(1P10E15.6)") rho,dens,epsilon,temp,temp0,ye
        call CCTK_WARN(failinfomode,warnline)
        write(warnline,"(A7,i8)") "code: ",keyerr(1)
        call CCTK_WARN(failinfomode,warnline)
        write(warnline,"(A10,i5)") "reflevel: ", GRHydro_reflevel
        call CCTK_WARN(failinfomode,warnline)
        call CCTK_WARN(failwarnmode,"Aborting!!!")
        !$OMP END CRITICAL
     endif
  endif

  f = pold - xpress

  if (f .ne. f) then 
     ! Ok, this yielded nonsense, let's enforce bisection!
     mustbisect = .true.
  endif

!!$Find the root
  
  count = 0
  pnew = pold
  do while ( ((abs(pnew - pold)/abs(pnew) .gt. local_perc_ptol) .and. &
       (abs(pnew - pold) .gt. GRHydro_del_ptol))  .or. &
       (count .lt. GRHydro_countmin))
    count = count + 1

#if 0
    if(myproc.eq.3.and.&
         cctk_iteration.eq.11929.and.ii.eq.64.and.jj.eq.19.and.kk.eq.29) then
       !$OMP CRITICAL
       write(warnline,"(i5,1P10E18.9)") count, (abs(pnew - pold)/abs(pnew)), epsilon, temp, press
       call CCTK_WARN(1,warnline)
       !$OMP END CRITICAL
    endif
#endif

    if (count > GRHydro_countmax) then
      ! non-convergence is now handled outside of this
      ! routine
      GRHydro_C2P_failed = 2.0d0
      temp = temp0
      return
    end if

    call EOS_Omni_dpderho_dpdrhoe(handle,keytemp,GRHydro_eos_rf_prec,&
         n,rho,epsilon,temp,ye,dpressbydeps,dpressbydrho,keyerr,anyerr)

    temp1 = (utau+udens+pnew)**2 - s2
    drhobydpress = udens * s2 / (sqrt(temp1)*(udens+utau+pnew)**2)
    depsbydpress = pnew * s2 / (rho * (udens + utau + pnew) * temp1)
    df = 1.0d0 - dpressbydrho*drhobydpress - &
         dpressbydeps*depsbydpress

    pold = pnew
    
    ! Try to obtain new pressure via Newton-Raphson.
    
    pnew = pold - f/df 
    
    ! Check if Newton-Raphson resulted in something reasonable!
    if (c2p_resort_to_bisection.ne.0) then 
    
      tmp = (utau + pnew + udens)**2 - s2
      plow = max(pminl, sqrt(s2) - utau - udens)
      
      if (pnew .lt. plow .or. tmp .le. 0.0d0 .or. mustbisect) then
      
         ! Ok, Newton-Raphson ended up finding something unphysical.
         ! Let's try to find our root via bisection (which converges slower but is more robust)
      
         pnew = (plow + pold) / 2
         tmp = (utau + pnew + udens)**2 - s2
         
         mustbisect = .false.
      end if
    
    else
      
      ! This is the standard algorithm without resorting to bisection.
    
      pnew = max(pminl, pnew)
      tmp = (utau + pnew + udens)**2 - s2
    
    endif
    
    ! Check if we are still in the physical range now.
    ! If not, we set the C2P failed mask, and set pnew = pold (which will exit the loop).
    
    if ((tmp .le. 0.0d0) .and. GRHydro_C2P_failed .eq. 0) then
      GRHydro_C2P_failed = 1.0d0
      pnew = pold
    endif
    
!!$    Recalculate primitive variables and function
       
    rho = udens * sqrt( (utau + pnew + udens)**2 - s2)/(utau + pnew + udens)
    w_lorentz = (utau + pnew + udens) / sqrt( (utau + pnew + udens)**2 - &
         s2)
    epsilon = (sqrt( (utau + pnew + udens)**2 - s2) - pnew * w_lorentz - &
         udens)/udens

    call EOS_Omni_press(handle,keytemp,GRHydro_eos_rf_prec,n,&
         rho,epsilon,temp,ye,xpress,keyerr,anyerr)
    ! error handling
    if(anyerr.ne.0) then
       if(keyerr(1).eq.668) then
          continue
!          GRHydro_C2P_failed = 3.0d0
!          write(warnline,"(A3,3i3,1P10E15.6)") "*3", ii,jj,kk,rho,epsilon,temp,ye,xpress
!          call CCTK_WARN(1,warnline)
       else 
          !$OMP CRITICAL
          call CCTK_WARN(failinfomode,"EOS error in c2p 2")
          write(warnline,"(4i5,1P10E15.6)") GRHydro_Reflevel,ii,jj,kk,x,y,z
          call CCTK_WARN(failinfomode,warnline)
          write(warnline,"(1P10E15.6)") rho,dens,epsilon,temp,temp0,ye
          call CCTK_WARN(failinfomode,warnline)
          write(warnline,"(A7,i8)") "code: ",keyerr(1)
          call CCTK_WARN(failinfomode,warnline)
          write(warnline,"(A10,i5)") "reflevel: ", GRHydro_reflevel
          call CCTK_WARN(failinfomode,warnline)
          call CCTK_WARN(failwarnmode,"Aborting!!!")
          !$OMP END CRITICAL
       endif
    endif

    f = pnew - xpress

    if (f .ne. f) then 
      ! Ok, this yielded nonsense, let's enforce bisection!
      mustbisect = .true.
    endif

  enddo

  
!!$  Polish the root

  do i=1,GRHydro_polish

    call EOS_Omni_DPressByDRho(handle,keytemp,GRHydro_eos_rf_prec,n,&
         rho,epsilon,temp,ye,dpressbydrho,keyerr,anyerr)

    call EOS_Omni_DPressByDEps(handle,keytemp,GRHydro_eos_rf_prec,n,&
         rho,epsilon,temp,ye,dpressbydeps,keyerr,anyerr)

    temp1 = (utau+udens+pnew)**2 - s2
    drhobydpress = udens * s2 / (sqrt(temp1)*(udens+utau+pnew)**2)
    depsbydpress = pnew * s2 / (rho * (udens + utau + pnew) * temp1)
    df = 1.0d0 - dpressbydrho*drhobydpress - &
         dpressbydeps*depsbydpress
    pold = pnew
    pnew = pold - f/df
    
!!$    Recalculate primitive variables and function

    rho = udens * sqrt( (utau + pnew + udens)**2 - s2)/(utau + pnew + udens)
    w_lorentz = (utau + pnew + udens) / sqrt( (utau + pnew + udens)**2 - &
         s2)
    epsilon = (sqrt( (utau + pnew + udens)**2 - s2) - pnew * w_lorentz - &
         udens)/udens

    call EOS_Omni_press(handle,keytemp,GRHydro_eos_rf_prec,n, &
         rho,epsilon,temp,ye,xpress,keyerr,anyerr)

    ! error handling
    if(anyerr.ne.0) then
       if(keyerr(1).eq.668) then
          continue
!          GRHydro_C2P_failed = 3.0d0
!          write(warnline,"(A3,3i3,1P10E15.6)") "*4", ii,jj,kk,rho,epsilon,temp,ye,xpress
!          call CCTK_WARN(1,warnline)
       else 
          !$OMP CRITICAL
          call CCTK_WARN(failinfomode,"EOS error in c2p 3")
          write(warnline,"(4i5,1P10E15.6)") GRHydro_Reflevel,ii,jj,kk,x,y,z
          call CCTK_WARN(failinfomode,warnline)
          write(warnline,"(1P10E15.6)") rho,dens,epsilon,temp,temp0,ye
          call CCTK_WARN(failinfomode,warnline)
          write(warnline,"(A7,i8)") "code: ",keyerr(1)
          call CCTK_WARN(failinfomode,warnline)
          write(warnline,"(A10,i5)") "reflevel: ", GRHydro_reflevel
          call CCTK_WARN(failinfomode,warnline)
          call CCTK_WARN(failwarnmode,"Aborting!!!")
          !$OMP END CRITICAL
       endif
    endif
    
    f = pold - xpress

  enddo

!!$  Calculate primitive variables from root

  !if (rho .le. GRHydro_rho_min*(1.d0+GRHydro_atmo_tolerance) ) then
  IF_BELOW_ATMO(rho, GRHydro_rho_min, GRHydro_atmo_tolerance, r) then
    SET_ATMO_MIN(rho, GRHydro_rho_min, r) !GRHydro_rho_min
    udens = rho
    dens = sdet * rho
    temp = GRHydro_hot_atmo_temp
    ye = GRHydro_hot_atmo_Y_e
    ye_con = dens * ye
    keytemp=1
    call EOS_Omni_press(handle,keytemp,GRHydro_eos_rf_prec,n, &
         rho,epsilon,temp,ye,xpress,keyerr,anyerr)
    keytemp=0
    ! w_lorentz=1, so the expression for utau reduces to:
    utau  = rho + rho*epsilon - udens
    sx = 0.d0
    sy = 0.d0
    sz = 0.d0
    s2 = 0.d0
    usx = 0.d0
    usy = 0.d0
    usz = 0.d0
    w_lorentz = 1.d0
  end if

  press = pnew
  vlowx = usx / ( (rho + rho*epsilon + press) * w_lorentz**2)
  vlowy = usy / ( (rho + rho*epsilon + press) * w_lorentz**2)
  vlowz = usz / ( (rho + rho*epsilon + press) * w_lorentz**2)
  velx = uxx * vlowx + uxy * vlowy + uxz * vlowz
  vely = uxy * vlowx + uyy * vlowy + uyz * vlowz
  velz = uxz * vlowx + uyz * vlowy + uzz * vlowz

  ! indicate to wrapper routine that we dropped out of the
  ! EOS table
  if(keyerr(1).eq.668) then
     GRHydro_C2P_failed = 3.0d0
  endif

end subroutine Con2Prim_pt_hot3

