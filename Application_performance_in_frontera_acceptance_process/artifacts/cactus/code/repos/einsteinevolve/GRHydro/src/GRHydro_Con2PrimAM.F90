/*@@
   @file      GRHydro_Con2PrimAM.F90
   @date      Sep 3, 2010
   @author    Scott Noble, Joshua Faber, Bruno Mundim, Tanja Bode
   @desc 
   The routines for converting conservative to primitive variables.
   Vector-potential MHD.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "SpaceMask.h"
#include "GRHydro_Macros.h"

#define ITER_TOL (1.0e-8)
#define MAXITER (50)

 /*@@
   @routine    Conservative2PrimitiveAM
   @date       Sep 3, 2010
   @author     Scott Noble, Joshua Faber, Bruno Mundim, Ian Hawke
   @desc 
   Wrapper routine that converts from conserved to primitive variables
   at every grid cell centre.
   @enddesc 
   @calls     
   @calledby   
   @history 
   Trimmed and altered from the GR3D routines, original author Mark Miller.
   2007?: Bruno excluded the points in the atmosphere and excision region from the computation.
   Aug. 2008: Luca added a check on whether a failure at a given point may be disregarded, 
   because that point will then be restricted from a finer level. This should be completely 
   safe only if *regridding happens at times when all levels are evolved.*
   Feb. 2009: The above procedure proved to be wrong, so Luca implemented another one. 
   When a failure occurs, it is temporarily ignored, except for storing the location of where 
   it occured in a mask. At the end, after a Carpet restriction, the mask is checked and if 
   it still contains failures, the run is aborted with an error message. Only used routines 
   have been updated to use this procedure.
   @endhistory 

@@*/

subroutine Conservative2PrimitiveAM(CCTK_ARGUMENTS)

  use Con2PrimM_fortran_interfaces

  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  integer :: i, j, k, itracer, nx, ny, nz
  CCTK_REAL :: uxx, uxy, uxz, uyy, uyz, uzz, det, sdet, pmin(1), epsmin(1)
  CCTK_REAL :: oob, b2, d2, s2, bscon, bxhat, byhat, bzhat, bhatscon
  CCTK_REAL :: Wm, Wm0, Wm_resid, Wmold
  CCTK_REAL :: s2m, s2m0, s2m_resid, s2mold, s2max, taum
  CCTK_INT :: niter
  CCTK_INT :: epsnegative
  character(len=100) warnline
  
  CCTK_REAL :: local_min_tracer, local_gam(1), local_pgam,local_K, sc
  CCTK_REAL :: local_perc_ptol

! begin EOS Omni vars                                                                                       
  CCTK_INT  :: n,keytemp,anyerr,keyerr(1)
  CCTK_REAL :: xpress(1),xtemp(1),xye(1),xeps(1),xrho(1),one(1)=1.0d0
! end EOS Omni vars                 

  CCTK_INT :: GRHydro_UseGeneralCoordinates
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3) :: g11, g12, g13, g22, g23, g33
  pointer (pg11,g11), (pg12,g12), (pg13,g13), (pg22,g22), (pg23,g23), (pg33,g33)
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3,3) :: vup, Bprim
  pointer (pvup,vup), (pBprim,Bprim)

  logical :: posdef

  CCTK_REAL :: g11c, g12c, g13c, g22c, g23c, g33c
  CCTK_REAL :: tmp1

  ! Save the primitive variables to temporary functions before calling the
  ! con2prim pointwise routines:
  CCTK_REAL :: rho_tmp, press_tmp, eps_tmp
  CCTK_REAL :: velx_tmp, vely_tmp, velz_tmp, w_lorentz_tmp
  CCTK_REAL :: Bvecx_tmp, Bvecy_tmp, Bvecz_tmp
  CCTK_REAL :: Bconsx_tmp, Bconsy_tmp, Bconsz_tmp

  CCTK_REAL :: bdotv, magpress

  ! Assume 3-metric is positive definite. Check deep inside the horizon 
  ! if this is actually satisfied and if it is not then cast the metric 
  !as conformally flat only for con2prim inversion purposes.
  posdef = .true. 
 
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
#define gxx faulty_gxx
#define gxy faulty_gxy
#define gxz faulty_gxz
#define gyy faulty_gyy
#define gyz faulty_gyz
#define gzz faulty_gzz
#define betax faulty_betax
#define betay faulty_betay
#define betaz faulty_betaz
#define vel faulty_vel
#define Bvec faulty_Bvec

! begin EOS Omni vars                                                                                       
  n=1;keytemp=0;anyerr=0;keyerr(1)=0
  xpress(1)=0.0d0;xtemp(1)=0.0d0;xye(1)=0.0d0;xeps(1)=0.0d0
! end EOS Omni vars                 

  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)
  
  if (use_min_tracer .ne. 0) then
     local_min_tracer = min_tracer
  else
     local_min_tracer = 0d0
  end if

!  call EOS_Omni_press(GRHydro_polytrope_handle,keytemp,GRHydro_eos_rf_prec,n,&
!         GRHydro_rho_min,xeps,xtemp,xye,pmin,keyerr,anyerr)
!  call EOS_Omni_EpsFromPress(GRHydro_polytrope_handle,keytemp,GRHydro_eos_rf_prec,n,&
!       GRHydro_rho_min,xeps,xtemp,xye,pmin,epsmin,keyerr,anyerr)
  ! this is a poly call
  xrho(1)=GRHydro_rho_min
  call EOS_Omni_press(GRHydro_polytrope_handle,keytemp,GRHydro_eos_rf_prec,n,&
       xrho,xeps,xtemp,xye,pmin,keyerr,anyerr)
  call EOS_Omni_EpsFromPress(GRHydro_polytrope_handle,keytemp,GRHydro_eos_rf_prec,n,&
       xrho,xeps,xtemp,xye,pmin,epsmin,keyerr,anyerr)

  call EOS_Omni_press(GRHydro_eos_handle,keytemp,GRHydro_eos_rf_prec,n,&
         one,one,xtemp,xye,local_gam,keyerr,anyerr)
  local_gam = local_gam+1.d0

  !call CCTK_WARN(1,"In Con2PrimAM")

  !$OMP PARALLEL DO PRIVATE(i,j,k,itracer,&
  !$OMP uxx, uxy, uxz, uyy, uyz, uzz, det, epsnegative, &
  !$OMP b2,xrho,xeps,xpress,xtemp,local_K,local_pgam,sc,keyerr,anyerr,keytemp, &
  !$OMP local_perc_ptol,posdef,g11c,g12c,g13c,g22c,g23c,g33c,tmp1, &
  !$OMP sdet,d2,s2,oob,bscon,bxhat,byhat,bzhat, &
  !$OMP bhatscon,Wm,Wm0,Wm_resid,Wmold,s2m,s2m0,s2m_resid,s2mold,s2max, &
  !$OMP taum,niter,rho_tmp,press_tmp,eps_tmp,velx_tmp,vely_tmp,velz_tmp, &
  !$OMP Bconsx_tmp, Bconsy_tmp, Bconsz_tmp, &
  !$OMP w_lorentz_tmp,Bvecx_tmp,Bvecy_tmp,Bvecz_tmp,bdotv,magpress)
  do k = 1, nz 
     do j = 1, ny 
        do i = 1, nx
           
           !do not compute if in atmosphere or in excised region
           if ((atmosphere_mask(i,j,k) .ne. 0) .or. &
                (hydro_excision_mask(i,j,k) .ne. 0)) cycle
           
           epsnegative = 0
           
           sdet = sdetg(i,j,k)

           call UpperMetric(uxx,uxy,uxz,uyy,uyz,uzz,sdet*sdet,&
                g11(i,j,k),g12(i,j,k),g13(i,j,k),g22(i,j,k),&
                g23(i,j,k),g33(i,j,k))        
           
!!$ Tracers don't need an MHD treatment!
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

           
           b2=g11(i,j,k)*Bprim(i,j,k,1)**2+g22(i,j,k)*Bprim(i,j,k,2)**2+g33(i,j,k)*Bprim(i,j,k,3)**2+ &
           2.0*(g12(i,j,k)*Bprim(i,j,k,1)*Bprim(i,j,k,2)+g13(i,j,k)*Bprim(i,j,k,1)*Bprim(i,j,k,3)+ &
           g23(i,j,k)*Bprim(i,j,k,2)*Bprim(i,j,k,3))

       
           if ( dens(i,j,k) .le. sdet*GRHydro_rho_min*(1.d0+GRHydro_atmo_tolerance) ) then
              
              !call CCTK_WARN(1,"Con2Prim: Resetting to atmosphere")
              !write(warnline,"(3i5,1P10E15.6)") i,j,k,x(i,j,k),y(i,j,k),z(i,j,k)
              !call CCTK_WARN(1,warnline)
              !write(warnline,"(1P10E15.6)") rho(i,j,k),dens(i,j,k),eps(i,j,k),&
              !      temperature(i,j,k),y_e(i,j,k)
              !call CCTK_WARN(1,warnline)

              
              dens(i,j,k) = sdet*GRHydro_rho_min !/(1.d0+GRHydro_atmo_tolerance)
              rho(i,j,k) = GRHydro_rho_min
              scon(i,j,k,:) = 0.d0
              vup(i,j,k,:) = 0.d0
              w_lorentz(i,j,k) = 1.d0
              
              if(evolve_temper.ne.0) then
              ! set hot atmosphere values
                temperature(i,j,k) = grhydro_hot_atmo_temp
                y_e(i,j,k) = grhydro_hot_atmo_Y_e
                y_e_con(i,j,k) = y_e(i,j,k) * dens(i,j,k)
                keytemp = 1
                call EOS_Omni_press(GRHydro_eos_handle,keytemp,GRHydro_eos_rf_prec,n,&
                   rho(i,j,k),eps(i,j,k),temperature(i,j,k),y_e(i,j,k),&
                   press(i,j,k),keyerr,anyerr)
              else
                keytemp = 0
                call EOS_Omni_press(GRHydro_polytrope_handle,keytemp,GRHydro_eos_rf_prec,n,&
                   rho(i,j,k),eps(i,j,k),xtemp,xye,press(i,j,k),keyerr,anyerr)
                call EOS_Omni_EpsFromPress(GRHydro_polytrope_handle,keytemp,GRHydro_eos_rf_prec,n,&
                   rho(i,j,k),eps(i,j,k),xtemp,xye,press(i,j,k),eps(i,j,k),keyerr,anyerr)
              endif

              !call EOS_Omni_press(GRHydro_polytrope_handle,keytemp,GRHydro_eos_rf_prec,n,&
              !     rho(i,j,k),eps(i,j,k),xtemp,xye,press(i,j,k),keyerr,anyerr)
              ! 
              !call EOS_Omni_EpsFromPress(GRHydro_polytrope_handle,keytemp,GRHydro_eos_rf_prec,n,&
              !     rho(i,j,k),eps(i,j,k),xtemp,xye,press(i,j,k),eps(i,j,k),keyerr,anyerr)

              ! w_lorentz=1, so the expression for tau reduces to:
              
              !!$ tau does need to take into account the existing B-field
              !!$ with w_lorentz=1, we find tau = sqrtdet*(rho (1+eps+b^2/2)) - dens  [Press drops out]
              tau(i,j,k)  = sdet * (rho(i,j,k)*(1.0+eps(i,j,k)+b2/2.0)) - dens(i,j,k)
              
              if(tau(i,j,k).le.sdet*b2*0.5d0)then
                tau(i,j,k) = GRHydro_tau_min + sdet*b2*0.5d0
              endif

              cycle
              
           end if


           if(evolve_temper.eq.0) then
            

              if(sqrtdet_thr.gt.0d0 .and. sdet.ge.sqrtdet_thr) then
                    d2 = dens(i,j,k)**2
                    s2 = uxx*scon(i,j,k,1)**2 + uyy*scon(i,j,k,2)**2 &
                                              + uzz*scon(i,j,k,3)**2 &
                       + 2.0d0*uxy*scon(i,j,k,1)*scon(i,j,k,2) &
                       + 2.0d0*uxz*scon(i,j,k,1)*scon(i,j,k,3) &
                       + 2.0d0*uyz*scon(i,j,k,2)*scon(i,j,k,3) 
                   oob = 1.0d0/sqrt(b2)
                 bxhat = oob*Bprim(i,j,k,1)
                 byhat = oob*Bprim(i,j,k,2)
                 bzhat = oob*Bprim(i,j,k,3)
                 bhatscon = bxhat*scon(i,j,k,1)+byhat*scon(i,j,k,2) &
                                               +bzhat*scon(i,j,k,3)
                 bscon = Bprim(i,j,k,1)*scon(i,j,k,1) &
                       + Bprim(i,j,k,2)*scon(i,j,k,2) &
                       + Bprim(i,j,k,3)*scon(i,j,k,3)
                 ! Initial guesses for iterative procedure to find Wm:
                 Wm0 = sdet*sqrt(bhatscon**2+d2)
                 s2m0 = (Wm0**2*s2+bhatscon**2*(b2+2.0d0*Wm0)) &
                      / (Wm0+b2)**2 
                 Wm = sdet*sqrt(s2m0+d2)
                 s2m = (Wm**2*s2+bscon**2*(b2+2.0d0*Wm)) &
                     / (Wm+b2)**2 
                 s2m_resid = 1.0d60
                 Wm_resid = 1.0d60
                 niter = 0
                 do while((s2m_resid.ge.ITER_TOL.and.Wm_resid.ge.ITER_TOL).and.&
                          niter.le.MAXITER)
                   Wmold = Wm
                   s2mold = s2m
                   Wm = sdet*sqrt(s2m+d2)
                   s2m = (Wm**2*s2+bscon**2*(b2+2.0d0*Wm)) &
                       / (Wm+b2)**2 
                   Wm_resid = abs(Wmold-Wm)
                   s2m_resid = abs(s2mold-s2m)
                   niter = niter + 1
                 end do
                 !TODO: abort execution if niter .eq. MAXITER and warn user
                 taum = tau(i,j,k) - 0.5d0*sdet*b2 -0.5d0*(b2*s2-bscon**2)/ &
                                                         (sdet*(Wm+b2)**2) 
                 s2max = taum*(taum+2.0d0*dens(i,j,k))
                 if(taum.lt.GRHydro_tau_min)then
                   tau(i,j,k) = GRHydro_tau_min + 0.5d0*sdet*b2 + 0.5d0* &
                                (b2*s2-bscon**2)/(sdet*(Wm+b2)**2)
                 end if
                 if(s2.gt.s2max) then
                   scon(i,j,k,1) = scon(i,j,k,1)*sqrt(s2max/s2)
                   scon(i,j,k,2) = scon(i,j,k,2)*sqrt(s2max/s2)
                   scon(i,j,k,3) = scon(i,j,k,3)*sqrt(s2max/s2)
                 end if
              endif 

              rho_tmp = rho(i,j,k)
              press_tmp = press(i,j,k)
              eps_tmp = eps(i,j,k)
              velx_tmp = vup(i,j,k,1)
              vely_tmp = vup(i,j,k,2)
              velz_tmp = vup(i,j,k,3)
              w_lorentz_tmp = w_lorentz(i,j,k)
              Bvecx_tmp = Bprim(i,j,k,1)
              Bvecy_tmp = Bprim(i,j,k,2)
              Bvecz_tmp = Bprim(i,j,k,3)

              !! We've already create B from A, construct Bcons
              Bconsx_tmp = sdet*Bvecx_tmp 
              Bconsy_tmp = sdet*Bvecy_tmp 
              Bconsz_tmp = sdet*Bvecz_tmp 

              keytemp = 0
              !Watch out for the values returned to b2. Here b2 is the Bprim^2
              !while inside the point-wise con2prim routines it is the square
              !of the comoving B-field, b^{\mu} b_{\mu}. It is overwritten 
              !in this routine, but we may need to find a better notation 
              !avoid future confusions.
             if(GRHydro_eos_handle .eq. 1 .or. GRHydro_eos_handle .eq. 2) then
             
                call GRHydro_Con2PrimM_ptold(GRHydro_eos_handle, local_gam(1), dens(i,j,k), &
                   scon(i,j,k,1),scon(i,j,k,2),scon(i,j,k,3), tau(i,j,k), &
                   Bconsx_tmp,Bconsy_tmp,Bconsz_tmp,rho_tmp, & 
                   velx_tmp,vely_tmp,velz_tmp,eps_tmp,press_tmp, &
                   Bvecx_tmp,Bvecy_tmp,Bvecz_tmp,b2, w_lorentz_tmp,&
                   g11(i,j,k),g12(i,j,k),g13(i,j,k),g22(i,j,k),g23(i,j,k),g33(i,j,k), &
                   uxx,uxy,uxz,uyy,uyz,uzz,sdet, &
                   epsnegative,GRHydro_C2P_failed(i,j,k))
                   
             else
              call GRHydro_Con2PrimM_pt2(GRHydro_eos_handle, keytemp, GRHydro_eos_rf_prec, GRHydro_perc_ptol, local_gam(1), dens(i,j,k), &
                   scon(i,j,k,1),scon(i,j,k,2),scon(i,j,k,3), tau(i,j,k), &
                   Bconsx_tmp, Bconsy_tmp, Bconsz_tmp,xye(1), &
                   xtemp(1),rho_tmp,velx_tmp,vely_tmp,velz_tmp,&
                   eps_tmp,press_tmp,Bvecx_tmp,Bvecy_tmp,Bvecz_tmp,b2,&
                   w_lorentz_tmp,g11(i,j,k),g12(i,j,k),g13(i,j,k),&
                   g22(i,j,k),g23(i,j,k),g33(i,j,k), &
                   uxx,uxy,uxz,uyy,uyz,uzz,sdet, &
                   epsnegative,GRHydro_C2P_failed(i,j,k))
              endif

              if(evolve_entropy.ne.0) then
                if(GRHydro_C2P_failed(i,j,k).ne.0) then
                  !Use previous time step for rho:
                 entropy(i,j,k) = entropycons(i,j,k)/dens(i,j,k)*rho(i,j,k)
                else
                  !Use the current correct value of rho returned by con2prim:
                  entropy(i,j,k) = entropycons(i,j,k)/dens(i,j,k)*rho_tmp
                endif
              endif

              if(GRHydro_C2P_failed(i,j,k).ne.0) then
                xrho=1.0d0; xtemp=0.0d0; xeps=1.0d0
                call EOS_Omni_press(GRHydro_polytrope_handle,keytemp,GRHydro_eos_rf_prec,n,&
                     xrho,xeps,xtemp,xye,xpress,keyerr,anyerr)
                local_K = xpress(1); 

                xrho=10.0d0; xeps=1.0d0
                call EOS_Omni_press(GRHydro_polytrope_handle,keytemp,GRHydro_eos_rf_prec,n,&
                     xrho,xeps,xtemp,xye,xpress,keyerr,anyerr)
                local_pgam=log(xpress(1)/local_K)/log(xrho(1))
                sc = local_K*dens(i,j,k)

                if(sqrtdet_thr.gt.0d0 .and. sdet.ge.sqrtdet_thr) then
                  GRHydro_C2P_failed(i,j,k) = 0

                  rho_tmp = rho(i,j,k)
                  press_tmp = press(i,j,k)
                  eps_tmp = eps(i,j,k)
                  velx_tmp = vup(i,j,k,1)
                  vely_tmp = vup(i,j,k,2)
                  velz_tmp = vup(i,j,k,3)
                  w_lorentz_tmp = w_lorentz(i,j,k)
                  Bvecx_tmp = Bprim(i,j,k,1)
                  Bvecy_tmp = Bprim(i,j,k,2)
                  Bvecz_tmp = Bprim(i,j,k,3)
                  Bconsx_tmp = sdet*Bvecx_tmp 
                  Bconsy_tmp = sdet*Bvecy_tmp 
                  Bconsz_tmp = sdet*Bvecz_tmp 

                  call GRHydro_Con2PrimM_Polytype_pt(GRHydro_eos_handle, local_pgam, &
                       dens(i,j,k),scon(i,j,k,1),scon(i,j,k,2),scon(i,j,k,3), sc, &
                       Bconsx_tmp, Bconsy_tmp, Bconsz_tmp,rho_tmp,&
                       velx_tmp,vely_tmp,velz_tmp,eps_tmp,press_tmp,&
                       Bvecx_tmp,Bvecy_tmp,Bvecz_tmp,b2,w_lorentz_tmp,&
                       g11(i,j,k),g12(i,j,k),g13(i,j,k),g22(i,j,k),g23(i,j,k),g33(i,j,k), &
                       uxx,uxy,uxz,uyy,uyz,uzz,det, &
                       epsnegative,GRHydro_C2P_failed(i,j,k))

                  if(GRHydro_C2P_failed(i,j,k).ne.0) then
                    GRHydro_C2P_failed(i,j,k) = 0
                    call prim2conAM(GRHydro_eos_handle,g11(i,j,k),g12(i,j,k), &
                         g13(i,j,k),g22(i,j,k),g23(i,j,k),g33(i,j,k),det, &
                       dens(i,j,k),scon(i,j,k,1),scon(i,j,k,2),scon(i,j,k,3), &
                      tau(i,j,k), &
                          rho(i,j,k),vup(i,j,k,1),vup(i,j,k,2),vup(i,j,k,3), &
                          eps(i,j,k),press(i,j,k), &
                  Bprim(i,j,k,1),Bprim(i,j,k,2),Bprim(i,j,k,3),w_lorentz(i,j,k))
                    cycle
                  end if
                end if

                bdotv=g11(i,j,k)*Bprim(i,j,k,1)*vup(i,j,k,1)+ &
                      g22(i,j,k)*Bprim(i,j,k,2)*vup(i,j,k,2)+ &
                      g33(i,j,k)*Bprim(i,j,k,3)*vup(i,j,k,3)+ &
                 2.0*(g12(i,j,k)*Bprim(i,j,k,1)*vup(i,j,k,2)+ &
                      g13(i,j,k)*Bprim(i,j,k,1)*vup(i,j,k,3)+ &
                      g23(i,j,k)*Bprim(i,j,k,2)*vup(i,j,k,3))

                magpress = 0.5d0*(b2/w_lorentz(i,j,k)**2+bdotv**2)

                if(rho(i,j,k)*eps(i,j,k)*max_magnetic_to_gas_pressure_ratio.le.magpress) then
                  GRHydro_C2P_failed(i,j,k) = 0

                  if(evolve_entropy.ne.0) then
                    local_K = entropycons(i,j,k)/dens(i,j,k) 
                    xrho=10.0d0; xeps=1.0d0
                    call EOS_Omni_press(GRHydro_polytrope_handle,keytemp,GRHydro_eos_rf_prec,n,&
                         xrho,xeps,xtemp,xye,xpress,keyerr,anyerr)
                    local_pgam=log(xpress(1)/local_K)/log(xrho(1))
                    sc = local_K*dens(i,j,k)
                  end if

                  rho_tmp = rho(i,j,k)
                  press_tmp = press(i,j,k)
                  eps_tmp = eps(i,j,k)
                  velx_tmp = vup(i,j,k,1)
                  vely_tmp = vup(i,j,k,2)
                  velz_tmp = vup(i,j,k,3)
                  w_lorentz_tmp = w_lorentz(i,j,k)
                  Bvecx_tmp = Bprim(i,j,k,1)
                  Bvecy_tmp = Bprim(i,j,k,2)
                  Bvecz_tmp = Bprim(i,j,k,3)
                  Bconsx_tmp = sdet*Bvecx_tmp
                  Bconsy_tmp = sdet*Bvecy_tmp
                  Bconsz_tmp = sdet*Bvecz_tmp

                  if(evolve_entropy.ne.0) then
                    call GRHydro_Con2PrimM_ptee(GRHydro_eos_handle, keytemp, &
                         GRHydro_eos_rf_prec, local_gam(1), dens(i,j,k), &
                         scon(i,j,k,1),scon(i,j,k,2),scon(i,j,k,3), tau(i,j,k), &
                         Bconsx_tmp,Bconsy_tmp,Bconsz_tmp, &
                         entropycons(i,j,k), xye(1), &
                         xtemp(1),rho_tmp,velx_tmp,vely_tmp,velz_tmp,&
                         eps_tmp,press_tmp,Bvecx_tmp,Bvecy_tmp,Bvecz_tmp,b2,&
                         w_lorentz_tmp,g11(i,j,k),g12(i,j,k),g13(i,j,k),&
                         g22(i,j,k),g23(i,j,k),g33(i,j,k), &
                         uxx,uxy,uxz,uyy,uyz,uzz,sdet, &
                         epsnegative,GRHydro_C2P_failed(i,j,k))
                  else
                  call GRHydro_Con2PrimM_Polytype_pt(GRHydro_eos_handle, local_pgam, &
                       dens(i,j,k),scon(i,j,k,1),scon(i,j,k,2),scon(i,j,k,3), sc, &
                       Bconsx_tmp, Bconsy_tmp, Bconsz_tmp, rho_tmp,&
                       velx_tmp,vely_tmp,velz_tmp,eps_tmp,press_tmp,&
                       Bvecx_tmp,Bvecy_tmp,Bvecz_tmp,b2,w_lorentz_tmp,&
                       g11(i,j,k),g12(i,j,k),g13(i,j,k),g22(i,j,k),g23(i,j,k),g33(i,j,k), &
                       uxx,uxy,uxz,uyy,uyz,uzz,sdet, &
                       epsnegative,GRHydro_C2P_failed(i,j,k))
                  end if 

                  rho(i,j,k) = rho_tmp 
                  press(i,j,k) = press_tmp 
                  eps(i,j,k) = eps_tmp 
                  vup(i,j,k,1) = velx_tmp 
                  vup(i,j,k,2) = vely_tmp 
                  vup(i,j,k,3) = velz_tmp 
                  w_lorentz(i,j,k) = w_lorentz_tmp 
                  Bprim(i,j,k,1) = Bvecx_tmp 
                  Bprim(i,j,k,2) = Bvecy_tmp 
                  Bprim(i,j,k,3) = Bvecz_tmp 
                  cycle
                end if
              end if

           else    ! if(evolve_temper.eq.0) then

              rho_tmp = rho(i,j,k)
              press_tmp = press(i,j,k)
              eps_tmp = eps(i,j,k)
              velx_tmp = vup(i,j,k,1)
              vely_tmp = vup(i,j,k,2)
              velz_tmp = vup(i,j,k,3)
              w_lorentz_tmp = w_lorentz(i,j,k)
              Bvecx_tmp = Bprim(i,j,k,1)
              Bvecy_tmp = Bprim(i,j,k,2)
              Bvecz_tmp = Bprim(i,j,k,3)
              Bconsx_tmp = sdet*Bvecx_tmp
              Bconsy_tmp = sdet*Bvecy_tmp
              Bconsz_tmp = sdet*Bvecz_tmp

              keytemp = 0 

              call GRHydro_Con2PrimM_pt(GRHydro_eos_handle, GRHydro_reflevel, int(i,ik), int(j,ik), int(k,ik), x(i,j,k), y(i,j,k), z(i,j,k), keytemp, GRHydro_eos_rf_prec, GRHydro_perc_ptol, local_gam(1), dens(i,j,k), &
                   scon(i,j,k,1),scon(i,j,k,2),scon(i,j,k,3), tau(i,j,k), &
                   Bconsx_tmp, Bconsy_tmp, Bconsz_tmp, Y_e(i,j,k), &
                   temperature(i,j,k),rho_tmp,velx_tmp,vely_tmp,velz_tmp,&
                   eps_tmp,press_tmp,Bvecx_tmp,Bvecy_tmp,Bvecz_tmp,b2,&
                   w_lorentz_tmp,g11(i,j,k),g12(i,j,k),g13(i,j,k),&
                   g22(i,j,k),g23(i,j,k),g33(i,j,k), &
                   uxx,uxy,uxz,uyy,uyz,uzz,sdet, &
                   epsnegative,GRHydro_C2P_failed(i,j,k))
              if(GRHydro_C2P_failed(i,j,k).ne.0) then
                ! this means c2p did not converge.
                ! In this case, we attempt to call c2p with a reduced
                ! accuracy requirement; if it fails again, we abort
                GRHydro_C2P_failed(i,j,k) = 0
                local_perc_ptol = GRHydro_eos_rf_prec*100.0d0
                ! Use the previous primitive values as initial guesses
                rho_tmp = rho(i,j,k)
                press_tmp = press(i,j,k)
                eps_tmp = eps(i,j,k)
                velx_tmp = vup(i,j,k,1)
                vely_tmp = vup(i,j,k,2)
                velz_tmp = vup(i,j,k,3)
                w_lorentz_tmp = w_lorentz(i,j,k)
                Bvecx_tmp = Bprim(i,j,k,1)
                Bvecy_tmp = Bprim(i,j,k,2)
                Bvecz_tmp = Bprim(i,j,k,3)
                Bconsx_tmp = sdet*Bvecx_tmp
                Bconsy_tmp = sdet*Bvecy_tmp
                Bconsz_tmp = sdet*Bvecz_tmp
              call GRHydro_Con2PrimM_pt(GRHydro_eos_handle, GRHydro_reflevel, int(i,ik), int(j,ik), int(k,ik), x(i,j,k), y(i,j,k), z(i,j,k), keytemp, GRHydro_eos_rf_prec, GRHydro_perc_ptol, local_gam(1), dens(i,j,k), &
                     scon(i,j,k,1),scon(i,j,k,2),scon(i,j,k,3), tau(i,j,k), &
                     Bconsx_tmp, Bconsy_tmp, Bconsz_tmp, Y_e(i,j,k), &
                     temperature(i,j,k),rho_tmp,velx_tmp,vely_tmp,velz_tmp,&
                     eps_tmp,press_tmp,Bvecx_tmp,Bvecy_tmp,Bvecz_tmp,b2,&
                     w_lorentz_tmp,g11(i,j,k),g12(i,j,k),g13(i,j,k),&
                     g22(i,j,k),g23(i,j,k),g33(i,j,k), &
                     uxx,uxy,uxz,uyy,uyz,uzz,sdet, &
                     epsnegative,GRHydro_C2P_failed(i,j,k))
                if(GRHydro_C2P_failed(i,j,k).ne.0) then
                  !$OMP CRITICAL
                  if (GRHydro_reflevel.ge.GRHydro_c2p_warn_from_reflevel) then
                     call CCTK_WARN(1,"Convergence problem in c2p")
                     write(warnline,"(A10,i5)") "reflevel: ",GRHydro_reflevel
                     call CCTK_WARN(1,warnline)
                     write(warnline,"(3i5,1P10E15.6)") i,j,k,x(i,j,k),y(i,j,k),z(i,j,k)
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
           endif    ! if(evolve_temper.eq.0) then

           
           if (epsnegative .ne. 0) then  
              
#if 0
              ! cott 2010/03/30:
              ! Set point to atmosphere, but continue evolution -- this is better than using
              ! the poly EOS -- it will lead the code to crash if this happens inside a (neutron) star,
              ! but will allow the job to continue if it happens in the atmosphere or in a
              ! zone that contains garbage (i.e. boundary, buffer zones)
              ! Ultimately, we want this fixed via a new carpet mask presently under development
              !           GRHydro_C2P_failed(i,j,k) = 1
              
              !$OMP CRITICAL
              call CCTK_WARN(GRHydro_NaN_verbose+2, 'Specific internal energy just went below 0! ')
              write(warnline,'(a28,i2)') 'on carpet reflevel: ',GRHydro_reflevel
              call CCTK_WARN(GRHydro_NaN_verbose+2,warnline)
              write(warnline,'(a20,3g16.7)') 'xyz location: ',&
                   x(i,j,k),y(i,j,k),z(i,j,k)
              call CCTK_WARN(GRHydro_NaN_verbose+2,warnline)
              write(warnline,'(a20,g16.7)') 'radius: ',r(i,j,k)
              call CCTK_WARN(GRHydro_NaN_verbose+2,warnline)
              call CCTK_WARN(GRHydro_NaN_verbose+2,"Setting the point to atmosphere")
              !$OMP END CRITICAL
              
              ! for safety, let's set the point to atmosphere
              dens(i,j,k) = sdet*GRHydro_rho_min !/(1.d0+GRHydro_atmo_tolerance)
              rho(i,j,k) = GRHydro_rho_min
              scon(i,j,k,:) = 0.d0
              vup(i,j,k,:) = 0.d0
              w_lorentz(i,j,k) = 1.d0

              call EOS_Omni_press(GRHydro_polytrope_handle,keytemp,GRHydro_eos_rf_prec,n,&
                   rho(i,j,k),eps(i,j,k),xtemp,xye,press(i,j,k),keyerr,anyerr)
              
              call EOS_Omni_EpsFromPress(GRHydro_polytrope_handle,keytemp,GRHydro_eos_rf_prec,n,&
                   rho(i,j,k),eps(i,j,k),xtemp,xye,press(i,j,k),eps(i,j,k),keyerr,anyerr)

              b2=g11(i,j,k)*Bprim(i,j,k,1)**2+g22(i,j,k)*Bprim(i,j,k,2)**2+g33(i,j,k)*Bprim(i,j,k,3)**2+ &
                   2.0*(g12(i,j,k)*Bprim(i,j,k,1)*Bprim(i,j,k,2)+g13(i,j,k)*Bprim(i,j,k,1)*Bprim(i,j,k,3)+ &
                   g23(i,j,k)*Bprim(i,j,k,2)*Bprim(i,j,k,3))
              
              
              ! w_lorentz=1, so the expression for tau reduces to [see above]:
              tau(i,j,k)  = sdet * (rho(i,j,k)*(1.0+eps(i,j,k)+b2/2.0)) - dens(i,j,k)           
#else
              ! cott 2010/03/27:      
              ! Honestly, this should never happen. We need to flag the point where
              ! this happened as having led to failing con2prim.
              
              !$OMP CRITICAL
              call CCTK_WARN(GRHydro_NaN_verbose+2, 'Specific internal energy just went below 0, trying polytype.')
              !$OMP END CRITICAL
              
              xrho=1.0d0; xtemp=0.0d0; xeps=1.0d0
              call EOS_Omni_press(GRHydro_polytrope_handle,keytemp,GRHydro_eos_rf_prec,n,&
                   xrho,xeps,xtemp,xye,xpress,keyerr,anyerr)
              local_K = xpress(1); 

              xrho=10.0d0; xeps=1.0d0
              call EOS_Omni_press(GRHydro_polytrope_handle,keytemp,GRHydro_eos_rf_prec,n,&
                   xrho,xeps,xtemp,xye,xpress,keyerr,anyerr)
              local_pgam=log(xpress(1)/local_K)/log(xrho(1))
              sc = local_K*dens(i,j,k)

              rho_tmp = rho(i,j,k)
              press_tmp = press(i,j,k)
              eps_tmp = eps(i,j,k)
              velx_tmp = vup(i,j,k,1)
              vely_tmp = vup(i,j,k,2)
              velz_tmp = vup(i,j,k,3)
              w_lorentz_tmp = w_lorentz(i,j,k)
              Bvecx_tmp = Bprim(i,j,k,1)
              Bvecy_tmp = Bprim(i,j,k,2)
              Bvecz_tmp = Bprim(i,j,k,3)
              Bconsx_tmp = sdet*Bvecx_tmp
              Bconsy_tmp = sdet*Bvecy_tmp
              Bconsz_tmp = sdet*Bvecz_tmp

            call GRHydro_Con2PrimM_Polytype_pt(GRHydro_eos_handle, local_pgam, dens(i,j,k), &
                   scon(i,j,k,1),scon(i,j,k,2),scon(i,j,k,3), sc, &
                   Bconsx_tmp, Bconsy_tmp, Bconsz_tmp,rho_tmp,&
                   velx_tmp,vely_tmp,velz_tmp,eps_tmp,press_tmp,&
                   Bvecx_tmp,Bvecy_tmp,Bvecz_tmp,b2,w_lorentz_tmp,&
                   g11(i,j,k),g12(i,j,k),g13(i,j,k),g22(i,j,k),g23(i,j,k),g33(i,j,k), &
                   uxx,uxy,uxz,uyy,uyz,uzz,sdet, &
                   epsnegative,GRHydro_C2P_failed(i,j,k))

#endif          
              
           end if    ! if (epsnegative .ne. 0) then  

           rho(i,j,k) = rho_tmp 
           press(i,j,k) = press_tmp 
           eps(i,j,k) = eps_tmp 
           vup(i,j,k,1) = velx_tmp 
           vup(i,j,k,2) = vely_tmp 
           vup(i,j,k,3) = velz_tmp 
           w_lorentz(i,j,k) = w_lorentz_tmp 
           Bprim(i,j,k,1) = Bvecx_tmp 
           Bprim(i,j,k,2) = Bvecy_tmp 
           Bprim(i,j,k,3) = Bvecz_tmp 
          
           if ( rho(i,j,k) .le. GRHydro_rho_min*(1.d0+GRHydro_atmo_tolerance)) then
!           if ( rho(i,j,k) .le. GRHydro_rho_min*(1.d0+GRHydro_atmo_tolerance) .or. GRHydro_C2P_failed(i,j,k) .ge. 1) then
               
              b2=g11(i,j,k)*Bprim(i,j,k,1)**2+g22(i,j,k)*Bprim(i,j,k,2)**2+g33(i,j,k)*Bprim(i,j,k,3)**2+ &
                   2.0*(g12(i,j,k)*Bprim(i,j,k,1)*Bprim(i,j,k,2)+g13(i,j,k)*Bprim(i,j,k,1)*Bprim(i,j,k,3)+ &
                   g23(i,j,k)*Bprim(i,j,k,2)*Bprim(i,j,k,3))
              
              dens(i,j,k) = sdet*GRHydro_rho_min !/(1.d0+GRHydro_atmo_tolerance)
              rho(i,j,k) = GRHydro_rho_min
              scon(i,j,k,:) = 0.d0
              vup(i,j,k,:) = 0.d0
              w_lorentz(i,j,k) = 1.d0

              if(evolve_temper.ne.0) then
                ! set hot atmosphere values
                temperature(i,j,k) = grhydro_hot_atmo_temp
                y_e(i,j,k) = grhydro_hot_atmo_Y_e
                y_e_con(i,j,k) = y_e(i,j,k) * dens(i,j,k)
                keytemp = 1
                call EOS_Omni_press(GRHydro_eos_handle,keytemp,GRHydro_eos_rf_prec,n,&
                   rho(i,j,k),eps(i,j,k),temperature(i,j,k),y_e(i,j,k),&
                   press(i,j,k),keyerr,anyerr)
              else
                keytemp = 0
                call EOS_Omni_press(GRHydro_polytrope_handle,keytemp,GRHydro_eos_rf_prec,n,&
                   rho(i,j,k),eps(i,j,k),xtemp,xye,press(i,j,k),keyerr,anyerr)
                call EOS_Omni_EpsFromPress(GRHydro_polytrope_handle,keytemp,GRHydro_eos_rf_prec,n,&
                   rho(i,j,k),eps(i,j,k),xtemp,xye,press(i,j,k),eps(i,j,k),keyerr,anyerr)
              endif

              ! w_lorentz=1, so the expression for tau reduces to:

              !!$ tau does need to take into account the existing B-field
              !!$ with w_lorentz=1, we find tau = sqrtdet*(rho (1+eps+b^2/2)) - dens  [Press drops out]
              tau(i,j,k)  = sdet * (rho(i,j,k)*eps(i,j,k)+b2/2.0)

           end if

        end do
     end do
  end do
  !$OMP END PARALLEL DO
  
  return
#undef faulty_gxx
#undef faulty_gxy
#undef faulty_gxz
#undef faulty_gyy
#undef faulty_gyz
#undef faulty_gzz
#undef faulty_betax
#undef faulty_betay
#undef faulty_betaz
#undef faulty_vel
#undef faulty_Bvec
  
end subroutine Conservative2PrimitiveAM


 /*@@
   @routine    Conservative2PrimitiveBoundariesAM
   @date       Sep 15, 2010
   @author     Scott Noble, Joshua Faber, Bruno Mundim, The GRHydro Developers    
   @desc 
        This routine is used only if the reconstruction is performed on the conserved variables. 
        It computes the primitive variables on cell boundaries.
        Since reconstruction on conservative had not proved to be very successful, 
        some of the improvements to the C2P routines (e.g. the check about 
        whether a failure happens in a point that will be restriced anyway) 
        are not implemented here yet.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/


subroutine Conservative2PrimitiveBoundsAM(CCTK_ARGUMENTS)
  
  use Con2PrimM_fortran_interfaces
  
  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  integer :: i, j, k, itracer, nx, ny, nz
  CCTK_REAL :: uxxl, uxyl, uxzl, uyyl, uyzl, uzzl,&
       uxxr, uxyr, uxzr, uyyr, uyzr, uzzr, pmin(1), epsmin(1)
  CCTK_REAL :: gxxl,gxyl,gxzl,gyyl,gyzl,gzzl,avg_sdetl,&
       gxxr,gxyr,gxzr,gyyr,gyzr,gzzr,avg_sdetr
  CCTK_REAL :: b2minus, b2plus, local_gam(1), local_pgam,local_K,scminus,scplus
  CCTK_INT :: epsnegative
  CCTK_REAL :: Bconsx_tmp, Bconsy_tmp, Bconsz_tmp
  character(len=100) warnline
 
  CCTK_REAL :: local_min_tracer

! begin EOS Omni vars                                            
  CCTK_INT  :: n,keytemp,anyerr,keyerr(1)
  CCTK_REAL :: xpress(1),xtemp(1),xye(1),xeps(1),xrho(1),one(1)=1.0d0
! end EOS Omni vars                       

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
#define gxx faulty_gxx
#define gxy faulty_gxy
#define gxz faulty_gxz
#define gyy faulty_gyy
#define gyz faulty_gyz
#define gzz faulty_gzz
#define betax faulty_betax
#define betay faulty_betay
#define betaz faulty_betaz
#define vel faulty_vel
#define Bvec faulty_Bvec

! begin EOS Omni vars                                            
  n=1;keytemp=0;anyerr=0;keyerr(1)=0
  xpress(1)=0.0d0;xeps(1)=0.0d0;xtemp(1)=0.0d0;xye(1)=0.0d0
! end EOS Omni vars                       

  ! this is a poly call
  xrho(1)=GRHydro_rho_min 
  call EOS_Omni_press(GRHydro_polytrope_handle,keytemp,GRHydro_eos_rf_prec,n,&
       xrho,one,xtemp,xye,pmin,keyerr,anyerr)

  call EOS_Omni_EpsFromPress(GRHydro_polytrope_handle,keytemp,GRHydro_eos_rf_prec,n,&
       xrho,epsmin,xtemp,xye,pmin,epsmin,keyerr,anyerr)

  call EOS_Omni_press(GRHydro_eos_handle,keytemp,GRHydro_eos_rf_prec,n,&
         one,one,xtemp,xye,local_gam,keyerr,anyerr)
  local_gam=local_gam+1.0

  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3) 
  
  if (use_min_tracer .ne. 0) then
     local_min_tracer = min_tracer
  else
     local_min_tracer = 0d0
  end if

  do k = GRHydro_stencil, nz - GRHydro_stencil + 1
     do j = GRHydro_stencil, ny - GRHydro_stencil + 1
        do i = GRHydro_stencil, nx - GRHydro_stencil + 1
           
           !do not compute if in atmosphere or in an excised region
           if ((atmosphere_mask(i,j,k) .ne. 0) .or. &
                GRHydro_enable_internal_excision /= 0 .and. (hydro_excision_mask(i,j,k) .ne. 0)) cycle
           
           gxxl = 0.5d0 * (g11(i,j,k) + g11(i-xoffset,j-yoffset,k-zoffset))
           gxyl = 0.5d0 * (g12(i,j,k) + g12(i-xoffset,j-yoffset,k-zoffset))
           gxzl = 0.5d0 * (g13(i,j,k) + g13(i-xoffset,j-yoffset,k-zoffset))
           gyyl = 0.5d0 * (g22(i,j,k) + g22(i-xoffset,j-yoffset,k-zoffset))
           gyzl = 0.5d0 * (g23(i,j,k) + g23(i-xoffset,j-yoffset,k-zoffset))
           gzzl = 0.5d0 * (g33(i,j,k) + g33(i-xoffset,j-yoffset,k-zoffset))
           gxxr = 0.5d0 * (g11(i,j,k) + g11(i+xoffset,j+yoffset,k+zoffset))
           gxyr = 0.5d0 * (g12(i,j,k) + g12(i+xoffset,j+yoffset,k+zoffset))
           gxzr = 0.5d0 * (g13(i,j,k) + g13(i+xoffset,j+yoffset,k+zoffset))
           gyyr = 0.5d0 * (g22(i,j,k) + g22(i+xoffset,j+yoffset,k+zoffset))
           gyzr = 0.5d0 * (g23(i,j,k) + g23(i+xoffset,j+yoffset,k+zoffset))
           gzzr = 0.5d0 * (g33(i,j,k) + g33(i+xoffset,j+yoffset,k+zoffset))
           
           epsnegative = 0
           
           avg_sdetl = sqrt(SPATIAL_DETERMINANT(gxxl,gxyl,gxzl,gyyl, gyzl,gzzl))
           avg_sdetr = sqrt(SPATIAL_DETERMINANT(gxxr,gxyr,gxzr,gyyr, gyzr,gzzr))
           call UpperMetric(uxxl,uxyl,uxzl,uyyl,uyzl,uzzl,avg_sdetl*avg_sdetl,&
                gxxl,gxyl,gxzl,gyyl,gyzl,gzzl)        
           call UpperMetric(uxxr,uxyr,uxzr,uyyr,uyzr,uzzr,avg_sdetr*avg_sdetr,&
                gxxr,gxyr,gxzr,gyyr,gyzr,gzzr)        
           
!!$ Tracers get no update for MHD!
           if (evolve_tracer .ne. 0) then
              do itracer=1,number_of_tracers
                 call Con2Prim_ptTracer(cons_tracer(i,j,k,itracer), &
                      tracer(i,j,k,itracer), dens(i,j,k))
              enddo
              
              if (use_min_tracer .ne. 0) then
                 if (tracer(i,j,k,itracer) .le. local_min_tracer) then
                    tracer(i,j,k,itracer) = local_min_tracer
                 end if
              end if
              
           endif

           if(evolve_Y_e.ne.0) then
              Y_e(i,j,k) = Y_e_con(i,j,k) / dens(i,j,k)
           endif

           Bconsx_tmp = avg_sdetl*Bvecxminus(i,j,k)
           Bconsy_tmp = avg_sdetl*Bvecyminus(i,j,k)
           Bconsz_tmp = avg_sdetl*Bveczminus(i,j,k)
           call GRHydro_Con2PrimM_pt2(GRHydro_eos_handle, keytemp, GRHydro_eos_rf_prec, GRHydro_perc_ptol, local_gam(1), densminus(i,j,k), &
                sxminus(i,j,k),syminus(i,j,k),szminus(i,j,k), tauminus(i,j,k), &
                Bconsx_tmp, Bconsy_tmp, Bconsz_tmp,  xye(1), xtemp(1), rhominus(i,j,k),&
                velxminus(i,j,k),velyminus(i,j,k),velzminus(i,j,k),epsminus(i,j,k),pressminus(i,j,k),&
                Bvecxminus(i,j,k), Bvecyminus(i,j,k), Bveczminus(i,j,k),b2minus, w_lorentzminus(i,j,k),&
                gxxl,gxyl,gxzl,gyyl,gyzl,gzzl, &
                uxxl,uxyl,uxzl,uyyl,uyzl,uzzl,avg_sdetl, &
                epsnegative,GRHydro_C2P_failed(i,j,k))
           
           if (epsnegative .ne. 0) then
              !$OMP CRITICAL
              call CCTK_WARN(GRHydro_NaN_verbose+2, 'Specific internal energy just went below 0, trying polytype!')
              !$OMP END CRITICAL
              
              xrho=10.0d0
              call EOS_Omni_press(GRHydro_polytrope_handle,keytemp,GRHydro_eos_rf_prec,n,&
                   one,one,xtemp,xye,xpress,keyerr,anyerr)
              local_K = xpress(1)

              call EOS_Omni_press(GRHydro_polytrope_handle,keytemp,GRHydro_eos_rf_prec,n,&
                   xrho,one,xtemp,xye,xpress,keyerr,anyerr)
              local_pgam=log(xpress(1)/local_K)/log(xrho(1))
              scminus = local_K*densminus(i,j,k)

              call GRHydro_Con2PrimM_Polytype_pt(GRHydro_eos_handle, local_pgam, densminus(i,j,k), &
                   sxminus(i,j,k),syminus(i,j,k),szminus(i,j,k), scminus, &
                   Bconsx_tmp, Bconsy_tmp, Bconsz_tmp, rhominus(i,j,k),&
                   velxminus(i,j,k),velyminus(i,j,k),velzminus(i,j,k),epsminus(i,j,k),pressminus(i,j,k),&
                   Bvecxminus(i,j,k), Bvecyminus(i,j,k), Bveczminus(i,j,k),b2minus, w_lorentzminus(i,j,k),&
                   gxxl,gxyl,gxzl,gyyl,gyzl,gzzl, &
                   uxxl,uxyl,uxzl,uyyl,uyzl,uzzl,avg_sdetl, &
                   epsnegative,GRHydro_C2P_failed(i,j,k))
          
           end if
           
           if (epsminus(i,j,k) .lt. 0.0d0) then
              if (GRHydro_reflevel.ge.GRHydro_c2p_warn_from_reflevel) then
                 !$OMP CRITICAL
                 call CCTK_WARN(1,'Con2Prim: stopping the code.')
                 call CCTK_WARN(1, '   specific internal energy just went below 0! ')
                 write(warnline,'(a28,i2)') 'on carpet reflevel: ',GRHydro_reflevel
                 call CCTK_WARN(1,warnline)
                 write(warnline,'(a20,3g16.7)') 'xyz location: ',&
                      x(i,j,k),y(i,j,k),z(i,j,k)
                 call CCTK_WARN(1,warnline)
                 write(warnline,'(a20,g16.7)') 'radius: ',r(i,j,k)
                 call CCTK_WARN(1,warnline)
                 write(warnline,'(a20,3g16.7)') 'velocities: ',&
                      velxminus(i,j,k),velyminus(i,j,k),velzminus(i,j,k)
                 call CCTK_WARN(1,warnline)
                 call CCTK_WARN(GRHydro_c2p_warnlevel, "Specific internal energy negative")
                 !$OMP END CRITICAL
                 exit
              endif
           endif
           
           epsnegative = 0

           Bconsx_tmp = avg_sdetr*Bvecxplus(i,j,k)
           Bconsy_tmp = avg_sdetr*Bvecyplus(i,j,k)
           Bconsz_tmp = avg_sdetr*Bveczplus(i,j,k)
           call GRHydro_Con2PrimM_pt2(GRHydro_eos_handle, keytemp, GRHydro_eos_rf_prec, GRHydro_perc_ptol, local_gam(1), densplus(i,j,k), &
                sxplus(i,j,k),syplus(i,j,k),szplus(i,j,k), tauplus(i,j,k), &
                Bconsx_tmp, Bconsy_tmp, Bconsz_tmp,  xye(1), xtemp(1), rhoplus(i,j,k),&
                velxplus(i,j,k),velyplus(i,j,k),velzplus(i,j,k),epsplus(i,j,k),pressplus(i,j,k),&
                Bvecxplus(i,j,k), Bvecyplus(i,j,k), Bveczplus(i,j,k),b2plus, w_lorentzplus(i,j,k),&
                gxxr,gxyr,gxzr,gyyr,gyzr,gzzr, &
                uxxr,uxyr,uxzr,uyyr,uyzr,uzzr,avg_sdetr, &
                epsnegative,GRHydro_C2P_failed(i,j,k))
 
           if (epsnegative .ne. 0) then
              !$OMP CRITICAL
              call CCTK_WARN(GRHydro_NaN_verbose+2, 'Specific internal energy just went below 0, trying polytype!!')
              !$OMP END CRITICAL

              xrho=10.0d0
              call EOS_Omni_press(GRHydro_polytrope_handle,keytemp,GRHydro_eos_rf_prec,n,&
                   one,one,xtemp,xye,xpress,keyerr,anyerr)
              local_K = xpress(1)

              call EOS_Omni_press(GRHydro_polytrope_handle,keytemp,GRHydro_eos_rf_prec,n,&
                   xrho,one,xtemp,xye,xpress,keyerr,anyerr)
              local_pgam=log(xpress(1)/local_K)/log(xrho(1))
              scplus = local_K*densplus(i,j,k)

              call GRHydro_Con2PrimM_Polytype_pt(GRHydro_eos_handle, local_pgam, densplus(i,j,k), &
                   sxplus(i,j,k),syplus(i,j,k),szplus(i,j,k), scplus,&
                   Bconsx_tmp, Bconsy_tmp, Bconsz_tmp, rhoplus(i,j,k),&
                   velxplus(i,j,k),velyplus(i,j,k),velzplus(i,j,k),epsplus(i,j,k),pressplus(i,j,k),&
                   Bvecxplus(i,j,k), Bvecyplus(i,j,k), Bveczplus(i,j,k),b2plus, w_lorentzplus(i,j,k),&
                   gxxr,gxyr,gxzr,gyyr,gyzr,gzzr, &
                   uxxr,uxyr,uxzr,uyyr,uyzr,uzzr,avg_sdetr, &
                   epsnegative,GRHydro_C2P_failed(i,j,k))
           end if
           
           if (epsplus(i,j,k) .lt. 0.0d0) then
              if (GRHydro_reflevel.ge.GRHydro_c2p_warn_from_reflevel) then
                 !$OMP CRITICAL
                 call CCTK_WARN(1,'Con2Prim: stopping the code.')
                 call CCTK_WARN(1, '   specific internal energy just went below 0! ')
                 write(warnline,'(a28,i2)') 'on carpet reflevel: ',GRHydro_reflevel
                 call CCTK_WARN(1,warnline)
                 write(warnline,'(a20,3g16.7)') 'xyz location: ',&
                      x(i,j,k),y(i,j,k),z(i,j,k)
                 call CCTK_WARN(1,warnline)
                 write(warnline,'(a20,g16.7)') 'radius: ',r(i,j,k)
                 call CCTK_WARN(1,warnline)
                 write(warnline,'(a20,3g16.7)') 'velocities: ',&
                      velxplus(i,j,k),velyplus(i,j,k),velzplus(i,j,k)
                 call CCTK_WARN(1,warnline)
                 call CCTK_WARN(GRHydro_c2p_warnlevel, "Specific internal energy negative")
                 write(warnline,'(a25,4g15.6)') 'coordinates: x,y,z,r:',&
                      x(i,j,k),y(i,j,k),z(i,j,k),r(i,j,k)
                 call CCTK_WARN(1,warnline)
                 !$OMP END CRITICAL
              endif
           endif
        end do
     end do
  end do
  
#undef faulty_gxx
#undef faulty_gxy
#undef faulty_gxz
#undef faulty_gyy
#undef faulty_gyz
#undef faulty_gzz
#undef faulty_betax
#undef faulty_betay
#undef faulty_betaz
#undef faulty_vel
#undef faulty_Bvec
  
end subroutine Conservative2PrimitiveBoundsAM

/*@@
@routine    Con2PrimPolytypeAM
@date       Sep 16, 2010
@author     SCott Noble, Joshua Faber, Bruno Mundim, Ian Hawke
@desc 
All routines below are identical to those above, just
specialised from polytropic type EOS.
@enddesc 
@calls     
@calledby   
@history 

@endhistory 

@@*/

subroutine Conservative2PrimitivePolytypeAM(CCTK_ARGUMENTS)
  
  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  !!DECLARE_CCTK_FUNCTIONS
  
  integer :: i, j, k, itracer, nx, ny, nz
  CCTK_REAL :: uxx, uxy, uxz, uyy, uyz, uzz, sdet,b2

  CCTK_INT :: epsnegative
  CCTK_REAL :: Bconsx_tmp, Bconsy_tmp, Bconsz_tmp
  
  CCTK_REAL :: local_min_tracer, local_pgam,local_K, sc
  !  character(len=400) :: warnline
  
! begin EOS Omni vars                                                                                       
  CCTK_INT  :: n,keytemp,anyerr,keyerr(1)
  CCTK_REAL :: xpress,xtemp,xye,xeps,xrho
! end EOS Omni vars                 

  ! save memory when MP is not used
  CCTK_INT :: GRHydro_UseGeneralCoordinates
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3) :: g11, g12, g13, g22, g23, g33
  pointer (pg11,g11), (pg12,g12), (pg13,g13), (pg22,g22), (pg23,g23), (pg33,g33)
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3,3) :: vup, Bprim
  pointer (pvup,vup), (pBprim,Bprim)

  call GRHydro_BvecfromAvec(CCTK_PASS_FTOF)

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
#define gxx faulty_gxx
#define gxy faulty_gxy
#define gxz faulty_gxz
#define gyy faulty_gyy
#define gyz faulty_gyz
#define gzz faulty_gzz
#define betax faulty_betax
#define betay faulty_betay
#define betaz faulty_betaz
#define vel faulty_vel
#define Bvec faulty_Bvec

! begin EOS Omni vars                                                                                       
  n=1;keytemp=0;anyerr=0;keyerr(1)=0
  xpress=0.0d0;xtemp=0.0d0;xye=0.0d0;xeps=0.0d0
! end EOS Omni vars                 
  
  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)
  
  if (use_min_tracer .ne. 0) then
     local_min_tracer = min_tracer
  else
     local_min_tracer = 0d0
  end if
  
!!$  do k = GRHydro_stencil + 1, nz - GRHydro_stencil
!!$    do j = GRHydro_stencil + 1, ny - GRHydro_stencil
!!$      do i = GRHydro_stencil + 1, nx - GRHydro_stencil


  !$OMP PARALLEL DO PRIVATE(i,j,k,itracer,&
  !$OMP uxx, uxy, uxz, uyy, uyz, uzz, sdet, epsnegative, &
  !$OMP Bconsx_tmp, Bconsy_tmp, Bconsz_tmp, &
  !$OMP b2, xrho, xpress, local_K, local_pgam, sc)
  do k = 1, nz
     do j = 1, ny
        do i = 1, nx
           
           !do not compute if in atmosphere or in an excised region
           if ((atmosphere_mask(i,j,k) .ne. 0) .or. &
                GRHydro_enable_internal_excision /= 0 .and. (hydro_excision_mask(i,j,k) .ne. 0)) cycle
           
           sdet = sdetg(i,j,k)
           call UpperMetric(uxx,uxy,uxz,uyy,uyz,uzz,sdet*sdet,&
                g11(i,j,k),g12(i,j,k),g13(i,j,k),g22(i,j,k),&
                g23(i,j,k),g33(i,j,k))        

!!$ No MHD changes to tracers
           if (evolve_tracer .ne. 0) then
              do itracer=1,number_of_tracers
                 call Con2Prim_ptTracer(cons_tracer(i,j,k,itracer), & 
                      tracer(i,j,k,itracer), dens(i,j,k))
              enddo
              
              if (use_min_tracer .ne. 0) then
                 if (tracer(i,j,k,itracer) .le. local_min_tracer) then
                    tracer(i,j,k,itracer) = local_min_tracer
                 end if
              end if
              
           endif

           if(evolve_Y_e.ne.0) then
              Y_e(i,j,k) = Y_e_con(i,j,k) / dens(i,j,k)
           endif
           
           xrho=10.0d0
           call EOS_Omni_press(GRHydro_polytrope_handle,keytemp,GRHydro_eos_rf_prec,n,&
                1.d0,1.0d0,xtemp,xye,xpress,keyerr,anyerr)
           local_K = xpress
           
           call EOS_Omni_press(GRHydro_polytrope_handle,keytemp,GRHydro_eos_rf_prec,n,&
                xrho,1.0d0,xtemp,xye,xpress,keyerr,anyerr)
           local_pgam=log(xpress/local_K)/log(xrho)
           sc = local_K*dens(i,j,k)
           Bconsx_tmp = sdet*Bprim(i,j,k,1)
           Bconsy_tmp = sdet*Bprim(i,j,k,2)
           Bconsz_tmp = sdet*Bprim(i,j,k,3)
           
           call GRHydro_Con2PrimM_Polytype_pt(GRHydro_eos_handle, local_pgam, dens(i,j,k), &
                   scon(i,j,k,1),scon(i,j,k,2),scon(i,j,k,3), sc, &
                   Bconsx_tmp, Bconsy_tmp, Bconsz_tmp, rho(i,j,k),&
                   vup(i,j,k,1),vup(i,j,k,2),vup(i,j,k,3),eps(i,j,k),press(i,j,k),&
                   Bprim(i,j,k,1), Bprim(i,j,k,2), Bprim(i,j,k,3),b2, w_lorentz(i,j,k),&
                   g11(i,j,k),g12(i,j,k),g13(i,j,k),g22(i,j,k),g23(i,j,k),g33(i,j,k), &
                   uxx,uxy,uxz,uyy,uyz,uzz,sdet, &
                   epsnegative,GRHydro_C2P_failed(i,j,k))
            
        end do
     end do
  end do
  
  !$OMP END PARALLEL DO
  
  return

#undef faulty_gxx
#undef faulty_gxy
#undef faulty_gxz
#undef faulty_gyy
#undef faulty_gyz
#undef faulty_gzz
#undef faulty_betax
#undef faulty_betay
#undef faulty_betaz
#undef faulty_vel
#undef faulty_Bvec
  
end subroutine Conservative2PrimitivePolytypeAM


 /*@@
   @routine    Cons2PrimBoundsPolytypeAM
   @date       Sep 16, 2010
   @author     Scott Noble, Joshua Faber, Bruno Mundim, The GRHydro Developers
   @desc 
        This routine is used only if the reconstruction is performed on the conserved variables. 
        It computes the primitive variables on cell boundaries.
        Since reconstruction on conservative had not proved to be very successful, 
        some of the improvements to the C2P routines (e.g. the check about 
        whether a failure happens in a point that will be restriced anyway) are not implemented here yet.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine Con2PrimBoundsPolytypeAM(CCTK_ARGUMENTS)
  
  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  !!DECLARE_CCTK_FUNCTIONS
  
  integer :: i, j, k, nx, ny, nz
  CCTK_REAL :: uxxl, uxyl, uxzl, uyyl, uyzl, uzzl,&
       uxxr, uxyr, uxzr, uyyr, uyzr, uzzr
  CCTK_REAL :: gxxl,gxyl,gxzl,gyyl,gyzl,gzzl,avg_sdetl,&
       gxxr,gxyr,gxzr,gyyr,gyzr,gzzr,avg_sdetr
  CCTK_REAL :: b2minus, b2plus
  CCTK_REAL :: Bconsx_tmp, Bconsy_tmp, Bconsz_tmp
  CCTK_INT :: epsnegative

  CCTK_REAL :: local_pgam,local_K,scplus,scminus

! begin EOS Omni vars                                                                                       
  CCTK_INT  :: n,keytemp,anyerr,keyerr(1)
  CCTK_REAL :: xpress,xtemp,xye,xeps,xrho
! end EOS Omni vars                       

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
#define gxx faulty_gxx
#define gxy faulty_gxy
#define gxz faulty_gxz
#define gyy faulty_gyy
#define gyz faulty_gyz
#define gzz faulty_gzz
#define betax faulty_betax
#define betay faulty_betay
#define betaz faulty_betaz
#define vel faulty_vel
#define Bvec faulty_Bvec

! begin EOS Omni vars                                            
  n=1;keytemp=0;anyerr=0;keyerr(1)=0
  xpress=0.0d0;xtemp=0.0d0;xye=0.0d0;xeps=0.0d0
! end EOS Omni vars                 
  
  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)
  
  do k = GRHydro_stencil, nz - GRHydro_stencil + 1
     do j = GRHydro_stencil, ny - GRHydro_stencil + 1
        do i = GRHydro_stencil, nx - GRHydro_stencil + 1
           
           !do not compute if in atmosphere or in an excised region
           if ((atmosphere_mask(i,j,k) .ne. 0) .or. &
                GRHydro_enable_internal_excision /= 0 .and. (hydro_excision_mask(i,j,k) .ne. 0)) cycle
           
           gxxl = 0.5d0 * (g11(i,j,k) + g11(i-xoffset,j-yoffset,k-zoffset))
           gxyl = 0.5d0 * (g12(i,j,k) + g12(i-xoffset,j-yoffset,k-zoffset))
           gxzl = 0.5d0 * (g13(i,j,k) + g13(i-xoffset,j-yoffset,k-zoffset))
           gyyl = 0.5d0 * (g22(i,j,k) + g22(i-xoffset,j-yoffset,k-zoffset))
           gyzl = 0.5d0 * (g23(i,j,k) + g23(i-xoffset,j-yoffset,k-zoffset))
           gzzl = 0.5d0 * (g33(i,j,k) + g33(i-xoffset,j-yoffset,k-zoffset))
           gxxr = 0.5d0 * (g11(i,j,k) + g11(i+xoffset,j+yoffset,k+zoffset))
           gxyr = 0.5d0 * (g12(i,j,k) + g12(i+xoffset,j+yoffset,k+zoffset))
           gxzr = 0.5d0 * (g13(i,j,k) + g13(i+xoffset,j+yoffset,k+zoffset))
           gyyr = 0.5d0 * (g22(i,j,k) + g22(i+xoffset,j+yoffset,k+zoffset))
           gyzr = 0.5d0 * (g23(i,j,k) + g23(i+xoffset,j+yoffset,k+zoffset))
           gzzr = 0.5d0 * (g33(i,j,k) + g33(i+xoffset,j+yoffset,k+zoffset))
           
           avg_sdetl = sqrt(SPATIAL_DETERMINANT(gxxl,gxyl,gxzl,gyyl, gyzl,gzzl))
           avg_sdetr = sqrt(SPATIAL_DETERMINANT(gxxr,gxyr,gxzr,gyyr, gyzr,gzzr))
           call UpperMetric(uxxl,uxyl,uxzl,uyyl,uyzl,uzzl,avg_sdetl*avg_sdetl,&
                gxxl,gxyl,gxzl,gyyl,gyzl,gzzl)        
           call UpperMetric(uxxr,uxyr,uxzr,uyyr,uyzr,uzzr,avg_sdetr*avg_sdetr,&
                gxxr,gxyr,gxzr,gyyr,gyzr,gzzr)        

           xrho=10.0d0
           call EOS_Omni_press(GRHydro_polytrope_handle,keytemp,GRHydro_eos_rf_prec,n,&
                1.d0,1.0d0,xtemp,xye,xpress,keyerr,anyerr)
           local_K = xpress
           
           call EOS_Omni_press(GRHydro_polytrope_handle,keytemp,GRHydro_eos_rf_prec,n,&
                xrho,1.0d0,xtemp,xye,xpress,keyerr,anyerr)
           local_pgam=log(xpress/local_K)/log(xrho)
           scminus = local_K*densminus(i,j,k)
           scplus = local_K*densplus(i,j,k)
           Bconsx_tmp = avg_sdetl*Bvecxminus(i,j,k)
           Bconsy_tmp = avg_sdetl*Bvecyminus(i,j,k)
           Bconsz_tmp = avg_sdetl*Bveczminus(i,j,k)

           call GRHydro_Con2PrimM_Polytype_pt(GRHydro_eos_handle, local_pgam,densminus(i,j,k), &
                sxminus(i,j,k),syminus(i,j,k),szminus(i,j,k), scminus,&
                Bconsx_tmp, Bconsy_tmp, Bconsz_tmp, rhominus(i,j,k),&
                velxminus(i,j,k),velyminus(i,j,k),velzminus(i,j,k),epsminus(i,j,k),pressminus(i,j,k),&
                Bvecxminus(i,j,k), Bvecyminus(i,j,k), Bveczminus(i,j,k),b2minus, w_lorentzminus(i,j,k),&
                gxxl,gxyl,gxzl,gyyl,gyzl,gzzl, &
                uxxl,uxyl,uxzl,uyyl,uyzl,uzzl,avg_sdetl, &
                epsnegative,GRHydro_C2P_failed(i,j,k))

           Bconsx_tmp = avg_sdetr*Bvecxplus(i,j,k)
           Bconsy_tmp = avg_sdetr*Bvecyplus(i,j,k)
           Bconsz_tmp = avg_sdetr*Bveczplus(i,j,k)
           call GRHydro_Con2PrimM_Polytype_pt(GRHydro_eos_handle, local_pgam,densplus(i,j,k), &
                sxplus(i,j,k),syplus(i,j,k),szplus(i,j,k), scplus,&
                Bconsx_tmp, Bconsy_tmp, Bconsz_tmp, rhoplus(i,j,k),&
                velxplus(i,j,k),velyplus(i,j,k),velzplus(i,j,k),epsplus(i,j,k),pressplus(i,j,k),&
                Bvecxplus(i,j,k), Bvecyplus(i,j,k),Bveczplus(i,j,k),b2plus,w_lorentzplus(i,j,k),&
                gxxr,gxyr,gxzr,gyyr,gyzr,gzzr, &
                uxxr,uxyr,uxzr,uyyr,uyzr,uzzr,avg_sdetr, &
                epsnegative,GRHydro_C2P_failed(i,j,k))
        end do
     end do
  end do
  
#undef faulty_gxx
#undef faulty_gxy
#undef faulty_gxz
#undef faulty_gyy
#undef faulty_gyz
#undef faulty_gzz
#undef faulty_betax
#undef faulty_betay
#undef faulty_betaz
#undef faulty_vel
#undef faulty_Bvec
  
end subroutine Con2PrimBoundsPolytypeAM

!!$ Con2Prim_ptTracer, Con2Prim_BoundsTracer, and Con2Prim_ptBoundsTracer need not be rewritten!

