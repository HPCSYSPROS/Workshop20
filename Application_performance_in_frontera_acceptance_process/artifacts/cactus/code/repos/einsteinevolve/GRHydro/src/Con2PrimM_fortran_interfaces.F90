#include "cctk.h" 

module Con2PrimM_fortran_interfaces
  implicit none
  
  CCTK_INT, parameter :: izero = 0
  integer, parameter :: ik = kind(izero)

  interface
    subroutine GRHydro_Con2PrimM_ptold(  handle, &
                                       local_gam, dens, &
                        sx,  sy,  sz, &
                        tau, &
                        Bconsx,  Bconsy,  Bconsz, &
                        rho, &
                        velx,  vely,  velz,&
                        epsilon,  pressure,&
                        Bx,  By,  Bz, &
                        bsq,&
                        w_lorentz, &
                        gxx,  gxy,  gxz, &
                        gyy,  gyz,  gzz, &
                        uxx, uxy, uxz,&
                        uyy, uyz, uzz,&
                        det,&
                        epsnegative, &
                        retval)

       implicit none       
       CCTK_INT handle
       CCTK_REAL local_gam
       CCTK_REAL dens 
       CCTK_REAL sx, sy, sz
       CCTK_REAL tau 
       CCTK_REAL Bconsx, Bconsy, Bconsz
       CCTK_REAL rho
       CCTK_REAL velx, vely, velz
       CCTK_REAL epsilon, pressure
       CCTK_REAL Bx, By, Bz
       CCTK_REAL bsq
       CCTK_REAL w_lorentz
       CCTK_REAL gxx, gxy, gxz 
       CCTK_REAL gyy, gyz, gzz
       CCTK_REAL uxx, uxy, uxz
       CCTK_REAL uyy, uyz, uzz
       CCTK_REAL det
       CCTK_INT epsnegative
       CCTK_REAL retval
     end subroutine  GRHydro_Con2PrimM_ptold

     subroutine GRHydro_Con2PrimM_pt2(  handle, keytemp, eos_prec, prec, &
                                       local_gam, dens, &
                        sx,  sy,  sz, &
                        tau, &
                        Bconsx,  Bconsy,  Bconsz, &
                        xye, xtemp, &
                        rho, &
                        velx,  vely,  velz,&
                        epsilon,  pressure,&
                        Bx,  By,  Bz, &
                        bsq,&
                        w_lorentz, &
                        gxx,  gxy,  gxz, &
                        gyy,  gyz,  gzz, &
                        uxx, uxy, uxz,&
                        uyy, uyz, uzz,&
                        det,&
                        epsnegative, &
                        retval)

       implicit none       
       CCTK_INT handle
       CCTK_INT keytemp
       CCTK_REAL eos_prec
       CCTK_REAL prec
       CCTK_REAL local_gam
       CCTK_REAL dens 
       CCTK_REAL sx, sy, sz
       CCTK_REAL tau 
       CCTK_REAL Bconsx, Bconsy, Bconsz
       CCTK_REAL xye
       CCTK_REAL xtemp 
       CCTK_REAL rho
       CCTK_REAL velx, vely, velz
       CCTK_REAL epsilon, pressure
       CCTK_REAL Bx, By, Bz
       CCTK_REAL bsq
       CCTK_REAL w_lorentz
       CCTK_REAL gxx, gxy, gxz 
       CCTK_REAL gyy, gyz, gzz
       CCTK_REAL uxx, uxy, uxz
       CCTK_REAL uyy, uyz, uzz
       CCTK_REAL det
       CCTK_INT epsnegative
       CCTK_REAL retval
     end subroutine  GRHydro_Con2PrimM_pt2


     subroutine GRHydro_Con2PrimM_pt(  handle, GRHydro_reflevel, i, j, k, x, y, z, keytemp, eos_prec, prec, &
                                       local_gam, dens, &
                        sx,  sy,  sz, &
                        tau, &
                        Bconsx,  Bconsy,  Bconsz, &
                        xye, xtemp, &
                        rho, &
                        velx,  vely,  velz,&
                        epsilon,  pressure,&
                        Bx,  By,  Bz, &
                        bsq,&
                        w_lorentz, &
                        gxx,  gxy,  gxz, &
                        gyy,  gyz,  gzz, &
                        uxx, uxy, uxz,&
                        uyy, uyz, uzz,&
                        det,&
                        epsnegative, &
                        retval)

       implicit none       
       CCTK_INT handle
       CCTK_INT GRHydro_reflevel 
       CCTK_INT i, j, k 
       CCTK_REAL x, y, z
       CCTK_INT keytemp
       CCTK_REAL eos_prec
       CCTK_REAL prec
       CCTK_REAL local_gam
       CCTK_REAL dens 
       CCTK_REAL sx, sy, sz
       CCTK_REAL tau 
       CCTK_REAL Bconsx, Bconsy, Bconsz
       CCTK_REAL xye
       CCTK_REAL xtemp 
       CCTK_REAL rho
       CCTK_REAL velx, vely, velz
       CCTK_REAL epsilon, pressure
       CCTK_REAL Bx, By, Bz
       CCTK_REAL bsq
       CCTK_REAL w_lorentz
       CCTK_REAL gxx, gxy, gxz 
       CCTK_REAL gyy, gyz, gzz
       CCTK_REAL uxx, uxy, uxz
       CCTK_REAL uyy, uyz, uzz
       CCTK_REAL det
       CCTK_INT epsnegative
       CCTK_REAL retval
     end subroutine  GRHydro_Con2PrimM_pt

  subroutine GRHydro_Con2PrimM_Polytype_pt(  handle, local_gam,&
                        dens, &
                        sx,  sy,  sz, &
                        sc, &
                        Bconsx,  Bconsy,  Bconsz, &
                        rho, &
                        velx,  vely,  velz,&
                        epsilon,  pressure,&
                        Bx,  By,  Bz, &
                        bsq,&
                        w_lorentz, &
                        gxx,  gxy,  gxz, &
                        gyy,  gyz,  gzz, &
                        uxx, uxy, uxz,&
                        uyy, uyz, uzz,&
                        det,&
                        epsnegative, &
                        retval)

       implicit none       
       CCTK_INT handle
       CCTK_REAL local_gam
       CCTK_REAL dens 
       CCTK_REAL sx, sy, sz
       CCTK_REAL sc 
       CCTK_REAL Bconsx, Bconsy, Bconsz
       CCTK_REAL rho
       CCTK_REAL velx, vely, velz
       CCTK_REAL epsilon, pressure
       CCTK_REAL Bx, By, Bz
       CCTK_REAL bsq
       CCTK_REAL w_lorentz
       CCTK_REAL gxx, gxy, gxz 
       CCTK_REAL gyy, gyz, gzz
       CCTK_REAL uxx, uxy, uxz
       CCTK_REAL uyy, uyz, uzz
       CCTK_REAL det
       CCTK_INT epsnegative
       CCTK_REAL retval
     end subroutine  GRHydro_Con2PrimM_Polytype_pt

  end interface

end module Con2PrimM_fortran_interfaces
