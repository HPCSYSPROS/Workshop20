 /*@@
   @file      GRHydro_MP5Reconstruct.F90
   @date      Fri Jan 3 2013
   @author    Ian Hawke, Christian Reisswig
   @desc 
   Routines to set up the coefficient array and to perform one dimensional 
   M5 reconstruction of arbitrary order.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"




 /*@@
   @routine    GRHydro_MP5Reconstruct1d
   @date       Fri Jan 3 2013
   @author     Christian Reisswig
   @desc 
   Perform a one dimensional reconstruction of a given array using M5.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

#define SpaceMask_CheckStateBitsF90_1D(mask,i,type_bits,state_bits) \
  (iand(mask((i)),(type_bits)).eq.(state_bits))

subroutine GRHydro_MP5Reconstruct1d(nx, v, vminus, vplus, trivial_rp, &
     hydro_excision_mask)

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: nx, i
  CCTK_REAL, dimension(nx) :: v, vplus, vminus

  CCTK_INT, dimension(nx) :: hydro_excision_mask
  logical, dimension(nx) :: trivial_rp
  logical, dimension(nx) :: excise
  logical :: normal_m5

  CCTK_REAL :: vl, vmp, djm1, dj, djp1, dm4jph, dm4jmh, vul, vav, vmd, vlc, vmin, vmax, vnorm

  ! sign requires its arguments to be of identical KIND
  CCTK_REAL, parameter :: one = 1d0
  
  vminus = 0.d0
  vplus = 0.d0

  excise = .false.
  trivial_rp = .false.

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
    normal_m5 = .true.
    if (i < nx) then
     if (excise(i+1)) then
      vminus(i) = v(i)
      vplus(i) = v(i)
      normal_m5 = .false.
     end if
    end if
    if (i > 1) then
     if (excise(i-1)) then
      vminus(i) = v(i)
      vplus(i) = v(i)
      normal_m5 = .false.
     end if
    end if

    if (normal_m5) then
    
      vnorm = sqrt(v(i-2)**2 + v(i-1)**2 + v(i)**2 + v(i+1)**2 + v(i+2)**2)
    
#define MINMOD(x,y)  \
     (0.5d0*(sign(one,x) + sign(one,y)) * min(abs(x), abs(y)))

#define MINMOD4(w,x,y,z) \
     (0.125d0*( sign(one,w)+sign(one,x) )*abs( (sign(one,w)+sign(one,y)) * (sign(one,w)+sign(one,z)) )*min(abs(w), abs(x), abs(y), abs(z)))

#define MP5(am2, am1, a, ap1, ap2, arecon) \
      vl = (2.0d0*am2 - 13.0d0*am1 + 47.0d0*a + 27.0d0*ap1 - 3.0d0*ap2)/60.0d0  &&\
      vmp = a + MINMOD( ap1-a, mp5_alpha*(a-am1) )  &&\
      if ((vl-a)*(vl-vmp) .le. mp5_eps*vnorm) then &&\
         arecon = vl &&\
      else &&\
        djm1 = am2 -2.0d0*am1 + a &&\
        dj   = am1 -2.0d0*a + ap1 &&\
        djp1 = a -2.d0*ap1 + ap2 &&\
        dm4jph = MINMOD4(4.0d0*dj-djp1, 4.0d0*djp1-dj, dj, djp1) &&\
        dm4jmh = MINMOD4(4.0d0*dj-djm1, 4.0d0*djm1-dj, dj, djm1) &&\
        vul = a + mp5_alpha*(a-am1) &&\
        vav = 0.5d0*(a+ap1) &&\
        vmd = vav - 0.5d0*dm4jph &&\
        vlc = a + 0.5d0*(a-am1) + 4.0d0/3.0d0*dm4jmh &&\
        vmin = max(min(a,ap1,vmd), min(a,vul,vlc)) &&\
        vmax = min(max(a,ap1,vmd), max(a,vul,vlc)) &&\
        arecon = vl + MINMOD(vmin-vl, vmax-vl) &&\
      endif
      
      MP5(v(i-2), v(i-1), v(i), v(i+1), v(i+2), vplus(i))
      MP5(v(i+2), v(i+1), v(i), v(i-1), v(i-2), vminus(i))
      
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

end subroutine GRHydro_MP5Reconstruct1d
