 /*@@
   @file      Utils.F
   @date      Sat Jan 26 02:28:46 2002
   @author    
   @desc 
   Utility functions for other thorns. Calculation of the determinant
   of the spatial metric and the upper metric.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "GRHydro_Macros.h"

 /*@@
   @routine    GRHydro_RefinementLevel
   @date       July 2005
   @author     
   @desc 
   Calculates the current refinement level from flesh data
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 
@@*/

subroutine GRHydro_Debug(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS

  integer i,j,k
  integer nx, ny, nz

  GRHydro_reflevel = int(log10(dble(cctk_levfac(1)))/log10(2.0d0))
  
  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)

  if (GRHydro_reflevel .lt. 3) return
 
  do i=4,nx
     do j=4,ny
        do k=4,4
           if( r(i,j,k)-4.0d0 .lt. 96.0d0 .and. &
                r(i,j,k)+4.0d0 .gt. 96.0d0) then
              write(6,"(4i4,1P10E15.6)") GRHydro_reflevel, i,j,k,&
                   r(i,j,k),rho(i,j,k),eps(i,j,k),y_e(i,j,k),temperature(i,j,k)
           endif
        enddo
     enddo
  enddo

!  call flush(6)
!  if(GRHydro_reflevel.eq.4.and.temperature(1,1,1).lt.0.1d0) then
!     call CCTK_EROR("stop")
!     STOP
!  endif

end subroutine GRHydro_Debug


subroutine GRHydro_RefinementLevel(CCTK_ARGUMENTS)

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS

  GRHydro_reflevel = int(log10(dble(cctk_levfac(1)))/log10(2.0d0))

end subroutine GRHydro_RefinementLevel

subroutine GRHydro_SqrtSpatialDeterminant(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  integer i,j,k
  integer nx, ny, nz

  ! save memory when MP is not used
  CCTK_INT :: GRHydro_UseGeneralCoordinates
  CCTK_REAL, DIMENSION(cctk_ash1,cctk_ash2,cctk_ash3) :: g11, g12, g13, g22, g23, g33
  pointer (pg11,g11), (pg12,g12), (pg13,g13), (pg22,g22), (pg23,g23), (pg33,g33)

  ! save memory when MP is not used
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

  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)

  !$OMP PARALLEL DO PRIVATE(i,j,k)
  do k=1,nz
     do j=1,ny
        do i=1,nx
           sdetg(i,j,k) = -(g13(i,j,k)**2)*g22(i,j,k) + & 
                2.0d0*g12(i,j,k)*g13(i,j,k)*g23(i,j,k) - &
                g11(i,j,k)*(g23(i,j,k)**2) - &
                (g12(i,j,k)**2)*g33(i,j,k)  + &
                (g11(i,j,k)*g22(i,j,k))*g33(i,j,k)
           sdetg(i,j,k) = sqrt(sdetg(i,j,k))
        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO

end subroutine GRHydro_SqrtSpatialDeterminant


subroutine SpatialDeterminant(gxx,gxy,gxz,gyy,gyz,gzz,det)

  implicit none
  
  CCTK_REAL :: det
  CCTK_REAL :: gxx,gxy,gxz,gyy,gyz,gzz

  det = -gxz**2*gyy + 2*gxy*gxz*gyz - gxx*gyz**2 - gxy**2*gzz + gxx*gyy*gzz


!!$  Why is this weird order necessary? Search me. It just seemed 
!!$  to make a really odd NaN go away.

!!$  det =  -gxz**2*gyy + 2*gxy*gxz*gyz 
!!$  det = det - gxx*gyz**2 - gxy**2*gzz
!!$  det = det + gxx*gyy*gzz
  
  return
  
end subroutine SpatialDeterminant



/*@@
@routine    UpperMetric
@date       Sat Jan 26 02:32:26 2002
@author     
@desc 
Calculates the upper metric. The determinant is given, not
calculated.
@enddesc 
@calls     
@calledby   
@history 

@endhistory 
@@*/

subroutine UpperMetric(uxx, uxy, uxz, uyy, uyz, uzz, &
     det, gxx, gxy, gxz, gyy, gyz, gzz)
  
  implicit none
  
  CCTK_REAL :: uxx, uxy, uxz, uyy, uyz, uzz, det, &
       gxx, gxy, gxz, gyy, gyz, gzz, invdet
  
  invdet = 1.d0 / det
  uxx=(-gyz**2 + gyy*gzz)*invdet
  uxy=(gxz*gyz - gxy*gzz)*invdet
  uyy=(-gxz**2 + gxx*gzz)*invdet
  uxz=(-gxz*gyy + gxy*gyz)*invdet
  uyz=(gxy*gxz - gxx*gyz)*invdet
  uzz=(-gxy**2 + gxx*gyy)*invdet
  
  return
  
end subroutine UpperMetric


