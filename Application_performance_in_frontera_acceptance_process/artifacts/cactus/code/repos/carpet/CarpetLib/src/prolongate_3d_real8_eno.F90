#include "cctk.h"


!!$ This routine performs "ENO" prolongation. It is intended to be used 
!!$ with GFs that are not expected to be smooth, particularly those
!!$ that must also obey certain constraints. The obvious example is the 
!!$ density in hydrodynamics, which may be discontinuous yet must be
!!$ strictly positive.
!!$
!!$ To ensure that this prolongation method is used you should add the
!!$ tag
!!$
!!$      tags='Prolongation="ENO"'
!!$
!!$ to the interface.ccl on the appropriate group.
!!$
!!$ This applies ENO2 type limiting to the slope, checking over the
!!$ entire coarse grid cell for the least oscillatory quadratic in each 
!!$ direction. If the slope changes sign over the extrema, linear
!!$ interpolation is used instead.


function eno1d(q)
  
  implicit none

  CCTK_REAL8 :: eno1d
  CCTK_REAL8, INTENT(IN) :: q(4)
  CCTK_REAL8 :: zero, one, two, three, six, half, eighth
  parameter (zero = 0)
  parameter (two = 2)
  parameter (one = 1)
  parameter (three = 3)
  parameter (six = 6)
  parameter (eighth = one / 8)
  parameter (half = one / two)
  CCTK_REAL8 :: diffleft, diffright

!!$  Directly find the second undivided differences
!!$  We need to pick between discrete values at
!!$  1 2 3 4 for the interpolation between 2 and 3.

  diffleft  = q(1) + q(3) - two * q(2)
  diffright = q(2) + q(4) - two * q(3)

  if ( abs(diffleft) .lt. abs(diffright) ) then

!!$      Apply the left quadratic

    eno1d = eighth * (-q(1) + six * q(2) + three * q(3))
    
  else

!!$      Apply the right quadratic

    eno1d = eighth * (three * q(2) + six * q(3) - q(4))
    
  end if

!!$    Check that the quadratic is reasonable:
!!$    Check 1: interpolated value between
!!$             values at interpolation points
!!$    Check 2: sign of the curvature of the interpolating 
!!$             polynomial does not change.

  if ( ((eno1d-q(2)) * (q(3)-eno1d) .lt. zero) &
       .or. &
       (diffleft*diffright .le. zero) ) then
      
!!$      Not reasonable. Linear interpolation

    eno1d = half * (q(2) + q(3))

  end if
  
end function eno1d

subroutine prolongate_3d_real8_eno ( &
     src, srcipadext, srcjpadext, srckpadext, srciext, srcjext, srckext, &
     dst, dstipadext, dstjpadext, dstkpadext, dstiext, dstjext, dstkext, &
     srcbbox, dstbbox, regbbox)

  implicit none

  CCTK_REAL8 one
  parameter (one = 1)

  integer srcipadext, srcjpadext, srckpadext
  integer srciext, srcjext, srckext
  CCTK_REAL8 src(srcipadext,srcjpadext,srckpadext)
  integer dstipadext, dstjpadext, dstkpadext
  integer dstiext, dstjext, dstkext
  CCTK_REAL8 dst(dstipadext,dstjpadext,dstkpadext)
!!$     bbox(:,1) is lower boundary (inclusive)
!!$     bbox(:,2) is upper boundary (inclusive)
!!$     bbox(:,3) is stride
  integer srcbbox(3,3), dstbbox(3,3), regbbox(3,3)

  integer offsetlo, offsethi

  integer regiext, regjext, regkext

  integer dstifac, dstjfac, dstkfac

  integer srcioff, srcjoff, srckoff
  integer dstioff, dstjoff, dstkoff

  integer i, j, k
  integer i0, j0, k0
  integer fi, fj, fk
  integer ii, jj, kk
  integer d

  CCTK_REAL8, dimension(0:3,0:3) :: tmp1
  CCTK_REAL8, dimension(0:3) :: tmp2

  external eno1d
  CCTK_REAL8 eno1d

  CCTK_REAL8 half, zero
  parameter (half = 0.5)
  parameter (zero = 0)

  do d=1,3
    if (srcbbox(d,3).eq.0 .or. dstbbox(d,3).eq.0 &
         .or. regbbox(d,3).eq.0) then
      call CCTK_WARN (0, "Internal error: stride is zero")
    end if
    if (srcbbox(d,3).le.regbbox(d,3) &
         .or. dstbbox(d,3).ne.regbbox(d,3)) then
      call CCTK_WARN (0, "Internal error: strides disagree")
    end if
    if (mod(srcbbox(d,3), dstbbox(d,3)).ne.0) then
      call CCTK_WARN (0, "Internal error: destination strides are not integer multiples of the source strides")
    end if
    if (mod(srcbbox(d,1), srcbbox(d,3)).ne.0 &
         .or. mod(dstbbox(d,1), dstbbox(d,3)).ne.0 &
         .or. mod(regbbox(d,1), regbbox(d,3)).ne.0) then
      call CCTK_WARN (0, "Internal error: array origins are not integer multiples of the strides")
    end if
    if (regbbox(d,1).gt.regbbox(d,2)) then
!!$     This could be handled, but is likely to point to an error elsewhere
      call CCTK_WARN (0, "Internal error: region extent is empty")
    end if
    regkext = (regbbox(d,2) - regbbox(d,1)) / regbbox(d,3) + 1
    dstkfac = srcbbox(d,3) / dstbbox(d,3)
    srckoff = (regbbox(d,1) - srcbbox(d,1)) / dstbbox(d,3)
    offsetlo = regbbox(d,3)
    if (mod(srckoff + 0, dstkfac).eq.0) then
      offsetlo = 0
      if (regkext.gt.1) then
        offsetlo = regbbox(d,3)
      end if
    end if
    offsethi = regbbox(d,3)
    if (mod(srckoff + regkext-1, dstkfac).eq.0) then
      offsethi = 0
      if (regkext.gt.1) then
        offsethi = regbbox(d,3)
      end if
    end if
    if (regbbox(d,1)-offsetlo.lt.srcbbox(d,1) &
         .or. regbbox(d,2)+offsethi.gt.srcbbox(d,2) &
         .or. regbbox(d,1).lt.dstbbox(d,1) &
         .or. regbbox(d,2).gt.dstbbox(d,2)) then
      call CCTK_WARN (0, "Internal error: region extent is not contained in array extent")
    end if
  end do

  if (srciext.ne.(srcbbox(1,2)-srcbbox(1,1))/srcbbox(1,3)+1 &
       .or. srcjext.ne.(srcbbox(2,2)-srcbbox(2,1))/srcbbox(2,3)+1 &
       .or. srckext.ne.(srcbbox(3,2)-srcbbox(3,1))/srcbbox(3,3)+1 &
       .or. dstiext.ne.(dstbbox(1,2)-dstbbox(1,1))/dstbbox(1,3)+1 &
       .or. dstjext.ne.(dstbbox(2,2)-dstbbox(2,1))/dstbbox(2,3)+1 &
       .or. dstkext.ne.(dstbbox(3,2)-dstbbox(3,1))/dstbbox(3,3)+1) then
    call CCTK_WARN (0, "Internal error: array sizes don't agree with bounding boxes")
  end if

  regiext = (regbbox(1,2) - regbbox(1,1)) / regbbox(1,3) + 1
  regjext = (regbbox(2,2) - regbbox(2,1)) / regbbox(2,3) + 1
  regkext = (regbbox(3,2) - regbbox(3,1)) / regbbox(3,3) + 1

  dstifac = srcbbox(1,3) / dstbbox(1,3)
  dstjfac = srcbbox(2,3) / dstbbox(2,3)
  dstkfac = srcbbox(3,3) / dstbbox(3,3)

  srcioff = (regbbox(1,1) - srcbbox(1,1)) / dstbbox(1,3)
  srcjoff = (regbbox(2,1) - srcbbox(2,1)) / dstbbox(2,3)
  srckoff = (regbbox(3,1) - srcbbox(3,1)) / dstbbox(3,3)

  dstioff = (regbbox(1,1) - dstbbox(1,1)) / dstbbox(1,3)
  dstjoff = (regbbox(2,1) - dstbbox(2,1)) / dstbbox(2,3)
  dstkoff = (regbbox(3,1) - dstbbox(3,1)) / dstbbox(3,3)

!!$     Loop over fine region

  !$omp parallel do collapse(3) private(i,j,k, i0,fi,j0,fj,k0,fk, tmp1,tmp2, ii,jj,kk)
  do k = 0, regkext-1
    do j = 0, regjext-1
      do i = 0, regiext-1
         
        i0 = (srcioff + i) / dstifac
        fi = mod(srcioff + i, dstifac)
        j0 = (srcjoff + j) / dstjfac
        fj = mod(srcjoff + j, dstjfac)
        k0 = (srckoff + k) / dstkfac
        fk = mod(srckoff + k, dstkfac)

!!$        Where is the fine grid point w.r.t the coarse grid?

        select case (fi + 10*fj + 100*fk)
        case (0)
!!$            On a coarse grid point exactly!

          dst (dstioff+i+1, dstjoff+j+1, dstkoff+k+1) = &
               src(i0+1,j0+1,k0+1)

        case (1)
!!$          Interpolate only in x

          dst (dstioff+i+1, dstjoff+j+1, dstkoff+k+1) = &
               eno1d(src(i0:i0+3,j0+1,k0+1))

        case (10)
!!$          Interpolate only in y

          dst (dstioff+i+1, dstjoff+j+1, dstkoff+k+1) = &
               eno1d(src(i0+1,j0:j0+3,k0+1))

        case (11)
!!$          Interpolate only in x and y

          do jj = 0, 3
            tmp2(jj) = eno1d(src(i0:i0+3,j0+jj,k0+1))
          end do

          dst (dstioff+i+1, dstjoff+j+1, dstkoff+k+1) = &
               eno1d(tmp2(0:3))

        case (100)
!!$          Interpolate only in z

          dst (dstioff+i+1, dstjoff+j+1, dstkoff+k+1) = &
               eno1d(src(i0+1,j0+1,k0:k0+3))

        case (101)
!!$          Interpolate only in x and z

          do kk = 0, 3
            tmp2(kk) = eno1d(src(i0:i0+3,j0+1,k0+kk))
          end do

          dst (dstioff+i+1, dstjoff+j+1, dstkoff+k+1) = &
               eno1d(tmp2(0:3))

        case (110)
!!$          Interpolate only in y and z

          do kk = 0, 3
            tmp2(kk) = eno1d(src(i0+1,j0:j0+3,k0+kk))
          end do

          dst (dstioff+i+1, dstjoff+j+1, dstkoff+k+1) = &
               eno1d(tmp2(0:3))

        case (111)
!!$          Interpolate in all of x, y, and z

          do jj = 0, 3
            do kk = 0, 3
              tmp1(jj,kk) = eno1d(src(i0:i0+3,j0+jj,k0+kk))
            end do
          end do
          do ii = 0, 3
            tmp2(ii) = eno1d(tmp1(0:3,ii))
          end do
          
          dst (dstioff+i+1, dstjoff+j+1, dstkoff+k+1) = &
               eno1d(tmp2(0:3))

        case default
          call CCTK_ERROR("Internal error in ENO prolongation. Should only be used with refinement factor 2!")
        end select

      end do
    end do
  end do

end subroutine prolongate_3d_real8_eno
