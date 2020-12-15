#ifndef OMIT_F90
!!$     -*-Fortran-*-

#include "cctk.h"


!!$ This routine performs "WENO" prolongation. It is intended to be used 
!!$ with GFs that are not expected to be smooth, particularly those
!!$ that must also obey certain constraints. The obvious example is the 
!!$ density in hydrodynamics, which may be discontinuous yet must be
!!$ strictly positive.
!!$
!!$ To ensure that this prolongation method is used you should add the
!!$ tag
!!$
!!$      tags='Prolongation="WENO"'
!!$
!!$ to the interface.ccl on the appropriate group.
!!$
!!$ This applies WENO2 type limiting to the slope, checking over the
!!$ entire coarse grid cell for the least oscillatory quadratic in each 
!!$ direction. If the slope changes sign over the extrema, linear
!!$ interpolation is used instead.


function weno1d(q)
  
  implicit none

  CCTK_REAL8 :: weno1d
  CCTK_REAL8 :: q(5)
  CCTK_REAL8 :: zero, one, two, three, four, five, eight, ten, eleven, &
      thirteen, fifteen, nineteen, twentyfive, thirtyone, epsilon
  parameter (zero       = 0)
  parameter (one        = 1)
  parameter (two        = 2)
  parameter (three      = 3)
  parameter (four       = 4)
  parameter (five       = 5)
  parameter (eight      = 8)
  parameter (ten        = 10)
  parameter (eleven     = 11)
  parameter (thirteen   = 13)
  parameter (fifteen    = 15)
  parameter (nineteen   = 19)
  parameter (twentyfive = 25)
  parameter (thirtyone  = 31)
  parameter (epsilon    = 0.000000001)

  CCTK_REAL8, dimension(0:2,0:2) :: c
  CCTK_REAL8 :: wtildesum, qmax, qmin
  CCTK_REAL8, dimension(0:2) :: d, beta, wtilde, w, vr
  logical, dimension(0:2) :: hacked

  integer :: j, k

!!$  Linear weights

  d(0) = three / ten
  d(1) = three / five
  d(2) = one   / ten

  c(0,0) = three / eight
  c(0,1) = three / four
  c(0,2) = -one / eight
  c(1,0) = -one / eight
  c(1,1) = three / four
  c(1,2) = three / eight
  c(2,0) = three / eight
  c(2,1) = -five / four
  c(2,2) = fifteen / eight

!!$  Substencils

  do j = 0, 2
    vr(j) = 0.d0
    do k = 0, 2
      vr(j) = vr(j) + c(j, k) * q(3 - j + k)
    end do
  end do
  
!!$  Nonlinear weights

  beta(0) = (ten * q(3)**2 - &
       thirtyone * q(3) * q(4) + &
       twentyfive * q(4)**2 + &
       eleven * q(3) * q(5) - &
       nineteen * q(4) * q(5) + &
       four * q(5)**2) / three
  beta(1) = (four * q(2)**2 - &
       thirteen * q(2) * q(3) + &
       thirteen * q(3)**2 + & 
       five * q(2) * q(4) - & 
       thirteen * q(3) * q(4) + &
       four * q(4)**2) / three
  beta(2) = (four * q(1)**2 - &
       nineteen * q(1) * q(2) + &
       twentyfive * q(2)**2 + &
       eleven * q(1) * q(3) - &
       thirtyone * q(2) * q(3) + & 
       ten * q(3)**2) / three

  do j = 0, 2
    wtilde(j) = d(j) / (epsilon + beta(j))**2
  end do

!!$  Hack the weights if outside the range

  qmax = maxval(q)
  qmin = minval(q)

  do j = 0, 2
    hacked(j) = .false.
    if ( (qmax - vr(j)) * (vr(j) - qmin) < 0.d0 ) then
      wtilde(j) = 0.d0
      hacked(j) = .true.
    end if
  end do

!!$  If all weights were hacked we cannot get a good interpolant;
!!$  drop to linear interpolation
  
  if (hacked(0).and.hacked(1).and.hacked(2)) then

!!$    Linear interpolation
    weno1d = 0.5d0 * (q(3) + q(4))

  else
  
    wtildesum = wtilde(0) + wtilde(1) + wtilde(2)
    w = wtilde / wtildesum

    weno1d = 0.d0

    do j = 0, 2
      weno1d = weno1d + w(j) * vr(j)
    end do

  end if

!!$  if (.not.( (weno1d .ge. 0.d0 ).or.(weno1d .le. 0.d0) ) ) then
!!$    write(*,*) 'Error?', weno1d
!!$    write(*,*) 'Done weno1d', weno1d, hacked
!!$    write(*,*) 'Substencil',vr
!!$    write(*,*) 'Weights', w
!!$    write(*,*) 'Indicators', beta
!!$    write(*,*) 'Input', q
!!$  end if
  
!!$  write(*,*) 'Done weno1d', weno1d, hacked
!!$  write(*,*) 'Substencil',vr
!!$  write(*,*) 'Weights', w
!!$  write(*,*) 'Indicators', beta
!!$  write(*,*) 'Input', q
  
end function weno1d

subroutine prolongate_3d_real8_weno ( &
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

  CCTK_REAL8, dimension(0:4,0:4) :: tmp1
  CCTK_REAL8, dimension(0:4) :: tmp2

  external weno1d
  CCTK_REAL8 weno1d

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
               weno1d(src(i0-1:i0+3,j0+1,k0+1))

        case (10)
!!$          Interpolate only in y

          dst (dstioff+i+1, dstjoff+j+1, dstkoff+k+1) = &
               weno1d(src(i0+1,j0-1:j0+3,k0+1))

        case (11)
!!$          Interpolate only in x and y

          do jj = 0, 4
            tmp2(jj) = weno1d(src(i0-1:i0+3,j0+jj-1,k0+1))
          end do

          dst (dstioff+i+1, dstjoff+j+1, dstkoff+k+1) = &
               weno1d(tmp2(0:4))

        case (100)
!!$          Interpolate only in z

          dst (dstioff+i+1, dstjoff+j+1, dstkoff+k+1) = &
               weno1d(src(i0+1,j0+1,k0-1:k0+3))

        case (101)
!!$          Interpolate only in x and z

          do kk = 0, 4
            tmp2(kk) = weno1d(src(i0-1:i0+3,j0+1,k0+kk-1))
          end do

          dst (dstioff+i+1, dstjoff+j+1, dstkoff+k+1) = &
               weno1d(tmp2(0:4))

        case (110)
!!$          Interpolate only in y and z

          do kk = 0, 4
            tmp2(kk) = weno1d(src(i0+1,j0-1:j0+3,k0+kk-1))
          end do

          dst (dstioff+i+1, dstjoff+j+1, dstkoff+k+1) = &
               weno1d(tmp2(0:4))

        case (111)
!!$          Interpolate in all of x, y, and z

          do jj = 0, 4
            do kk = 0, 4
              tmp1(jj,kk) = weno1d(src(i0-1:i0+3,j0+jj-1,k0+kk-1))
            end do
          end do
          do ii = 0, 4
            tmp2(ii) = weno1d(tmp1(0:4,ii))
          end do
          
          dst (dstioff+i+1, dstjoff+j+1, dstkoff+k+1) = &
               weno1d(tmp2(0:4))

        case default
          call CCTK_ERROR("Internal error in WENO prolongation. Should only be used with refinement factor 2!")
        end select

      end do
    end do
  end do

end subroutine prolongate_3d_real8_weno
#endif	/* !OMIT_F90 */
