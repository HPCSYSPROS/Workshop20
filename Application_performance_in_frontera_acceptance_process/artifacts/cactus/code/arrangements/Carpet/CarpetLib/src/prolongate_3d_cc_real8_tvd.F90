#include "cctk.h"

! This routine performs "TVD" prolongation. It is intended to be used 
! with GFs that are not expected to be smooth, particularly those
! that must also obey certain constraints. The obvious example is the 
! density in hydrodynamics, which may be discontinuous yet must be
! strictly positive.
!
! To ensure that this prolongation method is used you should add the
! tag
!
!      tags='Prolongation="TVD"'
!
! to the interface.ccl on the appropriate group.
!
! This applies minmod type limiting to the slope, checking over the
! entire coarse grid cell for the minimum modulus in each direction.



! Grid point locations and their indices:
!
! global   0   4      12      20      28   |
! local    |   0       1       2       3   |
!          |   C       C       C       C   |
!          | f   f   f   f   f   f   f   f |
! local    | 0   1   2   3   4   5   6   7 |
! global   0 2   6  10  14  18  22  26  30 |
!
! Interpolation with even interpolation order:
!    offset zero (fine index = 2 * coarse index):
!       [1]
!    offset one (fine index = 2 * coarse index + 1):
!       [1]
! Interpolation with odd interpolation order:
!    offset zero:
!       [1 1 0]/2
!    offset one:
!       [0 1 1]/2
! The centres of these stencils are located at the coarse grid point
! corresponding to the fine grid point (the interpolation target)
! minus the offset. Example: fine grid 8 -> coarse grid 8, fine grid
! 12 -> also coarse grid 8.



#define CHKIDX(i,j,k, imax,jmax,kmax, where) \
   if ((i)<1 .or. (i)>(imax) .or. (j)<1 .or. (j)>(jmax) .or. (k)<1 .or. (k)>(kmax)) then &&\
      write (msg, '(a, " array index out of bounds: shape is (",i4,",",i4,",",i4,"), index is (",i4,",",i4,",",i4,")")') (where), (imax),(jmax),(kmax), (i),(j),(k) &&\
      call CCTK_WARN (CCTK_WARN_ABORT, msg) &&\
   end if



subroutine prolongate_3d_cc_real8_tvd ( &
     src, srcipadext, srcjpadext, srckpadext, srciext, srcjext, srckext, &
     dst, dstipadext, dstjpadext, dstkpadext, dstiext, dstjext, dstkext, &
     srcbbox, dstbbox, regbbox)

  implicit none

  CCTK_REAL8, parameter :: one=1, fourth=one/4

  integer srcipadext, srcjpadext, srckpadext
  integer srciext, srcjext, srckext
  CCTK_REAL8 src(srcipadext,srcjpadext,srckpadext)
  integer dstipadext, dstjpadext, dstkpadext
  integer dstiext, dstjext, dstkext
  CCTK_REAL8 dst(dstipadext,dstjpadext,dstkpadext)
  ! bbox(:,1) is lower boundary (inclusive)
  ! bbox(:,2) is upper boundary (inclusive)
  ! bbox(:,3) is stride
  integer srcbbox(3,3), dstbbox(3,3), regbbox(3,3)

  character*1000 msg

  integer offsetlo, offsethi

  integer regiext, regjext, regkext

  integer dstifac, dstjfac, dstkfac

  integer srcioff, srcjoff, srckoff
  integer dstioff, dstjoff, dstkoff

  integer i, j, k
  integer i0, j0, k0
  integer fi, fj, fk
  integer d

  CCTK_REAL8 :: dlo, dup
  CCTK_REAL8 :: slopex, slopey, slopez

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
!!$     if (mod(srcbbox(d,1), srcbbox(d,3)).ne.0 &
!!$          .or. mod(dstbbox(d,1), dstbbox(d,3)).ne.0 &
!!$          .or. mod(regbbox(d,1), regbbox(d,3)).ne.0) then
!!$        call CCTK_WARN (0, "Internal error: array origins are not integer multiples of the strides")
!!$     end if
     if (regbbox(d,1).gt.regbbox(d,2)) then
        ! This could be handled, but is likely to point to an error elsewhere
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

  if (any(srcbbox(:,3) / dstbbox(:,3) /= 2)) then
     call CCTK_WARN (0, "Internal error: refinement factor is not 2")
  end if

  if (any(mod(regbbox(:,3), 2) /= 0)) then
     call CCTK_WARN (0, "Internal error: region stride is not a multiple of 2")
  end if

  regiext = (regbbox(1,2) - regbbox(1,1)) / regbbox(1,3) + 1
  regjext = (regbbox(2,2) - regbbox(2,1)) / regbbox(2,3) + 1
  regkext = (regbbox(3,2) - regbbox(3,1)) / regbbox(3,3) + 1

  dstifac = 2 ! srcbbox(1,3) / dstbbox(1,3)
  dstjfac = 2 ! srcbbox(2,3) / dstbbox(2,3)
  dstkfac = 2 ! srcbbox(3,3) / dstbbox(3,3)

  srcioff = (regbbox(1,1) - srcbbox(1,1) + regbbox(1,3) / 2) / dstbbox(1,3)
  srcjoff = (regbbox(2,1) - srcbbox(2,1) + regbbox(1,3) / 2) / dstbbox(2,3)
  srckoff = (regbbox(3,1) - srcbbox(3,1) + regbbox(1,3) / 2) / dstbbox(3,3)

  dstioff = (regbbox(1,1) - dstbbox(1,1)) / dstbbox(1,3)
  dstjoff = (regbbox(2,1) - dstbbox(2,1)) / dstbbox(2,3)
  dstkoff = (regbbox(3,1) - dstbbox(3,1)) / dstbbox(3,3)

  if (srcioff<0 .or. srcjoff<0 .or. srckoff<0) then
     call CCTK_WARN (0, "Internal error: source array offset is negative")
  end if

  if (dstioff<0 .or. dstjoff<0 .or. dstkoff<0) then
     call CCTK_WARN (0, "Internal error: destination array offset is negative")
  end if

  ! Loop over fine region

  do k = 0, regkext-1
     k0 = (srckoff + k) / dstkfac
     fk = mod(srckoff + k, dstkfac)

     do j = 0, regjext-1
        j0 = (srcjoff + j) / dstjfac
        fj = mod(srcjoff + j, dstjfac)

        do i = 0, regiext-1
           i0 = (srcioff + i) / dstifac
           fi = mod(srcioff + i, dstifac)

           ! We consider the nearest coarse grid point. From this
           ! point, we use a slope to interpolate (or extrapolate, as
           ! it may be) towards the interpolation point. We consider
           ! both the left and the right slope starting from this
           ! coarse grid point, and choose the minmod of these. This
           ! is, loosely speaking, either the smaller of these slopes,
           ! or zero if they differ too much.

           CHKIDX (i0  ,j0  ,k0  , srciext,srcjext,srckext, "src")
           CHKIDX (i0+2,j0+2,k0+2, srciext,srcjext,srckext, "src")

           dlo = src(i0+1,j0+1,k0+1) - src(i0+0,j0+1,k0+1)
           dup = src(i0+2,j0+1,k0+1) - src(i0+1,j0+1,k0+1)
           slopex = (2*fi-1) * fourth * minmod(dlo,dup)
           
           dlo = src(i0+1,j0+1,k0+1) - src(i0+1,j0+0,k0+1)
           dup = src(i0+1,j0+2,k0+1) - src(i0+1,j0+1,k0+1)
           slopey = (2*fj-1) * fourth * minmod(dlo,dup)
           
           dlo = src(i0+1,j0+1,k0+1) - src(i0+1,j0+1,k0+0)
           dup = src(i0+1,j0+1,k0+2) - src(i0+1,j0+1,k0+1)
           slopez = (2*fk-1) * fourth * minmod(dlo,dup)
           
           CHKIDX (dstioff+i+1, dstjoff+j+1, dstkoff+k+1, dstiext,dstjext,dstkext, "dst")
           
           dst (dstioff+i+1, dstjoff+j+1, dstkoff+k+1) = &
                src(i0+1,j0+1,k0+1) + slopex + slopey + slopez
           
        end do
     end do
  end do
  
contains
  
  function minmod(a, b)
    
    CCTK_REAL8 minmod
    CCTK_REAL8 a, b
    
    if (a * b < 0) then
       minmod = 0
    else if (abs(a) < abs(b)) then
       minmod = a
    else
       minmod = b
    end if
    
  end function minmod
  
end subroutine prolongate_3d_cc_real8_tvd
