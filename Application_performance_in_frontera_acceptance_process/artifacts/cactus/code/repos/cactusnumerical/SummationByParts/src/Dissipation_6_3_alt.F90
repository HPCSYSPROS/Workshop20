! $Header$

#include "cctk.h"
#include "cctk_Parameters.h"

subroutine dissipation_6_3_alt (var, ni, nj, nk, bb, gsize, offset, delta, epsilon, rhs)

  implicit none

  DECLARE_CCTK_PARAMETERS

  integer ::  ni, nj, nk
  CCTK_REAL, dimension(ni,nj,nk), intent(in) :: var
  CCTK_REAL, dimension(ni,nj,nk), intent(inout) :: rhs
  CCTK_INT, dimension(6), intent(in) :: bb
  CCTK_INT, dimension(3), intent(in) :: gsize
  CCTK_INT, dimension(6), intent(in) :: offset
  CCTK_REAL, dimension(3), intent(in) :: delta
  CCTK_REAL, intent(in) :: epsilon 

  CCTK_REAL :: zero = 0.0
  integer, parameter :: wp = kind(zero)
  CCTK_REAL, dimension(10,6) :: a
  CCTK_REAL :: idel

  CCTK_INT :: il, ir, jl, jr, kl, kr, ol, or

  call set_coeff ( a )

  idel = epsilon / ( 256 * delta(1) )

  if ( bb(1) == 0 ) then
    il = 1 + gsize(1)
  else
    ol = offset(1)

!$OMP PARALLEL WORKSHARE
    rhs(1+ol,:,:) = rhs(1+ol,:,:) + &
                 ( a(1,1) * var(1+ol,:,:) + a(2,1) * var(2+ol,:,:) + &
                   a(3,1) * var(3+ol,:,:) + a(4,1) * var(4+ol,:,:) + &
                   a(5,1) * var(5+ol,:,:) ) * idel
    rhs(2+ol,:,:) = rhs(2+ol,:,:) + &
                 ( a(1,2) * var(1+ol,:,:) + a(2,2) * var(2+ol,:,:) + &
                   a(3,2) * var(3+ol,:,:) + a(4,2) * var(4+ol,:,:) + &
                   a(5,2) * var(5+ol,:,:) + a(6,2) * var(6+ol,:,:) ) * idel
    rhs(3+ol,:,:) = rhs(3+ol,:,:) + &
                 ( a(1,3) * var(1+ol,:,:) + a(2,3) * var(2+ol,:,:) + &
                   a(3,3) * var(3+ol,:,:) + a(4,3) * var(4+ol,:,:) + &
                   a(5,3) * var(5+ol,:,:) + a(6,3) * var(6+ol,:,:) + &
                   a(7,3) * var(7+ol,:,:) ) * idel
    rhs(4+ol,:,:) = rhs(4+ol,:,:) + &
                 ( a(1,4) * var(1+ol,:,:) + a(2,4) * var(2+ol,:,:) + &
                   a(3,4) * var(3+ol,:,:) + a(4,4) * var(4+ol,:,:) + &
                   a(5,4) * var(5+ol,:,:) + a(6,4) * var(6+ol,:,:) + &
                   a(7,4) * var(7+ol,:,:) + a(8,4) * var(8+ol,:,:) ) * idel
    rhs(5+ol,:,:) = rhs(5+ol,:,:) + &
                 ( a(1,5) * var(1+ol,:,:) + a(2,5) * var(2+ol,:,:) + &
                   a(3,5) * var(3+ol,:,:) + a(4,5) * var(4+ol,:,:) + &
                   a(5,5) * var(5+ol,:,:) + a(6,5) * var(6+ol,:,:) + &
                   a(7,5) * var(7+ol,:,:) + a(8,5) * var(8+ol,:,:) + &
                   a(9,5) * var(9+ol,:,:) ) * idel
    rhs(6+ol,:,:) = rhs(6+ol,:,:) + &
                 ( a(2,6) * var(2+ol,:,:) + a(3,6) * var(3+ol,:,:) + &
                   a(4,6) * var(4+ol,:,:) + a(5,6) * var(5+ol,:,:) + &
                   a(6,6) * var(6+ol,:,:) + a(7,6) * var(7+ol,:,:) + &
                   a(8,6) * var(8+ol,:,:) + a(9,6) * var(9+ol,:,:) + &
                   a(10,6) * var(10+ol,:,:) ) * idel
!$OMP END PARALLEL WORKSHARE

    il = 7 + ol
  end if
  if ( bb(2) == 0 ) then
    ir = ni - gsize(1)
  else
    or = ni - offset(2)

!$OMP PARALLEL WORKSHARE
    rhs(or-5,:,:) = rhs(or-5,:,:) + &
                 ( a(2,6) * var(or-1,:,:) + a(3,6) * var(or-2,:,:) + &
                   a(4,6) * var(or-3,:,:) + a(5,6) * var(or-4,:,:) + &
                   a(6,6) * var(or-5,:,:) + a(7,6) * var(or-6,:,:) + &
                   a(8,6) * var(or-7,:,:) + a(9,6) * var(or-8,:,:) + &
                   a(10,6) * var(or-9,:,:) ) * idel
    rhs(or-4,:,:) = rhs(or-4,:,:) + &
                 ( a(1,5) * var(or,:,:) + a(2,5) * var(or-1,:,:) + &
                   a(3,5) * var(or-2,:,:) + a(4,5) * var(or-3,:,:) + &
                   a(5,5) * var(or-4,:,:) + a(6,5) * var(or-5,:,:) + &
                   a(7,5) * var(or-6,:,:) + a(8,5) * var(or-7,:,:) + &
                   a(9,5) * var(or-8,:,:) ) * idel
    rhs(or-3,:,:) = rhs(or-3,:,:) + &
                 ( a(1,4) * var(or,:,:) + a(2,4) * var(or-1,:,:) + &
                   a(3,4) * var(or-2,:,:) + a(4,4) * var(or-3,:,:) + &
                   a(5,4) * var(or-4,:,:) + a(6,4) * var(or-5,:,:) + &
                   a(7,4) * var(or-6,:,:) + a(8,4) * var(or-7,:,:) ) * idel
    rhs(or-2,:,:) = rhs(or-2,:,:) + &
                 ( a(1,3) * var(or,:,:) + a(2,3) * var(or-1,:,:) + &
                   a(3,3) * var(or-2,:,:) + a(4,3) * var(or-3,:,:) + &
                   a(5,3) * var(or-4,:,:) + a(6,3) * var(or-5,:,:) + &
                   a(7,3) * var(or-6,:,:) ) * idel
    rhs(or-1,:,:) = rhs(or-1,:,:) + &
                 ( a(1,2) * var(or,:,:) + a(2,2) * var(or-1,:,:) + &
                   a(3,2) * var(or-2,:,:) + a(4,2) * var(or-3,:,:) + &
                   a(5,2) * var(or-4,:,:) + a(6,2) * var(or-5,:,:) ) * idel
    rhs(or,:,:)   = rhs(or,:,:) + &
                 ( a(1,1) * var(or,:,:) + a(2,1) * var(or-1,:,:) + &
                   a(3,1) * var(or-2,:,:) + a(4,1) * var(or-3,:,:) + &
                   a(5,1) * var(or-4,:,:) ) * idel
!$OMP END PARALLEL WORKSHARE

    ir = or - 6
  end if

!$OMP PARALLEL WORKSHARE
  rhs(il:ir,:,:) = rhs(il:ir,:,:) + &
                   ( -70.0_wp * var(il:ir,:,:) + &
                      56.0_wp * ( var(il-1:ir-1,:,:) + &
                                  var(il+1:ir+1,:,:) ) - &
                      28.0_wp * ( var(il-2:ir-2,:,:) + &
                                  var(il+2:ir+2,:,:) ) + &
                       8.0_wp * ( var(il-3:ir-3,:,:) + &
                                  var(il+3:ir+3,:,:) ) - &
                                ( var(il-4:ir-4,:,:) + &
                                  var(il+4:ir+4,:,:) ) ) * idel
!$OMP END PARALLEL WORKSHARE

  if ( zero_derivs_y == 0 ) then
    call set_coeff ( a )
  
    idel = epsilon / ( 256 * delta(2) )
  
    if ( bb(3) == 0 ) then
      jl = 1 + gsize(2)
    else
      ol = offset(3)

!$OMP PARALLEL WORKSHARE
      rhs(:,1+ol,:) = rhs(:,1+ol,:) + &
                   ( a(1,1) * var(:,1+ol,:) + a(2,1) * var(:,2+ol,:) + &
                     a(3,1) * var(:,3+ol,:) + a(4,1) * var(:,4+ol,:) + &
                     a(5,1) * var(:,5+ol,:) ) * idel
      rhs(:,2+ol,:) = rhs(:,2+ol,:) + &
                   ( a(1,2) * var(:,1+ol,:) + a(2,2) * var(:,2+ol,:) + &
                     a(3,2) * var(:,3+ol,:) + a(4,2) * var(:,4+ol,:) + &
                     a(5,2) * var(:,5+ol,:) + a(6,2) * var(:,6+ol,:) ) * idel
      rhs(:,3+ol,:) = rhs(:,3+ol,:) + &
                   ( a(1,3) * var(:,1+ol,:) + a(2,3) * var(:,2+ol,:) + &
                     a(3,3) * var(:,3+ol,:) + a(4,3) * var(:,4+ol,:) + &
                     a(5,3) * var(:,5+ol,:) + a(6,3) * var(:,6+ol,:) + &
                     a(7,3) * var(:,7+ol,:) ) * idel
      rhs(:,4+ol,:) = rhs(:,4+ol,:) + &
                   ( a(1,4) * var(:,1+ol,:) + a(2,4) * var(:,2+ol,:) + &
                     a(3,4) * var(:,3+ol,:) + a(4,4) * var(:,4+ol,:) + &
                     a(5,4) * var(:,5+ol,:) + a(6,4) * var(:,6+ol,:) + &
                     a(7,4) * var(:,7+ol,:) + a(8,4) * var(:,8+ol,:) ) * idel
      rhs(:,5+ol,:) = rhs(:,5+ol,:) + &
                   ( a(1,5) * var(:,1+ol,:) + a(2,5) * var(:,2+ol,:) + &
                     a(3,5) * var(:,3+ol,:) + a(4,5) * var(:,4+ol,:) + &
                     a(5,5) * var(:,5+ol,:) + a(6,5) * var(:,6+ol,:) + &
                     a(7,5) * var(:,7+ol,:) + a(8,5) * var(:,8+ol,:) + &
                     a(9,5) * var(:,9+ol,:) ) * idel
      rhs(:,6+ol,:) = rhs(:,6+ol,:) + &
                   ( a(2,6) * var(:,2+ol,:) + a(3,6) * var(:,3+ol,:) + &
                     a(4,6) * var(:,4+ol,:) + a(5,6) * var(:,5+ol,:) + &
                     a(6,6) * var(:,6+ol,:) + a(7,6) * var(:,7+ol,:) + &
                     a(8,6) * var(:,8+ol,:) + a(9,6) * var(:,9+ol,:) + &
                     a(10,6) * var(:,10+ol,:) ) * idel
!$OMP END PARALLEL WORKSHARE
  
      jl = 7 + ol
    end if
    if ( bb(4) == 0 ) then
      jr = nj - gsize(2)
    else
      or = nj - offset(4)

!$OMP PARALLEL WORKSHARE
      rhs(:,or-5,:) = rhs(:,or-5,:) + &
                   ( a(2,6) * var(:,or-1,:) + a(3,6) * var(:,or-2,:) + &
                     a(4,6) * var(:,or-3,:) + a(5,6) * var(:,or-4,:) + &
                     a(6,6) * var(:,or-5,:) + a(7,6) * var(:,or-6,:) + &
                     a(8,6) * var(:,or-7,:) + a(9,6) * var(:,or-8,:) + &
                     a(10,6) * var(:,or-9,:) ) * idel
      rhs(:,or-4,:) = rhs(:,or-4,:) + &
                   ( a(1,5) * var(:,or,:) + a(2,5) * var(:,or-1,:) + &
                     a(3,5) * var(:,or-2,:) + a(4,5) * var(:,or-3,:) + &
                     a(5,5) * var(:,or-4,:) + a(6,5) * var(:,or-5,:) + &
                     a(7,5) * var(:,or-6,:) + a(8,5) * var(:,or-7,:) + &
                     a(9,5) * var(:,or-8,:) ) * idel
      rhs(:,or-3,:) = rhs(:,or-3,:) + &
                   ( a(1,4) * var(:,or,:) + a(2,4) * var(:,or-1,:) + &
                     a(3,4) * var(:,or-2,:) + a(4,4) * var(:,or-3,:) + &
                     a(5,4) * var(:,or-4,:) + a(6,4) * var(:,or-5,:) + &
                     a(7,4) * var(:,or-6,:) + a(8,4) * var(:,or-7,:) ) * idel
      rhs(:,or-2,:) = rhs(:,or-2,:) + &
                   ( a(1,3) * var(:,or,:) + a(2,3) * var(:,or-1,:) + &
                     a(3,3) * var(:,or-2,:) + a(4,3) * var(:,or-3,:) + &
                     a(5,3) * var(:,or-4,:) + a(6,3) * var(:,or-5,:) + &
                     a(7,3) * var(:,or-6,:) ) * idel
      rhs(:,or-1,:) = rhs(:,or-1,:) + &
                   ( a(1,2) * var(:,or,:) + a(2,2) * var(:,or-1,:) + &
                     a(3,2) * var(:,or-2,:) + a(4,2) * var(:,or-3,:) + &
                     a(5,2) * var(:,or-4,:) + a(6,2) * var(:,or-5,:) ) * idel
      rhs(:,or,:)   = rhs(:,or,:) + &
                   ( a(1,1) * var(:,or,:) + a(2,1) * var(:,or-1,:) + &
                     a(3,1) * var(:,or-2,:) + a(4,1) * var(:,or-3,:) + &
                     a(5,1) * var(:,or-4,:) ) * idel
!$OMP END PARALLEL WORKSHARE
  
      jr = or - 6
    end if

!$OMP PARALLEL WORKSHARE
    rhs(:,jl:jr,:) = rhs(:,jl:jr,:) + &
                     ( -70.0_wp * var(:,jl:jr,:) + &
                        56.0_wp * ( var(:,jl-1:jr-1,:) + &
                                    var(:,jl+1:jr+1,:) ) - &
                        28.0_wp * ( var(:,jl-2:jr-2,:) + &
                                    var(:,jl+2:jr+2,:) ) + &
                         8.0_wp * ( var(:,jl-3:jr-3,:) + &
                                    var(:,jl+3:jr+3,:) ) - &
                                  ( var(:,jl-4:jr-4,:) + &
                                    var(:,jl+4:jr+4,:) ) ) * idel
!$OMP END PARALLEL WORKSHARE
  end if

  if ( zero_derivs_z == 0 ) then
    call set_coeff ( a )
  
    idel = epsilon / ( 256 * delta(3) )
  
    if ( bb(5) == 0 ) then
      kl = 1 + gsize(3)
    else
      ol = offset(5)

!$OMP PARALLEL WORKSHARE
      rhs(:,:,1+ol) = rhs(:,:,1+ol) + &
                   ( a(1,1) * var(:,:,1+ol) + a(2,1) * var(:,:,2+ol) + &
                     a(3,1) * var(:,:,3+ol) + a(4,1) * var(:,:,4+ol) + &
                     a(5,1) * var(:,:,5+ol) ) * idel
      rhs(:,:,2+ol) = rhs(:,:,2+ol) + &
                   ( a(1,2) * var(:,:,1+ol) + a(2,2) * var(:,:,2+ol) + &
                     a(3,2) * var(:,:,3+ol) + a(4,2) * var(:,:,4+ol) + &
                     a(5,2) * var(:,:,5+ol) + a(6,2) * var(:,:,6+ol) ) * idel
      rhs(:,:,3+ol) = rhs(:,:,3+ol) + &
                   ( a(1,3) * var(:,:,1+ol) + a(2,3) * var(:,:,2+ol) + &
                     a(3,3) * var(:,:,3+ol) + a(4,3) * var(:,:,4+ol) + &
                     a(5,3) * var(:,:,5+ol) + a(6,3) * var(:,:,6+ol) + &
                     a(7,3) * var(:,:,7+ol) ) * idel
      rhs(:,:,4+ol) = rhs(:,:,4+ol) + &
                   ( a(1,4) * var(:,:,1+ol) + a(2,4) * var(:,:,2+ol) + &
                     a(3,4) * var(:,:,3+ol) + a(4,4) * var(:,:,4+ol) + &
                     a(5,4) * var(:,:,5+ol) + a(6,4) * var(:,:,6+ol) + &
                     a(7,4) * var(:,:,7+ol) + a(8,4) * var(:,:,8+ol) ) * idel
      rhs(:,:,5+ol) = rhs(:,:,5+ol) + &
                   ( a(1,5) * var(:,:,1+ol) + a(2,5) * var(:,:,2+ol) + &
                     a(3,5) * var(:,:,3+ol) + a(4,5) * var(:,:,4+ol) + &
                     a(5,5) * var(:,:,5+ol) + a(6,5) * var(:,:,6+ol) + &
                     a(7,5) * var(:,:,7+ol) + a(8,5) * var(:,:,8+ol) + &
                     a(9,5) * var(:,:,9+ol) ) * idel
      rhs(:,:,6+ol) = rhs(:,:,6+ol) + &
                   ( a(2,6) * var(:,:,2+ol) + a(3,6) * var(:,:,3+ol) + &
                     a(4,6) * var(:,:,4+ol) + a(5,6) * var(:,:,5+ol) + &
                     a(6,6) * var(:,:,6+ol) + a(7,6) * var(:,:,7+ol) + &
                     a(8,6) * var(:,:,8+ol) + a(9,6) * var(:,:,9+ol) + &
                     a(10,6) * var(:,:,10+ol) ) * idel
!$OMP END PARALLEL WORKSHARE
  
      kl = 7 + ol
    end if
    if ( bb(6) == 0 ) then
      kr = nk - gsize(3)
    else
      or = nk - offset(6)

!$OMP PARALLEL WORKSHARE
      rhs(:,:,or-5) = rhs(:,:,or-5) + &
                   ( a(2,6) * var(:,:,or-1) + a(3,6) * var(:,:,or-2) + &
                     a(4,6) * var(:,:,or-3) + a(5,6) * var(:,:,or-4) + &
                     a(6,6) * var(:,:,or-5) + a(7,6) * var(:,:,or-6) + &
                     a(8,6) * var(:,:,or-7) + a(9,6) * var(:,:,or-8) + &
                     a(10,6) * var(:,:,or-9) ) * idel
      rhs(:,:,or-4) = rhs(:,:,or-4) + &
                   ( a(1,5) * var(:,:,or) + a(2,5) * var(:,:,or-1) + &
                     a(3,5) * var(:,:,or-2) + a(4,5) * var(:,:,or-3) + &
                     a(5,5) * var(:,:,or-4) + a(6,5) * var(:,:,or-5) + &
                     a(7,5) * var(:,:,or-6) + a(8,5) * var(:,:,or-7) + &
                     a(9,5) * var(:,:,or-8) ) * idel
      rhs(:,:,or-3) = rhs(:,:,or-3) + &
                   ( a(1,4) * var(:,:,or) + a(2,4) * var(:,:,or-1) + &
                     a(3,4) * var(:,:,or-2) + a(4,4) * var(:,:,or-3) + &
                     a(5,4) * var(:,:,or-4) + a(6,4) * var(:,:,or-5) + &
                     a(7,4) * var(:,:,or-6) + a(8,4) * var(:,:,or-7) ) * idel
      rhs(:,:,or-2) = rhs(:,:,or-2) + &
                   ( a(1,3) * var(:,:,or) + a(2,3) * var(:,:,or-1) + &
                     a(3,3) * var(:,:,or-2) + a(4,3) * var(:,:,or-3) + &
                     a(5,3) * var(:,:,or-4) + a(6,3) * var(:,:,or-5) + &
                     a(7,3) * var(:,:,or-6) ) * idel
      rhs(:,:,or-1) = rhs(:,:,or-1) + &
                   ( a(1,2) * var(:,:,or) + a(2,2) * var(:,:,or-1) + &
                     a(3,2) * var(:,:,or-2) + a(4,2) * var(:,:,or-3) + &
                     a(5,2) * var(:,:,or-4) + a(6,2) * var(:,:,or-5) ) * idel
      rhs(:,:,or)   = rhs(:,:,or) + &
                   ( a(1,1) * var(:,:,or) + a(2,1) * var(:,:,or-1) + &
                     a(3,1) * var(:,:,or-2) + a(4,1) * var(:,:,or-3) + &
                     a(5,1) * var(:,:,or-4) ) * idel
!$OMP END PARALLEL WORKSHARE

      kr = or - 6
    end if

!$OMP PARALLEL WORKSHARE
    rhs(:,:,kl:kr) = rhs(:,:,kl:kr) + &
                     ( -70.0_wp * var(:,:,kl:kr) + &
                        56.0_wp * ( var(:,:,kl-1:kr-1) + &
                                    var(:,:,kl+1:kr+1) ) - &
                        28.0_wp * ( var(:,:,kl-2:kr-2) + &
                                    var(:,:,kl+2:kr+2) ) + &
                         8.0_wp * ( var(:,:,kl-3:kr-3) + &
                                    var(:,:,kl+3:kr+3) ) - &
                                  ( var(:,:,kl-4:kr-4) + &
                                    var(:,:,kl+4:kr+4) ) ) * idel
!$OMP END PARALLEL WORKSHARE
  end if

contains

  subroutine set_coeff ( a )

    implicit none

    CCTK_REAL, dimension(10,6), intent(OUT) :: a
    CCTK_REAL :: zero = 0.0
    integer, parameter :: wp = kind(zero)

    a(1,1) = -3.1650670378782328375705179866656897941241116565316_wp
    a(2,1) = 12.660268151512931350282071946662759176496446626126_wp
    a(3,1) = -18.990402227269397025423107919994138764744669939190_wp
    a(4,1) = 12.660268151512931350282071946662759176496446626126_wp
    a(5,1) = -3.1650670378782328375705179866656897941241116565316_wp
    a(6,1) = 0.0_wp
    a(7,1) = 0.0_wp
    a(8,1) = 0.0_wp
    a(9,1) = 0.0_wp
    a(10,1) = 0.0_wp
    a(1,2) = 2.8768833763422958461666527928077915591442603845834_wp
    a(2,2) = -12.226754349454757346208274369433114126363106634479_wp
    a(3,2) = 20.138183634396070923166569549654540914009822692084_wp
    a(4,2) = -15.822858569882627153916590360442853575293432115209_wp
    a(5,2) = 5.7537667526845916923333055856155831182885207691667_wp
    a(6,2) = -0.71922084408557396154166319820194788978606509614584_wp
    a(7,2) = 0.0_wp
    a(8,2) = 0.0_wp
    a(9,2) = 0.0_wp
    a(10,2) = 0.0_wp
    a(1,3) = -9.5610475839173736628550350424197713021025451862781_wp
    a(2,3) = 44.618222058281077093323496864625599409811877535965_wp
    a(3,3) = -84.455920324603467355219476208041313168572482478790_wp
    a(4,3) = 82.862412393950571744743637034304684618222058281077_wp
    a(5,3) = -44.618222058281077093323496864625599409811877535965_wp
    a(6,3) = 12.748063445223164883806713389893028402803393581704_wp
    a(7,3) = -1.5935079306528956104758391737366285503504241977130_wp
    a(8,3) = 0.0_wp
    a(9,3) = 0.0_wp
    a(10,3) = 0.0_wp
    a(1,4) = 3.2244821795111028176898675125956335137152453816010_wp
    a(2,4) = -17.734651987311065497294271319275984325433849598806_wp
    a(3,4) = 41.918268333644336629968277663743235678298189960814_wp
    a(4,4) = -55.622317596566523605150214592274678111587982832618_wp
    a(5,4) = 45.142750513155439447658145176338869192013435342415_wp
    a(6,4) = -22.571375256577719723829072588169434596006717671207_wp
    a(7,4) = 6.4489643590222056353797350251912670274304907632021_wp
    a(8,4) = -0.80612054487777570442246687814890837842881134540026_wp
    a(9,4) = 0.0_wp
    a(10,4) = 0.0_wp
    a(1,5) = -1.0968642884346832550463374381109559476958232829758_wp
    a(2,5) = 8.7749143074774660403706995048876475815665862638060_wp
    a(3,5) = -30.712200076171131141297448267106766535483051923321_wp
    a(4,5) = 61.424400152342262282594896534213533070966103846642_wp
    a(5,5) = -76.780500190427827853243620667766916338707629808303_wp
    a(6,5) = 61.424400152342262282594896534213533070966103846642_wp
    a(7,5) = -30.712200076171131141297448267106766535483051923321_wp
    a(8,5) = 8.7749143074774660403706995048876475815665862638060_wp
    a(9,5) = -1.0968642884346832550463374381109559476958232829758_wp
    a(10,5) = 0.0_wp
    a(1,6) = 0.0_wp
    a(2,6) = -0.98627885208100271683294901942878016483641925983425_wp
    a(3,6) = 7.8902308166480217346635921554302413186913540786740_wp
    a(4,6) = -27.615807858268076071322572544005844615419739275359_wp
    a(5,6) = 55.231615716536152142645145088011689230839478550718_wp
    a(6,6) = -69.039519645670190178306431360014611538549348188398_wp
    a(7,6) = 55.231615716536152142645145088011689230839478550718_wp
    a(8,6) = -27.615807858268076071322572544005844615419739275359_wp
    a(9,6) = 7.8902308166480217346635921554302413186913540786740_wp
    a(10,6) = -0.98627885208100271683294901942878016483641925983425_wp

  end subroutine set_coeff

end subroutine dissipation_6_3_alt

subroutine dissipation_6_3_delta (var, ni, nj, nk, bb, gsize, offset, &
                                  dx, dy, dz, epsilon, rhs)

  implicit none

  DECLARE_CCTK_PARAMETERS

  integer ::  ni, nj, nk
  CCTK_REAL, dimension(ni,nj,nk), intent(in) :: var
  CCTK_REAL, dimension(ni,nj,nk), intent(inout) :: rhs
  CCTK_INT, dimension(6), intent(in) :: bb
  CCTK_INT, dimension(3), intent(in) :: gsize
  CCTK_INT, dimension(6), intent(in) :: offset
  CCTK_REAL, dimension(ni,nj,nk), intent(in) :: dx, dy, dz
  CCTK_REAL, intent(in) :: epsilon 

  CCTK_REAL :: zero = 0.0
  integer, parameter :: wp = kind(zero)
  CCTK_REAL, dimension(10,6) :: a
  CCTK_REAL :: idel

  CCTK_INT :: il, ir, jl, jr, kl, kr, ol, or

  call set_coeff ( a )

  idel = epsilon / 256

  if ( bb(1) == 0 ) then
    il = 1 + gsize(1)
  else
    ol = offset(1)

!$OMP PARALLEL WORKSHARE
    rhs(1+ol,:,:) = rhs(1+ol,:,:) + &
                 ( a(1,1) * var(1+ol,:,:) + a(2,1) * var(2+ol,:,:) + &
                   a(3,1) * var(3+ol,:,:) + a(4,1) * var(4+ol,:,:) + &
                   a(5,1) * var(5+ol,:,:) ) * idel / dx(1+ol,:,:)
    rhs(2+ol,:,:) = rhs(2+ol,:,:) + &
                 ( a(1,2) * var(1+ol,:,:) + a(2,2) * var(2+ol,:,:) + &
                   a(3,2) * var(3+ol,:,:) + a(4,2) * var(4+ol,:,:) + &
                   a(5,2) * var(5+ol,:,:) + a(6,2) * var(6+ol,:,:) ) * &
                 idel / dx(2+ol,:,:)
    rhs(3+ol,:,:) = rhs(3+ol,:,:) + &
                 ( a(1,3) * var(1+ol,:,:) + a(2,3) * var(2+ol,:,:) + &
                   a(3,3) * var(3+ol,:,:) + a(4,3) * var(4+ol,:,:) + &
                   a(5,3) * var(5+ol,:,:) + a(6,3) * var(6+ol,:,:) + &
                   a(7,3) * var(7+ol,:,:) ) * idel / dx(3+ol,:,:)
    rhs(4+ol,:,:) = rhs(4+ol,:,:) + &
                 ( a(1,4) * var(1+ol,:,:) + a(2,4) * var(2+ol,:,:) + &
                   a(3,4) * var(3+ol,:,:) + a(4,4) * var(4+ol,:,:) + &
                   a(5,4) * var(5+ol,:,:) + a(6,4) * var(6+ol,:,:) + &
                   a(7,4) * var(7+ol,:,:) + a(8,4) * var(8+ol,:,:) ) * &
                 idel / dx(4+ol,:,:)
    rhs(5+ol,:,:) = rhs(5+ol,:,:) + &
                 ( a(1,5) * var(1+ol,:,:) + a(2,5) * var(2+ol,:,:) + &
                   a(3,5) * var(3+ol,:,:) + a(4,5) * var(4+ol,:,:) + &
                   a(5,5) * var(5+ol,:,:) + a(6,5) * var(6+ol,:,:) + &
                   a(7,5) * var(7+ol,:,:) + a(8,5) * var(8+ol,:,:) + &
                   a(9,5) * var(9+ol,:,:) ) * idel / dx(5+ol,:,:)
    rhs(6+ol,:,:) = rhs(6+ol,:,:) + &
                 ( a(2,6) * var(2+ol,:,:) + a(3,6) * var(3+ol,:,:) + &
                   a(4,6) * var(4+ol,:,:) + a(5,6) * var(5+ol,:,:) + &
                   a(6,6) * var(6+ol,:,:) + a(7,6) * var(7+ol,:,:) + &
                   a(8,6) * var(8+ol,:,:) + a(9,6) * var(9+ol,:,:) + &
                   a(10,6) * var(10+ol,:,:) ) * idel / dx(6+ol,:,:)
!$OMP END PARALLEL WORKSHARE

    il = 7 + ol
  end if
  if ( bb(2) == 0 ) then
    ir = ni - gsize(1)
  else
    or = ni - offset(2)

!$OMP PARALLEL WORKSHARE
    rhs(or-5,:,:) = rhs(or-5,:,:) + &
                 ( a(2,6) * var(or-1,:,:) + a(3,6) * var(or-2,:,:) + &
                   a(4,6) * var(or-3,:,:) + a(5,6) * var(or-4,:,:) + &
                   a(6,6) * var(or-5,:,:) + a(7,6) * var(or-6,:,:) + &
                   a(8,6) * var(or-7,:,:) + a(9,6) * var(or-8,:,:) + &
                   a(10,6) * var(or-9,:,:) ) * idel / dx(or-5,:,:)
    rhs(or-4,:,:) = rhs(or-4,:,:) + &
                 ( a(1,5) * var(or,:,:) + a(2,5) * var(or-1,:,:) + &
                   a(3,5) * var(or-2,:,:) + a(4,5) * var(or-3,:,:) + &
                   a(5,5) * var(or-4,:,:) + a(6,5) * var(or-5,:,:) + &
                   a(7,5) * var(or-6,:,:) + a(8,5) * var(or-7,:,:) + &
                   a(9,5) * var(or-8,:,:) ) * idel / dx(or-4,:,:)
    rhs(or-3,:,:) = rhs(or-3,:,:) + &
                 ( a(1,4) * var(or,:,:) + a(2,4) * var(or-1,:,:) + &
                   a(3,4) * var(or-2,:,:) + a(4,4) * var(or-3,:,:) + &
                   a(5,4) * var(or-4,:,:) + a(6,4) * var(or-5,:,:) + &
                   a(7,4) * var(or-6,:,:) + a(8,4) * var(or-7,:,:) ) * &
                 idel / dx(or-3,:,:)
    rhs(or-2,:,:) = rhs(or-2,:,:) + &
                 ( a(1,3) * var(or,:,:) + a(2,3) * var(or-1,:,:) + &
                   a(3,3) * var(or-2,:,:) + a(4,3) * var(or-3,:,:) + &
                   a(5,3) * var(or-4,:,:) + a(6,3) * var(or-5,:,:) + &
                   a(7,3) * var(or-6,:,:) ) * idel / dx(or-2,:,:)
    rhs(or-1,:,:) = rhs(or-1,:,:) + &
                 ( a(1,2) * var(or,:,:) + a(2,2) * var(or-1,:,:) + &
                   a(3,2) * var(or-2,:,:) + a(4,2) * var(or-3,:,:) + &
                   a(5,2) * var(or-4,:,:) + a(6,2) * var(or-5,:,:) ) * &
                 idel / dx(or-1,:,:)
    rhs(or,:,:)   = rhs(or,:,:) + &
                 ( a(1,1) * var(or,:,:) + a(2,1) * var(or-1,:,:) + &
                   a(3,1) * var(or-2,:,:) + a(4,1) * var(or-3,:,:) + &
                   a(5,1) * var(or-4,:,:) ) * idel / dx(or,:,:)
!$OMP END PARALLEL WORKSHARE

    ir = or - 6
  end if

!$OMP PARALLEL WORKSHARE
  rhs(il:ir,:,:) = rhs(il:ir,:,:) + &
                   ( -70.0_wp * var(il:ir,:,:) + &
                      56.0_wp * ( var(il-1:ir-1,:,:) + &
                                  var(il+1:ir+1,:,:) ) - &
                      28.0_wp * ( var(il-2:ir-2,:,:) + &
                                  var(il+2:ir+2,:,:) ) + &
                       8.0_wp * ( var(il-3:ir-3,:,:) + &
                                  var(il+3:ir+3,:,:) ) - &
                                ( var(il-4:ir-4,:,:) + &
                                  var(il+4:ir+4,:,:) ) ) * idel / dx(il:ir,:,:)
!$OMP END PARALLEL WORKSHARE

  if ( zero_derivs_y == 0 ) then
    call set_coeff ( a )
  
    idel = epsilon /  256
  
    if ( bb(3) == 0 ) then
      jl = 1 + gsize(2)
    else
      ol = offset(3)

!$OMP PARALLEL WORKSHARE
      rhs(:,1+ol,:) = rhs(:,1+ol,:) + &
                   ( a(1,1) * var(:,1+ol,:) + a(2,1) * var(:,2+ol,:) + &
                     a(3,1) * var(:,3+ol,:) + a(4,1) * var(:,4+ol,:) + &
                     a(5,1) * var(:,5+ol,:) ) * idel / dy(:,1+ol,:)
      rhs(:,2+ol,:) = rhs(:,2+ol,:) + &
                   ( a(1,2) * var(:,1+ol,:) + a(2,2) * var(:,2+ol,:) + &
                     a(3,2) * var(:,3+ol,:) + a(4,2) * var(:,4+ol,:) + &
                     a(5,2) * var(:,5+ol,:) + a(6,2) * var(:,6+ol,:) ) * &
                   idel / dy(:,2+ol,:)
      rhs(:,3+ol,:) = rhs(:,3+ol,:) + &
                   ( a(1,3) * var(:,1+ol,:) + a(2,3) * var(:,2+ol,:) + &
                     a(3,3) * var(:,3+ol,:) + a(4,3) * var(:,4+ol,:) + &
                     a(5,3) * var(:,5+ol,:) + a(6,3) * var(:,6+ol,:) + &
                     a(7,3) * var(:,7+ol,:) ) * idel / dy(:,3+ol,:)
      rhs(:,4+ol,:) = rhs(:,4+ol,:) + &
                   ( a(1,4) * var(:,1+ol,:) + a(2,4) * var(:,2+ol,:) + &
                     a(3,4) * var(:,3+ol,:) + a(4,4) * var(:,4+ol,:) + &
                     a(5,4) * var(:,5+ol,:) + a(6,4) * var(:,6+ol,:) + &
                     a(7,4) * var(:,7+ol,:) + a(8,4) * var(:,8+ol,:) ) * &
                   idel / dy(:,4+ol,:)
      rhs(:,5+ol,:) = rhs(:,5+ol,:) + &
                   ( a(1,5) * var(:,1+ol,:) + a(2,5) * var(:,2+ol,:) + &
                     a(3,5) * var(:,3+ol,:) + a(4,5) * var(:,4+ol,:) + &
                     a(5,5) * var(:,5+ol,:) + a(6,5) * var(:,6+ol,:) + &
                     a(7,5) * var(:,7+ol,:) + a(8,5) * var(:,8+ol,:) + &
                     a(9,5) * var(:,9+ol,:) ) * idel / dy(:,5+ol,:)
      rhs(:,6+ol,:) = rhs(:,6+ol,:) + &
                   ( a(2,6) * var(:,2+ol,:) + a(3,6) * var(:,3+ol,:) + &
                     a(4,6) * var(:,4+ol,:) + a(5,6) * var(:,5+ol,:) + &
                     a(6,6) * var(:,6+ol,:) + a(7,6) * var(:,7+ol,:) + &
                     a(8,6) * var(:,8+ol,:) + a(9,6) * var(:,9+ol,:) + &
                     a(10,6) * var(:,10+ol,:) ) * idel / dy(:,6+ol,:)
!$OMP END PARALLEL WORKSHARE
  
      jl = 7 + ol
    end if
    if ( bb(4) == 0 ) then
      jr = nj - gsize(2)
    else
      or = nj - offset(4)

!$OMP PARALLEL WORKSHARE
      rhs(:,or-5,:) = rhs(:,or-5,:) + &
                   ( a(2,6) * var(:,or-1,:) + a(3,6) * var(:,or-2,:) + &
                     a(4,6) * var(:,or-3,:) + a(5,6) * var(:,or-4,:) + &
                     a(6,6) * var(:,or-5,:) + a(7,6) * var(:,or-6,:) + &
                     a(8,6) * var(:,or-7,:) + a(9,6) * var(:,or-8,:) + &
                     a(10,6) * var(:,or-9,:) ) * idel / dy(:,or-5,:)
      rhs(:,or-4,:) = rhs(:,or-4,:) + &
                   ( a(1,5) * var(:,or,:) + a(2,5) * var(:,or-1,:) + &
                     a(3,5) * var(:,or-2,:) + a(4,5) * var(:,or-3,:) + &
                     a(5,5) * var(:,or-4,:) + a(6,5) * var(:,or-5,:) + &
                     a(7,5) * var(:,or-6,:) + a(8,5) * var(:,or-7,:) + &
                     a(9,5) * var(:,or-8,:) ) * idel / dy(:,or-4,:)
      rhs(:,or-3,:) = rhs(:,or-3,:) + &
                   ( a(1,4) * var(:,or,:) + a(2,4) * var(:,or-1,:) + &
                     a(3,4) * var(:,or-2,:) + a(4,4) * var(:,or-3,:) + &
                     a(5,4) * var(:,or-4,:) + a(6,4) * var(:,or-5,:) + &
                     a(7,4) * var(:,or-6,:) + a(8,4) * var(:,or-7,:) ) * &
                   idel / dy(:,or-3,:)
      rhs(:,or-2,:) = rhs(:,or-2,:) + &
                   ( a(1,3) * var(:,or,:) + a(2,3) * var(:,or-1,:) + &
                     a(3,3) * var(:,or-2,:) + a(4,3) * var(:,or-3,:) + &
                     a(5,3) * var(:,or-4,:) + a(6,3) * var(:,or-5,:) + &
                     a(7,3) * var(:,or-6,:) ) * idel / dy(:,or-2,:)
      rhs(:,or-1,:) = rhs(:,or-1,:) + &
                   ( a(1,2) * var(:,or,:) + a(2,2) * var(:,or-1,:) + &
                     a(3,2) * var(:,or-2,:) + a(4,2) * var(:,or-3,:) + &
                     a(5,2) * var(:,or-4,:) + a(6,2) * var(:,or-5,:) ) * &
                   idel / dy(:,or-1,:)
      rhs(:,or,:)   = rhs(:,or,:) + &
                   ( a(1,1) * var(:,or,:) + a(2,1) * var(:,or-1,:) + &
                     a(3,1) * var(:,or-2,:) + a(4,1) * var(:,or-3,:) + &
                     a(5,1) * var(:,or-4,:) ) * idel / dy(:,or,:)
!$OMP END PARALLEL WORKSHARE
  
      jr = or - 6
    end if

!$OMP PARALLEL WORKSHARE
    rhs(:,jl:jr,:) = rhs(:,jl:jr,:) + &
                     ( -70.0_wp * var(:,jl:jr,:) + &
                        56.0_wp * ( var(:,jl-1:jr-1,:) + &
                                    var(:,jl+1:jr+1,:) ) - &
                        28.0_wp * ( var(:,jl-2:jr-2,:) + &
                                    var(:,jl+2:jr+2,:) ) + &
                         8.0_wp * ( var(:,jl-3:jr-3,:) + &
                                    var(:,jl+3:jr+3,:) ) - &
                                  ( var(:,jl-4:jr-4,:) + &
                                    var(:,jl+4:jr+4,:) ) ) * &
                     idel / dy(:,jl:jr,:)
!$OMP END PARALLEL WORKSHARE
  end if

  if ( zero_derivs_z == 0 ) then
    call set_coeff ( a )
  
    idel = epsilon / 256
  
    if ( bb(5) == 0 ) then
      kl = 1 + gsize(3)
    else
      ol = offset(5)

!$OMP PARALLEL WORKSHARE
      rhs(:,:,1+ol) = rhs(:,:,1+ol) + &
                   ( a(1,1) * var(:,:,1+ol) + a(2,1) * var(:,:,2+ol) + &
                     a(3,1) * var(:,:,3+ol) + a(4,1) * var(:,:,4+ol) + &
                     a(5,1) * var(:,:,5+ol) ) * idel / dz(:,:,1+ol)
      rhs(:,:,2+ol) = rhs(:,:,2+ol) + &
                   ( a(1,2) * var(:,:,1+ol) + a(2,2) * var(:,:,2+ol) + &
                     a(3,2) * var(:,:,3+ol) + a(4,2) * var(:,:,4+ol) + &
                     a(5,2) * var(:,:,5+ol) + a(6,2) * var(:,:,6+ol) ) * &
                   idel / dz(:,:,2+ol)
      rhs(:,:,3+ol) = rhs(:,:,3+ol) + &
                   ( a(1,3) * var(:,:,1+ol) + a(2,3) * var(:,:,2+ol) + &
                     a(3,3) * var(:,:,3+ol) + a(4,3) * var(:,:,4+ol) + &
                     a(5,3) * var(:,:,5+ol) + a(6,3) * var(:,:,6+ol) + &
                     a(7,3) * var(:,:,7+ol) ) * idel / dz(:,:,3+ol)
      rhs(:,:,4+ol) = rhs(:,:,4+ol) + &
                   ( a(1,4) * var(:,:,1+ol) + a(2,4) * var(:,:,2+ol) + &
                     a(3,4) * var(:,:,3+ol) + a(4,4) * var(:,:,4+ol) + &
                     a(5,4) * var(:,:,5+ol) + a(6,4) * var(:,:,6+ol) + &
                     a(7,4) * var(:,:,7+ol) + a(8,4) * var(:,:,8+ol) ) * &
                   idel / dz(:,:,4+ol)
      rhs(:,:,5+ol) = rhs(:,:,5+ol) + &
                   ( a(1,5) * var(:,:,1+ol) + a(2,5) * var(:,:,2+ol) + &
                     a(3,5) * var(:,:,3+ol) + a(4,5) * var(:,:,4+ol) + &
                     a(5,5) * var(:,:,5+ol) + a(6,5) * var(:,:,6+ol) + &
                     a(7,5) * var(:,:,7+ol) + a(8,5) * var(:,:,8+ol) + &
                     a(9,5) * var(:,:,9+ol) ) * idel / dz(:,:,5+ol)
      rhs(:,:,6+ol) = rhs(:,:,6+ol) + &
                   ( a(2,6) * var(:,:,2+ol) + a(3,6) * var(:,:,3+ol) + &
                     a(4,6) * var(:,:,4+ol) + a(5,6) * var(:,:,5+ol) + &
                     a(6,6) * var(:,:,6+ol) + a(7,6) * var(:,:,7+ol) + &
                     a(8,6) * var(:,:,8+ol) + a(9,6) * var(:,:,9+ol) + &
                     a(10,6) * var(:,:,10+ol) ) * idel / dz(:,:,6+ol)
!$OMP END PARALLEL WORKSHARE
  
      kl = 7 + ol
    end if
    if ( bb(6) == 0 ) then
      kr = nk - gsize(3)
    else
      or = nk - offset(6)

!$OMP PARALLEL WORKSHARE
      rhs(:,:,or-5) = rhs(:,:,or-5) + &
                   ( a(2,6) * var(:,:,or-1) + a(3,6) * var(:,:,or-2) + &
                     a(4,6) * var(:,:,or-3) + a(5,6) * var(:,:,or-4) + &
                     a(6,6) * var(:,:,or-5) + a(7,6) * var(:,:,or-6) + &
                     a(8,6) * var(:,:,or-7) + a(9,6) * var(:,:,or-8) + &
                     a(10,6) * var(:,:,or-9) ) * idel / dz(:,:,or-5)
      rhs(:,:,or-4) = rhs(:,:,or-4) + &
                   ( a(1,5) * var(:,:,or) + a(2,5) * var(:,:,or-1) + &
                     a(3,5) * var(:,:,or-2) + a(4,5) * var(:,:,or-3) + &
                     a(5,5) * var(:,:,or-4) + a(6,5) * var(:,:,or-5) + &
                     a(7,5) * var(:,:,or-6) + a(8,5) * var(:,:,or-7) + &
                     a(9,5) * var(:,:,or-8) ) * idel / dz(:,:,or-4)
      rhs(:,:,or-3) = rhs(:,:,or-3) + &
                   ( a(1,4) * var(:,:,or) + a(2,4) * var(:,:,or-1) + &
                     a(3,4) * var(:,:,or-2) + a(4,4) * var(:,:,or-3) + &
                     a(5,4) * var(:,:,or-4) + a(6,4) * var(:,:,or-5) + &
                     a(7,4) * var(:,:,or-6) + a(8,4) * var(:,:,or-7) ) * &
                   idel / dz(:,:,or-3)
      rhs(:,:,or-2) = rhs(:,:,or-2) + &
                   ( a(1,3) * var(:,:,or) + a(2,3) * var(:,:,or-1) + &
                     a(3,3) * var(:,:,or-2) + a(4,3) * var(:,:,or-3) + &
                     a(5,3) * var(:,:,or-4) + a(6,3) * var(:,:,or-5) + &
                     a(7,3) * var(:,:,or-6) ) * idel / dz(:,:,or-2)
      rhs(:,:,or-1) = rhs(:,:,or-1) + &
                   ( a(1,2) * var(:,:,or) + a(2,2) * var(:,:,or-1) + &
                     a(3,2) * var(:,:,or-2) + a(4,2) * var(:,:,or-3) + &
                     a(5,2) * var(:,:,or-4) + a(6,2) * var(:,:,or-5) ) * &
                   idel / dz(:,:,or-1)
      rhs(:,:,or)   = rhs(:,:,or) + &
                   ( a(1,1) * var(:,:,or) + a(2,1) * var(:,:,or-1) + &
                     a(3,1) * var(:,:,or-2) + a(4,1) * var(:,:,or-3) + &
                     a(5,1) * var(:,:,or-4) ) * idel / dz(:,:,or)
!$OMP END PARALLEL WORKSHARE

      kr = or - 6
    end if

!$OMP PARALLEL WORKSHARE
    rhs(:,:,kl:kr) = rhs(:,:,kl:kr) + &
                     ( -70.0_wp * var(:,:,kl:kr) + &
                        56.0_wp * ( var(:,:,kl-1:kr-1) + &
                                    var(:,:,kl+1:kr+1) ) - &
                        28.0_wp * ( var(:,:,kl-2:kr-2) + &
                                    var(:,:,kl+2:kr+2) ) + &
                         8.0_wp * ( var(:,:,kl-3:kr-3) + &
                                    var(:,:,kl+3:kr+3) ) - &
                                  ( var(:,:,kl-4:kr-4) + &
                                    var(:,:,kl+4:kr+4) ) ) * &
                     idel / dz(:,:,kl:kr)
!$OMP END PARALLEL WORKSHARE
  end if

contains

  subroutine set_coeff ( a )

    implicit none

    CCTK_REAL, dimension(10,6), intent(OUT) :: a
    CCTK_REAL :: zero = 0.0
    integer, parameter :: wp = kind(zero)

    a(1,1) = -3.1650670378782328375705179866656897941241116565316_wp
    a(2,1) = 12.660268151512931350282071946662759176496446626126_wp
    a(3,1) = -18.990402227269397025423107919994138764744669939190_wp
    a(4,1) = 12.660268151512931350282071946662759176496446626126_wp
    a(5,1) = -3.1650670378782328375705179866656897941241116565316_wp
    a(6,1) = 0.0_wp
    a(7,1) = 0.0_wp
    a(8,1) = 0.0_wp
    a(9,1) = 0.0_wp
    a(10,1) = 0.0_wp
    a(1,2) = 2.8768833763422958461666527928077915591442603845834_wp
    a(2,2) = -12.226754349454757346208274369433114126363106634479_wp
    a(3,2) = 20.138183634396070923166569549654540914009822692084_wp
    a(4,2) = -15.822858569882627153916590360442853575293432115209_wp
    a(5,2) = 5.7537667526845916923333055856155831182885207691667_wp
    a(6,2) = -0.71922084408557396154166319820194788978606509614584_wp
    a(7,2) = 0.0_wp
    a(8,2) = 0.0_wp
    a(9,2) = 0.0_wp
    a(10,2) = 0.0_wp
    a(1,3) = -9.5610475839173736628550350424197713021025451862781_wp
    a(2,3) = 44.618222058281077093323496864625599409811877535965_wp
    a(3,3) = -84.455920324603467355219476208041313168572482478790_wp
    a(4,3) = 82.862412393950571744743637034304684618222058281077_wp
    a(5,3) = -44.618222058281077093323496864625599409811877535965_wp
    a(6,3) = 12.748063445223164883806713389893028402803393581704_wp
    a(7,3) = -1.5935079306528956104758391737366285503504241977130_wp
    a(8,3) = 0.0_wp
    a(9,3) = 0.0_wp
    a(10,3) = 0.0_wp
    a(1,4) = 3.2244821795111028176898675125956335137152453816010_wp
    a(2,4) = -17.734651987311065497294271319275984325433849598806_wp
    a(3,4) = 41.918268333644336629968277663743235678298189960814_wp
    a(4,4) = -55.622317596566523605150214592274678111587982832618_wp
    a(5,4) = 45.142750513155439447658145176338869192013435342415_wp
    a(6,4) = -22.571375256577719723829072588169434596006717671207_wp
    a(7,4) = 6.4489643590222056353797350251912670274304907632021_wp
    a(8,4) = -0.80612054487777570442246687814890837842881134540026_wp
    a(9,4) = 0.0_wp
    a(10,4) = 0.0_wp
    a(1,5) = -1.0968642884346832550463374381109559476958232829758_wp
    a(2,5) = 8.7749143074774660403706995048876475815665862638060_wp
    a(3,5) = -30.712200076171131141297448267106766535483051923321_wp
    a(4,5) = 61.424400152342262282594896534213533070966103846642_wp
    a(5,5) = -76.780500190427827853243620667766916338707629808303_wp
    a(6,5) = 61.424400152342262282594896534213533070966103846642_wp
    a(7,5) = -30.712200076171131141297448267106766535483051923321_wp
    a(8,5) = 8.7749143074774660403706995048876475815665862638060_wp
    a(9,5) = -1.0968642884346832550463374381109559476958232829758_wp
    a(10,5) = 0.0_wp
    a(1,6) = 0.0_wp
    a(2,6) = -0.98627885208100271683294901942878016483641925983425_wp
    a(3,6) = 7.8902308166480217346635921554302413186913540786740_wp
    a(4,6) = -27.615807858268076071322572544005844615419739275359_wp
    a(5,6) = 55.231615716536152142645145088011689230839478550718_wp
    a(6,6) = -69.039519645670190178306431360014611538549348188398_wp
    a(7,6) = 55.231615716536152142645145088011689230839478550718_wp
    a(8,6) = -27.615807858268076071322572544005844615419739275359_wp
    a(9,6) = 7.8902308166480217346635921554302413186913540786740_wp
    a(10,6) = -0.98627885208100271683294901942878016483641925983425_wp

  end subroutine set_coeff

end subroutine dissipation_6_3_delta
