! $Header$

#include "cctk.h"
#include "cctk_Parameters.h"

subroutine dissipation_6_3 (var, ni, nj, nk, bb, gsize, offset, delta, epsilon, rhs)

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
  CCTK_REAL, dimension(9,6) :: a
  CCTK_REAL :: idel

  CCTK_INT :: il, ir, jl, jr, kl, kr, ol, or

  call set_coeff ( a )

  if ( scale_with_h > 0 ) then
    idel = epsilon / ( 64 * delta(1) )
  else
    idel = epsilon / 64
  end if

  if ( bb(1) == 0 ) then
    il = 1 + gsize(1)
  else
    ol = offset(1)
    rhs(1+ol,:,:) = rhs(1+ol,:,:) + &
                 ( a(1,1) * var(1+ol,:,:) + a(2,1) * var(2+ol,:,:) + &
                   a(3,1) * var(3+ol,:,:) + a(4,1) * var(4+ol,:,:) ) * idel
    rhs(2+ol,:,:) = rhs(2+ol,:,:) + &
                 ( a(1,2) * var(1+ol,:,:) + a(2,2) * var(2+ol,:,:) + &
                   a(3,2) * var(3+ol,:,:) + a(4,2) * var(4+ol,:,:) + &
                   a(5,2) * var(5+ol,:,:) ) * idel
    rhs(3+ol,:,:) = rhs(3+ol,:,:) + &
                 ( a(1,3) * var(1+ol,:,:) + a(2,3) * var(2+ol,:,:) + &
                   a(3,3) * var(3+ol,:,:) + a(4,3) * var(4+ol,:,:) + &
                   a(5,3) * var(5+ol,:,:) + a(6,3) * var(6+ol,:,:) ) * idel
    rhs(4+ol,:,:) = rhs(4+ol,:,:) + &
                 ( a(1,4) * var(1+ol,:,:) + a(2,4) * var(2+ol,:,:) + &
                   a(3,4) * var(3+ol,:,:) + a(4,4) * var(4+ol,:,:) + &
                   a(5,4) * var(5+ol,:,:) + a(6,4) * var(6+ol,:,:) + &
                   a(7,4) * var(7+ol,:,:) ) * idel
    rhs(5+ol,:,:) = rhs(5+ol,:,:) + &
                 ( a(2,5) * var(2+ol,:,:) + a(3,5) * var(3+ol,:,:) + &
                   a(4,5) * var(4+ol,:,:) + a(5,5) * var(5+ol,:,:) + &
                   a(6,5) * var(6+ol,:,:) + a(7,5) * var(7+ol,:,:) + &
                   a(8,5) * var(8+ol,:,:) ) * idel
    rhs(6+ol,:,:) = rhs(6+ol,:,:) + &
                 ( a(3,6) * var(3+ol,:,:) + a(4,6) * var(4+ol,:,:) + &
                   a(5,6) * var(5+ol,:,:) + a(6,6) * var(6+ol,:,:) + &
                   a(7,6) * var(7+ol,:,:) + a(8,6) * var(8+ol,:,:) + &
                   a(9,6) * var(9+ol,:,:) ) * idel

    il = 7 + ol
  end if
  if ( bb(2) == 0 ) then
    ir = ni - gsize(1)
  else
    or = ni - offset(2)
    rhs(or-5,:,:) = rhs(or-5,:,:) + &
                 ( a(3,6) * var(or-2,:,:) + a(4,6) * var(or-3,:,:) + &
                   a(5,6) * var(or-4,:,:) + a(6,6) * var(or-5,:,:) + &
                   a(7,6) * var(or-6,:,:) + a(8,6) * var(or-7,:,:) + &
                   a(9,6) * var(or-8,:,:) ) * idel
    rhs(or-4,:,:) = rhs(or-4,:,:) + &
                 ( a(2,5) * var(or-1,:,:) + a(3,5) * var(or-2,:,:) + &
                   a(4,5) * var(or-3,:,:) + a(5,5) * var(or-4,:,:) + &
                   a(6,5) * var(or-5,:,:) + a(7,5) * var(or-6,:,:) + &
                   a(8,5) * var(or-7,:,:) ) * idel
    rhs(or-3,:,:) = rhs(or-3,:,:) + &
                 ( a(1,4) * var(or,:,:) + a(2,4) * var(or-1,:,:) + &
                   a(3,4) * var(or-2,:,:) + a(4,4) * var(or-3,:,:) + &
                   a(5,4) * var(or-4,:,:) + a(6,4) * var(or-5,:,:) + &
                   a(7,4) * var(or-6,:,:) ) * idel
    rhs(or-2,:,:) = rhs(or-2,:,:) + &
                 ( a(1,3) * var(or,:,:) + a(2,3) * var(or-1,:,:) + &
                   a(3,3) * var(or-2,:,:) + a(4,3) * var(or-3,:,:) + &
                   a(5,3) * var(or-4,:,:) + a(6,3) * var(or-5,:,:) ) * idel
    rhs(or-1,:,:) = rhs(or-1,:,:) + &
                 ( a(1,2) * var(or,:,:) + a(2,2) * var(or-1,:,:) + &
                   a(3,2) * var(or-2,:,:) + a(4,2) * var(or-3,:,:) + &
                   a(5,2) * var(or-4,:,:) ) * idel
    rhs(or,:,:)   = rhs(or,:,:) + &
                 ( a(1,1) * var(or,:,:) + a(2,1) * var(or-1,:,:) + &
                   a(3,1) * var(or-2,:,:) + a(4,1) * var(or-3,:,:) ) * idel

    ir = or - 6
  end if
  rhs(il:ir,:,:) = rhs(il:ir,:,:) + &
                   ( -20.0_wp * var(il:ir,:,:) + &
                      15.0_wp * ( var(il-1:ir-1,:,:) + &
                                  var(il+1:ir+1,:,:) ) - &
                       6.0_wp * ( var(il-2:ir-2,:,:) + &
                                  var(il+2:ir+2,:,:) ) + &
                                ( var(il-3:ir-3,:,:) + &
                                  var(il+3:ir+3,:,:) ) ) * idel

  if ( zero_derivs_y == 0 ) then
    call set_coeff ( a )
  
    if ( scale_with_h > 0 ) then
      idel = epsilon / ( 64 * delta(2) )
    else
      idel = epsilon / 64
    end if
  
    if ( bb(3) == 0 ) then
      jl = 1 + gsize(2)
    else
      ol = offset(3)
      rhs(:,1+ol,:) = rhs(:,1+ol,:) + &
                   ( a(1,1) * var(:,1+ol,:) + a(2,1) * var(:,2+ol,:) + &
                     a(3,1) * var(:,3+ol,:) + a(4,1) * var(:,4+ol,:) ) * idel
      rhs(:,2+ol,:) = rhs(:,2+ol,:) + &
                   ( a(1,2) * var(:,1+ol,:) + a(2,2) * var(:,2+ol,:) + &
                     a(3,2) * var(:,3+ol,:) + a(4,2) * var(:,4+ol,:) + &
                     a(5,2) * var(:,5+ol,:) ) * idel
      rhs(:,3+ol,:) = rhs(:,3+ol,:) + &
                   ( a(1,3) * var(:,1+ol,:) + a(2,3) * var(:,2+ol,:) + &
                     a(3,3) * var(:,3+ol,:) + a(4,3) * var(:,4+ol,:) + &
                     a(5,3) * var(:,5+ol,:) + a(6,3) * var(:,6+ol,:) ) * idel
      rhs(:,4+ol,:) = rhs(:,4+ol,:) + &
                   ( a(1,4) * var(:,1+ol,:) + a(2,4) * var(:,2+ol,:) + &
                     a(3,4) * var(:,3+ol,:) + a(4,4) * var(:,4+ol,:) + &
                     a(5,4) * var(:,5+ol,:) + a(6,4) * var(:,6+ol,:) + &
                     a(7,4) * var(:,7+ol,:) ) * idel
      rhs(:,5+ol,:) = rhs(:,5+ol,:) + &
                   ( a(2,5) * var(:,2+ol,:) + a(3,5) * var(:,3+ol,:) + &
                     a(4,5) * var(:,4+ol,:) + a(5,5) * var(:,5+ol,:) + &
                     a(6,5) * var(:,6+ol,:) + a(7,5) * var(:,7+ol,:) + &
                     a(8,5) * var(:,8+ol,:) ) * idel
      rhs(:,6+ol,:) = rhs(:,6+ol,:) + &
                   ( a(3,6) * var(:,3+ol,:) + a(4,6) * var(:,4+ol,:) + &
                     a(5,6) * var(:,5+ol,:) + a(6,6) * var(:,6+ol,:) + &
                     a(7,6) * var(:,7+ol,:) + a(8,6) * var(:,8+ol,:) + &
                     a(9,6) * var(:,9+ol,:) ) * idel
  
      jl = 7 + ol
    end if
    if ( bb(4) == 0 ) then
      jr = nj - gsize(2)
    else
      or = nj - offset(4)
      rhs(:,or-5,:) = rhs(:,or-5,:) + &
                   ( a(3,6) * var(:,or-2,:) + a(4,6) * var(:,or-3,:) + &
                     a(5,6) * var(:,or-4,:) + a(6,6) * var(:,or-5,:) + &
                     a(7,6) * var(:,or-6,:) + a(8,6) * var(:,or-7,:) + &
                     a(9,6) * var(:,or-8,:) ) * idel
      rhs(:,or-4,:) = rhs(:,or-4,:) + &
                   ( a(2,5) * var(:,or-1,:) + a(3,5) * var(:,or-2,:) + &
                     a(4,5) * var(:,or-3,:) + a(5,5) * var(:,or-4,:) + &
                     a(6,5) * var(:,or-5,:) + a(7,5) * var(:,or-6,:) + &
                     a(8,5) * var(:,or-7,:) ) * idel
      rhs(:,or-3,:) = rhs(:,or-3,:) + &
                   ( a(1,4) * var(:,or,:) + a(2,4) * var(:,or-1,:) + &
                     a(3,4) * var(:,or-2,:) + a(4,4) * var(:,or-3,:) + &
                     a(5,4) * var(:,or-4,:) + a(6,4) * var(:,or-5,:) + &
                     a(7,4) * var(:,or-6,:) ) * idel
      rhs(:,or-2,:) = rhs(:,or-2,:) + &
                   ( a(1,3) * var(:,or,:) + a(2,3) * var(:,or-1,:) + &
                     a(3,3) * var(:,or-2,:) + a(4,3) * var(:,or-3,:) + &
                     a(5,3) * var(:,or-4,:) + a(6,3) * var(:,or-5,:) ) * idel
      rhs(:,or-1,:) = rhs(:,or-1,:) + &
                   ( a(1,2) * var(:,or,:) + a(2,2) * var(:,or-1,:) + &
                     a(3,2) * var(:,or-2,:) + a(4,2) * var(:,or-3,:) + &
                     a(5,2) * var(:,or-4,:) ) * idel
      rhs(:,or,:)   = rhs(:,or,:) + &
                   ( a(1,1) * var(:,or,:) + a(2,1) * var(:,or-1,:) + &
                     a(3,1) * var(:,or-2,:) + a(4,1) * var(:,or-3,:) ) * idel
  
      jr = or - 6
    end if
    rhs(:,jl:jr,:) = rhs(:,jl:jr,:) + &
                     ( -20.0_wp * var(:,jl:jr,:) + &
                        15.0_wp * ( var(:,jl-1:jr-1,:) + &
                                    var(:,jl+1:jr+1,:) ) - &
                         6.0_wp * ( var(:,jl-2:jr-2,:) + &
                                    var(:,jl+2:jr+2,:) ) + &
                                  ( var(:,jl-3:jr-3,:) + &
                                    var(:,jl+3:jr+3,:) ) ) * idel
  end if

  if ( zero_derivs_z == 0 ) then
    call set_coeff ( a )
  
    if ( scale_with_h > 0 ) then
      idel = epsilon / ( 64 * delta(3) )
    else
      idel = epsilon / 64
    end if
  
    if ( bb(5) == 0 ) then
      kl = 1 + gsize(3)
    else
      ol = offset(5)
      rhs(:,:,1+ol) = rhs(:,:,1+ol) + &
                   ( a(1,1) * var(:,:,1+ol) + a(2,1) * var(:,:,2+ol) + &
                     a(3,1) * var(:,:,3+ol) + a(4,1) * var(:,:,4+ol) ) * idel
      rhs(:,:,2+ol) = rhs(:,:,2+ol) + &
                   ( a(1,2) * var(:,:,1+ol) + a(2,2) * var(:,:,2+ol) + &
                     a(3,2) * var(:,:,3+ol) + a(4,2) * var(:,:,4+ol) + &
                     a(5,2) * var(:,:,5+ol) ) * idel
      rhs(:,:,3+ol) = rhs(:,:,3+ol) + &
                   ( a(1,3) * var(:,:,1+ol) + a(2,3) * var(:,:,2+ol) + &
                     a(3,3) * var(:,:,3+ol) + a(4,3) * var(:,:,4+ol) + &
                     a(5,3) * var(:,:,5+ol) + a(6,3) * var(:,:,6+ol) ) * idel
      rhs(:,:,4+ol) = rhs(:,:,4+ol) + &
                   ( a(1,4) * var(:,:,1+ol) + a(2,4) * var(:,:,2+ol) + &
                     a(3,4) * var(:,:,3+ol) + a(4,4) * var(:,:,4+ol) + &
                     a(5,4) * var(:,:,5+ol) + a(6,4) * var(:,:,6+ol) + &
                     a(7,4) * var(:,:,7+ol) ) * idel
      rhs(:,:,5+ol) = rhs(:,:,5+ol) + &
                   ( a(2,5) * var(:,:,2+ol) + a(3,5) * var(:,:,3+ol) + &
                     a(4,5) * var(:,:,4+ol) + a(5,5) * var(:,:,5+ol) + &
                     a(6,5) * var(:,:,6+ol) + a(7,5) * var(:,:,7+ol) + &
                     a(8,5) * var(:,:,8+ol) ) * idel
      rhs(:,:,6+ol) = rhs(:,:,6+ol) + &
                   ( a(3,6) * var(:,:,3+ol) + a(4,6) * var(:,:,4+ol) + &
                     a(5,6) * var(:,:,5+ol) + a(6,6) * var(:,:,6+ol) + &
                     a(7,6) * var(:,:,7+ol) + a(8,6) * var(:,:,8+ol) + &
                     a(9,6) * var(:,:,9+ol) ) * idel
  
      kl = 7 + ol
    end if
    if ( bb(6) == 0 ) then
      kr = nk - gsize(3)
    else
      or = nk - offset(6)
      rhs(:,:,or-5) = rhs(:,:,or-5) + &
                   ( a(3,6) * var(:,:,or-2) + a(4,6) * var(:,:,or-3) + &
                     a(5,6) * var(:,:,or-4) + a(6,6) * var(:,:,or-5) + &
                     a(7,6) * var(:,:,or-6) + a(8,6) * var(:,:,or-7) + &
                     a(9,6) * var(:,:,or-8) ) * idel
      rhs(:,:,or-4) = rhs(:,:,or-4) + &
                   ( a(2,5) * var(:,:,or-1) + a(3,5) * var(:,:,or-2) + &
                     a(4,5) * var(:,:,or-3) + a(5,5) * var(:,:,or-4) + &
                     a(6,5) * var(:,:,or-5) + a(7,5) * var(:,:,or-6) + &
                     a(8,5) * var(:,:,or-7) ) * idel
      rhs(:,:,or-3) = rhs(:,:,or-3) + &
                   ( a(1,4) * var(:,:,or) + a(2,4) * var(:,:,or-1) + &
                     a(3,4) * var(:,:,or-2) + a(4,4) * var(:,:,or-3) + &
                     a(5,4) * var(:,:,or-4) + a(6,4) * var(:,:,or-5) + &
                     a(7,4) * var(:,:,or-6) ) * idel
      rhs(:,:,or-2) = rhs(:,:,or-2) + &
                   ( a(1,3) * var(:,:,or) + a(2,3) * var(:,:,or-1) + &
                     a(3,3) * var(:,:,or-2) + a(4,3) * var(:,:,or-3) + &
                     a(5,3) * var(:,:,or-4) + a(6,3) * var(:,:,or-5) ) * idel
      rhs(:,:,or-1) = rhs(:,:,or-1) + &
                   ( a(1,2) * var(:,:,or) + a(2,2) * var(:,:,or-1) + &
                     a(3,2) * var(:,:,or-2) + a(4,2) * var(:,:,or-3) + &
                     a(5,2) * var(:,:,or-4) ) * idel
      rhs(:,:,or)   = rhs(:,:,or) + &
                   ( a(1,1) * var(:,:,or) + a(2,1) * var(:,:,or-1) + &
                     a(3,1) * var(:,:,or-2) + a(4,1) * var(:,:,or-3) ) * idel
  
      kr = or - 6
    end if
    rhs(:,:,kl:kr) = rhs(:,:,kl:kr) + &
                     ( -20.0_wp * var(:,:,kl:kr) + &
                        15.0_wp * ( var(:,:,kl-1:kr-1) + &
                                    var(:,:,kl+1:kr+1) ) - &
                         6.0_wp * ( var(:,:,kl-2:kr-2) + &
                                    var(:,:,kl+2:kr+2) ) + &
                                  ( var(:,:,kl-3:kr-3) + &
                                    var(:,:,kl+3:kr+3) ) ) * idel
  end if

contains

  subroutine set_coeff ( a )

    implicit none

    CCTK_REAL, dimension(9,6), intent(OUT) :: a
    CCTK_REAL :: zero = 0.0
    integer, parameter :: wp = kind(zero)

    a(1,1) = -3.1650670378782328375705179866656897941241116565316_wp
    a(2,1) = 9.4952011136346985127115539599970693823723349695948_wp
    a(3,1) = -9.4952011136346985127115539599970693823723349695948_wp
    a(4,1) = 3.1650670378782328375705179866656897941241116565316_wp
    a(5,1) = 0.0_wp
    a(6,1) = 0.0_wp
    a(7,1) = 0.0_wp
    a(8,1) = 0.0_wp
    a(9,1) = 0.0_wp
    a(1,2) = 2.1576625322567218846249895946058436693581952884375_wp
    a(2,2) = -7.1922084408557396154166319820194788978606509614584_wp
    a(3,2) = 8.6306501290268875384999583784233746774327811537501_wp
    a(4,2) = -4.3153250645134437692499791892116873387163905768751_wp
    a(5,2) = 0.71922084408557396154166319820194788978606509614584_wp
    a(6,2) = 0.0_wp
    a(7,2) = 0.0_wp
    a(8,2) = 0.0_wp
    a(9,2) = 0.0_wp
    a(1,3) = -4.7805237919586868314275175212098856510512725931391_wp
    a(2,3) = 19.122095167834747325710070084839542604205090372556_wp
    a(3,3) = -30.276650682405016599040944300995942456658059756547_wp
    a(4,3) = 23.902618959793434157137587606049428255256362965695_wp
    a(5,3) = -9.5610475839173736628550350424197713021025451862781_wp
    a(6,3) = 1.5935079306528956104758391737366285503504241977130_wp
    a(7,3) = 0.0_wp
    a(8,3) = 0.0_wp
    a(9,3) = 0.0_wp
    a(1,4) = 0.80612054487777570442246687814890837842881134540026_wp
    a(2,4) = -4.8367232692666542265348012688934502705728680724016_wp
    a(3,4) = 12.091808173166635566337003172233625676432170181004_wp
    a(4,4) = -16.122410897555514088449337562978167568576226908005_wp
    a(5,4) = 12.091808173166635566337003172233625676432170181004_wp
    a(6,4) = -4.8367232692666542265348012688934502705728680724016_wp
    a(7,4) = 0.80612054487777570442246687814890837842881134540026_wp
    a(8,4) = 0.0_wp
    a(9,4) = 0.0_wp
    a(1,5) = 0.0_wp
    a(2,5) = 1.0968642884346832550463374381109559476958232829758_wp
    a(3,5) = -6.5811857306080995302780246286657356861749396978545_wp
    a(4,5) = 16.452964326520248825695061571664339215437349244636_wp
    a(5,5) = -21.937285768693665100926748762219118953916465659515_wp
    a(6,5) = 16.452964326520248825695061571664339215437349244636_wp
    a(7,5) = -6.5811857306080995302780246286657356861749396978545_wp
    a(8,5) = 1.0968642884346832550463374381109559476958232829758_wp
    a(9,5) = 0.0_wp
    a(1,6) = 0.0_wp
    a(2,6) = 0.0_wp
    a(3,6) = 0.98627885208100271683294901942878016483641925983425_wp
    a(4,6) = -5.9176731124860163009976941165726809890185155590055_wp
    a(5,6) = 14.794182781215040752494235291431702472546288897514_wp
    a(6,6) = -19.725577041620054336658980388575603296728385196685_wp
    a(7,6) = 14.794182781215040752494235291431702472546288897514_wp
    a(8,6) = -5.9176731124860163009976941165726809890185155590055_wp
    a(9,6) = 0.98627885208100271683294901942878016483641925983425_wp

  end subroutine set_coeff

end subroutine dissipation_6_3
