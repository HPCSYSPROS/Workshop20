!!$/*@@
!!$  @file      Riemann1d.F90
!!$  @date      Tue Apr 29 11:42:30 2003
!!$  @author    Ian Hawke
!!$  @desc 
!!$  An old 1d Riemann solver for ideal relativistic hydro.
!!$  @enddesc 
!!$@@*/

subroutine Riemann1d

  implicit none
  
  integer :: mn, i, step, test
  real(kind=8) :: tol, pmin, pmax, dvel1, dvel2, check, gamma, t, &
       &gm1, mtol, xi, rad, vshockl, vshockr, csl, csls, csrs, &
       &csr, vels, ul, uls, urs, ur, hl, hls, hrs, hr, &
       &v12, sqrtgm1, v12ss, v12sr, eps, v1, v2, p1, p2
  real(kind=8), dimension(3) :: left, right, starl, starr, guess
  real(kind=8), dimension(5) :: waves, pwaves
  real(kind=8), dimension(4, 400) :: final
  integer :: method
  
!!$Parameters
  
  mn = 400
!!$  tol = 1.d-8
  
!!$Machine precision
  
  mtol = 1.d0
  do while (1.d0 + mtol > 1.d0)
     mtol = mtol / 2.d0
  end do
!!$  tol = mtol
  tol = 1.d-8
  eps = 1.d-6

!!$Initial Data

  t = 3.75d-1
  test = 3
  method = 3

  if (test == 1) then
    gamma = 5.d0 / 3.d0
    left(1) = 1.d0
    left(2) = 0.d0
    left(3) = 1.d2
    right(1) = 1.d0
    right(2) = 0.d0
    right(3) = 1.d0
  else if (test == 2) then
    gamma = 4.d0 / 3.d0
    left(1) = 1.d0
    left(2) = 0.9d0
    left(3) = 1.d0
    right(1) = 1.d0
    right(2) = 0.d0
    right(3) = 1.d1
  else if (test == 3) then
    gamma = 5.d0 / 3.d0
    left(1) = 1.d0
    left(2) = 0.d0
    left(3) = 1.d3
    right(1) = 1.d0
    right(2) = 0.d0
    right(3) = 1.d-2
  else if (test == 4) then
    gamma = 5.d0 / 3.d0
    left(1) = 1.d0
    left(2) = -0.6d0
    left(3) = 1.d1
    right(1) = 1.d1
    right(2) = 0.5d0
    right(3) = 2.d1
  else if (test == 5) then
    gamma = 5.d0 / 3.d0
    left(1) = 25.d0
    left(2) = 0.d0
    left(3) = 5.d1
    right(1) = 1.d0
    right(2) = 0.d0
    right(3) = 6.d0 
  else if (test == 6) then
    gamma = 5.d0 / 3.d0
    left(1) = 1.d0
    left(2) = 0.d0
    left(3) = 1.d0
    right(1) = 1.01d0
    right(2) = 0.01d0
    right(3) = 1.01d0 
  else if (test == 7) then
    gamma = 5.d0 / 3.d0
    left(1) = 1.0139158702550264d0
    left(2) = 6.57962037903012369d-3
    left(3) = 1.0309221446370552d-1
    right(1) = 1.d0
    right(2) = 0.d0
    right(3) = 1.d-1
  else if (test == 8) then !!$ shock_case = Simple
    gamma = 5.d0 / 3.d0
    left(1) = 1.d1
    left(2) = 0.d0
    left(3) = 13.333333333333333333333d0
    right(1) = 1.d0
    right(2) = 0.d0
    right(3) = 6.666666666666666d-7
  else 
    write(*,*) "Test unrecognized"
    stop
  end if
  
  gm1 = gamma - 1.d0
  
!!$Find interval in which intermediate pressure lies

!!$The loop was from when I was doing efficiency testing
  
  do i = 1, 1
    
    if (method == 1) then
      
!!$This is standard interval bisection. 
!!$First we find an interval that contains the root.

      guess = 0.5d0 * (left + right)
      
      pmin = guess(3)
      pmax = pmin
      check = 1.d0
      
      step = 0
      
      do while (check > 0.d0)
        
        pmin = 0.5d0 * max(pmin, 0.d0)
        pmax = 2.d0 * pmax
        
        call getdvel(gamma, left, right, guess, pmin, dvel1, step)
        call getdvel(gamma, left, right, guess, pmax, dvel2, step)
        
        check = dvel1 * dvel2
        
      end do
      
!!$Find the intermediate pressure - and hence velocity
      
      call getp(mtol, gamma, pmin, pmax, tol, left, right, starl, &
           &starr, guess, step)
      
    else if (method==2) then

!!$This is Newton-Raphson method
      
      guess = (0.5d0 - eps) * (left + right)
      
      call getp2(gamma, tol, left, right, starl, starr, guess, &
           &step)
      
    else

!!$This is Rezzolla-Zanotti approach
      
      guess = 0.d0
      starr = 0.d0
      starl = 0.d0
      vshockl = 0.d0
      vshockr = 0.d0
      
      sqrtgm1 = sqrt(gm1)
      
      if (left(3) > right(3)) then
        v1 = left(2)
        v2 = right(2)
        p1 = left(3)
        p2 = right(3)
        guess = left
        call getvel(gamma, right, guess, 1.d0)
        v12ss = (guess(2) - right(2)) / (1.d0 - guess(2) * right(2))
        guess = right
        call getvel(gamma, left, guess, -1.d0)
        v12sr = -(guess(2) - left(2)) / (1.d0 - guess(2) * left(2))
      else
        v2 = -left(2)
        v1 = -right(2)
        p2 = left(3)
        p1 = right(3)
        guess = right
        call getvel(gamma, left, guess, -1.d0)
        v12ss = -(guess(2) - left(2)) / (1.d0 - guess(2) * left(2))
        guess = left
        call getvel(gamma, right, guess, 1.d0)
        v12sr = (guess(2) - right(2)) / (1.d0 - guess(2) * right(2))
      end if
      
      v12 = (v1 - v2) / (1.d0 - v1 * v2)
      
      if (v12 > v12ss) then
        
        pmin = max(p1, p2)
        pmax = pmin
        check = 1.d0
        do while (check > 0.d0)
          
          pmax = 2.d0 * pmax
          
          call getdvel(gamma, left, right, guess, pmin, dvel1, step)
          call getdvel(gamma, left, right, guess, pmax, dvel2, step)
          
          check = dvel1 * dvel2
          
        end do
        
      else if (v12 > v12sr) then
        
        pmin = min(p1, p2)
        pmax = max(p1, p2)
        
      else
        
        pmin = 1.d-20
        pmax = min(p1, p2)
        
      end if
      
!!$Safety!
      
      pmin = pmin * (1.d0 - eps)
      pmax = pmax * (1.d0 + eps)
      
!!$Find the intermediate pressure - and hence velocity
      
      call getp(mtol, gamma, pmin, pmax, tol, left, right, &
           &starl, starr, guess, step)
      
    end if
    
  end do
  
  call vshock(gamma, left, starl, starr, right)
  
  ul = left(3) / gm1 / left(1)
  uls = starl(3) / gm1 / starl(1)
  urs = starr(3) / gm1 / starr(1)
  ur = right(3) / gm1 / right(1)
  
  hl = 1.d0 + ul + left(3) / left(1)
  hls = 1.d0 + uls + starl(3) / starl(1)
  hrs = 1.d0 + urs + starr(3) / starr(1)
  hr = 1.d0 + ur + right(3) / right(1)

  csl = sqrt(gamma * left(3) / left(1) / hl)
  csls = sqrt(gamma * starl(3) / starl(1) / hls)
  csrs = sqrt(gamma * starr(3) / starr(1) / hrs)
  csr = sqrt(gamma * right(3) / right(1) / hr)

  vels = 0.5d0 * (starl(2) + starr(2))

!!$Find the characteristic speeds

  if (left(3) > starl(3)) then

!!$Left rarefaction

    waves(1) = (left(2) - csl) / (1.d0 - left(2) * csl)
    waves(2) = (starl(2) - csls) / (1.d0 - starl(2) * csls)

  else

!!$Left shock

    waves(1) = vshockl
    waves(2) = waves(1)

  end if

  waves(3) = vels
  
  if (right(3) > starl(3)) then

!!$Right rarefaction

    waves(4) = (starr(2) + csrs) / (1.d0 + starr(2) * csrs)
    waves(5) = (right(2) + csr) / (1.d0 + right(2) * csr)

  else

!!$Right shock

    waves(4) = vshockr
    waves(5) = waves(4)

  end if

  do i = 1, 2
    if (abs(waves(i)) > 1.d10) waves(i) = -1.d10
  end do
  do i = 4, 5
    if (abs(waves(i)) > 1.d10) waves(i) = 1.d10
  end do

  pwaves = 0.5d0 + t * waves

  do i = 1, mn

    rad = real(i) / real(mn)
    final(1, i) = rad

    if (rad < pwaves(1)) then

      final(2, i) = left(1)
      final(3, i) = left(2)
      final(4, i) = left(3)
      
    else if (rad < pwaves(2)) then
      
      xi = (rad - 0.5d0) / t
      
      call raref(gamma, xi, left, final(:, i), 1)
      
    else if (rad < pwaves(3)) then
      
      final(2, i) = starl(1)
      final(3, i) = starl(2)
      final(4, i) = starl(3)
      
    else if (rad < pwaves(4)) then
      
      final(2, i) = starr(1)
      final(3, i) = starr(2)
      final(4, i) = starr(3)
      
    else if (rad < pwaves(5)) then
      
      xi = (rad - 0.5d0) / t
      
      call raref(gamma, xi, right, final(:, i), -1)
      
    else
      
      final(2, i) = right(1)
      final(3, i) = right(2)
      final(4, i) = right(3)
      
    end if
    
  end do
  
  open(10, file='Exact_rho.dat')
  
  do i = 1, mn
    write(10, '(2e20.12)') final(1, i), final(2, i)
  end do
  
  close(10)

  open(11, file='Exact_velx.dat')
  
  do i = 1, mn
    write(11, '(2e20.12)') final(1, i), final(3, i)
  end do
  
  close(11)

  open(12, file='Exact_press.dat')
  
  do i = 1, mn
    write(12, '(2e20.12)') final(1, i), final(4, i)
  end do
  
  close(12)

contains
  
!!$Given a right and left state and a guess of the intermediate pressure, 
!!$returns the difference between the intermediate velocities connecting 
!!$right and left states. 

  subroutine getdvel(gamma, left, right, guess, pguess, dvel, step)

    implicit none

    real(kind=8), dimension(3) :: left, right, guess
    real(kind=8) :: dvel, pguess, gamma
    integer :: step
    
    step = step + 1
    
    guess(3) = pguess

!!$Left wave
      
    call getvel(gamma, left, guess, -1.d0)

    dvel = guess(2)

!!$Right wave

    call getvel(gamma, right, guess, 1.d0)

    dvel = dvel - guess(2)

  end subroutine getdvel

!!$Given a prewave state (known) and a postwave state (guess, with pressure
!!$and density guessed) returns with the postwave velocity. Sign indicates
!!$left (-) or right (+) wave

  subroutine getvel(gamma, known, guess, sign)

    implicit none

    real(kind=8), dimension(3) :: known, guess
    real(kind=8) :: gamma, gm1, cs, sign, a, b, c, h, e, ea, ua,&
         &ha, csa, k, sqgl1, v12
    
    gm1 = gamma - 1.d0
    ua = known(3) / gm1 / known(1)
    ea = known(1) + known(3) / gm1  
    ha = 1.d0 + ua + known(3) / known(1)
    csa = sqrt(gamma * known(3) / known(1) / ha) 
    
    if (guess(3) > known(3)) then
      
!!$Shock
      
      a = 1.d0 + gm1 * (known(3) - guess(3)) / gamma / guess(3)
      b = 1.d0 - a
      c = ha * (known(3) - guess(3)) / known(1) - ha * ha
      
      if (c > b * b / 4.d0 / a) then
        write (*,*) "unphysical enthalpy"
        stop
      endif
      
      h = (-b + sqrt(b * b - 4.d0 * a * c)) / 2.d0 / a
      
      guess(1) = gamma * guess(3) / gm1 / (h - 1.d0)
      e = guess(1) + guess(3) / gm1
      
      v12 = -sign * sqrt((guess(3) - known(3)) * (e - ea) / &
           &(ea + guess(3)) / (e + known(3)))
      
!!$         pp = guess(3) - known(3)
!!$         n = pp * ( guess(3) * (ha - 1.d0) * (gm1 + h) - &
!!$              &known(3) * (h - 1.d0) * (gm1 + ha) )
!!$         d = (ha - 1.d0) * (ea + guess(3)) * ( (gm1 + h) * pp + &
!!$              &gamma * h * known(3) )
!!$
!!$         v12 = -sign * sqrt(n / d)

      guess(2) = (known(2) - v12) / (1.d0 - known(2) * v12)

    else

!!$Rarefaction

      k = known(3) / known(1)**gamma
      
      guess(1) = (guess(3) / k)**(1.d0 / gamma)
      cs = sqrt(gamma * guess(3) / (guess(1) + gamma * guess(3) / gm1))
      sqgl1 = sqrt(gm1)
      a = (1.d0 + known(2)) / (1.d0 - known(2)) * &
           &((sqgl1 + csa) / (sqgl1 - csa) * (sqgl1 - cs) / &
           &(sqgl1 + cs))**(-sign * 2.d0 / sqgl1)
      guess(2) = (a - 1.d0) / (a + 1.d0)
      
    end if
    
  end subroutine getvel

!!$Given an interval, computes the correct intermediate pressure and 
!!$returns the intermediate state in star
    
  subroutine getp(mtol, gamma, pmin, pmax, tol, left, right, starl, &
       &starr, guess, step)
    
    implicit none
    
    real(kind=8), dimension(3) :: left, right, guess, starl, starr
    real(kind=8) :: pmin, pmax, tol, gamma, mtol, a, b, c, d, e, &
         &fa, fb, fc, tol1, xm, p, q, r, s
    logical :: flag
    integer :: step
    
    a = pmin
    b = pmax
    call getdvel(gamma, left, right, guess, a, fa, step)
    call getdvel(gamma, left, right, guess, b, fb, step)
    
    do
      
      c = a
      fc = fa
      d = b - a
      e = d
      
      flag = .false.
      
      do
        
        if (abs(fc) < abs(fb)) then
          a = b
          b = c
          c = a
          fa = fb
          fb = fc
          fc = fa
        end if

!!$Convergence test
        
        tol1 = 2.d0 * mtol * abs(b) + 0.5d0 * tol
        xm = 0.5d0 * (c - b)
        
        if ((abs(xm).le.tol1).or.(fb == 0.d0)) then
          flag = .true.
          exit
        end if
        
!!$Is bisection necessary?
        
        if ( (abs(e) .ge. tol1).and.(abs(fa) > abs(fb)) ) then
          
!!$If quadratic interpolation is impossible (a =  c) then do linear
          
          if (a == c) then
            
            s = fb / fa
            p = 2.d0 * xm * s
            q = 1.d0 - s
            
          else
            
            q = fa / fc
            r = fb / fc
            s = fb / fa
            p = s * ( 2.d0 * xm * q * (q - r) - (b - a) * (r - 1.d0) )
            q = (q - 1.d0) * (r - 1.d0) * (s - 1.d0)
            
          end if

!!$Adjust the signs
            
          if (p > 0.d0) q = -q
          p = abs(p)
          
!!$Is the interpolation acceptable?

          if ( (2.d0 * p < 3.d0 * xm * q - abs(tol1 * q)).and.&
               &(p < abs(0.5d0 * e * q)) ) then
            
            e = d
            d = p / q
            
          else
            
!!$Bisection
                  
            d = xm
            e = d
                  
          end if

        else

          d = xm
          e = d
          
        end if

!!$Complete the step
            
        a = b
        fa = fb
        
        if (abs(d) > tol1) then
          
          b = b + d
          
        else
          
          b = b + sign(tol1, xm)
          
        end if
        
        call getdvel(gamma, left, right, guess, b, fb, step)
        
        if ( fb * fc / abs(fc) > 0.d0 ) exit
        
      end do
      
      if (flag) exit
      
    end do

    starl(3) = b
    starr(3) = b
      
    guess = starl
      
    call getvel(gamma, left, guess, -1.d0)
    
    starl(1:2) = guess(1:2)
    
    guess = starr

    call getvel(gamma, right, guess, 1.d0)

    starr(1:2) = guess(1:2)

  end subroutine getp

  subroutine raref(gamma, xi, known, out, sign)

    implicit none
    
    real(kind=8), intent(in) :: gamma, xi
    real(kind=8), dimension(3) :: known
    real(kind=8), dimension(4) :: out
    integer, intent(in) :: sign
    real(kind=8) :: gm1, b, c, d, k, l, v, ua, ha, csa, cs2, &
         &ocs2, fcs2, dfdcs2, tol
    
    gm1 = gamma - 1.d0
    tol = 1.d-8
    ua = known(3) / gm1 / known(1)
    ha = 1.d0 + ua + known(3) / known(1)
    csa = sqrt(gamma * known(3) / known(1) / ha)
    b = sqrt(gm1)
    c = (b + csa) / (b - csa)
    d = -sign * b / 2.d0
    k = (1.d0 + xi) / (1.d0 - xi)
    l = c * k**d
    v = ((1.d0 - known(2)) / (1.d0 + known(2)))**d
    
    ocs2 = csa
    
    do
      
      fcs2 = l * v * (1.d0 + sign * ocs2)**d * (ocs2 - b) + &
           &(1.d0 - sign * ocs2)**d * (ocs2 + b)
      dfdcs2 = l * v * (1.d0 + sign * ocs2)**d * &
           &(1.d0 + sign * d * (ocs2 - b) / (1.d0 + sign * ocs2)) + &
           &(1.d0 - sign * ocs2)**d * (1.d0 - sign * d * (ocs2 + b) / &
           &(1.d0 - sign * ocs2)) 
      
      cs2 = ocs2 - fcs2 / dfdcs2
      
      if (abs(cs2 - ocs2) / ocs2 < tol) then
        
        exit
        
      else
        
        ocs2 = cs2
        
      end if
      
    end do
    
    out(3) = (xi + sign * cs2) / (1.d0 + sign * xi * cs2)
    out(2) = known(1) * ((cs2 * cs2 * (gm1 - csa * csa)) / &
         &(csa * csa * (gm1 - cs2 * cs2)))**(1.d0/gm1)
    out(4) = cs2 * cs2 * gm1 * out(2) / (gm1 - cs2 * cs2) / gamma
    
  end subroutine raref
  
!!$Given a right and left state and a guess of the intermediate pressure, 
!!$returns the difference between the intermediate velocities connecting 
!!$right and left states. 

  subroutine getdvel2(gamma, left, right, guess, pguess, dvel, step, &
       &dfdp)
    
    implicit none
    
    real(kind=8), dimension(3) :: left, right, guess
    real(kind=8) :: dvel, pguess, gamma, &
         &dfdp, dfdp1
    integer :: step
    
    step = step + 1
    
    guess(3) = pguess
    
!!$Left wave
    
    call getvel2(gamma, left, guess, -1.d0, dfdp1)
    
    dvel = guess(2)
    dfdp = dfdp1
    
!!$Right wave
    
    call getvel2(gamma, right, guess, 1.d0, dfdp1)
    
    dvel = dvel - guess(2)
    dfdp = dfdp - dfdp1
    
  end subroutine getdvel2
  
!!$Given a prewave state (known) and a postwave state (guess, with pressure
!!$and density guessed) returns with the postwave velocity. Sign indicates
!!$left (-) or right (+) wave
  
  subroutine getvel2(gamma, known, guess, sign, dfdp)

    implicit none

    real(kind=8), dimension(3) :: known, guess
    real(kind=8) :: gamma, gm1, cs, sign, a, b, c, h, ua, &
         &ha, csa, wa, k, sqgl1, dfdp, dhdp, &
         &dadp, e, dbdp, dcdp, &
         &dv12dp, v12, dedp, ea

    gm1 = gamma - 1.d0
    ua = known(3) / gm1 / known(1)
    ea = known(1) + known(3) / gm1 
    ha = 1.d0 + ua + known(3) / known(1)
    csa = sqrt(gamma * known(3) / known(1) / ha)
    wa = 1.d0 / sqrt(1.d0 - known(2) * known(2)) 

    if (guess(3) > known(3)) then

!!$Shock

      a = 1.d0 + gm1 * (known(3) - guess(3)) / gamma / guess(3)
      b = 1.d0 - a
      c = ha * (known(3) - guess(3)) / known(1) - ha * ha

      if (c > b * b / 4.d0 / a) then
        write(*,*) "unphysical enthalpy"
        stop
      endif

      h = (-b + sqrt(b * b - 4.d0 * a * c)) / 2.d0 / a

      dadp = -gm1 / gamma * known(3) / guess(3) / guess(3)
      dbdp = - dadp
      dcdp = - ha / known(1)
      dhdp = - (h * h * dadp + h * dbdp + dcdp) / (2.d0 * a * h + b)

      guess(1) = gamma * guess(3) / gm1 / (h - 1.d0)
      e = guess(1) + guess(3) / gm1

      dedp = gamma / gm1 * (h - 1.d0 - guess(3) * dhdp) / &
           &(h - 1.d0) / (h - 1.d0) + 1.d0 / gm1

      v12 = -sign * sqrt((guess(3) - known(3)) * (e - ea) / &
           &(ea + guess(3)) / (e + known(3)))

      dv12dp = (ea + known(3)) / 2.d0 / (ea + guess(3)) / &
           &(e + known(3)) * ( (e + known(3)) / (guess(3) - known(3)) + &
           &dedp * (ea + guess(3)) / (e - ea) )

      guess(2) = (known(2) - v12) / (1.d0 - known(2) * v12)

      dfdp = dv12dp * (known(2) * known(2) - 1.d0) / &
           &(known(2) * v12 - 1.d0)

    else

!!$Rarefaction

      k = known(3) / known(1)**gamma

      guess(1) = (guess(3) / k)**(1.d0 / gamma)
      cs = sqrt(gamma * guess(3) / (guess(1) + gamma * guess(3) / gm1))
      sqgl1 = sqrt(gm1)
      a = (1.d0 + known(2)) / (1.d0 - known(2)) * &
           &((sqgl1 + csa) / (sqgl1 - csa) * (sqgl1 - cs) / &
           &(sqgl1 + cs))**(-sign * 2.d0 / sqgl1)
      guess(2) = (a - 1.d0) / (a + 1.d0)

      dfdp = sign * cs / guess(3) * (1.d0 - guess(2) * guess(2))

    end if

  end subroutine getvel2

!!$Given an interval, computes the correct intermediate pressure and 
!!$returns the intermediate state in star
    
  subroutine getp2(gamma, tol, left, right, starl, &
       &starr, guess, step)

    implicit none

    real(kind=8), dimension(3) :: left, right, guess, starl, starr
    real(kind=8) :: tol, gamma, f, dfdp
    integer :: step

    step = 0

    if (guess(3) > 1.d-8) guess(3) = guess(3) - 1.d-8

    do

      if (step > 200) stop

      call getdvel2(gamma, left, right, guess, guess(3), f, &
           &step, dfdp)

      if (abs(f) < tol) then

        exit

      else

        guess(3) = guess(3) - f / dfdp  

        if (guess(3) < 1.d-20) guess(3) = 1.d-20

      end if

    end do

    starl(3) = guess(3)
    starr(3) = guess(3)

    guess = starl

    call getvel2(gamma, left, guess, -1.d0, dfdp)

    starl(1:2) = guess(1:2)

    guess = starr

    call getvel2(gamma, right, guess, 1.d0, dfdp)

    starr(1:2) = guess(1:2)

  end subroutine getp2

!!$Quick routine to find the shock velocities

  subroutine vshock(gamma, l, sl, sr, r)

    implicit none

    real(kind=8), dimension(3) :: l, sl, sr, r
    real(kind=8) :: h, j, gamma, gm1, &
         &ha, wa2, a, b

    gm1 = gamma - 1.d0
    vshockl = 0.d0
    vshockr = 0.d0

    if (sl(3) > l(3)) then

      wa2 = 1.d0 / (1.d0 - l(2) * l(2))
      ha = 1.d0 + l(3) / l(1) * gamma / gm1
      h = 1.d0 + sl(3) / sl(1) * gamma / gm1
      j = - sqrt((sl(3) - l(3)) / (ha / l(1) - h / sl(1)))
      a = j * j + l(1) * l(1) * wa2
      b = - l(2) * l(1) * l(1) * wa2
      vshockl = (-b - j * j * sqrt(1.d0 + l(1) * l(1) / j / j)) / a

    end if

    if (sr(3) > r(3)) then

      wa2 = 1.d0 / (1.d0 - r(2) * r(2))
      ha = 1.d0 + r(3) / r(1) * gamma / gm1
      h = 1.d0 + sr(3) / sr(1) * gamma / gm1
      j = sqrt((sr(3) - r(3)) / (ha / r(1) - h / sr(1)))
      a = j * j + r(1) * r(1) * wa2
      b = - r(2) * r(1) * r(1) * wa2
      vshockr = (-b + j * j * sqrt(1.d0 + r(1) * r(1) / j / j)) / a

    end if

  end subroutine vshock
    
end subroutine Riemann1d
