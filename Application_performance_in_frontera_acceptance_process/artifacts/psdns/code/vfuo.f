        subroutine vfuo( tforce, epsfor, dt, nrproc,
     1          velinc1, suo, w1, w2, kranf1, kranf2, istart )
c
c  routine to determine velocity increments due to
c       forcing according to uhlenbeck-ornstein processes.
c
c  arguments:
c       tforce  - integral time scale of forcing acceleration
c       epsfor  - the variance of the forcing acceleration is
c                 epsfor/tforce.
c       dt      - time step
c       nrproc  - number of real-valued stochastic processes
c       velinc1  - velocity increment
c       suo     - uo process
c       w1,w2   - work arrays
c
        dimension velinc1(nrproc),suo(nrproc),w1(nrproc),w2(nrproc)
        real*8 exp1,exp2,rana,ranua,ranub,ranud,dzero
c
        data ist/0/, dzero/0./
c
c
c------ first call ------- initialize uo process ----------------------
c
        if(ist.eq.0) then
           ist = 1
           if( istart .eq. 0 ) then
              call ranseq(kranf1,0)
              call rann2( suo , 1 , nrproc , 1 , 1 , nrproc )
           endif
        endif
c
c  ------  determine u-o parameters  ---------------------------------
c
        exp1 = 0.
        if( tforce .gt. 0.05*dt ) exp1 = exp( -dt / tforce )
        exp2 = exp1**2
c
        rana  = sqrt(1.-exp2)
        ranua = sqrt( epsfor*tforce ) * (1.-exp1)**2 / sqrt(1.-exp2)
        ranub =     ( epsfor*(  2.*dt  +
     1              tforce*( - 3. + 4.*exp1 - exp2 )  ) - ranua**2   )
c       ranub =dsqrt (dmax1(ranub,dzero))
        ranub =sqrt (max(ranub,dzero))
        ranud = sqrt( epsfor*tforce ) * ( 1. - exp1 )
c
c  generate gaussian random numbers ----------------------------------
c
        call ranseq(kranf2,0)
        call rann2( w1 , 1 , nrproc , 1 , 1 , nrproc )
        call rann2( w2 , 1 , nrproc , 1 , 1 , nrproc )
c
        do 100 i = 1, nrproc
        velinc1(i) = ranud * suo(i) + ranua * w1(i) + ranub * w2(i)
100     suo(i) = exp1 * suo(i) + rana * w1(i)
c
        return
        end
