        subroutine vfiuo( tforce, tfiuo, epsfor, dt, nrproc,
     1          velinc1, suo, svo, w1, w2, kranf1, kranf2, istart )
c
c  routine to determine velocity increments due to
c       forcing according to inegrated uhlenbeck-
c       ornstein processes.
c
c  arguments:
c       tforce  - integral time scale of forcing acceleration
c                 (it is assumed that  dt/tforce  is small enough,
c                  so that exp(dt/force) is well represented by
c                  taking only the first term)
c       tfiuo   - normalized taylor time scale of forcing
c                 ( tfiuo must be strictly less than 0.5.)
c       epsfor  - the variance of the forcing acceleration is
c                 epsfor/tforce.
c       dt      - time step
c       nrproc  - number of real-valued stochastic processes
c       velinc1  - velocity increments
c       suo     - uo processes
c       svo     - iuo processes
c       w1,w2   - work arrays
c
        dimension velinc1(nrproc),suo(nrproc),svo(nrproc),
     1          w1(nrproc),w2(nrproc)
        save t1,t2,siga,ist
c
        data ist/0/

c
c first call ------- determine time scales and initialize  ---------
c
        if(ist.eq.0) then
           ist = 1
c
           t1 = 0.5 * tforce * ( 1. + sqrt( 1. - 4.*abs(tfiuo)**2 ) )
           t2 = 0.5 * tforce * ( 1. - sqrt( 1. - 4.*abs(tfiuo)**2 ) )
           siga = sqrt( epsfor / tforce )
c
c           if( istart .eq. 0. and. tfiuo.gt.0. ) then
           if( istart .le. 0. and. tfiuo.gt.0. ) then
              call ranseq(kranf1,0)
              call raniuo( suo, svo, nrproc, 1, t1, t2, -1., w1, w2 )
           endif
        endif
c	print * ,'vfiuo: siga,dt,t1,t2=',siga,epsfor,tforce,dt,t1,t2
c
c--------------  increment iuo process  -------------------------------
c
        call ranseq(kranf2,0)
        call raniuo( suo, svo, nrproc, 1, t1, t2, dt, w1, w2 )
c
c  -------------------------------------------------------------------
c
        fmult=siga*dt
        do 100 i = 1, nrproc
        velinc1(i) = svo(i) * fmult
 100	 continue
c
c        if (istep.eq.1) then
c        write (6,*) 'vfiuo: after raniuo: velinc1(3)=',velinc1(3)
c        end if


c
        return
        end
