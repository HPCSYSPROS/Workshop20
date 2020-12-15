        subroutine raniuo(u,v,n,istr,t1,t2,dt,w1,w2)
c
c  routine to initialize and increment a standardized integrated
c       uhlenbeck-ornstien process.
c
c   u   - vector containing uo processes
c   v   - vector containing iuo processes
c   n   - number of processes
c   istr- vector stride of u and v
c   t1,t2 -  time scales
c       t1 is the time-scale of the uo process
c       t2 is the integration time scale
c       the integral time scale of the iuo process is t1+t2
c       the taylor time scale of the iuo process is sqrt(t1*t2)
c       t1 and t2 must both be positive, and should not be equal
c       (if t1 and t2 are commuted, the statistics of v are unaltered,
c        but the joint statistics of u and v are.)
c   dt - time increment
c   w1,w2 - work arrays of dimension at least n
c
c  notes:
c   to initialize, set dt < 0.
c
c  other routines called: - rann2
c
#ifdef RANDOM_SEED
        use ranseed
#endif

        dimension u(1),v(1),w1(n),w2(n)
c
      data ist/0/
      ist=ist+1
c
c  check time scales
c
        if( t1 .le. 0. ) then
           write(6,*)' raniuo: t1=', t1, 'must be positive. stopped'
           stop
        endif
c
        if( t2 .le. 0. ) then
           write(6,*)' raniuo: t2=', t2, 'must be positive. stopped'
           stop
        endif
c
        diff = abs( t1/t2 -1. )
        if( diff .lt. 0.05 ) then
           write(6,*)' raniuo: t1=', t1, 'too close to t2=', t2
           write(6,*)' t1 decreased to 0.95 t2 '
           t1 = 0.95 * t2
        endif
c
c   initialization
c
        if( dt .lt. 0. ) then
c
c  determine initial correlation
c
           rr= t1 / ( t1 + t2 )
           be=sqrt(rr)
           ga=sqrt(1.-rr)
c
        write (6,*) 'raniuo: iseed=',iseed
           call rann2(w1,1,n,1,1,n)
           call rann2(w2,1,n,1,1,n)
c
           ii = 1 - istr
           do 10 i = 1 , n
           ii = ii + istr
           u(ii)=w1(i)
10         v(ii)=be*w1(i)+ga*w2(i)
c
           return
        endif
c
c evaluate parameters
c
        om1=1./t1
        om2=1./t2
        s=om1+om2
        d=om2-om1
c
        a=exp(-om1*dt)
        b=exp(-om2*dt)
        c=(a-b)/d
        sco=c*sqrt(s*om2)
        r=sqrt(om2/s)
c
        uu=1.-a*a
        vv=1.-b*b+4.*om1*om2*b*c/d-om2*s*(a*a-b*b)/(d*d)
        uv=r*(1.-a*b-s*a*c)
c
        al=sqrt(uu)
        be=uv/al
        ga=sqrt( amax1( vv - be*be ,0. ) )
c
c increment iuo process
c
        call rann2(w1,1,n,1,1,n)
        call rann2(w2,1,n,1,1,n)
c
      if (ist.eq.2) write (6,*) 't1,t2=',t1,t2
        ii = 1 - istr
        do 100 i = 1 , n
        ii = ii + istr
        v(ii)=b*v(ii)+sco*u(ii)+be*w1(i)+ga*w2(i)
        u(ii)=a*u(ii)+al*w1(i)
 100    continue
c
        return
        end
