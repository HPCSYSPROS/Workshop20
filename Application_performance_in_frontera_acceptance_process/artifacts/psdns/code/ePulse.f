      function e(k,kmax)
c
c specification of initial 3-d energy spectrum in waveno space
c
c pulse spectrum similar to that used for R128 (Rogers')
c
      real k,kp,kc,k0
c
      data k0/1./  
c
      kc=kmax
c
      e=0.0
c
c     if (k.ge.16.and.k.le.kmax) e=1.
c     if (k.ge.8.and.k.le.16.) e=1.
c     if (k.ge.8.and.k.le.12.) e=1.
      if (k.ge.4.and.k.le.8.) e=1.
c     if (k.ge.1.and.k.le.2.) e=1.
c     if (k.ge.2.and.k.le.4.) e=1.
c     if (k.ge.16.and.k.le.32.) e=1.
c
      return
      end
