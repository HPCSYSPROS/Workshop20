      function e(k,kmax)
c
c specification of initial 3-d energy spectrum in waveno space
c
c this version chosen to match Moon Lee's initial conditions
c (see Lee & Reynolds pp. 226-227 and Table 4.1, p.50)
c run on 256**3 grid
c
      real k,kp,kc,k0
c
      data k0/1./
      data kp/8.04/
      data gamma/6.001e-4/
c
      kc=kmax
c
      e=0.0
c
      if (k.ge.k0.and.k.le.kp) then
      e=gamma*k*k
c
      else if (k.gt.kp.and.k.le.kc) then
      e=gamma*kp**(11./3.)*k**(-5./3.)
c
      end if
c
      return
      end
