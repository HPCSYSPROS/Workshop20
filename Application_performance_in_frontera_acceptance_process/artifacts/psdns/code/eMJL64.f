      real function e(k,kmax)
c
c specification of initial 3-d energy spectrum in waveno space
c
c this version chosen to match Moon Lee's initial conditions
c (see Lee & Reynolds pp. 226-227 and Table 4.1, p.50)
c run on 64**3 grid
c
      real k,kp,kc,k0
      data k0/1./  
c     data kp/7.04/
c     data gamma/1.526e-3/
      data kp/8.08/
      data gamma/7.888e-4/

      kc=kmax

      e=0.0

      if (k.ge.k0.and.k.le.kp) then
         e=gamma*k*k 
      else if (k.gt.kp.and.k.le.kc) then
         e=gamma*kp**(11./3.)*k**(-5./3.)
      end if

      return
      end
