      FUNCTION E(k, kmax)
#ifdef MODEL_SPECTRUM
C
C E: Specification of initial 3D energy spectrum.
C
C The initial energy spectrum follows the model spectrum in "Turbulent
C Flow" (Pope):
C
C        E(k) = C*eps**(2/3)*k**(-5/3)*fL(k*L)*fn(k*eta),
C
C where k is the wavenumber, C is the Kolmogorov constant, eps is the
C dissipation rate, L is the integral scale, and eta is the Kolmogorov
C microscale. The large-scale function fL is
C
C        fL(k*L) = (k*L/((k*L)**2 + cL)**0.5)**(5/3 + p0),
C
C where cL and p0 are model parameters. The microscale function fn is
C
C        fn(k*eta) = EXP(-beta*(((k*eta)**4 + ceta**4)**0.25 - ceta)),
C
C where beta and ceta are model parameters.
C
C Authors:
C              Matthew Clay and PK Yeung
C              School of Aerospace Engineering
C              Georgia Institute of Technology
C
C References:
C
C 1. Pope, S. B., "Turbulent Flows," Cambridge University Press, 2000.
C
      USE spectrum,ONLY: p0, beta, ceta, cL, intlen, eta, eps, kol
#ifdef SPECTRAL_TRUNCATION
      USE spectrum,ONLY: kLower, kUpper
#endif
      USE com,ONLY: viscos
      IMPLICIT NONE
      REAL :: k, kc, E, kL, fkL, keta, arg, fketa
      INTEGER :: kmax
C
      kc = kmax
      E = 0.0
      IF (k .LE. kc) THEN
C Integral scale function.
         kL = k*intlen
         fkL = (kL/SQRT(kL*kL + cL))**(5./3. + p0)
C Microscale function.
         keta = k*eta
         arg = (keta**4 + ceta**4)**0.25 - ceta
         fketa = EXP(-beta*arg)
C Energy spectrum function.
         E = kol*eps**(2./3.)*k**(-5./3.)*fkL*fketa
      END IF
#ifdef SPECTRAL_TRUNCATION
C Truncate wavenumbers so each direction uses the same range of physical
C wavenumbers.
      IF (k .LT. kLower) THEN
         E = 0.0
      END IF
      IF (k .GT. kUpper) THEN
         E = 0.0
      END IF
#endif
#endif
      RETURN
      END
