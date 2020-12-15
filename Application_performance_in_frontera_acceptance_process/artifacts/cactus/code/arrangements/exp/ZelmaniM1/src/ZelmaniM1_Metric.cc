#include "ZelmaniM1_Metric.hh"
#include <math.h>

using namespace ThreeTensor;
  
Metric::Metric(double gxxi, double gxyi, double gxzi,
            double gyyi, double gyzi, double gzzi,
            double alp, double bx, double by, double bz,
	    double Kxx, double Kxy, double Kxz, 
	    double Kyy, double Kyz, double Kzz) :
            gxx(gxxi), gxy(gxyi), gxz(gxzi), gyy(gyyi), gyz(gyzi), gzz(gzzi), 
            alpha(alp), betax(bx), betay(by), betaz(bz), 
	    Kxx(Kxx), Kxy(Kxy), Kxz(Kxz), Kyy(Kyy), Kyz(Kyz), Kzz(Kzz) {
  detg = (gxx*gyy*gzz + 2.0*gxy*gxz*gyz 
    - gyy*gxz*gxz - gxx*gyz*gyz - gzz*gxy*gxy);
  const double idenom = 1.0/detg;
  guxx = (gyy*gzz - gyz*gyz)*idenom;
  guxy = (gxz*gyz - gxy*gzz)*idenom;
  guxz = (gxy*gyz - gxz*gyy)*idenom;
  guyy = (gxx*gzz - gxz*gxz)*idenom;
  guyz = (gxy*gxz - gxx*gyz)*idenom;
  guzz = (gxx*gyy - gxy*gxy)*idenom;
}

void Metric::loadLapseShift(double alphai, double betaxi){
  alpha = alphai; betax = betaxi;
  betay = 0.0;    betaz = 0.0;
}
      
void Metric::loadLapseShiftAll(double alphai, double betaxi, 
                       double betayi, double betazi){
  alpha = alphai; betax = betaxi;
  betay = betayi; betaz = betazi; 
}

void Metric::loadExtCurv(double Kxxi, double Kxyi, double Kxzi, 
                 double Kyyi, double Kyzi, double Kzzi) {
  Kxx = Kxxi; Kxy = Kxyi; Kxz = Kxzi;
  Kyy = Kyyi; Kyz = Kyzi;
  Kzz = Kzzi;
  return;
}
      
double Metric::contractKllTuuSym(double Txx, double Txy, double Txz,
                         double Tyy, double Tyz, double Tzz) const {
  return        Txx*Kxx + Tyy*Kyy + Tzz*Kzz 
         + 2.0*(Txy*Kxy + Txz*Kxz + Tyz*Kyz); 
}


// Contract the three metric with three vectors 
double Metric::contractAlBu(double Ax, double Ay, double Az,
    double Bx, double By, double Bz) const {
  return Ax*Bx + Ay*By + Az*Bz;
}
double Metric::contractGllAuBu(double Ax, double Ay, double Az,
    double Bx, double By, double Bz) const {
  return  gxx*Ax*Bx + gxy*Ax*By + gxz*Ax*Bz 
        + gxy*Ay*Bx + gyy*Ay*By + gyz*Ay*Bz 
        + gxz*Az*Bx + gyz*Az*By + gzz*Az*Bz;
}

double Metric::contractGuuAlBl(double Ax, double Ay, double Az,
    double Bx, double By, double Bz) const {
  return  guxx*Ax*Bx + guxy*Ax*By + guxz*Ax*Bz 
        + guxy*Ay*Bx + guyy*Ay*By + guyz*Ay*Bz 
        + guxz*Az*Bx + guyz*Az*By + guzz*Az*Bz;
}

// Raise and lower three vector indices 
void Metric::lowerAu(double Axu, double Ayu, double Azu,
    double *Axl, double *Ayl, double *Azl) const {
  *Axl = gxx*Axu + gxy*Ayu + gxz*Azu;
  *Ayl = gxy*Axu + gyy*Ayu + gyz*Azu;
  *Azl = gxz*Axu + gyz*Ayu + gzz*Azu;
}

void Metric::raiseAl(double Axl, double Ayl, double Azl,
    double *Axu, double *Ayu, double *Azu) const {
  *Axu = guxx*Axl + guxy*Ayl + guxz*Azl;
  *Ayu = guxy*Axl + guyy*Ayl + guyz*Azl;
  *Azu = guxz*Axl + guyz*Ayl + guzz*Azl;
}

void Metric::raiseSAll(double  Axxl, double  Axyl, double  Axzl, 
  double  Ayyl, double  Ayzl, double  Azzl,
  double *Axxu, double *Axyu, double *Axzu, 
  double *Ayyu, double *Ayzu, double *Azzu) const {
  
  *Axxu =   guxx*guxx*Axxl + 2.0*guxx*(guxy*Axyl + guxz*Axzl)
          + guxy*guxy*Ayyl + 2.0*guxy* guxz*Ayzl
          + guxz*guxz*Azzl; 

  *Axyu =   guxy*guxx*Axxl + 2.0*guxy*(guxy*Axyl + guxz*Axzl)
          + guyy*guxy*Ayyl + 2.0*guyy* guxz*Ayzl
          + guyz*guxz*Azzl; 
  
  *Axzu =   guxz*guxx*Axxl + 2.0*guxz*(guxy*Axyl + guxz*Axzl)
          + guxz*guxy*Ayyl + 2.0*guyz* guxz*Ayzl
          + guzz*guxz*Azzl; 

  *Ayyu =   guxy*guxy*Axxl + 2.0*guxy*(guyy*Axyl + guyz*Axzl)
          + guyy*guyy*Ayyl + 2.0*guyy* guyz*Ayzl
          + guyz*guyz*Azzl; 

  *Ayzu =   guxz*guxy*Axxl + 2.0*guxz*(guyy*Axyl + guyz*Axzl)
          + guyz*guyy*Ayyl + 2.0*guyz* guyz*Ayzl
          + guzz*guyz*Azzl; 

  *Azzu =   guxz*guxz*Axxl + 2.0*guxz*(guyz*Axyl + guzz*Axzl)
          + guyz*guyz*Ayyl + 2.0*guyz* guzz*Ayzl
    + guzz*guzz*Azzl; 

        return;
}

void Metric::lowerSAuu(double  Axxu, double  Axyu, double  Axzu, 
    double  Ayyu, double  Ayzu, double  Azzu,
    double *Axxl, double *Axyl, double *Axzl, 
    double *Ayyl, double *Ayzl, double *Azzl) const {
  
  *Axxl =   gxx*gxx*Axxu + 2.0*gxx*(gxy*Axyu + gxz*Axzu)
          + gxy*gxy*Ayyu + 2.0*gxy*gxz*Ayzu
    + gxz*gxz*Azzu; 

  *Axyl =   gxy*gxx*Axxu + 2.0*gxy*(gxy*Axyu + gxz*Axzu)
          + gyy*gxy*Ayyu + 2.0*gyy*gxz*Ayzu
    + gyz*gxz*Azzu; 
  
  *Axzl =   gxz*gxx*Axxu + 2.0*gxz*(gxy*Axyu + gxz*Axzu)
          + gxz*gxy*Ayyu + 2.0*gyz*gxz*Ayzu
    + gzz*gxz*Azzu; 

  *Ayyl =   gxy*gxy*Axxu + 2.0*gxy*(gyy*Axyu + gyz*Axzu)
          + gyy*gyy*Ayyu + 2.0*gyy*gyz*Ayzu
    + gyz*gyz*Azzu; 

  *Ayzl =   gxz*gxy*Axxu + 2.0*gxz*(gyy*Axyu + gyz*Axzu)
          + gyz*gyy*Ayyu + 2.0*gyz*gyz*Ayzu
    + gzz*gyz*Azzu; 

  *Azzl =   gxz*gxz*Axxu + 2.0*gxz*(gyz*Axyu + gzz*Axzu)
          + gyz*gyz*Ayyu + 2.0*gyz*gzz*Ayzu
    + gzz*gzz*Azzu; 

        return;
}
  
