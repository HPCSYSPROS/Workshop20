#ifndef _THREEMETRICCLASS_
#define _THREEMETRICCLASS_

#include <math.h>

namespace ThreeTensor {

class Vector {
public:  
  double x,y,z;
  bool upper; 
  Vector(double x, double y, double z, bool up):x(x),y(y),z(z),upper(up){};
  double operator[](int i){
    int idx = i % 3; 
    if (idx==0) return x;
    if (idx==1) return y;
    return z;
  }
};

class Sym2Tensor { 
public:
  double xx,xy,xz;
  double yy,yz,zz;
  bool upper; 
  Sym2Tensor(double xx, double xy, double xz, 
             double yy, double yz, double zz,bool up): 
             xx(xx),xy(xy),xz(xz),yy(yy),yz(yz),zz(zz),upper(up){};
};

// This is a storage capsule for local properties of the spacetime that are useful for radiative transfer
class Metric {

protected:  
  // Input quantities
  double detg; 
   
public:
  double gxx,gxy,gxz,gyy,gyz,gzz; //gamma_{ij}
  double alpha,betax,betay,betaz; //Lapse and shift (raised index)
  double guxx,guxy,guxz,guyy,guyz,guzz; //gamma^{ij}
  double Kxx,Kxy,Kxz,Kyy,Kyz,Kzz; //K_{ij}
  
  Metric(double gxx = 1.0, double gxy = 0.0, double gxz = 0.0, 
      double gyy = 1.0, double gyz = 0.0, double gzz = 1.0,
      double alp = 1.0, 
      double bx = 0.0, double by = 0.0, double bz = 0.0,
      double Kxx = 0.0, double Kxy = 0.0, double Kxz = 0.0,
      double Kyy = 0.0, double Kyz = 0.0, double Kzz = 0.0);
 
  void loadLapseShift(double alphai, double betaxi);
        
  void loadLapseShiftAll(double alphai, double betaxi, 
                         double betayi, double betazi);
  
  void loadExtCurv(double Kxxi, double Kxyi, double Kxzi, 
                   double Kyyi, double Kyzi, double Kzzi);
        
  double contractKllTuuSym(double Txx, double Txy, double Txz,
                           double Tyy, double Tyz, double Tzz) const;

  double getDetg()  const { return detg;  } 
  double getAlpha() const { return alpha; } 
  double getBetax() const { return betax; } 
  double getBetay() const { return betay; } 
  double getBetaz() const { return betaz; } 

  // Contract the three metric with three vectors 
  double contractAlBu(double Ax, double Ay, double Az,
      double Bx, double By, double Bz) const;

  double contractGllAuBu(double Ax, double Ay, double Az,
      double Bx, double By, double Bz) const;

  double contractGuuAlBl(double Ax, double Ay, double Az,
      double Bx, double By, double Bz) const;

  // Raise and lower three vector indices 
  void lowerAu(double Axu, double Ayu, double Azu,
      double *Axl, double *Ayl, double *Azl) const;
  
  void raiseAl(double Axl, double Ayl, double Azl,
      double *Axu, double *Ayu, double *Azu) const;
  
  void raiseSAll(double  Axxl, double  Axyl, double  Axzl, 
    double  Ayyl, double  Ayzl, double  Azzl,
    double *Axxu, double *Axyu, double *Axzu, 
    double *Ayyu, double *Ayzu, double *Azzu) const;
  
  void lowerSAuu(double  Axxu, double  Axyu, double  Axzu, 
      double  Ayyu, double  Ayzu, double  Azzu,
      double *Axxl, double *Axyl, double *Axzl, 
      double *Ayyl, double *Ayzl, double *Azzl) const;
};
}
#endif
