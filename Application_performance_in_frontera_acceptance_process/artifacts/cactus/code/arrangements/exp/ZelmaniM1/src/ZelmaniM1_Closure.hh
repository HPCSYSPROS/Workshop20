#ifndef _CLOSURECLASS_
#define _CLOSURECLASS_
#include <math.h>
#include <vector>
#include "ZelmaniM1_Metric.hh"
using namespace std;

class Closure {
public:
  Closure(double vx=0.0, double vy=0.0, double vz=0.0, 
      ThreeTensor::Metric gamma=ThreeTensor::Metric(),
      double lambda=1.0, double fhxi=0.0, double fhyi=0.0, double fhzi=0.0);

  double setClosure(double E, double Fx, double Fy, double Fz, 
      double xi0=-1.0, bool Force_closure=false);  
  
  std::vector<double> getWaveSpeeds(int direction, bool fixed_closure=false) const;
  std::vector<double> getAsymptoticWaveSpeeds(int direction) const;
 
  // Return the fluxes of radiation quantities 
  void getFluxes(int direction, double *fenu, 
      double *fFx, double *fFy, double *fFz) const {
    double fnnu;
    getFluxes(direction,0.0,fenu,&fnnu,fFx,fFy,fFz);
  }
  void getFluxes(int direction, double nuNin, double *fenu, double *fnnu, 
      double *fFx, double *fFy, double *fFz) const;
  void getAsymptoticAdvectiveFluxes(int direction, double *fenu, 
      double *fFx, double *fFy, double *fFz) const;
  friend double getAsymptoticDiffusiveFlux(double kappa, double ida, 
      const Closure& cL, const Closure& Cr, const ThreeTensor::Metric& gmetC);
 
  // Return the momentum space fluxes 
  void getMomentumFluxes(double dlnadx, double dlnady, double dlnadz,
      double ndW, double dWx, double dWy, double dWz,
      double aa,  double bx,  double by,  double bz, 
      double cxx, double cxy, double cxz, double cyy, 
      double cyz, double czz, double *fE, double *fN, 
      double *fFx, double *fFy, double *fFz) const;
        
  void getMomentumFluxes(double dlnadx, double dlnady, double dlnadz,
      double *fE, double *fN, double *fFx, double *fFy, double *fFz) const;
        
  // Invert the source terms given explicit fluxes 
  void getSourceUpdate(double gE,   double gFx,  double gFy,  double gFz, 
      double kapa, double kaps, double eta,  double dtau, 
      double *Ep,  double *Fxp, double *Fyp, double *Fzp, 
      double *dJ,  double *dTau,double *dSx, double *dSy,
      double *dSz) const {
    double Jout,Np,dN;
    getSourceUpdate(gE,0.0,gFx,gFy,gFz,kapa,kaps,eta,dtau,Ep,&Np,Fxp,Fyp,
        Fzp,&Jout,dJ,&dN,dTau,dSx,dSy,dSz);  
  }
  void getSourceUpdate(double gE, double gN, double gFx, double gFy, double gFz,
      double kapa, double kaps, double eta, double dtau, 
      double *Ep, double *Np, double *Fxp, double *Fyp, 
      double *Fzp, double *Jout, double *dJ, double *dN, 
      double *dTau, double *dSx, double *dSy, double *dSz) const;

  // Functions for dealing with pressure tensor 
  void getP(double *Pxxo, double *Pxyo, double *Pxzo, 
      double *Pyyo, double *Pyzo, double *Pzzo) const {
    *Pxxo = Pxx; *Pxyo = Pxy; *Pxzo = Pxz;
    *Pyyo = Pyy; *Pyzo = Pyz; *Pzzo = Pzz;
  }
        
  void getPll(double *Pxxo, double *Pxyo, double *Pxzo, 
      double *Pyyo, double *Pyzo, double *Pzzo) const {
    gmet.lowerSAuu(Pxx, Pxy, Pxz, Pyy, Pyz, Pzz,
                   Pxxo,Pxyo,Pxzo,Pyyo,Pyzo,Pzzo);
  }
        
  double contractPuuTllSymmetric(double Txx, double Txy, double Txz, 
      double Tyy, double Tyz, double Tzz) const { 
    return Pxx*Txx + Pyy*Tyy + Pzz*Tzz + 2.0*Txy*Pxy 
         + 2.0*Txz*Pxz + 2.0*Tyz*Pyz;
  }
  
  // Debug string so one can determine where errors are coming from 
  void setDebugInfo(std::string debinfo) { debout = debinfo; }  

protected:  
  const double TINY;
  
  // Quantities describing the background (always available)
  double vx,vy,vz;    // v^i
  double vxl,vyl,vzl; // v_l 
  ThreeTensor::Metric gmet;
  double v2;
  double W,W2,W3; // Lorentz factor
  
  // Quantities describing radiation field
  
  double lambda;         // Interpolation for pressure tensor, from constructor
  double fhxi,fhyi,fhzi; // Externally specified unit vector, from constructor
  
  double E,Fx,Fy,Fz;     // E and F_i
  double Fxu,Fyu,Fzu;    // F^i 
  double fhx,fhy,fhz;    // \hat f_i
  double fhxu,fhyu,fhzu; // \hat f^i 
  double F,F2,vdotF;
  double vdotfh,Fdotfh;
  
  double xi; // Calculated xi
  double JJ,tHx,tHy,tHz,nH;  
  double Pxx,Pxy,Pxz,Pyy,Pyz,Pzz;
  
  // String for debugging purposes
  std::string debout;
  
  // Stuff for finding the closure
  double JJ0;
  double JJthin;
  double JJthick; 
  double HH0;
  double HHthickthick;
  double HHthinthin; 
  double HHthin;
  double HHthick;
  double HHthickthin;
  double Eclosure;
  double zFunc(double xi); 
  double zFuncWithD(double xi, double *dfdxi);
  
  // The actual analytic closure functions
  double closureFunc(double xi) const;
  double closureFunc(double xi, double *dchidxi) const;
     
};


#endif
