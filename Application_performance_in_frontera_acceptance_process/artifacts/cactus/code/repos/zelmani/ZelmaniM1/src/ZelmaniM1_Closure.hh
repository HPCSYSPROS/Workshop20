#ifndef _CLOSURECLASS_
#define _CLOSURECLASS_
#include <math.h>
#include "ZelmaniM1_Metric.hh"
using namespace std;

class Closure {

protected:	
	// Input quantities
	ThreeMetric gmet;
	double E,Fx,Fy,Fz; // E and F_i
	double vx,vy,vz;   // v^i
	double fhx,fhy,fhz; // \hat f_i
	 
	// Derived quantities
	double JJ0,JJthick,JJthin;
	double HH0,HHthick,HHthin,HHthickthick,HHthinthin,HHthickthin;
	double W,W2,W3; // Lorentz factor
	double v2,F,F2,vdotF;
        double vdotfh,Fdotfh;
	
	double JJ,tHx,tHy,tHz,nH;	 
	double xi;
	double Pxx,Pxy,Pxz,Pyy,Pyz,Pzz;
	double Fxu,Fyu,Fzu; // F^i 
	double fhxu,fhyu,fhzu; // \hat f^i 
	double vxl,vyl,vzl; // v_l 
        
	// Function to be zeroed to find the rest frame 
	// Eddington factor
	inline double zeroFunc(double xi);
	inline double zeroFunc(double xi, double *dfdxi);
  
	// The closure functions we would like to use
	inline double closureFunc(double xi){
	  	//return (3.0 + 4.0*xi*xi)/(5.0 + 2.0*sqrt(4.0 - 3.0*xi*xi));
		return 1.0/3.0 + xi*xi*(6.0 - 2.0*xi + 6.0*xi*xi)/15.0;
	}
	
	inline double closureFunc(double xi, double *dchidxi){
		// Levermore 84 closure
		//double rt = sqrt(4.0 - 3.0*xi*xi);
		//*dchidxi = 2.0*xi/rt;
		//return (3.0 + 4.0*xi)/(5.0 + 2.0*rt);
		// Maximum Entropy closure
		*dchidxi = 2.0*xi*(6.0 - 2.0*xi + 6.0*xi*xi)/15.0
		  + xi*xi*(12.0*xi - 2.0)/15.0;
		return 1.0/3.0 + xi*xi*(6.0 - 2.0*xi + 6.0*xi*xi)/15.0;
	}
		
public:
	Closure() { vx = 0.0; vy = 0.0; vz = 0.0;}
	
	
	// Some inline functions for loading the class
	// Need to load the metric before loading radiation quantities and velocity
	void loadBackground(double vxi, double vyi, double vzi, ThreeMetric gin) {
		gmet = gin;	
		loadV(vxi,vyi,vzi);
		return;
	}
			
	void loadAll(double Ei, double Fxi, double Fyi, double Fzi, 
		    double vxi, double vyi, double vzi, ThreeMetric gin) {
		loadBackground(vxi,vyi,vzi,gin);
		loadRad(Ei,Fxi,Fyi,Fzi);
		return;	
	}
	 
	void loadRad(double Ei, double Fxi,  double Fyi, double Fzi) {
		E = Ei; Fx = Fxi; Fy = Fyi; Fz = Fzi; 
		gmet.raiseAl(Fx,Fy,Fz,&Fxu,&Fyu,&Fzu);
	        F2 = gmet.contractAlBu(Fx,Fy,Fz,Fxu,Fyu,Fzu) + 1.e-50;
		F  = sqrt(F2);
	        vdotF = gmet.contractAlBu(Fx,Fy,Fz,vx,vy,vz) + 1.e-50;	
		fhx = Fx/F; fhy = Fy/F; fhz = Fz/F;
		gmet.raiseAl(fhx,fhy,fhz,&fhxu,&fhyu,&fhzu);
	        vdotfh = gmet.contractAlBu(fhx,fhy,fhz,vx,vy,vz) + 1.e-50;	
	        Fdotfh = gmet.contractAlBu(Fx,Fy,Fz,fhxu,fhyu,fhzu) + 1.e-50;
		return;
	}
	
	void loadV(double vxi,  double vyi,  double vzi) { 
		vx = vxi; vy = vyi; vz = vzi;
		gmet.lowerAu(vx,vy,vz,&vxl,&vyl,&vzl);	
		v2 = gmet.contractAlBu(vxl,vyl,vzl,vx,vy,vz);
		W = 1.0/sqrt(1.0-v2);
		W2 = W*W;
		W3 = W*W2;
		return;
	}
	
	void loadG(double gxxi, double gxyi, double gxzi,
		   double gyyi, double gyzi, double gzzi) {
		gmet.loadG(gxxi,gxyi,gxzi,gyyi,gyzi,gzzi);
	}
	
	void loadLapseShift(double alphai, double betaxi){
		gmet.loadLapseShift(alphai,betaxi);
	}
	
	void setFhat(double fhxi, double fhyi, double fhzi) {
		fhx = fhxi; fhy = fhyi; fhz = fhzi;
		double fh2 = gmet.contractGuuAlBl(fhx,fhy,fhz,fhx,fhy,fhz);
		fhx = fhx/sqrt(fh2); fhy = fhy/sqrt(fh2); fhz = fhz/sqrt(fh2);
		gmet.raiseAl(fhx,fhy,fhz,&fhxu,&fhyu,&fhzu);
		vdotfh = gmet.contractAlBu(fhx,fhy,fhz,vx,vy,vz) + 1.e-50;	
	        Fdotfh = gmet.contractAlBu(Fx,Fy,Fz,fhxu,fhyu,fhzu) + 1.e-50;
	}
	
	// Interface to the closure functions 	
	double setClosure(){ return setClosure(-1.0); }	
	double setClosure(double xi0);	
	
	void getWaveSpeeds(int direction, double *lam);
	void getWaveSpeeds(double *lam) { getWaveSpeeds(0,lam); } 
	void getAsymptoticWaveSpeeds(int direction, double *lam);
	void getAsymptoticWaveSpeeds(double *lam) { getAsymptoticWaveSpeeds(0,lam); } 
	
	void getSourceUpdate(double gE,   double gFx,  double gFy,  double gFz, 
	                     double kapa, double kaps, double eta,  double dtau, 
			     double *Ep,  double *Fxp, double *Fyp, double *Fzp, 
			     double *dJ,  double *dTau,double *dSx, double *dSy,
			     double *dSz){
	     double Jout,Np,dN;
	     getSourceUpdate(gE,0.0,gFx,gFy,gFz,kapa,kaps,eta,dtau,Ep,&Np,Fxp,Fyp,
	                     Fzp,&Jout,dJ,&dN,dTau,dSx,dSy,dSz);	
	     return;
	}
        void getSourceUpdate(double gE, double gN, double gFx, double gFy, double gFz,
                             double kapa, double kaps, double eta, double dtau, 
			     double *Ep, double *Np, double *Fxp, double *Fyp, double *Fzp, double *Jout,
			     double *dJ, double *dN, double *dTau, double *dSx, double *dSy, double *dSz);
        void getMomentumFluxes(double dlnadx, double dlnady, double dlnadz,
                               double ndW, double dWx, double dWy, double dWz,
                               double aa,  double bx,  double by,  double bz, 
                               double cxx, double cxy, double cxz, double cyy, 
                               double cyz, double czz,
			       double *fE, double *fN, double *fFx, double *fFy, double *fFz);
        
	int  getP(double *Pxxo, double *Pxyo, double *Pxzo, 
	          double *Pyyo, double *Pyzo, double *Pzzo){
		*Pxxo = Pxx; *Pxyo = Pxy; *Pxzo = Pxz;
		*Pyyo = Pyy; *Pyzo = Pyz; *Pzzo = Pzz;
		return 0;
	}
        
	int  getPll(double *Pxxo, double *Pxyo, double *Pxzo, 
	            double *Pyyo, double *Pyzo, double *Pzzo){
		
		gmet.lowerSAuu(Pxx, Pxy, Pxz, Pyy, Pyz, Pzz,
		               Pxxo,Pxyo,Pxzo,Pyyo,Pyzo,Pzzo);
		return 0;
	
	}
        
	double contractPuuTllSymmetric(double Txx, double Txy, double Txz, 
	                               double Tyy, double Tyz, double Tzz){	
		return Pxx*Txx + Pyy*Tyy + Pzz*Tzz + 2.0*Txy*Pxy 
		     + 2.0*Txz*Pxz + 2.0*Tyz*Pyz;
	}

	void getFluxes(double *fenu, double *fFx, double *fFy, double *fFz){
	
		//double sdetg = sqrt(gmet.getDetg());
		double Pxxl,Pxyl,Pxzl;
		gmet.lowerAu(Pxx,Pxy,Pxz,&Pxxl,&Pxyl,&Pxzl);
		double alpha = gmet.getAlpha();
		double betax = gmet.getBetax();
		
		*fenu = (alpha*Fxu  - betax*E );

		*fFx  = (alpha*Pxxl - betax*Fx);	
		*fFy  = (alpha*Pxyl - betax*Fy);	
		*fFz  = (alpha*Pxzl - betax*Fz);	
	
		return;
	}
	
	void getFluxes(double nuNin, double *fenu, double *fnnu, double *fFx, double *fFy, double *fFz){
	  getFluxes(0,nuNin,fenu,fnnu,fFx,fFy,fFz);
	}

	void getFluxes(int direction, double nuNin, double *fenu, double *fnnu, double *fFx, double *fFy, double *fFz){
	
		//double sdetg = sqrt(gmet.getDetg());
		double Paxl,Payl,Pazl;
		double alpha = gmet.getAlpha();
		double beta,Fau,va,tHa;
		if (direction==0) {
		  gmet.lowerAu(Pxx,Pxy,Pxz,&Paxl,&Payl,&Pazl);
		  beta = gmet.getBetax();
		  Fau  = Fxu; tHa = tHx;
		  va   = vx;
		} else if (direction==1) { 
		  gmet.lowerAu(Pxy,Pyy,Pyz,&Paxl,&Payl,&Pazl);
		  beta = gmet.getBetay();
		  Fau  = Fyu; tHa = tHy;
		  va   = vy;
		} else if (direction==2) {
		  gmet.lowerAu(Pxz,Pyz,Pzz,&Paxl,&Payl,&Pazl);
		  beta = gmet.getBetaz();
		  Fau  = Fzu; tHa = tHz;
		  va   = vz;
		} else {
		  Paxl = 0.0; Payl = 0.0; Pazl = 0.0;
		  beta = 0.0;
		  Fau  = 0.0; tHa = 0.0;
		  va   = 0.0;
		} 

		*fenu = (alpha*Fau  - beta*E );
		*fnnu = (alpha*W*va*JJ + alpha*tHa - beta*nuNin);

		*fFx  = (alpha*Paxl - beta*Fx);	
		*fFy  = (alpha*Payl - beta*Fy);	
		*fFz  = (alpha*Pazl - beta*Fz);	
	
		return;
	}
	
	
	void getAsymptoticFluxes(double nuNin, double *fenu, double *fnnu, double *fFx, double *fFy, double *fFz){
	  getAsymptoticFluxes(0,nuNin,fenu,fnnu,fFx,fFy,fFz);
	}

	void getAsymptoticFluxes(int direction, double nuNin, double *fenu, double *fnnu, double *fFx, double *fFy, double *fFz){
		
		//double sdetg = sqrt(gmet.getDetg());
		double Paxl,Payl,Pazl;
		double alpha = gmet.getAlpha();
		double beta,Fau,va,tHa;
		if (direction==0) {
		  gmet.lowerAu(Pxx,Pxy,Pxz,&Paxl,&Payl,&Pazl);
		  beta = gmet.getBetax();
		  Fau  = Fxu; tHa = tHx;
		  va   = vx;
		} else if (direction==1) { 
		  gmet.lowerAu(Pxy,Pyy,Pyz,&Paxl,&Payl,&Pazl);
		  beta = gmet.getBetay();
		  Fau  = Fyu; tHa = tHy;
		  va   = vy;
		} else if (direction==2) {
		  gmet.lowerAu(Pxz,Pyz,Pzz,&Paxl,&Payl,&Pazl);
		  beta = gmet.getBetaz();
		  Fau  = Fzu; tHa = tHz;
		  va   = vz;
		} else {
		  Paxl = 0.0; Payl = 0.0; Pazl = 0.0;
		  beta = 0.0;
		  Fau  = 0.0; tHa = 0.0;
		  va   = 0.0;
		} 
	
		double Jthick = 3.0/(2.0*W2 + 1.0) 
		               *((2.0*W2 - 1.0)*E - 2.0*W2*vdotF);
	        
		Jthick = 3.0/(4.0*W2-1.0)*E;
			
		*fenu = (alpha*4.0*W2/3.0*va*Jthick  - beta*E);
		*fnnu = (alpha*W*va*Jthick - beta*nuNin);

		*fFx  = (alpha*Paxl - beta*Fx);	
		*fFy  = (alpha*Payl - beta*Fy);	
		*fFz  = (alpha*Pazl - beta*Fz);	
	
		return;
	}

};


// Zero velocity version of the closure calculation class
class ClosureZV: public Closure {
protected:
public:
	// This is the only function we need to overide for the v^i = 0 case
	double setClosure();
	double setClosure(double xi0){ return setClosure(); }	
};

#endif
