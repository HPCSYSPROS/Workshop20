#ifndef _ASYMPTOTICFCLASS_
#define _ASYMPTOTICFCLASS_

#include <math.h>
#include "ZelmaniM1_Metric.hh"
using namespace std;

class AsymptoticF {

protected:
	
	ThreeMetric *gmetL;
	ThreeMetric *gmetR;	
	ThreeMetric *gmetC;	
	double vxR,vyR,vzR; // v^i 
	double vxL,vyL,vzL;
	double FxR,FyR,FzR; // Zone centered fluxes F_i
	double FxL,FyL,FzL;  
	double EL,ER; 

public:

	void loadAllR(double Ei, double Fxi, double Fyi, double Fzi, 
		      double vxi, double vyi, double vzi, ThreeMetric *gin) {
		gmetR = gin;	
		ER = Ei; FxR = Fxi; FyR = Fyi; FzR = Fzi; 
		vxR = vxi; vyR = vyi; vzR = vzi;
	}
	 
	void loadAllL(double Ei, double Fxi, double Fyi, double Fzi, 
		      double vxi, double vyi, double vzi, ThreeMetric *gin) {
		gmetL = gin;	
		EL = Ei; FxL = Fxi; FyL = Fyi; FzL = Fzi; 
		vxL = vxi; vyL = vyi; vzL = vzi;
	}

	void loadMetricC(ThreeMetric *gin) { gmetC = gin; } 

	double getAsymptoticFx(double kappa, double idx, double *nFx) {
	  return getAsymptoticFa(0,kappa,idx,nFx);	
	}
		
	double getAsymptoticFa(int direction, double kappa, double ida, double *nFa);	
};

double AsymptoticF::getAsymptoticFa(int direction, double kappa, double ida, double *nFa) {
	
	double vxC,vyC,vzC;
	double FxC,FyC,FzC;
	double EC;
        
	
	// Need to do a better reconstruction 
	vxC = 0.5*(vxL + vxR);
	vyC = 0.5*(vyL + vyR);
	vzC = 0.5*(vzL + vzR);
	FxC = 0.5*(FxL + FxR);
	FyC = 0.5*(FyL + FyR);
	FzC = 0.5*(FzL + FzR);
	EC  = 0.5*(EL  + ER );
        
	double v2L = gmetL->contractGllAuBu(vxL,vyL,vzL,vxL,vyL,vzL); 
	double v2R = gmetR->contractGllAuBu(vxR,vyR,vzR,vxR,vyR,vzR); 
	double v2C = gmetC->contractGllAuBu(vxC,vyC,vzC,vxC,vyC,vzC); 
	
	double vdotFL = gmetL->contractAlBu(FxL,FyL,FzL,vxL,vyL,vzL); 
	double vdotFR = gmetR->contractAlBu(FxR,FyR,FzR,vxR,vyR,vzR); 
	double vdotFC = gmetC->contractAlBu(FxC,FyC,FzC,vxC,vyC,vzC); 
	
	double WL = 1.0/sqrt(1.0 - v2L); 
	double WR = 1.0/sqrt(1.0 - v2R); 
	double WC = 1.0/sqrt(1.0 - v2C); 
	
	double JthickL = 3.0/(2.0*WL*WL + 1.0) * ((2.0*WL*WL - 1.0)*EL - 2.0*WL*WL*vdotFL);
	double JthickR = 3.0/(2.0*WR*WR + 1.0) * ((2.0*WR*WR - 1.0)*ER - 2.0*WR*WR*vdotFR);
	JthickL = 3.0/(4.0*WL*WL - 1.0)*EL;
	JthickR = 3.0/(4.0*WR*WR - 1.0)*ER;

	//double JthickC = 3.0/(2.0*WC*WC + 1.0) * ((2.0*WC*WC - 1.0)*EC - 2.0*WC*WC*vdotFC);
        
	//double JthickL = 3.0/(4.0*WL*WL-1.0)*EL;
        //double JthickR = 3.0/(4.0*WR*WR-1.0)*ER;

	double sdetgL = sqrt(gmetL->getDetg());		
	double sdetgR = sqrt(gmetR->getDetg());		
	double sdetgC = sqrt(gmetC->getDetg());		
	
	double dJda = (JthickR/sdetgR - JthickL/sdetgL)*ida*sdetgC;
	double dJdb = 0.0; // !! Think about how to calculate this, may be important !!
	double dJdc = 0.0; // !! Think about how to calculate this, may be important !!
			   // Only enters at order (v/c)^2	
	double vxl,vyl,vzl;
	gmetC->lowerAu(vxC,vyC,vzC,&vxl,&vyl,&vzl);
	//double dJdxu,dJdyu,dJdzu;
	//gmetC->raiseAl(dJdx,dJdy,dJdz,&dJdxu,&dJdyu,&dJdzu);
	
	double dJdotv, beta, vac;	
	double sdetg = sqrt(gmetC->getDetg());
	double alpha = gmetC->getAlpha();
	if (direction==0) {
	  dJdotv = gmetC->contractAlBu(vxC,vyC,vzC,dJda,dJdb,dJdc);
	  beta   = gmetC->getBetax();
	  vac    = vxC;
        } else if (direction==1) {
	  dJdotv = gmetC->contractAlBu(vxC,vyC,vzC,dJdc,dJda,dJdb);
	  beta   = gmetC->getBetay();
	  vac    = vyC;
        } else if (direction==2) {
	  dJdotv = gmetC->contractAlBu(vxC,vyC,vzC,dJdb,dJdc,dJda);
	  beta   = gmetC->getBetaz();
	  vac    = vzC;
	}
			
	//double Fxu = 4.0*WC*WC/3.0*JthickC*vxC - WC/(3.0*kappa)*(dJdx + vxC*dJdotv);	
	double Fau = - WC/(3.0*kappa)*(dJda + vac*dJdotv);	
        *nFa = -1.0/(3.0*kappa)*dJda; 
	return alpha*Fau;
}

#endif


