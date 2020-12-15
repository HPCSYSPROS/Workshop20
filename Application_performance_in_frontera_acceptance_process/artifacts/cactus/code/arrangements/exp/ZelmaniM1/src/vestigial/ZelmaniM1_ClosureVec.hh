#ifndef _CLOSUREVECCLASS_
#define _CLOSUREVECCLASS_
#include <math.h>
#include "ZelmaniM1_Metric.hh"

#define RESTRICT __restrict

#define CONTRACTAUBU(Ax,Ay,Az,Bx,By,Bz,i) (  gxx[i]*Ax*Bx + gxy[i]*Ax*By + gxz[i]*Ax*Bz \
                                           + gxy[i]*Ay*Bx + gyy[i]*Ay*By + gyz[i]*Ay*Bz \
                                           + gxz[i]*Az*Bx + gyz[i]*Az*By + gzz[i]*Az*Bz)
#define CONTRACTALBL(Ax,Ay,Az,Bx,By,Bz,i) (  gxxu[i]*Ax*Bx + gxyu[i]*Ax*By + gxzu[i]*Ax*Bz \
                                           + gxyu[i]*Ay*Bx + gyyu[i]*Ay*By + gyzu[i]*Ay*Bz \
                                           + gxzu[i]*Az*Bx + gyzu[i]*Az*By + gzzu[i]*Az*Bz)
#define DET3(axx,axy,axz,ayy,ayz,azz) (axx*ayy*azz + 2.0*axy*axz*ayz - ayy*axz*axz - axx*ayz*ayz - azz*axy*axy)

using namespace std;

class ClosureVec {

protected:	
        int size; 

	// Input quantities
	double *RESTRICT E,  *RESTRICT Fx, *RESTRICT Fy, *RESTRICT Fz; // E and F_i
	double *RESTRICT vx, *RESTRICT vy, *RESTRICT vz;               // v^i
	double *RESTRICT fhx,*RESTRICT fhy,*RESTRICT fhz;              // \hat f_i
        double *RESTRICT F;

        // Metric quantities
	double *RESTRICT alpha, *RESTRICT betax, *RESTRICT betay, *RESTRICT betaz; 
	double *RESTRICT gxx,   *RESTRICT gxy,   *RESTRICT gxz;
	double *RESTRICT gyy,   *RESTRICT gyz,   *RESTRICT gzz;
	double *RESTRICT guxx,  *RESTRICT guxy,  *RESTRICT guxz;
	double *RESTRICT guyy,  *RESTRICT guyz,  *RESTRICT guzz;
        double *RESTRICT sdetg;
			 
        double *RESTRICT JJ0, *RESTRICT JJthin,	*RESTRICT JJthick;
        double *RESTRICT HH0, *RESTRICT HHthin,	*RESTRICT HHthick;
        double *RESTRICT HHthickthick, *RESTRICT HHthinthin, *RESTRICT HHthickthin;
	
	//double *RESTRICT JJ,*RESTRICT tHx,*RESTRICT tHy,*RESTRICT tHz,*RESTRICT nH;	 
	double *RESTRICT xiarr;
	double *RESTRICT Pxx,*RESTRICT Pxy,*RESTRICT Pxz,*RESTRICT Pyy,*RESTRICT Pyz,*RESTRICT Pzz;
	double *RESTRICT Fxu,*RESTRICT Fyu,*RESTRICT Fzu; // F^ i
	double *RESTRICT fhxu,*RESTRICT fhyu,*RESTRICT fhzu; // \hat f^i 
	double *RESTRICT vxl,*RESTRICT vyl,*RESTRICT vzl; // v_l 
        
	double *RESTRICT conv;
	int *RESTRICT reiterate;
	double *RESTRICT v2, *RESTRICT W, *RESTRICT invW2fac; 
	// Function to be zeroed to find the rest frame 
	// Eddington factor
        inline double zeroFunc(double xi,double *dfdxi,double *JJ0,double *JJthin, double *JJthick,
                               double *HH0,double *HHthick,double *HHthin,double *HHthickthick,
			       double *HHthinthin, double *HHthickthin, double *Escale);
  
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
	ClosureVec();
	ClosureVec(int nsize);
	~ClosureVec();
	
        void loadBackground(double *vxi, double *vyi, double *vzi, double *ain, 
                            double *bx, double *by, double *bz, 
		            double *gxxi, double *gxyi, double *gxzi,
		            double *gyyi, double *gyzi, double *gzzi);
        void loadRad(double *Ein, double *Fxin, double *Fyin, double *Fzin);
	
	// Interface to the closure functions 	
//	double setClosureVec(){ return setClosureVec(-1.0); }	
	void setClosureVec(double *RESTRICT xi0);	
	void setClosureVecZV(double *RESTRICT xi0);	
	
	void setStressTensor();
        void getSourceUpdate(double *RESTRICT gE,    double *RESTRICT gN,
                             double *RESTRICT gFx,   double *RESTRICT gFy,   double *RESTRICT gFz,
                             double *RESTRICT kapaa, double *RESTRICT kapsa, double *RESTRICT etaa, double dtau,
                             double *RESTRICT Ep,    double *RESTRICT Np,
                             double *RESTRICT Fxp,   double *RESTRICT Fyp,   double *RESTRICT Fzp,
                             double *RESTRICT Jout,  double *RESTRICT dJ,    double *RESTRICT dN,   double *RESTRICT dTau,
                             double *RESTRICT dSx,   double *RESTRICT dSy,   double *RESTRICT dSz);
	//void getWaveSpeeds(int direction, double *lam);
	//void getWaveSpeeds(double *lam) { getWaveSpeeds(0,lam); } 
        void getWaveSpeeds(int direction, double *RESTRICT lammin, double *RESTRICT lammax);
	//void getAsymptoticWaveSpeeds(int direction, double *lam);
	//void getAsymptoticWaveSpeeds(double *lam) { getAsymptoticWaveSpeeds(0,lam); } 
	
	void getFluxes(int direction, double *RESTRICT nuNin, double *RESTRICT fenu, double *RESTRICT fnnu, double *RESTRICT fFx, double *RESTRICT fFy, double *RESTRICT fFz){
	
                if (direction==0) {
#pragma ivdep
		  for (int i=0;i<size;i++){
		    double alphaL = alpha[i];
		    double betaxL = betax[i];
		    double FxuL = Fxu[i];
		    double EL = E[i];
		    double gxxL = gxx[i];
		    double gxyL = gxy[i];
		    double gxzL = gxz[i];
		    double gyyL = gyy[i];
		    double gyzL = gyz[i];
		    double gzzL = gzz[i];
		    double PxxL = Pxx[i];
		    double PxyL = Pxy[i];
		    double PxzL = Pxz[i];
		    double FxL = Fx[i];
		    double FyL = Fy[i];
		    double FzL = Fz[i];

		    double fenuL = (alphaL*FxuL  - betaxL*EL);
		    double fnnuL = 0.0; //(alpha[i]*W*vx[i]*JJ + alpha*tHa - beta*nuNin);
                    double Pxxl = gxxL*PxxL + gxyL*PxyL + gxzL*PxzL; 
                    double Pxyl = gxyL*PxxL + gyyL*PxyL + gyzL*PxzL; 
                    double Pxzl = gxzL*PxxL + gyzL*PxyL + gzzL*PxzL; 

		    double fFxL  = (alphaL*Pxxl - betaxL*FxL);	
		    double fFyL  = (alphaL*Pxyl - betaxL*FyL);	
		    double fFzL  = (alphaL*Pxzl - betaxL*FzL);	

		    fenu[i] = fenuL;
		    fnnu[i] = fnnuL; //(alpha[i]*W*vx[i]*JJ + alpha*tHa - beta*nuNin);
		    fFx[i]  = fFxL;	
		    fFy[i]  = fFyL;	
		    fFz[i]  = fFzL;	
	          }
		} else if (direction==1) {
#pragma ivdep
		  for (int i=0;i<size;i++){
		    double alphaL = alpha[i];
		    double betayL = betay[i];
		    double FyuL = Fyu[i];
		    double EL = E[i];
		    double gxxL = gxx[i];
		    double gxyL = gxy[i];
		    double gxzL = gxz[i];
		    double gyyL = gyy[i];
		    double gyzL = gyz[i];
		    double gzzL = gzz[i];
		    double PxyL = Pxy[i];
		    double PyyL = Pyy[i];
		    double PyzL = Pyz[i];
		    double FxL = Fx[i];
		    double FyL = Fy[i];
		    double FzL = Fz[i];

		    double fenuL = (alphaL*FyuL  - betayL*EL);
		    double fnnuL = 0.0; //(alphaL*W*vxL*JJ + alpha*tHa - beta*nuNin);

                    double Pyxl = gxxL*PxyL + gxyL*PyyL + gxzL*PyzL; 
                    double Pyyl = gxyL*PxyL + gyyL*PyyL + gyzL*PyzL; 
                    double Pyzl = gxzL*PxyL + gyzL*PyyL + gzzL*PyzL; 

		    double fFxL  = (alphaL*Pyxl - betayL*FxL);	
		    double fFyL  = (alphaL*Pyyl - betayL*FyL);	
		    double fFzL  = (alphaL*Pyzl - betayL*FzL);	/**/

		    fenu[i] = fenuL;
		    fnnu[i] = fnnuL; //(alpha[i]*W*vx[i]*JJ + alpha*tHa - beta*nuNin);
		    fFx[i]  = fFxL;	
		    fFy[i]  = fFyL;	
		    fFz[i]  = fFzL;	
	          }
		} else if (direction==2) {
#pragma ivdep
		  for (int i=0;i<size;i++){
		    double alphaL = alpha[i];
		    double betazL = betaz[i];
		    double FzuL = Fzu[i];
		    double EL = E[i];
		    double gxxL = gxx[i];
		    double gxyL = gxy[i];
		    double gxzL = gxz[i];
		    double gyyL = gyy[i];
		    double gyzL = gyz[i];
		    double gzzL = gzz[i];
		    double PxzL = Pxz[i];
		    double PyzL = Pyz[i];
		    double PzzL = Pzz[i];
		    double FxL = Fx[i];
		    double FyL = Fy[i];
		    double FzL = Fz[i];

		    double fenuL = (alphaL*FzuL  - betazL*EL);
		    double fnnuL = 0.0; //(alphaL*W*vxL*JJ + alpha*tHa - beta*nuNin);
                    double Pzxl = gxxL*PxzL + gxyL*PyzL + gxzL*PzzL; 
                    double Pzyl = gxyL*PxzL + gyyL*PyzL + gyzL*PzzL; 
                    double Pzzl = gxzL*PxzL + gyzL*PyzL + gzzL*PzzL; 
		    double fFxL  = (alphaL*Pzxl - betazL*FxL);	
		    double fFyL  = (alphaL*Pzyl - betazL*FyL);	
		    double fFzL  = (alphaL*Pzzl - betazL*FzL);	

		    fenu[i] = fenuL;
		    fnnu[i] = fnnuL; //(alpha[i]*W*vx[i]*JJ + alpha*tHa - beta*nuNin);
		    fFx[i]  = fFxL;	
		    fFy[i]  = fFyL;	
		    fFz[i]  = fFzL;	
	          }
		}
		return;
	}
	
};

#endif
