#include <cmath>
#include <math.h>
#include <iostream>
#include "ZelmaniM1_ClosureVec.hh"
// These need to go to make this general
//#include "cctk.h"
//#include "cctk_Parameters.h"
//#include "cctk_Arguments.h"
//#include "cctk_Functions.h"

ClosureVec::ClosureVec() {
}

ClosureVec::ClosureVec(int nsize) {
        size = nsize;
        
	// Background 
	guxx = new double[size];
        guxy = new double[size];
        guxz = new double[size];
        guyy = new double[size];
        guyz = new double[size];
        guzz = new double[size];
        
	sdetg = new double[size];
	 
        vxl  = new double[size];
        vyl  = new double[size];
        vzl  = new double[size];
        v2   = new double[size];
        W    = new double[size];
	invW2fac = new double[size];

        
	fhx  = new double[size];
        fhy  = new double[size];
        fhz  = new double[size];
        fhxu = new double[size];
        fhyu = new double[size];
        fhzu = new double[size];

	xiarr = new double[size];
        
        F    = new double[size];
	Fxu  = new double[size];
        Fyu  = new double[size];
        Fzu  = new double[size];
        Pxx  = new double[size];
        Pxy  = new double[size];
        Pxz  = new double[size];
        Pyy  = new double[size];
        Pyz  = new double[size];
        Pzz  = new double[size];
        
	JJ0     = new double[size];
	JJthick = new double[size];
	JJthin  = new double[size];
	HH0     = new double[size];
	HHthick = new double[size];
	HHthin  = new double[size];
	HHthickthick = new double[size];
	HHthinthin  = new double[size];
	HHthickthin  = new double[size];
	conv = new double[size];
	reiterate = new int [size];

#pragma ivdep
        for(int i=0; i<size; i++) {
          Fxu[i] = 0.0; Fyu[i] = 0.0; Fzu[i] = 0.0; F[i] = 0.0;
          fhx[i]  = 0.0; fhy[i]  = 0.0; fhz[i]  = 0.0;
          fhxu[i] = 0.0; fhyu[i] = 0.0; fhzu[i] = 0.0;
          vxl[i] = 0.0; vyl[i] = 0.0; vzl[i] = 0.0;
          W[i] = 1.0; v2[i] = 0.0;
	  xiarr[i] = 0.0;
          Pxx[i] = 0.0; Pxy[i] = 0.0; Pxz[i] = 0.0;
          Pyy[i] = 0.0; Pyz[i] = 0.0; Pzz[i] = 0.0;
        }
      }

ClosureVec::~ClosureVec() {
        delete [] guxx;
        delete [] guxy;
        delete [] guxz;
        delete [] guyy;
        delete [] guyz;
        delete [] guzz;
	
	delete [] fhx;
        delete [] fhy;
        delete [] fhz;

        delete [] xiarr;
        delete [] Pxx;
        delete [] Pxy;
        delete [] Pxz;
        delete [] Pyy;
        delete [] Pyz;
        delete [] Pzz;
        delete [] F;
        delete [] Fxu;
        delete [] Fyu;
        delete [] Fzu;
        delete [] fhxu;
        delete [] fhyu;
        delete [] fhzu;
        delete [] vxl;
        delete [] vyl;
        delete [] vzl;
	delete [] W;
	delete [] invW2fac;
	delete [] v2;
	delete [] JJ0;
	delete [] JJthick;
	delete [] JJthin;
	delete [] HH0;
	delete [] HHthick;
	delete [] HHthin;
	delete [] HHthickthick;
	delete [] HHthinthin;
	delete [] HHthickthin;
	delete [] conv;
	delete [] reiterate;
} 

void ClosureVec::loadBackground(double *RESTRICT vxi, double *RESTRICT vyi, double *RESTRICT vzi, double *RESTRICT ain, 
                                double *RESTRICT bx, double *RESTRICT by, double *RESTRICT bz, 
		                double *RESTRICT gxxi, double *RESTRICT gxyi, double *RESTRICT gxzi,
		                double *RESTRICT gyyi, double *RESTRICT gyzi, double *RESTRICT gzzi) {
        
	vx = vxi; vy = vyi; vz = vzi; 
        alpha = ain;
        betax = bx; betay = by; betaz = bz; 
        gxx = gxxi; gxy = gxyi; gxz = gxzi; 
        gyy = gyyi; gyz = gyzi; gzz = gzzi; 
	
#pragma ivdep
	for (int i=0;i<size;i++) {
	  double gxxL = gxx[i];
	  double gxyL = gxy[i];
	  double gxzL = gxz[i];
	  double gyyL = gyy[i];
	  double gyzL = gyz[i];
	  double gzzL = gzz[i];

	  double vxL = vx[i];
	  double vyL = vy[i];
	  double vzL = vz[i];

	  double detgL = DET3(gxxL,gxyL,gxzL,gyyL,gyzL,gzzL);
//	  double sdetgL = sqrt(DET3(gxxL,gxyL,gxzL,gyyL,gyzL,gzzL));
	  double sdetgL = sqrt(detgL);
	  double invdetg = 1.0/detgL;
//	  double invsdetgsq = 1.0/(sdetgL*sdetgL);
	  double guxxL = (gyyL*gzzL - gyzL*gyzL)*invdetg;
	  double guxyL = (gxzL*gyzL - gxyL*gzzL)*invdetg;
	  double guxzL = (gxyL*gyzL - gxzL*gyyL)*invdetg;
	  double guyyL = (gxxL*gzzL - gxzL*gxzL)*invdetg;
	  double guyzL = (gxyL*gxzL - gxxL*gyzL)*invdetg;
	  double guzzL = (gxxL*gyyL - gxyL*gxyL)*invdetg;
	  
	  double vxlL  = gxxL*vxL + gxyL*vyL + gxzL*vzL;  
	  double vylL  = gxyL*vxL + gyyL*vyL + gyzL*vzL;  
	  double vzlL  = gxzL*vxL + gyzL*vyL + gzzL*vzL;
	  
	  double v2L = vxL*vxlL + vyL*vylL + vzL*vzlL; 
	  double WL  = 1.0/sqrt(1.0-v2L);
	  double W2 = WL*WL;
	  double invW2facL = 1.0/(1.0+2.0*W2);
          
	  v2[i] = v2L;
	  W[i]  = WL;
	  invW2fac[i]  = invW2facL;

	  guxx[i] = guxxL;
          guxy[i] = guxyL;
          guxz[i] = guxzL;
          guyy[i] = guyyL;
          guyz[i] = guyzL;
          guzz[i] = guzzL;
	  sdetg[i] = sdetgL;
	  vxl[i]  = vxlL;  
	  vyl[i]  = vylL;  
	  vzl[i]  = vzlL;
	}
} 

void ClosureVec::loadRad(double *RESTRICT Ein, double *RESTRICT Fxin, double *RESTRICT Fyin, double *RESTRICT Fzin) {
	E = Ein; Fx = Fxin; Fy = Fyin; Fz = Fzin;
#pragma ivdep
	for (int i=0;i<size;i++) {
	  double guxxL = guxx[i];
	  double guxyL = guxy[i];
	  double guxzL = guxz[i];
	  double guyyL = guyy[i];
	  double guyzL = guyz[i];
	  double guzzL = guzz[i];
	  
	  double FxL = Fx[i];
	  double FyL = Fy[i];
	  double FzL = Fz[i];

	  double FxuL = guxxL*FxL + guxyL*FyL + guxzL*FzL;
	  double FyuL = guxyL*FxL + guyyL*FyL + guyzL*FzL;
	  double FzuL = guxzL*FxL + guyzL*FyL + guzzL*FzL;
            
	  double FL = sqrt(FxL*FxuL + FyL*FyuL + FzL*FzuL)+1.e-50;
	  double invF = 1.0/FL;

	  double fhxL = FxL*invF;
	  double fhyL = FyL*invF;
	  double fhzL = FzL*invF;

	  double fhxuL = guxxL*fhxL + guxyL*fhyL + guxzL*fhzL;
	  double fhyuL = guxyL*fhxL + guyyL*fhyL + guyzL*fhzL;
	  double fhzuL = guxzL*fhxL + guyzL*fhyL + guzzL*fhzL;

	  Fxu[i] = FxuL;
	  Fyu[i] = FyuL;
	  Fzu[i] = FzuL;
          F[i]   = FL;
	  fhx[i] = fhxL;
	  fhy[i] = fhyL;
	  fhz[i] = fhzL;

	  fhxu[i] = fhxuL;
	  fhyu[i] = fhyuL;
	  fhzu[i] = fhzuL;
	}
}

void ClosureVec::setClosureVecZV(double *RESTRICT xi0) {
		
      double TINY = 1.e-45;

      // Given the rest frame E and F_i, as well as the local three 
      // velocity v^i and the three metric gamma_{ij}, this routine 
      // solves the closure relation for the radiation stress tensor P^{ij}
      // in the lab frame and returns the Eddington factor and components 
      // of the stress tensor. 
      
#pragma ivdep
      for (int i=0;i<size;i++) {
	  double xi;
	  	
	  // Calculate some quantities we are going to use a lot	
          double F2 = Fx[i]*Fxu[i] + Fy[i]*Fyu[i] + Fz[i]*Fzu[i];
	  double F  = sqrt(F2);
          double EL = E[i];  

	  // The velocity is basically zero and the closure is easy to determine
	  xi = sqrt(F2/(E[i]*E[i]+TINY));
	  if (F2 >= E[i]*E[i]) xi = 1.0;
	  xiarr[i] = xi;
	  xi0[i] = xi;
	} 
	return; 
}

void ClosureVec::setClosureVec(double *RESTRICT xi0) {
		
      double TINY = 1.e-45;

      // Given the rest frame E and F_i, as well as the local three 
      // velocity v^i and the three metric gamma_{ij}, this routine 
      // solves the closure relation for the radiation stress tensor P^{ij}
      // in the lab frame and returns the Eddington factor and components 
      // of the stress tensor. 
            
      for (int i=0;i<size;i++) {
        double xi0L = xi0[i];
        double FxL = Fx[i];
        double FyL = Fy[i];
        double FzL = Fz[i];
        double FxuL = Fxu[i];
        double FyuL = Fyu[i];
        double FzuL = Fzu[i];
        double EL = E[i];
        double xiarrL;

//	if (xi0<1.e-5 || xi0>1.0) {
        if (xi0L<0.0 || xi0L>1.0) {
          double F2 = FxL*FxuL + FyL*FyuL + FzL*FzuL;
          xiarrL = sqrt(F2/(EL*EL + TINY));
          if (F2<=TINY) xiarrL = 0.0;
          xiarrL = min(xiarrL,1.0);
        } else {
          xiarrL = xi0L;
        }

        xiarr[i] = xiarrL;
      }
      return;

#pragma ivdep
      for (int i=0;i<size;i++) {
		
	// Calculate some quantities we are going to use a lot	
        
	double v2L = v2[i]; 
	double WL  = W[i];						
	double W2  = WL*WL;
	double W3  = WL*W2;
	double invW2facL = invW2fac[i];
	double vxL = vx[i];
	double vyL = vy[i];
	double vzL = vz[i];
	double FxL = Fx[i];
	double FyL = Fy[i];
	double FzL = Fz[i];
	double fhxL = fhx[i];
	double fhyL = fhy[i];
	double fhzL = fhz[i];
	double FxuL = Fxu[i];
	double FyuL = Fyu[i];
	double FzuL = Fzu[i];
        double vdotF  = vxL*FxL  + vyL*FyL  + vzL*FzL;
        double vdotfh = vxL*fhxL + vyL*fhyL + vzL*fhzL;
        double F2 = FxL*FxuL + FyL*FyuL + FzL*FzuL;
//	double F = sqrt(F2);
        double Fdotfh = FxuL*fhxL + FyuL*fhyL + FzuL*fhzL; 	 
	double EL = E[i];

	// Calculate pieces of J 	

  	double JJ0L     = W2*(EL - 2.0*vdotF); 
	double JJthinL  = W2*EL*vdotfh*vdotfh;
        double JJthickL = (W2-1.0)*invW2facL*(4.0*W2*vdotF + (3.0-2.0*W2)*EL);

	// Calculate pieces of H^alpha H_alpha 
	double cn = WL*JJ0L + WL*(vdotF - EL);
	double cv = WL*JJ0L;
	double cF =-WL;
        
	double Cthickn = WL*JJthickL;
	double Cthickv = WL*JJthickL + WL*invW2facL*((3.0-2.0*W2)*EL + (2.0*W2-1.0)*vdotF);
	double CthickF = WL*v2L;		
	
	double Cthinn  = WL*JJthinL;
	double Cthinv  = Cthinn;
	double Cthinfh = WL*EL*vdotfh;

	
	HH0[i] = cv*cv*v2L + cF*cF*F2 + 2.0*cv*cF*vdotF - cn*cn;
	HHthickthick[i] = Cthickv*Cthickv*v2L + CthickF*CthickF*F2 + 2.0*CthickF*Cthickv*vdotF - Cthickn*Cthickn;
	HHthinthin[i]   = Cthinv* Cthinv* v2L + Cthinfh*Cthinfh    + 2.0*Cthinfh* Cthinv*vdotfh- Cthinn* Cthinn ;
        HHthin[i]       = 2.0*(cv*    Cthinv* v2L + cF*    Cthinfh*Fdotfh + Cthinfh* cv*     vdotfh + Cthinv* cF*     vdotF - Cthinn* cn     );	
        HHthick[i]      = 2.0*(cv*    Cthickv*v2L + cF*    CthickF*F2 + CthickF*cv*     vdotF + Cthickv*cF*     vdotF - Cthickn*cn     );	
        HHthickthin[i]  = 2.0*(Cthinv*Cthickv*v2L + Cthinfh*CthickF*Fdotfh + Cthinfh* Cthickv*vdotfh + Cthinv* CthickF*vdotF - Cthinn* Cthickn);	
  	JJ0[i]     = W2*(EL - 2.0*vdotF); 
	JJthin[i]  = W2*EL*vdotfh*vdotfh;
        JJthick[i] = (W2-1.0)/(1.0+2.0*W2)*(4.0*W2*vdotF + (3.0-2.0*W2)*EL);
      }
      
      // Do the initial NR iteration in a vectorized way 
      int NRITERINIT = 3; 
      for (int j=0;j<=NRITERINIT;j++) {
#pragma simd
        for (int i=0;i<size;i++) {
	  // Find xi that satisfies our non-linear equation using NR iteration
	  // Compare the zero velocity limit and the optically thick limit for xi
	  // and take the larger value for the starting guess
	  double xi = xiarr[i];
	  double f,dfdxi;
          
	  double EL = E[i];
	  double EsqL = EL*EL;
	  double invEsqL = 1.0/EsqL;
	  double JJ0L = JJ0[i], JJthinL = JJthin[i], JJthickL = JJthick[i];
	  double HH0L = HH0[i], HHthinL = HHthin[i], HHthickL = HHthick[i];
	  double HHthickthickL = HHthickthick[i], HHthinthinL = HHthinthin[i], HHthickthinL = HHthickthin[i];
          
          f = zeroFunc(xi,&dfdxi,&JJ0L,&JJthinL,&JJthickL,&HH0L,&HHthickL,&HHthinL,&HHthickthickL,&HHthinthinL,&HHthickthinL,&invEsqL);
	  double fodfdxi = f/dfdxi;
          xi = xi - fodfdxi;
	  xiarr[i] = xi;
	  conv[i] = abs(fodfdxi/xi);
	}
      }       
     
      // Check if there are bad points 
      int nbad = 0;
      for (int i=0;i<size;i++) {
        if (conv[i]>1.e-3) {
	  reiterate[nbad] = i; 
	  nbad++;
	}
      }
//      cout << "nbad = " << nbad << "\n";
      
      for (int ibad=0;ibad<nbad;ibad++) {
        int i = reiterate[ibad]; 
	
	double xi = xiarr[i];
	double f,dfdxi,dftemp;
	
	double EL = E[i];
        double EsqL = EL*EL;
	double invEsqL = 1.0/EsqL;
	double JJ0L = JJ0[i], JJthinL = JJthin[i], JJthickL = JJthick[i];
	double HH0L = HH0[i], HHthinL = HHthin[i], HHthickL = HHthick[i];
	double HHthickthickL = HHthickthick[i], HHthinthinL = HHthinthin[i], HHthickthinL = HHthickthin[i];

	double xil = xi*0.95;
	double xiu = min(1.0,xi*1.05);
	double fl = zeroFunc(xil,&dftemp,&JJ0L,&JJthinL,&JJthickL,&HH0L,&HHthickL,
	                     &HHthinL,&HHthickthickL,&HHthinthinL,&HHthickthinL,&invEsqL);
        double fu = zeroFunc(xiu,&dftemp,&JJ0L,&JJthinL,&JJthickL,&HH0L,&HHthickL,
	                     &HHthinL,&HHthickthickL,&HHthinthinL,&HHthickthinL,&invEsqL);
        
	int NRITER = 7; 
	int j;
	for (j=0;j<=NRITER;j++) {	
          f = zeroFunc(xi,&dfdxi,&JJ0L,&JJthinL,&JJthickL,&HH0L,&HHthickL,
	               &HHthinL,&HHthickthickL,&HHthinthinL,&HHthickthinL,&invEsqL);
	  if (abs(f/dfdxi)<1.e-3 && abs(f)<1.e-2) break;
	  if (abs(f)<1.e-3) break;
	  if ((xi - f/dfdxi) > 1.0 || (xi - f/dfdxi) < 0.0) {
	    // Bisect if we are out of the interval
	    if (f*fl>0.0) xil = xi; 
            else xiu = xi;
	    xi = 0.5*(xiu + xil);
	  } else {
	    // Take NR step if we are within the interval
            if (f*fl>0.0) xiu = xi; 
            else xil = xi;
	    xi = xi - f/dfdxi;
	  }
        }

	if (j>=NRITER && xi >1.0-1.e-5) {xi = 1.0; j = 0; }
	
	if (j>=NRITER){ // NR Failed, so do some good old bisection 
	  xil = 0.0;
	  xiu = 1.0;
	  fl = zeroFunc(xil,&dftemp,&JJ0L,&JJthinL,&JJthickL,&HH0L,&HHthickL,
	                &HHthinL,&HHthickthickL,&HHthinthinL,&HHthickthinL,&invEsqL);
          fu = zeroFunc(xiu,&dftemp,&JJ0L,&JJthinL,&JJthickL,&HH0L,&HHthickL,
	                &HHthinL,&HHthickthickL,&HHthinthinL,&HHthickthinL,&invEsqL);
        	
	  if ((xi<0.0) || (xi>1.0)) xi = 0.5;
	    
	  if (fl*fu>0.0) { // Doesn't go to zero in interval
            // Forge ahead wildly
            if (abs(fl)<abs(fu) && abs(fl)<2.e-1) xi = xil;
            else if (abs(fu)<abs(fl) && abs(fu)<2.e-1) xi = xiu;
            else xi = 1.0;
            //else cout << "Doesn't appear to be a zero in the interval. " << fl << " " << fu << " \n";
	  } else {	
	    for (int j=0;j<100;j++) {
	      if (abs(xiu-xil)<1.e-3) break;
	      f = zeroFunc(xi,&dftemp,&JJ0L,&JJthinL,&JJthickL,&HH0L,&HHthickL,
	                   &HHthinL,&HHthickthickL,&HHthinthinL,&HHthickthinL,&invEsqL);
	      if (f*fl>0.0){ fl = f; xil = xi; }
	      else { fu = f; xiu = xi; }
	      xi = (xiu + xil)*0.5;
            }
          }
	}
        xiarr[i] = xi;
      }
      

// #pragma nofusion
//       for (int i=0;i<size;i++) {
// 	if (isnan(Pxx[i])) {
// 	  double xi = xiarr[i]; 
// 	  double chi = closureFunc(xi);
// 	  double dthin  = 1.5*chi - 0.5;
// 	  double v2 = vx[i]*vxl[i] + vy[i]*vyl[i] + vz[i]*vzl[i];
// #pragma omp critical
// 	  cerr << "Pxx is a NaN. " << dthin << " " 
// 	       << chi << " " << xi << " " << Fx[i] << " " << E[i] << " " << v2 << " " << vx[i] << " " << vxl[i] << "\n" ;
// //	  CCTK_WARN(0,"aborting!");
// 	}	
//       }

	// Return the closure factor if anyone is interested
//	return xiarr[0]; 
        for (int i=0;i<size;i++) {
          xi0[i] = xiarr[i];
        }
	return;

}

void ClosureVec::setStressTensor() {

      // Fill the stress tensor 
#pragma ivdep
      for (int i=0;i<size;i++) {
	double xi = xiarr[i]; 
	double chi = closureFunc(xi);
	double dthick = 1.5 - 1.5*chi;
	double dthin  = 1.5*chi - 0.5;
        //CCTK_VInfo(CCTK_THORNSTRING,"ClosureVec: xi = %10.2e E = %10.2e Fx = %10.2e F2 = %10.2e chi = %10.2e",
	//		xi,E,Fx,F2,chi);
	
	double v2L = v2[i]; 
        double WL  = W[i];						
	double W2  = WL*WL;
	double invW2facL = invW2fac[i];
	double vxL = vx[i];
	double vyL = vy[i];
	double vzL = vz[i];
	double EL = E[i];
	double guxxL = guxx[i];
	double guxyL = guxy[i];
	double guxzL = guxz[i];
	double guyyL = guyy[i];
	double guyzL = guyz[i];
	double guzzL = guzz[i];
	double FxL = Fx[i];
	double FyL = Fy[i];
	double FzL = Fz[i];
	double FxuL = Fxu[i];
	double FyuL = Fyu[i];
	double FzuL = Fzu[i];
	double fhxuL = fhxu[i];
	double fhyuL = fhyu[i];
	double fhzuL = fhxu[i];

        double vdotF = vxL*FxL  + vyL*FyL  + vzL*FzL;
        
	// Rest frame energy density and projected fluxes in the optically thick
	// limit.
	double Jthick = 3.0*invW2facL * ((2.0*W2 - 1.0)*EL - 2.0*W2*vdotF);
	
	double tHxt = FxuL/WL + vxL*WL*invW2facL * ((4.0*W2 + 1.0)*vdotF - 4.0*W2*EL);
	double tHyt = FyuL/WL + vyL*WL*invW2facL * ((4.0*W2 + 1.0)*vdotF - 4.0*W2*EL);
	double tHzt = FzuL/WL + vzL*WL*invW2facL * ((4.0*W2 + 1.0)*vdotF - 4.0*W2*EL);
  
	double PxxL = dthick*(Jthick/3.0*(4.0*W2*vxL*vxL + guxxL) + WL*(tHxt*vxL + vxL*tHxt)) + dthin*W2*EL*fhxuL*fhxuL;
	double PxyL = dthick*(Jthick/3.0*(4.0*W2*vxL*vyL + guxyL) + WL*(tHxt*vyL + vxL*tHyt)) + dthin*W2*EL*fhxuL*fhyuL;
	double PxzL = dthick*(Jthick/3.0*(4.0*W2*vxL*vzL + guxzL) + WL*(tHxt*vzL + vxL*tHzt)) + dthin*W2*EL*fhxuL*fhzuL;
	double PyyL = dthick*(Jthick/3.0*(4.0*W2*vyL*vyL + guyyL) + WL*(tHyt*vyL + vyL*tHyt)) + dthin*W2*EL*fhyuL*fhyuL;
	double PyzL = dthick*(Jthick/3.0*(4.0*W2*vyL*vzL + guyzL) + WL*(tHyt*vzL + vyL*tHzt)) + dthin*W2*EL*fhyuL*fhzuL;
	double PzzL = dthick*(Jthick/3.0*(4.0*W2*vzL*vzL + guzzL) + WL*(tHzt*vzL + vzL*tHzt)) + dthin*W2*EL*fhzuL*fhzuL;
        
//	PxxL = PxxL + dthin*W2*EL*fhxuL*fhxuL;
//	PxyL = PxyL + dthin*W2*EL*fhxuL*fhyuL;
//	PxzL = PxzL + dthin*W2*EL*fhxuL*fhzuL;
//	PyyL = PyyL + dthin*W2*EL*fhyuL*fhyuL;
//	PyzL = PyzL + dthin*W2*EL*fhyuL*fhzuL;
//	PzzL = PzzL + dthin*W2*EL*fhzuL*fhzuL;
        
	// Calculate the rest frame quantities
        //JJ  = JJ0 + JJthick*dthick + JJthin*dthin;
	//nH  = cn  + Cthickn*dthick + Cthinn*dthin;
	//tHx = -(cv + dthin*Cthinv + dthick*Cthickv)*vx 
	//      -(cF + dthick*CthickF)*Fxu
	//      -Cthinfh*vdotfh*fhxu;
	//tHy = -(cv + dthin*Cthinv + dthick*Cthickv)*vy
	//      -(cF + dthick*CthickF)*Fyu
	//      -Cthinfh*vdotfh*fhyu;
	//tHz = -(cv + dthin*Cthinv + dthick*Cthickv)*vz
	//      -(cF + dthick*CthickF)*Fzu
	//      -Cthinfh*vdotfh*fhzu;
        
	Pxx[i] = PxxL; 
	Pxy[i] = PxyL;
	Pxz[i] = PxzL;
	Pyy[i] = PyyL;
	Pyz[i] = PyzL;
	Pzz[i] = PzzL;
      }
      return;
}

void ClosureVec::getWaveSpeeds(int direction, double *RESTRICT lammin, double *RESTRICT lammax) {
      if (direction == 0) {
#pragma ivdep
        for (int i=0;i<size;i++) {
	  double alpL = alpha[i]; 
	  double beta = betax[i]; 
	  double gaau = guxx[i]; 
	  double Fau  = Fxu[i]; 
	  double vel  = vx[i];
	  double FLi  = 1.0/F[i];
	  double F2i  = FLi*FLi;
	  double WL   = W[i];

          //double lamthinmin = min(-beta - alpL*abs(Fau)*FLi,
	  //      		  -beta + alpL*E*Fau*F2i);
          //double lamthinmax = max(-beta + alpL*abs(Fau)*FLi,
	  //			  -beta + alpL*E*Fau*F2i);	
          // Guarantee that the speed of light is not exceeded
	  double lamthinmin = -beta - alpL*abs(Fau)*FLi;
          double lamthinmax = -beta + alpL*abs(Fau)*FLi;

          double p = alpL*vel/WL;
	  double rt = sqrt(alpL*alpL*gaau*(2.0*WL*WL + 1.0) - 2.0*WL*WL*p*p);
          double lamthickmin = min(-beta+(2.0*WL*WL*p-rt)/(2.0*WL*WL+1.0),-beta+p);	
          double lamthickmax = max(-beta+(2.0*WL*WL*p+rt)/(2.0*WL*WL+1.0),-beta+p);	
          
	  double xi = xiarr[i]; 
	  double chi = closureFunc(xi); 
	  double dthick = 1.5 - 1.5*chi;
	  double dthin  = 1.5*chi - 0.5;
	  double lamminL = dthin*lamthinmin + dthick*lamthickmin;
          double lammaxL = dthin*lamthinmax + dthick*lamthickmax;
	  
	  lammin[i] = lamminL; 
	  lammax[i] = lammaxL; 
        } 
      } else if (direction == 1) {
#pragma ivdep
        for (int i=0;i<size;i++) {
	  double alpL = alpha[i]; 
	  double beta = betay[i]; 
	  double gaau = guyy[i]; 
	  double Fau  = Fyu[i]; 
	  double vel  = vy[i];
	  double FLi  = 1.0/F[i];
	  double F2i  = FLi*FLi;
	  double WL   = W[i];

          //double lamthinmin = min(-beta - alpL*abs(Fau)*FLi,
	  //      		  -beta + alpL*E*Fau*F2i);
          //double lamthinmax = max(-beta + alpL*abs(Fau)*FLi,
	  //			  -beta + alpL*E*Fau*F2i);	
          // Guarantee that the speed of light is not exceeded
	  double lamthinmin = -beta - alpL*abs(Fau)*FLi;
          double lamthinmax = -beta + alpL*abs(Fau)*FLi;

          double p = alpL*vel/WL;
	  double rt = sqrt(alpL*alpL*gaau*(2.0*WL*WL + 1.0) - 2.0*WL*WL*p*p);
          double lamthickmin = min(-beta+(2.0*WL*WL*p-rt)/(2.0*WL*WL+1.0),-beta+p);	
          double lamthickmax = max(-beta+(2.0*WL*WL*p+rt)/(2.0*WL*WL+1.0),-beta+p);	
          
	  double xi = xiarr[i]; 
	  double chi = closureFunc(xi); 
	  double dthick = 1.5 - 1.5*chi;
	  double dthin  = 1.5*chi - 0.5;
	  double lamminL = dthin*lamthinmin + dthick*lamthickmin;
          double lammaxL = dthin*lamthinmax + dthick*lamthickmax;
	  
	  lammin[i] = lamminL; 
	  lammax[i] = lammaxL; 
        } 
	   
      } else if (direction == 2) {
#pragma ivdep
        for (int i=0;i<size;i++) {
	  double alpL = alpha[i]; 
	  double beta = betaz[i]; 
	  double gaau = guzz[i]; 
	  double Fau  = Fzu[i]; 
	  double vel  = vz[i];
	  double FLi  = 1.0/F[i];
	  double F2i  = FLi*FLi;
	  double WL   = W[i];

          //double lamthinmin = min(-beta - alpL*abs(Fau)*FLi,
	  //      		  -beta + alpL*E*Fau*F2i);
          //double lamthinmax = max(-beta + alpL*abs(Fau)*FLi,
	  //			  -beta + alpL*E*Fau*F2i);	
          // Guarantee that the speed of light is not exceeded
	  double lamthinmin = -beta - alpL*abs(Fau)*FLi;
          double lamthinmax = -beta + alpL*abs(Fau)*FLi;

          double p = alpL*vel/WL;
	  double rt = sqrt(alpL*alpL*gaau*(2.0*WL*WL + 1.0) - 2.0*WL*WL*p*p);
          double lamthickmin = min(-beta+(2.0*WL*WL*p-rt)/(2.0*WL*WL+1.0),-beta+p);	
          double lamthickmax = max(-beta+(2.0*WL*WL*p+rt)/(2.0*WL*WL+1.0),-beta+p);	
          
	  double xi = xiarr[i]; 
	  double chi = closureFunc(xi); 
	  double dthick = 1.5 - 1.5*chi;
	  double dthin  = 1.5*chi - 0.5;
	  double lamminL = dthin*lamthinmin + dthick*lamthickmin;
          double lammaxL = dthin*lamthinmax + dthick*lamthickmax;
	  
	  lammin[i] = lamminL; 
	  lammax[i] = lammaxL; 
        } 
      } else {
#pragma ivdep
        for (int i=0;i<size;i++) { 
	  lammin[i] =-1.0;
	  lammax[i] = 1.0;
	}
      }
    return;  
}

void ClosureVec::getSourceUpdate(double *RESTRICT gE, double *RESTRICT gN,
                                 double *RESTRICT gFx, double *RESTRICT gFy, double *RESTRICT gFz,
                                 double *RESTRICT kapaa, double *RESTRICT kapsa, double *RESTRICT etaa, double dtau,
                                 double *RESTRICT Ep, double *RESTRICT Np,
                                 double *RESTRICT Fxp, double *RESTRICT Fyp, double *RESTRICT Fzp,
                                 double *RESTRICT Jout, double *RESTRICT dJ, double *RESTRICT dN, double *RESTRICT dTau,
                                 double *RESTRICT dSx, double *RESTRICT dSy, double *RESTRICT dSz) {

#pragma ivdep
      for (int i=0;i<size;i++) {

        double kapa = kapaa[i];
        double kaps = kapsa[i];
        double eta  = etaa[i];
        double vxL = vx[i], vyL = vy[i], vzL = vz[i];
        double vxlL = vxl[i],  vylL = vyl[i],  vzlL = vzl[i];
        double FxL = Fx[i],  FyL = Fy[i],  FzL = Fz[i];
        double FxuL = Fxu[i],  FyuL = Fyu[i],  FzuL = Fzu[i];
        double fhxL = fhx[i],  fhyL = fhy[i],  fhzL = fhz[i];
        double fhxuL = fhxu[i],  fhyuL = fhyu[i],  fhzuL = fhzu[i];
        double gEL = gE[i], gNL = gN[i], gFxL = gFx[i],  gFyL = gFy[i],  gFzL = gFz[i];
	double xiarrL = xiarr[i];
        double WL  = W[i];
	double invW2facL = invW2fac[i];

        double v2L = v2[i];
//        double W  = 1.0/sqrt(1.0-v2);
        double W2 = WL*WL;
        double W3 = WL*W2;
        double vdotF  = vxL*FxL  + vyL*FyL  + vzL*FzL;
        double vdotfh = vxL*fhxL + vyL*fhyL + vzL*fhzL;
//        double F2 = FxL*FxuL + FyL*FyuL + FzL*FzuL;
//        double F = sqrt(F2);
        double Fdotfh = FxuL*fhxL + FyuL*fhyL + FzuL*fhzL;

        // Given the explicit parts of the update (gE,gFx,etc.), calculate the implicit update to E and F_i
        // Assumes an explicit closure factor           
        double chi = closureFunc(xiarrL);
        double dthin  = 1.5*chi - 0.5;
        double dthick = 1.5 - 1.5*chi;

        // Build the components of the rest frame energy density
        double JdE = W2 + dthin*W2*vdotfh*vdotfh + dthick*(3.0-2.0*W2)*(W2-1.0)*invW2facL;
        double JdF = (4.0*dthick*W2*(W2-1.0)*invW2facL - 2.0*W2); // Derivative of wrt to any F_i component is JdF*v^i       

        // Build the components of the E source term
        double SE0  =  dtau*WL*eta;
        double SEdE =  dtau*WL*kaps*(1.0 - JdE); // Diagonal portion that goes to zero with v
        double SEdF = -dtau*WL*(kapa + kaps*(1.0 + JdF)); // Derivative of source term wrt to any F_i component is SEdF*v^i
        double dE   = 1.0 + dtau*WL*kapa; // Diagonal term that is finite in zero velocity limit

        // Build the components of the F_i source terms
        // To get any SF_idx take SFdx*v_i
        double SF0  =  dtau*WL*eta;
        double SFdE = -dtau*(WL*kaps*JdE + WL*invW2facL*dthick*(kaps+kapa)*(3.0-2.0*W2));
        double SFdF = -dtau*(WL*kaps*JdF + WL*invW2facL*dthick*(kaps+kapa)*(2.0*W2-1.0)); // Derivative of source term j wrt to 

          // any F_i component is SFdF*v^i*v_j
        //double dF    = 1.0 + dtau*W*(kapa+kaps)*(1.0 - dthick*v2 - dthin*E/sqrt(F2)*vdotfh);
        double dF = 1.0 + dtau*WL*(kapa+kaps)*(1.0 - dthick*v2L - dthin*vdotfh); // Assume E = sqrt(F) when last term matters

        // Build the RHS
        double bE = gEL  + SE0;
        double bx = gFxL + SF0*vxlL;
        double by = gFyL + SF0*vylL;
        double bz = gFzL + SF0*vzlL;

        // Solve Matrix equation directly
        // This is relatively simple because of the high degree of symmetry of the linear equations 
        // All of this is Mathematica generated
        double detA = dF*dF*( (dF + SFdF*v2L)*(dE + SEdE) - SEdF*SFdE*v2L );

        double bdotv = bx*vxL + by*vyL + bz*vzL;
	double invdE = 1.0/dE;
        double aE = -dF*dF*invdE*( bdotv*dE*SEdF + bE   *(dF*SEdE - (SEdF*SFdE - SEdE*SFdF)*v2L) );
        double aF = -dF*      ( bE   *dF*SFdE + bdotv*(dE*SFdF - (SEdF*SFdE - SEdE*SFdF)   ) );

	double invdetA = 1.0/detA;
        double EpL  = bE*invdE +      aE*invdetA; // The first term is just the diagonal in the zero velocity limit
	double invdF = 1.0/dF;
        double FxpL = bx*invdF + vxlL*aF*invdetA;
        double FypL = by*invdF + vylL*aF*invdetA;
        double FzpL = bz*invdF + vzlL*aF*invdetA;

        // Calculate the change in energy in the rest frame
        double Fdotv = (vxL*FxpL + vyL*FypL + vzL*FzpL);
        double J = JdE*EpL + JdF*Fdotv;
        double JoutL = J;
        double dJL = dtau*(eta-kapa*J);
        double dTauL = -(SEdE*EpL + SEdF*Fdotv - SE0)      - (dE - 1.0)*EpL ;
        double dSxL  = -(SFdE*EpL + SFdF*Fdotv - SF0)*vxlL - (dF - 1.0)*FxpL; //Lowered index
        double dSyL  = -(SFdE*EpL + SFdF*Fdotv - SF0)*vylL - (dF - 1.0)*FypL;
        double dSzL  = -(SFdE*EpL + SFdF*Fdotv - SF0)*vzlL - (dF - 1.0)*FzpL;

        // Calculate things required for number source term     
        double Jf      = W2*(EpL - 2.0*Fdotv);
        double Jfthin  = W2*EpL*vdotfh*vdotfh;
        double Jfthick = (W2-1.0)*invW2facL*(4.0*W2*Fdotv + (3.0-2.0*W2)*EpL);
        double cn = WL*Jf + WL*(Fdotv - EpL);
        double Cthickn = WL*Jfthick;
        double Cthinn  = WL*Jfthin;

        double nHf = cn + Cthickn*dthick + Cthinn*dthin;
        double NpL = (gNL + dtau*eta + dtau*kapa*nHf)/(1.0 + dtau*WL*kapa);
        double dN  = dtau*(eta - kapa*WL*NpL + kapa*nHf);

        Ep[i]  = EpL;
        Np[i]  = NpL;
        Fxp[i] = FxpL;
        Fyp[i] = FypL;
        Fzp[i] = FzpL;

        Jout[i] = JoutL;
        dJ[i]   = dJL;
        dTau[i] = dTauL;
        dSx[i]  = dSxL;
        dSy[i]  = dSyL;
        dSz[i]  = dSzL;

      }
      return;

}


//void ClosureVec::getAsymptoticWaveSpeeds(int direction, double * lam) {
//	
//	double alpha = gmet.getAlpha();
//	double beta;
//	double vel;
//	double gaau;
//	if (direction==0) {
//	  beta = gmet.getBetax();
//	  vel  = vx;
//	  gaau = gmet.guxx;
//	} else if (direction==1) {
//	  beta = gmet.getBetay();
//	  vel  = vy;
//	  gaau = gmet.guyy;
//        } else if (direction==2) {
//	  beta = gmet.getBetaz();
//	  vel  = vz;
//	  gaau = gmet.guzz;
//	} else {
//	  beta = 0.0;
//	  vel  = 0.0;
//	  gaau = 1.0;
//	} 
//
//	lam[0] = -beta + 4.0*alpha*W2*vel/(2.0*W2+1.0);
//	         //-2.0*alpha/(2.0*W2+1.0)*sqrt(2.0*W2*(W2-1.0)*(2.0*W2-3.0+2.0*gaau*pow(W2-1.0,2.0)));
//	lam[1] = -beta + 4.0*alpha*W2*vel/(2.0*W2+1.0);
//	         //+ 2.0*alpha/(2.0*W2+1.0)*sqrt(2.0*W2*(W2-1.0)*(2.0*W2-3.0+2.0*gaau*pow(W2-1.0,2.0)));
//	lam[2] = lam[0];
//	lam[3] = lam[1];
//}

//void ClosureVec::getWaveSpeeds(int direction,double * lam) {
// 	//double TINY = 1.e-45;
//						
//	double alpha = gmet.getAlpha();
//	double beta,vel,gaau,Fau;
//	if (direction==0) {
//	  beta = gmet.getBetax();
//	  vel  = vx;
//	  gaau = gmet.guxx;
//	  Fau  = Fxu;
//	} else if (direction==1) {
//	  beta = gmet.getBetay();
//	  vel  = vy;
//	  gaau = gmet.guyy;
//	  Fau  = Fyu;
//        } else if (direction==2) {
//	  beta = gmet.getBetaz();
//	  vel  = vz;
//	  gaau = gmet.guzz;
//	  Fau  = Fzu;
//	} else {
//	  beta = 0.0;
//	  vel  = 0.0;
//	  gaau = 1.0;
//	  Fau  = 0.0;
//	} 
//
//        double lamthinmin = min(-beta - alpha*abs(Fau)/F,
//				-beta + alpha*E*Fau/F2);
//        double lamthinmax = max(-beta + alpha*abs(Fau)/F,
//				-beta + alpha*E*Fau/F2);	
//        // Guarantee that the speed of light is not exceeded
//	lamthinmin = -beta - alpha*abs(Fau)/F;
//        lamthinmax = -beta + alpha*abs(Fau)/F;
//
//        double p = alpha*vel/W;
//	double rt = sqrt(alpha*alpha*gaau*(2.0*W*W + 1.0) - 2.0*W*W*p*p);
//        double lamthickmin = min(-beta+(2.0*W*W*p-rt)/(2.0*W*W+1.0),-beta+p);	
//        double lamthickmax = max(-beta+(2.0*W*W*p+rt)/(2.0*W*W+1.0),-beta+p);	
//
//	double chi = closureFunc(xi); 
//	double dthick = 1.5 - 1.5*chi;
//	double dthin  = 1.5*chi - 0.5;
//	
//	// Interpolate the wave speeds as in Shibata 
//        lam[0] = dthin*lamthinmin + dthick*lamthickmin;	
//        lam[1] = dthin*lamthinmax + dthick*lamthickmax;
//        lam[2] = lam[0];
//        lam[3] = lam[1];
//	
//	return;
//
//}

inline double ClosureVec::zeroFunc(double xi,double *dfdxi,double *JJ0,double *JJthin, double *JJthick,
                                   double *HH0,double *HHthick,double *HHthin,double *HHthickthick,
				   double *HHthinthin, double *HHthickthin, double *invEscale2){	
	double dchidxi;
  	double chi = closureFunc(xi,&dchidxi);
	double athin  = 1.5*chi - 0.5;
	double athick = 1.5 - 1.5*chi;
	double dathindxi  =  1.5*dchidxi;
	double dathickdxi = -1.5*dchidxi; 
		
	double J = (*JJ0) + athin*(*JJthin) + athick*(*JJthick);
	double dJdxi = (*JJthin)*dathindxi + (*JJthick)*dathickdxi;

	double g = (*HH0) + athick*(*HHthick) + athin*(*HHthin) + athick*athin*(*HHthickthin) 
	          + athick*athick*(*HHthickthick) + athin*athin*(*HHthinthin);
	double dgdathin  = (*HHthin)  + (*HHthickthin)*athick + 2.0*(*HHthinthin)*athin;
	double dgdathick = (*HHthick) + (*HHthickthin)*athin  + 2.0*(*HHthickthick)*athick;
	
	*dfdxi = 2.0*J*dJdxi*xi*xi + 2.0*J*J*xi - dathindxi*dgdathin - dathickdxi*dgdathick;
	*dfdxi = *dfdxi*(*invEscale2);
	return (J*J*xi*xi - g)*(*invEscale2); 
}


//inline double ClosureVec::getEquations() {
//	double fE,fFx,fFy,fFz,fClose;
//	
//	return;
//}

