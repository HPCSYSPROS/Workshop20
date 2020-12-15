#include <cmath>
#include <math.h>
#include <iostream>
#include "ZelmaniM1_Closure.hh"
// These need to go to make this general
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

double Closure::setClosure(double xi0) {
		
	// Given the rest frame E and F_i, as well as the local three 
	// velocity v^i and the three metric gamma_{ij}, this routine 
	// solves the closure relation for the radiation stress tensor P^{ij}
	// in the lab frame and returns the Eddington factor and components 
	// of the stress tensor. 
	
	// Calculate some quantities we are going to use a lot	
	double TINY = 1.e-45;
						
	// Calculate pieces of J 	
  	JJ0 = W2*(E - 2.0*vdotF); 
	JJthin  = W2*E*vdotfh*vdotfh;
        JJthick = (W2-1.0)/(1.0+2.0*W2)*(4.0*W2*vdotF + (3.0-2.0*W2)*E);

	// Calculate pieces of H^alpha H_alpha 
	double cn = W*JJ0 + W*(vdotF - E);
	double cv = W*JJ0;
	double cF =-W;
        
	double Cthickn = W*JJthick;
	double Cthickv = W*JJthick + W/(2.0*W2+1.0)*((3.0-2.0*W2)*E + (2.0*W2-1.0)*vdotF);
	double CthickF = W*v2;		
	
	double Cthinn  = W*JJthin;
	double Cthinv  = Cthinn;
	double Cthinfh = W*E*vdotfh;

	if (v2>1.e-15) { // Use the finite velocity closure 
	
	  HH0 = cv*cv*v2 + cF*cF*F2 + 2.0*cv*cF*vdotF - cn*cn;
	  HHthickthick = Cthickv*Cthickv*v2 + CthickF*CthickF*F2 + 2.0*CthickF*Cthickv*vdotF - Cthickn*Cthickn;
	  HHthinthin   = Cthinv* Cthinv* v2 + Cthinfh*Cthinfh    + 2.0*Cthinfh* Cthinv*vdotfh- Cthinn* Cthinn ;
          HHthin       = 2.0*(cv*    Cthinv* v2 + cF*    Cthinfh*Fdotfh + Cthinfh* cv*     vdotfh + Cthinv* cF*     vdotF - Cthinn* cn     );	
          HHthick      = 2.0*(cv*    Cthickv*v2 + cF*    CthickF*F2 + CthickF*cv*     vdotF + Cthickv*cF*     vdotF - Cthickn*cn     );	
          HHthickthin  = 2.0*(Cthinv*Cthickv*v2 + Cthinfh*CthickF*Fdotfh + Cthinfh* Cthickv*vdotfh + Cthinv* CthickF*vdotF - Cthinn* Cthickn);	

	  // Find xi that satisfies our non-linear equation using NR iteration
	  // Compare the zero velocity limit and the optically thick limit for xi
	  // and take the larger value for the starting guess
	  double HHt = HH0 + HHthick + HHthickthick; 
	  double JJt = JJ0 + JJthick;
	  if (xi0<0.0 || xi0>1.0) {
	  	xi = sqrt(F2/(E*E + TINY));
          	if (F2<=TINY) xi = 0.0;
		xi = min(xi,1.0);
	  } else if (xi0<1.e-5) { 
		HHt = max(0.0,HHt);
		xi = sqrt(HHt/(JJt*JJt)); 
	  } else {
	  	xi = xi0;
	  }
          
	  double dfdxi = 1.0;
	  double f     = 1.0;
	  int i=0;
	  int NRITER = 10;
          
	  double xil = xi*0.95;
	  double xiu = min(1.0,xi*1.05);
          double fl = zeroFunc(xil);
	  double fu = zeroFunc(xiu); 

	  //if (fl*fu>0.0) cout << "This interval must have either multiple or no roots " << fl << " " << fu << "\n.";
	  for (i=0;i<=NRITER;i++) {
	  	f = zeroFunc(xi,&dfdxi);

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
          
          
	  if (i>=NRITER && xi >1.0-1.e-5) {xi = 1.0; i = 0; }
	  if (HHt<0.0) { 
	    xi = 0.0; 
	  } else if (i>=NRITER){ // NR Failed, so do some good old bisection 
          	
	  	//CCTK_VInfo(CCTK_THORNSTRING,"Bisecting: fNR = %10.2e xiNR = %10.2e xi = %10.2e xi0 = %10.2e E2 = %10.2e F2 = %10.2e",
	  	//	f,xi,sqrt(F2/(E*E+TINY)),xi0,E*E,F2); 

	  	fl = zeroFunc(0.0);
	  	fu = zeroFunc(1.0);
	  	xil = 0.0;
	  	xiu = 1.0;
	  	
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
	  			f = zeroFunc(xi);
	  			if (f*fl>0.0){ fl = f; xil = xi; }
	  			else { fu = f; xiu = xi; }
	  			xi = (xiu + xil)*0.5;
	  		}
	  	}
	  }
	} else {
	  // The velocity is basically zero and the closure is easy to determine
	  xi = sqrt(F2/(E*E+TINY));
	  if (F2 >= E*E) xi = 1.0;
	}
        
	// Fill the stress tensor 
	double chi = closureFunc(xi);
	double dthick = 1.5 - 1.5*chi;
	double dthin  = 1.5*chi - 0.5;
        //CCTK_VInfo(CCTK_THORNSTRING,"Closure: xi = %10.2e E = %10.2e Fx = %10.2e F2 = %10.2e chi = %10.2e",
	//		xi,E,Fx,F2,chi);
        
	// Rest frame energy density and projected fluxes in the optically thick
	// limit.
	double Jthick = 3.0/(2.0*W2 + 1.0) * ((2.0*W2 - 1.0)*E - 2.0*W2*vdotF);
	
	double tHxt = Fxu/W + vx*W/(2.0*W2 + 1.0) * ((4.0*W2 + 1.0)*vdotF - 4.0*W2*E);
	double tHyt = Fyu/W + vy*W/(2.0*W2 + 1.0) * ((4.0*W2 + 1.0)*vdotF - 4.0*W2*E);
	double tHzt = Fzu/W + vz*W/(2.0*W2 + 1.0) * ((4.0*W2 + 1.0)*vdotF - 4.0*W2*E);
  
	Pxx = dthick*(Jthick/3.0*(4.0*W2*vx*vx + gmet.guxx) + W*(tHxt*vx + vx*tHxt));
	Pxy = dthick*(Jthick/3.0*(4.0*W2*vx*vy + gmet.guxy) + W*(tHxt*vy + vx*tHyt));
	Pxz = dthick*(Jthick/3.0*(4.0*W2*vx*vz + gmet.guxz) + W*(tHxt*vz + vx*tHzt));
	Pyy = dthick*(Jthick/3.0*(4.0*W2*vy*vy + gmet.guyy) + W*(tHyt*vy + vy*tHyt));
	Pyz = dthick*(Jthick/3.0*(4.0*W2*vy*vz + gmet.guyz) + W*(tHyt*vz + vy*tHzt));
	Pzz = dthick*(Jthick/3.0*(4.0*W2*vz*vz + gmet.guzz) + W*(tHzt*vz + vz*tHzt));
        
	Pxx = Pxx + dthin*W2*E*fhxu*fhxu;
	Pxy = Pxy + dthin*W2*E*fhxu*fhyu;
	Pxz = Pxz + dthin*W2*E*fhxu*fhzu;
	Pyy = Pyy + dthin*W2*E*fhyu*fhyu;
	Pyz = Pyz + dthin*W2*E*fhyu*fhzu;
	Pzz = Pzz + dthin*W2*E*fhzu*fhzu;
        
	// Calculate the rest frame quantities
        JJ  = JJ0 + JJthick*dthick + JJthin*dthin;
	nH  = cn  + Cthickn*dthick + Cthinn*dthin;
	tHx = -(cv + dthin*Cthinv + dthick*Cthickv)*vx 
	      -(cF + dthick*CthickF)*Fxu
	      -Cthinfh*vdotfh*fhxu;
	tHy = -(cv + dthin*Cthinv + dthick*Cthickv)*vy
	      -(cF + dthick*CthickF)*Fyu
	      -Cthinfh*vdotfh*fhyu;
	tHz = -(cv + dthin*Cthinv + dthick*Cthickv)*vz
	      -(cF + dthick*CthickF)*Fzu
	      -Cthinfh*vdotfh*fhzu;
        
	if (isnan(Pxx)) {
#pragma omp critical
	  cerr << "Pxx is a NaN. " << dthin << " " 
	       << chi << " " << xi << " " << Fx << " " << E << " " << v2 << " " << vx << " " << vxl << "\n" ;
	  CCTK_WARN(0,"aborting!");
	}	

	// Return the closure factor if anyone is interested
	return xi; 

}


void Closure::getSourceUpdate(double gE, double gN, double gFx, double gFy, double gFz,
                              double kapa, double kaps, double eta, double dtau, 
			      double *Ep, double *Np, double *Fxp, double *Fyp, double *Fzp, double *Jout, 
			      double *dJ, double *dN, double *dTau, double *dSx, double *dSy, double *dSz) {
        // Given the explicit parts of the update (gE,gFx,etc.), calculate the implicit update to E and F_i
        // Assumes an explicit closure factor		
	double chi = closureFunc(xi);
	double dthin  = 1.5*chi - 0.5;
	double dthick = 1.5 - 1.5*chi;
	
	// Build the components of the rest frame energy density
	double JdE = W2 + dthin*W2*vdotfh*vdotfh + dthick*(3.0-2.0*W2)*(W2-1.0)/(2.0*W2+1.0);	
	double JdF = (4.0*dthick*W2*(W2-1.0)/(2.0*W2+1.0) - 2.0*W2); // Derivative of wrt to any F_i component is JdF*v^i	

	// Build the components of the E source term
	double SE0  =  dtau*W*eta;
	double SEdE =  dtau*W*kaps*(1.0 - JdE); // Diagonal portion that goes to zero with v
	double SEdF = -dtau*W*(kapa + kaps*(1.0 + JdF)); // Derivative of source term wrt to any F_i component is SEdF*v^i
        double dE   = 1.0 + dtau*W*kapa; // Diagonal term that is finite in zero velocity limit

	// Build the components of the F_i source terms
	// To get any SF_idx take SFdx*v_i
	double SF0  =  dtau*W*eta; 
	double SFdE = -dtau*(W*kaps*JdE + W/(2.0*W2+1.0)*dthick*(kaps+kapa)*(3.0-2.0*W2)); 
	double SFdF = -dtau*(W*kaps*JdF + W/(2.0*W2+1.0)*dthick*(kaps+kapa)*(2.0*W2-1.0)); // Derivative of source term j wrt to 
		   	       						                  // any F_i component is SFdF*v^i*v_j
	//double dF    = 1.0 + dtau*W*(kapa+kaps)*(1.0 - dthick*v2 - dthin*E/sqrt(F2)*vdotfh);
	double dF = 1.0 + dtau*W*(kapa+kaps)*(1.0 - dthick*v2 - dthin*vdotfh); // Assume E = sqrt(F) when last term matters
	
	// Build the RHS
	double bE = gE  + SE0;
	double bx = gFx + SF0*vxl;
	double by = gFy + SF0*vyl;
	double bz = gFz + SF0*vzl;
	
	// Solve Matrix equation directly
	// This is relatively simple because of the high degree of symmetry of the linear equations 
	// All of this is Mathematica generated
	double detA = dF*dF*( (dF + SFdF*v2)*(dE + SEdE) - SEdF*SFdE*v2 );
        
	double bdotv = bx*vx + by*vy + bz*vz;
	double aE = -dF*dF/dE*( bdotv*dE*SEdF + bE   *(dF*SEdE - (SEdF*SFdE - SEdE*SFdF)*v2) );
	double aF = -dF*      ( bE   *dF*SFdE + bdotv*(dE*SFdF - (SEdF*SFdE - SEdE*SFdF)   ) );
	
	*Ep  = bE/dE +     aE/detA; // The first term is just the diagonal in the zero velocity limit
	*Fxp = bx/dF + vxl*aF/detA;
	*Fyp = by/dF + vyl*aF/detA;
	*Fzp = bz/dF + vzl*aF/detA;
        
	// Calculate the change in energy in the rest frame
	double Fdotv = (vx*(*Fxp) + vy*(*Fyp) + vz*(*Fzp));
	double J = JdE*(*Ep) + JdF*Fdotv;
	*Jout = J;
        *dJ = dtau*(eta-kapa*J);
        *dTau = -(SEdE*(*Ep) + SEdF*Fdotv - SE0)     - (dE - 1.0)*(*Ep) ; 
        *dSx  = -(SFdE*(*Ep) + SFdF*Fdotv - SF0)*vxl - (dF - 1.0)*(*Fxp); //Lowered index
        *dSy  = -(SFdE*(*Ep) + SFdF*Fdotv - SF0)*vyl - (dF - 1.0)*(*Fyp);
        *dSz  = -(SFdE*(*Ep) + SFdF*Fdotv - SF0)*vzl - (dF - 1.0)*(*Fzp);
	
	// Calculate things required for number source term 	
  	double Jf      = W2*(*Ep - 2.0*Fdotv); 
	double Jfthin  = W2*E*vdotfh*vdotfh;
        double Jfthick = (W2-1.0)/(1.0+2.0*W2)*(4.0*W2*Fdotv + (3.0-2.0*W2)*(*Ep));
	double cn = W*Jf + W*(Fdotv - *Ep);
	double Cthickn = W*Jfthick;
	double Cthinn  = W*Jfthin;
        
	double nHf = cn + Cthickn*dthick + Cthinn*dthin;
	*Np = (gN + dtau*eta + dtau*kapa*nHf)/(1.0 + dtau*W*kapa);
	*dN = dtau*(eta - kapa*W*(*Np) + kapa*nHf);
        	
	return;

}

void Closure::getAsymptoticWaveSpeeds(int direction, double * lam) {
	
	double alpha = gmet.getAlpha();
	double beta;
	double vel;
	double gaau;
	if (direction==0) {
	  beta = gmet.getBetax();
	  vel  = vx;
	  gaau = gmet.guxx;
	} else if (direction==1) {
	  beta = gmet.getBetay();
	  vel  = vy;
	  gaau = gmet.guyy;
        } else if (direction==2) {
	  beta = gmet.getBetaz();
	  vel  = vz;
	  gaau = gmet.guzz;
	} else {
	  beta = 0.0;
	  vel  = 0.0;
	  gaau = 1.0;
	} 

	lam[0] = -beta + 4.0*alpha*W2*vel/(2.0*W2+1.0);
	         //-2.0*alpha/(2.0*W2+1.0)*sqrt(2.0*W2*(W2-1.0)*(2.0*W2-3.0+2.0*gaau*pow(W2-1.0,2.0)));
	lam[1] = -beta + 4.0*alpha*W2*vel/(2.0*W2+1.0);
	         //+ 2.0*alpha/(2.0*W2+1.0)*sqrt(2.0*W2*(W2-1.0)*(2.0*W2-3.0+2.0*gaau*pow(W2-1.0,2.0)));
	lam[2] = lam[0];
	lam[3] = lam[1];
}

void Closure::getWaveSpeeds(int direction,double * lam) {
        // follows: Shibata, Kiuchi, Sekiguchi, Suwa Prog. Th. Phys. 1222 vol 125 (2011)
        // Section 6.4ff
 	//double TINY = 1.e-45;
						
	double alpha = gmet.getAlpha();
	double beta,vel,gaau,Fau;
	if (direction==0) {
	  beta = gmet.getBetax();
	  vel  = vx;
	  gaau = gmet.guxx;
	  Fau  = Fxu;
	} else if (direction==1) {
	  beta = gmet.getBetay();
	  vel  = vy;
	  gaau = gmet.guyy;
	  Fau  = Fyu;
        } else if (direction==2) {
	  beta = gmet.getBetaz();
	  vel  = vz;
	  gaau = gmet.guzz;
	  Fau  = Fzu;
	} else {
	  beta = 0.0;
	  vel  = 0.0;
	  gaau = 1.0;
	  Fau  = 0.0;
	} 

        double lamthinmin = min(-beta - alpha*abs(Fau)/F,
				-beta + alpha*E*Fau/F2);
        double lamthinmax = max(-beta + alpha*abs(Fau)/F,
				-beta + alpha*E*Fau/F2);	
        // Guarantee that the speed of light is not exceeded
	lamthinmin = -beta - alpha*abs(Fau)/F;
        lamthinmax = -beta + alpha*abs(Fau)/F;

        double p = alpha*vel/W;
	double rt = sqrt(alpha*alpha*gaau*(2.0*W*W + 1.0) - 2.0*W*W*p*p);
        double lamthickmin = min(-beta+(2.0*W*W*p-rt)/(2.0*W*W+1.0),-beta+p);	
        double lamthickmax = max(-beta+(2.0*W*W*p+rt)/(2.0*W*W+1.0),-beta+p);	

	double chi = closureFunc(xi); 
	double dthick = 1.5 - 1.5*chi;
	double dthin  = 1.5*chi - 0.5;
	
	// Interpolate the wave speeds as in Shibata 
        lam[0] = dthin*lamthinmin + dthick*lamthickmin;	
        lam[1] = dthin*lamthinmax + dthick*lamthickmax;
        lam[2] = lam[0];
        lam[3] = lam[1];
	
	return;

}

void Closure::getMomentumFluxes(double dlnadx, double dlnady, double dlnadz,
                                double ndW, double dWx, double dWy, double dWz,
                                double aa,  double bx,  double by,  double bz, 
                                double cxx, double cxy, double cxz, double cyy, 
                                double cyz, double czz,
				double *fE, double *fN, double *fFx, double *fFy, double *fFz) {
        
	double TINY = 1.e-45;
	
	// Calculate the optically thin/thick mixing
	double chi = closureFunc(xi); 
	double dthick = 1.5 - 1.5*chi;
	double dthin  = 1.5*chi - 0.5;
	
	// Calculate pieces of J 	
  	//JJ0 = W2*(E - 2.0*vdotF); 
	//JJthin  = W2*E*vdotfh*vdotfh;
        //JJthick = (W2-1.0)/(1.0+2.0*W2)*(4.0*W2*vdotF + (3.0-2.0*W2)*E);

	// Calculate pieces of H^alpha H_alpha 
	//double cv = W*JJ0;
	//double cF =-W;
        //
	//double Cthickv = W*JJthick + W/(2.0*W2+1.0)*((3.0-2.0*W2)*E + (2.0*W2-1.0)*vdotF);
	//double CthickF = W*v2;		
	//
	//double Cthinv  = W*JJthin;
	//double Cthinfh = W*E*vdotfh;

	//// Rest Frame energy density 
	//double J   = JJ0 + dthin*JJthin + dthick*JJthick;
	//
	// Rest frame flux and contractions with velocity 	
	//double tHx = -(cv + dthin*Cthinv + dthick*Cthickv)*vx 
	//             -(cF + dthick*CthickF)*Fxu
	//             -Cthinfh*vdotfh*fhxu;
	//double tHy = -(cv + dthin*Cthinv + dthick*Cthickv)*vy
	//             -(cF + dthick*CthickF)*Fyu
	//             -Cthinfh*vdotfh*fhyu;
	//double tHz = -(cv + dthin*Cthinv + dthick*Cthickv)*vz
	//             -(cF + dthick*CthickF)*Fzu
	//             -Cthinfh*vdotfh*fhzu;
        double J = JJ; 
	double vtH =  vxl*tHx + vyl*tHy + vzl*tHz; 
        
	// Spatially projected rest frame projection operator
	double thxx = gmet.guxx + W2*vx*vx;
	double thxy = gmet.guxy + W2*vx*vy;
	double thxz = gmet.guxz + W2*vx*vz;
	double thyy = gmet.guyy + W2*vy*vy;
	double thyz = gmet.guyz + W2*vy*vz;
	double thzz = gmet.guzz + W2*vz*vz;
	
	// Rest frame pressure tensor and contractions with velocity
	// This is a little bit of a fudge, since the closure relation 
	// we are using is not consistent
	// Also note that Kxx is *not* the extrinsic curvature here  
	double H2 = gmet.contractGllAuBu(tHx,tHy,tHz,tHx,tHy,tHz) + TINY;
	double Kxx = dthick*J/3.0*thxx + dthin*J*tHx*tHx/H2;      
	double Kxy = dthick*J/3.0*thxy + dthin*J*tHx*tHy/H2;      
	double Kxz = dthick*J/3.0*thxz + dthin*J*tHx*tHz/H2;      
	double Kyy = dthick*J/3.0*thyy + dthin*J*tHy*tHy/H2;      
	double Kyz = dthick*J/3.0*thyz + dthin*J*tHy*tHz/H2;      
	double Kzz = dthick*J/3.0*thzz + dthin*J*tHz*tHz/H2;     
	
	double vKx = vxl*Kxx + vyl*Kxy + vzl*Kxz; 
	double vKy = vxl*Kxy + vyl*Kyy + vzl*Kyz; 
	double vKz = vxl*Kxz + vyl*Kyz + vzl*Kzz; 
	
	double vvK = vxl*vKx + vyl*vKy + vzl*vKz;

	// Third moment in the rest frame 
	double Lxxx = dthin*tHx*tHx*tHx/H2 + dthick*3*tHx*thxx/5;
	double Lyyy = dthin*tHy*tHy*tHy/H2 + dthick*3*tHy*thyy/5;
	double Lzzz = dthin*tHz*tHz*tHz/H2 + dthick*3*tHz*thzz/5;
	double Lxxy = dthin*tHx*tHx*tHy/H2 + dthick*(2*tHx*thxy + tHy*thxx)/5;
	double Lxxz = dthin*tHx*tHx*tHz/H2 + dthick*(2*tHx*thxz + tHz*thxx)/5;
	double Lxyy = dthin*tHx*tHy*tHy/H2 + dthick*(2*tHy*thxy + tHx*thyy)/5;
	double Lyyz = dthin*tHy*tHy*tHz/H2 + dthick*(2*tHy*thyz + tHz*thyy)/5;
	double Lxzz = dthin*tHx*tHz*tHz/H2 + dthick*(tHx*thzz + 2*tHz*thxz)/5;
	double Lyzz = dthin*tHy*tHz*tHz/H2 + dthick*(tHy*thzz + 2*tHz*thyz)/5;
	double Lxyz = dthin*tHx*tHy*tHz/H2 + dthick*(tHx*thyz + tHy*thxz + tHz*thxy)/5;

        double vLxx = vxl*Lxxx + vyl*Lxxy + vzl*Lxxz;	
        double vLxy = vxl*Lxxy + vyl*Lxyy + vzl*Lxyz;	
        double vLxz = vxl*Lxxz + vyl*Lxyz + vzl*Lxzz;	
        double vLyy = vxl*Lxyy + vyl*Lyyy + vzl*Lyyz;	
        double vLyz = vxl*Lxyz + vyl*Lyyz + vzl*Lyzz;	
        double vLzz = vxl*Lxzz + vyl*Lyzz + vzl*Lzzz;
	
	double vvLx = vxl*vLxx + vyl*vLxy + vzl*vLxz;
	double vvLy = vxl*vLxy + vyl*vLyy + vzl*vLyz;
	double vvLz = vxl*vLxz + vyl*vLyz + vzl*vLzz;
		
        double vvvL = vxl*vvLx + vyl*vvLy + vzl*vvLz;
		
	// Projections in the lab frame 
	double Q  = W2*(W*J + 3*vtH) + 3*W*vvK + vvvL;

	double Rx = W2*(W*J + 2*vtH)*vx + W2*tHx + W*vvK*vx + 2*W*vKx + vvLx;	
	double Ry = W2*(W*J + 2*vtH)*vy + W2*tHy + W*vvK*vy + 2*W*vKy + vvLy;	
	double Rz = W2*(W*J + 2*vtH)*vz + W2*tHz + W*vvK*vz + 2*W*vKz + vvLz;	
	
	double Sxx = W2*vx*vx*(W*J + vtH) + W2*(tHx*vx + vx*tHx) + W*(vKx*vx + vx*vKx) + W*Kxx + vLxx; 
	double Sxy = W2*vx*vy*(W*J + vtH) + W2*(tHx*vy + vx*tHy) + W*(vKx*vy + vx*vKy) + W*Kxy + vLxy; 
	double Sxz = W2*vx*vz*(W*J + vtH) + W2*(tHx*vz + vx*tHz) + W*(vKx*vz + vx*vKz) + W*Kxz + vLxz; 
	double Syy = W2*vy*vy*(W*J + vtH) + W2*(tHy*vy + vy*tHy) + W*(vKy*vy + vy*vKy) + W*Kyy + vLyy; 
	double Syz = W2*vy*vz*(W*J + vtH) + W2*(tHy*vz + vy*tHz) + W*(vKy*vz + vy*vKz) + W*Kyz + vLyz; 
	double Szz = W2*vz*vz*(W*J + vtH) + W2*(tHz*vz + vz*tHz) + W*(vKz*vz + vz*vKz) + W*Kzz + vLzz; 
        
	double Txxx = W3*J*vx*vx*vx + 3*W2*tHx*vx*vx + 3*W*Kxx*vx + Lxxx;		
	double Tyyy = W3*J*vy*vy*vy + 3*W2*tHy*vy*vy + 3*W*Kyy*vy + Lyyy;		
	double Tzzz = W3*J*vz*vz*vz + 3*W2*tHz*vz*vz + 3*W*Kzz*vz + Lzzz;		
	
	double Txxy = W3*J*vx*vx*vy + W2*(2*tHx*vx*vy + vx*vx*tHy) + W*(Kxx*vy + 2*Kxy*vx) + Lxxy;
	double Txxz = W3*J*vx*vx*vz + W2*(2*tHx*vx*vz + vx*vx*tHz) + W*(Kxx*vz + 2*Kxz*vx) + Lxxz;

	double Txyy = W3*J*vx*vy*vy + W2*(tHx*vy*vy + 2*vx*tHy*vy) + W*(2*Kxy*vy + Kyy*vx) + Lxyy;
	double Txzz = W3*J*vx*vz*vz + W2*(tHx*vz*vz + 2*vx*tHz*vz) + W*(2*Kxz*vz + Kzz*vx) + Lxzz;
	
	double Tyyz = W3*J*vy*vy*vz + W2*(2*tHy*vy*vz + vy*vy*tHz) + W*(2*Kyz*vy + Kyy*vz) + Lyyz;
	
	double Tyzz = W3*J*vy*vz*vz + W2*(tHy*vz*vz + 2*vy*tHz*vz) + W*(2*Kyz*vz + Kzz*vy) + Lyzz;
	double Txyz = W3*J*vx*vy*vz + W2*(tHx*vy*vz + vx*tHy*vz + vx*vy*tHz) 
	            + W*(Kxy*vz + Kxz*vy + vz*Kxy) + Lxyz;		
       
       // Now calculate the momentum space fluxes
       // (NOTE: These fluxes are just the contraction of \gamma_{i\delta} and 
       //  n_\delta with M^{\alpha \beta \delta}u_{\alpha;\beta}, so the radiation 
       //  fluxes eventually require a minus sign)
       double fFxu,fFyu,fFzu;
       *fE = 0; fFxu = 0; fFyu = 0; fFzu = 0;
       *fFx = 0; *fFy = 0; *fFz = 0;
       *fN = 0;

       // Velocity dependent fluxes
       *fE  += Q*ndW + Rx*dWx + Ry*dWy + Rz*dWz;
       *fN  += E*ndW + Fxu*dWx + Fyu*dWy + Fzu*dWz;
       fFxu += Sxx*dWx + Sxy*dWy + Sxz*dWz + Rx*ndW;
       fFyu += Sxy*dWx + Syy*dWy + Syz*dWz + Ry*ndW;
       fFzu += Sxz*dWx + Syz*dWy + Szz*dWz + Rz*ndW;
       
       *fE  +=  Rx*bx + Ry*by + Rz*bz - aa*Q - Sxx*cxx - Syy*cyy - Szz*czz 
               - 2*(Sxy*cxy + Sxz*cxz + Syz*cyz);
       *fN  += -aa*E + Fxu*bx + Fyu*by + Fzu*bz - Pxx*cxx - Pyy*cyy - Pzz*czz 
               - 2*(Pxy*cxy + Pxz*cxz + Pyz*cyz);
       fFxu += - aa*Rx + 2*(bx*Sxx + by*Sxy + bz*Sxz) -   (Txxx*cxx + Txyy*cyy + Txzz*czz)
                                                      - 2*(Txxy*cxy + Txxz*cxz + Txyz*cyz);
       fFyu += - aa*Ry + 2*(bx*Sxy + by*Syy + bz*Syz) -   (Txxy*cxx + Tyyy*cyy + Tyzz*czz)
                                                      - 2*(Txyy*cxy + Txyz*cxz + Tyyz*cyz);
       fFzu += - aa*Rz + 2*(bx*Sxz + by*Syz + bz*Szz) -   (Txxz*cxx + Tyyz*cyy + Tzzz*czz)
                                                      - 2*(Txyz*cxy + Txzz*cxz + Tyzz*cyz);
       
       // Curvature fluxes
       *fE  +=   W*gmet.contractKllTuuSym(Sxx,Sxy,Sxz,Syy,Syz,Szz) 
               - W*(Rx*dlnadx + Ry*dlnady + Rz*dlnadz);	  
       *fN  += W*( gmet.contractKllTuuSym(Pxx,Pxy,Pxz,Pyy,Pyz,Pzz)
                 - gmet.contractAlBu(dlnadx,dlnady,dlnadz,Fxu,Fyu,Fzu));

       fFxu += - W*(Sxx*dlnadx + Sxy*dlnady + Sxz*dlnadz)
               + W*gmet.contractKllTuuSym(Txxx,Txxy,Txxz,Txyy,Txyz,Txzz);	        
       fFyu += - W*(Sxy*dlnadx + Syy*dlnady + Syz*dlnadz) 
               + W*gmet.contractKllTuuSym(Txxy,Txyy,Txyz,Tyyy,Tyyz,Tyzz);	        
       fFzu += - W*(Sxz*dlnadx + Syz*dlnady + Szz*dlnadz) 
               + W*gmet.contractKllTuuSym(Txxz,Txyz,Txzz,Tyyz,Tyzz,Tzzz);	         
       
       gmet.lowerAu(fFxu,fFyu,fFzu,fFx,fFy,fFz);  
       
       return;
}

inline double Closure::zeroFunc(double xi){	
	double chi = closureFunc(xi);
	double athin  = 1.5*chi - 0.5;
	double athick = 1.5 - 1.5*chi;
		
	double J = JJ0 + athin*JJthin + athick*JJthick;

	double g = HH0 + athick*HHthick + athin*HHthin + athick*athin*HHthickthin 
	          + athick*athick*HHthickthick + athin*athin*HHthinthin;
						
	return (J*J*xi*xi - g)/(E*E); 
}

inline double Closure::zeroFunc(double xi,double *dfdxi){	
	double dchidxi;
  	double chi = closureFunc(xi,&dchidxi);
	double athin  = 1.5*chi - 0.5;
	double athick = 1.5 - 1.5*chi;
	double dathindxi  =  1.5*dchidxi;
	double dathickdxi = -1.5*dchidxi; 
		
	double J = JJ0 + athin*JJthin + athick*JJthick;
	double dJdxi = JJthin*dathindxi + JJthick*dathickdxi;

	double g = HH0 + athick*HHthick + athin*HHthin + athick*athin*HHthickthin 
	          + athick*athick*HHthickthick + athin*athin*HHthinthin;
	double dgdathin  = (HHthin  + HHthickthin*athick + 2.0*HHthinthin*athin   );
	double dgdathick = (HHthick + HHthickthin*athin  + 2.0*HHthickthick*athick);
	
	*dfdxi = 2.0*J*dJdxi*xi*xi + 2.0*J*J*xi - dathindxi*dgdathin - dathickdxi*dgdathick;
	*dfdxi = *dfdxi/(E*E);
	return (J*J*xi*xi - g)/(E*E); 
}


//inline double Closure::getEquations() {
//	double fE,fFx,fFy,fFz,fClose;
//	
//	return;
//}

