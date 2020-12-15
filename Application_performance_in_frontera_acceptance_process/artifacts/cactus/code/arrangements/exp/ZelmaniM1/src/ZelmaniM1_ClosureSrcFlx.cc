#include <cmath>
#include <math.h>
#include <iostream>
#include "ZelmaniM1_Closure.hh"


void Closure::getFluxes (int direction, double nuNin, double *fenu, 
    double *fnnu, double *fFx, double *fFy, double *fFz) const {
  double Paxl,Payl,Pazl;
  double alpha = gmet.getAlpha();
  double beta,Fau,va,tHa;
  direction = direction%3;
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
  } else {
    gmet.lowerAu(Pxz,Pyz,Pzz,&Paxl,&Payl,&Pazl);
    beta = gmet.getBetaz();
    Fau  = Fzu; tHa = tHz;
    va   = vz;
  } 

  *fenu = (alpha*Fau  - beta*E );
  *fnnu = (alpha*W*va*JJ + alpha*tHa - beta*nuNin);

  *fFx  = (alpha*Paxl - beta*Fx); 
  *fFy  = (alpha*Payl - beta*Fy); 
  *fFz  = (alpha*Pazl - beta*Fz); 
}

double getAsymptoticDiffusiveFlux(double kappa, double ida, 
    const Closure& cL, const Closure& cR, const ThreeTensor::Metric& gmetC) { 
  double sdetgR = sqrt(cR.gmet.getDetg());
  double sdetgL = sqrt(cL.gmet.getDetg());
  double sdetgC = sqrt(gmetC.getDetg());
  return -sdetgC*gmetC.getAlpha()*sqrt(cL.W*cR.W)/(3.0*kappa) 
      * (cR.JJ/sdetgR - cL.JJ/sdetgL)*ida;
}

void Closure::getAsymptoticAdvectiveFluxes(int direction, double *fenu,
    double *fFx, double *fFy, double *fFz) const {
  //double Paxl,Payl,Pazl;
  double alpha = gmet.getAlpha();
  double beta,va;
  direction = direction%3;
  if (direction==0) {
    //gmet.lowerAu(Pxx,Pxy,Pxz,&Paxl,&Payl,&Pazl);
    beta = gmet.getBetax();
    va   = vx;
  } else if (direction==1) { 
    //gmet.lowerAu(Pxy,Pyy,Pyz,&Paxl,&Payl,&Pazl);
    beta = gmet.getBetay();
    va   = vy;
  } else {
    //gmet.lowerAu(Pxz,Pyz,Pzz,&Paxl,&Payl,&Pazl);
    beta = gmet.getBetaz();
    va   = vz;
  } 
  
  //double Jthick = 3.0/(2.0*W2 + 1.0)*((2.0*W2 - 1.0)*E - 2.0*W2*vdotF);
  //Jthick = 3.0/(4.0*W2-1.0)*E;
    
  *fenu = (alpha*4.0*W2/3.0*va*JJ  - beta*E);

  *fFx  = 0.0;//(alpha*Paxl - beta*Fx); 
  *fFy  = 0.0;//(alpha*Payl - beta*Fy); 
  *fFz  = 0.0;//(alpha*Pazl - beta*Fz); 
}
  
void Closure::getSourceUpdate(double gE, double gN, double gFx, double gFy, double gFz,
    double kapa, double kaps, double eta, double dtau, 
    double *Ep, double *Np, double *Fxp, double *Fyp, double *Fzp, double *Jout, 
    double *dJ, double *dN, double *dTau, double *dSx, double *dSy, double *dSz) const {
  // Given the explicit parts of the update (gE,gFx,etc.), calculate the 
  // implicit update to E and F_i
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

std::vector<double> Closure::getAsymptoticWaveSpeeds(int direction) const {
  
  double alpha = gmet.getAlpha();
  double beta;
  double vel;
  //double gaau;
  direction = direction%3;
  if (direction==0) {
    beta = gmet.getBetax();
    vel  = vx;
    //gaau = gmet.guxx;
  } else if (direction==1) {
    beta = gmet.getBetay();
    vel  = vy;
    //gaau = gmet.guyy;
  } else {
    beta = gmet.getBetaz();
    vel  = vz;
    //gaau = gmet.guzz;
  } 

  std::vector<double> lam;
  lam.push_back(-beta + 4.0*alpha*W2*vel/(2.0*W2+1.0));
           //-2.0*alpha/(2.0*W2+1.0)*sqrt(2.0*W2*(W2-1.0)*(2.0*W2-3.0+2.0*gaau*pow(W2-1.0,2.0)));
  lam.push_back(-beta + 4.0*alpha*W2*vel/(2.0*W2+1.0));
           //+ 2.0*alpha/(2.0*W2+1.0)*sqrt(2.0*W2*(W2-1.0)*(2.0*W2-3.0+2.0*gaau*pow(W2-1.0,2.0)));
  return lam;
}

std::vector<double> Closure::getWaveSpeeds(int direction,bool fixed_closure) const {
        // follows: Shibata, Kiuchi, Sekiguchi, Suwa Prog. Th. Phys. 1222 vol 125 (2011)
        // Section 6.4ff
 	//double TINY = 1.e-45;
            
  double alpha = gmet.getAlpha();
  double beta,vel,gaau,Fau,fha;
  direction = direction%3;
  if (direction==0) {
    beta = gmet.getBetax();
    vel  = vx;
    gaau = gmet.guxx;
    Fau  = Fxu;
    fha = fhx;
  } else if (direction==1) {
    beta = gmet.getBetay();
    vel  = vy;
    gaau = gmet.guyy;
    Fau  = Fyu;
    fha = fhy;
  } else {
    beta = gmet.getBetaz();
    vel  = vz;
    gaau = gmet.guzz;
    Fau  = Fzu;
    fha = fhz;
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
  std::vector<double> lam;
  if (!fixed_closure) {
    lam.push_back(dthin*lamthinmin + dthick*lamthickmin); 
    lam.push_back(dthin*lamthinmax + dthick*lamthickmax);
  } else {
    lam.push_back(-sqrt(dthick/3.0 + dthin*fha*fha)); 
    lam.push_back( sqrt(dthick/3.0 + dthin*fha*fha)); 
  }
  return lam;
}

void Closure::getMomentumFluxes(double dlnadx, double dlnady, double dlnadz,
    double ndW, double dWx, double dWy, double dWz, double aa,  double bx,  
    double by,  double bz, double cxx, double cxy, double cxz, double cyy, 
    double cyz, double czz, double *fE, double *fN, double *fFx, double *fFy, 
    double *fFz) const {
        
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
  
}



void Closure::getMomentumFluxes(double dlnadx, double dlnady, double dlnadz,
    double *fE, double *fN, double *fFx, double *fFy, double *fFz) const {
        
  // Calculate the optically thin/thick mixing
  double chi = closureFunc(xi); 
  double dthick = 1.5 - 1.5*chi;
  double dthin  = 1.5*chi - 0.5;
  
  // Third projected moment 
  double Lxxx = dthin*Fxu*Fxu*Fxu/F2 + dthick*3*Fxu*gmet.guxx/5;
  double Lyyy = dthin*Fyu*Fyu*Fyu/F2 + dthick*3*Fyu*gmet.guyy/5;
  double Lzzz = dthin*Fzu*Fzu*Fzu/F2 + dthick*3*Fzu*gmet.guzz/5;
  double Lxxy = dthin*Fxu*Fxu*Fyu/F2 + dthick*(2*Fxu*gmet.guxy+Fyu*gmet.guxx)/5;
  double Lxxz = dthin*Fxu*Fxu*Fzu/F2 + dthick*(2*Fxu*gmet.guxz+Fzu*gmet.guxx)/5;
  double Lxyy = dthin*Fxu*Fyu*Fyu/F2 + dthick*(2*Fyu*gmet.guxy+Fxu*gmet.guyy)/5;
  double Lyyz = dthin*Fyu*Fyu*Fzu/F2 + dthick*(2*Fyu*gmet.guyz+Fzu*gmet.guyy)/5;
  double Lxzz = dthin*Fxu*Fzu*Fzu/F2 + dthick*(2*Fzu*gmet.guxz+Fxu*gmet.guzz)/5;
  double Lyzz = dthin*Fyu*Fzu*Fzu/F2 + dthick*(2*Fzu*gmet.guyz+Fyu*gmet.guzz)/5;
  double Lxyz = dthin*Fxu*Fyu*Fzu/F2 
      + dthick*(Fxu*gmet.guyz + Fyu*gmet.guxz + Fzu*gmet.guxy)/5;
       
  // Now calculate the momentum space fluxes
  // (NOTE: These fluxes are just the contraction of \gamma_{i\delta} and 
  //  n_\delta with M^{\alpha \beta \delta}u_{\alpha;\beta}, so the radiation 
  //  fluxes eventually require a minus sign)
  
  // Curvature fluxes
  *fE = gmet.contractKllTuuSym(Pxx,Pxy,Pxz,Pyy,Pyz,Pzz) 
      - (Fxu*dlnadx + Fyu*dlnady + Fzu*dlnadz);   
  *fN = gmet.contractKllTuuSym(Pxx,Pxy,Pxz,Pyy,Pyz,Pzz)
      - (Fxu*dlnadx + Fyu*dlnady + Fzu*dlnadz);

  double fFxu,fFyu,fFzu;
  fFxu = - (Pxx*dlnadx + Pxy*dlnady + Pxz*dlnadz)
      + gmet.contractKllTuuSym(Lxxx,Lxxy,Lxxz,Lxyy,Lxyz,Lxzz);         
  fFyu = - (Pxy*dlnadx + Pyy*dlnady + Pyz*dlnadz) 
      + gmet.contractKllTuuSym(Lxxy,Lxyy,Lxyz,Lyyy,Lyyz,Lyzz);         
  fFzu = - (Pxz*dlnadx + Pyz*dlnady + Pzz*dlnadz) 
      + gmet.contractKllTuuSym(Lxxz,Lxyz,Lxzz,Lyyz,Lyzz,Lzzz);          
  
  gmet.lowerAu(fFxu,fFyu,fFzu,fFx,fFy,fFz);  
  
}



