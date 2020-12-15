#include <cmath>
#include <math.h>
#include <iostream>
#include "ZelmaniM1_Closure.hh"
// These need to go to make this general
//#include "cctk.h"
//#include "cctk_Parameters.h"
//#include "cctk_Arguments.h"
//#include "cctk_Functions.h"
//#include "carpet.hh"

double Closure::zFunc(double xi) {
  const double chi = closureFunc(xi);
  const double athin  = 1.5*chi - 0.5;
  const double athick = 1.5 - 1.5*chi;
  const double J = JJ0 + athin*JJthin + athick*JJthick;
  const double g = HH0 + athick*HHthick + athin*HHthin 
      + athick*athin*HHthickthin + athick*athick*HHthickthick 
      + athin*athin*HHthinthin;
  return (J*J*xi*xi - g)/(Eclosure*Eclosure); 
}

double Closure::zFuncWithD(double xi, double *dfdxi) {
  double dchidxi;
  const double chi = closureFunc(xi,&dchidxi);
  const double athin  = 1.5*chi - 0.5;
  const double athick = 1.5 - 1.5*chi;
  const double dathindxi  =  1.5*dchidxi;
  const double dathickdxi = -1.5*dchidxi; 
  const double J = JJ0 + athin*JJthin + athick*JJthick;
  const double dJdxi = JJthin*dathindxi + JJthick*dathickdxi;
  const double g = HH0 + athick*HHthick + athin*HHthin 
    + athick*athin*HHthickthin + athick*athick*HHthickthick 
    + athin*athin*HHthinthin;
  const double dgdathin  = HHthin  + HHthickthin*athick 
      + 2.0*HHthinthin*athin;
  const double dgdathick = HHthick + HHthickthin*athin 
      + 2.0*HHthickthick*athick;
  *dfdxi = 2.0*J*dJdxi*xi*xi + 2.0*J*J*xi - dathindxi*dgdathin 
      - dathickdxi*dgdathick;
  *dfdxi = *dfdxi/(Eclosure*Eclosure);
  return (J*J*xi*xi - g)/(Eclosure*Eclosure); 
}

Closure::Closure(double vx, double vy, double vz, ThreeTensor::Metric gamma,
    double lambda, double fhxi, double fhyi, double fhzi) :
  vx(vx), vy(vy), vz(vz), gmet(gamma), lambda(lambda), 
  fhxi(fhxi), fhyi(fhyi), fhzi(fhzi), TINY(1.e-45) {
  
  lambda = max(min(lambda,1.0),0.0);
   
  //Make sure this is a unit vector
  const double fhm = 
      sqrt(gmet.contractGuuAlBl(fhxi,fhyi,fhzi,fhxi,fhyi,fhzi) + TINY);
  fhxi = fhxi/fhm;
  fhyi = fhyi/fhm;
  fhzi = fhzi/fhm;
  
  gmet.lowerAu(vx,vy,vz,&vxl,&vyl,&vzl);  
  v2 = gmet.contractAlBu(vxl,vyl,vzl,vx,vy,vz);
  W = 1.0/sqrt(1.0-v2);
  W2 = W*W;
  W3 = W*W2;
  debout="";
}

double Closure::closureFunc(double xi) const {
  //return (3.0 + 4.0*xi*xi)/(5.0 + 2.0*sqrt(4.0 - 3.0*xi*xi));
  return 1.0/3.0 + xi*xi*(6.0 - 2.0*xi + 6.0*xi*xi)/15.0;
}

double Closure::closureFunc(double xi, double *dchidxi) const {
  // Levermore 84 closure
  //double rt = sqrt(4.0 - 3.0*xi*xi);
  //*dchidxi = 2.0*xi/rt;
  //return (3.0 + 4.0*xi)/(5.0 + 2.0*rt);
  // Maximum Entropy closure
  *dchidxi = 2.0*xi*(6.0 - 2.0*xi + 6.0*xi*xi)/15.0
    + xi*xi*(12.0*xi - 2.0)/15.0;
  return 1.0/3.0 + xi*xi*(6.0 - 2.0*xi + 6.0*xi*xi)/15.0;
} 

double Closure::setClosure(double Ei, double Fxi, double Fyi, double Fzi,
  double xi0, bool Force_closure) { 
  // Given the rest frame E and F_i, as well as the local three 
  // velocity v^i and the three metric gamma_{ij}, this routine 
  // solves the closure relation for the radiation stress tensor P^{ij}
  // in the lab frame and returns the Eddington factor and components 
  // of the stress tensor. 
  
  E = Ei; 
  Fx = Fxi; 
  Fy = Fyi; 
  Fz = Fzi; 
  gmet.raiseAl(Fx,Fy,Fz,&Fxu,&Fyu,&Fzu);
  F2 = gmet.contractAlBu(Fx,Fy,Fz,Fxu,Fyu,Fzu) + TINY;
  F  = sqrt(F2);
  vdotF = gmet.contractAlBu(Fx,Fy,Fz,vx,vy,vz) + TINY; 
   
  fhx = (1.0-lambda)*fhxi + lambda*Fx/F; 
  fhy = (1.0-lambda)*fhyi + lambda*Fy/F; 
  fhz = (1.0-lambda)*fhzi + lambda*Fz/F;
  gmet.raiseAl(fhx,fhy,fhz,&fhxu,&fhyu,&fhzu);
  
  vdotfh = gmet.contractAlBu(fhx,fhy,fhz,vx,vy,vz) + TINY;  
  Fdotfh = gmet.contractAlBu(Fx,Fy,Fz,fhxu,fhyu,fhzu) + TINY;

  // Calculate some quantities we are going to use a lot  
  // Calculate pieces of J  
  JJ0     = W2*(E - 2.0*vdotF); 
  JJthin  = W2*E*vdotfh*vdotfh;
  JJthick = (W2-1.0)/(1.0+2.0*W2)*(4.0*W2*vdotF + (3.0-2.0*W2)*E);

  // Calculate pieces of H^alpha H_alpha 
  const double cn = W*JJ0 + W*(vdotF - E);
  const double cv = W*JJ0;
  const double cF =-W;
         
  const double Cthickn = W*JJthick;
  const double Cthickv = W*JJthick + W/(2.0*W2+1.0)*((3.0-2.0*W2)*E 
      + (2.0*W2-1.0)*vdotF);
  const double CthickF = W*v2;    
   
  const double Cthinn  = W*JJthin;
  const double Cthinv  = Cthinn;
  const double Cthinfh = W*E*vdotfh;
  
  if (Force_closure) {
    xi = xi0;
  } else if (v2>1.e-15) { // Use the finite velocity closure 
    Eclosure = E; 
    HH0 = cv*cv*v2 + cF*cF*F2 + 2.0*cv*cF*vdotF - cn*cn;
    HHthickthick = Cthickv*Cthickv*v2 + CthickF*CthickF*F2 
        + 2.0*CthickF*Cthickv*vdotF - Cthickn*Cthickn;
    HHthinthin   = Cthinv* Cthinv* v2 + Cthinfh*Cthinfh    
        + 2.0*Cthinfh*Cthinv*vdotfh - Cthinn*Cthinn ;
    HHthin       = 2.0*(cv*Cthinv* v2 + cF*Cthinfh*Fdotfh 
        + Cthinfh*cv*vdotfh + Cthinv*cF*vdotF - Cthinn*cn); 
    HHthick      = 2.0*(cv*Cthickv*v2 + cF*CthickF*F2 
        + CthickF*cv*vdotF + Cthickv*cF*vdotF - Cthickn*cn); 
    HHthickthin  = 2.0*(Cthinv*Cthickv*v2 
        + Cthinfh*CthickF*Fdotfh + Cthinfh*Cthickv*vdotfh 
        + Cthinv*CthickF*vdotF - Cthinn* Cthickn);  
    
    // Find xi that satisfies our non-linear equation using NR iteration
    // Compare the zero velocity limit and the optically thick limit for xi
    // and take the larger value for the starting guess
    if (xi0<1.e-5 || xi0>1.0) {
      xi = sqrt(F2/(E*E + TINY));
      if (F2<=TINY) xi = 0.0;
      xi = min(xi,1.0);
    } else {
      xi = xi0;
    }
          
    double dfdxi = 1.0;
    double f     = 1.0;
    int i=0;
    int NRITER = 10;
          
    double xil = xi*0.95;
    double xiu = min(1.0,xi*1.05);
    double fl = zFunc(xil);
    double fu = zFunc(xiu); 

    for (i=0;i<=NRITER;i++) {
      f = zFuncWithD(xi,&dfdxi);

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
    if (i>=NRITER){ // NR Failed, so do some good old bisection 
      fl = zFunc(0.0);
      fu = zFunc(1.0);
      xil = 0.0;
      xiu = 1.0;
      
      if ((xi<0.0) || (xi>1.0)) xi = 0.5;
      
      if (fl*fu>0.0) { // Doesn't go to zero in interval
        // Forge ahead wildly
        if (abs(fl)<abs(fu) && abs(fl)<2.e-1) xi = xil;
        else if (abs(fu)<abs(fl) && abs(fu)<2.e-1) xi = xiu;
        else xi = 1.0;
      } else {  
        for (int j=0;j<100;j++) {
          if (abs(xiu-xil)<1.e-3) break;
          f = zFunc(xi);
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
  const double chi = closureFunc(xi);
  const double dthick = 1.5 - 1.5*chi;
  const double dthin  = 1.5*chi - 0.5;
        
  // Rest frame energy density and projected fluxes in the optically thick
  // limit.
  const double Jthick = 3.0/(2.0*W2 + 1.0) * ((2.0*W2 - 1.0)*E - 2.0*W2*vdotF);
   
  const double tHxt = Fxu/W + vx*W/(2.0*W2 + 1.0) 
      * ((4.0*W2 + 1.0)*vdotF - 4.0*W2*E);
  const double tHyt = Fyu/W + vy*W/(2.0*W2 + 1.0) 
      * ((4.0*W2 + 1.0)*vdotF - 4.0*W2*E);
  const double tHzt = Fzu/W + vz*W/(2.0*W2 + 1.0) 
      * ((4.0*W2 + 1.0)*vdotF - 4.0*W2*E);
  
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
        
//  if (isnan(Pxx)) {
//#pragma omp critical
//    cerr << "Pxx is a NaN. " << debout << " " << dthin << " " 
//         << chi << " " << xi << " " << Fx << " " << E << " " << v2 << " " 
//         << vx << " " << vxl << " " << Carpet::reflevel << "\n" ;
//    //    if(Carpet::reflevel > 3) CCTK_WARN(0,"aborting!");
//    //    CCTK_WARN(0,"aborting!");
//  } 
  
  // Return the closure factor if anyone is interested
  return xi; 

}

