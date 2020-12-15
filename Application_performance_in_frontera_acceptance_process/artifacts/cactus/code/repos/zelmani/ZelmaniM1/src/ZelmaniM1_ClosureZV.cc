#include <cmath>
#include <math.h>
#include <iostream>
#include "ZelmaniM1_Closure.hh"
// These need to go to make this general
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

double ClosureZV::setClosure() {
		
	// Given the rest frame E and F_i, as well as the local three 
	// velocity v^i and the three metric gamma_{ij}, this routine 
	// solves the closure relation for the radiation stress tensor P^{ij}
	// in the lab frame and returns the Eddington factor and components 
	// of the stress tensor. 
	
	// Calculate some quantities we are going to use a lot	
	double TINY = 1.e-45;
  
	xi = sqrt(F2/(E*E + TINY));
	if (xi>1.0) xi = 1.0;
	if (xi<0.0) xi = 0.0;	
        
	// Fill the stress tensor 
	double chi = closureFunc(xi);
		
	double dthick = 1.5 - 1.5*chi;
	double dthin  = 1.5*chi - 0.5;

	Pxx = dthick*E/3.0*gmet.guxx;
	Pxy = dthick*E/3.0*gmet.guxy;
	Pxz = dthick*E/3.0*gmet.guxz;
	Pyy = dthick*E/3.0*gmet.guyy;
	Pyz = dthick*E/3.0*gmet.guyz;
	Pzz = dthick*E/3.0*gmet.guzz;
	
	Pxx = Pxx + dthin*E*Fxu*Fxu/F2;
	Pxy = Pxy + dthin*E*Fxu*Fyu/F2;
	Pxz = Pxz + dthin*E*Fxu*Fzu/F2;
	Pyy = Pyy + dthin*E*Fyu*Fyu/F2;
	Pyz = Pyz + dthin*E*Fyu*Fzu/F2;
	Pzz = Pzz + dthin*E*Fzu*Fzu/F2;

	if (isnan(Pxx)) cout << "PxxZV is a NaN. " << xi << " " 
		<< Fy << " " << Fx << " " << Fz << " " << E << "\n" ;
				
	// Return the closure factor if anyone is interested
	return xi; 

}


