#include <iostream>
#include <algorithm>
#include <time.h>
#include <gsl/gsl_sf_erf.h>
#include "ZelmaniM1.hh"
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"


using namespace std;


namespace ZelmaniM1 {

  extern "C" 
  void zm1_setup_tests(CCTK_ARGUMENTS) {

    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    // so far have only homogeneous sphere test coded

    if(CCTK_EQUALS(zm1_test_type,"Homogeneous Sphere")) {

      double dx = CCTK_DELTA_SPACE(1);
      int nx = cctk_lsh[0];
      int ny = cctk_lsh[1];
      int nz = cctk_lsh[2];
      
      double kappa = 4.0;
      double b = 0.8;
    for(int ig=0;ig<ngroups*nspecies;ig++)
      for(int k=zm1_ghost-2;k<nz-zm1_ghost+2;k++)
        for(int j=zm1_ghost-2;j<ny-zm1_ghost+2;j++)
          for(int i=zm1_ghost-2;i<nx-zm1_ghost+2;i++) {
	    
	    // get indices
	    int index3D   = CCTK_GFINDEX3D(cctkGH,i,j,k);
	    int index4D   = CCTK_VECTGFINDEX3D(cctkGH,i,j,k,ig);
	    int indopac4D = ig + ngroups*nspecies*CCTK_GFINDEX3D(cctkGH,i,j,k); 
            indopac4D = index4D; 
	    double fnu0 = 0.0; //1.0e-10; 
	    double rad   = sqrt(x[index3D]*x[index3D] + z[index3D]*z[index3D] + y[index3D]*y[index3D]);
	    double theta = acos(z[index3D]/rad);
	    if (rad>0.0) fnuz[index4D] = fnu0*cos(theta);
	    else fnuz[index4D] = 0.0;

            double rho   = sqrt(x[index3D]*x[index3D] + y[index3D]*y[index3D]);
	    double phi   = acos(x[index3D]/rho);
            if (rho>0.0) {
	    	fnux[index4D] = fnu0*sin(theta)*cos(phi);
	    	fnuy[index4D] = fnu0*sin(theta)*sin(phi);
	    } else {
	    	fnux[index4D] = 0.0;
	    	fnuy[index4D] = 0.0;
	    }   	    

	    enu[index4D]  = 1.e-10;
            
	    // Roll off the opacity by tanh
	    double srad = zm1_spheretest_radius; 
	    double fac = 0.5 + 0.5*tanh(1.0*(srad-r[index3D])/srad);
	    scat[indopac4D]   = 1.0e-40;
	    absorb[indopac4D] = fac*zm1_spheretest_kappa;
	    emis[indopac4D]   = fac*zm1_spheretest_b*zm1_spheretest_kappa;

	    //if(r[index3D] <= zm1_spheretest_radius) {
	    //    scat[indopac4D]   = 1.0e-40;
	    //    absorb[indopac4D] = zm1_spheretest_kappa;
	    //    emis[indopac4D]   = zm1_spheretest_b*zm1_spheretest_kappa;
	    //} // r < zm1_sphere_radius
	    //else {
	    //    scat[indopac4D]   = 1.0e-40;
	    //    absorb[indopac4D] = 1.0e-10*zm1_spheretest_kappa;
	    //    emis[indopac4D]   = 1.0e-10*zm1_spheretest_b*zm1_spheretest_kappa;
	    //} 
	  } // for...
    } else if(CCTK_EQUALS(zm1_test_type,"Diffusion Pulse")) {
      
      double dx = CCTK_DELTA_SPACE(1);
      int nx = cctk_lsh[0];
      int ny = cctk_lsh[1];
      int nz = cctk_lsh[2];
      
      for(int ig=0;ig<ngroups*nspecies;ig++)
        for(int k=zm1_ghost-2;k<nz-zm1_ghost+2;k++)
          for(int j=zm1_ghost-2;j<ny-zm1_ghost+2;j++)
            for(int i=zm1_ghost-2;i<nx-zm1_ghost+2;i++) {
	      
	      int index3D   = CCTK_GFINDEX3D(cctkGH,i,j,k);
	      int index4D   = CCTK_VECTGFINDEX3D(cctkGH,i,j,k,ig);
	      
	      scat[index4D]   = zm1_spheretest_kappa;
	      absorb[index4D] = 0.0;
	      emis[index4D]   = 0.0;

	      fnux[index4D]   = 0.0;
	      fnuy[index4D]   = 0.0;
	      fnuz[index4D]   = 0.0;
	       
	      double srad = zm1_spheretest_radius; 
	      enu[index4D] = exp(-r[index3D]*r[index3D]/(srad*srad));

      }

    } else if(CCTK_EQUALS(zm1_test_type,"Advect Radiation")) {

      int nx = cctk_lsh[0];
      int ny = cctk_lsh[1];
      int nz = cctk_lsh[2];
      
      double dx = CCTK_DELTA_SPACE(0);
    
    for(int ig=0;ig<ngroups*nspecies;ig++)
      for(int k=zm1_ghost-2;k<nz-zm1_ghost+2;k++)
        for(int j=zm1_ghost-2;j<ny-zm1_ghost+2;j++)
          for(int i=zm1_ghost-2;i<nx-zm1_ghost+2;i++) {
	    
	    // get indices
	    int index3D   = CCTK_GFINDEX3D(cctkGH,i,j,k);
	    int index4D   = CCTK_VECTGFINDEX3D(cctkGH,i,j,k,ig);
	    int indopac4D = ig + ngroups*nspecies*CCTK_GFINDEX3D(cctkGH,i,j,k); 
            indopac4D = index4D;

	    enu[index4D]  = 1.e-40; 
	    fnux[index4D] = 1.e-40; 
	    fnuy[index4D] = 1.e-40; 
	    fnuz[index4D] = 1.e-40; 

	    scat[indopac4D]   = 0.0;
	    absorb[indopac4D] = 0.0;
	    emis[indopac4D]   = 0.0;
	    
	    if(abs(x[index3D]) - 0.5*dx <= zm1_spheretest_radius) {
	      double xl =  max((x[index3D] - 0.5*dx)/zm1_spheretest_radius,-zm1_spheretest_radius);
	      double xu =  min((x[index3D] + 0.5*dx)/zm1_spheretest_radius, zm1_spheretest_radius);
	      double lo = gsl_sf_erf(xl);
	      double up = gsl_sf_erf(xu);
	      //enu[index4D]  = 0.8*exp(-x[index3D]*x[index3D]/pow(zm1_spheretest_radius,2));
	      enu[index4D]  = 0.88627*(up-lo)*zm1_spheretest_radius/dx;
	      fnux[index4D] = enu[index4D]; //enu[indopac4D]; 
	    } // x < zm1_sphere_radius
	  } // for...

    } else if(CCTK_EQUALS(zm1_test_type,"Absorbing Sphere")) {

      int nx = cctk_lsh[0];
      int ny = cctk_lsh[1];
      int nz = cctk_lsh[2];
      
      double kappa = 4.0;
      double b = 0.8;
    for(int ig=0;ig<ngroups*nspecies;ig++)
      for(int k=zm1_ghost-2;k<nz-zm1_ghost+2;k++)
        for(int j=zm1_ghost-2;j<ny-zm1_ghost+2;j++)
          for(int i=zm1_ghost-2;i<nx-zm1_ghost+2;i++) {
	    
	    // get indices
	    int index3D   = CCTK_GFINDEX3D(cctkGH,i,j,k);
	    int index4D   = CCTK_VECTGFINDEX3D(cctkGH,i,j,k,ig);
	    int indopac4D = ig + ngroups*nspecies*CCTK_GFINDEX3D(cctkGH,i,j,k); 
            indopac4D = index4D; 
	    enu[indopac4D]  = 1.e-10;

	    // Set up an absorbing sphere at the origin 
	    // and a planar source
	    if(r[index3D] <= zm1_spheretest_radius) {
		scat[indopac4D]   = 0.0;
		absorb[indopac4D] = zm1_spheretest_kappa;
		emis[indopac4D]   = 0.0;
	        // r < zm1_sphere_radius
	    } else if (x[index3D]<-3.5*zm1_spheretest_radius) {
		scat[indopac4D]   = 0.0;
		absorb[indopac4D] = zm1_spheretest_kappa;
		emis[indopac4D]   = zm1_spheretest_b*zm1_spheretest_kappa;
	    } else {
		scat[indopac4D]   = 0.0;
		absorb[indopac4D] = 0.0;
		emis[indopac4D]   = 0.0;
	    } 
	  } // for...
    } else {
      CCTK_WARN(0,"Test type not coded in yet!");
    }

    if(use_nuTmunu) {
      zm1_ComputeNuTmunu(CCTK_PASS_CTOC);
    }

    return;
  }

  extern "C" 
  void zm1_tests_opacity(CCTK_ARGUMENTS) {

    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    // so far have only homogeneous sphere test coded
    CCTK_INFO("Calculating test opacities.");

    if(CCTK_EQUALS(zm1_test_type,"Homogeneous Sphere")) {

      double dx = CCTK_DELTA_SPACE(1);
      int nx = cctk_lsh[0];
      int ny = cctk_lsh[1];
      int nz = cctk_lsh[2];
      
      double kappa = 4.0;
      double b = 0.8;
    for(int ig=0;ig<ngroups*nspecies;ig++)
      for(int k=zm1_ghost-2;k<nz-zm1_ghost+2;k++)
        for(int j=zm1_ghost-2;j<ny-zm1_ghost+2;j++)
          for(int i=zm1_ghost-2;i<nx-zm1_ghost+2;i++) {
	    
	    // get indices
	    int index3D   = CCTK_GFINDEX3D(cctkGH,i,j,k);
	    int index4D   = CCTK_VECTGFINDEX3D(cctkGH,i,j,k,ig);
	    int indopac4D = ig + ngroups*nspecies*CCTK_GFINDEX3D(cctkGH,i,j,k); 
            indopac4D = index4D; 
	    
	    double rad   = sqrt(x[index3D]*x[index3D] + z[index3D]*z[index3D] + y[index3D]*y[index3D]);
	    double theta = acos(z[index3D]/rad);

            double rho   = sqrt(x[index3D]*x[index3D] + y[index3D]*y[index3D]);
	    double phi   = acos(x[index3D]/rho);
            
	    // Roll off the opacity by tanh
	    double srad = zm1_spheretest_radius; 
	    double fac = 0.5 + 0.5*tanh(1.0*(srad-r[index3D])/srad);
	    scat[indopac4D]   = 1.0e-40;
	    absorb[indopac4D] = fac*zm1_spheretest_kappa;
	    emis[indopac4D]   = fac*zm1_spheretest_b*zm1_spheretest_kappa;

	    //if(r[index3D] <= zm1_spheretest_radius) {
	    //    scat[indopac4D]   = 1.0e-40;
	    //    absorb[indopac4D] = zm1_spheretest_kappa;
	    //    emis[indopac4D]   = zm1_spheretest_b*zm1_spheretest_kappa;
	    //} // r < zm1_sphere_radius
	    //else {
	    //    scat[indopac4D]   = 1.0e-40;
	    //    absorb[indopac4D] = 1.0e-10*zm1_spheretest_kappa;
	    //    emis[indopac4D]   = 1.0e-10*zm1_spheretest_b*zm1_spheretest_kappa;
	    //} 
	  } // for...
    } else if(CCTK_EQUALS(zm1_test_type,"Diffusion Pulse")) {
      
      double dx = CCTK_DELTA_SPACE(1);
      int nx = cctk_lsh[0];
      int ny = cctk_lsh[1];
      int nz = cctk_lsh[2];
      
      for(int ig=0;ig<ngroups*nspecies;ig++)
        for(int k=zm1_ghost-2;k<nz-zm1_ghost+2;k++)
          for(int j=zm1_ghost-2;j<ny-zm1_ghost+2;j++)
            for(int i=zm1_ghost-2;i<nx-zm1_ghost+2;i++) {
	      
	      int index3D   = CCTK_GFINDEX3D(cctkGH,i,j,k);
	      int index4D   = CCTK_VECTGFINDEX3D(cctkGH,i,j,k,ig);
	      
	      scat[index4D]   = zm1_spheretest_kappa;
	      absorb[index4D] = 0.0;
	      emis[index4D]   = 0.0;

      }

    } else if(CCTK_EQUALS(zm1_test_type,"Advect Radiation")) {

      int nx = cctk_lsh[0];
      int ny = cctk_lsh[1];
      int nz = cctk_lsh[2];
      
      double dx = CCTK_DELTA_SPACE(0);
    
    for(int ig=0;ig<ngroups*nspecies;ig++)
      for(int k=zm1_ghost-2;k<nz-zm1_ghost+2;k++)
        for(int j=zm1_ghost-2;j<ny-zm1_ghost+2;j++)
          for(int i=zm1_ghost-2;i<nx-zm1_ghost+2;i++) {
	    
	    // get indices
	    int index3D   = CCTK_GFINDEX3D(cctkGH,i,j,k);
	    int index4D   = CCTK_VECTGFINDEX3D(cctkGH,i,j,k,ig);
	    int indopac4D = ig + ngroups*nspecies*CCTK_GFINDEX3D(cctkGH,i,j,k); 
            indopac4D = index4D;

	    scat[indopac4D]   = 0.0;
	    absorb[indopac4D] = 0.0;
	    emis[indopac4D]   = 0.0;
	    
	  } // for...

    } else if(CCTK_EQUALS(zm1_test_type,"Absorbing Sphere")) {

      int nx = cctk_lsh[0];
      int ny = cctk_lsh[1];
      int nz = cctk_lsh[2];
      
      double kappa = 4.0;
      double b = 0.8;
    for(int ig=0;ig<ngroups*nspecies;ig++)
      for(int k=zm1_ghost-2;k<nz-zm1_ghost+2;k++)
        for(int j=zm1_ghost-2;j<ny-zm1_ghost+2;j++)
          for(int i=zm1_ghost-2;i<nx-zm1_ghost+2;i++) {
	    
	    // get indices
	    int index3D   = CCTK_GFINDEX3D(cctkGH,i,j,k);
	    int index4D   = CCTK_VECTGFINDEX3D(cctkGH,i,j,k,ig);
	    int indopac4D = ig + ngroups*nspecies*CCTK_GFINDEX3D(cctkGH,i,j,k); 
            indopac4D = index4D; 

	    // Set up an absorbing sphere at the origin 
	    // and a planar source
	    if(r[index3D] <= zm1_spheretest_radius) {
		scat[indopac4D]   = 0.0;
		absorb[indopac4D] = zm1_spheretest_kappa;
		emis[indopac4D]   = 0.0;
	        // r < zm1_sphere_radius
	    } else if (x[index3D]<-3.5*zm1_spheretest_radius) {
		scat[indopac4D]   = 0.0;
		absorb[indopac4D] = zm1_spheretest_kappa;
		emis[indopac4D]   = zm1_spheretest_b*zm1_spheretest_kappa;
	    } else {
		scat[indopac4D]   = 0.0;
		absorb[indopac4D] = 0.0;
		emis[indopac4D]   = 0.0;
	    } 
	  } // for...
    }
    return;
  }
}

