#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include <math.h>
#include <stdio.h>

#include "sYlm.h"

#define RHOCGSTOCACTUS (1.0e0/6.1755e17)

#define sx (&scon[0*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define sy (&scon[1*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define sz (&scon[2*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])

void CCCCGlobalModes_LocalStuff(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;


  int nx = cctk_lsh[0];
  int ny = cctk_lsh[1];
  int nz = cctk_lsh[2];

  int i,j,k;
  int index;

  if (*volume_form_state == 0)
  {
     CCTK_WARN(0, "Tell your 'Coordinates' implementation to activate storage for the volume form!");
  }
  if (*jacobian_state == 0)
  {
     CCTK_WARN(0, "Tell your 'Coordinates' implementation to activate storage for the Jacobian!");
  }

  CCTK_REAL localweight = 1.0e0;

  CCTK_REAL tiny,rt,rt2,rt4,x2,y2;

  CCTK_REAL cutrho,cutradius;
  CCTK_REAL cutrhoCoM,cutradiusCoM;

  cutrho = 0.0e0;
  cutradius = 1.0e10;

  if(CCTK_EQUALS(cut_method,"abs density")) {
    cutrho = density_cut*RHOCGSTOCACTUS;
    cutradius = 1.0e20;
  } else if (CCTK_EQUALS(cut_method,"rel density")) {
    cutrho = (*global_rho_max)/density_cut_factor*RHOCGSTOCACTUS;
    cutradius = 1.0e20;
  } else if (CCTK_EQUALS(cut_method,"radius")) {
    cutrho = 0.0e0;
    cutradius = radius_cut;
  }

  if(use_CoM_radius) {
    cutrhoCoM = 1.0e0*RHOCGSTOCACTUS;
    cutradiusCoM = CoM_radius;
  } else {
    cutrhoCoM = cutrho;
    cutradiusCoM = cutradius;
  }

  tiny = 1.e-50;

  // We don`t want to do anything until shortly before our start time...

  if (!(cctk_time >= (start_time*2.03e2-2.0e0))){
    return;
  }

  CCTK_REAL densfac = 1.0;

  if(do_saijo) {

#pragma omp parallel for private(k,j,i,index,localweight,x2,y2,rt,rt2,rt4,densfac)
    for (k=0;k<nz;k++)
      for (j=0;j<ny;j++)
	for (i=0;i<nx;i++) {

	  index = CCTK_GFINDEX3D(cctkGH,i,j,k);

	  localweight = 1.0e0 * volume_form[index];
	  
	  // transformation factor for cons. density which is a scalar tensor density of weight 1
	  densfac = fabs((  J11[index] * J22[index] * J33[index]
                          + J12[index] * J23[index] * J31[index]
                          + J13[index] * J21[index] * J32[index]
                          - J11[index] * J23[index] * J32[index]
                          - J12[index] * J21[index] * J33[index]
                          - J13[index] * J22[index] * J31[index]));

	  if(rho[index] <= cutrho || r[index] > cutradius) localweight = 0.0e0;

	  x2 = x[index]*x[index];
	  y2 = y[index]*y[index];


	  rt = sqrt(x2 + y2);
	  rt2 = (x2+y2);
	  rt4 = (x2+y2)*(x2+y2);


	  //exclude the origin. This gives a small systematic error.
	  if(rt > 1.0e-4) {

	    ddi_re[index] = localweight*densfac*dens[index]*x[index] / rt;
	    ddi_im[index] = localweight*densfac*dens[index]*y[index] / rt;

	    dquad_re[index] = localweight*densfac*dens[index]*(x2 - y2) / rt2;
	    dquad_im[index] = localweight*densfac*dens[index]*(2.0e0*x[index]*y[index]) / rt2;
	    
	    dsextu_re[index] = localweight*densfac*dens[index]*( x2*x2+y2*y2 - 6.0e0*x2*y2 ) /
	      rt4;

	    dsextu_im[index] = localweight*densfac*dens[index]*( 4.0e0*x[index]*x2 - 4.0e0*x[index]*y[index]*y2) / 
	      rt4;

	  } else {
	    
	    ddi_re[index]    = 0.0e0;
	    ddi_im[index]    = 0.0e0;
	    dquad_re[index]  = 0.0e0;
	    dquad_im[index]  = 0.0e0;
	    dsextu_re[index] = 0.0e0;
	    dsextu_im[index] = 0.0e0;
	    
	  }


	}

  }

  if(do_shibata) {
#pragma omp parallel for private(k,j,i,index,localweight,densfac)
    for (k=0;k<nz;k++)
      for (j=0;j<ny;j++)
	for (i=0;i<nx;i++) {

	  index = CCTK_GFINDEX3D(cctkGH,i,j,k);

	  localweight = 1.0e0 * volume_form[index];

	  if(rho[index] <= cutrho || r[index] > cutradius) localweight = 0.0e0;

          // transformation factor for cons. density which is a scalar tensor density of weight 1
	  densfac = fabs((  J11[index] * J22[index] * J33[index]
                          + J12[index] * J23[index] * J31[index]
                          + J13[index] * J21[index] * J32[index]
                          - J11[index] * J23[index] * J32[index]
                          - J12[index] * J21[index] * J33[index]
                          - J13[index] * J22[index] * J31[index]));

	  dIxx[index] = localweight*densfac*dens[index]*x[index]*x[index];
	  dIyy[index] = localweight*densfac*dens[index]*y[index]*y[index];
	  dIxy[index] = localweight*densfac*dens[index]*x[index]*y[index];

	}


  }

  //  cutradiusCoM=1.0e9;

  if(do_CoM) {
#pragma omp parallel for private(k,j,i,index,localweight,densfac)
    for (k=0;k<nz;k++)
      for (j=0;j<ny;j++)
	for (i=0;i<nx;i++) {

	  index = CCTK_GFINDEX3D(cctkGH,i,j,k);


	  localweight = 1.0e0 * volume_form[index];

	  if( (rho[index] <= cutrhoCoM) || (r[index] > cutradiusCoM) ) localweight = 0.0e0;

          // transformation factor for cons. density which is a scalar tensor density of weight 1
	  densfac = fabs((  J11[index] * J22[index] * J33[index]
                          + J12[index] * J23[index] * J31[index]
                          + J13[index] * J21[index] * J32[index]
                          - J11[index] * J23[index] * J32[index]
                          - J12[index] * J21[index] * J33[index]
                          - J13[index] * J22[index] * J31[index]));

	  dMass[index] = localweight*densfac*dens[index];
	  dMx[index] = localweight*densfac*dens[index]*x[index];
	  dMy[index] = localweight*densfac*dens[index]*y[index];
	  dMz[index] = localweight*densfac*dens[index]*z[index];

	}


  }


  if(do_P) {
#pragma omp parallel for private(k,j,i,index,localweight)
    for (k=0;k<nz;k++)
      for (j=0;j<ny;j++)
	for (i=0;i<nx;i++) {
	  index = CCTK_GFINDEX3D(cctkGH,i,j,k);

	  localweight = 1.0e0 * volume_form[index];

	  if(rho[index] <= cutrho || r[index] > cutradius) localweight = 0.0e0;
#warning "Need to ransform to global basis!"
	  dPx[index] = localweight*sx[index];
	  dPy[index] = localweight*sy[index];
	  dPz[index] = localweight*sz[index];

	}


  }
  
  

  
  if(do_qlm) {
    CCTK_REAL theta = 0.0;
    CCTK_REAL phi = 0.0;
    CCTK_REAL re = 0;
    CCTK_REAL im = 0;

#pragma omp parallel for private(k,j,i,index,localweight,re,im,theta,phi)
    for (k=0;k<nz;k++)
      for (j=0;j<ny;j++)
	for (i=0;i<nx;i++) {
	  index = CCTK_GFINDEX3D(cctkGH,i,j,k);

	  localweight = 1.0e0 * volume_form[index];

	  if(rho[index] <= cutrho || r[index] > cutradius) localweight = 0.0e0;

          if (r[index] > 0.0) {
             theta = acos(z[index]/r[index]);
             phi = atan2(y[index], x[index]);

             sYlm(0, qlm_l, qlm_m, theta, phi, &re, &im);
          }
          else
          {
             re = im = 0;
          }
          
	  qlm_integrand_re[index] = localweight * rho[index] * pow(r[index], qlm_l) * re;
          qlm_integrand_im[index] = localweight * rho[index] * pow(r[index], qlm_l) * im;
          

	}


  }


  return;
}
