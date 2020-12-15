#ifndef ZelmaniM1_HH
#define ZelmaniM1_HH

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

#ifdef DEFINE_GLOBALS
#define EXTERN
#else
#define EXTERN extern
#endif

#define RHO_GF 1.61620075314614e-18
#define INV_RHO_GF 6.18735016707159e17

namespace ZelmaniM1 {

  // Global opacity table storage
  EXTERN int ntables;
  EXTERN int nrho,ntemp,nye,ngtab;
  EXTERN int bound_table;
  EXTERN double  *rho_points, *temp_points, *ye_points;
  EXTERN double  *neutrino_energies;
  EXTERN double  *bin_bottom, *bin_top, *bin_widths;
  EXTERN double  *alltables;

  void linterp_many(double x, double y, double z,
		    double* f, double* ft, 
		    int nx, int ny, int nz, int nvars,
		    double* xt, double* yt, double* zt);

  extern "C"
  void zm1_RegisterVars(CCTK_ARGUMENTS);

  extern "C"
  void zm1_SetSymmetries(CCTK_ARGUMENTS);

  extern "C"
  void zm1_Boundaries(CCTK_ARGUMENTS);

  extern "C"
  void zm1_readtable(CCTK_ARGUMENTS);

  extern "C"
  void zm1_Initialize(CCTK_ARGUMENTS);

  extern "C"
  void zm1_getopacs(CCTK_ARGUMENTS);
  
  extern "C"
  void zm1_getopacsi_gray(CCTK_ARGUMENTS);

  extern "C"
  void zm1_UpdateEoS(CCTK_ARGUMENTS);
  
  extern "C"
  void zm1_UpdateTmunu(CCTK_ARGUMENTS);
  
  extern "C"
  void zm1_StartRKLoop(CCTK_ARGUMENTS);
  
  extern "C"
  void zm1_AdvanceRKLoop(CCTK_ARGUMENTS);

  extern "C"
  void zm1_StartLoop(CCTK_ARGUMENTS);

  extern "C"
  void zm1_AdvanceLoop(CCTK_ARGUMENTS);

  extern "C"
  void zm1_Reconstruct(CCTK_ARGUMENTS);

  extern "C"
  void zm1_RiemannHLLE(CCTK_ARGUMENTS);
  
  extern "C"
  void zm1_Redshift(CCTK_ARGUMENTS);
  
  extern "C"
  void zm1_FluxesDG(CCTK_ARGUMENTS);

  extern "C"
  void zm1_CalcUpdate(CCTK_ARGUMENTS);

  extern "C"
  void zm1_setup_tests(CCTK_ARGUMENTS);
  
  extern "C"
  void zm1_luminosity_output(CCTK_ARGUMENTS);
  
  // Approximate relativistic fermi integrals a la Takahashi 1978
  double inline Fermi_1(double eta){
  	if (eta>1.e-3){
  		return (eta*eta/2.0 + 1.6449)/(1.0+exp(-1.6855*eta));
  	}else{
  		return exp(eta)/(1.0 + 0.2159*exp(0.8857*eta));
  	}
  }
  
  double inline Fermi_2(double eta){
  	if (eta>1.e-3){
  		return (pow(eta,3)/3.0 + 3.2899*eta)/(1.0-exp(-1.8246*eta));
  	}else{
  		return 2.0*exp(eta)/(1.0 + 0.1092*exp(0.8908*eta));
  	}
  }
  
  double inline Fermi_3(double eta){
  	if (eta>1.e-3){
  		return (pow(eta,4)/4.0 + 4.9348*eta*eta + 11.3644)/(1.0+exp(-1.9039*eta));
  	}else{
  		return 6.0*exp(eta)/(1.0 + 0.0559*exp(0.9069*eta));
  	}
  }
  
  double inline Fermi_4(double eta){
  	if (eta>1.e-3){
  		return (pow(eta,5)/5.0 + 6.5797*pow(eta,3) + 45.4576*eta)/(1.0-exp(-1.9484*eta));
  	}else{
  		return 24.0*exp(eta)/(1.0 + 0.0287*exp(0.9257*eta));
  	}
  }
  
  double inline Fermi_5(double eta){
  	if (eta>1.e-3){
  		return (pow(eta,6)/6.0 + 8.2247*pow(eta,4) + 113.6439*eta*eta + 236.5323)/(1.0+exp(-1.9727*eta));
  	}else{
  		return 120.0*exp(eta)/(1.0 + 0.0147*exp(0.9431*eta));
  	}
  }

}

#endif
