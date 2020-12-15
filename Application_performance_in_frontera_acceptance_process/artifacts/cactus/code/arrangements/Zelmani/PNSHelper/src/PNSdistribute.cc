#include <cassert>
#include <cmath>
#include <cstring>
#include <vector>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

using namespace std;



extern "C" {
  void PNSdistribute_init(CCTK_ARGUMENTS);
  void PNSdistribute_local(CCTK_ARGUMENTS);
  void PNSdistribute_combine(CCTK_ARGUMENTS);
  void PNSdistribute_transform(CCTK_ARGUMENTS);
}



namespace {
  CCTK_REAL det(CCTK_REAL const& gxx,
                CCTK_REAL const& gxy,
                CCTK_REAL const& gxz,
                CCTK_REAL const& gyy,
                CCTK_REAL const& gyz,
                CCTK_REAL const& gzz)
  {
    return
      -gxz*gxz*gyy + 2*gxy*gxz*gyz - gxx*gyz*gyz - gxy*gxy*gzz + gxx*gyy*gzz;
  }
}



// global mode
void PNSdistribute_init(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
   
  for (int n=0; n<nrad; ++n) {
    pns_vol[n]          = 0.0; 
    pns_rho_psi[n]      = 0.0; 
    pns_rho_alp[n]      = 0.0; 
    pns_vol[n]          = 0.0; 
    pns_mass[n]         = 0.0; 
    pns_av_rho[n]       = 0.0; 
    pns_av_eps[n]       = 0.0; 
    pns_av_press[n]     = 0.0; 
    pns_av_temp[n]      = 0.0; 
    pns_av_ye[n]        = 0.0; 
    pns_av_w_lorentz[n] = 0.0; 
    pns_av_gxx[n]       = 0.0; 
    pns_av_alp[n]       = 0.0; 
  }
}



// local mode
void PNSdistribute_local(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if ((cctk_iteration%update_GR_every!=0) && cctk_iteration > 16) return;
  
  if (verbose)  
    CCTK_Info(CCTK_THORNSTRING, "Doing distributed sum arrays");
  
  // Grid cell volume
  CCTK_REAL dV = 1.0;
  for (int d=0; d<3; ++d) {
    dV *= CCTK_DELTA_SPACE(d);
  }
  
  // Grid spacing of 1d array
  //CCTK_REAL const dr = rad_max / (nrad-1);
  CCTK_REAL const dr = *drad;
  
  // Loop bounds
  int imin[3], imax[3];
  for (int d=0; d<3; ++d) {
    imin[d] = cctk_nghostzones[d];
    imax[d] = cctk_lsh[d] - cctk_nghostzones[d];
  }
  
  // Weight function
  CCTK_REAL const *restrict const weight =
    static_cast<CCTK_REAL const*>
    (CCTK_VarDataPtr(cctkGH, 0, "CarpetReduce::weight"));
  if (not weight) {
    CCTK_WARN(CCTK_WARN_ABORT,
              "Grid function 'CarpetReduce::weight' does not have storage");
  }
  
  // This loop is only parallel if the reduction operations are
  // declared correctly. Alternatively, we could allocate one 1d array
  // per thread.
//#pragma omp parallel for 
  for (int k=imin[2]; k<imax[2]; ++k) {
    for (int j=imin[1]; j<imax[1]; ++j) {
      for (int i=imin[0]; i<imax[0]; ++i) {
        int const ind3d = CCTK_GFINDEX3D(cctkGH, i,j,k);
        CCTK_REAL const w = weight[ind3d];
        if (w == 0.0) continue;
        
        // Radius, and index into 1d array
        CCTK_REAL const rL = r[ind3d];
        int const n = rL / dr;
        if (n >= nrad) continue; // ignore far away grid points
        
        // Calculate volume in current cell
        CCTK_REAL const detg = det(gxx[ind3d], gxy[ind3d], gxz[ind3d],
                                   gyy[ind3d], gyz[ind3d], gzz[ind3d]);
        CCTK_REAL const sqrt_detg = sqrt(detg);
        CCTK_REAL const vol = sqrt_detg * w * dV;
        
	double rhos,S;
	
        // Calculate sources using HydroBase quantities
        //double mass = rho[ind3d]         * vol;
        //double ee   = eps[ind3d]         * mass;
        //double pp   = press[ind3d]       * vol;
	//double W    = w_lorentz[ind3d];
        //rhos = W*W*(mass+ee+pp) - pp;
	//S   = (W*W-1.0)*(mass+ee+pp) + 3.0*pp;
        
	// Calculate sources using TmunuBase
	double guxx = (gyy[ind3d]*gzz[ind3d] - gyz[ind3d]*gyz[ind3d])/detg;
	double guxy = (gxz[ind3d]*gyz[ind3d] - gxy[ind3d]*gzz[ind3d])/detg;
	double guxz = (gxy[ind3d]*gyz[ind3d] - gxz[ind3d]*gyy[ind3d])/detg;
	double guyy = (gxx[ind3d]*gzz[ind3d] - gxz[ind3d]*gxz[ind3d])/detg;
	double guyz = (gxy[ind3d]*gxz[ind3d] - gxx[ind3d]*gyz[ind3d])/detg;
	double guzz = (gxx[ind3d]*gyy[ind3d] - gxy[ind3d]*gxy[ind3d])/detg;
        
	rhos = (     eTtt[ind3d] 
	       - 2.0*eTtx[ind3d]*betax[ind3d] 
	       - 2.0*eTty[ind3d]*betay[ind3d] 
	       - 2.0*eTtz[ind3d]*betaz[ind3d] 
	       +     eTxx[ind3d]*betax[ind3d]*betax[ind3d] 
	       +     eTyy[ind3d]*betay[ind3d]*betay[ind3d] 
	       +     eTzz[ind3d]*betaz[ind3d]*betaz[ind3d] 
	       + 2.0*eTxy[ind3d]*betax[ind3d]*betay[ind3d] 
	       + 2.0*eTxz[ind3d]*betax[ind3d]*betaz[ind3d] 
	       + 2.0*eTyz[ind3d]*betay[ind3d]*betaz[ind3d]
	       )/(alp[ind3d]*alp[ind3d]);
	
	rhos = rhos * vol;      
	
	S =  guxx*eTxx[ind3d] + guyy*eTyy[ind3d] + guzz*eTzz[ind3d] 
	   + 2.0*(guxy*eTxy[ind3d] + guxz*eTxz[ind3d] + guyz*eTyz[ind3d]);
        
	S = S * vol;

        // Add to 1d arrays
        pns_vol[n]     += vol;
	pns_rho_psi[n] += rhos;
	pns_rho_alp[n] += rhos + 2.0*S;
        
	// Calculate mass in current cell
        CCTK_REAL const mass   = rho[ind3d]         * vol;
        CCTK_REAL const ee     = eps[ind3d]         * mass;
        CCTK_REAL const pp     = press[ind3d]       * vol;
        CCTK_REAL const tt     = temperature[ind3d] * vol;
        CCTK_REAL const ye     = Y_e[ind3d]         * mass;
        CCTK_REAL const ww     = w_lorentz[ind3d]   * vol;
        CCTK_REAL const alphao = alp[ind3d]         * vol;
        CCTK_REAL const gxxo   = gxx[ind3d]         * vol;
        
        // Add to 1d array
        pns_mass[n]     += mass;
	pns_av_eps[n]   += ee;
	pns_av_press[n] += pp;
	pns_av_temp[n]  += tt;
	pns_av_ye[n]    += ye;
	pns_av_w_lorentz[n] += ww;
	pns_av_gxx[n]   += gxxo;
	pns_av_alp[n]   += alphao;
      
      }
    }
  }
}



// global mode
void PNSdistribute_combine(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if ((cctk_iteration%update_GR_every!=0) && cctk_iteration > 16) return;
  
  if (verbose) 
    CCTK_Info(CCTK_THORNSTRING, "Summing arrays");
  
  int const sum = CCTK_ReductionArrayHandle("sum");
  if (sum<0) {
    CCTK_WARN(CCTK_WARN_ABORT, "'sum' reduction handle not defined");
  }
  
  vector<CCTK_REAL> tmp(nrad);
  
  // Sum the total volume in each radial bin 
  {
    int const ierr =
      CCTK_ReduceLocArrayToArray1D(cctkGH, -1, sum,
                                   pns_vol, &tmp[0], nrad, CCTK_VARIABLE_REAL);
    if (ierr<0) {
      CCTK_WARN(CCTK_WARN_ABORT, "Could not reduce pns_vol");
    }
    memcpy(pns_vol, &tmp[0], nrad*sizeof *pns_vol);
  }
  
  // Sum the total mass in each radial bin and get the average density 
  {
    int const ierr =
      CCTK_ReduceLocArrayToArray1D(cctkGH, -1, sum,
                                   pns_mass, &tmp[0], nrad, CCTK_VARIABLE_REAL);
    if (ierr<0) {
      CCTK_WARN(CCTK_WARN_ABORT, "Could not reduce pns_mass");
    }
    memcpy(pns_mass, &tmp[0], nrad*sizeof *pns_mass);
  }

  // Sum and average the total energy in each radial bin 
  {
    int const ierr =
      CCTK_ReduceLocArrayToArray1D(cctkGH, -1, sum,
                                   pns_av_eps, &tmp[0], nrad, CCTK_VARIABLE_REAL);
    if (ierr<0) {
      CCTK_WARN(CCTK_WARN_ABORT, "Could not reduce pns_av_eps");
    }
    memcpy(pns_av_eps, &tmp[0], nrad*sizeof *pns_av_eps);
  }

  // Sum the total pressure in each radial bin 
  {
    int const ierr =
      CCTK_ReduceLocArrayToArray1D(cctkGH, -1, sum,
                                   pns_av_press, &tmp[0], nrad, CCTK_VARIABLE_REAL);
    if (ierr<0) {
      CCTK_WARN(CCTK_WARN_ABORT, "Could not reduce pns_av_press");
    }
    memcpy(pns_av_press, &tmp[0], nrad*sizeof *pns_av_press);
  }

  // Sum the temperature in each radial bin 
  {
    int const ierr =
      CCTK_ReduceLocArrayToArray1D(cctkGH, -1, sum,
                                   pns_av_temp, &tmp[0], nrad, CCTK_VARIABLE_REAL);
    if (ierr<0) {
      CCTK_WARN(CCTK_WARN_ABORT, "Could not reduce pns_av_temp");
    }
    memcpy(pns_av_temp, &tmp[0], nrad*sizeof *pns_av_temp);
  }

  // Sum ye in each radial bin 
  {
    int const ierr =
      CCTK_ReduceLocArrayToArray1D(cctkGH, -1, sum,
                                   pns_av_ye, &tmp[0], nrad, CCTK_VARIABLE_REAL);
    if (ierr<0) {
      CCTK_WARN(CCTK_WARN_ABORT, "Could not reduce pns_av_ye");
    }
    memcpy(pns_av_ye, &tmp[0], nrad*sizeof *pns_av_ye);
  }

  // Sum the lorentz factor in each radial bin 
  {
    int const ierr =
      CCTK_ReduceLocArrayToArray1D(cctkGH, -1, sum,
                                   pns_av_w_lorentz, &tmp[0], nrad, CCTK_VARIABLE_REAL);
    if (ierr<0) {
      CCTK_WARN(CCTK_WARN_ABORT, "Could not reduce pns_av_w_lorentz");
    }
    memcpy(pns_av_w_lorentz, &tmp[0], nrad*sizeof *pns_av_w_lorentz);
  }

  // Sum the xx metric component in each radial bin 
  {
    int const ierr =
      CCTK_ReduceLocArrayToArray1D(cctkGH, -1, sum,
                                   pns_av_gxx, &tmp[0], nrad, CCTK_VARIABLE_REAL);
    if (ierr<0) {
      CCTK_WARN(CCTK_WARN_ABORT, "Could not reduce pns_av_gxx");
    }
    memcpy(pns_av_gxx, &tmp[0], nrad*sizeof *pns_av_gxx);
  }

  // Sum the lapse in each radial bin 
  {
    int const ierr =
      CCTK_ReduceLocArrayToArray1D(cctkGH, -1, sum,
                                   pns_av_alp, &tmp[0], nrad, CCTK_VARIABLE_REAL);
    if (ierr<0) {
      CCTK_WARN(CCTK_WARN_ABORT, "Could not reduce pns_av_alp");
    }
    memcpy(pns_av_alp, &tmp[0], nrad*sizeof *pns_av_alp);
  }
  
  // Sum the psi source in each radial bin 
  {
    int const ierr =
      CCTK_ReduceLocArrayToArray1D(cctkGH, -1, sum,
                                   pns_rho_psi, &tmp[0], nrad, CCTK_VARIABLE_REAL);
    if (ierr<0) {
      CCTK_WARN(CCTK_WARN_ABORT, "Could not reduce pns_rho_psi");
    }
    memcpy(pns_rho_psi, &tmp[0], nrad*sizeof *pns_rho_psi);
  }
  
  // Sum the alp psi source in each radial bin 
  {
    int const ierr =
      CCTK_ReduceLocArrayToArray1D(cctkGH, -1, sum,
                                   pns_rho_alp, &tmp[0], nrad, CCTK_VARIABLE_REAL);
    if (ierr<0) {
      CCTK_WARN(CCTK_WARN_ABORT, "Could not reduce pns_rho_alp");
    }
    memcpy(pns_rho_alp, &tmp[0], nrad*sizeof *pns_rho_alp);
  }
  
}

// global mode
void PNSdistribute_transform(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if ((cctk_iteration%update_GR_every!=0) && cctk_iteration > 16) return;
  
  
  if (verbose)  
  CCTK_Info(CCTK_THORNSTRING, "Transforming Arrays");
  
  for (int i=0; i<nrad; ++i) {
    
    double mas = pns_mass[i];
    if (mas>0.0) {
      pns_av_eps[i] = pns_av_eps[i]/mas;
      pns_av_ye[i]  = pns_av_ye[i]/mas;
    }else{
      pns_av_eps[i] = 0.0; 
      pns_av_ye[i]  = 0.0; 
    } 
    
    double vol = pns_vol[i];
    if (vol>0.0) {
      pns_av_rho[i]       = mas/vol;
      pns_av_press[i]     = pns_av_press[i]/vol;
      pns_av_temp[i]      = pns_av_temp[i]/vol;
      pns_av_w_lorentz[i] = pns_av_w_lorentz[i]/vol;
      pns_av_gxx[i]       = pns_av_gxx[i]/vol;
      pns_av_alp[i]       = pns_av_alp[i]/vol;
      
      pns_rho_psi[i]      = pns_rho_psi[i]/vol;
      pns_rho_alp[i]      = pns_rho_alp[i]/vol;
    }else{
      pns_rho_psi[i]      = 0.0;
      pns_rho_alp[i]      = 0.0;
      
      pns_av_rho[i]       = 0.0; 
      pns_av_press[i]     = 0.0;
      pns_av_temp[i]      = 0.0;
      pns_av_w_lorentz[i] = 1.0;
      pns_av_gxx[i]       = 1.0;
      pns_av_alp[i]       = 1.0;

    } 
  }
  
  *have_interp_data = 1;
  
}
