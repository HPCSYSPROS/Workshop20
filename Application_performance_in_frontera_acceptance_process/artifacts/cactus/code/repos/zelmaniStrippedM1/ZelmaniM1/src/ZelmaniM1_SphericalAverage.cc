#include <iostream>
#include <algorithm>
#include <time.h>
#include <math.h>
#include <vector>
#include <cstring>
#include "ZelmaniM1.hh"
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "carpet.hh"

using namespace std;


namespace ZelmaniM1 {

  const int navrad = 500;

  // global mode
  extern "C"
  void zm1_Spherical_Average_Init(CCTK_ARGUMENTS)
  {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    assert(navrad > 0);
    assert(avradmax > 0.0);
    const double dr = avradmax / (navrad - 1.0);

    for (int n=0; n<nrad; ++n) {
      // zm1_av_rad is the zone outer radius
      zm1_av_rad[n]  = (n+1)*dr;
      zm1_av_mass[n] = 0.0;
      zm1_av_vol[n]  = 0.0;
      zm1_av_rho[n]  = 0.0;
      zm1_av_temp[n] = 0.0;
      zm1_av_ye[n]   = 0.0;
    }

    return;
  }

  static inline CCTK_REAL det(CCTK_REAL const& gxx,
			      CCTK_REAL const& gxy,
			      CCTK_REAL const& gxz,
			      CCTK_REAL const& gyy,
			      CCTK_REAL const& gyz,
			      CCTK_REAL const& gzz)
  {
    return
      -gxz*gxz*gyy + 2*gxy*gxz*gyz - gxx*gyz*gyz - gxy*gxy*gzz + gxx*gyy*gzz;
  }



// local mode
  extern "C" 
  void zm1_SphericalAverage(CCTK_ARGUMENTS) {

    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;
    
    if ( (cctk_iteration%update_av_closure_every !=0)  ) return;

    if (zm1_verbose)
      CCTK_Info(CCTK_THORNSTRING, "Doing averaged background!");

    // Grid cell volume and radial spacing
    const double dr = avradmax / (navrad - 1.0);
    CCTK_REAL dV = 1.0;
    for (int d=0; d<3; ++d) {
      dV *= CCTK_DELTA_SPACE(d);
    }

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
	  if (w < 1.0e-10) continue;

	  // Radius, and index into 1d array
	  CCTK_REAL const rL = r[ind3d];
	  int const n = rL / dr;
	  if (n >= navrad) continue; // ignore far away grid points   

	  // Calculate volume in current cell
	  CCTK_REAL const detg = det(gxx[ind3d], gxy[ind3d], gxz[ind3d],
				     gyy[ind3d], gyz[ind3d], gzz[ind3d]);
	  CCTK_REAL const sqrt_detg = sqrt(detg);
	  CCTK_REAL const vol = sqrt_detg * w * dV;

	  // add to 1d arrays:
	  zm1_av_vol[n]  += vol;
	  zm1_av_mass[n] += rho[ind3d]*vol;
	  zm1_av_temp[n] += temperature[ind3d]*vol;
	  zm1_av_ye[n]   += Y_e[ind3d]*zm1_av_mass[n];
	  
	}
      }
    }
  } // zm1_SphericalAverage

  // global mode
  void zm1_SphericalAverage_combine(CCTK_ARGUMENTS)
  {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;
    
    if ((cctk_iteration%update_av_closure_every!=0)) return;

    if (zm1_verbose)
      CCTK_Info(CCTK_THORNSTRING, "Summing arrays");

    int const sum = CCTK_ReductionArrayHandle("sum");
    if (sum<0) {
      CCTK_WARN(CCTK_WARN_ABORT, "'sum' reduction handle not defined");
    }

    vector<CCTK_REAL> tmp(navrad);
    
    {
    int const ierr =
      CCTK_ReduceLocArrayToArray1D(cctkGH, -1, sum,
                                   zm1_av_mass, &tmp[0], navrad, CCTK_VARIABLE_REAL);
    if (ierr<0) {
      CCTK_WARN(CCTK_WARN_ABORT, "Could not reduce zm1_av_rho");
    }
    memcpy(zm1_av_mass, &tmp[0], nrad*sizeof *zm1_av_mass);
    }


    {
    int const ierr =
      CCTK_ReduceLocArrayToArray1D(cctkGH, -1, sum,
                                   zm1_av_temp, &tmp[0], navrad, CCTK_VARIABLE_REAL);
    if (ierr<0) {
      CCTK_WARN(CCTK_WARN_ABORT, "Could not reduce zm1_av_rho");
    }
    memcpy(zm1_av_temp, &tmp[0], nrad*sizeof *zm1_av_temp);
    }

    {
    int const ierr =
      CCTK_ReduceLocArrayToArray1D(cctkGH, -1, sum,
                                   zm1_av_ye, &tmp[0], navrad, CCTK_VARIABLE_REAL);
    if (ierr<0) {
      CCTK_WARN(CCTK_WARN_ABORT, "Could not reduce zm1_av_rho");
    }
    memcpy(zm1_av_ye, &tmp[0], nrad*sizeof *zm1_av_ye);
    }

    {
    int const ierr =
      CCTK_ReduceLocArrayToArray1D(cctkGH, -1, sum,
                                   zm1_av_vol, &tmp[0], navrad, CCTK_VARIABLE_REAL);
    if (ierr<0) {
      CCTK_WARN(CCTK_WARN_ABORT, "Could not reduce zm1_av_rho");
    }
    memcpy(zm1_av_vol, &tmp[0], nrad*sizeof *zm1_av_vol);
    }

  } // zm1_SphericalAverage_combine

  // global mode                                                                                                    
  void zm1_SphericalAverage_transform(CCTK_ARGUMENTS)
  {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    if ((cctk_iteration%update_av_closure_every!=0)) return;

    if (zm1_verbose)
      CCTK_Info(CCTK_THORNSTRING, "Transforming Arrays");

    for (int i=0; i<nrad; ++i) {

      double mas = zm1_av_mass[i];
      if (mas>0.0) {
	zm1_av_ye[i]  = zm1_av_ye[i]/mas;
      }else{
	zm1_av_ye[i]  = 0.0;
      }

      double vol = zm1_av_vol[i];
      if (vol>0.0) {
	zm1_av_rho[i]       = mas/vol;
	zm1_av_temp[i]      = zm1_av_temp[i]/vol;
      }else{
	zm1_av_rho[i]       = mas/vol;
	zm1_av_temp[i]      = zm1_av_temp[i]/vol;

      }
    }

  }

} // namespace ZelmaniM1

