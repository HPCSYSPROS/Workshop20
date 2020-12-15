#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

#include "loopcontrol.h"
#include "HydroBase.h"

#include <math.h>

// given a (theta, phi), return minimum radius of spherical shape of all neighboring points on sphere
// This algorithm assumes that the spherical surface encloses a star-shaped region.
CCTK_REAL SetMask_MinRad(const cGH* cctkGH, const CCTK_REAL theta, const CCTK_REAL phi, const int si)
{
   DECLARE_CCTK_PARAMETERS;
   DECLARE_CCTK_ARGUMENTS;

   // get neighboring points
   const int i = (theta-sf_origin_theta[si]) / sf_delta_theta[si];
   const int j = (phi-sf_origin_phi[si]) / sf_delta_phi[si];
   
   // handle poles (theta < sf_origin_theta && theta > sf_delta_theta * sf_ntheta)
   if (theta-sf_origin_theta[si] < 0) {
      CCTK_REAL min = 1e99;
      for (int j=0; j < sf_nphi[si]; ++j) {
         int const ij = 0 + maxntheta * (j + maxnphi * si);
         if (min > sf_radius[ij])
            min = sf_radius[ij];
      }
      return min;
   }
   if (i >= sf_ntheta[si]-1) {
      CCTK_REAL min = 1e99;
      for (int j=0; j < sf_nphi[si]; ++j) {
         int const ij = sf_ntheta[si]-1 + maxntheta * (j + maxnphi * si);
         if (min > sf_radius[ij])
            min = sf_radius[ij];
      }
      return min;
   }
   
   int const j0 = phi-sf_origin_phi[si] < 0 ? sf_nphi[si]-1 : j;
   int const j1 = j+1 >= sf_nphi[si]-1 ? 0 : j+1;
   
   int const ij00 = i + maxntheta * (j0 + maxnphi * si);
   int const ij01 = i + maxntheta * (j1 + maxnphi * si);
   int const ij10 = (i+1) + maxntheta * (j0 + maxnphi * si);
   int const ij11 = (i+1) + maxntheta * (j1 + maxnphi * si);
   
   const CCTK_REAL a = sf_radius[ij00] < sf_radius[ij01] ? sf_radius[ij00] : sf_radius[ij01];
   const CCTK_REAL b = sf_radius[ij10] < sf_radius[ij11] ? sf_radius[ij10] : sf_radius[ij11];

   return a < b ? a : b;
}


void SetMask_SphericalSurface (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT *mask = (CCTK_INT*) CCTK_VarDataPtr(cctkGH, 0, SetMask_MaskName);
  if (!mask)
    CCTK_WARN(0, "No such variable, or no storage enabled");

  /* Delete mask first!  */
  #pragma omp parallel
      {
        LC_LOOP3(mask_zero, i,j,k, 0,0,0,
                 cctk_lsh[0], cctk_lsh[1], cctk_lsh[2],
                 cctk_lsh[0], cctk_lsh[1], cctk_lsh[2])
        {
          CCTK_INT i3D = CCTK_GFINDEX3D(cctkGH, i, j, k);
          mask[i3D] = HYDRO_EXCISION_NORMAL;
        }
        LC_ENDLOOP3(mask_zero);
      }

  /* Now set excision! */
  for (int smi = 0; smi < 10; smi++)
  {
    CCTK_INT sfi = sf_IdFromName(SetMask_SurfaceIndex[smi], SetMask_SurfaceName[smi]);
    if (sfi >= 0 && sf_active[sfi] && sf_valid[sfi] >= 0)
    {
    
      if (!SetMask_TrueShape[smi]) {
      // classical crude algorithm that just looks at the minimum radius
      #pragma omp parallel
      {
        LC_LOOP3(setsurface, i,j,k, 0,0,0,
                 cctk_lsh[0], cctk_lsh[1], cctk_lsh[2],
                 cctk_lsh[0], cctk_lsh[1], cctk_lsh[2])
        {
          CCTK_INT i3D = CCTK_GFINDEX3D(cctkGH, i, j, k);
          CCTK_REAL dist2 = (sf_centroid_x[sfi]-x[i3D]) * (sf_centroid_x[sfi]-x[i3D]) +
                            (sf_centroid_y[sfi]-y[i3D]) * (sf_centroid_y[sfi]-y[i3D]) +
                            (sf_centroid_z[sfi]-z[i3D]) * (sf_centroid_z[sfi]-z[i3D]);
          if (dist2 < SetMask_RadiusFactor[smi] * sf_min_radius[sfi] *
                      SetMask_RadiusFactor[smi] * sf_min_radius[sfi])
          {
            mask[i3D] = HYDRO_EXCISION_EXCISED;
          }
        }
        LC_ENDLOOP3(setsurface);
      }
      
      } else {
      
      // Fast algorithm that approx. captures the true shape.
      // It assumes that the spherical surface encloses a star-shaped region.
      #pragma omp parallel
      {
        LC_LOOP3(setsurface, i,j,k, 0,0,0,
                 cctk_lsh[0], cctk_lsh[1], cctk_lsh[2],
                 cctk_lsh[0], cctk_lsh[1], cctk_lsh[2])
        {
          CCTK_INT i3D = CCTK_GFINDEX3D(cctkGH, i, j, k);
          CCTK_REAL dist2 = (sf_centroid_x[sfi]-x[i3D]) * (sf_centroid_x[sfi]-x[i3D]) +
                            (sf_centroid_y[sfi]-y[i3D]) * (sf_centroid_y[sfi]-y[i3D]) +
                            (sf_centroid_z[sfi]-z[i3D]) * (sf_centroid_z[sfi]-z[i3D]);
          
          // get (theta, phi) of current point
          const CCTK_REAL theta = acos(z[i3D] / r[i3D]);
          const CCTK_REAL tmp = atan2(y[i3D], x[i3D]);  // angle between (-pi,pi]
          const CCTK_REAL phi = fmod(tmp < 0 ? tmp + 2*M_PI : tmp, 2*M_PI);  // angle between [0, 2pi)
          
          // get min radius of neighboring (theta,phi) points on sphere
          const CCTK_REAL minRad = SetMask_MinRad(cctkGH, theta, phi, sfi);
          
          if (dist2 < SetMask_RadiusFactor[smi] * minRad *
                      SetMask_RadiusFactor[smi] * minRad)
          {
            mask[i3D] = HYDRO_EXCISION_EXCISED;
          }
        }
        LC_ENDLOOP3(setsurface);
      }
      
      
      }
      
    }
  }
}

