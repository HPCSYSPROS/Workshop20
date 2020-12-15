/* $Header$ */

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

#include "util_ErrorCodes.h"
#include "util_Table.h"

#define MAXDIM 3
#define REFLEVEL ((int)(0.1 + log10((CCTK_REAL)(cctk_levfac[0]))/log10(2.0)))

void
setup_epsdis (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int ai,aj,ak;
  int ni,nj,nk;
  int i,j,k,s,l,m;
  int index,indexP;
  int ierr;
  int npts;
  int *inds;
  CCTK_REAL *xa,*ya,*za,*rads;
  CCTK_REAL maxrad;
  CCTK_REAL xmin,xmax, ymin,ymax, zmin,zmax;
  CCTK_REAL odx,ody,odz;
  CCTK_REAL radp;
  const CCTK_INT  MAXSURFNUM=100; /* XXX hard limit */
  CCTK_INT doBC[2*MAXDIM],symbnd[2*MAXDIM];  
  CCTK_INT symtable;
  int reflvl = REFLEVEL;

  if (verbose) {
    CCTK_VInfo(CCTK_THORNSTRING,"Setting up spatially varying dissipation at T=%g",
               (double)cctk_time);
  }

  ai=cctk_ash[0];
  aj=cctk_ash[1];
  ak=cctk_ash[2];
  ni=cctk_lsh[0];
  nj=cctk_lsh[1];
  nk=cctk_lsh[2];

  if (epsdis_for_level[reflvl] > 0.0) {
    for (i=0; i<ai*aj*ak; ++i) {
      epsdisA[i] = epsdis_for_level[reflvl];
    }
  }
  else {
    for (i=0; i<ai*aj*ak; ++i) {
      epsdisA[i]=epsdis;
    }
  } 
  if (epsdis_for_level2[reflvl] > 0.0) {
      for (k=0;k<ak;k++) {
        for (j=0;j<aj;j++) {
          for (i=0;i<ai;i++) {
            index = CCTK_GFINDEX3D(cctkGH,i,j,k);
            CCTK_REAL radius = sqrt(x[index]*x[index] + y[index]*y[index] + z[index]*z[index]);
            CCTK_REAL omega = sqrt(x[index]*x[index] + y[index]*y[index]);
            if (radius > 60.0 && omega > 60.0) { 
              epsdisB[index]=epsdis_for_level2[reflvl];
            }
            else { 
              epsdisB[index] = 0.0;
            }
        }
      }
    }
  }
  else {
    for (i=0; i<ai*aj*ak; ++i) {
      epsdisB[i]=epsdis2;
    }
  } 
  if (extra_dissipation_at_outerbound) 
  {
    symtable = SymmetryTableHandleForGrid (cctkGH);
    if (symtable < 0) CCTK_WARN (1, "unable to get symmetry table");
    ierr=Util_TableGetIntArray(symtable, 6, symbnd, "symmetry_handle");
    if (ierr != 6) CCTK_WARN (1, "unable to get symmetry handle");
    for (i = 0; i < 6; i++) {
      doBC[i] = cctk_bbox[i]!=0 && symbnd[i] < 0;
    }

    if(doBC[0]) {
      for (k=0;k<nk;k++) {
        for (j=0;j<nj;j++) {
          for (i=0;i<outer_bound_npoints;i++) {
            index = CCTK_GFINDEX3D(cctkGH,i,j,k);
            indexP= CCTK_GFINDEX3D(cctkGH,outer_bound_npoints,j,k); 
            epsdisA[index]=epsdis+ob_slope*fabs(x[index]-x[indexP]);
	    if (epsdisA[index] > outer_boundary_max_epsdis) {
	      epsdisA[index] = outer_boundary_max_epsdis;
	    }
          }
        }
      }
    }
    if(doBC[1]) {
      for (k=0;k<nk;k++) {
        for (j=0;j<nj;j++) {
          for (i=ni-1;i>=ni-outer_bound_npoints;i--) {
            index = CCTK_GFINDEX3D(cctkGH,i,j,k);
            indexP= CCTK_GFINDEX3D(cctkGH,ni-outer_bound_npoints-1,j,k); 
	    epsdisA[index]=epsdis+ob_slope*fabs(x[index]-x[indexP]);
	    if (epsdisA[index] > outer_boundary_max_epsdis) {
	      epsdisA[index] = outer_boundary_max_epsdis;
	    }
          }
        }
      }
    }
    if(doBC[2]) {
      for (k=0;k<nk;k++) {
        for (j=0;j<outer_bound_npoints;j++) {
          for (i=0;i<ni;i++) {
            index = CCTK_GFINDEX3D(cctkGH,i,j,k);
            indexP= CCTK_GFINDEX3D(cctkGH,i,outer_bound_npoints,k); 
            epsdisA[index]=epsdis+ob_slope*fabs(y[index]-y[indexP]);
	    if (epsdisA[index] > outer_boundary_max_epsdis) {
	      epsdisA[index] = outer_boundary_max_epsdis;
	    }
          }
        }
      }
    }
    if(doBC[3]) {
      for (k=0;k<nk;k++) {
        for (j=nj-1;j>=nj-outer_bound_npoints;j--) {
          for (i=0;i<ni;i++) {
            index = CCTK_GFINDEX3D(cctkGH,i,j,k);
            indexP= CCTK_GFINDEX3D(cctkGH,i,nj-outer_bound_npoints-1,k); 
            epsdisA[index]=epsdis+ob_slope*fabs(y[index]-y[indexP]);
	    if (epsdisA[index] > outer_boundary_max_epsdis) {
	      epsdisA[index] = outer_boundary_max_epsdis;
	    }
	  }
	}
      }
    }

    if(doBC[4]) {
      for (k=0;k<outer_bound_npoints;k++) {
        for (j=0;j<nj;j++) {
          for (i=0;i<ni;i++) {
            index = CCTK_GFINDEX3D(cctkGH,i,j,k);
            indexP= CCTK_GFINDEX3D(cctkGH,i,j,outer_bound_npoints); 
            epsdisA[index]=epsdis+ob_slope*fabs(z[index]-z[indexP]);
	    if (epsdisA[index] > outer_boundary_max_epsdis) {
	      epsdisA[index] = outer_boundary_max_epsdis;
	    }
          }
        }
      }
    }
    if(doBC[5]) {
      for (k=nk-1;k>=nk-outer_bound_npoints;k--) {
        for (j=0;j<nj;j++) {
          for (i=0;i<ni;i++) {
            index = CCTK_GFINDEX3D(cctkGH,i,j,k);
            indexP= CCTK_GFINDEX3D(cctkGH,i,j,nk-outer_bound_npoints-1); 
            epsdisA[index]=epsdis+ob_slope*fabs(z[index]-z[indexP]);
	    if (epsdisA[index] > outer_boundary_max_epsdis) {
	      epsdisA[index] = outer_boundary_max_epsdis;
	    }
          }
        }
      }
    }
  }

  if (extra_dissipation_in_horizons && cctk_iteration%update_ah_every == 0)
  {
    if (verbose) {
      CCTK_INFO("Linear Interpolation into AH surfaces");
    }

    for (s=0;s<MAXSURFNUM;s++)
    {
      if (surface_number[s]==-1 && horizon_number[s]==-1) {
        continue;
      }
      
      if (surface_number[s]<0 || horizon_number[s]<=0) {
        CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Invalid specification for horizon %d", s);
        continue;
      }
      
      if (! sf_valid[surface_number[s]]) {
        if (verbose) {
          CCTK_VInfo(CCTK_THORNSTRING,
                     "Invalid Surface: s=%d, sf_va=%d, surf_no=%d, ah_no=%d",
                     s,(int)sf_valid[surface_number[s]],
                     (int)surface_number[s],(int)horizon_number[s]);
        }
        continue;
      }
      
      if (! CCTK_IsFunctionAliased ("HorizonRadiusInDirection"))
      {
        CCTK_WARN (1, "The aliased function \"HorizonRadiusInDirection\" must be defined when the parameter \"extra_dissipation_in_horizons\" is set and one of the sources is AHFinderDirect");
        continue;
      }

      maxrad=sf_max_radius[surface_number[s]];
      odx=sf_origin_x[surface_number[s]];
      ody=sf_origin_y[surface_number[s]];
      odz=sf_origin_z[surface_number[s]];
      xmin=odx-maxrad;
      xmax=odx+maxrad;
      ymin=ody-maxrad;
      ymax=ody+maxrad;
      zmin=odz-maxrad;
      zmax=odz+maxrad;

      assert (cctk_nghostzones[0]>=2);
      assert (cctk_nghostzones[1]>=2);
      assert (cctk_nghostzones[2]>=2);

      npts=0;
      for (i=2;i<ni-2;i++)
       for (j=2;j<nj-2;j++)
        for (k=2;k<nk-2;k++) {
          m=CCTK_GFINDEX3D(cctkGH,i,j,k);
          if ((  x[m]<=xmax&&x[m]>=xmin
               &&y[m]<=ymax&&y[m]>=ymin
               &&z[m]<=zmax&&z[m]>=zmin) &&
              (!respect_emask || (
               (emask[CCTK_GFINDEX3D(cctkGH,i-1,j,k)]>0.4) &&
               (emask[CCTK_GFINDEX3D(cctkGH,i-2,j,k)]>0.4) &&
               (emask[CCTK_GFINDEX3D(cctkGH,i+1,j,k)]>0.4) &&
               (emask[CCTK_GFINDEX3D(cctkGH,i+2,j,k)]>0.4) &&
               (emask[CCTK_GFINDEX3D(cctkGH,i,j-1,k)]>0.4) &&
               (emask[CCTK_GFINDEX3D(cctkGH,i,j-2,k)]>0.4) &&
               (emask[CCTK_GFINDEX3D(cctkGH,i,j+1,k)]>0.4) &&
               (emask[CCTK_GFINDEX3D(cctkGH,i,j+2,k)]>0.4) &&
               (emask[CCTK_GFINDEX3D(cctkGH,i,j,k-1)]>0.4) &&
               (emask[CCTK_GFINDEX3D(cctkGH,i,j,k-2)]>0.4) &&
               (emask[CCTK_GFINDEX3D(cctkGH,i,j,k+1)]>0.4) &&
               (emask[CCTK_GFINDEX3D(cctkGH,i,j,k+2)]>0.4)
              ))
             )
          {
            npts++;
          }
        }

      xa=(CCTK_REAL *)   malloc(npts*sizeof(CCTK_REAL));
      ya=(CCTK_REAL *)   malloc(npts*sizeof(CCTK_REAL));
      za=(CCTK_REAL *)   malloc(npts*sizeof(CCTK_REAL));
      rads=(CCTK_REAL *) malloc(npts*sizeof(CCTK_REAL));
      inds=(CCTK_INT *)  malloc(npts*sizeof(CCTK_INT));

      for (i=0;i<npts;i++) {
        rads[i]=0;
        inds[i]=0;
      }

      l=0;
      for (i=2;i<ni-2;i++)
       for (j=2;j<nj-2;j++)
        for (k=2;k<nk-2;k++) {
          m=CCTK_GFINDEX3D(cctkGH,i,j,k);
          if ((  x[m]<=xmax&&x[m]>=xmin
               &&y[m]<=ymax&&y[m]>=ymin
               &&z[m]<=zmax&&z[m]>=zmin) &&
              (!respect_emask || (
               (emask[CCTK_GFINDEX3D(cctkGH,i-1,j,k)]>0.4) &&
               (emask[CCTK_GFINDEX3D(cctkGH,i-2,j,k)]>0.4) &&
               (emask[CCTK_GFINDEX3D(cctkGH,i+1,j,k)]>0.4) &&
               (emask[CCTK_GFINDEX3D(cctkGH,i+2,j,k)]>0.4) &&
               (emask[CCTK_GFINDEX3D(cctkGH,i,j-1,k)]>0.4) &&
               (emask[CCTK_GFINDEX3D(cctkGH,i,j-2,k)]>0.4) &&
               (emask[CCTK_GFINDEX3D(cctkGH,i,j+1,k)]>0.4) &&
               (emask[CCTK_GFINDEX3D(cctkGH,i,j+2,k)]>0.4) &&
               (emask[CCTK_GFINDEX3D(cctkGH,i,j,k-1)]>0.4) &&
               (emask[CCTK_GFINDEX3D(cctkGH,i,j,k-2)]>0.4) &&
               (emask[CCTK_GFINDEX3D(cctkGH,i,j,k+1)]>0.4) &&
               (emask[CCTK_GFINDEX3D(cctkGH,i,j,k+2)]>0.4)
              ))
             )
          {
            xa[l]=x[m];
            ya[l]=y[m];
            za[l]=z[m];
            inds[l]=m;
            l++;
          }
        }

      ierr=HorizonRadiusInDirection(horizon_number[s],
                                    npts,
                                    xa, ya, za,  rads);
      assert(!ierr);

      for (i=0;i<npts;i++) {
        radp=sqrt((xa[i]-odx)*(xa[i]-odx)+(ya[i]-ody)*(ya[i]-ody)+
                  (za[i]-odz)*(za[i]-odz));
        if (radp<=rads[i]+ah_radius_offset) {
          epsdisA[inds[i]]=epsdis+ ah_slope*(rads[i]+ah_radius_offset-radp);
          if (epsdisA[inds[i]] > ah_max_epsdis) {
            epsdisA[inds[i]] = ah_max_epsdis;
          }
        }
      }

      free(xa);
      free(ya);
      free(za);
      free(rads);
      free(inds);
    }
  }

}
