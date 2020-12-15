#include <stdio.h>
#include <cassert>
#include <cmath>
#include <mpi.h>
#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <util_Table.h>
#include <cctk_Functions.h>
#include <cctk_Faces.h>
#include <Symmetry.h>
#include <loopcontrol.h>
#include <time.h>
#include "zst_find.hh"

using namespace ZST;

extern "C" {
  void zst_setup(CCTK_ARGUMENTS);
  void ZelmaniShockTracker2_Find(CCTK_ARGUMENTS);
  void zst_set_origin_init(CCTK_ARGUMENTS);
}

double alpha_grid_zst(double rmin, double rmax, int nzones,
		      double drmin, double prec);

void zst_setup(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;  

  if(CCTK_EQUALS(symm,"octant")) {
    if(ntheta > 1)
      *dtheta =  M_PI / (ntheta-1) / 2.0;
    else
      *dtheta = 0.0;
    if(nphi > 1) 
      *dphi =   M_PI / (nphi  -1) / 2.0;
    else
      *dphi = 0.0;
  } else if CCTK_EQUALS(symm,"bitant") {
    if(ntheta > 1)
      *dtheta =   M_PI / (ntheta-1) / 2.0;
    else
      *dtheta = 0.0;
    if(nphi > 1)
      *dphi   = 2*M_PI / (nphi  -1);
    else
      *dphi = 0.0;
  } else if CCTK_EQUALS(symm,"full") {
    if(ntheta > 1) 
      *dtheta =   M_PI / (ntheta-1);
    else
      *dtheta = 0.0;
    if(nphi > 1)
      *dphi = 2*M_PI / (nphi  -1);
    else
      *dphi = 0.0;
  } else {
    CCTK_WARN(0,"This symmetry is unknown");
  }
  
  // inner drad
  *drad   =  rad_max / (nrad - 1);
  // outer drad is tougher, since it's nonequidistant
  double prec = 1.0e-8;
  double alpha = alpha_grid_zst(rad_max,rad_max_outer,nrad_outer,*drad,prec);
  if(alpha < 0.0) {
    CCTK_WARN(0,"Grid setup failed!");
  }

  // set up radial coordinates
  for(int i=0;i<nrad;i++) {
    rad[i] = i*(*drad);
    //fprintf(stdout,"rad: %5d %15.6E\n",i,rad[i]);
  }
  // outer, non-equidistant part
  double dr2 = *drad;
  for(int i=nrad;i<nrad+nrad_outer;i++) {
    rad[i] = rad[i-1] + dr2;
    dr2 = dr2 * alpha;
    //  fprintf(stdout,"rad: %5d %15.6E\n",i,rad[i]);
  }
  // sanity check; tolerate error of not more than 1%
  if( abs(rad[nrad+nrad_outer-1] - rad_max_outer) > 0.01e0*rad_max_outer ) {
    CCTK_VInfo(CCTK_THORNSTRING,"rad[nrad+nrad_outer-1]: %18.9E    rad_max_outer: %18.9E",
	       rad[nrad+nrad_outer-1], rad_max_outer);
    CCTK_WARN(0,"Setting up radius array failed! Check your settings!");
  }

  for(int i=0;i<nphi;i++) {
    phi[i] = i*(*dphi);
    // fprintf(stdout,"phi: %5d %15.6E\n",i,phi[i]);
  }
  for(int i=0;i<ntheta;i++) {
    theta[i] = i*(*dtheta);
    //fprintf(stdout,"theta: %5d %15.6E\n",i,theta[i]);
  }

  // Calculate coordinates
  //  LC_LOOP3 (ZLtau_get_rays,
  //          ii,jj,kk, 0,0,0, nrad,ntheta,nphi, nrad,ntheta,nphi)
  #pragma omp parallel for
  for(int ii=0;ii<nrad+nrad_outer;ii++)
    for(int jj=0;jj<ntheta;jj++)
      for(int kk=0;kk<nphi;kk++)
  {
    int const iind3d = ii+(nrad+nrad_outer)*(jj+ntheta*kk);
    CCTK_REAL const xrad   = rad[ii];
    CCTK_REAL const xtheta = jj * (*dtheta);
    CCTK_REAL const xphi   = kk * (*dphi);
    CCTK_REAL xx, yy, zz;
    
    spher2cart (xrad,xtheta,xphi, &xx,&yy,&zz, x0, y0, z0);
    //    fprintf(stdout,"%5d %5d %5d %15.6E %15.6E %15.6E\n",
    //	    ii,jj,kk,xx,yy,zz);
    zst_x[iind3d] = xx;
    zst_y[iind3d] = yy;
    zst_z[iind3d] = zz;
  } //LC_ENDLOOP3 (ZLtau_get_rays);
      //} LC_ENDLOOP3 (ZLtau_get_rays);

  // some informative output
  CCTK_VInfo (CCTK_THORNSTRING,"Radii:");
  for(int ii=0;ii<nrad+nrad_outer;ii++) {
    if(ii<1)
      fprintf(stdout,"%4d %15.6E M\n",ii,rad[ii]);
    else
      fprintf(stdout,"%4d %15.6E M dr= %15.6E\n",ii,rad[ii],rad[ii]-rad[ii-1]);
  }
  CCTK_VInfo (CCTK_THORNSTRING,"Theta angles:");
  for(int jj=0;jj<ntheta;jj++) {
    fprintf(stdout,"%4d %15.6E Pi\n",jj,theta[jj]/M_PI);
  }
  CCTK_VInfo (CCTK_THORNSTRING,"Phi angles:");
  for(int kk=0;kk<nphi;kk++) {
    fprintf(stdout,"%4d %15.6E Pi\n",kk,phi[kk]/M_PI);
  }

  //  *have_interp_data = 0;

}

void ZelmaniShockTracker2_Find(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // don't do anything if we are not interested in doing anything
  if(!*dostuff) return;

  clock_t start, end;
  double elapsed;
  start = clock();
  CCTK_Info(CCTK_THORNSTRING, "Start interpolating for shock radius");
  // Interpolate
  {
    int const interp_handle = CCTK_InterpHandle (interpolator);
    assert (interp_handle >= 0);
    int const options_handle =
      Util_TableCreateFromString (interpolator_options);
    assert (options_handle >= 0);
    int const coords_handle = CCTK_CoordSystemHandle (coordinate_system);
    assert (coords_handle >= 0);
    
    // interpolate ONLY points on proc 0, then
    // do an MPI broadcast of the results
    int npoints = 0;
    if(CCTK_MyProc(cctkGH) == 0) {
      npoints = (nrad+nrad_outer)*ntheta*nphi;
    }

    void const *const interp_coords[] = { zst_x, zst_y, zst_z };
    
    CCTK_INT const input_array_indices[] = {
      CCTK_VarIndex ("HydroBase::vel[0]"),
      CCTK_VarIndex ("HydroBase::vel[1]"),
      CCTK_VarIndex ("HydroBase::vel[2]"),
      CCTK_VarIndex ("HydroBase::entropy")
    };
    int ninputs =
      sizeof input_array_indices / sizeof *input_array_indices;

    CCTK_INT const output_array_types[] ={
      CCTK_VARIABLE_REAL,
      CCTK_VARIABLE_REAL,
      CCTK_VARIABLE_REAL,
      CCTK_VARIABLE_REAL
    };
    assert (sizeof output_array_types / sizeof *output_array_types == ninputs);

    void *const output_arrays[] = { 
      zst_vx, zst_vy, zst_vz, zst_ent
    };

    assert (sizeof output_arrays / sizeof *output_arrays == ninputs);

    // skip interpolation of entropy when mode == simple!
    if (CCTK_EQUALS(mode, "simple"))
       --ninputs;

    int const ierr =
      CCTK_InterpGridArrays (cctkGH, 3,
                             interp_handle, options_handle, coords_handle,
                             npoints, CCTK_VARIABLE_REAL, interp_coords,
                             ninputs, input_array_indices,
                             ninputs, output_array_types, output_arrays);
    assert (not ierr);

    Util_TableDestroy (options_handle);
  }

  // now proc 0 broadcasts the interpolation results
  {
    int npoints = (nrad+nrad_outer)*ntheta*nphi;
    int err = MPI_Bcast((void*)zst_vx,npoints,MPI_DOUBLE,0,MPI_COMM_WORLD);
    err += MPI_Bcast((void*)zst_vy,npoints,MPI_DOUBLE,0,MPI_COMM_WORLD);
    err += MPI_Bcast((void*)zst_vz,npoints,MPI_DOUBLE,0,MPI_COMM_WORLD);
    if (CCTK_EQUALS(mode, "full"))
       err += MPI_Bcast((void*)zst_ent,npoints,MPI_DOUBLE,0,MPI_COMM_WORLD);
    assert(err==0);
  }

  end = clock();
  elapsed = ((double) (end - start)) / CLOCKS_PER_SEC;
  CCTK_VInfo (CCTK_THORNSTRING,"time used in proc 0: %5.2f",elapsed);
  CCTK_Info(CCTK_THORNSTRING, "Done interpolating for shock radius");
  //  *have_interp_data = 1;


  if (CCTK_EQUALS(mode, "full")) {

  // now compute the radial velocity
#pragma omp parallel for
  for (int kk=0; kk<nphi; ++kk) {
    for (int jj=0; jj<ntheta; ++jj) {
      for (int ii=0; ii<nrad+nrad_outer;ii++) {
        int const iind3d = ii+(nrad+nrad_outer)*(jj+ntheta*kk);

	zst_vr[iind3d] = zst_vx[iind3d]*cos(phi[kk])*sin(theta[jj]) 
	  + zst_vy[iind3d]*sin(phi[kk])*sin(theta[jj])
	  + zst_vz[iind3d]*cos(theta[jj]);
	
	if(ii>0 && ii<(nrad+nrad_outer-1)) {
	  int const ip = ii+(nrad+nrad_outer)*(jj+ntheta*kk) + 1;
	  int const im = ii+(nrad+nrad_outer)*(jj+ntheta*kk) - 1;
	  zst_dent[iind3d] = zst_ent[ip] - zst_ent[ii]; // forward difference
	}

      }
    }
  }

  } else if (CCTK_EQUALS(mode, "simple")) {
  
  // now compute the radial velocity
#pragma omp parallel for
  for (int kk=0; kk<nphi; ++kk) {
    for (int jj=0; jj<ntheta; ++jj) {
      for (int ii=0; ii<nrad+nrad_outer;ii++) {
        int const iind3d = ii+(nrad+nrad_outer)*(jj+ntheta*kk);

	zst_vr[iind3d] = zst_vx[iind3d]*cos(phi[kk])*sin(theta[jj]) 
	  + zst_vy[iind3d]*sin(phi[kk])*sin(theta[jj])
	  + zst_vz[iind3d]*cos(theta[jj]);
	
      }
    }
  }
  
  }

  // now find the shock
  // walk from the outside in, check for first
  // point at which entropy > 5 (or some other value); check for min velocity
  // if there is no entropy > 5 (or some other value), use min velocity point
  // else use the maximum radius of the two

  double minrad;
  if(cctk_time - (*bouncetime + start_time_inner_rad) < 0)
    minrad = 0.0;
  else
    minrad = inner_min_radius;

  if (CCTK_EQUALS(mode, "full")) {

#pragma omp parallel for
  for (int kk=0; kk<nphi; ++kk) {
    for (int jj=0; jj<ntheta; ++jj) {

      // going from inside out, finding innermost entropy jump
      // outside some radius
      int ishock = 0;
      int ientjump = 0;
      if(minrad > 0.0) {
	{
	  int ii = 0;
	  while(ii < nrad+nrad_outer && ientjump == 0) {
	    if(rad[ii] > 40.0) {
	      int const iind3d = ii+(nrad+nrad_outer)*(jj+ntheta*kk);
	      if (zst_ent[iind3d] < entropy_value) {
		ientjump = ii;
	      } 
	    }
	    ii++;
	  }
	}
      }
      ishock = ientjump;

      // if something went wrong, walk down and look for lowest
      // velocity
      if(ientjump >= nrad+nrad_outer - 1 ) {
	int ivelmin = nrad+nrad_outer-1 ;
	double velmin = 1.0;
	int ii = (nrad+nrad_outer-1);
	while(ii > 0) {
	  int const iind3d = ii+(nrad+nrad_outer)*(jj+ntheta*kk);
	  if(zst_vr[iind3d] < velmin && rad[ii] > minrad) {
	    velmin = zst_vr[iind3d];
	    ivelmin = ii;
	  }
	  ii--;
	}
      ishock = ivelmin;
      }

// hack
      // going outside in, find outermost entropy jump
      //int ishock = 0;
      //int ientjump = 0;
//      ishock = 0;
//      ientjump = 0;
//      {
//	int ii = (nrad+nrad_outer-1);
//	while(ii > 0 && ientjump == 0) {
//	  int const iind3d = ii+(nrad+nrad_outer)*(jj+ntheta*kk);
//	  if (zst_ent[iind3d] >= entropy_value ) {
//	    ientjump = ii;
//	  }
//	  ii--;
//	}
//      }
//
//      double minvel = -1.0e10;
//      int iminvel = 0;
//      for (int ii=(nrad+nrad_outer-1); ii >= 0;--ii) {
//	int const iind3d = ii+(nrad+nrad_outer)*(jj+ntheta*kk);
//	if(zst_vr[iind3d] < minvel && rad[ii] > minrad) {
//	  minvel = zst_vr[iind3d];
//	  iminvel = ii;
//	}
//      }
//
//      if(ientjump == 0) {
//	ishock = iminvel;
//      } else {
//	ishock = max(ientjump,iminvel);
//      }
//
// hack
      int iind2d = jj + ntheta*kk;
      shockpos[iind2d] = rad[ishock];

    }
  }

  } else if (CCTK_EQUALS(mode, "simple")) {

#pragma omp parallel for
  for (int kk=0; kk<nphi; ++kk) {
    for (int jj=0; jj<ntheta; ++jj) {

      // walk down and look for lowest velocity
        int ishock = 0;
	int ivelmin = nrad+nrad_outer-1 ;
	double velmin = 1.0;
	int ii = (nrad+nrad_outer-1);
	while(ii > 0) {
	  int const iind3d = ii+(nrad+nrad_outer)*(jj+ntheta*kk);
	  if(zst_vr[iind3d] < velmin && rad[ii] > minrad) {
	    velmin = zst_vr[iind3d];
	    ivelmin = ii;
	  }
	  ii--;
	}
        ishock = ivelmin;
      
        int iind2d = jj + ntheta*kk;
        shockpos[iind2d] = rad[ishock];
      
    }
  }

  }

  // now get min, max, average radii
  *shockmax = 0.0;
  *shockmin = 1.0e20;
  *shockav = 0.0;
  for (int kk=0; kk<nphi; ++kk) {
    for (int jj=0; jj<ntheta; ++jj) {
      int iind2d = jj + ntheta*kk;
      *shockav = *shockav + shockpos[iind2d];
      *shockmax = max(*shockmax,shockpos[iind2d]);
      *shockmin = min(*shockmin,shockpos[iind2d]);
    }
  }
  *shockav = *shockav / (ntheta*nphi);


  if(verbose) {
    CCTK_VInfo(CCTK_THORNSTRING,"average shock radius: %15.6E",*shockav);
    CCTK_VInfo(CCTK_THORNSTRING,"maximum shock radius: %15.6E",*shockmax);
    CCTK_VInfo(CCTK_THORNSTRING,"minimum shock radius: %15.6E",*shockmin);
  }

  return;

}

void zst_set_origin_init(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  *x0 = 0.0;
  *y0 = 0.0;
  *z0 = 0.0;


  for(int i=0;i<(nrad+nrad_outer)*nphi*ntheta;i++) {
    zst_vx[i]  = 0.0;
    zst_vy[i]  = 0.0;
    zst_vz[i]  = 0.0;
    zst_ent[i]  = 0.0;
    zst_dent[i]  = 0.0;
  }

  return;
}
