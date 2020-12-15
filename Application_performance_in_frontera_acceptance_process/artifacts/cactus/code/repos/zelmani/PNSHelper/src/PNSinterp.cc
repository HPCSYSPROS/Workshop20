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
#include "PNSinterp.hh"


using namespace PNSinterp;

extern "C" {
  void PNSinterp_setup(CCTK_ARGUMENTS);
  void PNSinterp_setup_local(CCTK_ARGUMENTS);
  void PNSinterp_get_rays(CCTK_ARGUMENTS);
  void PNSinterp_register(CCTK_ARGUMENTS);
  void PNSinterp_set_origin_init(CCTK_ARGUMENTS);
}

double alpha_grid_pns(double rmin, double rmax, int nzones,
		  double drmin, double prec);

void PNSinterp_setup_local(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

}

void PNSinterp_register(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;  
  
}

void PNSinterp_setup(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;  

  CCTK_VInfo (CCTK_THORNSTRING,"Setup");
  
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
  double alpha = alpha_grid_pns(rad_max,rad_max_outer,nrad_outer,*drad,prec);
  if(alpha < 0.0) {
    CCTK_WARN(0,"Grid setup failed!");
  }

  // set up radial coordinates, shift by half a zone
  for(int i=0;i<nrad;i++) {
    rad[i] = i*(*drad) + 0.5*(*drad);
    //fprintf(stdout,"rad: %5d %15.6E\n",i,rad[i]);
  }
  // outer, non-equidistant part
  if(nrad_outer > 0) {
    CCTK_WARN(0,"non-equidistant outer grid currently not supported");
  }
  double dr2 = *drad;
  for(int i=nrad;i<nrad+nrad_outer;i++) {
    rad[i] = rad[i-1] + dr2;
    dr2 = dr2 * alpha;
    //  fprintf(stdout,"rad: %5d %15.6E\n",i,rad[i]);
  }

  for(int i=0;i<nphi;i++) {
    phi[i] = i*(*dphi);
    // fprintf(stdout,"phi: %5d %15.6E\n",i,phi[i]);
  }
  for(int i=0;i<ntheta;i++) {
    theta[i] = i*(*dtheta);
    //fprintf(stdout,"theta: %5d %15.6E\n",i,theta[i]);
  }

  if (CCTK_EQUALS(collect, "interpolate")) { 
    // Calculate coordinates
    //  LC_LOOP3 (PNSinterp_get_rays,
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
      pns_x[iind3d] = xx;
      pns_y[iind3d] = yy;
      pns_z[iind3d] = zz;
    } //LC_ENDLOOP3 (PNSinterp_get_rays);
        //} LC_ENDLOOP3 (PNSinterp_get_rays);
  }

  // some informative output
  if (verbose) {
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
  }
   
  *have_interp_data = 0;
}

void PNSinterp_get_rays(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if(update_GR_every < 0) return;
  
  if (verbose) 
  fprintf(stdout,"cctk_iteration: %d force_interp: %d\n",cctk_iteration,*force_interp);

  if (cctk_time < update_GR_switch_time) {
    if( ((cctk_iteration-1) % update_GR_every_start != 0) && !*force_interp) return;
  } else {
    if( ((cctk_iteration-1) % update_GR_every != 0) && !*force_interp) return;
  }



  *force_interp = 0;
  if (verbose)
  CCTK_Info(CCTK_THORNSTRING, "Interpolating for GR");

  clock_t start, end;
  double elapsed;
  start = clock();

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

    void const *const interp_coords[] = { pns_x, pns_y, pns_z };
    
    CCTK_INT const input_array_indices[] = {
      CCTK_VarIndex ("HydroBase::rho"),
      CCTK_VarIndex ("HydroBase::press"),
      CCTK_VarIndex ("HydroBase::temperature"),
      CCTK_VarIndex ("HydroBase::Y_e"),
      CCTK_VarIndex ("HydroBase::eps"),
      CCTK_VarIndex ("HydroBase::w_lorentz"),
      CCTK_VarIndex ("ADMBase::alp"),
      CCTK_VarIndex ("ADMBase::gxx")
    };
    int ninputs =
      sizeof input_array_indices / sizeof *input_array_indices;

    CCTK_INT const output_array_types[] = {
      CCTK_VARIABLE_REAL,
      CCTK_VARIABLE_REAL,
      CCTK_VARIABLE_REAL,
      CCTK_VARIABLE_REAL,
      CCTK_VARIABLE_REAL,
      CCTK_VARIABLE_REAL,
      CCTK_VARIABLE_REAL,
      CCTK_VARIABLE_REAL
    };
    assert (sizeof output_array_types / sizeof *output_array_types == ninputs);

    void *const output_arrays[] = { 
      pns_rho, pns_press, pns_temp, pns_ye, pns_eps, pns_w_lorentz, pns_old_alp, pns_old_gxx
    };

    assert (sizeof output_arrays / sizeof *output_arrays == ninputs);

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
    int err = MPI_Bcast((void*)pns_rho,npoints,MPI_DOUBLE,0,MPI_COMM_WORLD);
    err += MPI_Bcast((void*)pns_press,npoints,MPI_DOUBLE,0,MPI_COMM_WORLD);
    err += MPI_Bcast((void*)pns_temp,npoints,MPI_DOUBLE,0,MPI_COMM_WORLD);
    err += MPI_Bcast((void*)pns_ye,npoints,MPI_DOUBLE,0,MPI_COMM_WORLD);
    err += MPI_Bcast((void*)pns_eps,npoints,MPI_DOUBLE,0,MPI_COMM_WORLD);
    err += MPI_Bcast((void*)pns_w_lorentz,npoints,MPI_DOUBLE,0,MPI_COMM_WORLD);
    err += MPI_Bcast((void*)pns_old_alp,npoints,MPI_DOUBLE,0,MPI_COMM_WORLD);
    err += MPI_Bcast((void*)pns_old_gxx,npoints,MPI_DOUBLE,0,MPI_COMM_WORLD);
    assert(err==0);
  }

  end = clock();
  elapsed = ((double) (end - start)) / CLOCKS_PER_SEC;
  if (verbose) {
    CCTK_VInfo (CCTK_THORNSTRING,"time used in proc 0: %5.2f",elapsed);
    CCTK_Info(CCTK_THORNSTRING, "Done interpolating for GR");
  }
  *have_interp_data = 1;

  return;
}

void PNSinterp_set_origin_init(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  if (verbose) 
  CCTK_Info(CCTK_THORNSTRING,"PNSinterp_set_origin_init");

  *x0 = 0.0;
  *y0 = 0.0;
  *z0 = 0.0;

  *pns_mgrav = 0.0;
  
  if (CCTK_EQUALS(collect, "interpolate")) { 
    for(int i=0;i<(nrad+nrad_outer)*nphi*ntheta;i++) {
      pns_rho[i]  = 0.0;
      pns_eps[i]  = 0.0;
      pns_press[i]  = 0.0;
      pns_temp[i]  = 0.0;
      pns_ye[i]  = 0.0;
      pns_w_lorentz[i]  = 0.0;
      pns_old_alp[i] = 0.0;
      pns_old_gxx[i] = 0.0;
    }
  }

  for(int i=0;i<(nrad+nrad_outer);i++) {
    pns_av_rho[i]  = 0.0;
    pns_av_eps[i]  = 0.0;
    pns_av_press[i]  = 0.0;
    pns_av_temp[i]  = 0.0;
    pns_av_ye[i]  = 0.0;
    pns_av_w_lorentz[i]  = 0.0;
    pns_rho_star[i] = 0.0;
    pns_real_rho_star[i] = 0.0;
    pns_av_alp[i]  = 0.0;
    pns_av_gxx[i]  = 0.0;
    pns_alp[i]  = 1.0;
    pns_psi[i]  = 1.0;
  }
  metric_times[2] = 0.0;
  metric_times[1] = 0.0;
  metric_times[0] = 0.0;
  pns_mgrav_store[0] = 0.0;
  pns_mgrav_store[1] = 0.0;
  pns_mgrav_store[2] = 0.0;

  for(int i=0;i<((nrad+nrad_outer)*3);i++) {
    pns_psi_store[i] = 1.0;
    pns_alp_store[i] = 1.0;
  }


  *force_interp = 1;
  
  return;
}

double alpha_grid_pns(double rmin, double rmax, int nzones,
		  double drmin, double prec) {

  double alpha,alpha_p;
  double f,fp,dfda;
  double rad,rad_p,dr;
  int it;

  rad = rmin;
  alpha = 1.0e0;
  dr = drmin;
  for(int i=0;i<nzones;i++) {
    rad = rad + dr;
    dr = alpha * dr;
  }
  rad_p = rad;
  alpha_p = alpha;

  rad = rmin;
  alpha = 1.01e0;
  dr = drmin;
  for(int i=0;i<nzones;i++) {
    rad = rad + dr;
    dr = alpha * dr;
  }

  it = 0;
  f = rmax - rad;
  while( (fabs(f/rmax) > prec) && (it < 100) ) {
    dfda =  ( (rmax - rad) - (rmax - rad_p) ) / (alpha - alpha_p);
    rad_p = rad;
    alpha_p = alpha;
    alpha = alpha - f / dfda;
    rad = rmin;
    dr = drmin;
    for(int i=0;i<nzones;i++) {
      rad = rad + dr;
      dr = alpha * dr;
    }
    f = rmax - rad;
    //    printf("%d %15.6E %15.6E %15.6E %15.6E\n",it,alpha,f,fabs(f/rmax),prec);
    it++;
  }
  if (it >= 200) {
    alpha = -1.0e0;
  }

  return alpha;
}
