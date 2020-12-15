
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#define velx (&vel[0*cctk_ash[0]*cctk_ash[1]*cctk_ash[2]])
#define vely (&vel[1*cctk_ash[0]*cctk_ash[1]*cctk_ash[2]])
#define velz (&vel[2*cctk_ash[0]*cctk_ash[1]*cctk_ash[2]])

static void SpatialDeterminant(CCTK_REAL gxx,
                               CCTK_REAL gxy,
                               CCTK_REAL gxz,
                               CCTK_REAL gyy,
                               CCTK_REAL gyz,
                               CCTK_REAL gzz,
                               CCTK_REAL *detg)
{

  *detg = -gxz*gxz*gyy + 2.0*gxy*gxz*gyz 
    - gxx*gyz*gyz - gxy*gxy*gzz + gxx*gyy*gzz;

  return;
}


static void UpperMetric(CCTK_REAL gxx, 
                        CCTK_REAL gxy, 
                        CCTK_REAL gxz, 
                        CCTK_REAL gyy, 
                        CCTK_REAL gyz, 
                        CCTK_REAL gzz, 
                        CCTK_REAL det,
                        CCTK_REAL *uxx, 
                        CCTK_REAL *uxy, 
                        CCTK_REAL *uxz, 
                        CCTK_REAL *uyy, 
                        CCTK_REAL *uyz, 
                        CCTK_REAL *uzz)
{
  
  *uxx=(-gyz*gyz + gyy*gzz)/det;
  *uxy=(gxz*gyz - gxy*gzz)/det;
  *uyy=(-gxz*gxz + gxx*gzz)/det;
  *uxz=(-gxz*gyy + gxy*gyz)/det;
  *uyz=(gxy*gxz - gxx*gyz)/det;
  *uzz=(-gxy*gxy + gxx*gyy)/det;
  
  return;
  
}


void ZelmaniAnalysis_Hydro_Local(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int nx = cctkGH->cctk_lsh[0];
  int ny = cctkGH->cctk_lsh[1];
  int nz = cctkGH->cctk_lsh[2];
  
  int ghost_size_x = 0; //cctkGH->cctk_nghostzones[0];
  int ghost_size_y = 0; //cctkGH->cctk_nghostzones[1];
  int ghost_size_z = 0; //cctkGH->cctk_nghostzones[2];

  CCTK_REAL dx = CCTK_DELTA_SPACE(0);
  CCTK_REAL dy = CCTK_DELTA_SPACE(1);
  CCTK_REAL dz = CCTK_DELTA_SPACE(2);

  /* helpers */
  CCTK_REAL h;
  CCTK_REAL detg;
  int index,i,j,k;
  CCTK_REAL gtx,gty,gtz;
  CCTK_REAL uxx, uxy, uxz, uyy, uyz, uzz;
  CCTK_REAL Jx, Jy, Omega, cyl_rad_squ;
  CCTK_REAL ulx, uly, ulz;
  CCTK_REAL vlx, vly, vlz;
  CCTK_REAL* is_inner_core = NULL;
  CCTK_REAL densfac;

  if (*volume_form_state == 0)
  {
     CCTK_WARN(0, "Tell your 'Coordinates' implementation to activate storage for the volume form!");
  }

  if (*jacobian_state == 0)
  {
     CCTK_WARN(0, "Tell your 'Coordinates' implementation to activate storage for Jacobians!");
  }
  

  // Check if InnerCore is active; if so, get pointer
  // to is_inner_core GF that tells us if we are in the
  // inner core
  if(do_ic_analysis) {
    if(CCTK_IsThornActive("InnerCore")) {
      int varindex = CCTK_VarIndex("InnerCore::is_inner_core");
      is_inner_core = (CCTK_REAL*) CCTK_VarDataPtrI (cctkGH, 0, varindex);
      assert(is_inner_core);
    } else {
      CCTK_WARN(0,"You are trying to do inner core analysis, but Thorn InnerCore is not active!");
    }
  }

  int do_enu = CCTK_IsThornActive("ZelmaniM1");
  int *ng = NULL,*ns = NULL;
  CCTK_REAL *Eg = NULL;
  if(do_enu) {
    ng = (int*)CCTK_ParameterGet("ngroups", "ZelmaniM1",NULL);
    ns = (int*)CCTK_ParameterGet("nspecies","ZelmaniM1",NULL);
    int vin = CCTK_VarIndex("ZelmaniM1::enu[0]");
    Eg = (CCTK_REAL *)(CCTK_VarDataPtrI(cctkGH,0,vin));
  }


#pragma omp parallel for private(k,j,i,index,gtx,gty,gtz,detg,uxx,uxy,uxz,uyy,uyz,uzz, \
				 ulx,uly,ulz,vlx,vly,vlz,h,Jx,Jy,Omega,cyl_rad_squ,densfac)
  for(k=ghost_size_z; k < nz-ghost_size_z;k++)
    for(j=ghost_size_y; j < ny-ghost_size_y;j++)
      for(i=ghost_size_x; i < nx-ghost_size_x;i++) {

	index = CCTK_GFINDEX3D(cctkGH,i,j,k);

	/* determinant of the 3-metric */
	SpatialDeterminant(gxx[index],gxy[index],gxz[index],
			   gyy[index],gyz[index],gzz[index],
			   &detg);	


	/* Contravariant components of ADM 3-metric. */
	UpperMetric(gxx[index],gxy[index],gxz[index],gyy[index],
		    gyz[index],gzz[index],detg,&uxx,&uxy,&uxz,
		    &uyy,&uyz,&uzz);


	/* Covariant components of ADM 4-metric. */
	gtx = gxx[index] * betax[index] + 
	  gxy[index] * betay[index] + 
	  gxz[index] * betaz[index];

	gty = gxy[index] * betax[index] + 
	  gyy[index] * betay[index] + 
	  gyz[index] * betaz[index];

	gtz = gxz[index] * betax[index] + 
	  gyz[index] * betay[index] + 
	  gzz[index] * betaz[index];

	/* Covariant components of four-velocity. */
	ulx = w_lorentz[index] * 
	  (gtx / alp[index] + 
	   gxx[index] * (velx[index] - betax[index] / alp[index]) + 
	   gxy[index] * (vely[index] - betay[index] / alp[index]) + 
	   gxz[index] * (velz[index] - betaz[index] / alp[index]));

	uly = w_lorentz[index] * 
	  (gty / alp[index] + 
	   gxy[index] * (velx[index] - betax[index] / alp[index]) + 
	   gyy[index] * (vely[index] - betay[index] / alp[index]) + 
	   gyz[index] * (velz[index] - betaz[index] / alp[index]));

	ulz = w_lorentz[index] * 
	  (gtz / alp[index] + 
	   gxz[index] * (velx[index] - betax[index] / alp[index]) + 
	   gyz[index] * (vely[index] - betay[index] / alp[index]) + 
	   gzz[index] * (velz[index] - betaz[index] / alp[index]));

	vlx = gxx[index]*velx[index] + gxy[index]*vely[index] + gxz[index]*velz[index];
	vly = gxy[index]*velx[index] + gyy[index]*vely[index] + gyz[index]*velz[index];
	vlz = gxz[index]*velx[index] + gyz[index]*vely[index] + gzz[index]*velz[index];


	h = 1.0 + eps[index] + press[index] / rho[index];
	/* Momentum constraint matter sources. */
	Jx = rho[index] * w_lorentz[index] * h * ulx;
	Jy = rho[index] * w_lorentz[index] * h * uly;

	/* Cylindrical radius squared. */
	cyl_rad_squ = x[index]*x[index] + y[index]*y[index];

	/* Angular velocity. */
	/* (assumes del_phi to be the direction defining rotation) */
	if (cyl_rad_squ > 1e-10) {
	  Omega = (x[index]*(alp[index]*vely[index]-betay[index]) - 
		   y[index]*(alp[index]*velx[index]-betax[index])) /
		   cyl_rad_squ;
	} else {
	  Omega = 0.0;
	}
	if(do_omega_3D) {
	  angular_velocity[index] = Omega;
	}


	/* -- Saijo's variant of the total angular momentum -- */
	angular_momentum_local[index] = sqrt(detg)*
	  (x[index]*Jy - y[index]*Jx) * volume_form[index];

	kinetic_energy_local[index] = 0.5 * Omega * 
	  angular_momentum_local[index];

	/* -- total kinetic energy -- */
	total_kinetic_energy_local[index] = (w_lorentz[index]*w_lorentz[index] - 1) * 
	  sqrt(detg) * rho[index] * volume_form[index];

	/* Masses */
        // transformation factor for cons. density which is a scalar tensor density of weight 1
        densfac = fabs((  J11[index] * J22[index] * J33[index]
                        + J12[index] * J23[index] * J31[index]
                        + J13[index] * J21[index] * J32[index]
                        - J11[index] * J23[index] * J32[index]
                        - J12[index] * J21[index] * J33[index]
                        - J13[index] * J22[index] * J31[index]));

	baryonic_mass_local[index] = densfac*dens[index] * volume_form[index];
	
	proper_mass_local[index] = densfac*dens[index] * (1.0 + eps[index]) * volume_form[index];
        
	// Include neutrino energy density if M1 thorn is active
	double Enutot = 0.0;
        if (do_enu) {
	  for (int ig=0;ig<(*ng)*(*ns);ig++) {
              int index = CCTK_VectGFIndex3D(cctkGH,i,j,k,ig);
	      Enutot += Eg[index];
	  }
#if 0          
	  if (i==10&&j==10&&k==10) 	
	    CCTK_VInfo(CCTK_THORNSTRING, "Analysis: %3d %3d %18.9E",*ng,*ns,Enutot);
#endif 
	}
        
	if (sqrt(detg)>1.e-16) Enutot = Enutot/sqrt(detg);
	else Enutot = Enutot/1.e-16;

	gravitational_mass_local[index] = (rho[index] * h * w_lorentz[index]*
					   w_lorentz[index]  - press[index] + Enutot)
	  *pow(detg,5.0e0/12.0e0) * volume_form[index];
       
	if(number_of_spheres>0) {
	  for(int m=0;m<number_of_spheres;m++) {
	    int i4D = CCTK_VECTGFINDEX3D(cctkGH,i,j,k,m);
	    if(r[index] <= radii[m]) {
	      adm_mass_volume_local[i4D] = gravitational_mass_local[index];
	      baryonic_mass_interior_local[i4D] = baryonic_mass_local[index];
	      if(do_enu_interior) {
		enu_interior_local[i4D] = Enutot*pow(detg,5.0e0/12.0e0) * volume_form[index];
	      }
	    } else {
	      adm_mass_volume_local[i4D] = 0.0;
	      baryonic_mass_interior_local[i4D] = 0.0;
	      if(do_enu_interior) {
		enu_interior_local[i4D] = 0.0;
	      }
	    }
	  }
	}


	//	gravitational_mass_local[index] = (-2.0e0 * rho[index] * h * w_lorentz[index]*
	/// w_lorentz[index] * (-alp[index] + betax[index]*vlx + 
	//  betay[index]*vly + betaz[index]*vlz)
	//  - alp[index]*rho[index] * h + alp[index]*2.0e0*press[index]) * sqrt(detg);

	if(rho[index] >= 1.61930347e-06) {
	  gravitational_mass_local_in1e12[index] = (rho[index] * h * w_lorentz[index]*
			w_lorentz[index]  - press[index] + Enutot)
	                *pow(detg,5.0e0/12.0e0) * volume_form[index];
	  
	  baryonic_mass_local_in1e12[index] = densfac*dens[index] * volume_form[index];

	} else {
	  gravitational_mass_local_in1e12[index] = 0.0e0 * volume_form[index];
	  baryonic_mass_local_in1e12[index] = 0.0e0 * volume_form[index];
	}
	
	// compute the same stuff in the inner core
	if(do_ic_analysis) {
	  if(is_inner_core[index] > 0.0e0) {
	    proper_mass_ic_local[index] = proper_mass_local[index];
	    kinetic_energy_ic_local[index] = kinetic_energy_local[index];
	    angular_momentum_ic_local[index] = angular_momentum_local[index];
	    gravitational_mass_ic_local[index] = gravitational_mass_local[index];
	  } else {
	    proper_mass_ic_local[index] = 0.0e0;
	    kinetic_energy_ic_local[index] = 0.0e0;
	    angular_momentum_ic_local[index] = 0.0e0;
	    gravitational_mass_ic_local[index] = 0.0e0;
	  }
	}
      }
}

void ZelmaniAnalysis_Hydro_Global(CCTK_ARGUMENTS)
{
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    CCTK_INT reduction_handle;
    CCTK_REAL radius; //, grid_spacing_product;
    
    reduction_handle = CCTK_ReductionHandle("sum");
    if (reduction_handle < 0)
        CCTK_WARN(0, "Unable to get reduction handle.");

    if (CCTK_Reduce(cctkGH, -1, reduction_handle, 1,
                    CCTK_VARIABLE_REAL,
                    baryonic_mass, 1,
                    //CCTK_VarIndex("GRHydro::dens")))
                    CCTK_VarIndex("ZelmaniAnalysis::baryonic_mass_local")))
        CCTK_WARN(0, "Error while reducing GRHydro::dens");

    if (CCTK_Reduce(cctkGH, -1, reduction_handle, 1,
                    CCTK_VARIABLE_REAL,
                    proper_mass, 1,
                    CCTK_VarIndex("ZelmaniAnalysis::proper_mass_local")))
        CCTK_WARN(0, "Error while reducing ZelmaniAnalysis::proper_mass_local");

    if (CCTK_Reduce(cctkGH, -1, reduction_handle, 1,
                    CCTK_VARIABLE_REAL,
                    gravitational_mass, 1,
                    CCTK_VarIndex("ZelmaniAnalysis::gravitational_mass_local")))
        CCTK_WARN(0, "Error while reducing ZelmaniAnalysis::gravitational_mass_local");

    if (CCTK_Reduce(cctkGH, -1, reduction_handle, 1,
                    CCTK_VARIABLE_REAL,
                    gravitational_mass_in1e12, 1,
                    CCTK_VarIndex("ZelmaniAnalysis::gravitational_mass_local_in1e12")))
        CCTK_WARN(0, "Error while reducing ZelmaniAnalysis::gravitational_mass_local_in1e12");

    if (CCTK_Reduce(cctkGH, -1, reduction_handle, 1,
                    CCTK_VARIABLE_REAL,
                    baryonic_mass_in1e12, 1,
                    CCTK_VarIndex("ZelmaniAnalysis::baryonic_mass_local_in1e12")))
        CCTK_WARN(0, "Error while reducing ZelmaniAnalysis::baryonic_mass_local_in1e12");

    if (CCTK_Reduce(cctkGH, -1, reduction_handle, 1,
                    CCTK_VARIABLE_REAL,
                    angular_momentum, 1,
                    CCTK_VarIndex("ZelmaniAnalysis::angular_momentum_local")))
        CCTK_WARN(0, "Error while reducing ZelmaniAnalysis::angular_momentum_local_local");

    if (CCTK_Reduce(cctkGH, -1, reduction_handle, 1,
                    CCTK_VARIABLE_REAL,
                    kinetic_energy, 1,
                    CCTK_VarIndex("ZelmaniAnalysis::kinetic_energy_local")))
        CCTK_WARN(0, "Error while reducing ZelmaniAnalysis::kinetic_energy_local_local");


    if (CCTK_Reduce(cctkGH, -1, reduction_handle, 1,
                    CCTK_VARIABLE_REAL,
                    total_kinetic_energy, 1,
                    CCTK_VarIndex("ZelmaniAnalysis::total_kinetic_energy_local")))
        CCTK_WARN(0, "Error while reducing ZelmaniAnalysis::total_kinetic_energy_local");

    if (CCTK_IsThornActive("ML_ADMQuantities")) {
       if (CCTK_Reduce(cctkGH, -1, reduction_handle, 1,
                       CCTK_VARIABLE_REAL,
                       Madm_red, 1,
                       CCTK_VarIndex("ML_ADMQuantities::Madm")))
           CCTK_WARN(0, "Error while reducing ML_ADMQuantities::Madm");

       if (CCTK_Reduce(cctkGH, -1, reduction_handle, 1,
                       CCTK_VARIABLE_REAL,
                       Jadm_z_red, 1,
                       CCTK_VarIndex("ML_ADMQuantities::Jadm3")))
           CCTK_WARN(0, "Error while reducing ML_ADMQuantities::Jadm3");
    }

//    grid_spacing_product = cctk_delta_space[0]*cctk_delta_space[1]*cctk_delta_space[2];

    if(do_ic_analysis) {
      // inner core reductions
      if (CCTK_Reduce(cctkGH, -1, reduction_handle, 1,
		      CCTK_VARIABLE_REAL,
                    kinetic_energy_ic, 1,
		      CCTK_VarIndex("ZelmaniAnalysis::kinetic_energy_ic_local")))
        CCTK_WARN(0, "Error while reducing ZelmaniAnalysis::kinetic_energy_local_local");
      
      if (CCTK_Reduce(cctkGH, -1, reduction_handle, 1,
		      CCTK_VARIABLE_REAL,
		      angular_momentum_ic, 1,
		      CCTK_VarIndex("ZelmaniAnalysis::angular_momentum_ic_local")))
        CCTK_WARN(0, "Error while reducing ZelmaniAnalysis::angular_momentum_ic_local");

      if (CCTK_Reduce(cctkGH, -1, reduction_handle, 1,
		      CCTK_VARIABLE_REAL,
		      gravitational_mass_ic, 1,
		      CCTK_VarIndex("ZelmaniAnalysis::gravitational_mass_ic_local")))
        CCTK_WARN(0, "Error while reducing ZelmaniAnalysis::gravitational_mass_ic_local");

      if (CCTK_Reduce(cctkGH, -1, reduction_handle, 1,
		      CCTK_VARIABLE_REAL,
		      proper_mass_ic, 1,
		      CCTK_VarIndex("ZelmaniAnalysis::proper_mass_ic_local")))
        CCTK_WARN(0, "Error while reducing ZelmaniAnalysis::proper_mass_ic_local");
/*
      *proper_mass_ic *= grid_spacing_product;
      *gravitational_mass_ic *= grid_spacing_product;
      *kinetic_energy_ic *= grid_spacing_product;
      *angular_momentum_ic *= grid_spacing_product;
*/
      *W_ic = (*proper_mass_ic) + (*kinetic_energy_ic) - (*gravitational_mass_ic);
      *T_over_W_ic = *kinetic_energy_ic / fabs(*W_ic);

    }

    if(number_of_spheres>0 && CCTK_QueryGroupStorage(cctkGH,"ZelmaniAnalysis::adm_mass_volume_local")) {
      for(int is=0;is<number_of_spheres;is++) {
	char varname[256];
	sprintf(varname,"%s[%d]","ZelmaniAnalysis::adm_mass_volume_local",is);
	if (CCTK_Reduce(cctkGH, -1, reduction_handle, 1,
			CCTK_VARIABLE_REAL,
			(void*)&(adm_mass_volume[is]), 1,
			CCTK_VarIndex(varname)))
	  CCTK_WARN(0, "Error while reducing ZelmaniAnalysis::adm_mass_volume_local");

	sprintf(varname,"%s[%d]","ZelmaniAnalysis::baryonic_mass_interior_local",is);
	if (CCTK_Reduce(cctkGH, -1, reduction_handle, 1,
			CCTK_VARIABLE_REAL,
			(void*)&(baryonic_mass_interior[is]), 1,
			CCTK_VarIndex(varname)))
	  CCTK_WARN(0, "Error while reducing ZelmaniAnalysis::baryonic_mass_interior_local");
	if(do_enu_interior) {
	  sprintf(varname,"%s[%d]","ZelmaniAnalysis::enu_interior_local",is);
	  if (CCTK_Reduce(cctkGH, -1, reduction_handle, 1,
			  CCTK_VARIABLE_REAL,
			  (void*)&(enu_interior[is]), 1,
			  CCTK_VarIndex(varname)))
	    CCTK_WARN(0, "Error while reducing ZelmaniAnalysis::enu_interior_local");
	}
      }
    }




/*
    *baryonic_mass *= grid_spacing_product;
    *baryonic_mass_in1e12 *= grid_spacing_product;
    *proper_mass *= grid_spacing_product;
    *gravitational_mass *= grid_spacing_product;
    *gravitational_mass_in1e12 *= grid_spacing_product;
    *angular_momentum *= grid_spacing_product;
    *kinetic_energy *= grid_spacing_product;
    *total_kinetic_energy *= grid_spacing_product;
*/

    *W = (*proper_mass) + (*kinetic_energy) - (*gravitational_mass);
    *T_over_W = *kinetic_energy / fabs(*W);

    //    CCTK_VInfo(CCTK_THORNSTRING, "T/|W|: %15.6E", *T_over_W);
/*
    *Madm_red *= grid_spacing_product;
    *Jadm_z_red *= grid_spacing_product;
*/
}

