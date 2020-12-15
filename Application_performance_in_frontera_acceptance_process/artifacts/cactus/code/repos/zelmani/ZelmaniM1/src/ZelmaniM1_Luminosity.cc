#include <cassert>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <cerrno>
#include <vector>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include "ZelmaniM1.hh"

#define PI 3.14159265358979
#define G_GRAV 6.673e-8
#define C_L 2.99792e10
#define M_SUN 1.98892e33
#define M_TO_L 1.47672e5 // Convert from solar masses to cm
#define Mevfm3_TO_Ergcm3 1.60218e33 
#define INV_PRESS_GF 5.56082181777535e38

using namespace std;


extern "C" {
  void zm1_luminosity_init(CCTK_ARGUMENTS);
  void zm1_luminosity_setup(CCTK_ARGUMENTS);
  void zm1_luminosity_local(CCTK_ARGUMENTS);
  void zm1_luminosity_combine(CCTK_ARGUMENTS);
  void zm1_luminosity_transform(CCTK_ARGUMENTS);
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


// open file and report errors and abort if failed
#define fopen_check(fn, mode) fopen_check_fun(fn, mode, __LINE__)
static FILE *fopen_check_fun(const char *fn, const char *mode, int line);

void zm1_luminosity_setup(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;  

  CCTK_VInfo (CCTK_THORNSTRING,"Setup");
  
  // Set up the radial grid, shifted so no point at origin 
  double const drad =  rad_max / (nrad - 1);
  for (int i=0;i<nrad;i++) rad[i] = i*drad + 0.5*drad;
  if (zm1_verbose) 
    CCTK_VInfo(CCTK_THORNSTRING, "Radii: %d %e %e",nrad,drad,rad[nrad-1]);

}

// global mode
void zm1_luminosity_init(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  if (zm1_verbose)  
    CCTK_Info(CCTK_THORNSTRING, "Zeroing arrays");
  
  //CCTK_WARN(CCTK_WARN_ABORT,
  //          "Need to distribute other quantities as well, not just density");
  //assert(0);

  #pragma omp parallel for
  for (int n=0; n<nrad; ++n) {
    tot_vol[n] = 0.0;
    for (int l=0;l<nspecies;l++)
      for (int m=0;m<ngroups; m++) {
        int const idx = (nspecies*n + l)*ngroups + m; 
	Luminosity[idx] = 0.0;
      }
  }
}

// local mode
void zm1_luminosity_local(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (zm1_verbose) 
    CCTK_Info(CCTK_THORNSTRING, "Doing distributed sum arrays");
  
  // Grid cell volume
  CCTK_REAL dV = 1.0;
  for (int d=0; d<3; ++d) {
    dV *= CCTK_DELTA_SPACE(d);
  }
  
  // Grid spacing of 1d array
  CCTK_REAL const dr = rad_max / (nrad-1.0);
  
  int const nx = cctk_lsh[0]; 
  int const ny = cctk_lsh[1]; 
  int const nz = cctk_lsh[2]; 

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
  // #pragma omp parallel for reduction(+:*Luminosity)
  for (int k=zm1_ghost; k<nz-zm1_ghost; ++k) {
    for (int j=zm1_ghost; j<ny-zm1_ghost; ++j) {
      for (int i=zm1_ghost; i<nx-zm1_ghost; ++i) {
        int const ind3d = CCTK_GFINDEX3D(cctkGH, i,j,k);
        CCTK_REAL const w = weight[ind3d];

        if (w == 0.0) continue;
        
        // Radius, and index into 1d array
        CCTK_REAL const rL = r[ind3d];
        int const n = rL / dr;
        if (n >= nrad) continue; // ignore far away grid points
        
        // Calculate surface area of current cell
        CCTK_REAL const detg = det(gxx[ind3d], gxy[ind3d], gxz[ind3d],
                                   gyy[ind3d], gyz[ind3d], gzz[ind3d]);
        CCTK_REAL const sqrt_detg = sqrt(detg);
        CCTK_REAL const vol = sqrt_detg * w * dV;

	CCTK_REAL const tiny = 1.0e-20; // a tiny radius
        CCTK_REAL const a1 = z[ind3d]/rL;
	CCTK_REAL const rhoxy = sqrt(x[ind3d]*x[ind3d] + y[ind3d]*y[ind3d]); 
        //CCTK_REAL const a2 = y[ind3d]/(x[ind3d]+tiny); 
	//CCTK_REAL const rhx = sqrt(1.0-a1*a1)/sqrt(1.0+a2*a2);
	//CCTK_REAL const rhy = sqrt(1.0-a1*a1)/sqrt(1.0+a2*a2)*a2;
	CCTK_REAL const a2 = x[ind3d]/rhoxy;
	CCTK_REAL rhx = sqrt(1.0-a1*a1)*a2;
	CCTK_REAL sinphi = sqrt(1.0-a2*a2);
	if (y[ind3d]<0.0) sinphi = -sinphi;
	CCTK_REAL rhy = sqrt(1.0-a1*a1)*sinphi;
	CCTK_REAL rhz = a1;
	
	if (rhoxy <= 0.0) {
		rhx = 0.0;
		rhy = 0.0;
		rhz = 1.0;
	}

	if (rL <= 0.0) {
		rhx = 0.0;
		rhy = 0.0;
		rhz = 0.0;
	}
	
	for (int l=0;l<nspecies;l++)
	  for (int m=0;m<ngroups;m++){
            int const ig = l*ngroups + m;
	    int const index4D = CCTK_VectGFIndex3D(cctkGH,i,j,k,ig);
            int const idxr = (nspecies*n + l)*ngroups + m; 
	    
	    // This is a very kludgy way to calculate the luminosity  
	    if (l==2 && spec_idx[l]<1 || spec_idx[l] == 3) {
	      Luminosity[idxr] += vol*fnux[index4D]*rhx*rL*rL/4.0; 
	      Luminosity[idxr] += vol*fnuy[index4D]*rhy*rL*rL/4.0; 
	      Luminosity[idxr] += vol*fnuz[index4D]*rhz*rL*rL/4.0; 
              //Luminosity[idxr] += vol*fnumag[index4D];
	    } else {
	      Luminosity[idxr] += vol*fnux[index4D]*rhx*rL*rL; 
	      Luminosity[idxr] += vol*fnuy[index4D]*rhy*rL*rL; 
	      Luminosity[idxr] += vol*fnuz[index4D]*rhz*rL*rL; 
              //Luminosity[idxr] += vol*fnumag[index4D]/4.0;
	    }
#if 0
	    if(rad[n] > 19.61 && rad[n] < 19.62) {
	      CCTK_VInfo(CCTK_THORNSTRING, "Lum: %18.9E %18.9E",Luminosity[idxr],vol);
	    }
#endif
#if 0
	    // DEBUGGING; remove for performance
	    if( isnan(Luminosity[idxr]) || isinf(Luminosity[idxr])){
	      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
			 "%15.6E %15.6E %15.6E %15.6E %15.6E",
			 vol,fnux[index4D],fnuy[index4D],fnuz[index4D],rhx);
	      CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
			 "NAN in Luminosity %d %d %d r: %15.6E!!!",i,j,k,r[ind3d]);
	    }
#endif
	     
	} 

        // Add to 1d array
        tot_vol[n] += vol;
      }
    }
  }
}



// global mode
void zm1_luminosity_combine(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (zm1_verbose) 
    CCTK_Info(CCTK_THORNSTRING, "Summing arrays");
  
  int const sum = CCTK_ReductionArrayHandle("sum");
  if (sum<0) {
    CCTK_WARN(CCTK_WARN_ABORT, "'sum' reduction handle not defined");
  }
  

  
  // Sum the total volume in each radial bin 
  {
    vector<CCTK_REAL> tmp(nrad);
    int const ierr =
      CCTK_ReduceLocArrayToArray1D(cctkGH, -1, sum,
                                   tot_vol, &tmp[0], nrad, CCTK_VARIABLE_REAL);
    if (ierr<0) {
      CCTK_WARN(CCTK_WARN_ABORT, "Could not reduce tot_vol");
    }
    memcpy(tot_vol, &tmp[0], nrad*sizeof *tot_vol);
  }
  
  // Sum the Luminosity in each radial bin 
  {
    vector<CCTK_REAL> tmp(nrad*ngroups*nspecies);
    int const ierr =
      CCTK_ReduceLocArrayToArray1D(cctkGH, -1, sum, Luminosity,
                    &tmp[0], nrad*ngroups*nspecies, CCTK_VARIABLE_REAL);
    if (ierr<0) {
      CCTK_WARN(CCTK_WARN_ABORT, "Could not reduce total luminosity");
    }
    memcpy(Luminosity, &tmp[0], nrad*ngroups*nspecies*sizeof *Luminosity);
  }
  
}

// global mode
void zm1_luminosity_transform(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (zm1_verbose) 
    CCTK_Info(CCTK_THORNSTRING, "Transforming Arrays");
  for (int i=0; i<nrad; ++i) {
    
    double vol = tot_vol[i];
    for (int l=0;l<nspecies;l++)
      for (int m=0;m<ngroups;m++){
        int const idxr = (nspecies*i + l)*ngroups + m; 
        if (vol>0.0) {
          Luminosity[idxr] = Luminosity[idxr]/vol*4.0*PI;
          if (!do_M1_testing){
            // Scale to ergs/s
	    Luminosity[idxr] = Luminosity[idxr]*M_TO_L*M_TO_L*INV_PRESS_GF*C_L;
	  }
	}else{
          Luminosity[idxr] = 0.0;
        } 
    }
  } 
  if (zm1_verbose)
    CCTK_Info(CCTK_THORNSTRING, "Done transforming Arrays");

}

// global mode
namespace ZelmaniM1 {
void zm1_luminosity_output(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  int do_stuff = ( (cctk_iteration % lum_out_every) == 0 );
  if(!do_stuff) return;


  
  if (CCTK_MyProc(cctkGH)==0) {

    {
      // output the full thing
      ostringstream outdirbuf;

      outdirbuf << out_dir;
      CCTK_CreateDirectory (0755, outdirbuf.str().c_str());

      ostringstream outfilename;

      outfilename << outdirbuf.str().c_str() << "/" << "ZelmaniM1_lums.dat";

      FILE* lumoutfile=fopen_check(outfilename.str().c_str(),"a");
    
      for (int i=0; i<nrad; ++i) {
	fprintf(lumoutfile,"%10d %18.9E %18.9E ",cctk_iteration,cctk_time,rad[i]);
	for (int l=0;l<nspecies;l++) {
	  double Lumtot = 0.0;		
	  for (int m=0;m<ngroups;m++){
	    int const idxr = (nspecies*i + l)*ngroups + m;
	    Lumtot += Luminosity[idxr]; 
	  }
	  fprintf(lumoutfile,"%18.9E ",Lumtot);
	}
	fprintf(lumoutfile,"\n");
      }

      fprintf(lumoutfile,"\n");
      fprintf(lumoutfile,"\n");

      fclose(lumoutfile);
    } // done with full output

    // output just at a few radii:
    {
      int in = 0;
      while(lum_out_radii[in]>0.0 && lum_out_radii[in] < rad[nrad-1]) {
	ostringstream outdirbuf;
	outdirbuf << out_dir;
	CCTK_CreateDirectory (0755, outdirbuf.str().c_str());

	ostringstream outfilename,sumoutfilename;
	char tmpstring[20];
	sprintf(tmpstring,"%08.3f",lum_out_radii[in]);
	outfilename << outdirbuf.str().c_str() << "/" 
		    << "ZelmaniM1_lums_r" << tmpstring << ".dat";
	sumoutfilename << outdirbuf.str().c_str() << "/" 
		    << "ZelmaniM1_sum_lums_r" << tmpstring << ".dat";

	double dr=rad[1]-rad[0];
	int ir = lum_out_radii[in]/dr;
        
	if (zm1_verbose) 	
	  CCTK_VInfo(CCTK_THORNSTRING, "Radius: %d %d %e %e",in,ir,dr,lum_out_radii[in]);
         
	// Output Spectra 	
	FILE* lumoutfile=fopen_check(outfilename.str().c_str(),"a");
	for(int m=0;m<ngroups;m++) {
	  fprintf(lumoutfile,"%10d %18.9E %18.9E ",cctk_iteration,cctk_time,neutrino_energies[sgroup+m]);
	  for(int l=0;l<nspecies;l++) {
	    int const idxr = (nspecies*ir + l)*ngroups + m; 
	    fprintf(lumoutfile,"%18.9E ",Luminosity[idxr]/bin_widths[sgroup+m]);	    
	  }
	  fprintf(lumoutfile,"\n");
	}
	fprintf(lumoutfile,"\n");
	fclose(lumoutfile);
	
	// Output luminosities
	FILE* sumlumoutfile=fopen_check(sumoutfilename.str().c_str(),"a");
	fprintf(sumlumoutfile,"%10d %18.9E ",cctk_iteration,cctk_time);
	
	for(int l=0;l<nspecies;l++) {
	  double Lumtot = 0.0;
	  for(int m=0;m<ngroups;m++) {
	    int const idxr = (nspecies*ir + l)*ngroups + m; 
	    Lumtot += Luminosity[idxr];
	  }
	  fprintf(sumlumoutfile,"%18.9E ",Lumtot);	    

	}
	
	for(int l=0;l<nspecies;l++) {
	  double Lumtot = 0.0;
	  double eavgden = 0.0;
	  for(int m=0;m<ngroups;m++) {
	    int const idxr = (nspecies*ir + l)*ngroups + m; 
	    Lumtot += Luminosity[idxr];
	    eavgden += Luminosity[idxr]/neutrino_energies[sgroup+m];
	  }
	  fprintf(sumlumoutfile,"%18.9E ",Lumtot/eavgden);	    

	}

	fprintf(sumlumoutfile,"\n");
	fclose(sumlumoutfile);
	in++;
      }

    }



  }
}
}

static FILE *fopen_check_fun(const char *fn, const char *mode, int line)
{
  FILE *fh = fopen(fn, mode);
  if (! fh) {
    CCTK_VWarn(CCTK_WARN_ABORT, line, __FILE__, CCTK_THORNSTRING,
               "Could not open file '%s': %s", fn, strerror(errno));
  }
  return fh;
}
