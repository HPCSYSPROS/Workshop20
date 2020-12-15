#include <iostream>
#include <algorithm>
#define H5_USE_16_API 1

#include "hdf5.h"

#define DEFINE_GLOBALS 1
#include "ZelmaniM1.hh"
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "mpi.h"
#include <sys/stat.h>

#define INV_R_GF 1.4768e5
#define INV_PRESS_GF 5.56082181777535e38
#define INV_EPS_GF 8.98755175549085e20
#define EPS_GF 1.11265006e-21
#define C_L 2.99792458e10
#define G_GRAV 6.673e-8
#define M_SUN 1.98892e33
#define PI 3.14159265359

using namespace std;


// mini NoMPI
#ifdef HAVE_CAPABILITY_MPI
#include <mpi.h>
#define BCAST(buffer, size) MPI_Bcast(buffer, size, MPI_BYTE, my_reader_process, MPI_COMM_WORLD)
#else
#define BCAST(buffer, size) do { /* do nothing */ } while(0)
#endif

// Catch HDF5 errors
#define HDF5_ERROR(fn_call)                                              \
  if(doIO) {                                                             \
    int _error_code = fn_call;                                           \
    if (_error_code < 0) {                                               \
      CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,                  \
                  "HDF5 call '%s' returned error code %d",               \
                  #fn_call, _error_code);                                \
    }                                                                    \
  }


namespace ZelmaniM1 {

  extern "C"
  void zm1_readtable(CCTK_ARGUMENTS) {

    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;
   
    int me,mp;
    MPI_Comm_rank(MPI_COMM_WORLD,&me);
    MPI_Comm_size(MPI_COMM_WORLD,&mp);

    CCTK_Info(CCTK_THORNSTRING,"*******************************");
    CCTK_Info(CCTK_THORNSTRING,"Reading NuLib table file:");
    CCTK_Info(CCTK_THORNSTRING,nulib_table_name);

    CCTK_INT my_reader_process = reader_process;
    if (my_reader_process < 0 || my_reader_process >= CCTK_nProcs(cctkGH))
    {
      CCTK_VWarn(CCTK_WARN_COMPLAIN, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Requested IO process %d out of range. Reverting to process 0.", my_reader_process);
      my_reader_process = 0;
    }
    const int doIO = !read_table_on_single_process || CCTK_MyProc(cctkGH) == my_reader_process;

    // check that file is readable
    if (false) {
      FILE* fp = NULL;
      fp = fopen(nulib_table_name,"r");
      if(fp != NULL) {
        fclose(fp);
      } else {
        CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
       "Could not read nulib_table_name'%s'",
       nulib_table_name);
      }
    } else {
      struct stat fstat;
      if (doIO && stat(nulib_table_name,&fstat)==-1) {
        CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Could not read nulib_table_name'%s'",
                   nulib_table_name);
      }
    }
     
    hid_t file;
    HDF5_ERROR(file = H5Fopen(nulib_table_name, H5F_ACC_RDONLY, H5P_DEFAULT));

    // crazy macro definition for an inline call to some shit
    // to read in data
#define READ_SHIT_HDF5(NAME,VAR,TYPE,NELMS)                                     \
    do {                                                                        \
      hid_t dataset;                                                            \
      HDF5_ERROR(dataset = H5Dopen(file, NAME));                                \
      HDF5_ERROR(H5Dread(dataset, TYPE, H5S_ALL, H5S_ALL, H5P_DEFAULT, VAR));   \
      HDF5_ERROR(H5Dclose(dataset));                                            \
      if (read_table_on_single_process)                                         \
        BCAST (VAR, sizeof(*(VAR))*(NELMS));                                    \
    } while (0)

    READ_SHIT_HDF5("nrho",&nrho, H5T_NATIVE_INT, 1);
    READ_SHIT_HDF5("ntemp",&ntemp, H5T_NATIVE_INT, 1);
    READ_SHIT_HDF5("nye",&nye, H5T_NATIVE_INT, 1);

    if(do_gray && do_opac) {
      CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Can't do spectral and gray transport together, bailing!");
    }

    int ng;
    int ns;
    READ_SHIT_HDF5("number_groups",&ng, H5T_NATIVE_INT, 1);
    READ_SHIT_HDF5("number_species",&ns, H5T_NATIVE_INT, 1);
    if (do_gray){
      ngtab = ng;
    } else {
      ngtab = ngroups;
    }

    CCTK_VInfo(CCTK_THORNSTRING,"NuLib:     nrho: %d",nrho);
    CCTK_VInfo(CCTK_THORNSTRING,"NuLib:    ntemp: %d",ntemp);
    CCTK_VInfo(CCTK_THORNSTRING,"NuLib:      nye: %d",nye);
    CCTK_VInfo(CCTK_THORNSTRING,"NuLib:  ngroups: %d",ng);
    CCTK_VInfo(CCTK_THORNSTRING,"NuLib: nspecies: %d",ns);


    if(ng < ngroups && do_opac) {
      CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Number of groups in NuLib table %d is too few for parameter ngtab %d",
                 ng,ngroups);
    }

    if(ns < nspecies) {
      CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Number of species in NuLib table %d is too few for parameter nspecies %d",
                 ns,nspecies);
    }

    rho_points = (double*) malloc(sizeof(double)*nrho);
    temp_points = (double*) malloc(sizeof(double)*ntemp);
    ye_points = (double*) malloc(sizeof(double)*nye);

    READ_SHIT_HDF5("rho_points",rho_points, H5T_NATIVE_DOUBLE, nrho);
    READ_SHIT_HDF5("temp_points",temp_points, H5T_NATIVE_DOUBLE, ntemp);
    READ_SHIT_HDF5("ye_points",ye_points, H5T_NATIVE_DOUBLE, nye);
    
    // Rescale temperature and density grid to log spacin 
    for (int m=0;m<nrho; m++) rho_points[m]  = log10(rho_points[m] );
    for (int m=0;m<ntemp;m++) temp_points[m] = log10(temp_points[m]);
    
    bin_bottom = (double*) malloc(sizeof(double)*ng);
    bin_top = (double*) malloc(sizeof(double)*ng);
    bin_widths = (double*) malloc(sizeof(double)*ng);
    neutrino_energies = (double*) malloc(sizeof(double)*ng);

    READ_SHIT_HDF5("bin_bottom",bin_bottom, H5T_NATIVE_DOUBLE, ng);
    READ_SHIT_HDF5("bin_top",bin_top, H5T_NATIVE_DOUBLE, ng);
    READ_SHIT_HDF5("bin_widths",bin_widths, H5T_NATIVE_DOUBLE, ng);
    READ_SHIT_HDF5("neutrino_energies",neutrino_energies, H5T_NATIVE_DOUBLE, ng);

    double *absorb_temp, *emis_temp, *scat_temp;
    const size_t table_size = ng*ns*nrho*ntemp*nye;
    absorb_temp = (double*) malloc(sizeof(double)*table_size);
    emis_temp = (double*) malloc(sizeof(double)*table_size);
    scat_temp = (double*) malloc(sizeof(double)*table_size);
    
    READ_SHIT_HDF5("absorption_opacity",absorb_temp,H5T_NATIVE_DOUBLE,table_size);
    READ_SHIT_HDF5("scattering_opacity",scat_temp,H5T_NATIVE_DOUBLE,table_size);
    READ_SHIT_HDF5("emissivities",emis_temp,H5T_NATIVE_DOUBLE,table_size);

    // mangle all these tables into one big table
    // indices in order of variation speed (fastest first)
    // igroup, ispecies, ikind, irho, itemp, iye

    ntables = 3*ngtab*nspecies;
    alltables = (double *) malloc(sizeof(double)*ntables*nrho*ntemp*nye);

    // tables read by EOS are in the following order (fastest index first)
    // irho, itemp, iye, ispecies, igroup
    double TINY = 1.e-90;
    
    
    for(int irho=0;irho < nrho;irho++)
      for(int itemp=0;itemp < ntemp;itemp++)
        for(int iye=0;iye < nye;iye++)
          for(int is=0;is < nspecies;is++)
            for(int ig=0;ig < ngtab;ig++) {
                
              int isold = is;
              if (spec_idx[is]>0) isold = spec_idx[is]-1;
              if (isold>=ns) CCTK_WARN(0,"Bad species index.");
              if (sgroup+ngtab>ng) CCTK_WARN(0,"Bad number of groups.");
              int indold = irho + nrho*(itemp + ntemp*(iye + nye*(isold + ns*(ig+sgroup))));
              // absorb is ikind = 1
              // emis is ikind = 2
              // scat is ikind = 3
              int nkind = 3;
        
              // absorb
              int ikind = 0;
              int indnew = ig + ngtab*
                (is + nspecies*(ikind + nkind*(irho + nrho*(itemp + ntemp*iye))));
                alltables[indnew] = log10(absorb_temp[indold]*INV_R_GF + TINY); // In geometrized units
        
              // emis
              ikind = 1;
              indnew = ig + ngtab*
                (is + nspecies*(ikind + nkind*(irho + nrho*(itemp + ntemp*iye))));
                alltables[indnew] = log10(4.0*PI*emis_temp[indold]*bin_widths[ig]*INV_R_GF/INV_PRESS_GF/C_L + TINY); // convert from erg/cc/MeV/s/Sr to geometrized units [ENERGY/L^4] 
        
              // scat
              ikind = 2;
              indnew = ig + ngtab*
                (is + nspecies*(ikind + nkind*(irho + nrho*(itemp + ntemp*iye))));
                alltables[indnew] = log10(scat_temp[indold]*INV_R_GF + TINY); // In geometrized units
        
            }

    free(absorb_temp);
    free(emis_temp);
    free(scat_temp);
    
    HDF5_ERROR(H5Fclose(file));


    CCTK_Info(CCTK_THORNSTRING,"Done preparing opacity tables!");
    CCTK_Info(CCTK_THORNSTRING,"*******************************");
    return;
  }
}

