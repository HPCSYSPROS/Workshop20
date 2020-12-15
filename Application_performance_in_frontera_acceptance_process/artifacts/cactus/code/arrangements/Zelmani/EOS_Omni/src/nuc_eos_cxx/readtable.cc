#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define H5_USE_16_API 1
#include "hdf5.h"
#include "nuc_eos.hh"

#include <cctk.h>
#include <cctk_Parameters.h>

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
    if (_error_code < 0) {	       				         \
      CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,                  \
                  "HDF5 call '%s' returned error code %d",               \
                  #fn_call, _error_code);                                \
    }                                                                    \
  }

static int file_is_readable(const char* filename);
static int file_is_readable(const char* filename)
{
    FILE* fp = NULL;
    fp = fopen(filename, "r");
    if(fp != NULL)
    {
        fclose(fp);
        return 1;
    }
    return 0;
}


// define the variables
namespace nuc_eos {
  double temp0, temp1;
  double energy_shift;

  double eos_rhomax, eos_rhomin;
  double eos_tempmin, eos_tempmax;
  double eos_yemin, eos_yemax;
  
  double c2p_tempmin;
  double c2p_tempmax;

}
namespace nuc_eos_private {
  int nrho;
  int ntemp;
  int nye;

  double * restrict alltables;
  double * restrict epstable;
  double * restrict logrho;
  double * restrict logtemp;
  double dlintemp, dlintempi;
  double drholintempi;
  double dlintempyei;
  double drholintempyei;
  double * restrict yes;
  double dtemp, dtempi;
  double drho, drhoi;
  double dye, dyei;
  double drhotempi;
  double drhoyei;
  double dtempyei;
  double drhotempyei;
}



namespace nuc_eos {
// TODO: replace with version in ET EOS_Omni. NOTE: table arrangement changed.

// Cactus calls this function. It reads in the table and calls a fortran
// function to setup values for the fortran eos module
extern "C"
void nuc_eos_C_ReadTable(char* , const cGH* cctkGH)
{
  DECLARE_CCTK_PARAMETERS;

  using namespace nuc_eos;
  using namespace nuc_eos_private;

  CCTK_VInfo(CCTK_THORNSTRING,"*******************************");
  CCTK_VInfo(CCTK_THORNSTRING,"Reading nuc_eos table file:");
  CCTK_VInfo(CCTK_THORNSTRING,"%s",nuceos_table_name);
  CCTK_VInfo(CCTK_THORNSTRING,"*******************************");

  CCTK_INT my_reader_process = reader_process;
  if (my_reader_process < 0 || my_reader_process >= CCTK_nProcs(cctkGH))
  {
    CCTK_VWarn(CCTK_WARN_COMPLAIN, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Requested IO process %d out of range. Reverting to process 0.", my_reader_process);
    my_reader_process = 0;
  }
  const int doIO = !read_table_on_single_process || CCTK_MyProc(cctkGH) == my_reader_process;

  hid_t file;
  if (doIO && !file_is_readable(nuceos_table_name)) {
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
               "Could not read nuceos_table_name %s \n",
               nuceos_table_name);
  }
  HDF5_ERROR(file = H5Fopen(nuceos_table_name, H5F_ACC_RDONLY, H5P_DEFAULT));

// Use these two defines to easily read in a lot of variables in the same way
// The first reads in one variable of a given type completely
#define READ_BCAST_EOS_HDF5(NAME,VAR,TYPE,MEM,NELEMS)                                      \
  do {                                                                        \
    hid_t dataset;                                                            \
    HDF5_ERROR(dataset = H5Dopen(file, NAME));                                \
    HDF5_ERROR(H5Dread(dataset, TYPE, MEM, H5S_ALL, H5P_DEFAULT, VAR));       \
    if (read_table_on_single_process)                                         \
      BCAST (VAR, sizeof(*(VAR))*(NELEMS));                                   \
    HDF5_ERROR(H5Dclose(dataset));                                            \
  } while (0)
// The second reads a given variable into a hyperslab of the alltables_temp array
#define READ_BCAST_EOSTABLE_HDF5(NAME,OFF,DIMS)                          \
  do {                                                                   \
    READ_BCAST_EOS_HDF5(NAME,&alltables_temp[(OFF)*(DIMS)[1]],H5T_NATIVE_DOUBLE,H5S_ALL,(DIMS)[1]); \
  } while (0)

  // Read size of tables
  READ_BCAST_EOS_HDF5("pointsrho",  &nrho,  H5T_NATIVE_INT, H5S_ALL, 1);
  READ_BCAST_EOS_HDF5("pointstemp", &ntemp, H5T_NATIVE_INT, H5S_ALL, 1);
  READ_BCAST_EOS_HDF5("pointsye",   &nye,   H5T_NATIVE_INT, H5S_ALL, 1);


  // Allocate memory for tables
  double* alltables_temp;
  if (!(alltables_temp = (double*)malloc(nrho * ntemp * nye * NTABLES * sizeof(double)))) {
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                "Cannot allocate memory for EOS table");
  }
  if (!(logrho = (double*)malloc(nrho * sizeof(double)))) {
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                "Cannot allocate memory for EOS table");
  }
  if (!(logtemp = (double*)malloc(ntemp * sizeof(double)))) {
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                "Cannot allocate memory for EOS table");
  }
  if (!(yes = (double*)malloc(nye * sizeof(double)))) {
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                "Cannot allocate memory for EOS table");
  }

  // Prepare HDF5 to read hyperslabs into alltables_temp
  hsize_t table_dims[2] = {NTABLES, (hsize_t)nrho * ntemp * nye};
  hsize_t var3[2]       = { 1, (hsize_t)nrho * ntemp * nye};
  hid_t mem3 =  H5Screate_simple(2, table_dims, NULL);

  // Read alltables_temp
  READ_BCAST_EOSTABLE_HDF5("logpress",  0, table_dims);
  READ_BCAST_EOSTABLE_HDF5("logenergy", 1, table_dims);
  READ_BCAST_EOSTABLE_HDF5("entropy",   2, table_dims);
  READ_BCAST_EOSTABLE_HDF5("munu",      3, table_dims);
  READ_BCAST_EOSTABLE_HDF5("cs2",       4, table_dims);
  READ_BCAST_EOSTABLE_HDF5("dedt",      5, table_dims);
  READ_BCAST_EOSTABLE_HDF5("dpdrhoe",   6, table_dims);
  READ_BCAST_EOSTABLE_HDF5("dpderho",   7, table_dims);
  // chemical potentials
  READ_BCAST_EOSTABLE_HDF5("muhat",     8, table_dims);
  READ_BCAST_EOSTABLE_HDF5("mu_e",      9, table_dims);
  READ_BCAST_EOSTABLE_HDF5("mu_p",     10, table_dims);
  READ_BCAST_EOSTABLE_HDF5("mu_n",     11, table_dims);
  // compositions
  READ_BCAST_EOSTABLE_HDF5("Xa",       12, table_dims);
  READ_BCAST_EOSTABLE_HDF5("Xh",       13, table_dims);
  READ_BCAST_EOSTABLE_HDF5("Xn",       14, table_dims);
  READ_BCAST_EOSTABLE_HDF5("Xp",       15, table_dims);
  // average nucleus
  READ_BCAST_EOSTABLE_HDF5("Abar",     16, table_dims);
  READ_BCAST_EOSTABLE_HDF5("Zbar",     17, table_dims);
  // Gamma
  READ_BCAST_EOSTABLE_HDF5("gamma",    18, table_dims);

  // Read additional tables and variables
  READ_BCAST_EOS_HDF5("logrho",       logrho,        H5T_NATIVE_DOUBLE, H5S_ALL, nrho);
  READ_BCAST_EOS_HDF5("logtemp",      logtemp,       H5T_NATIVE_DOUBLE, H5S_ALL, ntemp);
  READ_BCAST_EOS_HDF5("ye",           yes,            H5T_NATIVE_DOUBLE, H5S_ALL, nye);
  READ_BCAST_EOS_HDF5("energy_shift", &energy_shift, H5T_NATIVE_DOUBLE, H5S_ALL, 1);

  HDF5_ERROR(H5Sclose(mem3));
  HDF5_ERROR(H5Fclose(file));

  // change ordering of alltables array so that
  // the table kind is the fastest changing index
  if (!(alltables = (double*)malloc(nrho * ntemp * nye * NTABLES 
				    * sizeof(double)))) {
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                "Cannot allocate memory for EOS table");
  }
  for(int iv = 0;iv<NTABLES;iv++) 
    for(int k = 0; k<nye;k++) 
      for(int j = 0; j<ntemp; j++) 
	for(int i = 0; i<nrho; i++) {
	  int indold = i + nrho*(j + ntemp*(k + nye*iv));
	  int indnew = iv + NTABLES*(i + nrho*(j + ntemp*k));
	  alltables[indnew] = alltables_temp[indold];
	}

  // free memory of temporary array
  free(alltables_temp);

  // convert units, convert logs to natural log
  // The latter is great, because exp() is way faster than pow()
  // pressure
  energy_shift = energy_shift * EPSGF;
  for(int i=0;i<nrho;i++) {
    // rewrite:
    //logrho[i] = log(pow(10.0,logrho[i]) * RHOGF);
    // by using log(a^b*c) = b*log(a)+log(c)
    logrho[i] = logrho[i] * log(10.) + log(RHOGF);
  }

  for(int i=0;i<ntemp;i++) {
    //logtemp[i] = log(pow(10.0,logtemp[i]));
    logtemp[i] = logtemp[i]*log(10.0);
  }

  // allocate epstable; a linear-scale eps table
  // that allows us to extrapolate to negative eps
  if (!(epstable = (double*)malloc(nrho * ntemp * nye  
				    * sizeof(double)))) {
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                "Cannot allocate memory for eps table\n");
  }

  // convert units
  for(int i=0;i<nrho*ntemp*nye;i++) {

    { // pressure
      int idx = 0 + NTABLES*i;
      alltables[idx] = alltables[idx] * log(10.0) + log(PRESSGF);
    }

    { // eps
      int idx = 1 + NTABLES*i;
      alltables[idx] = alltables[idx] * log(10.0) + log(EPSGF);
      epstable[i] = exp(alltables[idx]);
    }

    { // cs2
      int idx = 4 + NTABLES*i;
      alltables[idx] *= LENGTHGF*LENGTHGF/TIMEGF/TIMEGF;
    }

    { // dedT
      int idx = 5 + NTABLES*i;
      alltables[idx] *= EPSGF;
    }

    { // dpdrhoe
      int idx = 6 + NTABLES*i;
      alltables[idx] *= PRESSGF/RHOGF;
    }

    { // dpderho
      int idx = 7 + NTABLES*i;
      alltables[idx] *= PRESSGF/EPSGF;
    }

  }

  temp0 = exp(logtemp[0]);
  temp1 = exp(logtemp[1]);

  // set up some vars
  dtemp  = (logtemp[ntemp-1] - logtemp[0]) / (1.0*(ntemp-1));
  dtempi = 1.0/dtemp;

  dlintemp = temp1-temp0;
  dlintempi = 1.0/dlintemp;

  drho  = (logrho[nrho-1] - logrho[0]) / (1.0*(nrho-1));
  drhoi = 1.0/drho;

  dye  = (yes[nye-1] - yes[0]) / (1.0*(nye-1));
  dyei = 1.0/dye;

  drhotempi   = drhoi * dtempi;
  drholintempi = drhoi * dlintempi;
  drhoyei     = drhoi * dyei;
  dtempyei    = dtempi * dyei;
  dlintempyei = dlintempi * dyei;
  drhotempyei = drhoi * dtempi * dyei;
  drholintempyei = drhoi * dlintempi * dyei;

  eos_rhomax = exp(logrho[nrho-1]);
  eos_rhomin = exp(logrho[0]);
  
  eos_tempmax = exp(logtemp[ntemp-1]);
  eos_tempmin = exp(logtemp[0]);

  eos_yemax = yes[nye-1];
  eos_yemin = yes[0];

}

extern "C"
void CCTK_FNAME(nuc_eos_c_get_energy_shift)(double *energy_shift_fortran,
					    double *eos_tempmin_fortran,
					    double *eos_tempmax_fortran,
					    double *eos_yemin_fortran,
					    double * eos_yemax_fortran) {

  *energy_shift_fortran = energy_shift;
  *eos_tempmin_fortran = eos_tempmin;
  *eos_tempmax_fortran = eos_tempmax;
  *eos_yemin_fortran = eos_yemin;
  *eos_yemax_fortran = eos_yemax;

}

} // namespace nuc_eos


