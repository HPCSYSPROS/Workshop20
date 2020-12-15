#ifndef CARPETIOHDF5_HH
#define CARPETIOHDF5_HH

#define H5_USE_16_API 1
#include <hdf5.h>

#include <vector>

#include "cctk_Arguments.h"
#include "CactusBase/IOUtil/src/ioutil_Utils.h"
#include "carpet.hh"

// some macros for HDF5 group names
#define METADATA_GROUP "Parameters and Global Attributes"
#define ALL_PARAMETERS "All Parameters"
#define GRID_STRUCTURE "Grid Structure v5"

// atomic HDF5 datatypes for the generic CCTK datatypes
// (the one for CCTK_COMPLEX is created at startup as a compound HDF5 datatype)
#define HDF5_CHAR H5T_NATIVE_CHAR

#ifdef CCTK_REAL_PRECISION_16
#define HDF5_REAL H5T_NATIVE_LDOUBLE
#elif CCTK_REAL_PRECISION_8
#define HDF5_REAL H5T_NATIVE_DOUBLE
#elif CCTK_REAL_PRECISION_4
#define HDF5_REAL H5T_NATIVE_FLOAT
#endif

#ifdef CCTK_INTEGER_PRECISION_8
#define HDF5_INT H5T_NATIVE_LLONG
#elif CCTK_INTEGER_PRECISION_4
#define HDF5_INT H5T_NATIVE_INT
#elif CCTK_INTEGER_PRECISION_2
#define HDF5_INT H5T_NATIVE_SHORT
#elif CCTK_INTEGER_PRECISION_1
#define HDF5_INT H5T_NATIVE_CHAR
#endif

// check return code of HDF5 call and print a warning in case of an error
#define HDF5_ERROR(fn_call)                                                    \
  do {                                                                         \
    int _error_code = fn_call;                                                 \
                                                                               \
    if (_error_code < 0) {                                                     \
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,                      \
                 "HDF5 call '%s' returned error code %d", #fn_call,            \
                 _error_code);                                                 \
      error_count++;                                                           \
    }                                                                          \
  } while (0)

// datatype of the start[] parameter in a call to H5Sselect_hyperslab()
// (the HDF5 API has changed in this respect after version 1.6.3)
#if (H5_VERS_MAJOR == 1 &&                                                     \
     (H5_VERS_MINOR < 6 || (H5_VERS_MINOR == 6 && H5_VERS_RELEASE < 4)))
#define slice_start_size_t hssize_t
#else
#define slice_start_size_t hsize_t
#endif

// CarpetIOHDF5 GH extension structure
typedef struct {
  // default number of times to output
  int out_every_default;

  // list of variables to output
  char *out_vars;

  // stop on I/O parameter parsing errors ?
  int stop_on_parse_errors;

  // I/O request description list (for all variables)
  vector<ioRequest *> requests;

  // directory in which to output
  char *out_dir;

  // ring buffer for list of successfully created cp files
  int checkpoint_keep;
  int cp_filename_index;
  char **cp_filename_list;

  // list of recovery files to remove
  int recovery_num_filenames;
  char **recovery_filename_list;

  // iteration number of the last checkpoint
  int last_checkpoint_iteration;

  // hdf5 datatype for complex variables; to be set at run time
  hid_t HDF5_COMPLEX, HDF5_COMPLEX8, HDF5_COMPLEX16, HDF5_COMPLEX32;

} CarpetIOHDF5GH;

namespace CarpetIOHDF5 {
// callback routine registered for recovery/filereader
int Recover(cGH *cctkGH, const char *basefilename, int called_from);

// worker routines to write a single variable
int WriteVarUnchunked(const cGH *const cctkGH, hid_t file, CCTK_REAL &io_bytes,
                      const ioRequest *const request,
                      bool called_from_checkpoint);
int WriteVarChunkedSequential(const cGH *const cctkGH, hid_t file,
                              CCTK_REAL &io_bytes,
                              const ioRequest *const request,
                              bool called_from_checkpoint, hid_t index = -1);
int WriteVarChunkedParallel(const cGH *const cctkGH, hid_t file,
                            CCTK_REAL &io_bytes, const ioRequest *const request,
                            bool called_from_checkpoint, hid_t index = -1);

int WriteMetadata(const cGH *const cctkGH, int const nioprocs,
                  int const firstvar, int const numvars,
                  bool const called_from_checkpoint, hid_t const file);

int AddSliceAttributes(const cGH *const cctkGH, const char *const fullname,
                       const int refinementlevel, const int multigridlevel,
                       const int map, const int timelevel,
                       const vector<double> &origin,
                       const vector<double> &delta, const vector<int> &iorigin,
                       const vector<int> &ioffset,
                       const vector<int> &ioffsetdenom, const vector<int> &bbox,
                       const vector<int> &nghostzones, const string &active,
                       hid_t &dataset, const vector<hsize_t> &shape,
                       const bool is_index);

int WriteAttribute(hid_t const group, char const *const name, int const ivalue);
int WriteAttribute(hid_t const group, char const *const name,
                   double const dvalue);
int WriteAttribute(hid_t const group, char const *const name,
                   char const *const svalue);
int WriteAttribute(hid_t const group, char const *const name,
                   int const *const ivalues, int const nvalues);
int WriteAttribute(hid_t const group, char const *const name,
                   double const *const dvalues, int const nvalues);
int WriteAttribute(hid_t const group, char const *const name,
                   char const *const *const svalues, int const nvalues);
int WriteAttribute(hid_t const group, char const *const name,
                   hsize_t const *const svalues, int const nvalues);
int WriteLargeAttribute(hid_t const group, char const *const name,
                        char const *const svalue);

// returns an HDF5 datatype corresponding to the given CCTK datatype
hid_t CCTKtoHDF5_Datatype(const cGH *const cctkGH, int cctk_type,
                          bool single_precision);

// Everything is a class template, so that it can easily be
// instantiated for all output dimensions

// computes the active region the way dh::regrid does
void GetAllActive(const dh *dd, const gh *hh, int ml, int rl, ibset &allactive);

template <int outdim> struct IOHDF5 {

  // name of the output directory
  static char *my_out_slice_dir;

  // list of variables to output
  static char *my_out_slice_vars;

  // I/O request description list (for all variables)
  static vector<ioRequest *> slice_requests;

  // number of I/O processors to use, which processor does I/O for me
  static int nioprocs;
  static int ioproc;
  static int ioproc_every;

  // Scheduled functions
  static int Startup();

  // Registered functions
  static void *SetupGH(tFleshConfig *fc, int convLevel, cGH *cctkGH);

  static int OutputGH(const cGH *cctkGH);
  static int OutputVarAs(const cGH *cctkGH, const char *varname,
                         const char *alias);
  static int TimeToOutput(const cGH *cctkGH, int vindex);
  static int TriggerOutput(const cGH *cctkGH, int vindex);

  // Other functions
  static void CheckSteerableParameters(const cGH *cctkGH);

  static bool DidOutput(const cGH *cctkGH, int vindex, string basefilename,
                        bool &is_new_file, bool &truncate_file);

  static bool DirectionIsRequested(const vect<int, outdim> &dirs);

  static void OutputDirection(const cGH *cctkGH, int vindex, string alias,
                              string basefilename,
                              const vect<int, outdim> &dirs, bool is_new_file,
                              bool truncate_file);

  static int OpenFile(const cGH *cctkGH, int m, int vindex, int numvars,
                      string alias, string basefilename,
                      const vect<int, outdim> &dirs, bool is_new_file,
                      bool truncate_file, hid_t &file, hid_t &index_file);

  static int WriteHDF5(const cGH *cctkGH, hid_t &file, hid_t &index_file,
                       vector<gdata *> const gfdatas,
                       const bbox<int, dim> &gfext, const int vi,
                       const vect<int, dim> &org, const vect<int, outdim> &dirs,
                       const int rl, const int ml, const int m, const int c,
                       const int output_component, const int tl,
                       const CCTK_REAL coord_time,
                       const vect<CCTK_REAL, dim> &coord_lower,
                       const vect<CCTK_REAL, dim> &coord_upper);

  static int CloseFile(const cGH *cctkGH, hid_t &file, hid_t &index_file);

  static ivect GetOutputOffset(const cGH *cctkGH, int m,
                               const vect<int, outdim> &dirs);
  static int IOProcForProc(int proc);

}; // struct IOHDF5

// scheduled and aliased routines (must be declared as C according
// to schedule.ccl)
extern "C" {

int CarpetIOHDF5_RecoverParameters(void);
int CarpetIOHDF5_SetNumRefinementLevels(void);
int CarpetIOHDF5_Startup(void);
void CarpetIOHDF5_Init(CCTK_ARGUMENTS);
void CarpetIOHDF5_InitCheckpointingIntervals(CCTK_ARGUMENTS);
void CarpetIOHDF5_RecoverGridStructure(CCTK_ARGUMENTS);
void CarpetIOHDF5_CloseFiles(CCTK_ARGUMENTS);
void CarpetIOHDF5_InitialDataCheckpoint(CCTK_ARGUMENTS);
void CarpetIOHDF5_EvolutionCheckpoint(CCTK_ARGUMENTS);
void CarpetIOHDF5_TerminationCheckpoint(CCTK_ARGUMENTS);

CCTK_INT CarpetIOHDF5_SetCheckpointGroups(CCTK_INT const *groups,
                                          CCTK_INT ngroups);

} // extern "C"

// Which groups should be checkpointed. If empty, all variables
// should be checkpointed (default).
extern vector<bool> groups_to_checkpoint;

} // namespace CarpetIOHDF5

#endif // !defined(CARPETIOHDF5_HH)
