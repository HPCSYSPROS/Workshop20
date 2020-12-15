#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <cctk_Version.h>
#include <util_Table.h>

#include "CactusBase/IOUtil/src/ioGH.h"
#include "CactusBase/IOUtil/src/ioutil_CheckpointRecovery.h"

#include <Timer.hh>

#include "CarpetIOHDF5.hh"

#include "defs.hh"
#include "gh.hh"

namespace CarpetIOHDF5 {

using namespace std;
using namespace Carpet;

// Variable definitions

// when was the last checkpoint written ?
static int last_checkpoint_iteration = -1;
static CCTK_REAL last_checkpoint_walltime;

// registered GH extension setup routine
static void *SetupGH(tFleshConfig *const fleshconfig, const int convLevel,
                     cGH *const cctkGH);

// callbacks for CarpetIOHDF5's I/O method
static int OutputGH(const cGH *const cctkGH);
static int OutputVarAs(const cGH *const cctkGH, const char *const varname,
                       const char *const alias);
static int TimeToOutput(const cGH *const cctkGH, const int vindex);
static int TriggerOutput(const cGH *const cctkGH, const int vindex);

// general checkpoint routine
static void Checkpoint(const cGH *const cctkGH, int called_from);

// callback for I/O parameter parsing routine
static void GetVarIndex(int vindex, const char *optstring, void *arg);

static void CheckSteerableParameters(const cGH *const cctkGH,
                                     CarpetIOHDF5GH *myGH);

//////////////////////////////////////////////////////////////////////////////
// public routines
//////////////////////////////////////////////////////////////////////////////

int CarpetIOHDF5_Startup(void) {
  CCTK_RegisterBanner("AMR HDF5 I/O provided by CarpetIOHDF5");

  const int GHExtension = CCTK_RegisterGHExtension(CCTK_THORNSTRING);
  CCTK_RegisterGHExtensionSetupGH(GHExtension, SetupGH);

  IOHDF5<0>::Startup();
  IOHDF5<1>::Startup();
  IOHDF5<2>::Startup();
  IOHDF5<3>::Startup();

  return (0);
}

// Called at basegrid during regular startup
void CarpetIOHDF5_Init(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;

  *this_iteration = -1;
  *next_output_iteration = 0;
  *next_output_time = cctk_time;

  for (int d = 0; d < 4; ++d) {
    this_iteration_slice[d] = 0;
    last_output_iteration_slice[d] = 0;
    last_output_time_slice[d] = cctk_time;
  }

  last_checkpoint_iteration = cctk_iteration;
  last_checkpoint_walltime = CCTK_RunTime() / 3600.0;
}

// Called after recovering
void CarpetIOHDF5_InitCheckpointingIntervals(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;

  last_checkpoint_iteration = cctk_iteration;
  last_checkpoint_walltime = CCTK_RunTime() / 3600.0;
}

void CarpetIOHDF5_InitialDataCheckpoint(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (checkpoint and checkpoint_ID) {
    if (not CCTK_Equals(verbose, "none")) {
      CCTK_INFO("---------------------------------------------------------");
      CCTK_VInfo(CCTK_THORNSTRING, "Dumping initial checkpoint at "
                                   "iteration %d, simulation time %g",
                 cctk_iteration, double(cctk_time));
      CCTK_INFO("---------------------------------------------------------");
    }
    Checkpoint(cctkGH, CP_INITIAL_DATA);
  }
}

void CarpetIOHDF5_EvolutionCheckpoint(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT const iteration = cctk_iteration;
  CCTK_REAL const walltime = CCTK_RunTime() / 3600.0;

  bool const checkpoint_by_iteration =
      checkpoint_every > 0 and
      iteration >= last_checkpoint_iteration + checkpoint_every;
  bool const checkpoint_by_iteration_divisor =
      checkpoint_every_divisor > 0 and
      iteration % checkpoint_every_divisor == 0;
  bool const checkpoint_by_walltime =
      checkpoint_every_walltime_hours > 0 and
      walltime >= last_checkpoint_walltime + checkpoint_every_walltime_hours;

  int do_checkpoint =
      checkpoint and
      (checkpoint_by_iteration or checkpoint_by_iteration_divisor or
       checkpoint_by_walltime or checkpoint_next);
  if (checkpoint_every_walltime_hours > 0) {
    // broadcast the decision since comparing wall times may differ on
    // different processors
    MPI_Bcast(&do_checkpoint, 1, MPI_INT, 0, dist::comm());
  }

  if (do_checkpoint) {
    if (not CCTK_Equals(verbose, "none")) {
      CCTK_INFO("---------------------------------------------------------");
      CCTK_VInfo(CCTK_THORNSTRING, "Dumping periodic checkpoint at "
                                   "iteration %d, simulation time %g",
                 cctk_iteration, double(cctk_time));
      CCTK_INFO("---------------------------------------------------------");
    }

    Checkpoint(cctkGH, CP_EVOLUTION_DATA);
  }
}

void CarpetIOHDF5_TerminationCheckpoint(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (checkpoint and checkpoint_on_terminate) {

    if (last_checkpoint_iteration < cctk_iteration or
        (cctk_iteration == 0 and not checkpoint_ID)) {
      if (not CCTK_Equals(verbose, "none")) {
        CCTK_INFO("---------------------------------------------------------");
        CCTK_VInfo(CCTK_THORNSTRING, "Dumping termination checkpoint at "
                                     "iteration %d, simulation time %g",
                   cctk_iteration, double(cctk_time));
        CCTK_INFO("---------------------------------------------------------");
      }
      Checkpoint(cctkGH, CP_EVOLUTION_DATA);

    } else {
      if (not CCTK_Equals(verbose, "none")) {
        CCTK_INFO("---------------------------------------------------------");
        CCTK_VInfo(CCTK_THORNSTRING,
                   "Termination checkpoint already dumped "
                   "as last evolution checkpoint at iteration %d",
                   last_checkpoint_iteration);
        CCTK_INFO("---------------------------------------------------------");
      }
    }
  }
}

CCTK_INT CarpetIOHDF5_SetCheckpointGroups(CCTK_INT const *const groups,
                                          CCTK_INT const ngroups) {
  if (ngroups == -1) {
    // Checkpoint all groups
    groups_to_checkpoint.clear();
  } else {
    assert(ngroups >= 0);
    groups_to_checkpoint.resize(CCTK_NumGroups());
    for (int n = 0; n < CCTK_NumGroups(); ++n) {
      groups_to_checkpoint.at(n) = false;
    }
    for (int n = 0; n < ngroups; ++n) {
      groups_to_checkpoint.at(groups[n]) = true;
    }
  }
  return 0;
}

vector<bool> groups_to_checkpoint;

hid_t CCTKtoHDF5_Datatype(const cGH *const cctkGH, int cctk_type,
                          bool single_precision) {
  hid_t retval;

  const CarpetIOHDF5GH *myGH =
      (CarpetIOHDF5GH *)CCTK_GHExtension(cctkGH, CCTK_THORNSTRING);

  if (single_precision) {
    switch (cctk_type) {
    case CCTK_VARIABLE_REAL:
      cctk_type = CCTK_VARIABLE_REAL4;
      break;
    case CCTK_VARIABLE_COMPLEX:
      cctk_type = CCTK_VARIABLE_COMPLEX8;
      break;
    default:
      break; // do nothing
    }
  }

  // Normalize the type
  switch (cctk_type) {
  case CCTK_VARIABLE_INT:
#ifdef CCTK_INTEGER_PRECISION_1
    cctk_type = CCTK_VARIABLE_INT1;
#endif
#ifdef CCTK_INTEGER_PRECISION_2
    cctk_type = CCTK_VARIABLE_INT2;
#endif
#ifdef CCTK_INTEGER_PRECISION_4
    cctk_type = CCTK_VARIABLE_INT4;
#endif
#ifdef CCTK_INTEGER_PRECISION_8
    cctk_type = CCTK_VARIABLE_INT8;
#endif
#ifdef CCTK_INTEGER_PRECISION_16
    cctk_type = CCTK_VARIABLE_INT16;
#endif
    break;
  case CCTK_VARIABLE_REAL:
#ifdef CCTK_REAL_PRECISION_4
    cctk_type = CCTK_VARIABLE_REAL4;
#endif
#ifdef CCTK_REAL_PRECISION_8
    cctk_type = CCTK_VARIABLE_REAL8;
#endif
#ifdef CCTK_REAL_PRECISION_16
    cctk_type = CCTK_VARIABLE_REAL16;
#endif
    break;
  case CCTK_VARIABLE_COMPLEX:
#ifdef CCTK_REAL_PRECISION_4
    cctk_type = CCTK_VARIABLE_COMPLEX8;
#endif
#ifdef CCTK_REAL_PRECISION_8
    cctk_type = CCTK_VARIABLE_COMPLEX16;
#endif
#ifdef CCTK_REAL_PRECISION_16
    cctk_type = CCTK_VARIABLE_COMPLEX32;
#endif
    break;
  default:
    // do nothing
    break;
  }

  switch (cctk_type) {

  case CCTK_VARIABLE_CHAR:
    retval = HDF5_CHAR;
    break;

  case CCTK_VARIABLE_INT1:
    retval = H5T_NATIVE_INT8;
    break;
  case CCTK_VARIABLE_INT2:
    retval = H5T_NATIVE_INT16;
    break;
  case CCTK_VARIABLE_INT4:
    retval = H5T_NATIVE_INT32;
    break;
  case CCTK_VARIABLE_INT8:
    retval = H5T_NATIVE_INT64;
    break;
  // case CCTK_VARIABLE_INT16:     retval = H5T_NATIVE_INT128; break;

  case CCTK_VARIABLE_REAL4:
    retval = H5T_NATIVE_FLOAT;
    break;
  case CCTK_VARIABLE_REAL8:
    retval = H5T_NATIVE_DOUBLE;
    break;
  case CCTK_VARIABLE_REAL16:
    retval = H5T_NATIVE_LDOUBLE;
    break;

  case CCTK_VARIABLE_COMPLEX8:
    retval = myGH->HDF5_COMPLEX8;
    break;
  case CCTK_VARIABLE_COMPLEX16:
    retval = myGH->HDF5_COMPLEX16;
    break;
  case CCTK_VARIABLE_COMPLEX32:
    retval = myGH->HDF5_COMPLEX32;
    break;

  default:
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Unsupported CCTK variable datatype %d", cctk_type);
    retval = -1;
  }

  return (retval);
}

// add attributes to an HDF5 slice dataset
int AddSliceAttributes(const cGH *const cctkGH, const char *const fullname,
                       const int refinementlevel, const int multigridlevel,
                       const int map, const int timelevel,
                       const vector<double> &origin,
                       const vector<double> &delta, const vector<int> &iorigin,
                       const vector<int> &ioffset,
                       const vector<int> &ioffsetdenom, const vector<int> &bbox,
                       const vector<int> &nghostzones, const string &active,
                       hid_t &dataset, const vector<hsize_t> &shape,
                       const bool is_index) {
  int error_count = 0;

  error_count += WriteAttribute(dataset, "time", cctkGH->cctk_time);
  error_count += WriteAttribute(dataset, "timestep", cctkGH->cctk_iteration);
  error_count += WriteAttribute(dataset, "name", fullname);
  error_count += WriteAttribute(dataset, "level", refinementlevel);
  error_count += WriteAttribute(dataset, "carpet_mglevel", multigridlevel);
  error_count += WriteAttribute(dataset, "group_timelevel", timelevel);
  error_count += WriteAttribute(dataset, "origin", &origin[0], origin.size());
  error_count += WriteAttribute(dataset, "delta", &delta[0], delta.size());
  error_count +=
      WriteAttribute(dataset, "iorigin", &iorigin[0], iorigin.size());
  error_count +=
      WriteAttribute(dataset, "ioffset", &ioffset[0], ioffset.size());
  error_count += WriteAttribute(dataset, "ioffsetdenom", &ioffsetdenom[0],
                                ioffsetdenom.size());
  // TODO: Add "active" only for grid functions since it's trivial otherwise
  error_count += WriteAttribute(dataset, "active", active.c_str());
  // bbox and nghostzones are only used for grid functions and grid arrays
  if (bbox.size() > 0) {
    error_count += WriteAttribute(dataset, "cctk_bbox", &bbox[0], bbox.size());
  }
  if (nghostzones.size() > 0) {
    error_count += WriteAttribute(dataset, "cctk_nghostzones", &nghostzones[0],
                                  nghostzones.size());
  }
  // Specify whether the coordinate system is Cartesian or not
  if (CCTK_IsFunctionAliased("MultiPatch_MapIsCartesian")) {
    int const map_is_cartesian = MultiPatch_MapIsCartesian(map);
    error_count += WriteAttribute(dataset, "MapIsCartesian", map_is_cartesian);
  }
  if (is_index) {
    error_count += WriteAttribute(dataset, "h5shape", &shape[0], shape.size());
  }

  return error_count;
}

// Write an int attribute
int WriteAttribute(hid_t const group, char const *const name,
                   int const ivalue) {
  hid_t dataspace, attribute;
  int error_count = 0;

  HDF5_ERROR(dataspace = H5Screate(H5S_SCALAR));
  HDF5_ERROR(attribute = H5Acreate(group, name, H5T_NATIVE_INT, dataspace,
                                   H5P_DEFAULT));
  HDF5_ERROR(H5Awrite(attribute, H5T_NATIVE_INT, &ivalue));
  HDF5_ERROR(H5Aclose(attribute));
  HDF5_ERROR(H5Sclose(dataspace));

  return error_count;
}

// Write a double attribute
int WriteAttribute(hid_t const group, char const *const name,
                   double const dvalue) {
  hid_t dataspace, attribute;
  int error_count = 0;

  HDF5_ERROR(dataspace = H5Screate(H5S_SCALAR));
  HDF5_ERROR(attribute = H5Acreate(group, name, H5T_NATIVE_DOUBLE, dataspace,
                                   H5P_DEFAULT));
  HDF5_ERROR(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &dvalue));
  HDF5_ERROR(H5Aclose(attribute));
  HDF5_ERROR(H5Sclose(dataspace));

  return error_count;
}

// Write a string attribute
int WriteAttribute(hid_t const group, char const *const name,
                   char const *const svalue) {
  hid_t datatype, dataspace, attribute;
  int error_count = 0;

  HDF5_ERROR(datatype = H5Tcopy(H5T_C_S1));
  HDF5_ERROR(H5Tset_size(datatype, strlen(svalue) + 1));
  HDF5_ERROR(dataspace = H5Screate(H5S_SCALAR));
  HDF5_ERROR(attribute =
                 H5Acreate(group, name, datatype, dataspace, H5P_DEFAULT));
  HDF5_ERROR(H5Awrite(attribute, datatype, svalue));
  HDF5_ERROR(H5Aclose(attribute));
  HDF5_ERROR(H5Sclose(dataspace));
  HDF5_ERROR(H5Tclose(datatype));

  return error_count;
}

// Write an array of int attributes
int WriteAttribute(hid_t const group, char const *const name,
                   int const *const ivalues, int const nvalues) {
  hid_t dataspace, attribute;
  int error_count = 0;

  hsize_t const size = nvalues;
  HDF5_ERROR(dataspace = H5Screate_simple(1, &size, NULL));
  HDF5_ERROR(attribute = H5Acreate(group, name, H5T_NATIVE_INT, dataspace,
                                   H5P_DEFAULT));
  HDF5_ERROR(H5Awrite(attribute, H5T_NATIVE_INT, ivalues));
  HDF5_ERROR(H5Aclose(attribute));
  HDF5_ERROR(H5Sclose(dataspace));

  return error_count;
}

// Write an array of double attributes
int WriteAttribute(hid_t const group, char const *const name,
                   double const *const dvalues, int const nvalues) {
  hid_t dataspace, attribute;
  int error_count = 0;

  hsize_t const size = nvalues;
  HDF5_ERROR(dataspace = H5Screate_simple(1, &size, NULL));
  HDF5_ERROR(attribute = H5Acreate(group, name, H5T_NATIVE_DOUBLE, dataspace,
                                   H5P_DEFAULT));
  HDF5_ERROR(H5Awrite(attribute, H5T_NATIVE_DOUBLE, dvalues));
  HDF5_ERROR(H5Aclose(attribute));
  HDF5_ERROR(H5Sclose(dataspace));

  return error_count;
}

// Write an array of hsize_t attributes
int WriteAttribute(hid_t const group, char const *const name,
                   hsize_t const *const svalues, int const nvalues) {
  hid_t dataspace, attribute;
  int error_count = 0;

  hsize_t const size = nvalues;
  HDF5_ERROR(dataspace = H5Screate_simple(1, &size, NULL));
  HDF5_ERROR(attribute = H5Acreate(group, name, H5T_NATIVE_HSIZE, dataspace,
                                   H5P_DEFAULT));
  HDF5_ERROR(H5Awrite(attribute, H5T_NATIVE_HSIZE, svalues));
  HDF5_ERROR(H5Aclose(attribute));
  HDF5_ERROR(H5Sclose(dataspace));

  return error_count;
}

// Write an array of string attributes
int WriteAttribute(hid_t const group, char const *const name,
                   char const *const *const svalues, int const nvalues) {
  hid_t datatype, dataspace, attribute;
  int error_count = 0;

  size_t maxstrlen = 0;
  for (int i = 0; i < nvalues; ++i) {
    maxstrlen = max(maxstrlen, strlen(svalues[i]));
  }
  vector<char> svalue(nvalues * (maxstrlen + 1));
  for (int i = 0; i < nvalues; ++i) {
    strncpy(&svalue.at(i * maxstrlen), svalues[i], maxstrlen + 1);
  }

  HDF5_ERROR(datatype = H5Tcopy(H5T_C_S1));
  HDF5_ERROR(H5Tset_size(datatype, maxstrlen + 1));
  hsize_t const size = nvalues;
  HDF5_ERROR(dataspace = H5Screate_simple(1, &size, NULL));
  HDF5_ERROR(attribute =
                 H5Acreate(group, name, datatype, dataspace, H5P_DEFAULT));
  HDF5_ERROR(H5Awrite(attribute, datatype, &svalue.front()));
  HDF5_ERROR(H5Aclose(attribute));
  HDF5_ERROR(H5Sclose(dataspace));
  HDF5_ERROR(H5Tclose(datatype));

  return error_count;
}

// Write a large string attribute
int WriteLargeAttribute(hid_t const group, char const *const name,
                        char const *const svalue) {
  hid_t dataspace, dataset, datatype;
  int error_count = 0;

  // Create a dataset, since the data may not fit into an attribute
  hsize_t const size = strlen(svalue) + 1;
  HDF5_ERROR(datatype = H5Tcopy(H5T_C_S1));
  HDF5_ERROR(H5Tset_size(datatype, size));

  HDF5_ERROR(dataspace = H5Screate(H5S_SCALAR));
  HDF5_ERROR(dataset =
                 H5Dcreate(group, name, datatype, dataspace, H5P_DEFAULT));
  HDF5_ERROR(
      H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, svalue));
  HDF5_ERROR(H5Dclose(dataset));
  HDF5_ERROR(H5Sclose(dataspace));
  HDF5_ERROR(H5Tclose(datatype));

  return error_count;
}

//////////////////////////////////////////////////////////////////////////////
// private routines
//////////////////////////////////////////////////////////////////////////////
static void *SetupGH(tFleshConfig *const fleshconfig, const int convLevel,
                     cGH *const cctkGH) {
  DECLARE_CCTK_PARAMETERS;

  int error_count = 0;

  // register CarpetIOHDF5's routines as a new I/O method
  const int IOMethod = CCTK_RegisterIOMethod("IOHDF5");
  CCTK_RegisterIOMethodOutputGH(IOMethod, OutputGH);
  CCTK_RegisterIOMethodOutputVarAs(IOMethod, OutputVarAs);
  CCTK_RegisterIOMethodTimeToOutput(IOMethod, TimeToOutput);
  CCTK_RegisterIOMethodTriggerOutput(IOMethod, TriggerOutput);

  if (not CCTK_Equals(verbose, "none")) {
    CCTK_INFO("I/O Method 'IOHDF5' registered: AMR output of grid variables "
              "to HDF5 files");
  }

  // register CarpetIOHDF5's recovery routine with IOUtil
  if (IOUtil_RegisterRecover("CarpetIOHDF5 recovery", Recover) < 0) {
    CCTK_WARN(1, "Failed to register " CCTK_THORNSTRING " recovery routine");
  }

  const int numvars = CCTK_NumVars();

  // allocate a new GH extension structure
  CarpetIOHDF5GH *myGH = new CarpetIOHDF5GH;

  myGH->requests.resize(numvars);
  myGH->cp_filename_index = 0;
  myGH->checkpoint_keep = abs(checkpoint_keep);
  myGH->cp_filename_list =
      (char **)calloc(myGH->checkpoint_keep, sizeof(char *));
  myGH->recovery_num_filenames = 0;
  myGH->recovery_filename_list = NULL;
  myGH->out_vars = strdup("");
  myGH->out_every_default = out_every - 1;

  // initial I/O parameter check
  myGH->out_dir = 0;
  myGH->stop_on_parse_errors = strict_io_parameter_check;
  CheckSteerableParameters(cctkGH, myGH);
  myGH->stop_on_parse_errors = 0;

  // check parallel I/O parameters for chunked output
  if (not(CCTK_EQUALS(out_mode, "onefile") or CCTK_EQUALS(out_mode, "proc"))) {
    CCTK_VWarn(CCTK_WARN_COMPLAIN, __LINE__, __FILE__, CCTK_THORNSTRING,
               "IO::out_mode = '%s' is not implemented in %s. "
               "Defaulting to one output file per processor...",
               out_mode, CCTK_THORNSTRING);
  }
  if (CCTK_EQUALS(out_mode, "proc") and io_out_unchunked) {
    CCTK_WARN(CCTK_WARN_COMPLAIN,
              "IO::out_unchunked = 'yes' is incompatible with IO::out_mode = "
              "'proc'. Ignoring setting for IO::out_unchunked...");
  }

  if (output_index && CCTK_EQUALS(out_mode, "onefile"))
    CCTK_VWarn(CCTK_WARN_COMPLAIN, __LINE__, __FILE__, CCTK_THORNSTRING,
               "CarpetIOHDF5::output_index with IO::out_mode = '%s' is not "
               "implemented in %s. "
               "Index will not be output",
               out_mode, CCTK_THORNSTRING);

  // Now set hdf5 complex datatypes
  HDF5_ERROR(myGH->HDF5_COMPLEX =
                 H5Tcreate(H5T_COMPOUND, sizeof(CCTK_COMPLEX)));
  HDF5_ERROR(H5Tinsert(myGH->HDF5_COMPLEX, "real", 0, HDF5_REAL));
  HDF5_ERROR(
      H5Tinsert(myGH->HDF5_COMPLEX, "imag", sizeof(CCTK_REAL), HDF5_REAL));
#ifdef HAVE_CCTK_REAL4
  HDF5_ERROR(myGH->HDF5_COMPLEX8 =
                 H5Tcreate(H5T_COMPOUND, sizeof(CCTK_COMPLEX8)));
  HDF5_ERROR(H5Tinsert(myGH->HDF5_COMPLEX8, "real", 0, H5T_NATIVE_FLOAT));
  HDF5_ERROR(H5Tinsert(myGH->HDF5_COMPLEX8, "imag", sizeof(CCTK_REAL4),
                       H5T_NATIVE_FLOAT));
#endif
#ifdef HAVE_CCTK_REAL8
  HDF5_ERROR(myGH->HDF5_COMPLEX16 =
                 H5Tcreate(H5T_COMPOUND, sizeof(CCTK_COMPLEX16)));
  HDF5_ERROR(H5Tinsert(myGH->HDF5_COMPLEX16, "real", 0, H5T_NATIVE_DOUBLE));
  HDF5_ERROR(H5Tinsert(myGH->HDF5_COMPLEX16, "imag", sizeof(CCTK_REAL8),
                       H5T_NATIVE_DOUBLE));
#endif
#ifdef HAVE_CCTK_REAL16
  HDF5_ERROR(myGH->HDF5_COMPLEX32 =
                 H5Tcreate(H5T_COMPOUND, sizeof(CCTK_COMPLEX32)));
  HDF5_ERROR(H5Tinsert(myGH->HDF5_COMPLEX32, "real", 0, H5T_NATIVE_LDOUBLE));
  HDF5_ERROR(H5Tinsert(myGH->HDF5_COMPLEX32, "imag", sizeof(CCTK_REAL16),
                       H5T_NATIVE_LDOUBLE));
#endif

  return (myGH);
}

static void CheckSteerableParameters(const cGH *const cctkGH,
                                     CarpetIOHDF5GH *myGH) {
  DECLARE_CCTK_PARAMETERS;

  // re-parse the 'IOHDF5::out_dir' parameter if it has changed
  const char *my_out_dir = *out_dir ? out_dir : io_out_dir;
  char *the_out_dir;
  if (strcmp(my_out_dir, ".")) {
    int i = strlen(my_out_dir);
    if (CCTK_Equals(out_mode, "onefile") or not strstr(my_out_dir, "%u")) {
      the_out_dir = (char *)malloc(i + 2);
      strcpy(the_out_dir, my_out_dir);
      the_out_dir[i] = '/';
      the_out_dir[i + 1] = 0;
    } else {
      // TODO: ensure that there is exactly one "%u" and no other "%"
      // substrings, except possibly "%%".
      the_out_dir = (char *)malloc(i + 20);
      snprintf(the_out_dir, i + 19, my_out_dir, dist::rank());
      strcat(the_out_dir, "/");
    }
  } else {
    the_out_dir = strdup("");
  }

  if (not myGH->out_dir or strcmp(the_out_dir, myGH->out_dir)) {
    free(myGH->out_dir);
    myGH->out_dir = the_out_dir;

    // create the output directory
    // const ioGH* const ioUtilGH = (const ioGH*) CCTK_GHExtension (cctkGH,
    // "IO");
    int result = IOUtil_CreateDirectory(cctkGH, myGH->out_dir,
                                        not CCTK_Equals(out_mode, "onefile"),
                                        dist::rank());
    if (result < 0) {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Problem creating HDF5 output directory '%s'", myGH->out_dir);
    } else if (result > 0 and CCTK_Equals(verbose, "full")) {
      CCTK_VInfo(CCTK_THORNSTRING, "HDF5 output directory '%s' already exists",
                 myGH->out_dir);
    }
  } else {
    free(the_out_dir);
  }

  // re-parse the 'IOHDF5::out_vars' parameter if it has changed
  if (strcmp(out_vars, myGH->out_vars)) {
#ifdef IOUTIL_PARSER_HAS_OUT_DT
    IOUtil_ParseVarsForOutput(cctkGH, CCTK_THORNSTRING, "IOHDF5::out_vars",
                              myGH->stop_on_parse_errors, out_vars, -1, -1.0,
                              &myGH->requests[0]);
#else
    IOUtil_ParseVarsForOutput(cctkGH, CCTK_THORNSTRING, "IOHDF5::out_vars",
                              myGH->stop_on_parse_errors, out_vars, -1,
                              &myGH->requests[0]);
#endif

    // notify the user about the new setting
    if (not CCTK_Equals(verbose, "none")) {
      int count = 0;
      ostringstream msg;
      msg << "Periodic AMR output requested for:";
      for (int vi = 0; vi < CCTK_NumVars(); ++vi) {
        if (myGH->requests[vi]) {
          ++count;
          char *const fullname = CCTK_FullName(vi);
          msg << eol << "   " << fullname;
          free(fullname);
        }
      }
      if (count > 0) {
        CCTK_INFO(msg.str().c_str());
      }
    }

    // save the last setting of 'IOHDF5::out_vars' parameter
    free(myGH->out_vars);
    myGH->out_vars = strdup(out_vars);
  }
}

static int OutputGH(const cGH *const cctkGH) {
  static Timers::Timer timer("OutputGH");
  timer.start();

  CarpetIOHDF5GH *myGH =
      (CarpetIOHDF5GH *)CCTK_GHExtension(cctkGH, CCTK_THORNSTRING);
  CheckSteerableParameters(cctkGH, myGH);

  if (strcmp(myGH->out_vars, "")) {
    for (int vindex = CCTK_NumVars() - 1; vindex >= 0; vindex--) {
      if (TimeToOutput(cctkGH, vindex)) {
        TriggerOutput(cctkGH, vindex);
      }
    }
  }
  timer.stop();

  return (0);
}

static int TimeToOutput(const cGH *const cctkGH, const int vindex) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const int numvars = CCTK_NumVars();
  assert(vindex >= 0 and vindex < numvars);

  if (CCTK_GroupTypeFromVarI(vindex) != CCTK_GF and not do_global_mode) {
    return 0;
  }

  CarpetIOHDF5GH *myGH =
      (CarpetIOHDF5GH *)CCTK_GHExtension(cctkGH, CCTK_THORNSTRING);

  // check if output for this variable was requested
  if (not myGH->requests[vindex]) {
    return (0);
  }

  // check whether this refinement level should be output
  if (not(myGH->requests[vindex]->refinement_levels & (1 << reflevel))) {
    return (0);
  }

  // check if output for this variable was requested individually
  // by a "<varname>{ out_every = <number> }" option string
  // this will overwrite the output criterion setting
  const char *myoutcriterion =
      CCTK_EQUALS(out_criterion, "default") ? io_out_criterion : out_criterion;
  if (myGH->requests[vindex]->out_every >= 0) {
    myoutcriterion = "divisor";
  }

  if (CCTK_EQUALS(myoutcriterion, "never")) {
    return (0);
  }

  // check whether to output at this iteration
  bool output_this_iteration = false;

  if (CCTK_EQUALS(myoutcriterion, "iteration")) {
    int myoutevery = out_every == -2 ? io_out_every : out_every;
    if (myoutevery > 0) {
      if (*this_iteration == cctk_iteration) {
        // we already decided to output this iteration
        output_this_iteration = true;
      } else if (cctk_iteration >= *next_output_iteration) {
        // it is time for the next output
        output_this_iteration = true;
        *this_iteration = cctk_iteration;
        *next_output_iteration = cctk_iteration + myoutevery;
      }
    }
  } else if (CCTK_EQUALS(myoutcriterion, "divisor")) {
    int myoutevery = out_every == -2 ? io_out_every : out_every;
    if (myGH->requests[vindex]->out_every >= 0) {
      myoutevery = myGH->requests[vindex]->out_every;
    }
    if (myoutevery > 0 and (cctk_iteration % myoutevery) == 0) {
      // we already decided to output this iteration
      output_this_iteration = true;
    }
  } else if (CCTK_EQUALS(myoutcriterion, "time")) {
    CCTK_REAL myoutdt = out_dt == -2 ? io_out_dt : out_dt;
    if (myoutdt == 0 or *this_iteration == cctk_iteration) {
      output_this_iteration = true;
    } else if (myoutdt > 0) {
      int do_output = (cctk_time / cctk_delta_time >=
                       *next_output_time / cctk_delta_time - 1.0e-12);
      MPI_Bcast(&do_output, 1, MPI_INT, 0, dist::comm());
      if (do_output) {
        // it is time for the next output
        output_this_iteration = true;
        *this_iteration = cctk_iteration;
        *next_output_time = cctk_time + myoutdt;
      }
    }
  }

  return output_this_iteration ? 1 : 0;
}

static int TriggerOutput(const cGH *const cctkGH, const int vindex) {
  DECLARE_CCTK_PARAMETERS;
  int retval;

  char *const fullname = CCTK_FullName(vindex);
  if (one_file_per_group) {
    const int gindex = CCTK_GroupIndexFromVarI(vindex);
    char *const groupname_c = CCTK_GroupName(gindex);
    string groupname(groupname_c);
    free(groupname_c);
    transform(groupname.begin(), groupname.end(), groupname.begin(), ::tolower);
    string const oldsep("::");
    size_t const oldseppos = groupname.find(oldsep);
    assert(oldseppos != string::npos);
    groupname.replace(oldseppos, oldsep.size(), out_group_separator);
    retval = OutputVarAs(cctkGH, fullname, groupname.c_str());
  } else {
    const char *varname = CCTK_VarName(vindex);
    retval = OutputVarAs(cctkGH, fullname, varname);
  }
  free(fullname);

  return (retval);
}

static void GetVarIndex(int vindex, const char *optstring, void *arg) {
  if (optstring) {
    char *fullname = CCTK_FullName(vindex);
    CCTK_VWarn(2, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Option string '%s' will be ignored for HDF5 output of "
               "variable '%s'",
               optstring, fullname);
    free(fullname);
  }

  *((int *)arg) = vindex;
}

static int OutputVarAs(const cGH *const cctkGH, const char *const fullname,
                       const char *const alias) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int error_count = 0;
  int vindex = -1;

  if (CCTK_TraverseString(fullname, GetVarIndex, &vindex, CCTK_VAR) < 0) {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "error while parsing variable name '%s' (alias name '%s')",
               fullname, alias);
    return (-1);
  }

  if (vindex < 0) {
    return (-1);
  }

  const int group = CCTK_GroupIndexFromVarI(vindex);
  assert(group >= 0);
  cGroup groupdata;
  CCTK_GroupData(group, &groupdata);
  if (groupdata.grouptype == CCTK_SCALAR or groupdata.grouptype == CCTK_ARRAY) {
    assert(do_global_mode);
  }

  // get the default I/O request for this variable
  const CarpetIOHDF5GH *myGH =
      (CarpetIOHDF5GH *)CCTK_GHExtension(cctkGH, CCTK_THORNSTRING);
  ioRequest *request = myGH->requests[vindex];
  if (not request) {
#ifdef IOUTIL_PARSER_HAS_OUT_DT
    request = IOUtil_DefaultIORequest(cctkGH, vindex, 1, -1.0);
#else
    request = IOUtil_DefaultIORequest(cctkGH, vindex, 1);
#endif
  }

  // Get grid hierarchy extentsion from IOUtil
  const ioGH *const iogh = (const ioGH *)CCTK_GHExtension(cctkGH, "IO");
  assert(iogh);

  // Invent a file name
  int ioproc = 0, nioprocs = 1;
  string filename;
  filename.append(myGH->out_dir);
  filename.append(alias);
  if (out_timesteps_per_file > 0) {
    // Round down to nearest multiple of out_timesteps_per_file
    int const iter =
        cctk_iteration / out_timesteps_per_file * out_timesteps_per_file;
    char buffer[32];
    snprintf(buffer, sizeof(buffer), ".iter_%d", iter);
    filename.append(buffer);
  }
  if (not(CCTK_EQUALS(out_mode, "onefile") or request->out_unchunked or
          groupdata.disttype == CCTK_DISTRIB_CONSTANT or dist::size() == 1)) {
    char buffer[32];
    ioproc = dist::rank();
    nioprocs = dist::size();
    snprintf(buffer, sizeof(buffer), ".file_%d", ioproc);
    filename.append(buffer);
  }

  string base_filename(filename);
  filename.append(out_extension);

  string index_filename(base_filename + ".idx" + out_extension);

  const char *const c_filename = filename.c_str();

  // check if the file has been created already
  typedef std::map<string, vector<vector<vector<int> > > > filelist;
  static filelist created_files;
  filelist::iterator thisfile = created_files.find(filename);
  bool is_new_file = thisfile == created_files.end();
  if (is_new_file) {
    int const numvars = CCTK_NumVars();
    vector<vector<vector<int> > > last_outputs; // [ml][rl][var]
    last_outputs.resize(mglevels);
    for (int ml = 0; ml < mglevels; ++ml) {
      last_outputs[ml].resize(maxreflevels);
      for (int rl = 0; rl < maxreflevels; ++rl) {
        last_outputs[ml][rl].resize(numvars, cctk_iteration - 1);
      }
    }
    thisfile = created_files.insert(
        thisfile, filelist::value_type(filename, last_outputs));
    assert(thisfile != created_files.end());
  }

  const int firstvar = one_file_per_group ? CCTK_FirstVarIndexI(group) : vindex;
  const int numvars = one_file_per_group ? CCTK_NumVarsInGroupI(group) : 1;

  // check if this variable has been output already during this iteration
  int &last_output = thisfile->second.at(mglevel).at(reflevel).at(vindex);
  if (last_output == cctk_iteration) {
    // Has already been output during this iteration
    if (not one_file_per_group or vindex == firstvar + numvars - 1) {
      char *varname = CCTK_FullName(vindex);
      CCTK_VWarn(5, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Skipping output for variable \"%s\", because this variable "
                 "has already been output during the current iteration -- "
                 "probably via a trigger during the analysis stage",
                 varname);
      free(varname);
    }
    return (0);
  }
  assert(last_output < cctk_iteration);
  last_output = cctk_iteration;

  // Check for storage
  if (not CCTK_QueryGroupStorageI(cctkGH, group)) {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Cannot output variable '%s' because it has no storage",
               fullname);
    return (0);
  }

  // Open the output file if this is a designated I/O processor
  hid_t file = -1;
  hid_t index_file = -1;
  CCTK_REAL io_files = 0;
  CCTK_REAL io_bytes = 0;
  BeginTimingIO(cctkGH);
  if (dist::rank() == ioproc) {

    if (is_new_file and not IO_TruncateOutputFiles(cctkGH)) {
      H5E_BEGIN_TRY { is_new_file = H5Fis_hdf5(c_filename) <= 0; }
      H5E_END_TRY;
    }

    if (is_new_file) {
      hid_t fapl_id;
      HDF5_ERROR(fapl_id = H5Pcreate(H5P_FILE_ACCESS));
      HDF5_ERROR(H5Pset_fclose_degree(fapl_id, H5F_CLOSE_STRONG));
      HDF5_ERROR(
          file = H5Fcreate(c_filename, H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id));
      if (output_index) {
        HDF5_ERROR(index_file = H5Fcreate(index_filename.c_str(), H5F_ACC_TRUNC,
                                          H5P_DEFAULT, fapl_id));
      }
      HDF5_ERROR(H5Pclose(fapl_id));
      // write metadata information
      error_count +=
          WriteMetadata(cctkGH, nioprocs, firstvar, numvars, false, file);

      if (output_index) {
        error_count += WriteMetadata(cctkGH, nioprocs, firstvar, numvars, false,
                                     index_file);
      }
    } else {
      hid_t fapl_id;
      HDF5_ERROR(fapl_id = H5Pcreate(H5P_FILE_ACCESS));
      HDF5_ERROR(H5Pset_fclose_degree(fapl_id, H5F_CLOSE_STRONG));
      HDF5_ERROR(file = H5Fopen(c_filename, H5F_ACC_RDWR, fapl_id));
      if (output_index)
        HDF5_ERROR(index_file =
                       H5Fopen(index_filename.c_str(), H5F_ACC_RDWR, fapl_id));
      HDF5_ERROR(H5Pclose(fapl_id));
    }
    io_files += 1;
  }

  if (CCTK_Equals(verbose, "full")) {
    CCTK_VInfo(CCTK_THORNSTRING,
               "Writing variable '%s' on mglevel %d reflevel %d", fullname,
               mglevel, reflevel);
  }
  for (int var = firstvar; var < firstvar + numvars; var++) {
    ioRequest *r = myGH->requests[var];
    if (not r) {
#ifdef IOUTIL_PARSER_HAS_OUT_DT
      r = IOUtil_DefaultIORequest(cctkGH, var, 1, -1.0);
#else
      r = IOUtil_DefaultIORequest(cctkGH, var, 1);
#endif
    }
    if ((CCTK_EQUALS(out_mode, "onefile") and io_out_unchunked) or
        r->out_unchunked or groupdata.disttype == CCTK_DISTRIB_CONSTANT) {
      error_count += WriteVarUnchunked(cctkGH, file, io_bytes, r, false);
    } else if (CCTK_EQUALS(out_mode, "onefile")) {
      error_count += WriteVarChunkedSequential(cctkGH, file, io_bytes, r, false,
                                               index_file);
    } else {
      error_count +=
          WriteVarChunkedParallel(cctkGH, file, io_bytes, r, false, index_file);
    }
    if (r != myGH->requests[var])
      IOUtil_FreeIORequest(&r);

    // mark this variable to have been output at this iteration
    thisfile->second.at(mglevel).at(reflevel).at(var) = cctk_iteration;
  }

  // free I/O request structure
  if (request != myGH->requests[vindex]) {
    IOUtil_FreeIORequest(&request);
  }

  // Close the file
  if (file >= 0) {
    HDF5_ERROR(H5Fclose(file));
    if (output_index)
      HDF5_ERROR(H5Fclose(index_file));
  }
  HDF5_ERROR(H5garbage_collect());
  {
    CCTK_REAL local[2], global[2];
    local[0] = io_files;
    local[1] = io_bytes;
    MPI_Allreduce(local, global, 2, dist::mpi_datatype(local[0]), MPI_SUM,
                  dist::comm());
    io_files = global[0];
    io_bytes = global[1];
  }
  EndTimingIO(cctkGH, io_files, io_bytes, true);

  if (error_count > 0 and abort_on_io_errors) {
    CCTK_WARN(0, "Aborting simulation due to previous I/O errors");
  }

  return (0);
}

static void Checkpoint(const cGH *const cctkGH, int called_from) {
  int error_count = 0;
  DECLARE_CCTK_PARAMETERS;

  /* get the filenames for both the temporary and real checkpoint file */
  int ioproc = 0, nioprocs = 1;
  int parallel_io = 0;
  if (not(CCTK_EQUALS(out_mode, "onefile") or dist::size() == 1)) {
    ioproc = dist::rank();
    nioprocs = dist::size();
    parallel_io = 1;
  }
  char *filename = IOUtil_AssembleFilename(cctkGH, NULL, "", ".h5", called_from,
                                           ioproc, not parallel_io);
  char *tempname = IOUtil_AssembleFilename(
      cctkGH, NULL, ".tmp", ".h5", called_from, ioproc, not parallel_io);

  char *index_tempname = IOUtil_AssembleFilename(
      cctkGH, NULL, ".tmp", ".idx.h5", called_from, ioproc, not parallel_io);
  char *index_filename = IOUtil_AssembleFilename(
      cctkGH, NULL, "", ".idx.h5", called_from, ioproc, not parallel_io);

  hid_t file = -1;
  hid_t index_file = -1;
  if (dist::rank() == ioproc) {
    if (CCTK_Equals(verbose, "full")) {
      CCTK_VInfo(CCTK_THORNSTRING, "Creating temporary checkpoint file '%s'",
                 tempname);
    }

    hid_t fapl_id;
    HDF5_ERROR(fapl_id = H5Pcreate(H5P_FILE_ACCESS));
    HDF5_ERROR(H5Pset_fclose_degree(fapl_id, H5F_CLOSE_STRONG));
    HDF5_ERROR(file = H5Fcreate(tempname, H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id));

    // write metadata information
    error_count += WriteMetadata(cctkGH, nioprocs, -1, -1, true, file);

    if (output_index) {

      HDF5_ERROR(index_file = H5Fcreate(index_tempname, H5F_ACC_TRUNC,
                                        H5P_DEFAULT, fapl_id));
      error_count += WriteMetadata(cctkGH, nioprocs, -1, -1, true, index_file);
    }
    HDF5_ERROR(H5Pclose(fapl_id));
  }

  // remember the current wall time
  CCTK_REAL const walltime = CCTK_RunTime() / 3600.0;

  // now dump the grid variables on all mglevels, reflevels, maps and components
  BEGIN_MGLEVEL_LOOP(cctkGH) {

    CCTK_REAL io_files = 1;
    CCTK_REAL io_bytes = 0;
    BeginTimingIO(cctkGH);

    BEGIN_REFLEVEL_LOOP(cctkGH) {

      if (CCTK_Equals(verbose, "full")) {
        CCTK_VInfo(CCTK_THORNSTRING,
                   "Dumping grid variables on mglevel %d reflevel %d ...",
                   mglevel, reflevel);
      }

      for (int group = CCTK_NumGroups() - 1; group >= 0; group--) {
        /* skip variables which have been disabled for checkpointing */
        if (not groups_to_checkpoint.empty() and
            not groups_to_checkpoint.at(group)) {
          continue;
        }
        /* only dump groups which have storage assigned */
        if (CCTK_QueryGroupStorageI(cctkGH, group) <= 0 or
            CCTK_NumVarsInGroupI(group) == 0) {
          continue;
        }

        cGroup gdata;
        CCTK_GroupData(group, &gdata);
        assert(gdata.grouptype == CCTK_ARRAY or gdata.grouptype == CCTK_GF or
               gdata.grouptype == CCTK_SCALAR);

        // scalars and grid arrays only have one reflevel
        if (gdata.grouptype != CCTK_GF and reflevel > 0) {
          continue;
        }

        int const len =
            Util_TableGetString(gdata.tagstable, 0, NULL, "checkpoint");
        if (len > 0) {
          vector<char> value_buf(len + 1);
          char *value = &value_buf[0];
          Util_TableGetString(gdata.tagstable, len + 1, value, "checkpoint");
          if (len == sizeof("no") - 1 and CCTK_Equals(value, "no")) {
            continue;
          } else if (not CCTK_Equals(value, "yes")) {
            char *groupname = CCTK_GroupName(group);
            CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                       "Ignoring unknown checkpoint tag '%s' for group '%s'",
                       value, groupname);
            free(groupname);
          }
        }

        /* get the number of active timelevels */
        gdata.numtimelevels = CCTK_ActiveTimeLevelsGI(cctkGH, group);

        int first_vindex = CCTK_FirstVarIndexI(group);

/* get the default I/O request for this group */
#ifdef IOUTIL_PARSER_HAS_OUT_DT
        ioRequest *request =
            IOUtil_DefaultIORequest(cctkGH, first_vindex, 1, -1.0);
#else
        ioRequest *request = IOUtil_DefaultIORequest(cctkGH, first_vindex, 1);
#endif

        /* disable checking for old data objects, disable datatype conversion
           and downsampling */
        request->check_exist = 0;
        request->hdatatype = gdata.vartype;
        for (request->hdim = 0; request->hdim < request->vdim;
             request->hdim++) {
          request->downsample[request->hdim] = 1;
        }

        /* loop over all variables in this group */
        for (request->vindex = first_vindex;
             request->vindex < first_vindex + gdata.numvars;
             request->vindex++) {
          char *fullname = CCTK_FullName(request->vindex);
          assert(fullname);

          /* loop over all timelevels of this variable */
          for (request->timelevel = 0; request->timelevel < gdata.numtimelevels;
               request->timelevel++) {
            if (CCTK_Equals(verbose, "full")) {
              CCTK_VInfo(CCTK_THORNSTRING, "  %s (timelevel %d)", fullname,
                         request->timelevel);
            }

            // write the var
            error_count +=
                parallel_io
                    ? WriteVarChunkedParallel(cctkGH, file, io_bytes, request,
                                              true, index_file)
                    : WriteVarChunkedSequential(cctkGH, file, io_bytes, request,
                                                true, index_file);
          }
          free(fullname);

        } /* end of loop over all variables */

        // free I/O request structure
        IOUtil_FreeIORequest(&request);

      } /* end of loop over all groups */
    }
    END_REFLEVEL_LOOP;

    {
      CCTK_REAL local[2], global[2];
      local[0] = io_files;
      local[1] = io_bytes;
      MPI_Allreduce(local, global, 2, dist::mpi_datatype(local[0]), MPI_SUM,
                    dist::comm());
      io_files = global[0];
      io_bytes = global[1];
    }
    EndTimingIO(cctkGH, io_files, io_bytes, true);
  }
  END_MGLEVEL_LOOP;

  // Close the file
  if (file >= 0) {
    HDF5_ERROR(H5Fclose(file));
  }
  if (index_file >= 0) {
    HDF5_ERROR(H5Fclose(index_file));
  }
  HDF5_ERROR(H5garbage_collect());

  // get global error count
  int temp = error_count;
  MPI_Allreduce(&temp, &error_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if (error_count == 0) {
    if (file >= 0) {
      if (rename(tempname, filename)) {
        CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Could not rename temporary checkpoint file '%s' to '%s'",
                   tempname, filename);
        error_count = -1;
      } else if (called_from == CP_EVOLUTION_DATA and checkpoint_keep > 0) {
        CarpetIOHDF5GH *myGH =
            (CarpetIOHDF5GH *)CCTK_GHExtension(cctkGH, CCTK_THORNSTRING);

        // should an older checkpoint file be removed ?
        if (myGH->cp_filename_list[myGH->cp_filename_index]) {
          // check whether the recovery checkpoint (which can be a list of
          // several chunked files) or a checkpoint file should be removed
          if (myGH->recovery_filename_list) {
            for (int i = 0; i < myGH->recovery_num_filenames; i++) {
              if (myGH->recovery_filename_list[i]) {
                // remove possible index file
                string old_index_filename = myGH->recovery_filename_list[i];
                size_t const basenamelen =
                    old_index_filename.rfind(out_extension);
                if (basenamelen != string::npos) {
                  old_index_filename.insert(basenamelen, ".idx");
                  remove(old_index_filename.c_str());
                }

                remove(myGH->recovery_filename_list[i]);
                free(myGH->recovery_filename_list[i]);
              }
            }
            free(myGH->recovery_filename_list);
            myGH->recovery_filename_list = NULL;
          } else {
            // remove possible index file
            string old_index_filename =
                myGH->cp_filename_list[myGH->cp_filename_index];
            size_t const basenamelen = old_index_filename.rfind(out_extension);
            if (basenamelen != string::npos) {
              old_index_filename.insert(basenamelen, ".idx");
              remove(old_index_filename.c_str());
            }
            remove(myGH->cp_filename_list[myGH->cp_filename_index]);
            free(myGH->cp_filename_list[myGH->cp_filename_index]);
          }
        }

        // add this checkpoint to the checkpoint filename ring buffer
        myGH->cp_filename_list[myGH->cp_filename_index] = strdup(filename);
        myGH->cp_filename_index =
            (myGH->cp_filename_index + 1) % checkpoint_keep;

        // since the 'checkpoint_keep' parameter is steerable,
        // we may need to resize the ring buffer
        if (myGH->checkpoint_keep != checkpoint_keep) {
          char **cp_filename_list =
              (char **)calloc(checkpoint_keep, sizeof(char *));
          int min = myGH->checkpoint_keep < checkpoint_keep
                        ? myGH->checkpoint_keep
                        : checkpoint_keep;
          while (min-- > 0) {
            cp_filename_list[min] = myGH->cp_filename_list[min];
          }
          free(myGH->cp_filename_list);
          myGH->cp_filename_list = cp_filename_list;
          myGH->checkpoint_keep = checkpoint_keep;
        }
      }
    }
    if (index_file >= 0) {
      if (rename(index_tempname, index_filename)) {
        CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Could not rename temporary index file '%s' to '%s'",
                   index_tempname, index_filename);
      }
    }
  } else {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Failed to create checkpoint at iteration %d",
               cctkGH->cctk_iteration);
  }

  // save the iteration number of this checkpoint
  last_checkpoint_iteration = cctkGH->cctk_iteration;
  last_checkpoint_walltime = walltime;
  if (checkpoint_next) {
    CCTK_ParameterSet("checkpoint_next", CCTK_THORNSTRING, "no");
  }

  // free allocated resources
  free(index_tempname);
  free(index_filename);
  free(tempname);
  free(filename);

  if (error_count > 0 and abort_on_io_errors) {
    CCTK_WARN(0, "Aborting simulation due to previous I/O errors");
  }

} // Checkpoint

int WriteMetadata(const cGH *const cctkGH, int const nioprocs,
                  int const firstvar, int const numvars,
                  bool const called_from_checkpoint, hid_t const file) {
  DECLARE_CCTK_PARAMETERS;
  hid_t group;
  int error_count = 0;

  if (CCTK_Equals(verbose, "full")) {
    CCTK_INFO("Writing simulation metadata...");
  }

  // const ioGH *ioUtilGH = (const ioGH *) CCTK_GHExtension (cctkGH, "IO");
  HDF5_ERROR(group = H5Gcreate(file, METADATA_GROUP, 0));

  error_count += WriteAttribute(group, "main loop index", CCTK_MainLoopIndex());
  error_count += WriteAttribute(group, "GH$iteration", cctkGH->cctk_iteration);
  error_count += WriteAttribute(group, "nioprocs", nioprocs);
  error_count += WriteAttribute(group, "carpet_reflevels", reflevels);
  error_count += WriteAttribute(group, "carpet_global_time", global_time);
  error_count += WriteAttribute(group, "carpet_delta_time", delta_time);
  error_count += WriteAttribute(group, "Cactus version", CCTK_FullVersion());

  // all times on all refinement levels
  // error_count += WriteAttribute (group,
  //                                "numberofmgtimes", mglevels);
  // for (int i = 0; i < mglevels; i++) {
  //   char name[100];
  //   snprintf (name, sizeof (name), "mgleveltimes %d", i);
  //   error_count += WriteAttribute
  //     (group, name, &leveltimes.at(i).front(), leveltimes.at(i).size());
  // }

  // unique configuration identifier
  if (CCTK_IsFunctionAliased("UniqueConfigID")) {
    error_count += WriteAttribute(
        group, "config id", static_cast<char const *>(UniqueConfigID(cctkGH)));
  }

  // unique build identifier
  if (CCTK_IsFunctionAliased("UniqueBuildID")) {
    error_count += WriteAttribute(
        group, "build id", static_cast<char const *>(UniqueBuildID(cctkGH)));
  }

  // unique simulation identifier
  if (CCTK_IsFunctionAliased("UniqueSimulationID")) {
    error_count +=
        WriteAttribute(group, "simulation id",
                       static_cast<char const *>(UniqueSimulationID(cctkGH)));
  }

  // unique run identifier
  if (CCTK_IsFunctionAliased("UniqueRunID")) {
    error_count += WriteAttribute(
        group, "run id", static_cast<char const *>(UniqueRunID(cctkGH)));
  }

  // list all datasets in this file
  {
    ostringstream buf;
    if (firstvar >= 0) {
      // list the selected variables
      assert(firstvar >= 0 and numvars >= 0);
      for (int vi = 0; vi < numvars; ++vi) {
        char *const fullname = CCTK_FullName(firstvar + vi);
        buf << fullname << "\n";
        free(fullname);
      }
    } else {
      // list all variables with storage and without checkpoint="no"
      // TODO: check this only once, then pass the list of variables around
      ENTER_GLOBAL_MODE(cctkGH, 0) {
        for (int gi = 0; gi < CCTK_NumGroups(); ++gi) {
          int const numvars_group = CCTK_NumVarsInGroupI(gi);
          if (numvars_group > 0) {
            if (CCTK_QueryGroupStorageI(cctkGH, gi) > 0) {
              cGroup gdata;
              CCTK_GroupData(gi, &gdata);
              bool do_checkpoint = true;
              int const len =
                  Util_TableGetString(gdata.tagstable, 0, NULL, "checkpoint");
              if (len >= 0) {
                char value[len + 1];
                Util_TableGetString(gdata.tagstable, len + 1, value,
                                    "checkpoint");
                if (CCTK_Equals(value, "no")) {
                  do_checkpoint = false;
                } else if (CCTK_Equals(value, "yes")) {
                  do_checkpoint = true;
                } else {
                  char *const groupname = CCTK_GroupName(gi);
                  CCTK_VWarn(
                      1, __LINE__, __FILE__, CCTK_THORNSTRING,
                      "Ignoring unknown checkpoint tag '%s' for group '%s'",
                      value, groupname);
                  free(groupname);
                }
              }
              if (do_checkpoint) {
                int const firstvar_group = CCTK_FirstVarIndexI(gi);
                for (int vi = 0; vi < numvars_group; ++vi) {
                  char *const fullname = CCTK_FullName(firstvar_group + vi);
                  buf << fullname << "\n";
                  free(fullname);
                }
              }
            }
          }
        } // for vi
      }
      LEAVE_GLOBAL_MODE;
    } // if firstvar < 0
    string const str = buf.str();
    error_count += WriteLargeAttribute(group, "Datasets", str.c_str());
  }

  // save parameters in a separate dataset (may be too big for an attribute)
  if (called_from_checkpoint or not CCTK_Equals(out_save_parameters, "no")) {
    const int get_all =
        called_from_checkpoint or CCTK_Equals(out_save_parameters, "all");
    char *parameters = IOUtil_GetAllParameters(cctkGH, get_all);
    assert(parameters);
    error_count += WriteLargeAttribute(group, ALL_PARAMETERS, parameters);
    free(parameters);
  }

  // Save grid structure
  if (called_from_checkpoint or not CCTK_Equals(out_save_parameters, "no")) {
    vector<vector<vector<region_t> > > grid_superstructure(maps);
    vector<vector<vector<region_t> > > grid_structure(maps);
    vector<vector<i2vect> > grid_ghosts(maps);
    vector<vector<i2vect> > grid_buffers(maps);
    vector<vector<int> > grid_prolongation_orders(maps);
    for (int m = 0; m < maps; ++m) {
      grid_superstructure.at(m) = vhh.at(m)->superregions;
      grid_structure.at(m) = vhh.at(m)->regions.at(0);
      grid_ghosts.at(m) = vdd.at(m)->ghost_widths;
      grid_buffers.at(m) = vdd.at(m)->buffer_widths;
      grid_prolongation_orders.at(m) = vdd.at(m)->prolongation_orders_space;
    }
    vector<vector<vector<CCTK_REAL> > > grid_times(mglevels);
    for (int ml = 0; ml < mglevels; ++ml) {
      grid_times.at(ml).resize(vhh.at(0)->reflevels());
      for (int rl = 0; rl < vhh.at(0)->reflevels(); ++rl) {
        grid_times.at(ml).at(rl).resize(tt->timelevels);
        for (int tl = 0; tl < tt->timelevels; ++tl) {
          grid_times.at(ml).at(rl).at(tl) = tt->get_time(ml, rl, tl);
        }
      }
    }
    vector<vector<CCTK_REAL> > grid_delta_times(mglevels);
    for (int ml = 0; ml < mglevels; ++ml) {
      grid_delta_times.at(ml).resize(vhh.at(0)->reflevels());
      for (int rl = 0; rl < vhh.at(0)->reflevels(); ++rl) {
        grid_delta_times.at(ml).at(rl) = tt->get_delta(ml, rl);
      }
    }
    ostringstream gs_buf;
    gs_buf << setprecision(17);
    // We could write this information only into one of the checkpoint
    // files (to save space), or write it into a separate metadata
    // file
    gs_buf << "grid_superstructure:" << grid_superstructure << ",";
    // We could omit the grid structure (to save space), or write it
    // only into one of the checkpoint files
    gs_buf << "grid_structure:" << grid_structure << ",";
    gs_buf << "grid_times:" << grid_times << ",";
    gs_buf << "grid_delta_times:" << grid_delta_times << ",";
    // gs_buf << "grid_leveltimes:" << leveltimes << ",";
    gs_buf << "grid_ghosts:" << grid_ghosts << ",";
    gs_buf << "grid_buffers:" << grid_buffers << ",";
    gs_buf << "grid_prolongation_orders:" << grid_prolongation_orders << ".";
    string const gs_str = gs_buf.str();
    error_count += WriteLargeAttribute(group, GRID_STRUCTURE, gs_str.c_str());
  }

  HDF5_ERROR(H5Gclose(group));

  return (error_count);
}

} // namespace CarpetIOHDF5
