#include <cassert>
#include <cctype>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <list>
#include <sstream>
#include <string>
#include <vector>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <util_Network.h>

#include <Timer.hh>

#include "CactusBase/IOUtil/src/ioGH.h"
#include "CactusBase/IOUtil/src/ioutil_Utils.h"

#include "carpet.hh"

#include "typeprops.hh"

// That's a hack
namespace Carpet {
void UnsupportedVarType(const int vindex);
}

namespace CarpetIOScalar {

using namespace std;
using namespace Carpet;

// Definition of local types
struct info {
  string reduction;
  int handle;
};

// Begin a new line without flushing the output buffer
char const *const eol = "\n";

// Registered functions
static void *SetupGH(tFleshConfig *fc, int convLevel, cGH *cctkGH);
static int OutputGH(const cGH *cctkGH);
static int OutputVarAs(const cGH *cctkGH, const char *varname,
                       const char *alias);
static int OutputVarAs(const cGH *cctkGH, const char *varname,
                       const char *alias, const char *out_reductions);
static int TimeToOutput(const cGH *cctkGH, int vindex);
static int TriggerOutput(const cGH *cctkGH, int vindex);

// Internal functions
#if 0
  static void SetFlag (int index, const char* optstring, void* arg);
#endif
static void CheckSteerableParameters(const cGH *const cctkGH,
                                     bool first_time = false);

// Definition of static members
vector<bool> do_truncate;
vector<bool> reductions_changed;
vector<int> last_output;

/* CarpetScalar GH extension structure */
static struct {
  /* list of variables to output */
  char *out_vars;

  /* reductions to apply */
  char *out_reductions;
  char **var_reductions;

  /* stop on I/O parameter parsing errors ? */
  int stop_on_parse_errors;

  /* I/O request description list (for all variables) */
  ioRequest **requests;
} IOparameters;

extern "C" int CarpetIOScalarStartup() {
  CCTK_RegisterBanner("AMR scalar I/O provided by CarpetIOScalar");

  int GHExtension = CCTK_RegisterGHExtension("CarpetIOScalar");
  CCTK_RegisterGHExtensionSetupGH(GHExtension, SetupGH);

  int IOMethod = CCTK_RegisterIOMethod("CarpetIOScalar");
  CCTK_RegisterIOMethodOutputGH(IOMethod, OutputGH);
  CCTK_RegisterIOMethodOutputVarAs(IOMethod, OutputVarAs);
  CCTK_RegisterIOMethodTimeToOutput(IOMethod, TimeToOutput);
  CCTK_RegisterIOMethodTriggerOutput(IOMethod, TriggerOutput);

  return 0;
}

extern "C" void CarpetIOScalarInit(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;

  *this_iteration = 0;
  *last_output_iteration = 0;
  *last_output_time = cctkGH->cctk_time;
}

void *SetupGH(tFleshConfig *const fc, int const convLevel, cGH *const cctkGH) {
  DECLARE_CCTK_PARAMETERS;
  const void *dummy;

  dummy = &fc;
  dummy = &convLevel;
  dummy = &cctkGH;
  dummy = &dummy;

  // Truncate all files if this is not a restart
  const int numvars = CCTK_NumVars();
  do_truncate.resize(numvars, true);
  reductions_changed.resize(numvars, false);

  // No iterations have yet been output
  last_output.resize(numvars, -1);

  IOparameters.requests = (ioRequest **)calloc(numvars, sizeof(ioRequest *));
  IOparameters.out_vars = strdup("");
  IOparameters.out_reductions = strdup("");
  IOparameters.var_reductions = (char **)calloc(numvars, sizeof(char *));

  // initial I/O parameter check
  IOparameters.stop_on_parse_errors = strict_io_parameter_check;
  CheckSteerableParameters(cctkGH, true);
  IOparameters.stop_on_parse_errors = 0;

  // We register only once, ergo we get only one handle.  We store
  // that statically, so there is no need to pass anything to
  // Cactus.
  return NULL;
}

int OutputGH(const cGH *const cctkGH) {
  static Timers::Timer timer("OutputGH");
  timer.start();
  CheckSteerableParameters(cctkGH);
  if (strcmp(IOparameters.out_vars, "")) {
    for (int vindex = 0; vindex < CCTK_NumVars(); ++vindex) {
      if (TimeToOutput(cctkGH, vindex)) {
        TriggerOutput(cctkGH, vindex);
      }
    }
  }
  timer.stop();
  return 0;
}

int OutputVarAs(const cGH *const cctkGH, const char *const varname,
                const char *const alias) {
  DECLARE_CCTK_PARAMETERS;

  int const retval = OutputVarAs(cctkGH, varname, alias, outScalar_reductions);

  return retval;
}

int OutputVarAs(const cGH *const cctkGH, const char *const varname,
                const char *const alias, const char *out_reductions) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  assert(is_level_mode() or (is_singlemap_mode() and maps == 1) or
         (is_local_mode() and maps == 1 and
          vhh.at(Carpet::map)->local_components(reflevel) == 1));
  BEGIN_LEVEL_MODE(cctkGH) {

    const int n = CCTK_VarIndex(varname);
    if (n < 0) {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Variable \"%s\" does not exist", varname);
      return -1;
    }
    assert(n >= 0 and n < CCTK_NumVars());
    const int group = CCTK_GroupIndexFromVarI(n);
    assert(group >= 0 and group < (int)Carpet::arrdata.size());
    const int n0 = CCTK_FirstVarIndexI(group);
    assert(n0 >= 0 and n0 < CCTK_NumVars());
    const int var = n - n0;
    assert(var >= 0 and var < CCTK_NumVarsInGroupI(group));
    const int num_tl = CCTK_MaxActiveTimeLevelsVI(cctkGH, n);
    assert(num_tl >= 1);

    // Check for storage
    if (not CCTK_QueryGroupStorageI(cctkGH, group)) {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Cannot output variable \"%s\" because it has no storage",
                 varname);
      return 0;
    }

    assert(do_global_mode);

    const int vartype = CCTK_VarTypeI(n);
    assert(vartype >= 0);

    // Get grid hierarchy extentsion from IOUtil
    const ioGH *const iogh = (const ioGH *)CCTK_GHExtension(cctkGH, "IO");
    assert(iogh);

    // Create the output directory
    const char *myoutdir = outScalar_dir;
    if (CCTK_EQUALS(myoutdir, "")) {
      myoutdir = out_dir;
    }
    int const iret = IOUtil_CreateDirectory(cctkGH, myoutdir, 0, 0);
    if (iret < 0) {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Could not create output directory \"%s\"", myoutdir);
    } else if (CCTK_Equals(verbose, "full")) {
      static bool firsttime = true;
      if (firsttime and iret > 0) {
        CCTK_VInfo(CCTK_THORNSTRING, "Output directory \"%s\" exists already",
                   myoutdir);
      } else if (not firsttime and iret == 0) {
        CCTK_VInfo(CCTK_THORNSTRING, "Created output directory \"%s\"",
                   myoutdir);
      }
      firsttime = false;
    }

    // Find the set of desired reductions
    list<info> reductions;
    string const redlist(out_reductions);
    string::const_iterator p = redlist.begin();
    while (p != redlist.end()) {
      while (p != redlist.end() and isspace(*p))
        ++p;
      if (p == redlist.end())
        break;
      string::const_iterator const start = p;
      while (p != redlist.end() and not isspace(*p))
        ++p;
      string::const_iterator const end = p;
      string const reduction(start, end);
      int const handle = CCTK_ReductionHandle(reduction.c_str());
      if (handle < 0) {
        CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Reduction operator \"%s\" does not exist (maybe there is "
                   "no reduction thorn active?)",
                   reduction.c_str());
      } else {
        info i;
        i.reduction = reduction;
        i.handle = handle;
        reductions.push_back(i);
      }
    }

    // Output in global mode
    BEGIN_GLOBAL_MODE(cctkGH) {

      // single fstreams object used for all output files.
      // This violates resource-allocation-is-initialization but is required
      // when outputting all reductions into a single file (in which case there
      // is only a single stream)
      fstream file;
      CCTK_REAL io_files = 0;
      CCTK_REAL io_bytes_begin = 0, io_bytes_end = 0;
      if (all_reductions_in_one_file) {
        BeginTimingIO(cctkGH);
      }

      for (list<info>::const_iterator ireduction = reductions.begin();
           ireduction != reductions.end(); ++ireduction) {
        string const reduction = ireduction->reduction;

        if (not all_reductions_in_one_file) {
          BeginTimingIO(cctkGH);
          io_files = 0;
          io_bytes_begin = 0;
          io_bytes_end = 0;
        }
        if (CCTK_MyProc(cctkGH) == 0) {

          bool created_file = false;

          if (not all_reductions_in_one_file || not file.is_open()) {
            // Invent a file name
            ostringstream filenamebuf;
            filenamebuf << myoutdir << "/" << alias;
            if (not all_reductions_in_one_file)
              filenamebuf << "." << reduction;
            else
              filenamebuf << ".scalars";
            filenamebuf << ".asc";
            // we need a persistent temporary here
            string filenamestr = filenamebuf.str();
            const char *const filename = filenamestr.c_str();

            if (do_truncate.at(n) and IO_TruncateOutputFiles(cctkGH)) {
              file.open(filename, ios::out | ios::trunc);
            } else {
              file.open(filename, ios::out | ios::app);
            }
            if (not file.good()) {
              CCTK_VWarn(
                  0, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Could not open output file \"%s\" for variable \"%s\"",
                  filename, varname);
            }
            assert(file.is_open());

            // first time we open this file if do_truncate is true
            // TODO: add the logic of DidOoutput from IOASCII
            created_file = do_truncate.at(n);

            io_files += 1;
          }
          io_bytes_begin = file.tellg();

          // write header if reductions changed or we create a new file
          if (created_file or reductions_changed.at(n)) {
            bool want_labels = false;
            bool want_date = false;
            bool want_parfilename = false;
            bool want_other = false;
            if (CCTK_EQUALS(out_fileinfo, "none")) {
              // do nothing
            } else if (CCTK_EQUALS(out_fileinfo, "axis labels")) {
              want_labels = true;
            } else if (CCTK_EQUALS(out_fileinfo, "creation date")) {
              want_date = true;
            } else if (CCTK_EQUALS(out_fileinfo, "parameter filename")) {
              want_parfilename = true;
            } else if (CCTK_EQUALS(out_fileinfo, "all")) {
              want_labels = true;
              want_date = true;
              want_parfilename = true;
              want_other = true;
            } else {
              CCTK_WARN(0, "internal error");
            }
            file << "# Scalar ASCII output created by CarpetIOScalar" << eol;
            if (want_date) {
              char run_host[1000];
              Util_GetHostName(run_host, sizeof run_host);
              char const *run_user = getenv("USER");
              if (not run_user) {
                run_user = "";
              }
              char run_date[1000];
              Util_CurrentDate(sizeof run_date, run_date);
              char run_time[1000];
              Util_CurrentTime(sizeof run_time, run_time);
              file << "# created on " << run_host << " by " << run_user
                   << " on " << run_date << " at " << run_time << eol;
            }
            if (want_parfilename) {
              char parameter_filename[10000];
              CCTK_ParameterFilename(sizeof parameter_filename,
                                     parameter_filename);
              file << "# parameter filename: \"" << parameter_filename << "\""
                   << eol;
            }
            if (want_other) {
              if (CCTK_IsFunctionAliased("UniqueBuildID")) {
                char const *const build_id =
                    (char const *)UniqueBuildID(cctkGH);
                file << "# Build ID: " << build_id << eol;
              }
              if (CCTK_IsFunctionAliased("UniqueSimulationID")) {
                char const *const job_id =
                    static_cast<char const *>(UniqueSimulationID(cctkGH));
                file << "# Simulation ID: " << job_id << eol;
              }
              if (CCTK_IsFunctionAliased("UniqueRunID")) {
                char const *const job_id =
                    static_cast<char const *>(UniqueRunID(cctkGH));
                file << "# Run ID: " << job_id << eol;
              }
            }
            file << "#" << eol;
            if (want_labels) {
              file << "# " << varname << " (" << alias << ")" << eol;
              file << "# 1:iteration 2:time 3:data" << eol;
              int col = 3;
              if (one_file_per_group or all_reductions_in_one_file) {
                file << "# data columns:";
                int const firstvar =
                    one_file_per_group ? CCTK_FirstVarIndexI(group) : n;
                int const numvars =
                    one_file_per_group ? CCTK_NumVarsInGroupI(group) : 1;
                list<info>::const_iterator first_reduction =
                    all_reductions_in_one_file ? reductions.begin()
                                               : ireduction;
                list<info>::const_iterator end_reduction =
                    all_reductions_in_one_file
                        ? reductions.end()
                        : ++list<info>::const_iterator(ireduction);
                for (list<info>::const_iterator jreduction = first_reduction;
                     jreduction != end_reduction; ++jreduction) {
                  for (int n = firstvar; n < firstvar + numvars; ++n) {
                    file << " " << col << ":" << CCTK_VarName(n);
                    if (all_reductions_in_one_file)
                      file << "(" << jreduction->reduction << ")";
                    col += CarpetSimpleMPIDatatypeLength(vartype);
                  }
                }
                file << eol;
              }
            }

            // Don't write header again (unless the reductions change)
            reductions_changed.at(n) = false;
          }

          file << setprecision(out_precision);
          assert(file.good());

        } // if on the root processor

        if (CCTK_MyProc(cctkGH) == 0) {
          if (not all_reductions_in_one_file or
              ireduction == reductions.begin()) {
            file << cctk_iteration << " " << cctk_time;
          }
        }

        int const handle = ireduction->handle;

        char result[100]; // assuming no type is larger

        int const firstvar =
            one_file_per_group ? CCTK_FirstVarIndexI(group) : n;
        int const numvars =
            one_file_per_group ? CCTK_NumVarsInGroupI(group) : 1;

        for (int n = firstvar; n < firstvar + numvars; ++n) {

          int const ierr =
              CCTK_Reduce(cctkGH, 0, handle, 1, vartype, result, 1, n);
          if (ierr) {
            char *const fullname = CCTK_FullName(n);
            CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                       "Error during reduction for variable \"%s\"", fullname);
            free(fullname);
            memset(result, 0, sizeof result);
          }

          if (CCTK_MyProc(cctkGH) == 0) {
            file << " ";

            switch (specific_cactus_type(vartype)) {
#define CARPET_NO_COMPLEX
#define TYPECASE(N, T)                                                         \
  case N:                                                                      \
    file << *(T const *)result;                                                \
    break;
#include "typecase.hh"
#undef TYPECASE
#undef CARPET_NO_COMPLEX
#define CARPET_COMPLEX
#define TYPECASE(N, T)                                                         \
  case N:                                                                      \
    file << real(*(T const *)result) << " " << imag(*(T const *)result);       \
    break;
#include "typecase.hh"
#undef TYPECASE
#undef CARPET_COMPLEX
            default:
              UnsupportedVarType(n);
            }
          }

        } // for n

        if (not all_reductions_in_one_file) {
          if (CCTK_MyProc(cctkGH) == 0) {
            file << eol;
            assert(file.good());

            io_bytes_end = file.tellg();
            file.close();
            assert(file.good());
          }

          assert(not file.is_open());

          CCTK_REAL const io_bytes = io_bytes_end - io_bytes_begin;
          EndTimingIO(cctkGH, io_files, io_bytes, false);
        }

      } // for reductions

      if (all_reductions_in_one_file) {
        if (CCTK_MyProc(cctkGH) == 0) {
          file << eol;
          assert(file.good());

          io_bytes_end = file.tellg();
          file.close();
          assert(file.good());
        }

        assert(not file.is_open());

        CCTK_REAL const io_bytes = io_bytes_end - io_bytes_begin;
        EndTimingIO(cctkGH, io_files, io_bytes, false);
      }

      assert(not file.is_open());
    }
    END_GLOBAL_MODE;

    // Don't truncate again, this is a per-variable setting, not a
    // per-file setting so it cannot be moved inside of the loop over
    // reductions
    do_truncate.at(n) = false;
  }
  END_LEVEL_MODE;

  return 0;
}

int TimeToOutput(const cGH *const cctkGH, int const vindex) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  assert(vindex >= 0 and vindex < CCTK_NumVars());

  if (not do_global_mode)
    return 0;

  // check if output for this variable was requested
  if (not IOparameters.requests[vindex]) {
    return (0);
  }

  // check whether to output at this iteration
  bool output_this_iteration;

  const char *myoutcriterion = outScalar_criterion;
  if (CCTK_EQUALS(myoutcriterion, "default")) {
    myoutcriterion = out_criterion;
  }

  if (CCTK_EQUALS(myoutcriterion, "never")) {

    // Never output
    output_this_iteration = false;

  } else if (CCTK_EQUALS(myoutcriterion, "iteration")) {

    int myoutevery = outScalar_every;
    if (myoutevery == -2) {
      myoutevery = out_every;
    }
    if (myoutevery <= 0) {
      // output is disabled
      output_this_iteration = false;
    } else if (cctk_iteration == *this_iteration) {
      // we already decided to output this iteration
      output_this_iteration = true;
    } else if (cctk_iteration >= *last_output_iteration + myoutevery) {
      // it is time for the next output
      output_this_iteration = true;
      *last_output_iteration = cctk_iteration;
      *this_iteration = cctk_iteration;
    } else {
      // we want no output at this iteration
      output_this_iteration = false;
    }

  } else if (CCTK_EQUALS(myoutcriterion, "divisor")) {

    int myoutevery = outScalar_every;
    if (myoutevery == -2) {
      myoutevery = out_every;
    }
    if (myoutevery <= 0) {
      // output is disabled
      output_this_iteration = false;
    } else if ((cctk_iteration % myoutevery) == 0) {
      // we already decided to output this iteration
      output_this_iteration = true;
    } else {
      // we want no output at this iteration
      output_this_iteration = false;
    }

  } else if (CCTK_EQUALS(myoutcriterion, "time")) {

    CCTK_REAL myoutdt = outScalar_dt;
    if (myoutdt == -2) {
      myoutdt = out_dt;
    }
    if (myoutdt < 0) {
      // output is disabled
      output_this_iteration = false;
    } else if (myoutdt == 0) {
      // output all iterations
      output_this_iteration = true;
    } else if (cctk_iteration == *this_iteration) {
      // we already decided to output this iteration
      output_this_iteration = true;
    } else {
      int do_output =
          (cctk_time / cctk_delta_time >=
           (*last_output_time + myoutdt) / cctk_delta_time - 1.0e-12);
      MPI_Bcast(&do_output, 1, MPI_INT, 0, dist::comm());
      if (do_output) {
        // it is time for the next output
        output_this_iteration = true;
        *last_output_time = cctk_time;
        *this_iteration = cctk_iteration;
      } else {
        // we want no output at this iteration
        output_this_iteration = false;
      }
    }

  } else {

    assert(0);

  } // select output criterion

  if (not output_this_iteration)
    return 0;

#if 0
    // check which variables to output
    static vector<bool> output_variables;
    static int output_variables_iteration = -1;

    if (cctk_iteration > output_variables_iteration) {
      output_variables.resize (CCTK_NumVars());

      const char* const varlist = outScalar_vars;
      if (CCTK_TraverseString (varlist, SetFlag, &output_variables,
                               CCTK_GROUP_OR_VAR) < 0)
      {
        CCTK_WARN (output_variables_iteration < 0 and strict_io_parameter_check ?
                   0 : 1,
                   "error while parsing parameter 'IOScalar::outScalar_vars'");
      }

      output_variables_iteration = cctk_iteration;
    }

    if (not output_variables.at(vindex)) return 0;
#endif

  if (last_output.at(vindex) == cctk_iteration) {
    // Has already been output during this iteration
    char *const varname = CCTK_FullName(vindex);
    CCTK_VWarn(5, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Skipping output for variable \"%s\", because this variable "
               "has already been output during the current iteration -- "
               "probably via a trigger during the analysis stage",
               varname);
    free(varname);
    return 0;
  }

  assert(last_output.at(vindex) < cctk_iteration);

  // Should be output during this iteration
  return 1;
}

int TriggerOutput(const cGH *const cctkGH, int const vindex) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  assert(vindex >= 0 and vindex < CCTK_NumVars());

  // use individual reductions list for this variable when given
  // otherwise IOScalar::outScalar_reductions
  const char *out_reductions = IOparameters.requests[vindex]->reductions;
  if (not out_reductions)
    out_reductions = outScalar_reductions;

  int retval;

  if (one_file_per_group) {

    char *const fullname = CCTK_FullName(vindex);
    int const gindex = CCTK_GroupIndexFromVarI(vindex);
    char *const groupname_c = CCTK_GroupName(gindex);
    string groupname(groupname_c);
    transform(groupname.begin(), groupname.end(), groupname.begin(), ::tolower);
    string const oldsep("::");
    size_t const oldseppos = groupname.find(oldsep);
    assert(oldseppos != string::npos);
    groupname.replace(oldseppos, oldsep.size(), out_group_separator);
    retval = OutputVarAs(cctkGH, fullname, groupname.c_str(), out_reductions);
    free(fullname);

    int const firstvar = CCTK_FirstVarIndexI(gindex);
    int const numvars = CCTK_NumVarsInGroupI(gindex);
    for (int n = firstvar; n < firstvar + numvars; ++n) {
      last_output.at(n) = cctk_iteration;
    }

  } else {

    char *const fullname = CCTK_FullName(vindex);
    char const *const varname = CCTK_VarName(vindex);
    retval = OutputVarAs(cctkGH, fullname, varname, out_reductions);
    free(fullname);

    last_output.at(vindex) = cctk_iteration;
  }

  return retval;
}

static void CheckSteerableParameters(const cGH *const cctkGH, bool first_time) {
  DECLARE_CCTK_PARAMETERS;

  const int numvars = CCTK_NumVars();

  // re-parse the 'IOScalar::outScalar_vars' parameter if it has changed
  if (strcmp(outScalar_vars, IOparameters.out_vars)) {
#ifdef IOUTIL_PARSER_HAS_OUT_DT
    IOUtil_ParseVarsForOutput(cctkGH, CCTK_THORNSTRING,
                              "IOScalar::outScalar_vars",
                              IOparameters.stop_on_parse_errors, outScalar_vars,
                              -1, -1.0, IOparameters.requests);
#else
    IOUtil_ParseVarsForOutput(cctkGH, CCTK_THORNSTRING,
                              "IOScalar::outScalar_vars",
                              IOparameters.stop_on_parse_errors, outScalar_vars,
                              -1, IOparameters.requests);
#endif

    // notify the user about the new setting
    if (not CCTK_Equals(verbose, "none")) {
      int count = 0;
      ostringstream msg;
      msg << "Periodic scalar output requested for:";
      for (int vi = 0; vi < CCTK_NumVars(); ++vi) {
        if (IOparameters.requests[vi]) {
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

    // save the last setting of 'IOScalar::outScalar_vars' parameter
    free(IOparameters.out_vars);
    IOparameters.out_vars = strdup(outScalar_vars);

    for (int vi = 0; vi < numvars; ++vi) {
      if (not IOparameters.requests[vi])
        continue;

      // if the setting differ copy the setting (NULL or string) into
      // IOparameters
      if ((IOparameters.requests[vi]->reductions and
           not IOparameters.var_reductions[vi]) or
          (not IOparameters.requests[vi]->reductions and
           IOparameters.var_reductions[vi]) or
          (IOparameters.requests[vi]->reductions and
           IOparameters.var_reductions[vi] and
           strcmp(IOparameters.requests[vi]->reductions,
                  IOparameters.var_reductions[vi]))) {
        if (not IOparameters.var_reductions[vi]) {
          free(IOparameters.var_reductions[vi]);
        }
        IOparameters.var_reductions[vi] =
            IOparameters.requests[vi]->reductions
                ? strdup(IOparameters.requests[vi]->reductions)
                : NULL;
        // only all_reductions_in_one_file actually mentions the reduction in
        // the header
        // we must ignore the initial change when we first learn of
        // requested reductions. This means no header is generated if
        // a checkpoint-recovery changes the reductions.
        // TODO: re-write header in one_file_per_group mode if variables change
        reductions_changed[vi] = all_reductions_in_one_file and not first_time;
      }
    }
  }

  if (strcmp(outScalar_reductions, IOparameters.out_reductions)) {
    // save the last seeting of 'IOScalar::outScalar_reductions' parameter
    free(IOparameters.out_reductions);
    IOparameters.out_reductions = strdup(outScalar_reductions);
    // bit of an overkill. We ask for new headers for all variables, even though
    // the ones using per-variable reductions will not have differing headers
    for (int vi = 0; vi < numvars; ++vi) {
      // only all_reductions_in_one_file actually mentions the reduction in the
      // header
      // we must ignore the initial change when we first learn of
      // requested reductions. This means no header is generated if
      // a checkpoint-recovery changes the reductions.
      reductions_changed[vi] = all_reductions_in_one_file and not first_time;
    }
  }
}

#if 0
  void
  SetFlag (int const index, const char * const optstring, void * const arg)
  {
    const void *dummy;

    dummy = &optstring;
    dummy = &dummy;

    vector<bool>& flags = *(vector<bool>*)arg;
    flags.at(index) = true;
  }
#endif

} // namespace CarpetIOScalar
