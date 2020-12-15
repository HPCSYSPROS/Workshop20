#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <list>
#include <sstream>
#include <string>
#include <vector>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include "CactusBase/IOUtil/src/ioGH.h"
#include "CactusBase/IOUtil/src/ioutil_Utils.h"

#include <Timer.hh>

#include "carpet.hh"

#include "typeprops.hh"

using namespace CarpetLib;

// That's a hack
namespace Carpet {
void UnsupportedVarType(const int vindex);
}

namespace CarpetIOBasic {

using namespace std;
using namespace Carpet;

// Definition of local types
struct info {
  string reduction;
  int handle;
};

// Registered functions
void *SetupGH(tFleshConfig *fc, int convLevel, cGH *cctkGH);
int OutputGH(const cGH *cctkGH);

// Internal functions
bool TimeToOutput(const cGH *cctkGH);
void OutputHeader(const cGH *cctkGH);
void OutputVar(const cGH *cctkGH, int vindex);
void OutputVar(const cGH *cctkGH, int vindex, const char *out_reductions);
void CheckSteerableParameters(const cGH *cctkGH);

void ExamineVariable(int vindex, bool &isint, int &numcomps, bool &isscalar);
vector<string> ParseReductions(char const *credlist);

template <typename T> bool UseScientificNotation(T const &x) {
  return false; // default
}

template <> bool UseScientificNotation(CCTK_REAL const &x) {
  DECLARE_CCTK_PARAMETERS;

  CCTK_REAL const xa = fabs(x);
  return xa != 0 and (xa < real_min or xa >= real_max);
}

template <> bool UseScientificNotation(CCTK_COMPLEX const &x) {
  return UseScientificNotation(x.real()) or UseScientificNotation(x.imag());
}

// Definition of members
int output_count;
int last_output;

/* CarpetBasic GH extension structure */
static struct {
  /* list of variables to output */
  char *out_vars;

  /* stop on I/O parameter parsing errors ? */
  bool stop_on_parse_errors;

  /* I/O request description list (for all variables) */
  ioRequest **requests;
} IOparameters;

extern "C" int CarpetIOBasicStartup() {
  CCTK_RegisterBanner("AMR info I/O provided by CarpetIOBasic");

  int GHExtension = CCTK_RegisterGHExtension("CarpetIOBasic");
  CCTK_RegisterGHExtensionSetupGH(GHExtension, SetupGH);

  int IOMethod = CCTK_RegisterIOMethod("CarpetIOBasic");
  CCTK_RegisterIOMethodOutputGH(IOMethod, OutputGH);

  return 0;
}

extern "C" void CarpetIOBasicInit(CCTK_ARGUMENTS) {
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
  output_count = 0;

  // No iterations have yet been output
  last_output = -1;

  IOparameters.requests = (ioRequest **)calloc(numvars, sizeof(ioRequest *));
  IOparameters.out_vars = strdup("");

  // initial I/O parameter check
  IOparameters.stop_on_parse_errors = strict_io_parameter_check;
  CheckSteerableParameters(cctkGH);
  IOparameters.stop_on_parse_errors = false;

  // We register only once, ergo we get only one handle.  We store
  // that statically, so there is no need to pass anything to
  // Cactus.
  return NULL;
}

int OutputGH(const cGH *const cctkGH) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  static Timers::Timer timer("OutputGH");
  timer.start();

  if (TimeToOutput(cctkGH)) {

    int const oldprec = cout.precision();
    ios_base::fmtflags const oldflags = cout.flags();

    if (output_count++ % outHeader_every == 0 && outHeader_every != -1) {
      // Print the header
      OutputHeader(cctkGH);
    }

    if (CCTK_MyProc(cctkGH) == 0) {
      cout << setw(iter_width) << setfill(' ') << cctk_iteration << " "
           << setw(time_width) << setfill(' ') << fixed
           << setprecision(time_prec) << cctk_time;
    }

    int const numvars = CCTK_NumVars();
    for (int vindex = 0; vindex < numvars; ++vindex) {
      if (IOparameters.requests[vindex]) {
        // use individual reductions list for this variable when given
        // otherwise IOScalar::outInfo_reductions
        const char *out_reductions = IOparameters.requests[vindex]->reductions;
        if (not out_reductions)
          out_reductions = outInfo_reductions;

        OutputVar(cctkGH, vindex, out_reductions);
      }
    }

    if (CCTK_MyProc(cctkGH) == 0) {
      cout << endl << flush;
    }

    last_output = cctk_iteration;

    cout.precision(oldprec);
    cout.setf(oldflags);

  } // if time to output

  timer.stop();

  return 0;
}

void OutputHeader(const cGH *const cctkGH) {
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_MyProc(cctkGH) == 0) {

    int const numvars = CCTK_NumVars();

    // The header consists of four lines:
    //    pass 0: print separator
    //    pass 1: print variable name
    //    pass 2: print reduction operation
    //    pass 3: print separator
    for (int pass = 0; pass < 4; ++pass) {

      // Print iteration and time
      switch (pass) {
      case 0:
      case 3:
        cout << "-------------------";
        break;
      case 1:
        cout << "Iteration      Time";
        break;
      case 2:
        cout << "                   ";
        break;
      default:
        assert(0);
      }

      // Loop over all variables that should be output
      for (int vindex = 0; vindex < numvars; ++vindex) {
        if (IOparameters.requests[vindex]) {

          bool isint;
          int numcomps;
          bool isscalar;
          ExamineVariable(vindex, isint, numcomps, isscalar);

          char *cfullname = CCTK_FullName(vindex);
          string const fullname(cfullname);
          free(cfullname);

          // Print a vertical separator
          switch (pass) {
          case 0:
          case 3:
            cout << "--";
            break;
          case 1:
          case 2:
            cout << " |";
            break;
          default:
            assert(0);
          }

          // use individual reductions list for this variable when given
          // otherwise IOScalar::outInfo_reductions
          const char *out_reductions =
              IOparameters.requests[vindex]->reductions;
          if (not out_reductions)
            out_reductions = outInfo_reductions;

          // Find the set of desired reductions
          vector<string> const reductions = ParseReductions(out_reductions);
          int const numreds = reductions.size();

          int const width = isint ? int_width : real_width;

          int const mynumreds = isscalar ? 1 : numreds;

          // Print the entry
          switch (pass) {
          case 0:
          case 3: {
            size_t const numchars = (width + 1) * numcomps * mynumreds;
            cout << setw(numchars) << setfill('-') << "";
          } break;
          case 1: {
            size_t const numchars = (width + 1) * numcomps * mynumreds - 1;
            if (fullname.length() > numchars) {
              int begin = fullname.length() - (numchars - 1);
              cout << " *" << fullname.substr(begin);
            } else {
              cout << " " << setw(numchars) << setfill(' ') << fullname;
            }
          } break;
          case 2:
            if (isscalar) {
              int const numchars = (width + 1) * numcomps;
              cout << setw(numchars) << setfill(' ') << "";
            } else {
              for (int red = 0; red < mynumreds; ++red) {
                int const numchars = (width + 1) * numcomps - 1;
                cout << " " << setw(numchars) << setfill(' ')
                     << reductions.at(red).substr(0, numchars);
              }
            }
            break;
          default:
            assert(0);
          }

        } // if want output
      }   // for vindex

      cout << endl;

    } // pass

  } // if on root processor
}

void OutputVar(const cGH *const cctkGH, int const n) {
  DECLARE_CCTK_PARAMETERS;

  OutputVar(cctkGH, n, outInfo_reductions);
}

void OutputVar(const cGH *const cctkGH, int n, const char *out_reductions) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  assert(is_level_mode());

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
    char *fullname = CCTK_FullName(n);
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Cannot output variable \"%s\" because it has no storage",
               fullname);
    free(fullname);
    return;
  }

  assert(do_global_mode);

  const int vartype = CCTK_VarTypeI(n);
  assert(vartype >= 0);

  // Get grid hierarchy extentsion from IOUtil
  const ioGH *const iogh = (const ioGH *)CCTK_GHExtension(cctkGH, "IO");
  assert(iogh);

  // Find the set of desired reductions
  vector<string> const reductionstrings = ParseReductions(out_reductions);
  list<info> reductions;
  for (vector<string>::const_iterator ireduction = reductionstrings.begin();
       ireduction != reductionstrings.end(); ++ireduction) {
    int const handle = CCTK_ReductionHandle((*ireduction).c_str());
    if (handle < 0) {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Reduction operator \"%s\" does not exist (maybe there is no "
                 "reduction thorn active?)",
                 (*ireduction).c_str());
    } else {
      info i;
      i.reduction = *ireduction;
      i.handle = handle;
      reductions.push_back(i);
    }
  }

  bool isint;
  int numcomps;
  bool isscalar;
  ExamineVariable(n, isint, numcomps, isscalar);

  // Output in global mode
  BEGIN_GLOBAL_MODE(cctkGH) {

    // Remember cout state
    int const oldprec = cout.precision();
    ios_base::fmtflags const oldflags = cout.flags();

    // Print vertical separator
    if (CCTK_MyProc(cctkGH) == 0) {
      cout << " |";
    }

    int const width = isint ? int_width : real_width;

    if (isscalar) {

      if (CCTK_MyProc(cctkGH) == 0) {

        cout << " " << setw(width);

        void const *const vardataptr = CCTK_VarDataPtrI(cctkGH, 0, n);
        assert(vardataptr);

        switch (specific_cactus_type(vartype)) {
#define CARPET_NO_COMPLEX
#define TYPECASE(N, T)                                                         \
  case N: {                                                                    \
    T const val = *static_cast<T const *>(vardataptr);                         \
    if (not isint) {                                                           \
      if (UseScientificNotation(val)) {                                        \
        cout << scientific << setprecision(real_prec_sci);                     \
      } else {                                                                 \
        cout << fixed << setprecision(real_prec);                              \
      }                                                                        \
    }                                                                          \
    cout << val;                                                               \
  } break;
#include "typecase.hh"
#undef TYPECASE
#undef CARPET_NO_COMPLEX
#define CARPET_COMPLEX
#define TYPECASE(N, T)                                                         \
  case N: {                                                                    \
    T const val = *static_cast<T const *>(vardataptr);                         \
    if (not isint) {                                                           \
      if (UseScientificNotation(val)) {                                        \
        cout << scientific << setprecision(real_prec_sci);                     \
      } else {                                                                 \
        cout << fixed << setprecision(real_prec);                              \
      }                                                                        \
    }                                                                          \
    cout << real(val) << " " << imag(val);                                     \
  } break;
#include "typecase.hh"
#undef TYPECASE
#undef CARPET_COMPLEX
        default:
          UnsupportedVarType(n);
        }

      } // if on root processor

    } else {

      for (list<info>::const_iterator ireduction = reductions.begin();
           ireduction != reductions.end(); ++ireduction) {
        string const reduction = ireduction->reduction;

        int const handle = ireduction->handle;

        char result[100]; // assuming no type is larger than this

        int const ierr =
            CCTK_Reduce(cctkGH, 0, handle, 1, vartype, result, 1, n);
        assert(not ierr);

        if (CCTK_MyProc(cctkGH) == 0) {

          cout << " " << setw(width);

          switch (specific_cactus_type(vartype)) {
#define TYPECASE(N, T)                                                         \
  case N: {                                                                    \
    T const &val = *(T const *)result;                                         \
    if (not isint) {                                                           \
      if (UseScientificNotation(val)) {                                        \
        cout << scientific << setprecision(real_prec_sci);                     \
      } else {                                                                 \
        cout << fixed << setprecision(real_prec);                              \
      }                                                                        \
    }                                                                          \
    cout << val;                                                               \
  } break;
#include "typecase.hh"
#undef TYPECASE
          default:
            UnsupportedVarType(n);
          }

        } // if on root processor

      } // for reductions

    } // not isscalar

    // Restore cout state
    cout.precision(oldprec);
    cout.setf(oldflags);
  }
  END_GLOBAL_MODE;
}

bool TimeToOutput(const cGH *const cctkGH) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (not do_global_mode)
    return false;

  CheckSteerableParameters(cctkGH);

  // check whether to output at this iteration
  bool output_this_iteration;

  const char *myoutcriterion = outInfo_criterion;
  if (CCTK_EQUALS(myoutcriterion, "default")) {
    myoutcriterion = out_criterion;
  }

  if (CCTK_EQUALS(myoutcriterion, "never")) {

    // Never output
    output_this_iteration = false;

  } else if (CCTK_EQUALS(myoutcriterion, "iteration")) {

    int myoutevery = outInfo_every;
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

    int myoutevery = outInfo_every;
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

    CCTK_REAL myoutdt = outInfo_dt;
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
    return false;

  // These should be true in general, but may be false if
  // CCTK_OutputGH is called explicitly:
  // assert (last_output != cctk_iteration);
  // assert (last_output < cctk_iteration);

  // Should be output during this iteration
  return true;
}

void CheckSteerableParameters(const cGH *const cctkGH) {
  DECLARE_CCTK_PARAMETERS;

  // re-parse the 'IOBasic::outInfo_vars' parameter if it has changed
  if (strcmp(outInfo_vars, IOparameters.out_vars)) {
#ifdef IOUTIL_PARSER_HAS_OUT_DT
    IOUtil_ParseVarsForOutput(cctkGH, CCTK_THORNSTRING, "IOBasic::outInfo_vars",
                              IOparameters.stop_on_parse_errors, outInfo_vars,
                              -1, -1.0, IOparameters.requests);
#else
    IOUtil_ParseVarsForOutput(cctkGH, CCTK_THORNSTRING, "IOBasic::outInfo_vars",
                              IOparameters.stop_on_parse_errors, outInfo_vars,
                              -1, IOparameters.requests);
#endif

    // save the last setting of 'IOBasic::outInfo_vars' parameter
    free(IOparameters.out_vars);
    IOparameters.out_vars = strdup(outInfo_vars);
  }
}

// Parse the set of reductions
vector<string> ParseReductions(char const *const credlist) {
  string const redlist(credlist);
  vector<string> reductions;
  string::const_iterator p = redlist.begin();
  while (p != redlist.end()) {
    // Skip white space
    while (p != redlist.end() and isspace(*p))
      ++p;
    // Exit if end of reductions is reached
    if (p == redlist.end())
      break;
    // Mark beginning of reduction entry
    string::const_iterator const start = p;
    // Walk over reduction entry
    while (p != redlist.end() and not isspace(*p))
      ++p;
    // Mark end of reduction entry
    string::const_iterator const end = p;
    // Remember reduction
    string const reduction(start, end);
    reductions.push_back(reduction);
  }
  return reductions;
}

void ExamineVariable(int const vindex, bool &isint, int &numcomps,
                     bool &isscalar) {
  switch (CCTK_VarTypeI(vindex)) {
  case CCTK_VARIABLE_BYTE:
  case CCTK_VARIABLE_INT:
  case CCTK_VARIABLE_INT1:
  case CCTK_VARIABLE_INT2:
  case CCTK_VARIABLE_INT4:
  case CCTK_VARIABLE_INT8:
  case CCTK_VARIABLE_INT16:
    isint = true;
    numcomps = 1;
    break;
  case CCTK_VARIABLE_REAL:
  case CCTK_VARIABLE_REAL4:
  case CCTK_VARIABLE_REAL8:
  case CCTK_VARIABLE_REAL16:
    isint = false;
    numcomps = 1;
    break;
  case CCTK_VARIABLE_COMPLEX:
  case CCTK_VARIABLE_COMPLEX8:
  case CCTK_VARIABLE_COMPLEX16:
  case CCTK_VARIABLE_COMPLEX32:
    isint = false;
    numcomps = 2;
    break;
  default:
    assert(0);
  }

  switch (CCTK_GroupTypeFromVarI(vindex)) {
  case CCTK_SCALAR:
    isscalar = true;
    break;
  case CCTK_ARRAY:
  case CCTK_GF:
    isscalar = false;
    break;
  default:
    assert(0);
  }
}

} // namespace CarpetIOBasic
