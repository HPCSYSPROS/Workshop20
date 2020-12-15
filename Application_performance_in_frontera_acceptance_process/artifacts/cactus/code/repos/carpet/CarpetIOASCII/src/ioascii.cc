#include <cassert>
#include <cctype>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <map>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

#include <cctk.h>
#include <cctk_Parameters.h>
#include <util_Network.h>
#include <util_Table.h>

#include <Timer.hh>

#include "CactusBase/IOUtil/src/ioGH.h"
#include "CactusBase/IOUtil/src/ioutil_Utils.h"

#include "carpet.hh"

#include "typeprops.hh"

#include "ioascii.hh"

// That's a hack
namespace Carpet {
void UnsupportedVarType(const int vindex);
}

#define GetParameter(parameter)                                                \
  (outdim == 0 ? out0D_##parameter : outdim == 1                               \
                                         ? out1D_##parameter                   \
                                         : outdim == 2 ? out2D_##parameter     \
                                                       : out3D_##parameter)

namespace CarpetIOASCII {

using namespace std;
using namespace Carpet;

// Begin a new line without flushing the output buffer
const char *const eol = "\n";

// IO processor
const int ioproc = 0;

// Global configuration parameters
bool stop_on_parse_errors = false;

int CarpetIOASCIIStartup() {
  IOASCII<0>::Startup();
  IOASCII<1>::Startup();
  IOASCII<2>::Startup();
  IOASCII<3>::Startup();
  return 0;
}

void CarpetIOASCIIInit(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;

  for (int d = 0; d < 4; ++d) {
    this_iteration[d] = 0;
    last_output_iteration[d] = 0;
    last_output_time[d] = cctk_time;
  }
}

// Definition of static members
template <int outdim> char *IOASCII<outdim>::my_out_dir;
template <int outdim> char *IOASCII<outdim>::my_out_vars;
template <int outdim> vector<ioRequest *> IOASCII<outdim>::requests;

template <int outdim> int IOASCII<outdim>::Startup() {
  ostringstream msg;
  msg << "AMR " << outdim << "D ASCII I/O provided by CarpetIOASCII";
  CCTK_RegisterBanner(msg.str().c_str());

  ostringstream extension_name;
  extension_name << "CarpetIOASCII_" << outdim << "D";
  const int GHExtension =
      CCTK_RegisterGHExtension(extension_name.str().c_str());
  CCTK_RegisterGHExtensionSetupGH(GHExtension, SetupGH);

  ostringstream method_name;
  method_name << "IOASCII_" << outdim << "D";
  const int IOMethod = CCTK_RegisterIOMethod(method_name.str().c_str());
  CCTK_RegisterIOMethodOutputGH(IOMethod, OutputGH);
  CCTK_RegisterIOMethodOutputVarAs(IOMethod, OutputVarAs);
  CCTK_RegisterIOMethodTimeToOutput(IOMethod, TimeToOutput);
  CCTK_RegisterIOMethodTriggerOutput(IOMethod, TriggerOutput);

  return 0;
}

template <int outdim>
void *IOASCII<outdim>::SetupGH(tFleshConfig *const fc, const int convLevel,
                               cGH *const cctkGH) {
  DECLARE_CCTK_PARAMETERS;
  const void *dummy;

  dummy = &fc;
  dummy = &convLevel;
  dummy = &cctkGH;
  dummy = &dummy;

  if (not CCTK_Equals(verbose, "none")) {
    CCTK_VInfo(CCTK_THORNSTRING,
               "I/O Method 'IOASCII_%dD' registered: "
               "%dD AMR output of grid variables to ASCII files",
               outdim, outdim);
  }

  const int numvars = CCTK_NumVars();
  requests.resize(numvars);

  // initial I/O parameter check
  my_out_dir = 0;
  my_out_vars = strdup("");
  stop_on_parse_errors = strict_io_parameter_check != 0;
  CheckSteerableParameters(cctkGH);
  stop_on_parse_errors = false;

  // We register only once, ergo we get only one handle.  We store
  // that statically, so there is no need to pass anything to
  // Cactus.
  return NULL;
}

template <int outdim>
void IOASCII<outdim>::CheckSteerableParameters(const cGH *const cctkGH) {
  DECLARE_CCTK_PARAMETERS;

  // re-parse the 'IOASCII::out%dD_dir' parameter if it has changed
  const char *the_out_dir = GetParameter(dir);
  if (CCTK_EQUALS(the_out_dir, "")) {
    the_out_dir = out_dir;
  }

  if (not my_out_dir or strcmp(the_out_dir, my_out_dir)) {
    free(my_out_dir);
    my_out_dir = strdup(the_out_dir);

    // create the output directory
    const int result = IOUtil_CreateDirectory(cctkGH, my_out_dir, 0, 0);
    if (result < 0) {
      CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Problem creating %dD-output directory '%s'", outdim,
                 my_out_dir);
    } else if (result > 0 and CCTK_Equals(verbose, "full")) {
      CCTK_VInfo(CCTK_THORNSTRING, "%dD-output directory '%s' already exists",
                 outdim, my_out_dir);
    }
  }

  // re-parse the 'IOASCII::out%d_vars' parameter if it has changed
  const char *const out_vars = GetParameter(vars);
  if (strcmp(out_vars, my_out_vars)) {
    ostringstream parameter_name;
    parameter_name << "IOASCII::out" << outdim << "D_vars";
#ifdef IOUTIL_PARSER_HAS_OUT_DT
    IOUtil_ParseVarsForOutput(
        cctkGH, CCTK_THORNSTRING, parameter_name.str().c_str(),
        stop_on_parse_errors, out_vars, -1, -1.0, &requests[0]);
#else
    IOUtil_ParseVarsForOutput(cctkGH, CCTK_THORNSTRING,
                              parameter_name.str().c_str(),
                              stop_on_parse_errors, out_vars, -1, &requests[0]);
#endif

    // notify the user about the new setting
    if (not CCTK_Equals(verbose, "none")) {
      int count = 0;
      ostringstream msg;
      msg << "Periodic " << outdim << "D AMR output requested for:";
      for (int vi = 0; vi < CCTK_NumVars(); ++vi) {
        if (requests.at(vi)) {
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

    // save the last setting of 'IOASCII::out%d_vars' parameter
    free(my_out_vars);
    my_out_vars = strdup(out_vars);
  }
}

template <int outdim> int IOASCII<outdim>::OutputGH(const cGH *const cctkGH) {
  DECLARE_CCTK_PARAMETERS;

  ostringstream timer_name;
  timer_name << "OutputGH<" << outdim << ">";
  Timers::Timer timer(timer_name.str());

  timer.start();
  CheckSteerableParameters(cctkGH);
  if (strcmp(my_out_vars, "")) {
    for (int vi = 0; vi < CCTK_NumVars(); ++vi) {
      if (TimeToOutput(cctkGH, vi)) {
        TriggerOutput(cctkGH, vi);
      }
    }
  }
  timer.stop();

  return 0;
}

template <int outdim>
int IOASCII<outdim>::TimeToOutput(const cGH *const cctkGH, const int vindex) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  assert(vindex >= 0 and vindex < CCTK_NumVars());

  if (CCTK_GroupTypeFromVarI(vindex) != CCTK_GF and not do_global_mode) {
    return 0;
  }

  // check if output for this variable was requested
  if (not requests.at(vindex)) {
    return 0;
  }

  // check whether this refinement level should be output
  if (not(requests.at(vindex)->refinement_levels & (1 << reflevel))) {
    return 0;
  }

  // check if output for this variable was requested individually by
  // a "<varname>{ out_every = <number> }" option string
  // this will overwrite the output criterion setting
  const char *myoutcriterion = GetParameter(criterion);
  if (CCTK_EQUALS(myoutcriterion, "default")) {
    myoutcriterion = out_criterion;
  }
  if (requests.at(vindex)->out_every >= 0) {
    myoutcriterion = "divisor";
  }

  if (CCTK_EQUALS(myoutcriterion, "never")) {
    return 0;
  }

  // check whether to output at this iteration
  bool output_this_iteration = false;

  if (CCTK_EQUALS(myoutcriterion, "iteration")) {
    int myoutevery = GetParameter(every);
    if (myoutevery == -2) {
      myoutevery = out_every;
    }
    if (myoutevery > 0) {
      if (cctk_iteration == this_iteration[outdim]) {
        // we already decided to output this iteration
        output_this_iteration = true;
      } else if (cctk_iteration >= last_output_iteration[outdim] + myoutevery) {
        // it is time for the next output
        output_this_iteration = true;
        last_output_iteration[outdim] = cctk_iteration;
        this_iteration[outdim] = cctk_iteration;
      }
    }
  } else if (CCTK_EQUALS(myoutcriterion, "divisor")) {
    int myoutevery = GetParameter(every);
    if (myoutevery == -2) {
      myoutevery = out_every;
    }
    if (requests[vindex]->out_every >= 0) {
      myoutevery = requests[vindex]->out_every;
    }
    if (myoutevery > 0 and (cctk_iteration % myoutevery) == 0) {
      // we already decided to output this iteration
      output_this_iteration = true;
    }
  } else if (CCTK_EQUALS(myoutcriterion, "time")) {
    CCTK_REAL myoutdt = GetParameter(dt);
    if (myoutdt == -2) {
      myoutdt = out_dt;
    }
    if (myoutdt == 0 or cctk_iteration == this_iteration[outdim]) {
      output_this_iteration = true;
    } else if (myoutdt > 0) {
      int do_output =
          cctk_time / cctk_delta_time >=
          (last_output_time[outdim] + myoutdt) / cctk_delta_time - 1.0e-12;
      MPI_Bcast(&do_output, 1, MPI_INT, 0, dist::comm());
      if (do_output) {
        // it is time for the next output
        output_this_iteration = true;
        last_output_time[outdim] = cctk_time;
        this_iteration[outdim] = cctk_iteration;
      }
    }
  } // select output criterion

  return output_this_iteration ? 1 : 0;
}

template <int outdim>
int IOASCII<outdim>::TriggerOutput(const cGH *const cctkGH, const int vindex) {
  DECLARE_CCTK_PARAMETERS;

  assert(vindex >= 0 and vindex < CCTK_NumVars());

  char *const fullname = CCTK_FullName(vindex);

  int retval;
  if (one_file_per_group) {
    char *const alias_c = CCTK_GroupNameFromVarI(vindex);
    string alias(alias_c);
    free(alias_c);
    transform(alias.begin(), alias.end(), alias.begin(), ::tolower);
    string const oldsep("::");
    size_t const oldseppos = alias.find(oldsep);
    assert(oldseppos != string::npos);
    alias.replace(oldseppos, oldsep.size(), out_group_separator);
    retval = OutputVarAs(cctkGH, fullname, alias.c_str());
  } else {
    const char *const alias = CCTK_VarName(vindex);
    retval = OutputVarAs(cctkGH, fullname, alias);
  }

  free(fullname);

  return retval;
}

static void GetVarIndex(const int vindex, const char *const optstring,
                        void *const arg) {
  if (optstring) {
    char *const fullname = CCTK_FullName(vindex);
    CCTK_VWarn(CCTK_WARN_COMPLAIN, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Option string '%s' will be ignored for ASCII output of "
               "variable '%s'",
               optstring, fullname);
    free(fullname);
  }

  *static_cast<int *>(arg) = vindex;
}

template <int outdim>
int IOASCII<outdim>::OutputVarAs(const cGH *const cctkGH,
                                 const char *const varname,
                                 const char *const alias) {
  DECLARE_CCTK_PARAMETERS;

  int vindex = -1;

  if (CCTK_TraverseString(varname, GetVarIndex, &vindex, CCTK_VAR) < 0) {
    CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,
               "error while parsing variable name '%s' (alias name '%s')",
               varname, alias);
    return -1;
  }

  if (vindex < 0) {
    return -1;
  }

  if (not(is_level_mode() or (is_singlemap_mode() and maps == 1) or
          (is_local_mode() and maps == 1 and
           vhh.at(Carpet::map)->local_components(reflevel) == 1))) {
    CCTK_WARN(CCTK_WARN_ALERT, "OutputVarAs must be called in level mode");
    return -1;
  }

  BEGIN_LEVEL_MODE(cctkGH) {

    // Get information
    const int group = CCTK_GroupIndexFromVarI(vindex);
    assert(group >= 0);
    cGroup groupdata;
    {
      int const ierr = CCTK_GroupData(group, &groupdata);
      assert(not ierr);
    }

    // Check information
    if (groupdata.grouptype != CCTK_GF) {
      assert(do_global_mode);
    }

    if (outdim > groupdata.dim) {
      CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Cannot produce %dD ASCII output file '%s' for variable '%s' "
                 "because it has only %d dimensions",
                 outdim, alias, varname, groupdata.dim);
      return -1;
    }

    // Check for storage
    if (not CCTK_QueryGroupStorageI(cctkGH, group)) {
      // This may be okay if storage is e.g. scheduled only in the
      // analysis bin
      CCTK_VWarn(CCTK_WARN_DEBUG, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Cannot output variable '%s' because it has no storage",
                 varname);
      return 0;
    }

    ostringstream basefilenamebuf;
    basefilenamebuf << my_out_dir << "/" << alias;
    const string basefilename = basefilenamebuf.str();

    // Check if the file has been created already
    bool is_new_file, truncate_file;
    const bool did_output =
        DidOutput(cctkGH, vindex, basefilename, is_new_file, truncate_file);
    if (did_output) {
      return 0;
    }

    // Loop over all direction combinations
    vect<int, outdim> dirs(0);
    bool done;
    do {

      // Output each combination only once
      bool ascending = true;
      for (int d1 = 0; d1 < outdim; ++d1) {
        for (int d2 = d1 + 1; d2 < outdim; ++d2) {
          ascending = ascending and dirs[d1] < dirs[d2];
        }
      }

      // Skip output if the dimensions are not ascending
      if (ascending) {

        // Skip output if not requested
        if (DirectionIsRequested(dirs)) {
          OutputDirection(cctkGH, vindex, alias, basefilename, dirs,
                          is_new_file, truncate_file);
        }

      } // if ascending

      // Next direction combination
      done = true;
      for (int d = 0; d < outdim; ++d) {
        ++dirs[d];
        if (dirs[d] < groupdata.dim + (outdim == 1 ? 1 : 0)) {
          done = false;
          break;
        }
        dirs[d] = 0;
      }

    } while (not done); // output all directions
  }
  END_LEVEL_MODE;

  return 0;
}

// Traverse all maps and components on this refinement and multigrid
// level
template <int outdim>
void IOASCII<outdim>::OutputDirection(const cGH *const cctkGH, const int vindex,
                                      const string alias,
                                      const string basefilename,
                                      const vect<int, outdim> &dirs,
                                      const bool is_new_file,
                                      const bool truncate_file) {
  DECLARE_CCTK_PARAMETERS;

  // Get information
  const int group = CCTK_GroupIndexFromVarI(vindex);
  assert(group >= 0);
  const int vindex0 = CCTK_FirstVarIndexI(group);
  assert(vindex0 >= 0 and vindex >= vindex0);
  const int var = vindex - vindex0;
  cGroup groupdata;
  {
    int const ierr = CCTK_GroupData(group, &groupdata);
    assert(not ierr);
  }

  const int ml = groupdata.grouptype == CCTK_GF ? mglevel : 0;
  const int rl = groupdata.grouptype == CCTK_GF ? reflevel : 0;

  // const int num_tl = CCTK_NumTimeLevelsFromVarI (vindex);
  const int num_tl = CCTK_ActiveTimeLevelsVI(cctkGH, vindex);
  assert(num_tl >= 1);

  int const coord_group = CCTK_GroupIndex("grid::coordinates");

  // Loop over all maps
  const int m_min = 0;
  const int m_max = groupdata.grouptype == CCTK_GF ? Carpet::maps : 1;
  for (int m = m_min; m < m_max; ++m) {

    fstream file;
    OpenFile(cctkGH, m, vindex, alias, basefilename, dirs, is_new_file,
             truncate_file, file);

    // Find the output offset
    const ivect offset =
        groupdata.grouptype == CCTK_GF ? GetOutputOffset(cctkGH, m, dirs) : 0;

    const gh *const hh = arrdata.at(group).at(m).hh;
    const dh *const dd = arrdata.at(group).at(m).dd;

    // Traverse all components on this multigrid level, refinement
    // level, and map
    const int c_min = 0;
    const int c_max = groupdata.grouptype == CCTK_GF
                          ? hh->components(reflevel)
                          : groupdata.disttype != CCTK_DISTRIB_CONSTANT
                                ? CCTK_nProcs(cctkGH)
                                : 1;
    for (int c = c_min; c < c_max; ++c) {
      int const lc = hh->get_local_component(rl, c);
      int const proc = hh->processor(rl, c);
      if (dist::rank() == proc or dist::rank() == ioproc) {

        const ibbox &data_ext = dd->light_boxes.at(ml).at(rl).at(c).exterior;
        const ibbox ext = GetOutputBBox(cctkGH, group, rl, m, c, data_ext);

        CCTK_REAL coord_time;
        rvect coord_lower, coord_upper;
        GetCoordinates(cctkGH, m, groupdata, ext, coord_time, coord_lower,
                       coord_upper);

        // Apply offset
        ivect offset1;
        if (groupdata.grouptype == CCTK_GF) {
          const ibbox &baseext = hh->baseextents.at(ml).at(rl);
          offset1 = baseext.lower() + offset * ext.stride();
        } else {
          offset1 = offset * ext.stride();
        }
        for (int d = 0; d < outdim; ++d) {
          if (dirs[d] < 3) {
            offset1[dirs[d]] = ext.lower()[dirs[d]];
          }
        }

        const int tl_min = 0;
        const int tl_max = output_all_timelevels ? num_tl : 1;
        for (int tl = tl_min; tl < tl_max; ++tl) {

          mempool pool;

          const int n_min = one_file_per_group ? 0 : var;
          const int n_max =
              one_file_per_group ? CCTK_NumVarsInGroupI(group) : var + 1;
          vector<const gdata *> datas(n_max - n_min);
          for (size_t n = 0; n < datas.size(); ++n) {
            if (dist::rank() == proc) {
              const ggf *const ff = arrdata.at(group).at(m).data.at(n + n_min);
              datas.at(n) = ff->data_pointer(tl, rl, lc, ml);
            } else {
              datas.at(n) = NULL;
            }
          }

          vector<const gdata *> coords;
          if (use_grid_coordinates and groupdata.grouptype == CCTK_GF) {
            coords.resize(dim);
            for (int d = 0; d < dim; ++d) {
              if (dist::rank() == proc) {
                const ggf *const ff = arrdata.at(coord_group).at(m).data.at(d);
                coords.at(d) = ff->data_pointer(0, rl, lc, ml);
              } else {
                coords.at(d) = NULL;
              }
            }
          }

          vector<gdata *> tmpdatas(datas.size());
          vector<gdata *> tmpcoords(coords.size());

          if (proc != ioproc) {

            for (size_t n = 0; n < datas.size(); ++n) {
              const ggf *const ff = arrdata.at(group).at(m).data.at(n + n_min);
              tmpdatas.at(n) = ff->new_typed_data();
              size_t const memsize =
                  tmpdatas.at(n)->allocsize(data_ext, ioproc);
              void *const memptr = pool.alloc(memsize);
              tmpdatas.at(n)->allocate(data_ext, ioproc, memptr, memsize);
            } // for n
            for (size_t n = 0; n < coords.size(); ++n) {
              const ggf *const ff = arrdata.at(coord_group).at(m).data.at(n);
              tmpcoords.at(n) = ff->new_typed_data();
              size_t const memsize =
                  tmpcoords.at(n)->allocsize(data_ext, ioproc);
              void *const memptr = pool.alloc(memsize);
              tmpcoords.at(n)->allocate(data_ext, ioproc, memptr, memsize);
            } // for n

            for (comm_state state; not state.done(); state.step()) {
              for (size_t n = 0; n < datas.size(); ++n) {
                tmpdatas.at(n)->copy_from(state, datas.at(n), data_ext,
                                          data_ext, NULL, ioproc, proc);
              }
              for (size_t n = 0; n < coords.size(); ++n) {
                tmpcoords.at(n)->copy_from(state, coords.at(n), data_ext,
                                           data_ext, NULL, ioproc, proc);
              }
            }

          } else {

            for (size_t n = 0; n < datas.size(); ++n) {
              tmpdatas.at(n) = const_cast<gdata *>(datas.at(n));
            }
            for (size_t n = 0; n < coords.size(); ++n) {
              tmpcoords.at(n) = const_cast<gdata *>(coords.at(n));
            }
          }

          if (dist::rank() == ioproc) {
            WriteASCII(file, tmpdatas, ext, vindex, cctkGH->cctk_iteration,
                       offset1, dirs, rl, ml, m, c, tl, coord_time, coord_lower,
                       coord_upper, tmpcoords);
          }

          if (proc != ioproc) {
            for (size_t n = 0; n < tmpdatas.size(); ++n) {
              delete tmpdatas.at(n);
            }
            for (size_t n = 0; n < tmpcoords.size(); ++n) {
              delete tmpcoords.at(n);
            }
          }

          // Append EOL after every component
          if (dist::rank() == ioproc) {
            if (separate_components) {
              assert(file.good());
              file << eol;
            }
            assert(file.good());
          }

        } // for tl
      }
    } // for c

    // Append EOL after every complete set of components
    if (dist::rank() == ioproc) {
      if (separate_grids and not compact_format) {
        assert(file.good());
        file << eol;
      }
      assert(file.good());
    }

    CloseFile(cctkGH, file);

  } // for m
}

template <int outdim>
bool IOASCII<outdim>::DidOutput(const cGH *const cctkGH, const int vindex,
                                const string basefilename, bool &is_new_file,
                                bool &truncate_file) {
  DECLARE_CCTK_PARAMETERS;

  typedef std::map<string, vector<vector<vector<int> > > > filelist;
  static filelist created_files;

  filelist::iterator thisfile = created_files.find(basefilename);
  is_new_file = thisfile == created_files.end();
  truncate_file = is_new_file and IO_TruncateOutputFiles(cctkGH);

  if (is_new_file) {
    const int numelems = one_file_per_group ? CCTK_NumGroups() : CCTK_NumVars();
    vector<vector<vector<int> > > last_outputs; // [ml][rl][var]
    last_outputs.resize(mglevels);
    for (int ml = 0; ml < mglevels; ++ml) {
      last_outputs.at(ml).resize(maxreflevels);
      for (int rl = 0; rl < maxreflevels; ++rl) {
        last_outputs.at(ml).at(rl).resize(numelems, cctkGH->cctk_iteration - 1);
      }
    }
    // TODO: this makes a copy of last_outputs, which is expensive;
    // change this to use a reference instead
    thisfile = created_files.insert(
        thisfile, filelist::value_type(basefilename, last_outputs));
    assert(thisfile != created_files.end());
  }

  // Check if this variable has been output already during this
  // iteration
  const int elem =
      one_file_per_group ? CCTK_GroupIndexFromVarI(vindex) : vindex;
  int &last_output = thisfile->second.at(mglevel).at(reflevel).at(elem);
  if (last_output == cctkGH->cctk_iteration) {
    // has already been output during this iteration
    char *const fullname = CCTK_FullName(vindex);
    CCTK_VWarn(5, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Skipping output for variable '%s', because this variable "
               "has already been output during the current iteration -- "
               "probably via a trigger during the analysis stage",
               fullname);
    free(fullname);
    return true;
  }
  assert(last_output < cctkGH->cctk_iteration);
  last_output = cctkGH->cctk_iteration;

  return false;
}

CCTK_REAL io_files;
CCTK_REAL io_bytes_begin, io_bytes_end;

template <int outdim>
void IOASCII<outdim>::OpenFile(const cGH *const cctkGH, const int m,
                               const int vindex, const string alias,
                               const string basefilename,
                               const vect<int, outdim> &dirs,
                               const bool is_new_file, const bool truncate_file,
                               fstream &file) {
  DECLARE_CCTK_PARAMETERS;

  BeginTimingIO(cctkGH);
  io_files = 0;
  io_bytes_begin = 0;
  io_bytes_end = 0;

  if (dist::rank() == ioproc) {

    const int grouptype = CCTK_GroupTypeFromVarI(vindex);
    assert(grouptype >= 0);

    // Invent a file name
    ostringstream filenamebuf;
    filenamebuf << basefilename;
    if (maps > 1 and grouptype == CCTK_GF) {
      filenamebuf << "." << m;
    }
    filenamebuf << ".";
    if (new_filename_scheme) {
      for (int d = 0; d < outdim; ++d) {
        const char *const coords = "xyzd";
        filenamebuf << coords[dirs[d]];
      }
      filenamebuf << ".asc";
    } else {
      for (int d = 0; d < outdim; ++d) {
        assert(dirs[d] >= 0 and dirs[d] < 4);
        const char *const coords = "xyzd";
        filenamebuf << coords[dirs[d]];
      }
      const char *const suffixes = "plpv";
      filenamebuf << suffixes[outdim];
    }
    // we need a persistent temporary here
    const string filenamestr = filenamebuf.str();
    const char *const filename = filenamestr.c_str();

    // Open the file
    file.open(filename, ios::out | (truncate_file ? ios::trunc : ios::app));
    if (not file.good()) {
      char *const fullname = CCTK_FullName(vindex);
      CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                  "Could not open output file '%s' for variable '%s'", filename,
                  fullname);
      free(fullname);
    }
    io_files += 1;
    io_bytes_begin = file.tellg();

    // If this is the first time, then write a nice header
    if (is_new_file) {

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
        CCTK_ERROR("internal error");
      }

      file << "# " << outdim << "D ASCII output created by CarpetIOASCII"
           << eol;

      if (want_date) {
        char run_host[1000];
        Util_GetHostName(run_host, sizeof run_host);
        const char *run_user = getenv("USER");
        if (not run_user) {
          run_user = "";
        }
        char run_date[1000];
        Util_CurrentDate(sizeof run_date, run_date);
        char run_time[1000];
        Util_CurrentTime(sizeof run_time, run_time);
        file << "# created on " << run_host << " by " << run_user << " on "
             << run_date << " at " << run_time << eol;
        assert(file.good());
      }

      if (want_parfilename) {
        char parameter_filename[10000];
        CCTK_ParameterFilename(sizeof parameter_filename, parameter_filename);
        file << "# parameter filename: \"" << parameter_filename << "\"" << eol;
      }

      if (want_other) {
        if (CCTK_IsFunctionAliased("UniqueBuildID")) {
          const char *const build_id = (const char *)UniqueBuildID(cctkGH);
          file << "# Build ID: " << build_id << eol;
        }
        if (CCTK_IsFunctionAliased("UniqueSimulationID")) {
          const char *const job_id = (const char *)UniqueSimulationID(cctkGH);
          file << "# Simulation ID: " << job_id << eol;
        }
        if (CCTK_IsFunctionAliased("UniqueRunID")) {
          const char *const job_id = (const char *)UniqueRunID(cctkGH);
          file << "# Run ID: " << job_id << eol;
        }
      }

      file << "#" << eol;

      if (want_labels) {
        if (one_file_per_group) {
          char *const groupname = CCTK_GroupNameFromVarI(vindex);
          file << "# " << groupname;
          free(groupname);
        } else {
          const char *const varname = CCTK_VarName(vindex);
          file << "# " << varname;
        }
        for (int d = 0; d < outdim; ++d) {
          file << " "
               << "xyzd"[dirs[d]];
        }
        file << " (" << alias << ")" << eol;
        file << "#" << eol;
      }

    } // if is_new_file

    file << setprecision(out_precision);
    assert(file.good());

  } // if on the I/O processor
}

template <int outdim>
void IOASCII<outdim>::CloseFile(const cGH *const cctkGH, fstream &file) {
  DECLARE_CCTK_PARAMETERS;

  if (dist::rank() == ioproc) {
    io_bytes_end = file.tellg();
    file.close();
    assert(file.good());
  }

  assert(not file.is_open());

  CCTK_REAL const io_bytes = io_bytes_end - io_bytes_begin;
  EndTimingIO(cctkGH, io_files, io_bytes, false);
}

// Check whether this output direction has been requested
template <int outdim>
bool IOASCII<outdim>::DirectionIsRequested(const vect<int, outdim> &dirs) {
  DECLARE_CCTK_PARAMETERS;

  switch (outdim) {

  case 0:
    // Output is always requested (if switched on)
    return true;

  case 1:
    switch (dirs[0]) {
    case 0:
      return out1D_x;
    case 1:
      return out1D_y;
    case 2:
      return out1D_z;
    case 3:
      return out1D_d;
    }
    assert(0);

  case 2:
    if (dirs[0] == 0 and dirs[1] == 1)
      return out2D_xy;
    if (dirs[0] == 0 and dirs[1] == 2)
      return out2D_xz;
    if (dirs[0] == 1 and dirs[1] == 2)
      return out2D_yz;
    assert(0);

  case 3:
    // Output is always requested (if switched on)
    return true;
  }
  assert(0);
}

// Get the region that should be output, in terms of grid points;
// this is the offset perpendicular to the output hyperslab
template <int outdim>
ivect IOASCII<outdim>::GetOutputOffset(const cGH *const cctkGH, const int m,
                                       const vect<int, outdim> &dirs) {
  DECLARE_CCTK_PARAMETERS;

  // Default is zero
  ivect offset(0);

  switch (outdim) {

  case 0:
    // 0D output
    offset[0] =
        GetGridOffset(cctkGH, m, 1, "out0D_point_xi", /*"out_point_xi"*/ NULL,
                      "out0D_point_x", /*"out_point_x"*/ NULL,
                      /*out_point_x*/ 0.0);
    offset[1] =
        GetGridOffset(cctkGH, m, 2, "out0D_point_yi", /*"out_point_yi"*/ NULL,
                      "out0D_point_y", /*"out_point_y"*/ NULL,
                      /*out_point_y*/ 0.0);
    offset[2] =
        GetGridOffset(cctkGH, m, 3, "out0D_point_zi", /*"out_point_zi"*/ NULL,
                      "out0D_point_z", /*"out_point_z"*/ NULL,
                      /*out_point_z*/ 0.0);
    break;

  case 1:
    // 1D output
    switch (dirs[0]) {
    case 0:
      offset[1] = GetGridOffset(cctkGH, m, 2, "out1D_xline_yi", "out_xline_yi",
                                "out1D_xline_y", "out_xline_y", out_xline_y);
      offset[2] = GetGridOffset(cctkGH, m, 3, "out1D_xline_zi", "out_xline_zi",
                                "out1D_xline_z", "out_xline_z", out_xline_z);
      break;
    case 1:
      offset[0] = GetGridOffset(cctkGH, m, 1, "out1D_yline_xi", "out_yline_xi",
                                "out1D_yline_x", "out_yline_x", out_yline_x);
      offset[2] = GetGridOffset(cctkGH, m, 3, "out1D_yline_zi", "out_yline_zi",
                                "out1D_yline_z", "out_yline_z", out_yline_z);
      break;
    case 2:
      offset[0] = GetGridOffset(cctkGH, m, 1, "out1D_zline_xi", "out_zline_xi",
                                "out1D_zline_x", "out_zline_x", out_zline_x);
      offset[1] = GetGridOffset(cctkGH, m, 2, "out1D_zline_yi", "out_zline_yi",
                                "out1D_zline_y", "out_zline_y", out_zline_y);
      break;
    case 3:
      // the diagonal: we don't care about the offset
      break;
    default:
      assert(0);
    }
    break;

  case 2:
    // 2D output
    if (dirs[0] == 0 and dirs[1] == 1) {
      offset[2] =
          GetGridOffset(cctkGH, m, 3, "out2D_xyplane_zi", "out_xyplane_zi",
                        "out2D_xyplane_z", "out_xyplane_z", out_xyplane_z);
    } else if (dirs[0] == 0 and dirs[1] == 2) {
      offset[1] =
          GetGridOffset(cctkGH, m, 2, "out2D_xzplane_yi", "out_xzplane_yi",
                        "out2D_xzplane_y", "out_xzplane_y", out_xzplane_y);
    } else if (dirs[0] == 1 and dirs[1] == 2) {
      offset[0] =
          GetGridOffset(cctkGH, m, 1, "out2D_yzplane_xi", "out_yzplane_xi",
                        "out2D_yzplane_x", "out_yzplane_x", out_yzplane_x);
    } else {
      assert(0);
    }
    break;

  case 3:
    // 3D output: the offset does not matter
    break;

  default:
    assert(0);
  }

  return offset;
}

// Omit symmetry and ghost zones if requested
ibbox GetOutputBBox(const cGH *const cctkGH, const int group, const int rl,
                    const int m, const int c, const ibbox &ext) {
  DECLARE_CCTK_PARAMETERS;

  const int groupdim = CCTK_GroupDimI(group);
  assert(groupdim >= 0);
  const int grouptype = CCTK_GroupTypeI(group);
  assert(grouptype >= 0);

  // TODO: This is a bit ad hoc
  CCTK_INT symtable;
  if (grouptype == CCTK_GF and groupdim == cctkGH->cctk_dim) {
    symtable = SymmetryTableHandleForGrid(cctkGH);
    if (symtable < 0)
      CCTK_ERROR("internal error");
  } else {
    symtable = SymmetryTableHandleForGI(cctkGH, group);
    if (symtable < 0)
      CCTK_ERROR("internal error");
  }

  CCTK_INT symbnd[2 * dim];
  int const ierr =
      Util_TableGetIntArray(symtable, 2 * groupdim, symbnd, "symmetry_handle");
  if (ierr != 2 * groupdim)
    CCTK_ERROR("internal error");

  bool is_symbnd[2 * dim];
  for (int d = 0; d < 2 * groupdim; ++d) {
    is_symbnd[d] = symbnd[d] >= 0;
  }

  ivect lo = ext.lower();
  ivect hi = ext.upper();
  const ivect str = ext.stride();

  const b2vect obnds = vhh.at(m)->outer_boundaries(rl, c);
  const i2vect ghost_width = arrdata.at(group).at(m).dd->ghost_widths.at(rl);

  for (int d = 0; d < groupdim; ++d) {
    bool const output_lower_ghosts =
        obnds[0][d] ? (is_symbnd[2 * d]
                           ? output_symmetry_points
                           : (output_boundary_points and out3D_outer_ghosts))
                    : (output_ghost_points and out3D_ghosts);
    bool const output_upper_ghosts =
        obnds[1][d] ? (is_symbnd[2 * d + 1]
                           ? output_symmetry_points
                           : (output_boundary_points and out3D_outer_ghosts))
                    : (output_ghost_points and out3D_ghosts);

    if (not output_lower_ghosts) {
      lo[d] += ghost_width[0][d] * str[d];
    }
    if (not output_upper_ghosts) {
      hi[d] -= ghost_width[1][d] * str[d];
    }
  }

  return ibbox(lo, hi, str);
}

// Determine coordinates
void GetCoordinates(const cGH *const cctkGH, const int m,
                    const cGroup &groupdata, const ibbox &ext,
                    CCTK_REAL &coord_time, rvect &coord_lower,
                    rvect &coord_upper) {
  coord_time = cctkGH->cctk_time;

  rvect global_lower;
  rvect coord_delta;
  if (groupdata.grouptype == CCTK_GF) {
    rvect const cctk_origin_space = origin_space.at(m).at(mglevel);
    rvect const cctk_delta_space = delta_space.at(m) * rvect(mglevelfact);
    for (int d = 0; d < dim; ++d) {
      // lower boundary of Carpet's integer indexing
      global_lower[d] = cctk_origin_space[d];
      // grid spacing of Carpet's integer indexing
      coord_delta[d] = (cctk_delta_space[d] /
                        vhh.at(m)->baseextents.at(0).at(0).stride()[d]);
    }
  } else {
    for (int d = 0; d < dim; ++d) {
      global_lower[d] = 0.0;
      coord_delta[d] = 1.0;
    }
  }

  coord_lower = global_lower + coord_delta * rvect(ext.lower());
  coord_upper = global_lower + coord_delta * rvect(ext.upper());
}

int GetGridOffset(const cGH *const cctkGH, const int m, const int dir,
                  const char *const iparam, const char *const iglobal,
                  const char *const cparam, const char *const cglobal,
                  const CCTK_REAL cfallback) {
  // First choice: explicit coordinate
  const int ncparam = CCTK_ParameterQueryTimesSet(cparam, CCTK_THORNSTRING);
  assert(ncparam >= 0);
  if (ncparam > 0) {
    int ptype;
    const CCTK_REAL *const pcoord = ((const CCTK_REAL *)CCTK_ParameterGet(
        cparam, CCTK_THORNSTRING, &ptype));
    assert(pcoord);
    const CCTK_REAL coord = *pcoord;
    assert(ptype == PARAMETER_REAL);
    return CoordToOffset(cctkGH, m, dir, coord, 0);
  }

  // Second choice: explicit index
  const int niparam = CCTK_ParameterQueryTimesSet(iparam, CCTK_THORNSTRING);
  assert(niparam >= 0);
  if (niparam > 0) {
    int ptype;
    const int *const pindex =
        (const int *)CCTK_ParameterGet(iparam, CCTK_THORNSTRING, &ptype);
    assert(pindex);
    const int index = *pindex;
    assert(ptype == PARAMETER_INT);
    return index;
  }

  // Third choice: explicit global coordinate
  const char *iothorn = CCTK_ImplementationThorn("IO");
  assert(iothorn);
  if (cglobal) {
    const int ncglobal = CCTK_ParameterQueryTimesSet(cglobal, iothorn);
    assert(ncglobal >= 0);
    if (ncglobal > 0) {
      int ptype;
      const CCTK_REAL *const pcoord =
          (const CCTK_REAL *)CCTK_ParameterGet(cglobal, iothorn, &ptype);
      assert(pcoord);
      const CCTK_REAL coord = *pcoord;
      assert(ptype == PARAMETER_REAL);
      return CoordToOffset(cctkGH, m, dir, coord, 0);
    }
  }

  // Fourth choice: explicit global index
  if (iglobal) {
    const int niglobal = CCTK_ParameterQueryTimesSet(iglobal, iothorn);
    assert(niglobal >= 0);
    if (niglobal > 0) {
      int ptype;
      const int *const pindex =
          (const int *)CCTK_ParameterGet(iglobal, iothorn, &ptype);
      assert(pindex);
      const int index = *pindex;
      assert(ptype == PARAMETER_INT);
      return index;
    }
  }

  // Fifth choice: default coordinate
  return CoordToOffset(cctkGH, m, dir, cfallback, 0);
}

int CoordToOffset(const cGH *cctkGH, const int m, const int dir,
                  const CCTK_REAL coord, const int ifallback) {
  assert(m >= 0 and m < Carpet::maps and dir >= 1 and dir <= dim);

  assert(mglevel != -1 and reflevel != -1 and Carpet::map == -1);

  rvect const cctk_origin_space = origin_space.at(m).at(mglevel);
  rvect const cctk_delta_space = delta_space.at(m) * rvect(mglevelfact);
  ivect const cctk_levfac = spacereffacts.at(reflevel);
  ibbox const &coarseext = vhh.at(m)->baseextents.at(mglevel).at(0);
  ibbox const &baseext = vhh.at(m)->baseextents.at(mglevel).at(reflevel);
  ivect const cctk_levoff = baseext.lower() - coarseext.lower();
  ivect const cctk_levoffdenom = baseext.stride();

  const CCTK_REAL delta = cctk_delta_space[dir - 1] / cctk_levfac[dir - 1];
  const CCTK_REAL lower = cctk_origin_space[dir - 1] +
                          cctk_delta_space[dir - 1] / cctk_levfac[dir - 1] *
                              cctk_levoff[dir - 1] / cctk_levoffdenom[dir - 1];

  const CCTK_REAL rindex = (coord - lower) / delta;
  int cindex = (int)floor(rindex + 0.75);

  return cindex;
}

CCTK_REAL nicelooking(const CCTK_REAL val, const CCTK_REAL base) {
  return floor(val / base + 0.5) * base;
}

// Output
template <int outdim>
void WriteASCII(ostream &os, vector<gdata *> const &gfdatas,
                const bbox<int, dim> &gfext, const int vi, const int time,
                const vect<int, dim> &org, const vect<int, outdim> &dirs,
                const int rl, const int ml, const int m, const int c,
                const int tl, const CCTK_REAL coord_time,
                const vect<CCTK_REAL, dim> &coord_lower,
                const vect<CCTK_REAL, dim> &coord_upper,
                vector<gdata *> const &gfcoords) {
  DECLARE_CCTK_PARAMETERS;

  assert(outdim <= dim);

  const int vartype = CCTK_VarTypeI(vi);
  const int grouptype = CCTK_GroupTypeFromVarI(vi);
  const int groupdim = CCTK_GroupDimFromVarI(vi);

  if (CCTK_EQUALS(out_fileinfo, "axis labels") or
      CCTK_EQUALS(out_fileinfo, "all")) {

    assert(os.good());

    if (not compact_format or grouptype == CCTK_GF) {
      // Don't output a comment with iteration number and time,
      // because these are always output with the real data anyway
      os << "# iteration " << time << "   time " << coord_time << eol;
      os << "# time level " << tl << eol;
      os << "# refinement level " << rl << "   multigrid level " << ml
         << "   map " << m << "   component " << c << eol;
    }
    static vector<bool> did_output_format;
    if (did_output_format.empty()) {
      did_output_format.resize(CCTK_NumVars());
    }
    if (not compact_format or not did_output_format.AT(vi)) {
      did_output_format.AT(vi) = true;
      os << "# column format: 1:it";
      int col = 2;
      if (not compact_format or output_all_timelevels) {
        os << "\t" << col++ << ":tl";
      }
      if (not compact_format) {
        os << "\t" << col++ << ":rl";
        os << " " << col++ << ":c";
        os << " " << col++ << ":ml";
      }
      assert(dim >= 0 and dim <= 3);
      const char *const coords = "xyz";
      if (not compact_format) {
        for (int d = 0; d < dim; ++d) {
          os << (d == 0 ? "\t" : " ") << col++ << ":i" << coords[d];
        }
      } else {
        for (int d = 0; d < min(dim, groupdim); ++d) {
          os << (d == 0 ? "\t" : " ") << col++ << ":i" << coords[d];
        }
      }
      os << "\t" << col++ << ":time";
      if (not compact_format or grouptype == CCTK_GF) {
        for (int d = 0; d < dim; ++d) {
          os << (d == 0 ? "\t" : " ") << col++ << ":" << coords[d];
        }
      }
      os << "\t" << col << ":data" << eol;
      if (one_file_per_group) {
        os << "# data columns:";
        int const gindex = CCTK_GroupIndexFromVarI(vi);
        int const firstvar = CCTK_FirstVarIndexI(gindex);
        int const numvars = CCTK_NumVarsInGroupI(gindex);
        for (int n = firstvar; n < firstvar + numvars; ++n) {
          os << " " << col << ":" << CCTK_VarName(n);
          col += CarpetSimpleMPIDatatypeLength(vartype);
        }
        os << eol;
      }
    }

  } // if out_fileinfo

  // boolean that says if we are doing 1D-diagonal output
  // This is not beautiful, but works for the moment
  bool const diagonal_output = outdim == 1 and dirs[0] == 3;

  if (not diagonal_output) { // not outputting the diagonal

    const vect<int, outdim> lo = gfext.lower()[dirs];
    const vect<int, outdim> up = gfext.upper()[dirs];
    const vect<int, outdim> str = gfext.stride()[dirs];
    const bbox<int, outdim> ext(lo, up, str);

    // check whether the output origin is contained in the extent of
    // the data that should be output
    ivect org1(org);
    for (int d = 0; d < outdim; ++d)
      org1[dirs[d]] = ext.lower()[d];
    if (gfext.contains(org1)) {

      typename bbox<int, outdim>::iterator it = ext.begin();
      do {

        ivect index(org);
        for (int d = 0; d < outdim; ++d)
          index[dirs[d]] = (*it)[d];
        os << time;
        if (not compact_format or output_all_timelevels) {
          os << "\t" << tl;
        }
        if (not compact_format) {
          // Don't output the grid structure in compact format (it
          // is still output in the comments above every component)
          os << "\t" << rl << " " << c << " " << ml;
        }
        if (not compact_format) {
          for (int d = 0; d < dim; ++d) {
            os << (d == 0 ? "\t" : " ") << index[d];
          }
        } else {
          // In the compact format, output one column for each of
          // the group's dimension, don't always output three
          // columns
          for (int d = 0; d < min(dim, groupdim); ++d) {
            os << (d == 0 ? "\t" : " ") << index[d];
          }
        }
        os << "\t" << coord_time;
        if (not compact_format or grouptype == CCTK_GF) {
          for (int d = 0; d < dim; ++d) {
            os << (d == 0 ? "\t" : " ");
            assert(gfext.upper()[d] - gfext.lower()[d] >= 0);
            if (gfcoords.empty()) {
              // Calculate coordinates
              if (gfext.upper()[d] - gfext.lower()[d] == 0) {
                os << coord_lower[d];
              } else {
                CCTK_REAL const dx = ((coord_upper[d] - coord_lower[d]) /
                                      (gfext.upper()[d] - gfext.lower()[d]));
                os << (nicelooking(coord_lower[d] +
                                       (index[d] - gfext.lower()[d]) * dx,
                                   dx * 1.0e-8));
              }
            } else {
              // Use coordinate grid functions
              const gdata *gfcoord = gfcoords.at(d);
              os << (*(const data<CCTK_REAL> *)gfcoord)[index];
            }
          }
        }
        for (size_t n = 0; n < gfdatas.size(); ++n) {
          const gdata *gfdata = gfdatas.at(n);
          os << (n == 0 ? "\t" : " ");
          switch (specific_cactus_type(vartype)) {
#define CARPET_NO_COMPLEX
#define TYPECASE(N, T)                                                         \
  case N:                                                                      \
    os << (*(const data<T> *)gfdata)[index];                                   \
    break;
#include "typecase.hh"
#undef TYPECASE
#undef CARPET_NO_COMPLEX
#define CARPET_COMPLEX
#define TYPECASE(N, T)                                                         \
  case N:                                                                      \
    os << real((*(const data<T> *)gfdata)[index]) << " "                       \
       << imag((*(const data<T> *)gfdata)[index]);                             \
    break;
#include "typecase.hh"
#undef TYPECASE
#undef CARPET_COMPLEX
          default:
            UnsupportedVarType(vi);
          }
        } // for n
        os << eol;

        ++it;

        for (int d = 0; d < outdim; ++d) {
          if ((*it)[d] != (*ext.end())[d])
            break;
          if (not compact_format or ext.shape()[d] > ext.stride()[d]) {
            // In the compact format, don't separate outputs that
            // consist of a single lines only
            os << eol;
          }
        }

      } while (it != ext.end());

    } else {

      if (not compact_format) {
        os << "#" << eol;
      }

    } // if not ext contains org

    assert(os.good());

  } else { // taking care of the diagonal

    const ivect lo = gfext.lower();
    const ivect up = gfext.upper();
    const ivect str = gfext.stride();
    const ibbox ext(lo, up, str);

    gh const &hh = *vhh.at(m);
    ibbox const &base = hh.baseextents.at(mglevel).at(reflevel);

    assert(base.stride()[0] == base.stride()[1] and
           base.stride()[0] == base.stride()[2]);

    // output the data on the diagonal
    for (int i = maxval(base.lower()); i <= minval(base.upper());
         i += base.stride()[0]) {

      ivect const pos = ivect(i, i, i);

      // check if the point in question is in our gf's extent
      if (gfext.contains(pos)) {
        os << time;
        if (not compact_format) {
          // Don't output the grid structure in compact format (it
          // is still output in the comments above every component)
          os << "\t" << tl << " " << rl << " " << c << " " << ml;
        }
        for (int d = 0; d < dim; ++d) {
          os << (d == 0 ? "\t" : " ") << pos[d];
        }
        os << "\t" << coord_time;
        for (int d = 0; d < dim; ++d) {
          os << (d == 0 ? "\t" : " ");
          assert(gfext.upper()[d] - gfext.lower()[d] >= 0);
          if (gfext.upper()[d] - gfext.lower()[d] == 0) {
            os << coord_lower[d];
          } else {
            CCTK_REAL const dx = ((coord_upper[d] - coord_lower[d]) /
                                  (gfext.upper()[d] - gfext.lower()[d]));
            os << (nicelooking(coord_lower[d] +
                                   (pos[d] - gfext.lower()[d]) * dx,
                               dx * 1.0e-8));
          }
        }
        for (size_t n = 0; n < gfdatas.size(); ++n) {
          const gdata *gfdata = gfdatas.at(n);
          os << (n == 0 ? "\t" : " ");
          switch (specific_cactus_type(vartype)) {
#define TYPECASE(N, T)                                                         \
  case N:                                                                      \
    os << (*(const data<T> *)gfdata)[pos];                                     \
    break;
#include "typecase.hh"
#undef TYPECASE
          default:
            UnsupportedVarType(vi);
          }
        } // for n
        os << eol;

      } else {

        if (not compact_format) {
          os << "#" << eol;
        }

      } // if not ext contains org

    } // end for loop

    if (not compact_format or maxval(base.lower()) > minval(base.upper())) {
      // In the compact format, don't separate outputs that consist
      // of a single lines only
      os << eol;
    }

    assert(os.good());

  } // if diagonal_output
}

// Explicit instantiation for all output dimensions
template struct IOASCII<0>;
template struct IOASCII<1>;
template struct IOASCII<2>;
template struct IOASCII<3>;

template void WriteASCII(ostream &os, vector<gdata *> const &gfdatas,
                         const bbox<int, dim> &gfext, const int vi,
                         const int time, const vect<int, dim> &org,
                         const vect<int, 0> &dirs, const int rl, const int ml,
                         const int m, const int c, const int tl,
                         const CCTK_REAL coord_time,
                         const vect<CCTK_REAL, dim> &coord_lower,
                         const vect<CCTK_REAL, dim> &coord_upper,
                         vector<gdata *> const &gfcoords);

template void WriteASCII(ostream &os, vector<gdata *> const &gfdatas,
                         const bbox<int, dim> &gfext, const int vi,
                         const int time, const vect<int, dim> &org,
                         const vect<int, 1> &dirs, const int rl, const int ml,
                         const int m, const int c, const int tl,
                         const CCTK_REAL coord_time,
                         const vect<CCTK_REAL, dim> &coord_lower,
                         const vect<CCTK_REAL, dim> &coord_upper,
                         vector<gdata *> const &gfcoords);

template void WriteASCII(ostream &os, vector<gdata *> const &gfdatas,
                         const bbox<int, dim> &gfext, const int vi,
                         const int time, const vect<int, dim> &org,
                         const vect<int, 2> &dirs, const int rl, const int ml,
                         const int m, const int c, const int tl,
                         const CCTK_REAL coord_time,
                         const vect<CCTK_REAL, dim> &coord_lower,
                         const vect<CCTK_REAL, dim> &coord_upper,
                         vector<gdata *> const &gfcoords);

template void WriteASCII(ostream &os, vector<gdata *> const &gfdatas,
                         const bbox<int, dim> &gfext, const int vi,
                         const int time, const vect<int, dim> &org,
                         const vect<int, 3> &dirs, const int rl, const int ml,
                         const int m, const int c, const int tl,
                         const CCTK_REAL coord_time,
                         const vect<CCTK_REAL, dim> &coord_lower,
                         const vect<CCTK_REAL, dim> &coord_upper,
                         vector<gdata *> const &gfcoords);

} // namespace CarpetIOASCII
