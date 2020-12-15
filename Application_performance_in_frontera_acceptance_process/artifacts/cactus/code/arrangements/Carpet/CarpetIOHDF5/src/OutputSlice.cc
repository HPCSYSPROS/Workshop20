#include <cassert>
#include <cctype>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <cstring>

#include <map>
#include <string>

#include <cctk.h>
#include <util_Table.h>

#include "CactusBase/IOUtil/src/ioGH.h"

#include <Timer.hh>

#include "typeprops.hh"

#include "CarpetIOHDF5.hh"

// That's a hack
namespace Carpet {
void UnsupportedVarType(const int vindex);
}

#define GetParameter(parameter)                                                \
  outdim == 0 ? out0D_##parameter : outdim == 1                                \
                                        ? out1D_##parameter                    \
                                        : outdim == 2 ? out2D_##parameter      \
                                                      : out3D_##parameter

namespace CarpetIOHDF5 {

using namespace std;
using namespace Carpet;

// routines which are independent of the output dimension
static ibset GetOutputBBoxes(const cGH *cctkGH, int group, int rl, int m, int c,
                             const ibbox &ext, const ibset &allactive);

static void GetCoordinates(const cGH *cctkGH, int m, const cGroup &groupdata,
                           const ibset &exts, CCTK_REAL &coord_time,
                           vector<rvect> &coord_lower,
                           vector<rvect> &coord_upper);

static int GetGridOffset(const cGH *cctkGH, int m, int dir, const char *itempl,
                         const char *iglobal, const char *ctempl,
                         const char *cglobal, CCTK_REAL cfallback);
static int CoordToOffset(const cGH *cctkGH, int m, int dir, CCTK_REAL coord,
                         int ifallback);

// IO processor
template <int outdim> int IOHDF5<outdim>::ioproc;
template <int outdim> int IOHDF5<outdim>::nioprocs;
template <int outdim> int IOHDF5<outdim>::ioproc_every;

// Global configuration parameters
bool stop_on_parse_errors = false;

// Definition of static members
template <int outdim> char *IOHDF5<outdim>::my_out_slice_dir;
template <int outdim> char *IOHDF5<outdim>::my_out_slice_vars;
template <int outdim> vector<ioRequest *> IOHDF5<outdim>::slice_requests;

template <int outdim> int IOHDF5<outdim>::Startup() {
  ostringstream msg;
  msg << "AMR " << outdim << "D HDF5 I/O provided by CarpetIOHDF5";
  CCTK_RegisterBanner(msg.str().c_str());

  ostringstream extension_name;
  extension_name << "CarpetIOHDF5_" << outdim << "D";
  const int GHExtension =
      CCTK_RegisterGHExtension(extension_name.str().c_str());
  CCTK_RegisterGHExtensionSetupGH(GHExtension, SetupGH);

  ostringstream method_name;
  method_name << "IOHDF5_" << outdim << "D";
  const int IOMethod = CCTK_RegisterIOMethod(method_name.str().c_str());
  CCTK_RegisterIOMethodOutputGH(IOMethod, OutputGH);
  CCTK_RegisterIOMethodOutputVarAs(IOMethod, OutputVarAs);
  CCTK_RegisterIOMethodTimeToOutput(IOMethod, TimeToOutput);
  CCTK_RegisterIOMethodTriggerOutput(IOMethod, TriggerOutput);

  return 0;
}

template <int outdim>
void *IOHDF5<outdim>::SetupGH(tFleshConfig *const fc, const int convLevel,
                              cGH *const cctkGH) {
  DECLARE_CCTK_PARAMETERS;
  const void *dummy;

  dummy = &fc;
  dummy = &convLevel;
  dummy = &cctkGH;
  dummy = &dummy;

  if (not CCTK_Equals(verbose, "none")) {
    CCTK_VInfo(CCTK_THORNSTRING,
               "I/O Method 'IOHDF5_%dD' registered: "
               "%dD AMR output of grid variables to HDF5 files",
               outdim, outdim);
  }

  const int numvars = CCTK_NumVars();
  slice_requests.resize(numvars);

  // initial I/O parameter check
  my_out_slice_dir = 0;
  my_out_slice_vars = strdup("");
  stop_on_parse_errors = strict_io_parameter_check != 0;
  CheckSteerableParameters(cctkGH);
  stop_on_parse_errors = false;

  // We register only once, ergo we get only one handle.  We store
  // that statically, so there is no need to pass anything to
  // Cactus.
  return NULL;
}

template <int outdim>
void IOHDF5<outdim>::CheckSteerableParameters(const cGH *const cctkGH) {
  DECLARE_CCTK_PARAMETERS;

  // re-parse the 'IOHDF5::out%dD_dir' parameter if it has changed
  const char *the_out_dir = GetParameter(dir);
  if (CCTK_EQUALS(the_out_dir, "")) {
    the_out_dir = io_out_dir;
  }

  if (not my_out_slice_dir or strcmp(the_out_dir, my_out_slice_dir)) {
    free(my_out_slice_dir);
    my_out_slice_dir = strdup(the_out_dir);

    // create the output directory
    const int result = IOUtil_CreateDirectory(cctkGH, my_out_slice_dir, 0, 0);
    if (result < 0) {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Problem creating %dD-output directory '%s'", outdim,
                 my_out_slice_dir);
    } else if (result > 0 and CCTK_Equals(verbose, "full")) {
      CCTK_VInfo(CCTK_THORNSTRING, "%dD-output directory '%s' already exists",
                 outdim, my_out_slice_dir);
    }
  }

  // re-parse the 'IOHDF5::out%d_vars' parameter if it has changed
  const char *const out_slice_vars = GetParameter(vars);
  if (strcmp(out_slice_vars, my_out_slice_vars)) {
    ostringstream parameter_name;
    parameter_name << "IOHDF5::out" << outdim << "D_vars";
#ifdef IOUTIL_PARSER_HAS_OUT_DT
    IOUtil_ParseVarsForOutput(
        cctkGH, CCTK_THORNSTRING, parameter_name.str().c_str(),
        stop_on_parse_errors, out_slice_vars, -1, -1.0, &slice_requests[0]);
#else
    IOUtil_ParseVarsForOutput(
        cctkGH, CCTK_THORNSTRING, parameter_name.str().c_str(),
        stop_on_parse_errors, out_slice_vars, -1, &slice_requests[0]);
#endif

    // notify the user about the new setting
    if (not CCTK_Equals(verbose, "none")) {
      int count = 0;
      ostringstream msg;
      msg << "Periodic " << outdim << "D AMR output requested for:";
      for (int vi = 0; vi < CCTK_NumVars(); ++vi) {
        if (slice_requests.at(vi)) {
          ++count;
          char *const fullname = CCTK_FullName(vi);
          msg << "\n"
              << "   " << fullname;
          free(fullname);
        }
      }
      if (count > 0) {
        CCTK_INFO(msg.str().c_str());
      }
    }

    // save the last setting of 'IOHDF5::out%d_vars' parameter
    free(my_out_slice_vars);
    my_out_slice_vars = strdup(out_slice_vars);
  }

  // copy ioprocs and ioproc ot
  if (outdim == 3) { // only 3D output splits files
    const ioGH *IO = (ioGH *)CCTK_GHExtension(cctkGH, "IO");
    assert(IO);
    nioprocs = IO->nioprocs;
    ioproc = IO->ioproc;
    ioproc_every = IO->ioproc_every;
    // cout << "nioprocs: " << nioprocs << " ioproc: " << ioproc << endl;
  } else {
    nioprocs = 1;
    ioproc = 0;
    ioproc_every = dist::size();
  }
}

template <int outdim> int IOHDF5<outdim>::OutputGH(const cGH *const cctkGH) {
  ostringstream timer_name;
  timer_name << "OutputGH<" << outdim << ">";

  Timers::Timer timer(timer_name.str());

  timer.start();

  CheckSteerableParameters(cctkGH);
  if (strcmp(my_out_slice_vars, "")) {
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
int IOHDF5<outdim>::TimeToOutput(const cGH *const cctkGH, const int vindex) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  assert(vindex >= 0 and vindex < CCTK_NumVars());

  if (CCTK_GroupTypeFromVarI(vindex) != CCTK_GF and not do_global_mode) {
    return 0;
  }

  // check if output for this variable was requested
  if (not slice_requests.at(vindex)) {
    return 0;
  }

  // check whether this refinement level should be output
  if (not(slice_requests.at(vindex)->refinement_levels & (1 << reflevel))) {
    return 0;
  }

  // check if output for this variable was requested individually by
  // a "<varname>{ out_every = <number> }" option string
  // this will overwrite the output criterion setting
  const char *myoutcriterion = GetParameter(criterion);
  if (CCTK_EQUALS(myoutcriterion, "default")) {
    myoutcriterion = io_out_criterion;
  }
  if (slice_requests.at(vindex)->out_every >= 0) {
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
      myoutevery = io_out_every;
    }
    if (myoutevery > 0) {
      if (cctk_iteration == this_iteration_slice[outdim]) {
        // we already decided to output this iteration
        output_this_iteration = true;
      } else if (cctk_iteration >=
                 last_output_iteration_slice[outdim] + myoutevery) {
        // it is time for the next output
        output_this_iteration = true;
        last_output_iteration_slice[outdim] = cctk_iteration;
        this_iteration_slice[outdim] = cctk_iteration;
      }
    }
  } else if (CCTK_EQUALS(myoutcriterion, "divisor")) {
    int myoutevery = GetParameter(every);
    if (myoutevery == -2) {
      myoutevery = io_out_every;
    }
    if (slice_requests.at(vindex)->out_every >= 0) {
      myoutevery = slice_requests.at(vindex)->out_every;
    }
    if (myoutevery > 0 and (cctk_iteration % myoutevery) == 0) {
      // we already decided to output this iteration
      output_this_iteration = true;
    }
  } else if (CCTK_EQUALS(myoutcriterion, "time")) {
    CCTK_REAL myoutdt = GetParameter(dt);
    if (myoutdt == -2) {
      myoutdt = io_out_dt;
    }
    if (myoutdt == 0 or cctk_iteration == this_iteration_slice[outdim]) {
      output_this_iteration = true;
    } else if (myoutdt > 0) {
      int do_output =
          cctk_time / cctk_delta_time >=
          (last_output_time_slice[outdim] + myoutdt) / cctk_delta_time -
              1.0e-12;
      MPI_Bcast(&do_output, 1, MPI_INT, 0, dist::comm());
      if (do_output) {
        // it is time for the next output
        output_this_iteration = true;
        last_output_time_slice[outdim] = cctk_time;
        this_iteration_slice[outdim] = cctk_iteration;
      }
    }
  } // select output criterion

  return output_this_iteration ? 1 : 0;
}

template <int outdim>
int IOHDF5<outdim>::TriggerOutput(const cGH *const cctkGH, const int vindex) {
  DECLARE_CCTK_PARAMETERS;

  assert(vindex >= 0 and vindex < CCTK_NumVars());

  char *const fullname = CCTK_FullName(vindex);

  int retval;
  if (one_file_per_proc) {
    char path[500];
    CCTK_ParameterFilename(500, path);
    char *value = strrchr (path, '/');
    if (value == NULL) {
      value = path;
    } else {
      value++;
    }
    char *dot = strrchr(value, '.');
    if (dot != NULL && strcmp(dot, ".par") == 0) {
      *dot = '\0';
    }
    retval = OutputVarAs(cctkGH, fullname, value);
  } else if (one_file_per_group) {
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
    CCTK_VWarn(2, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Option string '%s' will be ignored for HDF5 output of "
               "variable '%s'",
               optstring, fullname);
    free(fullname);
  }

  *static_cast<int *>(arg) = vindex;
}

template <int outdim>
int IOHDF5<outdim>::OutputVarAs(const cGH *const cctkGH,
                                const char *const varname,
                                const char *const alias) {
  DECLARE_CCTK_PARAMETERS;

  int vindex = -1;

  if (CCTK_TraverseString(varname, GetVarIndex, &vindex, CCTK_VAR) < 0) {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
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
    CCTK_WARN(1, "OutputVarAs must be called in level mode");
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
      CCTK_VWarn(
          1, __LINE__, __FILE__, CCTK_THORNSTRING,
          "Cannot produce %dD slice HDF5 output file '%s' for variable '%s' "
          "because it has only %d dimensions",
          outdim, alias, varname, groupdata.dim);
      return -1;
    }

    // Check for storage
    if (not CCTK_QueryGroupStorageI(cctkGH, group)) {
      // This may be okay if storage is e.g. scheduled only in the
      // analysis bin
      CCTK_VWarn(4, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Cannot output variable '%s' because it has no storage",
                 varname);
      return 0;
    }

    ostringstream basefilenamebuf;
    basefilenamebuf << my_out_slice_dir << "/" << alias;
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
void IOHDF5<outdim>::OutputDirection(const cGH *const cctkGH, const int vindex,
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

  const int num_tl = CCTK_MaxActiveTimeLevelsVI(cctkGH, vindex);
  assert(num_tl >= 1);

  const int numvars = one_file_per_group ? CCTK_NumVarsInGroupI(group) : 1;

  // Loop over all maps
  const int m_min = 0;
  const int m_max = groupdata.grouptype == CCTK_GF ? Carpet::maps : 1;
  for (int m = m_min; m < m_max; ++m) {

    hid_t file = -1, index_file = -1;
    int error_count = 0;
    error_count += OpenFile(cctkGH, m, vindex, numvars, alias, basefilename,
                            dirs, is_new_file, truncate_file, file, index_file);

    // Find the output offset
    const ivect offset =
        groupdata.grouptype == CCTK_GF ? GetOutputOffset(cctkGH, m, dirs) : 0;

    const gh *const hh = arrdata.at(group).at(m).hh;
    const dh *const dd = arrdata.at(group).at(m).dd;

    // re-compute the active (non-buffered) region
    ibset allactive;
    GetAllActive(dd, hh, ml, rl, allactive);

    // Traverse all components on this multigrid level, refinement
    // level, and map
    const int c_min = 0;
    const int c_max = groupdata.grouptype == CCTK_GF
                          ? vhh.at(m)->components(reflevel)
                          : groupdata.disttype != CCTK_DISTRIB_CONSTANT
                                ? CCTK_nProcs(cctkGH)
                                : 1;
    for (int c = c_min, c_base = c_min; c < c_max; ++c) {
      int const lc = hh->get_local_component(rl, c);
      int const proc = hh->processor(rl, c);
      const ibbox &data_ext = dd->light_boxes.at(ml).at(rl).at(c).exterior;
      const ibset exts =
          GetOutputBBoxes(cctkGH, group, rl, m, c, data_ext, allactive);
      // we have to take part in this output if we either own data to be
      // output or are the ioproc in the same group of processors as the
      // dataholder
      if (dist::rank() == proc or dist::rank() == IOProcForProc(proc)) {

        CCTK_REAL coord_time;
        vector<rvect> coord_lower, coord_upper;
        GetCoordinates(cctkGH, m, groupdata, exts, coord_time, coord_lower,
                       coord_upper);

        // Apply offset
        vector<ivect> offsets1;
        offsets1.reserve(exts.setsize());
        for (ibset::const_iterator ext = exts.begin(); ext != exts.end();
             ++ext) {
          ivect offset1;
          if (groupdata.grouptype == CCTK_GF) {
            const ibbox &baseext = hh->baseextents.at(ml).at(rl);
            offset1 = baseext.lower() + offset * ext->stride();
          } else {
            offset1 = offset * ext->stride();
          }
          for (int d = 0; d < outdim; ++d) {
            if (dirs[d] < 3) {
              offset1[dirs[d]] = ext->lower()[dirs[d]];
            }
          }
          offsets1.push_back(offset1);
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

          vector<gdata *> tmpdatas(datas.size());

          // tranfer data through the interconnect to ioproc if necessary
          if (proc != ioproc) {

            for (size_t n = 0; n < datas.size(); ++n) {
              const ggf *const ff = arrdata.at(group).at(m).data.at(n + n_min);
              tmpdatas.at(n) = ff->new_typed_data();
              size_t const memsize =
                  tmpdatas.at(n)->allocsize(data_ext, ioproc);
              void *const memptr = pool.alloc(memsize);
              tmpdatas.at(n)->allocate(data_ext, ioproc, memptr, memsize);
            } // for n

            for (comm_state state; not state.done(); state.step()) {
              for (size_t n = 0; n < datas.size(); ++n) {
                tmpdatas.at(n)->copy_from(state, datas.at(n), data_ext,
                                          data_ext, NULL, ioproc, proc);
              }
            }

          } else {

            for (size_t n = 0; n < datas.size(); ++n) {
              tmpdatas.at(n) = const_cast<gdata *>(datas.at(n));
            }
          }

          if (dist::rank() == IOProcForProc(proc)) {
            int c_offset = 0;
            for (ibset::const_iterator ext = exts.begin(); ext != exts.end();
                 ++ext, ++c_offset) {
              error_count += WriteHDF5(
                  cctkGH, file, index_file, tmpdatas, *ext, vindex,
                  offsets1[c_offset], dirs, rl, ml, m, c, c_base + c_offset, tl,
                  coord_time, coord_lower[c_offset], coord_upper[c_offset]);
            }
          }

          if (proc != ioproc) {
            for (size_t n = 0; n < tmpdatas.size(); ++n) {
              delete tmpdatas.at(n);
            }
          }

        } // for tl
      }

      c_base += exts.setsize();
    } // for c

    error_count += CloseFile(cctkGH, file, index_file);
    if (error_count > 0 and abort_on_io_errors) {
      CCTK_ERROR("Aborting simulation due to previous I/O errors");
    }

  } // for m
}

template <int outdim>
bool IOHDF5<outdim>::DidOutput(const cGH *const cctkGH, const int vindex,
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
CCTK_REAL io_bytes;

template <int outdim>
int IOHDF5<outdim>::OpenFile(const cGH *const cctkGH, const int m,
                             const int vindex, const int numvars,
                             const string alias, const string basefilename,
                             const vect<int, outdim> &dirs,
                             const bool is_new_file, const bool truncate_file,
                             hid_t &file, hid_t &index_file) {
  DECLARE_CCTK_PARAMETERS;

  int error_count = 0;

  BeginTimingIO(cctkGH);
  io_files = 0;
  io_bytes = 0;

  if (dist::rank() == ioproc) {

    const int grouptype = CCTK_GroupTypeFromVarI(vindex);
    assert(grouptype >= 0);

    // Invent a file name
    ostringstream filenamebuf;
    filenamebuf << basefilename;
    if (maps > 1 and grouptype == CCTK_GF) {
      filenamebuf << "." << m;
    }
    // historically 3d output files do not carry a label so the files created
    // here do not conflict with the old-style output files
    filenamebuf << ".";
    for (int d = 0; d < outdim; ++d) {
      const char *const coords = "xyzd";
      filenamebuf << coords[dirs[d]];
    }
    if (nioprocs > 1) {
      filenamebuf << ".file_" << dist::rank();
    }
    string index_filename(filenamebuf.str() + ".idx" + out_extension);
    filenamebuf << out_extension;

    // we need a persistent temporary here
    const string filenamestr = filenamebuf.str();
    const char *const filename = filenamestr.c_str();

    // Open the file
    bool file_exists = false;
    if (not truncate_file) {
      H5E_BEGIN_TRY { file_exists = H5Fis_hdf5(filename) > 0; }
      H5E_END_TRY;
    }

    if (truncate_file or not file_exists) {
      hid_t fapl_id;
      HDF5_ERROR(fapl_id = H5Pcreate(H5P_FILE_ACCESS));
      HDF5_ERROR(H5Pset_fclose_degree(fapl_id, H5F_CLOSE_STRONG));
      HDF5_ERROR(file =
                     H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id));
      if (output_index) {
        HDF5_ERROR(index_file = H5Fcreate(index_filename.c_str(), H5F_ACC_TRUNC,
                                          H5P_DEFAULT, fapl_id));
      }
      HDF5_ERROR(H5Pclose(fapl_id));
      // write metadata information
      error_count +=
          WriteMetadata(cctkGH, nioprocs, vindex, numvars, false, file);

      if (output_index) {
        error_count +=
            WriteMetadata(cctkGH, nioprocs, vindex, numvars, false, index_file);
      }
    } else {
      hid_t fapl_id;
      HDF5_ERROR(fapl_id = H5Pcreate(H5P_FILE_ACCESS));
      HDF5_ERROR(H5Pset_fclose_degree(fapl_id, H5F_CLOSE_STRONG));
      HDF5_ERROR(file = H5Fopen(filename, H5F_ACC_RDWR, fapl_id));
      if (output_index)
        HDF5_ERROR(index_file =
                       H5Fopen(index_filename.c_str(), H5F_ACC_RDWR, fapl_id));
      HDF5_ERROR(H5Pclose(fapl_id));
    }
    io_files += 1;

  } // if on the I/O processor

  return error_count;
}

template <int outdim>
int IOHDF5<outdim>::CloseFile(const cGH *const cctkGH, hid_t &file,
                              hid_t &index_file) {
  DECLARE_CCTK_PARAMETERS;

  int error_count = 0;

  if (dist::rank() == ioproc) {
    if (file >= 0) {
      HDF5_ERROR(H5Fclose(file));
    }
    if (output_index and index_file >= 0) {
      HDF5_ERROR(H5Fclose(index_file));
    }
    HDF5_ERROR(H5garbage_collect());
  }
  if (nioprocs > 1) {
    CCTK_REAL local[2], global[2];
    local[0] = io_files;
    local[1] = io_bytes;
    MPI_Allreduce(local, global, 2, dist::mpi_datatype(local[0]), MPI_SUM,
                  dist::comm());
    io_files = global[0];
    io_bytes = global[1];
  }

  EndTimingIO(cctkGH, io_files, io_bytes, true);

  return error_count;
}

// Check whether this output direction has been requested
template <int outdim>
bool IOHDF5<outdim>::DirectionIsRequested(const vect<int, outdim> &dirs) {
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
    default:
      assert(0);
    }

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

  default:
    assert(0);
    // Prevent compiler warning about missing return statement
    return false;
  }
}

// Get the region that should be output, in terms of grid points;
// this is the offset perpendicular to the output hyperslab
template <int outdim>
ivect IOHDF5<outdim>::GetOutputOffset(const cGH *const cctkGH, const int m,
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

// Omit symmetry, ghost and buffer zones if requested
ibset GetOutputBBoxes(const cGH *const cctkGH, const int group, const int rl,
                      const int m, const int c, const ibbox &ext,
                      const ibset &allactive) {
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
      CCTK_WARN(0, "internal error");
  } else {
    symtable = SymmetryTableHandleForGI(cctkGH, group);
    if (symtable < 0)
      CCTK_WARN(0, "internal error");
  }

  CCTK_INT symbnd[2 * dim];
  int const ierr =
      Util_TableGetIntArray(symtable, 2 * groupdim, symbnd, "symmetry_handle");
  if (ierr != 2 * groupdim)
    CCTK_WARN(0, "internal error");

  bool is_symbnd[2 * dim];
  for (int d = 0; d < 2 * groupdim; ++d) {
    is_symbnd[d] = symbnd[d] >= 0;
  }

  ivect lo = ext.lower();
  ivect hi = ext.upper();
  const ivect str = ext.stride();

  const b2vect obnds = vhh.at(m)->outer_boundaries(rl, c);
  const i2vect ghost_width = arrdata.at(group).at(m).dd->ghost_widths.AT(rl);

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

  ibset exts(ibbox(lo, hi, str));
  // do grid arrays have buffer zones?
  if (not output_buffer_points) {
    ivect loghosts(0);
    ivect highosts(0);
    // possibly re-add ghost points and symmetry points
    // TODO: come up with a nicer way to do this
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

      if (output_lower_ghosts) {
        loghosts[d] = ghost_width[0][d];
      }
      if (output_upper_ghosts) {
        highosts[d] = ghost_width[1][d];
      }
    }
    exts = exts & allactive.expand(loghosts, highosts);
  }

  return exts;
}

// Determine coordinates
void GetCoordinates(const cGH *const cctkGH, const int m,
                    const cGroup &groupdata, const ibset &exts,
                    CCTK_REAL &coord_time, vector<rvect> &coord_lower,
                    vector<rvect> &coord_upper) {
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

  coord_lower.reserve(exts.setsize());
  coord_upper.reserve(exts.setsize());

  for (ibset::const_iterator ext = exts.begin(); ext != exts.end(); ++ext) {
    coord_lower.push_back(global_lower + coord_delta * rvect(ext->lower()));
    coord_upper.push_back(global_lower + coord_delta * rvect(ext->upper()));
  }
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

// Output
template <int outdim>
int IOHDF5<outdim>::WriteHDF5(const cGH *cctkGH, hid_t &file, hid_t &indexfile,
                              vector<gdata *> const gfdatas,
                              const bbox<int, dim> &gfext, const int vi,
                              const vect<int, dim> &org,
                              const vect<int, outdim> &dirs, const int rl,
                              const int ml, const int m, const int c,
                              const int output_component, const int tl,
                              const CCTK_REAL coord_time,
                              const vect<CCTK_REAL, dim> &coord_lower,
                              const vect<CCTK_REAL, dim> &coord_upper) {
  DECLARE_CCTK_PARAMETERS;

  assert(outdim <= dim);

  cGroup groupdata;
  int const gi = CCTK_GroupIndexFromVarI(vi);
  assert(gi >= 0);
  int const ierr = CCTK_GroupData(gi, &groupdata);
  assert(not ierr);

  // boolean that says if we are doing 1D-diagonal output
  // This is not beautiful, but works for the moment
  bool const diagonal_output = outdim == 1 and dirs[0] == 3;

  // Check whether the output bbox overlaps
  // with the extent of the data to be output
  // FIXME: move this check up in the call stack
  bool output_bbox_overlaps_data_extent;
  if (not diagonal_output) {

    const vect<int, outdim> lo = gfext.lower()[dirs];
    const vect<int, outdim> up = gfext.upper()[dirs];
    const vect<int, outdim> str = gfext.stride()[dirs];
    const bbox<int, outdim> ext(lo, up, str);

    // Check whether the output origin is contained in the extent
    // of the data that should be output
    ivect org1(org);
    for (int d = 0; d < outdim; ++d)
      org1[dirs[d]] = ext.lower()[d];
    output_bbox_overlaps_data_extent = gfext.contains(org1);

  } else {

    gh const &hh = *vhh.at(m);
    ibbox const &base = hh.baseextents.at(mglevel).at(reflevel);

    assert(base.stride()[0] == base.stride()[1] and
           base.stride()[0] == base.stride()[2]);

    // Check if any point on the diagonal is in our gf's extent
    output_bbox_overlaps_data_extent = false;
    for (int i = maxval(base.lower()); i <= minval(base.upper());
         i += base.stride()[0]) {

      ivect const pos = ivect(i, i, i);
      output_bbox_overlaps_data_extent |= gfext.contains(pos);
    }
  }
  // Shortcut if there is nothing to output
  if (not output_bbox_overlaps_data_extent) {
    return 0;
  }

  int error_count = 0;

  ostringstream datasetname_suffix;
  datasetname_suffix << " it=" << cctkGH->cctk_iteration << " tl=" << tl;
  if (mglevels > 1)
    datasetname_suffix << " ml=" << ml;
  if (groupdata.grouptype == CCTK_GF) {
    if (maps > 1)
      datasetname_suffix << " m=" << m;
    datasetname_suffix << " rl=" << rl;
  }
  if (arrdata.at(gi).at(m).dd->light_boxes.at(ml).at(rl).size() > 1 and
      groupdata.disttype != CCTK_DISTRIB_CONSTANT) {
    datasetname_suffix << " c=" << output_component;
  }

  // enable compression and checksums if requested
  hid_t plist;
  HDF5_ERROR(plist = H5Pcreate(H5P_DATASET_CREATE));
  // enable checksums if requested
  if (use_checksums) {
    HDF5_ERROR(H5Pset_filter(plist, H5Z_FILTER_FLETCHER32, 0, 0, 0));
  }

  // enable datatype conversion if requested
  const hid_t mem_type = CCTKtoHDF5_Datatype(cctkGH, groupdata.vartype, false);
  const hid_t slice_type =
      CCTKtoHDF5_Datatype(cctkGH, groupdata.vartype, out_single_precision);

  if (not diagonal_output) { // not outputting the diagonal

    // there are currently three "extent" variables in use:
    // data_ext - the extent of the component in memory
    // gfext    - the bounding box of data to be written
    // ext      - gfext in the output dimension(s)
    // data_ext and gfext are used to construct the hyperslab location and
    // size, ext is just a shorthand.
    // TODO: maybe transfer only ext instead of data_ext in copy_from
    const vect<int, outdim> lo = gfext.lower()[dirs];
    const vect<int, outdim> up = gfext.upper()[dirs];
    const vect<int, outdim> str = gfext.stride()[dirs];
    const bbox<int, outdim> ext(lo, up, str);
    const dh *const dd = arrdata.at(gi).at(m).dd;
    const ibbox &data_ext = dd->light_boxes.at(ml).at(rl).at(c).exterior;

    // Check whether the output origin is contained in the extent of
    // the data that should be output
    ivect org1(org);
    for (int d = 0; d < outdim; ++d)
      org1[dirs[d]] = ext.lower()[d];
    assert(gfext.contains(org1));

    // HDF5 wants ranks to be >= 1
    const int rank = outdim > 0 ? outdim : 1;
    vector<hsize_t> mem_shape(dim);
    vector<hsize_t> slice_shape(rank, 1);
    hsize_t num_elems = 1;
    for (int d = 0; d < dim; d++) {
      mem_shape[dim - 1 - d] = data_ext.shape()[d] / data_ext.stride()[d];
      if (d < outdim) {
        slice_shape[outdim - 1 - d] = ext.shape()[d] / ext.stride()[d];
        num_elems *= slice_shape[outdim - 1 - d];
      }
    }

    ivect slice_lower(org - data_ext.lower());
    for (int d = 0; d < outdim; d++) {
      slice_lower[dirs[d]] = ext.lower()[d] - data_ext.lower()[dirs[d]];
    }
    ivect slice_upper(slice_lower);
    for (int d = 0; d < outdim; d++) {
      slice_upper[dirs[d]] = ext.upper()[d] - data_ext.lower()[dirs[d]];
    }
    slice_lower /= gfext.stride();
    slice_upper /= gfext.stride();

    slice_start_size_t slice_start[dim];
    hsize_t slice_count[dim];
    for (int d = 0; d < dim; d++) {
      slice_start[dim - 1 - d] = slice_lower[d];
      slice_count[dim - 1 - d] = slice_upper[d] - slice_lower[d] + 1;
    }
    const hsize_t size = num_elems * H5Tget_size(slice_type);
    if (compression_level and size > hsize_t(minimum_size_for_compression)) {
      HDF5_ERROR(H5Pset_shuffle(plist));
      HDF5_ERROR(H5Pset_deflate(plist, compression_level));
    }
    if (compression_level or use_checksums) {
      HDF5_ERROR(H5Pset_chunk(plist, slice_shape.size(), &slice_shape[0]));
    }
    hid_t slice_space, mem_space;
    HDF5_ERROR(slice_space =
                   H5Screate_simple(slice_shape.size(), &slice_shape[0], NULL));
    HDF5_ERROR(mem_space =
                   H5Screate_simple(mem_shape.size(), &mem_shape[0], NULL));
    HDF5_ERROR(H5Sselect_hyperslab(mem_space, H5S_SELECT_SET, slice_start, NULL,
                                   slice_count, NULL));

    vector<int> iorigin(rank, 0);
    vector<int> ioffset(rank, 0);
    vector<int> ioffsetdenom(rank, 1);
    vector<double> delta(rank, 0), origin(rank, 0);
    vector<int> bbox(2 * rank, 0), nghostzones(rank, 0);
    for (int d = 0; d < outdim; d++) {
      assert(gfext.upper()[dirs[d]] - gfext.lower()[dirs[d]] >= 0);
      iorigin[d] = ext.lower()[d];
      delta[d] = 0;
      origin[d] = coord_lower[dirs[d]];
      if (gfext.upper()[dirs[d]] - gfext.lower()[dirs[d]] > 0) {
        delta[d] = (coord_upper[dirs[d]] - coord_lower[dirs[d]]) /
                   (gfext.upper()[dirs[d]] - gfext.lower()[dirs[d]]) *
                   gfext.stride()[dirs[d]];
        origin[d] += (org1[dirs[d]] - gfext.lower()[dirs[d]]) * delta[d];
        ioffsetdenom[d] = gfext.stride()[dirs[d]];
        ioffset[d] =
            (iorigin[d] % gfext.stride()[dirs[d]] + gfext.stride()[dirs[d]]) %
            gfext.stride()[dirs[d]];
        assert((iorigin[d] - ioffset[d]) % gfext.stride()[dirs[d]] == 0);
        iorigin[d] = (iorigin[d] - ioffset[d]) / gfext.stride()[dirs[d]];
      }
    }
    string active;
    {
      // Determine extent of hyperslab that is output
      ivect lo = gfext.lower();
      ivect up = gfext.upper();
      ivect str = gfext.stride();
      for (int d = 0; d < dim; ++d) {
        bool isoutdir = false;
        for (int e = 0; e < outdim; ++e)
          isoutdir |= d == dirs[e];
        if (!isoutdir) {
          lo[d] = org[d];
          up[d] = org[d];
        }
      }
      const ibbox outputslab(lo, up, str);
      // Intersect active region with this hyperslab
      const int lc = vhh.at(m)->get_local_component(rl, c);
      const ibset &active0 = vdd.at(m)->level_boxes.at(ml).at(rl).active;
      const ibset active1 = active0 & outputslab;
      // Reduce dimensionality of active region
      bboxset<int, outdim> active2;
      for (ibset::const_iterator bi = active1.begin(), be = active1.end();
           bi != be; ++bi) {
        const ibbox &box0 = *bi;
        const vect<int, outdim> lo = box0.lower()[dirs];
        const vect<int, outdim> up = box0.upper()[dirs];
        const vect<int, outdim> str = box0.stride()[dirs];
        const ::bbox<int, outdim> box(lo, up, str);
        active2 += box;
      }
      ostringstream buf;
      buf << active2;
      active = buf.str();
    }

    // store cctk_bbox and cctk_nghostzones (for grid arrays only)
    if (groupdata.grouptype != CCTK_SCALAR) {
      const b2vect obnds = vhh.at(m)->outer_boundaries(rl, c);
      const i2vect ghost_width = arrdata.at(gi).at(m).dd->ghost_widths.AT(rl);
      for (int d = 0; d < outdim; d++) {
        nghostzones[d] = output_ghost_points ? ghost_width[0][dirs[d]] : 0;
        assert(all(ghost_width[0] == ghost_width[1]));

        bbox[2 * d] = obnds[0][dirs[d]];
        bbox[2 * d + 1] = obnds[1][dirs[d]];
      }
    }

    // now loop over all variables
    for (size_t n = 0; n < gfdatas.size(); n++) {

      // create a unique name for this variable's dataset
      char *fullname = CCTK_FullName(vi + n);
      string datasetname(fullname);
      datasetname.append(datasetname_suffix.str());

      // remove an already existing dataset of the same name
      ioRequest *request = slice_requests.at(vi + n);
      if (not request) {
#ifdef IOUTIL_PARSER_HAS_OUT_DT
        request = IOUtil_DefaultIORequest(cctkGH, vi + n, 1, -1.0);
#else
        request = IOUtil_DefaultIORequest(cctkGH, vi + n, 1);
#endif
      }
      if (request->check_exist) {
        H5E_BEGIN_TRY {
          H5Gunlink(file, datasetname.c_str());
          if (indexfile != -1)
            H5Gunlink(indexfile, datasetname.c_str());
        }
        H5E_END_TRY;
      }
      // free I/O request structure
      if (request != slice_requests.at(vi + n)) {
        IOUtil_FreeIORequest(&request);
      }

      // write the dataset
      hid_t dataset, index_dataset;
      HDF5_ERROR(dataset = H5Dcreate(file, datasetname.c_str(), slice_type,
                                     slice_space, plist));
      if (indexfile != -1) {
        HDF5_ERROR(index_dataset = H5Dcreate(indexfile, datasetname.c_str(),
                                             slice_type, slice_space, plist));
      }

      HDF5_ERROR(H5Dwrite(dataset, mem_type, mem_space, H5S_ALL, H5P_DEFAULT,
                          gfdatas[n]->storage()));
      error_count += AddSliceAttributes(
          cctkGH, fullname, rl, ml, m, tl, origin, delta, iorigin, ioffset,
          ioffsetdenom, bbox, nghostzones, active, dataset, slice_shape, false);
      HDF5_ERROR(H5Dclose(dataset));

      if (indexfile != -1) {
        error_count += AddSliceAttributes(
            cctkGH, fullname, rl, ml, m, tl, origin, delta, iorigin, ioffset,
            ioffsetdenom, bbox, nghostzones, active, index_dataset, slice_shape,
            true);
        HDF5_ERROR(H5Dclose(index_dataset));
      }
      free(fullname);

      io_bytes +=
          H5Sget_simple_extent_npoints(slice_space) * H5Tget_size(slice_type);
    } // for n

    HDF5_ERROR(H5Sclose(mem_space));
    HDF5_ERROR(H5Sclose(slice_space));

  } else { // taking care of the diagonal

    const ivect lo = gfext.lower();
    const ivect up = gfext.upper();
    const ivect str = gfext.stride();
    const ibbox ext(lo, up, str);

    gh const &hh = *vhh.at(m);
    ibbox const &base = hh.baseextents.at(mglevel).at(reflevel);

    assert(base.stride()[0] == base.stride()[1] and
           base.stride()[0] == base.stride()[2]);

    // count the number of points on the diagonal
    hsize_t npoints = 0;
    for (int i = maxval(base.lower()); i <= minval(base.upper());
         i += base.stride()[0]) {
      if (gfext.contains(i)) {
        ++npoints;
      }
    }
    assert(npoints > 0);

    // allocate a contiguous buffer for the diagonal points
    vector<char> buffer(CCTK_VarTypeSize(groupdata.vartype) * npoints *
                        gfdatas.size());

    // copy diagonal points into contiguous buffer
    hsize_t offset = 0;
    for (int i = maxval(base.lower()); i <= minval(base.upper());
         i += base.stride()[0]) {
      ivect const pos = ivect(i, i, i);
      if (gfext.contains(pos)) {
        for (size_t n = 0; n < gfdatas.size(); n++) {
          switch (specific_cactus_type(groupdata.vartype)) {
#define TYPECASE(N, T)                                                         \
  case N: {                                                                    \
    T *typed_buffer = (T *)&buffer.front();                                    \
    typed_buffer[offset + n * npoints] =                                       \
        (*(const data<T> *)gfdatas.at(n))[pos];                                \
    break;                                                                     \
  }
#include "typecase.hh"
#undef TYPECASE
          }
        }
        ++offset;
      }
    }
    assert(offset == npoints);

    const hsize_t size = npoints * H5Tget_size(slice_type);
    if (compression_level and size > hsize_t(minimum_size_for_compression)) {
      HDF5_ERROR(H5Pset_shuffle(plist));
      HDF5_ERROR(H5Pset_deflate(plist, compression_level));
    }
    if (compression_level or use_checksums) {
      HDF5_ERROR(H5Pset_chunk(plist, 1, &npoints));
    }
    hid_t slice_space;
    HDF5_ERROR(slice_space = H5Screate_simple(1, &npoints, NULL));

    // loop over all variables and write out diagonals
    for (size_t n = 0; n < gfdatas.size(); n++) {

      // create a unique name for this variable's dataset
      char *fullname = CCTK_FullName(vi + n);
      string datasetname(fullname);
      free(fullname);
      datasetname.append(datasetname_suffix.str());

      // remove an already existing dataset of the same name
      ioRequest *request = slice_requests.at(vi + n);
      if (not request) {
#ifdef IOUTIL_PARSER_HAS_OUT_DT
        request = IOUtil_DefaultIORequest(cctkGH, vi + n, 1, -1.0);
#else
        request = IOUtil_DefaultIORequest(cctkGH, vi + n, 1);
#endif
      }
      if (request->check_exist) {
        H5E_BEGIN_TRY { H5Gunlink(file, datasetname.c_str()); }
        H5E_END_TRY;
      }
      // free I/O request structure
      if (request != slice_requests.at(vi + n)) {
        IOUtil_FreeIORequest(&request);
      }

      // write the dataset
      hid_t dataset;
      HDF5_ERROR(dataset = H5Dcreate(file, datasetname.c_str(), slice_type,
                                     slice_space, plist));
      HDF5_ERROR(H5Dwrite(dataset, mem_type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                          &buffer.front() + n * npoints * gfdatas.size()));
      HDF5_ERROR(H5Dclose(dataset));

      io_bytes +=
          H5Sget_simple_extent_npoints(slice_space) * H5Tget_size(slice_type);
    }

    HDF5_ERROR(H5Sclose(slice_space));

  } // if(not diagonal_output)

  HDF5_ERROR(H5Pclose(plist));

  return error_count;
}

// which processor serves as IO processor for this group
// TODO: see if IOUtil offers and official way to get this information
template <int outdim> int IOHDF5<outdim>::IOProcForProc(int proc) {
  // according to IOUtil::SetupGH a proc with rank % nioprocs == 0 is an
  // ioproc
  return (proc / ioproc_every) * ioproc_every;
}

// Explicit instantiation for all slice output dimensions
template class IOHDF5<0>;
template class IOHDF5<1>;
template class IOHDF5<2>;
template class IOHDF5<3>;

} // namespace CarpetIOHDF5
