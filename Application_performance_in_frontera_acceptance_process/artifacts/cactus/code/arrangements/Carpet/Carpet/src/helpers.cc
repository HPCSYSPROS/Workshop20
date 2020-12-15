#include <cctk.h>
#include <cctk_FortranString.h>
#include <cctk_Parameters.h>
#include <util_Table.h>

#include <algorithm>
#include <cassert>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#ifdef CCTK_MPI
#include <mpi.h>
#else
#include "nompi.h"
#endif

#include <defs.hh>
#include <dist.hh>
#include <ggf.hh>
#include <typeprops.hh>

#include <carpet.hh>

namespace Carpet {

using namespace std;

// Get Carpet's GH extension
CarpetGH const *GetCarpetGH(const cGH *const cgh) {
  assert(cgh);
  return &carpetGH;
}

// Get current regridding epoch
extern "C" CCTK_INT
Carpet_GetRegriddingEpoch(CCTK_POINTER_TO_CONST const cctkGH) {
  return regridding_epoch;
}

// Get current regridding epochs
extern "C" CCTK_INT
Carpet_GetRegriddingEpochs(CCTK_POINTER_TO_CONST const cctkGH,
                           CCTK_INT const size,
                           CCTK_INT *restrict const epochs) {
  int const maxrl = min(size, CCTK_INT(level_regridding_epochs.size()));
  for (int rl = 0; rl < maxrl; ++rl) {
    epochs[rl] = level_regridding_epochs.AT(rl);
  }
  return CCTK_INT(level_regridding_epochs.size());
}

// Get current refinement level
extern "C" CCTK_INT
Carpet_GetRefinementLevel(CCTK_POINTER_TO_CONST const cctkGH) {
  return reflevel;
}

// Get number of refinement levels
extern "C" CCTK_INT
Carpet_GetRefinementLevels(CCTK_POINTER_TO_CONST const cctkGH) {
  return reflevels;
}

// Get maximum number of refinement levels
extern "C" CCTK_INT
Carpet_GetMaxRefinementLevels(CCTK_POINTER_TO_CONST const cctkGH) {
  return maxreflevels;
}

// Get current local component
extern "C" CCTK_INT
Carpet_GetLocalComponent(CCTK_POINTER_TO_CONST const cctkGH) {
  assert(reflevel >= 0);
  assert(map >= 0);
  return local_component;
}

// Get number of local components
extern "C" CCTK_INT
Carpet_GetLocalComponents(CCTK_POINTER_TO_CONST const cctkGH) {
  assert(reflevel >= 0);
  assert(map >= 0);
  return vhh.AT(map)->local_components(reflevel);
}

// Get current map level
extern "C" CCTK_INT Carpet_GetMap(CCTK_POINTER_TO_CONST const cctkGH) {
  return map;
}

// Get number of maps
extern "C" CCTK_INT Carpet_GetMaps(CCTK_POINTER_TO_CONST const cctkGH) {
  return maps;
}

// Get current map level
extern "C" CCTK_INT Carpet_GetTimeLevel(CCTK_POINTER_TO_CONST const cctkGH) {
  return timelevel;
}

// Get current map level
extern "C" CCTK_INT
Carpet_GetTimeLevelOffset(CCTK_POINTER_TO_CONST const cctkGH) {
  return timelevel_offset;
}

// Enable or disable prolongating
CCTK_INT CarpetEnableProlongating(const CCTK_INT flag) {
  assert(flag == 0 or flag == 1);
  do_prolongate = flag;
  if (do_prolongate) {
    Checkpoint("Prolongating enabled");
  } else {
    Checkpoint("Prolongating disabled");
  }
  return 0;
}

CCTK_INT CarpetQueryProlongating() { return do_prolongate; }

// Get pointer to grid variable for a specific map and refinement level
CCTK_POINTER
Carpet_VarDataPtrI(CCTK_POINTER_TO_CONST const cctkGH_, CCTK_INT const m,
                   CCTK_INT const rl, CCTK_INT const c, CCTK_INT const tl,
                   CCTK_INT const varindex) {
  assert(cctkGH_);
  cGH const *const cctkGH = static_cast<cGH const *>(cctkGH_);
  assert(varindex >= 0 and varindex < CCTK_NumVars());
  int const groupindex = CCTK_GroupIndexFromVarI(varindex);
  assert(groupindex >= 0);
  int const grouptype = CCTK_GroupTypeI(groupindex);
  assert(mglevel >= 0);
  if (grouptype == CCTK_GF) {
    assert(m >= 0 and m < maps);
    assert(rl >= 0 and rl < reflevels);
    assert(c >= 0 and
           c < arrdata.AT(groupindex).AT(m).hh->components(reflevel));
    assert(arrdata.AT(groupindex).AT(m).hh->is_local(rl, c));
  } else {
    assert(m == 0);
    assert(rl == 0);
    assert(c == CCTK_MyProc(NULL));
  }
  int const maxtls = CCTK_MaxActiveTimeLevelsGI(cctkGH, groupindex);
  assert(tl >= 0 and tl < maxtls);

  int const activetimelevels =
      groupdata.AT(groupindex).activetimelevels.AT(mglevel).AT(rl);
  if (tl < activetimelevels) {
    int const var = varindex - CCTK_FirstVarIndexI(groupindex);
    ggf *const ff = arrdata.AT(groupindex).AT(m).data.AT(var);
    int const lc = arrdata.AT(groupindex).AT(m).hh->get_local_component(rl, c);
    assert(lc >= 0 and
           lc < arrdata.AT(groupindex).AT(m).hh->local_components(rl));
    gdata *const data = ff->data_pointer(tl, rl, lc, mglevel);
    return data->storage();
  } else {
    return NULL;
  }
}

// Multi-Model
CCTK_POINTER_TO_CONST
Carpet_GetMPICommUniverse(CCTK_POINTER_TO_CONST const cctkGH) {
  assert(comm_universe != MPI_COMM_NULL);
  return &comm_universe;
}

CCTK_POINTER_TO_CONST
Carpet_GetMPICommWorld(CCTK_POINTER_TO_CONST const cctkGH) {
  assert(comm_world != MPI_COMM_NULL);
  return &comm_world;
}

// Hosts
extern "C" CCTK_INT Carpet_MyHost(CCTK_POINTER_TO_CONST const cctkGH_) {
  return HostId(dist::rank());
}

extern "C" CCTK_INT Carpet_nHosts(CCTK_POINTER_TO_CONST const cctkGH_) {
  return HostNames().size();
}

extern "C" CCTK_INT Carpet_nProcsOnHost(CCTK_POINTER_TO_CONST const cctkGH_,
                                        CCTK_INT const host) {
  return HostProcs(host).size();
}

extern "C" CCTK_INT Carpet_ProcsOnHost(CCTK_POINTER_TO_CONST const cctkGH_,
                                       CCTK_INT const host, CCTK_INT procs[],
                                       CCTK_INT const nprocs) {
  int const nprocs1 = HostProcs(host).size();
  for (int p = 0; p < nprocs1; ++p) {
    if (p >= nprocs)
      break;
    procs[p] = HostProcs(host).AT(p);
  }
  return nprocs1;
}

// Coordinates

CCTK_INT
Carpet_GetCoordRange(CCTK_POINTER_TO_CONST const cctkGH_, CCTK_INT const m,
                     CCTK_INT const ml, CCTK_INT const size,
                     CCTK_INT *const gsh, CCTK_REAL *const lower,
                     CCTK_REAL *const upper, CCTK_REAL *const delta) {
  cGH const *const cctkGH = static_cast<cGH const *>(cctkGH_);
  assert(cctkGH);
  assert(m >= 0 and m < maps);
  assert(ml >= 0 and ml < mglevels);
  assert(size >= 0);
  assert(not size or gsh);
  assert(not size or lower);
  assert(not size or upper);
  assert(not size or delta);

  assert(size == dim);

  ibbox const &baseext = vhh.AT(m)->baseextents.AT(ml).AT(0);
  ivect const igsh = baseext.shape() / baseext.stride();

  for (int d = 0; d < dim; ++d) {
    gsh[d] = igsh[d];
    lower[d] = origin_space.AT(m).AT(ml)[d];
    delta[d] = delta_space.AT(m)[d] * mglevelfact;
    upper[d] = lower[d] + delta[d] * (gsh[d] - 1);
  }

  return 0;
}

// Communication

int Barrier(const cGH *cgh) {
  const void *dummy = &dummy;
  dummy = &cgh;

  MPI_Barrier(dist::comm());
  return 0;
}

int NamedBarrier(const cGH *const cgh, const unsigned int id,
                 const char *name) {
  const void *dummy = &dummy;
  dummy = &cgh;

  // Prevent recursive calls
  static bool inside = false;
  if (inside)
    return 0;
  inside = true;

  Checkpoint("About to Barrier %ud \"%s\"", id, name);
  dist::barrier(dist::comm(), id, name);
  Checkpoint("Finished Barrier");

  inside = false;
  return 0;
}

int Exit(const cGH *cgh, int retval) {
  NamedBarrier(cgh, 792686462, "Carpet::Exit");
  dist::finalize();
  exit(retval);
  return -1;
}

int Abort(const cGH *cgh, int retval) {
  void *dummy = &dummy;
  dummy = &cgh;

  assert(0);

  MPI_Abort(MPI_COMM_WORLD, retval);
  exit(retval);
  return -1;
}

int MyProc(const cGH *cgh) {
  // This may be called very early, before dist:comm() is valid
  int rank;
  MPI_Comm_rank(dist::goodcomm(), &rank);
  return rank;
}

int nProcs(const cGH *cgh) {
  // This may be called very early, before dist:comm() is valid
  int size;
  MPI_Comm_size(dist::goodcomm(), &size);
  return size;
}

MPI_Comm CarpetMPIComm() { return dist::comm(); }

// Datatypes

MPI_Datatype CarpetMPIDatatype(const int vartype) {
  switch (specific_cactus_type(vartype)) {
#define TYPECASE(N, T)                                                         \
  case N: {                                                                    \
    T dummy;                                                                   \
    return dist::mpi_datatype(dummy);                                          \
  }
#include "typecase.hh"
#undef TYPECASE
  }
  CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
             "Carpet does not support the variable type %d.", vartype);
  // notreached
  return MPI_CHAR;
}

MPI_Datatype CarpetSimpleMPIDatatype(const int vartype) {
  switch (vartype) {
  case CCTK_VARIABLE_COMPLEX:
    return CarpetMPIDatatype(CCTK_VARIABLE_REAL);
#ifdef HAVE_CCTK_COMPLEX8
  case CCTK_VARIABLE_COMPLEX8:
    return CarpetMPIDatatype(CCTK_VARIABLE_REAL4);
#endif
#ifdef HAVE_CCTK_COMPLEX16
  case CCTK_VARIABLE_COMPLEX16:
    return CarpetMPIDatatype(CCTK_VARIABLE_REAL8);
#endif
#ifdef HAVE_CCTK_COMPLEX32
  case CCTK_VARIABLE_COMPLEX32:
    return CarpetMPIDatatype(CCTK_VARIABLE_REAL16);
#endif
  }
  // default
  return CarpetMPIDatatype(vartype);
}

int CarpetSimpleMPIDatatypeLength(const int vartype) {
  switch (vartype) {
  case CCTK_VARIABLE_COMPLEX:
#ifdef HAVE_CCTK_COMPLEX8
  case CCTK_VARIABLE_COMPLEX8:
#endif
#ifdef HAVE_CCTK_COMPLEX16
  case CCTK_VARIABLE_COMPLEX16:
#endif
#ifdef HAVE_CCTK_COMPLEX32
  case CCTK_VARIABLE_COMPLEX32:
#endif
    return 2;
  }
  // default
  return 1;
}

// Timelevels

int min_timelevel(const checktimes where, const int num_tl,
                  const bool persistent) {
  assert(num_tl > 0);
  switch (where) {
  case currenttime:
    return 0;
  case currenttimebutnotifonly:
    // don't include current time if there is only one (persistent)
    // time level
    return (not persistent or num_tl > 1) ? 0 : 1;
  case previoustime:
    return 1;
  case allbutlasttime:
    // do include current time if there is only one time level
    return 0;
  case allbutcurrenttime:
    return 1;
  case alltimes:
    return 0;
  default:
    assert(0);
  }
  return -999;
}

int max_timelevel(const checktimes where, const int num_tl,
                  const bool persistent) {
  assert(num_tl > 0);
  switch (where) {
  case currenttime:
    return 0;
  case currenttimebutnotifonly:
    return 0;
  case previoustime:
    return num_tl > 1 ? 1 : 0;
  case allbutlasttime:
    return num_tl - 2;
  case allbutcurrenttime:
    return num_tl - 1;
  case alltimes:
    return num_tl - 1;
  default:
    assert(0);
  }
  return -999;
}

// Diagnostic output

static void prepend_id(char *const msg, size_t const msglen) {
  bool needspace = false;
  if (mglevel != -1) {
    snprintf(msg + strlen(msg), msglen - strlen(msg), "[ml=%d]", mglevel);
    needspace = true;
  }
  if (reflevel != -1) {
    snprintf(msg + strlen(msg), msglen - strlen(msg), "[rl=%d]", reflevel);
    needspace = true;
  }
  if (map != -1) {
    snprintf(msg + strlen(msg), msglen - strlen(msg), "[m=%d]", map);
    needspace = true;
  }
  if (component != -1) {
    snprintf(msg + strlen(msg), msglen - strlen(msg), "[c=%d,lc=%d]", component,
             local_component);
    needspace = true;
  }
  if (timelevel != -1) {
    snprintf(msg + strlen(msg), msglen - strlen(msg), "[tl=%d]", timelevel);
    needspace = true;
  }
  if (needspace) {
    snprintf(msg + strlen(msg), msglen - strlen(msg), " ");
  }
}

void Output(const char *fmt, ...) {
  DECLARE_CCTK_PARAMETERS;
  va_list args;
  char msg[1000];
  snprintf(msg, sizeof msg, "%s", "");
  prepend_id(msg + strlen(msg), sizeof msg - strlen(msg));
  va_start(args, fmt);
  vsnprintf(msg + strlen(msg), sizeof msg - strlen(msg), fmt, args);
  va_end(args);
  CCTK_INFO(msg);
  if (barriers) {
    NamedBarrier(NULL, 988121033, "Carpet::Output");
  }
}

void Waypoint(const char *fmt, ...) {
  DECLARE_CCTK_PARAMETERS;
  if (verbose or veryverbose) {
    va_list args;
    char msg[1000];
    snprintf(msg, sizeof msg, "%s", "");
    prepend_id(msg + strlen(msg), sizeof msg - strlen(msg));
    va_start(args, fmt);
    vsnprintf(msg + strlen(msg), sizeof msg - strlen(msg), fmt, args);
    va_end(args);
    CCTK_INFO(msg);
  }
  if (barriers) {
    NamedBarrier(NULL, 471561432, "Carpet::Waypoint");
  }
}

void Checkpoint(const char *fmt, ...) {
  DECLARE_CCTK_PARAMETERS;
  if (veryverbose) {
    va_list args;
    char msg[1000];
    snprintf(msg, sizeof msg, "%s", "");
    prepend_id(msg + strlen(msg), sizeof msg - strlen(msg));
    va_start(args, fmt);
    vsnprintf(msg + strlen(msg), sizeof msg - strlen(msg), fmt, args);
    va_end(args);
    CCTK_INFO(msg);
  }
  if (barriers) {
    NamedBarrier(NULL, 572893143, "Carpet::Checkpoint");
  }
}

void UnsupportedVarType(const int vindex) {
  assert(vindex >= 0 and vindex < CCTK_NumVars());
  CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
              "Carpet does not support the type of the variable \"%s\".\n"
              "Either enable support for this type, "
              "or change the type of this variable.",
              CCTK_FullName(vindex));
}

bool IsMap0Group(const int gindex) {
  assert(gindex >= 0 and gindex < CCTK_NumGroups());
  const int table = CCTK_GroupTagsTableI(gindex);
  assert(table >= 0);
  CCTK_INT map0group;
  int status = Util_TableGetInt(table, &map0group, "map0group");
  if (status == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
    map0group = false;
    status = 1;
  }
  if (status != 1) {
    CCTK_ERROR("illegal table value");
  }
  return map0group;
}

} // namespace Carpet
