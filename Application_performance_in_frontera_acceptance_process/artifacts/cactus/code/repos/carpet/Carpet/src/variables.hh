// It is assumed that each group has at least one map.  All arrays
// have exactly one map.  All maps have the same number of refinement
// levels.

// It is assumed that each group has at least one component.

// It is assumed that the number of components of all arrays is equal
// to the number of components of the grid functions, and that their
// distribution onto the processors is the same, and that all
// processors own the same number of components.

#ifndef VARIABLES_HH
#define VARIABLES_HH

#include <cctk.h>

#include <vector>

#ifdef CCTK_MPI
#include <mpi.h>
#else
#include "nompi.h"
#endif

#include <bbox.hh>
#include <data.hh>
#include <dh.hh>
#include <ggf.hh>
#include <gh.hh>
#include <operators.hh>
#include <th.hh>
#include <vect.hh>

#include "carpet_public.h"

namespace Carpet {

using namespace std;

// Handle from CCTK_RegisterGHExtension
extern int GHExtension;

// Maximum number of refinement levels
extern int maxreflevels;

// Maximum number of time levels
extern int maxtimelevels;

// Timelevels
extern int timelevels;

// Refinement levels
extern int reflevels;

#define CARPET_NEW_REFFACT
// Temporal refinement factors over the coarsest grid
extern vector<int> timereffacts;

// Spatial refinement factors over the coarsest grid
extern vector<vect<int, dim> > spacereffacts;

// Maximum refinement factors on finest possible grid
extern int maxtimereflevelfact;
extern vect<int, dim> maxspacereflevelfact;

// Base multigrid level
extern int basemglevel;

// Multigrid levels
extern int mglevels;

// Multigrid factor
extern int mgfact;

// Multigrid factor on coarsest grid
extern int maxmglevelfact;

// Maps
extern int maps;

// Current position on the grid hierarchy
extern int reflevel;
extern int mglevel;
extern int mc_grouptype; // -1, CCTK_SCALAR/CCTK_ARRAY, CCTK_GF
extern int map;
extern int component;
extern int local_component; // -1 for non-local
extern int timelevel, timelevel_offset;

// Current refinement factors
extern int timereflevelfact;
extern vect<int, dim> spacereflevelfact;

// Current multigrid factor
extern int mglevelfact;

// Carpet's GH
extern CarpetGH carpetGH;

// Times and spaces on the refinement levels
extern CCTK_REAL global_time;
// extern vector<vector<CCTK_REAL> > leveltimes; // [mglevel][reflevel]
extern CCTK_REAL delta_time;

extern vector<vector<vect<CCTK_REAL, dim> > > origin_space; // [map][mglevel]
extern vector<vect<CCTK_REAL, dim> > delta_space;           // [map]

// Domain extent, as used for the original grid setup
// TODO: Unify this with origin_space and delta_space,
//       possibly asserting constency at run time
struct domainspec {
  vect<CCTK_REAL, dim> exterior_min, exterior_max;
  vect<int, dim> npoints;
};
extern vector<domainspec> domainspecs; // [map]

// Is this the time for a global mode call?
extern bool do_meta_mode;
extern bool do_early_meta_mode;
extern bool do_late_meta_mode;
extern bool do_global_mode;
extern bool do_early_global_mode;
extern bool do_late_global_mode;

// Can past time levels be accessed?
extern bool do_allow_past_timelevels;

// Is prolongation enabled?
// (This flag disables prolongation during MoL integration
// substeps.)
extern bool do_prolongate;

// Is tapering enabled?
// (This flag disables prolongation while the current refinement
// level is not aligned with the parent.)
extern bool do_taper;

// Should we warn about groups with insufficiently many time levels?
extern bool do_warn_about_storage;

// Are we in the analysis bin?
extern bool in_analysis_bin;

// Data for grid functions

// The grid hierarchy
extern vector<gh *> vhh; // [map]
extern vector<dh *> vdd; // [map]
extern th *tt;
extern int regridding_epoch; // increases with each regridding
extern vector<int> level_regridding_epochs;

// Data for the groups
struct groupdesc {
  cGroupDynamicData info;
  operator_type transport_operator;      // prolongation and restriction
  vector<vector<int> > activetimelevels; // [mglevel][reflevel]
};
extern vector<groupdesc> groupdata; // [group]

// Data for everything
struct arrdesc {
  // points to hh etc. for GF, and is unique for SCALAR and ARRAY
  gh *hh;
  dh *dd;
  th *tt;
  vector<ggf *> data; // [var]
};
extern vector<vector<arrdesc> > arrdata; // [group][map]

// MPI Communicators
extern MPI_Comm comm_universe;
extern MPI_Comm comm_world;

} // namespace Carpet

#endif // !defined(VARIABLES_HH)
