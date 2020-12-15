#include <variables.hh>

namespace Carpet {

using namespace std;

// Handle from CCTK_RegisterGHExtension
int GHExtension;

// Maximum number of refinement levels
int maxreflevels;

// Maximum number of time levels
int maxtimelevels;

// Timelevels
int timelevels;

// Refinement levels
int reflevels;

// Temporal refinement factors over the coarsest grid
vector<int> timereffacts;

// Spatial refinement factors over the coarsest grid
vector<vect<int, dim> > spacereffacts;

// Maximum refinement factors on finest possible grid
int maxtimereflevelfact;
vect<int, dim> maxspacereflevelfact;

// Base multigrid level
int basemglevel;

// Multigrid levels
int mglevels;

// Multigrid factor
int mgfact;

// Multigrid factor on coarsest grid
int maxmglevelfact;

// Maps
int maps;

// Current position on the grid hierarchy
int reflevel;
int mglevel;
int mc_grouptype;
int map;
int component;
int local_component;
int timelevel, timelevel_offset;

// Current refinement factors
int timereflevelfact;
vect<int, dim> spacereflevelfact;

// Current multigrid factor
int mglevelfact;

// Carpet's GH
CarpetGH carpetGH;

// Times and spaces on the refinement levels
CCTK_REAL global_time;
// vector<vector<CCTK_REAL> > leveltimes; // [mglevel][reflevel]
CCTK_REAL delta_time;

vector<vector<vect<CCTK_REAL, dim> > > origin_space; // [map][mglevel]
vector<vect<CCTK_REAL, dim> > delta_space;           // [map]

vector<domainspec> domainspecs; // [map]

// Is this the time for a global mode call?
bool do_meta_mode;
bool do_early_meta_mode;
bool do_late_meta_mode;
bool do_global_mode;
bool do_early_global_mode;
bool do_late_global_mode;

// Can past time levels be accessed?
bool do_allow_past_timelevels;

// Is prolongation enabled?
// (This flag disables prolongation during MoL integration
// substeps.)
bool do_prolongate;

// Is tapering enabled?
// (This flag disables prolongation while the current refinement
// level is not aligned with the parent.)
bool do_taper;

// Should we warn about groups with insufficiently many time levels?
bool do_warn_about_storage;

// Are we in the analysis bin?
bool in_analysis_bin;

// Data for grid functions

// The grid hierarchy
vector<gh *> vhh; // [map]
vector<dh *> vdd; // [map]
th *tt;
int regridding_epoch;
vector<int> level_regridding_epochs;

// Data for the groups
vector<groupdesc> groupdata; // [group]

// Data for everything
vector<vector<arrdesc> > arrdata; // [group][map]

// MPI Communicators
MPI_Comm comm_universe = MPI_COMM_NULL;
MPI_Comm comm_world = MPI_COMM_NULL;

} // namespace Carpet
