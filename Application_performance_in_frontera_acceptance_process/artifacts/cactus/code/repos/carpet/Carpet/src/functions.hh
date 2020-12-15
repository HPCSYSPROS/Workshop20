#ifndef FUNCTIONS_HH
#define FUNCTIONS_HH

#include <cctk.h>
#include <cctk_Schedule.h>

#include <map>
#include <string>
#include <vector>

#ifdef CCTK_MPI
#include <mpi.h>
#else
#include "nompi.h"
#endif

#include <bbox.hh>
#include <dh.hh>
#include <gh.hh>
#include <vect.hh>

namespace Carpet {

using namespace std;

int SyncGroupsByDirI(const cGH *cctkGH, int num_groups, const int *groups,
                     const int *directions);
int EnableGroupComm(const cGH *cgh, const char *groupname);
int DisableGroupComm(const cGH *cgh, const char *groupname);
int EnableGroupStorage(const cGH *cgh, const char *groupname);
int DisableGroupStorage(const cGH *cgh, const char *groupname);
int GroupStorageIncrease(const cGH *cgh, int n_groups, const int *groups,
                         const int *timelevels, int *status);
int GroupStorageDecrease(const cGH *cgh, int n_groups, const int *groups,
                         const int *timelevels, int *status);
int QueryMaxTimeLevels(const cGH *cgh, int n_groups, const int *groups,
                       int *status);
int Barrier(const cGH *cgh);
int NamedBarrier(const cGH *cgh, unsigned int id, const char *name);
int Exit(const cGH *cgh, int retval);
int Abort(const cGH *cgh, int retval);
int MyProc(const cGH *cgh) CCTK_ATTRIBUTE_PURE;
int nProcs(const cGH *cgh) CCTK_ATTRIBUTE_PURE;
const int *ArrayGroupSizeB(const cGH *cgh, int dir, int group,
                           const char *groupname) CCTK_ATTRIBUTE_PURE;
int QueryGroupStorageB(const cGH *cgh, int group,
                       const char *groupname) CCTK_ATTRIBUTE_PURE;
int GroupDynamicData(const cGH *cgh, int group, cGroupDynamicData *data);

void Restrict(const cGH *cgh);

// Multi-Model
void SplitUniverse(MPI_Comm const world, string const model, MPI_Comm &comm,
                   bool verbose);

// Model id to model name
vector<string> const &ModelNames();
string ModelName(int id);

// Model name to model id
std::map<string, int> const &ModelMap();
int ModelId(string name);

// Processor to model id
vector<int> const &ModelIds();
int ModelId(int proc);

// Model id to processes
vector<vector<int> > const &ModelProcs();
vector<int> const &ModelProcs(int id);

// Host mapping
void DetermineHosts(string host, bool verbose);

// Host id to host name
vector<string> const &HostNames();
string HostName(int id);

// Host name to host id
std::map<string, int> const &HostMap();
int HostId(string name);

// Process to host id
vector<int> const &HostIds();
int HostId(int proc);

// Host id to processes
vector<vector<int> > const &HostProcs();
vector<int> const &HostProcs(int id);

extern "C" {
CCTK_POINTER_TO_CONST
Carpet_GetMPICommUniverse(CCTK_POINTER_TO_CONST cctkGH);
CCTK_POINTER_TO_CONST
Carpet_GetMPICommWorld(CCTK_POINTER_TO_CONST cctkGH);
CCTK_INT
Carpet_GetCoordRange(CCTK_POINTER_TO_CONST const cctkGH_, CCTK_INT const m,
                     CCTK_INT const ml, CCTK_INT const size,
                     CCTK_INT *const gsh, CCTK_REAL *const lower,
                     CCTK_REAL *const upper, CCTK_REAL *const delta);
}

// Helpers for storage
void GroupsStorageCheck(cGH const *const cctkGH);

// Helpers for recomposing the grid hierarchy
void RegridMap(cGH const *cctkGH, int m, gh::rregs const &supeerregss,
               gh::mregs const &regsss, bool do_init);
void PostRegrid(cGH const *cctkGH);
bool Recompose(cGH const *cctkGH, int rl, bool do_init);
void RegridFree(cGH const *cctkGH, bool do_init);

void CheckRegions(gh::mregs const &regsss);

void OutputSuperregions(cGH const *cctkGH, int m, gh const &hh,
                        gh::rregs const &superregss);
void OutputGrids(cGH const *cctkGH, int const m, gh const &hh, dh const &dd);

void OutputGridStructure(cGH const *cctkGH, int const m,
                         gh::mregs const &regsss);

void OutputGridCoordinates(cGH const *cctkGH, int const m,
                           gh::mregs const &regsss);

void OutputGridStatistics(cGH const *cctkGH);

// Functions for recomposing the grid hierarchy
void SplitRegions(cGH const *cctkGH, vector<region_t> &superregs,
                  vector<region_t> &regs);
void SplitRegions_AlongZ(cGH const *cctkGH, vector<region_t> &superregs,
                         vector<region_t> &regs);
void SplitRegions_AlongDir(cGH const *cctkGH, vector<region_t> &superregs,
                           vector<region_t> &regs, int dir);
void SplitRegions_Automatic(cGH const *cctkGH, vector<region_t> &superregs,
                            vector<region_t> &regs,
                            bvect const &no_split_dims = false);
void SplitRegionsMaps_Recursively(cGH const *cctkGH,
                                  vector<vector<region_t> > &superregss,
                                  vector<vector<region_t> > &regss);

void SplitRegionsMaps(cGH const *cctkGH, vector<vector<region_t> > &superregss,
                      vector<vector<region_t> > &regss);
void SplitRegionsMaps_Automatic(cGH const *cctkGH,
                                vector<vector<region_t> > &superregss,
                                vector<vector<region_t> > &regss,
                                bvect const &no_split_dims = false);
void SplitRegionsMaps_Recursively(cGH const *cctkGH,
                                  vector<vector<region_t> > &superregss,
                                  vector<vector<region_t> > &regss);
void SplitRegionsMaps_Balanced(cGH const *cctkGH,
                               vector<vector<region_t> > &superregss,
                               vector<vector<region_t> > &regss);

void MakeMultigridBoxes(cGH const *cctkGH, int m, gh::rregs const &regss,
                        gh::mregs &regsss);

void MakeMultigridBoxesMaps(cGH const *cctkGH, vector<gh::rregs> const &regsss,
                            vector<gh::mregs> &regssss);

// Timing statistics functions
void InitTimingStats(cGH const *cctkGH);
void BeginTimingEvolution(cGH const *cctkGH);
void StepTimingEvolution(cGH const *cctkGH);
void BeginTimingLevel(cGH const *cctkGH);
void EndTimingLevel(cGH const *cctkGH);
void BeginTimingIO(cGH const *cctkGH);
void EndTimingIO(cGH const *cctkGH, CCTK_REAL files, CCTK_REAL bytes,
                 bool is_binary);
void BeginTimingCommunication(cGH const *cctkGH);
void EndTimingCommunication(cGH const *cctkGH, CCTK_REAL messages,
                            CCTK_REAL bytes);
void UpdateTimingStats(cGH const *cctkGH);
void PrintTimingStats(cGH const *cctkGH);

} // namespace Carpet

#endif // !defined(FUNCTIONS_HH)
