#ifndef CARPET_HH
#define CARPET_HH

#include <vector>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Functions.h>
#include <cctk_Schedule.h>

#include <gh.hh>

#include "carpet_public.hh"

namespace Carpet {

using namespace std;

// Scheduled functions
extern "C" {
int CarpetStartup(void);
int CarpetMultiModelStartup(void);
void CarpetParamCheck(CCTK_ARGUMENTS);
void CarpetRefineTimeStep(CCTK_ARGUMENTS);
void CarpetUnusedMask(CCTK_ARGUMENTS);
}

// Registered functions
void *SetupGH(tFleshConfig *fc, int convLevel, cGH *cgh);

int Initialise(tFleshConfig *config);
int Evolve(tFleshConfig *config);
int Shutdown(tFleshConfig *config);
int OutputGH(const cGH *cgh);

int CallFunction(void *function, cFunctionData *attribute, void *data);

// Other functions
bool Regrid(cGH const *cctkGH, bool force_recompose, bool do_init);

void CycleTimeLevels(cGH *cgh);
void UncycleTimeLevels(cGH *cgh);
void FlipTimeLevels(cGH *cgh);
void FillTimeLevels(const cGH *cgh);
void SyncGroups(const cGH *cgh, const vector<int> &groups);
int SyncProlongateGroups(const cGH *cgh, const vector<int> &groups,
                         cFunctionData const *function_data = NULL);

// Sanity checks
enum checktimes {
  currenttime,
  currenttimebutnotifonly,
  previoustime,
  allbutlasttime,
  allbutcurrenttime,
  alltimes
};

int min_timelevel(checktimes where, int num_tl,
                  bool persistent) CCTK_ATTRIBUTE_CONST;
int max_timelevel(checktimes where, int num_tl,
                  bool persistent) CCTK_ATTRIBUTE_CONST;

void Poison(const cGH *cgh, checktimes where, int what = 0);
void PoisonGroup(const cGH *cgh, int group, checktimes where);
void PoisonCheck(const cGH *cgh, checktimes where);

void CalculateChecksums(const cGH *cgh, checktimes where);
void CheckChecksums(const cGH *cgh, checktimes where);

// Schedule
int CallBeforeRoutines(cGH const *cctkGH, void *const function,
                       cFunctionData *const attribute, void *const data);
int CallAfterRoutines(cGH const *cctkGH, void *const function,
                      cFunctionData *const attribute, void *const data);

// Debugging output
void Output(const char *fmt, ...);
void Waypoint(const char *fmt, ...);
void Checkpoint(const char *fmt, ...);

// Error output
void UnsupportedVarType(int vindex);

// Check for a map0group
bool IsMap0Group(int gindex);

} // namespace Carpet

#endif // !defined(CARPET_HH)
