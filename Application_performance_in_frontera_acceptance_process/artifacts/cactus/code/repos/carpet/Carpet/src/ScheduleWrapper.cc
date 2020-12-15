#include <cassert>
#include <list>

#include <cctk.h>
#include <cctk_Arguments.h>

#include <carpet.hh>

namespace Carpet {

using namespace std;

typedef CCTK_INT (*func)(CCTK_POINTER_TO_CONST cctkGH, CCTK_POINTER function,
                         CCTK_POINTER attribute, CCTK_POINTER data);
typedef list<func> flist;

static flist func_befores, func_afters;

extern "C" CCTK_INT Carpet_RegisterScheduleWrapper(func const func_before,
                                                   func const func_after) {
  // Add functions
  if (func_before)
    func_befores.push_back(func_before);
  if (func_after)
    func_afters.push_front(func_after);
  return 0;
}

extern "C" CCTK_INT Carpet_UnRegisterScheduleWrapper(func const func_before,
                                                     func const func_after) {
  // Remove functions
  if (func_before)
    func_befores.remove(func_before);
  if (func_after)
    func_afters.remove(func_after);
  return 0;
}

int CallBeforeRoutines(cGH const *const cctkGH, void *const function,
                       cFunctionData *const attribute, void *const data) {
  int skip = 0;
  for (flist::const_iterator fli = func_befores.begin();
       fli != func_befores.end(); ++fli) {
    skip |= (*fli)(cctkGH, function, attribute, data);
  }
  return skip;
}

int CallAfterRoutines(cGH const *const cctkGH, void *const function,
                      cFunctionData *const attribute, void *const data) {
  int res = 0;
  for (flist::const_iterator fli = func_afters.begin();
       fli != func_afters.end(); ++fli) {
    res |= (*fli)(cctkGH, function, attribute, data);
  }
  return res;
}

} // namespace Carpet
