#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <sstream>

#include <cctk.h>
#include <cctk_Parameters.h>

#include <Timer.hh>

#include <dist.hh>

#include <carpet.hh>

namespace Carpet {

using namespace std;

int OutputGH(cGH const *const cctkGH) {
  DECLARE_CCTK_PARAMETERS;

  static Timers::Timer timer("OutputGH");
  timer.start();

  Checkpoint("OutputGH");

  int const num_methods = CCTK_NumIOMethods();
  if (num_methods == 0) {
    timer.stop();
    return -1;
  }

  static vector<Timers::Timer *> timers;
  timers.resize(num_methods, NULL);

  int num_vars = 0;
  for (int handle = 0; handle < num_methods; ++handle) {

    IOMethod const *const method = CCTK_IOMethod(handle);
    assert(method);

    if (not timers.AT(handle)) {
      ostringstream buf;
      buf << method->implementation << "::" << method->name << " [" << handle
          << "]";
      timers.AT(handle) = new Timers::Timer(buf.str().c_str());
    }

    timers.AT(handle)->start();
    num_vars += method->OutputGH(cctkGH);
    timers.AT(handle)->stop();

  } // for handle

  timer.stop();

  return num_vars;
}

} // namespace Carpet
