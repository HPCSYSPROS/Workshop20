#include <cassert>

#include <cctk.h>

#include "indexing.hh"

namespace CarpetRegrid2 {

// Get indexing information for a vector grid array
void getvectorindex2(cGH const *const cctkGH, char const *const groupname,
                     int *const lsh) {
  assert(groupname);
  assert(lsh);

  int const gi = CCTK_GroupIndex(groupname);
  assert(gi >= 0);

  {
    int const ierr = CCTK_GrouplshGI(cctkGH, 1, lsh, gi);
    assert(not ierr);
  }

  cGroup groupdata;
  {
    int const ierr = CCTK_GroupData(gi, &groupdata);
    assert(not ierr);
  }
  assert(groupdata.vectorgroup);
  assert(groupdata.vectorlength >= 0);
  lsh[1] = groupdata.vectorlength;
}

} // namespace CarpetRegrid2
