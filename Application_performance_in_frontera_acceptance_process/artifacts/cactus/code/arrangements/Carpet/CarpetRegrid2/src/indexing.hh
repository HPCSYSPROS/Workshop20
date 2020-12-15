#ifndef INDEXING_HH
#define INDEXING_HH

#include <cassert>

#include <cctk.h>

namespace CarpetRegrid2 {

// Get indexing information for a vector grid array
void getvectorindex2(cGH const *cctkGH, char const *groupname, int *lsh);

static inline int index2(int const *const lsh, int const i, int const j) {
  assert(lsh);
  assert(lsh[0] >= 0);
  assert(lsh[1] >= 0);
  assert(i >= 0 and i < lsh[0]);
  assert(j >= 0 and j < lsh[1]);
  return i + lsh[0] * j;
}

} // namespace CarpetRegrid2

#endif // #ifndef INDEXING_HH
