#include <cctk.h>
#include <cctk_Arguments.h>

namespace CarpetMask {

extern "C" {
void CarpetExcludedSetup(CCTK_ARGUMENTS);
}

} // namespace CarpetMask
