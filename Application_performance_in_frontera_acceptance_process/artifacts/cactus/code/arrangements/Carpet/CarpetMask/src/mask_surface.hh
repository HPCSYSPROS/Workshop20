#include <cctk.h>
#include <cctk_Arguments.h>

namespace CarpetMask {

extern "C" {
void CarpetSurfaceParamCheck(CCTK_ARGUMENTS);
void CarpetSurfaceSetup(CCTK_ARGUMENTS);
}

} // namespace CarpetMask
