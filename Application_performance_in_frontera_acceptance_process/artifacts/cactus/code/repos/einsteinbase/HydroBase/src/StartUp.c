#include <cctk.h>

int HydroBase_StartUp(void) {
  CCTK_RegisterBanner("HydroBase: Let it flow.");

  return 0;
}
