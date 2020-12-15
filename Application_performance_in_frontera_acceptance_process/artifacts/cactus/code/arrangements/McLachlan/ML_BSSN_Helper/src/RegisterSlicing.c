#include <cctk.h>
#include <Slicing.h>

int ML_BSSN_RegisterSlicing(void) {
  Einstein_RegisterSlicing("ML_BSSN");
  return 0;
}
