#include <cctk.h>
#include <stdlib.h>

#define XSTR(x) #x
#define STR(x) XSTR(x)

int hwloc_lstopo(void) {
  if (CCTK_MyProc(NULL) == 0) {
    CCTK_INFO("Hardware locality information:");
    CCTK_INFO(STR(HWLOC_UTILPATH) "/lstopo -v");
    system(STR(HWLOC_UTILPATH) "/lstopo -v");
  }
  return 0;
}
