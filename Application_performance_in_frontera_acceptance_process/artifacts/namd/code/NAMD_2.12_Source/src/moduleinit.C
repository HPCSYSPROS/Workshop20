#include "charm++.h"
extern void _registerNeighborLB(void);
void _registerExternalModules(void) {
  _registerNeighborLB();
}
