#ifndef BALANCE_HH
#define BALANCE_HH

#include <cctk.h>

#include <vector>

using namespace std;

namespace CarpetLib {

// This routine splits N items over P workers.  It can split the
// items if necessary to ensure a maximum imbalance, and it can
// ensure (artificially) that each worker receives the same number
// of items.

template <typename T>
void balance(vector<T> const &items_, vector<vector<T> > &split_items_,
             int nworkers, CCTK_REAL max_imbalance, bool ensure_same_size);

} // namespace CarpetLib

#endif // #ifndef BALANCE_HH
