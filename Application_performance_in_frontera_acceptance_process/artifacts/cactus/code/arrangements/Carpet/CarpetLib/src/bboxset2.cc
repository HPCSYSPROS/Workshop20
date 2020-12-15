#include <cctk.h>

#ifdef CARPET_ENABLE_BBOXSET2

#include "bboxset2.hh"

namespace bboxset2 {

template class bboxset<int, 0>;

template class bboxset<int, 1>;
template void bboxset<int, 1>::serialise(set<bbox> &) const;

template class bboxset<int, 2>;
template void bboxset<int, 2>::serialise(set<bbox> &) const;

template class bboxset<int, 3>;
template void bboxset<int, 3>::serialise(set<bbox> &) const;
template void bboxset<int, 3>::serialise(vector<bbox> &) const;

} // namespace bboxset2

#include "dh.hh"
#include "region.hh"

namespace bboxset2 {

template bboxset<int, 3>::bboxset(const vector<dh::full_dboxes> &,
                                  const bbox dh::full_dboxes::*);
template bboxset<int, 3>::bboxset(const vector<dh::full_dboxes> &,
                                  const bboxset dh::full_dboxes::*);
template bboxset<int, 3>::bboxset(const vector<region_t> &,
                                  const bbox region_t::*);

} // namespace bboxset2

#endif // #ifdef CARPET_ENABLE_BBOXSET2
