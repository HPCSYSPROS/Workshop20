#include <algorithm>
#include <cassert>
#include <iostream>
#include <limits>
#include <set>
#include <stack>
#include <vector>

#include "defs.hh"

#include "bboxset1.hh"

using namespace std;

namespace bboxset1 {

template class bboxset<int, 1>;
template void bboxset<int, 1>::serialise(set<bbox<int, 1> > &s) const;
template void bboxset<int, 1>::serialise(vector<bbox<int, 1> > &v) const;
template size_t memoryof(const bboxset<int, 1> &s);
template istream &operator>>(istream &is, bboxset<int, 1> &s);
template ostream &operator<<(ostream &os, const bboxset<int, 1> &s);

template class bboxset<int, 2>;
template void bboxset<int, 2>::serialise(set<bbox<int, 2> > &s) const;
template void bboxset<int, 2>::serialise(vector<bbox<int, 2> > &v) const;
template size_t memoryof(const bboxset<int, 2> &s);
template istream &operator>>(istream &is, bboxset<int, 2> &s);
template ostream &operator<<(ostream &os, const bboxset<int, 2> &s);

template class bboxset<int, 3>;
template void bboxset<int, 3>::serialise(set<bbox<int, 3> > &s) const;
template void bboxset<int, 3>::serialise(vector<bbox<int, 3> > &f) const;
template size_t memoryof(const bboxset<int, 3> &s);
template istream &operator>>(istream &is, bboxset<int, 3> &s);
template ostream &operator<<(ostream &os, const bboxset<int, 3> &s);

} // namespace bboxset1

#include "dh.hh"
#include "region.hh"

namespace bboxset1 {

template bboxset<int, dim>::bboxset(
    const vector<dh::full_dboxes> &vb,
    const bbox<int, dim> dh::full_dboxes::*const v);
template bboxset<int, dim>::bboxset(
    const vector<dh::full_dboxes> &vb,
    const bboxset<int, dim> dh::full_dboxes::*const v);
template bboxset<int, dim>::bboxset(const vector<region_t> &vb,
                                    const bbox<int, dim> region_t::*const v);

} // namespace bboxset1
