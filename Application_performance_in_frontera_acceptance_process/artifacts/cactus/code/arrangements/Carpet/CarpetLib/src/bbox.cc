#include <cctk.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <cassert>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <typeinfo>

#include "defs.hh"
#include "vect.hh"

#include "bbox.hh"

// intel-17.0.0 fails with internal error 0_76 if vectorization is not disabled                 
#if __INTEL_COMPILER >= 1700
#pragma GCC optimization_level 1
#endif

using namespace std;

// Consistency checks
template <typename T, int D> void bbox<T, D>::assert_bbox_limits() const {
  ASSERT_BBOX(all(_stride > T(0)));
  ASSERT_BBOX(all((_upper - _lower) % _stride == T(0)));
  if (numeric_limits<T>::is_integer) {
    // prevent accidental wrap-around
    if (any(_lower >= numeric_limits<T>::max() / 2) or
        any(_lower <= numeric_limits<T>::min() / 2) or
        any(_upper >= numeric_limits<T>::max() / 2) or
        any(_upper <= numeric_limits<T>::min() / 2)) {
      ostringstream buf;
      T dummy;
      buf << "Tried to create a very large bbox [" << _lower << "," << _upper
          << "," << _stride << "] for the type " << typeid(dummy).name()
          << " -- it is likely that this would lead to an integer overflow";
      CCTK_WARN(CCTK_WARN_ABORT, buf.str().c_str());
    }
  }
}

// Poison
template <typename T, int D> bbox<T, D> bbox<T, D>::poison() {
  DECLARE_CCTK_PARAMETERS;

  vect<T, D> const v(deadbeef);
  return bbox(v, v, v);
}

template <typename T, int D> bool bbox<T, D>::is_poison() const {
  return D > 0 and equal_to<bbox>()(*this, poison());
}

// Accessors
template <typename T, int D> size_type bbox<T, D>::size() const {
  if (empty())
    return 0;
  const vect<T, D> sh(shape() / stride());
#ifndef CARPET_DEBUG
  return prod(vect<size_type, D>(sh));
#else
  size_type sz = 1, max = numeric_limits<size_type>::max();
  for (int d = 0; d < D; ++d) {
    if (sh[d] > max) {
      CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "size of bbox of type %s is too large -- integer overflow",
                 typeid(*this).name());
    }
    sz *= sh[d];
    max /= sh[d];
  }
  return sz;
#endif
}

// Queries

// Containment
template <typename T, int D>
bool bbox<T, D>::contains(const vect<T, D> &x) const {
  if (empty())
    return false;
  // no alignment check
  return all(x >= lower() and x <= upper());
}

template <typename T, int D>
bool bbox<T, D>::is_contained_in(const bbox &b) const {
  if (empty())
    return true;
  // no alignment check
  return all(lower() >= b.lower() and upper() <= b.upper());
}

// Intersection
template <typename T, int D> bool bbox<T, D>::intersects(const bbox &b) const {
  if (empty())
    return false;
  if (b.empty())
    return false;
  // no alignment check
  return all(upper() >= b.lower() and lower() <= b.upper());
}

// Alignment check
template <typename T, int D>
bool bbox<T, D>::is_aligned_with(const bbox &b) const {
  return all(stride() == b.stride() and
             (lower() - b.lower()) % stride() == T(0));
}

// Operators
template <typename T, int D> bool bbox<T, D>::operator==(const bbox &b) const {
  if (empty() and b.empty())
    return true;
  ASSERT_BBOX(all(stride() == b.stride()));
  return all(lower() == b.lower() and upper() == b.upper());
}

template <typename T, int D> bool bbox<T, D>::operator!=(const bbox &b) const {
  return not(*this == b);
}

#if 0
// Introduce an ordering on bboxes
template<typename T, int D>
bool bbox<T,D>::operator< (const bbox& b) const {
  // An arbitraty order: empty boxes come first, then sorted by lower
  // bound, then by upper bound, then by coarseness
  if (b.empty()) return false;
  if (empty()) return true;
  for (int d=D-1; d>=0; --d) {
    if (lower()[d] < b.lower()[d]) return true;
    if (lower()[d] > b.lower()[d]) return false;
  }
  for (int d=D-1; d>=0; --d) {
    if (upper()[d] < b.upper()[d]) return true;
    if (upper()[d] > b.upper()[d]) return false;
  }
  for (int d=D-1; d>=0; --d) {
    if (stride()[d] > b.stride()[d]) return true;
    if (stride()[d] < b.stride()[d]) return false;
  }
  return false;
}
#endif

template <typename T, int D> bool bbox<T, D>::operator<=(const bbox &b) const {
  return is_contained_in(b);
}

template <typename T, int D> bool bbox<T, D>::operator>=(const bbox &b) const {
  return b <= *this;
}

template <typename T, int D> bool bbox<T, D>::operator<(const bbox &b) const {
  return *this <= b and *this != b;
}

template <typename T, int D> bool bbox<T, D>::operator>(const bbox &b) const {
  return b < *this;
}

// Expand the bbox a little by multiples of the stride
template <typename T, int D>
bbox<T, D> bbox<T, D>::expand(const vect<T, D> &lo,
                              const vect<T, D> &hi) const {
  // Allow expansion only into directions where the extent is not negative
  // ASSERT_BBOX (all(lower()<=upper() or (lo==T(0) and hi==T(0))));
  ASSERT_BBOX(all(shape() >= vect<T, D>(0) or (lo == T(0) and hi == T(0))));
  const vect<T, D> str = stride();
  const vect<T, D> lb = lower() - lo * str;
  const vect<T, D> ub = upper() + hi * str;
  return bbox(lb, ub, str);
}

// Expand the bbox a little by multiples of a fraction of the stride
template <typename T, int D>
bbox<T, D> bbox<T, D>::expand(const vect<T, D> &lo, const vect<T, D> &hi,
                              const vect<T, D> &denom) const {
  // Allow expansion only into directions where the extent is not negative
  // ASSERT_BBOX (all(lower()<=upper() or (lo==T(0) and hi==T(0))));
  ASSERT_BBOX(all(shape() >= vect<T, D>(0) or (lo == T(0) and hi == T(0))));
  ASSERT_BBOX(all(denom > vect<T, D>(0)));
  const vect<T, D> str = stride();
  ASSERT_BBOX(all(str % denom == vect<T, D>(0) or (lo == T(0) and hi == T(0))));
  const vect<T, D> lb = lower() - lo * str / denom;
  const vect<T, D> ub = upper() + hi * str / denom;
  return bbox(lb, ub, str);
}

// Find the smallest b-compatible box around *this
template <typename T, int D>
bbox<T, D> bbox<T, D>::expanded_for(const bbox &b) const {
  if (empty())
    return bbox(b.lower(), b.lower() - b.stride(), b.stride());
  const vect<T, D> str = b.stride();
  const vect<T, D> loff = imod(lower() - b.lower(), str);
  const vect<T, D> uoff = imod(upper() - b.lower(), str);
  const vect<T, D> lo = lower() - loff; // go outwards
  const vect<T, D> up = upper() + (str - uoff) % str;
  return bbox(lo, up, str);
}

// Find the largest b-compatible box inside *this
template <typename T, int D>
bbox<T, D> bbox<T, D>::contracted_for(const bbox &b) const {
  if (empty())
    return bbox(b.lower(), b.lower() - b.stride(), b.stride());
  const vect<T, D> str = b.stride();
  const vect<T, D> loff = imod(lower() - b.lower(), str);
  const vect<T, D> uoff = imod(upper() - b.lower(), str);
  const vect<T, D> lo = lower() + (str - loff) % str; // go inwards
  const vect<T, D> up = upper() - uoff;
  return bbox(lo, up, str);
}

// Find the smallest open b-compatible box around *this:
// This routine is similar to expanded_for. However, it returns a box
// that is possibly larger than the box returned by expanded_for.
template <typename T, int D>
bbox<T, D> bbox<T, D>::anti_contracted_for(const bbox &b) const {
  if (empty())
    return bbox(b.lower(), b.lower() - b.stride(), b.stride());
  return expand(1, 1).expanded_for(b).expand(-1, -1);
  // if (all(stride() <= b.stride())) {
  //   // New stride is larger or equal to current stride: forward call
  //   // to expanded_for
  //   return expanded_for(b);
  // }
  // if (all(stride() > b.stride())) {
  //   // New stride is smaller than current stride: expand box by one
  //   // stride, call expanded_for, and shrink box by on stride again
  //   return expand(1, 1).expanded_for(b).expand(-1, -1);
  // }
  // CCTK_BUILTIN_UNREACHABLE();
}

// Smallest bbox containing both boxes
template <typename T, int D>
bbox<T, D> bbox<T, D>::expanded_containing(const bbox &b) const {
  if (empty())
    return b;
  if (b.empty())
    return *this;
  ASSERT_BBOX(is_aligned_with(b));
  const vect<T, D> lo = min(lower(), b.lower());
  const vect<T, D> up = max(upper(), b.upper());
  const vect<T, D> str = min(stride(), b.stride());
  return bbox(lo, up, str);
}

// Iterators
template <typename T, int D>
bbox<T, D>::iterator::iterator(const bbox &box_, const vect<T, D> &pos_)
    : box(box_), pos(pos_) {
  if (box.empty())
    pos = box.upper();
}

template <typename T, int D>
bool bbox<T, D>::iterator::operator!=(const iterator &i) const {
  return any(pos != i.pos);
}

template <typename T, int D>
typename bbox<T, D>::iterator &bbox<T, D>::iterator::operator++() {
  for (int d = 0; d < D; ++d) {
    pos[d] += box.stride()[d];
    if (pos[d] <= box.upper()[d])
      break;
    pos[d] = box.lower()[d];
  }
  return *this;
}

template <typename T, int D>
typename bbox<T, D>::iterator bbox<T, D>::begin() const {
  return iterator(*this, lower());
}

template <typename T, int D>
typename bbox<T, D>::iterator bbox<T, D>::end() const {
  return iterator(*this, lower());
}

// Input
template <typename T, int D> void bbox<T, D>::input(istream &is) {
  try {
    skipws(is);
    consume(is, '(');
    is >> _lower;
    skipws(is);
    consume(is, ':');
    is >> _upper;
    skipws(is);
    consume(is, ':');
    is >> _stride;
    skipws(is);
    if (is.peek() == '/') {
      consume(is, '/');
      vect<T, D> lower_dummy;
      is >> lower_dummy;
      skipws(is);
      consume(is, ':');
      vect<T, D> upper_dummy;
      is >> upper_dummy;
      skipws(is);
      consume(is, '/');
      vect<T, D> shape_dummy;
      is >> shape_dummy;
      skipws(is);
      consume(is, '/');
      size_type size_dummy;
      is >> size_dummy;
      ASSERT_BBOX(is.good());
      skipws(is);
    }
    consume(is, ')');
  } catch (input_error &err) {
    T Tdummy;
    cout << "Input error while reading a bbox<" << typestring(Tdummy) << ","
         << D << ">" << endl;
    throw err;
  }
  if (any(_stride <= T(0))) {
    cout << "While reading the bbox " << *this << ":" << endl
         << "   The stride is not positive." << endl;
    throw input_error();
  }
  if (any((_upper - _lower) % _stride != T(0))) {
    cout << "While reading the bbox " << *this << ":" << endl
         << "   The stride does not evenly divide the extent." << endl;
    throw input_error();
  }
  ASSERT_BBOX(all(_stride > T(0)));
  ASSERT_BBOX(all((_upper - _lower) % _stride == T(0)));
}

// Output
template <typename T, int D> void bbox<T, D>::output(ostream &os) const {
  os << "(" << lower() << ":" << upper() << ":" << stride() << "/"
     << idiv(lower(), stride()) << ":" << idiv(upper(), stride()) << "/"
     << shape() / stride() << "/" << size() << ")";
}

// Comparison

namespace std {
// ==
template <typename T, int D>
bool equal_to<bbox<T, D> >::operator()(const bbox<T, D> &x,
                                       const bbox<T, D> &y) const {
  if (x.empty() and y.empty())
    return true;
  if (x.empty() or y.empty())
    return false;
  /*const*/ equal_to<vect<T, D> > vect_equal_to;
  if (not vect_equal_to(x.stride(), y.stride()))
    return false;
  if (not vect_equal_to(x.lower(), y.lower()))
    return false;
  if (not vect_equal_to(x.upper(), y.upper()))
    return false;
  return true;
}

// <
template <typename T, int D>
bool less<bbox<T, D> >::operator()(const bbox<T, D> &x,
                                   const bbox<T, D> &y) const {
  // Empty bboxes compare less than any non-empty bbox
  if (y.empty())
    return false;
  if (x.empty())
    return true;
  /*const*/ less<vect<T, D> > vect_less;
  if (vect_less(x.stride(), y.stride()))
    return true;
  if (vect_less(y.stride(), x.stride()))
    return false;
  if (vect_less(x.lower(), y.lower()))
    return true;
  if (vect_less(y.lower(), x.lower()))
    return false;
  if (vect_less(x.upper(), y.upper()))
    return true;
  if (vect_less(y.upper(), x.upper()))
    return false;
  return false;
}
}

// Note: We need all dimensions all the time.
template class bbox<int, 0>;
template class bbox<int, 1>;
template class bbox<int, 2>;
template class bbox<int, 3>;
template class bbox<int, 4>;
template class bbox<CCTK_REAL, dim>;

namespace std {
template struct less<bbox<int, 0> >;
template struct less<bbox<int, 1> >;
template struct less<bbox<int, 2> >;
template struct less<bbox<int, 3> >;
template struct less<bbox<int, 4> >;
}
