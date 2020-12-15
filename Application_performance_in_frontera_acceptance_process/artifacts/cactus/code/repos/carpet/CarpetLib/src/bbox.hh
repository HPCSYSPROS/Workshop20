#ifndef BBOX_HH
#define BBOX_HH

#include <cassert>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <limits>

#include "defs.hh"
#include "vect.hh"

using namespace std;

#ifdef CARPET_DEBUG
#define ASSERT_BBOX(x) assert(x)
#else
#define ASSERT_BBOX(x)
#endif

// Forward declaration
template <typename T, int D> class bbox;

// Input/Output
template <typename T, int D> istream &operator>>(istream &is, bbox<T, D> &b);
template <typename T, int D>
ostream &operator<<(ostream &os, const bbox<T, D> &b);

/**
 * A bounding box, i.e., a rectangle with lower and upper bound and a
 * stride.
 */
template <typename T, int D> class bbox {

  // Fields

  /** Bounding box bounds and stride.  The bounds are inclusive.  */
  vect<T, D> _lower, _upper, _stride;

  // Consistency checks

  void assert_bbox_limits() const CCTK_MEMBER_ATTRIBUTE_PURE;

public:
  // Constructors

  /** Construct an empty bbox.  */
  bbox() : _lower(T(1)), _upper(T(0)), _stride(T(1)) {}

  /** Copy constructor.  */
  bbox(const bbox &b)
      : _lower(b._lower), _upper(b._upper), _stride(b._stride) {}

  /** Assignment operator.  */
  bbox &operator=(const bbox &b) {
    _lower = b._lower;
    _upper = b._upper;
    _stride = b._stride;
    return *this;
  }

  /** Create a bbox from bounds and stride.  */
  bbox(const vect<T, D> &lower_, const vect<T, D> &upper_,
       const vect<T, D> &stride_)
      : _lower(lower_), _upper(upper_), _stride(stride_) {
#ifdef CARPET_DEBUG
    assert_bbox_limits();
#endif
  }

  // Poison
  static bbox poison() CCTK_ATTRIBUTE_PURE;

  bool is_poison() const CCTK_MEMBER_ATTRIBUTE_PURE;

  // Accessors
  // (Don't return references; *this might be a temporary)

  /** Get lower bound.  */
  vect<T, D> lower() const { return _lower; }

  /** Get upper bound.  */
  vect<T, D> upper() const { return _upper; }

  /** Get bounds.  */
  vect<vect<T, D>, 2> bounds() const {
    return vect<vect<T, D>, 2>(_lower, _upper);
  }

  /** Get stride.  */
  vect<T, D> stride() const { return _stride; }

  /** Get offset.  */
  vect<T, D> offset() const { return imod(_lower, _stride); }

  /** Get the shape (or extent).  */
  vect<T, D> shape() const { return _upper - _lower + _stride; }

  /** Determine whether the bbox is empty.  */
  bool empty() const { return any(lower() > upper()); }

  /** Return the size, which is the number of contained points.  */
  // T size () const;
  size_type size() const CCTK_MEMBER_ATTRIBUTE_PURE;

  // Queries

  /** Find out whether the bbox contains the point x.  */
  bool contains(const vect<T, D> &x) const CCTK_MEMBER_ATTRIBUTE_PURE;

  /** Find out whether this bbox is contained in the bbox b.  */
  bool is_contained_in(const bbox &b) const CCTK_MEMBER_ATTRIBUTE_PURE;

  /** Find out whether this bbox intersects the bbox b.  */
  bool intersects(const bbox &b) const CCTK_MEMBER_ATTRIBUTE_PURE;

  /** Find out whether this bbox is aligned with the bbox b.
      ("aligned" means that both bboxes have the same stride and that
      their boundaries are commensurate.)  */
  bool is_aligned_with(const bbox &b) const CCTK_MEMBER_ATTRIBUTE_PURE;

  // Operators
  bool operator==(const bbox &b) const CCTK_MEMBER_ATTRIBUTE_PURE;
  bool operator!=(const bbox &b) const CCTK_MEMBER_ATTRIBUTE_PURE;
  bool operator<=(const bbox &b) const CCTK_MEMBER_ATTRIBUTE_PURE;
  bool operator>=(const bbox &b) const CCTK_MEMBER_ATTRIBUTE_PURE;
  bool operator<(const bbox &b) const CCTK_MEMBER_ATTRIBUTE_PURE;
  bool operator>(const bbox &b) const CCTK_MEMBER_ATTRIBUTE_PURE;

  /** Calculate the intersection (the set of common points) with the
      bbox b.  */
  bbox operator&(const bbox &b) const {
    ASSERT_BBOX(all(stride() == b.stride()));
    ASSERT_BBOX(is_aligned_with(b));
    vect<T, D> lo = max(lower(), b.lower());
    vect<T, D> up = min(upper(), b.upper());
    return bbox(lo, up, stride());
  }

  /** Expand (enlarge) the bbox by multiples of the stride.  */
  bbox expand(const vect<T, D> &lo,
              const vect<T, D> &hi) const CCTK_MEMBER_ATTRIBUTE_PURE;
  bbox expand(const vect<vect<T, D>, 2> &lohi) const {
    return expand(lohi[0], lohi[1]);
  }

  /** Shift the bbox by multiples of the stride.  */
  bbox shift(const vect<T, D> &v) const { return expand(-v, v); }

  /** Expand (enlarge) the bbox by multiples of a fraction of the
      stride.  */
  bbox expand(const vect<T, D> &lo, const vect<T, D> &hi,
              const vect<T, D> &denom) const CCTK_MEMBER_ATTRIBUTE_PURE;
  bbox expand(const vect<vect<T, D>, 2> &lohi, const vect<T, D> &denom) const {
    return expand(lohi[0], lohi[1], denom);
  }

  /** Shift the bbox by multiples of a fraction of the stride.  */
  bbox shift(const vect<T, D> &v, const vect<T, D> &denom) const {
    return expand(-v, v, denom);
  }

  /** Find the smallest b-compatible box around this bbox.
      ("compatible" means having the same stride.)  */
  bbox expanded_for(const bbox &b) const CCTK_MEMBER_ATTRIBUTE_PURE;

  /** Find the largest b-compatible box inside this bbox.  */
  bbox contracted_for(const bbox &b) const CCTK_MEMBER_ATTRIBUTE_PURE;

  /** Find the smallest open b-compatible box around *this */
  bbox anti_contracted_for(const bbox &b) const CCTK_MEMBER_ATTRIBUTE_PURE;

  /** Find the smallest bbox containing both boxes.  */
  bbox
  expanded_containing(const bbox<T, D> &b) const CCTK_MEMBER_ATTRIBUTE_PURE;

  // Iterators

  /** An iterator over all points in a bbox.  */
  class iterator {
  protected:
    /** The bbox over which we iterate.  */
    const bbox &box;
    /** Current position.  */
    vect<T, D> pos;

  public:
    /** Constructor.  */
    iterator(const bbox &box, const vect<T, D> &pos);
    /** Accessor.  */
    const vect<T, D> &operator*() const { return pos; }
    /** Check whether the position is the same.  */
    bool operator!=(const iterator &i) const;
    /** Advance.  */
    iterator &operator++();
  };

  /** Create an iterator that points to the first point in a bbox.  */
  iterator begin() const CCTK_MEMBER_ATTRIBUTE_PURE;
  /** Create an iterator that points "after the last point" in a bbox,
      which means that it also points to the first point.  */
  iterator end() const CCTK_MEMBER_ATTRIBUTE_PURE;

  size_type index(const vect<T, D> &pos) const {
    assert(not empty());
    bbox const posbox(pos, pos, stride());
    assert(is_aligned_with(posbox));
    assert(contains(pos));
    size_type const i = ::index(shape(), (pos - lower()) / stride());
    assert(i >= 0 and i < size());
    return i;
  }

  // Memory usage
  size_t memory() const {
    return memoryof(_lower) + memoryof(_upper) + memoryof(_stride);
  }

  // Input/Output helpers
  void input(istream &is);
  void output(ostream &os) const;
};

// Memory usage

template <typename T, int D> inline size_t memoryof(bbox<T, D> const &b) {
  return b.memory();
}

// Input

/** Read a formatted bbox from a stream.  */
template <typename T, int D>
inline istream &operator>>(istream &is, bbox<T, D> &b) {
  b.input(is);
  return is;
}

// Output

/** Write a bbox formatted to a stream.  */
template <typename T, int D>
inline ostream &operator<<(ostream &os, const bbox<T, D> &b) {
  b.output(os);
  return os;
}

// Comparison

namespace std {
// ==
template <typename T, int D>
struct equal_to<bbox<T, D> > : binary_function<bbox<T, D>, bbox<T, D>, bool> {
  bool operator()(const bbox<T, D> &x, const bbox<T, D> &y) const;
};

// <
template <typename T, int D>
struct less<bbox<T, D> > : binary_function<bbox<T, D>, bbox<T, D>, bool> {
  bool operator()(const bbox<T, D> &x, const bbox<T, D> &y) const;
};

// >
template <typename T, int D>
struct greater<bbox<T, D> > : binary_function<bbox<T, D>, bbox<T, D>, bool> {
  bool operator()(const bbox<T, D> &x, const bbox<T, D> &y) const {
    return less<bbox<T, D> >()(y, x);
  }
};

// >=
template <typename T, int D>
struct greater_equal<bbox<T, D> >
    : binary_function<bbox<T, D>, bbox<T, D>, bool> {
  bool operator()(const bbox<T, D> &x, const bbox<T, D> &y) const {
    return not less<bbox<T, D> >()(x, y);
  }
};

// <=
template <typename T, int D>
struct less_equal<bbox<T, D> > : binary_function<bbox<T, D>, bbox<T, D>, bool> {
  bool operator()(const bbox<T, D> &x, const bbox<T, D> &y) const {
    return not greater<bbox<T, D> >()(x, y);
  }
};

// !=
template <typename T, int D>
struct not_equal_to<bbox<T, D> >
    : binary_function<bbox<T, D>, bbox<T, D>, bool> {
  bool operator()(const bbox<T, D> &x, const bbox<T, D> &y) const {
    return not equal_to<bbox<T, D> >()(x, y);
  }
};
}

#endif // BBOX_HH
