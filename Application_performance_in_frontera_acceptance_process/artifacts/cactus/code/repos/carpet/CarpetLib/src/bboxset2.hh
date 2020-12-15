#ifndef BBOXSET2_HH
#define BBOXSET2_HH

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <utility>
#include <vector>

#include <cctk.h>

#include "bbox.hh"
#include "defs.hh"
#include "vect.hh"

using namespace std;

/* A note on CARPET_AVOID_LAMBDA:

   When this macro is defined, then we cannot use C++11 lambdas.
   Instead, we use macros, both for the lambdas that we would
   otherwise create, and for the functions to which the lambdas would
   otherwise be passed. The functions taking lambdas as arguments
   could also be re-written as iterators, but we don't do this because
   this would be significantly more complex. Changing these functions
   and the lambdas to macros is a straightforward syntactical
   replacement.

   The macros can behave as expressions or as statements. Macros that
   behave like an expression can easily be created via the ({ ... })
   GNU extension, which is supported by almost all compilers.
   Unfortunately, e.g. the Blue Gene/Q IBM compiler does not seem to
   handle these correctly when they contain variable declarations
   requiring a non-trivial constructor. Therefore, we use macros that
   behave like statements.

   The functions that would otherwise take a lambda as argument are
   turned into a macro in the standard manner:
   - Name becomes all upper case
   - Name may have a suffix added since macros cannot be overloaded
   - Macro is wrapped in do { ... } while (0)
   - Macro arguments are safely accessed by defining local variables,
     and assigning the macro argument to them. Note that this is
     slightly non-trivial: the macro argument name must not be used
     elsewere (since it would otherwise be replaced as well), so we
     append an underscore to the original function argument name. If
     the macro argument that is passed has the name name as the local
     variable in which it is captured, then there is a conflict, so we
     choose a local variable name that has two underscores appended.
     Since this looks ugly, we then define a second local variable
     (without underscores) and assign in the double-underscore local
     variable.
   - Calling a lambda is also non-trivial, since the respective macro
     argument will be a sequence of statements, not a function-like
     object to which arguments can be passed. Instead, we store the
     arguments that the lambda would receive in local variables
     numbered _1, _2, ..., and define (if necessary) a local variable
     _0 that will hold the result. To keep our namespace clean, we put
     the "lambda call" into curly braces.

   The lambda definitions are then always located inside a macro call
   as just described. The actual lambda definition becomes a sequence
   of statements that is one argument to the macro call. (Note that
   one cannot have top-level commas in macro calls, since these would
   separate macro arguments. Luckily, top-level commas do not often
   appear in a sequence of statements.) Since the arguments will be
   passed in variables called _1, _2, ..., it is convenient to define
   local variables and assign these variables to them. Note also that
   the return value needs to be assigned to a variable called _0.

   Alternative approaches (that did not work) include:
   - Define local classes that can be called like a function. This
     does not work because these are local classes, and hence cannot
     be used as template argument.
   - Define local classes that can be called like a function, and make
     them a sub-class of an abstract base class. Instead of using
     templates, change the functions that accept lambda arguments from
     templates to using virtual function calls to evaluate the lambda.
     This does not work because ... (why did this not work?)
   - Use Boost lambdas. This does not work because these are intended
     for short expressions, and do not work for longer sections of
     code.

   Alternative approaches that could work:
   - Change the functions that accept lambda arguments into C++
     iterators. This would work, but will probably be d tedious.
   - Change the lambdas into function-like classes, but define these
     classes outside their containing function. This would work, but
     would be an invasive change to the code.
 */

#ifdef CARPET_ENABLE_BBOXSET2

namespace bboxset2 {

template <typename T, int D> class bboxset {
  template <typename, int> friend class bboxset;
  typedef ::vect<T, D> vect;
  typedef ::bbox<T, D> bbox;
  typedef ::vect<T, D - 1> vect1;
  typedef ::bbox<T, D - 1> bbox1;
  typedef bboxset<T, D - 1> bboxset1;

#if 0
  // We can't use auto_ptr or unique_ptr because we make read-only
  // copies of pointers. We could use weak_ptr for these copies, but
  // we don't do that yet.
  template<typename U>
  using smart_ptr = auto_ptr<U>;
#elif 0
  template <typename U> using smart_ptr = unique_ptr<U>;
#else
// This "using" declaration is not accepted by older C++ compilers
// template<typename U>
// using smart_ptr = shared_ptr<U>;
#define smart_ptr shared_ptr
#endif

  typedef map<T, smart_ptr<bboxset1> > subsets_t;
  subsets_t subsets;

  vect stride, offset;

  bool is_poison_;

  template <typename F>
  void traverse_subsets(const F &f) const
  // void traverse_subsets(function<void(int, const bboxset1&)>& f) const
  {
    bboxset1 decoded_subset;
    assert(decoded_subset.empty());
    forall(const auto &pos_subset, subsets) {
      const T &pos = pos_subset.first;
      const bboxset1 *const subsetp = pos_subset.second.get();
      decoded_subset ^= *subsetp;
      f(pos, decoded_subset);
    }
    assert(decoded_subset.empty());
  }
#define TRAVERSE_SUBSETS1(f)                                                   \
  do {                                                                         \
    bboxset1 decoded_subset;                                                   \
    assert(decoded_subset.empty());                                            \
    forall(const auto &pos_subset, subsets) {                                  \
      const T &pos = pos_subset.first;                                         \
      const bboxset1 *const subsetp = pos_subset.second.get();                 \
      decoded_subset ^= *subsetp;                                              \
      {                                                                        \
        const auto &_1(pos);                                                   \
        const auto &_2(decoded_subset);                                        \
        f;                                                                     \
      }                                                                        \
    }                                                                          \
    assert(decoded_subset.empty());                                            \
  } while (0)

  template <typename F>
  void traverse_subsets(const F &f, const bboxset &other) const {
    bboxset1 decoded_subset0;
    bboxset1 decoded_subset1;
    assert(decoded_subset0.empty());
    assert(decoded_subset1.empty());

    typedef typename subsets_t::const_iterator subsets_iter_t;
    subsets_iter_t iter0 = subsets.begin();
    subsets_iter_t iter1 = other.subsets.begin();
    subsets_iter_t const end0 = subsets.end();
    subsets_iter_t const end1 = other.subsets.end();
    while (iter0 != end0 or iter1 != end1) {
      const T next_pos0 =
          iter0 != end0 ? iter0->first : numeric_limits<T>::max();
      const T next_pos1 =
          iter1 != end1 ? iter1->first : numeric_limits<T>::max();
      const T pos = min(next_pos0, next_pos1);
      const bool active0 = next_pos0 == pos;
      const bool active1 = next_pos1 == pos;
      const bboxset1 *const subset0p = active0 ? iter0->second.get() : 0;
      const bboxset1 *const subset1p = active1 ? iter1->second.get() : 0;
      if (active0)
        decoded_subset0 ^= *subset0p;
      if (active1)
        decoded_subset1 ^= *subset1p;

      f(pos, decoded_subset0, decoded_subset1);

      if (active0)
        ++iter0;
      if (active1)
        ++iter1;
    }
    assert(decoded_subset0.empty());
    assert(decoded_subset1.empty());
  }
#define TRAVERSE_SUBSETS2(f, other_)                                           \
  do {                                                                         \
    const auto &other__(other_);                                               \
    const bboxset &other(other__);                                             \
                                                                               \
    bboxset1 decoded_subset0;                                                  \
    bboxset1 decoded_subset1;                                                  \
    assert(decoded_subset0.empty());                                           \
    assert(decoded_subset1.empty());                                           \
                                                                               \
    typedef typename subsets_t::const_iterator subsets_iter_t;                 \
    subsets_iter_t iter0 = subsets.begin();                                    \
    subsets_iter_t iter1 = other.subsets.begin();                              \
    subsets_iter_t const end0 = subsets.end();                                 \
    subsets_iter_t const end1 = other.subsets.end();                           \
    while (iter0 != end0 or iter1 != end1) {                                   \
      const T next_pos0 =                                                      \
          iter0 != end0 ? iter0->first : numeric_limits<T>::max();             \
      const T next_pos1 =                                                      \
          iter1 != end1 ? iter1->first : numeric_limits<T>::max();             \
      const T pos = min(next_pos0, next_pos1);                                 \
      const bool active0 = next_pos0 == pos;                                   \
      const bool active1 = next_pos1 == pos;                                   \
      const bboxset1 *const subset0p = active0 ? iter0->second.get() : 0;      \
      const bboxset1 *const subset1p = active1 ? iter1->second.get() : 0;      \
      if (active0)                                                             \
        decoded_subset0 ^= *subset0p;                                          \
      if (active1)                                                             \
        decoded_subset1 ^= *subset1p;                                          \
                                                                               \
      {                                                                        \
        const auto &_1(pos);                                                   \
        const auto &_2(decoded_subset0);                                       \
        const auto &_3(decoded_subset1);                                       \
        f;                                                                     \
      }                                                                        \
                                                                               \
      if (active0)                                                             \
        ++iter0;                                                               \
      if (active1)                                                             \
        ++iter1;                                                               \
    }                                                                          \
    assert(decoded_subset0.empty());                                           \
    assert(decoded_subset1.empty());                                           \
  } while (0)

  template <typename F>
  bboxset binary_operator(const F &op, const bboxset &other) const {
    assert(not is_poison() and not other.is_poison());

#if 0
    // TODO: This assumes that the empty set is a neutral element for
    // the operator
    if (other.empty()) return *this;
    if (empty()) return other;
    assert(all(other.stride == stride));
    assert(all(offset == other.offset));
#endif

    if (not empty() and not other.empty()) {
      assert(all(stride == other.stride));
      assert(all(offset == other.offset));
    }

    bboxset res;
    if (empty()) {
      res.stride = other.stride;
      res.offset = other.offset;
    } else {
      res.stride = stride;
      res.offset = offset;
    }

    bboxset1 old_decoded_subsetr;
#ifndef CARPET_AVOID_LAMBDA
    traverse_subsets(
        [&](const T &pos, const bboxset1 &decoded_subset0,
            const bboxset1 &decoded_subset1) {
          bboxset1 decoded_subsetr = op(decoded_subset0, decoded_subset1);
          const auto subsetrp =
              make_shared<bboxset1>(decoded_subsetr ^ old_decoded_subsetr);
          if (not subsetrp->empty()) {
            res.subsets.insert(res.subsets.end(), make_pair(pos, subsetrp));
          }
          swap(old_decoded_subsetr, decoded_subsetr);
        },
        other);
#else
    assert(0);
#endif
    assert(old_decoded_subsetr.empty());

    return res;
  }
#define BINARY_OPERATOR(op, other_)                                            \
  do {                                                                         \
    const auto &other__(other_);                                               \
    const bboxset &other(other__);                                             \
                                                                               \
    assert(not is_poison() and not other.is_poison());                         \
                                                                               \
    if (not empty() and not other.empty()) {                                   \
      assert(all(stride == other.stride));                                     \
      assert(all(offset == other.offset));                                     \
    }                                                                          \
                                                                               \
    bboxset res;                                                               \
    if (empty()) {                                                             \
      res.stride = other.stride;                                               \
      res.offset = other.offset;                                               \
    } else {                                                                   \
      res.stride = stride;                                                     \
      res.offset = offset;                                                     \
    }                                                                          \
                                                                               \
    bboxset1 old_decoded_subsetr;                                              \
    TRAVERSE_SUBSETS2(const T &pos(_1); const bboxset1 &decoded_subset0(_2);   \
                      const bboxset1 &decoded_subset1(_3); {                   \
                        bboxset1 decoded_subsetr;                              \
                        {                                                      \
                          auto &_0(decoded_subsetr);                           \
                          const auto &_1(decoded_subset0);                     \
                          const auto &_2(decoded_subset1);                     \
                          op;                                                  \
                        }                                                      \
                        const auto subsetrp = make_shared<bboxset1>(           \
                            decoded_subsetr ^ old_decoded_subsetr);            \
                        if (not subsetrp->empty()) {                           \
                          res.subsets.insert(res.subsets.end(),                \
                                             make_pair(pos, subsetrp));        \
                        }                                                      \
                        swap(old_decoded_subsetr, decoded_subsetr);            \
                      }, other);                                               \
    assert(old_decoded_subsetr.empty());                                       \
                                                                               \
    _0 = res;                                                                  \
  } while (0)

public:
  /** Invariant */
  bool invariant() const {
    if (is_poison() and not empty())
      return false;
    if (any(stride <= vect(0)))
      return false;
    if (any(offset < vect(0) or offset >= stride))
      return false;
    forall(const auto &pos_subset, subsets) {
      if (pos_subset.second.get()->empty())
        return false;
      if (any(pos_subset.second.get()->stride != init(stride)))
        return false;
      if (any(pos_subset.second.get()->offset != init(offset)))
        return false;
    }
    if (chi_size() % 2 != 0)
      return false;
    return true;
  }

  /** Create empty set */
  bboxset() : stride(vect(1)), offset(vect(0)), is_poison_(false) {}

  /** Return a is_poisoned set */
  static bboxset poison() {
    bboxset bs;
    bs.is_poison_ = true;
    return bs;
  }

  /** Copy constructor */
  bboxset(const bboxset &other);

  /** Copy constructor */
  // bboxset(bboxset&& other);

  /** Assignment */
  bboxset &operator=(const bboxset &other);

  /** Assignment */
  // bboxset& operator=(bboxset&& other);

  /** Create set from bbox */
  bboxset(const bbox &b);

  /** Create set from container of bboxes, bboxsets, or container
      thereof */
  template <typename C> bboxset(const C &elts);

  /** Create set from container of structs containing a bbox, bboxset,
      or container thereof */
  template <typename C, typename S, typename B>
  bboxset(const C &elts, const B S::*const mptr);

  /** Test for is_poison */
  bool is_poison() const { return is_poison_; }

  /** Test for emptiness */
  bool empty() const {
    assert(not is_poison());
    return subsets.empty();
  }

  /** Find a point contained in this set */
  vect front() const;

  /** Number of characteristic points */
  T chi_size() const;

  /** Number of elements */
  size_type size() const;

  /** Container (min and max) */
  bbox container() const;

  /** Test for equality */
  bool operator==(const bboxset &other) const;

  /** Test for is-subset-of */
  bool operator<=(const bboxset &other) const;

  /** Test for is-strict-subset-of */
  bool operator<(const bboxset &other) const;

  bool operator!=(const bboxset &other) const { return not(*this == other); }

  bool operator>=(const bboxset &other) const { return other <= *this; }

  bool operator>(const bboxset &other) const { return other < *this; }

  bool contains(const vect &v) const { return bbox(v, v, stride) <= *this; }

  bool intersects(const bbox &b) const { return not(b & *this).empty(); }

  /** Symmetric set difference */
  bboxset operator^(const bboxset &other) const;

  /** Symmetric set difference */
  bboxset &operator^=(const bboxset &other);

  /** Set intersection */
  bboxset operator&(const bboxset &other) const;

  /** Set intersection */
  bboxset &operator&=(const bboxset &other);

  /** Set Union */
  bboxset operator|(const bboxset &other) const;

  /** Set union */
  bboxset &operator|=(const bboxset &other);

  /** Symmetric set union */
  bboxset operator+(const bboxset &other) const;

  /** Symmetric set union */
  bboxset &operator+=(const bboxset &other);

  /** Set difference */
  bboxset operator-(const bboxset &other) const;

  /** Set difference */
  bboxset &operator-=(const bboxset &other);

  /** Operators taking bboxes as arguments */
  bboxset &operator&=(const bbox &b) { return *this &= bboxset(b); }
  bboxset &operator^=(const bbox &b) { return *this ^= bboxset(b); }
  bboxset &operator|=(const bbox &b) { return *this |= bboxset(b); }
  bboxset &operator+=(const bbox &b) { return *this += bboxset(b); }
  bboxset &operator-=(const bbox &b) { return *this -= bboxset(b); }

  /** Shift all points */
  bboxset shift(const vect &dist, const vect &dist_denom = vect(1)) const;

  /** Expand the set (convolute with a bbox) */
  bboxset expand(const vect &lo, const vect &hi) const;

  /** Expand the set (convolute with a bbox) */
  bboxset expand(const ::vect<vect, 2> &lohi) const {
    return expand(lohi[0], lohi[1]);
  }

private:
  static bool strides_are_compatible(const vect &str1, const vect &str2);

public:
  /** Expande the set (changing the stride): find the smallest
      b-compatible bboxset around this bboxset */
  // Possible implementation:
  // - convert bboxset into set of bboxes
  // - expand each bbox
  // - create new bboxset as union of these bboxes
  // Alternative, possibly more efficient implementation:
  // - for each dimension:
  //   - BS1 := shifted left to align with bbox
  //   - BS2 := shifted right to align with bbox
  //   - create union of BS1 and BS2
  bboxset expanded_for(const bbox &target) const;

  /** Shrink the set (changing the stride): find the largest
      b-compatible bboxset inside this bboxset */
  bboxset contracted_for(const bbox &target) const;

  /** Expande the set (changing the stride): find the smallest open
      b-compatible bboxset around this bboxset */
  bboxset anti_contracted_for(const bbox &target) const;

  /** Serialise the set */
  template <typename C> void serialise(C &out) const;

  /** Iterate over a serialised set */
private:
  typedef vector<bbox> iter_memo_t;

public:
  typedef typename iter_memo_t::const_iterator const_iterator;

private:
  mutable iter_memo_t iter_memo;

public:
  int setsize() const {
    iter_memo_t im;
    serialise(im);
    return im.size();
  }
  const_iterator begin() const {
    iter_memo.clear();
    serialise(iter_memo);
    return iter_memo.begin();
  }
  const_iterator end() const { return iter_memo.end(); }

  /** Memory usage */
  size_t memory() const;

  /** Input */
  istream &input(istream &is);

  /** Output */
  ostream &debug_output(ostream &os) const;
  ostream &output(ostream &os) const;
};

template <typename T> class bboxset<T, 0> {
  template <typename, int> friend class bboxset;
  typedef ::vect<T, 0> vect;
  typedef ::bbox<T, 0> bbox;

  bool state;
  vect stride, offset;
  bool is_poison_;

public:
  bool invariant() const {
    if (is_poison() and state)
      return false;
    return true;
  }

  bboxset() : state(false), is_poison_(false) {}
  static bboxset poison() {
    bboxset bs;
    bs.is_poison_ = true;
    return bs;
  }
  bboxset(const bbox &b) : state(true), is_poison_(false) {
    assert(not b.is_poison());
  }
  bboxset(const bboxset &other)
      : state(other.state), is_poison_(other.is_poison_) {}
  bboxset &operator=(const bboxset &other) {
    state = other.state;
    is_poison_ = other.is_poison_;
    return *this;
  }

  bool is_poison() const { return is_poison_; }
  bool empty() const {
    assert(not is_poison());
    return not state;
  }
  vect front() const {
    assert(not is_poison());
    assert(not empty());
    return vect();
  }
  T chi_size() const {
    assert(not is_poison());
    return 1;
  }
  T size() const {
    assert(not is_poison());
    return state;
  }
  bbox container() const {
    assert(not is_poison());
    return bbox();
  }

  bool operator==(const bboxset &other) const {
    assert(not is_poison() and not other.is_poison());
    return state == other.state;
  }
  bool operator<=(const bboxset &other) const {
    assert(not is_poison() and not other.is_poison());
    return state <= other.state;
  }
  bool operator!=(const bboxset &other) const { return not(*this == other); }
  bool operator>=(const bboxset &other) const { return other <= *this; }
  bool operator<(const bboxset &other) const { return not(*this >= other); }
  bool operator>(const bboxset &other) const { return not(*this <= other); }

  bboxset &operator&=(const bboxset &other) {
    assert(not is_poison() and not other.is_poison());
    state &= other.state;
    return *this;
  }
  bboxset operator&(const bboxset &other) const {
    bboxset res = *this;
    return res &= other;
  }
  bboxset &operator^=(const bboxset &other) {
    assert(not is_poison() and not other.is_poison());
    state ^= other.state;
    return *this;
  }
  bboxset operator^(const bboxset &other) const {
    bboxset res = *this;
    return res ^= other;
  }
  bboxset &operator|=(const bboxset &other) {
    assert(not is_poison() and not other.is_poison());
    state |= other.state;
    return *this;
  }
  bboxset operator|(const bboxset &other) const {
    bboxset res = *this;
    return res |= other;
  }
  bboxset &operator+=(const bboxset &other) {
    assert(not is_poison() and not other.is_poison());
    assert(not(state and other.state));
    return *this |= other;
  }
  bboxset operator+(const bboxset &other) const {
    bboxset res = *this;
    return res += other;
  }
  bboxset &operator-=(const bboxset &other) {
    assert(not is_poison() and not other.is_poison());
    state &= not other.state;
    return *this;
  }
  bboxset operator-(const bboxset &other) const {
    bboxset res = *this;
    return res -= other;
  }

  bboxset &operator&=(const bbox &b) { return *this &= bboxset(b); }
  bboxset &operator^=(const bbox &b) { return *this ^= bboxset(b); }
  bboxset &operator|=(const bbox &b) { return *this |= bboxset(b); }
  bboxset &operator+=(const bbox &b) { return *this += bboxset(b); }
  bboxset &operator-=(const bbox &b) { return *this -= bboxset(b); }

  bboxset shift(const vect &dist, const vect &dist_denom = vect(1)) const {
    assert(not is_poison());
    return *this;
  }
  bboxset expand(const vect &lo, const vect &hi) const {
    assert(not is_poison());
    return *this;
  }
  bboxset expand(const ::vect<vect, 2> &lohi) const {
    return expand(lohi[0], lohi[1]);
  }

  template <typename C> void serialise(C &out) const {
    assert(not is_poison());
    if (state)
      out.insert(out.end(), bbox());
  }

  size_t memory() const { return sizeof *this; }

  istream &input(istream &is);

  ostream &debug_output(ostream &os) const;
  ostream &output(ostream &os) const;
};

/** Operators involving bboxes */
template <typename T, int D>
inline bboxset<T, D> operator&(const bboxset<T, D> &bs, const bbox<T, D> &b) {
  return bs & bboxset<T, D>(b);
}
template <typename T, int D>
inline bboxset<T, D> operator&(const bbox<T, D> &b, const bboxset<T, D> &bs) {
  return bboxset<T, D>(b) & bs;
}
// Note: bbox & bbox -> bbox

template <typename T, int D>
inline bboxset<T, D> operator^(const bboxset<T, D> &bs, const bbox<T, D> &b) {
  return bs ^ bboxset<T, D>(b);
} template <typename T, int D>
inline bboxset<T, D> operator^(const bbox<T, D> &b, const bboxset<T, D> &bs) {
  return bboxset<T, D>(b) ^ bs;
} template <typename T, int D>
inline bboxset<T, D> operator^(const bbox<T, D> &b1, const bbox<T, D> &b2) {
  return bboxset<T, D>(b1) ^ bboxset<T, D>(b2);
}

template <typename T, int D>
inline bboxset<T, D> operator|(const bboxset<T, D> &bs, const bbox<T, D> &b) {
  return bs | bboxset<T, D>(b);
}
template <typename T, int D>
inline bboxset<T, D> operator|(const bbox<T, D> &b, const bboxset<T, D> &bs) {
  return bboxset<T, D>(b) | bs;
}
template <typename T, int D>
inline bboxset<T, D> operator|(const bbox<T, D> &b1, const bbox<T, D> &b2) {
  return bboxset<T, D>(b1) | bboxset<T, D>(b2);
}

template <typename T, int D>
inline bboxset<T, D> operator+(const bboxset<T, D> &bs, const bbox<T, D> &b) {
  return bs + bboxset<T, D>(b);
}
template <typename T, int D>
inline bboxset<T, D> operator+(const bbox<T, D> &b, const bboxset<T, D> &bs) {
  return bboxset<T, D>(b) + bs;
}
template <typename T, int D>
inline bboxset<T, D> operator+(const bbox<T, D> &b1, const bbox<T, D> &b2) {
  return bboxset<T, D>(b1) + bboxset<T, D>(b2);
}

template <typename T, int D>
inline bboxset<T, D> operator-(const bboxset<T, D> &bs, const bbox<T, D> &b) {
  return bs - bboxset<T, D>(b);
}
template <typename T, int D>
inline bboxset<T, D> operator-(const bbox<T, D> &b, const bboxset<T, D> &bs) {
  return bboxset<T, D>(b) - bs;
}
template <typename T, int D>
inline bboxset<T, D> operator-(const bbox<T, D> &b1, const bbox<T, D> &b2) {
  return bboxset<T, D>(b1) - bboxset<T, D>(b2);
}

template <typename T, int D>
inline bool operator==(const bbox<T, D> &b, const bboxset<T, D> &bs) {
  return bboxset<T, D>(b) == bs;
}
template <typename T, int D>
inline bool operator==(const bboxset<T, D> &bs, const bbox<T, D> &b) {
  return bs == bboxset<T, D>(b);
}

template <typename T, int D>
inline bool operator!=(const bbox<T, D> &b, const bboxset<T, D> &bs) {
  return bboxset<T, D>(b) != bs;
}
template <typename T, int D>
inline bool operator!=(const bboxset<T, D> &bs, const bbox<T, D> &b) {
  return bs != bboxset<T, D>(b);
}

template <typename T, int D>
inline bool operator<=(const bbox<T, D> &b, const bboxset<T, D> &bs) {
  return bboxset<T, D>(b) <= bs;
}
template <typename T, int D>
inline bool operator<=(const bboxset<T, D> &bs, const bbox<T, D> &b) {
  return bs <= bboxset<T, D>(b);
}

template <typename T, int D>
inline bool operator<(const bbox<T, D> &b, const bboxset<T, D> &bs) {
  return bboxset<T, D>(b) < bs;
}
template <typename T, int D>
inline bool operator<(const bboxset<T, D> &bs, const bbox<T, D> &b) {
  return bs < bboxset<T, D>(b);
}

template <typename T, int D>
inline bool operator>=(const bbox<T, D> &b, const bboxset<T, D> &bs) {
  return bboxset<T, D>(b) >= bs;
}
template <typename T, int D>
inline bool operator>=(const bboxset<T, D> &bs, const bbox<T, D> &b) {
  return bs >= bboxset<T, D>(b);
}

template <typename T, int D>
inline bool operator>(const bbox<T, D> &b, const bboxset<T, D> &bs) {
  return bboxset<T, D>(b) > bs;
}
template <typename T, int D>
inline bool operator>(const bboxset<T, D> &bs, const bbox<T, D> &b) {
  return bs > bboxset<T, D>(b);
}

/** Memory usage */
template <typename T, int D> inline size_t memoryof(bboxset<T, D> const &bs) {
  return bs.memory();
}

/** Input */
template <typename T, int D>
inline istream &operator>>(istream &is, bboxset<T, D> &bs) {
  return bs.input(is);
}

/** Output */
template <typename T, int D>
inline ostream &operator<<(ostream &os, const bboxset<T, D> &bs) {
  return bs.output(os);
}

//////////////////////////////////////////////////////////////////////////////
// bboxset<T,D>
//////////////////////////////////////////////////////////////////////////////

/** Copy constructor */
template <typename T, int D>
bboxset<T, D>::bboxset(const bboxset &other)
    : stride(other.stride), offset(other.offset), is_poison_(other.is_poison_) {
  forall(const auto &pos_subset, other.subsets) {
    const T &pos = pos_subset.first;
    const auto new_subsetp = make_shared<bboxset1>(*pos_subset.second.get());
    subsets.insert(subsets.end(), make_pair(pos, new_subsetp));
  }
}

/** Copy constructor */
// bboxset(bboxset&& other): stride(other.stride), offset(other.offset)
// {
//   swap(subsets, other.subsets);
// }

/** Assignment */
template <typename T, int D>
bboxset<T, D> &bboxset<T, D>::operator=(const bboxset &other) {
  if (&other == this)
    return *this;
  subsets.clear();
  forall(const auto &pos_subset, other.subsets) {
    const T &pos = pos_subset.first;
    const auto new_subsetp = make_shared<bboxset1>(*pos_subset.second.get());
    subsets.insert(subsets.end(), make_pair(pos, new_subsetp));
  }
  stride = other.stride;
  offset = other.offset;
  is_poison_ = other.is_poison_;
  return *this;
}

/** Assignment */
// bboxset& operator=(bboxset&& other)
// {
//   if (&other == this) return *this;
//   swap(subsets, other.subsets);
//   stride = other.stride;
//   offset = other.offset;
//   return *this;
// }

/** Create set from bbox */
template <typename T, int D>
bboxset<T, D>::bboxset(const bbox &b)
    : stride(b.stride()), offset(imod(b.lower(), b.stride())),
      is_poison_(false) {
  assert(not b.is_poison());
  if (b.empty())
    return;
  const T lo = last(b.lower());
  const T hi = last(b.upper()) + last(b.stride());
  const bbox1 b1(init(b.lower()), init(b.upper()), init(b.stride()));
  const auto lo_subsetp = make_shared<bboxset1>(b1);
  const auto hi_subsetp = make_shared<bboxset1>(b1);
  // subsets.emplace_hint(subsets.end(), lo, lo_subsetp);
  // subsets.emplace_hint(subsets.end(), hi, hi_subsetp);
  subsets.insert(subsets.end(), make_pair(lo, lo_subsetp));
  subsets.insert(subsets.end(), make_pair(hi, hi_subsetp));
}

/** Create set from container of bboxes, bboxsets, or container
    thereof */
template <typename T, int D>
template <typename C>
bboxset<T, D>::bboxset(const C &elts) {
  *this = bboxset();
  forall(const auto &elt, elts) { *this += bboxset(elt); }
}

/** Create set from container of structs containing a bbox, bboxset,
    or container thereof */
template <typename T, int D>
template <typename C, typename S, typename B>
bboxset<T, D>::bboxset(const C &elts, const B S::*const mptr) {
  *this = bboxset();
  forall(const auto &elt, elts) { *this += bboxset(elt.*mptr); }
}

/** Find a point contained in this set */
template <typename T, int D> vect<T, D> bboxset<T, D>::front() const {
  assert(not is_poison());
  assert(not empty());
  const auto &pos_subset = *subsets.begin();
  const auto &pos = pos_subset.first;
  const auto &subset = *pos_subset.second.get();
  const vect1 f1 = subset.front();
  return vect(f1, pos);
}

/** Number of characteristic points */
template <typename T, int D> T bboxset<T, D>::chi_size() const {
  assert(not is_poison());
  T sum = 0;
  forall(const auto &pos_subset, subsets) {
    sum += pos_subset.second->chi_size();
  }
  return sum;
}

/** Number of elements */
template <typename T, int D> size_type bboxset<T, D>::size() const {
  assert(not is_poison());
  size_type total_size = 0;             // accumulated total number of points
  T old_pos = numeric_limits<T>::min(); // location of last subset
  size_t old_subset_size = 0;           // number of points in the last subset
#ifndef CARPET_AVOID_LAMBDA
  traverse_subsets([&](const T &pos, const bboxset1 &subset) {
    const size_type subset_size = subset.size();
    if (old_subset_size > 0) {
      assert((pos - old_pos) % stride[D - 1] == 0);
      total_size += (pos - old_pos) / stride[D - 1] * old_subset_size;
    }
    old_pos = pos;
    old_subset_size = subset_size;
  });
#else
  TRAVERSE_SUBSETS1(const T &pos(_1); const bboxset1 &subset(_2); {
    const size_type subset_size = subset.size();
    if (old_subset_size > 0) {
      assert((pos - old_pos) % stride[D - 1] == 0);
      total_size += (pos - old_pos) / stride[D - 1] * old_subset_size;
    }
    old_pos = pos;
    old_subset_size = subset_size;
  });
#endif
  assert(old_subset_size == 0);
  return total_size;
}

/** Container (min and max) */
template <typename T, int D> bbox<T, D> bboxset<T, D>::container() const {
  assert(not is_poison());
  if (empty())
    return bbox();
  const T lo = subsets.begin()->first;
  const T hi = subsets.rbegin()->first;
  bbox1 container1;
#ifndef CARPET_AVOID_LAMBDA
  traverse_subsets([&](const T &pos, const bboxset1 &subset) {
    container1 = container1.expanded_containing(subset.container());
  });
#else
  TRAVERSE_SUBSETS1(const T &pos(_1); const bboxset1 &subset(_2); {
    container1 = container1.expanded_containing(subset.container());
  });
#endif
  return bbox(vect(container1.lower(), lo),
              vect(container1.upper(), hi - last(stride)),
              vect(container1.stride(), last(stride)));
}

/** Test for equality */
template <typename T, int D>
bool bboxset<T, D>::operator==(const bboxset &other) const {
  return (*this ^ other).empty();
}

/** Test for is-subset-of */
template <typename T, int D>
bool bboxset<T, D>::operator<=(const bboxset &other) const {
  return (*this | other) == other;
}

/** Test for is-strict-subset-of */
template <typename T, int D>
bool bboxset<T, D>::operator<(const bboxset &other) const {
  return *this != other and *this <= other;
}

/** Symmetric set difference */
template <typename T, int D>
bboxset<T, D> bboxset<T, D>::operator^(const bboxset &other) const {
// TODO: If other is much smaller than this, direct insertion may
// be faster
#ifndef CARPET_AVOID_LAMBDA
  return binary_operator([](const bboxset1 &set0, const bboxset1 &set1) {
    return set0 ^ set1;
  }, other);
#else
  bboxset _0;
  BINARY_OPERATOR(const bboxset1 &set0(_1); const bboxset1 &set1(_2);
                  { _0 = set0 ^ set1; }, other);
  return _0;
#endif
}

/** Symmetric set difference */
template <typename T, int D>
bboxset<T, D> &bboxset<T, D>::operator^=(const bboxset &other) {
  bboxset res = *this;
  return *this = *this ^ other;
}

/** Set intersection */
template <typename T, int D>
bboxset<T, D> bboxset<T, D>::operator&(const bboxset &other) const {
#ifndef CARPET_AVOID_LAMBDA
  return binary_operator([](const bboxset1 &set0, const bboxset1 &set1) {
    return set0 & set1;
  }, other);
#else
  bboxset _0;
  BINARY_OPERATOR(const bboxset1 &set0(_1); const bboxset1 &set1(_2);
                  { _0 = set0 & set1; }, other);
  return _0;
#endif
}

/** Set intersection */
template <typename T, int D>
bboxset<T, D> &bboxset<T, D>::operator&=(const bboxset &other) {
  return *this = *this & other;
}

/** Set Union */
template <typename T, int D>
bboxset<T, D> bboxset<T, D>::operator|(const bboxset &other) const {
#ifndef CARPET_AVOID_LAMBDA
  return binary_operator([](const bboxset1 &set0, const bboxset1 &set1) {
    return set0 | set1;
  }, other);
#else
  bboxset _0;
  BINARY_OPERATOR(const bboxset1 &set0(_1); const bboxset1 &set1(_2);
                  { _0 = set0 | set1; }, other);
  return _0;
#endif
}

/** Set union */
template <typename T, int D>
bboxset<T, D> &bboxset<T, D>::operator|=(const bboxset &other) {
  return *this = *this | other;
}

/** Symmetric set union */
template <typename T, int D>
bboxset<T, D> bboxset<T, D>::operator+(const bboxset &other) const {
// return binary_operator
//   ([](const bboxset1& set0, const bboxset1& set1) { return set0 + set1; },
//    other);
#ifdef CARPET_DEBUG
  assert((*this & other).empty());
#endif
  // Since the sets are disjoint, their symmetric set union is
  // identical to their symmetric difference.
  return *this ^ other;
}

/** Symmetric set union */
template <typename T, int D>
bboxset<T, D> &bboxset<T, D>::operator+=(const bboxset &other) {
  return *this = *this + other;
}

/** Set difference */
template <typename T, int D>
bboxset<T, D> bboxset<T, D>::operator-(const bboxset &other) const {
#ifndef CARPET_AVOID_LAMBDA
  return binary_operator([](const bboxset1 &set0, const bboxset1 &set1) {
    return set0 - set1;
  }, other);
#else
  bboxset _0;
  BINARY_OPERATOR(const bboxset1 &set0(_1); const bboxset1 &set1(_2);
                  { _0 = set0 - set1; }, other);
  return _0;
#endif
}

/** Set difference */
template <typename T, int D>
bboxset<T, D> &bboxset<T, D>::operator-=(const bboxset &other) {
  return *this = *this - other;
}

/** Shift all points */
template <typename T, int D>
bboxset<T, D> bboxset<T, D>::shift(const vect &dist,
                                   const vect &dist_denom) const {
  assert(not is_poison());
  assert(all(stride % dist_denom == 0));
  if (all(dist == 0))
    return *this;
  bboxset res;
  res.stride = stride;
  res.offset = imod(offset + dist * stride / dist_denom, res.stride);
  forall(const auto &pos_subset, subsets) {
    const T &pos = pos_subset.first;
    const bboxset1 &subset = *pos_subset.second.get();
    const T new_pos = pos + last(dist * stride / dist_denom);
    const auto new_subsetp =
        make_shared<bboxset1>(subset.shift(init(dist), init(dist_denom)));
    res.subsets.insert(res.subsets.end(), make_pair(new_pos, new_subsetp));
  }
  return res;
}

/** Expand the set (convolute with a bbox) */
template <typename T, int D>
bboxset<T, D> bboxset<T, D>::expand(const vect &lo, const vect &hi) const {
  assert(not is_poison());
  assert(all(lo >= 0 and hi >= 0));
  bboxset res = shift(-lo);
  for (int d = 0; d < D; ++d) {
    T to_expand = (hi + lo)[d];
    T current_size = 1;
    while (to_expand > 0) {
      const T this_expand = min(to_expand, current_size);
      res |= res.shift(vect::dir(d) * this_expand);
      current_size += this_expand;
      to_expand -= this_expand;
    }
    assert(to_expand == 0);
  }
  return res;
}

template <typename T, int D>
bool bboxset<T, D>::strides_are_compatible(const vect &str1, const vect &str2) {
  if (all(str1 >= str2)) {
    return all(imod(str1, str2) == 0);
  } else if (all(str1 <= str2)) {
    return all(imod(str2, str1) == 0);
  }
  return false;
}

/** Expande the set (changing the stride): find the smallest
    b-compatible bboxset around this bboxset */
template <typename T, int D>
bboxset<T, D> bboxset<T, D>::expanded_for(const bbox &target) const {
  // Check preconditions
  assert(not is_poison() and not target.is_poison());
  assert(not target.empty());
  assert(strides_are_compatible(stride, target.stride()));
  vector<bbox> bs;
  serialise(bs);
  bboxset res;
  forall(const auto &b, bs) { res |= b.expanded_for(target); }
  return res;
}

/** Shrink the set (changing the stride): find the largest
    b-compatible bboxset inside this bboxset */
template <typename T, int D>
bboxset<T, D> bboxset<T, D>::contracted_for(const bbox &target) const {
  // Check preconditions
  assert(not is_poison() and not target.is_poison());
  assert(not target.empty());
  assert(strides_are_compatible(stride, target.stride()));
  const bbox cont = container();
  const vect safety = 10;
  const bbox good_world = cont.anti_contracted_for(target).expand(safety);
  const bbox world1 = good_world.anti_contracted_for(cont).expand(safety);
  const bbox world2 = world1.anti_contracted_for(target).expand(safety);
  return (world2 ^ (world1 ^ *this).anti_contracted_for(target)) & good_world;
}

/** Expande the set (changing the stride): find the smallest open
    b-compatible bboxset around this bboxset */
template <typename T, int D>
bboxset<T, D> bboxset<T, D>::anti_contracted_for(const bbox &target) const {
  // Check preconditions
  assert(not is_poison() and not target.is_poison());
  assert(not target.empty());
  assert(strides_are_compatible(stride, target.stride()));
  vector<bbox> bs;
  serialise(bs);
  bboxset res;
  forall(const auto &b, bs) { res |= b.anti_contracted_for(target); }
  return res;
}

/** Serialise the set */
template <typename T, int D>
template <typename C>
void bboxset<T, D>::serialise(C &out) const {
  assert(not is_poison());
  typedef map<bbox1, T> subboxes_t;
  typedef set<bbox1> subboxes1_t;
  // TODO: Instead of copying from old_subboxes to subboxes,
  // maintain subboxes via a non-const iterator
  subboxes_t old_subboxes;
#ifndef CARPET_AVOID_LAMBDA
  traverse_subsets([&](const T &pos, const bboxset1 &subset) {
    // Convert subset to bboxes
    subboxes1_t subboxes1;
    subset.serialise(subboxes1);

    const bbox1 dummy1;
    subboxes_t subboxes;
    // subboxes.reserve(old_subboxes.size() + subboxes1.size());
    typedef typename subboxes_t::const_iterator subboxes_iter_t;
    typedef typename subboxes1_t::const_iterator subboxes1_iter_t;
    subboxes_iter_t iter0 = old_subboxes.begin();
    subboxes1_iter_t iter1 = subboxes1.begin();
    subboxes_iter_t const end0 = old_subboxes.end();
    subboxes1_iter_t const end1 = subboxes1.end();
    while (iter0 != end0 or iter1 != end1) {
      bool active0 = iter0 != end0;
      bool active1 = iter1 != end1;
      const bbox1 &subbox0 = active0 ? iter0->first : dummy1;
      const bbox1 &subbox1 = active1 ? *iter1 : dummy1;
      // When both subboxes are active, keep only the first (as
      // determined by less<>)
      if (active0 and active1) {
        if (subbox0 != subbox1) {
          /*const*/ less<bbox1> bbox1_less;
          if (bbox1_less(subbox0, subbox1)) {
            active1 = false;
          } else {
            active0 = false;
          }
        }
      }
      const T old_pos = active0 ? iter0->second : T();

      if ((active0 and active1) and (subbox0 == subbox1)) {
        // The current bbox continues unchanged -- keep it
        subboxes.insert(subboxes.end(), *iter0);
      } else {
        if (active0) {
          // The current box changed; finalize it
          const bbox new_box(vect(subbox0.lower(), old_pos),
                             vect(subbox0.upper(), pos - last(stride)), stride);
          out.insert(out.end(), new_box);
        }
        if (active1) {
          // There is a new box; add it
          subboxes.insert(subboxes.end(), make_pair(subbox1, pos));
        }
      }

      if (active0)
        ++iter0;
      if (active1)
        ++iter1;
    }
    swap(old_subboxes, subboxes);
  });
#else
  TRAVERSE_SUBSETS1(const T &pos(_1); const bboxset1 &subset(_2); {
    // Convert subset to bboxes
    subboxes1_t subboxes1;
    subset.serialise(subboxes1);

    const bbox1 dummy1;
    subboxes_t subboxes;
    // subboxes.reserve(old_subboxes.size() + subboxes1.size());
    typedef typename subboxes_t::const_iterator subboxes_iter_t;
    typedef typename subboxes1_t::const_iterator subboxes1_iter_t;
    subboxes_iter_t iter0 = old_subboxes.begin();
    subboxes1_iter_t iter1 = subboxes1.begin();
    subboxes_iter_t const end0 = old_subboxes.end();
    subboxes1_iter_t const end1 = subboxes1.end();
    while (iter0 != end0 or iter1 != end1) {
      bool active0 = iter0 != end0;
      bool active1 = iter1 != end1;
      const bbox1 &subbox0 = active0 ? iter0->first : dummy1;
      const bbox1 &subbox1 = active1 ? *iter1 : dummy1;
      // When both subboxes are active, keep only the first (as
      // determined by less<>)
      if (active0 and active1) {
        if (subbox0 != subbox1) {
          /*const*/ less<bbox1> bbox1_less;
          if (bbox1_less(subbox0, subbox1)) {
            active1 = false;
          } else {
            active0 = false;
          }
        }
      }
      const T old_pos = active0 ? iter0->second : T();

      if ((active0 and active1) and (subbox0 == subbox1)) {
        // The current bbox continues unchanged -- keep it
        subboxes.insert(subboxes.end(), *iter0);
      } else {
        if (active0) {
          // The current box changed; finalize it
          const bbox new_box(vect(subbox0.lower(), old_pos),
                             vect(subbox0.upper(), pos - last(stride)), stride);
          out.insert(out.end(), new_box);
        }
        if (active1) {
          // There is a new box; add it
          subboxes.insert(subboxes.end(), make_pair(subbox1, pos));
        }
      }

      if (active0)
        ++iter0;
      if (active1)
        ++iter1;
    }
    swap(old_subboxes, subboxes);
  });
#endif
  assert(old_subboxes.empty());
}

/** Memory usage */
template <typename T, int D> size_t bboxset<T, D>::memory() const {
  size_t s = sizeof *this;
  forall(const auto &pos_subset, subsets) {
    s += sizeof pos_subset;
    const auto *const subsetp = pos_subset.second.get();
    s += memoryof(*subsetp);
  }
  return s;
}

/** Input */
template <typename T, int D> istream &bboxset<T, D>::input(istream &is) {
  T Tdummy;
  try {
    skipws(is);
    consume(is, "bboxset<");
    consume(is, typestring(Tdummy));
    consume(is, ",");
    int D_;
    is >> D_;
    if (D_ != D) {
      ostringstream msg;
      msg << "Input error: Wrong bboxset dimension " << D_ << ", expected "
          << D;
      CCTK_WARN(CCTK_WARN_ALERT, msg.str().c_str());
      throw input_error();
    }
    consume(is, ">(");
    consume(is, "set<bbox>:");
    set<bbox> bs;
    is >> bs;
    *this = bboxset(bs);
    consume(is, ",");
    consume(is, "stride:");
    is >> stride;
    consume(is, ",");
    consume(is, "offset:");
    is >> offset;
    consume(is, ")");
  } catch (input_error &err) {
    ostringstream msg;
    msg << "Input error while reading a bboxset<" << typestring(Tdummy)
        << ",0>";
    CCTK_WARN(CCTK_WARN_ALERT, msg.str().c_str());
    throw err;
  }
  is_poison_ = false;
  return is;
}

/** Output */
template <typename T, int D>
ostream &bboxset<T, D>::debug_output(ostream &os) const {
  T Tdummy;
  os << "bboxset[debug]<" << typestring(Tdummy) << "," << D << ">("
     << "subsets:{";
  bool first = true;
  forall(const auto &pos_subset, subsets) {
    if (not first)
      os << ",";
    first = false;
    os << pos_subset.first << ":";
    pos_subset.second.get()->debug_output(os);
  }
  os << "},"
     << "stride:" << stride << ","
     << "offset:" << offset << ","
     << "is_poison_:" << is_poison_ << ")";
  return os;
}

template <typename T, int D> ostream &bboxset<T, D>::output(ostream &os) const {
  assert(not is_poison());
  T Tdummy;
  set<bbox> bs;
  serialise(bs);
  return os << "bboxset<" << typestring(Tdummy) << "," << D << ">("
            << "set<bbox>:" << bs << ","
            << "stride:" << stride << ","
            << "offset:" << offset << ")";
}

//////////////////////////////////////////////////////////////////////////////
// bboxset<T,0>
//////////////////////////////////////////////////////////////////////////////

template <typename T> istream &bboxset<T, 0>::input(istream &is) {
  T Tdummy;
  try {
    skipws(is);
    consume(is, "bboxset<");
    consume(is, typestring(Tdummy));
    consume(is, ",0>(");
    consume(is, "state:");
    is >> state;
    consume(is, ",");
    consume(is, "stride:");
    is >> stride;
    consume(is, ",");
    consume(is, "offset:");
    is >> offset;
    consume(is, ")");
  } catch (input_error &err) {
    ostringstream msg;
    msg << "Input error while reading a bboxset<" << typestring(Tdummy)
        << ",0>";
    CCTK_WARN(CCTK_WARN_ALERT, msg.str().c_str());
    throw err;
  }
  is_poison_ = false;
  return is;
}

template <typename T> ostream &bboxset<T, 0>::debug_output(ostream &os) const {
  T Tdummy;
  return os << "bboxset[debug]<" << typestring(Tdummy) << ",0>("
            << "state:" << state << ","
            << "stride:" << stride << ","
            << "offset:" << offset << ","
            << "is_poison_:" << is_poison_ << ")";
}

template <typename T> ostream &bboxset<T, 0>::output(ostream &os) const {
  assert(not is_poison());
  T Tdummy;
  return os << "bboxset<" << typestring(Tdummy) << ",0>("
            << "state:" << state << ","
            << "stride:" << stride << ","
            << "offset:" << offset << ")";
}

} // namespace bboxset2

#endif // #ifdef CARPET_ENABLE_BBOXSET2

#endif // #ifndef BBOXSET2_HH
