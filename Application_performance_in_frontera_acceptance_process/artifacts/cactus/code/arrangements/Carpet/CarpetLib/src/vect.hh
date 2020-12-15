#ifndef VECT_HH
#define VECT_HH

#include <cctk.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <iostream>

#include "defs.hh"
#include "dist.hh"
#include "vect_helpers.hh"

using namespace std;

#ifdef CARPET_DEBUG
#define ASSERT_VECT(x) assert(x)
#else
#define ASSERT_VECT(x)
#endif

// Forward definition
template <typename T, int D> class vect;

// Input/Output
template <typename T, int D> istream &operator>>(istream &is, vect<T, D> &a);
template <typename T, int D>
ostream &operator<<(ostream &os, const vect<T, D> &a);

/**
 * A short vector with a size that is specified at compile time.
 */
template <typename T, int D> class vect {

  // Fields

  /** Vector elements.  */
  T elt[D == 0 ? 1 : D];

public:
  // Constructors

  /** Explicit empty constructor.  */
  explicit vect() {}

  /** Create a vector from a lower-dimensional vector.  */
  vect(const vect<T, D == 0 ? 0 : D - 1> &x, const T &a) {
    assert(D > 0);
    for (int d = 0; d < D - 1; ++d)
      elt[d] = x[d];
    elt[D - 1] = a;
  }

  /** Create a vector from a lower-dimensional vector.  */
  vect(const T &a, const vect<T, D == 0 ? 0 : D - 1> &x) {
    assert(D > 0);
    elt[0] = a;
    for (int d = 0; d < D - 1; ++d)
      elt[d + 1] = x[d];
  }

  /** Copy constructor.  */
  vect(const vect &a) {
    for (int d = 0; d < D; ++d)
      elt[d] = a.elt[d];
  }

  /** Constructor from a single element.  This constructor might be
      confusing, but it is very convenient.  */
  vect(const T &x) {
    for (int d = 0; d < D; ++d)
      elt[d] = x;
  }

  /** Constructor for 2-element vectors from 2 elements.  */
  vect(const T &x, const T &y) {
    ASSERT_VECT(D == 2);
    // Note: this statement may give "index out of range" warnings.
    // You can safely ignore these.
    elt[0] = x;
    elt[1] = y;
  }

  /** Constructor for 3-element vectors from 3 elements.  */
  vect(const T &x, const T &y, const T &z) {
    ASSERT_VECT(D == 3);
    // Note: this statement may give "index out of range" warnings.
    // You can safely ignore these.
    elt[0] = x;
    elt[1] = y;
    elt[2] = z;
  }

  /** Constructor for 4-element vectors from 4 elements.  */
  vect(const T &x, const T &y, const T &z, const T &t) {
    ASSERT_VECT(D == 4);
    // Note: this statement may give "index out of range" warnings.
    // You can safely ignore these.
    elt[0] = x;
    elt[1] = y;
    elt[2] = z;
    elt[3] = t;
  }

#if 0
  // This creates confusion
  /** Constructor from a pointer, i.e., a C array.  */
  explicit vect (const T* const x)
  {
    for (int d=0; d<D; ++d) elt[d]=x[d];
  }
#endif

#if 0
  // This leads to an ICE on AIX
  template<int E>
  operator vect<vect<T,D>,E> ()
  {
    vect<vect<T,D>,E> r;
    for (int e=0; e<E; ++e) r[e]=*this;
    return r;
  }
#endif

  /** Constructor from a vector with a different type.  */
  template <typename S>
  /*explicit*/ vect(const vect<S, D> &a) {
    for (int d = 0; d < D; ++d)
      elt[d] = (T)a[d];
  }

  /** Create a new 0-element vector with a specific type.  */
  static vect make() {
    ASSERT_VECT(D == 0);
    return vect();
  }

  /** Create a new 1-element vector with a specific type.  */
  static vect make(const T &x) {
    ASSERT_VECT(D == 1);
    return vect(x);
  }

  /** Create a new 2-element vector with a specific type.  */
  static vect make(const T &x, const T &y) {
    ASSERT_VECT(D == 2);
    return vect(x, y);
  }

  /** Create a new 3-element vector with a specific type.  */
  static vect make(const T &x, const T &y, const T &z) {
    ASSERT_VECT(D == 3);
    return vect(x, y, z);
  }

  /** Create a new 4-element vector with a specific type.  */
  static vect make(const T &x, const T &y, const T &z, const T &t) {
    ASSERT_VECT(D == 4);
    return vect(x, y, z, t);
  }

  /** Treat a constant pointer as a reference to a constant vector.  */
  static const vect &ref(const T *const x) { return *(const vect *)x; }

  /** Treat a pointer as a reference to a vector.  */
  static vect &ref(T *const x) { return *(vect *)x; }

  /** Create a vector with one element set to 1 and all other elements
      set to zero.  */
  static vect dir(const int d) {
    vect r = (T)0;
    r[d] = 1;
    return r;
  }

  /** Create a vector with e[i] = i.  */
  static vect seq() {
    vect r;
    for (int d = 0; d < D; ++d)
      r[d] = d;
    return r;
  }

  /** Create a vector with e[i] = n + i.  */
  static vect seq(const int n) {
    vect r;
    for (int d = 0; d < D; ++d)
      r[d] = n + d;
    return r;
  }

  /** Create a vector with e[i] = n + s * i.  */
  static vect seq(const int n, const int s) {
    vect r;
    for (int d = 0; d < D; ++d)
      r[d] = n + s * d;
    return r;
  }

  // Accessors

  /** Return a non-writable element of a vector.  */
  // (Don't return a reference; *this might be a temporary)
  // Do return a reference, so that a vector can be accessed as array
  const T &operator[](const int d) const {
    ASSERT_VECT(d >= 0 && d < D);
    return elt[d];
  }

  /** Return a writable element of a vector as reference.  */
  T &operator[](const int d) {
    ASSERT_VECT(d >= 0 && d < D);
    return elt[d];
  }

#if 0
  // This creates confusion
  /** Return a pointer to a vector.  */
  operator const T* () const
  {
    return this;
  }
#endif

  /** Return a combination of the vector elements e[a[i]].  The
      element combination is selected by another vector.  */
  template <typename TT, int DD>
  vect<T, DD> operator[](const vect<TT, DD> &a) const {
    vect<T, DD> r;
    // (*this)[] performs index checking
    for (int d = 0; d < DD; ++d)
      r[d] = (*this)[a[d]];
    return r;
  }

  // Modifying operators
  DECLARE_MEMBER_OPERATOR_1_REF(operator+=, +=);
  DECLARE_MEMBER_OPERATOR_1_REF(operator-=, -=);
  DECLARE_MEMBER_OPERATOR_1_REF(operator*=, *=);
  DECLARE_MEMBER_OPERATOR_1_REF(operator/=, /=);
  DECLARE_MEMBER_OPERATOR_1_REF(operator%=, %=);
  DECLARE_MEMBER_OPERATOR_1_REF(operator&=, &=);
  DECLARE_MEMBER_OPERATOR_1_REF(operator|=, |=);
  DECLARE_MEMBER_OPERATOR_1_REF(operator^=, ^=);

  // Non-modifying operators

  /** Return a new vector where one element has been replaced.  */
  vect replace(const int d, const T &x) const {
    ASSERT_VECT(d >= 0 && d < D);
    vect r;
    for (int dd = 0; dd < D; ++dd)
      r[dd] = dd == d ? x : elt[dd];
    return r;
  }

  vect reverse() const {
    vect r;
    for (int d = 0; d < D; ++d)
      r[d] = elt[D - 1 - d];
    return r;
  }

  DECLARE_MEMBER_OPERATOR_0(operator+, +)
  DECLARE_MEMBER_OPERATOR_0(operator-, -)
  DECLARE_MEMBER_OPERATOR_0(operator~, ~)
// DECLARE_MEMBER_OPERATOR_0_RET (operator!, !, bool)

#if 0
  /** This corresponds to the ?: operator.  Return a vector with the
      elements set to either a[i] or b[i], depending on whether
      (*this)[i] is true or not.  */
  template<typename TT>
  vect<TT,D> ifthen (const vect<TT,D>& a, const vect<TT,D>& b)
  const
  {
    vect<TT,D> r;
    for (int d=0; d<D; ++d) r[d]=elt[d]?a[d]:b[d];
    return r;
  }
#endif

// Iterators
#if 0
  // This is non-standard
  class iter {
    vect &vec;
    int d;
  public:
    iter (vect &a): vec(a), d(0) { }
    iter& operator++ () { ASSERT_VECT(d<D); ++d; return *this; }
    bool operator bool () const { return d==D; }
    T& operator* () { return vec[d]; }
  };
#endif

  // Memory usage
  size_t memory() const {
    size_t mem = 0;
    for (int d = 0; d < D; ++d)
      mem += memoryof(elt[d]);
    return mem;
  }

  // Input/Output helpers
  void input(istream &is);
  void output(ostream &os) const;

  // MPI
  MPI_Datatype mpi_datatype() const;
};

// Operators

/** This corresponds to the ?: operator.  Return a vector with the
    elements set to either b[i] or c[i], depending on whether a[i] is
    true or not.  */
template <typename S, typename T, int D>
inline vect<T, D> either(const vect<S, D> &a, const vect<T, D> &b,
                         const vect<T, D> &c) {
  vect<T, D> r;
  for (int d = 0; d < D; ++d)
    r[d] = a[d] ? b[d] : c[d];
  return r;
}

template <typename S, typename T, int D>
inline vect<T, D> either(const vect<S, D> &a, const vect<T, D> &b, const T &c) {
  return either(a, b, vect<T, D>(c));
}

template <typename S, typename T, int D>
inline vect<T, D> either(const vect<S, D> &a, const T &b, const vect<T, D> &c) {
  return either(a, vect<T, D>(b), c);
}

template <typename S, typename T, int D>
inline vect<T, D> either(const vect<S, D> &a, const T &b, const T &c) {
  return either(a, vect<T, D>(b), vect<T, D>(c));
}

/** Transpose a vector of a vector */
template <typename T, int D, int DD>
inline vect<vect<T, D>, DD> xpose(vect<vect<T, DD>, D> const &a) {
  vect<vect<T, D>, DD> r;
  for (int dd = 0; dd < DD; ++dd)
    for (int d = 0; d < D; ++d)
      r[dd][d] = a[d][dd];
  return r;
}

/** Return the element-wise integer power of two vectors.  */
template <typename T, int D>
inline vect<T, D> ipow(const vect<T, D> &a, const vect<int, D> &b) {
  vect<T, D> r;
  for (int d = 0; d < D; ++d)
    r[d] = ipow(a[d], b[d]);
  return r;
}

DECLARE_FUNCTION_1(abs)
DECLARE_FUNCTION_1(ceil)
DECLARE_FUNCTION_1(fabs)
DECLARE_FUNCTION_1(floor)
DECLARE_FUNCTION_1(rint)
DECLARE_FUNCTION_1(round)
DECLARE_FUNCTION_1(sqrt)
DECLARE_FUNCTION_1(trunc)

DECLARE_FUNCTION_1_RET(lrint, int)

namespace std {
namespace Cactus {
DECLARE_FUNCTION_1_RET(good_fpclassify, int)
DECLARE_FUNCTION_1_RET(good_isfinite, int)
DECLARE_FUNCTION_1_RET(good_isinf, int)
DECLARE_FUNCTION_1_RET(good_isnan, int)
DECLARE_FUNCTION_1_RET(good_isnormal, int)
}
}

DECLARE_OPERATOR_1_RET(operator!, !, bool)

DECLARE_FUNCTION_2(max)
DECLARE_FUNCTION_2(min)
DECLARE_FUNCTION_2(pow)
DECLARE_FUNCTION_2(idiv)
DECLARE_FUNCTION_2(imod)

DECLARE_OPERATOR_2(operator+, +)
DECLARE_OPERATOR_2(operator-, -)
DECLARE_OPERATOR_2(operator*, *)
DECLARE_OPERATOR_2(operator/, /)
DECLARE_OPERATOR_2(operator%, %)
DECLARE_OPERATOR_2(operator&, &)
DECLARE_OPERATOR_2(operator|, |)
DECLARE_OPERATOR_2(operator^, ^)

DECLARE_OPERATOR_2_RET(operator&&, &&, bool)
DECLARE_OPERATOR_2_RET(operator||, ||, bool)
DECLARE_OPERATOR_2_RET(operator==, ==, bool)
DECLARE_OPERATOR_2_RET(operator!=, !=, bool)
DECLARE_OPERATOR_2_RET(operator<, <, bool)
DECLARE_OPERATOR_2_RET(operator<=, <=, bool)
DECLARE_OPERATOR_2_RET(operator>, >, bool)
DECLARE_OPERATOR_2_RET(operator>=, >=, bool)

// Reduction operators

// Identity
#define ID(x) (x)

DECLARE_REDUCTION_OPERATOR_1_T_RET(all, true, &=, ID, bool, bool)
DECLARE_REDUCTION_OPERATOR_1_T_RET(any, false, |=, ID, bool, bool)

DECLARE_REDUCTION_FUNCTION_1(maxval, a[0], max, ID)
DECLARE_REDUCTION_FUNCTION_1(minval, a[0], min, ID)
DECLARE_REDUCTION_OPERATOR_1(prod, 1, *=, ID)
DECLARE_REDUCTION_OPERATOR_1(sum, 0, +=, ID)

DECLARE_REDUCTION_OPERATOR_2(dot, 0, +=, *, ID)
DECLARE_REDUCTION_OPERATOR_2(hypot, 0, +=, *, sqrt)

#undef ID

/** Count the number of true (non-zero) elements in the vector.  */
template <typename T, int D> inline int count(const vect<T, D> &a) {
  return sum(ivect(a != T(0)));
}

/** Return the size (number of elements) of the vector.  */
template <typename T, int D> inline int size(const vect<T, D> &a) { return D; }

/** Return the first element.  */
template <typename T, int D> inline T head(const vect<T, D> &a) { return a[0]; }

/** Return all but the first element.  */
template <typename T, int D> inline vect<T, D - 1> tail(const vect<T, D> &a) {
  vect<T, D - 1> r;
  for (int d = 0; d < D - 1; ++d)
    r[d] = a[d + 1];
  return r;
}

/** Return the last element.  */
template <typename T, int D> inline T last(const vect<T, D> &a) {
  return a[D - 1];
}

/** Return all but the last element.  */
template <typename T, int D> inline vect<T, D - 1> init(const vect<T, D> &a) {
  vect<T, D - 1> r;
  for (int d = 0; d < D - 1; ++d)
    r[d] = a[d];
  return r;
}

/** Return the index of the first maximum element.  */
template <typename T, int D> inline int maxloc(const vect<T, D> &a) {
  ASSERT_VECT(D > 0);
  int r(0);
  for (int d = 1; d < D; ++d)
    if (a[d] > a[r])
      r = d;
  return r;
}

/** Return the index of the last maximum element.  */
template <typename T, int D> inline int maxloc1(const vect<T, D> &a) {
  ASSERT_VECT(D > 0);
  int r(D - 1);
  for (int d = D - 2; d >= 0; --d)
    if (a[d] > a[r])
      r = d;
  return r;
}

/** Return the index of the first minimum element.  */
template <typename T, int D> inline int minloc(const vect<T, D> &a) {
  ASSERT_VECT(D > 0);
  int r(0);
  for (int d = 1; d < D; ++d)
    if (a[d] < a[r])
      r = d;
  return r;
}

/** Return the index of the last minimum element.  */
template <typename T, int D> inline int minloc1(const vect<T, D> &a) {
  ASSERT_VECT(D > 0);
  int r(D - 1);
  for (int d = D - 2; d >= 0; --d)
    if (a[d] < a[r])
      r = d;
  return r;
}

/** Return the n-dimensional linear array index.  */
template <typename T, int D>
inline T index(const vect<T, D> &ash, const vect<T, D> &ind) {
  T r(0);
  for (int d = D - 1; d >= 0; --d) {
    ASSERT_VECT(ash[d] >= 0);
    // Be generous, and allow relative indices which may be negtive
    // ASSERT_VECT (ind[d]>=0 and ind[d]<ash[d]);
    ASSERT_VECT(abs(ind[d]) <= ash[d]);
    r = r * ash[d] + ind[d];
  }
  return r;
}

// Higher order functions

// This rarely used (and conflicts with Carpet::map), so it is renamed

/** Return a new vector where the function func() has been applied to
    all elements.  */
template <typename T, typename U, int D>
inline vect<U, D> vmap(U (*const func)(T x), const vect<T, D> &a);
template <typename T, typename U, int D>
inline vect<U, D> vmap(U (*const func)(T x), const vect<T, D> &a) {
  vect<U, D> r;
  for (int d = 0; d < D; ++d)
    r[d] = func(a[d]);
  return r;
}

/** Return a new vector where the function func() has been used
    element-wise to combine a and b.  */
template <typename S, typename T, typename U, int D>
inline vect<U, D> vzip(U (*const func)(S x, T y), const vect<S, D> &a,
                       const vect<T, D> &b);
template <typename S, typename T, typename U, int D>
inline vect<U, D> vzip(U (*const func)(S x, T y), const vect<S, D> &a,
                       const vect<T, D> &b) {
  vect<U, D> r;
  for (int d = 0; d < D; ++d)
    r[d] = func(a[d], b[d]);
  return r;
}

/** Return a scalar where the function func() has been used to reduce
    the vector a, starting with the scalar value val.  */
template <typename T, typename U, int D>
inline U fold(U (*const func)(U val, T x), U val, const vect<T, D> &a);
template <typename T, typename U, int D>
inline U fold(U (*const func)(U val, T x), U val, const vect<T, D> &a) {
  for (int d = 0; d < D; ++d)
    val = func(val, a[d]);
  return val;
}

/** Return a scalar where the function func() has been used to reduce
    the vector a, starting with element 0.  */
template <typename T, typename U, int D>
inline U fold1(U (*const func)(U val, T x), const vect<T, D> &a);
template <typename T, typename U, int D>
inline U fold1(U (*const func)(U val, T x), const vect<T, D> &a) {
  ASSERT_VECT(D >= 1);
  U val = a[0];
  for (int d = 1; d < D; ++d)
    val = func(val, a[d]);
  return val;
}

/** Return a vector where the function func() has been used to scan
    the vector a, starting with the scalar value val.  */
template <typename T, typename U, int D>
inline vect<U, D> scan0(U (*const func)(U val, T x), U val,
                        const vect<T, D> &a);
template <typename T, typename U, int D>
inline vect<U, D> scan0(U (*const func)(U val, T x), U val,
                        const vect<T, D> &a) {
  vect<U, D> r;
  for (int d = 0; d < D; ++d) {
    r[d] = val;
    val = func(val, a[d]);
  }
  return r;
}

/** Return a vector where the function func() has been used to scan
    the vector a, starting with element 0.  */
template <typename T, typename U, int D>
inline vect<U, D> scan1(U (*const func)(U val, T x), U val,
                        const vect<T, D> &a);
template <typename T, typename U, int D>
inline vect<U, D> scan1(U (*const func)(U val, T x), U val,
                        const vect<T, D> &a) {
  vect<U, D> r;
  for (int d = 0; d < D; ++d) {
    val = func(val, a[d]);
    r[d] = val;
  }
  return r;
}

// Memory usage

template <typename T, int D> inline size_t memoryof(vect<T, D> const &a);
template <typename T, int D> inline size_t memoryof(vect<T, D> const &a) {
  return a.memory();
}

// Input

/** Read a formatted vector from a stream.  */
template <typename T, int D>
inline istream &operator>>(istream &is, vect<T, D> &a) {
  a.input(is);
  return is;
}

// Output

/** Write a vector formatted to a stream.  */
template <typename T, int D>
inline ostream &operator<<(ostream &os, const vect<T, D> &a) {
  a.output(os);
  return os;
}

// Comparison

namespace std {
// ==
template <typename T, int D>
struct equal_to<vect<T, D> > : binary_function<vect<T, D>, vect<T, D>, bool> {
  bool operator()(const vect<T, D> &x, const vect<T, D> &y) const {
    /*const*/ equal_to<T> T_equal_to;
    for (int d = 0; d < D; ++d) {
      if (not T_equal_to(x[d], y[d]))
        return false;
    }
    return true;
  }
};

// <
template <typename T, int D>
struct less<vect<T, D> > : binary_function<vect<T, D>, vect<T, D>, bool> {
  bool operator()(const vect<T, D> &x, const vect<T, D> &y) const {
    /*const*/ less<T> T_less;
    for (int d = D - 1; d >= 0; --d) {
      if (T_less(x[d], y[d]))
        return true;
      if (T_less(y[d], x[d]))
        return false;
    }
    return false;
  }
};

// >
template <typename T, int D>
struct greater<vect<T, D> > : binary_function<vect<T, D>, vect<T, D>, bool> {
  bool operator()(const vect<T, D> &x, const vect<T, D> &y) const {
    return less<vect<T, D> >()(y, x);
  }
};

// >=
template <typename T, int D>
struct greater_equal<vect<T, D> >
    : binary_function<vect<T, D>, vect<T, D>, bool> {
  bool operator()(const vect<T, D> &x, const vect<T, D> &y) const {
    return not less<vect<T, D> >()(x, y);
  }
};

// <=
template <typename T, int D>
struct less_equal<vect<T, D> > : binary_function<vect<T, D>, vect<T, D>, bool> {
  bool operator()(const vect<T, D> &x, const vect<T, D> &y) const {
    return not greater<vect<T, D> >()(x, y);
  }
};

// !=
template <typename T, int D>
struct not_equal_to<vect<T, D> >
    : binary_function<vect<T, D>, vect<T, D>, bool> {
  bool operator()(const vect<T, D> &x, const vect<T, D> &y) const {
    return not equal_to<vect<T, D> >()(x, y);
  }
};
}

// MPI

template <typename T, int D>
inline MPI_Datatype mpi_datatype(vect<T, D> const &a) {
  return a.mpi_datatype();
}
namespace dist {
template <> inline MPI_Datatype mpi_datatype<ivect>() {
  ivect dummy;
  return mpi_datatype(dummy);
}
}

#if 0  
// Specialise explicit constructors

/** Constructor for 2-element vectors from 2 elements.  */
template<typename T>
inline vect<T,2>::vect<T,2> (const T& x, const T& y);
template<typename T>
inline vect<T,2>::vect<T,2> (const T& x, const T& y)
{
  elt[0]=x; elt[1]=y;
}

/** Constructor for 3-element vectors from 3 elements.  */
vect (const T& x, const T& y, const T& z);
vect (const T& x, const T& y, const T& z)
{
  ASSERT_VECT (D==3);
  elt[0]=x; elt[1]=y; elt[2]=z;
}

/** Constructor for 4-element vectors from 4 elements.  */
vect (const T& x, const T& y, const T& z, const T& t);
vect (const T& x, const T& y, const T& z, const T& t)
{
  ASSERT_VECT (D==4);
  elt[0]=x; elt[1]=y; elt[2]=z; elt[3]=t;
}
#endif

////////////////////////////////////////////////////////////////////////////////

// Specialise some constructors for lower dimensions
// These functions are declared, but must not be used.

template <> vect<int, 0>::vect(const int &x, const int &y);
template <> vect<int, 1>::vect(const int &x, const int &y);
template <> vect<int, 3>::vect(const int &x, const int &y);
template <> vect<int, 4>::vect(const int &x, const int &y);

template <> vect<int, 0>::vect(const int &x, const int &y, const int &z);
template <> vect<int, 1>::vect(const int &x, const int &y, const int &z);
template <> vect<int, 2>::vect(const int &x, const int &y, const int &z);
template <> vect<int, 4>::vect(const int &x, const int &y, const int &z);

template <>
vect<int, 0>::vect(const int &x, const int &y, const int &z, const int &t);
template <>
vect<int, 1>::vect(const int &x, const int &y, const int &z, const int &t);
template <>
vect<int, 2>::vect(const int &x, const int &y, const int &z, const int &t);
template <>
vect<int, 3>::vect(const int &x, const int &y, const int &z, const int &t);

// Specialise for CCTK_REAL

template <>
inline vect<CCTK_REAL, dim> &vect<CCTK_REAL, dim>::
operator%=(const vect<CCTK_REAL, dim> &a) {
  for (int d = 0; d < dim; ++d) {
    elt[d] = fmod(elt[d], a[d]);
    if (elt[d] > a[d] * (CCTK_REAL)(1.0 - 1.0e-10))
      elt[d] = (CCTK_REAL)0;
    if (elt[d] < a[d] * (CCTK_REAL)(1.0e-10))
      elt[d] = (CCTK_REAL)0;
  }
  return *this;
}

template <>
inline vect<CCTK_REAL, dim> operator%(const vect<CCTK_REAL, dim> &a,
                                      const vect<CCTK_REAL, dim> &b) {
  vect<CCTK_REAL, dim> r;
  for (int d = 0; d < dim; ++d) {
    r[d] = fmod(a[d], b[d]);
    if (r[d] > b[d] * (CCTK_REAL)(1.0 - 1.0e-10))
      r[d] = (CCTK_REAL)0;
    if (r[d] < b[d] * (CCTK_REAL)(1.0e-10))
      r[d] = (CCTK_REAL)0;
  }
  return r;
}

template <>
inline vect<CCTK_REAL, dim> idiv(const vect<CCTK_REAL, dim> &a,
                                 const vect<CCTK_REAL, dim> &b) {
  vect<CCTK_REAL, dim> r;
  for (int d = 0; d < dim; ++d) {
    r[d] = floor(a[d] / b[d]);
    if (r[d] > b[d] * (CCTK_REAL)(1.0 - 1.0e-10))
      r[d] = (CCTK_REAL)0;
    if (r[d] < b[d] * (CCTK_REAL)(1.0e-10))
      r[d] = (CCTK_REAL)0;
  }
  return r;
}

template <>
inline vect<CCTK_REAL, dim> imod(const vect<CCTK_REAL, dim> &a,
                                 const vect<CCTK_REAL, dim> &b) {
  vect<CCTK_REAL, dim> r;
  for (int d = 0; d < dim; ++d) {
    r[d] = a[d] / b[d];
    r[d] = b[d] * (r[d] - floor(r[d]));
    if (r[d] > b[d] * (CCTK_REAL)(1.0 - 1.0e-10))
      r[d] = (CCTK_REAL)0;
    if (r[d] < b[d] * (CCTK_REAL)(1.0e-10))
      r[d] = (CCTK_REAL)0;
  }
  return r;
}

#endif // VECT_HH
