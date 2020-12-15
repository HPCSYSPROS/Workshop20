#ifndef DEFS_HH
#define DEFS_HH

#include <cctk.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <list>
#include <map>
#include <memory>
#include <set>
#include <stack>
#include <vector>

#include "typeprops.hh"

// Disable bboxset2 if C++11 is not supported
#if !defined(HAVE_CCTK_CXX_AUTO_SPECIFIER) ||                                  \
    !defined(HAVE_CCTK_CXX_LAMBDA) || !defined(HAVE_CCTK_CXX_RANGE_BASED_FOR)
#ifndef CARPET_DISABLE_BBOXSET2
#define CARPET_WARN_DISABLE_BBOXSET2
#endif
#undef CARPET_DISABLE_BBOXSET2
#define CARPET_DISABLE_BBOXSET2
#endif

#ifndef CARPET_DISABLE_BBOXSET2
#define CARPET_ENABLE_BBOXSET2
#define CARPET_USE_BBOXSET2
#endif

using namespace std;

// TODO: autoconf these

#ifdef CARPET_USE_BOOST_FOREACH
#include <boost/foreach.hpp>
// We call the macro "forall", not "foreach", since the name "foreach"
// is taken by Boost and cannot be used here. An alternative name
// would be "foreach_" with a trailing underscore.
#define forall BOOST_FOREACH
// #  define forall_r BOOST_REVERSE_FOREACH
#else
#define forall(var, expr) for (var : expr)
#endif

#ifdef CARPET_USE_BOOST_SHARED_PTR
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#define shared_ptr boost::shared_ptr
#define make_shared boost::make_shared
#endif

// Stringify
#define STRINGIFY1(x) #x
#define STRINGIFY(x) STRINGIFY1(x)

// Structure member offsets
#undef offsetof
#define offsetof(TYPE, MEMBER) ((size_t) & ((TYPE *)0)->MEMBER)
#undef __offsetof__
#define __offsetof__ offsetof

// Number of dimensions
#ifndef CARPET_DIM
#define CARPET_DIM 3
#endif
const int dim = CARPET_DIM;

// A type for counting grid points (should be at least 64 bits)
typedef long long int size_type;

// Begin a new line without flushing the output buffer
char const *const eol = "\n";

// Check a return value
#define check(_expr)                                                           \
  do {                                                                         \
    bool const _val = (_expr);                                                 \
    assert(_val);                                                              \
  } while (0)

// Use this macro AT instead of vector's operator[] or at().
// Depending on the macro CARPET_OPTIMISE, this macro AT either checks
// for valid indices or not.
#if !defined(CARPET_OPTIMISE)
#define AT(index) at(index)
#else
#define AT(index) operator[](index)
#endif

// Some shortcuts for type names
template <typename T, int D> class vect;
template <typename T, int D> class bbox;
namespace bboxset1 {
template <typename T, int D> class bboxset;
}
#ifdef CARPET_ENABLE_BBOXSET2
namespace bboxset2 {
template <typename T, int D> class bboxset;
}
#endif
template <typename T, int D, typename P> class fulltree;

typedef vect<bool, dim> bvect;
typedef vect<int, dim> ivect;
typedef vect<CCTK_INT, dim> jvect;
typedef vect<CCTK_REAL, dim> rvect;
typedef bbox<int, dim> ibbox;
typedef bbox<CCTK_INT, dim> jbbox;
typedef bbox<CCTK_REAL, dim> rbbox;
namespace bboxset1 {
typedef bboxset<int, dim> ibset;
}
#ifdef CARPET_ENABLE_BBOXSET2
namespace bboxset2 {
typedef bboxset<int, dim> ibset;
}
#endif

// (Try to replace these by b2vect and i2vect)
typedef vect<vect<bool, 2>, dim> bbvect;
typedef vect<vect<int, 2>, dim> iivect;
typedef vect<vect<CCTK_INT, 2>, dim> jjvect;

typedef vect<vect<bool, dim>, 2> b2vect;
typedef vect<vect<int, dim>, 2> i2vect;
typedef vect<vect<CCTK_INT, dim>, 2> j2vect;
typedef vect<vect<CCTK_REAL, dim>, 2> r2vect;

struct pseudoregion_t;
struct region_t;

typedef fulltree<int, dim, pseudoregion_t> ipfulltree;

// A general type
enum centering { error_centered, vertex_centered, cell_centered };

// Divide, rounding to minus infinity
template <typename T> inline T idiv(T const x, T const y) {
  // round down manually if the result is negative
  return (x ^ y) >= T(0) ? x / y : (x - y + 1) / y;
}
// Modulo, rounding to minus infinity
template <typename T> inline T imod(T const x, T const y) {
  // return x - idiv(x,y)*y;
  return (x ^ y) >= T(0) ? x % y : (x - y + 1) % y + y - 1;
}

template <typename T> inline T div_down(T const x, T const y) {
  assert(x >= 0);
  assert(y > 0);
  return x / y;
}

template <typename T> inline T div_up(T const x, T const y) {
  assert(x >= 0);
  assert(y > 0);
  return (x + y - 1) / y;
}

template <typename T> inline T align_down(T const x, T const align) {
  assert(x >= 0);
  assert(align > 0);
  return div_down(x, align) * align;
}

template <typename T> inline T align_up(T const x, T const align) {
  assert(x >= 0);
  assert(align > 0);
  return div_up(x, align) * align;
}

// Useful helpers
template <class T> inline T square(const T x) { return x * x; }

template <class T> T ipow(T x, int y) CCTK_ATTRIBUTE_CONST;

template <class T> int ilog(T b, T x) CCTK_ATTRIBUTE_CONST;

// Access to CarpetLib parameters
CCTK_INT get_poison_value() CCTK_ATTRIBUTE_PURE;
CCTK_INT get_deadbeef() CCTK_ATTRIBUTE_PURE;

// Input streams
struct input_error {};
void skipws(istream &is);
void expect(istream &is, char c);
void consume(istream &is, char c);
void consume(istream &is, char const *c);

// Names for types

#ifdef HAVE_CCTK_INT1
inline const char *typestring(const CCTK_INT1 &) { return "CCTK_INT1"; }
#endif

#ifdef HAVE_CCTK_INT2
inline const char *typestring(const CCTK_INT2 &) { return "CCTK_INT2"; }
#endif

#ifdef HAVE_CCTK_INT4
inline const char *typestring(const CCTK_INT4 &) { return "CCTK_INT4"; }
#endif

#ifdef HAVE_CCTK_INT8
inline const char *typestring(const CCTK_INT8 &) { return "CCTK_INT8"; }
#endif

#ifdef HAVE_CCTK_INT16
inline const char *typestring(const CCTK_INT16 &) { return "CCTK_INT16"; }
#endif

#ifdef HAVE_CCTK_REAL4
inline const char *typestring(const CCTK_REAL4 &) { return "CCTK_REAL4"; }
#endif

#ifdef HAVE_CCTK_REAL8
inline const char *typestring(const CCTK_REAL8 &) { return "CCTK_REAL8"; }
#endif

#ifdef HAVE_CCTK_REAL16
inline const char *typestring(const CCTK_REAL16 &) { return "CCTK_REAL16"; }
#endif

#ifdef HAVE_CCTK_REAL4
inline const char *typestring(const CCTK_COMPLEX8 &) { return "CCTK_COMPLEX8"; }
#endif

#ifdef HAVE_CCTK_REAL8
inline const char *typestring(const CCTK_COMPLEX16 &) {
  return "CCTK_COMPLEX16";
}
#endif

#ifdef HAVE_CCTK_REAL16
inline const char *typestring(const CCTK_COMPLEX32 &) {
  return "CCTK_COMPLEX32";
}
#endif

// Provide implementations for some functions for complex numbers

#define IMPLEMENT_FUNCTIONS(T)                                                 \
                                                                               \
  inline int good_isfinite(T const &x) {                                       \
    return isfinite(x.real()) and isfinite(x.imag());                          \
  }                                                                            \
                                                                               \
  inline int good_isinf(T const &x) {                                          \
    return isinf(x.real()) or isinf(x.imag());                                 \
  }                                                                            \
                                                                               \
  inline int good_isnan(T const &x) {                                          \
    return isnan(x.real()) or isnan(x.imag());                                 \
  }                                                                            \
                                                                               \
  inline int good_isnormal(T const &x) {                                       \
    return isnormal(x.real()) and isnormal(x.imag());                          \
  }

namespace std {
namespace Cactus {

#ifdef HAVE_CCTK_COMPLEX8
IMPLEMENT_FUNCTIONS(CCTK_COMPLEX8)
#endif
#ifdef HAVE_CCTK_COMPLEX16
IMPLEMENT_FUNCTIONS(CCTK_COMPLEX16)
#endif
#ifdef HAVE_CCTK_COMPLEX32
IMPLEMENT_FUNCTIONS(CCTK_COMPLEX32)
#endif
}
}

#undef IMPLEMENT_FUNCTIONS

// Container memory usage
inline size_t memoryof(char const &e) { return sizeof e; }
inline size_t memoryof(short const &e) { return sizeof e; }
inline size_t memoryof(int const &e) { return sizeof e; }
inline size_t memoryof(long const &e) { return sizeof e; }
inline size_t memoryof(long long const &e) { return sizeof e; }
inline size_t memoryof(unsigned char const &e) { return sizeof e; }
inline size_t memoryof(unsigned short const &e) { return sizeof e; }
inline size_t memoryof(unsigned int const &e) { return sizeof e; }
inline size_t memoryof(unsigned long const &e) { return sizeof e; }
inline size_t memoryof(unsigned long long const &e) { return sizeof e; }
inline size_t memoryof(float const &e) { return sizeof e; }
inline size_t memoryof(double const &e) { return sizeof e; }
inline size_t memoryof(long double const &e) { return sizeof e; }
inline size_t memoryof(void const *const &e) { return sizeof e; }
template <class T> inline size_t memoryof(T const *const &e) {
  return sizeof e;
}
template <class T> inline size_t memoryof(typename list<T>::iterator const &i) {
  return sizeof i;
}
template <class T>
inline size_t memoryof(typename list<T>::const_iterator const &i) {
  return sizeof i;
}

template <class T> size_t memoryof(list<T> const &c) CCTK_ATTRIBUTE_PURE;
template <class T> size_t memoryof(set<T> const &c) CCTK_ATTRIBUTE_PURE;
template <class S, class T>
size_t memoryof(map<S, T> const &c) CCTK_ATTRIBUTE_PURE;
template <class T> size_t memoryof(stack<T> const &c) CCTK_ATTRIBUTE_PURE;
template <class T> size_t memoryof(vector<T> const &c) CCTK_ATTRIBUTE_PURE;

// Container input
template <class T> istream &input(istream &is, list<T> &l);
template <class T> istream &input(istream &is, set<T> &s);
template <class T> istream &input(istream &is, vector<T> &v);

template <class T> inline istream &operator>>(istream &is, list<T> &l) {
  return input(is, l);
}

template <class T> inline istream &operator>>(istream &is, set<T> &s) {
  return input(is, s);
}

template <class T> inline istream &operator>>(istream &is, vector<T> &v) {
  return input(is, v);
}

// Container output
template <class T> ostream &output(ostream &os, const list<T> &l);
template <class S, class T> ostream &output(ostream &os, const map<S, T> &m);
template <class S, class T> ostream &output(ostream &os, const pair<S, T> &p);
template <class T> ostream &output(ostream &os, const set<T> &s);
#ifdef CARPET_ENABLE_BBOXSET2
template <class T> ostream &output(ostream &os, const shared_ptr<T> &s);
#endif
template <class T> ostream &output(ostream &os, const stack<T> &s);
template <class T> ostream &output(ostream &os, const vector<T> &v);

template <class T> inline ostream &operator<<(ostream &os, const list<T> &l) {
  return output(os, l);
}

template <class S, class T>
inline ostream &operator<<(ostream &os, const map<S, T> &m) {
  return output(os, m);
}

template <class S, class T>
inline ostream &operator<<(ostream &os, const pair<S, T> &s) {
  return output(os, s);
}

template <class T> inline ostream &operator<<(ostream &os, const set<T> &s) {
  return output(os, s);
}

#ifdef CARPET_ENABLE_BBOXSET2
template <class T>
inline ostream &operator<<(ostream &os, const shared_ptr<T> &s) {
  return output(os, s);
}
#endif

template <class T> inline ostream &operator<<(ostream &os, const stack<T> &s) {
  return output(os, s);
}

template <class T> inline ostream &operator<<(ostream &os, const vector<T> &v) {
  return output(os, v);
}

#endif // DEFS_HH
