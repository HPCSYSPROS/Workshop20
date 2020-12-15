#include <cctk.h>

#include <cassert>
#include <iostream>
#include <typeinfo>

#include "defs.hh"
#include "bboxset.hh"

#include "vect.hh"

using namespace std;

// Input
template <typename T, int D> void vect<T, D>::input(istream &is) {
  skipws(is);
  consume(is, '[');
  for (int d = 0; d < D; ++d) {
    is >> (*this)[d];
    assert(is.good());
    if (d < D - 1) {
      skipws(is);
      consume(is, ',');
    }
  }
  skipws(is);
  consume(is, ']');
}

// Output
template <typename T, int D> void vect<T, D>::output(ostream &os) const {
  os << "[";
  for (int d = 0; d < D; ++d) {
    os << (*this)[d];
    if (d < D - 1)
      os << ",";
  }
  os << "]";
}

// MPI
template <typename T, int D> MPI_Datatype vect<T, D>::mpi_datatype() const {
  static bool initialised = false;
  static MPI_Datatype newtype;
  if (not initialised) {
    vect<T, D> const &s = *this;
#define ENTRY(type, name)                                                      \
  {                                                                            \
    sizeof s.name / sizeof(type), /* count elements */                         \
        (const char *) & s.name - (const char *) &                             \
            s,                      /* offsetof doesn't work (why?) */         \
        dist::mpi_datatype<type>(), /* find MPI datatype */                    \
        STRINGIFY(name),            /* field name */                           \
        STRINGIFY(type),            /* type name */                            \
  }
    dist::mpi_struct_descr_t const descr[] = {
        ENTRY(T, elt), {1, sizeof s, MPI_UB, "MPI_UB", "MPI_UB"}};
#undef ENTRY
    ostringstream buf;
    buf << "vect<" << typeid(T).name() << "," << D << ">";
    newtype = dist::create_mpi_datatype(sizeof descr / sizeof descr[0], descr,
                                        buf.str().c_str(), sizeof s);
    initialised = true;
  }
  return newtype;
}

// Specialise some constructors for lower dimensions
// These functions are declared, but must not be used.

// NOTE: __builtin_unreachable() triggers a bug in the OSX linker; we
// therefore use assert(0) instead.

template <> vect<int, 0>::vect(const int &x, const int &y) { assert(0); }
template <> vect<int, 1>::vect(const int &x, const int &y) { assert(0); }
template <> vect<int, 3>::vect(const int &x, const int &y) { assert(0); }
template <> vect<int, 4>::vect(const int &x, const int &y) { assert(0); }

template <> vect<int, 0>::vect(const int &x, const int &y, const int &z) {
  assert(0);
}
template <> vect<int, 1>::vect(const int &x, const int &y, const int &z) {
  assert(0);
}
template <> vect<int, 2>::vect(const int &x, const int &y, const int &z) {
  assert(0);
}
template <> vect<int, 4>::vect(const int &x, const int &y, const int &z) {
  assert(0);
}

template <>
vect<int, 0>::vect(const int &x, const int &y, const int &z, const int &t) {
  assert(0);
}
template <>
vect<int, 1>::vect(const int &x, const int &y, const int &z, const int &t) {
  assert(0);
}
template <>
vect<int, 2>::vect(const int &x, const int &y, const int &z, const int &t) {
  assert(0);
}
template <>
vect<int, 3>::vect(const int &x, const int &y, const int &z, const int &t) {
  assert(0);
}

// Note: We need all dimensions all the time.
template class vect<int, 0>;
template class vect<int, 1>;
template class vect<int, 2>;
template class vect<int, 3>;
template class vect<int, 4>;

template void vect<unsigned long long, dim>::input(istream &is);
template void vect<CCTK_REAL, dim>::input(istream &is);
template void vect<vect<bool, 2>, dim>::input(istream &is);
template void vect<vect<bool, dim>, 2>::input(istream &is);
template void vect<vect<int, dim>, 2>::input(istream &is);

template void vect<bool, 2>::output(ostream &os) const;
template void vect<bool, dim>::output(ostream &os) const;
template void vect<unsigned long long, dim>::output(ostream &os) const;
template void vect<CCTK_REAL, 2>::output(ostream &os) const;
template void vect<CCTK_REAL, dim>::output(ostream &os) const;
template void vect<vect<bool, 2>, dim>::output(ostream &os) const;
template void vect<vect<int, 2>, dim>::output(ostream &os) const;
template void vect<vect<long, 2>, dim>::output(ostream &os) const;
template void vect<vect<bool, dim>, 2>::output(ostream &os) const;
template void vect<vect<int, dim>, 2>::output(ostream &os) const;
template void vect<vect<CCTK_REAL, dim>, 2>::output(ostream &os) const;
template void vect<vect<vector<int>, 2>, dim>::output(ostream &os) const;

// Instantiate for bboxset class

#define DEFINE_FAKE_VECT_OPERATIONS(T, D)                                      \
  template <> vect<T, D> vect<T, D>::dir(const int d) { assert(0); }           \
  template <> vect<T, D> vect<T, D>::seq() { assert(0); }                      \
  template <> vect<T, D> vect<T, D>::seq(const int n) { assert(0); }           \
  template <> vect<T, D> vect<T, D>::seq(const int n, const int s) {           \
    assert(0);                                                                 \
  }                                                                            \
  template <> vect<T, D> &vect<T, D>::operator*=(const vect<T, D> &) {         \
    assert(0);                                                                 \
  }                                                                            \
  template <> vect<T, D> &vect<T, D>::operator*=(const T &) { assert(0); }     \
  template <> vect<T, D> &vect<T, D>::operator/=(const vect<T, D> &) {         \
    assert(0);                                                                 \
  }                                                                            \
  template <> vect<T, D> &vect<T, D>::operator/=(const T &) { assert(0); }     \
  template <> vect<T, D> &vect<T, D>::operator%=(const vect<T, D> &) {         \
    assert(0);                                                                 \
  }                                                                            \
  template <> vect<T, D> &vect<T, D>::operator%=(const T &) { assert(0); }     \
  template <> vect<T, D> &vect<T, D>::operator^=(const vect<T, D> &) {         \
    assert(0);                                                                 \
  }                                                                            \
  template <> vect<T, D> &vect<T, D>::operator^=(const T &) { assert(0); }     \
  template <> vect<T, D> vect<T, D>::operator+() const { assert(0); }          \
  template <> vect<T, D> vect<T, D>::operator-() const { assert(0); }          \
  template <> vect<T, D> vect<T, D>::operator~() const { assert(0); }          \
  template class vect<T, D>;                                                   \
  template size_t memoryof(const vect<T, D> &);                                \
  template istream &operator>>(istream &is, vect<T, D> &);                     \
  template ostream &operator<<(ostream &os, const vect<T, D> &);

// need typedefs since macro calls can't have commas in arguments
typedef bboxset<int, dim> T1;
typedef vect<bboxset<int, dim>, 2> T2;

DEFINE_FAKE_VECT_OPERATIONS(T1, dim)
DEFINE_FAKE_VECT_OPERATIONS(T2, dim)

#undef DEFINE_FAKE_VECT_OPERATIONS
