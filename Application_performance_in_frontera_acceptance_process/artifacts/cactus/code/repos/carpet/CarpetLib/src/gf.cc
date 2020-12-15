#include <cctk.h>

#include <cassert>

#include "defs.hh"

#include "gf.hh"

using namespace std;

// Constructors
template <typename T>
gf<T>::gf(const int varindex_, const operator_type transport_operator_, th &t_,
          dh &d_, const int prolongation_order_time_, const int vectorlength_,
          const int vectorindex_, ggf *const vectorleader_)
    : ggf(varindex_, transport_operator_, t_, d_, prolongation_order_time_,
          vectorlength_, vectorindex_, vectorleader_) {}

// Destructors
template <typename T> gf<T>::~gf() {}

// Memory usage
template <typename T> size_t gf<T>::memory() const { return ggf::memory(); }

// Output
template <typename T> ostream &gf<T>::output(ostream &os) const {
  T Tdummy;
  os << "gf<" << typestring(Tdummy) << ">:" << varindex << "["
     << CCTK_VarName(varindex) << "],"
     << "tls=" << timelevels_;
  return os;
}

#define TYPECASE(N, T) template class gf<T>;

#include "typecase.hh"

#undef TYPECASE
