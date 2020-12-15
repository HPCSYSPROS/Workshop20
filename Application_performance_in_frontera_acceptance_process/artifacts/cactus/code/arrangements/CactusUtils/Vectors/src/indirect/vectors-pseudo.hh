#include <cmath>
#include <cstdlib>
using namespace std;

template<typename T>
struct vec_t {
  size_t const D = 2;
  typedef T scalar_t;
  typedef T impl_t;
  impl_t v[D];
  
  static inline size_t size()
  {
    return D;
  }
  
  inline vec_t ()
  {
  }
  inline vec_t (scalar_t const& a)
  {
    for (size_t d=0; d<D; ++d) v[d]=a;
  }
  
  inline operator impl_t ()
  {
    return v;
  }
  
  inline scalar_t operator[] (size_t const d) const
  {
    return v[d];
  }
  
  static inline vec_t load (scalar_t const& a)
  {
    vec_t r;
    for (size_t d=0; d<D; ++d) r.v[d]=(&a)[d];
    return r;
  }
  static inline vec_t loadu (scalar_t const& a)
  {
    return load(a);
  }
  // Load a vector from memory that may or may not be aligned, as
  // decided by the offset and the vector size
  static inline vec_t loadu_maybe (int const off, scalar_t const& p)
  {
    return load(a);
  }
  static inline vec_t loadu_maybe3 (int const off0, int const off1,
                                    int const off2,
                                    scalar_t const& p)
  {
    return load(a);
  }
  inline void store (scalar_t& p) const
  {
    for (size_t d=0; d<D; ++d) (&p)[d]=v[d];
  }
  inline void storeu (scalar_t& p) const
  {
    store(p);
  }
  inline void store_nta (scalar_t& p) const
  {
    store(p);
  }
  inline void store_nta_partial_lo (scalar_t& p, size_t const cnt) const
  {
    for (size_t d=0; d<cnt; ++d) (&p)[d]=v[d];
  }
  inline void store_nta_partial_hi (scalar_t& p, size_t const cnt) const
  {
    for (size_t d=D-cnt; d<D; ++d) (&p)[d]=v[d];
  }
  
  inline vec_t operator+ () const
  {
    vec_t r=*this;
    for (size_t d=0; d<D; ++d) r.v[d]=+v[d];
    return r;
  }
  inline vec_t operator- () const
  {
    vec_t r=*this;
    for (size_t d=0; d<D; ++d) r.v[d]=-v[d];
    return r;
  }
  inline vec_t operator+ (vec_t const& x) const
  {
    vec_t r=*this;
    return r+=x;
  }
  inline vec_t operator- (vec_t const& x) const
  {
    vec_t r=*this;
    return r-=x;
  }
  inline vec_t operator* (vec_t const& x) const
  {
    vec_t r=*this;
    return r*=x;
  }
  inline vec_t operator/ (vec_t const& x) const
  {
    vec_t r=*this;
    return r/=x;
  }
  inline vec_t& operator+= (vec_t const& x)
  {
    for (size_t d=0; d<D; ++d) v[d]+=x.v[d];
    return *this;
  }
  inline vec_t& operator-= (vec_t const& x)
  {
    for (size_t d=0; d<D; ++d) v[d]-=x.v[d];
    return *this;
  }
  inline vec_t& operator*= (vec_t const& x)
  {
    for (size_t d=0; d<D; ++d) v[d]*=x.v[d];
    return *this;
  }
  inline vec_t& operator/= (vec_t const& x)
  {
    for (size_t d=0; d<D; ++d) v[d]/=x.v[d];
    return *this;
  }
};

template<typename T>
vec_t<T> exp (vec_t<T> const& x)
{
  vec_t r;
  for (size_t d=0; d<D; ++d) r.v[d]=exp(x.v[d]);
  return r;
}
template<typename T>
vec_t<T> fabs (vec_t<T> const& x)
{
  vec_t r;
  for (size_t d=0; d<D; ++d) r.v[d]=fabs(x.v[d]);
  return r;
}
template<typename T>
vec_t<T> fmax (vec_t<T> const& x, vec_t<T> const& y)
{
  vec_t r;
  for (size_t d=0; d<D; ++d) r.v[d]=fmax(x.v[d],y.v[d]);
  return r;
}
template<typename T>
vec_t<T> fmin (vec_t<T> const& x, vec_t<T> const& y)
{
  vec_t r;
  for (size_t d=0; d<D; ++d) r.v[d]=fmin(x.v[d],y.v[d]);
  return r;
}
template<typename T>
vec_t<T> ifthen (bool const b, vec_t<T> const& x, vec_t<T> const& y)
{
  return b ? x : y;
}
vec_t<T> log (vec_t<T> const& x)
{
  vec_t r;
  for (size_t d=0; d<D; ++d) r.v[d]=log(x.v[d]);
  return r;
}
template<typename T>
vec_t<T> pow (vec_t<T> const& x, typename vec_t<T>::scalar_t const& a)
{
  vec_t r;
  for (size_t d=0; d<D; ++d) r.v[d]=pow(x.v[d],a);
  return r;
}
vec_t<T> sqrt (vec_t<T> const& x)
{
  vec_t r;
  for (size_t d=0; d<D; ++d) r.v[d]=sqrt(x.v[d]);
  return r;
}
