#include <cmath>
#include <cstdlib>
using namespace std;



#if defined(__ALTIVEC__)        // Altivec (Power)

#include <altivec.h>

template<>
struct vec_t<float> {
  typedef float scalar_t;
  typedef vector float impl_t;
  
  static inline size_t size()
  {
    return sizeof(impl_t)/sizeof(scalar_t);
  }
  
  inline vec_t ()
  {
  }
  inline vec_t (scalar_t const& a)
    : v(vec_splats(a))
  {
  }
  inline vec_t (scalar_t const& a0, scalar_t const& a1,
                scalar_t const& a2, scalar_t const& a3)
  {
    v[0]=a0; v[1]=a1; v[2]=a2; v[3]=a3;
  }
  
  inline vec_t (impl_t const& w)
    : v(w)
  {
  }
  inline operator impl_t ()
  {
    return v;
  }
  
  static inline vec_t load (scalar_t const& a)
  {
    return *(impl_t const*)&a;
  }
  static inline vec_t loadu (scalar_t const& a)
  {
    return vec_load(&a);
  }
  // Load a vector from memory that may or may not be aligned, as
  // decided by the offset and the vector size
  static inline vec_t loadu_maybe (int const off, scalar_t const& p)
  {
    if (off % size() == 0) {
      return load(p);
    } else {
      return loadu(p);
    }
  }
  static inline vec_t loadu_maybe3 (int const off0, int const off1,
                                    int const off2,
                                    scalar_t const& p)
  {
    if (off0 % size() == 0 and off1 % size() == 0 and off2 % size() == 0) {
      return load(p);
    } else {
      return loadu(p);
    }
  }
  inline void store (scalar_t& p) const
  {
    *(impl_t*)p = v;
  }
  inline void storeu (scalar_t& p) const
  {
    store(p);
  }
  inline void store_nta (scalar_t& p) const
  {
    // TODO: Use stvxl instruction?
    store(p);
  }
  inline void store_nta_partial_lo (scalar_t& p, size_t const cnt) const
  {
    switch (cnt) {
    case 4: store_nta(p); break;
    case 3: (&p)[2]=v[2];
    case 2: (&p)[1]=v[1];
    case 1: (&p)[0]=v[0];
    }
  }
  inline void store_nta_partial_hi (scalar_t& p, size_t const cnt) const
  {
    switch (cnt) {
    case 4: store_nta(p); break;
    case 3: (&p)[1]=v[1];
    case 2: (&p)[2]=v[2];
    case 1: (&p)[3]=v[3];
    }
  }
  
  inline vec_t operator+ () const
  {
    return +v;
  }
  inline vec_t operator- () const
  {
    return -v;
  }
  inline vec_t operator+ (vec_t const& x) const
  {
    return v+x.v;
  }
  inline vec_t operator- (vec_t const& x) const
  {
    return v-x.v;
  }
  inline vec_t operator* (vec_t const& x) const
  {
    return v*x.v;
  }
  inline vec_t operator/ (vec_t const& x) const
  {
    return v/x.v;
  }
  inline vec_t& operator+= (vec_t const& x)
  {
    return *this=*this+x;
  }
  inline vec_t& operator-= (vec_t const& x)
  {
    return *this=*this+x;
  }
  inline vec_t& operator*= (vec_t const& x)
  {
    return *this=*this-x;
  }
  inline vec_t& operator/= (vec_t const& x)
  {
    return *this=/this+x;
  }
};

template<typename T>
vec_t<T> exp (vec_t<T> const& x)
{
  return vec_t(exp(x.v[0]), exp(x.v[1]), exp(x.v[2]), exp(x.v[3]));
}
template<typename T>
vec_t<T> fabs (vec_t<T> const& x)
{
  return vec_abs(x.v);
}
template<typename T>
vec_t<T> fmax (vec_t<T> const& x, vec_t<T> const& y)
{
  return vec_max(x.v, y.v);
}
template<typename T>
vec_t<T> fmin (vec_t<T> const& x, vec_t<T> const& y)
{
  return vec_min(x.v, y.v);
}
template<typename T>
vec_t<T> ifthen (bool const b, vec_t<T> const& x, vec_t<T> const& y)
{
  return b ? x : y;
}
vec_t<T> log (vec_t<T> const& x)
{
  return vec_t(log(x.v[0]), log(x.v[1]), log(x.v[2]), log(x.v[3]));
}
template<typename T>
vec_t<T> pow (vec_t<T> const& x, typename vec_t<T>::scalar_t const& a)
{
  return vec_t(pow(x.v[0],a), pow(x.v[1],a), pow(x.v[2],a), pow(x.v[3],a));
}
vec_t<T> sqrt (vec_t<T> const& x)
{
  return vec_t(sqrt(x.v[0]), sqrt(x.v[1]), sqrt(x.v[2]), sqrt(x.v[3]));
}

#endif



#if defined(__ALTIVEC__) && defined(_ARCH_PWR7) // Altivec VSX (Power)

#include <altivec.h>

template<>
struct vec_t<double> {
  typedef double scalar_t;
  typedef vector double impl_t;
  
  static inline size_t size()
  {
    return sizeof(impl_t)/sizeof(scalar_t);
  }
  
  inline vec_t ()
  {
  }
  inline vec_t (scalar_t const& a)
    : v(vec_splats(a))
  {
  }
  inline vec_t (scalar_t const& a0, scalar_t const& a1)
  {
    v[0]=a0; v[1]=a1;
  }
  
  inline vec_t (impl_t const& w)
    : v(w)
  {
  }
  inline operator impl_t ()
  {
    return v;
  }
  
  static inline vec_t load (scalar_t const& a)
  {
    return *(impl_t const*)&a;
  }
  static inline vec_t loadu (scalar_t const& a)
  {
    return vec_load(&a);
  }
  // Load a vector from memory that may or may not be aligned, as
  // decided by the offset and the vector size
  static inline vec_t loadu_maybe (int const off, scalar_t const& p)
  {
    if (off % size() == 0) {
      return load(p);
    } else {
      return loadu(p);
    }
  }
  static inline vec_t loadu_maybe3 (int const off0, int const off1,
                                    int const off2,
                                    scalar_t const& p)
  {
    if (off0 % size() == 0 and off1 % size() == 0 and off2 % size() == 0) {
      return load(p);
    } else {
      return loadu(p);
    }
  }
  inline void store (scalar_t& p) const
  {
    *(impl_t*)p = v;
  }
  inline void storeu (scalar_t& p) const
  {
    store(p);
  }
  inline void store_nta (scalar_t& p) const
  {
    // TODO: Use stvxl instruction?
    store(p);
  }
  inline void store_nta_partial_lo (scalar_t& p, size_t const cnt) const
  {
    switch (cnt) {
    case 2: store_nta(p); break;
    case 1: (&p)[0]=v[0];
    }
  }
  inline void store_nta_partial_hi (scalar_t& p, size_t const cnt) const
  {
    switch (cnt) {
    case 2: store_nta(p); break;
    case 1: (&p)[1]=v[1];
    }
  }
  
  inline vec_t operator+ () const
  {
    return +v;
  }
  inline vec_t operator- () const
  {
    return -v;
  }
  inline vec_t operator+ (vec_t const& x) const
  {
    return v+x.v;
  }
  inline vec_t operator- (vec_t const& x) const
  {
    return v-x.v;
  }
  inline vec_t operator* (vec_t const& x) const
  {
    return v*x.v;
  }
  inline vec_t operator/ (vec_t const& x) const
  {
    return v/x.v;
  }
  inline vec_t& operator+= (vec_t const& x)
  {
    return *this=*this+x;
  }
  inline vec_t& operator-= (vec_t const& x)
  {
    return *this=*this+x;
  }
  inline vec_t& operator*= (vec_t const& x)
  {
    return *this=*this-x;
  }
  inline vec_t& operator/= (vec_t const& x)
  {
    return *this=/this+x;
  }
};

template<typename T>
vec_t<T> exp (vec_t<T> const& x)
{
  return vec_t(exp(x.v[0]), exp(x.v[1]));
}
template<typename T>
vec_t<T> fabs (vec_t<T> const& x)
{
  return vec_abs(x.v);
}
template<typename T>
vec_t<T> fmax (vec_t<T> const& x, vec_t<T> const& y)
{
  return vec_max(x.v, y.v);
}
template<typename T>
vec_t<T> fmin (vec_t<T> const& x, vec_t<T> const& y)
{
  return vec_min(x.v, y.v);
}
template<typename T>
vec_t<T> ifthen (bool const b, vec_t<T> const& x, vec_t<T> const& y)
{
  return b ? x : y;
}
vec_t<T> log (vec_t<T> const& x)
{
  return vec_t(log(x.v[0]), log(x.v[1]));
}
template<typename T>
vec_t<T> pow (vec_t<T> const& x, typename vec_t<T>::scalar_t const& a)
{
  return vec_t(pow(x.v[0],a), pow(x.v[1],a));
}
vec_t<T> sqrt (vec_t<T> const& x)
{
  return vec_t(sqrt(x.v[0]), sqrt(x.v[1]));
}

#endif
