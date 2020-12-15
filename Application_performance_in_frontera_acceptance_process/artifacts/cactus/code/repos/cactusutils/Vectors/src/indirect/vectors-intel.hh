#include <cmath>
#include <cstdlib>
using namespace std;



#if defined(__SSE__)            // SSE (Intel)

#include <xmmintrin.h>

template<>
struct vec_T<float> {
  typedef float scalar_t;
  typedef __m128 impl_t;
  impl_t v;
  
  static inline size_t size()
  {
    return sizeof(impl_t)/sizeof(scalar_t);
  }
  
  inline vec_t ()
  {
  }
  inline vec_t (scalar_t const& a)
    : v(_mm_set1_ps(a))
  {
  }
  inline vec_t (scalar_t const& a0, scalar_t const& a1,
                scalar_t const& a2, scalar_t const& a3)
    : v(_mm_set_ps(a3,a2,a1,a0)) // reverse order!
  {
  }
  
  inline vec_t (impl_t const& w)
    : v(w)
  {
  }
  inline operator impl_t ()
  {
    return v;
  }
  
private:
  static inline scalar_t elt0 (impl_t const& v)
  {
    return _mm_cvtss_f32(v);    // this is a no-op
  }
public:
  inline scalar_t operator[] (size_t const d)
  {
    switch (d) {
    case 0: return elt0(v);
    case 1: return elt0(_mm_shuffle_ps(v,v,_MM_SHUFFLE(1,0,3,2)));
    case 2: return elt0(_mm_unpackhi_ps(v,v));
    case 3: return elt0(_mm_shuffle_ps(v,v,_MM_SHUFFLE(3,2,1,0)));
    }
  }
  
  static inline vec_t load (scalar_t const& a)
  {
    return _mm_load_ps(&a);
  }
  static inline vec_t loadu (scalar_t const& a)
  {
    return _mm_loadu_ps(&a);
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
    _mm_store_ps(&p,v);
  }
  inline void storeu (scalar_t& p) const
  {
    _mm_storeu_ps(&p,v);
  }
  inline void store_nta (scalar_t& p) const
  {
    _mm_stream_ps(&p,v);
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
  return _mm_and_ps(v,_mm_set1_pi32(0x7fffffffU));
}
template<typename T>
vec_t<T> fmax (vec_t<T> const& x, vec_t<T> const& y)
{
  return _mm_max_ps(x.v, y.v);
}
template<typename T>
vec_t<T> fmin (vec_t<T> const& x, vec_t<T> const& y)
{
  return _mm_min_ps(x.v, y.v);
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
  return _mm_sqrt_ps(x.v);
}

#endif



#if defined(__SSE2__)            // SSE2 (Intel)

#include <emmintrin.h>

template<>
struct vec_T<double> {
  typedef double scalar_t;
  typedef __m128d impl_t;
  impl_t v;
  
  static inline size_t size()
  {
    return sizeof(impl_t)/sizeof(scalar_t);
  }
  
  inline vec_t ()
  {
  }
  inline vec_t (scalar_t const& a)
    : v(_mm_set1_pd(a))
  {
  }
  inline vec_t (scalar_t const& a0, scalar_t const& a1)
    : v(_mm_set_pd(a1,a0))      // reverse order!
  {
  }
  
  inline vec_t (impl_t const& w)
    : v(w)
  {
  }
  inline operator impl_t ()
  {
    return v;
  }
  
private:
  static inline scalar_t elt0 (impl_t const& v)
  {
    return _mm_cvtss_f64(v);    // this is a no-op
  }
public:
  inline scalar_t operator[] (size_t const d)
  {
    switch (d) {
    case 0: return elt0(v);
    case 1: return elt0(_mm_unpackhi_pd(v,v));
    }
  }
  
  static inline vec_t load (scalar_t const& a)
  {
    return _mm_load_pd(&a);
  }
  static inline vec_t loadu (scalar_t const& a)
  {
    return _mm_loadu_pd(&a);
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
    _mm_store_pd(&p,v);
  }
  inline void storeu (scalar_t& p) const
  {
    _mm_storeu_pd(&p,v);
  }
  inline void store_nta (scalar_t& p) const
  {
    _mm_stream_pd(&p,v);
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
  return _mm_and_pd(v,_mm_set1_epi64(0x7fffffffffffffffULL));
}
template<typename T>
vec_t<T> fmax (vec_t<T> const& x, vec_t<T> const& y)
{
  return _mm_max_pd(x.v, y.v);
}
template<typename T>
vec_t<T> fmin (vec_t<T> const& x, vec_t<T> const& y)
{
  return _mm_min_pd(x.v, y.v);
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
  return _mm_sqrt_pd(x.v);
}

#endif
