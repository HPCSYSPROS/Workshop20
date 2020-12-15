using namespace std;



// A class template that provides a vectorised type for the underlying
// scalar type T. This implementation does "nothing", i.e. just
// provides a "vector" class with a vector size of 1, forwarding all
// operations to the scalar type.
//
// This implementation uses small integers in several places, e.g. in
// the [] operator. This is efficient only if these integers are
// compile time constants, so that the compiler can remove the
// corresponding if and switch statements.
template<typename T>
struct vec_t {
  
  // Names for the underlying scalar type, and for the vector type
  // used to implement this class. For example, with SSE2, it would be
  // scalar_t=double, and impl_t=__m128d.
  typedef T scalar_t;
  typedef T impl_t;
  
  // The payload -- the actual vector content
  impl_t v;
  
  // Vector size (number of elements)
  static inline size_t size()
  {
    return sizeof(impl_t)/sizeof(scalar_t);
  }
  
  // Constructors
  inline vec_t ()
  {
  }
  inline vec_t (scalar_t const& a)
    : v(a)
  {
  }
  
  // Convert to the implementation vector type
  inline operator impl_t ()
  {
    return v;
  }
  
  // Access individual vector elements
  inline scalar_t operator[] (size_t const d) const
  {
    return v;
  }
  
  // Load vectors from memory. For convenience when using this class,
  // these accept references to the scalar type instead of pointers to
  // the vector type. These routines are static members of the class,
  // so that they can be used as VEC::load(p); if they were
  // stand-alone functions, one would have to write load<SCALAR>(p)
  // instead.
  //
  // Aligned load
  static inline vec_t load (scalar_t const& p)
  {
    return p;
  }
  // Unaligned load
  static inline vec_t loadu (scalar_t const& p)
  {
    return p;
  }
  // Load a vector from memory that may or may not be aligned, as
  // decided by the offset and the vector size. These functions are
  // useful e.g. for loading neightbouring grid points while
  // evaluating finite differencing stencils.
  static inline vec_t loadu_maybe (int const off, scalar_t const& p)
  {
    return p;
  }
  static inline vec_t loadu_maybe3 (int const off0, int const off1,
                                    int const off2,
                                    scalar_t const& p)
  {
    return p;
  }
};

// Store vectors to memory. These routines are stand-alone functions,
// so that they can be used as vec_store(p,x); if they were class
// members, one would have to write x.store(p) instead, or possibly
// VEC::store(p,x).
//
// Aligned store
template<typename T>
inline void vec_store (typename vec_t<T>::scalar_t& p, vec_t<T> const& x)
{
  p=x.v;
}
// Unaligned store
template<typename T>
inline void vec_storeu (typename vec_t<T>::scalar_t& p, vec_t<T> const& x)
{
  p=x.v;
}
// Non-temporal store, i.e. a store that bypasses the cache
template<typename T>
inline void vec_store_nta (typename vec_t<T>::scalar_t& p, vec_t<T> const& x)
{
  p=x.v;
}
// Store the cnt lower elements of a vector, bypassing the cache if
// possible
template<typename T>
inline void vec_store_nta_partial_lo (typename vec_t<T>::scalar_t& p,
                                      vec_t<T> const& x,
                                      size_t const cnt)
{
  assert(0);
}
// Store the cnt higher elements of a vector, bypassing the cache if
// possible. This stores the vector elements into memory locations as
// if element 0 were stored at p.
template<typename T>
inline void vec_store_nta_partial_hi (typename vec_t<T>::scalar_t& p,
                                      vec_t<T> const& x,
                                      size_t const cnt)
{
  assert(0);
}

template<typename T>
inline vec_t<T> operator+ (vec_t<T> const& x)
{
  return +x.v;
}
template<typename T>
inline vec_t<T> operator- (vec_t<T> const& x)
{
  return -x.v;
}

template<typename T>
inline vec_t<T> operator+ (vec_t<T> const& x, vec_t<T> const& y)
{
  return x.v + y.v;
}
template<typename T>
inline vec_t<T> operator- (vec_t<T> const& x, vec_t<T> const& y)
{
  return x.v - y.v;
}
template<typename T>
inline vec_t<T> operator* (vec_t<T> const& x, vec_t<T> const& y)
{
  return x.v * y.v;
}
template<typename T>
inline vec_t<T> operator/ (vec_t<T> const& x, vec_t<T> const& y)
{
  return x.v / y.v;
}

template<typename T>
inline vec_t<T>& operator+= (vec_t<T>& x, vec_t<T> const& y)
{
  x.v += y.v;
  return x;
}
template<typename T>
inline vec_t<T>& operator-= (vec_t<T>& x, vec_t<T> const& y)
{
  x.v -= y.v;
  return x;
}
template<typename T>
inline vec_t<T>& operator*= (vec_t<T>& x, vec_t<T> const& y)
{
  x.v *= y.v;
  return x;
}
template<typename T>
inline vec_t<T>& operator/= (vec_t<T>& x, vec_t<T> const& y)
{
  x.v /= y.v;
  return x;
}

template<typename T>
inline vec_t<T> exp (vec_t<T> const& x)
{
  return exp(x.v);
}
template<typename T>
inline vec_t<T> fabs (vec_t<T> const& x)
{
  return fabs(x.v);
}
template<typename T>
inline vec_t<T> log (vec_t<T> const& x)
{
  return log(x.v);
}
template<typename T>
inline vec_t<T> sqrt (vec_t<T> const& x)
{
  return sqrt(x.v);
}

template<typename T>
inline vec_t<T> fmax (vec_t<T> const& x, vec_t<T> const& y)
{
  return fmax(x.v, y.v);
}
template<typename T>
inline vec_t<T> fmin (vec_t<T> const& x, vec_t<T> const& y)
{
  return fmin(x.v, y.v);
}
template<typename T>
inline vec_t<T> pow (vec_t<T> const& x, typename vec_t<T>::scalar_t const& a)
{
  return pow(x.v, a);
}
template<typename T>
inline vec_t<T> pow (vec_t<T> const& x, int const& i)
{
  return pow(x.v, i);
}
