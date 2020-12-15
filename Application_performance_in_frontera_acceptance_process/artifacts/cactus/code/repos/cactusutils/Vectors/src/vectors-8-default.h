// -*-C++-*-
// Fallback vectorisation implementation: Do not vectorise



#include <cassert>
#include <cmath>



#define vec8_architecture "scalar (no vectorisation, 64-bit precision)"

// Use CCTK_REAL8
#define CCTK_REAL8_VEC CCTK_REAL8

// Number of vector elements in a vector
#define CCTK_REAL8_VEC_SIZE 1

vec_static_assert(sizeof(CCTK_REAL8_VEC) ==
                  sizeof(CCTK_REAL8) * CCTK_REAL8_VEC_SIZE);

// Integer and boolean types corresponding to this real type
#define CCTK_INTEGER8     CCTK_INT8
#define CCTK_BOOLEAN8     CCTK_INT8
#define CCTK_INTEGER8_VEC CCTK_INT8
#define CCTK_BOOLEAN8_VEC CCTK_INT8



// Create a vector replicating a scalar
#define vec8_set1(a) ((CCTK_REAL8)(a))
// Create a vector from N scalars
#define vec8_set(a) ((CCTK_REAL8)(a))

// Access vectors elements
#define vec8_elt0(x) (x)
#define vec8_elt(x,d) (x)
#define vec8_elti(x,d) (x)
#define vec8_eltb(x,d) (x)



// Load an aligned vector from memory
#define vec8_load(p) (p)
// Load an unaligned vector from memory
#define vec8_loadu(p) (p)

// Load a vector from memory that may or may not be aligned, as
// decided by the offset and the vector size. These functions are
// useful e.g. for loading neightbouring grid points while evaluating
// finite differencing stencils.
#define vec8_loadu_maybe(off,p) (p)
#define vec8_loadu_maybe3(off1,off2,off3,p) (p)

// Aligned store
#define vec8_store(p,x) ((p)=(x))
// Unaligned store
#define vec8_store_nta(p,x) ((p)=(x))

#define vec8_store_partial_prepare(i,imin,imax) ((void)0)
#define vec8_store_nta_partial(p,x) (vec8_store_nta(p,x))
// Store the n lower elements of a vector to memory
#define vec8_store_nta_partial_lo(p,x,n) (CCTK_BUILTIN_UNREACHABLE())
// Store the n higher elements of a vector into memory. This stores
// the vector elements into memory locations as if element 0 were
// stored at p.
#define vec8_store_nta_partial_hi(p,x,n) (CCTK_BUILTIN_UNREACHABLE())
#define vec8_store_nta_partial_mid(p,x,nlo,nhi) (CCTK_BUILTIN_UNREACHABLE())



// Operators
#define k8neg(x) (-(x))

#define k8add(x,y) ((x)+(y))
#define k8sub(x,y) ((x)-(y))
#define k8mul(x,y) ((x)*(y))
#define k8div(x,y) ((x)/(y))

// Fused multiply-add, defined as [+-]x*y[+-]z
#define k8madd(x,y,z)  (+(x)*(y)+(z))
#define k8msub(x,y,z)  (+(x)*(y)-(z))
#define k8nmadd(x,y,z) (-(x)*(y)-(z))
#define k8nmsub(x,y,z) (-(x)*(y)+(z))

// Functions
#define k8acos(x)       (acos(x))
#define k8acosh(x)      (acosh(x))
#define k8asin(x)       (asin(x))
#define k8asinh(x)      (asinh(x))
#define k8atan(x)       (atan(x))
#define k8atan2(x,y)    (atan2(x,y))
#define k8atanh(x)      (atanh(x))
#define k8copysign(x,y) (std::copysign(x,y))
#define k8cos(x)        (cos(x))
#define k8cosh(x)       (cosh(x))
#define k8exp(x)        (exp(x))
#define k8fabs(x)       (fabs(x))
#define k8fmax(x,y)     (fmax(x,y))
#define k8fmin(x,y)     (fmin(x,y))
#define k8fmod(x,y)     (fmod(x,y))
#define k8fnabs(x)      (-fabs(x))
#define k8log(x)        (log(x))
#define k8pow(x,a)      (pow(x,a))
#define k8sin(x)        (sin(x))
#define k8sinh(x)       (sinh(x))
#define k8sqrt(x)       (sqrt(x))
#define k8tan(x)        (tan(x))
#define k8tanh(x)       (tanh(x))
#define k8signbit(x)    (std::signbit(x))

#define k8lfalse 0
#define k8ltrue  1
#define k8lnot(x)   (!(x))
#define k8land(x,y) ((x) && (y))
#define k8lor(x,y)  ((x) || (y))
#define k8lxor(x,y) (!(x) != !(y))

#define k8ifthen(x,y,z) ((x)?(y):(z))

#define k8cmpeq(x,y) ((x)==(y))
#define k8cmpne(x,y) ((x)!=(y))
#define k8cmpgt(x,y) ((x)>(y))
#define k8cmpge(x,y) ((x)>=(y))
#define k8cmplt(x,y) ((x)<(y))
#define k8cmple(x,y) ((x)<=(y))

static inline CCTK_REAL8_VEC k8sgn(CCTK_REAL8_VEC const x)
{
  return x==(CCTK_REAL8)0.0 ? (CCTK_REAL8)0.0 : k8copysign((CCTK_REAL8)1.0, x);
}
