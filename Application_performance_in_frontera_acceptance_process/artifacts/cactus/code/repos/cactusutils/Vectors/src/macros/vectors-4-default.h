// Fallback vectorisation implementation: Do not vectorise

// We use macros here, so that we are not surprised by compilers which
// don't like to inline functions. This should also make debug builds
// (which may not inline) more efficient.



#include <assert.h>
#include <math.h>



#define vec4_architecture "scalar (no vectorisation, 32-bit precision)"

// Use CCTK_REAL4
#define CCTK_REAL4_VEC CCTK_REAL4

// Number of vector elements in a vector
#define CCTK_REAL4_VEC_SIZE 1

// Integer and boolean types corresponding to this real type
#define CCTK_INTEGER4     CCTK_REAL4
#define CCTK_BOOLEAN4     CCTK_REAL4
#define CCTK_INTEGER4_VEC CCTK_REAL4_VEC
#define CCTK_BOOLEAN4_VEC CCTK_REAL4_VEC



// Create a vector replicating a scalar
#define vec4_set1(a) (a)
// Create a vector from N scalars
#define vec4_set(a) (a)

// Access vectors elements
#define vec4_elt0(x) (x)
#define vec4_elt(x,d) (x)



// Load an aligned vector from memory
#define vec4_load(p) (p)
// Load an unaligned vector from memory
#define vec4_loadu(p) (p)

// Load a vector from memory that may or may not be aligned, as
// decided by the offset and the vector size. These functions are
// useful e.g. for loading neightbouring grid points while evaluating
// finite differencing stencils.
#define vec4_loadu_maybe(off,p) (p)
#define vec4_loadu_maybe3(off1,off2,off3,p) (p)

// Aligned store
#define vec4_store(p,x) ((p)=(x))
#define vec4_storeu(p,x) ((p)=(x))

// Unaligned store
#define vec4_store_nta(p,x) ((p)=(x))

#define vec4_store_partial_prepare(i,imin,imax) (0)
#define vec4_store_nta_partial(p,x) (vec4_store_nta(p,x))
// Store the n lower elements of a vector to memory
#define vec4_store_nta_partial_lo(p,x,n) (assert(0))
// Store the n higher elements of a vector into memory. This stores
// the vector elements into memory locations as if element 0 were
// stored at p.
#define vec4_store_nta_partial_hi(p,x,n) (assert(0))
#define vec4_store_nta_partial_mid(p,x,nlo,nhi) (assert(0))



// Operators
#define k4neg(x) (-(x))

#define k4add(x,y) ((x)+(y))
#define k4sub(x,y) ((x)-(y))
#define k4mul(x,y) ((x)*(y))
#define k4div(x,y) ((x)/(y))

// Fused multiply-add, defined as [+-]x*y[+-]z
#define k4madd(x,y,z)  (+(x)*(y)+(z))
#define k4msub(x,y,z)  (+(x)*(y)-(z))
#define k4nmadd(x,y,z) (-(x)*(y)-(z))
#define k4nmsub(x,y,z) (-(x)*(y)+(z))

// Functions
#define k4acos(x)       (acosf(x))
#define k4acosh(x)      (acoshf(x))
#define k4asin(x)       (asinf(x))
#define k4asinh(x)      (asinhf(x))
#define k4atan(x)       (atanf(x))
#define k4atan2(x,y)    (atan2f(x,y))
#define k4atanh(x)      (atanhf(x))
#define k4copysign(x,y) (copysign(x,y))
#define k4cos(x)        (cosf(x))
#define k4cosh(x)       (coshf(x))
#define k4exp(x)        (expf(x))
#define k4fabs(x)       (fabsf(x))
#define k4fmax(x,y)     (fmaxf(x,y))
#define k4fmin(x,y)     (fminf(x,y))
#define k4fnabs(x)      (-fabsf(x))
#define k4log(x)        (logf(x))
#define k4pow(x,a)      (powf(x,a))
#define k4sin(x)        (sinf(x))
#define k4sinh(x)       (sinhf(x))
#define k4sqrt(x)       (sqrtf(x))
#define k4tan(x)        (tanf(x))
#define k4tanh(x)       (tanhf(x))

#define k4sgn(x_)                                                       \
  ({                                                                    \
    CCTK_REAL x__=(x_);                                                 \
    CCTK_REAL x=x__;                                                    \
    x==(CCTK_REAL)0.0 ? (CCTK_REAL)0.0 : std::copysign((CCTK_REAL)1.0, x); \
  })
#define k4signbit(x) (std::signbit(x))

#define k4l2r(x_) ({ CCTK_INT4 x__=(x_); CCTK_INT4 x=x__; *(CCTK_REAL4*)&x; })
#define k4r2l(x_) ({ CCTK_REAL4 x__=(x_); CCTK_REAL4 x=x__; *(CCTK_INT4*)&x; })
#define k4lfalse k4l2r(0)
#define k4ltrue  k4l2r(1)
#define k4lnot(x)   k4l2r(!k4r2l(x))
#define k4land(x,y) k4l2r(k4r2l(x) && k4r2l(y))
#define k4lor(x,y)  k4l2r(k4r2l(x) || k4r2l(y))
#define k4lxor(x,y) k4l2r(!k4r2l(x) != !k4r2l(y))

#define k4ifthen(x,y,z) (k4r2l(x)?(y):(z))

#define k4cmpeq(x,y) k4l2r((x)==(y))
#define k4cmpne(x,y) k4l2r((x)!=(y))
#define k4cmpgt(x,y) k4l2r((x)>(y))
#define k4cmpge(x,y) k4l2r((x)>=(y))
#define k4cmplt(x,y) k4l2r((x)<(y))
#define k4cmple(x,y) k4l2r((x)<=(y))
