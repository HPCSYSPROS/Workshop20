// -*-C++-*-
// Vectorise using IBM's Altivec (Power)

// Use the type vector double directly, without introducing a wrapper class
// Use macros instead of inline functions



#include <altivec.h>



#define vec4_architecture "Altivec"

// Vector type corresponding to CCTK_REAL
#define CCTK_REAL4_VEC vector float

// Number of vector elements in a CCTK_REAL_VEC
#define CCTK_REAL4_VEC_SIZE 4

// Integer and boolean types corresponding to this real type
#define CCTK_INTEGER4     CCTK_REAL4
#define CCTK_BOOLEAN4     CCTK_REAL4
#define CCTK_INTEGER4_VEC CCTK_REAL4_VEC
#define CCTK_BOOLEAN4_VEC CCTK_REAL4_VEC



// Create vectors, extract vector elements

#define vec4_set1(a)  (vec_splats(a))
#define vec4_set(a,b,c,d)                       \
  ({                                            \
    CCTK_REAL4_VEC x;                           \
    x[0]=(a);                                   \
    x[1]=(b);                                   \
    x[2]=(c);                                   \
    x[3]=(d);                                   \
    x;                                          \
  })

#define vec4_elt0(x) ((x)[0])
#define vec4_elt1(x) ((x)[1])
#define vec4_elt2(x) ((x)[2])
#define vec4_elt3(x) ((x)[3])
#define vec4_elt(x,d) ((x)[d])



// Load and store vectors

// Load a vector from memory (aligned and unaligned); this loads from
// a reference to a scalar
#define vec4_load(p)  (*(CCTK_REAL4_VEC const*)&(p))
#define vec4_loadu(p) (*(CCTK_REAL4_VEC const*)&(p))

// Load a vector from memory that may or may not be aligned, as
// decided by the offset and the vector size
#define vec4_loadu_maybe(off,p)             (vec4_loadu(p))
#define vec4_loadu_maybe3(off1,off2,off3,p) (vec4_loadu(p))

// Store a vector to memory (aligned and non-temporal); this stores to
// a reference to a scalar
#define vec4_store(p,x)  (*(CCTK_REAL4_VEC*)&(p)=(x))
#define vec4_storeu(p,x) (*(CCTK_REAL4_VEC*)&(p)=(x))
#if ! VECTORISE_STREAMING_STORES
#  define vec4_store_nta(p,x) (vec4_store(p,x))
#else
// use stvxl instruction
#  define vec4_store_nta(p,x) (vec_stl(x,0,(CCTK_REAL4_VEC*)&(p)))
#endif

// Store a lower or higher partial vector (aligned and non-temporal);
// the non-temporal hint is probably ignored
#define vec4_store_nta_partial_lo(p_,x_,n)      \
  ({                                            \
    CCTK_REAL4     const& p__=(p_);             \
    CCTK_REAL4_VEC const  x__=(x_);             \
    CCTK_REAL4     const& p=p__;                \
    CCTK_REAL4_VEC const  x=x__;                \
    switch (n) {                                \
    case 3: (&p)[2]=x[2];                       \
    case 2: (&p)[1]=x[1];                       \
    case 1: (&p)[0]=x[0];                       \
    }                                           \
  })
#define vec4_store_nta_partial_hi(p_,x_,n)      \
  ({                                            \
    CCTK_REAL4     const& p__=(p_);             \
    CCTK_REAL4_VEC const  x__=(x_);             \
    CCTK_REAL4     const& p=p__;                \
    CCTK_REAL4_VEC const  x=x__;                \
    switch (n) {                                \
    case 3: (&p)[1]=x[1];                       \
    case 2: (&p)[2]=x[2];                       \
    case 1: (&p)[3]=x[3];                       \
    }                                           \
  })
#define vec4_store_nta_partial_mid(p_,x_,nlo_,nhi_)     \
  ({                                                    \
    CCTK_REAL4     const& p__  =(p_);                   \
    CCTK_REAL4_VEC const  x__  =(x_);                   \
    int            const  nlo__=(nlo_);                 \
    int            const  nhi__=(nhi_);                 \
    CCTK_REAL4     const& p  =p__;                      \
    CCTK_REAL4_VEC const  x  =x__;                      \
    int            const  nlo=nlo__;                    \
    int            const  nhi=nhi__;                    \
    if (nlo==3 and nhi==3) {                            \
      (&p)[1]=x[1];                                     \
      (&p)[2]=x[2];                                     \
    } else if (nlo==2 and nhi==3) {                     \
      (&p)[1]=x[1];                                     \
    } else if (nlo==3 and nhi==2) {                     \
      (&p)[2]=x[2];                                     \
    }                                                   \
  })



// Functions and operators

// Operators
#define k4neg(x) (-(x))

#define k4add(x,y) ((x)+(y))
#define k4sub(x,y) ((x)-(y))
#define k4mul(x,y) ((x)*(y))
#define k4div(x,y) ((x)/(y))

// Fused multiply-add, defined as [+-]x*y[+-]z
#define k4madd(x,y,z)  (vec_madd(x,y,z))
#define k4msub(x,y,z)  (vec_msub(x,y,z))
#define k4nmadd(x,y,z) (vec_nmadd(x,y,z))
#define k4nmsub(x,y,z) (vec_nmsub(x,y,z))

// Cheap functions
#define k4fabs(x)   (vec_abs(x))
#define k4fmax(x,y) (vec_max(x,y))
#define k4fmin(x,y) (vec_min(x,y))
#define k4fnabs(x)  (vec_nabs(x))

// Expensive functions
#define K4REPL(f,x_)                            \
  ({                                            \
    CCTK_REAL4_VEC const x__=(x_);              \
    CCTK_REAL4_VEC const x=x__;                 \
    vec4_set(f(vec4_elt0(x)),                   \
             f(vec4_elt1(x)),                   \
             f(vec4_elt2(x)),                   \
             f(vec4_elt3(x)));                  \
  })
#define K4REPL2S(f,x_,a_)                       \
  ({                                            \
    CCTK_REAL4_VEC const x__=(x_);              \
    CCTK_REAL4     const a__=(a_);              \
    CCTK_REAL4_VEC const x=x__;                 \
    CCTK_REAL4     const a=a__;                 \
    vec4_set(f(vec4_elt0(x),a),                 \
             f(vec4_elt1(x),a),                 \
             f(vec4_elt2(x),a),                 \
             f(vec4_elt3(x),a));                \
  })
#define K4REPL2(f,x_,y_)                        \
  ({                                            \
    CCTK_REAL4_VEC const x__=(x_);              \
    CCTK_REAL4_VEC const y__=(y_);              \
    CCTK_REAL4_VEC const x=x__;                 \
    CCTK_REAL4_VEC const y=y__;                 \
    vec4_set(f(vec4_elt0(x),vec4_elt0(y)),      \
             f(vec4_elt1(x),vec4_elt1(y)),      \
             f(vec4_elt2(x),vec4_elt2(y)),      \
             f(vec4_elt3(x),vec4_elt3(y)));     \
  })

#define k4acos(x)    K4REPL(acosf,x)
#define k4acosh(x)   K4REPL(acoshf,x)
#define k4asin(x)    K4REPL(asinf,x)
#define k4asinh(x)   K4REPL(asinhf,x)
#define k4atan(x)    K4REPL(atanf,x)
#define k4atan2(x,y) K4REPL2(atan2f,x,y)
#define k4atanh(x)   K4REPL(atanhf,x)
#define k4cos(x)     K4REPL(cosf,x)
#define k4cosh(x)    K4REPL(coshf,x)
#define k4exp(x)     K4REPL(expf,x)
#define k4fmod(x,y)  K4REPL2(fmodf,x,y)
#define k4log(x)     K4REPL(logf,x)
#define k4pow(x,a)   K4REPL2S(powf,x,a)
#define k4sin(x)     K4REPL(sinf,x)
#define k4sinh(x)    K4REPL(sinhf,x)
#define k4tan(x)     K4REPL(tanf,x)
#define k4tanh(x)    K4REPL(tanhf,x)

#define k4ifmsb(x,y,z)                                                  \
  (vec_sel((z), (y), vec_sra(vec_convert((x), &(vector int*)0), 31)))
