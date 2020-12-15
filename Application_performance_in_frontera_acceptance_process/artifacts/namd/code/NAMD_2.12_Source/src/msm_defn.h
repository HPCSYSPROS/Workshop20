/* msm_defn.h */

#ifndef NL_MSM_DEFN_H
#define NL_MSM_DEFN_H

/* #define TEST_INLINING */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#ifdef WIN32
#define strcasecmp(s,t) stricmp(s,t)
#define strncasecmp(s,t,n) strnicmp(s,t,n)
#else
#include <strings.h>
#endif

#include "msm.h"

#include "wkfutils.h"  /* timers from John Stone */

/* avoid parameter name collisions with AIX5 "hz" macro */
#undef hz

/** Index vectors for x,y,z,q coordinates */
#undef  X
#define X  0
#undef  Y
#define Y  1
#undef  Z
#define Z  2
#undef  Q
#define Q  3

/** Return number of elements in static array */
#undef  NELEMS
#define NELEMS(arr)  (sizeof(arr)/sizeof(arr[0]))

/** Template to create a 3D grid to access data of given TYPE */
#undef GRID_TEMPLATE
#define GRID_TEMPLATE(TYPE) \
  typedef struct NL_Msmgrid_##TYPE##_t { \
    TYPE *buffer;       /* raw buffer */ \
    TYPE *data;         /* data access offset from buffer */ \
    size_t numbytes;    /* number of bytes in use by buffer */ \
    size_t maxbytes;    /* actual allocation for buffer */ \
    int i0, j0, k0;     /* starting index value for each dimension */ \
    int ni, nj, nk;     /* number of elements in each dimension */ \
  } NL_Msmgrid_##TYPE
  /* expect closing ';' */

/** Initialize grid to empty */
#define GRID_INIT(_p) \
  ((_p)->buffer=NULL, (_p)->data=NULL, (_p)->numbytes=0, (_p)->maxbytes=0, \
   (_p)->i0=0, (_p)->j0=0, (_p)->k0=0, (_p)->ni=0, (_p)->nj=0, (_p)->nk=0)
  /* expect closing ';' */

/** Finished with grid, free its memory */
#define GRID_DONE(_p) \
  free((_p)->buffer)
  /* expect closing ';' */

/** Determine the signed flattened index for 3D grid datum */
#define GRID_INDEX(_p, _i, _j, _k) \
  (((_k)*((_p)->nj) + (_j))*((_p)->ni) + (_i))
  /* expect closing ';' */

/** Obtain pointer to 3D grid datum */
#define GRID_POINTER(_p, _i, _j, _k) \
  ((_p)->data + GRID_INDEX(_p, _i, _j, _k))
  /* expect closing ';' */

/** Resize 3D grid buffer, setup its indexing.  Grab more memory
 * when needed (must be used within function returning int). */
#define GRID_RESIZE(_p, TYPE, __i0, __ni, __j0, __nj, __k0, __nk) \
  do { \
    int _i0=(__i0), _ni=(__ni); \
    int _j0=(__j0), _nj=(__nj); \
    int _k0=(__k0), _nk=(__nk); \
    size_t _numbytes = (_nk * _nj) * (size_t) _ni * sizeof((_p)->buffer[0]); \
    if ((_p)->maxbytes < _numbytes) { \
      void *_t = realloc((_p)->buffer, _numbytes); \
      if (NULL == _t) return NL_MSM_ERROR_MALLOC; \
      (_p)->buffer = (TYPE *) _t; \
      (_p)->maxbytes = _numbytes; \
    } \
    (_p)->numbytes = _numbytes; \
    (_p)->i0 = _i0, (_p)->ni = _ni; \
    (_p)->j0 = _j0, (_p)->nj = _nj; \
    (_p)->k0 = _k0, (_p)->nk = _nk; \
    (_p)->data = (_p)->buffer + GRID_INDEX((_p), -_i0, -_j0, -_k0); \
  } while (0)
  /* expect closing ';' */

/* reset 3D grid data to 0 */
#define GRID_ZERO(_p) \
  memset((_p)->buffer, 0, (_p)->numbytes)  /* ; */

/* check 3D grid index range when debugging
 * (must be used within function returning int) */
#ifdef MSM_DEBUG
#define ASSERT(expr)
#define GRID_INDEX_CHECK(a, _i, _j, _k) \
  do { \
    ASSERT((a)->i0 <= (_i) && (_i) < (a)->ni + (a)->i0); \
    ASSERT((a)->j0 <= (_j) && (_j) < (a)->nj + (a)->j0); \
    ASSERT((a)->k0 <= (_k) && (_k) < (a)->nk + (a)->k0); \
  } while (0)
  /* expect closing ';' */
#else
#define ASSERT(expr)
#define GRID_INDEX_CHECK(a, _i, _j, _k)
#endif


/** Default MSM parameters */
#define DEFAULT_GRIDSPACING  2.5
#define DEFAULT_APPROX       NL_MSM_APPROX_CUBIC
#define DEFAULT_SPLIT        NL_MSM_SPLIT_TAYLOR2
#define DEFAULT_NLEVELS      0  /* set number of levels as needed */

#define DEFAULT_DENSITY      0.1
#define DEFAULT_BINFILL      0.8
#define DEFAULT_NBINSLOTS    8


/** SPOLY() calculates the polynomial part of the
 * normalized smoothing of 1/r.
 *
 * Returns g(R), where R=r/a, and (d/dR)g(R).
 *
 *   pg - float*, points to variable to receive g(R)
 *   pdg - float*, points to variable to receive (d/dR)g(R)
 *   ra - (r/a), assumed to be between 0 and 1
 *   split - identify the type of smoothing used to split the potential
 */
#undef  SPOLY
#define SPOLY(pg, pdg, ra, split) \
  do { \
    double _r = ra;  /* _r=(r/a) */ \
    double _s = _r*_r;  /* _s=(r/a)^2 */ \
    double _g = 0, _dg = 0; \
    ASSERT(0 <= _s && _s <= 1); \
    switch (split) { \
      /* from section 5.1 of thesis */ \
      case NL_MSM_SPLIT_TAYLOR2: \
        _g = 1 + (_s-1)*(-1./2 + (_s-1)*(3./8));                               \
        _dg = (2*_r)*(-1./2 + (_s-1)*(3./4)); \
        break;                                                                 \
      case NL_MSM_SPLIT_TAYLOR3: \
        _g = 1 + (_s-1)*(-1./2 + (_s-1)*(3./8 + (_s-1)*(-5./16)));             \
        _dg = (2*_r)*(-1./2 + (_s-1)*(3./4 + (_s-1)*(-15./16))); \
        break;                                                                 \
      case NL_MSM_SPLIT_TAYLOR4: \
        _g = 1 + (_s-1)*(-1./2 + (_s-1)*(3./8 + (_s-1)*(-5./16                 \
                + (_s-1)*(35./128))));                                         \
        _dg = (2*_r)*(-1./2 + (_s-1)*(3./4 + (_s-1)*(-15./16 \
                + (_s-1)*(35./32))));      \
        break;                                                                 \
      case NL_MSM_SPLIT_TAYLOR5: \
        _g = 1 + (_s-1)*(-1./2 + (_s-1)*(3./8 + (_s-1)*(-5./16                 \
                + (_s-1)*(35./128 + (_s-1)*(-63./256)))));                     \
        _dg = (2*_r)*(-1./2 + (_s-1)*(3./4 + (_s-1)*(-15./16 + (_s-1)*(35./32  \
                + (_s-1)*(-315./256)))));                                      \
        break;                                                                 \
      case NL_MSM_SPLIT_TAYLOR6: \
        _g = 1 + (_s-1)*(-1./2 + (_s-1)*(3./8 + (_s-1)*(-5./16                 \
                + (_s-1)*(35./128 + (_s-1)*(-63./256                           \
                    + (_s-1)*(231./1024))))));                                 \
        _dg = (2*_r)*(-1./2 + (_s-1)*(3./4 + (_s-1)*(-15./16 + (_s-1)*(35./32  \
                + (_s-1)*(-315./256 + (_s-1)*(693./512))))));                  \
        break;                                                                 \
      case NL_MSM_SPLIT_TAYLOR7: \
        _g = 1 + (_s-1)*(-1./2 + (_s-1)*(3./8 + (_s-1)*(-5./16                 \
                + (_s-1)*(35./128 + (_s-1)*(-63./256                           \
                    + (_s-1)*(231./1024 + (_s-1)*(-429./2048)))))));           \
        _dg = (2*_r)*(-1./2 + (_s-1)*(3./4 + (_s-1)*(-15./16 + (_s-1)*(35./32  \
                + (_s-1)*(-315./256 + (_s-1)*(693./512                         \
                    + (_s-1)*(-3003./2048)))))));                              \
        break;                                                                 \
      case NL_MSM_SPLIT_TAYLOR8: \
        _g = 1 + (_s-1)*(-1./2 + (_s-1)*(3./8 + (_s-1)*(-5./16                 \
                + (_s-1)*(35./128 + (_s-1)*(-63./256                           \
                    + (_s-1)*(231./1024 + (_s-1)*(-429./2048                   \
                        + (_s-1)*(6435./32768))))))));                         \
        _dg = (2*_r)*(-1./2 + (_s-1)*(3./4 + (_s-1)*(-15./16 + (_s-1)*(35./32  \
                + (_s-1)*(-315./256 + (_s-1)*(693./512                         \
                    + (_s-1)*(-3003./2048 + (_s-1)*(6435./4096))))))));        \
        break;                                                                 \
      case NL_MSM_SPLIT_TAYLOR1: \
        _g = 1 + (_s-1)*(-1./2); \
        _dg = (2*_r)*(-1./2); \
        break; \
      /* from section 5.2 of thesis */ \
      case NL_MSM_SPLIT_SIGMA2_3:  /* the "perfect" smoothing */ \
        _g = 2 + _s*(-2 + _r);                                                \
        _dg = _r*(-4 + _r*3);                                                  \
        break;                                                                 \
      case NL_MSM_SPLIT_SIGMA3_5: \
        _g = 9./4 + _s*(-5./2 + _s*(9./4 - _r));                             \
        _dg = _r*(-5 + _s*(9 + _r*(-5)));                                     \
        break;                                                                 \
      case NL_MSM_SPLIT_SIGMA4_6: \
        _g = 21./8 + _s*(-35./8 + _s*(63./8 + _r*(-7 + _r*(15./8))));        \
        _dg = _r*(-35./4 + _s*(63./2 + _r*(-35 + _r*(45./4))));               \
        break;                                                                 \
      case NL_MSM_SPLIT_SIGMA4_7: \
        _g = 5./2 + _s*(-7./2 + _s*(7./2 + _s*(-5./2 + _r)));               \
        _dg = _r*(-7 + _s*(14 + _s*(-15 + _r*(7))));                         \
        break;                                                                 \
      case NL_MSM_SPLIT_SIGMA5_8: \
        _g = 45./16 + _s*(-21./4 + _s*(63./8 + _s*(-45./4                   \
                + _r*(9 + _r*(-35./16)))));                                    \
        _dg = _r*(-21./2 + _s*(63./2 + _s*(-135./2                           \
                + _r*(63 + _r*(-35./2)))));                                    \
        break;                                                                 \
      case NL_MSM_SPLIT_SIGMA5_9: \
        _g = 175./64 + _s*(-75./16 + _s*(189./32 + _s*(-75./16              \
                + _s*(175./64 - _r))));                                       \
        _dg = _r*(-75./8 + _s*(189./8 + _s*(-225./8 + _s*(175./8            \
                  + _r*(-9)))));                                               \
        break;                                                                 \
      case NL_MSM_SPLIT_SIGMA6_9: \
        _g = 25./8 + _s*(-15./2 + _s*(63./4 + _s*(-75./2                    \
                + _r*(45 + _r*(-175./8 + _r*4)))));                            \
        _dg = _r*(-15 + _s*(63 + _s*(-225                                    \
                + _r*(315 + _r*(-175 + _r*36)))));                             \
        break;                                                                 \
      case NL_MSM_SPLIT_SIGMA6_10: \
        _g = 385./128 + _s*(-825./128 + _s*(693./64 + _s*(-825./64          \
                + _s*(1925./128 + _r*(-11 + _r*(315./128))))));               \
        _dg = _r*(-825./64 + _s*(693./16 + _s*(-2475./32                     \
                + _s*(1925./16 + _r*(-99 + _r*(1575./64))))));                \
        break;                                                                 \
      case NL_MSM_SPLIT_SIGMA6_11: \
        _g = 189./64 + _s*(-385./64 + _s*(297./32 + _s*(-297./32            \
                + _s*(385./64 + _s*(-189./64 + _r)))));                      \
        _dg = _r*(-385./32 + _s*(297./8 + _s*(-891./16 + _s*(385./8         \
                  + _s*(-945./32 + _r*(11))))));                              \
        break;                                                                 \
      case NL_MSM_SPLIT_SIGMA7_11: \
        _g = 105./32 + _s*(-275./32 + _s*(297./16 + _s*(-495./16            \
                + _s*(1925./32 + _r*(-66 + _r*(945./32 + _r*(-5)))))));       \
        _dg = _r*(-275./16 + _s*(297./4 + _s*(-1485./8                       \
                + _s*(1925./4 + _r*(-594 + _r*(4725./16 + _r*(-55)))))));     \
        break;                                                                 \
      case NL_MSM_SPLIT_SIGMA7_12: \
        _g = 819./256 + _s*(-1001./128 + _s*(3861./256                       \
              + _s*(-1287./64 + _s*(5005./256 + _s*(-2457./128              \
                    + _r*(13 + _r*(-693./256)))))));                           \
        _dg = _r*(-1001./64 + _s*(3861./64 + _s*(-3861./32                   \
                + _s*(5005./32 + _s*(-12285./64 + _r*(143                    \
                      + _r*(-2079./64)))))));                                  \
        break;                                                                 \
      case NL_MSM_SPLIT_SIGMA7_13: \
        _g = 1617./512 + _s*(-1911./256 + _s*(7007./512 + _s*(-2145./128    \
                + _s*(7007./512 + _s*(-1911./256 + _s*(1617./512 - _r))))));\
        _dg = _r*(-1911./128 + _s*(7007./128 + _s*(-6435./64 + _s*(7007./64 \
                  + _s*(-9555./128 + _s*(4851./128 + _r*(-13)))))));         \
        break;                                                                 \
      case NL_MSM_SPLIT_SIGMA8_12: \
        _g = 455./128 + _s*(-715./64 + _s*(3861./128 + _s*(-2145./32        \
                + _s*(25025./128 + _r*(-286 + _r*(12285./64 + _r*(-65         \
                        + _r*(1155./128))))))));                               \
        _dg = _r*(-715./32 + _s*(3861./32 + _s*(-6435./16                    \
                + _s*(25025./16 + _r*(-2574 + _r*(61425./32 + _r*(-715        \
                        + _r*(3465./32))))))));                                \
        break;                                                                 \
      case NL_MSM_SPLIT_SIGMA8_13: \
        _g = 441./128 + _s*(-637./64 + _s*(3003./128                         \
              + _s*(-1287./32 + _s*(7007./128 + _s*(-5733./64               \
                    + _r*(91 + _r*(-4851./128 + _r*(6))))))));                 \
        _dg = _r*(-637./32 + _s*(3003./32 + _s*(-3861./16                    \
                + _s*(7007./16 + _s*(-28665./32 + _r*(1001                   \
                      + _r*(-14553./32 + _r*(78))))))));                       \
        break;                                                                 \
      case NL_MSM_SPLIT_SIGMA8_14: \
        _g = 3465./1024 + _s*(-9555./1024 + _s*(21021./1024                  \
              + _s*(-32175./1024 + _s*(35035./1024 + _s*(-28665./1024       \
                    + _s*(24255./1024 + _r*(-15 + _r*(3003./1024))))))));     \
        _dg = _r*(-9555./512 + _s*(21021./256 + _s*(-96525./512              \
                + _s*(35035./128 + _s*(-143325./512 + _s*(72765./256        \
                      + _r*(-195 + _r*(21021./512))))))));                     \
        break;                                                                 \
      case NL_MSM_SPLIT_SIGMA8_15: \
        _g = 429./128 + _s*(-1155./128 + _s*(2457./128 + _s*(-3575./128     \
                + _s*(3575./128 + _s*(-2457./128 + _s*(1155./128            \
                      + _s*(-429./128 + _r)))))));                            \
        _dg = _r*(-1155./64 + _s*(2457./32 + _s*(-10725./64                  \
                + _s*(3575./16 + _s*(-12285./64 + _s*(3465./32              \
                      + _s*(-3003./64 + _r*(15))))))));                       \
        break;                                                                 \
      /* from section 5.3 of thesis */ \
      case NL_MSM_SPLIT_SIGMA2_6: \
        _g = (31./16) + _s*(-23./16 + _s*(9./16 + _s*(-1./16))); \
        _dg = (2*_r)*(-23./16 + _s*(9./8 + _s*(-3./16))); \
        break; \
      /* from section 5.4 of thesis */ \
      case NL_MSM_SPLIT_SWITCH1_2:                                             \
        if (_r > 1./2) {                                                       \
          _g = 5./3 + _r + _s*(-3 + _r*(4./3));                               \
          _dg = 1 + _r*(-6 + _r*(4));                                          \
        }                                                                      \
        else {                                                                 \
          _g = 11./6 - _s;                                                    \
          _dg = _r*(-2);                                                       \
        }                                                                      \
        break;                                                                 \
      case NL_MSM_SPLIT_SWITCH3_4:                                             \
        if (_r > 3./4) {                                                       \
          _g = 5./7 + _r*(27./7 + _r*(-41./7 + _r*(16./7)));                   \
          _dg = 27./7 + _r*(-82./7 + _r*(48./7));                              \
        }                                                                      \
        else {                                                                 \
          _g = 47./28 + _s*(-5./7);                                           \
          _dg = _r*(-10./7);                                                   \
        }                                                                      \
        break;                                                                 \
      case NL_MSM_SPLIT_SWITCH7_8:                                             \
        if (_r > 7./8) {                                                       \
          _g = -19./15 + _r*(49./5 + _r*(-59./5 + _r*(64./15)));               \
          _dg = 49./5 + _r*(-118./5 + _r*(64./5));                             \
        }                                                                      \
        else {                                                                 \
          _g = 191./120 + _s*(-3./5);                                         \
          _dg = _r*(-6./5);                                                    \
        }                                                                      \
        break;                                                                 \
      default: \
        return NL_MSM_ERROR_SUPPORT; \
    } \
    *(pg) = _g; \
    *(pdg) = _dg; \
  } while (0)
  /* closing ';' from use as function call */


/** SPOLY_SPREC() calculates the polynomial part of the
 * normalized smoothing of 1/r.  Single precision version.
 *
 * Returns g(R), where R=r/a, and (d/dR)g(R).
 *
 *   pg - float*, points to variable to receive g(R)
 *   pdg - float*, points to variable to receive (d/dR)g(R)
 *   ra - (r/a), assumed to be between 0 and 1
 *   split - identify the type of smoothing used to split the potential
 */
#undef  SPOLY_SPREC
#define SPOLY_SPREC(pg, pdg, ra, split) \
  do { \
    float _r = ra;  /* _r=(r/a) */ \
    float _s = _r*_r;  /* _s=(r/a)^2 */ \
    float _g = 0, _dg = 0; \
    ASSERT(0 <= _s && _s <= 1); \
    switch (split) { \
      /* from section 5.1 of thesis */ \
      case NL_MSM_SPLIT_TAYLOR2: \
        _g = 1 + (_s-1)*(-1.f/2 + (_s-1)*(3.f/8));                             \
        _dg = (2*_r)*(-1.f/2 + (_s-1)*(3.f/4)); \
        break;                                                                 \
      case NL_MSM_SPLIT_TAYLOR3: \
        _g = 1 + (_s-1)*(-1.f/2 + (_s-1)*(3.f/8 + (_s-1)*(-5.f/16)));          \
        _dg = (2*_r)*(-1.f/2 + (_s-1)*(3.f/4 + (_s-1)*(-15.f/16))); \
        break;                                                                 \
      case NL_MSM_SPLIT_TAYLOR4: \
        _g = 1 + (_s-1)*(-1.f/2 + (_s-1)*(3.f/8 + (_s-1)*(-5.f/16              \
                + (_s-1)*(35.f/128))));                                        \
        _dg = (2*_r)*(-1.f/2 + (_s-1)*(3.f/4 + (_s-1)*(-15.f/16 \
                + (_s-1)*(35.f/32))));      \
        break;                                                                 \
      case NL_MSM_SPLIT_TAYLOR5: \
        _g = 1 + (_s-1)*(-1.f/2 + (_s-1)*(3.f/8 + (_s-1)*(-5.f/16              \
                + (_s-1)*(35.f/128 + (_s-1)*(-63.f/256)))));                   \
        _dg = (2*_r)*(-1.f/2 + (_s-1)*(3.f/4 + (_s-1)*(-15.f/16 \
                + (_s-1)*(35.f/32 + (_s-1)*(-315.f/256)))));                   \
        break;                                                                 \
      case NL_MSM_SPLIT_TAYLOR6: \
        _g = 1 + (_s-1)*(-1.f/2 + (_s-1)*(3.f/8 + (_s-1)*(-5.f/16              \
                + (_s-1)*(35.f/128 + (_s-1)*(-63.f/256                         \
                    + (_s-1)*(231.f/1024))))));                                \
        _dg = (2*_r)*(-1.f/2 + (_s-1)*(3.f/4 + (_s-1)*(-15.f/16 \
                + (_s-1)*(35.f/32 + (_s-1)*(-315.f/256 \
                    + (_s-1)*(693.f/512))))));                  \
        break;                                                                 \
      case NL_MSM_SPLIT_TAYLOR7: \
        _g = 1 + (_s-1)*(-1.f/2 + (_s-1)*(3.f/8 + (_s-1)*(-5.f/16              \
                + (_s-1)*(35.f/128 + (_s-1)*(-63.f/256                         \
                    + (_s-1)*(231.f/1024 + (_s-1)*(-429.f/2048)))))));         \
        _dg = (2*_r)*(-1.f/2 + (_s-1)*(3.f/4 + (_s-1)*(-15.f/16 \
                + (_s-1)*(35.f/32 + (_s-1)*(-315.f/256 + (_s-1)*(693.f/512     \
                    + (_s-1)*(-3003.f/2048)))))));                             \
        break;                                                                 \
      case NL_MSM_SPLIT_TAYLOR8: \
        _g = 1 + (_s-1)*(-1.f/2 + (_s-1)*(3.f/8 + (_s-1)*(-5.f/16              \
                + (_s-1)*(35.f/128 + (_s-1)*(-63.f/256                         \
                    + (_s-1)*(231.f/1024 + (_s-1)*(-429.f/2048                 \
                        + (_s-1)*(6435.f/32768))))))));                        \
        _dg = (2*_r)*(-1.f/2 + (_s-1)*(3.f/4 + (_s-1)*(-15.f/16 \
                + (_s-1)*(35.f/32 + (_s-1)*(-315.f/256 + (_s-1)*(693.f/512     \
                    + (_s-1)*(-3003.f/2048 + (_s-1)*(6435.f/4096))))))));      \
        break;                                                                 \
      case NL_MSM_SPLIT_TAYLOR1: \
        _g = 1 + (_s-1)*(-1.f/2); \
        _dg = (2*_r)*(-1.f/2); \
        break; \
      /* from section 5.2 of thesis */ \
      case NL_MSM_SPLIT_SIGMA2_3:  /* the "perfect" smoothing */ \
        _g = 2 + _s*(-2 + _r);                                                \
        _dg = _r*(-4 + _r*3);                                                  \
        break;                                                                 \
      case NL_MSM_SPLIT_SIGMA3_5: \
        _g = 9.f/4 + _s*(-5.f/2 + _s*(9.f/4 - _r));                            \
        _dg = _r*(-5 + _s*(9 + _r*(-5)));                                     \
        break;                                                                 \
      case NL_MSM_SPLIT_SIGMA4_6: \
        _g = 21.f/8 + _s*(-35.f/8 + _s*(63.f/8 + _r*(-7 + _r*(15.f/8))));      \
        _dg = _r*(-35.f/4 + _s*(63.f/2 + _r*(-35 + _r*(45.f/4))));             \
        break;                                                                 \
      case NL_MSM_SPLIT_SIGMA4_7: \
        _g = 5.f/2 + _s*(-7.f/2 + _s*(7.f/2 + _s*(-5.f/2 + _r)));              \
        _dg = _r*(-7 + _s*(14 + _s*(-15 + _r*(7))));                         \
        break;                                                                 \
      case NL_MSM_SPLIT_SIGMA5_8: \
        _g = 45.f/16 + _s*(-21.f/4 + _s*(63.f/8 + _s*(-45.f/4                  \
                + _r*(9 + _r*(-35.f/16)))));                                   \
        _dg = _r*(-21.f/2 + _s*(63.f/2 + _s*(-135.f/2                          \
                + _r*(63 + _r*(-35.f/2)))));                                   \
        break;                                                                 \
      case NL_MSM_SPLIT_SIGMA5_9: \
        _g = 175.f/64 + _s*(-75.f/16 + _s*(189.f/32 + _s*(-75.f/16             \
                + _s*(175.f/64 - _r))));                                       \
        _dg = _r*(-75.f/8 + _s*(189.f/8 + _s*(-225.f/8 + _s*(175.f/8           \
                  + _r*(-9)))));                                               \
        break;                                                                 \
      case NL_MSM_SPLIT_SIGMA6_9: \
        _g = 25.f/8 + _s*(-15.f/2 + _s*(63.f/4 + _s*(-75.f/2                   \
                + _r*(45 + _r*(-175.f/8 + _r*4)))));                           \
        _dg = _r*(-15 + _s*(63 + _s*(-225                                    \
                + _r*(315 + _r*(-175 + _r*36)))));                             \
        break;                                                                 \
      case NL_MSM_SPLIT_SIGMA6_10: \
        _g = 385.f/128 + _s*(-825.f/128 + _s*(693.f/64 + _s*(-825.f/64         \
                + _s*(1925.f/128 + _r*(-11 + _r*(315.f/128))))));              \
        _dg = _r*(-825.f/64 + _s*(693.f/16 + _s*(-2475.f/32                    \
                + _s*(1925.f/16 + _r*(-99 + _r*(1575.f/64))))));               \
        break;                                                                 \
      case NL_MSM_SPLIT_SIGMA6_11: \
        _g = 189.f/64 + _s*(-385.f/64 + _s*(297.f/32 + _s*(-297.f/32           \
                + _s*(385.f/64 + _s*(-189.f/64 + _r)))));                      \
        _dg = _r*(-385.f/32 + _s*(297.f/8 + _s*(-891.f/16 + _s*(385.f/8        \
                  + _s*(-945.f/32 + _r*(11))))));                              \
        break;                                                                 \
      case NL_MSM_SPLIT_SIGMA7_11: \
        _g = 105.f/32 + _s*(-275.f/32 + _s*(297.f/16 + _s*(-495.f/16           \
                + _s*(1925.f/32 + _r*(-66 + _r*(945.f/32 + _r*(-5)))))));      \
        _dg = _r*(-275.f/16 + _s*(297.f/4 + _s*(-1485.f/8                      \
                + _s*(1925.f/4 + _r*(-594 + _r*(4725.f/16 + _r*(-55)))))));    \
        break;                                                                 \
      case NL_MSM_SPLIT_SIGMA7_12: \
        _g = 819.f/256 + _s*(-1001.f/128 + _s*(3861.f/256                      \
              + _s*(-1287.f/64 + _s*(5005.f/256 + _s*(-2457.f/128              \
                    + _r*(13 + _r*(-693.f/256)))))));                          \
        _dg = _r*(-1001.f/64 + _s*(3861.f/64 + _s*(-3861.f/32                  \
                + _s*(5005.f/32 + _s*(-12285.f/64 + _r*(143                    \
                      + _r*(-2079.f/64)))))));                                 \
        break;                                                                 \
      case NL_MSM_SPLIT_SIGMA7_13: \
        _g = 1617.f/512 + _s*(-1911.f/256 + _s*(7007.f/512 + _s*(-2145.f/128   \
                + _s*(7007.f/512 + _s*(-1911.f/256 + _s*(1617.f/512 - _r))))));\
        _dg = _r*(-1911.f/128 + _s*(7007.f/128 + _s*(-6435.f/64 + _s*(7007.f/64\
                  + _s*(-9555.f/128 + _s*(4851.f/128 + _r*(-13)))))));         \
        break;                                                                 \
      case NL_MSM_SPLIT_SIGMA8_12: \
        _g = 455.f/128 + _s*(-715.f/64 + _s*(3861.f/128 + _s*(-2145.f/32       \
                + _s*(25025.f/128 + _r*(-286 + _r*(12285.f/64 + _r*(-65        \
                        + _r*(1155.f/128))))))));                              \
        _dg = _r*(-715.f/32 + _s*(3861.f/32 + _s*(-6435.f/16                   \
                + _s*(25025.f/16 + _r*(-2574 + _r*(61425.f/32 + _r*(-715       \
                        + _r*(3465.f/32))))))));                               \
        break;                                                                 \
      case NL_MSM_SPLIT_SIGMA8_13: \
        _g = 441.f/128 + _s*(-637.f/64 + _s*(3003.f/128                        \
              + _s*(-1287.f/32 + _s*(7007.f/128 + _s*(-5733.f/64               \
                    + _r*(91 + _r*(-4851.f/128 + _r*(6))))))));                \
        _dg = _r*(-637.f/32 + _s*(3003.f/32 + _s*(-3861.f/16                   \
                + _s*(7007.f/16 + _s*(-28665.f/32 + _r*(1001                   \
                      + _r*(-14553.f/32 + _r*(78))))))));                      \
        break;                                                                 \
      case NL_MSM_SPLIT_SIGMA8_14: \
        _g = 3465.f/1024 + _s*(-9555.f/1024 + _s*(21021.f/1024                 \
              + _s*(-32175.f/1024 + _s*(35035.f/1024 + _s*(-28665.f/1024       \
                    + _s*(24255.f/1024 + _r*(-15 + _r*(3003.f/1024))))))));    \
        _dg = _r*(-9555.f/512 + _s*(21021.f/256 + _s*(-96525.f/512             \
                + _s*(35035.f/128 + _s*(-143325.f/512 + _s*(72765.f/256        \
                      + _r*(-195 + _r*(21021.f/512))))))));                    \
        break;                                                                 \
      case NL_MSM_SPLIT_SIGMA8_15: \
        _g = 429.f/128 + _s*(-1155.f/128 + _s*(2457.f/128 + _s*(-3575.f/128    \
                + _s*(3575.f/128 + _s*(-2457.f/128 + _s*(1155.f/128            \
                      + _s*(-429.f/128 + _r)))))));                            \
        _dg = _r*(-1155.f/64 + _s*(2457.f/32 + _s*(-10725.f/64                 \
                + _s*(3575.f/16 + _s*(-12285.f/64 + _s*(3465.f/32              \
                      + _s*(-3003.f/64 + _r*(15))))))));                       \
        break;                                                                 \
      /* from section 5.3 of thesis */ \
      case NL_MSM_SPLIT_SIGMA2_6: \
        _g = (31.f/16) + _s*(-23.f/16 + _s*(9.f/16 + _s*(-1.f/16))); \
        _dg = (2*_r)*(-23.f/16 + _s*(9.f/8 + _s*(-3.f/16))); \
        break; \
      /* from section 5.4 of thesis */ \
      case NL_MSM_SPLIT_SWITCH1_2:                                             \
        if (_r > 1.f/2) {                                                      \
          _g = 5.f/3 + _r + _s*(-3 + _r*(4.f/3));                              \
          _dg = 1 + _r*(-6 + _r*(4));                                          \
        }                                                                      \
        else {                                                                 \
          _g = 11.f/6 - _s;                                                    \
          _dg = _r*(-2);                                                       \
        }                                                                      \
        break;                                                                 \
      case NL_MSM_SPLIT_SWITCH3_4:                                             \
        if (_r > 3.f/4) {                                                      \
          _g = 5.f/7 + _r*(27.f/7 + _r*(-41.f/7 + _r*(16.f/7)));               \
          _dg = 27.f/7 + _r*(-82.f/7 + _r*(48.f/7));                           \
        }                                                                      \
        else {                                                                 \
          _g = 47.f/28 + _s*(-5.f/7);                                          \
          _dg = _r*(-10.f/7);                                                  \
        }                                                                      \
        break;                                                                 \
      case NL_MSM_SPLIT_SWITCH7_8:                                             \
        if (_r > 7.f/8) {                                                      \
          _g = -19.f/15 + _r*(49.f/5 + _r*(-59.f/5 + _r*(64.f/15)));           \
          _dg = 49.f/5 + _r*(-118.f/5 + _r*(64.f/5));                          \
        }                                                                      \
        else {                                                                 \
          _g = 191.f/120 + _s*(-3.f/5);                                        \
          _dg = _r*(-6.f/5);                                                   \
        }                                                                      \
        break;                                                                 \
      default: \
        return NL_MSM_ERROR_SUPPORT; \
    } \
    *(pg) = _g; \
    *(pdg) = _dg; \
  } while (0)
  /* closing ';' from use as function call */


#ifdef __cplusplus
extern "C" {
#endif


  /* creates type NL_Msmgrid_double */
  GRID_TEMPLATE(double);   /* for MSM charge and potential grids */

  /* creates type NL_Msmgrid_float */
  GRID_TEMPLATE(float);    /* single prec MSM charge and potential grids */


  /** Container and data structures for MSM solver. */
  struct NL_Msm_t {
    double uelec;          /* calculate electrostatic potential energy */

    double *felec;         /* calculate electrostatic forces x/y/z/... */
    const double *atom;    /* the original atom array stored x/y/z/q/... */

    float *felec_f;        /* single prec elec forces x/y/z/... */
    const float *atom_f;   /* single prec atom array stored x/y/z/q/... */

    double cellvec1[3];    /* cell basis vector (along x) */
    double cellvec2[3];    /* cell basis vector (along y) */
    double cellvec3[3];    /* cell basis vector (along z) */
    double cellcenter[3];  /* cell center */
    double recipvec1[3];   /* reciprocal space row vector */
    double recipvec2[3];   /* reciprocal space row vector */
    double recipvec3[3];   /* reciprocal space row vector */

    float cellvec1_f[3];   /* single prec cell basis vector (along x) */
    float cellvec2_f[3];   /* single prec cell basis vector (along y) */
    float cellvec3_f[3];   /* single prec cell basis vector (along z) */
    float cellcenter_f[3]; /* single prec cell center */
    float recipvec1_f[3];  /* single prec reciprocal space row vector */
    float recipvec2_f[3];  /* single prec reciprocal space row vector */
    float recipvec3_f[3];  /* single prec reciprocal space row vector */

    int msmflags;          /* bitwise flags from enum in msm.h */

    int numbins;           /* number of bins (nbx * nby * nbz) */
    int maxbins;           /* maximum number of bins allocated */
    int nbx;               /* number of bins along first cell basis vector */
    int nby;               /* number of bins along second cell basis vector */
    int nbz;               /* number of bins along third cell basis vector */
    int *bin;              /* array length maxbins, index of first atom */

    int numatoms;          /* number of atoms */
    int maxatoms;          /* maximum number of atoms allocated */
    int *next;             /* array length maxatoms, index next atom in bin */

    /* additional data members needed for k-away neighborhood */
    double binfill;        /* ratio of number of atoms per bin (default 0.8) */
    double density;        /* expected density of atoms (default 0.1) */
    int nbinslots;         /* how many atoms per bin (default 8) */

    double bu[3];          /* bin basis vector along u */
    double bv[3];          /* bin basis vector along v */
    double bw[3];          /* bin basis vector along w */

    int *nbrhd;
    int nbrhdlen;
    int nbrhdmax;
    int nradius;

    /*
     * Fundamental MSM parameters:
     *
     * Default grid spacing hmin = 2.5 A,
     * C1 cubic interpolation, and C2 Taylor splitting,
     * together give reasonable accuracy for atomistic systems
     * using a cutoff between 8 and 12 A.
     *
     * Find grid spacings in range hmin <= hx, hy, hz <= (1.5 * hmin)
     * such that along periodic dimensions the number of grid points is
     * some power of 2 times one or zero powers of 3.
     *
     * For reasonable accuracy, maintain ratio 4 <= a/hx, a/hy, a/hz <= 6
     * (and grid cutoff stencils fit into GPU constant memory cache).
     */

    double gridspacing;    /* smallest MSM grid spacing, hmax = 1.5 * hmin */
    double hx, hy, hz;     /* the finest level grid spacings */
    double a;              /* the MSM cutoff between short- and long-range */

    float hx_f, hy_f, hz_f;/* single prec finest level grid spacings */
    float a_f;             /* single prec MSM cutoff */

    int nx, ny, nz;        /* count number of h spacings that cover cell */

    int approx;            /* the approximation/interpolation: "MSM_APPROX_" */
    int split;             /* the splitting: "MSM_SPLIT_" */

    int nlevels;           /* number of grid levels */

    double gx, gy, gz;     /* coordinate of h-grid (0,0,0) point */
    float gx_f, gy_f, gz_f;/* single prec coordinate of h-grid (0,0,0) point */

    double gzero;          /* self energy factor from chosen splitting */
    float gzero_f;         /* single prec self energy factor */

    /*
     * Grids for calculating long-range part:
     * q[0] is finest-spaced grid (hx, hy, hz),
     * grid level k has spacing 2^k * (hx, hy, hz).
     *
     * Use domain origin (px0, py0, pz0) for each grid origin, ao that
     * q[0](i,j,k) has position (i*hx + px0, j*hy + py0, k*hz + pz0)
     * the finest level grid is 0, e.g. q0 = qh[0].
     *
     * Grid indexes can be negative for non-periodic dimensions,
     * the periodic dimensions wrap around edges of grid.
     */

    NL_Msmgrid_double *qh; /* grids of charge, 0==finest */
    NL_Msmgrid_double *eh; /* grids of potential, 0==finest */
    NL_Msmgrid_double *gc; /* fixed-size grids of weights for convolution */

    NL_Msmgrid_float *qh_f;/* single prec grids of charge, 0==finest */
    NL_Msmgrid_float *eh_f;/* single prec grids of potential, 0==finest */
    NL_Msmgrid_float *gc_f;/* single prec fixed-size grids of weights */

    int maxlevels;         /* maximum number of grid levels allocated */

    int max_lzd, max_lyzd; /* alloc length of factor temp buffer space */
    double *lzd;           /* factor temp row buffer, length z-dim of h-level */
    double *lyzd;          /* factor temp row buffer, (y*z)-dim of h-level */

    float *lzd_f;          /* single prec factor temp row buffer, z-dim */
    float *lyzd_f;         /* single prec factor temp row buffer, (y*z)-dim */

    int report_timings;    /* Do we report timings? */
    wkf_timerhandle timer; /* timer from John Stone */
    wkf_timerhandle timer_longrng; /* time individual long-range parts */

    /* CUDA atom cutoff kernel */

    /* CUDA grid cutoff kernel */
    int   lk_nlevels;      /* number of levels for grid cutoff kernel */
    int   lk_srad;         /* subcube radius for grid cutoff kernel */
    int   lk_padding;      /* padding around internal array of subcubes */
    int   subcube_total;   /* total number of subcubes for compressed grids */
    int   block_total;     /* total number of thread blocks */
    /*
     * host_   -->  memory allocated on host
     * device_ -->  global memory allocated on device
     */
    int   *host_sinfo;     /* subcube info copy to device const mem */
    float *host_lfac;      /* level factor copy to device const mem */
    int lk_maxlevels;

    float *host_wt;        /* weights copy to device const mem */
    int lk_maxwts;

    float *host_qgrids;    /* q-grid subcubes copy to device global mem */
    float *host_egrids;    /* e-grid subcubes copy to device global mem */
    float *device_qgrids;  /* q-grid subcubes allocate on device */
    float *device_egrids;  /* e-grid subcubes allocate on device */
    int lk_maxgridpts;   
  };


  /* prototypes for internal use */
  void NL_msm_cleanup(NL_Msm *pm);
  int NL_msm_compute_short_range(NL_Msm *pm);
  int NL_msm_compute_long_range(NL_Msm *pm);

  int NL_msm_compute_short_range_sprec(NL_Msm *pm);
  int NL_msm_compute_long_range_sprec(NL_Msm *pm);

#ifdef NL_USE_CUDA
  /* CUDA routines for grid cutoff computation */
  int NL_msm_cuda_setup_gridcutoff(NL_Msm *);
  void NL_msm_cuda_cleanup_gridcutoff(NL_Msm *);

  int NL_msm_cuda_compute_gridcutoff(NL_Msm *);
  int NL_msm_cuda_condense_qgrids(NL_Msm *);
  int NL_msm_cuda_expand_egrids(NL_Msm *);
#endif


#ifdef __cplusplus
}
#endif

#endif /* NL_MSM_DEFN_H */
