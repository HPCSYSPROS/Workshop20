#include <assert.h>
#include "vectors.h"

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder21(u) (kmul(p1o2dx,ksub(KRANC_GFOFFSET3D(u,1,0,0),KRANC_GFOFFSET3D(u,-1,0,0))))
#else
#  define PDstandardNthfdOrder21(u) (PDstandardNthfdOrder21_impl(u,p1o2dx,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder21_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o2dx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder21_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o2dx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o2dx,ksub(KRANC_GFOFFSET3D(u,1,0,0),KRANC_GFOFFSET3D(u,-1,0,0)));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder22(u) (kmul(p1o2dy,ksub(KRANC_GFOFFSET3D(u,0,1,0),KRANC_GFOFFSET3D(u,0,-1,0))))
#else
#  define PDstandardNthfdOrder22(u) (PDstandardNthfdOrder22_impl(u,p1o2dy,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder22_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o2dy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder22_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o2dy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o2dy,ksub(KRANC_GFOFFSET3D(u,0,1,0),KRANC_GFOFFSET3D(u,0,-1,0)));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder23(u) (kmul(p1o2dz,ksub(KRANC_GFOFFSET3D(u,0,0,1),KRANC_GFOFFSET3D(u,0,0,-1))))
#else
#  define PDstandardNthfdOrder23(u) (PDstandardNthfdOrder23_impl(u,p1o2dz,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder23_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o2dz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder23_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o2dz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardNthfdOrder22_impl(u, p1o2dz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder41(u) (kmul(p1o12dx,kmadd(ToReal(-8),KRANC_GFOFFSET3D(u,-1,0,0),kmadd(ToReal(8),KRANC_GFOFFSET3D(u,1,0,0),ksub(KRANC_GFOFFSET3D(u,-2,0,0),KRANC_GFOFFSET3D(u,2,0,0))))))
#else
#  define PDstandardNthfdOrder41(u) (PDstandardNthfdOrder41_impl(u,p1o12dx,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder41_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o12dx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder41_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o12dx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o12dx,kmadd(ToReal(-8),KRANC_GFOFFSET3D(u,-1,0,0),kmadd(ToReal(8),KRANC_GFOFFSET3D(u,1,0,0),ksub(KRANC_GFOFFSET3D(u,-2,0,0),KRANC_GFOFFSET3D(u,2,0,0)))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder42(u) (kmul(p1o12dy,kmadd(ToReal(-8),KRANC_GFOFFSET3D(u,0,-1,0),kmadd(ToReal(8),KRANC_GFOFFSET3D(u,0,1,0),ksub(KRANC_GFOFFSET3D(u,0,-2,0),KRANC_GFOFFSET3D(u,0,2,0))))))
#else
#  define PDstandardNthfdOrder42(u) (PDstandardNthfdOrder42_impl(u,p1o12dy,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder42_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o12dy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder42_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o12dy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o12dy,kmadd(ToReal(-8),KRANC_GFOFFSET3D(u,0,-1,0),kmadd(ToReal(8),KRANC_GFOFFSET3D(u,0,1,0),ksub(KRANC_GFOFFSET3D(u,0,-2,0),KRANC_GFOFFSET3D(u,0,2,0)))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder43(u) (kmul(p1o12dz,kmadd(ToReal(-8),KRANC_GFOFFSET3D(u,0,0,-1),kmadd(ToReal(8),KRANC_GFOFFSET3D(u,0,0,1),ksub(KRANC_GFOFFSET3D(u,0,0,-2),KRANC_GFOFFSET3D(u,0,0,2))))))
#else
#  define PDstandardNthfdOrder43(u) (PDstandardNthfdOrder43_impl(u,p1o12dz,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder43_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o12dz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder43_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o12dz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardNthfdOrder42_impl(u, p1o12dz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder61(u) (kmul(p1o60dx,kmadd(ToReal(-45),KRANC_GFOFFSET3D(u,-1,0,0),kmadd(ToReal(45),KRANC_GFOFFSET3D(u,1,0,0),kmadd(ToReal(9),KRANC_GFOFFSET3D(u,-2,0,0),kmadd(ToReal(-9),KRANC_GFOFFSET3D(u,2,0,0),ksub(KRANC_GFOFFSET3D(u,3,0,0),KRANC_GFOFFSET3D(u,-3,0,0))))))))
#else
#  define PDstandardNthfdOrder61(u) (PDstandardNthfdOrder61_impl(u,p1o60dx,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder61_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o60dx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder61_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o60dx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o60dx,kmadd(ToReal(-45),KRANC_GFOFFSET3D(u,-1,0,0),kmadd(ToReal(45),KRANC_GFOFFSET3D(u,1,0,0),kmadd(ToReal(9),KRANC_GFOFFSET3D(u,-2,0,0),kmadd(ToReal(-9),KRANC_GFOFFSET3D(u,2,0,0),ksub(KRANC_GFOFFSET3D(u,3,0,0),KRANC_GFOFFSET3D(u,-3,0,0)))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder62(u) (kmul(p1o60dy,kmadd(ToReal(-45),KRANC_GFOFFSET3D(u,0,-1,0),kmadd(ToReal(45),KRANC_GFOFFSET3D(u,0,1,0),kmadd(ToReal(9),KRANC_GFOFFSET3D(u,0,-2,0),kmadd(ToReal(-9),KRANC_GFOFFSET3D(u,0,2,0),ksub(KRANC_GFOFFSET3D(u,0,3,0),KRANC_GFOFFSET3D(u,0,-3,0))))))))
#else
#  define PDstandardNthfdOrder62(u) (PDstandardNthfdOrder62_impl(u,p1o60dy,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder62_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o60dy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder62_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o60dy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o60dy,kmadd(ToReal(-45),KRANC_GFOFFSET3D(u,0,-1,0),kmadd(ToReal(45),KRANC_GFOFFSET3D(u,0,1,0),kmadd(ToReal(9),KRANC_GFOFFSET3D(u,0,-2,0),kmadd(ToReal(-9),KRANC_GFOFFSET3D(u,0,2,0),ksub(KRANC_GFOFFSET3D(u,0,3,0),KRANC_GFOFFSET3D(u,0,-3,0)))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder63(u) (kmul(p1o60dz,kmadd(ToReal(-45),KRANC_GFOFFSET3D(u,0,0,-1),kmadd(ToReal(45),KRANC_GFOFFSET3D(u,0,0,1),kmadd(ToReal(9),KRANC_GFOFFSET3D(u,0,0,-2),kmadd(ToReal(-9),KRANC_GFOFFSET3D(u,0,0,2),ksub(KRANC_GFOFFSET3D(u,0,0,3),KRANC_GFOFFSET3D(u,0,0,-3))))))))
#else
#  define PDstandardNthfdOrder63(u) (PDstandardNthfdOrder63_impl(u,p1o60dz,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder63_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o60dz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder63_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o60dz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardNthfdOrder62_impl(u, p1o60dz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder81(u) (kmul(p1o840dx,kmadd(ToReal(-672),KRANC_GFOFFSET3D(u,-1,0,0),kmadd(ToReal(672),KRANC_GFOFFSET3D(u,1,0,0),kmadd(ToReal(168),KRANC_GFOFFSET3D(u,-2,0,0),kmadd(ToReal(-168),KRANC_GFOFFSET3D(u,2,0,0),kmadd(ToReal(-32),KRANC_GFOFFSET3D(u,-3,0,0),kmadd(ToReal(32),KRANC_GFOFFSET3D(u,3,0,0),kmadd(ToReal(3),KRANC_GFOFFSET3D(u,-4,0,0),kmul(ToReal(-3),KRANC_GFOFFSET3D(u,4,0,0)))))))))))
#else
#  define PDstandardNthfdOrder81(u) (PDstandardNthfdOrder81_impl(u,p1o840dx,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder81_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o840dx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder81_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o840dx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o840dx,kmadd(ToReal(-672),KRANC_GFOFFSET3D(u,-1,0,0),kmadd(ToReal(672),KRANC_GFOFFSET3D(u,1,0,0),kmadd(ToReal(168),KRANC_GFOFFSET3D(u,-2,0,0),kmadd(ToReal(-168),KRANC_GFOFFSET3D(u,2,0,0),kmadd(ToReal(-32),KRANC_GFOFFSET3D(u,-3,0,0),kmadd(ToReal(32),KRANC_GFOFFSET3D(u,3,0,0),kmadd(ToReal(3),KRANC_GFOFFSET3D(u,-4,0,0),kmul(ToReal(-3),KRANC_GFOFFSET3D(u,4,0,0))))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder82(u) (kmul(p1o840dy,kmadd(ToReal(-672),KRANC_GFOFFSET3D(u,0,-1,0),kmadd(ToReal(672),KRANC_GFOFFSET3D(u,0,1,0),kmadd(ToReal(168),KRANC_GFOFFSET3D(u,0,-2,0),kmadd(ToReal(-168),KRANC_GFOFFSET3D(u,0,2,0),kmadd(ToReal(-32),KRANC_GFOFFSET3D(u,0,-3,0),kmadd(ToReal(32),KRANC_GFOFFSET3D(u,0,3,0),kmadd(ToReal(3),KRANC_GFOFFSET3D(u,0,-4,0),kmul(ToReal(-3),KRANC_GFOFFSET3D(u,0,4,0)))))))))))
#else
#  define PDstandardNthfdOrder82(u) (PDstandardNthfdOrder82_impl(u,p1o840dy,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder82_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o840dy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder82_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o840dy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o840dy,kmadd(ToReal(-672),KRANC_GFOFFSET3D(u,0,-1,0),kmadd(ToReal(672),KRANC_GFOFFSET3D(u,0,1,0),kmadd(ToReal(168),KRANC_GFOFFSET3D(u,0,-2,0),kmadd(ToReal(-168),KRANC_GFOFFSET3D(u,0,2,0),kmadd(ToReal(-32),KRANC_GFOFFSET3D(u,0,-3,0),kmadd(ToReal(32),KRANC_GFOFFSET3D(u,0,3,0),kmadd(ToReal(3),KRANC_GFOFFSET3D(u,0,-4,0),kmul(ToReal(-3),KRANC_GFOFFSET3D(u,0,4,0))))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder83(u) (kmul(p1o840dz,kmadd(ToReal(-672),KRANC_GFOFFSET3D(u,0,0,-1),kmadd(ToReal(672),KRANC_GFOFFSET3D(u,0,0,1),kmadd(ToReal(168),KRANC_GFOFFSET3D(u,0,0,-2),kmadd(ToReal(-168),KRANC_GFOFFSET3D(u,0,0,2),kmadd(ToReal(-32),KRANC_GFOFFSET3D(u,0,0,-3),kmadd(ToReal(32),KRANC_GFOFFSET3D(u,0,0,3),kmadd(ToReal(3),KRANC_GFOFFSET3D(u,0,0,-4),kmul(ToReal(-3),KRANC_GFOFFSET3D(u,0,0,4)))))))))))
#else
#  define PDstandardNthfdOrder83(u) (PDstandardNthfdOrder83_impl(u,p1o840dz,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder83_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o840dz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder83_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o840dz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardNthfdOrder82_impl(u, p1o840dz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder211(u) (kmul(p1odx2,kmadd(ToReal(-2),KRANC_GFOFFSET3D(u,0,0,0),kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)))))
#else
#  define PDstandardNthfdOrder211(u) (PDstandardNthfdOrder211_impl(u,p1odx2,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder211_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1odx2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder211_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1odx2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1odx2,kmadd(ToReal(-2),KRANC_GFOFFSET3D(u,0,0,0),kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder222(u) (kmul(p1ody2,kmadd(ToReal(-2),KRANC_GFOFFSET3D(u,0,0,0),kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)))))
#else
#  define PDstandardNthfdOrder222(u) (PDstandardNthfdOrder222_impl(u,p1ody2,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder222_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1ody2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder222_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1ody2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1ody2,kmadd(ToReal(-2),KRANC_GFOFFSET3D(u,0,0,0),kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder233(u) (kmul(p1odz2,kmadd(ToReal(-2),KRANC_GFOFFSET3D(u,0,0,0),kadd(KRANC_GFOFFSET3D(u,0,0,-1),KRANC_GFOFFSET3D(u,0,0,1)))))
#else
#  define PDstandardNthfdOrder233(u) (PDstandardNthfdOrder233_impl(u,p1odz2,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder233_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1odz2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder233_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1odz2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardNthfdOrder222_impl(u, p1odz2, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder411(u) (kmul(pm1o12dx2,kmadd(ToReal(30),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(-16),kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),kadd(KRANC_GFOFFSET3D(u,-2,0,0),KRANC_GFOFFSET3D(u,2,0,0))))))
#else
#  define PDstandardNthfdOrder411(u) (PDstandardNthfdOrder411_impl(u,pm1o12dx2,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder411_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o12dx2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder411_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o12dx2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(pm1o12dx2,kmadd(ToReal(30),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(-16),kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),kadd(KRANC_GFOFFSET3D(u,-2,0,0),KRANC_GFOFFSET3D(u,2,0,0)))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder422(u) (kmul(pm1o12dy2,kmadd(ToReal(30),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(-16),kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),kadd(KRANC_GFOFFSET3D(u,0,-2,0),KRANC_GFOFFSET3D(u,0,2,0))))))
#else
#  define PDstandardNthfdOrder422(u) (PDstandardNthfdOrder422_impl(u,pm1o12dy2,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder422_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o12dy2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder422_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o12dy2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(pm1o12dy2,kmadd(ToReal(30),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(-16),kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),kadd(KRANC_GFOFFSET3D(u,0,-2,0),KRANC_GFOFFSET3D(u,0,2,0)))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder433(u) (kmul(pm1o12dz2,kmadd(ToReal(30),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(-16),kadd(KRANC_GFOFFSET3D(u,0,0,-1),KRANC_GFOFFSET3D(u,0,0,1)),kadd(KRANC_GFOFFSET3D(u,0,0,-2),KRANC_GFOFFSET3D(u,0,0,2))))))
#else
#  define PDstandardNthfdOrder433(u) (PDstandardNthfdOrder433_impl(u,pm1o12dz2,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder433_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o12dz2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder433_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o12dz2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardNthfdOrder422_impl(u, pm1o12dz2, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder611(u) (kmul(p1o180dx2,kmadd(ToReal(-490),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(270),kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),kmadd(ToReal(-27),kadd(KRANC_GFOFFSET3D(u,-2,0,0),KRANC_GFOFFSET3D(u,2,0,0)),kmul(ToReal(2),kadd(KRANC_GFOFFSET3D(u,-3,0,0),KRANC_GFOFFSET3D(u,3,0,0))))))))
#else
#  define PDstandardNthfdOrder611(u) (PDstandardNthfdOrder611_impl(u,p1o180dx2,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder611_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o180dx2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder611_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o180dx2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o180dx2,kmadd(ToReal(-490),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(270),kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),kmadd(ToReal(-27),kadd(KRANC_GFOFFSET3D(u,-2,0,0),KRANC_GFOFFSET3D(u,2,0,0)),kmul(ToReal(2),kadd(KRANC_GFOFFSET3D(u,-3,0,0),KRANC_GFOFFSET3D(u,3,0,0)))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder622(u) (kmul(p1o180dy2,kmadd(ToReal(-490),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(270),kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),kmadd(ToReal(-27),kadd(KRANC_GFOFFSET3D(u,0,-2,0),KRANC_GFOFFSET3D(u,0,2,0)),kmul(ToReal(2),kadd(KRANC_GFOFFSET3D(u,0,-3,0),KRANC_GFOFFSET3D(u,0,3,0))))))))
#else
#  define PDstandardNthfdOrder622(u) (PDstandardNthfdOrder622_impl(u,p1o180dy2,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder622_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o180dy2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder622_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o180dy2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o180dy2,kmadd(ToReal(-490),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(270),kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),kmadd(ToReal(-27),kadd(KRANC_GFOFFSET3D(u,0,-2,0),KRANC_GFOFFSET3D(u,0,2,0)),kmul(ToReal(2),kadd(KRANC_GFOFFSET3D(u,0,-3,0),KRANC_GFOFFSET3D(u,0,3,0)))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder633(u) (kmul(p1o180dz2,kmadd(ToReal(-490),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(270),kadd(KRANC_GFOFFSET3D(u,0,0,-1),KRANC_GFOFFSET3D(u,0,0,1)),kmadd(ToReal(-27),kadd(KRANC_GFOFFSET3D(u,0,0,-2),KRANC_GFOFFSET3D(u,0,0,2)),kmul(ToReal(2),kadd(KRANC_GFOFFSET3D(u,0,0,-3),KRANC_GFOFFSET3D(u,0,0,3))))))))
#else
#  define PDstandardNthfdOrder633(u) (PDstandardNthfdOrder633_impl(u,p1o180dz2,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder633_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o180dz2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder633_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o180dz2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardNthfdOrder622_impl(u, p1o180dz2, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder811(u) (kmul(p1o5040dx2,kmadd(ToReal(-14350),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(8064),kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),kmadd(ToReal(-1008),kadd(KRANC_GFOFFSET3D(u,-2,0,0),KRANC_GFOFFSET3D(u,2,0,0)),kmadd(ToReal(128),kadd(KRANC_GFOFFSET3D(u,-3,0,0),KRANC_GFOFFSET3D(u,3,0,0)),kmul(ToReal(-9),kadd(KRANC_GFOFFSET3D(u,-4,0,0),KRANC_GFOFFSET3D(u,4,0,0)))))))))
#else
#  define PDstandardNthfdOrder811(u) (PDstandardNthfdOrder811_impl(u,p1o5040dx2,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder811_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o5040dx2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder811_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o5040dx2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o5040dx2,kmadd(ToReal(-14350),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(8064),kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),kmadd(ToReal(-1008),kadd(KRANC_GFOFFSET3D(u,-2,0,0),KRANC_GFOFFSET3D(u,2,0,0)),kmadd(ToReal(128),kadd(KRANC_GFOFFSET3D(u,-3,0,0),KRANC_GFOFFSET3D(u,3,0,0)),kmul(ToReal(-9),kadd(KRANC_GFOFFSET3D(u,-4,0,0),KRANC_GFOFFSET3D(u,4,0,0))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder822(u) (kmul(p1o5040dy2,kmadd(ToReal(-14350),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(8064),kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),kmadd(ToReal(-1008),kadd(KRANC_GFOFFSET3D(u,0,-2,0),KRANC_GFOFFSET3D(u,0,2,0)),kmadd(ToReal(128),kadd(KRANC_GFOFFSET3D(u,0,-3,0),KRANC_GFOFFSET3D(u,0,3,0)),kmul(ToReal(-9),kadd(KRANC_GFOFFSET3D(u,0,-4,0),KRANC_GFOFFSET3D(u,0,4,0)))))))))
#else
#  define PDstandardNthfdOrder822(u) (PDstandardNthfdOrder822_impl(u,p1o5040dy2,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder822_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o5040dy2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder822_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o5040dy2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o5040dy2,kmadd(ToReal(-14350),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(8064),kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),kmadd(ToReal(-1008),kadd(KRANC_GFOFFSET3D(u,0,-2,0),KRANC_GFOFFSET3D(u,0,2,0)),kmadd(ToReal(128),kadd(KRANC_GFOFFSET3D(u,0,-3,0),KRANC_GFOFFSET3D(u,0,3,0)),kmul(ToReal(-9),kadd(KRANC_GFOFFSET3D(u,0,-4,0),KRANC_GFOFFSET3D(u,0,4,0))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder833(u) (kmul(p1o5040dz2,kmadd(ToReal(-14350),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(8064),kadd(KRANC_GFOFFSET3D(u,0,0,-1),KRANC_GFOFFSET3D(u,0,0,1)),kmadd(ToReal(-1008),kadd(KRANC_GFOFFSET3D(u,0,0,-2),KRANC_GFOFFSET3D(u,0,0,2)),kmadd(ToReal(128),kadd(KRANC_GFOFFSET3D(u,0,0,-3),KRANC_GFOFFSET3D(u,0,0,3)),kmul(ToReal(-9),kadd(KRANC_GFOFFSET3D(u,0,0,-4),KRANC_GFOFFSET3D(u,0,0,4)))))))))
#else
#  define PDstandardNthfdOrder833(u) (PDstandardNthfdOrder833_impl(u,p1o5040dz2,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder833_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o5040dz2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder833_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o5040dz2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardNthfdOrder822_impl(u, p1o5040dz2, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder212(u) (kmul(p1o4dxdy,kadd(KRANC_GFOFFSET3D(u,-1,-1,0),ksub(ksub(KRANC_GFOFFSET3D(u,1,1,0),KRANC_GFOFFSET3D(u,1,-1,0)),KRANC_GFOFFSET3D(u,-1,1,0)))))
#else
#  define PDstandardNthfdOrder212(u) (PDstandardNthfdOrder212_impl(u,p1o4dxdy,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder212_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o4dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder212_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o4dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o4dxdy,kadd(KRANC_GFOFFSET3D(u,-1,-1,0),ksub(ksub(KRANC_GFOFFSET3D(u,1,1,0),KRANC_GFOFFSET3D(u,1,-1,0)),KRANC_GFOFFSET3D(u,-1,1,0))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder213(u) (kmul(p1o4dxdz,kadd(KRANC_GFOFFSET3D(u,-1,0,-1),ksub(ksub(KRANC_GFOFFSET3D(u,1,0,1),KRANC_GFOFFSET3D(u,1,0,-1)),KRANC_GFOFFSET3D(u,-1,0,1)))))
#else
#  define PDstandardNthfdOrder213(u) (PDstandardNthfdOrder213_impl(u,p1o4dxdz,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder213_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o4dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder213_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o4dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardNthfdOrder212_impl(u, p1o4dxdz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder221(u) (kmul(p1o4dxdy,kadd(KRANC_GFOFFSET3D(u,-1,-1,0),ksub(ksub(KRANC_GFOFFSET3D(u,1,1,0),KRANC_GFOFFSET3D(u,1,-1,0)),KRANC_GFOFFSET3D(u,-1,1,0)))))
#else
#  define PDstandardNthfdOrder221(u) (PDstandardNthfdOrder221_impl(u,p1o4dxdy,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder221_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o4dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder221_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o4dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardNthfdOrder212_impl(u, p1o4dxdy, cdj, cdk);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder223(u) (kmul(p1o4dydz,kadd(KRANC_GFOFFSET3D(u,0,-1,-1),ksub(ksub(KRANC_GFOFFSET3D(u,0,1,1),KRANC_GFOFFSET3D(u,0,1,-1)),KRANC_GFOFFSET3D(u,0,-1,1)))))
#else
#  define PDstandardNthfdOrder223(u) (PDstandardNthfdOrder223_impl(u,p1o4dydz,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder223_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o4dydz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder223_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o4dydz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o4dydz,kadd(KRANC_GFOFFSET3D(u,0,-1,-1),ksub(ksub(KRANC_GFOFFSET3D(u,0,1,1),KRANC_GFOFFSET3D(u,0,1,-1)),KRANC_GFOFFSET3D(u,0,-1,1))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder231(u) (kmul(p1o4dxdz,kadd(KRANC_GFOFFSET3D(u,-1,0,-1),ksub(ksub(KRANC_GFOFFSET3D(u,1,0,1),KRANC_GFOFFSET3D(u,1,0,-1)),KRANC_GFOFFSET3D(u,-1,0,1)))))
#else
#  define PDstandardNthfdOrder231(u) (PDstandardNthfdOrder231_impl(u,p1o4dxdz,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder231_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o4dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder231_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o4dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardNthfdOrder212_impl(u, p1o4dxdz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder232(u) (kmul(p1o4dydz,kadd(KRANC_GFOFFSET3D(u,0,-1,-1),ksub(ksub(KRANC_GFOFFSET3D(u,0,1,1),KRANC_GFOFFSET3D(u,0,1,-1)),KRANC_GFOFFSET3D(u,0,-1,1)))))
#else
#  define PDstandardNthfdOrder232(u) (PDstandardNthfdOrder232_impl(u,p1o4dydz,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder232_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o4dydz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder232_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o4dydz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardNthfdOrder223_impl(u, p1o4dydz, cdj, cdk);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder412(u) (kmul(p1o144dxdy,kmadd(ToReal(-64),kadd(KRANC_GFOFFSET3D(u,-1,1,0),KRANC_GFOFFSET3D(u,1,-1,0)),kmadd(ToReal(64),kadd(KRANC_GFOFFSET3D(u,-1,-1,0),KRANC_GFOFFSET3D(u,1,1,0)),kmadd(ToReal(8),kadd(KRANC_GFOFFSET3D(u,-1,2,0),kadd(KRANC_GFOFFSET3D(u,1,-2,0),kadd(KRANC_GFOFFSET3D(u,-2,1,0),KRANC_GFOFFSET3D(u,2,-1,0)))),kmadd(ToReal(-8),kadd(KRANC_GFOFFSET3D(u,-1,-2,0),kadd(KRANC_GFOFFSET3D(u,1,2,0),kadd(KRANC_GFOFFSET3D(u,-2,-1,0),KRANC_GFOFFSET3D(u,2,1,0)))),kadd(KRANC_GFOFFSET3D(u,-2,-2,0),ksub(ksub(KRANC_GFOFFSET3D(u,2,2,0),KRANC_GFOFFSET3D(u,2,-2,0)),KRANC_GFOFFSET3D(u,-2,2,0)))))))))
#else
#  define PDstandardNthfdOrder412(u) (PDstandardNthfdOrder412_impl(u,p1o144dxdy,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder412_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o144dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder412_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o144dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o144dxdy,kmadd(ToReal(-64),kadd(KRANC_GFOFFSET3D(u,-1,1,0),KRANC_GFOFFSET3D(u,1,-1,0)),kmadd(ToReal(64),kadd(KRANC_GFOFFSET3D(u,-1,-1,0),KRANC_GFOFFSET3D(u,1,1,0)),kmadd(ToReal(8),kadd(KRANC_GFOFFSET3D(u,-1,2,0),kadd(KRANC_GFOFFSET3D(u,1,-2,0),kadd(KRANC_GFOFFSET3D(u,-2,1,0),KRANC_GFOFFSET3D(u,2,-1,0)))),kmadd(ToReal(-8),kadd(KRANC_GFOFFSET3D(u,-1,-2,0),kadd(KRANC_GFOFFSET3D(u,1,2,0),kadd(KRANC_GFOFFSET3D(u,-2,-1,0),KRANC_GFOFFSET3D(u,2,1,0)))),kadd(KRANC_GFOFFSET3D(u,-2,-2,0),ksub(ksub(KRANC_GFOFFSET3D(u,2,2,0),KRANC_GFOFFSET3D(u,2,-2,0)),KRANC_GFOFFSET3D(u,-2,2,0))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder413(u) (kmul(p1o144dxdz,kmadd(ToReal(-64),kadd(KRANC_GFOFFSET3D(u,-1,0,1),KRANC_GFOFFSET3D(u,1,0,-1)),kmadd(ToReal(64),kadd(KRANC_GFOFFSET3D(u,-1,0,-1),KRANC_GFOFFSET3D(u,1,0,1)),kmadd(ToReal(8),kadd(KRANC_GFOFFSET3D(u,-1,0,2),kadd(KRANC_GFOFFSET3D(u,1,0,-2),kadd(KRANC_GFOFFSET3D(u,-2,0,1),KRANC_GFOFFSET3D(u,2,0,-1)))),kmadd(ToReal(-8),kadd(KRANC_GFOFFSET3D(u,-1,0,-2),kadd(KRANC_GFOFFSET3D(u,1,0,2),kadd(KRANC_GFOFFSET3D(u,-2,0,-1),KRANC_GFOFFSET3D(u,2,0,1)))),kadd(KRANC_GFOFFSET3D(u,-2,0,-2),ksub(ksub(KRANC_GFOFFSET3D(u,2,0,2),KRANC_GFOFFSET3D(u,2,0,-2)),KRANC_GFOFFSET3D(u,-2,0,2)))))))))
#else
#  define PDstandardNthfdOrder413(u) (PDstandardNthfdOrder413_impl(u,p1o144dxdz,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder413_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o144dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder413_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o144dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardNthfdOrder412_impl(u, p1o144dxdz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder421(u) (kmul(p1o144dxdy,kmadd(ToReal(-64),kadd(KRANC_GFOFFSET3D(u,-1,1,0),KRANC_GFOFFSET3D(u,1,-1,0)),kmadd(ToReal(64),kadd(KRANC_GFOFFSET3D(u,-1,-1,0),KRANC_GFOFFSET3D(u,1,1,0)),kmadd(ToReal(8),kadd(KRANC_GFOFFSET3D(u,-1,2,0),kadd(KRANC_GFOFFSET3D(u,1,-2,0),kadd(KRANC_GFOFFSET3D(u,-2,1,0),KRANC_GFOFFSET3D(u,2,-1,0)))),kmadd(ToReal(-8),kadd(KRANC_GFOFFSET3D(u,-1,-2,0),kadd(KRANC_GFOFFSET3D(u,1,2,0),kadd(KRANC_GFOFFSET3D(u,-2,-1,0),KRANC_GFOFFSET3D(u,2,1,0)))),kadd(KRANC_GFOFFSET3D(u,-2,-2,0),ksub(ksub(KRANC_GFOFFSET3D(u,2,2,0),KRANC_GFOFFSET3D(u,2,-2,0)),KRANC_GFOFFSET3D(u,-2,2,0)))))))))
#else
#  define PDstandardNthfdOrder421(u) (PDstandardNthfdOrder421_impl(u,p1o144dxdy,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder421_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o144dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder421_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o144dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardNthfdOrder412_impl(u, p1o144dxdy, cdj, cdk);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder423(u) (kmul(p1o144dydz,kmadd(ToReal(-64),kadd(KRANC_GFOFFSET3D(u,0,-1,1),KRANC_GFOFFSET3D(u,0,1,-1)),kmadd(ToReal(64),kadd(KRANC_GFOFFSET3D(u,0,-1,-1),KRANC_GFOFFSET3D(u,0,1,1)),kmadd(ToReal(8),kadd(KRANC_GFOFFSET3D(u,0,-1,2),kadd(KRANC_GFOFFSET3D(u,0,1,-2),kadd(KRANC_GFOFFSET3D(u,0,-2,1),KRANC_GFOFFSET3D(u,0,2,-1)))),kmadd(ToReal(-8),kadd(KRANC_GFOFFSET3D(u,0,-1,-2),kadd(KRANC_GFOFFSET3D(u,0,1,2),kadd(KRANC_GFOFFSET3D(u,0,-2,-1),KRANC_GFOFFSET3D(u,0,2,1)))),kadd(KRANC_GFOFFSET3D(u,0,-2,-2),ksub(ksub(KRANC_GFOFFSET3D(u,0,2,2),KRANC_GFOFFSET3D(u,0,2,-2)),KRANC_GFOFFSET3D(u,0,-2,2)))))))))
#else
#  define PDstandardNthfdOrder423(u) (PDstandardNthfdOrder423_impl(u,p1o144dydz,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder423_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o144dydz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder423_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o144dydz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o144dydz,kmadd(ToReal(-64),kadd(KRANC_GFOFFSET3D(u,0,-1,1),KRANC_GFOFFSET3D(u,0,1,-1)),kmadd(ToReal(64),kadd(KRANC_GFOFFSET3D(u,0,-1,-1),KRANC_GFOFFSET3D(u,0,1,1)),kmadd(ToReal(8),kadd(KRANC_GFOFFSET3D(u,0,-1,2),kadd(KRANC_GFOFFSET3D(u,0,1,-2),kadd(KRANC_GFOFFSET3D(u,0,-2,1),KRANC_GFOFFSET3D(u,0,2,-1)))),kmadd(ToReal(-8),kadd(KRANC_GFOFFSET3D(u,0,-1,-2),kadd(KRANC_GFOFFSET3D(u,0,1,2),kadd(KRANC_GFOFFSET3D(u,0,-2,-1),KRANC_GFOFFSET3D(u,0,2,1)))),kadd(KRANC_GFOFFSET3D(u,0,-2,-2),ksub(ksub(KRANC_GFOFFSET3D(u,0,2,2),KRANC_GFOFFSET3D(u,0,2,-2)),KRANC_GFOFFSET3D(u,0,-2,2))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder431(u) (kmul(p1o144dxdz,kmadd(ToReal(-64),kadd(KRANC_GFOFFSET3D(u,-1,0,1),KRANC_GFOFFSET3D(u,1,0,-1)),kmadd(ToReal(64),kadd(KRANC_GFOFFSET3D(u,-1,0,-1),KRANC_GFOFFSET3D(u,1,0,1)),kmadd(ToReal(8),kadd(KRANC_GFOFFSET3D(u,-1,0,2),kadd(KRANC_GFOFFSET3D(u,1,0,-2),kadd(KRANC_GFOFFSET3D(u,-2,0,1),KRANC_GFOFFSET3D(u,2,0,-1)))),kmadd(ToReal(-8),kadd(KRANC_GFOFFSET3D(u,-1,0,-2),kadd(KRANC_GFOFFSET3D(u,1,0,2),kadd(KRANC_GFOFFSET3D(u,-2,0,-1),KRANC_GFOFFSET3D(u,2,0,1)))),kadd(KRANC_GFOFFSET3D(u,-2,0,-2),ksub(ksub(KRANC_GFOFFSET3D(u,2,0,2),KRANC_GFOFFSET3D(u,2,0,-2)),KRANC_GFOFFSET3D(u,-2,0,2)))))))))
#else
#  define PDstandardNthfdOrder431(u) (PDstandardNthfdOrder431_impl(u,p1o144dxdz,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder431_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o144dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder431_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o144dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardNthfdOrder412_impl(u, p1o144dxdz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder432(u) (kmul(p1o144dydz,kmadd(ToReal(-64),kadd(KRANC_GFOFFSET3D(u,0,-1,1),KRANC_GFOFFSET3D(u,0,1,-1)),kmadd(ToReal(64),kadd(KRANC_GFOFFSET3D(u,0,-1,-1),KRANC_GFOFFSET3D(u,0,1,1)),kmadd(ToReal(8),kadd(KRANC_GFOFFSET3D(u,0,-1,2),kadd(KRANC_GFOFFSET3D(u,0,1,-2),kadd(KRANC_GFOFFSET3D(u,0,-2,1),KRANC_GFOFFSET3D(u,0,2,-1)))),kmadd(ToReal(-8),kadd(KRANC_GFOFFSET3D(u,0,-1,-2),kadd(KRANC_GFOFFSET3D(u,0,1,2),kadd(KRANC_GFOFFSET3D(u,0,-2,-1),KRANC_GFOFFSET3D(u,0,2,1)))),kadd(KRANC_GFOFFSET3D(u,0,-2,-2),ksub(ksub(KRANC_GFOFFSET3D(u,0,2,2),KRANC_GFOFFSET3D(u,0,2,-2)),KRANC_GFOFFSET3D(u,0,-2,2)))))))))
#else
#  define PDstandardNthfdOrder432(u) (PDstandardNthfdOrder432_impl(u,p1o144dydz,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder432_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o144dydz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder432_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o144dydz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardNthfdOrder423_impl(u, p1o144dydz, cdj, cdk);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder612(u) (kmul(p1o3600dxdy,kmadd(ToReal(-2025),kadd(KRANC_GFOFFSET3D(u,-1,1,0),KRANC_GFOFFSET3D(u,1,-1,0)),kmadd(ToReal(2025),kadd(KRANC_GFOFFSET3D(u,-1,-1,0),KRANC_GFOFFSET3D(u,1,1,0)),kmadd(ToReal(405),kadd(KRANC_GFOFFSET3D(u,-1,2,0),kadd(KRANC_GFOFFSET3D(u,1,-2,0),kadd(KRANC_GFOFFSET3D(u,-2,1,0),KRANC_GFOFFSET3D(u,2,-1,0)))),kmadd(ToReal(-405),kadd(KRANC_GFOFFSET3D(u,-1,-2,0),kadd(KRANC_GFOFFSET3D(u,1,2,0),kadd(KRANC_GFOFFSET3D(u,-2,-1,0),KRANC_GFOFFSET3D(u,2,1,0)))),kmadd(ToReal(-81),kadd(KRANC_GFOFFSET3D(u,-2,2,0),KRANC_GFOFFSET3D(u,2,-2,0)),kmadd(ToReal(81),kadd(KRANC_GFOFFSET3D(u,-2,-2,0),KRANC_GFOFFSET3D(u,2,2,0)),kmadd(ToReal(-45),kadd(KRANC_GFOFFSET3D(u,-1,3,0),kadd(KRANC_GFOFFSET3D(u,1,-3,0),kadd(KRANC_GFOFFSET3D(u,-3,1,0),KRANC_GFOFFSET3D(u,3,-1,0)))),kmadd(ToReal(45),kadd(KRANC_GFOFFSET3D(u,-1,-3,0),kadd(KRANC_GFOFFSET3D(u,1,3,0),kadd(KRANC_GFOFFSET3D(u,-3,-1,0),KRANC_GFOFFSET3D(u,3,1,0)))),kmadd(ToReal(9),kadd(KRANC_GFOFFSET3D(u,-2,3,0),kadd(KRANC_GFOFFSET3D(u,2,-3,0),kadd(KRANC_GFOFFSET3D(u,-3,2,0),KRANC_GFOFFSET3D(u,3,-2,0)))),kmadd(ToReal(-9),kadd(KRANC_GFOFFSET3D(u,-2,-3,0),kadd(KRANC_GFOFFSET3D(u,2,3,0),kadd(KRANC_GFOFFSET3D(u,-3,-2,0),KRANC_GFOFFSET3D(u,3,2,0)))),kadd(KRANC_GFOFFSET3D(u,-3,-3,0),ksub(ksub(KRANC_GFOFFSET3D(u,3,3,0),KRANC_GFOFFSET3D(u,3,-3,0)),KRANC_GFOFFSET3D(u,-3,3,0)))))))))))))))
#else
#  define PDstandardNthfdOrder612(u) (PDstandardNthfdOrder612_impl(u,p1o3600dxdy,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder612_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o3600dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder612_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o3600dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o3600dxdy,kmadd(ToReal(-2025),kadd(KRANC_GFOFFSET3D(u,-1,1,0),KRANC_GFOFFSET3D(u,1,-1,0)),kmadd(ToReal(2025),kadd(KRANC_GFOFFSET3D(u,-1,-1,0),KRANC_GFOFFSET3D(u,1,1,0)),kmadd(ToReal(405),kadd(KRANC_GFOFFSET3D(u,-1,2,0),kadd(KRANC_GFOFFSET3D(u,1,-2,0),kadd(KRANC_GFOFFSET3D(u,-2,1,0),KRANC_GFOFFSET3D(u,2,-1,0)))),kmadd(ToReal(-405),kadd(KRANC_GFOFFSET3D(u,-1,-2,0),kadd(KRANC_GFOFFSET3D(u,1,2,0),kadd(KRANC_GFOFFSET3D(u,-2,-1,0),KRANC_GFOFFSET3D(u,2,1,0)))),kmadd(ToReal(-81),kadd(KRANC_GFOFFSET3D(u,-2,2,0),KRANC_GFOFFSET3D(u,2,-2,0)),kmadd(ToReal(81),kadd(KRANC_GFOFFSET3D(u,-2,-2,0),KRANC_GFOFFSET3D(u,2,2,0)),kmadd(ToReal(-45),kadd(KRANC_GFOFFSET3D(u,-1,3,0),kadd(KRANC_GFOFFSET3D(u,1,-3,0),kadd(KRANC_GFOFFSET3D(u,-3,1,0),KRANC_GFOFFSET3D(u,3,-1,0)))),kmadd(ToReal(45),kadd(KRANC_GFOFFSET3D(u,-1,-3,0),kadd(KRANC_GFOFFSET3D(u,1,3,0),kadd(KRANC_GFOFFSET3D(u,-3,-1,0),KRANC_GFOFFSET3D(u,3,1,0)))),kmadd(ToReal(9),kadd(KRANC_GFOFFSET3D(u,-2,3,0),kadd(KRANC_GFOFFSET3D(u,2,-3,0),kadd(KRANC_GFOFFSET3D(u,-3,2,0),KRANC_GFOFFSET3D(u,3,-2,0)))),kmadd(ToReal(-9),kadd(KRANC_GFOFFSET3D(u,-2,-3,0),kadd(KRANC_GFOFFSET3D(u,2,3,0),kadd(KRANC_GFOFFSET3D(u,-3,-2,0),KRANC_GFOFFSET3D(u,3,2,0)))),kadd(KRANC_GFOFFSET3D(u,-3,-3,0),ksub(ksub(KRANC_GFOFFSET3D(u,3,3,0),KRANC_GFOFFSET3D(u,3,-3,0)),KRANC_GFOFFSET3D(u,-3,3,0))))))))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder613(u) (kmul(p1o3600dxdz,kmadd(ToReal(-2025),kadd(KRANC_GFOFFSET3D(u,-1,0,1),KRANC_GFOFFSET3D(u,1,0,-1)),kmadd(ToReal(2025),kadd(KRANC_GFOFFSET3D(u,-1,0,-1),KRANC_GFOFFSET3D(u,1,0,1)),kmadd(ToReal(405),kadd(KRANC_GFOFFSET3D(u,-1,0,2),kadd(KRANC_GFOFFSET3D(u,1,0,-2),kadd(KRANC_GFOFFSET3D(u,-2,0,1),KRANC_GFOFFSET3D(u,2,0,-1)))),kmadd(ToReal(-405),kadd(KRANC_GFOFFSET3D(u,-1,0,-2),kadd(KRANC_GFOFFSET3D(u,1,0,2),kadd(KRANC_GFOFFSET3D(u,-2,0,-1),KRANC_GFOFFSET3D(u,2,0,1)))),kmadd(ToReal(-81),kadd(KRANC_GFOFFSET3D(u,-2,0,2),KRANC_GFOFFSET3D(u,2,0,-2)),kmadd(ToReal(81),kadd(KRANC_GFOFFSET3D(u,-2,0,-2),KRANC_GFOFFSET3D(u,2,0,2)),kmadd(ToReal(-45),kadd(KRANC_GFOFFSET3D(u,-1,0,3),kadd(KRANC_GFOFFSET3D(u,1,0,-3),kadd(KRANC_GFOFFSET3D(u,-3,0,1),KRANC_GFOFFSET3D(u,3,0,-1)))),kmadd(ToReal(45),kadd(KRANC_GFOFFSET3D(u,-1,0,-3),kadd(KRANC_GFOFFSET3D(u,1,0,3),kadd(KRANC_GFOFFSET3D(u,-3,0,-1),KRANC_GFOFFSET3D(u,3,0,1)))),kmadd(ToReal(9),kadd(KRANC_GFOFFSET3D(u,-2,0,3),kadd(KRANC_GFOFFSET3D(u,2,0,-3),kadd(KRANC_GFOFFSET3D(u,-3,0,2),KRANC_GFOFFSET3D(u,3,0,-2)))),kmadd(ToReal(-9),kadd(KRANC_GFOFFSET3D(u,-2,0,-3),kadd(KRANC_GFOFFSET3D(u,2,0,3),kadd(KRANC_GFOFFSET3D(u,-3,0,-2),KRANC_GFOFFSET3D(u,3,0,2)))),kadd(KRANC_GFOFFSET3D(u,-3,0,-3),ksub(ksub(KRANC_GFOFFSET3D(u,3,0,3),KRANC_GFOFFSET3D(u,3,0,-3)),KRANC_GFOFFSET3D(u,-3,0,3)))))))))))))))
#else
#  define PDstandardNthfdOrder613(u) (PDstandardNthfdOrder613_impl(u,p1o3600dxdz,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder613_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o3600dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder613_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o3600dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardNthfdOrder612_impl(u, p1o3600dxdz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder621(u) (kmul(p1o3600dxdy,kmadd(ToReal(-2025),kadd(KRANC_GFOFFSET3D(u,-1,1,0),KRANC_GFOFFSET3D(u,1,-1,0)),kmadd(ToReal(2025),kadd(KRANC_GFOFFSET3D(u,-1,-1,0),KRANC_GFOFFSET3D(u,1,1,0)),kmadd(ToReal(405),kadd(KRANC_GFOFFSET3D(u,-1,2,0),kadd(KRANC_GFOFFSET3D(u,1,-2,0),kadd(KRANC_GFOFFSET3D(u,-2,1,0),KRANC_GFOFFSET3D(u,2,-1,0)))),kmadd(ToReal(-405),kadd(KRANC_GFOFFSET3D(u,-1,-2,0),kadd(KRANC_GFOFFSET3D(u,1,2,0),kadd(KRANC_GFOFFSET3D(u,-2,-1,0),KRANC_GFOFFSET3D(u,2,1,0)))),kmadd(ToReal(-81),kadd(KRANC_GFOFFSET3D(u,-2,2,0),KRANC_GFOFFSET3D(u,2,-2,0)),kmadd(ToReal(81),kadd(KRANC_GFOFFSET3D(u,-2,-2,0),KRANC_GFOFFSET3D(u,2,2,0)),kmadd(ToReal(-45),kadd(KRANC_GFOFFSET3D(u,-1,3,0),kadd(KRANC_GFOFFSET3D(u,1,-3,0),kadd(KRANC_GFOFFSET3D(u,-3,1,0),KRANC_GFOFFSET3D(u,3,-1,0)))),kmadd(ToReal(45),kadd(KRANC_GFOFFSET3D(u,-1,-3,0),kadd(KRANC_GFOFFSET3D(u,1,3,0),kadd(KRANC_GFOFFSET3D(u,-3,-1,0),KRANC_GFOFFSET3D(u,3,1,0)))),kmadd(ToReal(9),kadd(KRANC_GFOFFSET3D(u,-2,3,0),kadd(KRANC_GFOFFSET3D(u,2,-3,0),kadd(KRANC_GFOFFSET3D(u,-3,2,0),KRANC_GFOFFSET3D(u,3,-2,0)))),kmadd(ToReal(-9),kadd(KRANC_GFOFFSET3D(u,-2,-3,0),kadd(KRANC_GFOFFSET3D(u,2,3,0),kadd(KRANC_GFOFFSET3D(u,-3,-2,0),KRANC_GFOFFSET3D(u,3,2,0)))),kadd(KRANC_GFOFFSET3D(u,-3,-3,0),ksub(ksub(KRANC_GFOFFSET3D(u,3,3,0),KRANC_GFOFFSET3D(u,3,-3,0)),KRANC_GFOFFSET3D(u,-3,3,0)))))))))))))))
#else
#  define PDstandardNthfdOrder621(u) (PDstandardNthfdOrder621_impl(u,p1o3600dxdy,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder621_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o3600dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder621_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o3600dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardNthfdOrder612_impl(u, p1o3600dxdy, cdj, cdk);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder623(u) (kmul(p1o3600dydz,kmadd(ToReal(-2025),kadd(KRANC_GFOFFSET3D(u,0,-1,1),KRANC_GFOFFSET3D(u,0,1,-1)),kmadd(ToReal(2025),kadd(KRANC_GFOFFSET3D(u,0,-1,-1),KRANC_GFOFFSET3D(u,0,1,1)),kmadd(ToReal(405),kadd(KRANC_GFOFFSET3D(u,0,-1,2),kadd(KRANC_GFOFFSET3D(u,0,1,-2),kadd(KRANC_GFOFFSET3D(u,0,-2,1),KRANC_GFOFFSET3D(u,0,2,-1)))),kmadd(ToReal(-405),kadd(KRANC_GFOFFSET3D(u,0,-1,-2),kadd(KRANC_GFOFFSET3D(u,0,1,2),kadd(KRANC_GFOFFSET3D(u,0,-2,-1),KRANC_GFOFFSET3D(u,0,2,1)))),kmadd(ToReal(-81),kadd(KRANC_GFOFFSET3D(u,0,-2,2),KRANC_GFOFFSET3D(u,0,2,-2)),kmadd(ToReal(81),kadd(KRANC_GFOFFSET3D(u,0,-2,-2),KRANC_GFOFFSET3D(u,0,2,2)),kmadd(ToReal(-45),kadd(KRANC_GFOFFSET3D(u,0,-1,3),kadd(KRANC_GFOFFSET3D(u,0,1,-3),kadd(KRANC_GFOFFSET3D(u,0,-3,1),KRANC_GFOFFSET3D(u,0,3,-1)))),kmadd(ToReal(45),kadd(KRANC_GFOFFSET3D(u,0,-1,-3),kadd(KRANC_GFOFFSET3D(u,0,1,3),kadd(KRANC_GFOFFSET3D(u,0,-3,-1),KRANC_GFOFFSET3D(u,0,3,1)))),kmadd(ToReal(9),kadd(KRANC_GFOFFSET3D(u,0,-2,3),kadd(KRANC_GFOFFSET3D(u,0,2,-3),kadd(KRANC_GFOFFSET3D(u,0,-3,2),KRANC_GFOFFSET3D(u,0,3,-2)))),kmadd(ToReal(-9),kadd(KRANC_GFOFFSET3D(u,0,-2,-3),kadd(KRANC_GFOFFSET3D(u,0,2,3),kadd(KRANC_GFOFFSET3D(u,0,-3,-2),KRANC_GFOFFSET3D(u,0,3,2)))),kadd(KRANC_GFOFFSET3D(u,0,-3,-3),ksub(ksub(KRANC_GFOFFSET3D(u,0,3,3),KRANC_GFOFFSET3D(u,0,3,-3)),KRANC_GFOFFSET3D(u,0,-3,3)))))))))))))))
#else
#  define PDstandardNthfdOrder623(u) (PDstandardNthfdOrder623_impl(u,p1o3600dydz,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder623_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o3600dydz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder623_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o3600dydz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o3600dydz,kmadd(ToReal(-2025),kadd(KRANC_GFOFFSET3D(u,0,-1,1),KRANC_GFOFFSET3D(u,0,1,-1)),kmadd(ToReal(2025),kadd(KRANC_GFOFFSET3D(u,0,-1,-1),KRANC_GFOFFSET3D(u,0,1,1)),kmadd(ToReal(405),kadd(KRANC_GFOFFSET3D(u,0,-1,2),kadd(KRANC_GFOFFSET3D(u,0,1,-2),kadd(KRANC_GFOFFSET3D(u,0,-2,1),KRANC_GFOFFSET3D(u,0,2,-1)))),kmadd(ToReal(-405),kadd(KRANC_GFOFFSET3D(u,0,-1,-2),kadd(KRANC_GFOFFSET3D(u,0,1,2),kadd(KRANC_GFOFFSET3D(u,0,-2,-1),KRANC_GFOFFSET3D(u,0,2,1)))),kmadd(ToReal(-81),kadd(KRANC_GFOFFSET3D(u,0,-2,2),KRANC_GFOFFSET3D(u,0,2,-2)),kmadd(ToReal(81),kadd(KRANC_GFOFFSET3D(u,0,-2,-2),KRANC_GFOFFSET3D(u,0,2,2)),kmadd(ToReal(-45),kadd(KRANC_GFOFFSET3D(u,0,-1,3),kadd(KRANC_GFOFFSET3D(u,0,1,-3),kadd(KRANC_GFOFFSET3D(u,0,-3,1),KRANC_GFOFFSET3D(u,0,3,-1)))),kmadd(ToReal(45),kadd(KRANC_GFOFFSET3D(u,0,-1,-3),kadd(KRANC_GFOFFSET3D(u,0,1,3),kadd(KRANC_GFOFFSET3D(u,0,-3,-1),KRANC_GFOFFSET3D(u,0,3,1)))),kmadd(ToReal(9),kadd(KRANC_GFOFFSET3D(u,0,-2,3),kadd(KRANC_GFOFFSET3D(u,0,2,-3),kadd(KRANC_GFOFFSET3D(u,0,-3,2),KRANC_GFOFFSET3D(u,0,3,-2)))),kmadd(ToReal(-9),kadd(KRANC_GFOFFSET3D(u,0,-2,-3),kadd(KRANC_GFOFFSET3D(u,0,2,3),kadd(KRANC_GFOFFSET3D(u,0,-3,-2),KRANC_GFOFFSET3D(u,0,3,2)))),kadd(KRANC_GFOFFSET3D(u,0,-3,-3),ksub(ksub(KRANC_GFOFFSET3D(u,0,3,3),KRANC_GFOFFSET3D(u,0,3,-3)),KRANC_GFOFFSET3D(u,0,-3,3))))))))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder631(u) (kmul(p1o3600dxdz,kmadd(ToReal(-2025),kadd(KRANC_GFOFFSET3D(u,-1,0,1),KRANC_GFOFFSET3D(u,1,0,-1)),kmadd(ToReal(2025),kadd(KRANC_GFOFFSET3D(u,-1,0,-1),KRANC_GFOFFSET3D(u,1,0,1)),kmadd(ToReal(405),kadd(KRANC_GFOFFSET3D(u,-1,0,2),kadd(KRANC_GFOFFSET3D(u,1,0,-2),kadd(KRANC_GFOFFSET3D(u,-2,0,1),KRANC_GFOFFSET3D(u,2,0,-1)))),kmadd(ToReal(-405),kadd(KRANC_GFOFFSET3D(u,-1,0,-2),kadd(KRANC_GFOFFSET3D(u,1,0,2),kadd(KRANC_GFOFFSET3D(u,-2,0,-1),KRANC_GFOFFSET3D(u,2,0,1)))),kmadd(ToReal(-81),kadd(KRANC_GFOFFSET3D(u,-2,0,2),KRANC_GFOFFSET3D(u,2,0,-2)),kmadd(ToReal(81),kadd(KRANC_GFOFFSET3D(u,-2,0,-2),KRANC_GFOFFSET3D(u,2,0,2)),kmadd(ToReal(-45),kadd(KRANC_GFOFFSET3D(u,-1,0,3),kadd(KRANC_GFOFFSET3D(u,1,0,-3),kadd(KRANC_GFOFFSET3D(u,-3,0,1),KRANC_GFOFFSET3D(u,3,0,-1)))),kmadd(ToReal(45),kadd(KRANC_GFOFFSET3D(u,-1,0,-3),kadd(KRANC_GFOFFSET3D(u,1,0,3),kadd(KRANC_GFOFFSET3D(u,-3,0,-1),KRANC_GFOFFSET3D(u,3,0,1)))),kmadd(ToReal(9),kadd(KRANC_GFOFFSET3D(u,-2,0,3),kadd(KRANC_GFOFFSET3D(u,2,0,-3),kadd(KRANC_GFOFFSET3D(u,-3,0,2),KRANC_GFOFFSET3D(u,3,0,-2)))),kmadd(ToReal(-9),kadd(KRANC_GFOFFSET3D(u,-2,0,-3),kadd(KRANC_GFOFFSET3D(u,2,0,3),kadd(KRANC_GFOFFSET3D(u,-3,0,-2),KRANC_GFOFFSET3D(u,3,0,2)))),kadd(KRANC_GFOFFSET3D(u,-3,0,-3),ksub(ksub(KRANC_GFOFFSET3D(u,3,0,3),KRANC_GFOFFSET3D(u,3,0,-3)),KRANC_GFOFFSET3D(u,-3,0,3)))))))))))))))
#else
#  define PDstandardNthfdOrder631(u) (PDstandardNthfdOrder631_impl(u,p1o3600dxdz,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder631_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o3600dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder631_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o3600dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardNthfdOrder612_impl(u, p1o3600dxdz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder632(u) (kmul(p1o3600dydz,kmadd(ToReal(-2025),kadd(KRANC_GFOFFSET3D(u,0,-1,1),KRANC_GFOFFSET3D(u,0,1,-1)),kmadd(ToReal(2025),kadd(KRANC_GFOFFSET3D(u,0,-1,-1),KRANC_GFOFFSET3D(u,0,1,1)),kmadd(ToReal(405),kadd(KRANC_GFOFFSET3D(u,0,-1,2),kadd(KRANC_GFOFFSET3D(u,0,1,-2),kadd(KRANC_GFOFFSET3D(u,0,-2,1),KRANC_GFOFFSET3D(u,0,2,-1)))),kmadd(ToReal(-405),kadd(KRANC_GFOFFSET3D(u,0,-1,-2),kadd(KRANC_GFOFFSET3D(u,0,1,2),kadd(KRANC_GFOFFSET3D(u,0,-2,-1),KRANC_GFOFFSET3D(u,0,2,1)))),kmadd(ToReal(-81),kadd(KRANC_GFOFFSET3D(u,0,-2,2),KRANC_GFOFFSET3D(u,0,2,-2)),kmadd(ToReal(81),kadd(KRANC_GFOFFSET3D(u,0,-2,-2),KRANC_GFOFFSET3D(u,0,2,2)),kmadd(ToReal(-45),kadd(KRANC_GFOFFSET3D(u,0,-1,3),kadd(KRANC_GFOFFSET3D(u,0,1,-3),kadd(KRANC_GFOFFSET3D(u,0,-3,1),KRANC_GFOFFSET3D(u,0,3,-1)))),kmadd(ToReal(45),kadd(KRANC_GFOFFSET3D(u,0,-1,-3),kadd(KRANC_GFOFFSET3D(u,0,1,3),kadd(KRANC_GFOFFSET3D(u,0,-3,-1),KRANC_GFOFFSET3D(u,0,3,1)))),kmadd(ToReal(9),kadd(KRANC_GFOFFSET3D(u,0,-2,3),kadd(KRANC_GFOFFSET3D(u,0,2,-3),kadd(KRANC_GFOFFSET3D(u,0,-3,2),KRANC_GFOFFSET3D(u,0,3,-2)))),kmadd(ToReal(-9),kadd(KRANC_GFOFFSET3D(u,0,-2,-3),kadd(KRANC_GFOFFSET3D(u,0,2,3),kadd(KRANC_GFOFFSET3D(u,0,-3,-2),KRANC_GFOFFSET3D(u,0,3,2)))),kadd(KRANC_GFOFFSET3D(u,0,-3,-3),ksub(ksub(KRANC_GFOFFSET3D(u,0,3,3),KRANC_GFOFFSET3D(u,0,3,-3)),KRANC_GFOFFSET3D(u,0,-3,3)))))))))))))))
#else
#  define PDstandardNthfdOrder632(u) (PDstandardNthfdOrder632_impl(u,p1o3600dydz,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder632_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o3600dydz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder632_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o3600dydz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardNthfdOrder623_impl(u, p1o3600dydz, cdj, cdk);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder812(u) (kmul(p1o705600dxdy,kmadd(ToReal(-451584),kadd(KRANC_GFOFFSET3D(u,-1,1,0),KRANC_GFOFFSET3D(u,1,-1,0)),kmadd(ToReal(451584),kadd(KRANC_GFOFFSET3D(u,-1,-1,0),KRANC_GFOFFSET3D(u,1,1,0)),kmadd(ToReal(112896),kadd(KRANC_GFOFFSET3D(u,-1,2,0),kadd(KRANC_GFOFFSET3D(u,1,-2,0),kadd(KRANC_GFOFFSET3D(u,-2,1,0),KRANC_GFOFFSET3D(u,2,-1,0)))),kmadd(ToReal(-112896),kadd(KRANC_GFOFFSET3D(u,-1,-2,0),kadd(KRANC_GFOFFSET3D(u,1,2,0),kadd(KRANC_GFOFFSET3D(u,-2,-1,0),KRANC_GFOFFSET3D(u,2,1,0)))),kmadd(ToReal(-28224),kadd(KRANC_GFOFFSET3D(u,-2,2,0),KRANC_GFOFFSET3D(u,2,-2,0)),kmadd(ToReal(28224),kadd(KRANC_GFOFFSET3D(u,-2,-2,0),KRANC_GFOFFSET3D(u,2,2,0)),kmadd(ToReal(-21504),kadd(KRANC_GFOFFSET3D(u,-1,3,0),kadd(KRANC_GFOFFSET3D(u,1,-3,0),kadd(KRANC_GFOFFSET3D(u,-3,1,0),KRANC_GFOFFSET3D(u,3,-1,0)))),kmadd(ToReal(21504),kadd(KRANC_GFOFFSET3D(u,-1,-3,0),kadd(KRANC_GFOFFSET3D(u,1,3,0),kadd(KRANC_GFOFFSET3D(u,-3,-1,0),KRANC_GFOFFSET3D(u,3,1,0)))),kmadd(ToReal(5376),kadd(KRANC_GFOFFSET3D(u,-2,3,0),kadd(KRANC_GFOFFSET3D(u,2,-3,0),kadd(KRANC_GFOFFSET3D(u,-3,2,0),KRANC_GFOFFSET3D(u,3,-2,0)))),kmadd(ToReal(-5376),kadd(KRANC_GFOFFSET3D(u,-2,-3,0),kadd(KRANC_GFOFFSET3D(u,2,3,0),kadd(KRANC_GFOFFSET3D(u,-3,-2,0),KRANC_GFOFFSET3D(u,3,2,0)))),kmadd(ToReal(-1024),kadd(KRANC_GFOFFSET3D(u,-3,3,0),KRANC_GFOFFSET3D(u,3,-3,0)),kmadd(ToReal(1024),kadd(KRANC_GFOFFSET3D(u,-3,-3,0),KRANC_GFOFFSET3D(u,3,3,0)),kmadd(ToReal(2016),kadd(KRANC_GFOFFSET3D(u,-1,4,0),kadd(KRANC_GFOFFSET3D(u,1,-4,0),kadd(KRANC_GFOFFSET3D(u,-4,1,0),KRANC_GFOFFSET3D(u,4,-1,0)))),kmadd(ToReal(-2016),kadd(KRANC_GFOFFSET3D(u,-1,-4,0),kadd(KRANC_GFOFFSET3D(u,1,4,0),kadd(KRANC_GFOFFSET3D(u,-4,-1,0),KRANC_GFOFFSET3D(u,4,1,0)))),kmadd(ToReal(-504),kadd(KRANC_GFOFFSET3D(u,-2,4,0),kadd(KRANC_GFOFFSET3D(u,2,-4,0),kadd(KRANC_GFOFFSET3D(u,-4,2,0),KRANC_GFOFFSET3D(u,4,-2,0)))),kmadd(ToReal(504),kadd(KRANC_GFOFFSET3D(u,-2,-4,0),kadd(KRANC_GFOFFSET3D(u,2,4,0),kadd(KRANC_GFOFFSET3D(u,-4,-2,0),KRANC_GFOFFSET3D(u,4,2,0)))),kmadd(ToReal(96),kadd(KRANC_GFOFFSET3D(u,-3,4,0),kadd(KRANC_GFOFFSET3D(u,3,-4,0),kadd(KRANC_GFOFFSET3D(u,-4,3,0),KRANC_GFOFFSET3D(u,4,-3,0)))),kmadd(ToReal(-96),kadd(KRANC_GFOFFSET3D(u,-3,-4,0),kadd(KRANC_GFOFFSET3D(u,3,4,0),kadd(KRANC_GFOFFSET3D(u,-4,-3,0),KRANC_GFOFFSET3D(u,4,3,0)))),kmadd(ToReal(-9),kadd(KRANC_GFOFFSET3D(u,-4,4,0),KRANC_GFOFFSET3D(u,4,-4,0)),kmul(ToReal(9),kadd(KRANC_GFOFFSET3D(u,-4,-4,0),KRANC_GFOFFSET3D(u,4,4,0))))))))))))))))))))))))
#else
#  define PDstandardNthfdOrder812(u) (PDstandardNthfdOrder812_impl(u,p1o705600dxdy,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder812_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o705600dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder812_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o705600dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o705600dxdy,kmadd(ToReal(-451584),kadd(KRANC_GFOFFSET3D(u,-1,1,0),KRANC_GFOFFSET3D(u,1,-1,0)),kmadd(ToReal(451584),kadd(KRANC_GFOFFSET3D(u,-1,-1,0),KRANC_GFOFFSET3D(u,1,1,0)),kmadd(ToReal(112896),kadd(KRANC_GFOFFSET3D(u,-1,2,0),kadd(KRANC_GFOFFSET3D(u,1,-2,0),kadd(KRANC_GFOFFSET3D(u,-2,1,0),KRANC_GFOFFSET3D(u,2,-1,0)))),kmadd(ToReal(-112896),kadd(KRANC_GFOFFSET3D(u,-1,-2,0),kadd(KRANC_GFOFFSET3D(u,1,2,0),kadd(KRANC_GFOFFSET3D(u,-2,-1,0),KRANC_GFOFFSET3D(u,2,1,0)))),kmadd(ToReal(-28224),kadd(KRANC_GFOFFSET3D(u,-2,2,0),KRANC_GFOFFSET3D(u,2,-2,0)),kmadd(ToReal(28224),kadd(KRANC_GFOFFSET3D(u,-2,-2,0),KRANC_GFOFFSET3D(u,2,2,0)),kmadd(ToReal(-21504),kadd(KRANC_GFOFFSET3D(u,-1,3,0),kadd(KRANC_GFOFFSET3D(u,1,-3,0),kadd(KRANC_GFOFFSET3D(u,-3,1,0),KRANC_GFOFFSET3D(u,3,-1,0)))),kmadd(ToReal(21504),kadd(KRANC_GFOFFSET3D(u,-1,-3,0),kadd(KRANC_GFOFFSET3D(u,1,3,0),kadd(KRANC_GFOFFSET3D(u,-3,-1,0),KRANC_GFOFFSET3D(u,3,1,0)))),kmadd(ToReal(5376),kadd(KRANC_GFOFFSET3D(u,-2,3,0),kadd(KRANC_GFOFFSET3D(u,2,-3,0),kadd(KRANC_GFOFFSET3D(u,-3,2,0),KRANC_GFOFFSET3D(u,3,-2,0)))),kmadd(ToReal(-5376),kadd(KRANC_GFOFFSET3D(u,-2,-3,0),kadd(KRANC_GFOFFSET3D(u,2,3,0),kadd(KRANC_GFOFFSET3D(u,-3,-2,0),KRANC_GFOFFSET3D(u,3,2,0)))),kmadd(ToReal(-1024),kadd(KRANC_GFOFFSET3D(u,-3,3,0),KRANC_GFOFFSET3D(u,3,-3,0)),kmadd(ToReal(1024),kadd(KRANC_GFOFFSET3D(u,-3,-3,0),KRANC_GFOFFSET3D(u,3,3,0)),kmadd(ToReal(2016),kadd(KRANC_GFOFFSET3D(u,-1,4,0),kadd(KRANC_GFOFFSET3D(u,1,-4,0),kadd(KRANC_GFOFFSET3D(u,-4,1,0),KRANC_GFOFFSET3D(u,4,-1,0)))),kmadd(ToReal(-2016),kadd(KRANC_GFOFFSET3D(u,-1,-4,0),kadd(KRANC_GFOFFSET3D(u,1,4,0),kadd(KRANC_GFOFFSET3D(u,-4,-1,0),KRANC_GFOFFSET3D(u,4,1,0)))),kmadd(ToReal(-504),kadd(KRANC_GFOFFSET3D(u,-2,4,0),kadd(KRANC_GFOFFSET3D(u,2,-4,0),kadd(KRANC_GFOFFSET3D(u,-4,2,0),KRANC_GFOFFSET3D(u,4,-2,0)))),kmadd(ToReal(504),kadd(KRANC_GFOFFSET3D(u,-2,-4,0),kadd(KRANC_GFOFFSET3D(u,2,4,0),kadd(KRANC_GFOFFSET3D(u,-4,-2,0),KRANC_GFOFFSET3D(u,4,2,0)))),kmadd(ToReal(96),kadd(KRANC_GFOFFSET3D(u,-3,4,0),kadd(KRANC_GFOFFSET3D(u,3,-4,0),kadd(KRANC_GFOFFSET3D(u,-4,3,0),KRANC_GFOFFSET3D(u,4,-3,0)))),kmadd(ToReal(-96),kadd(KRANC_GFOFFSET3D(u,-3,-4,0),kadd(KRANC_GFOFFSET3D(u,3,4,0),kadd(KRANC_GFOFFSET3D(u,-4,-3,0),KRANC_GFOFFSET3D(u,4,3,0)))),kmadd(ToReal(-9),kadd(KRANC_GFOFFSET3D(u,-4,4,0),KRANC_GFOFFSET3D(u,4,-4,0)),kmul(ToReal(9),kadd(KRANC_GFOFFSET3D(u,-4,-4,0),KRANC_GFOFFSET3D(u,4,4,0)))))))))))))))))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder813(u) (kmul(p1o705600dxdz,kmadd(ToReal(-451584),kadd(KRANC_GFOFFSET3D(u,-1,0,1),KRANC_GFOFFSET3D(u,1,0,-1)),kmadd(ToReal(451584),kadd(KRANC_GFOFFSET3D(u,-1,0,-1),KRANC_GFOFFSET3D(u,1,0,1)),kmadd(ToReal(112896),kadd(KRANC_GFOFFSET3D(u,-1,0,2),kadd(KRANC_GFOFFSET3D(u,1,0,-2),kadd(KRANC_GFOFFSET3D(u,-2,0,1),KRANC_GFOFFSET3D(u,2,0,-1)))),kmadd(ToReal(-112896),kadd(KRANC_GFOFFSET3D(u,-1,0,-2),kadd(KRANC_GFOFFSET3D(u,1,0,2),kadd(KRANC_GFOFFSET3D(u,-2,0,-1),KRANC_GFOFFSET3D(u,2,0,1)))),kmadd(ToReal(-28224),kadd(KRANC_GFOFFSET3D(u,-2,0,2),KRANC_GFOFFSET3D(u,2,0,-2)),kmadd(ToReal(28224),kadd(KRANC_GFOFFSET3D(u,-2,0,-2),KRANC_GFOFFSET3D(u,2,0,2)),kmadd(ToReal(-21504),kadd(KRANC_GFOFFSET3D(u,-1,0,3),kadd(KRANC_GFOFFSET3D(u,1,0,-3),kadd(KRANC_GFOFFSET3D(u,-3,0,1),KRANC_GFOFFSET3D(u,3,0,-1)))),kmadd(ToReal(21504),kadd(KRANC_GFOFFSET3D(u,-1,0,-3),kadd(KRANC_GFOFFSET3D(u,1,0,3),kadd(KRANC_GFOFFSET3D(u,-3,0,-1),KRANC_GFOFFSET3D(u,3,0,1)))),kmadd(ToReal(5376),kadd(KRANC_GFOFFSET3D(u,-2,0,3),kadd(KRANC_GFOFFSET3D(u,2,0,-3),kadd(KRANC_GFOFFSET3D(u,-3,0,2),KRANC_GFOFFSET3D(u,3,0,-2)))),kmadd(ToReal(-5376),kadd(KRANC_GFOFFSET3D(u,-2,0,-3),kadd(KRANC_GFOFFSET3D(u,2,0,3),kadd(KRANC_GFOFFSET3D(u,-3,0,-2),KRANC_GFOFFSET3D(u,3,0,2)))),kmadd(ToReal(-1024),kadd(KRANC_GFOFFSET3D(u,-3,0,3),KRANC_GFOFFSET3D(u,3,0,-3)),kmadd(ToReal(1024),kadd(KRANC_GFOFFSET3D(u,-3,0,-3),KRANC_GFOFFSET3D(u,3,0,3)),kmadd(ToReal(2016),kadd(KRANC_GFOFFSET3D(u,-1,0,4),kadd(KRANC_GFOFFSET3D(u,1,0,-4),kadd(KRANC_GFOFFSET3D(u,-4,0,1),KRANC_GFOFFSET3D(u,4,0,-1)))),kmadd(ToReal(-2016),kadd(KRANC_GFOFFSET3D(u,-1,0,-4),kadd(KRANC_GFOFFSET3D(u,1,0,4),kadd(KRANC_GFOFFSET3D(u,-4,0,-1),KRANC_GFOFFSET3D(u,4,0,1)))),kmadd(ToReal(-504),kadd(KRANC_GFOFFSET3D(u,-2,0,4),kadd(KRANC_GFOFFSET3D(u,2,0,-4),kadd(KRANC_GFOFFSET3D(u,-4,0,2),KRANC_GFOFFSET3D(u,4,0,-2)))),kmadd(ToReal(504),kadd(KRANC_GFOFFSET3D(u,-2,0,-4),kadd(KRANC_GFOFFSET3D(u,2,0,4),kadd(KRANC_GFOFFSET3D(u,-4,0,-2),KRANC_GFOFFSET3D(u,4,0,2)))),kmadd(ToReal(96),kadd(KRANC_GFOFFSET3D(u,-3,0,4),kadd(KRANC_GFOFFSET3D(u,3,0,-4),kadd(KRANC_GFOFFSET3D(u,-4,0,3),KRANC_GFOFFSET3D(u,4,0,-3)))),kmadd(ToReal(-96),kadd(KRANC_GFOFFSET3D(u,-3,0,-4),kadd(KRANC_GFOFFSET3D(u,3,0,4),kadd(KRANC_GFOFFSET3D(u,-4,0,-3),KRANC_GFOFFSET3D(u,4,0,3)))),kmadd(ToReal(-9),kadd(KRANC_GFOFFSET3D(u,-4,0,4),KRANC_GFOFFSET3D(u,4,0,-4)),kmul(ToReal(9),kadd(KRANC_GFOFFSET3D(u,-4,0,-4),KRANC_GFOFFSET3D(u,4,0,4))))))))))))))))))))))))
#else
#  define PDstandardNthfdOrder813(u) (PDstandardNthfdOrder813_impl(u,p1o705600dxdz,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder813_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o705600dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder813_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o705600dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardNthfdOrder812_impl(u, p1o705600dxdz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder821(u) (kmul(p1o705600dxdy,kmadd(ToReal(-451584),kadd(KRANC_GFOFFSET3D(u,-1,1,0),KRANC_GFOFFSET3D(u,1,-1,0)),kmadd(ToReal(451584),kadd(KRANC_GFOFFSET3D(u,-1,-1,0),KRANC_GFOFFSET3D(u,1,1,0)),kmadd(ToReal(112896),kadd(KRANC_GFOFFSET3D(u,-1,2,0),kadd(KRANC_GFOFFSET3D(u,1,-2,0),kadd(KRANC_GFOFFSET3D(u,-2,1,0),KRANC_GFOFFSET3D(u,2,-1,0)))),kmadd(ToReal(-112896),kadd(KRANC_GFOFFSET3D(u,-1,-2,0),kadd(KRANC_GFOFFSET3D(u,1,2,0),kadd(KRANC_GFOFFSET3D(u,-2,-1,0),KRANC_GFOFFSET3D(u,2,1,0)))),kmadd(ToReal(-28224),kadd(KRANC_GFOFFSET3D(u,-2,2,0),KRANC_GFOFFSET3D(u,2,-2,0)),kmadd(ToReal(28224),kadd(KRANC_GFOFFSET3D(u,-2,-2,0),KRANC_GFOFFSET3D(u,2,2,0)),kmadd(ToReal(-21504),kadd(KRANC_GFOFFSET3D(u,-1,3,0),kadd(KRANC_GFOFFSET3D(u,1,-3,0),kadd(KRANC_GFOFFSET3D(u,-3,1,0),KRANC_GFOFFSET3D(u,3,-1,0)))),kmadd(ToReal(21504),kadd(KRANC_GFOFFSET3D(u,-1,-3,0),kadd(KRANC_GFOFFSET3D(u,1,3,0),kadd(KRANC_GFOFFSET3D(u,-3,-1,0),KRANC_GFOFFSET3D(u,3,1,0)))),kmadd(ToReal(5376),kadd(KRANC_GFOFFSET3D(u,-2,3,0),kadd(KRANC_GFOFFSET3D(u,2,-3,0),kadd(KRANC_GFOFFSET3D(u,-3,2,0),KRANC_GFOFFSET3D(u,3,-2,0)))),kmadd(ToReal(-5376),kadd(KRANC_GFOFFSET3D(u,-2,-3,0),kadd(KRANC_GFOFFSET3D(u,2,3,0),kadd(KRANC_GFOFFSET3D(u,-3,-2,0),KRANC_GFOFFSET3D(u,3,2,0)))),kmadd(ToReal(-1024),kadd(KRANC_GFOFFSET3D(u,-3,3,0),KRANC_GFOFFSET3D(u,3,-3,0)),kmadd(ToReal(1024),kadd(KRANC_GFOFFSET3D(u,-3,-3,0),KRANC_GFOFFSET3D(u,3,3,0)),kmadd(ToReal(2016),kadd(KRANC_GFOFFSET3D(u,-1,4,0),kadd(KRANC_GFOFFSET3D(u,1,-4,0),kadd(KRANC_GFOFFSET3D(u,-4,1,0),KRANC_GFOFFSET3D(u,4,-1,0)))),kmadd(ToReal(-2016),kadd(KRANC_GFOFFSET3D(u,-1,-4,0),kadd(KRANC_GFOFFSET3D(u,1,4,0),kadd(KRANC_GFOFFSET3D(u,-4,-1,0),KRANC_GFOFFSET3D(u,4,1,0)))),kmadd(ToReal(-504),kadd(KRANC_GFOFFSET3D(u,-2,4,0),kadd(KRANC_GFOFFSET3D(u,2,-4,0),kadd(KRANC_GFOFFSET3D(u,-4,2,0),KRANC_GFOFFSET3D(u,4,-2,0)))),kmadd(ToReal(504),kadd(KRANC_GFOFFSET3D(u,-2,-4,0),kadd(KRANC_GFOFFSET3D(u,2,4,0),kadd(KRANC_GFOFFSET3D(u,-4,-2,0),KRANC_GFOFFSET3D(u,4,2,0)))),kmadd(ToReal(96),kadd(KRANC_GFOFFSET3D(u,-3,4,0),kadd(KRANC_GFOFFSET3D(u,3,-4,0),kadd(KRANC_GFOFFSET3D(u,-4,3,0),KRANC_GFOFFSET3D(u,4,-3,0)))),kmadd(ToReal(-96),kadd(KRANC_GFOFFSET3D(u,-3,-4,0),kadd(KRANC_GFOFFSET3D(u,3,4,0),kadd(KRANC_GFOFFSET3D(u,-4,-3,0),KRANC_GFOFFSET3D(u,4,3,0)))),kmadd(ToReal(-9),kadd(KRANC_GFOFFSET3D(u,-4,4,0),KRANC_GFOFFSET3D(u,4,-4,0)),kmul(ToReal(9),kadd(KRANC_GFOFFSET3D(u,-4,-4,0),KRANC_GFOFFSET3D(u,4,4,0))))))))))))))))))))))))
#else
#  define PDstandardNthfdOrder821(u) (PDstandardNthfdOrder821_impl(u,p1o705600dxdy,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder821_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o705600dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder821_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o705600dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardNthfdOrder812_impl(u, p1o705600dxdy, cdj, cdk);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder823(u) (kmul(p1o705600dydz,kmadd(ToReal(-451584),kadd(KRANC_GFOFFSET3D(u,0,-1,1),KRANC_GFOFFSET3D(u,0,1,-1)),kmadd(ToReal(451584),kadd(KRANC_GFOFFSET3D(u,0,-1,-1),KRANC_GFOFFSET3D(u,0,1,1)),kmadd(ToReal(112896),kadd(KRANC_GFOFFSET3D(u,0,-1,2),kadd(KRANC_GFOFFSET3D(u,0,1,-2),kadd(KRANC_GFOFFSET3D(u,0,-2,1),KRANC_GFOFFSET3D(u,0,2,-1)))),kmadd(ToReal(-112896),kadd(KRANC_GFOFFSET3D(u,0,-1,-2),kadd(KRANC_GFOFFSET3D(u,0,1,2),kadd(KRANC_GFOFFSET3D(u,0,-2,-1),KRANC_GFOFFSET3D(u,0,2,1)))),kmadd(ToReal(-28224),kadd(KRANC_GFOFFSET3D(u,0,-2,2),KRANC_GFOFFSET3D(u,0,2,-2)),kmadd(ToReal(28224),kadd(KRANC_GFOFFSET3D(u,0,-2,-2),KRANC_GFOFFSET3D(u,0,2,2)),kmadd(ToReal(-21504),kadd(KRANC_GFOFFSET3D(u,0,-1,3),kadd(KRANC_GFOFFSET3D(u,0,1,-3),kadd(KRANC_GFOFFSET3D(u,0,-3,1),KRANC_GFOFFSET3D(u,0,3,-1)))),kmadd(ToReal(21504),kadd(KRANC_GFOFFSET3D(u,0,-1,-3),kadd(KRANC_GFOFFSET3D(u,0,1,3),kadd(KRANC_GFOFFSET3D(u,0,-3,-1),KRANC_GFOFFSET3D(u,0,3,1)))),kmadd(ToReal(5376),kadd(KRANC_GFOFFSET3D(u,0,-2,3),kadd(KRANC_GFOFFSET3D(u,0,2,-3),kadd(KRANC_GFOFFSET3D(u,0,-3,2),KRANC_GFOFFSET3D(u,0,3,-2)))),kmadd(ToReal(-5376),kadd(KRANC_GFOFFSET3D(u,0,-2,-3),kadd(KRANC_GFOFFSET3D(u,0,2,3),kadd(KRANC_GFOFFSET3D(u,0,-3,-2),KRANC_GFOFFSET3D(u,0,3,2)))),kmadd(ToReal(-1024),kadd(KRANC_GFOFFSET3D(u,0,-3,3),KRANC_GFOFFSET3D(u,0,3,-3)),kmadd(ToReal(1024),kadd(KRANC_GFOFFSET3D(u,0,-3,-3),KRANC_GFOFFSET3D(u,0,3,3)),kmadd(ToReal(2016),kadd(KRANC_GFOFFSET3D(u,0,-1,4),kadd(KRANC_GFOFFSET3D(u,0,1,-4),kadd(KRANC_GFOFFSET3D(u,0,-4,1),KRANC_GFOFFSET3D(u,0,4,-1)))),kmadd(ToReal(-2016),kadd(KRANC_GFOFFSET3D(u,0,-1,-4),kadd(KRANC_GFOFFSET3D(u,0,1,4),kadd(KRANC_GFOFFSET3D(u,0,-4,-1),KRANC_GFOFFSET3D(u,0,4,1)))),kmadd(ToReal(-504),kadd(KRANC_GFOFFSET3D(u,0,-2,4),kadd(KRANC_GFOFFSET3D(u,0,2,-4),kadd(KRANC_GFOFFSET3D(u,0,-4,2),KRANC_GFOFFSET3D(u,0,4,-2)))),kmadd(ToReal(504),kadd(KRANC_GFOFFSET3D(u,0,-2,-4),kadd(KRANC_GFOFFSET3D(u,0,2,4),kadd(KRANC_GFOFFSET3D(u,0,-4,-2),KRANC_GFOFFSET3D(u,0,4,2)))),kmadd(ToReal(96),kadd(KRANC_GFOFFSET3D(u,0,-3,4),kadd(KRANC_GFOFFSET3D(u,0,3,-4),kadd(KRANC_GFOFFSET3D(u,0,-4,3),KRANC_GFOFFSET3D(u,0,4,-3)))),kmadd(ToReal(-96),kadd(KRANC_GFOFFSET3D(u,0,-3,-4),kadd(KRANC_GFOFFSET3D(u,0,3,4),kadd(KRANC_GFOFFSET3D(u,0,-4,-3),KRANC_GFOFFSET3D(u,0,4,3)))),kmadd(ToReal(-9),kadd(KRANC_GFOFFSET3D(u,0,-4,4),KRANC_GFOFFSET3D(u,0,4,-4)),kmul(ToReal(9),kadd(KRANC_GFOFFSET3D(u,0,-4,-4),KRANC_GFOFFSET3D(u,0,4,4))))))))))))))))))))))))
#else
#  define PDstandardNthfdOrder823(u) (PDstandardNthfdOrder823_impl(u,p1o705600dydz,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder823_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o705600dydz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder823_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o705600dydz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o705600dydz,kmadd(ToReal(-451584),kadd(KRANC_GFOFFSET3D(u,0,-1,1),KRANC_GFOFFSET3D(u,0,1,-1)),kmadd(ToReal(451584),kadd(KRANC_GFOFFSET3D(u,0,-1,-1),KRANC_GFOFFSET3D(u,0,1,1)),kmadd(ToReal(112896),kadd(KRANC_GFOFFSET3D(u,0,-1,2),kadd(KRANC_GFOFFSET3D(u,0,1,-2),kadd(KRANC_GFOFFSET3D(u,0,-2,1),KRANC_GFOFFSET3D(u,0,2,-1)))),kmadd(ToReal(-112896),kadd(KRANC_GFOFFSET3D(u,0,-1,-2),kadd(KRANC_GFOFFSET3D(u,0,1,2),kadd(KRANC_GFOFFSET3D(u,0,-2,-1),KRANC_GFOFFSET3D(u,0,2,1)))),kmadd(ToReal(-28224),kadd(KRANC_GFOFFSET3D(u,0,-2,2),KRANC_GFOFFSET3D(u,0,2,-2)),kmadd(ToReal(28224),kadd(KRANC_GFOFFSET3D(u,0,-2,-2),KRANC_GFOFFSET3D(u,0,2,2)),kmadd(ToReal(-21504),kadd(KRANC_GFOFFSET3D(u,0,-1,3),kadd(KRANC_GFOFFSET3D(u,0,1,-3),kadd(KRANC_GFOFFSET3D(u,0,-3,1),KRANC_GFOFFSET3D(u,0,3,-1)))),kmadd(ToReal(21504),kadd(KRANC_GFOFFSET3D(u,0,-1,-3),kadd(KRANC_GFOFFSET3D(u,0,1,3),kadd(KRANC_GFOFFSET3D(u,0,-3,-1),KRANC_GFOFFSET3D(u,0,3,1)))),kmadd(ToReal(5376),kadd(KRANC_GFOFFSET3D(u,0,-2,3),kadd(KRANC_GFOFFSET3D(u,0,2,-3),kadd(KRANC_GFOFFSET3D(u,0,-3,2),KRANC_GFOFFSET3D(u,0,3,-2)))),kmadd(ToReal(-5376),kadd(KRANC_GFOFFSET3D(u,0,-2,-3),kadd(KRANC_GFOFFSET3D(u,0,2,3),kadd(KRANC_GFOFFSET3D(u,0,-3,-2),KRANC_GFOFFSET3D(u,0,3,2)))),kmadd(ToReal(-1024),kadd(KRANC_GFOFFSET3D(u,0,-3,3),KRANC_GFOFFSET3D(u,0,3,-3)),kmadd(ToReal(1024),kadd(KRANC_GFOFFSET3D(u,0,-3,-3),KRANC_GFOFFSET3D(u,0,3,3)),kmadd(ToReal(2016),kadd(KRANC_GFOFFSET3D(u,0,-1,4),kadd(KRANC_GFOFFSET3D(u,0,1,-4),kadd(KRANC_GFOFFSET3D(u,0,-4,1),KRANC_GFOFFSET3D(u,0,4,-1)))),kmadd(ToReal(-2016),kadd(KRANC_GFOFFSET3D(u,0,-1,-4),kadd(KRANC_GFOFFSET3D(u,0,1,4),kadd(KRANC_GFOFFSET3D(u,0,-4,-1),KRANC_GFOFFSET3D(u,0,4,1)))),kmadd(ToReal(-504),kadd(KRANC_GFOFFSET3D(u,0,-2,4),kadd(KRANC_GFOFFSET3D(u,0,2,-4),kadd(KRANC_GFOFFSET3D(u,0,-4,2),KRANC_GFOFFSET3D(u,0,4,-2)))),kmadd(ToReal(504),kadd(KRANC_GFOFFSET3D(u,0,-2,-4),kadd(KRANC_GFOFFSET3D(u,0,2,4),kadd(KRANC_GFOFFSET3D(u,0,-4,-2),KRANC_GFOFFSET3D(u,0,4,2)))),kmadd(ToReal(96),kadd(KRANC_GFOFFSET3D(u,0,-3,4),kadd(KRANC_GFOFFSET3D(u,0,3,-4),kadd(KRANC_GFOFFSET3D(u,0,-4,3),KRANC_GFOFFSET3D(u,0,4,-3)))),kmadd(ToReal(-96),kadd(KRANC_GFOFFSET3D(u,0,-3,-4),kadd(KRANC_GFOFFSET3D(u,0,3,4),kadd(KRANC_GFOFFSET3D(u,0,-4,-3),KRANC_GFOFFSET3D(u,0,4,3)))),kmadd(ToReal(-9),kadd(KRANC_GFOFFSET3D(u,0,-4,4),KRANC_GFOFFSET3D(u,0,4,-4)),kmul(ToReal(9),kadd(KRANC_GFOFFSET3D(u,0,-4,-4),KRANC_GFOFFSET3D(u,0,4,4)))))))))))))))))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder831(u) (kmul(p1o705600dxdz,kmadd(ToReal(-451584),kadd(KRANC_GFOFFSET3D(u,-1,0,1),KRANC_GFOFFSET3D(u,1,0,-1)),kmadd(ToReal(451584),kadd(KRANC_GFOFFSET3D(u,-1,0,-1),KRANC_GFOFFSET3D(u,1,0,1)),kmadd(ToReal(112896),kadd(KRANC_GFOFFSET3D(u,-1,0,2),kadd(KRANC_GFOFFSET3D(u,1,0,-2),kadd(KRANC_GFOFFSET3D(u,-2,0,1),KRANC_GFOFFSET3D(u,2,0,-1)))),kmadd(ToReal(-112896),kadd(KRANC_GFOFFSET3D(u,-1,0,-2),kadd(KRANC_GFOFFSET3D(u,1,0,2),kadd(KRANC_GFOFFSET3D(u,-2,0,-1),KRANC_GFOFFSET3D(u,2,0,1)))),kmadd(ToReal(-28224),kadd(KRANC_GFOFFSET3D(u,-2,0,2),KRANC_GFOFFSET3D(u,2,0,-2)),kmadd(ToReal(28224),kadd(KRANC_GFOFFSET3D(u,-2,0,-2),KRANC_GFOFFSET3D(u,2,0,2)),kmadd(ToReal(-21504),kadd(KRANC_GFOFFSET3D(u,-1,0,3),kadd(KRANC_GFOFFSET3D(u,1,0,-3),kadd(KRANC_GFOFFSET3D(u,-3,0,1),KRANC_GFOFFSET3D(u,3,0,-1)))),kmadd(ToReal(21504),kadd(KRANC_GFOFFSET3D(u,-1,0,-3),kadd(KRANC_GFOFFSET3D(u,1,0,3),kadd(KRANC_GFOFFSET3D(u,-3,0,-1),KRANC_GFOFFSET3D(u,3,0,1)))),kmadd(ToReal(5376),kadd(KRANC_GFOFFSET3D(u,-2,0,3),kadd(KRANC_GFOFFSET3D(u,2,0,-3),kadd(KRANC_GFOFFSET3D(u,-3,0,2),KRANC_GFOFFSET3D(u,3,0,-2)))),kmadd(ToReal(-5376),kadd(KRANC_GFOFFSET3D(u,-2,0,-3),kadd(KRANC_GFOFFSET3D(u,2,0,3),kadd(KRANC_GFOFFSET3D(u,-3,0,-2),KRANC_GFOFFSET3D(u,3,0,2)))),kmadd(ToReal(-1024),kadd(KRANC_GFOFFSET3D(u,-3,0,3),KRANC_GFOFFSET3D(u,3,0,-3)),kmadd(ToReal(1024),kadd(KRANC_GFOFFSET3D(u,-3,0,-3),KRANC_GFOFFSET3D(u,3,0,3)),kmadd(ToReal(2016),kadd(KRANC_GFOFFSET3D(u,-1,0,4),kadd(KRANC_GFOFFSET3D(u,1,0,-4),kadd(KRANC_GFOFFSET3D(u,-4,0,1),KRANC_GFOFFSET3D(u,4,0,-1)))),kmadd(ToReal(-2016),kadd(KRANC_GFOFFSET3D(u,-1,0,-4),kadd(KRANC_GFOFFSET3D(u,1,0,4),kadd(KRANC_GFOFFSET3D(u,-4,0,-1),KRANC_GFOFFSET3D(u,4,0,1)))),kmadd(ToReal(-504),kadd(KRANC_GFOFFSET3D(u,-2,0,4),kadd(KRANC_GFOFFSET3D(u,2,0,-4),kadd(KRANC_GFOFFSET3D(u,-4,0,2),KRANC_GFOFFSET3D(u,4,0,-2)))),kmadd(ToReal(504),kadd(KRANC_GFOFFSET3D(u,-2,0,-4),kadd(KRANC_GFOFFSET3D(u,2,0,4),kadd(KRANC_GFOFFSET3D(u,-4,0,-2),KRANC_GFOFFSET3D(u,4,0,2)))),kmadd(ToReal(96),kadd(KRANC_GFOFFSET3D(u,-3,0,4),kadd(KRANC_GFOFFSET3D(u,3,0,-4),kadd(KRANC_GFOFFSET3D(u,-4,0,3),KRANC_GFOFFSET3D(u,4,0,-3)))),kmadd(ToReal(-96),kadd(KRANC_GFOFFSET3D(u,-3,0,-4),kadd(KRANC_GFOFFSET3D(u,3,0,4),kadd(KRANC_GFOFFSET3D(u,-4,0,-3),KRANC_GFOFFSET3D(u,4,0,3)))),kmadd(ToReal(-9),kadd(KRANC_GFOFFSET3D(u,-4,0,4),KRANC_GFOFFSET3D(u,4,0,-4)),kmul(ToReal(9),kadd(KRANC_GFOFFSET3D(u,-4,0,-4),KRANC_GFOFFSET3D(u,4,0,4))))))))))))))))))))))))
#else
#  define PDstandardNthfdOrder831(u) (PDstandardNthfdOrder831_impl(u,p1o705600dxdz,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder831_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o705600dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder831_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o705600dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardNthfdOrder812_impl(u, p1o705600dxdz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder832(u) (kmul(p1o705600dydz,kmadd(ToReal(-451584),kadd(KRANC_GFOFFSET3D(u,0,-1,1),KRANC_GFOFFSET3D(u,0,1,-1)),kmadd(ToReal(451584),kadd(KRANC_GFOFFSET3D(u,0,-1,-1),KRANC_GFOFFSET3D(u,0,1,1)),kmadd(ToReal(112896),kadd(KRANC_GFOFFSET3D(u,0,-1,2),kadd(KRANC_GFOFFSET3D(u,0,1,-2),kadd(KRANC_GFOFFSET3D(u,0,-2,1),KRANC_GFOFFSET3D(u,0,2,-1)))),kmadd(ToReal(-112896),kadd(KRANC_GFOFFSET3D(u,0,-1,-2),kadd(KRANC_GFOFFSET3D(u,0,1,2),kadd(KRANC_GFOFFSET3D(u,0,-2,-1),KRANC_GFOFFSET3D(u,0,2,1)))),kmadd(ToReal(-28224),kadd(KRANC_GFOFFSET3D(u,0,-2,2),KRANC_GFOFFSET3D(u,0,2,-2)),kmadd(ToReal(28224),kadd(KRANC_GFOFFSET3D(u,0,-2,-2),KRANC_GFOFFSET3D(u,0,2,2)),kmadd(ToReal(-21504),kadd(KRANC_GFOFFSET3D(u,0,-1,3),kadd(KRANC_GFOFFSET3D(u,0,1,-3),kadd(KRANC_GFOFFSET3D(u,0,-3,1),KRANC_GFOFFSET3D(u,0,3,-1)))),kmadd(ToReal(21504),kadd(KRANC_GFOFFSET3D(u,0,-1,-3),kadd(KRANC_GFOFFSET3D(u,0,1,3),kadd(KRANC_GFOFFSET3D(u,0,-3,-1),KRANC_GFOFFSET3D(u,0,3,1)))),kmadd(ToReal(5376),kadd(KRANC_GFOFFSET3D(u,0,-2,3),kadd(KRANC_GFOFFSET3D(u,0,2,-3),kadd(KRANC_GFOFFSET3D(u,0,-3,2),KRANC_GFOFFSET3D(u,0,3,-2)))),kmadd(ToReal(-5376),kadd(KRANC_GFOFFSET3D(u,0,-2,-3),kadd(KRANC_GFOFFSET3D(u,0,2,3),kadd(KRANC_GFOFFSET3D(u,0,-3,-2),KRANC_GFOFFSET3D(u,0,3,2)))),kmadd(ToReal(-1024),kadd(KRANC_GFOFFSET3D(u,0,-3,3),KRANC_GFOFFSET3D(u,0,3,-3)),kmadd(ToReal(1024),kadd(KRANC_GFOFFSET3D(u,0,-3,-3),KRANC_GFOFFSET3D(u,0,3,3)),kmadd(ToReal(2016),kadd(KRANC_GFOFFSET3D(u,0,-1,4),kadd(KRANC_GFOFFSET3D(u,0,1,-4),kadd(KRANC_GFOFFSET3D(u,0,-4,1),KRANC_GFOFFSET3D(u,0,4,-1)))),kmadd(ToReal(-2016),kadd(KRANC_GFOFFSET3D(u,0,-1,-4),kadd(KRANC_GFOFFSET3D(u,0,1,4),kadd(KRANC_GFOFFSET3D(u,0,-4,-1),KRANC_GFOFFSET3D(u,0,4,1)))),kmadd(ToReal(-504),kadd(KRANC_GFOFFSET3D(u,0,-2,4),kadd(KRANC_GFOFFSET3D(u,0,2,-4),kadd(KRANC_GFOFFSET3D(u,0,-4,2),KRANC_GFOFFSET3D(u,0,4,-2)))),kmadd(ToReal(504),kadd(KRANC_GFOFFSET3D(u,0,-2,-4),kadd(KRANC_GFOFFSET3D(u,0,2,4),kadd(KRANC_GFOFFSET3D(u,0,-4,-2),KRANC_GFOFFSET3D(u,0,4,2)))),kmadd(ToReal(96),kadd(KRANC_GFOFFSET3D(u,0,-3,4),kadd(KRANC_GFOFFSET3D(u,0,3,-4),kadd(KRANC_GFOFFSET3D(u,0,-4,3),KRANC_GFOFFSET3D(u,0,4,-3)))),kmadd(ToReal(-96),kadd(KRANC_GFOFFSET3D(u,0,-3,-4),kadd(KRANC_GFOFFSET3D(u,0,3,4),kadd(KRANC_GFOFFSET3D(u,0,-4,-3),KRANC_GFOFFSET3D(u,0,4,3)))),kmadd(ToReal(-9),kadd(KRANC_GFOFFSET3D(u,0,-4,4),KRANC_GFOFFSET3D(u,0,4,-4)),kmul(ToReal(9),kadd(KRANC_GFOFFSET3D(u,0,-4,-4),KRANC_GFOFFSET3D(u,0,4,4))))))))))))))))))))))))
#else
#  define PDstandardNthfdOrder832(u) (PDstandardNthfdOrder832_impl(u,p1o705600dydz,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder832_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o705600dydz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder832_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o705600dydz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardNthfdOrder823_impl(u, p1o705600dydz, cdj, cdk);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthfdOrder21(u) (kmul(kmadd(ToReal(3),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(-4),KRANC_GFOFFSET3D(u,1,0,0),KRANC_GFOFFSET3D(u,2,0,0))),kmul(pm1o2dx,dir1)))
#else
#  define PDupwindNthfdOrder21(u) (PDupwindNthfdOrder21_impl(u,pm1o2dx,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthfdOrder21_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o2dx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthfdOrder21_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o2dx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{ assert(0); return ToReal(1e30); /* ERROR */ }
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthfdOrder22(u) (kmul(kmadd(ToReal(3),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(-4),KRANC_GFOFFSET3D(u,0,1,0),KRANC_GFOFFSET3D(u,0,2,0))),kmul(pm1o2dy,dir2)))
#else
#  define PDupwindNthfdOrder22(u) (PDupwindNthfdOrder22_impl(u,pm1o2dy,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthfdOrder22_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o2dy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthfdOrder22_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o2dy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{ assert(0); return ToReal(1e30); /* ERROR */ }
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthfdOrder23(u) (kmul(kmadd(ToReal(3),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(-4),KRANC_GFOFFSET3D(u,0,0,1),KRANC_GFOFFSET3D(u,0,0,2))),kmul(pm1o2dz,dir3)))
#else
#  define PDupwindNthfdOrder23(u) (PDupwindNthfdOrder23_impl(u,pm1o2dz,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthfdOrder23_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o2dz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthfdOrder23_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o2dz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{ assert(0); return ToReal(1e30); /* ERROR */ }
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthfdOrder41(u) (kmul(kmadd(ToReal(-10),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(-3),KRANC_GFOFFSET3D(u,-1,0,0),kmadd(ToReal(18),KRANC_GFOFFSET3D(u,1,0,0),kmadd(ToReal(-6),KRANC_GFOFFSET3D(u,2,0,0),KRANC_GFOFFSET3D(u,3,0,0))))),kmul(p1o12dx,dir1)))
#else
#  define PDupwindNthfdOrder41(u) (PDupwindNthfdOrder41_impl(u,p1o12dx,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthfdOrder41_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o12dx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthfdOrder41_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o12dx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{ assert(0); return ToReal(1e30); /* ERROR */ }
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthfdOrder42(u) (kmul(kmadd(ToReal(-10),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(-3),KRANC_GFOFFSET3D(u,0,-1,0),kmadd(ToReal(18),KRANC_GFOFFSET3D(u,0,1,0),kmadd(ToReal(-6),KRANC_GFOFFSET3D(u,0,2,0),KRANC_GFOFFSET3D(u,0,3,0))))),kmul(p1o12dy,dir2)))
#else
#  define PDupwindNthfdOrder42(u) (PDupwindNthfdOrder42_impl(u,p1o12dy,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthfdOrder42_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o12dy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthfdOrder42_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o12dy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{ assert(0); return ToReal(1e30); /* ERROR */ }
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthfdOrder43(u) (kmul(kmadd(ToReal(-10),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(-3),KRANC_GFOFFSET3D(u,0,0,-1),kmadd(ToReal(18),KRANC_GFOFFSET3D(u,0,0,1),kmadd(ToReal(-6),KRANC_GFOFFSET3D(u,0,0,2),KRANC_GFOFFSET3D(u,0,0,3))))),kmul(p1o12dz,dir3)))
#else
#  define PDupwindNthfdOrder43(u) (PDupwindNthfdOrder43_impl(u,p1o12dz,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthfdOrder43_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o12dz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthfdOrder43_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o12dz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{ assert(0); return ToReal(1e30); /* ERROR */ }
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthfdOrder61(u) (kmul(kmadd(ToReal(35),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(24),KRANC_GFOFFSET3D(u,-1,0,0),kmadd(ToReal(-80),KRANC_GFOFFSET3D(u,1,0,0),kmadd(ToReal(-2),KRANC_GFOFFSET3D(u,-2,0,0),kmadd(ToReal(30),KRANC_GFOFFSET3D(u,2,0,0),kmadd(ToReal(-8),KRANC_GFOFFSET3D(u,3,0,0),KRANC_GFOFFSET3D(u,4,0,0))))))),kmul(pm1o60dx,dir1)))
#else
#  define PDupwindNthfdOrder61(u) (PDupwindNthfdOrder61_impl(u,pm1o60dx,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthfdOrder61_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o60dx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthfdOrder61_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o60dx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{ assert(0); return ToReal(1e30); /* ERROR */ }
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthfdOrder62(u) (kmul(kmadd(ToReal(35),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(24),KRANC_GFOFFSET3D(u,0,-1,0),kmadd(ToReal(-80),KRANC_GFOFFSET3D(u,0,1,0),kmadd(ToReal(-2),KRANC_GFOFFSET3D(u,0,-2,0),kmadd(ToReal(30),KRANC_GFOFFSET3D(u,0,2,0),kmadd(ToReal(-8),KRANC_GFOFFSET3D(u,0,3,0),KRANC_GFOFFSET3D(u,0,4,0))))))),kmul(pm1o60dy,dir2)))
#else
#  define PDupwindNthfdOrder62(u) (PDupwindNthfdOrder62_impl(u,pm1o60dy,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthfdOrder62_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o60dy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthfdOrder62_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o60dy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{ assert(0); return ToReal(1e30); /* ERROR */ }
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthfdOrder63(u) (kmul(kmadd(ToReal(35),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(24),KRANC_GFOFFSET3D(u,0,0,-1),kmadd(ToReal(-80),KRANC_GFOFFSET3D(u,0,0,1),kmadd(ToReal(-2),KRANC_GFOFFSET3D(u,0,0,-2),kmadd(ToReal(30),KRANC_GFOFFSET3D(u,0,0,2),kmadd(ToReal(-8),KRANC_GFOFFSET3D(u,0,0,3),KRANC_GFOFFSET3D(u,0,0,4))))))),kmul(pm1o60dz,dir3)))
#else
#  define PDupwindNthfdOrder63(u) (PDupwindNthfdOrder63_impl(u,pm1o60dz,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthfdOrder63_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o60dz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthfdOrder63_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o60dz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{ assert(0); return ToReal(1e30); /* ERROR */ }
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthfdOrder81(u) (kmul(kmadd(ToReal(-378),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(-5),kmadd(ToReal(-210),KRANC_GFOFFSET3D(u,1,0,0),kmadd(ToReal(-12),KRANC_GFOFFSET3D(u,-2,0,0),kmadd(ToReal(84),kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,2,0,0)),kadd(KRANC_GFOFFSET3D(u,-3,0,0),kmadd(ToReal(-28),KRANC_GFOFFSET3D(u,3,0,0),kmul(KRANC_GFOFFSET3D(u,4,0,0),ToReal(6))))))),kmul(KRANC_GFOFFSET3D(u,5,0,0),ToReal(3)))),kmul(p1o840dx,dir1)))
#else
#  define PDupwindNthfdOrder81(u) (PDupwindNthfdOrder81_impl(u,p1o840dx,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthfdOrder81_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o840dx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthfdOrder81_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o840dx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{ assert(0); return ToReal(1e30); /* ERROR */ }
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthfdOrder82(u) (kmul(kmadd(ToReal(-378),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(-5),kmadd(ToReal(-210),KRANC_GFOFFSET3D(u,0,1,0),kmadd(ToReal(-12),KRANC_GFOFFSET3D(u,0,-2,0),kmadd(ToReal(84),kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,2,0)),kadd(KRANC_GFOFFSET3D(u,0,-3,0),kmadd(ToReal(-28),KRANC_GFOFFSET3D(u,0,3,0),kmul(KRANC_GFOFFSET3D(u,0,4,0),ToReal(6))))))),kmul(KRANC_GFOFFSET3D(u,0,5,0),ToReal(3)))),kmul(p1o840dy,dir2)))
#else
#  define PDupwindNthfdOrder82(u) (PDupwindNthfdOrder82_impl(u,p1o840dy,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthfdOrder82_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o840dy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthfdOrder82_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o840dy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{ assert(0); return ToReal(1e30); /* ERROR */ }
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthfdOrder83(u) (kmul(kmadd(ToReal(-378),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(-5),kmadd(ToReal(-210),KRANC_GFOFFSET3D(u,0,0,1),kmadd(ToReal(-12),KRANC_GFOFFSET3D(u,0,0,-2),kmadd(ToReal(84),kadd(KRANC_GFOFFSET3D(u,0,0,-1),KRANC_GFOFFSET3D(u,0,0,2)),kadd(KRANC_GFOFFSET3D(u,0,0,-3),kmadd(ToReal(-28),KRANC_GFOFFSET3D(u,0,0,3),kmul(KRANC_GFOFFSET3D(u,0,0,4),ToReal(6))))))),kmul(KRANC_GFOFFSET3D(u,0,0,5),ToReal(3)))),kmul(p1o840dz,dir3)))
#else
#  define PDupwindNthfdOrder83(u) (PDupwindNthfdOrder83_impl(u,p1o840dz,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthfdOrder83_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o840dz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthfdOrder83_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o840dz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{ assert(0); return ToReal(1e30); /* ERROR */ }
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthSymmfdOrder21(u) (kmul(pm1o4dx,kmadd(ToReal(6),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(-4),kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),kadd(KRANC_GFOFFSET3D(u,-2,0,0),KRANC_GFOFFSET3D(u,2,0,0))))))
#else
#  define PDupwindNthSymmfdOrder21(u) (PDupwindNthSymmfdOrder21_impl(u,pm1o4dx,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthSymmfdOrder21_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o4dx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthSymmfdOrder21_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o4dx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(pm1o4dx,kmadd(ToReal(6),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(-4),kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),kadd(KRANC_GFOFFSET3D(u,-2,0,0),KRANC_GFOFFSET3D(u,2,0,0)))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthSymmfdOrder22(u) (kmul(pm1o4dy,kmadd(ToReal(6),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(-4),kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),kadd(KRANC_GFOFFSET3D(u,0,-2,0),KRANC_GFOFFSET3D(u,0,2,0))))))
#else
#  define PDupwindNthSymmfdOrder22(u) (PDupwindNthSymmfdOrder22_impl(u,pm1o4dy,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthSymmfdOrder22_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o4dy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthSymmfdOrder22_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o4dy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(pm1o4dy,kmadd(ToReal(6),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(-4),kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),kadd(KRANC_GFOFFSET3D(u,0,-2,0),KRANC_GFOFFSET3D(u,0,2,0)))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthSymmfdOrder23(u) (kmul(pm1o4dz,kmadd(ToReal(6),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(-4),kadd(KRANC_GFOFFSET3D(u,0,0,-1),KRANC_GFOFFSET3D(u,0,0,1)),kadd(KRANC_GFOFFSET3D(u,0,0,-2),KRANC_GFOFFSET3D(u,0,0,2))))))
#else
#  define PDupwindNthSymmfdOrder23(u) (PDupwindNthSymmfdOrder23_impl(u,pm1o4dz,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthSymmfdOrder23_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o4dz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthSymmfdOrder23_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o4dz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDupwindNthSymmfdOrder22_impl(u, pm1o4dz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthSymmfdOrder41(u) (kmul(p1o24dx,kmadd(ToReal(-20),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(15),kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),kmadd(ToReal(-6),kadd(KRANC_GFOFFSET3D(u,-2,0,0),KRANC_GFOFFSET3D(u,2,0,0)),kadd(KRANC_GFOFFSET3D(u,-3,0,0),KRANC_GFOFFSET3D(u,3,0,0)))))))
#else
#  define PDupwindNthSymmfdOrder41(u) (PDupwindNthSymmfdOrder41_impl(u,p1o24dx,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthSymmfdOrder41_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o24dx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthSymmfdOrder41_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o24dx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o24dx,kmadd(ToReal(-20),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(15),kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),kmadd(ToReal(-6),kadd(KRANC_GFOFFSET3D(u,-2,0,0),KRANC_GFOFFSET3D(u,2,0,0)),kadd(KRANC_GFOFFSET3D(u,-3,0,0),KRANC_GFOFFSET3D(u,3,0,0))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthSymmfdOrder42(u) (kmul(p1o24dy,kmadd(ToReal(-20),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(15),kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),kmadd(ToReal(-6),kadd(KRANC_GFOFFSET3D(u,0,-2,0),KRANC_GFOFFSET3D(u,0,2,0)),kadd(KRANC_GFOFFSET3D(u,0,-3,0),KRANC_GFOFFSET3D(u,0,3,0)))))))
#else
#  define PDupwindNthSymmfdOrder42(u) (PDupwindNthSymmfdOrder42_impl(u,p1o24dy,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthSymmfdOrder42_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o24dy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthSymmfdOrder42_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o24dy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o24dy,kmadd(ToReal(-20),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(15),kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),kmadd(ToReal(-6),kadd(KRANC_GFOFFSET3D(u,0,-2,0),KRANC_GFOFFSET3D(u,0,2,0)),kadd(KRANC_GFOFFSET3D(u,0,-3,0),KRANC_GFOFFSET3D(u,0,3,0))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthSymmfdOrder43(u) (kmul(p1o24dz,kmadd(ToReal(-20),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(15),kadd(KRANC_GFOFFSET3D(u,0,0,-1),KRANC_GFOFFSET3D(u,0,0,1)),kmadd(ToReal(-6),kadd(KRANC_GFOFFSET3D(u,0,0,-2),KRANC_GFOFFSET3D(u,0,0,2)),kadd(KRANC_GFOFFSET3D(u,0,0,-3),KRANC_GFOFFSET3D(u,0,0,3)))))))
#else
#  define PDupwindNthSymmfdOrder43(u) (PDupwindNthSymmfdOrder43_impl(u,p1o24dz,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthSymmfdOrder43_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o24dz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthSymmfdOrder43_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o24dz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDupwindNthSymmfdOrder42_impl(u, p1o24dz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthSymmfdOrder61(u) (kmul(pm1o120dx,kmadd(ToReal(70),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(-56),kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),kmadd(ToReal(28),kadd(KRANC_GFOFFSET3D(u,-2,0,0),KRANC_GFOFFSET3D(u,2,0,0)),kmadd(ToReal(-8),kadd(KRANC_GFOFFSET3D(u,-3,0,0),KRANC_GFOFFSET3D(u,3,0,0)),kadd(KRANC_GFOFFSET3D(u,-4,0,0),KRANC_GFOFFSET3D(u,4,0,0))))))))
#else
#  define PDupwindNthSymmfdOrder61(u) (PDupwindNthSymmfdOrder61_impl(u,pm1o120dx,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthSymmfdOrder61_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o120dx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthSymmfdOrder61_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o120dx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(pm1o120dx,kmadd(ToReal(70),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(-56),kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),kmadd(ToReal(28),kadd(KRANC_GFOFFSET3D(u,-2,0,0),KRANC_GFOFFSET3D(u,2,0,0)),kmadd(ToReal(-8),kadd(KRANC_GFOFFSET3D(u,-3,0,0),KRANC_GFOFFSET3D(u,3,0,0)),kadd(KRANC_GFOFFSET3D(u,-4,0,0),KRANC_GFOFFSET3D(u,4,0,0)))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthSymmfdOrder62(u) (kmul(pm1o120dy,kmadd(ToReal(70),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(-56),kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),kmadd(ToReal(28),kadd(KRANC_GFOFFSET3D(u,0,-2,0),KRANC_GFOFFSET3D(u,0,2,0)),kmadd(ToReal(-8),kadd(KRANC_GFOFFSET3D(u,0,-3,0),KRANC_GFOFFSET3D(u,0,3,0)),kadd(KRANC_GFOFFSET3D(u,0,-4,0),KRANC_GFOFFSET3D(u,0,4,0))))))))
#else
#  define PDupwindNthSymmfdOrder62(u) (PDupwindNthSymmfdOrder62_impl(u,pm1o120dy,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthSymmfdOrder62_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o120dy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthSymmfdOrder62_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o120dy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(pm1o120dy,kmadd(ToReal(70),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(-56),kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),kmadd(ToReal(28),kadd(KRANC_GFOFFSET3D(u,0,-2,0),KRANC_GFOFFSET3D(u,0,2,0)),kmadd(ToReal(-8),kadd(KRANC_GFOFFSET3D(u,0,-3,0),KRANC_GFOFFSET3D(u,0,3,0)),kadd(KRANC_GFOFFSET3D(u,0,-4,0),KRANC_GFOFFSET3D(u,0,4,0)))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthSymmfdOrder63(u) (kmul(pm1o120dz,kmadd(ToReal(70),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(-56),kadd(KRANC_GFOFFSET3D(u,0,0,-1),KRANC_GFOFFSET3D(u,0,0,1)),kmadd(ToReal(28),kadd(KRANC_GFOFFSET3D(u,0,0,-2),KRANC_GFOFFSET3D(u,0,0,2)),kmadd(ToReal(-8),kadd(KRANC_GFOFFSET3D(u,0,0,-3),KRANC_GFOFFSET3D(u,0,0,3)),kadd(KRANC_GFOFFSET3D(u,0,0,-4),KRANC_GFOFFSET3D(u,0,0,4))))))))
#else
#  define PDupwindNthSymmfdOrder63(u) (PDupwindNthSymmfdOrder63_impl(u,pm1o120dz,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthSymmfdOrder63_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o120dz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthSymmfdOrder63_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o120dz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDupwindNthSymmfdOrder62_impl(u, pm1o120dz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthSymmfdOrder81(u) (kmul(p1o560dx,kmadd(ToReal(-252),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(210),kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),kmadd(ToReal(-120),kadd(KRANC_GFOFFSET3D(u,-2,0,0),KRANC_GFOFFSET3D(u,2,0,0)),kmadd(ToReal(45),kadd(KRANC_GFOFFSET3D(u,-3,0,0),KRANC_GFOFFSET3D(u,3,0,0)),kmadd(ToReal(-10),kadd(KRANC_GFOFFSET3D(u,-4,0,0),KRANC_GFOFFSET3D(u,4,0,0)),kadd(KRANC_GFOFFSET3D(u,-5,0,0),KRANC_GFOFFSET3D(u,5,0,0)))))))))
#else
#  define PDupwindNthSymmfdOrder81(u) (PDupwindNthSymmfdOrder81_impl(u,p1o560dx,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthSymmfdOrder81_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o560dx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthSymmfdOrder81_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o560dx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o560dx,kmadd(ToReal(-252),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(210),kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),kmadd(ToReal(-120),kadd(KRANC_GFOFFSET3D(u,-2,0,0),KRANC_GFOFFSET3D(u,2,0,0)),kmadd(ToReal(45),kadd(KRANC_GFOFFSET3D(u,-3,0,0),KRANC_GFOFFSET3D(u,3,0,0)),kmadd(ToReal(-10),kadd(KRANC_GFOFFSET3D(u,-4,0,0),KRANC_GFOFFSET3D(u,4,0,0)),kadd(KRANC_GFOFFSET3D(u,-5,0,0),KRANC_GFOFFSET3D(u,5,0,0))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthSymmfdOrder82(u) (kmul(p1o560dy,kmadd(ToReal(-252),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(210),kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),kmadd(ToReal(-120),kadd(KRANC_GFOFFSET3D(u,0,-2,0),KRANC_GFOFFSET3D(u,0,2,0)),kmadd(ToReal(45),kadd(KRANC_GFOFFSET3D(u,0,-3,0),KRANC_GFOFFSET3D(u,0,3,0)),kmadd(ToReal(-10),kadd(KRANC_GFOFFSET3D(u,0,-4,0),KRANC_GFOFFSET3D(u,0,4,0)),kadd(KRANC_GFOFFSET3D(u,0,-5,0),KRANC_GFOFFSET3D(u,0,5,0)))))))))
#else
#  define PDupwindNthSymmfdOrder82(u) (PDupwindNthSymmfdOrder82_impl(u,p1o560dy,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthSymmfdOrder82_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o560dy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthSymmfdOrder82_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o560dy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o560dy,kmadd(ToReal(-252),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(210),kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),kmadd(ToReal(-120),kadd(KRANC_GFOFFSET3D(u,0,-2,0),KRANC_GFOFFSET3D(u,0,2,0)),kmadd(ToReal(45),kadd(KRANC_GFOFFSET3D(u,0,-3,0),KRANC_GFOFFSET3D(u,0,3,0)),kmadd(ToReal(-10),kadd(KRANC_GFOFFSET3D(u,0,-4,0),KRANC_GFOFFSET3D(u,0,4,0)),kadd(KRANC_GFOFFSET3D(u,0,-5,0),KRANC_GFOFFSET3D(u,0,5,0))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthSymmfdOrder83(u) (kmul(p1o560dz,kmadd(ToReal(-252),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(210),kadd(KRANC_GFOFFSET3D(u,0,0,-1),KRANC_GFOFFSET3D(u,0,0,1)),kmadd(ToReal(-120),kadd(KRANC_GFOFFSET3D(u,0,0,-2),KRANC_GFOFFSET3D(u,0,0,2)),kmadd(ToReal(45),kadd(KRANC_GFOFFSET3D(u,0,0,-3),KRANC_GFOFFSET3D(u,0,0,3)),kmadd(ToReal(-10),kadd(KRANC_GFOFFSET3D(u,0,0,-4),KRANC_GFOFFSET3D(u,0,0,4)),kadd(KRANC_GFOFFSET3D(u,0,0,-5),KRANC_GFOFFSET3D(u,0,0,5)))))))))
#else
#  define PDupwindNthSymmfdOrder83(u) (PDupwindNthSymmfdOrder83_impl(u,p1o560dz,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthSymmfdOrder83_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o560dz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthSymmfdOrder83_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o560dz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDupwindNthSymmfdOrder82_impl(u, p1o560dz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthAntifdOrder21(u) (kmul(p1o4dx,kmadd(ToReal(-4),KRANC_GFOFFSET3D(u,-1,0,0),kmadd(ToReal(4),KRANC_GFOFFSET3D(u,1,0,0),ksub(KRANC_GFOFFSET3D(u,-2,0,0),KRANC_GFOFFSET3D(u,2,0,0))))))
#else
#  define PDupwindNthAntifdOrder21(u) (PDupwindNthAntifdOrder21_impl(u,p1o4dx,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthAntifdOrder21_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o4dx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthAntifdOrder21_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o4dx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o4dx,kmadd(ToReal(-4),KRANC_GFOFFSET3D(u,-1,0,0),kmadd(ToReal(4),KRANC_GFOFFSET3D(u,1,0,0),ksub(KRANC_GFOFFSET3D(u,-2,0,0),KRANC_GFOFFSET3D(u,2,0,0)))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthAntifdOrder22(u) (kmul(p1o4dy,kmadd(ToReal(-4),KRANC_GFOFFSET3D(u,0,-1,0),kmadd(ToReal(4),KRANC_GFOFFSET3D(u,0,1,0),ksub(KRANC_GFOFFSET3D(u,0,-2,0),KRANC_GFOFFSET3D(u,0,2,0))))))
#else
#  define PDupwindNthAntifdOrder22(u) (PDupwindNthAntifdOrder22_impl(u,p1o4dy,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthAntifdOrder22_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o4dy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthAntifdOrder22_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o4dy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o4dy,kmadd(ToReal(-4),KRANC_GFOFFSET3D(u,0,-1,0),kmadd(ToReal(4),KRANC_GFOFFSET3D(u,0,1,0),ksub(KRANC_GFOFFSET3D(u,0,-2,0),KRANC_GFOFFSET3D(u,0,2,0)))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthAntifdOrder23(u) (kmul(p1o4dz,kmadd(ToReal(-4),KRANC_GFOFFSET3D(u,0,0,-1),kmadd(ToReal(4),KRANC_GFOFFSET3D(u,0,0,1),ksub(KRANC_GFOFFSET3D(u,0,0,-2),KRANC_GFOFFSET3D(u,0,0,2))))))
#else
#  define PDupwindNthAntifdOrder23(u) (PDupwindNthAntifdOrder23_impl(u,p1o4dz,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthAntifdOrder23_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o4dz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthAntifdOrder23_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o4dz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDupwindNthAntifdOrder22_impl(u, p1o4dz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthAntifdOrder41(u) (kmul(p1o24dx,kmadd(ToReal(-21),KRANC_GFOFFSET3D(u,-1,0,0),kmadd(ToReal(21),KRANC_GFOFFSET3D(u,1,0,0),kmadd(ToReal(6),KRANC_GFOFFSET3D(u,-2,0,0),kmadd(ToReal(-6),KRANC_GFOFFSET3D(u,2,0,0),ksub(KRANC_GFOFFSET3D(u,3,0,0),KRANC_GFOFFSET3D(u,-3,0,0))))))))
#else
#  define PDupwindNthAntifdOrder41(u) (PDupwindNthAntifdOrder41_impl(u,p1o24dx,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthAntifdOrder41_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o24dx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthAntifdOrder41_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o24dx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o24dx,kmadd(ToReal(-21),KRANC_GFOFFSET3D(u,-1,0,0),kmadd(ToReal(21),KRANC_GFOFFSET3D(u,1,0,0),kmadd(ToReal(6),KRANC_GFOFFSET3D(u,-2,0,0),kmadd(ToReal(-6),KRANC_GFOFFSET3D(u,2,0,0),ksub(KRANC_GFOFFSET3D(u,3,0,0),KRANC_GFOFFSET3D(u,-3,0,0)))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthAntifdOrder42(u) (kmul(p1o24dy,kmadd(ToReal(-21),KRANC_GFOFFSET3D(u,0,-1,0),kmadd(ToReal(21),KRANC_GFOFFSET3D(u,0,1,0),kmadd(ToReal(6),KRANC_GFOFFSET3D(u,0,-2,0),kmadd(ToReal(-6),KRANC_GFOFFSET3D(u,0,2,0),ksub(KRANC_GFOFFSET3D(u,0,3,0),KRANC_GFOFFSET3D(u,0,-3,0))))))))
#else
#  define PDupwindNthAntifdOrder42(u) (PDupwindNthAntifdOrder42_impl(u,p1o24dy,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthAntifdOrder42_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o24dy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthAntifdOrder42_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o24dy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o24dy,kmadd(ToReal(-21),KRANC_GFOFFSET3D(u,0,-1,0),kmadd(ToReal(21),KRANC_GFOFFSET3D(u,0,1,0),kmadd(ToReal(6),KRANC_GFOFFSET3D(u,0,-2,0),kmadd(ToReal(-6),KRANC_GFOFFSET3D(u,0,2,0),ksub(KRANC_GFOFFSET3D(u,0,3,0),KRANC_GFOFFSET3D(u,0,-3,0)))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthAntifdOrder43(u) (kmul(p1o24dz,kmadd(ToReal(-21),KRANC_GFOFFSET3D(u,0,0,-1),kmadd(ToReal(21),KRANC_GFOFFSET3D(u,0,0,1),kmadd(ToReal(6),KRANC_GFOFFSET3D(u,0,0,-2),kmadd(ToReal(-6),KRANC_GFOFFSET3D(u,0,0,2),ksub(KRANC_GFOFFSET3D(u,0,0,3),KRANC_GFOFFSET3D(u,0,0,-3))))))))
#else
#  define PDupwindNthAntifdOrder43(u) (PDupwindNthAntifdOrder43_impl(u,p1o24dz,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthAntifdOrder43_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o24dz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthAntifdOrder43_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o24dz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDupwindNthAntifdOrder42_impl(u, p1o24dz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthAntifdOrder61(u) (kmul(p1o120dx,kmadd(ToReal(-104),KRANC_GFOFFSET3D(u,-1,0,0),kmadd(ToReal(104),KRANC_GFOFFSET3D(u,1,0,0),kmadd(ToReal(32),KRANC_GFOFFSET3D(u,-2,0,0),kmadd(ToReal(-32),KRANC_GFOFFSET3D(u,2,0,0),kmadd(ToReal(-8),KRANC_GFOFFSET3D(u,-3,0,0),kmadd(ToReal(8),KRANC_GFOFFSET3D(u,3,0,0),ksub(KRANC_GFOFFSET3D(u,-4,0,0),KRANC_GFOFFSET3D(u,4,0,0))))))))))
#else
#  define PDupwindNthAntifdOrder61(u) (PDupwindNthAntifdOrder61_impl(u,p1o120dx,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthAntifdOrder61_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o120dx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthAntifdOrder61_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o120dx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o120dx,kmadd(ToReal(-104),KRANC_GFOFFSET3D(u,-1,0,0),kmadd(ToReal(104),KRANC_GFOFFSET3D(u,1,0,0),kmadd(ToReal(32),KRANC_GFOFFSET3D(u,-2,0,0),kmadd(ToReal(-32),KRANC_GFOFFSET3D(u,2,0,0),kmadd(ToReal(-8),KRANC_GFOFFSET3D(u,-3,0,0),kmadd(ToReal(8),KRANC_GFOFFSET3D(u,3,0,0),ksub(KRANC_GFOFFSET3D(u,-4,0,0),KRANC_GFOFFSET3D(u,4,0,0)))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthAntifdOrder62(u) (kmul(p1o120dy,kmadd(ToReal(-104),KRANC_GFOFFSET3D(u,0,-1,0),kmadd(ToReal(104),KRANC_GFOFFSET3D(u,0,1,0),kmadd(ToReal(32),KRANC_GFOFFSET3D(u,0,-2,0),kmadd(ToReal(-32),KRANC_GFOFFSET3D(u,0,2,0),kmadd(ToReal(-8),KRANC_GFOFFSET3D(u,0,-3,0),kmadd(ToReal(8),KRANC_GFOFFSET3D(u,0,3,0),ksub(KRANC_GFOFFSET3D(u,0,-4,0),KRANC_GFOFFSET3D(u,0,4,0))))))))))
#else
#  define PDupwindNthAntifdOrder62(u) (PDupwindNthAntifdOrder62_impl(u,p1o120dy,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthAntifdOrder62_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o120dy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthAntifdOrder62_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o120dy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o120dy,kmadd(ToReal(-104),KRANC_GFOFFSET3D(u,0,-1,0),kmadd(ToReal(104),KRANC_GFOFFSET3D(u,0,1,0),kmadd(ToReal(32),KRANC_GFOFFSET3D(u,0,-2,0),kmadd(ToReal(-32),KRANC_GFOFFSET3D(u,0,2,0),kmadd(ToReal(-8),KRANC_GFOFFSET3D(u,0,-3,0),kmadd(ToReal(8),KRANC_GFOFFSET3D(u,0,3,0),ksub(KRANC_GFOFFSET3D(u,0,-4,0),KRANC_GFOFFSET3D(u,0,4,0)))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthAntifdOrder63(u) (kmul(p1o120dz,kmadd(ToReal(-104),KRANC_GFOFFSET3D(u,0,0,-1),kmadd(ToReal(104),KRANC_GFOFFSET3D(u,0,0,1),kmadd(ToReal(32),KRANC_GFOFFSET3D(u,0,0,-2),kmadd(ToReal(-32),KRANC_GFOFFSET3D(u,0,0,2),kmadd(ToReal(-8),KRANC_GFOFFSET3D(u,0,0,-3),kmadd(ToReal(8),KRANC_GFOFFSET3D(u,0,0,3),ksub(KRANC_GFOFFSET3D(u,0,0,-4),KRANC_GFOFFSET3D(u,0,0,4))))))))))
#else
#  define PDupwindNthAntifdOrder63(u) (PDupwindNthAntifdOrder63_impl(u,p1o120dz,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthAntifdOrder63_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o120dz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthAntifdOrder63_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o120dz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDupwindNthAntifdOrder62_impl(u, p1o120dz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthAntifdOrder81(u) (kmul(p1o1680dx,kmadd(ToReal(-1470),KRANC_GFOFFSET3D(u,-1,0,0),kmadd(ToReal(1470),KRANC_GFOFFSET3D(u,1,0,0),kmadd(ToReal(480),KRANC_GFOFFSET3D(u,-2,0,0),kmadd(ToReal(-480),KRANC_GFOFFSET3D(u,2,0,0),kmadd(ToReal(-145),KRANC_GFOFFSET3D(u,-3,0,0),kmadd(ToReal(145),KRANC_GFOFFSET3D(u,3,0,0),kmadd(ToReal(30),KRANC_GFOFFSET3D(u,-4,0,0),kmadd(ToReal(-30),KRANC_GFOFFSET3D(u,4,0,0),kmadd(ToReal(-3),KRANC_GFOFFSET3D(u,-5,0,0),kmul(ToReal(3),KRANC_GFOFFSET3D(u,5,0,0)))))))))))))
#else
#  define PDupwindNthAntifdOrder81(u) (PDupwindNthAntifdOrder81_impl(u,p1o1680dx,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthAntifdOrder81_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o1680dx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthAntifdOrder81_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o1680dx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o1680dx,kmadd(ToReal(-1470),KRANC_GFOFFSET3D(u,-1,0,0),kmadd(ToReal(1470),KRANC_GFOFFSET3D(u,1,0,0),kmadd(ToReal(480),KRANC_GFOFFSET3D(u,-2,0,0),kmadd(ToReal(-480),KRANC_GFOFFSET3D(u,2,0,0),kmadd(ToReal(-145),KRANC_GFOFFSET3D(u,-3,0,0),kmadd(ToReal(145),KRANC_GFOFFSET3D(u,3,0,0),kmadd(ToReal(30),KRANC_GFOFFSET3D(u,-4,0,0),kmadd(ToReal(-30),KRANC_GFOFFSET3D(u,4,0,0),kmadd(ToReal(-3),KRANC_GFOFFSET3D(u,-5,0,0),kmul(ToReal(3),KRANC_GFOFFSET3D(u,5,0,0))))))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthAntifdOrder82(u) (kmul(p1o1680dy,kmadd(ToReal(-1470),KRANC_GFOFFSET3D(u,0,-1,0),kmadd(ToReal(1470),KRANC_GFOFFSET3D(u,0,1,0),kmadd(ToReal(480),KRANC_GFOFFSET3D(u,0,-2,0),kmadd(ToReal(-480),KRANC_GFOFFSET3D(u,0,2,0),kmadd(ToReal(-145),KRANC_GFOFFSET3D(u,0,-3,0),kmadd(ToReal(145),KRANC_GFOFFSET3D(u,0,3,0),kmadd(ToReal(30),KRANC_GFOFFSET3D(u,0,-4,0),kmadd(ToReal(-30),KRANC_GFOFFSET3D(u,0,4,0),kmadd(ToReal(-3),KRANC_GFOFFSET3D(u,0,-5,0),kmul(ToReal(3),KRANC_GFOFFSET3D(u,0,5,0)))))))))))))
#else
#  define PDupwindNthAntifdOrder82(u) (PDupwindNthAntifdOrder82_impl(u,p1o1680dy,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthAntifdOrder82_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o1680dy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthAntifdOrder82_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o1680dy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o1680dy,kmadd(ToReal(-1470),KRANC_GFOFFSET3D(u,0,-1,0),kmadd(ToReal(1470),KRANC_GFOFFSET3D(u,0,1,0),kmadd(ToReal(480),KRANC_GFOFFSET3D(u,0,-2,0),kmadd(ToReal(-480),KRANC_GFOFFSET3D(u,0,2,0),kmadd(ToReal(-145),KRANC_GFOFFSET3D(u,0,-3,0),kmadd(ToReal(145),KRANC_GFOFFSET3D(u,0,3,0),kmadd(ToReal(30),KRANC_GFOFFSET3D(u,0,-4,0),kmadd(ToReal(-30),KRANC_GFOFFSET3D(u,0,4,0),kmadd(ToReal(-3),KRANC_GFOFFSET3D(u,0,-5,0),kmul(ToReal(3),KRANC_GFOFFSET3D(u,0,5,0))))))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthAntifdOrder83(u) (kmul(p1o1680dz,kmadd(ToReal(-1470),KRANC_GFOFFSET3D(u,0,0,-1),kmadd(ToReal(1470),KRANC_GFOFFSET3D(u,0,0,1),kmadd(ToReal(480),KRANC_GFOFFSET3D(u,0,0,-2),kmadd(ToReal(-480),KRANC_GFOFFSET3D(u,0,0,2),kmadd(ToReal(-145),KRANC_GFOFFSET3D(u,0,0,-3),kmadd(ToReal(145),KRANC_GFOFFSET3D(u,0,0,3),kmadd(ToReal(30),KRANC_GFOFFSET3D(u,0,0,-4),kmadd(ToReal(-30),KRANC_GFOFFSET3D(u,0,0,4),kmadd(ToReal(-3),KRANC_GFOFFSET3D(u,0,0,-5),kmul(ToReal(3),KRANC_GFOFFSET3D(u,0,0,5)))))))))))))
#else
#  define PDupwindNthAntifdOrder83(u) (PDupwindNthAntifdOrder83_impl(u,p1o1680dz,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthAntifdOrder83_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o1680dz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthAntifdOrder83_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o1680dz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDupwindNthAntifdOrder82_impl(u, p1o1680dz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDonesidedfdOrder21(u) (kmul(kmadd(ToReal(3),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(-4),KRANC_GFOFFSET3D(u,1,0,0),KRANC_GFOFFSET3D(u,2,0,0))),kmul(pm1o2dx,dir1)))
#else
#  define PDonesidedfdOrder21(u) (PDonesidedfdOrder21_impl(u,pm1o2dx,cdj,cdk))
static CCTK_REAL_VEC PDonesidedfdOrder21_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o2dx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDonesidedfdOrder21_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o2dx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{ assert(0); return ToReal(1e30); /* ERROR */ }
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDonesidedfdOrder22(u) (kmul(kmadd(ToReal(3),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(-4),KRANC_GFOFFSET3D(u,0,1,0),KRANC_GFOFFSET3D(u,0,2,0))),kmul(pm1o2dy,dir2)))
#else
#  define PDonesidedfdOrder22(u) (PDonesidedfdOrder22_impl(u,pm1o2dy,cdj,cdk))
static CCTK_REAL_VEC PDonesidedfdOrder22_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o2dy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDonesidedfdOrder22_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o2dy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{ assert(0); return ToReal(1e30); /* ERROR */ }
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDonesidedfdOrder23(u) (kmul(kmadd(ToReal(3),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(-4),KRANC_GFOFFSET3D(u,0,0,1),KRANC_GFOFFSET3D(u,0,0,2))),kmul(pm1o2dz,dir3)))
#else
#  define PDonesidedfdOrder23(u) (PDonesidedfdOrder23_impl(u,pm1o2dz,cdj,cdk))
static CCTK_REAL_VEC PDonesidedfdOrder23_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o2dz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDonesidedfdOrder23_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o2dz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{ assert(0); return ToReal(1e30); /* ERROR */ }
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDonesidedfdOrder41(u) (kmul(kmadd(ToReal(-11),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(18),KRANC_GFOFFSET3D(u,1,0,0),kmadd(ToReal(-9),KRANC_GFOFFSET3D(u,2,0,0),kmul(KRANC_GFOFFSET3D(u,3,0,0),ToReal(2))))),kmul(p1o6dx,dir1)))
#else
#  define PDonesidedfdOrder41(u) (PDonesidedfdOrder41_impl(u,p1o6dx,cdj,cdk))
static CCTK_REAL_VEC PDonesidedfdOrder41_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o6dx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDonesidedfdOrder41_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o6dx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{ assert(0); return ToReal(1e30); /* ERROR */ }
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDonesidedfdOrder42(u) (kmul(kmadd(ToReal(-11),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(18),KRANC_GFOFFSET3D(u,0,1,0),kmadd(ToReal(-9),KRANC_GFOFFSET3D(u,0,2,0),kmul(KRANC_GFOFFSET3D(u,0,3,0),ToReal(2))))),kmul(p1o6dy,dir2)))
#else
#  define PDonesidedfdOrder42(u) (PDonesidedfdOrder42_impl(u,p1o6dy,cdj,cdk))
static CCTK_REAL_VEC PDonesidedfdOrder42_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o6dy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDonesidedfdOrder42_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o6dy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{ assert(0); return ToReal(1e30); /* ERROR */ }
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDonesidedfdOrder43(u) (kmul(kmadd(ToReal(-11),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(18),KRANC_GFOFFSET3D(u,0,0,1),kmadd(ToReal(-9),KRANC_GFOFFSET3D(u,0,0,2),kmul(KRANC_GFOFFSET3D(u,0,0,3),ToReal(2))))),kmul(p1o6dz,dir3)))
#else
#  define PDonesidedfdOrder43(u) (PDonesidedfdOrder43_impl(u,p1o6dz,cdj,cdk))
static CCTK_REAL_VEC PDonesidedfdOrder43_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o6dz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDonesidedfdOrder43_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o6dz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{ assert(0); return ToReal(1e30); /* ERROR */ }
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDonesidedfdOrder61(u) (kmul(kmadd(ToReal(25),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(-48),KRANC_GFOFFSET3D(u,1,0,0),kmadd(ToReal(36),KRANC_GFOFFSET3D(u,2,0,0),kmadd(ToReal(-16),KRANC_GFOFFSET3D(u,3,0,0),kmul(KRANC_GFOFFSET3D(u,4,0,0),ToReal(3)))))),kmul(pm1o12dx,dir1)))
#else
#  define PDonesidedfdOrder61(u) (PDonesidedfdOrder61_impl(u,pm1o12dx,cdj,cdk))
static CCTK_REAL_VEC PDonesidedfdOrder61_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o12dx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDonesidedfdOrder61_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o12dx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{ assert(0); return ToReal(1e30); /* ERROR */ }
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDonesidedfdOrder62(u) (kmul(kmadd(ToReal(25),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(-48),KRANC_GFOFFSET3D(u,0,1,0),kmadd(ToReal(36),KRANC_GFOFFSET3D(u,0,2,0),kmadd(ToReal(-16),KRANC_GFOFFSET3D(u,0,3,0),kmul(KRANC_GFOFFSET3D(u,0,4,0),ToReal(3)))))),kmul(pm1o12dy,dir2)))
#else
#  define PDonesidedfdOrder62(u) (PDonesidedfdOrder62_impl(u,pm1o12dy,cdj,cdk))
static CCTK_REAL_VEC PDonesidedfdOrder62_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o12dy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDonesidedfdOrder62_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o12dy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{ assert(0); return ToReal(1e30); /* ERROR */ }
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDonesidedfdOrder63(u) (kmul(kmadd(ToReal(25),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(-48),KRANC_GFOFFSET3D(u,0,0,1),kmadd(ToReal(36),KRANC_GFOFFSET3D(u,0,0,2),kmadd(ToReal(-16),KRANC_GFOFFSET3D(u,0,0,3),kmul(KRANC_GFOFFSET3D(u,0,0,4),ToReal(3)))))),kmul(pm1o12dz,dir3)))
#else
#  define PDonesidedfdOrder63(u) (PDonesidedfdOrder63_impl(u,pm1o12dz,cdj,cdk))
static CCTK_REAL_VEC PDonesidedfdOrder63_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o12dz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDonesidedfdOrder63_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o12dz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{ assert(0); return ToReal(1e30); /* ERROR */ }
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDonesidedfdOrder81(u) (kmul(kmadd(ToReal(-137),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(25),kmadd(ToReal(12),KRANC_GFOFFSET3D(u,1,0,0),kmadd(ToReal(-12),KRANC_GFOFFSET3D(u,2,0,0),kmadd(ToReal(8),KRANC_GFOFFSET3D(u,3,0,0),kmul(KRANC_GFOFFSET3D(u,4,0,0),ToReal(-3))))),kmul(KRANC_GFOFFSET3D(u,5,0,0),ToReal(12)))),kmul(p1o60dx,dir1)))
#else
#  define PDonesidedfdOrder81(u) (PDonesidedfdOrder81_impl(u,p1o60dx,cdj,cdk))
static CCTK_REAL_VEC PDonesidedfdOrder81_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o60dx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDonesidedfdOrder81_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o60dx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{ assert(0); return ToReal(1e30); /* ERROR */ }
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDonesidedfdOrder82(u) (kmul(kmadd(ToReal(-137),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(25),kmadd(ToReal(12),KRANC_GFOFFSET3D(u,0,1,0),kmadd(ToReal(-12),KRANC_GFOFFSET3D(u,0,2,0),kmadd(ToReal(8),KRANC_GFOFFSET3D(u,0,3,0),kmul(KRANC_GFOFFSET3D(u,0,4,0),ToReal(-3))))),kmul(KRANC_GFOFFSET3D(u,0,5,0),ToReal(12)))),kmul(p1o60dy,dir2)))
#else
#  define PDonesidedfdOrder82(u) (PDonesidedfdOrder82_impl(u,p1o60dy,cdj,cdk))
static CCTK_REAL_VEC PDonesidedfdOrder82_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o60dy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDonesidedfdOrder82_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o60dy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{ assert(0); return ToReal(1e30); /* ERROR */ }
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDonesidedfdOrder83(u) (kmul(kmadd(ToReal(-137),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(25),kmadd(ToReal(12),KRANC_GFOFFSET3D(u,0,0,1),kmadd(ToReal(-12),KRANC_GFOFFSET3D(u,0,0,2),kmadd(ToReal(8),KRANC_GFOFFSET3D(u,0,0,3),kmul(KRANC_GFOFFSET3D(u,0,0,4),ToReal(-3))))),kmul(KRANC_GFOFFSET3D(u,0,0,5),ToReal(12)))),kmul(p1o60dz,dir3)))
#else
#  define PDonesidedfdOrder83(u) (PDonesidedfdOrder83_impl(u,p1o60dz,cdj,cdk))
static CCTK_REAL_VEC PDonesidedfdOrder83_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o60dz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDonesidedfdOrder83_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o60dz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{ assert(0); return ToReal(1e30); /* ERROR */ }
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDdissipationNthfdOrder21(u) (kmul(pm1o16dx,kmadd(ToReal(6),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(-4),kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),kadd(KRANC_GFOFFSET3D(u,-2,0,0),KRANC_GFOFFSET3D(u,2,0,0))))))
#else
#  define PDdissipationNthfdOrder21(u) (PDdissipationNthfdOrder21_impl(u,pm1o16dx,cdj,cdk))
static CCTK_REAL_VEC PDdissipationNthfdOrder21_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o16dx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDdissipationNthfdOrder21_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o16dx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(pm1o16dx,kmadd(ToReal(6),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(-4),kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),kadd(KRANC_GFOFFSET3D(u,-2,0,0),KRANC_GFOFFSET3D(u,2,0,0)))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDdissipationNthfdOrder22(u) (kmul(pm1o16dy,kmadd(ToReal(6),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(-4),kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),kadd(KRANC_GFOFFSET3D(u,0,-2,0),KRANC_GFOFFSET3D(u,0,2,0))))))
#else
#  define PDdissipationNthfdOrder22(u) (PDdissipationNthfdOrder22_impl(u,pm1o16dy,cdj,cdk))
static CCTK_REAL_VEC PDdissipationNthfdOrder22_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o16dy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDdissipationNthfdOrder22_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o16dy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(pm1o16dy,kmadd(ToReal(6),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(-4),kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),kadd(KRANC_GFOFFSET3D(u,0,-2,0),KRANC_GFOFFSET3D(u,0,2,0)))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDdissipationNthfdOrder23(u) (kmul(pm1o16dz,kmadd(ToReal(6),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(-4),kadd(KRANC_GFOFFSET3D(u,0,0,-1),KRANC_GFOFFSET3D(u,0,0,1)),kadd(KRANC_GFOFFSET3D(u,0,0,-2),KRANC_GFOFFSET3D(u,0,0,2))))))
#else
#  define PDdissipationNthfdOrder23(u) (PDdissipationNthfdOrder23_impl(u,pm1o16dz,cdj,cdk))
static CCTK_REAL_VEC PDdissipationNthfdOrder23_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o16dz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDdissipationNthfdOrder23_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o16dz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDdissipationNthfdOrder22_impl(u, pm1o16dz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDdissipationNthfdOrder41(u) (kmul(p1o64dx,kmadd(ToReal(-20),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(15),kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),kmadd(ToReal(-6),kadd(KRANC_GFOFFSET3D(u,-2,0,0),KRANC_GFOFFSET3D(u,2,0,0)),kadd(KRANC_GFOFFSET3D(u,-3,0,0),KRANC_GFOFFSET3D(u,3,0,0)))))))
#else
#  define PDdissipationNthfdOrder41(u) (PDdissipationNthfdOrder41_impl(u,p1o64dx,cdj,cdk))
static CCTK_REAL_VEC PDdissipationNthfdOrder41_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o64dx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDdissipationNthfdOrder41_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o64dx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o64dx,kmadd(ToReal(-20),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(15),kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),kmadd(ToReal(-6),kadd(KRANC_GFOFFSET3D(u,-2,0,0),KRANC_GFOFFSET3D(u,2,0,0)),kadd(KRANC_GFOFFSET3D(u,-3,0,0),KRANC_GFOFFSET3D(u,3,0,0))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDdissipationNthfdOrder42(u) (kmul(p1o64dy,kmadd(ToReal(-20),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(15),kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),kmadd(ToReal(-6),kadd(KRANC_GFOFFSET3D(u,0,-2,0),KRANC_GFOFFSET3D(u,0,2,0)),kadd(KRANC_GFOFFSET3D(u,0,-3,0),KRANC_GFOFFSET3D(u,0,3,0)))))))
#else
#  define PDdissipationNthfdOrder42(u) (PDdissipationNthfdOrder42_impl(u,p1o64dy,cdj,cdk))
static CCTK_REAL_VEC PDdissipationNthfdOrder42_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o64dy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDdissipationNthfdOrder42_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o64dy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o64dy,kmadd(ToReal(-20),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(15),kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),kmadd(ToReal(-6),kadd(KRANC_GFOFFSET3D(u,0,-2,0),KRANC_GFOFFSET3D(u,0,2,0)),kadd(KRANC_GFOFFSET3D(u,0,-3,0),KRANC_GFOFFSET3D(u,0,3,0))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDdissipationNthfdOrder43(u) (kmul(p1o64dz,kmadd(ToReal(-20),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(15),kadd(KRANC_GFOFFSET3D(u,0,0,-1),KRANC_GFOFFSET3D(u,0,0,1)),kmadd(ToReal(-6),kadd(KRANC_GFOFFSET3D(u,0,0,-2),KRANC_GFOFFSET3D(u,0,0,2)),kadd(KRANC_GFOFFSET3D(u,0,0,-3),KRANC_GFOFFSET3D(u,0,0,3)))))))
#else
#  define PDdissipationNthfdOrder43(u) (PDdissipationNthfdOrder43_impl(u,p1o64dz,cdj,cdk))
static CCTK_REAL_VEC PDdissipationNthfdOrder43_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o64dz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDdissipationNthfdOrder43_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o64dz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDdissipationNthfdOrder42_impl(u, p1o64dz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDdissipationNthfdOrder61(u) (kmul(pm1o256dx,kmadd(ToReal(70),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(-56),kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),kmadd(ToReal(28),kadd(KRANC_GFOFFSET3D(u,-2,0,0),KRANC_GFOFFSET3D(u,2,0,0)),kmadd(ToReal(-8),kadd(KRANC_GFOFFSET3D(u,-3,0,0),KRANC_GFOFFSET3D(u,3,0,0)),kadd(KRANC_GFOFFSET3D(u,-4,0,0),KRANC_GFOFFSET3D(u,4,0,0))))))))
#else
#  define PDdissipationNthfdOrder61(u) (PDdissipationNthfdOrder61_impl(u,pm1o256dx,cdj,cdk))
static CCTK_REAL_VEC PDdissipationNthfdOrder61_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o256dx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDdissipationNthfdOrder61_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o256dx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(pm1o256dx,kmadd(ToReal(70),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(-56),kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),kmadd(ToReal(28),kadd(KRANC_GFOFFSET3D(u,-2,0,0),KRANC_GFOFFSET3D(u,2,0,0)),kmadd(ToReal(-8),kadd(KRANC_GFOFFSET3D(u,-3,0,0),KRANC_GFOFFSET3D(u,3,0,0)),kadd(KRANC_GFOFFSET3D(u,-4,0,0),KRANC_GFOFFSET3D(u,4,0,0)))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDdissipationNthfdOrder62(u) (kmul(pm1o256dy,kmadd(ToReal(70),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(-56),kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),kmadd(ToReal(28),kadd(KRANC_GFOFFSET3D(u,0,-2,0),KRANC_GFOFFSET3D(u,0,2,0)),kmadd(ToReal(-8),kadd(KRANC_GFOFFSET3D(u,0,-3,0),KRANC_GFOFFSET3D(u,0,3,0)),kadd(KRANC_GFOFFSET3D(u,0,-4,0),KRANC_GFOFFSET3D(u,0,4,0))))))))
#else
#  define PDdissipationNthfdOrder62(u) (PDdissipationNthfdOrder62_impl(u,pm1o256dy,cdj,cdk))
static CCTK_REAL_VEC PDdissipationNthfdOrder62_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o256dy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDdissipationNthfdOrder62_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o256dy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(pm1o256dy,kmadd(ToReal(70),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(-56),kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),kmadd(ToReal(28),kadd(KRANC_GFOFFSET3D(u,0,-2,0),KRANC_GFOFFSET3D(u,0,2,0)),kmadd(ToReal(-8),kadd(KRANC_GFOFFSET3D(u,0,-3,0),KRANC_GFOFFSET3D(u,0,3,0)),kadd(KRANC_GFOFFSET3D(u,0,-4,0),KRANC_GFOFFSET3D(u,0,4,0)))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDdissipationNthfdOrder63(u) (kmul(pm1o256dz,kmadd(ToReal(70),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(-56),kadd(KRANC_GFOFFSET3D(u,0,0,-1),KRANC_GFOFFSET3D(u,0,0,1)),kmadd(ToReal(28),kadd(KRANC_GFOFFSET3D(u,0,0,-2),KRANC_GFOFFSET3D(u,0,0,2)),kmadd(ToReal(-8),kadd(KRANC_GFOFFSET3D(u,0,0,-3),KRANC_GFOFFSET3D(u,0,0,3)),kadd(KRANC_GFOFFSET3D(u,0,0,-4),KRANC_GFOFFSET3D(u,0,0,4))))))))
#else
#  define PDdissipationNthfdOrder63(u) (PDdissipationNthfdOrder63_impl(u,pm1o256dz,cdj,cdk))
static CCTK_REAL_VEC PDdissipationNthfdOrder63_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o256dz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDdissipationNthfdOrder63_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o256dz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDdissipationNthfdOrder62_impl(u, pm1o256dz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDdissipationNthfdOrder81(u) (kmul(p1o1024dx,kmadd(ToReal(-252),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(210),kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),kmadd(ToReal(-120),kadd(KRANC_GFOFFSET3D(u,-2,0,0),KRANC_GFOFFSET3D(u,2,0,0)),kmadd(ToReal(45),kadd(KRANC_GFOFFSET3D(u,-3,0,0),KRANC_GFOFFSET3D(u,3,0,0)),kmadd(ToReal(-10),kadd(KRANC_GFOFFSET3D(u,-4,0,0),KRANC_GFOFFSET3D(u,4,0,0)),kadd(KRANC_GFOFFSET3D(u,-5,0,0),KRANC_GFOFFSET3D(u,5,0,0)))))))))
#else
#  define PDdissipationNthfdOrder81(u) (PDdissipationNthfdOrder81_impl(u,p1o1024dx,cdj,cdk))
static CCTK_REAL_VEC PDdissipationNthfdOrder81_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o1024dx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDdissipationNthfdOrder81_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o1024dx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o1024dx,kmadd(ToReal(-252),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(210),kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),kmadd(ToReal(-120),kadd(KRANC_GFOFFSET3D(u,-2,0,0),KRANC_GFOFFSET3D(u,2,0,0)),kmadd(ToReal(45),kadd(KRANC_GFOFFSET3D(u,-3,0,0),KRANC_GFOFFSET3D(u,3,0,0)),kmadd(ToReal(-10),kadd(KRANC_GFOFFSET3D(u,-4,0,0),KRANC_GFOFFSET3D(u,4,0,0)),kadd(KRANC_GFOFFSET3D(u,-5,0,0),KRANC_GFOFFSET3D(u,5,0,0))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDdissipationNthfdOrder82(u) (kmul(p1o1024dy,kmadd(ToReal(-252),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(210),kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),kmadd(ToReal(-120),kadd(KRANC_GFOFFSET3D(u,0,-2,0),KRANC_GFOFFSET3D(u,0,2,0)),kmadd(ToReal(45),kadd(KRANC_GFOFFSET3D(u,0,-3,0),KRANC_GFOFFSET3D(u,0,3,0)),kmadd(ToReal(-10),kadd(KRANC_GFOFFSET3D(u,0,-4,0),KRANC_GFOFFSET3D(u,0,4,0)),kadd(KRANC_GFOFFSET3D(u,0,-5,0),KRANC_GFOFFSET3D(u,0,5,0)))))))))
#else
#  define PDdissipationNthfdOrder82(u) (PDdissipationNthfdOrder82_impl(u,p1o1024dy,cdj,cdk))
static CCTK_REAL_VEC PDdissipationNthfdOrder82_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o1024dy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDdissipationNthfdOrder82_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o1024dy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o1024dy,kmadd(ToReal(-252),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(210),kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),kmadd(ToReal(-120),kadd(KRANC_GFOFFSET3D(u,0,-2,0),KRANC_GFOFFSET3D(u,0,2,0)),kmadd(ToReal(45),kadd(KRANC_GFOFFSET3D(u,0,-3,0),KRANC_GFOFFSET3D(u,0,3,0)),kmadd(ToReal(-10),kadd(KRANC_GFOFFSET3D(u,0,-4,0),KRANC_GFOFFSET3D(u,0,4,0)),kadd(KRANC_GFOFFSET3D(u,0,-5,0),KRANC_GFOFFSET3D(u,0,5,0))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDdissipationNthfdOrder83(u) (kmul(p1o1024dz,kmadd(ToReal(-252),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(210),kadd(KRANC_GFOFFSET3D(u,0,0,-1),KRANC_GFOFFSET3D(u,0,0,1)),kmadd(ToReal(-120),kadd(KRANC_GFOFFSET3D(u,0,0,-2),KRANC_GFOFFSET3D(u,0,0,2)),kmadd(ToReal(45),kadd(KRANC_GFOFFSET3D(u,0,0,-3),KRANC_GFOFFSET3D(u,0,0,3)),kmadd(ToReal(-10),kadd(KRANC_GFOFFSET3D(u,0,0,-4),KRANC_GFOFFSET3D(u,0,0,4)),kadd(KRANC_GFOFFSET3D(u,0,0,-5),KRANC_GFOFFSET3D(u,0,0,5)))))))))
#else
#  define PDdissipationNthfdOrder83(u) (PDdissipationNthfdOrder83_impl(u,p1o1024dz,cdj,cdk))
static CCTK_REAL_VEC PDdissipationNthfdOrder83_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o1024dz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDdissipationNthfdOrder83_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o1024dz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDdissipationNthfdOrder82_impl(u, p1o1024dz, cdk, cdj);
}
#endif

