#include <assert.h>
#include "vectors.h"

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard2nd1(u) (kmul(p1o2dx,ksub(KRANC_GFOFFSET3D(u,1,0,0),KRANC_GFOFFSET3D(u,-1,0,0))))
#else
#  define PDstandard2nd1(u) (PDstandard2nd1_impl(u,p1o2dx,cdj,cdk))
static CCTK_REAL_VEC PDstandard2nd1_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o2dx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandard2nd1_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o2dx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o2dx,ksub(KRANC_GFOFFSET3D(u,1,0,0),KRANC_GFOFFSET3D(u,-1,0,0)));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard2nd2(u) (kmul(p1o2dy,ksub(KRANC_GFOFFSET3D(u,0,1,0),KRANC_GFOFFSET3D(u,0,-1,0))))
#else
#  define PDstandard2nd2(u) (PDstandard2nd2_impl(u,p1o2dy,cdj,cdk))
static CCTK_REAL_VEC PDstandard2nd2_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o2dy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandard2nd2_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o2dy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o2dy,ksub(KRANC_GFOFFSET3D(u,0,1,0),KRANC_GFOFFSET3D(u,0,-1,0)));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard2nd3(u) (kmul(p1o2dz,ksub(KRANC_GFOFFSET3D(u,0,0,1),KRANC_GFOFFSET3D(u,0,0,-1))))
#else
#  define PDstandard2nd3(u) (PDstandard2nd3_impl(u,p1o2dz,cdj,cdk))
static CCTK_REAL_VEC PDstandard2nd3_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o2dz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandard2nd3_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o2dz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandard2nd2_impl(u, p1o2dz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard2nd11(u) (kmul(p1odx2,kmadd(ToReal(-2),KRANC_GFOFFSET3D(u,0,0,0),kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)))))
#else
#  define PDstandard2nd11(u) (PDstandard2nd11_impl(u,p1odx2,cdj,cdk))
static CCTK_REAL_VEC PDstandard2nd11_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1odx2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandard2nd11_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1odx2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1odx2,kmadd(ToReal(-2),KRANC_GFOFFSET3D(u,0,0,0),kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard2nd22(u) (kmul(p1ody2,kmadd(ToReal(-2),KRANC_GFOFFSET3D(u,0,0,0),kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)))))
#else
#  define PDstandard2nd22(u) (PDstandard2nd22_impl(u,p1ody2,cdj,cdk))
static CCTK_REAL_VEC PDstandard2nd22_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1ody2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandard2nd22_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1ody2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1ody2,kmadd(ToReal(-2),KRANC_GFOFFSET3D(u,0,0,0),kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard2nd33(u) (kmul(p1odz2,kmadd(ToReal(-2),KRANC_GFOFFSET3D(u,0,0,0),kadd(KRANC_GFOFFSET3D(u,0,0,-1),KRANC_GFOFFSET3D(u,0,0,1)))))
#else
#  define PDstandard2nd33(u) (PDstandard2nd33_impl(u,p1odz2,cdj,cdk))
static CCTK_REAL_VEC PDstandard2nd33_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1odz2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandard2nd33_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1odz2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandard2nd22_impl(u, p1odz2, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard2nd12(u) (kmul(p1o4dxdy,kadd(KRANC_GFOFFSET3D(u,-1,-1,0),ksub(ksub(KRANC_GFOFFSET3D(u,1,1,0),KRANC_GFOFFSET3D(u,1,-1,0)),KRANC_GFOFFSET3D(u,-1,1,0)))))
#else
#  define PDstandard2nd12(u) (PDstandard2nd12_impl(u,p1o4dxdy,cdj,cdk))
static CCTK_REAL_VEC PDstandard2nd12_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o4dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandard2nd12_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o4dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o4dxdy,kadd(KRANC_GFOFFSET3D(u,-1,-1,0),ksub(ksub(KRANC_GFOFFSET3D(u,1,1,0),KRANC_GFOFFSET3D(u,1,-1,0)),KRANC_GFOFFSET3D(u,-1,1,0))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard2nd13(u) (kmul(p1o4dxdz,kadd(KRANC_GFOFFSET3D(u,-1,0,-1),ksub(ksub(KRANC_GFOFFSET3D(u,1,0,1),KRANC_GFOFFSET3D(u,1,0,-1)),KRANC_GFOFFSET3D(u,-1,0,1)))))
#else
#  define PDstandard2nd13(u) (PDstandard2nd13_impl(u,p1o4dxdz,cdj,cdk))
static CCTK_REAL_VEC PDstandard2nd13_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o4dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandard2nd13_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o4dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandard2nd12_impl(u, p1o4dxdz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard2nd21(u) (kmul(p1o4dxdy,kadd(KRANC_GFOFFSET3D(u,-1,-1,0),ksub(ksub(KRANC_GFOFFSET3D(u,1,1,0),KRANC_GFOFFSET3D(u,1,-1,0)),KRANC_GFOFFSET3D(u,-1,1,0)))))
#else
#  define PDstandard2nd21(u) (PDstandard2nd21_impl(u,p1o4dxdy,cdj,cdk))
static CCTK_REAL_VEC PDstandard2nd21_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o4dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandard2nd21_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o4dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandard2nd12_impl(u, p1o4dxdy, cdj, cdk);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard2nd23(u) (kmul(p1o4dydz,kadd(KRANC_GFOFFSET3D(u,0,-1,-1),ksub(ksub(KRANC_GFOFFSET3D(u,0,1,1),KRANC_GFOFFSET3D(u,0,1,-1)),KRANC_GFOFFSET3D(u,0,-1,1)))))
#else
#  define PDstandard2nd23(u) (PDstandard2nd23_impl(u,p1o4dydz,cdj,cdk))
static CCTK_REAL_VEC PDstandard2nd23_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o4dydz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandard2nd23_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o4dydz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o4dydz,kadd(KRANC_GFOFFSET3D(u,0,-1,-1),ksub(ksub(KRANC_GFOFFSET3D(u,0,1,1),KRANC_GFOFFSET3D(u,0,1,-1)),KRANC_GFOFFSET3D(u,0,-1,1))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard2nd31(u) (kmul(p1o4dxdz,kadd(KRANC_GFOFFSET3D(u,-1,0,-1),ksub(ksub(KRANC_GFOFFSET3D(u,1,0,1),KRANC_GFOFFSET3D(u,1,0,-1)),KRANC_GFOFFSET3D(u,-1,0,1)))))
#else
#  define PDstandard2nd31(u) (PDstandard2nd31_impl(u,p1o4dxdz,cdj,cdk))
static CCTK_REAL_VEC PDstandard2nd31_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o4dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandard2nd31_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o4dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandard2nd12_impl(u, p1o4dxdz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard2nd32(u) (kmul(p1o4dydz,kadd(KRANC_GFOFFSET3D(u,0,-1,-1),ksub(ksub(KRANC_GFOFFSET3D(u,0,1,1),KRANC_GFOFFSET3D(u,0,1,-1)),KRANC_GFOFFSET3D(u,0,-1,1)))))
#else
#  define PDstandard2nd32(u) (PDstandard2nd32_impl(u,p1o4dydz,cdj,cdk))
static CCTK_REAL_VEC PDstandard2nd32_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o4dydz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandard2nd32_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o4dydz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandard2nd23_impl(u, p1o4dydz, cdj, cdk);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard4th1(u) (kmul(p1o12dx,kmadd(ToReal(-8),KRANC_GFOFFSET3D(u,-1,0,0),kmadd(ToReal(8),KRANC_GFOFFSET3D(u,1,0,0),ksub(KRANC_GFOFFSET3D(u,-2,0,0),KRANC_GFOFFSET3D(u,2,0,0))))))
#else
#  define PDstandard4th1(u) (PDstandard4th1_impl(u,p1o12dx,cdj,cdk))
static CCTK_REAL_VEC PDstandard4th1_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o12dx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandard4th1_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o12dx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o12dx,kmadd(ToReal(-8),KRANC_GFOFFSET3D(u,-1,0,0),kmadd(ToReal(8),KRANC_GFOFFSET3D(u,1,0,0),ksub(KRANC_GFOFFSET3D(u,-2,0,0),KRANC_GFOFFSET3D(u,2,0,0)))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard4th2(u) (kmul(p1o12dy,kmadd(ToReal(-8),KRANC_GFOFFSET3D(u,0,-1,0),kmadd(ToReal(8),KRANC_GFOFFSET3D(u,0,1,0),ksub(KRANC_GFOFFSET3D(u,0,-2,0),KRANC_GFOFFSET3D(u,0,2,0))))))
#else
#  define PDstandard4th2(u) (PDstandard4th2_impl(u,p1o12dy,cdj,cdk))
static CCTK_REAL_VEC PDstandard4th2_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o12dy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandard4th2_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o12dy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o12dy,kmadd(ToReal(-8),KRANC_GFOFFSET3D(u,0,-1,0),kmadd(ToReal(8),KRANC_GFOFFSET3D(u,0,1,0),ksub(KRANC_GFOFFSET3D(u,0,-2,0),KRANC_GFOFFSET3D(u,0,2,0)))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard4th3(u) (kmul(p1o12dz,kmadd(ToReal(-8),KRANC_GFOFFSET3D(u,0,0,-1),kmadd(ToReal(8),KRANC_GFOFFSET3D(u,0,0,1),ksub(KRANC_GFOFFSET3D(u,0,0,-2),KRANC_GFOFFSET3D(u,0,0,2))))))
#else
#  define PDstandard4th3(u) (PDstandard4th3_impl(u,p1o12dz,cdj,cdk))
static CCTK_REAL_VEC PDstandard4th3_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o12dz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandard4th3_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o12dz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandard4th2_impl(u, p1o12dz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard4th11(u) (kmul(pm1o12dx2,kmadd(ToReal(30),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(-16),kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),kadd(KRANC_GFOFFSET3D(u,-2,0,0),KRANC_GFOFFSET3D(u,2,0,0))))))
#else
#  define PDstandard4th11(u) (PDstandard4th11_impl(u,pm1o12dx2,cdj,cdk))
static CCTK_REAL_VEC PDstandard4th11_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o12dx2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandard4th11_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o12dx2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(pm1o12dx2,kmadd(ToReal(30),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(-16),kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),kadd(KRANC_GFOFFSET3D(u,-2,0,0),KRANC_GFOFFSET3D(u,2,0,0)))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard4th22(u) (kmul(pm1o12dy2,kmadd(ToReal(30),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(-16),kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),kadd(KRANC_GFOFFSET3D(u,0,-2,0),KRANC_GFOFFSET3D(u,0,2,0))))))
#else
#  define PDstandard4th22(u) (PDstandard4th22_impl(u,pm1o12dy2,cdj,cdk))
static CCTK_REAL_VEC PDstandard4th22_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o12dy2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandard4th22_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o12dy2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(pm1o12dy2,kmadd(ToReal(30),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(-16),kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),kadd(KRANC_GFOFFSET3D(u,0,-2,0),KRANC_GFOFFSET3D(u,0,2,0)))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard4th33(u) (kmul(pm1o12dz2,kmadd(ToReal(30),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(-16),kadd(KRANC_GFOFFSET3D(u,0,0,-1),KRANC_GFOFFSET3D(u,0,0,1)),kadd(KRANC_GFOFFSET3D(u,0,0,-2),KRANC_GFOFFSET3D(u,0,0,2))))))
#else
#  define PDstandard4th33(u) (PDstandard4th33_impl(u,pm1o12dz2,cdj,cdk))
static CCTK_REAL_VEC PDstandard4th33_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o12dz2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandard4th33_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o12dz2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandard4th22_impl(u, pm1o12dz2, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard4th12(u) (kmul(p1o144dxdy,kmadd(ToReal(-64),kadd(KRANC_GFOFFSET3D(u,-1,1,0),KRANC_GFOFFSET3D(u,1,-1,0)),kmadd(ToReal(64),kadd(KRANC_GFOFFSET3D(u,-1,-1,0),KRANC_GFOFFSET3D(u,1,1,0)),kmadd(ToReal(8),kadd(KRANC_GFOFFSET3D(u,-1,2,0),kadd(KRANC_GFOFFSET3D(u,1,-2,0),kadd(KRANC_GFOFFSET3D(u,-2,1,0),KRANC_GFOFFSET3D(u,2,-1,0)))),kmadd(ToReal(-8),kadd(KRANC_GFOFFSET3D(u,-1,-2,0),kadd(KRANC_GFOFFSET3D(u,1,2,0),kadd(KRANC_GFOFFSET3D(u,-2,-1,0),KRANC_GFOFFSET3D(u,2,1,0)))),kadd(KRANC_GFOFFSET3D(u,-2,-2,0),ksub(ksub(KRANC_GFOFFSET3D(u,2,2,0),KRANC_GFOFFSET3D(u,2,-2,0)),KRANC_GFOFFSET3D(u,-2,2,0)))))))))
#else
#  define PDstandard4th12(u) (PDstandard4th12_impl(u,p1o144dxdy,cdj,cdk))
static CCTK_REAL_VEC PDstandard4th12_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o144dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandard4th12_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o144dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o144dxdy,kmadd(ToReal(-64),kadd(KRANC_GFOFFSET3D(u,-1,1,0),KRANC_GFOFFSET3D(u,1,-1,0)),kmadd(ToReal(64),kadd(KRANC_GFOFFSET3D(u,-1,-1,0),KRANC_GFOFFSET3D(u,1,1,0)),kmadd(ToReal(8),kadd(KRANC_GFOFFSET3D(u,-1,2,0),kadd(KRANC_GFOFFSET3D(u,1,-2,0),kadd(KRANC_GFOFFSET3D(u,-2,1,0),KRANC_GFOFFSET3D(u,2,-1,0)))),kmadd(ToReal(-8),kadd(KRANC_GFOFFSET3D(u,-1,-2,0),kadd(KRANC_GFOFFSET3D(u,1,2,0),kadd(KRANC_GFOFFSET3D(u,-2,-1,0),KRANC_GFOFFSET3D(u,2,1,0)))),kadd(KRANC_GFOFFSET3D(u,-2,-2,0),ksub(ksub(KRANC_GFOFFSET3D(u,2,2,0),KRANC_GFOFFSET3D(u,2,-2,0)),KRANC_GFOFFSET3D(u,-2,2,0))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard4th13(u) (kmul(p1o144dxdz,kmadd(ToReal(-64),kadd(KRANC_GFOFFSET3D(u,-1,0,1),KRANC_GFOFFSET3D(u,1,0,-1)),kmadd(ToReal(64),kadd(KRANC_GFOFFSET3D(u,-1,0,-1),KRANC_GFOFFSET3D(u,1,0,1)),kmadd(ToReal(8),kadd(KRANC_GFOFFSET3D(u,-1,0,2),kadd(KRANC_GFOFFSET3D(u,1,0,-2),kadd(KRANC_GFOFFSET3D(u,-2,0,1),KRANC_GFOFFSET3D(u,2,0,-1)))),kmadd(ToReal(-8),kadd(KRANC_GFOFFSET3D(u,-1,0,-2),kadd(KRANC_GFOFFSET3D(u,1,0,2),kadd(KRANC_GFOFFSET3D(u,-2,0,-1),KRANC_GFOFFSET3D(u,2,0,1)))),kadd(KRANC_GFOFFSET3D(u,-2,0,-2),ksub(ksub(KRANC_GFOFFSET3D(u,2,0,2),KRANC_GFOFFSET3D(u,2,0,-2)),KRANC_GFOFFSET3D(u,-2,0,2)))))))))
#else
#  define PDstandard4th13(u) (PDstandard4th13_impl(u,p1o144dxdz,cdj,cdk))
static CCTK_REAL_VEC PDstandard4th13_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o144dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandard4th13_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o144dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandard4th12_impl(u, p1o144dxdz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard4th21(u) (kmul(p1o144dxdy,kmadd(ToReal(-64),kadd(KRANC_GFOFFSET3D(u,-1,1,0),KRANC_GFOFFSET3D(u,1,-1,0)),kmadd(ToReal(64),kadd(KRANC_GFOFFSET3D(u,-1,-1,0),KRANC_GFOFFSET3D(u,1,1,0)),kmadd(ToReal(8),kadd(KRANC_GFOFFSET3D(u,-1,2,0),kadd(KRANC_GFOFFSET3D(u,1,-2,0),kadd(KRANC_GFOFFSET3D(u,-2,1,0),KRANC_GFOFFSET3D(u,2,-1,0)))),kmadd(ToReal(-8),kadd(KRANC_GFOFFSET3D(u,-1,-2,0),kadd(KRANC_GFOFFSET3D(u,1,2,0),kadd(KRANC_GFOFFSET3D(u,-2,-1,0),KRANC_GFOFFSET3D(u,2,1,0)))),kadd(KRANC_GFOFFSET3D(u,-2,-2,0),ksub(ksub(KRANC_GFOFFSET3D(u,2,2,0),KRANC_GFOFFSET3D(u,2,-2,0)),KRANC_GFOFFSET3D(u,-2,2,0)))))))))
#else
#  define PDstandard4th21(u) (PDstandard4th21_impl(u,p1o144dxdy,cdj,cdk))
static CCTK_REAL_VEC PDstandard4th21_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o144dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandard4th21_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o144dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandard4th12_impl(u, p1o144dxdy, cdj, cdk);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard4th23(u) (kmul(p1o144dydz,kmadd(ToReal(-64),kadd(KRANC_GFOFFSET3D(u,0,-1,1),KRANC_GFOFFSET3D(u,0,1,-1)),kmadd(ToReal(64),kadd(KRANC_GFOFFSET3D(u,0,-1,-1),KRANC_GFOFFSET3D(u,0,1,1)),kmadd(ToReal(8),kadd(KRANC_GFOFFSET3D(u,0,-1,2),kadd(KRANC_GFOFFSET3D(u,0,1,-2),kadd(KRANC_GFOFFSET3D(u,0,-2,1),KRANC_GFOFFSET3D(u,0,2,-1)))),kmadd(ToReal(-8),kadd(KRANC_GFOFFSET3D(u,0,-1,-2),kadd(KRANC_GFOFFSET3D(u,0,1,2),kadd(KRANC_GFOFFSET3D(u,0,-2,-1),KRANC_GFOFFSET3D(u,0,2,1)))),kadd(KRANC_GFOFFSET3D(u,0,-2,-2),ksub(ksub(KRANC_GFOFFSET3D(u,0,2,2),KRANC_GFOFFSET3D(u,0,2,-2)),KRANC_GFOFFSET3D(u,0,-2,2)))))))))
#else
#  define PDstandard4th23(u) (PDstandard4th23_impl(u,p1o144dydz,cdj,cdk))
static CCTK_REAL_VEC PDstandard4th23_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o144dydz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandard4th23_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o144dydz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o144dydz,kmadd(ToReal(-64),kadd(KRANC_GFOFFSET3D(u,0,-1,1),KRANC_GFOFFSET3D(u,0,1,-1)),kmadd(ToReal(64),kadd(KRANC_GFOFFSET3D(u,0,-1,-1),KRANC_GFOFFSET3D(u,0,1,1)),kmadd(ToReal(8),kadd(KRANC_GFOFFSET3D(u,0,-1,2),kadd(KRANC_GFOFFSET3D(u,0,1,-2),kadd(KRANC_GFOFFSET3D(u,0,-2,1),KRANC_GFOFFSET3D(u,0,2,-1)))),kmadd(ToReal(-8),kadd(KRANC_GFOFFSET3D(u,0,-1,-2),kadd(KRANC_GFOFFSET3D(u,0,1,2),kadd(KRANC_GFOFFSET3D(u,0,-2,-1),KRANC_GFOFFSET3D(u,0,2,1)))),kadd(KRANC_GFOFFSET3D(u,0,-2,-2),ksub(ksub(KRANC_GFOFFSET3D(u,0,2,2),KRANC_GFOFFSET3D(u,0,2,-2)),KRANC_GFOFFSET3D(u,0,-2,2))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard4th31(u) (kmul(p1o144dxdz,kmadd(ToReal(-64),kadd(KRANC_GFOFFSET3D(u,-1,0,1),KRANC_GFOFFSET3D(u,1,0,-1)),kmadd(ToReal(64),kadd(KRANC_GFOFFSET3D(u,-1,0,-1),KRANC_GFOFFSET3D(u,1,0,1)),kmadd(ToReal(8),kadd(KRANC_GFOFFSET3D(u,-1,0,2),kadd(KRANC_GFOFFSET3D(u,1,0,-2),kadd(KRANC_GFOFFSET3D(u,-2,0,1),KRANC_GFOFFSET3D(u,2,0,-1)))),kmadd(ToReal(-8),kadd(KRANC_GFOFFSET3D(u,-1,0,-2),kadd(KRANC_GFOFFSET3D(u,1,0,2),kadd(KRANC_GFOFFSET3D(u,-2,0,-1),KRANC_GFOFFSET3D(u,2,0,1)))),kadd(KRANC_GFOFFSET3D(u,-2,0,-2),ksub(ksub(KRANC_GFOFFSET3D(u,2,0,2),KRANC_GFOFFSET3D(u,2,0,-2)),KRANC_GFOFFSET3D(u,-2,0,2)))))))))
#else
#  define PDstandard4th31(u) (PDstandard4th31_impl(u,p1o144dxdz,cdj,cdk))
static CCTK_REAL_VEC PDstandard4th31_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o144dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandard4th31_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o144dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandard4th12_impl(u, p1o144dxdz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard4th32(u) (kmul(p1o144dydz,kmadd(ToReal(-64),kadd(KRANC_GFOFFSET3D(u,0,-1,1),KRANC_GFOFFSET3D(u,0,1,-1)),kmadd(ToReal(64),kadd(KRANC_GFOFFSET3D(u,0,-1,-1),KRANC_GFOFFSET3D(u,0,1,1)),kmadd(ToReal(8),kadd(KRANC_GFOFFSET3D(u,0,-1,2),kadd(KRANC_GFOFFSET3D(u,0,1,-2),kadd(KRANC_GFOFFSET3D(u,0,-2,1),KRANC_GFOFFSET3D(u,0,2,-1)))),kmadd(ToReal(-8),kadd(KRANC_GFOFFSET3D(u,0,-1,-2),kadd(KRANC_GFOFFSET3D(u,0,1,2),kadd(KRANC_GFOFFSET3D(u,0,-2,-1),KRANC_GFOFFSET3D(u,0,2,1)))),kadd(KRANC_GFOFFSET3D(u,0,-2,-2),ksub(ksub(KRANC_GFOFFSET3D(u,0,2,2),KRANC_GFOFFSET3D(u,0,2,-2)),KRANC_GFOFFSET3D(u,0,-2,2)))))))))
#else
#  define PDstandard4th32(u) (PDstandard4th32_impl(u,p1o144dydz,cdj,cdk))
static CCTK_REAL_VEC PDstandard4th32_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o144dydz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandard4th32_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o144dydz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandard4th23_impl(u, p1o144dydz, cdj, cdk);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder21(u) (kmul(p1o2dx,ksub(KRANC_GFOFFSET3D(u,1,0,0),KRANC_GFOFFSET3D(u,-1,0,0))))
#else
#  define PDstandardfdOrder21(u) (PDstandardfdOrder21_impl(u,p1o2dx,cdj,cdk))
static CCTK_REAL_VEC PDstandardfdOrder21_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o2dx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardfdOrder21_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o2dx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o2dx,ksub(KRANC_GFOFFSET3D(u,1,0,0),KRANC_GFOFFSET3D(u,-1,0,0)));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder22(u) (kmul(p1o2dy,ksub(KRANC_GFOFFSET3D(u,0,1,0),KRANC_GFOFFSET3D(u,0,-1,0))))
#else
#  define PDstandardfdOrder22(u) (PDstandardfdOrder22_impl(u,p1o2dy,cdj,cdk))
static CCTK_REAL_VEC PDstandardfdOrder22_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o2dy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardfdOrder22_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o2dy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o2dy,ksub(KRANC_GFOFFSET3D(u,0,1,0),KRANC_GFOFFSET3D(u,0,-1,0)));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder23(u) (kmul(p1o2dz,ksub(KRANC_GFOFFSET3D(u,0,0,1),KRANC_GFOFFSET3D(u,0,0,-1))))
#else
#  define PDstandardfdOrder23(u) (PDstandardfdOrder23_impl(u,p1o2dz,cdj,cdk))
static CCTK_REAL_VEC PDstandardfdOrder23_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o2dz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardfdOrder23_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o2dz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardfdOrder22_impl(u, p1o2dz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder41(u) (kmul(p1o12dx,kmadd(ToReal(-8),KRANC_GFOFFSET3D(u,-1,0,0),kmadd(ToReal(8),KRANC_GFOFFSET3D(u,1,0,0),ksub(KRANC_GFOFFSET3D(u,-2,0,0),KRANC_GFOFFSET3D(u,2,0,0))))))
#else
#  define PDstandardfdOrder41(u) (PDstandardfdOrder41_impl(u,p1o12dx,cdj,cdk))
static CCTK_REAL_VEC PDstandardfdOrder41_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o12dx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardfdOrder41_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o12dx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o12dx,kmadd(ToReal(-8),KRANC_GFOFFSET3D(u,-1,0,0),kmadd(ToReal(8),KRANC_GFOFFSET3D(u,1,0,0),ksub(KRANC_GFOFFSET3D(u,-2,0,0),KRANC_GFOFFSET3D(u,2,0,0)))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder42(u) (kmul(p1o12dy,kmadd(ToReal(-8),KRANC_GFOFFSET3D(u,0,-1,0),kmadd(ToReal(8),KRANC_GFOFFSET3D(u,0,1,0),ksub(KRANC_GFOFFSET3D(u,0,-2,0),KRANC_GFOFFSET3D(u,0,2,0))))))
#else
#  define PDstandardfdOrder42(u) (PDstandardfdOrder42_impl(u,p1o12dy,cdj,cdk))
static CCTK_REAL_VEC PDstandardfdOrder42_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o12dy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardfdOrder42_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o12dy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o12dy,kmadd(ToReal(-8),KRANC_GFOFFSET3D(u,0,-1,0),kmadd(ToReal(8),KRANC_GFOFFSET3D(u,0,1,0),ksub(KRANC_GFOFFSET3D(u,0,-2,0),KRANC_GFOFFSET3D(u,0,2,0)))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder43(u) (kmul(p1o12dz,kmadd(ToReal(-8),KRANC_GFOFFSET3D(u,0,0,-1),kmadd(ToReal(8),KRANC_GFOFFSET3D(u,0,0,1),ksub(KRANC_GFOFFSET3D(u,0,0,-2),KRANC_GFOFFSET3D(u,0,0,2))))))
#else
#  define PDstandardfdOrder43(u) (PDstandardfdOrder43_impl(u,p1o12dz,cdj,cdk))
static CCTK_REAL_VEC PDstandardfdOrder43_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o12dz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardfdOrder43_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o12dz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardfdOrder42_impl(u, p1o12dz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder61(u) (kmul(p1o60dx,kmadd(ToReal(-45),KRANC_GFOFFSET3D(u,-1,0,0),kmadd(ToReal(45),KRANC_GFOFFSET3D(u,1,0,0),kmadd(ToReal(9),KRANC_GFOFFSET3D(u,-2,0,0),kmadd(ToReal(-9),KRANC_GFOFFSET3D(u,2,0,0),ksub(KRANC_GFOFFSET3D(u,3,0,0),KRANC_GFOFFSET3D(u,-3,0,0))))))))
#else
#  define PDstandardfdOrder61(u) (PDstandardfdOrder61_impl(u,p1o60dx,cdj,cdk))
static CCTK_REAL_VEC PDstandardfdOrder61_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o60dx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardfdOrder61_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o60dx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o60dx,kmadd(ToReal(-45),KRANC_GFOFFSET3D(u,-1,0,0),kmadd(ToReal(45),KRANC_GFOFFSET3D(u,1,0,0),kmadd(ToReal(9),KRANC_GFOFFSET3D(u,-2,0,0),kmadd(ToReal(-9),KRANC_GFOFFSET3D(u,2,0,0),ksub(KRANC_GFOFFSET3D(u,3,0,0),KRANC_GFOFFSET3D(u,-3,0,0)))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder62(u) (kmul(p1o60dy,kmadd(ToReal(-45),KRANC_GFOFFSET3D(u,0,-1,0),kmadd(ToReal(45),KRANC_GFOFFSET3D(u,0,1,0),kmadd(ToReal(9),KRANC_GFOFFSET3D(u,0,-2,0),kmadd(ToReal(-9),KRANC_GFOFFSET3D(u,0,2,0),ksub(KRANC_GFOFFSET3D(u,0,3,0),KRANC_GFOFFSET3D(u,0,-3,0))))))))
#else
#  define PDstandardfdOrder62(u) (PDstandardfdOrder62_impl(u,p1o60dy,cdj,cdk))
static CCTK_REAL_VEC PDstandardfdOrder62_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o60dy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardfdOrder62_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o60dy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o60dy,kmadd(ToReal(-45),KRANC_GFOFFSET3D(u,0,-1,0),kmadd(ToReal(45),KRANC_GFOFFSET3D(u,0,1,0),kmadd(ToReal(9),KRANC_GFOFFSET3D(u,0,-2,0),kmadd(ToReal(-9),KRANC_GFOFFSET3D(u,0,2,0),ksub(KRANC_GFOFFSET3D(u,0,3,0),KRANC_GFOFFSET3D(u,0,-3,0)))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder63(u) (kmul(p1o60dz,kmadd(ToReal(-45),KRANC_GFOFFSET3D(u,0,0,-1),kmadd(ToReal(45),KRANC_GFOFFSET3D(u,0,0,1),kmadd(ToReal(9),KRANC_GFOFFSET3D(u,0,0,-2),kmadd(ToReal(-9),KRANC_GFOFFSET3D(u,0,0,2),ksub(KRANC_GFOFFSET3D(u,0,0,3),KRANC_GFOFFSET3D(u,0,0,-3))))))))
#else
#  define PDstandardfdOrder63(u) (PDstandardfdOrder63_impl(u,p1o60dz,cdj,cdk))
static CCTK_REAL_VEC PDstandardfdOrder63_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o60dz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardfdOrder63_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o60dz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardfdOrder62_impl(u, p1o60dz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder81(u) (kmul(p1o840dx,kmadd(ToReal(-672),KRANC_GFOFFSET3D(u,-1,0,0),kmadd(ToReal(672),KRANC_GFOFFSET3D(u,1,0,0),kmadd(ToReal(168),KRANC_GFOFFSET3D(u,-2,0,0),kmadd(ToReal(-168),KRANC_GFOFFSET3D(u,2,0,0),kmadd(ToReal(-32),KRANC_GFOFFSET3D(u,-3,0,0),kmadd(ToReal(32),KRANC_GFOFFSET3D(u,3,0,0),kmadd(ToReal(3),KRANC_GFOFFSET3D(u,-4,0,0),kmul(ToReal(-3),KRANC_GFOFFSET3D(u,4,0,0)))))))))))
#else
#  define PDstandardfdOrder81(u) (PDstandardfdOrder81_impl(u,p1o840dx,cdj,cdk))
static CCTK_REAL_VEC PDstandardfdOrder81_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o840dx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardfdOrder81_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o840dx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o840dx,kmadd(ToReal(-672),KRANC_GFOFFSET3D(u,-1,0,0),kmadd(ToReal(672),KRANC_GFOFFSET3D(u,1,0,0),kmadd(ToReal(168),KRANC_GFOFFSET3D(u,-2,0,0),kmadd(ToReal(-168),KRANC_GFOFFSET3D(u,2,0,0),kmadd(ToReal(-32),KRANC_GFOFFSET3D(u,-3,0,0),kmadd(ToReal(32),KRANC_GFOFFSET3D(u,3,0,0),kmadd(ToReal(3),KRANC_GFOFFSET3D(u,-4,0,0),kmul(ToReal(-3),KRANC_GFOFFSET3D(u,4,0,0))))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder82(u) (kmul(p1o840dy,kmadd(ToReal(-672),KRANC_GFOFFSET3D(u,0,-1,0),kmadd(ToReal(672),KRANC_GFOFFSET3D(u,0,1,0),kmadd(ToReal(168),KRANC_GFOFFSET3D(u,0,-2,0),kmadd(ToReal(-168),KRANC_GFOFFSET3D(u,0,2,0),kmadd(ToReal(-32),KRANC_GFOFFSET3D(u,0,-3,0),kmadd(ToReal(32),KRANC_GFOFFSET3D(u,0,3,0),kmadd(ToReal(3),KRANC_GFOFFSET3D(u,0,-4,0),kmul(ToReal(-3),KRANC_GFOFFSET3D(u,0,4,0)))))))))))
#else
#  define PDstandardfdOrder82(u) (PDstandardfdOrder82_impl(u,p1o840dy,cdj,cdk))
static CCTK_REAL_VEC PDstandardfdOrder82_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o840dy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardfdOrder82_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o840dy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o840dy,kmadd(ToReal(-672),KRANC_GFOFFSET3D(u,0,-1,0),kmadd(ToReal(672),KRANC_GFOFFSET3D(u,0,1,0),kmadd(ToReal(168),KRANC_GFOFFSET3D(u,0,-2,0),kmadd(ToReal(-168),KRANC_GFOFFSET3D(u,0,2,0),kmadd(ToReal(-32),KRANC_GFOFFSET3D(u,0,-3,0),kmadd(ToReal(32),KRANC_GFOFFSET3D(u,0,3,0),kmadd(ToReal(3),KRANC_GFOFFSET3D(u,0,-4,0),kmul(ToReal(-3),KRANC_GFOFFSET3D(u,0,4,0))))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder83(u) (kmul(p1o840dz,kmadd(ToReal(-672),KRANC_GFOFFSET3D(u,0,0,-1),kmadd(ToReal(672),KRANC_GFOFFSET3D(u,0,0,1),kmadd(ToReal(168),KRANC_GFOFFSET3D(u,0,0,-2),kmadd(ToReal(-168),KRANC_GFOFFSET3D(u,0,0,2),kmadd(ToReal(-32),KRANC_GFOFFSET3D(u,0,0,-3),kmadd(ToReal(32),KRANC_GFOFFSET3D(u,0,0,3),kmadd(ToReal(3),KRANC_GFOFFSET3D(u,0,0,-4),kmul(ToReal(-3),KRANC_GFOFFSET3D(u,0,0,4)))))))))))
#else
#  define PDstandardfdOrder83(u) (PDstandardfdOrder83_impl(u,p1o840dz,cdj,cdk))
static CCTK_REAL_VEC PDstandardfdOrder83_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o840dz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardfdOrder83_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o840dz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardfdOrder82_impl(u, p1o840dz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder211(u) (kmul(p1odx2,kmadd(ToReal(-2),KRANC_GFOFFSET3D(u,0,0,0),kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)))))
#else
#  define PDstandardfdOrder211(u) (PDstandardfdOrder211_impl(u,p1odx2,cdj,cdk))
static CCTK_REAL_VEC PDstandardfdOrder211_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1odx2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardfdOrder211_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1odx2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1odx2,kmadd(ToReal(-2),KRANC_GFOFFSET3D(u,0,0,0),kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder222(u) (kmul(p1ody2,kmadd(ToReal(-2),KRANC_GFOFFSET3D(u,0,0,0),kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)))))
#else
#  define PDstandardfdOrder222(u) (PDstandardfdOrder222_impl(u,p1ody2,cdj,cdk))
static CCTK_REAL_VEC PDstandardfdOrder222_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1ody2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardfdOrder222_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1ody2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1ody2,kmadd(ToReal(-2),KRANC_GFOFFSET3D(u,0,0,0),kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder233(u) (kmul(p1odz2,kmadd(ToReal(-2),KRANC_GFOFFSET3D(u,0,0,0),kadd(KRANC_GFOFFSET3D(u,0,0,-1),KRANC_GFOFFSET3D(u,0,0,1)))))
#else
#  define PDstandardfdOrder233(u) (PDstandardfdOrder233_impl(u,p1odz2,cdj,cdk))
static CCTK_REAL_VEC PDstandardfdOrder233_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1odz2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardfdOrder233_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1odz2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardfdOrder222_impl(u, p1odz2, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder411(u) (kmul(pm1o12dx2,kmadd(ToReal(30),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(-16),kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),kadd(KRANC_GFOFFSET3D(u,-2,0,0),KRANC_GFOFFSET3D(u,2,0,0))))))
#else
#  define PDstandardfdOrder411(u) (PDstandardfdOrder411_impl(u,pm1o12dx2,cdj,cdk))
static CCTK_REAL_VEC PDstandardfdOrder411_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o12dx2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardfdOrder411_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o12dx2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(pm1o12dx2,kmadd(ToReal(30),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(-16),kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),kadd(KRANC_GFOFFSET3D(u,-2,0,0),KRANC_GFOFFSET3D(u,2,0,0)))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder422(u) (kmul(pm1o12dy2,kmadd(ToReal(30),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(-16),kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),kadd(KRANC_GFOFFSET3D(u,0,-2,0),KRANC_GFOFFSET3D(u,0,2,0))))))
#else
#  define PDstandardfdOrder422(u) (PDstandardfdOrder422_impl(u,pm1o12dy2,cdj,cdk))
static CCTK_REAL_VEC PDstandardfdOrder422_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o12dy2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardfdOrder422_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o12dy2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(pm1o12dy2,kmadd(ToReal(30),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(-16),kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),kadd(KRANC_GFOFFSET3D(u,0,-2,0),KRANC_GFOFFSET3D(u,0,2,0)))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder433(u) (kmul(pm1o12dz2,kmadd(ToReal(30),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(-16),kadd(KRANC_GFOFFSET3D(u,0,0,-1),KRANC_GFOFFSET3D(u,0,0,1)),kadd(KRANC_GFOFFSET3D(u,0,0,-2),KRANC_GFOFFSET3D(u,0,0,2))))))
#else
#  define PDstandardfdOrder433(u) (PDstandardfdOrder433_impl(u,pm1o12dz2,cdj,cdk))
static CCTK_REAL_VEC PDstandardfdOrder433_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o12dz2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardfdOrder433_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o12dz2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardfdOrder422_impl(u, pm1o12dz2, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder611(u) (kmul(p1o180dx2,kmadd(ToReal(-490),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(270),kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),kmadd(ToReal(-27),kadd(KRANC_GFOFFSET3D(u,-2,0,0),KRANC_GFOFFSET3D(u,2,0,0)),kmul(ToReal(2),kadd(KRANC_GFOFFSET3D(u,-3,0,0),KRANC_GFOFFSET3D(u,3,0,0))))))))
#else
#  define PDstandardfdOrder611(u) (PDstandardfdOrder611_impl(u,p1o180dx2,cdj,cdk))
static CCTK_REAL_VEC PDstandardfdOrder611_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o180dx2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardfdOrder611_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o180dx2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o180dx2,kmadd(ToReal(-490),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(270),kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),kmadd(ToReal(-27),kadd(KRANC_GFOFFSET3D(u,-2,0,0),KRANC_GFOFFSET3D(u,2,0,0)),kmul(ToReal(2),kadd(KRANC_GFOFFSET3D(u,-3,0,0),KRANC_GFOFFSET3D(u,3,0,0)))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder622(u) (kmul(p1o180dy2,kmadd(ToReal(-490),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(270),kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),kmadd(ToReal(-27),kadd(KRANC_GFOFFSET3D(u,0,-2,0),KRANC_GFOFFSET3D(u,0,2,0)),kmul(ToReal(2),kadd(KRANC_GFOFFSET3D(u,0,-3,0),KRANC_GFOFFSET3D(u,0,3,0))))))))
#else
#  define PDstandardfdOrder622(u) (PDstandardfdOrder622_impl(u,p1o180dy2,cdj,cdk))
static CCTK_REAL_VEC PDstandardfdOrder622_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o180dy2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardfdOrder622_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o180dy2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o180dy2,kmadd(ToReal(-490),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(270),kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),kmadd(ToReal(-27),kadd(KRANC_GFOFFSET3D(u,0,-2,0),KRANC_GFOFFSET3D(u,0,2,0)),kmul(ToReal(2),kadd(KRANC_GFOFFSET3D(u,0,-3,0),KRANC_GFOFFSET3D(u,0,3,0)))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder633(u) (kmul(p1o180dz2,kmadd(ToReal(-490),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(270),kadd(KRANC_GFOFFSET3D(u,0,0,-1),KRANC_GFOFFSET3D(u,0,0,1)),kmadd(ToReal(-27),kadd(KRANC_GFOFFSET3D(u,0,0,-2),KRANC_GFOFFSET3D(u,0,0,2)),kmul(ToReal(2),kadd(KRANC_GFOFFSET3D(u,0,0,-3),KRANC_GFOFFSET3D(u,0,0,3))))))))
#else
#  define PDstandardfdOrder633(u) (PDstandardfdOrder633_impl(u,p1o180dz2,cdj,cdk))
static CCTK_REAL_VEC PDstandardfdOrder633_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o180dz2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardfdOrder633_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o180dz2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardfdOrder622_impl(u, p1o180dz2, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder811(u) (kmul(p1o5040dx2,kmadd(ToReal(-14350),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(8064),kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),kmadd(ToReal(-1008),kadd(KRANC_GFOFFSET3D(u,-2,0,0),KRANC_GFOFFSET3D(u,2,0,0)),kmadd(ToReal(128),kadd(KRANC_GFOFFSET3D(u,-3,0,0),KRANC_GFOFFSET3D(u,3,0,0)),kmul(ToReal(-9),kadd(KRANC_GFOFFSET3D(u,-4,0,0),KRANC_GFOFFSET3D(u,4,0,0)))))))))
#else
#  define PDstandardfdOrder811(u) (PDstandardfdOrder811_impl(u,p1o5040dx2,cdj,cdk))
static CCTK_REAL_VEC PDstandardfdOrder811_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o5040dx2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardfdOrder811_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o5040dx2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o5040dx2,kmadd(ToReal(-14350),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(8064),kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),kmadd(ToReal(-1008),kadd(KRANC_GFOFFSET3D(u,-2,0,0),KRANC_GFOFFSET3D(u,2,0,0)),kmadd(ToReal(128),kadd(KRANC_GFOFFSET3D(u,-3,0,0),KRANC_GFOFFSET3D(u,3,0,0)),kmul(ToReal(-9),kadd(KRANC_GFOFFSET3D(u,-4,0,0),KRANC_GFOFFSET3D(u,4,0,0))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder822(u) (kmul(p1o5040dy2,kmadd(ToReal(-14350),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(8064),kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),kmadd(ToReal(-1008),kadd(KRANC_GFOFFSET3D(u,0,-2,0),KRANC_GFOFFSET3D(u,0,2,0)),kmadd(ToReal(128),kadd(KRANC_GFOFFSET3D(u,0,-3,0),KRANC_GFOFFSET3D(u,0,3,0)),kmul(ToReal(-9),kadd(KRANC_GFOFFSET3D(u,0,-4,0),KRANC_GFOFFSET3D(u,0,4,0)))))))))
#else
#  define PDstandardfdOrder822(u) (PDstandardfdOrder822_impl(u,p1o5040dy2,cdj,cdk))
static CCTK_REAL_VEC PDstandardfdOrder822_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o5040dy2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardfdOrder822_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o5040dy2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o5040dy2,kmadd(ToReal(-14350),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(8064),kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),kmadd(ToReal(-1008),kadd(KRANC_GFOFFSET3D(u,0,-2,0),KRANC_GFOFFSET3D(u,0,2,0)),kmadd(ToReal(128),kadd(KRANC_GFOFFSET3D(u,0,-3,0),KRANC_GFOFFSET3D(u,0,3,0)),kmul(ToReal(-9),kadd(KRANC_GFOFFSET3D(u,0,-4,0),KRANC_GFOFFSET3D(u,0,4,0))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder833(u) (kmul(p1o5040dz2,kmadd(ToReal(-14350),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(8064),kadd(KRANC_GFOFFSET3D(u,0,0,-1),KRANC_GFOFFSET3D(u,0,0,1)),kmadd(ToReal(-1008),kadd(KRANC_GFOFFSET3D(u,0,0,-2),KRANC_GFOFFSET3D(u,0,0,2)),kmadd(ToReal(128),kadd(KRANC_GFOFFSET3D(u,0,0,-3),KRANC_GFOFFSET3D(u,0,0,3)),kmul(ToReal(-9),kadd(KRANC_GFOFFSET3D(u,0,0,-4),KRANC_GFOFFSET3D(u,0,0,4)))))))))
#else
#  define PDstandardfdOrder833(u) (PDstandardfdOrder833_impl(u,p1o5040dz2,cdj,cdk))
static CCTK_REAL_VEC PDstandardfdOrder833_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o5040dz2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardfdOrder833_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o5040dz2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardfdOrder822_impl(u, p1o5040dz2, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder212(u) (kmul(p1o4dxdy,kadd(KRANC_GFOFFSET3D(u,-1,-1,0),ksub(ksub(KRANC_GFOFFSET3D(u,1,1,0),KRANC_GFOFFSET3D(u,1,-1,0)),KRANC_GFOFFSET3D(u,-1,1,0)))))
#else
#  define PDstandardfdOrder212(u) (PDstandardfdOrder212_impl(u,p1o4dxdy,cdj,cdk))
static CCTK_REAL_VEC PDstandardfdOrder212_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o4dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardfdOrder212_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o4dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o4dxdy,kadd(KRANC_GFOFFSET3D(u,-1,-1,0),ksub(ksub(KRANC_GFOFFSET3D(u,1,1,0),KRANC_GFOFFSET3D(u,1,-1,0)),KRANC_GFOFFSET3D(u,-1,1,0))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder213(u) (kmul(p1o4dxdz,kadd(KRANC_GFOFFSET3D(u,-1,0,-1),ksub(ksub(KRANC_GFOFFSET3D(u,1,0,1),KRANC_GFOFFSET3D(u,1,0,-1)),KRANC_GFOFFSET3D(u,-1,0,1)))))
#else
#  define PDstandardfdOrder213(u) (PDstandardfdOrder213_impl(u,p1o4dxdz,cdj,cdk))
static CCTK_REAL_VEC PDstandardfdOrder213_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o4dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardfdOrder213_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o4dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardfdOrder212_impl(u, p1o4dxdz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder221(u) (kmul(p1o4dxdy,kadd(KRANC_GFOFFSET3D(u,-1,-1,0),ksub(ksub(KRANC_GFOFFSET3D(u,1,1,0),KRANC_GFOFFSET3D(u,1,-1,0)),KRANC_GFOFFSET3D(u,-1,1,0)))))
#else
#  define PDstandardfdOrder221(u) (PDstandardfdOrder221_impl(u,p1o4dxdy,cdj,cdk))
static CCTK_REAL_VEC PDstandardfdOrder221_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o4dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardfdOrder221_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o4dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardfdOrder212_impl(u, p1o4dxdy, cdj, cdk);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder223(u) (kmul(p1o4dydz,kadd(KRANC_GFOFFSET3D(u,0,-1,-1),ksub(ksub(KRANC_GFOFFSET3D(u,0,1,1),KRANC_GFOFFSET3D(u,0,1,-1)),KRANC_GFOFFSET3D(u,0,-1,1)))))
#else
#  define PDstandardfdOrder223(u) (PDstandardfdOrder223_impl(u,p1o4dydz,cdj,cdk))
static CCTK_REAL_VEC PDstandardfdOrder223_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o4dydz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardfdOrder223_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o4dydz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o4dydz,kadd(KRANC_GFOFFSET3D(u,0,-1,-1),ksub(ksub(KRANC_GFOFFSET3D(u,0,1,1),KRANC_GFOFFSET3D(u,0,1,-1)),KRANC_GFOFFSET3D(u,0,-1,1))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder231(u) (kmul(p1o4dxdz,kadd(KRANC_GFOFFSET3D(u,-1,0,-1),ksub(ksub(KRANC_GFOFFSET3D(u,1,0,1),KRANC_GFOFFSET3D(u,1,0,-1)),KRANC_GFOFFSET3D(u,-1,0,1)))))
#else
#  define PDstandardfdOrder231(u) (PDstandardfdOrder231_impl(u,p1o4dxdz,cdj,cdk))
static CCTK_REAL_VEC PDstandardfdOrder231_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o4dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardfdOrder231_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o4dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardfdOrder212_impl(u, p1o4dxdz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder232(u) (kmul(p1o4dydz,kadd(KRANC_GFOFFSET3D(u,0,-1,-1),ksub(ksub(KRANC_GFOFFSET3D(u,0,1,1),KRANC_GFOFFSET3D(u,0,1,-1)),KRANC_GFOFFSET3D(u,0,-1,1)))))
#else
#  define PDstandardfdOrder232(u) (PDstandardfdOrder232_impl(u,p1o4dydz,cdj,cdk))
static CCTK_REAL_VEC PDstandardfdOrder232_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o4dydz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardfdOrder232_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o4dydz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardfdOrder223_impl(u, p1o4dydz, cdj, cdk);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder412(u) (kmul(p1o144dxdy,kmadd(ToReal(-64),kadd(KRANC_GFOFFSET3D(u,-1,1,0),KRANC_GFOFFSET3D(u,1,-1,0)),kmadd(ToReal(64),kadd(KRANC_GFOFFSET3D(u,-1,-1,0),KRANC_GFOFFSET3D(u,1,1,0)),kmadd(ToReal(8),kadd(KRANC_GFOFFSET3D(u,-1,2,0),kadd(KRANC_GFOFFSET3D(u,1,-2,0),kadd(KRANC_GFOFFSET3D(u,-2,1,0),KRANC_GFOFFSET3D(u,2,-1,0)))),kmadd(ToReal(-8),kadd(KRANC_GFOFFSET3D(u,-1,-2,0),kadd(KRANC_GFOFFSET3D(u,1,2,0),kadd(KRANC_GFOFFSET3D(u,-2,-1,0),KRANC_GFOFFSET3D(u,2,1,0)))),kadd(KRANC_GFOFFSET3D(u,-2,-2,0),ksub(ksub(KRANC_GFOFFSET3D(u,2,2,0),KRANC_GFOFFSET3D(u,2,-2,0)),KRANC_GFOFFSET3D(u,-2,2,0)))))))))
#else
#  define PDstandardfdOrder412(u) (PDstandardfdOrder412_impl(u,p1o144dxdy,cdj,cdk))
static CCTK_REAL_VEC PDstandardfdOrder412_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o144dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardfdOrder412_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o144dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o144dxdy,kmadd(ToReal(-64),kadd(KRANC_GFOFFSET3D(u,-1,1,0),KRANC_GFOFFSET3D(u,1,-1,0)),kmadd(ToReal(64),kadd(KRANC_GFOFFSET3D(u,-1,-1,0),KRANC_GFOFFSET3D(u,1,1,0)),kmadd(ToReal(8),kadd(KRANC_GFOFFSET3D(u,-1,2,0),kadd(KRANC_GFOFFSET3D(u,1,-2,0),kadd(KRANC_GFOFFSET3D(u,-2,1,0),KRANC_GFOFFSET3D(u,2,-1,0)))),kmadd(ToReal(-8),kadd(KRANC_GFOFFSET3D(u,-1,-2,0),kadd(KRANC_GFOFFSET3D(u,1,2,0),kadd(KRANC_GFOFFSET3D(u,-2,-1,0),KRANC_GFOFFSET3D(u,2,1,0)))),kadd(KRANC_GFOFFSET3D(u,-2,-2,0),ksub(ksub(KRANC_GFOFFSET3D(u,2,2,0),KRANC_GFOFFSET3D(u,2,-2,0)),KRANC_GFOFFSET3D(u,-2,2,0))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder413(u) (kmul(p1o144dxdz,kmadd(ToReal(-64),kadd(KRANC_GFOFFSET3D(u,-1,0,1),KRANC_GFOFFSET3D(u,1,0,-1)),kmadd(ToReal(64),kadd(KRANC_GFOFFSET3D(u,-1,0,-1),KRANC_GFOFFSET3D(u,1,0,1)),kmadd(ToReal(8),kadd(KRANC_GFOFFSET3D(u,-1,0,2),kadd(KRANC_GFOFFSET3D(u,1,0,-2),kadd(KRANC_GFOFFSET3D(u,-2,0,1),KRANC_GFOFFSET3D(u,2,0,-1)))),kmadd(ToReal(-8),kadd(KRANC_GFOFFSET3D(u,-1,0,-2),kadd(KRANC_GFOFFSET3D(u,1,0,2),kadd(KRANC_GFOFFSET3D(u,-2,0,-1),KRANC_GFOFFSET3D(u,2,0,1)))),kadd(KRANC_GFOFFSET3D(u,-2,0,-2),ksub(ksub(KRANC_GFOFFSET3D(u,2,0,2),KRANC_GFOFFSET3D(u,2,0,-2)),KRANC_GFOFFSET3D(u,-2,0,2)))))))))
#else
#  define PDstandardfdOrder413(u) (PDstandardfdOrder413_impl(u,p1o144dxdz,cdj,cdk))
static CCTK_REAL_VEC PDstandardfdOrder413_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o144dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardfdOrder413_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o144dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardfdOrder412_impl(u, p1o144dxdz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder421(u) (kmul(p1o144dxdy,kmadd(ToReal(-64),kadd(KRANC_GFOFFSET3D(u,-1,1,0),KRANC_GFOFFSET3D(u,1,-1,0)),kmadd(ToReal(64),kadd(KRANC_GFOFFSET3D(u,-1,-1,0),KRANC_GFOFFSET3D(u,1,1,0)),kmadd(ToReal(8),kadd(KRANC_GFOFFSET3D(u,-1,2,0),kadd(KRANC_GFOFFSET3D(u,1,-2,0),kadd(KRANC_GFOFFSET3D(u,-2,1,0),KRANC_GFOFFSET3D(u,2,-1,0)))),kmadd(ToReal(-8),kadd(KRANC_GFOFFSET3D(u,-1,-2,0),kadd(KRANC_GFOFFSET3D(u,1,2,0),kadd(KRANC_GFOFFSET3D(u,-2,-1,0),KRANC_GFOFFSET3D(u,2,1,0)))),kadd(KRANC_GFOFFSET3D(u,-2,-2,0),ksub(ksub(KRANC_GFOFFSET3D(u,2,2,0),KRANC_GFOFFSET3D(u,2,-2,0)),KRANC_GFOFFSET3D(u,-2,2,0)))))))))
#else
#  define PDstandardfdOrder421(u) (PDstandardfdOrder421_impl(u,p1o144dxdy,cdj,cdk))
static CCTK_REAL_VEC PDstandardfdOrder421_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o144dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardfdOrder421_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o144dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardfdOrder412_impl(u, p1o144dxdy, cdj, cdk);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder423(u) (kmul(p1o144dydz,kmadd(ToReal(-64),kadd(KRANC_GFOFFSET3D(u,0,-1,1),KRANC_GFOFFSET3D(u,0,1,-1)),kmadd(ToReal(64),kadd(KRANC_GFOFFSET3D(u,0,-1,-1),KRANC_GFOFFSET3D(u,0,1,1)),kmadd(ToReal(8),kadd(KRANC_GFOFFSET3D(u,0,-1,2),kadd(KRANC_GFOFFSET3D(u,0,1,-2),kadd(KRANC_GFOFFSET3D(u,0,-2,1),KRANC_GFOFFSET3D(u,0,2,-1)))),kmadd(ToReal(-8),kadd(KRANC_GFOFFSET3D(u,0,-1,-2),kadd(KRANC_GFOFFSET3D(u,0,1,2),kadd(KRANC_GFOFFSET3D(u,0,-2,-1),KRANC_GFOFFSET3D(u,0,2,1)))),kadd(KRANC_GFOFFSET3D(u,0,-2,-2),ksub(ksub(KRANC_GFOFFSET3D(u,0,2,2),KRANC_GFOFFSET3D(u,0,2,-2)),KRANC_GFOFFSET3D(u,0,-2,2)))))))))
#else
#  define PDstandardfdOrder423(u) (PDstandardfdOrder423_impl(u,p1o144dydz,cdj,cdk))
static CCTK_REAL_VEC PDstandardfdOrder423_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o144dydz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardfdOrder423_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o144dydz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o144dydz,kmadd(ToReal(-64),kadd(KRANC_GFOFFSET3D(u,0,-1,1),KRANC_GFOFFSET3D(u,0,1,-1)),kmadd(ToReal(64),kadd(KRANC_GFOFFSET3D(u,0,-1,-1),KRANC_GFOFFSET3D(u,0,1,1)),kmadd(ToReal(8),kadd(KRANC_GFOFFSET3D(u,0,-1,2),kadd(KRANC_GFOFFSET3D(u,0,1,-2),kadd(KRANC_GFOFFSET3D(u,0,-2,1),KRANC_GFOFFSET3D(u,0,2,-1)))),kmadd(ToReal(-8),kadd(KRANC_GFOFFSET3D(u,0,-1,-2),kadd(KRANC_GFOFFSET3D(u,0,1,2),kadd(KRANC_GFOFFSET3D(u,0,-2,-1),KRANC_GFOFFSET3D(u,0,2,1)))),kadd(KRANC_GFOFFSET3D(u,0,-2,-2),ksub(ksub(KRANC_GFOFFSET3D(u,0,2,2),KRANC_GFOFFSET3D(u,0,2,-2)),KRANC_GFOFFSET3D(u,0,-2,2))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder431(u) (kmul(p1o144dxdz,kmadd(ToReal(-64),kadd(KRANC_GFOFFSET3D(u,-1,0,1),KRANC_GFOFFSET3D(u,1,0,-1)),kmadd(ToReal(64),kadd(KRANC_GFOFFSET3D(u,-1,0,-1),KRANC_GFOFFSET3D(u,1,0,1)),kmadd(ToReal(8),kadd(KRANC_GFOFFSET3D(u,-1,0,2),kadd(KRANC_GFOFFSET3D(u,1,0,-2),kadd(KRANC_GFOFFSET3D(u,-2,0,1),KRANC_GFOFFSET3D(u,2,0,-1)))),kmadd(ToReal(-8),kadd(KRANC_GFOFFSET3D(u,-1,0,-2),kadd(KRANC_GFOFFSET3D(u,1,0,2),kadd(KRANC_GFOFFSET3D(u,-2,0,-1),KRANC_GFOFFSET3D(u,2,0,1)))),kadd(KRANC_GFOFFSET3D(u,-2,0,-2),ksub(ksub(KRANC_GFOFFSET3D(u,2,0,2),KRANC_GFOFFSET3D(u,2,0,-2)),KRANC_GFOFFSET3D(u,-2,0,2)))))))))
#else
#  define PDstandardfdOrder431(u) (PDstandardfdOrder431_impl(u,p1o144dxdz,cdj,cdk))
static CCTK_REAL_VEC PDstandardfdOrder431_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o144dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardfdOrder431_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o144dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardfdOrder412_impl(u, p1o144dxdz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder432(u) (kmul(p1o144dydz,kmadd(ToReal(-64),kadd(KRANC_GFOFFSET3D(u,0,-1,1),KRANC_GFOFFSET3D(u,0,1,-1)),kmadd(ToReal(64),kadd(KRANC_GFOFFSET3D(u,0,-1,-1),KRANC_GFOFFSET3D(u,0,1,1)),kmadd(ToReal(8),kadd(KRANC_GFOFFSET3D(u,0,-1,2),kadd(KRANC_GFOFFSET3D(u,0,1,-2),kadd(KRANC_GFOFFSET3D(u,0,-2,1),KRANC_GFOFFSET3D(u,0,2,-1)))),kmadd(ToReal(-8),kadd(KRANC_GFOFFSET3D(u,0,-1,-2),kadd(KRANC_GFOFFSET3D(u,0,1,2),kadd(KRANC_GFOFFSET3D(u,0,-2,-1),KRANC_GFOFFSET3D(u,0,2,1)))),kadd(KRANC_GFOFFSET3D(u,0,-2,-2),ksub(ksub(KRANC_GFOFFSET3D(u,0,2,2),KRANC_GFOFFSET3D(u,0,2,-2)),KRANC_GFOFFSET3D(u,0,-2,2)))))))))
#else
#  define PDstandardfdOrder432(u) (PDstandardfdOrder432_impl(u,p1o144dydz,cdj,cdk))
static CCTK_REAL_VEC PDstandardfdOrder432_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o144dydz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardfdOrder432_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o144dydz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardfdOrder423_impl(u, p1o144dydz, cdj, cdk);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder612(u) (kmul(p1o3600dxdy,kmadd(ToReal(-2025),kadd(KRANC_GFOFFSET3D(u,-1,1,0),KRANC_GFOFFSET3D(u,1,-1,0)),kmadd(ToReal(2025),kadd(KRANC_GFOFFSET3D(u,-1,-1,0),KRANC_GFOFFSET3D(u,1,1,0)),kmadd(ToReal(405),kadd(KRANC_GFOFFSET3D(u,-1,2,0),kadd(KRANC_GFOFFSET3D(u,1,-2,0),kadd(KRANC_GFOFFSET3D(u,-2,1,0),KRANC_GFOFFSET3D(u,2,-1,0)))),kmadd(ToReal(-405),kadd(KRANC_GFOFFSET3D(u,-1,-2,0),kadd(KRANC_GFOFFSET3D(u,1,2,0),kadd(KRANC_GFOFFSET3D(u,-2,-1,0),KRANC_GFOFFSET3D(u,2,1,0)))),kmadd(ToReal(-81),kadd(KRANC_GFOFFSET3D(u,-2,2,0),KRANC_GFOFFSET3D(u,2,-2,0)),kmadd(ToReal(81),kadd(KRANC_GFOFFSET3D(u,-2,-2,0),KRANC_GFOFFSET3D(u,2,2,0)),kmadd(ToReal(-45),kadd(KRANC_GFOFFSET3D(u,-1,3,0),kadd(KRANC_GFOFFSET3D(u,1,-3,0),kadd(KRANC_GFOFFSET3D(u,-3,1,0),KRANC_GFOFFSET3D(u,3,-1,0)))),kmadd(ToReal(45),kadd(KRANC_GFOFFSET3D(u,-1,-3,0),kadd(KRANC_GFOFFSET3D(u,1,3,0),kadd(KRANC_GFOFFSET3D(u,-3,-1,0),KRANC_GFOFFSET3D(u,3,1,0)))),kmadd(ToReal(9),kadd(KRANC_GFOFFSET3D(u,-2,3,0),kadd(KRANC_GFOFFSET3D(u,2,-3,0),kadd(KRANC_GFOFFSET3D(u,-3,2,0),KRANC_GFOFFSET3D(u,3,-2,0)))),kmadd(ToReal(-9),kadd(KRANC_GFOFFSET3D(u,-2,-3,0),kadd(KRANC_GFOFFSET3D(u,2,3,0),kadd(KRANC_GFOFFSET3D(u,-3,-2,0),KRANC_GFOFFSET3D(u,3,2,0)))),kadd(KRANC_GFOFFSET3D(u,-3,-3,0),ksub(ksub(KRANC_GFOFFSET3D(u,3,3,0),KRANC_GFOFFSET3D(u,3,-3,0)),KRANC_GFOFFSET3D(u,-3,3,0)))))))))))))))
#else
#  define PDstandardfdOrder612(u) (PDstandardfdOrder612_impl(u,p1o3600dxdy,cdj,cdk))
static CCTK_REAL_VEC PDstandardfdOrder612_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o3600dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardfdOrder612_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o3600dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o3600dxdy,kmadd(ToReal(-2025),kadd(KRANC_GFOFFSET3D(u,-1,1,0),KRANC_GFOFFSET3D(u,1,-1,0)),kmadd(ToReal(2025),kadd(KRANC_GFOFFSET3D(u,-1,-1,0),KRANC_GFOFFSET3D(u,1,1,0)),kmadd(ToReal(405),kadd(KRANC_GFOFFSET3D(u,-1,2,0),kadd(KRANC_GFOFFSET3D(u,1,-2,0),kadd(KRANC_GFOFFSET3D(u,-2,1,0),KRANC_GFOFFSET3D(u,2,-1,0)))),kmadd(ToReal(-405),kadd(KRANC_GFOFFSET3D(u,-1,-2,0),kadd(KRANC_GFOFFSET3D(u,1,2,0),kadd(KRANC_GFOFFSET3D(u,-2,-1,0),KRANC_GFOFFSET3D(u,2,1,0)))),kmadd(ToReal(-81),kadd(KRANC_GFOFFSET3D(u,-2,2,0),KRANC_GFOFFSET3D(u,2,-2,0)),kmadd(ToReal(81),kadd(KRANC_GFOFFSET3D(u,-2,-2,0),KRANC_GFOFFSET3D(u,2,2,0)),kmadd(ToReal(-45),kadd(KRANC_GFOFFSET3D(u,-1,3,0),kadd(KRANC_GFOFFSET3D(u,1,-3,0),kadd(KRANC_GFOFFSET3D(u,-3,1,0),KRANC_GFOFFSET3D(u,3,-1,0)))),kmadd(ToReal(45),kadd(KRANC_GFOFFSET3D(u,-1,-3,0),kadd(KRANC_GFOFFSET3D(u,1,3,0),kadd(KRANC_GFOFFSET3D(u,-3,-1,0),KRANC_GFOFFSET3D(u,3,1,0)))),kmadd(ToReal(9),kadd(KRANC_GFOFFSET3D(u,-2,3,0),kadd(KRANC_GFOFFSET3D(u,2,-3,0),kadd(KRANC_GFOFFSET3D(u,-3,2,0),KRANC_GFOFFSET3D(u,3,-2,0)))),kmadd(ToReal(-9),kadd(KRANC_GFOFFSET3D(u,-2,-3,0),kadd(KRANC_GFOFFSET3D(u,2,3,0),kadd(KRANC_GFOFFSET3D(u,-3,-2,0),KRANC_GFOFFSET3D(u,3,2,0)))),kadd(KRANC_GFOFFSET3D(u,-3,-3,0),ksub(ksub(KRANC_GFOFFSET3D(u,3,3,0),KRANC_GFOFFSET3D(u,3,-3,0)),KRANC_GFOFFSET3D(u,-3,3,0))))))))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder613(u) (kmul(p1o3600dxdz,kmadd(ToReal(-2025),kadd(KRANC_GFOFFSET3D(u,-1,0,1),KRANC_GFOFFSET3D(u,1,0,-1)),kmadd(ToReal(2025),kadd(KRANC_GFOFFSET3D(u,-1,0,-1),KRANC_GFOFFSET3D(u,1,0,1)),kmadd(ToReal(405),kadd(KRANC_GFOFFSET3D(u,-1,0,2),kadd(KRANC_GFOFFSET3D(u,1,0,-2),kadd(KRANC_GFOFFSET3D(u,-2,0,1),KRANC_GFOFFSET3D(u,2,0,-1)))),kmadd(ToReal(-405),kadd(KRANC_GFOFFSET3D(u,-1,0,-2),kadd(KRANC_GFOFFSET3D(u,1,0,2),kadd(KRANC_GFOFFSET3D(u,-2,0,-1),KRANC_GFOFFSET3D(u,2,0,1)))),kmadd(ToReal(-81),kadd(KRANC_GFOFFSET3D(u,-2,0,2),KRANC_GFOFFSET3D(u,2,0,-2)),kmadd(ToReal(81),kadd(KRANC_GFOFFSET3D(u,-2,0,-2),KRANC_GFOFFSET3D(u,2,0,2)),kmadd(ToReal(-45),kadd(KRANC_GFOFFSET3D(u,-1,0,3),kadd(KRANC_GFOFFSET3D(u,1,0,-3),kadd(KRANC_GFOFFSET3D(u,-3,0,1),KRANC_GFOFFSET3D(u,3,0,-1)))),kmadd(ToReal(45),kadd(KRANC_GFOFFSET3D(u,-1,0,-3),kadd(KRANC_GFOFFSET3D(u,1,0,3),kadd(KRANC_GFOFFSET3D(u,-3,0,-1),KRANC_GFOFFSET3D(u,3,0,1)))),kmadd(ToReal(9),kadd(KRANC_GFOFFSET3D(u,-2,0,3),kadd(KRANC_GFOFFSET3D(u,2,0,-3),kadd(KRANC_GFOFFSET3D(u,-3,0,2),KRANC_GFOFFSET3D(u,3,0,-2)))),kmadd(ToReal(-9),kadd(KRANC_GFOFFSET3D(u,-2,0,-3),kadd(KRANC_GFOFFSET3D(u,2,0,3),kadd(KRANC_GFOFFSET3D(u,-3,0,-2),KRANC_GFOFFSET3D(u,3,0,2)))),kadd(KRANC_GFOFFSET3D(u,-3,0,-3),ksub(ksub(KRANC_GFOFFSET3D(u,3,0,3),KRANC_GFOFFSET3D(u,3,0,-3)),KRANC_GFOFFSET3D(u,-3,0,3)))))))))))))))
#else
#  define PDstandardfdOrder613(u) (PDstandardfdOrder613_impl(u,p1o3600dxdz,cdj,cdk))
static CCTK_REAL_VEC PDstandardfdOrder613_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o3600dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardfdOrder613_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o3600dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardfdOrder612_impl(u, p1o3600dxdz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder621(u) (kmul(p1o3600dxdy,kmadd(ToReal(-2025),kadd(KRANC_GFOFFSET3D(u,-1,1,0),KRANC_GFOFFSET3D(u,1,-1,0)),kmadd(ToReal(2025),kadd(KRANC_GFOFFSET3D(u,-1,-1,0),KRANC_GFOFFSET3D(u,1,1,0)),kmadd(ToReal(405),kadd(KRANC_GFOFFSET3D(u,-1,2,0),kadd(KRANC_GFOFFSET3D(u,1,-2,0),kadd(KRANC_GFOFFSET3D(u,-2,1,0),KRANC_GFOFFSET3D(u,2,-1,0)))),kmadd(ToReal(-405),kadd(KRANC_GFOFFSET3D(u,-1,-2,0),kadd(KRANC_GFOFFSET3D(u,1,2,0),kadd(KRANC_GFOFFSET3D(u,-2,-1,0),KRANC_GFOFFSET3D(u,2,1,0)))),kmadd(ToReal(-81),kadd(KRANC_GFOFFSET3D(u,-2,2,0),KRANC_GFOFFSET3D(u,2,-2,0)),kmadd(ToReal(81),kadd(KRANC_GFOFFSET3D(u,-2,-2,0),KRANC_GFOFFSET3D(u,2,2,0)),kmadd(ToReal(-45),kadd(KRANC_GFOFFSET3D(u,-1,3,0),kadd(KRANC_GFOFFSET3D(u,1,-3,0),kadd(KRANC_GFOFFSET3D(u,-3,1,0),KRANC_GFOFFSET3D(u,3,-1,0)))),kmadd(ToReal(45),kadd(KRANC_GFOFFSET3D(u,-1,-3,0),kadd(KRANC_GFOFFSET3D(u,1,3,0),kadd(KRANC_GFOFFSET3D(u,-3,-1,0),KRANC_GFOFFSET3D(u,3,1,0)))),kmadd(ToReal(9),kadd(KRANC_GFOFFSET3D(u,-2,3,0),kadd(KRANC_GFOFFSET3D(u,2,-3,0),kadd(KRANC_GFOFFSET3D(u,-3,2,0),KRANC_GFOFFSET3D(u,3,-2,0)))),kmadd(ToReal(-9),kadd(KRANC_GFOFFSET3D(u,-2,-3,0),kadd(KRANC_GFOFFSET3D(u,2,3,0),kadd(KRANC_GFOFFSET3D(u,-3,-2,0),KRANC_GFOFFSET3D(u,3,2,0)))),kadd(KRANC_GFOFFSET3D(u,-3,-3,0),ksub(ksub(KRANC_GFOFFSET3D(u,3,3,0),KRANC_GFOFFSET3D(u,3,-3,0)),KRANC_GFOFFSET3D(u,-3,3,0)))))))))))))))
#else
#  define PDstandardfdOrder621(u) (PDstandardfdOrder621_impl(u,p1o3600dxdy,cdj,cdk))
static CCTK_REAL_VEC PDstandardfdOrder621_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o3600dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardfdOrder621_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o3600dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardfdOrder612_impl(u, p1o3600dxdy, cdj, cdk);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder623(u) (kmul(p1o3600dydz,kmadd(ToReal(-2025),kadd(KRANC_GFOFFSET3D(u,0,-1,1),KRANC_GFOFFSET3D(u,0,1,-1)),kmadd(ToReal(2025),kadd(KRANC_GFOFFSET3D(u,0,-1,-1),KRANC_GFOFFSET3D(u,0,1,1)),kmadd(ToReal(405),kadd(KRANC_GFOFFSET3D(u,0,-1,2),kadd(KRANC_GFOFFSET3D(u,0,1,-2),kadd(KRANC_GFOFFSET3D(u,0,-2,1),KRANC_GFOFFSET3D(u,0,2,-1)))),kmadd(ToReal(-405),kadd(KRANC_GFOFFSET3D(u,0,-1,-2),kadd(KRANC_GFOFFSET3D(u,0,1,2),kadd(KRANC_GFOFFSET3D(u,0,-2,-1),KRANC_GFOFFSET3D(u,0,2,1)))),kmadd(ToReal(-81),kadd(KRANC_GFOFFSET3D(u,0,-2,2),KRANC_GFOFFSET3D(u,0,2,-2)),kmadd(ToReal(81),kadd(KRANC_GFOFFSET3D(u,0,-2,-2),KRANC_GFOFFSET3D(u,0,2,2)),kmadd(ToReal(-45),kadd(KRANC_GFOFFSET3D(u,0,-1,3),kadd(KRANC_GFOFFSET3D(u,0,1,-3),kadd(KRANC_GFOFFSET3D(u,0,-3,1),KRANC_GFOFFSET3D(u,0,3,-1)))),kmadd(ToReal(45),kadd(KRANC_GFOFFSET3D(u,0,-1,-3),kadd(KRANC_GFOFFSET3D(u,0,1,3),kadd(KRANC_GFOFFSET3D(u,0,-3,-1),KRANC_GFOFFSET3D(u,0,3,1)))),kmadd(ToReal(9),kadd(KRANC_GFOFFSET3D(u,0,-2,3),kadd(KRANC_GFOFFSET3D(u,0,2,-3),kadd(KRANC_GFOFFSET3D(u,0,-3,2),KRANC_GFOFFSET3D(u,0,3,-2)))),kmadd(ToReal(-9),kadd(KRANC_GFOFFSET3D(u,0,-2,-3),kadd(KRANC_GFOFFSET3D(u,0,2,3),kadd(KRANC_GFOFFSET3D(u,0,-3,-2),KRANC_GFOFFSET3D(u,0,3,2)))),kadd(KRANC_GFOFFSET3D(u,0,-3,-3),ksub(ksub(KRANC_GFOFFSET3D(u,0,3,3),KRANC_GFOFFSET3D(u,0,3,-3)),KRANC_GFOFFSET3D(u,0,-3,3)))))))))))))))
#else
#  define PDstandardfdOrder623(u) (PDstandardfdOrder623_impl(u,p1o3600dydz,cdj,cdk))
static CCTK_REAL_VEC PDstandardfdOrder623_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o3600dydz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardfdOrder623_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o3600dydz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o3600dydz,kmadd(ToReal(-2025),kadd(KRANC_GFOFFSET3D(u,0,-1,1),KRANC_GFOFFSET3D(u,0,1,-1)),kmadd(ToReal(2025),kadd(KRANC_GFOFFSET3D(u,0,-1,-1),KRANC_GFOFFSET3D(u,0,1,1)),kmadd(ToReal(405),kadd(KRANC_GFOFFSET3D(u,0,-1,2),kadd(KRANC_GFOFFSET3D(u,0,1,-2),kadd(KRANC_GFOFFSET3D(u,0,-2,1),KRANC_GFOFFSET3D(u,0,2,-1)))),kmadd(ToReal(-405),kadd(KRANC_GFOFFSET3D(u,0,-1,-2),kadd(KRANC_GFOFFSET3D(u,0,1,2),kadd(KRANC_GFOFFSET3D(u,0,-2,-1),KRANC_GFOFFSET3D(u,0,2,1)))),kmadd(ToReal(-81),kadd(KRANC_GFOFFSET3D(u,0,-2,2),KRANC_GFOFFSET3D(u,0,2,-2)),kmadd(ToReal(81),kadd(KRANC_GFOFFSET3D(u,0,-2,-2),KRANC_GFOFFSET3D(u,0,2,2)),kmadd(ToReal(-45),kadd(KRANC_GFOFFSET3D(u,0,-1,3),kadd(KRANC_GFOFFSET3D(u,0,1,-3),kadd(KRANC_GFOFFSET3D(u,0,-3,1),KRANC_GFOFFSET3D(u,0,3,-1)))),kmadd(ToReal(45),kadd(KRANC_GFOFFSET3D(u,0,-1,-3),kadd(KRANC_GFOFFSET3D(u,0,1,3),kadd(KRANC_GFOFFSET3D(u,0,-3,-1),KRANC_GFOFFSET3D(u,0,3,1)))),kmadd(ToReal(9),kadd(KRANC_GFOFFSET3D(u,0,-2,3),kadd(KRANC_GFOFFSET3D(u,0,2,-3),kadd(KRANC_GFOFFSET3D(u,0,-3,2),KRANC_GFOFFSET3D(u,0,3,-2)))),kmadd(ToReal(-9),kadd(KRANC_GFOFFSET3D(u,0,-2,-3),kadd(KRANC_GFOFFSET3D(u,0,2,3),kadd(KRANC_GFOFFSET3D(u,0,-3,-2),KRANC_GFOFFSET3D(u,0,3,2)))),kadd(KRANC_GFOFFSET3D(u,0,-3,-3),ksub(ksub(KRANC_GFOFFSET3D(u,0,3,3),KRANC_GFOFFSET3D(u,0,3,-3)),KRANC_GFOFFSET3D(u,0,-3,3))))))))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder631(u) (kmul(p1o3600dxdz,kmadd(ToReal(-2025),kadd(KRANC_GFOFFSET3D(u,-1,0,1),KRANC_GFOFFSET3D(u,1,0,-1)),kmadd(ToReal(2025),kadd(KRANC_GFOFFSET3D(u,-1,0,-1),KRANC_GFOFFSET3D(u,1,0,1)),kmadd(ToReal(405),kadd(KRANC_GFOFFSET3D(u,-1,0,2),kadd(KRANC_GFOFFSET3D(u,1,0,-2),kadd(KRANC_GFOFFSET3D(u,-2,0,1),KRANC_GFOFFSET3D(u,2,0,-1)))),kmadd(ToReal(-405),kadd(KRANC_GFOFFSET3D(u,-1,0,-2),kadd(KRANC_GFOFFSET3D(u,1,0,2),kadd(KRANC_GFOFFSET3D(u,-2,0,-1),KRANC_GFOFFSET3D(u,2,0,1)))),kmadd(ToReal(-81),kadd(KRANC_GFOFFSET3D(u,-2,0,2),KRANC_GFOFFSET3D(u,2,0,-2)),kmadd(ToReal(81),kadd(KRANC_GFOFFSET3D(u,-2,0,-2),KRANC_GFOFFSET3D(u,2,0,2)),kmadd(ToReal(-45),kadd(KRANC_GFOFFSET3D(u,-1,0,3),kadd(KRANC_GFOFFSET3D(u,1,0,-3),kadd(KRANC_GFOFFSET3D(u,-3,0,1),KRANC_GFOFFSET3D(u,3,0,-1)))),kmadd(ToReal(45),kadd(KRANC_GFOFFSET3D(u,-1,0,-3),kadd(KRANC_GFOFFSET3D(u,1,0,3),kadd(KRANC_GFOFFSET3D(u,-3,0,-1),KRANC_GFOFFSET3D(u,3,0,1)))),kmadd(ToReal(9),kadd(KRANC_GFOFFSET3D(u,-2,0,3),kadd(KRANC_GFOFFSET3D(u,2,0,-3),kadd(KRANC_GFOFFSET3D(u,-3,0,2),KRANC_GFOFFSET3D(u,3,0,-2)))),kmadd(ToReal(-9),kadd(KRANC_GFOFFSET3D(u,-2,0,-3),kadd(KRANC_GFOFFSET3D(u,2,0,3),kadd(KRANC_GFOFFSET3D(u,-3,0,-2),KRANC_GFOFFSET3D(u,3,0,2)))),kadd(KRANC_GFOFFSET3D(u,-3,0,-3),ksub(ksub(KRANC_GFOFFSET3D(u,3,0,3),KRANC_GFOFFSET3D(u,3,0,-3)),KRANC_GFOFFSET3D(u,-3,0,3)))))))))))))))
#else
#  define PDstandardfdOrder631(u) (PDstandardfdOrder631_impl(u,p1o3600dxdz,cdj,cdk))
static CCTK_REAL_VEC PDstandardfdOrder631_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o3600dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardfdOrder631_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o3600dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardfdOrder612_impl(u, p1o3600dxdz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder632(u) (kmul(p1o3600dydz,kmadd(ToReal(-2025),kadd(KRANC_GFOFFSET3D(u,0,-1,1),KRANC_GFOFFSET3D(u,0,1,-1)),kmadd(ToReal(2025),kadd(KRANC_GFOFFSET3D(u,0,-1,-1),KRANC_GFOFFSET3D(u,0,1,1)),kmadd(ToReal(405),kadd(KRANC_GFOFFSET3D(u,0,-1,2),kadd(KRANC_GFOFFSET3D(u,0,1,-2),kadd(KRANC_GFOFFSET3D(u,0,-2,1),KRANC_GFOFFSET3D(u,0,2,-1)))),kmadd(ToReal(-405),kadd(KRANC_GFOFFSET3D(u,0,-1,-2),kadd(KRANC_GFOFFSET3D(u,0,1,2),kadd(KRANC_GFOFFSET3D(u,0,-2,-1),KRANC_GFOFFSET3D(u,0,2,1)))),kmadd(ToReal(-81),kadd(KRANC_GFOFFSET3D(u,0,-2,2),KRANC_GFOFFSET3D(u,0,2,-2)),kmadd(ToReal(81),kadd(KRANC_GFOFFSET3D(u,0,-2,-2),KRANC_GFOFFSET3D(u,0,2,2)),kmadd(ToReal(-45),kadd(KRANC_GFOFFSET3D(u,0,-1,3),kadd(KRANC_GFOFFSET3D(u,0,1,-3),kadd(KRANC_GFOFFSET3D(u,0,-3,1),KRANC_GFOFFSET3D(u,0,3,-1)))),kmadd(ToReal(45),kadd(KRANC_GFOFFSET3D(u,0,-1,-3),kadd(KRANC_GFOFFSET3D(u,0,1,3),kadd(KRANC_GFOFFSET3D(u,0,-3,-1),KRANC_GFOFFSET3D(u,0,3,1)))),kmadd(ToReal(9),kadd(KRANC_GFOFFSET3D(u,0,-2,3),kadd(KRANC_GFOFFSET3D(u,0,2,-3),kadd(KRANC_GFOFFSET3D(u,0,-3,2),KRANC_GFOFFSET3D(u,0,3,-2)))),kmadd(ToReal(-9),kadd(KRANC_GFOFFSET3D(u,0,-2,-3),kadd(KRANC_GFOFFSET3D(u,0,2,3),kadd(KRANC_GFOFFSET3D(u,0,-3,-2),KRANC_GFOFFSET3D(u,0,3,2)))),kadd(KRANC_GFOFFSET3D(u,0,-3,-3),ksub(ksub(KRANC_GFOFFSET3D(u,0,3,3),KRANC_GFOFFSET3D(u,0,3,-3)),KRANC_GFOFFSET3D(u,0,-3,3)))))))))))))))
#else
#  define PDstandardfdOrder632(u) (PDstandardfdOrder632_impl(u,p1o3600dydz,cdj,cdk))
static CCTK_REAL_VEC PDstandardfdOrder632_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o3600dydz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardfdOrder632_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o3600dydz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardfdOrder623_impl(u, p1o3600dydz, cdj, cdk);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder812(u) (kmul(p1o705600dxdy,kmadd(ToReal(-451584),kadd(KRANC_GFOFFSET3D(u,-1,1,0),KRANC_GFOFFSET3D(u,1,-1,0)),kmadd(ToReal(451584),kadd(KRANC_GFOFFSET3D(u,-1,-1,0),KRANC_GFOFFSET3D(u,1,1,0)),kmadd(ToReal(112896),kadd(KRANC_GFOFFSET3D(u,-1,2,0),kadd(KRANC_GFOFFSET3D(u,1,-2,0),kadd(KRANC_GFOFFSET3D(u,-2,1,0),KRANC_GFOFFSET3D(u,2,-1,0)))),kmadd(ToReal(-112896),kadd(KRANC_GFOFFSET3D(u,-1,-2,0),kadd(KRANC_GFOFFSET3D(u,1,2,0),kadd(KRANC_GFOFFSET3D(u,-2,-1,0),KRANC_GFOFFSET3D(u,2,1,0)))),kmadd(ToReal(-28224),kadd(KRANC_GFOFFSET3D(u,-2,2,0),KRANC_GFOFFSET3D(u,2,-2,0)),kmadd(ToReal(28224),kadd(KRANC_GFOFFSET3D(u,-2,-2,0),KRANC_GFOFFSET3D(u,2,2,0)),kmadd(ToReal(-21504),kadd(KRANC_GFOFFSET3D(u,-1,3,0),kadd(KRANC_GFOFFSET3D(u,1,-3,0),kadd(KRANC_GFOFFSET3D(u,-3,1,0),KRANC_GFOFFSET3D(u,3,-1,0)))),kmadd(ToReal(21504),kadd(KRANC_GFOFFSET3D(u,-1,-3,0),kadd(KRANC_GFOFFSET3D(u,1,3,0),kadd(KRANC_GFOFFSET3D(u,-3,-1,0),KRANC_GFOFFSET3D(u,3,1,0)))),kmadd(ToReal(5376),kadd(KRANC_GFOFFSET3D(u,-2,3,0),kadd(KRANC_GFOFFSET3D(u,2,-3,0),kadd(KRANC_GFOFFSET3D(u,-3,2,0),KRANC_GFOFFSET3D(u,3,-2,0)))),kmadd(ToReal(-5376),kadd(KRANC_GFOFFSET3D(u,-2,-3,0),kadd(KRANC_GFOFFSET3D(u,2,3,0),kadd(KRANC_GFOFFSET3D(u,-3,-2,0),KRANC_GFOFFSET3D(u,3,2,0)))),kmadd(ToReal(-1024),kadd(KRANC_GFOFFSET3D(u,-3,3,0),KRANC_GFOFFSET3D(u,3,-3,0)),kmadd(ToReal(1024),kadd(KRANC_GFOFFSET3D(u,-3,-3,0),KRANC_GFOFFSET3D(u,3,3,0)),kmadd(ToReal(2016),kadd(KRANC_GFOFFSET3D(u,-1,4,0),kadd(KRANC_GFOFFSET3D(u,1,-4,0),kadd(KRANC_GFOFFSET3D(u,-4,1,0),KRANC_GFOFFSET3D(u,4,-1,0)))),kmadd(ToReal(-2016),kadd(KRANC_GFOFFSET3D(u,-1,-4,0),kadd(KRANC_GFOFFSET3D(u,1,4,0),kadd(KRANC_GFOFFSET3D(u,-4,-1,0),KRANC_GFOFFSET3D(u,4,1,0)))),kmadd(ToReal(-504),kadd(KRANC_GFOFFSET3D(u,-2,4,0),kadd(KRANC_GFOFFSET3D(u,2,-4,0),kadd(KRANC_GFOFFSET3D(u,-4,2,0),KRANC_GFOFFSET3D(u,4,-2,0)))),kmadd(ToReal(504),kadd(KRANC_GFOFFSET3D(u,-2,-4,0),kadd(KRANC_GFOFFSET3D(u,2,4,0),kadd(KRANC_GFOFFSET3D(u,-4,-2,0),KRANC_GFOFFSET3D(u,4,2,0)))),kmadd(ToReal(96),kadd(KRANC_GFOFFSET3D(u,-3,4,0),kadd(KRANC_GFOFFSET3D(u,3,-4,0),kadd(KRANC_GFOFFSET3D(u,-4,3,0),KRANC_GFOFFSET3D(u,4,-3,0)))),kmadd(ToReal(-96),kadd(KRANC_GFOFFSET3D(u,-3,-4,0),kadd(KRANC_GFOFFSET3D(u,3,4,0),kadd(KRANC_GFOFFSET3D(u,-4,-3,0),KRANC_GFOFFSET3D(u,4,3,0)))),kmadd(ToReal(-9),kadd(KRANC_GFOFFSET3D(u,-4,4,0),KRANC_GFOFFSET3D(u,4,-4,0)),kmul(ToReal(9),kadd(KRANC_GFOFFSET3D(u,-4,-4,0),KRANC_GFOFFSET3D(u,4,4,0))))))))))))))))))))))))
#else
#  define PDstandardfdOrder812(u) (PDstandardfdOrder812_impl(u,p1o705600dxdy,cdj,cdk))
static CCTK_REAL_VEC PDstandardfdOrder812_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o705600dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardfdOrder812_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o705600dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o705600dxdy,kmadd(ToReal(-451584),kadd(KRANC_GFOFFSET3D(u,-1,1,0),KRANC_GFOFFSET3D(u,1,-1,0)),kmadd(ToReal(451584),kadd(KRANC_GFOFFSET3D(u,-1,-1,0),KRANC_GFOFFSET3D(u,1,1,0)),kmadd(ToReal(112896),kadd(KRANC_GFOFFSET3D(u,-1,2,0),kadd(KRANC_GFOFFSET3D(u,1,-2,0),kadd(KRANC_GFOFFSET3D(u,-2,1,0),KRANC_GFOFFSET3D(u,2,-1,0)))),kmadd(ToReal(-112896),kadd(KRANC_GFOFFSET3D(u,-1,-2,0),kadd(KRANC_GFOFFSET3D(u,1,2,0),kadd(KRANC_GFOFFSET3D(u,-2,-1,0),KRANC_GFOFFSET3D(u,2,1,0)))),kmadd(ToReal(-28224),kadd(KRANC_GFOFFSET3D(u,-2,2,0),KRANC_GFOFFSET3D(u,2,-2,0)),kmadd(ToReal(28224),kadd(KRANC_GFOFFSET3D(u,-2,-2,0),KRANC_GFOFFSET3D(u,2,2,0)),kmadd(ToReal(-21504),kadd(KRANC_GFOFFSET3D(u,-1,3,0),kadd(KRANC_GFOFFSET3D(u,1,-3,0),kadd(KRANC_GFOFFSET3D(u,-3,1,0),KRANC_GFOFFSET3D(u,3,-1,0)))),kmadd(ToReal(21504),kadd(KRANC_GFOFFSET3D(u,-1,-3,0),kadd(KRANC_GFOFFSET3D(u,1,3,0),kadd(KRANC_GFOFFSET3D(u,-3,-1,0),KRANC_GFOFFSET3D(u,3,1,0)))),kmadd(ToReal(5376),kadd(KRANC_GFOFFSET3D(u,-2,3,0),kadd(KRANC_GFOFFSET3D(u,2,-3,0),kadd(KRANC_GFOFFSET3D(u,-3,2,0),KRANC_GFOFFSET3D(u,3,-2,0)))),kmadd(ToReal(-5376),kadd(KRANC_GFOFFSET3D(u,-2,-3,0),kadd(KRANC_GFOFFSET3D(u,2,3,0),kadd(KRANC_GFOFFSET3D(u,-3,-2,0),KRANC_GFOFFSET3D(u,3,2,0)))),kmadd(ToReal(-1024),kadd(KRANC_GFOFFSET3D(u,-3,3,0),KRANC_GFOFFSET3D(u,3,-3,0)),kmadd(ToReal(1024),kadd(KRANC_GFOFFSET3D(u,-3,-3,0),KRANC_GFOFFSET3D(u,3,3,0)),kmadd(ToReal(2016),kadd(KRANC_GFOFFSET3D(u,-1,4,0),kadd(KRANC_GFOFFSET3D(u,1,-4,0),kadd(KRANC_GFOFFSET3D(u,-4,1,0),KRANC_GFOFFSET3D(u,4,-1,0)))),kmadd(ToReal(-2016),kadd(KRANC_GFOFFSET3D(u,-1,-4,0),kadd(KRANC_GFOFFSET3D(u,1,4,0),kadd(KRANC_GFOFFSET3D(u,-4,-1,0),KRANC_GFOFFSET3D(u,4,1,0)))),kmadd(ToReal(-504),kadd(KRANC_GFOFFSET3D(u,-2,4,0),kadd(KRANC_GFOFFSET3D(u,2,-4,0),kadd(KRANC_GFOFFSET3D(u,-4,2,0),KRANC_GFOFFSET3D(u,4,-2,0)))),kmadd(ToReal(504),kadd(KRANC_GFOFFSET3D(u,-2,-4,0),kadd(KRANC_GFOFFSET3D(u,2,4,0),kadd(KRANC_GFOFFSET3D(u,-4,-2,0),KRANC_GFOFFSET3D(u,4,2,0)))),kmadd(ToReal(96),kadd(KRANC_GFOFFSET3D(u,-3,4,0),kadd(KRANC_GFOFFSET3D(u,3,-4,0),kadd(KRANC_GFOFFSET3D(u,-4,3,0),KRANC_GFOFFSET3D(u,4,-3,0)))),kmadd(ToReal(-96),kadd(KRANC_GFOFFSET3D(u,-3,-4,0),kadd(KRANC_GFOFFSET3D(u,3,4,0),kadd(KRANC_GFOFFSET3D(u,-4,-3,0),KRANC_GFOFFSET3D(u,4,3,0)))),kmadd(ToReal(-9),kadd(KRANC_GFOFFSET3D(u,-4,4,0),KRANC_GFOFFSET3D(u,4,-4,0)),kmul(ToReal(9),kadd(KRANC_GFOFFSET3D(u,-4,-4,0),KRANC_GFOFFSET3D(u,4,4,0)))))))))))))))))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder813(u) (kmul(p1o705600dxdz,kmadd(ToReal(-451584),kadd(KRANC_GFOFFSET3D(u,-1,0,1),KRANC_GFOFFSET3D(u,1,0,-1)),kmadd(ToReal(451584),kadd(KRANC_GFOFFSET3D(u,-1,0,-1),KRANC_GFOFFSET3D(u,1,0,1)),kmadd(ToReal(112896),kadd(KRANC_GFOFFSET3D(u,-1,0,2),kadd(KRANC_GFOFFSET3D(u,1,0,-2),kadd(KRANC_GFOFFSET3D(u,-2,0,1),KRANC_GFOFFSET3D(u,2,0,-1)))),kmadd(ToReal(-112896),kadd(KRANC_GFOFFSET3D(u,-1,0,-2),kadd(KRANC_GFOFFSET3D(u,1,0,2),kadd(KRANC_GFOFFSET3D(u,-2,0,-1),KRANC_GFOFFSET3D(u,2,0,1)))),kmadd(ToReal(-28224),kadd(KRANC_GFOFFSET3D(u,-2,0,2),KRANC_GFOFFSET3D(u,2,0,-2)),kmadd(ToReal(28224),kadd(KRANC_GFOFFSET3D(u,-2,0,-2),KRANC_GFOFFSET3D(u,2,0,2)),kmadd(ToReal(-21504),kadd(KRANC_GFOFFSET3D(u,-1,0,3),kadd(KRANC_GFOFFSET3D(u,1,0,-3),kadd(KRANC_GFOFFSET3D(u,-3,0,1),KRANC_GFOFFSET3D(u,3,0,-1)))),kmadd(ToReal(21504),kadd(KRANC_GFOFFSET3D(u,-1,0,-3),kadd(KRANC_GFOFFSET3D(u,1,0,3),kadd(KRANC_GFOFFSET3D(u,-3,0,-1),KRANC_GFOFFSET3D(u,3,0,1)))),kmadd(ToReal(5376),kadd(KRANC_GFOFFSET3D(u,-2,0,3),kadd(KRANC_GFOFFSET3D(u,2,0,-3),kadd(KRANC_GFOFFSET3D(u,-3,0,2),KRANC_GFOFFSET3D(u,3,0,-2)))),kmadd(ToReal(-5376),kadd(KRANC_GFOFFSET3D(u,-2,0,-3),kadd(KRANC_GFOFFSET3D(u,2,0,3),kadd(KRANC_GFOFFSET3D(u,-3,0,-2),KRANC_GFOFFSET3D(u,3,0,2)))),kmadd(ToReal(-1024),kadd(KRANC_GFOFFSET3D(u,-3,0,3),KRANC_GFOFFSET3D(u,3,0,-3)),kmadd(ToReal(1024),kadd(KRANC_GFOFFSET3D(u,-3,0,-3),KRANC_GFOFFSET3D(u,3,0,3)),kmadd(ToReal(2016),kadd(KRANC_GFOFFSET3D(u,-1,0,4),kadd(KRANC_GFOFFSET3D(u,1,0,-4),kadd(KRANC_GFOFFSET3D(u,-4,0,1),KRANC_GFOFFSET3D(u,4,0,-1)))),kmadd(ToReal(-2016),kadd(KRANC_GFOFFSET3D(u,-1,0,-4),kadd(KRANC_GFOFFSET3D(u,1,0,4),kadd(KRANC_GFOFFSET3D(u,-4,0,-1),KRANC_GFOFFSET3D(u,4,0,1)))),kmadd(ToReal(-504),kadd(KRANC_GFOFFSET3D(u,-2,0,4),kadd(KRANC_GFOFFSET3D(u,2,0,-4),kadd(KRANC_GFOFFSET3D(u,-4,0,2),KRANC_GFOFFSET3D(u,4,0,-2)))),kmadd(ToReal(504),kadd(KRANC_GFOFFSET3D(u,-2,0,-4),kadd(KRANC_GFOFFSET3D(u,2,0,4),kadd(KRANC_GFOFFSET3D(u,-4,0,-2),KRANC_GFOFFSET3D(u,4,0,2)))),kmadd(ToReal(96),kadd(KRANC_GFOFFSET3D(u,-3,0,4),kadd(KRANC_GFOFFSET3D(u,3,0,-4),kadd(KRANC_GFOFFSET3D(u,-4,0,3),KRANC_GFOFFSET3D(u,4,0,-3)))),kmadd(ToReal(-96),kadd(KRANC_GFOFFSET3D(u,-3,0,-4),kadd(KRANC_GFOFFSET3D(u,3,0,4),kadd(KRANC_GFOFFSET3D(u,-4,0,-3),KRANC_GFOFFSET3D(u,4,0,3)))),kmadd(ToReal(-9),kadd(KRANC_GFOFFSET3D(u,-4,0,4),KRANC_GFOFFSET3D(u,4,0,-4)),kmul(ToReal(9),kadd(KRANC_GFOFFSET3D(u,-4,0,-4),KRANC_GFOFFSET3D(u,4,0,4))))))))))))))))))))))))
#else
#  define PDstandardfdOrder813(u) (PDstandardfdOrder813_impl(u,p1o705600dxdz,cdj,cdk))
static CCTK_REAL_VEC PDstandardfdOrder813_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o705600dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardfdOrder813_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o705600dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardfdOrder812_impl(u, p1o705600dxdz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder821(u) (kmul(p1o705600dxdy,kmadd(ToReal(-451584),kadd(KRANC_GFOFFSET3D(u,-1,1,0),KRANC_GFOFFSET3D(u,1,-1,0)),kmadd(ToReal(451584),kadd(KRANC_GFOFFSET3D(u,-1,-1,0),KRANC_GFOFFSET3D(u,1,1,0)),kmadd(ToReal(112896),kadd(KRANC_GFOFFSET3D(u,-1,2,0),kadd(KRANC_GFOFFSET3D(u,1,-2,0),kadd(KRANC_GFOFFSET3D(u,-2,1,0),KRANC_GFOFFSET3D(u,2,-1,0)))),kmadd(ToReal(-112896),kadd(KRANC_GFOFFSET3D(u,-1,-2,0),kadd(KRANC_GFOFFSET3D(u,1,2,0),kadd(KRANC_GFOFFSET3D(u,-2,-1,0),KRANC_GFOFFSET3D(u,2,1,0)))),kmadd(ToReal(-28224),kadd(KRANC_GFOFFSET3D(u,-2,2,0),KRANC_GFOFFSET3D(u,2,-2,0)),kmadd(ToReal(28224),kadd(KRANC_GFOFFSET3D(u,-2,-2,0),KRANC_GFOFFSET3D(u,2,2,0)),kmadd(ToReal(-21504),kadd(KRANC_GFOFFSET3D(u,-1,3,0),kadd(KRANC_GFOFFSET3D(u,1,-3,0),kadd(KRANC_GFOFFSET3D(u,-3,1,0),KRANC_GFOFFSET3D(u,3,-1,0)))),kmadd(ToReal(21504),kadd(KRANC_GFOFFSET3D(u,-1,-3,0),kadd(KRANC_GFOFFSET3D(u,1,3,0),kadd(KRANC_GFOFFSET3D(u,-3,-1,0),KRANC_GFOFFSET3D(u,3,1,0)))),kmadd(ToReal(5376),kadd(KRANC_GFOFFSET3D(u,-2,3,0),kadd(KRANC_GFOFFSET3D(u,2,-3,0),kadd(KRANC_GFOFFSET3D(u,-3,2,0),KRANC_GFOFFSET3D(u,3,-2,0)))),kmadd(ToReal(-5376),kadd(KRANC_GFOFFSET3D(u,-2,-3,0),kadd(KRANC_GFOFFSET3D(u,2,3,0),kadd(KRANC_GFOFFSET3D(u,-3,-2,0),KRANC_GFOFFSET3D(u,3,2,0)))),kmadd(ToReal(-1024),kadd(KRANC_GFOFFSET3D(u,-3,3,0),KRANC_GFOFFSET3D(u,3,-3,0)),kmadd(ToReal(1024),kadd(KRANC_GFOFFSET3D(u,-3,-3,0),KRANC_GFOFFSET3D(u,3,3,0)),kmadd(ToReal(2016),kadd(KRANC_GFOFFSET3D(u,-1,4,0),kadd(KRANC_GFOFFSET3D(u,1,-4,0),kadd(KRANC_GFOFFSET3D(u,-4,1,0),KRANC_GFOFFSET3D(u,4,-1,0)))),kmadd(ToReal(-2016),kadd(KRANC_GFOFFSET3D(u,-1,-4,0),kadd(KRANC_GFOFFSET3D(u,1,4,0),kadd(KRANC_GFOFFSET3D(u,-4,-1,0),KRANC_GFOFFSET3D(u,4,1,0)))),kmadd(ToReal(-504),kadd(KRANC_GFOFFSET3D(u,-2,4,0),kadd(KRANC_GFOFFSET3D(u,2,-4,0),kadd(KRANC_GFOFFSET3D(u,-4,2,0),KRANC_GFOFFSET3D(u,4,-2,0)))),kmadd(ToReal(504),kadd(KRANC_GFOFFSET3D(u,-2,-4,0),kadd(KRANC_GFOFFSET3D(u,2,4,0),kadd(KRANC_GFOFFSET3D(u,-4,-2,0),KRANC_GFOFFSET3D(u,4,2,0)))),kmadd(ToReal(96),kadd(KRANC_GFOFFSET3D(u,-3,4,0),kadd(KRANC_GFOFFSET3D(u,3,-4,0),kadd(KRANC_GFOFFSET3D(u,-4,3,0),KRANC_GFOFFSET3D(u,4,-3,0)))),kmadd(ToReal(-96),kadd(KRANC_GFOFFSET3D(u,-3,-4,0),kadd(KRANC_GFOFFSET3D(u,3,4,0),kadd(KRANC_GFOFFSET3D(u,-4,-3,0),KRANC_GFOFFSET3D(u,4,3,0)))),kmadd(ToReal(-9),kadd(KRANC_GFOFFSET3D(u,-4,4,0),KRANC_GFOFFSET3D(u,4,-4,0)),kmul(ToReal(9),kadd(KRANC_GFOFFSET3D(u,-4,-4,0),KRANC_GFOFFSET3D(u,4,4,0))))))))))))))))))))))))
#else
#  define PDstandardfdOrder821(u) (PDstandardfdOrder821_impl(u,p1o705600dxdy,cdj,cdk))
static CCTK_REAL_VEC PDstandardfdOrder821_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o705600dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardfdOrder821_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o705600dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardfdOrder812_impl(u, p1o705600dxdy, cdj, cdk);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder823(u) (kmul(p1o705600dydz,kmadd(ToReal(-451584),kadd(KRANC_GFOFFSET3D(u,0,-1,1),KRANC_GFOFFSET3D(u,0,1,-1)),kmadd(ToReal(451584),kadd(KRANC_GFOFFSET3D(u,0,-1,-1),KRANC_GFOFFSET3D(u,0,1,1)),kmadd(ToReal(112896),kadd(KRANC_GFOFFSET3D(u,0,-1,2),kadd(KRANC_GFOFFSET3D(u,0,1,-2),kadd(KRANC_GFOFFSET3D(u,0,-2,1),KRANC_GFOFFSET3D(u,0,2,-1)))),kmadd(ToReal(-112896),kadd(KRANC_GFOFFSET3D(u,0,-1,-2),kadd(KRANC_GFOFFSET3D(u,0,1,2),kadd(KRANC_GFOFFSET3D(u,0,-2,-1),KRANC_GFOFFSET3D(u,0,2,1)))),kmadd(ToReal(-28224),kadd(KRANC_GFOFFSET3D(u,0,-2,2),KRANC_GFOFFSET3D(u,0,2,-2)),kmadd(ToReal(28224),kadd(KRANC_GFOFFSET3D(u,0,-2,-2),KRANC_GFOFFSET3D(u,0,2,2)),kmadd(ToReal(-21504),kadd(KRANC_GFOFFSET3D(u,0,-1,3),kadd(KRANC_GFOFFSET3D(u,0,1,-3),kadd(KRANC_GFOFFSET3D(u,0,-3,1),KRANC_GFOFFSET3D(u,0,3,-1)))),kmadd(ToReal(21504),kadd(KRANC_GFOFFSET3D(u,0,-1,-3),kadd(KRANC_GFOFFSET3D(u,0,1,3),kadd(KRANC_GFOFFSET3D(u,0,-3,-1),KRANC_GFOFFSET3D(u,0,3,1)))),kmadd(ToReal(5376),kadd(KRANC_GFOFFSET3D(u,0,-2,3),kadd(KRANC_GFOFFSET3D(u,0,2,-3),kadd(KRANC_GFOFFSET3D(u,0,-3,2),KRANC_GFOFFSET3D(u,0,3,-2)))),kmadd(ToReal(-5376),kadd(KRANC_GFOFFSET3D(u,0,-2,-3),kadd(KRANC_GFOFFSET3D(u,0,2,3),kadd(KRANC_GFOFFSET3D(u,0,-3,-2),KRANC_GFOFFSET3D(u,0,3,2)))),kmadd(ToReal(-1024),kadd(KRANC_GFOFFSET3D(u,0,-3,3),KRANC_GFOFFSET3D(u,0,3,-3)),kmadd(ToReal(1024),kadd(KRANC_GFOFFSET3D(u,0,-3,-3),KRANC_GFOFFSET3D(u,0,3,3)),kmadd(ToReal(2016),kadd(KRANC_GFOFFSET3D(u,0,-1,4),kadd(KRANC_GFOFFSET3D(u,0,1,-4),kadd(KRANC_GFOFFSET3D(u,0,-4,1),KRANC_GFOFFSET3D(u,0,4,-1)))),kmadd(ToReal(-2016),kadd(KRANC_GFOFFSET3D(u,0,-1,-4),kadd(KRANC_GFOFFSET3D(u,0,1,4),kadd(KRANC_GFOFFSET3D(u,0,-4,-1),KRANC_GFOFFSET3D(u,0,4,1)))),kmadd(ToReal(-504),kadd(KRANC_GFOFFSET3D(u,0,-2,4),kadd(KRANC_GFOFFSET3D(u,0,2,-4),kadd(KRANC_GFOFFSET3D(u,0,-4,2),KRANC_GFOFFSET3D(u,0,4,-2)))),kmadd(ToReal(504),kadd(KRANC_GFOFFSET3D(u,0,-2,-4),kadd(KRANC_GFOFFSET3D(u,0,2,4),kadd(KRANC_GFOFFSET3D(u,0,-4,-2),KRANC_GFOFFSET3D(u,0,4,2)))),kmadd(ToReal(96),kadd(KRANC_GFOFFSET3D(u,0,-3,4),kadd(KRANC_GFOFFSET3D(u,0,3,-4),kadd(KRANC_GFOFFSET3D(u,0,-4,3),KRANC_GFOFFSET3D(u,0,4,-3)))),kmadd(ToReal(-96),kadd(KRANC_GFOFFSET3D(u,0,-3,-4),kadd(KRANC_GFOFFSET3D(u,0,3,4),kadd(KRANC_GFOFFSET3D(u,0,-4,-3),KRANC_GFOFFSET3D(u,0,4,3)))),kmadd(ToReal(-9),kadd(KRANC_GFOFFSET3D(u,0,-4,4),KRANC_GFOFFSET3D(u,0,4,-4)),kmul(ToReal(9),kadd(KRANC_GFOFFSET3D(u,0,-4,-4),KRANC_GFOFFSET3D(u,0,4,4))))))))))))))))))))))))
#else
#  define PDstandardfdOrder823(u) (PDstandardfdOrder823_impl(u,p1o705600dydz,cdj,cdk))
static CCTK_REAL_VEC PDstandardfdOrder823_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o705600dydz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardfdOrder823_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o705600dydz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o705600dydz,kmadd(ToReal(-451584),kadd(KRANC_GFOFFSET3D(u,0,-1,1),KRANC_GFOFFSET3D(u,0,1,-1)),kmadd(ToReal(451584),kadd(KRANC_GFOFFSET3D(u,0,-1,-1),KRANC_GFOFFSET3D(u,0,1,1)),kmadd(ToReal(112896),kadd(KRANC_GFOFFSET3D(u,0,-1,2),kadd(KRANC_GFOFFSET3D(u,0,1,-2),kadd(KRANC_GFOFFSET3D(u,0,-2,1),KRANC_GFOFFSET3D(u,0,2,-1)))),kmadd(ToReal(-112896),kadd(KRANC_GFOFFSET3D(u,0,-1,-2),kadd(KRANC_GFOFFSET3D(u,0,1,2),kadd(KRANC_GFOFFSET3D(u,0,-2,-1),KRANC_GFOFFSET3D(u,0,2,1)))),kmadd(ToReal(-28224),kadd(KRANC_GFOFFSET3D(u,0,-2,2),KRANC_GFOFFSET3D(u,0,2,-2)),kmadd(ToReal(28224),kadd(KRANC_GFOFFSET3D(u,0,-2,-2),KRANC_GFOFFSET3D(u,0,2,2)),kmadd(ToReal(-21504),kadd(KRANC_GFOFFSET3D(u,0,-1,3),kadd(KRANC_GFOFFSET3D(u,0,1,-3),kadd(KRANC_GFOFFSET3D(u,0,-3,1),KRANC_GFOFFSET3D(u,0,3,-1)))),kmadd(ToReal(21504),kadd(KRANC_GFOFFSET3D(u,0,-1,-3),kadd(KRANC_GFOFFSET3D(u,0,1,3),kadd(KRANC_GFOFFSET3D(u,0,-3,-1),KRANC_GFOFFSET3D(u,0,3,1)))),kmadd(ToReal(5376),kadd(KRANC_GFOFFSET3D(u,0,-2,3),kadd(KRANC_GFOFFSET3D(u,0,2,-3),kadd(KRANC_GFOFFSET3D(u,0,-3,2),KRANC_GFOFFSET3D(u,0,3,-2)))),kmadd(ToReal(-5376),kadd(KRANC_GFOFFSET3D(u,0,-2,-3),kadd(KRANC_GFOFFSET3D(u,0,2,3),kadd(KRANC_GFOFFSET3D(u,0,-3,-2),KRANC_GFOFFSET3D(u,0,3,2)))),kmadd(ToReal(-1024),kadd(KRANC_GFOFFSET3D(u,0,-3,3),KRANC_GFOFFSET3D(u,0,3,-3)),kmadd(ToReal(1024),kadd(KRANC_GFOFFSET3D(u,0,-3,-3),KRANC_GFOFFSET3D(u,0,3,3)),kmadd(ToReal(2016),kadd(KRANC_GFOFFSET3D(u,0,-1,4),kadd(KRANC_GFOFFSET3D(u,0,1,-4),kadd(KRANC_GFOFFSET3D(u,0,-4,1),KRANC_GFOFFSET3D(u,0,4,-1)))),kmadd(ToReal(-2016),kadd(KRANC_GFOFFSET3D(u,0,-1,-4),kadd(KRANC_GFOFFSET3D(u,0,1,4),kadd(KRANC_GFOFFSET3D(u,0,-4,-1),KRANC_GFOFFSET3D(u,0,4,1)))),kmadd(ToReal(-504),kadd(KRANC_GFOFFSET3D(u,0,-2,4),kadd(KRANC_GFOFFSET3D(u,0,2,-4),kadd(KRANC_GFOFFSET3D(u,0,-4,2),KRANC_GFOFFSET3D(u,0,4,-2)))),kmadd(ToReal(504),kadd(KRANC_GFOFFSET3D(u,0,-2,-4),kadd(KRANC_GFOFFSET3D(u,0,2,4),kadd(KRANC_GFOFFSET3D(u,0,-4,-2),KRANC_GFOFFSET3D(u,0,4,2)))),kmadd(ToReal(96),kadd(KRANC_GFOFFSET3D(u,0,-3,4),kadd(KRANC_GFOFFSET3D(u,0,3,-4),kadd(KRANC_GFOFFSET3D(u,0,-4,3),KRANC_GFOFFSET3D(u,0,4,-3)))),kmadd(ToReal(-96),kadd(KRANC_GFOFFSET3D(u,0,-3,-4),kadd(KRANC_GFOFFSET3D(u,0,3,4),kadd(KRANC_GFOFFSET3D(u,0,-4,-3),KRANC_GFOFFSET3D(u,0,4,3)))),kmadd(ToReal(-9),kadd(KRANC_GFOFFSET3D(u,0,-4,4),KRANC_GFOFFSET3D(u,0,4,-4)),kmul(ToReal(9),kadd(KRANC_GFOFFSET3D(u,0,-4,-4),KRANC_GFOFFSET3D(u,0,4,4)))))))))))))))))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder831(u) (kmul(p1o705600dxdz,kmadd(ToReal(-451584),kadd(KRANC_GFOFFSET3D(u,-1,0,1),KRANC_GFOFFSET3D(u,1,0,-1)),kmadd(ToReal(451584),kadd(KRANC_GFOFFSET3D(u,-1,0,-1),KRANC_GFOFFSET3D(u,1,0,1)),kmadd(ToReal(112896),kadd(KRANC_GFOFFSET3D(u,-1,0,2),kadd(KRANC_GFOFFSET3D(u,1,0,-2),kadd(KRANC_GFOFFSET3D(u,-2,0,1),KRANC_GFOFFSET3D(u,2,0,-1)))),kmadd(ToReal(-112896),kadd(KRANC_GFOFFSET3D(u,-1,0,-2),kadd(KRANC_GFOFFSET3D(u,1,0,2),kadd(KRANC_GFOFFSET3D(u,-2,0,-1),KRANC_GFOFFSET3D(u,2,0,1)))),kmadd(ToReal(-28224),kadd(KRANC_GFOFFSET3D(u,-2,0,2),KRANC_GFOFFSET3D(u,2,0,-2)),kmadd(ToReal(28224),kadd(KRANC_GFOFFSET3D(u,-2,0,-2),KRANC_GFOFFSET3D(u,2,0,2)),kmadd(ToReal(-21504),kadd(KRANC_GFOFFSET3D(u,-1,0,3),kadd(KRANC_GFOFFSET3D(u,1,0,-3),kadd(KRANC_GFOFFSET3D(u,-3,0,1),KRANC_GFOFFSET3D(u,3,0,-1)))),kmadd(ToReal(21504),kadd(KRANC_GFOFFSET3D(u,-1,0,-3),kadd(KRANC_GFOFFSET3D(u,1,0,3),kadd(KRANC_GFOFFSET3D(u,-3,0,-1),KRANC_GFOFFSET3D(u,3,0,1)))),kmadd(ToReal(5376),kadd(KRANC_GFOFFSET3D(u,-2,0,3),kadd(KRANC_GFOFFSET3D(u,2,0,-3),kadd(KRANC_GFOFFSET3D(u,-3,0,2),KRANC_GFOFFSET3D(u,3,0,-2)))),kmadd(ToReal(-5376),kadd(KRANC_GFOFFSET3D(u,-2,0,-3),kadd(KRANC_GFOFFSET3D(u,2,0,3),kadd(KRANC_GFOFFSET3D(u,-3,0,-2),KRANC_GFOFFSET3D(u,3,0,2)))),kmadd(ToReal(-1024),kadd(KRANC_GFOFFSET3D(u,-3,0,3),KRANC_GFOFFSET3D(u,3,0,-3)),kmadd(ToReal(1024),kadd(KRANC_GFOFFSET3D(u,-3,0,-3),KRANC_GFOFFSET3D(u,3,0,3)),kmadd(ToReal(2016),kadd(KRANC_GFOFFSET3D(u,-1,0,4),kadd(KRANC_GFOFFSET3D(u,1,0,-4),kadd(KRANC_GFOFFSET3D(u,-4,0,1),KRANC_GFOFFSET3D(u,4,0,-1)))),kmadd(ToReal(-2016),kadd(KRANC_GFOFFSET3D(u,-1,0,-4),kadd(KRANC_GFOFFSET3D(u,1,0,4),kadd(KRANC_GFOFFSET3D(u,-4,0,-1),KRANC_GFOFFSET3D(u,4,0,1)))),kmadd(ToReal(-504),kadd(KRANC_GFOFFSET3D(u,-2,0,4),kadd(KRANC_GFOFFSET3D(u,2,0,-4),kadd(KRANC_GFOFFSET3D(u,-4,0,2),KRANC_GFOFFSET3D(u,4,0,-2)))),kmadd(ToReal(504),kadd(KRANC_GFOFFSET3D(u,-2,0,-4),kadd(KRANC_GFOFFSET3D(u,2,0,4),kadd(KRANC_GFOFFSET3D(u,-4,0,-2),KRANC_GFOFFSET3D(u,4,0,2)))),kmadd(ToReal(96),kadd(KRANC_GFOFFSET3D(u,-3,0,4),kadd(KRANC_GFOFFSET3D(u,3,0,-4),kadd(KRANC_GFOFFSET3D(u,-4,0,3),KRANC_GFOFFSET3D(u,4,0,-3)))),kmadd(ToReal(-96),kadd(KRANC_GFOFFSET3D(u,-3,0,-4),kadd(KRANC_GFOFFSET3D(u,3,0,4),kadd(KRANC_GFOFFSET3D(u,-4,0,-3),KRANC_GFOFFSET3D(u,4,0,3)))),kmadd(ToReal(-9),kadd(KRANC_GFOFFSET3D(u,-4,0,4),KRANC_GFOFFSET3D(u,4,0,-4)),kmul(ToReal(9),kadd(KRANC_GFOFFSET3D(u,-4,0,-4),KRANC_GFOFFSET3D(u,4,0,4))))))))))))))))))))))))
#else
#  define PDstandardfdOrder831(u) (PDstandardfdOrder831_impl(u,p1o705600dxdz,cdj,cdk))
static CCTK_REAL_VEC PDstandardfdOrder831_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o705600dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardfdOrder831_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o705600dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardfdOrder812_impl(u, p1o705600dxdz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardfdOrder832(u) (kmul(p1o705600dydz,kmadd(ToReal(-451584),kadd(KRANC_GFOFFSET3D(u,0,-1,1),KRANC_GFOFFSET3D(u,0,1,-1)),kmadd(ToReal(451584),kadd(KRANC_GFOFFSET3D(u,0,-1,-1),KRANC_GFOFFSET3D(u,0,1,1)),kmadd(ToReal(112896),kadd(KRANC_GFOFFSET3D(u,0,-1,2),kadd(KRANC_GFOFFSET3D(u,0,1,-2),kadd(KRANC_GFOFFSET3D(u,0,-2,1),KRANC_GFOFFSET3D(u,0,2,-1)))),kmadd(ToReal(-112896),kadd(KRANC_GFOFFSET3D(u,0,-1,-2),kadd(KRANC_GFOFFSET3D(u,0,1,2),kadd(KRANC_GFOFFSET3D(u,0,-2,-1),KRANC_GFOFFSET3D(u,0,2,1)))),kmadd(ToReal(-28224),kadd(KRANC_GFOFFSET3D(u,0,-2,2),KRANC_GFOFFSET3D(u,0,2,-2)),kmadd(ToReal(28224),kadd(KRANC_GFOFFSET3D(u,0,-2,-2),KRANC_GFOFFSET3D(u,0,2,2)),kmadd(ToReal(-21504),kadd(KRANC_GFOFFSET3D(u,0,-1,3),kadd(KRANC_GFOFFSET3D(u,0,1,-3),kadd(KRANC_GFOFFSET3D(u,0,-3,1),KRANC_GFOFFSET3D(u,0,3,-1)))),kmadd(ToReal(21504),kadd(KRANC_GFOFFSET3D(u,0,-1,-3),kadd(KRANC_GFOFFSET3D(u,0,1,3),kadd(KRANC_GFOFFSET3D(u,0,-3,-1),KRANC_GFOFFSET3D(u,0,3,1)))),kmadd(ToReal(5376),kadd(KRANC_GFOFFSET3D(u,0,-2,3),kadd(KRANC_GFOFFSET3D(u,0,2,-3),kadd(KRANC_GFOFFSET3D(u,0,-3,2),KRANC_GFOFFSET3D(u,0,3,-2)))),kmadd(ToReal(-5376),kadd(KRANC_GFOFFSET3D(u,0,-2,-3),kadd(KRANC_GFOFFSET3D(u,0,2,3),kadd(KRANC_GFOFFSET3D(u,0,-3,-2),KRANC_GFOFFSET3D(u,0,3,2)))),kmadd(ToReal(-1024),kadd(KRANC_GFOFFSET3D(u,0,-3,3),KRANC_GFOFFSET3D(u,0,3,-3)),kmadd(ToReal(1024),kadd(KRANC_GFOFFSET3D(u,0,-3,-3),KRANC_GFOFFSET3D(u,0,3,3)),kmadd(ToReal(2016),kadd(KRANC_GFOFFSET3D(u,0,-1,4),kadd(KRANC_GFOFFSET3D(u,0,1,-4),kadd(KRANC_GFOFFSET3D(u,0,-4,1),KRANC_GFOFFSET3D(u,0,4,-1)))),kmadd(ToReal(-2016),kadd(KRANC_GFOFFSET3D(u,0,-1,-4),kadd(KRANC_GFOFFSET3D(u,0,1,4),kadd(KRANC_GFOFFSET3D(u,0,-4,-1),KRANC_GFOFFSET3D(u,0,4,1)))),kmadd(ToReal(-504),kadd(KRANC_GFOFFSET3D(u,0,-2,4),kadd(KRANC_GFOFFSET3D(u,0,2,-4),kadd(KRANC_GFOFFSET3D(u,0,-4,2),KRANC_GFOFFSET3D(u,0,4,-2)))),kmadd(ToReal(504),kadd(KRANC_GFOFFSET3D(u,0,-2,-4),kadd(KRANC_GFOFFSET3D(u,0,2,4),kadd(KRANC_GFOFFSET3D(u,0,-4,-2),KRANC_GFOFFSET3D(u,0,4,2)))),kmadd(ToReal(96),kadd(KRANC_GFOFFSET3D(u,0,-3,4),kadd(KRANC_GFOFFSET3D(u,0,3,-4),kadd(KRANC_GFOFFSET3D(u,0,-4,3),KRANC_GFOFFSET3D(u,0,4,-3)))),kmadd(ToReal(-96),kadd(KRANC_GFOFFSET3D(u,0,-3,-4),kadd(KRANC_GFOFFSET3D(u,0,3,4),kadd(KRANC_GFOFFSET3D(u,0,-4,-3),KRANC_GFOFFSET3D(u,0,4,3)))),kmadd(ToReal(-9),kadd(KRANC_GFOFFSET3D(u,0,-4,4),KRANC_GFOFFSET3D(u,0,4,-4)),kmul(ToReal(9),kadd(KRANC_GFOFFSET3D(u,0,-4,-4),KRANC_GFOFFSET3D(u,0,4,4))))))))))))))))))))))))
#else
#  define PDstandardfdOrder832(u) (PDstandardfdOrder832_impl(u,p1o705600dydz,cdj,cdk))
static CCTK_REAL_VEC PDstandardfdOrder832_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o705600dydz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardfdOrder832_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o705600dydz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardfdOrder823_impl(u, p1o705600dydz, cdj, cdk);
}
#endif

