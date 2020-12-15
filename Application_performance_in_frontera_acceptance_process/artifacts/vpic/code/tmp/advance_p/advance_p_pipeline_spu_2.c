#ifdef CELL_SPU_BUILD

#define IN_particle_pipeline
#define HAS_SPU_PIPELINE
#include <particle_pipelines.h>
#include <v4c_spu.h>
#include <spu_mfcio.h>

#ifdef IN_HARNESS
#include <profile.h>
#endif

// DMA tag usage:
//  0: 2 - Particle buffer 0 (read, write, mover write)
//  3: 5 - Particle buffer 1 (read, write, mover write)
//  6: 8 - Particle buffer 2 (read, write, mover write)
//  9:16 - Interpolator cache
// 17:25 - Accumulator cache
// 26:29 - Neighbor cache
// 30:30 - Pipeline input arguments
// 31:31 - Pipeline return values

///////////////////////////////////////////////////////////////////////////////
// Pipeline Input / Output arguments

DECLARE_ALIGNED_ARRAY( advance_p_pipeline_args_t, 128, args, 1 );
DECLARE_ALIGNED_ARRAY( particle_mover_seg_t,      128, seg,  1 );
int pipeline_rank;
int n_pipeline;

///////////////////////////////////////////////////////////////////////////////
// Pipeline buffers

// NP_BLOCK_TARGET is such that a 16Kb of particle data is located in
// a particle block (the maximum amount of data that can be moved in a
// single DMA transfer).  Triple buffering requires 3 such particle
// blocks.  So local particle storage takes up 48Kb.  Each local
// particle has a mover associated with it such that a SPU pipeline
// cannot run out movers during particle processing.  This costs 24Kb.

#define NP_BLOCK_TARGET 512
DECLARE_ALIGNED_ARRAY( particle_t,       128, local_p,  3*NP_BLOCK_TARGET );
DECLARE_ALIGNED_ARRAY( particle_mover_t, 128, local_pm, 3*NP_BLOCK_TARGET );

///////////////////////////////////////////////////////////////////////////////
// External memory caches

// #define CACHE_STATS

// Since DMA transfers seem optimized for 128-bytes at aligned
// addresses, set up the cache lines to use 128-byte aligned lines
// 128-bytes in size.

#if 0
// Interpolator cache: 512 cached interpolators (64Kb). Roughly four
// second nearest neighborhoods (125 voxels) around a particle can exist
// within the cache.

#undef CACHE_NAME
#undef CACHED_TYPE
#undef CACHE_TYPE
#undef CACHELINE_LOG2SIZE
#undef CACHE_LOG2NWAY
#undef CACHE_LOG2NSETS
#undef CACHE_SET_TAGID
#undef CACHE_READ_X4
#undef CACHE_BASE_EA

#define CACHE_NAME           interpolator_cache
#define CACHED_TYPE          interpolator_t
#define CACHE_TYPE           0             /* r/o */
#define CACHELINE_LOG2SIZE   7             /* 1 per line - 128 byte lines */
#define CACHE_LOG2NWAY       2             /* 4 way */
#define CACHE_LOG2NSETS      7             /* 128 lines per way */
#define CACHE_SET_TAGID(set) (9+(set)&0x7) /* tags 9:16 */
#define CACHE_BASE_EA        args->f0
#include "cache-api.h"

#define PTR_INTERPOLATOR(v) cache_rw( interpolator_cache,               \
                                      (v)*sizeof(interpolator_t) )

#else
DECLARE_ALIGNED_ARRAY( interpolator_t, 128, interpolator_cache, 1 );
#define PTR_INTERPOLATOR(v) interpolator_cache
#endif

#if 0

// Accumulator cache: 512 cached accumulators (32Kb).  Roughly four
// second nearest neighborhood of accumulators can exist within the
// cache.

#undef CACHE_NAME
#undef CACHED_TYPE
#undef CACHE_TYPE
#undef CACHELINE_LOG2SIZE
#undef CACHE_LOG2NWAY
#undef CACHE_LOG2NSETS
#undef CACHE_SET_TAGID
#undef CACHE_READ_X4
#undef CACHE_BASE_EA

#define CACHE_NAME           accumulator_cache
#define CACHED_TYPE          accumulator_t
#define CACHE_TYPE           1              /* r/w */
#define CACHELINE_LOG2SIZE   7              /* 2 per line - 128 byte lines */
#define CACHE_LOG2NWAY       2              /* 4 way */
#define CACHE_LOG2NSETS      6              /* 64 lines per way */
#define CACHE_SET_TAGID(set) (17+(set)&0x7) /* tags 17:25 */
#define CACHE_BASE_EA        args->a0
#include "cache-api.h"

#define PTR_ACCUMULATOR(v) cache_rw( accumulator_cache,                 \
                                     (v)*sizeof(accumulator_t) )

#else
DECLARE_ALIGNED_ARRAY( accumulator_t, 128, accumulator_cache, 1 );
#define PTR_ACCUMULATOR(v) accumulator_cache
#endif

// Neighbor cache: 2048 cached cell adjacencies (16K).  Roughly three
// second nearest neighborhoods of voxel adjacencies can exist within the
// cache.

#undef CACHE_NAME
#undef CACHED_TYPE
#undef CACHE_TYPE
#undef CACHELINE_LOG2SIZE
#undef CACHE_LOG2NWAY
#undef CACHE_LOG2NSETS
#undef CACHE_SET_TAGID
#undef CACHE_READ_X4
#undef CACHE_BASE_EA

#define CACHE_NAME           neighbor_cache
#define CACHED_TYPE          int64_t
#define CACHE_TYPE           0              /* r/o */
#define CACHELINE_LOG2SIZE   7              /* 16 per line - 128 byte lines */
#define CACHE_LOG2NWAY       2              /* 4 way */
#define CACHE_LOG2NSETS      5              /* 32 lines per way */
#define CACHE_SET_TAGID(set) (25+(set)&0x3) /* tags 25:29 */
#define CACHE_BASE_EA        args->neighbor
#include "cache-api.h"

#define NEIGHBOR(v,face) cache_rd( neighbor_cache,                      \
                                   6*sizeof(int64_t)*(v) +              \
                                   (face)*sizeof(int64_t) )

///////////////////////////////////////////////////////////////////////////////
// Computational kernel

#if 0
// move_p moves the particle m->p by m->dispx, m->dispy, m->dispz
// depositing particle current as it goes. If the particle was moved
// sucessfully (particle mover is no longer in use), this returns 0.
// If the particle interacted with something this routine could not
// handle, this routine returns 1 (particle mover is still in use).
// On a successful move, the particle position is updated and dispx,
// dispy and dispz are zerod.  On a partial move, the particle
// position is updated to the point where the particle interacted and
// dispx, dispy, dispz contain the remaining particle displacement.
// The displacements are the physical displacments normalized current
// cell size.
//
// Because move_p is internal use only and frequently called, it does
// not check its input arguments.  Callers are responsible for
// insuring valid arguments.

int
move_p_spu( particle_t       * ALIGNED(32) p,
            particle_mover_t * ALIGNED(16) pm ) {
  float s_midx, s_midy, s_midz;
  float s_dispx, s_dispy, s_dispz;
  float s_dir[3];
  float v0, v1, v2, v3, v4, v5;
  int type;
  int64_t neighbor;
  float * ALIGNED(16) a;

  // FIXME: THIS CAN BE PARTIALLY HORIZONTAL SIMD ACCELERATED

  for(;;) {
    s_midx = p->dx;
    s_midy = p->dy;
    s_midz = p->dz;

    s_dispx = pm->dispx;
    s_dispy = pm->dispy;
    s_dispz = pm->dispz;

    s_dir[0] = (s_dispx>0) ? 1 : -1;
    s_dir[1] = (s_dispy>0) ? 1 : -1;
    s_dir[2] = (s_dispz>0) ? 1 : -1;
    
    // Compute the twice the fractional distance to each potential
    // streak/cell face intersection.

    v0 = (s_dispx==0) ? 3.4e38 : (s_dir[0]-s_midx)/s_dispx;
    v1 = (s_dispy==0) ? 3.4e38 : (s_dir[1]-s_midy)/s_dispy;
    v2 = (s_dispz==0) ? 3.4e38 : (s_dir[2]-s_midz)/s_dispz;

    // Determine the fractional length and type of current streak. The
    // streak ends on either the first face intersected by the
    // particle track or at the end of the particle track.
    // 
    //   type 0,1 or 2 ... streak ends on a x,y or z-face respectively
    //   type 3        ... streak ends at end of the particle track

    /**/      v3=2,  type=3;
    if(v0<v3) v3=v0, type=0;
    if(v1<v3) v3=v1, type=1;
    if(v2<v3) v3=v2, type=2;
    v3 *= 0.5;

    // Compute the midpoint and the normalized displacement of the streak

    s_dispx *= v3;
    s_dispy *= v3;
    s_dispz *= v3;

    s_midx += s_dispx;
    s_midy += s_dispy;
    s_midz += s_dispz;

    // Accumulate the streak.  Note: accumulator values are 4 times
    // the total physical charge that passed through the appropriate
    // current quadrant in a time-step

    v5 = p->q*s_dispx*s_dispy*s_dispz*(1./3.);
    a = (float * ALIGNED(64))PTR_ACCUMULATOR(p->i); // No need to lock

#   define ACCUMULATE_J(X,Y,Z,offset)                                   \
    v4  = p->q*s_disp##X; /* v2 = q ux                            */    \
    v1  = v4*s_mid##Y;    /* v1 = q ux dy                         */    \
    v0  = v4-v1;          /* v0 = q ux (1-dy)                     */    \
    v1 += v4;             /* v1 = q ux (1+dy)                     */    \
    v4  = 1+s_mid##Z;     /* v4 = 1+dz                            */    \
    v2  = v0*v4;          /* v2 = q ux (1-dy)(1+dz)               */    \
    v3  = v1*v4;          /* v3 = q ux (1+dy)(1+dz)               */    \
    v4  = 1-s_mid##Z;     /* v4 = 1-dz                            */    \
    v0 *= v4;             /* v0 = q ux (1-dy)(1-dz)               */    \
    v1 *= v4;             /* v1 = q ux (1+dy)(1-dz)               */    \
    v0 += v5;             /* v0 = q ux [ (1-dy)(1-dz) + uy*uz/3 ] */    \
    v1 -= v5;             /* v1 = q ux [ (1+dy)(1-dz) - uy*uz/3 ] */    \
    v2 -= v5;             /* v2 = q ux [ (1-dy)(1+dz) - uy*uz/3 ] */    \
    v3 += v5;             /* v3 = q ux [ (1+dy)(1+dz) + uy*uz/3 ] */    \
    a[offset+0] += v0;                                                  \
    a[offset+1] += v1;                                                  \
    a[offset+2] += v2;                                                  \
    a[offset+3] += v3

    ACCUMULATE_J( x,y,z, 0 );
    ACCUMULATE_J( y,z,x, 4 );
    ACCUMULATE_J( z,x,y, 8 );

#   undef ACCUMULATE_J

    // Compute the remaining particle displacment

    pm->dispx -= s_dispx;
    pm->dispy -= s_dispy;
    pm->dispz -= s_dispz;

    // Compute the new particle offset

    p->dx += s_dispx+s_dispx;
    p->dy += s_dispy+s_dispy;
    p->dz += s_dispz+s_dispz;

    // If an end streak, return success (should be ~50% of the time)

    if( type==3 ) return 0;

    // Determine if the cell crossed into a local cell or if it hit a
    // boundary.  Convert the coordinate system accordingly.  Note:
    // Crossing into a local cell should happen the other ~50% of
    // time; hitting a structure and parallel domain boundary should
    // usually be a rare event.  Note: the entry / exit coordinate for
    // the particle is guaranteed to be +/-1 _exactly_ for the
    // particle.

    // FIXME: WOULD SWITCHING ON TYPE BE FASTER ON THE SPUS?

    v0 = s_dir[type];
    neighbor = NEIGHBOR( p->i, ((v0>0)?3:0) + type );
    if( neighbor<args->rangel || neighbor>args->rangeh ) { // Hit a boundary
      (&(p->dx))[type] = v0;                               // Put on boundary
      if( neighbor!=reflect_particles ) return 1;          // Cannot handle it
      (&(p->ux))[type] = -(&(p->ux))[type];
      (&(pm->dispx))[type] = -(&(pm->dispx))[type];
    } else {
      p->i = neighbor - args->rangel; // Compute local index of neighbor
      /**/                            // Note: neighbor-args->rangel < 2^31 / 6
      (&(p->dx))[type] = -v0;         // Convert coordinate system
    }
  }
  return 0; // Never get here ... avoid compiler warning
}
#endif

// FIXME: Using restricted pointers makes this worse(!) on gcc??
// FIXME: Branch hints makes this worse(!) on gcc?? (No change on xlc.)

static int                             // Return number of movers used
advance_p_pipeline_spu( particle_t       * ALIGNED(128) p,  // Particle array
                        particle_mover_t * ALIGNED(16)  pm, // Mover array
                        int idx,       // Index of first particle
                        int ndq ) {    // Number of particle double quads
  USING_V4C;

  const vec_float4 qdt_2mc        = VEC_FLOAT4( args->qdt_2mc );
  const vec_float4 cdt_dx         = VEC_FLOAT4( args->cdt_dx  );
  const vec_float4 cdt_dy         = VEC_FLOAT4( args->cdt_dy  );
  const vec_float4 cdt_dz         = VEC_FLOAT4( args->cdt_dz  );
  const vec_float4 one            = VEC_FLOAT4(  1.           );
  const vec_float4 one_third      = VEC_FLOAT4(  1./3.        );
  const vec_float4 one_half       = VEC_FLOAT4(  1./2.        );
  const vec_float4 two_fifteenths = VEC_FLOAT4(  2./15.       );
  const vec_float4 neg_one        = VEC_FLOAT4( -1.           );

  vec_float4 p0r, p1r, p2r, p3r;                                     vec_float4 p4r, p5r, p6r, p7r;
  vec_float4 p0u, p1u, p2u, p3u;                                     vec_float4 p4u, p5u, p6u, p7u;

  vec_float4 dx, dy, dz;    vec_int4 i;                              vec_float4 dx_, dy_, dz_;    vec_int4 i_;
  vec_float4 ux, uy, uz, q;                                          vec_float4 ux_, uy_, uz_, q_;           

  vec_float4 ex0,  dexdy,  dexdz, d2exdydz;                          vec_float4 ex0_,  dexdy_,  dexdz_, d2exdydz_;
  vec_float4 ey0,  deydz,  deydx, d2eydzdx;                          vec_float4 ey0_,  deydz_,  deydx_, d2eydzdx_;
  vec_float4 ez0,  dezdx,  dezdy, d2ezdxdy;                          vec_float4 ez0_,  dezdx_,  dezdy_, d2ezdxdy_;
  vec_float4 cbx0, dcbxdx, cby0,  dcbydy;                            vec_float4 cbx0_, dcbxdx_, cby0_,  dcbydy_;
  vec_float4 cbz0, dcbzdz, v12,   v13;                               vec_float4 cbz0_, dcbzdz_, v12_,   v13_; 

  vec_float4 hax, hay, haz;                                          vec_float4 hax_, hay_, haz_;
  vec_float4 cbx, cby, cbz;                                          vec_float4 cbx_, cby_, cbz_;

  vec_float4 ux0, uy0, uz0;                                          vec_float4 ux0_, uy0_, uz0_;
  vec_float4 v14, cbs, ths, v15, v16, v17;                           vec_float4 v14_, cbs_, ths_, v15_, v16_, v17_;
  vec_float4 wx0, wy0, wz0;                                          vec_float4 wx0_, wy0_, wz0_;
  vec_float4 uxh, uyh, uzh;                                          vec_float4 uxh_, uyh_, uzh_;

  vec_float4 rgamma;                                                 vec_float4 rgamma_; 
  vec_float4 ddx, ddy, ddz;                                          vec_float4 ddx_, ddy_, ddz_;
  vec_float4 dxh, dyh, dzh;                                          vec_float4 dxh_, dyh_, dzh_;
  vec_float4 dx1, dy1, dz1;                                          vec_float4 dx1_, dy1_, dz1_;
  vec_uint4  outbnd;                                                 vec_uint4  outbnd_; 

  vec_float4 qa, ccc;                                                vec_float4 qa_, ccc_;
  vec_float4 a0x, a1x, a2x, a3x, a4x;                                vec_float4 a0x_, a1x_, a2x_, a3x_, a4x_;
  vec_float4 a0y, a1y, a2y, a3y, a4y;                                vec_float4 a0y_, a1y_, a2y_, a3y_, a4y_;
  vec_float4 a0z, a1z, a2z, a3z, a4z;                                vec_float4 a0z_, a1z_, a2z_, a3z_, a4z_;

  int i0, i1, i2, i3;                                                int i0_, i1_, i2_, i3_;

  const interpolator_t * ALIGNED(128) fi;
  accumulator_t * ALIGNED(64) ja;
  int nm = 0;

  // Process the particle quads for this pipeline

  for( ; ndq; ndq--, p+=8, idx+=8 ) {

    // Load the particle quad positions

    LOAD_4x1( &p[0].dx, p0r );                                       LOAD_4x1( &p[4].dx, p4r );
    LOAD_4x1( &p[1].dx, p1r );                                       LOAD_4x1( &p[5].dx, p5r );
    LOAD_4x1( &p[2].dx, p2r );                                       LOAD_4x1( &p[6].dx, p6r );
    LOAD_4x1( &p[3].dx, p3r );                                       LOAD_4x1( &p[7].dx, p7r );
    TRANSPOSE( p0r, p1r, p2r, p3r );                                 TRANSPOSE( p4r, p5r, p6r, p7r );
    dx = p0r;                                                        dx_ = p4r;
    dy = p1r;                                                        dy_ = p5r;
    dz = p2r;                                                        dz_ = p6r;
    i  = (vec_int4)p3r;                                              i_ = (vec_int4)p7r;

    // Interpolate fields

    i0 = EXTRACT( i, 0 );                                            i0_ = EXTRACT( i_, 0 );
    fi = PTR_INTERPOLATOR(i0);
    LOAD_4x1( &fi->ex,  ex0  );
    LOAD_4x1( &fi->ey,  ey0  );
    LOAD_4x1( &fi->ez,  ez0  );
    LOAD_4x1( &fi->cbx, cbx0 );
    LOAD_4x1( &fi->cbz, cbz0 );
    /**/                                                             fi = PTR_INTERPOLATOR(i0_);
    /**/                                                             LOAD_4x1( &fi->ex,  ex0_  );
    /**/                                                             LOAD_4x1( &fi->ey,  ey0_  );
    /**/                                                             LOAD_4x1( &fi->ez,  ez0_  );
    /**/                                                             LOAD_4x1( &fi->cbx, cbx0_ );
    /**/                                                             LOAD_4x1( &fi->cbz, cbz0_ );

    i1 = EXTRACT( i, 1 );                                            i1_ = EXTRACT( i_, 1 );
    fi = PTR_INTERPOLATOR(i1);
    LOAD_4x1( &fi->ex,  dexdy  );
    LOAD_4x1( &fi->ey,  deydz  );
    LOAD_4x1( &fi->ez,  dezdx  );
    LOAD_4x1( &fi->cbx, dcbxdx );
    LOAD_4x1( &fi->cbz, dcbzdz );
    /**/                                                             fi = PTR_INTERPOLATOR(i1_);
    /**/                                                             LOAD_4x1( &fi->ex,  dexdy_  );
    /**/                                                             LOAD_4x1( &fi->ey,  deydz_  );
    /**/                                                             LOAD_4x1( &fi->ez,  dezdx_  );
    /**/                                                             LOAD_4x1( &fi->cbx, dcbxdx_ );
    /**/                                                             LOAD_4x1( &fi->cbz, dcbzdz_ );

    i2 = EXTRACT( i, 2 );                                            i2_ = EXTRACT( i_, 2 );
    fi = PTR_INTERPOLATOR(i2);
    LOAD_4x1( &fi->ex,  dexdz );
    LOAD_4x1( &fi->ey,  deydx );
    LOAD_4x1( &fi->ez,  dezdy );
    LOAD_4x1( &fi->cbx, cby0  );
    LOAD_4x1( &fi->cbz, v12   );
    /**/                                                             fi = PTR_INTERPOLATOR(i2_);
    /**/                                                             LOAD_4x1( &fi->ex,  dexdz_ );
    /**/                                                             LOAD_4x1( &fi->ey,  deydx_ );
    /**/                                                             LOAD_4x1( &fi->ez,  dezdy_ );
    /**/                                                             LOAD_4x1( &fi->cbx, cby0_  );
    /**/                                                             LOAD_4x1( &fi->cbz, v12_   );

    i3 = EXTRACT( i, 3 );                                            i3_ = EXTRACT( i_, 3 );
    fi = PTR_INTERPOLATOR(i3);
    LOAD_4x1( &fi->ex,  d2exdydz );
    LOAD_4x1( &fi->ey,  d2eydzdx );
    LOAD_4x1( &fi->ez,  d2ezdxdy );
    LOAD_4x1( &fi->cbx, dcbydy   );
    LOAD_4x1( &fi->cbz, v13      );
    /**/                                                             fi = PTR_INTERPOLATOR(i3_);
    /**/                                                             LOAD_4x1( &fi->ex,  d2exdydz_ );
    /**/                                                             LOAD_4x1( &fi->ey,  d2eydzdx_ );
    /**/                                                             LOAD_4x1( &fi->ez,  d2ezdxdy_ );
    /**/                                                             LOAD_4x1( &fi->cbx, dcbydy_   );
    /**/                                                             LOAD_4x1( &fi->cbz, v13_      );

    TRANSPOSE( ex0, dexdy, dexdz, d2exdydz );                        TRANSPOSE( ex0_, dexdy_, dexdz_, d2exdydz_ );
    hax = MUL( qdt_2mc, FMA( FMA( d2exdydz, dy, dexdz ), dz,
                             FMA( dexdy,    dy, ex0   ) ) );         hax_ = MUL( qdt_2mc, FMA( FMA( d2exdydz_, dy_, dexdz_ ), dz_,
                                                                                               FMA( dexdy_,    dy_, ex0_   ) ) );

    TRANSPOSE( ey0, deydz, deydx, d2eydzdx );                        TRANSPOSE( ey0_, deydz_, deydx_, d2eydzdx_ );
    hay = MUL( qdt_2mc, FMA( FMA( d2eydzdx, dz, deydx ), dx,
                             FMA( deydz,    dz, ey0 ) ) );           hay_ = MUL( qdt_2mc, FMA( FMA( d2eydzdx_, dz_, deydx_ ), dx_,
                                                                                               FMA( deydz_,    dz_, ey0_ ) ) );

    TRANSPOSE( ez0, dezdx, dezdy, d2ezdxdy );                        TRANSPOSE( ez0_, dezdx_, dezdy_, d2ezdxdy_ );
    haz = MUL( qdt_2mc, FMA( FMA( d2ezdxdy, dx, dezdy ), dy,
                             FMA( dezdx,    dx, ez0 ) ) );           haz_ = MUL( qdt_2mc, FMA( FMA( d2ezdxdy_, dx_, dezdy_ ), dy_,
                                                                                               FMA( dezdx_,    dx_, ez0_ ) ) );

    TRANSPOSE( cbx0, dcbxdx, cby0, dcbydy );                         TRANSPOSE( cbx0_, dcbxdx_, cby0_, dcbydy_ );
    cbx = FMA( dcbxdx, dx, cbx0 );                                   cbx_ = FMA( dcbxdx_, dx_, cbx0_ );
    cby = FMA( dcbydy, dy, cby0 );                                   cby_ = FMA( dcbydy_, dy_, cby0_ );

    HALF_TRANSPOSE( cbz0, dcbzdz, v12, v13 );                        HALF_TRANSPOSE( cbz0_, dcbzdz_, v12_, v13_ );
    cbz = FMA( dcbzdz, dz, cbz0 );                                   cbz_ = FMA( dcbzdz_, dz_, cbz0_ );

    // Update momentum.  Note: Could eliminate a dependency in v14 calc
    // if willing to play fast and loose with numerics (saves about a spu
    // clock per particle).

    LOAD_4x1( &p[0].ux, p0u );                                       LOAD_4x1( &p[4].ux, p4u );
    LOAD_4x1( &p[1].ux, p1u );                                       LOAD_4x1( &p[5].ux, p5u );
    LOAD_4x1( &p[2].ux, p2u );                                       LOAD_4x1( &p[6].ux, p6u );
    LOAD_4x1( &p[3].ux, p3u );                                       LOAD_4x1( &p[7].ux, p7u );
    TRANSPOSE( p0u, p1u, p2u, p3u );                                 TRANSPOSE( p4u, p5u, p6u, p7u );
    ux  = p0u;                                                       ux_  = p4u;
    uy  = p1u;                                                       uy_  = p5u;
    uz  = p2u;                                                       uz_  = p6u;
    q   = p3u;                                                       q_   = p7u;
    ux0 = ADD( ux, hax );                                            ux0_ = ADD( ux_, hax_ );
    uy0 = ADD( uy, hay );                                            uy0_ = ADD( uy_, hay_ );
    uz0 = ADD( uz, haz );                                            uz0_ = ADD( uz_, haz_ );
    v14 = MUL( qdt_2mc, RSQRT( ADD( one, FMA( ux0,ux0, FMA( uy0,uy0, MUL( uz0,uz0 ) ) ) ) ) );
    /**/                                                             v14_ = MUL( qdt_2mc, RSQRT( ADD( one, FMA( ux0_,ux0_, FMA( uy0_,uy0_, MUL( uz0_,uz0_ ) ) ) ) ) );
    cbs = FMA( cbx,cbx, FMA( cby,cby, MUL(cbz,cbz) ) );              cbs_ = FMA( cbx_,cbx_, FMA( cby_,cby_, MUL(cbz_,cbz_) ) ); 
    ths = MUL( MUL( v14,v14 ), cbs );                                ths_ = MUL( MUL( v14_,v14_ ), cbs_ );
    v15 = MUL( v14, FMA( FMA( two_fifteenths, ths, one_third ), ths, one ) );
    /**/                                                             v15_ = MUL( v14_, FMA( FMA( two_fifteenths, ths_, one_third ), ths_, one ) );
    v16 = MUL( v15, RCP( FMA( MUL(v15,v15), cbs, one ) ) );          v16_ = MUL( v15_, RCP( FMA( MUL(v15_,v15_), cbs_, one ) ) );
    v17 = ADD( v16, v16 );                                           v17_ = ADD( v16_, v16_ );
    wx0 =      FMA( FMS( uy0,cbz, MUL(uz0,cby) ), v15, ux0 );        wx0_ =      FMA( FMS( uy0_,cbz_, MUL(uz0_,cby_) ), v15_, ux0_ );
    wy0 =      FMA( FMS( uz0,cbx, MUL(ux0,cbz) ), v15, uy0 );        wy0_ =      FMA( FMS( uz0_,cbx_, MUL(ux0_,cbz_) ), v15_, uy0_ );
    wz0 =      FMA( FMS( ux0,cby, MUL(uy0,cbx) ), v15, uz0 );        wz0_ =      FMA( FMS( ux0_,cby_, MUL(uy0_,cbx_) ), v15_, uz0_ );
    uxh = ADD( FMA( FMS( wy0,cbz, MUL(wz0,cby) ), v17, ux0 ), hax ); uxh_ = ADD( FMA( FMS( wy0_,cbz_, MUL(wz0_,cby_) ), v17_, ux0_ ), hax_ );
    uyh = ADD( FMA( FMS( wz0,cbx, MUL(wx0,cbz) ), v17, uy0 ), hay ); uyh_ = ADD( FMA( FMS( wz0_,cbx_, MUL(wx0_,cbz_) ), v17_, uy0_ ), hay_ );
    uzh = ADD( FMA( FMS( wx0,cby, MUL(wy0,cbx) ), v17, uz0 ), haz ); uzh_ = ADD( FMA( FMS( wx0_,cby_, MUL(wy0_,cbx_) ), v17_, uz0_ ), haz_ );
    p0u = uxh;                                                       p4u  = uxh_;
    p1u = uyh;                                                       p5u  = uyh_;
    p2u = uzh;                                                       p6u  = uzh_;
    /* p3u is unchanged */                                           /* p7u is unchanged */
    TRANSPOSE( p0u, p1u, p2u, p3u );                                 TRANSPOSE( p4u, p5u, p6u, p7u );
    STORE_4x1( p0u, &p[0].ux );                                      STORE_4x1( p4u, &p[4].ux );
    STORE_4x1( p1u, &p[1].ux );                                      STORE_4x1( p5u, &p[5].ux );
    STORE_4x1( p2u, &p[2].ux );                                      STORE_4x1( p6u, &p[6].ux );
    STORE_4x1( p3u, &p[3].ux );                                      STORE_4x1( p7u, &p[7].ux );
    
    // Update the position of inbnd particles

    rgamma = RSQRT( ADD( one, FMA( uxh,uxh, FMA( uyh,uyh, MUL(uzh,uzh) ) ) ) );
    /**/                                                             rgamma_ = RSQRT( ADD( one, FMA( uxh_,uxh_, FMA( uyh_,uyh_, MUL(uzh_,uzh_) ) ) ) );
    ddx    = MUL( MUL( uxh, cdt_dx ), rgamma );                      ddx_    = MUL( MUL( uxh_, cdt_dx ), rgamma_ );
    ddy    = MUL( MUL( uyh, cdt_dy ), rgamma );                      ddy_    = MUL( MUL( uyh_, cdt_dy ), rgamma_ );
    ddz    = MUL( MUL( uzh, cdt_dz ), rgamma );                      ddz_    = MUL( MUL( uzh_, cdt_dz ), rgamma_ );
    dxh    = ADD( dx,  ddx );                                        dxh_    = ADD( dx_,  ddx_ );
    dyh    = ADD( dy,  ddy );                                        dyh_    = ADD( dy_,  ddy_ );
    dzh    = ADD( dz,  ddz );                                        dzh_    = ADD( dz_,  ddz_ );        
    dx1    = ADD( dxh, ddx );                                        dx1_    = ADD( dxh_, ddx_ );
    dy1    = ADD( dyh, ddy );                                        dy1_    = ADD( dyh_, ddy_ );
    dz1    = ADD( dzh, ddz );                                        dz1_    = ADD( dzh_, ddz_ );
    outbnd = OR( OR( OR( CMPLT(dx1,neg_one), CMPGT(dx1,one) ),
                     OR( CMPLT(dy1,neg_one), CMPGT(dy1,one) ) ),
                     OR( CMPLT(dz1,neg_one), CMPGT(dz1,one) ) );     outbnd_ = OR( OR( OR( CMPLT(dx1_,neg_one), CMPGT(dx1_,one) ),
                                                                                       OR( CMPLT(dy1_,neg_one), CMPGT(dy1_,one) ) ),
                                                                                       OR( CMPLT(dz1_,neg_one), CMPGT(dz1_,one) ) );
    p0r    = MERGE( outbnd, dx, dx1 );                               p4r     = MERGE( outbnd_, dx_, dx1_ );
    p1r    = MERGE( outbnd, dy, dy1 );                               p5r     = MERGE( outbnd_, dy_, dy1_ );
    p2r    = MERGE( outbnd, dz, dz1 );                               p6r     = MERGE( outbnd_, dz_, dz1_ );
    /* p3r is unchanged */                                           /* p7r is unchanged */
    TRANSPOSE( p0r, p1r, p2r, p3r );                                 TRANSPOSE( p4r, p5r, p6r, p7r );
    STORE_4x1( p0r, &p[0].dx );                                      STORE_4x1( p4r, &p[4].dx );
    STORE_4x1( p1r, &p[1].dx );                                      STORE_4x1( p5r, &p[5].dx );
    STORE_4x1( p2r, &p[2].dx );                                      STORE_4x1( p6r, &p[6].dx );
    STORE_4x1( p3r, &p[3].dx );                                      STORE_4x1( p7r, &p[7].dx );
   
    // Accumulate current of inbnd particles
    // Note: accumulator values are 4 times the total physical charge that
    // passed through the appropriate current quadrant in a time-step

    qa  = CZERO( outbnd, q );                                        qa_  = CZERO( outbnd_, q_ );
    ccc = MUL( MUL( one_third, qa ), MUL( ddx, MUL( ddy, ddz ) ) );  ccc_ = MUL( MUL( one_third, qa_ ), MUL( ddx_, MUL( ddy_, ddz_ ) ) );

#   define ACCUMULATE_J(X,Y,Z)                                                                          \
    a4##X = MUL(qa,   dd##X);                                        a4##X##_ = MUL(qa_,     dd##X##_); \
    a1##X = MUL(a4##X,d##Y##h);                                      a1##X##_ = MUL(a4##X##_,d##Y##h_); \
    a0##X = SUB(a4##X,a1##X);                                        a0##X##_ = SUB(a4##X##_,a1##X##_); \
    a1##X = ADD(a1##X,a4##X);                                        a1##X##_ = ADD(a1##X##_,a4##X##_); \
    a4##X = ADD(one,  d##Z##h);                                      a4##X##_ = ADD(one,     d##Z##h_); \
    a2##X = MUL(a0##X,a4##X);                                        a2##X##_ = MUL(a0##X##_,a4##X##_); \
    a3##X = MUL(a1##X,a4##X);                                        a3##X##_ = MUL(a1##X##_,a4##X##_); \
    a4##X = SUB(one,  d##Z##h);                                      a4##X##_ = SUB(one,     d##Z##h_); \
    a0##X = MUL(a0##X,a4##X);                                        a0##X##_ = MUL(a0##X##_,a4##X##_); \
    a1##X = MUL(a1##X,a4##X);                                        a1##X##_ = MUL(a1##X##_,a4##X##_); \
    a0##X = ADD(a0##X,ccc);                                          a0##X##_ = ADD(a0##X##_,ccc_);     \
    a1##X = SUB(a1##X,ccc);                                          a1##X##_ = SUB(a1##X##_,ccc_);     \
    a2##X = SUB(a2##X,ccc);                                          a2##X##_ = SUB(a2##X##_,ccc_);     \
    a3##X = ADD(a3##X,ccc);                                          a3##X##_ = ADD(a3##X##_,ccc_);     \
    TRANSPOSE(a0##X,a1##X,a2##X,a3##X);                              TRANSPOSE(a0##X##_,a1##X##_,a2##X##_,a3##X##_)

    ACCUMULATE_J( x,y,z );
    ACCUMULATE_J( y,z,x );
    ACCUMULATE_J( z,x,y );

#   undef ACCUMULATE_J

    ja = PTR_ACCUMULATOR(i0);
    INCREMENT_4x1( ja->jx, a0x );
    INCREMENT_4x1( ja->jy, a0y );
    INCREMENT_4x1( ja->jz, a0z );

    ja = PTR_ACCUMULATOR(i1);
    INCREMENT_4x1( ja->jx, a1x );
    INCREMENT_4x1( ja->jy, a1y );
    INCREMENT_4x1( ja->jz, a1z );

    ja = PTR_ACCUMULATOR(i2);
    INCREMENT_4x1( ja->jx, a2x );
    INCREMENT_4x1( ja->jy, a2y );
    INCREMENT_4x1( ja->jz, a2z );

    ja = PTR_ACCUMULATOR(i3);
    INCREMENT_4x1( ja->jx, a3x );
    INCREMENT_4x1( ja->jy, a3y );
    INCREMENT_4x1( ja->jz, a3z );

    /**/                                                             ja = PTR_ACCUMULATOR(i0_);
    /**/                                                             INCREMENT_4x1( ja->jx, a0x_ );
    /**/                                                             INCREMENT_4x1( ja->jy, a0y_ );
    /**/                                                             INCREMENT_4x1( ja->jz, a0z_ );

    /**/                                                             ja = PTR_ACCUMULATOR(i1_);
    /**/                                                             INCREMENT_4x1( ja->jx, a1x_ );
    /**/                                                             INCREMENT_4x1( ja->jy, a1y_ );
    /**/                                                             INCREMENT_4x1( ja->jz, a1z_ );

    /**/                                                             ja = PTR_ACCUMULATOR(i2_);
    /**/                                                             INCREMENT_4x1( ja->jx, a2x_ );
    /**/                                                             INCREMENT_4x1( ja->jy, a2y_ );
    /**/                                                             INCREMENT_4x1( ja->jz, a2z_ );

    /**/                                                             ja = PTR_ACCUMULATOR(i3_);
    /**/                                                             INCREMENT_4x1( ja->jx, a3x_ );
    /**/                                                             INCREMENT_4x1( ja->jy, a3y_ );
    /**/                                                             INCREMENT_4x1( ja->jz, a3z_ );

    // Update position and accumulate outbnd

#   if 0
#   define MOVE_OUTBND( N )                   \
    if( unlikely( EXTRACT( outbnd, N ) ) ) {  \
      pm->dispx = EXTRACT( ddx, N );          \
      pm->dispy = EXTRACT( ddy, N );          \
      pm->dispz = EXTRACT( ddz, N );          \
      pm->i     = idx+N;                      \
      if( unlikely( move_p_spu( p+N, pm ) ) ) pm++, nm++; \
    }
#   define MOVE_OUTBND_( N )                                                                                     \
                                                                     if( unlikely( EXTRACT( outbnd_, N ) ) ) {   \
                                                                       pm->dispx = EXTRACT( ddx_, N );           \
                                                                       pm->dispy = EXTRACT( ddy_, N );           \
                                                                       pm->dispz = EXTRACT( ddz_, N );           \
                                                                       pm->i     = idx+4+N;                      \
                                                                       if( unlikely( move_p_spu( p+4+N, pm ) ) ) pm++, nm++; \
                                                                     }

    MOVE_OUTBND(0);
    MOVE_OUTBND(1);
    MOVE_OUTBND(2);
    MOVE_OUTBND(3);
    /**/                                                             MOVE_OUTBND_(0);
    /**/                                                             MOVE_OUTBND_(1);
    /**/                                                             MOVE_OUTBND_(2);
    /**/                                                             MOVE_OUTBND_(3);

#   undef MOVE_OUTBND
#   undef MOVE_OUTBND_
#   endif

  }

  return nm;
}

///////////////////////////////////////////////////////////////////////////////
// main (workload distribution and data buffering)

// FIXME: util functionality is not compiled for the spu
// FIXME: UNIQUE EXTENSIONS FOR SPU OBJS AND EXECS!

int
main( uint64_t spu_id,
      uint64_t argp,
      uint64_t envp ) {
  particle_t       * ALIGNED(128) p_block[3];
  int idx[3];
  int np_block[3];

  particle_mover_t * ALIGNED(128) pm_block[3];
  int nm_block[3];

  int buffer;
  const int next[3] = { 1, 2, 0 };
  const int prev[3] = { 2, 0, 1 };

  int np, next_idx, itmp;

# ifdef IN_HARNESS
  prof_clear();
  prof_start();
# endif

  // Get the pipeline arguments from the dispatcher

  mfc_get( args,
           argp,
           sizeof(*args),
           30, 0, 0 );
  mfc_write_tag_mask( (1<<30) );
  mfc_read_tag_status_all();

  pipeline_rank = envp & 0xffffffff;
  n_pipeline    = envp >> 32; // Note: pipeline_rank<n_pipeline

  // Determine which particle quads this pipeline processes

  DISTRIBUTE( args->np, 16, pipeline_rank, n_pipeline, next_idx, np );

  // Determine which movers are reserved for this pipeline
  // Movers (16 bytes) are reserved for pipelines in multiples of 8
  // such that the set of particle movers reserved for a pipeline is
  // 128-bit aligned and a multiple of 128-bits in size. 

  args->max_nm -= args->np&15; // Insure host gets enough
  if( args->max_nm<0 ) args->max_nm = 0;
  DISTRIBUTE( args->max_nm, 8, pipeline_rank, n_pipeline, itmp, seg->max_nm );
  seg->pm        = args->pm + itmp*sizeof(particle_mover_t);
  seg->nm        = 0;
  seg->n_ignored = 0;

  // Determine which accumulator array is reserved for this pipeline

  args->a0 += sizeof(accumulator_t)*(1+pipeline_rank)*
    POW2_CEIL((args->nx+2)*(args->ny+2)*(args->nz+2),2);

  // Process the particles assigned to this pipeline with triple buffering

# define BEGIN_GET_PBLOCK(buffer) do {                          \
                                                                \
    /* Determine the start and size of the block */             \
    idx[buffer]      = next_idx;                                \
    np_block[buffer] = NP_BLOCK_TARGET;                         \
    if( np_block[buffer]>np ) np_block[buffer] = np;            \
    next_idx += np_block[buffer];                               \
    np       -= np_block[buffer];                               \
                                                                \
    /* If we have a block, start reading it into the buffer */  \
    if( np_block[buffer] )                                      \
      mfc_get( p_block[buffer],                                 \
               args->p0 + idx[buffer]*sizeof(particle_t),       \
               np_block[buffer]*sizeof(particle_t),             \
               3*(buffer)+0, 0, 0 );                            \
  } while(0)
  
# define END_GET_PBLOCK(buffer) do {                                    \
    /* If we have a block, stop reading it into the buffer */           \
    if( np_block[buffer] ) {                                            \
      mfc_write_tag_mask( (np_block[buffer]?1:0)<<(3*(buffer)+0) );     \
      mfc_read_tag_status_all();                                        \
    }                                                                   \
  } while(0)
  
# define PROCESS_PBLOCK(buffer)                                         \
  nm_block[buffer] = advance_p_pipeline_spu( p_block[buffer],           \
                                             pm_block[buffer],          \
                                             idx[buffer],               \
                                             np_block[buffer]>>3 )

  // FIXME: mfc list for this??
# define BEGIN_PUT_PBLOCK(buffer) do {                                  \
                                                                        \
    /* If we have a block, begin writing the buffer into the block */   \
    if( np_block[buffer] )                                              \
      mfc_put( p_block[buffer],                                         \
               args->p0 + idx[buffer]*sizeof(particle_t),               \
               np_block[buffer]*sizeof(particle_t),                     \
               3*(buffer)+1, 0, 0 );                                    \
                                                                        \
    /* Begin writing the movers corresponding to this block */          \
    /* Ignore movers that would overflow the mover segment */           \
    itmp = seg->nm + nm_block[buffer] - seg->max_nm;                    \
    if( itmp > 0 ) {                                                    \
      seg->n_ignored += itmp;                                           \
      nm_block[buffer] = seg->max_nm - seg->nm;                         \
    }                                                                   \
    if( nm_block[buffer] ) {                                            \
      mfc_put( pm_block[buffer],                                        \
               seg->pm + seg->nm*sizeof(particle_mover_t),              \
               nm_block[buffer]*sizeof(particle_mover_t),               \
               3*(buffer)+2, 0, 0 );                                    \
      seg->nm += nm_block[buffer];                                      \
    }                                                                   \
                                                                        \
  } while(0)

# define END_PUT_PBLOCK(buffer) do {                                    \
    if( np_block[buffer] || nm_block[buffer] ) {                        \
      mfc_write_tag_mask( (np_block[buffer]?1:0)<<(3*(buffer)+1) |      \
                          (nm_block[buffer]?1:0)<<(3*(buffer)+2) );     \
      mfc_read_tag_status_all();                                        \
    }                                                                   \
    np_block[buffer] = 0;                                               \
    nm_block[buffer] = 0;                                               \
  } while(0)

  p_block[0]  = local_p;                      np_block[0] = 0;
  p_block[1]  = local_p  + NP_BLOCK_TARGET;   np_block[1] = 0;
  p_block[2]  = local_p  + NP_BLOCK_TARGET*2; np_block[2] = 0;

  pm_block[0] = local_pm;                     nm_block[0] = 0;
  pm_block[1] = local_pm + NP_BLOCK_TARGET;   nm_block[1] = 0;
  pm_block[2] = local_pm + NP_BLOCK_TARGET*2; nm_block[2] = 0;

  BEGIN_GET_PBLOCK(0);
  for( buffer=0; np_block[buffer]>0; buffer=next[buffer] ) {
    BEGIN_GET_PBLOCK( next[buffer] );
    END_GET_PBLOCK( buffer );
    PROCESS_PBLOCK( buffer );
    BEGIN_PUT_PBLOCK( buffer );
    END_PUT_PBLOCK( prev[buffer] );
  }
  END_PUT_PBLOCK( prev[buffer] );

  // Flush the caches. Only accumulator_cache needs to be flushed as
  // the other caches are read-only.

# if 0
  cache_flush( accumulator_cache );
# endif

  // Write the pipeline return values back

  mfc_put( seg,
           args->seg + pipeline_rank*sizeof(particle_mover_seg_t),
           sizeof(particle_mover_seg_t),
           31, 0, 0 );
  mfc_write_tag_mask( (1<<31) );
  mfc_read_tag_status_all();
  
# ifdef IN_HARNESS
  prof_stop();
# endif

  return 0;
}

#endif // CELL_SPU_BUILD
