#define IN_particle_pipeline
#include <particle_pipelines.h>

#ifdef CELL_SPU_BUILD

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

using namespace v4;

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

// Since DMA transfers seem optimized for 128-bytes at aligned
// addresses, set up the cache lines to use 128-byte aligned lines
// 128-bytes in size.

// Note: XLC cannot do caches this large!  However, XLC either hates
// this loop or completely miscompiles it (it was over 4X slower than
// on GCC).

// Interpolator cache: 512 cached interpolators (64Kb).  Roughly four
// second nearest neighborhoods (125 voxels) around the current
// particle can exist within the cache.

#undef CACHE_NAME
#undef CACHED_TYPE
#undef CACHE_TYPE
#undef CACHELINE_LOG2SIZE
#undef CACHE_LOG2NWAY
#undef CACHE_LOG2NSETS
#undef CACHE_SET_TAGID
#undef CACHE_READ_X4
#undef CACHE_STATS

#define CACHE_NAME           interpolator_cache
#define CACHED_TYPE          interpolator_t
#define CACHE_TYPE           0             /* r/o */
#define CACHELINE_LOG2SIZE   7             /* 1 per line - 128 byte lines */
#define CACHE_LOG2NWAY       2             /* 4 way */
#define CACHE_LOG2NSETS      7             /* 128 lines per way */
#define CACHE_SET_TAGID(set) (9+(set)&0x7) /* tags 9:16 */
#include <cache-api.h>

#define PTR_INTERPOLATOR(v) cache_rw( interpolator_cache,               \
                                      args->f0 +                        \
                                      (v)*sizeof(interpolator_t) )

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
#undef CACHE_STATS

#define CACHE_NAME           accumulator_cache
#define CACHED_TYPE          accumulator_t
#define CACHE_TYPE           1              /* r/w */
#define CACHELINE_LOG2SIZE   7              /* 2 per line - 128 byte lines */
#define CACHE_LOG2NWAY       2              /* 4 way */
#define CACHE_LOG2NSETS      6              /* 64 lines per way */
#define CACHE_SET_TAGID(set) (17+(set)&0x7) /* tags 17:25 */
#include <cache-api.h>

#define PTR_ACCUMULATOR(v) cache_rw( accumulator_cache,                 \
                                     args->a0 +                         \
                                     (v)*sizeof(accumulator_t) )

// Neighbor cache: 2048 cached cell adjacencies (16K).  Roughly three
// complete second nearest neighborhoods of voxel adjacencies can
// exist within the cache.

#undef CACHE_NAME
#undef CACHED_TYPE
#undef CACHE_TYPE
#undef CACHELINE_LOG2SIZE
#undef CACHE_LOG2NWAY
#undef CACHE_LOG2NSETS
#undef CACHE_SET_TAGID
#undef CACHE_READ_X4
#undef CACHE_STATS

#define CACHE_NAME           neighbor_cache
#define CACHED_TYPE          int64_t
#define CACHE_TYPE           0              /* r/o */
#define CACHELINE_LOG2SIZE   7              /* 16 per line - 128 byte lines */
#define CACHE_LOG2NWAY       2              /* 4 way */
#define CACHE_LOG2NSETS      5              /* 32 lines per way */
#define CACHE_SET_TAGID(set) (25+(set)&0x3) /* tags 25:29 */
#include <cache-api.h>

#define NEIGHBOR(v,face) cache_rd( neighbor_cache,                      \
                                   args->neighbor +                     \
                                   6*sizeof(int64_t)*(v) +              \
                                   (face)*sizeof(int64_t) )

///////////////////////////////////////////////////////////////////////////////
// Computational kernel

// move_p moves the particle m->p by m->dispx, m->dispy, m->dispz
// depositing particle current as it goes. If the particle was moved
// sucessfully (particle mover is no longer in use) returns 0. If the
// particle interacted with something this routine could not handle,
// this routine returns 1 (particle mover is still in use). On a
// successful move, the particle position is updated and m->dispx,
// m->dispy and m->dispz are zerod. On a partial move, the particle
// position is updated to the point where the particle interacted and
// m->dispx, m->dispy, m->dispz contains the remaining particle
// displacement. The displacements are the physical displacments
// normalized current cell cell size.
//
// Because move_p is internal use only and frequently called, it does
// not check its input arguments. Higher level routines are
// responsible for insuring valid arguments.

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

  // FIXME: THIS CAN BE HORIZONAL V4 ACCELERATED (AT LEAST PARTIALLY)!

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
    a = (float * ALIGNED(16))PTR_ACCUMULATOR(p->i); // No need to lock

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

static int // Return number of movers used
advance_p_pipeline_spu( particle_t       * ALIGNED(128) p,  // Particle array
                        particle_mover_t * ALIGNED(16)  pm, // Mover array
                        int idx,   // Index of first particle
                        int nq ) { // Number of particle quads

  const v4float qdt_2mc( args->qdt_2mc );
  const v4float cdt_dx(  args->cdt_dx  );
  const v4float cdt_dy(  args->cdt_dy  );
  const v4float cdt_dz(  args->cdt_dz  );
  const v4float one(1.);
  const v4float one_third(1./3.);
  const v4float two_fifteenths(2./15.);
  const v4float neg_one(-1.);

  v4float dx, dy, dz; v4int   ii; // particle position
  v4float ux, uy, uz; v4float q;  // particle momentum and charge

  v4float ex0,  dexdy,  dexdz, d2exdydz; // ex interp coeff
  v4float ey0,  deydz,  deydx, d2eydzdx; // ey interp coeff
  v4float ez0,  dezdx,  dezdy, d2ezdxdy; // ez interp coeff
  v4float cbx0, dcbxdx, cby0,  dcbydy;   // bx and by interp coeff
  v4float cbz0, dcbzdz, v12,   v13;      // bz interp coeff and scratch

  v4float hax, hay, haz; // half particle electric field acceleration
  v4float cbx, cby, cbz; // particle magnetic field

  v4float ux0, uy0, uz0; // pre Boris rotation momentum
  v4float v14, cbs, ths, v15, v16, v17; // Boris rotation scalars
  v4float wx0, wy0, wz0; // intermediate Boris rotation momentum
  v4float uxh, uyh, uzh; // new particle momentum

  v4float rgamma;        // new inverse gamma
  v4float ddx, ddy, ddz; // cell-normalized particle displacment
  v4float dxh, dyh, dzh; // intermediate position relative to starting cell
  v4float dx1, dy1, dz1; // new position relative to starting cell
  v4int   outbnd;        // Boolean (true if particle exited cell)

  v4float qa, ccc;
  v4float a0x, a1x, a2x, a3x, a4x;
  v4float a0y, a1y, a2y, a3y, a4y;
  v4float a0z, a1z, a2z, a3z, a4z;

  const interpolator_t * ALIGNED(128) fi;
  accumulator_t * ALIGNED(64) ja;
  int nm = 0, ii0, ii1, ii2, ii3;

  // Process the particle quads for this pipeline

  for( ; nq; nq--, p+=4, idx+=4 ) {

    // Load the particle quad positions

    load_4x4_tr( &p[0].dx, &p[1].dx, &p[2].dx, &p[3].dx, dx, dy, dz, ii );

    // Interpolate fields

    ii0 = ii(0);
    fi  = PTR_INTERPOLATOR(ii0);
    load_4x1( &fi->ex,  ex0  );
    load_4x1( &fi->ey,  ey0  );
    load_4x1( &fi->ez,  ez0  );
    load_4x1( &fi->cbx, cbx0 );
    load_4x1( &fi->cbz, cbz0 );

    ii1 = ii(1);
    fi  = PTR_INTERPOLATOR(ii1);
    load_4x1( &fi->ex,  dexdy  );
    load_4x1( &fi->ey,  deydz  );
    load_4x1( &fi->ez,  dezdx  );
    load_4x1( &fi->cbx, dcbxdx );
    load_4x1( &fi->cbz, dcbzdz );

    ii2 = ii(2);
    fi  = PTR_INTERPOLATOR(ii2);
    load_4x1( &fi->ex,  dexdz );
    load_4x1( &fi->ey,  deydx );
    load_4x1( &fi->ez,  dezdy );
    load_4x1( &fi->cbx, cby0  );
    load_4x1( &fi->cbz, v12   );

    ii3 = ii(3);
    fi  = PTR_INTERPOLATOR(ii3);
    load_4x1( &fi->ex,  d2exdydz );
    load_4x1( &fi->ey,  d2eydzdx );
    load_4x1( &fi->ez,  d2ezdxdy );
    load_4x1( &fi->cbx, dcbydy   );
    load_4x1( &fi->cbz, v13      );

    transpose( ex0, dexdy, dexdz, d2exdydz );
    hax = qdt_2mc*fma( fma( d2exdydz, dy, dexdz ), dz, fma( dexdy, dy, ex0 ) );

    transpose( ey0, deydz, deydx, d2eydzdx );
    hay = qdt_2mc*fma( fma( d2eydzdx, dz, deydx ), dx, fma( deydz, dz, ey0 ) );

    transpose( ez0, dezdx, dezdy, d2ezdxdy );
    haz = qdt_2mc*fma( fma( d2ezdxdy, dx, dezdy ), dy, fma( dezdx, dx, ez0 ) );

    transpose( cbx0, dcbxdx, cby0, dcbydy );
    cbx = fma( dcbxdx, dx, cbx0 );
    cby = fma( dcbydy, dy, cby0 );

    transpose( cbz0, dcbzdz, v12, v13 ); // FIXME: Could be optimized a bit
    cbz = fma( dcbzdz, dz, cbz0 );

    // Update momentum.  Note: Could eliminate a dependency in v14 calc
    // if willing to play extra fast and loose with numerics (saves
    // about a spu clock per particle).

    load_4x4_tr( &p[0].ux, &p[1].ux, &p[2].ux, &p[3].ux, ux, uy, uz, q  );
    ux0 = ux + hax;
    uy0 = uy + hay;
    uz0 = uz + haz;
    v14 = qdt_2mc*rsqrt( one + fma( ux0,ux0, fma( uy0,uy0, uz0*uz0 ) ) );
    cbs = fma( cbx,cbx, fma( cby,cby, cbz*cbz ) ); // |cB|^2
    ths = (v14*v14)*cbs;                           // |theta|^2 = |wc dt/2|^2
    v15 = v14*fma( fma( two_fifteenths, ths, one_third ), ths, one );
    v16 = v15*rcp( fma( v15*v15, cbs, one ) );
    v17 = v16 + v16;
    wx0 = fma( fms( uy0,cbz, uz0*cby ), v15, ux0 );
    wy0 = fma( fms( uz0,cbx, ux0*cbz ), v15, uy0 );
    wz0 = fma( fms( ux0,cby, uy0*cbx ), v15, uz0 );
    uxh = fma( fms( wy0,cbz, wz0*cby ), v17, ux0 ) + hax;
    uyh = fma( fms( wz0,cbx, wx0*cbz ), v17, uy0 ) + hay;
    uzh = fma( fms( wx0,cby, wy0*cbx ), v17, uz0 ) + haz;
    store_4x4_tr( uxh, uyh, uzh,  q, &p[0].ux, &p[1].ux, &p[2].ux, &p[3].ux );
    
    // Update the position of inbnd particles

    rgamma = rsqrt( one + fma( uxh,uxh, fma( uyh,uyh, uzh*uzh ) ) );
    ddx    = (uxh*cdt_dx)*rgamma;
    ddy    = (uyh*cdt_dy)*rgamma;
    ddz    = (uzh*cdt_dz)*rgamma; // ddx,ddy,ddz are normalized displacement
    dxh    = dx  + ddx;
    dyh    = dy  + ddy;
    dzh    = dz  + ddz;           // Half step position
    dx1    = dxh + ddx;
    dy1    = dyh + ddy;
    dz1    = dzh + ddz;           // New particle position
    outbnd = (dx1>one) | (dx1<neg_one) |
             (dy1>one) | (dy1<neg_one) |
             (dz1>one) | (dz1<neg_one);
    store_4x4_tr( merge(outbnd,dx,dx1), // Do not update outbnd particles
                  merge(outbnd,dy,dy1),
                  merge(outbnd,dz,dz1),
                  ii,
                  &p[0].dx, &p[1].dx, &p[2].dx, &p[3].dx );

    // Accumulate current of inbnd particles
    // Note: accumulator values are 4 times the total physical charge that
    // passed through the appropriate current quadrant in a time-step

    qa  = czero(outbnd,q);          // Do not accumulate outbnd particles
    /**/                            // The half step position is not the
    /**/                            // streak midpoint for them!
    ccc = qa*ddx*ddy*ddz*one_third; // Charge conservation correction

#   define ACCUMULATE_J(X,Y,Z)                                          \
    a4##X  = qa*dd##X;      /* a4 = qa ddx                                */ \
    a1##X  = a4##X*d##Y##h; /* a1 = qa ddx dyh                            */ \
    a0##X  = a4##X-a1##X;   /* a0 = qa ddx (1-dyh)                        */ \
    a1##X += a4##X;         /* a1 = qa ddx (1+dyh)                        */ \
    a4##X  = one+d##Z##h;   /* a4 = 1+dzh                                 */ \
    a2##X  = a0##X*a4##X;   /* a2 = qa ddx (1-dyh)(1+dzh)                 */ \
    a3##X  = a1##X*a4##X;   /* a3 = qa ddx (1+dyh)(1+dzh)                 */ \
    a4##X  = one-d##Z##h;   /* a4 = 1-dzh                                 */ \
    a0##X *= a4##X;         /* a0 = qa ddx (1-dyh)(1-dzh)                 */ \
    a1##X *= a4##X;         /* a1 = qa ddx (1+dyh)(1-dzh)                 */ \
    a0##X += ccc;           /* a0 = qa ddx [ (1-dyh)(1-dzh) + ddy*ddz/3 ] */ \
    a1##X -= ccc;           /* a1 = qa ddx [ (1+dyh)(1-dzh) - ddy*ddz/3 ] */ \
    a2##X -= ccc;           /* a2 = qa ddx [ (1-dyh)(1+dzh) - ddy*ddz/3 ] */ \
    a3##X += ccc;           /* a3 = qa ddx [ (1+dyh)(1+dzh) + ddy*ddz/3 ] */ \
    transpose(a0##X,a1##X,a2##X,a3##X)

    ACCUMULATE_J( x,y,z );
    ACCUMULATE_J( y,z,x );
    ACCUMULATE_J( z,x,y );

#   undef ACCUMULATE_J

    ja = PTR_ACCUMULATOR(ii0);
    increment_4x1( ja->jx, a0x );
    increment_4x1( ja->jy, a0y );
    increment_4x1( ja->jz, a0z );

    ja = PTR_ACCUMULATOR(ii1);
    increment_4x1( ja->jx, a1x );
    increment_4x1( ja->jy, a1y );
    increment_4x1( ja->jz, a1z );

    ja = PTR_ACCUMULATOR(ii2);
    increment_4x1( ja->jx, a2x );
    increment_4x1( ja->jy, a2y );
    increment_4x1( ja->jz, a2z );

    ja = PTR_ACCUMULATOR(ii3);
    increment_4x1( ja->jx, a3x );
    increment_4x1( ja->jy, a3y );
    increment_4x1( ja->jz, a3z );

    // Update position and accumulate outbnd

#   define MOVE_OUTBND(N)                                      \
    if( outbnd(N) ) {                                          \
      pm->dispx = ddx(N);                                      \
      pm->dispy = ddy(N);                                      \
      pm->dispz = ddz(N);                                      \
      pm->i     = idx + N;                                     \
      if( move_p_spu( p+N, pm ) ) pm++, nm++;                  \
    }

    MOVE_OUTBND(0);
    MOVE_OUTBND(1);
    MOVE_OUTBND(2);
    MOVE_OUTBND(3);

#   undef MOVE_OUTBND
    
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

  // Determine which quads of particle quads this pipeline processes

  DISTRIBUTE( args->np, 16, pipeline_rank, n_pipeline, next_idx, np );

  // Determine which movers are reserved for this pipeline
  // Movers (16 bytes) should be reserved for pipelines in at least
  // multiples of 8 such that the set of particle movers reserved for a
  // pipeline is 128-bit aligned and a multiple of 128-bits in size. 

  args->max_nm -= args->np&15; // Insure host gets enough
  if( args->max_nm<0 ) args->max_nm = 0;
  DISTRIBUTE( args->max_nm, 8, pipeline_rank, n_pipeline, itmp, seg->max_nm );
  seg->pm        = args->pm + itmp*sizeof(particle_mover_t);
  seg->nm        = 0;
  seg->n_ignored = 0;

  // Determine which accumulator array is reserved for this pipeline

  args->a0 += (1+pipeline_rank)*(args->nx+2)*(args->ny+2)*(args->nz+2)*
    sizeof(accumulator_t);
  
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
                                             np_block[buffer]>>2 )

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

  cache_flush( accumulator_cache );

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
