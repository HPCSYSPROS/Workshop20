// FIXME: This function assumes that the accumlator ghost values are
// zero.  Further, assumes that the ghost values of jfx, jfy, jfz are
// meaningless.  This might be changed to a more robust but slightly
// slower implementation in the near future.

#define IN_sf_interface
#include "sf_interface_private.h"

#define f(x,y,z) f[INDEX_FORTRAN_3(x,y,z,0,nx+1,0,ny+1,0,nz+1)]
#define a(x,y,z) a[INDEX_FORTRAN_3(x,y,z,0,nx+1,0,ny+1,0,nz+1)]

void
unload_accumulator_pipeline( unload_accumulator_pipeline_args_t * args,
			     int pipeline_rank,
                             int n_pipeline ) {
  field_t             * ALIGNED(128) f = args->f;
  const accumulator_t * ALIGNED(128) a = args->a;
  const grid_t        *              g = args->g;
  
  const accumulator_t * ALIGNED(16) a0;
  const accumulator_t * ALIGNED(16) ax,  * ALIGNED(16) ay,  * ALIGNED(16) az;
  const accumulator_t * ALIGNED(16) ayz, * ALIGNED(16) azx, * ALIGNED(16) axy;
  field_t * ALIGNED(16) f0;
  int x, y, z, n_voxel;
  
  const int nx = g->nx;
  const int ny = g->ny;
  const int nz = g->nz;

  const float cx = 0.25*g->rdy*g->rdz/g->dt;
  const float cy = 0.25*g->rdz*g->rdx/g->dt;
  const float cz = 0.25*g->rdx*g->rdy/g->dt;

  // Process the voxels assigned to this pipeline
  
  n_voxel = distribute_voxels( 1,nx+1, 1,ny+1, 1,nz+1, 16,
                               pipeline_rank, n_pipeline,
                               &x, &y, &z );

# define LOAD_STENCIL()                                                 \
  f0  = &f(x,  y,  z  );                                                \
  a0  = &a(x,  y,  z  );                                                \
  ax  = &a(x-1,y,  z  ); ay  = &a(x,  y-1,z  ); az  = &a(x,  y,  z-1);  \
  ayz = &a(x,  y-1,z-1); azx = &a(x-1,y,  z-1); axy = &a(x-1,y-1,z  )

  LOAD_STENCIL();

  for( ; n_voxel; n_voxel-- ) {

    f0->jfx += cx*( a0->jx[0] + ay->jx[1] + az->jx[2] + ayz->jx[3] );
    f0->jfy += cy*( a0->jy[0] + az->jy[1] + ax->jy[2] + azx->jy[3] );
    f0->jfz += cz*( a0->jz[0] + ax->jz[1] + ay->jz[2] + axy->jz[3] );

    f0++; a0++; ax++; ay++; az++; ayz++; azx++; axy++;

    x++;
    if( x>nx+1 ) {
      x=1, y++;
      if( y>ny+1 ) y=1, z++;
      LOAD_STENCIL();
    }

  }

# undef LOAD_STENCIL

}

#if defined(CELL_PPU_BUILD) && defined(USE_CELL_SPUS) && \
    defined(HAS_SPU_PIPELINE)

#error "SPU version not hooked up yet!"

#elif defined(V4_ACCELERATION) && defined(HAS_V4_PIPELINE)

#error "V4 version not hooked up yet!"

#endif

void
unload_accumulator( field_t             * ALIGNED(128) f,
                    const accumulator_t * ALIGNED(128) a,
                    const grid_t        * g ) {
  unload_accumulator_pipeline_args_t args[1];

  if( f==NULL ) ERROR(("Bad field"));
  if( a==NULL ) ERROR(("Bad accumulator"));
  if( g==NULL ) ERROR(("Bad grid"));

# if 0 // Original non-pipelined version

  for( z=1; z<=nz+1; z++ ) {
    for( y=1; y<=ny+1; y++ ) {

      x   = 1;
      f0  = &f(x,  y,  z  );
      a0  = &a(x,  y,  z  );
      ax  = &a(x-1,y,  z  ); ay  = &a(x,  y-1,z  ); az  = &a(x,  y,  z-1);
      ayz = &a(x,  y-1,z-1); azx = &a(x-1,y,  z-1); axy = &a(x-1,y-1,z  );

      for( x=1; x<=nx+1; x++ ) {

        f0->jfx += cx*( a0->jx[0] + ay->jx[1] + az->jx[2] + ayz->jx[3] );
        f0->jfy += cy*( a0->jy[0] + az->jy[1] + ax->jy[2] + azx->jy[3] );
        f0->jfz += cz*( a0->jz[0] + ax->jz[1] + ay->jz[2] + axy->jz[3] );

        f0++; a0++; ax++; ay++; az++; ayz++; azx++; axy++;

      }
    }
  }

# endif

  args->f = f;
  args->a = a;
  args->g = g;

  EXEC_PIPELINES( unload_accumulator, args, 0 );
  WAIT_PIPELINES();
}
