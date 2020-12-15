
#ifdef NAMD_CUDA

#define NAME(X) SLOWNAME( X )

#undef SLOW
#undef SLOWNAME
#ifdef DO_SLOW
#define SLOW(X) X
#define SLOWNAME(X) ENERGYNAME( X ## _slow )
#else
#define SLOW(X)
#define SLOWNAME(X) ENERGYNAME( X )
#endif

#undef ENERGY
#undef ENERGYNAME
#ifdef DO_ENERGY
#define ENERGY(X) X
#define ENERGYNAME(X) PAIRLISTNAME( X ## _energy )
#else
#define ENERGY(X)
#define ENERGYNAME(X) PAIRLISTNAME( X )
#endif

#undef GENPAIRLIST
#undef USEPAIRLIST
#undef PAIRLISTNAME
#ifdef MAKE_PAIRLIST
#define GENPAIRLIST(X) X
#define USEPAIRLIST(X)
#define PAIRLISTNAME(X) LAST( X ## _pairlist )
#else
#define GENPAIRLIST(X)
#define USEPAIRLIST(X) X
#define PAIRLISTNAME(X) LAST( X )
#endif

#define LAST(X) X

#undef KEPLER_SHUFFLE
#ifdef __CUDA_ARCH__
#define KEPLER_SHUFFLE
#if __CUDA_ARCH__ < 300
#undef KEPLER_SHUFFLE
#endif
#endif

#undef REG_JFORCE
#ifdef KEPLER_SHUFFLE
#ifndef MAKE_PAIRLIST
#define REG_JFORCE
#endif
#endif

#ifdef KEPLER_SHUFFLE
__device__ __forceinline__ static void NAME(shfl_reduction)(
		float *reg,
		float *smem
		)
{
  *reg += __shfl_xor(*reg, 16, 32);
  *reg += __shfl_xor(*reg, 8, 32);
  *reg += __shfl_xor(*reg, 4, 32);
  *reg += __shfl_xor(*reg, 2, 32);
  *reg += __shfl_xor(*reg, 1, 32);

  if ( threadIdx.x % 32 == 0 ) {
      atomicAdd(smem,*reg);
  }
  return;
}
#endif /* KEPLER_SHUFFLE */

__device__ __forceinline__ 
static void NAME(finish_forces_virials)(const int start, const int size, const int patch_ind,
					const atom* atoms,
					volatile float* sh_buf,
#ifndef KEPLER_SHUFFLE
					volatile float* sh_slow_buf, volatile float* sh_vcc,
#endif
					float4* tmpforces, float4* slow_tmpforces,
					float4* forces, float4* slow_forces,
					float* tmpvirials, float* slow_tmpvirials, 
					float* virials, float* slow_virials);

//
// Reduces up to three variables into global memory locations dst[0], dst[1], dst[2]
//
template<typename T, int n, int sh_buf_size>
__device__ __forceinline__
void NAME(reduceVariables)(volatile T* sh_buf, T* dst, T val1, T val2, T val3) {
	// Sanity check
	cuda_static_assert(n > 0 && n <= NUM_WARP);
#ifdef KEPLER_SHUFFLE
  // Requires NUM_WARP*n*sizeof(float) shared memory
  cuda_static_assert(sh_buf_size >= NUM_WARP*n*sizeof(T));
  // Reduce within warp
  for (int i=WARPSIZE/2;i >= 1;i/=2) {
    if (n >= 1) val1 += __shfl_xor(val1, i);
    if (n >= 2) val2 += __shfl_xor(val2, i);
    if (n >= 3) val3 += __shfl_xor(val3, i);
  }
  if (threadIdx.x == 0) {
    if (n >= 1) sh_buf[threadIdx.y*n + 0] = val1;
    if (n >= 2) sh_buf[threadIdx.y*n + 1] = val2;
    if (n >= 3) sh_buf[threadIdx.y*n + 2] = val3;
  }
  __syncthreads();
  if (threadIdx.x < n && threadIdx.y == 0) {
    T finalval = (T)0;
#pragma unroll
    for (int i=0;i < NUM_WARP;++i) {
      finalval += sh_buf[i*n + threadIdx.x];
    }
    atomicAdd(&dst[threadIdx.x], finalval);
  }
#else // ! KEPLER_SHUFFLE
  // Requires NUM_WARP*n*WARPSIZE*sizeof(float) shared memory
  cuda_static_assert(sh_buf_size >= NUM_WARP*n*WARPSIZE*sizeof(T));
  volatile T* sh_bufy = &sh_buf[threadIdx.y*n*WARPSIZE];
  if (n >= 1) sh_bufy[threadIdx.x*n + 0] = val1;
  if (n >= 2) sh_bufy[threadIdx.x*n + 1] = val2;
  if (n >= 3) sh_bufy[threadIdx.x*n + 2] = val3;
  // Reducue within warp
  for (int d=1;d < WARPSIZE;d*=2) {
    int pos = threadIdx.x + d;
    T val1t, val2t, val3t;
    if (n >= 1) val1t = (pos < WARPSIZE) ? sh_bufy[pos*n + 0] : (T)0;
    if (n >= 2) val2t = (pos < WARPSIZE) ? sh_bufy[pos*n + 1] : (T)0;
    if (n >= 3) val3t = (pos < WARPSIZE) ? sh_bufy[pos*n + 2] : (T)0;
    if (n >= 1) sh_bufy[threadIdx.x*n + 0] += val1t;
    if (n >= 2) sh_bufy[threadIdx.x*n + 1] += val2t;
    if (n >= 3) sh_bufy[threadIdx.x*n + 2] += val3t;
  }
  __syncthreads();
  if (threadIdx.x < n && threadIdx.y == 0) {
    T finalval = (T)0;
#pragma unroll
    for (int i=0;i < NUM_WARP;++i) {
      finalval += sh_buf[i*n*WARPSIZE + threadIdx.x];
    }
    atomicAdd(&dst[threadIdx.x], finalval);
  }
#endif // KEPLER_SHUFFLE
}

//
// Called with 2d thread setting:
// x-threadblock = warpSize = 32
// y-threadblock = NUM_WARP = 4
//
__global__ static void
GENPAIRLIST(__launch_bounds__(NUM_WARP*WARPSIZE, 10) )
USEPAIRLIST(__launch_bounds__(NUM_WARP*WARPSIZE, 12) )
NAME(dev_nonbonded)
     (const patch_pair* patch_pairs,
      const atom *atoms, const atom_param *atom_params,
      const int* vdw_types, unsigned int* plist,
      float4 *tmpforces, float4 *slow_tmpforces,
      float4 *forces, float4 *slow_forces,
      float* tmpvirials, float* slow_tmpvirials,
      float* virials, float* slow_virials,
      unsigned int* global_counters, int* force_ready_queue,
      const unsigned int *overflow_exclusions,
      const int npatches,
      const int block_begin, const int total_block_count, int* block_order,
      exclmask* exclmasks, const int lj_table_size,
      const float3 lata, const float3 latb, const float3 latc,
      const float cutoff2, const float plcutoff2, const int doSlow) {

  // Local structure definitions
  GENPAIRLIST(struct vdw_index {
    int vdw_type;
    int index;
  };)

  // Shared memory
  __shared__ patch_pair sh_patch_pair;
#ifndef REG_JFORCE
  __shared__ float3 sh_jforce_2d[NUM_WARP][WARPSIZE];
  SLOW(__shared__ float3 sh_jforce_slow_2d[NUM_WARP][WARPSIZE];)
#endif
#ifndef KEPLER_SHUFFLE
  __shared__ atom sh_jpq_2d[NUM_WARP][WARPSIZE];
#endif
  __shared__ float3 sh_iforcesum[SLOW(NUM_WARP+) NUM_WARP];

  ENERGY(
	 float totalev = 0.f;
	 float totalee = 0.f;
	 SLOW( float totales = 0.f; )
	)

  GENPAIRLIST(int nexcluded=0;);

  {
#ifndef KEPLER_SHUFFLE
  GENPAIRLIST(__shared__ atom_param sh_jap_2d[NUM_WARP][WARPSIZE];)
  USEPAIRLIST(__shared__ int sh_jap_vdw_type_2d[NUM_WARP][WARPSIZE];)
#endif
  USEPAIRLIST(__shared__ int sh_plist_ind[NUM_WARP];
              __shared__ unsigned int sh_plist_val[NUM_WARP];);

  // Load patch_pair -data into shared memory
  {
    const int t = threadIdx.x + threadIdx.y*WARPSIZE;

    if (t < 3*(SLOW(NUM_WARP+) NUM_WARP)) {
      float *p = (float *)sh_iforcesum;
      p[threadIdx.x] = 0.0f;
    }

    if (t < PATCH_PAIR_SIZE) {
      int* src = (int *)&patch_pairs[block_begin + blockIdx.x];
      int* dst = (int *)&sh_patch_pair;
      dst[t] = src[t];
    }
    // Need to sync here to make sure sh_patch_pair is ready
    __syncthreads();

    // Initialize pairlist index to impossible value
    USEPAIRLIST(if (threadIdx.x == 0) sh_plist_ind[threadIdx.y] = -1;);

    // Initialize pair list to "no interactions"
    GENPAIRLIST({
	if (t < sh_patch_pair.plist_size)
	  plist[sh_patch_pair.plist_start + t] = 0;
      })

    // convert scaled offset with current lattice and write into shared memory
    if (t == 0) {
      float offx = sh_patch_pair.offset.x * lata.x
	+ sh_patch_pair.offset.y * latb.x
	+ sh_patch_pair.offset.z * latc.x;
      float offy = sh_patch_pair.offset.x * lata.y
	+ sh_patch_pair.offset.y * latb.y
	+ sh_patch_pair.offset.z * latc.y;
      float offz = sh_patch_pair.offset.x * lata.z
	+ sh_patch_pair.offset.y * latb.z
	+ sh_patch_pair.offset.z * latc.z;
      sh_patch_pair.offset.x = offx;
      sh_patch_pair.offset.y = offy;
      sh_patch_pair.offset.z = offz;
    }

    __syncthreads();
  }

  // Compute pointers to shared memory to avoid point computation later on
#ifndef REG_JFORCE
  volatile float3* sh_jforce      = &sh_jforce_2d[threadIdx.y][0];
  SLOW(volatile float3* sh_jforce_slow = &sh_jforce_slow_2d[threadIdx.y][0];)
#endif

#ifndef KEPLER_SHUFFLE
  atom* sh_jpq       = &sh_jpq_2d[threadIdx.y][0];
  GENPAIRLIST(atom_param* sh_jap = &sh_jap_2d[threadIdx.y][0];);
  USEPAIRLIST(int* sh_jap_vdw_type = &sh_jap_vdw_type_2d[threadIdx.y][0];);
#endif

  for (int blocki = threadIdx.y*WARPSIZE;blocki < sh_patch_pair.patch1_size;blocki += WARPSIZE*NUM_WARP) {

    atom ipq;
    GENPAIRLIST(vdw_index iap;);
    USEPAIRLIST(int iap_vdw_type;);
    // Load i atom data
    if (blocki + threadIdx.x < sh_patch_pair.patch1_size) {
      int i = sh_patch_pair.patch1_start + blocki + threadIdx.x;
      float4 tmpa = ((float4*)atoms)[i];
      ipq.position.x = tmpa.x + sh_patch_pair.offset.x;
      ipq.position.y = tmpa.y + sh_patch_pair.offset.y;
      ipq.position.z = tmpa.z + sh_patch_pair.offset.z;
      ipq.charge = tmpa.w;
      GENPAIRLIST(uint4 tmpap = ((uint4*)atom_params)[i];
		  iap.vdw_type = tmpap.x*lj_table_size;
		  iap.index = tmpap.y;);
      USEPAIRLIST(iap_vdw_type = vdw_types[i]*lj_table_size;);
    }

    // i-forces in registers
    float3 iforce;
    iforce.x = 0.0f;
    iforce.y = 0.0f;
    iforce.z = 0.0f;
    SLOW(float3 iforce_slow;
	 iforce_slow.x = 0.0f;
	 iforce_slow.y = 0.0f;
	 iforce_slow.z = 0.0f;)

		const bool diag_patch_pair = (sh_patch_pair.patch1_start == sh_patch_pair.patch2_start) && 
		(sh_patch_pair.offset.x == 0.0f && sh_patch_pair.offset.y == 0.0f && sh_patch_pair.offset.z == 0.0f);
    int blockj = (diag_patch_pair) ? blocki : 0;
    for (;blockj < sh_patch_pair.patch2_size;blockj += WARPSIZE) {

      USEPAIRLIST({
	  const int size2 = (sh_patch_pair.patch2_size-1)/WARPSIZE+1;
	  int pos = (blockj/WARPSIZE) + (blocki/WARPSIZE)*size2;
	  int plist_ind = pos/32;
	  unsigned int plist_bit = 1 << (pos % 32);
	  // Check if we need to load next entry in the pairlist
	  if (plist_ind != sh_plist_ind[threadIdx.y]) {
	  	sh_plist_val[threadIdx.y] = plist[sh_patch_pair.plist_start + plist_ind];
	  	sh_plist_ind[threadIdx.y] = plist_ind;
	  }
	  if ((sh_plist_val[threadIdx.y] & plist_bit) == 0) continue;
	})

      // Load j atom data
#ifdef KEPLER_SHUFFLE
      atom jpq;
      GENPAIRLIST(atom_param jap;);
      USEPAIRLIST(int jap_vdw_type;);
#endif

      GENPAIRLIST(
  		// Avoid calculating pairs of blocks where all atoms on both blocks are fixed
   		if (blocki >= sh_patch_pair.patch1_free_size && blockj >= sh_patch_pair.patch2_free_size) continue;
		  int nfreej = sh_patch_pair.patch2_free_size - blockj;
		  int nloopj = min(sh_patch_pair.patch2_size - blockj, WARPSIZE);
		  );

      //GENPAIRLIST(bool inside_plcutoff = false;)
      if (blockj + threadIdx.x < sh_patch_pair.patch2_size) {
	int j = sh_patch_pair.patch2_start + blockj + threadIdx.x;
	float4 tmpa = ((float4*)atoms)[j];
#ifdef KEPLER_SHUFFLE
	jpq.position.x = tmpa.x;
	jpq.position.y = tmpa.y;
	jpq.position.z = tmpa.z;
	jpq.charge = tmpa.w;
#else
	sh_jpq[threadIdx.x].position.x = tmpa.x;
	sh_jpq[threadIdx.x].position.y = tmpa.y;
	sh_jpq[threadIdx.x].position.z = tmpa.z;
	sh_jpq[threadIdx.x].charge = tmpa.w;
#endif

#ifdef KEPLER_SHUFFLE
	GENPAIRLIST(jap = atom_params[j];)
        USEPAIRLIST(jap_vdw_type = vdw_types[j];)
#else
	GENPAIRLIST(sh_jap[threadIdx.x] = atom_params[j];)
        USEPAIRLIST(sh_jap_vdw_type[threadIdx.x] = vdw_types[j];)
#endif
      }

      // j-forces in shared memory
#ifdef REG_JFORCE
      float3 jforce;
      jforce.x = 0.0f;
      jforce.y = 0.0f;
      jforce.z = 0.0f;
      SLOW(float3 jforce_slow;
	   jforce_slow.x = 0.0f;
	   jforce_slow.y = 0.0f;
	   jforce_slow.z = 0.0f;
	   );
#else
      sh_jforce[threadIdx.x].x = 0.0f;
      sh_jforce[threadIdx.x].y = 0.0f;
      sh_jforce[threadIdx.x].z = 0.0f;
      SLOW(sh_jforce_slow[threadIdx.x].x = 0.0f;
	   sh_jforce_slow[threadIdx.x].y = 0.0f;
	   sh_jforce_slow[threadIdx.x].z = 0.0f;)
#endif

      GENPAIRLIST(unsigned int excl = 0;)
      USEPAIRLIST(
		  const int size2 = (sh_patch_pair.patch2_size-1)/WARPSIZE+1;
		  const int pos = (blockj/WARPSIZE) + (blocki/WARPSIZE)*size2;
		  unsigned int excl = exclmasks[sh_patch_pair.exclmask_start+pos].excl[threadIdx.x];
		  );
      GENPAIRLIST(
		  int nloopi = sh_patch_pair.patch1_size - blocki;
		  if (nloopi > WARPSIZE) nloopi = WARPSIZE;
		  // NOTE: We must truncate nfreei to be non-negative number since we're comparing to threadIdx.x (unsigned int) later on
		  int nfreei = max(sh_patch_pair.patch1_free_size - blocki, 0);
		  )
      const bool diag_tile = diag_patch_pair && (blocki == blockj);
      // Loop through tile diagonals. Local tile indices are:
      // i = threadIdx.x % WARPSIZE = constant
      // j = (t + threadIdx.x) % WARPSIZE
      const int modval = (diag_tile) ? 2*WARPSIZE-1 : WARPSIZE-1;
      int t = (diag_tile) ? 1 : 0;
      if (diag_tile) {
	USEPAIRLIST(excl >>= 1;);
#ifdef KEPLER_SHUFFLE
	jpq.charge = __shfl(jpq.charge, (threadIdx.x+1) & (WARPSIZE-1) );
	USEPAIRLIST(jap_vdw_type = __shfl(jap_vdw_type, (threadIdx.x+1) & (WARPSIZE-1) ););
	GENPAIRLIST(jap.vdw_type     = __shfl(jap.vdw_type, (threadIdx.x+1) & (WARPSIZE-1) );
		    jap.index        = __shfl(jap.index, (threadIdx.x+1) & (WARPSIZE-1) );
		    jap.excl_maxdiff = __shfl(jap.excl_maxdiff, (threadIdx.x+1) & (WARPSIZE-1) );
		    jap.excl_index   = __shfl(jap.excl_index, (threadIdx.x+1) & (WARPSIZE-1) );
		    );
#endif
      }

      for (; t < WARPSIZE; ++t) {
      	USEPAIRLIST(if (__any(excl & 1)))
      	{
	GENPAIRLIST(excl >>= 1;);
	int j = (t + threadIdx.x) & modval;
#ifdef KEPLER_SHUFFLE
	float tmpx = __shfl(jpq.position.x,j) - ipq.position.x;
	float tmpy = __shfl(jpq.position.y,j) - ipq.position.y;
	float tmpz = __shfl(jpq.position.z,j) - ipq.position.z;
	GENPAIRLIST(
		    int j_vdw_type     = jap.vdw_type;
		    int j_index        = jap.index;
		    int j_excl_maxdiff = jap.excl_maxdiff;
		    int j_excl_index   = jap.excl_index;
		    );
	float j_charge = jpq.charge;
	USEPAIRLIST(
		    int j_vdw_type = jap_vdw_type;
		    );
#endif
	GENPAIRLIST(if (j < nloopj && threadIdx.x < nloopi && (j < nfreej || threadIdx.x < nfreei) ))
	  {

#ifndef KEPLER_SHUFFLE
	  float tmpx = sh_jpq[j].position.x - ipq.position.x;
	  float tmpy = sh_jpq[j].position.y - ipq.position.y;
	  float tmpz = sh_jpq[j].position.z - ipq.position.z;
	  GENPAIRLIST(
		      int j_vdw_type     = sh_jap[j].vdw_type;
		      int j_index        = sh_jap[j].index;
		      int j_excl_maxdiff = sh_jap[j].excl_maxdiff;
		      int j_excl_index   = sh_jap[j].excl_index;
		      );
	  float j_charge = sh_jpq[j].charge;
	  USEPAIRLIST(
		      int j_vdw_type = sh_jap_vdw_type[j];
		      );
#endif
	  float r2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;
	  GENPAIRLIST(if (r2 < plcutoff2))
	  USEPAIRLIST(if ((excl & 1) && r2 < cutoff2))
	    {
	    GENPAIRLIST(
			bool excluded = false;
			int indexdiff = (int)(iap.index) - j_index;
			if ( abs(indexdiff) <= j_excl_maxdiff) {
			    indexdiff += j_excl_index;
			    int indexword = ((unsigned int) indexdiff) >> 5;
			    //indexword = tex1Dfetch(tex_exclusions, indexword);
			    if ( indexword < MAX_CONST_EXCLUSIONS )
			      indexword = const_exclusions[indexword];
			    else {
			      indexword = overflow_exclusions[indexword];
			    }
			    excluded = ((indexword & (1<<(indexdiff&31))) != 0);
			    if (excluded) nexcluded++;
			}
			if (!excluded) excl |= 0x80000000;
			)
	    GENPAIRLIST(if ( ! excluded && r2 < cutoff2))
	    {
	    ENERGY( float rsqrtfr2; );
	    float4 fi = tex1D(force_table, ENERGY(rsqrtfr2 =) rsqrtf(r2));
	    ENERGY( float4 ei = tex1D(energy_table, rsqrtfr2); );
	    GENPAIRLIST(float2 ljab = tex1Dfetch(lj_table, j_vdw_type + iap.vdw_type););
	    USEPAIRLIST(float2 ljab = tex1Dfetch(lj_table, j_vdw_type + iap_vdw_type););

	      float f_slow = ipq.charge * j_charge;
	      float f = ljab.x * fi.z + ljab.y * fi.y + f_slow * fi.x;
	      ENERGY(
		     float ev = ljab.x * ei.z + ljab.y * ei.y;
		     float ee = f_slow * ei.x;
		     SLOW( float es = f_slow * ei.w; )
		     )
	      SLOW( f_slow *= fi.w; )
	      ENERGY(
		     totalev += ev;
		     totalee += ee;
		     SLOW( totales += es; )
		     )
	      float fx = tmpx * f;
	      float fy = tmpy * f;
	      float fz = tmpz * f;
	      iforce.x += fx;
	      iforce.y += fy;
	      iforce.z += fz;
#ifdef REG_JFORCE
	      jforce.x -= fx;
	      jforce.y -= fy;
	      jforce.z -= fz;
#else
	      sh_jforce[j].x -= fx;
	      sh_jforce[j].y -= fy;
	      sh_jforce[j].z -= fz;
#endif
	      SLOW(
		   float fx_slow = tmpx * f_slow;
		   float fy_slow = tmpy * f_slow;
		   float fz_slow = tmpz * f_slow;
		   iforce_slow.x += fx_slow;
		   iforce_slow.y += fy_slow;
		   iforce_slow.z += fz_slow;
		   )
#ifdef REG_JFORCE
	      SLOW(
		   jforce_slow.x -= fx_slow;
		   jforce_slow.y -= fy_slow;
		   jforce_slow.z -= fz_slow;
		   )
#else
	      SLOW(
		   sh_jforce_slow[j].x -= fx_slow;
		   sh_jforce_slow[j].y -= fy_slow;
		   sh_jforce_slow[j].z -= fz_slow;
		   )
#endif
	    }
	    } // cutoff
	} // if (j < nloopj...)
}
	USEPAIRLIST(excl >>= 1;);
#ifdef KEPLER_SHUFFLE
	jpq.charge = __shfl(jpq.charge, (threadIdx.x+1) & (WARPSIZE-1) );
	USEPAIRLIST(jap_vdw_type = __shfl(jap_vdw_type, (threadIdx.x+1) & (WARPSIZE-1) ););
	GENPAIRLIST(jap.vdw_type     = __shfl(jap.vdw_type, (threadIdx.x+1) & (WARPSIZE-1) );
		    jap.index        = __shfl(jap.index, (threadIdx.x+1) & (WARPSIZE-1) );
		    jap.excl_maxdiff = __shfl(jap.excl_maxdiff, (threadIdx.x+1) & (WARPSIZE-1) );
		    jap.excl_index   = __shfl(jap.excl_index, (threadIdx.x+1) & (WARPSIZE-1) );
		    );
#ifdef REG_JFORCE
	jforce.x = __shfl(jforce.x, (threadIdx.x+1)&(WARPSIZE-1));
	jforce.y = __shfl(jforce.y, (threadIdx.x+1)&(WARPSIZE-1));
	jforce.z = __shfl(jforce.z, (threadIdx.x+1)&(WARPSIZE-1));
	SLOW(
	     jforce_slow.x = __shfl(jforce_slow.x, (threadIdx.x+1)&(WARPSIZE-1));
	     jforce_slow.y = __shfl(jforce_slow.y, (threadIdx.x+1)&(WARPSIZE-1));
	     jforce_slow.z = __shfl(jforce_slow.z, (threadIdx.x+1)&(WARPSIZE-1));
	     );
#endif
#endif
      } // t

      // Write j-forces
      GENPAIRLIST(if (__any(excl != 0))) {
	if ( blockj + threadIdx.x < sh_patch_pair.patch2_size ) {
	  int jforce_pos = sh_patch_pair.patch2_start + blockj + threadIdx.x;
#ifdef REG_JFORCE
	  atomicAdd(&tmpforces[jforce_pos].x, jforce.x);
	  atomicAdd(&tmpforces[jforce_pos].y, jforce.y);
	  atomicAdd(&tmpforces[jforce_pos].z, jforce.z);
	  SLOW(atomicAdd(&slow_tmpforces[jforce_pos].x, jforce_slow.x);
	       atomicAdd(&slow_tmpforces[jforce_pos].y, jforce_slow.y);
	       atomicAdd(&slow_tmpforces[jforce_pos].z, jforce_slow.z););
#else
	  atomicAdd(&tmpforces[jforce_pos].x, sh_jforce[threadIdx.x].x);
	  atomicAdd(&tmpforces[jforce_pos].y, sh_jforce[threadIdx.x].y);
	  atomicAdd(&tmpforces[jforce_pos].z, sh_jforce[threadIdx.x].z);
	  SLOW(atomicAdd(&slow_tmpforces[jforce_pos].x, sh_jforce_slow[threadIdx.x].x);
	       atomicAdd(&slow_tmpforces[jforce_pos].y, sh_jforce_slow[threadIdx.x].y);
	       atomicAdd(&slow_tmpforces[jforce_pos].z, sh_jforce_slow[threadIdx.x].z););
#endif
	}

      GENPAIRLIST(
		  const int size2 = (sh_patch_pair.patch2_size-1)/WARPSIZE+1;
		  int pos = (blockj/WARPSIZE) + (blocki/WARPSIZE)*size2;
		  exclmasks[sh_patch_pair.exclmask_start+pos].excl[threadIdx.x] = excl;
		  if (threadIdx.x == 0) {
		    int plist_ind = pos/32;
		    unsigned int plist_bit = 1 << (pos % 32);
		    atomicOr(&plist[sh_patch_pair.plist_start + plist_ind], plist_bit);
		  }
		  );
      }

    } // for (blockj)

    // Write i-forces
    if (blocki + threadIdx.x < sh_patch_pair.patch1_size) {
      int iforce_pos = sh_patch_pair.patch1_start + blocki + threadIdx.x;
      atomicAdd(&tmpforces[iforce_pos].x, iforce.x);
      atomicAdd(&tmpforces[iforce_pos].y, iforce.y);
      atomicAdd(&tmpforces[iforce_pos].z, iforce.z);
      SLOW(atomicAdd(&slow_tmpforces[iforce_pos].x, iforce_slow.x);
	   atomicAdd(&slow_tmpforces[iforce_pos].y, iforce_slow.y);
	   atomicAdd(&slow_tmpforces[iforce_pos].z, iforce_slow.z););
    }
    // Accumulate total forces for virial (warp synchronous)
#ifdef KEPLER_SHUFFLE
    for (int i=WARPSIZE/2;i >= 1;i/=2) {
      iforce.x += __shfl_xor(iforce.x, i);
      iforce.y += __shfl_xor(iforce.y, i);
      iforce.z += __shfl_xor(iforce.z, i);
      SLOW(
	   iforce_slow.x += __shfl_xor(iforce_slow.x, i);
	   iforce_slow.y += __shfl_xor(iforce_slow.y, i);
	   iforce_slow.z += __shfl_xor(iforce_slow.z, i);
	   );
    }
    if (threadIdx.x == 0) {
      sh_iforcesum[threadIdx.y].x += iforce.x;
      sh_iforcesum[threadIdx.y].y += iforce.y;
      sh_iforcesum[threadIdx.y].z += iforce.z;
      SLOW(
	   sh_iforcesum[threadIdx.y+NUM_WARP].x += iforce_slow.x;
	   sh_iforcesum[threadIdx.y+NUM_WARP].y += iforce_slow.y;
	   sh_iforcesum[threadIdx.y+NUM_WARP].z += iforce_slow.z;
	   );
    }
#else
    sh_jforce[threadIdx.x].x = iforce.x;
    sh_jforce[threadIdx.x].y = iforce.y;
    sh_jforce[threadIdx.x].z = iforce.z;
    SLOW(
	 sh_jforce_slow[threadIdx.x].x = iforce_slow.x;
	 sh_jforce_slow[threadIdx.x].y = iforce_slow.y;
	 sh_jforce_slow[threadIdx.x].z = iforce_slow.z;
	 );
    for (int d=1;d < WARPSIZE;d*=2) {
      int pos = threadIdx.x + d;
      float valx = (pos < WARPSIZE) ? sh_jforce[pos].x : 0.0f;
      float valy = (pos < WARPSIZE) ? sh_jforce[pos].y : 0.0f;
      float valz = (pos < WARPSIZE) ? sh_jforce[pos].z : 0.0f;
      SLOW(
	   float slow_valx = (pos < WARPSIZE) ? sh_jforce_slow[pos].x : 0.0f;
	   float slow_valy = (pos < WARPSIZE) ? sh_jforce_slow[pos].y : 0.0f;
	   float slow_valz = (pos < WARPSIZE) ? sh_jforce_slow[pos].z : 0.0f;
	   );
      sh_jforce[threadIdx.x].x += valx;
      sh_jforce[threadIdx.x].y += valy;
      sh_jforce[threadIdx.x].z += valz;
      SLOW(
	   sh_jforce_slow[threadIdx.x].x += slow_valx;
	   sh_jforce_slow[threadIdx.x].y += slow_valy;
	   sh_jforce_slow[threadIdx.x].z += slow_valz;
	   );
    }
    if (threadIdx.x == 0) {
      sh_iforcesum[threadIdx.y].x += sh_jforce[threadIdx.x].x;
      sh_iforcesum[threadIdx.y].y += sh_jforce[threadIdx.x].y;
      sh_iforcesum[threadIdx.y].z += sh_jforce[threadIdx.x].z;
      SLOW(
	   sh_iforcesum[threadIdx.y+NUM_WARP].x += sh_jforce_slow[threadIdx.x].x;
	   sh_iforcesum[threadIdx.y+NUM_WARP].y += sh_jforce_slow[threadIdx.x].y;
	   sh_iforcesum[threadIdx.y+NUM_WARP].z += sh_jforce_slow[threadIdx.x].z;
	   );
    }
#endif

  } // for (blocki)

  }

  {

#ifdef REG_JFORCE
#undef SH_BUF_SIZE
#define SH_BUF_SIZE NUM_WARP*(SLOW(9)+9)*sizeof(float)
    __shared__ float sh_buf[NUM_WARP*(SLOW(9)+9)];
#else // ! REG_JFORCE
#undef SH_BUF_SIZE
#define SH_BUF_SIZE NUM_WARP*WARPSIZE*3*sizeof(float)
    volatile float* sh_buf = (float *)&sh_jforce_2d[0][0];
  	// Sync here to make sure we can write into shared memory (sh_jforce_2d)
   	__syncthreads();
#endif

#if ENERGY(1+)0
   	NAME(reduceVariables)<float, SLOW(1+)2, SH_BUF_SIZE>(sh_buf, &tmpvirials[sh_patch_pair.patch1_ind*16 + 9], totalev, totalee, SLOW(totales+)0.0f);
#endif

#if GENPAIRLIST(1+)0
   	ENERGY(__syncthreads());
    NAME(reduceVariables)<int, 1, SH_BUF_SIZE>((int *)sh_buf, (int *)&tmpvirials[sh_patch_pair.patch1_ind*16 + 12], nexcluded, 0, 0);
#endif

  // Virials
  __syncthreads();
  if (threadIdx.x < SLOW(3+)3 && threadIdx.y == 0) {
    float* sh_virials = (float *)sh_iforcesum + (threadIdx.x % 3) + (threadIdx.x/3)*3*NUM_WARP;
    float iforcesum = 0.0f;
#pragma unroll
    for (int i=0;i < 3*NUM_WARP;i+=3) iforcesum += sh_virials[i];
    float vx = iforcesum*sh_patch_pair.offset.x;
    float vy = iforcesum*sh_patch_pair.offset.y;
    float vz = iforcesum*sh_patch_pair.offset.z;
    sh_iforcesum[threadIdx.x].x = vx;
    sh_iforcesum[threadIdx.x].y = vy;
    sh_iforcesum[threadIdx.x].z = vz;
  }
  if (threadIdx.x < SLOW(9+)9 && threadIdx.y == 0) {
    // virials are in sh_virials[0...8] and slow virials in sh_virials[9...17]
    float* sh_virials = (float *)sh_iforcesum;
    int patch1_ind = sh_patch_pair.patch1_ind;
    float *dst = (threadIdx.x < 9) ? tmpvirials : slow_tmpvirials;
    atomicAdd(&dst[patch1_ind*16 + (threadIdx.x % 9)], sh_virials[threadIdx.x]);
  }

  // Make sure forces are up-to-date in device global memory
  __threadfence();
  __syncthreads();

  // Mark patch pair (patch1_ind, patch2_ind) as "done"
  int patch1_ind = sh_patch_pair.patch1_ind;
  int patch2_ind = sh_patch_pair.patch2_ind;
  if (threadIdx.x == 0 && threadIdx.y == 0) {
    sh_patch_pair.patch_done[0] = false;
    sh_patch_pair.patch_done[1] = false;
    //
    // global_counters[0]: force_ready_queue
    // global_counters[1]: block_order
    // global_counters[2...npatches+1]: number of pairs finished for patch i+2
    //
    unsigned int patch1_num_pairs = sh_patch_pair.patch1_num_pairs;
    int patch1_old = atomicInc(&global_counters[patch1_ind+2], patch1_num_pairs-1);
    if (patch1_old+1 == patch1_num_pairs) sh_patch_pair.patch_done[0] = true;
    if (patch1_ind != patch2_ind) {
      unsigned int patch2_num_pairs = sh_patch_pair.patch2_num_pairs;
      int patch2_old = atomicInc(&global_counters[patch2_ind+2], patch2_num_pairs-1);
      if (patch2_old+1 == patch2_num_pairs) sh_patch_pair.patch_done[1] = true;
    }
  }
  // sync threads so that patch1_done and patch2_done are visible to all threads
  __syncthreads();

  if (sh_patch_pair.patch_done[0]) {

// #ifndef REG_JFORCE
//     volatile float* sh_buf = (float *)&sh_jforce_2d[0][0];
// #endif
#ifndef KEPLER_SHUFFLE
    volatile float* sh_vcc = (volatile float*)&sh_jpq_2d[0][0];
    volatile float* sh_slow_buf = NULL;
    SLOW(sh_slow_buf = (volatile float*)&sh_jforce_slow_2d[0][0];)
#endif
    NAME(finish_forces_virials)(sh_patch_pair.patch1_start, sh_patch_pair.patch1_size,
				patch1_ind, atoms, sh_buf,
#ifndef KEPLER_SHUFFLE
				sh_slow_buf, sh_vcc,
#endif
				tmpforces, slow_tmpforces, forces, slow_forces,
				tmpvirials, slow_tmpvirials, virials, slow_virials);

  }

  if (sh_patch_pair.patch_done[1]) {
// #ifndef REG_JFORCE
//     volatile float* sh_buf = (float *)&sh_jforce_2d[0][0];
// #endif
#ifndef KEPLER_SHUFFLE
    volatile float* sh_vcc = (volatile float*)&sh_jpq_2d[0][0];
    volatile float* sh_slow_buf = NULL;
    SLOW(sh_slow_buf = (volatile float*)&sh_jforce_slow_2d[0][0];)
#endif
    NAME(finish_forces_virials)(sh_patch_pair.patch2_start, sh_patch_pair.patch2_size,
				patch2_ind, atoms, sh_buf,
#ifndef KEPLER_SHUFFLE
				sh_slow_buf, sh_vcc,
#endif
				tmpforces, slow_tmpforces, forces, slow_forces,
				tmpvirials, slow_tmpvirials, virials, slow_virials);
  }

  if (force_ready_queue != NULL && (sh_patch_pair.patch_done[0] || sh_patch_pair.patch_done[1])) {
  	// Make sure page-locked host forces are up-to-date
#if __CUDA_ARCH__ < 200
    __threadfence();
#else
    __threadfence_system();
#endif
    __syncthreads();
    // Add patch into "force_ready_queue"
    if (threadIdx.x == 0 && threadIdx.y == 0) {
    	if (sh_patch_pair.patch_done[0]) {
        int ind = atomicInc(&global_counters[0], npatches-1);
        force_ready_queue[ind] = patch1_ind;
      }
      if (sh_patch_pair.patch_done[1]) {
        int ind = atomicInc(&global_counters[0], npatches-1);
        force_ready_queue[ind] = patch2_ind;
      }
      // Make sure "force_ready_queue" is visible in page-locked host memory
#if __CUDA_ARCH__ < 200
      __threadfence();
#else
      __threadfence_system();
#endif
    }
  }

  if (threadIdx.x == 0 && threadIdx.y == 0 && block_order != NULL) {
    int old = atomicInc(&global_counters[1], total_block_count-1);
    block_order[old] = block_begin + blockIdx.x;
  }
  }

}

//
// Copy patch forces to their final resting place in page-locked host memory
// and reduce virials
//
__device__ __forceinline__ 
static void NAME(finish_forces_virials)(const int start, const int size, const int patch_ind,
					const atom* atoms,
					volatile float* sh_buf,
#ifndef KEPLER_SHUFFLE
					volatile float* sh_slow_buf, volatile float* sh_vcc,
#endif
					float4* tmpforces, float4* slow_tmpforces,
					float4* forces, float4* slow_forces,
					float* tmpvirials, float* slow_tmpvirials, 
					float* virials, float* slow_virials) {
  float vxx = 0.f;
  float vxy = 0.f;
  float vxz = 0.f;
  float vyx = 0.f;
  float vyy = 0.f;
  float vyz = 0.f;
  float vzx = 0.f;
  float vzy = 0.f;
  float vzz = 0.f;
  SLOW(
       float slow_vxx = 0.f;
       float slow_vxy = 0.f;
       float slow_vxz = 0.f;
       float slow_vyx = 0.f;
       float slow_vyy = 0.f;
       float slow_vyz = 0.f;
       float slow_vzx = 0.f;
       float slow_vzy = 0.f;
       float slow_vzz = 0.f;
       )
  for (int i=threadIdx.x+threadIdx.y*WARPSIZE;i < size;i+=NUM_WARP*WARPSIZE) {
    const int p = start+i;
    float4 f = tmpforces[p];
    forces[p] = f;
    float4 pos = ((float4*)atoms)[p];
    vxx += f.x * pos.x;
    vxy += f.x * pos.y;
    vxz += f.x * pos.z;
    vyx += f.y * pos.x;
    vyy += f.y * pos.y;
    vyz += f.y * pos.z;
    vzx += f.z * pos.x;
    vzy += f.z * pos.y;
    vzz += f.z * pos.z;
    SLOW(
	 float4 slow_f = slow_tmpforces[p];
	 slow_forces[p] = slow_f;
	 slow_vxx += slow_f.x * pos.x;
	 slow_vxy += slow_f.x * pos.y;
	 slow_vxz += slow_f.x * pos.z;
	 slow_vyx += slow_f.y * pos.x;
	 slow_vyy += slow_f.y * pos.y;
	 slow_vyz += slow_f.y * pos.z;
	 slow_vzx += slow_f.z * pos.x;
	 slow_vzy += slow_f.z * pos.y;
	 slow_vzz += slow_f.z * pos.z;
 	 )
  }
#ifdef KEPLER_SHUFFLE
  // Reduce within warps
  for (int i=WARPSIZE/2;i >= 1;i/=2) {
    vxx += __shfl_xor(vxx, i);
    vxy += __shfl_xor(vxy, i);
    vxz += __shfl_xor(vxz, i);
    vyx += __shfl_xor(vyx, i);
    vyy += __shfl_xor(vyy, i);
    vyz += __shfl_xor(vyz, i);
    vzx += __shfl_xor(vzx, i);
    vzy += __shfl_xor(vzy, i);
    vzz += __shfl_xor(vzz, i);
    SLOW(
	 slow_vxx += __shfl_xor(slow_vxx, i);
	 slow_vxy += __shfl_xor(slow_vxy, i);
	 slow_vxz += __shfl_xor(slow_vxz, i);
	 slow_vyx += __shfl_xor(slow_vyx, i);
	 slow_vyy += __shfl_xor(slow_vyy, i);
	 slow_vyz += __shfl_xor(slow_vyz, i);
	 slow_vzx += __shfl_xor(slow_vzx, i);
	 slow_vzy += __shfl_xor(slow_vzy, i);
	 slow_vzz += __shfl_xor(slow_vzz, i);
	 )
  }
  // Reduce between warps
  // Requires NUM_WARP*(SLOW(9)+9)*sizeof(float) amount of shared memory
  if (threadIdx.x == 0) {
    sh_buf[threadIdx.y*(SLOW(9)+9) + 0] = vxx;
    sh_buf[threadIdx.y*(SLOW(9)+9) + 1] = vxy;
    sh_buf[threadIdx.y*(SLOW(9)+9) + 2] = vxz;
    sh_buf[threadIdx.y*(SLOW(9)+9) + 3] = vyx;
    sh_buf[threadIdx.y*(SLOW(9)+9) + 4] = vyy;
    sh_buf[threadIdx.y*(SLOW(9)+9) + 5] = vyz;
    sh_buf[threadIdx.y*(SLOW(9)+9) + 6] = vzx;
    sh_buf[threadIdx.y*(SLOW(9)+9) + 7] = vzy;
    sh_buf[threadIdx.y*(SLOW(9)+9) + 8] = vzz;
    SLOW(
	 sh_buf[threadIdx.y*(SLOW(9)+9) + 9]  = slow_vxx;
	 sh_buf[threadIdx.y*(SLOW(9)+9) + 10] = slow_vxy;
	 sh_buf[threadIdx.y*(SLOW(9)+9) + 11] = slow_vxz;
	 sh_buf[threadIdx.y*(SLOW(9)+9) + 12] = slow_vyx;
	 sh_buf[threadIdx.y*(SLOW(9)+9) + 13] = slow_vyy;
	 sh_buf[threadIdx.y*(SLOW(9)+9) + 14] = slow_vyz;
	 sh_buf[threadIdx.y*(SLOW(9)+9) + 15] = slow_vzx;
	 sh_buf[threadIdx.y*(SLOW(9)+9) + 16] = slow_vzy;
	 sh_buf[threadIdx.y*(SLOW(9)+9) + 17] = slow_vzz;
	 )
  }
  __syncthreads();
  // Write final virials into global memory
  if (threadIdx.x < SLOW(9+)9 && threadIdx.y == 0) {
    float v = 0.0f;
#pragma unroll
    for (int i=0;i < NUM_WARP;i++) v += sh_buf[i*(SLOW(9)+9) + threadIdx.x];
    float* dst = (threadIdx.x < 9) ? virials : slow_virials;
    const float* src = (threadIdx.x < 9) ? tmpvirials : slow_tmpvirials;
    int pos = patch_ind*16 + (threadIdx.x % 9);
    dst[pos] = v + src[pos];
  }
#else // ! KEPLER_SHUFFLE
  // We have total of NUM_WARP*WARPSIZE*3 floats, reduce in sets of three
  // (NOTE: we do have more shared memory available, so this could be optimized further
  //        for pre-Kepler architectures.)
  const int t = threadIdx.x + threadIdx.y*WARPSIZE;
  volatile float* sh_v1 = &sh_buf[0];
  volatile float* sh_v2 = &sh_buf[NUM_WARP*WARPSIZE];
  volatile float* sh_v3 = &sh_buf[2*NUM_WARP*WARPSIZE];
  SLOW(
       volatile float* sh_slow_v1 = &sh_slow_buf[0];
       volatile float* sh_slow_v2 = &sh_slow_buf[NUM_WARP*WARPSIZE];
       volatile float* sh_slow_v3 = &sh_slow_buf[2*NUM_WARP*WARPSIZE];
       )

  // vxx, vxy, vxz
  sh_v1[t] = vxx;
  sh_v2[t] = vxy;
  sh_v3[t] = vxz;
  SLOW(
       sh_slow_v1[t] = slow_vxx;
       sh_slow_v2[t] = slow_vxy;
       sh_slow_v3[t] = slow_vxz;
       )
  for (int d=1;d < NUM_WARP*WARPSIZE;d*=2) {
    int pos = t + d;
    float v1 = (pos < NUM_WARP*WARPSIZE) ? sh_v1[pos] : 0.0f;
    float v2 = (pos < NUM_WARP*WARPSIZE) ? sh_v2[pos] : 0.0f;
    float v3 = (pos < NUM_WARP*WARPSIZE) ? sh_v3[pos] : 0.0f;
    SLOW(
	 float slow_v1 = (pos < NUM_WARP*WARPSIZE) ? sh_slow_v1[pos] : 0.0f;
	 float slow_v2 = (pos < NUM_WARP*WARPSIZE) ? sh_slow_v2[pos] : 0.0f;
	 float slow_v3 = (pos < NUM_WARP*WARPSIZE) ? sh_slow_v3[pos] : 0.0f;
	 )
    __syncthreads();
    sh_v1[t] += v1;
    sh_v2[t] += v2;
    sh_v3[t] += v3;
    SLOW(
	 sh_slow_v1[t] += slow_v1;
	 sh_slow_v2[t] += slow_v2;
	 sh_slow_v3[t] += slow_v3;
	 )
    __syncthreads();
  }
  if (threadIdx.x == 0 && threadIdx.y == 0) {
    sh_vcc[0] = sh_v1[0];
    sh_vcc[1] = sh_v2[0];
    sh_vcc[2] = sh_v3[0];
    SLOW(
	 sh_vcc[9+0] = sh_slow_v1[0];
	 sh_vcc[9+1] = sh_slow_v2[0];
	 sh_vcc[9+2] = sh_slow_v3[0];
	 )
  }
  // vyx, vyy, vyz
  sh_v1[t] = vyx;
  sh_v2[t] = vyy;
  sh_v3[t] = vyz;
  SLOW(
       sh_slow_v1[t] = slow_vyx;
       sh_slow_v2[t] = slow_vyy;
       sh_slow_v3[t] = slow_vyz;
       )
  for (int d=1;d < NUM_WARP*WARPSIZE;d*=2) {
    int pos = t + d;
    float v1 = (pos < NUM_WARP*WARPSIZE) ? sh_v1[pos] : 0.0f;
    float v2 = (pos < NUM_WARP*WARPSIZE) ? sh_v2[pos] : 0.0f;
    float v3 = (pos < NUM_WARP*WARPSIZE) ? sh_v3[pos] : 0.0f;
    SLOW(
	 float slow_v1 = (pos < NUM_WARP*WARPSIZE) ? sh_slow_v1[pos] : 0.0f;
	 float slow_v2 = (pos < NUM_WARP*WARPSIZE) ? sh_slow_v2[pos] : 0.0f;
	 float slow_v3 = (pos < NUM_WARP*WARPSIZE) ? sh_slow_v3[pos] : 0.0f;
	 )
    __syncthreads();
    sh_v1[t] += v1;
    sh_v2[t] += v2;
    sh_v3[t] += v3;
    SLOW(
	 sh_slow_v1[t] += slow_v1;
	 sh_slow_v2[t] += slow_v2;
	 sh_slow_v3[t] += slow_v3;
	 )
    __syncthreads();
  }
  if (threadIdx.x == 0 && threadIdx.y == 0) {
    sh_vcc[3] = sh_v1[0];
    sh_vcc[4] = sh_v2[0];
    sh_vcc[5] = sh_v3[0];
    SLOW(
	 sh_vcc[9+3] = sh_slow_v1[0];
	 sh_vcc[9+4] = sh_slow_v2[0];
	 sh_vcc[9+5] = sh_slow_v3[0];
	 )
  }
  // vzx, vzy, vzz
  sh_v1[t] = vzx;
  sh_v2[t] = vzy;
  sh_v3[t] = vzz;
  SLOW(
       sh_slow_v1[t] = slow_vzx;
       sh_slow_v2[t] = slow_vzy;
       sh_slow_v3[t] = slow_vzz;
       )
  for (int d=1;d < NUM_WARP*WARPSIZE;d*=2) {
    int pos = t + d;
    float v1 = (pos < NUM_WARP*WARPSIZE) ? sh_v1[pos] : 0.0f;
    float v2 = (pos < NUM_WARP*WARPSIZE) ? sh_v2[pos] : 0.0f;
    float v3 = (pos < NUM_WARP*WARPSIZE) ? sh_v3[pos] : 0.0f;
    SLOW(
	 float slow_v1 = (pos < NUM_WARP*WARPSIZE) ? sh_slow_v1[pos] : 0.0f;
	 float slow_v2 = (pos < NUM_WARP*WARPSIZE) ? sh_slow_v2[pos] : 0.0f;
	 float slow_v3 = (pos < NUM_WARP*WARPSIZE) ? sh_slow_v3[pos] : 0.0f;
	 )
    __syncthreads();
    sh_v1[t] += v1;
    sh_v2[t] += v2;
    sh_v3[t] += v3;
    SLOW(
	 sh_slow_v1[t] += slow_v1;
	 sh_slow_v2[t] += slow_v2;
	 sh_slow_v3[t] += slow_v3;
	 )
    __syncthreads();
  }
  if (threadIdx.x == 0 && threadIdx.y == 0) {
    sh_vcc[6] = sh_v1[0];
    sh_vcc[7] = sh_v2[0];
    sh_vcc[8] = sh_v3[0];
    SLOW(
	 sh_vcc[9+6] = sh_slow_v1[0];
	 sh_vcc[9+7] = sh_slow_v2[0];
	 sh_vcc[9+8] = sh_slow_v3[0];
	 )
  }
  // Write final virials and energies into global memory
  if (threadIdx.x < SLOW(9+)9 && threadIdx.y == 0) {
    float* dst = (threadIdx.x < 9) ? virials : slow_virials;
    const float* src = (threadIdx.x < 9) ? tmpvirials : slow_tmpvirials;
    int pos = patch_ind*16 + (threadIdx.x % 9);
    dst[pos] = sh_vcc[threadIdx.x] + src[pos];
  }
#endif // KEPLER_SHUFFLE
  ENERGY(
	 // Write final energies into global memory
	 if (threadIdx.x < 3 && threadIdx.y == 0) {
	   int pos = patch_ind*16 + 9 + threadIdx.x;
	   virials[pos] = tmpvirials[pos];
	 }
	 );
  GENPAIRLIST(
  	if (threadIdx.x == 0 && threadIdx.y == 0) {
	  	int pos = patch_ind*16 + 12;
	   	virials[pos] = tmpvirials[pos];  		
  	}
  	);
} 

#endif // NAMD_CUDA

