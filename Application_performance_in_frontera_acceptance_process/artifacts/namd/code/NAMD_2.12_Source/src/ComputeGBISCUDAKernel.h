#include <stdio.h>
#define GBIS_CUDA
#include "ComputeGBIS.inl"

#undef KEPLER_SHUFFLE
#ifdef __CUDA_ARCH__
#define KEPLER_SHUFFLE
#if __CUDA_ARCH__ < 300
#undef KEPLER_SHUFFLE
#endif
#endif

//1111111111111111111111111111111111111111111111111111111111
//
// GBIS Phase 1 CUDA Kernal
//
//1111111111111111111111111111111111111111111111111111111111
__global__ static void GBIS_P1_Kernel (
  const patch_pair *patch_pairs, // atoms pointers and such
  const atom *atoms,             // position & charge
  const float *intRad0,          // read in intrinsic radius
  const float *intRadS,          // read in intrinsic radius
  GBReal *tmp_psiSum,            // temporary device memory
  GBReal *psiSum,                // host-mapped memory
  const float a_cut,             // P1 interaction cutoff
  const float rho_0,             // H(i,j) parameter
  float3 lata,
  float3 latb,
  float3 latc,
  unsigned int *P1_counters
) {

  // shared memory
  __shared__ GBReal sh_psiSumJ_2d[NUM_WARP][WARPSIZE];
  __shared__ patch_pair sh_patch_pair;
#ifndef KEPLER_SHUFFLE
  __shared__ atom sh_jpq_2d[NUM_WARP][WARPSIZE];
  __shared__ float sh_intRad0j_2d[NUM_WARP][WARPSIZE];
#endif

  volatile GBReal* sh_psiSumJ = sh_psiSumJ_2d[threadIdx.y];
#ifndef KEPLER_SHUFFLE
  volatile atom* sh_jpq = sh_jpq_2d[threadIdx.y];
  volatile float* sh_intRad0j = sh_intRad0j_2d[threadIdx.y];
#endif

  // Load data into shared memory
  {
    const int t = threadIdx.x + threadIdx.y*WARPSIZE;
    if (t < PATCH_PAIR_SIZE) {
      int* src = (int *)&patch_pairs[blockIdx.x];
      int* dst = (int *)&sh_patch_pair;
      dst[t] = src[t];
    }
    __syncthreads();

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

  //iterate over chunks of atoms within Patch 1
  for (int blocki = threadIdx.y*WARPSIZE; blocki < sh_patch_pair.patch1_size; blocki += NUM_WARP*WARPSIZE) {

    int nloopi = sh_patch_pair.patch1_size - blocki;
    nloopi = min(nloopi, WARPSIZE);

    //this thread calculates only the force on atomi; iterating over js
    atom atomi;
    float intRadSi;
    int i;

    // load BLOCK of Patch i atoms
    if ( blocki + threadIdx.x < sh_patch_pair.patch1_size ) {
      i = sh_patch_pair.patch1_start + blocki + threadIdx.x;
      float4 tmpa = ((float4*)atoms)[i];
      atomi.position.x = tmpa.x + sh_patch_pair.offset.x;
      atomi.position.y = tmpa.y + sh_patch_pair.offset.y;
      atomi.position.z = tmpa.z + sh_patch_pair.offset.z;
      atomi.charge = intRad0[i]; // overwrite charge with radius
      intRadSi = intRadS[i];
    } // load patch 1

    //init intermediate variables
    GBReal psiSumI = 0.f; // each thread accumulating single psi

    const bool diag_patch_pair = (sh_patch_pair.patch1_start == sh_patch_pair.patch2_start) && 
    (sh_patch_pair.offset.x == 0.0f && sh_patch_pair.offset.y == 0.0f && sh_patch_pair.offset.z == 0.0f);
    int blockj = (diag_patch_pair) ? blocki : 0;

    //iterate over chunks of atoms within Patch 2
    for (; blockj < sh_patch_pair.patch2_size; blockj += WARPSIZE ) {

      int nloopj = min(sh_patch_pair.patch2_size - blockj, WARPSIZE);

#ifdef KEPLER_SHUFFLE
      float xj;
      float yj;
      float zj;
      float chargej;
      float intRad0j_val;
#endif

      //load smaller chunk of j atoms into shared memory: coordinates
      if (blockj + threadIdx.x < sh_patch_pair.patch2_size) {
        int j = sh_patch_pair.patch2_start + blockj + threadIdx.x;
        float4 tmpa = ((float4*)atoms)[j];
#ifdef KEPLER_SHUFFLE
        xj = tmpa.x;
        yj = tmpa.y;
        zj = tmpa.z;
        chargej = intRadS[j];
        intRad0j_val = intRad0[j];
#else
        sh_jpq[threadIdx.x].position.x = tmpa.x;
        sh_jpq[threadIdx.x].position.y = tmpa.y;
        sh_jpq[threadIdx.x].position.z = tmpa.z;
        sh_jpq[threadIdx.x].charge = intRadS[j];
        sh_intRad0j[threadIdx.x] = intRad0[j];
#endif
      }

      const bool diag_tile = diag_patch_pair && (blocki == blockj);
      const int modval = diag_tile ? 2*WARPSIZE : WARPSIZE;

      sh_psiSumJ[threadIdx.x] = 0.f;

      //each thread loop over shared atoms
      int t = diag_tile ? 1 : 0;
#ifdef KEPLER_SHUFFLE
      if (diag_tile) {
        xj = __shfl(xj, (threadIdx.x+1) & (WARPSIZE-1) );
        yj = __shfl(yj, (threadIdx.x+1) & (WARPSIZE-1) );
        zj = __shfl(zj, (threadIdx.x+1) & (WARPSIZE-1) );
        chargej = __shfl(chargej, (threadIdx.x+1) & (WARPSIZE-1) );
        intRad0j_val = __shfl(intRad0j_val, (threadIdx.x+1) & (WARPSIZE-1) );
      }
#endif

      for (; t < WARPSIZE; ++t ) {
        int j = (t + threadIdx.x) % modval;
#ifndef KEPLER_SHUFFLE
        float xj = sh_jpq[j].position.x;
        float yj = sh_jpq[j].position.y;
        float zj = sh_jpq[j].position.z;
        float chargej = sh_jpq[j].charge;
        float intRad0j_val = sh_intRad0j[j];
#endif
        if (j < nloopj && threadIdx.x < nloopi)
        {
          float dx = atomi.position.x - xj;
          float dy = atomi.position.y - yj;
          float dz = atomi.position.z - zj;
          float r2 = dx*dx + dy*dy + dz*dz;

          // within cutoff        different atoms
          if (r2 < (a_cut+FS_MAX)*(a_cut+FS_MAX) && r2 > 0.01f) {
            // calculate H(i,j) [and not H(j,i)]
            float r_i = 1.f / sqrt(r2);
            float r  = r2 * r_i;
            float hij;
            int dij;
            CalcH(r,r2,r_i,a_cut,atomi.charge,chargej,hij,dij);
            psiSumI += hij;
            float hji;
            int dji;
            CalcH(r,r2,r_i,a_cut,intRad0j_val,intRadSi,hji,dji);
            sh_psiSumJ[j] += hji;
          } // cutoff
        } // if (j < nloopj)
#ifdef KEPLER_SHUFFLE        
        xj = __shfl(xj, (threadIdx.x+1) & (WARPSIZE-1) );
        yj = __shfl(yj, (threadIdx.x+1) & (WARPSIZE-1) );
        zj = __shfl(zj, (threadIdx.x+1) & (WARPSIZE-1) );
        chargej = __shfl(chargej, (threadIdx.x+1) & (WARPSIZE-1) );
        intRad0j_val = __shfl(intRad0j_val, (threadIdx.x+1) & (WARPSIZE-1) );
#endif
      } // for t

      if (blockj + threadIdx.x < sh_patch_pair.patch2_size) {
        int i_out = sh_patch_pair.patch2_start + blockj + threadIdx.x;
        atomicAdd(&tmp_psiSum[i_out], sh_psiSumJ[threadIdx.x]);
      }

    } // for block j
    //psiSumI now contains contributions from all j in 2nd patch

    // write psiSum to global memory buffer; to be accumulated later
    if ( blocki + threadIdx.x < sh_patch_pair.patch1_size) {
      int i_out = sh_patch_pair.patch1_start + blocki + threadIdx.x;
      atomicAdd(&tmp_psiSum[i_out], psiSumI);
    }
  } // for block i

  { // start of force sum
    // make sure psiSums are visible in global memory
    __threadfence();
    __syncthreads();

    // Mark patch pair (patch1_ind, patch2_ind) as "done"
    int patch1_ind = sh_patch_pair.patch1_ind;
    int patch2_ind = sh_patch_pair.patch2_ind;

    if (threadIdx.x == 0 && threadIdx.y == 0) {
      sh_patch_pair.patch_done[0] = false;
      sh_patch_pair.patch_done[1] = false;

      unsigned int patch1_num_pairs = sh_patch_pair.patch1_num_pairs;
      int patch1_old = atomicInc(&P1_counters[patch1_ind], patch1_num_pairs-1);
      if (patch1_old+1 == patch1_num_pairs) sh_patch_pair.patch_done[0] = true;
      if (patch1_ind != patch2_ind) {
        unsigned int patch2_num_pairs = sh_patch_pair.patch2_num_pairs;
        int patch2_old = atomicInc(&P1_counters[patch2_ind], patch2_num_pairs-1);
        if (patch2_old+1 == patch2_num_pairs) sh_patch_pair.patch_done[1] = true;
      }
    }
    // sync threads so that patch1_done and patch2_done are visible to all threads
    __syncthreads();

    if (sh_patch_pair.patch_done[0]) {
      const int start = sh_patch_pair.patch1_start;
      for (int i=threadIdx.x+threadIdx.y*WARPSIZE;i < sh_patch_pair.patch1_size;i+=NUM_WARP*WARPSIZE) {
        psiSum[start+i] = tmp_psiSum[start+i];
      }
    }

    if (sh_patch_pair.patch_done[1]) {
      const int start = sh_patch_pair.patch2_start;
      for (int i=threadIdx.x+threadIdx.y*WARPSIZE;i < sh_patch_pair.patch2_size;i+=NUM_WARP*WARPSIZE) {
        psiSum[start+i] = tmp_psiSum[start+i];
      }
    }

    if (sh_patch_pair.patch_done[0] || sh_patch_pair.patch_done[1]) {
      // Make sure page-locked host forces are up-to-date
#if __CUDA_ARCH__ < 200
      __threadfence();
#else
      __threadfence_system();
#endif
    }

  } // end of force sum
} //GBIS_P1

//2222222222222222222222222222222222222222222222222222222222
//
// GBIS Phase 2 CUDA Kernal
//
//2222222222222222222222222222222222222222222222222222222222
__global__ static void GBIS_P2_Kernel (
  const patch_pair *patch_pairs,// atoms pointers and such
  const atom *atoms,            // position & charge
  const float *bornRad,         // read in Born radius
  GBReal *tmp_dEdaSum,          // temporary device memory
  GBReal *dEdaSum,              // host-mapped memory
  const float a_cut,            // P1 interaction cutoff
  const float r_cut,            // P1 interaction cutoff
  const float scaling,          // scale nonbonded
  const float kappa,
  const float smoothDist,       // use interaction cutoff smoothing?
  const float epsilon_p,        // protein dielectric
  const float epsilon_s,        // solvent dielectric
  float3 lata,
  float3 latb,
  float3 latc,
  const int doEnergy,           // calculate energy too?
  const int doFullElec,         // calc dEdaSum for P3 full electrostatics
  float4 *tmp_forces,           // temporary  device memory
  float4 *forces,               // host-mapped memory
  float *tmp_energy,            // temporary device memory
  float *energy,                // host-mapped memory
  unsigned int *P2_counters
) {

  // Shared memory
  __shared__ patch_pair sh_patch_pair;
#ifndef KEPLER_SHUFFLE
  __shared__ atom sh_jpq_2d[NUM_WARP][WARPSIZE];
  __shared__ float sh_jBornRad_2d[NUM_WARP][WARPSIZE];
#endif
  __shared__ float3 sh_forceJ_2d[NUM_WARP][WARPSIZE];
  __shared__ float sh_dEdaSumJ_2d[NUM_WARP][WARPSIZE];

#ifndef KEPLER_SHUFFLE
  volatile atom* sh_jpq = sh_jpq_2d[threadIdx.y];
  volatile float* sh_jBornRad = sh_jBornRad_2d[threadIdx.y];
#endif
  volatile float3* sh_forceJ = sh_forceJ_2d[threadIdx.y];
  volatile float* sh_dEdaSumJ = sh_dEdaSumJ_2d[threadIdx.y];

  // Load data into shared memory
  {
    const int t = threadIdx.x + threadIdx.y*WARPSIZE;
    if (t < PATCH_PAIR_SIZE) {
      int* src = (int *)&patch_pairs[blockIdx.x];
      int* dst = (int *)&sh_patch_pair;
      dst[t] = src[t];
    }
    __syncthreads();

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

  float energyT = 0.f; // total energy for this thread; to be reduced

  //values used in loop
  float r_cut2 = r_cut*r_cut;
  float r_cut_2 = 1.f / r_cut2;
  float r_cut_4 = 4.f*r_cut_2*r_cut_2;
  float epsilon_s_i = 1.f / epsilon_s;
  float epsilon_p_i = 1.f / epsilon_p;

  //iterate over chunks of atoms within Patch 1
  for ( int blocki = threadIdx.y*WARPSIZE; blocki < sh_patch_pair.patch1_size; blocki += BLOCK_SIZE ) {

    int nloopi = sh_patch_pair.patch1_size - blocki;
    nloopi = min(nloopi, WARPSIZE);

    //this thread calculates only the force on atomi; iterating over js
    atom atomi;
    float bornRadI;
    int i;

    // load BLOCK of Patch i atoms
    if ( blocki + threadIdx.x < sh_patch_pair.patch1_size ) {
      i = sh_patch_pair.patch1_start + blocki + threadIdx.x;
      float4 tmpa = ((float4*)atoms)[i];
      atomi.position.x = tmpa.x + sh_patch_pair.offset.x;
      atomi.position.y = tmpa.y + sh_patch_pair.offset.y;
      atomi.position.z = tmpa.z + sh_patch_pair.offset.z;
      atomi.charge = - tmpa.w * scaling;
      bornRadI = bornRad[i];
    } // load patch 1

    //init intermediate variables
    GBReal dEdaSumI = 0.f; // each thread accumulating single psi
    float3 forceI;
    forceI.x = 0.f;
    forceI.y = 0.f;
    forceI.z = 0.f;

    const bool diag_patch_pair = (sh_patch_pair.patch1_start == sh_patch_pair.patch2_start) && 
    (sh_patch_pair.offset.x == 0.0f && sh_patch_pair.offset.y == 0.0f && sh_patch_pair.offset.z == 0.0f);
    int blockj = (diag_patch_pair) ? blocki : 0;
    //iterate over chunks of atoms within Patch 2
    for (; blockj < sh_patch_pair.patch2_size; blockj += WARPSIZE) {

      int nloopj = min(sh_patch_pair.patch2_size - blockj, WARPSIZE);

#ifdef KEPLER_SHUFFLE
      float xj;
      float yj;
      float zj;
      float chargej;
      float bornRadJ;
#endif

      //load smaller chunk of j atoms into shared memory: coordinates
      if (blockj + threadIdx.x < sh_patch_pair.patch2_size) {
        int j = sh_patch_pair.patch2_start + blockj + threadIdx.x;
        float4 tmpa = ((float4*)atoms)[j];
#ifdef KEPLER_SHUFFLE
        xj = tmpa.x;
        yj = tmpa.y;
        zj = tmpa.z;
        chargej = tmpa.w;
        bornRadJ = bornRad[j];
#else
        sh_jpq[threadIdx.x].position.x = tmpa.x;
        sh_jpq[threadIdx.x].position.y = tmpa.y;
        sh_jpq[threadIdx.x].position.z = tmpa.z;
        sh_jpq[threadIdx.x].charge = tmpa.w;
        sh_jBornRad[threadIdx.x] = bornRad[j];
#endif
      }

      sh_forceJ[threadIdx.x].x = 0.f;
      sh_forceJ[threadIdx.x].y = 0.f;
      sh_forceJ[threadIdx.x].z = 0.f;
      sh_dEdaSumJ[threadIdx.x] = 0.f;

      const bool diag_tile = diag_patch_pair && (blocki == blockj);
      const int modval = diag_tile ? 2*WARPSIZE : WARPSIZE;
      for (int t=0; t < WARPSIZE; ++t ) {
        int j = (t + threadIdx.x) % modval;
#ifndef KEPLER_SHUFFLE
        float xj = sh_jpq[j].position.x;
        float yj = sh_jpq[j].position.y;
        float zj = sh_jpq[j].position.z;
        float chargej = sh_jpq[j].charge;
        float bornRadJ = sh_jBornRad[j];
#endif
        if (j < nloopj && threadIdx.x < nloopi)
        {
          float dx = atomi.position.x - xj;
          float dy = atomi.position.y - yj;
          float dz = atomi.position.z - zj;
          float r2 = dx*dx + dy*dy + dz*dz;

          // within cutoff different atoms
          if (r2 < r_cut2 && r2 > 0.01f) {

            float r_i = 1.f / sqrt(r2);
            float r  = r2 * r_i;
            //float bornRadJ = sh_jBornRad[j];

            //calculate GB energy
            float qiqj = atomi.charge*chargej;
            float aiaj = bornRadI*bornRadJ;
            float aiaj4 = 4*aiaj;
            float expr2aiaj4 = exp(-r2/aiaj4);
            float fij = sqrt(r2+aiaj*expr2aiaj4);
            float f_i = 1/fij;
            float expkappa = exp(-kappa*fij);
            float Dij = epsilon_p_i - expkappa*epsilon_s_i;
            float gbEij = qiqj*Dij*f_i;

            //calculate energy derivatives
            float ddrfij = r*f_i*(1.f - 0.25f*expr2aiaj4);
            float ddrf_i = -ddrfij*f_i*f_i;
            float ddrDij = kappa*expkappa*ddrfij*epsilon_s_i;
            float ddrGbEij = qiqj*(ddrDij*f_i+Dij*ddrf_i);

            //NAMD smoothing function
            float scale = 1.f;
            float ddrScale = 0.f;
            float forcedEdr;
            if (smoothDist > 0.f) {
              scale = r2 * r_cut_2 - 1.f;
              scale *= scale;
              ddrScale = r*(r2-r_cut2)*r_cut_4;
              energyT += gbEij * scale;
              forcedEdr = -(ddrGbEij)*scale-(gbEij)*ddrScale;
            } else {
              energyT += gbEij;
              forcedEdr = -ddrGbEij;
            }

            //add dEda
            if (doFullElec) {
              float dEdai = 0.5f*qiqj*f_i*f_i
                        *(kappa*epsilon_s_i*expkappa-Dij*f_i)
                        *(aiaj+0.25f*r2)*expr2aiaj4/bornRadI*scale;//0
              dEdaSumI += dEdai;
              float dEdaj = 0.5f*qiqj*f_i*f_i
                        *(kappa*epsilon_s_i*expkappa-Dij*f_i)
                        *(aiaj+0.25f*r2)*expr2aiaj4/bornRadJ*scale;//0
              sh_dEdaSumJ[j] += dEdaj;
            }

            forcedEdr *= r_i;
            float tmpx = dx*forcedEdr;
            float tmpy = dy*forcedEdr;
            float tmpz = dz*forcedEdr;
            forceI.x += tmpx;
            forceI.y += tmpy;
            forceI.z += tmpz;
            sh_forceJ[j].x -= tmpx;
            sh_forceJ[j].y -= tmpy;
            sh_forceJ[j].z -= tmpz;
          } // within cutoff
          if (r2 < 0.01f) {
            // GB Self Energy
            if (doEnergy) {
              float fij = bornRadI;//inf
              float expkappa = exp(-kappa*fij);//0
              float Dij = epsilon_p_i - expkappa*epsilon_s_i;
              float gbEij = atomi.charge*(atomi.charge / (-scaling) )*Dij/fij;
              energyT += 0.5f*gbEij;
            }
          } //same atom or within cutoff
        } // if (j < nloopj)
#ifdef KEPLER_SHUFFLE
        xj = __shfl(xj, (threadIdx.x+1) & (WARPSIZE-1) );
        yj = __shfl(yj, (threadIdx.x+1) & (WARPSIZE-1) );
        zj = __shfl(zj, (threadIdx.x+1) & (WARPSIZE-1) );
        chargej = __shfl(chargej, (threadIdx.x+1) & (WARPSIZE-1) );
        bornRadJ = __shfl(bornRadJ, (threadIdx.x+1) & (WARPSIZE-1) );
#endif
      } // for t
      if ( blockj + threadIdx.x < sh_patch_pair.patch2_size) {
        int i_out = sh_patch_pair.patch2_start + blockj + threadIdx.x;
        atomicAdd(&tmp_dEdaSum[i_out], sh_dEdaSumJ[threadIdx.x]);
        atomicAdd(&tmp_forces[i_out].x, sh_forceJ[threadIdx.x].x);
        atomicAdd(&tmp_forces[i_out].y, sh_forceJ[threadIdx.x].y);
        atomicAdd(&tmp_forces[i_out].z, sh_forceJ[threadIdx.x].z);
      }
    } // for block j
    //psiSumI now contains contributions from all j in 2nd patch

    // write psiSum to global memory buffer; to be accumulated later
    if ( blocki + threadIdx.x < sh_patch_pair.patch1_size) {
      int i_out = sh_patch_pair.patch1_start + blocki + threadIdx.x;
      atomicAdd(&tmp_dEdaSum[i_out], dEdaSumI);
      atomicAdd(&tmp_forces[i_out].x, forceI.x);
      atomicAdd(&tmp_forces[i_out].y, forceI.y);
      atomicAdd(&tmp_forces[i_out].z, forceI.z);
    }
  } // for block i

  //Energy Reduction
  if (doEnergy) {
    // Do not have to sync here because each warp writes to the same
    // portion of sh_jforce_2d as in the above force computation loop
    volatile float* sh_energy = (float *)&sh_forceJ_2d[threadIdx.y][0].x;
    // Reduce within warps
    sh_energy[threadIdx.x] = (float)energyT;
    for (int d=1;d < WARPSIZE;d*=2) {
      int pos = threadIdx.x + d;
      float val = (pos < WARPSIZE) ? sh_energy[pos] : 0.0f;
      sh_energy[threadIdx.x] += val;
    }
    __syncthreads();
    // Reduce among warps
    if (threadIdx.x == 0 && threadIdx.y == 0) {
      float tot_energy = 0.0f;
#pragma unroll
      for (int i=0;i < NUM_WARP;++i) {
        tot_energy += ((float *)&sh_forceJ_2d[i][0].x)[0];
      }
      int patch1_ind = sh_patch_pair.patch1_ind;
      atomicAdd(&tmp_energy[patch1_ind], (float)tot_energy);
    }
  } //end Energy Reduction

  { // start of reduction
    // make sure tmp_forces and tmp_dEdaSum are visible in global memory
    __threadfence();
    __syncthreads();

    // Mark patch pair (patch1_ind, patch2_ind) as "done"
    int patch1_ind = sh_patch_pair.patch1_ind;
    int patch2_ind = sh_patch_pair.patch2_ind;

    if (threadIdx.x == 0 && threadIdx.y == 0) {
      sh_patch_pair.patch_done[0] = false;
      sh_patch_pair.patch_done[1] = false;

      unsigned int patch1_num_pairs = sh_patch_pair.patch1_num_pairs;
      int patch1_old = atomicInc(&P2_counters[patch1_ind], patch1_num_pairs-1);
      if (patch1_old+1 == patch1_num_pairs) sh_patch_pair.patch_done[0] = true;
      if (patch1_ind != patch2_ind) {
        unsigned int patch2_num_pairs = sh_patch_pair.patch2_num_pairs;
        int patch2_old = atomicInc(&P2_counters[patch2_ind], patch2_num_pairs-1);
        if (patch2_old+1 == patch2_num_pairs) sh_patch_pair.patch_done[1] = true;
      }
    }
    // sync threads so that patch1_done and patch2_done are visible to all threads
    __syncthreads();

    if (sh_patch_pair.patch_done[0]) {
      const int start = sh_patch_pair.patch1_start;
      for (int i=threadIdx.x+threadIdx.y*WARPSIZE;i < sh_patch_pair.patch1_size;i+=NUM_WARP*WARPSIZE) {
        forces[start+i] = tmp_forces[start+i];
        dEdaSum[start+i] = tmp_dEdaSum[start+i];
      }
      energy[patch1_ind] = tmp_energy[patch1_ind];
    }

    if (sh_patch_pair.patch_done[1]) {
      const int start = sh_patch_pair.patch2_start;
      for (int i=threadIdx.x+threadIdx.y*WARPSIZE;i < sh_patch_pair.patch2_size;i+=NUM_WARP*WARPSIZE) {
        forces[start+i] = tmp_forces[start+i];
        dEdaSum[start+i] = tmp_dEdaSum[start+i];
      }
      energy[patch2_ind] = tmp_energy[patch2_ind];
    }

    if (sh_patch_pair.patch_done[0] || sh_patch_pair.patch_done[1]) {
      // Make sure page-locked host arrays are up-to-date
#if __CUDA_ARCH__ < 200
      __threadfence();
#else
      __threadfence_system();
#endif
    }

  } // end of sum
} //GBIS_P2

//3333333333333333333333333333333333333333333333333333333333
//
// GBIS Phase 3 CUDA Kernal
//
//3333333333333333333333333333333333333333333333333333333333
__global__ static void GBIS_P3_Kernel (
  const patch_pair *patch_pairs,  // atoms pointers and such
  const atom *atoms,              // position & charge
  const float *intRad0,           // read in intrinsic radius
  const float *intRadS,           // read in intrinsic radius
  const float *dHdrPrefix,        // read in prefix
  const float a_cut,              // P1 interaction cutoff
  const float rho_0,              // H(i,j) parameter
  const float scaling,            // scale nonbonded
  float3 lata,
  float3 latb,
  float3 latc,
  float4 *tmp_forces,             // temporary device memory
  float4 *forces,                 // host-mapped memory
  unsigned int *P3_counters
) {

  // Shared memory
  __shared__ patch_pair sh_patch_pair;
#ifndef KEPLER_SHUFFLE
  __shared__ atom sh_jpq_2d[NUM_WARP][WARPSIZE];
  __shared__ float sh_intRadJ0_2d[NUM_WARP][WARPSIZE];
  __shared__ float sh_jDHdrPrefix_2d[NUM_WARP][WARPSIZE];
#endif
  __shared__ float3 sh_forceJ_2d[NUM_WARP][WARPSIZE];

#ifndef KEPLER_SHUFFLE
  volatile atom* sh_jpq = sh_jpq_2d[threadIdx.y];
  volatile float* sh_intRadJ0 = sh_intRadJ0_2d[threadIdx.y];
  volatile float* sh_jDHdrPrefix = sh_jDHdrPrefix_2d[threadIdx.y];
#endif
  volatile float3* sh_forceJ = sh_forceJ_2d[threadIdx.y];

  // Load data into shared memory
  {
    const int t = threadIdx.x + threadIdx.y*WARPSIZE;
    if (t < PATCH_PAIR_SIZE) {
      int* src = (int *)&patch_pairs[blockIdx.x];
      int* dst = (int *)&sh_patch_pair;
      dst[t] = src[t];
    }
    __syncthreads();

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

  //iterate over chunks of atoms within Patch 1
  for ( int blocki = threadIdx.y*WARPSIZE; blocki < sh_patch_pair.patch1_size; blocki += NUM_WARP*WARPSIZE ) {

    int nloopi = sh_patch_pair.patch1_size - blocki;
    nloopi = min(nloopi, WARPSIZE);

    //this thread calculates only the force on atomi; iterating over js
    atom atomi;
    float intRadIS;
    int i;
    float dHdrPrefixI;

    // load BLOCK of Patch i atoms
    if ( blocki + threadIdx.x < sh_patch_pair.patch1_size ) {
      i = sh_patch_pair.patch1_start + blocki + threadIdx.x;
      float4 tmpa = ((float4*)atoms)[i];
      atomi.position.x = tmpa.x + sh_patch_pair.offset.x;
      atomi.position.y = tmpa.y + sh_patch_pair.offset.y;
      atomi.position.z = tmpa.z + sh_patch_pair.offset.z;
      atomi.charge = intRad0[i]; // overwrite charge with radius
      intRadIS = intRadS[i];
      dHdrPrefixI = dHdrPrefix[i];
    } // load patch 1

    //init intermediate variables
    float3 forceI;
    forceI.x = 0.f;
    forceI.y = 0.f;
    forceI.z = 0.f;

    const bool diag_patch_pair = (sh_patch_pair.patch1_start == sh_patch_pair.patch2_start) && 
    (sh_patch_pair.offset.x == 0.0f && sh_patch_pair.offset.y == 0.0f && sh_patch_pair.offset.z == 0.0f);
    int blockj = (diag_patch_pair) ? blocki : 0;
    //iterate over chunks of atoms within Patch 2
    for (; blockj < sh_patch_pair.patch2_size; blockj += WARPSIZE ) {

      int nloopj = min(sh_patch_pair.patch2_size - blockj, WARPSIZE);

#ifdef KEPLER_SHUFFLE
      float xj;
      float yj;
      float zj;
      float intRadSJ;
      float dHdrPrefixJ;
      float intRadJ0;
#endif

      //load smaller chunk of j atoms into shared memory: coordinates
      if (blockj + threadIdx.x < sh_patch_pair.patch2_size) {
        int j = sh_patch_pair.patch2_start + blockj + threadIdx.x;
        float4 tmpa = ((float4*)atoms)[j];
#ifdef KEPLER_SHUFFLE
        xj = tmpa.x;
        yj = tmpa.y;
        zj = tmpa.z;
        intRadSJ = intRadS[j];
        dHdrPrefixJ = dHdrPrefix[j];
        intRadJ0 = intRad0[j];
#else
        sh_jpq[threadIdx.x].position.x = tmpa.x;
        sh_jpq[threadIdx.x].position.y = tmpa.y;
        sh_jpq[threadIdx.x].position.z = tmpa.z;
        sh_jpq[threadIdx.x].charge = intRadS[j];
        sh_jDHdrPrefix[threadIdx.x] = dHdrPrefix[j]; // load dHdrPrefix into shared
        sh_intRadJ0[threadIdx.x] = intRad0[j];
#endif
      }

      sh_forceJ[threadIdx.x].x = 0.f;
      sh_forceJ[threadIdx.x].y = 0.f;
      sh_forceJ[threadIdx.x].z = 0.f;

      const bool diag_tile = diag_patch_pair && (blocki == blockj);
      const int modval = diag_tile ? 2*WARPSIZE : WARPSIZE;
#ifdef KEPLER_SHUFFLE
      if (diag_tile) {
        xj = __shfl(xj, (threadIdx.x+1) & (WARPSIZE-1) );
        yj = __shfl(yj, (threadIdx.x+1) & (WARPSIZE-1) );
        zj = __shfl(zj, (threadIdx.x+1) & (WARPSIZE-1) );
        intRadSJ = __shfl(intRadSJ, (threadIdx.x+1) & (WARPSIZE-1) );
        dHdrPrefixJ = __shfl(dHdrPrefixJ, (threadIdx.x+1) & (WARPSIZE-1) );
        intRadJ0 = __shfl(intRadJ0, (threadIdx.x+1) & (WARPSIZE-1) );
      }
#endif
      int t = diag_tile ? 1 : 0;
      for (; t < WARPSIZE; ++t ) {
        int j = (t + threadIdx.x) % modval;
#ifndef KEPLER_SHUFFLE
        float xj = sh_jpq[j].position.x;
        float yj = sh_jpq[j].position.y;
        float zj = sh_jpq[j].position.z;
        float intRadSJ = sh_jpq[j].charge;
        float dHdrPrefixJ = sh_jDHdrPrefix[j];
        float intRadJ0 = sh_intRadJ0[j];
#endif
        if (j < nloopj && threadIdx.x < nloopi)
        {
          float dx = atomi.position.x - xj;
          float dy = atomi.position.y - yj;
          float dz = atomi.position.z - zj;
          float r2 = dx*dx + dy*dy + dz*dz;

          // within cutoff        different atoms
          if (r2 < (a_cut+FS_MAX)*(a_cut+FS_MAX) && r2 > 0.01f) {

            float r_i = 1.f / sqrt(r2);
            float r  = r2 * r_i;
            float dhij, dhji;
            int dij, dji;
            CalcDH(r,r2,r_i,a_cut,atomi.charge,intRadSJ,dhij,dij);
            CalcDH(r,r2,r_i,a_cut,intRadJ0,intRadIS,dhji,dji);

            float forceAlpha = -r_i*(dHdrPrefixI*dhij+dHdrPrefixJ*dhji);
            float tmpx = dx * forceAlpha;
            float tmpy = dy * forceAlpha;
            float tmpz = dz * forceAlpha;
            forceI.x += tmpx;
            forceI.y += tmpy;
            forceI.z += tmpz;
            sh_forceJ[j].x -= tmpx;
            sh_forceJ[j].y -= tmpy;
            sh_forceJ[j].z -= tmpz;
          } // cutoff
        } // if (j < nloopj...)
#ifdef KEPLER_SHUFFLE
        xj = __shfl(xj, (threadIdx.x+1) & (WARPSIZE-1) );
        yj = __shfl(yj, (threadIdx.x+1) & (WARPSIZE-1) );
        zj = __shfl(zj, (threadIdx.x+1) & (WARPSIZE-1) );
        intRadSJ = __shfl(intRadSJ, (threadIdx.x+1) & (WARPSIZE-1) );
        dHdrPrefixJ = __shfl(dHdrPrefixJ, (threadIdx.x+1) & (WARPSIZE-1) );
        intRadJ0 = __shfl(intRadJ0, (threadIdx.x+1) & (WARPSIZE-1) );
#endif
      } // for t
      if ( blockj + threadIdx.x < sh_patch_pair.patch2_size) {
        int i_out = sh_patch_pair.patch2_start + blockj + threadIdx.x;
        atomicAdd(&tmp_forces[i_out].x, sh_forceJ[threadIdx.x].x);
        atomicAdd(&tmp_forces[i_out].y, sh_forceJ[threadIdx.x].y);
        atomicAdd(&tmp_forces[i_out].z, sh_forceJ[threadIdx.x].z);
      }
    } // for block j

    // write psiSum to global memory buffer; to be accumulated later
    if ( blocki + threadIdx.x < sh_patch_pair.patch1_size) {
      int i_out = sh_patch_pair.patch1_start + blocki + threadIdx.x;
      atomicAdd(&tmp_forces[i_out].x, forceI.x);
      atomicAdd(&tmp_forces[i_out].y, forceI.y);
      atomicAdd(&tmp_forces[i_out].z, forceI.z);
    }
  } // for block i

  { // start of force sum
    // make sure forces are visible in global memory
    __threadfence();
    __syncthreads();

    // Mark patch pair (patch1_ind, patch2_ind) as "done"
    int patch1_ind = sh_patch_pair.patch1_ind;
    int patch2_ind = sh_patch_pair.patch2_ind;

    if (threadIdx.x == 0 && threadIdx.y == 0) {
      sh_patch_pair.patch_done[0] = false;
      sh_patch_pair.patch_done[1] = false;

      unsigned int patch1_num_pairs = sh_patch_pair.patch1_num_pairs;
      int patch1_old = atomicInc(&P3_counters[patch1_ind], patch1_num_pairs-1);
      if (patch1_old+1 == patch1_num_pairs) sh_patch_pair.patch_done[0] = true;
      if (patch1_ind != patch2_ind) {
        unsigned int patch2_num_pairs = sh_patch_pair.patch2_num_pairs;
        int patch2_old = atomicInc(&P3_counters[patch2_ind], patch2_num_pairs-1);
        if (patch2_old+1 == patch2_num_pairs) sh_patch_pair.patch_done[1] = true;
      }
    }
    // sync threads so that patch1_done and patch2_done are visible to all threads
    __syncthreads();

    if (sh_patch_pair.patch_done[0]) {
      const int start = sh_patch_pair.patch1_start;
      for (int i=threadIdx.x+threadIdx.y*WARPSIZE;i < sh_patch_pair.patch1_size;i+=NUM_WARP*WARPSIZE) {
        forces[start+i] = tmp_forces[start+i];
      }
    }

    if (sh_patch_pair.patch_done[1]) {
      const int start = sh_patch_pair.patch2_start;
      for (int i=threadIdx.x+threadIdx.y*WARPSIZE;i < sh_patch_pair.patch2_size;i+=NUM_WARP*WARPSIZE) {
        forces[start+i] = tmp_forces[start+i];
      }
    }

    if (sh_patch_pair.patch_done[0] || sh_patch_pair.patch_done[1]) {
      // Make sure page-locked host arrays are up-to-date
#if __CUDA_ARCH__ < 200
      __threadfence();
#else
      __threadfence_system();
#endif
    }

 } // end of force sum

} //GBIS_P3

