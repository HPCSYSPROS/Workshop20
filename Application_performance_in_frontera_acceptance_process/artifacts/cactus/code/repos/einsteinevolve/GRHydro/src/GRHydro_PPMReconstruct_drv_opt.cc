#include <cmath>
#include <algorithm>
#include <vector>
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

#include "SpaceMask.h"
#include "GRHydro_PPM_opt.h"
#include "GRHydro_Reconstruct_drv_impl.hh" 

#undef velx
#undef vely
#undef velz

using namespace std;

extern "C" CCTK_INT GRHydro_UseGeneralCoordinates(const cGH * cctkGH);

static inline void oned_slice(const int npt, const int nstart, const int nstride, CCTK_REAL const * const restrict vec3d, vector<CCTK_REAL>& vec1d)
{
  for (int i=0; i<npt; ++i)
     vec1d[i] = vec3d[nstart+i*nstride];
}

static inline void oned_unslice(const int npt, const int nstart, const int nstride, const vector<CCTK_REAL>& vec1d, CCTK_REAL * const restrict vec3d)
{
  for (int i=0; i<npt; ++i)
     vec3d[nstart+i*nstride]=vec1d[i];
}

/*
  Cases that must be considered:
  * basic hydro
  * hydro + temperature + ye
  * hydro + ye
  * basic mhd
  * mhd + temperature + ye 
  * mhd + ye 
  * mppm (not supported right now)
  * not supporting trivial_rp
  * with or without divergence cleaning
 */

extern "C"
void CCTK_FNAME(GRHydro_PPMReconstruct_drv_opt)(cGH const * const restrict & /*restrict*/ cctkGH)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  //Multipatch related pointers
  const CCTK_REAL * restrict vup;
  const CCTK_REAL * restrict Bprim;
  const CCTK_REAL * restrict g11;
  const CCTK_REAL * restrict g12;
  const CCTK_REAL * restrict g13;
  const CCTK_REAL * restrict g22;
  const CCTK_REAL * restrict g23;
  const CCTK_REAL * restrict g33;
  
  if(GRHydro_UseGeneralCoordinates(cctkGH)) {
    vup=lvel;
    Bprim=lBvec;
    g11=gaa; g12=gab; g13=gac;
    g22=gbb; g23=gbc;
    g33=gcc;
  } else {
    vup=vel;
    Bprim=Bvec;
    g11=gxx; g12=gxy; g13=gxz;
    g22=gyy; g23=gyz;
    g33=gzz;
  }

  const int nx=cctk_ash[0];
  const int ny=cctk_ash[1];
  const int nz=cctk_ash[2];
  const int nxy=nx*ny;
  const int nxyz=nxy*nz;

  const int nmax_xyz=MAX(nx,MAX(ny,nz));

  // bit masks and patterns for trivial Riemann Problem spacemask
  const int type_bitsx = SpaceMask_GetTypeBits("Hydro_RiemannProblemX");
  //  trivialx = SpaceMask_GetStateBits("Hydro_RiemannProblemX","trivial");
  const int not_trivialx = SpaceMask_GetStateBits("Hydro_RiemannProblemX",
                                                  "not_trivial");

  const int type_bitsy = SpaceMask_GetTypeBits("Hydro_RiemannProblemY");
  //  trivialy = SpaceMask_GetStateBits("Hydro_RiemannProblemY","trivial");
  const int not_trivialy = SpaceMask_GetStateBits("Hydro_RiemannProblemY",
                                                  "not_trivial");

  const int type_bitsz = SpaceMask_GetTypeBits("Hydro_RiemannProblemZ");
  //  trivialz = SpaceMask_GetStateBits("Hydro_RiemannProblemZ","trivial");
  const int not_trivialz = SpaceMask_GetStateBits("Hydro_RiemannProblemZ",
                                                  "not_trivial");
    
  //   if use_enhanced_ppm, allow old PPM on one level
  bool apply_enhanced_ppm = false;
  if (GRHydro_oppm_reflevel == (-1) ||
      *GRHydro_reflevel != GRHydro_oppm_reflevel)
      apply_enhanced_ppm = use_enhanced_ppm != 0;
  else
      apply_enhanced_ppm = false;
  
  
  typedef
  void (*ppm1d_cxx_ptr_t)(const int nx, \
                    const double dx,              \
                    const double* restrict rho,   \
                    const double* restrict velx,  \
                    const double* restrict vely,  \
                    const double* restrict velz,  \
                    const double* restrict eps,   \
                    const double* restrict press, \
                    const double* restrict temp,  \
                    const double* restrict ye,    \
                    const double* restrict Bvcx,  \
                    const double* restrict Bvcy,  \
                    const double* restrict Bvcz,  \
                    const double* restrict psidc, \
                    double* restrict rhominus,    \
                    double* restrict velxminus,   \
                    double* restrict velyminus,   \
                    double* restrict velzminus,   \
                    double* restrict epsminus,    \
                    double* restrict tempminus,   \
                    double* restrict yeminus,     \
                    double* restrict Bvcxminus,   \
                    double* restrict Bvcyminus,   \
                    double* restrict Bvczminus,   \
                    double* restrict psidcminus,  \
                    double* restrict rhoplus,     \
                    double* restrict velxplus,    \
                    double* restrict velyplus,    \
                    double* restrict velzplus,    \
                    double* restrict epsplus,     \
                    double* restrict tempplus,    \
                    double* restrict yeplus,      \
                    double* restrict Bvcxplus,    \
                    double* restrict Bvcyplus,    \
                    double* restrict Bvczplus,    \
                    double* restrict psidcplus);
  static ppm1d_cxx_ptr_t ppm1d_cxx_funcs[64];
  ppm1d_cxx_ptr_t ppm1d_cxx_func;

  #define MKINDEX(do_temp,do_ye,do_mhd,dc_flag,do_ppm_detect,do_eppm) ((do_eppm<<5)|(do_temp<<4)|(do_ye<<3)|(do_mhd<<2)|(dc_flag<<1)|(do_ppm_detect))
  #define SETPOINTER(do_temp,do_ye,do_mhd,dc_flag,do_ppm_detect,do_eppm) \
                     if (!do_eppm) \
                        ppm1d_cxx_funcs[MKINDEX(do_temp,do_ye,do_mhd,dc_flag,do_ppm_detect, do_eppm)] = \
                        GRHydro_ppm1d_cxx<do_temp,do_ye,do_mhd,dc_flag,do_ppm_detect>; \
                     else \
                        ppm1d_cxx_funcs[MKINDEX(do_temp,do_ye,do_mhd,dc_flag,do_ppm_detect, do_eppm)] = \
                        GRHydro_eppm1d_cxx<do_temp,do_ye,do_mhd,dc_flag,do_ppm_detect>;
                        
  //  --- original PPM
  // first stuff without MHD
  SETPOINTER(false,false,false,false,false,false);
  SETPOINTER(false,false,false,false,true,false);
  // do_temp can only be true of do_ye is also true
  SETPOINTER(true,true,false,false,false,false);
  SETPOINTER(true,true,false,false,true,false);
  // but do_ye can be true and do_temp can be false
  SETPOINTER(false,true,false,false,false,false);
  SETPOINTER(false,true,false,false,true,false);

  // with MHD, but dc_flag false
  SETPOINTER(false,false,true,false,false,false);
  SETPOINTER(false,false,true,false,true,false);
  // do_temp can only be true of do_ye is also true
  SETPOINTER(true,true,true,false,false,false);
  SETPOINTER(true,true,true,false,true,false);
  // but do_ye can be true and do_temp can be false
  SETPOINTER(false,true,true,false,false,false);
  SETPOINTER(false,true,true,false,true,false);

  // with MHD, dc_flag true
  SETPOINTER(false,false,true,true,false,false);
  SETPOINTER(false,false,true,true,true,false);
  // do_temp can only be true of do_ye is also true
  SETPOINTER(true,true,true,true,false,false);
  SETPOINTER(true,true,true,true,true,false);
  // but do_ye can be true and do_temp can be false
  SETPOINTER(false,true,true,true,false,false);
  SETPOINTER(false,true,true,true,true,false);

  // ---- enhanced PPM
  // first stuff without MHD
  SETPOINTER(false,false,false,false,false,true);
  SETPOINTER(false,false,false,false,true,true);
  // do_temp can only be true of do_ye is also true
  SETPOINTER(true,true,false,false,false,true);
  SETPOINTER(true,true,false,false,true,true);
  // but do_ye can be true and do_temp can be false
  SETPOINTER(false,true,false,false,false,true);
  SETPOINTER(false,true,false,false,true,true);

  // with MHD, but dc_flag false
  SETPOINTER(false,false,true,false,false,true);
  SETPOINTER(false,false,true,false,true,true);
  // do_temp can only be true of do_ye is also true
  SETPOINTER(true,true,true,false,false,true);
  SETPOINTER(true,true,true,false,true,true);
  // but do_ye can be true and do_temp can be false
  SETPOINTER(false,true,true,false,false,true);
  SETPOINTER(false,true,true,false,true,true);

  // with MHD, dc_flag true
  SETPOINTER(false,false,true,true,false,true);
  SETPOINTER(false,false,true,true,true,true);
  // do_temp can only be true of do_ye is also true
  SETPOINTER(true,true,true,true,false,true);
  SETPOINTER(true,true,true,true,true,true);
  // but do_ye can be true and do_temp can be false
  SETPOINTER(false,true,true,true,false,true);
  SETPOINTER(false,true,true,true,true,true);
  
  
  const bool do_temp = *evolve_temper;
  const bool do_ye = *evolve_Y_e;
  const bool do_mhd = *evolve_MHD;
  const bool dc_flag = clean_divergence;
  const bool do_ppm_detect = ppm_detect;
  const bool do_Avec = CCTK_Equals(Bvec_evolution_method,"GRHydro_Avec");
  ppm1d_cxx_func = ppm1d_cxx_funcs[MKINDEX(do_temp,do_ye,do_mhd,dc_flag,do_ppm_detect,apply_enhanced_ppm)];
  if(ppm1d_cxx_func == NULL) {
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                "Combination do_temp = %d do_ye = %d do_mhd = %d dc_flag = %d do_ppm_detect = %d do_eppm = %d is not supported for PPM. Please instantiate the respective template",
                do_temp, do_ye, do_mhd, dc_flag, do_ppm_detect, apply_enhanced_ppm);
  }

  //index = i + j*nx + k*nx*ny;

  vector<CCTK_REAL> Wvup(0);
  // allocate temp memory for Wv reconstruction
  if (reconstruct_Wv) {
     Wvup.resize(3*nxyz, 0);
     #pragma omp parallel for
     for (int i=0; i < nxyz; ++i) {
        Wvup[i]        = w_lorentz[i]*vup[i];
        Wvup[i+nxyz]   = w_lorentz[i]*vup[i+nxyz];
        Wvup[i+2*nxyz] = w_lorentz[i]*vup[i+2*nxyz];
     }
  }
  
// PPM starts:
  if (*flux_direction == 1) {
    #pragma omp parallel
    {
    vector<CCTK_REAL> rho1d(nmax_xyz), eps1d(nmax_xyz), press1d(nmax_xyz);
    vector<CCTK_REAL> velx1d(nmax_xyz), vely1d(nmax_xyz), velz1d(nmax_xyz);
    vector<CCTK_REAL> temp1d(nmax_xyz), ye1d(nmax_xyz);
    vector<CCTK_REAL> Bvcx1d(nmax_xyz), Bvcy1d(nmax_xyz), Bvcz1d(nmax_xyz);
    vector<CCTK_REAL> psidc1d(nmax_xyz);

    vector<CCTK_REAL> rho1d_plus(nmax_xyz), eps1d_plus(nmax_xyz), 
                      press1d_plus(nmax_xyz);
    vector<CCTK_REAL> velx1d_plus(nmax_xyz), vely1d_plus(nmax_xyz),
                      velz1d_plus(nmax_xyz);
    vector<CCTK_REAL> temp1d_plus(nmax_xyz), ye1d_plus(nmax_xyz);
    vector<CCTK_REAL> Bvcx1d_plus(nmax_xyz), Bvcy1d_plus(nmax_xyz),
                      Bvcz1d_plus(nmax_xyz);
    vector<CCTK_REAL> psidc1d_plus(nmax_xyz);

    vector<CCTK_REAL> rho1d_minus(nmax_xyz), eps1d_minus(nmax_xyz), 
                      press1d_minus(nmax_xyz);
    vector<CCTK_REAL> velx1d_minus(nmax_xyz), vely1d_minus(nmax_xyz), 
                      velz1d_minus(nmax_xyz);
    vector<CCTK_REAL> temp1d_minus(nmax_xyz), ye1d_minus(nmax_xyz);
    vector<CCTK_REAL> Bvcx1d_minus(nmax_xyz), Bvcy1d_minus(nmax_xyz),
                      Bvcz1d_minus(nmax_xyz);
    vector<CCTK_REAL> psidc1d_minus(nmax_xyz);

    #pragma omp for
    for(int k = GRHydro_stencil-1; k < nz - GRHydro_stencil + 1 + transport_constraints; ++k) {
      for(int j = GRHydro_stencil-1; j < ny - GRHydro_stencil + 1 + transport_constraints; ++j) {
        
        // TODO: there is no need to slice in the x direction, we can just point into the array

        //For reference, the slicing call is:
        //void oned_slice(int npt,int nstart,int nstride,CCTK_REAL* vec3d,CCTK_REAL* vec1d)
        const int nstart=nx*(j+k*ny);
        
        oned_slice(nx,nstart,1,rho,rho1d);
        
        // TODO: reconstruct_Wv is not yet templated!
        if (!reconstruct_Wv) {
          oned_slice(nx,nstart,       1,vup,velx1d);
          oned_slice(nx,nstart+nxyz,  1,vup,vely1d);
          oned_slice(nx,nstart+2*nxyz,1,vup,velz1d);
        } else {
          oned_slice(nx,nstart,       1,&Wvup.front(),velx1d);
          oned_slice(nx,nstart+nxyz,  1,&Wvup.front(),vely1d);
          oned_slice(nx,nstart+2*nxyz,1,&Wvup.front(),velz1d);
        }
        oned_slice(nx,nstart,1,eps,eps1d);
        oned_slice(nx,nstart,1,press,press1d);
        
        if(do_ye) {
          oned_slice(nx,nstart,1,Y_e,ye1d);
          if(do_temp){
            oned_slice(nx,nstart,1,temperature,temp1d);
          } 
        }
        
        if(do_mhd) {
          oned_slice(nx,nstart,       1,Bprim,Bvcx1d);
          oned_slice(nx,nstart+nxyz,  1,Bprim,Bvcy1d);
          oned_slice(nx,nstart+2*nxyz,1,Bprim,Bvcz1d);
          if(clean_divergence) {
            oned_slice(nx,nstart,1,psidc,psidc1d);
          }
          if (do_Avec)
            CCTK_WARN(0, "Someone needs to figure out which variables to pass for Avec");
        }
        
        ppm1d_cxx_func(cctk_lsh[0],CCTK_DELTA_SPACE(0),
           &rho1d[0],&velx1d[0],&vely1d[0],&velz1d[0],&eps1d[0],&press1d[0],&temp1d[0],
           &ye1d[0],&Bvcx1d[0],&Bvcy1d[0],&Bvcz1d[0],&psidc1d[0],
           &rho1d_minus[0],&velx1d_minus[0],&vely1d_minus[0],&velz1d_minus[0],&eps1d_minus[0],&temp1d_minus[0],
           &ye1d_minus[0],&Bvcx1d_minus[0],&Bvcy1d_minus[0],&Bvcz1d_minus[0],&psidc1d_minus[0],
           &rho1d_plus[0],&velx1d_plus[0],&vely1d_plus[0],&velz1d_plus[0],&eps1d_plus[0],&temp1d_plus[0],
           &ye1d_plus[0],&Bvcx1d_plus[0],&Bvcy1d_plus[0],&Bvcz1d_plus[0],&psidc1d_plus[0]);
        
        
        oned_unslice(nx,nstart,1,rho1d_minus,rhominus);
        oned_unslice(nx,nstart,1,rho1d_plus,rhoplus);
        
        oned_unslice(nx,nstart,1,velx1d_minus,velxminus);
        oned_unslice(nx,nstart,1,vely1d_minus,velyminus);
        oned_unslice(nx,nstart,1,velz1d_minus,velzminus);
        oned_unslice(nx,nstart,1,velx1d_plus,velxplus);
        oned_unslice(nx,nstart,1,vely1d_plus,velyplus);
        oned_unslice(nx,nstart,1,velz1d_plus,velzplus);
        if (reconstruct_Wv)
           undo_Wv<0>(nx, velxminus, velyminus, velzminus,
                          velxplus, velyplus, velzplus,
                          g11, g12, g13, g22, g23, g33,
                          cctkGH, j, k); 
        
        oned_unslice(nx,nstart,1,eps1d_minus,epsminus);
        oned_unslice(nx,nstart,1,eps1d_plus,epsplus);
        
        if(do_ye) {
          oned_unslice(nx,nstart,1,ye1d_minus,Y_e_minus);
          oned_unslice(nx,nstart,1,ye1d_plus,Y_e_plus);
          if(do_temp) {
            oned_unslice(nx,nstart,1,temp1d_minus,tempminus);
            oned_unslice(nx,nstart,1,temp1d_plus,tempplus);
          }
        }
        
        if(do_mhd) {
          oned_unslice(nx,nstart,1,Bvcx1d_minus,Bvecxminus);
          oned_unslice(nx,nstart,1,Bvcy1d_minus,Bvecyminus);
          oned_unslice(nx,nstart,1,Bvcz1d_minus,Bveczminus);
          oned_unslice(nx,nstart,1,Bvcx1d_plus,Bvecxplus);
          oned_unslice(nx,nstart,1,Bvcy1d_plus,Bvecyplus);
          oned_unslice(nx,nstart,1,Bvcz1d_plus,Bveczplus);
          if(clean_divergence) {
            oned_unslice(nx,nstart,1,psidc1d_minus,psidcminus);
            oned_unslice(nx,nstart,1,psidc1d_plus,psidcplus);
          }
        }
        for (int i=0; i<nx; ++i) {
          SpaceMask_SetStateBits(space_mask, CCTK_GFINDEX3D(cctkGH, i,j,k), type_bitsx, not_trivialx);
        }
      }
    }  
    } // pragma omp parallel

  } else if (*flux_direction==2) {

    #pragma omp parallel
    {
    vector<CCTK_REAL> rho1d(nmax_xyz), eps1d(nmax_xyz), press1d(nmax_xyz);
    vector<CCTK_REAL> velx1d(nmax_xyz), vely1d(nmax_xyz), velz1d(nmax_xyz);
    vector<CCTK_REAL> temp1d(nmax_xyz), ye1d(nmax_xyz);
    vector<CCTK_REAL> Bvcx1d(nmax_xyz), Bvcy1d(nmax_xyz), Bvcz1d(nmax_xyz);
    vector<CCTK_REAL> psidc1d(nmax_xyz);

    vector<CCTK_REAL> rho1d_plus(nmax_xyz), eps1d_plus(nmax_xyz), 
                      press1d_plus(nmax_xyz);
    vector<CCTK_REAL> velx1d_plus(nmax_xyz), vely1d_plus(nmax_xyz),
                      velz1d_plus(nmax_xyz);
    vector<CCTK_REAL> temp1d_plus(nmax_xyz), ye1d_plus(nmax_xyz);
    vector<CCTK_REAL> Bvcx1d_plus(nmax_xyz), Bvcy1d_plus(nmax_xyz),
                      Bvcz1d_plus(nmax_xyz);
    vector<CCTK_REAL> psidc1d_plus(nmax_xyz);

    vector<CCTK_REAL> rho1d_minus(nmax_xyz), eps1d_minus(nmax_xyz), 
                      press1d_minus(nmax_xyz);
    vector<CCTK_REAL> velx1d_minus(nmax_xyz), vely1d_minus(nmax_xyz), 
                      velz1d_minus(nmax_xyz);
    vector<CCTK_REAL> temp1d_minus(nmax_xyz), ye1d_minus(nmax_xyz);
    vector<CCTK_REAL> Bvcx1d_minus(nmax_xyz), Bvcy1d_minus(nmax_xyz),
                      Bvcz1d_minus(nmax_xyz);
    vector<CCTK_REAL> psidc1d_minus(nmax_xyz);

    #pragma omp for
    for (int k = GRHydro_stencil-1; k < nz - GRHydro_stencil + 1 + transport_constraints; ++k) {
      for (int j = GRHydro_stencil-1; j < nx - GRHydro_stencil + 1 + transport_constraints; ++j) {
        
        const int nstart=j+k*nx*ny;
        
        oned_slice(ny,nstart,nx,rho,rho1d);
        
        if (!reconstruct_Wv) {
          oned_slice(ny,nstart,       nx,vup,velx1d);
          oned_slice(ny,nstart+nxyz,  nx,vup,vely1d);
          oned_slice(ny,nstart+2*nxyz,nx,vup,velz1d);
        } else {
          oned_slice(ny,nstart,       nx,&Wvup.front(),velx1d);
          oned_slice(ny,nstart+nxyz,  nx,&Wvup.front(),vely1d);
          oned_slice(ny,nstart+2*nxyz,nx,&Wvup.front(),velz1d);
        }
        oned_slice(ny,nstart,nx,eps,eps1d);
        oned_slice(ny,nstart,nx,press,press1d);
        
        if(do_ye) {
          oned_slice(ny,nstart,nx,Y_e,ye1d);
          if(do_temp){
            oned_slice(ny,nstart,nx,temperature,temp1d);
          } 
        }
        
        if(do_mhd) {
          oned_slice(ny,nstart,       nx,Bprim,Bvcx1d);
          oned_slice(ny,nstart+nxyz,  nx,Bprim,Bvcy1d);
          oned_slice(ny,nstart+2*nxyz,nx,Bprim,Bvcz1d);
          if(clean_divergence) {
            oned_slice(ny,nstart,nx,psidc,psidc1d);
          }
          if (do_Avec)
            CCTK_WARN(0, "Someone needs to figure out which variables to pass for Avec");
        }
        
        ppm1d_cxx_func
          (cctk_lsh[1],CCTK_DELTA_SPACE(1),
           &rho1d[0],&vely1d[0],&velz1d[0],&velx1d[0],&eps1d[0],&press1d[0],&temp1d[0],
           &ye1d[0],&Bvcy1d[0],&Bvcz1d[0],&Bvcx1d[0],&psidc1d[0],
           &rho1d_minus[0],&vely1d_minus[0],&velz1d_minus[0],&velx1d_minus[0],&eps1d_minus[0],&temp1d_minus[0],
           &ye1d_minus[0],&Bvcy1d_minus[0],&Bvcz1d_minus[0],&Bvcx1d_minus[0],&psidc1d_minus[0],
           &rho1d_plus[0],&vely1d_plus[0],&velz1d_plus[0],&velx1d_plus[0],&eps1d_plus[0],&temp1d_plus[0],
           &ye1d_plus[0],&Bvcy1d_plus[0],&Bvcz1d_plus[0],&Bvcx1d_plus[0],&psidc1d_plus[0]);
        
        
        oned_unslice(ny,nstart,nx,rho1d_minus,rhominus);
        oned_unslice(ny,nstart,nx,rho1d_plus,rhoplus);
        
        oned_unslice(ny,nstart,nx,velx1d_minus,velxminus);
        oned_unslice(ny,nstart,nx,vely1d_minus,velyminus);
        oned_unslice(ny,nstart,nx,velz1d_minus,velzminus);
        oned_unslice(ny,nstart,nx,velx1d_plus,velxplus);
        oned_unslice(ny,nstart,nx,vely1d_plus,velyplus);
        oned_unslice(ny,nstart,nx,velz1d_plus,velzplus);
        if (reconstruct_Wv)
           undo_Wv<1>(ny, velyminus, velzminus, velxminus,
                          velyplus, velzplus, velxplus,
                          g22, g23, g12, g33, g13, g11,
                          cctkGH, j, k);
        
        oned_unslice(ny,nstart,nx,eps1d_minus,epsminus);
        oned_unslice(ny,nstart,nx,eps1d_plus,epsplus);
        
        if(do_ye) {
          oned_unslice(ny,nstart,nx,ye1d_minus,Y_e_minus);
          oned_unslice(ny,nstart,nx,ye1d_plus,Y_e_plus);
          if(do_temp) {
            oned_unslice(ny,nstart,nx,temp1d_minus,tempminus);
            oned_unslice(ny,nstart,nx,temp1d_plus,tempplus);
          }
        }
        
        if(do_mhd) {
          oned_unslice(ny,nstart,nx,Bvcx1d_minus,Bvecxminus);
          oned_unslice(ny,nstart,nx,Bvcy1d_minus,Bvecyminus);
          oned_unslice(ny,nstart,nx,Bvcz1d_minus,Bveczminus);
          oned_unslice(ny,nstart,nx,Bvcx1d_plus,Bvecxplus);
          oned_unslice(ny,nstart,nx,Bvcy1d_plus,Bvecyplus);
          oned_unslice(ny,nstart,nx,Bvcz1d_plus,Bveczplus);
          if(clean_divergence) {
            oned_unslice(ny,nstart,nx,psidc1d_minus,psidcminus);
            oned_unslice(ny,nstart,nx,psidc1d_plus,psidcplus);
          }
        }
        for (int i=0; i<ny; ++i) {
          SpaceMask_SetStateBits(space_mask, CCTK_GFINDEX3D(cctkGH, j,i,k), type_bitsy, not_trivialy);
        }
        
      }
    }
    } // pragma omp parallel
    
  } else if (*flux_direction==3) {

    #pragma omp parallel
    {
    vector<CCTK_REAL> rho1d(nmax_xyz), eps1d(nmax_xyz), press1d(nmax_xyz);
    vector<CCTK_REAL> velx1d(nmax_xyz), vely1d(nmax_xyz), velz1d(nmax_xyz);
    vector<CCTK_REAL> temp1d(nmax_xyz), ye1d(nmax_xyz);
    vector<CCTK_REAL> Bvcx1d(nmax_xyz), Bvcy1d(nmax_xyz), Bvcz1d(nmax_xyz);
    vector<CCTK_REAL> psidc1d(nmax_xyz);

    vector<CCTK_REAL> rho1d_plus(nmax_xyz), eps1d_plus(nmax_xyz), 
                      press1d_plus(nmax_xyz);
    vector<CCTK_REAL> velx1d_plus(nmax_xyz), vely1d_plus(nmax_xyz),
                      velz1d_plus(nmax_xyz);
    vector<CCTK_REAL> temp1d_plus(nmax_xyz), ye1d_plus(nmax_xyz);
    vector<CCTK_REAL> Bvcx1d_plus(nmax_xyz), Bvcy1d_plus(nmax_xyz),
                      Bvcz1d_plus(nmax_xyz);
    vector<CCTK_REAL> psidc1d_plus(nmax_xyz);

    vector<CCTK_REAL> rho1d_minus(nmax_xyz), eps1d_minus(nmax_xyz), 
                      press1d_minus(nmax_xyz);
    vector<CCTK_REAL> velx1d_minus(nmax_xyz), vely1d_minus(nmax_xyz), 
                      velz1d_minus(nmax_xyz);
    vector<CCTK_REAL> temp1d_minus(nmax_xyz), ye1d_minus(nmax_xyz);
    vector<CCTK_REAL> Bvcx1d_minus(nmax_xyz), Bvcy1d_minus(nmax_xyz),
                      Bvcz1d_minus(nmax_xyz);
    vector<CCTK_REAL> psidc1d_minus(nmax_xyz);

    #pragma omp for
    for(int k = GRHydro_stencil-1; k < ny - GRHydro_stencil + 1 + transport_constraints; ++k) {
      for(int j = GRHydro_stencil-1; j < nx - GRHydro_stencil + 1 + transport_constraints; ++j) {
        
        const int nstart=j+k*nx;

        oned_slice(nz,nstart,nxy,rho,rho1d);
        
        if (!reconstruct_Wv) {
          oned_slice(nz,nstart,       nxy,vup,velx1d);
          oned_slice(nz,nstart+nxyz,  nxy,vup,vely1d);
          oned_slice(nz,nstart+2*nxyz,nxy,vup,velz1d);
        } else {
          oned_slice(nz,nstart,       nxy,&Wvup.front(),velx1d);
          oned_slice(nz,nstart+nxyz,  nxy,&Wvup.front(),vely1d);
          oned_slice(nz,nstart+2*nxyz,nxy,&Wvup.front(),velz1d);
        }
        
        oned_slice(nz,nstart,nxy,eps,eps1d);
        oned_slice(nz,nstart,nxy,press,press1d);
        
        if(do_ye) {
          oned_slice(nz,nstart,nxy,Y_e,ye1d);
          if(do_temp){
            oned_slice(nz,nstart,nxy,temperature,temp1d);
          } 
        }
        
        if(do_mhd) {
          oned_slice(nz,nstart,       nxy,Bprim,Bvcx1d);
          oned_slice(nz,nstart+nxyz,  nxy,Bprim,Bvcy1d);
          oned_slice(nz,nstart+2*nxyz,nxy,Bprim,Bvcz1d);
          if(clean_divergence) {
            oned_slice(nz,nstart,nxy,psidc,psidc1d);
          }
          if (do_Avec)
            CCTK_WARN(0, "Someone needs to figure out which variables to pass for Avec");
        }
        
        ppm1d_cxx_func
          (cctk_lsh[2],CCTK_DELTA_SPACE(2),
           &rho1d[0],&velz1d[0],&velx1d[0],&vely1d[0],&eps1d[0],&press1d[0],&temp1d[0],
           &ye1d[0],&Bvcz1d[0],&Bvcx1d[0],&Bvcy1d[0],&psidc1d[0],
           &rho1d_minus[0],&velz1d_minus[0],&velx1d_minus[0],&vely1d_minus[0],&eps1d_minus[0],&temp1d_minus[0],
           &ye1d_minus[0],&Bvcz1d_minus[0],&Bvcx1d_minus[0],&Bvcy1d_minus[0],&psidc1d_minus[0],
           &rho1d_plus[0],&velz1d_plus[0],&velx1d_plus[0],&vely1d_plus[0],&eps1d_plus[0],&temp1d_plus[0],
           &ye1d_plus[0],&Bvcz1d_plus[0],&Bvcx1d_plus[0],&Bvcy1d_plus[0],&psidc1d_plus[0]);
        
        
        oned_unslice(nz,nstart,nxy,rho1d_minus,rhominus);
        oned_unslice(nz,nstart,nxy,rho1d_plus,rhoplus);
        
        oned_unslice(nz,nstart,nxy,velx1d_minus,velxminus);
        oned_unslice(nz,nstart,nxy,vely1d_minus,velyminus);
        oned_unslice(nz,nstart,nxy,velz1d_minus,velzminus);
        oned_unslice(nz,nstart,nxy,velx1d_plus,velxplus);
        oned_unslice(nz,nstart,nxy,vely1d_plus,velyplus);
        oned_unslice(nz,nstart,nxy,velz1d_plus,velzplus);
        if (reconstruct_Wv)
           undo_Wv<2>(nz, velzminus, velxminus, velyminus,
                          velzplus, velxplus, velyplus,
                          g33, g13, g23, g11, g12, g22,
                          cctkGH, j, k);
        
        oned_unslice(nz,nstart,nxy,eps1d_minus,epsminus);
        oned_unslice(nz,nstart,nxy,eps1d_plus,epsplus);
        
        if(do_ye) {
          oned_unslice(nz,nstart,nxy,ye1d_minus,Y_e_minus);
          oned_unslice(nz,nstart,nxy,ye1d_plus,Y_e_plus);
          if(do_temp) {
            oned_unslice(nz,nstart,nxy,temp1d_minus,tempminus);
            oned_unslice(nz,nstart,nxy,temp1d_plus,tempplus);
          }
        }
        
        if(do_mhd) {
          oned_unslice(nz,nstart,nxy,Bvcx1d_minus,Bvecxminus);
          oned_unslice(nz,nstart,nxy,Bvcy1d_minus,Bvecyminus);
          oned_unslice(nz,nstart,nxy,Bvcz1d_minus,Bveczminus);
          oned_unslice(nz,nstart,nxy,Bvcx1d_plus,Bvecxplus);
          oned_unslice(nz,nstart,nxy,Bvcy1d_plus,Bvecyplus);
          oned_unslice(nz,nstart,nxy,Bvcz1d_plus,Bveczplus);
          if(clean_divergence) {
            oned_unslice(nz,nstart,nxy,psidc1d_minus,psidcminus);
            oned_unslice(nz,nstart,nxy,psidc1d_plus,psidcplus);
          }
        }
        for (int i=0; i<nz; ++i) {
          SpaceMask_SetStateBits(space_mask, CCTK_GFINDEX3D(cctkGH, j,k,i), type_bitsz, not_trivialz);
        }
        
      }
    }
    } // pragma omp parallel

  }
}
