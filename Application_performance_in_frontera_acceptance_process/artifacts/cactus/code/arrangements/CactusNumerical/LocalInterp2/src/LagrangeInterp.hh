#ifndef LOCALINTERP_LAGRANGEINTERP_HH
#define LOCALINTERP_LAGRANGEINTERP_HH

#include <cassert>
#include <cmath>
#include <cstring>

#include "cctk.h"

// If this is uncommented always use symmetric operators
#ifndef LOCALINTERP_SYMMETRIC
#define LOCALINTERP_SYMMETRIC
#endif

template<int order>
class LagrangeInterp1D {
  public:
    LagrangeInterp1D(
        //! [in] Grid origin
        CCTK_REAL const origin,
        //! [in] Grid spacing
        CCTK_REAL const delta,
        //! [in] Number of grid points
        CCTK_INT siz,
        //! [in] Interpolation point
        CCTK_REAL const coord):
      m_origin(origin),
      m_delta(delta),
      m_delta_inv(1.0/m_delta),
      m_siz(siz),
      m_coord(coord),
      m_mid_flag(false),
      m_npoint(order + 1),
      m_out_of_bounds(false),
      point(m_point),
      npoint(m_npoint),
      out_of_bounds(m_out_of_bounds) {
      // Check if we are in the middle between two grid points
#ifdef LOCALINTERP_SYMMETRIC
      if(0 == order % 2) {
        int idx = std::floor((m_coord - m_origin)*m_delta_inv);
        if(m_coord - (idx*m_delta + m_origin) == 0.5 * m_delta) {
          m_mid_flag = true;
          ++m_npoint;
        }
        else {
          m_mid_flag = false;
        }
      }
#endif
      // First point (from the left) of the interpolation stencil
      m_point = std::floor((m_coord - m_origin)*m_delta_inv
          - 0.5*(order - 1));
#ifdef LOCALINTERP_SYMMETRIC
      if(m_mid_flag) {
          m_point -= 1;
      }
#endif

      // Shift the interpolation stencil if out of the grid
      CCTK_INT shift = m_point;
      if(shift < 0) {
        m_point -= shift;
        m_out_of_bounds = true;
      }
      shift = m_point + order - (m_siz - 1) + m_mid_flag;
      if(shift > 0) {
        m_point -= shift;
        m_out_of_bounds = true;
      }

      CCTK_REAL xp[order+2];
      for(int i = 0; i <= order; ++i) {
        xp[i] = (m_point + i) * m_delta + m_origin;
      }
#ifdef LOCALINTERP_SYMMETRIC
      if(0 == order % 2 && m_mid_flag) {
        xp[order+1] = (m_point + order + 1) * m_delta + m_origin;
      }
#endif
      m_calc_coeff_lr(xp, &m_coeff_lr[0]);
#ifdef LOCALINTERP_SYMMETRIC
      if(0 == order % 2 && m_mid_flag) {
        m_calc_coeff_rl(&xp[1], &m_coeff_rl[0]);
      }
      else {
        std::memcpy(m_coeff_rl, m_coeff_lr, sizeof(m_coeff_lr));
      }
#endif
    }

    //! Evaluates the interpolator
    template<typename T>
    T eval(
        //! [in] must be offset so that vals[0] = vals["point"]
        T const * const vals,
        //! [in] stride used to access vals
        CCTK_INT const stride
        ) const {
      T out_lr = 0;
      for(int i = 0; i <= order; ++i) {
        out_lr += static_cast<T>(m_coeff_lr[i]) * vals[i*stride];
      }
#ifdef LOCALINTERP_SYMMETRIC
      T out_rl = 0;
      int shift = static_cast<int>(0 == order % 2 && m_mid_flag);
      for(int i = order; i >= 0; --i) {
        out_rl += static_cast<T>(m_coeff_rl[i]) * vals[(i + shift)*stride];
      }
      return T(0.5)*(out_lr + out_rl);
#else
      return out_lr;
#endif
    }
  private:
    // Compute the Lagrange interpolation coefficients on a given stencil
    void m_calc_coeff_lr(
        CCTK_REAL const * const xp,
        CCTK_REAL * const coeff
        ) const {
#define TYPECASE(I0, TEST, I1, OP)                                            \
      for(int j = 0; j <= order; ++j) {                                       \
        CCTK_REAL num = 1.0;                                                  \
        CCTK_REAL den = 1.0;                                                  \
        for(int i = I0; i TEST I1; OP i) {                                    \
          if(i == j) {                                                        \
            continue;                                                         \
          }                                                                   \
          num = num * (m_coord - xp[i]);                                      \
          den = den * (xp[j] - xp[i]);                                        \
        }                                                                     \
        coeff[j] = num/den;                                                   \
      }
      TYPECASE(0, <=, order, ++)
    }
    void m_calc_coeff_rl(
        CCTK_REAL const * const xp,
        CCTK_REAL * const coeff
        ) const {
      TYPECASE(order, >=, 0, --)
#undef TYPECASE
    }
  private:
    CCTK_REAL m_origin;
    CCTK_REAL m_delta;
    CCTK_REAL m_delta_inv;
    CCTK_INT m_siz;

    CCTK_REAL m_coord;

    // If true we have an asymmetric stencil, but the interpolation point is
    // exactly in the middle between two grid points. In this case we need to
    // average the results obtained from the interpolation on two separate
    // stencils.
    bool m_mid_flag;
    // First point (going from left to right) of the stencil
    CCTK_INT m_point;
    // Number of points needed for the interpolation
    int m_npoint;
    // The stencil was shifted to avoid going out of bounds
    bool m_out_of_bounds;

    // Interpolation coefficients for interpolation from left to right
    CCTK_REAL m_coeff_lr[order+1];
#ifdef LOCALINTERP_SYMMETRIC
    // Interpolation coefficients for interpolation from right to left
    CCTK_REAL m_coeff_rl[order+1];
#endif
  public:
    //! Index of the first point of the interpolation stencil
    CCTK_INT const & point;
    //! Number of points needed for the interpolation
    int const & npoint;
    //! The stencil
    bool const & out_of_bounds;
};

template<int ndim, int D>
class NextStencil {
  public:
    enum { value = D - 1 };
};

template<int ndim>
class NextStencil<ndim, 0> {
  public:
    enum { value = 0 };
};

template<int order, int ndim>
class LagrangeInterpND {
  public:
    LagrangeInterpND(
        //! [in] Grid origin
        CCTK_REAL const origin[ndim],
        //! [in] Grid spacing
        CCTK_REAL const delta[ndim],
        //! [in] Number of grid points
        CCTK_INT const siz[ndim],
        //! [in] Interpolation point
        CCTK_REAL const coord[ndim]):
        m_out_of_bounds(false),
        out_of_bounds(m_out_of_bounds) {
      for(int d = 0; d < ndim; ++d) {
        m_origin[d] = origin[d];
        m_delta[d] = delta[d];
        m_siz[d] = siz[d];
        m_coord[d] = coord[d];
        mp_interp[d] = new (&m_interp_scratch[d][0])
          LagrangeInterp1D<order>(m_origin[d], m_delta[d],
              m_siz[d], m_coord[d]);
        m_out_of_bounds = m_out_of_bounds || mp_interp[d]->out_of_bounds;
      }
    }

    template<typename T>
    T eval(
        //! [in] Grid function to interpolate
        T const * const gf
        ) const {
      T vals[ndim][order+2];
      int pos[ndim];
      m_fill_stencil<T, ndim-1>(gf, pos, vals);
      return mp_interp[ndim-1]->eval(vals[ndim-1], 1);
    }
  private:
    // Recursively fill the stencil used for the interpolation
    template<typename T, int D>
    void m_fill_stencil(
        T const * const gf,
        int pos[ndim],
        T vals[ndim][order+2]
        ) const {
      assert(D >= 0 && D < ndim);
      if(D == 0) {
        CCTK_INT gidx = mp_interp[0]->point;
        CCTK_INT stride = 1;
        for(int d = 1; d < ndim; ++d) {
          stride *= m_siz[d-1];
          gidx += stride * (mp_interp[d]->point + pos[d]);
        }
        std::memcpy(&vals[0][0], &gf[gidx], mp_interp[0]->npoint*sizeof(T));
      }
      else {
        for(pos[D] = 0; pos[D] < mp_interp[D]->npoint; ++pos[D]) {
          m_fill_stencil<T, NextStencil<ndim, D>::value>(gf, pos, vals);
          vals[D][pos[D]] = mp_interp[D-1]->eval(vals[D-1], 1);
        }
      }
    }
  private:
    CCTK_REAL m_origin[ndim];
    CCTK_REAL m_delta[ndim];
    CCTK_INT m_siz[ndim];

    CCTK_REAL m_coord[ndim];
    bool m_out_of_bounds;

    // 1D interpolators (will be placed on the stack)
    LagrangeInterp1D<order> * mp_interp[ndim];
    // Scratch space used for placement new
    char m_interp_scratch[ndim][sizeof(LagrangeInterp1D<order>)];
  public:
    bool const & out_of_bounds;
};

#endif
