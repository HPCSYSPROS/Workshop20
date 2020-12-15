#include <cctk.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <sstream>

#ifdef CCTK_MPI
#include <mpi.h>
#else
#include "nompi.h"
#endif

#include <cacheinfo.hh>
#include <carpet.hh>
#include <vect.hh>
#include <Timer.hh>

#include "fasterp.hh"

#ifdef CARPETINTERP2_CHECK
#define CI2C(x, y) x, y
#else
#define CI2C(x, y)
#endif

namespace CarpetInterp2 {

int const ipoison = -1234567890;
CCTK_REAL const poison = -1.0e+12;

int get_poison(int const &) CCTK_ATTRIBUTE_CONST;
int get_poison(int const &) { return ipoison; }
CCTK_REAL get_poison(CCTK_REAL const &) CCTK_ATTRIBUTE_CONST;
CCTK_REAL get_poison(CCTK_REAL const &) { return poison; }

template <typename T> void fill(vector<T> &v, T const &val) {
// fill (v.begin(), v.end(), val);
#pragma omp parallel for
  for (int n = 0; n < int(v.size()); ++n) {
    v.AT(n) = val;
  }
}

template <typename T> void fill_with_poison(vector<T> &v) {
#ifndef NDEBUG
  T dummy;
  fill(v, get_poison(dummy));
// fill (v.begin(), v.end(), get_poison(dummy));
#endif
}

template <typename T> void check_for_poison(vector<T> const &v) {
#pragma omp parallel for
  for (int n = 0; n < int(v.size()); ++n) {
    T dummy;
    assert(v.AT(n) != get_poison(dummy));
  }
}

// Create an MPI datatype for location information
MPI_Datatype fasterp_iloc_t::mpi_datatype() {
  static bool initialised = false;
  static MPI_Datatype newtype;
  if (not initialised) {
    static fasterp_iloc_t s;
#define ENTRY(type, name)                                                      \
  {                                                                            \
    sizeof s.name / sizeof(type),         /* count elements */                 \
        (char *) & s.name - (char *) & s, /* offsetof doesn't work (why?) */   \
        dist::mpi_datatype<type>(),       /* find MPI datatype */              \
        STRINGIFY(name),                  /* field name */                     \
        STRINGIFY(type),                  /* type name */                      \
  }
    dist::mpi_struct_descr_t const descr[] = {
        ENTRY(int, mrc),
#ifdef CARPETINTERP2_CHECK
        ENTRY(int, pn),
        ENTRY(int, ipos),
        ENTRY(int, ind),
#endif
        ENTRY(int, ind3d),
        ENTRY(CCTK_REAL, offset),
        {1, sizeof(s), MPI_UB, "MPI_UB", "MPI_UB"}};
#undef ENTRY
    newtype =
        dist::create_mpi_datatype(sizeof descr / sizeof descr[0], descr,
                                  "fasterp_iloc_t::mpi_datatype", sizeof s);
    initialised = true;
  }
  return newtype;
}

void fasterp_iloc_t::output(ostream &os) const {
  os << "fasterp_iloc_t{"
     << "mrc=" << mrc << ","
#ifdef CARPETINTERP2_CHECK
     << "pn=" << pn << ","
     << "ipos=" << ipos << ","
     << "ind=" << ind << ","
#endif
     << "ind3d=" << ind3d << ","
     << "offset=" << offset << "}";
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

int mrc_t::maps = -1;
int mrc_t::reflevels = -1;
int mrc_t::components = -1;

void mrc_t::determine_mrc_info() {
  maps = Carpet::maps;
  reflevels = Carpet::reflevels;
  components = 0;
  for (int m = 0; m < maps; ++m) {
    for (int rl = 0; rl < reflevels; ++rl) {
      int const ncomps = Carpet::vhh.AT(m)->components(rl);
      components = max(components, ncomps);
    }
  }
}

mrc_t::mrc_t(int ind) {
  assert(ind >= 0);
  int const ind0 = ind;
  c = ind % components;
  ind /= components;
  rl = ind % reflevels;
  ind /= reflevels;
  m = ind % maps;
  ind /= maps;
  assert(ind == 0);
  assert(get_ind() == ind0);
}

void mrc_t::output(ostream &os) const {
  os << "mrc{m=" << m << ",rl=" << rl << ",c=" << c << "}";
}

void pn_t::output(ostream &os) const {
  os << "pn{p=" << p << ",n=" << n << "}";
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

// Calculate the coefficients for one interpolation stencil
// TODO: Could templatify this function on the order to improve
// efficiency
int fasterp_src_loc_t::calc_stencil(fasterp_iloc_t const &iloc,
                                    ivect const &ash,
#ifdef CARPETINTERP2_CHECK
                                    ivect const &lsh,
#endif
                                    int const order) {
  assert(order <= max_order);
  CCTK_REAL const eps = 1.0e-12;

#ifdef CARPETINTERP2_CHECK
  mrc = iloc.mrc;
  pn = iloc.pn;
  ipos = iloc.ipos;

  // Save ash
  saved_ash = ash;
#endif

#ifndef NDEBUG
  // Poison all coefficients
  for (int d = 0; d < dim; ++d) {
    for (int n = 0; n <= max_order; ++n) {
      coeffs[d][n] = poison;
    }
  }
#endif

  if (order == 0) {
    // Special case
    for (int d = 0; d < dim; ++d) {
      exact[d] = true;
    }

#ifdef CARPETINTERP2_CHECK
    ind = iloc.ind;
#endif
    ind3d = iloc.ind3d;
    return 0;
  }

  // Choose stencil anchor (first stencil point)
  ivect iorigin;
  if (order % 2 == 0) {
    iorigin = -ivect(order / 2);
  } else {
    // Potentially shift the stencil anchor for odd interpolation
    // orders (i.e., for even numbers of stencil points)
    ivect const ioffset(iloc.offset < CCTK_REAL(0.0));
    iorigin = -ivect((order - 1) / 2) - ioffset;
  }
  rvect const offset = iloc.offset - rvect(iorigin);
  // Ensure that the stencil is centred
  assert(all(offset >= CCTK_REAL(0.5) * (order - 1) and
             offset < CCTK_REAL(0.5) * (order + 1)));

  for (int d = 0; d < dim; ++d) {
    // C_n = PRODUCT_m,m!=n [(x - x_m) / (x_n - x_m)]
    CCTK_REAL const x = offset[d];
    // round is not available with PGI compilers
    // CCTK_REAL const rx = round(x);
    CCTK_REAL const rx = floor(x + 0.5);
    if (abs(x - rx) < eps * (1.0 + abs(x))) {
      // The interpolation point coincides with a grid point; no
      // interpolation is necessary (this is a special case)
      iorigin[d] += int(rx);
      exact[d] = true;
    } else {
      for (int n = 0; n <= order; ++n) {
        CCTK_REAL const xn = n;
        CCTK_REAL c = 1.0;
        for (int m = 0; m <= order; ++m) {
          if (m != n) {
            CCTK_REAL const xm = m;
            c *= (x - xm) / (xn - xm);
          }
        }
        coeffs[d][n] = c;
      }
      exact[d] = false;
      // The sum of the coefficients must be one
      CCTK_REAL s = 0.0;
      for (int n = 0; n <= order; ++n) {
        s += coeffs[d][n];
      }
      assert(abs(s - 1.0) <= eps);
    }
  }

// Set 3D location of stencil anchor
#ifdef CARPETINTERP2_CHECK
  ind = iloc.ind + iorigin;
  if (not(all(ind >= 0 and ind + either(exact, 0, order) < lsh))) {
    stringstream buf;
    buf << "*this=" << *this << " iloc=" << iloc << " "
        << "lsh=" << lsh << " order=" << order;
    CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Array access out of bounds: %s. Most likely, the interpolation "
               "point is too close to an outer or to a symmetry/inter-patch "
               "boundary. More information follows...",
               buf.str().c_str());
    return -1;
  }
#endif
  ind3d = iloc.ind3d + index(ash, iorigin);
#ifdef CARPETINTERP2_CHECK
  assert(index(ash, ind) == ind3d);
#endif
  return 0;
}

// Interpolate a set of variables at the same location, with a given
// set of interpolation coefficients.  The template mechanism allows
// the order in each direction to be different; in particular, the
// order can be 0 if the interpolation location coincides with a
// source grid point.  For interpolation to order O, this function
// should be instantiated eight times, with On=O and On=0, for
// n=[0,1,2].  We hope that the compiler optimises the for loops
// away if On=0.
template <int O0, int O1, int O2>
void fasterp_src_loc_t::interpolate(ivect const &ash,
#ifdef CARPETINTERP2_CHECK
                                    ivect const &lsh,
#endif
                                    vector<CCTK_REAL const *> const &varptrs,
                                    CCTK_REAL *restrict const vals) const {
#ifdef CARPETINTERP2_CHECK
  assert(all(ash == saved_ash));
#endif
  assert(O0 == 0 or not exact[0]);
  assert(O1 == 0 or not exact[1]);
  assert(O2 == 0 or not exact[2]);

  size_t const di = 1;
  size_t const dj = di * ash[0];
  size_t const dk = dj * ash[1];

#ifdef CARPETINTERP2_CHECK
  assert(all(ind >= 0 and ind + ivect(O0, O1, O2) < lsh));
  assert(ind3d == index(ash, ind));
  assert(int(di * ind[0] + dj * ind[1] + dk * ind[2]) == index(ash, ind));
#endif

  for (size_t v = 0; v < varptrs.size(); ++v) {
    CCTK_REAL tmp = 0.0;
    const CCTK_REAL *restrict const varptr = &varptrs.AT(v)[ind3d];

    for (size_t k = 0; k <= O2; ++k) {
      CCTK_REAL const coeff_k = O2 == 0 ? 1.0 : coeffs[2][k];
      for (size_t j = 0; j <= O1; ++j) {
        CCTK_REAL const coeff_jk = coeff_k * (O1 == 0 ? 1.0 : coeffs[1][j]);
#ifndef __BIGGEST_ALIGNMENT__    // this is a GCC extension
#define __BIGGEST_ALIGNMENT__ 64 // a guess
#endif
        CCTK_REAL buffer[O0 + 1]
            __attribute__((__aligned__(__BIGGEST_ALIGNMENT__)));
#ifdef __INTEL_COMPILER
#pragma vector always
#else
#pragma omp simd
#endif
        for (size_t i = 0; i <= O0; ++i) {
          buffer[i] = varptr[i * di + j * dj + k * dk];
        }
#ifdef __INTEL_COMPILER
#pragma vector always
#else
#pragma omp simd reduction(+ : tmp)
#endif
        for (size_t i = 0; i <= O0; ++i) {
          CCTK_REAL const coeff_ijk = coeff_jk * (O0 == 0 ? 1.0 : coeffs[0][i]);
          tmp += coeff_ijk * buffer[i];
        }
      }
    }
    vals[v] = tmp;
  } // for v

#ifdef CARPETINTERP2_CHECK
#if 0
    for (size_t v=0; v<varptrs.size(); ++v) {
      vals[v] = ind[0] + 1000 * (ind[1] + 1000 * ind[2]);
    } // for v
#endif
#if 0
    for (size_t v=0; v<varptrs.size(); ++v) {
      vals[v] = pn.p * 1000000 + pn.n;
    } // for v
#endif
#endif
}

// Interpolate a set of variables at the same location, with a given
// set of interpolation coefficients.  This calls the specialised
// interpolation function, depending on whether interpolation is
// exact in each of the directions.
template <int O>
void fasterp_src_loc_t::interpolate(ivect const &ash,
#ifdef CARPETINTERP2_CHECK
                                    ivect const &lsh,
#endif
                                    vector<CCTK_REAL const *> const &varptrs,
                                    CCTK_REAL *restrict const vals) const {
  int const Z = 0;
  if (exact[2]) {
    if (exact[1]) {
      if (exact[0]) {
        interpolate<Z, Z, Z>(ash, CI2C(lsh, ) varptrs, vals);
      } else {
        interpolate<O, Z, Z>(ash, CI2C(lsh, ) varptrs, vals);
      }
    } else {
      if (exact[0]) {
        interpolate<Z, O, Z>(ash, CI2C(lsh, ) varptrs, vals);
      } else {
        interpolate<O, O, Z>(ash, CI2C(lsh, ) varptrs, vals);
      }
    }
  } else {
    if (exact[1]) {
      if (exact[0]) {
        interpolate<Z, Z, O>(ash, CI2C(lsh, ) varptrs, vals);
      } else {
        interpolate<O, Z, O>(ash, CI2C(lsh, ) varptrs, vals);
      }
    } else {
      if (exact[0]) {
        interpolate<Z, O, O>(ash, CI2C(lsh, ) varptrs, vals);
      } else {
        interpolate<O, O, O>(ash, CI2C(lsh, ) varptrs, vals);
      }
    }
  }
}

// Interpolate a set of variables at the same location, with a given
// set of interpolation coefficients.  This calls a specialised
// interpolation function, depending on whether interpolation is
// exact in each of the directions.
void fasterp_src_loc_t::interpolate(ivect const &ash,
#ifdef CARPETINTERP2_CHECK
                                    ivect const &lsh,
#endif
                                    int const order,
                                    vector<CCTK_REAL const *> const &varptrs,
                                    CCTK_REAL *restrict const vals) const {
  switch (order) {
  case 0:
    interpolate<0>(ash, CI2C(lsh, ) varptrs, vals);
    break;
  case 1:
    interpolate<1>(ash, CI2C(lsh, ) varptrs, vals);
    break;
  case 2:
    interpolate<2>(ash, CI2C(lsh, ) varptrs, vals);
    break;
  case 3:
    interpolate<3>(ash, CI2C(lsh, ) varptrs, vals);
    break;
  case 4:
    interpolate<4>(ash, CI2C(lsh, ) varptrs, vals);
    break;
  case 5:
    interpolate<5>(ash, CI2C(lsh, ) varptrs, vals);
    break;
  case 6:
    interpolate<6>(ash, CI2C(lsh, ) varptrs, vals);
    break;
  case 7:
    interpolate<7>(ash, CI2C(lsh, ) varptrs, vals);
    break;
  case 8:
    interpolate<8>(ash, CI2C(lsh, ) varptrs, vals);
    break;
  case 9:
    interpolate<9>(ash, CI2C(lsh, ) varptrs, vals);
    break;
  case 10:
    interpolate<10>(ash, CI2C(lsh, ) varptrs, vals);
    break;
  case 11:
    interpolate<11>(ash, CI2C(lsh, ) varptrs, vals);
    break;
  default:
    // Add higher orders here as desired
    CCTK_ERROR("Interpolation orders larger than 11 are not yet implemented");
    assert(0);
  }
}

void fasterp_src_loc_t::output(ostream &os) const {
  os << "fasterp_src_loc_t{";
  os << "coeffs=[";
  for (int d = 0; d < dim; ++d) {
    if (d > 0)
      os << ",";
    os << "[";
    for (int n = 0; n <= max_order; ++n) {
      if (n > 0)
        os << ",";
      os << coeffs[d][n];
    }
    os << "]";
  }
  os << "],";
  os << "exact=" << exact << ",";
#ifdef CARPETINTERP2_CHECK
  os << "pn=" << pn << ",";
  os << "mrc=" << mrc << ",";
  os << "ipos=" << ipos << ",";
  os << "ind=" << ind << ",";
#endif
  os << "ind3d=" << ind3d;
#ifdef CARPETINTERP2_CHECK
  os << ","
     << "saved_ash=" << saved_ash;
#endif
  os << "}";
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

// Calculate the coefficients for one interpolation stencil
// TODO: Could templatify this function on the order to improve
// efficiency
int fasterp_eno2_src_loc_t::calc_stencil(fasterp_iloc_t const &iloc,
                                         ivect const &ash,
#ifdef CARPETINTERP2_CHECK
                                         ivect const &lsh,
#endif
                                         int /*order*/) {
  // assert (order <= max_order);
  CCTK_REAL const eps = 1.0e-12;

#ifdef CARPETINTERP2_CHECK
  mrc = iloc.mrc;
  pn = iloc.pn;
  ipos = iloc.ipos;

  // Save ash
  saved_ash = ash;
#endif

  // first we setup 1st order stencil!
  {
    const int order = 1;

#ifndef NDEBUG
    // Poison all coefficients
    for (int d = 0; d < dim; ++d) {
      for (int n = 0; n <= order; ++n) {
        coeffs[d][n] = poison;
      }
    }
#endif

    // Choose stencil anchor (first stencil point)
    ivect iorigin;

    // Potentially shift the stencil anchor for odd interpolation
    // orders (i.e., for even numbers of stencil points)
    ivect const ioffset(iloc.offset < CCTK_REAL(0.0));
    iorigin = -ivect((order - 1) / 2) - ioffset;
    rvect const offset = iloc.offset - rvect(iorigin);
    // Ensure that the stencil is centred
    assert(all(offset >= CCTK_REAL(0.5) * (order - 1) and
               offset < CCTK_REAL(0.5) * (order + 1)));

    for (int d = 0; d < dim; ++d) {
      // C_n = PRODUCT_m,m!=n [(x - x_m) / (x_n - x_m)]
      CCTK_REAL const x = offset[d];
      // round is not available with PGI compilers
      // CCTK_REAL const rx = round(x);
      CCTK_REAL const rx = floor(x + 0.5);
      if (abs(x - rx) < eps * (1.0 + abs(x))) {
        // The interpolation point coincides with a grid point; no
        // interpolation is necessary (this is a special case)
        iorigin[d] += int(rx);
        exact[d] = true;
      } else {
        for (int n = 0; n <= order; ++n) {
          CCTK_REAL const xn = n;
          CCTK_REAL c = 1.0;
          for (int m = 0; m <= order; ++m) {
            if (m != n) {
              CCTK_REAL const xm = m;
              c *= (x - xm) / (xn - xm);
            }
          }
          coeffs[d][n] = c;
        }
        exact[d] = false;
        // The sum of the coefficients must be one
        CCTK_REAL s = 0.0;
        for (int n = 0; n <= order; ++n) {
          s += coeffs[d][n];
        }
        assert(abs(s - 1.0) <= eps);
      }
    }
  }

  // second, we setup left 2nd order stencil!
  {
    int order = 2;

#ifndef NDEBUG
    // Poison all coefficients
    for (int d = 0; d < dim; ++d) {
      for (int n = 0; n <= order; ++n) {
        coeffsLeft[d][n] = poison;
      }
    }
#endif

    // Choose stencil anchor (first stencil point)
    ivect iorigin;
    iorigin = -ivect(order / 2) - ivect((iloc.offset < CCTK_REAL(0.0)));
    rvect const offset = iloc.offset - rvect(iorigin);
    // cout << "Left " << iorigin << " " << iloc.offset << " " << offset <<
    // endl;
    // Ensure that interpolation point is between second and third point
    assert(all(offset >= CCTK_REAL(1.0) and offset <= CCTK_REAL(2.0)));

    for (int d = 0; d < dim; ++d) {
      // C_n = PRODUCT_m,m!=n [(x - x_m) / (x_n - x_m)]
      CCTK_REAL const x = offset[d];
      // round is not available with PGI compilers
      // CCTK_REAL const rx = round(x);
      CCTK_REAL const rx = floor(x + 0.5);
      if (abs(x - rx) < eps * (1.0 + abs(x))) {
        // The interpolation point coincides with a grid point; no
        // interpolation is necessary (this is a special case)
        iorigin[d] += int(rx);
        exact[d] = true;
      } else {
        for (int n = 0; n <= order; ++n) {
          CCTK_REAL const xn = n;
          CCTK_REAL c = 1.0;
          for (int m = 0; m <= order; ++m) {
            if (m != n) {
              CCTK_REAL const xm = m;
              c *= (x - xm) / (xn - xm);
            }
          }
          coeffsLeft[d][n] = c;
        }
        exact[d] = false;
        // The sum of the coefficients must be one
        CCTK_REAL s = 0.0;
        for (int n = 0; n <= order; ++n) {
          s += coeffsLeft[d][n];
        }
        assert(abs(s - 1.0) <= eps);
      }
    }

// Set 3D location of stencil anchor.
// Since left stencil extends farthest to the left, we use
// this to compute stencil anchor
#ifdef CARPETINTERP2_CHECK
    ind = iloc.ind + iorigin;
    if (not(all(ind >= 0 and ind + either(exact, 0, order) < lsh))) {
      stringstream buf;
      buf << "*this=" << *this << " iloc=" << iloc << " "
          << "lsh=" << lsh << " order=" << order;
      CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Array access out of bounds: %s. Most likely, the "
                 "interpolation point is too close to an outer or to a "
                 "symmetry/inter-patch boundary. More information follows...",
                 buf.str().c_str());
      return -1;
    }
#endif
    ind3d = iloc.ind3d + index(ash, iorigin);
#ifdef CARPETINTERP2_CHECK
    assert(index(ash, ind) == ind3d);
#endif
  }

  // third, we setup right 2nd order stencil!
  {
    int order = 2;

#ifndef NDEBUG
    // Poison all coefficients
    for (int d = 0; d < dim; ++d) {
      for (int n = 0; n <= order; ++n) {
        coeffsLeft[d][n] = poison;
      }
    }
#endif

    // Choose stencil anchor (first stencil point)
    ivect iorigin;
    iorigin =
        -ivect(order / 2) + ivect(1.0) - ivect((iloc.offset < CCTK_REAL(0.0)));
    rvect const offset = iloc.offset - rvect(iorigin);
    // cout << "Right " << iorigin << " " << iloc.offset << " " << offset <<
    // endl;
    // Ensure that the interpolation point is between first and second point
    assert(all(offset >= CCTK_REAL(0.0) and offset <= CCTK_REAL(1.0)));

    for (int d = 0; d < dim; ++d) {
      // C_n = PRODUCT_m,m!=n [(x - x_m) / (x_n - x_m)]
      CCTK_REAL const x = offset[d];
      // round is not available with PGI compilers
      // CCTK_REAL const rx = round(x);
      CCTK_REAL const rx = floor(x + 0.5);
      if (abs(x - rx) < eps * (1.0 + abs(x))) {
        // The interpolation point coincides with a grid point; no
        // interpolation is necessary (this is a special case)
        iorigin[d] += int(rx);
        exact[d] = true;
      } else {
        for (int n = 0; n <= order; ++n) {
          CCTK_REAL const xn = n;
          CCTK_REAL c = 1.0;
          for (int m = 0; m <= order; ++m) {
            if (m != n) {
              CCTK_REAL const xm = m;
              c *= (x - xm) / (xn - xm);
            }
          }
          coeffsRight[d][n] = c;
        }
        exact[d] = false;
        // The sum of the coefficients must be one
        CCTK_REAL s = 0.0;
        for (int n = 0; n <= order; ++n) {
          s += coeffsRight[d][n];
        }
        assert(abs(s - 1.0) <= eps);
      }
    }
  }

  return 0;
}

// Interpolate a set of variables at the same location, with a given
// set of interpolation coefficients.  The template mechanism allows
// the order in each direction to be different; in particular, the
// order can be 0 if the interpolation location coincides with a
// source grid point.  For interpolation to order O, this function
// should be instantiated eight times, with On=O and On=0, for
// n=[0,1,2].  We hope that the compiler optimises the for loops
// away if On=0.
template <int O0, int O1, int O2>
void fasterp_eno2_src_loc_t::interpolate(
    ivect const &ash,
#ifdef CARPETINTERP2_CHECK
    ivect const &lsh,
#endif
    vector<CCTK_REAL const *> const &varptrs,
    CCTK_REAL *restrict const vals) const {
#ifdef CARPETINTERP2_CHECK
  assert(all(ash == saved_ash));
#endif
  assert(O0 == 0 or not exact[0]);
  assert(O1 == 0 or not exact[1]);
  assert(O2 == 0 or not exact[2]);

  size_t const di = 1;
  size_t const dj = di * ash[0];
  size_t const dk = dj * ash[1];

#ifdef CARPETINTERP2_CHECK
  assert(all(ind >= 0 and ind + ivect(O0, O1, O2) < lsh));
  assert(ind3d == index(ash, ind));
  assert(int(di * ind[0] + dj * ind[1] + dk * ind[2]) == index(ash, ind));
#endif

  for (size_t v = 0; v < varptrs.size(); ++v) {
    vals[v] = 0.0;
  }

  // Must construct interpolation stencils for each variable
  // independently!
  for (size_t v = 0; v < varptrs.size(); ++v) {

    vector<CCTK_REAL> qz(4);
    for (size_t k = 0; k <= 3; ++k) {

      vector<CCTK_REAL> qy(4);
      for (size_t j = 0; j <= 3; ++j) {

        vector<CCTK_REAL> qx(4);
        for (size_t i = 0; i <= 3; ++i) {

          qx[i] = varptrs.AT(v)[ind3d + i * di + j * dj + k * dk];
        }

        switch (O0) {
        case 0:
          qy[j] = qx[0];
          break;

        default: {
          const CCTK_REAL diffleft = qx[0] + qx[2] - 2.0 * qx[1];
          const CCTK_REAL diffright = qx[1] + qx[3] - 2.0 * qx[2];

          CCTK_REAL res;
          if (fabs(diffleft) < fabs(diffright))
            res = coeffsLeft[0][0] * qx[0] + coeffsLeft[0][1] * qx[1] +
                  coeffsLeft[0][2] * qx[2];
          else
            res = coeffsRight[0][0] * qx[1] + coeffsRight[0][1] * qx[2] +
                  coeffsRight[0][2] * qx[3];

          if (diffleft * diffright <= 0.0 ||
              (res - qx[1]) * (qx[2] - res) < 0.0) {
            res = coeffs[0][0] * qx[1] + coeffs[0][1] * qx[2];
          }

          qy[j] = res;
          break;
        }
        }
      }

      switch (O1) {
      case 0:
        qz[k] = qy[0];
        break;
      default: {
        const CCTK_REAL diffleft = qy[0] + qy[2] - 2.0 * qy[1];
        const CCTK_REAL diffright = qy[1] + qy[3] - 2.0 * qy[2];

        CCTK_REAL res;
        if (fabs(diffleft) < fabs(diffright))
          res = coeffsLeft[1][0] * qy[0] + coeffsLeft[1][1] * qy[1] +
                coeffsLeft[1][2] * qy[2];
        else
          res = coeffsRight[1][0] * qy[1] + coeffsRight[1][1] * qy[2] +
                coeffsRight[1][2] * qy[3];

        if (diffleft * diffright <= 0.0 ||
            (res - qy[1]) * (qy[2] - res) < 0.0) {
          res = coeffs[1][0] * qy[1] + coeffs[1][1] * qy[2];
        }

        qz[k] = res;
        break;
      }
      }
    }

    switch (O2) {
    case 0:
      vals[v] = qz[0];
      break;
    default: {
      const CCTK_REAL diffleft = qz[0] + qz[2] - 2.0 * qz[1];
      const CCTK_REAL diffright = qz[1] + qz[3] - 2.0 * qz[2];

      CCTK_REAL res;
      if (fabs(diffleft) < fabs(diffright))
        res = coeffsLeft[2][0] * qz[0] + coeffsLeft[2][1] * qz[1] +
              coeffsLeft[2][2] * qz[2];
      else
        res = coeffsRight[2][0] * qz[1] + coeffsRight[2][1] * qz[2] +
              coeffsRight[2][2] * qz[3];

      if (diffleft * diffright <= 0.0 || (res - qz[1]) * (qz[2] - res) < 0.0) {
        res = coeffs[2][0] * qz[1] + coeffs[2][1] * qz[2];
      }

      vals[v] = res;
      break;
    }
    }

  } // for v

#ifdef CARPETINTERP2_CHECK
#if 0
    for (size_t v=0; v<varptrs.size(); ++v) {
      vals[v] = ind[0] + 1000 * (ind[1] + 1000 * ind[2]);
    } // for v
#endif
#if 0
    for (size_t v=0; v<varptrs.size(); ++v) {
      vals[v] = pn.p * 1000000 + pn.n;
    } // for v
#endif
#endif
}

// Interpolate a set of variables at the same location, with a given
// set of interpolation coefficients.  This calls the specialised
// interpolation function, depending on whether interpolation is
// exact in each of the directions.
template <int O>
void fasterp_eno2_src_loc_t::interpolate(
    ivect const &ash,
#ifdef CARPETINTERP2_CHECK
    ivect const &lsh,
#endif
    vector<CCTK_REAL const *> const &varptrs,
    CCTK_REAL *restrict const vals) const {
  int const Z = 0;
  if (exact[2]) {
    if (exact[1]) {
      if (exact[0]) {
        interpolate<Z, Z, Z>(ash, CI2C(lsh, ) varptrs, vals);
      } else {
        interpolate<O, Z, Z>(ash, CI2C(lsh, ) varptrs, vals);
      }
    } else {
      if (exact[0]) {
        interpolate<Z, O, Z>(ash, CI2C(lsh, ) varptrs, vals);
      } else {
        interpolate<O, O, Z>(ash, CI2C(lsh, ) varptrs, vals);
      }
    }
  } else {
    if (exact[1]) {
      if (exact[0]) {
        interpolate<Z, Z, O>(ash, CI2C(lsh, ) varptrs, vals);
      } else {
        interpolate<O, Z, O>(ash, CI2C(lsh, ) varptrs, vals);
      }
    } else {
      if (exact[0]) {
        interpolate<Z, O, O>(ash, CI2C(lsh, ) varptrs, vals);
      } else {
        interpolate<O, O, O>(ash, CI2C(lsh, ) varptrs, vals);
      }
    }
  }
}

// Interpolate a set of variables at the same location, with a given
// set of interpolation coefficients.  This calls a specialised
// interpolation function, depending on whether interpolation is
// exact in each of the directions.
void fasterp_eno2_src_loc_t::interpolate(
    ivect const &ash,
#ifdef CARPETINTERP2_CHECK
    ivect const &lsh,
#endif
    int const order, vector<CCTK_REAL const *> const &varptrs,
    CCTK_REAL *restrict const vals) const {
  switch (order) {
  case 2:
    interpolate<2>(ash, CI2C(lsh, ) varptrs, vals);
    break;
  default:
    // Add higher orders here as desired
    CCTK_ERROR("Interpolation orders other than 2 don't make sense for eno2");
    assert(0);
  }
}

void fasterp_eno2_src_loc_t::output(ostream &os) const {
  os << "fasterp_src_loc_t{";
  os << "coeffs=[";
  for (int d = 0; d < dim; ++d) {
    if (d > 0)
      os << ",";
    os << "[";
    for (int n = 0; n <= 1; ++n) {
      if (n > 0)
        os << ",";
      os << coeffs[d][n];
    }
    os << "]";
  }
  os << "],";
  os << "exact=" << exact << ",";
#ifdef CARPETINTERP2_CHECK
  os << "pn=" << pn << ",";
  os << "mrc=" << mrc << ",";
  os << "ipos=" << ipos << ",";
  os << "ind=" << ind << ",";
#endif
  os << "ind3d=" << ind3d;
#ifdef CARPETINTERP2_CHECK
  os << ","
     << "saved_ash=" << saved_ash;
#endif
  os << "}";
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

// Set up an interpolation starting from global coordinates
template <typename FASTERP>
fasterp_setup_gen_t<FASTERP>::fasterp_setup_gen_t(
    cGH const *restrict const cctkGH, fasterp_glocs_t const &locations,
    int const order_)
    : order(order_), reflevel(Carpet::reflevel),
      regridding_epoch(reflevel == -1
                           ? Carpet::regridding_epoch
                           : Carpet::level_regridding_epochs.AT(reflevel)) {
  // Some global properties
  int const npoints = locations.size();

  // Calculate patch numbers and local coordinates
  fasterp_llocs_t local_locations(npoints);

  if (CCTK_IsFunctionAliased("MultiPatch_GlobalToLocal")) {
    // This is a multi-patch simulation: convert global to local
    // coordinates

    CCTK_REAL const *coords[dim];
    CCTK_REAL *local_coords[dim];
    for (int d = 0; d < dim; ++d) {
      coords[d] = &locations.coords[d].front();
      local_coords[d] = &local_locations.coords[d].front();
    }

    int const ierr = MultiPatch_GlobalToLocal(cctkGH, dim, npoints, coords,
                                              &local_locations.maps.front(),
                                              local_coords, NULL, NULL);
    assert(not ierr);

  } else {
    // This is a single-patch simulation

    // TODO: Optimise this: don't copy coordinates, and don't store
    // map numbers if there is only one map.
    for (int n = 0; n < npoints; ++n) {
      int const m = 0;
      local_locations.maps.AT(n) = m;
      for (int d = 0; d < dim; ++d) {
        local_locations.coords[d].AT(n) = locations.coords[d].AT(n);
      }
    }

  } // if not multi-patch

  setup(cctkGH, local_locations);
}

// Set up an interpolation starting from local coordinates
template <typename FASTERP>
fasterp_setup_gen_t<FASTERP>::fasterp_setup_gen_t(
    cGH const *restrict const cctkGH, fasterp_llocs_t const &locations,
    int const order_)
    : order(order_) {
  setup(cctkGH, locations);
}

// Helper for setting up an interpolation
template <typename FASTERP>
void fasterp_setup_gen_t<FASTERP>::setup(cGH const *restrict const cctkGH,
                                         fasterp_llocs_t const &locations) {
  DECLARE_CCTK_PARAMETERS;

  if (verbose)
    CCTK_VInfo(CCTK_THORNSTRING, "Setting up interpolation for %d grid points",
               int(locations.size()));

  if (order < 0) {
    CCTK_ERROR("Interpolation order must be non-negative");
  }
  if (order > max_order) {
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                "Interpolation order cannot be larger than max_order=%d; "
                "order=%d was requested.  "
                "(You can increase the compile time constant max_order "
                "in thorn CarpetInterp2.)",
                max_order, order);
  }

  // Some global properties
  int const npoints = locations.size();
  int const nprocs = CCTK_nProcs(cctkGH);
  // int const myproc = CCTK_MyProc (cctkGH);

  mrc_t::determine_mrc_info();
  int const maxmrc = mrc_t::get_max_ind();

  MPI_Comm const &comm_world = *(MPI_Comm const *)GetMPICommWorld(cctkGH);

  // Obtain the coordinate ranges for all patches
  vector<rvect> lower(Carpet::maps);
  vector<rvect> upper(Carpet::maps);
  vector<rvect> delta(Carpet::maps);  // spacing on finest possible grid
  vector<rvect> idelta(Carpet::maps); // inverse spacing
  for (int m = 0; m < Carpet::maps; ++m) {
    jvect gsh;
    int const ierr =
        GetCoordRange(cctkGH, m, Carpet::mglevel, dim, &gsh[0], &lower.AT(m)[0],
                      &upper.AT(m)[0], &delta.AT(m)[0]);
    assert(not ierr);
    // delta.AT(m) /= Carpet::maxspacereflevelfact;
    gh const *const hh = Carpet::vhh.AT(m);
    ibbox const &baseext = hh->baseextent(Carpet::mglevel, 0);
    delta.AT(m) /= baseext.stride();
    idelta.AT(m) = CCTK_REAL(1.0) / delta.AT(m);
    if (veryverbose) {
      cout << "GetCoordRange[" << m << "]: lower=" << lower.AT(m)
           << " upper=" << upper.AT(m) << " delta=" << delta.AT(m) << endl;
    }
  }

  // Calculate refinement levels, components, and integer grid point
  // indices
  if (verbose)
    CCTK_INFO("Mapping points onto components");
  vector<fasterp_iloc_t> ilocs(npoints);
  vector<int> proc(npoints);
  fill_with_poison(proc);
  vector<int> nlocs(nprocs, 0);
  int min_rl, max_rl;
  if (Carpet::is_level_mode()) {
    min_rl = Carpet::reflevel;
    max_rl = Carpet::reflevel + 1;
  } else if (Carpet::is_global_mode()) {
    min_rl = 0;
    max_rl = Carpet::reflevels;
  }
#pragma omp parallel for
  for (int n = 0; n < npoints; ++n) {
    int const m = locations.maps.AT(n);
    rvect const pos(locations.coords[0].AT(n), locations.coords[1].AT(n),
                    locations.coords[2].AT(n));

    gh const *const hh = Carpet::vhh.AT(m);
    // ibbox const & baseext = hh->baseextent(Carpet::mglevel, 0);
    // rvect const rpos =
    //   (pos - lower.AT(m)) / (upper.AT(m) - lower.AT(m)) *
    //   rvect (baseext.upper() - baseext.lower());
    rvect const rpos = (pos - lower.AT(m)) * idelta.AT(m);

    // Find refinement level and component
    int rl, c;
    ivect ipos;
    hh->locate_position(rpos, Carpet::mglevel, min_rl, max_rl, rl, c, ipos);
    if (not(rl >= 0 and c >= 0)) {
#pragma omp critical
      {
        ostringstream msg;
        msg << "Interpolation point " << n << " on map " << m << " "
            << "at " << pos << " is outside of the grid hierarchy";
        msg << "\n"
            << "rl=" << rl << " c=" << c << "\n"
            << "rpos=" << rpos << "\n"
            << "ipos=" << ipos << "\n"
            << "lower=" << lower << "\n"
            << "upper=" << upper << "\n"
            << "delta=" << delta << "\n"
            << "idelta=" << idelta << "\n"
            << "hh=" << *hh << "\n";
        CCTK_ERROR(msg.str().c_str());
      }
    }
    assert(rl >= 0 and c >= 0);

    ibbox const &ext = Carpet::vdd.AT(m)
                           ->light_boxes.AT(Carpet::mglevel)
                           .AT(rl)
                           .AT(c)
                           .exterior;
    rvect dpos = rpos - rvect(ipos);

    // Convert from Carpet indexing to grid point indexing
    assert(all(ipos % ext.stride() == ivect(0)));
    ipos /= ext.stride();
    dpos /= rvect(ext.stride());
    if (not(all(fabs(dpos) <= rvect(0.5)))) {
      cout << "fasterp.cc:659\n"
           << "   dpos=" << dpos << "\n"
           << "   ext=" << ext << "\n";
    }
    assert(all(fabs(dpos) <= rvect(0.5)));

    ivect const ind = ipos - ext.lower() / ext.stride();
    ivect const ash = pad_shape(ext);
    int const ind3d = index(ash, ind);
#if 0
      ENTER_SINGLEMAP_MODE (cctkGH, m, CCTK_GF) {
        ENTER_LOCAL_MODE (cctkGH, c, CCTK_GF) {
          CCTK_REAL const * restrict const xptr = (CCTK_REAL const *) CCTK_VarDataPtr (cctkGH, 0, "grid::x");
          CCTK_REAL const * restrict const yptr = (CCTK_REAL const *) CCTK_VarDataPtr (cctkGH, 0, "grid::y");
          CCTK_REAL const * restrict const zptr = (CCTK_REAL const *) CCTK_VarDataPtr (cctkGH, 0, "grid::z");
          assert (xptr);
          assert (yptr);
          assert (zptr);
          cout << "CI2 map=" << m << " pos=" << pos << " ind=" << ind << " x=" << xptr[ind3d] << " y=" << yptr[ind3d] << " z=" << zptr[ind3d] << endl;
        } LEAVE_LOCAL_MODE;
      } LEAVE_SINGLEMAP_MODE;
#endif

    // TODO: assert that there are enough ghost zones

    // TODO: store for every face/direction/component/map/reflevel
    // how wide the boundaries are in every direction. Then check
    // against this when the stencils are set up.

    // Store result
    fasterp_iloc_t &iloc = ilocs.AT(n);
    iloc.mrc = mrc_t(m, rl, c);
#ifdef CARPETINTERP2_CHECK
    iloc.pn.p = dist::rank();
    iloc.pn.n = n;
    iloc.ipos = ipos * ext.stride();
    iloc.ind = ind;
#endif
    iloc.ind3d = ind3d;
    iloc.offset = dpos;

    // Find source processor
    int const p = Carpet::vhh.AT(m)->processor(rl, c);
    proc.AT(n) = p;
    int &nlocs_p = nlocs.AT(p);
#pragma omp atomic
    ++nlocs_p;

    // Output
    if (veryverbose) {
#pragma omp critical
      { cout << "Point #" << n << " at " << pos << ": iloc " << iloc << endl; }
    }
  }

  // Find mapping from processors to "processor indices": It may be
  // too expensive to store data for all processors, so we store
  // only data about those processors with which we actually
  // communicate.
  if (verbose)
    CCTK_INFO("Determine set of communicating processors");
  {
    int n_nz_nlocs = 0;
    for (int p = 0; p < nprocs; ++p) {
      if (nlocs.AT(p) > 0) {
        ++n_nz_nlocs;
      }
    }
    recv_descr.procs.resize(n_nz_nlocs);
    recv_descr.procinds.resize(nprocs, -1);
    int pp = 0;
    int offset = 0;
    for (int p = 0; p < nprocs; ++p) {
      if (nlocs.AT(p) > 0) {
        recv_descr.procs.AT(pp).p = p;
        recv_descr.procs.AT(pp).offset = offset;
        recv_descr.procs.AT(pp).npoints = nlocs.AT(p);
        recv_descr.procinds.AT(p) = pp;
        ++pp;
        offset += nlocs.AT(p);
      }
    }
    assert(pp == n_nz_nlocs);
    assert(offset == npoints);
    recv_descr.npoints = npoints;
  }

  // Create a mapping "index" from location index, as specified by
  // the user, to the index order in which the data are received
  // from the other processors.
  if (verbose) {
    CCTK_INFO("Determine inter-processor gather index permutation");
  }
  {
    vector<int> index(recv_descr.procs.size());
    fill_with_poison(index);
    for (int pp = 0; pp < int(recv_descr.procs.size()); ++pp) {
      index.AT(pp) = recv_descr.procs.AT(pp).offset;
    }
    recv_descr.index.resize(npoints);
    // TODO: parallelise with OpenMP
    for (int n = 0; n < npoints; ++n) {
      int const p = proc.AT(n);
      int const pp = recv_descr.procinds.AT(p);
      assert(pp >= 0);
      recv_descr.index.AT(n) = index.AT(pp);
      ++index.AT(pp);
    }
    for (int pp = 0; pp < int(recv_descr.procs.size()); ++pp) {
      int const recv_npoints = index.AT(pp) - recv_descr.procs.AT(pp).offset;
      assert(recv_npoints == recv_descr.procs.AT(pp).npoints);
    }
#ifndef NDEBUG
    vector<bool> received(npoints);
    fill(received, false);
    // NOTE: Can't use OMP parallel here -- vector<bool> has the
    // wrong memory layout for this
    for (int n = 0; n < npoints; ++n) {
      assert(not received.AT(n));
      received.AT(n) = true;
    }
    for (int n = 0; n < npoints; ++n) {
      assert(received.AT(n));
    }
#endif
  }

  // Count the number of points which have to be sent to other
  // processors, and exchange this information with MPI
  if (verbose) {
    CCTK_INFO("Count and exchange number of communicated grid points");
  }
  vector<int> recv_npoints(nprocs), recv_offsets(nprocs);
  fill(recv_npoints, 0);
  fill_with_poison(recv_offsets);
  {
    int offset = 0;
    for (int pp = 0; pp < int(recv_descr.procs.size()); ++pp) {
      int const p = recv_descr.procs.AT(pp).p;
      recv_npoints.AT(p) = recv_descr.procs.AT(pp).npoints;
      recv_offsets.AT(p) = offset;
      offset += recv_npoints.AT(p);
    }
    assert(offset == npoints);
  }
  vector<int> send_npoints(nprocs), send_offsets(nprocs);
  fill_with_poison(send_npoints);
  fill_with_poison(send_offsets);
  MPI_Alltoall(&recv_npoints.front(), 1, MPI_INT, &send_npoints.front(), 1,
               MPI_INT, comm_world);
  int npoints_send = 0;
  for (int p = 0; p < nprocs; ++p) {
    send_offsets.AT(p) = npoints_send;
    npoints_send += send_npoints.AT(p);
  }

  // Scatter the location information into a send-array, and
  // exchange it with MPI
  vector<fasterp_iloc_t> scattered_ilocs(npoints);
#pragma omp parallel for
  for (int n = 0; n < npoints; ++n) {
    scattered_ilocs.AT(recv_descr.index.AT(n)) = ilocs.AT(n);
  }
  vector<fasterp_iloc_t> gathered_ilocs(npoints_send);
  MPI_Alltoallv(&scattered_ilocs.front(), &recv_npoints.front(),
                &recv_offsets.front(), fasterp_iloc_t::mpi_datatype(),
                &gathered_ilocs.front(), &send_npoints.front(),
                &send_offsets.front(), fasterp_iloc_t::mpi_datatype(),
                comm_world);

#ifdef CARPETINTERP2_CHECK
  // Ensure that the ilocs were sent to the right processor
  for (int p = 0; p < nprocs; ++p) {
#pragma omp parallel for
    for (int n = 0; n < send_npoints.at(p); ++n) {
      fasterp_iloc_t const &iloc = gathered_ilocs.at(send_offsets.at(p) + n);
      mrc_t const &mrc = iloc.mrc;
      int const m = mrc.m;
      int const rl = mrc.rl;
      int const c = mrc.c;
      assert(Carpet::vhh.AT(m)->is_local(rl, c));
    }
  }
#endif

  // Fill in send descriptors
  send_descr.npoints = npoints_send;

  {
    int n_nz_send_npoints = 0;
    for (int p = 0; p < nprocs; ++p) {
      if (send_npoints.AT(p) > 0)
        ++n_nz_send_npoints;
    }
    send_descr.procs.resize(n_nz_send_npoints);
    int pp = 0;
    for (int p = 0; p < nprocs; ++p) {
      if (send_npoints.AT(p) > 0) {
        send_proc_t<FASTERP> &send_proc = send_descr.procs.AT(pp);
        send_proc.p = p;
        send_proc.offset = send_offsets.AT(p);
        send_proc.npoints = send_npoints.AT(p);
        ++pp;
      }
    }
    assert(pp == n_nz_send_npoints);
  }

  // Calculate stencil coefficients
  if (verbose)
    CCTK_INFO("Calculate stencil coefficients");

  for (size_t pp = 0; pp < send_descr.procs.size(); ++pp) {
    send_proc_t<FASTERP> &send_proc = send_descr.procs.AT(pp);

    vector<int> mrc2comp(maxmrc, -1);
    vector<int> comp2mrc(maxmrc);
    fill_with_poison(comp2mrc);
    int ncomps = 0;

    vector<int> npoints_comp(maxmrc);
    fill(npoints_comp, 0);

    // TODO: parallelise with OpenMP
    for (int n = 0; n < int(send_proc.npoints); ++n) {
      fasterp_iloc_t const &iloc = gathered_ilocs.AT(send_proc.offset + n);
      int const mrc = iloc.mrc.get_ind();
      if (mrc2comp.AT(mrc) == -1) {
        mrc2comp.AT(mrc) = ncomps;
        comp2mrc.AT(ncomps) = mrc;
        ++ncomps;
      }
      ++npoints_comp.AT(mrc);
    }
    assert(ncomps <= maxmrc);
    send_proc.comps.resize(ncomps);

    int offset = 0;
    for (int comp = 0; comp < ncomps; ++comp) {
      send_comp_t<FASTERP> &send_comp = send_proc.comps.AT(comp);
      int const mrc = comp2mrc.AT(comp);
      send_comp.mrc = mrc;
      send_comp.locs.reserve(npoints_comp.AT(mrc));

      mrc_t const themrc(mrc);
      int const m = themrc.m;
      int const rl = themrc.rl;
      int const c = themrc.c;
      assert(Carpet::vhh.AT(m)->is_local(rl, c));
      ibbox const &ext = Carpet::vdd.AT(m)
                             ->light_boxes.AT(Carpet::mglevel)
                             .AT(rl)
                             .AT(c)
                             .exterior;
      send_comp.ash = pad_shape(ext);
#ifdef CARPETINTERP2_CHECK
      send_comp.lsh = ext.shape() / ext.stride();
#endif

      send_comp.offset = offset;
      send_comp.npoints = npoints_comp.AT(mrc);
      offset += send_comp.npoints;
    }
    assert(offset == send_proc.npoints);

    send_proc.index.resize(send_proc.npoints);
    fill_with_poison(send_proc.index);
    // TODO: This is not parallel!  Loop over comps instead?
    // #pragma omp parallel for
    for (int n = 0; n < int(send_proc.npoints); ++n) {
      fasterp_iloc_t const &iloc = gathered_ilocs.AT(send_proc.offset + n);
      int const mrc = iloc.mrc.get_ind();
      int const comp = mrc2comp.AT(mrc);
      send_comp_t<FASTERP> &send_comp = send_proc.comps.AT(comp);

      send_proc.index.AT(n) = send_comp.offset + send_comp.locs.size();

      // fasterp_src_loc_t sloc;
      FASTERP sloc;
      int const ierr =
          sloc.calc_stencil(iloc, send_comp.ash, CI2C(send_comp.lsh, ) order);
      if (ierr) {
        CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                    "Could not determine valid interpolation stencil for point "
                    "%d on map %d, refinement level %d, component %d",
                    n, iloc.mrc.m, iloc.mrc.rl, iloc.mrc.c);
      }
      send_comp.locs.push_back(sloc);
    }

    for (int comp = 0; comp < ncomps; ++comp) {
      send_comp_t<FASTERP> &send_comp = send_proc.comps.AT(comp);
      assert(int(send_comp.locs.size()) == send_comp.npoints);
    }

#ifndef NDEBUG
    {
      vector<bool> used(send_proc.npoints, false);
      for (int n = 0; n < int(send_proc.npoints); ++n) {
        assert(not used.AT(send_proc.index.AT(n)));
        used.AT(send_proc.index.AT(n)) = true;
      }
      for (int n = 0; n < int(send_proc.npoints); ++n) {
        assert(used.AT(send_proc.index.AT(n)));
      }
    }
#endif

#ifdef CARPETINTERP2_CHECK
    for (int comp = 0; comp < ncomps; ++comp) {
      send_comp_t<FASTERP> const &send_comp = send_proc.comps.at(comp);
      assert(int(send_comp.locs.size()) == send_comp.npoints);
      for (int n = 0; n < send_comp.npoints; ++n) {
        FASTERP const &sloc = send_comp.locs.at(n);
        assert(sloc.mrc == send_comp.mrc);
        assert(all(sloc.saved_ash == send_comp.ash));
      }
    }
#endif

  } // for pp

#ifdef CARPETINTERP2_CHECK
  {
    if (verbose)
      CCTK_INFO("Compare send and receive counts");

    vector<int> recv_count(dist::size());
    fill(recv_count, 0);
    for (size_t pp = 0; pp < recv_descr.procs.size(); ++pp) {
      recv_proc_t const &recv_proc = recv_descr.procs.AT(pp);
      assert(recv_count.AT(recv_proc.p) == 0);
      assert(recv_proc.npoints > 0);
      recv_count.AT(recv_proc.p) = recv_proc.npoints;
    }

    vector<int> send_count(dist::size());
    fill(send_count, 0);
    for (size_t pp = 0; pp < send_descr.procs.size(); ++pp) {
      send_proc_t<FASTERP> const &send_proc = send_descr.procs.AT(pp);
      assert(send_count.AT(send_proc.p) == 0);
      assert(send_proc.npoints > 0);
      send_count.AT(send_proc.p) = send_proc.npoints;
    }

    {
      vector<int> tmp_count(dist::size());
      MPI_Alltoall(&send_count.front(), 1, MPI_INT, &tmp_count.front(), 1,
                   MPI_INT, dist::comm());
      swap(send_count, tmp_count);
    }
    bool error = false;
    for (int p = 0; p < dist::size(); ++p) {
      if (not(send_count.AT(p) == recv_count.AT(p))) {
        ostringstream buf;
        buf << "Error: processor " << p << " "
            << "sends me " << send_count.AT(p) << " points, "
            << "while I receive " << recv_count.AT(p) << " from it";
        CCTK_WARN(CCTK_WARN_ALERT, buf.str().c_str());
        error = true;
      }
    }
    if (error) {
      CCTK_ERROR("Internal error in interpolation setup");
    }
  }
#endif

  if (verbose)
    CCTK_INFO("Done.");
}

// Free the setup for one interpolation
template <typename FASTERP>
fasterp_setup_gen_t<FASTERP>::~fasterp_setup_gen_t() {
  // do nothing -- C++ calls destructors automatically
}

// Interpolate
template <typename FASTERP>
void fasterp_setup_gen_t<FASTERP>::interpolate(
    cGH const *restrict const cctkGH, vector<int> const &varinds,
    vector<CCTK_REAL *> &values) const {
  DECLARE_CCTK_PARAMETERS;

  // Check regridding epoch
  if (outofdate()) {
    if (reflevel == -1) {
      CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                  "The Carpet grid structure was changed since this "
                  "fasterp_setup was created");
    } else {
      CCTK_VError(
          __LINE__, __FILE__, CCTK_THORNSTRING,
          "The Carpet grid structure on level %d was changed since this "
          "fasterp_setup was created",
          reflevel);
    }
  }

  // Desired time level
  int const tl = 0; // current time

  size_t const nvars = varinds.size();
  if (verbose)
    CCTK_VInfo(CCTK_THORNSTRING, "Interpolating %d variables", int(nvars));
  assert(values.size() == nvars);

  if (nvars == 0)
    return;

  if (interp_barrier) {
    static Timers::Timer barrier_timer("PreBarrier");
    barrier_timer.start();
    CCTK_Barrier(cctkGH);
    barrier_timer.stop();
  }

  for (size_t v = 0; v < values.size(); ++v) {
    int const vi = varinds.AT(v);
    assert(vi >= 0 and vi <= CCTK_NumVars());
    int const gi = CCTK_GroupIndexFromVarI(vi);
    assert(gi >= 0);
    cGroup group;
    check(not CCTK_GroupData(gi, &group));
    assert(group.grouptype == CCTK_GF);
    assert(group.vartype == CCTK_VARIABLE_REAL);
    assert(group.dim == dim);
#if 0
      // Cannot be called in level mode
      cGroupDynamicData dyndata;
      {
        int const ierr = CCTK_GroupDynamicData (cctkGH, gi, &dyndata);
        assert (not ierr);
      }
      assert (dyndata.activetimelevels > tl);
#endif
    assert(recv_descr.npoints == 0 or values.AT(v) != NULL);
  }

  MPI_Comm const &comm_world = *(MPI_Comm const *)GetMPICommWorld(cctkGH);
  int const mpi_tag = 0;

  // Post Irecvs
  if (verbose)
    CCTK_INFO("Posting MPI_Irecvs");

  static Timers::Timer irecvs_timer("PostIrecvs");
  irecvs_timer.start();

  vector<CCTK_REAL> recv_points(recv_descr.npoints * nvars);
  fill_with_poison(recv_points);
  vector<MPI_Request> recv_reqs(recv_descr.procs.size());
#ifdef CARPETINTERP2_CHECK
  vector<pn_t> recv_pn(recv_descr.npoints);
  vector<MPI_Request> recv_reqs_pn(recv_descr.procs.size());
#endif
  for (size_t pp = 0; pp < recv_descr.procs.size(); ++pp) {
    recv_proc_t const &recv_proc = recv_descr.procs.AT(pp);

    MPI_Irecv(&recv_points.AT(recv_proc.offset * nvars),
              recv_proc.npoints * nvars, dist::mpi_datatype<CCTK_REAL>(),
              recv_proc.p, mpi_tag, comm_world, &recv_reqs.AT(pp));
#ifdef CARPETINTERP2_CHECK
    MPI_Irecv(&recv_pn.AT(recv_proc.offset), recv_proc.npoints * 2,
              dist::mpi_datatype<int>(), recv_proc.p, mpi_tag, comm_world,
              &recv_reqs_pn.AT(pp));
#endif
  }
  irecvs_timer.stop();

  // Interpolate data and post Isends
  if (verbose)
    CCTK_INFO("Interpolating and posting MPI_Isends");
  static Timers::Timer interpolate_timer("Interpolate");

  interpolate_timer.instantiate();

  // TODO: Use one array per processor?
  vector<CCTK_REAL> send_points(send_descr.npoints * nvars);
  fill_with_poison(send_points);
  vector<MPI_Request> send_reqs(send_descr.procs.size());
#ifdef CARPETINTERP2_CHECK
  vector<pn_t> send_pn(send_descr.npoints);
  vector<MPI_Request> send_reqs_pn(send_descr.procs.size());
#endif
  for (size_t pp = 0; pp < send_descr.procs.size(); ++pp) {
    send_proc_t<FASTERP> const &send_proc = send_descr.procs.AT(pp);

    vector<CCTK_REAL> computed_points(send_proc.npoints * nvars);
    fill_with_poison(computed_points);
#ifdef CARPETINTERP2_CHECK
    vector<pn_t> computed_pn(send_descr.npoints);
#endif
    for (size_t comp = 0; comp < send_proc.comps.size(); ++comp) {
      send_comp_t<FASTERP> const &send_comp = send_proc.comps.AT(comp);

      int const m = send_comp.mrc.m;
      int const rl = send_comp.mrc.rl;
      int const c = send_comp.mrc.c;

      int const lc = Carpet::vhh.AT(m)->get_local_component(rl, c);

      // Get pointers to 3D data
      vector<CCTK_REAL const *> varptrs(nvars);
      for (size_t v = 0; v < nvars; ++v) {
        int const gi = CCTK_GroupIndexFromVarI(varinds.AT(v));
        assert(gi >= 0);
        int const vi = varinds.AT(v) - CCTK_FirstVarIndexI(gi);
        assert(vi >= 0);
        ggf const *const ff = Carpet::arrdata.AT(gi).AT(m).data.AT(vi);
        void const *const ptr =
            ff->data_pointer(tl, rl, lc, Carpet::mglevel)->storage();
        varptrs.AT(v) = (CCTK_REAL const *)ptr;
        assert(varptrs.AT(v));
      }

      // TODO: This loops seems unbalanced.  Maybe the different
      // interpolations have different costs.
      interpolate_timer.start();
#pragma omp parallel for schedule(dynamic, 1000)
      for (int n = 0; n < int(send_comp.locs.size()); ++n) {
        size_t const ind = (send_comp.offset + n) * nvars;
        send_comp.locs.AT(n).interpolate(send_comp.ash,
                                         CI2C(send_comp.lsh, ) order, varptrs,
                                         &computed_points.AT(ind));
#ifdef CARPETINTERP2_CHECK
        computed_pn.AT(send_comp.offset + n) = send_comp.locs.AT(n).pn;
#endif
      }
      interpolate_timer.stop();

    } // for comp

// Gather send points
#pragma omp parallel for
    for (int n = 0; n < int(send_proc.npoints); ++n) {
      size_t const nn = send_proc.index.AT(n);
      for (size_t v = 0; v < nvars; ++v) {
        send_points.AT((send_proc.offset + n) * nvars + v) =
            computed_points.AT(nn * nvars + v);
      }
#ifdef CARPETINTERP2_CHECK
      send_pn.AT(send_proc.offset + n) = computed_pn.AT(nn);
#endif
    }
// TODO: Use a scatter index list instead of a gather index
// list, and combine the scattering with the calculation

#ifdef CARPETINTERP2_CHECK
#pragma omp parallel for
    for (int n = 0; n < int(send_proc.npoints); ++n) {
      assert(send_pn.AT(send_proc.offset + n).p == send_proc.p);
      if (n > 0) {
        assert(send_pn.AT(send_proc.offset + n).n >
               send_pn.AT(send_proc.offset + n - 1).n);
      }
    }
#endif

    MPI_Isend(&send_points.AT(send_proc.offset * nvars),
              send_proc.npoints * nvars, dist::mpi_datatype<CCTK_REAL>(),
              send_proc.p, mpi_tag, comm_world, &send_reqs.AT(pp));
#ifdef CARPETINTERP2_CHECK
    MPI_Isend(&send_pn.AT(send_proc.offset), send_proc.npoints * 2,
              dist::mpi_datatype<int>(), send_proc.p, mpi_tag, comm_world,
              &send_reqs_pn.AT(pp));
#endif
  } // for pp

  // Wait for Irecvs to complete
  if (verbose)
    CCTK_INFO("Waiting for MPI_Irevcs to complete");

  static Timers::Timer waitall_ir_timer("WaitAll_Irecvs");
  waitall_ir_timer.start();
  MPI_Waitall(recv_reqs.size(), &recv_reqs.front(), MPI_STATUSES_IGNORE);
#ifdef CARPETINTERP2_CHECK
  MPI_Waitall(recv_reqs.size(), &recv_reqs_pn.front(), MPI_STATUSES_IGNORE);
#endif

  waitall_ir_timer.stop();
  // Gather data
  if (verbose)
    CCTK_INFO("Gathering data");
  static Timers::Timer gather_timer("Gather");
  gather_timer.start();

#pragma omp parallel for
  for (int n = 0; n < recv_descr.npoints; ++n) {
    size_t const nn = recv_descr.index.AT(n);
    for (size_t v = 0; v < nvars; ++v) {
      values.AT(v)[n] = recv_points.AT(nn * nvars + v);
#ifndef NDEBUG
      recv_points.AT(nn * nvars + v) = poison;
#endif
    }
#ifdef CARPETINTERP2_CHECK
    assert(recv_pn.AT(nn).p == dist::rank());
    assert(recv_pn.AT(nn).n == n);
#endif
  }

  gather_timer.stop();

  // Wait for Isends to complete
  if (verbose)
    CCTK_INFO("Waiting for MPI_Isends to complete");
  static Timers::Timer waitall_is_timer("WaitAll_Isend");
  waitall_is_timer.start();
  MPI_Waitall(send_reqs.size(), &send_reqs.front(), MPI_STATUSES_IGNORE);
#ifdef CARPETINTERP2_CHECK
  MPI_Waitall(send_reqs.size(), &send_reqs_pn.front(), MPI_STATUSES_IGNORE);
#endif

  waitall_is_timer.stop();

  if (interp_barrier) {
    static Timers::Timer barrier_timer("PostBarrier");
    barrier_timer.start();
    CCTK_Barrier(cctkGH);
    barrier_timer.stop();
  }

  if (verbose)
    CCTK_INFO("Done.");
}

// instantiate template
template class fasterp_setup_gen_t<fasterp_src_loc_t>;

template class fasterp_setup_gen_t<fasterp_eno2_src_loc_t>;

} // namespace CarpetInterp2
