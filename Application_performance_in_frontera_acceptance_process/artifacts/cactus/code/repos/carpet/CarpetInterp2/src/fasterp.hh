#ifndef FASTERP_HH
#define FASTERP_HH

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <vector>

#include <cctk.h>

#include <defs.hh>
#include <dist.hh>
#include <typeprops.hh>
#include <vect.hh>

#include <carpet.hh>

// Define this at all times, because otherwise out-of-bounds
// interpolations may not be detected early enough
#define CARPETINTERP2_CHECK

namespace CarpetInterp2 {

using namespace std;
using namespace Carpet;

// int const dim = 3;

// Each interpolation point descriptor requires
//    (dim * (max_order+1) + 1)
// CCTK_REAL values of memory
int const max_order = 11;

// Keep a map, refinement level, and component index together, and
// map the triple to small integers
struct mrc_t {
  static int maps;
  static int reflevels;
  static int components;

  static void determine_mrc_info();
  static int get_max_ind() CCTK_MEMBER_ATTRIBUTE_PURE {
    return maps * reflevels * components;
  }

  int m;
  int rl;
  int c;

  mrc_t() {}

  mrc_t(int const m_, int const rl_, int const c_) : m(m_), rl(rl_), c(c_) {
    assert(m >= 0 and m < maps);
    assert(rl >= 0 and rl < reflevels);
    assert(c >= 0 and c < components);
  }

  // Create a mrc from an integer index
  mrc_t(int ind);

  // Convert a mrc into an integer index
  int get_ind() const CCTK_MEMBER_ATTRIBUTE_PURE {
    return c + components * (rl + reflevels * m);
  }

  bool operator==(mrc_t const &a) const CCTK_MEMBER_ATTRIBUTE_PURE {
    return m == a.m and rl == a.rl and c == a.c;
  }

  void output(ostream &os) const;
};

inline ostream &operator<<(ostream &os, mrc_t const &mrc) {
  mrc.output(os);
  return os;
}

// Keep a processor and an integer index together
struct pn_t {
  int p;
  int n;

  pn_t() {}

  pn_t(int const p_, int const n_) : p(p_), n(n_) {
    assert(p >= 0 and p < dist::size());
    assert(n >= 0);
  }

  void output(ostream &os) const;
};

inline ostream &operator<<(ostream &os, pn_t const &pn) {
  pn.output(os);
  return os;
}

// A global location, given by its global coordinates
struct fasterp_glocs_t {
  // TODO: use vector<rvect> instead
  vector<CCTK_REAL> coords[dim];
  fasterp_glocs_t(size_t const n) {
    for (int d = 0; d < dim; ++d) {
      coords[d].resize(n);
    }
  }
  size_t size() const CCTK_MEMBER_ATTRIBUTE_PURE { return coords[0].size(); }
};

// A local location, given by map and local coordinates
struct fasterp_llocs_t {
  vector<CCTK_INT> maps;
  vector<CCTK_REAL> coords[dim];
  fasterp_llocs_t(size_t const n) {
    maps.resize(n);
    for (int d = 0; d < dim; ++d) {
      coords[d].resize(n);
    }
  }
  size_t size() const CCTK_MEMBER_ATTRIBUTE_PURE { return maps.size(); }
};

// An integer location, given by map, refinement level, and
// component
struct fasterp_iloc_t {
  mrc_t mrc; // map, refinement level, component

#ifdef CARPETINTERP2_CHECK
  pn_t pn;    // origin of this point
  ivect ipos; // closest grid point (Carpet indexing)
  ivect ind;  // closest grid point (local indexing)
#endif
  int ind3d;    // closest grid point
  rvect offset; // in terms of grid points

  static MPI_Datatype mpi_datatype() CCTK_MEMBER_ATTRIBUTE_PURE;

  void output(ostream &os) const;
};

inline ostream &operator<<(ostream &os, fasterp_iloc_t const &iloc) {
  iloc.output(os);
  return os;
}

struct fasterp_dest_loc_t {
  // ivect idx;
  int ind3d; // destination grid point index
};

/**
   This setup is tailored for standard Lagrange interpolation.
*/
class fasterp_src_loc_t {
  CCTK_REAL coeffs[dim][max_order + 1]; // interpolation coefficients
  bvect exact;

#ifdef CARPETINTERP2_CHECK
public:
  pn_t pn;    // origin of this point
  mrc_t mrc;  // map, refinement level, component
  ivect ipos; // closest grid point (Carpet indexing)
  ivect ind;  // source grid point offset
private:
#endif
  int ind3d; // source grid point offset

#ifdef CARPETINTERP2_CHECK
public:
  ivect saved_ash; // copy of ash
private:
#endif

public:
  int calc_stencil(fasterp_iloc_t const &iloc, ivect const &ash,
#ifdef CARPETINTERP2_CHECK
                   ivect const &lsh,
#endif
                   int order);
  void interpolate(ivect const &ash,
#ifdef CARPETINTERP2_CHECK
                   ivect const &lsh,
#endif
                   int order, vector<CCTK_REAL const *> const &varptrs,
                   CCTK_REAL *restrict vals) const;

private:
  template <int O>
  void interpolate(ivect const &ash,
#ifdef CARPETINTERP2_CHECK
                   ivect const &lsh,
#endif
                   vector<CCTK_REAL const *> const &varptrs,
                   CCTK_REAL *restrict vals) const;
  template <int O0, int O1, int O2>
  void interpolate(ivect const &ash,
#ifdef CARPETINTERP2_CHECK
                   ivect const &lsh,
#endif
                   vector<CCTK_REAL const *> const &varptrs,
                   CCTK_REAL *restrict vals) const;

public:
  void output(ostream &os) const;
};

inline ostream &operator<<(ostream &os, fasterp_src_loc_t const &sloc) {
  sloc.output(os);
  return os;
}

/**
   This setup is tailored for eno2 interpolation.
*/
class fasterp_eno2_src_loc_t {
  CCTK_REAL
      coeffs[dim][2]; // interpolation coefficients for first-order stencil
  CCTK_REAL coeffsLeft[dim][3];  // interpolation coefficients for left stencil
  CCTK_REAL coeffsRight[dim][3]; // interpolation coefficients for right stencil
  bvect exact;

#ifdef CARPETINTERP2_CHECK
public:
  pn_t pn;    // origin of this point
  mrc_t mrc;  // map, refinement level, component
  ivect ipos; // closest grid point (Carpet indexing)
  ivect ind;  // source grid point offset
private:
#endif
  int ind3d; // source grid point offset (computed from left stencil)

#ifdef CARPETINTERP2_CHECK
public:
  ivect saved_ash; // copy of ash
private:
#endif

public:
  int calc_stencil(fasterp_iloc_t const &iloc, ivect const &ash,
#ifdef CARPETINTERP2_CHECK
                   ivect const &lsh,
#endif
                   int /*order*/);
  void interpolate(ivect const &ash,
#ifdef CARPETINTERP2_CHECK
                   ivect const &lsh,
#endif
                   int order, vector<CCTK_REAL const *> const &varptrs,
                   CCTK_REAL *restrict vals) const;

private:
  template <int O>
  void interpolate(ivect const &ash,
#ifdef CARPETINTERP2_CHECK
                   ivect const &lsh,
#endif
                   vector<CCTK_REAL const *> const &varptrs,
                   CCTK_REAL *restrict vals) const;
  template <int O0, int O1, int O2>
  void interpolate(ivect const &ash,
#ifdef CARPETINTERP2_CHECK
                   ivect const &lsh,
#endif
                   vector<CCTK_REAL const *> const &varptrs,
                   CCTK_REAL *restrict vals) const;

public:
  void output(ostream &os) const;
};

inline ostream &operator<<(ostream &os, fasterp_eno2_src_loc_t const &sloc) {
  sloc.output(os);
  return os;
}

// A receive descriptor, describing what is received from other
// processors
struct recv_proc_t {
  int p; // sending processor
  int offset;
  int npoints; // total number of received points
};

struct recv_descr_t {
  vector<recv_proc_t> procs;
  vector<int> procinds;
  int npoints; // total number of received points

  vector<int> index; // gather index list: index[user-location] =
                     // mpi-location
};

// A send descriptor; describing what to send to other processors
template <typename FASTERP> struct send_comp_t {
  // This structure does not exist for all components -- components
  // which are not accessed are not described, making this a sparse
  // data structure.  The fields m, rl, and c identify the
  // component.
  vector<FASTERP> locs;

  mrc_t mrc; // source map, refinement level, component
  ivect ash;
#ifdef CARPETINTERP2_CHECK
  ivect lsh;
#endif
  int offset;
  int npoints;
};

template <typename FASTERP> struct send_proc_t {
  // This structure does not exist for all processors -- processors
  // with which there is no communication are not described, making
  // this a sparse data structure.  The field p contains the
  // processor number.
  vector<send_comp_t<FASTERP> > comps;

  int p; // receiving processor
  int offset;
  int npoints; // total number of sent points

  // gather index list: index[mpi-location] = calculation-location
  vector<int> index;
};

template <typename FASTERP> struct send_descr_t {
  vector<send_proc_t<FASTERP> > procs;
  // vector<int> procinds;
  int npoints; // total number of sent points
};

template <typename FASTERP> class fasterp_setup_gen_t {
  recv_descr_t recv_descr;
  send_descr_t<FASTERP> send_descr;
  int order;

  int reflevel;
  int regridding_epoch;

  void setup(cGH const *restrict cctkGH, fasterp_llocs_t const &locations);

public:
  fasterp_setup_gen_t(cGH const *restrict cctkGH,
                      fasterp_glocs_t const &locations, int order);

  fasterp_setup_gen_t(cGH const *restrict cctkGH,
                      fasterp_llocs_t const &locations, int order);

  ~fasterp_setup_gen_t();

  void interpolate(cGH const *restrict cctkGH, vector<int> const &varinds,
                   vector<CCTK_REAL *> &values) const;

  size_t get_npoints() const CCTK_MEMBER_ATTRIBUTE_PURE {
    return recv_descr.npoints;
  }

  bool outofdate() const {
    return regridding_epoch !=
           (reflevel == -1 ? Carpet::regridding_epoch
                           : Carpet::level_regridding_epochs.AT(reflevel));
  }
};

typedef fasterp_setup_gen_t<fasterp_src_loc_t> fasterp_setup_t;

typedef fasterp_setup_gen_t<fasterp_eno2_src_loc_t> fasterp_eno2_setup_t;

} // namespace CarpetInterp2

#endif // #define FASTERP_HH
