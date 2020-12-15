#ifndef BOUNDARY_HH
#define BOUNDARY_HH

#include <defs.hh>
#include <gh.hh>
#include <vect.hh>

#include <ostream>

namespace CarpetRegrid2 {

// Convert a coordinate location to an index location.  For cell
// centring, shift upwards.
ivect rpos2ipos(rvect const &rpos, rvect const &origin, rvect const &scale,
                gh const &hh, int const rl);

// Convert a coordinate location to an index location, rounding in
// the opposite manner as rpos2ipos.  For cell centring, shift
// downwards instead of upwards.
ivect rpos2ipos1(rvect const &rpos, rvect const &origin, rvect const &scale,
                 gh const &hh, int const rl);

// Convert an index location to a coordinate location
rvect ipos2rpos(ivect const &ipos, rvect const &origin, rvect const &scale,
                gh const &hh, int const rl);

// Convert an index bbox to a coordinate bbox
rbbox ibbox2rbbox(ibbox const &ib, rvect const &origin, rvect const &scale,
                  gh const &hh, int const rl);

// Snap (enlarge) a bbox to the next coarser level, if desired
ibbox snap_ibbox(ibbox const &ib, gh const &hh, int const rl);

void get_boundary_specification(jjvect &nboundaryzones, jjvect &is_internal,
                                jjvect &is_staggered, jjvect &shiftout);

void get_physical_boundary(rvect &physical_lower, rvect &physical_upper,
                           rvect &spacing);

void calculate_exterior_boundary(rvect const &physical_lower,
                                 rvect const &physical_upper,
                                 rvect &exterior_lower, rvect &exterior_upper,
                                 rvect const &spacing);

// Location and description of the outer boundary
struct domain_boundary {
  jjvect nboundaryzones, is_internal;
  jjvect is_staggered, shiftout;

  b2vect boundary_staggering_mismatch;

  rvect physical_lower, physical_upper;
  rvect spacing;

  rvect exterior_lower, exterior_upper;

  // The physical boundary
  rvect origin;
  rvect scale;

  // The location of the outermost grid points. For cell centring,
  // these are 1/2 grid spacing inside of the boundary.
  ivect physical_ilower, physical_iupper;

  domain_boundary(gh const &hh, dh const &dd);
};

struct level_boundary : public domain_boundary {
  rvect level_physical_lower;
  rvect level_physical_upper;
  rvect level_spacing;

  rvect level_exterior_lower, level_exterior_upper;

  // The location of the outermost grid points. For cell centring,
  // these are 1/2 grid spacing inside of the boundary.
  ivect level_physical_ilower;
  ivect level_physical_iupper;

  ivect level_exterior_ilower;
  ivect level_exterior_iupper;

  // The minimum necessary distance away from the outer boundary due
  // to buffer and ghost zones. This is e.g. the distance that the
  // lower boundary of a bbox has to have from the lower boundary.
  // This is in terms of grid points.
  i2vect min_bnd_dist_away;
  // The minimum necessary distance from the outer boundary due to
  // buffer and ghost zones. This is e.g. the distance that the
  // upper boundary of a bbox has to have from the lower boundary.
  // This is in terms of grid points.
  i2vect min_bnd_dist_incl;

  level_boundary(gh const &hh, dh const &dd, int rl);
};

} // namespace CarpetRegrid2

ostream &operator<<(ostream &os, CarpetRegrid2::domain_boundary const &bnd);
ostream &operator<<(ostream &os, CarpetRegrid2::level_boundary const &bnd);

#endif // #ifndef BOUNDARY_HH
