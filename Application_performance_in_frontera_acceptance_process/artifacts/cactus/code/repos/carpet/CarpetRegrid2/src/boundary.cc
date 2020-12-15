#include <cctk.h>
#include <cctk_Parameters.h>

#include <carpet.hh>

#include "boundary.hh"

namespace CarpetRegrid2 {

using namespace Carpet;

// Convert a coordinate location to an index location. For cell
// centring, shift upwards.
ivect rpos2ipos(rvect const &rpos, rvect const &origin, rvect const &scale,
                gh const &hh, int const rl) {
  ivect const istride = hh.baseextent(0, rl).stride();
  ivect const bistride = hh.baseextent(0, 0).stride();

  if (hh.refcent == cell_centered) {
    assert(all(istride % 2 == 0));
  }

  CCTK_REAL const eps = 1.0e-10; // prevent rounding errors

#if 1
  ivect const ipos =
      hh.refcent == vertex_centered
          ? ivect(floor(((rpos - origin) * scale) / rvect(istride) +
                        rvect(0.5) + rvect(eps))) *
                istride
          : ivect(floor(((rpos - origin) * scale - rvect(bistride / 2)) /
                            rvect(istride) +
                        rvect(eps))) *
                    istride +
                istride / 2 + bistride / 2;
#else
  ivect const ipos =
      hh.refcent == vertex_centered
          ? ivect(floor(((rpos - origin) * scale) / rvect(istride) +
                        rvect(0.5))) *
                istride
          : ivect(floor(((rpos - origin) * scale - rvect(bistride / 2)) /
                            rvect(istride) +
                        rvect(0.5))) *
                    istride +
                istride / 2 + bistride / 2;
#endif

  return ipos;
}

// Convert a coordinate location to an index location, rounding in
// the opposite manner as rpos2ipos. For cell centring, shift
// downwards instead of upwards.
ivect rpos2ipos1(rvect const &rpos, rvect const &origin, rvect const &scale,
                 gh const &hh, int const rl) {
  ivect const istride = hh.baseextent(0, rl).stride();
  ivect const bistride = hh.baseextent(0, 0).stride();

  if (hh.refcent == cell_centered) {
    assert(all(istride % 2 == 0));
  }

  CCTK_REAL const eps = 1.0e-10; // prevent rounding errors

#if 1
  ivect const ipos =
      hh.refcent == vertex_centered
          ? ivect(ceil(((rpos - origin) * scale) / rvect(istride) - rvect(0.5) -
                       rvect(eps))) *
                istride
          : ivect(ceil(((rpos - origin) * scale - rvect(bistride / 2)) /
                           rvect(istride) -
                       rvect(eps))) *
                    istride -
                istride / 2 + bistride / 2;
#else
  ivect const ipos =
      hh.refcent == vertex_centered
          ? ivect(
                ceil(((rpos - origin) * scale) / rvect(istride) - rvect(0.5))) *
                istride
          : ivect(ceil(((rpos - origin) * scale - rvect(bistride / 2)) /
                           rvect(istride) -
                       rvect(0.5))) *
                    istride -
                istride / 2 + bistride / 2;
#endif

  return ipos;
}

// Convert an index location to a coordinate location
rvect ipos2rpos(ivect const &ipos, rvect const &origin, rvect const &scale,
                gh const &hh, int const rl) {
  return rvect(ipos) / scale + origin;
}

// Convert an index bbox to a coordinate bbox
rbbox ibbox2rbbox(ibbox const &ib, rvect const &origin, rvect const &scale,
                  gh const &hh, int const rl) {
  rvect const zero(0);
  return rbbox(ipos2rpos(ib.lower(), origin, scale, hh, rl),
               ipos2rpos(ib.upper(), origin, scale, hh, rl),
               ipos2rpos(ib.stride(), zero, scale, hh, rl));
}

void get_boundary_specification(jjvect &nboundaryzones, jjvect &is_internal,
                                jjvect &is_staggered, jjvect &shiftout) {
  if (CCTK_IsFunctionAliased("MultiPatch_GetBoundarySpecification")) {
    assert(Carpet::map >= 0);
    CCTK_INT const ierr = MultiPatch_GetBoundarySpecification(
        Carpet::map, 2 * dim, &nboundaryzones[0][0], &is_internal[0][0],
        &is_staggered[0][0], &shiftout[0][0]);
    assert(not ierr);
  } else if (CCTK_IsFunctionAliased("GetBoundarySpecification")) {
    CCTK_INT const ierr = GetBoundarySpecification(
        2 * dim, &nboundaryzones[0][0], &is_internal[0][0], &is_staggered[0][0],
        &shiftout[0][0]);
    assert(not ierr);
  } else {
    assert(0);
  }
}

void get_physical_boundary(rvect &physical_lower, rvect &physical_upper,
                           rvect &spacing) {
  rvect interior_lower, interior_upper;
  rvect exterior_lower, exterior_upper;
  if (CCTK_IsFunctionAliased("MultiPatch_GetDomainSpecification")) {
    assert(Carpet::map >= 0);
    CCTK_INT const ierr = MultiPatch_GetDomainSpecification(
        Carpet::map, dim, &physical_lower[0], &physical_upper[0],
        &interior_lower[0], &interior_upper[0], &exterior_lower[0],
        &exterior_upper[0], &spacing[0]);
    assert(not ierr);
  } else if (CCTK_IsFunctionAliased("GetDomainSpecification")) {
    CCTK_INT const ierr = GetDomainSpecification(
        dim, &physical_lower[0], &physical_upper[0], &interior_lower[0],
        &interior_upper[0], &exterior_lower[0], &exterior_upper[0],
        &spacing[0]);
    assert(not ierr);
  } else {
    assert(0);
  }
}

void calculate_exterior_boundary(rvect const &physical_lower,
                                 rvect const &physical_upper,
                                 rvect &exterior_lower, rvect &exterior_upper,
                                 rvect const &spacing) {
  rvect interior_lower, interior_upper;
  if (CCTK_IsFunctionAliased("MultiPatch_ConvertFromPhysicalBoundary")) {
    assert(Carpet::map >= 0);
    CCTK_INT const ierr = MultiPatch_ConvertFromPhysicalBoundary(
        Carpet::map, dim, &physical_lower[0], &physical_upper[0],
        &interior_lower[0], &interior_upper[0], &exterior_lower[0],
        &exterior_upper[0], &spacing[0]);
    assert(not ierr);
  } else if (CCTK_IsFunctionAliased("ConvertFromPhysicalBoundary")) {
    CCTK_INT const ierr = ConvertFromPhysicalBoundary(
        dim, &physical_lower[0], &physical_upper[0], &interior_lower[0],
        &interior_upper[0], &exterior_lower[0], &exterior_upper[0],
        &spacing[0]);
    assert(not ierr);
  } else {
    assert(0);
  }
}

// Location and description of the outer boundary
domain_boundary::domain_boundary(gh const &hh, dh const &dd) {
  DECLARE_CCTK_PARAMETERS;

  if (veryverbose) {
    cout << "Determining domain outer boundary...\n";
  }

  get_boundary_specification(nboundaryzones, is_internal, is_staggered,
                             shiftout);

  boundary_staggering_mismatch =
      xpose((hh.refcent == vertex_centered) != (iivect(is_staggered) == 0));
  // TODO: This is too strict
  assert(all(all(not boundary_staggering_mismatch)));

  get_physical_boundary(physical_lower, physical_upper, spacing);

  // Adapt spacing for convergence level
  spacing *= ipow(CCTK_REAL(mgfact), basemglevel);

  calculate_exterior_boundary(physical_lower, physical_upper, exterior_lower,
                              exterior_upper, spacing);

  // The physical boundary
  origin = exterior_lower;
  scale = rvect(hh.baseextent(0, 0).stride()) / spacing;

  // The location of the outermost grid points. For cell centring,
  // these are 1/2 grid spacing inside of the boundary.
  physical_ilower = rpos2ipos(physical_lower, origin, scale, hh, 0);
  physical_iupper = rpos2ipos1(physical_upper, origin, scale, hh, 0);
}

level_boundary::level_boundary(gh const &hh, dh const &dd, int const rl)
    : domain_boundary(hh, dd) {
  DECLARE_CCTK_PARAMETERS;

  if (veryverbose) {
    cout << "Refinement level " << rl << ": determining outer boundary...\n";
  }

  level_physical_lower = physical_lower;
  level_physical_upper = physical_upper;
  level_spacing = spacing / rvect(hh.reffacts.at(rl));
  if (veryverbose) {
    cout << "Refinement level " << rl << ": physical coordinate boundary is at "
         << r2vect(level_physical_lower, level_physical_upper) << "\n";
    cout << "Refinement level " << rl << ": spacing is " << level_spacing
         << "\n";
  }

  calculate_exterior_boundary(level_physical_lower, level_physical_upper,
                              level_exterior_lower, level_exterior_upper,
                              level_spacing);
  if (veryverbose) {
    cout << "Refinement level " << rl << ": exterior coordinate boundary is at "
         << r2vect(level_exterior_lower, level_exterior_upper) << "\n";
  }

  ibbox const &baseextent = hh.baseextent(0, rl);
  ivect const istride = baseextent.stride();
  if (veryverbose) {
    cout << "Refinement level " << rl << ": stride is " << istride << "\n";
  }

  // This is the location of the outermost grid points. For cell
  // centring, these are 1/2 grid spacing inside of the boundary.
  level_physical_ilower = rpos2ipos(physical_lower, origin, scale, hh, rl);
  level_physical_iupper = rpos2ipos1(physical_upper, origin, scale, hh, rl);
  if (veryverbose) {
    cout << "Refinement level " << rl << ": physical boundary is at "
         << i2vect(level_physical_ilower, level_physical_iupper) << "\n";
    cout << "Refinement level " << rl
         << ": reconstructed physical coordinate boundary is at "
         << r2vect(
                ipos2rpos(level_physical_ilower -
                              (hh.refcent == cell_centered ? istride / 2 : 0),
                          origin, scale, hh, rl),
                ipos2rpos(level_physical_iupper +
                              (hh.refcent == cell_centered ? istride / 2 : 0),
                          origin, scale, hh, rl))
         << "\n";
  }

  level_exterior_ilower =
      rpos2ipos(level_exterior_lower, origin, scale, hh, rl);
  level_exterior_iupper =
      rpos2ipos1(level_exterior_upper, origin, scale, hh, rl);
  assert(all(level_exterior_ilower >= baseextent.lower()));
  assert(all(level_exterior_iupper <= baseextent.upper()));
  if (veryverbose) {
    cout << "Refinement level " << rl << ": exterior boundary is at "
         << i2vect(level_exterior_ilower, level_exterior_iupper) << "\n";
    cout << "Refinement level " << rl
         << ": reconstructed exterior coordinate boundary is at "
         << r2vect(ipos2rpos(level_exterior_ilower, origin, scale, hh, rl),
                   ipos2rpos(level_exterior_iupper, origin, scale, hh, rl))
         << "\n";
  }

  // Find the minimum necessary distance away from the outer
  // boundary due to buffer and ghost zones. This is e.g. the
  // distance that the lower boundary of a bbox has to have from the
  // lower boundary. This is in terms of grid points.
  min_bnd_dist_away = dd.ghost_widths.at(rl);
  // Find the minimum necessary distance from the outer boundary due
  // to buffer and ghost zones. This is e.g. the distance that the
  // upper boundary of a bbox has to have from the lower boundary.
  // This is in terms of grid points.
  min_bnd_dist_incl = dd.ghost_widths.at(rl);
  // TODO: The above is required only near symmetry boundaries.
}

} // namespace CarpetRegrid2

ostream &operator<<(ostream &os, CarpetRegrid2::domain_boundary const &bnd) {
  return os << "domain_boundary:{"
            << "nboundaryzones=" << bnd.nboundaryzones << ","
            << "is_internal=" << bnd.is_internal << ","
            << "is_staggered=" << bnd.is_staggered << ","
            << "shiftout=" << bnd.shiftout << ","
            << "boundary_staggering_mismatch="
            << bnd.boundary_staggering_mismatch << ","
            << "physical_lower=" << bnd.physical_lower << ","
            << "physical_upper=" << bnd.physical_upper << ","
            << "spacing=" << bnd.spacing << ","
            << "exterior_lower=" << bnd.exterior_lower << ","
            << "exterior_upper=" << bnd.exterior_upper << ","
            << "origin=" << bnd.origin << ","
            << "scale=" << bnd.scale << ","
            << "physical_ilower=" << bnd.physical_ilower << ","
            << "physical_iupper=" << bnd.physical_iupper << "}";
}

ostream &operator<<(ostream &os, CarpetRegrid2::level_boundary const &bnd) {
  return os << "level_boundary:{"
            << *static_cast<CarpetRegrid2::domain_boundary const *>(&bnd) << ","
            << "level_physical_lower=" << bnd.level_physical_lower << ","
            << "level_physical_upper=" << bnd.level_physical_upper << ","
            << "level_spacing=" << bnd.level_spacing << ","
            << "level_exterior_lower=" << bnd.level_exterior_lower << ","
            << "level_exterior_upper=" << bnd.level_exterior_upper << ","
            << "level_physical_ilower=" << bnd.level_physical_ilower << ","
            << "level_physical_iupper=" << bnd.level_physical_iupper << ","
            << "level_exterior_ilower=" << bnd.level_exterior_ilower << ","
            << "level_exterior_iupper=" << bnd.level_exterior_iupper << ","
            << "min_bnd_dist_away=" << bnd.min_bnd_dist_away << ","
            << "min_bnd_dist_incl=" << bnd.min_bnd_dist_incl << "}";
}
