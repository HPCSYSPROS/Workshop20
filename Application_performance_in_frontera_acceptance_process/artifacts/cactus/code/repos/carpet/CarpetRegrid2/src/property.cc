#include <cctk.h>
#include <cctk_Parameters.h>

#include <carpet.hh>

#include "boundary.hh"
#include "property.hh"

#include <typeinfo>

// Consistency properties for the grid structure

namespace CarpetRegrid2 {

using namespace std;
using namespace Carpet;

// Each property consists of a test, which returns true or false
// depending on whether the property is satisfied, and an action
// that enforces the property.

bool property::test(gh const &hh, dh const &dd, level_boundary const &bnd,
                    vector<ibset> const &regions, int const rl) {
  assert(rl >= 0 and rl < int(regions.size()));
  return test_impl(hh, dd, bnd, regions, rl);
}

void property::enforce(gh const &hh, dh const &dd, level_boundary const &bnd,
                       vector<ibset> &regions, int const rl) {
  DECLARE_CCTK_PARAMETERS;
  enforce_impl(hh, dd, bnd, regions, rl);
  if (not test(hh, dd, bnd, regions, rl)) {
    cout << "Property " << typeid(*this).name() << "\n";
    CCTK_WARN(CCTK_WARN_ABORT, "Property does not hold after being enforced");
  }
}

//////////////////////////////////////////////////////////////////////////////
// Ensure that this grid contains the next finer grid
//////////////////////////////////////////////////////////////////////////////

ibset proper_nesting::enlarged_fine_grid(gh const &hh, dh const &dd,
                                         level_boundary const &bnd,
                                         vector<ibset> const &regions,
                                         int const rl) {
  DECLARE_CCTK_PARAMETERS;

  assert(rl + 1 < int(regions.size()));

  // The minimum amount of space required between the boundaries of
  // this and the next finer grid. We need a certain amount of space
  // on the coarse and a certain amount on the fine grid.
  i2vect const fdistance = dd.ghost_widths.at(rl + 1);
  i2vect const cdistance =
      i2vect(min_distance + dd.prolongation_stencil_size(rl));

  ibset enlarged;

  // Loop over all bboxes that make up the next finer level
  for (ibset::const_iterator ibb = regions.at(rl + 1).begin();
       ibb != regions.at(rl + 1).end(); ++ibb) {
    ibbox const &fbb = *ibb;

    // Find out which faces are on a boundary
    bvect const lower_is_outer = fbb.lower() <= bnd.level_physical_ilower;
    bvect const upper_is_outer = fbb.upper() >= bnd.level_physical_iupper;
    b2vect const ob(lower_is_outer, upper_is_outer);

    ibbox const domext = hh.baseextent(0, rl);

    // Enlarge the bbox, first on the fine grid, then transfer it to
    // the coarse grid, then enlarge it again
    ibbox const ebb = fbb.expand(i2vect(not ob) * fdistance);
    ibbox const cbb = ebb.expanded_for(domext);
    ibbox const ecbb = cbb.expand(i2vect(not ob) * cdistance);

    // Add it
    enlarged |= ecbb;
  }

  return enlarged;
}

bool proper_nesting::test_impl(gh const &hh, const dh &dd,
                               level_boundary const &bnd,
                               vector<ibset> const &regions, int const rl) {
  // This should not be tested because it has to be applied
  // unconditionally and only once
  return true;
}

void proper_nesting::enforce_impl(gh const &hh, dh const &dd,
                                  level_boundary const &bnd,
                                  vector<ibset> &regions, int const rl) {
  DECLARE_CCTK_PARAMETERS;

  if (not ensure_proper_nesting)
    return;
  if (rl == int(regions.size()) - 1)
    return;

  if (veryverbose) {
    cout << "Refinement level " << rl << ": ensuring proper nesting...\n";
  }

  // Enlarge the level
  regions.AT(rl) |= enlarged_fine_grid(hh, dd, bnd, regions, rl);

  if (veryverbose) {
    cout << "   New regions are " << regions.at(rl) << "\n";
  }
}

//////////////////////////////////////////////////////////////////////////////
// Add buffer zones (do this only once)
//////////////////////////////////////////////////////////////////////////////

ibset add_buffers::buffered_regions(gh const &hh, dh const &dd,
                                    level_boundary const &bnd,
                                    vector<ibset> const &regions,
                                    int const rl) {
  return regions.at(rl)
      .expand(dd.buffer_widths.at(rl) + dd.overlap_widths.at(rl));
}

bool add_buffers::test_impl(gh const &hh, dh const &dd,
                            level_boundary const &bnd,
                            vector<ibset> const &regions, int const rl) {
  // This should not be tested because it has to be applied
  // unconditionally and only once
  return true;
}

void add_buffers::enforce_impl(gh const &hh, dh const &dd,
                               level_boundary const &bnd,
                               vector<ibset> &regions, int const rl) {
  DECLARE_CCTK_PARAMETERS;

  if (veryverbose) {
    cout << "Refinement level " << rl << ": adding buffer zones...\n";
  }

  regions.at(rl) = buffered_regions(hh, dd, bnd, regions, rl);

  if (veryverbose) {
    cout << "   New regions are " << regions.at(rl) << "\n";
  }
}

//////////////////////////////////////////////////////////////////////////////
// Combine all regions into a single region, if this is worthwhile
//////////////////////////////////////////////////////////////////////////////

ibbox combine_regions::combined_regions(gh const &hh, dh const &dd,
                                        level_boundary const &bnd,
                                        vector<ibset> const &regions,
                                        int const rl) {
  return regions.at(rl).container();
}

bool combine_regions::test_impl(gh const &hh, dh const &dd,
                                level_boundary const &bnd,
                                vector<ibset> const &regions, int const rl) {
  DECLARE_CCTK_PARAMETERS;

  if (min_fraction == 1.0)
    return true;

  ibbox const combined = combined_regions(hh, dd, bnd, regions, rl);

  CCTK_REAL const regions_size = static_cast<CCTK_REAL>(regions.at(rl).size());
  CCTK_REAL const combined_size = static_cast<CCTK_REAL>(combined.size());

  // Is the current setup "simple enough"? (It is "simple enough" if
  // either it already consists of a single box, or if using a
  // single bbox would be too inefficient.)
  // TODO: Check this also for pairs of regions
  return (regions.at(rl).setsize() == 1 or
          min_fraction * combined_size > regions_size);
}

void combine_regions::enforce_impl(gh const &hh, dh const &dd,
                                   level_boundary const &bnd,
                                   vector<ibset> &regions, int const rl) {
  DECLARE_CCTK_PARAMETERS;

  if (veryverbose) {
    cout << "Refinement level " << rl << ": combining regions...\n";
  }

  regions.at(rl) = combined_regions(hh, dd, bnd, regions, rl);

  if (veryverbose) {
    cout << "   New regions are " << regions.at(rl) << "\n";
  }
}

//////////////////////////////////////////////////////////////////////////////
// Align the boxes with the next coarser grid
//////////////////////////////////////////////////////////////////////////////

ibset snap_coarse::snapped_regions(gh const &hh, dh const &dd,
                                   level_boundary const &bnd,
                                   vector<ibset> const &regions, int const rl) {
  assert(rl > 0);

  ibbox const &base = hh.baseextent(0, rl);
  ibbox const &cbase = hh.baseextent(0, rl - 1);
  assert(all(cbase.stride() % base.stride() == 0));
  ivect const reffact = cbase.stride() / base.stride();
  i2vect const &buffers = dd.buffer_widths.at(rl);

  ibset snapped;

  for (ibset::const_iterator ibb = regions.at(rl).begin();
       ibb != regions.at(rl).end(); ++ibb) {
    ibbox const &bb = *ibb;

    // We want to align the current level (without its buffer zones)
    // with the next coarser level. Conceptually, we therefore
    // subtract the buffer zones, align, and add the buffer zones
    // again.

    // In detail, we first expand by reffact + reffact-1 - N points,
    // then contract to the next coarser grid, then expand back to
    // the current grid, and expand by N - reffact points again.
    // This sequence is correct for both vertex and cell centred
    // grids, and N is determined by the number of buffer zones.

    // N is the number of buffer zones modulo the refinement factor.
    // We cannot shrink the domain (since we cannot shrink
    // bboxsets). For alignment, only operations modulo the
    // refinement factor are relevant.

    // The additional constant reffact is necessary to ensure that
    // the coarsened boxes, which may be smaller than the fine
    // boxes, cannot become empty.

    ibbox aa = bb;
    assert(not aa.empty());
    aa = aa.expand(reffact + reffact - 1 - buffers % reffact);
    assert(not aa.empty());
    aa = aa.contracted_for(cbase);
    assert(not aa.empty());
    aa = aa.expanded_for(base);
    assert(not aa.empty());
    aa = aa.expand(buffers % reffact - reffact);
    assert(not aa.empty());

    snapped |= aa;
  }

  // We don't want to remove any points
  ibset const &original = regions.at(rl);
  // return snapped | original;
  if (not(snapped >= original)) {
    cout << "snapped=" << snapped << "\n"
         << "original=" << original << "\n"
         << "original-snapped=" << (original - snapped) << "\n";
  }
  assert(snapped >= original);

  return snapped;
}

bool snap_coarse::test_impl(gh const &hh, dh const &dd,
                            level_boundary const &bnd,
                            vector<ibset> const &regions, int const rl) {
  DECLARE_CCTK_PARAMETERS;

  if (not snap_to_coarse)
    return true;

  ibset const snapped = snapped_regions(hh, dd, bnd, regions, rl);

  // We cannot test for equality, since the difference may be
  // outside of the domain (and hence irrelevant)
  // return regions.AT(rl) == snapped;

  // Test whether any part of the difference (i.e. that part of the
  // level that would be added by snapping) is inside the domain. If
  // the difference is outside, we can safely ignore it.
  ibbox const &baseextent = hh.baseextent(0, rl);
  ibset const difference = snapped - regions.AT(rl);
  return (difference & baseextent).empty();
}

void snap_coarse::enforce_impl(gh const &hh, dh const &dd,
                               level_boundary const &bnd,
                               vector<ibset> &regions, int const rl) {
  DECLARE_CCTK_PARAMETERS;

  assert(snap_to_coarse);

  if (veryverbose) {
    cout << "Refinement level " << rl
         << ": aligning regions with next coarser grid...\n";
  }

  regions.AT(rl) = snapped_regions(hh, dd, bnd, regions, rl);

  if (veryverbose) {
    cout << "   New regions are " << regions.at(rl) << "\n";
  }
}

//////////////////////////////////////////////////////////////////////////////
// Make the boxes rotating-90 symmetric
//////////////////////////////////////////////////////////////////////////////

ibset rotsym90::symmetrised_regions(gh const &hh, dh const &dd,
                                    level_boundary const &bnd,
                                    vector<ibset> const &regions,
                                    int const rl) {
  // ibbox const& baseextent = hh.baseextent(0,rl);

  ibset symmetrised = regions.at(rl);
  for (ibset::const_iterator ibb = regions.at(rl).begin();
       ibb != regions.at(rl).end(); ++ibb) {
    ibbox const &bb = *ibb;

    bvect const lower_is_outside_lower =
        bb.lower() - bnd.min_bnd_dist_away[0] * bb.stride() <=
        bnd.level_physical_ilower;

    // Treat both x and y directions
    for (int dir = 0; dir <= 1; ++dir) {
      if (lower_is_outside_lower[dir]) {
        ivect const ilo = bb.lower();
        ivect const iup = bb.upper();
        ivect const istr = bb.stride();

        // Origin
        rvect const axis(bnd.physical_lower[0], bnd.physical_lower[1],
                         bnd.physical_lower[2]); // z component is unused
        ivect const iaxis0 = rpos2ipos(axis, bnd.origin, bnd.scale, hh, rl);
        assert(all(iaxis0 % istr == 0));
        ivect const iaxis1 = rpos2ipos1(axis, bnd.origin, bnd.scale, hh, rl);
        assert(all(iaxis1 % istr == 0));
        ivect const offset = iaxis1 - iaxis0;
        assert(all(offset % istr == 0));
        assert(all(offset >= 0 and offset < 2 * istr));
        assert(all((iaxis0 + iaxis1 - offset) % (2 * istr) == 0));
        ivect const iaxis = (iaxis0 + iaxis1 - offset) / 2;
        // negated (reflected) domain boundaries
        ivect const neg_ilo = (2 * iaxis + offset) - ilo;
        ivect const neg_iup = (2 * iaxis + offset) - iup;
        // offset to add when permuting directions
        ivect const permute01(-iaxis[0] + iaxis[1], -iaxis[1] + iaxis[0], 0);

        // Rotate 90 degrees about z axis
        ivect new_ilo, new_iup;
        if (dir == 0) {
          // rotate clockwise
          new_ilo = ivect(ilo[1], neg_iup[0], ilo[2]) + permute01;
          new_iup = ivect(iup[1], neg_ilo[0], iup[2]) + permute01;
        }
        if (dir == 1) {
          // rotate counterclockwise
          new_ilo = ivect(neg_iup[1], ilo[0], ilo[2]) + permute01;
          new_iup = ivect(neg_ilo[1], iup[0], iup[2]) + permute01;
        }
        ivect const new_istr(istr);

        ibbox const new_bb(new_ilo, new_iup, new_istr);
        // Will be clipped later
        // assert (new_bb.is_contained_in (baseextent));

        // symmetrised |= new_bb & baseextent;
        symmetrised |= new_bb;
      }
    }
  }

  return symmetrised;
}

bool rotsym90::test_impl(gh const &hh, dh const &dd, level_boundary const &bnd,
                         vector<ibset> const &regions, int const rl) {
  DECLARE_CCTK_PARAMETERS;

  if (not symmetry_rotating90)
    return true;

  ibset const symmetrised = symmetrised_regions(hh, dd, bnd, regions, rl);

  // We cannot test for equality, since the difference may be
  // outside of the domain (and hence irrelevant)
  // return regions.AT(rl) == symmetrised;

  // Test whether any part of the difference (i.e. that part of the
  // level that would be added by symmetrising) is inside the
  // domain. If the difference is outside, we can safely ignore it.
  ibbox const &baseextent = hh.baseextent(0, rl);
  ibset const difference = symmetrised - regions.AT(rl);
  return (difference & baseextent).empty();
}

void rotsym90::enforce_impl(gh const &hh, dh const &dd,
                            level_boundary const &bnd, vector<ibset> &regions,
                            int const rl) {
  DECLARE_CCTK_PARAMETERS;

  assert(symmetry_rotating90);

  if (veryverbose) {
    cout << "Refinement level " << rl
         << ": making regions rotating-90 symmetric...\n";
  }

  regions.AT(rl) = symmetrised_regions(hh, dd, bnd, regions, rl);

  if (veryverbose) {
    cout << "   New regions are " << regions.at(rl) << "\n";
  }
}

//////////////////////////////////////////////////////////////////////////////
// Make the boxes rotating-180 symmetric
//////////////////////////////////////////////////////////////////////////////

ibset rotsym180::symmetrised_regions(gh const &hh, dh const &dd,
                                     level_boundary const &bnd,
                                     vector<ibset> const &regions,
                                     int const rl) {
  ibbox const &baseextent = hh.baseextent(0, rl);

  ibset symmetrised = regions.at(rl);
  for (ibset::const_iterator ibb = regions.at(rl).begin();
       ibb != regions.at(rl).end(); ++ibb) {
    ibbox const &bb = *ibb;

    bvect const lower_is_outside_lower =
        bb.lower() - bnd.min_bnd_dist_away[0] * bb.stride() <=
        bnd.level_physical_ilower;

    // Treat x direction
    int const dir = 0;
    if (lower_is_outside_lower[dir]) {
      ivect const ilo = bb.lower();
      ivect const iup = bb.upper();
      ivect const istr = bb.stride();
      assert(istr[0] == istr[1]);

      // Origin
      assert(hh.refcent == vertex_centered or all(istr % 2 == 0));
      rvect const axis(bnd.physical_lower[0],
                       (bnd.physical_lower[1] + bnd.physical_upper[1]) / 2,
                       bnd.physical_lower[2]); // z component is unused
      ivect const iaxis0 = rpos2ipos(axis, bnd.origin, bnd.scale, hh, rl);
      assert(all((iaxis0 - baseextent.lower()) % istr == 0));
      ivect const iaxis1 = rpos2ipos1(axis, bnd.origin, bnd.scale, hh, rl);
      assert(all((iaxis1 - baseextent.lower()) % istr == 0));
      ivect const offset = iaxis1 - iaxis0;
      assert(all(offset % istr == 0));
      if (hh.refcent == vertex_centered) {
        assert(all(offset >= 0 and offset < 2 * istr));
        assert(all((iaxis0 + iaxis1 - offset) % (2 * istr) == 0));
      } else {
        // The offset may be negative because both boundaries are
        // shifted inwards by 1/2 grid spacing, and therefore iaxis0
        // < iaxis1 + istr
        assert(all(offset >= -istr and offset < istr));
        assert(all((iaxis0 + iaxis1 - offset) % (2 * istr) == istr));
        assert(all(istr % 2 == 0));
      }
      ivect const iaxis = (iaxis0 + iaxis1 - offset) / 2;
      ivect const neg_ilo = (2 * iaxis + offset) - ilo;
      ivect const neg_iup = (2 * iaxis + offset) - iup;

      // Rotate 180 degrees about z axis
      ivect const new_ilo(neg_iup[0], neg_iup[1], ilo[2]);
      ivect const new_iup(neg_ilo[0], neg_ilo[1], iup[2]);
      ivect const new_istr(istr);

      ibbox const new_bb(new_ilo, new_iup, new_istr);
      // Will be clipped later
      // assert (new_bb.is_contained_in (baseextent));

      // symmetrised |= new_bb & baseextent;
      symmetrised |= new_bb;
    }
  }

  return symmetrised;
}

bool rotsym180::test_impl(gh const &hh, dh const &dd, level_boundary const &bnd,
                          vector<ibset> const &regions, int const rl) {
  DECLARE_CCTK_PARAMETERS;

  if (not symmetry_rotating180)
    return true;

  ibset const symmetrised = symmetrised_regions(hh, dd, bnd, regions, rl);

  // We cannot test for equality, since the difference may be
  // outside of the domain (and hence irrelevant)
  // return regions.AT(rl) == symmetrised;

  // Test whether any part of the difference (i.e. that part of the
  // level that would be added by symmetrising) is inside the
  // domain. If the difference is outside, we can safely ignore it.
  ibbox const &baseextent = hh.baseextent(0, rl);
  ibset const difference = symmetrised - regions.AT(rl);
  return (difference & baseextent).empty();
}

void rotsym180::enforce_impl(gh const &hh, dh const &dd,
                             level_boundary const &bnd, vector<ibset> &regions,
                             int const rl) {
  DECLARE_CCTK_PARAMETERS;

  assert(symmetry_rotating180);

  if (veryverbose) {
    cout << "Refinement level " << rl
         << ": making regions rotating-180 symmetric...\n";
  }

  regions.AT(rl) = symmetrised_regions(hh, dd, bnd, regions, rl);

  if (veryverbose) {
    cout << "   New regions are " << regions.at(rl) << "\n";
  }
}

//////////////////////////////////////////////////////////////////////////////
// Make the boxes  parity symmetric
//////////////////////////////////////////////////////////////////////////////

ibset parsym::symmetrised_regions(gh const &hh, dh const &dd,
                                  level_boundary const &bnd,
                                  vector<ibset> const &regions, int const rl) {
  ibbox const &baseextent = hh.baseextent(0, rl);

  ibset symmetrised = regions.at(rl);
  for (ibset::const_iterator ibb = regions.at(rl).begin();
       ibb != regions.at(rl).end(); ++ibb) {
    ibbox const &bb = *ibb;

    bvect const lower_is_outside_lower =
        bb.lower() - bnd.min_bnd_dist_away[0] * bb.stride() <=
        bnd.level_physical_ilower;

    // Treat z direction
    int const dir = 2;
    if (lower_is_outside_lower[dir]) {
      ivect const ilo = bb.lower();
      ivect const iup = bb.upper();
      ivect const istr = bb.stride();
      assert(istr[0] == istr[1]);

      // Origin
      assert(hh.refcent == vertex_centered or all(istr % 2 == 0));
      rvect const axis((bnd.physical_lower[0] + bnd.physical_upper[0]) / 2,
                       (bnd.physical_lower[1] + bnd.physical_upper[1]) / 2,
                       bnd.physical_lower[2]);
      ivect const iaxis0 = rpos2ipos(axis, bnd.origin, bnd.scale, hh, rl);
      assert(all((iaxis0 - baseextent.lower()) % istr == 0));
      ivect const iaxis1 = rpos2ipos1(axis, bnd.origin, bnd.scale, hh, rl);
      assert(all((iaxis1 - baseextent.lower()) % istr == 0));
      ivect const offset = iaxis1 - iaxis0;
      assert(all(offset % istr == 0));
      if (hh.refcent == vertex_centered) {
        assert(all(offset >= 0 and offset < 2 * istr));
        assert(all((iaxis0 + iaxis1 - offset) % (2 * istr) == 0));
      } else {
        // The offset may be negative because both boundaries are
        // shifted inwards by 1/2 grid spacing, and therefore iaxis0
        // < iaxis1 + istr
        assert(all(offset >= -istr and offset < istr));
        assert(all((iaxis0 + iaxis1 - offset) % (2 * istr) == istr));
        assert(all(istr % 2 == 0));
      }
      ivect const iaxis = (iaxis0 + iaxis1 - offset) / 2;
      ivect const neg_ilo = (2 * iaxis + offset) - ilo;
      ivect const neg_iup = (2 * iaxis + offset) - iup;

      // Rotate 180 degrees about z axis
      ivect const new_ilo(neg_iup[0], neg_iup[1], neg_iup[2]);
      ivect const new_iup(neg_ilo[0], neg_ilo[1], neg_ilo[2]);
      ivect const new_istr(istr);

      ibbox const new_bb(new_ilo, new_iup, new_istr);
      // Will be clipped later
      // assert (new_bb.is_contained_in (baseextent));

      // symmetrised |= new_bb & baseextent;
      symmetrised |= new_bb;
    }
  }

  return symmetrised;
}

bool parsym::test_impl(gh const &hh, dh const &dd, level_boundary const &bnd,
                       vector<ibset> const &regions, int const rl) {
  DECLARE_CCTK_PARAMETERS;

  if (not symmetry_parity)
    return true;

  ibset const symmetrised = symmetrised_regions(hh, dd, bnd, regions, rl);

  // We cannot test for equality, since the difference may be
  // outside of the domain (and hence irrelevant)
  // return regions.AT(rl) == symmetrised;

  // Test whether any part of the difference (i.e. that part of the
  // level that would be added by symmetrising) is inside the
  // domain. If the difference is outside, we can safely ignore it.
  ibbox const &baseextent = hh.baseextent(0, rl);
  ibset const difference = symmetrised - regions.AT(rl);
  return (difference & baseextent).empty();
}

void parsym::enforce_impl(gh const &hh, dh const &dd, level_boundary const &bnd,
                          vector<ibset> &regions, int const rl) {
  DECLARE_CCTK_PARAMETERS;

  assert(symmetry_parity);

  if (veryverbose) {
    cout << "Refinement level " << rl
         << ": making regions parity symmetric...\n";
  }

  regions.AT(rl) = symmetrised_regions(hh, dd, bnd, regions, rl);

  if (veryverbose) {
    cout << "   New regions are " << regions.at(rl) << "\n";
  }
}

//////////////////////////////////////////////////////////////////////////////
// Make the boxes periodic in one direction
//////////////////////////////////////////////////////////////////////////////

template <int dir>
ibset periodic<dir>::symmetrised_regions(gh const &hh, dh const &dd,
                                         level_boundary const &bnd,
                                         vector<ibset> const &regions,
                                         int const rl) {
  ibbox const &baseextent = hh.baseextent(0, rl);

  // We are not using level_physical_ilower and
  // level_physical_iupper here, because these are rounded in
  // opposite directions for cell centring, so that their difference
  // is smaller than the domain size
  ivect const ilower =
      rpos2ipos(bnd.physical_lower, bnd.origin, bnd.scale, hh, rl);
  ivect const iupper =
      rpos2ipos(bnd.physical_upper, bnd.origin, bnd.scale, hh, rl);
  ivect const ioffset = ivect::dir(dir) * (iupper - ilower);
  assert(all(ioffset % baseextent.stride() == 0));

  ibset symmetrised = regions.at(rl);
  for (ibset::const_iterator ibb = regions.at(rl).begin();
       ibb != regions.at(rl).end(); ++ibb) {
    ibbox const &bb = *ibb;

    // Shift boxes upwards and downwards by one period
    symmetrised |= bb.shift(+ioffset / bb.stride());
    symmetrised |= bb.shift(-ioffset / bb.stride());
  }

  return symmetrised;
}

template <int dir>
bool periodic<dir>::test_impl(gh const &hh, dh const &dd,
                              level_boundary const &bnd,
                              vector<ibset> const &regions, int const rl) {
  DECLARE_CCTK_PARAMETERS;

  ivect const symmetry_periodic(symmetry_periodic_x, symmetry_periodic_y,
                                symmetry_periodic_z);
  if (not symmetry_periodic[dir])
    return true;

  ibset const symmetrised = symmetrised_regions(hh, dd, bnd, regions, rl);

  // We cannot test for equality, since the difference may be
  // outside of the domain (and hence irrelevant)
  // return regions.AT(rl) == symmetrised;

  // Test whether any part of the difference (i.e. that part of the
  // level that would be added by symmetrising) is inside the
  // domain. If the difference is outside, we can safely ignore it.
  ibbox const &baseextent = hh.baseextent(0, rl);
  ibset const difference = symmetrised - regions.AT(rl);
  return (difference & baseextent).empty();
}

template <int dir>
void periodic<dir>::enforce_impl(gh const &hh, dh const &dd,
                                 level_boundary const &bnd,
                                 vector<ibset> &regions, int const rl) {
  DECLARE_CCTK_PARAMETERS;

  ivect const symmetry_periodic(symmetry_periodic_x, symmetry_periodic_y,
                                symmetry_periodic_z);
  assert(symmetry_periodic[dir]);

  if (veryverbose) {
    cout << "Refinement level " << rl << ": making regions periodic in the "
         << "xyz"[dir] << " direction...\n";
  }

  regions.AT(rl) = symmetrised_regions(hh, dd, bnd, regions, rl);

  if (veryverbose) {
    cout << "   New regions are " << regions.at(rl) << "\n";
  }
}

template class periodic<0>;
template class periodic<1>;
template class periodic<2>;

//////////////////////////////////////////////////////////////////////////////
// Clip at the outer boundary
//////////////////////////////////////////////////////////////////////////////

ibset boundary_clip::clipped_regions(gh const &hh, dh const &dd,
                                     level_boundary const &bnd,
                                     vector<ibset> const &regions,
                                     int const rl) {
  ibbox const &baseextent = hh.baseextent(0, rl);

  ibset clipped;
  for (ibset::const_iterator ibb = regions.at(rl).begin();
       ibb != regions.at(rl).end(); ++ibb) {
    ibbox const &bb = *ibb;

    // Clip boxes that extend outside the boundary. Enlarge boxes
    // that are inside but too close to the outer boundary.
    bvect const lower_is_outside_lower =
        bb.lower() - bnd.min_bnd_dist_away[0] * bb.stride() <=
        bnd.level_physical_ilower;
    // Remove bboxes that are completely outside.
    bvect const upper_is_outside_lower = bb.upper() < bnd.level_physical_ilower;
    // Enlarge bboxes that extend not far enough inwards.
    bvect const upper_is_almost_outside_lower =
        bb.upper() <
        bnd.level_physical_ilower + bnd.min_bnd_dist_incl[0] * bb.stride();

    // Ditto for the upper boundary.
    bvect const upper_is_outside_upper =
        bb.upper() + bnd.min_bnd_dist_away[1] * bb.stride() >=
        bnd.level_physical_iupper;
    bvect const lower_is_outside_upper = bb.lower() > bnd.level_physical_iupper;
    bvect const lower_is_almost_outside_upper =
        bb.lower() >
        bnd.level_physical_iupper - bnd.min_bnd_dist_incl[1] * bb.stride();

    // Don't assert for trivial dimensions ( e.g., 2D: (x,y,1) )
    bvect const domain_is_trivial =
        bnd.level_physical_ilower == bnd.level_physical_iupper;

    assert(not any(lower_is_almost_outside_upper and lower_is_outside_lower and
                   not domain_is_trivial));
    assert(not any(upper_is_almost_outside_lower and upper_is_outside_upper and
                   not domain_is_trivial));

    if (any(upper_is_outside_lower or lower_is_outside_upper)) {
      // The box is completely outside. Ignore it.
      continue;
    }

    if (any((lower_is_outside_lower and bnd.boundary_staggering_mismatch[0]) or
            (upper_is_outside_upper and bnd.boundary_staggering_mismatch[1]))) {
      ostringstream msg;
      msg << "Level " << rl << " of the refinement hierarchy has inconsistent "
                               "bountary staggering."
          << "  The refined region extends up to the boundary, but the "
             "staggering of the boundary is different from the staggering of "
             "the mesh refinement."
          << "  lower_is_outside_lower=" << lower_is_outside_lower
          << "  upper_is_outside_upper=" << upper_is_outside_upper
          << "  boundary_staggering_mismatch="
          << bnd.boundary_staggering_mismatch
          << "  level_physical_ilower=" << bnd.level_physical_ilower
          << "  level_physical_iupper=" << bnd.level_physical_iupper
          << "  baseextent=" << baseextent;
      CCTK_WARN(CCTK_WARN_ABORT, msg.str().c_str());
    }

    ibbox const clipped_bb(
        either(lower_is_outside_lower, bnd.level_exterior_ilower,
               either(lower_is_almost_outside_upper,
                      (bnd.level_physical_iupper -
                       bnd.min_bnd_dist_incl[1] * bb.stride()),
                      bb.lower())),
        either(upper_is_outside_upper, bnd.level_exterior_iupper,
               either(upper_is_almost_outside_lower,
                      (bnd.level_physical_ilower +
                       bnd.min_bnd_dist_incl[0] * bb.stride()),
                      bb.upper())),
        bb.stride());
    if (not clipped_bb.is_contained_in(baseextent)) {
      ostringstream msg;
      msg << "Level " << rl << " of the refinement hierarchy is not contained "
                               "in the simulation domain."
          << "  (There may be too many ghost or buffer zones.)"
          << "  One bbox is " << clipped_bb << "."
          << "  lower_is_outside_lower=" << lower_is_outside_lower
          << "  upper_is_outside_upper=" << upper_is_outside_upper
          << "  lower_is_almost_outside_upper=" << lower_is_almost_outside_upper
          << "  upper_is_almost_outside_lower=" << upper_is_almost_outside_lower
          << "  level_exterior_ilower=" << bnd.level_exterior_ilower
          << "  level_exterior_iupper=" << bnd.level_exterior_iupper
          << "  level_physical_ilower=" << bnd.level_physical_ilower
          << "  level_physical_iupper=" << bnd.level_physical_iupper
          << "  baseextent=" << baseextent;
      CCTK_WARN(CCTK_WARN_ABORT, msg.str().c_str());
    }
    assert(clipped_bb.is_contained_in(baseextent));

    clipped |= clipped_bb;
  }

  return clipped;
}

bool boundary_clip::test_impl(gh const &hh, dh const &dd,
                              level_boundary const &bnd,
                              vector<ibset> const &regions, int const rl) {
  DECLARE_CCTK_PARAMETERS;

  ibset const clipped = clipped_regions(hh, dd, bnd, regions, rl);
  return regions.AT(rl) == clipped;
}

void boundary_clip::enforce_impl(gh const &hh, dh const &dd,
                                 level_boundary const &bnd,
                                 vector<ibset> &regions, int const rl) {
  DECLARE_CCTK_PARAMETERS;

  if (veryverbose) {
    cout << "Refinement level " << rl << ": clipping at outer boundary...\n";
  }

  regions.AT(rl) = clipped_regions(hh, dd, bnd, regions, rl);

  if (veryverbose) {
    cout << "   New regions are " << regions.at(rl) << "\n";
  }
}

//////////////////////////////////////////////////////////////////////////////
// Ensure that this grid is contained in the domain
//////////////////////////////////////////////////////////////////////////////

bool in_domain::test_impl(gh const &hh, dh const &dd, level_boundary const &bnd,
                          vector<ibset> const &regions, int const rl) {
  // TODO: baseextent is the exterior of the domain (i.e. it
  // includes the boundary); we need to use the interior of the
  // domain here (i.e. baseextent minus the boundary). Such a
  // calculation is already performed in the beginning of dh.cc's
  // regrid.
  return regions.at(rl) <= hh.baseextent(0, rl);
}

void in_domain::enforce_impl(gh const &hh, dh const &dd,
                             level_boundary const &bnd, vector<ibset> &regions,
                             int const rl) {
  // There is nothing we can do here, since we can't enlarge the
  // domain
  CCTK_WARN(CCTK_WARN_ABORT, "internal error");
}

//////////////////////////////////////////////////////////////////////////////
// Ensure that this grid is in the domain, if desired
//////////////////////////////////////////////////////////////////////////////

ibset is_symmetric::symmetrised_regions(gh const &hh, dh const &dd,
                                        level_boundary const &bnd,
                                        vector<ibset> const &regions,
                                        int const rl) {
  ibset symmetrised = regions.at(rl);
  for (ibset::const_iterator ibb = regions.at(rl).begin();
       ibb != regions.at(rl).end(); ++ibb) {
    ibbox const &bb = *ibb;

    ivect const ilo = bb.lower();
    ivect const iup = bb.upper();
    ivect const istr = bb.stride();

    // Origin
    rvect const axis(bnd.physical_lower[0], bnd.physical_lower[1],
                     bnd.physical_lower[2]);
    ivect const iaxis0 = rpos2ipos(axis, bnd.origin, bnd.scale, hh, rl);
    assert(all(iaxis0 % istr == 0));
    ivect const iaxis1 = rpos2ipos1(axis, bnd.origin, bnd.scale, hh, rl);
    assert(all(iaxis1 % istr == 0));
    ivect const offset = iaxis1 - iaxis0;
    assert(all(offset % istr == 0));
    assert(all(offset >= 0 and offset < 2 * istr));
    assert(all((iaxis0 + iaxis1 - offset) % (2 * istr) == 0));
    ivect const iaxis = (iaxis0 + iaxis1 - offset) / 2;
    // negated (reflected) domain boundaries
    ivect const neg_ilo = (2 * iaxis + offset) - ilo;
    ivect const neg_iup = (2 * iaxis + offset) - iup;

    // Mirror
    ivect const new_ilo = neg_iup;
    ivect const new_iup = neg_ilo;
    ivect const new_istr(istr);

    ibbox const new_bb(new_ilo, new_iup, new_istr);

    symmetrised |= new_bb;
  }

  return symmetrised;
}

bool is_symmetric::test_impl(gh const &hh, dh const &dd,
                             level_boundary const &bnd,
                             vector<ibset> const &regions, int const rl) {
  DECLARE_CCTK_PARAMETERS;

  if (not expect_symmetric_grids)
    return true;

  ibset const symmetrised = symmetrised_regions(hh, dd, bnd, regions, rl);
  return regions.AT(rl) == symmetrised;
}

void is_symmetric::enforce_impl(gh const &hh, dh const &dd,
                                level_boundary const &bnd,
                                vector<ibset> &regions, int const rl) {
  // There is nothing we want to do here
  CCTK_WARN(CCTK_WARN_ABORT, "internal error");
}

} // namespace CarpetRegrid2
