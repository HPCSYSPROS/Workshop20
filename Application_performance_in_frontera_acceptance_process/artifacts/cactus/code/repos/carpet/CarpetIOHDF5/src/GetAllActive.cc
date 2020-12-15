/*@@
  @file      GetAllActive.cc
  @date      Sun Sep 19 20:34:27 EDT 2010
  @author    Erik Schnetter, Roland Haas
  @desc
             a piece of dh::regrid to compute the active region
  @enddesc
@@*/

#include <cassert>

#include <vector>

#include "cctk.h"

#include "CarpetIOHDF5.hh"

namespace CarpetIOHDF5 {

/********************************************************************
 *********************     Internal Routines     ********************
 ********************************************************************/

static void assert_error(char const *restrict const checkstring,
                         char const *restrict const file, int const line,
                         int const ml, int const rl,
                         char const *restrict const message);
static void assert_error(char const *restrict const checkstring,
                         char const *restrict const file, int const line,
                         int const ml, int const rl, int const c,
                         char const *restrict const message);
static void assert_error(char const *restrict const checkstring,
                         char const *restrict const file, int const line,
                         int const ml, int const rl, int const c, int const cc,
                         char const *restrict const message);

/********************************************************************
 *********************     Local Data   *****************************
 ********************************************************************/

// accumulates any consistency errors
static bool there_was_an_error = false;

/*@@
  @routine
  @date       Sun Sep 19 20:46:35 EDT 2010
  @author     Erik Schnetter, Roland Haas
  @desc
              Output warning for consistency checks for a refinement level.

  @enddesc
  @history
  @endhistory
  @var        checkstring
  @vdesc      the check that failed
  @vtype      char const * restrict const
  @vio        in
  @endvar
  @var        file
  @vdesc      name of source file where error occured
  @vtype      char const * restrict const
  @vio        in
  @endvar
  @var        line
  @vdesc      source line in file where error occured
  @vtype      int const
  @vio        in
  @endvar
  @var        ml
  @vdesc      meta level where error occured
  @vtype      int const
  @vio        in
  @endvar
  @var        rl
  @vdesc      refinement level where error occured
  @vtype      int const
  @vio        in
  @endvar
  @var        message
  @vdesc      error message to display
  @vtype      char const * restrict const
  @vio        in
  @endvar
  @returntype void
@@*/
static void assert_error(char const *restrict const checkstring,
                         char const *restrict const file, int const line,
                         int const ml, int const rl,
                         char const *restrict const message) {
  CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,
             "\n%s:%d:\n   [ml=%d rl=%d] The following grid structure "
             "consistency check failed:\n   %s\n   %s",
             file, line, ml, rl, message, checkstring);
  there_was_an_error = true;
}

/*@@
  @routine
  @date       Sun Sep 19 20:46:35 EDT 2010
  @author     Erik Schnetter, Roland Haas
  @desc
              Output warning for consistency checks for a single component.

  @enddesc
  @history
  @endhistory
  @var        checkstring
  @vdesc      the check that failed
  @vtype      char const * restrict const
  @vio        in
  @endvar
  @var        file
  @vdesc      name of source file where error occured
  @vtype      char const * restrict const
  @vio        in
  @endvar
  @var        line
  @vdesc      source line in file where error occured
  @vtype      int const
  @vio        in
  @endvar
  @var        ml
  @vdesc      meta level where error occured
  @vtype      int const
  @vio        in
  @endvar
  @var        rl
  @vdesc      refinement level where error occured
  @vtype      int const
  @vio        in
  @endvar
  @var        c
  @vdesc      component where error occured
  @vtype      int const
  @vio        in
  @endvar
  @var        message
  @vdesc      error message to display
  @vtype      char const * restrict const
  @vio        in
  @endvar
  @returntype void
@@*/
static void assert_error(char const *restrict const checkstring,
                         char const *restrict const file, int const line,
                         int const ml, int const rl, int const c,
                         char const *restrict const message) {
  CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,
             "\n%s:%d:\n   [ml=%d rl=%d c=%d] The following grid structure "
             "consistency check failed:\n   %s\n   %s",
             file, line, ml, rl, c, message, checkstring);
  there_was_an_error = true;
}

/*@@
  @routine
  @date       Sun Sep 19 20:46:35 EDT 2010
  @author     Erik Schnetter, Roland Haas
  @desc
              Output warning for consistency checks between two components.

  @enddesc
  @history
  @endhistory
  @var        checkstring
  @vdesc      the check that failed
  @vtype      char const * restrict const
  @vio        in
  @endvar
  @var        file
  @vdesc      name of source file where error occured
  @vtype      char const * restrict const
  @vio        in
  @endvar
  @var        line
  @vdesc      source line in file where error occured
  @vtype      int const
  @vio        in
  @endvar
  @var        ml
  @vdesc      meta level where error occured
  @vtype      int const
  @vio        in
  @endvar
  @var        rl
  @vdesc      refinement level where error occured
  @vtype      int const
  @vio        in
  @endvar
  @var        c
  @vdesc      first component where error occured
  @vtype      int const
  @vio        in
  @endvar
  @var        cc
  @vdesc      second component where error occured
  @vtype      int const
  @vio        in
  @endvar
  @var        message
  @vdesc      error message to display
  @vtype      char const * restrict const
  @vio        in
  @endvar
  @returntype void
@@*/
static void assert_error(char const *restrict const checkstring,
                         char const *restrict const file, int const line,
                         int const ml, int const rl, int const c, int const cc,
                         char const *restrict const message) {
  CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,
             "\n%s:%d:\n   [ml=%d rl=%d c=%d cc=%d] The following grid "
             "structure consistency check failed:\n   %s\n   %s",
             file, line, ml, rl, c, cc, message, checkstring);
  there_was_an_error = true;
}

#ifdef CARPET_OPTIMISE

// For highest efficiency, omit all self-checks
#define ASSERT_rl(check, message)
#define ASSERT_c(check, message)
#define ASSERT_cc(check, message)

#else

#define ASSERT_rl(check, message)                                              \
  do {                                                                         \
    if (not(check)) {                                                          \
      assert_error(#check, __FILE__, __LINE__, ml, rl, message);               \
    }                                                                          \
  } while (false)

#define ASSERT_c(check, message)                                               \
  do {                                                                         \
    if (not(check)) {                                                          \
      assert_error(#check, __FILE__, __LINE__, ml, rl, c, message);            \
    }                                                                          \
  } while (false)

#define ASSERT_cc(check, message)                                              \
  do {                                                                         \
    if (not(check)) {                                                          \
      assert_error(#check, __FILE__, __LINE__, ml, rl, c, cc, message);        \
    }                                                                          \
  } while (false)

#endif

/*@@
  @routine
  @date       Sun Sep 19 20:46:35 EDT 2010
  @author     Erik Schnetter, Roland Haas
  @desc
              Recompute the union of all active regions. This routine is a
              stripped down copy of dh::regrid.
  @enddesc
  @history
  @endhistory
  @var        dh
  @vdesc      data hirarchy for this group
  @vtype      const dh *
  @vio        in
  @endvar
  @var        hh
  @vdesc      grid hirarchy
  @vtype      const gh *
  @vio        in
  @endvar
  @var        ml
  @vdesc      meta level to compute active region on
  @vtype      int
  @vio        in
  @endvar
  @var        rl
  @vdesc      refinement level to compute active region on
  @vtype      int
  @vio        in
  @endvar
  @var        allactive
  @vdesc      receives the list of active boxes
  @vtype      ibset &
  @vio        out
  @endvar
  @returntype void
@@*/
void GetAllActive(const dh *dd, const gh *hh, int ml, int rl,
                  ibset &allactive) {
  DECLARE_CCTK_PARAMETERS;

  // All owned regions
  ibset allowned;

  // Domain:

  ibbox const &domain_exterior = hh->baseextent(ml, rl);
  // Variables may have size zero
  // ASSERT_rl (not domain_exterior.empty(),
  //            "The exterior of the domain must not be empty");

  i2vect const &buffer_width = dd->buffer_widths.AT(rl);
  i2vect const &boundary_width = hh->boundary_width;
  ASSERT_rl(all(all(boundary_width >= 0)),
            "The gh boundary widths must not be negative");

  ibbox const domain_active = domain_exterior.expand(-boundary_width);
  // Variables may have size zero
  // ASSERT_rl (not domain_active.empty(),
  //            "The active part of the domain must not be empty");
  ASSERT_rl(domain_active <= domain_exterior, "The active part of the domain "
                                              "must be contained in the "
                                              "exterior part of the domain");

  for (int c = 0; c < hh->components(rl); ++c) {

    // Interior:

    ibbox intr;
    intr = ibbox::poison();

    // The interior of the grid has the extent as specified by the
    // regridding thorn
    intr = hh->extent(ml, rl, c);

    // (The interior must not be empty)
    // Variables may have size zero
    // ASSERT_c (not intr.empty(),
    //           "The interior must not be empty");

    // The interior must be contained in the domain
    ASSERT_c(intr <= hh->baseextent(ml, rl),
             "The interior must be contained in the domain");

    // Outer boundary faces:

    b2vect is_outer_boundary;

    // The outer boundary faces are where the interior extends up
    // to the outer boundary of the domain.  It is not possible to
    // check whether it extends past the active part of the
    // domain, since this would be wrong when the outer boundary
    // width is zero.
    is_outer_boundary[0] = intr.lower() == domain_exterior.lower();
    is_outer_boundary[1] = intr.upper() == domain_exterior.upper();

    // Owned region:

    ibbox owned;
    owned = ibbox::poison();

    owned = intr.expand(i2vect(is_outer_boundary) * (-boundary_width));

    // (The owned region must not be empty)
    // Variables may have size zero
    // ASSERT_c (not owned.empty(),
    //           "The owned region must not be empty");

    // The owned region must be contained in the active part of
    // the domain
    ASSERT_c(
        owned <= domain_active,
        "The owned region must be contained in the active part of the domain");

    allowned |= owned;
    ASSERT_rl(
        allowned <= domain_active,
        "The owned regions must be contained in the active part of the domain");

  } // for c

  // Enlarge active part of domain
  i2vect const safedist = i2vect(0);
  ibbox const domain_enlarged = domain_active.expand(safedist);

  // All not-owned regions
  ibset const notowned = domain_enlarged - allowned;

  // All not-active points
  ibset const notactive = notowned.expand(buffer_width);

  // All active points
  allactive = allowned - notactive;

  // Fail if any of the consistency checks failed
  assert(there_was_an_error == false);
}

} // namespace CarpetIOHDF5
