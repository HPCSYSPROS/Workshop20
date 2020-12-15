#include <cctk.h>
#include <cctk_Parameters.h>

#include <cassert>
#include <cstddef>
#include <sstream>

#include <Timer.hh>

#include "bbox.hh"
#include "bboxset.hh"
#include "defs.hh"
#include "dist.hh"
#include "ggf.hh"
#include "mpi_string.hh"
#include "timestat.hh"
#include "vect.hh"

#include "dh.hh"

using namespace std;
using namespace CarpetLib;

set<dh *> dh::alldh;

// Indexing over neighbours

static int num_nbs() { return ipow(3, dim); }

static ivect ind2nb(int ind) {
  ivect nb;
  for (int d = 0; d < dim; ++d) {
    nb[d] = ind % 3 - 1;
    ind /= 3;
  }
  assert(ind == 0);
  return nb;
}

vect<vect<dh::srpvect dh::fast_dboxes::*, 2>, dim>
    dh::fast_dboxes::fast_ref_refl_sendrecv;

void dh::fast_dboxes::init_fast_ref_refl_sendrecv() {
  static bool initialised = false;
  if (initialised)
    return;
  initialised = true;
  fast_ref_refl_sendrecv[0][0] = &fast_dboxes::fast_ref_refl_sendrecv_0_0;
  fast_ref_refl_sendrecv[0][1] = &fast_dboxes::fast_ref_refl_sendrecv_0_1;
  fast_ref_refl_sendrecv[1][0] = &fast_dboxes::fast_ref_refl_sendrecv_1_0;
  fast_ref_refl_sendrecv[1][1] = &fast_dboxes::fast_ref_refl_sendrecv_1_1;
  fast_ref_refl_sendrecv[2][0] = &fast_dboxes::fast_ref_refl_sendrecv_2_0;
  fast_ref_refl_sendrecv[2][1] = &fast_dboxes::fast_ref_refl_sendrecv_2_1;
}

vect<vect<dh::srpvect dh::fast_dboxes::*, 2>, dim>
    dh::fast_dboxes::fast_ref_refl_prol_sendrecv;

void dh::fast_dboxes::init_fast_ref_refl_prol_sendrecv() {
  static bool initialised = false;
  if (initialised)
    return;
  initialised = true;
  fast_ref_refl_prol_sendrecv[0][0] =
      &fast_dboxes::fast_ref_refl_prol_sendrecv_0_0;
  fast_ref_refl_prol_sendrecv[0][1] =
      &fast_dboxes::fast_ref_refl_prol_sendrecv_0_1;
  fast_ref_refl_prol_sendrecv[1][0] =
      &fast_dboxes::fast_ref_refl_prol_sendrecv_1_0;
  fast_ref_refl_prol_sendrecv[1][1] =
      &fast_dboxes::fast_ref_refl_prol_sendrecv_1_1;
  fast_ref_refl_prol_sendrecv[2][0] =
      &fast_dboxes::fast_ref_refl_prol_sendrecv_2_0;
  fast_ref_refl_prol_sendrecv[2][1] =
      &fast_dboxes::fast_ref_refl_prol_sendrecv_2_1;
}

// Constructors
dh::dh(gh &h_, vector<i2vect> const &ghost_widths_,
       vector<i2vect> const &buffer_widths_,
       vector<i2vect> const &overlap_widths_,
       vector<int> const &prolongation_orders_space_)
    : h(h_), ghost_widths(ghost_widths_), buffer_widths(buffer_widths_),
      overlap_widths(overlap_widths_),
      prolongation_orders_space(prolongation_orders_space_) {
  fast_dboxes::init_fast_ref_refl_sendrecv();
  fast_dboxes::init_fast_ref_refl_prol_sendrecv();
  size_t const maxreflevels = h.reffacts.size();
  assert(ghost_widths.size() >= maxreflevels);
  assert(buffer_widths.size() >= maxreflevels);
  assert(overlap_widths.size() >= maxreflevels);
  assert(prolongation_orders_space.size() >= maxreflevels);
  for (size_t rl = 0; rl < maxreflevels; ++rl) {
    assert(all(all(ghost_widths.AT(rl) >= 0)));
    assert(all(all(buffer_widths.AT(rl) >= 0)));
    assert(all(all(overlap_widths.AT(rl) >= 0)));
    assert(prolongation_orders_space.AT(rl) >= 0);
  }

  alldh.insert(this);
  h.insert(this);
  CHECKPOINT;
  regrid(false);
  for (int rl = 0; rl < h.reflevels(); ++rl) {
    recompose(rl, false);
  }
  regrid_free(false);
}

// Destructors
dh::~dh() {
  CHECKPOINT;
  assert(gfs.empty());
  h.erase(this);
  alldh.erase(this);
}

// Helpers
int dh::prolongation_stencil_size(int const rl) const {
  assert(prolongation_orders_space.AT(rl) >= 0);
  return prolongation_orders_space.AT(rl) / 2;
}

// Modifiers

// Calculate this quantity on this process? It does not need to be
// calculated if it won't be used later on.
inline int dh::this_proc(int const rl, int const c) const {
  return h.processor(rl, c);
}

inline bool dh::on_this_proc(int const rl, int const c) const {
  return this_proc(rl, c) == dist::rank();
}

inline int dh::this_oldproc(int const rl, int const c) const {
  return h.old_processor(rl, c);
}

inline bool dh::on_this_oldproc(int const rl, int const c) const {
  return this_oldproc(rl, c) == dist::rank();
}

bool there_was_an_error = false;

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

void dh::regrid(bool const do_init) {
  DECLARE_CCTK_PARAMETERS;

  static Timers::Timer timer("CarpetLib::dh::regrid");
  timer.start();

  CHECKPOINT;

  static Timer total("CarpetLib::dh::regrid");
  total.start();

  light_mboxes old_light_boxes;
  swap(light_boxes, old_light_boxes);

  full_mboxes full_boxes;

  fast_boxes.clear();
  local_boxes.clear();

  light_boxes.resize(h.mglevels());
  local_boxes.resize(h.mglevels());
  level_boxes.resize(h.mglevels());
  full_boxes.resize(h.mglevels());
  fast_boxes.resize(h.mglevels());
  for (int ml = 0; ml < h.mglevels(); ++ml) {
    light_boxes.AT(ml).resize(h.reflevels());
    local_boxes.AT(ml).resize(h.reflevels());
    level_boxes.AT(ml).resize(h.reflevels());
    full_boxes.AT(ml).resize(h.reflevels());
    fast_boxes.AT(ml).resize(h.reflevels());
    for (int rl = 0; rl < h.reflevels(); ++rl) {
      light_boxes.AT(ml).AT(rl).resize(h.components(rl));
      local_boxes.AT(ml).AT(rl).resize(h.local_components(rl));
      full_boxes.AT(ml).AT(rl).resize(h.components(rl));

      light_cboxes &light_level = light_boxes.AT(ml).AT(rl);
      local_cboxes &local_level = local_boxes.AT(ml).AT(rl);
      level_dboxes &level_level = level_boxes.AT(ml).AT(rl);
      full_cboxes &full_level = full_boxes.AT(ml).AT(rl);
      fast_dboxes &fast_level = fast_boxes.AT(ml).AT(rl);

      vector<fast_dboxes> fast_level_otherprocs(dist::size());

      i2vect const &boundary_width = h.boundary_width;
      i2vect const &ghost_width = ghost_widths.AT(rl);
      i2vect const &buffer_width = buffer_widths.AT(rl);
      i2vect const &overlap_width = overlap_widths.AT(rl);

      bool const boundaries_are_trivial =
          all(all(boundary_width == 0)) && all(all(ghost_width == 0)) &&
          all(all(buffer_width == 0)) && all(all(overlap_width == 0));

      // Domain:

      static Timers::Timer timer_domain("domain");
      timer_domain.start();

      ibbox const &domain_exterior = h.baseextent(ml, rl);
      // Variables may have size zero
      // ASSERT_rl (not domain_exterior.empty(),
      //            "The exterior of the domain must not be empty");

      ASSERT_rl(all(all(boundary_width >= 0)),
                "The gh boundary widths must not be negative");

      ibbox const domain_active = domain_exterior.expand(-boundary_width);
      // Variables may have size zero
      // ASSERT_rl (not domain_active.empty(),
      //            "The active part of the domain must not be empty");
      ASSERT_rl(domain_active <= domain_exterior,
                "The active part of the domain must be contained in the "
                "exterior part of the domain");

      ibset const domain_boundary = domain_exterior - domain_active;

      timer_domain.stop();

      static Timers::Timer timer_region("region");
      timer_region.start();

      for (int c = 0; c < h.components(rl); ++c) {

        full_dboxes &box = full_level.AT(c);

        // Interior:

        ibbox &intr = box.interior;
        intr = ibbox::poison();

        // The interior of the grid has the extent as specified by the
        // regridding thorn
        intr = h.extent(ml, rl, c);

        // (The interior must not be empty)
        // Variables may have size zero
        // ASSERT_c (not intr.empty(),
        //           "The interior must not be empty");

        // The interior must be contained in the domain
        ASSERT_c(intr <= h.baseextent(ml, rl),
                 "The interior must be contained in the domain");

// All interiors must be disjunct
#ifdef CARPET_DEBUG
        for (int cc = 0; cc < c; ++cc) {
          ASSERT_cc(not intr.intersects(full_level.AT(cc).interior),
                    "All interiors must be disjunct");
        }
#endif

        // Outer boundary faces:

        b2vect &is_outer_boundary = box.is_outer_boundary;

        // The outer boundary faces are where the interior extends up
        // to the outer boundary of the domain.  It is not possible to
        // check whether it extends past the active part of the
        // domain, since this would be wrong when the outer boundary
        // width is zero.
        is_outer_boundary[0] = intr.lower() == domain_exterior.lower();
        is_outer_boundary[1] = intr.upper() == domain_exterior.upper();

        // Exterior:

        ibbox &extr = box.exterior;
        extr = ibbox::poison();

        ASSERT_c(all(all(ghost_width >= 0)),
                 "The gh ghost widths must not be negative");

        if (boundaries_are_trivial) {

          extr = intr;

        } else {

          extr = intr.expand(i2vect(not is_outer_boundary) * ghost_width);
        }

        // (The exterior must not be empty)
        // Variables may have size zero
        // ASSERT_c (not extr.empty(),
        //           "The exterior must not be empty");

        // The exterior must be contained in the domain
        ASSERT_c(extr <= domain_exterior,
                 "The exterior must be contained in the domain");

        // Cactus ghost zones (which include outer boundaries):

        ibset &ghosts = box.ghosts;
        ghosts = ibset::poison();

        if (boundaries_are_trivial) {

          ghosts = ibset();

        } else {

          ghosts = extr - intr;
        }

        // The ghosts must be contained in the domain.  Different from
        // the boundaries, the ghost can include part of the outer
        // boundary of the domain.
        ASSERT_c(ghosts <= domain_exterior,
                 "The ghost zones must be contained in the domain");

        // Communicated region:

        ibbox &comm = box.communicated;
        comm = ibbox::poison();

        if (boundaries_are_trivial) {

          comm = extr;

        } else {

          comm = extr.expand(i2vect(is_outer_boundary) * (-boundary_width));
        }

        // (The communicated region must not be empty)
        // Variables may have size zero
        // ASSERT_c (not comm.empty(),
        //           "The communicated region must not be empty");

        // The communicated region must be contained in the active
        // part of the domain
        ASSERT_c(comm <= domain_active, "The communicated region must be "
                                        "contained in the active part of the "
                                        "domain");

        // Outer boundary:

        ibset &outer_boundaries = box.outer_boundaries;
        outer_boundaries = ibset::poison();

        if (boundaries_are_trivial) {

          outer_boundaries = ibset();

        } else {

          outer_boundaries = extr - comm;
        }

        // The outer boundary must be contained in the outer boundary
        // of the domain
        ASSERT_c(outer_boundaries <= domain_boundary,
                 "The outer boundary must be contained in the outer boundary "
                 "of the domain");

        // Owned region:

        ibbox &owned = box.owned;
        owned = ibbox::poison();

        if (boundaries_are_trivial) {

          owned = intr;

        } else {

          owned = intr.expand(i2vect(is_outer_boundary) * (-boundary_width));
        }

        // (The owned region must not be empty)
        // Variables may have size zero
        // ASSERT_c (not owned.empty(),
        //           "The owned region must not be empty");

        // The owned region must be contained in the active part of
        // the domain
        ASSERT_c(owned <= domain_active, "The owned region must be contained "
                                         "in the active part of the domain");

// All owned regions must be disjunct
#ifdef CARPET_DEBUG
        for (int cc = 0; cc < c; ++cc) {
          ASSERT_cc(not owned.intersects(full_level.AT(cc).owned),
                    "All owned regions must be disjunct");
        }
#endif

        // Boundary (Carpet ghost zones, which do not include outer
        // boundaries):

        ibset &boundaries = box.boundaries;
        boundaries = ibset::poison();

        if (boundaries_are_trivial) {

          boundaries = ibset();

        } else {

          boundaries = comm - owned;
        }

        // The boundary must be contained in the active part of the
        // domain.  This prevents that a region is too close to the
        // outer boundary, so that it has ghost zones overlapping with
        // the outer boundary.
        ASSERT_c(
            boundaries <= domain_active,
            "The boundary must be contained in the active part of the domain");

      } // for c

      timer_region.stop();

      // Conjunction of all buffer zones:

      static Timers::Timer timer_buffers("buffers");
      timer_buffers.start();

      // Enlarge active part of domain
      i2vect const safedist = i2vect(0);
      ibbox const domain_enlarged = domain_active.expand(safedist);

      // All owned regions
      ibset const allowned(full_level, &full_dboxes::owned);
      ASSERT_rl(allowned <= domain_active, "The owned regions must be "
                                           "contained in the active part of "
                                           "the domain");

      // All not-owned regions
      ibset const notowned = domain_enlarged - allowned;

      // All not-active points
      ibset const notactive = notowned.expand(buffer_width + overlap_width);

      // All not-active points, in stages
      int const num_substeps = any(any(ghost_width == 0))
                                   ? 0
                                   : minval(minval(buffer_width / ghost_width));
      if (not all(all(buffer_width == num_substeps * ghost_width))) {
        ostringstream buf;
        buf << "The buffer width " << buffer_width
            << " is not a multiple of the ghost width " << ghost_width
            << " on level " << rl;
        CCTK_WARN(CCTK_WARN_COMPLAIN, buf.str().c_str());
      }

#ifdef CARPET_DEBUG
      vector<ibset> notactive_stepped(num_substeps + 1);
      notactive_stepped.AT(0) = notowned;
      for (int substep = 1; substep <= num_substeps; ++substep) {
        notactive_stepped.AT(substep) =
            notactive_stepped.AT(substep - 1).expand(ghost_width);
      }
      ibset const notactive_overlaps =
          notactive_stepped.AT(num_substeps).expand(overlap_width);
      if (all(all(buffer_width == num_substeps * ghost_width))) {
        ASSERT_rl(notactive_overlaps == notactive,
                  "The stepped not-owned region including overlaps must be "
                  "equal to the not-active region");
      }
#endif

      // All buffer zones
      // ibset const allbuffers = allowned & notowned.expand (buffer_width);
      ibset &allbuffers = level_level.buffers;
      allbuffers = allowned & notowned.expand(buffer_width);

      // All overlap zones
      ibset const alloverlaps = (allowned & notactive) - allbuffers;

      // All active points
      ibset &allactive = level_level.active;
      allactive = allowned - notactive;

#ifdef CARPET_DEBUG
      // All stepped buffer zones
      vector<ibset> allbuffers_stepped(num_substeps);
      ibset allbuffers_stepped_combined;
      for (int substep = 0; substep < num_substeps; ++substep) {
        allbuffers_stepped.AT(substep) =
            allowned &
            (notactive_stepped.AT(substep + 1) - notactive_stepped.AT(substep));
        allbuffers_stepped_combined += allbuffers_stepped.AT(substep);
      }
      if (all(all(buffer_width == num_substeps * ghost_width))) {
        ASSERT_rl(allbuffers_stepped_combined == allbuffers,
                  "The stepped buffer zones must be equal to the buffer zones");
      }
#endif

      // Overlap zones and buffer zones must be in the active part of
      // the domain
      ASSERT_rl(allactive <= domain_active,
                "The active region must be in the active part of the domain");
      ASSERT_rl(alloverlaps <= domain_active,
                "The overlap zones must be in the active part of the domain");
      ASSERT_rl(allbuffers <= domain_active,
                "The buffer zones must be in the active part of the domain");
      ASSERT_rl((allactive & alloverlaps).empty(),
                "The active points and the overlap zones cannot overlap");
      ASSERT_rl((allactive & allbuffers).empty(),
                "The active points and the buffer zones cannot overlap");
      ASSERT_rl((alloverlaps & allbuffers).empty(),
                "The overlap zones and the buffer zones cannot overlap");
      ASSERT_rl(allactive + alloverlaps + allbuffers == allowned,
                "The active points, the overlap points, and buffer points "
                "together must be exactly the owned region");

      for (int c = 0; c < h.components(rl); ++c) {
        full_dboxes &box = full_level.AT(c);

        if (boundaries_are_trivial) {

          // Buffer zones:
          box.buffers = ibset();

          // Overlap zones:
          box.overlaps = ibset();

          // Active region:
          box.active = box.owned;

        } else {

          // Buffer zones:
          box.buffers = box.owned & allbuffers;

          // Overlap zones:
          box.overlaps = box.owned & alloverlaps;

          // Active region:
          box.active = box.owned & allactive;
        }

        ASSERT_c(box.active == box.owned - (box.buffers + box.overlaps),
                 "The active region must equal the owned region minus the "
                 "(buffer zones plus overlap zones)");
      } // for c

      for (int lc = 0; lc < h.local_components(rl); ++lc) {
        int const c = h.get_component(rl, lc);
        local_dboxes &local_box = local_level.AT(lc);
        full_dboxes const &box = full_level.AT(c);

        local_box.buffers = box.buffers;
        local_box.overlaps = box.overlaps;

#if 0
        // Stepped buffer zones:
        local_box.buffers_stepped.resize (num_substeps);
        for (int substep = 0; substep < num_substeps; ++ substep) {
          local_box.buffers_stepped.AT(substep) =
            box.owned & allbuffers_stepped.AT(substep);
        }
#endif

        local_box.active = box.active;
      } // for lc

      // The conjunction of all buffer zones must equal allbuffers
      ibset const allbuffers1(full_level, &full_dboxes::buffers);
      ASSERT_rl(allbuffers1 == allbuffers, "Buffer zone consistency check");

      timer_buffers.stop();

      // Test constituency relations:

      static Timers::Timer timer_test("test");
      timer_test.start();

      for (int c = 0; c < h.components(rl); ++c) {
        full_dboxes const &box = full_level.AT(c);

        ASSERT_c((box.active & box.buffers).empty(), "Consistency check");
        ASSERT_c((box.active & box.overlaps).empty(), "Consistency check");
        ASSERT_c((box.overlaps & box.buffers).empty(), "Consistency check");
        ASSERT_c((box.active | box.overlaps | box.buffers) == box.owned,
                 "Consistency check");

        ASSERT_c((box.owned & box.boundaries).empty(), "Consistency check");
        ASSERT_c((box.owned | box.boundaries) == box.communicated,
                 "Consistency check");

        ASSERT_c((box.communicated & box.outer_boundaries).empty(),
                 "Consistency check");
        ASSERT_c((box.communicated | box.outer_boundaries) == box.exterior,
                 "Consistency check");

        ASSERT_c(box.boundaries <= box.ghosts, "Consistency check");

        ASSERT_c((box.interior & box.ghosts).empty(), "Consistency check");
        ASSERT_c((box.interior | box.ghosts) == box.exterior,
                 "Consistency check");

      } // for c

      timer_test.stop();

      // Communication schedule:

      static Timers::Timer timer_comm("comm");
      timer_comm.start();

      static Timers::Timer timer_comm_mgrest("mgrest");
      static Timers::Timer timer_comm_mgprol("mgprol");
      static Timers::Timer timer_comm_refprol("refprol");
      static Timers::Timer timer_comm_sync("sync");
      static Timers::Timer timer_comm_refbndprol("refbndprol");
      timer_comm_mgrest.instantiate();
      timer_comm_mgprol.instantiate();
      timer_comm_refprol.instantiate();
      timer_comm_sync.instantiate();
      timer_comm_refbndprol.instantiate();

      for (int lc = 0; lc < h.local_components(rl); ++lc) {
        int const c = h.get_component(rl, lc);

        full_dboxes &box = full_level.AT(c);

        // Multigrid restriction:

        timer_comm_mgrest.start();

        if (ml > 0) {
          int const oml = ml - 1;

          // Multigrid restriction must fill all active points

          full_dboxes const &obox = full_boxes.AT(oml).AT(rl).AT(c);

          ibset needrecv = box.active;

          ibset const contracted_oactive(
              obox.active.contracted_for(box.interior));
          ibset const ovlp = needrecv & contracted_oactive;

          for (ibset::const_iterator ri = ovlp.begin(); ri != ovlp.end();
               ++ri) {
            ibbox const &recv = *ri;
            ibbox const send = recv.expanded_for(obox.interior);
            ASSERT_c(send <= obox.exterior, "Multigrid restriction: Send "
                                            "region must be contained in "
                                            "exterior");
            fast_level.fast_mg_rest_sendrecv.push_back(
                sendrecv_pseudoregion_t(send, c, recv, c));
          }

          needrecv -= ovlp;

          // All points must have been received
          ASSERT_c(needrecv.empty(),
                   "Multigrid restriction: All points must have been received");

        } // if ml > 0

        timer_comm_mgrest.stop();

        // Multigrid prolongation:

        timer_comm_mgprol.start();

        if (ml > 0) {
          int const oml = ml - 1;

          // Multigrid prolongation must fill all active points
          // (this could probably be relaxed)

          full_dboxes const &obox = full_boxes.AT(oml).AT(rl).AT(c);

          ibset oneedrecv = obox.active;

          i2vect const stencil_size = i2vect(prolongation_stencil_size(rl));

          ibset const expanded_active(box.active.expanded_for(obox.interior));
          ibset const ovlp = oneedrecv & expanded_active;

          for (ibset::const_iterator ri = ovlp.begin(); ri != ovlp.end();
               ++ri) {
            ibbox const &recv = *ri;
            ibbox const send =
                recv.expanded_for(box.interior).expand(stencil_size);
            ASSERT_c(send <= box.exterior, "Multigrid prolongation: Send "
                                           "region must be contained in "
                                           "exterior");
            fast_level.fast_mg_prol_sendrecv.push_back(
                sendrecv_pseudoregion_t(send, c, recv, c));
          }

          oneedrecv -= ovlp;

          // All points must have been received
          ASSERT_c(
              oneedrecv.empty(),
              "Multigrid prolongation: All points must have been received");

        } // if ml > 0

        timer_comm_mgprol.stop();

        // Refinement prolongation:

        timer_comm_refprol.start();

        if (rl > 0) {
          int const orl = rl - 1;

          // Refinement prolongation must fill all active points

          ibset needrecv = box.active + box.overlaps;

          i2vect const stencil_size = i2vect(prolongation_stencil_size(rl));

          ASSERT_c(
              all(h.reffacts.at(rl) % h.reffacts.at(orl) == 0),
              "Refinement factors must be integer multiples of each other");
          i2vect const reffact = i2vect(h.reffacts.at(rl) / h.reffacts.at(orl));

          for (int cc = 0; cc < h.components(orl); ++cc) {
            if (needrecv.empty())
              break;

            full_dboxes const &obox = full_boxes.AT(ml).AT(orl).AT(cc);

#if 0
            // This does not work for cell centering; if the domain is
            // 1 cell wide, the contraction disappears it
            ibset const expanded_oactive
              (obox.active.contracted_for (box.interior).expand (reffact));
#else
            ibset const expanded_oactive(
                obox.active.expanded_for(box.interior)
                    .expand(h.refcent == vertex_centered ? reffact
                                                         : reffact - 1));
#endif
            ibset const ovlp = needrecv & expanded_oactive;

            for (ibset::const_iterator ri = ovlp.begin(); ri != ovlp.end();
                 ++ri) {
              ibbox const &recv = *ri;
              ibbox const send =
                  recv.expanded_for(obox.interior).expand(stencil_size);
              ASSERT_c(send <= obox.exterior, "Refinement prolongation: Send "
                                              "region must be contained in "
                                              "exterior");
              fast_level.fast_ref_prol_sendrecv.push_back(
                  sendrecv_pseudoregion_t(send, cc, recv, c));
              if (not on_this_proc(orl, cc)) {
                fast_dboxes &fast_level_otherproc =
                    fast_level_otherprocs.AT(this_proc(orl, cc));
                fast_level_otherproc.fast_ref_prol_sendrecv.push_back(
                    sendrecv_pseudoregion_t(send, cc, recv, c));
              }
            }

            needrecv -= ovlp;

          } // for cc

          // All points must have been received
          if (not needrecv.empty()) {
            cerr << "box.active=" << box.active << "\n"
                 << "needrecv=" << needrecv << "\n";
          }
          ASSERT_c(
              needrecv.empty(),
              "Refinement prolongation: All points must have been received");

        } // if rl > 0

        timer_comm_refprol.stop();

        // Synchronisation:

        timer_comm_sync.start();

        {

// Synchronisation should fill as many boundary points as
// possible

#if 0
          // Outer boundaries are not synchronised, since they cannot
          // be filled by boundary prolongation either, and therefore
          // the user code must set them anyway.
          ibset needrecv = box.boundaries;
#else
          // Outer boundaries are synchronised for backward
          // compatibility.
          ibset needrecv = box.ghosts;
#endif

          ibset &sync = box.sync;

          for (int cc = 0; cc < h.components(rl); ++cc) {
            if (needrecv.empty())
              break;

            full_dboxes const &obox = full_level.AT(cc);

#if 0
            ibset const ovlp = needrecv & obox.owned;
#else
            ibset const ovlp = needrecv & obox.interior;
#endif

            if (cc == c) {
              ASSERT_cc(ovlp.empty(),
                        "A region may not synchronise from itself");
            }

            for (ibset::const_iterator ri = ovlp.begin(); ri != ovlp.end();
                 ++ri) {
              ibbox const &recv = *ri;
              ibbox const &send = recv;
              fast_level.fast_sync_sendrecv.push_back(
                  sendrecv_pseudoregion_t(send, cc, recv, c));
              if (not on_this_proc(rl, cc)) {
                fast_dboxes &fast_level_otherproc =
                    fast_level_otherprocs.AT(this_proc(rl, cc));
                fast_level_otherproc.fast_sync_sendrecv.push_back(
                    sendrecv_pseudoregion_t(send, cc, recv, c));
              }
            }

            needrecv -= ovlp;
            sync += ovlp;

          } // for cc
        }

        timer_comm_sync.stop();

        // Boundary prolongation:

        timer_comm_refbndprol.start();

        if (rl > 0) {
          int const orl = rl - 1;

#if 0
          // Outer boundaries are not synchronised, since they cannot
          // be filled by boundary prolongation either, and therefore
          // the user code must set them anyway.
          ibset needrecv = box.boundaries;
#else
          // Outer boundaries are synchronised for backward
          // compatibility.
          ibset needrecv = box.ghosts;
#endif

          // Points which are synchronised need not be boundary
          // prolongated
          needrecv -= box.sync;

          // Outer boundary points cannot be boundary prolongated
          needrecv &= box.communicated;

          // Prolongation must fill what cannot be synchronised, and
          // also all buffer zones
          needrecv += box.buffers;

          ibset &bndref = box.bndref;

          i2vect const stencil_size = i2vect(prolongation_stencil_size(rl));

          ASSERT_c(
              all(h.reffacts.at(rl) % h.reffacts.at(orl) == 0),
              "Refinement factors must be integer multiples of each other");
          i2vect const reffact = i2vect(h.reffacts.at(rl) / h.reffacts.at(orl));

          for (int cc = 0; cc < h.components(orl); ++cc) {
            if (needrecv.empty())
              break;

            full_dboxes const &obox = full_boxes.AT(ml).AT(orl).AT(cc);

            // See "refinement prolongation"
            ibset const expanded_oactive(
                obox.active.expanded_for(box.interior)
                    .expand(h.refcent == vertex_centered ? reffact
                                                         : reffact - 1));
            ibset const ovlp = needrecv & expanded_oactive;

            for (ibset::const_iterator ri = ovlp.begin(); ri != ovlp.end();
                 ++ri) {
              ibbox const &recv = *ri;
              ibbox const send =
                  recv.expanded_for(obox.interior).expand(stencil_size);
              ASSERT_c(send <= obox.exterior, "Boundary prolongation: Send "
                                              "region must be contained in "
                                              "exterior");
              fast_level.fast_ref_bnd_prol_sendrecv.push_back(
                  sendrecv_pseudoregion_t(send, cc, recv, c));
              if (not on_this_proc(orl, cc)) {
                fast_dboxes &fast_level_otherproc =
                    fast_level_otherprocs.AT(this_proc(orl, cc));
                fast_level_otherproc.fast_ref_bnd_prol_sendrecv.push_back(
                    sendrecv_pseudoregion_t(send, cc, recv, c));
              }
            }

            needrecv -= ovlp;
            bndref += ovlp;

          } // for cc

          // All points must now have been received, either through
          // synchronisation or through boundary prolongation
          ASSERT_c(needrecv.empty(), "Synchronisation and boundary "
                                     "prolongation: All points must have been "
                                     "received");

        } // if rl > 0

        timer_comm_refbndprol.stop();

      } // for lc

      // Refinement restriction:

      static Timers::Timer timer_comm_refrest("refrest");
      timer_comm_refrest.start();

      if (rl > 0) {
        int const orl = rl - 1;
        full_cboxes const &full_olevel = full_boxes.AT(ml).AT(orl);
        fast_dboxes &fast_olevel = fast_boxes.AT(ml).AT(orl);
        ibbox const &odomext = h.baseextent(ml, orl);

        // Refinement restriction may fill all coarse interior points,
        // and must use all fine active points

        ibset allrestricted;
        switch (h.refcent) {
        case vertex_centered:
          allrestricted = allactive.contracted_for(odomext);
          break;
        case cell_centered: {
          ibset const &source = allactive;
          ibbox const &target = odomext;
          ibbox const all_source = allactive.container().expand(10, 10);
          ibbox const all_target = all_source.contracted_for(target);
          ibset const tmp0 = source;
          ibset const tmp1 = tmp0.expanded_for(target);
          ibset const tmp2 = all_target - tmp1;
          ibset const tmp3 = tmp2.expand(1, 1);
          ibset const tmp4 = all_target - tmp3;
          allrestricted = tmp4;
          break;
        }
        default:
          assert(0);
        }

        for (int olc = 0; olc < h.local_components(orl); ++olc) {
          int const oc = h.get_component(orl, olc);
          full_dboxes const &obox = full_olevel.AT(oc);

          ibset needrecv;
          if (not use_higher_order_restriction and
              not support_staggered_operators) {
            // Ghost zones: We can restrict into ghost zones if not
            // using higher order restriction, which is probably much
            // cheaper than performing a sync after restriction
            needrecv = allrestricted & obox.exterior;
          } else {
            // We do not restrict into boundaries or ghost zones, but
            // we require them. This means that we need to synchronize
            // and apply boundary conditions after restricting on each
            // level.
            needrecv = allrestricted & obox.interior;
          }
          // Cannot restrict into buffer zones
          assert((allrestricted & obox.buffers).empty());

          for (int c = 0; c < h.components(rl); ++c) {
            if (needrecv.empty())
              break;

            full_dboxes const &box = full_level.AT(c);

            // If we make a mistake expanding the domain of dependence here, it
            // _should_ be caught be by the expand()ed is_contained_in(srcbbox)
            // test in the actual operator.
            int shrink_by;
            if (use_higher_order_restriction and h.refcent == cell_centered)
              shrink_by = restriction_order_space / 2;
            else if (support_staggered_operators)
              shrink_by = 2; // RH: not sure why 2 is required
            else
              shrink_by = 0;
            ibbox const contracted_exterior =
                box.exterior.expand(ivect(-shrink_by)).contracted_for(odomext);
            ibset const ovlp = needrecv & contracted_exterior;

            for (ibset::const_iterator ri = ovlp.begin(); ri != ovlp.end();
                 ++ri) {
              ibbox const &recv = *ri;
              ibbox const send =
                  recv.expanded_for(box.exterior).expand(ivect(shrink_by));
              ASSERT_c(send <= box.exterior, "Refinement restriction: Send "
                                             "region must be contained in "
                                             "exterior");

              sendrecv_pseudoregion_t const preg(send, c, recv, oc);
              fast_olevel.fast_ref_rest_sendrecv.push_back(preg);
              if (not on_this_proc(rl, c)) {
                fast_dboxes &fast_level_otherproc =
                    fast_level_otherprocs.AT(this_proc(rl, c));
                fast_level_otherproc.fast_ref_rest_sendrecv.push_back(preg);
              }
            }

            needrecv -= ovlp;

          } // for c

          // All points must have been received
          ASSERT_rl(
              needrecv.empty(),
              "Refinement restriction: All points must have been received");

        } // for olc

      } // if rl > 0

      timer_comm_refrest.stop();

      // Refinement refluxing:

      static Timers::Timer timer_comm_reflux("reflux");
      timer_comm_reflux.start();

      // If there is no coarser level, do nothing
      if (rl > 0) {
        int const orl = rl - 1;
        local_cboxes &local_olevel = local_boxes.AT(ml).AT(orl);
        full_cboxes const &full_olevel = full_boxes.AT(ml).AT(orl);
        fast_dboxes &fast_olevel = fast_boxes.AT(ml).AT(orl);

        // This works only with cell centred grids
        if (h.refcent == cell_centered) {

          // Fine grids
          ibset const &fine_level = allactive;

          assert(all(h.reffacts.AT(rl) % h.reffacts.AT(orl) == 0));

          vect<vect<ibset, 2>, dim> all_fine_boundary;

          for (int dir = 0; dir < dim; ++dir) {
            // Unit vector
            ivect const idir = ivect::dir(dir);

            // Left and right faces
            ibset const fine_level_minus = fine_level.shift(-idir, 2);
            ibset const fine_level_plus = fine_level.shift(+idir, 2);

            // Fine boundaries
            all_fine_boundary[dir][0] = fine_level_minus - fine_level_plus;
            all_fine_boundary[dir][1] = fine_level_plus - fine_level_minus;
          } // for dir

          // Coarse grids
          ibbox const &coarse_domain_exterior = h.baseextent(ml, orl);
          ibbox const &coarse_ext = coarse_domain_exterior;
          ibbox const &fine_domain_exterior = h.baseextent(ml, rl);
          ibbox const &fine_ext = fine_domain_exterior;

          // Coarse equivalent of fine boundary
          vect<vect<ibset, 2>, dim> all_coarse_boundary;
          for (int dir = 0; dir < 3; ++dir) {
            // Unit vector
            ivect const idir = ivect::dir(dir);
            ibbox const coarse_faces = coarse_ext.shift(-idir, 2);
            ibbox const fine_faces = fine_ext.shift(-idir, 2);
            for (int face = 0; face < 2; ++face) {
              if (verbose) {
                cout << "REFREF rl=" << rl << " dir=" << dir << " face=" << face
                     << "\n"
                     << "   all_fine_boundary=" << all_fine_boundary[dir][face]
                     << "\n"
                     << "   coarse_ext=" << coarse_ext << "\n"
                     << "   coarse_faces=" << coarse_faces << "\n";
              }
              all_coarse_boundary[dir][face] =
                  all_fine_boundary[dir][face].contracted_for(coarse_faces);
              if (verbose) {
                cout << "   all_coarse_boundary="
                     << all_coarse_boundary[dir][face] << "\n";
              }
              ASSERT_rl(all_coarse_boundary[dir][face].expanded_for(
                            fine_faces) == all_fine_boundary[dir][face],
                        "Refluxing: Coarse grid boundaries must be consistent "
                        "with fine grid boundaries");
            } // for face
          }   // for dir

          // Split onto components
          for (int lc = 0; lc < h.local_components(rl); ++lc) {
            int const c = h.get_component(rl, lc);
            full_dboxes &box = full_level.AT(c);
            local_dboxes &local_box = local_level.AT(lc);

            for (int dir = 0; dir < dim; ++dir) {
              // Unit vector
              ivect const idir = ivect::dir(dir);
              for (int face = 0; face < 2; ++face) {

                // This is not really used (only for debugging)
                local_box.fine_boundary[dir][face] =
                    box.exterior.shift(-idir, 2) & all_fine_boundary[dir][face];
                if (verbose) {
                  cout << "REFREF rl=" << rl << " lc=" << lc << " dir=" << dir
                       << " face=" << face << "\n"
                       << "   local.fine_boundary="
                       << local_box.fine_boundary[dir][face] << "\n";
                }

#if 0 // OFFSETS,SIZE
                ibset const& boxes =
                  local_box.fine_boundary[dir][face];
                vector<int>& offsets =
                  local_box.fine_boundary_offsets[dir][face];
                assert(offsets.empty());
                offsets.reserve(boxes.setsize());
                int offset = 0;
                for (ibset::const_iterator
                       bi=boxes.begin(); bi!=boxes.end(); ++bi)
                {
                  offsets.push_back(offset);
                  offset += (*bi).size();
                }
                local_box.fine_boundary_size[dir][face] = offset;
#endif

              } // for face
            }   // for dir

          } // for lc

#ifdef CARPET_DEBUG
          for (int dir = 0; dir < dim; ++dir) {
            ivect const idir = ivect::dir(dir); // Unit vector
            for (int face = 0; face < 2; ++face) {
              ibset all_fine_boundary_combined;
              for (int c = 0; c < h.components(rl); ++c) {
                full_dboxes const &box = full_level.AT(c);
                all_fine_boundary_combined |=
                    box.exterior.shift(-idir, 2) & all_fine_boundary[dir][face];
              } // for c
              ASSERT_rl(all_fine_boundary_combined ==
                            all_fine_boundary[dir][face],
                        "Refluxing: Union of all fine grid boundaries must be "
                        "the total fine grid boundary ");
            }
          }
#endif

          for (int lc = 0; lc < h.local_components(orl); ++lc) {
            int const c = h.get_component(orl, lc);
            full_dboxes const &obox = full_olevel.AT(c);
            local_dboxes &local_obox = local_olevel.AT(lc);

            for (int dir = 0; dir < dim; ++dir) {
              // Unit vector
              ivect const idir = ivect::dir(dir);
              for (int face = 0; face < 2; ++face) {

                // This is used for post-processing the fluxes
                // (calculating the difference between coarse and fine
                // fluxes, adjusting the state)
                local_obox.coarse_boundary[dir][face] =
                    obox.owned.shift(-idir, 2) & all_coarse_boundary[dir][face];
                // Cannot reflux into buffer zones
                if (verbose) {
                  cout << "REFREF orl=" << orl << " lc=" << lc << " dir=" << dir
                       << " face=" << face << "\n"
                       << "   local.coarse_boundary="
                       << local_obox.coarse_boundary[dir][face] << "\n";
                }
                assert((obox.buffers.shift(-idir, 2) &
                        all_coarse_boundary[dir][face])
                           .empty());

#if 0 // OFFSETS,SIZE
                ibset const& boxes =
                  local_obox.coarse_boundary[dir][face];
                vector<int>& offsets =
                  local_obox.coarse_boundary_offsets[dir][face];
                assert(offsets.empty());
                offsets.reserve(boxes.setsize());
                int offset = 0;
                for (ibset::const_iterator
                       bi=boxes.begin(); bi!=boxes.end(); ++bi)
                {
                  offsets.push_back(offset);
                  offset += (*bi).size();
                }
                local_obox.coarse_boundary_size[dir][face] = offset;
#endif

              } // for face
            }   // for dir

          } // for lc

#ifdef CARPET_DEBUG
          for (int dir = 0; dir < dim; ++dir) {
            ivect const idir = ivect::dir(dir); // Unit vector
            for (int face = 0; face < 2; ++face) {
              ibset all_coarse_boundary_combined;
              for (int c = 0; c < h.components(orl); ++c) {
                full_dboxes const &obox = full_olevel.AT(c);
                all_coarse_boundary_combined +=
                    obox.owned.shift(face ? +idir : -idir, 2) &
                    all_coarse_boundary[dir][face];
              } // for c
              ASSERT_rl(all_coarse_boundary_combined ==
                            all_coarse_boundary[dir][face],
                        "Refluxing: Union of all coarse grid boundaries must "
                        "be the total coarse grid boundary ");
            }
          }
#endif

          // Set up communication schedule
          for (int olc = 0; olc < h.local_components(orl); ++olc) {
            int const oc = h.get_component(orl, olc);
            local_dboxes &local_obox = local_olevel.AT(olc);
            if (verbose) {
              cout << "REF ref_refl ml=" << ml << " rl=" << rl << " olc=" << olc
                   << " oc=" << oc << "\n";
            }

            for (int dir = 0; dir < dim; ++dir) {
              // Unit vector
              ivect const idir = ivect::dir(dir);
              ibbox const coarse_faces = coarse_ext.shift(-idir, 2);
              for (int face = 0; face < 2; ++face) {

                srpvect fast_dboxes::*const fast_ref_refl_sendrecv =
                    fast_dboxes::fast_ref_refl_sendrecv[dir][face];
                srpvect fast_dboxes::*const fast_ref_refl_prol_sendrecv =
                    fast_dboxes::fast_ref_refl_prol_sendrecv[dir][face];

                // Refluxing must fill all coarse refluxing boundary
                // points, and may use all fine points

                ibset needrecv(local_obox.coarse_boundary[dir][face]);

                for (int c = 0; c < h.components(rl); ++c) {
                  if (needrecv.empty())
                    break;

                  full_dboxes const &box = full_level.AT(c);
                  if (verbose) {
                    cout << "REF ref_refl ml=" << ml << " rl=" << rl
                         << " olc=" << olc << " oc=" << oc << " c=" << c
                         << " dir=" << dir << " face=" << face << "\n";
                  }

                  ibbox const contracted_exterior =
                      box.exterior.shift(-idir, 2).contracted_for(coarse_faces);
                  if (verbose) {
                    cout << "   exterior=" << box.exterior.shift(-idir, 2)
                         << "\n"
                         << "   contracted=" << contracted_exterior << "\n";
                  }
                  ibset const ovlp = needrecv & contracted_exterior;
                  if (verbose) {
                    cout << "   ovlp=" << ovlp << "\n";
                  }

                  for (ibset::const_iterator ri = ovlp.begin();
                       ri != ovlp.end(); ++ri) {
                    ibbox const &recv = *ri;
                    ibbox const send =
                        recv.expanded_for(box.exterior.shift(-idir, 2));
                    ASSERT_c(send <= box.exterior.shift(-idir, 2),
                             "Refinement restriction: Send region must be "
                             "contained in exterior");
                    ibbox const shifted_recv = recv.shift(idir, 2);
                    ibbox const shifted_send = send.shift(idir, 2);
                    sendrecv_pseudoregion_t const preg(shifted_send, c,
                                                       shifted_recv, oc);
                    if (verbose) {
                      cout << "REF ref_refl ml=" << ml << " rl=" << rl
                           << " olc=" << olc << " c=" << c << " oc=" << oc
                           << " dir=" << dir << " face=" << face << "\n"
                           << "   preg=" << preg << "\n";
                    }
                    (fast_olevel.*fast_ref_refl_sendrecv).push_back(preg);
                    if (not on_this_proc(rl, c)) {
                      fast_dboxes &fast_level_otherproc =
                          fast_level_otherprocs.AT(this_proc(rl, c));
                      (fast_level_otherproc.*fast_ref_refl_sendrecv)
                          .push_back(preg);
                    }

                    // reflux-prolongation
                    sendrecv_pseudoregion_t const preg_prol(shifted_recv, oc,
                                                            shifted_send, c);
                    if (verbose) {
                      cout << "REF ref_refL_prol ml=" << ml << " rl=" << rl
                           << " olc=" << olc << " c=" << c << " oc=" << oc
                           << " dir=" << dir << " face=" << face << "\n"
                           << "   preg=" << preg_prol << "\n";
                    }
                    (fast_level.*fast_ref_refl_prol_sendrecv)
                        .push_back(preg_prol);
                    if (not on_this_proc(orl, oc)) {
                      fast_dboxes &fast_level_otherproc =
                          fast_level_otherprocs.AT(this_proc(orl, oc));
                      (fast_level_otherproc.*fast_ref_refl_prol_sendrecv)
                          .push_back(preg_prol);
                    }
                  }

                  needrecv -= ovlp;

                } // for c

                // All points must have been received
                if (not needrecv.empty()) {
                  cout << "needrecv=" << needrecv << endl;
                }
                ASSERT_rl(
                    needrecv.empty(),
                    "Refinement refluxing: All points must have been received");

              } // for face
            }   // for dir

          } // for olc

        } // if cell centered

      } // if rl > 0

      timer_comm_reflux.stop();

      timer_comm.stop();

      // Reduction mask:

      static Timers::Timer timer_mask("mask");
      timer_mask.start();

      for (int lc = 0; lc < h.local_components(rl); ++lc) {
        local_dboxes &local_box = local_level.AT(lc);
        local_box.prolongation_boundary.resize(num_nbs());
        local_box.restriction_boundary.resize(num_nbs());
      }

      if (rl > 0) {
        int const orl = rl - 1;

        // Set the weight in the interior of the notactive and the
        // allactive regions to zero, and set the weight on the
        // boundary of the notactive and allactive regions to 1/2.
        //
        // For the prolongation region, the "boundary" is the first
        // layer outside of the region, and the "region" consists of
        // the inactive points. This corresponds to the outermost
        // layer of active grid points.
        //
        // For the restricted region, the "boundary" is the outermost
        // layer of grid points if this layer is aligned with the next
        // coarser (i.e. the current) grid; otherwise, the boundary is
        // empty. The "region" consists of the active points.

        // Note: We calculate the prolongation information for the
        // current level, and the restriction information for the next
        // coarser level. We do it this way because both calculations
        // start from the current level's set of active grid points.

        full_cboxes const &full_olevel = full_boxes.AT(ml).AT(orl);
        // Local boxes are not communicated or post-processed, and
        // thus can be modified even for coarser levels
        local_cboxes &local_olevel = local_boxes.AT(ml).AT(orl);

        if (verbose) {
          ostringstream buf;
          buf << "Setting prolongation region " << notactive << " on level "
              << rl;
          CCTK_INFO(buf.str().c_str());
        }
        if (verbose) {
          ostringstream buf;
          buf << "Setting restriction region " << allactive << " on level "
              << orl;
          CCTK_INFO(buf.str().c_str());
        }

        for (int neighbour = 0; neighbour < num_nbs(); ++neighbour) {
          ivect const shift = ind2nb(neighbour);

          // In this loop, shift [1,1,1] denotes a convex corner of
          // the region which should be masked out, i.e. a region
          // where only a small part (1/8) of the region should be
          // masked out. Concave edges are treated implicitly
          // (sequentially), i.e. larger chunks are cut out multiple
          // times: a concave xy edge counts as both an x face and a y
          // face.

          // Prolongation boundary of this level: start with inactive
          // (prolongated from the next coarser level) points
          ibset boxes = notactive;
          // Restriction boundary of the next coarser level: start
          // with this level's active (restricted to the next coarser)
          // points
          ibset fboxes = allactive;

          switch (h.refcent) {
          case vertex_centered:
            fboxes |= domain_boundary;
            for (int d = 0; d < dim; ++d) {
              ivect const dir = ivect::dir(d);
              fboxes = fboxes.shift(-dir) & fboxes & fboxes.shift(+dir);
            }
            fboxes -= domain_boundary;
            for (int d = 0; d < dim; ++d) {
              // Calculate the boundary in direction d
              ivect const dir = ivect::dir(d);
              switch (shift[d]) {
              case -1: {
                // left boundary
                boxes = boxes.shift(-dir) - boxes;
                fboxes = fboxes.shift(-dir) - fboxes;
                break;
              }
              case 0: {
                // interior
                // do nothing
                break;
              }
              case +1: {
                // right boundary
                boxes = boxes.shift(+dir) - boxes;
                fboxes = fboxes.shift(+dir) - fboxes;
                break;
              }
              default:
                assert(0);
              }
            }
            break;
          case cell_centered:
            // For cell-centred grids we assume that all cell
            // boundaries are aligned. This may not actually be the
            // case.
            // TODO: assert this.
            if (all(shift == 0)) {
              // do nothing
            } else {
              boxes = ibset();
              fboxes = ibset();
            }
            break;
          default:
            assert(0);
          } // switch centering

          ibbox const &odomext = h.baseextent(ml, orl);
          ibset const cfboxes = fboxes.contracted_for(odomext);

          if (verbose) {
            ostringstream buf;
            buf << "Setting boundary " << shift << ": prolongation region "
                << boxes;
            CCTK_INFO(buf.str().c_str());
          }
          if (verbose) {
            ostringstream buf;
            buf << "Setting boundary " << shift << ": restriction region "
                << fboxes;
            CCTK_INFO(buf.str().c_str());
          }

          // For the current level, the light boxes do not (yet) exist
          // here, so we use the full boxes
          for (int lc = 0; lc < h.local_components(rl); ++lc) {
            int const c = h.get_component(rl, lc);
            full_dboxes const &full_box = full_level.AT(c);
            local_dboxes &local_box = local_level.AT(lc);
            local_box.prolongation_boundary.AT(neighbour) =
                boxes & full_box.exterior;
          }
          // For the next coarser level, the full boxes do not exist
          // (any more), so we use the light boxes
          for (int olc = 0; olc < h.local_components(orl); ++olc) {
            int const oc = h.get_component(orl, olc);
            full_dboxes const &full_obox = full_olevel.AT(oc);
            local_dboxes &local_obox = local_olevel.AT(olc);
            local_obox.restriction_boundary.AT(neighbour) =
                cfboxes & full_obox.exterior;
          }

        } // for neighbour

      } // if not coarsest level

      timer_mask.stop();

      // Mask for unused points on coarser level (which do not
      // influence the future evolution provided regridding is done at
      // the right times):
      static Timers::Timer timer_overwrittenmask("unusedpoints_mask");
      timer_overwrittenmask.start();

      // Declare this here to save it for later use. Contains all the boxes
      // which are active minus the boundary
      ibset all_refined;

      if (rl > 0) {
        int const orl = rl - 1;
        full_cboxes const &full_olevel = full_boxes.AT(ml).AT(orl);
        // Local boxes are not communicated or post-processed, and
        // thus can be modified even for coarser levels
        local_cboxes &local_olevel = local_boxes.AT(ml).AT(orl);

        ivect const reffact = h.reffacts.AT(rl) / h.reffacts.AT(orl);

        // This works only when the refinement factor is 2
        assert(all(reffact == 2));

        ibbox const &base = domain_exterior;
        ibbox const &obase = h.baseextent(ml, orl);

        // Calculate the union of all coarse regions
        ibset const allointr(full_olevel, &full_dboxes::interior);

        // Project to current level
        ivect const rf(reffact);
        ibset const parent(allointr.expanded_for(base));

        // Subtract the active region
        ibset const notrefined = parent - allactive;

        // Enlarge this set
        vect<vect<ibset, 2>, dim> enlarged;
        for (int d = 0; d < dim; ++d) {
          for (int f = 0; f < 2; ++f) {
            switch (h.refcent) {
            case vertex_centered: {
              ivect const dir = ivect::dir(d);
              enlarged[d][f] =
                  ibset(notrefined.expand(f == 1 ? dir : 0, f == 0 ? dir : 0));
              break;
            }
            case cell_centered: {
              enlarged[d][f] = notrefined;
              // TODO: restriction boundaries are wrong (they are
              // empty, but should not be) with cell centring when
              // fine cell cut coarse cells
              bool const old_there_was_an_error = there_was_an_error;
              ASSERT_rl(notrefined.contracted_for(obase).expanded_for(base) ==
                            notrefined,
                        "Refinement mask: Fine grid boundaries must be aligned "
                        "with coarse grid cells");
              // Do not abort because of this problem
              there_was_an_error = old_there_was_an_error;
              break;
            }
            default:
              assert(0);
            }
          }
        }

        // Intersect with the active region
        vect<vect<ibset, 2>, dim> all_boundaries;
        for (int d = 0; d < dim; ++d) {
          for (int f = 0; f < 2; ++f) {
            all_boundaries[d][f] = allactive & enlarged[d][f];
          }
        }

        // TODO: Ensure that the prolongation boundaries
        // all_boundaries are contained in the boundary prolongated
        // region

        // TODO: Ensure that the restriction boundaries and the
        // restricted region are contained in the restricted region

        // Subtract the boundaries from the refined region
        all_refined = allactive;
        for (int d = 0; d < dim; ++d) {
          for (int f = 0; f < 2; ++f) {
            all_refined -= all_boundaries[d][f];
          }
        }

        for (int lc = 0; lc < h.local_components(orl); ++lc) {
          int const c = h.get_component(orl, lc);
          full_dboxes const &obox = full_olevel.AT(c);
          // Local boxes are not communicated or post-processed, and
          // thus can be modified even for coarser levels
          local_dboxes &local_obox = local_olevel.AT(lc);

          // Set restriction information for next coarser level
          local_obox.restricted_region =
              all_refined.contracted_for(obox.exterior) & obox.owned;
        } // for lc

      } // if not coarsest level

      if (rl > 0) {
        int const orl = rl - 1;
        full_cboxes const &full_olevel = full_boxes.AT(ml).AT(orl);
        // Local boxes are not communicated or post-processed, and
        // thus can be modified even for coarser levels
        local_cboxes &local_olevel = local_boxes.AT(ml).AT(orl);

        // This works only when the refinement factor is 2
        ivect const reffact = h.reffacts.AT(rl) / h.reffacts.AT(orl);
        if (all(reffact == 2)) {
          // use the already computed 'all_refined' to get region from where
          // no information will be used later (overwritten)
          // First: get the region which will get restricted, on the coarse
          // level
          ibset restricted_region =
              all_refined.contracted_for(h.baseextent(ml, orl));
          // This is too big - during MoL-substeps information within this
          // region will be used to update points outside -> need to
          // shrink it by some points
          // The way we shrink it is to invert it, expand that, and invert
          // again. To invert we define an enclosing box and subtract it from
          // that.
          i2vect to_shrink = buffer_widths[orl] + ghost_widths[orl];
          ibbox enclosing =
              restricted_region.container().expand(ivect(1) + to_shrink);
          ibset unused_region =
              enclosing - (enclosing - restricted_region).expand(to_shrink);
          // Now we have the interesting region in 'unused_region' and need to
          // store
          // the union of this and the local regions
          for (int lc = 0; lc < h.local_components(orl); ++lc) {
            int const c = h.get_component(orl, lc);
            full_dboxes const &obox = full_olevel.AT(c);
            // Local boxes are not communicated or post-processed, and
            // thus can be modified even for coarser levels
            local_dboxes &local_obox = local_olevel.AT(lc);
            // Set unused information for next coarser level
            local_obox.unused_region = unused_region & obox.owned;
          } // for lc
        }   // if reffact != 2
      }     // if not coarsest level

      timer_overwrittenmask.stop();

      // Regridding schedule:

      fast_level.do_init = do_init;
      if (do_init) {

        static Timers::Timer timer_regrid("regrid");
        timer_regrid.start();

        static Timers::Timer timer_regrid_sync("sync");
        static Timers::Timer timer_regrid_prolongate("prolongate");

        for (int lc = 0; lc < h.local_components(rl); ++lc) {
          int const c = h.get_component(rl, lc);

          full_dboxes &box = full_level.AT(c);

          ibset needrecv = box.active + box.overlaps;

          // Regridding synchronisation:

          timer_regrid_sync.start();

          if (int(old_light_boxes.size()) > ml and
              int(old_light_boxes.AT(ml).size()) > rl) {

            int const oldcomponents = old_light_boxes.AT(ml).AT(rl).size();

            // Synchronisation copies from the same level of the old
            // grid structure.  It should fill as many active points
            // as possible.

            for (int cc = 0; cc < oldcomponents; ++cc) {
              if (needrecv.empty())
                break;

              light_dboxes const &obox = old_light_boxes.AT(ml).AT(rl).AT(cc);

              ibset const ovlp = needrecv & obox.owned;

              for (ibset::const_iterator ri = ovlp.begin(); ri != ovlp.end();
                   ++ri) {
                ibbox const &recv = *ri;
                ibbox const &send = recv;
                fast_level.fast_old2new_sync_sendrecv.push_back(
                    sendrecv_pseudoregion_t(send, cc, recv, c));
                if (not on_this_oldproc(rl, cc)) {
                  fast_dboxes &fast_level_otherproc =
                      fast_level_otherprocs.AT(this_oldproc(rl, cc));
                  fast_level_otherproc.fast_old2new_sync_sendrecv.push_back(
                      sendrecv_pseudoregion_t(send, cc, recv, c));
                }
              }

              needrecv -= ovlp;

            } // for cc

          } // if not old_light_boxes.empty

          timer_regrid_sync.stop();

          // Regridding prolongation:

          timer_regrid_prolongate.start();

          if (rl > 0) {
            int const orl = rl - 1;

            // Prolongation interpolates from the next coarser level
            // of the new grid structure.  It must fill what cannot be
            // synchronised.

            i2vect const stencil_size = i2vect(prolongation_stencil_size(rl));

            ASSERT_c(
                all(h.reffacts.at(rl) % h.reffacts.at(orl) == 0),
                "Refinement factors must be integer multiples of each other");
            i2vect const reffact =
                i2vect(h.reffacts.at(rl) / h.reffacts.at(orl));

            for (int cc = 0; cc < h.components(orl); ++cc) {
              if (needrecv.empty())
                break;

              full_dboxes const &obox = full_boxes.AT(ml).AT(orl).AT(cc);

              // See "refinement prolongation"
              ibset const expanded_oactive(
                  obox.active.expanded_for(box.interior)
                      .expand(h.refcent == vertex_centered ? reffact
                                                           : reffact - 1));
              ibset const ovlp = needrecv & expanded_oactive;

              for (ibset::const_iterator ri = ovlp.begin(); ri != ovlp.end();
                   ++ri) {
                ibbox const &recv = *ri;
                ibbox const send =
                    recv.expanded_for(obox.interior).expand(stencil_size);
                if (not(send <= obox.exterior)) {
                  cerr << "stencil_size=" << stencil_size << "\n"
                       << "obox=" << obox << "recv=" << recv << "\n"
                       << "send=" << send << "\n";
                }
                ASSERT_c(send <= obox.exterior, "Regridding prolongation: Send "
                                                "region must be contained in "
                                                "exterior");
                fast_level.fast_old2new_ref_prol_sendrecv.push_back(
                    sendrecv_pseudoregion_t(send, cc, recv, c));
                if (not on_this_proc(orl, cc)) {
                  fast_dboxes &fast_level_otherproc =
                      fast_level_otherprocs.AT(this_proc(orl, cc));
                  fast_level_otherproc.fast_old2new_ref_prol_sendrecv.push_back(
                      sendrecv_pseudoregion_t(send, cc, recv, c));
                }
              }

              needrecv -= ovlp;

            } // for cc

          } // if rl > 0

          if (int(old_light_boxes.size()) > ml and
              int(old_light_boxes.AT(ml).size()) > 0) {
            // All points must now have been received, either through
            // synchronisation or through prolongation
            ASSERT_c(
                needrecv.empty(),
                "Regridding prolongation: All points must have been received");
          }

          timer_regrid_prolongate.stop();

        } // for lc

        timer_regrid.stop();

      } // if do_init

      for (int lc = 0; lc < h.local_components(rl); ++lc) {
        int const c = h.get_component(rl, lc);

        light_level.AT(c).exterior = full_level.AT(c).exterior;
        light_level.AT(c).owned = full_level.AT(c).owned;
        light_level.AT(c).interior = full_level.AT(c).interior;

        light_level.AT(c).exterior_size = full_level.AT(c).exterior.size();
        light_level.AT(c).owned_size = full_level.AT(c).owned.size();
        light_level.AT(c).active_size = full_level.AT(c).active.size();

      } // for lc

      // Broadcast grid structure and communication schedule

      {

        static Timers::Timer timer_bcast_boxes("bcast_boxes");
        timer_bcast_boxes.start();

        int const count_send = h.local_components(rl);
        vector<light_dboxes> level_send(count_send);
        for (int lc = 0; lc < h.local_components(rl); ++lc) {
          int const c = h.get_component(rl, lc);
          level_send.AT(lc) = light_level.AT(c);
        }
        vector<vector<light_dboxes> > const level_recv =
            allgatherv(dist::comm(), level_send);
        vector<int> count_recv(dist::size(), 0);
        for (int c = 0; c < h.components(rl); ++c) {
          int const p = this_proc(rl, c);
          if (p != dist::rank()) {
            light_level.AT(c) = level_recv.AT(p).AT(count_recv.AT(p));
            ++count_recv.AT(p);
          }
        }
        for (int p = 0; p < dist::size(); ++p) {
          if (p != dist::rank()) {
            assert(count_recv.AT(p) == int(level_recv.AT(p).size()));
          }
        }

        timer_bcast_boxes.stop();
      }

      {

        static Timers::Timer timer_bcast_comm("bcast_comm");
        timer_bcast_comm.start();

        static Timers::Timer timer_bcast_comm_ref_prol("ref_prol");
        timer_bcast_comm_ref_prol.start();
        broadcast_schedule(fast_level_otherprocs, fast_level,
                           &fast_dboxes::fast_ref_prol_sendrecv);
        timer_bcast_comm_ref_prol.stop();

        static Timers::Timer timer_bcast_comm_sync("sync");
        timer_bcast_comm_sync.start();
        broadcast_schedule(fast_level_otherprocs, fast_level,
                           &fast_dboxes::fast_sync_sendrecv);
        timer_bcast_comm_sync.stop();

        static Timers::Timer timer_bcast_comm_ref_bnd_prol("ref_bnd_prol");
        timer_bcast_comm_ref_bnd_prol.start();
        broadcast_schedule(fast_level_otherprocs, fast_level,
                           &fast_dboxes::fast_ref_bnd_prol_sendrecv);
        timer_bcast_comm_ref_bnd_prol.stop();

        if (rl > 0) {
          int const orl = rl - 1;
          fast_dboxes &fast_olevel = fast_boxes.AT(ml).AT(orl);
          static Timers::Timer timer_bcast_comm_ref_rest("ref_rest");
          timer_bcast_comm_ref_rest.start();
          broadcast_schedule(fast_level_otherprocs, fast_olevel,
                             &fast_dboxes::fast_ref_rest_sendrecv);
          timer_bcast_comm_ref_rest.stop();
        }

        if (rl > 0) {
          int const orl = rl - 1;
          fast_dboxes &fast_olevel = fast_boxes.AT(ml).AT(orl);
          static Timers::Timer timer_bcast_comm_ref_refl("ref_refl");
          timer_bcast_comm_ref_refl.start();
          for (int dir = 0; dir < dim; ++dir) {
            for (int face = 0; face < 2; ++face) {
              srpvect fast_dboxes::*const fast_ref_refl_sendrecv =
                  fast_dboxes::fast_ref_refl_sendrecv[dir][face];
              broadcast_schedule(fast_level_otherprocs, fast_olevel,
                                 fast_ref_refl_sendrecv);
            }
          }
          timer_bcast_comm_ref_refl.stop();
        }

        static Timers::Timer timer_bcast_comm_ref_refl_prol("ref_refl_prol");
        timer_bcast_comm_ref_refl_prol.start();
        for (int dir = 0; dir < dim; ++dir) {
          for (int face = 0; face < 2; ++face) {
            srpvect fast_dboxes::*const fast_ref_refl_prol_sendrecv =
                fast_dboxes::fast_ref_refl_prol_sendrecv[dir][face];
            broadcast_schedule(fast_level_otherprocs, fast_level,
                               fast_ref_refl_prol_sendrecv);
          }
        }
        timer_bcast_comm_ref_refl_prol.stop();

        // TODO: Maybe broadcast old2new schedule only if do_init is
        // set
        static Timers::Timer timer_bcast_comm_old2new_sync("old2new_sync");
        timer_bcast_comm_old2new_sync.start();
        broadcast_schedule(fast_level_otherprocs, fast_level,
                           &fast_dboxes::fast_old2new_sync_sendrecv);
        timer_bcast_comm_old2new_sync.stop();

        static Timers::Timer timer_bcast_comm_old2new_ref_prol(
            "old2new_ref_prol");
        timer_bcast_comm_old2new_ref_prol.start();
        broadcast_schedule(fast_level_otherprocs, fast_level,
                           &fast_dboxes::fast_old2new_ref_prol_sendrecv);
        timer_bcast_comm_old2new_ref_prol.stop();

        timer_bcast_comm.stop();
      }

      // Output:
      if (output_bboxes or there_was_an_error) {

        cout << eol;
        cout << "ml=" << ml << " rl=" << rl << eol;
        cout << "baseextent=" << h.baseextent(ml, rl) << eol;

        for (int c = 0; c < h.components(rl); ++c) {
          cout << eol;
          cout << "ml=" << ml << " rl=" << rl << " c=" << c << eol;
          cout << "extent=" << h.extent(ml, rl, c) << eol;
          cout << "outer_boundaries=" << h.outer_boundaries(rl, c) << eol;
          cout << "processor=" << h.processor(rl, c) << eol;
        } // for c

        for (int c = 0; c < h.components(rl); ++c) {
          full_dboxes const &box = full_boxes.AT(ml).AT(rl).AT(c);
          cout << eol;
          cout << "ml=" << ml << " rl=" << rl << " c=" << c << eol;
          cout << box;
        } // for c

        // light, local, and fast boxes are output later

      } // if output_bboxes

      // Free memory early to save space
      if (int(old_light_boxes.size()) > ml and
          int(old_light_boxes.AT(ml).size()) > rl) {
        old_light_boxes.AT(ml).AT(rl).clear();
      }

      if (ml > 0) {
        if (rl > 0) {
          full_boxes.AT(ml - 1).AT(rl - 1).clear();
        }
        if (rl == h.reflevels() - 1) {
          full_boxes.AT(ml - 1).AT(rl).clear();
        }
      }
      if (ml == h.mglevels() - 1) {
        if (rl > 0) {
          full_boxes.AT(ml).AT(rl - 1).clear();
        }
        if (rl == h.reflevels() - 1) {
          full_boxes.AT(ml).AT(rl).clear();
        }
      }

    } // for rl

    if (ml > 0) {
      full_boxes.AT(ml - 1).clear();
    }
    if (ml == h.mglevels() - 1) {
      full_boxes.AT(ml).clear();
    }

  } // for ml

  // Output:
  if (output_bboxes or there_was_an_error) {

    for (int ml = 0; ml < h.mglevels(); ++ml) {
      for (int rl = 0; rl < h.reflevels(); ++rl) {

        cout << eol;
        cout << "ml=" << ml << " rl=" << rl << eol;
        cout << "baseextent=" << h.baseextent(ml, rl) << eol;

        for (int c = 0; c < h.components(rl); ++c) {
          cout << eol;
          cout << "ml=" << ml << " rl=" << rl << " c=" << c << eol;
          cout << "extent=" << h.extent(ml, rl, c) << eol;
          cout << "outer_boundaries=" << h.outer_boundaries(rl, c) << eol;
          cout << "processor=" << h.processor(rl, c) << eol;
        } // for c

        // full boxes have already been output (and deallocated)

        for (int c = 0; c < h.components(rl); ++c) {
          light_dboxes const &box = light_boxes.AT(ml).AT(rl).AT(c);
          cout << eol;
          cout << "ml=" << ml << " rl=" << rl << " c=" << c << eol;
          cout << box;
        } // for c

        for (int lc = 0; lc < h.local_components(rl); ++lc) {
          int const c = h.get_component(rl, lc);
          local_dboxes const &box = local_boxes.AT(ml).AT(rl).AT(lc);
          cout << eol;
          cout << "ml=" << ml << " rl=" << rl << " lc=" << lc << " c=" << c
               << eol;
          cout << box;
        } // for lc

        {
          level_dboxes const &box = level_boxes.AT(ml).AT(rl);
          cout << eol;
          cout << "ml=" << ml << " rl=" << rl << eol;
          cout << box;
        }

        fast_dboxes const &fast_box = fast_boxes.AT(ml).AT(rl);
        cout << eol;
        cout << "ml=" << ml << " rl=" << rl << eol;
        cout << fast_box;

      } // for rl
    }   // for ml

    cout << eol;
    cout << "memoryof(gh)=" << memoryof(h) << eol;
    cout << "memoryof(dh)=" << memoryof(*this) << eol;
    cout << "memoryof(dh.light_boxes)=" << memoryof(light_boxes) << eol;
    cout << "memoryof(dh.local_boxes)=" << memoryof(local_boxes) << eol;
    cout << "memoryof(dh.level_boxes)=" << memoryof(level_boxes) << eol;
    cout << "memoryof(dh.fast_boxes)=" << memoryof(fast_boxes) << eol;
    int gfcount = 0;
    size_t gfmemory = 0;
    for (map<int, ggf *>::const_iterator gfi = gfs.begin(); gfi != gfs.end();
         ++gfi) {
      ++gfcount;
      gfmemory += memoryof(*gfi->second);
    }
    cout << "#gfs=" << gfcount << eol;
    cout << "memoryof(gfs)=" << gfmemory << eol;

  } // if output_bboxes

  if (there_was_an_error) {
    CCTK_WARN(
        CCTK_WARN_ABORT,
        "The grid structure is inconsistent.  It is impossible to continue.");
  }

  total.stop(0);
  timer.stop();
}

void dh::broadcast_schedule(vector<fast_dboxes> &fast_level_otherprocs,
                            fast_dboxes &fast_level,
                            srpvect fast_dboxes::*const schedule_item) {
  static Timers::Timer timer_bs1("CarpetLib::dh::bs1");
  timer_bs1.start();
  vector<srpvect> send(dist::size());
  for (int p = 0; p < dist::size(); ++p) {
    swap(send.AT(p), fast_level_otherprocs.AT(p).*schedule_item);
  }
  timer_bs1.stop();

  static Timers::Timer timer_bs2("CarpetLib::dh::bs2");
  timer_bs2.start();
  srpvect const recv = alltoallv1(dist::comm(), send);
  timer_bs2.stop();

  static Timers::Timer timer_bs3("CarpetLib::dh::bs3");
  timer_bs3.start();
  (fast_level.*schedule_item)
      .insert((fast_level.*schedule_item).end(), recv.begin(), recv.end());
  timer_bs3.stop();
}

void dh::regrid_free(bool const do_init) {
  if (do_init) {
    for (int ml = 0; ml < h.mglevels(); ++ml) {
      for (int rl = 0; rl < h.reflevels(); ++rl) {
        fast_boxes.AT(ml).AT(rl).fast_old2new_sync_sendrecv.clear();
        fast_boxes.AT(ml).AT(rl).fast_old2new_ref_prol_sendrecv.clear();
      }
    }
  } else {
    for (int ml = 0; ml < h.mglevels(); ++ml) {
      for (int rl = 0; rl < h.reflevels(); ++rl) {
        assert(fast_boxes.AT(ml).AT(rl).fast_old2new_sync_sendrecv.empty());
        assert(fast_boxes.AT(ml).AT(rl).fast_old2new_ref_prol_sendrecv.empty());
      }
    }
  }
}

void dh::recompose(int const rl, bool const do_prolongate) {
  DECLARE_CCTK_PARAMETERS;

  assert(rl >= 0 and rl < h.reflevels());

  static Timers::Timer timer("CarpetLib::dh::recompose");
  timer.start();

  for (map<int, ggf *>::iterator f = gfs.begin(); f != gfs.end(); ++f) {
    f->second->recompose_crop();
  }

  if (combine_recompose) {
    // Recompose all grid functions of this refinement levels at once.
    // This may be faster, but requires more memory. This is the default.
    for (map<int, ggf *>::iterator f = gfs.begin(); f != gfs.end(); ++f) {
      f->second->recompose_allocate(rl);
    }
// TODO: If this works, rename do_prolongate to do_init here, and
// remove the do_prolongate parameter from ggf::recompose_fill
#if 0
    for (comm_state state; not state.done(); state.step()) {
      for (map<int,ggf*>::iterator f=gfs.begin(); f!=gfs.end(); ++f) {
        f->second->recompose_fill (state, rl, do_prolongate);
      }
    }
#endif
    if (do_prolongate) {
      for (comm_state state; not state.done(); state.step()) {
        for (map<int, ggf *>::iterator f = gfs.begin(); f != gfs.end(); ++f) {
          f->second->recompose_fill(state, rl, true);
        }
      }
    }
    for (map<int, ggf *>::iterator f = gfs.begin(); f != gfs.end(); ++f) {
      f->second->recompose_free_old(rl);
    }
  } else {
    // Recompose the grid functions sequentially.  This may be slower,
    // but requires less memory.
    for (map<int, ggf *>::iterator f = gfs.begin(); f != gfs.end(); ++f) {
      f->second->recompose_allocate(rl);
#if 0
      for (comm_state state; not state.done(); state.step()) {
        f->second->recompose_fill (state, rl, do_prolongate);
      }
#endif
      if (do_prolongate) {
        for (comm_state state; not state.done(); state.step()) {
          f->second->recompose_fill(state, rl, true);
        }
      }
      f->second->recompose_free_old(rl);
    }
  }

  timer.stop();
}

// Grid function management
void dh::insert(ggf *const f) {
  CHECKPOINT;
  assert(f);
  assert(f->varindex >= 0);
  assert(f->varindex < CCTK_NumVars());
  bool const inserted = gfs.insert(pair<int, ggf *>(f->varindex, f)).second;
  assert(inserted);
}

void dh::erase(ggf *const f) {
  CHECKPOINT;
  assert(f);
  assert(f->varindex >= 0);
  assert(f->varindex < CCTK_NumVars());
  size_t const erased = gfs.erase(f->varindex);
  assert(erased == 1);
}

// Equality

bool dh::full_dboxes::operator==(full_dboxes const &b) const {
  return exterior == b.exterior and
         all(all(is_outer_boundary == b.is_outer_boundary)) and
         outer_boundaries == b.outer_boundaries and
         communicated == b.communicated and boundaries == b.boundaries and
         owned == b.owned and buffers == b.buffers and
         overlaps == b.overlaps and active == b.active and sync == b.sync and
         bndref == b.bndref and ghosts == b.ghosts and interior == b.interior;
}

// MPI datatypes

MPI_Datatype mpi_datatype(dh::light_dboxes const &) {
  static bool initialised = false;
  static MPI_Datatype newtype;
  if (not initialised) {
    static dh::light_dboxes s;
#define ENTRY(type, name)                                                      \
  {                                                                            \
    sizeof s.name / sizeof(type),     /* count elements */                     \
        (char *)&s.name - (char *)&s, /* offsetof doesn't work (why?) */       \
        dist::mpi_datatype<type>(),   /* find MPI datatype */                  \
        STRINGIFY(name),              /* field name */                         \
        STRINGIFY(type),              /* type name */                          \
  }
    dist::mpi_struct_descr_t const descr[] = {
        ENTRY(int, exterior),
        ENTRY(int, owned),
        ENTRY(int, interior),
        ENTRY(size_type, exterior_size),
        ENTRY(size_type, owned_size),
        ENTRY(size_type, active_size),
        {1, sizeof s, MPI_UB, "MPI_UB", "MPI_UB"}};
#undef ENTRY
    newtype = dist::create_mpi_datatype(sizeof descr / sizeof descr[0], descr,
                                        "dh::light::dboxes", sizeof s);
#if 0
    int type_size;
    MPI_Type_size (newtype, & type_size);
    assert (type_size <= sizeof s);
    MPI_Aint type_lb, type_ub;
    MPI_Type_lb (newtype, & type_lb);
    MPI_Type_ub (newtype, & type_ub);
    assert (type_ub - type_lb == sizeof s);
#endif
    initialised = true;
  }
  return newtype;
}

MPI_Datatype mpi_datatype(dh::fast_dboxes const &) {
  static bool initialised = false;
  static MPI_Datatype newtype;
  if (not initialised) {
    static dh::fast_dboxes s;
#define ENTRY(type, name)                                                      \
  {                                                                            \
    sizeof s.name / sizeof(type),     /* count elements */                     \
        (char *)&s.name - (char *)&s, /* offsetof doesn't work (why?) */       \
        dist::mpi_datatype<type>(),   /* find MPI datatype */                  \
        STRINGIFY(name),              /* field name */                         \
        STRINGIFY(type),              /* type name */                          \
  }
    dist::mpi_struct_descr_t const descr[] = {
        ENTRY(dh::srpvect, fast_mg_rest_sendrecv),
        ENTRY(dh::srpvect, fast_mg_prol_sendrecv),
        ENTRY(dh::srpvect, fast_ref_prol_sendrecv),
        ENTRY(dh::srpvect, fast_ref_rest_sendrecv),
        ENTRY(dh::srpvect, fast_sync_sendrecv),
        ENTRY(dh::srpvect, fast_ref_bnd_prol_sendrecv),
        ENTRY(dh::srpvect, fast_old2new_sync_sendrecv),
        ENTRY(dh::srpvect, fast_old2new_ref_prol_sendrecv),
        ENTRY(dh::srpvect, fast_ref_refl_sendrecv_0_0),
        ENTRY(dh::srpvect, fast_ref_refl_sendrecv_0_1),
        ENTRY(dh::srpvect, fast_ref_refl_sendrecv_1_0),
        ENTRY(dh::srpvect, fast_ref_refl_sendrecv_1_1),
        ENTRY(dh::srpvect, fast_ref_refl_sendrecv_2_0),
        ENTRY(dh::srpvect, fast_ref_refl_sendrecv_2_1),
        ENTRY(dh::srpvect, fast_ref_refl_prol_sendrecv_0_0),
        ENTRY(dh::srpvect, fast_ref_refl_prol_sendrecv_0_1),
        ENTRY(dh::srpvect, fast_ref_refl_prol_sendrecv_1_0),
        ENTRY(dh::srpvect, fast_ref_refl_prol_sendrecv_1_1),
        ENTRY(dh::srpvect, fast_ref_refl_prol_sendrecv_2_0),
        ENTRY(dh::srpvect, fast_ref_refl_prol_sendrecv_2_1),
        {1, sizeof s, MPI_UB, "MPI_UB", "MPI_UB"}};
#undef ENTRY
    newtype = dist::create_mpi_datatype(sizeof descr / sizeof descr[0], descr,
                                        "dh::fast_dboxes", sizeof s);
    initialised = true;
  }
  return newtype;
}

// Memory usage

size_t dh::memory() const {
  return memoryof(ghost_widths) + memoryof(buffer_widths) +
         memoryof(overlap_widths) + memoryof(prolongation_orders_space) +
         memoryof(light_boxes) + memoryof(local_boxes) + memoryof(level_boxes) +
         memoryof(fast_boxes) + memoryof(gfs);
}

size_t dh::allmemory() {
  size_t mem = memoryof(alldh);
  for (set<dh *>::const_iterator dhi = alldh.begin(); dhi != alldh.end();
       ++dhi) {
    mem += memoryof(**dhi);
  }
  return mem;
}

size_t dh::light_dboxes::memory() const {
  return memoryof(exterior) + memoryof(owned) + memoryof(interior) +
         memoryof(exterior_size) + memoryof(owned_size) + memoryof(active_size);
}

size_t dh::local_dboxes::memory() const {
  return memoryof(buffers) + memoryof(overlaps) + memoryof(active) +
#if 0
    memoryof (buffers_stepped) +
#endif
#if 0
    memoryof (fine_active) +
#endif
         memoryof(prolongation_boundary) + memoryof(restriction_boundary) +
         memoryof(restricted_region) + memoryof(unused_region) +
         memoryof(coarse_boundary) + memoryof(fine_boundary);
#if 0 // OFFSETS,SIZE
    memoryof (coarse_boundary_offsets) +
    memoryof (fine_boundary_offsets) +
    memoryof (coarse_boundary_size) +
    memoryof (fine_boundary_size);
#endif
}

size_t dh::level_dboxes::memory() const { return memoryof(active); }

size_t dh::full_dboxes::memory() const {
  return memoryof(exterior) + memoryof(is_outer_boundary) +
         memoryof(outer_boundaries) + memoryof(communicated) +
         memoryof(boundaries) + memoryof(owned) + memoryof(buffers) +
         memoryof(overlaps) + memoryof(active) + memoryof(sync) +
         memoryof(bndref) + memoryof(ghosts) + memoryof(interior);
}

size_t dh::fast_dboxes::memory() const {
  return memoryof(fast_mg_rest_sendrecv) + memoryof(fast_mg_prol_sendrecv) +
         memoryof(fast_ref_prol_sendrecv) + memoryof(fast_ref_rest_sendrecv) +
         memoryof(fast_sync_sendrecv) + memoryof(fast_ref_bnd_prol_sendrecv) +
         memoryof(fast_ref_refl_sendrecv_0_0) +
         memoryof(fast_ref_refl_sendrecv_0_1) +
         memoryof(fast_ref_refl_sendrecv_1_0) +
         memoryof(fast_ref_refl_sendrecv_1_1) +
         memoryof(fast_ref_refl_sendrecv_2_0) +
         memoryof(fast_ref_refl_sendrecv_2_1) +
         memoryof(fast_ref_refl_prol_sendrecv_0_0) +
         memoryof(fast_ref_refl_prol_sendrecv_0_1) +
         memoryof(fast_ref_refl_prol_sendrecv_1_0) +
         memoryof(fast_ref_refl_prol_sendrecv_1_1) +
         memoryof(fast_ref_refl_prol_sendrecv_2_0) +
         memoryof(fast_ref_refl_prol_sendrecv_2_1) + memoryof(do_init) +
         memoryof(fast_old2new_sync_sendrecv) +
         memoryof(fast_old2new_ref_prol_sendrecv);
}

// Input

istream &dh::light_dboxes::input(istream &is) {
  // Regions:
  try {
    skipws(is);
    consume(is, "dh::light_dboxes:{");
    skipws(is);
    consume(is, "exterior:");
    is >> exterior;
    exterior_size = exterior.size();
    skipws(is);
    consume(is, "owned:");
    is >> owned;
    owned_size = owned.size();
    skipws(is);
    consume(is, "interior:");
    is >> interior;
    skipws(is);
    consume(is, "active_size:");
    is >> active_size;
    skipws(is);
    consume(is, "}");
  } catch (input_error &err) {
    cout << "Input error while reading a dh::light_dboxes" << endl;
    throw err;
  }
  return is;
}

istream &dh::local_dboxes::input(istream &is) {
  // Regions:
  try {
    skipws(is);
    consume(is, "dh::local_dboxes:{");
    skipws(is);
    consume(is, "buffers:");
    is >> buffers;
    skipws(is);
    consume(is, "overlaps:");
    is >> overlaps;
    skipws(is);
    consume(is, "active:");
    is >> active;
#if 0
    skipws (is);
    consume (is, "buffers_stepped:");
    is >> buffers_stepped;
#endif
#if 0
    skipws (is);
    consume (is, "fine_active:");
    is >> fine_active;
#endif
    skipws(is);
    consume(is, "prolongation_boundary:");
    is >> prolongation_boundary;
    skipws(is);
    consume(is, "restriction_boundary:");
    is >> restriction_boundary;
    skipws(is);
    consume(is, "restricted_region:");
    is >> restricted_region;
    skipws(is);
    consume(is, "unused_region:");
    is >> unused_region;
    skipws(is);
    consume(is, "coarse_boundary:");
    is >> coarse_boundary;
    skipws(is);
    consume(is, "fine_boundary:");
    is >> fine_boundary;
    skipws(is);
    // TODO: read boundary sizes and boundary offsets
    consume(is, "}");
  } catch (input_error &err) {
    cout << "Input error while reading a dh::local_dboxes" << endl;
    throw err;
  }
  return is;
}

istream &dh::level_dboxes::input(istream &is) {
  // Regions:
  try {
    skipws(is);
    consume(is, "dh::level_dboxes:{");
    skipws(is);
    consume(is, "active:");
    is >> active;
    skipws(is);
    consume(is, "}");
  } catch (input_error &err) {
    cout << "Input error while reading a dh::level_dboxes" << endl;
    throw err;
  }
  return is;
}

istream &dh::full_dboxes::input(istream &is) {
  // Regions:
  try {
    skipws(is);
    consume(is, "dh::full_dboxes:{");
    skipws(is);
    consume(is, "exterior:");
    is >> exterior;
    skipws(is);
    consume(is, "is_outer_boundary:");
    is >> is_outer_boundary;
    skipws(is);
    consume(is, "outer_boundaries:");
    is >> outer_boundaries;
    skipws(is);
    consume(is, "communicated:");
    is >> communicated;
    skipws(is);
    consume(is, "boundaries:");
    is >> boundaries;
    skipws(is);
    consume(is, "owned:");
    is >> owned;
    skipws(is);
    consume(is, "buffers:");
    is >> buffers;
    skipws(is);
    consume(is, "overlaps:");
    is >> overlaps;
    skipws(is);
    consume(is, "active:");
    is >> active;
    skipws(is);
    consume(is, "sync:");
    is >> sync;
    skipws(is);
    consume(is, "bndref:");
    is >> bndref;
    skipws(is);
    consume(is, "ghosts:");
    is >> ghosts;
    skipws(is);
    consume(is, "interior:");
    is >> interior;
    skipws(is);
    consume(is, "}");
  } catch (input_error &err) {
    cout << "Input error while reading a dh::full_dboxes" << endl;
    throw err;
  }
  return is;
}

istream &dh::fast_dboxes::input(istream &is) {
  // Communication schedule:
  try {
    skipws(is);
    consume(is, "dh::fast_dboxes:{");
    skipws(is);
    consume(is, "fast_mg_rest_sendrecv:");
    is >> fast_mg_rest_sendrecv;
    skipws(is);
    consume(is, "fast_mg_prol_sendrecv:");
    is >> fast_mg_prol_sendrecv;
    skipws(is);
    consume(is, "fast_ref_prol_sendrecv:");
    is >> fast_ref_prol_sendrecv;
    skipws(is);
    consume(is, "fast_ref_rest_sendrecv:");
    is >> fast_ref_rest_sendrecv;
    skipws(is);
    consume(is, "fast_sync_sendrecv:");
    is >> fast_sync_sendrecv;
    skipws(is);
    consume(is, "fast_ref_bnd_prol_sendrecv:");
    is >> fast_ref_bnd_prol_sendrecv;
    skipws(is);
    consume(is, "fast_old2new_sync_sendrecv:");
    is >> fast_old2new_sync_sendrecv;
    skipws(is);
    consume(is, "fast_old2new_ref_prol_sendrecv:");
    is >> fast_old2new_ref_prol_sendrecv;
    skipws(is);
    consume(is, "}");
  } catch (input_error &err) {
    cout << "Input error while reading a dh::fast_dboxes" << endl;
    throw err;
  }
  return is;
}

// Output

ostream &dh::output(ostream &os) const {
  os << "dh:"
     << "ghost_widths=" << ghost_widths << ","
     << "buffer_widths=" << buffer_widths << ","
     << "overlap_widths=" << overlap_widths << ","
     << "prolongation_orders_space=" << prolongation_orders_space << ","
     << "light_boxes=" << light_boxes << ","
     << "local_boxes=" << local_boxes << ","
     << "level_boxes=" << level_boxes << ","
     << "fast_boxes=" << fast_boxes << ","
     << "gfs={";
  {
    bool isfirst = true;
    for (map<int, ggf *>::const_iterator f = gfs.begin(); f != gfs.end();
         ++f, isfirst = false) {
      if (not isfirst)
        os << ",";
      os << *f;
    }
  }
  os << "}";
  return os;
}

ostream &dh::light_dboxes::output(ostream &os) const {
  // Regions:
  os << "dh::light_dboxes:{" << eol << "   exterior: " << exterior << eol
     << "   owned: " << owned << eol << "   interior: " << interior << eol
     << "   active_size: " << active_size << eol << "}" << eol;
  return os;
}

ostream &dh::local_dboxes::output(ostream &os) const {
  // Regions:
  os << "dh::local_dboxes:{" << eol << "   buffers: " << buffers << eol
     << "   overlaps: " << overlaps << eol << "   active: " << active << eol
#if 0
     << "   buffers_stepped: " << buffers_stepped << eol
#endif
#if 0
     << "   fine_active: " << fine_active << eol
#endif
     << "   prolongation_boundary: " << prolongation_boundary << eol
     << "   restriction_boundary: " << restriction_boundary << eol
     << "   restricted_region: " << restricted_region << eol
     << "   unused_region: " << unused_region << eol
     << "   coarse_boundary: " << coarse_boundary << eol
     << "   fine_boundary: " << fine_boundary << eol
#if 0 // OFFSETS,SIZE
     << "   coarse_boundary_offsets: " << coarse_boundary_offsets << eol
     << "   fine_boundary_offsets: " << fine_boundary_offsets << eol
     << "   coarse_boundary_size: " << coarse_boundary_size << eol
     << "   fine_boundary_size: " << fine_boundary_size << eol
#endif
     << "}" << eol;
  return os;
}

ostream &dh::level_dboxes::output(ostream &os) const {
  // Regions:
  os << "dh::level_dboxes:{" << eol << "   active: " << active << eol << "}"
     << eol;
  return os;
}

ostream &dh::full_dboxes::output(ostream &os) const {
  // Regions:
  os << "dh::full_dboxes:{" << eol << "   exterior: " << exterior << eol
     << "   is_outer_boundary: " << is_outer_boundary << eol
     << "   outer_boundaries: " << outer_boundaries << eol
     << "   communicated: " << communicated << eol
     << "   boundaries: " << boundaries << eol << "   owned: " << owned << eol
     << "   buffers: " << buffers << eol << "   overlaps: " << overlaps << eol
     << "   active: " << active << eol << "   sync: " << sync << eol
     << "   bndref: " << bndref << eol << "   ghosts: " << ghosts << eol
     << "   interior: " << interior << eol << "}" << eol;
  return os;
}

ostream &dh::fast_dboxes::output(ostream &os) const {
  // Communication schedule:
  os << "dh::fast_dboxes:{" << eol
     << "   fast_mg_rest_sendrecv: " << fast_mg_rest_sendrecv << eol
     << "   fast_mg_prol_sendrecv: " << fast_mg_prol_sendrecv << eol
     << "   fast_ref_prol_sendrecv: " << fast_ref_prol_sendrecv << eol
     << "   fast_ref_rest_sendrecv: " << fast_ref_rest_sendrecv << eol
     << "   fast_sync_sendrecv: " << fast_sync_sendrecv << eol
     << "   fast_ref_bnd_prol_sendrecv: " << fast_ref_bnd_prol_sendrecv << eol
     << "   fast_ref_refl_sendrecv_0_0: " << fast_ref_refl_sendrecv_0_0 << eol
     << "   fast_ref_refl_sendrecv_0_1: " << fast_ref_refl_sendrecv_0_1 << eol
     << "   fast_ref_refl_sendrecv_1_0: " << fast_ref_refl_sendrecv_1_0 << eol
     << "   fast_ref_refl_sendrecv_1_1: " << fast_ref_refl_sendrecv_1_1 << eol
     << "   fast_ref_refl_sendrecv_2_0: " << fast_ref_refl_sendrecv_2_0 << eol
     << "   fast_ref_refl_sendrecv_2_1: " << fast_ref_refl_sendrecv_2_1 << eol
     << "   fast_ref_refl_prol_sendrecv_0_0: "
     << fast_ref_refl_prol_sendrecv_0_0 << eol
     << "   fast_ref_refl_prol_sendrecv_0_1: "
     << fast_ref_refl_prol_sendrecv_0_1 << eol
     << "   fast_ref_refl_prol_sendrecv_1_0: "
     << fast_ref_refl_prol_sendrecv_1_0 << eol
     << "   fast_ref_refl_prol_sendrecv_1_1: "
     << fast_ref_refl_prol_sendrecv_1_1 << eol
     << "   fast_ref_refl_prol_sendrecv_2_0: "
     << fast_ref_refl_prol_sendrecv_2_0 << eol
     << "   fast_ref_refl_prol_sendrecv_2_1: "
     << fast_ref_refl_prol_sendrecv_2_1 << eol << "   do_init: " << do_init
     << eol << "   fast_old2new_sync_sendrecv: " << fast_old2new_sync_sendrecv
     << eol
     << "   fast_old2new_ref_prol_sendrecv: " << fast_old2new_ref_prol_sendrecv
     << eol << "}" << eol;
  return os;
}
