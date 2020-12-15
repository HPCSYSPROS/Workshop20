#include <cctk.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <list>
#include <sstream>
#include <string>
#include <vector>

#include <sys/stat.h>
#include <sys/types.h>

#ifdef CCTK_MPI
#include <mpi.h>
#else
#include "nompi.h"
#endif

#include <Timer.hh>

#include <loopcontrol.h>

#include <balance.hh>
#include <bbox.hh>
#include <bboxset.hh>
#include <defs.hh>
#include <dh.hh>
#include <gh.hh>
#include <region.hh>
#include <vect.hh>

#include <carpet.hh>
#include <modes.hh>
#include <variables.hh>

namespace Carpet {

using namespace std;

// Helper routines for spliting regions automatically

// The cost for a region, assuming a cost of 1 per interior point
static rvect cost(region_t const &reg) {
  DECLARE_CCTK_PARAMETERS;
  static rvect costfactor;
  static bool initialised = false;
  if (not initialised) {
    costfactor = rvect(1.0);
    if (dim > 0)
      costfactor[0] = 1.0 / aspect_ratio_x;
    if (dim > 1)
      costfactor[1] = 1.0 / aspect_ratio_y;
    if (dim > 2)
      costfactor[2] = 1.0 / aspect_ratio_z;
  }
  if (reg.extent.empty())
    return rvect(0);
  return rvect(reg.extent.shape() / reg.extent.stride()) * costfactor;
}

static void ClassifyPoints(cGH const *cctkGH, int rl);

static void SplitRegionsMaps_Automatic_Recursively(bvect const &dims,
                                                   int const firstproc,
                                                   int const nprocs,
                                                   region_t &superreg,
                                                   vector<region_t> &newregs);
static void SplitRegionsMaps_AsSpecified(cGH const *cctkGH,
                                         vector<vector<region_t> > &superregss,
                                         vector<vector<region_t> > &regss);
static void SplitRegions_AsSpecified(cGH const *cctkGH,
                                     vector<region_t> &superregs,
                                     vector<region_t> &regs);

void CheckRegions(gh::mregs const &regsss) {
  char const *const where = "CheckRegions";
  static Timers::Timer timer(where);
  timer.start();

  // At least one multigrid level
  if (regsss.size() == 0) {
    CCTK_WARN(0,
              "I cannot set up a grid hierarchy with zero multigrid levels.");
  }
  assert(regsss.size() > 0);
  // At least one refinement level
  if (regsss.AT(0).size() == 0) {
    CCTK_WARN(0,
              "I cannot set up a grid hierarchy with zero refinement levels.");
  }
  assert(regsss.AT(0).size() > 0);
  // At most maxreflevels levels
  if ((int)regsss.AT(0).size() > maxreflevels) {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "I cannot set up a grid hierarchy with more than "
               "Carpet::max_refinement_levels refinement levels.  I found "
               "Carpet::max_refinement_levels=%d, while %d levels were "
               "requested.",
               (int)maxreflevels, (int)regsss.AT(0).size());
  }
  assert((int)regsss.AT(0).size() <= maxreflevels);
  for (int ml = 0; ml < (int)regsss.size(); ++ml) {
    int num_regions = 0;
    for (int rl = 0; rl < (int)regsss.AT(0).size(); ++rl) {
      // No empty levels
      // (but allow some empty maps)
      // assert (regsss.AT(ml).AT(rl).size() > 0);
      num_regions += regsss.AT(ml).AT(rl).size();
      for (int c = 0; c < (int)regsss.AT(ml).AT(rl).size(); ++c) {
        // Check sizes
        // (but allow processors with zero grid points)
        // assert (all(regsss.AT(rl).AT(c).AT(ml).extent.lower() <=
        //             regsss.AT(rl).AT(c).AT(ml).extent.upper()));
        // Check strides
        ivect const str =
            (maxspacereflevelfact / spacereffacts.AT(rl) * ipow(mgfact, ml));
        assert(all(regsss.AT(ml).AT(rl).AT(c).extent.stride() % str == 0));
        // Check alignments
        assert(all(regsss.AT(ml).AT(rl).AT(c).extent.lower() % str == 0));
        assert(all(regsss.AT(ml).AT(rl).AT(c).extent.upper() % str == 0));
      }
    }
    // No empty levels
    assert(num_regions > 0);
  }

  timer.stop();
}

bool Regrid(cGH const *const cctkGH, bool const force_recompose,
            bool const do_init) {
  DECLARE_CCTK_PARAMETERS;

  char const *const where = "Regrid";
  static Timers::Timer timer(where);
  timer.start();

  Checkpoint("Regridding level %d...", reflevel);

  assert(is_level_mode());

  bool const have_regrid = CCTK_IsFunctionAliased("Carpet_Regrid");
  bool const have_regridmaps = CCTK_IsFunctionAliased("Carpet_RegridMaps");
  bool const use_regridmaps = regrid_in_level_mode and have_regridmaps;

  if (not use_regridmaps and not have_regrid) {
    static bool didtell = false;
    if (maxreflevels > 1 and not didtell) {
      CCTK_WARN(2, "The regridding routine Carpet_Regrid has not been "
                   "provided.  There will be no regridding.  Maybe you forgot "
                   "to activate a regridding thorn?");
      didtell = true;
    }
    timer.stop();
    return false;
  }

  if (regrid_in_level_mode and not have_regridmaps) {
    static bool didtell = false;
    if (maps > 1 and maxreflevels > 1 and not didtell) {
      CCTK_WARN(2, "The regridding routine Carpet_RegridMaps has not been "
                   "provided.  Regridding will be performed in singlemap mode "
                   "instead of level mode.");
      didtell = true;
    }
  }

  bool did_change = false;

  if (not use_regridmaps) {

    BEGIN_MAP_LOOP(cctkGH, CCTK_GF) {

      gh::rregs superregss = vhh.AT(map)->superregions;
      gh::mregs regsss;

      // Check whether to recompose
      CCTK_INT const do_recompose =
          Carpet_Regrid(cctkGH, &superregss, &regsss, force_recompose);
      assert(do_recompose >= 0);
      did_change = did_change or do_recompose;

      if (do_recompose) {
        RegridMap(cctkGH, map, superregss, regsss, do_init);
      }
    }
    END_MAP_LOOP;

  } else {

    vector<gh::rregs> superregsss(maps);
    vector<gh::mregs> regssss(maps);
    BEGIN_MAP_LOOP(cctkGH, CCTK_GF) {
      superregsss.AT(map) = vhh.AT(map)->superregions;
    }
    END_MAP_LOOP;

    // Check whether to recompose
    CCTK_INT const do_recompose =
        Carpet_RegridMaps(cctkGH, &superregsss, &regssss, force_recompose);
    assert(do_recompose >= 0);
// TODO
#if 1 // #ifdef CARPET_DEBUG
    {
      int ival = do_recompose;
      MPI_Bcast(&ival, 1, MPI_INT, 0, dist::comm());
      assert(ival == do_recompose);
    }
#endif
    did_change = did_change or do_recompose;

    if (do_recompose) {
      BEGIN_MAP_LOOP(cctkGH, CCTK_GF) {
        gh::rregs const &superregss = superregsss.AT(map);
        gh::mregs const &regsss = regssss.AT(map);
        RegridMap(cctkGH, map, superregss, regsss, do_init);
      }
      END_MAP_LOOP;
    }

  } // if regrid in level mode

  if (did_change) {

    PostRegrid(cctkGH);

  } // if did change

  timer.stop();

  return did_change;
}

void RegridMap(cGH const *const cctkGH, int const m,
               gh::rregs const &superregss, gh::mregs const &regsss,
               bool const do_init) {
  DECLARE_CCTK_PARAMETERS;

  char const *const where = "RegridMap";
  static Timers::Timer timer(where);
  timer.start();

  Waypoint("Regridding map %d...", m);

// TODO: keep levels fixed here
#if 0
    //
    // Keep this level fixed if it is not evolved
    //
    if (regrid_only_evolved_levels) {
      int type;
      bool const use_tapered_grids =
        * static_cast<CCTK_INT const *>
        (CCTK_ParameterGet ("use_tapered_grids", "Carpet", &type));
      assert (type == PARAMETER_BOOLEAN);
      
      int const do_every =
        use_tapered_grids ?
        maxtimereflevelfact / timereffacts.AT(max(0,rl-1)):
        maxtimereflevelfact / timereffacts.AT(      rl   );
      
      bool const regrid_this_level =
        (cctkGH->cctk_iteration - 1) % do_every == 0;
      
      if (not regrid_this_level) {
        // Set regions from current grid structure
        regions.AT(rl) = ...;
      }
    }
#endif

  // Check the regions
  CheckRegions(regsss);
  // TODO: check also that the current and all coarser levels did
  // not change

  // Regrid
  vhh.AT(m)->regrid(superregss, regsss, do_init);

  if (print_grid_info) {
    // Output grid structure to screen
    OutputSuperregions(cctkGH, m, *vhh.AT(m), superregss);
    OutputGrids(cctkGH, m, *vhh.AT(m), *vdd.AT(m));
  }

  // Write grid structure to file
  OutputGridStructure(cctkGH, m, regsss);
  OutputGridCoordinates(cctkGH, m, regsss);

  Waypoint("Done regridding map %d.", m);
  timer.stop();
}

void PostRegrid(cGH const *const cctkGH) {
  DECLARE_CCTK_PARAMETERS;

  char const *const where = "PostRegrid";
  static Timers::Timer timer(where);
  timer.start();

  // Calculate new number of levels
  // int const oldreflevels = reflevels;
  reflevels = vhh.AT(0)->reflevels();
  for (int m = 0; m < maps; ++m) {
    assert(vhh.AT(m)->reflevels() == reflevels);
  }

// TODO
#if 1 // #ifdef CARPET_DEBUG
  {
    // All processes must use the same number of levels
    int ival = reflevels;
    MPI_Bcast(&ival, 1, MPI_INT, 0, dist::comm());
    assert(ival == reflevels);
  }
#endif

  // One cannot switch off the current level
  assert(reflevels > reflevel);

  // Set new number of active time levels
  for (int n = 0; n < CCTK_NumGroups(); ++n) {
    int const grouptype = CCTK_GroupTypeI(n);
    if (grouptype == CCTK_GF) {
      for (int ml = 0; ml < mglevels; ++ml) {
        groupdata.AT(n).activetimelevels.AT(ml).resize(
            reflevels, groupdata.AT(n).activetimelevels.AT(ml).AT(0));
      }
    }
  }

  // Set new level times
  // for (int ml=0; ml<mglevels; ++ml) {
  //   leveltimes.AT(ml).resize
  //     (reflevels, leveltimes.AT(ml).AT(oldreflevels-1));
  // }

  ++regridding_epoch;
  // Remove levels that vanished
  if (int(level_regridding_epochs.size()) > reflevels) {
    level_regridding_epochs.resize(reflevels);
  }
  // Insert entries for new levels
  if (int(level_regridding_epochs.size()) < reflevels) {
    level_regridding_epochs.resize(reflevels, 0);
  }

  OutputGridStatistics(cctkGH);

  timer.stop();
}

bool Recompose(cGH const *const cctkGH, int const rl, bool const do_init) {
  char const *const where = "Recompose";
  static Timers::Timer timer(where);
  timer.start();

  bool did_recompose = false;

  for (int m = 0; m < maps; ++m) {
    Waypoint("Recomposing the grid hierarchy for map %d level %d...", m, rl);

    assert(rl >= 0 and rl < vhh.AT(m)->reflevels());
    did_recompose |= vhh.AT(m)->recompose(rl, do_init);

    Waypoint("Done recomposing the grid hierarchy for map %d level %d.", m, rl);
  }

  ClassifyPoints(cctkGH, rl);

  if (did_recompose) {
    ++level_regridding_epochs.at(rl);
  }

  timer.stop();
  return did_recompose;
}

void RegridFree(cGH const *const cctkGH, bool const do_init) {
  char const *const where = "RegridFree";
  static Timers::Timer timer(where);
  timer.start();

  Checkpoint("Freeing after regridding level %d...", reflevel);

  assert(is_level_mode());

  // Free old grid structure
  for (int m = 0; m < maps; ++m) {
    vhh.AT(m)->regrid_free(do_init);
  }

  timer.stop();
}

void OutputSuperregions(cGH const *const cctkGH, int const m, gh const &hh,
                        gh::rregs const &superregss) {
  CCTK_INFO("Grid structure (superregions, grid points):");
  for (int rl = 0; rl < (int)superregss.size(); ++rl) {
    assert(rl < hh.reflevels());
    for (int c = 0; c < (int)superregss.AT(rl).size(); ++c) {
      const ivect &levfact = spacereffacts.AT(rl);
      const ibbox &ext = superregss.AT(rl).AT(c).extent;
      const ivect &lower = ext.lower();
      const ivect &upper = ext.upper();
      const ivect &stride = ext.stride();
      assert(all(lower * levfact % maxspacereflevelfact == 0));
      assert(all(upper * levfact % maxspacereflevelfact == 0));
      cout << "   [" << rl << "][" << m << "][" << c << "]"
           << "   exterior: " << lower / stride << " : " << upper / stride
           << "   (" << (upper - lower) / stride + 1 << " + PADDING) "
           << prod((upper - lower) / stride + 1) << eol;
    }
  }

  CCTK_INFO("Grid structure (superregions, coordinates):");
  for (int rl = 0; rl < (int)superregss.size(); ++rl) {
    assert(rl < hh.reflevels());
    for (int c = 0; c < (int)superregss.AT(rl).size(); ++c) {
      const rvect origin = domainspecs.AT(m).exterior_min;
      const rvect delta =
          (domainspecs.AT(m).exterior_max - domainspecs.AT(m).exterior_min) /
          rvect(domainspecs.AT(m).npoints - 1);
      const ibbox &ext = superregss.AT(rl).AT(c).extent;
      const ivect &ilower = ext.lower();
      const ivect &iupper = ext.upper();
      const ivect &levfact = spacereffacts.AT(rl);
      const ibbox &base = hh.baseextent(0, 0);
      const ivect &bstride = base.stride();
      cout << "   [" << rl << "][" << m << "][" << c << "]"
           << "   exterior: " << origin + delta * rvect(ilower) / rvect(bstride)
           << " : " << origin + delta * rvect(iupper) / rvect(bstride) << " : "
           << delta / rvect(levfact) << eol;
    }
  }

  fflush(stdout);
}

void OutputGrids(cGH const *const cctkGH, int const m, gh const &hh,
                 dh const &dd) {
  DECLARE_CCTK_PARAMETERS;

  if (verbose or veryverbose) {

    CCTK_INFO("Grid structure (grid points):");
    for (int ml = 0; ml < hh.mglevels(); ++ml) {
      for (int rl = 0; rl < hh.reflevels(); ++rl) {
        for (int c = 0; c < hh.components(rl); ++c) {
          const int convfact = ipow(mgfact, ml);
          const ivect levfact = spacereffacts.AT(rl);
          const ibbox ext = hh.extent(ml, rl, c);
          const ivect &lower = ext.lower();
          const ivect &upper = ext.upper();
          const ivect &stride = ext.stride();
          assert(all(lower * levfact % maxspacereflevelfact == 0));
          assert(all(upper * levfact % maxspacereflevelfact == 0));
          assert(all(((upper - lower) * levfact / maxspacereflevelfact) %
                         convfact ==
                     0));
          cout << "   [" << ml << "][" << rl << "][" << m << "][" << c << "]"
               << "   exterior: "
               << "proc " << hh.processor(rl, c) << "   " << lower / stride
               << " : " << upper / stride << "   ("
               << (upper - lower) / stride + 1 << ") "
               << prod((upper - lower) / stride + 1) << eol;
        }
      }
    }

    CCTK_INFO("Grid structure (boundaries):");
    for (int rl = 0; rl < hh.reflevels(); ++rl) {
      for (int c = 0; c < hh.components(rl); ++c) {
        cout << "   [" << rl << "][" << m << "][" << c << "]"
             << "   bbox: " << hh.outer_boundaries(rl, c) << eol;
      }
    }

    CCTK_INFO("Grid structure (coordinates):");
    for (int ml = 0; ml < hh.mglevels(); ++ml) {
      for (int rl = 0; rl < hh.reflevels(); ++rl) {
        for (int c = 0; c < hh.components(rl); ++c) {
          const rvect origin = domainspecs.AT(m).exterior_min;
          const rvect delta = (domainspecs.AT(m).exterior_max -
                               domainspecs.AT(m).exterior_min) /
                              rvect(domainspecs.AT(m).npoints - 1);
          const ibbox ext = hh.extent(ml, rl, c);
          const ivect &ilower = ext.lower();
          const ivect &iupper = ext.upper();
          const int convfact = ipow(mgfact, ml);
          const ivect levfact = spacereffacts.AT(rl);
          const ibbox &base = hh.baseextent(ml, 0);
          const ivect &bstride = base.stride();
          cout << "   [" << ml << "][" << rl << "][" << m << "][" << c << "]"
               << "   exterior: "
               << origin + delta * rvect(ilower) / rvect(bstride) << " : "
               << origin + delta * rvect(iupper) / rvect(bstride) << " : "
               << delta * rvect(convfact) / rvect(levfact) << eol;
        }
      }
    }

    CCTK_INFO("Grid structure (coordinates, including ghosts):");
    for (int ml = 0; ml < hh.mglevels(); ++ml) {
      for (int rl = 0; rl < hh.reflevels(); ++rl) {
        for (int c = 0; c < hh.components(rl); ++c) {
          const rvect origin = domainspecs.AT(m).exterior_min;
          const rvect delta = (domainspecs.AT(m).exterior_max -
                               domainspecs.AT(m).exterior_min) /
                              rvect(domainspecs.AT(m).npoints - 1);
          const ivect lower =
              dd.light_boxes.AT(ml).AT(rl).AT(c).exterior.lower();
          const ivect upper =
              dd.light_boxes.AT(ml).AT(rl).AT(c).exterior.upper();
          const int convfact = ipow(mgfact, ml);
          const ivect levfact = spacereffacts.AT(rl);
          const ibbox &base = hh.baseextent(ml, 0);
          const ivect &bstride = base.stride();
          cout << "   [" << ml << "][" << rl << "][" << m << "][" << c << "]"
               << "   exterior: "
               << origin + delta * rvect(lower) / rvect(bstride) << " : "
               << origin + delta * rvect(upper) / rvect(bstride) << " : "
               << delta * rvect(convfact) / rvect(levfact) << eol;
        }
      }
    }

    CCTK_INFO("Grid statistics:");
    const streamsize oldprecision = cout.precision();
    const ios_base::fmtflags oldflags = cout.flags();
    cout.setf(ios::fixed);
    for (int ml = 0; ml < hh.mglevels(); ++ml) {
      CCTK_REAL CCTK_ATTRIBUTE_UNUSED coarsevolume = 0;
      for (int rl = 0; rl < hh.reflevels(); ++rl) {

        const CCTK_REAL basevolume = hh.baseextents.AT(0).AT(0).size();
        CCTK_REAL countvolume = 0;
        CCTK_REAL totalvolume = 0;
        CCTK_REAL totalvolume2 = 0;

        for (int c = 0; c < hh.components(rl); ++c) {
          const CCTK_REAL volume = hh.extent(ml, rl, c).size();
          ++countvolume;
          totalvolume += volume;
          totalvolume2 += ipow(volume, 2);
        }

        const CCTK_REAL avgvolume = totalvolume / countvolume;
        const CCTK_REAL stddevvolume = sqrt(
            max(CCTK_REAL(0), totalvolume2 / countvolume - ipow(avgvolume, 2)));

        for (int c = 0; c < hh.components(rl); ++c) {
          const CCTK_REAL volume = hh.extent(ml, rl, c).size();
          cout << "   [" << ml << "][" << rl << "][" << m << "][" << c << "]"
               << "   volume: " << setprecision(0) << volume
               << "   of parent: " << setprecision(1)
               << 100 * volume / totalvolume << "%"
               << "   of domain: " << setprecision(1)
               << 100 * volume / basevolume << "%" << eol;
        }

        cout << "   [" << ml << "][" << rl << "][" << m << "]"
             << "   average volume: " << setprecision(0) << avgvolume
             << "   of parent: " << setprecision(1)
             << 100 * avgvolume / totalvolume << "%"
             << "   of domain: " << setprecision(1)
             << 100 * avgvolume / basevolume << "%" << eol;
        cout << "   [" << ml << "][" << rl << "][" << m << "]"
             << "   standard deviation: " << setprecision(0) << stddevvolume
             << "   of parent: " << setprecision(1)
             << 100 * stddevvolume / totalvolume << "%"
             << "   of domain: " << setprecision(1)
             << 100 * stddevvolume / basevolume << "%" << eol;

        // TODO: ghost points vs. interior points (and boundary
        // points)

        CCTK_REAL countquadrupole = 0;
        CCTK_REAL minquadrupole = 1;
        CCTK_REAL totalquadrupole = 0;
        CCTK_REAL totalquadrupole2 = 0;

        for (int c = 0; c < hh.components(rl); ++c) {
          const ibbox ext = hh.extent(ml, rl, c);
          const ivect shape = ext.shape();
          const ivect stride = ext.stride();
          const CCTK_REAL minlength = minval(rvect(shape) / rvect(stride));
          const CCTK_REAL maxlength = maxval(rvect(shape) / rvect(stride));
          const CCTK_REAL quadrupole = minlength / maxlength;
          ++countquadrupole;
          minquadrupole = min(minquadrupole, quadrupole);
          totalquadrupole += quadrupole;
          totalquadrupole2 += ipow(quadrupole, 2);
          cout << "   [" << ml << "][" << rl << "][" << m << "][" << c << "]"
               << "   length ratio: " << setprecision(3) << quadrupole << eol;
        }

        const CCTK_REAL avgquadrupole = totalquadrupole / countquadrupole;
        const CCTK_REAL stddevquadrupole =
            sqrt(max(CCTK_REAL(0), (totalquadrupole2 / countquadrupole -
                                    ipow(avgquadrupole, 2))));

        cout << "   [" << ml << "][" << rl << "][" << m << "]"
             << "   average length ratio: " << setprecision(3) << avgquadrupole
             << "   standard deviation: " << setprecision(3) << stddevquadrupole
             << eol;

        // TODO: processor distribution, average load, std deviation

        coarsevolume = totalvolume * prod(rvect(spacereflevelfact));
        // TODO: output coarsevolume
      }
    }
    cout.precision(oldprecision);
    cout.setf(oldflags);

    fflush(stdout);
  }
}

void OutputGridStructure(cGH const *const cctkGH, int const m,
                         gh::mregs const &regsss) {
  DECLARE_CCTK_PARAMETERS;

  // Output only on the root processor
  if (CCTK_MyProc(cctkGH) != 0)
    return;

  // Output only if output is desired
  if (CCTK_EQUALS(grid_structure_filename, ""))
    return;

  // Create the output directory
  CCTK_CreateDirectory(0755, out_dir);

  ostringstream filenamebuf;
  filenamebuf << out_dir << "/" << grid_structure_filename;
  // we need a persistent temporary here
  string filenamestr = filenamebuf.str();
  const char *filename = filenamestr.c_str();

  ofstream file;

  static bool do_truncate = true;

  if (do_truncate) {
    do_truncate = false;
    struct stat fileinfo;
    if (IO_TruncateOutputFiles(cctkGH) or stat(filename, &fileinfo) != 0) {
      file.open(filename, ios::out | ios::trunc);
      assert(file.good());
      file << "# grid structure" << eol
           << "# format: map reflevel component mglevel   processor "
              "bounding-box is-outer-boundary"
           << eol;
      assert(file.good());
    }
  }
  if (not file.is_open()) {
    file.open(filename, ios::app);
    assert(file.good());
  }

  if (m == 0) {
    // Output a separator
    file << eol;
    file << "iteration " << cctkGH->cctk_iteration << eol;
    file << "maps " << maps << eol;
  }
  file << m << " mglevels " << regsss.size() << eol;
  for (int ml = 0; ml < (int)regsss.size(); ++ml) {
    file << m << " " << ml << " reflevels " << regsss.AT(ml).size() << eol;
    for (int rl = 0; rl < (int)regsss.AT(ml).size(); ++rl) {
      file << m << " " << ml << " " << rl << " components "
           << regsss.AT(ml).AT(rl).size() << eol;
      for (int c = 0; c < (int)regsss.AT(ml).AT(rl).size(); ++c) {
        file << m << " " << ml << " " << rl << " " << c << "   "
             << regsss.AT(ml).AT(rl).AT(c).processor << " "
             << regsss.AT(ml).AT(rl).AT(c).extent << " "
             << regsss.AT(ml).AT(rl).AT(c).outer_boundaries << eol;
      }
    }
  }

  file.close();
  assert(file.good());
}

void OutputGridCoordinates(cGH const *const cctkGH, int const m,
                           gh::mregs const &regsss) {
  DECLARE_CCTK_PARAMETERS;

  // Output only on the root processor
  if (CCTK_MyProc(cctkGH) != 0)
    return;

  // Output only if output is desired
  if (CCTK_EQUALS(grid_coordinates_filename, ""))
    return;

  // Create the output directory
  CCTK_CreateDirectory(0755, out_dir);

  ostringstream filenamebuf;
  filenamebuf << out_dir << "/" << grid_coordinates_filename;
  // we need a persistent temporary here
  string filenamestr = filenamebuf.str();
  const char *filename = filenamestr.c_str();

  ofstream file;

  static bool do_truncate = true;

  if (do_truncate) {
    do_truncate = false;
    struct stat fileinfo;
    if (IO_TruncateOutputFiles(cctkGH) or stat(filename, &fileinfo) != 0) {
      file.open(filename, ios::out | ios::trunc);
      assert(file.good());
      file << "# grid coordinates" << eol
           << "# format: map reflevel region mglevel   bounding-box" << eol;
      assert(file.good());
    }
  }
  if (not file.is_open()) {
    file.open(filename, ios::app);
    assert(file.good());
  }

  // Affine transformation between index space and coordinate space
  rvect const origin = domainspecs.AT(m).exterior_min;
  rvect const spacing =
      (domainspecs.AT(m).exterior_max - domainspecs.AT(m).exterior_min) /
      rvect(domainspecs.AT(m).npoints - 1);
  rvect const scale =
      rvect(vhh.AT(m)->baseextents.AT(0).AT(0).stride()) / spacing;

  file << "iteration " << cctkGH->cctk_iteration << eol;
  file << "maps " << maps << eol;
  file << m << " mglevels " << regsss.size() << eol;
  for (int ml = 0; ml < (int)regsss.size(); ++ml) {
    file << m << " " << ml << " reflevels " << regsss.AT(ml).size() << eol;
    for (int rl = 0; rl < (int)regsss.AT(ml).size(); ++rl) {
      ibset const extents(regsss.AT(ml).AT(rl), &region_t::extent);
      file << m << " " << ml << " " << rl << " regions " << extents.setsize()
           << eol;
      int c = 0;
      for (ibset::const_iterator bi = extents.begin(); bi != extents.end();
           ++bi) {
#if 0
          ibbox const & ext       = * bi;
          ibbox const & baseext   = vhh.AT(m)->baseextents.AT(ml).AT(rl);
          ibbox const & coarseext = vhh.AT(m)->baseextents.AT(ml).AT(0 );
          
          // This is nice, but is wrong since CartGrid3D has not yet
          // initialised the coordinates
          ivect const cctk_levfac       = spacereffacts.AT(rl);
          ivect const cctk_levoff       = baseext.lower() - coarseext.lower();
          ivect const cctk_levoffdenom  = baseext.stride();
          
          ivect const cctk_lbnd =
            (ext.lower() - baseext.lower()) / ext.stride();
          ivect const cctk_ubnd =
            (ext.upper() - baseext.lower()) / ext.stride();
          
          rvect const cctk_origin_space =
            origin_space.AT(m).AT(ml);
          rvect const cctk_delta_space  =
            delta_space.AT(m) * rvect (ipow (mgfact, ml));
          
          rvect const CCTK_ORIGIN_SPACE =
            cctk_origin_space +
            cctk_delta_space / rvect (cctk_levfac) *
            rvect (cctk_levoff) / rvect (cctk_levoffdenom);
          rvect const CCTK_DELTA_SPACE =
            cctk_delta_space / rvect (cctk_levfac);
          
          rvect const rlower =
            CCTK_ORIGIN_SPACE + rvect (cctk_lbnd) * CCTK_DELTA_SPACE;
          rvect const rupper =
            CCTK_ORIGIN_SPACE + rvect (cctk_ubnd) * CCTK_DELTA_SPACE;
          rvect const rdelta =
            CCTK_DELTA_SPACE;
#endif

        ibbox const &ext = *bi;

        rvect const rlower = origin + rvect(ext.lower()) / scale;
        rvect const rupper = origin + rvect(ext.upper()) / scale;
        rvect const rdelta = rvect(ext.stride()) / scale;

        rbbox const rb(rlower, rupper, rdelta);
        const std::streamsize oldprecision = file.precision(19);
        file << m << " " << ml << " " << rl << " " << c << "   " << rb << eol;
        file.precision(oldprecision);
        ++c;
      }
    }
  }
  file << eol;

  file.close();
  assert(file.good());
}

// TODO: this routine should go into CarpetRegrid (except maybe
// SplitRegions_AlongZ for grid arrays)
void SplitRegions(cGH const *const cctkGH, vector<region_t> &superregs,
                  vector<region_t> &regs) {
  DECLARE_CCTK_PARAMETERS;

  assert(regs.empty());
  if (CCTK_EQUALS(processor_topology, "along-z")) {
    SplitRegions_AlongZ(cctkGH, superregs, regs);
  } else if (CCTK_EQUALS(processor_topology, "along-dir")) {
    SplitRegions_AlongDir(cctkGH, superregs, regs, split_direction);
  } else if (CCTK_EQUALS(processor_topology, "automatic")) {
    SplitRegions_Automatic(cctkGH, superregs, regs);
  } else if (CCTK_EQUALS(processor_topology, "manual")) {
    SplitRegions_AsSpecified(cctkGH, superregs, regs);
  } else {
    assert(0);
  }
}

void OutputGridStatistics(cGH const *const cctkGH) {
  // Grid array statistics
  int num_arrays = 0;
  int num_gfs = 0;
  int size_gfs = 0;
  CCTK_REAL num_active_array_points = 0;
  CCTK_REAL num_total_array_points = 0;
  CCTK_REAL size_total_array_points = 0;
  for (int g = 0; g < CCTK_NumGroups(); ++g) {
    cGroup gdata;
    check(not CCTK_GroupData(g, &gdata));
    int const num_tl = CCTK_ActiveTimeLevelsGI(cctkGH, g);
    int const num_vars = gdata.numvars;
    int const size_vars = gdata.numvars * CCTK_VarTypeSize(gdata.vartype);
    switch (gdata.grouptype) {
    case CCTK_SCALAR:
    case CCTK_ARRAY: {
      num_arrays += num_tl * num_vars;
      gh const *const hh = arrdata.AT(g).AT(0).hh;
      dh const *const dd = arrdata.AT(g).AT(0).dd;
      for (int c = 0; c < hh->components(0); ++c) {
        dh::light_dboxes const &b = dd->light_boxes.AT(0).AT(0).AT(c);
        num_active_array_points += num_tl * num_vars * b.active_size;
        num_total_array_points += num_tl * num_vars * b.exterior_size;
        size_total_array_points += num_tl * size_vars * b.exterior_size;
      }
      break;
    }
    case CCTK_GF:
      num_gfs += num_tl * num_vars;
      size_gfs += num_tl * size_vars;
      break;
    default:
      assert(0);
    }
  }

  // Grid function statistics
  int num_comps = 0;
  int num_steps = 0;
  CCTK_REAL num_active_mem_points = 0;
  CCTK_REAL num_owned_mem_points = 0;
  CCTK_REAL num_total_mem_points = 0;
  CCTK_REAL size_total_mem_points = 0;
  CCTK_REAL num_active_cpu_points = 0;
  CCTK_REAL num_owned_cpu_points = 0;
  CCTK_REAL num_total_cpu_points = 0;
  for (int m = 0; m < maps; ++m) {
    gh const *const hh = vhh.AT(m);
    dh const *const dd = vdd.AT(m);
    for (int ml = 0; ml < mglevels; ++ml) {
      for (int rl = 0; rl < reflevels; ++rl) {
        int const trf = timereffacts.AT(rl);
        if (m == 0 and ml == 0) {
          num_steps += trf;
        }
        for (int c = 0; c < hh->components(rl); ++c) {
          ++num_comps;
          dh::light_dboxes const &b = dd->light_boxes.AT(ml).AT(rl).AT(c);
          num_active_mem_points += num_gfs * b.active_size;
          num_owned_mem_points += num_gfs * b.owned_size;
          num_total_mem_points += num_gfs * b.exterior_size;
          size_total_mem_points += size_gfs * b.exterior_size;
          num_active_cpu_points += trf * b.active_size;
          num_owned_cpu_points += trf * b.owned_size;
          num_total_cpu_points += trf * b.exterior_size;
        }
      }
    }
  }
  num_active_cpu_points /= num_steps * delta_time;
  num_owned_cpu_points /= num_steps * delta_time;
  num_total_cpu_points /= num_steps * delta_time;

  // Load balance statistics
  // We count owned points, since these are evolved. Owned points
  // include active points, overlap points, and buffer points.
  CCTK_REAL const rzero = 0;
  CCTK_REAL const rmax = numeric_limits<CCTK_REAL>::max();
  vector<CCTK_REAL> min_owned(reflevels, rmax);
  vector<CCTK_REAL> max_owned(reflevels, 0);
  vector<CCTK_REAL> avg_owned(reflevels, 0);
  vector<CCTK_REAL> sdv_owned(reflevels, 0);
  for (int rl = 0; rl < reflevels; ++rl) {
    vector<CCTK_REAL> num_owned_per_proc(dist::size(), 0);
    for (int m = 0; m < maps; ++m) {
      gh const *const hh = vhh.AT(m);
      dh const *const dd = vdd.AT(m);
      for (int ml = 0; ml < mglevels; ++ml) {
        for (int c = 0; c < hh->components(rl); ++c) {
          int const p = hh->processor(rl, c);
          dh::light_dboxes const &b = dd->light_boxes.AT(ml).AT(rl).AT(c);
          num_owned_per_proc.AT(p) += num_gfs * b.owned_size;
        }
      }
    }
    for (int p = 0; p < dist::size(); ++p) {
      min_owned.AT(rl) = min(min_owned.AT(rl), num_owned_per_proc.AT(p));
      max_owned.AT(rl) = max(max_owned.AT(rl), num_owned_per_proc.AT(p));
      avg_owned.AT(rl) += num_owned_per_proc.AT(p);
      sdv_owned.AT(rl) += ipow(num_owned_per_proc.AT(p), 2);
    }
    avg_owned.AT(rl) /= dist::size();
    sdv_owned.AT(rl) /= dist::size();
    sdv_owned.AT(rl) =
        sqrt(max(rzero, sdv_owned.AT(rl) - ipow(avg_owned.AT(rl), 2)));
  } // for rl

  // Output
  CCTK_VInfo(CCTK_THORNSTRING, "Global grid structure statistics:");
  CCTK_VInfo(
      CCTK_THORNSTRING, "GF: rhs: %.0fk active, %.0fk owned (+%.0f%%), %.0fk "
                        "total (+%.0f%%), %.3g steps/time",
      double(num_active_cpu_points / 1e+3), double(num_owned_cpu_points / 1e+3),
      double(num_active_cpu_points == 0
                 ? 0
                 : num_owned_cpu_points / num_active_cpu_points * 100 - 100),
      double(num_total_cpu_points / 1e+3),
      double(num_owned_cpu_points == 0
                 ? 0
                 : num_total_cpu_points / num_owned_cpu_points * 100 - 100),
      double(num_steps / delta_time));
  CCTK_VInfo(
      CCTK_THORNSTRING, "GF: vars: %d, pts: %.0fM active, %.0fM owned "
                        "(+%.0f%%), %.0fM total (+%.0f%%), %.1f comp/proc",
      num_gfs, double(num_active_mem_points / 1e+6),
      double(num_owned_mem_points / 1e+6),
      double(num_active_mem_points == 0
                 ? 0
                 : num_owned_mem_points / num_active_mem_points * 100 - 100),
      double(num_total_mem_points / 1e+6),
      double(num_owned_mem_points == 0
                 ? 0
                 : num_total_mem_points / num_owned_mem_points * 100 - 100),
      double(num_comps / (reflevels * CCTK_nProcs(cctkGH))));
  CCTK_VInfo(
      CCTK_THORNSTRING,
      "GA: vars: %d, pts: %.0fM active, %.0fM total (+%.0f%%)", num_arrays,
      double(num_active_array_points / 1e+6),
      double(num_total_array_points / 1e+6),
      double(num_active_array_points == 0
                 ? 0
                 : num_total_array_points / num_active_array_points * 100 -
                       100));
  CCTK_VInfo(
      CCTK_THORNSTRING,
      "Total required memory: %.3f GByte (for GAs and currently active GFs)",
      double((size_total_array_points + size_total_mem_points) / 1e+9));
  CCTK_VInfo(CCTK_THORNSTRING,
             "Load balance:  min     avg     max     sdv     max/avg-1");
  for (int rl = 0; rl < reflevels; ++rl) {
    CCTK_VInfo(
        CCTK_THORNSTRING,
        "Level %2d:    %4.0fM   %4.0fM   %4.0fM   %4.0fM owned   %4.0f%%", rl,
        double(min_owned.AT(rl) / 1e+6), double(avg_owned.AT(rl) / 1e+6),
        double(max_owned.AT(rl) / 1e+6), double(sdv_owned.AT(rl) / 1e+6),
        double(100 * (max_owned.AT(rl) / avg_owned.AT(rl) - 1)));
  }

  // After this, we will begin to allocate memory for the grid
  // structure.  If we run out of memory, ensure that this output
  // still makes it to disk.
  fflush(stdout);
  fflush(stderr);
}

void SplitRegions_AlongZ(cGH const *const cctkGH, vector<region_t> &superregs,
                         vector<region_t> &regs) {
  assert(regs.empty());
  SplitRegions_AlongDir(cctkGH, superregs, regs, dim - 1);
}

void SplitRegions_AlongDir(cGH const *const cctkGH, vector<region_t> &superregs,
                           vector<region_t> &regs, int const dir) {
  assert(regs.empty());

  // Something to do?
  if (superregs.size() == 0) {
    return;
  }

  const int nprocs = CCTK_nProcs(cctkGH);

  assert(superregs.size() == 1);
  regs = superregs;

  if (nprocs == 1) {
    regs.AT(0).processor = 0;
    pseudoregion_t const preg(regs.AT(0).extent, regs.AT(0).processor);
    assert(superregs.AT(0).processors == NULL);
    superregs.AT(0).processors = new ipfulltree(preg);
    return;
  }

  assert(dir >= 0 and dir < dim);

  region_t const &reg0 = regs.AT(0);
  const ivect rstr0 = reg0.extent.stride();
  const ivect rlb0 = reg0.extent.lower();
  const ivect rub0 = reg0.extent.upper() + rstr0;
  const b2vect obnd0 = reg0.outer_boundaries;

  regs.resize(nprocs);
  vector<int> bounds(nprocs + 1);
  vector<ipfulltree *> subtrees(nprocs);
  for (int c = 0; c < nprocs; ++c) {
    ivect cstr = rstr0;
    ivect clb = rlb0;
    ivect cub = rub0;
    const int glonpz = (rub0[dir] - rlb0[dir]) / cstr[dir];
    const int locnpz = (glonpz + nprocs - 1) / nprocs;
    const int zstep = locnpz * cstr[dir];
    clb[dir] = rlb0[dir] + zstep * c;
    cub[dir] = rlb0[dir] + zstep * (c + 1);
    if (clb[dir] > rub0[dir])
      clb[dir] = rub0[dir];
    if (cub[dir] > rub0[dir])
      cub[dir] = rub0[dir];
    assert(clb[dir] <= cub[dir]);
    assert(cub[dir] <= rub0[dir]);
    region_t &reg = regs.AT(c);
    ibbox &ext = reg.extent;
    b2vect &obnd = reg.outer_boundaries;
    int &proc = reg.processor;
    ext = ibbox(clb, cub - cstr, cstr);
    obnd = obnd0;
    obnd[0][dir] &= clb[dir] == rlb0[dir];
    obnd[1][dir] &= cub[dir] == rub0[dir];
    proc = c;
    pseudoregion_t const preg(reg.extent, c);
    bounds.AT(c) = reg.extent.lower()[dir];
    subtrees.AT(c) = new ipfulltree(preg);
  }
  bounds.AT(nprocs) = rub0[dir] + rstr0[dir];

  assert(superregs.AT(0).processors == NULL);
  superregs.AT(0).processors = new ipfulltree(dir, bounds, subtrees);

  for (int c = 0; c < (int)regs.size(); ++c) {
    assert(regs.AT(c).processor == c);
  }
}

void SplitRegions_Automatic(cGH const *const cctkGH,
                            vector<region_t> &superregs, vector<region_t> &regs,
                            bvect const &no_split_dims) {
  assert(regs.empty());
  vector<vector<region_t> > superregss(1, superregs);
  vector<vector<region_t> > regss(1);
  SplitRegionsMaps_Automatic(cctkGH, superregss, regss, no_split_dims);
  assert(superregss.size() == 1);
  superregs = superregss.AT(0);
  assert(regss.size() == 1);
  regs = regss.AT(0);
}

static void SplitRegions_AsSpecified(cGH const *const cctkGH,
                                     vector<region_t> &superregs,
                                     vector<region_t> &regs) {
  DECLARE_CCTK_PARAMETERS;

  assert(regs.empty());

  // Something to do?
  if (superregs.size() == 0) {
    return;
  }

  const int nprocs = CCTK_nProcs(cctkGH);

  assert(superregs.size() == 1);

  region_t const &reg0 = superregs.AT(0);
  const ivect rstr0 = reg0.extent.stride();
  const ivect rlb0 = reg0.extent.lower();
  const ivect rub0 = reg0.extent.upper() + rstr0;
  const b2vect obnd0 = reg0.outer_boundaries;

  const ivect nprocs_dir(processor_topology_3d_x, processor_topology_3d_y,
                         processor_topology_3d_z);
  assert(all(nprocs_dir > 0));
  if (prod(nprocs_dir) != nprocs) {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "The specified processor topology [%d,%d,%d] requires %d "
               "processors, but there are %d processors",
               nprocs_dir[0], nprocs_dir[1], nprocs_dir[2], prod(nprocs_dir),
               nprocs);
  }
  assert(prod(nprocs_dir) == nprocs);

  regs.resize(nprocs);
  const ivect cstr = rstr0;
  const ivect glonp = (rub0 - rlb0) / cstr;
  //     const ivect locnp = (glonp + nprocs_dir - 1) / nprocs_dir;
  const ivect locnp = glonp / nprocs_dir;
  const ivect rem = glonp % nprocs_dir;
  const ivect step = locnp * cstr;
  assert(dim == 3);

  vector<int> boundsz(nprocs_dir[2] + 1);
  vector<ipfulltree *> subtreesz(nprocs_dir[2]);
  for (int k = 0; k < nprocs_dir[2]; ++k) {
    vector<int> boundsy(nprocs_dir[1] + 1);
    vector<ipfulltree *> subtreesy(nprocs_dir[1]);
    for (int j = 0; j < nprocs_dir[1]; ++j) {
      vector<int> boundsx(nprocs_dir[0] + 1);
      vector<ipfulltree *> subtreesx(nprocs_dir[0]);
      for (int i = 0; i < nprocs_dir[0]; ++i) {

        const int c = i + nprocs_dir[0] * (j + nprocs_dir[1] * k);
        const ivect ipos(i, j, k);
        ivect clb = rlb0 + step * ipos;
        ivect cub = rlb0 + step * (ipos + 1);
        // 	  clb = min (clb, rub);
        // 	  cub = min (cub, rub);
        for (int d = 0; d < dim; ++d) {
          if (ipos[d] < rem[d]) {
            clb[d] += cstr[d] * ipos[d];
            cub[d] += cstr[d] * (ipos[d] + 1);
          } else {
            clb[d] += cstr[d] * rem[d];
            cub[d] += cstr[d] * rem[d];
          }
        }
        assert(all(clb >= 0));
        assert(all(clb <= cub));
        assert(all(cub <= rub0));
        assert(all(not(ipos == 0) or clb == rlb0));
        assert(all(not(ipos == nprocs_dir - 1) or cub == rub0));
        region_t &reg = regs.AT(c);
        ibbox &ext = reg.extent;
        b2vect &obnd = reg.outer_boundaries;
        int &proc = reg.processor;
        ext = ibbox(clb, cub - cstr, cstr);
        obnd = obnd0;
        obnd[0] &= clb == rlb0;
        obnd[1] &= cub == rub0;
        proc = c;

        pseudoregion_t preg(reg.extent, c);
        subtreesx.AT(i) = new ipfulltree(preg);
      }
      boundsx.AT(nprocs_dir[0]) = rub0[0] + rstr0[0];
      subtreesy.AT(j) = new ipfulltree(0, boundsx, subtreesx);
    }
    boundsy.AT(nprocs_dir[1]) = rub0[1] + rstr0[1];
    subtreesz.AT(k) = new ipfulltree(1, boundsy, subtreesy);
  }
  boundsz.AT(nprocs_dir[2]) = rub0[2] + rstr0[2];

  assert(superregs.AT(0).processors == NULL);
  superregs.AT(0).processors = new ipfulltree(2, boundsz, subtreesz);

  for (int c = 0; c < (int)regs.size(); ++c) {
    assert(regs.AT(c).processor == c);
  }
}

// TODO: this routine should go into CarpetRegrid (except maybe
// SplitRegions_AlongZ for grid arrays)
void SplitRegionsMaps(cGH const *const cctkGH,
                      vector<vector<region_t> > &superregss,
                      vector<vector<region_t> > &regss) {
  DECLARE_CCTK_PARAMETERS;

  assert(regss.size() == superregss.size());
  for (size_t m = 0; m < regss.size(); ++m) {
    assert(regss.AT(m).empty());
  }
  if (CCTK_EQUALS(processor_topology, "along-z")) {
    assert(0);
    //       SplitRegionsMaps_AlongZ (cctkGH, superregss, regss);
  } else if (CCTK_EQUALS(processor_topology, "along-dir")) {
    assert(0);
    //       SplitRegionsMaps_AlongDir (cctkGH, superregss, regss,
    //       split_direction);
  } else if (CCTK_EQUALS(processor_topology, "automatic")) {
    SplitRegionsMaps_Automatic(cctkGH, superregss, regss);
  } else if (CCTK_EQUALS(processor_topology, "recursive")) {
    SplitRegionsMaps_Recursively(cctkGH, superregss, regss);
  } else if (CCTK_EQUALS(processor_topology, "balanced")) {
    SplitRegionsMaps_Balanced(cctkGH, superregss, regss);
  } else if (CCTK_EQUALS(processor_topology, "manual")) {
    SplitRegionsMaps_AsSpecified(cctkGH, superregss, regss);
  } else {
    assert(0);
  }

  for (size_t m = 0; m < superregss.size(); ++m) {
    for (size_t r = 0; r < superregss.AT(m).size(); ++r) {
      bool const good_superregs = superregss.AT(m).AT(r).check_region(true);
      if (not good_superregs) {
        cout << "superregs[" << m << "][" << r << "]:\n"
             << superregss.AT(m).AT(r) << "\n";
        cout << "all superregions:\n" << superregss << "\n";
        CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Superregions failed self-check for map %d superregion %d",
                   int(m), int(r));
      }
    }
  }
  assert(regss.size() == superregss.size());
  for (size_t m = 0; m < regss.size(); ++m) {
    for (size_t r = 0; r < regss.AT(m).size(); ++r) {
      bool const good_regs = regss.AT(m).AT(r).check_region(false);
      if (not good_regs) {
        cout << "regs[" << m << "][" << r << "]:\n" << regss.AT(m).AT(r)
             << "\n";
        cout << "all superregions:\n" << superregss << "\n";
        cout << "all regions:\n" << regss << "\n";
        CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Regions failed self-check for map %d superregion %d",
                   int(m), int(r));
      }
    }
  }
}

static void SplitRegionsMaps_Automatic_Recursively(bvect const &dims,
                                                   int const firstproc,
                                                   int const nprocs,
                                                   region_t &superreg,
                                                   vector<region_t> &newregs) {
  DECLARE_CCTK_PARAMETERS;

  if (recompose_verbose)
    cout << "SRMAR enter" << endl;
  // Check preconditions
  assert(nprocs >= 1);

  // Remember old size
  size_t const oldsize = newregs.size();

  // Are we done?
  if (all(dims)) {
    if (recompose_verbose)
      cout << "SRMAR bottom" << endl;

    // Check precondition
    assert(nprocs == 1);

    // Return argument
    region_t newreg = superreg;
    newreg.processor = firstproc;
    newregs.push_back(newreg);
    pseudoregion_t const preg(newreg.extent, newreg.processor);
    superreg.processors = new ipfulltree(preg);

    // Check postcondition
    assert(newregs.size() == oldsize + nprocs);

    if (recompose_verbose)
      cout << "SRMAR exit" << endl;
    return;
  }

  // Special case
  if (superreg.extent.empty()) {
    if (recompose_verbose)
      cout << "SRMAR empty" << endl;

    // Choose a direction
    int mydim = -1;
    for (int d = 0; d < dim; ++d) {
      if (not dims[d]) {
        mydim = d;
        break;
      }
    }
    assert(mydim >= 0 and mydim < dim);

    // Create a new region
    region_t newreg(superreg);
    newreg.outer_boundaries = b2vect(false);
    if (recompose_verbose)
      cout << "SRMAR newreg " << newreg << endl;

#if 0
      // Store
      for (int p=0; p<nprocs; ++p) {
        newreg.processor = firstproc + p;
        newregs.push_back (newreg);
      }
      assert (superreg.processors == NULL);
      superreg.processors = new ipfulltree ();
#endif

    // Store
    vector<int> bounds(nprocs + 1);
    vector<ipfulltree *> subtrees(nprocs);
    for (int p = 0; p < nprocs; ++p) {
      newreg.processor = firstproc + p;
      newregs.push_back(newreg);
      bounds.AT(p) = newreg.extent.lower()[mydim];
      pseudoregion_t const preg(newreg.extent, newreg.processor);
      subtrees.AT(p) = new ipfulltree(preg);
    }
    bounds.AT(nprocs) = newreg.extent.lower()[mydim];
    assert(superreg.processors == NULL);
    superreg.processors = new ipfulltree(mydim, bounds, subtrees);

    // Check postcondition
    assert(newregs.size() == oldsize + nprocs);

    if (recompose_verbose)
      cout << "SRMAR exit" << endl;
    return;
  }

  // Calculate cost factors
  rvect rcost = cost(superreg);
  CCTK_REAL const rfact = pow(nprocs / prod(rcost), CCTK_REAL(1) / dim);
  rcost *= rfact;
  assert(fabs(prod(rcost) - nprocs) <= 1.0e-6);
  if (recompose_verbose)
    cout << "SRMA shapes " << rcost << endl;

  // Choose a direction
  int mydim = -1;
  int nslices = -1;
  if (no_split_direction != -1 and not dims[no_split_direction]) {
    mydim = no_split_direction;
    nslices = 1;
  } else {
    int alldims = 0;
    CCTK_REAL mycost = 0;
    CCTK_REAL totalcost = 1;
    for (int d = 0; d < dim; ++d) {
      if (not dims[d]) {
        ++alldims;
        CCTK_REAL const thiscost = rcost[d];
        // Prefer to split in the z direction
        // (Choose this dimension even if it is slightly worse than
        // a previous dimension)
        if (thiscost >= 0.999999 * mycost) {
          mydim = d;
          mycost = thiscost;
        }
        totalcost *= thiscost;
      }
    }
    assert(mydim >= 0 and mydim < dim);
    assert(mycost >= 0);
    if (recompose_verbose)
      cout << "SRMAR mydim " << mydim << endl;
    if (recompose_verbose)
      cout << "SRMAR mycost " << mycost << endl;

    // Choose a number of slices for this direction
    CCTK_REAL const mycost1 =
        mycost * pow(nprocs / totalcost, CCTK_REAL(1) / alldims);
    // Round up, so that we use more slices in this direction than
    // in the left-over directions
    nslices = min(nprocs, int(lrint(ceil(mycost1))));
  }
  assert(nslices <= nprocs);
  if (recompose_verbose)
    cout << "SRMAR " << mydim << " nprocs " << nprocs << endl;
  if (recompose_verbose)
    cout << "SRMAR " << mydim << " nslices " << nslices << endl;

  // Mark this direction as done
  assert(not dims[mydim]);
  bvect const newdims = dims.replace(mydim, true);

  // Split the remaining processors
  vector<int> mynprocs(nslices);
  int const mynprocs_base = nprocs / nslices;
  int const mynprocs_left = nprocs - nslices * mynprocs_base;
  int sum_mynprocs = 0;
  for (int n = 0; n < nslices; ++n) {
    mynprocs.AT(n) = mynprocs_base + int(n < mynprocs_left);
    sum_mynprocs += mynprocs.AT(n);
  }
  assert(sum_mynprocs == nprocs);
  if (recompose_verbose)
    cout << "SRMAR " << mydim << " mynprocs " << mynprocs << endl;

  // Split the region
  vector<int> mynpoints(nslices);
  int const npoints =
      (superreg.extent.shape() / superreg.extent.stride())[mydim];
  int const npoints_bnd_lo =
      superreg.outer_boundaries[0][mydim] ? granularity_boundary : 0;
  int const npoints_bnd_hi =
      superreg.outer_boundaries[1][mydim] ? granularity_boundary : 0;

  // Keep track of how many points and processors we have left to
  // distribute
  int npoints_left = npoints;
  int nprocs_left = nprocs;
  for (int n = 0; n < nslices; ++n) {
    assert(nprocs_left > 0);
    CCTK_REAL const npoints1 =
        CCTK_REAL(1) * npoints_left * mynprocs.AT(n) / nprocs_left;
    mynpoints.AT(n) =
        int(floor(npoints1 / granularity + CCTK_REAL(0.5))) * granularity;
    if (n == 0)
      mynpoints.AT(n) += npoints_bnd_lo;
    if (npoints_left - mynpoints.AT(n) <= npoints_bnd_hi)
      mynpoints.AT(n) += npoints_bnd_hi;
    mynpoints.AT(n) = min(npoints_left, mynpoints.AT(n));
    assert(mynpoints.AT(n) >= 0 and mynpoints.AT(n) <= npoints_left);
    npoints_left -= mynpoints.AT(n);
    nprocs_left -= mynprocs.AT(n);
    assert(npoints_left >= 0);
    assert(nprocs_left >= 0);
  }
  assert(npoints_left == 0);
  assert(nprocs_left == 0);
  if (recompose_verbose)
    cout << "SRMAR " << mydim << " mynpoints " << mynpoints << endl;

  // Create the regions and recurse
  if (recompose_verbose)
    cout << "SRMAR " << mydim << ": create bboxes" << endl;
  ivect lo = superreg.extent.lower();
  ivect up = superreg.extent.upper();
  ivect const str = superreg.extent.stride();
  int p = firstproc;
  vector<int> bounds(nslices + 1);
  vector<ipfulltree *> subtrees(nslices);
  for (int n = 0; n < nslices; ++n) {
    if (recompose_verbose)
      cout << "SRMAR " << mydim << " n " << n << endl;

    // Create a new region
    up[mydim] = lo[mydim] + (mynpoints.AT(n) - 1) * str[mydim];
    b2vect newob(superreg.outer_boundaries);
    if (n > 0) {
      newob[0][mydim] = false;
    }
    if (n < nslices - 1) {
      up[mydim] = lo[mydim] + (mynpoints.AT(n) - 1) * str[mydim];
      newob[1][mydim] = false;
    }
    region_t newreg(superreg);
    newreg.extent = ibbox(lo, up, str);
    newreg.outer_boundaries = newob;
    if (recompose_verbose)
      cout << "SRMAR " << mydim << " newreg " << newreg << endl;

    // Recurse
    bounds.AT(n) = lo[mydim];
    SplitRegionsMaps_Automatic_Recursively(newdims, p, mynprocs.AT(n), newreg,
                                           newregs);
    if (recompose_verbose)
      cout << "SRMAR " << mydim << " newregs " << newregs << endl;

    assert(newreg.processors != NULL);
    subtrees.AT(n) = newreg.processors;
    newreg.processors = NULL;
    newreg.processor = p;

    // Next slice
    lo[mydim] = up[mydim] + str[mydim];
    p += mynprocs.AT(n);
  }
  assert(up[mydim] == superreg.extent.upper()[mydim]);
  assert(p == firstproc + nprocs);
  bounds.AT(nslices) = up[mydim] + str[mydim];

  // Create tree
  assert(superreg.processors == NULL);
  superreg.processors = new ipfulltree(mydim, bounds, subtrees);

  // Check postcondition
  assert(newregs.size() == oldsize + nprocs);

  if (recompose_verbose)
    cout << "SRMAR exit" << endl;
}

void SplitRegionsMaps_Automatic(cGH const *const cctkGH,
                                vector<vector<region_t> > &superregss,
                                vector<vector<region_t> > &regss,
                                bvect const &no_split_dims) {
  DECLARE_CCTK_PARAMETERS;

  if (recompose_verbose)
    cout << "SRMA enter" << endl;

  int const nmaps = superregss.size();
  int minmap = 1000000000;
  for (int m = 0; m < nmaps; ++m) {
    for (int r = 0; r < int(superregss.AT(m).size()); ++r) {
      minmap = min(minmap, superregss.AT(m).AT(r).map);
    }
  }
  int nregs = 0;
  for (int m = 0; m < nmaps; ++m) {
    nregs += superregss.AT(m).size();
  }
  if (recompose_verbose)
    cout << "SRMA nregs " << nregs << endl;

  // Something to do?
  if (nregs == 0) {
    return;
  }

  // Collect slices
  vector<region_t> superregs;
  {
    for (int m = 0; m < nmaps; ++m) {
      combine_regions(superregss.AT(m), superregs);
    }
    nregs = superregs.size();

    // If the last region was removed, add a new empty region again.
    // A set of regions (corresponding to a refinement level or a
    // grid array) cannot be empty.
    if (nregs == 0) {
      assert(nmaps == 1); // we should only be here for grid
                          // arrays
      region_t reg;
      reg.extent = ibbox(ivect(0), ivect(-1), ivect(1));
      reg.outer_boundaries = b2vect(bvect(true), bvect(true));
      reg.map = 0;
      superregs.push_back(reg);
      nregs = superregs.size();
    }
  }

  int const real_nprocs = CCTK_nProcs(cctkGH);
  if (recompose_verbose)
    cout << "SRMA real_nprocs " << real_nprocs << endl;

  // Deactivate some processors if there are too many
  int nprocs;
  if (min_points_per_proc == 0) {
    nprocs = real_nprocs;
  } else {
    CCTK_REAL mycost = 0;
    for (int r = 0; r < nregs; ++r) {
      mycost += prod(cost(superregs.AT(r)));
    }
    int const goodnprocs = int(floor(mycost / min_points_per_proc));
    nprocs = max(1, min(real_nprocs, goodnprocs));
  }
  if (recompose_verbose)
    cout << "SRMA nprocs " << nprocs << endl;

  // ncomps: number of components per processor
  int const ncomps = (nregs + nprocs - 1) / nprocs;
  if (recompose_verbose)
    cout << "SRMA ncomps " << ncomps << endl;
  assert(ncomps > 0);

  int const newnregs = nprocs * ncomps;
  if (recompose_verbose)
    cout << "SRMA newnregs " << newnregs << endl;

  if (recompose_verbose)
    cout << "SRMA: distributing processors to regions" << endl;
  vector<CCTK_REAL> mycosts(nregs);
  for (int r = 0; r < nregs; ++r) {
    mycosts.AT(r) = prod(cost(superregs.AT(r)));
  }
  int nregs_left = newnregs;
  vector<int> mynprocs(nregs);
  for (int r = 0; r < nregs; ++r) {
    mynprocs.AT(r) = 1;
    --nregs_left;
  }
  // TODO: split regions if necessary
  while (nregs_left > 0) {
    if (recompose_verbose)
      cout << "SRMA nregs_left " << nregs_left << endl;
    int maxr = -1;
    CCTK_REAL maxratio = -1;
    for (int r = 0; r < nregs; ++r) {
      CCTK_REAL const ratio = mycosts.AT(r) / mynprocs.AT(r);
      if (ratio > maxratio) {
        maxr = r;
        maxratio = ratio;
      }
    }
    assert(maxr >= 0 and maxr < nregs);
    ++mynprocs.AT(maxr);
    if (recompose_verbose)
      cout << "SRMA maxr " << maxr << endl;
    if (recompose_verbose)
      cout << "SRMA mynprocs[maxr] " << mynprocs.AT(maxr) << endl;
    --nregs_left;
  }
  assert(nregs_left == 0);
  int sum_nprocs = 0;
  for (int r = 0; r < nregs; ++r) {
    sum_nprocs += mynprocs.AT(r);
  }
  assert(sum_nprocs == newnregs);
  if (recompose_verbose)
    cout << "SRMA mynprocs " << mynprocs << endl;

  if (recompose_verbose)
    cout << "SRMA: splitting work units" << endl;
  // TODO: rename newregs to regs
  // TODO: rename nregs to nsuperregs
  // TODO: rename newnregs to nregs
  vector<region_t> newregs;
  newregs.reserve(newnregs);
  for (int r = 0, p = 0; r < nregs; p += mynprocs.AT(r), ++r) {
    if (recompose_verbose)
      cout << "SRMA superreg[" << r << "] " << superregs.AT(r) << endl;
    // bvect const dims = false;
    bvect const dims = not split_components || no_split_dims;
    SplitRegionsMaps_Automatic_Recursively(dims, p, mynprocs.AT(r),
                                           superregs.AT(r), newregs);
  } // for r
  if (recompose_verbose)
    cout << "SRMA newregs " << newregs << endl;
  assert(int(newregs.size()) == newnregs);

  // Count components per map
  vector<int> myncomps(nmaps, 0);
#if 0
    vector<int> empty_comps(nmaps, 0);
#endif
  for (int r = 0; r < newnregs; ++r) {
    int const m = newregs.AT(r).map - minmap;
    assert(m >= 0 and m < nmaps);
#if 0
      if (not newregs.AT(r).extent.empty()) {
        // Ignore empty regions, which may come from empty grid arrays
        ++ myncomps.AT(m);
      } else {
        ++ empty_comps.AT(m);
      }
#endif
    ++myncomps.AT(m);
  }
  vector<int> mynregs(nmaps, 0);
#if 0
    vector<int> empty_regs(nmaps, 0);
#endif
  for (int r = 0; r < nregs; ++r) {
    int const m = superregs.AT(r).map - minmap;
    assert(m >= 0 and m < nmaps);
    ++mynregs.AT(m);
#if 0
      if (not superregs.AT(r).extent.empty()) {
        // Ignore empty regions, which may come from empty grid arrays
        ++ mynregs.AT(m);
      } else {
        ++ empty_regs.AT(m);
      }
#endif
  }

  // Convert virtual to real processors
  for (int r = 0; r < newnregs; ++r) {
    newregs.AT(r).processor /= ncomps;
    assert(newregs.AT(r).processor >= 0 and newregs.AT(r).processor < nprocs);
  }
  {
    vector<int> tmpncomps(nmaps, 0);
    for (int r = 0; r < nregs; ++r) {
      int const m = superregs.AT(r).map - minmap;
      ipfulltree *const regf = superregs.AT(r).processors;
      assert(regf != NULL);
      for (ipfulltree::iterator fti(*regf); not fti.done(); ++fti) {
        pseudoregion_t &preg = (*fti).payload();
        preg.component = tmpncomps.AT(m)++;
      }
    }
    for (int m = 0; m < nmaps; ++m) {
      if (not(tmpncomps.AT(m) == myncomps.AT(m))) {
        cout << "Recompose.cc" << endl
             << "superregss=" << superregss << endl
             << "regss=" << regss << endl
             << "nregs=" << nregs << endl
             << "newregs=" << newregs << endl
             << "newnregs=" << newnregs << endl
             << "m=" << m << endl
             << "tmpncomps.AT(m)=" << tmpncomps.AT(m) << endl
             << "myncomps.AT(m)=" << myncomps.AT(m) << endl;
        cout << "newregs:" << endl;
        for (int r = 0; r < newnregs; ++r) {
          int const m = newregs.AT(r).map - minmap;
          assert(m >= 0 and m < nmaps);
          cout << "r=" << r << endl
               << "newregs=" << newregs.at(r) << endl
               << "newregs.extent.size()=" << newregs.at(r).extent.size()
               << endl
               << "newregs.extent.empty()=" << newregs.at(r).extent.empty()
               << endl;
        }
        cout << "*regf:" << endl;
        for (int r = 0; r < nregs; ++r) {
          int const m = superregs.AT(r).map - minmap;
          ipfulltree *const regf = superregs.AT(r).processors;
          cout << "r=" << r << endl << "m=" << m << endl;
          cout << "*regf=" << *regf << endl;
        }
        cout.flush();
      }
      assert(tmpncomps.AT(m) == myncomps.AT(m));
    }
  }
  // Count componets per process
  {
    vector<int> ncompsperproc(real_nprocs, 0);
    for (int r = 0; r < newnregs; ++r) {
      int const p = newregs.AT(r).processor;
      ++ncompsperproc.at(p);
    }
    if (same_number_of_components_on_each_process) {
      for (int p = 0; p < real_nprocs; ++p) {
        assert(ncompsperproc.at(p) == ncomps);
      }
    }
  }

  // Distribute regions
  // Allocate regions
  assert((int)regss.size() == nmaps);
  for (int m = 0; m < nmaps; ++m) {
    assert(regss.AT(m).empty());
    regss.AT(m).reserve(myncomps.AT(m));
    superregss.AT(m).clear();
    superregss.AT(m).reserve(mynregs.AT(m));
  }
  // Assign regions
  for (int r = 0; r < newnregs; ++r) {
    int const m = newregs.AT(r).map - minmap;
    assert(m >= 0 and m < nmaps);
    regss.AT(m).push_back(newregs.AT(r));
  }
  for (int r = 0; r < nregs; ++r) {
    int const m = superregs.AT(r).map - minmap;
    assert(m >= 0 and m < nmaps);
    superregss.AT(m).push_back(superregs.AT(r));
  }
  // Output regions
  if (recompose_verbose) {
    cout << "SRMA superregss " << superregss << endl;
    cout << "SRMA regss " << regss << endl;
  }
  // Check sizes
  for (int m = 0; m < nmaps; ++m) {
#if 0
      assert (int(regss.AT(m).size()) == myncomps.AT(m) + empty_comps.AT(m));
#endif
    assert(int(regss.AT(m).size()) == myncomps.AT(m));
    assert(int(superregss.AT(m).size()) == mynregs.AT(m));
  }

  if (recompose_verbose)
    cout << "SRMA exit" << endl;
}

void SplitRegionsMaps_Balanced(cGH const *const cctkGH,
                               vector<vector<region_t> > &superregss,
                               vector<vector<region_t> > &regss) {
  DECLARE_CCTK_PARAMETERS;

  // Find map numbers
  int const nmaps = superregss.size();
  vector<int> map_index(maps, -1);
  for (int m = 0; m < nmaps; ++m) {
    assert(not superregss.AT(m).empty());
    if (not superregss.AT(m).empty()) {
      assert(map_index.AT(superregss.AT(m).AT(0).map) == -1);
      map_index.AT(superregss.AT(m).AT(0).map) = m;
    }
  }

  // Combine regions from all maps
  vector<region_t> regs;
  for (int m = 0; m < nmaps; ++m) {
    combine_regions(superregss.AT(m), regs);
  }

  // Create superregion tree structure
  {
    int const ncomps = regs.size();
    for (int c = 0; c < ncomps; ++c) {
      assert(regs.AT(c).processor == -1);
      assert(not regs.AT(c).processors);
      pseudoregion_t const preg(regs.AT(c).extent, -1);
      regs.AT(c).processors = new ipfulltree(preg);
    }
  }

  // Set up new superregions
  // Note: The tree structure is now reachable from both the
  // superregions and the regions
  for (int m = 0; m < nmaps; ++m) {
    superregss.AT(m).clear();
  }
  {
    int const ncomps = regs.size();
    for (int c = 0; c < ncomps; ++c) {
      int const m = map_index.AT(regs.AT(c).map);
      superregss.AT(m).push_back(regs.AT(c));
    }
  }

  // Split onto processes
  vector<vector<region_t> > split_regss;
  int const nprocs = CCTK_nProcs(cctkGH);
  int const nworkers = nprocs;
  balance(regs, split_regss, nworkers, maximum_imbalance,
          same_number_of_components_on_each_process);

  // Assign process numbers, and make the tree structure
  // inaccessible from the regions
  for (int p = 0; p < nprocs; ++p) {
    int const ncomps = split_regss.AT(p).size();
    for (int c = 0; c < ncomps; ++c) {
      // Set process number in tree structure
      assert(split_regss.AT(p).AT(c).processors);
      assert(split_regss.AT(p).AT(c).processors->is_leaf());
      assert(split_regss.AT(p).AT(c).processors->payload().component == -1);
      split_regss.AT(p).AT(c).processors->payload().component = p;
      // Set process number in region
      split_regss.AT(p).AT(c).processors = NULL;
      assert(split_regss.AT(p).AT(c).processor == -1);
      split_regss.AT(p).AT(c).processor = p;
    }
  }

  // Distribute by map, implicitly assigning component numbers
  for (int p = 0; p < nprocs; ++p) {
    int const ncomps = split_regss.AT(p).size();
    for (int c = 0; c < ncomps; ++c) {
      int const m = map_index.AT(split_regss.AT(p).AT(c).map);
      regss.AT(m).push_back(split_regss.AT(p).AT(c));
    }
  }
}

static void SplitRegionsMaps_AsSpecified(cGH const *cctkGH,
                                         vector<vector<region_t> > &superregss,
                                         vector<vector<region_t> > &regss) {
  DECLARE_CCTK_PARAMETERS;

  int const nmaps = superregss.size();
  assert(nmaps == 1); // for now
  assert((int)regss.size() == nmaps);

  SplitRegions_AsSpecified(cctkGH, superregss.AT(0), regss.AT(0));
}

//////////////////////////////////////////////////////////////////////////////

static void MakeMultigridBoxes(cGH const *const cctkGH, int const m,
                               ibbox const &base, region_t const &reg,
                               vector<region_t> &regs) {
  DECLARE_CCTK_PARAMETERS;

  regs.resize(mglevels, reg);
  if (mglevels > 1) {
    // boundary offsets
    jjvect nboundaryzones, is_internal, is_staggered, shiftout;
    if (domain_from_multipatch and
        CCTK_IsFunctionAliased("MultiPatch_GetBoundarySpecification")) {
      check(not MultiPatch_GetBoundarySpecification(
          m, 2 * dim, &nboundaryzones[0][0], &is_internal[0][0],
          &is_staggered[0][0], &shiftout[0][0]));
    } else {
      check(not GetBoundarySpecification(2 * dim, &nboundaryzones[0][0],
                                         &is_internal[0][0],
                                         &is_staggered[0][0], &shiftout[0][0]));
    }
    // (distance in grid points between the exterior and the physical boundary)
    iivect offset;
    for (int d = 0; d < dim; ++d) {
      for (int f = 0; f < 2; ++f) {
        assert(not is_staggered[d][f]);
        offset[d][f] = (+(is_internal[d][f] ? 0 : nboundaryzones[d][f] - 1) +
                        shiftout[d][f]);
      }
    }
    vector<ibbox> bases(mglevels);
    bases.AT(0) = base;
    for (int ml = 1; ml < mglevels; ++ml) {
      // next finer base
      ivect const fbaselo = bases.AT(ml - 1).lower();
      ivect const fbasehi = bases.AT(ml - 1).upper();
      ivect const fbasestr = bases.AT(ml - 1).stride();
      // this base
      ivect const basestr = fbasestr * mgfact;
      ivect const baselo =
          fbaselo +
          (xpose(offset)[0] - ivect(mgfact) * xpose(offset)[0]) * fbasestr;
      ivect const basehi =
          fbasehi +
          (xpose(offset)[1] - ivect(mgfact) * xpose(offset)[1]) * fbasestr;
      ivect const baselo1 = baselo;
      ivect const basehi1 = baselo1 + (basehi - baselo1) / basestr * basestr;
      bases.AT(ml) = ibbox(baselo1, basehi1, basestr);
      // next finer grid
      ivect const flo = regs.AT(ml - 1).extent.lower();
      ivect const fhi = regs.AT(ml - 1).extent.upper();
      ivect const fstr = regs.AT(ml - 1).extent.stride();
      // this grid
      ivect const str = fstr * mgfact;
      ivect const lo =
          flo +
          either(reg.outer_boundaries[0],
                 (xpose(offset)[0] - ivect(mgfact) * xpose(offset)[0]) * fstr,
                 ivect(0));
      ivect const hi =
          fhi +
          either(reg.outer_boundaries[1],
                 -(xpose(offset)[1] - ivect(mgfact) * xpose(offset)[1]) * fstr,
                 ivect(0));
      ivect const lo1 = baselo1 + (lo - baselo1 + str - 1) / str * str;
      ivect const hi1 = lo1 + (hi - lo1) / str * str;
      regs.AT(ml).extent = ibbox(lo1, hi1, str);
    }
  }
}

void MakeMultigridBoxes(cGH const *const cctkGH, int const m,
                        vector<vector<region_t> > const &regss,
                        vector<vector<vector<region_t> > > &regsss) {
  regsss.resize(mglevels);
  for (int ml = 0; ml < mglevels; ++ml) {
    regsss.AT(ml).resize(regss.size());
    for (int rl = 0; rl < (int)regss.size(); ++rl) {
      regsss.AT(ml).AT(rl).resize(regss.AT(rl).size());
    }
  }

  for (int rl = 0; rl < (int)regss.size(); ++rl) {

    ibbox base;
    for (int c = 0; c < (int)regss.AT(rl).size(); ++c) {
      base = base.expanded_containing(regss.AT(rl).AT(c).extent);
    }

    for (int c = 0; c < (int)regss.AT(rl).size(); ++c) {
      vector<region_t> mg_regs;
      MakeMultigridBoxes(cctkGH, m, base, regss.AT(rl).AT(c), mg_regs);
      assert((int)mg_regs.size() == mglevels);
      for (int ml = 0; ml < mglevels; ++ml) {
        regsss.AT(ml).AT(rl).AT(c) = mg_regs.AT(ml);
      }

    } // for c

  } // for rl
}

void MakeMultigridBoxesMaps(
    cGH const *const cctkGH, vector<vector<vector<region_t> > > const &regsss,
    vector<vector<vector<vector<region_t> > > > &regssss) {
  regssss.resize(maps);
  for (int m = 0; m < maps; ++m) {
    MakeMultigridBoxes(cctkGH, m, regsss.AT(m), regssss.AT(m));
  } // for m
}

static void ClassifyPoints(cGH const *const cctkGH, int const rl) {
  // negative: needs to be set explicitly (e.g. boundary)
  // zero:     unused (e.g. ghost)
  // positive: needs to be evolved
  // -1:      boundary point (needs to be set explicitly)
  //  0:      unused (e.g. ghost point, or restriction target)
  //  n=1..N: evolved, used for integrator substeps i<=n
  //          (i=N..1, counting backwards; see MoL documentation)
  //          i.e.: n=1: used between time steps (i.e., should be visualised)
  //                n>1: used only while time stepping (e.g. buffer zones)

  BEGIN_META_MODE(cctkGH) {
    BEGIN_MGLEVEL_LOOP(cctkGH) {
      ENTER_LEVEL_MODE(cctkGH, rl) {
        BEGIN_LOCAL_MAP_LOOP(cctkGH, CCTK_GF) {
          BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, CCTK_GF) {
            DECLARE_CCTK_ARGUMENTS;
            if (point_class) {
#pragma omp parallel
              CCTK_LOOP3_ALL(CarpetClassifyPoints, cctkGH, i, j, k) {
                int const ind = CCTK_GFINDEX3D(cctkGH, i, j, k);
                point_class[ind] = 1;
              }
              CCTK_ENDLOOP3_ALL(CarpetClassifyPoints);
            }
          }
          END_LOCAL_COMPONENT_LOOP;
        }
        END_LOCAL_MAP_LOOP;
      }
      LEAVE_LEVEL_MODE;
    }
    END_MGLEVEL_LOOP;
  }
  END_META_MODE;
}

} // namespace Carpet
