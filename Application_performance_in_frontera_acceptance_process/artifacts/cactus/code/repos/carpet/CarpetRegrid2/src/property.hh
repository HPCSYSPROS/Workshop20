#ifndef PROPERTY_HH
#define PROPERTY_HH

// Consistency properties for the grid structure

#include <vector>

#include <bboxset.hh>
#include <defs.hh>
#include <dh.hh>
#include <gh.hh>

namespace CarpetRegrid2 {

// Each property consists of a test, which returns true or false
// depending on whether the property is satisfied, and an action
// that enforces the property.
class property {
protected:
  virtual bool test_impl(gh const &hh, dh const &dd, level_boundary const &bnd,
                         vector<ibset> const &regions, int rl) = 0;
  virtual void enforce_impl(gh const &hh, dh const &dd,
                            level_boundary const &bnd, vector<ibset> &regions,
                            int rl) = 0;

public:
  virtual ~property() {}
  virtual const char *name() = 0;
  bool test(gh const &hh, dh const &dd, level_boundary const &bnd,
            vector<ibset> const &regions, int rl);
  void enforce(gh const &hh, dh const &dd, level_boundary const &bnd,
               vector<ibset> &regions, int rl);
};

// Ensure that this grid contains the next finer grid
class proper_nesting : public property {
  ibset enlarged_fine_grid(gh const &hh, dh const &dd,
                           level_boundary const &bnd,
                           vector<ibset> const &regions, int rl);
  bool test_impl(gh const &hh, dh const &dd, level_boundary const &bnd,
                 vector<ibset> const &regions, int rl);
  void enforce_impl(gh const &hh, dh const &dd, level_boundary const &bnd,
                    vector<ibset> &regions, int rl);
  const char* name() { return "proper_nesting"; };
};

// Add buffer zones (do this only once)
class add_buffers : public property {
  ibset buffered_regions(gh const &hh, dh const &dd, level_boundary const &bnd,
                         vector<ibset> const &regions, int rl);
  bool test_impl(gh const &hh, dh const &dd, level_boundary const &bnd,
                 vector<ibset> const &regions, int rl);
  void enforce_impl(gh const &hh, dh const &dd, level_boundary const &bnd,
                    vector<ibset> &regions, int rl);
  const char* name() { return "add_buffers"; };
};

// Combine all regions into a single region, if this is worthwhile
class combine_regions : public property {
  ibbox combined_regions(gh const &hh, dh const &dd, level_boundary const &bnd,
                         vector<ibset> const &regions, int rl);
  bool test_impl(gh const &hh, dh const &dd, level_boundary const &bnd,
                 vector<ibset> const &regions, int rl);
  void enforce_impl(gh const &hh, dh const &dd, level_boundary const &bnd,
                    vector<ibset> &regions, int rl);
  const char* name() { return "combine_regions"; };
};

// Align the boxes with the next coarser grid
class snap_coarse : public property {
  ibset snapped_regions(gh const &hh, dh const &dd, level_boundary const &bnd,
                        vector<ibset> const &regions, int rl);
  bool test_impl(gh const &hh, dh const &dd, level_boundary const &bnd,
                 vector<ibset> const &regions, int rl);
  void enforce_impl(gh const &hh, dh const &dd, level_boundary const &bnd,
                    vector<ibset> &regions, int rl);
  const char* name() { return "snap_coarse"; };
};

// Make the boxes rotating-90 symmetric
class rotsym90 : public property {
  ibset symmetrised_regions(gh const &hh, dh const &dd,
                            level_boundary const &bnd,
                            vector<ibset> const &regions, int rl);
  bool test_impl(gh const &hh, dh const &dd, level_boundary const &bnd,
                 vector<ibset> const &regions, int rl);
  void enforce_impl(gh const &hh, dh const &dd, level_boundary const &bnd,
                    vector<ibset> &regions, int rl);
  const char* name() { return "rotsym90"; };
};

// Make the boxes parity symmetric
class parsym : public property {
  ibset symmetrised_regions(gh const &hh, dh const &dd,
                            level_boundary const &bnd,
                            vector<ibset> const &regions, int rl);
  bool test_impl(gh const &hh, dh const &dd, level_boundary const &bnd,
                 vector<ibset> const &regions, int rl);
  void enforce_impl(gh const &hh, dh const &dd, level_boundary const &bnd,
                    vector<ibset> &regions, int rl);
  const char* name() { return "parsym"; };
};

// Make the boxes rotating-180 symmetric
class rotsym180 : public property {
  ibset symmetrised_regions(gh const &hh, dh const &dd,
                            level_boundary const &bnd,
                            vector<ibset> const &regions, int rl);
  bool test_impl(gh const &hh, dh const &dd, level_boundary const &bnd,
                 vector<ibset> const &regions, int rl);
  void enforce_impl(gh const &hh, dh const &dd, level_boundary const &bnd,
                    vector<ibset> &regions, int rl);
  const char* name() { return "rotsym180"; };
};

// Make the boxes periodic in one direction
template <int dir> class periodic : public property {
  ibset symmetrised_regions(gh const &hh, dh const &dd,
                            level_boundary const &bnd,
                            vector<ibset> const &regions, int rl);
  bool test_impl(gh const &hh, dh const &dd, level_boundary const &bnd,
                 vector<ibset> const &regions, int rl);
  void enforce_impl(gh const &hh, dh const &dd, level_boundary const &bnd,
                    vector<ibset> &regions, int rl);
  const char* name() { return "periodic"; };
};

// Clip at the outer boundary
class boundary_clip : public property {
  ibset clipped_regions(gh const &hh, dh const &dd, level_boundary const &bnd,
                        vector<ibset> const &regions, int rl);
  bool test_impl(gh const &hh, dh const &dd, level_boundary const &bnd,
                 vector<ibset> const &regions, int rl);
  void enforce_impl(gh const &hh, dh const &dd, level_boundary const &bnd,
                    vector<ibset> &regions, int rl);
  const char* name() { return "boundary_clip"; };
};

// Ensure that this grid is contained in the domain
class in_domain : public property {
  bool test_impl(gh const &hh, dh const &dd, level_boundary const &bnd,
                 vector<ibset> const &regions, int rl);
  void enforce_impl(gh const &hh, dh const &dd, level_boundary const &bnd,
                    vector<ibset> &regions, int rl);
  const char* name() { return "in_domain"; };
};

// Ensure that this grid is symmetric, if desired
class is_symmetric : public property {
  ibset symmetrised_regions(gh const &hh, dh const &dd,
                            level_boundary const &bnd,
                            vector<ibset> const &regions, int rl);
  bool test_impl(gh const &hh, dh const &dd, level_boundary const &bnd,
                 vector<ibset> const &regions, int rl);
  void enforce_impl(gh const &hh, dh const &dd, level_boundary const &bnd,
                    vector<ibset> &regions, int rl);
  const char* name() { return "is_symmetric"; };
};

} // namespace CarpetRegrid2

#endif // #ifndef PROPERTY_HH
