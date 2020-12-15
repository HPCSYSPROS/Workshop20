#ifndef GGF_HH
#define GGF_HH

#include <cctk.h>

#include <cassert>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include "defs.hh"
#include "dh.hh"
#include "gdata.hh"
#include "gh.hh"
#include "th.hh"

using namespace std;

// Forward declaration
class ggf;

// Output
ostream &operator<<(ostream &os, const ggf &f);

// A generic grid function without type information
class ggf {

  static set<ggf *> allggf;

  // Types
  typedef vector<pseudoregion_t> pvect;
  typedef vector<sendrecv_pseudoregion_t> srpvect;

  typedef gdata *tdata;        // data ...
  typedef vector<tdata> fdata; // ... for each time level
  typedef vector<fdata> cdata; // ... for each local component
  typedef vector<cdata> rdata; // ... for each refinement level
  typedef vector<rdata> mdata; // ... for each multigrid level

public: // should be readonly
  // Fields
  const int varindex; // Cactus variable index
  const operator_type transport_operator;

  const th &t;                       // time hierarchy
  const int prolongation_order_time; // order of temporal prolongation operator

  const gh &h; // grid hierarchy
  dh &d;       // data hierarchy

protected:
  vector<vector<int> > timelevels_; // time levels [ml][rl]

  mdata storage; // storage

public:
  const int vectorlength;  // vector length
  const int vectorindex;   // index of *this
  ggf *const vectorleader; // first vector element

private:
  mdata oldstorage; // temporary storage

public:
  // Constructors
  ggf(const int varindex, const operator_type transport_operator, th &t, dh &d,
      const int prolongation_order_time, const int vectorlength,
      const int vectorindex, ggf *const vectorleader);

  // Destructors
  virtual ~ggf();

  // Comparison
  bool operator==(const ggf &f) const CCTK_MEMBER_ATTRIBUTE_PURE;

  // Querying
  int timelevels(int const ml, int const rl) const {
    return timelevels_.AT(ml).AT(rl);
  }

  // Modifiers
  void set_timelevels(int ml, int rl, int new_timelevels);

  void recompose_crop();
  void recompose_allocate(int rl);
  void recompose_fill(comm_state &state, int rl, bool do_prolongate);
  void recompose_free_old(int rl);
  void recompose_free(int rl);

  // Cycle the time levels by rotating the data sets
  void cycle_all(int rl, int ml);

  // Uncycle the time levels by rotating the data sets
  void uncycle_all(int rl, int ml);

  // Flip the time levels by exchanging the data sets
  void flip_all(int rl, int ml);

  // Fill all time levels from the current time level
  void fill_all(int rl, int ml);

  // The grid boundaries have to be updated after calling mg_restrict,
  // mg_prolongate, ref_restrict, or ref_prolongate.

  // "Updating" means here that the boundaries have to be
  // synchronised. They don't need to be prolongated.

  // Synchronise the boundaries of a component
  void sync_all(comm_state &state, int tl, int rl, int ml);

  // Prolongate the boundaries of a component
  void ref_bnd_prolongate_all(comm_state &state, int tl, int rl, int ml,
                              CCTK_REAL time);

  // Restrict a multigrid level
  void mg_restrict_all(comm_state &state, int tl, int rl, int ml,
                       CCTK_REAL time);

  // Prolongate a multigrid level
  void mg_prolongate_all(comm_state &state, int tl, int rl, int ml,
                         CCTK_REAL time);

  // Restrict a refinement level
  void ref_restrict_all(comm_state &state, int tl, int rl, int ml);

  // Prolongate a refinement level
  void ref_prolongate_all(comm_state &state, int tl, int rl, int ml,
                          CCTK_REAL time);

  // Reflux a refinement level
  void ref_reflux_all(comm_state &state, int tl, int rl, int ml, int dir,
                      int face);

  // Reflux a refinement level
  void ref_reflux_prolongate_all(comm_state &state, int tl, int rl, int ml,
                                 int dir, int face);

  // Helpers

  virtual gdata *typed_data(int tl, int rl, int lc, int ml) const = 0;
  virtual gdata *new_typed_data() const = 0;

protected:
  // Transfer regions
  void transfer_from_all(comm_state &state, int tl1, int rl1, int ml1,
                         srpvect const dh::fast_dboxes::*sendrecvs,
                         vector<int> const &tl2s, int rl2, int ml2,
                         CCTK_REAL const &time, bool use_old_storage = false,
                         bool flip_send_recv = false,
                         islab const *restrict slabinfo = NULL);

  void transfer_from_all(comm_state &state, int tl1, int rl1, int ml1,
                         srpvect const dh::fast_dboxes::*sendrecvs, int tl2,
                         int rl2, int ml2, bool use_old_storage = false,
                         bool flip_send_recv = false,
                         islab const *restrict slabinfo = NULL) {
    vector<int> tl2s(1);
    tl2s.AT(0) = tl2;
    CCTK_REAL const time = t.get_time(ml2, rl2, tl2);
    transfer_from_all(state, tl1, rl1, ml1, sendrecvs, tl2s, rl2, ml2, time,
                      use_old_storage, flip_send_recv, slabinfo);
  }

public:
  // Access to the data
  const gdata *data_pointer(int const tl, int const rl, int const lc,
                            int const ml) const {
    assert(rl >= 0 and rl < h.reflevels());
    assert(lc >= 0 and lc < h.local_components(rl));
    assert(ml >= 0 and ml < h.mglevels());
    assert(tl >= 0 and tl < timelevels(ml, rl));
    return storage.AT(ml).AT(rl).AT(lc).AT(tl);
  }
  gdata *data_pointer(int const tl, int const rl, int const lc, int const ml) {
    assert(rl >= 0 and rl < h.reflevels());
    assert(lc >= 0 and lc < h.local_components(rl));
    assert(ml >= 0 and ml < h.mglevels());
    assert(tl >= 0 and tl < timelevels(ml, rl));
    return storage.AT(ml).AT(rl).AT(lc).AT(tl);
  }

  // Output
  virtual size_t memory() const CCTK_MEMBER_ATTRIBUTE_PURE = 0;
  static size_t allmemory() CCTK_MEMBER_ATTRIBUTE_PURE;
  virtual ostream &output(ostream &os) const = 0;

private:
  ggf();                       // canonical default construtor
  ggf(const ggf &);            // canonical copy construtor
  ggf &operator=(const ggf &); // canonical copy
};

inline size_t memoryof(ggf const &f) { return f.memory(); }

inline ostream &operator<<(ostream &os, const ggf &f) { return f.output(os); }

#endif // GGF_HH
