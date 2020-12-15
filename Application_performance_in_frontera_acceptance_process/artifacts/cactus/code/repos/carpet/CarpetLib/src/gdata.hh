#ifndef GDATA_HH
#define GDATA_HH

#include <cctk.h>

#include <cassert>
#include <cstdlib>
#include <queue>
#include <iostream>
#include <string>
#include <vector>

#include "bbox.hh"
#include "commstate.hh"
#include "defs.hh"
#include "dist.hh"
#include "operators.hh"
#include "timestat.hh"
#include "vect.hh"

using namespace std;

// Slabbing description
template <typename T, int D> struct slab {
  // For periodic boundary conditions:
  vect<T, D> offset; // dst[ipos] = src[ipos + offset * box.stride]
  // For refluxing:
  vect<T, D> is_centered;
  slab() : offset(0), is_centered(1) {}
};
typedef slab<int, dim> islab;

template <typename T, int D>
ostream &operator<<(ostream &os, slab<T, D> const &slabinfo);

template <typename T, int D>
MPI_Datatype mpi_datatype(slab<T, D> const &) CCTK_ATTRIBUTE_CONST;
namespace dist {
template <> inline MPI_Datatype mpi_datatype<islab>() CCTK_ATTRIBUTE_CONST;
template <> inline MPI_Datatype mpi_datatype<islab>() {
  islab dummy;
  return mpi_datatype(dummy);
}
}

// A generic data storage without type information
class gdata {

  static set<gdata *> allgdata;

protected: // should be readonly
  // Fields
  void *_storage; // A copy of the storage pointer

public:
  const int varindex; // Cactus variable index, or -1

protected:
  centering cent;
  operator_type transport_operator;

  bool _has_storage; // has storage associated (on some process)
  int _size;         // size (number of elements including padding)

  int _proc; // stored on process

  ivect _shape;                 // shape
  ivect _padded_shape, _stride; // allocated shape and index order

  ibbox _extent; // bbox for all data

  bool comm_active;    // a communication is going on
  MPI_Request request; // outstanding MPI request

private:
  // Forbid copying and passing by value
  gdata(gdata const &);
  gdata &operator=(gdata const &);

public:
  // Constructors
  gdata(const int varindex, const centering cent = error_centered,
        const operator_type transport_operator = op_error);

  // Destructors
  virtual ~gdata();

  // Pseudo constructors
  virtual gdata *
  make_typed(const int varindex, const centering cent = error_centered,
             const operator_type transport_operator = op_error) const = 0;

  // Storage management
  virtual void allocate(const ibbox &extent, const int proc,
                        void *const memptr = NULL,
                        size_t const memsize = 0) = 0;
  virtual void free() = 0;
  virtual size_t allocsize(const ibbox &extent, const int proc) const = 0;

  // true if fence is intact
  virtual bool check_fence(const int upperlower) const = 0;
  static bool fence_is_energized();

  // Accessors
  bool has_storage() const { return _has_storage; }

  void const *storage() const {
    assert(_has_storage);
    return _storage;
  }

  void *storage() {
    assert(_has_storage);
    return _storage;
  }

  int size() const {
    assert(_has_storage);
    return _size;
  }

  int proc() const {
    assert(_has_storage);
    return _proc;
  }

  const ivect &shape() const {
    assert(_has_storage);
    return _shape;
  }

  const ivect &padded_shape() const {
    assert(_has_storage);
    return _padded_shape;
  }

  const ivect &stride() const {
    assert(_has_storage);
    return _stride;
  }

  const ibbox &extent() const {
    assert(_has_storage);
    return _extent;
  }

  int elementsize() const { return c_datatype_size(); }

  // Data accessors
  int offset(const ivect &index) const {
    assert(_has_storage);
    assert(all((index - extent().lower()) % extent().stride() == 0));
    ivect const ind = (index - extent().lower()) / extent().stride();
    assert(all(ind >= 0 and ind < shape()));
    int const off = dot(ind, stride());
    assert(off >= 0 and off < size());
    return off;
  }

private:
  // Datatype accessors
  // maps the C datatype of a data class object to a 0-based index
  virtual unsigned int c_datatype() const CCTK_MEMBER_ATTRIBUTE_PURE = 0;
  virtual size_t c_datatype_size() const CCTK_MEMBER_ATTRIBUTE_PURE = 0;

  // Data manipulators

public:
  void copy_from(comm_state &state, gdata const *src, ibbox const &dstbox,
                 ibbox const &srcbox, islab const *restrict const slabinfo,
                 int dstproc, int srcproc);

  void transfer_from(comm_state &state, vector<gdata const *> const &srcs,
                     vector<CCTK_REAL> const &times, ibbox const &dstbox,
                     ibbox const &srcbox, islab const *restrict const slabinfo,
                     int dstproc, int srcproc, CCTK_REAL time, int order_space,
                     int order_time);

protected:
  void find_source_timelevel(vector<CCTK_REAL> const &times, CCTK_REAL time,
                             int order_time, operator_type transport_operator,
                             int &timelevel0, int &ntimelevels);

private:
  virtual void copy_from_innerloop(gdata const *gsrc, ibbox const &dstbox,
                                   ibbox const &srcbox,
                                   islab const *slabinfo) = 0;

  virtual void transfer_from_innerloop(vector<gdata const *> const &gsrcs,
                                       vector<CCTK_REAL> const &times,
                                       ibbox const &dstbox, ibbox const &srcbox,
                                       islab const *restrict const slabinfo,
                                       CCTK_REAL time, int order_space,
                                       int order_time) = 0;

public:
  virtual size_t memory() const CCTK_MEMBER_ATTRIBUTE_PURE = 0;
  static size_t allmemory() CCTK_MEMBER_ATTRIBUTE_PURE;
  virtual ostream &output(ostream &os) const = 0;
};

inline size_t memoryof(gdata const &d) { return d.memory(); }

inline ostream &operator<<(ostream &os, const gdata &d) { return d.output(os); }

#endif // GDATA_HH
