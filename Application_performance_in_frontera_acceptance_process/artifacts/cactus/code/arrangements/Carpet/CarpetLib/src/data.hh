#ifndef DATA_HH
#define DATA_HH

#include <cctk.h>

#include <cassert>
#include <iostream>
#include <string>
#include <vector>

#include "defs.hh"
#include "dist.hh"
#include "bbox.hh"
#include "gdata.hh"
#include "mem.hh"
#include "vect.hh"

using namespace std;

template <typename T> class data;

// A distributed multi-dimensional array
template <typename T> class data : public gdata {

  // Fields
  mem<T> *_memory; // the data (if located on this process)

  // For vector groups with contiguous storage
  int vectorlength;   // number of vector elements
  int vectorindex;    // index of this vector element
  data *vectorleader; // if index!=0: first vector element

private:
  // Forbid copying and passing by value
  data(data const &);
  data &operator=(data const &);

public:
  // Constructors
  data(const int varindex = -1, const centering cent = error_centered,
       const operator_type transport_operator = op_error,
       const int vectorlength = 1, const int vectorindex = 0,
       data *const vectorleader = NULL);
  data(const int varindex, const centering cent,
       const operator_type transport_operator, const int vectorlength,
       const int vectorindex, data *const vectorleader, const ibbox &extent,
       const int proc);

  // Destructors
  virtual ~data();

  // Pseudo constructors
  virtual data *make_typed(const int varindex, const centering cent,
                           const operator_type transport_operator) const;

  // Storage management
  virtual void allocate(const ibbox &extent, const int proc,
                        void *const memptr = NULL, size_t const memsize = 0);
  virtual void free();
  virtual size_t allocsize(const ibbox &extent, const int proc) const;
  virtual bool check_fence(const int upperlower) const;

public:
  // Data accessors
  const T &operator[](const ivect &index) const {
    assert(_has_storage);
    assert(_memory);
    return _memory->storage(vectorindex)[offset(index)];
  }

  T &operator[](const ivect &index) {
    assert(_has_storage);
    assert(_memory);
    return _memory->storage(vectorindex)[offset(index)];
  }

  // Datatype accessors
private:
  // maps the C datatype of a data class object to a 0-based index
  virtual unsigned int c_datatype() const {
    T dummy;
    return dist::c_datatype(dummy);
  }
  // size of the C datatype
  virtual size_t c_datatype_size() const { return sizeof(T); }

  // Data manipulators

private:
  void copy_from_innerloop(gdata const *gsrc, ibbox const &dstbox,
                           ibbox const &srcbox, islab const *slabinfo);

  void transfer_from_innerloop(vector<gdata const *> const &gsrcs,
                               vector<CCTK_REAL> const &times,
                               ibbox const &dstbox, ibbox const &srcbox,
                               islab const *restrict const slabinfo,
                               CCTK_REAL time, int order_space, int order_time);

  void transfer_time(vector<gdata const *> const &gsrcs,
                     vector<CCTK_REAL> const &times, ibbox const &dstbox,
                     ibbox const &srcbox, islab const *restrict const slabinfo,
                     CCTK_REAL time, int order_space, int order_time);

  void transfer_p_r(data const *const src, ibbox const &dstbox,
                    ibbox const &srcbox, islab const *restrict const slabinfo,
                    int order_space);

  void transfer_p_vc_cc(data const *const src, ibbox const &dstbox,
                        ibbox const &srcbox, int order_space);

  void transfer_prolongate(data const *const src, ibbox const &dstbox,
                           ibbox const &srcbox, int order_space);

  void transfer_restrict(data const *const src, ibbox const &dstbox,
                         ibbox const &srcbox,
                         islab const *restrict const slabinfo, int order_space);

  void time_interpolate(vector<data *> const &srcs, ibbox const &dstbox,
                        ibbox const &srcbox, vector<CCTK_REAL> const &times,
                        CCTK_REAL time, int order_time);

public:
  // Memory usage
  virtual size_t memory() const CCTK_MEMBER_ATTRIBUTE_PURE;

  // Output
  virtual ostream &output(ostream &os) const;
};

#endif // DATA_HH
