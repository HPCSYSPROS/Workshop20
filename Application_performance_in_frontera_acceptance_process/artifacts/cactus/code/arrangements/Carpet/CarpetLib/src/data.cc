#include <cctk.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <cassert>
#include <climits>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "bbox.hh"
#include "bboxset.hh"
#include "cacheinfo.hh"
#include "defs.hh"
#include "dist.hh"
#include "timestat.hh"
#include "vect.hh"

#include "data.hh"
#include "operator_prototypes_3d.hh"
#include "operator_prototypes_4d.hh"

using namespace std;
using namespace CarpetLib;

template <typename T>
static void
call_operator(void (*the_operator)(
                  T const *restrict const src, ivect3 const &restrict srcpadext,
                  ivect3 const &restrict srcext, T *restrict const dst,
                  ivect3 const &restrict dstpadext,
                  ivect3 const &restrict dstext, ibbox3 const &restrict srcbbox,
                  ibbox3 const &restrict dstbbox,
                  ibbox3 const &restrict srcregbbox,
                  ibbox3 const &restrict dstregbbox, void *const extraargs),
              T const *restrict const src, ivect3 const &restrict srcpadext,
              ivect3 const &restrict srcext, T *restrict const dst,
              ivect3 const &restrict dstpadext, ivect3 const &restrict dstext,
              ibbox3 const &restrict srcbbox, ibbox3 const &restrict dstbbox,
              ibbox3 const &restrict srcregbbox,
              ibbox3 const &restrict dstregbbox, void *const extraargs) {
  DECLARE_CCTK_PARAMETERS;

  if (use_loopcontrol_in_operators) {
    (*the_operator)(src, srcpadext, srcext, dst, dstpadext, dstext, srcbbox,
                    dstbbox, srcregbbox, dstregbbox, extraargs);
  } else if (use_openmp) {
#ifdef CARPET_DEBUG
    ibset alldstregbboxes;
#endif
#pragma omp parallel
    {
      int const num_threads = omp_get_num_threads();
      int const thread_num = omp_get_thread_num();
      // Parallelise in z direction
      // int const dir = 2;
      // Parallelise along longest extent
      int const dir = maxloc(dstregbbox.shape());
      int const stride = dstregbbox.stride()[dir];
      int const first_point = dstregbbox.lower()[dir];
      int const last_point = dstregbbox.upper()[dir] + stride;
      int const num_points = last_point - first_point;
      assert(num_points >= 0);
      assert(num_points % stride == 0);
      int const my_num_points =
          (num_points / stride + num_threads - 1) / num_threads * stride;
      int const my_first_point =
          min(last_point, first_point + thread_num * my_num_points);
      int const my_last_point = min(last_point, my_first_point + my_num_points);
      assert(my_last_point >= my_first_point);
      ibbox3 const mydstregbbox(
          dstregbbox.lower().replace(dir, my_first_point),
          dstregbbox.upper().replace(dir, my_last_point - stride),
          dstregbbox.stride());
      if (not mydstregbbox.empty()) {
        (*the_operator)(src, srcpadext, srcext, dst, dstpadext, dstext, srcbbox,
                        dstbbox, srcregbbox, mydstregbbox, extraargs);
#ifdef CARPET_DEBUG
#pragma omp critical
        alldstregbboxes += mydstregbbox;
#endif
      }
    }
#ifdef CARPET_DEBUG
    if (not(alldstregbboxes == ibset(dstregbbox))) {
      cout << "alldstregbboxes=" << alldstregbboxes << endl
           << "dstregbbox=" << dstregbbox << endl;
    }
    assert(alldstregbboxes == ibset(dstregbbox));
#endif
  }
}

template <typename T>
static void
call_operator(void (*the_operator)(
                  T const *restrict const src, ivect4 const &restrict srcpadext,
                  ivect4 const &restrict srcext, T *restrict const dst,
                  ivect4 const &restrict dstpadext,
                  ivect4 const &restrict dstext, ibbox4 const &restrict srcbbox,
                  ibbox4 const &restrict dstbbox,
                  ibbox4 const &restrict srcregbbox,
                  ibbox4 const &restrict dstregbbox, void *const extraargs),
              T const *restrict const src, ivect4 const &restrict srcpadext,
              ivect4 const &restrict srcext, T *restrict const dst,
              ivect4 const &restrict dstpadext, ivect4 const &restrict dstext,
              ibbox4 const &restrict srcbbox, ibbox4 const &restrict dstbbox,
              ibbox4 const &restrict srcregbbox,
              ibbox4 const &restrict dstregbbox, void *const extraargs) {
  DECLARE_CCTK_PARAMETERS;

  if (use_loopcontrol_in_operators) {
    (*the_operator)(src, srcpadext, srcext, dst, dstpadext, dstext, srcbbox,
                    dstbbox, srcregbbox, dstregbbox, extraargs);
  } else if (use_openmp) {
#ifdef CARPET_DEBUG
    ibset4 alldstregbboxes;
#endif
#pragma omp parallel
    {
      int const num_threads = omp_get_num_threads();
      int const thread_num = omp_get_thread_num();
      // Parallelise in z direction
      // int const dir = 2;
      // Parallelise along longest extent
      int const dir = maxloc(dstregbbox.shape());
      int const stride = dstregbbox.stride()[dir];
      int const first_point = dstregbbox.lower()[dir];
      int const last_point = dstregbbox.upper()[dir] + stride;
      int const num_points = last_point - first_point;
      assert(num_points >= 0);
      assert(num_points % stride == 0);
      int const my_num_points =
          (num_points / stride + num_threads - 1) / num_threads * stride;
      int const my_first_point =
          min(last_point, first_point + thread_num * my_num_points);
      int const my_last_point = min(last_point, my_first_point + my_num_points);
      assert(my_last_point >= my_first_point);
      ibbox4 const mydstregbbox(
          dstregbbox.lower().replace(dir, my_first_point),
          dstregbbox.upper().replace(dir, my_last_point - stride),
          dstregbbox.stride());
      if (not mydstregbbox.empty()) {
        (*the_operator)(src, srcpadext, srcext, dst, dstpadext, dstext, srcbbox,
                        dstbbox, srcregbbox, mydstregbbox, extraargs);
#ifdef CARPET_DEBUG
#pragma omp critical
        alldstregbboxes += mydstregbbox;
#endif
      }
    }
#ifdef CARPET_DEBUG
    if (not(alldstregbboxes == ibset4(dstregbbox))) {
      cout << "alldstregbboxes=" << alldstregbboxes << endl
           << "dstregbbox=" << dstregbbox << endl;
    }
    assert(alldstregbboxes == ibset4(dstregbbox));
#endif
  } else { // serial
    (*the_operator)(src, srcpadext, srcext, dst, dstpadext, dstext, srcbbox,
                    dstbbox, srcregbbox, dstregbbox, extraargs);
  }
}

// Fortran wrappers

template <typename T>
void prolongate_3d_eno(T const *restrict const /*src*/,
                       ivect3 const & /*srcpadext*/, ivect3 const & /*srcext*/,
                       T *restrict const /*dst*/, ivect3 const & /*dstpadext*/,
                       ivect3 const & /*dstext*/, ibbox3 const & /*srcbbox*/,
                       ibbox3 const & /*dstbbox*/,
                       ibbox3 const & /*srcregbbox*/,
                       ibbox3 const & /*dstregbbox*/, void * /*extraargs*/) {
  CCTK_ERROR("Data type not supported");
}

extern "C" void CCTK_FCALL CCTK_FNAME(prolongate_3d_real8_eno)(
    const CCTK_REAL8 *src, const int &srcipadext, const int &srcjpadext,
    const int &srckpadext, const int &srciext, const int &srcjext,
    const int &srckext, CCTK_REAL8 *dst, const int &dstipadext,
    const int &dstjpadext, const int &dstkpadext, const int &dstiext,
    const int &dstjext, const int &dstkext, const int srcbbox[3][3],
    const int dstbbox[3][3], const int dstregbbox[3][3]);

template <>
void prolongate_3d_eno(CCTK_REAL8 const *restrict const src,
                       ivect3 const &srcpadext, ivect3 const &srcext,
                       CCTK_REAL8 *restrict const dst, ivect3 const &dstpadext,
                       ivect3 const &dstext, ibbox3 const &srcbbox,
                       ibbox3 const &dstbbox, ibbox3 const &srcregbbox,
                       ibbox3 const &dstregbbox, void *extraargs) {
  assert(not extraargs);
  CCTK_FNAME(prolongate_3d_real8_eno)
  (src, srcpadext[0], srcpadext[1], srcpadext[2], srcext[0], srcext[1],
   srcext[2], dst, dstpadext[0], dstpadext[1], dstpadext[2], dstext[0],
   dstext[1], dstext[2], reinterpret_cast<int const(*)[3]>(&srcbbox),
   reinterpret_cast<int const(*)[3]>(&dstbbox),
   reinterpret_cast<int const(*)[3]>(&dstregbbox));
}

template <typename T>
void prolongate_3d_weno(T const *restrict const /*src*/,
                        ivect3 const & /*srcpadext*/, ivect3 const & /*srcext*/,
                        T *restrict const /*dst*/, ivect3 const & /*dstpadext*/,
                        ivect3 const & /*dstext*/, ibbox3 const & /*srcbbox*/,
                        ibbox3 const & /*dstbbox*/,
                        ibbox3 const & /*srcregbbox*/,
                        ibbox3 const & /*dstregbbox*/, void *extraargs) {
  CCTK_ERROR("Data type not supported");
}

extern "C" void CCTK_FCALL CCTK_FNAME(prolongate_3d_real8_weno)(
    const CCTK_REAL8 *src, const int &srcipadext, const int &srcjpadext,
    const int &srckpadext, const int &srciext, const int &srcjext,
    const int &srckext, CCTK_REAL8 *dst, const int &dstipadext,
    const int &dstjpadext, const int &dstkpadext, const int &dstiext,
    const int &dstjext, const int &dstkext, const int srcbbox[3][3],
    const int dstbbox[3][3], const int regbbox[3][3]);

template <>
void prolongate_3d_weno(CCTK_REAL8 const *restrict const src,
                        ivect3 const &srcpadext, ivect3 const &srcext,
                        CCTK_REAL8 *restrict const dst, ivect3 const &dstpadext,
                        ivect3 const &dstext, ibbox3 const &srcbbox,
                        ibbox3 const &dstbbox, ibbox3 const &srcregbbox,
                        ibbox3 const &dstregbbox, void *extraargs) {
  CCTK_FNAME(prolongate_3d_real8_weno)
  (src, srcpadext[0], srcpadext[1], srcpadext[2], srcext[0], srcext[1],
   srcext[2], dst, dstpadext[0], dstpadext[1], dstpadext[2], dstext[0],
   dstext[1], dstext[2], reinterpret_cast<int const(*)[3]>(&srcbbox),
   reinterpret_cast<int const(*)[3]>(&dstbbox),
   reinterpret_cast<int const(*)[3]>(&dstregbbox));
}

template <typename T>
void prolongate_3d_tvd(T const *restrict const /*src*/,
                       ivect3 const & /*srcpadext*/, ivect3 const & /*srcext*/,
                       T *restrict const /*dst*/, ivect3 const & /*dstpadext*/,
                       ivect3 const & /*dstext*/, ibbox3 const & /*srcbbox*/,
                       ibbox3 const & /*dstbbox*/,
                       ibbox3 const & /*srcregbbox*/,
                       ibbox3 const & /*dstregbbox*/, void * /*extraargs*/) {
  CCTK_ERROR("Data type not supported");
}

#ifndef OMIT_F90
extern "C" void CCTK_FCALL CCTK_FNAME(prolongate_3d_real8_tvd)(
    const CCTK_REAL8 *src, const int &srcipadext, const int &srcjpadext,
    const int &srckpadext, const int &srciext, const int &srcjext,
    const int &srckext, CCTK_REAL8 *dst, const int &dstipadext,
    const int &dstjpadext, const int &dstkpadext, const int &dstiext,
    const int &dstjext, const int &dstkext, const int srcbbox[3][3],
    const int dstbbox[3][3], const int regbbox[3][3]);

template <>
void prolongate_3d_tvd(CCTK_REAL8 const *restrict const src,
                       ivect3 const &srcpadext, ivect3 const &srcext,
                       CCTK_REAL8 *restrict const dst, ivect3 const &dstpadext,
                       ivect3 const &dstext, ibbox3 const &srcbbox,
                       ibbox3 const &dstbbox, ibbox3 const &srcregbbox,
                       ibbox3 const &dstregbbox, void *const extraargs) {
  assert(not extraargs);
  CCTK_FNAME(prolongate_3d_real8_tvd)
  (src, srcpadext[0], srcpadext[1], srcpadext[2], srcext[0], srcext[1],
   srcext[2], dst, dstpadext[0], dstpadext[1], dstpadext[2], dstext[0],
   dstext[1], dstext[2], reinterpret_cast<int const(*)[3]>(&srcbbox),
   reinterpret_cast<int const(*)[3]>(&dstbbox),
   reinterpret_cast<int const(*)[3]>(&dstregbbox));
}
#endif

template <typename T>
void prolongate_3d_cc_tvd(T const *restrict const /*src*/,
                          ivect3 const & /*srcpadext*/,
                          ivect3 const & /*srcext*/, T *restrict const /*dst*/,
                          ivect3 const & /*dstpadext*/,
                          ivect3 const & /*dstext*/, ibbox3 const & /*srcbbox*/,
                          ibbox3 const & /*dstbbox*/,
                          ibbox3 const & /*srcregbbox*/,
                          ibbox3 const & /*dstregbbox*/, void * /*extraargs*/) {
  CCTK_ERROR("Data type not supported");
}

#ifndef OMIT_F90
extern "C" void CCTK_FCALL CCTK_FNAME(prolongate_3d_cc_real8_tvd)(
    const CCTK_REAL8 *src, const int &srcipadext, const int &srcjpadext,
    const int &srckpadext, const int &srciext, const int &srcjext,
    const int &srckext, CCTK_REAL8 *dst, const int &dstipadext,
    const int &dstjpadext, const int &dstkpadext, const int &dstiext,
    const int &dstjext, const int &dstkext, const int srcbbox[3][3],
    const int dstbbox[3][3], const int regbbox[3][3]);

template <>
void prolongate_3d_cc_tvd(CCTK_REAL8 const *restrict const src,
                          ivect3 const &srcpadext, ivect3 const &srcext,
                          CCTK_REAL8 *restrict const dst,
                          ivect3 const &dstpadext, ivect3 const &dstext,
                          ibbox3 const &srcbbox, ibbox3 const &dstbbox,
                          ibbox3 const &srcregbbox, ibbox3 const &dstregbbox,
                          void *const extraargs) {
  assert(not extraargs);
  CCTK_FNAME(prolongate_3d_cc_real8_tvd)
  (src, srcpadext[0], srcpadext[1], srcpadext[2], srcext[0], srcext[1],
   srcext[2], dst, dstpadext[0], dstpadext[1], dstpadext[2], dstext[0],
   dstext[1], dstext[2], reinterpret_cast<int const(*)[3]>(&srcbbox),
   reinterpret_cast<int const(*)[3]>(&dstbbox),
   reinterpret_cast<int const(*)[3]>(&dstregbbox));
}
#endif

// Constructors
template <typename T>
data<T>::data(const int varindex_, const centering cent_,
              const operator_type transport_operator_, const int vectorlength_,
              const int vectorindex_, data *const vectorleader_)
    : gdata(varindex_, cent_, transport_operator_), _memory(NULL),
      vectorlength(vectorlength_), vectorindex(vectorindex_),
      vectorleader(vectorleader_) {
  assert(vectorlength >= 1);
  assert(vectorindex >= 0 and vectorindex < vectorlength);
  assert((vectorindex == 0 and not vectorleader) or
         (vectorindex != 0 and vectorleader));
}

template <typename T>
data<T>::data(const int varindex_, const centering cent_,
              const operator_type transport_operator_, const int vectorlength_,
              const int vectorindex_, data *const vectorleader_,
              const ibbox &extent_, const int proc_)
    : gdata(varindex_, cent_, transport_operator_), _memory(NULL),
      vectorlength(vectorlength_), vectorindex(vectorindex_),
      vectorleader(vectorleader_) {
  assert(vectorlength >= 1);
  assert(vectorindex >= 0 and vectorindex < vectorlength);
  assert((vectorindex == 0 and not vectorleader) or
         (vectorindex != 0 and vectorleader));
  allocate(extent_, proc_);
}

// Destructors
template <typename T> data<T>::~data() {
  if (_memory)
    free();
}

// Pseudo constructors
template <typename T>
data<T> *data<T>::make_typed(const int varindex_, const centering cent_,
                             const operator_type transport_operator_) const {
  return new data(varindex_, cent_, transport_operator_, 1, 0, NULL);
}

// Storage management
template <typename T>
void data<T>::allocate(const ibbox &extent_, const int proc_,
                       void *const memptr, size_t const memsize) {
  assert(not _has_storage);
  _has_storage = true;
  // prevent accidental wrap-around
  assert(all(extent_.lower() < numeric_limits<int>::max() / 2));
  assert(all(extent_.lower() > numeric_limits<int>::min() / 2));
  assert(all(extent_.upper() < numeric_limits<int>::max() / 2));
  assert(all(extent_.upper() > numeric_limits<int>::min() / 2));

  // data
  _extent = extent_;
  _shape = max(ivect(0), _extent.shape() / _extent.stride());
  _padded_shape = pad_shape(_shape);
  _stride[0] = 1;
  for (int d = 1; d < dim; ++d) {
    _stride[d] = _stride[d - 1] * _padded_shape[d - 1];
  }
  _size = _stride[dim - 1] * _padded_shape[dim - 1];

  _proc = proc_;
  if (dist::rank() == _proc) {
    if (vectorindex == 0) {
      assert(not vectorleader);
      _memory = new mem<T>(vectorlength, _size, (T *)memptr, memsize);
    } else {
      assert(vectorleader);
      _memory = vectorleader->_memory;
      assert(_memory);
    }
    _memory->register_client(vectorindex);
    _storage = _memory->storage(vectorindex);
  } else {
    assert(not memptr);
  }
}

template <typename T> void data<T>::free() {
  assert(_has_storage);
  assert(_memory);
  _memory->unregister_client(vectorindex);
  if (not _memory->has_clients())
    delete _memory;
  _memory = NULL;
  _has_storage = false;
  _storage = NULL;
}

template <typename T>
size_t data<T>::allocsize(const ibbox &extent_, const int proc_) const {
  if (dist::rank() != proc_)
    return 0;
  if (vectorindex != 0)
    return 0;
  assert(not vectorleader);
  return vectorlength * prod(pad_shape(extent_)) * sizeof(T);
}

template <typename T> bool data<T>::check_fence(const int upperlower) const {
  assert((vectorleader != NULL) xor (vectorindex == 0));
  bool retval = true;
  // since vectors share _memory we only check the first/last vector element
  if ((vectorindex == 0 and upperlower == 0) or
      (vectorindex == vectorlength - 1 and upperlower == 1)) {
    retval = _memory->is_fence_intact(upperlower);
  }
  return retval;
}

// Data manipulators

template <typename T>
void data<T>::copy_from_innerloop(gdata const *const gsrc,
                                  ibbox const &dstregbox,
                                  ibbox const &srcregbox,
                                  islab const *restrict const slabinfo) {
  // data const * const src = dynamic_cast <data const *> (gsrc);
  data const *const src = (data const *)gsrc;
  assert(has_storage() and src->has_storage());

  assert(proc() == src->proc());
  assert(dist::rank() == proc());

  ibbox const &srcbox = src->extent();
  ibbox const &dstbox = this->extent();

#if CARPET_DIM == 3
  // Don't use call_operator, because we parallelise ourselves
  copy_3d(static_cast<T const *>(src->storage()), src->padded_shape(),
          src->shape(), static_cast<T *>(this->storage()), this->padded_shape(),
          this->shape(), srcbox, dstbox, srcregbox, dstregbox,
          (void *)slabinfo);
#elif CARPET_DIM == 4
  // Don't use call_operator, because we parallelise ourselves
  copy_4d(static_cast<T const *>(src->storage()), src->padded_shape(),
          src->shape(), static_cast<T *>(this->storage()), this->padded_shape(),
          this->shape(), srcbox, dstbox, srcregbox, dstregbox,
          (void *)slabinfo);
#else
#error "Value for CARPET_DIM not supported"
#endif
}

template <typename T>
void data<T>::transfer_from_innerloop(vector<gdata const *> const &gsrcs,
                                      vector<CCTK_REAL> const &times,
                                      ibbox const &dstbox, ibbox const &srcbox,
                                      islab const *restrict const slabinfo,
                                      CCTK_REAL const time,
                                      int const order_space,
                                      int const order_time) {
  assert(has_storage());
  assert(dist::rank() == proc());
  for (size_t tl = 0; tl < gsrcs.size(); ++tl) {
    if (gsrcs.AT(tl)) {
      assert(gsrcs.AT(tl)->has_storage());
      assert(gsrcs.AT(tl)->proc() == proc());
    }
  }

  transfer_time(gsrcs, times, dstbox, srcbox, slabinfo, time, order_space,
                order_time);
}

template <typename T>
void data<T>::transfer_time(vector<gdata const *> const &gsrcs,
                            vector<CCTK_REAL> const &times, ibbox const &dstbox,
                            ibbox const &srcbox,
                            islab const *restrict const slabinfo,
                            CCTK_REAL const time, int const order_space,
                            int const order_time) {
  // Use this timelevel, or interpolate in time if set to -1
  int timelevel0, ntimelevels;
  find_source_timelevel(times, time, order_time, transport_operator, timelevel0,
                        ntimelevels);

  if (ntimelevels > 1) {
    // Time interpolation is necessary
    assert(timelevel0 == 0);

    assert((int)gsrcs.size() >= ntimelevels);
    assert((int)times.size() >= ntimelevels);

    vector<data *> tmps(timelevel0 + ntimelevels, (data *)0);

    for (int tl = timelevel0; tl < timelevel0 + ntimelevels; ++tl) {
      tmps.AT(tl) =
          new data(this->varindex, this->cent, this->transport_operator);
      tmps.AT(tl)->allocate(dstbox, this->proc());

      assert(gsrcs.AT(tl));
      // data const * const src = dynamic_cast <data const *> (gsrcs.AT(tl));
      data const *const src = (data const *)gsrcs.AT(tl);
      tmps.AT(tl)->transfer_p_r(src, dstbox, srcbox, slabinfo, order_space);
    }

    time_interpolate(tmps, dstbox, dstbox, times, time, order_time);

    for (int tl = timelevel0; tl < timelevel0 + ntimelevels; ++tl) {
      delete tmps.AT(tl);
    }

  } else {
    // No time interpolation

    assert((int)gsrcs.size() > timelevel0);
    assert((int)times.size() > timelevel0);

    // data const * const src = dynamic_cast <data const *>
    // (gsrcs.AT(timelevel0));
    data const *const src = (data const *)gsrcs.AT(timelevel0);

    transfer_p_r(src, dstbox, srcbox, slabinfo, order_space);
  } // if
}

template <typename T>
void data<T>::transfer_p_r(data const *const src, ibbox const &dstbox,
                           ibbox const &srcbox,
                           islab const *restrict const slabinfo,
                           int const order_space) {
  if (all(src->extent().stride() == this->extent().stride())) {
    // Copy
    copy_from_innerloop(src, dstbox, srcbox, slabinfo);
  } else if (all(src->extent().stride() > this->extent().stride())) {
    // Prolongate
    assert(transport_operator != op_sync and transport_operator != op_restrict);
    assert(not slabinfo);
    transfer_p_vc_cc(src, dstbox, srcbox, order_space);
  } else if (all(src->extent().stride() < this->extent().stride())) {
    // Restrict
    assert(transport_operator != op_sync);
    transfer_restrict(src, dstbox, srcbox, slabinfo, order_space);
  } else {
    assert(0);
  }
}

template <typename T>
void data<T>::transfer_p_vc_cc(data const *const src, ibbox const &dstbox,
                               ibbox const &srcbox, int const order_space) {
  transfer_prolongate(src, dstbox, srcbox, order_space);
}

template <>
void data<CCTK_INT>::transfer_p_vc_cc(data const *const /*src*/,
                                      ibbox const & /*dstbox*/,
                                      ibbox const & /*srcbox*/,
                                      int const /*order_space*/) {
  CCTK_ERROR("Data type not supported");
}

template <typename T>
void data<T>::transfer_prolongate(data const *const src, ibbox const &dstbox,
                                  ibbox const &srcbox, int const order_space) {
  DECLARE_CCTK_PARAMETERS;

  static Timer total("prolongate");
  total.start();

#if CARPET_DIM == 3

  switch (transport_operator) {

  case op_copy:
  case op_Lagrange: {
    static Timer timer("prolongate_Lagrange");
    timer.start();
    // enum centering { vertex_centered, cell_centered };
    switch (cent) {
    case vertex_centered: {
      static void (*the_operators[])(
          T const *restrict const src, ivect3 const &restrict srcpadext,
          ivect3 const &restrict srcext, T *restrict const dst,
          ivect3 const &restrict dstpadext, ivect3 const &restrict dstext,
          ibbox3 const &restrict srcbbox, ibbox3 const &restrict dstbbox,
          ibbox3 const &restrict srcregbbox, ibbox3 const &restrict dstregbbox,
          void *const extraargs) = {
          NULL, &prolongate_3d_rf2<T, 1>, NULL, &prolongate_3d_rf2<T, 3>,
          NULL, &prolongate_3d_rf2<T, 5>, NULL, &prolongate_3d_rf2<T, 7>,
          NULL, &prolongate_3d_rf2<T, 9>, NULL, &prolongate_3d_rf2<T, 11>,
      };
#if defined(CCTK_REAL_PRECISION_4)
      if (order_space == 11) {
        CCTK_ERROR("There is no single precision vertex-centred stencil for "
                   "op=\"LAGRANGE\" or op=\"COPY\" with order_space=11");
      }
#endif
      if (order_space < 0 or order_space > 11 or
          not the_operators[order_space]) {
        CCTK_ERROR("There is no vertex-centred stencil for op=\"LAGRANGE\" or "
                   "op=\"COPY\" with order_space not in {1, 3, 5, 7, 9, 11}");
      }
      call_operator<T>(
          the_operators[order_space], static_cast<T const *>(src->storage()),
          src->padded_shape(), src->shape(), static_cast<T *>(this->storage()),
          this->padded_shape(), this->shape(), src->extent(), this->extent(),
          srcbox, dstbox, NULL);
      break;
    }
    case cell_centered: {
      if (use_dgfe) {
        // Don't use call_operator, because we parallelise ourselves
        prolongate_3d_dgfe_rf2<T, 5>(
            static_cast<T const *>(src->storage()), src->padded_shape(),
            src->shape(), static_cast<T *>(this->storage()),
            this->padded_shape(), this->shape(), src->extent(), this->extent(),
            srcbox, dstbox, NULL);
        break;
      }
      static void (*the_operators[])(
          T const *restrict const src, ivect3 const &restrict srcpadext,
          ivect3 const &restrict srcext, T *restrict const dst,
          ivect3 const &restrict dstpadext, ivect3 const &restrict dstext,
          ibbox3 const &restrict srcbbox, ibbox3 const &restrict dstbbox,
          ibbox3 const &restrict srcregbbox, ibbox3 const &restrict dstregbbox,
          void *const extraargs) = {
          &prolongate_3d_cc_rf2<T, 0>, &prolongate_3d_cc_rf2<T, 1>,
          &prolongate_3d_cc_rf2<T, 2>, &prolongate_3d_cc_rf2<T, 3>,
          &prolongate_3d_cc_rf2<T, 4>, &prolongate_3d_cc_rf2<T, 5>};
      if (order_space < 0 or order_space > 5) {
        CCTK_ERROR("There is no cell-centred stencil for op=\"LAGRANGE\" or "
                   "op=\"COPY\" with order_space not in {0, 1, 2, 3, 4, 5}");
      }
      call_operator<T>(
          the_operators[order_space], static_cast<T const *>(src->storage()),
          src->padded_shape(), src->shape(), static_cast<T *>(this->storage()),
          this->padded_shape(), this->shape(), src->extent(), this->extent(),
          srcbox, dstbox, NULL);
      break;
    }
    default:
      assert(0);
    }
    timer.stop(0);
    break;
  }

  case op_ENO: {
    static Timer timer("prolongate_ENO");
    timer.start();
    // enum centering { vertex_centered, cell_centered };
    switch (cent) {
    case vertex_centered: {
      switch (order_space) {
      case 1:
        CCTK_ERROR("There is no stencil for op=\"ENO\" with order_space=1");
        break;
      case 3:
        call_operator<T>(
            &prolongate_3d_eno, static_cast<T const *>(src->storage()),
            src->padded_shape(), src->shape(),
            static_cast<T *>(this->storage()), this->padded_shape(),
            this->shape(), src->extent(), this->extent(), srcbox, dstbox, NULL);
        break;
      case 5:
        // There is only one parameter for the prolongation order, but
        // Whisky may want 5th order for spacetime and 3rd order for
        // hydro, so we cheat here.
        call_operator<T>(
            &prolongate_3d_eno, static_cast<T const *>(src->storage()),
            src->padded_shape(), src->shape(),
            static_cast<T *>(this->storage()), this->padded_shape(),
            this->shape(), src->extent(), this->extent(), srcbox, dstbox, NULL);
        break;
      default:
        CCTK_ERROR("There is no stencil for op=\"ENO\" with order_space!=3");
        break;
      }
      break;
    }
    case cell_centered: {
      static void (*the_operators[])(
          T const *restrict const src, ivect3 const &restrict srcpadext,
          ivect3 const &restrict srcext, T *restrict const dst,
          ivect3 const &restrict dstpadext, ivect3 const &restrict dstext,
          ibbox3 const &restrict srcbbox, ibbox3 const &restrict dstbbox,
          ibbox3 const &restrict srcregbbox, ibbox3 const &restrict dstregbbox,
          void *const extraargs) = {
          &prolongate_3d_cc_eno_rf2<T, 2>,
          &prolongate_3d_cc_eno_rf2<T, 2>, // note that we cheat here: order is
                                           // still 2 even though 3 was
                                           // requested!
          &prolongate_3d_cc_eno_rf2<T, 2>, // note that we cheat here: order is
                                           // 2 even though 4 was requested!
          &prolongate_3d_cc_eno_rf2<T, 3> // note that we cheat here: order is 3
                                          // even though 5 was requested!
          // have cheated here for two reasons: first, the ENO prolongation
          // operator stencil radius is larger than Lagrange (and dh.cc assumes
          // that the stencil goes as order_space/2!),
          // and second, we want to allow spacetime interpolation to be of
          // higher order while keeping the implemeneted ENO order!
      };
      if (order_space < 2 or order_space > 5) {
        CCTK_ERROR("There is no cell-centred stencil for op=\"ENO\" with "
                   "order_space not in {2,3,4,5}");
      }

      call_operator<T>(
          the_operators[order_space - 2],
          static_cast<T const *>(src->storage()), src->padded_shape(),
          src->shape(), static_cast<T *>(this->storage()), this->padded_shape(),
          this->shape(), src->extent(), this->extent(), srcbox, dstbox, NULL);
    } break;
    default:
      assert(0);
    }
    timer.stop(0);
  } break;

  case op_WENO: {
    static Timer timer("prolongate_WENO");
    timer.start();
    // enum centering { vertex_centered, cell_centered };
    switch (cent) {
    case vertex_centered: {
      switch (order_space) {
      case 1:
        CCTK_ERROR("There is no stencil for op=\"WENO\" with order_space=1");
        break;
      case 3:
        CCTK_ERROR("There is no stencil for op=\"WENO\" with order_space=3");
        break;
      case 5:
        call_operator<T>(
            &prolongate_3d_eno, static_cast<T const *>(src->storage()),
            src->padded_shape(), src->shape(),
            static_cast<T *>(this->storage()), this->padded_shape(),
            this->shape(), src->extent(), this->extent(), srcbox, dstbox, NULL);
        break;
      default:
        CCTK_ERROR("There is no stencil for op=\"WENO\" with order_space!=5");
        break;
      }
      break;
    }
    case cell_centered: {
      CCTK_ERROR(
          "There are currently no cell-centred stencils for op=\"WENO\"");
      break;
    }
    default:
      assert(0);
    }
    timer.stop(0);
  } break;
  case op_TVD: {
    static Timer timer("prolongate_TVD");
    timer.start();
    // enum centering { vertex_centered, cell_centered };
    switch (cent) {
    case vertex_centered: {
      switch (order_space) {
      case 1:
        call_operator<T>(
            &prolongate_3d_tvd, static_cast<T const *>(src->storage()),
            src->padded_shape(), src->shape(),
            static_cast<T *>(this->storage()), this->padded_shape(),
            this->shape(), src->extent(), this->extent(), srcbox, dstbox, NULL);
        break;
      default:
        CCTK_ERROR("There is no stencil for op=\"TVD\" with order_space!=1");
        break;
      }
      break;
    }
    case cell_centered: {
      switch (order_space) {
      case 1:
        call_operator<T>(
            &prolongate_3d_cc_tvd, static_cast<T const *>(src->storage()),
            src->padded_shape(), src->shape(),
            static_cast<T *>(this->storage()), this->padded_shape(),
            this->shape(), src->extent(), this->extent(), srcbox, dstbox, NULL);
        break;
      default:
        CCTK_ERROR("There is no stencil for op=\"TVD\" with order_space!=1");
        break;
      }
      break;
    }
    default:
      assert(0);
    }
    timer.stop(0);
    break;
  } break;
  case op_Lagrange_monotone: {
    static Timer timer("prolongate_Lagrange_monotone");
    timer.start();
    switch (order_space) {
    case 1:
      CCTK_ERROR("There is no stencil for op=\"Lagrange_monotone\" with "
                 "order_space=1");
      break;
    case 3:
      CCTK_ERROR("There is no stencil for op=\"Lagrange_monotone\" with "
                 "order_space=3");
      break;
    case 5:
      call_operator<T>(
          &prolongate_3d_o5_monotone_rf2,
          static_cast<T const *>(src->storage()), src->padded_shape(),
          src->shape(), static_cast<T *>(this->storage()), this->padded_shape(),
          this->shape(), src->extent(), this->extent(), srcbox, dstbox, NULL);
      break;
    default:
      CCTK_ERROR("There is no stencil for op=\"Lagrange_monotone\" with "
                 "order_space!=5");
      break;
    }
    timer.stop(0);
    break;
  }
  case op_STAGGER011: {
    static void (*the_operators[])(
        T const *restrict const src, ivect3 const &restrict srcpadext,
        ivect3 const &restrict srcext, T *restrict const dst,
        ivect3 const &restrict dstpadext, ivect3 const &restrict dstext,
        ibbox3 const &restrict srcbbox, ibbox3 const &restrict dstbbox,
        ibbox3 const &restrict srcregbbox, ibbox3 const &restrict dstregbbox,
        void *const extraargs) = {
        &prolongate_3d_stagger011<T, 2>, &prolongate_3d_stagger011<T, 3>,
        &prolongate_3d_stagger011<T, 3>, &prolongate_3d_stagger011<T, 3>};

    static Timer timer("prolongate_STAGGER011");
    timer.start();
    call_operator<T>(
        the_operators[order_space - 2], static_cast<T const *>(src->storage()),
        src->padded_shape(), src->shape(), static_cast<T *>(this->storage()),
        this->padded_shape(), this->shape(), src->extent(), this->extent(),
        srcbox, dstbox, NULL);
    timer.stop(0);
    break;
  }

  case op_STAGGER101: {
    static void (*the_operators[])(
        T const *restrict const src, ivect3 const &restrict srcpadext,
        ivect3 const &restrict srcext, T *restrict const dst,
        ivect3 const &restrict dstpadext, ivect3 const &restrict dstext,
        ibbox3 const &restrict srcbbox, ibbox3 const &restrict dstbbox,
        ibbox3 const &restrict srcregbbox, ibbox3 const &restrict dstregbbox,
        void *const extraargs) = {
        &prolongate_3d_stagger101<T, 2>, &prolongate_3d_stagger101<T, 3>,
        &prolongate_3d_stagger101<T, 3>, &prolongate_3d_stagger101<T, 3>};

    static Timer timer("prolongate_STAGGER101");
    timer.start();
    call_operator<T>(
        the_operators[order_space - 2], static_cast<T const *>(src->storage()),
        src->padded_shape(), src->shape(), static_cast<T *>(this->storage()),
        this->padded_shape(), this->shape(), src->extent(), this->extent(),
        srcbox, dstbox, NULL);
    timer.stop(0);
    break;
  }

  case op_STAGGER110: {
    static void (*the_operators[])(
        T const *restrict const src, ivect3 const &restrict srcpadext,
        ivect3 const &restrict srcext, T *restrict const dst,
        ivect3 const &restrict dstpadext, ivect3 const &restrict dstext,
        ibbox3 const &restrict srcbbox, ibbox3 const &restrict dstbbox,
        ibbox3 const &restrict srcregbbox, ibbox3 const &restrict dstregbbox,
        void *const extraargs) = {
        &prolongate_3d_stagger110<T, 2>, &prolongate_3d_stagger110<T, 3>,
        &prolongate_3d_stagger110<T, 3>, &prolongate_3d_stagger110<T, 3>};

    static Timer timer("prolongate_STAGGER110");
    timer.start();
    call_operator<T>(
        the_operators[order_space - 2], static_cast<T const *>(src->storage()),
        src->padded_shape(), src->shape(), static_cast<T *>(this->storage()),
        this->padded_shape(), this->shape(), src->extent(), this->extent(),
        srcbox, dstbox, NULL);
    timer.stop(0);
    break;
  }

  case op_STAGGER111: {
    static void (*the_operators[])(
        T const *restrict const src, ivect3 const &restrict srcpadext,
        ivect3 const &restrict srcext, T *restrict const dst,
        ivect3 const &restrict dstpadext, ivect3 const &restrict dstext,
        ibbox3 const &restrict srcbbox, ibbox3 const &restrict dstbbox,
        ibbox3 const &restrict srcregbbox, ibbox3 const &restrict dstregbbox,
        void *const extraargs) = {
        &prolongate_3d_stagger111<T, 2>, &prolongate_3d_stagger111<T, 3>,
        &prolongate_3d_stagger111<T, 3>, &prolongate_3d_stagger111<T, 3>};

    static Timer timer("prolongate_STAGGER111");
    timer.start();
    call_operator<T>(
        the_operators[order_space - 2], static_cast<T const *>(src->storage()),
        src->padded_shape(), src->shape(), static_cast<T *>(this->storage()),
        this->padded_shape(), this->shape(), src->extent(), this->extent(),
        srcbox, dstbox, NULL);
    timer.stop(0);
    break;
  }

  default:
    assert(0);
  } // switch (transport_operator)

#elif CARPET_DIM == 4

  switch (transport_operator) {

  case op_copy:
  case op_Lagrange: {
    static Timer timer("prolongate_Lagrange");
    timer.start();
    // enum centering { vertex_centered, cell_centered };
    switch (cent) {
    case vertex_centered:
      switch (order_space) {
      case 1:
        call_operator<T>(
            &prolongate_4d_o1_rf2, static_cast<T const *>(src->storage()),
            src->padded_shape(), src->shape(),
            static_cast<T *>(this->storage()), this->padded_shape(),
            this->shape(), src->extent(), this->extent(), srcbox, dstbox, NULL);
        break;
      default:
        CCTK_ERROR("There is no vertex-centred stencil for op=\"LAGRANGE\" "
                   "with order_space not in {1}");
        break;
      }
      break;
    default:
      assert(0);
    }
    timer.stop(0);
    break;
  }
  default:
    assert(0);
  } // switch (transport_operator)

#else
#error "Value for CARPET_DIM not supported"
#endif

  total.stop(0);
}

template <>
void data<CCTK_INT>::transfer_prolongate(data const *const /*src*/,
                                         ibbox const & /*dstbox*/,
                                         ibbox const & /*srcbox*/,
                                         int const /*order_space*/) {
  CCTK_ERROR("Data type not supported");
}

template <typename T>
void data<T>::transfer_restrict(data const *const src, ibbox const &dstregbox,
                                ibbox const &srcregbox,
                                islab const *restrict const slabinfo,
                                int const /*order_space*/) {
  DECLARE_CCTK_PARAMETERS;

  static Timer total("restrict");
  total.start();

#if CARPET_DIM == 3

  switch (transport_operator) {

  case op_STAGGER011:
    assert(not slabinfo);
    call_operator<T>(
        &restrict_3d_stagger011, static_cast<T const *>(src->storage()),
        src->padded_shape(), src->shape(), static_cast<T *>(this->storage()),
        this->padded_shape(), this->shape(), src->extent(), this->extent(),
        srcregbox, dstregbox, NULL);
    break;

  case op_STAGGER101:
    assert(not slabinfo);
    call_operator<T>(
        &restrict_3d_stagger101, static_cast<T const *>(src->storage()),
        src->padded_shape(), src->shape(), static_cast<T *>(this->storage()),
        this->padded_shape(), this->shape(), src->extent(), this->extent(),
        srcregbox, dstregbox, NULL);
    break;

  case op_STAGGER110:
    assert(not slabinfo);
    call_operator<T>(
        &restrict_3d_stagger110, static_cast<T const *>(src->storage()),
        src->padded_shape(), src->shape(), static_cast<T *>(this->storage()),
        this->padded_shape(), this->shape(), src->extent(), this->extent(),
        srcregbox, dstregbox, NULL);
    break;

  case op_STAGGER111:
    assert(not slabinfo);
    call_operator<T>(
        &restrict_3d_stagger111, static_cast<T const *>(src->storage()),
        src->padded_shape(), src->shape(), static_cast<T *>(this->storage()),
        this->padded_shape(), this->shape(), src->extent(), this->extent(),
        srcregbox, dstregbox, NULL);
    break;

  case op_copy:
  case op_Lagrange:
  case op_ENO:
  case op_WENO:
  case op_TVD:
  case op_Lagrange_monotone:
  case op_restrict:
    // enum centering { vertex_centered, cell_centered };
    switch (cent) {
    case vertex_centered:
      assert(not slabinfo);
      call_operator<T>(&restrict_3d_rf2, static_cast<T const *>(src->storage()),
                       src->padded_shape(), src->shape(),
                       static_cast<T *>(this->storage()), this->padded_shape(),
                       this->shape(), src->extent(), this->extent(), srcregbox,
                       dstregbox, NULL);
      break;
    case cell_centered: {
      assert(all(dstregbox.stride() == this->extent().stride()));
      ivect const is_centered = slabinfo ? slabinfo->is_centered : 1;

      ibbox const &srcbox = src->extent();
      ibbox const &dstbox = this->extent();

      if (all(is_centered == ivect(1, 1, 1))) {
        if (use_dgfe) {
          // Don't use call_operator, because we parallelise ourselves
          restrict_3d_dgfe_rf2<T, 5>(
              static_cast<T const *>(src->storage()), src->padded_shape(),
              src->shape(), static_cast<T *>(this->storage()),
              this->padded_shape(), this->shape(), srcbox, dstbox, srcregbox,
              dstregbox, NULL);
          break;
        }
        if (use_higher_order_restriction and transport_operator != op_WENO and
            transport_operator != op_ENO) { // HACK
          switch (restriction_order_space) {
          case 1:
            // Don't use call_operator, because we parallelise ourselves
            restrict_3d_cc_rf2(static_cast<T const *>(src->storage()),
                               src->padded_shape(), src->shape(),
                               static_cast<T *>(this->storage()),
                               this->padded_shape(), this->shape(), srcbox,
                               dstbox, srcregbox, dstregbox, NULL);
            break;
          case 3:
            // Don't use call_operator, because we parallelise ourselves
            restrict_3d_cc_o3_rf2(static_cast<T const *>(src->storage()),
                                  src->padded_shape(), src->shape(),
                                  static_cast<T *>(this->storage()),
                                  this->padded_shape(), this->shape(), srcbox,
                                  dstbox, srcregbox, dstregbox, NULL);
            break;
          case 5:
            // Don't use call_operator, because we parallelise ourselves
            restrict_3d_cc_o5_rf2(static_cast<T const *>(src->storage()),
                                  src->padded_shape(), src->shape(),
                                  static_cast<T *>(this->storage()),
                                  this->padded_shape(), this->shape(), srcbox,
                                  dstbox, srcregbox, dstregbox, NULL);
            break;
          default:
            CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                        "There is no restriction stencil with "
                        "restriction_order_space==%d",
                        int(restriction_order_space));
            break;
          }
          break;
        }
        // Don't use call_operator, because we parallelise ourselves
        restrict_3d_cc_rf2(static_cast<T const *>(src->storage()),
                           src->padded_shape(), src->shape(),
                           static_cast<T *>(this->storage()),
                           this->padded_shape(), this->shape(), srcbox, dstbox,
                           srcregbox, dstregbox, NULL);
      } else if (all(is_centered == ivect(0, 1, 1))) {
        // Don't use call_operator, because we parallelise ourselves
        restrict_3d_vc_rf2<T, 0, 1, 1>(
            static_cast<T const *>(src->storage()), src->padded_shape(),
            src->shape(), static_cast<T *>(this->storage()),
            this->padded_shape(), this->shape(), srcbox, dstbox, srcregbox,
            dstregbox, NULL);
      } else if (all(is_centered == ivect(1, 0, 1))) {
        // Don't use call_operator, because we parallelise ourselves
        restrict_3d_vc_rf2<T, 1, 0, 1>(
            static_cast<T const *>(src->storage()), src->padded_shape(),
            src->shape(), static_cast<T *>(this->storage()),
            this->padded_shape(), this->shape(), srcbox, dstbox, srcregbox,
            dstregbox, NULL);
      } else if (all(is_centered == ivect(1, 1, 0))) {
        // Don't use call_operator, because we parallelise ourselves
        restrict_3d_vc_rf2<T, 1, 1, 0>(
            static_cast<T const *>(src->storage()), src->padded_shape(),
            src->shape(), static_cast<T *>(this->storage()),
            this->padded_shape(), this->shape(), srcbox, dstbox, srcregbox,
            dstregbox, NULL);
      } else {
        assert(0);
      }
      break;
    }
    default:
      assert(0);
    }
    break;

  default:
    assert(0);
  }

#elif CARPET_DIM == 4

  switch (transport_operator) {

  case op_copy:
  case op_Lagrange:
    // enum centering { vertex_centered, cell_centered };
    switch (cent) {
    case vertex_centered:
      // Don't use call_operator, because we parallelise ourselves
      restrict_4d_rf2(static_cast<T const *>(src->storage()),
                      src->padded_shape(), src->shape(),
                      static_cast<T *>(this->storage()), this->padded_shape(),
                      this->shape(), src->extent(), this->extent(), srcregbox,
                      dstregbox, NULL);
      break;
    default:
      assert(0);
    }
    break;

  default:
    assert(0);
  }

#else
#error "Value for CARPET_DIM not supported"
#endif

  total.stop(0);
}

template <>
void data<CCTK_INT>::transfer_restrict(data const *const /*src*/,
                                       ibbox const & /*dstbox*/,
                                       ibbox const & /*srcbox*/,
                                       islab const *restrict const /*slabinfo*/,
                                       int const /*order_space*/) {
  CCTK_ERROR("Data type not supported");
}

template <typename T>
void data<T>::time_interpolate(vector<data *> const &srcs, ibbox const &dstbox,
                               ibbox const &srcbox,
                               vector<CCTK_REAL> const &times,
                               CCTK_REAL const time, int const order_time) {
  static Timer total("time_interpolate");
  total.start();

#if CARPET_DIM == 3

  switch (transport_operator) {

  case op_copy:
  case op_STAGGER011:
  case op_STAGGER101:
  case op_STAGGER110:
  case op_STAGGER111:
  case op_Lagrange: {
    static Timer timer("time_interpolate_Lagrange");
    timer.start();
    switch (order_time) {

    case 1:
      assert(times.size() >= 2);
      interpolate_3d_2tl(
          static_cast<T const *>(srcs.AT(0)->storage()), times.AT(0),
          static_cast<T const *>(srcs.AT(1)->storage()), times.AT(1),
          srcs.AT(0)->padded_shape(), srcs.AT(0)->shape(),
          static_cast<T *>(this->storage()), time, this->padded_shape(),
          this->shape(), srcs.AT(0)->extent(), this->extent(), srcbox, dstbox,
          NULL);
      break;

    case 2:
      assert(times.size() >= 3);
      interpolate_3d_3tl(
          static_cast<T const *>(srcs.AT(0)->storage()), times.AT(0),
          static_cast<T const *>(srcs.AT(1)->storage()), times.AT(1),
          static_cast<T const *>(srcs.AT(2)->storage()), times.AT(2),
          srcs.AT(0)->padded_shape(), srcs.AT(0)->shape(),
          static_cast<T *>(this->storage()), time, this->padded_shape(),
          this->shape(), srcs.AT(0)->extent(), this->extent(), srcbox, dstbox,
          NULL);
      break;

    case 3:
      assert(times.size() >= 4);
      interpolate_3d_4tl(
          static_cast<T const *>(srcs.AT(0)->storage()), times.AT(0),
          static_cast<T const *>(srcs.AT(1)->storage()), times.AT(1),
          static_cast<T const *>(srcs.AT(2)->storage()), times.AT(2),
          static_cast<T const *>(srcs.AT(3)->storage()), times.AT(3),
          srcs.AT(0)->padded_shape(), srcs.AT(0)->shape(),
          static_cast<T *>(this->storage()), time, this->padded_shape(),
          this->shape(), srcs.AT(0)->extent(), this->extent(), srcbox, dstbox,
          NULL);
      break;

    case 4:
      assert(times.size() >= 5);
      interpolate_3d_5tl(
          static_cast<T const *>(srcs.AT(0)->storage()), times.AT(0),
          static_cast<T const *>(srcs.AT(1)->storage()), times.AT(1),
          static_cast<T const *>(srcs.AT(2)->storage()), times.AT(2),
          static_cast<T const *>(srcs.AT(3)->storage()), times.AT(3),
          static_cast<T const *>(srcs.AT(4)->storage()), times.AT(4),
          srcs.AT(0)->padded_shape(), srcs.AT(0)->shape(),
          static_cast<T *>(this->storage()), time, this->padded_shape(),
          this->shape(), srcs.AT(0)->extent(), this->extent(), srcbox, dstbox,
          NULL);
      break;

    default:
      assert(0);
    }
    timer.stop(0);
    break;
  }

  case op_ENO:
  case op_WENO:
  case op_TVD:
  case op_Lagrange_monotone: {
    // ENO, WENO, TVD, and Lagrange_monotone time interpolation is the
    // same for order_time <= 2
    static Timer timer("time_interpolate_ENO");
    timer.start();
    switch (order_time) {

    case 1:
      assert(times.size() >= 2);
      interpolate_3d_2tl(
          static_cast<T const *>(srcs.AT(0)->storage()), times.AT(0),
          static_cast<T const *>(srcs.AT(1)->storage()), times.AT(1),
          srcs.AT(0)->padded_shape(), srcs.AT(0)->shape(),
          static_cast<T *>(this->storage()), time, this->padded_shape(),
          this->shape(), srcs.AT(0)->extent(), this->extent(), srcbox, dstbox,
          NULL);
      break;

    case 2:
      assert(times.size() >= 3);
      interpolate_eno_3d_3tl(
          static_cast<T const *>(srcs.AT(0)->storage()), times.AT(0),
          static_cast<T const *>(srcs.AT(1)->storage()), times.AT(1),
          static_cast<T const *>(srcs.AT(2)->storage()), times.AT(2),
          srcs.AT(0)->padded_shape(), srcs.AT(0)->shape(),
          static_cast<T *>(this->storage()), time, this->padded_shape(),
          this->shape(), srcs.AT(0)->extent(), this->extent(), srcbox, dstbox,
          NULL);
      break;

    default:
      assert(0);
    }
    timer.stop(0);
    break;
  }

  default:
    assert(0);
  } // switch (transport_operator)

#elif CARPET_DIM == 4

  assert(0);

#else
#error "Value for CARPET_DIM not supported"
#endif

  total.stop(0);
}

template <>
void data<CCTK_INT>::time_interpolate(vector<data *> const & /*srcs*/,
                                      ibbox const & /*dstbox*/,
                                      ibbox const & /*srcbox*/,
                                      vector<CCTK_REAL> const & /*times*/,
                                      CCTK_REAL const /*time*/,
                                      int const /*order_time*/) {
  CCTK_ERROR("Data type not supported");
}

// Memory usage
template <typename T> size_t data<T>::memory() const {
  return memoryof(_memory) + memoryof(vectorlength) + memoryof(vectorindex) +
         memoryof(vectorleader);
}

// Output
template <typename T> ostream &data<T>::output(ostream &os) const {
  T Tdummy;
  os << "data<" << typestring(Tdummy) << ">:"
     << "extent=" << extent() << ","
     << "shape=" << shape() << ","
     << "padded_shape=" << padded_shape() << ","
     << "stride=" << stride() << ",size=" << size();
  return os;
}

#define TYPECASE(N, T) template class data<T>;
#include "typecase.hh"
#undef TYPECASE
