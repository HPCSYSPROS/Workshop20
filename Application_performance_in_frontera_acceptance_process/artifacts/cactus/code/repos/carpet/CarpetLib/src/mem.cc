#include <cctk.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include <vectors.h>

#include "defs.hh"
#include "mem.hh"

using namespace std;

double const gmem::KILO = 1000.0;
double const gmem::MEGA = 1000.0 * 1000.0;
double const gmem::GIGA = 1000.0 * 1000.0 * 1000.0;
double const gmem::TERA = 1000.0 * 1000.0 * 1000.0 * 1000.0;
double const gmem::PETA = 1000.0 * 1000.0 * 1000.0 * 1000.0 * 1000.0;
double const gmem::EXA = 1000.0 * 1000.0 * 1000.0 * 1000.0 * 1000.0 * 1000.0;

// Total number of currently allocated bytes and objects
double gmem::total_allocated_bytes = 0;
double gmem::total_allocated_objects = 0;

// Maximum of the above (over time)
double gmem::max_allocated_bytes = 0;
double gmem::max_allocated_objects = 0;

namespace {
size_t get_max_cache_linesize() {
  static size_t max_cache_linesize = 0;
  if (CCTK_BUILTIN_EXPECT(max_cache_linesize == 0, false)) {
#pragma omp barrier
#pragma omp master
    {
      max_cache_linesize = 1;
      if (CCTK_IsFunctionAliased("GetCacheInfo1")) {
        int const num_levels =
            GetCacheInfo1(NULL, NULL, NULL, NULL, NULL, NULL, 0);
        vector<CCTK_INT> types(num_levels);
        vector<CCTK_INT> linesizes(num_levels);
        vector<CCTK_INT> strides(num_levels);
        GetCacheInfo1(NULL, &types[0], NULL, &linesizes[0], &strides[0], NULL,
                      num_levels);
        for (int level = 0; level < num_levels; ++level) {
          if (types[level] == 0) { // if this is a cache
            max_cache_linesize =
                max(max_cache_linesize, size_t(linesizes[level]));
          }
        }
      }
    }
#pragma omp barrier
  }
  assert(max_cache_linesize > 0);
  return max_cache_linesize;
}
}

// TODO: Make this a plain class instead of a template

template <typename T>
mem<T>::mem(size_t const vectorlength, size_t const nelems, T *const memptr,
            size_t const memsize)
    : storage_base_(memptr), storage_(memptr), nelems_(nelems),
      vectorlength_(vectorlength), owns_storage_(false),
      clients_(vectorlength, false), num_clients_(0) {
  DECLARE_CCTK_PARAMETERS;
  if (memptr == NULL) {

#if VECTORISE
    size_t const vector_size = CCTK_REAL_VEC_SIZE;
#else
    size_t const vector_size = 1;
#endif
    size_t const canary = electric_fence ? 2 * fence_width : 0;
    size_t const final_padding = vector_size - 1;
    size_t const max_cache_linesize = get_max_cache_linesize();
    size_t const alignment =
        align_up(max_cache_linesize, vector_size * sizeof(T));
    assert(alignment >= 1);
    // Safety check
    assert(alignment <= 1024);

    const size_t nbytes =
        (vectorlength * nelems + canary + final_padding) * sizeof(T);
    if (max_allowed_memory_MB > 0 and
        (total_allocated_bytes + nbytes > MEGA * max_allowed_memory_MB)) {
      T Tdummy;
      CCTK_VError(
          __LINE__, __FILE__, CCTK_THORNSTRING,
          "Refusing to allocate %.0f bytes (%.3f MB) of memory for type %s.  "
          "%.0f bytes (%.3f MB) are currently allocated in %d objects.  The "
          "parameter file specifies a maximum of %d MB",
          double(nbytes), double(nbytes / MEGA), typestring(Tdummy),
          double(total_allocated_bytes), double(total_allocated_bytes / MEGA),
          int(total_allocated_objects), int(max_allowed_memory_MB));
    }

    // void* ptr;
    // const int ierr = posix_memalign(&ptr, alignment, nbytes);
    void *ptr = malloc(nbytes + alignment - 1);
    if (not ptr) {
      T Tdummy;
      CCTK_VError(
          __LINE__, __FILE__, CCTK_THORNSTRING,
          "Failed to allocate %.0f bytes (%.3f MB) of memory for type %s.  "
          "%.0f bytes (%.3f MB) are currently allocated in %d objects",
          double(nbytes), double(nbytes / MEGA), typestring(Tdummy),
          double(total_allocated_bytes), double(total_allocated_bytes / MEGA),
          int(total_allocated_objects));
    }

    storage_base_ = (T *)ptr;
    storage_ = (T *)align_up(size_t(storage_base_ + canary / 2), alignment);
    assert(size_t(storage_) % alignment == 0);
    owns_storage_ = true;

    total_allocated_bytes += nbytes;
    max_allocated_bytes = max(max_allocated_bytes, total_allocated_bytes);
    if (poison_new_memory) {
      memset(storage_, poison_value, nbytes);
    }
    if (electric_fence) {
      // put poison just before and just after payload region
      // FIXME: this will not work with alignment for vectorizing. Not sure how
      // to support that as well as protect grid scalars.
      memset(storage_ - fence_width, poison_value, fence_width * sizeof(T));
      memset(storage_ + vectorlength_ * nelems_, poison_value,
             fence_width * sizeof(T));
    }

  } else {
    assert(memsize >= vectorlength * nelems * sizeof(T));
    // Don't poison the memory.  Passing in a pointer allows the
    // pointer to be re-interpreted as a mem object, keeping the
    // previous content.  This is e.g. used to turn communication
    // buffers into mem objects.
  }
  ++total_allocated_objects;
  max_allocated_objects = max(max_allocated_objects, total_allocated_objects);
}

template <typename T> mem<T>::~mem() {
  DECLARE_CCTK_PARAMETERS;
  assert(not has_clients());
  if (owns_storage_) {
    // do we really want this automated check in the destructor?
    // what if we are already terminating to to a failed fence check?
    if (electric_fence)
      assert(is_fence_intact(0) && is_fence_intact(1));
    free(storage_base_);

#if VECTORISE
    size_t const vector_size = CCTK_REAL_VEC_SIZE;
#else
    size_t const vector_size = 1;
#endif
    size_t const canary = electric_fence ? 2 * fence_width : 0;
    size_t const final_padding = vector_size - 1;

    const size_t nbytes =
        (vectorlength_ * nelems_ + canary + final_padding) * sizeof(T);
    total_allocated_bytes -= nbytes;
  }
  --total_allocated_objects;
}

template <typename T> bool mem<T>::is_fence_intact(const int upperlower) const {
  DECLARE_CCTK_PARAMETERS;
  bool retval = true;

  if (owns_storage_) {
    assert(storage_ and storage_base_);
    if (electric_fence) {
      T worm;
      memset(&worm, poison_value, sizeof(T));
      if (upperlower) {
        for (int i = 0; i < fence_width; ++i) {
          retval =
              retval && (memcmp(&worm, storage_ + vectorlength_ * nelems_ + i,
                                sizeof(T)) == 0);
        }
      } else {
        for (int i = 0; i < fence_width; ++i) {
          retval = retval && (memcmp(&worm, storage_ - 1 - i, sizeof(T)) == 0);
        }
      }
    }
  }

  return retval;
}

template <typename T> void mem<T>::register_client(size_t const vectorindex) {
  assert(vectorindex < vectorlength_);
  assert(not clients_.AT(vectorindex));
  clients_.AT(vectorindex) = true;
  ++num_clients_;
}

template <typename T> void mem<T>::unregister_client(size_t const vectorindex) {
  assert(vectorindex < vectorlength_);
  assert(clients_.AT(vectorindex));
  clients_.AT(vectorindex) = false;
  assert(num_clients_ > 0);
  --num_clients_;
}

template <typename T> bool mem<T>::has_clients() const {
  // return find (clients_.begin(), clients_.end(), true) != clients_.end();
  return num_clients_ > 0;
}

// Memory usage
template <typename T> size_t mem<T>::memory() const {
  return memoryof(storage_base_) + memoryof(storage_) + memoryof(nelems_) +
         memoryof(vectorlength_) + memoryof(owns_storage_) +
         memoryof(clients_) + memoryof(num_clients_) +
         (owns_storage_ ? (vectorlength_ * nelems_ + storage_ - storage_base_)
                        : 0) *
             sizeof(T);
}

size_t const mempool::chunksize;
size_t const mempool::align;

mempool::mempool() : allocated(0), freeptr(0), freesize(0) {}

mempool::~mempool() {
  while (not chunks.empty()) {
    free(chunks.top());
    chunks.pop();
  }
}

// TODO: add electric fence
void *mempool::alloc(size_t nbytes) {
  // Take a shortcut for silly requests
  if (nbytes == 0)
    return 0;

  // Round up request size
  nbytes = (nbytes + align - 1) / align * align;

  // If there is not enough memory left, allocate a new chunk.  Ignore
  // whatever is left in the old chunk.
  if (nbytes > freesize) {
    // Allocate the usual chunk size, or more if more is requested
    freesize = max(chunksize, nbytes);
    freeptr = malloc(freesize);
    allocated += freesize;
    if (not freeptr) {
      CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Failed to allocate %.3f MB of memory",
                 double(freesize / gmem::MEGA));
    }
    // Remember the pointer so that it can be freed
    chunks.push(freeptr);
  }

  // Allocate a piece from the current chunk
  void *const ptr = freeptr;
  assert(freesize >= nbytes);
  freesize -= nbytes;
  assert(freeptr);
  freeptr = static_cast<char *>(freeptr) + nbytes;

  return ptr;
}

// Memory usage
size_t mempool::memory() const {
  return memoryof(chunks) + memoryof(freeptr) + memoryof(freesize) +
         memoryof(allocated);
}

#define TYPECASE(N, T) template class mem<T>;

#include "typecase.hh"

#undef TYPECASE
