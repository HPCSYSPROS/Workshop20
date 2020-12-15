#ifndef MEM_HH
#define MEM_HH

#include <cctk.h>

#include <cassert>
#include <cstdlib>
#include <stack>
#include <vector>

#include "defs.hh"

using namespace std;

// A chunk of memory, possibly shared between some clients
class gmem {
public:
  static double const KILO;
  static double const MEGA;
  static double const GIGA;
  static double const TERA;
  static double const PETA;
  static double const EXA;

  // Total number of currently allocated bytes and objects
  static double total_allocated_bytes;
  static double total_allocated_objects;

  // Maximum of the above (over time)
  static double max_allocated_bytes;
  static double max_allocated_objects;
};

template <typename T> class mem : public gmem {
  T *storage_base_;
  T *storage_;
  size_t nelems_;
  size_t vectorlength_;
  bool owns_storage_;

  vector<bool> clients_;
  size_t num_clients_;

public:
  mem(size_t vectorlength, size_t nelems, T *memptr = NULL, size_t memsize = 0);
  ~mem();

  T *storage(size_t vectorindex) const {
    assert(vectorindex < vectorlength_);
    assert(clients_.AT(vectorindex));
    return &storage_[vectorindex * nelems_];
  }

  // return true if fence is intact, ie poison value is still there
  bool is_fence_intact(const int upperlower) const;

  void register_client(size_t vectorindex);
  void unregister_client(size_t vectorindex);
  bool has_clients() const CCTK_MEMBER_ATTRIBUTE_PURE;

  // Memory usage
  size_t memory() const CCTK_MEMBER_ATTRIBUTE_PURE;
};

template <typename T> inline size_t memoryof(mem<T> const &m) {
  return m.memory();
}

// A mempool (memory pool) is a large chunk of memory.  You can
// allocate pieces of it.  In order to simplify things there is no way
// to free a piece again.  If the mempool is destroyed, then all its
// memory is freed.  This is dangerous: you have to make sure that no
// one continues to use that memory afterwards.  Using a memory pool
// for short-lived objects can reduce memory fragmentation.
class mempool {
  // The minimum chunk size which is requested via malloc.  If a
  // larger piece is required, then a larger chunk is allocated.
  static size_t const chunksize = 10 * 1024 * 1024;
  // The alignment of the returned memory.  All requests are rounded
  // up to the next multiple of this alignment.
  static size_t const align = 32;

  // List of all allocated chunks.  When the mempool is destroyed,
  // these pointers need to be freed.
  stack<void *> chunks;
  // Total size of all allocated chunks
  size_t allocated;

  // Pointer to the beginning of some unused memory
  void *freeptr;
  // Size of that unused memory
  size_t freesize;

private:
  // Forbid copying
  mempool(mempool const &);

public:
  // Create and destroy a memory pool
  mempool();
  ~mempool();

  // Allocate some memory and return a pointer to it.  This cannot
  // fail.
  void *alloc(size_t nbytes);

  // Memory usage
  size_t memory() const CCTK_MEMBER_ATTRIBUTE_PURE;
};

inline size_t memoryof(mempool const &m) { return m.memory(); }

#endif // ifndef MEM_HH
