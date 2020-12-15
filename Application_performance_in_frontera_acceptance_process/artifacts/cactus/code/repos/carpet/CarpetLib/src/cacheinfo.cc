#include "cacheinfo.hh"

#include <iostream>
#include <sstream>

#include <cctk.h>
#include <cctk_Parameters.h>

#include <vectors.h>

template <int D> vect<int, D> pad_shape(bbox<int, D> const &extent) {
  assert(all(extent.shape() >= 0));
  return pad_shape(extent.shape() / extent.stride());
}

namespace {
struct cache_info_t {
  int type;
  int linesize;
  int stride;
};
bool have_cache_info = false;
vector<cache_info_t> cache_info;
}

template <int D> vect<int, D> pad_shape(vect<int, D> const &shape) {
  DECLARE_CCTK_PARAMETERS;

  assert(all(shape >= 0));

  // Don't pad empty arrays; we don't want to handle all the special
  // cases for this below
  if (any(shape == 0))
    return shape;

  if (CCTK_BUILTIN_EXPECT(not have_cache_info, false)) {
#pragma omp barrier
#pragma omp master
    {
      if (CCTK_IsFunctionAliased("GetCacheInfo1")) {
        int const num_levels =
            GetCacheInfo1(NULL, NULL, NULL, NULL, NULL, NULL, 0);
        vector<CCTK_INT> types(num_levels);
        vector<CCTK_INT> linesizes(num_levels);
        vector<CCTK_INT> strides(num_levels);
        GetCacheInfo1(NULL, &types[0], NULL, &linesizes[0], &strides[0], NULL,
                      num_levels);
        cache_info.resize(num_levels);
        for (int level = 0; level < num_levels; ++level) {
          cache_info[level].type = types[level];
          cache_info[level].linesize = linesizes[level];
          cache_info[level].stride = strides[level];
        }
      }
      have_cache_info = true;
    }
#pragma omp barrier
  }

  vect<int, D> padded_shape;
  size_type accumulated_npoints = 1;
  for (int d = 0; d < D; ++d) {
    size_type npoints = shape[d];

#if VECTORISE && VECTORISE_ALIGNED_ARRAYS
    if (d == 0) {
      // Pad array to a multiple of the vector size. Note that this is
      // a hard requirement, so that we can emit aligned load/store
      // operations.
      npoints = align_up(npoints, size_type(CCTK_REAL_VEC_SIZE));
    }

    if (pad_to_cachelines) {
      for (size_t cache_level = 0; cache_level < cache_info.size();
           ++cache_level) {
        if (cache_info[cache_level].type == 0) {
          // Pad array in this direction to a multiple of this cache
          // line size
          int const cache_linesize = cache_info[cache_level].linesize;
          int const cache_stride = cache_info[cache_level].stride;

          assert(cache_linesize % sizeof(CCTK_REAL) == 0);
          int const linesize = cache_linesize / sizeof(CCTK_REAL);
          if (npoints * accumulated_npoints < linesize) {
            // The extent is less than one cache line long: Ensure
            // that the array size divides the cache line size evenly
            // by rounding to the next power of 2
            assert(is_power_of_2(linesize));
            npoints = next_power_of_2(npoints);
          } else {
            // The extent is at least one cache line long: round up to
            // the next full cache line
            size_type total_npoints = npoints * accumulated_npoints;
            total_npoints = align_up(total_npoints, size_type(linesize));
            assert(total_npoints % accumulated_npoints == 0);
            npoints = total_npoints / accumulated_npoints;
          }

          // Avoid multiples of the cache stride
          if (cache_stride > 0) {
            assert(cache_stride % sizeof(CCTK_REAL) == 0);
            int const stride = cache_stride / sizeof(CCTK_REAL);
            if (npoints * accumulated_npoints % stride == 0) {
              assert(stride > linesize);
              size_type total_npoints = npoints * accumulated_npoints;
              total_npoints += max(size_type(linesize), accumulated_npoints);
              assert(total_npoints % accumulated_npoints == 0);
              npoints = total_npoints / accumulated_npoints;
            }
          }
        } // if is cache
      }   // for cache_level
    }     // if pad_to_cachelines
#endif

    padded_shape[d] = npoints;
    accumulated_npoints *= npoints;
  }
  assert(prod(vect<size_type, D>(padded_shape)) == accumulated_npoints);

  // self-check
  for (int d = 0; d < D; ++d) {
    assert(padded_shape[d] >= shape[d]);
#if VECTORISE && VECTORISE_ALIGNED_ARRAYS
    if (d == 0) {
      assert(padded_shape[d] % CCTK_REAL_VEC_SIZE == 0);
    }
#endif
  }

  // Safety check
  if (not(prod(vect<size_type, D>(padded_shape)) <=
          2 * prod(vect<size_type, D>(shape)) + 1000)) {
    cerr << "shape=" << shape << "   "
         << "prod(shape)=" << prod(vect<size_type, D>(shape)) << "\n"
         << "padded_shape=" << padded_shape << "   "
         << "prod(padded_shape)=" << prod(vect<size_type, D>(padded_shape))
         << "\n";
  }
  assert(prod(vect<size_type, D>(padded_shape)) <=
         2 * prod(vect<size_type, D>(shape)) + 1000);

  if (verbose) {
    ostringstream buf;
    buf << "padding " << shape << " to " << padded_shape;
    CCTK_INFO(buf.str().c_str());
  }

  return padded_shape;
}

template vect<int, 3> pad_shape(bbox<int, 3> const &extent);
template vect<int, 3> pad_shape(vect<int, 3> const &shape);

template vect<int, 4> pad_shape(bbox<int, 4> const &extent);
template vect<int, 4> pad_shape(vect<int, 4> const &shape);
