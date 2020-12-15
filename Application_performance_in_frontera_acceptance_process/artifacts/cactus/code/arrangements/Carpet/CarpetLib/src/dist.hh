#ifndef DIST_HH
#define DIST_HH

#include <cctk.h>

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <iostream>

#ifdef CCTK_MPI
#include <mpi.h>
#else
#include "nompi.h"
#endif
#ifdef _OPENMP
#include <omp.h>
#else
static inline int omp_get_max_threads() { return 1; }
static inline int omp_get_thread_num() { return 0; }
static inline int omp_get_num_threads() { return 1; }
#endif

#include "defs.hh"

using namespace std;

namespace dist {

extern MPI_Comm comm_;

extern MPI_Datatype mpi_complex8;
extern MPI_Datatype mpi_complex16;
extern MPI_Datatype mpi_complex32;
extern MPI_Op mpi_max;
extern MPI_Op mpi_min;
extern MPI_Op mpi_prod;
extern MPI_Op mpi_sum;

extern int total_num_threads_;

void init(int &argc, char **&argv);
void pseudoinit(MPI_Comm const c);
void finalize();
void barrier(MPI_Comm const c, int const id, char const *const errormsg);

// Create MPI datatypes from C structures

struct mpi_struct_descr_t {
  int blocklength;
  MPI_Aint displacement;
  MPI_Datatype type;
  char const *field_name;
  char const *type_name;
};

ostream &operator<<(ostream &os, mpi_struct_descr_t const &descr);

MPI_Datatype create_mpi_datatype(size_t const count,
                                 mpi_struct_descr_t const descr[],
                                 char const *name, size_t size);
#if 0
  
  class generic_mpi_datatype_t {
    
    string const type_name;
    virtual size_t type_size() const = 0;
    
    struct field_t {
      size_t       offset;
      size_t       count;
      MPI_Datatype mpi_datatype;
      string       field_name;
      string       type_name;
      field_t (size_t const offset_,
               size_t const count_,
               MPI_Datatype const mpi_datatype_,
               string const field_name_,
               string const type_name_)
        : offset(offset_),
          count(count_),
          mpi_datatype(mpi_datatype_),
          field_name(field_name_),
          type_name(type_name_)
      {
      }
      ostream& output (ostream& os) const;
    };
    friend ostream& operator<< (ostream& os,
                                generic_mpi_datatype_t::field_t const& field);
    
    list<field_t> entries;
    
    bool         type_is_committed;
    MPI_Datatype mpi_datatype;
    
  public:
    
    generic_mpi_datatype_t (string const type_name_);
    
    template <typename U>
    void add_field (size_t offset, size_t count, string field_name);
    
    void commit ();
    
    MPI_Datatype get () const
    {
      assert (type_is_committed);
      return mpi_datatype;
    }
    
    ostream& output (ostream& os) const;
  };

  template <typename T>
  class mpi_datatype_t: public generic_mpi_datatype_t {
    virtual size_t type_size() const
    {
      return sizeof(T);
    }
  };

  inline ostream& operator<< (ostream& os,
                              generic_mpi_datatype_t::field_t const& field)
  {
    return field.output(os);
  }
  
  inline ostream& operator<< (ostream& os, generic_mpi_datatype_t const& type)
  {
    return type.output(os);
  }

#endif

// Debugging output
#define CHECKPOINT dist::checkpoint(__FILE__, __LINE__)
void checkpoint(const char *file, int line);

// Information about the communicator

// Return the communicator
inline MPI_Comm comm() { return comm_; }

// Always return a good communicator
inline MPI_Comm goodcomm() {
  return comm_ != MPI_COMM_NULL ? comm_ : MPI_COMM_WORLD;
}

// Rank in the communicator (this process's number, 0 .. size-1)
inline int rank() {
  static int rank_ = -1;
  if (rank_ == -1)
    MPI_Comm_rank(comm(), &rank_);
  return rank_;
}

// Size of the communicator
inline int size() {
  static int size_ = -1;
  if (size_ == -1)
    MPI_Comm_size(comm(), &size_);
  return size_;
}

// Set number of threads
void set_num_threads(int num_threads);

// Local number of threads
inline int num_threads() { return omp_get_max_threads(); }

// My own thread number
inline int thread_num() { return omp_get_thread_num(); }

// Global number of threads
void collect_total_num_threads();
inline int total_num_threads() { return total_num_threads_; }

/////////////////////////////////////////////////////////////////////////
// C Datatype helpers
// Map a C datatype to a 0-based index running up to c_ndatatypes().
/////////////////////////////////////////////////////////////////////////
inline unsigned int c_datatype(const char &) { return 0; }

inline unsigned int c_datatype(const signed char &) { return 1; }

inline unsigned int c_datatype(const unsigned char &) { return 2; }

inline unsigned int c_datatype(const short &) { return 3; }

inline unsigned int c_datatype(const unsigned short &) { return 4; }

inline unsigned int c_datatype(const int &) { return 5; }

inline unsigned int c_datatype(const unsigned int &) { return 6; }

inline unsigned int c_datatype(const long &) { return 7; }

inline unsigned int c_datatype(const unsigned long &) { return 8; }

inline unsigned int c_datatype(const long long &) { return 9; }

inline unsigned int c_datatype(const unsigned long long &) { return 10; }

inline unsigned int c_datatype(const float &) { return 11; }

inline unsigned int c_datatype(const double &) { return 12; }

inline unsigned int c_datatype(const long double &) { return 13; }

#ifdef HAVE_CCTK_COMPLEX8
inline unsigned int c_datatype(const CCTK_COMPLEX8 &) { return 14; }
#endif

#ifdef HAVE_CCTK_COMPLEX16
inline unsigned int c_datatype(const CCTK_COMPLEX16 &) { return 15; }
#endif

#ifdef HAVE_CCTK_COMPLEX32
inline unsigned int c_datatype(const CCTK_COMPLEX32 &) { return 16; }
#endif

// keep this function's return code consistent with functions above
inline unsigned int c_ndatatypes() { return 17; }

template <typename T> unsigned int c_datatype() { abort(); }
template <> inline unsigned int c_datatype<char>() { return 0; }
template <> inline unsigned int c_datatype<signed char>() { return 1; }
template <> inline unsigned int c_datatype<unsigned char>() { return 2; }
template <> inline unsigned int c_datatype<short>() { return 3; }
template <> inline unsigned int c_datatype<unsigned short>() { return 4; }
template <> inline unsigned int c_datatype<int>() { return 5; }
template <> inline unsigned int c_datatype<unsigned int>() { return 6; }
template <> inline unsigned int c_datatype<long>() { return 7; }
template <> inline unsigned int c_datatype<unsigned long>() { return 8; }
template <> inline unsigned int c_datatype<long long>() { return 9; }
template <> inline unsigned int c_datatype<unsigned long long>() { return 10; }
template <> inline unsigned int c_datatype<float>() { return 11; }
template <> inline unsigned int c_datatype<double>() { return 12; }
template <> inline unsigned int c_datatype<long double>() { return 13; }
#ifdef HAVE_CCTK_COMPLEX8
template <> inline unsigned int c_datatype<CCTK_COMPLEX8>() { return 14; }
#endif
#ifdef HAVE_CCTK_COMPLEX16
template <> inline unsigned int c_datatype<CCTK_COMPLEX16>() { return 15; }
#endif
#ifdef HAVE_CCTK_COMPLEX32
template <> inline unsigned int c_datatype<CCTK_COMPLEX32>() { return 16; }
#endif

// Map a C datatype index to a string
char const *c_datatype_name(unsigned type) CCTK_ATTRIBUTE_CONST;

/////////////////////////////////////////////////////////////////
// MPI Datatype helpers
// Map a C datatype to its corresponding MPI datatype.
/////////////////////////////////////////////////////////////////
inline MPI_Datatype mpi_datatype(const char &) { return MPI_CHAR; }

inline MPI_Datatype mpi_datatype(const signed char &) { return MPI_CHAR; }

inline MPI_Datatype mpi_datatype(const unsigned char &) {
  return MPI_UNSIGNED_CHAR;
}

inline MPI_Datatype mpi_datatype(const short &) { return MPI_SHORT; }

inline MPI_Datatype mpi_datatype(const unsigned short &) {
  return MPI_UNSIGNED_SHORT;
}

inline MPI_Datatype mpi_datatype(const int &) { return MPI_INT; }

inline MPI_Datatype mpi_datatype(const unsigned int &) { return MPI_UNSIGNED; }

inline MPI_Datatype mpi_datatype(const long &) { return MPI_LONG; }

inline MPI_Datatype mpi_datatype(const unsigned long &) {
  return MPI_UNSIGNED_LONG;
}

inline MPI_Datatype mpi_datatype(const long long &) {
  return MPI_LONG_LONG_INT;
}

inline MPI_Datatype mpi_datatype(const unsigned long long &) {
  return MPI_LONG_LONG_INT;
} // should be unsigned, but this doesn't exist

inline MPI_Datatype mpi_datatype(const float &) { return MPI_FLOAT; }

inline MPI_Datatype mpi_datatype(const double &) { return MPI_DOUBLE; }

inline MPI_Datatype mpi_datatype(const long double &) {
  return MPI_LONG_DOUBLE;
}

#ifdef HAVE_CCTK_COMPLEX8
inline MPI_Datatype mpi_datatype(const CCTK_COMPLEX8 &) { return mpi_complex8; }
#endif

#ifdef HAVE_CCTK_COMPLEX16
inline MPI_Datatype mpi_datatype(const CCTK_COMPLEX16 &) {
  return mpi_complex16;
}
#endif

#ifdef HAVE_CCTK_COMPLEX32
inline MPI_Datatype mpi_datatype(const CCTK_COMPLEX32 &) {
  return mpi_complex32;
}
#endif

template <typename T> MPI_Datatype mpi_datatype() { abort(); }
template <> inline MPI_Datatype mpi_datatype<char>() { return MPI_CHAR; }
template <> inline MPI_Datatype mpi_datatype<signed char>() { return MPI_CHAR; }
template <> inline MPI_Datatype mpi_datatype<unsigned char>() {
  return MPI_UNSIGNED_CHAR;
}
template <> inline MPI_Datatype mpi_datatype<short>() { return MPI_SHORT; }
template <> inline MPI_Datatype mpi_datatype<unsigned short>() {
  return MPI_UNSIGNED_SHORT;
}
template <> inline MPI_Datatype mpi_datatype<int>() { return MPI_INT; }
template <> inline MPI_Datatype mpi_datatype<unsigned int>() {
  return MPI_UNSIGNED;
}
template <> inline MPI_Datatype mpi_datatype<long>() { return MPI_LONG; }
template <> inline MPI_Datatype mpi_datatype<unsigned long>() {
  return MPI_UNSIGNED_LONG;
}
template <> inline MPI_Datatype mpi_datatype<long long>() {
  return MPI_LONG_LONG_INT;
}
template <> inline MPI_Datatype mpi_datatype<unsigned long long>() {
  return MPI_LONG_LONG_INT;
} // should be unsigned, but this doesn't exist
template <> inline MPI_Datatype mpi_datatype<float>() { return MPI_FLOAT; }
template <> inline MPI_Datatype mpi_datatype<double>() { return MPI_DOUBLE; }
template <> inline MPI_Datatype mpi_datatype<long double>() {
  return MPI_LONG_DOUBLE;
}
#ifdef HAVE_CCTK_COMPLEX8
template <> inline MPI_Datatype mpi_datatype<CCTK_COMPLEX8>() {
  return mpi_complex8;
}
#endif
#ifdef HAVE_CCTK_COMPLEX16
template <> inline MPI_Datatype mpi_datatype<CCTK_COMPLEX16>() {
  return mpi_complex16;
}
#endif
#ifdef HAVE_CCTK_COMPLEX32
template <> inline MPI_Datatype mpi_datatype<CCTK_COMPLEX32>() {
  return mpi_complex32;
}
#endif

} // namespace dist

#endif // DIST_HH
