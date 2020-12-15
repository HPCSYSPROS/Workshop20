#include <cctk.h>

#include <algorithm>
#include <cassert>
#include <cstring>
#include <string>
#include <vector>

#ifdef CCTK_MPI
#include <mpi.h>
#else
#include "nompi.h"
#endif

#include "dh.hh"
#include "mpi_string.hh"
#include "region.hh"

namespace CarpetLib {

using namespace std;

vector<string> gather_string(MPI_Comm const comm, int const root,
                             string const &data) {
  // Get my rank
  int rank;
  MPI_Comm_rank(comm, &rank);

  if (rank == root) {

    // Get the total number of processes
    int num_procs;
    MPI_Comm_size(comm, &num_procs);

    // Gather the lengths of the data strings
    int const length = data.length();
    vector<int> lengths(num_procs);

    MPI_Gather(const_cast<int *>(&length), 1, MPI_INT, &lengths.front(), 1,
               MPI_INT, root, comm);

    // Allocate space for all data strings
    vector<int> offsets(num_procs + 1);
    offsets.AT(0) = 0;
    for (int n = 0; n < num_procs; ++n) {
      offsets.AT(n + 1) = offsets.AT(n) + lengths.AT(n);
    }
    int const total_length = offsets.AT(num_procs);
    vector<char> alldata_buffer(total_length);

    // Gather all data strings
    MPI_Gatherv(const_cast<char *>(data.c_str()), length, MPI_CHAR,
                &alldata_buffer.front(), const_cast<int *>(&lengths.front()),
                const_cast<int *>(&offsets.front()), MPI_CHAR, root, comm);

    // Convert data buffer with C strings to C++ strings
    vector<string> alldata(num_procs);
    for (int n = 0; n < num_procs; ++n) {
      alldata.AT(n) = string(&alldata_buffer.AT(offsets.AT(n)), lengths.AT(n));
    }

    return alldata;

  } else {

    // Gather the lengths of the data strings
    int const length = data.length();

    MPI_Gather(const_cast<int *>(&length), 1, MPI_INT, NULL, 1, MPI_INT, root,
               comm);

    // Gather all data strings
    MPI_Gatherv(const_cast<char *>(data.c_str()), length, MPI_CHAR, NULL, NULL,
                NULL, MPI_CHAR, root, comm);

    // Convert data buffer with C strings to C++ strings
    vector<string> alldata;

    return alldata;
  }
}

vector<string> allgather_string(MPI_Comm const comm, string const &data) {
  // Get the total number of processes
  int num_procs;
  MPI_Comm_size(comm, &num_procs);

  // Gather the lengths of the data strings
  int const length = data.length();
  vector<int> lengths(num_procs);

  MPI_Allgather(const_cast<int *>(&length), 1, MPI_INT, &lengths.front(), 1,
                MPI_INT, comm);

  // Allocate space for all data strings
  vector<int> offsets(num_procs + 1);
  offsets.AT(0) = 0;
  for (int n = 0; n < num_procs; ++n) {
    offsets.AT(n + 1) = offsets.AT(n) + lengths.AT(n);
  }
  int const total_length = offsets.AT(num_procs);
  vector<char> alldata_buffer(total_length);

  // Gather all data strings
  MPI_Allgatherv(const_cast<char *>(data.c_str()), length, MPI_CHAR,
                 &alldata_buffer.front(), const_cast<int *>(&lengths.front()),
                 const_cast<int *>(&offsets.front()), MPI_CHAR, comm);

  // Convert data buffer with C strings to C++ strings
  vector<string> alldata(num_procs);
  for (int n = 0; n < num_procs; ++n) {
    alldata.AT(n) = string(&alldata_buffer.AT(offsets.AT(n)), lengths.AT(n));
  }

  return alldata;
}

vector<string> alltoallv_string(MPI_Comm const comm,
                                vector<string> const &data) {
  // Get the total number of processes
  int num_procs;
  MPI_Comm_size(comm, &num_procs);

  // Exchange the lengths of the data strings
  vector<int> lengths_in(num_procs);
  for (int n = 0; n < num_procs; ++n) {
    lengths_in.AT(n) = data.AT(n).length();
  }
  vector<int> lengths(num_procs);
  MPI_Alltoall(&lengths_in.front(), 1, MPI_INT, &lengths.front(), 1, MPI_INT,
               comm);

  // Allocate space for all data strings
  vector<int> offsets_in(num_procs + 1);
  offsets_in.AT(0) = 0;
  for (int n = 0; n < num_procs; ++n) {
    offsets_in.AT(n + 1) = offsets_in.AT(n) + lengths_in.AT(n);
  }
  int const total_length_in = offsets_in.AT(num_procs);
  vector<char> alldata_buffer_in(total_length_in);

  vector<int> offsets(num_procs + 1);
  offsets.AT(0) = 0;
  for (int n = 0; n < num_procs; ++n) {
    offsets.AT(n + 1) = offsets.AT(n) + lengths.AT(n);
  }
  int const total_length = offsets.AT(num_procs);
  vector<char> alldata_buffer(total_length);

  // Convert C++ strings to data buffer with C strings
  for (int n = 0; n < num_procs; ++n) {
    memcpy(&alldata_buffer_in.AT(offsets_in.AT(n)), data.AT(n).c_str(),
           lengths_in.AT(n));
  }

  // Exchange all data strings
  MPI_Alltoallv(&alldata_buffer_in.front(), &lengths_in.front(),
                &offsets_in.front(), MPI_CHAR, &alldata_buffer.front(),
                &lengths.front(), &offsets.front(), MPI_CHAR, comm);

  // Convert data buffer with C strings to C++ strings
  vector<string> alldata(num_procs);
  for (int n = 0; n < num_procs; ++n) {
    alldata.AT(n) = string(&alldata_buffer.AT(offsets.AT(n)), lengths.AT(n));
  }

  return alldata;
}

string broadcast_string(MPI_Comm const comm, int const root,
                        string const &data) {
  // Get my rank
  int rank;
  MPI_Comm_rank(comm, &rank);

  if (rank == root) {

    // Broadcast the length of the data string
    int const length = data.length();
    MPI_Bcast(const_cast<int *>(&length), 1, MPI_INT, root, comm);

    // Broadcast data string
    char const *const buf = data.c_str();
    MPI_Bcast(const_cast<char *>(buf), length, MPI_CHAR, root, comm);

    // Return original string
    return data;

  } else {

    // Broadcast the length of the data string
    int length;
    MPI_Bcast(&length, 1, MPI_INT, root, comm);

    // Allocate space for data string
    vector<char> data_buffer(length);

    // Broadcast data string
    char *const buf = &data_buffer.front();
    MPI_Bcast(buf, length, MPI_CHAR, root, comm);

    // Convert data buffer with C strings to C++ strings
    string const result = string(&data_buffer.front(), length);

    return result;
  }
}

//////////////////////////////////////////////////////////////////////////////

template vector<vector<dh::light_dboxes> >
allgatherv(MPI_Comm comm, vector<dh::light_dboxes> const &data);

template vector<ivect> allgatherv1(MPI_Comm comm, vector<ivect> const &data);

} // namespace CarpetLib
