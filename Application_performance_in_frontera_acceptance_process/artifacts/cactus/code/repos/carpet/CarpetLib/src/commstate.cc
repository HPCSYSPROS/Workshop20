#include <cctk.h>
#include <cctk_Parameters.h>

#include <cassert>
#include <cstdlib>
#include <cstring>
#include <iostream>

#include "bbox.hh"
#include "defs.hh"
#include "dist.hh"
#include "vect.hh"

#include "commstate.hh"
#include "timestat.hh"

using namespace std;
using namespace CarpetLib;

char const *tostring(astate const &thestate) {
  switch (thestate) {
  case state_get_buffer_sizes:
    return "state_get_buffer_sizes";
  case state_fill_send_buffers:
    return "state_fill_send_buffers";
  case state_do_some_work:
    return "state_do_some_work";
  case state_empty_recv_buffers:
    return "state_empty_recv_buffers";
  case state_done:
    return "state_done";
  default:
    assert(0);
  }
  return NULL;
}

comm_state::procbufdesc::procbufdesc()
    : sendbufsize(0), recvbufsize(0), sendbuf(NULL), recvbuf(NULL),
      did_post_send(false), did_post_recv(false) {}

void comm_state::procbufdesc::reinitialize() {
  // Note: calling resize(0) instead of clear() ensures that the
  // vector capacity does not change
  sendbufbase.resize(0);
  recvbufbase.resize(0);
  // Re-initialize the procbuf
  sendbufsize = 0;
  recvbufsize = 0;
  sendbuf = NULL;
  recvbuf = NULL;
  did_post_send = false;
  did_post_recv = false;
}

comm_state::typebufdesc::typebufdesc()
    : in_use(false), mpi_datatype(MPI_DATATYPE_NULL), datatypesize(0) {}

// Define static class members
bool comm_state::typebufs_busy = false;
vector<comm_state::typebufdesc> comm_state::typebufs;
vector<MPI_Request> comm_state::srequests;
vector<MPI_Request> comm_state::rrequests;

// Communication state control
comm_state::comm_state() {
  DECLARE_CCTK_PARAMETERS;

  static Timer timer("commstate::create");
  timer.start();
  thestate = state_get_buffer_sizes;

  assert(not typebufs_busy);
  typebufs_busy = true;
  if (typebufs.empty()) {
    typebufs.resize(dist::c_ndatatypes());
#define TYPECASE(N, T)                                                         \
  {                                                                            \
    T dummy;                                                                   \
    unsigned const type = dist::c_datatype(dummy);                             \
    typebufs.AT(type).mpi_datatype = dist::mpi_datatype(dummy);                \
    typebufs.AT(type).datatypesize = sizeof dummy;                             \
  }
#include "typecase.hh"
#undef TYPECASE
  }

  srequests.reserve(dist::c_ndatatypes() * dist::size());
  rrequests.reserve(dist::c_ndatatypes() * dist::size());
  assert(srequests.empty());
  assert(rrequests.empty());

  timer.stop(0);
}

void comm_state::step() {
  DECLARE_CCTK_PARAMETERS;
  static Timer total("commstate::step");
  total.start();

  if (barrier_between_stages) {
    // Add a barrier, ensuring e.g. that all Irecvs are posted before
    // the first Isends are made
    if (commstate_verbose) {
      CCTK_VInfo(CCTK_THORNSTRING, "before MPI_Barrier; state=%s",
                 tostring(thestate));
    }
    dist::barrier(dist::comm(), 404924393, "CarpetLib::comm_state::step");
    if (commstate_verbose) {
      CCTK_INFO("after MPI_Barrier");
    }
  }

  switch (thestate) {

  case state_get_buffer_sizes: {

    if (check_communication_schedule) {
      vector<int> sendcount(dist::size() * dist::c_ndatatypes());
      for (unsigned type = 0; type < dist::c_ndatatypes(); ++type) {
        for (int proc = 0; proc < dist::size(); ++proc) {
          sendcount.AT(proc * dist::c_ndatatypes() + type) =
              typebufs.AT(type).in_use
                  ? typebufs.AT(type).procbufs.AT(proc).sendbufsize
                  : 0;
        }
        assert(sendcount.AT(dist::rank() * dist::c_ndatatypes() + type) == 0);
      }
      vector<int> recvcount(dist::size() * dist::c_ndatatypes());
      if (commstate_verbose) {
        CCTK_INFO("before MPI_Alltoall");
      }
      MPI_Alltoall(&sendcount.front(), dist::c_ndatatypes(), MPI_INT,
                   &recvcount.front(), dist::c_ndatatypes(), MPI_INT,
                   dist::comm());
      if (commstate_verbose) {
        CCTK_INFO("after MPI_Alltoall");
      }
      for (unsigned type = 0; type < dist::c_ndatatypes(); ++type) {
        for (int proc = 0; proc < dist::size(); ++proc) {
          assert(recvcount.AT(proc * dist::c_ndatatypes() + type) ==
                 (typebufs.AT(type).in_use
                      ? int(typebufs.AT(type).procbufs.AT(proc).recvbufsize)
                      : 0));
        }
        assert(recvcount.AT(dist::rank() * dist::c_ndatatypes() + type) == 0);
      }
    }

    // The sizes of the collective communication buffers are known so
    // now allocate them.
    // The receive operations are also posted here already (a clever
    // MPI layer may take advantage of such early posting).

    for (unsigned type = 0; type < dist::c_ndatatypes(); ++type) {
      if (typebufs.AT(type).in_use) {

        for (int proc1 = 0; proc1 < dist::size(); ++proc1) {
          int const proc = interleave_communications
                               ? (proc1 + dist::rank()) % dist::size()
                               : proc1;

          int const datatypesize = typebufs.AT(type).datatypesize;
          procbufdesc &procbuf = typebufs.AT(type).procbufs.AT(proc);

          assert(procbuf.sendbufbase.empty());
          assert(procbuf.recvbufbase.empty());
          procbuf.sendbufbase.resize(procbuf.sendbufsize * datatypesize *
                                     message_size_multiplier);
          procbuf.recvbufbase.resize(procbuf.recvbufsize * datatypesize *
                                     message_size_multiplier);
          // TODO: this may be a bit extreme, and it is only for
          // internal consistency checking
          if (poison_new_memory) {
            memset(&procbuf.sendbufbase.front(), poison_value,
                   procbuf.sendbufsize * datatypesize *
                       message_size_multiplier);
            memset(&procbuf.recvbufbase.front(), poison_value,
                   procbuf.recvbufsize * datatypesize *
                       message_size_multiplier);
          }
          procbuf.sendbuf = &procbuf.sendbufbase.front();
          procbuf.recvbuf = &procbuf.recvbufbase.front();

          if (procbuf.recvbufsize > 0) {
            static Timer timer("commstate::sizes_irecv");
            timer.start();
            int const tag = type;
            if (commstate_verbose) {
              CCTK_VInfo(CCTK_THORNSTRING,
                         "About to MPI_Irecv from process %d for type %s", proc,
                         dist::c_datatype_name(type));
            }
            MPI_Irecv(&procbuf.recvbufbase.front(),
                      procbuf.recvbufsize * message_size_multiplier,
                      typebufs.AT(type).mpi_datatype, proc, tag, dist::comm(),
                      &push_back(rrequests));
            if (commstate_verbose) {
              CCTK_INFO("Finished MPI_Irecv");
            }
            assert(not procbuf.did_post_recv);
            procbuf.did_post_recv = true;
            timer.stop(procbuf.recvbufsize * datatypesize);
          }

        } // for proc
      }
    } // for type

    if (check_communication_schedule) {
      for (unsigned type = 0; type < dist::c_ndatatypes(); ++type) {
        if (typebufs.AT(type).in_use) {
          for (int proc = 0; proc < dist::size(); ++proc) {
            procbufdesc const &procbuf = typebufs.AT(type).procbufs.AT(proc);
            assert(procbuf.did_post_recv == (procbuf.recvbufsize > 0));
          }
        }
      }
    }

    thestate = state_fill_send_buffers;
    break;
  }

  case state_fill_send_buffers: {
    if (combine_sends) {
      for (unsigned type = 0; type < dist::c_ndatatypes(); ++type) {
        if (typebufs.AT(type).in_use) {

          for (int proc1 = 0; proc1 < dist::size(); ++proc1) {
            int const proc =
                interleave_communications
                    ? (proc1 + dist::size() - dist::rank()) % dist::size()
                    : proc1;

            procbufdesc &procbuf = typebufs.AT(type).procbufs.AT(proc);
            if (procbuf.sendbufsize > 0) {

              int const datatypesize = typebufs.AT(type).datatypesize;

              size_t const fillstate =
                  procbuf.sendbuf - &procbuf.sendbufbase.front();
              assert(fillstate == procbuf.sendbufsize * datatypesize);

              // Enlarge messages for performance testing
              if (message_size_multiplier > 1) {
                size_t const nbytes = procbuf.sendbufsize * datatypesize *
                                      (message_size_multiplier - 1);
                memset(procbuf.sendbuf, poison_value, nbytes);
              }

              int const tag = type;
              if (use_mpi_send) {
                // use MPI_Send
                static Timer timer("commstate::send");
                timer.start();
                if (commstate_verbose) {
                  CCTK_VInfo(CCTK_THORNSTRING,
                             "About to MPI_Send to process %d for type %s",
                             proc, dist::c_datatype_name(type));
                }
                MPI_Send(const_cast<char *>(&procbuf.sendbufbase.front()),
                         procbuf.sendbufsize * message_size_multiplier,
                         typebufs.AT(type).mpi_datatype, proc, tag,
                         dist::comm());
                assert(not procbuf.did_post_send);
                procbuf.did_post_send = true;
                if (commstate_verbose) {
                  CCTK_INFO("Finished MPI_Send");
                }
                timer.stop(procbuf.sendbufsize * datatypesize);
              } else if (use_mpi_ssend) {
                // use MPI_Ssend
                static Timer timer("commstate::ssend");
                timer.start();
                if (commstate_verbose) {
                  CCTK_VInfo(CCTK_THORNSTRING,
                             "About to MPI_Ssend to process %d for type %s",
                             proc, dist::c_datatype_name(type));
                }
                MPI_Ssend(const_cast<char *>(&procbuf.sendbufbase.front()),
                          procbuf.sendbufsize * message_size_multiplier,
                          typebufs.AT(type).mpi_datatype, proc, tag,
                          dist::comm());
                assert(not procbuf.did_post_send);
                procbuf.did_post_send = true;
                if (commstate_verbose) {
                  CCTK_INFO("Finished MPI_Ssend");
                }
                timer.stop(procbuf.sendbufsize * datatypesize);
              } else {
                // use MPI_Isend
                static Timer timer("commstate::isend");
                timer.start();
                if (commstate_verbose) {
                  CCTK_VInfo(CCTK_THORNSTRING,
                             "About to MPI_Isend to process %d for type %s",
                             proc, dist::c_datatype_name(type));
                }
                MPI_Isend(const_cast<char *>(&procbuf.sendbufbase.front()),
                          procbuf.sendbufsize * message_size_multiplier,
                          typebufs.AT(type).mpi_datatype, proc, tag,
                          dist::comm(), &push_back(srequests));
                assert(not procbuf.did_post_send);
                procbuf.did_post_send = true;
                if (commstate_verbose) {
                  CCTK_INFO("Finished MPI_Isend");
                }
                timer.stop(procbuf.sendbufsize * datatypesize);
              }
            }
          } // for proc
        }
      } // for type
    }   // if combine_sends

    if (check_communication_schedule) {
      for (unsigned type = 0; type < dist::c_ndatatypes(); ++type) {
        if (typebufs.AT(type).in_use) {
          for (int proc = 0; proc < dist::size(); ++proc) {
            procbufdesc const &procbuf = typebufs.AT(type).procbufs.AT(proc);
            assert(procbuf.did_post_send == (procbuf.sendbufsize > 0));
          }
        }
      }
    }

    thestate = state_do_some_work;
    break;
  }

  case state_do_some_work: {
    static Timer timer("commstate::do_some_work::waitall");
    timer.start();
    if (commstate_verbose) {
      CCTK_INFO("About to MPI_Waitall");
    }
    MPI_Waitall(rrequests.size(), &rrequests.front(), MPI_STATUSES_IGNORE);
    if (commstate_verbose) {
      CCTK_INFO("Finished MPI_Waitall");
    }
    timer.stop(0);

    thestate = state_empty_recv_buffers;
    break;
  }

  case state_empty_recv_buffers: {
    static Timer timer("commstate::empty_recv_buffers::waitall");
    timer.start();
    if (commstate_verbose) {
      CCTK_INFO("About to MPI_Waitall");
    }
    MPI_Waitall(srequests.size(), &srequests.front(), MPI_STATUSES_IGNORE);
    if (commstate_verbose) {
      CCTK_INFO("Finished MPI_Waitall");
    }
    timer.stop(0);

    // Transfer messages again for performance testing
    for (int n = 1; n < message_count_multiplier; ++n) {

      srequests.resize(0);
      srequests.reserve(dist::c_ndatatypes() * dist::size());
      rrequests.resize(0);
      rrequests.reserve(dist::c_ndatatypes() * dist::size());

      // Irecv
      for (unsigned type = 0; type < dist::c_ndatatypes(); ++type) {
        if (typebufs.AT(type).in_use) {

          for (int proc1 = 0; proc1 < dist::size(); ++proc1) {
            int const proc = interleave_communications
                                 ? (proc1 + dist::rank()) % dist::size()
                                 : proc1;

            procbufdesc &procbuf = typebufs.AT(type).procbufs.AT(proc);

            if (procbuf.recvbufsize > 0) {
              static Timer timer("commstate::message_count_multiplier::irecv");
              timer.start();
              int const tag = type;
              if (commstate_verbose) {
                CCTK_VInfo(CCTK_THORNSTRING,
                           "About to MPI_Irecv from process %d for type %s",
                           proc, dist::c_datatype_name(type));
              }
              MPI_Irecv(&procbuf.recvbufbase.front(),
                        procbuf.recvbufsize * message_size_multiplier,
                        typebufs.AT(type).mpi_datatype, proc, tag, dist::comm(),
                        &push_back(rrequests));
              if (commstate_verbose) {
                CCTK_INFO("Finished MPI_Irecv");
              }
              timer.stop(procbuf.recvbufsize * typebufs.AT(type).datatypesize);
            }

          } // for proc
        }
      } // for type

      // Isend
      for (unsigned type = 0; type < dist::c_ndatatypes(); ++type) {
        if (typebufs.AT(type).in_use) {

          for (int proc1 = 0; proc1 < dist::size(); ++proc1) {
            int const proc =
                interleave_communications
                    ? (proc1 + dist::size() - dist::rank()) % dist::size()
                    : proc1;

            procbufdesc &procbuf = typebufs.AT(type).procbufs.AT(proc);

            if (procbuf.sendbufsize > 0) {
              int const tag = type;
              assert(not use_mpi_send);
              assert(not use_mpi_ssend);
              static Timer timer("commstate::message_count_multiplier::isend");
              timer.start();
              if (commstate_verbose) {
                CCTK_VInfo(CCTK_THORNSTRING,
                           "About to MPI_Isend to process %d for type %s", proc,
                           dist::c_datatype_name(type));
              }
              MPI_Isend(const_cast<char *>(&procbuf.sendbufbase.front()),
                        procbuf.sendbufsize * message_size_multiplier,
                        typebufs.AT(type).mpi_datatype, proc, tag, dist::comm(),
                        &push_back(srequests));
              if (commstate_verbose) {
                CCTK_INFO("Finished MPI_Isend");
              }
              timer.stop(procbuf.sendbufsize * typebufs.AT(type).datatypesize);
            }

          } // for proc
        }
      } // for type

      // Waitall
      {
        static Timer timer(
            "commstate::message_count_multiplier::waitall(irecv)");
        timer.start();
        if (commstate_verbose) {
          CCTK_INFO("About to MPI_Waitall");
        }
        MPI_Waitall(rrequests.size(), &rrequests.front(), MPI_STATUSES_IGNORE);
        if (commstate_verbose) {
          CCTK_INFO("Finished MPI_Waitall");
        }
        timer.stop(0);
      }

      // Waitall
      {
        static Timer timer(
            "commstate::message_count_multiplier::waitall(isend)");
        timer.start();
        if (commstate_verbose) {
          CCTK_INFO("About to MPI_Waitall");
        }
        MPI_Waitall(srequests.size(), &srequests.front(), MPI_STATUSES_IGNORE);
        if (commstate_verbose) {
          CCTK_INFO("Finished MPI_Waitall");
        }
        timer.stop(0);
      }

    } // for n

    thestate = state_done;
    break;
  }

  case state_done: {
    assert(0);
  }

  default:
    assert(0);
  }

  total.stop(0);
}

bool comm_state::done() const { return thestate == state_done; }

comm_state::~comm_state() {
  DECLARE_CCTK_PARAMETERS;

  // Note: calling resize(0) instead of clear() ensures that the
  // vector capacity does not change
  srequests.resize(0);
  rrequests.resize(0);

  for (size_t type = 0; type < typebufs.size(); ++type) {
    typebufdesc &typebuf = typebufs.AT(type);
    for (size_t proc = 0; proc < typebuf.procbufs.size(); ++proc) {
      procbufdesc &procbuf = typebuf.procbufs.AT(proc);
      procbuf.reinitialize();
    }
  }

  assert(typebufs_busy);
  typebufs_busy = false;

  assert(thestate == state_done or thestate == state_get_buffer_sizes);
}

void comm_state::reserve_send_space(unsigned const type, int const proc,
                                    int const npoints) {
  assert(type < dist::c_ndatatypes());
  assert(proc >= 0 and proc < dist::size());
  assert(npoints >= 0);
  typebufdesc &typebuf = typebufs.AT(type);
  if (not typebuf.in_use) {
    typebuf.procbufs.resize(dist::size());
    typebuf.in_use = true;
  }
  procbufdesc &procbuf = typebuf.procbufs.AT(proc);
  procbuf.sendbufsize += npoints;
}

void comm_state::reserve_recv_space(unsigned const type, int const proc,
                                    int const npoints) {
  assert(type < dist::c_ndatatypes());
  assert(proc >= 0 and proc < dist::size());
  assert(npoints >= 0);
  typebufdesc &typebuf = typebufs.AT(type);
  if (not typebuf.in_use) {
    typebuf.procbufs.resize(dist::size());
    typebuf.in_use = true;
  }
  procbufdesc &procbuf = typebuf.procbufs.AT(proc);
  procbuf.recvbufsize += npoints;
}

void *comm_state::send_buffer(unsigned const type, int const proc,
                              int const npoints) {
  assert(type < dist::c_ndatatypes());
  assert(proc >= 0 and proc < dist::size());
  assert(npoints > 0);
  typebufdesc const &typebuf = typebufs.AT(type);
  procbufdesc const &procbuf = typebuf.procbufs.AT(proc);

  assert(procbuf.sendbuf + npoints * typebuf.datatypesize <=
         &procbuf.sendbufbase.front() +
             procbuf.sendbufsize * typebuf.datatypesize);

  return procbuf.sendbuf;
}

void *comm_state::recv_buffer(unsigned const type, int const proc,
                              int const npoints) {
  assert(type < dist::c_ndatatypes());
  assert(proc >= 0 and proc < dist::size());
  assert(npoints > 0);
  typebufdesc const &typebuf = typebufs.AT(type);
  procbufdesc const &procbuf = typebuf.procbufs.AT(proc);

  assert(procbuf.recvbuf + npoints * typebuf.datatypesize <=
         &procbuf.recvbufbase.front() +
             procbuf.recvbufsize * typebuf.datatypesize);

  return procbuf.recvbuf;
}

void comm_state::commit_send_space(unsigned const type, int const proc,
                                   int const npoints) {
  DECLARE_CCTK_PARAMETERS;

  assert(type < dist::c_ndatatypes());
  assert(proc >= 0 and proc < dist::size());
  assert(npoints >= 0);
  assert(npoints > 0);
  typebufdesc &typebuf = typebufs.AT(type);
  procbufdesc &procbuf = typebuf.procbufs.AT(proc);
  procbuf.sendbuf += npoints * typebuf.datatypesize;
  assert(procbuf.sendbuf <= &procbuf.sendbufbase.front() +
                                procbuf.sendbufsize * typebuf.datatypesize);

  if (not combine_sends) {
    // post the send if the buffer is full
    if (procbuf.sendbuf ==
        &procbuf.sendbufbase.front() +
            procbuf.sendbufsize * typebuf.datatypesize) {
      if (message_size_multiplier > 1) {
        size_t const nbytes = procbuf.sendbufsize * typebuf.datatypesize *
                              (message_size_multiplier - 1);
        memset(procbuf.sendbuf, poison_value, nbytes);
      }

      static Timer timer("commit_send_space::isend");
      timer.start();
      if (commstate_verbose) {
        CCTK_VInfo(CCTK_THORNSTRING,
                   "About to MPI_Isend to process %d for type %s", proc,
                   dist::c_datatype_name(type));
      }
      int const tag = type;
      assert(procbuf.sendbufsize > 0);
      assert(not use_mpi_send);
      assert(not use_mpi_ssend);
      MPI_Isend(&procbuf.sendbufbase.front(),
                procbuf.sendbufsize * message_size_multiplier,
                typebuf.mpi_datatype, proc, tag, dist::comm(),
                &push_back(srequests));
      assert(not procbuf.did_post_send);
      procbuf.did_post_send = true;
      if (commstate_verbose) {
        CCTK_INFO("Finished MPI_Isend");
      }
      timer.stop(procbuf.sendbufsize * typebuf.datatypesize);
    }
  }
}

void comm_state::commit_recv_space(unsigned const type, int const proc,
                                   int const npoints) {
  assert(type < dist::c_ndatatypes());
  assert(proc >= 0 and proc < dist::size());
  assert(npoints >= 0);
  assert(npoints > 0);
  typebufdesc &typebuf = typebufs.AT(type);
  procbufdesc &procbuf = typebuf.procbufs.AT(proc);
  procbuf.recvbuf += npoints * typebuf.datatypesize;
  assert(procbuf.recvbuf <= &procbuf.recvbufbase.front() +
                                procbuf.recvbufsize * typebuf.datatypesize);
}
