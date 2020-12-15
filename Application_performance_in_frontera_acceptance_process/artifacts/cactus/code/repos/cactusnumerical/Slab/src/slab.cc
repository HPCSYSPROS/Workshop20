// TODO:
// Provide facilities for dim > 3
// Set up the slab exchange information in advance
// Slab in several stages
// Test slabbing without MPI
// Allow not setting the ghost zones
// Allow not using / not setting the boundaries
// Allow different numbers of ghost zones at the lower and upper boundary



// Print debug information?
#undef DEBUG

// Perform expensive self-checks?
#undef CHECK

// Omit all self-checks?  (Overrides CHECK)
#undef NDEBUG

// Byte value for poison checks: use 255 for nan, or e.g. 113 for a
// large value
#define POISON_VALUE 254



#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <cctk_DefineThorn.h>
#include <util_ErrorCodes.h>
#include <util_Table.h>

#include <loopcontrol.h>

#ifdef CCTK_MPI
#  include <mpi.h>
#  define HAVE_MPI 1
#else
#  define HAVE_MPI 0
#endif

#include "slab.h"

using namespace std;



#ifdef DEBUG
#  define ifdebug
#else
#  define ifdebug while (0)
#endif

#ifdef CHECK
#  define ifcheck
#else
#  define ifcheck while (0)
#endif

#ifdef NDEBUG
#  define check(x) ((x) ? 0 : CCTK_WARN (0, "internal error"))
#else
#  define check(x) assert (x)
#endif



static int timer_init                  = -1;
static int timer_apply                 = -1;
static int timer_copy_in               = -1;
static int timer_copy_in_noxpose       = -1;
static int timer_copy_in_xposexy       = -1;
static int timer_copy_in_xposegeneral  = -1;
static int timer_copy_in_general       = -1;
static int timer_xfer                  = -1;
static int timer_copy_back             = -1;
static int timer_copy_back_noflip      = -1;
static int timer_copy_back_flipx       = -1;
static int timer_copy_back_flipy       = -1;
static int timer_copy_back_flipxy      = -1;
static int timer_copy_back_flipgeneral = -1;
static int timer_copy_back_general     = -1;



extern "C"
void
Slab_InitTimers (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  timer_init                  = CCTK_TimerCreate ("Slab/init"                 );
  timer_apply                 = CCTK_TimerCreate ("Slab/apply"                );
  timer_copy_in               = CCTK_TimerCreate ("Slab/copy_in"              );
  timer_copy_in_noxpose       = CCTK_TimerCreate ("Slab/copy_in_noxpose"      );
  timer_copy_in_xposexy       = CCTK_TimerCreate ("Slab/copy_in_xposexy"      );
  timer_copy_in_xposegeneral  = CCTK_TimerCreate ("Slab/copy_in_xposegeneral" );
  timer_copy_in_general       = CCTK_TimerCreate ("Slab/copy_in_general"      );
  timer_xfer                  = CCTK_TimerCreate ("Slab/xfer"                 );
  timer_copy_back             = CCTK_TimerCreate ("Slab/copy_back"            );
  timer_copy_back_noflip      = CCTK_TimerCreate ("Slab/copy_back_noflip"     );
  timer_copy_back_flipx       = CCTK_TimerCreate ("Slab/copy_back_flipx"      );
  timer_copy_back_flipy       = CCTK_TimerCreate ("Slab/copy_back_flipy"      );
  timer_copy_back_flipxy      = CCTK_TimerCreate ("Slab/copy_back_flipxy"     );
  timer_copy_back_flipgeneral = CCTK_TimerCreate ("Slab/copy_back_flipgeneral");
  timer_copy_back_general     = CCTK_TimerCreate ("Slab/copy_back_general"    );
}

extern "C"
int
Slab_PrintTimers ()
{
  CCTK_TimerPrintDataI (timer_init                 , -1);
  CCTK_TimerPrintDataI (timer_apply                , -1);
  CCTK_TimerPrintDataI (timer_copy_in              , -1);
  CCTK_TimerPrintDataI (timer_copy_in_noxpose      , -1);
  CCTK_TimerPrintDataI (timer_copy_in_xposexy      , -1);
  CCTK_TimerPrintDataI (timer_copy_in_xposegeneral , -1);
  CCTK_TimerPrintDataI (timer_copy_in_general      , -1);
  CCTK_TimerPrintDataI (timer_xfer                 , -1);
  CCTK_TimerPrintDataI (timer_copy_back            , -1);
  CCTK_TimerPrintDataI (timer_copy_back_noflip     , -1);
  CCTK_TimerPrintDataI (timer_copy_back_flipx      , -1);
  CCTK_TimerPrintDataI (timer_copy_back_flipy      , -1);
  CCTK_TimerPrintDataI (timer_copy_back_flipxy     , -1);
  CCTK_TimerPrintDataI (timer_copy_back_flipgeneral, -1);
  CCTK_TimerPrintDataI (timer_copy_back_general    , -1);
  return 0;
}



// Find out which driver to use
#ifdef CCTK_MPI
#  if defined CARPET_CARPET
#    include "Carpet/Carpet/src/carpet_public.h"
#  endif
#  if defined CACTUSPUGH_PUGH
#    include "CactusPUGH/PUGH/src/include/pugh.h"
#  endif
#endif



#ifdef CCTK_MPI
// Determine MPI type sizes

#  define CACTUS_MPI_BYTE MPI_CHAR

#  define CACTUS_MPI_INT1 MPI_CHAR

#  if SIZEOF_SHORT_INT == 2
#    define CACTUS_MPI_INT2 MPI_SHORT
#  elif SIZEOF_INT == 2
#    define CACTUS_MPI_INT2 MPI_INT
#  elif SIZEOF_LONG_INT == 2
#    define CACTUS_MPI_INT2 MPI_LONG
#  elif SIZEOF_LONG_LONG == 2
#    define CACTUS_MPI_INT2 MPI_LONG_LONG_INT
#  endif

#  if SIZEOF_SHORT_INT == 4
#    define CACTUS_MPI_INT4 MPI_SHORT
#  elif SIZEOF_INT == 4
#    define CACTUS_MPI_INT4 MPI_INT
#  elif SIZEOF_LONG_INT == 4
#    define CACTUS_MPI_INT4 MPI_LONG
#  elif SIZEOF_LONG_LONG == 4
#    define CACTUS_MPI_INT4 MPI_LONG_LONG_INT
#  endif

#  if SIZEOF_SHORT_INT == 8
#    define CACTUS_MPI_INT8 MPI_SHORT
#  elif SIZEOF_INT == 8
#    define CACTUS_MPI_INT8 MPI_INT
#  elif SIZEOF_LONG_INT == 8
#    define CACTUS_MPI_INT8 MPI_LONG
#  elif SIZEOF_LONG_LONG == 8
#    define CACTUS_MPI_INT8 MPI_LONG_LONG_INT
#  endif

#  if SIZEOF_SHORT_INT == 16
#    define CACTUS_MPI_INT16 MPI_SHORT
#  elif SIZEOF_INT == 16
#    define CACTUS_MPI_INT16 MPI_INT
#  elif SIZEOF_LONG_INT == 16
#    define CACTUS_MPI_INT16 MPI_LONG
#  elif SIZEOF_LONG_LONG == 16
#    define CACTUS_MPI_INT16 MPI_LONG_LONG_INT
#  endif

#  if SIZEOF_FLOAT == 4
#    define CACTUS_MPI_REAL4 MPI_FLOAT
#  elif SIZEOF_DOUBLE == 4
#    define CACTUS_MPI_REAL4 MPI_DOUBLE
#  elif SIZEOF_LONG_DOUBLE == 4
#    define CACTUS_MPI_REAL4 MPI_LONG_DOUBLE
#  endif

#  if SIZEOF_FLOAT == 8
#    define CACTUS_MPI_REAL8 MPI_FLOAT
#  elif SIZEOF_DOUBLE == 8
#    define CACTUS_MPI_REAL8 MPI_DOUBLE
#  elif SIZEOF_LONG_DOUBLE == 8
#    define CACTUS_MPI_REAL8 MPI_LONG_DOUBLE
#  endif

#  if SIZEOF_FLOAT == 16
#    define CACTUS_MPI_REAL16 MPI_FLOAT
#  elif SIZEOF_DOUBLE == 16
#    define CACTUS_MPI_REAL16 MPI_DOUBLE
#  elif SIZEOF_LONG_DOUBLE == 16
#    define CACTUS_MPI_REAL16 MPI_LONG_DOUBLE
#  endif

static MPI_Datatype CACTUS_MPI_COMPLEX8;
static MPI_Datatype CACTUS_MPI_COMPLEX16;
static MPI_Datatype CACTUS_MPI_COMPLEX32;

#endif



// Replace MPI functions if MPI is disabled
#ifndef CCTK_MPI

typedef int MPI_Comm;

typedef enum {
  CACTUS_MPI_BYTE      = CCTK_VARIABLE_BYTE,
  CACTUS_MPI_INT       = CCTK_VARIABLE_INT,
  CACTUS_MPI_INT1      = CCTK_VARIABLE_INT1,
  CACTUS_MPI_INT2      = CCTK_VARIABLE_INT2,
  CACTUS_MPI_INT4      = CCTK_VARIABLE_INT4,
  CACTUS_MPI_INT8      = CCTK_VARIABLE_INT8,
  CACTUS_MPI_INT16     = CCTK_VARIABLE_INT16,
  CACTUS_MPI_REAL      = CCTK_VARIABLE_REAL,
  CACTUS_MPI_REAL4     = CCTK_VARIABLE_REAL4,
  CACTUS_MPI_REAL8     = CCTK_VARIABLE_REAL8,
  CACTUS_MPI_REAL16    = CCTK_VARIABLE_REAL16,
  CACTUS_MPI_COMPLEX   = CCTK_VARIABLE_COMPLEX,
  CACTUS_MPI_COMPLEX8  = CCTK_VARIABLE_COMPLEX8,
  CACTUS_MPI_COMPLEX16 = CCTK_VARIABLE_COMPLEX16,
  CACTUS_MPI_COMPLEX32 = CCTK_VARIABLE_COMPLEX32
} MPI_Datatype;

static MPI_Datatype MPI_INT;

typedef enum { MPI_MIN, MPI_MAX } MPI_Op;

static int
MPI_Barrier (MPI_Comm comm)
{
  return 0;
}

static int
MPI_Bcast (void * buffer, int count,
           MPI_Datatype datatype, int root, MPI_Comm comm)
{
  assert (buffer);
  assert (count >= 0);
  assert (root == 0);
  return 0;
}

static int
MPI_Comm_size (MPI_Comm comm, int * size)
{
  *size = 1;
  return 0;
}

static int
MPI_Comm_rank (MPI_Comm comm, int * rank)
{
  *rank = 0;
  return 0;
}

static int
MPI_Allgather (void * sendbuf, int sendcnt, int sendtype,
	       void * recvbuf, int recvcnt, int recvtype,
	       MPI_Comm comm)
{
  int recvsize;
  assert (sendbuf);
  assert (recvbuf);
  assert (sendcnt == recvcnt);
  assert (recvcnt >= 0);
  assert (sendtype == recvtype);
  recvsize = CCTK_VarTypeSize (recvtype);
  assert (recvsize > 0);
  memcpy (recvbuf, sendbuf, recvcnt * recvsize);
  return 0;
}

static int
MPI_Alltoall (void * sendbuf, int sendcnt, int sendtype,
	      void * recvbuf, int recvcnt, int recvtype,
	      MPI_Comm comm)
{
  int recvsize;
  assert (sendbuf);
  assert (recvbuf);
  assert (sendcnt == recvcnt);
  assert (recvcnt >= 0);
  assert (sendtype == recvtype);
  recvsize = CCTK_VarTypeSize (recvtype);
  assert (recvsize > 0);
  memcpy (recvbuf, sendbuf, recvcnt * recvsize);
  return 0;
}

static int
MPI_Alltoallv (void * sendbuf, int * sendcnt, int * sendoff, int sendtype,
	       void * recvbuf, int * recvcnt, int * recvoff, int recvtype,
	       MPI_Comm comm)
{
  int recvsize;
  assert (sendbuf);
  assert (recvbuf);
  assert (sendcnt);
  assert (recvcnt);
  assert (*sendcnt == *recvcnt);
  assert (*recvcnt >= 0);
  assert (sendoff);
  assert (recvoff);
  assert (*sendoff == 0);
  assert (*recvoff == 0);
  assert (sendtype == recvtype);
  recvsize = CCTK_VarTypeSize (recvtype);
  assert (recvsize > 0);
  memcpy (recvbuf, sendbuf, *recvcnt * recvsize);
  return 0;
}

static int
MPI_Allreduce (void * sendbuf, void * recvbuf, int count,
               MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{
  int recvsize;
  assert (sendbuf);
  assert (recvbuf);
  assert (count >= 0);
  recvsize = CCTK_VarTypeSize (datatype);
  assert (recvsize > 0);
  memcpy (recvbuf, sendbuf, count * recvsize);
  return 0;
}

typedef int MPI_Request;        // dummy
typedef int MPI_Status;         // dummy
#define MPI_REQUEST_NULL 0      // dummy
#define MPI_STATUSES_IGNORE 0   // dummy

static int
MPI_Irecv (void *buf, int count, MPI_Datatype datatype,
           int source, int tag, MPI_Comm comm, MPI_Request *request)
{
  abort();
}

static int
MPI_Isend (void *buf, int count, MPI_Datatype datatype, int dest,
           int tag, MPI_Comm comm, MPI_Request *request)
{
  abort();
}

static int
MPI_Waitall (int count, MPI_Request *array_of_requests,
             MPI_Status *array_of_statuses)
{
  abort();
}

#endif



// Get the MPI COMM_WOLRD communicator from the driver
static MPI_Comm
get_mpi_comm (cGH const * restrict const cctkGH)
{
#ifdef CCTK_MPI
  if (CCTK_IsFunctionAliased ("GetMPICommWorld")) {
    return * (MPI_Comm const *) GetMPICommWorld (cctkGH);
  }
#  if defined CACTUSPUGH_PUGH
  {
    static int PUGH_active = -1;
    if (PUGH_active == -1) PUGH_active = CCTK_IsThornActive ("PUGH");
    assert (PUGH_active >= 0);
    if (PUGH_active) return PUGH_pGH(cctkGH)->PUGH_COMM_WORLD;
  }
#  endif
  return MPI_COMM_WORLD;
#else
  return 0;
#endif
}



// Initialise the MPI datatypes for complex variables
extern "C"
int
Slab_InitMPIDatatypes ()
{
#ifdef CCTK_MPI
#  ifdef HAVE_CCTK_REAL4
  assert (CACTUS_MPI_REAL4 != MPI_DATATYPE_NULL);
  MPI_Type_contiguous (2, CACTUS_MPI_REAL4, &CACTUS_MPI_COMPLEX8);
  MPI_Type_commit (&CACTUS_MPI_COMPLEX8);
#  endif
#  ifdef HAVE_CCTK_REAL8
  assert (CACTUS_MPI_REAL8 != MPI_DATATYPE_NULL);
  MPI_Type_contiguous (2, CACTUS_MPI_REAL8, &CACTUS_MPI_COMPLEX16);
  MPI_Type_commit (&CACTUS_MPI_COMPLEX16);
#  endif
#  ifdef HAVE_CCTK_REAL16
  if (CACTUS_MPI_REAL16 != MPI_DATATYPE_NULL) {
    MPI_Type_contiguous (2, CACTUS_MPI_REAL16, &CACTUS_MPI_COMPLEX32);
    MPI_Type_commit (&CACTUS_MPI_COMPLEX32);
  } else {
    // CCTK_REAL16 is not supported by MPI
    CCTK_WARN (CCTK_WARN_ALERT,
               "CCTK_REAL16 support is enabled in Cactus, but is not supported by MPI. All MPI operations with this datatype will fail.");
    CACTUS_MPI_COMPLEX32 = MPI_DATATYPE_NULL;
  }
#  endif
#endif
  
#ifndef CCTK_MPI
  switch (sizeof(int)) {
#ifdef HAVE_CCTK_INT1
  case sizeof(CCTK_INT1): MPI_INT = CACTUS_MPI_INT1; break;
#endif
#ifdef HAVE_CCTK_INT2
  case sizeof(CCTK_INT2): MPI_INT = CACTUS_MPI_INT2; break;
#endif
#ifdef HAVE_CCTK_INT4
  case sizeof(CCTK_INT4): MPI_INT = CACTUS_MPI_INT4; break;
#endif
#ifdef HAVE_CCTK_INT8
  case sizeof(CCTK_INT8): MPI_INT = CACTUS_MPI_INT8; break;
#endif
  default: CCTK_BUILTIN_UNREACHABLE();
  }
#endif
  
  return 0;
}



// Normalise a Cactus datatype
static int
normal_type (int cactustype)
{
  switch (cactustype) {
  case CCTK_VARIABLE_INT:
#ifdef CCTK_INTEGER_PRECISION_1
    return CCTK_VARIABLE_INT1;
#endif
#ifdef CCTK_INTEGER_PRECISION_2
    return CCTK_VARIABLE_INT2;
#endif
#ifdef CCTK_INTEGER_PRECISION_4
    return CCTK_VARIABLE_INT4;
#endif
#ifdef CCTK_INTEGER_PRECISION_8
    return CCTK_VARIABLE_INT8;
#endif
#ifdef CCTK_INTEGER_PRECISION_16
    return CCTK_VARIABLE_INT16;
#endif
    CCTK_BUILTIN_UNREACHABLE();
  case CCTK_VARIABLE_REAL:
#ifdef CCTK_REAL_PRECISION_4
    return CCTK_VARIABLE_REAL4;
#endif
#ifdef CCTK_REAL_PRECISION_8
    return CCTK_VARIABLE_REAL8;
#endif
#ifdef CCTK_REAL_PRECISION_16
    return CCTK_VARIABLE_REAL16;
#endif
    CCTK_BUILTIN_UNREACHABLE();
  case CCTK_VARIABLE_COMPLEX:
#ifdef CCTK_REAL_PRECISION_4
    return CCTK_VARIABLE_COMPLEX8;
#endif
#ifdef CCTK_REAL_PRECISION_8
    return CCTK_VARIABLE_COMPLEX16;
#endif
#ifdef CCTK_REAL_PRECISION_16
    return CCTK_VARIABLE_COMPLEX32;
#endif
    CCTK_BUILTIN_UNREACHABLE();
  }
  return cactustype;
}



// Find the MPI datatype corresponding to a Cactus datatype
static MPI_Datatype
mpi_type (int const cactustype)
{
  int const normaltype = normal_type (cactustype);
  switch (normaltype) {
  case CCTK_VARIABLE_BYTE: return CACTUS_MPI_BYTE;
#ifdef HAVE_CCTK_INT1
  case CCTK_VARIABLE_INT1: return CACTUS_MPI_INT1;
#endif
#ifdef HAVE_CCTK_INT2
  case CCTK_VARIABLE_INT2: return CACTUS_MPI_INT2;
#endif
#ifdef HAVE_CCTK_INT4
  case CCTK_VARIABLE_INT4: return CACTUS_MPI_INT4;
#endif
#ifdef HAVE_CCTK_INT8
  case CCTK_VARIABLE_INT8: return CACTUS_MPI_INT8;
#endif
// #ifdef HAVE_CCTK_INT16
//   case CCTK_VARIABLE_INT16: return CACTUS_MPI_INT16;
// #endif
#ifdef HAVE_CCTK_REAL4
  case CCTK_VARIABLE_REAL4: return CACTUS_MPI_REAL4;
  case CCTK_VARIABLE_COMPLEX8: return CACTUS_MPI_COMPLEX8;
#endif
#ifdef HAVE_CCTK_REAL8
  case CCTK_VARIABLE_REAL8: return CACTUS_MPI_REAL8;
  case CCTK_VARIABLE_COMPLEX16: return CACTUS_MPI_COMPLEX16;
#endif
#ifdef HAVE_CCTK_REAL16
  case CCTK_VARIABLE_REAL16: return CACTUS_MPI_REAL16;
  case CCTK_VARIABLE_COMPLEX32: return CACTUS_MPI_COMPLEX32;
#endif
  }
  CCTK_ERROR("Unsupported datatype");
}



////////////////////////////////////////////////////////////////////////////////



struct bbox {
  int off;                      // offset
  int len;                      // length
  int alen;                     // allocated length
  int str;                      // stride
};

struct arrays {
  bbox global;                  // global region (all processors)
  bbox local;                   // this processor's region
  bbox active;                  // ???
  bbox slab;                    // the slab that should be transferred
};

struct xfer {
  arrays src;                   // source region
  arrays dst;                   // destination region
  int xpose;                    // exchange axes
  int flip;                     // change directions of axes
};



static inline int
roundup (int const x, int const y)
{
  assert (x >= 0);
  assert (y > 0);
  return (x + y - 1) / y * y;
}



static void
bbox_print (bbox const * restrict const bbox)
{
  assert (bbox);
  printf
    ("[%d:%d:%d]",
     bbox->off, bbox->off + (bbox->len - 1) * bbox->str, bbox->str);
}

static void
bbox_check (bbox const * restrict const bbox)
{
  assert (bbox);
  assert (bbox->len >= 0);
  assert (bbox->alen >= bbox->len);
  assert (bbox->str > 0);
}

static void
global2bbox (slabinfo const * restrict const slab,
             bbox           * restrict const bbox)
{
  assert (slab);
  assert (bbox);
  assert (slab->gsh >= 0);
  bbox->off = 0;
  bbox->len = slab->gsh;
  bbox->alen = slab->gsh;
  bbox->str = 1;
  bbox_check (bbox);
}

static void
local2bbox (slabinfo const * restrict const slab,
            bbox           * restrict const bbox)
{
  assert (slab);
  assert (bbox);
  assert (slab->lbnd >= 0);
  assert (slab->lsh >= 0);
  assert (slab->ash >= slab->lsh);
  assert (slab->lbnd + slab->lsh <= slab->gsh);
  bbox->off = slab->lbnd;
  bbox->len = slab->lsh;
  bbox->alen = slab->ash;
  bbox->str = 1;
  bbox_check (bbox);
}

static void
active2bbox (slabinfo const * restrict const slab,
             bbox           * restrict const bbox,
             int                       const useghosts)
{
  assert (slab);
  assert (bbox);
  assert (useghosts == 0 or useghosts == 1);
  assert (slab->lbnd >= 0);
  assert (slab->lsh >= 0);
  assert (slab->ash >= slab->lsh);
  assert (slab->lbnd + slab->lsh <= slab->gsh);
  assert (slab->lbbox == 0 or slab->lbbox == 1);
  assert (slab->ubbox == 0 or slab->ubbox == 1);
  assert (slab->nghostzones >= 0);
  int const nlghostzones = slab->lbbox or useghosts ? 0 : slab->nghostzones;
  int const nughostzones = slab->ubbox or useghosts ? 0 : slab->nghostzones;
  bbox->off = slab->lbnd + nlghostzones;
  bbox->len = slab->lsh - nlghostzones - nughostzones;
  bbox->alen = slab->ash;
  bbox->str = 1;
  bbox_check (bbox);
}

static void
slab2bbox (slabinfo const * restrict const slab,
           bbox           * restrict const bbox)
{
  assert (slab);
  assert (bbox);
  bbox->off = slab->off;
  bbox->len = slab->len;
  bbox->alen = slab->len;
  bbox->str = slab->str;
  bbox_check (bbox);
}

static int
bbox_iscontained (bbox const * restrict const inner,
                  bbox const * restrict const outer)
{
  bbox_check (inner);
  bbox_check (outer);
  int const inner_last = inner->off + (inner->len - 1) * inner->str;
  int const outer_last = outer->off + (outer->len - 1) * outer->str;
  return inner->off >= outer->off and inner_last <= outer_last;
}

static void
bbox_clip (bbox       * restrict const inner,
           bbox const * restrict const outer)
{
  bbox_check (inner);
  bbox_check (outer);
  int inner_last = inner->off + (inner->len - 1) * inner->str;
  int const outer_last = outer->off + (outer->len - 1) * outer->str;
  if (inner->off < outer->off) {
    inner->off += roundup (outer->off - inner->off, inner->str);
  }
  if (inner_last > outer_last) {
    inner_last -= roundup (inner_last - outer_last, inner->str);
  }
  assert ((inner_last - inner->off) % inner->str == 0);
  if (inner_last >= inner->off) {
    inner->len = (inner_last - inner->off + inner->str) / inner->str;
  } else {
    inner->len = 0;
  }
  bbox_check (inner);
}

// ydst = xdst + Flip (ysrc - xsrc)]
// This function is its own inverse.
static void
bbox_xform (bbox       * restrict const ydst,
            bbox const * restrict const ysrc,
            bbox const * restrict const xdst,
            bbox const * restrict const xsrc,
            int                   const flip)
{
  assert (ydst);
  bbox_check (ysrc);
  bbox_check (xdst);
  bbox_check (xsrc);
  assert (ysrc->str == xsrc->str);
  /* int const xsrc_last = xsrc->off + (xsrc->len - 1) * xsrc->str; */
  int const xdst_last = xdst->off + (xdst->len - 1) * xdst->str;
  int const ysrc_last = ysrc->off + (ysrc->len - 1) * ysrc->str;
  ydst->str = xdst->str;
  assert ((ysrc->off - xsrc->off) % ysrc->str == 0);
  ydst->off = xdst->off + (ysrc->off - xsrc->off) / ysrc->str * ydst->str;
  int ydst_last =
    xdst->off + (ysrc_last - xsrc->off) / ysrc->str * ydst->str;
  if (flip) {
    int const off = ydst->off;
    int const last = ydst_last;
    ydst->off = xdst->off + xdst_last - last;
    ydst_last = xdst_last - (off - xdst->off);
  }
  assert ((ysrc_last - xsrc->off) % ysrc->str == 0);
  assert (ydst_last - ydst->off + ydst->str >= 0);
  ydst->len = (ydst_last - ydst->off + ydst->str) / ydst->str;
  ydst->alen = xdst->alen;
  bbox_check (ydst);
}



extern "C"
void
print_slabinfo (FILE           *          const out,
                slabinfo const * restrict const slabinfo)
{
  fprintf (out, "    gsh: %d\n", slabinfo->gsh);
  fprintf (out, "    lbnd: %d, lsh: %d\n", slabinfo->lbnd, slabinfo->lsh);
  fprintf (out, "    ash: %d\n", slabinfo->ash);
  fprintf (out, "    lbbox: %d, ubbox: %d, nghostzones: %d\n",
           slabinfo->lbbox, slabinfo->ubbox, slabinfo->nghostzones);
  fprintf (out, "    off: %d, str: %d, len: %d\n",
           slabinfo->off, slabinfo->str, slabinfo->len);
}

extern "C"
void
print_xferinfo (FILE           *          const out,
                xferinfo const * restrict const xferinfo)
{
  fprintf (out, "  src:\n");
  print_slabinfo (out, & xferinfo->src);
  fprintf (out, "  dst:\n");
  print_slabinfo (out, & xferinfo->dst);
  fprintf (out, "  xpose: %d\n", xferinfo->xpose);
  fprintf (out, "  flip: %d\n", xferinfo->flip);
}                

// workhorse routines for the actual copying transposing, and flipping
// of data
template<typename T>
inline void
copy_data (const vector<xfer> &info,
           const vector<bbox> &srcdetail,
           const vector<int> &srcoffset,
           const vector<int> &srcelems,
           vector<char> &srcdata,
           void  const * restrict const * restrict const srcptrs,
           const int n,
           const vector<int> &varis,
           const int nvaris,
           const int xpose_x,
           const int xpose_y,
           const int xpose_z)
{
  assert (srcptrs);
  
  int const srcoffi = info[0].src.local.off;
  int const srcoffj = info[1].src.local.off;
  int const srcoffk = info[2].src.local.off;
  
  int const srcleni = info[0].src.local.len;
  int const srclenj = info[1].src.local.len;
  int const srclenk = info[2].src.local.len;
  
  int const srcaleni = info[0].src.local.alen;
  int const srcalenj = info[1].src.local.alen;
  int const srcalenk = info[2].src.local.alen;
  
  int const srcdetailoffi = srcdetail[n*SLAB_MAXDIM+0].off;
  int const srcdetailoffj = srcdetail[n*SLAB_MAXDIM+1].off;
  int const srcdetailoffk = srcdetail[n*SLAB_MAXDIM+2].off;
  
  int const srcdetailleni = srcdetail[n*SLAB_MAXDIM+0].len;
  int const srcdetaillenj = srcdetail[n*SLAB_MAXDIM+1].len;
  int const srcdetaillenk = srcdetail[n*SLAB_MAXDIM+2].len;
  
  int const dstdetailleni = srcdetail[n*SLAB_MAXDIM+xpose_x].len;
  int const dstdetaillenj = srcdetail[n*SLAB_MAXDIM+xpose_y].len;
  int const dstdetaillenk = srcdetail[n*SLAB_MAXDIM+xpose_z].len;
  
  if (n==0) assert (srcoffset[n]==0);
  // TODO: This does not take nvaris into account
  // if (n<size-1) assert (srcoffset[n+1]==srcoffset[n]+srcdetailleni*srcdetaillenj*srcdetaillenk);
  
  ifcheck {
    const int xpose[SLAB_MAXDIM] = {xpose_x, xpose_y, xpose_z};
    for (int i = 0; i < SLAB_MAXDIM; ++i) {
      for (int j = i+1; j < SLAB_MAXDIM; ++j) {
        assert(xpose[i] != xpose[j]);
      }
    }
  }
  
  assert (dstdetailleni*dstdetaillenj*dstdetaillenk == srcelems[n]);
  if (srcelems[n] == 0) return;
  
  for (int vari=0; vari<nvaris; ++vari) {
    T * restrict const srcdataptr =
      (T *)&srcdata.front() + srcoffset[n] + vari * srcelems[n];
    T const * restrict const srcptr = (T const *)srcptrs[varis[vari]];
    assert(srcptr);
    
    if (xpose_x==0 and xpose_y==1 and xpose_z==2) {
      // no transposition
      
      CCTK_TimerStartI (timer_copy_in_noxpose);
#     pragma omp parallel
      CCTK_LOOP3(Slab_copy_in_noxpose, i,j,k,
                 0,0,0, srcdetailleni,srcdetaillenj,srcdetaillenk,
                 srcaleni,srcalenj,srcalenk)
      {
        int const srcindi = srcdetailoffi + i - srcoffi;
        int const srcindj = srcdetailoffj + j - srcoffj;
        int const srcindk = srcdetailoffk + k - srcoffk;
        ifcheck assert (srcindi>=0 and srcindi<srcleni);
        ifcheck assert (srcindj>=0 and srcindj<srclenj);
        ifcheck assert (srcindk>=0 and srcindk<srclenk);
        size_t const srcind =
          srcindi + srcaleni * (srcindj + srcalenj * srcindk);
        size_t const bufind =
          i + dstdetailleni * (j + dstdetaillenj * k);
        srcdataptr[bufind] = srcptr[srcind];
      } CCTK_ENDLOOP3(Slab_copy_in_noxpose);
      CCTK_TimerStopI (timer_copy_in_noxpose);
      
    } else if (xpose_x==1 and xpose_y==0 and xpose_z==2) {
      // transpose x and y
      
      CCTK_TimerStartI (timer_copy_in_xposexy);
#     pragma omp parallel
      // Interchange i and j loops
      CCTK_LOOP3(Slab_copy_in_xposexy, j,i,k,
                 0,0,0, srcdetaillenj,srcdetailleni,srcdetaillenk,
                 srcalenj,srcaleni,srcalenk)
      {
        int const srcindi = srcdetailoffi + i - srcoffi;
        int const srcindj = srcdetailoffj + j - srcoffj;
        int const srcindk = srcdetailoffk + k - srcoffk;
        ifcheck assert (srcindi>=0 and srcindi<srcleni);
        ifcheck assert (srcindj>=0 and srcindj<srclenj);
        ifcheck assert (srcindk>=0 and srcindk<srclenk);
        size_t const srcind =
          srcindi + srcaleni * (srcindj + srcalenj * srcindk);
        size_t const bufind =
          j + dstdetailleni * (i + dstdetaillenj * k);
        srcdataptr[bufind] = srcptr[srcind];
      } CCTK_ENDLOOP3(Slab_copy_in_xposexy);
      CCTK_TimerStopI (timer_copy_in_xposexy);
      
    } else {
      // general transposition
      
      CCTK_TimerStartI (timer_copy_in_xposegeneral);
#     pragma omp parallel
      CCTK_LOOP3(Slab_copy_in_xposegeneral, i,j,k,
                 0,0,0, srcdetailleni,srcdetaillenj,srcdetaillenk,
                 srcaleni,srcalenj,srcalenk)
      {
        int ipos[SLAB_MAXDIM];
        ipos[0] = i;
        ipos[1] = j;
        ipos[2] = k;
        int const srcindi = srcdetailoffi + i - srcoffi;
        int const srcindj = srcdetailoffj + j - srcoffj;
        int const srcindk = srcdetailoffk + k - srcoffk;
        ifcheck assert (srcindi>=0 and srcindi<srcleni);
        ifcheck assert (srcindj>=0 and srcindj<srclenj);
        ifcheck assert (srcindk>=0 and srcindk<srclenk);
        size_t const srcind =
          srcindi + srcaleni * (srcindj + srcalenj * srcindk);
        // TODO: rewrite this, moving the index onto dstdetaillen and
        // outside the loop; then collapse all loop specialisations
        size_t const bufind =
          (ipos[xpose_x] + dstdetailleni *
           (ipos[xpose_y] + dstdetaillenj * ipos[xpose_z]));
        srcdataptr[bufind] = srcptr[srcind];
      } CCTK_ENDLOOP3(Slab_copy_in_xposegeneral);
      CCTK_TimerStopI (timer_copy_in_xposegeneral);
      
    }
    
  } // for vari
}

// workhorse routine responsible for the actual copying/flipping of data
template<typename T> inline void
copy_data_back (const vector<xfer> &info,
                const vector<bbox> &dstdetail,
                const vector<int> &dstoffset,
                const vector<int> &dstelems,
                const vector<char> &dstdata,
                void * restrict const * restrict const dstptrs,
                const int n,
                const vector<int> &varis,
                const int nvaris,
                const bool flip_x,
                const bool flip_y,
                const bool flip_z)
{
  assert (dstptrs);
  
  int const dstoffi = info[0].dst.local.off;
  int const dstoffj = info[1].dst.local.off;
  int const dstoffk = info[2].dst.local.off;
  
  int const dstleni = info[0].dst.local.len;
  int const dstlenj = info[1].dst.local.len;
  int const dstlenk = info[2].dst.local.len;
  
  int const dstaleni = info[0].dst.local.alen;
  int const dstalenj = info[1].dst.local.alen;
  int const dstalenk = info[2].dst.local.alen;
  
  int const dstdetailoffi = dstdetail[n*SLAB_MAXDIM+0].off;
  int const dstdetailoffj = dstdetail[n*SLAB_MAXDIM+1].off;
  int const dstdetailoffk = dstdetail[n*SLAB_MAXDIM+2].off;
  
  int const dstdetailleni = dstdetail[n*SLAB_MAXDIM+0].len;
  int const dstdetaillenj = dstdetail[n*SLAB_MAXDIM+1].len;
  int const dstdetaillenk = dstdetail[n*SLAB_MAXDIM+2].len;
  
  assert (dstdetailleni*dstdetaillenj*dstdetaillenk == dstelems[n]);
  if (dstelems[n] == 0) return;
  
  for (int vari=0; vari<nvaris; ++vari) {
    T * restrict const dstptr = (T * restrict)dstptrs[varis[vari]];
    assert (dstptr);
    T const * restrict const dstdataptr =
      (T const *)&dstdata.front() + dstoffset[n] + vari * dstelems[n];
    
    if (not flip_x and not flip_y and not flip_z) {
      // no flipping
      
      CCTK_TimerStartI (timer_copy_back_noflip);
#     pragma omp parallel
      CCTK_LOOP3(Slab_copy_back_noflip, i,j,k,
                 0,0,0, dstdetailleni,dstdetaillenj,dstdetaillenk,
                 dstaleni,dstalenj,dstalenk)
      {
        int const dstindi = dstdetailoffi + i - dstoffi;
        int const dstindj = dstdetailoffj + j - dstoffj;
        int const dstindk = dstdetailoffk + k - dstoffk;
        ifcheck assert (dstindi>=0 and dstindi<dstleni);
        ifcheck assert (dstindj>=0 and dstindj<dstlenj);
        ifcheck assert (dstindk>=0 and dstindk<dstlenk);
        size_t const dstind =
          dstindi + dstaleni * (dstindj + dstalenj * dstindk);
        size_t const bufind =
          i + dstdetailleni * (j + dstdetaillenj * k);
        dstptr[dstind] = dstdataptr[bufind];
      } CCTK_ENDLOOP3(Slab_copy_back_noflip);
      CCTK_TimerStopI (timer_copy_back_noflip);
      
    } else if (flip_x and not flip_y and not flip_z) {
      // flip in x direction
      
      CCTK_TimerStartI (timer_copy_back_flipx);
#     pragma omp parallel
      CCTK_LOOP3(Slab_copy_back_flipx, i,j,k,
                 0,0,0, dstdetailleni,dstdetaillenj,dstdetaillenk,
                 dstaleni,dstalenj,dstalenk)
      {
        int const dstindi = dstdetailoffi + (dstdetailleni - 1 - i) - dstoffi;
        int const dstindj = dstdetailoffj + j - dstoffj;
        int const dstindk = dstdetailoffk + k - dstoffk;
        ifcheck assert (dstindi>=0 and dstindi<dstleni);
        ifcheck assert (dstindj>=0 and dstindj<dstlenj);
        ifcheck assert (dstindk>=0 and dstindk<dstlenk);
        size_t const dstind =
          dstindi + dstaleni * (dstindj + dstalenj * dstindk);
        size_t const bufind =
          i + dstdetailleni * (j + dstdetaillenj * k);
        dstptr[dstind] = dstdataptr[bufind];
      } CCTK_ENDLOOP3(Slab_copy_back_flipx);
      CCTK_TimerStopI (timer_copy_back_flipx);
      
    } else if (not flip_x and flip_y and not flip_z) {
      // flip in y direction
      
      CCTK_TimerStartI (timer_copy_back_flipy);
#     pragma omp parallel
      CCTK_LOOP3(Slab_copy_back_flipy, i,j,k,
                 0,0,0, dstdetailleni,dstdetaillenj,dstdetaillenk,
                 dstaleni,dstalenj,dstalenk)
      {
        int const dstindi = dstdetailoffi + i - dstoffi;
        int const dstindj = dstdetailoffj + (dstdetaillenj - 1 - j) - dstoffj;
        int const dstindk = dstdetailoffk + k - dstoffk;
        ifcheck assert (dstindi>=0 and dstindi<dstleni);
        ifcheck assert (dstindj>=0 and dstindj<dstlenj);
        ifcheck assert (dstindk>=0 and dstindk<dstlenk);
        size_t const dstind =
          dstindi + dstaleni * (dstindj + dstalenj * dstindk);
        size_t const bufind =
          i + dstdetailleni * (j + dstdetaillenj * k);
        dstptr[dstind] = dstdataptr[bufind];
      } CCTK_ENDLOOP3(Slab_copy_back_flipy);
      CCTK_TimerStopI (timer_copy_back_flipy);
      
    } else if (flip_x and flip_y and not flip_z) {
      // flip in both x and y direction
      
      CCTK_TimerStartI (timer_copy_back_flipxy);
#     pragma omp parallel
      CCTK_LOOP3(Slab_copy_back_flipxy, i,j,k,
                 0,0,0, dstdetailleni,dstdetaillenj,dstdetaillenk,
                 dstaleni,dstalenj,dstalenk)
      {
        int const dstindi = dstdetailoffi + (dstdetailleni - 1 - i) - dstoffi;
        int const dstindj = dstdetailoffj + (dstdetaillenj - 1 - j) - dstoffj;
        int const dstindk = dstdetailoffk + k - dstoffk;
        ifcheck assert (dstindi>=0 and dstindi<dstleni);
        ifcheck assert (dstindj>=0 and dstindj<dstlenj);
        ifcheck assert (dstindk>=0 and dstindk<dstlenk);
        size_t const dstind =
          dstindi + dstaleni * (dstindj + dstalenj * dstindk);
        size_t const bufind =
          i + dstdetailleni * (j + dstdetaillenj * k);
        dstptr[dstind] = dstdataptr[bufind];
      } CCTK_ENDLOOP3(Slab_copy_back_flipxy);
      CCTK_TimerStopI (timer_copy_back_flipxy);
      
    } else {
      // general flipping
      
      CCTK_TimerStartI (timer_copy_back_flipgeneral);
#     pragma omp parallel
      CCTK_LOOP3(Slab_copy_back_flipgeneral, i,j,k,
                 0,0,0, dstdetailleni,dstdetaillenj,dstdetaillenk,
                 dstaleni,dstalenj,dstalenk)
      {
        // TODO: rewrite this, moving the flipping onto dstlen and
        // outside the loop; then collapse all loop specialisations
        int const dstindi =
          dstdetailoffi + (flip_x ? dstdetailleni - 1 - i : i) - dstoffi;
        int const dstindj =
          dstdetailoffj + (flip_y ? dstdetaillenj - 1 - j : j) - dstoffj;
        int const dstindk =
          dstdetailoffk + (flip_z ? dstdetaillenk - 1 - k : k) - dstoffk;
        ifcheck assert (dstindi>=0 and dstindi<dstleni);
        ifcheck assert (dstindj>=0 and dstindj<dstlenj);
        ifcheck assert (dstindk>=0 and dstindk<dstlenk);
        size_t const dstind =
          dstindi + dstaleni * (dstindj + dstalenj * dstindk);
        size_t const bufind =
          i + dstdetailleni * (j + dstdetaillenj * k);
        dstptr[dstind] = dstdataptr[bufind];
      } CCTK_ENDLOOP3(Slab_copy_back_flipgeneral);
      CCTK_TimerStopI (timer_copy_back_flipgeneral);
      
    }
    
  } // for vari
}



struct slabsetup {
  MPI_Comm comm;
  vector<xfer> info;
  vector<xfer> allinfo;
  vector<bbox> srcdetail, dstdetail;
  size_t srclentot, dstlentot;
};



extern "C"
slabsetup *
Slab_MultiTransfer_Init
(cGH       const* restrict const cctkGH,
 int                       const dim,
 xferinfo  const* restrict const xferinfo,
 int                       const options)
{
  DECLARE_CCTK_PARAMETERS;
  
  // Check arguments
  check (cctkGH);
  check (dim >= 0);
  check (xferinfo);
  
  
  
  CCTK_TimerStartI (timer_init);
  
  slabsetup * restrict const slabsetup = new struct slabsetup;
  
  
  
  bool useghosts;
  {
    CCTK_INT tmp;
    int const iret = Util_TableGetInt (options, &tmp, "useghosts");
    if (iret == 1) {
      // There was an entry, use it
      useghosts = tmp;
    } else if (iret == UTIL_ERROR_BAD_HANDLE or
               iret == UTIL_ERROR_TABLE_NO_SUCH_KEY)
    {
      // There was no entry, use a default
      useghosts = false;
    } else {
      // Something went wrong, abort
      check (0);
    }
  }
  
  check (dim <= SLAB_MAXDIM);
  vector<xfer>& info = slabsetup->info;
  info.resize (SLAB_MAXDIM);
  for (int d=0; d<dim; ++d) {
    global2bbox (&xferinfo[d].src, &info[d].src.global);
    local2bbox  (&xferinfo[d].src, &info[d].src.local);
    active2bbox (&xferinfo[d].src, &info[d].src.active, useghosts);
    slab2bbox   (&xferinfo[d].src, &info[d].src.slab);
    check (bbox_iscontained (&info[d].src.active, &info[d].src.local));
    check (bbox_iscontained (&info[d].src.local, &info[d].src.global));
    
    global2bbox (&xferinfo[d].dst, &info[d].dst.global);
    local2bbox  (&xferinfo[d].dst, &info[d].dst.local);
    active2bbox (&xferinfo[d].dst, &info[d].dst.active, 1); // fill ghosts
    slab2bbox   (&xferinfo[d].dst, &info[d].dst.slab);
    check (bbox_iscontained (&info[d].dst.active, &info[d].dst.local));
    check (bbox_iscontained (&info[d].dst.local, &info[d].dst.global));
    
    info[d].xpose = xferinfo[d].xpose;
    check (info[d].xpose >= 0 and info[d].xpose < dim);
    info[d].flip = xferinfo[d].flip;
    check (info[d].flip == 0 or info[d].flip == 1);
  }
  for (int d=dim; d<SLAB_MAXDIM; ++d) {
    static bbox const fake_bbox = { 0, 1, 1 };
    static arrays const fake_arrays =
      { { 0, 1, 1 }, { 0, 1, 1 }, { 0, 1, 1 }, { 0, 1, 1 } };
    
    bbox_check (&fake_bbox);
    
    info[d].src = fake_arrays;
    check (bbox_iscontained (&info[d].src.active, &info[d].src.local));
    check (bbox_iscontained (&info[d].src.local, &info[d].src.global));
    
    info[d].dst = fake_arrays;
    check (bbox_iscontained (&info[d].dst.active, &info[d].dst.local));
    check (bbox_iscontained (&info[d].dst.local, &info[d].dst.global));
    
    info[d].xpose = d;
    check (info[d].xpose >= 0 and info[d].xpose < SLAB_MAXDIM);
    info[d].flip = 0;
    check (info[d].flip == 0 or info[d].flip == 1);
  }
  
  ifcheck {
    ifdebug printf ("srcinfo:\n");
    for (int d=0; d<SLAB_MAXDIM; ++d) {
      printf ("   src.global d=%d ", d);
      bbox_print (&info[d].src.global);
      printf ("\n");
      printf ("   src.local  d=%d ", d);
      bbox_print (&info[d].src.local);
      printf ("\n");
      printf ("   src.active d=%d ", d);
      bbox_print (&info[d].src.active);
      printf ("\n");
      printf ("   src.slab   d=%d ", d);
      bbox_print (&info[d].src.slab);
      printf ("\n");
    }
    ifdebug printf ("dstinfo:\n");
    for (int d=0; d<SLAB_MAXDIM; ++d) {
      printf ("   dst.global d=%d ", d);
      bbox_print (&info[d].dst.global);
      printf ("\n");
      printf ("   dst.local  d=%d ", d);
      bbox_print (&info[d].dst.local);
      printf ("\n");
      printf ("   dst.active d=%d ", d);
      bbox_print (&info[d].dst.active);
      printf ("\n");
      printf ("   dst.slab   d=%d ", d);
      bbox_print (&info[d].dst.slab);
      printf ("\n");
    }
    ifdebug printf ("info:\n");
    for (int d=0; d<SLAB_MAXDIM; ++d) {
      printf ("   xpose      d=%d %d\n", d, info[d].xpose);
      printf ("   flip       d=%d %d\n", d, info[d].flip);
    }
  }
  
  {
    bool iflag[SLAB_MAXDIM];
    for (int d=0; d<SLAB_MAXDIM; ++d) {
      iflag[d] = false;
    }
    for (int d=0; d<SLAB_MAXDIM; ++d) {
      assert (not iflag[info[d].xpose]);
      iflag[info[d].xpose] = true;
    }
    for (int d=0; d<SLAB_MAXDIM; ++d) {
      assert (iflag[d]);
    }
    // Allow non-contributing processors to be non-knowledgeable
    for (int d=0; d<SLAB_MAXDIM; ++d) {
      if (info[info[d].xpose].src.slab.len and info[d].dst.slab.len > 0) {
        assert (info[info[d].xpose].src.slab.len == info[d].dst.slab.len);
      }
    }
  }
  
  size_t& srclentot = slabsetup->srclentot;
  size_t& dstlentot = slabsetup->dstlentot;
  srclentot = 1;
  dstlentot = 1;
  for (int d=0; d<SLAB_MAXDIM; ++d) {
    srclentot *= info[d].src.local.len;
    dstlentot *= info[d].dst.local.len;
  }
  
  
  
  MPI_Comm& comm = slabsetup->comm;
  {
    CCTK_POINTER_TO_CONST tmp1;
    int const iret1 = Util_TableGetPointerToConst (options, &tmp1, "comm");
    if (iret1 == 1) {
      // There was an entry, use it
      comm = * (MPI_Comm const *) tmp1;
    } else if (iret1 == UTIL_ERROR_TABLE_WRONG_DATA_TYPE) {
      // Entry has wrong type, fall back
      CCTK_POINTER tmp2;
      int const iret2 = Util_TableGetPointer (options, &tmp2, "comm");
      if (iret2 == 1) {
        // There was an entry, use it
        comm = * (MPI_Comm const *) tmp2;
      } else {
        // Something went wrong, abort
        check (0);
      }
    } else if (iret1 == UTIL_ERROR_BAD_HANDLE or
               iret1 == UTIL_ERROR_TABLE_NO_SUCH_KEY)
    {
      // There was no entry, use a default
      comm = get_mpi_comm (cctkGH);
    } else {
      // Something went wrong, abort
      check (0);
    }
  }
  
  ifcheck {
    ifdebug fflush (stdout);
    MPI_Barrier (comm);
  }
  
  int size, rank;
  MPI_Comm_size (comm, &size);
  MPI_Comm_rank (comm, &rank);
  
  ifcheck {
    static int count = 424242;
    int mycount = count;
    ifdebug fflush (stdout);
    MPI_Bcast (&mycount, 1, MPI_INT, 0, comm);
    assert (mycount == count);
    ++ count;
  }
  
  
  
  vector<xfer>& allinfo = slabsetup->allinfo;
  allinfo.resize (size * SLAB_MAXDIM);
  {
    int const info_nints = sizeof(xfer) / sizeof(int);
    ifdebug fflush (stdout);
    MPI_Allgather
      (&info.front(),    SLAB_MAXDIM * info_nints, MPI_INT,
       &allinfo.front(), SLAB_MAXDIM * info_nints, MPI_INT, comm);
  }
  
  for (int n = 0; n < size; ++n) {
    for (int d=0; d<SLAB_MAXDIM; ++d) {
      // Allow non-contributing processors to be non-knowledgeable
      if (allinfo[n*SLAB_MAXDIM+d].src.slab.len > 0 and
          info[d].src.slab.len > 0)
      {
        assert
          (allinfo[n*SLAB_MAXDIM+d].src.global.off == info[d].src.global.off);
        assert
          (allinfo[n*SLAB_MAXDIM+d].src.global.len == info[d].src.global.len);
        assert
          (allinfo[n*SLAB_MAXDIM+d].src.global.str == info[d].src.global.str);
        assert
          (allinfo[n*SLAB_MAXDIM+d].src.local.str == info[d].src.local.str);
        assert
          (allinfo[n*SLAB_MAXDIM+d].src.active.str == info[d].src.active.str);
        // 2003-03-01 eschnett: I don't know why the following should
        // be necessary
        assert
          (allinfo[n*SLAB_MAXDIM+d].src.slab.str == info[d].src.slab.str);
      }
      if (allinfo[n*SLAB_MAXDIM+d].dst.slab.len > 0 and
          info[d].dst.slab.len > 0)
      {
        assert
          (allinfo[n*SLAB_MAXDIM+d].dst.global.off == info[d].dst.global.off);
        assert
          (allinfo[n*SLAB_MAXDIM+d].dst.global.len == info[d].dst.global.len);
        assert
          (allinfo[n*SLAB_MAXDIM+d].dst.global.str == info[d].dst.global.str);
        assert
          (allinfo[n*SLAB_MAXDIM+d].dst.local.str == info[d].dst.local.str);
        assert
          (allinfo[n*SLAB_MAXDIM+d].dst.active.str == info[d].dst.active.str);
        assert
          (allinfo[n*SLAB_MAXDIM+d].dst.slab.str == info[d].dst.slab.str);
      }
      assert (allinfo[n*SLAB_MAXDIM+d].xpose == info[d].xpose);
      assert (allinfo[n*SLAB_MAXDIM+d].flip == info[d].flip);
    }
  }
  
  
  
  vector<bbox>& srcdetail = slabsetup->srcdetail;
  srcdetail.resize (size * SLAB_MAXDIM);
  for (int n = 0; n < size; ++n) {
    ifdebug printf ("srcdetail n=%d:\n", n);
    for (int d=0; d<SLAB_MAXDIM; ++d) {
      srcdetail[n*SLAB_MAXDIM+d] = allinfo[n*SLAB_MAXDIM+d].src.slab;
      ifdebug printf ("   src.slab                d=%d ", d);
      ifdebug bbox_print (&srcdetail[n*SLAB_MAXDIM+d]);
      ifdebug printf ("\n");
      bbox_clip (&srcdetail[n*SLAB_MAXDIM+d], &info[d].src.active);
      ifdebug printf ("   clipped with src.active d=%d ", d);
      ifdebug bbox_print (&srcdetail[n*SLAB_MAXDIM+d]);
      ifdebug printf ("\n");
    }
    for (int d=0; d<SLAB_MAXDIM; ++d) {
      bbox whereto;
      bbox wherefrom;
      whereto = allinfo[n*SLAB_MAXDIM+d].dst.slab;
      ifdebug printf ("   dst.slab                d=%d ", info[d].xpose);
      ifdebug bbox_print (&whereto);
      ifdebug printf ("\n");
      bbox_clip (&whereto, &allinfo[n*SLAB_MAXDIM+d].dst.active);
      ifdebug printf ("   whereto                 d=%d ", info[d].xpose);
      ifdebug bbox_print (&whereto);
      ifdebug printf ("\n");
      bbox_xform
	(&wherefrom, &whereto,
	 &allinfo[n*SLAB_MAXDIM+info[d].xpose].src.slab,
         &allinfo[n*SLAB_MAXDIM+d].dst.slab,
	 info[d].flip);
      ifdebug printf ("   wherefrom               d=%d ", info[d].xpose);
      ifdebug bbox_print (&wherefrom);
      ifdebug printf ("\n");
      bbox_clip (&srcdetail[n*SLAB_MAXDIM+info[d].xpose], &wherefrom);
      ifdebug printf ("   clipped with wherefrom  d=%d ", info[d].xpose);
      ifdebug bbox_print (&srcdetail[n*SLAB_MAXDIM+info[d].xpose]);
      ifdebug printf ("\n");
    }
  }
  
  
  
  vector<bbox>& dstdetail = slabsetup->dstdetail;
  dstdetail.resize (size * SLAB_MAXDIM);
  for (int n = 0; n < size; ++n) {
    ifdebug printf ("dstdetail n=%d:\n", n);
    for (int d=0; d<SLAB_MAXDIM; ++d) {
      // dstdetail[n*SLAB_MAXDIM+d] = allinfo[n*SLAB_MAXDIM+d].dst.slab;
      dstdetail[n*SLAB_MAXDIM+d] = info[d].dst.slab;
      ifdebug printf ("   dst.slab                d=%d ", d);
      ifdebug bbox_print (&dstdetail[n*SLAB_MAXDIM+d]);
      ifdebug printf ("\n");
      bbox_clip (&dstdetail[n*SLAB_MAXDIM+d], &info[d].dst.active);
      ifdebug printf ("   clipped with dst.active d=%d ", d);
      ifdebug bbox_print (&dstdetail[n*SLAB_MAXDIM+d]);
      ifdebug printf ("\n");
    }
    for (int d=0; d<SLAB_MAXDIM; ++d) {
      bbox wherefrom;
      bbox whereto;
      // wherefrom = allinfo[n*SLAB_MAXDIM+info[d].xpose].src.slab;
      wherefrom = info[info[d].xpose].src.slab;
      ifdebug printf ("   src.slab                d=%d ", d);
      // ifdebug bbox_print (&dstdetail[n*SLAB_MAXDIM+d]);
      ifdebug bbox_print (&wherefrom);
      ifdebug printf ("\n");
      bbox_clip (&wherefrom, &allinfo[n*SLAB_MAXDIM+info[d].xpose].src.active);
      ifdebug printf ("   wherefrom               d=%d ", d);
      // ifdebug bbox_print (&dstdetail[n*SLAB_MAXDIM+d]);
      ifdebug bbox_print (&wherefrom);
      ifdebug printf ("\n");
      bbox_xform
	(&whereto, &wherefrom,
	 &allinfo[n*SLAB_MAXDIM+d].dst.slab,
         // &allinfo[n*SLAB_MAXDIM+info[d].xpose].src.slab,
         &info[info[d].xpose].src.slab,
	 info[d].flip);
      ifdebug printf ("   whereto                 d=%d ", d);
      // ifdebug bbox_print (&dstdetail[n*SLAB_MAXDIM+d]);
      ifdebug bbox_print (&whereto);
      ifdebug printf ("\n");
      bbox_clip (&dstdetail[n*SLAB_MAXDIM+d], &whereto);
      ifdebug printf ("   clipped with whereto    d=%d ", d);
      ifdebug bbox_print (&dstdetail[n*SLAB_MAXDIM+d]);
      ifdebug printf ("\n");
    }
  }
  
  
  
  CCTK_TimerStopI (timer_init);
  
  return slabsetup;
}



extern "C"
int
Slab_MultiTransfer_Apply
(cGH       const                  * restrict const cctkGH,
 slabsetup const                  * restrict const slabsetup,
 int                                         const nvars,
 int       const                  * restrict const srctypes,
 void      const * restrict const * restrict const srcptrs,
 int       const                  * restrict const dsttypes,
 void            * restrict const * restrict const dstptrs)
{
  DECLARE_CCTK_PARAMETERS;
  
  // Check arguments
  check (cctkGH);
  check (slabsetup);
  check (nvars >= 0);
  check (nvars==0 or srctypes);
  for (int var=0; var<nvars; ++var) check (srctypes[var] >= 0);
  check (nvars==0 or srcptrs);
  size_t const& srclentot = slabsetup->srclentot;
  for (int var=0; var<nvars; ++var) if (srclentot > 0) assert (srcptrs[var]);
  check (nvars==0 or dsttypes);
  for (int var=0; var<nvars; ++var) check (dsttypes[var] >= 0);
  check (nvars==0 or dstptrs);
  size_t const& dstlentot = slabsetup->dstlentot;
  for (int var=0; var<nvars; ++var) if (dstlentot > 0) assert (dstptrs[var]);
  
  if (nvars==0) return 0;
  
  
  
  CCTK_TimerStartI (timer_apply);
  
  
  
  MPI_Comm const& comm = slabsetup->comm;
  
  ifcheck {
    ifdebug fflush (stdout);
    MPI_Barrier (comm);
  }
  
  int size, rank;
  MPI_Comm_size (comm, &size);
  MPI_Comm_rank (comm, &rank);
  
  ifcheck {
    static int count = 424242;
    int mycount = count;
    ifdebug fflush (stdout);
    MPI_Bcast (&mycount, 1, MPI_INT, 0, comm);
    assert (mycount == count);
    ++ count;
  }
  
  
  
  vector<xfer> const& info = slabsetup->info;
  vector<xfer> const& allinfo = slabsetup->allinfo;
  vector<bbox> const& srcdetail = slabsetup->srcdetail;
  vector<bbox> const& dstdetail = slabsetup->dstdetail;
  
  vector<int> srcelems (size);
  vector<int> srccount (size);
  vector<int> srcoffset (size + 1);
  
  vector<int> dstelems (size);
  vector<int> dstcount (size);
  vector<int> dstoffset (size + 1);
  
  
  
  int nvartypes = 0;
  vector<int> vartypes (nvars);
  vector<int> vartypecount (nvars);
  for (int var=0; var<nvars; ++var) {
    int const srctype = srctypes[var];
    int const dsttype = dsttypes[var];
    check (srctype == dsttype);
    int vartypei;
    for (vartypei=0; vartypei<nvartypes; ++vartypei) {
      if (srctype == vartypes[vartypei]) break;
    }
    if (vartypei>=nvartypes) {
      vartypes[nvartypes] = srctype;
      vartypecount[nvartypes] = 0;
      ++nvartypes;
    }
    assert (vartypei<nvartypes);
    ++ vartypecount[vartypei];
  }
  
  for (int vartypei=0; vartypei<nvartypes; ++vartypei) {
    int const vartype = vartypes[vartypei];
    int nvaris = 0;
    vector<int> varis (nvars);
    for (int var=0; var<nvars; ++var) {
      if (srctypes[var] == vartype) {
        varis[nvaris] = var;
        ++nvaris;
      }
    }
    assert (nvaris == vartypecount[vartypei]);
    assert (nvaris > 0);
    
    int const vartypesize = CCTK_VarTypeSize (vartype);
    check (vartypesize > 0);
    MPI_Datatype const vardatatype = mpi_type (vartype);
    
    
    
    srcoffset[0] = 0;
    for (int n = 0; n < size; ++n) {
      srcelems[n] = 1;
      for (int d=0; d<SLAB_MAXDIM; ++d) {
        srcelems[n] *= srcdetail[n*SLAB_MAXDIM+d].len;
      }
      srccount[n] = nvaris * srcelems[n];
      ifdebug printf
        ("srccnt n=%d offset=%d count=%d\n", n, srcoffset[n], srccount[n]);
      srcoffset[n+1] = srcoffset[n] + srccount[n];
    }
    vector<char> srcdata (srcoffset[size] * vartypesize);
    ifcheck {
      if (vartype == CCTK_VARIABLE_REAL) {
        CCTK_REAL * restrict const srcdataptr = (CCTK_REAL *)&srcdata.front();
        CCTK_REAL marker;
        memset (&marker, POISON_VALUE, sizeof marker);
        for (int i = 0; i < srcoffset[size]; ++i) {
          memcpy (&srcdataptr[i], &marker, sizeof marker);
        }
      }
    }
    
    
    
    dstoffset[0] = 0;
    for (int n = 0; n < size; ++n) {
      dstelems[n] = 1;
      for (int d=0; d<SLAB_MAXDIM; ++d) {
        dstelems[n] *= dstdetail[n*SLAB_MAXDIM+d].len;
      }
      dstcount[n] = nvaris * dstelems[n];
      ifdebug printf
        ("dstcnt n=%d offset=%d count=%d\n", n, dstoffset[n], dstcount[n]);
      dstoffset[n+1] = dstoffset[n] + dstcount[n];
    }
    vector<char> dstdata (dstoffset[size] * vartypesize);
    ifcheck {
      if (vartype == CCTK_VARIABLE_REAL) {
        CCTK_REAL * restrict const dstdataptr = (CCTK_REAL *)&dstdata.front();
        CCTK_REAL marker;
        memset (&marker, POISON_VALUE, sizeof marker);
        for (int i = 0; i < dstoffset[size]; ++i) {
          memcpy (&dstdataptr[i], &marker, sizeof marker);
        }
      }
    }
    
    check (srccount[rank] == dstcount[rank]);
    
    
    
    ifcheck {
      vector<int> src2count (size);
      vector<int> dst2count (size);
      ifdebug fflush (stdout);
      MPI_Alltoall
        (&srccount.front(), 1, MPI_INT, &src2count.front(), 1, MPI_INT, comm);
      MPI_Alltoall
        (&dstcount.front(), 1, MPI_INT, &dst2count.front(), 1, MPI_INT, comm);
      for (int n = 0; n < size; ++n) {
        check (src2count[n] == dstcount[n]);
        check (dst2count[n] == srccount[n]);
      }
    }



    CCTK_TimerStartI (timer_copy_in);
    
    for (int n = 0; n < size; ++n) {
      check (SLAB_MAXDIM == 3);
      
      if (srcdetail[n*SLAB_MAXDIM  ].str==1 and
          srcdetail[n*SLAB_MAXDIM+1].str==1 and
          srcdetail[n*SLAB_MAXDIM+2].str==1 and
          vartype == CCTK_VARIABLE_REAL)
      {
        // Optimised for stride 1 and CCTK_REAL
        copy_data<CCTK_REAL>
          (info, srcdetail, srcoffset, srcelems, srcdata, srcptrs, 
           n, varis, nvaris, info[0].xpose, info[1].xpose, info[2].xpose);
      } else if (srcdetail[n*SLAB_MAXDIM  ].str==1 and
                 srcdetail[n*SLAB_MAXDIM+1].str==1 and
                 srcdetail[n*SLAB_MAXDIM+2].str==1 and
                 vartype == CCTK_VARIABLE_COMPLEX)
      {
        // Optimised for stride 1 and CCTK_COMPLEX
        copy_data<CCTK_COMPLEX>
          (info, srcdetail, srcoffset, srcelems, srcdata, srcptrs, 
           n, varis, nvaris, info[0].xpose, info[1].xpose, info[2].xpose);
      } else if (srcdetail[n*SLAB_MAXDIM  ].str==1 and
                 srcdetail[n*SLAB_MAXDIM+1].str==1 and
                 srcdetail[n*SLAB_MAXDIM+2].str==1 and
                 vartype == CCTK_VARIABLE_INT)
      {
        // Optimised for stride 1 and CCTK_INT
        copy_data<CCTK_INT>
          (info, srcdetail, srcoffset, srcelems, srcdata, srcptrs, 
           n, varis, nvaris, info[0].xpose, info[1].xpose, info[2].xpose);
      } else {
        // Generic, unoptimised version
        CCTK_TimerStartI (timer_copy_in_general);
        
        int const srcdetailleni = srcdetail[n*SLAB_MAXDIM+info[0].xpose].len;
        int const srcdetaillenj = srcdetail[n*SLAB_MAXDIM+info[1].xpose].len;
        int const srcdetaillenk = srcdetail[n*SLAB_MAXDIM+info[2].xpose].len;
        
        for (int vari=0; vari<nvaris; ++vari) {
          char * restrict const srcdataptr =
            &srcdata[vartypesize * (srcoffset[n] + vari * srcelems[n])];
          char const * restrict const srcptr =
            (char const *)srcptrs[varis[vari]];
          
#         pragma omp parallel for
          for (int k = 0; k < srcdetaillenk; ++k) {
            for (int j = 0; j < srcdetaillenj; ++j) {
              for (int i = 0; i < srcdetailleni; ++i) {
                int ipos[SLAB_MAXDIM];
                ipos[0] = i;
                ipos[1] = j;
                ipos[2] = k;
                int srcipos[SLAB_MAXDIM];
                int bufipos[SLAB_MAXDIM];
                for (int d=0; d<SLAB_MAXDIM; ++d) {
                  int const c = info[d].xpose;
                  srcipos[c] = srcdetail[n*SLAB_MAXDIM+c].off + ipos[d] * srcdetail[n*SLAB_MAXDIM+c].str;
                  assert (srcipos[c] >= info[c].src.local.off and
                          srcipos[c] < info[c].src.local.off + info[c].src.local.len);
                  assert (srcipos[c] >= allinfo[n*SLAB_MAXDIM+c].src.slab.off and
                          srcipos[c] <= allinfo[n*SLAB_MAXDIM+c].src.slab.off + (allinfo[n*SLAB_MAXDIM+c].src.slab.len - 1) * allinfo[n*SLAB_MAXDIM+c].src.slab.str);
                  assert ((srcipos[c] - allinfo[n*SLAB_MAXDIM+c].src.slab.off) % allinfo[n*SLAB_MAXDIM+c].src.slab.str == 0);
                  bufipos[d] = ipos[d];
                  assert (bufipos[d] >= 0 and bufipos[d] < srcdetail[n*SLAB_MAXDIM+c].len);
                }
                size_t srcind = 0;
                size_t bufind = 0;
                for (int d=SLAB_MAXDIM-1; d>=0; --d) {
                  int const c = info[d].xpose;
                  srcind = srcind * info[d].src.local.alen + srcipos[d] - info[d].src.local.off;
                  bufind = bufind * srcdetail[n*SLAB_MAXDIM+c].len + bufipos[d];
                }
                assert (srcind < srclentot);
                assert (bufind < (size_t)srccount[n]);
                memcpy (srcdataptr + vartypesize * bufind,
                        srcptr + vartypesize * srcind,
                        vartypesize);
              }
            }
          }
          
        } // for vari
        CCTK_TimerStopI (timer_copy_in_general);
        
      }
    } // for n
    
    ifcheck {
      if (vartype == CCTK_VARIABLE_REAL) {
        CCTK_REAL const * restrict const srcdataptr =
          (CCTK_REAL const *)&srcdata.front();
        CCTK_REAL marker;
        memset (&marker, POISON_VALUE, sizeof marker);
        for (int i = 0; i < srcoffset[size]; ++i) {
          assert (memcmp(&srcdataptr[i], &marker, sizeof marker) != 0);
        }
      }
    }
    CCTK_TimerStopI (timer_copy_in);
    
    
    
    CCTK_TimerStartI (timer_xfer);
    ifdebug fflush (stdout);
    if (not HAVE_MPI or use_alltoallv) {
      MPI_Alltoallv
        (&srcdata.front(), &srccount.front(), &srcoffset.front(), vardatatype,
         &dstdata.front(), &dstcount.front(), &dstoffset.front(), vardatatype,
         comm);
    } else {
      vector<MPI_Request> requests;
      requests.reserve (2 * size);
      // Start receive
      for (int n = 0; n < size; ++n) {
        if (n != rank and dstcount[n] > 0) {
          MPI_Request req;
          MPI_Irecv
            (&dstdata[vartypesize * dstoffset[n]], dstcount[n], vardatatype,
             n, 0, comm, &req);
          requests.push_back (req);
        }
      }
      // Start send
      for (int n = 0; n < size; ++n) {
        if (n != rank and srccount[n] > 0) {
          MPI_Request req;
          MPI_Isend
            (&srcdata[vartypesize * srcoffset[n]], srccount[n], vardatatype,
             n, 0, comm, &req);
          requests.push_back (req);
        }
      }
      // Self communication
      {
        int const n = rank;
        assert (dstcount[n] == srccount[n]);
        memcpy (&dstdata[vartypesize * dstoffset[n]],
                &srcdata[vartypesize * srcoffset[n]],
                dstcount[n] * vartypesize);
      }
      // Wait
      MPI_Waitall (requests.size(), &requests.front(), MPI_STATUSES_IGNORE);
    }
    
    ifcheck {
      if (vartype == CCTK_VARIABLE_REAL) {
        for (int vari=0; vari<nvaris; ++vari) {
          CCTK_REAL const * restrict const dstdataptr =
            (CCTK_REAL const *)&dstdata.front();
          CCTK_REAL marker;
          memset (&marker, POISON_VALUE, sizeof marker);
          for (int i = 0; i < dstoffset[size]; ++i) {
            assert (memcmp(&dstdataptr[i], &marker, sizeof marker) != 0);
          }
        }
      }
    }
    CCTK_TimerStopI (timer_xfer);
    
    
    
    CCTK_TimerStartI (timer_copy_back);
    for (int n = 0; n < size; ++n) {
      check (SLAB_MAXDIM == 3);
      
      if (dstdetail[n*SLAB_MAXDIM  ].str==1 and
          dstdetail[n*SLAB_MAXDIM+1].str==1 and
          dstdetail[n*SLAB_MAXDIM+2].str==1 and
          vartype == CCTK_VARIABLE_REAL)
      {
        // Optimised version for stride 1 and CCTK_REAL
        copy_data_back<CCTK_REAL>
          (info, dstdetail, dstoffset, dstelems, dstdata, dstptrs, 
           n, varis, nvaris, info[0].flip, info[1].flip, info[2].flip);
      } else if (dstdetail[n*SLAB_MAXDIM  ].str==1 and
                 dstdetail[n*SLAB_MAXDIM+1].str==1 and
                 dstdetail[n*SLAB_MAXDIM+2].str==1 and
                 vartype == CCTK_VARIABLE_COMPLEX)
      {
        // Optimised version for stride 1 and CCTK_COMPLEX
        copy_data_back<CCTK_COMPLEX>
          (info, dstdetail, dstoffset, dstelems, dstdata, dstptrs, 
           n, varis, nvaris, info[0].flip, info[1].flip, info[2].flip);
      } else if (dstdetail[n*SLAB_MAXDIM  ].str==1 and
                 dstdetail[n*SLAB_MAXDIM+1].str==1 and
                 dstdetail[n*SLAB_MAXDIM+2].str==1 and
                 vartype == CCTK_VARIABLE_INT)
      {
        // Optimised version for stride 1 and CCTK_INT
        copy_data_back<CCTK_INT>
          (info, dstdetail, dstoffset, dstelems, dstdata, dstptrs, 
           n, varis, nvaris, info[0].flip, info[1].flip, info[2].flip);
      } else {
        // Generic, unoptimised version
        CCTK_TimerStartI (timer_copy_back_general);
        
        int const dstdetailleni = dstdetail[n*SLAB_MAXDIM+0].len;
        int const dstdetaillenj = dstdetail[n*SLAB_MAXDIM+1].len;
        int const dstdetaillenk = dstdetail[n*SLAB_MAXDIM+2].len;
        
        for (int vari=0; vari<nvaris; ++vari) {
          char * restrict const dstptr = (char * restrict)dstptrs[varis[vari]];
          char const * restrict const dstdataptr =
            &dstdata[vartypesize * (dstoffset[n] + vari * dstelems[n])];
          
#         pragma omp parallel for
          for (int k = 0; k < dstdetaillenk; ++k) {
            for (int j = 0; j < dstdetaillenj; ++j) {
              for (int i = 0; i < dstdetailleni; ++i) {
                int ipos[SLAB_MAXDIM];
                ipos[0] = i;
                ipos[1] = j;
                ipos[2] = k;
                int bufipos[SLAB_MAXDIM];
                int dstipos[SLAB_MAXDIM];
                for (int d=0; d<SLAB_MAXDIM; ++d) {
                  if (not info[d].flip) {
                    bufipos[d] = ipos[d];
                  } else {
                    bufipos[d] = dstdetail[n*SLAB_MAXDIM+d].len - 1 - ipos[d];
                  }
                  ifcheck assert (bufipos[d] >= 0 and bufipos[d] < dstdetail[n*SLAB_MAXDIM+d].len);
                  dstipos[d] = dstdetail[n*SLAB_MAXDIM+d].off + ipos[d] * info[d].dst.slab.str;
                  ifcheck assert (dstipos[d] >= info[d].dst.local.off and
                                  dstipos[d] < info[d].dst.local.off + info[d].dst.local.len);
                  ifcheck assert (dstipos[d] >= info[d].dst.slab.off and
                                  dstipos[d] <= info[d].dst.slab.off + (info[d].dst.slab.len - 1) * info[d].dst.slab.str);
                  ifcheck assert ((dstipos[d] - info[d].dst.slab.off) % info[d].dst.slab.str == 0);
                }
                size_t bufind = 0;
                size_t dstind = 0;
                for (int d=SLAB_MAXDIM-1; d>=0; --d) {
                  bufind = bufind * dstdetail[n*SLAB_MAXDIM+d].len + bufipos[d];
                  dstind = dstind * info[d].dst.local.alen + dstipos[d] - info[d].dst.local.off;
                }
                ifcheck assert (bufind < (size_t)dstcount[n]);
                ifcheck assert (dstind < dstlentot);
                memcpy (dstptr + vartypesize * dstind,
                        dstdataptr + vartypesize * bufind,
                        vartypesize);
              }
            }
          }
          
        } // for vari
        CCTK_TimerStopI (timer_copy_back_general);
        
      }
      
    } // for n
    CCTK_TimerStopI (timer_copy_back);
    
  } // for vartypei
  
  
  
  ifcheck {
    ifdebug fflush (stdout);
    MPI_Barrier (comm);
  }
  
  CCTK_TimerStopI (timer_apply);
  
  return 0;
}



extern "C"
int
Slab_MultiTransfer_Finalize
(cGH       const * restrict const cctkGH,
 slabsetup       * restrict const slabsetup)
{
  DECLARE_CCTK_PARAMETERS;
  
  // Check arguments
  check (cctkGH);
  check (slabsetup);
  
  delete slabsetup;
  
  return 0;
}



// Interface for transferring a variable in one go
extern "C"
int
Slab_MultiTransfer (cGH      const                  * restrict const cctkGH,
                    int                                        const dim,
                    xferinfo const                  * restrict const xferinfo,
                    int                                        const options,
                    int                                        const nvars,
                    int      const                  * restrict const srctypes,
                    void     const * restrict const * restrict const srcptrs,
                    int      const                  * restrict const dsttypes,
                    void           * restrict const * restrict const dstptrs)
{
  slabsetup * restrict const slabsetup =
    Slab_MultiTransfer_Init (cctkGH, dim, xferinfo, options);
  Slab_MultiTransfer_Apply (cctkGH, slabsetup,
                            nvars, srctypes, srcptrs, dsttypes, dstptrs);
  Slab_MultiTransfer_Finalize (cctkGH, slabsetup);
  return 0;
}



// Old interface for transferring a single variable
extern "C"
int
Slab_Transfer (cGH      const * restrict const cctkGH,
               int                       const dim,
               xferinfo const * restrict const xferinfo,
               int                       const options,
               int                       const srctype,
               void     const * restrict const srcptr,
               int                       const dsttype,
               void           * restrict const dstptr)
{
  int const nvars = 1;
  int const srctypes[] = { srctype };
  void const * restrict const srcptrs[] = { srcptr };
  int const dsttypes[] = { dsttype };
  void * restrict const dstptrs[] = { dstptr };
  return Slab_MultiTransfer (cctkGH, dim, xferinfo, options,
                             nvars, srctypes, srcptrs, dsttypes, dstptrs);
}



// Fortran wrapper
extern "C"
void CCTK_FCALL
CCTK_FNAME(Slab_Transfer) (int                * restrict const ierr,
                           cGH  const * const * restrict const cctkGH,
                           int  const         * restrict const dim,
                           int  const         * restrict const src_gsh,
                           int  const         * restrict const src_lbnd,
                           int  const         * restrict const src_lsh,
                           int  const         * restrict const src_ash,
                           int  const         * restrict const src_lbbox,
                           int  const         * restrict const src_ubbox,
                           int  const         * restrict const src_nghostzones,
                           int  const         * restrict const src_off,
                           int  const         * restrict const src_str,
                           int  const         * restrict const src_len,
                           int  const         * restrict const dst_gsh,
                           int  const         * restrict const dst_lbnd,
                           int  const         * restrict const dst_lsh,
                           int  const         * restrict const dst_ash,
                           int  const         * restrict const dst_lbbox,
                           int  const         * restrict const dst_ubbox,
                           int  const         * restrict const dst_nghostzones,
                           int  const         * restrict const dst_off,
                           int  const         * restrict const dst_str,
                           int  const         * restrict const dst_len,
                           int  const         * restrict const xpose,
                           int  const         * restrict const flip,
                           int  const         * restrict const options,
                           int  const         * restrict const srctype,
                           void const         * restrict const srcptr,
                           int  const         * restrict const dsttype,
                           void               * restrict const dstptr)
{
  vector<xferinfo> xferinfo (*dim);
  
  for (int d=0; d<*dim; ++d) {
    xferinfo[d].src.gsh         = src_gsh[d];
    xferinfo[d].src.lbnd        = src_lbnd[d];
    xferinfo[d].src.lsh         = src_lsh[d];
    xferinfo[d].src.ash         = src_ash[d];
    xferinfo[d].src.lbbox       = src_lbbox[d];
    xferinfo[d].src.ubbox       = src_ubbox[d];
    xferinfo[d].src.nghostzones = src_nghostzones[d];
    xferinfo[d].src.off         = src_off[d];
    xferinfo[d].src.str         = src_str[d];
    xferinfo[d].src.len         = src_len[d];
    
    xferinfo[d].dst.gsh         = dst_gsh[d];
    xferinfo[d].dst.lbnd        = dst_lbnd[d];
    xferinfo[d].dst.lsh         = dst_lsh[d];
    xferinfo[d].dst.ash         = dst_ash[d];
    xferinfo[d].dst.lbbox       = dst_lbbox[d];
    xferinfo[d].dst.ubbox       = dst_ubbox[d];
    xferinfo[d].dst.nghostzones = dst_nghostzones[d];
    xferinfo[d].dst.off         = dst_off[d];
    xferinfo[d].dst.str         = dst_str[d];
    xferinfo[d].dst.len         = dst_len[d];
    
    xferinfo[d].xpose           = xpose[d];
    xferinfo[d].flip            = flip[d];
  }
  
  *ierr = Slab_Transfer (*cctkGH, *dim, &xferinfo.front(), *options,
                         *srctype, srcptr, *dsttype, dstptr);
}
