/*@@
  @file      hdf5_merge.c
  @date      Thu 10 Jan 2002
  @author    Thomas Radke
  @desc
             This utility program takes a list of Cactus HDF5 datafiles,
             merges them at the group hierarchy level and dumps the resulting
             tree to a new HDF5 file.
  @enddesc
  @version   $Id: hdf5_merge.c,v 1.3 2009/09/29 14:38:14 schnetter Exp $
@@*/

#include "cctk.h"

#define H5_USE_16_API 1
#include <hdf5.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* the rcs ID and its dummy function to use it */
static const char *rcsid =
    "$Header: /cactusdevcvs/CactusExternal/HDF5/src/util/hdf5_merge.c,v 1.3 "
    "2009/09/29 14:38:14 schnetter Exp $";
CCTK_FILEVERSION(CactusExternal_HDF5_util_hdf5_merge_c)

/*****************************************************************************/
/*                           macro definitions                               */
/*****************************************************************************/
/* macro to do an HDF5 call, check its return code, and print a warning
   in case of an error */
#define CHECK_ERROR(hdf5_call)                                                 \
  do {                                                                         \
    int _error_code = hdf5_call;                                               \
                                                                               \
    if (_error_code < 0) {                                                     \
      fprintf(stderr, "WARNING: line %d: HDF5 call '%s' returned "             \
                      "error code %d\n",                                       \
              __LINE__, #hdf5_call, _error_code);                              \
      nerrors++;                                                               \
    }                                                                          \
  } while (0)

/*****************************************************************************/
/*                           global variables                                */
/*****************************************************************************/
/* NOTE: although it isn't good programming practice
         we make these variables global for convenience
         since they are accessed from recursively or
         indirectly called routines which only get passed
         a single user-supplied argument */
static char *pathname = NULL;    /* pathname of the current object */
static unsigned int nerrors = 0; /* global error counter */

/*****************************************************************************/
/*                           local function prototypes                       */
/*****************************************************************************/
static herr_t CopyObject(hid_t copy_from, const char *objectname, void *arg);
static herr_t CopyAttribute(hid_t src, const char *attr_name, void *arg);

/*@@
  @routine    main
  @date       Sat 24 Feb 2001
  @author     Thomas Radke
  @desc
              Main routine of the HDF5 file merger
  @enddesc

  @calls      CopyObject

  @var        argc
  @vdesc      number of command line arguments
  @vtype      int
  @vio        in
  @endvar
  @var        argv
  @vdesc      command line arguments
  @vtype      char *[]
  @vio        in
  @endvar

  @returntype int
  @returndesc
              0 for success, negative return values indicate an error
  @endreturndesc
@@*/
int main(int argc, char *argv[]) {
  int i;
  hid_t infile, outfile;

  /* give some help if called with incorrect number of parameters */
  if (argc < 3) {
    fprintf(stderr, "Usage: %s <infile1> [<infile2> ...] <outfile>\n", argv[0]);
    fprintf(stderr, "   eg, %s alp.time*.h5 alp.h5"
                    "\n\n",
            argv[0]);
    return (0);
  }

  H5E_BEGIN_TRY {
    /* open the input files as early test that this works */
    for (i = 0; i < argc - 2; i++) {
      infile = H5Fopen(argv[i + 1], H5F_ACC_RDONLY, H5P_DEFAULT);
      if (infile < 0) {
        fprintf(stderr, "ERROR: Cannot open HDF5 input file '%s' !\n\n",
                argv[i + 1]);
        return (-1);
      }
      /* close them again, or we might have too many open files and fail */
      H5Fclose(infile);
    }

    /* try to open an existing outfile file in append mode,
       if this fails create it as a new file */
    outfile = H5Fopen(argv[argc - 1], H5F_ACC_RDWR, H5P_DEFAULT);
    if (outfile < 0) {
      outfile =
          H5Fcreate(argv[argc - 1], H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    }
    if (outfile < 0) {
      fprintf(stderr, "ERROR: Cannot open HDF5 output file '%s' !\n\n",
              argv[argc - 1]);
      return (-1);
    }
  }
  H5E_END_TRY

  printf("\n  -------------------------\n"
         "  Cactus 4 HDF5 File Merger\n"
         "  -------------------------\n");

  /* do the copying by iterating over all objects */
  for (i = 0; i < argc - 2; i++) {
    printf("\n  Merging objects from input file '%s' into output file '%s'\n",
           argv[i + 1], argv[argc - 1]);
    pathname = "";
    H5E_BEGIN_TRY {
      infile = H5Fopen(argv[i + 1], H5F_ACC_RDONLY, H5P_DEFAULT);
      if (infile < 0) {
        fprintf(stderr, "ERROR: Cannot open HDF5 input file '%s' !\n\n",
                argv[i + 1]);
        return (-1);
      }
    }
    H5E_END_TRY
    CHECK_ERROR(H5Giterate(infile, "/", NULL, CopyObject, &outfile));
    CHECK_ERROR(H5Fclose(infile));
  }

  CHECK_ERROR(H5Fclose(outfile));

  /* report status */
  if (nerrors == 0) {
    printf("\n\n   *** All input files successfully merged. ***\n\n");
  } else {
    fprintf(stderr, "\n\n   *** WARNING: %u errors occured during "
                    "file merging. ***\n\n",
            nerrors);
  }

  return (0);
}

/*****************************************************************************/
/*                           local routines                                  */
/*****************************************************************************/
/*@@
  @routine    CopyObject
  @date       Sat 24 Feb 2001
  @author     Thomas Radke
  @desc
              Iterator recursively called by H5Giterate() for every object
              in the input file
              It copies the current object to the output file if it didn't
              already exist there.
  @enddesc

  @calls      CopyAttribute

  @var        from
  @vdesc      identifier for the group the current object belongs to
  @vtype      hid_t
  @vio        in
  @endvar
  @var        objectname
  @vdesc      name of the current object
  @vtype      const char *
  @vio        in
  @endvar
  @var        _to
  @vdesc      user-supplied argument indicating the output object identifier
  @vtype      hid_t
  @vio        in
  @endvar

  @returntype int
  @returndesc
              0 - continue the iteration for following group objects
              1 - short-curcuit, no further iteration of this group
  @endreturndesc
@@*/
static herr_t CopyObject(hid_t from, const char *objectname, void *_to) {
  int i;
  hid_t to, datatype, dataspace;
  H5G_stat_t objectinfo;
  char *current_pathname;
  size_t objectsize;
  void *data;

  /* get the output object identifier */
  to = *(hid_t *)_to;

  /* check whether an object by that name already exists */
  H5E_BEGIN_TRY { i = H5Gget_objinfo(to, objectname, 0, &objectinfo) >= 0; }
  H5E_END_TRY
  if (i) {
    printf("   object '%s/%s' will not be copied (already exists)\n", pathname,
           objectname);
    return (0);
  }

  /* build the full pathname for the current to object to process */
  current_pathname = pathname;
  pathname = (char *)malloc(strlen(current_pathname) + strlen(objectname) + 2);
  sprintf(pathname, "%s/%s", current_pathname, objectname);

  /* check the type of the current object */
  CHECK_ERROR(H5Gget_objinfo(from, objectname, 0, &objectinfo));
  if (objectinfo.type == H5G_GROUP) {
    printf("   copying group '%s'\n", pathname);

    CHECK_ERROR(from = H5Gopen(from, objectname));
    CHECK_ERROR(to = H5Gcreate(to, objectname, 0));
    /* iterate over all objects in the (first) input file */
    CHECK_ERROR(H5Giterate(from, ".", NULL, CopyObject, &to));
    CHECK_ERROR(H5Aiterate(from, NULL, CopyAttribute, &to));
    CHECK_ERROR(H5Gclose(to));
    CHECK_ERROR(H5Gclose(from));
  } else if (objectinfo.type == H5G_DATASET) {
    printf("   copying dataset '%s'\n", pathname);

    CHECK_ERROR(from = H5Dopen(from, objectname));
    CHECK_ERROR(datatype = H5Dget_type(from));
    CHECK_ERROR(dataspace = H5Dget_space(from));
    CHECK_ERROR(
        to = H5Dcreate(to, objectname, datatype, dataspace, H5P_DEFAULT));
    objectsize = H5Sget_select_npoints(dataspace) * H5Tget_size(datatype);
    if (objectsize > 0) {
      data = malloc(objectsize);
      if (data == NULL) {
        fprintf(stderr, "failled to allocate %d bytes of memory, giving up\n",
                (int)objectsize);
        exit(-1);
      }
      CHECK_ERROR(H5Dread(from, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data));
      CHECK_ERROR(H5Dwrite(to, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data));
      free(data);
    }
    CHECK_ERROR(H5Aiterate(from, NULL, CopyAttribute, &to));
    CHECK_ERROR(H5Dclose(to));
    CHECK_ERROR(H5Dclose(from));
    CHECK_ERROR(H5Sclose(dataspace));
    CHECK_ERROR(H5Tclose(datatype));
  } else {
    fprintf(stderr, "WARNING: Found object '%s' which is not neither a "
                    "group nor a dataset ! Object will not be copied.\n",
            pathname);
    nerrors++;
  }

  /* reset the pathname */
  free(pathname);
  pathname = current_pathname;

  return (0);
}

/*@@
  @routine    CopyAttribute
  @date       Sat 24 Feb 2001
  @author     Thomas Radke
  @desc
              Iterator recursively called by H5Aiterate() for every attribute
              of an object (dataset or group)
  @enddesc

  @var        from
  @vdesc      identifier for the group or dataset to read the attribute from
  @vtype      hid_t
  @vio        in
  @endvar
  @var        attrname
  @vdesc      name of the current attribute
  @vtype      const char *
  @vio        in
  @endvar
  @var        _to
  @vdesc      user-supplied argument indicating the group or dataset
              to copy the attribute to
  @vtype      hid_t
  @vio        in
  @endvar

  @returntype int
  @returndesc
              0 - continue the iteration for following attributes
  @endreturndesc
@@*/
static herr_t CopyAttribute(hid_t from, const char *attrname, void *_to) {
  hid_t attr, datatype, dataspace, to;
  size_t attrsize;
  void *value;

  /* get the target group/dataset identifier */
  to = *(hid_t *)_to;

  /* open the attribute given by its name, get type, dataspace, and value
     and just copy it */
  CHECK_ERROR(attr = H5Aopen_name(from, attrname));
  CHECK_ERROR(datatype = H5Aget_type(attr));
  CHECK_ERROR(dataspace = H5Aget_space(attr));
  attrsize = H5Tget_size(datatype);
  if (H5Sis_simple(dataspace) > 0) {
    attrsize *= H5Sget_simple_extent_npoints(dataspace);
  }
  if (attrsize > 0) {
    value = malloc(attrsize);
    CHECK_ERROR(H5Aread(attr, datatype, value));
    CHECK_ERROR(H5Aclose(attr));
    CHECK_ERROR(attr =
                    H5Acreate(to, attrname, datatype, dataspace, H5P_DEFAULT));
    CHECK_ERROR(H5Awrite(attr, datatype, value));
    free(value);
  }
  CHECK_ERROR(H5Aclose(attr));
  CHECK_ERROR(H5Sclose(dataspace));
  CHECK_ERROR(H5Tclose(datatype));

  return (0);
}
