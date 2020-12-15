/*@@
  @header    ioGH.h
  @date      Tue 9th Jan 1999
  @author    Gabrielle Allen
  @desc
             The extensions to the GH structure from IOUtil.
  @enddesc
  @version   $Header$
@@*/

#ifndef _IOUTIL_IOGH_H_
#define _IOUTIL_IOGH_H_ 1

#ifdef __cplusplus
extern "C" {
#endif

/* advertise that this is the new version of this API, which has 'alias'
 * members in the grid extension */
#define IOUTIL_IOGH_HAS_ALIAS 1

/* IOUtil's GH extension structure */
typedef struct {
  /* for full-dimensional parallel output */
  int ioproc;       /* the I/O processor each proc belongs to */
  int nioprocs;     /* total number of I/O processors */
  int ioproc_every; /* output by every N'th processor */
  int unchunked;    /* if true generate unchunked output file */
  int *downsample;  /* downsampling parameters array of size cctk_maxdim */
  int out_single;   /* if true output 3D data in single precision */

  /* for recovery */
  int recovered; /* flag indicating restart after successful recovery */
  int stop_on_parse_errors; /* stop on I/O parameter parsing errors ? */

  /* for data file reader */
  CCTK_INT *do_inVars; /* flags indicating to read in variable i with
                          iteration number do_inVars[i] (or -1 to read
                          the last iteration */
  const char **alias;  /* name under which a variable appears in data files.
                          If NULL, use CCTK_FullName() */
} ioGH;

#ifdef __cplusplus
}
#endif

#endif /* _IOUTIL_IOGH_H_ */
