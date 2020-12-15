/*@@
  @header    ioutil_Utils.h
  @date      Tue 19 Sep 2000
  @author    Thomas Radke
  @desc
             Function prototypes for setting up slice centers.
  @enddesc
  @version   $Header$
@@*/

#ifndef _IOUTIL_IOUTIL_UTILS_H_
#define _IOUTIL_IOUTIL_UTILS_H_ 1

#ifdef __cplusplus
extern "C" {
#endif

/* structure describing an I/O request (including hyperslab parameters) */
typedef struct {
  /* output frequency */
  CCTK_INT out_every;
  CCTK_REAL out_dt;

  /* index and timelevel of the variable */
  int vindex, timelevel;

  /* dimensionality of the variable and the hyperslab */
  CCTK_INT vdim, hdim;

  /* CCTK datatype for the hyperslab */
  int hdatatype;

  /* flag indicating whether an object to be written already exists
     (and remove it in that case) */
  int check_exist;

  /* flag indicating whether to include ghostzones in the hyperslab mapping */
  int with_ghostzones;

  /* flag indicating whether to output in chunked or unchunked format */
  int out_unchunked;

  /* compression level */
  CCTK_INT compression_level;

  /* pointer to allocated buffers */
  CCTK_INT *vectors;

  /* hyperslab mapping parameters */
  CCTK_INT *origin, *direction, *extent, *downsample;

  /* offset and sizes of hyperslab into the variable's dataspace */
  CCTK_INT *hoffset, *hsize, *hsize_chunk;

  /* bitmask for refinement levels to output */
  CCTK_INT refinement_levels;

  /* list of reductions to output */
  char *reductions;

} ioRequest;

/* Advertise that this is the new version of this API, which has
   'out_dt' arguments in the IOUtil_* functions below  */
#define IOUTIL_PARSER_HAS_OUT_DT 1

/* parse a given 'out_vars' parameter string */
void IOUtil_ParseVarsForOutput(const cGH *GH, const char *method_name,
                               const char *parameter_name,
                               int stop_on_parse_errors, const char *out_vars,
                               int out_every_default, CCTK_REAL out_dt_default,
                               ioRequest *request_list[]);

/* parse a given I/O parameter option string for the 'out_every' and
   'out_dt' options */
void IOUtil_ParseOutputFrequency(const char *method_name,
                                 const char *parameter_name,
                                 int stop_on_parse_errors, int vindex,
                                 const char *optstring, CCTK_INT *out_every,
                                 CCTK_REAL *out_dt);

/* return the default I/O request description structure for a variable */
ioRequest *IOUtil_DefaultIORequest(const cGH *GH, int vindex,
                                   int out_every_default,
                                   CCTK_REAL out_dt_default);

/* free an I/O request description */
void IOUtil_FreeIORequest(ioRequest **request);

/* set the slice center for 1D lines */
int IOUtil_1DLines(const cGH *GH, int num_dims, int *const *const origin_index,
                   CCTK_REAL *const *const origin_phys,
                   int *const *const slice_center);

/* set the slice center for 2D planes */
int IOUtil_2DPlanes(const cGH *GH, int num_dims, const int *origin_index,
                    const CCTK_REAL *origin_phys, int *slice_center);

/* create an output directory on all I/O processors */
int IOUtil_CreateDirectory(const cGH *GH, const char *dirname,
                           int multiple_io_procs, int ioproc);

#ifdef __cplusplus
}
#endif

#endif /* _IOUTIL_IOUTIL_UTILS_H_ */
