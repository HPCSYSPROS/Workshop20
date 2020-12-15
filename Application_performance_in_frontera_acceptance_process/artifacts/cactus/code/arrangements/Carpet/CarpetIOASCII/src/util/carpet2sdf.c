
/*******************************************************************
 * Program: carpet2sdf
 * Description: Converts from Carpet ASCII file format to SDF format
 * Author:  Scott H. Hawley
 * Date:    June 17, 2002
 *
 * Similar to carpet2sdf, except that instead of sending output to
 * sdtout, it sends it to <infile>.sdf
 *******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <bbhutil.h>

/*******************************************************************
 * Function: read_next_set
 * Supposed to read one time step of data (for one level of refinement)
 *******************************************************************/
int read_next_set(FILE *infile, int rlev, int *numelems, double *time,
                  double **coord, double **data, int clip_data,
                  double clip_val) {
  char in_line[200];
  int it, tl, rl, c, ml, ix, iy, iz;
  double xval, yval, zval;
  double coordval, dataval;
  int hit_first_non_blank = 0;
  int stop_reading = 0;
  int retval = 0;
  int rc = 0;

  *numelems = 0;
  sprintf(in_line, " ");

  /* reads up to first blank line after non-blank lines, or up to
   * eof */
  while ((fgets(in_line, 200, infile) != NULL) && (!stop_reading)) {

    if (in_line[0] == '#') {
      hit_first_non_blank = 1;

    } else if ((strncmp("\n", in_line, 2) == 0) || (in_line[0] == ' ')) {
      if (hit_first_non_blank) {
        stop_reading = 1;
      }

    } else { /* Assume that we're dealing with a line of numbers */
      hit_first_non_blank = 1;
      retval = sscanf(in_line, "%d %d %d %d %d %d %d %d %lg %lg %lg %lg %lg",
                      &it, &tl, &rl, &c, &ml, &ix, &iy, &iz, time, &xval, &yval,
                      &zval, &dataval);
      coordval = xval;

      /* Only add input to arrays if it's for the right ref. level */
      if ((retval == 13) && (rl == rlev)) {
        (*numelems)++;

        if (coord == NULL) {
          *coord = (double *)malloc(sizeof(double) * (*numelems));
          *data = (double *)malloc(sizeof(double) * (*numelems));
        } else {
          *coord = (double *)realloc(*coord, sizeof(double) * (*numelems));
          *data = (double *)realloc(*data, sizeof(double) * (*numelems));
        }
        if (coord == NULL) {
          fprintf(stderr, "Error: malloc/realloc returned NULL\n");
          exit(1);
        }

        (*coord)[(*numelems) - 1] = coordval;

        /* If desired, "clip" data above a certain value */
        if (clip_data && (dataval > clip_val)) {
          (*data)[(*numelems) - 1] = clip_val;
        } else {
          (*data)[(*numelems) - 1] = dataval;
        }
      }
    }
  }
  rc = feof(infile);

  return rc;
}

/********************************************************************
 * Main part of program
 ********************************************************************/
int main(int argc, char **argv) {
  char *infilename;
  FILE *infile;

  char func_name[200];
  double time;
  int shape[1];
  double bounds[2];
  int rank = 1;
  int numelems_in_set = 0;
  int rlev = 0;
  double *coord, *data;
  int clip_data = 0;
  double clip_val;

  /*
   * Parse command-line arguments
   */
  if (argc <= 1) {
    fprintf(stderr, "usage: carpet2sdf <infile> [ref_lev] [clip_data_at]\n");
    exit(1);
  }
  if (argc >= 3) {
    sscanf(argv[2], "%d", &rlev);
  }
  if (argc >= 4) {
    sscanf(argv[3], "%lf", &clip_val);
    clip_data = 1;
  }

  /*
   * Open the input file
   */
  infilename = argv[1];
  infile = fopen(infilename, "r");
  if (infile == NULL) {
    fprintf(stderr, "Error opening file '%s'.\n", infilename);
    exit(1);
  }

  time = 0;
  coord = NULL;
  data = NULL;
  sprintf(func_name, "%s_%d", infilename, rlev);

  /*
   * Main loop for reading from input file and writing to output file
   */
  while (read_next_set(infile, rlev, &numelems_in_set, &time, &coord, &data,
                       clip_data, clip_val) == 0) {
    if (numelems_in_set > 0) {

      shape[0] = numelems_in_set;
      bounds[0] = coord[0];
      bounds[1] = coord[numelems_in_set - 1];

      gft_out_bbox(func_name, time, shape, rank, bounds, data);

      free(coord);
      free(data);
      coord = NULL;
      data = NULL;
    }
  }

  return 0;
}
