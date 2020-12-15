/* test_molecule_posn -- test driver for AEILocalInterp_molecule_posn() */
/* $Header$ */

/*
 * This program is a test driver for  AEILocalInterp_molecule_posn() .
 *
 * Usage:
 *	test_molecule_posn			# run a preset set of tests
 *	test_molecule_posn molecule_size		\
 *	                   boundary_off_centering_tolerance_min	\
 *	                   boundary_off_centering_tolerance_max	\
 *	                   boundary_extrapolation_tolerance_min	\
 *	                   boundary_extrapolation_tolerance_max	\
 *	                   boundary_off_centering_tolerance_max	\
 *	                   x			# do a single test as specified
 */

const char* help_msg =
"usage:\n"
"# run a preset series of tests:\n"
"   test_molecule_posn\n"
"# run a single test as specified:\n"
"   test_molecule_posn molecule_size \\\n"
"                      boundary_off_centering_tolerance_min \\\n"
"                      boundary_off_centering_tolerance_max \\\n"
"                      boundary_extrapolation_tolerance_min \\\n"
"                      boundary_extrapolation_tolerance_max \\\n"
"                      x\n"
;

#include <math.h>
#include <string.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef AEILOCALINTERP_STANDALONE_TEST
  #include "cctk.h"
#endif

#include "InterpLocalUniform.h"

#define fuzzy_EQ(x,y)	(fabs(x-y) <= 1.0e-10)

/* hard-wired arguments for  AEILocalInterp_molecule_posn() */
const fp grid_x0 = 3.1;
const fp grid_dx = 0.1;
const int grid_i_min = 42;	/* x_min =  7.3 */
const int grid_i_max = 105;	/* x_max = 13.6 */

/******************************************************************************/

/*
 * prototypes for functions local to this file
*/
static
  void run_interactive_test(int molecule_size,
			    fp boundary_off_centering_tolerance_min,
			    fp boundary_off_centering_tolerance_max,
			    fp boundary_extrapolation_tolerance_min,
			    fp boundary_extrapolation_tolerance_max,
			    fp x);
static
  int run_batch_tests(void);

/******************************************************************************/

/*
 * test data for batch tests
 */

/* arguments/results structure for AEILocalInterp_molecule_posn() calls */
struct	args_results
	{
	/* args */
	int molecule_size;
	fp boundary_off_centering_tolerance_min;
	fp boundary_off_centering_tolerance_max;
	fp boundary_extrapolation_tolerance_min;
	fp boundary_extrapolation_tolerance_max;
	fp x;
	/* results */
	int status;
	int i_center;
	fp x_rel;
	};

/* test data */
#define OK		0
#define X_LT_MIN	MOLECULE_POSN_ERROR_X_LT_MIN
#define X_GT_MAX	MOLECULE_POSN_ERROR_X_GT_MAX
static const struct args_results test_data[] =
  {
    /*** molecule size 2 ***/
    /*                                        status           */
    /*     off_centering     extrapolation        i_center     */
    /*      min     max     min     max    x             x_rel */
    { 2,  999.0,  999.0,    1.0,  999.0,  7.19, X_LT_MIN,   0,  0.0 },
    { 2,  999.0,  999.0,    1.0,  999.0,  7.21,       OK,  42, -0.9 },
    { 2,  999.0,  999.0,    1.0,  999.0,  7.24,       OK,  42, -0.6 },
    { 2,  999.0,  999.0,    1.0,  999.0,  7.26,       OK,  42, -0.4 },
    { 2,    0.0,    0.0,    1.0,  999.0,  7.26, X_LT_MIN,   0,  0.0 },
    { 2, 1.0e-6,  999.0,    0.2,  999.0,  7.29, X_LT_MIN,   0,  0.0 },
    { 2,    0.2,  999.0, 1.0e-6,  999.0,  7.29, X_LT_MIN,   0,  0.0 },
    { 2,    0.2,  999.0,    0.2,  999.0,  7.29,       OK,  42, -0.1 },
    { 2,    0.0,    0.0, 1.0e-6,  999.0,  7.31,       OK,  42, +0.1 },
    { 2,    0.0,    0.0, 1.0e-6,  999.0,  7.34,       OK,  42, +0.4 },
    { 2,    0.0,    0.0, 1.0e-6,  999.0,  7.36,       OK,  42, +0.6 },
    { 2,    0.0,    0.0, 1.0e-6,  999.0,  7.39,       OK,  42, +0.9 },
    { 2,    0.0,    0.0, 1.0e-6,  999.0,  7.41,       OK,  43, +0.1 },
    { 2,    0.0,    0.0, 1.0e-6, 1.0e-6,  9.81,       OK,  67, +0.1 },
    { 2,    0.0,    0.0, 1.0e-6, 1.0e-6,  9.84,       OK,  67, +0.4 },
    { 2,    0.0,    0.0, 1.0e-6, 1.0e-6,  9.85,       OK,  67, +0.5 },
    { 2,    0.0,    0.0, 1.0e-6, 1.0e-6,  9.86,       OK,  67, +0.6 },
    { 2,    0.0,    0.0, 1.0e-6, 1.0e-6,  9.89,       OK,  67, +0.9 },
    { 2,    0.0,    0.0, 1.0e-6, 1.0e-6, 13.45,       OK, 103, +0.5 },
    { 2,    0.0,    0.0, 1.0e-6, 1.0e-6, 13.51,       OK, 104, +0.1 },
    { 2,    0.0,    0.0, 1.0e-6, 1.0e-6, 13.59,       OK, 104, +0.9 },
    { 2,    0.0,    0.2, 1.0e-6, 1.0e-6, 13.61, X_GT_MAX,   0,  0.0 },
    { 2,    0.0, 1.0e-6, 1.0e-6, 1.0e-6, 13.61, X_GT_MAX,   0,  0.0 },
    { 2,    0.0,    0.2,    0.0,    0.2, 13.61,       OK, 104, +1.1 },

    /*** molecule size 3 ***/
    /*                                        status           */
    /*     off_centering     extrapolation        i_center     */
    /*      min     max     min     max    x             x_rel */
    { 3,    0.5,  999.0,  999.0,  999.0,  7.29, X_LT_MIN,   0,  0.0 },
    { 3,  999.0,  999.0,    0.0,  999.0,  7.29, X_LT_MIN,   0,  0.0 },
    { 3,    0.7,  999.0,    0.2,  999.0,  7.29,       OK,  43, -1.1 },
    { 3,    0.3,  999.0,    0.0,  999.0,  7.31, X_LT_MIN,   0,  0.0 },
    { 3,    0.5,  999.0,    0.0,  999.0,  7.31,       OK,  43, -0.9 },
    { 3,    0.0,    0.0,    0.0,    0.0,  7.36,       OK,  43, -0.4 },
    { 3,    0.0,    0.0,    0.0,    0.0,  7.44,       OK,  43, +0.4 },
    { 3,    0.0,    0.0,    0.0,    0.0,  7.46,       OK,  44, -0.4 },
    { 3,    0.0,    0.0,    0.0,    0.0,  9.81,       OK,  67, +0.1 },
    { 3,    0.0,    0.0,    0.0,    0.0,  9.84,       OK,  67, +0.4 },
    { 3,    0.0,    0.0,    0.0,    0.0,  9.86,       OK,  68, -0.4 },
    { 3,    0.0,    0.0,    0.0,    0.0,  9.89,       OK,  68, -0.1 },
    { 3,    0.0,    0.0,    0.0,    0.0, 13.44,       OK, 103, +0.4 },
    { 3,    0.0,    0.0,    0.0,    0.0, 13.46,       OK, 104, -0.4 },
    { 3,    0.0,    0.0,    0.0,    0.0, 13.54,       OK, 104, +0.4 },
    { 3,    0.0,    0.2,    0.0,    0.0, 13.56,       OK, 104, +0.6 },
    { 3,  999.0,    0.0,  999.0,  999.0, 13.56, X_GT_MAX,   0,  0.0 },
    { 3,    0.0,    0.5,    0.0,    0.0, 13.59,       OK, 104, +0.9 },
    { 3,  999.0,    0.3,  999.0,  999.0, 13.59, X_GT_MAX,   0,  0.0 },
    { 3,  999.0,    0.5,  999.0,  999.0, 13.61, X_GT_MAX,   0,  0.0 },
    { 3,  999.0,  999.0,  999.0,    0.0, 13.61, X_GT_MAX,   0,  0.0 },
    { 3,  999.0,    0.7,  999.0,    0.2, 13.61,       OK, 104, +1.1 },
    { 3,  999.0,  999.0,  999.0,    0.2, 13.63, X_GT_MAX,   0,  0.0 },
    { 3,  999.0,    0.7,  999.0,  999.0, 13.63, X_GT_MAX,   0,  0.0 },
    { 3,  999.0,    0.9,  999.0,    0.4, 13.63,       OK, 104, +1.3 },

    /*** molecule size 4 ***/
    /*                                        status           */
    /*     off_centering     extrapolation        i_center     */
    /*      min     max     min     max    x             x_rel */
    { 4,    1.0,  999.0,  999.0,  999.0,  7.29, X_LT_MIN,   0,  0.0 },
    { 4,  999.0,  999.0,    0.0,  999.0,  7.29, X_LT_MIN,   0,  0.0 },
    { 4,    1.2,  999.0,    0.2,  999.0,  7.29,       OK,  43, -1.1 },
    { 4,    0.5,  999.0,    0.0,  999.0,  7.31, X_LT_MIN,   0,  0.0 },
    { 4,    1.0,  999.0,    0.0,  999.0,  7.31,       OK,  43, -0.9 },
    { 4,    0.0,  999.0,    0.0,  999.0,  7.39, X_LT_MIN,   0,  0.0 },
    { 4,    0.2,  999.0,    0.0,  999.0,  7.39,       OK,  43, -0.1 },
    { 4,    0.0,    0.0,    0.0,  999.0,  7.41,       OK,  43, +0.1 },
    { 4,    0.0,    0.0,    0.0,  999.0,  7.49,       OK,  43, +0.9 },
    { 4,    0.0,    0.0,    0.0,  999.0,  7.51,       OK,  44, +0.1 },
    { 4,    0.0,    0.0,    0.0,    0.0,  9.81,       OK,  67, +0.1 },
    { 4,    0.0,    0.0,    0.0,    0.0,  9.84,       OK,  67, +0.4 },
    { 4,    0.0,    0.0,    0.0,    0.0,  9.85,       OK,  67, +0.5 },
    { 4,    0.0,    0.0,    0.0,    0.0,  9.86,       OK,  67, +0.6 },
    { 4,    0.0,    0.0,    0.0,    0.0,  9.89,       OK,  67, +0.9 },
    { 4,    0.0,    0.0,    0.0,    0.0, 13.44,       OK, 103, +0.4 },
    { 4,    0.0,    0.0,    0.0,    0.0, 13.49,       OK, 103, +0.9 },
    { 4,    0.0,    0.2,    0.0,    0.0, 13.51,       OK, 103, +1.1 },
    { 4,  999.0,    0.0,  999.0,  999.0, 13.51, X_GT_MAX,   0,  0.0 },
    { 4,  999.0,    0.8,  999.0,    0.0, 13.59, X_GT_MAX,   0,  0.0 },
    { 4,    0.0,    1.0,    0.0,    0.0, 13.59,       OK, 103, +1.9 },
    { 4,  999.0,  999.0,  999.0,    0.0, 13.61, X_GT_MAX,   0,  0.0 },
    { 4,    0.0,    1.0,  999.0,  999.0, 13.61, X_GT_MAX,   0,  0.0 },
    { 4,    0.0,    1.2,    0.0,    0.2, 13.61,       OK, 103, +2.1 },

    /*** molecule size 5 ***/
    /*                                        status           */
    /*     off_centering     extrapolation        i_center     */
    /*      min     max     min     max    x             x_rel */
    { 5,    1.5,  999.0,  999.0,  999.0,  7.29, X_LT_MIN,   0,  0.0 },
    { 5,  999.0,  999.0,    0.0,  999.0,  7.29, X_LT_MIN,   0,  0.0 },
    { 5,    1.7,  999.0,    0.2,  999.0,  7.29,       OK,  44, -2.1 },
    { 5,    0.9,  999.0,    0.0,  999.0,  7.35, X_LT_MIN,   0,  0.0 },
    { 5,    1.1,  999.0,    0.0,  999.0,  7.35,       OK,  44, -1.5 },
    { 5,    0.0,  999.0,  999.0,  999.0,  7.44, X_LT_MIN,   0,  0.0 },
    { 5,    0.2,  999.0,    0.0,  999.0,  7.44,       OK,  44, -0.6 },
    { 5,    0.0,    0.0,    0.0,    0.0,  7.46,       OK,  44, -0.4 },
    { 5,    0.0,    0.0,    0.0,    0.0,  9.81,       OK,  67, +0.1 },
    { 5,    0.0,    0.0,    0.0,    0.0,  9.84,       OK,  67, +0.4 },
    { 5,    0.0,    0.0,    0.0,    0.0,  9.86,       OK,  68, -0.4 },
    { 5,    0.0,    0.0,    0.0,    0.0,  9.89,       OK,  68, -0.1 },
    { 5,    0.0,    0.0,    0.0,    0.0, 13.44,       OK, 103, +0.4 },
    { 5,  999.0,    0.0,  999.0,  999.0, 13.46, X_GT_MAX,   0,  0.0 },
    { 5,    0.0,    0.2,    0.0,    0.0, 13.46,       OK, 103, +0.6 },
    { 5,    0.0,    1.1,    0.0,    0.0, 13.55,       OK, 103, +1.5 },
    { 5,  999.0,    0.9,  999.0,  999.0, 13.55, X_GT_MAX,   0,  0.0 },
    { 5,  999.0,  999.0,  999.0,    0.0, 13.61, X_GT_MAX,   0,  0.0 },
    { 5,  999.0,    1.5,  999.0,  999.0, 13.61, X_GT_MAX,   0,  0.0 },
    { 5,  999.0,    1.7,  999.0,    0.2, 13.61,       OK, 103, +2.1 },
  };
#define N_TESTS ((int) (sizeof(test_data)/sizeof(test_data[0])))


/******************************************************************************/

int main(int argc, const char *const argv[])
{
bool N_fail;
int molecule_size;
fp boundary_off_centering_tolerance_min;
fp boundary_off_centering_tolerance_max;
fp boundary_extrapolation_tolerance_min;
fp boundary_extrapolation_tolerance_max;
double x;

switch	(argc)
	{
case 1:
	/* run batch tests */
	N_fail = run_batch_tests();
	if (N_fail == 0)
	   then {
		printf("*** all %d tests successful ***\n", N_TESTS);
		return 0;
		}
	   else {
		printf("*** %d/%d test(s) failed ***\n", N_fail, N_TESTS);
		return 1;
		}

case 7:
	if (    (sscanf(argv[1], "%d",  &molecule_size) == 1)
	     && (sscanf(argv[2], "%lf", &boundary_off_centering_tolerance_min) == 1)
	     && (sscanf(argv[3], "%lf", &boundary_off_centering_tolerance_max) == 1)
	     && (sscanf(argv[4], "%lf", &boundary_extrapolation_tolerance_min) == 1)
	     && (sscanf(argv[5], "%lf", &boundary_extrapolation_tolerance_max) == 1)
	     && (sscanf(argv[6], "%lf", &x) == 1)    )
	   then {
		run_interactive_test(molecule_size,
				     boundary_off_centering_tolerance_min,
				     boundary_off_centering_tolerance_max,
				     boundary_extrapolation_tolerance_min,
				     boundary_extrapolation_tolerance_max,
				     x);
		return 0;					/*** NORMAL EXIT ***/
		}
	/* error ==> fall through */

default:
	fprintf(stderr, help_msg);
	return 1;
	}
}

/******************************************************************************/

/*
 * This function runs a single test as specified.
 */
static
  void run_interactive_test(int molecule_size,
			    fp boundary_off_centering_tolerance_min,
			    fp boundary_off_centering_tolerance_max,
			    fp boundary_extrapolation_tolerance_min,
			    fp boundary_extrapolation_tolerance_max,
			    fp x)
{
int i_center;
fp x_rel;

printf("testing with molecule_size=%d\n", molecule_size);
printf("             boundary_off_centering_tolerance_[min,max]=[%g,%g]\n",
       (double) boundary_off_centering_tolerance_min,
       (double) boundary_off_centering_tolerance_max);
printf("             boundary_extrapolation_tolerance_[min,max]=[%g,%g]\n",
       (double) boundary_extrapolation_tolerance_min,
       (double) boundary_extrapolation_tolerance_max);
printf("             x=%g\n", (double) x);

  {
const int debug = 0;
const int status
   = AEILocalInterp_molecule_posn(grid_x0, grid_dx,
				  grid_i_min, grid_i_max,
				  molecule_size,
				  boundary_off_centering_tolerance_min,
				  boundary_off_centering_tolerance_max,
				  boundary_extrapolation_tolerance_min,
				  boundary_extrapolation_tolerance_max,
				  x,
				  debug,
				  &i_center, &x_rel);

if (status < 0)
   then printf("status=%d\n", status);
   else printf("i_center=%d x_rel=%g\n", i_center, x_rel);
  }
}

/******************************************************************************/

/*
 * This function run the preset set of tests specified by the
 *  test_data  array
 *
 * Results:
 * This function returns the number of test failures.
 */
static
  int run_batch_tests(void)
{
int i;
int failure_count = 0;

	for (i = 0 ; i < N_TESTS ; ++i)
	{
	const struct args_results *p = & test_data[i];
	int i_center = 999;
	fp  x_rel = 999.999;
	const int debug = 0;
	const int status = AEILocalInterp_molecule_posn
				(grid_x0, grid_dx,
				 grid_i_min, grid_i_max,
				 p->molecule_size,
				 p->boundary_off_centering_tolerance_min,
				 p->boundary_off_centering_tolerance_max,
				 p->boundary_extrapolation_tolerance_min,
				 p->boundary_extrapolation_tolerance_max,
				 p->x,
				 debug,
				 &i_center, &x_rel);
	bool ok = (status == 0)
		  ? ( (status == p->status) && (i_center == p->i_center)
					    && fuzzy_EQ(x_rel, p->x_rel) )
		  :   (status == p->status);

	printf("size=%d off_centering_tol=[%g,%g] extrapolation_tol=[%g,%g] x=%g",
	       p->molecule_size,
	       (double) p->boundary_off_centering_tolerance_min,
	       (double) p->boundary_off_centering_tolerance_max,
	       (double) p->boundary_extrapolation_tolerance_min,
	       (double) p->boundary_extrapolation_tolerance_max,
	       (double) p->x);
	printf(": status=%d", status);
	if (status == 0)
	   then printf(" i_center=%d x_rel=%g: ",
		       i_center, (double) x_rel);

	if (ok)
	   then printf(": ok\n");
	   else {
		++failure_count;
		printf(": FAIL\n");
		if (p->status == 0)
		   then printf("*** expected status=%d i_center=%d x_rel=%g\n",
			       p->status, p->i_center, (double) p->x_rel);
		   else printf("*** expected status=%d\n", p->status);
		if (status == 0)
		   then printf("*** got      status=%d i_center=%d x_rel=%g\n",
			       status, i_center, (double) x_rel);
		   else printf("*** got      status=%d\n", status);
		}
	}

return failure_count;
}
