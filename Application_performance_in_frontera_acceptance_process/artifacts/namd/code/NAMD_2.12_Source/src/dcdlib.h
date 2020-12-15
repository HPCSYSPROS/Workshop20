/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   dcdlib contains C routines used for reading and writing binary
   dcd files.  The format of these files is from binary FORTRAN output,
   so its pretty ugly.  If you're squimish, don't look!!
*/

#ifndef DCDLIB_H
#define DCDLIB_H

#include "largefiles.h"  // must be first!
#include "common.h" // for int32 definition
#include "Vector.h"
#include <stdio.h>
#include <string.h>
#ifndef _NO_MALLOC_H
#include <malloc.h>
#endif
#ifndef WIN32
#include <unistd.h>
#endif
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#ifndef WIN32
#include <pwd.h>
#endif
#include <time.h>
#ifdef WIN32
#include <io.h>
#endif

#ifdef WIN32
#define OFF_T __int64
#else
#define OFF_T off_t
#endif


/*  DEFINE ERROR CODES THAT MAY BE RETURNED BY DCD ROUTINES		*/
#define DCD_DNE		-2	/*  DCD file does not exist		*/
#define DCD_OPENFAILED	-3	/*  Open of DCD file failed		*/
#define DCD_BADREAD 	-4	/*  read call on DCD file failed	*/
#define DCD_BADEOF	-5	/*  premature EOF found in DCD file	*/
#define DCD_BADFORMAT	-6	/*  format of DCD file is wrong		*/
#define DCD_FILEEXISTS  -7	/*  output file already exists		*/
#define DCD_BADMALLOC   -8	/*  malloc failed			*/

/*			FUNCTION ALLUSIONS				*/
int open_dcd_read(char *);      /*  Open a DCD file for reading 	*/
int read_dcdheader(int, int*, int*, int*, int*, double*, int*, int**);	
				/*  Read the DCD header			*/
int read_dcdstep(int, int, float*, float*, float*, int, int, int*);	
				/*  Read a timestep's values		*/
int open_dcd_write(const char *);     /*  Open a DCD file for writing		*/

int write_dcdstep(int, int, float *, float *, float *, double *unitcell);
				/*  Write out a timesteps values	*/
int write_dcdheader(int, const char*, int, int, int, int, int, double, int);	
				/*  Write a dcd header			*/
int get_dcdheader_size(); 
				/* Get the total size of the header */
void close_dcd_read(int, int, int *);
				/*  Close a dcd file open for reading   */
void close_dcd_write(int);	/*  Close a dcd file open for writing   */

int open_dcd_write_par_slave(char *dcdname);
/* Slaves open existing file created by master */
int write_dcdstep_par_cell(int fd, double *cell);
int write_dcdstep_par_XYZUnits(int fd, int N);
     /* Master writes unit cell and natom info */
int update_dcdstep_par_header(int fd); 
     /* Master updates header */

/* Write out a timesteps values partially in parallel for part [parL, parU] */
int write_dcdstep_par_slave(int fd, int parL, int parU, int N, float *X, float *Y, float *Z);
    
/* wrapper for seeking the dcd file */
OFF_T NAMD_seek(int file, OFF_T offset, int whence);

#endif /* ! DCDLIB_H */

