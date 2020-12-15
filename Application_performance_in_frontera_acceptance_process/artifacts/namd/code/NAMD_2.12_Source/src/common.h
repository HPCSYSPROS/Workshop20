/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   common definitions for namd.
*/

#ifndef COMMON_H
#define COMMON_H

#if !defined(WIN32) || defined(__CYGWIN__)
#include <unistd.h>
#endif
#include <stdio.h>
#include <limits.h>

#if ( INT_MAX == 2147483647L )
typedef	int	int32;
#elif ( SHRT_MAX == 2147483647L )
typedef	short	int32;
#endif

#ifdef _MSC_VER
typedef __int64 int64;
#else
#if ( INT_MAX == 9223372036854775807LL )
typedef int int64;
#elif ( LONG_MAX == 9223372036854775807LL )
typedef long int64;
#else
typedef long long int64;
#endif
#endif

#if defined(PLACEMENT_NEW)
void * ::operator new (size_t, void *p) { return p; }
#elif defined(PLACEMENT_NEW_GLOBAL)
void * operator new (size_t, void *p) { return p; }
#endif

#define COULOMB 332.0636
#define BOLTZMANN 0.001987191
#define TIMEFACTOR 48.88821
#define PRESSUREFACTOR 6.95E4
#define PDBVELFACTOR 20.45482706
#define PDBVELINVFACTOR (1.0/PDBVELFACTOR)
#define PNPERKCALMOL 69.479

/* Some plagtforms don't have nearbyint or round, so we'll define one */
/* that works everywhere */
#ifdef POWERPC_TANINT
extern "builtin" double __tanint(double); //IEEE round
#define mynearbyint(x)  __tanint(x)
#else
#define mynearbyint(x)  floor((x)+0.5)
#endif

#ifndef PI
#define PI	3.141592653589793
#endif

#ifndef TWOPI
#define TWOPI	2.0 * PI
#endif

#ifndef ONE
#define ONE	1.000000000000000
#endif

#ifndef ZERO
#define ZERO	0.000000000000000
#endif

#ifndef SMALLRAD
#define SMALLRAD      0.0005
#endif

#ifndef SMALLRAD2
#define SMALLRAD2     SMALLRAD*SMALLRAD
#endif

/* Define the size for Real and BigReal.  Real is usually mapped to float */
/* and BigReal to double.  To get BigReal mapped to float, use the 	  */
/* -DSHORTREALS compile time option					  */
typedef float	Real;

#ifdef SHORTREALS
typedef float	BigReal;
#else
typedef double  BigReal;
#endif

#ifndef FALSE
#define FALSE 0
#define TRUE 1
#endif

#ifndef NO
#define NO 0
#define YES 1
#endif

#ifndef STRINGNULL
#define STRINGNULL '\0'
#endif

#define MAX_NEIGHBORS 27

typedef int Bool;

class Communicate;

// global functions
void NAMD_quit(const char *);
void NAMD_die(const char *);
void NAMD_err(const char *);  // also prints strerror(errno)
void NAMD_bug(const char *);
int NAMD_file_exists(const char *filename);
void NAMD_backup_file(const char *filename, const char *extension = 0);
int NAMD_open_text(const char *fname, int append=0);
void NAMD_write(int fd, const char *buf, size_t count, const char *fname = "in NAMD_write()"); // NAMD_die on error
void NAMD_close(int fd, const char *fname);
char *NAMD_stringdup(const char *);
FILE *Fopen(const char *filename, const char *mode);
int  Fclose(FILE *fout);

// message tags
#define SIMPARAMSTAG	100	//  Tag for SimParameters class
#define STATICPARAMSTAG 101	//  Tag for Parameters class
#define MOLECULETAG	102	//  Tag for Molecule class
#define FULLTAG	104
#define FULLFORCETAG 105
#define DPMTATAG 106
#define GRIDFORCEGRIDTAG 107
#define COMPUTEMAPTAG 108

#define CYCLE_BARRIER   0
#define PME_BARRIER     0
#define STEP_BARRIER    0

#define USE_BARRIER   (CYCLE_BARRIER || PME_BARRIER || STEP_BARRIER)


// DMK - Atom Separation (water vs. non-water)
//   Setting this define to a non-zero value will cause the
//   HomePatches to separate the hydrogen groups in their
//   HomePatch::atom lists (all water molecules first, in arbitrary
//   order, followed by all non-waters, in arbitrary order).
#define NAMD_SeparateWaters    0

// DMK - Atom Sort
//   Setting this define to a non-zero value will cause the nonbonded compute
//   objects (pairs, not selfs) to sort the atoms along a line connecting the
//   center of masses of the two patches.  This is only done during timesteps
//   where the pairlists are being generated.  As the pairlist is being
//   generated, once an atom that is far enough away along the line is found,
//   the remaining atoms are automatically skipped (avoiding a distance
//   calculation/check for them).
// NOTE: The "less branches" flag toggles between two versions of merge sort.
//   When it is non-zero, a version that has fewer branches (but more integer
//   math) is used.  This version may or may not be faster or some architectures.
#define NAMD_ComputeNonbonded_SortAtoms                   1
  #define NAMD_ComputeNonbonded_SortAtoms_LessBranches    1

// plf -- alternate water models
#define WAT_TIP3 0
#define WAT_TIP4 1
#define WAT_SWM4 2  /* Drude model (5 charge sites) */


#include "converse.h"

#endif

