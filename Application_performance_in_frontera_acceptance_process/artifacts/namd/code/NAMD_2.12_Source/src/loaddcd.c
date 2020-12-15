/********************************************************************/
/*                                                                  */
/*   Copyright 1999, Jim Phillips and the University of Illinois.   */
/*                                                                  */
/********************************************************************/

#include "largefiles.h"  /* must be first! */

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <stdio.h>
#include <limits.h>

#ifndef MAP_FILE
#define MAP_FILE 0
#endif

#if ( INT_MAX == 2147483647 )
typedef int     int32;
#else
typedef short   int32;
#endif

int main(int argc, char *argv[]) {

int fd;
struct stat statbuf;
int i, j, isbig, itmp;
off_t n;
double delta;
float delta4;
int32 icntrl[20];
char *ccntrl;
char b[8];
char *d;

if ( argc != 2 ) {
  fprintf(stderr,"This program writes the ICNTRL array to DCD files.\n");
  fprintf(stderr,"Usage: %s <filename> < <data>\n",argv[0]);
  exit(-1);
}

if ( ( fd = open(argv[1], O_RDWR) ) < 0 ) {
  fprintf(stderr,"Can't open %s for updating.\n",argv[1]);
  exit(-1);
}

if ( fstat(fd,&statbuf) < 0 ) {
  fprintf(stderr,"Can't stat %s.\n",argv[1]);
  exit(-1);
}

n = statbuf.st_size;

if ( n <= 104 ) {
  fprintf(stderr,"%s is not in DCD format.\n",argv[1]);
  exit(-1);
}

if ( n % 4 ) {
  fprintf(stderr,"%s is not in DCD format.\n",argv[1]);
  exit(-1);
}

if ( ( d = mmap(0,n,PROT_READ|PROT_WRITE,MAP_FILE|MAP_SHARED,fd,0) )
							== (caddr_t) -1 ) {
  fprintf(stderr,"Can't mmap %s.\n",argv[1]);
  exit(-1);
}

#define SKIPFOUR {d+=4;n-=4;}
#define SKIP(X) {d+=(X);n-=(X);}
#define READINT(X) { X=0; if (isbig) { for(j=0;j<4;++j,X<<8) X+=d[j]; } \
	else { for(j=3;j>=0;--j,X<<8) X+=d[j]; } }

SKIPFOUR;  /* 84 */
SKIPFOUR;  /* "CORD" */

for(j=0;j<9;++j) {
  scanf("%d",&itmp);
  icntrl[j] = itmp;
}

scanf("%f\n",((float*)(icntrl+9)));

for(j=10;j<20;++j) {
  scanf("%d",&itmp);
  icntrl[j] = itmp;
}

ccntrl = (char*)(&(icntrl[0]));

for(j=0;j<80;++j) {
  d[j] = ccntrl[j];
}

}

