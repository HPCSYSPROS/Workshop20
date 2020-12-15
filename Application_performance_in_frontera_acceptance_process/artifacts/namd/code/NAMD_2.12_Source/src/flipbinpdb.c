/********************************************************************/
/*                                                                  */
/*   Copyright 1998, Jim Phillips and the University of Illinois.   */
/*                                                                  */
/********************************************************************/

#include "largefiles.h"  /* must be first! */

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <stdio.h>

#ifndef MAP_FILE
#define MAP_FILE 0
#endif

int main(int argc, char *argv[]) {

int fd;
struct stat statbuf;
off_t i, j, n;
char b[8];
char *d;

if ( argc != 2 ) {
  fprintf(stderr,"This program flips byte-ordering of 8-byte doubles.\n");
  fprintf(stderr,"Usage: %s <filename>\n",argv[0]);
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

if ( (n < 4) || ((n-4) % 24) ) {
  fprintf(stderr,"Size of %s is not 4 plus a multiple of 24.\n",argv[1]);
  exit(-1);
}

if ( ( d = mmap(0,n,PROT_READ|PROT_WRITE,MAP_FILE|MAP_SHARED,fd,0) )
							== (caddr_t) -1 ) {
  fprintf(stderr,"Can't mmap %s.\n",argv[1]);
  exit(-1);
}

for ( j = 0; j < 4; ++j ) b[j] = d[j];
for ( j = 3; j >= 0; --j, ++d ) *d = b[j];

for ( i = 4; i < n; i += 8 ) {
  for ( j = 0; j < 8; ++j ) b[j] = d[j];
  for ( j = 7; j >= 0; --j, ++d ) *d = b[j];
}

exit(0);

}

