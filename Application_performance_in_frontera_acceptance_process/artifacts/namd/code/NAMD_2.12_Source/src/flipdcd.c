/********************************************************************/
/*                                                                  */
/*   Copyright 1998, Jim Phillips and the University of Illinois.   */
/*                                                                  */
/********************************************************************/

/* Modified to include -s[tatus], -B[ig], -L[ittle] options,        */
/* and to process multiple files.                                   */
/* D. Barsky March 2002 LLNL                                        */

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
off_t n;
int i, j, isbig, itmp, argcount=0;
int status_only=0, make_big_only=0, make_little_only=0;
char b[8];
char *d;

if ( argc < 2 ) {
  fprintf(stderr,"This program flips byte-ordering of DCD files.\n");
  usage:
  fprintf(stderr,"Usage: %s [-s] [-B] [-L] file . . . \n",argv[0]);
  fprintf(stderr,"      The default behavior is to flip the byte ordering. Other options are:\n");
  fprintf(stderr,"      -s report the byte-order status of each <file> without changing it\n");
  fprintf(stderr,"      -B make/keep each <file> big-endian\n");
  fprintf(stderr,"      -L make/keep each <file> little-endian\n");
  fprintf(stderr,"      The options are mutually exclusive; the last one read is used.\n\n");
  exit(-1);
}

while (++argcount < argc){
   /* debug: printf("Current argument %d (out of %d): %s\n",argcount,argc-1,argv[argcount]); */ 
   if ((strncmp(argv[argcount],"-S",2) == 0) || (strncmp(argv[argcount],"-s",2) == 0)){
      status_only=1; make_big_only=0; make_little_only=0;
   }
   else if ((strncmp(argv[argcount],"-B",2) == 0) || (strncmp(argv[argcount],"-b",2) == 0)){
      make_big_only=1; status_only=0; make_little_only=0;
   }
   else if ((strncmp(argv[argcount],"-L",2) == 0) || (strncmp(argv[argcount],"-l",2) == 0)){
      make_little_only=1; make_big_only=0; status_only=0; 
   }
   else if (strncmp(argv[argcount],"-",1) == 0){
      printf("\n Error: %s not a valid option. \n\n",argv[argcount]);
      goto usage;
   }
   else{
      if ( ( fd = open(argv[argcount], O_RDWR) ) < 0 ) {
        fprintf(stderr,"Can't open %s for updating. (File must be read/writeable.)\n",argv[argcount]);
        goto end;
      }

      if ( fstat(fd,&statbuf) < 0 ) {
        fprintf(stderr,"Can't stat %s.\n",argv[1]);
        goto end;
      }
      
      n = statbuf.st_size;
      
      if ( n <= 104 ) {
        fprintf(stderr,"%s is not in DCD format.\n",argv[argcount]);
        goto end;
      }

      if ( ( sizeof(char*) < 8 ) && ( n >> 32 ) ) {
        fprintf(stderr,"%s is too large, 64-bit build required\n",argv[argcount]);
        goto end;
      }
      
      if ( n % 4 ) {
        fprintf(stderr,"%s is not in DCD format.\n",argv[argcount]);
        goto end;
      }
      if ( ( d = mmap(0,n,PROT_READ|PROT_WRITE,MAP_FILE|MAP_SHARED,fd,0) )
							      == (caddr_t) -1 ) {
        fprintf(stderr,"Can't mmap %s.\n",argv[argcount]);
        goto end;
      }
      
      if ( status_only ){
        if ( d[0] == 84 ) {
          isbig = 0;
          fprintf(stderr,"%s is little-endian.\n",argv[argcount]);
        }
        else if ( d[3] == 84 ) {
          isbig = 1;
          fprintf(stderr,"%s is big-endian.\n",argv[argcount]);
        }
        else {
          fprintf(stderr,"%s is not in DCD format.\n",argv[argcount]);
        }
        goto end;  /* Done if only status is requested */
        
      }
      else {
        if ( d[0] == 84 ) {
          if ( make_little_only ){
            fprintf(stderr,"%s is already little-endian. (No change made.)\n",argv[argcount]);
            goto end;
          }
          else {   
            isbig = 0;
            fprintf(stderr,"%s was little-endian, will be big-endian.\n",argv[argcount]);
          }
        }
        else if ( d[3] == 84 ) {
          if ( make_big_only ){
            fprintf(stderr,"%s is already big-endian. (No change made.)\n",argv[argcount]);
            goto end;
          }
          else {   
            isbig = 1;
            fprintf(stderr,"%s was big-endian, will be little-endian.\n",argv[argcount]);
          }
        }
        else {
          fprintf(stderr,"%s is not in DCD format.\n",argv[argcount]);
          goto end;
        }
      }
      
#define FLIPFOUR {for(j=0;j<4;++j)b[j]=d[j];for(j=3;j>=0;--j,++d)*d=b[j];n-=4;}
#define FLIPEIGHT {for(j=0;j<8;++j)b[j]=d[j];for(j=7;j>=0;--j,++d)*d=b[j];n-=8;}
#define SKIPFOUR {d+=4;n-=4;}
#define SKIP(X) {d+=(X);n-=(X);}
#define READINT(X) { X=0; if (isbig) { for(j=0;j<4;++j,X<<8) X+=d[j]; } \
	   else { for(j=3;j>=0;--j,X<<8) X+=d[j]; } }


      FLIPFOUR;  /* 84 */
      SKIPFOUR;  /* "CORD" */
      for ( i = 0; i < 20; ++i ) FLIPFOUR;
      FLIPFOUR;  /* 84 */
      FLIPFOUR;  /* TITLE SIZE */
      READINT(itmp); FLIPFOUR;  /* NTITLE */
      if ( n <= (80*(off_t)itmp + 4) ) {
        fprintf(stderr,"%s is too short for DCD format.\n",argv[argcount]);
        goto end;
      }
      SKIP(80*itmp);
      FLIPFOUR;  /* TITLE SIZE */
      
      while ( n ) FLIPFOUR;  /* WHATEVER UNTIL END */
        ;
      
      }
   end:
      ;
   }
}
