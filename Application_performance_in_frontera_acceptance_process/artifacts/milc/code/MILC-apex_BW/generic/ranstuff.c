/******************* ranstuff.c *********************/
/* MIMD version 7 */
/*  c language random number generator for parallel processors */
/*  exclusive or of feedback shift register and integer congruence
    generator.  Use a different multiplier on each generator, and make sure
    that fsr is initialized differently on each generator.  */

// updated for 144^3*288 lattice.  long long int in interior, so
// multipliers won't overflow, but still ordinary ints in external
// interface.  For this size lattice, the index won't overflow.
// Also increases the integer congruence generator to 64 bits

#include "generic_includes.h"
#include "../include/random.h"

/* usage:
   int seed;
   Real x;
   initialize_prn(prn_pt,seed,index)
	prn_pt is address of a struct double_prn, index selects algorithm
   x = myrand(prn_pt);
   ...
*/

void initialize_prn(double_prn *prn_pt, int seed, int index) {
    /* "index" selects which random number generator - which multiplier */
    if(sizeof(unsigned long long)<8) {
      node0_printf("Type long long less than 8 bytes on this machine. Rand problems?\n");
      exit(1);
    }
    seed = (69607+8*index)*seed+12345;
    prn_pt->r0 = (seed>>8) & 0xffffff;
    seed = (69607+8*index)*seed+12345;
    prn_pt->r1 = (seed>>8) & 0xffffff;
    seed = (69607+8*index)*seed+12345;
    prn_pt->r2 = (seed>>8) & 0xffffff;
    seed = (69607+8*index)*seed+12345;
    prn_pt->r3 = (seed>>8) & 0xffffff;
    seed = (69607+8*index)*seed+12345;
    prn_pt->r4 = (seed>>8) & 0xffffff;
    seed = (69607+8*index)*seed+12345;
    prn_pt->r5 = (seed>>8) & 0xffffff;
    seed = (69607+8*index)*seed+12345;
    prn_pt->r6 = (seed>>8) & 0xffffff;
    seed = (69607+8*index)*seed+12345;
    prn_pt->ic_state = seed;
    prn_pt->multiplier = 100000005 + 8*index;
    prn_pt->addend = 12345;
    prn_pt->scale = 1.0/((Real)0x1000000);
}

Real myrand(double_prn *prn_pt) {
  unsigned long long int s;
  int t;

    t = ( ((prn_pt->r5 >> 7) | (prn_pt->r6 << 17)) ^
	  ((prn_pt->r4 >> 1) | (prn_pt->r5 << 23)) ) & 0xffffff;
    prn_pt->r6 = prn_pt->r5;
    prn_pt->r5 = prn_pt->r4;
    prn_pt->r4 = prn_pt->r3;
    prn_pt->r3 = prn_pt->r2;
    prn_pt->r2 = prn_pt->r1;
    prn_pt->r1 = prn_pt->r0;
    prn_pt->r0 = t;
    s = prn_pt->ic_state * prn_pt->multiplier + prn_pt->addend;
    prn_pt->ic_state = s;
    return( prn_pt->scale*(t ^ ((s>>40)&0xffffff)) );
}

