/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "SortAtoms.h"
#include "NamdTypes.h"
#include <algorithm>

// #include "charm++.h"


struct sortop_base {
  const FullAtom * const a;
  sortop_base(const FullAtom* atoms) : a(atoms) { }
};

struct sortop_x : public sortop_base {
  sortop_x(const FullAtom* atoms) : sortop_base(atoms) { }
  bool operator() (int i, int j) const {
    return ( a[i].position.x < a[j].position.x );
  }
};

struct sortop_y : public sortop_base {
  sortop_y(const FullAtom* atoms) : sortop_base(atoms) { }
  bool operator() (int i, int j) const {
    return ( a[i].position.y < a[j].position.y );
  }
};

struct sortop_z : public sortop_base {
  sortop_z(const FullAtom* atoms) : sortop_base(atoms) { }
  bool operator() (int i, int j) const {
    return ( a[i].position.z < a[j].position.z );
  }
};


static void partition(int *order, const FullAtom *atoms, int begin, int end) {

  //  Applies orthogonal recursive bisection with splittings limited
  //  to multiples of 32 for warps and a final split on multiples of 16.

  int split;
  // must be a multiple of 32 or 16 between begin and end to split at
  if ( begin/32 < (end-1)/32 ) {
    // find a multiple of 32 near the median
    split = ((begin + end + 32) / 64) * 32;
  } else if ( begin/16 < (end-1)/16 ) {
    // find a multiple of 16 near the median
    split = ((begin + end + 16) / 32) * 16;
  } else {
    return;
  }

  BigReal xmin, ymin, zmin, xmax, ymax, zmax;
  {
    const Position &pos = atoms[order[begin]].position;
    xmin = pos.x;
    ymin = pos.y;
    zmin = pos.z;
    xmax = pos.x;
    ymax = pos.y;
    zmax = pos.z;
  }
  for ( int i=begin+1; i<end; ++i ) {
    const Position &pos = atoms[order[i]].position;
    if ( pos.x < xmin ) { xmin = pos.x; }
    if ( pos.y < ymin ) { ymin = pos.y; }
    if ( pos.z < zmin ) { zmin = pos.z; }
    if ( pos.x > xmax ) { xmax = pos.x; }
    if ( pos.y > ymax ) { ymax = pos.y; }
    if ( pos.z > zmax ) { zmax = pos.z; }
  }
  xmax -= xmin;
  ymax -= ymin;
  zmax -= zmin;

#define NTH_ELEMENT(BEGIN,SPLIT,END,OP) std::nth_element(BEGIN,SPLIT,END,OP)
#if defined(NAMD_CUDA) && defined(__GNUC_PATCHLEVEL__)
#if __GNUC__ == 4 && __GNUC_MINOR__ == 8 && __GNUC_PATCHLEVEL__ == 2
#define NTH_ELEMENT(BEGIN,SPLIT,END,OP) std::sort(BEGIN,END,OP)
#warning gcc 4.8.2 std::nth_element would segfault (see http://gcc.gnu.org/bugzilla/show_bug.cgi?id=58800)
#endif
#endif

  if ( xmax >= ymax && xmax >= zmax ) {
    NTH_ELEMENT(order+begin, order+split, order+end, sortop_x(atoms));
  } else if ( ymax >= xmax && ymax >= zmax ) {
    NTH_ELEMENT(order+begin, order+split, order+end, sortop_y(atoms));
  } else {
    NTH_ELEMENT(order+begin, order+split, order+end, sortop_z(atoms));
  }

  if ( split & 16 ) return;

  // recursively partition before and after split
  partition(order, atoms, begin, split);
  partition(order, atoms, split, end);

}

#if defined(NAMD_CUDA) && defined(__GNUC_PATCHLEVEL__)
#if __GNUC__ == 4 && __GNUC_MINOR__ == 8 && __GNUC_PATCHLEVEL__ == 2
// #error gcc 4.8.2 std::nth_element would segfault (see http://gcc.gnu.org/bugzilla/show_bug.cgi?id=58800)
#endif
#endif

void sortAtomsForCUDA(int *order, const FullAtom *atoms, int nfree, int n) {

  // partition free atoms
  // CkPrintf("%d %d\n", 0, nfree);
  partition(order, atoms, 0, nfree);

  // partition fixed atoms
  // CkPrintf("%d %d\n", nfree, n);
  partition(order, atoms, nfree, n);

}

void sortAtomsForPatches(int *order, int *breaks,
                         const FullAtom *atoms, int nmgrps, int natoms,
                         int ni, int nj, int nk) {

//CkPrintf("sorting %d atoms in %d groups to %d x %d x %d\n",
//    natoms, nmgrps, nk, nj, ni);
  std::sort(order, order+nmgrps, sortop_z(atoms));
  int pid = 0;
  int ibegin = 0;
  int nai = 0;
  for ( int ip=0; ip < ni; ++ip ) {
    int naj = nai;
    int targi = nai + (natoms - nai - 1) / (ni - ip) + 1;
    int iend;
    for ( iend=ibegin; iend<nmgrps; ++iend ) { 
      int mgs = atoms[order[iend]].migrationGroupSize;
      if (nai + mgs <= targi) nai += mgs;
      else break;
    }
//CkPrintf("  Z %d %d (%d) %d\n", ibegin, iend, iend-ibegin, nai);
    std::sort(order+ibegin, order+iend, sortop_y(atoms));
    int jbegin = ibegin;
    for ( int jp=0; jp < nj; ++jp ) {
      int nak = naj;
      int targj = naj + (nai - naj - 1) / (nj - jp) + 1;
      int jend;
      for ( jend=jbegin; jend<iend; ++jend ) { 
        int mgs = atoms[order[jend]].migrationGroupSize;
        if (naj + mgs <= targj) naj += mgs;
        else break;
      }

//CkPrintf("    Y %d %d (%d) %d\n", jbegin, jend, jend-jbegin, naj);
      std::sort(order+jbegin, order+jend, sortop_x(atoms));
      int kbegin = jbegin;
      for ( int kp=0; kp < nk; ++kp ) {
        int targk = nak + (naj - nak - 1) / (nk - kp) + 1;
        int kend;  
        for ( kend=kbegin; kend<jend; ++kend ) {
          int mgs = atoms[order[kend]].migrationGroupSize;
          if (nak + mgs <= targk) nak += mgs;
          else break;
//CkPrintf("        atom %d %d %.2f\n", atoms[order[kend]].id, mgs,
//                  atoms[order[kend]].position.x);
        }
//CkPrintf("      X %d %d (%d) %d\n", kbegin, kend, kend-kbegin, nak);
        breaks[pid++] = kend;
        kbegin = kend;
      }
      jbegin = jend;
    }
    ibegin = iend;
  }

}

