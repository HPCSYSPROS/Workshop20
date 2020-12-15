/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Common operations for ComputeNonbonded classes
*/

#ifndef COMPUTENONBONDEDINL_H
#define COMPUTENONBONDEDINL_H

#ifndef NAMD_RESTRICT
#define restrict
#endif

#include "ComputeNonbondedUtil.h"
#include "SimParameters.h"
#include "Molecule.h"
#include "LJTable.h"
#include "ReductionMgr.h"
#include "ReserveArray.h"

#include "PressureProfile.h"

// BEGIN LA
#include "Random.h"
// END LA

#include <algorithm>

#ifdef NAMD_KNL
inline int pairlist_from_pairlist_knl(float cutoff2,
                                  float p_i_x, float p_i_y, float p_i_z,
                                  const CompAtomFlt *p_j,
                                  const plint *list, int list_size, int *newlist,
                                  float *r2list, float *xlist, float *ylist, float *zlist) {

  int *nli = newlist;
  float *r2i = r2list;
  float *xli = xlist;
  float *yli = ylist;
  float *zli = zlist;

  if ( list_size <= 0) return 0;

  int g = 0;

  float p_j_x, p_j_y, p_j_z;
  float x2, y2, z2, r2;
#pragma vector aligned
#pragma ivdep
  for ( g = 0 ; g < list_size; ++g ) {
    int gi=list[g];
    p_j_x = p_j[ gi ].position.x;
    p_j_y = p_j[ gi ].position.y;
    p_j_z = p_j[ gi ].position.z;

    x2 = p_i_x - p_j_x;
    r2 = x2 * x2;
    y2 = p_i_y - p_j_y;
    r2 += y2 * y2;
    z2 = p_i_z - p_j_z;
    r2 += z2 * z2;

    if ( r2 <= cutoff2 ) {
      *nli = gi;  ++nli;
      *r2i = r2;  ++r2i;
      *xli = x2;  ++xli;
      *yli = y2;  ++yli;
      *zli = z2;  ++zli;
    }
  }

  return nli - newlist;
}
#endif// NAMD_KNL

inline int pairlist_from_pairlist(BigReal cutoff2,
				  BigReal p_i_x, BigReal p_i_y, BigReal p_i_z,
				  const CompAtom *p_j,
				  const plint *list, int list_size, int *newlist,
				  BigReal r2_delta, BigReal *r2list) {
  
  BigReal cutoff2_delta = cutoff2 + r2_delta;
  int *nli = newlist;
  BigReal *r2i = r2list;

  if ( list_size <= 0) return 0;

  int g = 0;

#ifdef NAMD_KNL

  BigReal p_j_x, p_j_y, p_j_z;
  BigReal x2, y2, z2, r2;
#pragma vector aligned
#pragma ivdep
  for ( g = 0 ; g < list_size; ++g ) {
    int gi=list[g];
    p_j_x = p_j[ gi ].position.x;
    p_j_y = p_j[ gi ].position.y;
    p_j_z = p_j[ gi ].position.z;

    x2 = p_i_x - p_j_x;
    r2 = x2 * x2 + r2_delta;
    y2 = p_i_y - p_j_y;
    r2 += y2 * y2;
    z2 = p_i_z - p_j_z;
    r2 += z2 * z2;

    if ( r2 <= cutoff2_delta ) {
      *nli = gi;  ++nli;
      *r2i = r2;  ++r2i;
    }
  }

#else // NAMD_KNL
#ifndef SIMPLE_PAIRLIST
#ifdef  A2_QPX
  //***************************************************************
  //* 4-way unrolled and software-pipelined and optiized via QPX
  //***************************************************************  

  int jout = 0;
  if ( list_size > 16) {
    // prefetch
    int jcur0 = list[g];
    int jcur1 = list[g + 1];
    int jcur2 = list[g + 2];
    int jcur3 = list[g + 3];
    
    int j0, j1, j2, j3;
    
    vector4double    pj_v_0, pj_v_1, pj_v_2, pj_v_3; 
    vector4double    v_0, v_1, v_2, v_3;
    register BigReal r2_0, r2_1, r2_2, r2_3;
    
    vector4double p_i_v = {p_i_x, p_i_y, p_i_z, 0.};
    vector4double r2_delta_v = {r2_delta};

    pj_v_0 = vec_ld(jcur0 * sizeof(CompAtom), (BigReal *)p_j);
    pj_v_1 = vec_ld(jcur1 * sizeof(CompAtom), (BigReal *)p_j);  
    pj_v_2 = vec_ld(jcur2 * sizeof(CompAtom), (BigReal *)p_j);  
    pj_v_3 = vec_ld(jcur3 * sizeof(CompAtom), (BigReal *)p_j);  
    
    for ( g = 4 ; g < list_size - 4; g += 4 ) {
      // compute 1d distance, 4-way parallel      
      //Save the previous iterations values, gives more flexibility 
      //to the compiler to schedule the loads and the computation
      j0   =   jcur0;           j1   =   jcur1;
      j2   =   jcur2;           j3   =   jcur3;
      
      jcur0  =  list[g    ];    jcur1  =  list[g + 1];
      jcur2  =  list[g + 2];    jcur3  =  list[g + 3];

      __dcbt((void*)(p_j + jcur0));

      v_0 = vec_sub (p_i_v, pj_v_0);
      v_1 = vec_sub (p_i_v, pj_v_1);
      v_2 = vec_sub (p_i_v, pj_v_2);
      v_3 = vec_sub (p_i_v, pj_v_3);
      
      v_0 = vec_madd (v_0, v_0, r2_delta_v);
      v_1 = vec_madd (v_1, v_1, r2_delta_v);
      v_2 = vec_madd (v_2, v_2, r2_delta_v);
      v_3 = vec_madd (v_3, v_3, r2_delta_v);

      pj_v_0 = vec_ld(jcur0 * sizeof(CompAtom), (BigReal *)p_j);
      pj_v_1 = vec_ld(jcur1 * sizeof(CompAtom), (BigReal *)p_j);  
      pj_v_2 = vec_ld(jcur2 * sizeof(CompAtom), (BigReal *)p_j);  
      pj_v_3 = vec_ld(jcur3 * sizeof(CompAtom), (BigReal *)p_j);  

      r2_0 = vec_extract(v_0, 0) + vec_extract(v_0, 1) + vec_extract(v_0, 2);
      r2_1 = vec_extract(v_1, 0) + vec_extract(v_1, 1) + vec_extract(v_1, 2);
      r2_2 = vec_extract(v_2, 0) + vec_extract(v_2, 1) + vec_extract(v_2, 2);
      r2_3 = vec_extract(v_3, 0) + vec_extract(v_3, 1) + vec_extract(v_3, 2);
      
      size_t test0, test1, test2, test3;
      size_t jout0, jout1, jout2, jout3;
      
      test0 = ( r2_0   <   cutoff2_delta );
      test1 = ( r2_1   <   cutoff2_delta );
      test2 = ( r2_2   <   cutoff2_delta );
      test3 = ( r2_3   <   cutoff2_delta );
      
      jout0 = jout;
      nli[ jout0 ]  = j0;         r2i[ jout0 ] = r2_0;
      jout += test0;              jout1 = jout;

      nli[ jout1 ]  = j1;         r2i[ jout1 ] = r2_1;
      jout += test1;              jout2 = jout;

      nli[ jout2 ]  = j2;         r2i[ jout2 ] = r2_2;
      jout += test2;              jout3 = jout;

      nli[ jout3 ]  = j3;         r2i[ jout3 ] = r2_3;
      jout += test3;
    }
    g -= 4;
  }

  nli += jout;
  r2i += jout;  
#else
  //***************************************************************
  //* 4-way unrolled and software-pipelined 
  //***************************************************************

  int jout = 0;
  if ( list_size > 16) {
    // prefetch
    int jcur0 = list[g];
    int jcur1 = list[g + 1];
    int jcur2 = list[g + 2];
    int jcur3 = list[g + 3];
    
    int j0, j1, j2, j3;
    
    register  BigReal pj_x_0, pj_x_1, pj_x_2, pj_x_3; 
    register  BigReal pj_y_0, pj_y_1, pj_y_2, pj_y_3; 
    register  BigReal pj_z_0, pj_z_1, pj_z_2, pj_z_3; 
    
    register BigReal t_0, t_1, t_2, t_3, r2_0, r2_1, r2_2, r2_3;
    
    pj_x_0 = p_j[jcur0].position.x;
    pj_x_1 = p_j[jcur1].position.x;  
    pj_x_2 = p_j[jcur2].position.x;  
    pj_x_3 = p_j[jcur3].position.x;  
    pj_y_0 = p_j[jcur0].position.y; 
    pj_y_1 = p_j[jcur1].position.y;  
    pj_y_2 = p_j[jcur2].position.y;  
    pj_y_3 = p_j[jcur3].position.y;  
    pj_z_0 = p_j[jcur0].position.z; 
    pj_z_1 = p_j[jcur1].position.z;
    pj_z_2 = p_j[jcur2].position.z; 
    pj_z_3 = p_j[jcur3].position.z;
    
    for ( g = 4 ; g < list_size - 4; g += 4 ) {
      // compute 1d distance, 4-way parallel
      
      //Save the previous iterations values, gives more flexibility 
      //to the compiler to schedule the loads and the computation
      j0   =   jcur0;           j1   =   jcur1;
      j2   =   jcur2;           j3   =   jcur3;
      
      jcur0  =  list[g    ];    jcur1  =  list[g + 1];
      jcur2  =  list[g + 2];    jcur3  =  list[g + 3];

#ifdef ARCH_POWERPC
      __dcbt ((void *) &p_j[jcur0]);
#endif      

      //Compute X distance
      t_0   =  p_i_x - pj_x_0;   t_1   =  p_i_x - pj_x_1;
      t_2   =  p_i_x - pj_x_2;   t_3   =  p_i_x - pj_x_3;
      
      r2_0  =  t_0 * t_0 + r2_delta; 
      r2_1  =  t_1 * t_1 + r2_delta;
      r2_2  =  t_2 * t_2 + r2_delta;
      r2_3  =  t_3 * t_3 + r2_delta;
      
      //Compute y distance
      t_0    =  p_i_y - pj_y_0;   t_1    =  p_i_y - pj_y_1;
      t_2    =  p_i_y - pj_y_2;   t_3    =  p_i_y - pj_y_3;
      r2_0  +=  t_0 * t_0;        r2_1  +=  t_1 * t_1;
      r2_2  +=  t_2 * t_2;        r2_3  +=  t_3 * t_3;
      
      //compute z distance
      t_0    =  p_i_z - pj_z_0;   t_1    =  p_i_z - pj_z_1;
      t_2    =  p_i_z - pj_z_2;   t_3    =  p_i_z - pj_z_3;
      r2_0  +=  t_0 * t_0;        r2_1  +=  t_1 * t_1;
      r2_2  +=  t_2 * t_2;        r2_3  +=  t_3 * t_3;
      
      pj_x_0 = p_j[jcur0].position.x;
      pj_x_1 = p_j[jcur1].position.x;  
      pj_x_2 = p_j[jcur2].position.x;  
      pj_x_3 = p_j[jcur3].position.x;  
      pj_y_0 = p_j[jcur0].position.y; 
      pj_y_1 = p_j[jcur1].position.y;  
      pj_y_2 = p_j[jcur2].position.y;  
      pj_y_3 = p_j[jcur3].position.y;  
      pj_z_0 = p_j[jcur0].position.z; 
      pj_z_1 = p_j[jcur1].position.z;
      pj_z_2 = p_j[jcur2].position.z; 
      pj_z_3 = p_j[jcur3].position.z;
      
      bool test0, test1, test2, test3;
      
      test0 = ( r2_0   <   cutoff2_delta );
      test1 = ( r2_1   <   cutoff2_delta );
      test2 = ( r2_2   <   cutoff2_delta );
      test3 = ( r2_3   <   cutoff2_delta );
      
      int jout0, jout1, jout2, jout3;

      jout0 = jout;
      nli[ jout0 ]  = j0;         r2i[ jout0 ] = r2_0;
      jout += test0;              jout1 = jout;
      nli[ jout1 ]  = j1;         r2i[ jout1 ] = r2_1;
      jout += test1;              jout2 = jout;
      nli[ jout2 ]  = j2;         r2i[ jout2 ] = r2_2;
      jout += test2;              jout3 = jout;
      nli[ jout3 ]  = j3;         r2i[ jout3 ] = r2_3;

      jout += test3;
    }
    g -= 4;
  }

  nli += jout;
  r2i += jout;  
#endif
#endif

  int j2 = list[g];
  BigReal p_j_x = p_j[j2].position.x;
  BigReal p_j_y = p_j[j2].position.y;
  BigReal p_j_z = p_j[j2].position.z;
  while ( g < list_size ) {
    int j = j2;
    j2 = list[++g];
    BigReal t2 = p_i_x - p_j_x;
    BigReal r2 = t2 * t2 + r2_delta;
    p_j_x = p_j[j2].position.x;
    t2 = p_i_y - p_j_y;
    r2 += t2 * t2;
    p_j_y = p_j[j2].position.y;
    t2 = p_i_z - p_j_z;
    r2 += t2 * t2;
    p_j_z = p_j[j2].position.z;
    if ( r2 <= cutoff2_delta ) {
      *nli= j; ++nli;
      *r2i = r2; ++r2i;
    }
  }

#endif // NAMD_KNL

  return nli - newlist;
}

// clear all
// define interaction type (pair or self)
#define NBPAIR	1
#define NBSELF	2



// Various atom sorting functions
#if NAMD_ComputeNonbonded_SortAtoms != 0

inline void sortEntries_selectionSort(SortEntry * const se, const int seLen) {

  register int i;

  for (i = 0; i < seLen; i++) {

    // Search through the remaining elements, finding the lowest
    //   value, and then swap it with the first remaining element.
    //   Start by assuming the first element is the smallest.
    register int smallestIndex = i;
    register BigReal smallestValue = se[i].sortValue;
    register int j;
    for (j = i + 1; j < seLen; j++) {
      register BigReal currentValue = se[j].sortValue;
      if (currentValue < smallestValue) {
        smallestIndex = j;
        smallestValue = currentValue;
      }
    }

    // Swap the first remaining element with the smallest element
    if (smallestIndex != i) {
      register SortEntry* entryA = se + i;
      register SortEntry* entryB = se + smallestIndex;
      register unsigned int tmpIndex = entryA->index;
      register BigReal tmpSortValue = entryA->sortValue;
      entryA->index = entryB->index;
      entryA->sortValue = entryB->sortValue;
      entryB->index = tmpIndex;
      entryB->sortValue = tmpSortValue;
    }
  }
}

inline void sortEntries_bubbleSort(SortEntry * const se, const int seLen) {

  register int keepSorting = 0;

  do {

    // Reset the keepSorting flag (assume no swaps will occur)
    keepSorting = 0;

    // Loop through the pairs and swap if needed
    register SortEntry* sortEntry1 = se;
    for (int i = 1; i < seLen; i++) {

      register SortEntry* sortEntry0 = sortEntry1;
      sortEntry1 = se + i;
      register BigReal sortEntry0_sortValue = sortEntry0->sortValue;
      register BigReal sortEntry1_sortValue = sortEntry1->sortValue;

      if (sortEntry0_sortValue > sortEntry1_sortValue) {
        register int sortEntry0_index = sortEntry0->index;
        register int sortEntry1_index = sortEntry1->index;
        sortEntry0->index = sortEntry1_index;
        sortEntry0->sortValue = sortEntry1_sortValue;
        sortEntry1->index = sortEntry0_index;
        sortEntry1->sortValue = sortEntry0_sortValue;
        keepSorting = 1;
      }
    }

  } while (keepSorting != 0);  // Loop again if at least one set of
                               //   elements was swapped.
}

// NOTE: The buf parameter should point to a buffer that this function
//   can use as a temp storage location (scratch pad).
// NOTE: This function may swap the values of se and buf, but will not
//   return any other value (nor will it set either to NULL).
inline void sortEntries_mergeSort_v1(SortEntry * &se, SortEntry * &buf, int seLen) {

  register SortEntry* srcArray = se;
  register SortEntry* dstArray = buf;

  // Start with each element being a separate list.  Start
  //   merging the "lists" into larger lists.
  register int subListSize = 1;
  while (subListSize < seLen) {

    // NOTE: This iteration consumes sublists of length
    //   subListSize and produces sublists of length
    //   (2*subListSize).  So, keep looping while the length of a
    //   single sorted sublist is not the size of the entire array.

    // Iterate through the lists, merging each consecutive pair of lists.
    register int firstListOffset = 0;
    while (firstListOffset < seLen) {

      /// Setup pointers and counts for sublists in the pair. ///

      register int numElements = std::min(2 * subListSize, seLen - firstListOffset);
      register int list0len;
      register int list1len;
      if (numElements > subListSize) {
        list0len = subListSize;                // First list full
        list1len = numElements - subListSize;  // 1+ elements in second list
      } else {
        list0len = numElements;                // 1+ elements in first list
        list1len = 0;                          // Zero elements in second list
      }

      register SortEntry* list0ptr = srcArray + firstListOffset;
      register SortEntry* list1ptr = list0ptr + subListSize;
      register SortEntry* dstptr = dstArray + firstListOffset;

      /// Merge the sublists ///

      // While there are elements in both lists, pick from one
      while (list0len > 0 && list1len > 0) {

        register BigReal sortValue0 = list0ptr->sortValue;
        register BigReal sortValue1 = list1ptr->sortValue;

        if (sortValue0 < sortValue1) {  // choose first list (list0)

          // Copy the value from srcArray to dstArray
          register int index0 = list0ptr->index;
          dstptr->sortValue = sortValue0;
          dstptr->index = index0;

          // Move the pointers forward for the sublists
          dstptr++;
          list0ptr++;
          list0len--;

        } else {                        // choose second list (list1)

          // Copy the value from srcArray to dstArray
          register int index1 = list1ptr->index;
          dstptr->sortValue = sortValue1;
          dstptr->index = index1;

          // Move the pointers forward for the sublists
          dstptr++;
          list1ptr++;
          list1len--;
        }

      } // end while (list0len > 0 && list1len > 0)

      // NOTE: Either list0len or list1len is zero at this point
      //   so only one of the following loops should execute.

      // Drain remaining elements from the first list (list0)
      while (list0len > 0) {

        // Copy the value from srcArray to dstArray
        register BigReal sortValue0 = list0ptr->sortValue;
        register int index0 = list0ptr->index;
        dstptr->sortValue = sortValue0;
        dstptr->index = index0;

        // Move the pointers forward for the sublists
        dstptr++;
        list0ptr++;
        list0len--;

      } // end while (list0len > 0)

      // Drain remaining elements from the first list (list1)
      while (list1len > 0) {

        // Copy the value from srcArray to dstArray
        register BigReal sortValue1 = list1ptr->sortValue;
        register int index1 = list1ptr->index;
        dstptr->sortValue = sortValue1;
        dstptr->index = index1;

        // Move the pointers forward for the sublists
        dstptr++;
        list1ptr++;
        list1len--;

      } // end while (list1len > 0)

      // Move forward to the next pair of sub-lists
      firstListOffset += (2 * subListSize);

    } // end while (firstListOffset < seLen) {

    // Swap the dstArray and srcArray pointers
    register SortEntry* tmpPtr = dstArray;
    dstArray = srcArray;
    srcArray = tmpPtr;

    // Double the subListSize
    subListSize <<= 1;

  }  // end while (subListSize < seLen)

  // Set the sort values pointers (NOTE: srcArray and dstArray are
  //   swapped at the end of each iteration of the merge sort outer-loop).
  buf = dstArray;
  se = srcArray;
}

// NOTE: The buf parameter should point to a buffer that this function
//   can use as a temp storage location (scratch pad).
// NOTE: This function may swap the values of se and buf, but will not
//   return any other value (nor will it set either to NULL).
inline void sortEntries_mergeSort_v2(SortEntry * &se, SortEntry * &buf, int seLen) {

  // NOTE: This macro "returns" either val0 (if test == 0) or val1 (if
  // test == 1).  It expects test to be either 0 or 1 (no other values).
  #define __TERNARY_ASSIGN(test, val0, val1)   ((test * val0) + ((1 - test) * val1))

  register SortEntry* srcArray = se;
  register SortEntry* dstArray = buf;

  // Start with each element being a separate list.  Start
  //   merging the "lists" into larger lists.
  register int subListSize = 1;
  while (subListSize < seLen) {

    // NOTE: This iteration consumes sublists of length
    //   subListSize and produces sublists of length
    //   (2*subListSize).  So, keep looping while the length of a
    //   single sorted sublist is not the size of the entire array.

    // Iterate through the lists, merging each consecutive pair of lists.
    register int firstListOffset = 0;
    while (firstListOffset < seLen) {

      /// Setup pointers and counts for sublists in the pair. ///

      // Calculate the number of elements for both sublists...
      //   min(2 * subListSize, seLen - firstListOffset);
      register int numElements;
      {
        register int numElements_val0 = 2 * subListSize;
        register int numElements_val1 = seLen - firstListOffset;
        register bool numElements_test = (numElements_val0 < numElements_val1);
        numElements = __TERNARY_ASSIGN(numElements_test, numElements_val0, numElements_val1);
      }

      // Setup the pointers for the source and destination arrays
      register SortEntry* dstptr = dstArray + firstListOffset;    // destination array pointer
      register SortEntry* list0ptr = srcArray + firstListOffset;  // source list 0 pointer
      register SortEntry* list1ptr = list0ptr + subListSize;      // source list 1 pointer
      register SortEntry* list0ptr_end;  // pointer to end of source list0's elements (element after last)
      register SortEntry* list1ptr_end;  // pointer to end of source list1's elements (element after last)
      {
        register bool lenTest = (numElements > subListSize);
        register int list0len_val0 = subListSize;
        register int list1len_val0 = numElements - subListSize;
        register int list0len_val1 = numElements;  // NOTE: list1len_val1 = 0
        register int list0len = __TERNARY_ASSIGN(lenTest, list0len_val0, list0len_val1);
        register int list1len = __TERNARY_ASSIGN(lenTest, list1len_val0, 0);
        // avoid pre-load of sortValue1 from past end of array
        if ( ! lenTest ) list1ptr = list0ptr;
        list0ptr_end = list0ptr + list0len;
        list1ptr_end = list1ptr + list1len;
      }

      // The firstListOffset variable won't be used again until the next
      //   iteration, so go ahead and update it now...
      //   Move forward to the next pair of sub-lists
      firstListOffset += (2 * subListSize);

      /// Merge the sublists ///

      // Pre-load values from both source arrays
      register BigReal sortValue0 = list0ptr->sortValue;
      register BigReal sortValue1 = list1ptr->sortValue;
      register int index0 = list0ptr->index;
      register int index1 = list1ptr->index;

      // While both lists have at least one element in them, compare the
      //   heads of each list and place the smaller of the two in the
      //   destination array.
      while (list0ptr < list0ptr_end && list1ptr < list1ptr_end) {

        // Compare the values
        register bool test = (sortValue0 < sortValue1);

        // Place the "winner" in the destination array
        dstptr->sortValue = __TERNARY_ASSIGN(test, sortValue0, sortValue1);
        dstptr->index = __TERNARY_ASSIGN(test, index0, index1);
        dstptr++;

        // Update the pointers
        list0ptr += __TERNARY_ASSIGN(test, 1, 0);
        list1ptr += __TERNARY_ASSIGN(test, 0, 1);

        // Refill the sortValue and index register
        // NOTE: These memory locations are likely to be in cache
        sortValue0 = list0ptr->sortValue;
        sortValue1 = list1ptr->sortValue;
        index0 = list0ptr->index;
        index1 = list1ptr->index;

      } // end while (list0ptr < list0ptr_end && list1ptr < list1ptr_end)

      // NOTE: At this point, at least one of the lists is empty so no
      //   more than one of the loops will be executed.

      // Drain the remaining elements from list0
      while (list0ptr < list0ptr_end) {

        // Place the value into the destination array
        dstptr->sortValue = sortValue0;
        dstptr->index = index0;
        dstptr++;

        // Load the next entry in list0
        list0ptr++;
        sortValue0 = list0ptr->sortValue;
        index0 = list0ptr->index;

      } // end while (list0ptr < list0ptr_end)

      // Drain the remaining elements from list1
      while (list1ptr < list1ptr_end) {

        // Place the value into the destination array
	dstptr->sortValue = sortValue1;
        dstptr->index = index1;
        dstptr++;

        // Load the next entry in list1
        list1ptr++;
        sortValue1 = list1ptr->sortValue;
        index1 = list1ptr->index;

      } // end while (list1ptr < list1ptr_end)

    } // end while (firstListOffset < seLen) {

    // Swap the dstArray and srcArray pointers
    register SortEntry* tmpPtr = dstArray;
    dstArray = srcArray;
    srcArray = tmpPtr;

    // Double the subListSize
    subListSize <<= 1;

  }  // end while (subListSize < seLen)

  // Set the sort values pointers (NOTE: srcArray and dstArray are
  //   swapped at the end of each iteration of the merge sort outer-loop).
  buf = dstArray;
  se = srcArray;

  #undef __TERNARY_ASSIGN
}

#endif  // NAMD_ComputeNonbonded_SortAtoms != 0


#endif // COMPUTENONBONDEDINL_H

