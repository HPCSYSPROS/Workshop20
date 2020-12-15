/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef RECBISECTION_H
#define RECBISECTION_H

#include "converse.h"

class PatchMap;

#if USE_TOPOMAP 
/******
       NAMD likes the X dimension to be the largest, followed by Y and
       then Z. This structure stores the relationsip between x,y,z and
       virtualx, virtualy and virtualz.
***************/

struct DimensionMap {
    int x;
    int y; 
    int z;
};


inline void findOptimalDimensions(int X, int Y, int Z, 
				 int & new_X, int & new_Y, int &new_Z,
				 DimensionMap &dm) {
    if(X == Y && Y == Z)
	return;
    
    if(X >= Y) {
	if(X >= Z) {
	    new_X = X;
	    dm.x = 0;
	    
	    if(Z >= Y) {
		new_Y = Z;
		new_Z = Y;
		
		dm.y = 2;
		dm.z = 1;
	    }
	    else {
		new_Y = Y;
		new_Z = Z;

		dm.y = 1;
		dm.z = 2;
	    }
	}
	else {
	    new_X = Z;
	    new_Y = X;
	    new_Z = Y;

	    dm.x = 1;
	    dm.y = 2;
	    dm.z = 0;
	}
    }
    else {
	if(Y >= Z) {
	    new_X = Y;
	    dm.y  = 0;

	    if(Z >= X) {
		new_Y = Z;
		new_Z = X;

		dm.x  = 2;
		dm.z  = 1;
	    }
	    else {
		new_Y = X;
		new_Z = Z;

		dm.x  = 1;
		dm.z  = 2;		
	    }
	}
	else {
	    new_X = Z;
	    new_Y = Y;
	    new_Z = X;

	    dm.x = 2;
	    dm.y = 1;
	    dm.z = 0;
	}
    }
}
#endif

/* *********************************************************************** */
/* This class performs a recursive coordinate bisection partitioning       */
/* together with some communication and computation refinements            */
/* *********************************************************************** */

#define MAXNEIGHBOUR 26
class RecBisection 
{
    private:

      typedef struct {                          // a rectangular prism
         float  load;                           // is represented with
         struct { int x,y,z; } origin;          // origin and corner coordinates
         struct { int x,y,z; } corner;
      } Partition;


      typedef struct {                          // the cost of a patch
         float total;                           // is represented here
         float local;
         float edge;
         float icompute[MAXNEIGHBOUR];
      } PatchLoad;

      enum {XDIR=0,YDIR,ZDIR};

      // cost parameters 
      float c_local0;       // fixed cost per patch
      float c_local1;       // cost per atom in the patch
      float c_edge0;        // fixed cost for having a neighbor patch
      float c_edge1;        // cost per atom of the neighboring patch
      float c_icompute0;    // fixed cost per calc. forces for a neighbor
      float c_icompute1;    // cost per atom of the neihgbor that I calc force.

 

      int          numPatches;
      int          npartition; 
      int          currentp;
      Partition    *partitions;

      PatchLoad    *patchload;
      PatchMap *patchMap;     
      Partition    top_partition;

      void compute_patch_load();              // determine cost of each patch
      void rec_divide(int, const Partition&); // recursive  partitioning 
      void assignNodes();                     // assign partitions to nodes
      void assign_nodes_arr(int *);           // assign partitions to array
      void refine_edges();                    
      void refine_boundaries();
      void refine_surface();
      int  prev_better(float,float,float);   

    public:

      RecBisection(int, PatchMap *);

      ~RecBisection();
      int partition(int *);                    // perform partitioning.
					       // if parameter=NULL, store
					       // results in patchDistrib,
					       // otherwise, store in array
      
#if USE_TOPOMAP 
      RecBisection(int, int , int, PatchMap *);  //Pass in a 3d
						 //processor grid
      void assignPatchesToProcGrid(int *dest_arr, int X, int Y, int Z, 
				   DimensionMap dm);
      int topogrid_rec_divide(Partition &proc_p, Partition &patch_p);     
      int partitionProcGrid(int X, int Y, int Z, int *dest_arr);
#endif
};

#endif

