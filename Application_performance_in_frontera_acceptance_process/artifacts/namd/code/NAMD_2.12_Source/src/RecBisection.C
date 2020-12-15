/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include <math.h>
#include <stdlib.h>

#include "RecBisection.h"
#include "PatchMap.inl"
#include "Patch.h"
#include "PatchMgr.h"

/* ********************************************************************* */
/* Constructor for the RecBisection Class                                */
/* ********************************************************************* */
RecBisection::RecBisection(int numpartitions, PatchMap *thePatchMap)
{
    patchMap    = thePatchMap;
    npartition  = numpartitions;
    numPatches  = patchMap->numPatches();
    partitions  = new Partition[npartition];
    patchload   = new PatchLoad[numPatches];
    currentp    = 0;

    if ( partitions == NULL ||
	 patchload == NULL )
    {
      NAMD_die("memory allocation failed in RecBisection::RecBisection");
    }


    // the cost coeffiencients that is used to compute the load introduced
    // to the processor by a patch

    c_local0    = 1;
    c_local1    = 0.015;
    c_edge0     = 0.001;
    c_edge1     = 0.001;
    c_icompute0 = 0.001;
    c_icompute1 = 0.000035;

    //c_local0    = 1.0;
    //c_local1    = 0.;
    //c_edge0     = 0.;
    //c_edge1     = 0.;
    //c_icompute0 = 0.;
    //c_icompute1 = 0.;
}


/* ********************************************************************* */
/* Destructor for the RecBisection Class                                */
/* ********************************************************************* */
RecBisection::~RecBisection()
{
    delete [] partitions; 
    delete [] patchload;
}



/* *********************************************************************** */
/* This is a recursive function to partition a 3-D mesh into n cubical     */
/* subpartitions with approximately equal loads.                           */
/*                                                                         */
/* Input  n : total number of subpartitions that the partition "p" has to  */
/*            be divided.  Its value is any integer > 0                    */
/*        p : the current partition                                        */
/*                                                                         */
/* The partition p is divided into two partitions (say p1 and p2) such     */
/* that a) p1 and p2 equally loaded b) p1 and p2 further to be partitioned */
/* ino n1 and n2 subpartitions.                                            */
/* If n is even then n1=n2=n/2. Otherwise n1 is one more than n2.          */
/*                                                                         */
/* Since the subpartitions are rectangualr prisms (not an artibrary shape),*/
/* it is not always possible to find a bisection point where the load      */
/* is equally divided.                                                     */
/* The following strategy is used to get a good partitioning:              */
/* We divide the initial partition along x, y and z directions tentatively */
/* and choose the one that gives the best division                          */
/* *********************************************************************** */

void RecBisection::rec_divide(int n, const Partition &p)
{
    int       i=0,j=0,k=0;              // general purpose index vars
    int       posi[3],posj[3],posk[3];  // division points along x,y,z direction
    int       mindir;                   // the best direction
    int       n1, n2;                   // number of subpartitions in p1 and p2
    int       p1_empty[3], p2_empty[3]; // true if a subpartition is empty
    float     load1;                   // desired load of the first subpartition
    float     prevload,currentload;    
    float     loadarray[3];            // actual loads of p1 (for each x,y,z
                                       // division)
    float     diff;                    // load diffenrence (actual - desired)
    float     mindiff;                 // minimum difference (among directions)
    Partition p1;                      // first subpartition
    Partition p2;                      // second subpartition


    if (n==1)
    {
       // no further subdivision
       // record teh partition p as a final partition
       partitions[currentp++] = p;
       return;
    }

    // calculate division ratio: 1/2 iff n is even, otherwise
    // first partition has more load

    n2 = n/2;
    n1 = n-n2;

    load1  = ( (float) n1/(float) n ) * p.load;

    for(i=XDIR; i<=ZDIR; i++) {p1_empty[i] = p2_empty[i] = 0;}

    p1 = p;
    p2 = p;

    // now try dividing along the x,y,z directions 
    
    // try x-axis 
    currentload = 0.0;
    i = p.origin.x;
    while(currentload < load1 && i<= p.corner.x) 
      {
        prevload = currentload;
        for(j=p.origin.y; j<=p.corner.y; j++)
          for(k=p.origin.z; k<=p.corner.z; k++)
            currentload += patchload[patchMap->pid(i,j,k)].total;
        if (currentload > load1)
           if ( prev_better(prevload,currentload,load1) ) 
               {currentload=prevload; break;}
        i++; 
      }
    posi[XDIR] = i; posj[XDIR] = j; posk[XDIR] = k;
    loadarray[XDIR] = currentload; 
    if (i == p.origin.x) p1_empty[XDIR] = 1;
    if (i > p.corner.x)  p2_empty[XDIR] = 1;


    // z axis
    currentload = 0.0;
    k = p.origin.z;
    while(currentload < load1 && k <= p.corner.z)
      {
        prevload = currentload;
        for(i=p.origin.x; i<=p.corner.x; i++)
          for(j=p.origin.y; j<=p.corner.y; j++)
            currentload += patchload[patchMap->pid(i,j,k)].total;
        if (currentload > load1)
           if ( prev_better(prevload,currentload,load1) )
               {currentload=prevload; break;}
        k++;
      }

    posi[ZDIR] = i; posj[ZDIR] = j; posk[ZDIR] = k;
    loadarray[ZDIR] = currentload; 
    if (k == p.origin.z) p1_empty[ZDIR] = 1;
    if (k > p.corner.z)  p2_empty[ZDIR] = 1;
    

    // y axis
    currentload = 0.0;
    j = p.origin.y;
    while(currentload < load1 && j <= p.corner.y) 
      {
        prevload = currentload;
        for(i=p.origin.x; i<=p.corner.x; i++)
          for(k=p.origin.z; k<=p.corner.z; k++)
            currentload += patchload[patchMap->pid(i,j,k)].total;
        if (currentload > load1)
           if ( prev_better(prevload,currentload,load1) )
               {currentload=prevload; break;}
        j++; 
      }
    posi[YDIR] = i; posj[YDIR] = j; posk[YDIR] = k;
    loadarray[YDIR] = currentload; 
    if (j == p.origin.y) p1_empty[YDIR] = 1;
    if (j > p.corner.y)  p2_empty[YDIR] = 1;

    // determine the best division direction
    mindiff = load1;
    mindir   = -1;
    for(i=XDIR; i<=ZDIR; i++) { 
       diff =  load1 - loadarray[i];
       if (diff < 0.0) diff = -diff;
       if (mindiff >= diff) {mindiff = diff; mindir = i;}
    }

    // always divide along x or y if possible
    int lx = p.corner.x - p.origin.x + 1;
    if ( n >= 2*lx || lx >= 2*n || n > 10 || lx > 10 ) {
      int n2x = n * (p.load - loadarray[XDIR]) / p.load + 0.5;
      if ( lx > 1 && n2x > 0 && n2x < n ) mindir = XDIR;
    }

    // revise n1 and n2 based on selected splitting dimension
    n2 = n * (p.load - loadarray[mindir]) / p.load + 0.5;
    if ( n2 < 1 ) n2 = 1;
    if ( n2 >= n ) n2 = n-1;
    n1 = n-n2;

    // divide along mindir
    switch (mindir) {
      case XDIR: p1.corner.x = posi[XDIR] - 1;
                 p2.origin.x = posi[XDIR];
                 break;
      case YDIR: p1.corner.y = posj[YDIR] - 1;
                 p2.origin.y = posj[YDIR];
                 break;
      case ZDIR: p1.corner.z = posk[ZDIR] - 1;
                 p2.origin.z = posk[ZDIR];
                 break;
      default:   NAMD_bug("RecBisection failing horribly!");
    }
    p1.load = loadarray[mindir];
    p2.load = p.load - p1.load;
    if (!p1_empty[mindir]) rec_divide(n1,p1); 
    if (!p2_empty[mindir]) rec_divide(n2,p2);
}


/* ************************************************************************ */
/* Compute the initial overhead/load of each patch to the processor. The    */
/* load of patches are computed as follows: Each patch has a fixed cost and */
/* a variable cost depending of number atoms it has; c_local0 and c_local1  */
/* respectively. Secondly, due to interaction with each patch, the patch    */
/* incurs some load determined by c_edge0 and c_edge1 cost parameters.      */
/* Finally, each patch has load due to computing forces between its certain */
/* neighbours, with cost coefficients c_icompute0 and c_icompute1           */
/* ************************************************************************ */

void RecBisection::compute_patch_load() 
{  
   int   i,nix,neighbour;
   int   numAtoms, numFixed, numNeighAtoms, numNeighFixed;
   float total_icompute;

   for(i=0; i<numPatches; i++) { 
#ifdef MEM_OPT_VERSION
     numAtoms = patchMap->numAtoms(i);
     numFixed = patchMap->numFixedAtoms(i);
#else
     numAtoms      = patchMap->patch(i)->getNumAtoms();
     numFixed      = patchMap->patch(i)->getNumFixedAtoms();
#endif
     patchload[i].total = 0.0;
     patchload[i].edge  = 0.0;

     patchload[i].local = c_local0 + c_local1 * numAtoms;
     patchload[i].local += c_icompute0 +
	c_icompute1*(numAtoms*numAtoms-numFixed*numFixed);


     total_icompute = 0.0;

     PatchID neighbors[PatchMap::MaxOneAway + PatchMap::MaxTwoAway];
     int nNeighbors = patchMap->oneAwayNeighbors(i,neighbors);
     
     for(nix=0; nix<nNeighbors; nix++) {
       neighbour = neighbors[nix];
#ifdef MEM_OPT_VERSION      
       numNeighAtoms = patchMap->numAtoms(neighbour);
       numNeighFixed = patchMap->numFixedAtoms(neighbour);
#else
       numNeighAtoms = patchMap->patch(neighbour)->getNumAtoms();
       numNeighFixed = patchMap->patch(neighbour)->getNumFixedAtoms();
#endif
       patchload[i].icompute[nix] = 
	 c_icompute0 +
	 c_icompute1*(numNeighAtoms*numAtoms-numNeighFixed*numFixed);
       total_icompute += patchload[i].icompute[nix]; 
       patchload[i].edge += c_edge0 + c_edge1 * numNeighAtoms; 
     }
     patchload[i].total+=patchload[i].local+total_icompute+patchload[i].edge; 
   }
}




/* *********************************************************************  */
/* Partitions a 3D space. First a recursive algorithm divides the initial */
/* space into rectangular prisms with approximately equal loads.          */
/* Then these rectangular prisms ar emodified to firther increase the     */
/* load balance and reduce communication cost                             */
/* *********************************************************************  */

int RecBisection::partition(int *dest_arr)
{
    int i;
  
    top_partition.origin.x = 0; 
    top_partition.origin.y = 0; 
    top_partition.origin.z = 0; 
    top_partition.corner.x  = patchMap->gridsize_a()-1;
    top_partition.corner.y  = patchMap->gridsize_b()-1;
    top_partition.corner.z  = patchMap->gridsize_c()-1;
    top_partition.load      = 0.0;

    // calculate estimated computational load due to each patch
    compute_patch_load();

    for(i=0; i<numPatches; i++) top_partition.load += patchload[i].total;

    // divide into rectangular prisms with load as equal as possible
    rec_divide(npartition,top_partition);

    if (currentp != npartition) 
          return 0;
    else  {
      if (dest_arr==NULL)
	  assignNodes();
      else
	  assign_nodes_arr(dest_arr);
    }

    return 1;
}


/* ********************************************************************* */
/* partitioning done, update PatchDistib data structure to assign the    */
/* patches to nodes. (patches in partition i is assigned to node i)      */
/* ********************************************************************* */

void RecBisection::assignNodes()
{
    int i,j,k,pix;
    Partition p;

    for(pix=0; pix<npartition; pix++)
    {
       p = partitions[pix];
       for(i=p.origin.x; i<=p.corner.x; i++)
          for(j=p.origin.y; j<=p.corner.y; j++)
             for(k=p.origin.z; k<=p.corner.z; k++)
                patchMap->assignNode(patchMap->pid(i,j,k),pix);  
    }

    int npatches = patchMap->numPatches();
    for (int pid=0; pid<npatches; ++pid ) 
      patchMap->assignBaseNode(pid);
}

/* ********************************************************************* */
/* partitioning done, save results to an array, rather than updating     */
/* the PatchDistib data structure.  Otherwise, thist is identical to     */
/* assignNodes()                                                         */
/* ********************************************************************* */

void RecBisection::assign_nodes_arr(int *dest_arr)
{
    int i,j,k,pix;
    Partition p;

    for(pix=0; pix<npartition; pix++)
    {
       p = partitions[pix];
       for(i=p.origin.x; i<=p.corner.x; i++)
          for(j=p.origin.y; j<=p.corner.y; j++)
             for(k=p.origin.z; k<=p.corner.z; k++)  {
                dest_arr[patchMap->pid(i,j,k)] = pix;  
	      }
    }
}


/* ************************************************************************ */
/* to be implemented:                                                       */
/* this functions revisits edge directions to refine the load distribution  */
/* ************************************************************************ */
void RecBisection::refine_edges()
{
}


/* ************************************************************************ */
/* to be implemented:                                                       */
/* this function refines the boundaries of subpartitions to redice the      */
/* communication across processors                                          */
/* ************************************************************************ */
void RecBisection::refine_boundaries()
{
}


/* ************************************************************************ */
/* to be implemented:                                                       */
/* refine boundries invokes this function for eeach surface                 */
/* ************************************************************************ */
void RecBisection::refine_surface()
{
}

/* ************************************************************************ */
/* return true if the difference between previous load (prev1) and desired  */
/* load (load1) is less than  teh difference between current an desired     */
/* ************************************************************************ */

int RecBisection::prev_better(float prev, float current, float load1)
{
   float diff1,diff2;

   diff1 = load1 - prev;
   diff2 = current - load1;

   if (diff1 < 0.0) diff1 = -diff1;
   if (diff2 < 0.0) diff2 = -diff2;

   return (diff1 <= diff2);
}

#if USE_TOPOMAP 

/* *********************************************************************  */
/* Partitions a 3D processor into rectangular prisms of different
   sizes corresponding to the patches. Then a processor in the prism
   is assigned to hold the patch.
*************************************/

int RecBisection::partitionProcGrid(int X, int Y, int Z, int *dest_arr) {
    
    int i = 0;
    top_partition.origin.x = 0; 
    top_partition.origin.y = 0; 
    top_partition.origin.z = 0; 
    top_partition.corner.x  = patchMap->gridsize_a()-1;
    top_partition.corner.y  = patchMap->gridsize_b()-1;
    top_partition.corner.z  = patchMap->gridsize_c()-1;
    top_partition.load      = 0.0;

    // calculate estimated computational load due to each patch
    compute_patch_load();
   
    for(i=0; i<numPatches; i++) top_partition.load += patchload[i].total;
    
    int new_X = X, new_Y = Y, new_Z = Z;
    DimensionMap dm;

    dm.x = 0; dm.y = 1; dm.z = 2;
    
    findOptimalDimensions(X, Y, Z, new_X, new_Y, new_Z, dm);
    
    Partition proc_p;
    
    proc_p.origin.x = 0;
    proc_p.corner.x = new_X - 1;
    
    proc_p.origin.y = 0;
    proc_p.corner.y = new_Y - 1;
    
    proc_p.origin.z = 0;
    proc_p.corner.z = new_Z - 1;
    
    iout << "\nLDB: Partitioning Proc Grid of size " << "[" << X << "] [" << Y << "] [" << Z << "]\n" << endi;

    // divide processor grid into rectangular prisms whose sizes
    // corrspond to the computation load of the patches. Moreover
    // neoghboring should also be on nearby processors on the grid
    int rc = topogrid_rec_divide(proc_p, top_partition);
    
    if (rc < 0)
      return rc;

    assignPatchesToProcGrid(dest_arr, X, Y, Z, dm);
    
    return 1;
}

//Partition both a processor grid and a patch grid. Only useful when
//number of processors is 2 times greater than number of patches. It
//returns a processor partition for each patch.

int RecBisection::topogrid_rec_divide(Partition &proc_p, Partition &patch_p) {
  
    Partition proc_p1, proc_p2, patch_p1, patch_p2;
    int i=0, j=0, k=0;
    int posi[3],posj[3],posk[3];  // division points along x,y,z direction
    int mindir;                   // the best direction
    int p1_empty[3], p2_empty[3]; // true if a subpartition is empty
    double loadarray[3];          // actual loads of p1 (for each x,y,z division)
    double diff;                  // load diffenrence (actual - desired)
    double mindiff;               // minimum difference (among directions)
    
    proc_p1 = proc_p2 = proc_p;
    patch_p1 = patch_p2 = patch_p;

    p1_empty[0] = p1_empty[1] = p1_empty[2] = 0;
    p2_empty[0] = p2_empty[1] = p2_empty[2] = 0;    
    /*
    CkPrintf("topo rec divide pe grid = (%d,%d,%d) to (%d,%d,%d) and patch grid = (%d,%d,%d) to (%d,%d,%d)\n", proc_p.origin.x, proc_p.origin.y, 
	     proc_p.origin.z, proc_p.corner.x, proc_p.corner.y, proc_p.corner.z, 
	     patch_p.origin.x, patch_p.origin.y, patch_p.origin.z, patch_p.corner.x, 
	     patch_p.corner.y, patch_p.corner.z);
    */

    //We have just one patch left in the partition
    if(patch_p.origin.x == patch_p.corner.x && 
       patch_p.origin.y == patch_p.corner.y && 
       patch_p.origin.z == patch_p.corner.z) {
	
        int pid = patchMap->pid(patch_p.origin.x,
				patch_p.origin.y, patch_p.origin.z);
      
	//This patch owns all processors in the partition proc_p
	partitions[pid] = proc_p;
	
	//CkPrintf("Assigning patch %d to partition with origin at %d,%d,%d\n", pid, 
	// proc_p.origin.x, 
	// proc_p.origin.y, proc_p.origin.z);

	return 1;
    }  
    
    
    if(proc_p.origin.x == proc_p.corner.x && 
       proc_p.origin.y == proc_p.corner.y && 
       proc_p.origin.z == proc_p.corner.z) {
      
      return -1;
      /*
	for(k = patch_p.origin.z; k < patch_p.corner.z; k++) 
	for(j = patch_p.origin.y; j < patch_p.corner.y; j++) 
	for(i = patch_p.origin.x; i < patch_p.corner.x; i++) {
	
	partitions[patchMap->pid(i,j,k)] = proc_p;
	}	
	return;
      */
    }
    
    double load_x, load_y, load_z;
    load_x = (int) (proc_p.corner.x - proc_p.origin.x + 1)/2;
    load_x /= proc_p.corner.x - proc_p.origin.x + 1;
    load_x *= patch_p.load;
    
    
    load_y = (int) (proc_p.corner.y - proc_p.origin.y + 1)/2;
    load_y /= proc_p.corner.y - proc_p.origin.y + 1;
    load_y *= patch_p.load;
    
    load_z = (int) (proc_p.corner.z - proc_p.origin.z + 1)/2;
    load_z /= proc_p.corner.z - proc_p.origin.z + 1;
    load_z *= patch_p.load;
    
    loadarray[XDIR] = loadarray[YDIR] = loadarray[ZDIR] = 10.0 * patch_p.load;
    
    double currentload = 0;
    double prevload = 0;
    if(load_x > 0 && patch_p.origin.x != patch_p.corner.x) { 
	//Split the processors in X dimension
	//Also split patches in X dimension
	
	currentload = 0.0;
	i = patch_p.origin.x;
	
	while(currentload < load_x && i <= patch_p.corner.x) {
	    prevload = currentload;
	    for(j=patch_p.origin.y; j<=patch_p.corner.y; j++)
		for(k=patch_p.origin.z; k<=patch_p.corner.z; k++)
		    currentload += patchload[patchMap->pid(i,j,k)].total;
	    
	    if (currentload > load_x)
		if ( prev_better(prevload,currentload,load_x) ) { 
		    currentload=prevload; 
		    break;
		}
	    i++; 
	}
	
	posi[XDIR] = i; posj[XDIR] = j; posk[XDIR] = k;
	loadarray[XDIR] = currentload; 
	if (i == patch_p.origin.x) p1_empty[XDIR] = 1;
	if (i > patch_p.corner.x)  p2_empty[XDIR] = 1;
    }

    if(load_z > 0 && patch_p.origin.z != patch_p.corner.z) { //Split patches in Z dimension
        //Also split patches in Z dimension
	
	// z axis
	currentload = 0.0;
	k = patch_p.origin.z;
	while(currentload < load_z && k <= patch_p.corner.z){
	    prevload = currentload;
	    for(i=patch_p.origin.x; i<=patch_p.corner.x; i++)
		for(j=patch_p.origin.y; j<=patch_p.corner.y; j++)
		    currentload += patchload[patchMap->pid(i,j,k)].total;
	    
	    if (currentload > load_z)
		if ( prev_better(prevload,currentload,load_z) )
		    {currentload=prevload; break;}
	    k++;
	}
	
	posi[ZDIR] = i; posj[ZDIR] = j; posk[ZDIR] = k;
	loadarray[ZDIR] = currentload; 
	if (k == patch_p.origin.z) p1_empty[ZDIR] = 1;
	if (k > patch_p.corner.z)  p2_empty[ZDIR] = 1;
    }
    
    if(load_y > 0 && patch_p.origin.y != patch_p.corner.y) { //Y dimension
	
	// y axis
	currentload = 0.0;
	j = patch_p.origin.y;
	while(currentload < load_y && j <= patch_p.corner.y) {
	    prevload = currentload;
	    for(i=patch_p.origin.x; i<=patch_p.corner.x; i++)
		for(k=patch_p.origin.z; k<=patch_p.corner.z; k++)
		    currentload += patchload[patchMap->pid(i,j,k)].total;
	    
	    if (currentload > load_y)
		if ( prev_better(prevload,currentload,load_y) )
		    {currentload=prevload; break;}
	    j++; 
	}
	posi[YDIR] = i; posj[YDIR] = j; posk[YDIR] = k;
	loadarray[YDIR] = currentload; 
	if (j == patch_p.origin.y) p1_empty[YDIR] = 1;
	if (j > patch_p.corner.y)  p2_empty[YDIR] = 1;
    }
    
    //Need to choose the best division. Need to be careful because the
    //processor grid may not be partitionable in all dimensions
    
    mindiff  = 10.0 * patch_p.load;
    mindir   = -1;
    for(i=ZDIR; i >= XDIR; i--) { 
	
	double load1 = 0;
	if(i == XDIR)
	    load1 = load_x;
	else if(i == YDIR)
	    load1 = load_y;
	else if(i == ZDIR)
	    load1= load_z;
	
	diff =  load1 - loadarray[i];
	if (diff < 0.0) diff = -diff;
	if (mindiff > diff) {mindiff = diff; mindir = i;}
    }

    int pdir = 0;

    //Give up on loads and divide it some way.  Later change to the
    //max of the differences along processors and patches
    if(mindiff >= patch_p.load) {
	if(patch_p.origin.x != patch_p.corner.x) {
	    mindir = XDIR;
	    posi[XDIR] = (patch_p.origin.x + patch_p.corner.x + 1)/2;
	}
	else if(patch_p.origin.y != patch_p.corner.y) {
	    mindir = YDIR;
	    posj[YDIR] = (patch_p.origin.y + patch_p.corner.y + 1)/2;
	}
	else {
	    mindir = ZDIR;
	    posk[ZDIR] = (patch_p.origin.z + patch_p.corner.z + 1)/2;
	}
	
	if(load_x > 0) 
	  pdir = XDIR;
	else if(load_y > 0) 
	  pdir = YDIR;
	else 	
	  pdir = ZDIR;
	    	
	loadarray[mindir] = 0.5 * patch_p.load;
	loadarray[pdir] = 0.5 * patch_p.load;
    }
    else {
	pdir = mindir;
    }

    int n1 = 0, n2 = 0;

    // divide processor along pdir
    switch (pdir) {
    case XDIR: 
        n1 = loadarray[XDIR] * (proc_p.corner.x - proc_p.origin.x) / patch_p.load;
	n2 = proc_p.corner.x - proc_p.origin.x - n1;
	
	proc_p1.corner.x = proc_p.origin.x + n1;
	proc_p2.origin.x = proc_p.origin.x + n1 + 1;
	
	break;
    case YDIR: 
	
	n1 = loadarray[YDIR] * (proc_p.corner.y - proc_p.origin.y) / patch_p.load;
	n2 = proc_p.corner.y - proc_p.origin.y - n1;
	
	proc_p1.corner.y = proc_p.origin.y + n1;
	proc_p2.origin.y = proc_p.origin.y + n1 + 1;	
	break;
	
    case ZDIR: 
	n1 = loadarray[ZDIR] * (proc_p.corner.z - proc_p.origin.z) / patch_p.load;
	n2 = proc_p.corner.z - proc_p.origin.z - n1;
	
	proc_p1.corner.z = proc_p.origin.z + n1;
	proc_p2.origin.z = proc_p.origin.z + n1 + 1;
	
	break;
    default:   NAMD_bug("RecBisection failing horribly!");
    }    

    // divide PATCH along mindir
    switch (mindir) {
    case XDIR: 
	patch_p1.corner.x = posi[XDIR] - 1;
	patch_p2.origin.x = posi[XDIR];
    
	break;
    case YDIR: 
	patch_p1.corner.y = posj[YDIR] - 1;
	patch_p2.origin.y = posj[YDIR];
	break;
	
    case ZDIR: 
	patch_p1.corner.z = posk[ZDIR] - 1;
	patch_p2.origin.z = posk[ZDIR];
	break;

    default:   NAMD_bug("RecBisection failing horribly!");
    }
  
    patch_p1.load = loadarray[mindir];
    patch_p2.load = patch_p.load - patch_p1.load;

    /*
    CkPrintf("Calling Topo Rec Divide along %d direction with pe grid = (%d,%d,%d) to (%d,%d,%d) and patch grid = (%d,%d,%d) to (%d,%d,%d)\n", mindir, proc_p1.origin.x, proc_p1.origin.y, 
	     proc_p1.origin.z, proc_p1.corner.x, proc_p1.corner.y, proc_p1.corner.z, 
	     patch_p1.origin.x, patch_p1.origin.y, patch_p1.origin.z, patch_p1.corner.x, 
	     patch_p1.corner.y, patch_p1.corner.z);
    
    CkPrintf("Calling Topo Rec Divide along %d direction with pe grid = (%d,%d,%d) to (%d,%d,%d) and patch grid = (%d,%d,%d) to (%d,%d,%d)\n\n\n", mindir, proc_p2.origin.x, proc_p2.origin.y, 
	     proc_p2.origin.z, proc_p2.corner.x, proc_p2.corner.y, proc_p2.corner.z, 
	     patch_p2.origin.x, patch_p2.origin.y, patch_p2.origin.z, patch_p2.corner.x, 
	     patch_p2.corner.y, patch_p2.corner.z);
    */

    int rc = 0;

    if (!p1_empty[mindir]) 
      rc = topogrid_rec_divide(proc_p1, patch_p1); 
    else
      return -1;

    if (rc < 0)
      return rc;

    if (!p2_empty[mindir]) 
      rc = topogrid_rec_divide(proc_p2, patch_p2);
    else
      return -1;

    if (rc < 0)
      return rc;    

    return 1;
}


//Assign all patches to their corresponding processors
void RecBisection::assignPatchesToProcGrid(int *dest_arr, int X, int Y, int Z,
					   DimensionMap dm)
{
    int pix;
    Partition p;
    TopoManager tmgr;
 
    srand(CkNumPes() * 1000);
    
    int coord[3];

    //Naive scheme where I just choose the origin of the proc grie
    for(pix=0; pix<npartition; pix++) {
        p = partitions[pix];
	//Get the actual origin of that processor prism and assign it
	//to this patch
	
	int xdiff, ydiff, zdiff;

	xdiff = p.corner.x - p.origin.x + 1;
	ydiff = p.corner.y - p.origin.y + 1;
	zdiff = p.corner.z - p.origin.z + 1;

	coord[0] = p.origin.x + rand() % xdiff;
	coord[1] = p.origin.y + rand() % ydiff;
	coord[2] = p.origin.z + rand() % zdiff;
	
	int pe = tmgr.coordinatesToRank(coord[dm.x], coord[dm.y], coord[dm.z], 0);
	patchMap->assignNode(pix, pe);  
	dest_arr[pix] = pe;  
    }

    int npatches = patchMap->numPatches();
    for (int pid=0; pid<npatches; ++pid ) 
      patchMap->assignBaseNode(pid);
}

#endif //USE_TOPOMAP
