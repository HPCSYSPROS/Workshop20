#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include <math.h>
#include <stdio.h>
#include <vector>
#include <deque>

#include "util_Table.h"
#include "loopcontrol.h"
#include <assert.h>
#include "mpi.h"


using namespace std;



// stores the global rho maxima per Carpet component
deque<vector<CCTK_REAL> > global_rho_max_per_comp(0);
deque<vector<CCTK_REAL> > global_rho_max_x_per_comp(0);
deque<vector<CCTK_REAL> > global_rho_max_y_per_comp(0);
deque<vector<CCTK_REAL> > global_rho_max_z_per_comp(0);


bool save_max(const CCTK_REAL rho, 
              const CCTK_REAL x,
              const CCTK_REAL y,
              const CCTK_REAL z,
              const vector<CCTK_REAL>& rho_max, const int n,
              const vector<CCTK_REAL>& x_max,
              const vector<CCTK_REAL>& y_max,
              const vector<CCTK_REAL>& z_max,
              const CCTK_REAL min_separation_between_maxima)
{
        if (rho > rho_max[n]) {
            // RH: this check is redundant as the for loop will shortircuit
            // for n==0 anyway
            if (n > 0) {
               {
                  // the maximum must be at least a given distance away to all other maxima!
                  bool sufficiently_separated = true;
                  for (int l=0; l < n; ++l) {
                     const CCTK_REAL dist = sqrt(pow(x-x_max[l], 2) 
                                               + pow(y-y_max[l], 2)
                                               + pow(z-z_max[l], 2));
                     if (dist < min_separation_between_maxima && rho_max[l] >= 0) {
                        sufficiently_separated = false;
                        break;
                     }
                  }
                  
                  if (sufficiently_separated)
                     return true;
               }
            }
            else {
               // this is the n=0 case, so the maximum we find can immediately be stored!
               return true;
            }
         }
         
         return false;
}

/** 
    Finds N regions with local density maxima.
    The N maxima must be separated by at least some given distance to avoid capturing
    close by values that effectively belong to the same local high density region.
    This routine only finds fragment maxima per Carpet component
*/
extern "C"
void CoreCollapseControl_Find_Fragments_Rho_Max(CCTK_ARGUMENTS) 
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (cctk_time < find_fragments_after_time) {
     return;
  }

  if (cctk_iteration % find_fragments_every != 0)
     return;

  if (find_fragments_only_once && *fragments_found)
     return;

  int const reflevel = GetRefinementLevel (cctkGH);

  if (reflevel > max_reflevel_to_search || reflevel < min_reflevel_to_search)
     return;

  if (verbose) {
     stringstream str;
     str << "Searching for fragments on level " << reflevel;
     CCTK_INFO(str.str().c_str());
  }

  // get new scope to avoid redeclaration of stuff that is
  // buried in the DECLARE macros above
  {
  
  // loop over the grid and collect pair (max,pos)
  
  vector<CCTK_REAL> local_rho_max(N_rho_maxima, -1.0);
  vector<CCTK_REAL> local_rho_max_x(N_rho_maxima, 0.0);
  vector<CCTK_REAL> local_rho_max_y(N_rho_maxima, 0.0);
  vector<CCTK_REAL> local_rho_max_z(N_rho_maxima, 0.0);

  // loop over the grid and collect local maxima
  
  for (int n=0; n < N_rho_maxima; ++n) {
    
    #pragma omp parallel
    {
      LC_LOOP3(CCC_Find_Local_Maxima, i,j,k,
               cctk_nghostzones[0], cctk_nghostzones[1], cctk_nghostzones[2],
               cctk_lsh[0] - cctk_nghostzones[0], cctk_lsh[1] - cctk_nghostzones[1],
               cctk_lsh[2] - cctk_nghostzones[2], cctk_lsh[0], cctk_lsh[1], cctk_lsh[2])
      {
         const CCTK_INT ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
         
         // ignore if below threshold density
         if (rho[ijk] < rho_max_threshold)
            continue;
         
         // ignore if not within search radius of old maximum
         const CCTK_REAL rad = sqrt(pow(x[ijk]-fragment_pos_x[n], 2)
                                  + pow(y[ijk]-fragment_pos_y[n], 2)
                                  + pow(z[ijk]-fragment_pos_z[n], 2));
         if (*fragments_found && rad > search_radius_around_old_maximum)
           continue;
         
         // check if we can add the current rho[ijk] to the local maxima list
         bool save = save_max(rho[ijk], x[ijk], y[ijk], z[ijk],
                              local_rho_max, n,
                              local_rho_max_x,
                              local_rho_max_y,
                              local_rho_max_z,
                              min_separation_between_maxima);
         
         if (save) {
            local_rho_max[n] = rho[ijk];
            local_rho_max_x[n] = x[ijk];
            local_rho_max_y[n] = y[ijk];
            local_rho_max_z[n] = z[ijk];
         }
      }
      LC_ENDLOOP3(CCC_Find_Local_Maxima);
    }
  }
  
  // Now we have to see what other processors on the current refinement level have found!
  // For this, proc 0 collects the data from all other procs
  // and produces a new list.
  
  // set up global rho maxima
  vector<CCTK_REAL> global_rho_max(N_rho_maxima, -1.0);
  vector<CCTK_REAL> global_rho_max_x(N_rho_maxima, 0.0);
  vector<CCTK_REAL> global_rho_max_y(N_rho_maxima, 0.0);
  vector<CCTK_REAL> global_rho_max_z(N_rho_maxima, 0.0);
  
  // first, arrange local data in a flat vector [rho1,x1,y1,z1, rho2,x2,y2,z2,...]
  // for MPI Gather
  vector<CCTK_REAL> sendbuf(4*N_rho_maxima, 0.0);
  for (int n=0; n < N_rho_maxima; ++n) {
     sendbuf[4*n  ] = local_rho_max[n];
     sendbuf[4*n+1] = local_rho_max_x[n];
     sendbuf[4*n+2] = local_rho_max_y[n];
     sendbuf[4*n+3] = local_rho_max_z[n];
  }
  
  
  if (CCTK_MyProc(cctkGH) == 0)
  {
     // this is a a "flat" receive buffer just like the send buffer
     // but now additionially with the local results from all processes!
     // [rho11,x11,y11,z11, rho21,x21,y21,z21, ..., rho12,x12,y12,z12, rho22,x22,y22,z22,...]
     vector<CCTK_REAL> recvbuf(4*N_rho_maxima*CCTK_nProcs(cctkGH), 0.0);
    
     // root process receives data...
     // RH: this will hang if the number of components differ between
     // processes
     // RH: this algorithm will also depend on the processor decomposition
     // since the number of candidates found depdends on the way the
     // components are laid out.
     MPI_Gather(&sendbuf.front(), sendbuf.size(), MPI_DOUBLE, &recvbuf.front(), sendbuf.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
     
     // Now process data from all procs!
     for (int n=0; n < N_rho_maxima; ++n) {
        for (int i=0; i < CCTK_nProcs(cctkGH); ++i) {
           for (int j=0; j < N_rho_maxima; ++j) {
              const CCTK_REAL lrho = recvbuf[4*N_rho_maxima*i + 4*j    ];
              const CCTK_REAL lx   = recvbuf[4*N_rho_maxima*i + 4*j + 1];
              const CCTK_REAL ly   = recvbuf[4*N_rho_maxima*i + 4*j + 2];
              const CCTK_REAL lz   = recvbuf[4*N_rho_maxima*i + 4*j + 3];
              
              // check if we can add the current lrho to the local maxima list
              bool save = save_max(lrho, lx, ly, lz,
                                   global_rho_max, n,
                                   global_rho_max_x,
                                   global_rho_max_y,
                                   global_rho_max_z,
                                   min_separation_between_maxima);
              
              if (save) {
                 global_rho_max[n] = lrho;
                 global_rho_max_x[n] = lx;
                 global_rho_max_y[n] = ly;
                 global_rho_max_z[n] = lz;
              }
           }
        }
     }
     
     // now the global_* vectors contain the maxima on the current reflevel/component!
     // save this until global-mode call to CoreCollapseControl_Set_Fragment_Props
     global_rho_max_per_comp.push_back(global_rho_max);
     global_rho_max_x_per_comp.push_back(global_rho_max_x);
     global_rho_max_y_per_comp.push_back(global_rho_max_y);
     global_rho_max_z_per_comp.push_back(global_rho_max_z);
     
  } else {
     // non-root processes just send data
     MPI_Gather(&sendbuf.front(), sendbuf.size(), MPI_DOUBLE, NULL, 0, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }
  
  // we executed this function at least once for all reflevels!
//  if (reflevel == max_reflevel_to_search)
//     *fragments_found = 1;
  
  } // end scope

}




/**
    Converts the maxima found per component into true global
    maxima.
*/
extern "C"
void CoreCollapseControl_Set_Fragment_Props(CCTK_ARGUMENTS) 
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  if (cctk_time < find_fragments_after_time)
     return;

  if (cctk_iteration % find_fragments_every != 0)
     return;

  if (find_fragments_only_once && *fragments_found)
     return;
  
  if (verbose) {
     stringstream str;
     str << "Getting fragment properties";
     CCTK_INFO(str.str().c_str());
  }

  if (CCTK_MyProc(cctkGH) == 0)
  {
     // root process does the search.
  
     vector<CCTK_REAL> global_rho_max(N_rho_maxima, -1.0);
     vector<CCTK_REAL> global_rho_max_x(N_rho_maxima, 0.0);
     vector<CCTK_REAL> global_rho_max_y(N_rho_maxima, 0.0);
     vector<CCTK_REAL> global_rho_max_z(N_rho_maxima, 0.0);
  
     // Now process data from all Carpet components!
     for (int n=0; n < N_rho_maxima; ++n) {
        for (int i=0; i < global_rho_max_per_comp.size(); ++i) {
           for (int j=0; j < N_rho_maxima; ++j) {
              const CCTK_REAL lrho = global_rho_max_per_comp[i][j];
              const CCTK_REAL lx   = global_rho_max_x_per_comp[i][j];
              const CCTK_REAL ly   = global_rho_max_y_per_comp[i][j];
              const CCTK_REAL lz   = global_rho_max_z_per_comp[i][j];
              
              // check if we can add the current lrho to the local maxima list
              bool save = save_max(lrho, lx, ly, lz,
                                   global_rho_max, n,
                                   global_rho_max_x,
                                   global_rho_max_y,
                                   global_rho_max_z,
                                   min_separation_between_maxima);
              
              if (save) {
                 global_rho_max[n] = lrho;
                 global_rho_max_x[n] = lx;
                 global_rho_max_y[n] = ly;
                 global_rho_max_z[n] = lz;
              }
           }
        }
     }
     
     // set new maxima.
     for (int n=0; n < N_rho_maxima; ++n) {
        // We need to make sure that we actually preserve the numbering of the fragements before!
        // Otherwise our AMR may be screwed up.
        // For this we check whether the respective maximum has not moved by more than a given threshold value.
        // RH for AMR the fragment number does not matter. All that Carpet
        // itself gets to see is the union of all CarpetRegrid2 boxes, so
        // swaping two centers has no effect on what Carpet sees.
        if (*fragments_found) {
           CCTK_REAL dist = 10000000000000;
           CCTK_INT fragment_no = n;
           for (int m=0; m < N_rho_maxima; ++m) {
              // find closest maximum (from previous search) that matches current maximum
              CCTK_REAL my_dist = sqrt(pow(fragment_pos_x[n]-global_rho_max_x[m], 2) 
                                     + pow(fragment_pos_y[n]-global_rho_max_y[m], 2)
                                     + pow(fragment_pos_z[n]-global_rho_max_z[m], 2));
              if (my_dist < dist) {
                 dist = my_dist;
                 fragment_no = m;
              }
           }
           fragment_rho_max[n] = global_rho_max[fragment_no];
           fragment_pos_x[n] = global_rho_max_x[fragment_no];
           fragment_pos_y[n] = global_rho_max_y[fragment_no];
           fragment_pos_z[n] = global_rho_max_z[fragment_no];
        }
        else { // first time we found fragments, we accept them as they are
           fragment_rho_max[n] = global_rho_max[n];
           fragment_pos_x[n] = global_rho_max_x[n];
           fragment_pos_y[n] = global_rho_max_y[n];
           fragment_pos_z[n] = global_rho_max_z[n];
        }
     }
     
     
     // RH: clear() apparently does not free memory, so this saves no memory.
     // Insted one has to swap() with an empty vector.
     global_rho_max_per_comp.clear();
     global_rho_max_x_per_comp.clear();
     global_rho_max_y_per_comp.clear();
     global_rho_max_z_per_comp.clear();
     
     // We finally need to send the result to all other procs.
     // First, we prepare a sendbuffer
     vector<CCTK_REAL> sendbuf(4*N_rho_maxima, 0.0);
     for (int n=0; n < N_rho_maxima; ++n) {
        sendbuf[4*n  ] = fragment_rho_max[n];
        sendbuf[4*n+1] = fragment_pos_x[n];
        sendbuf[4*n+2] = fragment_pos_y[n];
        sendbuf[4*n+3] = fragment_pos_z[n];
     }
     
     // Second, broadcast result to all other procs!
     MPI_Bcast(&sendbuf.front(), sendbuf.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
     
  } else {
     // non-root proc receive results.
     // First, prepare receive buffer
     vector<CCTK_REAL> recvbuf(4*N_rho_maxima, 0.0);
     
     // Second, receive broadcast from root process!
     MPI_Bcast(&recvbuf.front(), recvbuf.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
     
     // put result into "fragment_prop" variables
     for (int n=0; n < N_rho_maxima; ++n) {
        fragment_rho_max[n] = recvbuf[4*n  ];
        fragment_pos_x[n]   = recvbuf[4*n+1];
        fragment_pos_y[n]   = recvbuf[4*n+2];
        fragment_pos_z[n]   = recvbuf[4*n+3];
     }
  }
  
  // check if any fragments have *ever* been found
  for (int n=0; n < N_rho_maxima; ++n) {
     if (fragment_rho_max[n] > 0.0) {
        (*fragments_found) = 1;
        break;
     }
  }
  
  if (!(*fragments_found)) {
    if (track_fragments) {
       for (int n=0; n < N_rho_maxima; ++n) {
          if (which_surface_to_store_fragment_info[n] >= 0) {
             const int sn = which_surface_to_store_fragment_info[n];
             sf_active[sn] = 0;
             sf_valid[sn] = 0;
          }
       }
    }
    return;
  }
  
  // Print info and set SphericalSurface
  stringstream str;
  for (int n=0; n < N_rho_maxima; ++n) {
        str << endl;
        str << "   Fragment " << n << ":  max(rho) = " << fragment_rho_max[n] << endl;
        str << "             :         x = " << fragment_pos_x[n] << endl;
        str << "             :         y = " << fragment_pos_y[n] << endl;
        str << "             :         z = " << fragment_pos_z[n] << endl;
        if (track_fragments) {
          if (which_surface_to_store_fragment_info[n] >= 0) {
            const int sn = which_surface_to_store_fragment_info[n];
            if (fragment_rho_max[n] > 0.0) {
              sf_centroid_x[sn] = fragment_pos_x[n];
              sf_centroid_y[sn] = fragment_pos_y[n];
              if (!ignore_z_tracking[n]) {
                sf_centroid_z[sn] = fragment_pos_z[n];
              } else {
                sf_centroid_z[sn] = 0;
              }
              sf_active[sn] = 1;
              sf_valid[sn] = 1;
            }
          }
        }
  }
  CCTK_INFO(str.str().c_str());
  
  // Enforce symmetric tracking
  if (force_symmetric_tracking)
  {
     if (N_rho_maxima != 2)
        CCTK_WARN(1, "Symmetric tracking only for N_rho_maxima=2!");
     
     const int sn1 = which_surface_to_store_fragment_info[0];
     const int sn2 = which_surface_to_store_fragment_info[1];
     
     const double sx = fabs(sf_centroid_x[sn1] - sf_centroid_x[sn2]) / 2.0;
     const double sy = fabs(sf_centroid_y[sn1] - sf_centroid_y[sn2]) / 2.0;
     const double sz = fabs(sf_centroid_z[sn1] - sf_centroid_z[sn2]) / 2.0;
     sf_centroid_x[sn1] = sx * (sf_centroid_x[sn1] < 0 ? -1.0 : 1.0);
     sf_centroid_y[sn1] = sy * (sf_centroid_y[sn1] < 0 ? -1.0 : 1.0);
     sf_centroid_z[sn1] = sz * (sf_centroid_z[sn1] < 0 ? -1.0 : 1.0);
     
     sf_centroid_x[sn2] = sx * (sf_centroid_x[sn2] < 0 ? -1.0 : 1.0);
     sf_centroid_y[sn2] = sy * (sf_centroid_y[sn2] < 0 ? -1.0 : 1.0);
     sf_centroid_z[sn2] = sz * (sf_centroid_z[sn2] < 0 ? -1.0 : 1.0);
     
  }
  
}



