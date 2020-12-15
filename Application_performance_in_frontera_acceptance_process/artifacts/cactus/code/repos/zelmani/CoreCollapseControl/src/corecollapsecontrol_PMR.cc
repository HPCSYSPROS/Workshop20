#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "carpet.hh"
#define RHOGF 6.1755e17

int CoreCollapseControl_PMRStartup(void) 
{
  const char *banner = "CoreCollapseControl: Progressive Mesh Refinement";
  
  CCTK_RegisterBanner(banner);
  
  return 0;
  
}

void CoreCollapseControl_PMRInit(CCTK_ARGUMENTS) 
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  for (int i=0; i < 10; ++i)
     rho_max_index[i] = 0;
  *global_rho_max = 0.0e0;
  *global_entropy_max = 0.0e0;
  *rho_max = 0.0e0;
  *in_prebounce = 0;

  *bounce = 0;
  *bouncetime = 0.0e0;

  *global_alp_min = 10.0e0;
  *alp_min = 10.0e0;
  *in_preBH = 0;
  *force_check = 0;

  if (find_fragments)
     *fragments_found = 0;

  return;
}

void CoreCollapseControl_PMRGetRhoMax(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  int reduction_handle, varindex, error;
  CCTK_REAL current_rho_max;
  CCTK_REAL current_alp_min;

  if(  (( (cctk_iteration) % rho_max_every) != 0) && !*force_check ) 
    return;

  varindex = CCTK_VarIndex("hydrobase::rho");
  if(varindex < 0) {
    // this should never happen since we REQUIRE HydroBase
    CCTK_ERROR("Variable hydrobase::rho does not exist! Can't regrid");
  }    

  reduction_handle = CCTK_ReductionHandle("maximum");
  assert(reduction_handle >= 0);
  error = CCTK_Reduce(cctkGH,-1,reduction_handle,1,
		      CCTK_VARIABLE_REAL,&current_rho_max,1,varindex);

  if(error != 0) {
      CCTK_WARN(0,"Could not perform reduction on rho_max, no regridding.");
  }

  // RH: what is the difference between rho_max and global_rho_max?
  if(current_rho_max >= *rho_max) {
      *rho_max = current_rho_max;
  }

  if(*rho_max >= *global_rho_max) {
    *global_rho_max = *rho_max;
  }

  varindex = CCTK_VarIndex("admbase::alp");
  if(varindex < 0) {
    // this should never happen since we REQUIRE ADMBase
    CCTK_ERROR("Variable admbase::alp does not exist! Can't steer parameters based on its values.");
  }    
  reduction_handle = CCTK_ReductionHandle("minimum");
  error = CCTK_Reduce(cctkGH,-1,reduction_handle,1,
		      CCTK_VARIABLE_REAL,&current_alp_min,1,varindex);

  if(error != 0) {
      CCTK_WARN(0,"Could not perform reduction on alp_min, no changes to AH/output frequency.");
  }

  // RH: this treats alp_min differently than rho_max
  *alp_min = current_alp_min;
  
  if(*alp_min < *global_alp_min) {
    *global_alp_min = *alp_min;
  }


  return;
} 

void CoreCollapseControl_PMRCheckRegrid(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  CCTK_REAL current_rho_max;

  if( ( (cctk_iteration) % check_every) != 0 ) 
    return;

  CCTK_INT ncenters = 1;
  CCTK_INT do_it_the_usual_way = 1;

  if (CCTK_EQUALS(get_rho_max_from, "global rho max")) {
     // this just deals with one refinement center as it was originally implemented
     ncenters = 1;
     do_it_the_usual_way = 1;
  } else if (CCTK_EQUALS(get_rho_max_from, "fragments")) {
     ncenters = N_rho_maxima;
     do_it_the_usual_way = 0;
     if (!(*fragments_found))
        return;
  }
  
  
  for (int k=0; k < ncenters; ++k) {
    
    int n = PMR_ref_center_to_handle[k];  // refinement centre number of CarpetRegrid2
    if (do_it_the_usual_way) {
      // this just deals with one refinement center as it was originally implemented
      current_rho_max = *rho_max;
    }
    else {
      current_rho_max = fragment_rho_max[n];
    }

    // don't do anything if we don't have a density set at this rho_max_index
    // RH: this should simply test for 0. which does not suffer from the usual
    // issues with floating point comparison since 0. exactly representable and
    // also we want to check if someone sets precisely this value
    if(rho_max_list[rho_max_index[n]] < 1.0e-10) 
       continue;

    // don't do anything if current refinement center is not active
    if (!active[n]) {
       continue;
    }

    if (current_rho_max >= rho_max_list[rho_max_index[n]]/RHOGF) {
      // loop over all provided densities and check if they are above threshold
      // (this way we can switch on multiple levels at once!)
      while (current_rho_max >= rho_max_list[rho_max_index[n]]/RHOGF) {
	if(num_levels[n] < Carpet::maxreflevels) {
	  rho_max_index[n]++;
	  num_levels[n]++;

	  // regrid by changing CarpetRegrid2 parameters
	  CCTK_VInfo(CCTK_THORNSTRING,"CarpetRegrid2 has been told to regrid refinement center #%d (rho exceeded %5.5E g/cm^3)!!!!", n+1, rho_max_list[rho_max_index[n]-1]);
	  CCTK_VInfo(CCTK_THORNSTRING,"New number of refinement levels = %d!", num_levels[n]);
	  if(rho_max_list[rho_max_index[n]] > 0.0e0) {
	    CCTK_VInfo(CCTK_THORNSTRING,"Next regrid at rho_max = %5.5E, for refinement_center = #%d",
		      rho_max_list[rho_max_index[n]], n);
	  } else {
	    CCTK_Info(CCTK_THORNSTRING,"No further regridding planned");
	  }
	    //       fprintf(stderr,"max rls: %d, cur rls: %d\n", Carpet::maxreflevels, num_levels[0]);
	} else {
	  break;
	}
      }
    } else {
      CCTK_VInfo(CCTK_THORNSTRING,
		"Sorry, no regridding for refinement_center = %d! Current rho max: %5.5E g/cm^3",n+1,current_rho_max*RHOGF);
      CCTK_VInfo(CCTK_THORNSTRING,"Next regrid at rho_max = %5.5E g/cm^3 for refinement_center = #%d",
		rho_max_list[rho_max_index[n]], n+1);

      //    fprintf(stderr,"max rls: %d, cur rls: %d\n", Carpet::maxreflevels, num_levels[0]);
    }

    // in case we are beyond a given time, we may want to swictch off finer levels
    if (switch_off_refinement && cctk_time > switch_off_refinement_after_time)
    {
      num_levels[n] = 1;

      // regrid by changing CarpetRegrid2 parameters
      CCTK_VInfo(CCTK_THORNSTRING,"CarpetRegrid2 has been told to SWITCH OFF REFINEMENT!!!!");
    }
  }
  

  return;
} 


