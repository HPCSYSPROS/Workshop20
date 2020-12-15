#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

void ZelmaniShockTracker2_CheckStatus(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  *dostuff = ( (((cctk_iteration) % compute_every) == 0 ) &&
               (cctk_time >= (start_time)) || (cctk_iteration)==0);

  //bouncetime comes from CoreCollapseControl
  *dotrack = *bounce && (cctk_time-*bouncetime >= start_time_tracking) && 
    ((((cctk_iteration) % check_tracking_every) == 0) && check_tracking_every > 0 );

  // bounce comes from CoreCollapseControl
  *dostuff = *dostuff && *bounce || *dotrack;

  //  CCTK_VInfo(CCTK_THORNSTRING,"%d, %d, %d, %d, %d",*dotrack,*dostuff, cctk_iteration, compute_every, check_tracking_every);
  // CCTK_VInfo(CCTK_THORNSTRING,"%10.5f %10.5f",cctk_time-*bouncetime, start_time_tracking);

  return;
}


void ZelmaniShockTracker2_ForceFind(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // bounce comes from CoreCollapseControl
  *dostuff =  *bounce;

  return;
}
