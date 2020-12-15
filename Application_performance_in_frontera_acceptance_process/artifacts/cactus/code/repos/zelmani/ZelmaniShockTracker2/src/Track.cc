#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"


#include <carpet.hh>
#include <dh.hh>
#include <gh.hh>
#include <vect.hh>

#include <fstream>
#include <iomanip>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>


#include "assert.h"
#include "math.h"
#include <carpet.hh>

#include "util_Table.h"

extern "C" { 
  void ZelmaniShockTracker2_Track(CCTK_ARGUMENTS);
}



void ZelmaniShockTracker2_Track(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  using namespace std;
  using namespace Carpet;

  if(!*dostuff) return;
  if(!*dotrack) return;


  //int imax = maxloc(shockpos,3);
  //  int imin = minloc(shockpos,3);
  CCTK_REAL maxrad = *shockmax;
  CCTK_REAL mydxyz = 0.0;

  int do_regrid = 0;

  *dochangeregrid = 0;

  for (int i=0;i<track_n_levels;i++) {

    int do_regrid_local = 0;

    if(track_levels[i] == -1) continue;
    if(track_levels[i] == 0) continue;
    if(track_levels[i] > maxreflevels-1) continue;
    if(num_levels[0] < track_levels[i]) continue;

    CCTK_REAL newrad = radius[track_levels[i]];
    CCTK_REAL newradz = radius[track_levels[i]];
    // get grid spacing from Carpet
    enter_level_mode(cctkGH, track_levels[i]);
    enter_singlemap_mode(cctkGH,0,CCTK_GF);
    {
      DECLARE_CCTK_ARGUMENTS;
      mydxyz = CCTK_DELTA_SPACE(0);
    }
    leave_singlemap_mode(cctkGH);
    leave_level_mode(cctkGH);


    CCTK_VInfo(CCTK_THORNSTRING,"rl %d, dx=%15.6E",track_levels[i],mydxyz);

    if(do_regrid_z) {
       if(maxrad > radius_z[track_levels[i]] && radius_z[track_levels[i]] != -1.0) {
        CCTK_VInfo(CCTK_THORNSTRING,"rl %d, r=%15.6E, is inside shock",track_levels[i],radius_z[track_levels[i]]);
      } else {
        // reflevel boundary is outside shock
        CCTK_VInfo(CCTK_THORNSTRING,"rl %d, r=%15.6E, is outside shock",track_levels[i],radius[track_levels[i]]);
        CCTK_VInfo(CCTK_THORNSTRING,"rl %d, rz=%15.6E, is outside shock",track_levels[i],radius_z[track_levels[i]]);
        CCTK_REAL dist = radius_z[track_levels[i]] - maxrad;
        CCTK_INT  idist = (CCTK_INT)ceil(dist/mydxyz);
  
  //      CCTK_REAL distz = radius_z[track_levels[i]] - maxrad;
  //      CCTK_INT  idistz = (CCTK_INT)ceil(distz/mydxyz);
  
        if(dist < track_level_out_dist[i] || idist < track_level_out_zones[i]) {
  	CCTK_VInfo(CCTK_THORNSTRING,"rl %d, dist to shock =%15.6E, %d zones to shock: too small!!!",
  		   track_levels[i],dist,idist);
  	CCTK_INT inewzones = CCTK_INT(ceil( (maxrad + 
  					     track_level_dist_fac*max(track_level_out_dist[i], 
  								      track_level_out_zones[i]*mydxyz))/mydxyz));
  //	CCTK_INT inewzonesz = CCTK_INT(ceil( (maxrad + 
  //					     track_level_dist_fac*max(track_level_out_dist[i], 
  //								      track_level_out_zones[i]*mydxyz))/mydxyz));
  	//	inewzones = CCTK_INT(ceil(inewzones*track_level_dist_fac));
  
  //	newrad = inewzones*mydxyz;
  //	CCTK_VInfo(CCTK_THORNSTRING,"rl %d, newrad: %15.6E, oldrad: %15.6E, newzones: %d, oldzones: %d!!!",
  //		   track_levels[i],newrad,radius[track_levels[i]],inewzones,(CCTK_INT)ceil(radius[track_levels[i]]/mydxyz));
  	newradz = inewzones*mydxyz;
  	CCTK_VInfo(CCTK_THORNSTRING,"rl %d, newrad: %15.6E, oldrad: %15.6E, newzones: %d, oldzones: %d!!!",
  		   track_levels[i],newrad,radius[track_levels[i]],inewzones,(CCTK_INT)ceil(radius[track_levels[i]]/mydxyz));
  	CCTK_VInfo(CCTK_THORNSTRING,"rl %d, newradz: %15.6E, oldradz: %15.6E, newzones: %d, oldzonesz: %d!!!",
  		   track_levels[i],newradz,radius_z[track_levels[i]],inewzones,(CCTK_INT)ceil(radius_z[track_levels[i]]/mydxyz));
  
  	do_regrid = 1;
  	do_regrid_local = 1;
        }
      }
    } else {
      if(maxrad > radius[track_levels[i]]) {
        CCTK_VInfo(CCTK_THORNSTRING,"rl %d, r=%15.6E, is inside shock",track_levels[i],radius[track_levels[i]]);
      } else {
        // reflevel boundary is outside shock
        CCTK_VInfo(CCTK_THORNSTRING,"rl %d, r=%15.6E, is outside shock",track_levels[i],radius[track_levels[i]]);
        CCTK_VInfo(CCTK_THORNSTRING,"rl %d, rz=%15.6E, is outside shock",track_levels[i],radius_z[track_levels[i]]);
        CCTK_REAL dist = radius[track_levels[i]] - maxrad;
        CCTK_INT  idist = (CCTK_INT)ceil(dist/mydxyz);
  
  //      CCTK_REAL distz = radius_z[track_levels[i]] - maxrad;
  //      CCTK_INT  idistz = (CCTK_INT)ceil(distz/mydxyz);
  
        if(dist < track_level_out_dist[i] || idist < track_level_out_zones[i]) {
  	CCTK_VInfo(CCTK_THORNSTRING,"rl %d, dist to shock =%15.6E, %d zones to shock: too small!!!",
  		   track_levels[i],dist,idist);
  	CCTK_INT inewzones = CCTK_INT(ceil( (maxrad + 
  					     track_level_dist_fac*max(track_level_out_dist[i], 
  								      track_level_out_zones[i]*mydxyz))/mydxyz));
  //	CCTK_INT inewzonesz = CCTK_INT(ceil( (maxrad + 
  //					     track_level_dist_fac*max(track_level_out_dist[i], 
  //								      track_level_out_zones[i]*mydxyz))/mydxyz));
  	//	inewzones = CCTK_INT(ceil(inewzones*track_level_dist_fac));
  
  //	newrad = inewzones*mydxyz;
  //	CCTK_VInfo(CCTK_THORNSTRING,"rl %d, newrad: %15.6E, oldrad: %15.6E, newzones: %d, oldzones: %d!!!",
  //		   track_levels[i],newrad,radius[track_levels[i]],inewzones,(CCTK_INT)ceil(radius[track_levels[i]]/mydxyz));
  	newradz = inewzones*mydxyz;
  	CCTK_VInfo(CCTK_THORNSTRING,"rl %d, newrad: %15.6E, oldrad: %15.6E, newzones: %d, oldzones: %d!!!",
  		   track_levels[i],newrad,radius[track_levels[i]],inewzones,(CCTK_INT)ceil(radius[track_levels[i]]/mydxyz));
  	CCTK_VInfo(CCTK_THORNSTRING,"rl %d, newradz: %15.6E, oldradz: %15.6E, newzones: %d, oldzonesz: %d!!!",
  		   track_levels[i],newradz,radius_z[track_levels[i]],inewzones,(CCTK_INT)ceil(radius_z[track_levels[i]]/mydxyz));
  
  	do_regrid = 1;
  	do_regrid_local = 1;
        }
      }
    }
    if(do_regrid_local) {
      // change CarpetRegrid2 radius for this level
      if(do_regrid_z) {
        if (track_levels[i] ==6) {
          radius_x[track_levels[i]] = 120.0;
          radius_y[track_levels[i]] = 120.0;
          radius_z[track_levels[i]] = 120.0;
        }
        else {
          newrad = 300.0;
          radius_x[track_levels[i]] = newrad;
          radius_y[track_levels[i]] = newrad;
          radius_z[track_levels[i]] = newradz;
        }
      }
      else {
        radius[track_levels[i]] = newradz;
      }
    }
    
  } // end loop over all levels we need to track


  if(do_regrid) {

    CCTK_INFO ("Regridding");

    int pt;
    const CCTK_INT *temp_regrid_every;
    temp_regrid_every = (const CCTK_INT*)CCTK_ParameterGet("regrid_every","CarpetRegrid2",&pt);
    assert(temp_regrid_every);
    CCTK_VInfo(CCTK_THORNSTRING,"saved_regrid_every = %d !!!",*temp_regrid_every);
    *saved_regrid_every = *temp_regrid_every;
    //    CCTK_VInfo(CCTK_THORNSTRING,"saved_regrid_every = %d !!!",*saved_regrid_every);

    char regridstring[5];
    sprintf(regridstring,"%d",zst_regrid_every);
    int ierr =
      CCTK_ParameterSet ("regrid_every",
                         "CarpetRegrid2",
                         regridstring);
    ierr +=
      CCTK_ParameterSet ("radius_rel_change_threshold_1",
                         "CarpetRegrid2",
                         "0.000");

    if (ierr) {
      CCTK_WARN (CCTK_WARN_ABORT, "Could not enable regridding");
    }

    *dochangeregrid = 1;
  }

}

void ZelmaniShockTracker_StopRegridding (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INFO ("Disabling regridding again");

  ostringstream re;
  re << *saved_regrid_every;

  if(*dochangeregrid) {
    int ierr =
      CCTK_ParameterSet ("regrid_every",
                         "CarpetRegrid2",
                         re.str().c_str());

    ierr +=
      CCTK_ParameterSet ("radius_rel_change_threshold_1",
                         "CarpetRegrid2",
                         "0.01");

    if (ierr) {
      CCTK_WARN (CCTK_WARN_ABORT, "Could not disable regridding");
    }

    *dochangeregrid = 0;
  }

}

