#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#include "carpet.hh"
#include "defs.hh"
#include "vect.hh"


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
  void ZelmaniShockTracker2_Print(CCTK_ARGUMENTS);
}

void ZelmaniShockTracker2_Print(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  using namespace std;
  using namespace Carpet;

  if(!*dostuff) return;
  // if(!*dotrack) return;

  if(verbose_level >= 1) {
    for(int i=0;i<num_levels[0];i++) {
      if(do_regrid_z) {
        CCTK_VInfo(CCTK_THORNSTRING,"RL radius_x/y[%d] = %8.3f",i,radius_x[i]);
        CCTK_VInfo(CCTK_THORNSTRING,"RL radius_z[%d] = %8.3f",i,radius_z[i]);
      } else {
        CCTK_VInfo(CCTK_THORNSTRING,"RL radius[%d] = %8.3f",i,radius[i]);
      }
    }
  }

  //    int pt; const CCTK_INT *temp_regrid_every; temp_regrid_every =
  //(const
  //CCTK_INT*)CCTK_ParameterGet("regrid_every","CarpetRegrid2",&pt);
  //assert(temp_regrid_every);
  //CCTK_VInfo(CCTK_THORNSTRING,"saved_regrid_every = %d
  //!!!",*temp_regrid_every);

}
