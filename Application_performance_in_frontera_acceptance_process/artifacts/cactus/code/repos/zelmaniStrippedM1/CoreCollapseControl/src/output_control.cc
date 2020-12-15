#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "carpet.hh"
#define RHOGF 6.1755e17

void CoreCollapseControl_OutputControl(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  double inv_rho_gf = 1.0/RHOGF;

  if(force_postbounce && (!(*bounce) || *force_check)) {
    *in_prebounce = 0;
    *bounce = 1;
    *force_check = 0;
    *bouncetime = forced_bounce_time;
    CCTK_VInfo(CCTK_THORNSTRING,"CoreCollapseControl forcing bounce! bounce time: %10.2f",forced_bounce_time);
  }

  if(output_control) {
    // prebounce
    if( ( (*rho_max > prebounce_rho*inv_rho_gf && !(*in_prebounce))) 
	|| (force_postbounce && !(*in_prebounce))  ) {
      *in_prebounce = 1;

      // setting things to divisor
      int error=0;
      ostringstream outevery;
      outevery << "divisor";
      error = CCTK_ParameterSet("outScalar_criterion","CarpetIOScalar",
				outevery.str().c_str());
      assert(error==0);

      // change checking of corecollapse control 
      {
	int every=check_every_prebounce;
	outevery.str("");
	outevery << every;
	error = CCTK_ParameterSet("check_every","CoreCollapseControl",
				  outevery.str().c_str());
	assert(error==0);
	CCTK_VInfo(CCTK_THORNSTRING,"CoreCollapseControl check every %d time steps",every);
      }
      

      // setting output every every timestep
      if(preb_outscalar_every != -2) {
	int every=preb_outscalar_every;
	outevery.str("");
	outevery << every;
	
	error = CCTK_ParameterSet("outScalar_every","CarpetIOScalar",
				  outevery.str().c_str());
    
	//    fprintf(stderr,"%s %d \n",outevery.str().c_str(),error);
	assert(error==0);
	CCTK_VInfo(CCTK_THORNSTRING,"CarpetIOScalar scalar output every %d time steps",every);
      }  
      
      if(preb_out0D_every != -2) {
	int every=preb_out0D_every;
	outevery.str("");
	outevery << every;
	
	error = CCTK_ParameterSet("out0D_every","CarpetIOASCII",
				  outevery.str().c_str());

	assert(error==0);
	CCTK_VInfo(CCTK_THORNSTRING,"CarpetIOASCII 0D output every %d time steps",every);
      }
      
      if(preb_out2D_every != -2) {
	int every=preb_out2D_every;
	outevery.str("");
	outevery << every;
	
	error = CCTK_ParameterSet("out2D_every","CarpetIOASCII",
			      outevery.str().c_str());

	assert(error==0);
	CCTK_VInfo(CCTK_THORNSTRING,"CarpetIOASCII 2D output every %d time steps",every);
      }

      if(preb_out3Dhdf5_every != -2) {
	int every=preb_out3Dhdf5_every;
	outevery.str("");
	outevery << every;
	
	error = CCTK_ParameterSet("out_every","CarpetIOHDF5",
				  outevery.str().c_str());


	error += CCTK_ParameterSet("out3D_every","CarpetIOHDF5",
				  outevery.str().c_str());
	
	assert(error==0);
	CCTK_VInfo(CCTK_THORNSTRING,"CarpetIOHDF5 3D output every %d time steps",every);
      }

      if(preb_checkpoint_every != -2) {
	int every=preb_checkpoint_every;
	outevery.str("");
	outevery << every;
	
	error = CCTK_ParameterSet("checkpoint_every","IOUtil",
				  outevery.str().c_str());
	
	//      assert(error==0);
	if(error!=0) {
	  CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
		     "Failed to steer parameter io::checkpoint_every; tried to set it to %s, error was %d",
		     outevery.str().c_str(),error);
	} else {
	  CCTK_VInfo(CCTK_THORNSTRING,"Checkpoint output every %d time steps",every);
	}
      }

      // this part refers to many private thorns.
      if(preb_waves_every != -2) {
	int every=preb_waves_every;
	outevery.str("");
	outevery << every;
	
	if(CCTK_IsThornActive("InnerCore")) {
	  error = CCTK_ParameterSet("compute_every","InnerCore",
				    outevery.str().c_str());
     
	  if(error!=0) {
	    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
		       "Failed to steer parameter InnerCore::compute_every; tried to set it to %s, error was %d",
		       outevery.str().c_str(),error);
	  } else {
	    CCTK_VInfo(CCTK_THORNSTRING,"InnerCore output every %d time steps",every);
	  }
	}
	
	if(CCTK_IsThornActive("ZelmaniQuadWaveExtract")) {
	  error = CCTK_ParameterSet("compute_every","ZelmaniQuadWaveExtract",
				    outevery.str().c_str());
	  
	  if(error!=0) {
	    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
		       "Failed to steer parameter ZelmaniQuadWaveExtract::compute_every; tried to set it to %s, error was %d",
		       outevery.str().c_str(),error);
	  } else {
	    CCTK_VInfo(CCTK_THORNSTRING,"ZelmaniQuadWaveExtract output every %d time steps",every);
	  }
	}

	if(CCTK_IsThornActive("NPScalars")) {
	  error = CCTK_ParameterSet("NP_every","NPScalars",
				    outevery.str().c_str());
	  
	  if(error!=0) {
	    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
		       "Failed to steer parameter NPScalars::NP_every; tried to set it to %s, error was %d",
		       outevery.str().c_str(),error);
	  } else {
	    CCTK_VInfo(CCTK_THORNSTRING,"NPScalars output every %d time steps",every);
	  }
	}

	if(CCTK_IsThornActive("ZerilliIEF")) {
	  error = CCTK_ParameterSet("zerilli_every","ZerilliIEF",
				    outevery.str().c_str());
	  
	  if(error!=0) {
	    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
		       "Failed to steer parameter NPScalars::NP_every; tried to set it to %s, error was %d",
		       outevery.str().c_str(),error);
	  } else {
	    CCTK_VInfo(CCTK_THORNSTRING,"ZerilliIEF output every %d time steps",every);
	  }
	}
	
	if(CCTK_IsThornActive("WaveExtract")) {
	  error = CCTK_ParameterSet("out_every","WaveExtract",
				    outevery.str().c_str());
	  
	  if(error!=0) {
	    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
		       "Failed to steer parameter WaveExtract::out_every; tried to set it to %s, error was %d",
		       outevery.str().c_str(),error);
	  } else {
	    CCTK_VInfo(CCTK_THORNSTRING,"WaveExtract output every %d time steps",every);
	  }
	}
	

      }

    
      // CarpetIOASCII
      outevery.str("");
      outevery << "divisor";
      error = CCTK_ParameterSet("out1D_criterion","CarpetIOASCII",
				outevery.str().c_str());
      assert(!error);

      if(preb_out1D_every != -2) {
	// setting output every every timestep
	int every=preb_out1D_every;
	outevery.str("");
	outevery << every;
	
	error = CCTK_ParameterSet("out1D_every","CarpetIOASCII",
				  outevery.str().c_str());
	assert(!error);
	CCTK_VInfo(CCTK_THORNSTRING,"CarpetIOASCII 1D output every %d time steps",every);
      }
    }
  } // if(output_control)


  // bounce
  if(CCTK_EQUALS(bounce_criterion,"density")) 
    {
      if(*rho_max > bounce_rho*inv_rho_gf && !(*bounce)) {
	*bounce = 1;
	*bouncetime = cctk_time;
	CCTK_VInfo(CCTK_THORNSTRING,"max(rho) >= %15.6E g/cm^3 -- WE HAVE BOUNCE!!!",bounce_rho);
      }    
    } else {
    if(*global_entropy_max > bounce_entropy && !(*bounce)) {
	 *bounce = 1;
	 *bouncetime = cctk_time;
	 CCTK_VInfo(CCTK_THORNSTRING,"max(entropy in 100M) >= %5.2f -- WE HAVE BOUNCE!!!",bounce_entropy);
       }
     }

  if(*bounce && toggle_grhydro_eos_hot_eps_fix) {
      int error;
      error = CCTK_ParameterSet("GRHydro_eos_hot_eps_fix","GRHydro",
				"yes");
      assert(error==0);
  }

  if(*alp_min <= preBH_alpA && !(*in_preBH)) {
    if(preBH_force_cooling_off) {
      int error;
      error = CCTK_ParameterSet("force_off","ZelmaniHybridCool",
				"yes");
      assert(error==0);
      CCTK_VInfo(CCTK_THORNSTRING,"TURNING COOLING OFF!!!");
    }

    if(preBH_out3Dhdf5_every != -2) {
      int every=preBH_out3Dhdf5_every;
      int error;
      ostringstream outevery;
      outevery.str("");
      outevery << every;

      error = CCTK_ParameterSet("out_every","CarpetIOHDF5",
				outevery.str().c_str());

      assert(error==0);
      CCTK_VInfo(CCTK_THORNSTRING,"CarpetIOHDF5 3D output every %d time steps",every);
    }

  }


  if(*alp_min <= preBH_alpB && !(*in_preBH)) {
    *in_preBH = 1;
    if(preBH_AH_every != -2) {
      // setting AH finding every every timestep
      int every=preBH_AH_every;
      int error;
      ostringstream outevery;
      outevery.str("");
      outevery << every;
      
      error = CCTK_ParameterSet("find_every","AHFinderDirect",
				outevery.str().c_str());
      assert(!error);
      CCTK_VInfo(CCTK_THORNSTRING,"min(lapse) = %10.5f",*alp_min);
      CCTK_VInfo(CCTK_THORNSTRING,"AHFinderDirect::find_every to %d time steps",every);
    }


  } // preBH




  return;
}

void CoreCollapseControl_OutputControlRecover(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;
  // Force a new check after recovery
  *in_prebounce = 0;
  *in_preBH = 0;
  *force_check = 1;

  if(*bounce && toggle_grhydro_eos_hot_eps_fix) {
      int error;
      error = CCTK_ParameterSet("GRHydro_eos_hot_eps_fix","GRHydro",
				"yes");
      assert(error==0);
  }

  return;
}
