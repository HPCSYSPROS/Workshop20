/************************* control.c *******************************/
/* MIMD version 7 */
/* Main procedure for SU3 with dynamical staggered fermions        */
/* general quark action, general gauge action */

/* This file is for lattice generation with the RHMC algorithm */

#define CONTROL
#include "ks_imp_includes.h"	/* definitions files and prototypes */
#include "lattice_qdp.h"
#ifdef HAVE_QUDA
#include <quda_milc_interface.h>
#endif

#ifdef MILC_GLOBAL_DEBUG
#include "debug.h"
#endif /* MILC_GLOBAL_DEBUG */

/* For information */
#define NULL_FP -1

EXTERN gauge_header start_lat_hdr;	/* Input gauge field header */
int init_flag = 0;			/* let other functions know the stanza */

#if defined (NERSC_TIME) || defined (PROF_MARK)
#include <sys/time.h>
#include <sys/resource.h>
struct rusage usage;
double t_total, t_ks_ratinv, t_eo_fermion, t_restore_fermion; /* timing for update_h_fermion() */
#endif

int
main( int argc, char **argv )
{
  int i,meascount,traj_done, naik_index;
  int prompt;
  int s_iters, avs_iters, avbcorr_iters;
  double dtime, dclock();
  initialize_machine(&argc,&argv);

  /* Remap standard I/O */
  if(remap_stdio_from_args(argc, argv) == 1)terminate(1);
  
  g_sync();
  /* set up */
  prompt = setup();

  /* loop over input sets */
  while( readin(prompt) == 0) {
    
    /* perform warmup trajectories */
#ifdef MILC_GLOBAL_DEBUG
    global_current_time_step = 0;
#endif /* MILC_GLOBAL_DEBUG */

#if defined (NERSC_TIME) || defined (PROF_MARK)
    t_ks_ratinv = t_eo_fermion = t_restore_fermion = 0.0;
    t_total = -dclock();
#endif

    dtime = -dclock();


    for( traj_done=0; traj_done < warms; traj_done++ ){
      update();
    }
    node0_printf("WARMUPS COMPLETED\n"); fflush(stdout);
    
    /* perform measuring trajectories, reunitarizing and measuring 	*/
    meascount=0;		/* number of measurements 		*/
    avs_iters = avbcorr_iters = 0;

    for( traj_done=0; traj_done < trajecs; traj_done++ ){ 
#ifdef MILC_GLOBAL_DEBUG
#ifdef HISQ_REUNITARIZATION_DEBUG
  {
  int isite, idir;
  site *s;
  FORALLSITES(isite,s) {
    for( idir=XUP;idir<=TUP;idir++ ) {
      lattice[isite].on_step_Y[idir] = 0;
      lattice[isite].on_step_W[idir] = 0;
      lattice[isite].on_step_V[idir] = 0;
    }
  }
  }
#endif /* HISQ_REUNITARIZATION_DEBUG */
#endif /* MILC_GLOBAL_DEBUG */
      /* do the trajectories */
      s_iters=update();

      /* measure every "propinterval" trajectories */
      if( (traj_done%propinterval)==(propinterval-1) ){
	
	/* call gauge_variable fermion_variable measuring routines */
	/* results are printed in output file */
	rephase(OFF);
	g_measure( );
	rephase(ON);
#ifdef MILC_GLOBAL_DEBUG
#if ( FERM_ACTION == HISQ || FERM_ACTION == HYPISQ )
        g_measure_plaq( );
#endif
#ifdef MEASURE_AND_TUNE_HISQ
        g_measure_tune( );
#endif /* MEASURE_AND_TUNE_HISQ */
#endif /* MILC_GLOBAL_DEBUG */


	/**************************************************************/
	/* Compute chiral condensate and related quantities           */
	
	/* Make fermion links if not already done */
	
	restore_fermion_links_from_site(fn_links, par_buf.prec_pbp);
	for(i = 0; i < par_buf.num_pbp_masses; i++){
#if ( FERM_ACTION == HISQ || FERM_ACTION == HYPISQ )
	  naik_index = par_buf.ksp_pbp[i].naik_term_epsilon_index;
#else
	  naik_index = 0;
#endif
 	  f_meas_imp_field( par_buf.npbp_reps, &par_buf.qic_pbp[i], 
 			    par_buf.ksp_pbp[i].mass, naik_index, fn_links);

#ifdef D_CHEM_POT
	  Deriv_O6_field( par_buf.npbp_reps, &par_buf.qic_pbp[i],
			  par_buf.ksp_pbp[i].mass, naik_index, fn_links);
#endif
	}

	avs_iters += s_iters;
	++meascount;
	fflush(stdout);
      }
    }	/* end loop over trajectories */
    
    dtime += dclock();

    node0_printf("RUNNING COMPLETED\n"); fflush(stdout);
    if(meascount>0)  {
      node0_printf("average cg iters for step= %.3f\n",
		   (double)avs_iters/meascount);
    }

    node0_printf("total_iters = %d\n",total_iters);
    if( init_flag==0 ) {
      node0_printf("INIT_TIME = %.3f seconds\n",dtime);
      init_flag = 1;
    } else {
#ifdef NERSC_TIME
      node0_printf("TOTAL_TIME %.3f secs\n", dtime);
      node0_printf("CG_TIME    %.3f secs\n", t_ks_ratinv);
      node0_printf("FORCE_TIME %.3f secs\n", t_eo_fermion);
      node0_printf("LINK_TIME  %.3f secs\n", t_restore_fermion);
      if (getrusage(RUSAGE_SELF, &usage) == 0) {
        node0_printf("Approximate memory usage = %.3f MB\n", (float)usage.ru_maxrss*numnodes()/1024.0);
      }
#else
      node0_printf("TIME = %.3f seconds\n",dtime);
#endif
    }      

#ifdef HISQ_SVD_COUNTER
    node0_printf("hisq_svd_counter = %d\n",hisq_svd_counter);
#endif
#ifdef HYPISQ_SVD_COUNTER
    node0_printf("hypisq_svd_counter = %d\n",hypisq_svd_counter);
#endif
      
#ifdef HISQ_FORCE_FILTER_COUNTER
    node0_printf("hisq_force_filter_counter = %d\n",hisq_force_filter_counter);
#endif
#ifdef HYPISQ_FORCE_FILTER_COUNTER
    node0_printf("hypisq_force_filter_counter = %d\n",hypisq_force_filter_counter);
#endif
    fflush(stdout);
    
    /* save lattice if requested */
    if( saveflag != FORGET ){
      rephase( OFF );
      save_lattice( saveflag, savefile, stringLFN );
      rephase( ON );
    }

    /* Destroy fermion links (created in readin() */

#if FERM_ACTION == HISQ
    destroy_fermion_links_hisq(fn_links);
#elif FERM_ACTION == HYPISQ
    destroy_fermion_links_hypisq(fn_links);
#else
    destroy_fermion_links(fn_links);
#endif
    fn_links = NULL;
  }

#ifdef HAVE_QUDA
  qudaFinalize();
#endif

  normal_exit(0);
  return 0;
}

