/****** update_h_rhmc.c  -- ******************/
/* MIMD version 7 */
/* updates momentum matrices for improved action  with RHMC algorith*/
/* D.T. & J.H., naik term    8/96
*  D.T., fat link fermion term 5/97
*  D.T. general quark action 1/98
*  D.T. two types of quarks 3/99
*  T.D. and A.H. improved gauge updating spliced in 5/97
*  D.T. first try at RHMC version 12/05
*  D.T. 3/07 Gang together multiple pseudofermion terms in update_h_fermion
*/

#include "ks_imp_includes.h"	/* definitions files and prototypes */

#ifdef PROF_MARK
#  include <stdlib.h>
#  include <limits.h>
   extern int init_flag;
   int        profile_start=0;
   int        profile_stop=INT_MAX;
   double     t_profile;
   // Include any profile marking headers here
#  ifdef PROF_VTUNE
#    include <ittnotify.h>
#  endif
#endif

int update_h_rhmc( Real eps, su3_vector **multi_x ){
  int iters;
#ifdef FN
  invalidate_fermion_links(fn_links);
#endif
  /*  node0_printf("update_h_rhmc:\n"); */
  /* gauge field force */
  rephase(OFF);
  imp_gauge_force(eps,F_OFFSET(mom));
  rephase(ON);
  /* fermionic force */
  
  iters = update_h_fermion( eps,  multi_x );
  return iters;
} /* update_h_rhmc */

// gauge and fermion force parts separately, for algorithms that use
// different time steps for them
void update_h_gauge( Real eps ){
  /* node0_printf("update_h_gauge:\n");*/
  /* gauge field force */
  rephase(OFF);
  imp_gauge_force(eps,F_OFFSET(mom));
  rephase(ON);
} /* update_h_gauge */

// fermion force update grouping pseudofermions with the same path coeffs
int update_h_fermion( Real eps, su3_vector **multi_x ){
  int iphi,jphi;
  Real final_rsq;
  int i,j,n;
  int order, tmporder;
  Real *residues,*allresidues;
  Real *roots;
  int iters = 0;
  imp_ferm_links_t **fn;

#if defined(NERSC_TIME) || defined (PROF_MARK)
  extern struct rusage usage;
  extern double t_total, t_ks_ratinv, t_eo_fermion, t_restore_fermion;
  double t_tmp;
#endif

#ifdef PROF_MARK          
  if ((init_flag == 1) && (profile_stop == INT_MAX)) {
    char *s;
    profile_start = total_iters;
    if ((s = getenv("PROF_NUM_ITERS")) != NULL) {
      profile_stop = profile_start + atoi(s);
      node0_printf("------> PROFILING STARTED: start iters = %d, stop iters >= %d iters\n",
        profile_start, profile_stop);
      // Put profile start marking here
      #ifdef PROF_VTUNE
        __SSC_MARK(0x111);                   // start SDE tracing
        if (this_node == 0) __itt_resume();  // start Vtune
      #endif
      t_profile = -dclock();
    }
  }
#endif

  /* Algorithm sketch: assemble multi_x with all |X> fields,
     then call force routine for each part (so far we have to parts:
     zero correction to Naik and non-zero correction to Naik */

  allresidues = (Real *)malloc(n_order_naik_total*sizeof(Real));

  // Group the fermion force calculation according to sets of like
  // path coefficients.
  tmporder = 0;
  iphi = 0;
#if ( FERM_ACTION == HISQ || FERM_ACTION == HYPISQ )
  n = fermion_links_get_n_naiks(fn_links);
#else
  n = 1;
#endif
  for( i=0; i<n; i++ ) {
    for( jphi=0; jphi<n_pseudo_naik[i]; jphi++ ) {
#ifdef NERSC_TIME
      t_tmp = -dclock();
#endif
      restore_fermion_links_from_site(fn_links, prec_md[iphi]);
#ifdef NERSC_TIME
      t_restore_fermion += t_tmp + dclock();
#endif
      fn = get_fm_links(fn_links);

      // Add the current pseudofermion to the current set
      order = rparam[iphi].MD.order;
      residues = rparam[iphi].MD.res;
      roots = rparam[iphi].MD.pole;

      // Compute ( M^\dagger M)^{-1} in xxx_even
      // Then compute M*xxx in temporary vector xxx_odd 
      /* See long comment at end of file */
	/* The diagonal term in M doesn't matter */
#ifdef NERSC_TIME
      t_tmp = -dclock();
#endif
      iters += ks_ratinv( F_OFFSET(phi[iphi]), multi_x+tmporder, roots, order, 
			  niter_md[iphi], rsqmin_md[iphi], prec_md[iphi], EVEN, 
			  &final_rsq, fn[i], 
			  i, rparam[iphi].naik_term_epsilon );
#ifdef NERSC_TIME
      t_ks_ratinv += t_tmp + dclock();
#endif

      for(j=0;j<order;j++){
	dslash_field( multi_x[tmporder+j], multi_x[tmporder+j],  ODD,
		      fn[i]);
	allresidues[tmporder+j] = residues[j+1];
	// remember that residues[0] is constant, no force contribution.
      }
    tmporder += order;
    iphi++;
    }
  }

#ifdef MILC_GLOBAL_DEBUG
  node0_printf("update_h_rhmc: MULTI_X ASSEMBLED\n");fflush(stdout);
  node0_printf("update_h_rhmc: n_distinct_Naik=%d\n",n);
  for(j=0;j<n;j++)
    node0_printf("update_h_rhmc: orders[%d]=%d\n",j,n_orders_naik[j]);
#if ( FERM_ACTION == HISQ || FERM_ACTION == HYPISQ )
  for(j=0;j<n;j++)
    node0_printf("update_h_rhmc: masses_Naik[%d]=%f\n",j,fn_links.hl.eps_naik[j]);
#endif
  fflush(stdout);
#endif /* MILC_GLOBAL_DEBUG */

#ifdef NERSC_TIME
  t_tmp = -dclock();
#endif
  restore_fermion_links_from_site(fn_links, prec_ff);
#ifdef NERSC_TIME
  t_restore_fermion += t_tmp + dclock();
  t_tmp = -dclock();
#endif
  eo_fermion_force_multi( eps, allresidues, multi_x,
			  n_order_naik_total, prec_ff, fn_links );
#ifdef NERSC_TIME
  t_eo_fermion += t_tmp + dclock();
  node0_printf("step iters=%d: t_cg=%e, t_force=%e, t_link=%e\n", iters, t_ks_ratinv, t_eo_fermion, t_restore_fermion);
#endif

  free(allresidues);

#ifdef PROF_MARK
# include <sys/time.h>
# include <sys/resource.h>
  if (total_iters >= profile_stop) {
    t_tmp = dclock();
    t_profile += t_tmp;
    t_total += t_tmp;
    // Put profile stop marking here
    #ifdef PROF_VTUNE
      __SSC_MARK(0x222);                   // stop SDE tracing
      if (this_node == 0) __itt_pause();   // stop Vtune
    #endif
    node0_printf("------> PROFILING STOPPED: stop iters = %d profiled iters = %d\n",
      total_iters, total_iters - profile_start);
    node0_printf("TOTAL_TIME    %.3f secs\n", t_total);
    node0_printf("PROFILED_TIME %.3f secs\n", t_profile);
    node0_printf("CG_TIME       %.3f secs\n", t_ks_ratinv);
    node0_printf("FORCE_TIME    %.3f secs\n", t_eo_fermion);
    node0_printf("LINK_TIME     %.3f secs\n", t_restore_fermion);
    if (getrusage(RUSAGE_SELF, &usage) == 0) {
      node0_printf("Approximate memory usage = %.3f MiB\n", (float)usage.ru_maxrss*numnodes()/1024.0);
    }
    node0_printf("Exiting due to profiling iteration count exceeded\n");
    // call MILC exit routine and abort
    normal_exit(0); 
  }
#endif

  return iters;
} /* update_h_fermion */

