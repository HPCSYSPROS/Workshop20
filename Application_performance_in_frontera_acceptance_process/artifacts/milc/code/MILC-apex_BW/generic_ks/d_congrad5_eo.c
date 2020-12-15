/******* d_congrad5_eo.c - conjugate gradient for SU3/fermions ****/
/* MIMD version 7 */
/* Kogut-Susskind fermions  -- improved */
/* for general even-odd dslash, using dslash_eo(), no longlinks
   or fatlinks */

/* Jim Hetrick, Kari Rummukainen, Doug Toussaint */
/* 5/10/07   C. DeTar True residual stopping criterion now. */

/* This version looks at the initial vector every "niter" passes */
/* The source vector is in "src", and the initial guess and answer
   in "dest".  "resid" is the residual vector, and "cg_p" and "ttt" are
   working vectors for the conjugate gradient.
   niter = maximum number of iterations.
   rsqmin = desired rsq, quit when we reach rsq <= rsqmin*source_norm.

   reinitialize after niters iterations and try once more.
*/
#include "generic_ks_includes.h"	/* definitions files and prototypes */
#include "../include/fermion_links.h"
#include "../include/prefetch.h"
#define FETCH_UP 1

#ifdef CGTIME
static const char *prec_label[2] = {"F", "D"};
#endif

/*#define CG_DEBUG*/

#define LOOPEND
#include "../include/loopend.h"

/* The Fermilab relative residue */

static Real 
relative_residue(su3_vector *p, su3_vector *q, int parity)
{
  double residue, num, den;
  int i;
  site *s;
  
  residue = 0;
  FORSOMEPARITY(i,s,parity){
    num = (double)magsq_su3vec( &(p[i]) );
    den = (double)magsq_su3vec( &(q[i]) );
    residue += (den==0) ? 1.0 : (num/den);
  } END_LOOP

  g_doublesum(&residue);

  if(parity == EVENANDODD)
    return sqrt(residue/volume);
  else
    return sqrt(2*residue/volume);
}

static int
ks_congrad_eo_parity( su3_vector *t_src, su3_vector *t_dest, 
		      quark_invert_control *qic, Real mass,
		      eo_links_t *eo){
  register int i;
  register site *s;
  int iteration;	/* counter for iterations */
  Real a,b;           	/* Sugar's a,b */
#ifdef FEWSUMS
  double actual_rsq;      /* rsq from actual summation of resid */
  double c_tr,c_tt,tempsum[4];	/* Re<resid|ttt>, <ttt|ttt> */
#endif
  double rsq = 0,relrsq = 0; /* resid**2, rel resid*2 */
  double oldrsq,pkp;	/*last resid*2,pkp = cg_p.K.cg_p */
  Real msq_x4;	/* 4*mass*mass */
  double source_norm;	/* squared magnitude of source vector */
  int otherparity = 0; /* the other parity */
  int nrestart;  /* Restart counter */
  su3_vector *ttt, *cg_p, *resid;
  char myname[] = "ks_congrad_eo_parity";

  /* Unpack structure */
  int niter        = qic->max;      /* maximum number of iters per restart */
  int max_restarts = qic->nrestart; /* maximum restarts */
  Real rsqmin      = qic->resid * qic->resid;    /* desired squared residual - 
			 normalized as sqrt(r*r)/sqrt(src_e*src_e) */
  Real relrsqmin   = qic->relresid * qic->relresid; /* desired relative residual (FNAL)*/
  int parity       = qic->parity;   /* EVEN, ODD */

  int max_cg = max_restarts*niter; /* Maximum number of iterations */

  msq_x4 = 4.0*mass*mass;

  switch(parity){
  case(EVEN): otherparity=ODD; break;
  case(ODD):  otherparity=EVEN; break;
  }

  /* Allocate temporary variables */
  /* PAD may be used to avoid cache trashing */
#define PAD 0
  ttt = (su3_vector *) malloc((sites_on_node+PAD)*sizeof(su3_vector));
  cg_p = (su3_vector *) malloc((sites_on_node+PAD)*sizeof(su3_vector));
  resid = (su3_vector *) malloc((sites_on_node+PAD)*sizeof(su3_vector));

  if(ttt == NULL || cg_p == NULL || resid == NULL){
    printf("%s(%d): No room for temporaries\n",myname,this_node);
  }

  /* Source norm */
  source_norm = 0.0;
  FORSOMEPARITY(i,s,parity){
    source_norm += (double)magsq_su3vec( &t_src[i] );
  } END_LOOP
  g_doublesum( &source_norm );
#ifdef CG_DEBUG
  node0_printf("congrad: source_norm = %e\n", (double)source_norm);
#endif

  /* Start CG iterations */
  
  nrestart = 0;
  iteration = 0;
  qic->size_r = 0;
  qic->size_relr = 0;
  qic->final_iters   = 0;
  qic->final_restart = 0;
  qic->converged     = 1;

  while(1) {
    /* Check for completion */
    if( ( iteration % niter == 0 ) || 
	( ( rsqmin    <= 0 || rsqmin    > qic->size_r   ) &&
	  ( relrsqmin <= 0 || relrsqmin > qic->size_relr) ) ) 
      {
	
	/* (re)initialization process */
	
	/* Compute true residual and relative residual */

	/* ttt <-  (-1)*M_adjoint*M*dest
	   resid,cg_p <- src + ttt
	   rsq = |resid|^2
	   source_norm = |src|^2
	*/
        rsq = 0.0;
	dslash_eo_field(t_dest, ttt, otherparity, eo);
	dslash_eo_field(ttt, ttt, parity, eo);
	/* ttt  <- ttt - msq_x4*src	(msq = mass squared) */
	FORSOMEPARITY(i,s,parity){
	  if( i < loopend-FETCH_UP ){
	    prefetch_VVVV( &ttt[i+FETCH_UP], 
			   &t_dest[i+FETCH_UP],
			   &t_src[i+FETCH_UP],
			   &resid[i+FETCH_UP]);
	  }
	  scalar_mult_add_su3_vector( &ttt[i], &t_dest[i], -msq_x4, &ttt[i] );
	  add_su3_vector( &t_src[i], &ttt[i], &resid[i] );
	  /* remember ttt contains -M_adjoint*M*src */
	  cg_p[i] = resid[i];
	  rsq += (double)magsq_su3vec( &resid[i] );
	} END_LOOP
#ifdef FEWSUMS
	actual_rsq = rsq; /* not yet summed over nodes */
#endif
        g_doublesum( &rsq );

	if(relrsqmin > 0)
	  relrsq = relative_residue(resid, t_dest, parity);

	if(source_norm == 0){
	  qic->final_rsq = 0.;
	  qic->final_relrsq = 0.;
	} else {
	  qic->final_rsq    = (Real)rsq/source_norm;
	  qic->final_relrsq = (Real)relrsq;
	}

	iteration++ ;  /* iteration counts number of multiplications
			  by M_adjoint*M */
	total_iters++;

#ifdef CG_DEBUG
	if(this_node==0)printf("CONGRAD: (re)start %d rsq = %.10e relrsq %.10e\n",
			       nrestart, qic->final_rsq, qic->final_relrsq);
#endif
	/* Quit when true residual and true relative residual are within
	   tolerance or when we exhaust iterations or restarts */
	
	if( iteration >= max_cg || 
	    nrestart  >= max_restarts ||
	    ( ( rsqmin    <= 0 || rsqmin    > qic->final_rsq   ) &&
	      ( relrsqmin <= 0 || relrsqmin > qic->final_relrsq) ) ) break;
	
	nrestart++;
      }

    /*
      oldrsq <- rsq
      ttt <- (-1)*M_adjoint*M*cg_p
      pkp <- (-1)*cg_p.M_adjoint*M.cg_p
      a <- -rsq/pkp
      dest <- dest + a*cg_p
      resid <- resid + a*ttt
      rsq <- |resid|^2
      b <- rsq/oldrsq
      cg_p <- resid + b*cg_p
    */

#ifdef FEWSUMS
    oldrsq = actual_rsq;	/* not yet summed over nodes */
#else
    oldrsq = rsq;
#endif
    pkp = 0.0;
    /* sum of neighbors */
    
    dslash_eo_field( cg_p, ttt, otherparity, eo);
    dslash_eo_field( ttt, ttt, parity, eo);
    
    /* finish computation of M_adjoint*m*p and p*M_adjoint*m*Kp */
    /* ttt  <- ttt - msq_x4*cg_p	(msq = mass squared) */
    /* pkp  <- cg_p.(ttt - msq*cg_p) */
    pkp = 0.0;
#ifdef FEWSUMS
    c_tr=0.0; c_tt=0.0;
#endif
    FORSOMEPARITY(i,s,parity){
      if( i < loopend-FETCH_UP ){
	prefetch_VV( &ttt[i+FETCH_UP], &cg_p[i+FETCH_UP] );
      }
      scalar_mult_add_su3_vector( &ttt[i], &cg_p[i], -msq_x4,
				  &ttt[i] );
      pkp += (double)su3_rdot( &cg_p[i], &ttt[i] );
#ifdef FEWSUMS
      c_tr += (double)su3_rdot( &ttt[i], &resid[i] );
      c_tt += (double)su3_rdot( &ttt[i], &ttt[i] );
#endif
    } END_LOOP
#ifdef FEWSUMS
    /* finally sum oldrsq over nodes, also other sums */
    tempsum[0] = pkp; tempsum[1] = c_tr; 
    tempsum[2] = c_tt; tempsum[3] = oldrsq;
    g_vecdoublesum( tempsum, 4 );
    pkp = tempsum[0]; c_tr = tempsum[1]; 
    c_tt = tempsum[2]; oldrsq = tempsum[3];
#else
    g_doublesum( &pkp );
#endif
    iteration++;
    total_iters++;
    
    a = (Real) (-rsq/pkp);
    
    /* dest <- dest - a*cg_p */
    /* resid <- resid - a*ttt */
#ifdef FEWSUMS
    actual_rsq=0.0;
#else
    rsq=0.0;
#endif
    FORSOMEPARITY(i,s,parity){
      if( i < loopend-FETCH_UP ){
	prefetch_VVVV( &t_dest[i+FETCH_UP], 
		       &cg_p[i+FETCH_UP], 
		       &resid[i+FETCH_UP], 
		       &ttt[i+FETCH_UP] );
      }
      scalar_mult_add_su3_vector( &t_dest[i], &cg_p[i], a, &t_dest[i] );
      scalar_mult_add_su3_vector( &resid[i], &ttt[i], a, &resid[i]);
#ifdef FEWSUMS
      actual_rsq += (double)magsq_su3vec( &resid[i] );
#else
      rsq += (double)magsq_su3vec( &resid[i] );
#endif
    } END_LOOP
#ifdef FEWSUMS
    /**printf("XXX:  node %d\t%e\t%e\t%e\n",this_node,oldrsq,c_tr,c_tt);**/
    rsq = oldrsq + 2.0*a*c_tr + a*a*c_tt; /*TEST - should equal actual_rsq */
    /**c_tt = actual_rsq;**/ /* TEMP for test */
    /**g_doublesum(&c_tt);**/ /* TEMP true value for rsq */
    /**node0_printf("RSQTEST: %e\t%e\t%e\n",rsq,c_tt,rsq-c_tt);**/
#else
    g_doublesum(&rsq);
#endif	

    if(relrsqmin > 0)
      relrsq = relative_residue(resid, t_dest, parity);
    
    qic->size_r        = (Real)rsq/source_norm;
    qic->size_relr     = (Real)relrsq;
    qic->final_iters   = iteration;
    qic->final_restart = nrestart;
    qic->converged     = 1;

#ifdef CG_DEBUG
    if(mynode()==0){printf("iter=%d, rsq/src= %e, relrsq= %e, pkp=%e\n",
			   iteration,(double)qic->size_r,
			   (double)qic->size_relr,
			   (double)pkp);fflush(stdout);}
#endif
    
    b = (Real)rsq/oldrsq;
    /* cg_p  <- resid + b*cg_p */
    FORSOMEPARITY(i,s,parity){
      scalar_mult_add_su3_vector( &resid[i],
				  &cg_p[i] , b , &cg_p[i]);
    } END_LOOP
  }

  if(nrestart == max_restarts || iteration == max_cg){
    qic->converged = 0;
    node0_printf("ks_congrad_eo_parity: CG not converged after %d iterations, \n",
		 iteration);
    node0_printf("rsq. = %e wanted %e relrsq = %e wanted %e\n",
		 qic->final_rsq,rsqmin,qic->final_relrsq,relrsqmin);
    fflush(stdout);
  }

  cleanup_dslash_temps();

  free(ttt); free(cg_p); free(resid);
  return iteration;
}

/* API for field arguments */

int ks_congrad_field( su3_vector *src, su3_vector *dest, 
			 quark_invert_control *qic, Real mass,
			 eo_links_t *eo )
{
  int iters = 0;
  double dtimec;
  double nflop = 1187;
  int parity = qic->parity;

  if(parity==EVENANDODD)nflop *=2;

  dtimec = -dclock(); 

  if(parity == EVEN || parity == EVENANDODD){
    qic->parity = EVEN;
    iters += ks_congrad_eo_parity(src, dest, qic, mass, eo);
  }
  if(parity == ODD || parity == EVENANDODD){
    qic->parity = ODD;
    iters += ks_congrad_eo_parity(src, dest, qic, mass, eo);
  }

  qic->parity = parity;

  dtimec += dclock();
#ifdef CGTIME
  if(this_node==0){
    printf("CONGRAD5: time = %e (eo %s) masses = 1 iters = %d mflops = %e\n",
	   dtimec, prec_label[PRECISION-1], iters, 
	   (double)(nflop*volume*iters/(1.0e6*dtimec*numnodes())) );
    fflush(stdout);}
#endif

  report_status(qic);

  return iters;
}

/* Needed for consistency with d_congrad5_fn.c. */
int ks_congrad_field_cpu( su3_vector *src, su3_vector *dest, 
			  quark_invert_control *qic, Real mass,
			  eo_links_t *eo ){
  return ks_congrad_field(src, dest, qic, mass, eo);
}

/* New API for site arguments */

int ks_congrad_site( field_offset src, field_offset dest, 
			quark_invert_control *qic, Real mass,
			eo_links_t *eo )
{
  int i;
  site *s;
  int iters = 0;
  su3_vector *t_src, *t_dest;
  double dtimec;
  double nflop = 1187;
  int parity = qic->parity;

  if(parity==EVENANDODD)nflop *=2;

  dtimec = -dclock(); 

  /* Map src and dest from site to field of correct precision */
  
  t_src  = (su3_vector *)malloc(sizeof(su3_vector)*sites_on_node);
  t_dest = (su3_vector *)malloc(sizeof(su3_vector)*sites_on_node);
  if(t_src == NULL || t_dest == NULL){
    printf("ks_congrad_site(%d): No room for temporaries\n",this_node);
    terminate(1);
  }

  FORALLSITES(i,s){
    t_src[i]  = *((su3_vector *)F_PT(s,src) );
    t_dest[i] = *((su3_vector *)F_PT(s,dest));
  }

  if(parity == EVEN || parity == EVENANDODD){
    qic->parity = EVEN;
    iters += ks_congrad_eo_parity(t_src, t_dest, qic, mass, eo );
  }
  else if(parity == ODD || parity == EVENANDODD){
    qic->parity = ODD;
    iters += ks_congrad_eo_parity(t_src, t_dest, qic, mass, eo );
  }

  qic->parity = parity;
    
  /* Map solution to site structure */

  FORALLSITES(i,s){
    *((su3_vector *)F_PT(s,dest)) = t_dest[i];
  }

  free(t_src); free(t_dest);

  dtimec += dclock();
#ifdef CGTIME
  if(this_node==0){
    printf("CONGRAD5: time = %e (eo %s) masses = 1 iters = %d mflops = %e\n",
	   dtimec, prec_label[PRECISION-1], iters, 
	   (double)(nflop*volume*iters/(1.0e6*dtimec*numnodes())) );
    fflush(stdout);}
#endif

  report_status(qic);

  return iters;
}

/* Traditional MILC API for site arguments and no relative residual test */

int ks_congrad( field_offset src, field_offset dest, Real mass,
		   int niter, int nrestart, Real rsqmin, int prec,
		   int parity, Real *final_rsq,
		   eo_links_t *eo )
{
  int iters;
  quark_invert_control qic;

  /* Pack structure */
  qic.prec      = prec;  /* Currently ignored */
  qic.min       = 0;
  qic.max       = niter;
  qic.nrestart  = nrestart;
  qic.start_flag = 0;
  qic.nsrc      = 1;
  qic.parity    = parity;
  qic.resid     = sqrt(rsqmin);
  qic.relresid  = 0;     /* Suppresses this test */

  /* Solve the system */
  iters = ks_congrad_site( src, dest, &qic, mass, eo );

  /* Unpack the results */
  *final_rsq    = qic.final_rsq;

  return iters;
}
