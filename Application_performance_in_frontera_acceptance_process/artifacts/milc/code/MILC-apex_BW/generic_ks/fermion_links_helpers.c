/**************** fermion_links_helpers.c *****************************/
/* MILC Version 7 */
/* Routines used by all of the fermion_links*.c files */

#include "generic_ks_includes.h"	/* definitions files and prototypes */
#define IMP_QUARK_ACTION_INFO_ONLY
#include <quark_action.h>

#define GOES_FORWARDS(dir) (dir<=TUP)
#define GOES_BACKWARDS(dir) (dir>TUP)

#ifdef QCDOC
#define special_alloc qcdoc_alloc
#define special_free qfree
#else
#define special_alloc malloc
#define special_free free
#endif

#ifdef ASQ_OPTIMIZED_FATTENING
static void 
compute_gen_staple_field(su3_matrix *staple, int mu, int nu, 
			 su3_matrix *link, su3_matrix *fatlink, Real coef);

static void 
compute_gen_staple_site(su3_matrix *staple, int mu, int nu, 
			field_offset link, su3_matrix* fatlink, Real coef);
#endif

void load_longlinks(ferm_links_t *fn, ks_action_paths *ap) {
  su3_matrix **t_ll = &fn->lng;
  register int i;
  register site *s;
  int ipath,dir;
  int disp[4];
  int num_q_paths = ap->num_q_paths;
  Q_path *q_paths = ap->q_paths;
  register su3_matrix *long1;
  su3_matrix *staple = NULL, *tempmat1 = NULL;
  char myname[] = "load_longlinks";

#ifdef LLTIME
  int nflop = 1804;
  double dtime;
  dtime=-dclock();
#endif

  if( phases_in != 1){
    node0_printf("BOTCH: %s needs phases in\n",myname);
    terminate(0);
  }
  /* Allocate space for t_ll if NULL */
  if(*t_ll == NULL){
    *t_ll = (su3_matrix *)special_alloc(sites_on_node*4*sizeof(su3_matrix));
    if(*t_ll==NULL){
      printf("%s(%d): no room for t_ll\n",myname, this_node);
      terminate(1);
    }
  }
  
  staple = (su3_matrix *)special_alloc(sites_on_node*sizeof(su3_matrix));
  if(staple == NULL){
    printf("%s(%d): Can't malloc temporary\n",myname,this_node);
    terminate(1);
  }

  tempmat1 = (su3_matrix *)special_alloc(sites_on_node*sizeof(su3_matrix));
  if(tempmat1 == NULL){
    printf("%s(%d): Can't malloc temporary\n",myname,this_node);
    terminate(1);
  }

  for (dir=XUP; dir<=TUP; dir++){ /* loop over longlink directions */
    /* set longlink to zero */
    FORALLSITES(i,s){
      long1 = *t_ll + 4*i +dir;
      clear_su3mat( long1 );
    }

    /* loop over paths, checking for ones with total displacement 3*dir */
    for( ipath=0; ipath<num_q_paths; ipath++ ){  /* loop over paths */
	/* compute total displacement of path */
	for(i=XUP;i<=TUP;i++)disp[i]=0;
	for( i=0; i<q_paths[ipath].length; i++){
	  if( GOES_FORWARDS(q_paths[ipath].dir[i]) )
	    disp[        q_paths[ipath].dir[i]  ]++;
	  else
	    disp[OPP_DIR(q_paths[ipath].dir[i]) ]--;
	}
	for( disp[dir]+=3,i=XUP; i<=TUP; i++)if(disp[i]!=0)break;
	if( i<=TUP )continue;  /* skip if path doesn't go to right place */
/**printf("ipath = %d, found a path:  ",ipath);
for(j=0;j<q_paths[ipath].length;j++)printf("\t%d", q_paths[ipath].dir[j]);
printf("\n");**/

	path_product( q_paths[ipath].dir, q_paths[ipath].length, tempmat1 );
	FORALLSITES(i,s){
	  su3_adjoint( &tempmat1[i], &staple[i] );
	  long1 = *t_ll + 4*i + dir;
          scalar_mult_add_su3_matrix( long1,
	    &staple[i], -q_paths[ipath].coeff, long1 );
		/* minus sign in coeff. because we used backward path*/
	}
    } /* ipath */

  } /* loop over directions */


  special_free(staple); staple = NULL;
  special_free(tempmat1); tempmat1 = NULL;

#ifdef LLTIME
dtime += dclock();
node0_printf("LLTIME(long): time =  %e (Naik) mflops = %e\n",dtime,
	     (Real)nflop*volume/(1e6*dtime*numnodes()) );
#endif
}  /* load_longlinks() */

/*
 * 10/01/02, flopcount for ASQ_OPTIMIZED_FATTENING - C. DeTar
 * Fatlinks:       61632 for load_fatlinks
 */

/* KS phases and APBC must be in the links. See long comment at 
   end of fermion_force_general.c */
void load_fatlinks(ferm_links_t *fn, ks_action_paths *ap){
  su3_matrix **t_fl = &fn->fat;
  register int i;
  register site *s;
  int dir;
  register su3_matrix *fat1;
  su3_matrix *staple = NULL, *tempmat1 = NULL;
  char myname[] = "load_fatlinks";

#ifdef ASQ_OPTIMIZED_FATTENING
  int  nu,rho,sig ;
  Real one_link;
  Real *act_path_coeff = ap->act_path_coeff;
#ifdef LLTIME
  char method[] = "Asqtad opt";
#endif
#else
  int ipath;
  int disp[4];
  int num_q_paths = ap->num_q_paths;
  Q_path *q_paths = ap->q_paths;
#ifdef LLTIME
  char method[] = "FN nonopt";
#endif
#endif

#ifdef LLTIME
  int nflop = 61632;
double dtime;
dtime=-dclock();
#endif

  if( phases_in != 1){
    node0_printf("BOTCH: %s needs phases in\n",myname);
    terminate(1);
  }

  /* Allocate space for t_fl if NULL */
  if(*t_fl == NULL){
    *t_fl = (su3_matrix *)special_alloc(sites_on_node*4*sizeof(su3_matrix));
    if(*t_fl==NULL){
      printf("NODE %d: no room for t_fl\n",this_node);
      terminate(1);
    }
  }
  
  staple = (su3_matrix *)special_alloc(sites_on_node*sizeof(su3_matrix));
  if(staple == NULL){
    printf("%s: Can't malloc temporary\n",myname);
    terminate(1);
  }

  tempmat1 = (su3_matrix *)special_alloc(sites_on_node*sizeof(su3_matrix));
  if(tempmat1 == NULL){
    printf("%s: Can't malloc temporary\n",myname);
    terminate(1);
  }

#ifndef  ASQ_OPTIMIZED_FATTENING   /* general case code */
  for (dir=XUP; dir<=TUP; dir++){ /* loop over fatlink directions */
    /* set fatlink to zero */
    FORALLSITES(i,s){
      fat1 = (*t_fl) + 4*i + dir;
      clear_su3mat( fat1 );
    }
    
    /* loop over paths, checking for ones with total displacement 1*dir */
    for( ipath=0; ipath<num_q_paths; ipath++ ){  /* loop over paths */
	/* compute total displacement of path */
	for(i=XUP;i<=TUP;i++)disp[i]=0;
	for( i=0; i<q_paths[ipath].length; i++){
	  if( GOES_FORWARDS(q_paths[ipath].dir[i]) )
	    disp[        q_paths[ipath].dir[i]  ]++;
	  else
	    disp[OPP_DIR(q_paths[ipath].dir[i]) ]--;
	}
	for( disp[dir]+=1,i=XUP; i<=TUP; i++)if(disp[i]!=0)break;
	if( i<=TUP )continue;  /* skip if path doesn't go to right place */
/**printf("dir = %d, found a path:  ",dir);
for(j=0;j<q_paths.[ipath].length;j++)printf("\t%d", q_paths[ipath].dir[j]);
printf("\n");**/

	path_product( q_paths[ipath].dir, q_paths[ipath].length, tempmat1 );
	FORALLSITES(i,s){
	  su3_adjoint( &tempmat1[i], &staple[i] );
	  fat1 = (*t_fl) +  4*i + dir;
          scalar_mult_add_su3_matrix( fat1,
	    &staple[i], -q_paths[ipath].coeff, fat1 );
		/* minus sign in coeff. because we used backward path*/
	}
    } /* ipath */
  } /* loop over directions */
#else	/* ASQ_OPTIMIZED_FATTENING, for Asq and Asqtad actions */
/*  Optimized fattening code for the Asq and Asqtad actions.           *
 * I assume that path 0 is the one link path 2 the 3-staple            *
 * path 3 the 5-staple path 4 the 7-staple and path 5 the Lapage term. *
 * Path 1 is the Naik term.                                            */
 
 /* to fix up the Lepage term, included by a trick below */
 one_link = (act_path_coeff[0] - 6.0*act_path_coeff[5]);
 
 for (dir=XUP; dir<=TUP; dir++){
   FORALLSITES(i,s) /* Intialize fat links with c_1*U_\mu(x) */
     {
       fat1 = (*t_fl) +  4*i + dir;
       scalar_mult_su3_matrix(&(s->link[dir]), one_link,
			      fat1 );
     }
   for(nu=XUP; nu<=TUP; nu++) if(nu!=dir)
     {
       compute_gen_staple_site(staple,dir,nu,F_OFFSET(link[dir]),
			       *t_fl, act_path_coeff[2]);
       /* The Lepage term */
       /* Note this also involves modifying c_1 (above) */
       compute_gen_staple_field(NULL,dir,nu,staple,
				*t_fl, act_path_coeff[5]);
       for(rho=XUP; rho<=TUP; rho++) if((rho!=dir)&&(rho!=nu))
	 {
	   compute_gen_staple_field( tempmat1, dir, rho, staple,
				     *t_fl, act_path_coeff[3]);
	   for(sig=XUP; sig<=TUP; sig++)
	     if((sig!=dir)&&(sig!=nu)&&(sig!=rho))
	       {
		 compute_gen_staple_field(NULL,dir,sig,tempmat1,
					  *t_fl, act_path_coeff[4]);
	       } /* sig */
	 } /* rho */
     } /* nu */
 }/* dir */  
#endif

 special_free(staple);  staple = NULL;
 special_free(tempmat1); tempmat1 = NULL;
#ifdef LLTIME
dtime += dclock();
 node0_printf("LLTIME(Fat): time = %e (%s) mflops = %e\n",dtime,method,
	      (Real)nflop*volume/(1e6*dtime*numnodes()) );
#endif
}  /* load_fatlinks() */


#ifdef  ASQ_OPTIMIZED_FATTENING   /* Asqtad action only, "_fn" executables */

#ifndef FN
BOMB THE COMPILE
#endif

static void 
compute_gen_staple_site(su3_matrix *staple, int mu, int nu, 
			field_offset link, su3_matrix* fatlink, Real coef) {
  su3_matrix tmat1,tmat2;
  msg_tag *mtag0,*mtag1;
  su3_matrix *tempmat = NULL;
  register site *s ;
  register int i ;
  register su3_matrix *fat1;

  /* Computes the staple :
                   mu
                +-------+
          nu	|	|
		|	|
		X	X
    Where the mu link can be any su3_matrix. The result is saved in staple.
    if staple==NULL then the result is not saved.
    It also adds the computed staple to the fatlink[mu] with weight coef.
  */

  /* Upper staple */
  mtag0 = start_gather_site( link, sizeof(su3_matrix), nu, EVENANDODD, gen_pt[0] );
  mtag1 = start_gather_site( F_OFFSET(link[nu]), sizeof(su3_matrix), mu, 
			EVENANDODD, gen_pt[1] );
  wait_gather(mtag0);
  wait_gather(mtag1);
  
  if(staple!=NULL){/* Save the staple */
    FORALLSITES(i,s){
      mult_su3_na( (su3_matrix *)gen_pt[0][i],
		   (su3_matrix *)gen_pt[1][i], &tmat1 );
      mult_su3_nn( &(s->link[nu]), &tmat1, &staple[i] );
    }
  }
  else{ /* No need to save the staple. Add it to the fatlinks */
    FORALLSITES(i,s){
      mult_su3_na( (su3_matrix *)gen_pt[0][i],
		   (su3_matrix *)gen_pt[1][i], &tmat1 );
      mult_su3_nn( &(s->link[nu]), &tmat1, &tmat2 );
      fat1 = &(fatlink[4*i+mu]);
      scalar_mult_add_su3_matrix(fat1, &tmat2, coef,
				 fat1) ;
    }
  }
  cleanup_gather(mtag0);
  cleanup_gather(mtag1);

  /* lower staple */
  tempmat = (su3_matrix *)special_alloc( sites_on_node*sizeof(su3_matrix) );
  mtag0 = start_gather_site( F_OFFSET(link[nu]),
			sizeof(su3_matrix), mu, EVENANDODD, gen_pt[0] );
  wait_gather(mtag0);
  FORALLSITES(i,s){
    mult_su3_an( &(s->link[nu]),(su3_matrix *)F_PT(s,link), &tmat1 );
    mult_su3_nn( &(tmat1),(su3_matrix *)gen_pt[0][i], &(tempmat[i]) );
  }
  cleanup_gather(mtag0);
  mtag0 = start_gather_field( tempmat, sizeof(su3_matrix),
				  OPP_DIR(nu), EVENANDODD, gen_pt[0] );
  wait_gather(mtag0);

  if(staple!=NULL){/* Save the staple */
    FORALLSITES(i,s){
      add_su3_matrix( &staple[i],(su3_matrix *)gen_pt[0][i], 
		      &staple[i] );
      fat1 = &(fatlink[4*i+mu]);
      scalar_mult_add_su3_matrix( fat1,
				 &staple[i], coef, 
				 fat1 );
    }
  }
  else{ /* No need to save the staple. Add it to the fatlinks */
    FORALLSITES(i,s){
      fat1 = &(fatlink[4*i+mu]);
      scalar_mult_add_su3_matrix( fat1,
				 (su3_matrix *)gen_pt[0][i], coef, 
				 fat1 );
    }
  }

  special_free(tempmat); tempmat = NULL;
  cleanup_gather(mtag0);
} /* compute_gen_staple_site */

/* Asqtad action only, "_fn" executables */

#ifndef FN
BOMB THE COMPILE
#endif

static void 
compute_gen_staple_field(su3_matrix *staple, int mu, int nu, 
			 su3_matrix *link, su3_matrix *fatlink, Real coef) {
  su3_matrix tmat1,tmat2;
  msg_tag *mtag0,*mtag1;
  su3_matrix *tempmat = NULL;
  register site *s ;
  register int i ;
  register su3_matrix *fat1;

  /* Computes the staple :
                   mu
                +-------+
          nu	|	|
		|	|
		X	X
    Where the mu link can be any su3_matrix. The result is saved in staple.
    if staple==NULL then the result is not saved.
    It also adds the computed staple to fatlink[mu] with weight coef.
  */

  /* Upper staple */
  mtag0 = start_gather_field( link, sizeof(su3_matrix), nu, EVENANDODD, gen_pt[0] );
  mtag1 = start_gather_site( F_OFFSET(link[nu]), sizeof(su3_matrix), mu, 
			EVENANDODD, gen_pt[1] );
  wait_gather(mtag0);
  wait_gather(mtag1);
  
  if(staple!=NULL){/* Save the staple */
    FORALLSITES(i,s){
      mult_su3_na( (su3_matrix *)gen_pt[0][i],
		   (su3_matrix *)gen_pt[1][i], &tmat1 );
      mult_su3_nn( &(s->link[nu]), &tmat1, &staple[i] );
    }
  }
  else{ /* No need to save the staple. Add it to the fatlinks */
    FORALLSITES(i,s){
      mult_su3_na( (su3_matrix *)gen_pt[0][i],
		   (su3_matrix *)gen_pt[1][i], &tmat1 );
      mult_su3_nn( &(s->link[nu]), &tmat1, &tmat2 );
      fat1 = &(fatlink[4*i+mu]);
      scalar_mult_add_su3_matrix(fat1, &tmat2, coef,
				 fat1) ;
    }
  }
  cleanup_gather(mtag0);
  cleanup_gather(mtag1);

  /* lower staple */
  tempmat = (su3_matrix *)special_alloc( sites_on_node*sizeof(su3_matrix) );
  mtag0 = start_gather_site( F_OFFSET(link[nu]),
			sizeof(su3_matrix), mu, EVENANDODD, gen_pt[0] );
  wait_gather(mtag0);
  FORALLSITES(i,s){
    mult_su3_an( &(s->link[nu]),&link[i], &tmat1 );
    mult_su3_nn( &(tmat1),(su3_matrix *)gen_pt[0][i], &(tempmat[i]) );
  }
  cleanup_gather(mtag0);
  mtag0 = start_gather_field( tempmat, sizeof(su3_matrix),
				  OPP_DIR(nu), EVENANDODD, gen_pt[0] );
  wait_gather(mtag0);

  if(staple!=NULL){/* Save the staple */
    FORALLSITES(i,s){
      add_su3_matrix( &staple[i],(su3_matrix *)gen_pt[0][i], 
		      &staple[i] );
      fat1 = &(fatlink[4*i+mu]);
      scalar_mult_add_su3_matrix( fat1,
				 &staple[i], coef, 
				 fat1 );
    }
  }
  else{ /* No need to save the staple. Add it to the fatlinks */
    FORALLSITES(i,s){
      fat1 = &(fatlink[4*i+mu]);
      scalar_mult_add_su3_matrix( fat1,
				 (su3_matrix *)gen_pt[0][i], coef, 
				 fat1 );
    }
  }

  special_free(tempmat); tempmat = NULL;
  cleanup_gather(mtag0);
} /* compute_gen_staple_field */
#endif  /* ASQ_OPTIMIZED_FATTENING   */

/* Move up the backward longlinks.  Result in t_lbl */
void 
load_longbacklinks(ferm_links_t *fn){
  su3_matrix **t_lbl = &fn->lngback;
  su3_matrix *t_ll = fn->lng;
  register int i;
  register site *s;
  int dir;
  su3_matrix *tempmat1 = NULL;
  msg_tag *tag[4];
  char myname[] = "load_longbacklinks";

  /* Allocate space for t_lbl if NULL */
  if(*t_lbl == NULL){
    *t_lbl = (su3_matrix *)special_alloc(sites_on_node*4*sizeof(su3_matrix));
    if(*t_lbl==NULL){
      printf("%s(%d): no room for t_lbl\n",myname,this_node);
      terminate(1);
    }
  }

  tempmat1 = (su3_matrix *)special_alloc(sites_on_node*sizeof(su3_matrix));
  if(tempmat1 == NULL){
    printf("%s: Can't malloc temporary\n",myname);
    terminate(1);
  }

  /* gather backwards longlinks */
  for( dir=XUP; dir<=TUP; dir ++){
    FORALLSITES(i,s){
      tempmat1[i] = t_ll[dir+4*i];
    }
    tag[dir] = start_gather_field( tempmat1,
      sizeof(su3_matrix), OPP_3_DIR(DIR3(dir)), EVENANDODD, gen_pt[dir] );
    wait_gather( tag[dir] );
    FORALLSITES(i,s){
      su3_adjoint( (su3_matrix *)gen_pt[dir][i],
      (*t_lbl) + dir + 4*i );
    }
    cleanup_gather( tag[dir] );
  }

  special_free(tempmat1); tempmat1 = NULL;
}

/* Move up the backward fatlinks.  Result in t_fbl */
void 
load_fatbacklinks(ferm_links_t *fn){
  su3_matrix **t_fbl = &fn->fatback;
  su3_matrix *t_fl = fn->fat;
  register int i;
  register site *s;
  int dir;
  su3_matrix *tempmat1 = NULL;
  msg_tag *tag[4];
  char myname[] = "load_fatbacklinks";

  /* Allocate space for t_fbl if NULL */
  if(*t_fbl == NULL){
    *t_fbl = (su3_matrix *)special_alloc(sites_on_node*4*sizeof(su3_matrix));
    if(*t_fbl==NULL){
      printf("%s(%d): no room for t_fbl\n",myname,this_node);
      terminate(1);
    }
  }

  tempmat1 = (su3_matrix *)special_alloc(sites_on_node*sizeof(su3_matrix));
  if(tempmat1 == NULL){
    printf("%s: Can't malloc temporary\n",myname);
    terminate(1);
  }

  /* gather backwards fatlinks */
  for( dir=XUP; dir<=TUP; dir ++){
    FORALLSITES(i,s){
      tempmat1[i] = t_fl[dir+4*i];
    }
    tag[dir] = start_gather_field( tempmat1,
      sizeof(su3_matrix), OPP_DIR(dir), EVENANDODD, gen_pt[dir] );
    wait_gather( tag[dir] );
    FORALLSITES(i,s){
      su3_adjoint( (su3_matrix *)gen_pt[dir][i],
      (*t_fbl) + dir + 4*i );
    }
    cleanup_gather( tag[dir] );
  }

  special_free(tempmat1); tempmat1 = NULL;
}

static void 
free_t_links(su3_matrix **t_l){
  if(*t_l != NULL) special_free(*t_l);
  *t_l = NULL;
}

/* Wrappers for MILC call to QOP */
void 
free_fn_links(ferm_links_t *fn){
  free_t_links(&fn->fat);
  free_t_links(&fn->lng);
#ifdef DBLSTORE_FN
  free_t_links(&fn->fatback);
  free_t_links(&fn->lngback);
#endif
  invalidate_all_ferm_links(fn);
}

#ifdef DM_DU0
/* Routines for dDslash/du0 */

void free_fn_links_dmdu0(ferm_links_t *fn){
  free_t_links(&fn->fat);
  invalidate_all_ferm_links(fn);
}
#endif

