/****** eigen_stuff_Ritz.c  ******************/
/* Eigenvalue and Eigevector computation routines.
* K.O. 8/99 Started. 
* UH did some stuff to this I think
* EBG 2/15/02 changed phi -> phi1, mass -> mass1 so could be used 
*    w/ ks_imp_dyn2
* EBG 6/2004 fixed some memory leaks
* MIMD version 7
*
*  These routines are for the computation of the Eigenvalues and Eigevectors
* of the Kogut-Susskind dslash^2. 
*/
/**#define DEBUG**/
/* This Variable controls the tollerence of the JACOBI iteration */
#define JACOBI_TOL 1.110223e-16
#define MINITER 5

/*   If you use strict convergence then you stop when g = |M*v - l*v|     *
 *  becomes less than the eigenvalue tolerance. Other wise you stop       *
 *  when the estimated error of the eigenvalue gets smaller than          *
 *  the eigenvalue tolerance. The later option takes abvantage of         *
 *  the quadratic convergence of the algorithm. (See Kalkreuter's paper    * 
 *  for details). I prefer the strict convergenve which results some      *
 *  overkill.                                                             */
/**#define STRICT_CONVERGENCE	(taken out by UMH)**/


/* If you do not define the USE_DSLASH_SPECIAL then calls to the standard   *
 * dslash are made. If you do define USE_DSLASH_SPECIAL then DSLASH_SPECIAL *
 * is used                                                                  */
#ifdef FN
#define USE_DSLASH_SPECIAL
#endif


/* Include files */
#include "generic_ks_includes.h"
#include "../include/jacobi.h"
#include "../include/dslash_ks_redefine.h"
#include <string.h>

static void copy_Vector(su3_vector *src, su3_vector *res ) ; 
static void norm2(su3_vector *vec, double *norm, int parity); 
static void dot_product(su3_vector *vec1, su3_vector *vec2, 
		 double_complex *dot, int parity) ;
static void complex_vec_mult_sub(double_complex *cc, su3_vector *vec1, 
			  su3_vector *vec2, int parity) ;
static void complex_vec_mult_add(double_complex *cc, su3_vector *vec1, 
			  su3_vector *vec2, int parity) ;
static void double_vec_mult(double *a, su3_vector *vec1, 
		     su3_vector *vec2, int parity) ;
static void double_vec_mult_sub(double *rr, su3_vector *vec1,  
			 su3_vector *vec2, int parity) ;
static void dax_p_by(double *a, su3_vector *vec1, double *b, su3_vector *vec2, 
	      int parity) ; 
static void vec_plus_double_vec_mult(su3_vector *vec1, double *a, su3_vector *vec2, 
			      int parity) ;
static void normalize(su3_vector *vec,int parity) ;
static void project_out(su3_vector *vec, su3_vector **vector, int Num, int parity);
static void constructArray(su3_vector **eigVec, su3_vector **MeigVec, Matrix *A,
		    double *err, int *converged, int parity,
		    imp_ferm_links_t *fn);
static void RotateBasis(su3_vector **eigVec, Matrix *V, int parity) ;

/************************************************************************/

#if 0
/* Modified Gram-Schmidt orthonormalization pg. 219 Golub & Van Loan     *
 * Num is the number of vectors. vector are the vectors. They get        *
 * overwritten by the orthonormal vectors.                               *
 * parity is the parity on which we work on (EVEN,ODD,ENENANDODD).       */
static void GramSchmidt(su3_vector **vector, int Num, int parity)/* NEVER USED */{
  register int i,j ;
  double norm ;
  double_complex cc ;

  for(i=0;i<Num;i++){
    norm2(vector[i], &norm, parity);
    norm = 1.0/norm ;
    double_vec_mult(&norm, vector[i], vector[i], parity);
    for(j=i+1;j<Num;j++){
      dot_product(vector[i], vector[j], &cc, parity) ;
      complex_vec_mult_sub(&cc, vector[i], vector[j], parity);
    }
  }
}
#endif
/************************************************************************/

/*  Projects out the *vectors from the  vec. Num is the Number of vectors  *
 * and parity is the parity on which we work on.                           *
 * The vectors are assumed to be orthonormal.                              */
static void project_out(su3_vector *vec, su3_vector **vector, int Num, int parity){
  register int i ;
  double_complex cc ;
  for(i=Num-1;i>-1;i--){
    dot_product(vector[i], vec, &cc, parity) ;
    complex_vec_mult_sub(&cc, vector[i], vec, parity);
  }
}

/************************************************************************/
/* Copies scr to res */
static void copy_Vector(su3_vector *src, su3_vector *res){
  memcpy((void *)res, (void *)src, sites_on_node*sizeof(su3_vector)) ;
}

/************************************************************************/
/* Returns the 2-norm of a fermion vector */
static void norm2(su3_vector *vec, double *norm, int parity){
  register double n ;
  register site *s;
  register  int i;
  
  n=0 ; 
  FORSOMEPARITY(i,s,parity){
    n += magsq_su3vec(&(vec[i]));
  }
  *norm = n ;
  g_doublesum(norm);
  *norm = sqrt(*norm) ;
}
 
/*****************************************************************************/
/* Returns the dot product of two fermion vectors */
static void dot_product(su3_vector *vec1, su3_vector *vec2, 
		   double_complex *dot, int parity) {
  register double re,im ;
  register site *s;
  register  int i;
  complex cc ;
  
  re=im=0.0;
  FORSOMEPARITY(i,s,parity){
    cc = su3_dot( &(vec1[i]), &(vec2[i]) );
    re += cc.real ;
    im += cc.imag ;
  }
  dot->real = re ; 
  dot->imag = im ;
  g_dcomplexsum(dot);
}

/*****************************************************************************/
/* Returns vec2 = vec2 - cc*vec1   cc is a double complex   */
static void complex_vec_mult_sub(double_complex *cc, su3_vector *vec1, 
			  su3_vector *vec2, int parity){

  register site *s;
  register  int i;
  complex sc ;
  
  sc.real= (Real)(cc->real) ; 
  sc.imag= (Real)(cc->imag) ;

  FORSOMEPARITY(i,s,parity){
    c_scalar_mult_sub_su3vec(&(vec2[i]), (&sc), &(vec1[i])) ;
  }
}

/*****************************************************************************/
/* Returns vec2 = vec2 + cc*vec1   cc is a double complex   */
static void complex_vec_mult_add(double_complex *cc, su3_vector *vec1, 
			  su3_vector *vec2, int parity){
  register site *s;
  register  int i;
  complex sc ;
  
  sc.real= (Real)(cc->real) ; sc.imag= (Real)(cc->imag) ;
  FORSOMEPARITY(i,s,parity){
    c_scalar_mult_add_su3vec(&(vec2[i]), (&sc), &(vec1[i])) ;
  }
}

/*****************************************************************************/
/* Returns vec2 = vec2 - rr*vec1   rr is a double    */
static void double_vec_mult_sub(double *rr, su3_vector *vec1,  
			 su3_vector *vec2, int parity){
  register site *s;
  register  int i;
  
  FORSOMEPARITY(i,s,parity){
    scalar_mult_sub_su3_vector(&(vec2[i]), &(vec1[i]), (Real)*rr, &(vec2[i]));
  }
}

#if 0
/*****************************************************************************/
/* Returns vec2 = vec2 + rr*vec1   rr is a double    */
static void double_vec_mult_add(double *rr, su3_vector *vec1,  
			 su3_vector *vec2, int parity){
  register site *s;
  register  int i;
  
  FORSOMEPARITY(i,s,parity){
    scalar_mult_add_su3_vector(&(vec2[i]), &(vec1[i]), (Real)*rr, &(vec2[i]));
  }
}
#endif

/*****************************************************************************/
/* Returns vec2 = a*vec1   a is a double vec2 can be vec1*/
static void double_vec_mult(double *a, su3_vector *vec1, 
		     su3_vector *vec2, int parity){
  
  register site *s;
  register  int i;
  
  FORSOMEPARITY(i,s,parity){ 
    scalar_mult_su3_vector( &(vec1[i]),(Real)*a, &(vec2[i])) ;
  }
}

/*****************************************************************************/
/* Returns vec2 = vec1 + a*vec2 */
static void vec_plus_double_vec_mult(su3_vector *vec1, double *a, su3_vector *vec2, 
			      int parity){
  register site *s;
  register  int i;
  
  FORSOMEPARITY(i,s,parity){ 
    scalar_mult_su3_vector( &(vec2[i]), *a, &(vec2[i]) ) ;
    add_su3_vector( &(vec1[i]), &(vec2[i]), &(vec2[i]) ) ;
  }
}

/*****************************************************************************/
/* Returns vec1 = a*vec1 + b*vec2   a,b are double    */
static void dax_p_by(double *a, su3_vector *vec1, double *b, su3_vector *vec2, 
	      int parity) {

  register site *s;
  register  int i;
  
  FORSOMEPARITY(i,s,parity){
    scalar_mult_su3_vector(&(vec1[i]), (Real)(*a), &(vec1[i])) ;
    scalar_mult_add_su3_vector(&(vec1[i]), &(vec2[i]), (Real)(*b), 
			       &(vec1[i])) ;
  }
}

/*****************************************************************************/
/* normalizes the vecror vec. Work only on parity. */
static void normalize(su3_vector *vec,int parity){

  double norm ;
  norm2(vec,&norm,parity) ;
  norm = 1.0/norm ;
  double_vec_mult(&norm,vec,vec,parity) ;
}

/*****************************************************************************/
/* Run a CG iteration to minimize the Rayleigh-Ritz quotient over the
   span "eigVec".  The initial vector is "vec".  The CG iteration stops
   when the tolerance is reached or the norm of the residual drops
   faster than RelTol, or when the MaxIter and Restart are exhausted. */

int Rayleigh_min(su3_vector *vec, su3_vector **eigVec, Real Tolerance, 
		 Real RelTol, int Nvecs, int MaxIter, int Restart, 
		 int parity, imp_ferm_links_t *fn){

  int iter ;
  double beta, cos_theta, sin_theta ;
  double quot, old_quot, P_norm, theta, real_vecMp, pMp ;
#ifdef DEBUG
  double vec_norm ;
#endif
  double g_norm, old_g_norm, start_g_norm ;
  double_complex cc ;
  su3_vector *Mvec, *grad, *P, *MP ;
  
  Mvec     = (su3_vector *)malloc(sites_on_node*sizeof(su3_vector)) ;
  grad     = (su3_vector *)malloc(sites_on_node*sizeof(su3_vector)) ;
  P        = (su3_vector *)malloc(sites_on_node*sizeof(su3_vector)) ;
  MP       = (su3_vector *)malloc(sites_on_node*sizeof(su3_vector)) ;


  old_quot = 1.0e+16 ;
  project_out(vec, eigVec, Nvecs, parity);
  normalize(vec,parity) ; 
  Matrix_Vec_mult(vec,Mvec,parity, fn) ;
  project_out(Mvec, eigVec, Nvecs, parity);
  
  /* Compute the quotient quot=vev*M*vec */
  dot_product(vec, Mvec, &cc, parity) ;
  /* quot is real since M is hermitian. quot = vec*M*vec */
  quot = cc.real ;
#ifdef DEBUG
  node0_printf("Rayleigh_min: Start -- quot=%g,%g\n", quot,cc.imag) ;
#endif  
  /* Compute the grad=M*vec - quot*vec */
  copy_Vector(Mvec,grad) ;
  double_vec_mult_sub(&quot, vec, grad, parity) ;
  /* set P (the search direction) equal to grad */
  copy_Vector(grad,P) ;
  /* compute the norms of P and grad */
  norm2(P   , &P_norm, parity) ;
  g_norm = P_norm ;
  start_g_norm = g_norm ;
#ifdef DEBUG
  node0_printf("Rayleigh_min: Start -- g_norm=%g\n", g_norm) ;
#endif  

  iter = 0 ;
  while( (g_norm>Tolerance)&&
	 ( ((iter<MaxIter)&&(g_norm/start_g_norm>RelTol)) || (iter<MINITER) ) 
	 ){
    iter++ ;
    Matrix_Vec_mult(P,MP,parity, fn) ;
    dot_product(vec, MP, &cc, parity) ;
    real_vecMp = cc.real ;
    dot_product(P, MP, &cc, parity) ;
    pMp = cc.real ; /*pMp is real */
    theta = 0.5*atan(2.0*real_vecMp/(quot*P_norm - pMp/P_norm)) ;
    sin_theta = sin(theta) ;
    cos_theta = cos(theta) ;
    if(sin_theta*cos_theta*real_vecMp>0){
      theta = theta - 0.5*PI ; /* chose the minimum not the maximum */
      sin_theta = sin(theta) ; /* the sin,cos calls can be avoided */
      cos_theta = cos(theta) ;
    }
    sin_theta = sin_theta/P_norm ;
    /* vec = cos(theta)*vec +sin(theta)*P/p_norm */
    dax_p_by(&cos_theta, vec,&sin_theta, P,parity) ;
    /* Mvec = cos(theta)*Mvec +sin(theta)*MP/p_norm */
    dax_p_by(&cos_theta,Mvec,&sin_theta,MP,parity) ;
    /* renormalize vec ... */
    if(iter%Restart == 0 ) {
#ifdef DEBUG
      node0_printf("Renormalizing...");
      norm2(vec,&vec_norm,parity) ;
      node0_printf("  norm: %g\n",1.0/vec_norm);
#endif
      /* Project vec on the orthogonal complement of eigVec */
      project_out(vec, eigVec, Nvecs, parity);
      normalize(vec,parity) ;
      Matrix_Vec_mult(vec,Mvec,parity, fn);
      /* Recompute the quotient */
      dot_product(vec, Mvec, &cc, parity) ;
      /* quot is real since M is hermitian. quot = vec*M*vec */
      quot = cc.real ;
      /* Recompute the grad */
      copy_Vector(Mvec,grad) ;
      double_vec_mult_sub(&quot, vec, grad, parity) ;
      norm2(grad,&g_norm,parity) ; /* recompute the g_norm */
      /* Project P on the orthogonal complement of eigVec */
      project_out(P, eigVec, Nvecs, parity);
      /* make P orthogonal to vec */
      dot_product(vec, P, &cc, parity) ;
      complex_vec_mult_sub(&cc, vec, P, parity);
      /* make P orthogonal to grad */
      dot_product(grad, P, &cc, parity) ;
      complex_vec_mult_sub(&cc, grad, P, parity);
      norm2(P,&P_norm,parity) ; /* recompute the P_norm */
    }
    dot_product(vec, Mvec, &cc, parity) ;
    /* quot is real since M is hermitian. quot = vec*M*vec */
    quot = cc.real ;
#ifdef DEBUG
    node0_printf("Rayleigh_min: %i, quot=%8g g=%8g b=%6g P:%6g\n",
		 iter,cc.real,g_norm,beta,P_norm) ;
#endif      
    old_g_norm = g_norm ;
       
    copy_Vector(Mvec,grad) ;
    double_vec_mult_sub(&quot, vec, grad, parity) ;

    norm2(grad,&g_norm,parity) ;
    beta = cos_theta*g_norm*g_norm/(old_g_norm*old_g_norm) ;
    /* Cut off beta */
    if( beta>2.0 ) {
      beta = 2.0 ;
    }
    dot_product(vec, P, &cc, parity) ;
    cc.real *= beta ; cc.imag *= beta ;
    /*      P_norm = 1.0/P_norm ;
	    double_vec_mult(&P_norm,  P,  P, parity) ; */
    vec_plus_double_vec_mult(grad, &beta, P, parity) ;/* P = grad + beta*P */
    complex_vec_mult_sub(&cc, vec, P, parity) ; /* P = P - cc*vec */
    norm2(P   , &P_norm, parity) ;

    if(fabs(old_quot -quot)< Tolerance/100.0){
      /*break*/
      g_norm = Tolerance/10.0; 
    }

    old_quot = quot ;
  }

  project_out(vec, eigVec, Nvecs, parity);
  normalize(vec,parity) ;
  free(MP) ;
  free(P) ;
  free(grad) ;
  free(Mvec) ;
  iter++ ;
  return iter ;
}

/*****************************************************************************/
/* Returns the projected matrix A and the error of each eigenvector */
static void constructArray(su3_vector **eigVec, su3_vector **MeigVec, Matrix *A,
		    double *err, int *converged, int parity,
		    imp_ferm_links_t *fn){
  int i,j,Nvecs ;
  su3_vector *grad ;
  double_complex cc,Aij,Aji ;

  Nvecs = A->N ;
  grad = (su3_vector *)malloc(sites_on_node*sizeof(su3_vector));
  for(i=0;i<Nvecs;i++){
    if(converged[i]==0) Matrix_Vec_mult(eigVec[i],MeigVec[i],parity, fn) ;
    dot_product(MeigVec[i], eigVec[i], &cc, parity) ;
    A->M[i][i].real = cc.real ; A->M[i][i].imag = 0.0 ;
    copy_Vector(MeigVec[i], grad) ;
    double_vec_mult_sub(&cc.real, eigVec[i], grad, parity) ;
    norm2(grad, &err[i], parity) ;
    for(j=i+1;j<Nvecs;j++){
      dot_product(MeigVec[i], eigVec[j], &cc, parity) ;
      Aij=cc ;
      CONJG(cc, Aji) ;
      dot_product(eigVec[j], MeigVec[i], &cc, parity) ;
      CSUM(Aji,cc) ;
      CONJG(cc, cc) ;
      CSUM(Aij,cc) ;
      CMULREAL(Aij,0.5,A->M[i][j]) ;
      CMULREAL(Aji,0.5,A->M[j][i]) ;
    }
  }
  free(grad) ;
}


/*****************************************************************************/
static void RotateBasis(su3_vector **eigVec, Matrix *V, int parity){

  su3_vector **Tmp ;
  register site *s;
  register  int i,N ;
  int j,k ;
  
  N = V->N ;
  /* Allocate the temporary vectors needed */
  Tmp = (su3_vector **)malloc(N*sizeof(su3_vector *));
  for(j=0;j<N;j++){
    Tmp[j]=(su3_vector *)malloc(sites_on_node*sizeof(su3_vector));
  }

  for(j=0;j<N;j++){
    FORSOMEPARITY(i,s,parity){
      clearvec(&Tmp[j][i]) ;
    }
    for(k=0;k<N;k++){
      complex_vec_mult_add(&V->M[k][j],eigVec[k],Tmp[j],parity) ;
    }
  }
  
  /* Copy rotated basis to the eigVec and free temporaries */
  for(j=0;j<N;j++) {
    /* Copy only wanted parity (UMH) */
    /** copy_Vector(Tmp[j], eigVec[j]) ; **/
    FORSOMEPARITY(i,s,parity){
      su3vec_copy( &(Tmp[j][i]), &(eigVec[j][i]));
    }
    free(Tmp[j]) ;
  }
  free(Tmp) ;
}

/*****************************************************************************/
int Kalkreuter_Ritz(su3_vector **eigVec, double *eigVal, Real Tolerance, 
	       Real RelTol, int Nvecs, int MaxIter, 
	       int Restart, int Kiters ){

  int total_iters=0 ;
  int j;
  Matrix Array,V ;
  register site *s ;
  register  int i ;
  su3_vector *vec ;
  su3_vector **MeigVec ;
  double max_error = 1.0e+10 ;
  double *grad, *err ;
  int iter = 0 ;
  int *converged ;
  int parity = active_parity ;
  Real ToleranceG ;
#ifdef EIGTIME
  double dtimec;
#endif
  imp_ferm_links_t **fn;

  fn = get_fm_links(fn_links);


  ToleranceG = 10.0*Tolerance ;

  /** Allocate the array **/
  Array = AllocateMatrix(Nvecs) ;  
  /** Allocate the Eigenvector matrix **/
  V = AllocateMatrix(Nvecs) ;  

  vec = (su3_vector *)malloc(sites_on_node*sizeof(su3_vector));
  grad = (double *) malloc(Nvecs*sizeof(double)) ;
  err = (double *) malloc(Nvecs*sizeof(double)) ;
  converged = (int *)malloc(Nvecs*sizeof(int));
  MeigVec = (su3_vector **)malloc(Nvecs*sizeof(su3_vector*));
  for(i=0;i<Nvecs;i++){
    MeigVec[i] = (su3_vector *)malloc(sites_on_node*sizeof(su3_vector));
  }

  /* Initiallize all the eigenvectors to a random vector */
  for(j=0;j<Nvecs;j++) {
    grad[j] = 1.0e+10 ;
    grsource_plain(F_OFFSET(g_rand), parity);  
    FORSOMEPARITY(i,s,parity){
      su3vec_copy(&(s->g_rand),&(eigVec[j][i]));
    }
    eigVal[j] = 1.0e+16 ;
  }

#ifdef EIGTIME
  dtimec = -dclock();
#endif

  while((max_error>Tolerance)&&(iter<Kiters)) {
    iter++ ;
    /* Run through all current eigenvectors, minimizing the
       Rayleigh-Ritz quotient in the space spanned by the current
       eigenvectors.  Replace the eigenvectors with each
       improvement. Stop improving the eigenvector when it has
       converged. */
    for(j=0;j<Nvecs;j++){
      if(grad[j]>(ToleranceG)) {
	converged[j] = 0 ;
	copy_Vector(eigVec[j],vec) ;
	total_iters += Rayleigh_min(vec, eigVec, ToleranceG, RelTol,
				    j, MaxIter , Restart, parity, fn[0]) ;
	/* Copy only wanted parity (UMH) */
	/** copy_Vector(vec,eigVec[j]) ; **/
	FORSOMEPARITY(i,s,parity){
	  su3vec_copy( &(vec[i]), &(eigVec[j][i]));
	}
      }else{ 
	converged[j] = 1;
      }
    }

    /* Diagonalize the operator on the subspace of improved
       eigenvectors and rotate the eigenvectors to the new basis. */

    /* if you didn't act on eigVec[i] last time, converged[i]=1,
       and  MeigVec hasn't changed, so don't compute it */
    constructArray(eigVec, MeigVec, &Array, grad, converged, parity, fn[0]) ;

#ifdef DEBUG
    node0_printf("Eigenvalues before diagonalization\n");
    for(i=0;i<Nvecs;i++){
      node0_printf("quot(%i) = %g |grad|=%g\n",i,Array.M[i][i].real,grad[i]);
    }
#endif

    Jacobi(&Array, &V, JACOBI_TOL) ;
    sort_eigenvectors(&Array,&V) ;
    RotateBasis(eigVec,&V,parity) ;
    RotateBasis(MeigVec,&V,parity) ;

    /* find the maximum error */
    max_error = 0.0 ;
    for(i=0;i<Nvecs;i++){
      /* First recompute gradient--vec=MeigVec - eigVal*eigVec */
      /* recall we are recycling MeigVec's */
      Array.M[i][i].imag=0.0;
      copy_Vector(MeigVec[i], vec) ;
      double_vec_mult_sub(&Array.M[i][i].real, eigVec[i], vec, parity) ;
      norm2(vec, &grad[i], parity) ;
      err[i] = eigVal[i] ;
      eigVal[i] = Array.M[i][i].real ;
      err[i] = fabs(err[i] - eigVal[i])/(1.0 - RelTol*RelTol) ;
#ifndef STRICT_CONVERGENCE
      if(err[i]>max_error){
	max_error = err[i] ;
      }
#else	
      if(grad[i]>max_error){
	max_error = grad[i] ;
      }
#endif
    }

    /* this makes a LOT of output! Im  commenting it out for a bit. -EBG */

    /*
    node0_printf("\nEigenvalues after diagonalization at iteration %i\n",
		 iter);
    for(i=0;i<Nvecs;i++)
      node0_printf("quot(%i) = %g +/- %8e |grad|=%g\n",
		   i,eigVal[i],err[i],grad[i]);

    */
  }
  
#ifdef EIGTIME
  dtimec += dclock();
  node0_printf("KAULKRITER: time = %e iters = %d iters/vec = %e\n",
	   dtimec,total_iters, (double)(total_iters)/Nvecs);
#endif

  node0_printf("Eigenvalues after diagonalization at iteration %i\n",iter);

  node0_printf("BEGIN RESULTS\n");
  for(i=0;i<Nvecs;i++){
    node0_printf("Eigenvalue(%i) = %g +/- %8e\t cvg? %d  \n",
		 i,eigVal[i],err[i],converged[i]);
  }

  /** Deallocate the arrays **/
  deAllocate(&V) ;
  deAllocate(&Array) ;
  free(err) ;
  free(grad) ;
  free(vec) ;
  free(converged) ;
  for(i=0;i<Nvecs;i++){
    free(MeigVec[i]);
  }
  free(MeigVec);
  cleanup_Matrix();

  return total_iters ;
}



