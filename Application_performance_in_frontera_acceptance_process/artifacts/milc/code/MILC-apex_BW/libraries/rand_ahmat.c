/******************  rand_ahmat.c  (in su3.a) ***************************
*									*
* void random_anti_hermitian( anti_hermitmat *mat_antihermit, passthru *prn_pt)*
* Creates gaussian random anti-hermitian matrices			*
* Normalization is < |m01|^2 > = 1, or < m01.real*m01.real > = 1/2	*
* The argument "prn_pt" is a pointer to be passed to gaussian_rand_no() *
* RS6000 may choke on void *						*
*/
#include "../include/config.h"
#include <math.h>
#include "../include/complex.h"
#include "../include/su3.h"
#include <stdio.h>

void random_anti_hermitian( anti_hermitmat *mat_antihermit, double_prn *prn_pt) {
complex r38;
Real sqrt_third;

	sqrt_third = sqrt( (double)(1.0/3.0) );
	r38 = complex_gaussian_rand_no(prn_pt);
	mat_antihermit->m00im=r38.real+sqrt_third*r38.imag;
	mat_antihermit->m11im= -r38.real+sqrt_third*r38.imag;
	mat_antihermit->m22im= -2.0*sqrt_third*r38.imag;
	mat_antihermit->m01=complex_gaussian_rand_no(prn_pt);
	mat_antihermit->m02=complex_gaussian_rand_no(prn_pt);
	mat_antihermit->m12=complex_gaussian_rand_no(prn_pt);

}/*random_anti_hermitian_*/
