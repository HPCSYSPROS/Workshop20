      subroutine epfor_param (tl_hat,kmaxeta)
c

	use com

	implicit none
	real tl_hat,n_f,eta_t,kmaxeta,eps_t,efdot,c3,tl,c2,sigma2
	real re_hat,rlam,c1,c5,c4,beta_EP
c
#ifdef SBP_2010
#ifndef FEK_FORC
c
	if (tforce.lt.0.) then
	tforce=abs(tforce)
	return
	end if
c
	beta_EP=0.8
c
	pi=atan(1.)*4.
c
	c1=(3.*kforce/16./pi)**(1./3.)
	c2=(3./16./pi)**(1./3.)*kforce**(-5./3.)
	c3=1/beta_EP * (.25)**(1./3.)*kforce**(-2./3.)/c2
	c4=8.5*(4.*pi/3.)**(-2./9.)*kforce**(1./6.)
	c5=c4*c1**(-5./8.)

c	c1=0.55
c	c2=0.07
c	c3=5.7
c	c5=10.7
c
c
	if (taskid.eq.0) write (6,*) 'kforce,epsfor=',kforce,epsfor
c
c 6/10/11: use float(nx) instead of nx below, to avoid possible
c integer overflow for large nx,ny,nz
c
c	kmax=sqrt(2.)*(nx*ny*nz*beta1*beta2*beta3)**(1./3.)/3
	kmax=sqrt(2.)*(float(nx)*ny*nz*beta1*beta2*beta3)**(1./3.)/3
c
	n_f=4.*pi/3.*kforce**3/beta1/beta2/beta3
c
	eta_t=(kmaxeta)/kmax
c
	eps_t=viscos**3/eta_t**4
c 
	efdot=eps_t*(1+c3*tl_hat)

	tl=tl_hat/(c2*efdot**(1./3.)*kforce**(2./3.))
c
	tforce=tl
c
	sigma2=efdot/4./n_f/tl
c
 	epsfor=sigma2*tl
c
	re_hat=c1*(efdot)**(1./3.)/kforce**(4./3.)/viscos
c
	rlam=c5*(eps_t/efdot)**(5./24.)*re_hat**(5./8.)
c
	if (taskid.eq.0) then
        write (6,*) 'epfor_param: c1,c2,c3=',c1,c2,c3
        write (6,*) 'epfor_param: c4,c5=',c4,c5
	write (6,*) 'epfor_param: eta_t,eps_t=',eta_t,eps_t
	write (6,*) 'epfor_param: efdot,n_f=',efdot,n_f
	write (6,*) 'epfor_param: kmaxeta,kmax=',kmaxeta,kmax
	write (6,*) 'epfor_param: tforce,epsfor=',tforce,epsfor
	write (6,*) 'epfor_param: re_hat,rlam=',re_hat,rlam
	end if
#endif
#endif

 
      return
      end
