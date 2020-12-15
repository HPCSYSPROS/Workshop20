	subroutine noise (uny,unoy,m)
c this routine does not appear to be in use

c
c calculation of "horizontal spectrum" (in 2D kx-ky plane)
c
      use bouss_module
      use comsp
c
c     implicit none
	include 'intvars'
c
	complex(b8) :: uny(xisz,zjsz,ny,3)
	complex(b8) :: unoy(xisz,zjsz,ny,2)
	integer :: i,j,k,l,ik,m,lu
	real(b8) :: rk,rk2,term,term1,term2,enk,rksq,const
c
      common/rancom/iseed,rseed(2)
        integer iseed
        integer rseed

	integer :: xst
c	
	complex ctmp,twopii
	real(b8) ampsq, phase

	tfact(:)=2.
	tfact(1)=1.
	if (taskid.eq.0) write (6,*) 'enter noise'

	kx2(:)=b11(m)*kx(:)**2
	ky2(:)=b22(m)*ky(:)**2
	kz2(:)=b33(m)*kz(:)**2

	ky2(nyhp)=0.
	kz2(nzhp)=0.
c
	unoy(:,:,:,:)=0.
c
	iseed=1000000+taskid

	pi=atan(1.)*4.
	twopii=2.*pi*imagi
c
	sum=0.
	sum2=0.
	do 100 y=1,ny
	if (y.eq.nyhp) go to 100
	call masks (y)
	do 110 zp=1,zjsz
	z=zp+zjst-1
	if (z.eq.nzhp) go to 110
	do 120 xp=1,xisz
	x=xp+xist-1
 	rk2=kx2(x)+ky2(y)
	rk=sqrt(rk2+kz2(z))
	term1=uny(xp,zp,y,1)*conjg(uny(xp,zp,y,1))*b11(m)
	term2=uny(xp,zp,y,2)*conjg(uny(xp,zp,y,2))*b22(m)
	term3=uny(xp,zp,y,3)*conjg(uny(xp,zp,y,3))*b33(m)
	term=term1+term2+term3
	sum=sum+tfact(x)*(term1+term2+term3)
	sum2=sum2+tfact(x)*(term1+term2+term3)*2.*viscos*rk*rk
 120	continue
 110	continue
 100	continue
	call MPI_ALLREDUCE (sum,sum0,1,mpireal,MPI_SUM,
     1			MPI_COMM_WORLD,mpierr)
	if (taskid.eq.0) write (6,*) 'noise: sum0=',sum0
	tke0=sum0
	call MPI_ALLREDUCE (sum2,sum0,1,mpireal,MPI_SUM,
     1			MPI_COMM_WORLD,mpierr)
	if (taskid.eq.0) write (6,*) 'noise: sum0=',sum0
	diss0=sum0
	eta0=(viscos**3/diss0)**(.25)

	const=sqrt(beta1*beta2*beta3/2./pi)

      do 10 y=1,ny
	
	if (y.eq.nyhp) go to 10

	call masks(y)
	
      do 11 zp=1,zjsz
	z=zp+zjst-1
	if (z.eq.nzhp) go to 11
c
	xst=1
	if (y.eq.1.and.z.eq.1) xst=2
c
      do 12 xp=xst,xisz
	x=xp+xist-1

 	rk2=kx2(x)+ky2(y)
	rksq=rk2+kz2(z)
	rk=sqrt(rksq)
	
	if (x.eq.1.and.y.eq.1) then
	unoy(xp,zp,y,1)=1.
	unoy(xp,zp,y,2)=-1.
	else
	unoy(xp,zp,y,1)=beta1*ky(y)/sqrt(rk2)
	unoy(xp,zp,y,2)=-beta2*kx(x)/sqrt(rk2)
	end if
c
	phase=ranu(1)
	ampsq=enk(rk,kmax,eta0,noise_spconst,noise_spslope)
	
#ifdef OLD
c	unoy(xp,zp,y,1)=ctmp*unoy(xp,zp,y,1)/rk2/beta1
c	unoy(xp,zp,y,2)=ctmp*unoy(xp,zp,y,2)/rk2/beta2
 	unoy(xp,zp,y,1)=ctmp*unoy(xp,zp,y,1)/rksq/beta1
 	unoy(xp,zp,y,2)=ctmp*unoy(xp,zp,y,2)/rksq/beta2
#endif
 	ctmp=mask(xp,zp)*sqrt(ampsq)*cexp(twopii*phase)*const
  	unoy(xp,zp,y,1)=ctmp*unoy(xp,zp,y,1)/rk/beta1
  	unoy(xp,zp,y,2)=ctmp*unoy(xp,zp,y,2)/rk/beta2
 

 12	continue	
c

 11	continue
c
 10	continue
c
	sum=0.
	do 200 y=1,ny
	if (y.eq.nyhp) go to 200
	call masks (y)
	do 210 zp=1,zjsz
	z=zp+zjst-1
	if (z.eq.nzhp) go to 210
	do 220 xp=1,xisz
	x=xp+xist-1
	term1=unoy(xp,zp,y,1)*conjg(unoy(xp,zp,y,1))*b11(m)
	term2=unoy(xp,zp,y,2)*conjg(unoy(xp,zp,y,2))*b22(m)
	sum=sum+tfact(x)*(term1+term2)
 220	continue
 210	continue
 200	continue
	call MPI_ALLREDUCE (sum,sum0,1,mpireal,MPI_SUM,
     2			MPI_COMM_WORLD,mpierr)
c
	if (taskid.eq.0) write (6,*) 'noise: sum0=',sum0
	ratio=sqrt(noise_ratio*tke0/sum0)
c
	sum=0.
	do 300 y=1,ny
	if (y.eq.nyhp) go to 300
	call masks (y)
	do 310 zp=1,zjsz
	z=zp+zjst-1
	if (z.eq.nzhp) go to 310
	do 320 xp=1,xisz
	x=xp+xist-1
	unoy(xp,zp,y,1)=unoy(xp,zp,y,1)*ratio
	unoy(xp,zp,y,2)=unoy(xp,zp,y,2)*ratio
	term1=unoy(xp,zp,y,1)*conjg(unoy(xp,zp,y,1))*b11(m)
	term2=unoy(xp,zp,y,2)*conjg(unoy(xp,zp,y,2))*b22(m)
	sum=sum+tfact(x)*(term1+term2)
 320	continue
 310	continue
 300	continue
	call MPI_ALLREDUCE (sum,sum0,1,mpireal,MPI_SUM,
     3			MPI_COMM_WORLD,mpierr)
c
	if (taskid.eq.0) write (6,*) 'noise: sum0=',sum0
c
	do y=1,ny
	do zp=1,zjsz
	do xp=1,xisz
 	uny(xp,zp,y,1)=uny(xp,zp,y,1)+unoy(xp,zp,y,1)
 	uny(xp,zp,y,2)=uny(xp,zp,y,2)+unoy(xp,zp,y,2)
c	uny(xp,zp,y,1)=unoy(xp,zp,y,1)
c	uny(xp,zp,y,2)=unoy(xp,zp,y,2)
	end do
	end do
	end do

	if (taskid.eq.0) write (6,*) ' exit noise'

	return
	end
c
#ifdef TEST
	real function enk(rk,kmax,kscale,const,slope)
	real kc,kscale
      kc=kmax
	enk=0.
      if (rk.gt.kc) then
	enk=0.
	else
c	const=128.
c	const=512.
	term1=rk**4
	term2=(1.+const*(rk*kscale)**2)**slope
	enk=term1/term2
	end if
	return
	end
#endif
c
	real function enk(rk,kmax,kscale,const,slope)
	real kp,kscale
	theta=1.5
	eta=theta/kmax
	b=64.
c	kp=0.5
	kp=1.26
	if (rk.lt.kp) then
	enk=rk**4
	else
	enk=exp(-b*(rk*eta)**(4./3.))
	end if
	return
	end
#ifdef TEST
	real function enk(rk,kmax,kscale,const,slope)
	real kp,kscale
	theta=1.5
	eta=theta/kmax
	enk=0.
  	if (rk.ge.8.1.and.rk.le.12.1) enk=1.
	return
	end
#endif
