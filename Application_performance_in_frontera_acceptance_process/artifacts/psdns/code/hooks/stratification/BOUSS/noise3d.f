	subroutine noise3d (uny,unoy,m)
c
c
c 3D noise
c
	use stratification_module
      use comsp
c
      implicit none
	include 'intvars'
c
	complex(b8) :: uny(ny,zjsz*xisz,3)
	complex(b8) :: unoy(ny,zjsz*xisz,3)
	integer :: i,j,k,l,ik,m,lu
	real(b8) :: rk,rk2,term,term1,term2,enk,rksq,const
c
	real(b8) sum,sum2,term3,tke0,diss0,eta0,sum0,sum3,sum30
	integer a
	real(b8) ranu,ratio
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
	unoy(:,:,:)=0.
c
	iseed=1000000+taskid

	pi=atan(1.)*4.
	twopii=2.*pi*imagi
c
	sum=0.
	sum2=0.
	sum3=0.
c
      xp = mystart_x
      z = mystart_z
c
        do 100 a=1,num_al
           x=xp+xist-1

	if (z.eq.nzhp) go to 100
	do 120 y=1,ny
	if (.not.mask(y,a)) go to 120
 	rk2=kx2(x)+ky2(y)
	rk=sqrt(rk2+kz2(z))
	term1=uny(y,a,1)*conjg(uny(y,a,1))*b11(m)
	term2=uny(y,a,2)*conjg(uny(y,a,2))*b22(m)
	term3=uny(y,a,3)*conjg(uny(y,a,3))*b33(m)
	term=term1+term2+term3
	sum=sum+tfact(x)*(term1+term2+term3)
	sum2=sum2+tfact(x)*(term1+term2+term3)*2.*viscos*rk*rk
	sum3=sum3+tfact(x)*term3
 120	continue
	call next_xz (xp,z)
 100	continue
c
	call MPI_ALLREDUCE (sum,sum0,1,mpireal,MPI_SUM,
     1			MPI_COMM_WORLD,mpierr)
	if (taskid.eq.0) write (6,*) 'noise: sum0=',sum0
	tke0=sum0
	call MPI_ALLREDUCE (sum2,sum0,1,mpireal,MPI_SUM,
     1			MPI_COMM_WORLD,mpierr)
	if (taskid.eq.0) write (6,*) 'noise: sum0=',sum0
	diss0=sum0
	eta0=(viscos**3/diss0)**(.25)
	call MPI_ALLREDUCE (sum3,sum30,1,mpireal,MPI_SUM,
     1			MPI_COMM_WORLD,mpierr)
	if (taskid.eq.0) write (6,*) 'noise: sum30=',sum30
c
	if (taskid.eq.0) write (6,*) 'noise: before gauss3d'
	call gauss3d (unoy,m)
	if (taskid.eq.0) write (6,*) 'noise:  after gauss3d'

	sum=0.
      xp = mystart_x
      z = mystart_z
c
        do 200 a=1,num_al
           x=xp+xist-1
	if (z.eq.nzhp) go to 200
	
	do 220 y=1,ny
	if (.not.mask(y,a)) go to 220
	term1=unoy(y,a,1)*conjg(unoy(y,a,1))*b11(m)
	term2=unoy(y,a,2)*conjg(unoy(y,a,2))*b22(m)
	term3=unoy(y,a,3)*conjg(unoy(y,a,3))*b33(m)
	sum=sum+tfact(x)*(term1+term2+term3)
 220	continue
	call next_xz (xp,z)
c
 200	continue
	call MPI_ALLREDUCE (sum,sum0,1,mpireal,MPI_SUM,
     2			MPI_COMM_WORLD,mpierr)
c
	if (taskid.eq.1) write (6,*) 'noise: sum0=',sum0
	ratio=sqrt(noise_ratio*tke0/sum0)
	if (taskid.eq.1) write (6,*) 'noise:noise_ratio: tke0=',
     1         noise_ratio,tke0
c
ccccc
c	ratio=1.
ccccc
c
c
	sum=0.
      xp = mystart_x
      z = mystart_z
c
	sum3=0.
c
        do 300 a=1,num_al
           x=xp+xist-1
	if (z.eq.nzhp) go to 300
	
	do 320 y=1,ny
	if (.not.mask(y,a)) go to 320
	unoy(y,a,1)=unoy(y,a,1)*ratio
	unoy(y,a,2)=unoy(y,a,2)*ratio
	unoy(y,a,3)=unoy(y,a,3)*ratio
	term1=unoy(y,a,1)*conjg(unoy(y,a,1))*b11(m)
	term2=unoy(y,a,2)*conjg(unoy(y,a,2))*b22(m)
	term3=unoy(y,a,3)*conjg(unoy(y,a,3))*b33(m)
	sum=sum+tfact(x)*(term1+term2+term3)
	sum3=sum3+tfact(x)*term3
 320	continue
	call next_xz (xp,z)
c
 300	continue
	call MPI_ALLREDUCE (sum,sum0,1,mpireal,MPI_SUM,
     3			MPI_COMM_WORLD,mpierr)
	call MPI_ALLREDUCE (sum3,sum30,1,mpireal,MPI_SUM,
     3			MPI_COMM_WORLD,mpierr)
c
	if (taskid.eq.0) write (6,*) 'noise: sum0=',sum0
	if (taskid.eq.0) write (6,*) 'noise: sum30=',sum30
c
	do a=1,num_al
	do y=1,ny
 	uny(y,a,1)=uny(y,a,1)+unoy(y,a,1)
 	uny(y,a,2)=uny(y,a,2)+unoy(y,a,2)
 	uny(y,a,3)=uny(y,a,3)+unoy(y,a,3)
c	uny(y,a,1)=unoy(y,a,1)
c	uny(y,a,2)=unoy(y,a,2)
c	uny(y,a,3)=unoy(y,a,3)
	end do
	end do

	if (taskid.eq.0) write (6,*) ' exit noise3d'

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
