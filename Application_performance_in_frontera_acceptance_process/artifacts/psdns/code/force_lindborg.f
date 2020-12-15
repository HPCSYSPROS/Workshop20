!NPM converting to stride-1 (not complete)

	subroutine force_lindborg (uny,icall)
	
#ifdef LINDBORG
	use comsp
	
	implicit none
	include 'intvars'
	
	complex(b8) :: uny(ny,zjsz,xisz,3+nc)	  
	integer :: i,j,k,l,ik,m,icall
	real :: rk,rk2,term,term1,term2,enk_2d,rksq,const,rkh,kzabs
	real tke1,tke2,tke1_all,tke2_all,ranu,rate,kf_target
	real ratio1,ratio2,factor1,factor2
	real tke_before, tke_after,sum1,sum2,rate1
c
c	complex, allocatable :: ufh(:,:,:),ufv(:,:)
c	save ufh,ufv
c
      common/rancom/iseed,rseed(2)
        integer iseed
        integer rseed

	integer :: xst
	
	complex ctmp,twopii,cfactor
	real ampsq, phase

	go to (1,2), icall

 1	continue
c
c	allocate (ufh(xisz,ny,2),ufv(zisz,2))
c
c compute signal of given spectral characteristics and mean-square 1.0
c
	kf_target=ei_waveno

	tfact(:)=2.
	tfact(1)=1.

	m=2

	kx2(:)=b11(m)*kx(:)**2
	ky2(:)=b22(m)*ky(:)**2
c
	iseed=1000000+taskid

	twopii=2.*pi*(0,1)
c
	const=sqrt(beta1*beta2)/sqrt(pi)
c
	tke1=0.
	tke2=0.
	ufh(:,:,:)=0.
	ufv(:,:)=0.
c	
c horizontal modes
c
	if (zjst.eq.1) then
c
      do 10 y=1,ny
c
	if (y.eq.nyhp) go to 10
	call masks(y)
	xst=1
	if (y.eq.1) xst=2
c
      do 12 xp=xst,xisz
	x=xp+xist-1
 	rk2=kx2(x)+ky2(y)
	rkh=sqrt(rk2)
c
	ufh(xp,y,1)=beta1*ky(y)/rkh
	ufh(xp,y,2)=-beta2*kx(x)/rkh
c
	phase=ranu(1)
	ampsq=enk_2d(rkh,kf_target)
 	ctmp=sqrt(ampsq)*cexp(twopii*phase)*const
c
	cfactor=ctmp/sqrt(rkh)
   	ufh(xp,y,1)=cfactor*ufh(xp,y,1)/beta1
   	ufh(xp,y,2)=cfactor*ufh(xp,y,2)/beta2
c
	term1=ufh(xp,y,1)*conjg(ufh(xp,y,1))*b11(m)
        term2=ufh(xp,y,2)*conjg(ufh(xp,y,2))*b22(m)
        tke1=tke1+.5*tfact(x)*(term1+term2)
 12	continue	
c
 10     continue
c
	end if
c
c vertical modes
c
	if (xist.eq.1) then
	do zp=1,zjsz
	z=zp+zjst-1
c	kzabs=abs(kz(z))*beta3
	kzabs=abs(kz(z))
	if (xist.eq.1.and.(kzabs.gt.2.999.and.kzabs.lt.5.001)) then
	write (6,*) 'force_lindborg: z=',z
	ufv(zp,1)=1./beta1
	ufv(zp,2)=-1./beta2
	tke2=tke2+1.0
	end if
	end do
	end if
c
	call MPI_ALLREDUCE (tke1,tke1_all,1,MPI_REAL,MPI_SUM,
     1		            MPI_COMM_WORLD,mpierr)
	call MPI_ALLREDUCE (tke2,tke2_all,1,MPI_REAL,MPI_SUM,
     1		            MPI_COMM_WORLD,mpierr)
	if (taskid.eq.0) write (6,*) 'force_lindborg: tke1=',tke1_all
	if (taskid.eq.0) write (6,*) 'force_lindborg: tke2=',tke2_all
c
 	ratio1=sqrt(.99/tke1_all)
 	ratio2=sqrt(.01/tke2_all)
c	ratio1=sqrt(.5/tke1_all)
c	ratio2=sqrt(.5/tke2_all)
c
c scale the energies
c
	ufh(:,:,:)=ufh(:,:,:)*ratio1
	ufv(:,:)=ufv(:,:)*ratio2
c
c check the scaled energy
c
	tke1=0.
	tke2=0.
	if (zjst.eq.1) then
c
      do 110 y=1,ny
c
	if (y.eq.nyhp) go to 110
	call masks(y)
	xst=1
	if (y.eq.1) xst=2
c
      do 112 xp=xst,xisz
	x=xp+xist-1
	term1=ufh(xp,y,1)*conjg(ufh(xp,y,1))*b11(m)
        term2=ufh(xp,y,2)*conjg(ufh(xp,y,2))*b22(m)
        tke1=tke1+.5*tfact(x)*(term1+term2)
 112	continue	
c
 110     continue
c
	end if
c
	if (xist.eq.1) then
	do zp=1,zjsz
	z=zp+zjst-1
c	kzabs=abs(kz(z))*beta3
	kzabs=abs(kz(z))
	tke2=tke2+ufv(zp,1)*conjg(ufv(zp,1))/b11(m)
	end do
	end if
c
	call MPI_ALLREDUCE (tke1,tke1_all,1,MPI_REAL,MPI_SUM,
     1		            MPI_COMM_WORLD,mpierr)
	call MPI_ALLREDUCE (tke2,tke2_all,1,MPI_REAL,MPI_SUM,
     1		            MPI_COMM_WORLD,mpierr)
	if (taskid.eq.0) write (6,*) 'force_lindborg: tke1=',tke1_all
	if (taskid.eq.0) write (6,*) 'force_lindborg: tke2=',tke2_all
c
c
	return
c
 2	continue
c
c copy into 3D array
c recognize that Rogallo's algorithm stores u/beta1, v.beta2
c
c	if (taskid.eq.0) write (6,*) 'enter force_lindborg: icall=2'
c
	
	call force_eirate (uny,m,rate1)

 	factor1=dt
 	factor2=dt

	sum1=0.
	sum2=0.
	
	if (zjst.eq.1) then
	   zp=1
	   do 30 y=1,ny
	      if (y.eq.nyhp) go to 30
	      xst=1
	      if (y.eq.1) xst=2
	      do 32 xp=xst,xisz
		 x=xp+xist-1
		 term1=uny(y,zp,xp,1)*conjg(uny(y,zp,xp,1))*b11(m)
		 term2=uny(y,zp,xp,2)*conjg(uny(y,zp,xp,2))*b22(m)
		 sum1=sum1+.5*tfact(x)*(term1+term2)
#ifdef TEST
		 write (910+taskid,901) xp,y,ufh(xp,y,1),uny(y,zp,xp,1),factor1
		 write (910+taskid,901) xp,y,ufh(xp,y,2),uny(y,zp,xp,2)
#endif
		 uny(y,zp,xp,1)=uny(y,zp,xp,1)+ufh(xp,y,1)*factor1
		 uny(y,zp,xp,2)=uny(y,zp,xp,2)+ufh(xp,y,2)*factor1
		 term1=uny(y,zp,xp,1)*conjg(uny(y,zp,xp,1))*b11(m)
		 term2=uny(y,zp,xp,2)*conjg(uny(y,zp,xp,2))*b22(m)
		 sum2=sum2+.5*tfact(x)*(term1+term2)
#ifdef TEST
		 if (istep.eq.1) then
		    write (900+taskid,901) xp,y,ufh(xp,y,1),uny(y,zp,xp,1)
		    write (900+taskid,901) xp,y,ufh(xp,y,2),uny(y,zp,xp,2)
 901		    format (2i4,1p,5e12.4)
		 end if
#endif

 32	      continue	
 30	   continue
	end if
c       
	if (xist.eq.1) then
	   do zp=1,zjsz
	      z=zp+zjst-1
c	kzabs=abs(kz(z))*beta3
	      kzabs=abs(kz(z))
	      if ((kzabs.gt.2.999.and.kzabs.lt.5.001)) then
		 term1=uny(1,zp,1,1)*conjg(uny(1,zp,1,1))*b11(m)
		 term2=uny(1,zp,1,2)*conjg(uny(1,zp,1,2))*b22(m)
		 sum1=sum1+.5*tfact(x)*(term1+term2)
		 uny(1,zp,1,1)=uny(1,zp,1,1)+ufv(zp,1)*factor2
		 uny(1,zp,1,2)=uny(1,zp,1,2)+ufv(zp,2)*factor2
		 term1=uny(1,zp,1,1)*conjg(uny(1,zp,1,1))*b11(m)
		 term2=uny(1,zp,1,2)*conjg(uny(1,zp,1,2))*b22(m)
		 sum2=sum2+.5*tfact(x)*(term1+term2)
c	write (6,*) 'vertical mode: istep,z,uny=',istep,z,uny(1,zp,1,1)
	      end if
	   end do
	end if
c       
	call MPI_REDUCE (sum1,tke_before,1,mpireal,MPI_SUM,0,
	1  MPI_COMM_WORLD,mpierr)
	call MPI_REDUCE (sum2,tke_after, 1,mpireal,MPI_SUM,0,
	1  MPI_COMM_WORLD,mpierr)
	if (taskid.eq.0) then
	   write (670,701) istep,tke_after-tke_before
 701	   format ('istep,energy added in force_lindborg=',i6,1p,e13.5)
	end if
c       

c	if (taskid.eq.0) write (6,*) ' exit force_lindborg: icall=2'
#endif
c       
	return
	end

c       
#ifdef TEST
	real function enk_2d(rk,kmax)
	real kp,kscale
	enk_2d=0.
c 	if (rk.ge.12.1.and.rk.le.20.1) enk_2d=1.
  	if (rk.ge.22.1.and.rk.le.30.1) enk_2d=1.
c 	if (rk.le.1.1) enk_2d=1.
	return
	end
#endif
	real function enk_2d(rkh,kscale)
	real kscale
	enk_2d=rkh*exp(-(rkh-kscale)**2)
	return
	end
