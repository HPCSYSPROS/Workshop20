      subroutine press_stat

! Last deallocate statement in this routine added on 10/3/2012
c

      use comsp
      implicit none
      include 'intvars'

	real, allocatable :: eky(:),eks(:)
c
	real term,rk2,sum
	integer a,ik
c
	integer ithr,omp_get_thread_num,iy1,iy2,ii
c
	integer, allocatable :: ixp(:),iz(:)
        real, allocatable :: eky_thr(:,:)
	integer icall
	save icall
	data icall/0/
c
        real(b8) beta_min

#ifdef SHELL_DK
        beta_min=min(beta1,beta2,beta3)
#endif

#ifdef CVC_PRESS
	if (.not.flagc) return
c
	allocate (eky(nxh),eks(nxh))
	allocate (ixp(0:num_thr-1),iz(0:num_thr-1),eky_thr(nxh,0:num_thr-1))
c

	kx2(:)=kx(:)**2
	ky2(:)=ky(:)**2
	kz2(:)=kz(:)**2
c
	ky2(nyhp)=0.
	kz2(nzhp)=0.
c
	eky_thr=0.
	eky=0.
c
	ithr=0
c
!$OMP PARALLEL private (ithr,a,x,tfact_x,y,rk2,ik,term,xp,z,iy1,iy2) shared (mystart_x,mystart_z)

#ifdef OPENMP
      ithr = OMP_GET_THREAD_NUM()
#endif

      xp = mystart_x
      z = mystart_z
c
      do a=1,num_al
                                                                                
         x = xp + xist-1
                                                                                
         if(x .eq. 1) then
            tfact_x=1
         else
            tfact_x=2
         endif
                                                                                
	iy1=ithr*ny/num_thr+1
	iy2=(ithr+1)*ny/num_thr

         do 10 y=iy1,iy2
c
	if (.not.mask(y,a)) go to 10
c
	term=real(upy(y,a)*conjg(upy(y,a)))
c
	rk2=kx2(x)+ky2(y)+kz2(z)
#ifdef SHELL_DK
        ik=sqrt(rk2)/beta_min+1.5
#else
        ik=sqrt(rk2)+1.5
#endif

	eky_thr(ik,ithr)=eky_thr(ik,ithr)+term*tfact_x
c
 10	continue
c
#ifdef OPENMP
	ixp(ithr)=xp
	iz(ithr)=z
  	call next_xz(ixp(ithr),iz(ithr))
	xp=ixp(ithr)
	z=iz(ithr)
#else
  	call next_xz(xp,z)
#endif

	end do
c
!$OMP BARRIER
!$OMP MASTER
	do ii=0,num_thr-1
	do ik=1,nxh
	eky(ik)=eky(ik)+eky_thr(ik,ii)
	end do
	end do
	
	call MPI_REDUCE (eky,eks,nxh,mpireal,
     1            MPI_SUM,0,MPI_COMM_WORLD,mpierr)
!$OMP END MASTER
!$OMP BARRIER
c
!$OMP END PARALLEL

	if (taskid.eq.0) then
	sum=0.
	do ik=1,nxh
	sum=sum+eks(ik)
	end do
c
	icall=icall+1
  	if (icall.eq.1) open (200,file='press_sptr')
	write (200,201) istep,time
	write (200,202) sum,(2.*tke/3.)**2
	write (200,"('pressure spectrum')")
	write (200,203) (eks(ik),ik=1,nxh)
 201 	format ('istep=',i6,' time=',1p,e12.4) 
 202	format ('pressure variance and urms**4=',1p,2e12.4)
 203	format ((1p,10e13.5))
	end if
c
	deallocate (eky,eks)
	deallocate (ixp,iz,eky_thr)
c
#endif
	return
	end
