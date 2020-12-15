      subroutine sptvar (uny)
c
c hybrid version: note: slight round-off differences
c w.r.t. pure MPI version is possible when adding up
c contributions to spectrum from different threads

      use comsp
      implicit none
      include 'intvars'

      complex(b8) uny(ny,zjsz*xisz,3)
c
	real, allocatable :: eky(:),eks(:)
c
	real term,term1,term2,term3,rk2
	integer a,ik, icall
	data icall/0/
	save icall
c
	integer ithr,omp_get_thread_num,iy1,iy2,ii
c
	integer, allocatable :: ixp(:),iz(:)
        real, allocatable :: eky_thr(:,:)
c
	allocate (eky(nxh),eks(nxh))
c#ifdef OPENMP
	allocate (ixp(0:num_thr-1),iz(0:num_thr-1),eky_thr(nxh,0:num_thr-1))
c#endif
c
	icall=icall+1
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
!$OMP PARALLEL private (ithr,a,x,tfact_x,y,term1,term2,term3,rk2,ik,term,xp,z,iy1,iy2) shared (mystart_x,mystart_z)

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
	term1=real(uny(y,a,1)*conjg(uny(y,a,1)))
	term2=real(uny(y,a,2)*conjg(uny(y,a,2)))
	term3=real(uny(y,a,3)*conjg(uny(y,a,3)))
c
	rk2=kx2(x)+ky2(y)+kz2(z)
	ik=sqrt(rk2)+1.5

	term=tfact_x*(term1+term2+term3)*0.5
c
c	if (kinit(1).eq.0.and.istep.le.1) then
	if (kinit(1).eq.-2.and.term.gt.1.e-4) then
	write (6,"('sptvar,taskid,x,y,z,uny=',5i4,i6,1p,2e12.4)") taskid,x,y,z,ik,a,uny(y,a,1)

	end if
c	end if
c
	if (ik.gt.nxh.or.ik.lt.0) then 
	write (6,*) 'bad value of ik=',ik,x,y,z,taskid,ithr
	end if
c
	if (ik.eq.2) then
	write (6,"('sptvar: taskid, x,y,z,a,uny=',4i5,i7,1p,2e12.4)") taskid,x,y,z,a,uny(y,a,1)
	end if
c
	eky_thr(ik,ithr)=eky_thr(ik,ithr)+term
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

	call MPI_BARRIER (MPI_COMM_WORLD,mpierr)
	if (taskid.eq.0) then
	if (icall.eq.1) then
  	open (20,file='spectrum')
	else
 	open (20,file='spectrum',position='append')
	end if
	write (20,*) 'istep,kstep=',istep,kstep
	write (20,201) (eks(ik),ik=1,nxh)
 201	format ((1p,10e13.5))
	write (6,*) 'sptvar: sum of spectrum=',sum(eks(1:nxh))
 	close (20)
	end if
c
	deallocate (eky,eks)
	deallocate (ixp,iz,eky_thr)
c
	if (taskid.eq.0) write (6,*) ' exit sptvar,istep,kstep=',istep,kstep
	return
	end
