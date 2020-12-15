! reduced version of kxcomm2, used for the "ut" array
! (which does not have the y-dimension) in proc2a and proc2b 
!
! This routine reorders X,Z indices while packing the send buffer
! (using loop blocking),
! exchanges data to arrange X pencils, expandsX dimension to nxhppad
! and performs inverse complex-to-real FFT
! Input: Z pencils, complex
! Output: X-pencils, real 
!
! Multivariable version
!
! sy is parameter indicating which y-plane this is
!    

      subroutine kxcomm2t_trans(source,dest,npy,yplane,nv)

      use com
      use timers_comm
      use timers_comp
      use timers_tran
      implicit none
      include 'intvars'

#include "fft_stuff.f"

	integer npy
      real(b8) dest(nxpad,zisz,npy,nv)
      complex(b8) source(nzpad,xisz,npy,nv)

      complex(b8), allocatable :: recvbuf(:,:,:),buf12(:,:,:)
      complex(b8), allocatable :: sendbuf(:,:,:)

      integer, allocatable :: sndcnts(:),sndstrt(:),rcvcnts(:),rcvstrt(:)
      integer :: iarg
      integer position,pos0,pos1,pos2
      integer i,j,k,n,ix,iz,x2,ii,nv,yplane
        integer ithr,omp_get_thread_num,iz1,iz2

	integer ostep
	data ostep/0/
	real*8 rtime0,rtime1,rtime2
c
      real(b8) norm
c
        integer, allocatable :: jzp1(:,:),jzpi(:,:)
	integer mzdim,NB1Z

        allocate (jzp1(0:num_thr-1,0:iproc-1))
        allocate (jzpi(0:num_thr-1,0:iproc-1))
c
        call MPI_BARRIER (MPI_COMM_WORLD,mpierr)

        if (jstep.gt.ostep.and.kstep.eq.1) then
	ostep=jstep
c       t_kxcomm2t(:,2)=0.
        end if

      norm = 1.

	rtime0=MPI_WTIME()
	rtime1=MPI_WTIME()

	call divide_thr (zisz,izp1,izpi,'kxcomm2t: zisz')
	mzdim=0
	do ithr=0,num_thr-1
	mzdim=max0(mzdim,izpi(ithr))
	end do
	
c

#ifdef USE_EVEN
      iarg = KrCntMax/(b8*2)
      allocate (sendbuf(iarg,nv,0:iproc-1))
      allocate (recvbuf(iarg,nv,0:iproc-1))
#else
      allocate (recvbuf(zisz*nxh*nv,1,1))
      allocate (sendbuf(nzpad*xisz*nv,1,1))
#endif


c     allocate(buf12(nxhppad,zisz/num_thr,0:num_thr-1))
      allocate(buf12(nxhppad,mzdim,0:num_thr-1))
c
	rtime2=MPI_WTIME()
	t_kxcomm2t(3,2)=t_kxcomm2t(3,2)+(rtime2-rtime1)

	tcpu_other=tcpu_other+MPI_WTIME()-rtime1
c
! Pack send buffer
! Interchange indices, using loop blocking NB1

c
	ithr=0
c
	NB1Z=NB1_Z
	if (kisz(0)/NB1_Z.lt.num_thr) NB1Z=1
c	if (mod(kisz(0),NB1_Z).ne.0) NB1Z=kisz(0)/2
	if (mod(kisz(0),NB1_Z).ne.0) NB1Z=4
c
         do i=0,iproc-1
	ii=mymap(i)
        call divide_thr (kisz(ii)/NB1Z,jzp1(0,i),jzpi(0,i),'kxcomm2t: kisz/NB1Z')
        end do
c
!$OMP PARALLEL private (ithr,ii,z2,pos0,x2,pos1,position,z,x,iz,n,pos2,iz1,j,ix,iz2,i)

#ifdef OPENMP
      ithr = OMP_GET_THREAD_NUM()
#endif


#ifdef USE_EVEN
! Use MPI_Alltoall

!$OMP MASTER
        rtime1=MPI_WTIME()
!$OMP END MASTER

c
!$OMP DO
      do i=0,iproc-1
         ii = mymap(i)
         
         do z=kist(ii),kien(ii),NB1Z
            z2 = min(z+NB1Z-1,kien(ii))
            pos0 = (z-kist(ii))*xisz
            
            do x=1,xisz,NB1_X
               x2 = min(x+NB1_X-1,xisz)
               pos1 = pos0 + x
               
               do iz=z,z2
                  position = pos1
                  do ix=x,x2
                     sendbuf(position,1:nv,i) = source(iz,ix,yplane,1:nv)
                     position = position +1
                  enddo
                  pos1 = pos1 + xisz
               enddo
            enddo
         enddo
         n = KrCntMax/(b8*2) - kisz(i)*xisz
         do k=1,n
            sendbuf(position,1:nv,i) = 0.0
            position = position +1
         enddo
      enddo
!$OMP END DO

!$OMP MASTER
        rtime2=MPI_WTIME()
        tcpu_other=tcpu_other+(rtime2-rtime1)
	t_kxcomm2t(3,2)=t_kxcomm2t(3,2)+(rtime2-rtime1)
!$OMP END MASTER


! Exchange x-z buffers

!$OMP BARRIER
!$OMP MASTER
      t_alltoall = t_alltoall - MPI_Wtime()
      t4t_comm = t4t_comm - MPI_Wtime() 

	rtime1=MPI_WTIME()

      iarg = KrCntMax*nv
#ifdef CAF
      call compi_alltoall(sendbuf,recvbuf,iarg,mpi_comm_row)
#else
      call mpi_alltoall(sendbuf,iarg,mpi_byte,
     &      		recvbuf,iarg,mpi_byte,mpi_comm_row,ierr)
#endif

	rtime2=MPI_WTIME()

      t4t_comm = t4t_comm + MPI_Wtime() 
      t_alltoall = t_alltoall + MPI_Wtime()
	t_kxcomm2t(1,2)=t_kxcomm2t(1,2)+(rtime2-rtime1)
!$OMP END MASTER
!$OMP BARRIER

#else
! Not USE EVEN
! Use MPI_Alltoallv

! Pack send buffer
! Interchange indices, using loop blocking NB1

!$OMP MASTER
        rtime1=MPI_WTIME()
!$OMP END MASTER
c
	do j=1,nv
c
      do i=0,iproc-1
         ii = mymap(i)
c
         pos0 = KrSndStrt(i)/(b8*2) * nv + (j-1)*xisz*kisz(ii)

        do z1=jzp1(ithr,i),jzp1(ithr,i)+jzpi(ithr,i)-1
        z=kist(ii)+(z1-1)*NB1Z

            z2 = min(z+NB1Z-1,kien(ii))
            pos1 = pos0 + (z-kist(ii))*xisz 
            
            do x=1,xisz,NB1_X
               x2 = min(x+NB1_X-1,xisz)
               pos2 = pos1 + x
               
               do iz=z,z2
                  position = pos2
                  do ix=x,x2
                     sendbuf(position,1,1) = source(iz,ix,yplane,j)
                     position = position +1
                  enddo
                  pos2 = pos2 + xisz
               enddo
            enddo
         enddo

	end do


      enddo

!$OMP BARRIER
!$OMP MASTER

! Exchange x-z buffers
      allocate(sndcnts(0:iproc-1))
      allocate(sndstrt(0:iproc-1))
      allocate(rcvcnts(0:iproc-1))
      allocate(rcvstrt(0:iproc-1))

      sndcnts = KrSndCnts * nv
      sndstrt = KrSndStrt*nv
      rcvcnts = KrRcvCnts*nv
      rcvstrt = KrRcvStrt*nv

        rtime2=MPI_WTIME()
        tcpu_other=tcpu_other+(rtime2-rtime1)
	t_kxcomm2t(3,2)=t_kxcomm2t(3,2)+(rtime2-rtime1)
	t_alltoall = t_alltoall - MPI_Wtime()
      t4t_comm = t4t_comm - MPI_Wtime() 
        
        rtime1=MPI_WTIME()

#ifdef CAF
      call compi_alltoallv(sendbuf,SndCnts,SndStrt,KrCntMax*nv,
     &                     recvbuf,RcvCnts,RcvStrt,mpi_comm_row)
#else
      call mpi_alltoallv (sendbuf,SndCnts, SndStrt,mpi_byte, 
     &     recvbuf,RcvCnts,RcvStrt,mpi_byte,mpi_comm_row,ierr)
#endif

        rtime2=MPI_WTIME()

      t4t_comm = t4t_comm + MPI_Wtime() 
	t_alltoall = t_alltoall + MPI_Wtime()
	t_kxcomm2t(1,2)=t_kxcomm2t(1,2)+(rtime2-rtime1)
      deallocate(sndcnts,sndstrt,rcvcnts,rcvstrt)


!$OMP END MASTER
!$OMP BARRIER

#endif
c

! Unpack receive buffers 


#ifdef ESSL
         call crft(1,buf12,nxhppad,dest,nxpad,nxpad,zisz/num_thr,-1,norm,
     &      raux1(1,ithr), rnaux1,raux2(1,ithr),rnaux2,raux3,rnaux3)
#endif

	iz1=izp1(ithr)
	iz2=izp1(ithr)+izpi(ithr)-1
c
	do 300 j=1,nv
c
!$OMP MASTER
        rtime1=MPI_WTIME()
!$OMP END MASTER
c
         do 70 i=0,iproc-1
            ii = mymap(i)
#ifdef USE_EVEN
            position = 1 + (iz1-1)*iisz(ii)
#else
            position = KrRcvStrt(i)/(b8*2) * nv + 1 + (j-1)*zisz*iisz(ii)
     1               + (iz1-1)*iisz(ii)
#endif
            do z=1,iz2-iz1+1
               do x=iist(ii),iien(ii)
#ifdef USE_EVEN
                  buf12(x,z,ithr) = recvbuf(position,j,i)
#else
                  buf12(x,z,ithr) = recvbuf(position,1,1)
#endif
                  position = position+1
               enddo
            enddo
 70	continue
         
c!$OMP BARRIER
         
! Add and zero extra elements in X

         do z=1,iz2-iz1+1
            do x=nxhp,nxhppad
               buf12(x,z,ithr) = 0.0
            enddo
         enddo
c

c!$OMP BARRIER

!$OMP MASTER
        rtime2=MPI_WTIME()
        tcpu_other=tcpu_other+(rtime2-rtime1)
        t_kxcomm2t(3,2)=t_kxcomm2t(3,2)+(rtime2-rtime1)
        rtime1=MPI_WTIME()
!$OMP END MASTER
c

! C2R Transform    
#ifdef FFTW
         call fftw_execute_c2r(plan3_p2a(ithr),buf12(1,1,ithr),dest(1,iz1,yplane,j))
#elif defined ESSL
         call crft(0,buf12(1,1,ithr),nxhppad,dest(1,iz1,yplane,j),nxpad,nxpad,
     1           zisz/num_thr,-1,norm,
     &      raux1(1,ithr), rnaux1,raux2(1,ithr),rnaux2,raux3,rnaux3)
#endif
c
!$OMP MASTER
        rtime2=MPI_WTIME()
        tcpu_fft=tcpu_fft+(rtime2-rtime1)
        t_kxcomm2t(2,2)=t_kxcomm2t(2,2)+(rtime2-rtime1)
!$OMP END MASTER

 300	continue

c!$OMP BARRIER


!$OMP END PARALLEL
         
        rtime1=MPI_WTIME()
      deallocate(sendbuf)
      deallocate (recvbuf)
      deallocate (buf12)

!RAF one more to plug memory leak
      deallocate (jzp1,jzpi)

        rtime2=MPI_WTIME()
        tcpu_fft=tcpu_fft+(rtime2-rtime1)
        t_kxcomm2t(3,2)=t_kxcomm2t(3,2)+(rtime2-rtime1)


      i4t_comm = i4t_comm + 1


        t_kxcomm2t(4,2)=t_kxcomm2t(4,2)+(MPI_WTIME()-rtime0)

      return
      end


