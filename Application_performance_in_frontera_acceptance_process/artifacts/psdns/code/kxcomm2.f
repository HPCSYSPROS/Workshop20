! Update of kxcomm2, by PKY Jan 14, 2012, DP Jan 13, 2012
! all y-planes, but 1 alltoallv call per variable  being transposed

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

      subroutine kxcomm2_trans(source,dest,nv)

      use com
      use timers_comm
      use timers_comp
      use timers_tran
      implicit none
      include 'intvars'

#include "fft_stuff.f"

c        complex(b8) dest(nxpad/2,zisz,yjsz,nv)
        real(b8) dest(nxpad,zisz,yjsz,nv)
c
      complex(b8) source(nzpad,xisz,yjsz,nv)

      complex(b8), allocatable :: recvbuf(:,:,:),buf12(:,:,:)
      complex(b8), allocatable :: sendbuf(:,:,:)

#ifndef USE_EVEN
      integer, allocatable :: sndcnts(:),sndstrt(:),rcvcnts(:),rcvstrt(:)
#endif

      integer :: iarg
      integer position,pos0,pos1,pos2
      integer i,j,k,n,ix,iz,x2,ii,nv

        integer ostep
        data ostep/0/
	real*8 rtime1,rtime2,rtime3,rtime0
c
      real(b8) norm
c
        integer ithr,omp_get_thread_num,iz1,iz2
!        integer, allocatable :: jzp1(:,:),jzpi(:,:)
c
	rtime0=MPI_WTIME()
        rtime1=MPI_WTIME()
c
c already allocated in main program
c       allocate (iyp1(0:num_thr-1))
c        allocate (iypi(0:num_thr-1))

	call divide_thr (yjsz,iyp1,iypi,'kxcomm2: yjsz')
c
      norm = 1.

        if (jstep.gt.ostep.and.kstep.eq.1) then
	ostep=jstep
c        t_kxcomm2(:,2)=0.
        end if
c

#ifdef USE_EVEN

      iarg = KrCntMax*yjsz/(b8*2)
      allocate (sendbuf(iarg,nv,0:iproc-1))
      allocate (recvbuf(iarg,nv,0:iproc-1))

#else
      allocate (sendbuf(nzpad*xisz*yjsz,nv,1))
      allocate (recvbuf(zisz*nxh*yjsz,nv,1))

      allocate(sndcnts(0:iproc-1))
      allocate(sndstrt(0:iproc-1))
      allocate(rcvcnts(0:iproc-1))
      allocate(rcvstrt(0:iproc-1))
      sndcnts = KrSndCnts*yjsz
      sndstrt = KrSndStrt*yjsz
      rcvcnts = KrRcvCnts*yjsz
      rcvstrt = KrRcvStrt*yjsz
#endif

      allocate(buf12(nxhppad,zisz,0:num_thr-1))

        rtime2=MPI_WTIME()
        t_kxcomm2(3,2)=t_kxcomm2(3,2)+(rtime2-rtime1)

        ithr=0

!$OMP PARALLEL private (ithr,ii,z2,pos0,x2,pos1,position,z,x,iz,n,pos2,iz1,j,ix,iz2,i)

#ifdef OPENMP
      ithr = OMP_GET_THREAD_NUM()
#endif


! Pack send buffer
! Interchange indices, using loop blocking NB1

#ifdef USE_EVEN
! Use MPI_Alltoall
c
!$OMP MASTER
        rtime1=MPI_WTIME()
!$OMP END MASTER

      do 100 j=1,nv
c
        do yp=iyp1(ithr),iyp1(ithr)+iypi(ithr)-1

         do i=0,iproc-1
            ii = mymap(i)

                pos0 = xisz * kisz(ii) * (yp-1)

               do z=kist(ii),kien(ii),NB1_Z
                  z2 = min(z+NB1_Z-1,kien(ii))
                  pos1 = pos0 + (z-kist(ii))*xisz

                  do x=1,xisz,NB1_X
                     x2 = min(x+NB1_X-1,xisz)
                     pos2 = pos1 + x

                     do iz=z,z2
                        position = pos2
                        do ix=x,x2
                           sendbuf(position,j,i) = source(iz,ix,yp,j)
                           position = position +1
                        enddo
                        pos2 = pos2 + xisz
                     enddo
                  enddo
               enddo
            enddo
         enddo
 100  continue

!$OMP BARRIER
!$OMP MASTER


#ifdef TIMERS
        rtime2=MPI_WTIME()
        tcpu_other=tcpu_other+(rtime2-rtime1)
        t_kxcomm2(3,2)=t_kxcomm2(3,2)+(rtime2-rtime1)
#endif

! Exchange x-z buffers

#ifdef TIMERS	
      t_alltoall = t_alltoall - MPI_Wtime()
      t4_comm = t4_comm - MPI_Wtime() 
#endif	
	rtime3=MPI_WTIME()


      iarg = KrCntMax*nv*yjsz
#ifdef ALLTOALLV
        atar_count(:)=iarg
        atar_start(0)=0
        do i=1,iproc-1
        atar_start(i)=atar_start(i-1)+atar_count(i-1)
        end do
      call mpi_alltoallv (sendbuf,atar_count,atar_start,mpi_byte,
     &      recvbuf,atar_count,atar_start,mpi_byte,mpi_comm_row,ierr)
#else
      call mpi_alltoall(sendbuf,iarg,mpi_byte,
     &                  recvbuf,iarg,mpi_byte,mpi_comm_row,ierr)
#endif
#ifdef CAF
      call compi_alltoall(sendbuf,recvbuf,iarg,mpi_comm_row)
#endif

	rtime3=MPI_WTIME()-rtime3
	t_kxcomm2(1,2)=t_kxcomm2(1,2)+rtime3
#ifdef TIMERS	
      t4_comm = t4_comm + MPI_Wtime() 
      t_alltoall = t_alltoall + MPI_Wtime()
#endif	

!$OMP END MASTER
!$OMP BARRIER
#else
! Not USE EVEN
! Use MPI_Alltoallv

! Pack send buffer
! Interchange indices, using loop blocking NB1

c
      do 200 j=1,nv

!$OMP MASTER
        rtime1=MPI_WTIME()
!$OMP END MASTER

      pos0 = 0
      do i=0,iproc-1
         ii = mymap(i)

c           do yp=1,yjsz
            do yp=iyp1(ithr),iyp1(ithr)+iypi(ithr)-1
               pos0 = xisz * kisz(ii) * (yp-1) + KrSndStrt(i)/(b8*2) *yjsz
               do z=kist(ii),kien(ii),NB1_Z
                  z2 = min(z+NB1_Z-1,kien(ii))
                  pos1 = pos0 + (z-kist(ii))*xisz 
            
                  do x=1,xisz,NB1_X
                     x2 = min(x+NB1_X-1,xisz)
                     pos2 = pos1 + x
               
                     do iz=z,z2
                        position = pos2
                        do ix=x,x2
                           sendbuf(position,j,1) = source(iz,ix,yp,j)
                           position = position +1
                        enddo
                        pos2 = pos2 + xisz
                     enddo
                  enddo
               enddo
            enddo
         enddo
c
! Exchange x-z buffers

!$OMP BARRIER 
!$OMP MASTER
        rtime2=MPI_WTIME()
        tcpu_other=tcpu_other+(rtime2-rtime1)
        t_kxcomm2(3,2)=t_kxcomm2(3,2)+(rtime2-rtime1)
c
	t_alltoall = t_alltoall - MPI_Wtime()
      t4_comm = t4_comm - MPI_Wtime() 

	rtime3=MPI_WTIME()
c
#ifdef CAF
        call compi_alltoallv(sendbuf(1,j,1),SndCnts,SndStrt,KrCntMax*yjsz,
     &                       recvbuf(1,j,1),RcvCnts,RcvStrt,mpi_comm_row)
#else
        call mpi_alltoallv (sendbuf(1,j,1),SndCnts, SndStrt,mpi_byte, 
     &       recvbuf(1,j,1),RcvCnts,RcvStrt,mpi_byte,mpi_comm_row,ierr)
#endif
c
	rtime3=MPI_WTIME()-rtime3
	t_kxcomm2(1,2)=t_kxcomm2(1,2)+rtime3
      t4_comm = t4_comm + MPI_Wtime() 
	t_alltoall = t_alltoall + MPI_Wtime()
c
!$OMP END MASTER
c
 200	continue



! This OMP barrier is important: otherwise some threads go ahead
! to unpack elements of the recvbuf array before everything is ready
!$OMP BARRIER 


#endif


! Unpack receive buffers 


#ifdef ESSL
         call crft(1,buf12,nxhppad,dest,nxpad,nxpad,zisz,-1,norm,
     &      raux1(1,ithr), rnaux1,raux2(1,ithr),rnaux2,raux3,rnaux3)
#endif

      do 300 j=1,nv
            do yp=iyp1(ithr),iyp1(ithr)+iypi(ithr)-1

!$OMP MASTER
            rtime1=MPI_WTIME()
!$OMP END MASTER

            do i=0,iproc-1
               ii = mymap(i)
#ifdef USE_EVEN
               position = 1
#else
               position = KrRcvStrt(i)/(b8*2) *yjsz + 1
#endif
               position = position + (yp-1)*iisz(ii)*zisz 
               do z=1,zisz
                  do x=iist(ii),iien(ii)
#ifdef USE_EVEN
                     buf12(x,z,ithr) = recvbuf(position,j,i)
#else
                     buf12(x,z,ithr) = recvbuf(position,j,1)
#endif
                     position = position+1
                  enddo
               enddo
            enddo
         
         
! Add and zero extra elements in X
!
            do z=1,zisz
               do x=nxhp,nxhppad
                  buf12(x,z,ithr) = 0.0
               enddo
            enddo
!$OMP MASTER
            rtime2=MPI_WTIME()
            tcpu_other=tcpu_other+(rtime2-rtime1)
            t_kxcomm2(3,2)=t_kxcomm2(3,2)+(rtime2-rtime1)
! C2R Transform    
            rtime1=MPI_WTIME()
!$OMP END MASTER
#ifdef FFTW
             call fftw_execute_c2r(plan3_p2c(ithr),buf12(1,1,ithr),dest(1,1,yp,j))
#else
        call crft(0,buf12(1,1,ithr),nxhppad,dest(1,1,yp,j),nxpad,nxpad,zisz,-1,norm,
     &           raux1(1,ithr),rnaux1,raux2(1,ithr),rnaux2,raux3,rnaux3)

#endif
!$OMP MASTER
            rtime2=MPI_WTIME()
            tcpu_fft=tcpu_fft+(rtime2-rtime1)
            t_kxcomm2(2,2)=t_kxcomm2(2,2)+(rtime2-rtime1)
!$OMP END MASTER
         enddo
 300  continue
         
!$OMP END PARALLEL

        rtime1=MPI_WTIME()

#ifndef USE_EVEN
      deallocate(sndcnts,sndstrt,rcvcnts,rcvstrt)
#endif
      deallocate(sendbuf)
      deallocate (recvbuf)
      deallocate (buf12)

        rtime2=MPI_WTIME()
        t_kxcomm2(3,2)=t_kxcomm2(3,2)+(rtime2-rtime1)

      i4_comm = i4_comm + 1

	t_kxcomm2(4,2)=t_kxcomm2(4,2)+(MPI_WTIME()-rtime0)
c
      return
      end
