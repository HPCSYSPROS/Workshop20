! Input: X-pencils (real)
! Output: Z-pencils (complex)
!
! This routine performs real-to-complex FFT,
! truncates the X dimension from nxpad to nx, and
! transposes the data to arrange Z-pencils while 
! interchanging the order of X and Z indices (using loop blocking)
!
! Multivariable version
!

      subroutine xkcomm1_trans(source,dest,nv)

      use com
      use timers_comm
      use timers_comp
      use timers_tran
      implicit none
      include 'intvars'

#include "fft_stuff.f"

      real(b8) source(nxpad,zisz,yjsz,nv)
      complex(b8) dest(nzpad,xisz,yjsz,nv)
      complex(b8), allocatable :: buf(:,:,:) 

      integer position,pos0,pos1,pos2
      integer, allocatable :: pos(:)
      integer i,j,k,n,x2,ix,iz,ii,nv
      integer :: iarg
      real(b8) factor
      real xnorm,ynorm,znorm,tmpvmaxz,vel,vm
c
	real*8 rtime0,rtime1,rtime2
c
      complex(b8), allocatable :: recvbuf(:,:,:)
      complex(b8), allocatable :: sendbuf(:,:,:)

#ifndef USE_EVEN
      integer, allocatable :: sndcnts(:),sndstrt(:),rcvcnts(:),rcvstrt(:)
#endif
c
	integer ithr,omp_get_thread_num,iy1,yp1,yp2,i1,i2,iz1,iz2
	integer, allocatable :: iip1(:),iipi(:)
c
c next line moved to sub. transform, 12/24/2012
c	if (kstep.eq.1) t_xkcomm1(:,2)=0.

      rtime0=MPI_WTIME()

#ifdef USE_EVEN
      if(IfCntUneven) then
         allocate (pos(0:iproc-1))
      endif
      iarg = IfCntMax/(b8*2)
      allocate (sendbuf(iarg,nv,0:iproc-1))
      allocate (recvbuf(iarg,nv,0:iproc-1))
#endif

      allocate(buf(nxhppad,zisz,yjsz))
      factor=1./real(nxpad,b8)

      tp1_comm = tp1_comm - MPI_Wtime() 
c
	call divide_thr (yjsz,iyp1,iypi,'xkcomm1: yjsz')
	if (iproc.ge.num_thr) then
	allocate (iip1(0:num_thr-1),iipi(0:num_thr-1))
	call divide_thr (iproc,iip1,iipi,'xkcomm1: iproc')
	iip1(:)=iip1(:)-1
	end if
c
        t_xkcomm1(3,2)=t_xkcomm1(3,2)+(MPI_WTIME()-rtime0)

	ithr=0
!$OMP PARALLEL private (ithr,iy1,ii,pos0,z1,pos1,x2,pos2,position,n,i,j,k,iz,yp1,yp2,i1,i2,iz1,iz2,z2)
#ifdef OPENMP
      ithr = OMP_GET_THREAD_NUM()
#endif
	
#ifdef ESSL
          call rcft(1,source,nxpad,buf,nxhppad,nxpad,
     1      zisz*yjsz/num_thr,1,factor,
     &      raux1(1,ithr), rnaux1,raux2(1,ithr),rnaux2,raux3,1)           
#endif

      do 100 j=1,nv

!$OMP MASTER
	rtime1=MPI_WTIME()
!$OMP END MASTER

! Transform R2C

	iy1=iyp1(ithr)
c
#ifdef FFTW
         call fftw_execute_r2c(plan1_p2b(ithr),source(1,1,iy1,j),buf(1,1,iy1))
#elif defined ESSL
          call rcft(0,source(1,1,iy1,j),nxpad,buf(1,1,iy1),nxhppad,nxpad,
     1                zisz*yjsz/num_thr,1,factor,
     &      raux1(1,ithr),rnaux1,raux2(1,ithr),rnaux2,raux3,rnaux3)           
#endif
c
!$OMP MASTER
	rtime2=MPI_WTIME()
	tcpu_fft=tcpu_fft+(rtime2-rtime1)
        t_xkcomm1(2,2)=t_xkcomm1(2,2)+(rtime2-rtime1)
!$OMP END MASTER

c

c Pack the send buffer for exchanging z and x (within a given y plane ) into sendbuf
         
!$OMP MASTER
	rtime1=MPI_WTIME()
!$OMP END MASTER
c

#ifdef USE_EVEN
! If using MPI_Alltoall

!$OMP DO
      do yp=1,yjsz

          do i=0,iproc-1
             ii = mymap(i)
             position = (yp-1)*iisz(ii)*zisz+1
             do z=1,zisz
! Pack only those components with kx < nxh
                do x=iist(ii),iien(ii)
#ifdef ESSL

                   sendbuf(position,j,i) = buf(x,z,yp)
#else
                   sendbuf(position,j,i) = buf(x,z,yp) * factor
#endif
                position = position +1
                enddo
             enddo
             if(IfCntUneven) then
                pos(i) = position
             endif
          enddo

	end do
!$OMP END DO

#else
! If using MPI_Alltoallv
c
c Two possible ways of spreading workload among the threads
c
	if (yjsz.ge.num_thr.or.iproc.lt.num_thr) then
	yp1=iyp1(ithr)
	yp2=iyp1(ithr)+iypi(ithr)-1
 	i1=0
 	i2=iproc-1
	else
 	yp1=1
 	yp2=yjsz
 	i1=iip1(ithr)
 	i2=iip1(ithr)+iipi(ithr)-1
	end if
c
c      do yp=yp1,yp2
          do i=i1,i2
             ii = mymap(i)
      do yp=yp1,yp2
             pos0 = IfSndStrt(i)/(b8*2)*nv + (yp-1)*iisz(ii)*zisz +1
     &            + (j-1)*yjsz*zisz*iisz(ii)
	position=pos0
c
             do z=1,zisz
! Pack only those components with kx < nxh
                do x=iist(ii),iien(ii)
#ifdef ESSL
                   sendbuf(position,1,1) = buf(x,z,yp) 
#else
                   sendbuf(position,1,1) = buf(x,z,yp) * factor
#endif
                   position = position +1
                enddo
             enddo
c
          enddo
       enddo
#endif
c
!$OMP MASTER
	rtime2=MPI_WTIME()
	tcpu_other=tcpu_other+(rtime2-rtime1)
        t_xkcomm1(3,2)=t_xkcomm1(3,2)+(rtime2-rtime1)
!$OMP END MASTER
c
 100	continue
c
!$OMP BARRIER 
#ifdef USE_EVEN
! Not needed if sendbuf is preallocated
       if(IfCntUneven) then
!$OMP DO
          do i=0,iproc-1
             ii = mymap(i)
             position = pos(i)
             n = IfCntMax/(b8*2) - iisz(ii)*zisz*yjsz
             do j=1,nv
             do k=1,n
                sendbuf(position,j,i) = 0.0
                position = position +1
             enddo
          enddo
          enddo
!$OMP END DO
       endif
#endif       

!$OMP BARRIER

!$OMP MASTER


      tp1_comm = tp1_comm + MPI_Wtime() 
      ip_comm = ip_comm + 1

c Exchange the z-x buffers

#ifdef USE_EVEN

c
      t_alltoall = t_alltoall - MPI_Wtime()
      t1_comm = t1_comm - MPI_Wtime() 

	rtime1=MPI_WTIME()

      iarg = IfCntMax*nv
#ifdef CAF
      call compi_alltoall(sendbuf,recvbuf,iarg,mpi_comm_row)
#else
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

#endif

	rtime2=MPI_WTIME()

      t1_comm = t1_comm + MPI_Wtime() 
      t_alltoall = t_alltoall + MPI_Wtime()
      t_xkcomm1(1,2) = t_xkcomm1(1,2) + rtime2-rtime1

#else
c (not use_even)


      allocate(sndcnts(0:iproc-1))
      allocate(sndstrt(0:iproc-1))
      allocate(rcvcnts(0:iproc-1))
      allocate(rcvstrt(0:iproc-1))

      sndcnts = IfSndCnts * nv
      sndstrt = IfSndStrt*nv
      rcvcnts = IfRcvCnts*nv
      rcvstrt = IfRcvStrt*nv


! Exchange X-Z buffers

      t_alltoall = t_alltoall - MPI_Wtime()
      t1_comm = t1_comm - MPI_Wtime() 

        rtime1=MPI_WTIME()

#ifdef CAF
      call compi_alltoallv(sendbuf,SndCnts,SndStrt,IfCntMax*nv,
     &                     recvbuf,RcvCnts,RcvStrt,mpi_comm_row)
#else
      call mpi_alltoallv(sendbuf,SndCnts, SndStrt,mpi_byte, 
     &     recvbuf,RcvCnts, RcvStrt,mpi_byte,mpi_comm_row,ierr)
#endif

        rtime2=MPI_WTIME()
c
      t1_comm = MPI_Wtime() + t1_comm
      t_alltoall = t_alltoall + MPI_Wtime()

      t_xkcomm1(1,2) = t_xkcomm1(1,2) + rtime2-rtime1


      deallocate(sndcnts,sndstrt,rcvcnts,rcvstrt)


#endif

!$OMP END MASTER


!Next line (the barrier) appears necessary to enure correct results!
!$OMP BARRIER

!$OMP MASTER
	rtime1=MPI_WTIME()
!$OMP END MASTER
c

#ifdef USE_EVEN
! Unpack the receive buffer, interchanging indices, with loop blocking NB1

      do 200 j=1,nv
!$OMP DO
      do yp=1,yjsz

         do i =0,iproc-1
            ii = mymap(i)
            pos0 =(yp-1)*kisz(ii)*xisz

            do z=kist(ii),kien(ii),NB1_Z
               z2 = min(z+NB1_Z-1,kien(ii))
               pos1 = pos0 + (z-kist(ii))*xisz                

               do x=1,xisz,NB1_X
                  x2 = min(x+NB1_X-1,xisz)
                  
                  pos2 = pos1 +x
                  do iz=z,z2
                     position = pos2
                     do ix=x,x2
                        dest(iz,ix,yp,j) = recvbuf(position,j,i)
                        position = position +1
                     enddo
                     pos2 = pos2 + xisz
                  enddo
               enddo
            enddo
         enddo
      enddo

!$OMP END DO
 200	continue
c
!$OMP MASTER
	rtime2=MPI_WTIME()
	tcpu_other=tcpu_other+(rtime2-rtime1)
      t_xkcomm1(3,2) = t_xkcomm1(3,2) + rtime2-rtime1
!$OMP END MASTER
c
#else
! (not USE_EVEN)

! Unpack the receive buffer, interchanging indices, with loop blocking NB1

!$OMP MASTER
	rtime1=MPI_WTIME()
!$OMP END MASTER
c
 	if (iproc.ge.num_thr) then
 	yp1=1
 	yp2=yjsz
 	i1=iip1(ithr)
 	i2=iip1(ithr)+iipi(ithr)-1
 	else
 	yp1=iyp1(ithr)
 	yp2=iyp1(ithr)+iypi(ithr)-1
  	i1=0
  	i2=iproc-1
 	end if
c
      do 300 j=1,nv
c
      do i=i1,i2
         ii = mymap(i)
         

         do yp=yp1,yp2
c
          pos0 = IfRcvStrt(i)/(b8*2) *nv + (j-1)*yjsz*xisz*kisz(ii)
     1     + xisz*kisz(ii)*(yp-1)

            do z=kist(ii),kien(ii),NB1_Z
               z2 = min(z+NB1_Z-1,kien(ii))
            
               pos1 = pos0 + (z-kist(ii))*xisz 

               do x=1,xisz,NB1_X
                   x2 = min(x+NB1_X-1,xisz)
                  pos2 = pos1 + x 
                  
                  do iz=z,z2
                     position = pos2
                     do ix=x,x2
                        dest(iz,ix,yp,j) = recvbuf(position,1,1)
                        position = position +1
                     enddo
                     pos2 = pos2 + xisz
                  enddo
               enddo

            enddo
         enddo
      enddo
c
 300	continue

!$OMP MASTER
	rtime2=MPI_WTIME()
	tcpu_other=tcpu_other+(rtime2-rtime1)
      t_xkcomm1(3,2) = t_xkcomm1(3,2) + rtime2-rtime1
!$OMP END MASTER
c
#endif
!$OMP END PARALLEL

	rtime2=MPI_WTIME()

      deallocate(buf)
c
      deallocate(recvbuf)
      deallocate (sendbuf)

#ifdef USE_EVEN
      if(IfCntUneven) deallocate (pos)
#endif

	if (allocated(iipi)) deallocate (iip1,iipi)
c
      t_xkcomm1(3,2) = t_xkcomm1(3,2) + (MPI_WTIME()-rtime2)
c
      i1_comm = i1_comm + 1
c

      return
      end
