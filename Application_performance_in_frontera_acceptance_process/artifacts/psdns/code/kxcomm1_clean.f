
! Version with Cylindrical representation 
!
! Updates by PKY, 1/16/2012: assume jproc.gt.num_thr
!
! Starting with Y-cylinders, wavenumber-space in all three dimensions
! - Inverse-Transform in Y, 
! - transpose to get data in z-pencils with stride-1 in Z
!
! Multivariable version
!
!     subroutine kxcomm1_trans_cyl (source,dest,nv)
      subroutine kxcomm1_clean (source,dest,nv)

      use com
      use timers_comm
      use timers_comp
      use timers_tran
	use timers_rkstep

      implicit none
      include 'intvars'
#ifdef FFTW
#include "fft_stuff.f"
#endif

      complex(b8) source(ny,zjsz*xisz,nv)
      complex(b8) dest(nz,xisz,yjsz,nv)

      complex(b8), allocatable :: recvbuf(:,:,:)
      complex(b8), allocatable :: sendbuf(:,:,:)
      integer position,pos0,pos1,n,a
      integer, allocatable :: pos(:,:)
      integer i,j,iz,dnz,iy,ia,a1,a2,c,nv
      real vmax
	integer :: iarg
c
	real*8 rtime0,rtime1,rtime2
c
        integer ithr,omp_get_thread_num,iy1,iy2,num_fft,ia1,jthr,NB2Y
        integer, allocatable :: ja_st(:)

       integer, allocatable :: jpp1(:),jppi(:)

        if (kstep.eq.1) t_kxcomm1(:,2)=0.
c

	rtime1=MPI_WTIME()

      if(num_al_i(0) .gt. 0) then
c
        allocate (sendbuf(jjsz(0)*num_al_i(0),nv,0:jproc-1))
        allocate (recvbuf(jjsz(0)*num_al_i(0),nv,0:jproc-1))
c
c ASSUME jproc.ge.num_thr
c
        allocate (jpp1(0:num_thr-1),jppi(0:num_thr-1))
        call divide_thr (jproc,jpp1,jppi,'xkcomm2: jproc')
        jpp1(:)=jpp1(:)-1

! This loop does the following:
! - for each x, it executes zjsz transforms in Y, to reduce memory references
! - it packs the send buffer while rearranging the data 
!   so Z has stride 1, using loop blocking for cache optimization

        allocate(pos(0:jproc-1,0:num_thr-1))
	
        allocate (ja_st(0:num_thr-1))
        ia=1
        do jthr=0,num_thr-1
        ja_st(jthr)=ia
        ia=ia+num_fft0(jthr)
        end do
c
      end if

      rtime2=MPI_WTIME()
      tcpu_other=tcpu_other+(rtime2-rtime1)
      t_kxcomm1(3,2)=t_kxcomm1(3,2)+(rtime2-rtime1)

      call divide_thr (yjsz,iyp1,iypi,'kxcomm1: yjsz')

      ithr=0
c
!$OMP PARALLEL private (ithr,y,a2,pos1,ia,ia1,position,jthr,i,n,c,a1,iy1,iy2,x,z,a)

#ifdef OPENMP
      ithr = OMP_GET_THREAD_NUM()
#endif

      if(num_al_i(0) .gt. 0) then

!$OMP MASTER
        rtime1=MPI_WTIME()
!$OMP END MASTER
        do j=jpp1(ithr),jpp1(ithr)+jppi(ithr)-1
        sendbuf (:,:,j)=0.
        end do
!$OMP MASTER
        rtime2=MPI_WTIME()
        t_kxcomm1(3,2)=t_kxcomm1(3,2)+(rtime2-rtime1)
!$OMP END MASTER

        if(num_al .gt. 0) then
c
#ifdef ESSL
            call escpft (source,1,ny,ny,num_fft0(ithr),1,ithr)
#endif
c
         do 100 j=1,nv

!$OMP MASTER
          rtime1=MPI_WTIME()
!$OMP END MASTER

          ia1=ja_st(ithr)
#ifdef FFTW
            call fftw_execute(plan4_p1(ithr),source(1,ia1,j),source(1,ia1,j))
#elif defined ESSL
            call escfft (source(1,ia1,j),1,ny,ny,num_fft0(ithr),1,ithr)
#endif

!$OMP MASTER
        rtime2=MPI_WTIME()
        tcpu_fft=tcpu_fft+(rtime2-rtime1)
        t_kxcomm1(2,2)=t_kxcomm1(2,2)+(rtime2-rtime1)
!$OMP END MASTER


!$OMP BARRIER

!$OMP MASTER
        rtime1=MPI_WTIME()
!$OMP END MASTER

            do i=jpp1(ithr),jpp1(ithr)+jppi(ithr)-1

                do y=jjst(i),jjen(i),NB2_Y


                  y2 = min(y+NB2_Y-1,jjen(i))
               
                  do a=1,num_al,NB2_Z
                     a2 = min(a+NB2_Z-1,num_al)
                  
                     pos1 = a + (y-jjst(i))*num_al
               
                     do iy=y,y2
                        position = pos1
                        do ia=a,a2
                           sendbuf(position,j,i) = source(iy,ia,j)
                           position = position +1
                        enddo
                        pos1 = pos1 + num_al
                     enddo
                  enddo
               enddo
            
            enddo

!$OMP MASTER
        rtime2=MPI_WTIME()
        tcpu_other=tcpu_other+(rtime2-rtime1)
        t_kxcomm1(3,2)=t_kxcomm1(3,2)+(rtime2-rtime1)
!$OMP END MASTER

 100	continue
      
c
      endif

!$OMP BARRIER
!$OMP MASTER

      t_alltoall = t_alltoall - MPI_Wtime()
      t3_comm = t3_comm - MPI_Wtime() 

	iarg=jjsz(0)*num_al_i(0)*b8*nv*2
c
	rtime1=MPI_WTIME()

!      call mpi_alltoall(sendbuf,iarg,mpi_byte,
!     &      		recvbuf,iarg,mpi_byte,mpi_comm_col,ierr)

#ifdef CAF
      call compi_alltoall(sendbuf,recvbuf,iarg,mpi_comm_col)
#else
      call mpi_alltoall(sendbuf,iarg,mpi_byte,
     &      		recvbuf,iarg,mpi_byte,mpi_comm_col,ierr)
#endif

	rtime2=MPI_WTIME()

      t3_comm = t3_comm + MPI_Wtime() 
      t_alltoall = t_alltoall + MPI_Wtime()
        t_kxcomm1(1,2)=t_kxcomm1(1,2)+(rtime2-rtime1)
!$OMP END MASTER
!$OMP BARRIER


!$OMP MASTER
        rtime1=MPI_WTIME()
!$OMP END MASTER
c
	do j=1,nv
	do y=iyp1(ithr),iyp1(ithr)+iypi(ithr)-1
	do x=max_al_x+1,xisz
	do z=1,nz
	dest(z,x,y,j)=0.
	end do
	end do
	end do
	end do
c

	do y=iyp1(ithr),iyp1(ithr)+iypi(ithr)-1
	
	do i=0,jproc-1
	pos(i,ithr)=(y-1)*num_al_i(i)+1
	end do

            i = 0
            n=0
            do x=1,max_al_x
               c = cut_z(x)
               a1 = nzhp - c -1
               a2 = nzhp + c +1
               do z=1,a1
                  dest(z,x,y,1:nv) = recvbuf(pos(i,ithr),1:nv,i) 
                  pos(i,ithr) = pos(i,ithr)+1
                  n = n+1
                  if(n .eq. num_al_i(i)) then
                     n = 0
                     i = i+1
                  endif
               enddo
               
	do z=a1+1,a2-1
	dest(z,x,y,1:nv)=0.
	end do
	
               
               do z=a2,nz
                  dest(z,x,y,1:nv) = recvbuf(pos(i,ithr),1:nv,i) 

                  pos(i,ithr) = pos(i,ithr)+1
                  n = n+1
                  if(n .eq. num_al_i(i)) then
                     n = 0
                     i = i+1
                  endif
               enddo
               
            enddo
         enddo
c
         
!$OMP MASTER
      i3_comm = i3_comm + 1
        rtime2=MPI_WTIME()
        tcpu_other=tcpu_other+(rtime2-rtime1)
        t_kxcomm1(3,2)=t_kxcomm1(3,2)+(rtime2-rtime1)
!$OMP END MASTER

      else
! if num_al_i(0) = 0
c
!$OMP MASTER
        rtime1=MPI_WTIME()
!$OMP END MASTER
c
         do j=1,nv
	do y=iyp1(ithr),iyp1(ithr)+iypi(ithr)-1
               do x=1,xisz
                  do z=1,nz
                     dest(z,x,y,j) = 0.0
                  enddo
               enddo
            enddo
         enddo

!$OMP MASTER
        rtime2=MPI_WTIME()
        tcpu_other=tcpu_other+(rtime2-rtime1)
        t_kxcomm1(3,2)=t_kxcomm1(3,2)+(rtime2-rtime1)
!$OMP END MASTER
c
      endif


!$OMP END PARALLEL

c

      if(num_al_i(0) .gt. 0) then
c
      deallocate (pos)
c
	deallocate (ja_st)
      deallocate (recvbuf)
      deallocate (jpp1,jppi)
c
	end if

!RAF Must allocate coarrays on all images, so all ranks must deallocate them
      if (allocated(sendbuf)) deallocate(sendbuf)
c
      return
      end
