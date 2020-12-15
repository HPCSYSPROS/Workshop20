! Version with Cylindrical representation 
!
! Updates by PKY, 1/16/2012: assume jproc.gt.num_thr
!
! Starting with data in Z-pencils, X and Z in wavenumber space, 
! Y in real space:
! - transpose data into Y-cylinders, changing indices ordering as well
!   to arrange the data in stride-1 in Y
! - Then do a forward transform in Y
!
c This version has improvements made by Dmitry, Dec 22, 2011
c
! Multivariable version
!
!     subroutine xkcomm2_trans_cyl(source,dest,nv)
      subroutine xkcomm2_clean (source,dest,nv)

      use com
      use timers_comm
      use timers_comp
      use timers_tran
      implicit none
      include 'intvars'
#ifdef FFTW
#include "fft_stuff.f"
#endif

      complex(b8) source(nzpad,xisz,yjsz,nv)
      complex(b8) dest(nypad,xisz*zjsz,nv)

      integer position,pos0,pos1,pos2,n,a
      integer i,j,iz,dnz,iy,ia,a1,a2,c,nv
      real xnorm,ynorm,znorm,tmpvmaxz,vel,vm

      integer, allocatable :: pos(:,:)
      complex(b8), allocatable :: recvbuf(:,:,:)
      complex(b8), allocatable :: sendbuf(:,:,:)

	integer :: iarg
	integer :: ii,jj
c
	real*8 rtime0,rtime1,rtime2,rtime3,rtime4
c
	integer lu
c
        integer ithr,omp_get_thread_num,iy1,iy2,ia1
 	integer, allocatable :: ja_st(:)

	integer ip,index,jthr,jzthr,NB2Y
c
	integer, allocatable :: jpp1(:),jppi(:)
c
	if (kstep.eq.1) t_xkcomm2(:,2)=0.

	rtime1=MPI_WTIME()
c
      if(num_al_i(0) .gt. 0) then
        allocate (sendbuf(jjsz(0)*num_al_i(0),nv,0:jproc-1))
        allocate (recvbuf(jjsz(0)*num_al_i(0),nv,0:jproc-1))
      endif
        sendbuf = 0.
c

	rtime1=MPI_WTIME()-rtime1
	tcpu_other=tcpu_other+rtime1
        t_xkcomm2(3,2)=t_xkcomm2(3,2)+rtime1

        rtime1=MPI_WTIME()

        allocate(pos(0:jproc-1,0:num_thr-1))
c
c ASSUME jproc.ge.num_thr
c
	allocate (jpp1(0:num_thr-1),jppi(0:num_thr-1))
	call divide_thr (jproc,jpp1,jppi,'xkcomm2: jproc')
	jpp1(:)=jpp1(:)-1

	ithr=0

!$OMP PARALLEL private(ithr,i,j,x,y,z,n,c,a1,a2,iy1,iy2,ip,index,jthr,pos1,position,ia,ia1,y1,y2)
c
#ifdef OPENMP
	ithr = OMP_GET_THREAD_NUM()
#endif


      do i=0,jproc-1
         pos(i,ithr) = 1+num_al_i(i)*(iyp1(ithr)-1)
      enddo
	
	iy1=iyp1(ithr)
	iy2=iyp1(ithr)+iypi(ithr)-1

      do y=iy1,iy2
         i = 0
         n=0
         do x=1,max_al_x
            c = cut_z(x) 
	    a1 = nzhp - c -1
            a2 = nzhp + c +1
            do z=1,a1

	sendbuf(pos(i,ithr),1:nv,i) = source(z,x,y,1:nv)
               pos(i,ithr) = pos(i,ithr)+1
               n = n+1
               if(n .eq. num_al_i(i)) then
                  n = 0
                  i = i+1
               endif
            enddo
            
            do z=a2,nz
	sendbuf(pos(i,ithr),1:nv,i) = source(z,x,y,1:nv)
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
!$OMP BARRIER
!$OMP MASTER
c

        rtime1=MPI_WTIME()-rtime1
        tcpu_other=tcpu_other+rtime1
        t_xkcomm2(3,2)=t_xkcomm2(3,2)+rtime1

      t_alltoall = t_alltoall - MPI_Wtime()
      t2_comm = t2_comm - MPI_Wtime() 


	iarg=jjsz(0)*num_al_i(0)*b8*nv*2
c
	rtime1=MPI_WTIME()


#ifdef CAF
      call compi_alltoall(sendbuf,recvbuf,iarg,mpi_comm_col)
#else
#ifdef ALLTOALLV
        atac_count(:)=iarg
        atac_start(0)=0
        do j=1,jproc-1
        atac_start(j)=atac_start(j-1)+atac_count(j-1)
        end do
      call mpi_alltoallv (sendbuf,atac_count,atac_start,mpi_byte,
     &      recvbuf,atac_count,atac_start,mpi_byte,mpi_comm_col,ierr)
#else
       call mpi_alltoall(sendbuf,iarg,mpi_byte,
     1                 recvbuf,iarg,mpi_byte,mpi_comm_col,ierr)
#endif
#endif


        rtime1=MPI_WTIME()-rtime1
        t_xkcomm2(1,2)=t_xkcomm2(1,2)+rtime1
c
      t2_comm = t2_comm + MPI_Wtime() 
      t_alltoall = t_alltoall + MPI_Wtime()

      deallocate (pos)
c
!$OMP END MASTER
!$OMP BARRIER

c
      if(num_al .gt. 0) then

#ifdef ESSL		
            call escpft (dest,1,nypad,nypad,num_fft0(ithr),-1,ithr)
#endif

!$OMP MASTER
	allocate (ja_st(0:num_thr-1))
	ia=1
	do jthr=0,num_thr-1
	ja_st(jthr)=ia
	ia=ia+num_fft0(jthr)
	end do

!$OMP END MASTER
!$OMP BARRIER

c
           do j=1,nv

!$OMP MASTER
        rtime2=MPI_WTIME()
!$OMP END MASTER

         do i=jpp1(ithr),jpp1(ithr)+jppi(ithr)-1
c
            do y=jjst(i),jjen(i),NB2_Y
               y2 = min(y+NB2_Y-1,jjen(i))
               do a=1,num_al,NB2_Z
                  a2 = min(a+NB2_Z-1,num_al)
                  pos1 = a + (y-jjst(i))*num_al
                  do iy=y,y2
                     position = pos1
                     do ia=a,a2
                         dest(iy,ia,j) = recvbuf(position,j,i)
                        position = position +1
                     enddo
                     pos1 = pos1 + num_al
                  enddo
               enddo
            enddo
            
	end do
c	end if

!$OMP BARRIER

!$OMP MASTER
        rtime2=MPI_WTIME()-rtime2
        tcpu_other=tcpu_other+rtime2
        t_xkcomm2(3,2)=t_xkcomm2(3,2)+rtime2

        rtime3=MPI_WTIME()
!$OMP END MASTER

	ia1=ja_st(ithr)
#ifdef FFTW
            call fftw_execute(plan4_p3(ithr),dest(1,ia1,j),dest(1,ia1,j))
#elif defined ESSL
        call escfft (dest(1,ia1,j),1,nypad,nypad,num_fft0(ithr),-1,ithr)
#endif
        

!$OMP BARRIER
!$OMP MASTER
        rtime3=MPI_WTIME()-rtime3
        tcpu_fft=tcpu_fft+rtime3
        t_xkcomm2(2,2)=t_xkcomm2(2,2)+rtime3
!$OMP END MASTER

      enddo
c
c Bug fix on 4/8/2012
c This deallocate can give an error if placed outside this endif,
c since num_al could be zero for a minority of MPI tasks

!$OMP MASTER
 	deallocate (ja_st)
!$OMP END MASTER

      endif
         
!$OMP END PARALLEL

c
      deallocate (recvbuf)

      i2_comm = i2_comm + 1

	if (allocated(jppi)) deallocate (jpp1,jppi)
c
! num_al_i(0) > 0

!RAF Must allocate coarrays on all images, so all ranks must deallocate them
      if (allocated(sendbuf)) deallocate (sendbuf)

      return
      end
