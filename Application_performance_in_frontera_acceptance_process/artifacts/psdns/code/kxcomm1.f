
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
      subroutine kxcomm1_trans_cyl (source,dest,nv)

      use com
      use timers_comm
      use timers_comp
      use timers_tran

      implicit none
      include 'intvars'
#ifdef FFTW
#include "fft_stuff.f"
#endif

      integer i,j,iz,dnz,iy,ia,a1,a2,c,nv
      complex(b8) source(ny,zjsz*xisz,nv)
      complex(b8) dest(nz,xisz,yjsz,nv)

      complex(b8), allocatable :: recvbuf(:,:,:)
      complex(b8), allocatable :: sendbuf(:,:,:)

      integer position,pos0,pos1,n,a
      integer, allocatable :: pos(:,:)
      real vmax
	integer :: iarg
        integer icall
        data icall/0/
        save icall
c
        integer ostep
        data ostep/0/
	real*8 rtime0,rtime1,rtime2
c
        integer ithr,omp_get_thread_num,iy1,iy2,num_fft,ia1,jthr,NB2Y
        integer, allocatable :: ja_st(:)

       integer, allocatable :: jpp1(:),jppi(:)

        icall=icall+1
         if (jstep.gt.ostep.and.kstep.eq.1) then
         ostep=jstep
c        t_kxcomm1(:,2)=0.
         end if


c
	rtime1=MPI_WTIME()

      if(num_al_i(0) .gt. 0) then

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
!$OMP PARALLEL private (ithr,y,a2,pos1,ia,ia1,position,jthr,i,n,c,a1,iy1,iy2,x,z,a,y2)

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
        t_kxcomm1(2,2)=t_kxcomm1(2,2)+(rtime2-rtime1)
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
     1      recvbuf,atac_count,atac_start,mpi_byte,mpi_comm_col,ierr)
#else
       call mpi_alltoall(sendbuf,iarg,mpi_byte,
     1                 recvbuf,iarg,mpi_byte,mpi_comm_col,ierr)
#endif

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


      if(num_al_i(0) .gt. 0) then
c
      deallocate (pos)
c
	deallocate (ja_st)
      deallocate (sendbuf)
      deallocate (recvbuf)
      deallocate (jpp1,jppi)
c
	end if
c
      return
      end

!-----------------------------------------------------------------
! Square pencils vesion (used in inpen_read)
!
! Starting with Y-pencils, wavenumber-space in all three dimensions
! - transpose to get data in z-pencils with stride-1 in Z
! (only transpose, no transform)
!
! change made to ensure correct execution on Cray compiler, 6/1/2012
!-----------------------------------------------------------------

      subroutine kxcomm1_pen (source,dest,nv)

      use com
      use timers_comm
      implicit none
      include 'intvars'
      integer position,pos0,pos1,n,nv,cnt
      complex(b8) source(ny,zjsz,xisz,nv)
      complex(b8) dest(nz,xisz,yjsz,nv)

      complex(b8), allocatable :: sendbuf(:,:,:),recvbuf(:,:)
#ifdef CAF_TEST
      complex(b8), allocatable :: cafrbuf(:,:)
#endif
      integer, allocatable :: pos(:)
! SndCnts(:),RcvCnts(:),SndStrt(:),RcvStrt(:)
      integer i,j,iz,dnz,iy,a,ia,a1,a2,c
      real vmax
c
	integer position1

      allocate (sendbuf(xisz*zjsz*ny,nv,1))
c      allocate (recvbuf(xisz*yjsz*nz,nv,1))
      allocate (recvbuf(xisz*yjsz*nz,nv))
#ifdef CAF_TEST
      allocate (cafrbuf(xisz*yjsz*nz,nv))
#endif

      t3_comm = t3_comm - MPI_Wtime() 
         
! This loop does the following:
! - for each x, it executes zjsz transforms in Y, to reduce memory references
! - it packs the send buffer while rearranging the data 
!   so Z has stride 1, using loop blocking for cache optimization

      do j=1,nv

         do x=1,xisz

            pos0 = (x-1)*zjsz 

            do y=1,ny,NB2_Y
               y2 = min(y+NB2_Y-1,ny)
               
               do z=1,zjsz,NB2_Z
                  z2 = min(z+NB2_Z-1,zjsz)
                  pos1 = pos0 +z
                  
                  do iy=y,y2
                     position = pos1 
                     do iz=z,z2
                        sendbuf(position,j,1) = source(iy,iz,x,j)
                        position = position +1
                     enddo
                     pos1 = pos1 + zjsz * xisz
                  enddo
               enddo
               pos0 = pos0 + xisz*zjsz*NB2_Y
            enddo
         enddo
         
      enddo
      
      cnt = xisz * yjsz * zjsz * b8 * 2

#ifdef TIMERS	
      t_alltoall = t_alltoall - MPI_Wtime()
#endif	

#if defined (CAF_TEST) || !defined (CAF)
      call mpi_alltoall(sendbuf,cnt,mpi_byte,
     &     recvbuf,cnt, mpi_byte,mpi_comm_col,ierr)
#endif
#ifdef CAF
      call compi_alltoall(sendbuf,recvbuf,cnt,mpi_comm_col)
#endif
#ifdef CAF_TEST
      call compi_alltoall(sendbuf,cafrbuf,cnt,mpi_comm_col)
      if (any(cafrbuf.ne.recvbuf)) then
        print *,'kxcomm1_pen_taskid correctness failed',taskid
      else
        print *,'kxcomm1_pen_taskid correctness passed',taskid
      endif
#endif

#ifdef TIMERS	
	t_alltoall = t_alltoall + MPI_Wtime()
#endif	

! Unpack the receive buffer 

      do j=1,nv

! PKY, 6/1/2012: the original version of this loop
! (where 'position' was incremented by 1 each time)
! does not execute correctly when using the Cray compiler,
! presumably because of some automatic loop re-structuring issues
! (results were correct if an IF or I/O statement was placed
! inside the z-loop to break up the flow).

         position = 1
	position1=1
         do i=0,jproc-1

            do y=1,yjsz
               do x=1,xisz
                  do z=kjst(i),kjen(i)
		position=position1+z-kjst(i)
                     dest(z,x,y,j) = recvbuf(position,j)
                  enddo
	position1=position1+kjsz(i)
               enddo
            enddo
         enddo

      enddo
      
      t3_comm = t3_comm + MPI_Wtime() 
      i3_comm = i3_comm + 1

      deallocate(sendbuf)
      deallocate (recvbuf)
#ifdef CAF_TEST
      deallocate (cafrbuf)
#endif

      return
      end
!-----------------------------------------------------------------
! Version with Cylindrical representation (used in dwrtu2)
!
! Starting with Y-cylinders, wavenumber-space in all three dimensions,
! - transpose to get data in square z-pencils with X dimension leading
! (only transpose, no transform)
!-----------------------------------------------------------------
 
      subroutine kxcomm1_cyl2sq (source,dest,nv)

      use com
      use timers_comm
      implicit none
      include 'intvars'
      integer position,pos0,pos1,n,nv
      complex(b8) source(ny,zjsz*xisz,nv)
      complex(b8) dest(xisz,nz,yjsz,nv)

      complex(b8), allocatable :: sendbuf(:,:,:),recvbuf(:,:,:)
#ifdef CAF_TEST
      complex(b8), allocatable :: cafrbuf(:,:,:)
#endif
      integer, allocatable :: pos(:)
      integer i,j,iz,dnz,iy,a,ia,a1,a2,c
      real vmax

      if(num_al_i(0) .gt. 0) then

      t3_comm = t3_comm - MPI_Wtime() 

      allocate (sendbuf(jjsz(0)*num_al_i(0),nv,0:jproc-1))
      allocate (recvbuf(jjsz(0)*num_al_i(0),nv,0:jproc-1))
#ifdef CAF_TEST
      allocate (cafrbuf(jjsz(0)*num_al_i(0),nv,0:jproc-1))
#endif

      sendbuf = 0
      recvbuf = 0
#ifdef CAF_TEST
      cafrbuf = 0
#endif

! This loop does the following:
! - for each x, 
! - it packs the send buffer while rearranging the data 
!   so Z has stride 1, using loop blocking for cache optimization

      do j=1,nv
      
      if(num_al .gt. 0) then

         do i=0,jproc-1

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
      endif

      enddo
      
#ifdef TIMERS	
      t_alltoall = t_alltoall - MPI_Wtime()
#endif	

#if defined (CAF_TEST) || !defined (CAF)
      call mpi_alltoall(sendbuf,jjsz(0)*num_al_i(0)*b8*2*nv,mpi_byte,
     &      recvbuf,jjsz(0)*num_al_i(0)*b8*2*nv,mpi_byte,mpi_comm_col,ierr)
#endif
#ifdef CAF
      call compi_alltoall(sendbuf,recvbuf,jjsz(0)*num_al_i(0)*b8*2*nv,
     & mpi_comm_col)
#endif
#ifdef CAF_TEST
      call compi_alltoall(sendbuf,cafrbuf,jjsz(0)*num_al_i(0)*b8*2*nv,
     & mpi_comm_col)
      if (any(cafrbuf.ne.recvbuf)) then
        print *,'kxcomm1_cyl2sq taskid correctness failed',taskid
      else
        print *,'kxcomm1_cyl2sq taskid correctness passed',taskid
      endif
#endif

#ifdef TIMERS	
      t_alltoall = t_alltoall + MPI_Wtime()
#endif	

      allocate(pos(0:jproc-1))

      do j=1,nv

         do i=0,jproc-1
            pos(i) = 1
         enddo

         do y=1,yjsz
            i = 0
            n=0
            do x=1,max_al_x
               c = cut_z(x)
               a1 = nzhp - c -1
               a2 = nzhp + c +1
               do z=1,a1
                  dest(x,z,y,j) = recvbuf(pos(i),j,i) 
                  pos(i) = pos(i)+1
                  n = n+1
                  if(n .eq. num_al_i(i)) then
                     n = 0
                     i = i+1
                  endif
               enddo
               
               do z=a1+1,a2-1
                  dest(x,z,y,j) = 0.0
               enddo
               
               do z=a2,nz
                  dest(x,z,y,j) = recvbuf(pos(i),j,i) 
                  pos(i) = pos(i)+1
                  n = n+1
                  if(n .eq. num_al_i(i)) then
                     n = 0
                     i = i+1
                  endif
               enddo
               
            enddo
            do x=max_al_x+1,xisz
               do z=1,nz
                  dest(x,z,y,j) = 0.0
               enddo
            enddo
         enddo
      enddo
      deallocate (pos)
         
      t3_comm = t3_comm + MPI_Wtime() 
      i3_comm = i3_comm + 1

      deallocate(sendbuf)
      deallocate (recvbuf)
#ifdef CAF_TEST
      deallocate (cafrbuf)
#endif

      else
! if num_al_i(0) = 0
         do j=1,nv
            do y=1,yjsz
               do x=1,xisz
                  do z=1,nz
c                    dest(z,x,y,j) = 0.0
                     dest(x,z,y,j) = 0.0
                  enddo
               enddo
            enddo
         enddo

      endif

      return
      end

!-----------------------------------------------------------------
! Square pencils version (used in inhomogeneous version)
!
! Starting with Y-pencils, wavenumber-space in all three dimensions,
! padded in Y:
! - Inverse-Transform in Y, 
! - transpose to get data in z-pencils with stride-1 in Z
! - and optionally expand in Z by filling the center with zeros
! 
!-----------------------------------------------------------------
      subroutine kxcomm1_trans_pen (source,dest)
      use com
      use timers_comm
      implicit none
      include 'intvars'
#ifdef FFTW
#include "fft_stuff.f"
#endif

      complex(b8) source(nypad,zjsz,xisz)
      complex(b8) dest(nzpad,xisz,yjsz)

      complex(b8), allocatable :: sendbuf(:,:),recvbuf(:,:)
#ifdef CAF_TEST
      complex(b8), allocatable :: cafrbuf(:,:)
#endif
      integer position,pos0,pos1
      integer, allocatable :: pos(:)
      integer i,j,n,iz,dnz,iy

      t3_comm = t3_comm - MPI_Wtime() 


! This case is for using MPI_Alltoall 
#ifdef USE_EVEN

! Allocate these buffers in the beginning ? 

      allocate (recvbuf(JrCntMax/(b8*2),0:jproc-1))
#ifdef CAF_TEST
      allocate (cafrbuf(JrCntMax/(b8*2),0:jproc-1))
#endif
      allocate(sendbuf(JrCntMax/(b8*2),0:jproc-1))
      if(JrCntUneven) then
         allocate(pos(0:jproc-1))
      endif

! This loop does the following:
! - for each x, it executes zjsz transforms in Y, to reduce memory references
! - it packs the send buffer while rearranging the data 
!   so Z has stride 1, using loop blocking for cache optimization

      do x=1,xisz

         do z=1,zjsz
            source(nyhp,z,x)=0.
         enddo

! What about nzhp?

#ifdef FFTW
         call fftw_execute(plan2_p1,source(1,1,x),source(1,1,x))
#elif ESSL
         call escfft (source(1,1,x),1,nypad,nypad,zjsz,1)
#endif
         
         do i=0,jproc-1
           
            pos0 = (x-1)*zjsz 

            do y=jjst(i),jjen(i),NB2_Y
               y2 = min(y+NB2_Y-1,jjen(i))

               do z=1,zjsz,NB2_Z
                  z2 = min(z+NB2_Z-1,zjsz)

                  pos1 = pos0+z
                  do iy=y,y2
                     position = pos1 
                     do iz=z,z2
                        sendbuf(position,i) = source(iy,iz,x)
                        position = position +1
                     enddo
                     pos1 = pos1 + zjsz * xisz
                  enddo
               enddo
               pos0 = pos0 + xisz*zjsz*NB2_Y
            enddo
         enddo
         if(JrCntUneven) then
            do i=0,jproc-1
               pos(i) = position
            enddo
         endif
      enddo

! Fill the buffer with zeros (not necessary if preallocating)
      
      if(JrCntUneven) then
         do i=0,jproc-1
            position = pos(i)
            n = JrCntMax/(b8*2) - jjsz(i)*zjsz*xisz
            do j=1,n
               sendbuf(position,i) = 0.0
               position = position+1
            enddo
         enddo
         deallocate(pos)
      endif


! Exchange the data (send the buffers)

#ifdef TIMERS	
      t_alltoall = t_alltoall - MPI_Wtime()
#endif	

#if defined (CAF_TEST) || !defined (CAF)
      call mpi_alltoall(sendbuf,JrCntMax, mpi_byte,
     &     recvbuf,JrCntMax, mpi_byte,mpi_comm_col,ierr)
#endif
#ifdef CAF
      call compi_alltoall(sendbuf,recvbuf,JrCntMax,
     & mpi_comm_col)
#endif
#ifdef CAF_TEST
      call compi_alltoall(sendbuf,cafrbuf,JrCntMax,
     & mpi_comm_col)
      if (any(cafrbuf.ne.recvbuf)) then
        print *,'kxcomm1_trans_pen taskid correctness failed',taskid
      else
        print *,'kxcomm1_trans_pen taskid correctness passed',taskid
      endif
#endif

#ifdef TIMERS	
      t_alltoall = t_alltoall + MPI_Wtime()
#endif	

      deallocate(sendbuf)

! Unpack the receive buffer and expand the Z dimension by filling 
! the center with (nzpad - nz) zeroes

      dnz = nzpad - nz
      do i=0,jproc-1
         position = 1
! If clearly in the first half of nz
         if(kjen(i) .le. nzh) then
            do y=1,yjsz
               do x=1,xisz
                  do z=kjst(i),kjen(i)
                     dest(z,x,y) = recvbuf(position,i)
                     position = position +1
                  enddo
               enddo
            enddo

! If clearly in the second half of nz
         else if (kjst(i) .ge. nzh+1) then
            do y=1,yjsz
               do x=1,xisz
                  do z=kjst(i)+dnz,kjen(i)+dnz
                     dest(z,x,y) = recvbuf(position,i)
                     position = position +1
                  enddo
               enddo
            enddo         

! If spanning the first and second half of nz (i.e. iproc is odd)  
         else
            do y=1,yjsz
               do x=1,xisz
                  do z=kjst(i),nzh
                     dest(z,x,y) = recvbuf(position,i)
                     position = position +1
                  enddo
                  do z=nzpad-nzh+2,kjen(i)+dnz
                     dest(z,x,y) = recvbuf(position,i)
                     position = position +1
                  enddo
               enddo
            enddo
         endif

! Fill center with zeros
         do y=1,yjsz
            do x=1,xisz
               do z=nzhp,nzpad-nzh+1
                  dest(z,x,y) = 0.0
               enddo
            enddo
         enddo

      enddo

#else
! Not USE EVEN - using MPI_Alltoallv

      allocate (sendbuf(xisz*zjsz*nypad,1))
      allocate (recvbuf(xisz*yjsz*nz,1))
#ifdef CAF_TEST
      allocate (cafrbuf(xisz*yjsz*nz,1))
#endif

! This loop does the following:
! - for each x, it executes zjsz transforms in Y, to reduce memory references
! - it packs the send buffer while rearranging the data 
!   so Z has stride 1, using loop blocking for cache optimization

      do x=1,xisz

         do z=1,zjsz
            source(nyhp,z,x)=0.
         enddo

! What about nzhp?

#ifdef FFTW
         call fftw_execute(plan2_p1,source(1,1,x),source(1,1,x))
#elif ESSL
         call escfft (source(1,1,x),1,nypad,nypad,zjsz,1)
#endif

         pos0 = (x-1)*zjsz 

         do y=1,nypad,NB2_Y
            y2 = min(y+NB2_Y-1,nypad)

            do z=1,zjsz,NB2_Z
               z2 = min(z+NB2_Z-1,zjsz)
               pos1 = pos0 +z
               
               do iy=y,y2
                  position = pos1 
                  do iz=z,z2
                     sendbuf(position,1) = source(iy,iz,x)
                     position = position +1
                  enddo
                  pos1 = pos1 + zjsz * xisz
               enddo
            enddo
            pos0 = pos0 + xisz*zjsz*NB2_Y
         enddo
      enddo
         
#ifdef TIMERS	
	t_alltoall = t_alltoall - MPI_Wtime()
#endif	

#if defined (CAF_TEST) || !defined (CAF)
      call mpi_alltoallv(sendbuf,JrSndCnts, JrSndStrt,mpi_byte,
     &     recvbuf,JrRcvCnts, JrRcvStrt,mpi_byte,mpi_comm_col,ierr)
#endif
#ifdef CAF
      call compi_alltoallv(sendbuf,JrSndCnts,JrSndStrt,JrCntMax,
     &                     recvbuf,JrRcvCnts,JrRcvStrt,mpi_comm_col)
#endif
#ifdef CAF_TEST
      call compi_alltoallv(sendbuf,JrSndCnts,JrSndStrt,JrCntMax,
     &                     cafrbuf,JrRcvCnts,JrRcvStrt,mpi_comm_col)
      if (any(cafrbuf.ne.recvbuf)) then
        print *,'kxcomm1_trans_pen taskid correctness failed',taskid
      else
        print *,'kxcomm1_trans_pen taskid correctness passed',taskid
      endif
#endif

#ifdef TIMERS	
	t_alltoall = t_alltoall + MPI_Wtime()
#endif	

      deallocate(sendbuf)

! Unpack the receive buffer and expand the Z dimension by filling 
! the center with (nzpad - nz) zeroes

      dnz = nzpad - nz
      position = 1
      do i=0,jproc-1

! If clearly in the first half of nz
         if(kjen(i) .le. nzh) then
            do y=1,yjsz
               do x=1,xisz
                  do z=kjst(i),kjen(i)
                     dest(z,x,y) = recvbuf(position,1)
                     position = position +1
                  enddo
               enddo
            enddo

! If clearly in the second half of nz
         else if (kjst(i) .gt. nzh) then
            do y=1,yjsz
               do x=1,xisz
                  do z=kjst(i)+dnz,kjen(i)+dnz
                     dest(z,x,y) = recvbuf(position,1)
                     position = position +1
                  enddo
               enddo
            enddo         

! If spanning the first and second half of nz (i.e. iproc is odd)  
         else
            do y=1,yjsz
               do x=1,xisz
                  do z=kjst(i),nzh
                     dest(z,x,y) = recvbuf(position,1)
                     position = position +1
                  enddo
                  do z=nzpad-nzh+1,kjen(i)+dnz
                     dest(z,x,y) = recvbuf(position,1)
                     position = position +1
                  enddo
               enddo
            enddo
         endif

      enddo

! Fill center with zeros
      do y=1,yjsz
         do x=1,xisz
            do z=nzhp,nzpad-nzh+1
               dest(z,x,y) = 0.0
            enddo
         enddo
      enddo
      
#endif
      
      t3_comm = t3_comm + MPI_Wtime() 
      i3_comm = i3_comm + 1

      deallocate (recvbuf)
#ifdef CAF_TEST
      deallocate (cafrbuf)
#endif

      return
      end
