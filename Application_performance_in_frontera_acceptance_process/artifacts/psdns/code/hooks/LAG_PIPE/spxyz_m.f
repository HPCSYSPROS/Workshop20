      subroutine spxyz_m (f,bsy,nv)
c
#ifdef LAG
c
c Revised by PKY: Jul 11, 2011: replaced global indices
c (which seems to cause seg-faults sometimes) by local indices
c
c 2D code
c
	use com, only:istep,nsteps,jstep
	use mpilag
	use lag_timers
	use compart, only: iout
	implicit none
c

      real(b8)  f(nx,zisz,yjsz,nv)
      real(p8), allocatable :: bsz(:,:,:,:),buf1(:),buf2(:),buf3(:)
      real(p8)  bsy(bxisize,nby,bzjsize,nv)


        integer ithr, omp_get_thread_num,maxcount
	integer, allocatable :: bp1(:),bpi(:), ipp1(:),ippi(:)
	integer, allocatable :: jpp1(:),jppi(:)
	integer, allocatable :: zpp1(:),zppi(:)

	real*8 rtime0,rtime1,rtime2
	real*8 tcpu,avcpu,cpumin,cpumax
c
        real*8 cpu_alloc, cpu_buf1, cpu_buf2, cpu_buf3
        real*8 cpu_splx, cpu_sply, cpu_splz, cpu_buf1z
        real*8 cpu_comm1, cpu_comm2, cpu_spxyzin

	integer iv,nv,jj,ii,iz
	real mem
c
	integer icall
	data icall/0/
c
	integer, allocatable :: pos0(:),pos2(:),pos3(:)
c
      integer z,y,pos1,pos,z1,z2,jjy
      integer x,i,j
c
	integer xp,yp,zp
c
	real(p8), allocatable :: bsx(:,:,:,:)
c
	
	jj=1
	if (iout.gt.0) jj=2
c
	icall=icall+1
c	if (taskid.eq.0) write(6,*) 'spxyz: icall=',icall
c
	rtime0=MPI_WTIME()
	rtime1=MPI_WTIME()
	
        allocate (bp1(0:num_thr-1),bpi(0:num_thr-1))
        allocate (ipp1(0:num_thr-1),ippi(0:num_thr-1))
        allocate (jpp1(0:num_thr-1),jppi(0:num_thr-1))
        allocate (zpp1(0:num_thr-1),zppi(0:num_thr-1))
	allocate (pos0(0:iproc-1))
	allocate (pos2(0:jproc-1))
	allocate (pos3(0:jproc-1))

      allocate (bsx(nbx,zisz,yjsz,nv))
c
      allocate(buf1(nbx*zisz*nv),stat=ierr)
	if (ierr.gt.0) then
	if (taskid.eq.0) then
	mem=float(nbx)*zisz*nv*4/1048576
	write (6,*) 'spxyz_m: cannot allocate buf1, needs', mem, ' MB/core'
	end if
	stop
	end if
c
c     allocate(bsz(nv,bxistart:bxiend,nbz,yjsz),stat=ierr)
      allocate(bsz(nv,bxisize,nbz,yjsz),stat=ierr)
	bsz=0.
	if (ierr.gt.0) then
	if (taskid.eq.0) then
	mem=float(nv)*(bxiend-bxistart+1)*nbz*(yjen-yjst+1)*4/1048576
	write (6,*) 'spxyz_m: cannot allocate bsz, needs', mem, ' MB/core'
	end if
	stop
	end if
c
	cpu_buf1 = 0.
	cpu_comm1 = 0.

c form splines in x and z within each slab
c
        call divide_thr (yjsz, bp1, bpi, 'spxyz_m: yjsz')
        call divide_thr (iproc, ipp1, ippi, 'spxyz_m: iproc')
	ipp1(:)=ipp1(:)-1
        call divide_thr (jproc, jpp1, jppi, 'spxyz_m: jproc')
	jpp1(:)=jpp1(:)-1
        call divide_thr (bzjsize, zpp1, zppi, 'spxyz_m: bzjsize')
c
	pos0(0)=1
        do i=1,iproc-1
	ii=mymap(i)
	pos0(i)=pos0(i-1)+bxisz(i-1)*zisz*nv
        enddo


	ithr=0
!$OMP PARALLEL private (ithr,pos1)
#ifdef OPENMP
        ithr = omp_get_thread_num()
#endif

!$OMP MASTER
	if(jstep.ne.1) then
	cpu_spxyz(1,jj)=cpu_spxyz(1,jj) + MPI_WTIME()-rtime1
	endif
	rtime1=MPI_WTIME()
!$OMP END MASTER
c
	do 10 iv=1,nv
	do 10 y=bp1(ithr),bp1(ithr)+bpi(ithr)-1
          call splx (f(1,1,y,iv),bsx(1,1,y,iv))
 10   continue

!$OMP BARRIER
	if (ithr.eq.0) call splx_minmax ('after do 10',bsx,nv)

!$OMP MASTER
	if(jstep.ne.1) then
	cpu_spxyz(2,jj)=cpu_spxyz(2,jj) + MPI_WTIME()-rtime1
	endif
        rtime1 = MPI_WTIME()
!$OMP END MASTER

c pack and send
c


      do 20 yp=1,yjsz
!$OMP MASTER
	rtime1=MPI_WTIME()
!$OMP END MASTER
         pos1=1
c         do i=0,iproc-1
	do i=ipp1(ithr),ipp1(ithr)+ippi(ithr)-1
	ii=mymap(i)
         pos1=pos0(ii)
            do zp=1,zisz
               do x=bxist(ii),bxien(ii)
		do iv=1,nv
                  buf1(pos1) = bsx(x,zp,yp,iv)
                  pos1 = pos1+1
	        end do
               enddo
            enddo
         enddo

!$OMP BARRIER
!$OMP MASTER

	if(jstep.ne.1) then
	cpu_spxyz(3,jj)=cpu_spxyz(3,jj) + MPI_WTIME()-rtime1
	endif
	rtime1=MPI_WTIME()
c
	bcnts_zi(:)=bcnts_zi(:)*nv
	bdisp_zi(:)=bdisp_zi(:)*nv
	bcnts_xi(:)=bcnts_xi(:)*nv
	bdisp_xi(:)=bdisp_xi(:)*nv


! 12/15/12 : D. Buaria: Use of pure mpi seems to be faster and more
! reliable compared to CAF for this alltoallv
! But for the 2nd alltoallv CAF is faster
#ifdef CAF 

#ifdef CF_LAG
        maxcount=maxval(bcnts_zi)
	call MPI_BARRIER (MPI_COMM_WORLD, mpierr)

         call compi2_alltoallv(buf1,bcnts_zi,bdisp_zi,maxcount,
     &        bsz(1,1,1,yp),bcnts_xi,bdisp_xi,
     &        mpi_comm_row)

#else
         call mpi_alltoallv(buf1,bcnts_zi,bdisp_zi,pptype,
     &        bsz(1,1,1,yp),bcnts_xi,bdisp_xi,pptype,
     &        mpi_comm_row,ierr)

#endif

#else
         call mpi_alltoallv(buf1,bcnts_zi,bdisp_zi,pptype,
     &        bsz(1,1,1,yp),bcnts_xi,bdisp_xi,pptype,
     &        mpi_comm_row,ierr)
#endif


	bcnts_zi(:)=bcnts_zi(:)/nv
	bdisp_zi(:)=bdisp_zi(:)/nv
	bcnts_xi(:)=bcnts_xi(:)/nv
	bdisp_xi(:)=bdisp_xi(:)/nv

	if(jstep.ne.1) then
	cpu_spxyz(4,jj)=cpu_spxyz(4,jj) + MPI_WTIME()-rtime1
	endif

!$OMP END MASTER
!$OMP BARRIER

 20	continue

!$OMP MASTER
	rtime1=MPI_WTIME()
!$OMP END MASTER
	
!$OMP BARRIER
	if (ithr.eq.0) call splz_minmax ('after do 20',bsz,nv)
!$OMP BARRIER

	do yp=bp1(ithr),bp1(ithr)+bpi(ithr)-1
         call splz_m (bsz(1,1,1,yp),nv)
	end do

!$OMP BARRIER
	if (ithr.eq.0) call splz_minmax ('after splz',bsz,nv)
!$OMP BARRIER

!$OMP MASTER
	if(jstep.ne.1) then
	cpu_spxyz(5,jj)=cpu_spxyz(5,jj) + MPI_WTIME()-rtime1
	endif
!$OMP END MASTER

!$OMP END PARALLEL

	rtime1=MPI_WTIME()

      deallocate(buf1)
c
	deallocate (bsx)


      allocate(buf2(bxisize*yjsz*nbz*nv),stat=ierr)
	if (ierr.gt.0) then
	if (taskid.eq.0) then
	mem=float(nv)*(bxisize*yjsz*nbz)*4/1048576
	write (6,*) 'spxyz_m: cannot allocate buf2, needs', mem, ' MB/core'
	end if
	stop
	end if

	pos2(0)=1
        do j=1,jproc-1
	pos2(j)=pos2(j-1)+bzjsz(j-1)*yjsz*bxisize*nv
        enddo

	pos3(0)=1
        do j=1,jproc-1
	pos3(j)=pos3(j-1)+bzjsize*jjsz(j)*bxisize*nv
        enddo

	if(jstep.ne.1) then
	cpu_spxyz(1,jj)=cpu_spxyz(1,jj) + MPI_WTIME()-rtime1
	endif

	ithr=0
!$OMP PARALLEL private (ithr,pos1,z1,z2,jjy)
#ifdef OPENMP
        ithr = omp_get_thread_num()
#endif

!$OMP MASTER
	rtime1=MPI_WTIME()
!$OMP END MASTER

      pos1=1
c     do j=0,jproc-1
	do j=jpp1(ithr),jpp1(ithr)+jppi(ithr)-1
	pos1=pos2(j)
         do z=bzjst(j),bzjen(j)
            do yp=1,yjsz
               do x=1,bxisize
		do iv=1,nv
                  buf2(pos1) = bsz(iv,x,z,yp)
                  pos1 = pos1+1
               enddo
               enddo
            enddo
         enddo
      enddo

!$OMP MASTER
	if(jstep.ne.1) then
	cpu_spxyz(6,jj)=cpu_spxyz(6,jj) + MPI_WTIME()-rtime1
	endif
!$OMP END MASTER

!$OMP BARRIER
!$OMP MASTER

	rtime1=MPI_WTIME()

      deallocate(bsz)
      allocate(buf3(bxisize*nby*bzjsize*nv),stat=ierr)
	if (ierr.gt.0) then
	if (taskid.eq.0) then
	mem=float(nv)*(bxisize*nby*bzjsize)*4/1048576
	write (6,*) 'spxyz_m: cannot allocate buf3, needs', mem, ' MB/core'
	end if
	stop
	end if
c
	if(jstep.ne.1) then
	cpu_spxyz(1,jj)=cpu_spxyz(1,jj) + MPI_WTIME()-rtime1
	endif
	rtime1=MPI_WTIME()
c
	bcnts_yj(:)=bcnts_yj(:)*nv
	bdisp_yj(:)=bdisp_yj(:)*nv
	bcnts_zj(:)=bcnts_zj(:)*nv
	bdisp_zj(:)=bdisp_zj(:)*nv
c

#ifdef CAF

#ifdef CF_LAG

        maxcount=maxval(bcnts_yj)

      call compi2_alltoallv(buf2,bcnts_yj,bdisp_yj,maxcount,
     &     buf3,bcnts_zj,bdisp_zj,mpi_comm_col)

#else

      call mpi_alltoallv(buf2,bcnts_yj,bdisp_yj,pptype,
     &     buf3,bcnts_zj,bdisp_zj,pptype,
     &     mpi_comm_col,ierr)

#endif

#else

      call mpi_alltoallv(buf2,bcnts_yj,bdisp_yj,pptype,
     &     buf3,bcnts_zj,bdisp_zj,pptype,
     &     mpi_comm_col,ierr)

#endif

	bcnts_yj(:)=bcnts_yj(:)/nv
	bdisp_yj(:)=bdisp_yj(:)/nv
	bcnts_zj(:)=bcnts_zj(:)/nv
	bdisp_zj(:)=bdisp_zj(:)/nv

      deallocate(buf2)

	if(jstep.ne.1) then
	cpu_spxyz(7,jj)=cpu_spxyz(7,jj) + MPI_WTIME()-rtime1
	endif
	
!$OMP END MASTER
!$OMP BARRIER


!$OMP MASTER
	rtime1=MPI_WTIME()
!$OMP END MASTER

	z1=zpp1(ithr)
	z2=zpp1(ithr)+zppi(ithr)-1
c
#ifdef TEST
      pos1 = 1
       do j=0,jproc-1
c	do j=jpp1(ithr),jpp1(ithr)+jppi(ithr)-1
	pos1=pos3(j)
c         do zp=1,bzjsize
c	do zp=zpp1(ithr),zpp1(ithr)+zppi(ithr)-1
	do zp=z1,z2
            do y=jjst(j),jjen(j)
               do xp=1,bxisize
		do iv=1,nv
                  bsy(xp,y,zp,iv) = buf3(pos1)
                  pos1 = pos1 + 1
	        end do
               enddo
            enddo
         enddo
      enddo

#else
c this loop based on effort by D. Buaria
      do j=0,jproc-1
         jjy = (jjen(j) - jjst(j) + 1)*bxisize*nv
 	do zp=zpp1(ithr),zpp1(ithr)+zppi(ithr)-1
            do y=jjst(j),jjen(j)
               do xp=1,bxisize
                do iv=1,nv
                  pos1 = (j*bzjsize + zp-1)*jjy
     &                    + (y - jjst(j))*bxisize*nv + (xp-1)*nv + iv
                  bsy(xp,y,zp,iv) = buf3(pos1)
                end do
               enddo
            enddo
         enddo
      enddo
#endif

!$OMP MASTER
	if(jstep.ne.1) then
	cpu_spxyz(8,jj)=cpu_spxyz(8,jj) + MPI_WTIME()-rtime1
	endif
	rtime1 = MPI_WTIME()
!$OMP END MASTER

      do iv=1,nv
c         do zp=1,bzjsize
	do zp=zpp1(ithr),zpp1(ithr)+zppi(ithr)-1
            call sply (bsy(1,1,zp,iv))
      enddo
      enddo

!$OMP MASTER
	if(jstep.ne.1) then
	cpu_spxyz(9,jj)=cpu_spxyz(9,jj) + MPI_WTIME()-rtime1
	endif
!$OMP END MASTER

!$OMP END PARALLEL

	do iv=1,nv
c	call bsminmax ('exit spxyz_m',bsy(1,1,1,iv),1)
c	call bsminmaxv ('exit spxyz_m',bsy(1,1,1,iv),1)
	end do
      deallocate(buf3)
c
        deallocate (bp1,bpi,ipp1,ippi,jpp1,jppi,zpp1,zppi)
	deallocate (pos0,pos2,pos3)

	if(jstep.ne.1) then
	cpu_spxyz(10,jj)=cpu_spxyz(10,jj) + MPI_WTIME()-rtime0
	endif

c	if (taskid.eq.0.and.icall.eq.1) write (6,*) ' exit spxyz_m'

      return
c
 99   call MPI_ABORT (MPI_COMM_WORLD,ierr)

 703    format('spxyz (',a7,'):  min/ave/max=',3f9.3,'   secs.')
c
#endif
c
      end
