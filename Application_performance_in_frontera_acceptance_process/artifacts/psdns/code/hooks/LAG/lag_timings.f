	subroutine lag_timings
c
       	use comp
	use timers_comm
	use timers_comp
	use timers_io
	use timers_rkstep
	use timers_tran
	use lag_timers
	use compart, only:nop,nop2,nom

        implicit none

	integer itask,i,i2,jj
	real(8) sum1,comm2,sum2,sum3,sum4,sum5,sum6,sum7,sum8,sum9
        real(8), allocatable :: t_rks_all(:,:),t_itrans_all(:,:)
        real(8), allocatable :: t_spxyz_all(:,:)
        real(8), allocatable :: t_partic_all(:)
        real(8), allocatable :: t_partic2_all(:)
        real(8), allocatable :: t_partsp_all(:,:)
        real(8), allocatable :: t_rkpart_all(:,:)
        real(8), allocatable :: t_intbf_all(:,:,:)
        real(8), allocatable :: t_spcal_all(:,:,:)
        real(8), allocatable :: t_lagout_all(:)
c
	integer itask1,itask2,j
	character*20 string
	character*12 dir
c
	data dir/'lag_timings/'/
c
c detailed instrumentation timings, added by PKY, 2/3/2012
c
	if (ncpusteps.eq.0) return
c
	allocate (t_rks_all(10,0:numtasks-1))
	allocate (t_spxyz_all(0:10,0:numtasks-1))
	allocate (t_partic_all(0:numtasks-1))
	allocate (t_partic2_all(0:numtasks-1))
	allocate (t_partsp_all(0:5,0:numtasks-1))
	allocate (t_rkpart_all(2,0:numtasks-1))
	allocate (t_intbf_all(4,3,0:numtasks-1))
	allocate (t_spcal_all(4,3,0:numtasks-1))
	allocate (t_lagout_all(0:numtasks-1))

	cpu_spxyz(:,1)=cpu_spxyz(:,1)/(nspcal/2)
        call  MPI_ALLGATHER (cpu_spxyz(0,1),11,MPI_DOUBLE_PRECISION,
     1                       t_spxyz_all,11,MPI_DOUBLE_PRECISION,
     1                       MPI_COMM_WORLD,mpierr)
	cpu_partic2(1)=cpu_partic2(1)/(nspcal/2)
        call  MPI_ALLGATHER (cpu_partic2(1),1,MPI_DOUBLE_PRECISION,
     1                       t_partic2_all,1,MPI_DOUBLE_PRECISION,
     1                       MPI_COMM_WORLD,mpierr)
	cpu_partsp(:,1)=cpu_partsp(:,1)/(nspcal/2)
        call  MPI_ALLGATHER (cpu_partsp(0,1),6,MPI_DOUBLE_PRECISION,
     1                       t_partsp_all,6,MPI_DOUBLE_PRECISION,
     1                       MPI_COMM_WORLD,mpierr)
	cpu_rkpart(:,1)=cpu_rkpart(:,1)/(nspcal/2)
        call  MPI_ALLGATHER (cpu_rkpart(1,1),2,MPI_DOUBLE_PRECISION,
     1                       t_rkpart_all,2,MPI_DOUBLE_PRECISION,
     1                       MPI_COMM_WORLD,mpierr)
	cpu_intbf(:,:,1)=cpu_intbf(:,:,1)/nintbf
        call  MPI_ALLGATHER (cpu_intbf(1,1,1),4*3,MPI_DOUBLE_PRECISION,
     1                       t_intbf_all,4*3,MPI_DOUBLE_PRECISION,
     1                       MPI_COMM_WORLD,mpierr)
	cpu_spcal(:,:,1)=cpu_spcal(:,:,1)/nspcal
	if (taskid.eq.0) write (6,*) 'nintbf,nspcal=',nintbf,nspcal
        call  MPI_ALLGATHER (cpu_spcal(1,1,1),4*3,MPI_DOUBLE_PRECISION,
     1                       t_spcal_all,4*3,MPI_DOUBLE_PRECISION,
     1                       MPI_COMM_WORLD,mpierr)
!	cpu_lagout=cpu_lagout/ncpusteps
        call  MPI_ALLGATHER (cpu_lagout,1,MPI_DOUBLE_PRECISION,
     1                       t_lagout_all,1,MPI_DOUBLE_PRECISION,
     1                       MPI_COMM_WORLD,mpierr)

	t_rks(:,1)=t_rks(:,1)/ncpusteps
        call  MPI_ALLGATHER (t_rks(1,1),10,MPI_DOUBLE_PRECISION,
     1                       t_rks_all,10,MPI_DOUBLE_PRECISION,
     1                       MPI_COMM_WORLD,mpierr)
c
	if (taskid.eq.0) then
c
c

	itask1=2*iproc-1
	itask2=numtasks-2*iproc

	open (7,file=dir//'cpu_rk_all')
      write (7,711) nx,ny,nz,nc,iproc, jproc
 711  format ('Pencils code, nx,ny,nz,nc=',3i6,i3,'  iproc, jproc=',i4,
     1        ' x ',i4)
#ifndef OPENMP
	num_thr=0
#endif
#ifdef DOUBLE_PREC
        write (7,"('num_thr=',i3,'  double prec,  rkmethod=',i3)") num_thr,
     1    rkmethod
#else
        write (7,"('num_thr=',i3,'  single prec,  rkmethod=',i3)") num_thr,
     1    rkmethod
#endif
#ifndef OPENMP
	num_thr=1
#endif

	write (7,"('Detailed breakdown of CPU costs per time step')")
	write (7,"('averaged over', i4,'  non-output steps')") ncpusteps
	write (7,"('taskid ipid jpid           itransform realspace',
     1             ' transform           overall')")
        do itask=0,numtasks-1
	if (itask.le.itask1.or.itask.ge.itask2) then
	sum1=0.
	do i=1,5
	sum1=sum1+t_rks_all(i,itask)
	end do
        write (7,630) itask,ipid_all(itask),jpid_all(itask),
     1             (t_rks_all(i,itask),i=1,6)
        if (mod(itask+1,iproc).eq.0) then
        write (7,"('--------------------------------------------')")
        end if
	end if
	end do
 630    format (i6,2i5,1x,1p,6e10.3)
	close (7)
c

	open(unit=7,file=dir//'t_spxyz2_all')
      write (7,711) nx,ny,nz,nc,iproc, jproc
        write (7,"('num_thr=',i3,'  single prec,  rkmethod=',i3)") num_thr,
     1    rkmethod
	write(7,"('nop = ',i9,'  nop2 =',i9)")nop,nop2
        write (7,"('--------------------------------------------')")
	write(7,*)'taskid ipid jpid    splx      buf1      comm1     splz      buf2     comm2      buf3      sply',
     1        '     sum       overall '
	do itask=0,numtasks-1
	if (itask.le.iproc-1.or.itask.ge.numtasks-iproc) then
	write(7,634)itask,ipid_all(itask),jpid_all(itask),(t_spxyz_all(i,itask),
     1             i=2,9),sum(t_spxyz_all(1:9,itask)),t_spxyz_all(0,itask)
	end if
	enddo
 634    format (i6,2i5,1x,1p,10e10.3)
	close(7)

	open(unit=7,file=dir//'t_spxyz_all')
      write (7,711) nx,ny,nz,nc,iproc, jproc
        write (7,"('num_thr=',i3,'  single prec,  rkmethod=',i3)") num_thr,
     1    rkmethod
	write(7,"('nop = ',i9,'  nop2 =',i9)")nop,nop2
        write (7,"('--------------------------------------------')")
	write(7,*)'taskid ipid jpid    spls      bufs      comms     sum     overall '
	do itask=0,numtasks-1
	if (itask.le.iproc-1.or.itask.ge.numtasks-iproc) then
	sum1=t_spxyz_all(2,itask)+t_spxyz_all(5,itask)+t_spxyz_all(9,itask)
	sum2=t_spxyz_all(3,itask)+t_spxyz_all(6,itask)+t_spxyz_all(8,itask)
	sum3=t_spxyz_all(4,itask)+t_spxyz_all(7,itask)
	write(7,639)itask,ipid_all(itask),jpid_all(itask),sum1,
     1              sum2,sum3,sum1+sum2+sum3,t_spxyz_all(0,itask)         
	end if
	enddo
 639    format (i6,2i5,1x,1p,5e10.3)
	close(7)

c
	open(unit=7,file=dir//'t_partsp_all')
      write (7,711) nx,ny,nz,nc,iproc, jproc
        write (7,"('num_thr=',i3,'  single prec,  rkmethod=',i3)") num_thr,
     1    rkmethod
	write(7,"('nop = ',i9,'  nop2 =',i9)")nop,nop2
        write (7,"('--------------------------------------------')")
	write(7,*)'taskid ipid jpid  rkadvos1  rkpart1   rkadvos2  rkpart0     sum     overall    partic '
	do itask=0,numtasks-1
	if (itask.le.iproc-1.or.itask.ge.numtasks-iproc) then
	write(7,632)itask,ipid_all(itask),jpid_all(itask),
     1                (t_partsp_all(i,itask),i=2,5), sum(t_partsp_all(1:5,itask)),
     1                t_partsp_all(0,itask),t_partic2_all(itask)
	end if
	enddo
 632    format (i6,2i5,1x,1p,7e10.3)
	close(7)

	open(unit=7,file=dir//'t_intbf_all')
      write (7,711) nx,ny,nz,nc,iproc, jproc
        write (7,"('num_thr=',i3,'  single prec,  rkmethod=',i3)") num_thr,
     1    rkmethod
	write(7,"('nop = ',i9,'  nop2 =',i9,'  nom = ',i9)")nop,nop2,nom
        write (7,"('--------------------------------------------')")
	write(7,"('taskid ipid jpid    loop      comm     other      total    overall')")
	jj=3
	if (nom.eq.0) jj=2
	do j=1,jj
	write (7,"('--------------------------------------------')")
	if (j.eq.1) write (7,"('nop=',i11)") nop
	if (j.eq.2) write (7,"('nop2=',i10)") nop2
	if (j.eq.3) write (7,"('nom=',i11)") nom
	do itask=0,numtasks-1
	if (itask.le.iproc-1.or.itask.ge.numtasks-iproc) then
	write(7,635)itask,ipid_all(itask),jpid_all(itask),(t_intbf_all(i,j,itask),i=1,3),
     1             sum(t_intbf_all(1:3,j,itask)),t_intbf_all(4,j,itask)
	end if
	enddo
	end do
 635    format (i6,2i5,1x,1p,8e10.3)
	close(7)

	open(unit=7,file=dir//'t_spcal_all')
      write (7,711) nx,ny,nz,nc,iproc, jproc
        write (7,"('num_thr=',i3,'  single prec,  rkmethod=',i3)") num_thr,
     1    rkmethod
	write(7,"('nop = ',i9,'  nop2 =',i9,'  nom = ',i9)")nop,nop2,nom
        write (7,"('--------------------------------------------')")
	write(7,"('taskid ipid jpid    loop      comm     other      total    overall')")
	
	jj=3
	if (nom.eq.0) jj=2
	do j=1,jj
	write (7,"('--------------------------------------------')")
	if (j.eq.1) write (7,"('nop=',i11)") nop
	if (j.eq.2) write (7,"('nop2=',i10)") nop2
	if (j.eq.3) write (7,"('nom=',i11)") nom
	do itask=0,numtasks-1
	if (itask.le.iproc-1.or.itask.ge.numtasks-iproc) then
	write(7,645)itask,ipid_all(itask),jpid_all(itask),(t_spcal_all(i,j,itask),i=1,3),
     1             sum(t_spcal_all(1:3,j,itask)),t_spcal_all(4,j,itask)
	end if
	enddo
	end do
 645    format (i6,2i5,1x,1p,8e10.3)
	close(7)

	open(unit=7,file=dir//'t_rkpart_all')
      write (7,711) nx,ny,nz,nc,iproc, jproc
        write (7,"('num_thr=',i3,'  single prec,  rkmethod=',i3)") num_thr,
     1    rkmethod
	write(7,"('nop = ',i9,'  nop2 =',i9,'  nom = ',i9)")nop,nop2,nom
        write (7,"('--------------------------------------------')")
	write(7,"('taskid ipid jpid    spxyz  intbf+comm    total')")
	do itask=0,numtasks-1
	if (itask.le.iproc-1.or.itask.ge.numtasks-iproc) then
	write(7,633)itask,ipid_all(itask),jpid_all(itask),
     1             t_rkpart_all(1,itask),t_rkpart_all(2,itask),
     1             t_rkpart_all(1,itask)+t_rkpart_all(2,itask)
	end if
	enddo
 633    format (i6,2i5,1x,1p,3e10.3)
	close(7)


	end if

	return
	end
