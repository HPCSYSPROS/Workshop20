        subroutine write_timings 

       	use comp
	use timers_comm
	use timers_comp
	use timers_io
	use timers_rkstep
	use timers_tran

        implicit none

	integer itask,i,i2t,ii,jj
	real(8) sum,comm2
        real(8), allocatable :: t_rks_all(:,:),t_itrans_all(:,:)
        integer, allocatable :: numal_all(:)
        real(8), allocatable :: t_xkcomm1_all(:,:)
        real(8), allocatable :: t_xkcomm2_all(:,:)
        real(8), allocatable :: t_kxcomm1_all(:,:)
        real(8), allocatable :: t_kxcomm2_all(:,:)
        real(8), allocatable :: t_kxcomm2t_all(:,:)
        real(8), allocatable :: t_trans_all(:,:)
c
	integer itask1,itask2
	character*20 string
	character*6 filepos
c
c detailed instrumentation timings, added by PKY, 2/3/2012
c
	if (ncpusteps.eq.0) return
c
	if (taskid.eq.0) write (6,*) 'enter write_timings'
	allocate (t_rks_all(10,0:numtasks-1))
c       allocate (ipid_all(0:numtasks-1))
c       allocate (jpid_all(0:numtasks-1))
	allocate (t_itrans_all(4,0:numtasks-1))
	allocate (t_trans_all(4,0:numtasks-1))
	allocate (t_kxcomm1_all(4,0:numtasks-1))
	allocate (t_xkcomm1_all(4,0:numtasks-1))
	allocate (t_xkcomm2_all(4,0:numtasks-1))

c       call  MPI_ALLGATHER (ipid,1,MPI_INTEGER,ipid_all,1,MPI_INTEGER,
c    1                       MPI_COMM_WORLD,mpierr)
c       call  MPI_ALLGATHER (jpid,1,MPI_INTEGER,jpid_all,1,MPI_INTEGER,
c    1                       MPI_COMM_WORLD,mpierr)


	t_rks(:,1)=t_rks(:,1)/ncpusteps
        call  MPI_ALLGATHER (t_rks(1,1),10,MPI_DOUBLE_PRECISION,
     1                       t_rks_all,10,MPI_DOUBLE_PRECISION,
     1                       MPI_COMM_WORLD,mpierr)
	t_itrans(:,1)=t_itrans(:,1)/ncpusteps
        call  MPI_ALLGATHER (t_itrans,4,MPI_DOUBLE_PRECISION,
     1                       t_itrans_all,4,MPI_DOUBLE_PRECISION,
     1                       MPI_COMM_WORLD,mpierr)

	t_trans(:,1)=t_trans(:,1)/ncpusteps
        call  MPI_ALLGATHER (t_trans,4,MPI_DOUBLE_PRECISION,
     1                       t_trans_all,4,MPI_DOUBLE_PRECISION,
     1                       MPI_COMM_WORLD,mpierr)

	t_kxcomm1(:,1)=t_kxcomm1(:,1)/ncpusteps
        call  MPI_ALLGATHER (t_kxcomm1,4,MPI_DOUBLE_PRECISION,
     1                       t_kxcomm1_all,4,MPI_DOUBLE_PRECISION,
     1                       MPI_COMM_WORLD,mpierr)
c
c#ifdef KXCOMM2
	allocate (t_kxcomm2_all(4,0:numtasks-1))
	t_kxcomm2(:,1)=t_kxcomm2(:,1)/ncpusteps
        call  MPI_ALLGATHER (t_kxcomm2,4,MPI_DOUBLE_PRECISION,
     1                       t_kxcomm2_all,4,MPI_DOUBLE_PRECISION,
     1                       MPI_COMM_WORLD,mpierr)
c#else
	allocate (t_kxcomm2t_all(4,0:numtasks-1))
	t_kxcomm2t(:,1)=t_kxcomm2t(:,1)/ncpusteps
        call  MPI_ALLGATHER (t_kxcomm2t,4,MPI_DOUBLE_PRECISION,
     1                       t_kxcomm2t_all,4,MPI_DOUBLE_PRECISION,
     1                       MPI_COMM_WORLD,mpierr)
c#endif
c
	t_xkcomm1(:,1)=t_xkcomm1(:,1)/ncpusteps
        call  MPI_ALLGATHER (t_xkcomm1,4,MPI_DOUBLE_PRECISION,
     1                       t_xkcomm1_all,4,MPI_DOUBLE_PRECISION,
     1                       MPI_COMM_WORLD,mpierr)
c
	t_xkcomm2(:,1)=t_xkcomm2(:,1)/ncpusteps
        call  MPI_ALLGATHER (t_xkcomm2,4,MPI_DOUBLE_PRECISION,
     1                       t_xkcomm2_all,4,MPI_DOUBLE_PRECISION,
     1                       MPI_COMM_WORLD,mpierr)
c
	t_rks(:,1)=t_rks(:,1)*ncpusteps
	t_itrans(:,1)=t_itrans(:,1)*ncpusteps
	t_trans(:,1)=t_trans(:,1)*ncpusteps
	t_kxcomm1(:,1)=t_kxcomm1(:,1)*ncpusteps
	t_kxcomm2(:,1)=t_kxcomm2(:,1)*ncpusteps
	t_kxcomm2t(:,1)=t_kxcomm2t(:,1)*ncpusteps
	t_xkcomm1(:,1)=t_xkcomm1(:,1)*ncpusteps
	t_xkcomm2(:,1)=t_xkcomm2(:,1)*ncpusteps
c
	if (taskid.eq.0) then
c
#ifdef USE_EVEN
#ifdef KXCOMM2
        write (string,"('USE_EVEN, kxcomm2')")
#else
        write (string,"('USE_EVEN, kxcomm2t')")
#endif
#else
#ifdef KXCOMM2
        write (string,"('uneven, kxcomm2')")
#else
        write (string,"('uneven, kxcomm2t')")
#endif
#endif
c
	if (taskid.eq.0) write (6,*) 'write_timings, string=',string

	itask1=2*iproc-1
	itask2=numtasks-2*iproc
	write (6,*) 'write_timings, itask1,itask2=',itask1,itask2

	filepos='rewind'
	write (6,*) 'write_timings, ichkpt=',ichkpt
	if (ichkpt.gt.1) filepos='append'

 	open (77,file='cpu_rk_all',position=filepos)
      write (77,711) nx,ny,nz,nc,iproc, jproc
 711  format ('Pencils code, nx,ny,nz,nc=',3i6,i3,'  iproc, jproc=',i4,
     1        ' x ',i6)
#ifndef OPENMP
	num_thr=0
#endif
#ifdef DOUBLE_PREC
        write (77,"('num_thr=',i3,'  double prec,  rkmethod=',i3)") num_thr,
     1    rkmethod
#else
        write (77,"('num_thr=',i3,'  single prec,  rkmethod=',i3)") num_thr,
     1    rkmethod
#endif
#ifndef OPENMP
	num_thr=1
#endif

	write (77,"(a20)") string
	write (77,"('Detailed breakdown of CPU costs per time step')")
	write (77,"('averaged over', i4,'  non-output steps')") ncpusteps
	write (77,"('taskid ipid jpid           itransform realspace',
     1             ' transform           overall')")
        do itask=0,numtasks-1
	if (itask.le.itask1.or.itask.ge.itask2) then
	sum=0.
	do i=1,5
	sum=sum+t_rks_all(i,itask)
	end do
        write (77,630) itask,ipid_all(itask),jpid_all(itask),
     1             (t_rks_all(i,itask),i=1,6)
        if (mod(itask+1,iproc).eq.0) then
        write (77,"('--------------------------------------------')")
        end if
	end if
	end do
 630    format (i6,2i5,1x,1p,6e10.3)
	close (77)
c
	open (77,file='aftrans_all',position=filepos)
	write (77,"('Detailed breakdown of CPU costs, after sub. transform, per time step')")
	write (77,"('averaged over', i4,'  non-output steps')") ncpusteps
	write (77,"('taskid ipid jpid    step     wavespace   advanc     barrier    total')")
        do itask=0,numtasks-1
	if (itask.le.itask1.or.itask.ge.itask2) then
        write (77,631) itask,ipid_all(itask),jpid_all(itask),
     1             (t_rks_all(i,itask),i=7,10),t_rks_all(5,itask)
        if (mod(itask+1,iproc).eq.0) then
        write (77,"('--------------------------------------------')")
	end if
	end if
	end do
 631    format (i6,2i5,1p,5e11.3)
	close (77)
c
	open (77,file='t_itrans_all',position=filepos)
	write (77,"('Detailed breakdown of CPU costs, sub. itransform_vel, per time step')")
	write (77,"('averaged over', i4,'  non-output steps')") ncpusteps
#ifdef NOSCALAR
	write (77,"('taskid ipid jpid   kxcomm1    fft(z)     kxcomm2    total')")
#else
	if (nc.gt.0) then
	write (77,"('taskid ipid jpid   kxcomm1    fft(z)     kxcomm2    scalars    total')")
	else
	write (77,"('taskid ipid jpid   kxcomm1    fft(z)     kxcomm2    velgrad    total')")
	end if
#endif
        do itask=0,numtasks-1
	if (itask.le.itask1.or.itask.ge.itask2) then
#ifdef NOSCALAR
        write (77,631) itask,ipid_all(itask),jpid_all(itask),
     1             (t_itrans_all(i,itask),i=1,3),
     1             t_itrans_all(1,itask)+t_itrans_all(2,itask)
     1             +t_itrans_all(3,itask)
#else
        write (77,631) itask,ipid_all(itask),jpid_all(itask),
     1             (t_itrans_all(i,itask),i=1,4),
     1             t_itrans_all(1,itask)+t_itrans_all(2,itask)
     1             +t_itrans_all(3,itask)+t_itrans_all(4,itask)
#endif
        if (mod(itask+1,iproc).eq.0) then
        write (77,"('--------------------------------------------')")
        end if
        end if
	end do
	close (77)
c
	open (77,file='t_trans_all',position=filepos)
	write (77,"('Detailed breakdown of CPU costs, sub. transform_vel, per time step')")
	write (77,"('averaged over', i4,'  non-output steps')") ncpusteps
	write (77,"('taskid ipid jpid   xkcomm1    fft(z)     convec     xkcomm2    total')")
        do itask=0,numtasks-1
	if (itask.le.itask1.or.itask.ge.itask2) then
        write (77,631) itask,ipid_all(itask),jpid_all(itask),
     1             (t_trans_all(i,itask),i=1,4),
     1             t_trans_all(1,itask)+t_trans_all(2,itask)
     1             +t_trans_all(3,itask)+t_trans_all(4,itask)
        if (mod(itask+1,iproc).eq.0) then
        write (77,"('--------------------------------------------')")
        end if
        end if
	end do
	close (77)
c
	open (77,file='t_kxcomm1_all',position=filepos)
	write (77,"('Detailed breakdown of CPU costs, sub. kxcomm1_clean per time step')")
	write (77,"('averaged over', i4,'  non-output steps')") ncpusteps
	write (77,"('taskid ipid jpid  alltoall  comp_fft   comp_other   total')")
        do itask=0,numtasks-1
	if (itask.le.itask1.or.itask.ge.itask2) then
        write (77,631) itask,ipid_all(itask),jpid_all(itask),
     1             (t_kxcomm1_all(i,itask),i=1,3),
     1             t_kxcomm1_all(1,itask)+t_kxcomm1_all(2,itask)
     1             +t_kxcomm1_all(3,itask)
        if (mod(itask+1,iproc).eq.0) then
        write (77,"('--------------------------------------------')")
        end if
        end if
	end do
	close (77)
c
#ifdef KXCOMM2
	if (nc.eq.0) then
	open (77,file='t_kxcomm2_all',position=filepos)
	write (77,"('Detailed breakdown of CPU costs, sub. kxcomm2 per time step')")
	write (77,"('averaged over', i4,'  non-output steps')") ncpusteps
	write (77,"('taskid ipid jpid  alltoall   comp_fft  comp_other   total      overall')")
        do itask=0,numtasks-1
	if (itask.le.itask1.or.itask.ge.itask2) then
        write (77,631) itask,ipid_all(itask),jpid_all(itask),
     1             (t_kxcomm2_all(i,itask),i=1,3),
     1             t_kxcomm2_all(1,itask)+t_kxcomm2_all(2,itask)
     1             +t_kxcomm2_all(3,itask),t_kxcomm2_all(4,itask)
        if (mod(itask+1,iproc).eq.0) then
        write (77,"('--------------------------------------------')")
        end if
        end if
	end do
	else
	open (77,file='t_kxcomm2t_all',position=filepos)
	write (77,"('Detailed breakdown of CPU costs, sub. kxcomm2t per time step')")
	write (77,"('averaged over', i4,'  non-output steps')") ncpusteps
	write (77,"('taskid ipid jpid  alltoall   comp_fft  comp_other   total      overall')")
        do itask=0,numtasks-1
	if (itask.le.itask1.or.itask.ge.itask2) then
        write (77,631) itask,ipid_all(itask),jpid_all(itask),
     1             (t_kxcomm2t_all(i,itask),i=1,3),
     1             t_kxcomm2t_all(1,itask)+t_kxcomm2t_all(2,itask)
     1             +t_kxcomm2t_all(3,itask),t_kxcomm2t_all(4,itask)
        if (mod(itask+1,iproc).eq.0) then
        write (77,"('--------------------------------------------')")
        end if
        end if
	end do
	end if
	close (77)
#else
	open (77,file='t_kxcomm2t_all',position=filepos)
	write (77,"('Detailed breakdown of CPU costs, sub. kxcomm2t per time step')")
	write (77,"('averaged over', i4,'  non-output steps')") ncpusteps
	write (77,"('taskid ipid jpid  alltoall  comp_fft   comp_other   total')")
        do itask=0,numtasks-1
	if (itask.le.itask1.or.itask.ge.itask2) then
        write (77,631) itask,ipid_all(itask),jpid_all(itask),
     1             (t_kxcomm2t_all(i,itask),i=1,3),
     1             t_kxcomm2t_all(1,itask)+t_kxcomm2t_all(2,itask)
     1             +t_kxcomm2t_all(3,itask),t_kxcomm2t_all(4,itask)
        if (mod(itask+1,iproc).eq.0) then
        write (77,"('--------------------------------------------')")
        end if
        end if
	end do
	close (77)
#endif
c
	open (77,file='t_xkcomm1_all',position=filepos)
	write (77,"('Detailed breakdown of CPU costs, sub. xkcomm1_clean per time step')")
	write (77,"('averaged over', i4,'  non-output steps')") ncpusteps
	write (77,"('taskid ipid jpid  alltoall  comp_fft   comp_other   total')")
        do itask=0,numtasks-1
	if (itask.le.itask1.or.itask.ge.itask2) then
        write (77,631) itask,ipid_all(itask),jpid_all(itask),
     1             (t_xkcomm1_all(i,itask),i=1,3),
     1             t_xkcomm1_all(1,itask)+t_xkcomm1_all(2,itask)
     1             +t_xkcomm1_all(3,itask)
        if (mod(itask+1,iproc).eq.0) then
        write (77,"('--------------------------------------------')")
        end if
        end if
	end do
	close (77)
c
	open (77,file='t_xkcomm2_all',position=filepos)
	write (77,"('Detailed breakdown of CPU costs, sub. xkcomm2_clean per time step')")
	write (77,"('averaged over', i4,'  non-output steps')") ncpusteps
	write (77,"('taskid ipid jpid  alltoall   comp_fft  comp_other   total')")
        do itask=0,numtasks-1
	if (itask.le.itask1.or.itask.ge.itask2) then
        write (77,631) itask,ipid_all(itask),jpid_all(itask),
     1             (t_xkcomm2_all(i,itask),i=1,3),
     1             t_xkcomm2_all(1,itask)+t_xkcomm2_all(2,itask)
     1             +t_xkcomm2_all(3,itask)
        if (mod(itask+1,iproc).eq.0) then
        write (77,"('--------------------------------------------')")
        end if
	end if
	end do
	close (77)
c
	i2t=2
#ifdef KXCOMM2
	if (nc.eq.0) i2t=1
#endif

	open (77,file='t_comm_all',position=filepos)
	write (77,"('Detailed breakdown of alltoall communication costs, per time step')")
	write (77,"('averaged over', i4,'  non-output steps')") ncpusteps
	if (i2t.eq.1) then
	write (77,"('taskid ipid jpid    kxcomm1   kxcomm2  xkcomm1   xkcomm2    total      %Code')")
	else
	write (77,"('taskid ipid jpid    kxcomm1   kxcomm2t  xkcomm1   xkcomm2    total     %Code')")
	end if

        do itask=0,numtasks-1
	if (itask.le.itask1.or.itask.ge.itask2) then
	if (i2t.eq.1) then
	comm2=t_kxcomm2_all(1,itask)
	else
	comm2=t_kxcomm2t_all(1,itask)
	end if
        sum=t_kxcomm1_all(1,itask)+comm2+
     1             t_xkcomm1_all(1,itask)+t_xkcomm2_all(1,itask)
        write (77,632) itask,ipid_all(itask),jpid_all(itask),
     1             t_kxcomm1_all(1,itask),comm2,
     1             t_xkcomm1_all(1,itask),t_xkcomm2_all(1,itask),
     1		sum,sum/t_rks_all(6,itask)*100.
        if (mod(itask+1,iproc).eq.0) then
        write (77,"('--------------------------------------------')")
        end if
        end if
	end do
	close (77)
 632    format (i6,2i5,2x,1p,5e10.3,0p,f7.1,'%')
c
	open (77,file='t_loctsp_all',position=filepos)
	write (77,"('Detailed breakdown of local-transpose costs, per time step')")
	write (77,"('averaged over', i4,'  non-output steps')") ncpusteps
	if (i2t.eq.1) then
	write (77,"('taskid ipid jpid    kxcomm1   kxcomm2  xkcomm1   xkcomm2    total      %Code')")
	else
	write (77,"('taskid ipid jpid    kxcomm1   kxcomm2t  xkcomm1   xkcomm2    total     %Code')")
	end if
        do itask=0,numtasks-1
	if (itask.le.itask1.or.itask.ge.itask2) then
	if (i2t.eq.1) then
	comm2=t_kxcomm2_all(3,itask)
	else
	comm2=t_kxcomm2t_all(3,itask)
	end if
        sum=t_kxcomm1_all(3,itask)+comm2+
     1             t_xkcomm1_all(3,itask)+t_xkcomm2_all(3,itask)
        write (77,632) itask,ipid_all(itask),jpid_all(itask),
     1             t_kxcomm1_all(3,itask),comm2,
     1             t_xkcomm1_all(3,itask),t_xkcomm2_all(3,itask),
     1		sum,sum/t_rks_all(6,itask)*100.
        if (mod(itask+1,iproc).eq.0) then
        write (77,"('--------------------------------------------')")
        end if
        end if
	end do
	close (77)
c
	end if


	deallocate (t_rks_all)
c       deallocate (ipid_all)
c       deallocate (jpid_all)
	deallocate (t_itrans_all)
	deallocate (t_trans_all)
	deallocate (t_kxcomm1_all)
	deallocate (t_xkcomm1_all)
	deallocate (t_xkcomm2_all)
	deallocate (t_kxcomm2_all)
	deallocate (t_kxcomm2t_all)

	return
	end
