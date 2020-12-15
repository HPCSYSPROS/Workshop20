      subroutine readpos_comm_h5 (pdata, pind, npdim, nump, gid, molgrp)
c
! D. Buaria, Sep 3, 2015
! New restart using parallel hdf5

! Uses IO.F90 modules by M. Clay

#ifdef LAG
c

	use hdf5
	use IO
	use mpilag
!	use mpicom, only : numtasks
	use compart, only : nppad, xyzl, gx, gy, gz, indir_lag

	implicit none
c
	integer nn,numc,npr,indx
	integer i,i1,i2,k,nr,lu,j,mp
	character*12 caux
	character*2 numer 
	character*1 gchar, mchar
	character*6 dsetx, dsety, dsetz
	character*7 gname, fileg
!	character*20 fname
        character(LEN=FILE_NAME_LENGTH) :: fname

c
	integer npdim, nump, gid, molgrp
	integer npdim2, ntasks, iproc2, nrowall
	real(p8) pdata(nppad*npdim/numtasks,3)
	integer pind(nppad*npdim/numtasks) 

        integer, allocatable :: pindtemp(:)
	real(p8), allocatable :: ptemp(:)

        integer, allocatable :: num_row(:),idisp(:)

	integer itask
	integer ip,igm
	integer ndirs, jwrite

! detailed descriptions of these variables can found in details in 
! OutcommWriteHDF5.F90 by M. Clay
! in ..../trunk/homogeneous/ 
	integer :: info, insize, h5num
	integer(kind=HID_T) :: file_id, group_id, plist_id, h5prec
	integer(kind=HSIZE_T) :: offset0(1), offset(1), counter(1)
	integer(kind=HSIZE_T) :: dimt1(1), dim3(1), dimP(1)
	integer(kind=HSIZE_T) :: stride(1), block(1)

c


	if (taskid.eq.0) write (6,*) 'enter readpos_comm_h5', gid

	if(gid.eq.3) npdim = npdim*nppad

	if(gid.eq.1) gchar = '1'
	if(gid.eq.2) gchar = '2'
	if(gid.eq.3) gchar = 'm'

	gname = 'group-'//gchar
	dsetx = 'pp-'//gchar//'-x'
	dsety = 'pp-'//gchar//'-y'
	dsetz = 'pp-'//gchar//'-z'

	fileg = 'outpos'//gchar
	ndirs = 2**(pow2(jproc)/2)
	ndirs = 2*ndirs
	jwrite = mod(jpid,ndirs)

!	call FileName ('lagres', jwrite, fileg, jpid, ic, fname)
	call FILE_NAME (trim(indir_lag), jwrite, fileg, jpid, 0, fname)

	if(ipid.eq.0) write(6,*) 'reading ',trim(fname)

	allocate(num_row(iproc), idisp(iproc))

	if (npdim.gt.0) then


	dim3 = [INT(3, HSIZE_T)]
!	dimP = [INT(numtasks, HSIZE_T)]
	dimP = [INT(iproc, HSIZE_T)]
	offset0 = [INT(0, HSIZE_T)]

! initialize Fortran HDF5 interface
	call H5OPEN_F (ierr)

	h5prec = H5KIND_TO_TYPE (p8, H5_REAL_KIND)

	CALL H5FOPEN_F (trim(fname), H5F_ACC_RDONLY_F, file_id, ierr)

!	write(6,*) 'taskid, fname=', taskid, trim(fname)

	CALL H5GOPEN_F(file_id, gname, group_id, ierr)   
   
	CALL READ_ATTRIBUTE(group_id, 'npdim', npdim2)
	CALL READ_ATTRIBUTE(group_id, 'numtasks', ntasks)
	CALL READ_ATTRIBUTE(group_id, 'nrowall', nrowall)
	CALL READ_ATTRIBUTE(group_id, 'iproc', iproc2)


	if(npdim.ne.npdim2) then
	write(6,*) 'declared npdim and read npdim2 different', npdim, npdim2
	stop
	endif
	if(ntasks.ne.numtasks) then
	write(6,*) 'numtasks and ntasks dont match', numtasks, ntasks
	stop
	endif
	if(iproc2.ne.iproc) then
	write(6,*) 'iproc2 and iproc dont match', iproc2, iproc
	stop
	endif

	if(ipid.eq.0) then
	CALL READ_DATA(group_id, dimP, num_row, 'num_row', offset0)
	endif

        call MPI_BCAST (num_row, iproc, MPI_INTEGER, 0,
     &                       mpi_comm_row, mpierr)

	if(sum(num_row).ne.nrowall) then
	write(6,*) 'sum(num_row) and nrowall dont match', sum(num_row), nrowall
	stop
	endif

	idisp(1) = 0
	do itask=2,iproc
	idisp(itask) = idisp(itask-1) + num_row(itask-1)
	enddo 
	offset = [INT(idisp(ipid+1), HSIZE_T)]

	nump = num_row(ipid+1)

	allocate (ptemp(nump))
	allocate (pindtemp(nump))

	block = [INT(nump, HSIZE_T)]

	CALL READ_DATA(group_id, block, ptemp, dsetx, offset)
	pdata(1:nump,1) = ptemp(:)

	CALL READ_DATA(group_id, block, ptemp, dsety, offset)
	pdata(1:nump,2) = ptemp(:)

	CALL READ_DATA(group_id, block, ptemp, dsetz, offset)
	pdata(1:nump,3) = ptemp(:)

	if(gid.ne.3) then
	CALL READ_DATA(group_id, block, pindtemp, 'particle-index', offset)
	endif

	pind(1:nump) = pindtemp(:)


	deallocate (ptemp, pindtemp)


	call H5GCLOSE_F (group_id, ierr)
	call H5FCLOSE_F (file_id, ierr)

	call H5CLOSE_F (ierr)


	endif

	if (taskid.eq.0) write (6,*) ' exit readpos_h5'

#endif

	return

	end subroutine readpos_comm_h5
