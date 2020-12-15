      subroutine readpos_h5 (pdata, pind, npdim, nump, gid, molgrp)
c
! D. Buaria, Sep 3, 2015
! New restart using parallel hdf5

! Uses IO.F90 modules by M. Clay

#ifdef LAG
c

	use hdf5
	use IO
	use mpilag
	use mpicom, only : numtasks
	use compart, only : nppad, xyzl, gx, gy, gz

	implicit none
c
	integer nn,numc,npr,indx
	integer i,i1,i2,k,nr,lu,j,mp
	character*12 caux
	character*2 numer 
	character*1 gchar, mchar
	character*6 dsetx, dsety, dsetz
	character*7 gname
	character*20 fname
c
	integer npdim, nump, gid, molgrp
	integer npdim2, ntasks
	real(p8) pdata(nppad*npdim/numtasks,3)
	integer pind(nppad*npdim/numtasks) 

        integer, allocatable :: pindtemp(:)
	real(p8), allocatable :: ptemp(:)

        integer, allocatable :: num_all(:),idisp(:)

	integer itask
	integer ip,igm

! detailed descriptions of these variables can found in details in 
! OutcommWriteHDF5.F90 by M. Clay
! in ..../trunk/homogeneous/ 
	integer :: info, insize, h5num
	integer(kind=HID_T) :: file_id, group_id, plist_id, h5prec
	integer(kind=HSIZE_T) :: offset0(1), offset(1), counter(1)
	integer(kind=HSIZE_T) :: dimt1(1), dim3(1), dimP(1)
	integer(kind=HSIZE_T) :: stride(1), block(1)

c


	if (taskid.eq.0) write (6,*) 'enter readpos_h5', gid

	if(gid.eq.3) npdim = npdim*nppad

	if(gid.eq.1) gchar = '1'
	if(gid.eq.2) gchar = '2'
	if(gid.eq.3) gchar = 'm'

	gname = 'group-'//gchar
	dsetx = 'pp-'//gchar//'-x'
	dsety = 'pp-'//gchar//'-y'
	dsetz = 'pp-'//gchar//'-z'

	fname = 'in'//gchar//'pos.h5'

	if(molgrp.ne.0) then
	write(mchar,"(i1)") molgrp
	fname = 'in'//gchar//mchar//'pos.h5'
	endif

	if(taskid.eq.0) write(6,*) 'reading ',trim(fname)

	allocate(num_all(numtasks), idisp(numtasks))

	if (npdim.gt.0) then


	dim3 = [INT(3, HSIZE_T)]
	dimP = [INT(numtasks, HSIZE_T)]
	offset0 = [INT(0, HSIZE_T)]

! initialize Fortran HDF5 interface
	call H5OPEN_F (ierr)

	h5prec = H5KIND_TO_TYPE (p8, H5_REAL_KIND)

	CALL H5FOPEN_F (fname, H5F_ACC_RDONLY_F, file_id, ierr)

	CALL H5GOPEN_F(file_id, gname, group_id, ierr)   
   
	CALL READ_ATTRIBUTE(group_id, 'npdim', npdim2)


	if(npdim.ne.npdim2) then
	write(6,*) 'declared npdim and read npdim2 different', npdim, npdim2
	stop
	endif


	CALL READ_ATTRIBUTE(group_id, 'numtasks', ntasks)

	if(ntasks.ne.numtasks) then
	write(6,*) 'numtasks and ntasks dont match', numtasks, ntasks
	stop
	endif

	if(taskid.eq.0) then
	CALL READ_DATA(group_id, dimP, num_all, 'num_all', offset0)
	endif

        call MPI_BCAST (num_all, numtasks, MPI_INTEGER, 0,
     &                       MPI_COMM_WORLD, mpierr)

	idisp(1) = 0
	do itask=2,numtasks
	idisp(itask) = idisp(itask-1) + num_all(itask-1)
	enddo 
	offset = [INT(idisp(taskid+1), HSIZE_T)]

	nump = num_all(taskid+1)

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

	end subroutine readpos_h5
