      subroutine wrtpos_comm_h5 (pdata, pind, npdim, nump, gid, ic, molgrp)
c
! D. Buaria, Sep 3, 2015
! New checkpointing using parallel hdf5

! Uses IO.F90 modules by M. Clay


c task 0 for fluid particles
c tasks 2 and numtasks-3 for molecules
c
#ifdef LAG
c

	use hdf5
	use IO
	use mpilag
!	use mpicom, only : numtasks
	use compart, only : nppad, xyzl, gx, gy, gz

	implicit none
c
	integer nn,numc,npr,indx, pid
	integer i,i1,i2,k,nr,lu,j,mp
	character*12 caux
	character*2 numer 
	character*1 gchar, mchar
	character*6 dsetx, dsety, dsetz
	character*7 gname, fileg
	character(LEN=FILE_NAME_LENGTH) :: fname
c
	integer npdim, nump, gid, ic, molgrp
	real(p8) pdata(nppad*npdim/numtasks,3)
	integer pind(nppad*npdim/numtasks) 

        integer, allocatable :: num_row(:),idisp(:)

	integer itask,nrowall
	integer ip,igm
	integer ndirs, jwrite

! detailed descriptions of these variables can found in details in 
! OutcommWriteHDF5.F90 by M. Clay
! in ..../trunk/homogeneous/ 
	integer :: info, insize, h5num
	integer(kind=HID_T) :: file_id, group_id, plist_id, h5prec
	integer(kind=HSIZE_T) :: dimpp1(1), offset(1), counter(1)
	integer(kind=HSIZE_T) :: dimt1(1), dim3(1), dimP(1)
	integer(kind=HSIZE_T) :: stride(1), block(1)

c


	if (taskid.eq.0) write (6,*) 'enter wrtpos_comm_h5', gid

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

!	call FileName ('lagcomm', iwgid, 'outpos', jpid, ic, fname)
	call FILE_NAME ('lagres', jwrite, fileg, jpid, ic, fname)



	allocate(num_row(iproc), idisp(iproc))

	if (npdim.gt.0) then

        call MPI_ALLGATHER (nump, 1, MPI_INTEGER, num_row, 
     &                       1, MPI_INTEGER, mpi_comm_row, mpierr)

	idisp(1) = 0
	do itask=2,iproc
	idisp(itask) = idisp(itask-1) + num_row(itask-1)
	enddo 

	nrowall = sum(num_row)

	offset = [INT(idisp(ipid+1), HSIZE_T)]
	counter = [1_HSIZE_T]
	stride = [1_HSIZE_T]
	block = [INT(nump, HSIZE_T)]
	dimpp1 = [INT(nrowall, HSIZE_T)]
	dim3 = [INT(3, HSIZE_T)]
	dimP = [INT(iproc, HSIZE_T)]
	dimt1 = [INT(nppad*npdim/numtasks, HSIZE_T)] 


! initialize Fortran HDF5 interface
	call H5OPEN_F (ierr)

! set up file access property list with parallel IO access
	info = MPI_INFO_NULL
	call H5PCREATE_F (H5P_FILE_ACCESS_F, plist_id, ierr)
	call H5PSET_FAPL_MPIO_F (plist_id, mpi_comm_row, info, ierr)



	call H5FCREATE_F (trim(fname), H5F_ACC_TRUNC_F, file_id, ierr, 
     &               access_prp=plist_id)


   
! Create the 'Data' group.
	call H5GCREATE_F (file_id, gname, group_id, ierr)
   
! Determine the HDF5 precision type corresponding to the working precision of
! the DNS code.
      h5prec = H5KIND_TO_TYPE (p8, H5_REAL_KIND)
      pid = jpid*iproc
      call WRITE_ATTRIBUTE (group_id, 'npdim', npdim)
      call WRITE_ATTRIBUTE (group_id, 'nrowall', nrowall)
      call WRITE_DATA (dim3, 'xyzl', xyzl, group_id, tid=taskid, wid=pid)
      call WRITE_DATA (dim3, 'gx-gy-gz', [gx(1), gy(1), gz(1)], 
     &                 group_id, tid=taskid, wid=pid)
      call WRITE_ATTRIBUTE (group_id, 'numtasks', numtasks)
      call WRITE_ATTRIBUTE (group_id, 'iproc', iproc)
      call WRITE_DATA (dimP, 'num_row', num_row, group_id, tid=taskid, wid=pid)
      if(gid.ne.3) then
         call WRITE_DATA_PARALLEL (group_id, 'particle-index', dimpp1, dimt1,
     &                             offset, counter, stride, pind, dimsW=block)
      endif
      call WRITE_DATA_PARALLEL (group_id, dsetx, dimpp1, dimt1, offset,
     &                          counter, stride, pdata(:,1), dimsW=block)
      call WRITE_DATA_PARALLEL (group_id, dsety, dimpp1, dimt1, offset,
     &                          counter, stride, pdata(:,2), dimsW=block)
      call WRITE_DATA_PARALLEL (group_id, dsetz, dimpp1, dimt1, offset,
     &                          counter, stride, pdata(:,3), dimsW=block)
      call H5GCLOSE_F (group_id, ierr)
      call H5FCLOSE_F (file_id, ierr)
      call H5PCLOSE_F (plist_id, ierr)
      call H5CLOSE_F (ierr)


	endif

	if (taskid.eq.0) write (6,*) ' exit wrtpos_comm_h5'

#endif

	return

	end subroutine wrtpos_comm_h5


