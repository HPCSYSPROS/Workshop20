      subroutine lagout_comm_h5 (pdata, pind, nvar, npdim, nump, gid)
c
! D. Buaria, Sep 3, 2015
! New output routine using parallel hdf5

! Uses IO.F90 modules by M. Clay


c task 0 for fluid particles
c tasks 2 and numtasks-3 for molecules
c
#ifdef LAG
c

	use hdf5
	use IO
	use mpilag
	use compart, only : nppad, xyzl, gx, gy, gz,istep, time

	implicit none
c
	integer nn,numc,npr,indx, pid
	integer i,i1,i2,k,nr,lu,j,mp,iv
	character*12 caux
	character*2 numer 
	character*1 gchar, mchar
	character*6 dsetx, dsety, dsetz
	character*7 gname, fileg
	character*20 grpname, dname
	character(LEN=FILE_NAME_LENGTH) :: fname
c
	integer npdim, nump, nvar,gid, ic, molgrp
	real(p8) pdata(nppad*npdim/numtasks,nvar)
	integer pind(nppad*npdim/numtasks) 

        integer, allocatable :: num_all(:),idisp(:)

	integer itask,ntotal, isubset
	integer ip,igm
	integer ndirs, jwrite

	integer, save :: jcall(3) = (/ 0, 0, 0/) 

! detailed descriptions of these variables can found in details in 
! OutcommWriteHDF5.F90 by M. Clay
! in ..../trunk/homogeneous/ 
	integer :: info, insize, h5num
	integer(kind=HID_T) :: file_id, group_id, plist_id, h5prec
	integer(kind=HSIZE_T) :: dimpp1(1), offset(1), counter(1)
	integer(kind=HSIZE_T) :: dimt1(1), dim3(1), dimP(1)
	integer(kind=HSIZE_T) :: stride(1), block(1)
      integer :: itmp

c
	jcall(gid) = jcall(gid) + 1


	if (taskid.eq.0) write (6,*) 'enter lagout_comm_h5, gid, icall=', gid, jcall(gid)

	if(gid.eq.3) npdim = npdim*nppad

	if(gid.eq.1) gchar = '1'
	if(gid.eq.2) gchar = '2'
	if(gid.eq.3) gchar = 'm'

	gname = 'group-'//gchar
	fileg = 'outpos'//gchar

	ndirs = 2**(pow2(jproc)/2)
	ndirs = 2*ndirs
	jwrite = mod(jpid,ndirs)
	call FILE_NAME ('lagout', jwrite, fileg, jpid, fname)

	if(ipid.eq.0) write(6,*) 'taskid, fname =', taskid, trim(fname)

	allocate(num_all(iproc), idisp(iproc))

	if (npdim.gt.0) then

        call MPI_ALLGATHER (nump, 1, MPI_INTEGER, num_all, 
     &                       1, MPI_INTEGER, mpi_comm_row, mpierr)

	idisp(1) = 0
	do itask=2,iproc
	idisp(itask) = idisp(itask-1) + num_all(itask-1)
	enddo 

	ntotal = sum(num_all)

	offset = [INT(idisp(ipid+1), HSIZE_T)]
	counter = [1_HSIZE_T]
	stride = [1_HSIZE_T]
	block = [INT(nump, HSIZE_T)]
	dimpp1 = [INT(ntotal, HSIZE_T)]
	dim3 = [INT(3, HSIZE_T)]
	dimP = [INT(iproc, HSIZE_T)]
	dimt1 = [INT(nppad*npdim/numtasks, HSIZE_T)] 


! initialize Fortran HDF5 interface
	call H5OPEN_F (ierr)

! set up file access property list with parallel IO access
	info = MPI_INFO_NULL
	call H5PCREATE_F (H5P_FILE_ACCESS_F, plist_id, ierr)
	call H5PSET_FAPL_MPIO_F (plist_id, mpi_comm_row, info, ierr)


	if (jcall(gid).eq.1) then
	call H5FCREATE_F (trim(fname), H5F_ACC_TRUNC_F, file_id, ierr, 
     &               access_prp=plist_id)
	else
	call H5FOPEN_F (trim(fname), H5F_ACC_RDWR_F, file_id, ierr, 
     &               access_prp=plist_id)
	endif

! Create the 'Data' group.
	write(grpname,"(a6,i3.3)") 'icall=',jcall(gid)
	call H5GCREATE_F (file_id, trim(grpname), group_id, ierr)
   
! Determine the HDF5 precision type corresponding to the working precision of
! the DNS code.
      h5prec = H5KIND_TO_TYPE (p8, H5_REAL_KIND)
      pid = jpid*iproc
      if(jcall(gid).eq.1) then
         CALL WRITE_ATTRIBUTE (file_id, 'npdim', npdim)
         CALL WRITE_ATTRIBUTE (file_id, 'numtasks', numtasks)
         CALL WRITE_ATTRIBUTE (file_id, 'jproc', jproc)
      endif
      CALL WRITE_ATTRIBUTE (group_id, 'ntotal', ntotal)
      itmp = istep - 1
      CALL WRITE_ATTRIBUTE (group_id, 'istep', itmp)
      CALL WRITE_ATTRIBUTE (group_id, 'time', time)
!	call WRITE_DATA (dim3, 'xyzl', xyzl, h5prec, group_id, proc_id=pid)
!	call WRITE_DATA (dim3, 'gx-gy-gz', [gx(1), gy(1), gz(1)], 
!     &                     h5prec, group_id, proc_id=pid)
      CALL WRITE_DATA (dimP, 'num_all', num_all, group_id, tid=taskid, wid=pid)
      if(gid.ne.3) then
         CALL WRITE_DATA_PARALLEL (group_id, 'particle-index', dimpp1, dimt1,
     &                             offset, counter, stride, pind, dimsW=block)
      endif
      do iv = 1,nvar
         write(dname,"(a4,i2.2)") 'var=',iv
         CALL WRITE_DATA_PARALLEL (group_id, trim(dname), dimpp1, dimt1, offset,
     &                             counter, stride, pdata(:,iv), dimsW=block)
	   enddo
      CALL H5GCLOSE_F (group_id, ierr)
      CALL H5FCLOSE_F (file_id, ierr)
      CALL H5PCLOSE_F (plist_id, ierr)
      CALL H5CLOSE_F (ierr)


	endif

	if (taskid.eq.0) write (6,*) ' exit lagout_comm_h5'

#endif

	return

	end subroutine lagout_comm_h5


