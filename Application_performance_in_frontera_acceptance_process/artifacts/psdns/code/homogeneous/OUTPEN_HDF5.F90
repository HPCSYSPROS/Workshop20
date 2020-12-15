!> @file OUTPEN_HDF5.F90
!> @author Matthew Clay
!> @brief Write out pencil data in HDF5 format.
!!
!! This routine mimics the behavior of outpen_write, but with HDF5 output for
!! the data instead of Fortran binary output. Because HDF5 does not have a
!! native complex type, the inner-most stride of the complex data is doubled
!! while casting the array to real form.
!!
!> @param[in] buf Array with complex Fourier velocity field in real pencil form.
SUBROUTINE OUTPEN_HDF5(buf)
   ! Required modules.
   USE MPI
   USE HDF5
   USE IO
   USE param,ONLY: B8,nx,ny,nz,nc
   USE mpicom,ONLY: taskid,numtasks,xist,xien,xisz,zjst,zjen,zjsz, &
                    ipid,jpid,iproc,jproc
   USE com,ONLY: isave,dtsave,fninit,pr,viscos,beta1,beta2,beta3
   IMPLICIT NONE
   ! Calling arguments.
   REAL(KIND=B8),DIMENSION(2*xisz,ny,zjsz,3+nc),INTENT(IN) :: buf
   ! Local variables.
   ! HDF5 identifier for the file.
   INTEGER(KIND=HID_T) :: file_id
   ! Array to hold the dimensions of the data being written.
   INTEGER(KIND=HSIZE_T),DIMENSION(3) :: dims3
   ! Output file name.
   CHARACTER(LEN=FILE_NAME_LENGTH) :: fname
   ! String buffers for forming scalar names.
   CHARACTER(LEN=16) :: scnum,scname
   ! Character arrays to hold the date and time stamps.
   CHARACTER(LEN=8) :: date
   CHARACTER(LEN=10) :: time
   ! Parameters used to determine which directory this process writes to.
   INTEGER :: power,ndirs,dirnm
   ! Array for the dimensions of attribute arrays (Schmidt number).
   INTEGER(KIND=HSIZE_T),DIMENSION(1) :: dims_att
   ! Buffers used to gether process layout information for the VDS.
   INTEGER,DIMENSION(:),ALLOCATABLE :: sbuf
   INTEGER,DIMENSION(:,:),ALLOCATABLE :: rbuf
   ! Datasets, property lists, and dataspaces required for making the VDS.
   INTEGER(KIND=HID_T) :: dset_id,vspace_id,srcspc_id,plist_id
   ! Arrays used for sizing/offsets when forming the VDS.
   INTEGER(KIND=HSIZE_T),DIMENSION(3) :: dimT_vds,dimT_prc,oset_vds
   INTEGER(KIND=HSIZE_T),DIMENSION(2) :: dimT2D_vds,dimT2D_prc,oset2D_vds
   ! Looping variables and lists for linking datasets into the VDS.
   INTEGER :: nlist,n
   CHARACTER(LEN=16),DIMENSION(:),ALLOCATABLE :: vlist
   ! When using BWIO, the slice direction and index.
   INTEGER :: srank,sind
   ! Working precision of the code as an HDF5 datatype. Only used for VDS.
   INTEGER(KIND=HID_T) :: h5prec
   ! File IO unit for outdata.grid file.
   INTEGER :: iounit
   ! Looping index for scalars.
   INTEGER :: s
   ! Other looping indices, etc.
   INTEGER :: i
   ! Error handling.
   INTEGER :: ierr
   !
   ! Initialize the Fortran HDF5 interface.
   CALL H5OPEN_F(ierr)
   !
   ! Determine the HDF5 type corresponding to the working precision.
   h5prec = H5KIND_TO_TYPE(B8, H5_REAL_KIND)
   !
   ! Determine number of subdirectories and subdirectory for this process.
   power = POW2(numtasks)
   ndirs = 2**(power/2)
   dirnm = MOD(taskid,ndirs)
   !
   ! Open the file for writing.
   chkpt_num = chkpt_num + 1
   CALL FILE_NAME('outpen', dirnm, 'OUTPEN', taskid, chkpt_num, fname)
   CALL H5FCREATE_F(fname, H5F_ACC_TRUNC_F, file_id, ierr)
   !
   ! Write attributes summarizing the data.
   CALL WRITE_ATTRIBUTE(file_id, 'nc', nc)
   CALL WRITE_ATTRIBUTE(file_id, 'nx', nx)
   CALL WRITE_ATTRIBUTE(file_id, 'ny', ny)
   CALL WRITE_ATTRIBUTE(file_id, 'nz', nz)
   CALL WRITE_ATTRIBUTE(file_id, 'b1', beta1)
   CALL WRITE_ATTRIBUTE(file_id, 'b2', beta2)
   CALL WRITE_ATTRIBUTE(file_id, 'b3', beta3)
   CALL WRITE_ATTRIBUTE(file_id, 'xist', xist)
   CALL WRITE_ATTRIBUTE(file_id, 'xisz', xisz)
   CALL WRITE_ATTRIBUTE(file_id, 'zjst', zjst)
   CALL WRITE_ATTRIBUTE(file_id, 'zjsz', zjsz)
   CALL WRITE_ATTRIBUTE(file_id, 'ipid', ipid)
   CALL WRITE_ATTRIBUTE(file_id, 'jpid', jpid)
   CALL WRITE_ATTRIBUTE(file_id, 'iproc', iproc)
   CALL WRITE_ATTRIBUTE(file_id, 'jproc', jproc)
   CALL WRITE_ATTRIBUTE(file_id, 'bwio', hdf5_output_bwio)
   CALL WRITE_ATTRIBUTE(file_id, 'viscos', viscos)
   IF (nc .GT. 0) THEN
      dims_att(1) = INT(nc,HSIZE_T)
      CALL WRITE_ATTRIBUTE(file_id, 'schmidt', dims_att, pr)
   END IF
   !
   ! Record when the data was written for future reference.
   CALL DATE_AND_TIME(DATE=date, TIME=time)
   CALL WRITE_ATTRIBUTE(file_id, 'date', date)
   CALL WRITE_ATTRIBUTE(file_id, 'time', time)
   !
   ! Write out the velocity field.
   dims3 = [INT(2*xisz,HSIZE_T), INT(ny,HSIZE_T), INT(zjsz,HSIZE_T)]
   CALL WRITE_DATA(dims3, 'u', buf(:,:,:,1), file_id)
   IF (hdf5_output_bwio .EQ. BWIO_HDF5_ON) THEN
      srank = 2
      sind = 1
      CALL WRITE_SLICE(dims3, 'v-plane', buf(:,:,:,2), file_id, srank, sind)
   ELSE
      CALL WRITE_DATA(dims3, 'v', buf(:,:,:,2), file_id)
   END IF
   CALL WRITE_DATA(dims3, 'w', buf(:,:,:,3), file_id)
   !
   ! Write out any scalars.
   DO s = 1, nc
      WRITE(scnum,10) s
      WRITE(scname,20) TRIM(ADJUSTL(scnum))
      CALL WRITE_DATA(dims3, TRIM(scname), buf(:,:,:,3+s), file_id, &
                     real_att_name='schmidt', real_att_val=pr(s))
   END DO
   !
   ! Close the file.
   CALL H5FCLOSE_F(file_id, ierr)
   !
   ! Gather data required to form the VDS.
   ALLOCATE(sbuf(5))
   ALLOCATE(rbuf(5,numtasks))
   sbuf(1) = taskid
   sbuf(2) = xist
   sbuf(3) = xien
   sbuf(4) = zjst
   sbuf(5) = zjen
   CALL MPI_GATHER(sbuf,5,MPI_INTEGER,rbuf,5,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
   !
   ! Write out outdata.grid in case the VDS file gets lost or corrupted.
   IF (taskid .EQ. 0) THEN
      OPEN(NEWUNIT=iounit,FILE='outdata.grid',FORM='FORMATTED',STATUS='REPLACE')
      WRITE(iounit,*) numtasks,iproc,jproc
      DO i = 1,numtasks
         WRITE (iounit,90) rbuf(:,i)
      END DO
      CLOSE(UNIT=iounit)
   END IF
   !
   ! Root process makes a virtual dataset file to link all files together.
#ifdef NEW_HDF5
   IF (taskid .EQ. 0) THEN
      !
      ! Create the file for the virtual dataset.
      CALL FILE_NAME('outpen', 'OUTPEN', chkpt_num, fname)
      CALL H5FCREATE_F(TRIM(fname), H5F_ACC_TRUNC_F, file_id, ierr)
      !
      ! Write attributes summarizing the data.
      CALL WRITE_ATTRIBUTE(file_id, 'nc', nc)
      CALL WRITE_ATTRIBUTE(file_id, 'nx', nx)
      CALL WRITE_ATTRIBUTE(file_id, 'ny', ny)
      CALL WRITE_ATTRIBUTE(file_id, 'nz', nz)
      CALL WRITE_ATTRIBUTE(file_id, 'b1', beta1)
      CALL WRITE_ATTRIBUTE(file_id, 'b2', beta2)
      CALL WRITE_ATTRIBUTE(file_id, 'b3', beta3)
      CALL WRITE_ATTRIBUTE(file_id, 'iproc', iproc)
      CALL WRITE_ATTRIBUTE(file_id, 'jproc', jproc)
      CALL WRITE_ATTRIBUTE(file_id, 'bwio', hdf5_output_bwio)
      CALL WRITE_ATTRIBUTE(file_id, 'viscos', viscos)
      IF (nc .GT. 0) THEN
         dims_att(1) = INT(nc,HSIZE_T)
         CALL WRITE_ATTRIBUTE(file_id, 'schmidt', dims_att, pr)
      END IF
      !
      ! Record when data was written for future reference.
      CALL DATE_AND_TIME(DATE=date, TIME=time)
      CALL WRITE_ATTRIBUTE(file_id, 'date', date)
      CALL WRITE_ATTRIBUTE(file_id, 'time', time)
      !
      ! Create the dataspace for the virtual dataset.
      dimT_vds = [INT(nx,HSIZE_T), INT(ny,HSIZE_T), INT(nz,HSIZE_T)]
      CALL H5SCREATE_SIMPLE_F(3, dimT_vds, vspace_id, ierr)
      !
      ! The dataspace for all files (source dataspaces) is the same.
      dimT_prc = [INT(2*xisz,HSIZE_T), INT(ny,HSIZE_T), INT(zjsz,HSIZE_T)]
      CALL H5SCREATE_SIMPLE_F(3, dimT_prc, srcspc_id, ierr)
      !
      ! Make a list of all variables to link in the VDS.
      IF (hdf5_output_bwio .EQ. BWIO_HDF5_ON) THEN
         nlist = 2 + nc
         ALLOCATE(vlist(nlist))
         vlist(1) = 'u'
         vlist(2) = 'w'
         DO s = 1, nc
            WRITE(scnum,10) s
            WRITE(scname,20) TRIM(ADJUSTL(scnum))
            vlist(2+s) = TRIM(ADJUSTL(scname))
         END DO
      ELSE
         nlist = 3 + nc
         ALLOCATE(vlist(nlist))
         vlist(1) = 'u'
         vlist(2) = 'v'
         vlist(3) = 'w'
         DO s = 1, nc
            WRITE(scnum,10) s
            WRITE(scname,20) TRIM(ADJUSTL(scnum))
            vlist(3+s) = TRIM(ADJUSTL(scname))
         END DO
      END IF
      !
      ! Loop over each variable and create a VDS for it.
      DO n = 1,nlist
         !
         ! VDS creation property list for each variable.
         CALL H5PCREATE_F(H5P_DATASET_CREATE_F, plist_id, ierr)
         CALL H5PSET_FILL_VALUE_F(plist_id, h5prec, -999.999_B8, ierr)
         !
         ! Loop over each process-based file to link to the VDS.
         DO i = 1,numtasks
            !
            ! Figure out the name of the file for this dataset.
            dirnm = MOD(rbuf(1,i),ndirs)
            CALL FILE_NAME('inpen',dirnm,'OUTPEN',rbuf(1,i),chkpt_num,fname)
            !
            ! Form the offset array for this dataset in the VDS.
            oset_vds = [INT(2*(rbuf(2,i)-1),HSIZE_T), &
                        0_HSIZE_T, &
                        INT(rbuf(4,i)-1,HSIZE_T)]
            !
            ! Select the hyperslab from the VDS for this checkpoint.
            CALL H5SSELECT_HYPERSLAB_F(vspace_id,H5S_SELECT_SET_F, &
                                       oset_vds,dimT_prc,ierr)
            !
            ! Map the hyperslab to the dataset in the HDF5 file.
            CALL H5PSET_VIRTUAL_F(plist_id,vspace_id,fname, &
                                  TRIM(vlist(n)),srcspc_id,ierr)
         END DO
         !
         ! Now that the mapping is complete, create the VDS.
         CALL H5DCREATE_F(file_id,TRIM(vlist(n)),h5prec,vspace_id, &
                          dset_id,ierr,plist_id)
         !
         ! Close the virtual dataset.
         CALL H5DCLOSE_F(dset_id,ierr)
         !
         ! Close the property list.
         CALL H5PCLOSE_F(plist_id, ierr)
      END DO
      !
      ! Close the dataspaces before using them for BWIO.
      CALL H5SCLOSE_F(srcspc_id, ierr)
      CALL H5SCLOSE_F(vspace_id, ierr)
      !
      IF (hdf5_output_bwio .EQ. BWIO_HDF5_ON) THEN
         !
         ! Create the dataspace for the virtual dataset.
         dimT2D_vds = [INT(nx,HSIZE_T), INT(nz,HSIZE_T)]
         CALL H5SCREATE_SIMPLE_F(2, dimT2D_vds, vspace_id, ierr)
         !
         ! VDS creation property list.
         CALL H5PCREATE_F(H5P_DATASET_CREATE_F, plist_id, ierr)
         CALL H5PSET_FILL_VALUE_F(plist_id, h5prec, -999.999_B8, ierr)
         !
         ! The dataspace for all files (source dataspaces) is the same.
         dimT2D_prc = [INT(2*xisz,HSIZE_T), INT(zjsz,HSIZE_T)]
         CALL H5SCREATE_SIMPLE_F(2, dimT2D_prc, srcspc_id, ierr)
         !
         ! Loop over each process-based file to link to the VDS.
         DO i = 1,numtasks
            !
            ! Figure out the name of the file for this dataset.
            dirnm = MOD(rbuf(1,i),ndirs)
            CALL FILE_NAME('inpen',dirnm,'OUTPEN',rbuf(1,i),chkpt_num,fname)
            !
            ! Form the offset array for this dataset in the VDS.
            oset2D_vds = [INT(2*(rbuf(2,i)-1),HSIZE_T), &
                          INT(rbuf(4,i)-1,HSIZE_T)]
            !
            ! Select the hyperslab from the VDS for this checkpoint.
            CALL H5SSELECT_HYPERSLAB_F(vspace_id,H5S_SELECT_SET_F, &
                                       oset2D_vds,dimT2D_prc,ierr)
            !
            ! Map the hyperslab to the dataset in the HDF5 file.
            CALL H5PSET_VIRTUAL_F(plist_id,vspace_id,fname, &
                                  'v-plane',srcspc_id,ierr)
         END DO
         !
         ! Now that the mapping is complete, create the VDS.
         CALL H5DCREATE_F(file_id,'v-plane',h5prec,vspace_id, &
                          dset_id,ierr,plist_id)
         !
         ! Close the property list and dataspaces.
         CALL H5PCLOSE_F(plist_id, ierr)
         CALL H5SCLOSE_F(srcspc_id, ierr)
         CALL H5SCLOSE_F(vspace_id, ierr)
         !
         ! Close the virtual dataset.
         CALL H5DCLOSE_F(dset_id,ierr)
      END IF
      !
      ! Close the HDF5 file for the VDS.
      CALL H5FCLOSE_F(file_id, ierr)
      !
      ! Free memory on root process.
      DEALLOCATE(vlist)
   END IF
#endif
   !
   ! Close the Fortran-HDF5 interface.
   CALL H5CLOSE_F(ierr)
   !
   ! Free memory.
   DEALLOCATE(sbuf)
   DEALLOCATE(rbuf)
   !
   ! IO formats.
   10 FORMAT (I4)
   20 FORMAT ('s', A)
   90 FORMAT (5I6)
END SUBROUTINE OUTPEN_HDF5
