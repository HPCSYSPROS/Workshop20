!> @file RESTART_HDF5.F90
!> @author Matthew Clay
!> @brief Subroutine to restart a simulation from an HDF5 checkpoint.
!!
!! This subroutine assumes that we are reading data from either 1) an HDF5
!! virtual dataset (VDS) linking all checkpoint files into a single dataset, or
!! 2) a single HDF5 file containing all data for a checkpoint. The dimensions of
!! the dataset are obtained from the HDF5 file and are subsequently used to form
!! the HDF5 offsets, etc. that are used to put the data in the correct spots of
!! the working arrays for the current simulation. If the "bwio" attribute in the
!! HDF5 file is active, the ky=0 data for the v velocity is read, and the rest
!! of the field is computed with the continuity equation.
!!
!> @param[in] path Path to the checkpoint file.
!> @param[in] cnum Number of the checkpoint to read.
!> @param[in,out] buff Array to hold the velocity/scalar fields.
SUBROUTINE RESTART_HDF5(path,cnum,buff)
   ! Required modules.
   USE MPI
   USE HDF5
   USE IO
   USE param,ONLY: B8,nc,nx,ny,nz
   USE com,ONLY: kx,ky,kz,b11,b22,b33,imagi,kinit
   USE mpicom,ONLY: numtasks,taskid,xisz,zjsz,yjsz,xist,zjst
   IMPLICIT NONE
   ! Calling arguments.
   CHARACTER(LEN=*),INTENT(IN) :: path
   INTEGER,INTENT(IN) :: cnum
   COMPLEX(KIND=B8),DIMENSION(ny,zjsz,xisz,3+nc),INTENT(INOUT) :: buff
   ! Local variables.
   ! Identifier for the HDF5 file.
   INTEGER(KIND=HID_T) :: file_id
   ! Grid extents for the data in the HDF5 file.
   INTEGER :: nxi,nyi,nzi
   ! Whether the checkpoint used BWIO for the v velocity.
   INTEGER :: bwio_in
   ! Dimensions for the working arrays (for one variable) on this process.
   INTEGER(KIND=HSIZE_T),DIMENSION(3) :: dims_buff
   ! Dimensions for each read (+/- ky wavenumbers in separate reads).
   INTEGER(KIND=HSIZE_T),DIMENSION(3) :: dims_read_1,dims_read_2
   ! Offset arrays when filling data (put +/- ky wavenumber in correct spot).
   INTEGER(KIND=HSIZE_T),DIMENSION(3) :: oset_fill_1,oset_fill_2
   ! Offset arrays when reading data from the file.
   INTEGER(KIND=HSIZE_T),DIMENSION(3) :: oset_read_1,oset_read_2
   ! Dimensions of data for this process when using BWIO.
   INTEGER(KIND=HSIZE_T),DIMENSION(2) :: dims_buff_bwio
   ! Dimensions to read from the file when using BWIO.
   INTEGER(KIND=HSIZE_T),DIMENSION(2) :: dims_read_bwio
   ! Offset arrays for filling data after reading it from HDF5 with BWIO.
   INTEGER(KIND=HSIZE_T),DIMENSION(2) :: oset_fill_bwio
   ! Offset arrays when reading data from the HDF5 file for BWIO.
   INTEGER(KIND=HSIZE_T),DIMENSION(2) :: oset_read_bwio
   ! Wavenumber range in the checkpoint file.
   INTEGER :: kxi1,kxi2,kzi1,kzi2
   ! Wavenumber range for the current simulation.
   INTEGER :: kxs1,kxs2,kzs1,kzs2
   ! Wavenumber extents in x and z of data being read.
   INTEGER :: kxr1,kxr2,kzr1,kzr2
   ! How much is being read by a given process in x and z.
   INTEGER :: exr,ezr
   ! Logic indicating if we are reading in x/z and at all.
   LOGICAL :: read_x,read_z,read_bool
   ! Offset for z in when reading in data.
   INTEGER :: oset_z
   ! String for the file name.
   CHARACTER(LEN=FILE_NAME_LENGTH) :: fname
   ! Buffers to read data from the HDF5 files (2D array if BWIO used).
   REAL(KIND=B8),DIMENSION(:,:,:),ALLOCATABLE :: rbuf
   REAL(KIND=B8),DIMENSION(:,:),ALLOCATABLE :: rbuf_bwio
   ! Variables used to form a list of what to read and where to put it in buff.
   INTEGER :: nvar,nloc,ind
   CHARACTER(LEN=32),DIMENSION(:),ALLOCATABLE :: vlst
   INTEGER,DIMENSION(:),ALLOCATABLE :: vloc
   ! Looping indices.
   INTEGER :: i,j,k,n,x,y,z
   ! Used to access computational wavenumbers when using BWIO.
   INTEGER :: xp,zp
   ! Offsets when filling real data arrays into complex array buff.
   INTEGER :: x1,x2
   ! String buffers for forming scalar names.
   CHARACTER(LEN=16) :: scnum,scnam
   ! Error handling.
   INTEGER :: ierr

   ! Open the Fortran-HDF5 interface.
   CALL H5OPEN_F(ierr)

   ! File name based on pencil, row-comm., or col-comm. formats.
   IF (hdf5_input_type .EQ. PENCIL_HDF5_IO) THEN
      CALL FILE_NAME(TRIM(ADJUSTL(path)),'OUTPEN',cnum,fname)
   ELSE IF (hdf5_input_type .EQ. ROWCOM_HDF5_IO) THEN
      CALL FILE_NAME(TRIM(ADJUSTL(path)),'OUTROW',cnum,fname)
   END IF

   ! Open the checkpoint VDS file.
   CALL H5FOPEN_F(TRIM(ADJUSTL(fname)),H5F_ACC_RDONLY_F,file_id,ierr)

   ! Determine the size of the checkpoint data being read in.
   CALL READ_ATTRIBUTE(file_id,'nx',nxi)
   CALL READ_ATTRIBUTE(file_id,'ny',nyi)
   CALL READ_ATTRIBUTE(file_id,'nz',nzi)
   !
   ! Calculate integer wavenumber extents for the checkpoint. Used when
   ! determining which portions of the data each process must read.
   kxi1 = 0
   kxi2 = nxi/2 - 1
   !
   kzi1 = -nzi/2 + 1
   kzi2 = nzi/2

   ! Calculate integer wavenumber extents for the current process. Note that
   ! checkpoints contain all data in the ky direction. We will later perform two
   ! separate reads to put the -/+ ky wavenumbers in the right place.
   kxs1 = INT(kx(xist))
   kxs2 = INT(kx(xist+xisz-1))
   !
   kzs1 = INT(kz(zjst))
   kzs2 = INT(kz(zjst+zjsz-1))

   ! Figure out if this process has any data in the checkpoint file. Do this
   ! by comparing wavenumber range on process with wavenumber range in file.
   !
   ! X-direction. Simpler than Z direction because only kx>=0 stored in file.
   read_x = .FALSE.
   IF (kxs2 .LE. kxi2) THEN
      read_x = .TRUE.
      kxr1 = kxs1
      kxr2 = kxs2
   END IF
   !
   ! Z-direction.
   read_z = .FALSE.
   IF (((kzs1 .GE. kzi1) .AND. (kzs1 .LE. kzi2)) .OR. &
       ((kzs2 .GE. kzi1) .AND. (kzs2 .LE. kzi2))) THEN
      !
      ! This process overlaps the checkpoint in the z direction.
      read_z = .TRUE.
      !
      ! Lower index.
      IF (kzs1 .LT. kzi1) THEN
         kzr1 = kzi1
      ELSE
         kzr1 = kzs1
      END IF
      !
      ! Upper index.
      IF (kzs2 .GT. kzi2) THEN
         kzr2 = kzi2
      ELSE
         kzr2 = kzs2
      END IF
      !
      ! Offset for filling in z data. When restarting with the same resolution
      ! grid as the checkpoint, this offset will be zero. When starting on a
      ! higher resolution grid, a single process no longer "owns" the kzi/2
      ! wavenumber and the negative "high frequency" kz wavenumbers (-kzi/2+1,
      ! etc.) in the checkpoint. It is necessary to use an offset to make sure
      ! the negative kzi wavenumbers are put in the right spot of the arrays.
      oset_z = kzr1 - INT(kz(zjst))
   END IF
   !
   ! Must read in x and z to read data at all.
   read_bool = read_x .AND. read_z

   ! Only processes reading data should continue.
   IF (read_bool) THEN
      !
      ! How much data is being read by each process. Make sure to take care
      ! of situation when kzr1 is nzi/2 and kzr2 is negative!
      !
      exr = kxr2 - kxr1 + 1
      ezr = kzr2 - kzr1 + 1
      IF (ezr .LT. 0) ezr = ezr + nzi
      !
      ! Buffer used to read datasets from the HDF5 file. Note that this is
      ! allocated according to the size of the current simulation, NOT the size
      ! of the simulation stored in the HDF5 file. We will use HDF5 offsets,
      ! etc. to put data from low-resolution grids into this buffer, if needed.
      ALLOCATE(rbuf(2*xisz,ny,zjsz))
      rbuf(:,:,:) = 0.0_B8
      !
      ! Convert the wavenumber range being read into HDF5 offsets, etc. Note
      ! that we work with two offset arrays because the ky wavenumber must be
      ! read into different portions of the current processes array when using
      ! higher resolution grids. In the future we might be able to do this with
      ! a single read by setting an appropriate dataspace for the memory buffer
      ! used in the code, but for now this works fine.
      !
      IF (kzr1 .LT. 0) THEN
         oset_read_1 = [INT(2*kxr1,HSIZE_T), &
                        0_HSIZE_T, &
                        INT(kzr1+nzi,HSIZE_T)]
         oset_read_2 = [INT(2*kxr1,HSIZE_T), &
                        INT(nyi/2+1,HSIZE_T), &
                        INT(kzr1+nzi,HSIZE_T)]
      ELSE
         oset_read_1 = [INT(2*kxr1,HSIZE_T), &
                        0_HSIZE_T, &
                        INT(kzr1,HSIZE_T)]
         oset_read_2 = [INT(2*kxr1,HSIZE_T), &
                        INT(nyi/2+1,HSIZE_T), &
                        INT(kzr1,HSIZE_T)]
      END IF
      !
      oset_fill_1 = [0_HSIZE_T,0_HSIZE_T,INT(oset_z,HSIZE_T)]
      oset_fill_2 = [0_HSIZE_T,INT(ny-nyi/2+1,HSIZE_T),INT(oset_z,HSIZE_T)]
      !
      dims_read_1 = [INT(2*exr,HSIZE_T), &
                     INT(nyi/2+1,HSIZE_T), &
                     INT(ezr,HSIZE_T)]
      dims_read_2 = [INT(2*exr,HSIZE_T), &
                     INT(nyi/2-1,HSIZE_T), &
                     INT(ezr,HSIZE_T)]
      !
      dims_buff = [INT(2*xisz,HSIZE_T), &
                   INT(ny,HSIZE_T), &
                   INT(zjsz,HSIZE_T)]

      ! Form the list of variables to read from the HDF5 file.
      !
      ! First determine the number of variables.
      !
      CALL READ_ATTRIBUTE(file_id,'bwio',bwio_in)
      IF (bwio_in .EQ. 0) THEN
         nvar = 3
      ELSE IF (bwio_in .EQ. 1) THEN
         nvar = 2
      ELSE
         WRITE(*,*) 'BWIO VARIABLE IN HDF5 FILE TAKES ERRONEOUS VALUE',bwio_in
      END IF
      !
      DO i = 1,nc
         IF (kinit(3+i) .GT. 0) THEN
            nvar = nvar + 1
         END IF
      END DO
      !
      ! Form the list of variables to read from the HDF5 file and their
      ! corresponding indices in the working arrays for the code.
      !
      ALLOCATE(vlst(nvar))
      ALLOCATE(vloc(nvar))
      vlst(:) = ''
      vloc(:) = -1
      !
      IF (bwio_in .EQ. 0) THEN
         vlst(1) = 'u'
         vlst(2) = 'v'
         vlst(3) = 'w'
         vloc(1) = 1
         vloc(2) = 2
         vloc(3) = 3
         nloc = 4
      ELSE
         vlst(1) = 'u'
         vlst(2) = 'w'
         vloc(1) = 1
         vloc(2) = 3
         nloc = 3
      END IF
      !
      DO i = 1,nc
         IF (kinit(3+i) .GT. 0) THEN
            ! The number of the scalar in the checkpoint is stored in kinit.
            WRITE(scnum,10) kinit(3+i)
            WRITE(scnam,20) TRIM(ADJUSTL(scnum))
            !
            vlst(nloc) = scnam
            vloc(nloc) = 3+i
            !
            nloc = nloc + 1
         END IF
      END DO
      !
      ! Print this information to the user. This implicitly assumes taskid 0
      ! will always read data, which should be the case.
      IF (taskid .EQ. 0) THEN
         WRITE(*,30) 'READING FOLLOWING VARIABLES INTO UNY ARRAY'
         DO i = 1,nvar
            WRITE(*,40) TRIM(ADJUSTL(vlst(i))),vloc(i)
         END DO
         IF (bwio_in .EQ. 1) THEN
            WRITE(*,30) 'V VELOCITY READ WITH BWIO STRATEGY.'
         END IF
      END IF

      ! Read in the variables one at a time.
      DO n = 1,nvar
         !
         ! Zero out working buffer.
         rbuf(:,:,:) = 0.0_B8
         !
         ! Separate reads for +/- ky wavenumbers in case current grid is at a
         ! higher resolution than the checkpoint.
         CALL READ_DATA(file_id,dims_buff,rbuf,vlst(n), &
                        dims_read_1,oset_read_1,oset_fill_1)
         CALL READ_DATA(file_id,dims_buff,rbuf,vlst(n), &
                        dims_read_2,oset_read_2,oset_fill_2)
         !
         ! Unpack data into proper PSDNS memory layout.
         ind = vloc(n)
         DO i = 1, xisz
            x1 = 2*i - 1
            x2 = 2*i
            DO k = 1, zjsz
               DO j = 1, ny
                  buff(j,k,i,ind) = rbuf(x1,j,k) + imagi*rbuf(x2,j,k)
               END DO
            END DO
         END DO
      END DO

      ! Free temporary memory used for 3D dataset reading.
      DEALLOCATE(vlst)
      DEALLOCATE(vloc)
      DEALLOCATE(rbuf)

      ! When using BWIO, form the 2D array offsets, etc. and read the data for
      ! the v velocity at ky=0.
      IF (bwio_in .EQ. 1) THEN
         !
         ! These offset arrays are like those above without the ky components.
         IF (kzr1 .LT. 0) THEN
            oset_read_bwio = [INT(2*kxr1,HSIZE_T),INT(kzr1+nzi,HSIZE_T)]
         ELSE
            oset_read_bwio = [INT(2*kxr1,HSIZE_T),INT(kzr1,HSIZE_T)]
         END IF
         oset_fill_bwio = [0_HSIZE_T,INT(oset_z,HSIZE_T)]
         dims_read_bwio = [INT(2*exr,HSIZE_T),INT(ezr,HSIZE_T)]
         dims_buff_bwio = [INT(2*xisz,HSIZE_T),INT(zjsz,HSIZE_T)]
         !
         ! Buffer used to read BWIO data.
         ALLOCATE(rbuf_bwio(2*xisz,zjsz))
         rbuf_bwio(:,:) = 0.0_B8
         !
         ! Read in the v-plane dataset.
         CALL READ_DATA(file_id,dims_buff_bwio,rbuf_bwio,'v-plane', &
                        dims_read_bwio,oset_read_bwio,oset_fill_bwio)
         !
         ! Unpack the data into proper PSDNS memory layout.
         DO i = 1, xisz
            x1 = 2*i - 1
            x2 = 2*i
            DO k = 1, zjsz
               buff(1,k,i,2) = rbuf_bwio(x1,k) + imagi*rbuf_bwio(x2,k)
            END DO
         END DO
         !
         ! Calculate the rest of the v velocity with continuity. Grid metrics
         ! taken into account for non-cubic grids.
         DO x = 1, xisz
            xp = x + xist - 1
            DO z = 1, zjsz
               zp = z + zjst - 1
               DO y = 2, ny
                  buff(y,z,x,2) = -(b11(2)*kx(xp)*buff(y,z,x,1) + &
                                    b33(2)*kz(zp)*buff(y,z,x,3))/(b22(2)*ky(y))
               END DO
            END DO
         END DO
         !
         ! Free memory used for BWIO reading.
         DEALLOCATE(rbuf_bwio)
      END IF
   END IF

   ! Close the HDF5 file.
   CALL H5FCLOSE_F(file_id,ierr)

   ! Close the Fortran-HDF5 interface.
   CALL H5CLOSE_F(ierr)

   ! IO formats.
   10 FORMAT (I4)
   20 FORMAT ('s',A)
   30 FORMAT (A)
   40 FORMAT (4X,A4,I4)
END SUBROUTINE RESTART_HDF5
