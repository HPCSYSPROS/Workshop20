!> @file AXI_SPTR.F90
!> @author Matthew Clay
!> @brief Subroutine to calculate and write out axisymmetric spectra.
!!
!! Due to the size of the axisymmetric spectra, we will write one file per
!! iostep. The user should also consider skipping multiple iosteps of the
!! spectrum files are too big.
!!
!! Users can check the results by comparing the 1D spectra calculated in this
!! subroutine to those in sptr1d.
!!
!! TODO:
!!
!!    - Add OpenMP threading support.
!!
!> @param[in] uny Array for velocity field in cylinder form.
!> @param[in] axi Direction of axisymmetry (kx=1,ky=2,kz=3).
!> @param[in,optional] dir Subdirectory to store the files in.
!> @param[in,optional] nam Root name for the files.
SUBROUTINE AXI_SPTR(uny,axi,dir,nam)

! Required modules.
USE MPI
USE HDF5
USE IO
USE param,ONLY: B8,mpireal,nc,nx,ny,nz,nxhp,nyhp,nzhp
USE mpicom,ONLY: taskid,numtasks,xist,xisz,zjsz,mystart_x,mystart_z,num_al
USE com,ONLY: istep,time,beta1,beta2,beta3,viscos,kx,ky,kz
USE comp,ONLY: mask

IMPLICIT NONE

! Calling arguments.
COMPLEX(KIND=B8),DIMENSION(ny,zjsz*xisz,3+nc),INTENT(IN) :: uny
INTEGER,INTENT(IN) :: axi
CHARACTER(LEN=*),INTENT(IN) :: dir,nam
! Local variables.
! Switches for the direction of axisymmetry.
INTEGER,PARAMETER :: KX_AXI_SPEC = 1
INTEGER,PARAMETER :: KY_AXI_SPEC = 2
INTEGER,PARAMETER :: KZ_AXI_SPEC = 3
! Arrays for the component axisymmetric spectra and dissipation spectra.
REAL(KIND=B8),DIMENSION(:,:,:),ALLOCATABLE :: uijp,uijr,dijp,dijr
! Arrays for the energy spectrum and energy dissipation spectrum.
REAL(KIND=B8),DIMENSION(:,:),ALLOCATABLE :: uiip,uiir,diip,diir
! Used to calculate and bin the energy and dissipation spectra.
REAL(KIND=B8) :: uij,dij
! Size of the axisymmetric arrays in the direction of axisymmetry.
INTEGER :: n_parl
! Size of the axisymmetric arrays in the perpendicular direction.
INTEGER :: n_perp
! Used to set which bins a mode goes into.
INTEGER :: bin_perp,bin_parl
! Integer wavenumbers for each direction.
INTEGER :: k1i,k2i,k3i
! Physical wavenumbers for each direction.
REAL(KIND=B8) :: k1r,k2r,k3r
! Radial wavenumber magnitude.
REAL(KIND=B8) :: kr
! Wavenumber magnitude.
REAL(KIND=B8) :: kmag
! Overall square of the wavenumber magnitude.
REAL(KIND=B8) :: ksq
! Metric factor used to form the axisymmetric spectrum.
REAL(KIND=B8) :: kr_fac
! Metric factors in the parallel and perpendicular directions.
REAL(KIND=B8) :: beta_parl,beta_perp
! Used to form spectral density.
REAL(KIND=B8) :: beta0
! Factor for conjugate symmetry.
REAL(KIND=B8) :: tfact_x
! Factor used when binning modes. Takes into account con. sym. for ky and kz.
REAL(KIND=B8) :: fact
! Variables used for MPI calls.
INTEGER :: ntotal,ierr
! Arrays for 1D parallel spectra.
REAL(KIND=B8),DIMENSION(:,:),ALLOCATABLE :: uij_parl,dij_parl
! Arrays for 1D perpendicular spectra.
REAL(KIND=B8),DIMENSION(:,:),ALLOCATABLE :: uij_perp,dij_perp
! Used to calculate Reynolds stresses.
REAL(KIND=B8),DIMENSION(3) :: rs_parl,rs_perp
! Used to calculate dissipation of each RS component.
REAL(KIND=B8),DIMENSION(3) :: ds_parl,ds_perp
! File name for the output.
CHARACTER(LEN=256) :: fname
! File and group identifiers for HDF5.
INTEGER(KIND=HID_T) :: file_id,group_id
! Used to pass dimensions of datasets to HDF5 routines.
INTEGER(KIND=HSIZE_T),DIMENSION(1) :: dims1
INTEGER(KIND=HSIZE_T),DIMENSION(2) :: dims2
INTEGER(KIND=HSIZE_T),DIMENSION(3) :: dims3
! Normalization factor used when saving axisymmetric spectra.
REAL(KIND=B8) :: stheta
! Metric factors to recover the true velocities.
REAL(KIND=B8),DIMENSION(3) :: beta_arr
REAL(KIND=B8) :: beta_sqr
! Temporary when writing out attributes that need to be INOUT in IO.F90.
INTEGER :: tmp
! Strings to save date and time to file.
CHARACTER(LEN=8) :: date
CHARACTER(LEN=10) :: clock
! Looping indices for variables, cylinder data structure, etc.
INTEGER :: xp,x,y,z,n,a,i,j

! Inform the user of the main inputs.
IF (taskid .EQ. 0) THEN
   WRITE(*,300) 'CALL AXI_SPTR WITH AXI =',axi
END IF

! Metric factor to form the spectral density.
beta0 = beta1*beta2*beta3

! Metric factors to recover true velocities.
beta_arr(1) = beta1**2
beta_arr(2) = beta2**2
beta_arr(3) = beta3**2

! Based on the direction of axisymmetry, set the allocation sizes, etc.
SELECT CASE (axi)
   CASE (KX_AXI_SPEC)
      !
      ! Check that nz and ny are the same.
      IF ((ny .NE. nz) .OR. (ABS(beta2-beta3)/beta2 .GT. 1.0E-04_B8)) THEN
         IF (taskid .EQ. 0) THEN
            WRITE(*,200) 'KX AXI SPECTRA REQUIRE NY=NZ BETA2=BETA3. EXITING.'
         END IF
         RETURN
      END IF
      !
      n_parl = nxhp
      n_perp = nyhp
      !
      kr_fac = 2.0_B8*beta2/beta0
      !
      beta_parl = beta1
      beta_perp = beta2
   CASE (KY_AXI_SPEC)
      !
      ! Check that nx and nz are the same. NOTE: need to add metric check.
      IF ((nx .NE. nz) .OR. (ABS(beta1-beta3)/beta1 .GT. 1.0E-04_B8)) THEN
         IF (taskid .EQ. 0) THEN
            WRITE(*,200) 'KY AXI SPECTRA REQUIRE NX=NZ BETA1=BETA3. EXITING.'
         END IF
         RETURN
      END IF
      !
      n_parl = nyhp
      n_perp = nxhp
      !
      kr_fac = 2.0_B8*beta1/beta0
      !
      beta_parl = beta2
      beta_perp = beta1
   CASE (KZ_AXI_SPEC)
      !
      ! Check that nx and ny are the same. NOTE: need to add metric check.
      IF ((nx .NE. ny) .OR. (ABS(beta1-beta2)/beta1 .GT. 1.0E-04_B8)) THEN
         IF (taskid .EQ. 0) THEN
            WRITE(*,200) 'KZ AXI SPECTRA REQUIRE NX=NY BETA1=BETA2. EXITING.'
         END IF
         RETURN
      END IF
      !
      n_parl = nzhp
      n_perp = nxhp
      !
      kr_fac = 2.0_B8*beta1/beta0
      !
      beta_parl = beta3
      beta_perp = beta1
   CASE DEFAULT
      !
      ! Exit.
      IF (taskid .EQ. 0) THEN
         WRITE(*,200) 'INVALID OPTION FOR AXI SPECTRUM DIRECTION. EXITING.'
      END IF
      RETURN
END SELECT

! Allocate memory. Three entries for E11,E22,E33 and one for Eii/2.
ALLOCATE(uijp(n_perp,n_parl,3))
ALLOCATE(uijr(n_perp,n_parl,3))
ALLOCATE(dijp(n_perp,n_parl,3))
ALLOCATE(dijr(n_perp,n_parl,3))
ALLOCATE(uiip(n_perp,n_parl))
ALLOCATE(uiir(n_perp,n_parl))
ALLOCATE(diip(n_perp,n_parl))
ALLOCATE(diir(n_perp,n_parl))
uijp(:,:,:) = 0.0_B8
uijr(:,:,:) = 0.0_B8
dijp(:,:,:) = 0.0_B8
dijr(:,:,:) = 0.0_B8
uiip(:,:) = 0.0_B8
uiir(:,:) = 0.0_B8
diip(:,:) = 0.0_B8
diir(:,:) = 0.0_B8

! Loop over the cylinder data structure to bin the spectra.
DO n = 1,3
   !
   ! Metri factor to use true velocities in calculations.
   beta_sqr = beta_arr(n)
   !
   xp = mystart_x
   z = mystart_z
   aloop: DO a = 1,num_al
      x = xp + xist - 1
      !
      ! Conjugate symmetry factor.
      IF (x .EQ. 1) THEN
         tfact_x = 1.0_B8
      ELSE
         tfact_x = 2.0_B8
      END IF
      !
      ! Computational and physical wavenumbers in x and z.
      k1i = INT(kx(x))
      k3i = INT(kz(z))
      k1r = REAL(k1i,B8)*beta1
      k3r = REAL(k3i,B8)*beta3
      yloop: DO y = 1,ny
         IF (.NOT. mask(y,a)) THEN
            CYCLE yloop
         END IF
         !
         ! Computational and physical wavenumbers in y.
         k2i = ky(y)
         k2r = REAL(k2i,B8)*beta2
         !
         ! Square of wavenumber magnitude.
         ksq = k1r**2 + k2r**2 + k3r**2
         !
         ! Determine which bins to store the contributions in.
         SELECT CASE (axi)
            CASE (KX_AXI_SPEC)
               kr = SQRT(REAL(k2i,B8)**2 + REAL(k3i,B8)**2)
               bin_perp = INT(kr + 1.5_B8)
               bin_parl = x
               fact = 1.0_B8
            CASE (KY_AXI_SPEC)
               kr = SQRT(REAL(k1i,B8)**2 + REAL(k3i,B8)**2)
               bin_perp = INT(kr + 1.5_B8)
               IF (y .EQ. 1) THEN
                  bin_parl = y
                  fact = tfact_x
               ELSE IF (y .LT. nyhp) THEN
                  bin_parl = y
                  fact = 1.0_B8
               ELSE
                  bin_parl = ny + 2 - y
                  fact = tfact_x - 1.0_B8
               END IF
            CASE (KZ_AXI_SPEC)
               kr = SQRT(REAL(k1i,B8)**2 + REAL(k2i,B8)**2)
               bin_perp = INT(kr + 1.5_B8)
               IF (z .EQ. 1) THEN
                  bin_parl = z
                  fact = tfact_x
               ELSE IF (z .LT. nzhp) THEN
                  bin_parl = z
                  fact = 1.0_B8
               ELSE
                  bin_parl = nz + 2 - z
                  fact = tfact_x - 1.0_B8
               END IF
         END SELECT
         !
         ! Form the component of the velocity spectrum tensor and diss tensor.
         uij = fact*beta_sqr*REAL(CONJG(uny(y,a,n))*uny(y,a,n),B8)
         dij = ksq*uij
         !
         ! Bin the energy and dissipation spectra.
         uijp(bin_perp,bin_parl,n) = uijp(bin_perp,bin_parl,n) + uij
         dijp(bin_perp,bin_parl,n) = dijp(bin_perp,bin_parl,n) + dij
      END DO yloop
      CALL next_xz(xp,z)
   END DO aloop
END DO

! Take into account metric factors and physical constants.
uijp = kr_fac*uijp
dijp = kr_fac*2.0_B8*viscos*dijp

! Take the trace to get the energy spectrum and dissipation spectrum.
uiip = 0.5_B8*SUM(uijp,DIM=3)
diip = 0.5_B8*SUM(dijp,DIM=3)

! Reduce all contributions to root.
ntotal = n_parl*n_perp*3
CALL MPI_REDUCE(uijp,uijr,ntotal,mpireal,MPI_SUM,0,MPI_COMM_WORLD,ierr)
CALL MPI_REDUCE(dijp,dijr,ntotal,mpireal,MPI_SUM,0,MPI_COMM_WORLD,ierr)
ntotal = n_parl*n_perp
CALL MPI_REDUCE(uiip,uiir,ntotal,mpireal,MPI_SUM,0,MPI_COMM_WORLD,ierr)
CALL MPI_REDUCE(diip,diir,ntotal,mpireal,MPI_SUM,0,MPI_COMM_WORLD,ierr)

! Only the root process continues from here.
IF (taskid .EQ. 0) THEN
   ! Normalize first kr bin since it is smaller than the other bins.
   uijr(1,:,:) = 2.0_B8*uijr(1,:,:)
   dijr(1,:,:) = 2.0_B8*dijr(1,:,:)
   uiir(1,:) = 2.0_B8*uiir(1,:)
   diir(1,:) = 2.0_B8*diir(1,:)

   ! Allocate memory for the 1D spectra.
   ALLOCATE(uij_parl(n_parl,3))
   ALLOCATE(dij_parl(n_parl,3))
   ALLOCATE(uij_perp(n_perp,3))
   ALLOCATE(dij_perp(n_perp,3))

   ! Form the 1D parallel spectra by integrating over the perp. direction.
   uij_parl(:,:) = 0.0_B8
   dij_parl(:,:) = 0.0_B8
   DO n = 1,3
      DO i = 1,n_parl
         uij_parl(i,n) = 0.5_B8*uijr(1,i,n)
         dij_parl(i,n) = 0.5_B8*dijr(1,i,n)
         DO j = 2,n_perp
            uij_parl(i,n) = uij_parl(i,n) + uijr(j,i,n)
            dij_parl(i,n) = dij_parl(i,n) + dijr(j,i,n)
         END DO
      END DO
   END DO
   uij_parl = beta_perp*uij_parl
   dij_parl = beta_perp*dij_parl

   ! Form the 1D perp. spectra by integrating over the parallel direction.
   uij_perp(:,:) = 0.0_B8
   dij_perp(:,:) = 0.0_B8
   DO n = 1,3
      DO i = 1,n_perp
         uij_perp(i,n) = 0.5_B8*uijr(i,1,n)
         dij_perp(i,n) = 0.5_B8*dijr(i,1,n)
         DO j = 2,n_parl
            uij_perp(i,n) = uij_perp(i,n) + uijr(i,j,n)
            dij_perp(i,n) = dij_perp(i,n) + dijr(i,j,n)
         END DO
      END DO
   END DO
   uij_perp = beta_parl*uij_perp
   dij_perp = beta_parl*dij_perp

   ! Calculate the RS tensor components using parallel and perp. spectra.
   rs_parl(:) = 0.0_B8
   rs_perp(:) = 0.0_B8
   DO n = 1,3
      rs_parl(n) = 0.5_B8*uij_parl(1,n)
      DO i = 2,n_parl
         rs_parl(n) = rs_parl(n) + uij_parl(i,n)
      END DO
      rs_parl(n) = beta_parl*rs_parl(n)
      !
      rs_perp(n) = 0.5_B8*uij_perp(1,n)
      DO i = 2,n_perp
         rs_perp(n) = rs_perp(n) + uij_perp(i,n)
      END DO
      rs_perp(n) = beta_perp*rs_perp(n)
   END DO

   ! Calculate the dissipation components using parallel and perp. spectra.
   ds_parl(:) = 0.0_B8
   ds_perp(:) = 0.0_B8
   DO n = 1,3
      ds_parl(n) = 0.5_B8*dij_parl(1,n)
      DO i = 2,n_parl
         ds_parl(n) = ds_parl(n) + dij_parl(i,n)
      END DO
      ds_parl(n) = beta_parl*ds_parl(n)
      !
      ds_perp(n) = 0.5_B8*dij_perp(1,n)
      DO i = 2,n_perp
         ds_perp(n) = ds_perp(n) + dij_perp(i,n)
      END DO
      ds_perp(n) = beta_perp*ds_perp(n)
   END DO

   ! Calculate the normalized spectra. We store these results in the
   ! process-based array used by the root process before the MPI reduction.
   uijp(:,:,:) = 0.0_B8
   dijp(:,:,:) = 0.0_B8
   uiip(:,:) = 0.0_B8
   diip(:,:) = 0.0_B8
   DO n = 1,3
      DO i = 1,n_parl
         DO j = 2,n_perp
            kmag = SQRT((beta_parl*REAL(i-1,B8))**2 + &
                        (beta_perp*REAL(j-1,B8))**2)
            stheta = beta_perp*REAL(j-1,B8)/kmag
            uijp(j,i,n) = uijr(j,i,n)/stheta
            dijp(j,i,n) = dijr(j,i,n)/stheta
         END DO
      END DO
   END DO
   DO i = 1,n_parl
      DO j = 2,n_perp
         kmag = SQRT((beta_parl*REAL(i-1,B8))**2 + &
                     (beta_perp*REAL(j-1,B8))**2)
         stheta = beta_perp*REAL(j-1,B8)/kmag
         uiip(j,i) = uiir(j,i)/stheta
         diip(j,i) = diir(j,i)/stheta
      END DO
   END DO
   uijp(1,:,:) = uijr(1,:,:)
   dijp(1,:,:) = dijr(1,:,:)
   uiip(1,:) = uiir(1,:)
   diip(1,:) = diir(1,:)

   ! Write the data to HDF5.
   !
   ! Determine the file name.
   CALL FILE_NAME(dir,nam,istep,6,fname)
   !
   ! Open the Fortran-HDF5 interface.
   CALL H5OPEN_F(ierr)
   !
   ! Open the HDF5 file.
   CALL H5FCREATE_F(TRIM(fname),H5F_ACC_TRUNC_F,file_id,ierr)
   !
   ! Record when the data was created.
   CALL DATE_AND_TIME(DATE=date,TIME=clock)
   CALL WRITE_ATTRIBUTE(file_id,'DATE',date)
   CALL WRITE_ATTRIBUTE(file_id,'CLOCK',clock)
   !
   ! Write information about the simulation.
   CALL WRITE_ATTRIBUTE(file_id,'NX',nx)
   CALL WRITE_ATTRIBUTE(file_id,'NY',ny)
   CALL WRITE_ATTRIBUTE(file_id,'NZ',nz)
   CALL WRITE_ATTRIBUTE(file_id,'BETA1',beta1)
   CALL WRITE_ATTRIBUTE(file_id,'BETA2',beta2)
   CALL WRITE_ATTRIBUTE(file_id,'BETA3',beta3)
   CALL WRITE_ATTRIBUTE(file_id,'BETA_PERP',beta_perp)
   CALL WRITE_ATTRIBUTE(file_id,'BETA_PARL',beta_parl)
   CALL WRITE_ATTRIBUTE(file_id,'NU',viscos)
   CALL WRITE_ATTRIBUTE(file_id,'BIN_PERP',n_perp)
   CALL WRITE_ATTRIBUTE(file_id,'BIN_PARL',n_parl)
   CALL WRITE_ATTRIBUTE(file_id,'STEP',istep)
   CALL WRITE_ATTRIBUTE(file_id,'TIME',time)
   tmp = axi
   CALL WRITE_ATTRIBUTE(file_id,'AXI_DIR',tmp)
   !
   ! Write out the 1D spectra in the parallel direction.
   CALL H5GCREATE_F(file_id,'PARL_SPECTRA',group_id,ierr)
   dims2 = [INT(n_parl,HSIZE_T),3_HSIZE_T]
   CALL WRITE_DATA(dims2,'UIJ',uij_parl,group_id)
   CALL WRITE_DATA(dims2,'DIJ',uij_parl,group_id)
   dims1 = [3_HSIZE_T]
   CALL WRITE_DATA(dims1,'RS',rs_parl,group_id)
   CALL WRITE_DATA(dims1,'DS',ds_parl,group_id)
   CALL H5GCLOSE_F(group_id,ierr)
   !
   ! Write out the 1D spectra in the perpendicular direction.
   CALL H5GCREATE_F(file_id,'PERP_SPECTRA',group_id,ierr)
   dims2 = [INT(n_perp,HSIZE_T),3_HSIZE_T]
   CALL WRITE_DATA(dims2,'UIJ',uij_perp,group_id)
   CALL WRITE_DATA(dims2,'DIJ',dij_perp,group_id)
   dims1 = [3_HSIZE_T]
   CALL WRITE_DATA(dims1,'RS',rs_perp,group_id)
   CALL WRITE_DATA(dims1,'DS',ds_perp,group_id)
   CALL H5GCLOSE_F(group_id,ierr)
   !
   ! Write out the axisymmetric spectra.
   CALL H5GCREATE_F(file_id,'AXI_SPECTRA',group_id,ierr)
   dims2 = [INT(n_perp,HSIZE_T),INT(n_parl,HSIZE_T)]
   CALL WRITE_DATA(dims2,'EII',uiir,group_id)
   CALL WRITE_DATA(dims2,'EII_NORM',uiip,group_id)
   CALL WRITE_DATA(dims2,'DII',diir,group_id)
   CALL WRITE_DATA(dims2,'DII_NORM',diip,group_id)
   dims3 = [INT(n_perp,HSIZE_T),INT(n_parl,HSIZE_T),3_HSIZE_T]
   CALL WRITE_DATA(dims3,'UIJ',uijr,group_id)
   CALL WRITE_DATA(dims3,'UIJ_NORM',uijp,group_id)
   CALL WRITE_DATA(dims3,'DIJ',dijr,group_id)
   CALL WRITE_DATA(dims3,'DIJ_NORM',dijp,group_id)
   CALL H5GCLOSE_F(group_id,ierr)
   !
   ! Close the HDF5 file.
   CALL H5FCLOSE_F(file_id,ierr)
   !
   ! Close the Fortran-HDF5 interface.
   CALL H5CLOSE_F(ierr)

   ! Free memory for the 1D spectra.
   DEALLOCATE(uij_parl)
   DEALLOCATE(dij_parl)
   DEALLOCATE(uij_perp)
   DEALLOCATE(dij_perp)
END IF

! Free memory.
DEALLOCATE(uijp)
DEALLOCATE(uijr)
DEALLOCATE(dijp)
DEALLOCATE(dijr)
DEALLOCATE(uiip)
DEALLOCATE(uiir)
DEALLOCATE(diip)
DEALLOCATE(diir)

! IO formats.
200 FORMAT (A)
300 FORMAT (A,1X,I2)

END SUBROUTINE AXI_SPTR
