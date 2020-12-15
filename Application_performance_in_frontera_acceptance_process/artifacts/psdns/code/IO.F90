!> @file IO.F90
!> @author Matthew Clay
!> @brief General IO routies.
!!
!! This module provides procedures to form file names and write/read data
!! to/from HDF5 files. Since HDF5 does not natively support a complex datatype,
!! all arrays are written/read as reals. We also provide single and double
!! subroutines in the case that multiple precisions are used in the code.
!!
!! NOTES:
!!
!!    - 2016-05: Found bug in HDF5 1.10.0 H5AWRITE_F interface routines. They
!!      currently require INTENT(INOUT) for all variables instead of INTENT(IN).
!!      Recheck after first patch to 1.10.0 to see if it is fixed.
MODULE IO

   ! Required modules.
   USE HDF5

   IMPLICIT NONE

   !> Length of strings containing file names.
   INTEGER,PARAMETER,PUBLIC :: FILE_NAME_LENGTH = 256
   !> Parameter to define pencil-based HDF5 IO.
   INTEGER,PARAMETER,PUBLIC :: PENCIL_HDF5_IO = 1
   !> Parameter to define row-communicator based HDF5 IO.
   INTEGER,PARAMETER,PUBLIC :: ROWCOM_HDF5_IO = 2
   !> Parameter to define non-BWIO data format (write all velocity fields).
   INTEGER,PARAMETER,PUBLIC :: BWIO_HDF5_OFF = 0
   !> Parameter to defined BWIO format, which writes a reduced set of v.
   INTEGER,PARAMETER,PUBLIC :: BWIO_HDF5_ON = 1

   ! Parameters controlling HDF5 input.
   !
   !> Whether or not to use HDF5 input.
   LOGICAL,PUBLIC :: hdf5_input_flag = .FALSE.
   !> Whether the input is pencil or row-communicator based.
   INTEGER,PUBLIC :: hdf5_input_type = -1
   !> Integer for the initial checkpoint when restarting a simulation.
   INTEGER,PUBLIC :: hdf5_init = -1

   ! Parameters controlling HDF5 output.
   !
   !> Whether or not to use HDF5 output.
   LOGICAL,PUBLIC :: hdf5_output_flag = .FALSE.
   !> Whether the output should be in pencil form or row-comm. form.
   INTEGER,PUBLIC :: hdf5_output_type = -1
   !> Wheter or not to use BWIO format for the velocity field.
   INTEGER,PUBLIC :: hdf5_output_bwio = BWIO_HDF5_OFF
   !> Integer to keep track of the latest file number written out.
   INTEGER,PUBLIC :: chkpt_num = 0

   !> Interface for forming output file names.
   INTERFACE FILE_NAME
      MODULE PROCEDURE FILE_NAME_1
      MODULE PROCEDURE FILE_NAME_2
      MODULE PROCEDURE FILE_NAME_3
      MODULE PROCEDURE FILE_NAME_4
      MODULE PROCEDURE FILE_NAME_5
      MODULE PROCEDURE FILE_NAME_6
   END INTERFACE FILE_NAME

   !> Interface for writing attributes to HDF5 files.
   INTERFACE WRITE_ATTRIBUTE
      MODULE PROCEDURE WRITE_ATTRIBUTE_INT32
      MODULE PROCEDURE WRITE_ATTRIBUTE_INT32_RANK1
#ifdef NEW_HDF5
      MODULE PROCEDURE WRITE_ATTRIBUTE_INT64
      MODULE PROCEDURE WRITE_ATTRIBUTE_INT64_RANK1
#endif
      MODULE PROCEDURE WRITE_ATTRIBUTE_REAL32
      MODULE PROCEDURE WRITE_ATTRIBUTE_REAL32_RANK1
      MODULE PROCEDURE WRITE_ATTRIBUTE_REAL64
      MODULE PROCEDURE WRITE_ATTRIBUTE_REAL64_RANK1
      MODULE PROCEDURE WRITE_ATTRIBUTE_CHAR
   END INTERFACE WRITE_ATTRIBUTE

   !> Interface to read attributes from HDF5 files.
   INTERFACE READ_ATTRIBUTE
      MODULE PROCEDURE READ_ATTRIBUTE_INT32
      MODULE PROCEDURE READ_ATTRIBUTE_INT32_RANK1
#ifdef NEW_HDF5
      MODULE PROCEDURE READ_ATTRIBUTE_INT64
      MODULE PROCEDURE READ_ATTRIBUTE_INT64_RANK1
#endif
      MODULE PROCEDURE READ_ATTRIBUTE_REAL32
      MODULE PROCEDURE READ_ATTRIBUTE_REAL32_RANK1
      MODULE PROCEDURE READ_ATTRIBUTE_REAL64
      MODULE PROCEDURE READ_ATTRIBUTE_REAL64_RANK1
   END INTERFACE READ_ATTRIBUTE

   !> Interface for writing datasets to HDF5 files.
   INTERFACE WRITE_DATA
      MODULE PROCEDURE WRITE_DATA_INT32_RANK1
      MODULE PROCEDURE WRITE_DATA_INT32_RANK2
#ifdef NEW_HDF5
      MODULE PROCEDURE WRITE_DATA_INT64_RANK1
      MODULE PROCEDURE WRITE_DATA_INT64_RANK2
#endif
      MODULE PROCEDURE WRITE_DATA_REAL32_RANK1
      MODULE PROCEDURE WRITE_DATA_REAL32_RANK2
      MODULE PROCEDURE WRITE_DATA_REAL32_RANK3
      MODULE PROCEDURE WRITE_DATA_SUBSET_REAL32_RANK3
      MODULE PROCEDURE WRITE_DATA_REAL64_RANK1
      MODULE PROCEDURE WRITE_DATA_REAL64_RANK2
      MODULE PROCEDURE WRITE_DATA_REAL64_RANK3
      MODULE PROCEDURE WRITE_DATA_SUBSET_REAL64_RANK3
   END INTERFACE WRITE_DATA

   !> Interface for writing datasets to HDF5 files in parallel.
   INTERFACE WRITE_DATA_PARALLEL
      MODULE PROCEDURE WRITE_DATA_PARALLEL_INT32_RANK1
#ifdef NEW_HDF5
      MODULE PROCEDURE WRITE_DATA_PARALLEL_INT64_RANK1
#endif
      MODULE PROCEDURE WRITE_DATA_PARALLEL_REAL32_RANK1
      MODULE PROCEDURE WRITE_DATA_PARALLEL_REAL32_RANK3
      MODULE PROCEDURE WRITE_DATA_PARALLEL_REAL64_RANK1
      MODULE PROCEDURE WRITE_DATA_PARALLEL_REAL64_RANK3
   END INTERFACE WRITE_DATA_PARALLEL

   !> Interface to write a slice of data to an HDF5 file.
   INTERFACE WRITE_SLICE
      MODULE PROCEDURE WRITE_SLICE_REAL32_RANK3
      MODULE PROCEDURE WRITE_SLICE_REAL64_RANK3
   END INTERFACE WRITE_SLICE

   !> Interface to write a slice of data to an HDF5 file in parallel.
   INTERFACE WRITE_SLICE_PARALLEL
      MODULE PROCEDURE WRITE_SLICE_PARALLEL_REAL32_RANK3
      MODULE PROCEDURE WRITE_SLICE_PARALLEL_REAL64_RANK3
   END INTERFACE WRITE_SLICE_PARALLEL

   !> Interface for reading datasets from HDF5 files.
   INTERFACE READ_DATA
      MODULE PROCEDURE READ_DATA_INT32_RANK1
#ifdef NEW_HDF5
      MODULE PROCEDURE READ_DATA_INT64_RANK1
#endif
      MODULE PROCEDURE READ_DATA_REAL32_RANK1
      MODULE PROCEDURE READ_DATA_REAL32_RANK2
      MODULE PROCEDURE READ_DATA_PARTIAL_REAL32_RANK2
      MODULE PROCEDURE READ_DATA_REAL32_RANK3
      MODULE PROCEDURE READ_DATA_PARTIAL_REAL32_RANK3
      MODULE PROCEDURE READ_DATA_REAL64_RANK1
      MODULE PROCEDURE READ_DATA_REAL64_RANK2
      MODULE PROCEDURE READ_DATA_PARTIAL_REAL64_RANK2
      MODULE PROCEDURE READ_DATA_REAL64_RANK3
      MODULE PROCEDURE READ_DATA_PARTIAL_REAL64_RANK3
   END INTERFACE READ_DATA

   ! Other module procedures.
   PUBLIC :: GET_REAL_PREC, POW2

CONTAINS

   !> Procedure to form filename with root and subdirectories.
   !!
   !! The file name takes the following form:
   !!
   !!    rootdir/subdir/froot_num1_num2.h5,
   !!
   !! where num1 is written with 6 digits, and num2 is written with 2.
   !! Currently, num1 is intended for the process id and num2 for the checkpoint
   !! number.
   !!
   !> @param[in] rootdir Root directory for the file.
   !> @param[in] subdir Subdirectory for the file.
   !> @param[in] froot Root name for the output file.
   !> @param[in] num1 First output number for the file.
   !> @param[in] num2 Second output number for the file.
   !> @param[out] fout Output string with the file name.
   SUBROUTINE FILE_NAME_1(rootdir, subdir, froot, num1, num2, fout)
      IMPLICIT NONE
      ! Calling arguments.
      CHARACTER(LEN=*),INTENT(IN) :: rootdir, froot
      INTEGER,INTENT(IN) :: subdir, num1, num2
      CHARACTER(LEN=FILE_NAME_LENGTH),INTENT(OUT) :: fout
      ! Local variables.
      ! String used to form the subdirectory name as an integer with no padding.
      CHARACTER(LEN=6) :: numer
      !
      ! Form the subdirectory name keeping the existing format.
      WRITE(numer,10) subdir
      10 FORMAT (I6)
      !
      ! Form the file name.
      WRITE(fout,20) TRIM(rootdir),TRIM(ADJUSTL(numer)),TRIM(froot),num1,num2
      20 FORMAT (A,'/',A,'/',A,'_',I6.6,'_',I2.2,'.h5')
   END SUBROUTINE FILE_NAME_1

   !> Procedure to form filename with root and subdirectories.
   !!
   !! The file name takes the following form:
   !!
   !!    rootdir/subdir/froot_num.h5,
   !!
   !! where num is written with 4 digits.
   !!
   !> @param[in] rootdir Root directory for the file.
   !> @param[in] subdir Subdirectory for the file.
   !> @param[in] froot Root name for the output file.
   !> @param[in] num Output number for the file.
   !> @param[out] fout Output string with the file name.
   SUBROUTINE FILE_NAME_2(rootdir, subdir, froot, num, fout)
      IMPLICIT NONE
      ! Calling arguments.
      CHARACTER(LEN=*),INTENT(IN) :: rootdir, froot
      INTEGER,INTENT(IN) :: subdir, num
      CHARACTER(LEN=FILE_NAME_LENGTH),INTENT(OUT) :: fout
      ! Local variables.
      ! String used to form the subdirectory name as an integer with no padding.
      CHARACTER(LEN=6) :: numer
      !
      ! Form the subdirectory name keeping the existing format.
      WRITE(numer,10) subdir
      10 FORMAT (I6)
      !
      ! Form the file name.
      WRITE(fout,20) TRIM(rootdir),TRIM(ADJUSTL(numer)),TRIM(froot),num
      20 FORMAT (A,'/',A,'/',A,'_',I4.4,'.h5')
   END SUBROUTINE FILE_NAME_2

   !> Procedure to form filename with single root directory.
   !!
   !! The file name takes the following form:
   !!
   !!    rootdir/froot_num.h5
   !!
   !! where num is written with 2 digits.
   !!
   !> @param[in] rootdir Root directory for the file.
   !> @param[in] froot Root name for the output file.
   !> @param[in] num Output number for the file.
   !> @param[out] fout Output string with the file name.
   SUBROUTINE FILE_NAME_3(rootdir, froot, num, fout)
      IMPLICIT NONE
      ! Calling arguments.
      CHARACTER(LEN=*),INTENT(IN) :: rootdir, froot
      INTEGER,INTENT(IN) :: num
      CHARACTER(LEN=FILE_NAME_LENGTH),INTENT(OUT) :: fout
      !
      WRITE(fout,20) TRIM(rootdir),TRIM(froot),num
      20 FORMAT (A,'/',A,'_',I2.2,'.h5')
   END SUBROUTINE FILE_NAME_3

   !> Procedure to form a filename with root directory that is a number.
   !!
   !! The file name takes the following form:
   !!
   !!    subdir/froot_num1_num2.h5,
   !!
   !! where num1 is written with 6 digits, and num2 is written with 2.
   !! Currently, num1 is intended for the process id and num2 for the checkpoint
   !! number.
   !!
   !> @param[in] subdir Directory for the file (integer).
   !> @param[in] froot Root name for the output file.
   !> @param[in] num1 First number for the file (process ID).
   !> @param[in] num2 Second number for the file (checkpoint number).
   !> @param[out] fout Output string with the file name.
   SUBROUTINE FILE_NAME_4(subdir, froot, num1, num2, fout)
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER,INTENT(IN) :: subdir, num1, num2
      CHARACTER(LEN=*),INTENT(IN) :: froot
      CHARACTER(LEN=FILE_NAME_LENGTH),INTENT(OUT) :: fout
      ! Local variables.
      ! String used to form the subdirectory name as an integer with no padding.
      CHARACTER(LEN=6) :: numer
      !
      ! Translate the subdirectory to a string.
      WRITE(numer,10) subdir
      10 FORMAT (I6)
      !
      ! Form the file name.
      WRITE(fout,20) TRIM(ADJUSTL(numer)),TRIM(froot),num1,num2
      20 FORMAT (A,'/',A,'_',I6.6,'_',I2.2,'.h5')
   END SUBROUTINE FILE_NAME_4

   !> Procedure to form a filename with root path, sub path, and number.
   !!
   !! The file name takes the following form:
   !!
   !!    rootdir/subdir/froot_num1.h5
   !!
   !> @param[in] rootdir Root directory for the file.
   !> @param[in] subdir Subdirectory (relative to rootdir) for the file.
   !> @param[in] froot Root name for the file.
   !> @param[in] num Number for the file (filled to 2 integer width).
   !> @param[out] fout Output string with the file name.
   SUBROUTINE FILE_NAME_5(rootdir, subdir, froot, num, fout)
      IMPLICIT NONE
      ! Calling arguments.
      CHARACTER(LEN=*),INTENT(IN) :: rootdir, subdir, froot
      INTEGER,INTENT(IN) :: num
      CHARACTER(LEN=FILE_NAME_LENGTH),INTENT(OUT) :: fout
      !
      ! Form the file name.
      WRITE(fout,20) TRIM(ADJUSTL(rootdir)),TRIM(ADJUSTL(subdir)), &
                     TRIM(ADJUSTL(froot)),num
      20 FORMAT (A,'/',A,'/',A,'_',I2.2,'.h5')
   END SUBROUTINE FILE_NAME_5

   !> Procedure to form filename with single root directory.
   !!
   !! The file name takes the following form:
   !!
   !!    rootdir/froot_num.h5
   !!
   !! where num is written with a user-specified number of digits. This is a
   !! generalized version of FILE_NAME_3.
   !!
   !> @param[in] rootdir Root directory for the file.
   !> @param[in] froot Root name for the output file.
   !> @param[in] num Output number for the file.
   !> @param[in] dig Number of digits to use when writing.
   !> @param[out] fout Output string with the file name.
   SUBROUTINE FILE_NAME_6(rootdir, froot, num, dig, fout)
      IMPLICIT NONE
      ! Calling arguments.
      CHARACTER(LEN=*),INTENT(IN) :: rootdir, froot
      INTEGER,INTENT(IN) :: num,dig
      CHARACTER(LEN=FILE_NAME_LENGTH),INTENT(OUT) :: fout
      ! Local variables.
      ! Used to form the IO format for the number.
      CHARACTER(LEN=16) :: iotmp,iofmt
      ! Buffer to write num with specified number of digits.
      CHARACTER(LEN=32) :: buf
      !
      WRITE(iotmp,25) dig
      WRITE(iofmt,30) TRIM(iotmp),TRIM(iotmp)
      WRITE(buf,iofmt) num
      !
      WRITE(fout,20) TRIM(rootdir),TRIM(froot),TRIM(buf)
      !
      20 FORMAT (A,'/',A,'_',A,'.h5')
      25 FORMAT (I6)
      30 FORMAT ('(I',A,'.',A,')')
   END SUBROUTINE FILE_NAME_6

   !> Procedure to return LOG2 of a number. If the number is not a power of 2,
   !! the function returns zero.
   !!
   !> @param[in] n Number for evaluation.
   !> @return LOG2(n) if n is a power of two, else 0.
   INTEGER FUNCTION POW2(n)
      ! Required modules.
      USE ISO_FORTRAN_ENV,ONLY: REAL32
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER,INTENT(IN) :: n
      ! Local variables.
      ! LOG2 of n as a real.
      REAL(KIND=REAL32) :: p
      !
      ! Calculate the LOG base 2 of the number.
      p = LOG(REAL(n,REAL32))/LOG(2.0_REAL32)
      !
      ! Round to the nearest integer.
      POW2 = NINT(p)
      !
      ! Return 0 if n is not a power of 2.
      IF (2**POW2 .NE. n) THEN
         POW2 = 0
      END IF
   END FUNCTION POW2

   !> Procedure to determine the Fortran precision corresponding to the
   !! precision of a real data field stored in an HDF5 file.
   !!
   !> @param[in] dname Name of the dataset to query.
   !> @param[in] hdf5_id HDF5 ID for the object where the data resides.
   !> @param[out] sglbool Logical set to .TRUE. if the variable is single prec.
   !> @param[out] dblbool Logical set to .TRUE. if the variable is double prec.
   !> @param[in,optional] att_bool Set to .TRUE. if the field to be queried is
   !! an HDF5 attribute instead of a HDF5 dataset.
   SUBROUTINE GET_REAL_PREC(dname, hdf5_id, sglbool, dblbool, att_bool)
      ! Required modules.
      USE ISO_FORTRAN_ENV,ONLY: REAL32, REAL64
      IMPLICIT NONE
      ! Calling arguments.
      CHARACTER(LEN=*),INTENT(IN) :: dname
      INTEGER(KIND=HID_T),INTENT(IN) :: hdf5_id
      LOGICAL,INTENT(OUT) :: sglbool, dblbool
      LOGICAL,OPTIONAL,INTENT(IN) :: att_bool
      ! Local variables.
      ! HDF5 identifier for the dataset or attribute.
      INTEGER(KIND=HID_T) :: dsetid
      ! HDF5 datatype corresponding to the dataset.
      INTEGER(KIND=HID_T) :: typeid
      ! Native Fortran type corresponding to typeid.
      INTEGER(KIND=HID_T) :: fortid
      ! HDF5 precision types corresponding to single and double precision.
      INTEGER(KIND=HID_T) :: hsgl, hdbl
      ! Whether or not to use attribute or dataset calls in HDF5.
      LOGICAL :: att_check
      ! Error handling.
      INTEGER :: ierr
      !
      ! Handle the optional argument specifying if attributes or dataset calls
      ! should be made to HDF5.
      IF (PRESENT(att_bool)) THEN
         IF (att_bool) THEN
            att_check = .TRUE.
         ELSE
            att_check = .FALSE.
         END IF
      ELSE
         att_check = .FALSE.
      END IF
      !
      ! Determine the HDF5 precision types for single and double precision.
      hsgl = H5KIND_TO_TYPE(REAL32, H5_REAL_KIND)
      hdbl = H5KIND_TO_TYPE(REAL64, H5_REAL_KIND)
      !
      ! Open the HDF5 dataset or attribute.
      IF (att_check) THEN
         CALL H5AOPEN_F(hdf5_id, TRIM(dname), dsetid, ierr)
      ELSE
         CALL H5DOPEN_F(hdf5_id, TRIM(dname), dsetid, ierr)
      END IF
      !
      ! Get the datatype associated with the dataset.
      IF (att_check) THEN
         CALL H5AGET_TYPE_F(dsetid, typeid, ierr)
      ELSE
         CALL H5DGET_TYPE_F(dsetid, typeid, ierr)
      END IF
      !
      ! Get the native Fortran type for this datatype.
      CALL H5TGET_NATIVE_TYPE_F(typeid, H5T_DIR_ASCEND_F, fortid, ierr)
      !
      ! Check to see if the precision matches single or double in fortran.
      CALL H5TEQUAL_F(hsgl, fortid, sglbool, ierr)
      CALL H5TEQUAL_F(hdbl, fortid, dblbool, ierr)
      !
      ! Free the type resources.
      CALL H5TCLOSE_F(fortid, ierr)
      CALL H5TCLOSE_F(typeid, ierr)
      !
      ! Close the dataset or attribute.
      IF (att_check) THEN
         CALL H5ACLOSE_F(dsetid, ierr)
      ELSE
         CALL H5DCLOSE_F(dsetid, ierr)
      END IF
   END SUBROUTINE GET_REAL_PREC

   !> Procedure to write 32 bit integer attribute to an HDF5 entity.
   !!
   !> @param[in] hdf5_id HDF5 identifier for the HDF5 item.
   !> @param[in] att_name Name for the attribute.
   !> @param[in] var The data to be written.
   SUBROUTINE WRITE_ATTRIBUTE_INT32(hdf5_id, att_name, var)
      ! Required modules.
      USE ISO_FORTRAN_ENV,ONLY: INT32
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=HID_T),INTENT(IN) :: hdf5_id
      CHARACTER(LEN=*),INTENT(IN) :: att_name
      INTEGER(KIND=INT32),INTENT(INOUT) :: var
      ! Local variables.
      ! HDF5 identifier for the variable type.
      INTEGER(KIND=HID_T) :: h5prec
      ! HDF5 dataspace identifier for the attribute.
      INTEGER(KIND=HID_T) :: aspace_id
      ! HDF5 identifier for the attribute.
      INTEGER(KIND=HID_T) :: attr_id
      ! Dimensions of the data to be written.
      INTEGER(KIND=HSIZE_T),DIMENSION(1),PARAMETER :: dims = [1]
      ! Error handling.
      INTEGER :: ierr
      !
      ! Identify the HDF5 type.
      h5prec = H5KIND_TO_TYPE(INT32, H5_INTEGER_KIND)
      !
      ! Create the dataspace for the attribute.
      CALL H5SCREATE_F(H5S_SCALAR_F, aspace_id, ierr)
      !
      ! Create the attribute.
      CALL H5ACREATE_F(hdf5_id, att_name, h5prec, aspace_id, attr_id, ierr)
      !
      ! Write to the attribute.
      CALL H5AWRITE_F(attr_id, h5prec, var, dims, ierr)
      !
      ! Close the attribute.
      CALL H5ACLOSE_F(attr_id, ierr)
      !
      ! Close the dataspace for the attribute.
      CALL H5SCLOSE_F(aspace_id, ierr)
   END SUBROUTINE WRITE_ATTRIBUTE_INT32

   !> Procedure to write 1D 32 bit integer attribute array to file.
   !!
   !> @param[in] hdf5_id HDF5 identifier for the HDF5 item.
   !> @param[in] att_name Name for the attribute.
   !> @param[in] dims Dimensions of the data.
   !> @param[in] var The data to be written.
   SUBROUTINE WRITE_ATTRIBUTE_INT32_RANK1(hdf5_id, att_name, dims, var)
      ! Required modules.
      USE ISO_FORTRAN_ENV,ONLY: INT32
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=HID_T),INTENT(IN) :: hdf5_id
      CHARACTER(LEN=*),INTENT(IN) :: att_name
      INTEGER(KIND=HSIZE_T),DIMENSION(1),INTENT(IN) :: dims
      INTEGER(KIND=INT32),DIMENSION(dims(1)),INTENT(INOUT) :: var
      ! Local variables.
      ! HDF5 identifier for the variable type.
      INTEGER(KIND=HID_T) :: h5prec
      ! HDF5 dataspace identifier for the attribute.
      INTEGER(KIND=HID_T) :: aspace_id
      ! HDF5 identifier for the attribute.
      INTEGER(KIND=HID_T) :: attr_id
      ! Error handling.
      INTEGER :: ierr
      !
      ! Identify the HDF5 type.
      h5prec = H5KIND_TO_TYPE(INT32, H5_INTEGER_KIND)
      !
      ! Create the dataspace for the attribute.
      CALL H5SCREATE_SIMPLE_F(1, dims, aspace_id, ierr)
      !
      ! Create the attribute.
      CALL H5ACREATE_F(hdf5_id, att_name, h5prec, aspace_id, attr_id, ierr)
      !
      ! Write to the attribute.
      CALL H5AWRITE_F(attr_id, h5prec, var, dims, ierr)
      !
      ! Close the attribute.
      CALL H5ACLOSE_F(attr_id, ierr)
      !
      ! Close the dataspace for the attribute.
      CALL H5SCLOSE_F(aspace_id, ierr)
   END SUBROUTINE WRITE_ATTRIBUTE_INT32_RANK1

   !> Procedure to write 64 bit integer attribute to an HDF5 entity.
   !!
   !> @param[in] hdf5_id HDF5 identifier for the HDF5 item.
   !> @param[in] att_name Name for the attribute.
   !> @param[in] var The data to be written.
#ifdef NEW_HDF5
   SUBROUTINE WRITE_ATTRIBUTE_INT64(hdf5_id, att_name, var)
      ! Required modules.
      USE ISO_FORTRAN_ENV,ONLY: INT64
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=HID_T),INTENT(IN) :: hdf5_id
      CHARACTER(LEN=*),INTENT(IN) :: att_name
      INTEGER(KIND=INT64),INTENT(INOUT) :: var
      ! Local variables.
      ! HDF5 identifier for the variable type.
      INTEGER(KIND=HID_T) :: h5prec
      ! HDF5 dataspace identifier for the attribute.
      INTEGER(KIND=HID_T) :: aspace_id
      ! HDF5 identifier for the attribute.
      INTEGER(KIND=HID_T) :: attr_id
      ! Dimensions of the data to be written.
      INTEGER(KIND=HSIZE_T),DIMENSION(1),PARAMETER :: dims = [1]
      ! Error handling.
      INTEGER :: ierr
      !
      ! Identify the HDF5 type.
      h5prec = H5KIND_TO_TYPE(INT64, H5_INTEGER_KIND)
      !
      ! Create the dataspace for the attribute.
      CALL H5SCREATE_F(H5S_SCALAR_F, aspace_id, ierr)
      !
      ! Create the attribute.
      CALL H5ACREATE_F(hdf5_id, att_name, h5prec, aspace_id, attr_id, ierr)
      !
      ! Write to the attribute.
      CALL H5AWRITE_F(attr_id, h5prec, var, dims, ierr)
      !
      ! Close the attribute.
      CALL H5ACLOSE_F(attr_id, ierr)
      !
      ! Close the dataspace for the attribute.
      CALL H5SCLOSE_F(aspace_id, ierr)
   END SUBROUTINE WRITE_ATTRIBUTE_INT64
#endif

   !> Procedure to write 1D 64 bit integer attribute array to file.
   !!
   !> @param[in] hdf5_id HDF5 identifier for the HDF5 item.
   !> @param[in] att_name Name for the attribute.
   !> @param[in] dims Dimensions of the data.
   !> @param[in] var The data to be written.
#ifdef NEW_HDF5
   SUBROUTINE WRITE_ATTRIBUTE_INT64_RANK1(hdf5_id, att_name, dims, var)
      ! Required modules.
      USE ISO_FORTRAN_ENV,ONLY: INT64
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=HID_T),INTENT(IN) :: hdf5_id
      CHARACTER(LEN=*),INTENT(IN) :: att_name
      INTEGER(KIND=HSIZE_T),DIMENSION(1),INTENT(IN) :: dims
      INTEGER(KIND=INT64),DIMENSION(dims(1)),INTENT(INOUT) :: var
      ! Local variables.
      ! HDF5 identifier for the variable type.
      INTEGER(KIND=HID_T) :: h5prec
      ! HDF5 dataspace identifier for the attribute.
      INTEGER(KIND=HID_T) :: aspace_id
      ! HDF5 identifier for the attribute.
      INTEGER(KIND=HID_T) :: attr_id
      ! Error handling.
      INTEGER :: ierr
      !
      ! Identify the HDF5 type.
      h5prec = H5KIND_TO_TYPE(INT64, H5_INTEGER_KIND)
      !
      ! Create the dataspace for the attribute.
      CALL H5SCREATE_SIMPLE_F(1, dims, aspace_id, ierr)
      !
      ! Create the attribute.
      CALL H5ACREATE_F(hdf5_id, att_name, h5prec, aspace_id, attr_id, ierr)
      !
      ! Write to the attribute.
      CALL H5AWRITE_F(attr_id, h5prec, var, dims, ierr)
      !
      ! Close the attribute.
      CALL H5ACLOSE_F(attr_id, ierr)
      !
      ! Close the dataspace for the attribute.
      CALL H5SCLOSE_F(aspace_id, ierr)
   END SUBROUTINE WRITE_ATTRIBUTE_INT64_RANK1
#endif

   !> Procedure to write a single precision real attribute to an HDF5 entity.
   !!
   !> @param[in] hdf5_id HDF5 identifier for the HDF5 item.
   !> @param[in] att_name Name for the attribute.
   !> @param[in] var The data to be written.
   SUBROUTINE WRITE_ATTRIBUTE_REAL32(hdf5_id, att_name, var)
      ! Required modules.
      USE ISO_FORTRAN_ENV,ONLY: REAL32
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=HID_T),INTENT(IN) :: hdf5_id
      CHARACTER(LEN=*),INTENT(IN) :: att_name
      REAL(KIND=REAL32),INTENT(IN) :: var
      ! Local variables.
      ! HDF5 identifier for the variable type.
      INTEGER(KIND=HID_T) :: h5prec
      ! HDF5 dataspace identifier for the attribute.
      INTEGER(KIND=HID_T) :: aspace_id
      ! HDF5 identifier for the attribute.
      INTEGER(KIND=HID_T) :: attr_id
      ! Dimensions of the data to be written.
      INTEGER(KIND=HSIZE_T),DIMENSION(1),PARAMETER :: dims = [1]
      ! Error handling.
      INTEGER :: ierr
      !
      ! Identify the HDF5 type.
      h5prec = H5KIND_TO_TYPE(REAL32, H5_REAL_KIND)
      !
      ! Create the dataspace for the attribute.
      CALL H5SCREATE_F(H5S_SCALAR_F, aspace_id, ierr)
      !
      ! Create the attribute.
      CALL H5ACREATE_F(hdf5_id, att_name, h5prec, aspace_id, attr_id, ierr)
      !
      ! Write to the attribute.
      CALL H5AWRITE_F(attr_id, h5prec, var, dims, ierr)
      !
      ! Close the attribute.
      CALL H5ACLOSE_F(attr_id, ierr)
      !
      ! Close the dataspace for the attribute.
      CALL H5SCLOSE_F(aspace_id, ierr)
   END SUBROUTINE WRITE_ATTRIBUTE_REAL32

   !> Procedure to write 1D single precision real attribute to file.
   !!
   !> @param[in] hdf5_id HDF5 identifier for the HDF5 item.
   !> @param[in] att_name Name for the attribute.
   !> @param[in] dims Dimensions of the data.
   !> @param[in] var The data to be written.
   SUBROUTINE WRITE_ATTRIBUTE_REAL32_RANK1(hdf5_id, att_name, dims, var)
      ! Required modules.
      USE ISO_FORTRAN_ENV,ONLY: REAL32
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=HID_T),INTENT(IN) :: hdf5_id
      CHARACTER(LEN=*),INTENT(IN) :: att_name
      INTEGER(KIND=HSIZE_T),DIMENSION(1),INTENT(IN) :: dims
      REAL(KIND=REAL32),DIMENSION(dims(1)),INTENT(IN) :: var
      ! Local variables.
      ! HDF5 identifier for the variable type.
      INTEGER(KIND=HID_T) :: h5prec
      ! HDF5 dataspace identifier for the attribute.
      INTEGER(KIND=HID_T) :: aspace_id
      ! HDF5 identifier for the attribute.
      INTEGER(KIND=HID_T) :: attr_id
      ! Error handling.
      INTEGER :: ierr
      !
      ! Identify the HDF5 type.
      h5prec = H5KIND_TO_TYPE(REAL32, H5_REAL_KIND)
      !
      ! Create the dataspace for the attribute.
      CALL H5SCREATE_SIMPLE_F(1, dims, aspace_id, ierr)
      !
      ! Create the attribute.
      CALL H5ACREATE_F(hdf5_id, att_name, h5prec, aspace_id, attr_id, ierr)
      !
      ! Write to the attribute.
      CALL H5AWRITE_F(attr_id, h5prec, var, dims, ierr)
      !
      ! Close the attribute.
      CALL H5ACLOSE_F(attr_id, ierr)
      !
      ! Close the dataspace for the attribute.
      CALL H5SCLOSE_F(aspace_id, ierr)
   END SUBROUTINE WRITE_ATTRIBUTE_REAL32_RANK1

   !> Procedure to write a double precision real attribute to an HDF5 entity.
   !!
   !> @param[in] hdf5_id HDF5 identifier for the HDF5 item.
   !> @param[in] att_name Name for the attribute.
   !> @param[in] var The data to be written.
   SUBROUTINE WRITE_ATTRIBUTE_REAL64(hdf5_id, att_name, var)
      ! Required modules.
      USE ISO_FORTRAN_ENV,ONLY: REAL64
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=HID_T),INTENT(IN) :: hdf5_id
      CHARACTER(LEN=*),INTENT(IN) :: att_name
      REAL(KIND=REAL64),INTENT(IN) :: var
      ! Local variables.
      ! HDF5 identifier for the variable type.
      INTEGER(KIND=HID_T) :: h5prec
      ! HDF5 dataspace identifier for the attribute.
      INTEGER(KIND=HID_T) :: aspace_id
      ! HDF5 identifier for the attribute.
      INTEGER(KIND=HID_T) :: attr_id
      ! Dimensions of the data to be written.
      INTEGER(KIND=HSIZE_T),DIMENSION(1),PARAMETER :: dims = [1]
      ! Error handling.
      INTEGER :: ierr
      !
      ! Identify the HDF5 type.
      h5prec = H5KIND_TO_TYPE(REAL64, H5_REAL_KIND)
      !
      ! Create the dataspace for the attribute.
      CALL H5SCREATE_F(H5S_SCALAR_F, aspace_id, ierr)
      !
      ! Create the attribute.
      CALL H5ACREATE_F(hdf5_id, att_name, h5prec, aspace_id, attr_id, ierr)
      !
      ! Write to the attribute.
      CALL H5AWRITE_F(attr_id, h5prec, var, dims, ierr)
      !
      ! Close the attribute.
      CALL H5ACLOSE_F(attr_id, ierr)
      !
      ! Close the dataspace for the attribute.
      CALL H5SCLOSE_F(aspace_id, ierr)
   END SUBROUTINE WRITE_ATTRIBUTE_REAL64

   !> Procedure to write 1D double precision real attribute to file.
   !!
   !> @param[in] hdf5_id HDF5 identifier for the HDF5 item.
   !> @param[in] att_name Name for the attribute.
   !> @param[in] dims Dimensions of the data.
   !> @param[in] var The data to be written.
   SUBROUTINE WRITE_ATTRIBUTE_REAL64_RANK1(hdf5_id, att_name, dims, var)
      ! Required modules.
      USE ISO_FORTRAN_ENV,ONLY: REAL64
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=HID_T),INTENT(IN) :: hdf5_id
      CHARACTER(LEN=*),INTENT(IN) :: att_name
      INTEGER(KIND=HSIZE_T),DIMENSION(1),INTENT(IN) :: dims
      REAL(KIND=REAL64),DIMENSION(dims(1)),INTENT(IN) :: var
      ! Local variables.
      ! HDF5 identifier for the variable type.
      INTEGER(KIND=HID_T) :: h5prec
      ! HDF5 dataspace identifier for the attribute.
      INTEGER(KIND=HID_T) :: aspace_id
      ! HDF5 identifier for the attribute.
      INTEGER(KIND=HID_T) :: attr_id
      ! Error handling.
      INTEGER :: ierr
      !
      ! Identify the HDF5 type.
      h5prec = H5KIND_TO_TYPE(REAL64, H5_REAL_KIND)
      !
      ! Create the dataspace for the attribute.
      CALL H5SCREATE_SIMPLE_F(1, dims, aspace_id, ierr)
      !
      ! Create the attribute.
      CALL H5ACREATE_F(hdf5_id, att_name, h5prec, aspace_id, attr_id, ierr)
      !
      ! Write to the attribute.
      CALL H5AWRITE_F(attr_id, h5prec, var, dims, ierr)
      !
      ! Close the attribute.
      CALL H5ACLOSE_F(attr_id, ierr)
      !
      ! Close the dataspace for the attribute.
      CALL H5SCLOSE_F(aspace_id, ierr)
   END SUBROUTINE WRITE_ATTRIBUTE_REAL64_RANK1

   !> Procedure to write a character attribute to an HDF5 entity.
   !!
   !> @param[in] hdf5_id HDF5 identifier for the HDF5 item.
   !> @param[in] att_name Name for the attribute.
   !> @param[in] var The data to be written.
   SUBROUTINE WRITE_ATTRIBUTE_CHAR(hdf5_id, att_name, var)
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=HID_T),INTENT(IN) :: hdf5_id
      CHARACTER(LEN=*),INTENT(IN) :: att_name
      CHARACTER(LEN=*),INTENT(IN) :: var
      ! Local variables.
      ! HDF5 dataspace identifier for the attribute.
      INTEGER(KIND=HID_T) :: aspace_id
      ! HDF5 identifier for the attribute.
      INTEGER(KIND=HID_T) :: attr_id
      ! HDF5 identifier for the data type.
      INTEGER(KIND=HID_T) :: type_id
      ! Dimensions of the data to be written.
      INTEGER(KIND=HSIZE_T),DIMENSION(1),PARAMETER :: dims = [1]
      ! Length of the character string.
      INTEGER(KIND=SIZE_T) :: strlen
      ! Error handling.
      INTEGER :: ierr
      !
      ! Create the dataspace for the attribute.
      CALL H5SCREATE_F(H5S_SCALAR_F, aspace_id, ierr)
      !
      ! Get the length of the string figured out.
      strlen = INT(LEN_TRIM(var), HSIZE_T)
      CALL H5TCOPY_F(H5T_NATIVE_CHARACTER, type_id, ierr)
      CALL H5TSET_SIZE_F(type_id, strlen, ierr)
      !
      ! Create the attribute.
      CALL H5ACREATE_F(hdf5_id, att_name, type_id, aspace_id, attr_id, ierr)
      !
      ! Write to the attribute.
      CALL H5AWRITE_F(attr_id, type_id, TRIM(var), dims, ierr)
      !
      ! Close the attribute.
      CALL H5ACLOSE_F(attr_id, ierr)
      !
      ! Close the dataspace for the attribute.
      CALL H5SCLOSE_F(aspace_id, ierr)
   END SUBROUTINE WRITE_ATTRIBUTE_CHAR

   !> Procedure to read an integer attribute into a 32bit buffer.
   !!
   !> @param[in] hdf5_id HDF5 identifier for the HDF5 item.
   !> @param[in] att_name Name of the attribute.
   !> @param[out] var Variable to store the attribute.
   SUBROUTINE READ_ATTRIBUTE_INT32(hdf5_id, att_name, var)
      ! Required modules.
      USE ISO_FORTRAN_ENV,ONLY: INT32
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=HID_T),INTENT(IN) :: hdf5_id
      CHARACTER(LEN=*),INTENT(IN) :: att_name
      INTEGER(KIND=INT32),INTENT(OUT) :: var
      ! Local variables.
      ! HDF5 identifier for the variable type.
      INTEGER(KIND=HID_T) :: h5prec
      ! HDF5 identifier for the attribute.
      INTEGER(KIND=HID_T) :: att_id
      ! Dimensions of the variable being read.
      INTEGER(KIND=HSIZE_T),DIMENSION(1) :: dims = [1]
      ! Error handling.
      INTEGER :: ierr
      !
      ! Identify the HDF5 type.
      h5prec = H5KIND_TO_TYPE(INT32, H5_INTEGER_KIND)
      !
      ! Open the attribute.
      CALL H5AOPEN_F(hdf5_id, att_name, att_id, ierr)
      !
      ! Read the attribute.
      CALL H5AREAD_F(att_id, h5prec, var, dims, ierr)
      !
      ! Close the attribute.
      CALL H5ACLOSE_F(att_id, ierr)
   END SUBROUTINE READ_ATTRIBUTE_INT32

   !> Procedure to read an integer attribute array into 32bit integer buffer.
   !!
   !> @param[in] hdf5_id HDF5 identifier for the HDF5 item.
   !> @param[in] att_name Name of the attribute.
   !> @param[in] dims Size of the variable buffer to read the attribute.
   !> @param[out] var Variable to store the attribute.
   SUBROUTINE READ_ATTRIBUTE_INT32_RANK1(hdf5_id, att_name, dims, var)
      ! Required modules.
      USE ISO_FORTRAN_ENV,ONLY: INT32
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=HID_T),INTENT(IN) :: hdf5_id
      CHARACTER(LEN=*),INTENT(IN) :: att_name
      INTEGER(KIND=HSIZE_T),DIMENSION(1),INTENT(IN) :: dims
      INTEGER(KIND=INT32),DIMENSION(dims(1)),INTENT(OUT) :: var
      ! Local variables.
      ! HDF5 identifier for the variable type.
      INTEGER(KIND=HID_T) :: h5prec
      ! HDF5 identifier for the attribute.
      INTEGER(KIND=HID_T) :: att_id
      ! Error handling.
      INTEGER :: ierr
      !
      ! Identify the HDF5 type.
      h5prec = H5KIND_TO_TYPE(INT32, H5_INTEGER_KIND)
      !
      ! Open the attribute.
      CALL H5AOPEN_F(hdf5_id, att_name, att_id, ierr)
      !
      ! Read the attribute.
      CALL H5AREAD_F(att_id, h5prec, var, dims, ierr)
      !
      ! Close the attribute.
      CALL H5ACLOSE_F(att_id, ierr)
   END SUBROUTINE READ_ATTRIBUTE_INT32_RANK1

   !> Procedure to read an integer attribute into a 64bit buffer.
   !!
   !> @param[in] hdf5_id HDF5 identifier for the HDF5 item.
   !> @param[in] att_name Name of the attribute.
   !> @param[out] var Variable to store the attribute.
#ifdef NEW_HDF5
   SUBROUTINE READ_ATTRIBUTE_INT64(hdf5_id, att_name, var)
      ! Required modules.
      USE ISO_FORTRAN_ENV,ONLY: INT64
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=HID_T),INTENT(IN) :: hdf5_id
      CHARACTER(LEN=*),INTENT(IN) :: att_name
      INTEGER(KIND=INT64),INTENT(OUT) :: var
      ! Local variables.
      ! HDF5 identifier for the variable type.
      INTEGER(KIND=HID_T) :: h5prec
      ! HDF5 identifier for the attribute.
      INTEGER(KIND=HID_T) :: att_id
      ! Dimensions of the variable being read.
      INTEGER(KIND=HSIZE_T),DIMENSION(1) :: dims = [1]
      ! Error handling.
      INTEGER :: ierr
      !
      ! Identify the HDF5 type.
      h5prec = H5KIND_TO_TYPE(INT64, H5_INTEGER_KIND)
      !
      ! Open the attribute.
      CALL H5AOPEN_F(hdf5_id, att_name, att_id, ierr)
      !
      ! Read the attribute.
      CALL H5AREAD_F(att_id, h5prec, var, dims, ierr)
      !
      ! Close the attribute.
      CALL H5ACLOSE_F(att_id, ierr)
   END SUBROUTINE READ_ATTRIBUTE_INT64
#endif

   !> Procedure to read an integer attribute array into 64bit integer buffer.
   !!
   !> @param[in] hdf5_id HDF5 identifier for the HDF5 item.
   !> @param[in] att_name Name of the attribute.
   !> @param[in] dims Size of the variable buffer to read the attribute.
   !> @param[out] var Variable to store the attribute.
#ifdef NEW_HDF5
   SUBROUTINE READ_ATTRIBUTE_INT64_RANK1(hdf5_id, att_name, dims, var)
      ! Required modules.
      USE ISO_FORTRAN_ENV,ONLY: INT64
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=HID_T),INTENT(IN) :: hdf5_id
      CHARACTER(LEN=*),INTENT(IN) :: att_name
      INTEGER(KIND=HSIZE_T),DIMENSION(1),INTENT(IN) :: dims
      INTEGER(KIND=INT64),DIMENSION(dims(1)),INTENT(OUT) :: var
      ! Local variables.
      ! HDF5 identifier for the variable type.
      INTEGER(KIND=HID_T) :: h5prec
      ! HDF5 identifier for the attribute.
      INTEGER(KIND=HID_T) :: att_id
      ! Error handling.
      INTEGER :: ierr
      !
      ! Identify the HDF5 type.
      h5prec = H5KIND_TO_TYPE(INT64, H5_INTEGER_KIND)
      !
      ! Open the attribute.
      CALL H5AOPEN_F(hdf5_id, att_name, att_id, ierr)
      !
      ! Read the attribute.
      CALL H5AREAD_F(att_id, h5prec, var, dims, ierr)
      !
      ! Close the attribute.
      CALL H5ACLOSE_F(att_id, ierr)
   END SUBROUTINE READ_ATTRIBUTE_INT64_RANK1
#endif

   !> Procedure to read a single precision real attribute.
   !!
   !> @param[in] hdf5_id HDF5 identifier for the HDF5 item.
   !> @param[in] att_name Name of the attribute.
   !> @param[out] var Variable to store the attribute.
   SUBROUTINE READ_ATTRIBUTE_REAL32(hdf5_id, att_name, var)
      ! Required modules.
      USE ISO_FORTRAN_ENV,ONLY: REAL32
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=HID_T),INTENT(IN) :: hdf5_id
      CHARACTER(LEN=*),INTENT(IN) :: att_name
      REAL(KIND=REAL32),INTENT(OUT) :: var
      ! Local variables.
      ! HDF5 identifier for the variable type.
      INTEGER(KIND=HID_T) :: h5prec
      ! HDF5 identifier for the attribute.
      INTEGER(KIND=HID_T) :: att_id
      ! Dimensions of the variable being read.
      INTEGER(KIND=HSIZE_T),DIMENSION(1) :: dims = [1]
      ! Error handling.
      INTEGER :: ierr
      !
      ! Identify the HDF5 type.
      h5prec = H5KIND_TO_TYPE(REAL32, H5_REAL_KIND)
      !
      ! Open the attribute.
      CALL H5AOPEN_F(hdf5_id, att_name, att_id, ierr)
      !
      ! Read the attribute.
      CALL H5AREAD_F(att_id, h5prec, var, dims, ierr)
      !
      ! Close the attribute.
      CALL H5ACLOSE_F(att_id, ierr)
   END SUBROUTINE READ_ATTRIBUTE_REAL32

   !> Procedure to read a single precision real attribute array.
   !!
   !> @param[in] hdf5_id HDF5 identifier for the HDF5 item.
   !> @param[in] att_name Name of the attribute.
   !> @param[in] dims Size of the variable buffer to read the attribute.
   !> @param[out] var Variable to store the attribute.
   SUBROUTINE READ_ATTRIBUTE_REAL32_RANK1(hdf5_id, att_name, dims, var)
      ! Required modules.
      USE ISO_FORTRAN_ENV,ONLY: REAL32
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=HID_T),INTENT(IN) :: hdf5_id
      CHARACTER(LEN=*),INTENT(IN) :: att_name
      INTEGER(KIND=HSIZE_T),DIMENSION(1),INTENT(IN) :: dims
      REAL(KIND=REAL32),DIMENSION(dims(1)),INTENT(OUT) :: var
      ! Local variables.
      ! HDF5 identifier for the variable type.
      INTEGER(KIND=HID_T) :: h5prec
      ! HDF5 identifier for the attribute.
      INTEGER(KIND=HID_T) :: att_id
      ! Error handling.
      INTEGER :: ierr
      !
      ! Identify the HDF5 type.
      h5prec = H5KIND_TO_TYPE(REAL32, H5_REAL_KIND)
      !
      ! Open the attribute.
      CALL H5AOPEN_F(hdf5_id, att_name, att_id, ierr)
      !
      ! Read the attribute.
      CALL H5AREAD_F(att_id, h5prec, var, dims, ierr)
      !
      ! Close the attribute.
      CALL H5ACLOSE_F(att_id, ierr)
   END SUBROUTINE READ_ATTRIBUTE_REAL32_RANK1

   !> Procedure to read a double precision real attribute.
   !!
   !> @param[in] hdf5_id HDF5 identifier for the HDF5 item.
   !> @param[in] att_name Name of the attribute.
   !> @param[out] var Variable to store the attribute.
   SUBROUTINE READ_ATTRIBUTE_REAL64(hdf5_id, att_name, var)
      ! Required modules.
      USE ISO_FORTRAN_ENV,ONLY: REAL64
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=HID_T),INTENT(IN) :: hdf5_id
      CHARACTER(LEN=*),INTENT(IN) :: att_name
      REAL(KIND=REAL64),INTENT(OUT) :: var
      ! Local variables.
      ! HDF5 identifier for the variable type.
      INTEGER(KIND=HID_T) :: h5prec
      ! HDF5 identifier for the attribute.
      INTEGER(KIND=HID_T) :: att_id
      ! Dimensions of the variable being read.
      INTEGER(KIND=HSIZE_T),DIMENSION(1) :: dims = [1]
      ! Error handling.
      INTEGER :: ierr
      !
      ! Identify the HDF5 type.
      h5prec = H5KIND_TO_TYPE(REAL64, H5_REAL_KIND)
      !
      ! Open the attribute.
      CALL H5AOPEN_F(hdf5_id, att_name, att_id, ierr)
      !
      ! Read the attribute.
      CALL H5AREAD_F(att_id, h5prec, var, dims, ierr)
      !
      ! Close the attribute.
      CALL H5ACLOSE_F(att_id, ierr)
   END SUBROUTINE READ_ATTRIBUTE_REAL64

   !> Procedure to read a double precision real attribute array.
   !!
   !> @param[in] hdf5_id HDF5 identifier for the HDF5 item.
   !> @param[in] att_name Name of the attribute.
   !> @param[in] dims Size of the variable buffer to read the attribute.
   !> @param[out] var Variable to store the attribute.
   SUBROUTINE READ_ATTRIBUTE_REAL64_RANK1(hdf5_id, att_name, dims, var)
      ! Required modules.
      USE ISO_FORTRAN_ENV,ONLY: REAL64
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=HID_T),INTENT(IN) :: hdf5_id
      CHARACTER(LEN=*),INTENT(IN) :: att_name
      INTEGER(KIND=HSIZE_T),DIMENSION(1),INTENT(IN) :: dims
      REAL(KIND=REAL64),DIMENSION(dims(1)),INTENT(OUT) :: var
      ! Local variables.
      ! HDF5 identifier for the variable type.
      INTEGER(KIND=HID_T) :: h5prec
      ! HDF5 identifier for the attribute.
      INTEGER(KIND=HID_T) :: att_id
      ! Error handling.
      INTEGER :: ierr
      !
      ! Identify the HDF5 type.
      h5prec = H5KIND_TO_TYPE(REAL64, H5_REAL_KIND)
      !
      ! Open the attribute.
      CALL H5AOPEN_F(hdf5_id, att_name, att_id, ierr)
      !
      ! Read the attribute.
      CALL H5AREAD_F(att_id, h5prec, var, dims, ierr)
      !
      ! Close the attribute.
      CALL H5ACLOSE_F(att_id, ierr)
   END SUBROUTINE READ_ATTRIBUTE_REAL64_RANK1

   !> Procedure to write 1D 32bit integer data to file.
   !!
   !! This routine is designed to work in the case that parallel HDF5 is being
   !! used, but the user only wants one process to write to the dataset. Since
   !! dataset creation is collective, all processes participate in making the
   !! dataset, but only the desired process writes to it (assuming the optional
   !! arguments tid and wid are present).
   !!
   !> @param[in] dims Dimensions of the data.
   !> @param[in] dname Output name for the data.
   !> @param[in] dat Array with dimensions given by dims.
   !> @param[in] hdf5_id HDF5 ID for the object where data will be written.
   !> @param[in,optional] tid Process identifier (rank).
   !> @param[in,optional] wid Process ID to write out the data.
   SUBROUTINE WRITE_DATA_INT32_RANK1(dims, dname, dat, hdf5_id, tid, wid)
      ! Required modules.
      USE ISO_FORTRAN_ENV,ONLY: INT32
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=HSIZE_T),DIMENSION(1),INTENT(IN) :: dims
      CHARACTER(LEN=*),INTENT(IN) :: dname
      INTEGER(KIND=INT32),DIMENSION(dims(1)),INTENT(IN) :: dat
      INTEGER(KIND=HID_T),INTENT(IN) :: hdf5_id
      INTEGER,OPTIONAL,INTENT(IN) :: tid, wid
      ! Local variables.
      ! HDF5 identifier for the variable type.
      INTEGER(KIND=HID_T) :: h5prec
      ! HDF5 identifier for the dataspace and dataset.
      INTEGER(KIND=HID_T) :: dataspace_id, dataset_id
      ! Error handling.
      INTEGER :: ierr
      !
      ! Identify the HDF5 type.
      h5prec = H5KIND_TO_TYPE(INT32, H5_INTEGER_KIND)
      !
      ! Create the dataspace for the data set.
      CALL H5SCREATE_SIMPLE_F(1, dims, dataspace_id, ierr)
      !
      ! Create a dataset in hdf5_id with default properties.
      CALL H5DCREATE_F(hdf5_id, dname, h5prec, dataspace_id, dataset_id, ierr)
      !
      ! Write to the dataset.
      IF (PRESENT(tid) .AND. PRESENT(wid)) THEN
         IF (tid .EQ. wid) THEN
            CALL H5DWRITE_F(dataset_id, h5prec, dat, dims, ierr)
         END IF
      ELSE
         CALL H5DWRITE_F(dataset_id, h5prec, dat, dims, ierr)
      END IF
      !
      ! Close the dataset.
      CALL H5DCLOSE_F(dataset_id, ierr)
      !
      ! Close the dataspace.
      CALL H5SCLOSE_F(dataspace_id, ierr)
   END SUBROUTINE WRITE_DATA_INT32_RANK1

   !> Procedure to write 2D 32bit integer data to file.
   !!
   !> @param[in] dims Dimensions of the data.
   !> @param[in] dname Output name for the data.
   !> @param[in] dat Array with dimensions given by dims.
   !> @param[in] hdf5_id HDF5 ID for the object where data will be written.
   SUBROUTINE WRITE_DATA_INT32_RANK2(dims, dname, dat, hdf5_id)
      ! Required modules.
      USE ISO_FORTRAN_ENV,ONLY: INT32
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=HSIZE_T),DIMENSION(2),INTENT(IN) :: dims
      CHARACTER(LEN=*),INTENT(IN) :: dname
      INTEGER(KIND=INT32),DIMENSION(dims(1),dims(2)),INTENT(IN) :: dat
      INTEGER(KIND=HID_T),INTENT(IN) :: hdf5_id
      ! Local variables.
      ! HDF5 identifier for the variable type.
      INTEGER(KIND=HID_T) :: h5prec
      ! HDF5 identifier for the dataspace and dataset.
      INTEGER(KIND=HID_T) :: dataspace_id, dataset_id
      ! Error handling.
      INTEGER :: ierr
      !
      ! Identify the HDF5 type.
      h5prec = H5KIND_TO_TYPE(INT32, H5_INTEGER_KIND)
      !
      ! Create the dataspace for the data set.
      CALL H5SCREATE_SIMPLE_F(2, dims, dataspace_id, ierr)
      !
      ! Create a dataset in hdf5_id with default properties.
      CALL H5DCREATE_F(hdf5_id, dname, h5prec, dataspace_id, dataset_id, ierr)
      !
      ! Write to the dataset.
      CALL H5DWRITE_F(dataset_id, h5prec, dat, dims, ierr)
      !
      ! Close the dataset.
      CALL H5DCLOSE_F(dataset_id, ierr)
      !
      ! Close the dataspace.
      CALL H5SCLOSE_F(dataspace_id, ierr)
   END SUBROUTINE WRITE_DATA_INT32_RANK2

   !> Procedure to write 1D 64bit integer data to file.
   !!
   !! This routine is designed to work in the case that parallel HDF5 is being
   !! used, but the user only wants one process to write to the dataset. Since
   !! dataset creation is collective, all processes participate in making the
   !! dataset, but only the desired process writes to it (assuming the optional
   !! arguments tid and wid are present).
   !!
   !> @param[in] dims Dimensions of the data.
   !> @param[in] dname Output name for the data.
   !> @param[in] dat Array with dimensions given by dims.
   !> @param[in] hdf5_id HDF5 ID for the object where data will be written.
   !> @param[in,optional] tid Process identifier (rank).
   !> @param[in,optional] wid Process ID to write out the data.
#ifdef NEW_HDF5
   SUBROUTINE WRITE_DATA_INT64_RANK1(dims, dname, dat, hdf5_id, tid, wid)
      ! Required modules.
      USE ISO_FORTRAN_ENV,ONLY: INT64
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=HSIZE_T),DIMENSION(1),INTENT(IN) :: dims
      CHARACTER(LEN=*),INTENT(IN) :: dname
      INTEGER(KIND=INT64),DIMENSION(dims(1)),INTENT(IN) :: dat
      INTEGER(KIND=HID_T),INTENT(IN) :: hdf5_id
      INTEGER,OPTIONAL,INTENT(IN) :: tid, wid
      ! Local variables.
      ! HDF5 identifier for the variable type.
      INTEGER(KIND=HID_T) :: h5prec
      ! HDF5 identifier for the dataspace and dataset.
      INTEGER(KIND=HID_T) :: dataspace_id, dataset_id
      ! Error handling.
      INTEGER :: ierr
      !
      ! Identify the HDF5 type.
      h5prec = H5KIND_TO_TYPE(INT64, H5_INTEGER_KIND)
      !
      ! Create the dataspace for the data set.
      CALL H5SCREATE_SIMPLE_F(1, dims, dataspace_id, ierr)
      !
      ! Create a dataset in hdf5_id with default properties.
      CALL H5DCREATE_F(hdf5_id, dname, h5prec, dataspace_id, dataset_id, ierr)
      !
      ! Write to the dataset.
      IF (PRESENT(tid) .AND. PRESENT(wid)) THEN
         IF (tid .EQ. wid) THEN
            CALL H5DWRITE_F(dataset_id, h5prec, dat, dims, ierr)
         END IF
      ELSE
         CALL H5DWRITE_F(dataset_id, h5prec, dat, dims, ierr)
      END IF
      !
      ! Close the dataset.
      CALL H5DCLOSE_F(dataset_id, ierr)
      !
      ! Close the dataspace.
      CALL H5SCLOSE_F(dataspace_id, ierr)
   END SUBROUTINE WRITE_DATA_INT64_RANK1
#endif

   !> Procedure to write 2D 64bit integer data to file.
   !!
   !> @param[in] dims Dimensions of the data.
   !> @param[in] dname Output name for the data.
   !> @param[in] dat Array with dimensions given by dims.
   !> @param[in] hdf5_id HDF5 ID for the object where data will be written.
#ifdef NEW_HDF5
   SUBROUTINE WRITE_DATA_INT64_RANK2(dims, dname, dat, hdf5_id)
      ! Required modules.
      USE ISO_FORTRAN_ENV,ONLY: INT64
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=HSIZE_T),DIMENSION(2),INTENT(IN) :: dims
      CHARACTER(LEN=*),INTENT(IN) :: dname
      INTEGER(KIND=INT64),DIMENSION(dims(1),dims(2)),INTENT(IN) :: dat
      INTEGER(KIND=HID_T),INTENT(IN) :: hdf5_id
      ! Local variables.
      ! HDF5 identifier for the variable type.
      INTEGER(KIND=HID_T) :: h5prec
      ! HDF5 identifier for the dataspace and dataset.
      INTEGER(KIND=HID_T) :: dataspace_id, dataset_id
      ! Error handling.
      INTEGER :: ierr
      !
      ! Identify the HDF5 type.
      h5prec = H5KIND_TO_TYPE(INT64, H5_INTEGER_KIND)
      !
      ! Create the dataspace for the data set.
      CALL H5SCREATE_SIMPLE_F(2, dims, dataspace_id, ierr)
      !
      ! Create a dataset in hdf5_id with default properties.
      CALL H5DCREATE_F(hdf5_id, dname, h5prec, dataspace_id, dataset_id, ierr)
      !
      ! Write to the dataset.
      CALL H5DWRITE_F(dataset_id, h5prec, dat, dims, ierr)
      !
      ! Close the dataset.
      CALL H5DCLOSE_F(dataset_id, ierr)
      !
      ! Close the dataspace.
      CALL H5SCLOSE_F(dataspace_id, ierr)
   END SUBROUTINE WRITE_DATA_INT64_RANK2
#endif

   !> Procedure to write 1D single precision real data to file.
   !!
   !! This routine is designed to work in the case that parallel HDF5 is being
   !! used, but the user only wants one process to write to the dataset. Since
   !! dataset creation is collective, all processes participate in making the
   !! dataset, but only the desired process writes to it (assuming the optional
   !! arguments tid and wid are present).
   !!
   !> @param[in] dims Dimensions of the data.
   !> @param[in] dname Output name for the data.
   !> @param[in] dat Array with dimensions given by dims.
   !> @param[in] hdf5_id HDF5 ID for the object where data will be written.
   !> @param[in,optional] tid Process identifier (rank).
   !> @param[in,optional] wid Process ID to write out the data.
   SUBROUTINE WRITE_DATA_REAL32_RANK1(dims, dname, dat, hdf5_id, tid, wid)
      ! Required modules.
      USE ISO_FORTRAN_ENV,ONLY: REAL32
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=HSIZE_T),DIMENSION(1),INTENT(IN) :: dims
      CHARACTER(LEN=*),INTENT(IN) :: dname
      REAL(KIND=REAL32),DIMENSION(dims(1)),INTENT(IN) :: dat
      INTEGER(KIND=HID_T),INTENT(IN) :: hdf5_id
      INTEGER,OPTIONAL,INTENT(IN) :: tid, wid
      ! Local variables.
      ! HDF5 identifier for the variable type.
      INTEGER(KIND=HID_T) :: h5prec
      ! HDF5 identifier for the dataspace and dataset.
      INTEGER(KIND=HID_T) :: dataspace_id, dataset_id
      ! Error handling.
      INTEGER :: ierr
      !
      ! Identify the HDF5 type.
      h5prec = H5KIND_TO_TYPE(REAL32, H5_REAL_KIND)
      !
      ! Create the dataspace for the data set.
      CALL H5SCREATE_SIMPLE_F(1, dims, dataspace_id, ierr)
      !
      ! Create a dataset in hdf5_id with default properties.
      CALL H5DCREATE_F(hdf5_id, dname, h5prec, dataspace_id, dataset_id, ierr)
      !
      ! Write to the dataset.
      IF (PRESENT(tid) .AND. PRESENT(wid)) THEN
         IF (tid .EQ. wid) THEN
            CALL H5DWRITE_F(dataset_id, h5prec, dat, dims, ierr)
         END IF
      ELSE
         CALL H5DWRITE_F(dataset_id, h5prec, dat, dims, ierr)
      END IF
      !
      ! Close the dataset.
      CALL H5DCLOSE_F(dataset_id, ierr)
      !
      ! Close the dataspace.
      CALL H5SCLOSE_F(dataspace_id, ierr)
   END SUBROUTINE WRITE_DATA_REAL32_RANK1

   !> Procedure to write 2D single precision real data to file.
   !!
   !> @param[in] dims Dimensions of the data.
   !> @param[in] dname Output name for the data.
   !> @param[in] dat Array with dimensions given by dims.
   !> @param[in] hdf5_id HDF5 ID for the object where data will be written.
   SUBROUTINE WRITE_DATA_REAL32_RANK2(dims, dname, dat, hdf5_id)
      ! Required modules.
      USE ISO_FORTRAN_ENV,ONLY: REAL32
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=HSIZE_T),DIMENSION(2),INTENT(IN) :: dims
      CHARACTER(LEN=*),INTENT(IN) :: dname
      REAL(KIND=REAL32),DIMENSION(dims(1),dims(2)),INTENT(IN) :: dat
      INTEGER(KIND=HID_T),INTENT(IN) :: hdf5_id
      ! Local variables.
      ! HDF5 identifier for the variable type.
      INTEGER(KIND=HID_T) :: h5prec
      ! HDF5 identifier for the dataspace and dataset.
      INTEGER(KIND=HID_T) :: dataspace_id, dataset_id
      ! Error handling.
      INTEGER :: ierr
      !
      ! Identify the HDF5 type.
      h5prec = H5KIND_TO_TYPE(REAL32, H5_REAL_KIND)
      !
      ! Create the dataspace for the data set.
      CALL H5SCREATE_SIMPLE_F(2, dims, dataspace_id, ierr)
      !
      ! Create a dataset in hdf5_id with default properties.
      CALL H5DCREATE_F(hdf5_id, dname, h5prec, dataspace_id, dataset_id, ierr)
      !
      ! Write to the dataset.
      CALL H5DWRITE_F(dataset_id, h5prec, dat, dims, ierr)
      !
      ! Close the dataset.
      CALL H5DCLOSE_F(dataset_id, ierr)
      !
      ! Close the dataspace.
      CALL H5SCLOSE_F(dataspace_id, ierr)
   END SUBROUTINE WRITE_DATA_REAL32_RANK2

   !> Procedure to write 3D single precision real data to file.
   !!
   !! If you want to write a real attribute to the dataset, e.g. the schmidt
   !! number for a scalar field, then both real_att_name and real_att_val
   !! must be passed into the routine.
   !!
   !> @param[in] dims Dimensions of the data.
   !> @param[in] dname Output name for the data.
   !> @param[in] dat Array with dimensions given by dims.
   !> @param[in] hdf5_id HDF5 ID for the object where data will be written.
   !> @param[in,optional] real_att_name Name of a real attributes.
   !> @param[in,optional] real_att_val Value for a real attributes.
   !> @param[in,optional] gbool Whether to activate gzip compression.
   !> @param[in,optional] level Level of compression.
   !> @param[in,optional] shuff Whether to add the shuffle filter. Only used in
   !! conjunction with compression filters.
   !> @param[in,optional] chunk Chunking dimensions. Can be enabled even without
   !! using compression.
   SUBROUTINE WRITE_DATA_REAL32_RANK3(dims, dname, dat, hdf5_id, &
                                      real_att_name, real_att_val, &
                                      gbool, level, shuff, chunk)
      ! Required modules.
      USE ISO_FORTRAN_ENV,ONLY: REAL32
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=HSIZE_T),DIMENSION(3),INTENT(IN) :: dims
      CHARACTER(LEN=*),INTENT(IN) :: dname
      REAL(KIND=REAL32),DIMENSION(dims(1),dims(2),dims(3)),INTENT(IN) :: dat
      INTEGER(KIND=HID_T),INTENT(IN) :: hdf5_id
      CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: real_att_name
      REAL(KIND=REAL32),OPTIONAL,INTENT(IN) :: real_att_val
      LOGICAL,OPTIONAL,INTENT(IN) :: gbool, shuff
      INTEGER,OPTIONAL,INTENT(IN) :: level
      INTEGER(KIND=HSIZE_T),DIMENSION(3),OPTIONAL,INTENT(IN) :: chunk
      ! Local variables.
      ! HDF5 identifier for the variable type.
      INTEGER(KIND=HID_T) :: h5prec
      ! HDF5 identifier for the dataspace and dataset.
      INTEGER(KIND=HID_T) :: dataspace_id, dataset_id
      ! HDF5 identifier for the dataset creation property list.
      INTEGER(KIND=HID_T) :: plist_id
      ! Used to handle optional logical arguments.
      LOGICAL :: enable_compress, shuffle_filter, enable_chunking
      INTEGER :: compress_level
      INTEGER(KIND=HSIZE_T) :: chunk_size
      INTEGER(KIND=HSIZE_T),DIMENSION(3) :: chunk_dims
      ! Error handling.
      INTEGER :: ierr
      !
      ! Handle optional arguments. Compression & shuffle off by default.
      IF (PRESENT(gbool)) THEN
         enable_compress = gbool
         !
         ! Turn off shuffle filter by default.
         IF (PRESENT(shuff)) THEN
            shuffle_filter = shuff
         ELSE
            shuffle_filter = .FALSE.
         END IF
         !
         ! Set gzip compression level of 4 by default.
         IF (PRESENT(level)) THEN
            compress_level = level
         ELSE
            compress_level = 4
         END IF
      ELSE
         enable_compress = .FALSE.
         shuffle_filter = .FALSE.
      END IF
      !
      ! Check chunking options. Turned off by default.
      IF (PRESENT(chunk)) THEN
         enable_chunking = .TRUE.
         chunk_dims(:) = chunk(:)
      ELSE IF (enable_compress) THEN
         enable_chunking = .TRUE.
         chunk_size = 32_HSIZE_T
         DO WHILE (chunk_size .GT. MINVAL(dims))
            chunk_size = chunk_size/2_HSIZE_T
         END DO
         chunk_dims(:) = chunk_size
      ELSE
         enable_chunking = .FALSE.
      END IF
      !
      ! Identify the HDF5 type.
      h5prec = H5KIND_TO_TYPE(REAL32, H5_REAL_KIND)
      !
      ! Create the dataspace for the data set.
      CALL H5SCREATE_SIMPLE_F(3, dims, dataspace_id, ierr)
      !
      ! Use property list to control dataset creation, if needed.
      IF (enable_compress .OR. enable_chunking) THEN
         !
         ! Make the property list.
         CALL H5PCREATE_F(H5P_DATASET_CREATE_F, plist_id, ierr)
         !
         ! Add the chunking option (this if statement is redundant).
         IF (enable_chunking) THEN
            CALL H5PSET_CHUNK_F(plist_id, 3, chunk_dims, ierr)
         END IF
         !
         ! Add the shuffle filter.
         IF (shuffle_filter) THEN
            CALL H5PSET_SHUFFLE_F(plist_id, ierr)
         END IF
         !
         ! Add the compression options.
         IF (enable_compress) THEN
            CALL H5PSET_DEFLATE_F(plist_id, compress_level, ierr)
         END IF
         !
         ! Make the dataset with the property list.
         CALL H5DCREATE_F(hdf5_id, dname, h5prec, dataspace_id, dataset_id, &
                          ierr, dcpl_id=plist_id)
         !
         ! Destroy the property list.
         CALL H5PCLOSE_F(plist_id, ierr)
      ELSE
         !
         ! Create a dataset in hdf5_id with default properties.
         CALL H5DCREATE_F(hdf5_id,dname,h5prec,dataspace_id,dataset_id,ierr)
      END IF
      !
      ! Write to the dataset.
      CALL H5DWRITE_F(dataset_id, h5prec, dat, dims, ierr)
      !
      ! Write a real attribute if it is present.
      IF (PRESENT(real_att_name) .AND. PRESENT(real_att_val)) THEN
         CALL WRITE_ATTRIBUTE(dataset_id, real_att_name, real_att_val)
      END IF
      !
      ! Close the dataset.
      CALL H5DCLOSE_F(dataset_id, ierr)
      !
      ! Close the dataspace.
      CALL H5SCLOSE_F(dataspace_id, ierr)
   END SUBROUTINE WRITE_DATA_REAL32_RANK3

   !> Procedure to write a subset of 3D single precision data to file.
   !!
   !> @param[in] dimsT Dimensions of the data array.
   !> @param[in] dimsW Dimensions of the subset to write out.
   !> @param[in] offset Offset of the subset in the main data array.
   !> @param[in] dname Name for the dataset.
   !> @param[in] dat Array with dimensions given by dims.
   !> @param[in] hdf5_id HDF5 ID for the object where data will be written.
   !> @param[in,optional] gbool Whether to activate gzip compression.
   !> @param[in,optional] shuff Whether to add the shuffle filter. Only used in
   !! conjunction with compression filters.
   !> @param[in,optional] level Level of compression.
   !> @param[in,optional] chunk Chunking dimensions. Can be enabled even without
   !! using compression.
   SUBROUTINE WRITE_DATA_SUBSET_REAL32_RANK3(dimsT, dimsW, offset, dname, &
                                             dat, hdf5_id, gbool, shuff, &
                                             level, chunk)
      ! Required modules.
      USE ISO_FORTRAN_ENV,ONLY: REAL32
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=HSIZE_T),DIMENSION(3),INTENT(IN) :: dimsT, dimsW, offset
      CHARACTER(LEN=*),INTENT(IN) :: dname
      REAL(KIND=REAL32),DIMENSION(dimsT(1),dimsT(2),dimsT(3)),INTENT(IN) :: dat
      INTEGER(KIND=HID_T),INTENT(IN) :: hdf5_id
      LOGICAL,OPTIONAL,INTENT(IN) :: gbool, shuff
      INTEGER,OPTIONAL,INTENT(IN) :: level
      INTEGER(KIND=HSIZE_T),DIMENSION(3),OPTIONAL,INTENT(IN) :: chunk
      ! Local variables.
      ! HDF5 identifier for the variable type.
      INTEGER(KIND=HID_T) :: h5prec
      ! HDF5 identifiers for the dataspaces and dataset.
      INTEGER(KIND=HID_T) :: filespace_id, memspace_id, dset_id
      ! Dataset creation property list.
      INTEGER(KIND=HID_T) :: plist_id
      ! Arrays required when selecting the hyperslab from the mem dataspace.
      INTEGER(KIND=HSIZE_T),DIMENSION(3) :: counter, stride
      ! Used to handle optional logical arguments.
      LOGICAL :: enable_compress, shuffle_filter, enable_chunking
      INTEGER :: compress_level
      INTEGER(KIND=HSIZE_T) :: chunk_size
      INTEGER(KIND=HSIZE_T),DIMENSION(3) :: chunk_dims
      ! Error handling.
      INTEGER :: ierr
      !
      ! Handle optional arguments. Compression & shuffle off by default.
      IF (PRESENT(gbool)) THEN
         enable_compress = gbool
         !
         ! Turn off shuffle filter by default.
         IF (PRESENT(shuff)) THEN
            shuffle_filter = shuff
         ELSE
            shuffle_filter = .FALSE.
         END IF
         !
         ! Set gzip compression level of 4 by default.
         IF (PRESENT(level)) THEN
            compress_level = level
         ELSE
            compress_level = 4
         END IF
      ELSE
         enable_compress = .FALSE.
         shuffle_filter = .FALSE.
      END IF
      !
      ! Check chunking options. Turned off by default.
      IF (PRESENT(chunk)) THEN
         enable_chunking = .TRUE.
         chunk_dims(:) = chunk(:)
      ELSE IF (enable_compress) THEN
         enable_chunking = .TRUE.
         chunk_size = 32_HSIZE_T
         DO WHILE (chunk_size .GT. MINVAL(dimsW))
            chunk_size = chunk_size/2_HSIZE_T
         END DO
         chunk_dims(:) = chunk_size
      ELSE
         enable_chunking = .FALSE.
      END IF
      !
      ! Identify the HDF5 type.
      h5prec = H5KIND_TO_TYPE(REAL32, H5_REAL_KIND)
      !
      ! Create the dataspace for the file space.
      CALL H5SCREATE_SIMPLE_F(3, dimsW, filespace_id, ierr)
      !
      ! Use property list to control dataset creation, if needed.
      IF (enable_compress .OR. enable_chunking) THEN
         !
         ! Make the property list.
         CALL H5PCREATE_F(H5P_DATASET_CREATE_F, plist_id, ierr)
         !
         ! Add the chunking option (this if statement is redundant).
         IF (enable_chunking) THEN
            CALL H5PSET_CHUNK_F(plist_id, 3, chunk_dims, ierr)
         END IF
         !
         ! Add the shuffle filter.
         IF (shuffle_filter) THEN
            CALL H5PSET_SHUFFLE_F(plist_id, ierr)
         END IF
         !
         ! Add the compression options.
         IF (enable_compress) THEN
            CALL H5PSET_DEFLATE_F(plist_id, compress_level, ierr)
         END IF
         !
         ! Make the dataset with the property list.
         CALL H5DCREATE_F(hdf5_id, dname, h5prec, filespace_id, dset_id, &
                          ierr, dcpl_id=plist_id)
         !
         ! Destroy the property list.
         CALL H5PCLOSE_F(plist_id, ierr)
      ELSE
         !
         ! Create a dataset in hdf5_id with default properties.
         CALL H5DCREATE_F(hdf5_id, dname, h5prec, filespace_id, dset_id, ierr)
      END IF
      !
      ! Create the dataspace for the memory on this process.
      CALL H5SCREATE_SIMPLE_F(3, dimsT, memspace_id, ierr)
      !
      ! Select the hyperslab in the memory dataspace to write out.
      counter = [1_HSIZE_T, 1_HSIZE_T, 1_HSIZE_T]
      stride = [1_HSIZE_T, 1_HSIZE_T, 1_HSIZE_T]
      CALL H5SSELECT_HYPERSLAB_F(memspace_id, H5S_SELECT_SET_F, offset, &
                                 counter, ierr, stride, dimsW)
      !
      ! Write out the data to file.
      CALL H5DWRITE_F(dset_id, h5prec, dat, dimsT, ierr, &
                      file_space_id=filespace_id, mem_space_id=memspace_id)
      !
      ! Close the dataset.
      CALL H5DCLOSE_F(dset_id, ierr)
      !
      ! Close the dataspaces.
      CALL H5SCLOSE_F(filespace_id, ierr)
      CALL H5SCLOSE_F(memspace_id, ierr)
   END SUBROUTINE WRITE_DATA_SUBSET_REAL32_RANK3

   !> Procedure to write 1D double precision real data to file.
   !!
   !! This routine is designed to work in the case that parallel HDF5 is being
   !! used, but the user only wants one process to write to the dataset. Since
   !! dataset creation is collective, all processes participate in making the
   !! dataset, but only the desired process writes to it (assuming the optional
   !! arguments tid and wid are present).
   !!
   !> @param[in] dims Dimensions of the data.
   !> @param[in] dname Output name for the data.
   !> @param[in] dat Array with dimensions given by dims.
   !> @param[in] hdf5_id HDF5 ID for the object where data will be written.
   !> @param[in,optional] tid Process identifier (rank).
   !> @param[in,optional] wid Process ID to write out the data.
   SUBROUTINE WRITE_DATA_REAL64_RANK1(dims, dname, dat, hdf5_id, tid, wid)
      ! Required modules.
      USE ISO_FORTRAN_ENV,ONLY: REAL64
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=HSIZE_T),DIMENSION(1),INTENT(IN) :: dims
      CHARACTER(LEN=*),INTENT(IN) :: dname
      REAL(KIND=REAL64),DIMENSION(dims(1)),INTENT(IN) :: dat
      INTEGER(KIND=HID_T),INTENT(IN) :: hdf5_id
      INTEGER,OPTIONAL,INTENT(IN) :: tid, wid
      ! Local variables.
      ! HDF5 identifier for the variable type.
      INTEGER(KIND=HID_T) :: h5prec
      ! HDF5 identifier for the dataspace and dataset.
      INTEGER(KIND=HID_T) :: dataspace_id, dataset_id
      ! Error handling.
      INTEGER :: ierr
      !
      ! Identify the HDF5 type.
      h5prec = H5KIND_TO_TYPE(REAL64, H5_REAL_KIND)
      !
      ! Create the dataspace for the data set.
      CALL H5SCREATE_SIMPLE_F(1, dims, dataspace_id, ierr)
      !
      ! Create a dataset in hdf5_id with default properties.
      CALL H5DCREATE_F(hdf5_id, dname, h5prec, dataspace_id, dataset_id, ierr)
      !
      ! Write to the dataset.
      IF (PRESENT(tid) .AND. PRESENT(wid)) THEN
         IF (tid .EQ. wid) THEN
            CALL H5DWRITE_F(dataset_id, h5prec, dat, dims, ierr)
         END IF
      ELSE
         CALL H5DWRITE_F(dataset_id, h5prec, dat, dims, ierr)
      END IF
      !
      ! Close the dataset.
      CALL H5DCLOSE_F(dataset_id, ierr)
      !
      ! Close the dataspace.
      CALL H5SCLOSE_F(dataspace_id, ierr)
   END SUBROUTINE WRITE_DATA_REAL64_RANK1

   !> Procedure to write 2D double precision real data to file.
   !!
   !> @param[in] dims Dimensions of the data.
   !> @param[in] dname Output name for the data.
   !> @param[in] dat Array with dimensions given by dims.
   !> @param[in] hdf5_id HDF5 ID for the object where data will be written.
   SUBROUTINE WRITE_DATA_REAL64_RANK2(dims, dname, dat, hdf5_id)
      ! Required modules.
      USE ISO_FORTRAN_ENV,ONLY: REAL64
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=HSIZE_T),DIMENSION(2),INTENT(IN) :: dims
      CHARACTER(LEN=*),INTENT(IN) :: dname
      REAL(KIND=REAL64),DIMENSION(dims(1),dims(2)),INTENT(IN) :: dat
      INTEGER(KIND=HID_T),INTENT(IN) :: hdf5_id
      ! Local variables.
      ! HDF5 identifier for the variable type.
      INTEGER(KIND=HID_T) :: h5prec
      ! HDF5 identifier for the dataspace and dataset.
      INTEGER(KIND=HID_T) :: dataspace_id, dataset_id
      ! Error handling.
      INTEGER :: ierr
      !
      ! Identify the HDF5 type.
      h5prec = H5KIND_TO_TYPE(REAL64, H5_REAL_KIND)
      !
      ! Create the dataspace for the data set.
      CALL H5SCREATE_SIMPLE_F(2, dims, dataspace_id, ierr)
      !
      ! Create a dataset in hdf5_id with default properties.
      CALL H5DCREATE_F(hdf5_id, dname, h5prec, dataspace_id, dataset_id, ierr)
      !
      ! Write to the dataset.
      CALL H5DWRITE_F(dataset_id, h5prec, dat, dims, ierr)
      !
      ! Close the dataset.
      CALL H5DCLOSE_F(dataset_id, ierr)
      !
      ! Close the dataspace.
      CALL H5SCLOSE_F(dataspace_id, ierr)
   END SUBROUTINE WRITE_DATA_REAL64_RANK2

   !> Procedure to write 3D double precision real data to file.
   !!
   !! If you want to write a real attribute to the dataset, e.g. the schmidt
   !! number for a scalar field, then both real_att_name and real_att_val
   !! must be passed into the routine.
   !!
   !> @param[in] dims Dimensions of the data.
   !> @param[in] dname Output name for the data.
   !> @param[in] dat Array with dimensions given by dims.
   !> @param[in] hdf5_id HDF5 ID for the object where data will be written.
   !> @param[in,optional] real_att_name Name of a real attributes.
   !> @param[in,optional] real_att_val Value for a real attributes.
   !> @param[in,optional] gbool Whether to activate gzip compression.
   !> @param[in,optional] level Level of compression.
   !> @param[in,optional] shuff Whether to add the shuffle filter. Only used in
   !! conjunction with compression filters.
   !> @param[in,optional] chunk Chunking dimensions. Can be enabled even without
   !! using compression.
   SUBROUTINE WRITE_DATA_REAL64_RANK3(dims, dname, dat, hdf5_id, &
                                      real_att_name, real_att_val, &
                                      gbool, level, shuff, chunk)
      ! Required modules.
      USE ISO_FORTRAN_ENV,ONLY: REAL64
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=HSIZE_T),DIMENSION(3),INTENT(IN) :: dims
      CHARACTER(LEN=*),INTENT(IN) :: dname
      REAL(KIND=REAL64),DIMENSION(dims(1),dims(2),dims(3)),INTENT(IN) :: dat
      INTEGER(KIND=HID_T),INTENT(IN) :: hdf5_id
      CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: real_att_name
      REAL(KIND=REAL64),OPTIONAL,INTENT(IN) :: real_att_val
      LOGICAL,OPTIONAL,INTENT(IN) :: gbool, shuff
      INTEGER,OPTIONAL,INTENT(IN) :: level
      INTEGER(KIND=HSIZE_T),DIMENSION(3),OPTIONAL,INTENT(IN) :: chunk
      ! Local variables.
      ! HDF5 identifier for the variable type.
      INTEGER(KIND=HID_T) :: h5prec
      ! HDF5 identifier for the dataspace and dataset.
      INTEGER(KIND=HID_T) :: dataspace_id, dataset_id
      ! HDF5 identifier for the dataset creation property list.
      INTEGER(KIND=HID_T) :: plist_id
      ! Used to handle optional logical arguments.
      LOGICAL :: enable_compress, shuffle_filter, enable_chunking
      INTEGER :: compress_level
      INTEGER(KIND=HSIZE_T) :: chunk_size
      INTEGER(KIND=HSIZE_T),DIMENSION(3) :: chunk_dims
      ! Error handling.
      INTEGER :: ierr
      !
      ! Handle optional arguments. Compression & shuffle off by default.
      IF (PRESENT(gbool)) THEN
         enable_compress = gbool
         !
         ! Turn off shuffle filter by default.
         IF (PRESENT(shuff)) THEN
            shuffle_filter = shuff
         ELSE
            shuffle_filter = .FALSE.
         END IF
         !
         ! Set gzip compression level of 4 by default.
         IF (PRESENT(level)) THEN
            compress_level = level
         ELSE
            compress_level = 4
         END IF
      ELSE
         enable_compress = .FALSE.
         shuffle_filter = .FALSE.
      END IF
      !
      ! Check chunking options. Turned off by default.
      IF (PRESENT(chunk)) THEN
         enable_chunking = .TRUE.
         chunk_dims(:) = chunk(:)
      ELSE IF (enable_compress) THEN
         enable_chunking = .TRUE.
         chunk_size = 32_HSIZE_T
         DO WHILE (chunk_size .GT. MINVAL(dims))
            chunk_size = chunk_size/2_HSIZE_T
         END DO
         chunk_dims(:) = chunk_size
      ELSE
         enable_chunking = .FALSE.
      END IF
      !
      ! Identify the HDF5 type.
      h5prec = H5KIND_TO_TYPE(REAL64, H5_REAL_KIND)
      !
      ! Create the dataspace for the data set.
      CALL H5SCREATE_SIMPLE_F(3, dims, dataspace_id, ierr)
      !
      ! Use property list to control dataset creation, if needed.
      IF (enable_compress .OR. enable_chunking) THEN
         !
         ! Make the property list.
         CALL H5PCREATE_F(H5P_DATASET_CREATE_F, plist_id, ierr)
         !
         ! Add the chunking option (this if statement is redundant).
         IF (enable_chunking) THEN
            CALL H5PSET_CHUNK_F(plist_id, 3, chunk_dims, ierr)
         END IF
         !
         ! Add the shuffle filter.
         IF (shuffle_filter) THEN
            CALL H5PSET_SHUFFLE_F(plist_id, ierr)
         END IF
         !
         ! Add the compression options.
         IF (enable_compress) THEN
            CALL H5PSET_DEFLATE_F(plist_id, compress_level, ierr)
         END IF
         !
         ! Make the dataset with the property list.
         CALL H5DCREATE_F(hdf5_id, dname, h5prec, dataspace_id, dataset_id, &
                          ierr, dcpl_id=plist_id)
         !
         ! Destroy the property list.
         CALL H5PCLOSE_F(plist_id, ierr)
      ELSE
         !
         ! Create a dataset in hdf5_id with default properties.
         CALL H5DCREATE_F(hdf5_id,dname,h5prec,dataspace_id,dataset_id,ierr)
      END IF
      !
      ! Write to the dataset.
      CALL H5DWRITE_F(dataset_id, h5prec, dat, dims, ierr)
      !
      ! Write a real attribute if it is present.
      IF (PRESENT(real_att_name) .AND. PRESENT(real_att_val)) THEN
         CALL WRITE_ATTRIBUTE(dataset_id, real_att_name, real_att_val)
      END IF
      !
      ! Close the dataset.
      CALL H5DCLOSE_F(dataset_id, ierr)
      !
      ! Close the dataspace.
      CALL H5SCLOSE_F(dataspace_id, ierr)
   END SUBROUTINE WRITE_DATA_REAL64_RANK3

   !> Procedure to write a subset of 3D double precision data to file.
   !!
   !> @param[in] dimsT Dimensions of the data array.
   !> @param[in] dimsW Dimensions of the subset to write out.
   !> @param[in] offset Offset of the subset in the main data array.
   !> @param[in] dname Name for the dataset.
   !> @param[in] dat Array with dimensions given by dims.
   !> @param[in] hdf5_id HDF5 ID for the object where data will be written.
   !> @param[in,optional] gbool Whether to activate gzip compression.
   !> @param[in,optional] shuff Whether to add the shuffle filter. Only used in
   !! conjunction with compression filters.
   !> @param[in,optional] level Level of compression.
   !> @param[in,optional] chunk Chunking dimensions. Can be enabled even without
   !! using compression.
   SUBROUTINE WRITE_DATA_SUBSET_REAL64_RANK3(dimsT, dimsW, offset, dname, &
                                             dat, hdf5_id, gbool, shuff, &
                                             level, chunk)
      ! Required modules.
      USE ISO_FORTRAN_ENV,ONLY: REAL64
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=HSIZE_T),DIMENSION(3),INTENT(IN) :: dimsT, dimsW, offset
      CHARACTER(LEN=*),INTENT(IN) :: dname
      REAL(KIND=REAL64),DIMENSION(dimsT(1),dimsT(2),dimsT(3)),INTENT(IN) :: dat
      INTEGER(KIND=HID_T),INTENT(IN) :: hdf5_id
      LOGICAL,OPTIONAL,INTENT(IN) :: gbool, shuff
      INTEGER,OPTIONAL,INTENT(IN) :: level
      INTEGER(KIND=HSIZE_T),DIMENSION(3),OPTIONAL,INTENT(IN) :: chunk
      ! Local variables.
      ! HDF5 identifier for the variable type.
      INTEGER(KIND=HID_T) :: h5prec
      ! HDF5 identifiers for the dataspaces and dataset.
      INTEGER(KIND=HID_T) :: filespace_id, memspace_id, dset_id
      ! Dataset creation property list.
      INTEGER(KIND=HID_T) :: plist_id
      ! Arrays required when selecting the hyperslab from the mem dataspace.
      INTEGER(KIND=HSIZE_T),DIMENSION(3) :: counter, stride
      ! Used to handle optional logical arguments.
      LOGICAL :: enable_compress, shuffle_filter, enable_chunking
      INTEGER :: compress_level
      INTEGER(KIND=HSIZE_T) :: chunk_size
      INTEGER(KIND=HSIZE_T),DIMENSION(3) :: chunk_dims
      ! Error handling.
      INTEGER :: ierr
      !
      ! Handle optional arguments. Compression & shuffle off by default.
      IF (PRESENT(gbool)) THEN
         enable_compress = gbool
         !
         ! Turn off shuffle filter by default.
         IF (PRESENT(shuff)) THEN
            shuffle_filter = shuff
         ELSE
            shuffle_filter = .FALSE.
         END IF
         !
         ! Set gzip compression level of 4 by default.
         IF (PRESENT(level)) THEN
            compress_level = level
         ELSE
            compress_level = 4
         END IF
      ELSE
         enable_compress = .FALSE.
         shuffle_filter = .FALSE.
      END IF
      !
      ! Check chunking options. Turned off by default.
      IF (PRESENT(chunk)) THEN
         enable_chunking = .TRUE.
         chunk_dims(:) = chunk(:)
      ELSE IF (enable_compress) THEN
         enable_chunking = .TRUE.
         chunk_size = 32_HSIZE_T
         DO WHILE (chunk_size .GT. MINVAL(dimsW))
            chunk_size = chunk_size/2_HSIZE_T
         END DO
         chunk_dims(:) = chunk_size
      ELSE
         enable_chunking = .FALSE.
      END IF
      !
      ! Identify the HDF5 type.
      h5prec = H5KIND_TO_TYPE(REAL64, H5_REAL_KIND)
      !
      ! Create the dataspace for the file space.
      CALL H5SCREATE_SIMPLE_F(3, dimsW, filespace_id, ierr)
      !
      ! Use property list to control dataset creation, if needed.
      IF (enable_compress .OR. enable_chunking) THEN
         !
         ! Make the property list.
         CALL H5PCREATE_F(H5P_DATASET_CREATE_F, plist_id, ierr)
         !
         ! Add the chunking option (this if statement is redundant).
         IF (enable_chunking) THEN
            CALL H5PSET_CHUNK_F(plist_id, 3, chunk_dims, ierr)
         END IF
         !
         ! Add the shuffle filter.
         IF (shuffle_filter) THEN
            CALL H5PSET_SHUFFLE_F(plist_id, ierr)
         END IF
         !
         ! Add the compression options.
         IF (enable_compress) THEN
            CALL H5PSET_DEFLATE_F(plist_id, compress_level, ierr)
         END IF
         !
         ! Make the dataset with the property list.
         CALL H5DCREATE_F(hdf5_id, dname, h5prec, filespace_id, dset_id, &
                          ierr, dcpl_id=plist_id)
         !
         ! Destroy the property list.
         CALL H5PCLOSE_F(plist_id, ierr)
      ELSE
         !
         ! Create a dataset in hdf5_id with default properties.
         CALL H5DCREATE_F(hdf5_id, dname, h5prec, filespace_id, dset_id, ierr)
      END IF
      !
      ! Create the dataspace for the memory on this process.
      CALL H5SCREATE_SIMPLE_F(3, dimsT, memspace_id, ierr)
      !
      ! Select the hyperslab in the memory dataspace to write out.
      counter = [1_HSIZE_T, 1_HSIZE_T, 1_HSIZE_T]
      stride = [1_HSIZE_T, 1_HSIZE_T, 1_HSIZE_T]
      CALL H5SSELECT_HYPERSLAB_F(memspace_id, H5S_SELECT_SET_F, offset, &
                                 counter, ierr, stride, dimsW)
      !
      ! Write out the data to file.
      CALL H5DWRITE_F(dset_id, h5prec, dat, dimsT, ierr, &
                      file_space_id=filespace_id, mem_space_id=memspace_id)
      !
      ! Close the dataset.
      CALL H5DCLOSE_F(dset_id, ierr)
      !
      ! Close the dataspaces.
      CALL H5SCLOSE_F(filespace_id, ierr)
      CALL H5SCLOSE_F(memspace_id, ierr)
   END SUBROUTINE WRITE_DATA_SUBSET_REAL64_RANK3

   !> Procedure to write 1D integer data to file in parallel.
   !!
   !! Note that if the dimsW argument is passed, in the first 1:dimsW elements
   !! of the data array are written to file. The user must ensure the offset,
   !! etc. are properly set for each process when writing either the entire
   !! dimsP or reduced dimsW amount of data.
   !!
   !> @param[in] hdf5_id HDF5 ID for the object where data will be written.
   !> @param[in] dname Output name for the data.
   !> @param[in] dimsT Dimensions of the complete dataset across all procs.
   !> @param[in] dimsP Size of the data array on this process.
   !> @param[in] offset Offset for this process in the joined checkpoint.
   !> @param[in] counter Number of dimsP blocks for this process.
   !> @param[in] stride Strides to take in the HDF5 data in each direction.
   !> @param[in] dat Array with dimensions given by dimsP.
   !> @param[in,optional] chunk Chunking dimension, if desired.
   !> @param[in,optional] dimsW Section of data to write for this process.
   SUBROUTINE WRITE_DATA_PARALLEL_INT32_RANK1(hdf5_id, dname, dimsT, dimsP, &
                                              offset, counter, stride, dat, &
                                              chunk, dimsW)
      ! Required modules.
      USE ISO_FORTRAN_ENV,ONLY: INT32
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=HID_T),INTENT(IN) :: hdf5_id
      CHARACTER(LEN=*),INTENT(IN) :: dname
      INTEGER(KIND=HSIZE_T),DIMENSION(1),INTENT(IN) :: dimsT, dimsP, offset, &
                                                       counter, stride
      INTEGER(KIND=INT32),DIMENSION(dimsP(1)),INTENT(IN) :: dat
      INTEGER(KIND=HSIZE_T),DIMENSION(1),OPTIONAL,INTENT(IN) :: chunk, dimsW
      ! Local variables.
      ! HDF5 identifier for the variable type.
      INTEGER(KIND=HID_T) :: h5prec
      ! HDF5 dataspace identifier for the file space.
      INTEGER(KIND=HID_T) :: filespace
      ! HDF5 dataspace identifier for the memory on this process to be written.
      INTEGER(KIND=HID_T) :: memspace
      ! HDF5 property list identifier.
      INTEGER(KIND=HID_T) :: plist_id
      ! HDF5 dataset identifier.
      INTEGER(KIND=HID_T) :: dataset_id
      ! Local arrays to trim the process memory space if dimsW is present.
      INTEGER(KIND=HSIZE_T),DIMENSION(1) :: oset_tmp
      ! Error handling.
      INTEGER :: ierr
      !
      ! Identify the HDF5 type.
      h5prec = H5KIND_TO_TYPE(INT32, H5_INTEGER_KIND)
      !
      ! Create the dataspace for the file space.
      CALL H5SCREATE_SIMPLE_F(1, dimsT, filespace, ierr)
      !
      ! Create the dataspace for the memory on this process.
      CALL H5SCREATE_SIMPLE_F(1, dimsP, memspace, ierr)
      !
      ! Create the dataset with or without chunking, depending on whether or not
      ! the optional "chunk" argument was passed in.
      IF (PRESENT(chunk)) THEN
         CALL H5PCREATE_F(H5P_DATASET_CREATE_F, plist_id, ierr)
         CALL H5PSET_CHUNK_F(plist_id, 1, chunk, ierr)
         CALL H5DCREATE_F(hdf5_id, dname, h5prec, filespace, dataset_id, &
                          ierr, plist_id)
         CALL H5PCLOSE_F(plist_id, ierr)
      ELSE
         CALL H5DCREATE_F(hdf5_id, dname, h5prec, filespace, dataset_id, ierr)
      END IF
      !
      ! Close the dataspace for the file space.
      CALL H5SCLOSE_F(filespace, ierr)
      !
      ! Select the hyperslab in the file. If the dimsW optional argument is
      ! present, we write out the first dimsW elements of the memory on this
      ! process. If it is not present, we write out the entire dimsP memory.
      IF (PRESENT(dimsW)) THEN
         CALL H5DGET_SPACE_F(dataset_id, filespace, ierr)
         CALL H5SSELECT_HYPERSLAB_F(filespace, H5S_SELECT_SET_F, offset, &
                                    counter, ierr, stride, dimsW)
         !
         ! If dimsW is present, we also need to trim the memory space for the
         ! local process to only write out the first dimsW elements.
         oset_tmp(1) = 0_HSIZE_T
         CALL H5SSELECT_HYPERSLAB_F(memspace, H5S_SELECT_SET_F, oset_tmp, &
                                    dimsW, ierr)
      ELSE
         CALL H5DGET_SPACE_F(dataset_id, filespace, ierr)
         CALL H5SSELECT_HYPERSLAB_F(filespace, H5S_SELECT_SET_F, offset, &
                                    counter, ierr, stride, dimsP)
      END IF
      !
      ! Create the property list for collective dataset writing.
      CALL H5PCREATE_F(H5P_DATASET_XFER_F, plist_id, ierr)
      CALL H5PSET_DXPL_MPIO_F(plist_id, H5FD_MPIO_COLLECTIVE_F, ierr)
      !
      ! Write the dataset collectively.
      CALL H5DWRITE_F(dataset_id, h5prec, dat, dimsP, ierr, &
                      file_space_id=filespace, mem_space_id=memspace, &
                      xfer_prp=plist_id)
      !
      ! Close the dataset.
      CALL H5DCLOSE_F(dataset_id, ierr)
      !
      ! Close the property list.
      CALL H5PCLOSE_F(plist_id, ierr)
      !
      ! Close the dataspaces.
      CALL H5SCLOSE_F(filespace, ierr)
      CALL H5SCLOSE_F(memspace, ierr)
   END SUBROUTINE WRITE_DATA_PARALLEL_INT32_RANK1

   !> Procedure to write 1D integer data to file in parallel.
   !!
   !! Note that if the dimsW argument is passed, in the first 1:dimsW elements
   !! of the data array are written to file. The user must ensure the offset,
   !! etc. are properly set for each process when writing either the entire
   !! dimsP or reduced dimsW amount of data.
   !!
   !> @param[in] hdf5_id HDF5 ID for the object where data will be written.
   !> @param[in] dname Output name for the data.
   !> @param[in] dimsT Dimensions of the complete dataset across all procs.
   !> @param[in] dimsP Size of the data array on this process.
   !> @param[in] offset Offset for this process in the joined checkpoint.
   !> @param[in] counter Number of dimsP blocks for this process.
   !> @param[in] stride Strides to take in the HDF5 data in each direction.
   !> @param[in] dat Array with dimensions given by dimsP.
   !> @param[in,optional] chunk Chunking dimension, if desired.
   !> @param[in,optional] dimsW Section of data to write for this process.
#ifdef NEW_HDF5
   SUBROUTINE WRITE_DATA_PARALLEL_INT64_RANK1(hdf5_id, dname, dimsT, dimsP, &
                                              offset, counter, stride, dat, &
                                              chunk, dimsW)
      ! Required modules.
      USE ISO_FORTRAN_ENV,ONLY: INT64
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=HID_T),INTENT(IN) :: hdf5_id
      CHARACTER(LEN=*),INTENT(IN) :: dname
      INTEGER(KIND=HSIZE_T),DIMENSION(1),INTENT(IN) :: dimsT, dimsP, offset, &
                                                       counter, stride
      INTEGER(KIND=INT64),DIMENSION(dimsP(1)),INTENT(IN) :: dat
      INTEGER(KIND=HSIZE_T),DIMENSION(1),OPTIONAL,INTENT(IN) :: chunk, dimsW
      ! Local variables.
      ! HDF5 identifier for the variable type.
      INTEGER(KIND=HID_T) :: h5prec
      ! HDF5 dataspace identifier for the file space.
      INTEGER(KIND=HID_T) :: filespace
      ! HDF5 dataspace identifier for the memory on this process to be written.
      INTEGER(KIND=HID_T) :: memspace
      ! HDF5 property list identifier.
      INTEGER(KIND=HID_T) :: plist_id
      ! HDF5 dataset identifier.
      INTEGER(KIND=HID_T) :: dataset_id
      ! Local arrays to trim the process memory space if dimsW is present.
      INTEGER(KIND=HSIZE_T),DIMENSION(1) :: oset_tmp
      ! Error handling.
      INTEGER :: ierr
      !
      ! Identify the HDF5 type.
      h5prec = H5KIND_TO_TYPE(INT64, H5_INTEGER_KIND)
      !
      ! Create the dataspace for the file space.
      CALL H5SCREATE_SIMPLE_F(1, dimsT, filespace, ierr)
      !
      ! Create the dataspace for the memory on this process.
      CALL H5SCREATE_SIMPLE_F(1, dimsP, memspace, ierr)
      !
      ! Create the dataset with or without chunking, depending on whether or not
      ! the optional "chunk" argument was passed in.
      IF (PRESENT(chunk)) THEN
         CALL H5PCREATE_F(H5P_DATASET_CREATE_F, plist_id, ierr)
         CALL H5PSET_CHUNK_F(plist_id, 1, chunk, ierr)
         CALL H5DCREATE_F(hdf5_id, dname, h5prec, filespace, dataset_id, &
                          ierr, plist_id)
         CALL H5PCLOSE_F(plist_id, ierr)
      ELSE
         CALL H5DCREATE_F(hdf5_id, dname, h5prec, filespace, dataset_id, ierr)
      END IF
      !
      ! Close the dataspace for the file space.
      CALL H5SCLOSE_F(filespace, ierr)
      !
      ! Select the hyperslab in the file. If the dimsW optional argument is
      ! present, we write out the first dimsW elements of the memory on this
      ! process. If it is not present, we write out the entire dimsP memory.
      IF (PRESENT(dimsW)) THEN
         CALL H5DGET_SPACE_F(dataset_id, filespace, ierr)
         CALL H5SSELECT_HYPERSLAB_F(filespace, H5S_SELECT_SET_F, offset, &
                                    counter, ierr, stride, dimsW)
         !
         ! If dimsW is present, we also need to trim the memory space for the
         ! local process to only write out the first dimsW elements.
         oset_tmp(1) = 0_HSIZE_T
         CALL H5SSELECT_HYPERSLAB_F(memspace, H5S_SELECT_SET_F, oset_tmp, &
                                    dimsW, ierr)
      ELSE
         CALL H5DGET_SPACE_F(dataset_id, filespace, ierr)
         CALL H5SSELECT_HYPERSLAB_F(filespace, H5S_SELECT_SET_F, offset, &
                                    counter, ierr, stride, dimsP)
      END IF
      !
      ! Create the property list for collective dataset writing.
      CALL H5PCREATE_F(H5P_DATASET_XFER_F, plist_id, ierr)
      CALL H5PSET_DXPL_MPIO_F(plist_id, H5FD_MPIO_COLLECTIVE_F, ierr)
      !
      ! Write the dataset collectively.
      CALL H5DWRITE_F(dataset_id, h5prec, dat, dimsP, ierr, &
                      file_space_id=filespace, mem_space_id=memspace, &
                      xfer_prp=plist_id)
      !
      ! Close the dataset.
      CALL H5DCLOSE_F(dataset_id, ierr)
      !
      ! Close the property list.
      CALL H5PCLOSE_F(plist_id, ierr)
      !
      ! Close the dataspaces.
      CALL H5SCLOSE_F(filespace, ierr)
      CALL H5SCLOSE_F(memspace, ierr)
   END SUBROUTINE WRITE_DATA_PARALLEL_INT64_RANK1
#endif

   !> Procedure to write 1D single precision real data to file in parallel.
   !!
   !! Note that if the dimsW argument is passed, in the first 1:dimsW elements
   !! of the data array are written to file. The user must ensure the offset,
   !! etc. are properly set for each process when writing either the entire
   !! dimsP or reduced dimsW amount of data.
   !!
   !> @param[in] hdf5_id HDF5 ID for the object where data will be written.
   !> @param[in] dname Output name for the data.
   !> @param[in] dimsT Dimensions of the complete dataset across all procs.
   !> @param[in] dimsP Size of the data array on this process.
   !> @param[in] offset Offset for this process in the joined checkpoint.
   !> @param[in] counter Number of dimsP blocks for this process.
   !> @param[in] stride Strides to take in the HDF5 data in each direction.
   !> @param[in] dat Array with dimensions given by dimsP.
   !> @param[in,optional] chunk Chunking dimension, if desired.
   !> @param[in,optional] dimsW Section of data to write for this process.
   SUBROUTINE WRITE_DATA_PARALLEL_REAL32_RANK1(hdf5_id, dname, dimsT, dimsP, &
                                               offset, counter, stride, dat, &
                                               chunk, dimsW)
      ! Required modules.
      USE ISO_FORTRAN_ENV,ONLY: REAL32
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=HID_T),INTENT(IN) :: hdf5_id
      CHARACTER(LEN=*),INTENT(IN) :: dname
      INTEGER(KIND=HSIZE_T),DIMENSION(1),INTENT(IN) :: dimsT, dimsP, offset, &
                                                       counter, stride
      REAL(KIND=REAL32),DIMENSION(dimsP(1)),INTENT(IN) :: dat
      INTEGER(KIND=HSIZE_T),DIMENSION(1),OPTIONAL,INTENT(IN) :: chunk, dimsW
      ! Local variables.
      ! HDF5 identifier for the variable type.
      INTEGER(KIND=HID_T) :: h5prec
      ! HDF5 dataspace identifier for the file space.
      INTEGER(KIND=HID_T) :: filespace
      ! HDF5 dataspace identifier for the memory on this process to be written.
      INTEGER(KIND=HID_T) :: memspace
      ! HDF5 property list identifier.
      INTEGER(KIND=HID_T) :: plist_id
      ! HDF5 dataset identifier.
      INTEGER(KIND=HID_T) :: dataset_id
      ! Local arrays to trim the process memory space if dimsW is present.
      INTEGER(KIND=HSIZE_T),DIMENSION(1) :: oset_tmp
      ! Error handling.
      INTEGER :: ierr
      !
      ! Identify the HDF5 type.
      h5prec = H5KIND_TO_TYPE(REAL32, H5_REAL_KIND)
      !
      ! Create the dataspace for the file space.
      CALL H5SCREATE_SIMPLE_F(1, dimsT, filespace, ierr)
      !
      ! Create the dataspace for the memory on this process.
      CALL H5SCREATE_SIMPLE_F(1, dimsP, memspace, ierr)
      !
      ! Create the dataset with or without chunking, depending on whether or not
      ! the optional "chunk" argument was passed in.
      IF (PRESENT(chunk)) THEN
         CALL H5PCREATE_F(H5P_DATASET_CREATE_F, plist_id, ierr)
         CALL H5PSET_CHUNK_F(plist_id, 1, chunk, ierr)
         CALL H5DCREATE_F(hdf5_id, dname, h5prec, filespace, dataset_id, &
                          ierr, plist_id)
         CALL H5PCLOSE_F(plist_id, ierr)
      ELSE
         CALL H5DCREATE_F(hdf5_id, dname, h5prec, filespace, dataset_id, ierr)
      END IF
      !
      ! Close the dataspace for the file space.
      CALL H5SCLOSE_F(filespace, ierr)
      !
      ! Select the hyperslab in the file. If the dimsW optional argument is
      ! present, we write out the first dimsW elements of the memory on this
      ! process. If it is not present, we write out the entire dimsP memory.
      IF (PRESENT(dimsW)) THEN
         CALL H5DGET_SPACE_F(dataset_id, filespace, ierr)
         CALL H5SSELECT_HYPERSLAB_F(filespace, H5S_SELECT_SET_F, offset, &
                                    counter, ierr, stride, dimsW)
         !
         ! If dimsW is present, we also need to trim the memory space for the
         ! local process to only write out the first dimsW elements.
         oset_tmp(1) = 0_HSIZE_T
         CALL H5SSELECT_HYPERSLAB_F(memspace, H5S_SELECT_SET_F, oset_tmp, &
                                    dimsW, ierr)
      ELSE
         CALL H5DGET_SPACE_F(dataset_id, filespace, ierr)
         CALL H5SSELECT_HYPERSLAB_F(filespace, H5S_SELECT_SET_F, offset, &
                                    counter, ierr, stride, dimsP)
      END IF
      !
      ! Create the property list for collective dataset writing.
      CALL H5PCREATE_F(H5P_DATASET_XFER_F, plist_id, ierr)
      CALL H5PSET_DXPL_MPIO_F(plist_id, H5FD_MPIO_COLLECTIVE_F, ierr)
      !
      ! Write the dataset collectively.
      CALL H5DWRITE_F(dataset_id, h5prec, dat, dimsP, ierr, &
                      file_space_id=filespace, mem_space_id=memspace, &
                      xfer_prp=plist_id)
      !
      ! Close the dataset.
      CALL H5DCLOSE_F(dataset_id, ierr)
      !
      ! Close the property list.
      CALL H5PCLOSE_F(plist_id, ierr)
      !
      ! Close the dataspaces.
      CALL H5SCLOSE_F(filespace, ierr)
      CALL H5SCLOSE_F(memspace, ierr)
   END SUBROUTINE WRITE_DATA_PARALLEL_REAL32_RANK1

   !> Procedure to write 3D single precision real data to file in parallel.
   !!
   !! If you want to write a real attribute to the dataset, e.g. the schmidt
   !! number for a scalar field, then both real_att_name and real_att_val
   !! must be passed into the routine.
   !!
   !> @param[in] hdf5_id HDF5 ID for the object where data will be written.
   !> @param[in] dname Output name for the data.
   !> @param[in] dimT_proc Total array size for data on this processor.
   !> @param[in] dimW_proc Dimensions of data on this process to write out.
   !> @param[in] oset_proc Offset for the data to write from this process.
   !> @param[in] dimT_file Total size of the data to write out.
   !> @param[in] oset_file Offset for the data to be written from this process.
   !> @param[in] dat Array with dimensions given by dimT_proc.
   !> @param[in,optional] real_att_name Name of a real attributes.
   !> @param[in,optional] real_att_val Value for a real attributes.
   SUBROUTINE WRITE_DATA_PARALLEL_REAL32_RANK3(hdf5_id, dname, dimT_proc, &
                                               dimW_proc, oset_proc, &
                                               dimT_file, oset_file, dat, &
                                               real_att_name, real_att_val)
      ! Required modules.
      USE ISO_FORTRAN_ENV,ONLY: REAL32
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=HID_T),INTENT(IN) :: hdf5_id
      CHARACTER(LEN=*),INTENT(IN) :: dname
      INTEGER(KIND=HSIZE_T),DIMENSION(3),INTENT(IN) :: dimT_proc, dimW_proc, &
                                                       oset_proc, dimT_file, &
                                                       oset_file
      REAL(KIND=REAL32),DIMENSION(dimT_proc(1),dimT_proc(2),dimT_proc(3)), &
         INTENT(IN) :: dat
      CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: real_att_name
      REAL(KIND=REAL32),OPTIONAL,INTENT(IN) :: real_att_val
      ! Local variables.
      ! HDF5 identifier for the variable type.
      INTEGER(KIND=HID_T) :: h5prec
      ! HDF5 dataspace identifier for the file space.
      INTEGER(KIND=HID_T) :: fle_space
      ! HDF5 dataspace identifier for the memory on this process to be written.
      INTEGER(KIND=HID_T) :: mem_space
      ! HDF5 property list identifier.
      INTEGER(KIND=HID_T) :: plist_id
      ! HDF5 dataset identifier.
      INTEGER(KIND=HID_T) :: dset_id
      ! Error handling.
      INTEGER :: ierr

      ! Identify the HDF5 type to write.
      h5prec = H5KIND_TO_TYPE(REAL32,H5_REAL_KIND)

      ! Create the dataspace for the memory on this process.
      CALL H5SCREATE_SIMPLE_F(3,dimT_proc,mem_space,ierr)
      !
      ! Select the region of the memory space to write out.
      CALL H5SSELECT_HYPERSLAB_F(mem_space,H5S_SELECT_SET_F,oset_proc, &
                                 dimW_proc,ierr)

      ! Create the dataspace for the file.
      CALL H5SCREATE_SIMPLE_F(3,dimT_file,fle_space,ierr)
      !
      ! Create the dataset with chunking based on the process-based data size.
      CALL H5PCREATE_F(H5P_DATASET_CREATE_F,plist_id,ierr)
      CALL H5PSET_CHUNK_F(plist_id,3,dimW_proc,ierr)
      CALL H5DCREATE_F(hdf5_id,dname,h5prec,fle_space,dset_id,ierr,plist_id)
      CALL H5PCLOSE_F(plist_id,ierr)
      !
      ! Select the region of the filespace dataspace for this process.
      CALL H5SSELECT_HYPERSLAB_F(fle_space,H5S_SELECT_SET_F,oset_file, &
                                 dimW_proc,ierr)

      ! Create the property list for collective dataset writing.
      CALL H5PCREATE_F(H5P_DATASET_XFER_F,plist_id,ierr)
      CALL H5PSET_DXPL_MPIO_F(plist_id,H5FD_MPIO_COLLECTIVE_F,ierr)
      !
      ! Write the data collectively.
      CALL H5DWRITE_F(dset_id,h5prec,dat,dimT_proc,ierr, &
                      file_space_id=fle_space,mem_space_id=mem_space, &
                      xfer_prp=plist_id)

      ! Close the property list.
      CALL H5PCLOSE_F(plist_id,ierr)
      !
      ! Write a real attribute if it is present.
      IF (PRESENT(real_att_name) .AND. PRESENT(real_att_val)) THEN
         CALL WRITE_ATTRIBUTE(dset_id, real_att_name, real_att_val)
      END IF
      !
      ! Close the dataset.
      CALL H5DCLOSE_F(dset_id,ierr)
      !
      ! Close the dataspaces.
      CALL H5SCLOSE_F(fle_space,ierr)
      CALL H5SCLOSE_F(mem_space,ierr)
   END SUBROUTINE WRITE_DATA_PARALLEL_REAL32_RANK3

   !> Procedure to write 1D double precision real data to file in parallel.
   !!
   !! Note that if the dimsW argument is passed, in the first 1:dimsW elements
   !! of the data array are written to file. The user must ensure the offset,
   !! etc. are properly set for each process when writing either the entire
   !! dimsP or reduced dimsW amount of data.
   !!
   !> @param[in] hdf5_id HDF5 ID for the object where data will be written.
   !> @param[in] dname Output name for the data.
   !> @param[in] dimsT Dimensions of the complete dataset across all procs.
   !> @param[in] dimsP Size of the data array on this process.
   !> @param[in] offset Offset for this process in the joined checkpoint.
   !> @param[in] counter Number of dimsP blocks for this process.
   !> @param[in] stride Strides to take in the HDF5 data in each direction.
   !> @param[in] dat Array with dimensions given by dimsP.
   !> @param[in,optional] chunk Chunking dimension, if desired.
   !> @param[in,optional] dimsW Section of data to write for this process.
   SUBROUTINE WRITE_DATA_PARALLEL_REAL64_RANK1(hdf5_id, dname, dimsT, dimsP, &
                                               offset, counter, stride, dat, &
                                               chunk, dimsW)
      ! Required modules.
      USE ISO_FORTRAN_ENV,ONLY: REAL64
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=HID_T),INTENT(IN) :: hdf5_id
      CHARACTER(LEN=*),INTENT(IN) :: dname
      INTEGER(KIND=HSIZE_T),DIMENSION(1),INTENT(IN) :: dimsT, dimsP, offset, &
                                                       counter, stride
      REAL(KIND=REAL64),DIMENSION(dimsP(1)),INTENT(IN) :: dat
      INTEGER(KIND=HSIZE_T),DIMENSION(1),OPTIONAL,INTENT(IN) :: chunk, dimsW
      ! Local variables.
      ! HDF5 identifier for the variable type.
      INTEGER(KIND=HID_T) :: h5prec
      ! HDF5 dataspace identifier for the file space.
      INTEGER(KIND=HID_T) :: filespace
      ! HDF5 dataspace identifier for the memory on this process to be written.
      INTEGER(KIND=HID_T) :: memspace
      ! HDF5 property list identifier.
      INTEGER(KIND=HID_T) :: plist_id
      ! HDF5 dataset identifier.
      INTEGER(KIND=HID_T) :: dataset_id
      ! Local arrays to trim the process memory space if dimsW is present.
      INTEGER(KIND=HSIZE_T),DIMENSION(1) :: oset_tmp
      ! Error handling.
      INTEGER :: ierr
      !
      ! Identify the HDF5 type.
      h5prec = H5KIND_TO_TYPE(REAL64, H5_REAL_KIND)
      !
      ! Create the dataspace for the file space.
      CALL H5SCREATE_SIMPLE_F(1, dimsT, filespace, ierr)
      !
      ! Create the dataspace for the memory on this process.
      CALL H5SCREATE_SIMPLE_F(1, dimsP, memspace, ierr)
      !
      ! Create the dataset with or without chunking, depending on whether or not
      ! the optional "chunk" argument was passed in.
      IF (PRESENT(chunk)) THEN
         CALL H5PCREATE_F(H5P_DATASET_CREATE_F, plist_id, ierr)
         CALL H5PSET_CHUNK_F(plist_id, 1, chunk, ierr)
         CALL H5DCREATE_F(hdf5_id, dname, h5prec, filespace, dataset_id, &
                          ierr, plist_id)
         CALL H5PCLOSE_F(plist_id, ierr)
      ELSE
         CALL H5DCREATE_F(hdf5_id, dname, h5prec, filespace, dataset_id, ierr)
      END IF
      !
      ! Close the dataspace for the file space.
      CALL H5SCLOSE_F(filespace, ierr)
      !
      ! Select the hyperslab in the file. If the dimsW optional argument is
      ! present, we write out the first dimsW elements of the memory on this
      ! process. If it is not present, we write out the entire dimsP memory.
      IF (PRESENT(dimsW)) THEN
         CALL H5DGET_SPACE_F(dataset_id, filespace, ierr)
         CALL H5SSELECT_HYPERSLAB_F(filespace, H5S_SELECT_SET_F, offset, &
                                    counter, ierr, stride, dimsW)
         !
         ! If dimsW is present, we also need to trim the memory space for the
         ! local process to only write out the first dimsW elements.
         oset_tmp(1) = 0_HSIZE_T
         CALL H5SSELECT_HYPERSLAB_F(memspace, H5S_SELECT_SET_F, oset_tmp, &
                                    dimsW, ierr)
      ELSE
         CALL H5DGET_SPACE_F(dataset_id, filespace, ierr)
         CALL H5SSELECT_HYPERSLAB_F(filespace, H5S_SELECT_SET_F, offset, &
                                    counter, ierr, stride, dimsP)
      END IF
      !
      ! Create the property list for collective dataset writing.
      CALL H5PCREATE_F(H5P_DATASET_XFER_F, plist_id, ierr)
      CALL H5PSET_DXPL_MPIO_F(plist_id, H5FD_MPIO_COLLECTIVE_F, ierr)
      !
      ! Write the dataset collectively.
      CALL H5DWRITE_F(dataset_id, h5prec, dat, dimsP, ierr, &
                      file_space_id=filespace, mem_space_id=memspace, &
                      xfer_prp=plist_id)
      !
      ! Close the dataset.
      CALL H5DCLOSE_F(dataset_id, ierr)
      !
      ! Close the property list.
      CALL H5PCLOSE_F(plist_id, ierr)
      !
      ! Close the dataspaces.
      CALL H5SCLOSE_F(filespace, ierr)
      CALL H5SCLOSE_F(memspace, ierr)
   END SUBROUTINE WRITE_DATA_PARALLEL_REAL64_RANK1

   !> Procedure to write 3D double precision real data to file in parallel.
   !!
   !> @param[in] hdf5_id HDF5 ID for the object where data will be written.
   !> @param[in] dname Output name for the data.
   !> @param[in] dimT_proc Total array size for data on this processor.
   !> @param[in] dimW_proc Dimensions of data on this process to write out.
   !> @param[in] oset_proc Offset for the data to write from this process.
   !> @param[in] dimT_file Total size of the data to write out.
   !> @param[in] oset_file Offset for the data to be written from this process.
   !> @param[in] dat Array with dimensions given by dimT_proc.
   !> @param[in,optional] real_att_name Name of a real attributes.
   !> @param[in,optional] real_att_val Value for a real attributes.
   SUBROUTINE WRITE_DATA_PARALLEL_REAL64_RANK3(hdf5_id, dname, dimT_proc, &
                                               dimW_proc, oset_proc, &
                                               dimT_file, oset_file, dat, &
                                               real_att_name, real_att_val)
      ! Required modules.
      USE ISO_FORTRAN_ENV,ONLY: REAL64
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=HID_T),INTENT(IN) :: hdf5_id
      CHARACTER(LEN=*),INTENT(IN) :: dname
      INTEGER(KIND=HSIZE_T),DIMENSION(3),INTENT(IN) :: dimT_proc, dimW_proc, &
                                                       oset_proc, dimT_file, &
                                                       oset_file
      REAL(KIND=REAL64),DIMENSION(dimT_proc(1),dimT_proc(2),dimT_proc(3)), &
         INTENT(IN) :: dat
      CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: real_att_name
      REAL(KIND=REAL64),OPTIONAL,INTENT(IN) :: real_att_val
      ! Local variables.
      ! HDF5 identifier for the variable type.
      INTEGER(KIND=HID_T) :: h5prec
      ! HDF5 dataspace identifier for the file space.
      INTEGER(KIND=HID_T) :: fle_space
      ! HDF5 dataspace identifier for the memory on this process to be written.
      INTEGER(KIND=HID_T) :: mem_space
      ! HDF5 property list identifier.
      INTEGER(KIND=HID_T) :: plist_id
      ! HDF5 dataset identifier.
      INTEGER(KIND=HID_T) :: dset_id
      ! Error handling.
      INTEGER :: ierr

      ! Identify the HDF5 type to write.
      h5prec = H5KIND_TO_TYPE(REAL64,H5_REAL_KIND)

      ! Create the dataspace for the memory on this process.
      CALL H5SCREATE_SIMPLE_F(3,dimT_proc,mem_space,ierr)
      !
      ! Select the region of the memory space to write out.
      CALL H5SSELECT_HYPERSLAB_F(mem_space,H5S_SELECT_SET_F,oset_proc, &
                                 dimW_proc,ierr)

      ! Create the dataspace for the file.
      CALL H5SCREATE_SIMPLE_F(3,dimT_file,fle_space,ierr)
      !
      ! Create the dataset with chunking based on the process-based data size.
      CALL H5PCREATE_F(H5P_DATASET_CREATE_F,plist_id,ierr)
      CALL H5PSET_CHUNK_F(plist_id,3,dimW_proc,ierr)
      CALL H5DCREATE_F(hdf5_id,dname,h5prec,fle_space,dset_id,ierr,plist_id)
      CALL H5PCLOSE_F(plist_id,ierr)
      !
      ! Select the region of the filespace dataspace for this process.
      CALL H5SSELECT_HYPERSLAB_F(fle_space,H5S_SELECT_SET_F,oset_file, &
                                 dimW_proc,ierr)

      ! Create the property list for collective dataset writing.
      CALL H5PCREATE_F(H5P_DATASET_XFER_F,plist_id,ierr)
      CALL H5PSET_DXPL_MPIO_F(plist_id,H5FD_MPIO_COLLECTIVE_F,ierr)
      !
      ! Write the data collectively.
      CALL H5DWRITE_F(dset_id,h5prec,dat,dimT_proc,ierr, &
                      file_space_id=fle_space,mem_space_id=mem_space, &
                      xfer_prp=plist_id)

      ! Close the property list.
      CALL H5PCLOSE_F(plist_id,ierr)
      !
      ! Write a real attribute if it is present.
      IF (PRESENT(real_att_name) .AND. PRESENT(real_att_val)) THEN
         CALL WRITE_ATTRIBUTE(dset_id, real_att_name, real_att_val)
      END IF
      !
      ! Close the dataset.
      CALL H5DCLOSE_F(dset_id,ierr)
      !
      ! Close the dataspaces.
      CALL H5SCLOSE_F(fle_space,ierr)
      CALL H5SCLOSE_F(mem_space,ierr)
   END SUBROUTINE WRITE_DATA_PARALLEL_REAL64_RANK3

   !> Procedure to write a 2D slice from 3D single precision data.
   !!
   !> @param[in] dims Dimensions of the data.
   !> @param[in] dname Output name for the data.
   !> @param[in] dat Array with dimensions given by dims.
   !> @param[in] hdf5_id HDF5 ID for the object where data will be written.
   !> @param[in] srank Rank of the slice to take.
   !> @param[in] sind Index of the slice in the rank dimension.
   SUBROUTINE WRITE_SLICE_REAL32_RANK3(dims, dname, dat, hdf5_id, srank, sind)
      ! Required modules.
      USE ISO_FORTRAN_ENV,ONLY: REAL32
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=HSIZE_T),DIMENSION(3),INTENT(IN) :: dims
      CHARACTER(LEN=*),INTENT(IN) :: dname
      REAL(KIND=REAL32),DIMENSION(dims(1),dims(2),dims(3)),INTENT(IN) :: dat
      INTEGER(KIND=HID_T),INTENT(IN) :: hdf5_id
      INTEGER,INTENT(INOUT) :: srank, sind
      ! Local variables.
      ! HDF5 identifier for the variable type.
      INTEGER(KIND=HID_T) :: h5prec
      ! Dimensions for the slice in the filespace.
      INTEGER(KIND=HSIZE_T),DIMENSION(2) :: dims_out
      ! Dataspace in the file.
      INTEGER(KIND=HID_T) :: fspace_id
      ! Dataspace for the data being passed in.
      INTEGER(KIND=HID_T) :: memspace_id
      ! HDF5 identifier for the dataspace and dataset.
      INTEGER(KIND=HID_T) :: dset_id
      ! Offset in the incoming array to start writing data at.
      INTEGER(KIND=HSIZE_T),DIMENSION(3) :: offset
      ! Stride to take between elements in the data.
      INTEGER(KIND=HSIZE_T),DIMENSION(3) :: stride
      ! Number of blocks to select from the data array.
      INTEGER(KIND=HSIZE_T),DIMENSION(3) :: counter
      ! Dimensions of the data to select from the array.
      INTEGER(KIND=HSIZE_T),DIMENSION(3) :: dims_select
      ! Error handling.
      INTEGER :: ierr
      !
      ! Identify the HDF5 type.
      h5prec = H5KIND_TO_TYPE(REAL32, H5_REAL_KIND)
      !
      ! Determine the output dimensions in the HDF5 file based on the desired
      ! rank through which the slice is taken.
      IF (srank .EQ. 1) THEN
         dims_out = [dims(2), dims(3)]
      ELSE IF (srank .EQ. 2) THEN
         dims_out = [dims(1), dims(3)]
      ELSE IF (srank .EQ. 3) THEN
         dims_out = [dims(1), dims(2)]
      END IF
      !
      ! Create the dataspace corresponding to the file space.
      CALL H5SCREATE_SIMPLE_F(2, dims_out, fspace_id, ierr)
      !
      ! Create the dataset in hdf5_id to hold the output.
      CALL H5DCREATE_F(hdf5_id, dname, h5prec, fspace_id, dset_id, ierr)
      !
      ! Create a dataspace corresponding to the memory of the local 3D array.
      CALL H5SCREATE_SIMPLE_F(3, dims, memspace_id, ierr)
      !
      ! Select the hyperslab from memory that we actually want to write.
      IF (srank .EQ. 1) THEN
         offset = [INT(sind-1, HSIZE_T), 0_HSIZE_T, 0_HSIZE_T]
         counter = [1_HSIZE_T, 1_HSIZE_T, 1_HSIZE_T]
         stride = [1_HSIZE_T, 1_HSIZE_T, 1_HSIZE_T]
         dims_select = [1_HSIZE_T, dims(2), dims(3)]
      ELSE IF (srank .EQ. 2) THEN
         offset = [0_HSIZE_T, INT(sind-1, HSIZE_T), 0_HSIZE_T]
         counter = [1_HSIZE_T, 1_HSIZE_T, 1_HSIZE_T]
         stride = [1_HSIZE_T, 1_HSIZE_T, 1_HSIZE_T]
         dims_select = [dims(1), 1_HSIZE_T, dims(3)]
      ELSE IF (srank .EQ. 3) THEN
         offset = [0_HSIZE_T, 0_HSIZE_T, INT(sind-1, HSIZE_T)]
         counter = [1_HSIZE_T, 1_HSIZE_T, 1_HSIZE_T]
         stride = [1_HSIZE_T, 1_HSIZE_T, 1_HSIZE_T]
         dims_select = [dims(1), dims(2), 1_HSIZE_T]
      END IF
      !
      ! Select the hyperslab of memory to write out.
      CALL H5SSELECT_HYPERSLAB_F(memspace_id, H5S_SELECT_SET_F, offset, &
                                 counter, ierr, stride, dims_select)
      !
      ! Write out the dataset.
      CALL H5DWRITE_F(dset_id, h5prec, dat, dims, ierr, &
                      file_space_id=fspace_id, mem_space_id=memspace_id)
      !
      ! Write attributes to the dataset indicating what plane it is from.
      CALL WRITE_ATTRIBUTE(dset_id, 'srank', srank)
      CALL WRITE_ATTRIBUTE(dset_id, 'sind', sind)
      !
      ! Close the dataset.
      CALL H5DCLOSE_F(dset_id, ierr)
      !
      ! Close the dataspaces.
      CALL H5SCLOSE_F(fspace_id, ierr)
      CALL H5SCLOSE_F(memspace_id, ierr)
   END SUBROUTINE WRITE_SLICE_REAL32_RANK3

   !> Procedure to write a slice of a 3D single precision dataset in parallel.
   !!
   !! A current limitation of this routine is that dimW_proc must be the same
   !! on all processors writing data in parallel. This is because dimW_proc is
   !! used as the chunking dimensions for the dataset, which must be the same
   !! accross all processors. Change this in the future if the slice is to be
   !! written with different amounts of data from each process.
   !!
   !! The user must be aware of the current process layout when making a call to
   !! this routine, and make sure that only those processes that contain the
   !! desired slice are participating. For example, when the pencils are in
   !! y-aligned form, and we want the ky=0 slice, all processors in the 2D MPI
   !! process layout should participate in writing out the slice. On the other
   !! hand, consider writing out a kx=0 slice with y-aligned pencils. In this
   !! instance, only the processes in the first column communicator should
   !! be participating in IO, as only they own the data for the kx=0 slice.
   !! Currently, this routine has only been tested for writing out the ky=0
   !! plane in a y-aligned pencil layout.
   !!
   !> @param[in] hdf5_id HDF5 ID for the object where data will be written.
   !> @param[in] dname Output name for the data.
   !> @param[in] dimT_proc Dimensions of the complete dataset on this process.
   !> @param[in] dimW_proc Dimensions to write from this process.
   !> @param[in] oset_proc Offsets for the data to write from this process.
   !> @param[in] dimT_comm Total dims. of the 3D dataset for the communicator.
   !> @param[in] oset_comm 3D offsets for this process in the communicator.
   !> @param[in] dat Data array with dimensions given by dimT_proc
   !> @param[in] srank Rank of the slice in the data array.
   !> @param[in] sind Index for the slice in the rank dimension.
   SUBROUTINE WRITE_SLICE_PARALLEL_REAL32_RANK3(hdf5_id,dname, &
                                                dimT_proc,dimW_proc,oset_proc, &
                                                dimT_comm,oset_comm, &
                                                dat,srank,sind)
      ! Required modules.
      USE ISO_FORTRAN_ENV,ONLY: REAL32
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=HID_T),INTENT(IN) :: hdf5_id
      CHARACTER(LEN=*),INTENT(IN) :: dname
      INTEGER(KIND=HSIZE_T),DIMENSION(3),INTENT(IN) :: dimT_proc,dimW_proc
      INTEGER(KIND=HSIZE_T),DIMENSION(3),INTENT(IN) :: oset_proc
      INTEGER(KIND=HSIZE_T),DIMENSION(3),INTENT(IN) :: dimT_comm,oset_comm
      REAL(KIND=REAL32),DIMENSION(dimT_proc(1),dimT_proc(2),dimT_proc(3)), &
         INTENT(IN) :: dat
      INTEGER,INTENT(INOUT) :: srank,sind
      ! Local variables.
      ! HDF5 identifier for the variable type.
      INTEGER(KIND=HID_T) :: h5prec
      ! Dataspace in the file.
      INTEGER(KIND=HID_T) :: fspace_id
      ! Dataspace for the data being written from this process.
      INTEGER(KIND=HID_T) :: mspace_id
      ! HDF5 identifier for the dataset.
      INTEGER(KIND=HID_T) :: dset_id
      ! HDF5 property list identifier to control parallel IO access.
      INTEGER(KIND=HID_T) :: plist_id
      ! Arrays used to form the file space dataset, and to select the hyperslab
      ! for this particular process from the file space dataset.
      INTEGER(KIND=HSIZE_T),DIMENSION(2) :: dimT_fs,dimP_fs,oset_fs
      ! Arrays used to select the hyperslab from the data on this process.
      INTEGER(KIND=HSIZE_T),DIMENSION(3) :: oset_ms,dimW_ms
      ! Error handling.
      INTEGER :: ierr
      !
      ! Identify the HDF5 type.
      h5prec = H5KIND_TO_TYPE(REAL32, H5_REAL_KIND)
      !
      ! Determine the output dimensions, offsets, etc. in the HDF5 file space
      ! based on the desired rank through which the slice is taken. The file
      ! space is treated as one rank less than the data array passed in.
      IF (srank .EQ. 1) THEN
         dimT_fs = [dimT_comm(2), dimT_comm(3)]
         dimP_fs = [dimW_proc(2), dimW_proc(3)]
         oset_fs = [oset_comm(2), oset_comm(3)]
      ELSE IF (srank .EQ. 2) THEN
         dimT_fs = [dimT_comm(1), dimT_comm(3)]
         dimP_fs = [dimW_proc(1), dimW_proc(3)]
         oset_fs = [oset_comm(1), oset_comm(3)]
      ELSE IF (srank .EQ. 3) THEN
         dimT_fs = [dimT_comm(1), dimT_comm(2)]
         dimP_fs = [dimW_proc(1), dimW_proc(2)]
         oset_fs = [oset_comm(1), oset_comm(2)]
      END IF
      !
      ! Create the dataspace corresponding to the file space.
      CALL H5SCREATE_SIMPLE_F(2,dimT_fs,fspace_id,ierr)
      !
      ! Create the chunked dataset in the file space. The chunking dimensions
      ! are based on the process data size, which must be consistent across all
      ! processors writing out data.
      CALL H5PCREATE_F(H5P_DATASET_CREATE_F,plist_id,ierr)
      CALL H5PSET_CHUNK_F(plist_id,2,dimP_fs,ierr)
      CALL H5DCREATE_F(hdf5_id,dname,h5prec,fspace_id,dset_id,ierr,plist_id)
      CALL H5PCLOSE_F(plist_id,ierr)
      !
      ! Select the slab from the file space that this process will write to.
      CALL H5SSELECT_HYPERSLAB_F(fspace_id,H5S_SELECT_SET_F,oset_fs, &
                                 dimP_fs,ierr)
      !
      ! Create a dataspace corresponding to the memory of the local 3D array.
      CALL H5SCREATE_SIMPLE_F(3,dimT_proc,mspace_id,ierr)
      !
      ! Select the hyperslab from memory that we actually want to write.
      IF (srank .EQ. 1) THEN
         oset_ms = [INT(sind-1, HSIZE_T), 0_HSIZE_T, 0_HSIZE_T]
         dimW_ms = [1_HSIZE_T, dimW_proc(2), dimW_proc(3)]
      ELSE IF (srank .EQ. 2) THEN
         oset_ms = [0_HSIZE_T, INT(sind-1, HSIZE_T), 0_HSIZE_T]
         dimW_ms = [dimW_proc(1), 1_HSIZE_T, dimW_proc(3)]
      ELSE IF (srank .EQ. 3) THEN
         oset_ms = [0_HSIZE_T, 0_HSIZE_T, INT(sind-1, HSIZE_T)]
         dimW_ms = [dimW_proc(1), dimW_proc(2), 1_HSIZE_T]
      END IF
      !
      ! Select the hyperslab of memory to write out.
      CALL H5SSELECT_HYPERSLAB_F(mspace_id,H5S_SELECT_SET_F,oset_ms, &
                                 dimW_ms,ierr)
      !
      ! Create the property list for collective dataset writing.
      CALL H5PCREATE_F(H5P_DATASET_XFER_F,plist_id,ierr)
      CALL H5PSET_DXPL_MPIO_F(plist_id,H5FD_MPIO_COLLECTIVE_F,ierr)
      !
      ! Write the dataset collectively.
      CALL H5DWRITE_F(dset_id,h5prec,dat,dimT_proc,ierr, &
                      file_space_id=fspace_id,mem_space_id=mspace_id, &
                      xfer_prp=plist_id)
      !
      ! Write attributes to the dataset indicating what plane it is from.
      CALL WRITE_ATTRIBUTE(dset_id,'srank',srank)
      CALL WRITE_ATTRIBUTE(dset_id,'sind',sind)
      !
      ! Close the dataset.
      CALL H5DCLOSE_F(dset_id,ierr)
      !
      ! Close the property list.
      CALL H5PCLOSE_F(plist_id,ierr)
      !
      ! Close the dataspaces.
      CALL H5SCLOSE_F(fspace_id,ierr)
      CALL H5SCLOSE_F(mspace_id,ierr)
   END SUBROUTINE WRITE_SLICE_PARALLEL_REAL32_RANK3

   !> Procedure to write a 2D slice from 3D double precision data.
   !!
   !> @param[in] dims Dimensions of the data.
   !> @param[in] dname Output name for the data.
   !> @param[in] dat Array with dimensions given by dims.
   !> @param[in] hdf5_id HDF5 ID for the object where data will be written.
   !> @param[in] srank Rank of the slice to take.
   !> @param[in] sind Index of the slice in the rank dimension.
   SUBROUTINE WRITE_SLICE_REAL64_RANK3(dims, dname, dat, hdf5_id, srank, sind)
      ! Required modules.
      USE ISO_FORTRAN_ENV,ONLY: REAL64
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=HSIZE_T),DIMENSION(3),INTENT(IN) :: dims
      CHARACTER(LEN=*),INTENT(IN) :: dname
      REAL(KIND=REAL64),DIMENSION(dims(1),dims(2),dims(3)),INTENT(IN) :: dat
      INTEGER(KIND=HID_T),INTENT(IN) :: hdf5_id
      INTEGER,INTENT(INOUT) :: srank, sind
      ! Local variables.
      ! HDF5 identifier for the variable type.
      INTEGER(KIND=HID_T) :: h5prec
      ! Dimensions for the slice in the filespace.
      INTEGER(KIND=HSIZE_T),DIMENSION(2) :: dims_out
      ! Dataspace in the file.
      INTEGER(KIND=HID_T) :: fspace_id
      ! Dataspace for the data being passed in.
      INTEGER(KIND=HID_T) :: memspace_id
      ! HDF5 identifier for the dataspace and dataset.
      INTEGER(KIND=HID_T) :: dset_id
      ! Offset in the incoming array to start writing data at.
      INTEGER(KIND=HSIZE_T),DIMENSION(3) :: offset
      ! Stride to take between elements in the data.
      INTEGER(KIND=HSIZE_T),DIMENSION(3) :: stride
      ! Number of blocks to select from the data array.
      INTEGER(KIND=HSIZE_T),DIMENSION(3) :: counter
      ! Dimensions of the data to select from the array.
      INTEGER(KIND=HSIZE_T),DIMENSION(3) :: dims_select
      ! Error handling.
      INTEGER :: ierr
      !
      ! Identify the HDF5 type.
      h5prec = H5KIND_TO_TYPE(REAL64, H5_REAL_KIND)
      !
      ! Determine the output dimensions in the HDF5 file based on the desired
      ! rank through which the slice is taken.
      IF (srank .EQ. 1) THEN
         dims_out = [dims(2), dims(3)]
      ELSE IF (srank .EQ. 2) THEN
         dims_out = [dims(1), dims(3)]
      ELSE IF (srank .EQ. 3) THEN
         dims_out = [dims(1), dims(2)]
      END IF
      !
      ! Create the dataspace corresponding to the file space.
      CALL H5SCREATE_SIMPLE_F(2, dims_out, fspace_id, ierr)
      !
      ! Create the dataset in hdf5_id to hold the output.
      CALL H5DCREATE_F(hdf5_id, dname, h5prec, fspace_id, dset_id, ierr)
      !
      ! Create a dataspace corresponding to the memory of the local 3D array.
      CALL H5SCREATE_SIMPLE_F(3, dims, memspace_id, ierr)
      !
      ! Select the hyperslab from memory that we actually want to write.
      IF (srank .EQ. 1) THEN
         offset = [INT(sind-1, HSIZE_T), 0_HSIZE_T, 0_HSIZE_T]
         counter = [1_HSIZE_T, 1_HSIZE_T, 1_HSIZE_T]
         stride = [1_HSIZE_T, 1_HSIZE_T, 1_HSIZE_T]
         dims_select = [1_HSIZE_T, dims(2), dims(3)]
      ELSE IF (srank .EQ. 2) THEN
         offset = [0_HSIZE_T, INT(sind-1, HSIZE_T), 0_HSIZE_T]
         counter = [1_HSIZE_T, 1_HSIZE_T, 1_HSIZE_T]
         stride = [1_HSIZE_T, 1_HSIZE_T, 1_HSIZE_T]
         dims_select = [dims(1), 1_HSIZE_T, dims(3)]
      ELSE IF (srank .EQ. 3) THEN
         offset = [0_HSIZE_T, 0_HSIZE_T, INT(sind-1, HSIZE_T)]
         counter = [1_HSIZE_T, 1_HSIZE_T, 1_HSIZE_T]
         stride = [1_HSIZE_T, 1_HSIZE_T, 1_HSIZE_T]
         dims_select = [dims(1), dims(2), 1_HSIZE_T]
      END IF
      !
      ! Select the hyperslab of memory to write out.
      CALL H5SSELECT_HYPERSLAB_F(memspace_id, H5S_SELECT_SET_F, offset, &
                                 counter, ierr, stride, dims_select)
      !
      ! Write out the dataset.
      CALL H5DWRITE_F(dset_id, h5prec, dat, dims, ierr, &
                      file_space_id=fspace_id, mem_space_id=memspace_id)
      !
      ! Write attributes to the dataset indicating what plane it is from.
      CALL WRITE_ATTRIBUTE(dset_id, 'srank', srank)
      CALL WRITE_ATTRIBUTE(dset_id, 'sind', sind)
      !
      ! Close the dataset.
      CALL H5DCLOSE_F(dset_id, ierr)
      !
      ! Close the dataspaces.
      CALL H5SCLOSE_F(fspace_id, ierr)
      CALL H5SCLOSE_F(memspace_id, ierr)
   END SUBROUTINE WRITE_SLICE_REAL64_RANK3

   !> Procedure to write a slice of a 3D single precision dataset in parallel.
   !!
   !! A current limitation of this routine is that dimW_proc must be the same
   !! on all processors writing data in parallel. This is because dimW_proc is
   !! used as the chunking dimensions for the dataset, which must be the same
   !! accross all processors. Change this in the future if the slice is to be
   !! written with different amounts of data from each process.
   !!
   !! The user must be aware of the current process layout when making a call to
   !! this routine, and make sure that only those processes that contain the
   !! desired slice are participating. For example, when the pencils are in
   !! y-aligned form, and we want the ky=0 slice, all processors in the 2D MPI
   !! process layout should participate in writing out the slice. On the other
   !! hand, consider writing out a kx=0 slice with y-aligned pencils. In this
   !! instance, only the processes in the first column communicator should
   !! be participating in IO, as only they own the data for the kx=0 slice.
   !! Currently, this routine has only been tested for writing out the ky=0
   !! plane in a y-aligned pencil layout.
   !!
   !> @param[in] hdf5_id HDF5 ID for the object where data will be written.
   !> @param[in] dname Output name for the data.
   !> @param[in] dimT_proc Dimensions of the complete dataset on this process.
   !> @param[in] dimW_proc Dimensions to write from this process.
   !> @param[in] oset_proc Offsets for the data to write from this process.
   !> @param[in] dimT_comm Total dims. of the 3D dataset for the communicator.
   !> @param[in] oset_comm 3D offsets for this process in the communicator.
   !> @param[in] dat Data array with dimensions given by dimT_proc
   !> @param[in] srank Rank of the slice in the data array.
   !> @param[in] sind Index for the slice in the rank dimension.
   SUBROUTINE WRITE_SLICE_PARALLEL_REAL64_RANK3(hdf5_id,dname, &
                                                dimT_proc,dimW_proc,oset_proc, &
                                                dimT_comm,oset_comm, &
                                                dat,srank,sind)
      ! Required modules.
      USE ISO_FORTRAN_ENV,ONLY: REAL64
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=HID_T),INTENT(IN) :: hdf5_id
      CHARACTER(LEN=*),INTENT(IN) :: dname
      INTEGER(KIND=HSIZE_T),DIMENSION(3),INTENT(IN) :: dimT_proc,dimW_proc
      INTEGER(KIND=HSIZE_T),DIMENSION(3),INTENT(IN) :: oset_proc
      INTEGER(KIND=HSIZE_T),DIMENSION(3),INTENT(IN) :: dimT_comm,oset_comm
      REAL(KIND=REAL64),DIMENSION(dimT_proc(1),dimT_proc(2),dimT_proc(3)), &
         INTENT(IN) :: dat
      INTEGER,INTENT(INOUT) :: srank,sind
      ! Local variables.
      ! HDF5 identifier for the variable type.
      INTEGER(KIND=HID_T) :: h5prec
      ! Dataspace in the file.
      INTEGER(KIND=HID_T) :: fspace_id
      ! Dataspace for the data being written from this process.
      INTEGER(KIND=HID_T) :: mspace_id
      ! HDF5 identifier for the dataset.
      INTEGER(KIND=HID_T) :: dset_id
      ! HDF5 property list identifier to control parallel IO access.
      INTEGER(KIND=HID_T) :: plist_id
      ! Arrays used to form the file space dataset, and to select the hyperslab
      ! for this particular process from the file space dataset.
      INTEGER(KIND=HSIZE_T),DIMENSION(2) :: dimT_fs,dimP_fs,oset_fs
      ! Arrays used to select the hyperslab from the data on this process.
      INTEGER(KIND=HSIZE_T),DIMENSION(3) :: oset_ms,dimW_ms
      ! Error handling.
      INTEGER :: ierr
      !
      ! Identify the HDF5 type.
      h5prec = H5KIND_TO_TYPE(REAL64, H5_REAL_KIND)
      !
      ! Determine the output dimensions, offsets, etc. in the HDF5 file space
      ! based on the desired rank through which the slice is taken. The file
      ! space is treated as one rank less than the data array passed in.
      IF (srank .EQ. 1) THEN
         dimT_fs = [dimT_comm(2), dimT_comm(3)]
         dimP_fs = [dimW_proc(2), dimW_proc(3)]
         oset_fs = [oset_comm(2), oset_comm(3)]
      ELSE IF (srank .EQ. 2) THEN
         dimT_fs = [dimT_comm(1), dimT_comm(3)]
         dimP_fs = [dimW_proc(1), dimW_proc(3)]
         oset_fs = [oset_comm(1), oset_comm(3)]
      ELSE IF (srank .EQ. 3) THEN
         dimT_fs = [dimT_comm(1), dimT_comm(2)]
         dimP_fs = [dimW_proc(1), dimW_proc(2)]
         oset_fs = [oset_comm(1), oset_comm(2)]
      END IF
      !
      ! Create the dataspace corresponding to the file space.
      CALL H5SCREATE_SIMPLE_F(2,dimT_fs,fspace_id,ierr)
      !
      ! Create the chunked dataset in the file space. The chunking dimensions
      ! are based on the process data size, which must be consistent across all
      ! processors writing out data.
      CALL H5PCREATE_F(H5P_DATASET_CREATE_F,plist_id,ierr)
      CALL H5PSET_CHUNK_F(plist_id,2,dimP_fs,ierr)
      CALL H5DCREATE_F(hdf5_id,dname,h5prec,fspace_id,dset_id,ierr,plist_id)
      CALL H5PCLOSE_F(plist_id,ierr)
      !
      ! Select the slab from the file space that this process will write to.
      CALL H5SSELECT_HYPERSLAB_F(fspace_id,H5S_SELECT_SET_F,oset_fs, &
                                 dimP_fs,ierr)
      !
      ! Create a dataspace corresponding to the memory of the local 3D array.
      CALL H5SCREATE_SIMPLE_F(3,dimT_proc,mspace_id,ierr)
      !
      ! Select the hyperslab from memory that we actually want to write.
      IF (srank .EQ. 1) THEN
         oset_ms = [INT(sind-1, HSIZE_T), 0_HSIZE_T, 0_HSIZE_T]
         dimW_ms = [1_HSIZE_T, dimW_proc(2), dimW_proc(3)]
      ELSE IF (srank .EQ. 2) THEN
         oset_ms = [0_HSIZE_T, INT(sind-1, HSIZE_T), 0_HSIZE_T]
         dimW_ms = [dimW_proc(1), 1_HSIZE_T, dimW_proc(3)]
      ELSE IF (srank .EQ. 3) THEN
         oset_ms = [0_HSIZE_T, 0_HSIZE_T, INT(sind-1, HSIZE_T)]
         dimW_ms = [dimW_proc(1), dimW_proc(2), 1_HSIZE_T]
      END IF
      !
      ! Select the hyperslab of memory to write out.
      CALL H5SSELECT_HYPERSLAB_F(mspace_id,H5S_SELECT_SET_F,oset_ms, &
                                 dimW_ms,ierr)
      !
      ! Create the property list for collective dataset writing.
      CALL H5PCREATE_F(H5P_DATASET_XFER_F,plist_id,ierr)
      CALL H5PSET_DXPL_MPIO_F(plist_id,H5FD_MPIO_COLLECTIVE_F,ierr)
      !
      ! Write the dataset collectively.
      CALL H5DWRITE_F(dset_id,h5prec,dat,dimT_proc,ierr, &
                      file_space_id=fspace_id,mem_space_id=mspace_id, &
                      xfer_prp=plist_id)
      !
      ! Write attributes to the dataset indicating what plane it is from.
      CALL WRITE_ATTRIBUTE(dset_id,'srank',srank)
      CALL WRITE_ATTRIBUTE(dset_id,'sind',sind)
      !
      ! Close the dataset.
      CALL H5DCLOSE_F(dset_id,ierr)
      !
      ! Close the property list.
      CALL H5PCLOSE_F(plist_id,ierr)
      !
      ! Close the dataspaces.
      CALL H5SCLOSE_F(fspace_id,ierr)
      CALL H5SCLOSE_F(mspace_id,ierr)
   END SUBROUTINE WRITE_SLICE_PARALLEL_REAL64_RANK3

   !> Procedure to read hyperslabs of data into contiguous buffer arrays.
   !!
   !> @param[in] hdf5_id HDF5 identifier to the object containing the dataset.
   !> @param[in] dims Dimensions of the buffer array.
   !> @param[out] buf Buffer array to read the data into.
   !> @param[in] dname Name of the dataset in the HDF5 file.
   !> @param[in] offset Offset for the hyperslab in the checkpoint file.
   SUBROUTINE READ_DATA_INT32_RANK1(hdf5_id, dims, buf, dname, offset)
      ! Required modules.
      USE ISO_FORTRAN_ENV,ONLY: INT32
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=HID_T),INTENT(IN) :: hdf5_id
      INTEGER(KIND=HSIZE_T),DIMENSION(1),INTENT(IN) :: dims
      INTEGER(KIND=INT32),DIMENSION(dims(1)),INTENT(OUT) :: buf
      CHARACTER(LEN=*),INTENT(IN) :: dname
      INTEGER(KIND=HSIZE_T),DIMENSION(1),INTENT(IN) :: offset
      ! Local variables.
      ! HDF5 identifier for the variable type.
      INTEGER(KIND=HID_T) :: h5prec
      ! HDF5 identifiers for the dataspaces and data sets.
      INTEGER(KIND=HID_T) :: dset_id, dspace_id, memspace
      ! Rank of the memory dataspace.
      INTEGER :: memrank = 1
      ! Error handling.
      INTEGER :: ierr
      !
      ! Identify the HDF5 type.
      h5prec = H5KIND_TO_TYPE(INT32, H5_INTEGER_KIND)
      !
      ! Open the dataset.
      CALL H5DOPEN_F(hdf5_id, dname, dset_id, ierr)
      !
      ! Get the dataset's datspace identifier.
      CALL H5DGET_SPACE_F(dset_id, dspace_id, ierr)
      !
      ! Select the hyperslab in the dataset.
      CALL H5SSELECT_HYPERSLAB_F(dspace_id, H5S_SELECT_SET_F, offset, &
                                 dims, ierr)
      !
      ! Create the memory dataspace.
      CALL H5SCREATE_SIMPLE_F(memrank, dims, memspace, ierr)
      !
      ! Read data from hyperslab in the file into the hyperslab in memory.
      CALL H5DREAD_F(dset_id, h5prec, buf, dims, ierr, memspace, dspace_id)
      !
      ! Close the memory space.
      CALL H5SCLOSE_F(memspace, ierr)
      !
      ! Close the dataspace for the dataset.
      CALL H5SCLOSE_F(dspace_id, ierr)
      !
      ! Close the dataset.
      CALL H5DCLOSE_F(dset_id, ierr)
   END SUBROUTINE READ_DATA_INT32_RANK1

   !> Procedure to read hyperslabs of data into contiguous buffer arrays.
   !!
   !> @param[in] hdf5_id HDF5 identifier to the object containing the dataset.
   !> @param[in] dims Dimensions of the buffer array.
   !> @param[out] buf Buffer array to read the data into.
   !> @param[in] dname Name of the dataset in the HDF5 file.
   !> @param[in] offset Offset for the hyperslab in the checkpoint file.
#ifdef NEW_HDF5
   SUBROUTINE READ_DATA_INT64_RANK1(hdf5_id, dims, buf, dname, offset)
      ! Required modules.
      USE ISO_FORTRAN_ENV,ONLY: INT64
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=HID_T),INTENT(IN) :: hdf5_id
      INTEGER(KIND=HSIZE_T),DIMENSION(1),INTENT(IN) :: dims
      INTEGER(KIND=INT64),DIMENSION(dims(1)),INTENT(OUT) :: buf
      CHARACTER(LEN=*),INTENT(IN) :: dname
      INTEGER(KIND=HSIZE_T),DIMENSION(1),INTENT(IN) :: offset
      ! Local variables.
      ! HDF5 identifier for the variable type.
      INTEGER(KIND=HID_T) :: h5prec
      ! HDF5 identifiers for the dataspaces and data sets.
      INTEGER(KIND=HID_T) :: dset_id, dspace_id, memspace
      ! Rank of the memory dataspace.
      INTEGER :: memrank = 1
      ! Error handling.
      INTEGER :: ierr
      !
      ! Identify the HDF5 type.
      h5prec = H5KIND_TO_TYPE(INT64, H5_INTEGER_KIND)
      !
      ! Open the dataset.
      CALL H5DOPEN_F(hdf5_id, dname, dset_id, ierr)
      !
      ! Get the dataset's datspace identifier.
      CALL H5DGET_SPACE_F(dset_id, dspace_id, ierr)
      !
      ! Select the hyperslab in the dataset.
      CALL H5SSELECT_HYPERSLAB_F(dspace_id, H5S_SELECT_SET_F, offset, &
                                 dims, ierr)
      !
      ! Create the memory dataspace.
      CALL H5SCREATE_SIMPLE_F(memrank, dims, memspace, ierr)
      !
      ! Read data from hyperslab in the file into the hyperslab in memory.
      CALL H5DREAD_F(dset_id, h5prec, buf, dims, ierr, memspace, dspace_id)
      !
      ! Close the memory space.
      CALL H5SCLOSE_F(memspace, ierr)
      !
      ! Close the dataspace for the dataset.
      CALL H5SCLOSE_F(dspace_id, ierr)
      !
      ! Close the dataset.
      CALL H5DCLOSE_F(dset_id, ierr)
   END SUBROUTINE READ_DATA_INT64_RANK1
#endif

   !> Procedure to read hyperslabs of data into contiguous buffer arrays.
   !!
   !> @param[in] hdf5_id HDF5 identifier to the object containing the dataset.
   !> @param[in] dims Dimensions of the buffer array.
   !> @param[out] buf Buffer array to read the data into.
   !> @param[in] dname Name of the dataset in the HDF5 file.
   !> @param[in] offset Offset for the hyperslab in the checkpoint file.
   SUBROUTINE READ_DATA_REAL32_RANK1(hdf5_id, dims, buf, dname, offset)
      ! Required modules.
      USE ISO_FORTRAN_ENV,ONLY: REAL32
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=HID_T),INTENT(IN) :: hdf5_id
      INTEGER(KIND=HSIZE_T),DIMENSION(1),INTENT(IN) :: dims
      REAL(KIND=REAL32),DIMENSION(dims(1)),INTENT(OUT) :: buf
      CHARACTER(LEN=*),INTENT(IN) :: dname
      INTEGER(KIND=HSIZE_T),DIMENSION(1),INTENT(IN) :: offset
      ! Local variables.
      ! HDF5 identifier for the variable type.
      INTEGER(KIND=HID_T) :: h5prec
      ! HDF5 identifiers for the dataspaces and data sets.
      INTEGER(KIND=HID_T) :: dset_id, dspace_id, memspace
      ! Rank of the memory dataspace.
      INTEGER :: memrank = 1
      ! Error handling.
      INTEGER :: ierr
      !
      ! Identify the HDF5 type.
      h5prec = H5KIND_TO_TYPE(REAL32, H5_REAL_KIND)
      !
      ! Open the dataset.
      CALL H5DOPEN_F(hdf5_id, dname, dset_id, ierr)
      !
      ! Get the dataset's datspace identifier.
      CALL H5DGET_SPACE_F(dset_id, dspace_id, ierr)
      !
      ! Select the hyperslab in the dataset.
      CALL H5SSELECT_HYPERSLAB_F(dspace_id, H5S_SELECT_SET_F, offset, &
                                 dims, ierr)
      !
      ! Create the memory dataspace.
      CALL H5SCREATE_SIMPLE_F(memrank, dims, memspace, ierr)
      !
      ! Read data from hyperslab in the file into the hyperslab in memory.
      CALL H5DREAD_F(dset_id, h5prec, buf, dims, ierr, memspace, dspace_id)
      !
      ! Close the memory space.
      CALL H5SCLOSE_F(memspace, ierr)
      !
      ! Close the dataspace for the dataset.
      CALL H5SCLOSE_F(dspace_id, ierr)
      !
      ! Close the dataset.
      CALL H5DCLOSE_F(dset_id, ierr)
   END SUBROUTINE READ_DATA_REAL32_RANK1

   !> Procedure to read hyperslabs of data into contiguous buffer arrays.
   !!
   !> @param[in] hdf5_id HDF5 identifier to the object containing the dataset.
   !> @param[in] dims Dimensions of the buffer array.
   !> @param[out] buf Buffer array to read the data into.
   !> @param[in] dname Name of the dataset in the HDF5 file.
   !> @param[in] offset Offset for the hyperslab in the checkpoint file.
   SUBROUTINE READ_DATA_REAL32_RANK2(hdf5_id, dims, buf, dname, offset)
      ! Required modules.
      USE ISO_FORTRAN_ENV,ONLY: REAL32
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=HID_T),INTENT(IN) :: hdf5_id
      INTEGER(KIND=HSIZE_T),DIMENSION(2),INTENT(IN) :: dims
      REAL(KIND=REAL32),DIMENSION(dims(1),dims(2)),INTENT(OUT) :: buf
      CHARACTER(LEN=*),INTENT(IN) :: dname
      INTEGER(KIND=HSIZE_T),DIMENSION(2),INTENT(IN) :: offset
      ! Local variables.
      ! HDF5 identifier for the variable type.
      INTEGER(KIND=HID_T) :: h5prec
      ! HDF5 identifiers for the dataspaces and data sets.
      INTEGER(KIND=HID_T) :: dset_id, dspace_id, memspace
      ! Rank of the memory dataspace.
      INTEGER :: memrank = 2
      ! Error handling.
      INTEGER :: ierr
      !
      ! Identify the HDF5 type.
      h5prec = H5KIND_TO_TYPE(REAL32, H5_REAL_KIND)
      !
      ! Open the dataset.
      CALL H5DOPEN_F(hdf5_id, dname, dset_id, ierr)
      !
      ! Get the dataset's datspace identifier.
      CALL H5DGET_SPACE_F(dset_id, dspace_id, ierr)
      !
      ! Select the hyperslab in the dataset.
      CALL H5SSELECT_HYPERSLAB_F(dspace_id, H5S_SELECT_SET_F, offset, &
                                 dims, ierr)
      !
      ! Create the memory dataspace.
      CALL H5SCREATE_SIMPLE_F(memrank, dims, memspace, ierr)
      !
      ! Read data from hyperslab in the file into the hyperslab in memory.
      CALL H5DREAD_F(dset_id, h5prec, buf, dims, ierr, memspace, dspace_id)
      !
      ! Close the memory space.
      CALL H5SCLOSE_F(memspace, ierr)
      !
      ! Close the dataspace for the dataset.
      CALL H5SCLOSE_F(dspace_id, ierr)
      !
      ! Close the dataset.
      CALL H5DCLOSE_F(dset_id, ierr)
   END SUBROUTINE READ_DATA_REAL32_RANK2

   !> Procedure to read hyperslabs of data into hyperslabs of buffer arrays.
   !!
   !> @param[in] hdf5_id HDF5 identifier to the object containing the dataset.
   !> @param[in] dims Dimensions of the buffer array.
   !> @param[in,out] buf Buffer array to read the data into.
   !> @param[in] dname Name of the dataset in the HDF5 file.
   !> @param[in] dims_slab Dimensions of the hyperslab.
   !> @param[in] offset_read Offset for the hyperslab in the checkpoint file.
   !> @param[in] offset_fill Offset for the hyperslab in the buffer array.
   SUBROUTINE READ_DATA_PARTIAL_REAL32_RANK2(hdf5_id, dims, buf, dname, &
                                             dims_slab, offset_read, &
                                             offset_fill)
      ! Required modules.
      USE ISO_FORTRAN_ENV,ONLY: REAL32
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=HID_T),INTENT(IN) :: hdf5_id
      INTEGER(KIND=HSIZE_T),DIMENSION(2),INTENT(IN) :: dims
      REAL(KIND=REAL32),DIMENSION(dims(1),dims(2)),INTENT(INOUT) :: buf
      CHARACTER(LEN=*),INTENT(IN) :: dname
      INTEGER(KIND=HSIZE_T),DIMENSION(2),INTENT(IN) :: dims_slab
      INTEGER(KIND=HSIZE_T),DIMENSION(2),INTENT(IN) :: offset_read, offset_fill
      ! Local variables.
      ! HDF5 identifier for the variable type.
      INTEGER(KIND=HID_T) :: h5prec
      ! HDF5 identifiers for the dataspaces, memory spaces, and data sets.
      INTEGER(KIND=HID_T) :: dset_id, dspace_id, memspace
      ! Rank of the memory dataspace.
      INTEGER :: memrank = 2
      ! Error handling.
      INTEGER :: ierr
      !
      ! Identify the HDF5 type.
      h5prec = H5KIND_TO_TYPE(REAL32, H5_REAL_KIND)
      !
      ! Open the dataset.
      CALL H5DOPEN_F(hdf5_id, dname, dset_id, ierr)
      !
      ! Get the dataset's datspace identifier.
      CALL H5DGET_SPACE_F(dset_id, dspace_id, ierr)
      !
      ! Select the hyperslab in the dataset.
      CALL H5SSELECT_HYPERSLAB_F(dspace_id, H5S_SELECT_SET_F, offset_read, &
                                 dims_slab, ierr)
      !
      ! Create the memory dataspace.
      CALL H5SCREATE_SIMPLE_F(memrank, dims, memspace, ierr)
      !
      ! Select the hyperslab in memory.
      CALL H5SSELECT_HYPERSLAB_F(memspace, H5S_SELECT_SET_F, &
                                 offset_fill, dims_slab, ierr)
      !
      ! Read data from hyperslab in the file into the hyperslab in memory.
      CALL H5DREAD_F(dset_id, h5prec, buf, dims, ierr, memspace, dspace_id)
      !
      ! Close the memory space.
      CALL H5SCLOSE_F(memspace, ierr)
      !
      ! Close the dataspace for the dataset.
      CALL H5SCLOSE_F(dspace_id, ierr)
      !
      ! Close the dataset.
      CALL H5DCLOSE_F(dset_id, ierr)
   END SUBROUTINE READ_DATA_PARTIAL_REAL32_RANK2

   !> Procedure to read hyperslabs of data into contiguous buffer arrays.
   !!
   !> @param[in] hdf5_id HDF5 identifier to the object containing the dataset.
   !> @param[in] dims Dimensions of the buffer array.
   !> @param[in,out] buf Buffer array to read the data into.
   !> @param[in] dname Name of the dataset in the HDF5 file.
   !> @param[in] offset Offset for the hyperslab in the checkpoint file.
   SUBROUTINE READ_DATA_REAL32_RANK3(hdf5_id, dims, buf, dname, offset)
      ! Required modules.
      USE ISO_FORTRAN_ENV,ONLY: REAL32
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=HID_T),INTENT(IN) :: hdf5_id
      INTEGER(KIND=HSIZE_T),DIMENSION(3),INTENT(IN) :: dims
      REAL(KIND=REAL32),DIMENSION(dims(1),dims(2),dims(3)),INTENT(OUT) :: buf
      CHARACTER(LEN=*),INTENT(IN) :: dname
      INTEGER(KIND=HSIZE_T),DIMENSION(3),INTENT(IN) :: offset
      ! Local variables.
      ! HDF5 identifier for the variable type.
      INTEGER(KIND=HID_T) :: h5prec
      ! HDF5 identifiers for the dataspaces, memory spaces, and data sets.
      INTEGER(KIND=HID_T) :: dset_id, dspace_id, memspace
      ! Rank of the memory dataspace.
      INTEGER :: memrank = 3
      ! Error handling.
      INTEGER :: ierr
      !
      ! Identify the HDF5 type.
      h5prec = H5KIND_TO_TYPE(REAL32, H5_REAL_KIND)
      !
      ! Open the dataset.
      CALL H5DOPEN_F(hdf5_id, dname, dset_id, ierr)
      !
      ! Get the dataset's datspace identifier.
      CALL H5DGET_SPACE_F(dset_id, dspace_id, ierr)
      !
      ! Select the hyperslab in the dataset.
      CALL H5SSELECT_HYPERSLAB_F(dspace_id, H5S_SELECT_SET_F, offset, &
                                 dims, ierr)
      !
      ! Create the memory dataspace.
      CALL H5SCREATE_SIMPLE_F(memrank, dims, memspace, ierr)
      !
      ! Read data from hyperslab in the file into the hyperslab in memory.
      CALL H5DREAD_F(dset_id, h5prec, buf, dims, ierr, memspace, dspace_id)
      !
      ! Close the memory space.
      CALL H5SCLOSE_F(memspace, ierr)
      !
      ! Close the dataspace for the dataset.
      CALL H5SCLOSE_F(dspace_id, ierr)
      !
      ! Close the dataset.
      CALL H5DCLOSE_F(dset_id, ierr)
   END SUBROUTINE READ_DATA_REAL32_RANK3

   !> Procedure to read hyperslabs of data into hyperslabs of buffer arrays.
   !!
   !> @param[in] hdf5_id HDF5 identifier to the object containing the dataset.
   !> @param[in] dims Dimensions of the buffer array.
   !> @param[in,out] buf Buffer array to read the data into.
   !> @param[in] dname Name of the dataset in the HDF5 file.
   !> @param[in] dims_slab Dimensions of the hyperslab.
   !> @param[in] offset_read Offset for the hyperslab in the checkpoint file.
   !> @param[in] offset_fill Offset for the hyperslab in the buffer array.
   SUBROUTINE READ_DATA_PARTIAL_REAL32_RANK3(hdf5_id, dims, buf, dname, &
                                             dims_slab, offset_read, &
                                             offset_fill)
      ! Required modules.
      USE ISO_FORTRAN_ENV,ONLY: REAL32
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=HID_T),INTENT(IN) :: hdf5_id
      INTEGER(KIND=HSIZE_T),DIMENSION(3),INTENT(IN) :: dims
      REAL(KIND=REAL32),DIMENSION(dims(1),dims(2),dims(3)),INTENT(INOUT) :: buf
      CHARACTER(LEN=*),INTENT(IN) :: dname
      INTEGER(KIND=HSIZE_T),DIMENSION(3),INTENT(IN) :: dims_slab
      INTEGER(KIND=HSIZE_T),DIMENSION(3),INTENT(IN) :: offset_read, offset_fill
      ! Local variables.
      ! HDF5 identifier for the variable type.
      INTEGER(KIND=HID_T) :: h5prec
      ! HDF5 identifiers for the dataspaces, memory spaces, and data sets.
      INTEGER(KIND=HID_T) :: dset_id, dspace_id, memspace
      ! Rank of the memory dataspace.
      INTEGER :: memrank = 3
      ! Error handling.
      INTEGER :: ierr
      !
      ! Identify the HDF5 type.
      h5prec = H5KIND_TO_TYPE(REAL32, H5_REAL_KIND)
      !
      ! Open the dataset.
      CALL H5DOPEN_F(hdf5_id, dname, dset_id, ierr)
      !
      ! Get the dataset's datspace identifier.
      CALL H5DGET_SPACE_F(dset_id, dspace_id, ierr)
      !
      ! Select the hyperslab in the dataset.
      CALL H5SSELECT_HYPERSLAB_F(dspace_id, H5S_SELECT_SET_F, offset_read, &
                                 dims_slab, ierr)
      !
      ! Create the memory dataspace.
      CALL H5SCREATE_SIMPLE_F(memrank, dims, memspace, ierr)
      !
      ! Select the hyperslab in memory.
      CALL H5SSELECT_HYPERSLAB_F(memspace, H5S_SELECT_SET_F, &
                                 offset_fill, dims_slab, ierr)
      !
      ! Read data from hyperslab in the file into the hyperslab in memory.
      CALL H5DREAD_F(dset_id, h5prec, buf, dims, ierr, memspace, dspace_id)
      !
      ! Close the memory space.
      CALL H5SCLOSE_F(memspace, ierr)
      !
      ! Close the dataspace for the dataset.
      CALL H5SCLOSE_F(dspace_id, ierr)
      !
      ! Close the dataset.
      CALL H5DCLOSE_F(dset_id, ierr)
   END SUBROUTINE READ_DATA_PARTIAL_REAL32_RANK3

   !> Procedure to read hyperslabs of data into contiguous buffer arrays.
   !!
   !> @param[in] hdf5_id HDF5 identifier to the object containing the dataset.
   !> @param[in] dims Dimensions of the buffer array.
   !> @param[out] buf Buffer array to read the data into.
   !> @param[in] dname Name of the dataset in the HDF5 file.
   !> @param[in] offset Offset for the hyperslab in the checkpoint file.
   SUBROUTINE READ_DATA_REAL64_RANK1(hdf5_id, dims, buf, dname, offset)
      ! Required modules.
      USE ISO_FORTRAN_ENV,ONLY: REAL64
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=HID_T),INTENT(IN) :: hdf5_id
      INTEGER(KIND=HSIZE_T),DIMENSION(1),INTENT(IN) :: dims
      REAL(KIND=REAL64),DIMENSION(dims(1)),INTENT(OUT) :: buf
      CHARACTER(LEN=*),INTENT(IN) :: dname
      INTEGER(KIND=HSIZE_T),DIMENSION(1),INTENT(IN) :: offset
      ! Local variables.
      ! HDF5 identifier for the variable type.
      INTEGER(KIND=HID_T) :: h5prec
      ! HDF5 identifiers for the dataspaces and data sets.
      INTEGER(KIND=HID_T) :: dset_id, dspace_id, memspace
      ! Rank of the memory dataspace.
      INTEGER :: memrank = 1
      ! Error handling.
      INTEGER :: ierr
      !
      ! Identify the HDF5 type.
      h5prec = H5KIND_TO_TYPE(REAL64, H5_REAL_KIND)
      !
      ! Open the dataset.
      CALL H5DOPEN_F(hdf5_id, dname, dset_id, ierr)
      !
      ! Get the dataset's datspace identifier.
      CALL H5DGET_SPACE_F(dset_id, dspace_id, ierr)
      !
      ! Select the hyperslab in the dataset.
      CALL H5SSELECT_HYPERSLAB_F(dspace_id, H5S_SELECT_SET_F, offset, &
                                 dims, ierr)
      !
      ! Create the memory dataspace.
      CALL H5SCREATE_SIMPLE_F(memrank, dims, memspace, ierr)
      !
      ! Read data from hyperslab in the file into the hyperslab in memory.
      CALL H5DREAD_F(dset_id, h5prec, buf, dims, ierr, memspace, dspace_id)
      !
      ! Close the memory space.
      CALL H5SCLOSE_F(memspace, ierr)
      !
      ! Close the dataspace for the dataset.
      CALL H5SCLOSE_F(dspace_id, ierr)
      !
      ! Close the dataset.
      CALL H5DCLOSE_F(dset_id, ierr)
   END SUBROUTINE READ_DATA_REAL64_RANK1

   !> Procedure to read hyperslabs of data into contiguous buffer arrays.
   !!
   !> @param[in] hdf5_id HDF5 identifier to the object containing the dataset.
   !> @param[in] dims Dimensions of the buffer array.
   !> @param[out] buf Buffer array to read the data into.
   !> @param[in] dname Name of the dataset in the HDF5 file.
   !> @param[in] offset Offset for the hyperslab in the checkpoint file.
   SUBROUTINE READ_DATA_REAL64_RANK2(hdf5_id, dims, buf, dname, offset)
      ! Required modules.
      USE ISO_FORTRAN_ENV,ONLY: REAL64
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=HID_T),INTENT(IN) :: hdf5_id
      INTEGER(KIND=HSIZE_T),DIMENSION(2),INTENT(IN) :: dims
      REAL(KIND=REAL64),DIMENSION(dims(1),dims(2)),INTENT(OUT) :: buf
      CHARACTER(LEN=*),INTENT(IN) :: dname
      INTEGER(KIND=HSIZE_T),DIMENSION(2),INTENT(IN) :: offset
      ! Local variables.
      ! HDF5 identifier for the variable type.
      INTEGER(KIND=HID_T) :: h5prec
      ! HDF5 identifiers for the dataspaces and data sets.
      INTEGER(KIND=HID_T) :: dset_id, dspace_id, memspace
      ! Rank of the memory dataspace.
      INTEGER :: memrank = 2
      ! Error handling.
      INTEGER :: ierr
      !
      ! Identify the HDF5 type.
      h5prec = H5KIND_TO_TYPE(REAL64, H5_REAL_KIND)
      !
      ! Open the dataset.
      CALL H5DOPEN_F(hdf5_id, dname, dset_id, ierr)
      !
      ! Get the dataset's datspace identifier.
      CALL H5DGET_SPACE_F(dset_id, dspace_id, ierr)
      !
      ! Select the hyperslab in the dataset.
      CALL H5SSELECT_HYPERSLAB_F(dspace_id, H5S_SELECT_SET_F, offset, &
                                 dims, ierr)
      !
      ! Create the memory dataspace.
      CALL H5SCREATE_SIMPLE_F(memrank, dims, memspace, ierr)
      !
      ! Read data from hyperslab in the file into the hyperslab in memory.
      CALL H5DREAD_F(dset_id, h5prec, buf, dims, ierr, memspace, dspace_id)
      !
      ! Close the memory space.
      CALL H5SCLOSE_F(memspace, ierr)
      !
      ! Close the dataspace for the dataset.
      CALL H5SCLOSE_F(dspace_id, ierr)
      !
      ! Close the dataset.
      CALL H5DCLOSE_F(dset_id, ierr)
   END SUBROUTINE READ_DATA_REAL64_RANK2

   !> Procedure to read hyperslabs of data into hyperslabs of buffer arrays.
   !!
   !> @param[in] hdf5_id HDF5 identifier to the object containing the dataset.
   !> @param[in] dims Dimensions of the buffer array.
   !> @param[in,out] buf Buffer array to read the data into.
   !> @param[in] dname Name of the dataset in the HDF5 file.
   !> @param[in] dims_slab Dimensions of the hyperslab.
   !> @param[in] offset_read Offset for the hyperslab in the checkpoint file.
   !> @param[in] offset_fill Offset for the hyperslab in the buffer array.
   SUBROUTINE READ_DATA_PARTIAL_REAL64_RANK2(hdf5_id, dims, buf, dname, &
                                             dims_slab, offset_read, &
                                             offset_fill)
      ! Required modules.
      USE ISO_FORTRAN_ENV,ONLY: REAL64
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=HID_T),INTENT(IN) :: hdf5_id
      INTEGER(KIND=HSIZE_T),DIMENSION(2),INTENT(IN) :: dims
      REAL(KIND=REAL64),DIMENSION(dims(1),dims(2)),INTENT(INOUT) :: buf
      CHARACTER(LEN=*),INTENT(IN) :: dname
      INTEGER(KIND=HSIZE_T),DIMENSION(2),INTENT(IN) :: dims_slab
      INTEGER(KIND=HSIZE_T),DIMENSION(2),INTENT(IN) :: offset_read, offset_fill
      ! Local variables.
      ! HDF5 identifier for the variable type.
      INTEGER(KIND=HID_T) :: h5prec
      ! HDF5 identifiers for the dataspaces, memory spaces, and data sets.
      INTEGER(KIND=HID_T) :: dset_id, dspace_id, memspace
      ! Rank of the memory dataspace.
      INTEGER :: memrank = 2
      ! Error handling.
      INTEGER :: ierr
      !
      ! Identify the HDF5 type.
      h5prec = H5KIND_TO_TYPE(REAL64, H5_REAL_KIND)
      !
      ! Open the dataset.
      CALL H5DOPEN_F(hdf5_id, dname, dset_id, ierr)
      !
      ! Get the dataset's datspace identifier.
      CALL H5DGET_SPACE_F(dset_id, dspace_id, ierr)
      !
      ! Select the hyperslab in the dataset.
      CALL H5SSELECT_HYPERSLAB_F(dspace_id, H5S_SELECT_SET_F, offset_read, &
                                 dims_slab, ierr)
      !
      ! Create the memory dataspace.
      CALL H5SCREATE_SIMPLE_F(memrank, dims, memspace, ierr)
      !
      ! Select the hyperslab in memory.
      CALL H5SSELECT_HYPERSLAB_F(memspace, H5S_SELECT_SET_F, &
                                 offset_fill, dims_slab, ierr)
      !
      ! Read data from hyperslab in the file into the hyperslab in memory.
      CALL H5DREAD_F(dset_id, h5prec, buf, dims, ierr, memspace, dspace_id)
      !
      ! Close the memory space.
      CALL H5SCLOSE_F(memspace, ierr)
      !
      ! Close the dataspace for the dataset.
      CALL H5SCLOSE_F(dspace_id, ierr)
      !
      ! Close the dataset.
      CALL H5DCLOSE_F(dset_id, ierr)
   END SUBROUTINE READ_DATA_PARTIAL_REAL64_RANK2

   !> Procedure to read hyperslabs of data into contiguous buffer arrays.
   !!
   !> @param[in] hdf5_id HDF5 identifier to the object containing the dataset.
   !> @param[in] dims Dimensions of the buffer array.
   !> @param[in,out] buf Buffer array to read the data into.
   !> @param[in] dname Name of the dataset in the HDF5 file.
   !> @param[in] offset Offset for the hyperslab in the checkpoint file.
   SUBROUTINE READ_DATA_REAL64_RANK3(hdf5_id, dims, buf, dname, offset)
      ! Required modules.
      USE ISO_FORTRAN_ENV,ONLY: REAL64
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=HID_T),INTENT(IN) :: hdf5_id
      INTEGER(KIND=HSIZE_T),DIMENSION(3),INTENT(IN) :: dims
      REAL(KIND=REAL64),DIMENSION(dims(1),dims(2),dims(3)),INTENT(OUT) :: buf
      CHARACTER(LEN=*),INTENT(IN) :: dname
      INTEGER(KIND=HSIZE_T),DIMENSION(3),INTENT(IN) :: offset
      ! Local variables.
      ! HDF5 identifier for the variable type.
      INTEGER(KIND=HID_T) :: h5prec
      ! HDF5 identifiers for the dataspaces, memory spaces, and data sets.
      INTEGER(KIND=HID_T) :: dset_id, dspace_id, memspace
      ! Rank of the memory dataspace.
      INTEGER :: memrank = 3
      ! Error handling.
      INTEGER :: ierr
      !
      ! Identify the HDF5 type.
      h5prec = H5KIND_TO_TYPE(REAL64, H5_REAL_KIND)
      !
      ! Open the dataset.
      CALL H5DOPEN_F(hdf5_id, dname, dset_id, ierr)
      !
      ! Get the dataset's datspace identifier.
      CALL H5DGET_SPACE_F(dset_id, dspace_id, ierr)
      !
      ! Select the hyperslab in the dataset.
      CALL H5SSELECT_HYPERSLAB_F(dspace_id, H5S_SELECT_SET_F, offset, &
                                 dims, ierr)
      !
      ! Create the memory dataspace.
      CALL H5SCREATE_SIMPLE_F(memrank, dims, memspace, ierr)
      !
      ! Read data from hyperslab in the file into the hyperslab in memory.
      CALL H5DREAD_F(dset_id, h5prec, buf, dims, ierr, memspace, dspace_id)
      !
      ! Close the memory space.
      CALL H5SCLOSE_F(memspace, ierr)
      !
      ! Close the dataspace for the dataset.
      CALL H5SCLOSE_F(dspace_id, ierr)
      !
      ! Close the dataset.
      CALL H5DCLOSE_F(dset_id, ierr)
   END SUBROUTINE READ_DATA_REAL64_RANK3

   !> Procedure to read hyperslabs of data into hyperslabs of buffer arrays.
   !!
   !> @param[in] hdf5_id HDF5 identifier to the object containing the dataset.
   !> @param[in] dims Dimensions of the buffer array.
   !> @param[in,out] buf Buffer array to read the data into.
   !> @param[in] dname Name of the dataset in the HDF5 file.
   !> @param[in] dims_slab Dimensions of the hyperslab.
   !> @param[in] offset_read Offset for the hyperslab in the checkpoint file.
   !> @param[in] offset_fill Offset for the hyperslab in the buffer array.
   SUBROUTINE READ_DATA_PARTIAL_REAL64_RANK3(hdf5_id, dims, buf, dname, &
                                             dims_slab, offset_read, &
                                             offset_fill)
      ! Required modules.
      USE ISO_FORTRAN_ENV,ONLY: REAL64
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=HID_T),INTENT(IN) :: hdf5_id
      INTEGER(KIND=HSIZE_T),DIMENSION(3),INTENT(IN) :: dims
      REAL(KIND=REAL64),DIMENSION(dims(1),dims(2),dims(3)),INTENT(INOUT) :: buf
      CHARACTER(LEN=*),INTENT(IN) :: dname
      INTEGER(KIND=HSIZE_T),DIMENSION(3),INTENT(IN) :: dims_slab
      INTEGER(KIND=HSIZE_T),DIMENSION(3),INTENT(IN) :: offset_read, offset_fill
      ! Local variables.
      ! HDF5 identifier for the variable type.
      INTEGER(KIND=HID_T) :: h5prec
      ! HDF5 identifiers for the dataspaces, memory spaces, and data sets.
      INTEGER(KIND=HID_T) :: dset_id, dspace_id, memspace
      ! Rank of the memory dataspace.
      INTEGER :: memrank = 3
      ! Error handling.
      INTEGER :: ierr
      !
      ! Identify the HDF5 type.
      h5prec = H5KIND_TO_TYPE(REAL64, H5_REAL_KIND)
      !
      ! Open the dataset.
      CALL H5DOPEN_F(hdf5_id, dname, dset_id, ierr)
      !
      ! Get the dataset's datspace identifier.
      CALL H5DGET_SPACE_F(dset_id, dspace_id, ierr)
      !
      ! Select the hyperslab in the dataset.
      CALL H5SSELECT_HYPERSLAB_F(dspace_id, H5S_SELECT_SET_F, offset_read, &
                                 dims_slab, ierr)
      !
      ! Create the memory dataspace.
      CALL H5SCREATE_SIMPLE_F(memrank, dims, memspace, ierr)
      !
      ! Select the hyperslab in memory.
      CALL H5SSELECT_HYPERSLAB_F(memspace, H5S_SELECT_SET_F, &
                                 offset_fill, dims_slab, ierr)
      !
      ! Read data from hyperslab in the file into the hyperslab in memory.
      CALL H5DREAD_F(dset_id, h5prec, buf, dims, ierr, memspace, dspace_id)
      !
      ! Close the memory space.
      CALL H5SCLOSE_F(memspace, ierr)
      !
      ! Close the dataspace for the dataset.
      CALL H5SCLOSE_F(dspace_id, ierr)
      !
      ! Close the dataset.
      CALL H5DCLOSE_F(dset_id, ierr)
   END SUBROUTINE READ_DATA_PARTIAL_REAL64_RANK3

END MODULE IO
