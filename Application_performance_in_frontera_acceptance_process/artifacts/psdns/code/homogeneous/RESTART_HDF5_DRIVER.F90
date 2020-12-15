!> @file RESTART_HDF5_DRIVER.F90
!> @author Matthew Clay
!> @brief Driver subroutine to read HDF5 checkpoint data.
!!
!! This is a subroutine wrapper for the RESTART_HDF5.F90 subroutine. Here we add
!! a relay system to improve IO performance for large process counts, and call
!! the necessary subroutine to transpose to the cylinder data layout. We also
!! print out some basic timers to assess IO perfomance.
!!
!> @param[in] path Path to the checkpoint file.
!> @param[in] cnum Number of the checkpoint to read.
!> @param[in,out] buff Array to hold the velocity/scalar fields.
SUBROUTINE RESTART_HDF5_DRIVER(path,cnum,uny)
   ! Required modules.
   USE ISO_FORTRAN_ENV,ONLY: REAL64
   USE MPI
   USE com,ONLY: B8,nc,kinit
   USE param,ONLY: ny,nz
   USE mpicom,ONLY: numtasks,taskid,num_al_i,zjsz,yjsz,xisz
   IMPLICIT NONE
   ! Calling arguments.
   CHARACTER(LEN=*),INTENT(IN) :: path
   INTEGER,INTENT(IN) :: cnum
   COMPLEX(KIND=B8),DIMENSION(ny,zjsz,xisz,3+nc),INTENT(INOUT) :: uny
   ! Local variables.
   ! Used to determine number of processes to read at once.
   INTEGER :: mread,nread,iread
   ! Signal sent from one process to the next in the relay.
   INTEGER :: rid
   ! MPI status for RECV.
   INTEGER,DIMENSION(MPI_STATUS_SIZE) :: stat
   ! Timers for IO performance.
   REAL(KIND=REAL64) :: tovr,tprc,tmin,tmax,tavg
   ! Buffer to redistribute the data to cylinder layout.
   COMPLEX(KIND=B8),DIMENSION(:,:,:),ALLOCATABLE :: buf
   ! Looping index for cylinder redistribution.
   INTEGER :: i
   ! MPI error handling.
   INTEGER :: mpierr
   !
   ! Overall timer on the root process.
   tovr = MPI_WTIME()
   !
   ! Determine the number of processes that will read at once.
   IF (numtasks .LE. 4096) THEN
      mread = numtasks
   ELSE IF (numtasks .LT. 262144) THEN
      mread = 4096
   ELSE
      mread = 8192
   END IF
   !
   ! Determine number of batch reads and identify the starting processes.
   nread = numtasks/mread
   iread = taskid/mread
   IF (taskid .EQ. 0) THEN
      WRITE(*,50) 'READING WITH RELAY SYSTEM. MREAD/NREAD/IREAD=', &
                  mread,nread,iread
   END IF
   !
   IF (iread .EQ. 0) THEN
      tprc = MPI_WTIME()
      CALL RESTART_HDF5(path,cnum,uny)
      tprc = MPI_WTIME() - tprc
   ELSE
      ! Wait for a signal from other process that has finished reading data
      ! before beginning to read data.
      CALL MPI_RECV(rid,1,MPI_INTEGER,taskid-mread,taskid-mread, &
                    MPI_COMM_WORLD,stat,mpierr)
      !
      tprc = MPI_WTIME()
      CALL RESTART_HDF5(path,cnum,uny)
      tprc = MPI_WTIME() - tprc
   END IF
   !
   ! When a process finishes reading data send signal to next process to
   ! start reading data.
   IF (iread .LT. (nread-1)) THEN
      CALL MPI_SSEND(taskid,1,MPI_INTEGER,taskid+mread,taskid, &
                     MPI_COMM_WORLD,mpierr)
   END IF
   !
   ! Redistribute to cylinder layout.
   ALLOCATE(buf(nz,xisz,yjsz))
   iloop: DO i = 1,3+nc
      IF (kinit(i) .LE. 0) CYCLE iloop
      buf(:,:,:) = (0.0_B8,0.0_B8)
      CALL kxcomm1_pen(uny(1,1,1,i),buf,1)
      IF (num_al_i(0) .GT. 0) THEN
         CALL xkcomm2_sq2cyl(buf,uny(1,1,1,i),1)
      END IF
   END DO iloop
   DEALLOCATE(buf)
   !
   ! Reduce timers.
   tovr = MPI_WTIME() - tovr
   CALL MPI_REDUCE(tprc,tmin,1,MPI_REAL8,MPI_MIN,0,MPI_COMM_WORLD,mpierr)
   CALL MPI_REDUCE(tprc,tmax,1,MPI_REAL8,MPI_MAX,0,MPI_COMM_WORLD,mpierr)
   CALL MPI_REDUCE(tprc,tavg,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpierr)
   !
   ! Report timers.
   IF (taskid .EQ. 0) THEN
      tavg = tavg/REAL(numtasks,REAL64)
      WRITE(*,70) 'HDF5 IO TIMINGS. OALL/TMIN/TAVG/TMAX:',tovr,tmin,tavg,tmax
   END IF
   !
   ! IO formats.
   50 FORMAT (A,1X,3I7)
   60 FORMAT (A,1X,I6)
   70 FORMAT (A,4ES12.5)
END SUBROUTINE RESTART_HDF5_DRIVER
