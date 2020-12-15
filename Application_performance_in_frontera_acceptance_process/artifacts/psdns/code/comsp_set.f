      subroutine comsp_set
c
      use comsp
c
      implicit none
#ifdef MODEL_SPECTRUM
      LOGICAL :: exs
      INTEGER :: exsInt, stat, mxyzIn
#endif

#ifdef SHELL_DK
      real(b8) beta_min
      beta_min=min(beta1,beta2,beta3)
      mxyz=max(nx*beta1,ny*beta2,nz*beta3)/beta_min
      if (taskid.eq.0) write (6,*) 'comsp_set: mxyz=',mxyz
#else
#ifndef MODEL_SPECTRUM
      mxyz = max(nx,ny,nz)
#else
      real(b8) beta_min
      beta_min=min(beta1,beta2,beta3)
      mxyz=max(nx*beta1,ny*beta2,nz*beta3)/beta_min
      if (taskid.eq.0) write (6,*) 'comsp_set: mxyz=',mxyz
#endif
#endif

#ifdef MODEL_SPECTRUM
      ! The root process checks if the file "mxyz" is present. If it is, it
      ! attempts to read the specified value of mxyz, and broadcasts it to all
      ! other processes.
      IF (taskid .EQ. 0) THEN
         INQUIRE(FILE='mxyz',EXIST=exs)
         IF (exs) THEN
            ! Read the data.
            OPEN(UNIT=1115,FILE='mxyz',STATUS='OLD',ACTION='READ')
            READ(1115,*,IOSTAT=stat) mxyzIn
            CLOSE(UNIT=1115)
            !
            ! The parameter exsInt is zero if the file exists.
            exsInt = 0
         ELSE
            ! Set exsInt to a positive value to indicate it is not present.
            exsInt = 1
         END IF
      END IF
      !
      ! Broadcast the logical checks to see how to continue.
      CALL MPI_BCAST(exsInt, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST(stat, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      !
      ! Only proceed if the file existed in the first place.
      IF (exsInt .EQ. 0) THEN
         ! Check if there was an IO error, and if there was, use the default
         ! value of mxyz set in the above code.
         IF (stat .NE. 0) THEN
            IF (taskid .EQ. 0) THEN
               WRITE(6,*) 'Error reading mxyz file. Using default value.'
            END IF
         ELSE
            ! Broadcast the value of mxyz from the file and set the local value
            ! of mxyz to that value.
            CALL MPI_BCAST(mxyzIn, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
            mxyz = mxyzIn
            IF (taskid .EQ. 0) THEN
               WRITE(6,*) 'Overwriting mxyz to:', mxyz
            END IF
         END IF
      END IF
#endif


      allocate (lijk(ncp,3))
      allocate (tmij(nc+3,3))
      allocate (taylor(nc+3,nc+3,3))

	allocate (laak(nc+3,3),taak(nc+3,3),raak(3))
	allocate (tmre(3),rms(nc+3),ett(3))
	allocate (sk(mxyz),corr(ncp))

	allocate (ek(mxyz))
	allocate (dk(mxyz))
c
	allocate (kx2(nxh),ky2(ny),kz2(nz))
c
      allocate (tfact(nxh))
      tfact(:)=2.
      tfact(1)=1.
c
c
c change from nxh to mxyz in first dimension of these arrays
c is necessary for domains where ny or nz may be larger than nx

      if (iovor.ge.1) allocate (vijky(mxyz,6,0:num_thr-1),
     1                          vijk(mxyz,6),
     2                          vk(mxyz))

c
	allocate (umean(nu),varce(nu),skcof(nu),flat(nu))
        allocate (imean(nu))
c	
!RAF	allocate (kt(3+nc))
	allocate (kt(max(3+nc,4)))
c
#ifndef NOSCALAR
	ncgd=nc*(nc+1)/2
	allocate (scdiss(nc))
	allocate (scgmsq(nc,3))
	allocate (scgcor(ncgd,3),scgcov(ncgd,3),scgvar(ncgd,3))
#endif
c
      return
      end
