      subroutine print_hwm
      use mpicom
      implicit none
      integer, parameter :: lun=99
      character*80 :: line
      integer(kind=8):: hwmrss(2)
      real (kind=8) :: hr(2),hrmax(2),hrmin(2),hrsum(2)
!      integer :: ierr,taskid,ntasks

      open(lun,file='/proc/self/status')
      hwmrss=0
      do while(.true.)
        read(lun,'(a)',end=99) line
        if (line(1:6).eq.'VmHWM:') read(line(8:80),*) hwmrss(1)
        if (line(1:6).eq.'VmRSS:') read(line(8:80),*) hwmrss(2)
      enddo
  99  close(lun)
!      print *,'HWM = ',hwmrss(1),' KiB'
!      print *,'RSS = ',hwnrss(2),' KiB'

      hr = hwmrss/1024.
      call MPI_Reduce(hr,hrmax,2,MPI_REAL8,MPI_MAX,0,
     & MPI_COMM_WORLD,ierr)
      call MPI_Reduce(hr,hrmin,2,MPI_REAL8,MPI_MIN,0,
     & MPI_COMM_WORLD,ierr)
      call MPI_Reduce(hr,hrsum,2,MPI_REAL8,MPI_SUM,0,
     & MPI_COMM_WORLD,ierr)

      if (taskid.eq.0) then
        print *,'Memory high water mark (min/ave/max) = ',
     &  hrmin(1),hrsum(1)/numtasks,hrmax(1),' MiB'
        print *,'Resident Set size (min/ave/max) = ',
     &  hrmin(2),hrsum(2)/numtasks,hrmax(2),' MiB'
      endif
c
	return
	end
