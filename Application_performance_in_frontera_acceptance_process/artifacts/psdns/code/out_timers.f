	subroutine out_timers(acccpu2)

! writes file timers with communication and IO timers
	
	use com
	use timers_comm
	use timers_io
      implicit none
	real :: aux(3),aux2, acccpu2
#ifdef TIMERS_IO	  
	real :: aux3(3),aux4
#endif

#if defined TIMERS | defined TIMERS_IO
	if (taskid.eq.0) open (801,file='timers')
#endif

#ifdef TIMERS
	aux2=acccpu2
      call MPI_REDUCE (aux2,aux(1),1,mpireal,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE (aux2,aux(2),1,mpireal,MPI_MIN,0,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE (aux2,aux(3),1,mpireal,MPI_MAX,0,MPI_COMM_WORLD,ierr)
	aux(1)=aux(1)/numtasks

      if (taskid.eq.0) then
	  write (801,704) nx,numtasks,iproc,jproc, nu, nsteps
 704	  format ('nx=',i6,'  numtasks=',i6, ' ;  iproc, jproc=',i4, ' x ',i4,
     1          /'nu=',i3, '  repeat count=', i3)
        write (801,fmt="('total time       : min/ave/max=',3f11.3)") aux(2),aux(1),aux(3)
        write (801,fmt="('tot time / step  : min/ave/max=',3f11.3)") aux(2)/nsteps,aux(1)/nsteps,aux(3)/nsteps
	endif

	aux2=t_alltoall
      call MPI_REDUCE (aux2,aux(1),1,mpireal,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE (aux2,aux(2),1,mpireal,MPI_MIN,0,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE (aux2,aux(3),1,mpireal,MPI_MAX,0,MPI_COMM_WORLD,ierr)
	aux(1)=aux(1)/numtasks
      if (taskid.eq.0) then
        write (801,fmt="('alltoall time    : min/ave/max=',3f11.3)") aux(2),aux(1),aux(3)
        write (801,fmt="('alltoall / step  : min/ave/max=',3f11.3)") aux(2)/nsteps,aux(1)/nsteps,aux(3)/nsteps
	endif

	aux2=acccpu2- t_alltoall
      call MPI_REDUCE (aux2,aux(1),1,mpireal,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE (aux2,aux(2),1,mpireal,MPI_MIN,0,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE (aux2,aux(3),1,mpireal,MPI_MAX,0,MPI_COMM_WORLD,ierr)
	aux(1)=aux(1)/numtasks
      if (taskid.eq.0) then
        write (801,fmt="('tot-alltoall     : min/ave/max=',3f11.3)") aux(2),aux(1),aux(3)
        write (801,fmt="('tot-alltoall/step: min/ave/max=',3f11.3)") aux(2)/nsteps,aux(1)/nsteps,aux(3)/nsteps
	endif

	aux2=t_alltoall/acccpu2
      call MPI_REDUCE (aux2,aux(1),1,mpireal,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE (aux2,aux(2),1,mpireal,MPI_MIN,0,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE (aux2,aux(3),1,mpireal,MPI_MAX,0,MPI_COMM_WORLD,ierr)
	if (taskid.eq.0)  then 
	  aux=aux*100
        write (801,fmt="('alltoall/total % : min/ave/max=',3f11.0)") aux(2),aux(1)/numtasks,aux(3)
      end if
#endif	
#ifdef TIMERS_IO
	aux4=iread_io*(2*b8)/(1024.*1024.)
	call MPI_REDUCE(aux4,aux3(1),1,mpireal,MPI_SUM,0,MPI_COMM_WORLD,ierr)
	call MPI_REDUCE(aux4,aux3(2),1,mpireal,MPI_MIN,0,MPI_COMM_WORLD,ierr)
	call MPI_REDUCE(aux4,aux3(3),1,mpireal,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      if (taskid.eq.0) then
	  write (801,*) '----------- IO timers -------------'
	  write (801,705) nx,numtasks,iproc,jproc, nu
 705	  format ('nx=',i6,'  numtasks=',i6, ' ;  iproc, jproc=',i4, ' x ',i4,
     1          /'nu=',i3)
	  write (801,fmt="('Read, # of files    : ',i6)") abs(irz)
        write (801,fmt="('Read, Total GBytes  : ',1f11.3)") aux3(1)/1024.
	  aux3(1)=aux3(1)/numtasks
        write (801,fmt="('Read, MBytes     : min/ave/max=',3f11.3)") aux3(2),aux3(1),aux3(3)
	endif  

	aux4=tread_io
	call MPI_REDUCE(aux4,aux3(1),1,mpireal,MPI_SUM,0,MPI_COMM_WORLD,ierr)
	call MPI_REDUCE(aux4,aux3(2),1,mpireal,MPI_MIN,0,MPI_COMM_WORLD,ierr)
	call MPI_REDUCE(aux4,aux3(3),1,mpireal,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      if (taskid.eq.0) then
	  aux3(1)=aux3(1)/numtasks
        write (801,fmt="('Read, total time : min/ave/max=',3f11.3)") aux3(2),aux3(1),aux3(3)
	endif  

	aux4=iread_io*(2*b8)/(1024.*1024.)/tread_io
	call MPI_REDUCE(aux4,aux3(1),1,mpireal,MPI_SUM,0,MPI_COMM_WORLD,ierr)
	call MPI_REDUCE(aux4,aux3(2),1,mpireal,MPI_MIN,0,MPI_COMM_WORLD,ierr)
	call MPI_REDUCE(aux4,aux3(3),1,mpireal,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      if (taskid.eq.0) then
	  aux3(1)=aux3(1)/numtasks
        write (801,fmt="('Read, MB/secs    : min/ave/max=',3f11.3)") aux3(2),aux3(1),aux3(3)
	endif  

	aux4=iwrite_io*(2*b8)/(1024.*1024.)
	call MPI_REDUCE(aux4,aux3(1),1,mpireal,MPI_SUM,0,MPI_COMM_WORLD,ierr)
	call MPI_REDUCE(aux4,aux3(2),1,mpireal,MPI_MIN,0,MPI_COMM_WORLD,ierr)
	call MPI_REDUCE(aux4,aux3(3),1,mpireal,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      if (taskid.eq.0) then
	  write (801,fmt="('Write, # of files   : ',i6)") numtasks
        write (801,fmt="('Write, Total GBytes : ',1f11.3)") aux3(1)/1024.
	  aux3(1)=aux3(1)/numtasks
        write (801,fmt="('Write, MBytes    : min/ave/max=',3f11.3)") aux3(2),aux3(1),aux3(3)
	endif  
	aux4=twrite_io
	call MPI_REDUCE(aux4,aux3(1),1,mpireal,MPI_SUM,0,MPI_COMM_WORLD,ierr)
	call MPI_REDUCE(aux4,aux3(2),1,mpireal,MPI_MIN,0,MPI_COMM_WORLD,ierr)
	call MPI_REDUCE(aux4,aux3(3),1,mpireal,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      if (taskid.eq.0) then
	  aux3(1)=aux3(1)/numtasks
        write (801,fmt="('Write, total time: min/ave/max=',3f11.3)") aux3(2),aux3(1),aux3(3)
	endif  
	aux4=iwrite_io*(2*b8)/(1024.*1024.)/twrite_io
	call MPI_REDUCE(aux4,aux3(1),1,mpireal,MPI_SUM,0,MPI_COMM_WORLD,ierr)
	call MPI_REDUCE(aux4,aux3(2),1,mpireal,MPI_MIN,0,MPI_COMM_WORLD,ierr)
	call MPI_REDUCE(aux4,aux3(3),1,mpireal,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      if (taskid.eq.0) then
	  aux3(1)=aux3(1)/numtasks
        write (801,fmt="('Write, MB/secs   : min/ave/max=',3f11.3)") aux3(2),aux3(1),aux3(3)
	endif  

#endif
#if defined TIMERS | defined TIMERS_IO
	if (taskid.eq.0)  close (801)
#endif
	end subroutine out_timers
