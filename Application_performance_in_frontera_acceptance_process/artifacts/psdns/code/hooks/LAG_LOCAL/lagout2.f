       subroutine lagout2 (rtime0)
c
c routine to control output of particle properties,
c called from sub. partic
c
c
#ifdef LAG
c
	use compart
	use lag_timers
	implicit none
c
      integer dmy(3),hms(3),date_time(8)
c      common/timing/acccpu,totcpu,rwall0
	integer iendwt

	real totcpu
	integer, allocatable :: allmiseq(:)
	real(8) rtime0
c
	if (taskid.eq.0) write (6,*) 'enter lagout2'
c

      if (ngp.eq.1) then
      call wrtpop1 (0)
      else 
	if (shear.eq.0..and.itetra.eq.0) then
#ifdef NOTYET
	call wrtpop2 (0)
#endif
	else
!	if (nsubset.eq.1) call wrtpop3 (0)
!	if (nsubset.eq.1) call pop_file3 (0)
!	if (nsubset.gt.1) call wrtpop4 (0)
!	if (nsubset.gt.1) call pop_file4 (0)
!	call pop_file4 (0)
	call pop1_file_local (0)
	end if
      end if
c
#ifdef LAGSC2
 	if (nsubset.gt.1) call pop2_file5(0)
#endif

c
	cpu_partic2(2)=cpu_partic2(2)+MPI_WTIME()-rtime0
c
	call write_timings
	if (taskid.eq.0) write (6,*) 'after write_timings'
	call lag_timings
	call lag_iotimings
	if (taskid.eq.0) write (6,*) 'after lag_timings'

      call MPI_REDUCE (acccpu,totcpu,1,mpireal,MPI_SUM,0,
     1                 MPI_COMM_WORLD,mpierr)
c
      call date_and_time (values=date_time)
      dmy(1:3)=date_time(3:1:-1)
      hms(1:3)=date_time(5:7)
        if (taskid.eq.0.or.taskid.eq.numtasks-1) then
      write (6,605) taskid,dmy(2),dmy(1),dmy(3),hms
 605  format (' taskid=',i7, ' DONE:',' date & time is  ',i2,'/',i2,
     1        '/',i4,2x,i2,':',i2,':',i2)
        end if

      if (taskid.eq.0) then
      call rdate (dmy,hms,iendwt)
      write (6,604) numtasks,totcpu
        rwall0 = MPI_WTIME() - rwall0
	open (7,file='log',position='append')
      write (7,604) numtasks,totcpu,totcpu/3600.
      write (7,703) rwall0
 604  format ('total cpu over ',i6,' nodes=',f12.3,' secs.',
     1       f11.2, ' processor-hours')
 703  format ('overall wall time for task 0 =',f7.0, ' secs.')
      end if
c
c
      call MPI_FINALIZE (mpierr)
      stop
c
#endif
      return
      end
