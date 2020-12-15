       subroutine lagout_local (outid, rtime0)
c
c routine to control output of particle properties,
c called from sub. partic
c
c
#ifdef LAG
c
	use compart
	use lag_timers
	integer ic, outid

        integer dmy(3),hms(3),date_time(8)
c       common/timing/acccpu,totcpu,rwall0
	integer iendwt
	real totcpu
	integer, allocatable :: allmiseq(:)
	real(8) rtime0

        integer iread, mread, nread, relay, jsource
        integer mpistatus (MPI_STATUS_SIZE)

c
! write particle properties on file, including velocities if desired
! vel. gradients also (controlled by LGRAD)


	if(outid.ne.1.and.outid.ne.2) return

	if(taskid.eq.0) write(6,*) 'lagout_local, iout=', iout

	if(outid.eq.1) then

	call MPI_BARRIER (MPI_COMM_WORLD, mpierr)
	call time_stamp ('hdf5_io start')

#ifdef NOV15
#else
       if (jproc.le.1024) mread=jproc
       if (mod(jproc,1024).eq.0) mread=1024

        nread=jproc/mread
        iread=jpid/mread

        if (iread.eq.0) then

        relay=taskid

#ifdef LGRAD
	if (lpgrad(1).eq.1) then
	if(nop.gt.0) call lagout_comm_h5 (pp(1,1), pp1ind, 14, ndpart, nump1, 1)
	if(nop2.gt.0) call lagout_comm_h5 (pp2(1,1), pp2ind, 14, ndpart2, nump2, 2)
	else
	if(nop.gt.0) call lagout_comm_h5 (pp(1,1), pp1ind, 6, ndpart, nump1, 1)
	if(nop2.gt.0) call lagout_comm_h5 (pp2(1,1), pp2ind, 6, ndpart2, nump2, 2)
	endif
#endif

        else

        jsource=taskid-mread*iproc

        call MPI_RECV (relay,1,MPI_INTEGER,jsource,jsource,
     1                 MPI_COMM_WORLD,mpistatus,mpierr)

#ifdef LGRAD
	if (lpgrad(1).eq.1) then
	if(nop.gt.0) call lagout_comm_h5 (pp(1,1), pp1ind, 14, ndpart, nump1, 1)
	if(nop2.gt.0) call lagout_comm_h5 (pp2(1,1), pp2ind, 14, ndpart2, nump2, 2)
	else
	if(nop.gt.0) call lagout_comm_h5 (pp(1,1), pp1ind, 6, ndpart, nump1, 1)
	if(nop2.gt.0) call lagout_comm_h5 (pp2(1,1), pp2ind, 6, ndpart2, nump2, 2)
	endif
#endif

        end if

        if (iread.lt.nread-1) then
        call MPI_SSEND (relay,1,MPI_INTEGER,taskid+mread*iproc,taskid,
     1                MPI_COMM_WORLD,mpierr)
        end if
c

	call MPI_BARRIER (MPI_COMM_WORLD, mpierr)
	call time_stamp ('hdf5_io end')

#endif
!NOV15 directive


	elseif (outid.eq.2) then


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

	endif ! if(outid..)
c
#endif

      return
      end
