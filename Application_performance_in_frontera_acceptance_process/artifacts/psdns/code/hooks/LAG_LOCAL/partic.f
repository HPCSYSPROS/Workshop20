	subroutine partic (icall)
c
#ifdef LAG
c
	use compart
	use lag_timers
	implicit none
	include 'intvars'
	real*8 rtime1
	integer jj,ii,i
	integer icall
c
	real(8) rtime0


c
	if (nop.eq.0.and.nop2.eq.0) return 
c
	if (lpfunc(1).eq.0) return
c
	if (taskid.eq.0)
     1	write (6,*) 'enter partic: istep,kstep,icall=',istep,kstep,icall

c
	rtime0=MPI_WTIME()
	rtime1=MPI_WTIME()
c

 	call partsp (icall)
c
        call time_stamp ('partic: after partsp')
	jj=1
	if (iout.gt.0) jj=2
	if(icall.eq.2) then
	cpu_partsp(0,jj)=cpu_partsp(0,jj) + MPI_WTIME() - rtime1
	endif
c
c
	rtime1=MPI_WTIME()

#ifdef DEBUG
	write(2000+taskid,*) 'istep=',istep
	do i=1,nump2
	write(2000+taskid,*) pp2(i,0:3)
	enddo
#endif

! 	if (iout.eq.1) call lagout1
! 	if (iout.eq.2) call lagout2 (rtime0)
	call lagout_local (iout, rtime0)

	cpu_lagout = cpu_lagout + MPI_WTIME() - rtime1
c
!	if (taskid.eq.0)
!     1	write (6,*) ' exit partic: istep,kstep,icall=',istep,kstep,icall
c
#endif
c
	return
	end
