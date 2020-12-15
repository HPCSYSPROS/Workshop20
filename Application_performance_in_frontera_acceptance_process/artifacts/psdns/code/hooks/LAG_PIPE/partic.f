	subroutine partic (icall)
c
#ifdef LAG
c
	use compart
	use lag_timers
	implicit none
	include 'intvars'
	real*8 rtime1
	integer jj,ii
	integer icall
c
	real(8) rtime0
c
	if (nop.eq.0.and.nop2.eq.0) return 
c
	if (lpfunc(1).eq.0) return
c
c	if (taskid.eq.0)
c     1	write (6,*) 'enter partic: istep,kstep,icall=',istep,kstep,icall

c        allocate (bs(bxisize,nby,bzjsize,nbsdim))
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
c 	deallocate (bs)
c
	rtime1=MPI_WTIME()
 	if (iout.eq.1) call lagout1
 	if (iout.eq.2) call lagout2 (rtime0)
	cpu_lagout = cpu_lagout + MPI_WTIME() - rtime1
c
c	if (taskid.eq.0)
c     1	write (6,*) ' exit partic: istep,kstep,icall=',istep,kstep,icall
c
#endif
c
	return
	end
