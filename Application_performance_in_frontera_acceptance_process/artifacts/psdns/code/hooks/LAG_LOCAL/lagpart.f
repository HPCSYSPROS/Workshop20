	subroutine lagpart
c
          use comp
	use compart, only: nop,nop2,iout
	use lag_timers

          implicit none
          include 'intvars'

	integer icall,jj
	data icall/0/
	save icall

	real(8) rtime1,rtime2,rtime3,rtime4
c
	icall=icall+1
	
	if (jstep.gt.1.and.nop.eq.0.and.nop2.eq.0) return
c
	if (icall.eq.1) then
	cpu_lag = 0.
	cpu_lagout = 0.
	cpu_spxyz(:,:)=0.
	cpu_partic=0.
	cpu_partic2(:)=0.
	cpu_partsp(:,:)=0.
	cpu_intbf(:,:,:)=0.
	cpu_spcal(:,:,:)=0.
	cpu_rkpart(:,:)=0.
	itrack=0
	icpvrms=0
	nspxyz=0
	nspcal=0
	nintbf=0
	nspcal_io=0
	nintbf_io=0
	end if


c RK2
	
	rtime1=MPI_WTIME()
c
	if (rkmethod.eq.2) then
c
	if (kstep.eq.1) then
	if (jstep.eq.1) then
	call partic (4)
	else  
	rtime3 = MPI_WTIME()
	call partic (2)
	jj=1
	if (iout.gt.0) jj=2
	cpu_partic2(jj) = cpu_partic2(jj) + MPI_WTIME() - rtime3
c	if (taskid.eq.0) write (6,*) 'partic:',istep,jj,cpu_partsp(0,jj),cpu_partic2(jj)
	end if
	end if
c
c RK4
c
	else if (rkmethod.eq.4) then
c
	call partic_rk4
	
	end if

	rtime2 = MPI_WTIME()
	cpu_partic = cpu_partic + rtime2 - rtime1
	cpu_lag = cpu_lag + rtime2 - rtime1
c
	return
	end
