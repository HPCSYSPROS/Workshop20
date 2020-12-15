       subroutine check_pvrms
c
c more general particle tracking call, for a given value of kstep
c
#ifdef LAG
c
	use compart
	implicit none
	include 'intvars'
c
	integer i,igm, npt
	real sum, sum0, pvrms, pvrms2, mvrms, ratio
	real pvmax, pvmin
	real pvmaxe, pvmine, pvmax1, pvmin1, pvmax2, pvmin2
c
c extra call to rk_part(0) removed, 9/1/2012
c 
c check rms particle velocity
c



	if (ioflag.eq.1) then
c
      sum=0.
#ifdef NOV15
#else
	pvrms=0.
	if (nop.gt.0) then
!      do 390 i=1,nop/numtasks
      do 390 i=1,nump1
      sum=sum+pp(i,4)**2
 390  continue
        call MPI_ALLREDUCE (sum,sum0,1,mpireal,MPI_SUM,
     1                   MPI_COMM_WORLD,mpierr)
	sum=sum0
      pvrms=sqrt(sum/nop)
	endif
#endif
c
#ifdef MOL
	if (nom.gt.0) then
      sum=0.
        do igm=1,ngm
      do 391 i=1,nom/numtasks/ngm
      sum=sum+mpp(i,4,igm)**2
 391  continue
        end do
        call MPI_ALLREDUCE (sum,sum0,1,mpireal,MPI_SUM,
     1                   MPI_COMM_WORLD,mpierr)
	sum=sum0
      mvrms=sqrt(sum/nom)
      if(taskid.eq.0) then 
        write (6,"('partsp: istep-1, E/L/M urms =',i6,1p,4e14.6)")
     1             istep-1,rms(1),pvrms,mvrms
        end if
	end if
#else
#ifndef LAGSC2
      if (taskid.eq.0) write (6,601) istep-1,rms(1),pvrms
 601  format ('partsp: istep-1, Eul/Lag/rms x-vel.=',i6,1p,2e14.6)
#endif
#endif
c
#ifdef LAGSC2
      sum=0.
#ifdef NOV15
#else
	pvrms2=0.
	if (nop2.gt.0) then
!      do 392 i=1,nop2/numtasks
      do 392 i=1,nump2
      sum=sum+pp2(i,4)**2
 392  continue
        call MPI_ALLREDUCE (sum,sum0,1,mpireal,MPI_SUM,
     1                   MPI_COMM_WORLD,mpierr)
	sum=sum0
      pvrms2=sqrt(sum/nop2)
	endif
#endif
      if (taskid.eq.0) write (6,602) istep-1,rms(1),pvrms,pvrms2
 602  format ('partsp: istep-1, E/L urms=',i6,1p,3e14.6)
#endif
c
	pvmax = maxval (u(1:nx,1:zisz,1:yjsz,1:3))
	pvmin = minval (u(1:nx,1:zisz,1:yjsz,1:3))
        call MPI_ALLREDUCE (pvmax, pvmaxe, 1, mpireal, MPI_MAX,
     1                   MPI_COMM_WORLD,mpierr)
        call MPI_ALLREDUCE (pvmin, pvmine, 1, mpireal, MPI_MIN,
     1                   MPI_COMM_WORLD,mpierr)




	if (nop.gt.0) then
!	npt = nop/numtasks
	pvmax = maxval(pp(1:nump1,4:6))
	pvmin = minval(pp(1:nump1,4:6))
        call MPI_ALLREDUCE (pvmax, pvmax1, 1, mpireal, MPI_MAX,
     1                   MPI_COMM_WORLD,mpierr)
        call MPI_ALLREDUCE (pvmin, pvmin1, 1, mpireal, MPI_MIN,
     1                   MPI_COMM_WORLD,mpierr)
	endif

	if (nop2.gt.0) then
	pvmax = maxval(pp2(1:nump2,4:6))
	pvmin = minval(pp2(1:nump2,4:6))
        call MPI_ALLREDUCE (pvmax, pvmax2, 1, mpireal, MPI_MAX,
     1                   MPI_COMM_WORLD,mpierr)
        call MPI_ALLREDUCE (pvmin, pvmin2, 1, mpireal, MPI_MIN,
     1                   MPI_COMM_WORLD,mpierr)
	endif

      if (taskid.eq.0) write (6,603) istep-1,pvmaxe,pvmax1,pvmax2
 603  format ('partsp: istep-1, E/L umax=',i6,1p,3e14.6)
      if (taskid.eq.0) write (6,604) istep-1,pvmine,pvmin1,pvmin2
 604  format ('partsp: istep-1, E/L umin=',i6,1p,3e14.6)

#ifdef TEST
      ratio=pvrms/rms(1)
      if (ratio.lt.0.95.or.ratio.gt.1.05) then
      write (6,*) 'partsp: Eul/Lag rms differing by more than 5%'
      write (6,*) 'stop and investigate:'
      write (6,*) 'particle velocities dumped onto fort.99'
      do i=1,nop
      write (99,609) i,pp(i,4)
      end do
        call MPI_ABORT (MPI_COMM_WORLD,ierror)
      end if
#endif

	end if
#endif
c
      return
      end
