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
	integer i,igm
	real sum, sum0, pvrms, pvrms2, mvrms, ratio
c
c extra call to rk_part(0) removed, 9/1/2012
c 
c check rms particle velocity
c
	if (ioflag.eq.1) then
c
      sum=0.
      do 390 i=1,nop/numtasks
      sum=sum+pp(i,4)**2
 390  continue
        call MPI_ALLREDUCE (sum,sum0,1,mpireal,MPI_SUM,
     1                   MPI_COMM_WORLD,mpierr)
	sum=sum0
      pvrms=sqrt(sum/nop)
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
      do 392 i=1,nop2/numtasks
      sum=sum+pp2(4,i)**2
 392  continue
        call MPI_ALLREDUCE (sum,sum0,1,mpireal,MPI_SUM,
     1                   MPI_COMM_WORLD,mpierr)
	sum=sum0
      pvrms2=sqrt(sum/nop2)
      if (taskid.eq.0) write (6,602) istep-1,rms(1),pvrms,pvrms2
 602  format ('partsp: istep-1, E/L urms=',i6,1p,3e14.6)
#endif
c

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
