      subroutine update_part (pxv, pind, npdim, nump, nvar)

#ifdef LAG
#ifdef LAGSC2
c
c
	use compart
	implicit none
c
        real*8 rtime1,rtime2,rtime3,rtime0
        real cpu,term1,term2
        real avcpu,cpumin,cpumax,cpu_other
	real*8 rtimea,rtimeb,rtimec
	real cpu_loop,cpu_comm,cpu_spcal
	save cpu_other,cpu_loop,cpu_comm,cpu_spcal

	integer npdim, nump
	real (p8) pxv (nppad*npdim/numtasks,nvar)
	integer pind(nppad*npdim/numtasks)

        integer jj,nps,tid,nid,ip,newp,nvar
	integer npmax, npmin
        integer ips,ipin(8),ipout(8),msgid(16),msgid2(16)
        integer mpist(MPI_STATUS_SIZE,16)
        integer mpist2(MPI_STATUS_SIZE,16)

	integer, allocatable :: nump_all(:)
	integer i,j,itemp

c
	if (nump.eq.0) return

        ips = 0
        ipin(:) = 0 ! no particles coming in
        ipout(:) = 0  ! no of particles going out 

        do ip=1,nump


! gives taskid (tid) and neighbour id (nid) for a particle 
! using its x- and z- coordinate
! nid can be 0 to 8, 0 being for itself
        call pptoid(pxv(ip,1),pxv(ip,3),tid,nid,1000.)

!        if(nid.lt.0 .or. nid.gt.8) then
!        write(6,*) 'nid exceeds limits, taskid=',taskid,tid,nid
!        stop
!        endif

! if the particle doesnt change taskid copy it to pxv
! using integer "ips" as the counter
! else copy it to ppout array and keep track of particles going out 
! in variable "ipout"
        if (nid.eq.0) then
        ips = ips+1
        pxv(ips,1:nvar) = pxv(ip,1:nvar)
        pind(ips) = pind(ip)
        else
!        pout(nvar*ipout(nid)+1:nvar*ipout(nid)+nvar,nid) = pxv(ip,0:nvar-1)
        pout(nvar*ipout(nid)+1:nvar*ipout(nid)+nvar,nid) = pxv(ip,1:nvar)
        pindout(ipout(nid)+1,nid) = pind(ip)
        ipout(nid) = ipout(nid) + 1
        endif

        enddo


	if(maxval(ipout).gt.nump*nppad) then
	write(6,*) 'increase array size in update part, ',taskid, maxval(ipout)
	endif

! first send how many particles are to be sent to each neighbour
        do jj=1,8
        call MPI_IRECV (ipin(jj),1,MPI_INTEGER,nbtask(jj),9-jj,
     &                  MPI_COMM_WORLD,msgid(jj),ierr)
	enddo
        do jj=1,8
        call MPI_ISEND (ipout(jj),1,MPI_INTEGER,nbtask(jj),jj,
     &                  MPI_COMM_WORLD,msgid(8+jj),ierr)
        enddo
        call MPI_WAITALL (16,msgid,mpist,ierr)

! then send in the particle positions (and indices)
! array "pin" receives the particles on each neighbour
        do jj=1,8
        call MPI_IRECV (pin(1,jj),nvar*ipin(jj)+1,pptype,nbtask(jj),9-jj,
     &                  MPI_COMM_WORLD,msgid(jj),ierr)
        call MPI_IRECV (pindin(1,jj),ipin(jj)+1,MPI_INTEGER,nbtask(jj),19-jj,
     &                  MPI_COMM_WORLD,msgid2(jj),ierr)
	enddo
        do jj=1,8
        call MPI_ISEND (pout(1,jj),nvar*ipout(jj)+1,pptype,nbtask(jj),jj,
     &                  MPI_COMM_WORLD,msgid(8+jj),ierr)
        call MPI_ISEND (pindout(1,jj),ipout(jj)+1,MPI_INTEGER,nbtask(jj),10+jj,
     &                  MPI_COMM_WORLD,msgid2(8+jj),ierr)
        enddo
        call MPI_WAITALL (16,msgid,mpist,ierr)
        call MPI_WAITALL (16,msgid2,mpist2,ierr)


! append the particle list received in pin to the pxv array
! and update the number of particles now present on each task
        nump = ips
        do jj=1,8
	if(ipin(jj).gt.0) then
        do ip=1,ipin(jj)
!        pxv(nump+ip,0:nvar-1) = pin(nvar*(ip-1)+1:nvar*ip,jj)
        pxv(nump+ip,1:nvar) = pin(nvar*(ip-1)+1:nvar*ip,jj)
        pind(nump+ip) = pindin(ip,jj)
        enddo
        nump = nump + ipin(jj)
	endif
        enddo

        nps = ips + sum(ipin)
	if(nps.ne.nump) then
	if(taskid.eq.0) write(6,*) 'error in update_part, nps,nump=',nps,nump
	stop
	endif
	
!        write(6,*) 'new nump =',nump,nps,taskid


        call MPI_ALLREDUCE (nump,newp,1,MPI_INTEGER,MPI_SUM,
     1                   MPI_COMM_WORLD,mpierr)

	if(newp.ne.npdim) then
	if(taskid.eq.0) write(6,*) 'error in update_part,newp,npdim=',newp,npdim
	stop
	endif

        call MPI_ALLREDUCE (nump,npmax,1,MPI_INTEGER,MPI_MAX,
     1                   MPI_COMM_WORLD,mpierr)
        call MPI_ALLREDUCE (nump,npmin,1,MPI_INTEGER,MPI_MIN,
     1                   MPI_COMM_WORLD,mpierr)


#ifdef DIAG
	allocate (nump_all(0:numtasks-1))
	call MPI_ALLGATHER (nump,1,MPI_INTEGER,nump_all,1,MPI_INTEGER,
     1                      MPI_COMM_WORLD,mpierr)

	do j=0,numtasks-2
	do i=0,numtasks-2
	if (nump_all(i+1).gt.nump_all(i)) then
	itemp=nump_all(i)
	nump_all(i)=nump_all(i+1)
	nump_all(i+1)=itemp
	end if
	end do
	end do
#endif
	

	if(taskid.eq.0) then
	if(npdim.eq.nop) then
#ifdef DIAG
	write (702,"('distribution of pp across MPI tasks, istep=',i6)") istep
	write (702,"(12i6)") nump_all
#endif
	write(6,*) 'update_part g1: istep=', istep, 'max, min, avg=', npmax, npmin, npdim/numtasks
	elseif(npdim.eq.nop2) then
	write(6,*) 'update_part g2: istep=', istep, 'max, min, avg=', npmax, npmin, npdim/numtasks
	endif
	endif

	if(npmax.gt.nppad*npdim/numtasks) then
	if(taskid.eq.0) then
	write(6,*) 'npmax exceeds array dimensions in update part'
	write(6,*) 'increase nppad value'
	endif
	stop
	endif

#ifdef DIAG
	deallocate (nump_all)
#endif

      return


#endif
#endif
      end

