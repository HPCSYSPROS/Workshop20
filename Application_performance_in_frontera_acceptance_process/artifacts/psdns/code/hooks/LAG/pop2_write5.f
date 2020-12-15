      subroutine pop2_write5 (k,j1,j2)
c
c alternative to pop2_write4, 2/1/09
c 
c output action corresponding to wrt2pop4
c
c let every processor handle some of the output, 7/4/08
c (desirable if no. of particles and no. of processors are both very large)
c
c output of particle properties
c
c this routine is used if there is only 1 "group" of particles
c (otherwise, use wrtpop2 or wrtpop3)
c
c if iopt=0, simply close the files without writing new data
c
#ifdef LAG
#ifdef LAGSC2
c
	use compart
	implicit none
c
      character*15 fn
      character*6 name,numer
      character*11 unform
      logical secvel
      data unform/'unformatted'/
      character*110 string
	character*30 caux
c
	real(p8), allocatable :: pg2(:)
c
 	integer j,j1,j2,k,mp
c
	integer ipfst,iplst,nl,i,i1,i2,il,mm,npr,lu
c
	real*8 rtime1,rtime2,rtime0
	real cpumpi,cpuio,cpu
c
cccc 
	integer iraddr,itask,jtask,irtask1,irtask2,dest,msgids,imt1
        integer, allocatable :: msgidr(:),msgtype(:,:)
        integer, allocatable :: mpistr(:,:),mpists(:)
c
	integer isubset
c
	rtime0 = MPI_WTIME()

	allocate (pg2(ndpart2/nsubset))
c
c first and last particles in this group
c
	if (taskid.eq.0.and.jstep.eq.1) write (6,*) 'enter pop2_write5, istep,j1,j2=',istep,j1,j2
c

	isubset=(taskid+1)/(numtasks/nsubset)+1
	mm=nop2/nsubset
	ipfst=isubset*mm+1
	iplst=ipfst+mm-1
c
c npr = no. of particles per record
c
      npr=min0(8192,mm)
      nl=nop2/nsubset/npr
      if (mod(nop2/nsubset,npr).ne.0) nl=nl+1
c
	mp=nop2/numtasks
c
	cpumpi=0.
	cpuio=0.
c
        allocate (msgidr(numtasks/nsubset))
        allocate (msgtype(nsubset,numtasks/nsubset-1))
        allocate (mpists(MPI_STATUS_SIZE))
        allocate (mpistr(MPI_STATUS_SIZE,numtasks/nsubset))

c
	do 100 j=j1,j2
c
	rtime1=MPI_WTIME()
c
        if (mod(taskid,numtasks/nsubset).eq.0) then
c
c post the receives
c taskids 0, nsubset, 2*nsubset, 3*nsubset, etc to receive messages
c from the next numtasks/nsubset - 1 tasks
c
c
        irtask1=taskid+1
        irtask2=irtask1+numtasks/nsubset-2
c
        do itask=irtask1,irtask2
        jtask=itask-irtask1+1
        imt1=itask/(numtasks/nsubset)+1
        msgtype(imt1,jtask)=1000*imt1+jtask
        iraddr=jtask*mp+1
        call MPI_IRECV (pg2(iraddr),mp,pptype,itask,
     1                  msgtype(imt1,jtask),
     1                  MPI_COMM_WORLD,msgidr(jtask),mpierr)
        end do
c
	pg2(1:mp)=pp2(1:mp,j)
        call MPI_WAITALL (numtasks/nsubset-1, msgidr,mpistr,mpierr)
c
        else
c
c post the sends
c
        imt1=taskid/(numtasks/nsubset)+1
        dest=taskid/(numtasks/nsubset)*(numtasks/nsubset)
        jtask=taskid-dest
        msgtype(imt1,jtask)=1000*imt1+jtask
        call MPI_ISEND (pp2(1,j),mp,pptype,dest,
     1                  msgtype(imt1,jtask),
     1                  MPI_COMM_WORLD,msgids,mpierr)
        call MPI_WAITALL (1,msgids,mpists,mpierr)
c
	end if
c
c force synchronization
c
!        call MPI_BARRIER (MPI_COMM_WORLD,mpierr)
c
	rtime2=MPI_WTIME()
	cpumpi=cpumpi+rtime2-rtime1
c
cif (taskid.ge.nsubset) go to 100
	if (mod(taskid,numtasks/nsubset).ne.0) go to 100
c
	rtime1=MPI_WTIME()
      lu=lupop+k-1+1000
c
      if (j.eq.j1.and.jstep.eq.1) then 
	write (6,903) taskid,lu,istep-1,time,j1,j2
 903  format ('pop2_write5: taskid,lu,istep-1,time,j1,j2=', 3i6,1p,e12.4,2i3)
	end if
c
	if (j.eq.j1) write (lu) istep-1,time
	
c
c for velocity components
c
c
      do 15 il=1,nl
      i1=(il-1)*npr+1
c fix by PK Yeung, 2/14/09
c     i2=min0(i1+npr-1,iplst-ipfst)
      i2=min0(i1+npr-1,iplst)
      write (lu) (pg2(i),i=i1,i2)
 15   continue
c
	rtime2=MPI_WTIME()
	cpuio=cpuio+rtime2-rtime1
c
 100   continue
c
	deallocate (pg2)
	deallocate (msgidr,msgtype,mpists,mpistr)

	cpu = MPI_WTIME() -rtime0
	if (taskid.eq.0) then
	if (jstep.eq.1) open (71,file='lag_timings/timings_pop2_write5')
	write (71,"('istep,cpumpi,cpuio,cpu=',i6,2i3,1p,3e11.3)") istep,j1,j2,cpumpi,cpuio,cpu
	end if 


c
#endif
#endif
      end
