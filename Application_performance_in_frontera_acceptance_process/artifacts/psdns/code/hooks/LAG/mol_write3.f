      subroutine mol_write3 (k,j1,j2)
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
#ifdef MOL
c
	use compart
	implicit none
c
      character*15 fn
      character*6 name,numer
      character*11 unform
      data unform/'unformatted'/
	character*30 caux

c
	real(p8), allocatable :: mpg2(:)
c
 	integer j,j1,j2,k,mp
c
	integer ipfst,iplst,nl,i,i1,i2,il,mm,npr,lu
c
	real*8 rtime1,rtime2
	real cpumpi,cpuio
c
cccc 
	integer iraddr,itask,jtask,irtask1,irtask2,dest,msgids,imt1
        integer, allocatable :: msgidr(:),msgtype(:,:)
        integer, allocatable :: mpistr(:,:),mpists(:)
c
	integer isubset,igm
c
	allocate (mpg2(nom/nmsubset/ngm))
c
c first and last particles in this group
c
c

	isubset=(taskid+1)/(numtasks/nsubset)+1
	mm=nom/nmsubset/ngm
c
c npr = no. of particles per record
c
      npr=min0(8192,mm)
      nl=nom/nmsubset/npr/ngm
      if (mod(nom/nmsubset/ngm,npr).ne.0) nl=nl+1
cc
        nl=4
        npr=nom/nmsubset/ngm/nl
c
	mp=nom/numtasks/ngm
c
	cpumpi=0.
	cpuio=0.
c
        allocate (msgidr(numtasks/nmsubset))
        allocate (msgtype(nmsubset,numtasks/nmsubset-1))
        allocate (mpists(MPI_STATUS_SIZE))
        allocate (mpistr(MPI_STATUS_SIZE,numtasks/nmsubset))

        do 10 igm=1,ngm
c
      lu=lupop+k-1+3000+igm*10
c
	do 100 j=j1,j2
c
	rtime1=MPI_WTIME()
c
        if (mod(taskid,numtasks/nmsubset).eq.0) then
c
c post the receives
c taskids 0, nsubset, 2*nsubset, 3*nsubset, etc to receive messages
c from the next numtasks/nsubset - 1 tasks
c
c
        irtask1=taskid+1
        irtask2=irtask1+numtasks/nmsubset-2
c
        do itask=irtask1,irtask2
        jtask=itask-irtask1+1
        imt1=itask/(numtasks/nmsubset)+1
        msgtype(imt1,jtask)=1000*imt1+jtask
        iraddr=jtask*mp+1
        call MPI_IRECV (mpg2(iraddr),mp,pptype,itask,
     1                  msgtype(imt1,jtask),
     1                  MPI_COMM_WORLD,msgidr(jtask),mpierr)
        end do
c
	mpg2(1:mp)=mpp(1:mp,j,igm)
        call MPI_WAITALL (numtasks/nmsubset-1, msgidr,mpistr,mpierr)
c
        else
c
c post the sends
c
        imt1=taskid/(numtasks/nmsubset)+1
        dest=taskid/(numtasks/nmsubset)*(numtasks/nmsubset)
        jtask=taskid-dest
        msgtype(imt1,jtask)=1000*imt1+jtask
        call MPI_ISEND (mpp(1,j,igm),mp,pptype,dest,
     1                  msgtype(imt1,jtask),
     1                  MPI_COMM_WORLD,msgids,mpierr)
        call MPI_WAITALL (1,msgids,mpists,mpierr)
c
	end if
c
c force synchronization
c
        call MPI_BARRIER (MPI_COMM_WORLD,mpierr)
c
	rtime2=MPI_WTIME()
	cpumpi=cpumpi+rtime2-rtime1
c
cif (taskid.ge.nmsubset) go to 100
	if (mod(taskid,numtasks/nmsubset).ne.0) go to 100
c
	rtime1=MPI_WTIME()
c
	if (j.eq.j1) write (lu) istep-1,time
c
c
c for velocity components
c
      do 15 il=1,nl
      i1=(il-1)*npr+1
      i2=min0(i1+npr-1,mm)
      write (lu) (mpg2(i),i=i1,i2)
 15   continue
c
	rtime2=MPI_WTIME()
	cpuio=cpuio+rtime2-rtime1
c
 100   continue
c
 10     continue
c
	deallocate (mpg2)
	deallocate (msgidr,msgtype,mpists,mpistr)

c
#endif
#endif
      end
