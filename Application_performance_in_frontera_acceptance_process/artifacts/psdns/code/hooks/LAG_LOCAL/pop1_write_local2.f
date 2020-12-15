      subroutine pop1_write_local2 (k,j1,j2)
c
c alternative to pop1_write4, 2/1/09
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
	character*6 string2
	integer nchar
c
	real(p8), allocatable :: pptmp(:,:),rbuf(:)
	integer, allocatable :: ppindtmp(:,:),rindbuf(:)
c
 	integer j,j1,j2,k,mp,ic,indx,ic2
c
	integer ipfst,iplst,nl,i,i1,i2,il,mm,npr,lu,npset
	integer npuset, npu, ii, igp
c
	real*8 rtime1,rtime2,rtime0
	real cpumpi,cpuio,cpu,estmem
c
cccc 
	integer iraddr,itask,jtask,irtask1,irtask2,dest,msgids,imt1
        integer, allocatable :: msgidr(:),msgtype(:,:)
        integer, allocatable :: mpistr(:,:),mpists(:)
	integer, allocatable :: ipout(:),ipin(:),ipall(:),ipall2(:),idisp(:)
c
	integer isubset,nid,nvar,iset,jj,itemp
c

	if (taskid.eq.0.and.jstep.eq.1) write (6,*) 'enter pop1_write_local, 
     1                    istep,j1,j2=',istep,j1,j2
 	call time_stamp ('enter pop1_write_local2')

	rtime0 = MPI_WTIME()

	cpumpi=0.
	cpuio=0.

	npset = ndpart/nsubset 

	npu = ndpart/(1+3*ngp)
	npuset = npu/nsubset

! index of particle in pp is now used from pp1ind array
!	nvar = j2-j1+2
	nvar = j2-j1+1

	isubset=(taskid+1)/(numtasks/nsubset)+1
c
! npr = no. of particles per record
c
      npr=min0(8192,npset)
      nl=npset/npr
      if (mod(npset,npr).ne.0) nl=nl+1
c

	sync memory
	call MPI_BARRIER (MPI_COMM_WORLD,mpierr)


	if(istep.eq.1) then
	if (taskid.eq.0) then
	write(6,*) 'memory usage for pp1 group in IO'
	estmem=nump1*(nvar+1)*nsubset+(nump1+1)*nsubset
	estmem=estmem*p8/1024./1024.
	write (6,"('memory for pptmp and ppindtmp arrays:',f8.1,'  MB')") estmem
	end if

	if (taskid.eq.0) then
	estmem=(npset*nvar+numtasks)+(npset+numtasks)
	estmem=estmem*p8/1024./1024.
	write (6,"('memory for rbuf and rinbuf arrays:',f8.1,'  MB')") estmem
	end if
	endif

! sendbuffer arrays to sort particles on each MPI task according to the 
! destination MPI task  it will be going to for IO
	allocate (pptmp(nump1*nvar+1,nsubset))
	allocate (ppindtmp(nump1+1,nsubset))
! receive buffer arrays to get information on destination MPI tasks
	allocate (rbuf(npset*nvar + numtasks))
	allocate (rindbuf(npset + numtasks))
! integer arrays required for counts and displacements for MPI_Gather
	allocate (ipout(nsubset),ipall(numtasks))
	allocate (ipall2(numtasks),idisp(numtasks))


	pptmp = 0.
	ipout(:) = 0


! sort all particles according to the which MPI task
! they will go to for IO
	do i=1,nump1
	iset = (pp1ind(i) - 1)/npset + 1
	ppindtmp(ipout(iset) + 1,iset) = pp1ind(i)
	do j=j1,j2
	pptmp(nvar*ipout(iset) + j-j1+1,iset) = pp(i,j)
	enddo
	ipout(iset) = ipout(iset) + 1 
	enddo

	ipall(:)=0


	rtime1=MPI_WTIME()

	do jj=1,nsubset

! nid = MPI which does IO
! send info about how many particles each "nid" 
! will receive from various  MPI tasks

	nid = (jj-1)*numtasks/nsubset
        call MPI_GATHER (ipout(jj), 1, MPI_INTEGER, ipall, 1, MPI_INTEGER,
     1                   nid, MPI_COMM_WORLD,mpierr)

	enddo

	if(sum(ipall).ne.0) then
        if(sum(ipall).ne.npset) then
        write(6,*) 'error in pop1_write_local: nump1,sum ipall,taskid=',nump1,sum(ipall),taskid
        endif
	endif

#ifdef DEBUG
	write (string2,"(i6)") taskid
	call blanks (string2,nchar)
	caux='psort/nump1.'//string2(1:nchar)
	if (jstep.eq.1) then
	open (120,file=caux)
	else
	open (120,file=caux,position='append')
	end if
	write(120,*) 'in pop1_write_local: nump1,sum ipall,taskid=',nump1,sum(ipall),taskid
	close (120)
#endif

! send the variables from j1 to j2 to their respective "nid"
	idisp(1)=0
	ipall2(1) = ipall(1)*nvar + 1
	do i=2,numtasks
	ipall2(i) = ipall(i)*nvar + 1
	idisp(i) = idisp(i-1) + ipall2(i-1)
	enddo
	do jj=1,nsubset
	nid = (jj-1)*numtasks/nsubset
        call MPI_GATHERV (pptmp(1,jj), ipout(jj)*nvar+1, pptype, rbuf,
     1                   ipall2, idisp, pptype, nid, MPI_COMM_WORLD,mpierr)
	enddo

	call MPI_BARRIER (MPI_COMM_WORLD, mpierr)

! send the particle indices 
	idisp(1)=0
	ipall2(1) = ipall(1) + 1
	do i=2,numtasks
	ipall2(i) = ipall(i) + 1
	idisp(i) = idisp(i-1) + ipall2(i-1)
	enddo
	do jj=1,nsubset
	nid = (jj-1)*numtasks/nsubset
        call MPI_GATHERV (ppindtmp(1,jj), ipout(jj)+1, MPI_INTEGER, rindbuf,
     1                   ipall2, idisp, MPI_INTEGER, nid, MPI_COMM_WORLD,mpierr)
	enddo

	deallocate (pptmp)
	deallocate (ppindtmp)

	rtime2=MPI_WTIME()
	cpumpi=cpumpi+rtime2-rtime1


	allocate (pptmp(npset,nvar))
!	allocate (ppindtmp(npset,1))

        if (mod(taskid,numtasks/nsubset).eq.0) then

	ic = 1
	ic2 = 1

	if(j1.eq.1) write(8000+taskid,*) 'j1,j2,istep=',j1,j2,istep
	do i=1,numtasks

! sort the particles according to its index 
! this will result in output written out as old format

	if(ipall(i).ne.0) then
	do j=1,ipall(i)
	if(j1.eq.1) write(8000+taskid,*) rindbuf(ic2)
	indx = mod(rindbuf(ic2) - 1, npset) + 1
	pptmp(indx,1:nvar) = rbuf(ic:ic+nvar-1)
!	ppindtmp(indx,1) = rindbuf(ic2)
	ic = ic + nvar 
	ic2 = ic2 + 1
	enddo
	endif

	ic = ic + 1
	ic2 = ic2 + 1

	enddo

#ifdef TESTING
	if(taskid.lt.1000) then
	if(j1.eq.1) then
	write(8000+taskid,*) 'j1,j2,nvar=',j1,j2,nvar
	do i=1,npset
	write(8000+taskid,*) i, ppindtmp(i,1)
	enddo
	close(8000+taskid)
	endif
	endif
#endif

	rtime1=MPI_WTIME()

! write out the variables fron j1 to j2
	do 100 j=j1,j2

	do 200 igp=1,ngp

	lu=lupop+k-1+igp*10
c
      if (j.eq.j1.and.jstep.eq.1) then 
	write (6,903) taskid,lu,istep-1,time,j1,j2
	end if
c
	if (j.eq.j1) write (lu) istep-1,time

	write (lu) (pptmp(i,j-j1+1),i=1,npuset)
        
!	ipfst=nrefp+3*(m1-1)+1+(igp-1)*3*nrefp
!        i1=ipfst+(ii-1)*mrefp
!        i2=i1+mrefp-1
        ipfst = npuset + (igp-1)*3*npuset + 1
        do 15 ii=1,3
        i1 = ipfst + (ii-1)*npuset
        i2 = i1 + npuset - 1
        write (lu) (pptmp(i,j-j1+1),i=i1,i2)
 15     continue

!      do 15 il=1,nl
!      i1=(il-1)*npr+1
!      i2=min0(i1+npr-1,npset)
!      write (lu) (pptmp(i,j-j1+1),i=i1,i2)
! 15   continue
c
c
 200   continue
 100   continue

	rtime2=MPI_WTIME()
	cpuio=cpuio+rtime2-rtime1

	endif    ! if (mod(taskid,numtasks/nsubset).eq.0) then
c
 903  format ('pop1_write_local: taskid,lu,istep-1,time,j1,j2=', 3i6,1p,e12.4,2i3)

	call MPI_BARRIER (MPI_COMM_WORLD,mpierr)

	cpu = MPI_WTIME() -rtime0
	if (taskid.eq.0) then
	if (jstep.eq.1) open (71,file='lag_timings/timings_pop1_write_local')
	write (71,"('istep,cpumpi,cpuio,cpu=',i6,2i3,1p,3e11.3)") istep,j1,j2,cpumpi,cpuio,cpu
	end if 

	deallocate (pptmp,rbuf,rindbuf)
	deallocate (ipout,ipall,ipall2,idisp)

 	call time_stamp ('exit pop1_write_local2')

	return
c
#endif
#endif
      end
