      subroutine wrtpos_pipe (lu)
c
c PKY, 6/29/2014, for pp2 in pipelined algorithm
c
c save particle positions on unformatted disk file
c (for purpose of checkpointing)
c
c --- MPI version 
c
c task 0 for fluid particles
c tasks 2 and numtasks-3 for molecules
c
#ifdef LAG
c
	use compart
	implicit none
c
	integer nn,numc,npr
	integer i,i1,i2,k,nr,lu,j,mp
	character*12 caux
	character*2 numer 
c
	real(p8) avx,avx0
	real(p8), allocatable :: pptmp(:,:)
	real(p8), allocatable :: pptmp1(:,:)
	real(p8), allocatable :: pptmp2(:,:)
	integer itask
c
	integer sendtag,recvtag,isource,idest
	integer mpistatus(MPI_STATUS_SIZE)
	integer, allocatable :: tag(:)

	integer, allocatable :: pstatus(:,:)
c
	integer relay,np
c
	character*11 form,unform
      data form,unform/'formatted','unformatted'/


c npr = no of particles written per record
c
c      data npr/100/
      data npr/1024/
c
	if (taskid.eq.0) write (6,*) 'enter wrtpos_pipe'
	caux='savepos'
      if (taskid.eq.0) call fopen1 (lu,caux,unform)
c
c
	if (taskid.eq.numtasks-1) then
        write (numer,600) ichkpt
 600    format (i2)
        numc=1+floor(alog10(1.*ichkpt))
        caux='savepos.'//numer(3-numc:2)
        call blanks (caux,nn)
      call fopen1 (lu,caux(1:nn),unform)
	end if
c
c use pptmp as scratch global array, since (see partsp.f)
c new interpolation weights are to be calculated next
c
	mp=nop/numtasks
c
	allocate (pptmp(nop,3))
	do j=1,3
	call MPI_ALLGATHER (pp(1,j),mp,pptype,pptmp(1,j),mp,pptype,
     1                      MPI_COMM_WORLD,mpierr)
	end do
c

      if (taskid.eq.0.or.taskid.eq.numtasks-1) then  
c
	npr=min0(nop,npr)
      write (lu) nop,npr
c
      write (lu) xyzl
      write (lu) gx(1),gy(1),gz(1)
c
c determine no. of records (nr)
c
      nr=nop/npr
      if (mod(nop,npr).ne.0) nr=nr+1
c
      do 10 k=1,nr
      i1=(k-1)*npr+1
      i2=min0(i1+npr-1,nop)
      write (lu) (pptmp(i,1),pptmp(i,2),pptmp(i,3),i=i1,i2)
 10   continue
c
      call fclose1 (lu)
      write (6,*) 'particle positions saved, step=',istep-1,pptmp(4,1)
c
      end if
c
	deallocate (pptmp)
c
#ifdef LAGSC2

c
c
c
        caux='save2pos'
      if (taskid.eq.1) call fopen1 (lu,caux,unform)
c
	if (taskid.eq.numtasks-2) then
        write (numer,"(i2)") ichkpt
        numc=1+floor(alog10(1.*ichkpt))
        caux='save2pos.'//numer(3-numc:2)
        call blanks (caux,nn)
      call fopen1 (lu,caux(1:nn),unform)
	end if
c
	allocate (tag(0:numtasks-1))
	allocate (pstatus(MPI_STATUS_SIZE,0:numtasks-1))
	allocate (pptmp1(nop2/numtasks,3))
	allocate (pptmp2(nop2/numtasks,3))
c
c pack particle positions on each MPI task into the pptmp1 array
c
	np=nop2/numtasks
	do j=1,3
	do i=1,np
	pptmp1(i,j)=pp2(j,i)
	end do
	end do
c
c every MPI task to send data to task 1, one by one
c
	relay=0
c
	idest=taskid+1
	isource=taskid-1
c
	if (taskid.eq.0) idest=2
	if (taskid.eq.2) isource=0
	recvtag=isource
	sendtag=taskid
	
	do itask=0,numtasks-1
	tag(itask)=itask
	end do
c
	npr=np
c
	if (taskid.eq.1) then
      write (lu) nop2,npr
      write (lu) xyzl
      write (lu) gx(1),gy(1),gz(1)
	end if
c
	pptmp1(2,2)=taskid+.5
	avx0=0.
c
	do 20 itask=0,numtasks-1
c
	if (itask.eq.1.and.taskid.eq.1) then
	avx=0.
	do i=1,np
	avx=avx+pptmp1(i,1)
	end do
	avx=avx/np
	avx0=avx0+avx
      write (lu) (pptmp1(i,1),pptmp1(i,2),pptmp1(i,3),i=1,np)
	go to 20
	end if

	if (taskid.eq.1) then
	call MPI_RECV (pptmp2,np*3,pptype,itask,tag(itask),
     1                 MPI_COMM_WORLD,pstatus(1,0),mpierr)
	avx=0.
	do i=1,np
	avx=avx+pptmp2(i,1)
	end do
	avx=avx/np
	avx0=avx0+avx
      write (lu) (pptmp2(i,1),pptmp2(i,2),pptmp2(i,3),i=1,np)
	else if (taskid.eq.itask) then
	call MPI_SSEND (pptmp1,np*3,pptype,1,tag(itask),
     1                 MPI_COMM_WORLD,mpierr)
	end if
c
 20	continue
     
	if (taskid.eq.1) write (6,*) 'wrtpos_pipe,avx0=',avx0/numtasks
      if  (taskid.eq.1) call fclose1 (lu)
c
#ifdef TEST
	if (taskid.eq.0) then
	relay=-1
	write (6,*) 'wrtpos_pipe: taskid=',taskid,relay
c
	else 
 	if (taskid.ne.1) then
	call MPI_RECV (relay,1,MPI_INTEGER,isource,recvtag,
     1                 MPI_COMM_WORLD,mpistatus,mpierr)
 	end if
	if (taskid.ne.1) write (6,*) 'wrtpos_pipe: taskid=',taskid,relay

	end if
c
	if (taskid.lt.numtasks-1) then
	if (taskid.ne.1) then
	relay=relay+taskid
	call MPI_SSEND (relay,1,MPI_INTEGER,idest,sendtag,
     1                  MPI_COMM_WORLD,mpierr)
	end if
	end if
#endif
	
c
	deallocate (tag)
	deallocate (pptmp1,pptmp2)
c
	if (taskid.eq.0) write (6,*) ' exit wrtpos_pipe'

#endif

#endif
      return
      end
