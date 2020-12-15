      subroutine pop1_write_local (k,j1,j2)
c
c output action corresponding to wrtpop4
c
c Comment added on 2/2/09: replacement of MPI_ALLGATHER
c by discrete sends and receives (as in pop2_write5.f)
c for this routine is difficult unless the reference particles
c that are shared by tetrads of different initial size
c are stored separately in a different array.
c
c more processors participating in writing the output, 7/4/08
c
c output of particle properties
c
c MPI version, adapted for relative diffusion with more than
c one group of particle pairs, when (i) shear flow is simulated
c or (ii) particle pairs initialized using the tetrahedron scheme
c
c if iopt=0, simply close the files without writing new data
c
#ifdef LAG
c
	use compart
	implicit none
c
	integer j1,j2,nopg,nrefp,mp,igp,ifsttask,itask,mrefp,m1,m2,npr
	integer k,j,nl,lu,i,ipfst,ii,i1,i2
c
        real(p8), allocatable :: pg(:)
	real(p8), allocatable :: pptmp(:)
	integer, allocatable :: ppindtmp(:), num_all(:), idisp(:)
c
	real*8 rtime1,rtime2
	real cpu
c
      nopg=nop/(1+3*ngp)*4
      nrefp=nopg/4
c
	call MPI_BARRIER (MPI_COMM_WORLD,mpierr)
c
	call time_stamp ('enter pop_write4')
	mp=nop/numtasks
c
	if (taskid.eq.0) write (6,*) 'enter pop_write4,istep,j1,j2=',
     1 istep,j1,j2

	allocate (pg(ndpart))
	allocate (ppindtmp(ndpart), pptmp(ndpart))
	allocate (num_all(numtasks), idisp(numtasks))
c
        call MPI_ALLGATHER (nump1, 1 , MPI_INTEGER, num_all, 1, 
     1                      MPI_INTEGER, MPI_COMM_WORLD, mpierr)

        idisp(1) = 0
        do itask=2,numtasks
        idisp(itask) = idisp(itask-1) + num_all(itask-1)
        enddo

        call MPI_ALLGATHERV (pp1ind(1), nump1, MPI_INTEGER, ppindtmp(1), 
     1            num_all, idisp, MPI_INTEGER, MPI_COMM_WORLD,mpierr)

#ifdef DEBUG
	if(taskid.eq.0) then
	do i=1,ndpart
	write(8000,*) i, ppindtmp(i)
	enddo
	close(8000)
	endif
	call MPI_BARRIER (MPI_COMM_WORLD, mpierr)
#endif
c
	do 100 j=j1,j2
c
        call MPI_ALLGATHERV (pp(1,j), nump1, pptype, pptmp(1), 
     1            num_all, idisp, pptype, MPI_COMM_WORLD,mpierr)
!        call MPI_ALLGATHER (pp(1,j),mp,pptype,pg,mp,pptype,
!     1                      MPI_COMM_WORLD,mpierr)
c
	do i=1,ndpart
	pg(ppindtmp(i)) = pptmp(i)
	enddo

      do 200 igp=1,ngp
c
        ifsttask=numtasks-nsubset
        if (taskid.lt.ifsttask) go to 200
        itask=taskid-ifsttask
	mrefp=nrefp/nsubset
	m1=itask*mrefp+1
	m2=m1+mrefp-1


c
	npr=mrefp
	nl=4
      lu=lupop+k-1+igp*10

      if (j.eq.j1) write (lu) istep-1,time
c
c
      write (lu) (pg(i),i=m1,m2)
	
	ipfst=nrefp+3*(m1-1)+1+(igp-1)*3*nrefp
	do 15 ii=1,3
	i1=ipfst+(ii-1)*mrefp
	i2=i1+mrefp-1
	write (lu) (pg(i),i=i1,i2)
 15	continue
c
 200  	continue
c
 100	continue
c
	deallocate (pg,pptmp,ppindtmp)
c
	call time_stamp (' exit pop_write4')
c
#endif
      return
      end
