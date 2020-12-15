      subroutine pop_write4 (k,j1,j2)
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
c
c
	do 100 j=j1,j2
c
        call MPI_ALLGATHER (pp(1,j),mp,pptype,pg,mp,pptype,
     1                      MPI_COMM_WORLD,mpierr)
c
      do 200 igp=1,ngp
c
c note that, for the convenience of post-processing programs,
c reference particles are repeated for each group
c
c FEB26011: when this directive is activated a total 
c of ngp x nsubset MPI tasks will be writing the Lagrangian
c time series (for the main ensemble). Otherwise, 
c only nsubset tasks will do it, each writing multiple groups.
c Differences in timing may be small (not clear yet).
c
#ifdef FEB2611
        ifsttask=(igp-1)*nsubset+1
        if (taskid.eq.0.or.taskid.gt.ngp*nsubset) go to 200
        itask=mod(taskid,nsubset)
        ii=taskid/nsubset+1
        if (itask.eq.0) then
        itask=nsubset
        ii=ii-1
        end if
        if (ii.ne.igp) go to 200
	mrefp=nrefp/nsubset
	m1=(itask-1)*mrefp+1
	m2=m1+mrefp-1
#else
        ifsttask=numtasks-nsubset
        if (taskid.lt.ifsttask) go to 200
        itask=taskid-ifsttask
	mrefp=nrefp/nsubset
	m1=itask*mrefp+1
	m2=m1+mrefp-1
#endif


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
	deallocate (pg)
c
	call time_stamp (' exit pop_write4')
c
#endif
      return
      end
