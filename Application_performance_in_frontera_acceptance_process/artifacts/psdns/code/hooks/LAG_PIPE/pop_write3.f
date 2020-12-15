      subroutine pop_write3 (k,j1,j2)
c
c output action corresponding to wrtpop3
c
c output of particle properties
c
c
#ifdef LAG
c
	use compart
	implicit none
c
	integer i,j,ii,i1,i2,k,j1,j2,nopg,nrefp,mp,igp,igpw1,igpw2,lu,ipfst
	real(b8), allocatable :: pg(:)
c
	if (taskid.eq.0) write (6,*) 'enter pop_write3,k=',k
c
c note that, for the convenience of post-processing programs,
c reference particles are repeated for each group
c
      nopg=nop/(1+3*ngp)*4
      nrefp=nopg/4
c
c
c set no. of particles per record to the number of reference particles
c
c
	mp=nop/numtasks
c
	allocate (pg(ndpart))
c
c output of velocity components
c
c
	do 100 j=j1,j2
	write(6,*)'pop_write3: taskid, j1, j2=', taskid,j1,j2
c
	call MPI_ALLGATHER (pp(1,j),mp,mpireal,pg,mp,mpireal,
     1                      MPI_COMM_WORLD,mpierr)
c
	do 200 igp=1,ngp
c
      igpw1=igp/numtasks
      if (mod(igp,numtasks).eq.0) igpw1=igpw1-1
      igpw2=igp-igpw1*numtasks
c
	if (taskid+1.eq.igpw2) then
      lu=lupop+k-1+igpw1*10
c
	if (j.eq.j1) write (lu) istep-1,time
c
      write (lu) (pg(i),i=1,nrefp)
      ipfst=(igp-1)*3*nrefp+nrefp+1
      do ii=1,3
      i1=ipfst+(ii-1)*nrefp
      i2=i1+nrefp-1
      write (lu) (pg(i),i=i1,i2)
      end do
c
	end if
c
 200	continue
c
 100	continue
c
	deallocate (pg)

#endif
c
      end
