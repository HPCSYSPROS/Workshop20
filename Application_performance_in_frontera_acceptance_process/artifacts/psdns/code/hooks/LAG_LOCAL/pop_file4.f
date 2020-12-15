	subroutine pop_file4 (iopt)
c
c new routine to handle opening and closing of files
c corresponding to wrtpop3
c
c P.K Yeung, 12/28/08
c
c multiple processors participating in writing the output, 7/4/08
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
      character*20 fn
	character*6 numer
c
      character*6 name
      logical secvel
      character*11 unform,form
      data unform/'unformatted'/
      data form/'formatted'/
      character*110 string
c
	integer iopt,mm,ipfst,iplst,mpr,nl,k,lu,nchar,numc,ios,
     1          j,il,i1,i2,i,jp,ip
	integer igp,nopg,nrefp,ifsttask,itask,ii,mrefp,m1,m2,npr,
     1          number,icall

	real dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,
     1       s12,s13,s23,term1,term2,term3,term4,omegax,omegay,omegaz
c
      integer, allocatable :: ifile(:,:)
      data icall/0/
      save icall,ifile
	character*40 caux
	character*2 cgid
c
	real*8 rtime1,rtime2,rtime3,rtime4
	real cpuio1,cpuio2,cpuio3,cpu
	save cpuio1,cpuio2,cpuio3,cpu
	real, allocatable :: cpuio_all(:,:)
c
c	call time_stamp ('enter pop_file4')
c	call MPI_BARRIER (MPI_COMM_WORLD,mpierr)
	rtime1=MPI_WTIME()

 	call time_stamp ('enter pop_file4')
      if (icall.eq.0) then
      allocate (ifile(nplu,ngp))
      ifile(:,:)=0
      icall=icall+1
	cpuio1=0.
	cpuio2=0.
	cpuio3=0.
	cpu=0.
      end if
c
       do 1000 igp=1,ngp
c
c note that, for the convenience of post-processing programs,
c reference particles are repeated for each group
c
      nopg=nop/(1+3*ngp)*4
      nrefp=nopg/4
c
#ifdef FEB2611
	ifsttask=(igp-1)*nsubset+1
	if (taskid.eq.0.or.taskid.gt.ngp*nsubset) go to 1000
	itask=mod(taskid,nsubset)
	ii=taskid/nsubset+1
	if (itask.eq.0) then 
	itask=nsubset
	ii=ii-1
	end if
	if (ii.ne.igp) go to 1000
#else
	ifsttask=numtasks-nsubset
        if (taskid.lt.ifsttask) go to 1000
	itask=taskid-ifsttask
#endif

	mrefp=nrefp/nsubset
	m1=itask*mrefp+1
	m2=m1+mrefp-1
c
cccc
	npo(1)=3
	npo(nplu-1)=3
cccc
c
	npr=mrefp
	nl=4

      do 100 k=1,nplu-1
c
      if (npo(k).eq.0) go to 100
c
#ifdef LAGSC2
!	if (k.ne.1.and.k.ne.nplu-1) go to 100
#endif
c
      if (k.eq.1) then
      name='pvel'
      else if (k.eq.nplu-1) then
      name='ppos'
      else if (k.eq.2)then
      name='pvgrad'
      else if (k.ge.4.and.k.le.3+ncop) then
      write (name,603) k-3
 603  format ('pscal',i1)
      end if
c
      lu=lupop+k-1+igp*10
c
      if (ifile(k,igp).eq.0) then
      write (fn,602) name,ifile(k,igp)+1,igp
 602  format (a6,i2,'g',i2)
	call blanks (fn,nchar)
c
	write (numer,611) itask
 611	format (i6)
	if (itask.eq.0) then
	numc=1
	else
	numc=1+floor(alog10(1.*itask))
	end if
	fn=fn(1:nchar)//'_s'//numer(6-numc:6)
c
      call blanks (fn,nchar)
      ifile(k,igp)=ifile(k,igp)+1
	write (cgid,"(i1,'/')") igp
      caux='outpp/'//cgid//fn(1:nchar)
	write (6,*) 'pop_file4:',istep,taskid,caux
c
	rtime3=MPI_WTIME()
      open (lu,file=caux,form='unformatted')
	rtime4=MPI_WTIME()
	cpuio1=cpuio1+rtime4-rtime3
c
	number=nopg/nsubset
c
	rtime3=MPI_WTIME()
         write (lu) number,npr,nl,npo(k)
	rtime4=MPI_WTIME()
	cpuio2=cpuio2+rtime4-rtime3
      end if
c
c close the file and begin a new one
c if the current lagrangian output step is
c a checkpointing step
c
        if (iopt.eq.0.or.isflag.eq.1) then
c
      ifile(k,igp)=ifile(k,igp)+1
	rtime3=MPI_WTIME()
      close (lu)
	rtime4=MPI_WTIME()
	cpuio3=cpuio3+rtime4-rtime3
c
      if (iopt.eq.0) go to 100
c
      write (fn,602) name,ifile(k,igp),igp
      call blanks (fn,nchar)
        write (numer,611) itask
        if (itask.eq.0) then
        numc=1
        else
        numc=1+floor(alog10(1.*itask))
        end if
        fn=fn(1:nchar)//'_s'//numer(6-numc:6)
      call blanks (fn,nchar)
	write (cgid,"(i1,'/')") igp
      caux='outpp/'//cgid//fn(1:nchar)
	write (6,*) 'file4:',istep,taskid,caux
c
	rtime3=MPI_WTIME()
c      call fopen1 (lu,'outpp/'//fn(1:nchar),unform)
      open (lu,file=caux,form='unformatted')
	rtime4=MPI_WTIME()
	cpuio1=cpuio1+rtime4-rtime3

	rtime3=MPI_WTIME()
      write (lu) nopg/nsubset,npr,nl,npo(k)
	rtime4=MPI_WTIME()
	cpuio2=cpuio2+rtime4-rtime3
c
      end if
c
c
 100  continue
c

c

 1000 continue
c
	rtime2=MPI_WTIME()
	cpu=cpu+ rtime2-rtime1
c
	if (iopt.eq.0) then
c
	allocate (cpuio_all(0:numtasks-1,4))
	call MPI_ALLGATHER (cpuio1,1,mpireal,cpuio_all(0,1),1,mpireal,
     1                   MPI_COMM_WORLD,mpierr)
	call MPI_ALLGATHER (cpuio2,1,mpireal,cpuio_all(0,2),1,mpireal,
     1                   MPI_COMM_WORLD,mpierr)
	call MPI_ALLGATHER (cpuio3,1,mpireal,cpuio_all(0,3),1,mpireal,
     1                   MPI_COMM_WORLD,mpierr)
	call MPI_ALLGATHER (cpu,1,mpireal,cpuio_all(0,4),1,mpireal,
     1                   MPI_COMM_WORLD,mpierr)
	if (taskid.eq.numtasks-1) then
	caux='lag_timings/cpuio_pop_file4'
 	call fopen1 (7,caux,form)
	write (7,"('Time spent in opening and closing files, pop_file4',
     1               '  at istep=',i6)") istep
	do itask=0,numtasks-1
 	if (cpuio_all(itask,1).gt.0.) then
	write (7,"('taskid,cpuio=',i5,1p,4e12.4)") itask,cpuio_all(itask,1),
     1             cpuio_all(itask,3:4)
 	end if
 	end do
	close (7)
	end if

	deallocate (cpuio_all)
	end if

 	call time_stamp (' exit pop_file4')
c
      return
c
#endif
      end
