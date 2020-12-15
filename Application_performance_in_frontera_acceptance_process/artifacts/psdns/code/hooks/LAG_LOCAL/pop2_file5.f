      subroutine pop2_file5 (iopt)
c
c
c new routine to handle opening and closing of files
c corresponding to wrt2pop5. 2/1/09
c 
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
c "fn" needs to be larger! fix by PK Yeung, 3/27/09
c     character*15 fn
      character*20 fn
c
      character*6 name,numer
      character*11 unform,form
      data unform/'unformatted'/
      data form/'formatted'/
      character*110 string
	character*30 caux
c
c
	integer iopt,mm,mpr,nl,k,lu,nchar,numc,ios,
     1          j,il,i1,i2,i,jp,ip,npr,jj,itask
c
        real*8 rtime1,rtime2,rtime3,rtime4
        real cpuio1,cpuio2,cpuio3
        save cpuio1,cpuio2,cpuio3
        real, allocatable :: cpuio_all(:,:)
c
	integer isubset, npset
c
      integer, allocatable :: ifile(:)
      integer icall
      data icall/0/
      save icall,ifile
c
      if (icall.eq.0) then
      allocate (ifile(nplu))
      ifile(:)=0
      icall=icall+1
        cpuio1=0.
        cpuio2=0.
        cpuio3=0.
      end if
c
c
      if (taskid.eq.1) 
     1	write (6,*) 'enter pop2_file5: istep,isflag,iopt=',istep,isflag,iopt
      if (taskid.eq.1) 
     1	write (6,"('enter pop2_file5: npo=',12i4)") npo
      if (taskid.eq.1.and.iopt.eq.0)
     1    write (6,*) 'pop2_file5 called with iopt=0'
c
!	if (taskid.ge.nsubset) return
	if (mod(taskid,numtasks/nsubset).ne.0) go to 90
c
c
c first and last particles in this group
c
	mm=nop2/nsubset
	npset = ndpart2/nsubset
c
c npr = no. of particles per record
c
      npr=min0(8192,npset)
      nl=npset/npr
      if (mod(npset,npr).ne.0) nl=nl+1
!      nl=nop2/nsubset/npr
!      if (mod(nop2/nsubset,npr).ne.0) nl=nl+1
c
      do 100 k=1,nplu-1
c
      if (npo(k).eq.0) go to 100
c
      if (k.eq.1) then
      name='pvel'
      else if (k.eq.nplu-1) then
      name='ppos'
      else if (k.eq.nplu-2) then
      name='pvlap'
      else if (k.eq.2) then
      name='pvgrad'
      else if (k.eq.3) then
      name='pdiss'
      else if (k.ge.4.and.k.le.3+ncop) then
      write (name,603) k-3
 603  format ('pscal',i1)
      end if
c
c All old "iwmss=1" action removed
c
      lu=lupop+k-1
c     lu=lu+20
      lu=lu+1000
c
      if (ifile(k).eq.0) then
      ifile(k)=ifile(k)+1
      fn=name//'1g0'
      call blanks (fn,nchar)
c
c       if (numtasks.le.numtasks) then
        if (numtasks.le.nsubset) then
        isubset=taskid
        else
	isubset=(taskid+1)/(numtasks/nsubset)+1
        end if
	write (numer,611) isubset
 611	format (i6)
        numc=1+floor(alog10(1.*isubset))
        fn=fn(1:nchar)//'_s'//numer(6-numc:6)
c
      call blanks (fn,nchar)
        caux='outppg0/'//fn(1:nchar)
	write (6,*) 'pop2_file5:',taskid,caux

        rtime3=MPI_WTIME()
      call fopen1 (lu,caux,unform)
        rtime4=MPI_WTIME()
        cpuio1=cpuio1+rtime4-rtime3
        rtime3=MPI_WTIME()
      write (lu) nop2/nsubset,npr,nl,npo(k)
        rtime4=MPI_WTIME()
        cpuio2=cpuio2+rtime4-rtime3
      end if
c
c close the file and begin a new one
c if the current lagrangian output step is
c a checkpointing step
c
c
      if (iopt.eq.0.or.isflag.eq.1) then
c
        rtime3=MPI_WTIME()
      close (lu)
        rtime4=MPI_WTIME()
        cpuio3=cpuio3+rtime4-rtime3
c
      if (iopt.eq.0) go to 100
c
      ifile(k)=ifile(k)+1
      write (fn,601) name,ifile(k)
 601  format (a6,i2,'g0')
      call blanks (fn,nchar)
c
	isubset=(taskid+1)/(numtasks/nsubset)+1
        numc=1+floor(alog10(1.*isubset))
	write (numer,611) isubset
        fn=fn(1:nchar)//'_s'//numer(6-numc:6)
c
	call blanks (fn,nchar)
        caux='outppg0/'//fn(1:nchar)
	write (6,*) 'pop2_file5:',taskid,caux
        rtime3=MPI_WTIME()
      call fopen1 (lu,caux,unform)
        rtime4=MPI_WTIME()
        cpuio1=cpuio1+rtime4-rtime3
        rtime3=MPI_WTIME()
      write (lu) npset,npr,nl,npo(k)
!      write (lu) nop2/nsubset,npr,nl,npo(k)
!	if(npset.ne.524288 .and. npr.ne.8192) then
!	write(6,*) 'error in pop2_file5, taskid,  npset,npr=',taskid,npset,npr
!	endif
        rtime4=MPI_WTIME()
        cpuio2=cpuio2+rtime4-rtime3
      else
c
      write (fn,601) name,ifile(k)
      call blanks (fn,nchar)
	isubset=(taskid+1)/(numtasks/nsubset)+1
        numc=1+floor(alog10(1.*isubset))
	write (numer,611) isubset
        fn=fn(1:nchar)//'_s'//numer(6-numc:6)
      call blanks (fn,nchar)
      open(lu,file='outppg0/'//fn(1:nchar),form=unform,iostat=ios,position='append')
c
      end if
c
 100	continue
c
 90	continue
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
        if (taskid.eq.numtasks-1) then
        caux='lag_timings/cpuio_pop2_file5'
        call fopen1 (7,caux,form)
        write (7,"('Accumulated time spent in opening and closing files, pop2_file5',
     1               '  at istep=',i6)") istep
        do itask=0,numtasks-1
        if (cpuio_all(itask,1).gt.0.) then
        write (7,"('taskid,cpuio=',i5,1p,4e12.4)") itask,cpuio_all(itask,1:3),
     1            sum(cpuio_all(itask,1:3))
        end if
        end do
        close (7)
        end if

        deallocate (cpuio_all)
        end if

        call time_stamp (' exit pop2_file5')


#endif
#endif
      return
      end
