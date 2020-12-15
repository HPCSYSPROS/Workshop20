      subroutine wrtpop1 (iopt)
c
c updated for use in 2D code, 6/3/07
c
c output of particle properties
c
c this routine is used if there is only 1 "group" of particles
c (otherwise, use wrtpop2 or wrtpop3)
c
c if iopt=0, simply close the files without writing new data
c
#ifdef LAG
c
	use compart
c
      character*7 fn
      character*6 name
      character*11 unform
      logical secvel
      data unform/'unformatted'/
      character*110 string
c
      integer, allocatable :: ifile(:)
      data icall/0/
      save icall,ifile
c
      if (icall.eq.0) then
      allocate (ifile(nplu))
      ifile(:)=0
      icall=icall+1
      end if
c
	write (6,*) 'enter wrtpop1: istep,isflag,iopt=',istep,isflag,iopt
      if (taskid.eq.1.and.iopt.eq.0)
     1    write (6,*) 'wrtpop1 called with iopt=0'
c
c first and last particles in this group
c
      ipfst=1
      iplst=nop
c
c npr = no. of particles per record
c
c     npr=min0(20000,nop)
      npr=min0(8192,nop)
      nl=nop/npr
      if (mod(nop,npr).ne.0) nl=nl+1
c
c for now (8/10/96),
c task0 to output the velocities
c task1 to output the positions, mechanical and scalar dissipations
c task2 to output the scalars and scalar gradients
c task3 to output the velocity gradients
c
      do 100 k=1,nplu-1
c
      if (npo(k).eq.0) go to 100
c
      if (k.eq.1) then
      if (taskid.ne.0) go to 100
      name='pvel'
      else if (k.eq.nplu-1) then
      if (taskid.ne.1) go to 100
      name='ppos'
      else if (k.eq.nplu-2) then
      if ((numtasks.ge.4.and.taskid.ne.3)) go to 100
      if ((numtasks.lt.4.and.taskid.ne.1)) go to 100
      name='pvlap'
      else if (k.eq.2) then
      if ((numtasks.ge.4.and.taskid.ne.3)) go to 100
      if ((numtasks.lt.4.and.taskid.ne.1)) go to 100
      name='pvgrad'
      else if (k.eq.3) then
      if ((numtasks.ge.4.and.taskid.ne.2)) go to 100
      if ((numtasks.lt.4.and.taskid.ne.0)) go to 100
      name='pdiss'
      else if (k.ge.4.and.k.le.3+ncop) then
      if ((numtasks.ge.4.and.taskid.ne.2)) go to 100
      if ((numtasks.lt.4.and.taskid.ne.0)) go to 100
      write (name,603) k-3
 603  format ('pscal',i1)
      end if
c
c All old "iwmss=1" action removed
c
      lu=lupop+k-1
c
      if (ifile(k).eq.0) then
      ifile(k)=ifile(k)+1
      fn=name//'1'
      call blanks (fn,nchar)
      write (6,*)
     1     'wrtpop1: taskid,lu,k,fn(1:nchar)=',taskid,lu,k,fn(1:nchar)
      call fopen1 (lu,fn(1:nchar),unform)
      write (lu) nop,npr,nl,npo(k)
      end if
c
c close the file and begin a new one
c if the current lagrangian output step is
c a checkpointing step
c
#ifdef OLD
      if ((jstep-1.ge.iofreq.and.psave.gt.0.and.mod(jstep-1,psave).eq.0).
     1   .or.(iopt.eq.0.and.psave.ge.0).or.isflag.eq.1) then
#endif
c
      if (iopt.eq.0.or.isflag.eq.1) then
c
      close (lu)
c
      if (iopt.eq.0) go to 100
c
      ifile(k)=ifile(k)+1
      write (fn,601) name,ifile(k)
 601  format (a6,i1)
      call blanks (fn,nchar)
      write (6,*)'wrtpop1: taskid,lu,fn(1:nchar)=',taskid,lu,fn(1:nchar)
      call fopen1 (lu,fn(1:nchar),unform)
      write (lu) nop,npr,nl,npo(k)
c
      else
c
      write (fn,601) name,ifile(k)
      call blanks (fn,nchar)
      open(lu,file=fn(1:nchar),form=unform,iostat=ios,position='append')
c
      end if
c
      write (lu) istep-1,time
c
c for velocity components
c
      if (k.eq.1) then
c
      do 10 j=4,6
c
      do 15 il=1,nl
      i1=(il-1)*npr+ipfst
      i2=min0(i1+npr-1,iplst)
      write (lu) (pp(i,j),i=i1,i2)
 15   continue
c
 10   continue
c
c for (1st) velocity gradients
c
      else if (k.eq.2) then
c
      do 20 j=1,8
      if (ipvg(j).eq.0) go to 20
      jp=ipvg(j)
c
      do 25 il=1,nl
      i1=(il-1)*npr+ipfst
      i2=min0(i1+npr-1,iplst)
      write (lu) (pp(i,jp),i=i1,i2)
 25   continue
c
 20   continue
c
      else if (k.eq.3) then
c
c for mechanical and scalar dissipations
c
      do 30 j=1,1+ncop
      if (ipdis(j).eq.0) go to 30
      jp=ipdis(j)
c
      do 35 il=1,nl
      i1=(il-1)*npr+ipfst
      i2=min0(i1+npr-1,iplst)
      write (lu) (pp(i,jp),i=i1,i2)
 35   continue
c
 30   continue
c
c     else if (k.ge.4.and.k.lt.nplu-1) then
      else if (k.ge.4.and.k.lt.nplu-2) then
c
c for scalars and scalar gradients,
c similarly for velocity and scalar second derivatives
c
c in the case of velocity second derivatives, because of continuity
c the derivatives of dw/dz need not be written on the files
c if second gradients of u & v already are
c
      secvel=.false.
      if (k.eq.6+ncop.and.lpsecd(1).eq.1.and.lpsecd(2).eq.1.
     1    and.lpsecd(3).eq.1) secvel=.true.
c
      do 40 j=ipj1(k-3),ipj2(k-3)
c
      jj=j-ipj1(k-3)+1
      if (secvel.and.(jj.eq.3.or.jj.ge.5)) go to 40
c
      do 45 il=1,nl
      i1=(il-1)*npr+ipfst
      i2=min0(i1+npr-1,iplst)
      write (lu) (pp(i,j),i=i1,i2)
 45   continue
c
 40   continue
c
      else if (k.eq.nplu-2) then
c
c for velocity Laplacians
c
      do 50 j=ipj1(4+nc),ipj2(4+nc)
c
      do 55 il=1,nl
      i1=(il-1)*npr+ipfst
      i2=min0(i1+npr-1,iplst)
      write (lu) (pp(i,j),i=i1,i2)
 55   continue
c
 50   continue
c
      else if (k.eq.nplu-1) then
c
c for particle positions
c
      do 60 j=1,3
c
      do 65 il=1,nl
      i1=(il-1)*npr+ipfst
      i2=min0(i1+npr-1,iplst)
      write (lu) (pp(i,j),i=i1,i2)
 65   continue
c
 60   continue
c
      end if
c
 100  continue
c
      if (taskid.eq.1.and.iopt.ne.0)
     1    write (6,*) ' particle properties written, step=',istep-1
c

#ifdef MATP
	call wrtpomp2 (iopt)
#endif
c
#endif
      return
      end
