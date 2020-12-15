      subroutine mol_file3 (iopt)


! subroutine for opening and closing of files corresponding to molecules
! D. Buaria, 6/12/12, comments below need to be changed.
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
c if iopt=0, simply close the files without writing new data
c
#ifdef LAG
#ifdef MOL
c
	use compart
	implicit none
c
      character*20 fn
c
      character*6 name,numer
      character*11 unform
      data unform/'unformatted'/
      character*110 string
	character*30 caux
        character*1 cigm
        integer igm
c
c
	integer iopt,mm,mpr,nl,k,lu,nchar,numc,ios,
     1          j,il,i1,i2,i,jp,ip,npr,jj
c
	integer isubset
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
      end if
c
      if (taskid.eq.1) 
     1	write (6,*) 'enter mol_file3: istep,isflag,iopt=',istep,isflag,iopt
      if (taskid.eq.1.and.iopt.eq.0)
     1    write (6,*) 'mol_file3 called with iopt=0'
c
!	if (taskid.ge.nmsubset) return
	if (mod(taskid,numtasks/nmsubset).ne.0) go to 90
c
c
c first and last particles in this group
c
	mm=nom/nmsubset
c
c npr = no. of particles per record
c
      npr=min0(8192,mm)
      nl=nom/nmsubset/npr/ngm
      if (mod(nom/nmsubset/ngm,npr).ne.0) nl=nl+1
c
        nl=4
        npr=nom/nsubset/nl
c
      do 100 k=1,nplu-1
c
        if (k.ne.1.and.k.ne.nplu-1) go to 100
c
      if (npo(k).eq.0) go to 100
c
      if (k.eq.1) then
      name='mvel'
      else if (k.eq.nplu-1) then
      name='mpos'
      else if (k.eq.nplu-2) then
      name='mvlap'
      else if (k.eq.2) then
      name='mvgrad'
      else if (k.eq.3) then
      name='mdiss'
      else if (k.ge.4.and.k.le.3+ncop) then
      write (name,603) k-3
 603  format ('pscal',i1)
      end if
c
c
      if (ifile(k).eq.0) then
      ifile(k)=ifile(k)+1
      fn=name//'1g0'
      call blanks (fn,nchar)
c
	isubset=(taskid+1)/(numtasks/nmsubset)+1
	write (numer,611) isubset
 611	format (i6)
        numc=1+floor(alog10(1.*isubset))
        fn=fn(1:nchar)//'_s'//numer(6-numc:6)
      call blanks (fn,nchar)
c
        do igm=1,ngm
      lu=lupop+k-1+3000+igm*10
        write (cigm,"(i1)") igm
        caux='outmp/'//cigm//'/'//fn(1:nchar)
        write (6,*) 'mol_file3: caux=',caux
      call fopen1 (lu,caux,unform)
      write (lu) nom/nmsubset/ngm,npr,nl,npo(k)
        end do
      end if
c
c close the file and begin a new one
c if the current lagrangian output step is
c a checkpointing step
c
      if (iopt.eq.0.or.isflag.eq.1) then
c
        do igm=1,ngm
      lu=lupop+k-1+3000+igm*10
      call fclose (taskid,istep,lu)
        end do
c
      if (iopt.eq.0) go to 100
c
      ifile(k)=ifile(k)+1
      write (fn,601) name,ifile(k)
 601  format (a6,i2,'g0')
      call blanks (fn,nchar)
c
	isubset=(taskid+1)/(numtasks/nmsubset)+1
        numc=1+floor(alog10(1.*isubset))
	write (numer,611) isubset
        fn=fn(1:nchar)//'_s'//numer(6-numc:6)
	call blanks (fn,nchar)
c
        do igm=1,ngm
      lu=lupop+k-1+3000+igm*10
        write (cigm,"(i1)") igm
        caux='outmp/'//cigm//'/'//fn(1:nchar)
        write (6,*) 'mol_file3: caux=',caux
      call fopen1 (lu,caux,unform)
      write (lu) nom/nmsubset/ngm,npr,nl,npo(k)
        end do
c
c
      else
c
      write (fn,601) name,ifile(k)
      call blanks (fn,nchar)
	isubset=(taskid+1)/(numtasks/nmsubset)+1
        numc=1+floor(alog10(1.*isubset))
	write (numer,611) isubset
        fn=fn(1:nchar)//'_s'//numer(6-numc:6)
      call blanks (fn,nchar)
        do igm=1,ngm
      lu=lupop+k-1+3000+igm*10
      call fclose (taskid,istep,lu)
        write (cigm,"(i1)") igm
        caux='outmp/'//cigm//'/'//fn(1:nchar)
        write (6,*) 'mol_file3: caux=',caux
      open(lu,file=caux,form=unform,iostat=ios,position='append')
        end do
c
      end if
c
 100	continue
c
 90     continue
#endif
#endif
      return
      end
        subroutine fclose (taskid,istep,lu)
        integer taskid,lu,istep
        if (lu.gt.3000) then
        write (6,*) 'fclose: taskid,istep,lu=',taskid,istep,lu
        close (lu)
        end if
        return
        end
        
