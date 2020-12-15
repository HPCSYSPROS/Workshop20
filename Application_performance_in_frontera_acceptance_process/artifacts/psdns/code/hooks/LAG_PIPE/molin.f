      subroutine molin
c
#ifdef LAG
#ifdef MOL
c
	use compart
	implicit none

      real lx,ly,lz,lx1,ly1,lz1
	integer nompart
	real(p8), allocatable :: mppin(:,:)
	real(p8) px,py,pz

	integer npr,igm,i,j1,j2,j,lu,nchar,nn
	real(b8) gx1,gy1,gz1
      data npr/100/
      character*75 string
        real(b8) avx0(3)
        real(p8) mpmin,mpmax
c


c set initial molecular positions corresponding to fluid particles
c
c to be called after sub. partin
c
c nompp = no. of molecules per fluid particles PER GROUP
c nom   = nompp x nop x ngm = total no. of molecules
c
! this part now moved to compart_set
!
!      nom=nop*nompp*ngm
!
!	allocate(mpp(nom,3),bfmxyz(nom,3),ibmxyz(nom,3))
!
! all taskids now initialize their own molecules' positions 
!      if (taskid.ne.0) go to 1000
c
	if (taskid.eq.1) write (6,*) 'molin: nom=',nom
!
      if (pstart.le.0) then
c
      do 100 igm=1,ngm
c

#ifdef MOL2
      do 10 i=1,nop2/numtasks
      j1=(i-1)*nompp+1
      j2=j1+nompp-1
      do 20 j=j1,j2
      mpp(j,1,igm)=pp2(i,1)
      mpp(j,2,igm)=pp2(i,2)
      mpp(j,3,igm)=pp2(i,3)
 20   continue
 10   continue
#else
      do 10 i=1,nop/numtasks
      j1=(i-1)*nompp+(igm-1)*nm/ngm+1
      j2=j1+nompp-1
      if (i.eq.1) write (6,*) 'molin: igm, 1st=',igm,j1
      if (i.eq.nop/numtasks) write (6,*) 'molin: igm, last=',igm,j2
      do 20 j=j1,j2
      mpp(j,1)=pp(i,1)
      mpp(j,2)=pp(i,2)
      mpp(j,3)=pp(i,3)
 20   continue
 10   continue
#endif
c
 100  continue
c
      if (taskid.eq.0) write (6,*) 'molin: nop,nom=',nop,nom
c
c calculate average initial position coordinates
c
      avx0(1)=0.
      avx0(2)=0.
      avx0(3)=0.
c
      do j=1,3
        do igm=1,ngm
      do i=1,nom/numtasks/ngm
      avx0(j)=avx0(j)+mpp(i,j,igm)
      end do
        end do
      avx0(j)=avx0(j)/(nom/numtasks)
      end do
c
        if (taskid.le.1.or.taskid.ge.numtasks-2) then
      write (6,*) 'molin: taskid,avx0=',taskid,avx0
        end if

c option for checkpointing
c
	else if (pstart.gt.0) then
c
c if (pstart.gt.0), whole set of particle positions is to be
c read in from disk file
c
      lu=pstart
c
c get inpos in place
c
c
        if (taskid.eq.0) then
c
      open (lu,file='inmpos',form='unformatted',action='read')

c npr = no of particles written per record
c
      read (lu) nompart,npr
        write (6,*) 'nom=',nom,'nompart=',nompart
c
c check nom against the dimensioned limit
c
      if (nom.ne.nompart) then
        write (6,*) 'no of molecules are not consistent'
        write (6,*) 'nom=',nom,'nompart=',nompart
        write (6,*) 'check nop and/or nom'
        stop 'stopped in molin (pstart.gt.0)'
	call MPI_ABORT (MPI_COMM_WORLD,mpierr)
      end if
c
	allocate (mppin(nom/ngm,3))
        write (6,*) 'molin: after allocate'
c
      read (lu) lx1,ly1,lz1
      read (lu) gx1,gy1,gz1
      write (6,*) 'lx1,ly1,lz1=',lx1,ly1,lz1
      write (6,*) 'gx1,gy1,gz1=',gx1,gy1,gz1
c
        end if
c
c
c determine no. of records (nr)
c
        do 40 igm=1,ngm
c
        if (taskid.eq.0) then
        do j=1,3
      read (lu) mppin(:,j)
        end do
        end if
c
	do j=1,3
	call MPI_SCATTER (mppin(1,j),nom/ngm/numtasks,pptype,
     1                    mpp(1,j,igm),nom/ngm/numtasks,pptype,
     1                    0,MPI_COMM_WORLD,mpierr)
	enddo
c
 40     continue

c
      end if
c
 1000 continue
c
c
      if (taskid.eq.1) write (6,*) 'molin: nom,nm=',nom,nm
c
        mpmin=1.e10
        mpmax=-1.e10
        do igm=1,ngm
        do i=1,nom/numtasks/ngm
        mpmin=min(mpmin,mpp(i,2,igm))
        mpmax=max(mpmax,mpp(i,2,igm))
        end do
        end do
        write (6,"('molin: taskid,mpmin,mpmax=',i5,1p,2e12.4)") taskid,mpmin,mpmax
c
	if (taskid.eq.0) write (6,*) 'exit molin'
c
#endif
#endif
c
      return
      end
