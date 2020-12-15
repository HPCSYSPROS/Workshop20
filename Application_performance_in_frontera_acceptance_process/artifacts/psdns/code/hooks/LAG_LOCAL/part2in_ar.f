      subroutine part2in_ar
c

#ifdef LAG
#ifdef LAGSC2
#ifdef PARTIN_AR
! D. Buaria: partin_ar directive added such that
! each mpitask now initializes its own particle posns
! Furthermore, all posns now random in all 3 directions 
! spreading into entire domain

! D Buaria : Feb 4, 2015 : partin_ar still initializes
! its own posns randomly in all 3 directions  
! BUT on local domain of MPI task 


	use compart
	implicit none
	include 'intvars'
c
	common/rancom/iseed,rseed(2)
	integer iseed
	real(b8)  rseed
c
! initial particle positions are set by task 0, then broadcast to
! all other tasks
c
      real lx,ly,lz,lx1,ly1,lz1
        real(p8) px,py,pz
        integer npr,igm,i,j1,j2,j,lu,nchar,nn,nop22,nr,i1,i2,ic
	integer npy,niny,ni,ifst,ii,k,ntasks,itask,tid,nid
        real(b8) gx1,gy1,gz1
	real(p8) xstart,xend,zstart,zend

      real*8 avx0(3)
      data npr/100/
      character*80 string
c
	integer*8 nxnynz
c
	integer  itemp
	integer dmy(3),hms(3)
	integer date_time(8)
c
! Dec 2008: "ppin" is a global array used as temporary buffer
! in setting initial particle positions
c
	real(p8), allocatable :: pp2in(:)
	integer, allocatable :: pp2indin(:)
c
        real(b8), allocatable :: postmp(:,:)
c
        integer seed2
        integer, allocatable :: pseed(:),idisp(:),num_all(:)
c
! if pstart2 is 0, initialise particle positions by uniformly
c distributed random numbers, with (nearly) equal number of
c particles per y-slice
c
! if pstart2 is -ve, particles are initially assigned to be
c at a uniformly spaced subset of grid nodes, (cf. the choice
c of eulerian fixed-points in sub. eultd)
c
! if pstart2 is +ve integer,
c initial particle positions are to be obtained by reading data
c from prepared file (logical unit = pstart2), converting for different
c box sizes and grid shifts if necessary
c
! lx,ly,lz are scaled linear dimensions of the box
c
      if (taskid.eq.0) write (6,*) 'enter part2in_ar: pstart2=',pstart2
      lx=xyzl(1)
      ly=xyzl(2)
      lz=xyzl(3)
c
c
c
      if (pstart2.le.0) then
c
c set the random no. seed for kipran-th sequence
c
! note each taskid sets its own random seed
! make sure kipran+1 = 12
	if(kipran.ne.11) goto 469 

#ifdef LAGSC2_OLD
        call ranseq (kipran+1,0)
#else
        call random_seed (size=seed2)
        allocate(pseed(seed2))
        do i=1,seed2
        pseed(i) = (taskid+1)*123456*i
        enddo
        call random_seed (put=pseed)
#endif
c
	nump2 = ndpart2/numtasks
	allocate (postmp(nump2,3))
c
#ifdef LAGSC2_OLD
      call ranu2 (postmp(1,1),1,nump2,1,1,nump2)
      call ranu2 (postmp(1,2),1,nump2,1,1,nump2)
      call ranu2 (postmp(1,3),1,nump2,1,1,nump2)
#else
	if(taskid.eq.0) write(6,*) 'in part2in_ar, using local init'
        xstart = float(nx)*ipid/iproc
        xend = xstart + float(nx)/iproc
        zstart = float(nz)*jpid/jproc
        zend = zstart + float(nz)/jproc
        call random_number (postmp)
#endif
c
	pp2=0.

        do 120 i=1,nump2


! march 19 2015, D. Buaria
! be careful about loss of significant digits
! when converting itemp to pp2(i,0)
! if itemp is too large then pp2 needs to be double precision 

!	itemp = taskid*nump2 + i
!	pp2(i,0) = 1.d0*itemp

! jul 13 2015, D. Buaria
! introduced pp2ind (an integer array) to keep track of index
! single precision not sufficient for large particle count

	pp2ind(i) = taskid*nump2 + i
        pp2(i,1) = xstart + postmp(i,1)*(xend-xstart) + gx(1)
        pp2(i,2) = postmp(i,2)*ny + gy(1)
        pp2(i,3) = zstart + postmp(i,3)*(zend-zstart) + gz(1)
 120  continue
c
c
      avx0(:)=0.
	
      do j=1,3
      do i=1,nump2
      avx0(j)=avx0(j)+pp2(i,j)
      end do
      avx0(j)=avx0(j)/nump2
      end do
c
	
      if (taskid.eq.0)  write (6,603) pstart2,nop2,avx0
 603  format ('partin: pstart2,nop2,avx0=',i3,i10,3f12.4)
c
	deallocate (postmp)
c


      else if (pstart2.gt.0) then


	pp2=0.
	pp2ind = 0

	call time_stamp ('before readpos 2')
	call readpos_comm_h5 (pp2(1,1), pp2ind, ndpart2, nump2, 2, 0)
!	call readpos_h5 (pp2(1,1), pp2ind, ndpart2, nump2, 2, 0)
	call time_stamp ('after readpos 2')


#ifdef FORTBIO

        allocate (num_all(numtasks),idisp(numtasks),stat=ierr)

	num_all(:) = 0


	if (taskid.eq.0) then
c
	lu=pstart2
	open (lu,file='in2pos',form='unformatted', action='read')

! npr = no of particles written per record
c

      read (lu,err=459) nop22,npr
	if (nop2.ne.nop22) then
     	write (6,*) 'warning: number of particles in in2pos file',
     1              'does not match nop2 in lbdata.f'
	end if
c
c check nop against the dimensioned limit
c
      if (nop2.gt.ndpart2) then
        write (6,*) 'no of particles exceeds dimensioned limit'
        write (6,*) 'nop=',nop,'ndpart=',ndpart
        write (6,*) 'use less particles, or increase ndpart in'
        write (6,*) 'parameter statement, and re-compile'
        stop 'stopped in partin (pstart2.gt.0)'
	call MPI_ABORT (MPI_COMM_WORLD,mpierr)
      end if
c
c
      read (lu,err=459) lx1,ly1,lz1
      read (lu,err=459) gx1,gy1,gz1
c
! determine no. of records (nr)
      nr=nop2/npr
      if (mod(nop2,npr).ne.0) nr=nr+1

	read (lu) ntasks
	if (ntasks.ne.numtasks) then
	write(6,*) 'part2in_ar: restart with same number of MPI tasks'
	stop
	endif

	read (lu) (num_all(i),i=1,numtasks)

	endif

	call MPI_BARRIER (MPI_COMM_WORLD,mpierr)

      call MPI_BCAST (nop2,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST (num_all,numtasks,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

	nump2 = num_all (taskid+1)

        idisp(1) = 0
        do itask=2,numtasks
        idisp(itask) = idisp(itask-1) + num_all(itask-1)
        enddo

	allocate (pp2indin(ndpart2))

	if(taskid.eq.0) then

      do 44 k=1,nr
      i1=(k-1)*npr+1
      i2=min0(i1+npr-1,nop2)
      read (lu,err=459,end=459) (pp2indin(i),i=i1,i2)
 44   continue

	endif

	call MPI_BARRIER (MPI_COMM_WORLD,mpierr)
        call MPI_SCATTERV (pp2indin(1),num_all,idisp,MPI_INTEGER,
     1                    pp2ind(1),nump2,MPI_INTEGER,
     1                    0,MPI_COMM_WORLD,mpierr)


	deallocate (pp2indin)
	allocate (pp2in(ndpart2))

	do j=1,3

	if(taskid.eq.0) then

      do 45 k=1,nr
      i1=(k-1)*npr+1
      i2=min0(i1+npr-1,nop2)
      read (lu,err=459,end=459) (pp2in(i),i=i1,i2)
 45   continue


! check statistics of initial positions
      avx0=0.
      do i=1,nop2
      avx0=avx0+pp2in(i)
      end do
      avx0=avx0/nop2
	write (6,*) 'part2in: pstart2,nop2,j,avx=',pstart2,nop2,j,avx0

	endif

	call MPI_BARRIER (MPI_COMM_WORLD,mpierr)

        call MPI_SCATTERV (pp2in(1),num_all,idisp,pptype,
     1                    pp2(1,j),nump2,pptype,
     1                    0,MPI_COMM_WORLD,mpierr)


	enddo

        if(taskid.eq.0) call fclose1 (lu)
c
	deallocate (pp2in,num_all,idisp)
#endif
! ifdef OLD


 2000	continue

	call time_stamp ('part2in_ar at do 2000')
c

	endif   ! if (pstart....)


	return



 99   continue

        call MPI_ABORT (MPI_COMM_WORLD,mpierr)
	return
c
 459  write (6,*) 'abort: error in read in2pos file'
      go to 99

 469  write (6,*) 'abort: kipran should be 11'
      go to 99
      
c
#endif
#endif
#endif
      end
