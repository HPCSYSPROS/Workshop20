      subroutine part2in_ar
c

#ifdef LAG
#ifdef LAGSC2
#ifdef PARTIN_AR
! D. Buaria: partin_ar directive added such that
! each mpitask now initializes its own particle posns
! Furthermore, all posns now random in all 3 directions 
! spreading into entire domain


	use compart
	implicit none
	include 'intvars'
c
	common/rancom/iseed,rseed(2)
	integer iseed
	real(b8)  rseed
c
c initial particle positions are set by task 0, then broadcast to
c all other tasks
c
      real lx,ly,lz,lx1,ly1,lz1
        real(p8) px,py,pz
        integer npr,igm,i,j1,j2,j,lu,nchar,nn,nop22,nr,i1,i2
	integer npy,npe,niny,ni,ifst,ii,k
        real(b8) gx1,gy1,gz1


      real*8 avx0(3)
      data npr/100/
      character*80 string
c
	integer*8 nxnynz
c
	integer dmy(3),hms(3)
	integer date_time(8)
c
c Dec 2008: "ppin" is a global array used as temporary buffer
c in setting initial particle positions
c
	real(p8), allocatable :: ppin(:,:)
c
        real(b8), allocatable :: postmp(:,:)
c
        integer seed2
        integer, allocatable :: pseed(:)
c
c if pstart is 0, initialise particle positions by uniformly
c distributed random numbers, with (nearly) equal number of
c particles per y-slice
c
c if pstart is -ve, particles are initially assigned to be
c at a uniformly spaced subset of grid nodes, (cf. the choice
c of eulerian fixed-points in sub. eultd)
c
c if pstart is +ve integer,
c initial particle positions are to be obtained by reading data
c from prepared file (logical unit = pstart), converting for different
c box sizes and grid shifts if necessary
c
c lx,ly,lz are scaled linear dimensions of the box
c
      if (taskid.eq.0) write (6,*) 'enter part2in_ar: pstart=',pstart
      lx=xyzl(1)
      ly=xyzl(2)
      lz=xyzl(3)
c
c
c############### this section for 2nd sub-ensemble ########
c (note the ppin array is re-allocated and re-used)
c
c
!	allocate (ppin(ndpart2,3))
!	if (taskid.ne.0) go to 2000
c
c
      if (pstart.le.0) then
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
	npe = ndpart2/numtasks
	allocate (postmp(npe,3))
c
#ifdef LAGSC2_OLD
      call ranu2 (postmp(1,1),1,npe,1,1,npe)
      call ranu2 (postmp(1,2),1,npe,1,1,npe)
      call ranu2 (postmp(1,3),1,npe,1,1,npe)
#else
        call random_number (postmp)
#endif
c
      do 120 i=1,npe
      pp2(i,1)=postmp(i,1)*nx+gx(1)
      pp2(i,2)=postmp(i,2)*ny+gy(1)
      pp2(i,3)=postmp(i,3)*nz+gz(1)

 120  continue
c

c
      avx0(:)=0.
	
      do j=1,3
      do i=1,npe
      avx0(j)=avx0(j)+pp2(i,j)
      end do
      avx0(j)=avx0(j)/npe
      end do
c
	
      if (taskid.eq.0)  write (6,603) pstart,nop2,avx0
 603  format ('partin: pstart,nop2,avx0=',i3,i10,3f12.4)
c
	deallocate (postmp)
c

      else if (pstart.gt.0) then

	allocate (ppin(ndpart2,3))
	if (taskid.ne.0) go to 2000

c
	lu=pstart
c     call fopen1 (lu,'in2pos  ','unformatted')
	open (lu,file='in2pos',form='unformatted', action='read')

c npr = no of particles written per record
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
        stop 'stopped in partin (pstart.gt.0)'
	call MPI_ABORT (MPI_COMM_WORLD,mpierr)
      end if
c
c
      read (lu,err=459) lx1,ly1,lz1
      read (lu,err=459) gx1,gy1,gz1
c
c determine no. of records (nr)
c
      nr=nop2/npr
      if (mod(nop2,npr).ne.0) nr=nr+1
c
      do 45 k=1,nr
      i1=(k-1)*npr+1
      i2=min0(i1+npr-1,nop2)
      read (lu,err=459,end=459) (ppin(i,1),ppin(i,2),ppin(i,3),i=i1,i2)
 45   continue
c
      call fclose1 (lu)
c
!	end if
! moved to after scatter as for pstart <=0 all taskids initialize

c
c check statistics of initial positions
c
      avx0(:)=0.
	
      do j=1,3
      do i=1,nop2
      avx0(j)=avx0(j)+ppin(i,j)
      end do
      avx0(j)=avx0(j)/nop2
      end do
c
	
      if (taskid.eq.0)  write (6,602) pstart,nop2,avx0(1)
 602  format ('partin: pstart,nop2,avx0=',i3,i10,f12.4)
c
 2000	continue
c
	call time_stamp ('partin at do 2000')
      call MPI_BCAST (nop2,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
c
        do j=1,3
        call MPI_SCATTER (ppin(1,j),nop2/numtasks,pptype,
     1                    pp2(1,j),nop2/numtasks,pptype,
     1                    0,MPI_COMM_WORLD,mpierr)
        end do
c
	deallocate (ppin)

	endif

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
