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
c
	integer itask,np,no_err


      real(p8) avx0(3),avx1(3)
	real(p8) avx
c
        real(p8), allocatable :: pptmp(:,:)
        real(p8), allocatable :: pptmp1(:,:)
        real(p8), allocatable :: pptmp2(:,:)
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
      if (taskid.le.1) write (6,*) 'enter part2in_ar: pstart=',pstart
      lx=xyzl(1)
      ly=xyzl(2)
      lz=xyzl(3)
c
c
c############### this section for 2nd sub-ensemble ########
c
	if (pstart.eq.-9) then
c
        call random_seed (size=seed2)
        allocate(pseed(seed2))
        do i=1,seed2
        pseed(i) = (taskid+1)*123456*i
        enddo
        call random_seed (put=pseed)
c
	npe = ndpart2/numtasks
	allocate (pptmp(npe,3))
c
        call random_number (pptmp)
c
      do 120 i=1,npe
      pp2(1,i)=pptmp(i,1)*nx+gx(1)
      pp2(2,i)=pptmp(i,2)*ny+gy(1)
      pp2(3,i)=pptmp(i,3)*nz+gz(1)
 120  continue
c

c
      avx1(:)=0.
	
      do j=1,3
      do i=1,npe
      avx1(j)=avx1(j)+pp2(j,i)
      end do
      avx1(j)=avx1(j)/npe
      end do
c
	call MPI_ALLREDUCE (avx1,avx0,3,pptype,MPI_SUM,MPI_COMM_WORLD,mpierr)
	avx0=avx0/numtasks
c
	
      if (taskid.eq.0)  write (6,603) pstart,nop2,avx0
 603  format ('part2in_ar: pstart,nop2,avx0=',i3,i10,3f12.4)
c
	deallocate (pptmp)
c
	else if (pstart.gt.0) then
c
        allocate (pptmp1(3,nop2))
        allocate (pptmp2(3,nop2/numtasks))
c
	np=nop2/numtasks 
	if (taskid.eq.1) then
c
        lu=pstart
        open (lu,file='in2pos',form='unformatted', action='read')

c npr = no of particles written per record
c
      read (lu,err=459) nop22,npr
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
      read (lu,err=459,end=459) (pptmp1(:,i),i=i1,i2)
 45   continue
c
	close (lu)
c
c check statistics of initial positions
c
      avx0(:)=0.

      do j=1,3
      do i=1,nop2
      avx0(j)=avx0(j)+pptmp1(j,i)
      end do
      avx0(j)=avx0(j)/nop2
      end do
c
      write (6,602) pstart,nop2,avx0
 602  format ('partin2_ar: pstart,nop2,avx0=',i3,i10,3f12.4)
c
	no_err=0
	goto 450
 459 	no_err=1
	write (6,*) 'ABORT: part2in_ar: problems in reading in2pos file'
 450	continue
c
	end if
c
	call MPI_BCAST (no_err,1,MPI_INTEGER,1,MPI_COMM_WORLD,mpierr)
	if (no_err.gt.0) then
	call MPI_FINALIZE (mpierr)
	stop
	end if
	
c
	call MPI_SCATTER (pptmp1,np*3,pptype,pptmp2,np*3,pptype,
     1                    1,MPI_COMM_WORLD,mpierr)
c
	do i=1,np
	do j=1,3
	pp2(j,i)=pptmp2(j,i)
	end do
	end do

	deallocate (pptmp1,pptmp2)
c
	end if
c
c
#endif
#endif
#endif
      end
