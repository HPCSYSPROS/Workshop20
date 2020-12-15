      subroutine partin_ar
c
#ifdef LAG

	use compart
	implicit none
	include 'intvars'
c
	common/rancom/iseed,rseed(2)
	integer iseed
	real(b8) rseed
c
c initial particle positions are set by task 0, then broadcast to
c all other tasks
c
      real lx,ly,lz,lx1,ly1,lz1
      real*8 avx0(3), avxall(3)
	real(b8) gx1,gy1,gz1
	real(p8) px,py,pz
	integer i1,i2,nr,iy1,iy2,iz1,iz2,index,lu
	integer ix1,ix2,intx,inty,intz,k,ic, ip,igp
	integer npy,npe,niny,n1,ifst,ii,npu,npr,ni
	integer subsid, igrp, npset

      data npr/100/
      character*80 string
c
	integer*8 nxnynz
        integer no_err, ntotal, npmax, npmin
c
	integer dmy(3),hms(3)
	integer date_time(8)
c
c Dec 2008: "ppin" is a global array used as temporary buffer
c in setting initial particle positions
c
	real(p8), allocatable :: ppin(:,:), pptmp(:,:)
	integer, allocatable :: ppindin(:), ppindtmp(:)
c
	integer ntasks,itask

	integer, allocatable :: num_all(:), idisp(:)
	integer, allocatable :: scount(:), sdisp(:)
	integer, allocatable :: rcount(:), rdisp(:)
	integer, allocatable :: icount(:)
c
	integer i,j, tid, nid, niderr
	real(b8) pmin,pmax
	integer nerr

	integer iread, mread, nread, relay, jsource
	integer mpistatus (MPI_STATUS_SIZE)

c
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

      if (taskid.eq.0) write (6,*) 'enter partin_ar: pstart=',pstart

	if(pstart.le.0 .and. pstart.ne.-10) then
	write(6,*) 'pstart must be -10 for non restart cases'
	stop
	endif

c
	if (nop.gt.0) then



	if(pstart.eq.-10) then

      if (nop.gt.ndpart) then
      if (taskid.eq.0) then
        write (6,*) 'no of particles exceeds dimensioned limit'
        write (6,*) 'nop=',nop,'ndpart=',ndpart
        write (6,*) 'use fewer particles, or increase ndpart in'
        write (6,*) 'parameter statement, and re-compile'
        write (6,*) 'stopped in partin_ar (pstart=-1)'
	end if
	stop
	end if


      lx=xyzl(1)
      ly=xyzl(2)
      lz=xyzl(3)
c
	allocate (num_all(numtasks), idisp(numtasks))

	allocate (scount(numtasks), sdisp(numtasks))
	allocate (rcount(numtasks), rdisp(numtasks))
	allocate (icount(numtasks))
	allocate (ppin(ndpart/numtasks,3))
	allocate (pptmp(ndpart/numtasks,3))
	allocate (ppindin(ndpart/numtasks))
	allocate (ppindtmp(ndpart/numtasks))

c
	npu = ndpart/(1+3*ngp)
	ni = npu/numtasks
	npset = ndpart/nsubset

	nump1 = ndpart/numtasks
! this value of nump1 will change every step

       call tetra3_local (ppin,ndpart/numtasks)

! 3/18/11: force ndpart equal to nop
! (otherwise, because some arrays may be dimensioned as ndpart
! elements, there is risk of array trouble)
	ndpart = nop

! D Buaria, Sep 7, 2015,
! depending on which subroutine is used for IO,
! we need use different indexing for pp1 particles
! to ensure consisten IO for postprocessing codes
#ifdef ALLG
! this indexing is used when using pop1_write_local
	do i=1,ni
	ppindin(i) = taskid*ni + i
	enddo
	do igp=0,ngp-1
	do i=1,3*ni
	ppindin(ni + igp*3*ni + i) = npu + igp*3*npu + taskid*3*ni+i 
	enddo
	enddo
#else
! this indexing is used when using pop1_write_local2
	subsid = mod(taskid,numtasks/nsubset)
	igrp = taskid/(numtasks/nsubset)
	do i=1,ni
	ppindin(i) = igrp*npset + subsid*ni + i
	enddo
	do igp=0,ngp-1
	do i=1,3*ni
	ppindin(ni + igp*3*ni + i) = igrp*npset + npu/nsubset 
     &                 + igp*3*npu/nsubset + subsid*3*ni + i 
	enddo
	enddo
#endif

!	if(taskid.lt.1000) then
!	do i=1,nump1
!	write(7000+taskid,*) i, ppindin(i)
!	enddo
!	close(7000+taskid)
!	endif

#ifdef DXPR_INIT
! D. Buaria, Sep 15 2015
! If the dxpr for tetrads is less than grid spacing,
! the vertices of tetrads cannot go beyond neighboring tasks
! Thus can use the update_part subroutine to shuffle them around,
! else need to do alltoall 

	pp(1:nump1, 1:3) = ppin(1:nump1,1:3)
	pp1ind(1:nump1) = ppindin(1:nump1)
	call update_part (pp(1,1), pp1ind(1), ndpart, nump1, 3)
#else
	scount = 0
	rcount = 0
	icount = 0

	do i=1,ndpart/numtasks
	call pptoid( ppin(i,1), ppin(i,3), tid, nid, 1.)
	tid = tid + 1
	if (tid.lt.1.or.tid.gt.numtasks) then
	write(6,*) 'problem with tid=',tid, taskid, ppin(i,1:3)
	endif
	scount(tid) = scount(tid) + 1
	enddo

	sdisp(1) = 0
	do itask=2,numtasks
	sdisp(itask) = sdisp(itask-1) + scount(itask-1)
	enddo

	do i=1,ndpart/numtasks
	call pptoid( ppin(i,1), ppin(i,3), tid, nid, 1.)
	tid = tid+1
	if (tid.lt.1.or.tid.gt.numtasks) then
	write(6,*) 'problem with tid=',tid, taskid, ppin(i,1:3)
	endif
	icount(tid) = icount(tid) + 1
	pptmp(sdisp(tid)+icount(tid),1:3) = ppin(i,1:3)
	ppindtmp(sdisp(tid)+icount(tid)) = ppindin(i)
	enddo

	deallocate (ppin, ppindin)

      call MPI_ALLTOALL (scount, 1, MPI_INTEGER,
     &                  rcount, 1, MPI_INTEGER, MPI_COMM_WORLD, mpierr)
	rdisp(1) = 0
	do itask=2,numtasks
	rdisp(itask) = rdisp(itask-1) + rcount(itask-1)
	enddo


      call mpi_alltoallv (ppindtmp(1), scount, sdisp, MPI_INTEGER, 
     &      pp1ind(1), rcount, rdisp, MPI_INTEGER, MPI_COMM_WORLD, mpierr)

	do j=1,3

      call mpi_alltoallv (pptmp(1,j), scount, sdisp, pptype,
     &      pp(1,j), rcount, rdisp, pptype, MPI_COMM_WORLD, mpierr)

	enddo


	nump1 = sum(rcount)

#endif

        call MPI_ALLREDUCE (nump1, ntotal, 1, MPI_INTEGER, MPI_SUM,
     1                   MPI_COMM_WORLD, mpierr)
        call MPI_ALLREDUCE (nump1,npmax,1,MPI_INTEGER,MPI_MAX,
     1                   MPI_COMM_WORLD,mpierr)
        call MPI_ALLREDUCE (nump1,npmin,1,MPI_INTEGER,MPI_MIN,
     1                   MPI_COMM_WORLD,mpierr)


	if(taskid.eq.0) then
	write(6,*) 'partin_ar, max, min, avg nump1=',npmax, npmin, nop/numtasks
	endif

	if(ntotal.ne.nop) then
	write(6,*) 'partin_ar, ntotal .ne. nop', ntotal, nop
	stop
	endif        


	deallocate (pptmp, ppindtmp)

! check nop against the dimensioned limit
c
      if (nop.gt.ndpart) then
        write (6,*) 'no of particles exceeds dimensioned limit'
        write (6,*) 'nop=',nop,'ndpart=',ndpart
        write (6,*) 'use fewer particles, or increase ndpart in'
        write (6,*) 'parameter statement, and re-compile'
        write (6,*) 'stopped in partin_ar (pstart=-1)'
	stop
      end if
c
c
! for shear flow, force the mean of initial particle positions
! in the y-direction to be at the 1st y-plane
! (so that the mean velocity could be near zero when averaged
! over all particles)
c

      if (shear.gt.0..and.pstart.le.0) then
      do i=1,nump1
      pp(i,2)=pp(i,2)-(ny/2)/ps(2)
      end do
      end if
c

! calculate average initial particle position coordinates
c
	elseif (pstart.gt.0) then

	pp=0.
	pp1ind=0

	call time_stamp ('before readpos 1')
	call MPI_BARRIER (MPI_COMM_WORLD, mpierr)

! relay scheme 

       if (jproc.le.1024) mread=jproc
       if (mod(jproc,1024).eq.0) mread=1024

        nread=jproc/mread
        iread=jpid/mread

        if (iread.eq.0) then

        relay=taskid

	call readpos_comm_h5 (pp(1,1), pp1ind, ndpart, nump1, 1, 0)
        else
        jsource=taskid-mread*iproc

        call MPI_RECV (relay,1,MPI_INTEGER,jsource,jsource,
     1                 MPI_COMM_WORLD,mpistatus,mpierr)

	call readpos_comm_h5 (pp(1,1), pp1ind, ndpart, nump1, 1, 0)
        end if

        if (iread.lt.nread-1) then
        call MPI_SSEND (relay,1,MPI_INTEGER,taskid+mread*iproc,taskid,
     1                MPI_COMM_WORLD,mpierr)
        end if
c



	call time_stamp ('after readpos 1')
	call MPI_BARRIER (MPI_COMM_WORLD, mpierr)

	endif ! if (pstart.eq.-10)

      avx0(1)=0.
      avx0(2)=0.
      avx0(3)=0.
c
      do j=1,3
      do i=1,nump1
      avx0(j)=avx0(j)+pp(i,j)
      end do
      end do

        call MPI_ALLREDUCE (avx0, avxall, 3, MPI_REAL8, MPI_SUM,
     1                   MPI_COMM_WORLD, mpierr)

      avxall(:)=avxall(:)/nop

      if(taskid.eq.0) write (6,*) 'partin_ar: pstart,nop,avx0=',pstart,nop,avxall

c
!        call MPI_BCAST (no_err,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!        if (no_err.gt.0) then
!        if (taskid.eq.0) write (6,*) 'abort in partin_ar: no_err=',no_err
!        stop
!        end if 
        
      if (taskid.eq.0) then
      call date_and_time (values=date_time)
      dmy(1:3)=date_time(3:1:-1)
      hms(1:3)=date_time(5:7)
      write (6,701) dmy(2),dmy(1),dmy(3),hms
 701  format ('partin_ar at do 1000:  date & time is  ',i2,'/',i2,
     1        '/',i4,2x,i2,':',i2,':',i2)
      end if


	call time_stamp ('partin_ar: before LAGSC2')


	endif ! if (nop.gt.0) 


!############### this section for 2nd sub-ensemble ########

	if (nop2.gt.0) then

#ifdef PARTIN_AR
	call part2in_ar
#else
	call part2in 
#endif

	endif




c###################################################

      
#ifdef LAGSC2
c
      if (taskid.eq.1) write (6,*) 'partin_ar: nop,nop2=',nop,nop2
	call ppminmax ('exit partin_ar:',pp(1,1),nop)
	call ppminmax ('exit partin_ar:',pp2(1,1),nop2)
c
#else
      if (taskid.eq.1) write (6,*) 'partin_ar: nop=',nop
#endif
c
#ifdef MOL
	if (nom.gt.0) call molin
#endif
c
#ifdef MATP
	call mpin
#endif
c
      call date_and_time (values=date_time)
      dmy(1:3)=date_time(3:1:-1)
      hms(1:3)=date_time(5:7)
      if (taskid.eq.0)
     1   write (6,702) dmy(2),dmy(1),dmy(3),hms
 702  format ('partin_ar at    exit:  date & time is  ',i2,'/',i2,
     1        '/',i4,2x,i2,':',i2,':',i2)

      return
c
c
#endif
      end
