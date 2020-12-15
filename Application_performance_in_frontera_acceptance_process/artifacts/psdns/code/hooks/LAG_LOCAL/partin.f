      subroutine partin
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
      real*8 avx0(3)
	real(b8) gx1,gy1,gz1
	real(p8) px,py,pz
	integer i1,i2,nr,iy1,iy2,iz1,iz2,index,lu
	integer ix1,ix2,intx,inty,intz,k,ic
	integer npy,npe,niny,n1,ifst,ii,npu,npr,ni

      data npr/100/
      character*80 string
c
	integer*8 nxnynz
        integer no_err
c
	integer dmy(3),hms(3)
	integer date_time(8)
c
c Dec 2008: "ppin" is a global array used as temporary buffer
c in setting initial particle positions
c
	real(p8), allocatable :: ppin(:,:), pptmp(:,:)
c
	integer ntasks,itask
        real(b8), allocatable :: postmp(:,:)
	integer, allocatable :: num_all(:), ppindtmp(:,:)
	integer, allocatable :: ppindin(:), idisp(:)
c
	integer i,j, tid, nid, niderr
	real(b8) pmin,pmax
	integer nerr
c
	integer, allocatable :: sendcount(:),num_allb(:)
	integer ib,nb,ip1
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

      if (taskid.eq.0) write (6,*) 'enter partin: pstart=',pstart

      lx=xyzl(1)
      ly=xyzl(2)
      lz=xyzl(3)
c
	if (nop.gt.0) then

      if (nop.gt.ndpart) then
      if (taskid.eq.0) then
        write (6,*) 'no of particles exceeds dimensioned limit'
        write (6,*) 'nop=',nop,'ndpart=',ndpart
        write (6,*) 'use fewer particles, or increase ndpart in'
        write (6,*) 'parameter statement, and re-compile'
        write (6,*) 'stopped in partin (pstart=-1)'
	end if
	stop
	end if
c
	allocate (num_all(numtasks),idisp(numtasks))
	allocate (ppin(ndpart,3))
	allocate (ppindin(ndpart))

	allocate (num_allb(numtasks))
c
      if (taskid.ne.0) go to 1000
c
      if (pstart.eq.0) then
c
c subroutine ranu2 is called to generate uniformly distributed random no
c in (0,1) which are used to compute particle position coordinates
c
c in this case, no. of particles (nop) is specified in block data
c
c
      npy=nop/ny
      npe=nop-npy*ny
      niny=0
c
c if no. of particles is not a multiple of ny,
c then the first npe y-slices will have 1 particle more than npy
c
c set the random no. seed for kipran-th sequence
c
      call ranseq (kipran,0)
c
      do 10 y=1,ny
c
      ni=npy
      if (y.le.npe) ni=npy+1
      ifst=niny+1
c
      call ranu2 (ppin(ifst,1),1,ni,1,1,ni)
      call ranu2 (ppin(ifst,2),1,ni,1,1,ni)
      call ranu2 (ppin(ifst,3),1,ni,1,1,ni)
c
      do 20 ii=1,ni
      i=niny+ii
      ppin(i,1)=ppin(i,1)*lx+gx(1)
      ppin(i,2)=ppin(i,2)/ps(2)+gy(y)
      ppin(i,3)=ppin(i,3)*lz+gz(1)
 20   continue
c
      niny=niny+ni
c
 10   continue
c
      else if (pstart.eq.-1) then
c
      npu=nop
      nop=nop**3
c
c nop may not exceed the total number of grid points
c
	nxnynz=nx*ny*nz
      if (nop.gt.nxnynz) then
        write (6,*) 'no. of particles (initialized as a cubic lattice)'
        write (6,*) 'may not exceed total number of grid points'
        write (6,*) 'reduce nop in lbdata'
        no_err=1
        write (6,*) 'stopped in partin (pstart=-1)'
        go to 90
      end if
c
c define the grid locations
c
      intx=nx/npu
      inty=ny/npu
      intz=nz/npu
      if (mod(nx,npu).ne.0) intx=intx+1
      if (mod(ny,npu).ne.0) inty=inty+1
      if (mod(nz,npu).ne.0) intz=intz+1
c
      ix1=nx/2+1-npu/2*intx+intx/2
      ix2=nx/2+1+npu/2*intx-intx/2
      iy1=ny/2+1-npu/2*inty+inty/2
      iy2=ny/2+1+npu/2*inty-inty/2
      iz1=nz/2+1-npu/2*intz+intz/2
      iz2=nz/2+1+npu/2*intz-intz/2
c
c ensure that ix1,ix2, etc., are within the range from 1 to nx
c
      ix1=max0(ix1,1)
      ix2=min0(ix2,nx)
      iy1=max0(iy1,1)
      iy2=min0(iy2,ny)
      iz1=max0(iz1,1)
      iz2=min0(iz2,nz)
c
      index=0
      do 30 k=iz1,iz2,intz
      do 30 j=iy1,iy2,inty
      do 30 i=ix1,ix2,intx
      index=index+1
      ppin(index,1)=gx(i)
      ppin(index,2)=gy(j)
      ppin(index,3)=gz(k)
 30   continue
c
      if (index.ne.nop) then
      write (6,*) 'option pstart=-1 in sub. partin'
      write (6,*) 'index=',index,'  not equal to nop=',nop
      write (6,*) 'initialization mechanism fails'
      write (6,*) 'program stopped in partin'
	call MPI_ABORT (MPI_COMM_WORLD,mpierr)
      end if
c
      else if (pstart.eq.-2) then
c
      npu=nop
      nop=2.*nop**3
c
c check nop against the dimensioned limit
c
      if (nop.gt.ndpart) then
        write (6,*) 'no of particles exceeds dimensioned limit'
        write (6,*) 'nop=',nop,'ndpart=',ndpart
        write (6,*) 'use fewer particles, or increase ndpart in'
        write (6,*) 'parameter statement, and re-compile'
        write (6,*) 'stopped in partin (pstart=-1)'
	call MPI_ABORT (MPI_COMM_WORLD,mpierr)
      end if
c
c nop may not exceed the total number of grid points
c
	nxnynz=nx*ny*nz
      if (nop.gt.nxnynz) then
        write (6,*) 'no. of particles (initialized as a cubic lattice)'
        write (6,*) 'may not exceed total number of grid points'
        write (6,*) 'reduce nop in lbdata'
        write (6,*)  'stopped in partin (pstart=-2)'
	call MPI_ABORT (MPI_COMM_WORLD,mpierr)
      end if
c
c define the grid locations
c
      intx=nx/npu
      inty=ny/npu
      intz=nz/npu
      if (mod(nx,npu).ne.0) intx=intx+1
      if (mod(ny,npu).ne.0) inty=inty+1
      if (mod(nz,npu).ne.0) intz=intz+1
c
      ix1=nx/2+1-npu/2*intx+intx/2
      ix2=nx/2+1+npu/2*intx-intx/2
      iy1=ny/2+1-npu/2*inty+inty/2
      iy2=ny/2+1+npu/2*inty-inty/2
      iz1=nz/2+1-npu/2*intz+intz/2
      iz2=nz/2+1+npu/2*intz-intz/2
c
c ensure that ix1,ix2, etc., are within the range from 1 to nx
c
      ix1=max0(ix1,1)
      ix2=min0(ix2,nx)
      iy1=max0(iy1,1)
      iy2=min0(iy2,ny)
      iz1=max0(iz1,1)
      iz2=min0(iz2,nz)
c
c to place first grid (same as pstart=-1)
c
      index=0
      do 33 k=iz1,iz2,intz
      do 33 j=iy1,iy2,inty
      do 33 i=ix1,ix2,intx
      index=index+1
#ifdef DEc2108
      ppin(index,1)=gx(i)
      ppin(index,2)=gy(j)
      ppin(index,3)=gz(k)
#else
      pp(index,1)=gx(i)
      pp(index,2)=gy(j)
      pp(index,3)=gz(k)
#endif
 33   continue
c
c to place second grid (staggered relative to the first)
c
      ix1=ix1-intx/2
      ix2=ix2-intx/2
      iy1=iy1-inty/2
      iy2=iy2-inty/2
      iz1=iz1-intz/2
      iz2=iz2-intz/2
c
      do 35 k=iz1,iz2,intz
      do 35 j=iy1,iy2,inty
      do 35 i=ix1,ix2,intx
      index=index+1
      ppin(index,1)=gx(i)
      ppin(index,2)=gy(j)
      ppin(index,3)=gz(k)
 35   continue
c
      if (index.ne.nop) then
      write (6,*) 'option pstart=-2 in sub. partin'
      write (6,*) 'index=',index,'  not equal to nop=',nop
      write (6,*) 'initialization mechanism fails'
      write (6,*) 'program stopped in partin'
	call MPI_ABORT (MPI_COMM_WORLD,mpierr)
      end if
c
c if (pstart.gt.0), whole set of particle positions is to be
c read in from disk file
c
      else if (pstart.gt.0) then
c
      lu=pstart
c
c      call fopen1 (lu,'inpos  ','unformatted')
	open (lu,file='inpos',form='unformatted', action='read')

c npr = no of particles written per record
c
      read (lu) nop,npr
c
c because of dynamic array memory allocation, results will go wrong if
c nop in incoming inpos file is not the same as that in input.lag file
c
      if (nop.ne.ndpart) then
        write (6,*) 'PARTIN: Stop due to inconsistency in nop! (pstart.gt.0)'
        write (6,*) 'PARTIN: change nop in input.lag to ',nop
        no_err=2
	go to 90
      end if
c
c
      read (lu) lx1,ly1,lz1
      read (lu) gx1,gy1,gz1
      write (6,*) 'nop,npr=',nop,npr
      write (6,*) 'gx1,gy1,gz1=',gx1,gy1,gz1
c
c determine no. of records (nr)
c
      nr=nop/npr
      if (mod(nop,npr).ne.0) nr=nr+1
c
        read (lu) ntasks
        if (ntasks.ne.numtasks) then
        write(6,*) 'partin: restart with same number of MPI tasks'
        stop
        endif

        read (lu) (num_all(i),i=1,numtasks)


      do 40 k=1,nr
      i1=(k-1)*npr+1
      i2=min0(i1+npr-1,nop)
      read (lu) (ppindin(i),i=i1,i2)
 40   continue

	do j=1,3
      do 41 k=1,nr
      i1=(k-1)*npr+1
      i2=min0(i1+npr-1,nop)
      read (lu) (ppin(i,j),i=i1,i2)
 41   continue
	enddo
      call fclose1 (lu)
c
c determine whether conversion is necessary
c
      if (lx.ne.lx1) go to 55
      if (ly.ne.ly1) go to 55
      if (lz.ne.lz1) go to 55
c
c note there may be some round-off errors for restart from
c chkptou/ranout from another machine
c
      if (abs(gx(1)-gx1).gt.1.e-5) go to 55
      if (abs(gy(1)-gy1).gt.1.e-5) go to 55
      if (abs(gz(1)-gz1).gt.1.e-5) go to 55
c     if (gx(1).ne.gx1) go to 55
c     if (gy(1).ne.gy1) go to 55
c     if (gz(1).ne.gz1) go to 55
c
      go to 66
c
c for a checkpoint continuation run, the box size and grid shifts
c must match those of the incoming data
c
 55   if (istart.ne.0.and.pstart.gt.0) then
        write (6,*) 'no conversion of particle coordinates should'
        write (6,*) 'have been needed for checkpoint continuation run'
        write (6,6091) lx,ly,lz,lx1,ly1,lz1
        write (6,6092) gx(1),gy(1),gz(1),gx1,gy1,gz1
        write (6,*) 'stopped in partin (pstart.gt.0)'
	call MPI_ABORT (MPI_COMM_WORLD,mpierr)
 6091 format ('lx,lx1,=...',6f10.5)
 6092 format ('gx,gx1,=...',6f10.5)
      end if
c
      do 50 i=1,nop
c
c first convert to coordinates normalised by length of box
c (i.e. between 0 and 1)
c then express particle positions in new coord system
c
      px=(ppin(i,1)-gx1)/lx1
      py=(ppin(i,2)-gy1)/ly1
      pz=(ppin(i,3)-gz1)/lz1
      ppin(i,1)=px*lx+gx(1)
      ppin(i,2)=py*ly+gy(1)
      ppin(i,3)=pz*lz+gz(1)
c
 50   continue
c
 66   continue
c
c
c
      else if (pstart.eq.-7) then
      call tetra (1,ppin,ndpart)
      else if (pstart.eq.-8) then
      call tetra (2,ppin,ndpart)
      else if (pstart.eq.-9) then
      if (shear.gt.0.) then
        no_err=3
        go to 90
        end if
      call tetra (3,ppin,ndpart)
c

      end if

        if (pstart.eq.-10) call tetra3 (ppin,ndpart)
c
c check nop against the dimensioned limit
c
      if (nop.gt.ndpart) then
        write (6,*) 'no of particles exceeds dimensioned limit'
        write (6,*) 'nop=',nop,'ndpart=',ndpart
        write (6,*) 'use fewer particles, or increase ndpart in'
        write (6,*) 'parameter statement, and re-compile'
        write (6,*) 'stopped in partin (pstart=-1)'
        no_err=4
	go to 90
      end if
c
c
      write (6,*) 'pstart, no. of particles=',pstart,nop
c
c for shear flow, force the mean of initial particle positions
c in the y-direction to be at the 1st y-plane
c (so that the mean velocity could be near zero when averaged
c over all particles)
c

      if (shear.gt.0..and.pstart.le.0) then
      do i=1,nop
      ppin(i,2)=ppin(i,2)-(ny/2)/ps(2)
      end do
      end if
c

c calculate average initial particle position coordinates
c
      avx0(1)=0.
      avx0(2)=0.
      avx0(3)=0.
c
      do j=1,3
      do i=1,nop
      avx0(j)=avx0(j)+ppin(i,j)
      end do
      avx0(j)=avx0(j)/nop
      end do
c
      write (6,*) 'partin: pstart,nop,avx0=',pstart,nop,avx0
c
        no_err=0
        go to 1000
c
 90    continue
c
 1000 continue
c
        call MPI_BCAST (no_err,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        if (no_err.gt.0) then
        if (taskid.eq.0) write (6,*) 'abort in partin: no_err=',no_err
        stop
        end if 
        
      if (taskid.eq.0) then
      call date_and_time (values=date_time)
      dmy(1:3)=date_time(3:1:-1)
      hms(1:3)=date_time(5:7)
      write (6,701) dmy(2),dmy(1),dmy(3),hms
 701  format ('partin at do 1000:  date & time is  ',i2,'/',i2,
     1        '/',i4,2x,i2,':',i2,':',i2)
      end if
      
c some of the pstart options may have changed nop to its true value
c
      call MPI_BCAST (nop,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
c
c 3/18/11: force ndpart equal to nop
c (otherwise, because some arrays may be dimensioned as ndpart
c elements, there is risk of array trouble)
c
	ndpart=nop


	pp = 0.

	if(pstart.gt.0) then

	call MPI_BCAST (num_all,numtasks,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
	nump1 = num_all(taskid+1)

        idisp(1) = 0
        do itask=2,numtasks
        idisp(itask) = idisp(itask-1) + num_all(itask-1)
        enddo

        call MPI_SCATTERV (ppindin(1),num_all,idisp,MPI_INTEGER,
     1                    pp1ind(1),nump1,MPI_INTEGER,
     1                    0,MPI_COMM_WORLD,mpierr)


	do j=1,3
        call MPI_SCATTERV (ppin(1,j),num_all,idisp,pptype,
     1                    pp(1,j),nump1,pptype,
     1                    0,MPI_COMM_WORLD,mpierr)
	enddo

	deallocate (ppin, ppindin)


	else

	allocate (sendcount(numtasks))
	nb=2
	ip1=1
c
	allocate (pptmp(nppad*ndpart/numtasks/nb,numtasks))
	allocate (ppindtmp(nppad*ndpart/numtasks/nb,numtasks))

	pptmp = 0.
	ppindtmp = 0
	num_all = 0
	niderr = 0


	do 200 ib=1,nb

	i1=(ib-1)*nop/nb+1
	i2=ib*nop/nb
	num_allb = 0
	
	if (taskid.eq.0) then

	do i=i1,i2
	call pptoid( ppin(i,1), ppin(i,3), tid, nid, 1.)
	tid = tid + 1
	if (tid.lt.1.or.tid.gt.numtasks) then
	write(6,*) 'problem with tid=',tid, taskid, ppin(i,1:3)
	endif
	num_allb(tid) = num_allb(tid) + 1
	ppindtmp (num_allb(tid), tid) = i
	enddo

	do i=1,numtasks
	write(7003,*) i, num_allb(i)
	end do

	num_all = num_all + num_allb
	write(6,*) 'partin, nop, max, min, avg=',maxval(num_all),
     1                 minval(num_all), nop/numtasks

 	endif

	call MPI_BCAST (niderr,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
	if(niderr.eq.1) then
	if(taskid.eq.0) write(6,*)'error in partin, num_all out of bounds'
	stop
	endif

	call MPI_BCAST (num_allb,numtasks,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
	call MPI_BCAST (num_all,numtasks,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)

	do itask=1,numtasks
 	idisp(itask)=(itask-1)*nppad*nop/numtasks/nb
	sendcount(itask)=num_allb(itask)
	end do
	call MPI_SCATTERV (ppindtmp(1,1),sendcount,idisp,MPI_INTEGER,
     1                    pp1ind(ip1),num_allb(taskid+1),MPI_INTEGER,
     1                    0,MPI_COMM_WORLD,mpierr)

! 9/1/2015: to save on memory overhead of the buffer arrays involved,
! do the following 1 coordinate component at a time.

	do 100 j=1,3

	if(taskid.eq.0) then
	 do i=1,numtasks
	   do ii=1,num_allb(i)
	   pptmp(ii,i) = ppin(ppindtmp(ii,i),j)
	   enddo
	 enddo
	endif

c	call MPI_BARRIER (MPI_COMM_WORLD, mpierr)

	call MPI_SCATTERV (pptmp(1,1),sendcount,idisp,pptype,
     1                    pp(ip1,j),num_allb(taskid+1),pptype,
     1                    0,MPI_COMM_WORLD,mpierr)

c	call MPI_BARRIER (MPI_COMM_WORLD, mpierr)

 100	continue

	ip1=ip1+num_allb(taskid+1)

 200    continue

	if (taskid.eq.27) then
	do i=1,num_all(taskid+1)
	write (7000+taskid,*) i,pp1ind(i)
	end do
	end if
	nump1 = num_all(taskid+1)

	do i=1,num_all(taskid+1)
	if (pp1ind(i).eq.0) then
	write (6,*) 'pp1ind=0:',taskid,i
	end if
	end do

 	deallocate (ppin)
	deallocate (pptmp, ppindtmp)

	if (taskid.eq.0) then
	ic=0
	do i=1,numtasks
	write(8003,*) i, num_all(i)
	if(num_all(i).gt.nppad*ndpart/numtasks) niderr = 1
	if(num_all(i).le.0) niderr = 1
#ifdef DEBUG
	if(num_all(i).ne.nop/numtasks) then
	ic=ic+1
	write(8003+ic,*) 'task=',i
	do ii=1,num_all(i)
	write(8003+ic,*) ii, pptmp(ii,1:3,i)
	enddo
	endif
#endif
	enddo
	close(8003)
	endif

	endif

	call MPI_FINALIZE (mpierr)
	stop 'stop in partin'
	call time_stamp ('partin: before LAGSC2')


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
      if (taskid.eq.1) write (6,*) 'partin: nop,nop2=',nop,nop2
	call ppminmax ('exit partin:',pp(1,1),nop)
	call ppminmax ('exit partin:',pp2(1,1),nop2)
c
#else
      if (taskid.eq.1) write (6,*) 'partin: nop=',nop
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
 702  format ('partin at    exit:  date & time is  ',i2,'/',i2,
     1        '/',i4,2x,i2,':',i2,':',i2)

      return
c
c
#endif
      end
