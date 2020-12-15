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
	integer ix1,ix2,intx,inty,intz,k
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
	real(p8), allocatable :: ppin(:,:)
c
        real(b8), allocatable :: postmp(:,:)
c
	integer i,j
	real(b8) pmin,pmax
	integer nerr
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
      if (taskid.eq.0) write (6,*) 'enter partin: pstart=',pstart
      lx=xyzl(1)
      ly=xyzl(2)
      lz=xyzl(3)
c
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
	allocate (ppin(ndpart,3))
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

      do 40 k=1,nr
      i1=(k-1)*npr+1
      i2=min0(i1+npr-1,nop)
      read (lu) (ppin(i,1),ppin(i,2),ppin(i,3),i=i1,i2)
 40   continue
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
c
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

	do j=1,3
	call MPI_SCATTER (ppin(1,j),nop/numtasks,pptype,
     1                    pp(1,j),nop/numtasks,pptype,
     1                    0,MPI_COMM_WORLD,mpierr)
	end do
 	deallocate (ppin)
c
	call time_stamp ('partin: before LAGSC2')
c
#ifdef PARTIN_AR
	call part2in_ar
#else
	call part2in 
#endif

#ifdef TEST
c############### this section for 2nd sub-ensemble ########
c (note the ppin array is re-allocated and re-used)
c
#ifdef LAGSC2
c
	allocate (ppin(ndpart2,3))
	
	if (taskid.ne.0) go to 2000
c
c
      if (pstart.le.0) then
c
      npy=nop2/ny
      npe=nop2-npy*ny
      niny=0
c
c if no. of particles is not a multiple of ny,
c then the first npe y-slices will have 1 particle more than npy
c
c set the random no. seed for kipran-th sequence
c
      call ranseq (kipran+1,0)
	
	allocate (postmp(npy,3))
c
      do 110 y=1,ny
c
      ni=npy
      if (y.le.npe) ni=npy+1
      ifst=niny+1
c
      call ranu2 (postmp(1,1),1,ni,1,1,ni)
      call ranu2 (postmp(1,2),1,ni,1,1,ni)
      call ranu2 (postmp(1,3),1,ni,1,1,ni)
c
      do 120 ii=1,ni
      i=niny+ii
#ifdef FEB2411
      ppin(i,1)=postmp(ii,1)*nx+gx(1)
      ppin(i,2)=postmp(ii,2)+gy(y)
      ppin(i,3)=postmp(ii,3)*nz+gz(1)
#else
      ppin(i,1)=postmp(ii,1)*lx+gx(1)
      ppin(i,2)=postmp(ii,2)/ps(2)+gy(y)
      ppin(i,3)=postmp(ii,3)*lz+gz(1)
#endif
 120  continue
c
      niny=niny+ni
c
 110  continue
c
	deallocate (postmp)
c
      else if (pstart.gt.0) then
c
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
	end if
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
	nerr=0
 459  write (6,*) 'abort: error in read in2pos file'
	nerr=1
	call MPI_BCAST (nerr,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
	if (nerr.gt.0) then
	call MPI_ABORT (MPI_COMM_WORLD,mpierr)
	end if

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
#ifdef TEST
c 
c workaround on 3/18/11

      call MPI_BCAST (ppin,nop2*3,mpireal,0,MPI_COMM_WORLD,ierr)
c
	ip1=taskid*nop2/numtasks+1
	ip2=ip1+nop2/numtasks-1
	do j=1,3
	do ip=ip1,ip2
	i=ip-ip1+1
	pp2(i,j)=ppin(ip,j)
	end do
	end do
	
#endif

	

	deallocate (ppin)
#endif
#endif
c###################################################

      
#ifdef LAGSC2
c
      if (taskid.eq.1) write (6,*) 'partin: nop,nop2=',nop,nop2
	call ppminmax ('exit partin:',pp,nop)
	call ppminmax ('exit partin:',pp2,nop2)
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
