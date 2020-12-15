      subroutine tetra3_local (posn,npdim)
c
c
#ifdef LAG
c
        use compart
	implicit none
        include 'intvars'

	integer npdim,tid,nid
	real(p8) :: posn(npdim,3)
c
	real(p8), allocatable :: postmp(:,:)
        real(p8) xstart,xend,zstart,zend

c
      real lx,ly,lz,lx1,ly1,lz1
        real sum1(3),avx(3)
        real sumt(3),avxt(3)
c
      integer*8 ngrid
	integer ipfst,iplst,isgn,ipref,ipx,ipy,ipz,igp,ifst,i1,npu3
	integer i,j,ni,ix1,ix2,iy1,iy2,iz1,iz2
	integer npu,intx,inty,intz,k,indx
	real(b8) rx,ry,rz
        integer itask,ii,jj,seed2
	integer, allocatable :: pseed(:)
c
      integer isignx(8),isigny(8),isignz(8)
      data isignx/1, 1, 1, 1,-1,-1,-1,-1/
      data isigny/1, 1,-1,-1, 1, 1,-1,-1/
      data isignz/1,-1, 1,-1, 1,-1, 1,-1/
c
      ngrid=nx
      ngrid=ngrid*ny*nz
      if (nop.gt.ngrid) then
        write (6,*) 'no. of tetrahedra vertices (in sub. tetra)'
        write (6,*) 'may not exceed total number of grid points'
        write (6,*) 'reduce nop in lbdata'
	write (6,*) 'aborts in sub. tetra_local'
      call MPI_ABORT (MPI_COMM_WORLD,mpierr)
c
      end if
c
c distribute the reference particles (to be used as
c base vertices of the tetrahedra) randomly
c
      lx=xyzl(1)
      ly=xyzl(2)
      lz=xyzl(3)
c
      npu=nop
        ni=npu/numtasks
	allocate (postmp(ni,3))
c

        call random_seed (size=seed2)
        allocate(pseed(seed2))
        do i=1,seed2
        pseed(i) = (taskid+1)*9897434*i
        enddo
        call random_seed (put=pseed)

!        pseed(1)=17272710
!        pseed(2)=pseed(1)*123456
!        call random_seed(put=pseed)

        i=0
        sum1=0.
c
!        do 10 itask=0,numtasks-1
	itask = taskid
c
        call random_number (postmp)
c
        ii=ipid_all(itask)
        jj=jpid_all(itask)
!        ii = mod(itask,iproc)
!        jj = itask/iproc
        xstart = float(nx)*ii/iproc
        xend = xstart + float(nx)/iproc
        zstart = float(nz)*jj/jproc
        zend = zstart + float(nz)/jproc

        sumt=0.
c
        do 20 k=1,ni
        i=i+1
	if(i.ne.k) write(6,*) 'i and k should be equal'
!        posn(i,1)=postmp(k,1)*nx+gx(1)
!        posn(i,2)=postmp(k,2)*kisz(ii)+kist(ii)-1+gy(1)
!        posn(i,3)=postmp(k,3)*kjsz(jj)+kjst(jj)-1+gz(1)
!
! modfication by Dhawal: splines are in y-slabs, not x-slabs

        posn(i,1)= xstart + postmp(k,1)*(xend-xstart) + gx(1)
        posn(i,2)= postmp(k,2)*ny + gy(1)
        posn(i,3)= zstart + postmp(k,3)*(zend-zstart) + gz(1)

#ifdef DEBUG
        call pptoid( posn(i,1), posn(i,3), tid, nid, 1.)
        if(tid-itask.ne.0) then
        write(6,*) 'tetra3_local, itask, tid=',itask, tid, postmp(k,1:3)
        endif
#endif

        sum1(:)=sum1(:)+postmp(k,:)
        sumt(:)=sumt(:)+posn(i,:)
 20   continue
        avxt=sumt/ni
!        write (6,"('tetra3: after do 20: itask,avxt=',i4,1p,3e12.4)") itask,avxt
c
! 10   continue
c
        avx=sum1/npu
!        write (6,*) 'tetra3_local: after do 10: avx=',avx
c
      nop=(1+3*ngp)*nop
      if(taskid.eq.0) write (6,*) 'tetra3_local, nop=',nop
c
      if (nop.gt.ndpart) then
      write (6,*) 'no. of particles exceeding dimensioned limit'
      write (6,*) 'stop in tetra_local: nop,ndpart=',nop,ndpart
      call MPI_ABORT (MPI_COMM_WORLD,mpierr)
      end if
c
c define initial positions of other particles
c
      do 40 igp=1,ngp
c
	if(taskid.eq.0) then
	write (6,*) 'tetra_local, igp,dxpr=',igp,dxpr(igp)
	endif

      if (dypr(igp).eq.0.) dypr(igp)=dxpr(igp)
      if (dzpr(igp).eq.0.) dzpr(igp)=dxpr(igp)
c
      rx=dxpr(igp)*nx
      ry=dypr(igp)*ny
      rz=dzpr(igp)*nz
c
!      ipfst=(igp-1)*3*npu+npu+1
!      iplst=ipfst+3*npu-1
      ipfst = (igp-1)*3*ni + ni + 1
      iplst = ipfst + 3*ni - 1
	if(taskid.eq.0) then
      write (6,"('igp,ipfst,iplst=',i3,2i8,'  rx,ry,rz=',1p,3e12.4)")
     1             igp,ipfst,iplst,rx,ry,rz

	endif
c
! Change by Dhawal:
! note use of 'indx' and 'isgn' instead of index and isign
! helps ensure avoidance of conflicts with Fortran intrinsics
! called 'index' and 'isign'.

      isgn=1
      ipref=0
      do indx=ipfst,iplst-2,3
      isgn=-isgn
      ipref=ipref+1
      posn(indx,1)=posn(ipref,1)+isgn*rx
      posn(indx,2)=posn(ipref,2)
      posn(indx,3)=posn(ipref,3)
      posn(indx+1,1)=posn(ipref,1)
      posn(indx+1,2)=posn(ipref,2)+isgn*ry
      posn(indx+1,3)=posn(ipref,3)
      posn(indx+2,1)=posn(ipref,1)
      posn(indx+2,2)=posn(ipref,2)
      posn(indx+2,3)=posn(ipref,3)+isgn*rz
      end do
c
c
 40   continue
c
	if(taskid.eq.0) write (6,*) ' exit tetra3_local'
#endif
      return
      end
