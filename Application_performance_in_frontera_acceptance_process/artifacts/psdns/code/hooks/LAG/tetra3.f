      subroutine tetra3 (posn,npdim)
c
c Note: additional arguments added by PK Yeung, Dec 2008
c particle pair positions in the form of tetrahedra with the
c base vertex located at grid points
c
c called from sub. partin (if pstart=-7), task 0 only
c
#ifdef LAG
c
        use compart
	implicit none
        include 'intvars'

	integer npdim
	real(p8) :: posn(npdim,3)
c
	real(b8), allocatable :: postmp(:,:)
c
      real lx,ly,lz,lx1,ly1,lz1
        real sum(3),avx(3)
        real sumt(3),avxt(3)
c
      integer*8 ngrid
	integer ipfst,iplst,isign,ipref,ipx,ipy,ipz,igp,ifst,i1,npu3
	integer i,j,ni,ix1,ix2,iy1,iy2,iz1,iz2
	integer npu,intx,inty,intz,k,index
	real(b8) rx,ry,rz
        integer itask,ii,jj,pseed(2)
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
	write (6,*) 'aborts in sub. tetra'
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
        do itask=0,numtasks-1
        ii=ipid_all(itask)
        jj=jpid_all(itask)
c       write (700,"(7i6)") itask,ii,jj,iist(ii),iisz(ii),jjst(jj),jjsz(jj)
c       write (701,"(7i6)") itask,ii,jj,kist(ii),kisz(ii),kjst(jj),kjsz(jj)
        end do
c
        i=0
        sum=0.

        pseed(1)=17272710
        pseed(2)=pseed(1)*123456
        call random_seed(put=pseed)
c
        do 10 itask=0,numtasks-1
c
        call random_number (postmp)
c
        ii=ipid_all(itask)
        jj=jpid_all(itask)
        sumt=0.
c
        do 20 k=1,ni
        i=i+1
        posn(i,1)=postmp(k,1)*nx+gx(1)
        posn(i,2)=postmp(k,2)*kisz(ii)+kist(ii)-1+gy(1)
        posn(i,3)=postmp(k,3)*kjsz(jj)+kjst(jj)-1+gz(1)
        sum(:)=sum(:)+postmp(k,:)
        sumt(:)=sumt(:)+posn(i,:)
 20   continue
        avxt=sumt/ni
        write (6,"('tetra3: after do 20: itask,avxt=',i4,1p,3e12.4)") itask,avxt
c
 10   continue
c
        avx=sum/npu
        write (6,*) 'tetra3: after do 10: avx=',avx
c
      nop=(1+3*ngp)*nop
      write (6,*) 'tetra3, nop=',nop
c
      if (nop.gt.ndpart) then
      write (6,*) 'no. of particles exceeding dimensioned limit'
      write (6,*) 'stop in tetra: nop,ndpart=',nop,ndpart
      call MPI_ABORT (MPI_COMM_WORLD,mpierr)
      end if
c
c define initial positions of other particles
c
      do 40 igp=1,ngp
c
	write (6,*) 'tetra, igp,dxpr=',igp,dxpr(igp)
      if (dypr(igp).eq.0.) dypr(igp)=dxpr(igp)
      if (dzpr(igp).eq.0.) dzpr(igp)=dxpr(igp)
c
      rx=dxpr(igp)*nx
      ry=dypr(igp)*ny
      rz=dzpr(igp)*nz
c
      ipfst=(igp-1)*3*npu+npu+1
      iplst=ipfst+3*npu-1
      write (6,"('igp,ipfst,iplst=',i3,2i8,'  rx,ry,rz=',1p,3e12.4)")
     1             igp,ipfst,iplst,rx,ry,rz
c
      isign=1
      ipref=0
      do index=ipfst,iplst-2,3
      isign=-isign
      ipref=ipref+1
      posn(index,1)=posn(ipref,1)+isign*rx
      posn(index,2)=posn(ipref,2)
      posn(index,3)=posn(ipref,3)
      posn(index+1,1)=posn(ipref,1)
      posn(index+1,2)=posn(ipref,2)+isign*ry
      posn(index+1,3)=posn(ipref,3)
      posn(index+2,1)=posn(ipref,1)
      posn(index+2,2)=posn(ipref,2)
      posn(index+2,3)=posn(ipref,3)+isign*rz
      end do
c
c
c
 40   continue
c
 210  format ('particle clusters in the form of tetrahedra'/
     1        'no. of clusters(i4)=',i4)
 201  format ('iref,ip1,ip2(3i6)=',3i6)
c
	write (6,*) ' exit tetra3'
#endif
      return
      end
