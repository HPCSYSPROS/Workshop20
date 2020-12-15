      subroutine tetra (iopt,posn,npdim)
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
c
      integer*8 ngrid
	integer ipfst,iplst,isign,ipref,ipx,ipy,ipz,igp,ifst,i1,npu3
	integer i,j,npy,npe,niny,ni,ii,ix1,ix2,iy1,iy2,iz1,iz2
	integer npu,intx,inty,intz,k,iopt,index
	real(b8) rx,ry,rz
c
      character*7 fnx,fny,fnz
      character*1 numc(9)
      data numc/'1','2','3','4','5','6','7','8','9'/
c
      integer isignx(8),isigny(8),isignz(8)
      data isignx/1, 1, 1, 1,-1,-1,-1,-1/
      data isigny/1, 1,-1,-1, 1, 1,-1,-1/
      data isignz/1,-1, 1,-1, 1,-1, 1,-1/
c
	write (6,*) ' enter tetra: iopt=',iopt
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
      if (iopt.le.2) then
c
      npu=nop
      nop=iopt*(1+3*ngp)*npu**3
c
c define the grid locations for the vertices
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
c (they should be, even without the action below)
c
      ix1=max0(ix1,1)
      ix2=min0(ix2,nx)
      iy1=max0(iy1,1)
      iy2=min0(iy2,ny)
      iz1=max0(iz1,1)
      iz2=min0(iz2,nz)
c
      if (iopt.eq.2) then
      ix1=ix1-intx/4
      ix2=ix2-intx/4
      iy1=iy1-inty/4
      iy2=iy2-inty/4
      iz1=iz1-intz/4
      iz2=iz2-intz/4
      end if
c
c locate the reference group of particles at grid points
c
      write (6,*) 'tetra (iopt=2): ix1,ix2=',ix1,ix2
      write (6,*) 'tetra (iopt=2): iy1,iy2=',iy1,iy2
      write (6,*) 'tetra (iopt=2): iz1,iz2=',iz1,iz2
      index=0
      do 30 k=iz1,iz2,intz
      do 30 j=iy1,iy2,inty
      do 30 i=ix1,ix2,intx
      index=index+1
      posn(index,1)=gx(i)
      posn(index,2)=gy(j)
      posn(index,3)=gz(k)
 30   continue
c
c if iopt=2, place a second staggered grid of reference particles
c
      if (iopt.eq.2) then
c
      ix1=ix1+intx/2
      ix2=ix2+intx/2
      iy1=iy1+inty/2
      iy2=iy2+inty/2
      iz1=iz1+intz/2
      iz2=iz2+intz/2
      write (6,*) 'tetra (iopt=2): ix1,ix2=',ix1,ix2
      write (6,*) 'tetra (iopt=2): iy1,iy2=',iy1,iy2
      write (6,*) 'tetra (iopt=2): iz1,iz2=',iz1,iz2
c
      do 35 k=iz1,iz2,intz
      do 35 j=iy1,iy2,inty
      do 35 i=ix1,ix2,intx
      index=index+1
      posn(index,1)=gx(i)
      posn(index,2)=gy(j)
      posn(index,3)=gz(k)
 35   continue
c
      end if
c
      else if (iopt.eq.3) then
c
c iopt=3: distribute the reference particles (to be used as
c base vertices of the tetrahedra) randomly
c
      lx=xyzl(1)
      ly=xyzl(2)
      lz=xyzl(3)
c
      npu=nop
      npy=nop/ny
      npe=nop-npy*ny
      niny=0
	allocate (postmp(npy,3))
c
      call ranseq (kipran,0)
c
      do 10 y=1,ny
c
      ni=npy
      if (y.le.npe) ni=npy+1
      ifst=niny+1
c
      call ranu2 (postmp(1,1),1,ni,1,1,ni)
      call ranu2 (postmp(1,2),1,ni,1,1,ni)
      call ranu2 (postmp(1,3),1,ni,1,1,ni)
c
      do 20 ii=1,ni
      i=niny+ii
c     write (101,901) i,posn(i,1),posn(i,2),posn(i,3)
c901  format ('i,pp=',i6,1p,3e12.4)
#ifdef FEB2411
      posn(i,1)=postmp(ii,1)*nx+gx(1)
      posn(i,2)=postmp(ii,2)+gy(y)
      posn(i,3)=postmp(ii,3)*nz+gz(1)
#else
      posn(i,1)=postmp(ii,1)*lx+gx(1)
      posn(i,2)=postmp(ii,2)/ps(2)+gy(y)
      posn(i,3)=postmp(ii,3)*lz+gz(1)
#endif
 20   continue
c
      niny=niny+ni
c
 10   continue
c
      nop=(1+3*ngp)*nop
      write (6,*) 'tetra: iopt=3, nop=',nop
c
      end if
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
#ifdef FEB2411
      rx=dxpr(igp)*nx/beta1
      ry=dypr(igp)*ny/beta2
      rz=dzpr(igp)*nz/beta3
#else
      rx=dxpr(igp)*nx
      ry=dypr(igp)*ny
      rz=dzpr(igp)*nz
#endif
c
      if (shear.eq.0) then
c
      if (iopt.le.2) then
      npu3=iopt*npu**3
      ipfst=(igp-1)*3*npu3+npu3+1
      iplst=ipfst+3*npu3-1
      else if (iopt.eq.3) then
      ipfst=(igp-1)*3*npu+npu+1
      iplst=ipfst+3*npu-1
      end if
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
	if (ipref.eq.1) then
	write (8100,*) 'ipref,pp=',ipref,(posn(ipref,j),j=1,3)
	write (8100,*) 'index,pp=',index,(posn(index,j),j=1,3)
	write (8100,*) 'index+1,pp=',index+1,(posn(index+1,j),j=1,3)
	write (8100,*) 'index+2,pp=',index+2,(posn(index+2,j),j=1,3)
	end if
      end do
c
      else if (shear.gt.0.) then
c revision 2/16/97 for shear flow
c
      npu3=iopt*npu**3
      ipfst=1
      iplst=npu3
      write (6,"('igp,ipfst,iplst=',i3,2i8,'  rx,ry,rz=',1p,3e12.4)")
     1             igp,ipfst,iplst,rx,ry,rz
c
      do ipref=ipfst,iplst
c
      ii=mod(ipref,8)
      if (ii.eq.0) ii=8
c
      ipx=iplst+(ipref-ipfst)*3+1+(igp-1)*3*npu3
      ipy=ipx+1
      ipz=ipx+2
c
      posn(ipx,1)=posn(ipref,1)+isignx(ii)*rx
      posn(ipx,2)=posn(ipref,2)
      posn(ipx,3)=posn(ipref,3)
      posn(ipy,1)=posn(ipref,1)
      posn(ipy,2)=posn(ipref,2)+isigny(ii)*ry
      posn(ipy,3)=posn(ipref,3)
      posn(ipz,1)=posn(ipref,1)
      posn(ipz,2)=posn(ipref,2)
      posn(ipz,3)=posn(ipref,3)+isignz(ii)*rz
c
      end do
c
      end if
c
c write pairing information
c
c Comment this out: too time-consuming and not strictly necessary
#ifdef TEST
      if (shear.eq.0.) then
c
c 3 pairs for each cluster in isotropic turbulence
c
      if (ngp.eq.1) then
      call fopen1 (70,'tetra','formatted  ')
      else
      call fopen1 (70,'tetra'//numc(igp),'formatted  ')
      end if
c
      if (iopt.le.2) then
      write (70,601) 3*npu**3
      else
      write (70,601) 3*npu
      npu3=npu
      end if
c
      ipref=0
      ipair=0
      do index=ipfst,iplst-2,3
      ipref=ipref+1
      ii=index-(igp-1)*npu3*3
      ipair=ipair+1
      write (70,602) ipair,ipref,ii,index,(posn(ipref,j),j=1,3),
     1                              (posn(index,j),j=1,3)
      ipair=ipair+1
      write (70,602) ipair,ipref,ii+1,index+1,(posn(ipref,j),j=1,3),
     1                                (posn(index+1,j),j=1,3)
      ipair=ipair+1
      write (70,602) ipair,ipref,ii+2,index+2,(posn(ipref,j),j=1,3),
     1                                (posn(index+2,j),j=1,3)
      end do
c
      call fclose1 (70)
c
c different for shear flow
c
      else
c
      if (ngp.eq.1) then
      call fopen1 (71,'xtetra','formatted  ')
      call fopen1 (72,'ytetra','formatted  ')
      call fopen1 (73,'ztetra','formatted  ')
      else
      fnx='xtetra'//numc(igp)
      fny='ytetra'//numc(igp)
      fnz='ztetra'//numc(igp)
      call fopen1 (71,fnx,'formatted  ')
      call fopen1 (72,fny,'formatted  ')
      call fopen1 (73,fnz,'formatted  ')
      end if
c
      write (71,601) npu3
      write (72,601) npu3
      write (73,601) npu3
c601  format ('no. of pairs in this subset=',i6)
 601  format ('no. of pairs=',i6)
c
c revised for shear flow, 2/16/97
c
      npu3=iopt*npu**3
      ipfst=1
      iplst=npu3
      write (6,"('igp,ipfst,iplst=',i3,2i8,'  rx,ry,rz=',1p,3e12.4)")
     1             igp,ipfst,iplst,rx,ry,rz
c
      sum1=0.
      sum2=0.
      sum3=0.
      sum4=0.
c
      ii=(igp-1)*3*npu3
c
      do ipref=ipfst,iplst
c
      ipx=iplst+(ipref-ipfst)*3+1+ii
      ipy=ipx+1
      ipz=ipx+2
c
      ipair=ipref-ipfst+1
c
      write (71,602) ipair,ipref,ipx-ii,(posn(ipref,j),j=1,3),
     1                              (posn(ipx,j),j=1,3)
      write (72,602) ipair,ipref,ipy-ii,(posn(ipref,j),j=1,3),
     1                              (posn(ipy,j),j=1,3)
      write (73,602) ipair,ipref,ipz-ii,(posn(ipref,j),j=1,3),
     1                              (posn(ipz,j),j=1,3)
c
      end do
c
c602  format (3i6,2(2x,3f6.2))
c changed by DD. Numbers were too large for large problems.
c602  format (3i6,2(2x,3e12.5))
 602  format (4i6,2(2x,3e12.5))
c
      call fclose1 (71)
      call fclose1 (72)
      call fclose1 (73)
c
      end if
c
#endif
c
 40   continue
c
 210  format ('particle clusters in the form of tetrahedra'/
     1        'no. of clusters(i4)=',i4)
 201  format ('iref,ip1,ip2(3i6)=',3i6)
c
	write (6,*) ' exit tetra'
#endif
      return
      end
