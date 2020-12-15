      subroutine pcheck (unxr)
c
      use com
	include 'intvars'
c
c set initial velocity field in physical space
c
c sinusoidal velocity field (products of sines and cosines)
c
c  u = au * sin(x)*cos(y)*cos(z)
c  v = av * cos(x)*sin(y)*cos(z)
c  w = aw * cos(x)*cos(y)*sin(z)
c
c  continuity requires au + av + aw = 0.
c
      real(b8) :: unxr(nx,zisz,yjsz,3)
c
      real(b8), allocatable :: sinx(:),siny(:),cosx(:),cosy(:)
c
      allocate (sinx(nx),stat=ierr)
      allocate (cosx(nx),stat=ierr)
      allocate (siny(ny),stat=ierr)
      allocate (cosy(ny),stat=ierr)
c
      data au,av/1.,1./
c
      yg=taskid+1
c     write (6,*) 'enter pcheck, taskid=',taskid, nxhp2,zisz,yjsz,
c    1             Ntot
c
      aw=-(au+av)
c
      pi=atan(1.)*4.
      dx=2.*pi/nx
      dy=2.*pi/ny
      dz=2.*pi/nz
c
      do 1 x=1,nx
      sinx(x)=sin((x-1)*dx)
      cosx(x)=cos((x-1)*dx)
 1    continue
      do 2 y=1,ny
      siny(y)=sin((y-1)*dy)
      cosy(y)=cos((y-1)*dy)
 2    continue
c
      umax=0.
      difmax=0.
      difm1=0.
      difm2=0.
      difm3=0.
c
      do 10 yp=1,yjsz
c
      y=yjst+yp-1
c     siny=sin((y-1)*dy)
c     cosy=cos((y-1)*dy)
c
      do 20 zp=1,zisz
      z=zist+zp-1
      sinz=sin((z-1)*dz)
      cosz=cos((z-1)*dz)
c
      termu=cosz*cosy(y)
c     termv=sinz*cosy(y)
c     termw=cosz*siny(y)
      termv=siny(y)*cosz
      termw=cosy(y)*sinz
c
      do x=1,nx
      u1=au*sinx(x)*termu
      u2=av*cosx(x)*termv
      u3=aw*cosx(x)*termw
      if (abs(unxr(x,zp,yp,3)).gt.umax) then
      umax=abs(unxr(x,zp,yp,3)) 
      ixm=x
      iym=y
      izm=z
      end if
      dif1=abs(unxr(x,zp,yp,1)-u1)
      dif2=abs(unxr(x,zp,yp,2)-u2)
      dif3=abs(unxr(x,zp,yp,3)-u3)
      difm1=amax1(dif1,difm1)
      difm2=amax1(dif2,difm2)
      difm3=amax1(dif3,difm3)
c     difm=amax1(dif1,dif2,dif3)
      difm=dif1
      if (difm.gt.difmax) then
      difmax=difm
      exact=u1
      actual=unxr(x,zp,yp,1)
      ixmax=x
      iymax=y
      izmax=z
      end if
      end do
c
 20   continue
c
 10   continue
c
      write (6,*) 'pcheck: max.abs.value=',umax,ixm,iym,izm,
     1            'yjsz,zisz=',yjsz,zisz
      write (6,*) 'pcheck: max-diff=',difmax,ixmax,iymax,izmax,
     1             exact,actual
      write (6,*) 'pcheck: difm1,2,3=',difm1,difm2,difm3
c
      write (6,*) ' exit pcheck, taskid=',taskid
      return
      end
