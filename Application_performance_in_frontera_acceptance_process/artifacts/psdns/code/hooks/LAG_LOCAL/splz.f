      subroutine splz (bs)
c
#ifdef LAG
	use mpilag
      implicit none

      real bs(bxistart:bxiend,nbz)
      real uz(0:nz),vz(0:nz),wz(0:nbz)
c
      integer x,z
      integer xm1,xp1,nzp1
c
c
c form nbx z-splines
c
      uz(0)=0.
      vz(nz)=0.
      nzp1=nz+1
c
      do 200 x=bxistart,bxiend
c
      do 50 z=1,nz
         uz(z)=(bs(x,z)*216.-uz(z-1))*p(z,3)
 50   continue
c
      do 60 z=nz-1,1,-1
         vz(z)=q(z,3)*vz(z+1)+uz(z)
 60   continue
c
      wz(nzp1)=denom(3)*(bs(x,nz)*216.-vz(1)-vz(nz-1))
c
      do 70 z=nz-1,1,-1
         wz(z+1)=t(z,3)*wz(nzp1)+vz(z)
 70   continue
c
      do 80 z=2,nz
         bs(x,z)=wz(z)
 80   continue
      bs(x,nzp1)=wz(nzp1)
c
c     bs(x,nz+1)=wz(1)
c     bs(x,nz+2)=wz(2)
c     bs(x,nz+3)=wz(3)
c
 200  continue
c
      do x=bxistart,bxiend
         bs(x,1)=bs(x,nz+1)
         bs(x,nz+2)=bs(x,2)
         bs(x,nz+3)=bs(x,3)
      end do
c
c
#endif
      return
      end
