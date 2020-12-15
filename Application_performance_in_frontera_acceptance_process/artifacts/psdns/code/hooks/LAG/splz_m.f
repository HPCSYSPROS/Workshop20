      subroutine splz_m (bs,nv,uz,vz,wz)
c
#ifdef LAG
#ifdef SPXYZM
	use mpilag
      implicit none

c uz,vz,wz now allocated outside of this routine
c
      real(p8) bs(nv,bxisize,nbz)
      real(p8) uz(0:nz),vz(nz),wz(nz+1)
c
      integer x,z
      integer xm1,xp1,nzp1
c
	integer iv,nv
c
c
c form nbx z-splines
c
      uz(0)=0.
      vz(nz)=0.
      nzp1=nz+1
c
	do 200 iv=1,nv
      do 200 x=1,bxisize
c
      do 50 z=1,nz
         uz(z)=(bs(iv,x,z)*216.-uz(z-1))*p(z,3)
 50   continue
c
      do 60 z=nz-1,1,-1
         vz(z)=q(z,3)*vz(z+1)+uz(z)
 60   continue
c
      wz(nzp1)=denom(3)*(bs(iv,x,nz)*216.-vz(1)-vz(nz-1))
c
      do 70 z=nz-1,1,-1
         wz(z+1)=t(z,3)*wz(nzp1)+vz(z)
 70   continue
c
      do 80 z=2,nz
         bs(iv,x,z)=wz(z)
 80   continue
      bs(iv,x,nzp1)=wz(nzp1)
c
c
 200  continue
c
	do iv=1,nv
      do x=1,bxisize
         bs(iv,x,1)=bs(iv,x,nz+1)
         bs(iv,x,nz+2)=bs(iv,x,2)
         bs(iv,x,nz+3)=bs(iv,x,3)
      end do
      end do
c
#endif
#endif
      return
      end
