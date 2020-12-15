      subroutine splx (f,bsx,ux,vx)
c
#ifdef LAG
	use mpilag
      implicit none

      integer nxdim
      real f(nx,zisz)

c ux,vx now allocated outside of this routine

	real(p8) bsx(nbx,zisz)
	real(p8) ux(0:nx),vx(0:nx)
c
      integer x,z
      integer xm1,xp1
      integer nxp1,nxp2,nxp3
c
	integer zp
c
c
c form nz x-splines
c
      ux(0)=0.
      vx(nx)=0.
c
      nxp1=nx+1
      nxp2=nx+2
      nxp3=nx+3
c
      do 100 zp=1,zisz
c
      do 10 x=1,nx
         xm1=x-1
         ux(x)=(f(x,zp)-ux(xm1))*p(x,1)
 10   continue
c
      do 20 x=nx-1,1,-1
         xp1=x+1
         vx(x)=q(x,1)*vx(xp1)+ux(x)
 20   continue
c
      bsx(nxp1,zp)=denom(1)*(f(nx,zp)-vx(1)-vx(nx-1))
c
      do 30 x=nx-1,1,-1
         bsx(x+1,zp)=t(x,1)*bsx(nxp1,zp)+vx(x)
 30   continue
c
      bsx(1,zp)=bsx(nxp1,zp)
      bsx(nxp2,zp)=bsx(2,zp)
      bsx(nxp3,zp)=bsx(3,zp)
c     
 100  continue

#endif
      return
      end
