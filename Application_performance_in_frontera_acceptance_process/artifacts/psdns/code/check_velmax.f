	subroutine check_velmax (ux,char)
!
	use comp
	implicit none
	include 'intvars'
	character(*) char
!
	integer icall
c
 	real(b8)    :: ux(nx,zisz,yjsz,nu)

!
      real(b8) s1,s2,s3,xnorm,ynorm,znorm,vel,k3z,bk3z
c
	data icall/0/
c

	icall=icall+1

      xnorm=b11(1)*nx
      ynorm=b22(1)*ny
      znorm=b33(1)*nz
!
	velmax=0.
!
      do 120 yp=1,yjsz
      do 120 zp=1,zisz
!
      do 125 x=1,nx
!
         vel=xnorm*abs(ux(x,zp,yp,1))+ynorm*abs(ux(x,zp,yp,2))
     1      +znorm*abs(ux(x,zp,yp,3))
         velmax=max(vel,velmax)
c
 125   continue
!
 120   continue
!
	if (icall.eq.1) then
	write (4000+taskid,*) 'check_velmax, xnorm=',xnorm,ynorm,znorm
	end if
c
	write (6,900) char,taskid,velmax
	write (4000+taskid,900) char,taskid,velmax
 900	format ('check_velmax:char,taskid,velmax=',a20,i5,1p,e12.4)
c
	return
	end
