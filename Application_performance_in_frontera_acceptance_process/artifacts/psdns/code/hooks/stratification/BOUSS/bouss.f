      subroutine bouss (uy,uny,uby,icall)


#ifdef BOUSS
c
 
	use stratification_module
c
      use comp
      implicit none
	include 'intvars'
 
      complex(b8) :: uy(ny,xisz*zjsz,nu)
      complex(b8) :: uny(ny,xisz*zjsz,3+nc)
      complex(b8) :: uby(ny,xisz*zjsz,3:4)
c
	integer icall
	real const
	real bvf
c
	integer a
c
	real factor3,factor4
c

c	if (taskid.eq.0) write (6,*) 'enter bouss, icall=',icall
c
	pi=atan(1.)*4.
	if (froude.gt.1.e6) return
	const=(2.*pi/froude)**2
c
      xp = mystart_x
      z = mystart_z
c

	go to (1,2,3), icall
c
 1	continue
c
	factor3=dt/2.*const/beta3
	factor4=dt/2.*beta3
c
	do 10 a=1,num_al
	x=xp+xist-1
	if (z.eq.nzhp) go to 10
	do 12 y=1,ny
	if (.not.mask(y,a)) go to 12
 	uy(y,a,3)=uy(y,a,3)-uny(y,a,4)*factor3
   	uy(y,a,4)=uy(y,a,4)+uny(y,a,3)*factor4
 12	continue
	call next_xz (xp,z)
 10	continue
c
	return
c
 2	continue
c
	do 20 a=1,num_al
	x=xp+xist-1
	if (z.eq.nzhp) go to 20
	do 22 y=1,ny
	if (.not.mask(y,a)) go to 22
 	uby(y,a,3)=uy(y,a,3)
	uby(y,a,4)=uy(y,a,4)
 22	continue
	call next_xz (xp,z)
 20	continue
c
	return
c
 3	continue
c
c on second thought, perhaps it is more accurate to use uby than uny below
c for update of uy(y,a,3)
c	uy(y,a,3)=uy(y,a,3)-uny(y,a,4)*dt/2.*const/beta3
c
	factor3=dt/2.*const/beta3
	factor4=dt/2.*beta3
c
	do 30 a=1,num_al
	x=xp+xist-1
	if (z.eq.nzhp) go to 30
	do 32 y=1,ny
	if (.not.mask(y,a)) go to 32
	uy(y,a,3)=uy(y,a,3)-uby(y,a,4)*factor3
  	uy(y,a,4)=uy(y,a,4)+uby(y,a,3)*factor4
 32	continue
	call next_xz (xp,z)
 30	continue

      return
c
#endif
      end
