	subroutine param
c specified parameters
	use com

c set grid
 	nx=8  ; ny=8  ; nz=8
c	nx=64 ; ny=64 ; nz=64

	nu=6

c number of scalars	
	nc=0


c initial mesh metric
	beta1=1. ; beta2=1. ; beta3=1. 
	bet12=0.

c grid shift option (see sub. shifts)
c option 1 is original rogallo's scheme
c option 3 is used by yeung & pope (jfm, 1989)
c option 4 is zero shift
c
	kshift=4
	

	return
	end
