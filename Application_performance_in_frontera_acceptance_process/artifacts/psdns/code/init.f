	subroutine init()

	use com
	implicit none
	include 'intvars'
 
	integer :: i

! define some used constants
	imagi=cmplx(0.,1.)
        pi=4.*atan(1.)

! initialize wavenumbers

	do x=1,nxhpad
  	  kx(x)=x-1
 	end do

	do y=1,nypad
	  ky(y)=y-1
 	  if(y.gt.nyhppad) ky(y)=ky(y)-ny
 	end do

      do z=1,nzpad
        kz(z)=z-1
        if(z.gt.nzhppad) kz(z)=kz(z)-nz
	end do

	kmax=sqrt(2.)/3.*max(nxpad*beta1,nypad*beta2,nzpad*beta3)

!  define mesh:  b11=b(1,1)**2 etc., b12=b(1,2)/b(2,2).
 
	do 4 i=1,2
	b11(i)=beta1**2
	b22(i)=beta2**2
	b33(i)=beta3**2
 4    b12(i)=bet12
 
!  form scaled scalar gradients
 
#ifndef NOSCALAR
        do 5 i=1,nc
        cb(1,i)=beta1*grad(1,i)
        cb(2,i)=beta2*grad(2,i)
5       cb(3,i)=beta3*grad(3,i)
#endif

	return
	end
