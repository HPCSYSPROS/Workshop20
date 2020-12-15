	subroutine hsptr (uy,m)
c
c calculation of "horizontal spectrum" (in 2D kx-ky plane)
c
!!I dont believe this needs to call the bouss_module
      use comsp
c
      implicit none
	include 'intvars'

	complex(b8) :: uy(ny,zjsz,xisz,3)
	integer :: i,j,k,l,ik,m,lu
	real(b8) :: rk2,term,term1,term2,sum,beta1r,tfact_x
	real(b8), allocatable :: hekz(:,:),hek(:,:)
c
	integer :: iz,iy,nzp2,nyp2,xst
	complex(b8) :: czero

	czero=cmplx(0.,0.)
	nzp2=nz+2
	nyp2=ny+2
	
	allocate (hekz(nx,nz))
	allocate (hek(nx,nz))

c	tfact(:)=2.
c	tfact(1)=1.

	hekz(:,:)=0.

	kx2(:)=b11(m)*kx(:)**2
	ky2(:)=b22(m)*ky(:)**2
c	kz2(:)=b33(m)*kz(:)**2
c
c	ky2(nyhp)=0.
c	kz2(nzhp)=0.

	beta1r = 1.0/beta1

      do 10 xp=1,xisz
	x=xp+xist-1

	if(x .eq. 1) then
	   tfact_x=1
	else
	   tfact_x=2
	endif
	
      do 11 zp=1,zjsz
	z=zp+zjst-1
	if (z.eq.nzhp) go to 11

      do 12 y=1,ny
	
	if (.not. mask(y,zp,xp)) go to 12

c	call masks(y)

 	rk2=kx2(x)+ky2(y)
	ik=sqrt(rk2)*beta1r+1.5
c
	term1=uy(y,zp,xp,1)*conjg(uy(y,zp,xp,1))*b11(m)
	term2=uy(y,zp,xp,2)*conjg(uy(y,zp,xp,2))*b22(m)
	term=tfact_x*(term1+term2)
 	hekz(ik,z)=hekz(ik,z)+.5*term
c
 12	continue	

 11	continue
c
 10	continue
c

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c collect contributions from all processors

 
 	call MPI_ALLREDUCE (hekz,hek,nx*nz,mpireal,
     1   MPI_SUM,MPI_COMM_WORLD,ierr)
c
	if (taskid.eq.0) then
c
#ifdef TEST
	lu=61
	write (lu,600) istep, time
	do z=1,nz/2
	write (lu,601) z,kz(z)*beta3
	if (z.eq.1) then
	write (lu,602) (hek(k,z),k=1,nxh)
	else 
	iz=nz+2-z
	write (lu,602) (hek(k,z)+hek(k,iz),k=1,nxh)
	end if
	end do
#endif
 600	format ('horizontal spectrum, istep,time=',i6,1p,e12.4)
 601	format ('horizontal spectrum, z, kz=',i5,f8.1)
 602	format (1p,6e12.4)
c
	lu=62
	write (lu,600) istep, time
	do k=1,nxh
	sum=0.
	do z=1,nz
	sum=sum+hek(k,z)
	end do
	hek(k,1)=sum
	end do
	write (lu,602) (hek(k,1),k=1,nxh)
c
	end if
c
        deallocate (hekz,hek)
	if (taskid.eq.0) write (6,*) ' exit hsptr, istep=',istep

	return
	end
