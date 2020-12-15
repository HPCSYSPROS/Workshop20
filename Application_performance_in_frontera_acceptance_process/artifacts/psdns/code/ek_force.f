	subroutine ek_force  (uy,uny,m)
c
#ifdef FEK_FORC
c 
c test capability for freezing energy in forced wavenumber shells
c
c converted to cylindrical truncation scheme, by P.K.Y. and D.D., 11/29/09
c
      use comsp
c
      implicit none
	include 'intvars'

c
	real, allocatable :: ratio(:)
c
	integer m,ik
	real rk2
c
	integer a
	complex(b8) :: uny(ny,zjsz*xisz,3+nc)       
	complex(b8) :: uy(ny,zjsz*xisz,3+nc)

        real(b8) deltak

c
	allocate (ratio(kfor))
c
	call ek_shell (uny,m)
c
	ratio(1)=0.
	do ik=2,kfor
c	ratio(ik)=sqrt(ek_nf1(ik)/ek_nf2(ik))
	ratio(ik)=sqrt(ekinf(ik)/ek_nf2(ik))
	end do

	kx2(:)=kx(:)**2*b11(2)
	ky2(:)=ky(:)**2*b22(2)
	kz2(:)=kz(:)**2*b33(2)
c
#ifdef SHELL_DK
        deltak=amin1(beta1,beta2,beta3)
#else
        deltak=1.
#endif

      xp = mystart_x
      z = mystart_z

	do 10 a=1,num_al
	   x=xp+xist-1
	   if (z.eq.nzhp) go to 10
	   do 12 y=1,ny
	     if (.not.mask(y,a)) go to 12
 	     rk2=kx2(x)+ky2(y)+kz2(z)
 	     ik=sqrt(rk2)/deltak+1.5
	     if (ik.gt.kf_shell) go to 12
    	     uny(y,a,1)=uny(y,a,1)*ratio(ik)
    	     uny(y,a,2)=uny(y,a,2)*ratio(ik)
    	     uny(y,a,3)=uny(y,a,3)*ratio(ik)
    	     uy(y,a,1)=uy(y,a,1)*ratio(ik)
    	     uy(y,a,2)=uy(y,a,2)*ratio(ik)
    	     uy(y,a,3)=uy(y,a,3)*ratio(ik)
 12	   continue
	   call next_xz(xp,z)
 10	continue   
c
c
	call ek_shell (uny,m)
c
	deallocate (ratio)
c
c
#endif
	return
	end
