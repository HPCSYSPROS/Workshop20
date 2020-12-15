c 
c test capability for freezing energy in forced wavenumber shells
c
!NPM finishing stride-1 conversion 3/9/09
	subroutine ek_init  (uny)

c converted to cylindrical truncation scheme by P.K.Y. and D.D. 11/29/09

#ifdef FEK_FORC
	use comsp	
	implicit none
	include 'intvars'
	real, allocatable :: kden(:)
	logical iex	
	character*20 :: caux	
	real, allocatable :: ratio(:)
	integer m,ik,iflag,nnn
	real rk2

	integer a
	complex(b8) :: uny(ny,zjsz*xisz,3+nc)       
c
	real(b8) deltak
c
	allocate (kden(nxh))
	kden=1.

	if (taskid.eq.0) then

	ekinf(:)=ekinf_input(:)
	deallocate (ekinf_input)


	write (caux,fmt="('den',i0)") min0(nx,4096)
      call fopen1 (4,caux,'formatted  ')
c       
	if (nx.le.128) then
	   do ik=2,nxh
	      read (4,102) kden(ik)
	   end do
	else
	   do ik=2,min0(nxh,2048)
	      read (4,1021,err=888,end=888) kden(ik)
	   end do
	end if
	close(4)

 	go to 889
 888	call abrt('Error when reading denXXX file')	
 889	continue 
c       
	end if
c       
 	call MPI_BCAST (ekinf,nxh,mpireal,0,MPI_COMM_WORLD,mpierr)
	call MPI_BCAST (kden,nxh,mpireal,0,MPI_COMM_WORLD,mpierr)
c       
	if (taskid.eq.0) write (6,*) 'ek_init: ekinf,kden=',ekinf(2),kden(2)

c should we do this (or not) --- probably yes, assuming ekinp comes
c from eultav.x which divided the spectrum by kden(ik)
c
        do ik=2,kmax+1
        ekinf(ik)=ekinf(ik)*kden(ik)
        end do
c       
	allocate (ratio(kfor))
c       
	call ek_shell (uny,m)
c       
c       
	ratio(1)=0.
	do ik=2,kfor
	   ratio(ik)=sqrt(ekinf(ik)/ek_nf2(ik))
	if (taskid.eq.0) then
 	   write (6,"('ek_init: ik,ekinf,ek_nf2,ratio',i3,1p,3e12.4)")
     1          ik,ekinf(ik),ek_nf2(ik),ratio(ik)
	end if
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
c
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

 12	      continue
c
	call next_xz(xp,z)
c       
 10     continue
c       
c       
	deallocate (kden,ratio)
c       
 101	format (10e13.5)
 102	format (28x,f6.4)
 1021	format (29x,f6.4)
 501	format (1x)

c       
#endif
	return
	end subroutine ek_init
