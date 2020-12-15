c 
c test capability for freezing energy in forced wavenumber shells
c
!NPM finishing stride-1 convert 3/9/09
	subroutine ek_shell (uny,m)
c
c converted to cylindrical truncation scheme, by P.K.Y and D.D., 11/29/09
       
#ifdef FEK_FORC
	use comsp
       
	implicit none
	include 'intvars'

	integer a
	complex(b8) :: uny(ny,zjsz*xisz,3+nc)       
c
	real(b8) :: rk2,term
	real, allocatable :: eky_nf(:),eknf(:)
	real term1,term2,term3
	integer :: io,i,j,k,l,ij,ncount,ik,m
	integer icall
	character*20 caux
	data icall/0/
	save icall

	real(b8) deltak


	if (taskid.eq.0.and.icall.eq.0) then
	caux='fek_modes'
	call fopen1 (91,caux,'formatted  ')
        write (6,*) 'enter ek_shell,kf_shell=',kf_shell
        write (91,*) 'enter ek_shell,kf_shell=',kf_shell
	end if
	
	allocate (eky_nf(kfor),eknf(kfor))
	icall=icall+1

c note "tfact(:)" array is set in comsp_set.f
c
	eky_nf(:)=0.


	kx2(:)=kx(:)**2*b11(2)
	ky2(:)=ky(:)**2*b22(2)
	kz2(:)=kz(:)**2*b33(2)

#ifdef SHELL_DK
	deltak=amin1(beta1,beta2,beta3)
#else
	deltak=1.
#endif
c
      xp = mystart_x
      z = mystart_z

	do 10 a=1,num_al
	  x=xp+xist-1
	  if (z.eq.nzhp) go to 10

	  do 12 y=1,ny
	     if (.not. mask(y,a)) go to 12
		 
	     rk2=kx2(x)+ky2(y)+kz2(z)
	     ik=sqrt(rk2)/deltak+1.5
c     if (ik.gt.kf_shell) go to 12
	     if (ik.gt.kf_shell+1) go to 12
	     if (icall.eq.1) then 
	write (910+taskid,"('ekshell:,x,y,z,ik,k=',4i4,1p,e12.4)") 
     1          x,y,z,ik,sqrt(rk2)
		end if
	     term1=real(uny(y,a,1)*conjg(uny(y,a,1)))*b11(2)
	     term2=real(uny(y,a,2)*conjg(uny(y,a,2)))*b22(2)
	     term3=real(uny(y,a,3)*conjg(uny(y,a,3)))*b33(2)
	     term=(term1+term2+term3)/2.*tfact(x)
	     eky_nf(ik)=eky_nf(ik)+term
 12	  continue

	  call next_xz(xp,z)
 10   continue

	call MPI_ALLREDUCE (eky_nf,eknf,kfor,mpireal,MPI_SUM,MPI_COMM_WORLD,ierr)


	if (icall.eq.1) then
	   ek_nf1(:)=eknf(:)
	   ek_nf2(:)=eknf(:)
	else
	   ek_nf1(:)=ek_nf2(:)
	   ek_nf2(:)=eknf(:)
	end if
	
	deallocate (eky_nf,eknf)
	
	if (icall.eq.1.and.taskid.eq.0) close (91)
	
#endif
	return
	end subroutine ek_shell
