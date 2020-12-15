        subroutine covcof
!
	use comsp
!
! write covariance and correlation coefficient matrices
!
! routine to output the eulerian statistics computed by sub. spectr
!
	real(b8), allocatable :: crcof(:)
	real(b8) denom
!
	allocate (crcof(ncp),stat=ierr)
!
      if (shear.eq.0.) then
        write(30,600) istep,time,dt
	else
        write(30,6001) istep,time,dt,shear*time
	end if
!
!  covariance matrix
!
        ijf=0
        do 5 i=1,3+nc
        ijs=ijf+1
        ijf=ijs+(3+nc)-i
        write(30,651) (corr(ij),ij=ijs,ijf)
 5    continue
!
! correlation coefficients
!
        ijf=0
        do 6 i=1,3+nc
        ijs=ijf+1
        ijf=ijs+(3+nc)-i
      do 7 ij=ijs,ijf
      denom=corr(ijs)*corr(kt(i+ij-ijs))
	if (denom.eq.0.) then
	crcof(ij)=0.
	else
	crcof(ij)=corr(ij)/sqrt(denom)
	end if
c     crcof(ij)=corr(ij)/sqrt(corr(ijs)*corr(kt(i+ij-ijs)))
 7    continue
        write(30,651) (crcof(ij),ij=ijs,ijf)
 6    continue
!
	deallocate (crcof)
!
        return
!
600     format(/'output for step ',i6,'   t,dt=',1p,2e13.5)
6001    format(/'output for step ',i6,'   t,dt=',1p,2e13.5,
     1          '   st=',1p,e13.5)
651     format(1p,10e13.6)
!
        end
