        subroutine out1d (lu1,lu2,io)
!
#ifdef OUT1D
! routine to output 1-d eulerian statistics computed by sub. spectr
! (diagonal components only)
!
! lu1 for 1-dimensional spectra
! lu2 for 1-dimensional spatial correlation
!
	use comsp
!
	implicit none
	include 'intvars'
!
	integer i,j,l,lu1,lu2,io,ij
	real unit
	real, allocatable :: sp1x(:,:),sp1y(:,:),sp1z(:,:)
	real factor1
	complex, allocatable :: work(:)
        integer ndim
!
        
        ndim=max(nx,ny,nz)
	allocate (sp1x(nxhp,3+nc),sp1y(nyhp,3+nc),sp1z(nzhp,3+nc))
	allocate (qijx(nxhp,ncp),qijy(nyhp,ncp),qijz(nzhp,ncp))
	allocate (work(ndim*2))
!
	do 5 i=1,nc+3
	do x=1,nxhp
	sp1x(x,i)=real(xrij(x,kt(i),io))
	end do
	do y=1,nyhp
	sp1y(y,i)=real(yrij(y,kt(i),io))
	end do
	do z=1,nzhp
	sp1z(z,i)=real(zrij(z,kt(i),io))
	end do
 5 	continue
!
!
      write (lu1,210) istep,time,3+nc
      write (lu2,210) istep,time,3+nc
!
! only the necessary numbers are written on file, bearing in mind that:
!
! 1) diagonal terms of the spectrum tensors are real,
!    and nonzero only for wavenumbers 0 to n/2-1
!
! 2) spatial correlations at zero separation are always unity
!
! 3) since the 1-d correlations are even functions, they have
!    fourier cosine transforms in the 1-d spectrum arrays
!
!
      write(lu1,211)
      do 10 i=1,3+nc
      write (lu1,201) i, (sp1x(l,i),l=1,10)
      write (lu1,202)    (sp1x(l,i),l=11,nx/2)
 10   continue
!
      write(lu1,212)
      do 20 i=1,3+nc
      write (lu1,201) i, (sp1y(l,i),l=1,10)
      write (lu1,202)    (sp1y(l,i),l=11,ny/2)
 20   continue
!
      write(lu1,213)
      do 30 i=1,3+nc
      write (lu1,201) i, (sp1z(l,i),l=1,10)
      write (lu1,202)    (sp1z(l,i),l=11,nz/2)
 30   continue
!
! since xrij is 1-d spectrum,
! copy xrij into qijx, then take inverse fft to get two-point
! spatial correlations (which are real) in physical space.
! similarly for the other 2 dirns.
!
      do 50 ij=1,ncp
!
      do 51 x=1,nxhp
 51   qijx(x,ij)=xrij(x,ij,io)
      do 52 y=1,nyhp
 52   qijy(y,ij)=yrij(y,ij,io)
      do 53 l=1,nzhp
 53   qijz(l,ij)=zrij(l,ij,io)
!
      call rfft99 (qijx(1,ij),work,1,1,nx,1,1)
      call rfft99 (qijy(1,ij),work,1,1,ny,1,1)
      call rfft99 (qijz(1,ij),work,1,1,nz,1,1)
!
 50   continue
!
	call corr1d (qijx,qijy,qijz,lu2)
!
	deallocate (sp1x,sp1y,sp1z)
	deallocate (qijx,qijy,qijz)
	deallocate (work)
!
        return
!
210     format(/,1x,' istep,time,3+nc=',i5,1p,e14.5,2x,i4)
211     format(1x,'one dimensional x-spectra')
212     format(1x,'one dimensional y-spectra')
213     format(1x,'one dimensional z-spectra')
221     format(1x,'one dimensional x-correlation')
222     format(1x,'one dimensional y-correlation')
223     format(1x,'one dimensional z-correlation')
!
 201  format (1x,i2,1p,10e13.5)
 202  format (3x,1p,10e13.5)
!
#endif
        end
