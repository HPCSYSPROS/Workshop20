        subroutine corr1d (xcor,ycor,zcor,lu2)
c
#ifdef OUT1D
c
c routine to output 1-d eulerian statistics computed by sub. spectr
c (diagonal components only)
c
c lu1 for 1-dimensional spectra
c lu2 for 1-dimensional spatial correlation
c
	use comsp
c
	implicit none
	include 'intvars'
c
	integer i,lu2,ij
	real factor1
c
	real xcor(nx+2,ncp),ycor(ny+2,ncp),zcor(nz+2,ncp)

c
c normalisation to get correlation coefficients
c
      do 55 ij=1,ncp
c
      factor1=xcor(1,ij)
      if (factor1.ne.0.) then
      do x=nx+2,1,-1
      xcor(x,ij)=xcor(x,ij)/factor1
      end do
      end if
c
      factor1=ycor(1,ij)
      if (factor1.ne.0.) then
      do y=ny+2,1,-1
      ycor(y,ij)=ycor(y,ij)/factor1
      end do
      end if
c
      factor1=zcor(1,ij)
      if (factor1.ne.0.) then
      do z=nz+2,1,-1
      zcor(z,ij)=zcor(z,ij)/factor1
      end do
      end if
c
 55   continue
                         
      write(lu2,221)
      do 15 i=1,3+nc
      write (lu2,201) i, (xcor(x,kt(i)),x=2,11)
      write (lu2,202)    (xcor(x,kt(i)),x=12,nxhp)
 15   continue

      write(lu2,222)
      do 25 i=1,3+nc
      write (lu2,201) i, (ycor(y,kt(i)),y=2,11)
      write (lu2,202)    (ycor(y,kt(i)),y=12,nyhp)
 25   continue

      write(lu2,223)
      do 35 i=1,3+nc
      write (lu2,201) i, (zcor(z,kt(i)),z=2,11)
      write (lu2,202)    (zcor(z,kt(i)),z=12,nzhp)
 35   continue

c
c
        return
c
210     format(/,1x,' istep,time,3+nc=',i5,1p,e14.5,2x,i4)
211     format(1x,'one dimensional x-spectra')
212     format(1x,'one dimensional y-spectra')
213     format(1x,'one dimensional z-spectra')
221     format(1x,'one dimensional x-correlation')
222     format(1x,'one dimensional y-correlation')
223     format(1x,'one dimensional z-correlation')
c
 201  format (1x,i2,1p,10e13.5)
 202  format (3x,1p,10e13.5)
c
#endif
        end
