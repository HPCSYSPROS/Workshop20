      subroutine vorout (lu)
 
	use comsp
	include 'intvars'
 
! write out vorticity spectra, like sub. eulout
 
      character*2 cdirn(6)
      data cdirn/'xx','xy','xz','yy','yz','zz'/
 
      vvij(2,1)=vvij(1,2)
      vvij(3,1)=vvij(1,3)
      vvij(3,2)=vvij(2,3)
      vortsq=vvij(1,1)+vvij(2,2)+vvij(3,3)
 
      st=shear*time
 
      write (lu,*)
      if (shear.eq.0.) then
      write (lu,201) istep,time,min(beta1,beta2,beta3)
      else
      write (lu,2101) istep,time,st
      end if
 
      write (lu,202) vortsq
      write (lu,203) vij(1),vij(2),vij(3)
      write (lu,204) vij(4),vij(5),vij(6)
 
      do 30 ij=1,6
      write (lu,205) cdirn(ij),(vijk(k,ij),k=1,min0(10,mxyz))
      if (mxyz.gt.10) write (lu,206) (vijk(k,ij),k=11,mxyz)
 30   continue
 
      write (lu,212) min(beta1,beta2,beta3)
      write (lu,206) (vk(k),k=1,mxyz)
 
 201  format (' istep=',i5,3x,'time=',1p,e13.5,  'shell thickness=',
     1        1p,e12.4)
 2101 format (' istep=',i5,3x,'time=',1p,e13.5,3x,'st=',1p,e13.5)
 
 202  format ('vorticity-vorticity covariance tensor:  enstrophy=',
     1        1p,e12.4)
 203  format ('v1v1, v1v2, v1v3: ',1p,3e12.4)
 204  format ('v2v2, v2v3, v3v3: ',1p,3e12.4)
 205  format (1x,a2,1p,10e13.5)
 206  format ((3x,1p,10e13.5))
 211  format ('components of 3-d vorticity spectrum tensor, ',
     1        'shell thickness=',e12.4)
 212  format ('trace of 3-d vorticity spectrum tensor, ',
     1        'shell thickness=',e12.4)
 
      return
      end
