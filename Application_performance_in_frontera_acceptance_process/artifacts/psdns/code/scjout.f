      subroutine scjout (i1,i2,lu)
 
#ifndef NOSCALAR
 
	use comsp
 
! routine to print co-spectra of scalars i1 and i2
 
      ifst=min0(i1,i2)
      isec=max0(i1,i2)
      ij=kt(ifst+3)+isec-ifst
 
      if (jstep.eq.0) write (lu,211) viscos,pr(i1),pr(i2)
 
      fmult=.5*(1./pr(i1)+1./pr(i2))
      do 10 k=1,nxh
 10   dijk(k,ij)=dijk(k,ij)*fmult
 
      if (shear.eq.0.) then
      write (lu,201) istep,time
      else
      write (lu,2011) istep,time,shear*time
      end if
      write (lu,202) (eijk(k,ij),k=1,nxh)
      write (lu,202) (dijk(k,ij),k=1,nxh)
 
 201  format (/'co-spectra at step',i5,'   time=',1p,e12.4)
 2011 format (/'co-spectra at step',i5,'   time=',1p,e12.4,
     1         '   st=',1p,e13.5)
 202  format (1p,10e13.5)
 211  format ('viscosity and prandtl nos.=',1p,3e12.4)
 
#endif
 
      return
      end
