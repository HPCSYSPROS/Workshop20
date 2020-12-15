      subroutine scgout 
c
#ifndef NOSCALAR
c
c output of scalar gradient covariances
c
	use comsp
c
      data icall/0/
      save icall
c
	integer, allocatable :: indsc(:)
c
	if (nc.eq.0) return

	allocate (indsc(nc))
c
c     if (icall.eq.0.and.istart.eq.0) then
      if (icall.eq.0.and.istart.le.0) then
      icall=1
      write (38,601) nc, (viscos/pr(i),i=1,nc)
      end if
c
c obtain correlation coefficients by normalization
c
      do 51 i=1,nc
 51   indsc(i)=1+(i-1)*(nc+1-i/2.)
c
      do 50 k=1,3
      ij=0
      do 50 i=1,nc
      do 50 j=i,nc
      ij=ij+1
	if ((kinit(3+i).eq.0.or.kinit(3+j).eq.0).and.istep.eq.1) then
	scgcor(ij,k)=0.
	else 
      scgcor(ij,k)=scgcov(ij,k)/
     1      sqrt(scgcov(indsc(i),k)*scgcov(indsc(j),k))
	end if
 50   continue
c
c
      if (shear.eq.0.) then
      write (38,201) istep-1,time
      else
      write (38,2011) istep-1,time,shear*time
      end if
	if (nc.le.4) then
      write (38,202) 'x',(scgcov(ij,1),ij=1,ncgd)
      write (38,202) 'y',(scgcov(ij,2),ij=1,ncgd)
      write (38,202) 'z',(scgcov(ij,3),ij=1,ncgd)
      write (38,202) 'x',(scgcor(ij,1),ij=1,ncgd)
      write (38,202) 'y',(scgcor(ij,2),ij=1,ncgd)
      write (38,202) 'z',(scgcor(ij,3),ij=1,ncgd)
	end if
c
 201  format (/'istep,time=',i5,1p,e12.4)
 2011 format (/'istep,time=',i5,1p,e12.4,'   st=',1p,e12.4)
c 202  format (a1,1x,1p,6e12.5)
 202  format (a1,1x,1p,10e12.5)
c
 601  format ('joint covariances and corr. coeffs of scalar gradients'/
c    1        'number of scalars=',i2,17x,1p,3e11.3)
     1        'number of scalars=',i2,17x,1p,6e11.3)
c
	deallocate (indsc)
c
#endif
      return
      end
