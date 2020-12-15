      subroutine scdlpdf (ifc,ilc,i,yp,grx,gry,grz)
c
c get standardized PDF of logarithm of scalar dissipation
c
#ifndef NOSCALAR
c
	use comsp
	include 'intvars'
c
c values in physical space, called from sub. proc2a
c
c also skewness, flatness factor and superskewness
c
        real(b8) grx(nx,zisz),gry(nx,zisz),grz(nx,zisz)
!
	parameter (nddu=81)
c
	real(b8), allocatable :: hist(:,:),bldd(:),brdd(:)
	double precision, allocatable :: chisum(:,:)
	save hist,chisum,bldd,brdd
c
        integer icall
        save icall
        data icall/0/
!
	real(b8) hist0(nddu)
	real*8 scdmom(6)
	real(b8) nn
	real mean
c
      real scdlog(nx),chi(nx)
      real*8 scmom(6),chimom(6)
      real udum(1),wdum(1)
      character*7 head
      character*1 tail
      character*8 fn
      character*12 fmt
	character*40 name
c
      data head/'scdlpdf'/
      data fmt/'(1p,7e11.4)'/
c
      if (icall.eq.0) then
      icall=icall+1
        allocate (hist(nddu,nc))
        allocate (chisum(6,nc))
	allocate (bldd(nc),brdd(nc))
        if (taskid.eq.0) then
        do ic=1,nc
      lu=luscd+(ic-1)*6
        write (tail,"(i1)") ic
        fn=head//tail
        call fopen1 (lu,fn,'formatted  ')
        end do
        end if
      end if
c
c initialize on first plane of the current slab
c
      if (ifc.eq.1) then
	chisum(:,i)=0.
	bldd(i)=-11.
	brdd(i)=7.
      call hist1 (1,0,1,udum,wdum,0,1,nddu,bldd(i),brdd(i),
     1            hist(1,i),fmt)
      end if
c
      diffus=viscos/pr(i)
      const=2.*diffus
c
      do 10 zp=1,zisz
c
      do 20 x=1,nx
      scdlog(x)=grx(x,zp)**2+gry(x,zp)**2+grz(x,zp)**2
      scdlog(x)=alog(const*scdlog(x))
      do k=1,6
      chisum(k,i)=chisum(k,i)+scdlog(x)**k
      end do
 20   continue
c
      call hist1 (0,0,1,scdlog,wdum,nx,1,nddu,bldd(i),brdd(i),
     1            hist(1,i),fmt)
c
 10   continue
c
c collect tally over all tasks
c
      if (ilc.eq.1) then
      call MPI_REDUCE (hist(1,i),hist0,nddu,MPI_REAL,MPI_SUM,0,
     1                 MPI_COMM_WORLD,mpierr)
      call MPI_REDUCE (chisum(1,i),scdmom,6,MPI_DOUBLE_PRECISION,
     1                 MPI_SUM,0,MPI_COMM_WORLD,mpierr)
      end if
c
! output by task 0
                                                                                
      if (ilc.ne.1.or.taskid.ne.0) return
                                                                                
c
      nn=float(nx)*ny*nz
      do k=1,6
      scdmom(k)=scdmom(k)/nn
      end do
c
      mean=scdmom(1)
      varcl=scdmom(2)-mean**2
      skcl=(scdmom(3)-3.*mean*scdmom(2)+2.*mean**3)/varcl**1.5
      flatcl=(scdmom(4)-4.*mean*scdmom(3)+6.*mean**2*scdmom(2)-3.*mean**4)
     1     /varcl**2
      sskcl=(scdmom(6)-6.*mean*scdmom(5)+15.*mean**2*scdmom(4)
     1   -20.*mean**3*scdmom(3)+15.*mean**4*scdmom(2)-5.*mean**6)/varcl**3
c
      lu=luscd+(i-1)*6
      write (lu,211) 0,0,40,fmt,mean,varcl,skcl,flatcl
      write (lu,212) istep-1,time,sskcl
c
c re-scale to standardize logarithm
c
	bldd(i)=(bldd(i)-mean)/sqrt(varcl)
	brdd(i)=(brdd(i)-mean)/sqrt(varcl)
c
      call hist1 (0,lu,1,udum,wdum,0,1,nddu,bldd(i),brdd(i),hist0,fmt)
c
c
      return
c
 211  format (3i3,3x,a12,3x,'moments:',1p,4e10.3)
 212  format ('log scalar diss. step',i5,' time=',1p,e9.3,
     1         5x,' sskew=',1p,e10.3)
c
#endif
c
      end
