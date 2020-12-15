        subroutine hist1vip ( ifc,ilc,kwt, u,wt,n,stride,
     1                       nbin, bin, hist, fmt )
c
c 04/05: non-iterative algorithm suggested by SB Pope
c
c 11/97: version of hist1 with prescribed variable interval sizes
c
	use mpicom, only: taskid
	use precision, only: b8
        integer stride
        real(b8) u(1),wt(1)
c        real(b8) u(2),wt(1)
        real(b8) hist(nbin)
	real(b8) bin(nbin)
        character*(*) fmt
c
	real(b8), allocatable :: bl(:),br(:),bin2(:),bl2(:),br2(:)
	integer, allocatable :: jj(:,:),nj(:)
	
	save bin2,jj,nj,nbin2,bl,br,bl2,br2
c
c
        ilast = 1 + (n-1)*stride
c
c clean-up call
c
	if (ifc.eq.-1) then
	deallocate (bl,br)
	if (nbin2.le.2000) deallocate (bin2,bl2,br2,jj,nj)
	return
	end if
c
c  initialization
c
      if (ifc.eq.1) then
c
c
      do ku=1,nbin
      hist(ku)=0.
      end do
c
c find the minimum interval spacing
c
      dxmin=1000.
      do i=1,nbin-1
      dx=bin(i+1)-bin(i)
      dxmin=amin1(dx,dxmin)
      end do
c
	if (.not.allocated(bl)) then
	allocate (bl(nbin),br(nbin))
	end if
c
c
c fix on 3/15/06: relevant for cases where some intervals
c are of the same (minimum) size
c
        dxmin=dxmin*.99
c
c form the secondary "j-bins" with uniform spacing
c
c 04/11/10: raw value of nbin2 may exceed single-precision ceiling
c if range is very wide and number of bins is not very large
c
      nbin2=(bin(nbin)-bin(1))/dxmin+1
	if (nbin2.lt.0) nbin2=10000
c
c 10/18/09: if either of nbiu2 is really large
c do not use but switch back to the basic scheme
c
c	if (taskid.le.1) write (6,*) 'hist1vip1: nbin2=',nbin2,bin(nbin),bin(1),dxmin
c
        if (nbin2.gt.2000) go to 111

      dbin2=(bin(nbin)-bin(1))/(nbin2-1)
	if (.not.allocated(bin2)) then
	allocate (bin2(nbin2),bl2(nbin2),br2(nbin2))
	allocate (nj(nbin2),jj(2,nbin2))
	end if
c

      do j=1,nbin2
      bin2(j)=(j-1)*dbin2+bin(1)
      end do
c
c define the centered intervals
c
      do i=2,nbin-1
      bl(i)=.5*(bin(i-1)+bin(i))
      br(i)=.5*(bin(i+1)+bin(i))
      end do
      bl(1)=bin(1)
      br(1)=bl(2)
      bl(nbin)=br(nbin-1)
      br(nbin)=bin(nbin)
c
      do j=2,nbin2-1
      bl2(j)=.5*(bin2(j-1)+bin2(j))
      br2(j)=.5*(bin2(j+1)+bin2(j))
      end do
      bl2(1)=bin2(1)
      br2(1)=bl2(2)
      bl2(nbin2)=br2(nbin2-1)
      br2(nbin2)=bin2(nbin2)

      do j=1,nbin2
      nj(j)=0
      do i=1,nbin
      if (bl2(j).le.br(i).and.br2(j).ge.bl(i)) then
      nj(j)=nj(j)+1
      jj(nj(j),j)=i 
c     write (103,*) j,i,bl2(j),br2(j),bl(i),br(i)
      end if
      end do
c     write(102,900)j,bin2(j),bl2(j),br2(j),nj(j),
c    1                (jj(ij,j),ij=1,nj(j))
c900  format (i4,1p,3e12.4,4i4)
      end do
c
	go to 112
c
      end if
c
 111	continue
c
c define the intervals
c
      bl(1)=bin(1)
      br(1)=(bin(1)+bin(2))/2.
      do i=2,nbin-1
      bl(i)=(bin(i-1)+bin(i))/2.
      br(i)=(bin(i+1)+bin(i))/2.
      end do
      bl(nbin)=(bin(nbin-1)+bin(nbin))/2.
      br(nbin)=bin(nbin)
c 
 112	continue
c
        if (n.eq.0.and.ilc.eq.0) return
c
	if (nbin2.le.2000) then
	bldd=bin2(1)
	brdd=bin2(nbin2)
	umult=(nbin2-1)/(brdd-bldd)
	end if
c
c  form histogram
c
        if( kwt .eq. 0 ) then
c
           do 100 i=1,ilast,stride
      	do 10 m=1,nbin
      	if (u(i).ge.bl(m).and.u(i).lt.br(m)) go to 15
 10   	continue
 15    	hist(m) = hist(m) + wt(i)
 100	   continue 
c
        else
c
	if (nbin2.le.2000) then
c
           do 200 i=1,ilast,stride
	ku=(u(i)-bldd)*umult+1.5
	ku=max0(1,min0(ku,nbin2))
	m=jj(1,ku)
        if (nj(ku).eq.2) then
	if (u(i).ge.bl(jj(2,ku)).and.u(i).lt.br(jj(2,ku)))
     1  m=jj(2,ku)
	end if
       	hist(m) = hist(m) + 1.
  200	   continue 
c
	else
c
        do 220 i=1,ilast,stride
        do m=1,nbin
        if (u(i).ge.bl(m).and.u(i).lt.br(m)) go to 35
        end do
        m=nbin
 35     ku=m
        if (u(i).lt.bl(1)) ku=1
	hist(ku)=hist(ku)+1
 220	continue
c
	end if

c
        end if
c
c  on last call, normalize
c
        if( ilc .eq. 0 ) return
c
c
c add up total number of samples
c
      sum=0.
      do ku=1,nbin
      sum=sum+hist(ku)
      end do
c
c normalize by total number of samples and interval size,
c with special treatment for first and last nodes
c
           if( sum .ne. 0. ) then
       do 30 ku=2,nbin-1
       hist(ku)=hist(ku)/sum/(br(ku)-bl(ku))
 30    continue
       hist(1)=hist(1)/sum/(bin(2)-bin(1))
       hist(nbin)=hist(nbin)/sum/(bin(nbin)-bin(nbin-1))
           end if
c
c
        if( ilc .eq.-1 ) return
c
c  write out data
c
      lu=iabs(ilc)
c
      bldd=bin(1)
      brdd=bin(nbin)
c
      write (lu,201) 1,0
      write (lu,202) nbin,bldd,brdd
      write (lu,fmt) bin
      write (lu,fmt) hist
c	close (lu)
c
 201  format (2i3)
 202  format (i5,1p,2e12.4)
c
        return
        end
