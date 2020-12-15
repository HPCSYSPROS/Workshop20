      subroutine spcal( ibxyz, bfpos,bsxy,f,sf,np,npdim,iod)

#ifdef LAG
c version for 2D code
c
c routine to perform the basis function summations for
c all particles and thus store the particle properties
c
c this is a completely vectorised version,
c and the loop ordering (k,m,j,i) is found to run fastest,
c presumably because of reduced paging costs in the array bs
c
c this may be regarded as the analogue of sub. interp
c
c bfxyz:  array of basis functions determined by sub. intbf
c ibxyz:  indices of those basis functions
c
c f(m):   spline-interpolated property of mth particle
c sf:     scale factor, equals sqrt(b11(kstep)) etc for velcocities,
c         unity for all other properties
c np:     no. of particles
c npdim:  declared upper limit of np
c
c swork1 & swork2:  working variables in the summation procedure
c
c iod(n) (n=1,2,3) is 0 (function) or 1 (derivative) in direction n
c
c generalised to include use of spline derivatives in ns code
c
c
	use mpilag
c
      parameter (d6=1./6., tt=2./3.)
c
      real bsxy(bxisize,nby,bzjsize)
      integer iod(3)
      integer ibxyz
	integer ibx,iby,ibz
	real bfx,bfy,bfz
      real f(npdim)
      real bff(4,3)
      real bfpos(npdim,3)
c

      real, allocatable :: work(:)
c
      integer xp,zp
c
#ifdef DEC2008
	mp=np/numtasks
#else
	mp=np
#endif

      allocate (work(mp),stat=ierr)

c for function or 1st derivative of spline, just pick the
c correct basis functions to use in the summation, based
c on the order of differentiation
c
      ib1=1+iod(1)*3
      ib2=2+iod(2)*3
      ib3=3+iod(3)*3
c
	ix1=bxistart
	ix2=bxiend
	iz1=bzjstart
	iz2=bzjend
c
#ifdef DEC2008
c
 	do 100 itask=1,numtasks
c
 	mp1=(itask-1)*mp+1
	mp2=mp1+mp-1
c
	work(:)=0.
c
      do 15 m=mp1,mp2
c
c	ip=mp1-m+1
	ip=m-mp1+1
c
	ibx=bfpos(m,1)
	iby=bfpos(m,2)
	ibz=bfpos(m,3)
	bfx=bfpos(m,1)-ibx
	bfy=bfpos(m,2)-iby
	bfz=bfpos(m,3)-ibz
c
 	if (ibx.lt.ix1-4.or.ibx.gt.ix2-1) go to 15
 	if (ibz.lt.iz1-4.or.ibz.gt.iz2-1) go to 15
c
c
c     if (ibxyz(m,1).lt.ix1-4.or.ibxyz(m,1).gt.ix2-1) go to 15
c     if (ibxyz(m,3).lt.iz1-4.or.ibxyz(m,3).gt.iz2-1) go to 15
c
	
c compute and store basis functions
c
      pos = bfx
      ym = 1. - pos
      yms = ym*ym
      bff(1,1) = ym*yms*d6
      bff(2,1) = tt + 0.5*pos*(yms - 1.)
      bff(4,1) = pos*pos*pos*d6
      bff(3,1) = 1. - bff(1,1) - bff(2,1) - bff(4,1)
      pos = bfy
      ym = 1. - pos
      yms = ym*ym
      bff(1,2) = ym*yms*d6
      bff(2,2) = tt + 0.5*pos*(yms - 1.)
      bff(4,2) = pos*pos*pos*d6
      bff(3,2) = 1. - bff(1,2) - bff(2,2) - bff(4,2)
      pos = bfz
      ym = 1. - pos
      yms = ym*ym
      bff(1,3) = ym*yms*d6
      bff(2,3) = tt + 0.5*pos*(yms - 1.)
      bff(4,3) = pos*pos*pos*d6
      bff(3,3) = 1. - bff(1,3) - bff(2,3) - bff(4,3)
c
      do 30 k=1,4
c
      kk=ibz+k
      if (kk.lt.iz1.or.kk.gt.iz2) go to 30
      zp=kk-iz1+1
      swork1=0.
c
      do 10 i=1,4
c
      ii=ibx+i
      if (ii.lt.ix1.or.ii.gt.ix2) go to 10
      xp=ii-ix1+1
      swork2=0.
c
      do 20 j=1,4
      jj=iby+j
      swork2 = swork2 + bsxy(xp,jj,zp) * bff(j,ib2)
 20   continue
c
      swork1=swork1+swork2*bff(i,ib1)
c
 10   continue
c
      work(ip)=work(ip)+swork1*bff(k,ib3)
 30   continue
c
 15   continue
c
c this mp_combine is time-consuming: however it is implemented as
c an mp_reduce followed by mp_bcast
c
 	call MPI_ALLREDUCE (work,f(mp1),mp,MPI_REAL,MPI_SUM,
     1                   MPI_COMM_WORLD,mpierr)
c
 100	continue
c
#else
c
      do 15 m=1,mp
c
	ibx=bfpos(m,1)
	iby=bfpos(m,2)
	ibz=bfpos(m,3)
	bfx=bfpos(m,1)-ibx
	bfy=bfpos(m,2)-iby
	bfz=bfpos(m,3)-ibz
c
	if (ibx.lt.ix1-4.or.ibx.gt.ix2-1) go to 15
	if (ibz.lt.iz1-4.or.ibz.gt.iz2-1) go to 15
c
c     if (ibxyz(m,1).lt.ix1-4.or.ibxyz(m,1).gt.ix2-1) go to 15
c     if (ibxyz(m,3).lt.iz1-4.or.ibxyz(m,3).gt.iz2-1) go to 15
c
	
c compute and store basis functions
c
      pos = bfx
      ym = 1. - pos
      yms = ym*ym
      bff(1,1) = ym*yms*d6
      bff(2,1) = tt + 0.5*pos*(yms - 1.)
      bff(4,1) = pos*pos*pos*d6
      bff(3,1) = 1. - bff(1,1) - bff(2,1) - bff(4,1)
      pos = bfy
      ym = 1. - pos
      yms = ym*ym
      bff(1,2) = ym*yms*d6
      bff(2,2) = tt + 0.5*pos*(yms - 1.)
      bff(4,2) = pos*pos*pos*d6
      bff(3,2) = 1. - bff(1,2) - bff(2,2) - bff(4,2)
      pos = bfz
      ym = 1. - pos
      yms = ym*ym
      bff(1,3) = ym*yms*d6
      bff(2,3) = tt + 0.5*pos*(yms - 1.)
      bff(4,3) = pos*pos*pos*d6
      bff(3,3) = 1. - bff(1,3) - bff(2,3) - bff(4,3)
c
      do 30 k=1,4
c
      kk=ibz+k
      if (kk.lt.iz1.or.kk.gt.iz2) go to 30
      zp=kk-iz1+1
      swork1=0.
c
      do 10 i=1,4
c
      ii=ibx+i
      if (ii.lt.ix1.or.ii.gt.ix2) go to 10
      xp=ii-ix1+1
      swork2=0.
c
      do 20 j=1,4
      jj=iby+j
      swork2 = swork2 + bsxy(xp,jj,zp) * bff(j,ib2)
 20   continue
c
      swork1=swork1+swork2*bff(i,ib1)
c
 10   continue
c
      work(m)=work(m)+swork1*bff(k,ib3)
 30   continue
c
 15   continue
c
c this mp_combine is time-consuming: however it is implemented as
c an mp_reduce followed by mp_bcast
c
      call MPI_ALLREDUCE (work,f,np,MPI_REAL,MPI_SUM,
     1                   MPI_COMM_WORLD,mpierr)
c
#endif
 
c
      deallocate (work,stat=ierr)
c
      if (abs(sf-1.).gt.1.e-6) f(:)=f(:)*sf


      return


#endif
      end

