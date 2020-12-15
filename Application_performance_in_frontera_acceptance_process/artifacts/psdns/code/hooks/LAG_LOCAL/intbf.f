      subroutine intbf( xyzl,nxyz,shif, xp,np,npdim, ibxyz, bfpos,slant )

#ifdef LAG
c
c  routine to determine indicial locations of particles in the
c  grid and to compute the basis function coefficients
c  (for the function)
c
c  extended for slanted box occurring in homogeneous shear, 8/10/94
c
c  this may be regarded as the analogue of sub. intwt
c
c  input:
c      xyzl(3)   - x,y,z lengths of the solution domain.
c      nxyz(3)   - number of grid nodes in each direction.
c      shif(3)  - grid shift (i.e. the node (1,1,1) has
c                    coordinates shif(1)+dx,...etc)
c      xp(np,3)  - particle coordinates.
c      np        - number of particles.
c      npdim     - first dimension of xp.
c      ipbc      - 1 if particles are inside the box
c                - 2 if particles may be outside, and 
c                    periodicity to be accounted for
c
c  output:
c      bfpos(i,j) the 1st of 4 basis functions in the jth
c                   direction for the ith particle
c      ibxyz(i,j) index of the 1st of those basis functions
c      (note that ibxyz(i,j) gives the j-index of the node to
c       the immediate "left" of the ith particle, in jth direction)
c
	use mpicom
c
      parameter (d6=1./6., tt=2./3.)
c
      integer nxyz(3),ibxyz
      dimension xyzl(3),shif(3),xp(npdim,3)
      dimension bfpos(npdim,3)
c
c
c note that items in next line must be declared double precision:
c otherwise some of the ibxyz and bfpos's may not be correct
c
      real*8 gfn(3),denom,ploc,dxyz(3),box(3)
c
c dxyz are sizes of grid spacing in particle coordinates
c (this would be unity for dx=dy=dz with nx=ny=nz)
c
c gfn(i) is coord of first node in i-th dirn
c
      do k=1,3
      dxyz(k)=xyzl(k)/nxyz(k)
      gfn(k)=shif(k)+1.
      box(k)=nxyz(k)
      end do
c
      ipbc=2
c
	toler=.00001
c
c separate loop for x-direction, allowing for slanted box
c
      denom=1./dxyz(1)
c
      work=100*xyzl(1)-gfn(1)
c
c loop over particles
c
#ifdef DEC2008
	mp=np/numtasks
	mp1=taskid*mp+1
	mp2=mp1+mp-1
#else
	mp1=1
	mp2=np
#endif
c
      do 22 i=mp1,mp2
c
c locate the particle and
c calculate local normalised coord (between 0 and 1)
c
      r=slant*(xp(i,2)-gfn(2))
c
      px=xp(i,1)+work-r
      px=amod(px,xyzl(1))+gfn(1)+r
      ploc=(px-gfn(1)-r)*denom+1.
c
      jloc=ploc+toler
      bfpos(i,1)=ploc-jloc
c
c     ibxyz(i,1)=jloc-1
      bfpos(i,1)=bfpos(i,1)+jloc-1
c
 22   continue
c
c y- and z- directions
c
      do 20 j=2,3
c
      denom=1./dxyz(j)
c
      work=100*xyzl(j)-gfn(j)
c
c loop over particles
c
      do 25 i=mp1,mp2
c
c locate the particle and
c calculate local normalised coord (between 0 and 1)
c
      px=xp(i,j)+work
      px=amod(px,xyzl(j))+gfn(j)
      ploc=(px-gfn(j))*denom+1.
c
      jloc=ploc+toler
      bfpos(i,j)=ploc-jloc
c
c     ibxyz(i,j)=jloc-1
      bfpos(i,j)=bfpos(i,j)+jloc-1
c
 25   continue
c
 20   continue
c
#ifdef DEC2008
	do j=1,3
	call MPI_ALLGATHER (bfpos(mp1,j),mp,MPI_REAL,
     1                      bfpos(1,j),mp,MPI_REAL,
     1                      MPI_COMM_WORLD,mpierr)
	end do
#endif
     
c
#endif
      end
