      subroutine vorsp2 (ifc,ilc,m)
c
	use comsp
	include 'intvars'
c
	integer ithr
c
c routine to manage computation of vorticity component spectra,
c in a similar way to sub. sptr2
c
c all normalization factors relating to grid size are now
c applied in sub. vorsp1
c
c
      if (ifc.eq.1) then
c
c on 1st call, set sums to zero
c
      do 1 l=1,mxyz
 1    vk(l)=0.
c
      do 2 ij=1,6
      do 2 l=1,mxyz
 2    vijk(l,ij)=0.
c
      do 3 ithr=0,num_thr-1
      do 3 ij=1,6
      do 3 l=1,mxyz
 3    vijky(l,ij,ithr)=0.
c
      return
c
      end if
c
c on last call, sum up the 3 components, compute variances
c and produce output
c
      if (ilc.eq.1) then
c
	call sumthr (vijky,mxyz,6,1,num_thr)
c
c sum up the z-slabs
c
      nn=6*mxyz
      call MPI_REDUCE (vijky(1,1,0),vijk(1,1),nn,mpireal,MPI_SUM,
     1                 0,MPI_COMM_WORLD,mpierr)
c
      if (taskid.ne.0) return
c
      do 20 k=1,mxyz
      vk(k)=vijk(k,1)+vijk(k,4)+vijk(k,6)
 20   continue
c
      ij=0
      do 30 i=1,3
      do 30 j=i,3
      ij=ij+1
      vvij(i,j)=0.
      do 35 k=1,mxyz
      vvij(i,j)=vvij(i,j)+vijk(k,ij)
 35   continue
      vij(ij)=vvij(i,j)
 30   continue
c
      end if
c
      return
c
      end
