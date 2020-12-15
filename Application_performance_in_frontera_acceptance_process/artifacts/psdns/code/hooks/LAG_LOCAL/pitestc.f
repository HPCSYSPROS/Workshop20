      subroutine pitestc (bsy,k)
c
#ifdef LAG
#ifndef NOSCALAR
c
c adopted from new version of pitest.f, Jan 2008
c
c cubic-spline interpolation test for particles,
c performed at first time step if shear=0.,
c but every time step for shear flow
c
	use compart
c
	real bsy(bxisize,nby,bzjsize)
c
        real, allocatable :: ppnode(:,:),bfnode(:,:)
        integer, allocatable :: ibnode(:,:)
c
        integer, allocatable :: kcall(:)
	integer icall
	save icall,kcall
	data icall/0/
c
	if (icall.eq.0) then 
	allocate (kcall(3+2*nc))
	kcall(:)=0
	icall=1
	end if
c
c
	if (kcall(k).gt.0) go to 90
	kcall(k)=1
c
      ifail=0
c
c use one test particle per processor
c
        npdim=numtasks
        allocate (ppnode(npdim,4),bfnode(npdim,3),ibnode(npdim,3))
c
c set test particle positions
c
      toler=.00001
      ppwork1=5+gsh(1,1)+toler
      ppwork2=3+jjst(jpid)+gsh(2,1)+toler
      ppwork3=4+kist(ipid)+gsh(3,1)+toler
      call MPI_ALLGATHER (ppwork1,1,mpireal,ppnode(1,1),1,mpireal,
     1                 MPI_COMM_WORLD,mpierr)
      call MPI_ALLGATHER (ppwork2,1,mpireal,ppnode(1,2),1,mpireal,
     1                 MPI_COMM_WORLD,mpierr)
      call MPI_ALLGATHER (ppwork3,1,mpireal,ppnode(1,3),1,mpireal,
     1                 MPI_COMM_WORLD,mpierr)
c
c find indicial positions of particles, and compute basis functions
c at the particle locations
c
      call intbf( xyzl, nxyz, gsh(1,1), ppnode, npdim, npdim, ibnode,bfnode, 0.)
c
c obtain interpolated values by
c taking the triple tensor product between the basis functions
c and their coefficients
c
c
      call spcal ( ibnode, bfnode, bsy(1,1,1),
     1             ppnode(1,4), 1., npdim, npdim, iod )
c
c check against the grid point values
c
        do 10 ip=1,npdim
c
        ixloc=ibnode(1,ip)+1
        iyloc=ibnode(2,ip)+1
        izloc=ibnode(3,ip)+1
        if (bfnode(1,ip).gt.1.-toler) ixloc=ixloc+1
        if (bfnode(2,ip).gt.1.-toler) iyloc=iyloc+1
        if (bfnode(3,ip).gt.1.-toler) izloc=izloc+1
c
        if (iyloc.ge.yjst.and.iyloc.le.yjen.
     1  and.izloc.ge.zist.and.izloc.le.zien) then
c
        iy=iyloc-jjst(jpid)+1
        iz=izloc-kist(ipid)+1
	ugrid=u(ixloc,iz,iy,k)
       write (6,603) ixloc,iz,iy,taskid,k,ugrid,ppnode(ip,4)
c
	difabs=abs(ugrid-ppnode(ip,4))
	if (abs(difabs/ugrid).gt.0.02.or.difabs.gt.0.02) then
 	ifail=ifail+1
      write (6,*) 'pitestc: interpolation testc fails: more than',
     1            ' 2% error at nodal point'
      write (6,*) 'pitestc: ppnode=',ppnode(ip,1),ppnode(ip,2),ppnode(ip,3)
      write (6,*) 'pitestc: ibnode=',ibnode(ip,1),ibnode(ip,2),ibnode(ip,3)
      write (6,*) 'pitestc: bfnode=',bfnode(ip,1),bfnode(ip,2),bfnode(ip,3)
	end if


        end if
 10	continue

 603  format ('test particle: taskid,u(',i3,',',i3,',',i3,
     1         '), grid value=',2i3,1p,2e14.6)

	if (taskid.eq.0) write (6,*) 'exit pitestc'
c
        deallocate (ppnode,bfnode,ibnode)
      if (ifail.ge.1) call MPI_ABORT (MPI_COMM_WORLD,ierror)
c
#endif
#endif
 90   return
      end
