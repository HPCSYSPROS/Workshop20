      subroutine pitest (bsy,k)
c
#ifdef LAG
c
c cubic-spline interpolation test for particles,
c performed at first time step if shear=0.,
c but every time step for shear flow
c
	use compart
c
	real bsy(bxisize,nby,bzjsize)
	real(b8), allocatable  :: posnode(:,:)
c
        real, allocatable :: ppnode(:,:),bfnode(:,:)
	integer ibnode
c
c ntp= no. of test particles per MPI task
	integer itp,ntp
c
	real, allocatable :: bfw(:,:)
c
c
	ntp=1
c
      ifail=0
c
c use one test particle per processor
c
        npdim=numtasks*ntp
	allocate (posnode(ntp,3))
c
        allocate (ppnode(1,6),bfnode(3,npdim))
c
c set test particle positions
c
c      toler=.00001
c      toler=.0000001
	toler=0.
c     ppwork1=5+gsh(1,1)+toler
c     ppwork2=3+jjst(jpid)+gsh(2,1)+toler 
c     ppwork3=4+kist(ipid)+gsh(3,1)+toler
      ppwork1=1+gsh(1,1)+toler
      ppwork2=jjst(jpid)+gsh(2,1)+toler 
      ppwork3=kist(ipid)+gsh(3,1)+toler
	do itp=1,ntp
	posnode(itp,1)=ppwork1
	posnode(itp,2)=ppwork2
	posnode(itp,3)=ppwork3
	end do
c
c
c find indicial positions of particles, and compute basis functions
c at the particle locations
c
      call intbf_ars ( xyzl, nxyz, gsh(1,1), posnode, npdim, npdim, bfnode, 0.,
0)
c
c obtain interpolated values by
c taking the triple tensor product between the basis functions
c and their coefficients
c
c
      call spcal_ars ( bfnode, bsy(1,1,1),
     1             ppnode(1,3+k), 1, 1., npdim, npdim, iod, 0 )
c
c check against the grid point values
c
        do 10 ip=1,npdim
c
        ibx=bfnode(1,ip)
        iby=bfnode(2,ip)
        ibz=bfnode(3,ip)
        bfx=bfnode(1,ip)-ibx
        bfy=bfnode(2,ip)-iby
        bfz=bfnode(3,ip)-ibz
c
        ixloc=ibx+1
        iyloc=iby+1
        izloc=ibz+1
        if (bfx.gt.1.-toler) ixloc=ixloc+1
        if (bfy.gt.1.-toler) iyloc=iyloc+1
        if (bfz.gt.1.-toler) izloc=izloc+1
c
        if (iyloc.ge.yjst.and.iyloc.le.yjen.
     1  and.izloc.ge.zist.and.izloc.le.zien) then
c
        iy=iyloc-jjst(jpid)+1
        iz=izloc-kist(ipid)+1
	ugrid=u(ixloc,iz,iy,k)
c
       write (6,603) ixloc,iz,iy,taskid,ugrid,ppnode(1,3+k)

c	write (6,903) ip,(ppnode(1,j),j=1,3),ibx,iby,ibz
c903	format ('ip,ppnode,ibx=',i3,1p,3e12.4,3i4)
c	write (6,904) ip,(bfnode(ip,j),j=1,3),ibx,iby,ibz
c904	format ('ip,bfnode,ibx=',i3,1p,3e12.4,3i4)
c
	difabs=abs(ugrid-ppnode(1,3+k))
c	if (abs(difabs/ugrid).gt.0.02.or.difabs.gt.0.02) then
	if (abs(difabs/ugrid).gt.0.02.and.difabs.gt.0.02) then
 	ifail=ifail+1
      write (6,610) taskid,ipid,jpid,kist(ipid),jjst(jpid)
      write (6,*) 'pitest: ppnode=',(ppnode(1,j),j=1,3)
      write (6,*) 'pitest: ibnode=',ibx,iby,ibz
      write (6,*) 'pitest: bfnode=',bfx,bfy,bfz
 610  format ('interpolation test fails: >2% error at nodal point',
     1        ' taskid=',3i3,2i4)
	end if
c
        end if
 10	continue
c

c
 603  format ('test particle: taskid,u(',i3,',',i3,',',i3,
     1         '), vel=',i5,1p,2e14.6)

c
c
      if (ifail.ge.1) call MPI_ABORT (MPI_COMM_WORLD,ierror)
c
c	stop 'stop at exit of pitest'

c
#endif
      return
      end
