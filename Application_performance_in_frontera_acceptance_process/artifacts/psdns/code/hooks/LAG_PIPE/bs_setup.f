	subroutine bs_setup
c
#ifdef LAG

	use mpilag
	implicit none
	integer i,j,ii,jj
c
c! 
c! Next lines needed for cubic splines
c!
c!Mapping 3-D data arrays onto 2-D process grid
c! (nx+2,ny+2,nz) => (iproc,jproc)
c!
	if (taskid.eq.1) write (6,*) 'enter bs_setup'
       allocate (bxist(0:iproc-1))
       allocate (bxisz(0:iproc-1))
       allocate (bxien(0:iproc-1))
       allocate (bzjst(0:jproc-1))
       allocate (bzjsz(0:jproc-1))
       allocate (bzjen(0:jproc-1))
       allocate (bzist(0:iproc-1))
       allocate (bzisz(0:iproc-1))
       allocate (bzien(0:iproc-1))
c!

 	nbx=nx+3
 	nby=ny+3
 	nbz=nz+3
c
       call MapDataToProc(nbx,iproc,bxist,bxien,bxisz)
       call MapDataToProc(nbz,iproc,bzist,bzien,bzisz)
       call MapDataToProc(nbz,jproc,bzjst,bzjen,bzjsz)

c!
        bxistart = bxist(ipid)
        bzistart =bzist(ipid)
        bzjstart = bzjst(jpid)
 
        bxisize= bxisz(ipid)
        bzjsize= bzjsz(jpid)
        bzisize= bzisz(ipid)
 
        bxiend = bxien(ipid)
        bzjend = bzjen(jpid)
        bziend = bzien(ipid)
 
c
!  mpi derived data types
        allocate(bcnts_zi(0:iproc-1))
        allocate(bdisp_zi(0:iproc-1))
        allocate(bcnts_xi(0:iproc-1))
        allocate(bdisp_xi(0:iproc-1))
 
        allocate(bcnts_zj(0:jproc-1))
        allocate(bdisp_zj(0:jproc-1))
        allocate(bcnts_yj(0:jproc-1))
        allocate(bdisp_yj(0:jproc-1))
 
!        write (7000+taskid,*) 'ipid,bxisize=',ipid,bxisize
         do i=0,iproc-1
		ii = mymap(i)
           bcnts_zi(i) = bxisz(ii)* zisz
           bdisp_zi(i) = (bxist(ii)-1)* zisz
           bcnts_xi(i) = kisz(ii) * bxisize
           bdisp_xi(i) = (kist(ii)-1)* bxisize
!      write (7000+taskid,"(i4,4i8)") i,bcnts_zi(i),bdisp_zi(i),
!     1              bcnts_xi(i),bdisp_xi(i)
        enddo
!	close (7000+taskid)
c
        do j=0,jproc-1
           bcnts_yj(j) = bzjsz(j) * yjsz *bxisize
           bdisp_yj(j) = (bzjst(j)-1)*yjsz*bxisize
           bcnts_zj(j) = jjsz(j)* bzjsize*bxisize
           bdisp_zj(j) = (jjst(j)-1)* bzjsize*bxisize
        enddo
c
	if (taskid.eq.1) write (6,*) ' exit bs_setup'

#endif
	end subroutine
