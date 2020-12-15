      subroutine vgmomt (utx,yp,icall,ithr)
c
c hybridized version by PKY, 1/1/2012
c 
c
	use comp
	include 'intvars'
c
c up to 8th order moments added, 12/18/01
c
c  routine called by sub. proc2a, if idflag=1, to handle
c  calculations of moments of du/dx and dw/dz
c
	real(b8) utx(nxpad,zisz,2)
c
      real(8) dudxm(8),dwdzm(8)
      real(8), allocatable :: sydudx(:,:),sydwdz(:,:)
	real(b8), allocatable :: sqterm(:)
      save sydudx,sydwdz,sqterm
	real ncube
c
c considered making 'sqterm' a 2D way with ithr as second
c subscript, but tests showed it was unnecessary.
c
	integer ithr,jthr,m
	real(8) sum
c
c
c for efficient parallel implementation, broken into 2 different kinds
c of calls: tally contributions from current y-plane if icall=0      
c form global moments if icall=1
c
      go to (1,2,3), icall+1
c
 1   	if (.not.allocated(sydudx)) then
        allocate (sydudx(8,0:num_thr-1))
        allocate (sydwdz(8,0:num_thr-1))
      allocate (sqterm(nxpad))
	end if
    	sydudx(:,:)=0.
	sydwdz(:,:)=0.
	go to 90
c
 2    continue
c
c in the no-scalar case, utr(x,1,z) now contains du/dx or dw/dz in physical
c space: calculate the contributions to first, second and third moments by
c direct summation
c
      sumu1=0.
      sumu2=0.
      sumu3=0.
      sumu4=0.
      sumu5=0.
      sumu6=0.
      sumu7=0.
      sumu8=0.
c
      sumw1=0.
      sumw2=0.
      sumw3=0.
      sumw4=0.
      sumw5=0.
      sumw6=0.
      sumw7=0.
      sumw8=0.
c
      do 66 zp=izp1(ithr),izp1(ithr)+izpi(ithr)-1
      do 66 x=1,nxpad
      sqterm(x)=utx(x,zp,1)*utx(x,zp,1)
      sumu1=sumu1+utx(x,zp,1)
      sumu2=sumu2+sqterm(x)
      sumu3=sumu3+sqterm(x)*utx(x,zp,1)
      sumu4=sumu4+sqterm(x)*sqterm(x)
      sumu5=sumu5+utx(x,zp,1)**5
      sumu6=sumu6+utx(x,zp,1)**6
      sumu7=sumu7+utx(x,zp,1)**7
      sumu8=sumu8+utx(x,zp,1)**8
      sqterm(x)=utx(x,zp,2)*utx(x,zp,2)
      sumw1=sumw1+utx(x,zp,2)
      sumw2=sumw2+sqterm(x)
      sumw3=sumw3+sqterm(x)*utx(x,zp,2)
      sumw4=sumw4+sqterm(x)*sqterm(x)
      sumw5=sumw5+utx(x,zp,2)**5
      sumw6=sumw6+utx(x,zp,2)**6
      sumw7=sumw7+utx(x,zp,2)**7
      sumw8=sumw8+utx(x,zp,2)**8
 66   continue
c
      sydudx(1,ithr)=sydudx(1,ithr)+sumu1
      sydudx(2,ithr)=sydudx(2,ithr)+sumu2
      sydudx(3,ithr)=sydudx(3,ithr)+sumu3
      sydudx(4,ithr)=sydudx(4,ithr)+sumu4
      sydudx(5,ithr)=sydudx(5,ithr)+sumu5
      sydudx(6,ithr)=sydudx(6,ithr)+sumu6
      sydudx(7,ithr)=sydudx(7,ithr)+sumu7
      sydudx(8,ithr)=sydudx(8,ithr)+sumu8
c
      sydwdz(1,ithr)=sydwdz(1,ithr)+sumw1
      sydwdz(2,ithr)=sydwdz(2,ithr)+sumw2
      sydwdz(3,ithr)=sydwdz(3,ithr)+sumw3
      sydwdz(4,ithr)=sydwdz(4,ithr)+sumw4
      sydwdz(5,ithr)=sydwdz(5,ithr)+sumw5
      sydwdz(6,ithr)=sydwdz(6,ithr)+sumw6
      sydwdz(7,ithr)=sydwdz(7,ithr)+sumw7
      sydwdz(8,ithr)=sydwdz(8,ithr)+sumw8
c
c
      go to 90
c
 3    continue
c
c icall=1: form global moments by a reduction operation across all tasks,
c        placing the result in task 0, which proceeds to write the info.
c
	do m=1,8
	sum=0.
	do jthr=0,num_thr-1
	sum=sum+sydudx(m,jthr)
	end do
	sydudx(m,0)=sum
	end do
	do m=1,8
	sum=0.
	do jthr=0,num_thr-1
	sum=sum+sydwdz(m,jthr)
	end do
	sydwdz(m,0)=sum
	end do
c
      call MPI_REDUCE (sydudx,dudxm,8,MPI_DOUBLE_PRECISION,MPI_SUM,0,
     1                 MPI_COMM_WORLD,ierr)
      call MPI_REDUCE (sydwdz,dwdzm,8,MPI_DOUBLE_PRECISION,MPI_SUM,0,
     1                 MPI_COMM_WORLD,ierr)
c
      if (taskid.ne.0) go to 90
c
c normalize
c
      ncube=float(nxpad)*nypad*nzpad
      do m=1,8
      dudxm(m)=dudxm(m)/ncube
      dwdzm(m)=dwdzm(m)/ncube
      end do
c
      if (jstep.eq.1) write (190,611) 
 611  format ('un-normalized moments of du/dx and dw/dz, up to 8th')
      write (190,612) (dudxm(m),m=1,8)
 612  format ('du/dx: ',1p,8e12.4)
      write (190,613) (dwdzm(m),m=1,8)
 613  format ('dw/dz: ',1p,8e12.4)
      
      do m=3,8
      dudxm(m)=dudxm(m)/dudxm(2)**(m/2.)
      dwdzm(m)=dwdzm(m)/dwdzm(2)**(m/2.)
      end do
c
      if (jstep.eq.1) write (191,621)
 621  format ('normalized moments of du/dx and dw/dz, 3rd to 8th')
      write (191,622) (dudxm(m),m=3,8)
 622  format ('du/dx: ',1p,6e12.4)
      write (191,623) (dwdzm(m),m=3,8)
 623  format ('dw/dz: ',1p,6e12.4)
      
      if (jstep.eq.1) write (18,603)
      write (18,601) istep,(dudxm(m),m=1,4)
 601  format ('istep=',i5,' moms of dudx=',1p,4e12.4)
      write (18,602) istep,(dwdzm(m),m=1,4)
 602  format ('istep=',i5,' moms of dwdz=',1p,4e12.4)
c
      write (6,604) ' dudx:',istep,dudxm(3),dudxm(4)
      write (6,604) ' dwdz:',istep,dwdzm(3),dwdzm(4)
 603  format ('skewness and flatness factor of dudx and dwdz')
 604  format (a5,2x,i6,1p,2e14.4)
c
 90    continue
c
      return
      end
