      subroutine mpisetup()
        use com
#ifdef LAG
        use mpilag
#endif
        implicit none
        
        integer i,j,k,n1,ixp
	logical iex
        real rnx
c
	integer ithr
	integer balperiod
c
#ifdef OPENMP
        integer OMP_GET_NUM_THREADS,OMP_NUM_THREADS,
     1                OMP_GET_THREAD_NUM
#endif

#ifdef OPENMP
!$OMP PARALLEL private(ithr)
      num_thr = OMP_GET_NUM_THREADS()
      ithr = OMP_GET_THREAD_NUM()
      if (taskid.eq.0)   write (6,*) 'mpisetup: taskid,ithr,num_thr=',
     1               taskid,ithr,num_thr
!$OMP END PARALLEL
#else
      ithr=0
      num_thr = 1
#endif
c

      allocate( status(MPI_STATUS_SIZE,numtasks))

      dims(1) = 0
      dims(2) = 0

!    numtasks is devided into a iproc x jproc stencle
!

      call MPI_Dims_create(numtasks,2,dims,ierr)

#ifdef ONED 
       dims(2) = numtasks
       dims(1) = 1
#endif	

      if(dims(1) .gt. dims(2)) then
         dims(1) = dims(2)
         dims(2) = numtasks / dims(1)
      endif

	if (taskid.eq.0) then
	inquire(file='dims',exist=iex)
	if (iex) then
	  if (taskid.eq.0) print *, 'Reading grid from file dims'
	  open (999,file='dims')
	  read (999,*) dims(1), dims(2)
	  close (999)
	else
	  if (taskid.eq.0) print *, 'Creating grid with mpi_dims_create'
	endif  
	end if
	call MPI_BCAST (dims,2,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)


      iproc = dims(1)
      jproc = dims(2)

      if (iproc*jproc.ne.numtasks) then
!      if (taskid.eq.0) then
!        write (6,*)  'ABORT: invalid user-specified choice of iproc x jproc!'
!        write (6,*)  'Correct choices in the dims file'
!	  write (6,*) 'iproc,jproc,numtasks=',iproc,jproc,numtasks
!        call MPI_ABORT (MPI_COMM_WORLD,ierr)
!      end if
	jproc=numtasks/iproc
	dims(2)=jproc
      if (taskid.eq.0) then
        write (6,*)  'Warning: invalid user-specified choice of iproc x jproc'
        write (6,*)  'Proceed by adjusting jproc, now equal to',jproc
      end if
      end if

	if (jproc.gt.nz) then
	jproc=nz
	iproc=numtasks/jproc
      if (taskid.eq.0) then
        write (6,*)  'User-specified jproc was larger than nz'
        write (6,*)  'Proceed by adjusting iproc',iproc,jproc
      end if
	end if
	
#ifdef REVERSE_DIMS
       i = dims(1)  
       dims(1) = dims(2)
       dims(2) = i
#endif
          
	if (taskid.eq.0) write (6,*) 'mpisetup: iproc,jproc=',iproc,jproc
c
      periodic(1) = .false.
      periodic(2) = .false.
! creating cartesian processor grid
      call MPI_Cart_create(MPI_COMM_WORLD,2,dims,periodic,
     &     .false.,mpi_comm_cart,ierr)
	if (taskid.eq.0) write (6,*) 'mpisetup: after create'
! Obtaining process ids with in the cartesian grid
      call MPI_Cart_coords(mpi_comm_cart,taskid,2,cartid,ierr)
	if (taskid.eq.0) write (6,*) 'mpisetup: after coords'
! process with a linear id of 5 may have cartid of (3,1)
      ipid = cartid(1)
      jpid = cartid(2)
#ifdef REVERSE_DIMS
      ipid = cartid(2)
      jpid = cartid(1)
#endif
c	if (taskid.eq.0) write (6,*) 'mpisetup: a',taskid,ipid,jpid
c 	write (6,"('mpisetup: taskid,ipid,jpid=',4i3)") taskid,ipid,jpid

        allocate (ipid_all(0:numtasks-1))
        allocate (jpid_all(0:numtasks-1))
        call  MPI_ALLGATHER (ipid,1,MPI_INTEGER,ipid_all,1,MPI_INTEGER,
     1                       MPI_COMM_WORLD,mpierr)
        call  MPI_ALLGATHER (jpid,1,MPI_INTEGER,jpid_all,1,MPI_INTEGER,
     1                       MPI_COMM_WORLD,mpierr)

! here i is east-west j is north-south
! impid is west neighbour ippid is east neighbour and so on
      impid = ipid - 1
      ippid = ipid + 1
      jmpid = jpid - 1
      jppid = jpid + 1
!boundary processes
      if (ipid.eq.0) impid = MPI_PROC_NULL
      if (jpid.eq.0) jmpid = MPI_PROC_NULL
      if (ipid.eq.iproc-1) ippid = MPI_PROC_NULL
      if (jpid.eq.jproc-1) jppid = MPI_PROC_NULL
! using cart comworld create east-west(row) sub comworld
#ifdef REVERSE_DIMS
      remain_dims(1) = .true.
      remain_dims(2) = .false.
      call MPI_Cart_sub(mpi_comm_cart,remain_dims,mpi_comm_col,ierr)
	if (taskid.eq.0) write (6,*) 'mpisetup: b'
! using cart comworld create north-south(column) sub comworld
      remain_dims(1) = .false.
      remain_dims(2) = .true.
      call MPI_Cart_sub(mpi_comm_cart,remain_dims,mpi_comm_row,ierr)
	if (taskid.eq.0) write (6,*) 'mpisetup: c'
#else
      remain_dims(1) = .true.
      remain_dims(2) = .false.
      call MPI_Cart_sub(mpi_comm_cart,remain_dims,mpi_comm_row,ierr)
! using cart comworld create north-south(column) sub comworld
      remain_dims(1) = .false.
      remain_dims(2) = .true.
      call MPI_Cart_sub(mpi_comm_cart,remain_dims,mpi_comm_col,ierr)
#endif

	if (taskid.eq.0) write (6,*) 'mpisetup: A'

! mapping i onto iproc, i onto jproc, j onto iproc etc.
      allocate (iist(0:iproc-1))
      allocate (iisz(0:iproc-1))
      allocate (iien(0:iproc-1))
      allocate (jjst(0:jproc-1))
      allocate (jjsz(0:jproc-1))
      allocate (jjen(0:jproc-1))
      allocate (kist(0:iproc-1))
      allocate (kisz(0:iproc-1))
      allocate (kien(0:iproc-1))
      allocate (kjst(0:jproc-1))
      allocate (kjsz(0:jproc-1))
      allocate (kjen(0:jproc-1))
!
!Mapping 3-D data arrays onto 2-D process grid
! (nx+2,ny+2,nz) => (iproc,jproc)      
! 
      call MapDataToProc(nxh,iproc,iist,iien,iisz)
      call MapDataToProc(nypad,jproc,jjst,jjen,jjsz)
      call MapDataToProc(nzpad,iproc,kist,kien,kisz)
      call MapDataToProc(nz,jproc,kjst,kjen,kjsz)
!

      allocate(mymap(0:iproc-1))
      allocate(inverse_map(0:iproc-1))

#ifdef BAL_PERIOD
      bal_period = BAL_PERIOD
#else
      bal_period = 4
c     bal_period = 8
#endif

#ifdef LAG
! proposed by PKY, Jul 12, 2011
!       bal_period = iproc
! D. Buaria: removed Dec 2012 by making adjustments in bs_setup.f and spxyz_m.f

       bal_period = iproc
! Mar 2015, readded since it helps the new communication algorithm for particles
#endif


! Next 10 lines added by PKY, Dec 16, 2012
! allow user to override 'compiled' value of bal_period'
	balperiod=0
	if (taskid.eq.0) then
	inquire (file='balperiod.inp',exist=iex)
	if (iex) then
	open (1111,file='balperiod.inp')
	read (1111,*) balperiod
	end if
	end if
	call MPI_BCAST (balperiod,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	if (balperiod.ne.0) bal_period=balperiod
c
      if(taskid .le. 1) then
c         if(bal_period .lt. iproc) then
c            print *,'Error: period for load balancing is less than iproc, ',
c     &           bal_period,iproc
c         endif
         
         print *,'taskid,Balance Period is ',taskid,bal_period
      endif

	if(bal_period.eq.iproc) then
	if(taskid.eq.0) write(6,*) 'bal_period=iproc, might cause possible
     1                  conflict when restarting from bal_period=4'
	endif


	if (taskid.eq.0) write (6,*) 'mpisetup: B'

      call create_balance_map
	if (taskid.eq.0) write (6,*) 'mpisetup: C'

      ipid = mymap(ipid)

       xist = iist(ipid)
       yjst = jjst(jpid)
       zist = kist(ipid)
       zjst = kjst(jpid)
       xisz= iisz(ipid)
       yjsz= jjsz(jpid)
       zisz= kisz(ipid)
       zjsz= kjsz(jpid)
       xien = iien(ipid)
       yjen = jjen(jpid)
       zien = kien(ipid)
       zjen = kjen(jpid)

       

c         call spline_initialize !member function of module: mpilag
#ifdef LAG
         call bs_setup  !member function of module: mpilag
#endif

 	if (taskid.eq.0) write (6,*) ' exit mpisetup'
      end subroutine
!==================================================================       
      subroutine MapDataToProc (data,proc,st,en,szs)
!    
       implicit none 
       integer data,proc,st(0:proc-1),en(0:proc-1),szs(0:proc-1)
       integer i,size,nadd,size2
      size=data/proc
      nadd=mod(data,proc)
	size2=size
      if(nadd.ne.0) size2= size2+1
      st(0) = 1
      szs(0) = size2
      en(0) = size2
 	if (proc .gt. 1) then
      do i=1,proc-1
	   size2=size
 	   if (i.lt.nadd) size2=size2+1
         st(i) = st(i-1) + szs(i-1)
         szs(i) = size2
         en(i) = en(i-1) + size2
      enddo
      en(proc-1)= data 
      szs(proc-1)= data-st(i-1)+1
	endif
!
      end subroutine
!==================================================================       


!==================================================================       
      subroutine mpisetup2
!==================================================================       

        use comp
#ifdef LAG
        use mpilag
#endif
        implicit none
        
        integer ii,i,j,k,n1,ixp,x,c,n,xp,zp,z,xz(2)
        integer xz_tot
        real rnxi,rc
c
	integer ikmax

#ifdef HOMOGENEOUS

! Extent of local X dimension for a given task (the smaller of xisz 
! and cylinder radius)
	ikmax=nx*sqrt(2.)/3.
       max_al_x = min(xisz, ikmax +2-xist)
cccccccccx_al_x = min(xisz, (nx * sqrt(2.0))/3.0 +2-xist)
       if(max_al_x .le. 0) then
          max_al_x = 0
          xz_tot = 0
          num_al = 0
       else

! cut_z(x) is the number of positive Z wavenumbers for a given task 
! that fall outside the cylinder, for the array slice with a given x 
! (where x is the local X index)
       allocate(cut_z(max_al_x))
c       print *,taskid,': max_al_x=',max_al_x
       num_al = 0
       rnxi = 1.0/nx
       xz_tot = 0
       do ixp=1,max_al_x
          x = ixp + xist -1
! Compute how many positive Z wavevectors are outside the cylinder
          rc = nz*(1.0 - 2.0*sqrt(2.0/9.0 - (kx(x)*rnxi)**2))
          c = rc/2
          cut_z(ixp) = c

! For a given x, the number of Z points inside the cylinder = nz-1-c*2,
! arranged symmetrically around nzhp; that is nzhp point itself is always
! cut out, since it is outside the cylinder radius. The factor of 2 
! accounts for positive and negative wave numbers in Z. 
! xz_tot is the total number of points inside the cylinder. 
          xz_tot = xz_tot + nz - 1 - c*2
       enddo
       
       if(taskid .eq. 0) then
          print *,'xz_tot=',xz_tot
       endif

! Distribute the points inside the cylinder to jproc tasks
! (This is done only within columns of the proc. grid; rows are not affected)
! num_al is the number of XZ points on a given task
       num_al = xz_tot / jproc
       i = mod(xz_tot,jproc)
       if(i .ne. 0 .and. jpid .lt. i) then
          num_al = num_al+1
       endif

       allocate (mask(ny,num_al))

      endif

       if(jpid .eq. 0) then
          write(*,150) taskid,num_al,xz_tot
 150      format(i6,': Num_al,xz_tot =',i10,i10)
       endif

! Gather num_al from different tasks
       allocate(num_al_i(0:jproc-1))
       call mpi_allgather(num_al,1,MPI_INTEGER,num_al_i,1,MPI_INTEGER,
     &     mpi_comm_col,ierr)

! Subdivide the work in X and Z based on each processor's number of 
! points within the cylinder (num_al)
! Note: here x is local index, ranging from 1 to xisz, since X dimension 
! is decomposed among rows, while z is global index (from 1 to nz)
! The algorithm is as follows: start from x=z=1 for jpid=0 for each column, 
! then go one task at a time. Z index changes first. Add nz-1-cut_z(x)*2
! points for a given x, then increase x by 1 and repeat, until we reach 
! num_al total number of points allotted for the given task. Then we 
! stop, and continue from the nexxt z point with the next task.
       mystart_x = 1
       mystart_z = 1
       if(xz_tot .gt. 0) then
          do i=0,jpid-1
             call get_xz(num_al_i(i)+1,mystart_x,mystart_z,xz)
             mystart_x = xz(1)
             mystart_z = xz(2)
          enddo
       endif

c       print *,taskid,': mystart=',mystart_x,mystart_z

       if(mystart_x .gt. xisz) then
          print *,taskid,': Error: mystart_x =',mystart_x
       endif

       if(mystart_z .gt. nz) then
          print *,taskid,': Error: mystart_z =',mystart_z
       endif

#endif

      allocate (IfSndCnts(0:iproc-1))     
      allocate (IfRcvCnts(0:iproc-1))
      allocate (KfSndCnts(0:jproc-1))     
      allocate (KfRcvCnts(0:jproc-1))

      allocate (JrSndCnts(0:jproc-1))     
      allocate (JrRcvCnts(0:jproc-1))

      allocate (KrSndCnts(0:iproc-1))     
      allocate (KrRcvCnts(0:iproc-1))

      allocate (IfSndStrt(0:iproc-1))     
      allocate (IfRcvStrt(0:iproc-1))
      allocate (KfSndStrt(0:jproc-1))     
      allocate (KfRcvStrt(0:jproc-1))

      allocate (JrSndStrt(0:jproc-1))     
      allocate (JrRcvStrt(0:jproc-1))

      allocate (KrSndStrt(0:iproc-1))     
      allocate (KrRcvStrt(0:iproc-1))

#ifdef INPUT_PEN
      allocate (KfSndCnts_io(0:jproc-1))     
      allocate (KfRcvCnts_io(0:jproc-1))
      allocate (KfSndStrt_io(0:jproc-1))     
      allocate (KfRcvStrt_io(0:jproc-1))
#endif

#ifdef OUTPUT_PEN
      allocate (JrSndCnts_io(0:jproc-1))     
      allocate (JrRcvCnts_io(0:jproc-1))
      allocate (JrSndStrt_io(0:jproc-1))     
      allocate (JrRcvStrt_io(0:jproc-1))
#endif

#ifdef ALLTOALLV
        allocate (atac_start(0:jproc-1),atac_count(0:jproc-1))
        allocate (atar_start(0:iproc-1),atar_count(0:iproc-1))
#endif

!   start pointers and types of send  for the 1st forward transpose
      IfCntMax = iisz(0) * kisz(0) * b8 * 2 *yjsz
      if(mod(nxh,iproc) .ne. 0 .or. mod(nzpad,iproc) .ne. 0) then 
         IfCntUneven = .true.
      else
         IfCntUneven = .false.
      endif

      do i=0,iproc-1
         ii = mymap(i)
         IfSndStrt(i) = (iist(ii) -1)* zisz*b8*2*yjsz
         IfSndCnts(i) = iisz(ii) * zisz*b8*2*yjsz

!   start pointers and types of recv for the 1st forward transpose
         IfRcvStrt(i) = (kist(ii) -1) * xisz*b8*2*yjsz
         IfRcvCnts(i) = kisz(ii) * xisz*b8*2*yjsz
      end do

#ifdef HOMOGENEOUS
!   start pointers and types of send  for the 2nd forward transpose
      KfCntMax = num_al_i(0) *b8*2
      if(mod(xz_tot,jproc) .ne. 0) then
         KfCntUneven = .true.
      else
         KfCntUneven = .false.
      endif

      KfSndStrt(0) = 0
      KfSndCnts(0) = yjsz*b8*2 * num_al_i(0)
      KfRcvStrt(0) = 0
      KfRcvCnts(0) = num_al*jjsz(0)*b8*2
      do i=1,jproc-1
c         KfSndStrt(i) = (kjst(i) -1)*xisz*yjsz*b8*2
c         KfSndCnts(i) = xisz*yjsz*kjsz(i)*b8*2
         KfSndCnts(i) = yjsz*b8*2 * num_al_i(i)
         KfSndStrt(i) = KfSndStrt(i-1) + KfSndCnts(i-1)

!   start pointers and types of recv for the 2nd forward transpose
         KfRcvStrt(i) = (jjst(i) -1) * num_al*b8*2
         KfRcvCnts(i) = num_al*jjsz(i)*b8*2
      end do
#else
!   start pointers and types of send  for the 2nd forward transpose
      KfCntMax = xisz * jjsz(0) * kjsz(0) *b8*2
      if(mod(nypad,jproc) .ne. 0 .or. mod(nz,jproc) .ne. 0) then
         KfCntUneven = .true.
      else
         KfCntUneven = .false.
      endif

      do i=0,jproc-1
         KfSndStrt(i) = (kjst(i) -1)*xisz*yjsz*b8*2
         KfSndCnts(i) = xisz*yjsz*kjsz(i)*b8*2

!   start pointers and types of recv for the 2nd forward transpose
         KfRcvStrt(i) = (jjst(i) -1) * xisz * zjsz*b8*2
         KfRcvCnts(i) = xisz*zjsz*jjsz(i)*b8*2
      end do
#endif


!   start pointers and types of send  for the 1st inverse transpose
      JrCntMax = KfCntMax
      JrCntUneven = KfCntUneven
      do i=0,jproc-1
         JrSndStrt(i) = KfRcvStrt(i)
         JrSndCnts(i) = KfRcvCnts(i)

!   start pointers and types of recv for the 1st inverse transpose
          JrRcvStrt(i) = KfSndStrt(i)
          JrRcvCnts(i) = KfSndCnts(i)
      end do

!   start pointers and types of send  for the 2nd inverse transpose
      KrCntMax = iisz(0) * kisz(0) * b8 * 2
      if(mod(nxh,iproc) .ne. 0 .or. mod(nzpad,iproc) .ne. 0) then 
         KrCntUneven = .true.
      else
         KrCntUneven = .false.
      endif

      do i=0,iproc-1
         ii = mymap(i)
         KrSndStrt(i) = (kist(ii) -1) * xisz*b8*2
         KrSndCnts(i) = kisz(ii) * xisz*b8*2

!   start pointers and types of recv for the 2nd inverse transpose
         KrRcvStrt(i) = (iist(ii) -1) * zisz*b8*2
         KrRcvCnts(i) = zisz*iisz(ii)*b8*2
      enddo
!

#ifdef INPUT_PEN
      do i=0,jproc-1
         KfSndStrt_io(i) = (kjst(i) -1)*xisz*yjsz*b8*2
         KfSndCnts_io(i) = xisz*yjsz*kjsz(i)*b8*2

!   start pointers and types of recv for the 2nd forward transpose
         KfRcvStrt_io(i) = (jjst(i) -1) * xisz * zjsz*b8*2
         KfRcvCnts_io(i) = xisz*zjsz*jjsz(i)*b8*2
      enddo
#endif

#ifdef OUTPUT_PEN
      do i=0,jproc-1
         JrSndStrt_io(i) = KfRcvStrt_io(i)
         JrSndCnts_io(i) = KfRcvCnts_io(i)

!   start pointers and types of recv for the 1st inverse transpose
          JrRcvStrt_io(i) = KfSndStrt_io(i)
          JrRcvCnts_io(i) = KfSndCnts_io(i)
      end do
#endif

!      if(taskid .eq. 0) then
!         print *,taskid,': CntMax - If,Kf,Jr,Kr: ',IfCntMax,KfCntMax,JrCntMax,KrCntMax
!      endif

#ifdef CACHE_BL
        cache_l = CACHE_BL
#else
! A reasonable default
        cache_l = 32768
#endif

#ifdef NBL1_X
        NB1_X = NBL1_X
#else
! A reasonable default
        NB1_X = cache_l/(nzpad*b8)
#endif

#ifdef NBL1_Z
        NB1_Z = NBL1_Z
#else
! A reasonable default
        NB1_Z = cache_l/(xisz*b8)
#endif

#ifdef NBL2_Z
        NB2_Z = NBL2_Z
#else
! A reasonable default
        NB2_Z = cache_l/(nypad*b8)
#endif

#ifdef NBL2_Y
        NB2_Y = NBL2_Y
#else
! A reasonable default
        if(num_al .eq. 0) then
           NB2_Y = cache_l/(zjsz*b8)
        else
           NB2_Y = cache_l/(num_al*b8)
        endif
#endif


        if(NB1_X .gt. xisz) then
           NB1_X = xisz
        endif

        if(NB1_Z .gt. kisz(iproc-1)) then
           NB1_Z = kisz(iproc-1)
        endif
        
        if(NB2_Z .gt. num_al .and. num_al .gt. 0) then
           NB2_Z = num_al
        endif

        if(NB2_Y .gt. jjsz(jproc-1)) then
           NB2_Y = jjsz(jproc-1)
        endif

        !test that blocking factors are not zero -- if zero, set to 1
        if(NB1_X.eq.0)  NB1_X=1
        if(NB1_Z.eq.0)  NB1_Z=1 
        if(NB2_Z.eq.0)  NB2_Z=1
        if(NB2_Y.eq.0)  NB2_Y=1

c        if(taskid .eq. 0) then
c           print *,taskid,': Block factors NB1_x,NB1_z,NB2_z,NB2_y:',NB1_X,NB1_Z,NB2_Z,NB2_Y
c        endif

      end subroutine mpisetup2

! This subroutine takes starting counters, x and sz, and 'a' 
! which is the number of points by which to advance the counters. 
! It returns the xz(2) pair of the X and Z counters
! advanced by 'a' points. x is local (1 to xisz), z is global (1 to nz). 
!
! The algorithm is as follows: Z index changes first. Add nz-1-cut_z(x)*2
! points for a given x, then increase x by 1 and repeat, until we reach 
! a total number of points allotted for the given task. 

      subroutine get_xz(a1,x,start_z,xz)

      use comp
      implicit none
      integer a,z,a1
      integer x,z1,z2,c,start_z,xz(2)

      a = a1

 11   continue

      c = cut_z(x)
      z1 = nzhp - c - 1
      z2 = nzhp + c + 1

      z = a + start_z -1
      if(start_z .le. z1 .and. z .gt. z1) then
         z = z + c*2+1
      endif
            
      if(z .gt. nz) then
         if(x .eq. max_al_x .and. z .gt. nz+1) then
            write(6,100) taskid,a,start_z,x
 100        format(i5,': Error in get_xz:',i12,i5,i5)
         endif
         x = x+1
         start_z = 1
         a = z-nz
         goto 11

c         call get_xz(z - nz,x +1,1,xz)
      else
         xz(1) = x
         xz(2) = z
      endif

      return
      end subroutine

! Recursive version of get_xz - it did not work right, and is not used
      recursive subroutine get_xz2(l,x,zi,xz)

      use comp
      implicit none
      integer l,z
      integer x,z1,z2,c,zi,xz(2)

      c = cut_z(x)
      z1 = nzhp - c - 1
      z2 = nzhp + c + 1

      z = l + zi -1
      if(zi .le. z1 .and. z .gt. z1) then
         z = z + c*2+1
      else if(zi .lt. z2) then
         z = z + (z2-zi)
      endif
            
      if(z .gt. nz) then
         if(x .eq. max_al_x) then
            write(6,100) taskid,l,zi,x
 100        format(i5,': Error in get_xz: l,zi,x =max_al_x=',i10,i4,i4)
         endif
         call get_xz(z - nz,x +1,1,xz)
      else
         xz(1) = x
         xz(2) = z
      endif

      return
      end subroutine

! This is a subroutine to advance (x,z) counter by one. 
      subroutine next_xz(x,z)

      use comp
      implicit none
      integer x,z,c,z1,z2

      c = cut_z(x)
      z1 = nzhp - c -1
      z2 = nzhp + c +1

      if(z .lt. z1) then
         z = z+1
      else if(z .lt. z2) then
         z = z2
      else if(z .lt. nz) then
         z = z+1
      else
         z = 1
         x = x+1
      endif

      return
      end

      subroutine create_balance_map

      use comp
      implicit none
      integer i, id,st,en,n

      st = 0
      id = 0
      n = (iproc/bal_period)

 10   continue

      do i=st,iproc-1,bal_period
         mymap(id) = i
         inverse_map(i) = id
         id = id +1
         if(id .eq. iproc) goto 11
      enddo

      st = st + 1
      en = st + bal_period * n
      if(en .gt. iproc-1) then
         en = en - bal_period
      endif
      do i=en,st,-bal_period
         mymap(id) = i
         inverse_map(i) = id
         id = id +1
         if(id .eq. iproc) goto 11
      enddo

      st = st + 1

      goto 10

 11   continue
      return

      end

! This is a function mapping task id within a row to physical cores
! This is a load-balancing tool
! Currently we just do a round-robin with period 8 to have a good mix
! of loads for a given node with 8 cores (should work for 16 just as well)
c      function mymap(id)
c
c      use comp
c      integer myid,c,mymap
c
c      c = iproc / 4
c      mymap = id /c + mod(id,c) * 4
c
c      return
c      end

! Same as above, except this is the inverse map of cores to the row task id's
c      function inverse_map(id)
c
c      use comp
c      integer myid,c,inverse_map
c
c      c = iproc/4
c      inverse_map = mod((id * c),iproc) + id/4
c      
c      return
c      end

       
