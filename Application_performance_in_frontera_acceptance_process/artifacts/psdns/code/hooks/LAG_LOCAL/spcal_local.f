      subroutine spcal_local ( nxyz,shif,bsxy,posn,func,npar,nv,sf,npdim,iod, ityp)

#ifdef LAG

! Mar 20, 2015 D. Buaria
! new version to directly interpolate the required property 
! from particle position

! INPUT
! nxyz  = grid points in each direction
! shif  = grid shifting from pseudo-spectral scheme
! bsxy  = spline coefficients 
! posn  = particle positions 
! nv    = no. of variables (max allowed 3)
! npdim = total no. of particles 
! npar  = no. of particles locally on MPI task  (npar .ne. npdim/numtasks)
! sf    = 1.
! iod   = index related to vel. and vel. gradient
! ityp  = parameter for storing timing data


! OUTPUT
! func = property to be interpolated: velocity or vel. gradient


! OTHERS
! swork1 & swork2:  working variables in the summation procedure
! iod(n) (n=1,2,3) is 0 (function) or 1 (derivative) in direction n
c
c
        use com, only:istep,nsteps
        use compart, only:nop,nop2,nom,ngm,iout,nppad
        use mpilag
        use lag_timers
        implicit none
c
	integer n1,npdim,nprec,nptp,nrec,n2
	real d6,tt,sum1
c
	integer itid,jtid,tid,iaid,jaid
	integer xnew,znew,xpp,zpp
	integer ipx,jpx,ipmax,jloc
	real(p8) px,temp,px1,px2,px3,px4
c
      integer npar,nxyz(3)
      real(b8) shif(3),gfn(3)
	integer iv,nv
        integer ityp, ncomm, ncomm0
c
      integer iod(3)
	integer ibx,iby,ibz

#if defined(CAF) && defined(CF_LAG)
      real(p8) bsxy(xiind,nby,zjind,nv)[0:*]
#else
      real(p8) bsxy(xiind,nby,zjind,nv)
#endif

      real(p8) posn(nppad*npdim/numtasks,3)
      real(p8) func(nppad*npdim/numtasks,nv)
      real(p8), allocatable :: btemp(:,:),bt2(:,:)
      real(p8), allocatable :: bfp(:)
       real(p8) bfx,bfy,bfz
c
	integer ib1,ib2,ib3,ix1,ix2,iz1,iz2,mp,mp1,mp2,m,ip,k,kk,
     1      i,ii,j,jj,itask
	real(b8) sf
	real(p8) pos,ym,yms
c

      real(p8), allocatable :: swork1(:)
	real(p8) swork2,const
	real(p8), allocatable :: bff(:,:)
c
        real(8) cpu_other,cpu_loop,cpu_comm,cpuspcal
        save cpu_other,cpu_loop,cpu_comm,cpuspcal
c
      integer xp,zp
	integer, save :: icls = 0
c
        integer ithr,omp_get_thread_num,inc,m1,m2,im
	integer, allocatable :: bp1(:),bpi(:)
        real*8 rtime1,rtime2,rtime3,rtime0
        real term1,term2
	integer disp_int, win
	integer (kind=MPI_ADDRESS_KIND) lowerbound,bssize,realextent,disp_aint
c
	ncomm=0
	icls=icls+1
	if(taskid.eq.0) write(6,*) 'spcal_local, npdim=',npdim

	if(npar.eq.0) return

        rtime0=MPI_WTIME()
        rtime1=MPI_WTIME()
c
#if defined(CAF) && defined(CF_LAG)
#else
	call MPI_TYPE_GET_EXTENT (pptype, lowerbound, realextent, ierr)
	disp_int = realextent
	bssize = xiind*nby*zjind*nv*realextent
	call MPI_WIN_CREATE(bsxy, bssize, disp_int, MPI_INFO_NULL,
     1                        MPI_COMM_WORLD, win, ierr)

#endif



c for function or 1st derivative of spline, just pick the
c correct basis functions to use in the summation, based
c on the order of differentiation
c
	if (nv.gt.3) then
	write(6,*) 'nv should be always <= 3, nv=',nv
	stop
	endif

	do j=1,3
	gfn(j) = shif(j) + 1.
	enddo

	d6 = 1./6.
	tt = 2./3.

      ib1=1+iod(1)*3
      ib2=2+iod(2)*3
      ib3=3+iod(3)*3
c
	ix1=bxistart
	ix2=bxiend
	iz1=bzjstart
	iz2=bzjend
c
c
        cpu_loop=0.
        cpu_comm=0.
        cpu_other=0.
        cpuspcal=0.
c
	func=0.

c
        ithr=0

        cpu_other=cpu_other+MPI_WTIME()-rtime1
        rtime2=MPI_WTIME()
c
	allocate (bff(4,3),swork1(nv),bfp(3),btemp(nv,4))
	allocate (bt2(nv,4))
!	call divide_thr (nprec, bp1, bpi, 'spcal_ars: nprec')


!$OMP PARALLEL private (ithr,ibx,iby,ibz,bfx,bfy,bfz,pos,
!$OMP& ym,yms,ii,jj,kk,xp,zp,inc,m1,m2,im,term1,term2,bff,
!$OMP& swork1,swork2,const,px,px1,px2,ipx,jpx,jloc,bfp,btemp,temp,
!$OMP& jtid,itid,tid,znew,zpp,xnew,xpp,jaid,iaid,px3,px4,bt2) 
c
#ifdef OPENMP
        ithr = OMP_GET_THREAD_NUM()
#endif


	btemp=0.
	bt2=0.
	bff=0.
	bfp=0.

c
!$OMP DO
        do 15 m=1,npar


	do 20 j=1,3

        px1=posn(m,j)- gfn(j) 
	px = px1 - real(nxyz(j))
        ipx=floor(px)
        jpx=ipx
        if (ipx.lt.0) jpx=ipx+1000*nxyz(j)
        jloc=mod(jpx,nxyz(j))
        bfp(j)=real(jloc)+px-real(ipx)


!	if(icls.eq.1) then
!	write(7000+taskid,*) m,j,posn(m,j),nxyz(j),shif(j),gfn(j),px1,px2,px3,px4,px,bfp(j)
!	endif

        if (jloc.lt.0.or.jloc.gt.nxyz(j)-1) then
        write (6,"('spcal_local, bad jloc=',i7,i5,i3,1p,2e12.4,2i9,1p,e12.4,i12)")
     1           taskid,m,j,posn(m,j),px,ipx,jloc,bfp(j),npdim
	stop
        end if
c
 20   continue

c
	ibx=bfp(1)
	iby=bfp(2)
	ibz=bfp(3)

	bfx=bfp(1)-ibx
	bfy=bfp(2)-iby
	bfz=bfp(3)-ibz


!	if (ibz.lt.iz1-4.or.ibz.gt.iz2-1) go to 15
!	if (ibx.lt.ix1-4.or.ibx.gt.ix2-1) go to 15
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
c
      do 30 k=1,4
c
	kk = ibz+k
	if (kk.gt.nbz) kk=kk-nbz

! obtain jpid corresponding to kk and then 
! corresponding z index (znew) of bsxy array on jpid
        if(kk.gt.3*zjind) then
        jtid = (kk - 4)/(zjind-1)
        else
        jtid = (kk-1)/zjind
        endif
        znew = kk - bzjst(jtid) + 1

#ifdef DEBUG
	if(znew.lt.1.or.znew.gt.bzjen(jtid)) then
	write(6,*) 'znew out of bounds, taskid,jtid=',taskid,jtid,kk,znew
	stop
	endif
	if (jtid.lt.0.or.jtid.ge.jproc) then
	write(6,*) 'jtid out of bounds, taskid,jtid=',taskid,jtid,kk,znew
	stop
	endif
!      if (kk.lt.iz1.or.kk.gt.iz2) go to 30
!      zp=kk-iz1+1
#endif

	swork1(:)=0.

      do 10 i=1,4
c
	ii=ibx+i
	if (ii.gt.nbx) ii=ii-nbx

! obtain ipid corresponding to ii and then 
! corresponding x index (xnew) of bsxy array on jpid
        if(ii.gt.3*xiind) then
        itid = (ii - 4)/(xiind-1)
        else
        itid = (ii-1)/xiind
        endif
        xnew = ii - bxist(itid) + 1

! obtain taskid corresponding to jpid and ipid obtained above
        tid = jtid*iproc + itid

#ifdef NCOMM
	if (tid.ne.taskid) ncomm = ncomm+1
#endif

#ifdef DEBUG
	if(xnew.lt.1.or.xnew.gt.bxien(itid)) then
	write(6,*) 'xnew out of bounds, taskid,itid=',taskid,itid,ii,xnew
	stop
	endif
	if (itid.lt.0.or.itid.ge.iproc) then
	write(6,*) 'itid out of bounds, taskid,itid=',taskid,itid,ii,xnew
	stop
	endif
!      if (ii.lt.ix1.or.ii.gt.ix2) go to 10
!      xp=ii-ix1+1

	if (tid.lt.0.or.tid.ge.numtasks) then
	write(6,*) 'tid out of bounds, taskid,tid=',taskid,tid
	stop
	endif
#endif

! use CAF to obtain the spline coefficients
! if tid = taskid then coeff obtained from local memory

        rtime3=MPI_WTIME()

#if defined(CAF) && defined(CF_LAG)
        btemp(1:nv,1) = bsxy(xnew,iby+1,znew,1:nv)[tid]
        btemp(1:nv,2) = bsxy(xnew,iby+2,znew,1:nv)[tid]
        btemp(1:nv,3) = bsxy(xnew,iby+3,znew,1:nv)[tid]
        btemp(1:nv,4) = bsxy(xnew,iby+4,znew,1:nv)[tid]
#else



	do jj=1,4
	bt2(1:nv,jj) = bsxy(xnew,iby+jj,znew,1:nv)
	enddo

	bssize = xiind*nby*zjind
	call MPI_WIN_FENCE (0, win, ierr)


	if (tid.ne.taskid) then

	do jj=1,4
	do iv=1,nv
	disp_aint = xnew + (iby+jj-1)*xiind + (znew-1)*xiind*nby + (iv-1)*bssize
	call MPI_GET (bsxy(xnew,iby+jj,znew,iv), 1, pptype, tid, 
     1                  disp_aint-1, 1, pptype, win, ierr)
	enddo
	enddo

	endif

	call MPI_WIN_FENCE (0, win, ierr)

	do jj=1,4
	btemp(1:nv,jj) = bsxy(xnew,iby+jj,znew,1:nv)
	enddo
	if (tid.ne.taskid) then
	do jj=1,4
	bsxy(xnew,iby+jj,znew,1:nv) = bt2(1:nv,jj)
	enddo
	endif
#endif


        cpu_comm=cpu_comm+MPI_WTIME()-rtime3
c
	const=bff(i,ib1)
 	do iv=1,nv
       term1= btemp(iv,1)*bff(1,ib2) + btemp(iv,3)*bff(3,ib2)
       term2= btemp(iv,2)*bff(2,ib2) + btemp(iv,4)*bff(4,ib2) 
       swork2 = term1+term2
       swork1(iv)=swork1(iv)+swork2*const
 	end do
c
 10   continue
c
	const = bff(k,ib3)

	do iv=1,nv
      func(m,iv)=func(m,iv)+swork1(iv)*const
	end do

 30   continue
c
 15   continue
!$OMP END DO
c

!$OMP END PARALLEL


	do iv=1,nv      
	if (abs(sf-1.).gt.1.e-6) func(:,iv)=func(:,iv)*sf
	enddo

        cpu_loop=cpu_loop+MPI_WTIME()-rtime2
! Mar 25, 2015, D. Buaria: added for new algorithm because of how
! the loop/comm times are computed and reported
	cpu_loop = cpu_loop - cpu_comm
c


	deallocate (bff,swork1,bfp,btemp)


#if defined(CAF) && defined(CF_LAG)
	sync memory
	call MPI_BARRIER (MPI_COMM_WORLD,mpierr)
#else
	call MPI_WIN_FREE (win, ierr)
#endif

#ifdef NCOMM
! gives a measure of how much communication is done, however 
! affects the performance somewhat also. 
! use for benchmarking purposes but disabled for production
	call MPI_ALLREDUCE (ncomm,ncomm0,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
	if (taskid.eq.0) write(6,*) 'spcal, ncomm0=',ncomm0/numtasks
#endif


         cpuspcal=cpuspcal+MPI_WTIME()-rtime0


	if (iout.eq.0) then
	jj=1
	if (npdim.eq.nop) then 
	nspcal=nspcal+1
!	if (taskid.eq.0) write (6,*) 'istep,nspcal=',istep,nspcal,nv
	end if
	else
	jj=2
	if (npdim.eq.nop) then 
	nspcal_io=nspcal_io+1
!	if (taskid.eq.0) write (6,*) 'istep,nspcal_io=',istep,nspcal_io,nv
	end if
	end if
c
        if (ityp.ge.1.and.ityp.le.3) then
	cpu_spcal(1,ityp,jj)=cpu_spcal(1,ityp,jj)+cpu_loop
	cpu_spcal(2,ityp,jj)=cpu_spcal(2,ityp,jj)+cpu_comm
	cpu_spcal(3,ityp,jj)=cpu_spcal(3,ityp,jj)+cpu_other
	cpu_spcal(4,ityp,jj)=cpu_spcal(4,ityp,jj)+cpuspcal
        end if
c
 99   return


#endif
      end

