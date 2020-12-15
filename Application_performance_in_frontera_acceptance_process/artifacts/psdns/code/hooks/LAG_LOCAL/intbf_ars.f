      subroutine intbf_ars ( xyzl,nxyz,shif, xp,n1,nrec,npdim, bfp,slant,ityp )

#ifdef LAG


! MPI_IN_PLACE introduced in ALLGATHER call by D. Buaria, 8/23/12
! also extended to give correct results for 6 threads

c Feb 2011: "ars" stands for allgather/allgatherv (this routine)
c and reduction, scatter (in spcal_ars).
c Feb 6, 2009: this version switches the subscripts in bfp,
c and uses one large ALLGATHER covering all three coordinate
c components in one call.
c
c "Batches operation": assume each 1 task carry
c an equal-sized subset of the global population of particles
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
c      n1        - change this
c      npdim     - first dimension of xp.
c      ipbc      - 1 if particles are inside the box
c                - 2 if particles may be outside, and 
c                    periodicity to be accounted for
c
c  output:
c      bfp(i,j) the 1st of 4 basis functions in the jth
c                   direction for the ith particle
c
        use compart, only: nop,nop2,nom,ngm,iout
        use com, only:istep,nsteps
        use mpicom
        use mpilag, only: p8, pptype
        use lag_timers
c
	implicit none
c
	integer npdim,i,j,k,j1,mp,mp1,mp2,ipbc,m,jloc,n1,nrec,ii,jj
c
	real(b8)  d6,tt,slant,toler
c
	real(p8) px,xyzl8
c
      integer nxyz(3)
      real(b8) xyzl(3)
	real(b8) shif(3),gfn(3)
      real(p8) xp(npdim/numtasks,3)
      real(p8) bfp(3,npdim/nrec)
c
	integer ip,ipx,jpx,ipmax
!	real(p8), allocatable :: bfpm(:,:)

        integer ityp
	integer nbad,no_bad
	real(8) rtime0,rtime1,rtime2,rtime3
        real(8) cpu_other,cpu_loop,cpu_comm,cpuintbf
        save cpu_other,cpu_loop,cpu_comm,cpuintbf

c
c note that items in next line must be declared double precision:
c otherwise some of the bfps may not be correct
c
      real(b8) denom,ploc

	integer, save :: icls = 0
        integer ithr,omp_get_thread_num
        integer, allocatable :: nbad_thr(:)
        integer, allocatable :: bp1(:),bpi(:)
c
        integer, allocatable :: num_all(:),idisp(:)
        integer itask
	real*8 workj(3)
	real(p8) sizej(3)
c
	d6=1./6.
	tt=2./3.

	if(npdim.eq.16*1024*1024) icls=icls+1

	rtime0=MPI_WTIME()
	rtime1=MPI_WTIME()
c
        if (n1.eq.1) then
        cpu_loop=0.
        cpu_comm=0.
        cpu_other=0.
        cpuintbf=0.
	end if
c
 	if (npdim.gt.0) call ppminmax ('enter intbf_ars',xp,npdim)

c gfn(i) is coord of first node in i-th dirn
c
c EXPRESS PARTICLE COORDINATES IN TERMS OF GRID SPACINGS
c (even if on non-2pi^3 domains)
c
      do k=1,3
      gfn(k)=shif(k)+1.d0
      end do
c
c
        allocate (nbad_thr(0:num_thr-1))
        allocate (bp1(0:num_thr-1))
        allocate (bpi(0:num_thr-1))
c
	cpu_other = cpu_other + MPI_WTIME() -rtime1
	rtime2=MPI_WTIME()
c
	nbad=0
c
c loop over particles
c
	ithr = 0

	mp=npdim/nrec/numtasks
	mp1=taskid*mp+1
	mp2=mp1+mp-1
!	allocate (bfpm(3,mp1:mp2))
c
	call divide_thr (mp, bp1, bpi, 'intbf_ars: mp')
c

!$OMP PARALLEL private (ithr, i,j,m,ip,px,ipx,jpx,ploc,jloc,mp1,mp2,ipmax) shared (nbad_thr)
!$OMP& REDUCTION (+:nbad)
#ifdef OPENMP
        ithr = OMP_GET_THREAD_NUM()
#endif
c
        nbad_thr(ithr)=0
c
	mp1 = taskid*mp + bp1(ithr)
	mp2 = mp1 + bpi(ithr) - 1

	ipmax=0

      do 25 i=mp1,mp2
c
	m=i-taskid*mp
c
c locate the particle and
c calculate local normalised coord (between 0 and 1)
c
c	ip=m+n1
	ip=m+n1-1
	ipmax=max0(ip,ipmax)
c
      do 20 j=1,3
c
        px=xp(ip,j)-nxyz(j)-gfn(j)
ccc        ipx=px
c fix by PKY on 6/28/2014
        ipx=floor(px)
        jpx=ipx
c        if (ipx.lt.0) jpx=ipx+100*nxyz(j)
        if (ipx.lt.0) jpx=ipx+1000*nxyz(j)
        jloc=mod(jpx,nxyz(j))
!        bfpm(j,i)=jloc+px-ipx
        bfp(j,i)=jloc+px-ipx

        if (jloc.lt.0.or.jloc.gt.nxyz(j)-1) then
        write (6,"('intbf_ars, bad jloc=',i5,4i9,i3,1p,2e12.4,2i5,1p,e12.4,i9)")
     1           taskid,n1,i,m,ip,j,xp(ip,j),px,ipx,jloc,bfp(j,i),npdim
        nbad_thr(ithr)=nbad_thr(ithr)+1
        end if
c
 20   continue

c
 25   continue
c
c
        nbad=nbad+nbad_thr(ithr)
c
!$OMP END PARALLEL


        cpu_loop = cpu_loop + MPI_WTIME() - rtime2
	rtime3=MPI_WTIME()

	call MPI_ALLREDUCE (nbad,no_bad,1,MPI_INTEGER,MPI_SUM,
     1                      MPI_COMM_WORLD,mpierr)
	if (no_bad.gt.0) then
	call MPI_FINALIZE (mpierr)
	stop 'Stop in intbf_ars: bad indices encountered'
	end if

#ifdef ALLGATHERV
c
        allocate (num_all(numtasks),idisp(numtasks),stat=ierr)
        num_all(:)=mp
        do itask=1,numtasks
        idisp(itask)=(itask-1)*mp
        end do
        do j=1,3
        call MPI_ALLGATHERV (bfpm(mp1,j),mp,pptype,
     1                       bfpm(1,j),num_all,idisp,pptype,
     1                       MPI_COMM_WORLD,mpierr)
        end do
        deallocate (num_all,idisp,stat=ierr)
c
#else

	call MPI_ALLGATHER (MPI_IN_PLACE,mp*3,pptype,
     1                      bfp(1,1),mp*3,pptype,
     1                      MPI_COMM_WORLD,mpierr)
c
c
#endif
c
        cpu_comm = cpu_comm + MPI_WTIME() - rtime3

         cpuintbf=cpuintbf+MPI_WTIME()-rtime0


	mp=npdim/nrec/numtasks
        if(n1+mp-1.eq.npdim/numtasks) then
c
        if (iout.eq.0) then
        jj=1
        if (npdim.eq.nop) then
        nintbf=nintbf+1
        end if
        else
        jj=2
        if (npdim.eq.nop) then
        nintbf_io=nintbf_io+1
        end if
        end if
c
        if (ityp.ge.1.and.ityp.le.3) then
        cpu_intbf(1,ityp,jj)=cpu_intbf(1,ityp,jj)+cpu_loop
        cpu_intbf(2,ityp,jj)=cpu_intbf(2,ityp,jj)+cpu_comm
        cpu_intbf(3,ityp,jj)=cpu_intbf(3,ityp,jj)+cpu_other
        cpu_intbf(4,ityp,jj)=cpu_intbf(4,ityp,jj)+cpuintbf
        end if
c
        endif
c
!	deallocate (bfpm)

	deallocate (nbad_thr,bp1,bpi)

	call MPI_BARRIER (MPI_COMM_WORLD,mpierr)

#endif
c
      return
      end
