      subroutine spcal_pipe ( bsxy,ppi,iv1,ivv,nv,sf,npdim,iod, ityp)
c      subroutine spcal_ars ( bfp,bsxy,f,nv,sf,nrec,npdim,iod, ityp)

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
c npdim:  declared upper limit of nop
c
c swork1 & swork2:  working variables in the summation procedure
c
c iod(n) (n=1,2,3) is 0 (function) or 1 (derivative) in direction n
c
c generalised to include use of spline derivatives in ns code
c
c Fix by D. Buaria 9/1/2012: dimensioning of 'work' array
c is changed to work with cases where nprec is not 
c integer multiple of num_thr
c
        use mpilag
        use com, only:istep,nsteps
        use compart, only:nop,nop2,nom,ngm,iout,gsh
        use lag_timers
        implicit none
c
	integer npdim
	real d6,tt
c
      parameter (d6=1./6., tt=2./3.)
c
	integer iv,nv,iv1,iv2,ivv
        integer ityp
c
      integer iod(3)
	integer ibx,iby,ibz

	real(p8) ppi(nv,npdim/numtasks)
      real(p8) bsxy(bxisize,nby,bzjsize,ivv)
	real(p8) bfx,bfy,bfz
c
	integer ib1,ib2,ib3,ix1,ix2,iz1,iz2,mp,mp1,mp2,m,ip,k,kk,
     1      i,ii,j,jj,itask,minc
	real(b8) sf(ivv)
	real(p8) pos,ym,yms
c
        integer mpistatus(MPI_STATUS_SIZE)
        integer isource,idest,sendtag,recvtag,ipipe

	real(p8) swork2,const,tmp
	real(p8), allocatable :: bff(:,:),swork1(:)
c
        real(8) cpu_other,cpu_loop,cpu_comm,cpuspcal
        save cpu_other,cpu_loop,cpu_comm,cpuspcal
c
      integer xp,zp
c
        integer ithr,omp_get_thread_num,inc,m1,m2,im
        real*8 rtime1,rtime2,rtime3,rtime0
        real term1,term2
c
	real shif(3),gfn(3),fbox(3)

	if (taskid.eq.0) write (6,*) 'enter spcal_pipe'
c
        rtime0=MPI_WTIME()
        rtime1=MPI_WTIME()
c
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
	mp=npdim/numtasks
c
        shif(:)=gsh(:,1)
        gfn(:)=1.+shif(:)
        fbox(1)=float(nx)
        fbox(2)=float(ny)
        fbox(3)=float(nz)
c
        cpu_loop=0.
        cpu_comm=0.
        cpu_other=0.
        cpuspcal=0.
c
	iv2=iv1+ivv-1
c
        idest=taskid+1
        idest=mod(idest+numtasks,numtasks)
        isource=taskid-1
        isource=mod(isource+numtasks,numtasks)
c
        ithr=0
        minc=npdim/numtasks/num_thr

        cpu_other=cpu_other+MPI_WTIME()-rtime1
c
	ppi(iv1:iv2,:)=0.

	do 100 ipipe=1,numtasks
c
        rtime2=MPI_WTIME()
c
!$OMP PARALLEL private (ithr,ibx,iby,ibz,bfx,bfy,bfz,pos,ym,yms,ii,jj,kk,xp,zp,inc,m1,m2,im,term1,term2,bff,swork1,swork2,const) 
c
#ifdef OPENMP
        ithr = OMP_GET_THREAD_NUM()
#endif

	allocate (bff(4,3),swork1(iv1:iv2))

        m1=ithr*minc+1
        m2=m1+minc-1
c
      do 15 m=m1,m2
c

        tmp=mod(ppi(3,m)-gfn(3),fbox(3))
        bfz=tmp-floor(tmp/fbox(3))*fbox(3)
        ibz=bfz
        if (ibz.lt.iz1-4.or.ibz.gt.iz2-1) go to 15

        tmp=mod(ppi(1,m)-gfn(1),fbox(1))
        bfx=tmp-floor(tmp/fbox(1))*fbox(1)
        ibx=bfx
        if (ibx.lt.ix1-4.or.ibx.gt.ix2-1) go to 15
c
        tmp=mod(ppi(2,m)-gfn(2),fbox(2))
        bfy=tmp-floor(tmp/fbox(2))*fbox(2)
        iby=bfy

        bfy=bfy-iby
        bfx=bfx-ibx
        bfz=bfz-ibz

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
      kk=ibz+k
      if (kk.lt.iz1.or.kk.gt.iz2) go to 30
      zp=kk-iz1+1
      swork1(:)=0.
c
      do 10 i=1,4
c
      ii=ibx+i
      if (ii.lt.ix1.or.ii.gt.ix2) go to 10
      xp=ii-ix1+1
c
	const=bff(i,ib1)
 	do iv=1,ivv
       term1= bsxy(xp,iby+1,zp,iv) * bff(1,ib2)
       term2= bsxy(xp,iby+2,zp,iv) * bff(2,ib2)
       term1= term1 + bsxy(xp,iby+3,zp,iv) * bff(3,ib2)
       term2= term2 + bsxy(xp,iby+4,zp,iv) * bff(4,ib2)
       swork2 = term1+term2
       swork1(iv+iv1-1)=swork1(iv+iv1-1)+swork2*const
 	end do
c
 10   continue
c
        const=bff(k,ib3)
        do iv=iv1,iv2
        ppi(iv,m)=ppi(iv,m)+swork1(iv)*const
        end do

 30   continue
c
 15   continue
c
	deallocate (bff,swork1)

!$OMP END PARALLEL

        cpu_loop=cpu_loop+MPI_WTIME()-rtime2
c
	sendtag=0
	recvtag=0
        rtime3=MPI_WTIME()
c
        call mpi_sendrecv_replace (ppi,nv*mp,pptype,idest,sendtag,
     1         isource,recvtag,MPI_COMM_WORLD,mpistatus,ierr)
        cpu_comm=cpu_comm+MPI_WTIME()-rtime3
c
 100	continue

	do iv=iv1,iv2
	ppi(iv,:)=ppi(iv,:)*sf(iv-iv1+1)
	end do



         cpuspcal=cpuspcal+MPI_WTIME()-rtime0


	if (iout.eq.0) then
	jj=1
	if (npdim.eq.nop) then 
	nspcal=nspcal+1
c	if (taskid.eq.0) write (6,*) 'istep,nspcal=',istep,nspcal,nv
	end if
	else
	jj=2
	if (npdim.eq.nop) then 
	nspcal_io=nspcal_io+1
c	if (taskid.eq.0) write (6,*) 'istep,nspcal_io=',istep,nspcal_io,nv
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

	if (taskid.eq.0) write (6,*) ' exit spcal_pipe'

#endif
      end

