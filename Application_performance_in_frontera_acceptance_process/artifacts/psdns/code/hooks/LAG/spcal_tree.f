      subroutine spcal_tree ( bfp,bsxy,f,nv,sf,n1,nrec,npdim,iod)

#ifdef LAG
#ifdef CAF
#ifdef CF_LAG
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
c n1:     change this
c nprec:     change this
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
c Early March 2015 (D. Buaria): allocate arrays outside
c inside of inside parallel region(s)
c
        use com, only:istep,nsteps
        use compart, only:nop,nop2,nom,iout
        use mpilag
        use lag_timers
        implicit none
c
	integer n1,npdim,nprec,nptp,nrec,n2,np,sbb
	integer bit,mypal,numdim,imod,id
	real d6,tt
c
      parameter (d6=1./6., tt=2./3.)
c
	integer iv,nv
c
      integer iod(3)
	integer ibx,iby,ibz

      real(p8) bsxy(bxisize,nby,bzjsize,nv)
	real(p8) bfx,bfy,bfz
      real(p8) f(npdim/numtasks,nv)
      real(p8) bfp(3,npdim/nrec)
c
	integer ib1,ib2,ib3,ix1,ix2,iz1,iz2,mp,mp1,mp2,m,ip,k,kk,
     1      i,ii,j,jj,itask,i1,i2,ip1,ip2
	real(b8) sf
	real(p8) pos,ym,yms
c

      real(p8), allocatable :: work(:,:)[:],swork1(:),temp(:)
	real(p8) swork2,const
	real(p8), allocatable :: bff(:,:)
c
        real(8) cpu_other,cpu_loop,cpu_comm,cpuspcal
        save cpu_other,cpu_loop,cpu_comm,cpuspcal
	save work,temp
c
      integer xp,zp
c
        integer ithr,omp_get_thread_num,inc,m1,m2,im
	integer, allocatable :: bp1(:),bpi(:)
        real*8 rtime1,rtime2,rtime3,rtime0
        real term1,term2
c
        integer, save :: icallsp
        data icallsp/0/

        icallsp=icallsp+1

        if(taskid.eq.0.and.icallsp.le.3) write(6,*) 'icallsp_tree=',icallsp

        rtime0=MPI_WTIME()
        rtime1=MPI_WTIME()
c
c for function or 1st derivative of spline, just pick the
c correct basis functions to use in the summation, based
c on the order of differentiation
c
	nprec = npdim/nrec
	np = nprec/numtasks

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
	if (n1.eq.1) then
c
        cpu_loop=0.
        cpu_comm=0.
        cpu_other=0.
        cpuspcal=0.
c
	sbb=-1
        if(sbb.lt.0) allocate (temp(nprec/numtasks))
        if(sbb.gt.0) allocate (temp(sbb))

        allocate (work(nprec,nv)[0:*],stat=ierr)
        if (ierr.gt.0) then
        write (6,*) 'spcal_tree: error at allocate work'
        end if
c
	end if

	allocate (bp1(0:num_thr-1),bpi(0:num_thr-1))
c
        ithr=0

        cpu_other=cpu_other+MPI_WTIME()-rtime1
        rtime2=MPI_WTIME()
c
	work(:,:)=0.
	call divide_thr (nprec, bp1, bpi, 'spcal_tree: nprec')

	allocate (bff(4,3),swork1(nv))

!$OMP PARALLEL private (ithr,ibx,iby,ibz,bfx,bfy,bfz,pos,ym,yms,ii,jj,kk,xp,zp,inc,m1,m2,im,term1,term2,bff,swork1,swork2,const) 
c
#ifdef OPENMP
        ithr = OMP_GET_THREAD_NUM()
#endif

	m1 = bp1(ithr)
	m2 = m1 + bpi(ithr) - 1
c
        do 15 m=m1,m2
c
	ibz=bfp(3,m)
 	if (ibz.lt.iz1-4.or.ibz.gt.iz2-1) go to 15
	ibx=bfp(1,m)
 	if (ibx.lt.ix1-4.or.ibx.gt.ix2-1) go to 15
c
	iby=bfp(2,m)
	bfy=bfp(2,m)-iby
	bfx=bfp(1,m)-ibx
	bfz=bfp(3,m)-ibz
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
 	do iv=1,nv
       term1= bsxy(xp,iby+1,zp,iv) * bff(1,ib2)
       term2= bsxy(xp,iby+2,zp,iv) * bff(2,ib2)
       term1= term1 + bsxy(xp,iby+3,zp,iv) * bff(3,ib2)
       term2= term2 + bsxy(xp,iby+4,zp,iv) * bff(4,ib2)
       swork2 = term1+term2
       swork1(iv)=swork1(iv)+swork2*const
 	end do
c
 10   continue
c
	do iv=1,nv
      work(m,iv)=work(m,iv)+swork1(iv)*bff(k,ib3)
	end do

 30   continue
c
 15   continue
c
!$OMP END PARALLEL

	deallocate (bff,swork1)

        cpu_loop=cpu_loop+MPI_WTIME()-rtime2

        sync memory
        call MPI_BARRIER (MPI_COMM_WORLD,mpierr)

        rtime3=MPI_WTIME()
c
#ifdef TEST
	do 200 iv=1,nv
c
 	call MPI_REDUCE (work(1,iv),fg,nprec,pptype,MPI_SUM,0,
     1                   MPI_COMM_WORLD,mpierr)
c
	call MPI_SCATTER (fg,nprec/numtasks,pptype,
     1                    f(n1,iv),nprec/numtasks,pptype,
     1                    0,MPI_COMM_WORLD,mpierr)
c
 200	continue
#endif


        numdim = log(numtasks + 0.1)/log(2.)

        do iv=1,nv

        bit=1

        do id=1,numdim

        mypal = ieor(taskid,bit)
        bit=bit*2
        imod = mod(taskid,bit)


        do itask=imod,numtasks-1,bit
        mp1=itask*np+1
        mp2=mp1+np-1
        temp(:) = work(mp1:mp2,iv)[mypal]
        work(mp1:mp2,iv) = work(mp1:mp2,iv) + temp(:)
        enddo

        sync all

        enddo ! id=1,numdim

        mp1=np*taskid+1
        mp2=mp1+np-1
        ip1=n1
        ip2=n1+np-1
!        f(:,iv) = work(mp1:mp2,iv)
        f(ip1:ip2,iv)=work(mp1:mp2,iv)

        enddo  ! iv=1,nv
c

        cpu_comm=cpu_comm+MPI_WTIME()-rtime3
c
	if(n1+mp/nrec-1.eq.mp) then
	do iv=1,nv      
	if (abs(sf-1.).gt.1.e-6) f(:,iv)=f(:,iv)*sf
	enddo
      deallocate (temp,work,stat=ierr)
	endif


	deallocate (bp1,bpi)

         cpuspcal=cpuspcal+MPI_WTIME()-rtime0


	if(n1+mp/nrec-1.eq.mp) then

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
	if (npdim.eq.nop) then
	ii=1
	else if (npdim.eq.nop2) then
	ii=2
	else if (npdim.eq.nom) then
	ii=3
	end if
c
	cpu_spcal(1,ii,jj)=cpu_spcal(1,ii,jj)+cpu_loop
	cpu_spcal(2,ii,jj)=cpu_spcal(2,ii,jj)+cpu_comm
	cpu_spcal(3,ii,jj)=cpu_spcal(3,ii,jj)+cpu_other
	cpu_spcal(4,ii,jj)=cpu_spcal(4,ii,jj)+cpuspcal
c
	endif
	
 99   return

#endif
#endif
#endif
      end

