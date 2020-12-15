       subroutine partvg (indexp)
c
c  routine to perform cubic spline interpolation for 
c  velocity gradient fluctuations following fluid particles
c (by fitting a cubic spline for 8 of the gradients)
c
c  called by sub. partsp
c
#ifdef LAG
#ifdef LGRAD
c
	use compart
	implicit none
	include 'intvars'
c
	real slant
	integer indx, ir, nprec
	integer i,j,jj,indexp
	real(b8) sum,termu,termv,termw,term,sum0
c
	real(8) rtime0,rtime1
	real(b8) cpu(6)
	real(b8) cpu_overall

c
c
c in the current formulation, since each scalar is processed
c sequentially, all the scalars can share one storage position
c in the (bsxy) spline array
c
	if (lpgrad(1).eq.0) return
c
	if (taskid.eq.0) write (6,*) 'enter partvg, istep=',istep
c
	rtime0=MPI_WTIME()
	cpu=0.

	allocate (udxzr(nx,zisz,yjsz,3))

c
	indexp=7
c

      do 10 i=1,3
c
      if (lpgrad(i).ne.1) go to 10
c
	rtime1=MPI_WTIME()
 	call xktran (u(1,1,1,i),u(1,1,1,i),u(1,1,1,i),1)
	cpu(1)=cpu(1)+MPI_WTIME()-rtime1
c
      ipj1(i+nc)=indexp
c
      if (i.lt.3) then
      ipj2(i+nc)=indexp+3
      else
      ipj2(i+nc)=indexp+2
      end if
c
	jj=3
	if (i.eq.3) jj=2

c dw/dz can be skipped, since only 8 of the velocity gradients
c are independent in incompressible flow
c
	rtime1=MPI_WTIME()
   	call udxphy (u(1,1,1,i),udxzr,udxzr,jj)
	cpu(2)=cpu(2)+MPI_WTIME()-rtime1
c
c form basis spline coefficients
c
c  	call spxyz_m (udxzr,bs,jj)
c
	rtime1=MPI_WTIME()
	do j=1,jj
	call spxyz_m (udxzr(1,1,1,j),bs(1,1,1,j),1)
	end do
	cpu(3)=cpu(3)+MPI_WTIME()-rtime1

	slant=-b12(2)*ps(2)/ps(1)*nx/ny
c
c interpolate 
c
	rtime1=MPI_WTIME()
c
	pp(:,indexp)=0.
c
        nprec = ndpart/nrec1
        allocate (bfpos(3,nprec))

        do ir=1,nrec1
        indx = (ir-1)*nprec/numtasks + 1

	call intbf_ars ( xyzl,nxyz,gsh(1,1),pp(1,1),
     1              indx,nrec1,ndpart,bfpos,slant, 1 )

#ifdef CF_LAG
        call spcal_tree  ( bfpos, bs(1,1,1,1), pp(1,indexp), jj,
     1             1., indx, nrec1, ndpart, iod(1,1) )
#else
        call spcal_ars  ( bfpos, bs(1,1,1,1), pp(1,indexp), jj,
     1             1., indx, nrec1, ndpart, iod(1,1), 1 )
#endif

	enddo
c
	deallocate (bfpos)
c
	cpu(4)=cpu(4)+MPI_WTIME()-rtime1
c
#ifdef TEST
#ifdef LAGSC2
#ifdef SEP0411
	if(i.eq.1) then
	nprec = ndpart2/nrec2

	allocate (bfpos3(3,nprec,nrec2))

        do ir=1,nrec2
        indx = (ir-1)*nprec/numtasks + 1

	call intbf_ars ( xyzl,nxyz,gsh(1,1),pp2(1,1),indx,nrec2,ndpart2,
     1               bfpos3(1,1,ir),slant, 2)

	enddo
	endif
#endif
#endif
#endif

#ifdef LAGSC2
	rtime1=MPI_WTIME()
        pp2(:,indexp)=0.
c
#ifdef SEP0411
	nprec = ndpart2/nrec2

	do ir=1,nrec2
        indx = (ir-1)*nprec/numtasks + 1

#ifdef CF_LAG
        call spcal_tree  ( bfpos3(1,1,ir), bs(1,1,1,1), pp2(1,indexp), jj,
     1             1., indx, nrec2, ndpart2, iod(1,1) )
#else
        call spcal_ars  ( bfpos3(1,1,ir), bs(1,1,1,1), pp2(1,indexp), jj,
     1             1., indx, nrec2, ndpart2, iod(1,1), 2 )
#endif

	enddo


#else

        nprec = ndpart2/nrec2
        allocate(bfpos(3,nprec))

        do ir=1,nrec2
        indx = (ir-1)*nprec/numtasks + 1

	call intbf_ars ( xyzl,nxyz,gsh(1,1),pp2(1,1),indx,nrec2,ndpart2,
     1               bfpos,slant, 2)

#ifdef CF_LAG
        call spcal_tree  ( bfpos, bs(1,1,1,1), pp2(1,indexp), jj,
     1             1., indx, nrec2, ndpart2, iod(1,1) )
#else
        call spcal_ars  ( bfpos, bs(1,1,1,1), pp2(1,indexp), jj,
     1             1., indx, nrec2, ndpart2, iod(1,1), 2 )
#endif
c
	enddo
c
	deallocate (bfpos)
#endif

	cpu(5)=cpu(5)+MPI_WTIME()-rtime1
#endif
c
	indexp=indexp+jj
c
c
c transform velocity back to wavenumber space
c
	rtime1=MPI_WTIME()
 	call kxtran (u(1,1,1,i),u(1,1,1,i),u(1,1,1,i),1)
	cpu(1)=cpu(1)+MPI_WTIME()-rtime1
c
c
 10	continue
c
	deallocate (udxzr)
c
	rtime1=MPI_WTIME()
c
c
c	if (jstep.eq.1) then
c
      ipvg(1)=ipj1(1+nc)
      ipvg(2)=ipj1(1+nc)+1
      ipvg(3)=ipj1(1+nc)+2
      ipvg(4)=ipj1(2+nc)
      ipvg(5)=ipj1(2+nc)+1
      ipvg(6)=ipj1(2+nc)+2
      ipvg(7)=ipj1(3+nc)
      ipvg(8)=ipj1(3+nc)+1
      ipvg(9)=ipj1(3+nc)+2
c
c
	sum=0.
      do i=1,nop/numtasks
      termu=pp(i,ipvg(1))**2+pp(i,ipvg(2))**2+pp(i,ipvg(3))**2
      termv=pp(i,ipvg(4))**2+pp(i,ipvg(5))**2+pp(i,ipvg(6))**2
      pp(i,ipvg(9))=-pp(i,ipvg(1))-pp(i,ipvg(5))
      termw=pp(i,ipvg(7))**2+pp(i,ipvg(8))**2+pp(i,ipvg(9))**2
      term=viscos*(termu+termv+termw)
      sum=sum+term
	end do
c
	call MPI_REDUCE (sum,sum0,1,mpireal,MPI_SUM,0,
     1                      MPI_COMM_WORLD,mpierr)
c
	if (taskid.eq.0) write (6,601) istep,epslon,sum0/nop
 601  format ('partvg: istep,<eps>,ave. pseudo-diss=',i6,1p,2e12.4)
c
#ifdef LAGSC2
        sum=0.
      do i=1,nop2/numtasks
      termu=pp2(i,ipvg(1))**2+pp2(i,ipvg(2))**2+pp2(i,ipvg(3))**2
      termv=pp2(i,ipvg(4))**2+pp2(i,ipvg(5))**2+pp2(i,ipvg(6))**2
      pp2(i,ipvg(9))=-pp2(i,ipvg(1))-pp2(i,ipvg(5))
      termw=pp2(i,ipvg(7))**2+pp2(i,ipvg(8))**2+pp2(i,ipvg(9))**2
      term=viscos*(termu+termv+termw)
      sum=sum+term
        end do
c
        call MPI_REDUCE (sum,sum0,1,mpireal,MPI_SUM,0,
     1                      MPI_COMM_WORLD,mpierr)
c
        if (taskid.eq.0) write (6,602) istep,epslon,sum0/nop2
 602  format ('part2vg: istep,<eps>,ave. pseudo-diss=',i6,1p,2e12.4)
c


#endif

c
#ifdef NOTYET
c
c velocity Laplacians if desired
c
      do 20 i=1,3
c
      if (lplapl(i).eq.0) go to 20
c
      indexp=indexp+1
      if (i.eq.1) ipj1(4+nc)=indexp
      ipj2(4+nc)=indexp
c
      call udxphy (uxy(1,1,1,i),uxzr(1,1,1,i),1,4)
c
      yg=taskid+1
      call spxyz ( udxzr, bsxy(1,1,1,1), 1, 1, nx+2, yg )
c
      zg=taskid+1
      iz1=sz1(zg)
      iz2=sz1(zg)+nozp(zg)-1
      call spcal ( ibxyz, bfpos, bsxy(1,1,1,1), pp(1,indexp),
     1             1., nop, ndpart, iod(1,1), iz1, iz2 )
c
 20   continue
c 
c form the mechanical dissipation if desired
c
      if (ipdis(1).gt.0) then
c
      indexp=indexp+1
      ipdis(1)=indexp
c
      do 30 i=1,nop
c
      termu=pp(i,ipvg(1))**2+pp(i,ipvg(2))**2+pp(i,ipvg(3))**2
      termv=pp(i,ipvg(4))**2+pp(i,ipvg(5))**2+pp(i,ipvg(6))**2
      pp(i,ipvg(9))=-pp(i,ipvg(1))-pp(i,ipvg(5))
      termw=pp(i,ipvg(7))**2+pp(i,ipvg(8))**2+pp(i,ipvg(9))**2
      pp(i,ipdis(1))=viscos*(termu+termv+termw)
c
 30   continue
c
      end if
c
      do i=1,3
      if (lpgrad(i).eq.1.or.lplapl(i).eq.1)
     1       call kxtran (uxy(1,1,1,i),uxz(1,1,1,i))
      end do
c
      icall=1
#endif
	cpu(6)=cpu(6)+MPI_WTIME()-rtime1
c
#endif

	cpu_overall=MPI_WTIME()-rtime0
c
	if (taskid.eq.0) write (6,*) 'exit partvg, istep=',istep
	if (taskid.eq.1) then
	write (71,"(i6,1p,7e10.3)") istep,cpu,cpu_overall
	end if
#endif
      return
      end
