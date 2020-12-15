      subroutine part2sc (indexp)
c
c  routine to perform cubic spline interpolation for 
c  scalar fluctuations following fluid particles
c
c scalar gradients also implemented, 3/1/97
c
c  called by sub. partsp
c
#ifdef NOTYET
#ifndef NOSCALAR
#ifdef LAG
#ifdef LGRAD
#ifdef LAGSC2
c
        use compart
        include 'intvars'
c
	real sum,sm0
	real sg(3)
c
	integer icall
	data icall/0/
c
c
c in the current formulation, since each scalar is processed
c sequentially, all the scalars can share one storage position
c in the (bsxy) spline array
c
!	allocate (udxzr(nx,zisz,yjsz))
	allocate (udxzr(nx,zisz,yjsz,3))
c
	sg(1)=1./beta1
	sg(2)=1./beta2
	sg(3)=1./beta3
c
	icall=icall+1
c
      do 10 i=1,nc
c
      ipj1(i)=indexp+1
      ipj2(i)=indexp+1
c
	jc=i+3+nc
c
      if (lpfunc(i+3).eq.1) then
c
      indexp=indexp+1
c
c
      yg=taskid+1
!      call spxyz ( u(1,1,1,jc), bs(1,1,1,1), bs(1,1,1,1), nx )
      call spxyz ( u(1,1,1,jc), bs(1,1,1,1), 1 )
c
      zg=taskid+1
      call spcal_ars ( bfpos2, bs(1,1,1,1), pp2(1,indexp), 1,
     1             1., nop2, ndpart2, iod(1,1) )
c
c      if (.not.(jstep.eq.1.and.kinit(i+3).eq.0))
c    1   call pitestc (bs(1,1,1,1),jc)
c
      end if
c
      if (lpgrad(i+3).eq.1) then
c
      ipj2(i)=indexp+3
c
!      call xktran (u(1,1,1,jc),u(1,1,1,jc),u(1,1,1,jc))
      call xktran (u(1,1,1,jc),u(1,1,1,jc),u(1,1,1,jc),1)
c
c derivates in x,y,z directions
c
c
	do 15 k=1,3
c
c
      indexp=indexp+1
c
      call udxphy (u(1,1,1,jc),udxzr,udxzr,k)
c
      yg=taskid+1
!      call spxyz ( udxzr, bs(1,1,1,1), bs(1,1,1,1),nx )
      call spxyz ( udxzr, bs(1,1,1,1), 1 )
c
      zg=taskid+1
      call spcal_ars ( bfpos2, bs(1,1,1,1), pp2(1,indexp), 1,
     1             sg(k), nop2, ndpart2, iod(1,1) )
c
	if (jstep.eq.1) then
      sum=0.
      do 390 ii=1,nop2/numtasks
      sum=sum+pp2(ii,indexp)**2
 390  continue
        call MPI_REDUCE (sum,sum0,1,mpireal,MPI_SUM,0,
     1                   MPI_COMM_WORLD,mpierr)
        sum=sum0
	if (taskid.eq.0) write (6,*) 'from part2sc, mean-sq=',sum/nop2
	end if

 15     continue
c
c
c form scalar dissipation if desired
c
      if (ipdis(i+1).gt.0) then
c
      indexp=indexp+1
      ipdis(i+1)=indexp
      factor=viscos/pr(i)
c
      do 20 ii=1,nop2
      pp2(ii,indexp)=factor*(pp2(ii,indexp-3)**2+pp2(ii,indexp-2)**2
     1                      +pp2(ii,indexp-1)**2)
 20   continue
c
      end if
c
      end if
c
#ifdef NOTYET
c Laplacian of scalar if desired
c
      if (lplapl(j).eq.1) then
c
c
      ipj2(i)=indexp+1
c
      call udxphy (uxy(1,1,1,i+3+nc),uxzr(1,1,1,i+3+nc),1,4)
c
      indexp=indexp+1
c
      yg=taskid+1
      call spxyz ( udxzr, bsxy(1,1,1,1), 1, 1, nx+2, yg )
c
      zg=taskid+1
      iz1=sz1(zg)
      iz2=sz1(zg)+nozp(zg)-1
      call spcal_ars ( ibxyz2, bfpos2, bsxy(1,1,1,1), pp2(1,indexp),
     1             1., nop2, ndpart2, iod(1,1), iz1, iz2 )
c
      end if
c
#endif
c
      if (lpgrad(i+3).eq.1.or.lplapl(i+3).eq.1)  then
!      call kxtran (u(1,1,1,jc),u(1,1,1,jc),u(1,1,1,jc))
      call kxtran (u(1,1,1,jc),u(1,1,1,jc),u(1,1,1,jc),1)
      end if
c
c
 10   continue
c
	deallocate (udxzr)
#endif
#endif
#endif
#endif
#endif
      return
      end
