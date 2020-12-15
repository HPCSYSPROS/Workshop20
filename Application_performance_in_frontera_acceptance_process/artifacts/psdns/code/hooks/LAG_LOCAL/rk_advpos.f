       subroutine rk_advpos (pstep,cof1,cof2)
c
c
#ifdef LAG
c
	use compart
	implicit none
c
	real(4) cof1,cof2
	real(b8) fact1,fact2
	integer ispc, np,np2,i,k,pstep
	integer tid,nid,niderr,niderr0

	integer, save :: np_local
c
	ps(1)=sqrt(b11(2))
	ps(2)=sqrt(b22(2))
	ps(3)=sqrt(b33(2))
c
	np=nop/numtasks
!	np2=nop2/numtasks
c
        if (taskid.eq.0) write (10000,"('calling rk_advpos, istep,kstep,pstep=',3i4)") istep,kstep,pstep
c
	if (pstep.eq.1) then

	np_local = nump2
c
	do 10 k=1,3
c
	fact1 = cof1 * dt/2./pi*xyzl(k)*ps(k)
	fact2 = cof2 * dt/2./pi*xyzl(k)*ps(k)
c
	if (ndpart.gt.0) then
	do i=1,nump1
      pp(i,6+k) = pp(i,  k) + fact1*pp(i,3+k)
      pp(i,  k) = pp(i,6+k) + fact2*pp(i,3+k)
	end do
	endif
c
#ifdef LAGSC2
	if (ndpart2.gt.0) then
	do i=1,nump2
      pp2(i,6+k) = pp2(i,  k) + fact1*pp2(i,3+k)
      pp2(i,  k) = pp2(i,6+k) + fact2*pp2(i,3+k)
	end do
	endif
#endif
c
 10	continue
c
	else if (pstep.eq.rkmethod) then

	if(np_local.ne.nump2) then
	write (6,*) 'nump2 cannot change between 2 calls of rk_part, nump2,npl=',nump2,np_local
	stop
	endif
c
	if (ndpart.gt.0) then
	call ppminmax ('advpos: before do 20',pp(1,1),ndpart)
	call ppminmax ('advpos: before do 20',pp(1,4),ndpart)
	call ppminmax ('advpos: before do 20',pp(1,7),ndpart)
	endif

	do 20 k=1,3
	fact2 = cof2 * dt/2./pi*xyzl(k)*ps(k)

	if (ndpart.gt.0) then
	do i=1,nump1
      pp(i,  k) = pp(i,6+k) + fact2*pp(i,3+k)
	end do
	endif

 20     continue
	if (ndpart.gt.0) call ppminmax ('advpos:  after do 20',pp(1,1),ndpart)
c
#ifdef LAGSC2
	do 25 k=1,3
	fact2 = cof2 * dt/2./pi*xyzl(k)*ps(k)

	if(ndpart2.gt.0) then
	do i=1,nump2
      pp2(i,  k) = pp2(i,6+k) + fact2*pp2(i,3+k)
	end do
	endif

	niderr=0

	if(nump1.gt.0) then
	do i=1,nump1
	call pptoid (pp(i,1), pp(i,3), tid, nid, 1000.)

	if(nid.lt.0 .or. nid.gt.8) then
	write(6,*) 'niderr: i,nump1=',i,nump1
	write(6,*) 'niderr: taskid, tid, nid',taskid, tid, nid
	write(6,*) 'niderr: taskid, pp13=',taskid, pp(i,1:3)
	write(6,*) 'niderr: taskid, pp79=',taskid, pp(i,7:9)
	write(6,*) 'niderr: taskid, pp46=',taskid, pp(i,4:6)
	write(6,*) 'niderr: taskid, dt,=',taskid, dt, fact2
	write(6,*) 'niderr: taskid, nbtask,=',taskid, nbtask(1:8)
	niderr = niderr + 1
	endif
	enddo
	endif

	if(nump2.gt.0) then
        do i=1,nump2
        call pptoid (pp2(i,1), pp2(i,3), tid, nid, 1000.)

        if(nid.lt.0 .or. nid.gt.8) then
	write(6,*) 'niderr: i,nump2=',i,nump2
	write(6,*) 'niderr: taskid, tid, nid',taskid, tid, nid
        write(6,*) 'niderr: taskid, tid, nid',taskid, tid, nid
        write(6,*) 'niderr: taskid, pp13=',taskid, pp2(i,1:3)
        write(6,*) 'niderr: taskid, pp79=',taskid, pp2(i,7:9)
        write(6,*) 'niderr: taskid, pp46=',taskid, pp2(i,4:6)
        write(6,*) 'niderr: taskid, dt,=',taskid, dt, fact2
        write(6,*) 'niderr: taskid, nbtask,=',taskid, nbtask(1:8)
        niderr = niderr + 1
        endif
        enddo
	endif

        call MPI_ALLREDUCE (niderr,niderr0,1,MPI_INTEGER,MPI_SUM,
     1                   MPI_COMM_WORLD,mpierr)


	if(niderr0.gt.0) then
	if(taskid.eq.0) write(6,*) 'abort in rk_advpos, niderr0=',niderr0
	call MPI_ABORT (MPI_COMM_WORLD, mpierr)
	endif

 25	continue
#endif
c
c
	end if

c
#endif
c
      return
      end
