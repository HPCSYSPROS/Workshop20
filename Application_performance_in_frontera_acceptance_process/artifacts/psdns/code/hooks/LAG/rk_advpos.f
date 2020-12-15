       subroutine rk_advpos (pstep,cof1,cof2)
c
c
#ifdef LAG
c
	use compart
	implicit none
	include 'intvars'
c
	real(4) cof1,cof2
	real(b8) fact1,fact2
	integer ispc, np,np2,i,k,pstep
c
	ps(1)=sqrt(b11(2))
	ps(2)=sqrt(b22(2))
	ps(3)=sqrt(b33(2))
c
	np=nop/numtasks
	np2=nop2/numtasks
c
        if (taskid.eq.0) write (10000,"('calling rk_advpos, istep,kstep,pstep=',3i4)") istep,kstep,pstep
c
	if (pstep.eq.1) then
c
	do 10 k=1,3
c
	fact1 = cof1 * dt/2./pi*xyzl(k)*ps(k)
	fact2 = cof2 * dt/2./pi*xyzl(k)*ps(k)
c
	do i=1,np
      pp(i,6+k) = pp(i,  k) + fact1*pp(i,3+k)
      pp(i,  k) = pp(i,6+k) + fact2*pp(i,3+k)
	end do
c
#ifdef LAGSC2
	do i=1,np2
      pp2(i,6+k) = pp2(i,  k) + fact1*pp2(i,3+k)
      pp2(i,  k) = pp2(i,6+k) + fact2*pp2(i,3+k)
	end do
#endif
c
 10	continue
c
	else if (pstep.eq.rkmethod) then
c
	call ppminmax ('advpos: before do 20',pp(1,1),ndpart)
	call ppminmax ('advpos: before do 20',pp(1,4),ndpart)
	call ppminmax ('advpos: before do 20',pp(1,7),ndpart)

	do 20 k=1,3
	fact2 = cof2 * dt/2./pi*xyzl(k)*ps(k)
	do i=1,np
      pp(i,  k) = pp(i,6+k) + fact2*pp(i,3+k)
	end do
 20     continue
	call ppminmax ('advpos:  after do 20',pp(1,1),ndpart)
c
#ifdef LAGSC2
	do 25 k=1,3
	fact2 = cof2 * dt/2./pi*xyzl(k)*ps(k)
	do i=1,np2
      pp2(i,  k) = pp2(i,6+k) + fact2*pp2(i,3+k)
	end do
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
