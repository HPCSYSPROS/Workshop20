       subroutine rk4_advpos (pstep,cof1,cof2)
c
c
#ifdef LAG
c
	use compart
	implicit none
	include 'intvars'
c
	real cof1,cof2,fact1,fact2
	integer ispc, np,np2,i,k,pstep,kk
c
	ps(1)=sqrt(b11(2))
	ps(2)=sqrt(b22(2))
	ps(3)=sqrt(b33(2))
c
	np=nop/numtasks
	np2=nop2/numtasks
c
	if (taskid.eq.0) write (10000,*) 'calling rk_advpos',
     1              dt/2./pi*xyzl(1)*ps(1)
c
	if (pstep.lt.rkmethod) then
c
	do 10 k=1,3
c
	fact1 = cof1 * dt/2./pi*xyzl(k)*ps(k)
	fact2 = cof2 * dt/2./pi*xyzl(k)*ps(k)
c
	kk=ndprop-3+k
	do i=1,np
      pp(i,  k) = pp(i, kk) + fact1*pp(i,3+k)
      pp(i,6+k) = pp(i,6+k) + fact2*pp(i,3+k)
	end do
c
#ifdef LAGSC2
	kk=ndprop2-3+k
	do i=1,np2
      pp2(i,  k) = pp2(i, kk) + fact1*pp2(i,3+k)
      pp2(i,6+k) = pp2(i,6+k) + fact2*pp2(i,3+k)
	end do
#endif
c
 10	continue
c
	else if (pstep.eq.rkmethod) then
c
	do 20 k=1,3
c
	fact2 = cof2 * dt/2./pi*xyzl(k)*ps(k)
c
	do i=1,np
      pp(i,  k) = pp(i,6+k) + fact2*pp(i,3+k)
	end do
c
#ifdef LAGSC2
	do i=1,np2
      pp2(i,  k) = pp2(i,6+k) + fact2*pp2(i,3+k)
	end do
#endif
c
 20	continue
c
	end if

c
#endif
c
      return
      end
