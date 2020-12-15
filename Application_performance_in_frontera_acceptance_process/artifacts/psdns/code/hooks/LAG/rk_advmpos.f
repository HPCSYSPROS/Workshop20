       subroutine rk_advmpos (pstep)
c

c Several changes made by PKY, 2/1/2014
c 
c OpenMP now used on the loops in this routine --- may be advantageous
c when number of molecules gets really big?
c
c Also assume we are tracking molecules only in flow without mean shear
c
c 2/2/2014: each MPI task (and its master thread) handles 
c its own random number sequence for the Brownian motion.
c Quantities for molecules will depend on number of MPI tasks
c (but not number of threads per MPI task)
c
#ifdef LAG
#ifdef MOL
c
	use compart
	implicit none
c
        real moldif,mvrms,fact1,fact2,factor,hdx,hdt
	real work
        integer i,k,igm,im1,im2,inc,nsamp,j
	integer pstep,itemp,isz1,isz2,isum
        real(p8), allocatable :: mwork(:)
c
        real*8 rtime1,rtime2
        integer ithr,omp_get_thread_num,ii
c
	integer, allocatable :: imp1(:),impi(:)

c
#ifndef MOLOLD
!use the last saved mol seed 
        call random_seed(put=mseed)      
#endif
        rtime1=MPI_WTIME()

	if (nm.eq.0) go to 90
c
	ps(1)=sqrt(b11(2))
	ps(2)=sqrt(b22(2))
	ps(3)=sqrt(b33(2))
c
	allocate (imp1(0:num_thr-1),impi(0:num_thr-1))


!	call divide_thr (nm/ngm, imp1, impi, 'rk_advmpo: nm/ngm')

        itemp = nm/ngm/num_thr/nompp
        isz1 = (itemp+1)*nompp
        isz2 = nm/ngm - (num_thr-1)*isz1

        impi(:) = isz1
        impi(num_thr-1) = isz2

        imp1(0)=1
        isum=imp1(0)
        do i=1,num_thr-1
        isum=isum+impi(i-1)
        imp1(i)=isum
        end do

c
	if (taskid.eq.0) write (10000,"('calling rk_advmpos, istep,kstep,pstep=',3i4)") istep,kstep,pstep
c

	hdt = .5 * dt /2./pi

	if (pstep.eq.1) then
c
        nsamp=nm/ngm/nompp
        allocate (mwork(nsamp))
c
        ithr=0
c
!$OMP PARALLEL private (ithr,fact1,igm,im1,im2,moldif,factor,hdx,ii) shared(mwork)

#ifdef OPENMP
      ithr = OMP_GET_THREAD_NUM()
#endif

c
      do 220 igm=1,ngm
c
      moldif=viscos/molsc(igm)
#ifdef MOLOLD
!$OMP MASTER
      call ranseq (12+igm,0)
!$OMP END MASTER
#endif
c
	im1=imp1(ithr)
	im2=imp1(ithr)+impi(ithr)-1
c
      do k=1,3
c
      factor=sqrt(2.*moldif*dt)/(2.*pi/nxyz(k)/ps(k))
      fact1=hdt*xyzl(k)*ps(k)
c
        do j=1,nompp
c
! Apr 18 2014, mwork now declared to be double precision and 
! random_number intrinsic is used to generate the random numbers
!$OMP MASTER
#ifdef MOLOLD
       call rann2 (mwork,1,nsamp,1,1,nsamp)
#else
       call random_number(mwork)
	call boxmul (mwork,nsamp)
!       call molrann (mwork,1,nsamp,1,1,nsamp)
#endif
!$OMP END MASTER
!$OMP BARRIER

c
      do 225 i=im1+j-1,im2-nompp+j,nompp
c
        ii=(i-j)/nompp+1
      hdx = fact1 * mpp(i,3+k,igm)
c
      mpp(i,6+k,igm) = mpp(i,  k,igm) + hdx + factor * mwork(ii)
      mpp(i,  k,igm) = mpp(i,6+k,igm) + hdx
 225  continue
c
c!$OMP BARRIER
      end do
c
      end do
c
 220  continue

!$OMP END PARALLEL
c
        deallocate (mwork)
c
c
	else if (pstep.eq.rkmethod) then
c
c
c take care of the molecules
c
        ithr=0
c
!$OMP PARALLEL private (ithr,fact1,fact2,igm,im1,im2)

#ifdef OPENMP
      ithr = OMP_GET_THREAD_NUM()
#endif
c
	im1=imp1(ithr)
	im2=imp1(ithr)+impi(ithr)-1

        do igm=1,ngm
      do k=1,3
      fact1 = hdt * xyzl(k) * ps(k)
      do i=im1,im2
      mpp(i,k,igm)=mpp(i,6+k,igm) + fact1 * mpp(i,3+k,igm) 
      end do
      end do
        end do
c
!$OMP END PARALLEL

	end if

	deallocate (imp1,impi)

#ifndef MOLOLD
! save the current mol seed 
        call random_seed(get=mseed)
#endif
c
        rtime1=MPI_WTIME()-rtime1
        if (taskid.eq.0.and.istep.le.5) then
        write (6,"('rk_advmpos istep,pstep,cpu=',i6,i3,1p,e12.4)")
     1            istep,pstep,rtime1
        end if

 90    continue
c
#endif
#endif
c
      return
      end
