!######################################################################
!#       Transform Routines
!# Takes as input the velocity fields
!# Transforms to fourier space
!#
!# Only transforming / transposing two components of velocity (x,z)
!#
!######################################################################

       subroutine transform_vel (ux,uy,uz,m) 
	use comp
	use timers_rkstep
        use timers_tran
	implicit none
	include 'intvars'
#include "fft_stuff.f"
 	real(b8)    :: ux(nxpad,zisz,yjsz,nu)
 	complex(b8) :: uz(nzpad,xisz,yjsz,nu)
 	complex(b8) :: uy(nypad,zjsz*xisz,nu)
	integer :: m,i, is1,is2,novt,rkstep,i2f,l
        real(b8) s1,s2,s3
	complex(b8), allocatable :: bk1i(:,:),bk3i(:,:)
c
	integer ithr,omp_get_thread_num,yp1

	integer i1,ii
	real(8) rtime1,rtime2

        if (kstep.eq.1) then
        t_trans(:,2)=0.
        t_xkcomm1(:,2)=0.
        t_xkcomm2(:,2)=0.
        end if

        s1=sqrt(b11(m))
        s2=sqrt(b22(m))
        s3=sqrt(b33(m))
        is1=4+nc
        is2=5+nc

! transform to x-wavenumber space, and take transpose 
!

	novt=5+nc
#ifdef ROTD
        if (kstep.eq.1) novt=6+nc
#endif
#ifdef CVC_PRESS
        if (flagc) novt=6+nc
#endif


	rtime1=MPI_WTIME()
#ifdef MULTIVAR
           call xkcomm1_trans(ux,uz,novt)
#else
           do i=1,novt
             call xkcomm1_trans(ux(1,1,1,i),uz(1,1,1,i),1)
           end do
#endif
	t_trans(1,2)=t_trans(1,2)+(MPI_WTIME()-rtime1)
c
	call divide_thr (yjsz,iyp1,iypi,'transform: yjsz')

! these new arrays introduced by PKY, 4/9/2012
	allocate (bk1i(xist:xien,2))
	allocate (bk3i(nz,2))

	ithr=0
        
!$OMP PARALLEL private (ithr,yp1,x,bk1i,bk3i)
c
#ifdef OPENMP
      ithr = OMP_GET_THREAD_NUM()
#endif
#ifdef ESSL
        call escpft (uz,1,nzpad,nzpad,xisz*yjsz/num_thr,-1,ithr)
#endif
c
        yp1=iyp1(ithr)
!
!$OMP MASTER
	rtime1=MPI_WTIME()
!$OMP END MASTER

#ifdef FFTW
            do i=1,novt
         if (i.ne.1.and.i.ne.3.and.i.ne.is2) then
         call fftw_execute(plan2_xk(ithr),uz(1,1,yp1,i),uz(1,1,yp1,i))
         end if
         end do

         call fftw_execute(plan2_xk(ithr),uz(1,1,yp1,1),uz(1,1,yp1,1))
         call fftw_execute(plan2_xk(ithr),uz(1,1,yp1,3),uz(1,1,yp1,3))
         call fftw_execute(plan2_xk(ithr),uz(1,1,yp1,is2),uz(1,1,yp1,is2))
#elif defined ESSL
           do i=1,novt
         if (i.ne.1.and.i.ne.3.and.i.ne.is2) then
         call escfft (uz(1,1,yp1,i),1,nzpad,nzpad,xisz*yjsz/num_thr,-1,ithr)
         end if
         end do

         call escfft (uz(1,1,yp1,1),1,nzpad,nzpad,xisz*yjsz/num_thr,-1,ithr)
         call escfft (uz(1,1,yp1,3),1,nzpad,nzpad,xisz*yjsz/num_thr,-1,ithr)
         call escfft (uz(1,1,yp1,is2),1,nzpad,nzpad,xisz*yjsz/num_thr,-1,ithr)
#endif

!$OMP MASTER
	t_trans(2,2)=t_trans(2,2)+(MPI_WTIME()-rtime1)
	rtime1=MPI_WTIME()
!$OMP END MASTER


	do x=xist,xien
	bk1i(x,1)=imagi*kx(x)
	bk1i(x,2)=imagi*bk1(x,m)
	end do
c
	do z=1,nz
	bk3i(z,1)=imagi*bk3(z,m)
	bk3i(z,2)=imagi*kz(z)
	end do

!#################################################
! operations for convective terms
!#################################################
c
        do 70 yp=iyp1(ithr),iyp1(ithr)+iypi(ithr)-1
           
	y=yjst+yp-1
           do 71 xp=1,xisz
              x=xist+xp-1
#ifdef ROTD
              if (kstep.eq.1) then
                 do 75 z=1,nzpad
               if (z.eq.nzhppad) go to 75
               uz(z,xp,yp,1)=bk1i(x,1)*uz(z,xp,yp,1)+bk3i(z,1)*uz(z,xp,yp,5+nc)
               uz(z,xp,yp,3)=bk1i(x,2)*uz(z,xp,yp,5+nc)+bk3i(z,2)*uz(z,xp,yp,3)
               uz(z,xp,yp,5+nc)=uz(z,xp,yp,2)
               uz(z,xp,yp,2)=bk1i(x,1)*uz(z,xp,yp,2)+bk3i(z,1)*uz(z,xp,yp,4+nc)
 75              continue
              else if (kstep.eq.2) then
                 do 76 z=1,nzpad
                    if (z.eq.nzhppad) go to 76
                 uz(z,xp,yp,1)=bk1i(x,1)*uz(z,xp,yp,1)+bk3i(z,1)*uz(z,xp,yp,is2)
                 uz(z,xp,yp,3)=bk1i(x,2)*uz(z,xp,yp,is2)+bk3i(z,2)*uz(z,xp,yp,3)
 76              continue
              end if
#elif defined CVC_PRESS
        if (flagc) then
                 do 75 z=1,nzpad
                    if (z.eq.nzhppad) go to 75
                    uz(z,xp,yp,1)=bk1i(x,1)*uz(z,xp,yp,1)+bk3i(z,1)*uz(z,xp,yp,5+nc)
                    uz(z,xp,yp,3)=bk1i(x,2)*uz(z,xp,yp,5+nc)+bk3i(z,2)*uz(z,xp,yp,3)
                    uz(z,xp,yp,5+nc)=uz(z,xp,yp,2)
                    uz(z,xp,yp,2)=bk1i(x,1)*uz(z,xp,yp,2)+bk3i(z,1)*uz(z,xp,yp,4+nc)
 75              continue
        else
              do 76 z=1,nzpad
                 if (z.eq.nzhppad) go to 76
                 uz(z,xp,yp,1)=bk1i(x,1)*uz(z,xp,yp,1)+bk3i(z,1)*uz(z,xp,yp,is2)
                 uz(z,xp,yp,3)=bk1i(x,2)*uz(z,xp,yp,is2)+bk3i(z,2)*uz(z,xp,yp,3)
 76           continue
        end if
#else
              do 76 z=1,nzpad
                 if (z.eq.nzhppad) go to 76
                 uz(z,xp,yp,1)=bk1i(x,1)*uz(z,xp,yp,1)+bk3i(z,1)*uz(z,xp,yp,is2)
                 uz(z,xp,yp,3)=bk1i(x,2)*uz(z,xp,yp,is2)+bk3i(z,2)*uz(z,xp,yp,3)
 76           continue
#endif

 71        continue

 70   continue

!$OMP MASTER
	t_trans(3,2)=t_trans(3,2)+(MPI_WTIME()-rtime1)
!$OMP END MASTER


!$OMP END PARALLEL
	deallocate (bk1i,bk3i)

!###################################################
        
      i2f=4+nc
#ifdef ROTD
        if (kstep.eq.1) i2f=6+nc
#endif
#ifdef CVC_PRESS
        if (flagc) i2f=6+nc
#endif


! Truncate in Z if needed, Transpose ZXY -> YZX (Z- to Y-pencils), transform in Y 

	rtime1=MPI_WTIME()
!RAF: Always call xkcomm2_clean so that CAFI will work (similar to kxcomm1_clean).
!         if(num_al_i(0) .gt. 0) then
c            call xkcomm2_trans_cyl (uz,uy,i2f)
#ifdef MULTIVAR
            call xkcomm2_clean (uz,uy,i2f)
#else
            do i=1,i2f
              call xkcomm2_clean (uz(1,1,1,i),uy(1,1,i),1)
            end do
#endif
!	endif
	t_trans(4,2)=t_trans(4,2)+(MPI_WTIME()-rtime1)


c
      return
      end subroutine transform_vel
