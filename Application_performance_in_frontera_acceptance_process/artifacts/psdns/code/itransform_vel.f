

!################################################
!#   Subroutine ITRANSFORM
!#
!#   Performs an inverse transform on X and Z 
!#        (return to physical space)
!#
!#
!################################################
	subroutine itransform_vel (ux,uy,uz)
          use comsp
	use timers_rkstep
	use timers_tran
          implicit none
          include 'intvars'
#include "fft_stuff.f"

          real(b8) :: ux(nxpad,zisz,yjsz,nu)
          complex(b8) :: uz(nzpad,xisz,yjsz,nu)
          complex(b8) :: uy(nypad,zjsz*xisz,nu)
          integer :: m,i,j,zc,iz,int,ifc,ilc,idflag,x2,ix,dnz,iy,l

	real*8 rtime1,rtime2

        integer iy1,iy2,ithr,omp_get_thread_num,ix1
        real(b8) norm
          
	if (taskid.eq.0.and.istep.eq.1) write (6,*) 'enter itransform_vel, kstep=',kstep

	rtime1=MPI_WTIME()
c
        if (kstep.eq.1) then
        t_itrans(:,2)=0.
        t_kxcomm1(:,2)=0.
        t_kxcomm2(:,2)=0.
        t_kxcomm2t(:,2)=0.
        end if

! Transform in Y, transpose to get data in z-pencils, add extra elements if needed
c
	rtime2=MPI_WTIME()
#ifdef CLEAN
#ifdef MULTIVAR
	     call kxcomm1_clean (uy,uz,3)
#else
             do i=1,3
	       call kxcomm1_clean (uy(1,1,i),uz(1,1,1,i),1)
             end do
#endif
#else
#ifdef MULTIVAR
 	     call kxcomm1_trans_cyl (uy,uz,3)
#else
             do i=1,3
 	       call kxcomm1_trans_cyl (uy(1,1,i),uz(1,1,1,i),3)
             end do
#endif
#endif
	rtime2=MPI_WTIME()-rtime2
	t_itrans(1,2)=t_itrans(1,2)+rtime2


	rtime2=MPI_WTIME()

        call divide_thr (yjsz,iyp1,iypi,'itransform: yjsz')
c
        ithr=0
!$OMP PARALLEL private(ithr,iy1)

#ifdef OPENMP
      ithr = OMP_GET_THREAD_NUM()
#endif
#ifdef ESSL
c        call cft (1,uz,1,nzpad,uz,1,nzpad,nzpad,xisz,-1,norm,
c    1            caux1,cnaux1,caux2,cnaux2)
         call escpft (uz,1,nzpad,nzpad,xisz*yjsz/num_thr,1,ithr)
#endif
c
	     do i=1,3
	do yp=iyp1(ithr),iyp1(ithr)+iypi(ithr)-1
	do xp=1,xisz
	uz(nzhp,xp,yp,i)=0.
	end do
	end do
	iy1=iyp1(ithr)
#ifdef FFTW
       call fftw_execute (plan2_kx(ithr),uz(1,1,iy1,i),uz(1,1,iy1,i))
#elif defined ESSL
c       call cft (0,uz(1,1,iy1,i),1,nzpad,uz(1,1,iy1,i),1,nzpad,nzpad,
c    1        xisz,-1,norm,caux1,cnaux1,caux2,cnaux2)
         call escfft (uz(1,1,iy1,i),1,nzpad,nzpad,xisz*yjsz/num_thr,1,ithr)
#endif
	     enddo
c
!$OMP END PARALLEL
	rtime2=MPI_WTIME()-rtime2
	t_itrans(2,2)=t_itrans(2,2)+rtime2

	rtime2=MPI_WTIME()
#ifdef KXCOMM2
#ifdef MULTIVAR
 	call kxcomm2_trans (uz(1,1,1,1),ux(1,1,1,1),3)
#else
        do i=1,3
 	  call kxcomm2_trans (uz(1,1,1,i),ux(1,1,1,i),1)
       end do
#endif
#else
#ifdef MULTIVAR
	do yp=1,yjsz
 	  call kxcomm2t_trans (uz(1,1,1,1),ux(1,1,1,1),yjsz,yp,3)
        end do
#else
        do yp=1,yjsz
          do i=1,3
            call kxcomm2t_trans (uz(1,1,1,i),ux(1,1,1,i),yjsz,yp,1)
          end do
        end do
#endif
#endif
	rtime2=MPI_WTIME()-rtime2
	t_itrans(3,2)=t_itrans(3,2)+rtime2

	rtime1=MPI_WTIME()-rtime1
              
c
	return
        end subroutine itransform_vel

