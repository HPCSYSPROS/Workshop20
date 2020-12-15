

!################################################
!#   Subroutine ITRANSFORM
!#
!#   Performs an inverse transform on X and Z 
!#        (return to physical space)
!#  This version transforms both the velocity and vorticity
!#
!#
!################################################
	subroutine itransform_om (ux,uy,uz,omx,omy,omz)
#ifdef VORPHY
          use comsp
	use timers_rkstep
          implicit none
          include 'intvars'
#include "fft_stuff.f"

          real(b8) :: ux(nxpad,zisz,yjsz,nu)
          complex(b8) :: uz(nzpad,xisz,yjsz,nu)
          complex(b8) :: uy(nypad,zjsz*xisz,nu)
          real(b8) :: omx(nxpad,zisz,yjsz,3)
          complex(b8) :: omz(nzpad,xisz,yjsz,3)
          complex(b8) :: omy(nypad,zjsz*xisz,3)
          integer :: m,i,j,zc,iz,int,ifc,ilc,idflag,x2,ix,dnz,iy,l

	real*8 rtime1,rtime2

        integer iy1,iy2,ithr,omp_get_thread_num,ix1
          
	if (taskid.eq.0) write (6,*) 'enter itransform_om, istep,kstep=',istep,kstep

	rtime1=MPI_WTIME()
c
	if (kstep.eq.1) t_itrans(:,2)=0.

! Transform in Y, transpose to get data in z-pencils, add extra elements if needed
c
	rtime2=MPI_WTIME()

#ifdef MULTIVAR

#ifdef CLEAN
	     call kxcomm1_clean (uy,uz,3)
	     call kxcomm1_clean (omy,omz,3)
#else
 	     call kxcomm1_trans_cyl (uy,uz,3)
 	     call kxcomm1_trans_cyl (omy,omz,3)
#endif

#else

#ifdef CLEAN
             do i=1,3
	       call kxcomm1_clean (uy(1,1,i),uz(1,1,1,i),1)
	       call kxcomm1_clean (omy(1,1,i),omz(1,1,1,i),1)
             end do
#else
             do i=1,3
 	       call kxcomm1_trans_cyl (uy(1,1,i),uz(1,1,1,i),1)
 	       call kxcomm1_trans_cyl (omy(1,1,i),omz(1,1,1,i),1)
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
c
	     do i=1,3
	do yp=iyp1(ithr),iyp1(ithr)+iypi(ithr)-1
	do xp=1,xisz
	uz(nzhp,xp,yp,i)=0.
	end do
	end do
	iy1=iyp1(ithr)
       call fftw_execute (plan2_kx(ithr),uz(1,1,iy1,i),uz(1,1,iy1,i))
       call fftw_execute (plan2_kx(ithr),omz(1,1,iy1,i),omz(1,1,iy1,i))
	     enddo
c
!$OMP END PARALLEL
	rtime2=MPI_WTIME()-rtime2
	t_itrans(2,2)=t_itrans(2,2)+rtime2

	rtime2=MPI_WTIME()

#ifdef KXCOMM2

#ifdef MULTIVAR
 	call kxcomm2_trans (uz(1,1,1,1),ux(1,1,1,1),3)
 	call kxcomm2_trans (omz(1,1,1,1),omx(1,1,1,1),3)
#else
        do i=1,3
 	  call kxcomm2_trans (uz(1,1,1,i),ux(1,1,1,i),1)
 	  call kxcomm2_trans (omz(1,1,1,i),omx(1,1,1,i),1)
        end do
#endif

#else

	  do 10 yp=1,yjsz
#ifdef MULTIVAR
 	call kxcomm2t_trans (uz(1,1,1,1),ux(1,1,1,1),yjsz,yp,3)
 	call kxcomm2t_trans (omz(1,1,1,1),omx(1,1,1,1),yjsz,yp,3)
#else
        do i=1,3
 	  call kxcomm2t_trans (uz(1,1,1,i),ux(1,1,1,i),yjsz,yp,1)
 	  call kxcomm2t_trans (omz(1,1,1,i),omx(1,1,1,i),yjsz,yp,1)
        end do
#endif
 10	continue

#endif

	rtime2=MPI_WTIME()-rtime2
	t_itrans(3,2)=t_itrans(3,2)+rtime2

	rtime1=MPI_WTIME()-rtime1
              
c
#endif
	return
        end subroutine itransform_om

