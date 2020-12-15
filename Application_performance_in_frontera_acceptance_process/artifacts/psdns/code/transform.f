!######################################################################
!#       Transform Routines
!# Takes as input the velocity fields
!# Transforms to fourier space
!#
!# Only transforming / transposing two components of velocity (x,z)
!#
!######################################################################

       subroutine transform(ux,uy,uz,m) 
	use comp
	use timers_rkstep
	use timers_tran
	implicit none
	include 'intvars'
#include "fft_stuff.f"
	real(b8) factor
	complex(b8), allocatable :: buf(:,:),buf1(:,:,:)
 	real(b8)    :: ux(nxpad,zisz,yjsz,nu)
 	complex(b8) :: uz(nzpad,xisz,yjsz,nu)
#ifdef HOMOGENEOUS
 	complex(b8) :: uy(nypad,zjsz*xisz,nu)
#else
 	complex(b8) :: uy(nypad,zjsz,xisz,nu)
#endif
	integer :: m,i, is1,is2,novt,rkstep,i2f,x2,iz,ix,iz2,dnz,iy,l
        real(b8) s1,s2,s3
        complex(b8), allocatable :: bk1i(:,:),bk3i(:,:)

c 8/21/2015: fix for do 75 in ROTD mode which caused incorrect
c results on non-2pi^3 domains
c
	integer ithr,omp_get_thread_num,yp1,yp2

	integer i1,ii
	real(8) rtime1

        if (taskid.eq.0.and.istep.eq.1) write (6,*) 'call transform, kstep=',kstep
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
!factor=1./nx
	factor=1./real(nxpad,b8)

! transform to x-wavenumber space, and take transpose 
!
#ifdef ROTD
	if (kstep.eq.1) then
	novt=6+nc
	else
	novt=5+nc
	end if
#elif defined CVC_PRESS
	if (flagc) then
        novt=6+nc
	else
	novt=5+nc
	end if
#else
	novt=5+nc
#endif

	rtime1=MPI_WTIME()
c
! Transform R2C in X, transpose x <-> Z, truncate in X if needed
        if (iproc.gt.1) then
#ifdef MULTIVAR
#ifdef LVSTEP
	if (mod(jstep-1,ivstep).eq.0)
     1      call xkcomm1_trans(ux,uz,3)
         call xkcomm1_trans(ux(1,1,1,4),uz(1,1,1,4),novt-3)
#else
           call xkcomm1_trans(ux,uz,novt)
#endif
#else

! Watch for array overlap in case of uneven division

#ifdef LVSTEP
	i1=1
	if (mod(jstep-1,ivstep).ne.0) i1=4
           do i=i1,novt
#else
           do i=1,novt
#endif
              call xkcomm1_trans(ux(1,1,1,i),uz(1,1,1,i),1)
           enddo
#endif
        else
           call f_reorder_xz_trans_xz(ux,uz,novt)
        endif

        t_trans(1,2)=t_trans(1,2)+(MPI_WTIME()-rtime1)

	if (iproc .gt. 1) then
!
! transform in z
!
	i1=1
#ifdef LVSTEP
	if (mod(jstep-1,ivstep).ne.0) i1=4
#endif
c
	call divide_thr (yjsz,iyp1,iypi,'transform: yjsz')
c
	rtime1=MPI_WTIME()

	ithr=0

!$OMP PARALLEL private (ithr,yp1,yp2)
c
#ifdef OPENMP
      ithr = OMP_GET_THREAD_NUM()
#endif
#ifdef ESSL
        call cft(1,uz,1,nzpad,uz,1,nzpad,nzpad,xisz*yjsz/num_thr,
     1       1,1.,
     1        caux1(1,ithr),cnaux,caux2(1,ithr),cnaux)
#endif
c
           do 60 i=i1,novt
c
        yp1=iyp1(ithr)

#ifdef FFTW
        call fftw_execute(plan2_xk(ithr),uz(1,1,yp1,i),uz(1,1,yp1,i))
#else
        call escfft (uz(1,1,yp1,i),1,nzpad,nzpad,xisz*yjsz/num_thr,-1,ithr)
#endif
 60           continue

!$OMP END PARALLEL

        t_trans(2,2)=t_trans(2,2)+(MPI_WTIME()-rtime1)
        rtime1=MPI_WTIME()
c
 	endif


!#################################################
! operations for convective terms
!#################################################

#ifdef LVSTEP
	if (mod(jstep-1,ivstep).ne.0) go to 77
#endif
c
        allocate (bk1i(xist:xien,2))
        allocate (bk3i(nz,2))
c
	ithr=0
!$OMP PARALLEL private (ithr,y,x,bk1i,bk3i,yp1,yp2)
c
#ifdef OPENMP
      ithr = OMP_GET_THREAD_NUM()
#endif
c
!,kstep these new arrays introduced by PKY, 4/9/2012


        do x=xist,xien
        bk1i(x,1)=imagi*kx(x)
        bk1i(x,2)=imagi*bk1(x,m)
        end do
c
        do z=1,nz
        bk3i(z,1)=imagi*bk3(z,m)
        bk3i(z,2)=imagi*kz(z)
        end do

	yp1=ithr*yjsz/num_thr+1
	yp2=yp1+yjsz/num_thr-1
c
c       do 70 yp=yp1,yp2
        do 70 yp=iyp1(ithr),iyp1(ithr)+iypi(ithr)-1
           y=yjst+yp-1
           
           do 71 xp=1,xisz
              x=xist+xp-1

#ifdef ROTD
              if (kstep.eq.1) then
                 do 75 z=1,nzpad
             if (z.eq.nzhppad) go to 75
#ifdef ORIG
             uz(z,xp,yp,1)=bk1i(x,1)*uz(z,xp,yp,1)+bk3i(z,1)*uz(z,xp,yp,5+nc)
             uz(z,xp,yp,3)=bk1i(x,2)*uz(z,xp,yp,5+nc)+bk3i(z,2)*uz(z,xp,yp,3)
             uz(z,xp,yp,5+nc)=uz(z,xp,yp,2)
             uz(z,xp,yp,2)=bk1i(x,1)*uz(z,xp,yp,2)+bk3i(z,1)*uz(z,xp,yp,4+nc)
#else
             uz(z,xp,yp,1)=s1*bk1i(x,1)*uz(z,xp,yp,1)+s3*bk3i(z,2)*uz(z,xp,yp,5+nc)
             uz(z,xp,yp,3)=s1*bk1i(x,1)*uz(z,xp,yp,5+nc)+s3*bk3i(z,2)*uz(z,xp,yp,3)
             uz(z,xp,yp,5+nc)=uz(z,xp,yp,2)
             uz(z,xp,yp,2)=s1*bk1i(x,1)*uz(z,xp,yp,2)+s3*bk3i(z,2)*uz(z,xp,yp,4+nc)
#endif
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


!$OMP END PARALLEL
	deallocate (bk1i,bk3i)

        t_trans(3,2)=t_trans(3,2)+(MPI_WTIME()-rtime1)
        rtime1=MPI_WTIME()

 77	continue 
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

      if(jproc .gt. 1) then
#ifdef HOMOGENEOUS
         if(num_al_i(0) .gt. 0) then
#ifdef MULTIVAR
#ifdef LVSTEP
 	if (mod(jstep-1,ivstep).eq.0)
     1      call xkcomm2_trans_cyl (uz,uy,3)
             call xkcomm2_trans_cyl (uz(1,1,1,4),uy(1,1,4),i2f-3)
#else
            call xkcomm2_trans_cyl (uz,uy,i2f)
#endif
#else
            do i=1,i2f
               call xkcomm2_trans_cyl (uz(1,1,1,i),uy(1,1,i),1)
            end do
#endif
         endif

#else
#ifdef MULTIVAR
         call xkcomm2_trans_pen (uz,uy,i2f)
#else
         do i=1,i2f
            call xkcomm2_trans_pen (uz(1,1,1,i),uy(1,1,1,i),1)
         end do
#endif
#endif
      else
! 1D with jproc=1

         call f_reorder_yz_trans_y(uz,uy,i2f)

      endif

        t_trans(4,2)=t_trans(4,2)+(MPI_WTIME()-rtime1)


      return
      end subroutine transform

!-----------------------------------------------------------

      subroutine f_reorder_xz_trans_xz(ux,uz,int)

      use comp
      implicit none
      include 'intvars'

      complex(b8) uz(nzpad,xisz,yjsz,int)
      real(b8) ux(nxpad,zisz,yjsz,int)
      integer i,iz,ix,x2,int
      complex(b8), allocatable :: buf(:,:)
      real(b8) factor

      factor=1./real(nxpad,b8)

      allocate(buf(nxhppad,zisz))
      buf = 0

      do i=1,int
         do yp=1,yjsz
                 
#ifdef FFTW
            call sfftw_execute_dft_r2c(plan1_p2b,ux(1,1,yp,i),buf)
#else
            call rcft(0,ux(1,1,yp,i),nxpad,buf,nxhppad,nxpad,zisz,1,factor,
     &                raux1, rnaux1,raux2,rnaux2,raux3,rnaux3)           
#endif


            do z=1,nzpad,NB1_Z
               z2 = min(z+NB1_Z-1,nzpad)
               
               do x=1,nxh,NB1_X
                  x2 = min(x+NB1_X-1,nxh)
                  
                  do iz=z,z2
                     do ix=x,x2
#ifdef FFTW
                        uz(iz,ix,yp,i) = buf(ix,iz) * factor
#else
                        uz(iz,ix,yp,i) = buf(ix,iz)
#endif
                     enddo
                  enddo
               enddo
            enddo
            
#ifdef FFTW
            call fftw_execute(plan2_p2b,uz(1,1,yp,i),uz(1,1,yp,i))
#else
            call cft(0,uz(1,1,yp,i),1,nzpad,uz(1,1,yp,i),1,nzpad,nzpad,nxh,1,1.0, 
     1           caux1,cnaux1,caux2,cnaux2)                 
#endif
         enddo

      enddo

      deallocate(buf)

      return
      end subroutine

      
      subroutine f_reorder_yz_trans_y(uz,uy,int)

      use comp
      implicit none
      include 'intvars'

      integer i,iz,iy,dnz,int
      complex(b8) uz(nzpad,xisz,yjsz,int)

#ifdef HOMOGENEOUS
      complex(b8) uy(ny,zjsz*xisz,int)


#else
      complex(b8) uy(nypad,zjsz,xisz,int)
      complex(b8), allocatable :: buf(:,:,:)
      
      allocate(buf(nypad,nz,xisz))
      buf = 0

      dnz = nzpad - nz
      do i=1,int
         do x=1,xisz
            
            do z=1,nz,NB2_Z
               z2 = min(z+NB2_Z-1,nz)
               
               do y=1,nypad,NB2_Y
                  y2 = min(y+NB2_Y-1,nypad)
                  
                  do iz=z,min(z2,nzh)
                     do iy=y,y2
                        buf(iy,iz,x) = uz(iz,x,iy,i)
                     enddo
                  enddo
                  iz2 = iz
                  do iz=iz2,z2
                     do iy=y,y2
                        buf(iy,iz,x) = uz(iz+dnz,x,iy,i)
                     enddo
                  enddo
               enddo
            enddo
            
c #ifdef FFTW
c            call fftw_execute(plan2_p3,buf(1,1,x),buf(1,1,x))
c #else
c            call escfft (buf(1,1,x),1,nypad,nypad,zjsz,-1)
c #endif
         enddo
         
         do x=1,xisz
            do z=1,nz
               do y=1,nypad
                  uy(y,z,x,i) = buf(y,z,x)
               enddo
            enddo
         enddo
         
      enddo
      deallocate(buf)
#endif
      
      
      return
      end subroutine
