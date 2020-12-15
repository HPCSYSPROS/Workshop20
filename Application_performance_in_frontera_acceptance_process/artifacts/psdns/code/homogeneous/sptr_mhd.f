      subroutine sptr_mhd (uny,m)
c
#ifdef MHD

      use comsp
      implicit none
      include 'intvars'

        complex(b8) :: uny(ny,zjsz*xisz,3)
        integer :: i,j,k,l,ik,m
        real(b8) :: rk2,term

c
        integer :: iy,nzp2,nyp2,a

        real sum
c
        real(b8) beta_min
c
	integer ithr,omp_get_thread_num,iy1,iy2,ii
c
	integer, allocatable :: ixp(:),iz(:)
        real, allocatable :: jdky_thr(:,:,:)
c
        real factor
        integer nshell

        nshell=max(nxh,nyh,nzh)
c
        allocate (jdky(nshell,3),jdkc(nshell,3),jdk(nshell))
	allocate (ixp(0:num_thr-1),iz(0:num_thr-1))
	allocate (jdky_thr(nshell,3,0:num_thr-1))
c

        kx2(:)=b11(m)*kx(:)**2
        ky2(:)=b22(m)*ky(:)**2
        kz2(:)=b33(m)*kz(:)**2

	beta_min=amin1(beta1,beta2,beta3)
c
	ky2(nyhp)=0.
	kz2(nzhp)=0.
c
	jdky_thr=0.
	jdky=0.
c
	ithr=0
c
!$OMP PARALLEL private (ithr,a,x,tfact_x,y,rk2,ik,term,xp,z,iy1,iy2,factor) shared (mystart_x,mystart_z)

#ifdef OPENMP
      ithr = OMP_GET_THREAD_NUM()
#endif

      xp = mystart_x
      z = mystart_z
c
      do a=1,num_al
                                                                                
         x = xp + xist-1
                                                                                
         if(x .eq. 1) then
            tfact_x=1
         else
            tfact_x=2
         endif
                                                                                
	iy1=ithr*ny/num_thr+1
	iy2=(ithr+1)*ny/num_thr
	if (x.eq.1.and.z.eq.1.and.ithr.eq.0) iy1=2

         do 10 y=iy1,iy2
c
	if (.not.mask(y,a)) go to 10
c
	rk2=kx2(x)+ky2(y)+kz2(z)
	ik=sqrt(rk2)/beta_min+1.5

	if (ik.gt.nshell) go to 10
c
c fix spotted by X.M. Zhai on 1/3/2016
c to accommodate different choices of imaxis

        if (imaxis.eq.1) then
        factor=kx2(x)/rk2
        else if (imaxis.eq.2) then
        factor=ky2(y)/rk2
        else if (imaxis.eq.3) then
        factor=kz2(z)/rk2
        end if
        
        
        do i=1,3
        term=tfact_x*uny(y,a,i)*conjg(uny(y,a,i))
        jdky_thr(ik,i,ithr)=jdky_thr(ik,i,ithr)+term*factor
        end do
c
c
 10	continue
c
#ifdef OPENMP
	ixp(ithr)=xp
	iz(ithr)=z
  	call next_xz(ixp(ithr),iz(ithr))
	xp=ixp(ithr)
	z=iz(ithr)
#else
  	call next_xz(xp,z)
#endif

	end do
c
!$OMP BARRIER
!$OMP MASTER
c Bug fixed by X.M. Zhai, 1/3/2016:
c "ithr" should be "ii" in the loop below
	do ii=0,num_thr-1
	do i=1,3
	do ik=1,nshell
c	jdky(ik,i)=jdky(ik,i)+jdky_thr(ik,i,ithr)
	jdky(ik,i)=jdky(ik,i)+jdky_thr(ik,i,ii)
	end do
	end do
	end do
!$OMP END MASTER
!
!$OMP END PARALLEL

      call MPI_REDUCE (jdky,jdkc,nshell*3,mpireal,MPI_SUM,0,
     1                 MPI_COMM_WORLD,ierr)

        if (taskid.eq.0) then
        !jdk(:)=jdkc(:,1)+jdkc(:,2)+jdkc(:,3)
! ZHAI: 09/18/2016 consider different dirs + correct metric factors
        joule_diss(1:3)=0.
        do ik=1,nshell
        joule_diss(1)=joule_diss(1)+b11(m)*jdkc(ik,1)
        joule_diss(2)=joule_diss(2)+b22(m)*jdkc(ik,2)
        joule_diss(3)=joule_diss(3)+b33(m)*jdkc(ik,3)
        end do
        sum=joule_diss(1)+joule_diss(2)+joule_diss(3)
        joule_diss(:)=conduc*joule_diss(:)
	write (6,*) 'conduc,sum,joule_diss=',conduc,sum,conduc*sum
        end if

        deallocate (jdky,jdkc,jdk)
	deallocate (ixp,iz,jdky_thr)
c
#endif
	return
	end
