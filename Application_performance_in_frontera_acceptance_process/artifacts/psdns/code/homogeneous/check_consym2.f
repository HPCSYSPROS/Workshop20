      subroutine check_consym2 (uy,nv,title,shellk)
      use comsp
      implicit none
      include 'intvars'
      complex(b8) :: uy(ny,xisz*zjsz,nv)
      integer :: nzp2,nyp2,xst,k,ik,lu,nchar,a,shellk
      real :: term,rk2,cmag1,cmag2,cmag3,diff, difmax
      character(*) title
      character*30 caux
      integer,save :: icall
      data icall/0/
      complex(b8),dimension(:,:,:),allocatable :: modes,modesr
      integer :: k2i,k3i,bnd,cnt,iv,nv
      bnd=shellk-1
      allocate(modes(-bnd:bnd,-bnd:bnd,nv))
      allocate(modesr(-bnd:bnd,-bnd:bnd,nv))
      modes(:,:,:)=(0.0_b8,0.0_b8)
      modesr(:,:,:)=(0.0_b8,0.0_b8)
      kx2(:)=kx(:)**2
      ky2(:)=ky(:)**2
      kz2(:)=kz(:)**2
      ky2(nyhp)=0.
      kz2(nzhp)=0.
      caux=title
      call blanks (caux,nchar)

      xp=mystart_x
      z=mystart_z
      do 10 a=1,num_al
         x=xp+xist-1
         do 11 y=1,ny
            if(.not.mask(y,a)) go to 11
            rk2=kx2(x)+ky2(y)+kz2(z)
            ik=sqrt(rk2)+1.5
            if (ik.gt.nxh) go to 11
            if (ik.le.shellk.and.x.eq.1) then
               k2i = int(ky(y))
               k3i = int(kz(z))
               modes(k2i,k3i,:) = uy(y,a,:)
            end if
 11      continue
         call next_xz(xp,z)
 10   continue
c
c Reduce the modes at kx=0 to the root process for writing.
c
      cnt=(2*bnd+1)**2*nv
      call mpi_reduce(modes,modesr,cnt,mpicomplex,
     1                MPI_SUM,0,MPI_COMM_WORLD,mpierr)
      if (taskid.eq.0) then
c
c modes with +k2 & +k3 and -k2 & -k3.
c
	difmax=0.

	do 20 iv=1,nv

         do k3i=0,bnd
            do k2i=0,bnd
               ik = int(sqrt(real(k2i,b8)**2+real(k3i,b8)**2)+1.5)
               if (ik.gt.shellk) cycle
		cmag1=cabs(modesr(k2i,k3i,iv))
		cmag2=cabs(modesr(-k2i,-k3i,iv))
		cmag3=cabs(modesr(k2i,k3i,iv)-conjg(modesr(-k2i,-k3i,iv)))
		diff=cmag3/sqrt(cmag1*cmag2)
c              write(6,222)caux(1:nchar),k2i,k3i,modesr(k2i,k3i,iv)
		if (cmag1*cmag2.gt.0.) then
c              write(6,223)caux(1:nchar),-k2i,-k3i,modesr(-k2i,-k3i,iv),ik,diff
		difmax=max(difmax,diff)
		else
c              write(6,222)caux(1:nchar),-k2i,-k3i,modesr(-k2i,-k3i,iv)
		end if
 222           format('check_consym2: ',a10,2i4,1p,2e16.8)
 223           format('check_consym2: ',a7,2i4,1p,2e16.8,2x,i3,1p,2e12.4)
            end do
         end do
c
c modes with +k2 & -k3 and -k2 & +k3.
c
         do k3i=-1,-bnd,-1
            do k2i=1,bnd
               ik = int(sqrt(real(k2i,b8)**2+real(k3i,b8)**2)+1.5)
               if (ik.gt.shellk) cycle
		cmag1=cabs(modesr(k2i,k3i,iv))
		cmag2=cabs(modesr(-k2i,-k3i,iv))
		cmag3=cabs(modesr(k2i,k3i,iv)-conjg(modesr(-k2i,-k3i,iv)))
		diff=cmag3/sqrt(cmag1*cmag2)
c              write(6,222)caux(1:nchar),k2i,k3i,modesr(k2i,k3i,iv)
		if (cmag1*cmag2.gt.0.) then
c              write(6,223)caux(1:nchar),-k2i,-k3i,modesr(-k2i,-k3i,iv),ik,diff
		difmax=max(difmax,diff)
		else
c              write(6,222)caux(1:nchar),-k2i,-k3i,modesr(-k2i,-k3i,iv)
		end if
            end do
         end do

 20	continue
c
	write (6,*) 'check_consym2: istep, difmax=',istep,title,difmax
      endif
      deallocate(modes)
      deallocate(modesr)
      return
      end
