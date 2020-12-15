          subroutine dissenst (uny,unx,uy,ux)
c
            use comp
            implicit none
            include 'intvars'

        complex(b8) :: uny(nypad,zjsz*xisz,3)
        complex(b8) :: uy(nypad,zjsz*xisz,5)
        real(b8) :: unx(nx,zisz,yjsz,3)
        real(b8) :: ux(nx,zisz,yjsz,5)
        
            integer :: i,j,m
            complex(b8), allocatable :: shift(:)
            integer iy,a

        integer ithr,omp_get_thread_num,ia
	complex bk1i,bk3i
        real s11,s22,s33,s12,s13,s23,om1,om2,om3,ave1,ave2
c
        integer nhl
        real, allocatable :: binl(:),hist0(:)
        character*40 title
        real epsmax,epsmax0,enstmax,enstmax0

#ifdef ROTD
	real(b8) vortm(4,3),omskew,omflat
	real(b8), allocatable :: vort(:,:,:)
	allocate (vort(nx,zisz,yjsz))
#endif

	if (taskid.eq.0) write (6,*) 'enter dissenst'
        m=2
        uy(:,:,1:4)=0.
c
!$OMP PARALLEL private (ithr,x,xp,z,bk1i,bk2i,bk3i)

#ifdef OPENMP
      ithr = OMP_GET_THREAD_NUM()
#else
        ithr=0
#endif
c
      xp = ia_xst(ithr)
      z = ia_zst(ithr)
c
      do 10 a = ia_st(ithr), ia_st(ithr)+ia_sz(ithr)-1
c
         x = xp + xist-1
	bk1i=imagi*b11(m)*kx(x)
	bk3i=imagi*b33(m)*kz(z)
c
c differentiate in wavenumber space
c operate on already phase-shifted velocities)
!
c du/dx, du/dy, dv/dx, dw/dz
c
            do 15 y=1,nypad
               if(mask(y,a)) then
                  bk2i=imagi*b22(m)*ky(y)
 		uy(y,a,1)=bk1i*uny(y,a,1)
 		uy(y,a,2)=bk2i*uny(y,a,1)
 		uy(y,a,3)=bk1i*uny(y,a,2)
 		uy(y,a,4)=bk3i*uny(y,a,3)
		end if
 15         continue

            call next_xz(xp,z)
 10	continue

!$OMP END PARALLEL
!
        call kxtran (uy,ux,ux,4)
c
        do 20 yp=1,yjsz
        do 20 zp=1,zisz
        do 20 x=1,nx
        s11=ux(x,zp,yp,1)
        s12=.5*(ux(x,zp,yp,2)+ux(x,zp,yp,3))
        s33=ux(x,zp,yp,4)
        om3=(ux(x,zp,yp,3)-ux(x,zp,yp,2))
        s22=-s11-s33
        ux(x,zp,yp,4)=s11*s11+2.*s12*s12+s22*s22+s33*s33
        ux(x,zp,yp,5)=om3*om3
#ifdef ROTD
	vort(x,zp,yp)=om3
#endif
 20     continue
c
#ifdef ROTD
	call phymom4 (vort,1,1,1,vortm(1,3))
#endif

        uy(:,:,1)=0.
        uy(:,:,2)=0.
c
c
!$OMP PARALLEL private (ithr,x,xp,z,bk1i,bk2i,bk3i)

#ifdef OPENMP
      ithr = OMP_GET_THREAD_NUM()
#else
        ithr=0
#endif

      xp = ia_xst(ithr)
      z = ia_zst(ithr)
c
      do 30 a = ia_st(ithr), ia_st(ithr)+ia_sz(ithr)-1
c
         x = xp + xist-1
	bk1i=imagi*b11(m)*kx(x)
	bk3i=imagi*b33(m)*kz(z)
c
c differentiate in wavenumber space
c operate on already phase-shifted velocities)
!
c du/dz, dw/dx
c
            do y=1,nypad
               if(mask(y,a)) then
 		uy(y,a,1)=bk3i*uny(y,a,1)
 		uy(y,a,2)=bk1i*uny(y,a,3)
		end if
            end do

            call next_xz(xp,z)
 30	continue

!$OMP END PARALLEL
!
        call kxtran (uy,ux,ux,2)
c
        do 40 yp=1,yjsz
        do 40 zp=1,zisz
        do 40 x=1,nx
        s13=.5*(ux(x,zp,yp,1)+ux(x,zp,yp,2))
        om2=(ux(x,zp,yp,1)-ux(x,zp,yp,2))
        ux(x,zp,yp,4)=ux(x,zp,yp,4)+2.*s13*s13
        ux(x,zp,yp,5)=ux(x,zp,yp,5)+om2*om2
#ifdef ROTD
	vort(x,zp,yp)=om2
#endif
 40     continue
#ifdef ROTD
	call phymom4 (vort,1,1,1,vortm(1,2))
#endif
c
        uy(:,:,1)=0.
        uy(:,:,2)=0.
c
c
!$OMP PARALLEL private (ithr,x,xp,z,bk1i,bk2i,bk3i)

#ifdef OPENMP
      ithr = OMP_GET_THREAD_NUM()
#else
        ithr=0
#endif

      xp = ia_xst(ithr)
      z = ia_zst(ithr)
c
      do 50 a = ia_st(ithr), ia_st(ithr)+ia_sz(ithr)-1
c
         x = xp + xist-1
	bk3i=imagi*b33(m)*kz(z)
c
c differentiate in wavenumber space
c operate on already phase-shifted velocities)
!
c dv/dz, dw/dy
c
            do y=1,nypad
               if(mask(y,a)) then
                  bk2i=imagi*b22(m)*ky(y)
 		uy(y,a,1)=bk2i*uny(y,a,3)
 		uy(y,a,2)=bk3i*uny(y,a,2)
		end if
            end do

            call next_xz(xp,z)
 50	continue

!$OMP END PARALLEL
!
        call kxtran (uy,ux,ux,2)
c
        do 60 yp=1,yjsz
        do 60 zp=1,zisz
        do 60 x=1,nx
        s23=.5*(ux(x,zp,yp,1)+ux(x,zp,yp,2))
        om1=(ux(x,zp,yp,1)-ux(x,zp,yp,2))
        ux(x,zp,yp,4)=ux(x,zp,yp,4)+2.*s23*s23
        ux(x,zp,yp,5)=ux(x,zp,yp,5)+om1*om1
#ifdef ROTD
	vort(x,zp,yp)=om1
#endif
 60     continue
c
#ifdef ROTD
	call phymom4 (vort,1,1,1,vortm(1,1))
	if (taskid.eq.0) then
	if (istep.eq.0) then
	write (753,"(' istep    time         skew(omx)   flat(omx)   skew(omy)   flat(omy)   skew(omz)   flat(omz)')") 
	end if
	vortm(3,:)=vortm(3,:)/vortm(2,:)**1.5
	vortm(4,:)=vortm(4,:)/vortm(2,:)**2
	write (753,"(i6,1p,e12.4,2x,1p,6e12.4)") istep,time,(vortm(3,i),vortm(4,i),i=1,3)
	end if
#endif

        uy(:,:,1)=uny(:,:,1)
        uy(:,:,2)=uny(:,:,2)
        uy(:,:,3)=uny(:,:,3)
c
        call phyave (ux(1,1,1,4),ave1)
        call phyave (ux(1,1,1,5),ave2)
c
        ave1=2.*viscos*ave1
        ave2=viscos*ave2
        if (taskid.eq.0) write (6,*) 'dissenst: ave1,ave2=',ave1,ave2

        ave1=ave1/2./viscos
        ave2=ave2/viscos
c
        nhl=201
        allocate (binl(nhl),hist0(nhl))
        call logspace (-7.,7.,nhl,binl,1)
        write (title,"('PDF of eps/<eps>, istep=',i6)") istep
        call phypdfkrsp (ux(1,1,1,4),751,title,
     1                    0.,ave1,hist0,nhl,binl)
        write (title,"('PDF of enst/<enst>, istep=',i6)") istep
        call phypdfkrsp (ux(1,1,1,5),751,title,
     1                    0.,ave2,hist0,nhl,binl)
c
        epsmax=0.
        enstmax=0.
        do 70 yp=1,yjsz
        do 70 zp=1,zisz
        do 70 x=1,nx
        epsmax=max(epsmax,ux(x,zp,yp,4))
        enstmax=max(enstmax,ux(x,zp,yp,5))
 70     continue
        epsmax=epsmax/ave1
        enstmax=enstmax/ave2
        call MPI_REDUCE (epsmax,epsmax0,1,mpireal,MPI_MAX,0,
     1                   MPI_COMM_WORLD,mpierr)
        call MPI_REDUCE (enstmax,enstmax0,1,mpireal,MPI_MAX,0,
     1                   MPI_COMM_WORLD,mpierr)
        if (taskid.eq.0) then
        write (6,"('max eps/<eps>,enst/<enst>=',1p,2e11.4)") epsmax0,enstmax0
        write (752,"('istep,time=',i6,1p,e12.4,'  max eps/<eps>,enst/<enst>=',1p,2e11.4)") istep,time,epsmax0,enstmax0
        end if
        
        deallocate (binl,hist0)
        
#ifdef ROTD
	deallocate (vort)
#endif
        
        return 
      end
