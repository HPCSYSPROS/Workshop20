c     
c     Taylor-Green vortex, following Riley & deBruyn Kops (PF Jul 2003)
c     
      subroutine pfield (unxr)
      use com
      implicit none
      include 'intvars'

c     set initial velocity field in physical space
c     
c     u =  U cos(kappa.z)cos(kappa.x)sin(kappa.y)
c     v = -U cos(kappa.z)sin(kappa.x)cos(kappa.y)
c     w =  0
c     
#ifndef UNEV
      real(b8) :: unxr(nx,zisz,yjsz,3)
#else	
      real(b8) :: unxr(nx,zisz,yjsz+padx,3)
#endif
      
      real(b8) kappa, uscale, dx, dy, dz, twopi
      real(b8) sum,sum0,usqave
      real(b8), allocatable :: sinx(:),siny(:),cosx(:),cosy(:),cosz(:)
      
      allocate (sinx(nx),stat=ierr)
      allocate (cosx(nx),stat=ierr)
      allocate (siny(ny),stat=ierr)
      allocate (cosy(ny),stat=ierr)
      allocate (cosz(nz),stat=ierr)
      
      uscale=1.
      kappa=1.

!     print to standard output the initial conditions in use
      if (taskid.eq.0) then
         write (6,*) 'enter pfield, taskid=',taskid, nx,zisz,yjsz, Ntot         
         write (6,*) "initial conditions: using hook taylorgreen"
      endif

      pi=atan(1.)*4.
      twopi=2.*pi
      dx=twopi/beta1/nx
      dy=twopi/beta2/ny
      dz=twopi/beta3/nz
      
      do x=1,nx
         sinx(x)=sin(kappa*(x-1)*dx)
         cosx(x)=cos(kappa*(x-1)*dx)
      end do

      do y=1,ny
         siny(y)=sin(kappa*(y-1)*dy)
         cosy(y)=cos(kappa*(y-1)*dy)
      end do
      
      do z=1,nz
         cosz(z)=cos(2.*kappa*(z-1)*dz)
      end do
      
      sum=0.

      do yp=1,yjsz        
         y=yjst+yp-1         
         do zp=1,zisz
            z=zist+zp-1            
            do x=1,nx
               unxr(x,zp,yp,1)=uscale*cosz(z)*cosx(x)*siny(y)
               unxr(x,zp,yp,2)=-uscale*cosz(z)*sinx(x)*cosy(y)
               sum=sum+unxr(x,zp,yp,1)**2
               unxr(x,zp,yp,3)=0.
            enddo
         enddo
      enddo

      call MPI_ALLREDUCE(sum,sum0,1,mpireal,MPI_SUM,MPI_COMM_WORLD,mpierr)

      if (taskid.eq.0) then
         write (6,*) 'taylorgreen: usqave=',sum0/nx/ny/nz
      end if

c     
c     note Rogallo's algorithm stores u/b11, etc. 
c     
      do yp=1,yjsz
         do zp=1,zisz
            do x=1,nx
               unxr(x,zp,yp,1)=unxr(x,zp,yp,1)/beta1
               unxr(x,zp,yp,2)=unxr(x,zp,yp,2)/beta2
            end do
         end do
      end do
c
c 	call check_velmax (unxr,'taylorgreen')
      
      deallocate (sinx,cosx,siny,cosy,cosz)

      return
      end
