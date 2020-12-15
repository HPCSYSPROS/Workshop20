!#########################################
!#
!# pfield routine added here to link 
!# should the user not specify any initial
!# conditions (i.e. no pfield present in 
!# initial conditions hook folder
!#
!#   INPUT : field to be initialized
!#   OUTPUT: initialized field
!# 
!#########################################
      subroutine pfield (unxr)
      use com
      include 'intvars'
c     
c     set initial velocity field in physical space
c     
c     sinusoidal velocity field (products of sines and cosines)
c     
c     u = au * sin(x)*cos(y)*cos(z)
c     v = av * cos(x)*sin(y)*cos(z)
c     w = aw * cos(x)*cos(y)*sin(z)
c     
c     continuity requires au + av + aw = 0.
c     
#ifndef UNEV
      real(b8) :: unxr(nx,zisz,yjsz,3)
#else	
      real(b8) :: unxr(nx,zisz,yjsz+padx,3)
#endif

        integer ithr,OMP_GET_THREAD_NUM,yp1,yp2
     
      real(b8), allocatable :: sinx(:),siny(:),cosx(:),cosy(:)
c
	real(b8) sum,sum0
      
      allocate (sinx(nx),stat=ierr)
      allocate (cosx(nx),stat=ierr)
      allocate (siny(ny),stat=ierr)
      allocate (cosy(ny),stat=ierr)

      data au,av/1.,1./

      yg=taskid+1

!     print to standard output the initial conditions in use
      if (taskid.eq.0) then
         write (6,*) 'enter pfield, taskid=',taskid, nx,zisz,yjsz, Ntot         
         write (6,*) "using hook pfield: sinusoidal"
      endif
      
      aw=-(au+av)

      pi=atan(1.)*4.
      dx=2.*pi/nx
      dy=2.*pi/ny
      dz=2.*pi/nz

      do x=1,nx
         sinx(x)=sin((x-1)*dx)
         cosx(x)=cos((x-1)*dx)
      enddo
      do y=1,ny
         siny(y)=sin((y-1)*dy)
         cosy(y)=cos((y-1)*dy)
      enddo

      umax=0.

	ithr=0

!$OMP PARALLEL private(ithr,y,z,termu,termv,termw,sinz,cosz,yp1,yp2)

#ifdef OPENMP
        ithr=OMP_GET_THREAD_NUM()
#endif

        yp1=ithr*yjsz/num_thr+1
        yp2=(ithr+1)*yjsz/num_thr
      do yp=yp1,yp2


         y=yjst+yp-1

         do zp=1,zisz
            z=zist+zp-1
            sinz=sin((z-1)*dz)
            cosz=cos((z-1)*dz)
     
            termu=cosz*cosy(y)
            termv=siny(y)*cosz
            termw=cosy(y)*sinz
            
            do x=1,nx
               unxr(x,zp,yp,1)=au*sinx(x)*termu
               unxr(x,zp,yp,2)=av*cosx(x)*termv
               unxr(x,zp,yp,3)=aw*cosx(x)*termw
            enddo               !done with x
         enddo                  !done with zp
      enddo                     !done with yp
!$OMP END PARALLEL

	sum=0.
	do yp=1,yjsz
	do zp=1,zisz
	do x=1,nx
	sum=sum+unxr(x,zp,yp,1)**2
	end do
	end do
	end do

	call MPI_REDUCE (sum,sum0,1,mpireal,MPI_SUM,
     1                   0,MPI_COMM_WORLD,mpierr)
	if (taskid.eq.0) then
	sum0=sum0/nx/ny/nz
	write (6,*) 'pfield: rms x-vel=',sqrt(sum0)
	end if
c
      deallocate (sinx,cosx,siny,cosy)
	if (taskid.eq.0) write (6,*) 'exit pfield'
      

      return
      end
