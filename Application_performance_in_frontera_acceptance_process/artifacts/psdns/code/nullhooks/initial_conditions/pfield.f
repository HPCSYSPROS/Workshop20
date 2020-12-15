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
      
      real(b8), allocatable :: sinx(:),siny(:),cosx(:),cosy(:)
      
      allocate (sinx(nx),stat=ierr)
      allocate (cosx(nx),stat=ierr)
      allocate (siny(ny),stat=ierr)
      allocate (cosy(ny),stat=ierr)

      data au,av/1.,1./

      yg=taskid+1
      if (taskid.eq.0) 
 1    write (6,*) 'enter pfield, taskid=',taskid, nx,zisz,yjsz, Ntot
      
!     print to standard output the initial conditions in use
      if (taskid.eq.0) write (6,*) "using nullhook pfield: sinusoidal"

      aw=-(au+av)
      
      pi=atan(1.)*4.
      dx=2.*pi/nx/beta1
      dy=2.*pi/ny/beta2
      dz=2.*pi/nz/beta3
      
      do x=1,nx
         sinx(x)=sin(beta1*(x-1)*dx)
         cosx(x)=cos(beta1*(x-1)*dx)
      enddo

      do y=1,ny
         siny(y)=sin(beta2*(y-1)*dy)
         cosy(y)=cos(beta2*(y-1)*dy)
      enddo
      
      umax=0.
      
      do yp=1,yjsz
         
         y=yjst+yp-1
         
         do zp=1,zisz
            z=zist+zp-1
            sinz=sin(beta3*(z-1)*dz)
            cosz=cos(beta3*(z-1)*dz)
            
            termu=cosz*cosy(y)
            termv=siny(y)*cosz
            termw=cosy(y)*sinz
    
            do x=1,nx
               unxr(x,zp,yp,1)=au*sinx(x)*termu
               unxr(x,zp,yp,2)=av*cosx(x)*termv
               unxr(x,zp,yp,3)=aw*cosx(x)*termw
               if (abs(unxr(x,zp,yp,3)).gt.umax) then
                  umax=abs(unxr(x,zp,yp,3)) 
                  ixm=x
                  iym=y
                  izm=z
               end if
            enddo !done with x
         enddo  !done with zp
      enddo  !done with yp

      deallocate (sinx,cosx,siny,cosy)

     write (6,*) 'pfield: max.abs.value=',umax,ixm,iym,izm,
 1    'yjsz,zisz=',yjsz,zisz

      return
      end
