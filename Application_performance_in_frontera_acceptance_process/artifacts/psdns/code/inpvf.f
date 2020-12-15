!inputs are un, u
c
c routine for using a user-supplied definition to prescribe the
c initial velocity field in physical space, with subsequent
c transformation to wavenumber space (to fit with the structure
c of the main code).
c
!NPM converting to stride-1 3/16/09 
      subroutine inpvf(unxr,uxr)
	use comp
	implicit none
	include 'intvars'

#ifndef UNEV
      real(b8) :: unxr(nxpad,zisz,yjsz,3)
      real(b8) :: uxr(nxpad,zisz,yjsz,3)
#else	
      real(b8) :: unxr(nxpad,zisz,yjsz+padx,3)
      real(b8) :: uxr(nxpad,zisz,yjsz+padx,3)
#endif
	integer ii
c
c the specification is contained in a user-supplied
c subroutine, which should be called "pfield"
c (continuity must be satisfied)
c
!$OMP DO
	do yp=1,yjsz
           do zp=1,zisz
              do x=1,nxpad
                 unxr(x,zp,yp,:)=0.
              enddo
           enddo
	enddo
!$OMP END DO

        call pfield(unxr)
c
c transform to wavenumber space
c
#ifdef MULTIVAR
	call xktran (unxr,unxr,unxr,3)
#else
      do ii=1,3
         call xktran (unxr(1,1,1,ii),unxr(1,1,1,ii),unxr(1,1,1,ii),1)
      end do
#endif
c
c convert to u/beta1, v/beta2, w/beta3
c
c	unxr(:,:,:,1)=unxr(:,:,:,1)/beta1
c	unxr(:,:,:,2)=unxr(:,:,:,2)/beta2
c	unxr(:,:,:,3)=unxr(:,:,:,3)/beta3
c
      if (luinit(1).eq.0.and.xist.eq.1.and.zjst.eq.1) then
         write (6,*) 'mean velocities forced to zero in sub. inpvf'
         unxr(1,1,1,1)=0.
         unxr(1,1,1,2)=0.
         unxr(1,1,1,3)=0.
      end if


#ifdef BOUSS
	call noise3d (unxr,uxr,2)
	call enfsym (unxr)
c	call check_velmax (unxr,'exit of inpvf')
#endif

       if (taskid.eq.0) write (6,*) 'exit of inpvf'

c	write (6,*) 'taskid,unxr(2,2,2,2)=',taskid,unxr(2,2,2,2)
      return
      end
