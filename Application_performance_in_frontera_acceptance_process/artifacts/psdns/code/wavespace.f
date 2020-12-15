      subroutine wavespace(uy,m)
        !If one desired, this routine could be combined with 
        !yitransform
        !to be wavespace in wavespace_module
        use wavespace_module
        use comp
        implicit none
	include 'intvars'
        
#include "fft_stuff.f"
        
	integer :: i,j,m, is1, i2f
        complex(b8) :: uy(xisz,zjsz,ny,nu)
        complex(b8) :: uny(xisz,zjsz,ny,3+nu)
        logical flag

!     complex(b8) :: uy(ny,zjsz,xisz,nu)
!     complex(b8) :: uny(ny,zjsz,xisz,3+nc)                  

        !This is calling a routine from the wavespace_module
c        call wavespace_sub(uy,uny,m,flag)
        call wavespace_sub(uy,m)
        !this call is necessary because fortran 90 modules use 
        !an explicit interfact, which causes a datatype mismatch

        !For inhomogeneous, this will be an empty subroutine call
      
      return
      end
