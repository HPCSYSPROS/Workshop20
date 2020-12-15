      subroutine yitransform(uy,m)
!to be yitransform in wavespace
        use wavespace_module
        use comp
        implicit none

        complex(b8) :: uy(ny,zjsz,xisz,nu)
	integer :: i,j,m

        !This is calling a routine from the wavespace_module
        call yitransform_sub(uy,2)
        !this call is necessary because fortran 90 modules use 
        !an explicit interfact, which causes a datatype mismatch

        !For inhomogeneous, this will be an empty subroutine call

      return
      end
