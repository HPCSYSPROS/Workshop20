!#
!#       Null Hook Call
!#

	subroutine hsptr (uy,m)
!#
!# calculation of "horizontal spectrum" (in 2D kx-ky plane)
!#
       use comsp
       implicit none
       include 'intvars'

       complex(b8) :: uy(xisz,zjsz,ny,3)
       integer :: m

       return
       end subroutine hsptr
