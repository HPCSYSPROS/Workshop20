!################################################
!#  module stratification module
!#
!#  This null hook module
!#  contains all empty subroutine calls
!#
!################################################

      module stratification_module
        use comsp
       

!########################
      contains
!########################
        subroutine stratification(uy,uny,int)
          implicit none
          complex(b8) :: uy(xisz,zjsz,ny,nu)
          complex(b8) :: uny(xisz,zjsz,ny,3+nc)
          integer     :: int

          return
        end subroutine stratification
!########################
      end module stratification_module
