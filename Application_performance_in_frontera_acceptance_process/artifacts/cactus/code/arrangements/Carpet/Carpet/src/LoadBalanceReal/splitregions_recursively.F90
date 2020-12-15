#include "cctk.h"

subroutine splitregions_recursively ( &
     cxx_superregs, nsuperregs, &
     cxx_regs, &
     nprocs, &
     ghostsize_, alpha_, limit_size_, &
     procid_)
  use carpet_boxtypes
  implicit none
  
  integer,      intent(in) :: nsuperregs
  CCTK_POINTER, intent(in) :: cxx_superregs
  CCTK_POINTER, intent(in) :: cxx_regs
  integer,      intent(in) :: nprocs
  integer,      intent(in) :: ghostsize_
  CCTK_REAL,    intent(in) :: alpha_
  integer,      intent(in) :: limit_size_
  integer,      intent(in) :: procid_
  
  type(ptr), allocatable :: sregions(:)
  type(boundary)         :: outbound
  
  CCTK_POINTER :: cxx_superreg
  CCTK_POINTER :: cxx_tree
  type(bbox)   :: box
  integer      :: i
  
  
  
  ! Callback functions that are implemented in C++
  interface
     
     subroutine carpet_get_region (cxx_superregs, i, cxx_superreg)
       use carpet_boxtypes
       implicit none
       CCTK_POINTER, intent(in)  :: cxx_superregs
       integer,      intent(in)  :: i
       CCTK_POINTER, intent(out) :: cxx_superreg
     end subroutine carpet_get_region
     
     subroutine carpet_get_bbox (cxx_superreg, box, obound)
       use carpet_boxtypes
       implicit none
       CCTK_POINTER,   intent(in)  :: cxx_superreg
       type(bbox),     intent(out) :: box
       type(boundary), intent(out) :: obound
     end subroutine carpet_get_bbox
     
     subroutine carpet_insert_region (cxx_regs, reg)
       use carpet_boxtypes
       implicit none
       CCTK_POINTER,           intent(in) :: cxx_regs
       type(superregion2slim), intent(in) :: reg
     end subroutine carpet_insert_region
     
     subroutine carpet_create_tree_branch &
          (nch, dir, bounds, cxx_subtrees, cxx_tree)
       use carpet_boxtypes
       implicit none
       integer,      intent(in) :: nch
       integer,      intent(in) :: dir
       integer,      intent(in) :: bounds(nch+1)
       CCTK_POINTER, intent(in) :: cxx_subtrees(nch)
       CCTK_POINTER, intent(in) :: cxx_tree
     end subroutine carpet_create_tree_branch
     
     subroutine carpet_create_tree_leaf (sreg, cxx_tree)
       use carpet_boxtypes
       implicit none
       type(superregion2slim), intent(in)  :: sreg
       CCTK_POINTER,           intent(out) :: cxx_tree
     end subroutine carpet_create_tree_leaf
     
  end interface
  
  
  
  ! Set global parameters
  ghostsize  = ghostsize_
  alpha      = alpha_
  limit_size = limit_size_ /= 0
  procid     = procid_
  
  
  
  allocate (sregions(nsuperregs))
  do i=1, nsuperregs
     call carpet_get_region (cxx_superregs, i-1, cxx_superreg)
     call carpet_get_bbox (cxx_superreg, box, outbound)
     call create_sregion (box, outbound, i-1, sregions(i)%point)
  end do
  
  call SplitSuperRegions (sregions, nprocs)
  
  do i=1, nsuperregs
     call carpet_get_region (cxx_superregs, i-1, cxx_superreg)
     call insert_region (sregions(i)%point, cxx_tree, cxx_regs)
     call carpet_set_tree (cxx_superreg, cxx_tree)
  end do
  
  do i=1, nsuperregs
     call destroy_sregion (sregions(i)%point)
  end do
  deallocate (sregions)
  
contains
  
  recursive subroutine insert_region (sreg, cxx_tree, cxx_regs)
!   The intent has been removed to make it compile with gfortran 4.1.
!    type(superregion2), pointer, intent(in)  :: sreg
    type(superregion2), pointer  :: sreg
    CCTK_POINTER,                intent(in)  :: cxx_regs
    CCTK_POINTER,                intent(out) :: cxx_tree
    
    integer                   :: nch, ich
    integer                   :: dir
    integer                   :: mydir
    integer,      allocatable :: bounds(:)
    CCTK_POINTER, allocatable :: cxx_subtrees(:)
    type(superregion2slim)    :: sregslim
    
    ! TODO: insert tree dependencies into superregs
    
    if (associated(sreg%children)) then
       ! The region has children: traverse them recursively
       nch = size(sreg%children)
       allocate (bounds(nch+1))
       allocate (cxx_subtrees(nch))
       if (nch /= 2) then
          call CCTK_WARN (CCTK_WARN_ABORT, "number of children is not 2")
       end if
       mydir = -1
       do dir=1, 3
          if (sreg%children(1)%point%extent%dim(dir)%upper + &
               sreg%children(2)%point%extent%dim(dir)%stride == &
               sreg%children(2)%point%extent%dim(dir)%lower) then
             ! Found direction
             if (mydir > 0) then
                call CCTK_WARN (CCTK_WARN_ABORT, "could not determine direction")
             end if
             mydir = dir
          else
             if (sreg%children(1)%point%extent%dim(dir)%lower /= &
                  sreg%children(2)%point%extent%dim(dir)%lower .or. &
                  sreg%children(1)%point%extent%dim(dir)%upper /= &
                  sreg%children(2)%point%extent%dim(dir)%upper) then
                call CCTK_WARN (CCTK_WARN_ABORT, "children differ in unexpected ways")
             end if
          end if
       end do
       if (mydir < 0) then
          call CCTK_WARN (CCTK_WARN_ABORT, "could not determine direction")
       end if
       bounds(1) = sreg%children(1)%point%extent%dim(mydir)%lower
       bounds(2) = sreg%children(2)%point%extent%dim(mydir)%lower
       bounds(3) = sreg%children(2)%point%extent%dim(mydir)%upper + &
            sreg%children(2)%point%extent%dim(mydir)%stride
       do ich=1, nch
          call insert_region &
               (sreg%children(ich)%point, cxx_subtrees(ich), cxx_regs)
       end do
       call carpet_create_tree_branch &
            (nch, mydir-1, bounds, cxx_subtrees, cxx_tree)
    else
       ! The region is a leaf: insert it
       sregslim%extent           = sreg%extent
       sregslim%outer_boundaries = sreg%outer_boundaries
       sregslim%map              = sreg%map
       sregslim%processor        = sreg%processor
       call carpet_create_tree_leaf (sregslim, cxx_tree)
       call carpet_insert_region (cxx_regs, sregslim)
    end if
  end subroutine insert_region
  
end subroutine splitregions_recursively
