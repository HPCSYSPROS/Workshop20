module carpet_region 
  
  implicit none

  integer, parameter :: wp = selected_real_kind(12,99)
  
  ! These empty arrays are used to initialize variables to either the
  ! min or max possible number of kind wp or integer.
  real(wp), dimension(2:1) :: empty
  integer, dimension(2:1) :: iempty
  
  
  
  ! Note that using intent on pointer arguments requires fortran 2003.
  ! It is not allowed in fortran 90 or 95.
  
  ! The basic range structure with lower and upper bounds and the stride
  type range
     integer :: lower, upper, stride
  end type range
  
  ! The basic box structure: a length 3 vector of ranges.
  type bbox
     type(range), dimension(3) :: dim
  end type bbox
  
  ! The outer boundary information structure.
  type boundary
     integer, dimension(3,2) :: obound
  end type boundary
  
! A super region structure similar to the Carpet tree structure.
  type ptr
    type(superregion2), pointer :: point
  end type ptr
  type superregion2
    type(bbox) :: extent
    type(boundary) :: outer_boundaries
    integer :: map
    integer :: processor
    real(wp) :: frac
    type(ptr), pointer :: children(:)
  end type superregion2
  type superregion2slim
    type(bbox) :: extent
    type(boundary) :: outer_boundaries
    integer :: map
    integer :: processor
  end type superregion2slim
 
  contains

!   Routine to allocate memory for a super region. Input is the bbox
!   structure, outer boundary information structure and map number.
!   The processor id is initialized to -1 and the pointer to the children
!   is nullified.
    subroutine create_sregion ( box, outerbound, map, sregion )
      type(bbox), intent(in) :: box
      type(boundary), intent(in) :: outerbound
      integer, intent(in) :: map
!     The intent has been removed to make it compile with gfortran 4.1.
!      type(superregion2), pointer, intent(out) :: sregion
      type(superregion2), pointer :: sregion
!      gfortran dies with an internal compiler error on this initialization.
!      type(superregion2), pointer, intent(out) :: sregion => null()

      allocate ( sregion )
      sregion%extent = box
      sregion%outer_boundaries = outerbound
      sregion%map = map
      sregion%processor = -1
      sregion%frac = 0.0_wp
      nullify(sregion%children)

    end subroutine create_sregion
 

!   Routine to nullify the pointers of the children. This has to be used
!   before destroying a temporary super region and assumes another super
!   region has children with pointers to these super regions as well.
!   The use of intent(in) on the sregion just ensures that sregion itself
!   is not modified. Its childrens pointer can be nullified.
    subroutine disassociate ( sregion )
!     The intent has been removed to make it compile with gfortran 4.1.
!      type(superregion2), pointer, intent(in) :: sregion
      type(superregion2), pointer :: sregion
      integer :: n, i

!     Only do something if the children are associated.
      if ( associated(sregion%children) ) then
        n = size(sregion%children)
        do i = 1, n
          nullify ( sregion%children(i)%point )
        end do
      end if
    end subroutine disassociate


!   Routine to assign the pointers of the children of sregion1 to the same
!   super regions as sregion2.
    subroutine point_to_children ( sregion1, sregion2 )
!     The intent has been removed to make it compile with gfortran 4.1.
!      type(superregion2), pointer, intent(inout) :: sregion1
!      type(superregion2), pointer, intent(in) :: sregion2
      type(superregion2), pointer :: sregion1
      type(superregion2), pointer :: sregion2
      integer :: n1, n2, i

!     Only do something if sregion2 actually has children.
      if ( associated(sregion2%children) ) then
        n2 = size(sregion2%children)

!       If sregion1 has children as well make sure to deallocate their
!       sub tree (this probably shouldn't happen.
        if ( associated(sregion1%children) ) then
          n1 = size(sregion1%children)
          do i = 1, n1
            call destroy_sregion ( sregion1%children(i)%point )
          end do
        end if

!       Allocate the children and assign the pointers.
        allocate ( sregion1%children(n2) )
        do i = 1, n2
          sregion1%children(i)%point => sregion2%children(i)%point
        end do

!     Otherwise print a warning. This should not happen.
      else
        print*, 'Warning sregion2 has no children'
      end if  
    end subroutine point_to_children
      

!   Routine to recursively deallocate memory for a super region.
    recursive subroutine destroy_sregion ( sregion )
!     The intent has been removed to make it compile with gfortran 4.1.
!      type(superregion2), pointer, intent(inout) :: sregion
      type(superregion2), pointer :: sregion
      integer :: n, i

!     If the super region has children.
      if ( associated(sregion%children) ) then

!       Find out how many children.
        n = size(sregion%children)

!       Loop over all children.
        do i = 1, n

!         If the children has children call recursively.
          if ( associated ( sregion%children(i)%point ) ) then
            call destroy_sregion ( sregion%children(i)%point )
          end if
        end do

!       Then deallocate the storage for the children pointers.
        if ( n > 0 ) then
          deallocate ( sregion%children )
        end if
      end if

!     Finally deallocate storage for the super region itself.
      deallocate (sregion)

    end subroutine destroy_sregion


!   Routine to recursively print all information about a super region tree
!   structure to stdout. 
    recursive subroutine print_tree ( sregion )
!     The intent has been removed to make it compile with gfortran 4.1.
!      type(superregion2), pointer, intent(in) :: sregion
      type(superregion2), pointer :: sregion
      integer :: n, i
      
!     Only do something if the super region pointer is associated with a
!     valid target.
      if ( associated(sregion) ) then
        print*, 'bbox = ', sregion%extent
        print*, 'outer_boundaries = ', sregion%outer_boundaries
        print*, 'map = ', sregion%map
        print*, 'processor = ', sregion%processor
        print*, 'fraction = ', sregion%frac

!       If the super region has no children...
        if ( .not. associated(sregion%children) ) then
          print*, 'No children '
          print*

!       Otherwise call recursively for all the children.
        else
          n = size(sregion%children)
          print*, 'Number of children = ', n
          print*
          print*
          do i = 1, n
            call print_tree ( sregion%children(i)%point )
            print*
          end do
        end if
      else
        print*, "Pointer to super region not associated"
        print*
      end if

    end subroutine print_tree


!   Split a super region into an arbitrary number of pieces in direction dir.
!   The number of pieces is size(frac)+1. The integer array frac contains the
!   upper boundary of the n first regions to split and has to be ordered with
!   the largest element smaller than the size of the super region in direction
!   dir. There is currently no check that this is indeed the case.
    subroutine split_sregion ( sregion, dir, frac )
!     The intent has been removed to make it compile with gfortran 4.1.
!      type(superregion2), pointer, intent(in) :: sregion
      type(superregion2), pointer :: sregion
      integer, intent(in) :: dir
      integer, dimension(:), intent(in) :: frac
      type(bbox), dimension(size(frac)+1) :: newboxes
      type(boundary), dimension(size(frac)+1) :: newboundaries
      integer :: n, i, lo, up, str

!     Determine the size of frac.
      n = size(frac)

!     Only do something if n is larger than 0.
      if (n>0) then

!       Allocate storage for the the children pointers.
        allocate(sregion%children(n+1))

!       Determine the low and high boundaries and the stride of the super
!       region.
        lo = sregion%extent%dim(dir)%lower
        up = sregion%extent%dim(dir)%upper
        str = sregion%extent%dim(dir)%stride

!       copy the bbox and boundary information from the super region.
        do i = 1, n+1
          newboxes(i) = sregion%extent
          newboundaries(i) = sregion%outer_boundaries
        end do

!       Update the lower and upper boundaries of the children regions
!       according to the number of points passed in through frac.
!       Also set the outer boundary information to 0.
        do i = 1, n
          newboxes(i)%dim(dir)%upper = lo + ( frac(i) - 1 ) * str
          newboxes(i+1)%dim(dir)%lower = newboxes(i)%dim(dir)%upper + str
          newboundaries(i)%obound(dir,2) = 0
          newboundaries(i+1)%obound(dir,1) = 0
        end do

!       Allocate memory for the children and copy the bbox and outer boundary
!       information in. nullify the pointers of the children to ensure that
!       they are correctly disassociated.
        do i = 1, n+1
          allocate ( sregion%children(i)%point )
          sregion%children(i)%point%extent = newboxes(i)
          sregion%children(i)%point%outer_boundaries = newboundaries(i)
          sregion%children(i)%point%map = sregion%map
          sregion%children(i)%point%processor = sregion%processor
          sregion%children(i)%point%frac = 0.0_wp
          nullify(sregion%children(i)%point%children)
        end do
      end if

    end subroutine split_sregion
          
end module carpet_region
