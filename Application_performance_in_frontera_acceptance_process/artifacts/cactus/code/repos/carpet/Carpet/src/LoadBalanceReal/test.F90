program test

  use carpet_boxtypes
  implicit none

  type(bbox), dimension(100) :: inboxes
  integer, dimension(3,2) :: obound_init =  reshape (source = (/ 1, 1, 1, 1, 1, 1 /), shape = (/ 3, 2 /))
  type(boundary) :: outbound
  integer :: nprocs, i, j, nregs, &
             mode
  type(ptr), dimension(:), allocatable :: sregions
  integer, dimension(:), allocatable :: proc_load
  real(wp) :: totcost, idealcost

  namelist /input/ nregs
  namelist /superregions/ inboxes, nprocs, ghostsize, alpha, limit_size

  open(unit=1,file='boxes.par',status='old',action='read')
  read(1, nml=input)
  close(1)

  open(unit=1,file='boxes.par',status='old',action='read')
  read(1, nml=superregions)
  close(1)

  mode = 2

! In mode 1 the processor decomposition is obtained for all integer values
! from 1 to nprocs and only the processor number and load imbalance is printed.
  if ( mode == 1) then 
    do j = 1, nprocs

      allocate ( sregions(nregs) )
      outbound%obound = obound_init

      do i = 1, nregs
        call create_sregion ( inboxes(i), outbound, i-1, sregions(i)%point )
      end do

      call SplitSuperRegions ( sregions, j )

      allocate ( proc_load(j) )

      proc_load = 0

      call calc_proc_load ( sregions, proc_load )
      totcost = sum ( proc_load )
      idealcost = totcost / j

      do i = 1, nregs
        call destroy_sregion ( sregions(i)%point )
      end do

      print*, j, 100*(maxval(proc_load) - idealcost)/idealcost

      deallocate ( sregions, proc_load )

    end do

! In mode 2 the processor decomposition is done for nprocs only but more datail
! including the full tree information is printed.
  else

    allocate ( sregions(nregs) )
    outbound%obound = obound_init

    do i = 1, nregs
      call create_sregion ( inboxes(i), outbound, i-1, sregions(i)%point )
    end do

    call SplitSuperRegions ( sregions, nprocs )

    do i = 1, nregs
      call print_tree ( sregions(i)%point )
    end do

    allocate ( proc_load(nprocs) )

    proc_load = 0

    call calc_proc_load ( sregions, proc_load )

    print*,'proc_load = ', proc_load 

    totcost = sum ( proc_load )
    print*,'total load = ', totcost

    idealcost = totcost / nprocs
    print*,'ideal load = ', idealcost
    print*,'max load = ', maxval ( proc_load )
    print*,'load_imbalance (%) = ', 100*(maxval(proc_load) - idealcost)/idealcost

    do i = 1, nregs
      call destroy_sregion ( sregions(i)%point )
    end do

    deallocate ( sregions, proc_load )

  end if
  
end program test
