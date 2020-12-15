#include "cctk.h"

module carpet_boxtypes

  use carpet_region
  implicit none

  integer :: ghostsize
  real(wp) :: alpha
  logical :: limit_size
  integer :: procid

  contains

!   The basic calculation of the cost in each of the 3 directions.
    subroutine cost_box ( reg, cost )
      type(bbox), intent(in) :: reg
      integer, dimension(3), intent(out) :: cost
      integer :: i
    
      do i = 1, 3
        cost(i) = (reg%dim(i)%upper -  reg%dim(i)%lower)/reg%dim(i)%stride + 1
      end do
    end subroutine cost_box


!   Find a sorted list of the best options for splitting npoints in 2 chunks
!   among nprocs processors. The fractional remainder is passed in through 
!   fract. The maximum allowed number of options to return is nmax, and the
!   actual number of returned options is nact. The result is returned in
!   np (number of procs for the first chunk) and nc (a measure of the
!   imabalance of the split.
    subroutine choosenprocs ( npoints, nprocs, fract, nmax, nact, np, nc )
      integer, intent(in) :: npoints, nmax
      real(wp), intent(in) :: nprocs, fract
      integer, intent(out) :: nact
      integer, dimension(:), intent(out) :: np
      real(wp), dimension(:), intent(out) :: nc
      integer :: idealcost
      integer :: quot, i, j, round, ind
      real(wp) :: c1, c2, currentimbal
      logical :: exchange

!     The idealcost should be npoints/nprocs but I rescale with nprocs for
!     The purpose of avoiding some problems with floating point comparisons.
      idealcost = npoints

!     Initialize the cost to the maximum possible floating point number.
      nc = minval( empty )

!     initialize the actual returned number of options to 0.
      nact = 0

!     Initialize an integer index to zero. Used to keep track where the current
!     option should be inserted in the result vector.
      ind = 1

!     Calculate the starting number of processors to start the search with
!     (quot) as the largest integer smaller than nprocs.
      quot = floor(nprocs)
!     If the fractional part we have to leave unsplit is larger than 1
!     decrement quot by 1.
      if ( fract > 1.0_wp ) quot = quot-1
!     If nprocs has zero fractional part, we can start at quot/2 due to
!     symmetry.
      if ( .not. ( nprocs>quot ) ) then
        quot = quot / 2
      end if

!     Try all options until reaching a single processor.
      options: do i=quot, 1, -1

!       Calculate the number of points that will yield the smallest imbalance.
        round = nint(real(npoints,wp)/nprocs*i)

!       Calculate the imbalance for the current choice of processors.
!       This has been scaled by nprocs (like idealcost has been above).
        c1 = abs ( (round*nprocs)/i - idealcost )
        c2 = abs ( (npoints-round)/(nprocs-i)*nprocs - idealcost )
        currentimbal = max ( c1, c2 )

!       Initialize the boolean exchange that keeps track of whether the current
!       choice needs to replace an entry in the return vectors.
        exchange = .false.

!       Start from the current number of entries in the return vector and check
!       if the current choice is better than any of them. If so indicate that
!       an exchange has to happen and record the index.
        do j = nact, 1, -1
          if ( currentimbal < nc(j) ) then
            exchange = .true.
            ind = j
          end if
        end do

!       If an exchange has to happen, shift all elements with indices higher
!       than or equal to ind to the right.
        if (exchange) then
          np(ind:nmax) = eoshift(np(ind:nmax),-1)
          nc(ind:nmax) = eoshift(nc(ind:nmax),-1)
        end if

!       If we don't need to do an exchange, and if nact<nmax add the current
!       case after the previous entries.
        if (nact<nmax .and. .not.exchange) then
          nact = nact+1
          np(nact) = i
          nc(nact) = currentimbal
        end if

!       If the return vectors are empty or an exhange has to be done, insert
!       the current choice at the appropriate position (ind) and increment the
!       number of elements in the return vector (capped at nmax)
        if (nact == 0 .or. exchange) then
          nact = min(nact+1,nmax)
          np(ind) = i
          nc(ind) = currentimbal
        end if
      end do options

!     This is commented out for now. Used to return only the elements that are
!     as good as the best case.
!      nbad = nact + 1
!      do i = 2, nact
!        if ( abs(nc(i)-nc(i-1))>1.0d-10 ) nbad = i
!      end do
!      nact = nbad - 1
    end subroutine choosenprocs


!   Calculate a measure of the load imbalance of a super region.
!   The measure is a weighted sum taking into account the imbalance in the
!   number of points (computation) and the imbalance in the number of ghost
!   zones (communication). Ghost zones are added with a multiplication factor
!   of alpha.
    subroutine calc_imbalance ( sregion, imbalance )
!     The intent has been removed to make it compile with gfortran 4.1.
!      type(superregion2), pointer, intent(in) :: sregion
      type(superregion2), pointer :: sregion
      real(wp), intent(out) :: imbalance
      integer, dimension(3) :: cost
      real(wp) :: maxcost, maxgcost, n
      real(wp) :: idealcost, idealgcost

!     Initialize counter and accumulation variables
      n = 0.0_wp
      maxcost = 0.0_wp
      idealcost = 0.0_wp
      maxgcost = 0.0_wp
      idealgcost = 0.0_wp

!     Traverse the tree structure and add values for all leaves.
      call traverse_tree ( sregion )

!     Calculate the ideal cost
      idealcost = idealcost/n
      idealgcost = idealgcost/n

!     Calculate the imbalance
      imbalance = ( maxcost/idealcost - 1.0_wp ) &
                  + alpha * ( maxgcost/idealgcost - 1.0_wp )

      contains

!       Routine to traverse the tree and accumulate the cost of the leaves.
        recursive subroutine traverse_tree ( sreg )
!         The intent has been removed to make it compile with gfortran 4.1.
!          type(superregion2), pointer, intent(in) :: sreg
          type(superregion2), pointer :: sreg
          integer :: ich, nch, tcost, gcost

!         If the region has children recursively traverse the children.
          if ( associated(sreg%children) ) then
            nch = size(sreg%children)
            do ich = 1, nch
              call traverse_tree ( sreg%children(ich)%point )
            end do

!         Otherwise it's a leaf and we accumulate the cost and update the
!         maximum cost.
          else
            call cost_box ( sreg%extent, cost )
            tcost = product ( cost, mask = cost /= 0 )
            gcost = product ( cost+2*ghostsize, mask = cost /= 0 ) - tcost
            tcost = tcost
            gcost = gcost
            maxcost = max ( maxcost, tcost / sreg%frac )
            maxgcost = max ( maxgcost, gcost / sreg%frac )
            idealcost = idealcost + tcost
            idealgcost = idealgcost + gcost
            n = n + sreg%frac
          end if

        end subroutine traverse_tree

    end subroutine calc_imbalance


!   Assign processor numbers to the leaves in a tree starting from start_proc.
    subroutine AssignProcs ( minproc, maxproc, frac1, frac2, sregion )
      integer, intent(in) :: minproc, maxproc
      real(wp ), intent(in) :: frac1, frac2
      logical :: lower_is_outer, upper_is_outer
      integer :: boundarysize, minsize
!     The intent has been removed to make it compile with gfortran 4.1.
!      type(superregion2), pointer, intent(inout) :: sregion
      type(superregion2), pointer :: sregion
      character :: msg*200

!     Initialize.
      if ( frac1 == 0 ) then
        procid = minproc
      else
        procid = minproc + 1
      end if

!     Traverse the tree.
      call traverse_treep ( sregion )
      
      contains

!       Routine to traverse the tree and assign processor ids to the leaves.
        recursive subroutine traverse_treep ( sreg )
!          The intent has been removed to make it compile with gfortran 4.1.
!          type(superregion2), pointer, intent(inout) :: sreg
          type(superregion2), pointer :: sreg
          integer :: ich, nch, maxcost, np
          integer, dimension(3) :: cost
          integer, dimension(1) :: mydim

!         If the region has children recursively traverse the children.
          if ( associated(sreg%children) ) then
            nch = size(sreg%children)
            do ich = 1, nch
              call traverse_treep ( sreg%children(ich)%point )
            end do

!         Otherwise it's a leaf and we assign the processor id and increment it.
          else
            if ( sreg%frac == 1.0_wp ) then
              sreg%processor = procid
              procid = procid + 1
            else
              call cost_box ( sreg%extent, cost )
              maxcost = maxval(cost)
              mydim = maxloc (cost)
              lower_is_outer = sreg%outer_boundaries%obound(mydim(1),1)/=0
              upper_is_outer = sreg%outer_boundaries%obound(mydim(1),2)/=0
              if (frac1 == 0.0_wp) then
                np = 0
              else if (frac2 == 0.0_wp) then
                np = maxcost
              else
                np = nint(frac1/(frac1+frac2)*maxcost)
              end if
              if (limit_size) then
                boundarysize = ghostsize
                minsize = boundarysize+ghostsize
                if (lower_is_outer .and. np>0 .and. np<minsize) then
                  if (np<minsize/2) then
                    np = 0
                  else
                    np = minsize
                  end if
                end if
                if (upper_is_outer .and. &
                    np<maxcost .and. np>maxcost-minsize) then
                  if (np>maxcost-minsize/2) then
                    np = maxcost
                  else
                    np = maxcost-minsize
                  end if
                end if
              end if
              if ((lower_is_outer .and. np>0 .and. np<minsize) .or. &
                  (upper_is_outer .and. &
                   np<maxcost .and. np>maxcost-minsize)) then
                write (msg, '("Not enough grid points: np=",i4," gs=",i4, "maxcost=",i4)') &
                       np, boundarysize, maxcost
                call CCTK_WARN(CCTK_WARN_ABORT, msg)
              end if
              if (np>0 .and. np<maxcost) then
                call split_sregion ( sreg, mydim(1), (/ np /) )
                sreg%processor = -1
                sreg%children(1)%point%processor = minproc
                sreg%children(1)%point%frac = frac1
                sreg%children(2)%point%processor = maxproc
                sreg%children(2)%point%frac = frac2
              else if (np==0) then
                sreg%processor = maxproc
              else if (np==maxcost) then
                sreg%processor = minproc
              end if
              
            end if
          end if

        end subroutine traverse_treep

      end subroutine AssignProcs


!   Routine to sum up the load of an array of super regions. Return the result
!   in the pre-allocated array of length nprocs. The array should be
!   initialized to zero.
    subroutine calc_proc_load ( sregions, proc_load )
      type(ptr), dimension(:), intent(in) :: sregions
      integer, dimension(:), intent(inout) :: proc_load
      integer :: i

!     For each of the super regions in the array.
      do i = 1, size(sregions)

!       Traverse the tree and 
        call traverse_treel ( sregions(i)%point )
      end do

      contains

!       Routine to traverse the tree and add the load of each leaf to the
!       assigned processor.
        recursive subroutine traverse_treel ( sreg )
!         The intent has been removed to make it compile with gfortran 4.1.
!          type(superregion2), pointer, intent(in) :: sreg
          type(superregion2), pointer :: sreg
          integer, dimension(3) :: cost
          integer :: ich, nch

!         If the region has children recursively traverse the children.
          if ( associated(sreg%children) ) then
            nch = size(sreg%children)
            do ich = 1, nch
              call traverse_treel ( sreg%children(ich)%point )
            end do

!         Otherwise it's a leaf and we calculate the cost and accumulate it
!         into the correct spot in the proc_load array.

          else
            call cost_box ( sreg%extent, cost )
            proc_load(sreg%processor+1) = proc_load(sreg%processor+1) + &
                                          product ( cost )
          end if

        end subroutine traverse_treel

    end subroutine calc_proc_load


!   Split a superregion (newreg) recursively into smaller boxes while trying
!   to minimize the load imbalance. Here we allow for fractional processor
!   numbers, nprocs. The fraction that should be left unsplit is given by
!   fract (0.0<=fract<2). The inputs maxchoices and maxfactor are used to
!   determine how many options to examine at a given recursion level. The
!   variable clevel should be initialized to 0 and is used to keep track of
!   the recursion level. The number of options testes is counter and returned
!   in ncount.
    recursive subroutine SplitRegionsRecursive ( nprocs, fract, maxchoices, &
                                                 maxfactor, clevel, &
                                                 newreg, ncount )
      real(wp), intent(in) :: nprocs, fract
      integer, intent(in) :: maxchoices, maxfactor
      integer, intent(inout) :: clevel
!     The intent has been removed to make it compile with gfortran 4.1.
!      type(superregion2), pointer, intent(inout) :: newreg
      type(superregion2), pointer :: newreg
      integer, intent(inout) :: ncount
      type(bbox) :: reg
      integer, dimension(3) :: cost
      integer :: maxcost, ngood, i, np1, lo1, up1, str
      real(wp) :: p1, p2
      integer, dimension(1) :: mydim
      integer, dimension(maxchoices) :: npr
      real(wp), dimension(maxchoices) :: ncr
      type(ptr), dimension(maxchoices) :: newregarr
      real(wp), dimension(maxchoices) :: imbalancearr
      logical :: lower_is_outer, upper_is_outer
      integer :: boundarysize, minsize
      type(bbox) :: reg1, reg2
      integer, dimension(1) :: minimbalanceloc
      integer :: nlocal
      character :: msg*200

!     Extract the bounding box information
      reg = newreg%extent

!     Make sure the array with the test super regions are nullified.
      do i = 1, maxchoices
        nullify ( newregarr(i)%point )
      end do

!     Initialize the imbalance array to max floating point value.
      imbalancearr(:) = minval(empty)

!     Increment the recursion level variable.
      clevel = clevel + 1

!     Calculate the number of choices to examine at the current recursion
!     level based on the value of maxchoices and maxfactor.
      nlocal = max ( maxchoices - ( clevel - 1)  / maxfactor, 1 )

!     If we are called with nprocs==max(1,fract) mark the region to be a leaf
!     in the tree, increment the choice counter and decrement the recursion
!     level counter and return.
      if (nprocs <= max(1.0_wp,fract+0.000001_wp) ) then
        ncount = ncount + 1
        clevel = clevel - 1
        newreg%processor = 1
        newreg%frac = nprocs
        return
      end if

!     Calculate the cost of the current region.
      call cost_box ( reg, cost )

!     Find the dimension with the most points.
      maxcost = maxval(cost)
      mydim = maxloc (cost)

!     Get the list of chunk options to investigate.
      call choosenprocs ( cost(mydim(1)), nprocs, fract, nlocal, ngood, &
                                                                 npr, ncr )

!     Loop over the chunk options to investigate.
      chunks: do i=1, ngood

!       Initialize the super region.
        allocate ( newregarr(i)%point )
        newregarr(i)%point = newreg

!       Find the processor numbers of the two chunks (p1 always an integer).
        p1 = npr(i)
        p2 = nprocs - npr(i)

!       Find the number of points in each of the 2 chunks while making
!       sure the chunks are not too small compared to the ghostsize if
!       limit_size is true.
!       At the outer boundary, take also the boundary size into account,
!       assuming that the boundary size is equal to the ghost size.
!       (Boundary points cannot be split, and near a boundary, the
!        minimum number of interior points is the number of ghost points.)
        np1 = nint((real(maxcost,wp)*p1)/nprocs)
        if (limit_size) then
          lower_is_outer = &
               newregarr(i)%point%outer_boundaries%obound(mydim(1),1)/=0
          upper_is_outer = &
               newregarr(i)%point%outer_boundaries%obound(mydim(1),2)/=0
          boundarysize = ghostsize
          minsize = boundarysize+ghostsize
          if (lower_is_outer .and. np1>0 .and. np1<minsize) then
            np1 = minsize
          end if
          if (upper_is_outer .and. np1<maxcost .and. np1>maxcost-minsize) then
            np1 = maxcost-minsize
          end if
          if ((lower_is_outer .and. np1>0 .and. np1<minsize) .or. &
              (upper_is_outer .and. np1<maxcost .and. np1>maxcost-minsize)) then
            write (msg, '("Not enough grid points: np=",i4," gs=",i4, "maxcost=",i4)') &
                   np1, boundarysize, maxcost
            call CCTK_WARN(CCTK_WARN_ABORT, msg)
          end if
        end if

!       Split the region (reg) into two regions (reg1, reg2) taking the
!       stride into account.
        reg1 = reg
        reg2 = reg
        lo1 = reg%dim(mydim(1))%lower
        up1 = reg%dim(mydim(1))%upper
        str = reg%dim(mydim(1))%stride
        reg1%dim(mydim(1))%upper = lo1 + (np1 - 1)*str
        reg2%dim(mydim(1))%lower = reg1%dim(mydim(1))%upper + str

!       Add the new leaves to the tree.
        call split_sregion ( newregarr(i)%point, mydim(1), (/ np1 /) )

!       Split region reg1 among p1 processors and region reg2 among p2
!       processors adding the new regions to the super region.
        call SplitRegionsRecursive ( p1, fract, maxchoices, maxfactor, &
                                     clevel, &
                                     newregarr(i)%point%children(1)%point, &
                                     ncount )
        call SplitRegionsRecursive ( p2, fract, maxchoices, maxfactor, &
                                     clevel, &
                                     newregarr(i)%point%children(2)%point, &
                                     ncount )

!       Calculate the load imbalance of the resulting super region.
        call calc_imbalance ( newregarr(i)%point, imbalancearr(i) )

!       If the balance is perfect exit the do loop and don't do any further
!       tests.
        if ( imbalancearr(i) == 0.0_wp ) then
          ngood = i
          exit
        end if
      end do chunks

!     Find the location of the minimum load imbalance.
      minimbalanceloc = minloc ( imbalancearr(1:ngood) )      

!     Insert the best choice into the tree structure      
      call point_to_children ( newreg,  newregarr(minimbalanceloc(1))%point )

!     Cut the connection between the best choice and its children.
      call disassociate ( newregarr(minimbalanceloc(1))%point )

!     Deallocate all the tested options.
      do i = 1, ngood
        call destroy_sregion ( newregarr(i)%point )
      end do

!     Decrement the recursion level counter.
      clevel = clevel - 1

    end subroutine SplitRegionsRecursive


!   Routine to return the union of 2 pre-sorted arrays of integers.
!   only n1 elements of list1 and n2 elements of list2 are actually used.
!   The number of elements in the union is n. The union is returned in list.
    subroutine union (n1, list1, n2, list2, n, list )
      integer, intent(in) :: n1, n2
      integer, dimension(:), intent(in) :: list1, list2
      integer, intent(out) :: n
      integer, dimension(:), intent(out) :: list
      integer :: j, k, newl

!     Initialize n and list to 0. The integers j and k are used to march through
!     the 2 input arrays and are initialized to 1 if the lists have some
!     valid elements and 0 otherwise.
      n = 0; list = 0; j = min(1,n1); k = min(1,n2)

!     Initialize the output array to 0
      list = 0

!     If both input arrays are 0, so are the output array and we are done.
      if (n1==0 .and. n2==0) then
        return
      end if

!     Repeat until j>n1 and k>n2.
      loop: do

!       If there are still valid elements in both list1 and list2...
        lists: if ( (j>0 .and. j<=n1) .and. (k>0 .and. k<=n2) ) then

!         The next element to add to list is the minimum of the current values
!         in list1 and list2
          newl = min ( list1(j), list2(k) )

!         If the element is not already in the output list....
          both: if ( n > 0 ) then
            if ( newl > list(n) ) then

!             Add the element
              n = n + 1
              list(n) = newl

!             If the added element comes from list1 increment j
              if ( list1(j) == list(n) ) j = j+1

!             If the added element comes from list2 increment j
              if ( list2(k) == list(n) ) k = k+1

            end if

!         Or if the output list is empty
          else both

!           Add the element
            n = n + 1
            list(n) = newl

!           If the added element comes from list1 increment j
            if ( list1(j) == list(n) ) j = j+1

!           If the added element comes from list2 increment j
            if ( list2(k) == list(n) ) k = k+1

          end if both

!       If there there is only valid elements left in list1
        else if (j>0 .and. j<=n1) then lists

!         Add the current element of list1 to list if it is not already there
!         or list is empty and increment j.
          newl = list1(j)
          list_one: if ( n > 0 ) then
            if ( newl > list(n) ) then
              n = n + 1
              list(n) = newl
              j = j+1
            end if
          else list_one
            n = n + 1
            list(n) = newl
            j = j+1
          end if list_one

!       If there there is only valid elements left in list2
        else if (k>0 .and. k<=n2) then lists

!         Add the current element of list2 to list if it is not already there
!         or list is empty and increment k.
          newl = list2(k)
          list_two: if ( n > 0 ) then
            if ( newl > list(n) ) then
              n = n + 1
              list(n) = newl
              k = k+1
            end if
          else list_two
            n = n + 1
            list(n) = newl
            k = k+1
          end if list_two
        end if lists

!       If both list1 and list2 have been fully traversed then exit.
        if ( j > n1 .and. k > n2 ) exit

      end do loop

    end subroutine union


!   Routine to figure out how many processors to assign to a list of super
!   regions and if necessary split super regions into regions to maintain load
!   balance across all regions. Returns a lits of potentially split super
!   regions and 2 arrays containing the minimum and maximum processor id
!   assigned to all leaf regions in the list of super region tree structures.
    subroutine SplitSuperRegions ( sregion, nprocs )
      type(ptr), dimension(:), intent(inout) :: sregion
      integer, intent(in) :: nprocs
      integer :: i, nregs, totcost, &
                 clevel, ncount, maxcost, np
      logical :: lower_is_outer, upper_is_outer
      integer :: boundarysize, minsize
      integer, dimension(1) :: mydir
      integer, dimension(:), allocatable :: cost
      integer, dimension(3) :: mycost
      integer, dimension(:,:), allocatable :: lcost
      real(wp) :: idealcost, nrprocs
      real(wp), dimension(:), allocatable :: nsplit, accumsplit, frac1, frac2
      real(wp) :: fract
      integer, dimension(:), allocatable :: accumcost, minmap, maxmap, &
                                            splitnum
      character :: msg*200

!     Set nregs to the number of super regions in the input list.
      nregs = size(sregion)

!     Allocate helper arrays. 
      allocate ( cost(nregs), accumcost(nregs+1), minmap(nregs), &
                 maxmap(nregs), nsplit(nregs), accumsplit(nregs+1), &
                 frac1(nregs), frac2(nregs), lcost(nregs,3), splitnum(nregs) )

!     Set cost to the total cost of each super region and accumcost to the
!     accumulated cost.
      accumcost(1) = 0
      do i = 1, nregs
        call cost_box ( sregion(i)%point%extent, lcost(i,:) )
        cost(i) = product ( lcost(i,:) )
        accumcost(i+1) = accumcost(i) + cost(i)
      end do

!     Calculate the total cost (of all super regions) as well as the idealcost.
!     Note that totcost is an integer, while idealcost is real.
      totcost = sum ( cost )
      idealcost = real(totcost,wp) / nprocs

!     Calculate the ideal processor location for the beginning and end of each
!     super region, the fractional number of processors to be assigned to this
!     region as well as the accumulated fractional processor count.
      accumsplit(1) = 0.0_wp
      do i = 1, nregs
        minmap(i) = floor(accumcost(i)/idealcost) 
        maxmap(i) = ceiling(accumcost(i+1)/idealcost) - 1 
        nsplit(i) = cost(i)/idealcost
        accumsplit(i+1) = accumsplit(i)+nsplit(i)
      end do

!     Calculate the size of the fractional regions at the beginning and the
!     end of each super region.
      do i = 1, nregs
        if ( maxmap(i) > minmap(i) ) then
          frac1(i) = ceiling(accumsplit(i))-accumsplit(i)
          frac2(i) = accumsplit(i+1)-floor(accumsplit(i+1))
        else
          frac1(i) = 0.0_wp
          frac2(i) = 0.0_wp
        end if
      end do
 
!     Loop over the super regions.
      do i = 1, nregs
        nrprocs = nsplit(i)
!       If the region begins and ends on the same processor, don't split
!       but just assign the region to that processor.
        if ( maxmap(i) == minmap(i) ) then
          sregion(i)%point%processor = minmap(i)
          sregion(i)%point%frac = nrprocs
!       Otherwise split the region.
        else
!         Initialize counter variables used in the recursive routine.
          clevel = 0
          ncount = 0
!         Calculate the fractional part that has to be left unsplit. This is
!         the sum of the fractional leftovers at either end.
          fract = frac1(i)+frac2(i)
!         If nrprocs > fract (fuzzy logic) then split recursively on nrprocs
!         processors with fract left unsplit (0.0<fract<2.0) and then
!         assign processors to the resulting regions.
          if ( abs ( nrprocs -  fract ) > 0.000001_wp ) then
            call SplitRegionsRecursive ( nrprocs, fract, 4, 3, &
                    clevel, sregion(i)%point, ncount )
            call AssignProcs ( minmap(i), maxmap(i), frac1(i), frac2(i), &
                               sregion(i)%point )
!         Otherwise we just need to split the superregion into 2 regions
!         of size fract1 and fract2.
          else
!           Find the cost of the region.
            call cost_box ( sregion(i)%point%extent, mycost )
!           Find the direction of maximal cost.
            maxcost = maxval(mycost)
            mydir = maxloc (mycost)
!           Find the number of gridpoints in direction mydir of the first
!           split region.
            np = nint(frac1(i)/fract*mycost(mydir(1)))
            if (limit_size) then
              lower_is_outer = &
                   sregion(i)%point%outer_boundaries%obound(mydir(1),1)/=0
              upper_is_outer = &
                   sregion(i)%point%outer_boundaries%obound(mydir(1),2)/=0
              boundarysize = ghostsize
              minsize = boundarysize+ghostsize
              if (lower_is_outer .and. np>0 .and. np<minsize) then
                np = minsize
              end if
              if (upper_is_outer .and. np<maxcost .and. np>maxcost-minsize) then
                np = maxcost-minsize
              end if
              if ((lower_is_outer .and. np>0 .and. np<minsize) .or. &
                  (upper_is_outer .and. &
                   np<maxcost .and. np>maxcost-minsize)) then
                write (msg, '("Not enough grid points: np=",i4," gs=",i4, "maxcost=",i4)') &
                       np, boundarysize, maxcost
                call CCTK_WARN(CCTK_WARN_ABORT, msg)
              end if
            end if
!           Split the region.
            call split_sregion ( sregion(i)%point, mydir(1), (/ np /) )
!           Assign the processors to the resulting regions.
            sregion(i)%point%processor = -1
            sregion(i)%point%children(1)%point%processor = minmap(i)
            sregion(i)%point%children(1)%point%frac = frac1(i)
            sregion(i)%point%children(2)%point%processor = maxmap(i)
            sregion(i)%point%children(2)%point%frac = frac2(i)
          end if
        end if
      end do

!     Deallocate local helper arrays.
      deallocate ( cost, accumcost, minmap, &
                 maxmap, nsplit, accumsplit, &
                 frac1, frac2, lcost, splitnum )
 
    end subroutine SplitSuperRegions

end module carpet_boxtypes
