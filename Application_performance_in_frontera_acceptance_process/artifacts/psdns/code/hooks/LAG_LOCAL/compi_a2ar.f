
#if defined(CAF) && defined(CF_LAG)

      module mod_compi2_alltoall
      use precision
	use mpilag, only:p8
      implicit none
! Set the stripmine blocksize in bytes; aligned on ialignb byte boundaries
      integer,parameter::sbb=512,ialignb=64
! Set the size of the coarray buffer used during the transfers
      integer,parameter::buffersize=4*1024**2
!      integer,parameter::buffersize=4*1024
      integer, save :: first=1, maxtasks
      real(p8) :: recvbuck(buffersize)[0:*]
!     integer, parameter :: elemsz=p8*2
      integer, parameter :: elemsz=1
      integer, save :: sb=sbb/elemsz, ialign=ialignb/elemsz
      integer, allocatable, save :: myimage(:)
      end module



      subroutine init_compi2_alltoall()
! The first time compi_alltoall[v] is called, we need
! to set up some arrays and a unique random number seed for each image/rank.
! WARNING: PSDNS ALSO SETS ITS OWN RANDOM SEEDS; THERE COULD BE SIDE EFFETS!
      use mod_compi2_alltoall
      use mpicom
      implicit none
      integer seedsize, i
      integer,allocatable::seed(:)
      call random_seed(size=seedsize)    
      allocate(seed(seedsize))
      do i=1,seedsize
        seed(i)=this_image()*12345*i
      enddo
      call random_seed(put=seed)      
        if(taskid.eq.0) write(6,*) 'randomseed compi_lag'
      maxtasks = num_images()/min(iproc,jproc)
      allocate (myimage(0:maxtasks-1))
      myimage=this_image()-1
      first=0
      end subroutine init_compi2_alltoall
  


      subroutine compi2_alltoall(src,des,sendbytes,comm)
      use mod_compi2_alltoall
      implicit none
      
      real(p8), intent(out), dimension(*) :: des
      real(p8), intent(in) , dimension(*) :: src
      integer :: blksize,rri,ri,ti,im,sendcount,sendbytes
      integer :: comm,mycommid,ntasks,ierr
      integer :: i,i1,i2,ib,ie,ibucketbeg,ibucketend,bucketsize,i_co
      integer, allocatable, save :: randlist(:,:),image_map(:,:)
      integer, save :: comms(2)=0, icomms=0
      integer :: icomm, newcomm
      real(p8) x
      
      if(first==1)then
        call init_compi2_alltoall()
      endif

      call MPI_COMM_SIZE (comm,ntasks,ierr)
      call MPI_COMM_RANK (comm,mycommid,ierr)

! On the first call, allocate storage for the lists of ranks in communicators
! and the randomized ists

      if (icomms.eq.0) then
        allocate(image_map(0:maxtasks-1,2))
        allocate(randlist(0:maxtasks-1,2))
      endif

! Determine whether this communicator has been seen before.  
! If not, get and randomize the list of ranks.

      newcomm = 1
! Loop over known communicators
      do i=1,icomms
        if (comms(i).eq.comm) then
! We have seen this communicator before -- nothing to do
          newcomm = 0
          icomm = i
          exit
        endif
      enddo
      if (newcomm.eq.1) then
! New communicator; update bookkeeping
        icomms = icomms + 1
        icomm = icomms
        comms(icomm) = comm

! Generate random node list for icomm

! In general we need to get the global image number for all of 
! the ranks in this communicator.  The only 
! way to do this is via an alltoall
! This adds a lot of overhead for small message alltoalls
! but is easily amoritized for large alltoalls.
! Plus, we do it only on the first call with a given communicator.

        call mpi_alltoall (myimage,1, mpi_integer, 
     &          image_map(0,icomm),1,mpi_integer,comm,ierr)
!       write(6,*)this_image()-1,"image map=",image_map(0:ntasks-1)

        do i=0,ntasks-1
          randlist(i,icomm)=i
        enddo
        do i=0,ntasks-1
          call random_number(x)
          i1=max(0,min(floor(x*ntasks),ntasks-1))
          i2=randlist(i1,icomm)
          randlist(i1,icomm)=randlist(i,icomm)
          randlist(i,icomm)=i2
        enddo
      endif

      sendcount = sendbytes/elemsz

! Adjust the bucketsize for the number of tasks.  
! For relatively small numbers of tasks, the bucketsize 
! will be quite large.  For large numbers of tasks
! the bucketsize will be small.

!RAF ensure buckesize is no less than the stripmine block size in elements
      bucketsize = max((buffersize/(ialign*ntasks))*ialign,sb)

      do ibucketbeg=1,sendcount,bucketsize

! Stripmine loop to communicate only as much data as will
! fit into a bucket for each PE.
        
        ibucketend=min(ibucketbeg+bucketsize-1,sendcount)
        ib=ibucketbeg; ie = ibucketend;

! Copy into the symmetric buffer.  It is this 
! copy that allows us to use PGAS without
! requiring any of the subroutine arguments to be coarrays.

        do ri=0,ntasks-1
          im=mod(ri+mycommid,ntasks)
          recvbuck(im*bucketsize+1:im*bucketsize+(ie-ib)+1)= 
     &         src(im*sendcount+ib:im*sendcount+ie)
        enddo
        
!         sync all

! We need to make sure that all of the copies above are 
! complete before entering the get phase.  We can't use 
! sync all but we can use an mpi_barrier in this hybrid version.

        call mpi_barrier(comm,ierr)
        
! start the global data movement.

        do ib=ibucketbeg,ibucketend,sb

! Stripmine loop to make for very short (sb elemets) messages.  
! On XE and Cascade these short messages flow through the
! network more efficiently, especially when combined with 
! random target ordering.

          ie=min(ib+sb-1,ibucketend)
          do rri=0,ntasks-1

! Randomize the order of the images from which we are 
! pulling data.  This reduces contention significantly.

            ri=randlist(rri,icomm)  
            i_co = image_map(ri,icomm)
!             write(6,*)"image",this_image()-1,"getting from",image_map(ri)
!             !dir$ pgas defer_sync    ! defer_sync seems to be less efficient on XE.
! Get the data from the symmetric buffer directly into the final destination.
            des(ri*sendcount+ib:ri*sendcount+ie) =   
     &              recvbuck(mycommid*bucketsize+ib-ibucketbeg+1:
     &                       mycommid*bucketsize+ie-ibucketbeg+1) [i_co]

          enddo
        enddo
!         sync all 

! In this hybrid version we cannot use sync all, and sync images does not 
! scale well, but we can get the same effect by calling sync memory
! and mpi_barrier on the comm.

        sync memory
        call mpi_barrier(comm,ierr)

      enddo
!       sync all

      end subroutine compi2_alltoall




      subroutine compi2_alltoallv(src,scount,sdispls,maxscount,
     &                           des,rcount,rdispls,comm)
! Note that all the counts are in bytes.
      use mod_compi2_alltoall
      implicit none
      real(p8), intent(out), dimension(*) :: des
      real(p8), intent(in) , dimension(*) :: src
      integer blksize,rri,ri,ti,im,maxrcount
      integer, intent(in) :: scount(0:*),rcount(0:*)
      integer, intent(in) :: sdispls(0:*),rdispls(0:*)
      integer, intent(in) :: maxscount, comm
      integer :: i,i1,i2,ib,ie,ibucketbeg,ibucketend,bucketsize,i_co
      integer :: ntasks,mycommid,ierr
      integer, allocatable, save :: randlist(:,:),image_map(:,:)
      integer, save :: comms(2)=0, icomms=0
      integer :: icomm, newcomm
      real(p8) x

      if(first==1)then
        call init_compi2_alltoall()
      endif


      call MPI_COMM_SIZE (comm,ntasks,ierr)
      call MPI_COMM_RANK (comm,mycommid,ierr)


! On the first call, allocate storage for the lists of ranks in communicators
! and the randomized ists

      if (icomms.eq.0) then
        allocate(image_map(0:maxtasks-1,2))
        allocate(randlist(0:maxtasks-1,2))
      endif

! Determine whether this communicator has been seen before.  
! If not, get and randomize the list of ranks.

      newcomm = 1
! Loop over known communicators
      do i=1,icomms
        if (comms(i).eq.comm) then
! We have seen this communicator before -- nothing to do
          newcomm = 0
          icomm = i
          exit
        endif
      enddo
      if (newcomm.eq.1) then
! New communicator; update bookkeeping
        icomms = icomms + 1
        icomm = icomms
        comms(icomm) = comm

! Generate random node list for icomm

! In general we need to get the global image number for all of 
! the ranks in this communicator.  The only 
! way to do this is via an alltoall
! This adds a lot of overhead for small message alltoalls
! but is easily amoritized for large alltoalls.
! Plus, we do it only on the first call with a given communicator.

        call mpi_alltoall (myimage,1, mpi_integer, 
     &          image_map(0,icomm),1,mpi_integer,comm,ierr)
!       write(6,*)this_image()-1,"image map=",image_map(0:ntasks-1)

        do i=0,ntasks-1
          randlist(i,icomm)=i
        enddo
        do i=0,ntasks-1
          call random_number(x)
          i1=max(0,min(floor(x*ntasks),ntasks-1))
          i2=randlist(i1,icomm)
          randlist(i1,icomm)=randlist(i,icomm)
          randlist(i,icomm)=i2
        enddo
      endif

! Find the largest value of rcount

      maxrcount = rcount(0)
      do i=1,ntasks-1
        if (rcount(i).gt.maxrcount) maxrcount=rcount(i)
      enddo

! Eliminate the need for this AllReduce by passing the max send
! count over all ranks in the communicator from the calling program.
!      lmaxsc = maxval(scount(0:ntasks-1))
!      call mpi_allreduce(lmaxsc,maxscount,1,mpi_integer,mpi_max,comm,ierr)

! Copy as much data at a time that fits into the coarray bucket.
! There is a contribution for each outgoing message.

      bucketsize = max((buffersize/(ialign*ntasks))*ialign,sb)
      do ibucketbeg=1,maxscount/elemsz,bucketsize
        do ri=0,ntasks-1
          im=mod(ri+mycommid,ntasks)
          ibucketend=min(ibucketbeg+bucketsize-1,scount(im)/elemsz)
          ib=ibucketbeg; ie = ibucketend;
!	write(6,*) 'ib,ie,isize,im=',ib,ie,bucketsize,im,this_image()
            recvbuck(im*bucketsize+1:im*bucketsize+(ie-ib)+1)= 
     &          src(sdispls(im)/elemsz+ib:sdispls(im)/elemsz+ie)
        enddo

! Before going on, make sure all copies have completed.
!         sync all
        call mpi_barrier(comm,ierr)

! Now copy a coarray bucketfull of data from each image to the
! local receive buffer.  The ordering of the ranks from which we
! get the data is randomized to reduce contention significantly.

        do ib=ibucketbeg,min(ibucketbeg+bucketsize-1,maxrcount/elemsz),sb
          do rri=0,ntasks-1
            ri=randlist(rri,icomm)
            i_co = image_map(ri,icomm)
            ie=min(ib+sb-1,rcount(ri)/elemsz)
            des(rdispls(ri)/elemsz+ib:rdispls(ri)/elemsz+ie) =   
     &                recvbuck(mycommid*bucketsize+ib-ibucketbeg+1:
     &                         mycommid*bucketsize+ie-ibucketbeg+1) [i_co]
          enddo
        enddo

! Before going on to the next bucket-full, make sure all copies have 
! been completed.

!         sync all 
        sync memory
        call mpi_barrier(comm,ierr)
      enddo
!       sync all
      return
      end subroutine compi2_alltoallv

#endif
