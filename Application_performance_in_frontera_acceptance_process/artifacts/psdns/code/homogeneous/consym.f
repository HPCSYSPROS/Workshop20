c
c 7/13/97: to save some memory, use an external loop
c over the variables
c
c Fix by Matthew Clay:
c 8/5/15: use blocking receive commands to avoid segmentation
c fault associated with non-contiguous buffer arrays.
c
c  routine to enforce conjugate symmetry of initial fields
c  u(kx,ky,kz) contains x-y planes of initial data
c
      subroutine consym (uny,nv)
      use comp
      implicit none
      include 'intvars'
      integer i,j,zdest,zorig,ip,ircont,iscont,ik
      integer, allocatable :: compgrid(:,:)
      integer, allocatable :: sid(:),rid(:)
      integer, allocatable :: mpist(:,:)
      integer sndbuf(5)
      integer nv
      complex(b8) :: uny(ny,zjsz*xisz,nv)
      integer cz,zh,np1l,np1h,np2l,np2h,extra,a
      integer, allocatable :: num_send(:),num_recv(:),stz(:),
     1                        send_req(:,:),recv_req(:)
      integer, allocatable :: mpi_status(:,:,:)
      logical, allocatable :: fsend(:),frecv(:)
      complex(b8), allocatable :: ub(:,:,:), uc(:,:,:),buf(:,:)
      logical flg_first
      integer mpi_status1(MPI_STATUS_SIZE)
      integer :: pp
#ifdef FRESH
      complex(b8), allocatable :: u01(:),u02(:,:)
      if (taskid.eq.0) write (6,*) 'calling consym, nv=',nv
      if (xist.eq.1 .and. zjst.eq.1) then
         do i=1,nv
            uny(1,1,i)=cmplx( real( uny(1,1,i) ) , 0. )
         end do
         do i=1,nv
            do y=nyhp+1,ny
               uny(y,1,i)=conjg(uny(ny+2-y,1,i))
            end do
         end do
      end if
      allocate (compgrid(5,numtasks))
      sndbuf(1)=taskid
      sndbuf(2)=zjst
      sndbuf(3)=zjen
      sndbuf(4)=xist
      sndbuf(5)=xien
      call MPI_ALLGATHER (sndbuf,5,MPI_INTEGER,compgrid,5,
     1                    MPI_INTEGER,MPI_COMM_WORLD,ierr)
c
c Only processes with xist=1 need to worry about conjugate symmetry.
      if (xist.ne.1) go to 89
c
      if(num_al_i(0) .ge. nz) then
         if(taskid .eq. 0) then
c Check if the entire Z dimension is within the first proc.
            do i=1,nv
               do z=nzhp+cut_z(1)+1,nz
                  zorig = nz+2-z
c This translates z into the right index a
c zorig is the same as a_orig since it is before nzhp
                  a = z - 1 - 2*cut_z(1)
c zorig is the same as a_orig since it is before nzhp
c Conjugate/copy the uny array
                  uny(1,a,i) = conjg(uny(1,zorig,i))
                  do y=2,ny
                     uny(y,a,i) = conjg(uny(Ny+2-y,zorig,i))
                  enddo
               enddo
            enddo
         endif
c Need to communicate data
      else
         allocate(fsend(0:jproc-1))
         allocate(frecv(0:jproc-1))
         allocate(num_send(0:jproc-1))
         allocate(num_recv(0:jproc-1))
         allocate(stz(0:jproc-1))
         cz=0
         fsend = .false.
         frecv = .false.
         num_send = 0
         num_recv = 0
c We need lower zs up to this point
         zh = nzhp-cut_z(1) -1
         np1l = 1
         flg_first = .true.
c Determine the number of tasks that span zh
         do i=0,jproc-1
            if(flg_first) then
               cz = cz + num_al_i(i)
               fsend(i) = .true.
               num_send(i) = num_al_i(i)
               if(cz .ge. zh .and. flg_first) then
c Last task in the send list
                  np1h = i
                  extra = cz - zh
                  num_send(i) = num_send(i) - extra
                  flg_first = .false.
               endif
            endif
         enddo
c Determine the number of tasks that span the destination image of zh
         np2l = jproc
         flg_first = .true.
         do i=np1h,jproc-1
            if(flg_first) then
               if(cz .gt. zh) then
                  frecv(i) = .true.
                  if(cz-zh .le. num_al_i(i)) then
c first task in the receive list
                     num_recv(i) = cz-zh
                     stz(i) = zh+num_al_i(i) - cz +1
                     np2l = i
                  else
                     num_recv(i) = num_al_i(i)
                     stz(i) = 1
                  endif
                  if(cz .ge. zh*2-1 .and. flg_first) then
c Last task in the receive list
                     np2h = i
                     num_recv(i) = num_recv(i) - (cz-zh*2) -1
                     flg_first = .false.
                  endif
               endif
               cz = cz + num_al_i(i)
            endif
         enddo
         if(np1l .eq. 0) then
            np1l = 1
         endif
         if(taskid .eq. 0) then
            write(*,*) 'np1l,np1h:', np1l, np1h
            write(*,*) 'np2l,np2h:', np2l, np2h
c Receive data; do not receive from self
 500        format (a,i10)
 100        format (a,i4,a,1p,3i4)
            allocate(ub(ny,nz,nv))
            allocate(uc(ny,nz,nv))

            if(np1h .ge. np1l) then
               allocate(recv_req(np1l:np1h))
               recv_req = 0
            endif

            cz = num_send(0)+1
            do i=np1l,np1h
               write(*,99) 'RECV1: jpid=', jpid, ' from=', i,
     1                     ' cnt=', num_send(i), ' tag=', i
 99            format (a,i4,a,i4,a,i8,a,i4)
               CALL MPI_RECV(ub(1:ny,cz:cz+num_send(i)-1,1:nv),
     1                       num_send(i)*ny*nv,mpicomplex,i,i,
     2                       MPI_COMM_COL,mpi_status1,ierr)
               cz = cz + num_send(i)
            enddo
c Data from self
            do i=1,nv
               do z=2,num_send(0)
                  do y=1,ny
                     ub(y,z,i) = uny(y,z,i)
                  enddo
               enddo
            enddo
c Conjugate the data, prepare send buffer
            do i=1,nv
               do z=2,zh
                  zdest = zh+1-z
                  uc(1,zdest,i) = conjg(ub(1,z,i))
                  do y=2,ny
                     uc(y,zdest,i) = conjg(ub(ny+2-y,z,i))
                  enddo
               enddo
            enddo
c Send out image data; could be also sending to self
            if(np2h .ge. np2l) then
               allocate(send_req(np2l:np2h,nv))
               send_req = 0
            endif
            do j=1,nv
               cz = 1
               do i=np2l,np2h
                  write(*,99) 'SEND2: jpid=', jpid, ' dest=', i,
     1                        ' cnt=', num_recv(i)*ny, ' tag=', j
c                  CALL MPI_SEND(uc(1,cz,j),num_recv(i)*ny,mpicomplex,
c     1                          i,j,MPI_COMM_COL,ierr)
                  CALL MPI_ISEND(uc(1,cz,j),num_recv(i)*ny,mpicomplex,
     1                           i,j,MPI_COMM_COL,send_req(i,j),ierr)
                  cz = cz + num_recv(i)
               enddo
            enddo
c
         else if(fsend(jpid)) then
            write(*,99) 'SEND1: jpid=', jpid, ' dest=', 0,
     1                  ' cnt=', num_send(jpid), ' tag=', jpid
            call MPI_Send(uny(1:ny,1:num_send(jpid),1:nv),
     1                    ny*num_send(jpid)*(nv),mpicomplex,0,jpid,
     2                    MPI_COMM_COL,ierr)
         endif
c All processes to receive data post receives.
         if(frecv(jpid)) then
            allocate(buf(ny,num_recv(jpid)))
            do i=1,nv
               write(*,99) 'RECV2: jpid=', jpid, ' from=', 0,
     1                     ' cnt=', ny*num_recv(jpid), ' tag=', i
               call MPI_recv(buf,ny*num_recv(jpid),mpicomplex,0,i,
     1                       MPI_COMM_COL,mpi_status1,ierr)
               do a=1,num_recv(jpid)
                  do y=1,ny
                     uny(y,a+stz(jpid)-1,i) = buf(y,a)
                  enddo
               enddo
            enddo
            deallocate(buf)
         endif
c Free temporary memory.
         if (taskid.eq.0) then
            if(np2h .ge. np2l) then
               allocate(mpi_status(MPI_STATUS_SIZE,np2l:np2h,nv))
               call MPI_Waitall((np2h-np2l+1)*(nv),send_req(np2l,1),
     1                          mpi_status(1,np2l,1),ierr)
               deallocate(send_req)
               deallocate(mpi_status)
            endif
            deallocate(ub)
            deallocate(uc)
         endif
         deallocate(fsend)
         deallocate(frecv)
         deallocate(num_send)
         deallocate(num_recv)
         deallocate(stz)
      endif
#endif
 89   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
 90   return
      end
