c
c 7/13/97: to save some memory, use an external loop
c over the variables
c
c  routine to enforce conjugate symmetry of initial fields
c  u(kx,ky,kz) contains x-y planes of initial data
c
      subroutine enfsym (uny)
      use comp
      implicit none
      include 'intvars'

      integer i,j,zdest,zorig,ip,ircont,iscont,ik
      integer, allocatable :: compgrid(:,:)
      integer, allocatable :: sid(:),rid(:)
      integer, allocatable :: mpist(:,:)
      integer sndbuf(5)
      complex(b8) :: uny(ny,zjsz*xisz,3+nc)
      integer cz,zh,np1l,np1h,np2l,np2h,extra,a
      integer, allocatable :: num_send(:),num_recv(:),stz(:),send_req(:,:),recv_req(:)
      integer, allocatable :: mpi_status(:,:,:)
      logical, allocatable :: fsend(:),frecv(:)
      complex(b8), allocatable :: ub(:,:,:), uc(:,:,:),buf(:,:)
      logical flg_first
      integer mpi_status1(MPI_STATUS_SIZE)

#ifdef FRESH
      complex(b8), allocatable :: u01(:),u02(:,:)
c     common/sym/u01(ny,mz),u02(ny,nz)
c     complex u01,u02
c     
c     
c     kz=0
c     
        write (6,*) 'enter enfsym, taskid=',taskid
c     if (taskid.eq.0) then
      if (xist.eq.1 .and. zjst.eq.1) then
         do i=1,3+nc
            uny(1,1,i)=cmplx( real( uny(1,1,i) ) , 0. )
         end do
         do i=1,3+nc
            do y=nyhp+1,ny
               uny(y,1,i)=conjg(uny(ny+2-y,1,i))
            end do
         end do
      end if
c     
!      allocate (sndbuf(5))
      allocate (compgrid(5,numtasks))
      sndbuf(1)=taskid
      sndbuf(2)=zjst
      sndbuf(3)=zjen
      sndbuf(4)=xist
      sndbuf(5)=xien
      call MPI_ALLGATHER (sndbuf,5,MPI_INTEGER,compgrid,5,
     1	   MPI_INTEGER,MPI_COMM_WORLD,ierr)
c

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (xist.ne.1) go to 89
c

      if(num_al_i(0) .ge. nz) then
         if(taskid .eq. 0) then
c Check if the entire Z dimension is within the first proc.
c            print *,'Enfsym: no communication is required'
            do i=1,3+nc
               do z=nzhp+cut_z(1)+1,nz
                  zorig = nz+2-z
c     This translates z into the right index a
c     zorig is the same as a_orig since it is before nzhp
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
      else
c Need to communicate data
c
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
c We need lower z's up to this point
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
c     Last task in the receive list
                     np2h = i
                     num_recv(i) = num_recv(i) - (cz-zh*2) -1
                     flg_first = .false.
                  endif
               endif
               cz = cz + num_al_i(i)
            endif
         enddo

c         print *,'Enfsym: tasks ',np1l,':',np1h,' sending, tasks ',np2l,':',np2h,' receiving'
c         print *,'zh = ',zh
c         print *,taskid,': num_al = ',num_al
c         if(taskid .eq. 0) then
c            print *,'num_recv: ',num_recv
c            print *,'num_send: ',num_send
c         endif

         if(np1l .eq. 0) then
            np1l = 1
         endif

         if(taskid .eq. 0) then
c Receive data; don't receive from self
            allocate(ub(ny,nz,3+nc))
            allocate(uc(ny,nz,3+nc))

            if(np1h .ge. np1l) then
               allocate(recv_req(np1l:np1h))
               recv_req = 0
            endif

            cz = num_send(0)+1
            do i=np1l,np1h
               call MPI_irecv(ub(1:ny,cz:cz+num_send(i)-1,1:3+nc),
     &     num_send(i)*ny*(3+nc),mpicomplex,i,i,MPI_COMM_COL,recv_req(i),ierr)
               cz = cz + num_send(i)
            enddo
c Data from self
            do i=1,3+nc
               do z=2,num_send(0)
                  do y=1,ny
                     ub(y,z,i) = uny(y,z,i)
                  enddo
               enddo
            enddo

            if(np1h .ge. np1l) then
c               allocate(mpi_status(MPI_STATUS_SIZE,np1l:np1h,1))
               do i=np1l,np1h
                  call MPI_Wait(recv_req(i),mpi_status1,ierr)
               enddo
               deallocate(recv_req)
c               deallocate(mpi_status)
            endif

c Conjugate the data, prepare send buffer
            do i=1,3+nc
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
               allocate(send_req(np2l:np2h,3+nc))
               send_req = 0
            endif
            do j=1,3+nc
               cz = 1
               do i=np2l,np2h
                  call MPI_Isend(uc(1,cz,j),
     &    num_recv(i)*ny,mpicomplex,i,j,MPI_COMM_COL,send_req(i,j),ierr)
                  cz = cz + num_recv(i)
               enddo
            enddo

         else if(fsend(jpid)) then
            call MPI_Send(uny(1:ny,1:num_send(jpid),1:3+nc),
     &        ny*num_send(jpid)*(3+nc),mpicomplex,0,jpid, MPI_COMM_COL,ierr)
         endif
c
c
         if(frecv(jpid)) then
            allocate(buf(ny,num_recv(jpid)))
            do i=1,3+nc
               call MPI_recv(buf,ny*num_recv(jpid),mpicomplex,0,i, 
     1           MPI_COMM_COL,mpi_status1,ierr)
               do a=1,num_recv(jpid)
                  do y=1,ny
                     uny(y,a+stz(jpid)-1,i) = buf(y,a)
                  enddo
               enddo
            enddo
            deallocate(buf)
         endif
         if(taskid .eq. 0) then
            if(np2h .ge. np2l) then
               allocate(mpi_status(MPI_STATUS_SIZE,np2l:np2h,3+nc))
               call MPI_Waitall((np2h-np2l+1)*(3+nc),send_req(np2l,1),
     *           mpi_status(1,np2l,1),ierr)
               deallocate(send_req)
               deallocate(mpi_status)
            endif
            deallocate(ub)
            deallocate(uc)
         endif

c
         deallocate(fsend)
         deallocate(frecv)
         deallocate(num_send)
         deallocate(num_recv)
         deallocate(stz)

      endif

c     check conjugate symmerty	

c      print *,'Enfsym: checking conjugate symmetry'
c      xp = mystart_x
c      z=mystart_z
c      do a=1,num_al
c         if(xp .eq. 1) then
c            print *,'x,z=',xp,z
c            do y=1,ny
c               if(abs(uny(y,a,3)) .gt. 0.00001) then
c                  print *,y,uny(y,a,3)
c!                  ik = 1.5 + sqrt(kx(xp)**2+ky(y)**2+kz(z)**2)
c!                  if(ik .lt. 5) then
c!                     print *,'Warning in enfsym: ',ik,y,z,xp,uny(y,a,3)
c!                  endif
c               endif
c            enddo
c         endif
c         call next_xz(xp,z)
c      enddo

c     do i=1,3
c     write (i*2000+taskid,*) xist,xien,zjst,zjen
c     do y=1,ny
c     do zp=1,zjsz
c     z=zp+zjst-1
c     do xp=1,xisz
c     x=xp+xist-1
c     if (x.eq.1.)
c     1 write (i*2000+taskid,*) x,y,z,kx(x),ky(y),kz(z),uny(xp,zp,y,i) c     enddo
c     enddo
c     enddo
c     enddo
c     close(i*2000+taskid)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#endif
 89   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
 90   return
      end
