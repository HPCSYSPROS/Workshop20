      subroutine spxyz (f,bsy,bsx,nxdim)
c
#ifdef LAG
c
c 2D code
c
	use mpilag
c

      integer nxdim
      real f(nxdim,zist:zien,yjst:yjen)
      real bsx(nbx,zist:zien,yjst:yjen)
      real, allocatable :: bsz(:,:,:),buf1(:),buf2(:),buf3(:)
      real bsy(bxistart:bxiend,nby,bzjstart:bzjend)

      integer z,y,pos1
c	,yp,yg
      integer x,i,j
c
c form splines in x and z within each slab
c
	write (6,*) 'spxyz, before do 10, taskid=',taskid
      do 10 y=yjst,yjen
c         y=(yg-1)*yjsize+yp
         call splx (f(1,zist,y),nxdim,bsx(1,zist,y),y)
 10   continue
	write (6,*) 'spxyz,  after do 10, taskid=',taskid

c
c pack and send
c

      allocate(buf1(nbx*zisz))
      allocate(bsz(bxistart:bxiend,nbz,yjst:yjen))
	write (6,*) 'spxyz,  A, taskid=',taskid
      
      do y=yjst,yjen
         pos1=1
         do i=0,iproc-1
            do z=zist,zien
               do x=bxist(i),bxien(i)
                  buf1(pos1) = bsx(x,z,y)
                  pos1 = pos1+1
               enddo
            enddo
         enddo

         call mpi_alltoallv(buf1,bcnts_zi,bdisp_zi,mpi_real,
     &        bsz(bxistart,1,y),bcnts_xi,bdisp_xi,mpi_real,
     &        mpi_comm_row,ierr)
	write (6,*) 'spxyz,  A1, taskid,y=',taskid,y

         call splz(bsz(bxistart,1,y))
	write (6,*) 'spxyz,  A2, taskid,y=',taskid,y

      enddo
      deallocate(buf1)

	write (6,*) 'spxyz,  B, taskid=',taskid

      allocate(buf2(bxisize*yjsz*nbz))

      pos1=1
      do j=0,jproc-1
         do z=bzjst(j),bzjen(j)
            do y=yjst,yjen
               do x=bxistart,bxiend
                  buf2(pos1) = bsz(x,z,y)
                  pos1 = pos1+1
               enddo
            enddo
         enddo
      enddo

      deallocate(bsz)
      allocate(buf3(bxisize*nby*bzjsize))

	write (6,*) 'spxyz,  C, taskid=',taskid
      call mpi_alltoallv(buf2,bcnts_yj,bdisp_yj,mpi_real,
     &     buf3,bcnts_zj,bdisp_zj,mpi_real,
     &     mpi_comm_col,ierr)

      deallocate(buf2)

      pos1 = 1
      do j=0,jproc-1
         do z=bzjstart,bzjend
            do y=jjst(j),jjen(j)
               do x=bxistart,bxiend
                  bsy(x,y,z) = buf3(pos1)
                  pos1 = pos1 + 1
               enddo
            enddo
         enddo
      enddo
	write (6,*) 'spxyz,  D, taskid=',taskid


      do z=bzjstart,bzjend
c     z=sz1(zg)+zp-1
         call sply (bsy(bxistart,1,z))
      enddo

      deallocate(buf3)
      
	write (6,*) 'spxyz,  exit, taskid=',taskid


      return
c
 99   call MPI_ABORT (MPI_COMM_WORLD,ierr)
c
#endif
c
      end
