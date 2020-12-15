      subroutine outcomm_row (buf,iwgid,bwout,mys)
c
c Revised by PK Yeung, 3/24/09, for batches/relay mode
c called by sub. outpen
c
c Clean-up of previous directives:
c --- Implemented JUl1308 and TIMERS_IO
c --- Removed BGW, WRTALL
c Also removed old choice of iwfld (always operate with iwfld=2)
c
c
	use comp
c       use rewrite
	implicit none
	include 'intvars'
c
	complex(b8) :: buf(xisz,ny,zjsz,3+nc)
c
c  routine to write all u-fields on the appropriate logical units
c
 	integer maxsav,lu,i,icall
      parameter (maxsav=100)
      character*6 numer
      character*6 head
      character*16 fn
c
        character*40 longname
        character*4 cpid
        character*2 char_icall
        integer nchar_pid,nchar_icall
c
	character*8 fname
	integer nchar,nf
	integer iwgid,bwout
        integer iys,jys,mys

	real, allocatable :: ureal(:,:),uimag(:,:)
c
c
        integer ii,irow
        integer, allocatable :: xist_row(:),xisz_row(:)
        integer, allocatable :: zjst_row(:),zjsz_row(:)
        complex(b8), allocatable :: buf_row(:,:,:,:)
        complex(b8), allocatable :: buf_row1(:,:)
c
      save icall
      data icall/0/
c
        allocate (xist_row(iproc))
        allocate (xisz_row(iproc))
        allocate (zjst_row(iproc))
        allocate (zjsz_row(iproc))
        allocate (buf_row(xisz,ny/2/mys,iproc,mys))
        if (bwout.eq.1) allocate (buf_row1(xisz,iproc))
c
        call MPI_GATHER (xist,1,MPI_INTEGER,xist_row,1,MPI_INTEGER,0,mpi_comm_row,mpierr)
        call MPI_GATHER (xisz,1,MPI_INTEGER,xisz_row,1,MPI_INTEGER,0,mpi_comm_row,mpierr)
        call MPI_GATHER (zjst,1,MPI_INTEGER,zjst_row,1,MPI_INTEGER,0,mpi_comm_row,mpierr)
        call MPI_GATHER (zjsz,1,MPI_INTEGER,zjsz_row,1,MPI_INTEGER,0,mpi_comm_row,mpierr)

        if (ipid.ne.0) go to 41
c

        write (numer,"(i6)") jpid
        call blanks (numer,nchar)

        write (cpid,"(i3,'/')") iwgid
        call blanks (cpid,nchar_pid)
c
c With fnmu(i)=outpen for i=1 to 3+nc, 
c files to be called outpen.pXXXX if expecting only 1 checkpoint,
c but outpenYY.pXXXX if more than 1 checkpoint
c
      if (isave.eq.0.and.dtsave.eq.0..and.clock_chkpt.eq.0.) then
	fname=fnwu(1)
        else
	icall=icall+1
	head=fnwu(1)
	write (char_icall,"(i2)") icall
	call blanks (char_icall,nchar_icall)
	fname=head//char_icall(1:nchar_icall)
	end if
c
	call blanks (fname,nf)
      longname='outcom/'//cpid(1:nchar_pid)//fname(1:nf)//'.c'//numer(1:nchar)
	call blanks (longname,nchar)
c
      lu=luwu(1)
        write (6,*) 'longname=',longname(1:nchar)
        open (lu,file=longname(1:nchar),form='unformatted')
c
c	allocate (ureal(xisz,ny),uimag(xisz,ny))
c
 41     continue
c
      do 100 i = 1 , 3+nc
c
        if (ipid.eq.0) then
      write( lu ) i , nx , ny , nz
        do irow=1,iproc
      write( lu ) xist_row(irow),xisz,zjst,zjsz,1,ny
        end do
        end if
c
      do 120 zp=1,zjsz
      z=zjst+zp-1
      if( z .eq. nzhp ) go to 120
c
        if (i.eq.2.and.bwout.eq.1) then
        call MPI_GATHER (buf(1,1,zp,i),xisz,mpicomplex,
     1                    buf_row1,xisz,mpicomplex,
     1                    0,mpi_comm_row,mpierr)
        if (ipid.eq.0) then
        write (lu) (buf_row1(:,irow),irow=1,iproc)
        end if

        else
c

        do iys=1,mys
        jys=(iys-1)*ny/2/mys+1
        call MPI_GATHER (buf(1,jys,zp,i),xisz*ny/2/mys,mpicomplex,
     1                    buf_row(1,1,1,iys),xisz*ny/2/mys,mpicomplex,
     1                    0,mpi_comm_row,mpierr)
        end do
        if (ipid.eq.0) then
        write (lu) (((buf_row(:,y,irow,iys),irow=1,iproc),y=1,nyh/mys),iys=1,mys)
        end if
c
        do iys=1,mys
        jys=(iys-1)*ny/2/mys+1+ny/2
        call MPI_GATHER (buf(1,jys,zp,i),xisz*ny/2/mys,mpicomplex,
     1                    buf_row(1,1,1,iys),xisz*ny/2/mys,mpicomplex,
     1                    0,mpi_comm_row,mpierr)
        end do
        if (ipid.eq.0) then
        write (lu) ((buf_row(:,y,irow,1),irow=1,iproc),y=2,nyh/mys),
     1 (((buf_row(:,y,irow,iys),irow=1,iproc),y=1,nyh/mys),iys=2,mys)
        end if
c
c
        end if
c
120   continue
100   continue
c
        deallocate (xist_row,xisz_row,zjst_row,zjsz_row)
        deallocate (buf_row)
        if (bwout.eq.1) deallocate (buf_row1)
c
        if (ipid.ne.0) go to 90
c
c	deallocate (ureal,uimag)
	close (lu)
c
 90     continue
c
      return
      end
