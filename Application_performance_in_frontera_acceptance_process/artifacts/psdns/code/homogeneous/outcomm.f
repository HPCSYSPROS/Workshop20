      subroutine outcomm (uny)
c
c Revised by PK Yeung, 3/24/09, for batches/relay mode
c
c Clean-up of previous directives:
c --- Implemented JUl1308 and TIMERS_IO
c --- Removed BGW, WRTALL
c Also removed old choice of iwfld (always operate with iwfld=2)
c
c
      USE IO
	use comp
c       use rewrite
	implicit none
	include 'intvars'
c
	complex(b8) :: uny(ny,xisz*zjsz,3+nc)
c
c  routine to write all u-fields on the appropriate logical units
c
 	integer maxsav,lu,i,ndirs
      parameter (maxsav=100)
      character*6 numer,numer2
      character*6 head
      character*16 fn
c
        character*11 form,unform
      data form,unform/'formatted','unformatted'/
        character*20 caux
c
	integer,allocatable :: sendbuf(:),recvbuf(:,:)
	integer ip2
c
	character*8 fname
c
	character*25 string
c
        integer iwrite,mwrite,nwrite,relay,jwrite,jsource
        integer mpistatus (MPI_STATUS_SIZE)

#ifdef OUTPUT_PEN
        complex, allocatable :: buf(:,:,:,:),buf1(:,:,:)
#endif
	logical bwexs
	integer idummy,bwout,lubw
c
        integer mys
        logical flag_mys

c
	if (taskid.eq.0) write (6,*) 'enter outpen_comm'

	bwout=0
	if (taskid.eq.0) then
	inquire(file='bwio',exist=bwexs)
	if (bwexs) then
	lubw=819
	open(unit=lubw,file='bwio',action='read')
	read (lubw,'(1x)')
	read (lubw,*) idummy,bwout
	close(lubw)
	endif
	endif
	call MPI_BCAST (bwout,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
c
        fnwu(1)(4:6)='row'
	if (bwout.eq.1) fnwu(1)(4:6)='rwb'
        if (taskid.le.1) write (6,*) 'outcomm: bwout,fnwu(1)=',bwout,fnwu(1)

#ifdef OUTPUT_PEN
      allocate(buf1(xisz,nz,yjsz))
      allocate(buf(xisz,ny,zjsz,3+nc))

c First transpose the data from Y-cylinders to Z-pencils

      do i=1,3+nc
c
      call kxcomm1_cyl2sq(uny(1,1,i),buf1,1)

c Transpose from Z-pencils to Y-pencils

      call xkcomm2_pen(buf1,buf(1,1,1,i),1)
c
        end do
c
	deallocate (buf1)
#endif
c
c write computational grid in outdata.grid	
	if (isave.ge.0 .or. nsteps.eq.0) then
	   allocate (sendbuf(5))
	   sendbuf(1)=taskid
	   sendbuf(2)=xist
	   sendbuf(3)=xien
	   sendbuf(4)=zjst
	   sendbuf(5)=zjen
	     allocate (recvbuf(5,numtasks))
	   call MPI_GATHER (sendbuf,5,MPI_INTEGER,recvbuf,5,
     1	   MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	   if (taskid.eq.0) then
        caux='outdata.grid'
                call fopen1 (800,caux,form)
c	     open (800,file='outdata.grid')
	     write (800,*) numtasks,iproc,jproc
	     do i=1,numtasks
	       write (800,400) recvbuf(:,i)
	     enddo
	     close (800)
	   endif
 400    format(5i6)	   
	  deallocate (sendbuf)
	     deallocate(recvbuf)
	endif
c
        if (taskid.eq.0) then
        mys=1
        inquire (file='mys',exist=flag_mys)
        if (flag_mys) then
        open (300,file='mys')
        read (300,*) mys
        close (300)
        end if
        write (6,*) 'outcomm: value of mys=',mys
        end if
        call MPI_BCAST (mys,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)

        call MPI_BARRIER (MPI_COMM_WORLD,mpierr)
c relay scheme
c
      if (jproc.le.1024) mwrite=jproc
      if (mod(jproc,1024).eq.0) mwrite=1024
      nwrite=jproc/mwrite
      iwrite=jpid/mwrite
      ndirs=2**(ip2(jproc)/2)
      jwrite=mod(jpid,ndirs)
      if (taskid.eq.0) write (6,*) 'outcomm: jproc,ndirs,mwrite=',
     1   jproc,ndirs,mwrite
c
      if (iwrite.eq.0) then
         if (hdf5_output_flag) then
            CALL OUTROW_HDF5(buf)
         else
            call outcomm_row (buf,jwrite,bwout,mys)
         endif
         relay=taskid
      else
         jsource=taskid-mwrite*iproc
         call MPI_RECV (relay,1,MPI_INTEGER,jsource,jsource,
     1                  MPI_COMM_WORLD,mpistatus,mpierr)
         if (hdf5_output_flag) then
            CALL OUTROW_HDF5(buf)
         else
            call outcomm_row (buf,jwrite,bwout,mys)
         endif
      end if
c
      if (iwrite.lt.nwrite-1) then
         call MPI_SSEND (relay,1,MPI_INTEGER,taskid+mwrite*iproc,taskid,
     1                   MPI_COMM_WORLD,mpierr)
      end if
c
      call MPI_BARRIER (MPI_COMM_WORLD,mpierr)
c
#ifdef OUTPUT_PEN
      deallocate (buf)
#endif
c
      if (taskid.eq.0) then
         write (6,*) ' report from sub. outpen'
         write (6,600) istep
      end if
 600  format (' velocity and scalar fields saved at step no. ',i7)
c
      return
      end
#ifdef NOT_NEEDED
      function ip2(n)
c
c  function to determine if n is a positive power of two.
c  if n is an integer power of two, ip2 is that integer.
c  if not, ip2 is zero
c
      integer ip2,n
c
      ip2 = 0
      if( n .le. 0 ) return
      l = n
100   ip2 = ip2 + 1
      modl = mod(l,2)
      if( modl .eq. 0 ) then
         l = l/2
         go to 100
      else if ( modl .eq. 1 ) then
c corrected by p.k. 2/90
         ip2 = ip2 -1
         return
      else
         ip2 = 0
         return
      endif
c
      end
#endif
