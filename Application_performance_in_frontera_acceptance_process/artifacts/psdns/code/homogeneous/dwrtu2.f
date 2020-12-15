      subroutine dwrtu2(uny)
	use comp
	use timers_io
	implicit none
	include 'intvars'

        complex(b8) :: uny(ny,zjsz*xisz,3+nc)
c
c  routine to write all u-fields on the appropriate logical units
c
 	integer lu,i,numc,icall,numc2,iidir,ii,ndirs,j,k,l
      character*6 numer,numer1,numer2
      character*6 head
      integer relay,pos

	integer,allocatable :: sendbuf(:),recvbuf(:,:)
c
      character*8 fname
#ifdef KRAKEN_JUL11
	character*40 longname
	character*4 cpid
	integer nchar_pid
#endif
c
	integer nn1
c
      save icall
      data icall/0/
	integer ioz,nchar,nchar1
c
#ifdef OUTPUT_PEN
	complex, allocatable :: buf(:,:,:,:)
#endif
c
c if multiple checkpoints are expected the code will
c number the files sequentially
c
	call time_stamp ('enter dwrtu2')
c
#ifndef BGW
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c write computational grid in outdata.grid	
	if (isave.ge.0 .or. nsteps.eq.0) then
	   allocate (sendbuf(5))
	   sendbuf(1)=taskid
	   sendbuf(2)=xist
	   sendbuf(3)=xien
	   sendbuf(4)=zjst
	   sendbuf(5)=zjen
	   if (taskid.eq.0)  allocate (recvbuf(5,numtasks))
	   call MPI_GATHER (sendbuf,5,MPI_INTEGER,recvbuf,5,
     1	   MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	   if (taskid.eq.0) then
	     open (800,file='outdata.grid')
	     write (800,*) numtasks,iproc,jproc
	     do i=1,numtasks
	       write (800,400) recvbuf(:,i)
	     enddo
	     close (800)
	     deallocate(recvbuf)
	   endif
 400    format(5i6)	   
	  deallocate (sendbuf)
	endif
#endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#ifdef BGW
c different directories on gpfs on BGW
	ndirs=256  ! it divides the files among ndirs directories 0/,1/, etc
	ii=mod(taskid,ndirs)
      write (numer2,611) ii
c611  format (i6)
      if (ii.eq.0) then
      numc2=1
      else
      numc2=1+floor(alog10(1.*ii))
      end if
c	print *, 'this proc would write on ', numer2(7-numc2:6)
#endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	
c
	call time_stamp ('dwrtu2: before opening files')
c
	write (numer,"(i6)") taskid
	call blanks (numer,nchar)
c
#ifdef KRAKEN_JUL11
	write (cpid,"(i3,'/')") ipid
	call blanks (cpid,nchar_pid)
      if (isave.eq.0.and.dtsave.eq.0.) then
      longname='outpen/'//cpid(1:nchar_pid)//fnwu(1)//'.p'//numer(1:nchar)
	else
	end if
	open (luwu(1),file=longname,form='unformatted')
	go to 16
#endif
 
	

      if (isave.eq.0.and.dtsave.eq.0.) then
      do 5 i=1,3+nc
      if (i.gt.1.and.luwu(i).eq.luwu(i-1)) go to 5
      call fopen1
#ifdef BGW     
     1 (luwu(i),numer2(7-numc2:6)//'/'//fnwu(i)//'.p'//numer(1:nchar),'unformatted')
#else     
     1 (luwu(i),'outpen/'//fnwu(i)//'.p'//numer(1:nchar),'unformatted')
#endif     
 5    continue
c
      else
c
      icall=icall+1
c
      do 6 i=1,3+nc
      head=fnwu(i)
c
	write (numer1,"(i6)") icall
	call blanks (numer1,nchar1)
	fname=head//numer1(1:nchar1)
	
      if (i.gt.1.and.luwu(i).eq.luwu(i-1)) go to 6
	call blanks (fname,nn1)

      call fopen1
#ifdef BGW     
     1 (luwu(i),numer2(7-numc2:6)//'/'//fname(1:nn1)//'.p'//numer(7-numc:6),'unformatted')
#else     
     1 (luwu(i),'outpen/'//fname(1:nn1)//'.p'//numer(1:nchar),'unformatted')
#endif     
c
 6    continue
c
      end if
 16	call time_stamp ('dwrtu2:  after opening files')
c
	twrite_io = twrite_io - MPI_WTIME()

#ifdef OUTPUT_PEN
      allocate(buf(xisz,ny,zjsz,3+nc))
#endif

c First transpose the data from Y-cylinders to Z-pencils

      do i=1,3+nc
      call kxcomm1_cyl2sq(uny(1,1,i),buf(1,1,1,i),1)
	call time_stamp ('dwrtu2:  after kxcomm1_cyl2sq')

c Transpose from Z-pencils to Y-pencils

      call xkcomm2_pen(buf(1,1,1,i),buf(1,1,1,i),1)      
	call time_stamp ('dwrtu2:  after xkcomm2_pen')
	end do


      pos = 0
      do 100 i = 1 , 3+nc
      lu=luwu(i)
      write( lu ) i , nx , ny , nz
	call time_stamp ('dwrtu2:  after writing i')

#ifdef OUTPUT_PEN

      write( lu ) xist,xisz,zjst,zjsz,1,ny
	call time_stamp ('dwrtu2:  after writing xist')

	if (jpid.eq.0) then
      do 110 zp=1,zjsz
         z=zjst+zp-1
         if( z .eq. nzhp ) go to 110

         write (lu) ((buf(xp,y,zp,i),xp=1,xisz),y=1,nyh)
         write (lu) ((buf(xp,y,zp,i),xp=1,xisz),y=nyhp+1,ny)


 110	continue
	end if
	call time_stamp ('dwrtu2:  after do 110')
#else
c Cylindrical output
        
c Note: different control variables
        write(lu) xist,max_al_x,mystart_x,mystart_z,num_al,ny

        do a=1,num_al
           do y=1,ny
              if(mask(y,a)) then
                 write(lu) uny(y,a,i)
                 pos = pos +1
              endif
           enddo
        enddo

#endif

100   continue

#ifdef OUTPUT_PEN
        deallocate(buf)
#endif


#ifdef OUTPUT_PEN
	ioz=zjsz
	if (zjst.le.nzhp.and.nzhp.le.zjen) ioz=ioz-1
	iwrite_io = iwrite_io+ (ny-1)*ioz*xisz*(3+nc)
#else
        iwrite_io = iwrite_io + pos 
#endif

	twrite_io = twrite_io + MPI_WTIME()
c
c
      do 8 i=1,3+nc
      if (i.gt.1.and.luwu(i).eq.luwu(1)) go to 8
      call fclose1 (luwu(i))
 8    continue
c
	call time_stamp ('dwrtu2:  exit')
c
c
      if (taskid.eq.0) then
      write (6,*) ' report from sub. dwrtu2'
      write (6,600) istep
      end if
 600  format (' velocity and scalar fields saved at step no. ',i7)
c
      return
      end
