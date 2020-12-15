      subroutine incomm_row  (uny,iproc0,jproc0,ixfst,ixlst,jzfst,jzlst,ndirs,ind_prec,bwin,relay,itask, mys)
      use comp
      use timers_io
      implicit none
      include 'intvars'
c
c Changes (re: "mys") by PKY, 12/9/2015
      
c It doesn't really matter how we declare uny here; 
c Upon exit, uny will contain data in cylindrical form, that is 
c (ny,zjsz*xisz,3+nc)

      complex(b8) :: uny(ny,zjsz,xisz,3+nc)
        real(4), allocatable :: unpen_s_r(:,:)
        real(4), allocatable :: unpen_s_i(:,:)
        real(8), allocatable :: unpen_d_r(:,:)
        real(8), allocatable :: unpen_d_i(:,:)
c
      complex(b8), allocatable :: unpen(:,:),buf(:,:,:),utemp(:,:)

      integer nchar, nchar1
      integer lu,ii,jj,kk,iproc0,jproc0,ip
      integer ixfst(irz),ixlst(irz),jzfst(irz),jzlst(irz)
      integer i,iin,jjn,nxin,nyin,nzin,a,sty,eny,szy
      integer  xistin,xiszin,zjstin,zjszin,iaux
      integer jxp1,jxp2,jz
      integer dmy(3),hms(3),date_time(8)
c
        integer irow,bwin
        integer,allocatable :: xistin_row(:)
        integer,allocatable :: xiszin_row(:)
        complex(b8), allocatable :: buf_row1(:,:)
        complex(b8), allocatable :: buf_row(:,:,:,:)


      logical casex1,casex2,casez1,casez2          
      logical :: opn         
      character*20 string
      character*80 :: nam
      character*6 numer
      character*100 fn
      character*11 unform

      data unform/'unformatted'/

        integer ndirs,jwrite,nchar_pid,luout,ndirs0,relay,ip2,itask
        character*4 cpid
        real*8 rtime1,rtime2
        real cpuread
        character*30 fnout

        integer iys,jys,mys

        
	integer ind_prec

        luout=6
        if (itask.eq.taskid.and.relay.ne.numtasks) then
        if (relay.gt.0) luout=60
        else
        if (itask.ge.iproc) luout=60
        end if
c
        rtime1=MPI_WTIME()
c
        if (relay.gt.0) then
        ndirs0=2**(ip2(numtasks)/2)
        if (taskid.eq.1) write (6,*) 'ndirs0=',ndirs0
        jwrite=mod(taskid,ndirs0)
        write (cpid,"(i3,'/')") jwrite
        call blanks (cpid,nchar_pid)

        write (numer,"(i6)") taskid
        call blanks (numer,nchar1)
            fnout='readtimes'//'/'//cpid(1:nchar_pid)
            call blanks (fnout,nchar)
            fnout=fnout(1:nchar)//'time.'//numer(1:nchar1)
            call blanks (fnout,nchar)
        open (luout,file=fnout(1:nchar))
        end if


      call date_and_time (values=date_time)
      dmy(1:3)=date_time(3:1:-1)
      hms(1:3)=date_time(5:7)
      
	if (taskid.le.1) then
      string='enter incomm_row'
      write (luout,610) string,taskid,dmy(2),dmy(1),dmy(3),
     1            hms,ndirs,irz
 610  format (a20,2x,i5, ' date & time is  ',i2,'/',i2,
     1     '/',i4,2x,i2,':',i2,':',i2,i3,i5)
	end if
c


        allocate (xistin_row(iproc))
        allocate (xiszin_row(iproc))

      do 100 kk=1,irz
         
         casex1=.false.
         casex2=.false.
         casez1=.false.
         casez2=.false.
         
         ip=kk-1

         if (xien.ge.ixfst(kk).and.xien.le.ixlst(kk)) casex1=.true.
         if (ixlst(kk).ge.xist.and.ixlst(kk).le.xien) casex2=.true.
         if (zjen.ge.jzfst(kk).and.zjen.le.jzlst(kk)) casez1=.true.
         if (jzlst(kk).ge.zjst.and.jzlst(kk).le.zjen) casez2=.true.

         if ((casex1.or.casex2).and.(casez1.or.casez2)) then
            write (luout,601) ip, taskid,xist,xien,zjst,zjen,
     1           ixfst(kk),ixlst(kk),jzfst(kk),jzlst(kk),
     1           casex1,casex2,casez1,casez2
 601        format ('incomm_row: ip', i5,1x,9i4,4L3)

c
	if (ipid.eq.0) then
	write (numer,"(i5)") ip/iproc0
	call blanks (numer,nchar1)
 	write (luout,*) 'incomm_row: ndirs=',taskid,ndirs
	        jwrite=mod(ip/iproc0,ndirs)
        write (cpid,"(i3,'/')") jwrite
            call blanks (indir_fn,nchar)
        call blanks (cpid,nchar_pid)
            fn=indir_fn(1:nchar)//'/'//cpid(1:nchar_pid)
            call blanks (fn,nchar)
            fn=fn(1:nchar)//fninit(1)//'.c'//numer(1:nchar1)
            call blanks (fn,nchar)
 	write (luout,*) 'incomm_row: fn=',taskid,fn(1:nchar)
c
            lu=luinit(1)
c            call fopen1 (lu,fn(1:nchar),unform)
            open (lu,file=fn(1:nchar),form='unformatted',action='read')
	end if
c
            do 210 i=1,3+nc                   
               if (kinit(i).gt.0) then                      
c
                if (ipid.eq.0) then
c
                  iin=kinit(i)                      
                  read( lu , err=10 , end=20 ) jjn , nxin , nyin , nzin                      
                do ii=1,iproc0
                     read( lu , err=10 , end=20 )
     1               xistin_row(ii), xiszin_row(ii) ,zjstin,zjszin,iaux,nyin
                end do
c
                end if

                call MPI_SCATTER (xistin_row,1,MPI_INTEGER,
     1                            xistin,1,MPI_INTEGER,
     1                    0,mpi_comm_row,mpierr)
                call MPI_SCATTER (xiszin_row,1,MPI_INTEGER,
     1                            xiszin,1,MPI_INTEGER,
     1                    0,mpi_comm_row,mpierr)
c
        allocate (buf_row(xiszin,ny/2/mys,iproc0,mys))
        if (i.eq.2.and.bwin.eq.1) allocate (buf_row1(xiszin,iproc0))
        allocate (unpen(xiszin,ny))
c
c
                  do 220 jz=jzfst(kk),jzlst(kk)

                     if( jz .eq. nzhp ) go to 220

        if (i.eq.2.and.bwin.eq.1) then
c
        if (ipid.eq.0) then
        read (lu,err=30,end=40) (buf_row1(:,irow),irow=1,iproc0)
        end if
        call MPI_SCATTER (buf_row1,xiszin,mpicomplex,
     1                   unpen(1,1),xiszin,mpicomplex,
     1                    0,mpi_comm_row,mpierr)
c
        else
c
c 
c--------------------------
c 12/9/2015: replacing lines below by a new version which
c can use smaller messages in MPI_SCATTER.
c (The parameter nys can be 1, 2, 4 etc towards smaller messages)
c (Note in original version buf_row was a 3D array, with last
c dimension equal to ny)
c 
#ifdef TEST
        if (ipid.eq.0) then
        read (lu,err=30,end=40) ((buf_row(:,y,irow),irow=1,iproc0),y=1,nyin/2)
        read (lu,err=30,end=40) ((buf_row(:,y,irow),irow=1,iproc0),y=nyhp+1,ny)
        end if
        call MPI_SCATTER (buf_row,xiszin*ny,mpicomplex,
     1                   unpen,xiszin*ny,mpicomplex,
     1                    0,mpi_comm_row,mpierr)
#endif
c--------------------------
c 
c Fourier modes with 1<=y<=ny/2
        if (ipid.eq.0) then
        read (lu,err=30,end=40) (((buf_row(:,y,irow,iys),irow=1,iproc0),y=1,nyin/2/mys),iys=1,mys)
        end if
        do iys=1,mys
        jys=(iys-1)*ny/2/mys+1
        call MPI_SCATTER (buf_row(1,1,1,iys),xiszin*ny/2/mys,mpicomplex,
     1                   unpen(1,jys),xiszin*ny/2/mys,mpicomplex,
     1                    0,mpi_comm_row,mpierr)
        end do

c Fourier modes with y >= nyhp
        if (ipid.eq.0) then
        buf_row(:,1,:,1)=0.
        read (lu,err=30,end=40)
     1 ((buf_row(:,y,irow,1),irow=1,iproc0),y=2,nyin/2/mys),
     1 (((buf_row(:,y,irow,iys),irow=1,iproc0),y=1,nyin/2/mys),iys=2,mys)
        end if
        do iys=1,mys
        jys=(iys-1)*ny/2/mys+1+ny/2
        call MPI_SCATTER (buf_row(1,1,1,iys),xiszin*ny/2/mys,mpicomplex,
     1                   unpen(1,jys),xiszin*ny/2/mys,mpicomplex,
     1                    0,mpi_comm_row,mpierr)
        end do
c
c-----------------------------------------------

        end if
c

                     if (jz.lt.zjst.or.jz.gt.zjen) go to 220

                     zp=jz-zjst+1

                     if (iproc.ge.iproc0) then

                        jxp1=xist-ixfst(kk)+1
                        jxp2=jxp1+xisz-1
!     
                        do y=1,ny/2
                           uny(y,zp,1:xisz,i)=unpen(jxp1:jxp2,y)
                        end do
                        do y=ny/2+2,ny
                           uny(y,zp,1:xisz,i)=unpen(jxp1:jxp2,y)
                        end do

                     else

                        jxp1=ixfst(kk)-xist+1
                        jxp2=jxp1+xiszin-1

                        do y=1,ny/2
                           uny(y,zp,jxp1:jxp2,i)=unpen(1:xiszin,y)
                        end do
                        do y=ny/2+2,ny
                           uny(y,zp,jxp1:jxp2,i)=unpen(1:xiszin,y)
                        end do

                     end if
                     
 220              continue

                  deallocate (unpen,buf_row)
        if (i.eq.2.and.bwin.eq.1) deallocate (buf_row1)


               end if

 210        continue
c


 
         if (ipid.eq.0) close (lu)

c
         end if

 100  continue


        if (bwin.eq.1) then
        do x=1,xisz
        xp=x+xist-1
        do z=1,zjsz
        zp = z+zjst-1
        do y=2,ny/2
        uny(y,z,x,2) = -(kx(xp)*uny(y,z,x,1)*b11(2)+kz(zp)*uny(y,z,x,3)*b33(2))/ky(y)/b22(2)
        enddo
        do y=nyhp+1,ny
        uny(y,z,x,2) = -(kx(xp)*uny(y,z,x,1)*b11(2)+kz(zp)*uny(y,z,x,3)*b33(2))/ky(y)/b22(2)
        enddo
        enddo
        enddo
        end if
c
        rtime2=MPI_WTIME()
        cpuread=rtime2-rtime1
      call date_and_time (values=date_time)
      dmy(1:3)=date_time(3:1:-1)
      hms(1:3)=date_time(5:7)
         string=' exit incomm_row'
         write (luout,620)
     1  string(1:20),taskid,dmy(2),dmy(1),dmy(3),hms,cpuread,ipid,jpid
        if (cpuread.gt.100.) then
         write (6,620)
     1  string(1:20),taskid,dmy(2),dmy(1),dmy(3),hms,cpuread,ipid,jpid
        end if
 620  format (a20,2x,i6, ' date & time is  ',i2,'/',i2,
     1     '/',i4,2x,i2,':',i2,':',i2,f6.1,' secs',i3,i5)

c     endif

      inquire (60, opened=opn, name=nam)
        if (opn) close (60)


 90   return

 10   write(6,*)' error reading field header: jjn=',jjn
      write(6,*)' stopped in incomm_row, taskid= ',taskid,fn(1:nchar)
      inquire (lu, opened=opn, name=nam)
      write(6,*)' reading lu,filename=',lu,nam
      go to 99

 20   write(6,*)' end of file hit reading field header: jjn=',jjn
      write(6,*)' stopped in incomm_row, taskid= ',taskid,fn(1:nchar)
      inquire (lu, opened=opn, name=nam)
      write(6,*)' reading lu,filename=',lu,nam
      go to 99

 30   write(6,*)' error reading field:i,iin',i,iin,taskid
      write(6,*)' stopped in incomm_row, taskid= ',taskid,fn(1:nchar)
      inquire (lu, opened=opn, name=nam)
      write(6,*)' reading lu,filename=',lu,nam
      go to 99

 40   write(6,*)' end of file hit reading field:i,iin',i,iin,jz
      write(6,*)' stopped in incomm_row, taskid= ',taskid,fn(1:nchar)
      inquire (lu, opened=opn, name=nam)
      write(6,*)' reading lu,filename=',lu,nam
      go to 99

 99   call abrt('from incomm_row')

      return
      end subroutine incomm_row
