

!     NPM converting to stride-1 3/9/08
      subroutine inpen_read_bw  (uny,slabpenc,iproc0,ixfst,ixlst,jzfst,jzlst,ndirs,ind_prec,relay,xi1,zj1,itask)
      use comp
      use timers_io
      implicit none
      include 'intvars'
      
c It doesnt really matter how we declare uny here; 
c Upon exit, uny will contain data in cylindrical form, that is 
c (ny,zjsz*xisz,3+nc)


      complex(b8) :: uny(ny,zjsz,xisz,3+nc)
        real(4), allocatable :: unpen_s_r(:,:)
        real(4), allocatable :: unpen_s_i(:,:)
        real(8), allocatable :: unpen_d_r(:,:)
        real(8), allocatable :: unpen_d_i(:,:)



      complex(b8), allocatable :: unpen(:,:),buf(:,:,:),utemp(:,:)

      integer nchar, nchar1, nchar0
      integer lu,ii,jj,kk,iproc0,jproc0,ip
      integer ixfst(irz),ixlst(irz),jzfst(irz),jzlst(irz)
      integer i,iin,jjn,nxin,nyin,nzin,a,sty,eny,szy
      integer  xistin,xiszin,zjstin,zjszin,iaux
      integer xp1,xp2,jz,slabpenc          
      integer dmy(3),hms(3),date_time(8)
      logical casex1,casex2,casez1,casez2          
      logical :: opn         
      character*20 string
      character*100 :: nam
      character*6 numer
c suggestion from D. Whitaker: allow fn to be as long as 256 characters
      character*256 fn,fn0
      character*11 unform
	logical exs, exs0
c
	integer ip2,ip3,luout
	integer ndirs0
	character*30 fnout
	real*8 rtime1,rtime2,cpuread
	integer relay
	integer itask
	real term,term1

	integer xi1,xi2,zj1,zj2
c
      data unform/'unformatted'/
c

	integer ndirs,jwrite,nchar_pid
	character*4 cpid

        integer ind_prec
c
c this line added 4/11/12
	integer jxp1,jxp2


	if(taskid.eq.0) write(6,*) 'enter inpen_read_bw'

	luout=6
	if (itask.eq.taskid.and.relay.ne.numtasks) then
	if (relay.gt.0) luout=60
	else
	if (itask.ge.iproc) luout=60
	end if

c
c
c	xi1=xist
	xi2=xi1+xisz-1
c	zj1=zjst
	zj2=zj1+zjsz-1
c
	rtime1=MPI_WTIME()

	if (relay.gt.0) then
c	ndirs0=2**(ip2(numtasks)/2)
	ndirs0=2**(ip3(numtasks)/2)
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
      
      string='enter inpen_read'
      write (luout,610) string(1:20),taskid,dmy(2),dmy(1),dmy(3),hms
 610  format (a20,2x,i6, ' date & time is  ',i2,'/',i2,
     1     '/',i4,2x,i2,':',i2,':',i2)
c        call MPI_FINALIZE (MPI_COMM_WORLD,mpierr)
c        stop




      do 100 kk=1,irz
         
         casex1=.false.
         casex2=.false.
         casez1=.false.
         casez2=.false.
         
         ip=kk-1

         if (xi2.ge.ixfst(kk).and.xi2.le.ixlst(kk)) casex1=.true.
         if (ixlst(kk).ge.xi1.and.ixlst(kk).le.xi2) casex2=.true.
         if (zj2.ge.jzfst(kk).and.zj2.le.jzlst(kk)) casez1=.true.
         if (jzlst(kk).ge.zj1.and.jzlst(kk).le.zj2) casez2=.true.

         if ((casex1.or.casex2).and.(casez1.or.casez2)) then
            write (luout,601) ip, taskid,xi1,xi2,zj1,zj2,
     1           ixfst(kk),ixlst(kk),jzfst(kk),jzlst(kk),
     1           casex1,casex2,casez1,casez2
 601        format ('inpen: ip', 2i6,1x,8i4,4L3)

	write (numer,"(i6)") ip
	call blanks (numer,nchar1)
            call blanks (indir_fn,nchar)
c
	if (ndirs.eq.0) then
c
            fn=indir_fn(1:nchar)//'/'//fninit(1)//'.p'//numer(1:nchar1)
            call blanks (fn,nchar)
c
	else

	jwrite=mod(ip,ndirs)
        write (cpid,"(i3,'/')") jwrite
        call blanks (cpid,nchar_pid)

            fn=indir_fn(1:nchar)//'/'//cpid(1:nchar_pid)
            call blanks (fn,nchar)
	    fn=fn(1:nchar)//fninit(1)//'.p'//numer(1:nchar1)
            call blanks (fn,nchar)
c
            fn0=indir_fn(1:nchar)//'/'//'0/'
            call blanks (fn0,nchar0)
	    fn0=fn0(1:nchar0)//fninit(1)//'.p'//numer(1:nchar1)
            call blanks (fn0,nchar0)
c
	end if
c
            write (luout,*) 'inpen:taskid,',taskid,fn(1:nchar)

            lu=luinit(1)
	inquire (file=fn(1:nchar),exist=exs)
	if (exs) then
	open (lu,file=fn(1:nchar),form='unformatted',action='read')
	else
	inquire (file=fn0(1:nchar0),exist=exs0)
	if (exs0) then
	open (lu,file=fn0(1:nchar0),form='unformatted',action='read')
	else
	write (6,*) 'cannot locate file',fn(1:nchar)
	end if
	end if


c            call fopen1 (lu,fn(1:nchar),unform)
            do 210 i=1,3+nc                   
               if (kinit(i).gt.0) then                      
                  iin=kinit(i)                      
 215		continue
                  read( lu , err=10 , end=20 ) jjn , nxin , nyin , nzin

                  if (slabpenc.eq.2) then
                     read( lu , err=10 , end=20 )
     1                    xistin, xiszin,zjstin,zjszin,iaux,nyin
                  else
		     xistin=1
                     xiszin=nxin/2
                  end if

                  allocate (unpen(xiszin,ny))

                if (ind_prec.eq.1) then
                  allocate (unpen_s_r(xiszin,ny))
                  allocate (unpen_s_i(xiszin,ny))
		unpen_s_r=0.
		unpen_s_i=0.
                else
                  allocate (unpen_d_r(xiszin,ny))
                  allocate (unpen_d_i(xiszin,ny))
		unpen_d_r=0.
		unpen_d_i=0.
                end if

                  do 220 jz=jzfst(kk),jzlst(kk)

c                    if( jz .eq. nzhp ) go to 220
                     if( jz .eq. nz-nzin/2+1 ) go to 220

                     tread_io = tread_io - MPI_WTIME()

                if (ind_prec.eq.1) then

                     if (slabpenc.eq.2) then

	if(i.eq.2) then
          read (lu,err=30,end=40) (unpen_s_r(xp,1),unpen_s_i(xp,1),xp=1,xiszin)
	else
          read (lu,err=30,end=40) ((unpen_s_r(xp,y),unpen_s_i(xp,y),xp=1,xiszin),y=1,nyin/2)
          read (lu,err=30,end=40) ((unpen_s_r(xp,y),unpen_s_i(xp,y),xp=1,xiszin),y=ny-nyin/2+2,ny)

	endif
                     else




                        do 231 y=1,ny
                           if (y.gt.nyin/2.and.y.lt.ny-nyin/2+2) go to 231
                           read (lu,err=30,end=40) (unpen_s_r(xp,y),unpen_s_i(xp,y),xp=1,xiszin)
 231                    continue

                     end if
                unpen(:,:)=unpen_s_r(:,:)+unpen_s_i(:,:)*(0.,1.)

                else

                     if (slabpenc.eq.2) then
	if (i.eq.2) then
          read (lu,err=30,end=40) (unpen_d_r(xp,1),unpen_d_i(xp,1),xp=1,xiszin)
	else
          read (lu,err=30,end=40) ((unpen_d_r(xp,y),unpen_d_i(xp,y),xp=1,xiszin),y=1,nyin/2)
          read (lu,err=30,end=40) ((unpen_d_r(xp,y),unpen_d_i(xp,y),xp=1,xiszin),y=ny-nyin/2+2,ny)
	endif
                     else




                        do 232 y=1,ny
                           if (y.gt.nyin/2.and.y.lt.ny-nyin/2+2) go to 232
                           read (lu,err=30,end=40) (unpen_d_r(xp,y),unpen_d_i(xp,y),xp=1,xiszin)
 232                    continue

                     end if
                unpen(:,:)=unpen_d_r(:,:)+unpen_d_i(:,:)*(0.,1.)
c
                end if


                     tread_io = tread_io + MPI_WTIME()

                     if (jz.lt.zj1.or.jz.gt.zj2) go to 220

                     zp=jz-zj1+1

c ZHAI: 08/09/2016
                if(mod(nx,nxin).eq.0) then
                        xp1=ixfst(kk)-xi1+1
                        xp2=xp1+xiszin-1
                        jxp1=1
                        jxp2=xiszin
                else
c say from 256^3 to 384^3
c assume iproc=iproc0, jproc=jproc0
                if(xi1.eq.ixfst(kk)) then
                xp1=1
                xp2=xiszin
                jxp1=1
                jxp2=xiszin
                else if(xi2.eq.ixlst(kk)) then
                xp2=xisz
                xp1=xp2-xiszin+1
                jxp1=1
                jxp2=xiszin
                else
                        if(casex1.and..not.casex2) then
                        !
                        jxp1=1
                        jxp2=jxp1+(xien-ixfst(kk))
                        xp2=xisz
                        xp1=xp2-(jxp2-jxp1)
                        else if(.not.casex1.and.casex2) then
                        !
                        jxp2=xiszin
                        jxp1=jxp2-(ixlst(kk)-xist)
                        xp1=mod(xist,xisz)
                        xp2=xp1+(jxp2-jxp1)
                        end if
                end if
                end if



		if (slabpenc.eq.1) then
		xp1=1
		xp2=xisz
		jxp1=xi1
		jxp2=xi1+xisz-1
		end if
c
c fix on 5/28/2012 for cases where iproc > iproc0
c
	if (iproc.gt.iproc0) then
	xp1=1
	xp2=xisz
	jxp1=xi1-xistin+1
	jxp2=jxp1+xisz-1
	end if
c
	if (xp1.lt.1.or.xp2.gt.xisz) then
	write (luout,"('taskid,ipid,xi1,xistin,xisz,xp1,xp2=',7i5)") taskid,ipid,xi1,xistin,xisz,xp1,xp2
	call abrt ('inpen_read: illegal xp1 or xp2')
	end if

       if (jxp1.lt.1.or.jxp2.gt.xiszin) then
	write (luout,"('taskid,ipid,xi1,xistin,xisz,jxp1,jxp2=',7i5)") taskid,ipid,xi1,xistin,xisz,jxp1,jxp2
	call abrt ('inpen_read: illegal jxp1 or jxp2')
	end if
	

!	if(i.eq.2) then
!		uny(1,zp,xp1:xp2,i) = unpen(jxp1:jxp2,1)
!	else

                        do y=1,nyin/2
                           uny(y,zp,xp1:xp2,i)=unpen(jxp1:jxp2,y)
                        end do
                        do y=ny-nyin/2+2,ny
                           uny(y,zp,xp1:xp2,i)=unpen(jxp1:jxp2,y)
                        end do
!	endif


 220              continue

                  deallocate (unpen)
                  if (ind_prec.eq.1) then
                deallocate (unpen_s_r)
                deallocate (unpen_s_i)
                else
                deallocate (unpen_d_r)
                deallocate (unpen_d_i)
                end if
c
	if (i.ge.4.and.jjn.lt.kinit(i)) go to 215
c
               end if

 210        continue



	do x=1,xisz
	xp=x+xist-1
	do z=1,zjsz
	zp = z+zjst-1
	do y=2,ny
	uny(y,z,x,2) = -(kx(xp)*uny(y,z,x,1)*b11(2) + kz(zp)*uny(y,z,x,3)*b33(2))/ky(y)/b22(2)
	enddo
	enddo
	enddo

         close (lu)

         end if

 100  continue
	



c Now redistribute the data based on cylindrical data allocation

c Lines originally here have been moved to inpen.f, to be executed
c after all 1 tasks have read files assigned to them
c 

c     call MPI_Barrier(MPI_COMM_WORLD,ierr)
c     if(taskid .eq. 0) then
	rtime2=MPI_WTIME()
	cpuread=rtime2-rtime1
      call date_and_time (values=date_time)
      dmy(1:3)=date_time(3:1:-1)
      hms(1:3)=date_time(5:7)
         string=' exit inpen_read'
         write (luout,620) string(1:20),taskid,dmy(2),dmy(1),dmy(3),hms,cpuread,ipid,jpid
	if (cpuread.gt.100.) then 
         write (6,620) string(1:20),taskid,dmy(2),dmy(1),dmy(3),hms,cpuread,ipid,jpid
	end if
 620  format (a20,2x,i6, ' date & time is  ',i2,'/',i2,
     1     '/',i4,2x,i2,':',i2,':',i2,f6.1,' secs',i3,i5)
	
c     endif

      inquire (60, opened=opn, name=nam)
	if (opn) close (60)
	
c	write (6,*) 'inpen_read_bw: itask,uny(2,2,2,1)=',itask,uny(2,2,2,1)

 90   return

 10   write(6,*)' error reading field header: jjn=',jjn
      write(6,*)' stopped in inpen, taskid= ',taskid,fn(1:nchar)
      inquire (lu, opened=opn, name=nam)
      write(6,*)' reading lu,filename=',lu,nam
      go to 99

 20   write(6,*)' end of file hit reading field header: jjn=',jjn
      write(6,*)' stopped in inpen, taskid= ',taskid,fn(1:nchar)
      inquire (lu, opened=opn, name=nam)
      write(6,*)' reading lu,filename=',lu,nam
      go to 99

 30   write(6,*)' error reading field:i,iin',i,iin,taskid
      write(6,*)' stopped in inpen, taskid= ',taskid,fn(1:nchar)
      inquire (lu, opened=opn, name=nam)
      write(6,*)' reading lu,filename=',lu,nam
      go to 99

 40   write(6,*)' end of file hit reading field:i,iin',i,iin,y,xp,jz
      write(6,*)' stopped in inpen, taskid= ',taskid,fn(1:nchar)
      inquire (lu, opened=opn, name=nam)
      write(6,*)' reading lu,filename=',lu,nam
      go to 99

 99   call abrt('from inpen')

      return
      end subroutine inpen_read_bw
