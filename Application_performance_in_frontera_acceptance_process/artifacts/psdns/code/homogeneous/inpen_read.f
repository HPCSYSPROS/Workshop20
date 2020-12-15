!     NPM converting to stride-1 3/9/08
      subroutine inpen_read  (uny,slabpenc,iproc0,ixfst,ixlst,jzfst,jzlst,ndirs,ind_prec)
      use comp
      use timers_io
      implicit none
      include 'intvars'
      
c It doesn't really matter how we declare uny here; 
c Upon exit, uny will contain data in cylindrical form, that is 
c (ny,zjsz*xisz,3+nc)

#ifdef INPUT_PEN
      complex(b8) :: uny(ny,zjsz,xisz,3+nc)
        complex(4), allocatable :: unpen_s(:,:)
        complex(8), allocatable :: unpen_d(:,:)
#else
      complex(b8) :: uny(ny,zjsz*xisz,3+nc)
#endif
      complex(b8), allocatable :: unpen(:,:),buf(:,:,:),utemp(:,:)

      integer nchar, nchar1
      integer lu,ii,jj,kk,iproc0,jproc0,ip
      integer ixfst(irz),ixlst(irz),jzfst(irz),jzlst(irz)
      integer i,iin,jjn,nxin,nyin,nzin,a,sty,eny,szy
      integer  xistin,xiszin,zjstin,zjszin,iaux
      integer xp1,xp2,jz,slabpenc          
      integer dmy(3),hms(3),date_time(8)
      logical casex1,casex2,casez1,casez2          
      logical :: opn         
      character*20 string
      character*80 :: nam
      character*6 numer
c suggestion from D. Whitaker: allow fn to be as long as 256 characters
      character*256 fn
      character*11 unform

	character*20 fnout

      data unform/'unformatted'/
c
#ifdef SUBDIR
	integer ndirs,jwrite,nchar_pid
	character*4 cpid
#endif
        integer ind_prec
c
c this line added 4/11/12
	integer jxp1,jxp2
      
      call date_and_time (values=date_time)
      dmy(1:3)=date_time(3:1:-1)
      hms(1:3)=date_time(5:7)
      
	jwrite=mod(taskid,ndirs)
        write (cpid,"(i3,'/')") jwrite
        call blanks (cpid,nchar_pid)

	write (numer,"(i6)") taskid
	call blanks (numer,nchar1)
            fnout='readtimes'//'/'//cpid(1:nchar_pid)
            call blanks (fnout,nchar)
	    fnout=fnout(1:nchar)//'time.'//numer(1:nchar1)
            call blanks (fnout,nchar)
	open (60,file=fnout(1:nchar))
	write (6,*) 'taskid: fnoutout=',taskid,fnout(1:nchar)
c
      string='enter inpen_read'
      write (60,610) string(1:20),taskid,dmy(2),dmy(1),dmy(3),hms
 610  format (a20,2x,i6, ' date & time is  ',i2,'/',i2,
     1     '/',i4,2x,i2,':',i2,':',i2)
	call MPI_FINALIZE (mpierr)
	stop


#ifdef INPUT_PEN

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
            write (60,601) ip, taskid,xist,xien,zjst,zjen,
     1           ixfst(kk),ixlst(kk),jzfst(kk),jzlst(kk),
     1           casex1,casex2,casez1,casez2
 601        format ('inpen: ip', i6,1x,9i4,4L3)

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
	end if
c
            write (60,*) 'inpen:taskid,',taskid,fn(1:nchar)

            lu=luinit(1)
            call fopen1 (lu,fn(1:nchar),unform)
            do 210 i=1,3+nc                   
               if (kinit(i).gt.0) then                      
                  iin=kinit(i)                      
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
                  allocate (unpen_s(xiszin,ny))
                else
                  allocate (unpen_d(xiszin,ny))
                end if

                  do 220 jz=jzfst(kk),jzlst(kk)

c                    if( jz .eq. nzhp ) go to 220
                     if( jz .eq. nz-nzin/2+1 ) go to 220

                     tread_io = tread_io - MPI_WTIME()

                if (ind_prec.eq.1) then

                     if (slabpenc.eq.2) then
          read (lu,err=30,end=40) ((unpen_s(xp,y),xp=1,xiszin),y=1,nyin/2)
          read (lu,err=30,end=40) ((unpen_s(xp,y),xp=1,xiszin),y=ny-nyin/2+2,ny)
                     else
#ifdef BLOCK
          read (lu,err=30,end=40) ((unpen_s(xp,y),xp=1,xiszin),y=1,nyin/2)
          read (lu,err=30,end=40) ((unpen_s(xp,y),xp=1,xiszin),y=ny-nyin/2+2,ny)
#else
                        do 231 y=1,ny
                           if (y.gt.nyin/2.and.y.lt.ny-nyin/2+2) go to 231
                           read (lu,err=30,end=40) (unpen_s(xp,y),xp=1,xiszin)
 231                    continue
#endif
                     end if
                unpen(:,:)=unpen_s(:,:)

                else

                     if (slabpenc.eq.2) then
          read (lu,err=30,end=40) ((unpen_d(xp,y),xp=1,xiszin),y=1,nyin/2)
          read (lu,err=30,end=40) ((unpen_d(xp,y),xp=1,xiszin),y=ny-nyin/2+2,ny)
                     else
#ifdef BLOCK
          read (lu,err=30,end=40) ((unpen_d(xp,y),xp=1,xiszin),y=1,nyin/2)
          read (lu,err=30,end=40) ((unpen_d(xp,y),xp=1,xiszin),y=ny-nyin/2+2,ny)
#else
                        do 232 y=1,ny
                           if (y.gt.nyin/2.and.y.lt.ny-nyin/2+2) go to 232
                           read (lu,err=30,end=40) (unpen_d(xp,y),xp=1,xiszin)
 232                    continue
#endif
                     end if
                unpen(:,:)=unpen_d(:,:)
c
                end if


                     tread_io = tread_io + MPI_WTIME()

                     if (jz.lt.zjst.or.jz.gt.zjen) go to 220

                     zp=jz-zjst+1


                        xp1=ixfst(kk)-xist+1
                        xp2=xp1+xiszin-1
			jxp1=1
			jxp2=xiszin
c
		if (slabpenc.eq.1) then
		xp1=1
		xp2=xisz
		jxp1=xist
		jxp2=xist+xisz-1
		end if
c
c fix on 5/28/2012 for cases where iproc > iproc0
c
	if (iproc.gt.iproc0) then
	xp1=1
	xp2=xisz
	jxp1=xist-xistin+1
	jxp2=jxp1+xisz-1
	end if
c
	if (xp1.lt.1.or.xp2.gt.xisz) then
	write (60,"('taskid,ipid,xist,xistin,xisz,xp1,xp2=',i6,6i5)") taskid,ipid,xist,xistin,xisz,xp1,xp2
	call abrt ('inpen_read: illegal xp1 or xp2')
	end if

       if (jxp1.lt.1.or.jxp2.gt.xiszin) then
	write (60,"('taskid,ipid,xist,xistin,xisz,jxp1,jxp2=',i6,6i5)") taskid,ipid,xist,xistin,xisz,jxp1,jxp2
	call abrt ('inpen_read: illegal jxp1 or jxp2')
	end if
	

                        do y=1,nyin/2
                           uny(y,zp,xp1:xp2,i)=unpen(jxp1:jxp2,y)
                        end do
                        do y=ny-nyin/2+2,ny
                           uny(y,zp,xp1:xp2,i)=unpen(jxp1:jxp2,y)
                        end do

                     
 220              continue

                  deallocate (unpen)
                  if (ind_prec.eq.1) then
                deallocate (unpen_s)
                else
                deallocate (unpen_d)
                end if
c
               end if

 210        continue

         close (lu)

         end if

 100  continue

#ifdef TEST_CYLIO
      ucomp = uny
#endif

c Now redistribute the data based on cylindrical data allocation

c Lines originally here have been moved to inpen.f, to be executed
c after all MPI tasks have read files assigned to them
c 
#else

      do i=1,3+nc                   
         if (kinit(i).gt.0) then                      
            iin=kinit(i)                      
            read( lu , err=10 , end=20 ) jjn , nxin , nyin , nzin                      
            read( lu , err=10 , end=20 )
     1           xistin,max_al_xin,mystart_xin,mystart_zin,num_alin,nyin 
c     1           xistin, xiszin,zjstin,zjszin,iaux,nyin

            do a=1,num_al
               do y=1,ny
                  if(mask(y,a)) then
                     read(lu, err=30,end=40)  uny(y,a,i)
                  endif
               enddo
            enddo
         endif
      enddo
#endif

c     call MPI_Barrier(MPI_COMM_WORLD,ierr)
c     if(taskid .eq. 0) then
         string=' exit inpen_read'
         write (60,610) string(1:20),taskid,dmy(2),dmy(1),dmy(3),hms
c     endif

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
      end subroutine inpen_read
