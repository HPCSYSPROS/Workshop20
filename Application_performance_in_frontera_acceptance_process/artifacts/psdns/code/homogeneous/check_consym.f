	subroutine check_consym (uy,title)
c
      use comsp
c
      implicit none
	include 'intvars'

	complex(b8) :: uy(ny,xisz*zjsz)
c
	integer :: nzp2,nyp2,xst,k,ik,lu,nchar,a
	real :: term,rk2

	character(*) title
        character*30 caux
	integer itask
	integer, allocatable :: startx_all(:),startz_all(:)
	

	allocate (startx_all(0:numtasks-1))
	allocate (startz_all(0:numtasks-1))
	kx2(:)=kx(:)**2
	ky2(:)=ky(:)**2
	kz2(:)=kz(:)**2

	ky2(nyhp)=0.
	kz2(nzhp)=0.
c
        caux=title
        call blanks (caux,nchar)
c
        xp=mystart_x
        z=mystart_z
	if (taskid.le.1) write (801,*) taskid,mystart_x,mystart_z
        call MPI_GATHER (mystart_x,1,MPI_INTEGER,startx_all,1,
     1                  MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
        call MPI_GATHER (mystart_z,1,MPI_INTEGER,startz_all,1,
     1                  MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
	do itask=0,numtasks-1
	write (801,"('taskid,mystart_x,mystart_z',3i6)") itask,
     1             startx_all(itask),startz_all(itask)
	end do

        do 10 a=1,num_al
        x=xp+xist-1

        do 11 y=1,ny

        if(.not.mask(y,a)) go to 11

 	rk2=kx2(x)+ky2(y)+kz2(z)
 	ik=sqrt(rk2)+1.5
	if (ik.gt.nxh) go to 11
c	
        if (ik.le.2.and.x.eq.1) then
        write (6,221) caux(1:nchar),z,y,int(kz(z)),int(ky(y)),uy(y,a),taskid
221     format ('check_consym: ',a10,4i4,1p,2e16.8,i4)
        end if 
        
c
 11	continue		
        call next_xz(xp,z)
 10	continue
c
	deallocate (startx_all,startz_all)
c
	return
	end
