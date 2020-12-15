      subroutine wrtpos (lu)
c
c Revised by PKY. 3/10/2014
c
c save particle positions on unformatted disk file
c (for purpose of checkpointing)
c
c --- MPI version 
c
c task 0 for fluid particles
c tasks 2 and numtasks-3 for molecules
c
#ifdef LAG
c
	use compart
	implicit none
c
	integer nn,numc,npr,indx
	integer i,i1,i2,k,nr,lu,j,mp
	character*12 caux
	character*2 numer 
c
	real(p8), allocatable :: ppmtmp(:,:)
	real(p8), allocatable :: pptmp(:)
	real(p8), allocatable :: pp2tmp(:)
	integer, allocatable :: pp2indtmp(:)
	integer, allocatable :: ppindtmp(:)
        integer, allocatable :: num_all(:),idisp(:)
	integer itask
c
	integer ip,igm
c
      character*75 string
c
	character*11 form,unform
      data form,unform/'formatted','unformatted'/

c npr = no of particles written per record
c
!      data npr/100/
      data npr/1024/
c
	if (taskid.eq.0) write (6,*) 'enter wrtpos'

	if (nop.gt.0) then

	caux='savepos'
      if (taskid.eq.0) call fopen1 (lu,caux,unform)
c
c
	if (taskid.eq.numtasks-1) then
        write (numer,600) ichkpt
 600    format (i2)
        numc=1+floor(alog10(1.*ichkpt))
        caux='savepos.'//numer(3-numc:2)
        call blanks (caux,nn)
      call fopen1 (lu,caux(1:nn),unform)
	end if
c
c use pptmp as scratch global array, since (see partsp.f)
c new interpolation weights are to be calculated next
c
c
        allocate (num_all(numtasks),idisp(numtasks),stat=ierr)

        call MPI_ALLGATHER (nump1,1,MPI_INTEGER,
     1                       num_all,1,MPI_INTEGER,MPI_COMM_WORLD,mpierr)

        idisp(1) = 0
        do itask=2,numtasks
        idisp(itask) = idisp(itask-1) + num_all(itask-1)
        enddo



      if (taskid.eq.0.or.taskid.eq.numtasks-1) then  
c
	npr=min0(nop,npr)
        write (6,*) 'wrtpos: nop,npr=',nop,npr
      write (lu) nop,npr
c
      write (lu) xyzl
      write (lu) gx(1),gy(1),gz(1)
        write(lu) numtasks
        write(lu) num_all(1:numtasks)
c
c determine no. of records (nr)
c
      nr=nop/npr
      if (mod(nop,npr).ne.0) nr=nr+1
c
	endif

	allocate (ppindtmp(nop))


        call MPI_ALLGATHERV (pp1ind(1),nump1,MPI_INTEGER,ppindtmp(1),
     1                       num_all,idisp,MPI_INTEGER,MPI_COMM_WORLD,mpierr)



      if (taskid.eq.0.or.taskid.eq.numtasks-1) then

      do 10 k=1,nr
      i1=(k-1)*npr+1
      i2=min0(i1+npr-1,nop)
      write (lu) (ppindtmp(i),i=i1,i2)
 10   continue

	endif

	deallocate (ppindtmp)
	allocate (pptmp(nop))

        do j=1,3

        call MPI_ALLGATHERV (pp(1,j),nump1,pptype,
     1                       pptmp(1),num_all,idisp,pptype,MPI_COMM_WORLD,mpierr)


	if (taskid.eq.0.or.taskid.eq.numtasks-1) then

      do 11 k=1,nr
      i1=(k-1)*npr+1
      i2=min0(i1+npr-1,nop)
      write (lu) (pptmp(i),i=i1,i2)
 11   continue

	endif

	enddo

        if (taskid.eq.0.or.taskid.eq.numtasks-1) then 
	call fclose1 (lu)
	write (6,*) 'particle positions saved, step=',istep-1,pptmp(4)
	endif
c
	deallocate (pptmp)
c
	endif ! if (nop.gt.0) 
c

	call MPI_BARRIER (MPI_COMM_WORLD, mpierr)

#ifdef MOL
c
	if (nom.gt.0) then
c
        caux='savempos'
      if (taskid.eq.2) call fopen1 (lu,caux,unform)
c
        if (taskid.eq.numtasks-3) then
        write (numer,600) ichkpt
        numc=1+floor(alog10(1.*ichkpt))
        caux='savempos.'//numer(3-numc:2)
        call blanks (caux,nn)
      call fopen1 (lu,caux(1:nn),unform)
        end if
c
      if (taskid.eq.2.or.taskid.eq.numtasks-3) then
      write (lu) nom,npr
      write (lu) xyzl
      write (lu) gx(1),gy(1),gz(1)
      end if
c
        allocate (ppmtmp(nom/ngm,3))
c
        do igm=1,ngm
c
        do j=1,3
        call MPI_ALLGATHER (mpp(1,j,igm),nm/ngm,pptype,
     1                      ppmtmp(1,j),nm/ngm,pptype,
     1                      MPI_COMM_WORLD,mpierr)
        end do
c
      if (taskid.eq.2.or.taskid.eq.numtasks-3) then
        do j=1,3
      write (lu) ppmtmp(:,j)
        end do
        end if
c
        end do
c
      if (taskid.eq.2.or.taskid.eq.numtasks-3) then
      call fclose1 (lu)
      write (6,*) 'molecular positions saved, step=',istep-1
        end if
c
	deallocate (ppmtmp)

      end if

#endif

#ifdef LAGSC2

	if (nop2.gt.0) then
c
	if (taskid.eq.1) then
        caux='save2pos'
        call fopen1 (lu,caux,unform)
	endif
c
	if (taskid.eq.numtasks-2) then
        write (numer,600) ichkpt
        numc=1+floor(alog10(1.*ichkpt))
        caux='save2pos.'//numer(3-numc:2)
        call blanks (caux,nn)
        call fopen1 (lu,caux(1:nn),unform)
	end if
c
c
! checkpointing for new algorithm

	allocate (pp2indtmp(nop2))
!        allocate (num_all(numtasks),idisp(numtasks),stat=ierr)

 	call MPI_ALLGATHER (nump2,1,MPI_INTEGER,
     1                       num_all,1,MPI_INTEGER,MPI_COMM_WORLD,mpierr)


	idisp(1) = 0
	do itask=2,numtasks
	idisp(itask) = idisp(itask-1) + num_all(itask-1)
	enddo


      if (taskid.eq.1.or.taskid.eq.numtasks-2) then
	npr=min0(nop2,npr)
      write (lu) nop2,npr
c
      write (lu) xyzl
      write (lu) gx(1),gy(1),gz(1)
	write(lu) numtasks
	write(lu) num_all(1:numtasks)

! determine no. of records (nr)
      nr=nop2/npr
      if (mod(nop2,npr).ne.0) nr=nr+1

	endif
c

 	call MPI_ALLGATHERV (pp2ind(1),nump2,MPI_INTEGER,pp2indtmp(1),
     1                       num_all,idisp,MPI_INTEGER,MPI_COMM_WORLD,mpierr)

      if (taskid.eq.1.or.taskid.eq.numtasks-2) then
      do k=1,nr
      i1=(k-1)*npr+1
      i2=min0(i1+npr-1,nop2)
      write (lu) (pp2indtmp(i),i=i1,i2)
      enddo
      end if

	deallocate (pp2indtmp)
	allocate (pp2tmp(nop2))


	do j=1,3

 	call MPI_ALLGATHERV (pp2(1,j),nump2,pptype,
     1                       pp2tmp(1),num_all,idisp,pptype,MPI_COMM_WORLD,mpierr)

      if (taskid.eq.1.or.taskid.eq.numtasks-2) then
c
!	do ip=1,nop2
!	indx = int(ppindx(ip))
!	if (ppindx(ip)-1.d0*indx.gt.1.e-7) stop
!	ppt2 (indx) = pp2tmp(ip)
!	enddo

      do k=1,nr
      i1=(k-1)*npr+1
      i2=min0(i1+npr-1,nop2)
      write (lu) (pp2tmp(i),i=i1,i2)
	enddo

      end if

	enddo
c
      if (taskid.eq.1.or.taskid.eq.numtasks-2) call fclose1 (lu)

	deallocate (pp2tmp)


	endif
c
#endif

#ifdef MATP
c
      if (taskid.eq.2) then
c
	caux='savempos'
      call fopen1 (lu,caux,unform)
c
      write (lu) nomp,npr
c
      write (lu) xyzl
      write (lu) gx(1),gy(1),gz(1)
c
c determine no. of records (nr)
c
      nr=nomp/npr
      if (mod(nomp,npr).ne.0) nr=nr+1
c
      do 40 k=1,nr
      i1=(k-1)*npr+1
      i2=min0(i1+npr-1,nomp)
      write (lu) (mppv(i,1),mppv(i,2),mppv(i,3),i=i1,i2)
 40   continue
      write (6,*) 'material particles written to savempos:',
     1             mppv(2,2),mppv(3,3),mppv(1,1)
c
      call fclose1 (lu)
c
      end if
c
#endif


	if (taskid.eq.0) write (6,*) ' exit wrtpos'
#endif
      return
      end
