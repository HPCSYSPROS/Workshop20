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
	integer nn,numc,npr
	integer i,i1,i2,k,nr,lu,j,mp
	character*12 caux
	character*2 numer 
c
	real(p8), allocatable :: pptmp(:,:)
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
c      data npr/100/
      data npr/1024/
c
	if (taskid.eq.0) write (6,*) 'enter wrtpos'
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
	mp=nop/numtasks
c
	allocate (pptmp(nop,3))
	do j=1,3
	call MPI_ALLGATHER (pp(1,j),mp,pptype,pptmp(1,j),mp,pptype,
     1                      MPI_COMM_WORLD,mpierr)
	end do
c

      if (taskid.eq.0.or.taskid.eq.numtasks-1) then  
c
	npr=min0(nop,npr)
        write (6,*) 'wrtpos: nop,npr=',nop,npr
      write (lu) nop,npr
c
      write (lu) xyzl
      write (lu) gx(1),gy(1),gz(1)
c
c determine no. of records (nr)
c
      nr=nop/npr
      if (mod(nop,npr).ne.0) nr=nr+1
c
      do 10 k=1,nr
      i1=(k-1)*npr+1
      i2=min0(i1+npr-1,nop)
      write (lu) (pptmp(i,1),pptmp(i,2),pptmp(i,3),i=i1,i2)
 10   continue
c
      call fclose1 (lu)
      write (6,*) 'particle positions saved, step=',istep-1,pptmp(4,1)
c
      end if
c
c
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
        deallocate (pptmp)
        allocate (pptmp(nom/ngm,3))
c
        do igm=1,ngm
c
        do j=1,3
        call MPI_ALLGATHER (mpp(1,j,igm),nm/ngm,pptype,
     1                      pptmp(1,j),nm/ngm,pptype,
     1                      MPI_COMM_WORLD,mpierr)
        end do
c
      if (taskid.eq.2.or.taskid.eq.numtasks-3) then
        do j=1,3
      write (lu) pptmp(:,j)
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
      end if
#endif

#ifdef LAGSC2
c
        caux='save2pos'
      if (taskid.eq.1) call fopen1 (lu,caux,unform)
c
	if (taskid.eq.numtasks-2) then
        write (numer,600) ichkpt
        numc=1+floor(alog10(1.*ichkpt))
        caux='save2pos.'//numer(3-numc:2)
        call blanks (caux,nn)
      call fopen1 (lu,caux(1:nn),unform)
	end if
c
	deallocate (pptmp)
	allocate (pptmp(nop2,3))
c
#ifndef SEP1912
	mp=nop2/numtasks
	do j=1,3
 	call MPI_ALLGATHER (pp2(1,j),mp,pptype,pptmp(1,j),mp,pptype,
     1                      MPI_COMM_WORLD,mpierr)
	end do
#else
! Warning: the SEP1912 directive is not tested
        allocate (num_all(numtasks),idisp(numtasks),stat=ierr)
	mp=nop2/numtasks
        num_all(:)=mp/2
        do itask=1,numtasks
        idisp(itask)=(itask-1)*mp
        end do
	do j=1,3
 	call MPI_ALLGATHERV (pp2(1,j),mp/2,pptype,
     1                       pptmp(1,j),num_all,idisp,pptype,MPI_COMM_WORLD,mpierr)
 	call MPI_ALLGATHERV (pp2(mp/2+1,j),mp/2,pptype,
     1                       pptmp(mp/2+1,j),num_all,idisp,pptype,MPI_COMM_WORLD,mpierr)
	end do
#endif

#ifdef TEST
	if (taskid.eq.0) then
	do ip=1,nop2
	write (801,"('ip,pp2=',i6,1p,3e12.4)") ip,(pptmp2(ip,j),j=1,3)
	end do
	end if
c
	do ip=1,nop2/numtasks
	write (1000+taskid,"('ip,pp2=',i6,1p,3e12.4)") ip,(pp2(ip,j),j=1,3)
	end do
#endif
c
      if (taskid.eq.1.or.taskid.eq.numtasks-2) then
c
c
	npr=min0(nop2,npr)
      write (lu) nop2,npr
c
      write (lu) xyzl
      write (lu) gx(1),gy(1),gz(1)
c
c determine no. of records (nr)
c
      nr=nop2/npr
      if (mod(nop2,npr).ne.0) nr=nr+1
c
      do 30 k=1,nr
      i1=(k-1)*npr+1
      i2=min0(i1+npr-1,nop2)
      write (lu) (pptmp(i,1),pptmp(i,2),pptmp(i,3),i=i1,i2)
 30   continue
c
      call fclose1 (lu)
c
      end if
c
	deallocate (pptmp)
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
