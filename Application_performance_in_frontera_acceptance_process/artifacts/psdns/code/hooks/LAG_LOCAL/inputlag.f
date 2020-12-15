      subroutine inputlag
c
#ifdef LAG
	use compart
c
	implicit none
	integer lu,igp,ic,igm
	integer buf_int(50),i1,i2,n_int
	real dxp1
c
	logical iex
c
c
	lu=1112

c
	nop=0
	nop2=0
c
c if input.lag file does not exist, set nop and nop2 to 0
c (which would cause all other routines to skip Lagrangian action)
c
	if (taskid.eq.1) then
	write (6,*) 'inquiring'
	inquire (file='input.lag',exist=iex)
	end if
c
	call MPI_BCAST (iex,1,MPI_LOGICAL,1,MPI_COMM_WORLD,mpierr)
	if (taskid.eq.0) write (6,*) 'inputlag: after bcast 1'

	if (.not.iex) return
c
	if (taskid.eq.1) then
	
 	open(lu,file='input.lag')
c
	read (lu,501) 
	read (lu,*) pstart, pstart2

c note true nop is (1+3*ngp)*nop if pstart=-9
c
	read (lu,501) 
	read (lu,*) nop
	write (6,*) 'pstart,nop=',pstart,nop

	read (lu,501) 
	read (lu,*) ngp
c
	buf_int(1)=pstart
	buf_int(2)=nop
	buf_int(3)=ngp
c
	end if
c
	call MPI_BCAST (buf_int,3,MPI_INTEGER,1,MPI_COMM_WORLD,mpierr)
	call MPI_BCAST (pstart2, 1, MPI_INTEGER, 1, MPI_COMM_WORLD, mpierr)
	if (taskid.eq.0) write (6,*) 'inputlag: after bcast 2'
c
	pstart=buf_int(1)
	nop=buf_int(2)
	ngp=buf_int(3)
c
	if (ngp.gt.0) allocate (dxpr(ngp),dypr(ngp),dzpr(ngp))
c
c
	if (taskid.eq.1) then
c
	if (ngp.gt.0) then
	read (lu,501) 
	read (lu,*) dxp1
	read (lu,501) 
	read (lu,*) (dxpr(igp),igp=1,ngp)
	dxpr(:) = dxpr(:)*dxp1
	else
	read (lu,501) 
	read (lu,501) 
	read (lu,501) 
	read (lu,501) 
	end if
	dypr(:)=0.
	dzpr(:)=0.
c
	read (lu,501) 
	read (lu,*) itetra
c
	read (lu,501) 
	read (lu,*) lpgrad
c
c read flags for Lagrangian scalar time series
c
#ifdef NOSCALAR
	read (lu,501) 
	read (lu,501) 
#else
	read (lu,501) 
	read (lu,*) (lpfunc(ic+3),ic=1,nc)
#endif

	read (lu,501) 
	read (lu,*) lpsdis
#ifndef LGRAD
	lpsdis(:)=0
#endif
c
#ifdef LAGSC2
	read (lu,501)
	read (lu,*) nop2
	ndpart2=nop2
#endif
c
	nsubset=1
	read (lu,501)
	read (lu,*,end=6) nsubset
 6	continue
c
	numbatch=0
	pim_flag=1
	read (lu,501,end=7)
	read (lu,*) ibflag,pim_flag
 7	continue
c
	read (lu,501,end=8)
	read (lu,*) nrec1,nrec2,nmrec
 8	continue

	read (lu,501,end=9)
	read (lu,"(a120)") indir_lag
 9	continue
c
 501	format (1x)

 99	close(lu)
c

c
	buf_int(4)=itetra
	i1=5
	i2=i1+3+nc-1
	buf_int(i1:i2)=lpfunc(1:3+nc)
	i1=i2+1
	i2=i1+3+nc-1
	buf_int(i1:i2)=lpgrad(1:3+nc)
	i1=i2+1
	i2=i1+3+nc-1
	buf_int(i1:i2)=lpsdis(1:3+nc)
	buf_int(i2+1)=nop2
	buf_int(i2+2)=ndpart
	buf_int(i2+3)=nsubset
	buf_int(i2+4)=ibflag
	buf_int(i2+5)=pim_flag
	buf_int(i2+6)=nrec1
	buf_int(i2+7)=nrec2
	buf_int(i2+8)=nmrec
c
	end if
c
        call MPI_BCAST(indir_lag,120,MPI_CHARACTER,1,MPI_COMM_WORLD,mpierr)

!	n_int=9+3*(3+nc)
	n_int=12+3*(3+nc)
	call MPI_BCAST (buf_int(4),n_int-2,MPI_INTEGER,1,MPI_COMM_WORLD,mpierr)
	if (taskid.eq.0) write (6,*) 'inputlag: after bcast 3'
c
	itetra=buf_int(4)
	i1=5
	i2=i1+3+nc-1
	lpfunc(1:3+nc)=buf_int(i1:i2)
	i1=i2+1
	i2=i1+3+nc-1
	lpgrad(1:3+nc)=buf_int(i1:i2)
	i1=i2+1
	i2=i1+3+nc-1
	lpsdis(1:3+nc)=buf_int(i1:i2)
	nop2=buf_int(i2+1)
	ndpart=buf_int(i2+2)
	nsubset=buf_int(i2+3)
	ibflag=buf_int(i2+4)
	pim_flag=buf_int(i2+5)
	nrec1=buf_int(i2+6)
	nrec2=buf_int(i2+7)
	nmrec=buf_int(i2+8)
c
	if (mod(nop,numtasks).ne.0.or.mod(nop2,numtasks).ne.0) then
	if (taskid.eq.0) then
	write (6,*) 'inputlag: nop and nop2 must be multiples of numtasks'
	end if
	call MPI_ABORT (MPI_COMM_WORLD,mpierr)
	end if
c
#ifdef LAGSC2
	ndpart2=nop2
#endif
c
	ndprop=6
#ifdef LGRAD
	if (any(lpgrad(1:3).gt.0)) ndprop=max0(ndprop,15)
#endif
c
c reset ndpart if necessary
c
	ndpart=nop
	if (pstart.eq.-9.or.pstart.eq.-10) then
	ndpart=(1+3*ngp)*nop
	end if

	if (ngp.gt.0) then
	call MPI_BCAST (dxpr,ngp,mpireal,1,MPI_COMM_WORLD,mpierr)
	if (taskid.eq.0) write (6,*) 'inputlag: after bcast 4'
	dypr(:)=0.
	dzpr(:)=0.
	end if
c
	if (ibflag.eq.1) numbatch=iproc
	if (ibflag.eq.2) numbatch=jproc
c

#ifdef MOL
	if (taskid.eq.1) then
	write(6,*) 'inquiring input.mol'
	inquire (file='input.mol',exist=iex)
	endif

c
	call MPI_BCAST (iex,1,MPI_LOGICAL,1,MPI_COMM_WORLD,mpierr)
	if (taskid.eq.0) write (6,*) 'inputlag: after bcast mol'

	if (.not.iex) then
	if (taskid.eq.0) write (6,*) 'no molecules: input.mol file absent'
	nompp=0
	ngm=0
	nom=0
	return
	end if

	if (taskid.eq.1) then

	open(lu,file='input.mol')
c
	read (lu,501)
	read (lu,*) nompp

	read (lu,501)
	read (lu,*) ngm

	endif

	call MPI_BCAST (nompp,1,MPI_INTEGER,1,MPI_COMM_WORLD,mpierr)
	call MPI_BCAST (ngm,1,MPI_INTEGER,1,MPI_COMM_WORLD,mpierr)

	allocate (molsc(ngm))

	if (taskid.eq.1) then
	read (lu,501)
	read (lu,*) (molsc(igm),igm=1,ngm)

	read (lu,501)
	read (lu,*) nmsubset
#ifndef MOLOLD
	read (lu,501)
	read (lu,*) mstart, minitseed
#endif
	endif

	call MPI_BCAST (molsc,ngm,mpireal,1,MPI_COMM_WORLD,mpierr)
	call MPI_BCAST (nmsubset,1,MPI_INTEGER,1,MPI_COMM_WORLD,mpierr)
#ifndef MOLOLD
	call MPI_BCAST (mstart,1,MPI_INTEGER,1,MPI_COMM_WORLD,mpierr)
	call MPI_BCAST (minitseed,1,MPI_INTEGER,1,MPI_COMM_WORLD,mpierr)
#endif

	if(taskid.eq.0) write(6,"('input.mol read:nompp,ngm,molsc=',2i4,1p,7e12.4)") nompp,ngm,molsc
#endif

c
#endif
      return
      end
c
