	subroutine compart_set 

! D Buaria, Feb 05 2015,
! modified declaration of bs and pp2 arrays
! to fit the new interpolation algorithm
c	
#ifdef LAG
c
	use compart
	implicit none
	integer nipj
	integer jj,i,j,id,jd,itask
	integer numpmax
c
c read the Lagrangian parameters
c
	call inputlag
c
c next line added 1/25/09
c
	if (nop.eq.0.and.nop2.eq.0) return
c
c nbx, nby, nbz are set earlier in sub. mpisetup
c

        nbx=nx+3
        nby=ny+3
        nbz=nz+3
c       ndim=nbx
        ndim=max(nbx,nby,nbz)
       allocate (p(ndim,3),q(0:ndim,3),t(ndim,3),denom(3),ss(0:ndim))
c
	allocate (gx(nx),gy(ny),gz(nz))
c
	allocate (lpfunc(3+nc),lpgrad(3+nc),lplapl(3+nc),ipdis(3+nc))
	allocate (lpsecd(3+nc),lpsdis(3+nc))
c
	lpfunc(1:3)=1
	lpgrad(:)=0
	lplapl(:)=0
	lpsecd(:)=0
	ipdis(:)=0
c
	iovel=1
	ioposn=1
	
	nplu=9+nc
	allocate (npo(nplu))
	nipj=4+2*nc
	allocate (ipj1(nipj),ipj2(nipj))
	npo(:)=0
	ipj1(:)=0
	ipj2(:)=0
c
c initialization calls for spline routines
c
      call solsp1 (p(1,1),q(0,1),t(1,1),nx,ss,denom(1))
      call solsp1 (p(1,2),q(0,2),t(1,2),ny,ss,denom(2))
      call solsp1 (p(1,3),q(0,3),t(1,3),nz,ss,denom(3))
c
	call valid1
c

        jj=0
        do j=-1,1
        do i=-1,1
        if(i.eq.0 .and. j.eq.0) cycle
        jj=jj+1
        id = mod(ipid+i+iproc,iproc)
        jd = mod(jpid+j+jproc,jproc)
        itask = jd*iproc + id
        nbtask(jj) = itask
        enddo
        enddo

        xb = float(nx)/iproc
        zb = float(nz)/jproc

	xiind = nx/iproc + 1
	zjind = nz/jproc + 1
	nbsdim=3

! spline coefficients now a global coarray
#if defined(CAF) && defined(CF_LAG)
	 allocate (bs(xiind,nby,zjind,nbsdim)[0:*])
#else
	 allocate (bs(xiind,nby,zjind,nbsdim))
#endif
c
	nppad = 4

	nump1 = ndpart/numtasks
	if(taskid.eq.0) write(6,*) 'nump1=',nump1
! Note, starting with late-Dec-2008 version, the "pp" and 
! "pp2" arrays are now partitioned across the processors
	if (ndpart.gt.0) then
!	allocate (pp(nump1,ndprop))
! Aug 4 2015, D. Buaria
! pp group now also uses the new local algorithm
	allocate (pp(nppad*nump1,ndprop))
	allocate (pp1ind(nppad*nump1))
	endif
c
	if (taskid.eq.0) write (6,*) ' compart_set, ndprop=',ndprop
	allocate (ps(ndprop))
c
#ifdef LAGSC2
	if (taskid.eq.0) then 
        write (6,601) ndpart, ndprop, ndpart2, ndprop2
 601 	format ('compart_set: ndpart, ndprop,ndpart2,ndprop2=',2(i12,i4))
        end if

	nump2 = ndpart2/numtasks
! 2nd index of pp2 starts from 0, for tagging particles
	if(ndpart2.gt.0) then
! Jul 13 2015, The particle index is now stored in pp2ind 
! Hence 2nd index of pp2 doesnt start from 0 anymore 
!	allocate (pp2(2*nump2,0:ndprop2))
	allocate (pp2(nppad*nump2,ndprop2))
	allocate (pp2ind(nppad*nump2))

! mar 23 2015, D Buaria: The next 2 lines assume that the particle
! population on one task cannot change more than by a factor of 25%
	numpmax = max(nump1, nump2)
	allocate (pin(nppad*numpmax*4,8),pout(nppad*numpmax*4,8))
	allocate (pindin(nppad*numpmax*4,8),pindout(nppad*numpmax*4,8))
!	allocate (pin(100,8),pout(100,8))
	endif
#else
	if (taskid.eq.0) then 
        write (6,602) ndpart, ndprop
 602 	format ('compart_set: ndpart, ndprop=',i8,i3)
        end if
#endif
c

#ifdef MOL
c nompp = no. of molecules per fluid particles PER GROUP
c nom   = nompp x nop x ngm = total no. of molecules
c
#ifdef MOL2
	nom=ndpart2*nompp*ngm
#else
	nom=ndpart*nompp*ngm
#endif
	nm = nom/numtasks
	if (nm.gt.0) then

	allocate(mpp(nm/ngm,9,ngm))
	if (taskid.eq.0) write(6,*) 'compart_set: nop,nop2,nom=',nop,nop2,nom

! Sep 4, 2015, D. Buaria
! index array for molecules, might be useful in future
	allocate (mpind(nm/ngm,ngm))

	endif
#endif


	kipran=11
c
	lupop=101
c
#endif
c

#if defined(CAF) && defined(CF_LAG)
        if (taskid.eq.0) then
        write (6,"('CF_LAG directive: using CAF for particle tracking: Crays only!')") 
        end if
        call setrandlist
#endif

!	call lagcomm_setup
c
	return
	end
