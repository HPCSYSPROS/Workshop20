	subroutine compart_set 
c	
#ifdef LAG
c
	use compart
	implicit none
	integer nipj
c
c 9/13/2016: these 8 lines have to appear before calling inputlag

	allocate (lpfunc(3+nc),lpgrad(3+nc),lplapl(3+nc),ipdis(3+nc))
	allocate (lpsecd(3+nc),lpsdis(3+nc))
c
	lpfunc(1:3)=1
	lpgrad(:)=0
	lplapl(:)=0
	lpsecd(:)=0
	ipdis(:)=0
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
	nbsdim=3

	 allocate (bs(bxisize,nby,bzjsize,nbsdim))
c
c Note, starting with late-Dec-2008 version, the "pp" and 
c "pp2" arrays are now partitioned across the processors
c
	allocate (pp(ndpart/numtasks,ndprop))
c
	if (taskid.eq.0) write (6,*) ' compart_set, ndprop=',ndprop
	allocate (ps(ndprop))
c
#ifdef LAGSC2
	if (taskid.eq.0) then 
        write (6,601) ndpart, ndprop, ndpart2, ndprop2
 601 	format ('compart_set: ndpart, ndprop,ndpart2,ndprop2=',2(i8,i3))
        end if
	allocate (pp2(ndpart2/numtasks,ndprop2))
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
	if (nm.gt.0) allocate(mpp(nm/ngm,9,ngm))
	if (taskid.eq.0) write(6,*) 'compart_set: nop,nop2,nom=',nop,nop2,nom
#endif

ccc

	kipran=11
c
	lupop=101
c
#endif
c

#ifdef CF_LAG
        if (taskid.eq.0) then
        write (6,"('CF_LAG directive: using CAF for particle tracking: Crays only!')") 
        end if
        call setrandlist()
#endif


c
	return
	end
