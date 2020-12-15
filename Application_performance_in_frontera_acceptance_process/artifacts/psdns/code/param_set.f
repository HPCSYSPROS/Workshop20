      subroutine param_set
c
	use param
c
	use com, only: nrproc,kforce,beta1,beta2,beta3,taskid
c
c
c note that nx, ny, nz, nc, kfor are to be read in
c (set up a subroutine input.f)
c
      nxhvs=nx/2 ; nyvs=ny ; nzvs=nz
c
c  use this statement for shear flow              
        nxhvs = nx/2 ; nyvs =ny  ; nzvs = nz
c                                                 
c  use this statement if no shear                 
c       parameter ( nxhvs = 1 , nyvs =  1  , nzvs = 1 )                 
c                                                 
      if (nc.eq.0) then
      ncd=1
      ncop=0
      nut=2
      else
      ncd=nc
      ncop=nc
      nut=ncd+ncop
      end if
c
c kfor must be at least 1, even if no forcing used
c                                                 
c need larger kfor for runs on bigger domains
c
      kfor=kforce/min(beta1,beta2+beta3)+1.5
c	
	if (taskid.eq.0) write (6,*) 'param_set: kfor=',kfor
c
	k2fo=2*kfor-1
c                                                 
c  derived parameters                             
c                                                 
	nxh=nx/2 ; nxhp=nxh+1 ; nxhp2=2*nxhp
	nyh=ny/2 ; nyhp=nyh+1 ; nyhp2=2*nyhp
	nzh=nz/2 ; nzhp=nzh+1 ; nzhp2=2*nzhp

!! NPM edit: set dealiasing in input file

        if(nxpad.lt.nx) then !user disable dealiasing for x
           nxpad = nx
        endif
        
        if(nypad.lt.ny) then !user disable dealiasing for y
           nypad = ny
        endif
        
        if(nzpad.lt.nz) then !user disable dealiasing for z
           nzpad = nz     
        endif        

#ifdef HOMOGENEOUS
	nxpad=nx
	nypad=ny
	nzpad=nz
#endif

	nxhpad=nxpad/2 ; nxhppad=nxhpad+1 ; nxhp2pad=2*nxhppad
	nyhpad=nypad/2 ; nyhppad=nyhpad+1 
	nzhpad=nzpad/2 ; nzhppad=nzhpad+1 

       is1=4+nc ; is2=5+nc ; ipsi=4+nc 
                                                 
#ifdef ROTD
      nu=6+nc
#else
      if (nc.eq.1) then
      nu=is2
      else
      nu=3+2*ncd
      end if
#endif

#ifdef CVC_PRESS
        nu=6+nc
#endif
        ncpp=nc*nc+7*nc+12 ; ncp=ncpp/2                   
                                                 
        k2fo=2*kfor-1
        nrproc=4*kfor*k2fo*k2fo
                                                 
      return
      end
