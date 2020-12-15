	subroutine partic_rk4
c
#ifdef LAG
c
	use compart
	implicit none
	include 'intvars'
c
#ifdef FROZEN_SP
        if (jstep.eq.1.and.kstep.eq.1)
     1       allocate (bs(bxisize,nby,bzjsize,nbsdim)[0:*])
#else
        allocate (bs(bxisize,nby,bzjsize,nbsdim)[0:*])
#endif
c
c next line added 1/25/09
c
c
	if (nop.eq.0.and.nop2.eq.0) return 
c
	if (lpfunc(1).eq.0) return
c
         allocate (bfpos(3,ndpart))
#ifdef LAGSC2
        allocate (bfpos2(3,ndpart2))
#endif

 	if (kstep.eq.1) call partsp_1 
 	if (kstep.eq.2) call partsp_2 
 	if (kstep.eq.3) call partsp_3 
 	if (kstep.eq.4) call partsp_4 
c
#ifndef FROZEN_SP
 	deallocate (bs)
#endif
c
 	if (iout.eq.1) call lagout1
 	if (iout.eq.2) call lagout2
c
	deallocate (bfpos)
#ifdef LAGSC2
	deallocate (bfpos2)
#endif
c
#endif
c
	return
	end
