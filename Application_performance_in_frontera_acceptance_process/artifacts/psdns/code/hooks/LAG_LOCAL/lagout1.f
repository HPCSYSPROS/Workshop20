       subroutine lagout1
c
c routine to control output of particle properties,
c called from sub. partic
c
c
#ifdef LAG
c
	use compart
	implicit none
	integer ic
c
c write particle properties on file, including velocities if desired
c vel. gradients also (controlled by LGRAD)


#ifdef NOV15
#else


	if (nop.gt.0) then
!	call pop_file4 (1)
! 	call pop_write4 (1,4,6)
!	call pop_write4 (nplu-1,1,3)
	call pop1_file_local (1)
	call pop1_write_local2 (1,4,6)
	call pop1_write_local2 (nplu-1,1,3)
#ifdef LGRAD
	if (lpgrad(1).eq.1) then
!  	call pop_write4 (2,ipvg(1),ipvg(8))
  	call pop1_write_local2 (2,ipvg(1),ipvg(8))
	end if
#endif
	endif
c

#ifdef LAGSC2
	if (nop2.gt.0) then

c write positions and velocities and vel. gradients
	call pop2_file5 (1)
!  	call pop2_write5 (1,4,6)
! 	call pop2_write5 (nplu-1,1,3)
  	call pop2_write_local (1,4,6)
 	call pop2_write_local (nplu-1,1,3)
#ifdef LGRAD
	if (lpgrad(1).eq.1) then
!  	call pop2_write5 (2,ipvg(1),ipvg(8))
  	call pop2_write_local (2,ipvg(1),ipvg(8))
	end if
#endif
	endif
#endif



#endif
!NOV15 directive

#ifdef MOL
	if (nom.gt.0) then
	call mol_file3 (1)
	call mol_write3 (1,4,6)
	call mol_write3 (nplu-1,1,3)
	end if
#endif
c
#ifdef NOV15
#else
#ifndef NOSCALAR
	do 20 ic=1,nc
	if (lpfunc(ic+3).eq.0.and.lpgrad(ic+3).eq.0) go to 20
	call pop2_write5 (3+ic,ipj1(ic),ipj2(ic))
 20	continue
#endif
#endif

c
#endif

 90	continue

	if (taskid.eq.0) write (6,*) 'exit lagout1, istep=,istep'
      return
      end
