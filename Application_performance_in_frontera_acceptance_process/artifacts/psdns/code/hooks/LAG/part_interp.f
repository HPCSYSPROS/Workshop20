       subroutine part_interp (indexp)
c
c interpolate for other flow variables at particle positions
c
#ifdef LAG
c
	use compart
	implicit none
	include 'intvars'
c
	integer indexp
c
#ifdef NOV15
#else
c interpolate for scalars, scalar gradients and velocity gradients  if needed
c

#ifndef LAGSC2
        call partsc (indexp)
#else
        call part2sc (indexp)
#endif


#ifdef LGRAD
        if (any(lpgrad(1:3).gt.0)) then

        call partvg (indexp)

#ifdef LAGSC2
! 5/8/12, Dhawal Buaria : part2vg now merged with partvg
!        call part2vg (indexp)
#endif

        end if
#endif


#endif 
#endif 
!ifdef LAG
c
      return
      end
