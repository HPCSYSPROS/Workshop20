      subroutine forcdt
c
c
	use comsp
	implicit none

	integer i,k
c
c
#ifdef EPFOR

c     if (kforce.gt.0..and.mod(istep-1,iostep).eq.0) then
c     if (kforce.gt.0..and.(istep.eq.1.or.iopflag.eq.1)) then
      if (kforce.gt.0..and.(jstep.eq.1.or.iopflag.eq.1)) then
c
c add results for energy input from all slabs
c
      call MPI_REDUCE (efkz,efki,3*nxh,mpireal,MPI_SUM,0,
     1                 MPI_COMM_WORLD,ierr)
c
      if (taskid.eq.0) then
c
c divide by dt since in DNS code "for" is velocity increment over dt
c
      do i=1,3
      do k=1,kfor
      efki(k,i)=efki(k,i)/dt
      end do
      end do
c
      erate=0.
      eirate(1)=0.
      eirate(2)=0.
      eirate(3)=0.
      do k=1,kfor
      efk(k)=.5*(efki(k,1)+efki(k,2)+efki(k,3))
      eirate(1)=eirate(1)+efki(k,1)
      eirate(2)=eirate(2)+efki(k,2)
      eirate(3)=eirate(3)+efki(k,3)
      erate=erate+efk(k)
      end do
c
c forcing energy input information (added 6/27/96)
c
      if (istep.eq.1) then
      write (50,243) kfor,istep-1,time-dt
      else
      call clapf (50,'forcinp')
      write (50,210) istep-1,time-dt
      end if
      write (50,244) erate,epslon,erate-epslon
      write (50,245) eirate
      write (50,246) (efk(k),k=1,kfor)
c
      end if
c
      return
c
      end if
c
#endif
        return
c
 210  format (/'istep=',i6,2x,'time=',1p,e13.5)
 243  format ('forcing input: kfor=',i3/'istep=',i6,2x,'time=',1p,e13.5)
 244  format ('force/diss/balance=',1p,3e12.4)
 245  format ('input in each comp=',1p,3e12.4)
 246  format ('spectrum:',2x,1p,5e13.5)
c
        end
