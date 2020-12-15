      real function ranu (ii)
c
c  routine to generate independent random numbers
c  uniformly distributed in the exclusive interval (0,1).
c
c  the argument ii is a dummy
c
c
	use com
c
#ifdef RANDOM_SEED
	use ranseed
#else
      common/rancom/ix,rseed(2)
      real rseed
#endif
c
c
      save icount
      data icount/0/
c
      if (iseed.gt.0) then
c
c  generate random number
c
c     icount=icount+1
          iseed=iseed*65539
          if (iseed.lt.0) iseed=iseed+2147483647+1
          ranu=iseed*0.46566127e-9
          return
c
        else
           write(6,*)' seed iseed=',iseed,' must be positive',taskid
           write(6,*)' stopped in ranu '
      call MPI_ABORT (MPI_COMM_WORLD,ierror)
c
        endif
c
      end
