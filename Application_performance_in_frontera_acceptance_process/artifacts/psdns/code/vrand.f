      subroutine vrand(i,x,s,n)
c
c  routine to generate a vector of pseudo-random
c  numbers uniformly distributed on the exclusive
c  interval (0,1).
c  for use on 32 bit machines.
c
c  arguments:
c    i    integer seed ( i .ne. 0 ).
c    x    vector in which random numbers are returned.
c    n    dimension of vector x.
c    s    dummy argument - set to 1 in calling program for
c          consistency with ap164 vrand.
c
c  method:
c    for seed i (0 < i < 2**31-1), ix = mod(65539*i,2**31) is
c    a pseudo-random integer uniform on the exclusive interval
c    (0,2**31-1).  x = ix*0.46566127e-9 is then a pseudo-random
c    real number uniform on the exclusive interval (0.,1.).
c    i is replaced by ix on each call.
c    takes advantage of 32 bit machine's truncation of
c    integers >= 2**31.
c
c  notes:
c    in single precision (32 bit) 4.6566127e-10 <= x <= (1.0 - 5.96e-8).
c    if i = 0, seed locks at 0 and x(j) = 0. for all j.
c
c  other routines called - none.
c
      dimension x(n)
      integer s
c
c  check for i = 0 - invalid seed
c
      if( i .eq. 0 ) go to 300
c
c  generate random numbers
c
      do 100 j=1,n
       i = i*65539
       if( i .lt. 0 ) i=i+2147483647+1
       xi = i
       x(j) = xi*0.46566127e-9
100    continue
c
      return
c
300   write(6,600)
      stop
c
600   format(/,'error in vrand - seed = 0.')
c
      end
