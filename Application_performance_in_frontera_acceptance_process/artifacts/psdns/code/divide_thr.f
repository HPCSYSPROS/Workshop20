      subroutine divide_thr (ntotal,p1,ipi,label)
c
c This routine partitions 'ntotal' units of workload among 
c OpenMP threads. The distribution is made as close to even
c as possible: if ntotal is not an integer multiple of num_thr
c then some threads take on one more unit than the others.
c
c p1 is an array that will hold the ID of the first work unit
c assigned to each thread
c pi holds the number of work units assigned
c
c if ntotal < num_thr the code aborts. This may suggest a need
c for change of programming strategy in the calling routine
c or it may indicate a case of overkill in terms of numtasks
c and/or num_thr



      use com
      include 'intvars'
	character*(*) label
c     
	integer nchar,nn,ntotal
c
	integer ithr
c
	integer p1(0:num_thr-1),ipi(0:num_thr-1)
	integer l,psum
	
	nchar=len(label)
	nn=min0(nchar,20)

	if (ntotal.lt.num_thr) then
	write (6,"('divide_thr: ntotal,num_thr=',2i4,4x,a20)") ntotal,num_thr,
     1            label(1:nn)
	write (6,*) 'divide_thr: stop: invalid choice of parameters'
	stop
	end if
c
	ipi(:)=ntotal/num_thr
	l=mod(ntotal,num_thr) 
	if (l.gt.0) then
	do ithr=0,l-1
	ipi(ithr)=ipi(ithr)+1
	end do
	end if
	p1(0)=1
	psum=p1(0)
	do ithr=1,num_thr-1
	psum=psum+ipi(ithr-1)
	p1(ithr)=psum
	end do
c
      return
      end
