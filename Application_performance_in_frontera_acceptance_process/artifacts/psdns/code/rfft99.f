        subroutine rfft99(u,work,itrans,ivec,ntrans,nvec,isgn)
c
c       routine to perform a real fft on a plane of data.
c               transform is normalized
c
c       the arguments are essentially the same as in rogallo's code
c               except that his third and fourth arguments are omitted.
c
c       u (forward) real array containing data to be transformed.
c                       u(1) contains the first element of the first vec
c       u (inverse) complex array containing the first ntrans/2
c                       fourier coefficients.
c                       the transforms are returned in u
c       itrans  transform stride ( in real words ) (must be 1 ? )
c       ivec            vector stride
c       ntrans  length of transform (must be even)
c       nvec            number of vectors
c       isgn            = -1 for forward transform, =+1 for inverse
c       work            complex work array of dimension 2*ntrans
c
c       method: each vector is copied into a sequential one-dimensional
c                       array. subroutine fourt is then called to perfor
c                       the vector is the copied back into u.
c
c
        complex work(ntrans,2)
c       real u(2)
        real u(1)
	real norm
	integer ndim(1)
c
	ndim(1)=ntrans
c
c  forward transform
c
        if( isgn .eq. 1 ) go to 199
c
c  loop over vectors
c
        norm=1./float( ntrans )
        ik=1-ivec
        do 100 i=1,nvec
        ik=ik+ivec
c
c  copy u into (complex) vector and set imaginary part to zero
c
        k=ik-itrans
        do 40 j=1,ntrans
        k=k+itrans
40      work(j,1)=cmplx( u(k) , 0. )
c
c  transform
c
        call fourt(work,ndim(1),1,-1,0,work(1,2) )
c
c  copy first half of transform back into u. (second half is
c               the complex conjugate of the first)
c
        k=ik-itrans
        do 60 j=1,ntrans/2
        k=k+itrans
        u(k)=real( work(j,1) ) * norm
        k=k+itrans
60      u(k)=aimag( work(j,1) ) * norm
c
100     continue
c
        return
c
c
c  inverse transform
c
199     continue
c     write (6,*) ' rfft99: before do 200'
c
c  loop over vectors
c
        ik=1-ivec
        do 200 i=1,nvec
        ik=ik+ivec
c
c  copy complex elements form u into first half of work
c
        k=ik-itrans*2
        do 140 j=1,ntrans/2
        k=k+itrans*2
140     work(j,1)=cmplx( u(k) , u(k+1) )
c
c  use conjugate symmetry to form remainder of vector
c
        work(ntrans/2+1,1)=cmplx(0.,0.)
c
        do 150 j=ntrans/2+2,ntrans
        jj=ntrans+2-j
150     work(j,1)=conjg( work(jj,1) )
c
c***
c       work(1,1)=cmplx(real(work(1,1)),0.)
c***
c
c  transform
c
        call fourt(work,ndim(1),1,1,1,work(1,2) )
c
c  copy real part of transform into u
c
        k=ik-itrans
        do 160 j=1,ntrans
        k=k+itrans
160     u(k)=real( work(j,1) )
c
200     continue
c
        return
        end
