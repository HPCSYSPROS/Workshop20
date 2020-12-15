c fourt      from xlib                                     09/08/82
      subroutine fourt(data,nn,ndim,isign,iform,work)
c
c     the cooley-tukey fast fourier transform in usasi basic fortran
c
c     transform(j1,j2,,,,) = sum(data(i1,i2,,,,)*w1**((i2-1)*(j2-1))
c                                 *w2**((i2-1)*(j2-1))*,,,),
c     where i1 and j1 run from 1 to nn(1) and w1=exp(isign*2*pi=
c     sqrt(-1)/nn(1)), etc.  there is no limit on the dimensionality
c     (number of subscripts) of the data array.  if an inverse
c     transform (isign=+1) is performed upon an array of transformed
c     (isign=-1) data, the original data will reappear.
c     multiplied by nn(1)*nn(2)*,,,  the array of input data must be
c     in complex format.  however, if all imaginary parts are zero (i.e.
c     the data are disguised real) running time is cut up to forty per-
c     cent.  (for fastest transform of real data, nn(1) should be even.)
c     the transform values are always complex and are returned in the
c     original array of data, replacing the input data.  the length
c     of each dimension of the data array may be any integer.  the
c     program runs faster on composite integers than on primes, and is
c     particularly fast on numbers rich in factors of two.
c
c     timing is in fact given by the following formula.  let ntot be the
c     total number of points (real or complex) in the data array, that
c     is, ntot=nn(1)*nn(2)*...  decompose ntot into its prime factors,
c     such as 2**k2 * 3**k3 * 5**k5 * ...  let sum2 be the sum of all
c     the factors of two in ntot, that is, sum2 = 2*k2.  let sumf be
c     the sum of all other factors of ntot, that is, sumf = 3*k3*5*k5*..
c     the time taken by a multidimensional transform on these ntot data
c     is t = t0 + ntot*(t1+t2*sum2+t3*sumf).  on the cdc 3300 (floating
c     point add time = six microseconds), t = 3000 + ntot*(600+40*sum2+
c     175*sumf) microseconds on complex data.
c
c     implementation of the definition by summation will run in a time
c     proportional to ntot*(nn(1)+nn(2)+...).  for highly composite ntot
c     the savings offered by this program can be dramatic.  a one-dimen-
c     sional array 4000 in length will be transformed in 4000*(600+
c     40*(2+2+2+2+2)+175*(5+5+5)) = 14.5 seconds versus about 4000*
c     4000*175 = 2800 seconds for the straightforward technique.
c
c     the fast fourier transform places three restrictions upon the
c     data.
c     1.  the number of input data and the number of transform values
c     must be the same.
c     2.  both the input data and the transform values must represent
c     equispaced points in their respective domains of time and
c     frequency.  calling these spacings deltat and deltaf, it must be
c     true that deltaf=2*pi/(nn(i)*deltat).  of course, deltat need not
c     be the same for every dimension.
c     3.  conceptually at least, the input data and the transform output
c     represent single cycles of periodic functions.
c
c     the calling sequence is--
c     call fourt(data,nn,ndim,isign,iform,work)
c
c     data is the array used to hold the real and imaginary parts
c     of the data on input and the transform values on output.  it
c     is a multidimensional floating point array, with the real and
c     imaginary parts of a datum stored immediately adjacent in storage
c     (such as fortran iv places them).  normal fortran ordering is
c     expected, the first subscript changing fastest.  the dimensions
c     are given in the integer array nn, of length ndim.  isign is -1
c     to indicate a forward transform (exponential sign is -) and +1
c     for an inverse transform (sign is +).  iform is +1 if the data are
c     complex, 0 if the data are real.  if it is 0, the imaginary
c     parts of the data must be set to zero.  as explained above, the
c     transform values are always complex and are stored in array data.
c     work is an array used for working storage.  it is floating point
c     real, one dimensional of length equal to twice the largest array
c     dimension nn(i) that is not a power of two.  if all nn(i) are
c     powers of two, it is not needed and may be replaced by zero in the
c     calling sequence.  thus, for a one-dimensional array, nn(1) odd,
c     work occupies as many storage locations as data.  if supplied,
c     work must not be the same array as data.  all subscripts of all
c     arrays begin at one.
c
c     example 1.  three-dimensional forward fourier transform of a
c     complex array dimensioned 32 by 25 by 13 in fortran iv.
c     dimension data(32,25,13),work(50),nn(3)
c     complex data
c     data nn/32,25,13/
c     do 1 i=1,32
c     do 1 j=1,25
c     do 1 k=1,13
c  1  data(i,j,k)=complex value
c     call fourt(data,nn,3,-1,1,work)
c
c     example 2.  one-dimensional forward transform of a real array of
c     length 64 in fortran ii,
c     dimension data(2,64)
c     do 2 i=1,64
c     data(1,i)=real part
c  2  data(2,i)=0.
c     call fourt(data,64,1,-1,0,0)
c
c     there are no error messages or error halts in this program.  the
c     program returns immediately if ndim or any nn(i) is less than one.
c
c     program by norman brenner from the basic program by charles
c     rader,  june 1967.  the idea for the digit reversal was
c     suggested by ralph alter.
c
c     this is the fastest and most versatile version of the fft known
c     to the author.  a program called four2 is available that also
c     performs the fast fourier transform and is written in usasi basic
c     fortran.  it is about one third as long and restricts the
c     dimensions of the input array (which must be complex) to be powers
c     of two.  another program, called four1, is one tenth as long and
c     runs two thirds as fast on a one-dimensional complex array whose
c     length is a power of two.
c
c     reference--
c     ieee audio transactions (june 1967), special issue on the fft.
      dimension data(*),nn(1),ifact(32),work(1)
      data twopi/6.2831853071796/,rthlf/0.70710678118655/
      data nprev/0/,np0/0/
c the following call is for gathering statistics on library use at ncar
c     call q8qst4( 4hxlib      , 5hfourt     ,5hfourt  ,10hversion  9)
      if(ndim-1)920,1,1
1     ntot=2
      do 2 idim=1,ndim
      if(nn(idim))920,920,2
2     ntot=ntot*nn(idim)
c
c     main loop for each dimension
c
      np1=2
      do 910 idim=1,ndim
      n=nn(idim)
      np2=np1*n
      if(n-1)920,900,5
c
c     is n a power of two and if not, what are its factors
c
5     m=n
      ntwo=np1
      if=1
      idiv=2
10    iquot=m/idiv
      irem=m-idiv*iquot
      if(iquot-idiv)50,11,11
11    if(irem)20,12,20
12    ntwo=ntwo+ntwo
      ifact(if)=idiv
      if=if+1
      m=iquot
      go to 10
20    idiv=3
      inon2=if
30    iquot=m/idiv
      irem=m-idiv*iquot
      if(iquot-idiv)60,31,31
31    if(irem)40,32,40
32    ifact(if)=idiv
      if=if+1
      m=iquot
      go to 30
40    idiv=idiv+2
      go to 30
50    inon2=if
      if(irem)60,51,60
51    ntwo=ntwo+ntwo
      go to 70
60    ifact(if)=m
c
c     separate four cases--
c        1. complex transform or real transform for the 4th, 9th,etc.
c           dimensions.
c        2. real transform for the 2nd or 3rd dimension.  method--
c           transform half the data, supplying the other half by con-
c           jugate symmetry.
c        3. real transform for the 1st dimension, n odd.  method--
c           set the imaginary parts to zero.
c        4. real transform for the 1st dimension, n even.  method--
c           transform a complex array of length n/2 whose real parts
c           are the even numbered real values and whose imaginary parts
c           are the odd numbered real values.  separate and supply
c           the second half by conjugate symmetry.
c
70    icase=1
      ifmin=1
      i1rng=np1
      if(idim-4)71,100,100
71    if(iform)72,72,100
72    icase=2
      i1rng=np0*(1+nprev/2)
      if(idim-1)73,73,100
73    icase=3
      i1rng=np1
      if(ntwo-np1)100,100,74
74    icase=4
      ifmin=2
      ntwo=ntwo/2
      n=n/2
      np2=np2/2
      ntot=ntot/2
      i=1
      do 80 j=1,ntot
      data(j)=data(i)
80    i=i+2
c
c     shuffle data by bit reversal, since n=2**k.  as the shuffling
c     can be done by simple interchange, no working array is needed
c
100   if(ntwo-np2)200,110,110
110   np2hf=np2/2
      j=1
      do 150 i2=1,np2,np1
      if(j-i2)120,130,130
120   i1max=i2+np1-2
      do 125 i1=i2,i1max,2
      do 125 i3=i1,ntot,np2
      j3=j+i3-i2
      tempr=data(i3)
      tempi=data(i3+1)
      data(i3)=data(j3)
      data(i3+1)=data(j3+1)
      data(j3)=tempr
125   data(j3+1)=tempi
130   m=np2hf
140   if(j-m)150,150,145
145   j=j-m
      m=m/2
      if(m-np1)150,140,140
150   j=j+m
      go to 300
c
c     shuffle data by digit reversal for general n
c
200   nwork=2*n
      do 270 i1=1,np1,2
      do 270 i3=i1,ntot,np2
      j=i3
      do 260 i=1,nwork,2
      if(icase-3)210,220,210
210   work(i)=data(j)
      work(i+1)=data(j+1)
      go to 230
220   work(i)=data(j)
      work(i+1)=0.
230   ifp2=np2
      if=ifmin
240   ifp1=ifp2/ifact(if)
      j=j+ifp1
      if(j-i3-ifp2)260,250,250
250   j=j-ifp2
      ifp2=ifp1
      if=if+1
      if(ifp2-np1)260,260,240
260   continue
      i2max=i3+np2-np1
      i=1
      do 270 i2=i3,i2max,np1
      data(i2)=work(i)
      data(i2+1)=work(i+1)
270   i=i+2
c
c     main loop for factors of two.  perform fourier transforms of
c     length four, with one of length two if needed.  the twiddle factor
c     w=exp(isign*2*pi*sqrt(-1)*m/(4*mmax)).  check for w=isign*sqrt(-1)
c     and repeat for w=w*(1+isign*sqrt(-1))/sqrt(2).
c
300   if(ntwo-np1)600,600,305
305   np1tw=np1+np1
      ipar=ntwo/np1
310   if(ipar-2)350,330,320
320   ipar=ipar/4
      go to 310
330   do 340 i1=1,i1rng,2
      do 340 k1=i1,ntot,np1tw
      k2=k1+np1
      tempr=data(k2)
      tempi=data(k2+1)
      data(k2)=data(k1)-tempr
      data(k2+1)=data(k1+1)-tempi
      data(k1)=data(k1)+tempr
340   data(k1+1)=data(k1+1)+tempi
350   mmax=np1
360   if(mmax-ntwo/2)370,600,600
370   lmax=max0(np1tw,mmax/2)
      do 570 l=np1,lmax,np1tw
      m=l
      if(mmax-np1)420,420,380
380   theta=-twopi*float(l)/float(4*mmax)
      if(isign)400,390,390
390   theta=-theta
400   wr=cos(theta)
      wi=sin(theta)
410   w2r=wr*wr-wi*wi
      w2i=2.*wr*wi
      w3r=w2r*wr-w2i*wi
      w3i=w2r*wi+w2i*wr
420   do 530 i1=1,i1rng,2
      kmin=i1+ipar*m
      if(mmax-np1)430,430,440
430   kmin=i1
440   kdif=ipar*mmax
450   kstep=4*kdif
      if(kstep-ntwo)460,460,530
460   do 520 k1=kmin,ntot,kstep
      k2=k1+kdif
      k3=k2+kdif
      k4=k3+kdif
      if(mmax-np1)470,470,480
470   u1r=data(k1)+data(k2)
      u1i=data(k1+1)+data(k2+1)
      u2r=data(k3)+data(k4)
      u2i=data(k3+1)+data(k4+1)
      u3r=data(k1)-data(k2)
      u3i=data(k1+1)-data(k2+1)
      if(isign)471,472,472
471   u4r=data(k3+1)-data(k4+1)
      u4i=data(k4)-data(k3)
      go to 510
472   u4r=data(k4+1)-data(k3+1)
      u4i=data(k3)-data(k4)
      go to 510
480   t2r=w2r*data(k2)-w2i*data(k2+1)
      t2i=w2r*data(k2+1)+w2i*data(k2)
      t3r=wr*data(k3)-wi*data(k3+1)
      t3i=wr*data(k3+1)+wi*data(k3)
      t4r=w3r*data(k4)-w3i*data(k4+1)
      t4i=w3r*data(k4+1)+w3i*data(k4)
      u1r=data(k1)+t2r
      u1i=data(k1+1)+t2i
      u2r=t3r+t4r
      u2i=t3i+t4i
      u3r=data(k1)-t2r
      u3i=data(k1+1)-t2i
      if(isign)490,500,500
490   u4r=t3i-t4i
      u4i=t4r-t3r
      go to 510
500   u4r=t4i-t3i
      u4i=t3r-t4r
510   data(k1)=u1r+u2r
      data(k1+1)=u1i+u2i
      data(k2)=u3r+u4r
      data(k2+1)=u3i+u4i
      data(k3)=u1r-u2r
      data(k3+1)=u1i-u2i
      data(k4)=u3r-u4r
520   data(k4+1)=u3i-u4i
      kdif=kstep
      kmin=4*(kmin-i1)+i1
      go to 450
530   continue
      m=m+lmax
      if(m-mmax)540,540,570
540   if(isign)550,560,560
550   tempr=wr
      wr=(wr+wi)*rthlf
      wi=(wi-tempr)*rthlf
      go to 410
560   tempr=wr
      wr=(wr-wi)*rthlf
      wi=(tempr+wi)*rthlf
      go to 410
570   continue
      ipar=3-ipar
      mmax=mmax+mmax
      go to 360
c
c     main loop for factors not equal to two.  apply the twiddle factor
c     w=exp(isign*2*pi*sqrt(-1)*(j1-1)*(j2-j1)/(ifp1+ifp2)), then
c     perform a fourier transform of length ifact(if), making use of
c     conjugate symmetries.
c
600   if(ntwo-np2)605,700,700
605   ifp1=ntwo
      if=inon2
      np1hf=np1/2
610   ifp2=ifact(if)*ifp1
      j1min=np1+1
      if(j1min-ifp1)615,615,640
615   do 635 j1=j1min,ifp1,np1
      theta=-twopi*float(j1-1)/float(ifp2)
      if(isign)625,620,620
620   theta=-theta
625   wstpr=cos(theta)
      wstpi=sin(theta)
      wr=wstpr
      wi=wstpi
      j2min=j1+ifp1
      j2max=j1+ifp2-ifp1
      do 635 j2=j2min,j2max,ifp1
      i1max=j2+i1rng-2
      do 630 i1=j2,i1max,2
      do 630 j3=i1,ntot,ifp2
      tempr=data(j3)
      data(j3)=data(j3)*wr-data(j3+1)*wi
630   data(j3+1)=tempr*wi+data(j3+1)*wr
      tempr=wr
      wr=wr*wstpr-wi*wstpi
635   wi=tempr*wstpi+wi*wstpr
640   theta=-twopi/float(ifact(if))
      if(isign)650,645,645
645   theta=-theta
650   wstpr=cos(theta)
      wstpi=sin(theta)
      j2rng=ifp1*(1+ifact(if)/2)
      do 695 i1=1,i1rng,2
      do 695 i3=i1,ntot,np2
      j2max=i3+j2rng-ifp1
      do 690 j2=i3,j2max,ifp1
      j1max=j2+ifp1-np1
      do 680 j1=j2,j1max,np1
      j3max=j1+np2-ifp2
      do 680 j3=j1,j3max,ifp2
      jmin=j3-j2+i3
      jmax=jmin+ifp2-ifp1
      i=1+(j3-i3)/np1hf
      if(j2-i3)655,655,665
655   sumr=0.
      sumi=0.
      do 660 j=jmin,jmax,ifp1
 659  sumr=sumr+data(j)
660   sumi=sumi+data(j+1)
      work(i)=sumr
      work(i+1)=sumi
      go to 680
665   iconj=1+(ifp2-2*j2+i3+j3)/np1hf
      j=jmax
      sumr=data(j)
      sumi=data(j+1)
      oldsr=0.
      oldsi=0.
      j=j-ifp1
670   tempr=sumr
      tempi=sumi
      sumr=twowr*sumr-oldsr+data(j)
      sumi=twowr*sumi-oldsi+data(j+1)
      oldsr=tempr
      oldsi=tempi
      j=j-ifp1
      if(j-jmin)675,675,670
675   tempr=wr*sumr-oldsr+data(j)
      tempi=wi*sumi
      work(i)=tempr-tempi
      work(iconj)=tempr+tempi
      tempr=wr*sumi-oldsi+data(j+1)
      tempi=wi*sumr
      work(i+1)=tempr+tempi
      work(iconj+1)=tempr-tempi
680   continue
      if(j2-i3)685,685,686
685   wr=wstpr
      wi=wstpi
      go to 690
686   tempr=wr
      wr=wr*wstpr-wi*wstpi
      wi=tempr*wstpi+wi*wstpr
690   twowr=wr+wr
      i=1
      i2max=i3+np2-np1
      do 695 i2=i3,i2max,np1
      data(i2)=work(i)
      data(i2+1)=work(i+1)
695   i=i+2
      if=if+1
      ifp1=ifp2
      if(ifp1-np2)610,700,700
c
c     complete a real transform in the 1st dimension, n even, by con-
c     jugate symmetries.
c
700   go to (900,800,900,701),icase
701   nhalf=n
      n=n+n
      theta=-twopi/float(n)
      if(isign)703,702,702
702   theta=-theta
703   wstpr=cos(theta)
      wstpi=sin(theta)
      wr=wstpr
      wi=wstpi
      imin=3
      jmin=2*nhalf-1
      go to 725
710   j=jmin
      do 720 i=imin,ntot,np2
      sumr=(data(i)+data(j))/2.
      sumi=(data(i+1)+data(j+1))/2.
      difr=(data(i)-data(j))/2.
      difi=(data(i+1)-data(j+1))/2.
      tempr=wr*sumi+wi*difr
      tempi=wi*sumi-wr*difr
      data(i)=sumr+tempr
      data(i+1)=difi+tempi
      data(j)=sumr-tempr
      data(j+1)=-difi+tempi
720   j=j+np2
      imin=imin+2
      jmin=jmin-2
      tempr=wr
      wr=wr*wstpr-wi*wstpi
      wi=tempr*wstpi+wi*wstpr
725   if(imin-jmin)710,730,740
730   if(isign)731,740,740
731   do 735 i=imin,ntot,np2
735   data(i+1)=-data(i+1)
740   np2=np2+np2
      ntot=ntot+ntot
      j=ntot+1
      imax=ntot/2+1
745   imin=imax-2*nhalf
      i=imin
      go to 755
750   data(j)=data(i)
      data(j+1)=-data(i+1)
755   i=i+2
      j=j-2
      if(i-imax)750,760,760
760   data(j)=data(imin)-data(imin+1)
      data(j+1)=0.
      if(i-j)770,780,780
765   data(j)=data(i)
      data(j+1)=data(i+1)
770   i=i-2
      j=j-2
      if(i-imin)775,775,765
775   data(j)=data(imin)+data(imin+1)
      data(j+1)=0.
      imax=imin
      go to 745
780   data(1)=data(1)+data(2)
      data(2)=0.
      go to 900
c
c     complete a real transform for the 2nd or 3rd dimension by
c     conjugate symmetries.
c
800   if(i1rng-np1)805,900,900
805   do 860 i3=1,ntot,np2
      i2max=i3+np2-np1
      do 860 i2=i3,i2max,np1
      imin=i2+i1rng
      imax=i2+np1-2
      jmax=2*i3+np1-imin
      if(i2-i3)820,820,810
810   jmax=jmax+np2
820   if(idim-2)850,850,830
830   j=jmax+np0
      do 840 i=imin,imax,2
      data(i)=data(j)
      data(i+1)=-data(j+1)
840   j=j-2
850   j=jmax
      do 860 i=imin,imax,np0
      data(i)=data(j)
      data(i+1)=-data(j+1)
860   j=j-np0
c
c     end of loop on each dimension
c
900   np0=np1
      np1=np2
910   nprev=n
920   return
      end
