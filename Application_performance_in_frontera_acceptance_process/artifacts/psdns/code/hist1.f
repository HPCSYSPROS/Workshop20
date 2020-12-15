        subroutine hist1( ifc,ilc,kwt, u,wt,n,stride,
     1          nddu, bldd,brdd, hist, fmt )
!
!  routine to form a histogram from data samples.
!  each sample can be ascribed a weight wt.  the histogram of
!  u is formed.  the routine can be called
!  more than once to average over several sets of samples.
!
! --- modified for formatted output
!
!  input:
!       ifc = 1 on first call
!           = 0 on subsequent calls
!       ilc = 0 on all but last call
!           = -1 on last call ( normalization performed )
!           = lu > 0 on last call, normalization performed,
!               and data written on logical unit lu (in "f" format)
!       kwt = 0 if equal weights are to be used
!           = 1 if weights are to be taken from wt
!       u - vector containing samples of u
!       wt- vector containing weights ( if kwt=0, wt is not referenced )
!       n - number of samples
!       stride - vector stride (assumed the same for each vector)
!            e.g. the j-th element of u is u( 1+(j-1)*stride )
!       nddu - number of histogram nodes in u
!
!  output: - all output must be passed as input in subsequent calls
!       bldd,brdd - range of histogram
!       hist(nddu) - histogram
!
!  determination of range:
!       if, on the first call, the number of samples n is greater
!       than zero, then the range [bldd,brdd] is determined from
!       the samples.  the same range is used on subsequent calls
!       even if some samples fall outside.  such samples are shifted
!       in to the boundary.
!       if, on the other hand, you want to specify the range, then
!       in the first call (ifc=1) set n=0, and specify bldd and brdd.
!
!  nodes:
!       the nodes are equally spaced, the first and last being at
!       bldd and brdd.
!
!  normalization:
!       on the final call (ilc .ne. 0) the histogram is normalised to
!       approximate a probability density function.
!
!  diagnostics:
!       warnings are printed on unit lue if iwarn is set to 1
!
        integer stride
        dimension u(1),wt(1)
        dimension hist(nddu)
        character*(*) fmt
!
        data lue,iwarn/0,1/
!
        ilast = 1 + (n-1)*stride
!
!  initialization
!
        if( ifc .ne. 1 ) go to 100
!
        do 10 ku=1,nddu
10      hist(ku)=0.
!
        if( n .eq. 0 ) return
!
!  find bldd and brdd
!
        bldd= 1.e25
        brdd=-1.e25
!
        do 30 i=1,ilast,stride
        bldd = amin1( u(i) , bldd )
30      brdd = amax1( u(i) , brdd )
!
100     continue
!
!  check ranges
!
	if (ifc.eq.1) then
        if( bldd .gt. brdd ) then
           write(lue,*)' error in hist1: negative range '
           write(lue,*)' bldd,brdd=',bldd,brdd
           write(lue,*)' execution terminated '
        elseif( bldd .eq. brdd  .and. iwarn .eq. 1 ) then
           write(lue,*)' warning in hist1: zero range '
           write(lue,*)' bldd=brdd=',bldd
        endif
	end if
!
!  form histogram
!
        umult=(nddu-1)/amax1( 1.e-10 , brdd-bldd )
!
        if( kwt .eq. 0 ) then
           do 200 i=1,ilast,stride
           ku = ( u(i)-bldd )*umult + 1.5
           ku = max0( 1 , min0( ku,nddu )  )
200        hist(ku) = hist(ku) + wt(i)
        else
           do 210 i=1,ilast,stride
           ku = ( u(i)-bldd )*umult + 1.5
           ku = max0( 1 , min0( ku,nddu )  )
210        hist(ku) = hist(ku) + 1.
        endif
!
!  on last call, normalize
!
        if( ilc .eq. 0 ) return
!
           sum=0.
           do 310 ku=1,nddu
310        sum=sum+hist(ku)
!
           if( sum .ne. 0. ) then
              sum=umult/sum
              do 320 ku=1,nddu
320           hist(ku)=sum*hist(ku)
           endif
!
           hist(1   )=2.*hist(1   )
           hist(nddu)=2.*hist(nddu)
!
        if( ilc .eq.-1 ) return
!
!  write out data
!
      lu=iabs(ilc)
!
      write (lu,201) 1,0
      write (lu,202) nddu,bldd,brdd
      write (lu,fmt) hist
!
 201  format (2i3)
 202  format (i5,1p,2e12.4)
 203  format (10f7.4)
 204  format (1p,6e12.4)
!
        return
        end
