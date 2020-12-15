        subroutine escfft(u,itrans,ivec,ntrans,nvec,isgn,ithr)          
#ifdef ESSL	 
#include "fft_stuff.f" 
c                           
c essl library version of old sub. cfft99                               
c                  
c       routine to perform a complex fft on a plane of data.            
c               transform is unnormalised                               
c                                      
c       u               complex array containing data to be transformed.
c                       u(1) contains the first element of the first vec
c                       the transform is returned in u.                 
c       itrans  transform stride                    
c       ivec            vector stride               
c       ntrans  length of transform                 
c       nvec            number of vectors           
c       isgn            = -1 for forward transform, =+1 for inverse     
c       work            complex work array of dimension ntrans          
c                                                   
	use com
        complex(b8) u(1)                                
c                                                   
c isgn=-1 for forward transform                     
c isgn= 1 for inverse transform                     
c                                                   
        if (nvec.eq.0) return
        call cft(0,u,itrans,ivec,u,itrans,ivec,ntrans,nvec,-isgn,1.0,  
     1   caux1(1,ithr),cnaux,caux2(1,ithr),cnaux)                 
c                                                   
#endif
        return                                      
        end                                         
