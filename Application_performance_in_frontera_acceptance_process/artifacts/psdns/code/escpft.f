        subroutine escpft(u,itrans,ivec,ntrans,nvec0,isgn,ithr)               
#ifdef ESSL
#include "fft_stuff.f"
c                                                                       
c essl library version                                                  
c this call to set up auxiliary arrays for complex fft                  
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
c       nvec0            number of vectors                               
c       isgn            = -1 for forward transform, =+1 for inverse     
c                                                                       
	use com
        complex(b8) u(1)                                                    
c                                                                       
c isgn=-1 for forward transform                                         
c isgn= 1 for inverse transform                                         
c                                                                       
        if (nvec0.eq.0) return
        call cft(1,u,itrans,ivec,u,itrans,ivec,ntrans,nvec0,-isgn,1.0,  
     1   caux1(1,ithr),cnaux,caux2(1,ithr),cnaux)                              
c                                                                       
#endif
        return                                                          
        end                                                             
