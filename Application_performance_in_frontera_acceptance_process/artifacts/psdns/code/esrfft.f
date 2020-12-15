        subroutine esrfft(u,itrans,ivec,ntrans,nvec,isgn)          
#ifdef ESSL	  
c
c CURRENTLY NOT USED
c                                                                       
c essl version of old sub. fft991                                       
c                                                                       
c       routine to perform a real fft on a plane of data.               
c               transform is normalized                                 
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
c                                                                       
        use com
        real(b8) u(1),norm                                                  
c                                                                       
c                                                                       
c  forward transform                                                    
c                                                                       
        if( isgn .eq. 1 ) go to 199                                     
c                                                                       
        norm=1./float( ntrans )                                         
c                                                                       
c  transform                                                            
c                                                                       
c                                                                       
c                                                                       
         call srcft(0,u,ivec,u,ivec/2,ntrans,nvec,-isgn,norm,raux1,     
     1   rnaux1,raux2,rnaux2,raux3,rnaux3)                              
c                                                                       
c                                                                       
        if (nvec.eq.1)                                                  
     1  write (6,*) 'exit esrfft: isgn,nvec=',isgn,nvec                 
        return                                                          
c                                                                       
c                                                                       
c  inverse transform                                                    
c                                                                       
199     continue                                                        
c                                                                       
c  loop over vectors                                                    
c                                                                       
c  transform                                                            
c                                                                       
c                                                                       
          call scrft(0,u,ivec/2,u,ivec,ntrans,nvec,-isgn,1.,raux1,      
     1    rnaux1,raux2,rnaux2,raux3,rnaux3)                             
c                                                                       
c                                                                       
        if (nvec.eq.1)                                                  
     1  write (6,*) 'exit esrfft: isgn,nvec=',isgn,nvec                 
#endif
        return                                                          
        end                                                             
