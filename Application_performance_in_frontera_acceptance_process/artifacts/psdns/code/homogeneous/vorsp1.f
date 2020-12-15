! computation of 3-d vorticity spectrum tensor (like sub. sptr1)
! called from sub. sptr1
! note that in the DNS code we have u1/b11, etc. 
! this revised version takes these normalization factors into account
 
      subroutine vorsp1 (uy,m,ithr,bk2y) 
	use comsp
	implicit none
	include 'intvars'
 
	complex(b8) uy(ny,zjsz*xisz,3) 
!	complex(b8), allocatable :: bk2y(:)
	complex(b8) bk2y(ny)
	complex(b8) bk3i,bk1i,v1,v2,v3
        integer a,get_l
c
c the following declarations added, P.K Yeung 2/2/10
c also made fixes for non-2pi^3 domains
c
	integer m,ik
	real s1,s2,s3,rk2
	complex vt1,v3c
c
        real(b8) beta_min
c
	integer ithr

      beta_min=min(beta1,beta2,beta3)


!        allocate(bk2y(ny))
c
c

        kx2(:)=b11(m)*kx(:)**2
        ky2(:)=b22(m)*ky(:)**2
        kz2(:)=b33(m)*kz(:)**2
 
! define the metric normalization factors
 
        s1=sqrt(b11(m))
        s2=sqrt(b22(m))
        s3=sqrt(b33(m))
 
 
      do y=1,ny
         bk2y(y)=imagi*s2*ky(y)
      end do
c
      xp = ia_xst(ithr)
      z = ia_zst(ithr)
c
      do 5 a = ia_st(ithr), ia_st(ithr)+ia_sz(ithr)-1

         x = xp + xist-1

         bk1i=imagi*s1*kx(x)
         
         bk3i=imagi*s3*kz(z)            
         do  y=1,ny                  
            if (mask(y,a)) then
                     
               rk2=kx2(x)+ky2(y)+kz2(z)
               ik=sqrt(rk2)/beta_min+1.5
               
               v1=bk2y(y)*uy(y,a,3)*s3-bk3i*uy(y,a,2)*s2
               v2=bk3i*uy(y,a,1)*s1-bk1i*uy(y,a,3)*s3
               v3=bk1i*uy(y,a,2)*s2-bk2y(y)*uy(y,a,1)*s1
               
               v3c=conjg(v3)
               vt1=tfact(x)*v1
               
               vijky(ik,1,ithr)=vijky(ik,1,ithr)+vt1*conjg(v1)
               vijky(ik,2,ithr)=vijky(ik,2,ithr)+vt1*conjg(v2)
               vijky(ik,3,ithr)=vijky(ik,3,ithr)+vt1*v3c
               vijky(ik,4,ithr)=vijky(ik,4,ithr)+tfact(x)*v2*conjg(v2)
               vijky(ik,5,ithr)=vijky(ik,5,ithr)+tfact(x)*v2*v3c
               vijky(ik,6,ithr)=vijky(ik,6,ithr)+tfact(x)*v3*v3c
            endif
         enddo

         call next_xz(xp,z)
 5	continue
 
!      deallocate(bk2y)

      return
      end
