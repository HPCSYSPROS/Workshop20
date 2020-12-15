c     $Header$

#include "cctk.h"

      subroutine apply_dissipation (var, rhs, ni, nj, nk, dx, order, epsdis)
      implicit none
      
      integer   ni, nj, nk
      CCTK_REAL var(ni,nj,nk), rhs(ni,nj,nk)
      CCTK_REAL dx(3), idx(3)
      CCTK_INT  order
      CCTK_REAL epsdis(ni,nj,nk)
      
      integer   i, j, k
      
      idx(1) = 1 / dx(1)
      idx(2) = 1 / dx(2)
      idx(3) = 1 / dx(3)
      
      if (order .eq. 1) then
         
c$omp parallel do private(i,j,k)
         do k = 2, nk-1
            do j = 2, nj-1
               do i = 2, ni-1
                  rhs(i,j,k) = rhs(i,j,k) + epsdis(i,j,k)
     $                 * (+ (var(i-1,j,k) - 2*var(i,j,k) + var(i+1,j,k)) * idx(1)
     $                    + (var(i,j-1,k) - 2*var(i,j,k) + var(i,j+1,k)) * idx(2)
     $                    + (var(i,j,k-1) - 2*var(i,j,k) + var(i,j,k+1)) * idx(3))
                  
               end do
            end do
         end do
         
      else if (order .eq. 3) then
         
c$omp parallel do private(i,j,k)
         do k = 3, nk-2
            do j = 3, nj-2
               do i = 3, ni-2
                  
                  rhs(i,j,k) = rhs(i,j,k) - epsdis(i,j,k) / 16
     $                 * (+ (var(i-2,j,k) - 4*var(i-1,j,k) + 6*var(i,j,k) - 4*var(i+1,j,k) + var(i+2,j,k)) * idx(1)
     $                    + (var(i,j-2,k) - 4*var(i,j-1,k) + 6*var(i,j,k) - 4*var(i,j+1,k) + var(i,j+2,k)) * idx(2)
     $                    + (var(i,j,k-2) - 4*var(i,j,k-1) + 6*var(i,j,k) - 4*var(i,j,k+1) + var(i,j,k+2)) * idx(3))
                  
               end do
            end do
         end do
         
      else if (order .eq. 5) then
         
c$omp parallel do private(i,j,k)
        do k = 4, nk-3
          do j = 4, nj-3
            do i = 4, ni-3
               
              rhs(i,j,k) = rhs(i,j,k) + epsdis(i,j,k) / 64
     $           * (+ (var(i-3,j,k) - 6*var(i-2,j,k) + 15*var(i-1,j,k) - 20*var(i,j,k) + 15*var(i+1,j,k) - 6*var(i+2,j,k) + var(i+3,j,k)) * idx(1)
     $              + (var(i,j-3,k) - 6*var(i,j-2,k) + 15*var(i,j-1,k) - 20*var(i,j,k) + 15*var(i,j+1,k) - 6*var(i,j+2,k) + var(i,j+3,k) ) * idx(2)
     $              + (var(i,j,k-3) - 6*var(i,j,k-2) + 15*var(i,j,k-1) - 20*var(i,j,k) + 15*var(i,j,k+1) - 6*var(i,j,k+2) + var(i,j,k+3) ) * idx(3))
              
            end do
          end do
        end do
        
      else if (order .eq. 7) then
         
c$omp parallel do private(i,j,k)
        do k = 5, nk-4
          do j = 5, nj-4
            do i = 5, ni-4
               
              rhs(i,j,k) = rhs(i,j,k) - epsdis(i,j,k) / 256
     $           * (+ (var(i-4,j,k) - 8*var(i-3,j,k) + 28*var(i-2,j,k) - 56*var(i-1,j,k) + 70*var(i,j,k) - 56*var(i+1,j,k) + 28*var(i+2,j,k) - 8*var(i+3,j,k) + var(i+4,j,k)) * idx(1)
     $              + (var(i,j-4,k) - 8*var(i,j-3,k) + 28*var(i,j-2,k) - 56*var(i,j-1,k) + 70*var(i,j,k) - 56*var(i,j+1,k) + 28*var(i,j+2,k) - 8*var(i,j+3,k) + var(i,j+4,k)) * idx(2)
     $              + (var(i,j,k-4) - 8*var(i,j,k-3) + 28*var(i,j,k-2) - 56*var(i,j,k-1) + 70*var(i,j,k) - 56*var(i,j,k+1) + 28*var(i,j,k+2) - 8*var(i,j,k+3) + var(i,j,k+4)) * idx(3))
              
            end do
          end do
        end do
        
      else if (order .eq. 9) then
         
c$omp parallel do private(i,j,k)
        do k = 6, nk-5
          do j = 6, nj-5
            do i = 6, ni-5
               
              rhs(i,j,k) = rhs(i,j,k) + epsdis(i,j,k) / 1024
     $           * (+ (var(i-5,j,k) - 10*var(i-4,j,k) + 45*var(i-3,j,k) - 120*var(i-2,j,k) + 210*var(i-1,j,k) - 252*var(i,j,k) + 210*var(i+1,j,k) - 120*var(i+2,j,k) + 45*var(i+3,j,k) - 10*var(i+4,j,k) + var(i+5,j,k)) * idx(1)
     $              + (var(i,j-5,k) - 10*var(i,j-4,k) + 45*var(i,j-3,k) - 120*var(i,j-2,k) + 210*var(i,j-1,k) - 252*var(i,j,k) + 210*var(i,j+1,k) - 120*var(i,j+2,k) + 45*var(i,j+3,k) - 10*var(i,j+4,k) + var(i,j+5,k)) * idx(2)
     $              + (var(i,j,k-5) - 10*var(i,j,k-4) + 45*var(i,j,k-3) - 120*var(i,j,k-2) + 210*var(i,j,k-1) - 252*var(i,j,k) + 210*var(i,j,k+1) - 120*var(i,j,k+2) + 45*var(i,j,k+3) - 10*var(i,j,k+4) + var(i,j,k+5)) * idx(3))
              
            end do
          end do
        end do
        
      else
         
         call CCTK_WARN (0, "internal error")
         
      end if
      
      end
