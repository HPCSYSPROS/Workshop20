        subroutine phymsq (vx,meansq)
                                                                                
        use com
        implicit none
	include 'intvars'
                                                                                
        real(b8) :: vx(nx,zisz,yjsz)
	  real*8 :: sum1
        real(b8) :: sum2,ee,meansq
                                                                                
        sum1=0.
        do yp=1,yjsz
        do zp=1,zisz
        do x=1,nx
          sum1=sum1+vx(x,zp,yp)*vx(x,zp,yp)
        enddo
        enddo
        enddo
        sum1=sum1/nx/yjsz/zisz
	  sum2=sum1
                                                                                
        call MPI_ALLREDUCE (sum2,ee,1,mpireal,MPI_SUM,
     1 	            MPI_COMM_WORLD,ierr)
	meansq=ee/numtasks
c
        return
        end

