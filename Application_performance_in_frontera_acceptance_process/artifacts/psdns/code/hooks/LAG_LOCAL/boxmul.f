      subroutine boxmul(x,n)

	use mpilag, only : p8
        implicit none
        integer n,i
        real(p8) :: x(n), R, theta, u1,u2,pi

        pi=atan(1.)*4.
c
! n must be a multiple of 2
        if(mod(n,2).ne.0) then
	write(6,*) 'stop in boxmul, n must be a multiple of 2'
	stop
	endif

        do i=1,n,2

        u1 = x(i)
        u2 = x(i+1)

        R = sqrt(-2.*log(u1))
        theta = 2*pi*u2

        x(i) = R*cos(theta)
        x(i+1) = R*sin(theta)


        enddo

        return
        end

