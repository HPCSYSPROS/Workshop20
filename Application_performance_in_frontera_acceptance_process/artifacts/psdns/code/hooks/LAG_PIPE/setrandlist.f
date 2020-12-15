
        subroutine setrandlist()

#ifdef LAG
#ifdef CF_LAG

        use mpilag

        integer i
        real*8 x

        allocate (rantasks(0:numtasks-1))

        do i=0,numtasks-1
        rantasks(i) = i
        enddo
        do i=0,numtasks-1
        call random_number (x)
        i1 = max(0,min(floor(x*numtasks),numtasks-1))
        i2 = rantasks(i1)
        rantasks(i1)=rantasks(i)
        rantasks(i)=i2
        enddo

#endif
#endif

        return

        end

