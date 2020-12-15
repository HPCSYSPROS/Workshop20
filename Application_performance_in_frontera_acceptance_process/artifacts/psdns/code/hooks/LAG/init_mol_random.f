        subroutine init_mol_random

#ifdef LAG
#ifdef MOL
#ifndef MOLOLD

        use compart
        implicit none
        integer i, seedsize,tseed,ttask, inew
        integer, allocatable :: seed(:),allseeds(:)


        if(taskid.eq.0) write(6,*) 'inside init_mol_random'

        if(mstart.eq.0) then

        call random_seed(size=seedsize)
        allocate (mseed(seedsize))


        do i=1,seedsize
        mseed(i) = (taskid+1)*minitseed*i
        enddo
        call random_seed(put=mseed)


        elseif (mstart.gt.0) then

	inew=0
        call random_seed(size=seedsize)
        allocate (mseed(seedsize))

        if(taskid.eq.0) then

        open(2300,file='molseedin')
        read(2300,*) tseed,ttask

        if(tseed.ne.seedsize) then
        write(6,*) 'warning: seedsize in molseedin dont match'
        write(6,*) 'setting new seedsize and seeds'
        inew = 1
        endif
        if(ttask.ne.numtasks) then
        write(6,*) 'warning: numtasks in molseedin dont match'
        write(6,*) 'setting new seedsize and seeds'
        inew = 1
        endif

        endif

        call MPI_BCAST (inew,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)

	if(inew.eq.0) then

        allocate(allseeds(seedsize*numtasks))
	if(taskid.eq.0) read(2300,100) allseeds


        call MPI_SCATTER (allseeds,seedsize,MPI_INTEGER,mseed,seedsize,
     1                 MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)



        call random_seed(put=mseed)

        deallocate(allseeds)

	else

        do i=1,seedsize
        mseed(i) = (taskid+1)*(minitseed+10)*i 
        enddo
        call random_seed(put=mseed)

        endif

	endif
#endif
#endif
#endif

100     format( (5i15) )

        return
        end

