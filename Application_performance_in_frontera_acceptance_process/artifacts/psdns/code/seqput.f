        subroutine seqput(n,iseqpt)
c
c  routine to set seeds
c  arguments:
c       n      - number of seeds
c       iseqpt - seeds
c
        parameter (nk=20)
        common/seqcom/iseq(nk),klast
        dimension iseqpt(n)
        save ifst
        data ifst/0/
c
        if( n .gt. nk ) then
           write(6,*)' nk=',nk,' must be .ge. n=',n
           write(6,*)' stopped in seqput '
           stop
        else if( n .gt. nk ) then
           write(6,*)' only ',n,' out of ',nk,' seeds set in seqput '
        endif
c
        ifst = 1
        do 10 i=1,n
        if( iseqpt(i) .gt. 0 ) then
           iseq(i) = iseqpt(i)
        else
           write(6,*)' invalid seed iseqpt(',i,') =',iseqpt(i)
           write(6,*)' stopped in seqput '
           stop
        endif
10      continue
c
        return
c
c       entry seqset(n,iseqpt)
        entry seqset(n,iset)
c
c       iseqpt(1) = ifst
        iset = ifst
c
        return
        end
