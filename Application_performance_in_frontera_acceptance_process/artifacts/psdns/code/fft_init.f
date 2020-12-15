!!--------------------------------------------------------------
c One-time initialization
C Initialize work arrays for ESSL

      subroutine init_work
      
      use com
      implicit none

      integer err,ithr,omp_get_thread_num

#ifdef ESSL

#ifndef SINGLE_PREC
      integer nyz
      nyz = max(ny,nz)
      if(nyz .le. 2048) then
         cnaux = 20000
      else 
         cnaux = 20000+2.28*nyz
      endif
      if(nyz .ge. 252) then
         cnaux = cnaux+(2*nyz+256)*64
      endif
      
      if(nx .le. 4096) then
         rnaux1 = 22000
         rnaux2 = 20000
      else
         rnaux1 = 20000+1.64*nx
         rnaux2=20000+1.14*nx
      endif
#else
      integer nyz
      nyz = max(ny,nz)
      if(nyz .le. 8192) then
         cnaux = 20000
      else 
         cnaux = 20000+1.14*nyz
      endif
      if(nyz .ge. 252) then
         cnaux = cnaux +(nyz+256)*64
      endif

      if(nx .le. 16384) then
         rnaux1 = 25000
         rnaux2 = 20000
      else
         rnaux1 = 20000+0.82*nx
         rnaux2=20000+0.57*nx
      endif
#endif

      
c      print *,'In init_work, num_thr,cnaux=',num_thr,cnaux
      allocate(caux1(cnaux,0:num_thr-1),stat=err)
      if(err .ne. 0) then
         print *,taskid,': Error allocating caux1(',cnaux
      endif
c      print *,'In init_work, allocating caux2'
      allocate(caux2(cnaux,0:num_thr-1),stat=err)
      if(err .ne. 0) then
         print *,taskid,': Error allocating caux2(',cnaux
      endif
c      print *,'In init_work, allocating raux1'
      allocate(raux1(rnaux1,0:num_thr-1),stat=err)
      if(err .ne. 0) then
         print *,taskid,': Error allocating raux1(',rnaux1
      endif
c      print *,'In init_work, allocating raux2'
      allocate(raux2(rnaux2,0:num_thr-1),stat=err)
      if(err .ne. 0) then
         print *,taskid,': Error allocating raux2(',rnaux2
      endif

#endif

      return
      end

!!--------------------------------------------------------------
      
      subroutine free_work

#ifdef ESSL
c Release work arrays for ESSL

      use com

      deallocate(caux1)      
      deallocate(caux2)
      deallocate(raux1)
      deallocate(raux2)
#endif

      return
      end
