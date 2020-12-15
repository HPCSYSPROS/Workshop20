      subroutine lvalid
c
c routine to validate the choice of lagrangian parameters
c in the simulation: write messages and stop the program
c if appropriate
c
#ifdef LAG
	use compart
	implicit none
c
c
      if (taskid.eq.0) write (6,*) 'enter lvalid'
c for pstart.eq.1, check that the number of particles specified
c is within the dimensioned limit
c
      if (pstart.eq.1.and.nop.gt.ndpart) then
	if (taskid.eq.0) then
        write (6,*) 'no of particles exceeds dimensioned limit'
        write (6,*) 'nop=',nop,'ndpart=',ndpart
        write (6,*) 'use less particles, or increase ndpart in'
        write (6,*) 'parameter statement, and re-compile'
        write (6,*) 'stopped in lvalid (pstart=1)'
	end if
	call MPI_ABORT (MPI_COMM_WORLD,ierror)
      end if
c
      if (nop2.gt.ndpart2) then
	if (taskid.eq.0) then
        write (6,*) 'no of particles exceeds dimensioned limit'
        write (6,*) 'nop2=',nop2,'ndpart2=',ndpart2
        write (6,*) 'use less particles, or increase ndpart2 in'
        write (6,*) 'parameter statement, and re-compile'
        write (6,*) 'stopped in lvalid'
	end if
	call MPI_ABORT (MPI_COMM_WORLD,ierror)
      end if
c
c check that LGRAD is used if any of the velocity or scalar
c gradients (and Laplacians) are selected
c
#ifdef LGRAD
#else
      do i=1,3+nc
      if (lpgrad(i).gt.0.or.lpsecd(i).gt.0.or.lplapl(i).gt.0) then
      write (6,*) 'LGRAD option must be selected if gradients'
      write (6,*) 'are required following the fluid particles'
      write (6,*) 'i,lpgrad,lpsecd,lplapl=',
     1             i,lpgrad(i),lpsecd(i),lplapl(i)
	call MPI_ABORT (MPI_COMM_WORLD,ierror)
      end if
      end do
#endif
c
! same considerations for molecule groups
c
#ifdef MOL
	if (nom.gt.0) then
      if (ngm.gt.numtasks) then
      write (6,*) 'Abort: ngm should not exceed no. of MPI tasks'
      call MPI_ABORT (MPI_COMM_WORLD,ierror)
	end if
#endif
c
c determine actual no. of properties per particle
c
      npp=6
c
      do 10 i=1,3
      npp=npp+lpgrad(i)*3
 10   continue
      nppv=npp
c
c if all 9 velocity gradients are calculated, dw/dz need not be stored
c and may be overwritten by mechanical dissipation
c
      if (lpgrad(1)*lpgrad(2)*lpgrad(3).eq.1) npp=npp-1
c
      ldudx=7
      ldvdy=ldudx+2*lpgrad(1)+6*lpsecd(1)+2*lpgrad(2)
      ldwdz=ldvdy+1*lpgrad(2)+6*lpsecd(2)+3*lpgrad(3)
      npp=ldwdz
c
c for the case where lpmdis is set, but velocity gradients are not
c to be saved, second spatial derivatives are not allowed
c
      if (lpmdis.eq.1.and.npo(2).eq.0) then
c
      if (lpsecd(1).ne.0.or.lpsecd(2).ne.0.or.lpsecd(3).ne.0) then
        write (6,*) 'mech. diss. computed w/o saving velocity gradients'
        write (6,*) 'but some of lpsecd(1,2,3) are nonzero:',
     1              (lpsecd(i),i=1,3)
        write (6,*) 'this is not allowed'
        write (6,*) 're-specify correct values in lbdata'
        write (6,*) 'stopped in sub. lvalid'
	call MPI_ABORT (MPI_COMM_WORLD,ierror)
      else
      npp=ldudx
      end if
c
      end if
c
c second-derivatives of w (if selected) will follow dwdz
c (or mechanical dissipation)

      npp=npp+lpsecd(3)*6
      npp=npp+lplapl(1)+lplapl(2)+lplapl(3)
c
	write (6,*) 'lvalid, before do 20, npp,ncop=',npp,ncop
c#ifdef OLD
c allowance is made for storage savings if scalar gradients are
c not saved after scalar dissipation is formed
c
      if (ncop.gt.0) then
c
      do 20 i=4,3+ncop
c
c inc1 is no. of storage slots needed in evaluating particle
c scalar properties.      inc2 is  the no. needed to store them for
c later output to disk files
c
c if lpgrad(i)=0 & lpsdis(i)=1, then lpsecd(i)=1 is not allowed
c
      inc1=0
      inc2=lpfunc(i)+lpgrad(i)*3+lpsecd(i)*6+lpsdis(i)+lplapl(i)
      write (6,*) 'lvalid: i,inc2=',i,inc2
      if (.not.(lpsdis(i).eq.1.and.lpgrad(i).eq.0)) go to 25
c
      inc1=lpfunc(i)+3
      write (6,*) 'lvalid: i,inc1,inc2=',i,inc1,inc2
      if (lpsecd(i).ne.0) then
         write (6,*) 'the combination (lpgrad,lpsdis,lpsecd)=(0,1,1)'
         write (6,*) 'is not allowed for scalars'
         write (6,*) 're-specify correct values in lbdata'
         write (6,*) 'stopped in sub. lvalid'
	call MPI_ABORT (MPI_COMM_WORLD,ierror)
      end if
c
 25   npp=npp+inc2
c
 20   continue
c
      npp=max0(npp,npp+inc1-inc2)
c
      else
      npp=max0(npp,nppv)
c
      end if
c#endif
	write (6,*) 'lvalid,  after do 20, npp=',npp
c
c note that npp must be at least 9, to allow for temporary
c storage needed in updating particle displacements
c
      npp=max0(npp,9)
c
c if velocity gradients are computed, npp must be at least 15
c
      if (lpmdis.eq.1) npp=max0(npp,15)
c
c check whether npp is within the dimensioned limit
c
#ifndef LAGSC2
      if (npp.gt.ndprop) then
	if (taskid.eq.0) then
        write (6,*) 'no. of properties per particle is', npp, 'exceeds'
        write (6,*) 'the dimensioned limit of', ndprop, 'increase'
        write (6,*) 'ndprop in parameter statement and re-compile'
        write (6,*) 'stopped in lvalid'
	end if
	call MPI_ABORT (MPI_COMM_WORLD,ierror)
      else
      if (taskid.eq.0) 
     1  write (6,*) 'storage slots needed per particle is',npp,
     1              '(ok with ndprop=',ndprop,')'
      end if
#endif
c
#ifdef LAGSC2
c
      if (ndprop.lt.9) then
	if (taskid.eq.0) then
        write (6,*) 'no. of properties per particle is', 9, 'exceeds'
        write (6,*) 'the dimensioned limit of', ndprop, 'increase'
        write (6,*) 'ndprop in parameter statement and re-compile'
        write (6,*) 'stopped in lvalid'
	end if
	call MPI_ABORT (MPI_COMM_WORLD,ierror)
      else
      if (taskid.eq.0) 
     1  write (6,*) 'storage slots needed per particle is',9,
     1              '(ok with ndprop=',ndprop,')'
      end if
c
      if (npp.gt.ndprop2) then
	if (taskid.eq.0) then
        write (6,*) 'no. of properties per particle (2nd ensemble) is',
     1              npp, 'exceeds'
        write (6,*) 'the dimensioned limit of', ndprop2, 'increase'
        write (6,*) 'ndprop2 in parameter statement and re-compile'
        write (6,*) 'stopped in lvalid'
	end if
	call MPI_ABORT (MPI_COMM_WORLD,ierror)
      else
      if (taskid.eq.0) 
     1  write (6,*) 'storage slots needed per particle (2nd ensemble)',
     1              'is',npp,'(ok with ndprop2=',ndprop2,')'
      end if
c
#endif
c
#ifdef OLD
c for scalars, if lpfunc is 1, then so must kcps (otherwise
c the scalar would not be obtained in physical space in
c the eulerian code)
c
c also, kcps(i) should be zero for i.gt.ncop
c
      do 30 i=1,ncop
      j=i+3
c
      if (lpfunc(j).eq.1 .and. kcps(i).eq.0 ) then
        write (6,*) 'the',i,'-th scalar is required in physical space'
        write (6,*) 'by the lagrangian code: kcps(i) should be set to 1'
        write (6,*) 're-set kcps array in ebdata'
        write (6,*) 'also check the value of ncop in yparam'
        write (6,*) 'stopped in lvalid'
	call MPI_ABORT (MPI_COMM_WORLD,ierror)
      end if
c
      if (kcps(i).eq.1 .and. i.gt.ncop) then
        write (6,*) 'kcps array and value of ncop in conflict'
        write (6,*) 'kcps(i) should be zero for i.gt.ncop'
        write (6,*) 'stopped in lvalid'
	call MPI_ABORT (MPI_COMM_WORLD,ierror)
      end if
c
 30   continue
c
c check whether nbg and neg are within dimensioned limits
c
      if ( nbg.gt.nbgdim .or. neg.gt.negdim ) then
        write (6,*) 'nbg and/or neg exceeds dimensioned limit'
        write (6,*) 'increase nbgdim and/or negdim in ycompart'
        write (6,*) 'stopped in lvalid'
	call MPI_ABORT (MPI_COMM_WORLD,ierror)
      end if
c
#endif
c
c check that random no. sequence used for generating initial particle
c positions are distinct from sequences used in eulerian code
c kipran also should not be 1, since the 1st sequence is reserved
c for generating velocity fields in wavenumber space
c
      if (pstart.ne.0) go to 49
#ifndef NOSCALAR
      do 40 i=1,nc
      if (kipran.eq.kran(i)) go to 45
 40   continue
#endif
      if (kipran.eq.ksran.or.kipran.eq.kranf1.or.kipran.eq.kranf2)
     1    go to 45
      go to 49
c
 45   write (6,*) 'warning from sub. lvalid:'
      write (6,*) 'kipran is the same as one of'
      write (6,*) 'kran,ksran,kranf1,kranf2, or 1'
      write (6,*) 'this is not proper/recommended usage'
c
 49   continue
c
c check spline parameters
c
      if (meth.eq.2.or.meth.eq.3) then
      if (ndb.ne.max0(nx,ny,nz)) then
      write (6,*) 'ndb should be equal to max(mx,ny,nz) in spvcom'
      write (6,*) 'program aborts in sub. lvalid'
	call MPI_ABORT (MPI_COMM_WORLD,ierror)
      end if
      end if
c
c ngp should be 1 if pstart=0 is used to initialize particles
c
      if (ngp.gt.1.and.pstart.eq.0) then
      write (6,*) 'ngp should be set to 1 if pstart=0:', ngp
      write (6,*) 'modify lbdata.f and try again'
	call MPI_ABORT (MPI_COMM_WORLD,ierror)
      end if
c
c for shear flow, ngp at present should not exceed numtasks
c (since each node will write output for only a single
c  group of particles)
c
      if (shear.gt.0..and.ngp.gt.numtasks) then
      write (6,*) 'ngp should not exceed no. of nodes'
      write (6,*) 'modify param/lbdata.f and try again'
	call MPI_ABORT (MPI_COMM_WORLD,ierror)
      end if
c
c
c itetra must be >0 if pstart=-8
c
      if (pstart.eq.-8.and.itetra.eq.0) then
      write (6,*) 'lvalid: itetra must be >0 if pstart=-8'
	call MPI_ABORT (MPI_COMM_WORLD,ierror)
      end if

#endif
      if (taskid.eq.0) write (6,*) ' exit lvalid'
      return
      end
