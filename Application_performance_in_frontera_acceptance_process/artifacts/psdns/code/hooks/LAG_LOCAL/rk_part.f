       subroutine rk_part (ispc)
c
c more general particle tracking call, for a given value of kstep
c
#ifdef LAG
c
	use compart
	use lag_timers
	implicit none
	include 'intvars'
c
	real slant,umax
	integer ispc,k,jj,j
	real(8) rtime1,rtime2
	integer indx, ir, nprec,ip
        integer igm
c
	real*8 cpu1,cpu2
	save cpu1,cpu2
c
c generate cubic spline coefficients for velocity components
c
	if(iout.eq.0) then
	jj=1
	else
	jj=2
	end if
c
	rtime1 = MPI_WTIME()
	if (ispc.eq.1) then
c
	cpu1=0.
	cpu2=0.


#ifdef SPXYZM
       call spxyz_m (u(1,1,1,1),bs(1,1,1,1),3)
#else
       call spxyz_m (u(1,1,1,1),bs(1,1,1,1),1)
       call spxyz_m (u(1,1,1,2),bs(1,1,1,2),1)
       call spxyz_m (u(1,1,1,3),bs(1,1,1,3),1)
#endif
	call bsminmaxv ('rk_part, ispc=1',bs(1,1,1,1),1)
	call bsminmaxv ('rk_part, ispc=1',bs(1,1,1,2),1)
	call bsminmaxv ('rk_part, ispc=1',bs(1,1,1,3),1)
	else
	call bsminmaxv ('rk_part, ispc=0',bs(1,1,1,1),1)
	call bsminmaxv ('rk_part, ispc=0',bs(1,1,1,2),1)
	call bsminmaxv ('rk_part, ispc=0',bs(1,1,1,3),1)
	end if

	cpu1 = cpu1 + MPI_WTIME() - rtime1
	cpu_spxyz(0,jj) = cpu_spxyz(0,jj) + MPI_WTIME() - rtime1
	

	
	ps(1)=sqrt(b11(2))
	ps(2)=sqrt(b22(2))
	ps(3)=sqrt(b33(2))
c
c carry out actual interpolation operations
c
        slant=-b12(2)*ps(2)/ps(1)*nx/ny
c
	rtime2 = MPI_WTIME()

#ifdef NOV15
#else

	if (ndpart.gt.0) then

	call MPI_BARRIER (MPI_COMM_WORLD, mpierr)

      call ppminmaxv ('before local',pp(1,4),ndpart,nump1)
      call spcal_local (nxyz, gsh(1,1), bs, pp(1,1), pp(1,4),
     1                  nump1, 3, ps(1), ndpart, iod(1,1), 1)
      call ppminmaxv ('after local',pp(1,4),ndpart,nump1)

	endif  ! if (nop.gt.0)

#ifdef LAGSC2
	if (ndpart2.gt.0) then


	call MPI_BARRIER (MPI_COMM_WORLD, mpierr)

      call ppminmaxv ('before local',pp2(1,4),ndpart2,nump2)
      call spcal_local (nxyz, gsh(1,1), bs, pp2(1,1), pp2(1,4),
     1                  nump2, 3, ps(1), ndpart2, iod(1,1), 2)
      call ppminmaxv ('after local',pp2(1,4),ndpart2,nump2)


	endif
#endif
!LAGSC2 directive


#endif
! NOV15 directive



#ifdef MOL
	if (nom.gt.0) then

	nprec = nom/nmrec/ngm
	allocate (bfpos(3,nprec))

        do igm=1,ngm
	do ir=1,nmrec

	indx = (ir-1)*nprec/numtasks + 1

        call intbf_ars ( xyzl,nxyz,gsh(1,1),mpp(1,1,igm),indx,nmrec,nom/ngm,
     1               bfpos,slant, 3)

#if defined(CAF) && defined(CF_LAG)
        call spcal_tree ( bfpos, bs(1,1,1,1),  mpp(1,4,igm), 3, ps(1),
     1             indx, nmrec, nom/ngm, iod(1,1))
#else
         call spcal_ars ( bfpos, bs(1,1,1,1),  mpp(1,4,igm), 3, ps(1),
     1             indx, nmrec, nom/ngm, iod(1,1), 3)
#endif

	enddo
	enddo

	deallocate (bfpos)

	end if  ! if(nom.gt.0)
#endif
! MOL directive
c
	cpu2 = cpu2 + MPI_WTIME() - rtime2
c
	if (ispc.eq.0) then
	cpu_rkpart(1,jj) = cpu_rkpart(1,jj) + cpu1
	cpu_rkpart(2,jj) = cpu_rkpart(2,jj) + cpu2
	if (taskid.eq.0) then
	write (6,"('cpu_rkpart:',i6,i3,1p,3e12.4)") istep,jj,sum(cpu_intbf(4,1:2,jj)),sum(cpu_spcal(4,1:2,jj)),cpu_rkpart(2,jj)
	end if
	end if
	
c
#endif
c

      return
      end
