#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"



subroutine PNSMapper_ReadProfile(CCTK_ARGUMENTS)

  use EOS_Omni_Module
  use mpi

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  character(len=512) filename
  integer filenamelength
  integer n,i,myiostat,ibuffer
  real*8 buffer
  character(len=128) sbuffer
  
  real(8) :: pvel_temp(Num_Max_Radial)
  integer :: idx,j,nl,nu
  integer :: rank, ierr

  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

  ! initialize
  n = Num_Max_Radial
  pradius(:) = 0.0d0
  prho(:) = 0.0d0
  ptemp(:) = 0.0d0
  pvel(:) = 0.0d0
  pye(:) = 0.0d0
  ppsi(:) = 1.0d0


  filenamelength = 512
  ! Let's figure out the profile filename
  call CCTK_FortranString(filenamelength,Profile_File,filename)
  
  call CCTK_INFO("Reading PNS Profile")

  if(CCTK_EQUALS(Profile_Type,"GR1Dspecial")) then
     if(rank .eq. 0) then
        open(666,file=trim(filename),status='unknown',form='formatted',action='read')
        ! first line is the profile length
        read(666,*) ibuffer
        myiostat = 0
        i=0
        do while(myiostat .eq. 0) 
           i=i+1
           if(i.gt.Num_Max_Radial) then
              call CCTK_WARN(0,"PNS profile file longer than Num_Max_Radial parameter")
           endif
           read(666,*,iostat=myiostat) pradius(i),prho(i),ptemp(i),pye(i),pvel(i),ppsi(i),buffer
        enddo
        close(666)
     end if
     call MPI_Bcast(i, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(pradius, Num_Max_Radial, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(prho, Num_Max_Radial, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(ptemp, Num_Max_Radial, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(pye, Num_Max_Radial, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(pvel, Num_Max_Radial, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(ppsi, Num_Max_Radial, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     n = i-1
     pnszones = n


     ! radius is in cm
     ! density is in g/ccm
     ! temperature is in MeV
     ! velocity is in cm/s
     pradius = pradius * length_gf
     prho = prho * rho_gf
     pvel = pvel / clite

  else if(CCTK_EQUALS(Profile_Type,"GR1Dshort")) then
     if(rank .eq. 0) then
        open(666,file=trim(filename),status='unknown',form='formatted',action='read')
        myiostat = 0
        i=0
        do while(myiostat .eq. 0) 
           i=i+1
           if(i.gt.Num_Max_Radial) then
              call CCTK_WARN(0,"PNS profile file longer than Num_Max_Radial parameter")
           endif
           read(666,*,iostat=myiostat) buffer,buffer,buffer,ptemp(i),prho(i),pye(i),&
                pvel(i),buffer,buffer,buffer,buffer,pradius(i)
        enddo
        close(666)
     end if
     call MPI_Bcast(i, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(pradius, Num_Max_Radial, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(prho, Num_Max_Radial, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(ptemp, Num_Max_Radial, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(pye, Num_Max_Radial, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(pvel, Num_Max_Radial, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     n = i-1
     pnszones = n
     ! radius is in cm
     ! density is in g/ccm
     ! temperature is in MeV
     ! velocity is in cm/s
     pradius = pradius * length_gf
     prho = prho * rho_gf
     pvel = pvel / clite

  else if(CCTK_EQUALS(Profile_Type,"GR1Dspecial2")) then
     if(rank .eq. 0) then
        open(666,file=trim(filename),status='unknown',form='formatted',action='read')
        myiostat = 0
        i=0
        do while(myiostat .eq. 0) 
           i=i+1
           if(i.gt.Num_Max_Radial) then
              call CCTK_WARN(0,"PNS profile file longer than Num_Max_Radial parameter")
           endif
           read(666,*,iostat=myiostat) buffer,buffer,buffer,ptemp(i),prho(i),pvel(i),pye(i),&
                buffer,buffer,buffer,buffer,pradius(i)
        enddo
        close(666)
     end if
     call MPI_Bcast(i, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(pradius, Num_Max_Radial, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(prho, Num_Max_Radial, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(ptemp, Num_Max_Radial, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(pye, Num_Max_Radial, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(pvel, Num_Max_Radial, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     n = i-1
     pnszones = n
     ! radius is in cm
     ! density is in g/ccm
     ! temperature is in K
     ! velocity is in cm/s
     pradius = pradius * length_gf
     ptemp = ptemp * 8.6173324d-11
     prho = prho * rho_gf
     pvel = pvel / clite

  else if(CCTK_EQUALS(Profile_Type,"GR1Dformat2")) then
     if(rank .eq. 0) then
        open(666,file=trim(filename),status='unknown',form='formatted',action='read')
        myiostat = 0
        i=0
        ! read header
        read(666,*) sbuffer
        do while(myiostat .eq. 0) 
           i=i+1
           if(i.gt.Num_Max_Radial) then
              call CCTK_WARN(0,"PNS profile file longer than Num_Max_Radial parameter")
           endif
           read(666,*,iostat=myiostat) ibuffer,buffer,pradius(i),ptemp(i),prho(i),pvel(i),pye(i)
        enddo
        close(666)
     end if
     call MPI_Bcast(i, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(pradius, Num_Max_Radial, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(prho, Num_Max_Radial, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(ptemp, Num_Max_Radial, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(pye, Num_Max_Radial, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(pvel, Num_Max_Radial, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     n = i-1
     pnszones = n
     ! radius is in cm
     ! density is in g/ccm
     ! temperature is in K
     ! velocity is in cm/s
     pradius = pradius * length_gf
     ptemp = ptemp 
     prho = prho * rho_gf
     pvel = pvel / clite

  else
     call CCTK_WARN(0,"This profile type is not supported!")

  endif
  
  ! Smooth things out by averaging over nearest neighbors
  if (nz_vel_smooth>0) then
    do i=1,pnszones
      pvel_temp(i) = pvel(i)
      pvel(i) = 0.d0
    enddo   
    
    nl = -FLOOR(nz_vel_smooth/2.d0) 
    nu = CEILING(nz_vel_smooth/2.d0)-1
    do i=1,pnszones
      do j=nl,nu
        idx = i+j
        if (idx<1) idx = 1
        if (idx>pnszones) idx = pnszones
        pvel(i) = pvel(i) + pvel_temp(idx)/10.d0
      enddo
    enddo  
  endif


#if 0
  write(6,*) n
  do i=1,n
     write(6,"(i5,1P10E15.6)") i,pradius(i),prho(i),ptemp(i),pye(i),pvel(i)
  enddo
  call CCTK_WARN(0,"guender")
#endif


end subroutine PNSMapper_ReadProfile

