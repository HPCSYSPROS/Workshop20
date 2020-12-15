#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"


subroutine EOS_Omni_BarotropicReadTable(CCTK_ARGUMENTS)

  use EOS_Omni_Module
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  character(len=512) :: warnline
  character(len=256) :: eosfilename
  CCTK_INT :: slength, count
  LOGICAL :: tablethere
  CCTK_INT :: lnrho,lnye,lntemp
  CCTK_REAL :: temp_press, buffer1
  integer :: i, iostatus

  ! convert file name to fortran string
  call CCTK_FortranString ( slength, barotropiceos_table_name, &
       eosfilename )

  ! read and parse EOS table file
  inquire(file=trim(adjustl(eosfilename)),exist=tablethere)

  if(.not.tablethere) then
     write(warnline,"(A10,A,A15)") "EOS file ", trim(adjustl(eosfilename)), " does not exist!" 
     call CCTK_ERROR(warnline)
     STOP
  endif

  ! read header
  open(unit=667,file=trim(adjustl(eosfilename)),form="formatted")
  read(667,*) buffer1

  barotropiceos_energyshift = buffer1*eps_gf

  if(barotropiceos_use_thermal_gamma_law.ne.1) then
     barotropiceos_thfac = 0.0d0
  endif

  ! figure out how many entries we have
  count = 0
  do
     read(667,*,iostat=iostatus) buffer1
     if(iostatus .ne. 0) exit
     count = count + 1
  enddo
  close(667)
  barotropiceos_nrho = count

  allocate(barotropiceos_logrho(barotropiceos_nrho))
  allocate(barotropiceos_logpress(barotropiceos_nrho))
  allocate(barotropiceos_logeps(barotropiceos_nrho))
  allocate(barotropiceos_temp(barotropiceos_nrho))
  allocate(barotropiceos_ye(barotropiceos_nrho))

  ! now read everything
  open(unit=667,file=trim(adjustl(eosfilename)),form="formatted")
  read(667,*) buffer1
  do i=1,barotropiceos_nrho
     read(667,*) barotropiceos_logrho(i), &
          barotropiceos_logpress(i),barotropiceos_logeps(i),barotropiceos_ye(i),&
          barotropiceos_temp(i)
!     read(667,"(E23.14,E23.14,E23.14)") barotropiceos_logrho(i), &
!          barotropiceos_logpress(i),barotropiceos_logeps(i),barotropiceos_temp(i),&
!          barotropiceos_ye(i)
  enddo
  close(667)


  barotropiceos_logrho(:) = log10(10.0d0**barotropiceos_logrho(:) * rho_gf)

  barotropiceos_rhomax = 10.0**barotropiceos_logrho(barotropiceos_nrho)
  barotropiceos_rhomin = 10.0**barotropiceos_logrho(1)
  
  ! set up drhos
  barotropiceos_dlrho = barotropiceos_logrho(2)-barotropiceos_logrho(1)
  barotropiceos_dlrhoi = 1.0d0/barotropiceos_dlrho

  ! set up low_kappa to be able to extrapolate to very low densities
  ! we are using gamma=2
  barotropiceos_low_gamma = 2.0d0
  temp_press = 10.0**barotropiceos_logpress(1)
  barotropiceos_low_kappa = temp_press / &
       ((10.0d0**barotropiceos_logrho(1))**barotropiceos_low_gamma)


end subroutine EOS_Omni_BarotropicReadTable
