#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"


subroutine EOS_Omni_ReadColdTable(CCTK_ARGUMENTS)

  use EOS_Omni_Module
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  character(len=512) :: warnline
  character(len=256) :: eosfilename
  CCTK_INT :: slength
  LOGICAL :: tablethere
  character(len=256) :: buffer1,buffer2,buffer3
  CCTK_INT :: lnrho,lnye,lntemp
  CCTK_REAL :: temp_press
 integer :: i

  ! convert file name to fortran string
  call CCTK_FortranString ( slength, coldeos_table_name, &
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
  read(667,"(A)") buffer1
  ! check if okay
  if(.not. (trim(adjustl(buffer1)).eq."EoSType = Tabulated")) then
     call CCTK_ERROR("Wrong EOS table format -- check header")
     STOP
  endif
  ! now read number of rho entries
  read(667,"(A7,i7,A6,i6,A5,i6)") buffer1,lnrho,buffer2,lnye,buffer3,lntemp
  if(lnye.ne.1.or.lntemp.ne.1) then
     call CCTK_ERROR("Wrong EOS table format -- cannot handle Nye!=1 or NT!=`")
     STOP
  endif
  coldeos_nrho = lnrho


  ! allocate vars
  allocate(coldeos_logrho(coldeos_nrho))
  allocate(coldeos_eps(coldeos_nrho))
  allocate(coldeos_gamma(coldeos_nrho))
  allocate(coldeos_cs2(coldeos_nrho))

  ! read rho min and max
  read(667,"(A16,A16)") buffer1,buffer2
  read(buffer1(10:),*) coldeos_rhomin
  read(buffer2(10:),*) coldeos_rhomax

  ! read heat capacity? Not sure what exactly that is used for, just ignore it
  ! for now
  read(667,*) buffer1
  
  ! read gamma thermal
  read(667,"(A30)") buffer1
  read(buffer1(10:),*) coldeos_gammath

  if(coldeos_use_thermal_gamma_law.ne.1) then
     coldeos_thfac = 0.0d0
  endif

  ! read kappa 
  read(667,"(A30)") buffer1
  read(buffer1(9:),*) coldeos_kappa


  ! make sure density is in log spacing
  read(667,"(A30)") buffer1
  if(.not.(trim(adjustl(buffer1)).eq."RhoSpacing = Log")) then
     call CCTK_ERROR("Density spacing not log? Check table format!")
     STOP
  endif

  ! now read everything
  do i=1,coldeos_nrho
     read(667,"(E24.14,E24.14,E24.14)") coldeos_eps(i),coldeos_gamma(i),coldeos_cs2(i)
     coldeos_cs2(i) = coldeos_cs2(i)**2
  enddo
  close(667)

  ! set up rho
  coldeos_dlrho = (log10(coldeos_rhomax) - log10(coldeos_rhomin)) / (coldeos_nrho-1)
  coldeos_dlrhoi = 1.0d0/coldeos_dlrho
  do i=1,coldeos_nrho
     coldeos_logrho(i) = log10(coldeos_rhomin) + (i-1)*coldeos_dlrho
  enddo

  ! set up low_kappa to be able to extrapolate to very low densities
  ! we are using gamma=2
  coldeos_low_gamma = 2.0d0
  temp_press = (10.0d0**coldeos_logrho(1))**coldeos_gamma(1)
  coldeos_low_kappa = temp_press / &
       ((10.0d0**coldeos_logrho(1))**coldeos_low_gamma)


#if 0
  ! debug output
  do i=1,coldeos_nrho
     write(6,"(i5,1P10E15.6)") i, coldeos_logrho(i), &
          coldeos_eps(i),coldeos_gamma(i),sqrt(coldeos_cs2(i))
  enddo
#endif


end subroutine EOS_Omni_ReadColdTable
