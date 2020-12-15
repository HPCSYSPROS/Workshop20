#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine qlm_outputvtk(CCTK_ARGUMENTS,nhor,file_name,unit_nr)
  use cctk
  use constants
  use qlm_variables
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  
    integer, intent(in)          :: nhor
    character(len=*), intent(in) :: file_name
    integer, intent(in)          :: unit_nr 

    integer   :: i, j, nth, nph, ngth, ngph
    CCTK_REAL :: xx, yy, zz

    ngth = qlm_nghoststheta(nhor) 
    ngph = qlm_nghostsphi(nhor)
    nth = qlm_ntheta(nhor)-2*ngth
    nph = qlm_nphi(nhor)-2*ngph

    open (unit=unit_nr, file=file_name,action='write')
    
    write(unit_nr,'(A)') '# vtk DataFile Version 2.0' 
    write(unit_nr,'(A)') 'Horizon data' 
    write(unit_nr,'(A)') 'ASCII' 
    write(unit_nr,'(A)') 'DATASET POLYDATA'
    write(unit_nr,'(A,X,I5,X,A)') 'POINTS', nth*nph, 'float'

    do j=1, nth
       do i=1, nph
          xx = qlm_x(j+ngth,i+ngph,nhor) 
          yy = qlm_y(j+ngth,i+ngph,nhor) 
          zz = qlm_z(j+ngth,i+ngph,nhor) 
          write(unit_nr,*) xx, yy, zz 
       end do
    end do

    write(unit_nr,'(A)') ''
    write(unit_nr,'(A,X,I5,X,I10)') 'POLYGONS', nph*(nth-1), 5*nph*(nth-1)
    
    do j=0, nth-2
       do i=0, nph-2
          write(unit_nr,'(I1,4(I6))') 4, j*nph+i, (j+1)*nph+i, (j+1)*nph+i+1, j*nph+i+1
       end do
       write(unit_nr,'(I1,4(I6))') 4, j*nph+nph-1, (j+1)*nph+nph-1, (j+1)*nph, j*nph
    end do

    write(unit_nr,'(A)') ''
    write(unit_nr,'(A,X,I5)') 'POINT_DATA', nth*nph

    call qlm_writescalar(CCTK_PASS_FTOF,nth,nph,nhor,'shape',unit_nr)
    call qlm_writescalar(CCTK_PASS_FTOF,nth,nph,nhor,'l0',unit_nr)
    call qlm_writescalar(CCTK_PASS_FTOF,nth,nph,nhor,'l1',unit_nr)
    call qlm_writescalar(CCTK_PASS_FTOF,nth,nph,nhor,'l2',unit_nr)
    call qlm_writescalar(CCTK_PASS_FTOF,nth,nph,nhor,'l3',unit_nr)
    call qlm_writescalar(CCTK_PASS_FTOF,nth,nph,nhor,'n0',unit_nr)
    call qlm_writescalar(CCTK_PASS_FTOF,nth,nph,nhor,'n1',unit_nr)
    call qlm_writescalar(CCTK_PASS_FTOF,nth,nph,nhor,'n2',unit_nr)
    call qlm_writescalar(CCTK_PASS_FTOF,nth,nph,nhor,'n3',unit_nr)
    call qlm_writescalar(CCTK_PASS_FTOF,nth,nph,nhor,'rem0',unit_nr)
    call qlm_writescalar(CCTK_PASS_FTOF,nth,nph,nhor,'imm0',unit_nr)
    call qlm_writescalar(CCTK_PASS_FTOF,nth,nph,nhor,'rem1',unit_nr)
    call qlm_writescalar(CCTK_PASS_FTOF,nth,nph,nhor,'imm1',unit_nr)
    call qlm_writescalar(CCTK_PASS_FTOF,nth,nph,nhor,'rem2',unit_nr)
    call qlm_writescalar(CCTK_PASS_FTOF,nth,nph,nhor,'imm2',unit_nr)
    call qlm_writescalar(CCTK_PASS_FTOF,nth,nph,nhor,'rem3',unit_nr)
    call qlm_writescalar(CCTK_PASS_FTOF,nth,nph,nhor,'imm3',unit_nr)
    call qlm_writescalar(CCTK_PASS_FTOF,nth,nph,nhor,'renpkappa',unit_nr)
    call qlm_writescalar(CCTK_PASS_FTOF,nth,nph,nhor,'imnpkappa',unit_nr)
    call qlm_writescalar(CCTK_PASS_FTOF,nth,nph,nhor,'renptau',unit_nr)
    call qlm_writescalar(CCTK_PASS_FTOF,nth,nph,nhor,'imnptau',unit_nr)
    call qlm_writescalar(CCTK_PASS_FTOF,nth,nph,nhor,'renpsigma',unit_nr)
    call qlm_writescalar(CCTK_PASS_FTOF,nth,nph,nhor,'imnpsigma',unit_nr)
    call qlm_writescalar(CCTK_PASS_FTOF,nth,nph,nhor,'renprho',unit_nr)
    call qlm_writescalar(CCTK_PASS_FTOF,nth,nph,nhor,'imnprho',unit_nr)
    call qlm_writescalar(CCTK_PASS_FTOF,nth,nph,nhor,'renpepsilon',unit_nr)
    call qlm_writescalar(CCTK_PASS_FTOF,nth,nph,nhor,'imnpepsilon',unit_nr)
    call qlm_writescalar(CCTK_PASS_FTOF,nth,nph,nhor,'renpgamma',unit_nr)
    call qlm_writescalar(CCTK_PASS_FTOF,nth,nph,nhor,'imnpgamma',unit_nr)
    call qlm_writescalar(CCTK_PASS_FTOF,nth,nph,nhor,'renpbeta',unit_nr)
    call qlm_writescalar(CCTK_PASS_FTOF,nth,nph,nhor,'imnpbeta',unit_nr)
    call qlm_writescalar(CCTK_PASS_FTOF,nth,nph,nhor,'renpalpha',unit_nr)
    call qlm_writescalar(CCTK_PASS_FTOF,nth,nph,nhor,'imnpalpha',unit_nr)
    call qlm_writescalar(CCTK_PASS_FTOF,nth,nph,nhor,'renppi',unit_nr)
    call qlm_writescalar(CCTK_PASS_FTOF,nth,nph,nhor,'imnppi',unit_nr)
    call qlm_writescalar(CCTK_PASS_FTOF,nth,nph,nhor,'renpnu',unit_nr)
    call qlm_writescalar(CCTK_PASS_FTOF,nth,nph,nhor,'imnpnu',unit_nr)
    call qlm_writescalar(CCTK_PASS_FTOF,nth,nph,nhor,'renpmu',unit_nr)
    call qlm_writescalar(CCTK_PASS_FTOF,nth,nph,nhor,'imnpmu',unit_nr)
    call qlm_writescalar(CCTK_PASS_FTOF,nth,nph,nhor,'renplambda',unit_nr)
    call qlm_writescalar(CCTK_PASS_FTOF,nth,nph,nhor,'imnplambda',unit_nr)
    call qlm_writescalar(CCTK_PASS_FTOF,nth,nph,nhor,'repsi0',unit_nr)
    call qlm_writescalar(CCTK_PASS_FTOF,nth,nph,nhor,'impsi0',unit_nr)
    call qlm_writescalar(CCTK_PASS_FTOF,nth,nph,nhor,'repsi1',unit_nr)
    call qlm_writescalar(CCTK_PASS_FTOF,nth,nph,nhor,'impsi1',unit_nr)
    call qlm_writescalar(CCTK_PASS_FTOF,nth,nph,nhor,'repsi2',unit_nr)
    call qlm_writescalar(CCTK_PASS_FTOF,nth,nph,nhor,'impsi2',unit_nr)
    call qlm_writescalar(CCTK_PASS_FTOF,nth,nph,nhor,'repsi3',unit_nr)
    call qlm_writescalar(CCTK_PASS_FTOF,nth,nph,nhor,'impsi3',unit_nr)
    call qlm_writescalar(CCTK_PASS_FTOF,nth,nph,nhor,'repsi4',unit_nr)
    call qlm_writescalar(CCTK_PASS_FTOF,nth,nph,nhor,'impsi4',unit_nr)
    call qlm_writescalar(CCTK_PASS_FTOF,nth,nph,nhor,'xit',unit_nr)
    call qlm_writescalar(CCTK_PASS_FTOF,nth,nph,nhor,'xip',unit_nr)
    call qlm_writescalar(CCTK_PASS_FTOF,nth,nph,nhor,'chi',unit_nr)

    close(1)

end subroutine qlm_outputvtk

subroutine qlm_writescalar(CCTK_ARGUMENTS,nth,nph,nhor,array_name,unit_nr)
  use cctk
  use constants
  use qlm_variables
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  
    integer, intent(in)          :: nth, nph, nhor
    character(len=*), intent(in) :: array_name
    integer, intent(in)          :: unit_nr 
    CCTK_REAL                    :: array(1:nth,1:nph)

    integer   :: i, j

    select case (array_name)
        case('shape')
                array(1:nth,1:nph) = qlm_shape(qlm_nghoststheta(nhor):nth+qlm_nghoststheta(nhor),&
                                                qlm_nghostsphi(nhor):nph+qlm_nghostsphi(nhor),nhor)
        case('l0')
                array(1:nth,1:nph) = qlm_l0(qlm_nghoststheta(nhor):nth+qlm_nghoststheta(nhor),&
                                                qlm_nghostsphi(nhor):nph+qlm_nghostsphi(nhor),nhor)
        case('l1')
                array(1:nth,1:nph) = qlm_l1(qlm_nghoststheta(nhor):nth+qlm_nghoststheta(nhor),&
                                                qlm_nghostsphi(nhor):nph+qlm_nghostsphi(nhor),nhor)
        case('l2')
                array(1:nth,1:nph) = qlm_l2(qlm_nghoststheta(nhor):nth+qlm_nghoststheta(nhor),&
                                                qlm_nghostsphi(nhor):nph+qlm_nghostsphi(nhor),nhor)
        case('l3')
                array(1:nth,1:nph) = qlm_l3(qlm_nghoststheta(nhor):nth+qlm_nghoststheta(nhor),&
                                                qlm_nghostsphi(nhor):nph+qlm_nghostsphi(nhor),nhor)
        case('n0')
                array(1:nth,1:nph) = qlm_n0(qlm_nghoststheta(nhor):nth+qlm_nghoststheta(nhor),&
                                                qlm_nghostsphi(nhor):nph+qlm_nghostsphi(nhor),nhor)
        case('n1')
                array(1:nth,1:nph) = qlm_n1(qlm_nghoststheta(nhor):nth+qlm_nghoststheta(nhor),&
                                                qlm_nghostsphi(nhor):nph+qlm_nghostsphi(nhor),nhor)
        case('n2')
                array(1:nth,1:nph) = qlm_n2(qlm_nghoststheta(nhor):nth+qlm_nghoststheta(nhor),&
                                                qlm_nghostsphi(nhor):nph+qlm_nghostsphi(nhor),nhor)
        case('n3')
                array(1:nth,1:nph) = qlm_n3(qlm_nghoststheta(nhor):nth+qlm_nghoststheta(nhor),&
                                                qlm_nghostsphi(nhor):nph+qlm_nghostsphi(nhor),nhor)
        case('rem0')
                array(1:nth,1:nph) = real(qlm_m0(qlm_nghoststheta(nhor):nth+qlm_nghoststheta(nhor),&
                                                qlm_nghostsphi(nhor):nph+qlm_nghostsphi(nhor),nhor))
        case('imm0')
                array(1:nth,1:nph) = aimag(qlm_m0(qlm_nghoststheta(nhor):nth+qlm_nghoststheta(nhor),&
                                                qlm_nghostsphi(nhor):nph+qlm_nghostsphi(nhor),nhor))
        case('rem1')
                array(1:nth,1:nph) = real(qlm_m1(qlm_nghoststheta(nhor):nth+qlm_nghoststheta(nhor),&
                                                qlm_nghostsphi(nhor):nph+qlm_nghostsphi(nhor),nhor))
        case('imm1')
                array(1:nth,1:nph) = aimag(qlm_m1(qlm_nghoststheta(nhor):nth+qlm_nghoststheta(nhor),&
                                                qlm_nghostsphi(nhor):nph+qlm_nghostsphi(nhor),nhor))
        case('rem2')
                array(1:nth,1:nph) = real(qlm_m2(qlm_nghoststheta(nhor):nth+qlm_nghoststheta(nhor),&
                                                qlm_nghostsphi(nhor):nph+qlm_nghostsphi(nhor),nhor))
        case('imm2')
                array(1:nth,1:nph) = aimag(qlm_m2(qlm_nghoststheta(nhor):nth+qlm_nghoststheta(nhor),&
                                                qlm_nghostsphi(nhor):nph+qlm_nghostsphi(nhor),nhor))
        case('rem3')
                array(1:nth,1:nph) = real(qlm_m3(qlm_nghoststheta(nhor):nth+qlm_nghoststheta(nhor),&
                                                qlm_nghostsphi(nhor):nph+qlm_nghostsphi(nhor),nhor))
        case('imm3')
                array(1:nth,1:nph) = aimag(qlm_m3(qlm_nghoststheta(nhor):nth+qlm_nghoststheta(nhor),&
                                                qlm_nghostsphi(nhor):nph+qlm_nghostsphi(nhor),nhor))
        case('renpkappa')
                array(1:nth,1:nph) = real(qlm_npkappa(qlm_nghoststheta(nhor):nth+qlm_nghoststheta(nhor),&
                                                qlm_nghostsphi(nhor):nph+qlm_nghostsphi(nhor),nhor))
        case('imnpkappa')
                array(1:nth,1:nph) = aimag(qlm_npkappa(qlm_nghoststheta(nhor):nth+qlm_nghoststheta(nhor),&
                                                qlm_nghostsphi(nhor):nph+qlm_nghostsphi(nhor),nhor))
        case('renptau')
                array(1:nth,1:nph) = real(qlm_nptau(qlm_nghoststheta(nhor):nth+qlm_nghoststheta(nhor),&
                                                qlm_nghostsphi(nhor):nph+qlm_nghostsphi(nhor),nhor))
        case('imnptau')
                array(1:nth,1:nph) = aimag(qlm_nptau(qlm_nghoststheta(nhor):nth+qlm_nghoststheta(nhor),&
                                                qlm_nghostsphi(nhor):nph+qlm_nghostsphi(nhor),nhor))
        case('renpsigma')
                array(1:nth,1:nph) = real(qlm_npsigma(qlm_nghoststheta(nhor):nth+qlm_nghoststheta(nhor),&
                                                qlm_nghostsphi(nhor):nph+qlm_nghostsphi(nhor),nhor))
        case('imnpsigma')
                array(1:nth,1:nph) = aimag(qlm_npsigma(qlm_nghoststheta(nhor):nth+qlm_nghoststheta(nhor),&
                                                qlm_nghostsphi(nhor):nph+qlm_nghostsphi(nhor),nhor))
        case('renprho')
                array(1:nth,1:nph) = real(qlm_nprho(qlm_nghoststheta(nhor):nth+qlm_nghoststheta(nhor),&
                                                qlm_nghostsphi(nhor):nph+qlm_nghostsphi(nhor),nhor))
        case('imnprho')
                array(1:nth,1:nph) = aimag(qlm_nprho(qlm_nghoststheta(nhor):nth+qlm_nghoststheta(nhor),&
                                                qlm_nghostsphi(nhor):nph+qlm_nghostsphi(nhor),nhor))
        case('renpepsilon')
                array(1:nth,1:nph) = real(qlm_npepsilon(qlm_nghoststheta(nhor):nth+qlm_nghoststheta(nhor),&
                                                qlm_nghostsphi(nhor):nph+qlm_nghostsphi(nhor),nhor))
        case('imnpepsilon')
                array(1:nth,1:nph) = aimag(qlm_npepsilon(qlm_nghoststheta(nhor):nth+qlm_nghoststheta(nhor),&
                                                qlm_nghostsphi(nhor):nph+qlm_nghostsphi(nhor),nhor))
        case('renpgamma')
                array(1:nth,1:nph) = real(qlm_npgamma(qlm_nghoststheta(nhor):nth+qlm_nghoststheta(nhor),&
                                                qlm_nghostsphi(nhor):nph+qlm_nghostsphi(nhor),nhor))
        case('imnpgamma')
                array(1:nth,1:nph) = aimag(qlm_npgamma(qlm_nghoststheta(nhor):nth+qlm_nghoststheta(nhor),&
                                                qlm_nghostsphi(nhor):nph+qlm_nghostsphi(nhor),nhor))
        case('renpbeta')
                array(1:nth,1:nph) = real(qlm_npbeta(qlm_nghoststheta(nhor):nth+qlm_nghoststheta(nhor),&
                                                qlm_nghostsphi(nhor):nph+qlm_nghostsphi(nhor),nhor))
        case('imnpbeta')
                array(1:nth,1:nph) = aimag(qlm_npbeta(qlm_nghoststheta(nhor):nth+qlm_nghoststheta(nhor),&
                                                qlm_nghostsphi(nhor):nph+qlm_nghostsphi(nhor),nhor))
        case('renpalpha')
                array(1:nth,1:nph) = real(qlm_npalpha(qlm_nghoststheta(nhor):nth+qlm_nghoststheta(nhor),&
                                                qlm_nghostsphi(nhor):nph+qlm_nghostsphi(nhor),nhor))
        case('imnpalpha')
                array(1:nth,1:nph) = aimag(qlm_npalpha(qlm_nghoststheta(nhor):nth+qlm_nghoststheta(nhor),&
                                                qlm_nghostsphi(nhor):nph+qlm_nghostsphi(nhor),nhor))
        case('renppi')
                array(1:nth,1:nph) = real(qlm_nppi(qlm_nghoststheta(nhor):nth+qlm_nghoststheta(nhor),&
                                                qlm_nghostsphi(nhor):nph+qlm_nghostsphi(nhor),nhor))
        case('imnppi')
                array(1:nth,1:nph) = aimag(qlm_nppi(qlm_nghoststheta(nhor):nth+qlm_nghoststheta(nhor),&
                                                qlm_nghostsphi(nhor):nph+qlm_nghostsphi(nhor),nhor))
        case('renpnu')
                array(1:nth,1:nph) = real(qlm_npnu(qlm_nghoststheta(nhor):nth+qlm_nghoststheta(nhor),&
                                                qlm_nghostsphi(nhor):nph+qlm_nghostsphi(nhor),nhor))
        case('imnpnu')
                array(1:nth,1:nph) = aimag(qlm_npnu(qlm_nghoststheta(nhor):nth+qlm_nghoststheta(nhor),&
                                                qlm_nghostsphi(nhor):nph+qlm_nghostsphi(nhor),nhor))
        case('renpmu')
                array(1:nth,1:nph) = real(qlm_npmu(qlm_nghoststheta(nhor):nth+qlm_nghoststheta(nhor),&
                                                qlm_nghostsphi(nhor):nph+qlm_nghostsphi(nhor),nhor))
        case('imnpmu')
                array(1:nth,1:nph) = aimag(qlm_npmu(qlm_nghoststheta(nhor):nth+qlm_nghoststheta(nhor),&
                                                qlm_nghostsphi(nhor):nph+qlm_nghostsphi(nhor),nhor))
        case('renplambda')
                array(1:nth,1:nph) = real(qlm_nplambda(qlm_nghoststheta(nhor):nth+qlm_nghoststheta(nhor),&
                                                qlm_nghostsphi(nhor):nph+qlm_nghostsphi(nhor),nhor))
        case('imnplambda')
                array(1:nth,1:nph) = aimag(qlm_nplambda(qlm_nghoststheta(nhor):nth+qlm_nghoststheta(nhor),&
                                                qlm_nghostsphi(nhor):nph+qlm_nghostsphi(nhor),nhor))
        case('repsi0')
                array(1:nth,1:nph) = real(qlm_psi0(qlm_nghoststheta(nhor):nth+qlm_nghoststheta(nhor),&
                                                qlm_nghostsphi(nhor):nph+qlm_nghostsphi(nhor),nhor))
        case('impsi0')
                array(1:nth,1:nph) = aimag(qlm_psi0(qlm_nghoststheta(nhor):nth+qlm_nghoststheta(nhor),&
                                                qlm_nghostsphi(nhor):nph+qlm_nghostsphi(nhor),nhor))
        case('repsi1')
                array(1:nth,1:nph) = real(qlm_psi1(qlm_nghoststheta(nhor):nth+qlm_nghoststheta(nhor),&
                                                qlm_nghostsphi(nhor):nph+qlm_nghostsphi(nhor),nhor))
        case('impsi1')
                array(1:nth,1:nph) = aimag(qlm_psi1(qlm_nghoststheta(nhor):nth+qlm_nghoststheta(nhor),&
                                                qlm_nghostsphi(nhor):nph+qlm_nghostsphi(nhor),nhor))
        case('repsi2')
                array(1:nth,1:nph) = real(qlm_psi2(qlm_nghoststheta(nhor):nth+qlm_nghoststheta(nhor),&
                                                qlm_nghostsphi(nhor):nph+qlm_nghostsphi(nhor),nhor))
        case('impsi2')
                array(1:nth,1:nph) = aimag(qlm_psi2(qlm_nghoststheta(nhor):nth+qlm_nghoststheta(nhor),&
                                                qlm_nghostsphi(nhor):nph+qlm_nghostsphi(nhor),nhor))
        case('repsi3')
                array(1:nth,1:nph) = real(qlm_psi3(qlm_nghoststheta(nhor):nth+qlm_nghoststheta(nhor),&
                                                qlm_nghostsphi(nhor):nph+qlm_nghostsphi(nhor),nhor))
        case('impsi3')
                array(1:nth,1:nph) = aimag(qlm_psi3(qlm_nghoststheta(nhor):nth+qlm_nghoststheta(nhor),&
                                                qlm_nghostsphi(nhor):nph+qlm_nghostsphi(nhor),nhor))
        case('repsi4')
                array(1:nth,1:nph) = real(qlm_psi4(qlm_nghoststheta(nhor):nth+qlm_nghoststheta(nhor),&
                                                qlm_nghostsphi(nhor):nph+qlm_nghostsphi(nhor),nhor))
        case('impsi4')
                array(1:nth,1:nph) = aimag(qlm_psi4(qlm_nghoststheta(nhor):nth+qlm_nghoststheta(nhor),&
                                                qlm_nghostsphi(nhor):nph+qlm_nghostsphi(nhor),nhor))
        case('xit')
                array(1:nth,1:nph) = qlm_xi_t(qlm_nghoststheta(nhor):nth+qlm_nghoststheta(nhor),&
                                                qlm_nghostsphi(nhor):nph+qlm_nghostsphi(nhor),nhor)
        case('xip')
                array(1:nth,1:nph) = qlm_xi_p(qlm_nghoststheta(nhor):nth+qlm_nghoststheta(nhor),&
                                                qlm_nghostsphi(nhor):nph+qlm_nghostsphi(nhor),nhor)
        case('chi')
                array(1:nth,1:nph) = qlm_chi(qlm_nghoststheta(nhor):nth+qlm_nghoststheta(nhor),&
                                                qlm_nghostsphi(nhor):nph+qlm_nghostsphi(nhor),nhor)
    end select

    write(unit_nr,'(/A,X,A,X,A)') 'SCALARS', array_name, 'float 1'
    write(unit_nr,'(A)') 'LOOKUP_TABLE default'
    do j=1, nth
       do i=1, nph
          write(unit_nr,*) array(j,i)
       end do
    end do

end subroutine qlm_writescalar
