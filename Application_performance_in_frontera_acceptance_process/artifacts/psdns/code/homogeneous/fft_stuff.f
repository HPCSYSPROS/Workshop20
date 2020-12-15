#ifdef FFTW
      include "fftw3.f"

!     implicit none

      integer*8 plan1,plan2,plan3 
! next line added PK Yeung 2/4/07
	integer*8 plan4
      integer, parameter :: NULL = 0
!	integer, parameter :: fftw_flag = FFTW_ESTIMATE
       integer, parameter :: fftw_flag = FFTW_MEASURE
!     integer, parameter :: fftw_flag = FFTW_PATIENT
!     integer, parameter :: fftw_flag = FFTW_EXHAUSTIVE


#ifdef DOUBLE_PREC
#define fftw_plan_r2c dfftw_plan_many_dft_r2c
#define fftw_plan_c2r dfftw_plan_many_dft_c2r
#define fftw_plan dfftw_plan_many_dft
#define fftw_execute_r2c dfftw_execute_dft_r2c
#define fftw_execute_c2r dfftw_execute_dft_c2r
#define fftw_execute dfftw_execute_dft
#define fftw_destroy dfftw_destroy_plan
#else
#define fftw_plan_r2c sfftw_plan_many_dft_r2c
#define fftw_plan_c2r sfftw_plan_many_dft_c2r
#define fftw_plan sfftw_plan_many_dft
#define fftw_execute_r2c sfftw_execute_dft_r2c
#define fftw_execute_c2r sfftw_execute_dft_c2r
#define fftw_execute sfftw_execute_dft
#define fftw_destroy sfftw_destroy_plan
#endif

#else

#ifdef DOUBLE_PREC
#define rcft drcft
#define crft dcrft
#define cft dcft
#else
#define rcft srcft
#define crft scrft
#define cft scft
#endif

!  vec----------------------------------------------------------------  
! common blocks when using essllibs vector fft subroutines              
!                                                                       

!         real(8),allocatable,save:: caux1(:),caux2(:)
!         real(8),allocatable,save:: raux1(:),raux2(:),raux3(:)
    

!         common/vecfft/caux1(cnaux1),caux2(cnaux2),                     
!     1              raux1(rnaux1),raux2(rnaux2),raux3(rnaux3)           

!      use fft_ar
!      implicit none

#endif
