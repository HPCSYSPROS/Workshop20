#include "cctk.h"
  
#ifndef CCTK_BLAS_INT8
#  error "Internal error: CCTK_BLAS_INT8 not defined"
#endif
#if CCTK_BLAS_INT8
#  define IK 8
#else
#  define IK 4
#endif

module lapack
  implicit none

  integer, parameter :: lapack_integer_kind = IK
  
  interface geev
     SUBROUTINE SGEEV( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR, &
          LDVR, WORK, LWORK, INFO )
       IMPLICIT NONE
       CHARACTER          JOBVL, JOBVR
       INTEGER*IK         INFO, LDA, LDVL, LDVR, LWORK, N
       REAL               A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ), &
            WI( * ), WORK( * ), WR( * )
     END SUBROUTINE SGEEV
     
     SUBROUTINE DGEEV( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR, &
          LDVR, WORK, LWORK, INFO )
       IMPLICIT NONE
       CHARACTER          JOBVL, JOBVR
       INTEGER*IK         INFO, LDA, LDVL, LDVR, LWORK, N
       DOUBLE PRECISION   A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ), &
            WI( * ), WORK( * ), WR( * )
     END SUBROUTINE DGEEV
     
     SUBROUTINE CGEEV( JOBVL, JOBVR, N, A, LDA, W, VL, LDVL, VR, LDVR, &
          WORK, LWORK, RWORK, INFO )
       IMPLICIT NONE
       CHARACTER          JOBVL, JOBVR
       INTEGER*IK         INFO, LDA, LDVL, LDVR, LWORK, N
       REAL               RWORK( * )
       COMPLEX            A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ), &
            W( * ), WORK( * )
     END SUBROUTINE CGEEV
     
     SUBROUTINE ZGEEV( JOBVL, JOBVR, N, A, LDA, W, VL, LDVL, VR, LDVR, &
          WORK, LWORK, RWORK, INFO )
       IMPLICIT NONE
       CHARACTER          JOBVL, JOBVR
       INTEGER*IK         INFO, LDA, LDVL, LDVR, LWORK, N
       DOUBLE PRECISION   RWORK( * )
       COMPLEX*16         A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ), &
            W( * ), WORK( * )
     END SUBROUTINE ZGEEV
  end interface geev
  
  interface gesv
     SUBROUTINE SGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
       IMPLICIT NONE
       INTEGER*IK         INFO, LDA, LDB, N, NRHS
       INTEGER*IK         IPIV( * )
       REAL               A( LDA, * ), B( LDB, * )
     END SUBROUTINE SGESV
     
     SUBROUTINE DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
       IMPLICIT NONE
       INTEGER*IK         INFO, LDA, LDB, N, NRHS
       INTEGER*IK         IPIV( * )
       DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
     END SUBROUTINE DGESV
     
     SUBROUTINE CGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
       IMPLICIT NONE
       INTEGER*IK         INFO, LDA, LDB, N, NRHS
       INTEGER*IK         IPIV( * )
       COMPLEX            A( LDA, * ), B( LDB, * )
     END SUBROUTINE CGESV
     
      SUBROUTINE ZGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
        IMPLICIT NONE
        INTEGER*IK         INFO, LDA, LDB, N, NRHS
        INTEGER*IK         IPIV( * )
        COMPLEX*16         A( LDA, * ), B( LDB, * )
      END SUBROUTINE ZGESV
   end interface gesv
   
   interface getrf
      SUBROUTINE SGETRF( M, N, A, LDA, IPIV, INFO )
        IMPLICIT NONE
        INTEGER*IK         INFO, LDA, M, N
        INTEGER*IK         IPIV( * )
        REAL               A( LDA, * )
      END SUBROUTINE SGETRF
      
      SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
        IMPLICIT NONE
        INTEGER*IK         INFO, LDA, M, N
        INTEGER*IK         IPIV( * )
        DOUBLE PRECISION   A( LDA, * )
      END SUBROUTINE DGETRF
      
      SUBROUTINE CGETRF( M, N, A, LDA, IPIV, INFO )
        IMPLICIT NONE
        INTEGER*IK         INFO, LDA, M, N
        INTEGER*IK         IPIV( * )
        COMPLEX            A( LDA, * )
      END SUBROUTINE CGETRF
      
      SUBROUTINE ZGETRF( M, N, A, LDA, IPIV, INFO )
        IMPLICIT NONE
        INTEGER*IK         INFO, LDA, M, N
        INTEGER*IK         IPIV( * )
        COMPLEX*16         A( LDA, * )
      END SUBROUTINE ZGETRF
   end interface getrf
   
   interface posv
      SUBROUTINE SPOSV( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
        IMPLICIT NONE
        CHARACTER          UPLO
        INTEGER*IK         INFO, LDA, LDB, N, NRHS
        REAL               A( LDA, * ), B( LDB, * )
      END SUBROUTINE SPOSV
      
      SUBROUTINE DPOSV( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
        IMPLICIT NONE
        CHARACTER          UPLO
        INTEGER*IK         INFO, LDA, LDB, N, NRHS
        DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
      END SUBROUTINE DPOSV
      
      SUBROUTINE CPOSV( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
        IMPLICIT NONE
        CHARACTER          UPLO
        INTEGER*IK         INFO, LDA, LDB, N, NRHS
        COMPLEX            A( LDA, * ), B( LDB, * )
      END SUBROUTINE CPOSV
      
      SUBROUTINE ZPOSV( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
        IMPLICIT NONE
        CHARACTER          UPLO
        INTEGER*IK         INFO, LDA, LDB, N, NRHS
        COMPLEX*16         A( LDA, * ), B( LDB, * )
      END SUBROUTINE ZPOSV
   end interface posv
   
   interface sysv
      SUBROUTINE SSYSV( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, &
           LWORK, INFO )
        IMPLICIT NONE
        CHARACTER          UPLO
        INTEGER*IK         INFO, LDA, LDB, LWORK, N, NRHS
        INTEGER*IK         IPIV( * )
        REAL               A( LDA, * ), B( LDB, * ), WORK( * )
      END SUBROUTINE SSYSV
      
      SUBROUTINE DSYSV( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, &
           LWORK, INFO )
        IMPLICIT NONE
        CHARACTER          UPLO
        INTEGER*IK         INFO, LDA, LDB, LWORK, N, NRHS
        INTEGER*IK         IPIV( * )
        DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), WORK( * )
      END SUBROUTINE DSYSV
      
      SUBROUTINE CSYSV( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, &
           LWORK, INFO )
        IMPLICIT NONE
        CHARACTER          UPLO
        INTEGER*IK         INFO, LDA, LDB, LWORK, N, NRHS
        INTEGER*IK         IPIV( * )
        COMPLEX            A( LDA, * ), B( LDB, * ), WORK( * )
      END SUBROUTINE CSYSV
      
      SUBROUTINE ZSYSV( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, &
           LWORK, INFO )
        IMPLICIT NONE
        CHARACTER          UPLO
        INTEGER*IK         INFO, LDA, LDB, LWORK, N, NRHS
        INTEGER*IK         IPIV( * )
        COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( * )
      END SUBROUTINE ZSYSV
   end interface sysv
   
   interface sytrf
      SUBROUTINE SSYTRF( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
        IMPLICIT NONE
        CHARACTER          UPLO
        INTEGER*IK         INFO, LDA, LWORK, N
        INTEGER*IK         IPIV( * )
        REAL               A( LDA, * ), WORK( * )
      END SUBROUTINE SSYTRF
      
      SUBROUTINE DSYTRF( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
        IMPLICIT NONE
        CHARACTER          UPLO
        INTEGER*IK         INFO, LDA, LWORK, N
        INTEGER*IK         IPIV( * )
        DOUBLE PRECISION   A( LDA, * ), WORK( * )
      END SUBROUTINE DSYTRF
      
      SUBROUTINE CSYTRF( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
        IMPLICIT NONE
        CHARACTER          UPLO
        INTEGER*IK         INFO, LDA, LWORK, N
        INTEGER*IK         IPIV( * )
        COMPLEX            A( LDA, * ), WORK( * )
      END SUBROUTINE CSYTRF
      
      SUBROUTINE ZSYTRF( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
        IMPLICIT NONE
        CHARACTER          UPLO
        INTEGER*IK         INFO, LDA, LWORK, N
        INTEGER*IK         IPIV( * )
        COMPLEX*16         A( LDA, * ), WORK( * )
      END SUBROUTINE ZSYTRF
   end interface sytrf
 end module lapack
