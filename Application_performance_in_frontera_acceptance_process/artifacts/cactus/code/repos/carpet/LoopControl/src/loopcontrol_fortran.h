/* -*-f90-*- */

#ifndef LOOPCONTROL_FORTRAN_H
#define LOOPCONTROL_FORTRAN_H

#include "cctk.h"



#define LC_COARSE_DECLARE(name, D)                              \
   && integer :: name/**/_cmin/**/D, name/**/_cmax/**/D,        \
                 name/**/_cstep/**/D, name/**/_cpos/**/D
#define LC_COARSE_OMP_PRIVATE(name, D)                          \
   && !$omp private (name/**/_cmin/**/D, name/**/_cmax/**/D,    \
                     name/**/_cstep/**/D, name/**/_cpos/**/D)
#define LC_COARSE_SETUP(name, D)                                        \
   && name/**/_control%coarse_loop%min%v(D) =                           \
         name/**/_control%coarse_thread%pos%v(D)                        \
   && name/**/_control%coarse_loop%max%v(D) =                           \
         min(name/**/_control%coarse_thread%max%v(D),                   \
             name/**/_control%coarse_loop%min%v(D) +                    \
                name/**/_control%coarse_thread%step%v(D))               \
   && name/**/_cmin/**/D = name/**/_control%coarse_loop%min%v(D)        \
   && name/**/_cmax/**/D = name/**/_control%coarse_loop%max%v(D)        \
   && name/**/_cstep/**/D = name/**/_control%coarse_loop%step%v(D)
#define LC_COARSE_LOOP(name, D)                                         \
   && do name/**/_cpos/**/D = name/**/_cmin/**/D, name/**/_cmax/**/D,   \
                              name/**/_cstep/**/D

#define LC_FINE_DECLARE(name, D)                                \
   && integer :: name/**/_fmin/**/D, name/**/_fmax/**/D,        \
                 name/**/_fstep/**/D
#define LC_FINE_OMP_PRIVATE(name, I, NI, D)                     \
   && !$omp private (name/**/_fmin/**/D, name/**/_fmax/**/D,    \
                     name/**/_fstep/**/D, I, NI)
#define LC_FINE_SETUP(name, D)                                          \
   && name/**/_control%fine_loop%min%v(D) = name/**/_cpos/**/D          \
   && name/**/_control%fine_loop%max%v(D) =                             \
         min(name/**/_control%coarse_loop%max%v(D),                     \
             name/**/_control%fine_loop%min%v(D) +                      \
                name/**/_control%coarse_loop%step%v(D))                 \
   && name/**/_fmin/**/D = name/**/_control%fine_loop%min%v(D)          \
   && name/**/_fmax/**/D = name/**/_control%fine_loop%max%v(D)          \
   && name/**/_fstep/**/D = name/**/_control%fine_loop%step%v(D)
#define LC_FINE_LOOP(name, I, NI, D)                                    \
   && do I = name/**/_fmin/**/D, name/**/_fmax/**/D,                    \
             name/**/_fstep/**/D                                        \
   &&    NI = 0                                                         \
   &&    if (name/**/_dir/**/D<0) NI = I                                \
   &&    if (name/**/_dir/**/D>0) NI = name/**/_control%overall%max%v(D)+1-I



#define LC_LOOP3STR_NORMAL_DECLARE(name)                                \
   && integer :: name/**/_dir1, name/**/_dir2, name/**/_dir3            \
   && integer :: name/**/_ash1, name/**/_ash2, name/**/_ash3            \
   && integer :: name/**/_str1                                          \
   && CCTK_POINTER, save :: name/**/_descr = 0                          \
   && type(lc_control_t) :: name/**/_control                            \
      LC_COARSE_DECLARE(name,1)                                         \
      LC_COARSE_DECLARE(name,2)                                         \
      LC_COARSE_DECLARE(name,3)                                         \
      LC_FINE_DECLARE(name,1)                                           \
      LC_FINE_DECLARE(name,2)                                           \
      LC_FINE_DECLARE(name,3)

#define LC_LOOP3STR_NORMAL_OMP_PRIVATE(name, i,j,k)     \
   && !$omp private (name/**/_control)                  \
      LC_COARSE_OMP_PRIVATE(name,1)                     \
      LC_COARSE_OMP_PRIVATE(name,2)                     \
      LC_COARSE_OMP_PRIVATE(name,3)                     \
      LC_FINE_OMP_PRIVATE(name, i, ni,1)                \
      LC_FINE_OMP_PRIVATE(name, j, nj,2)                \
      LC_FINE_OMP_PRIVATE(name, k, nk,3)



#define LC_LOOP3STR_NORMAL(name, i,j,k, ni,nj,nk,                       \
                           idir_,jdir_,kdir_,                           \
                           imin_,jmin_,kmin_,                           \
                           imax_,jmax_,kmax_,                           \
                           iash_,jash_,kash_,                           \
                           vec_imin,vec_imax, istr_)                    \
   && name/**/_dir1 = (idir_)                                           \
   && name/**/_dir2 = (jdir_)                                           \
   && name/**/_dir3 = (kdir_)                                           \
   && name/**/_ash1 = (iash_)                                           \
   && name/**/_ash2 = (jash_)                                           \
   && name/**/_ash3 = (kash_)                                           \
   && name/**/_str1 = (istr_)                                           \
                                                                        \
   && call lc_descr_init(name/**/_descr, __LINE__, __FILE__, "name")    \
   && call lc_control_init(name/**/_control, name/**/_descr,            \
                           (imin_), (jmin_), (kmin_),                   \
                           (imax_), (jmax_), (kmax_),                   \
                           name/**/_ash1, name/**/_ash2, name/**/_ash3, \
                           name/**/_str1)                               \
                                                                        \
      /* Multithreading */                                              \
   && call lc_thread_init(name/**/_control)                             \
   && do while (lc_thread_done(name/**/_control) == 0)                  \
                                                                        \
         /* Coarse loops */                                             \
         LC_COARSE_SETUP(name,3)                                        \
         LC_COARSE_SETUP(name,2)                                        \
         LC_COARSE_SETUP(name,1)                                        \
         LC_COARSE_LOOP(name,3)                                         \
         LC_COARSE_LOOP(name,2)                                         \
         LC_COARSE_LOOP(name,1)                                         \
                                                                        \
            /* Fine loops */                                            \
            LC_FINE_SETUP(name,3)                                       \
            LC_FINE_SETUP(name,2)                                       \
            LC_FINE_SETUP(name,1)                                       \
            LC_FINE_LOOP(name, k, nk,3)                                 \
            LC_FINE_LOOP(name, j, nj,2)                                 \
            LC_FINE_LOOP(name, i, ni,1)

#define LC_ENDLOOP3STR_NORMAL(name)                             \
   &&       end do                                              \
   &&       end do                                              \
   &&       end do                                              \
   &&    end do                                                 \
   &&    end do                                                 \
   &&    end do                                                 \
   &&    call lc_thread_step(name/**/_control)                  \
   && end do                                                    \
   && call lc_control_finish(name/**/_control, name/**/_descr)



#define LC_LOOP3(name, i,j,k, imin,jmin,kmin, imax,jmax,kmax, iash,jash,kash) \
   LC_LOOP3STR(name, i,j,k, imin,jmin,kmin, imax,jmax,kmax, iash,jash,kash, 1)
#define LC_ENDLOOP3(name)                       \
   LC_ENDLOOP3STR(name)



/* Replace CCTK_LOOP macros */
#if (!defined CCTK_LOOP3STR_NORMAL_DECLARE ||           \
     !defined CCTK_LOOP3STR_NORMAL_OMP_PRIVATE ||       \
     !defined CCTK_LOOP3STR_NORMAL ||                   \
     !defined CCTK_ENDLOOP3STR_NORMAL)
#  error "internal error"
#endif
#undef CCTK_LOOP3STR_NORMAL_DECLARE
#undef CCTK_LOOP3STR_NORMAL_OMP_PRIVATE
#undef CCTK_LOOP3STR_NORMAL
#undef CCTK_ENDLOOP3STR_NORMAL
#define CCTK_LOOP3STR_NORMAL_DECLARE     LC_LOOP3STR_NORMAL_DECLARE
#define CCTK_LOOP3STR_NORMAL_OMP_PRIVATE LC_LOOP3STR_NORMAL_OMP_PRIVATE
#define CCTK_LOOP3STR_NORMAL             LC_LOOP3STR_NORMAL
#define CCTK_ENDLOOP3STR_NORMAL          LC_ENDLOOP3STR_NORMAL



#endif  /* #ifndef LOOPCONTROL_FORTRAN_H */
