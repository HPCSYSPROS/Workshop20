module cctk_Math
  implicit none

  interface

     double precision function CCTK_copysign (x, y)
       implicit none
       double precision, intent(in) :: x, y
     end function CCTK_copysign

     integer function CCTK_fpclassify (x)
       implicit none
       double precision, intent(in) :: x
     end function CCTK_fpclassify

     integer function CCTK_isfinite (x)
       implicit none
       double precision, intent(in) :: x
     end function CCTK_isfinite

     integer function CCTK_isinf (x)
       implicit none
       double precision, intent(in) :: x
     end function CCTK_isinf

     integer function CCTK_isnan (x)
       implicit none
       double precision, intent(in) :: x
     end function CCTK_isnan

     integer function CCTK_isnormal (x)
       implicit none
       double precision, intent(in) :: x
     end function CCTK_isnormal

     integer function CCTK_signbit (x)
       implicit none
       double precision, intent(in) :: x
     end function CCTK_signbit

  end interface

end module cctk_Math
