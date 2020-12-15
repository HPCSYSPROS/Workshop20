#include "cctk.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

CCTK_REAL function GetCoeff ()

  implicit none

  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS

  CCTK_REAL, parameter :: zero = 0.0
  integer, parameter :: wp = kind(zero)

  if ( CCTK_EQUALS(norm_type,'Diagonal') ) then
    select case (order)
    case (2)
      GetCoeff = 1.0_wp / 2.0_wp
    case (4)
      GetCoeff = 17.0_wp/48.0_wp
    case (6)
      GetCoeff = 13649.0_wp/43200.0_wp
    case (8)
      GetCoeff =  1498139.0_wp/5080320.0_wp
   case default
      call CCTK_WARN (0, "operator not implemented")
    end select
  else
    if ( CCTK_EQUALS(operator_type,'Minimal Bandwidth') ) then
      select case (order)
      case(4)
        GetCoeff = 3.0_wp/11.0_wp
      case(6)
        GetCoeff = 30.0_wp/137.0_wp
     case default
        call CCTK_WARN (0, "operator not implemented")
      end select
    else
      select case (order)
      case(4)
        GetCoeff = 0.2388575707774486064323157210922003533466_wp
      case(6)
        GetCoeff = 0.2028105550720356346665604029847379994496_wp
     case default
        call CCTK_WARN (0, "operator not implemented")
      end select
    end if
  end if
end function GetCoeff
