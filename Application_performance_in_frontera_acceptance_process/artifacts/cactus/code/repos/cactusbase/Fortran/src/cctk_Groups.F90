#include "cctk.h"

module cctk_Groups
  implicit none

  interface

     subroutine CCTK_DecomposeName (ierr, fullname, implementation, implementation_nchars, name, name_nchars)
       implicit none
       integer      ierr
       character(*) fullname
       character(*) implementation
       integer      implementation_nchars
       character(*) name
       integer      name_nchars
     end subroutine CCTK_DecomposeName

     subroutine CCTK_FirstVarIndex (index, group)
       implicit none
       integer      index
       character(*) group
     end subroutine CCTK_FirstVarIndex

     subroutine CCTK_FirstVarIndexI (index, group)
       implicit none
       integer index
       integer group
     end subroutine CCTK_FirstVarIndexI

     subroutine CCTK_FullName (nchars, var, fullname)
       implicit none
       integer      nchars
       integer      var
       character(*) fullname
     end subroutine CCTK_FullName

     ! CCTK_GroupData fills a structure and has no Fortran wrapper
     
     subroutine CCTK_GroupDimI (dim, group)
       implicit none
       integer dim
       integer group
     end subroutine CCTK_GroupDimI
     
     subroutine CCTK_GroupDimFromVarI (dim, var)
       implicit none
       integer dim
       integer var
     end subroutine CCTK_GroupDimFromVarI

     subroutine CCTK_GroupDistribNumber (number, distrib)
       implicit none
       integer      number
       character(*) distrib
     end subroutine CCTK_GroupDistribNumber
     
     ! CCTK_GroupGhostsizesI is a strange function and has no Fortran wrapper

     subroutine CCTK_ImplementationI (nchars, group, implementation)
       implicit none
       integer      nchars
       integer      group
       character(*) implementation
     end subroutine CCTK_ImplementationI

     subroutine CCTK_GroupIndex (index, group)
       implicit none
       integer      index
       character(*) group
     end subroutine CCTK_GroupIndex

     subroutine CCTK_GroupIndexFromVar (index, var)
       implicit none
       integer      index
       character(*) var
     end subroutine CCTK_GroupIndexFromVar

     subroutine CCTK_GroupIndexFromVarI (index, var)
       implicit none
       integer index
       integer var
     end subroutine CCTK_GroupIndexFromVarI

     subroutine CCTK_GroupName (nchars, group, groupname)
       implicit none
       integer      nchars
       integer      group
       character(*) groupname
     end subroutine CCTK_GroupName
     
     subroutine CCTK_GroupNameFromVarI (nchars, var, groupname)
       implicit none
       integer      nchars
       integer      var
       character(*) groupname
     end subroutine CCTK_GroupNameFromVarI

     subroutine CCTK_GroupScopeNumber (number, scope)
       implicit none
       integer      number
       character(*) scope
     end subroutine CCTK_GroupScopeNumber
     
     ! CCTK_GroupSizesI is a strange function and has no Fortran wrapper
     
     subroutine CCTK_GroupTypeFromVarI (type, var)
       implicit none
       integer type
       integer var
     end subroutine CCTK_GroupTypeFromVarI

     subroutine CCTK_GroupTypeNumber (number, type)
       implicit none
       integer      number
       character(*) type
     end subroutine CCTK_GroupTypeNumber
     
     subroutine CCTK_GroupTypeI (type, group)
       implicit none
       integer type
       integer group
     end subroutine CCTK_GroupTypeI
     
     subroutine CCTK_ImpFromVarI (nchars, var, imp)
       implicit none
       integer      nchars
       integer      var
       character(*) imp
     end subroutine CCTK_ImpFromVarI
     
     subroutine CCTK_MaxDim (maxdim)
       implicit none
       integer maxdim
     end subroutine CCTK_MaxDim
     
     subroutine CCTK_NumGroups (numgroups)
       implicit none
       integer numgroups
     end subroutine CCTK_NumGroups
     
     subroutine CCTK_NumTimeLevelsFromVar (numtimelevels, var)
       implicit none
       integer      numtimelevels
       character(*) var
     end subroutine CCTK_NumTimeLevelsFromVar
     
     subroutine CCTK_NumTimeLevelsFromVarI (numtimelevels, var)
       implicit none
       integer numtimelevels
       integer var
     end subroutine CCTK_NumTimeLevelsFromVarI
     
     subroutine CCTK_NumTimeLevels (numtimelevels, var)
       implicit none
       integer      numtimelevels
       character(*) var
     end subroutine CCTK_NumTimeLevels
     
     subroutine CCTK_NumTimeLevelsI (numtimelevels, var)
       implicit none
       integer numtimelevels
       integer var
     end subroutine CCTK_NumTimeLevelsI
     
     subroutine CCTK_MaxTimeLevels (maxtimelevels, group)
       implicit none
       integer      maxtimelevels
       character(*) group
     end subroutine CCTK_MaxTimeLevels
     
     subroutine CCTK_MaxTimeLevelsVN (maxtimelevels, var)
       implicit none
       integer      maxtimelevels
       character(*) var
     end subroutine CCTK_MaxTimeLevelsVN
     
     subroutine CCTK_MaxTimeLevelsVI (maxtimelevels, var)
       implicit none
       integer maxtimelevels
       integer var
     end subroutine CCTK_MaxTimeLevelsVI
     
     subroutine CCTK_MaxTimeLevelsGN (maxtimelevels, group)
       implicit none
       integer      maxtimelevels
       character(*) group
     end subroutine CCTK_MaxTimeLevelsGN
     
     subroutine CCTK_MaxTimeLevelsGI (maxtimelevels, group)
       implicit none
       integer maxtimelevels
       integer group
     end subroutine CCTK_MaxTimeLevelsGI
     
     subroutine CCTK_NumVars (numvars)
       implicit none
       integer numvars
     end subroutine CCTK_NumVars
     
     subroutine CCTK_NumVarsInGroup (numvars, group)
       implicit none
       integer      numvars
       character(*) group
     end subroutine CCTK_NumVarsInGroup
     
     subroutine CCTK_NumVarsInGroupI (numvars, group)
       implicit none
       integer numvars
       integer group
     end subroutine CCTK_NumVarsInGroupI
     
     subroutine CCTK_VarIndex (index, var)
       implicit none
       integer      index
       character(*) var
     end subroutine CCTK_VarIndex
     
     subroutine CCTK_VarName (nchars, var, varname)
       implicit none
       integer      nchars
       integer      var
       character(*) varname
     end subroutine CCTK_VarName
     
     subroutine CCTK_VarTypeI (type, var)
       implicit none
       integer type
       integer var
     end subroutine CCTK_VarTypeI
     
     subroutine CCTK_VarTypeNumber (number, type)
       implicit none
       integer      number
       character(*) type
     end subroutine CCTK_VarTypeNumber
     
     subroutine CCTK_VarTypeName (nchars, type, typename)
       implicit none
       integer      nchars
       integer      type
       character(*) typename
     end subroutine CCTK_VarTypeName
     
     subroutine CCTK_VarTypeSize (size, type)
       implicit none
       integer size
       integer type
     end subroutine CCTK_VarTypeSize
     
     ! CCTKi_GroupLengthAsPointer is a strange function and has no
     ! Fortran wrapper

     ! CCTK_TraverseString has no Fortran wrapper
     
     subroutine CCTK_GroupTagsTable (table, group)
       implicit none
       integer      table
       character(*) group
     end subroutine CCTK_GroupTagsTable
     
     subroutine CCTK_GroupTagsTableI (table, group)
       implicit none
       integer table
       integer group
     end subroutine CCTK_GroupTagsTableI
     
  end interface
  
end module cctk_Groups
