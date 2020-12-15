#include "cctk.h"

module util_Table
  implicit none

  interface

     ! create/destroy

     subroutine Util_TableCreate (handle, flags)
       implicit none
       integer               handle
       integer               flags
     end subroutine Util_TableCreate

     subroutine Util_TableClone (ierr, handle)
       implicit none
       integer               ierr
       integer               handle
     end subroutine Util_TableClone

     subroutine Util_TableDestroy (ierr, handle)
       implicit none
       integer               ierr
       integer               handle
     end subroutine Util_TableDestroy

     ! query

     subroutine Util_TableQueryFlags (flags, handle)
       implicit none
       integer               flags
       integer               handle
     end subroutine Util_TableQueryFlags

     subroutine Util_TableQueryNKeys (nkeys, handle)
       implicit none
       integer               nkeys
       integer               handle
     end subroutine Util_TableQueryNKeys

     subroutine Util_TableQueryMaxKeyLength (maxkeylength, handle)
       implicit none
       integer               maxkeylength
       integer               handle
     end subroutine Util_TableQueryMaxKeyLength

     subroutine Util_TableQueryValueInfo &
          (keyinfo, handle, type_code, N_elements, key)
       implicit none
       integer               keyinfo
       integer               handle
       CCTK_INT              type_code
       CCTK_INT              N_elements
       character(*)          key
     end subroutine Util_TableQueryValueInfo

     ! misc stuff

     subroutine Util_TableDeleteKey (ierr, handle, key)
       implicit none
       integer               ierr
       integer               handle
       character(*)          key
     end subroutine Util_TableDeleteKey

     ! convenience routines to create and/or set from a "parameter-file" string

     subroutine Util_TableCreateFromString (handle, string)
       implicit none
       integer               handle
       character(*)          string
     end subroutine Util_TableCreateFromString

     ! set/get a C-style null-terminated character string
     subroutine Util_TableSetString (info, handle, string, key)
       implicit none
       integer               info
       integer               handle
       character(*)          string
       character(*)          key
     end subroutine Util_TableSetString

     subroutine Util_TableGetString (length, handle, buffer, key)
       implicit none
       integer               length
       integer               handle
       character(*)          buffer
       character(*)          key
     end subroutine Util_TableGetString

     ! set/get generic types described by CCTK_VARIABLE_* type codes

     subroutine Util_TableSetGeneric (info, handle, type_code, value_ptr, key)
       implicit none
       integer               info
       integer               handle
       integer               type_code
       CCTK_POINTER_TO_CONST value_ptr
       character(*)          key
     end subroutine Util_TableSetGeneric

     subroutine Util_TableSetGenericArray &
          (info, handle, type_code, N_elements, array, key)
       implicit none
       integer               info
       integer               handle
       integer               type_code
       integer               N_elements
       CCTK_POINTER_TO_CONST array(*)
       character(*)          key
     end subroutine Util_TableSetGenericArray

     subroutine Util_TableGetGeneric &
          (length, handle, type_code, value_ptr, key)
       implicit none
       integer               length
       integer               handle
       integer               type_code
       CCTK_POINTER          value_ptr
       character(*)          key
     end subroutine Util_TableGetGeneric

     subroutine Util_TableGetGenericArray &
          (length, handle, type_code, N_elements, array, key)
       implicit none
       integer               length
       integer               handle
       integer               type_code
       integer               N_elements
       CCTK_POINTER          array(*)
       character(*)          key
     end subroutine Util_TableGetGenericArray

     ! set routines

     ! pointers

     subroutine Util_TableSetPointer (info, handle, value, key)
       implicit none
       integer               info
       integer               handle
       CCTK_POINTER          value
       character(*)          key
     end subroutine Util_TableSetPointer

     subroutine Util_TableSetFPointer (info, handle, value, key)
       implicit none
       integer               info
       integer               handle
       CCTK_FPOINTER         value
       character(*)          key
     end subroutine Util_TableSetFPointer

     ! integers

     subroutine Util_TableSetByte (info, handle, value, key)
       implicit none
       integer               info
       integer               handle
       CCTK_BYTE              value
       character(*)          key
     end subroutine Util_TableSetByte

     subroutine Util_TableSetInt (info, handle, value, key)
       implicit none
       integer               info
       integer               handle
       CCTK_INT              value
       character(*)          key
     end subroutine Util_TableSetInt

#ifdef CCTK_INT1
     subroutine Util_TableSetInt1 (info, handle, value, key)
       implicit none
       integer               info
       integer               handle
       CCTK_INT1             value
       character(*)          key
     end subroutine Util_TableSetInt1
#endif

#ifdef CCTK_INT2
     subroutine Util_TableSetInt2 (info, handle, value, key)
       implicit none
       integer               info
       integer               handle
       CCTK_INT2             value
       character(*)          key
     end subroutine Util_TableSetInt2
#endif

#ifdef CCTK_INT4
     subroutine Util_TableSetInt4 (info, handle, value, key)
       implicit none
       integer               info
       integer               handle
       CCTK_INT4             value
       character(*)          key
     end subroutine Util_TableSetInt4
#endif

#ifdef CCTK_INT8
     subroutine Util_TableSetInt8 (info, handle, value, key)
       implicit none
       integer               info
       integer               handle
       CCTK_INT8             value
       character(*)          key
     end subroutine Util_TableSetInt8
#endif

     ! real numbers

     subroutine Util_TableSetReal (info, handle, value, key)
       implicit none
       integer               info
       integer               handle
       CCTK_REAL             value
       character(*)          key
     end subroutine Util_TableSetReal

#ifdef HAVE_CCTK_REAL4
     subroutine Util_TableSetReal4 (info, handle, value, key)
       implicit none
       integer               info
       integer               handle
       CCTK_REAL4            value
       character(*)          key
     end subroutine Util_TableSetReal4
#endif

#ifdef HAVE_CCTK_REAL8
     subroutine Util_TableSetReal8 (info, handle, value, key)
       implicit none
       integer               info
       integer               handle
       CCTK_REAL8            value
       character(*)          key
     end subroutine Util_TableSetReal8
#endif

#ifdef HAVE_CCTK_REAL16
     subroutine Util_TableSetReal16 (info, handle, value, key)
       implicit none
       integer               info
       integer               handle
       CCTK_REAL16           value
       character(*)          key
     end subroutine Util_TableSetReal16
#endif

     ! complex numbers

     subroutine Util_TableSetComplex (info, handle, value, key)
       implicit none
       integer               info
       integer               handle
       CCTK_COMPLEX          value
       character(*)          key
     end subroutine Util_TableSetComplex

#ifdef HAVE_CCTK_REAL4
     subroutine Util_TableSetComplex8 (info, handle, value, key)
       implicit none
       integer               info
       integer               handle
       CCTK_COMPLEX8         value
       character(*)          key
     end subroutine Util_TableSetComplex8
#endif

#ifdef HAVE_CCTK_REAL8
     subroutine Util_TableSetComplex16 (info, handle, value, key)
       implicit none
       integer               info
       integer               handle
       CCTK_COMPLEX16        value
       character(*)          key
     end subroutine Util_TableSetComplex16
#endif

#ifdef HAVE_CCTK_REAL16
     subroutine Util_TableSetComplex32 (info, handle, value, key)
       implicit none
       integer               info
       integer               handle
       CCTK_COMPLEX32        value
       character(*)          key
     end subroutine Util_TableSetComplex32
#endif

     ! arrays of pointers

     subroutine Util_TableSetPointerArray &
          (info, handle, N_elements, array, key)
       implicit none
       integer               info
       integer               handle
       integer               N_elements
       CCTK_POINTER          array(*)
       character(*)          key
     end subroutine Util_TableSetPointerArray

     subroutine Util_TableSetFPointerArray &
          (info, handle, N_elements, array, key)
       implicit none
       integer               info
       integer               handle
       integer               N_elements
       CCTK_FPOINTER         array(*)
       character(*)          key
     end subroutine Util_TableSetFPointerArray

     ! arrays of integers

     subroutine Util_TableSetByteArray (info, handle, N_elements, array, key)
       implicit none
       integer               info
       integer               handle
       integer               N_elements
       CCTK_BYTE              array(*)
       character(*)          key
     end subroutine Util_TableSetByteArray

     subroutine Util_TableSetIntArray (info, handle, N_elements, array, key)
       implicit none
       integer               info
       integer               handle
       integer               N_elements
       CCTK_INT              array(*)
       character(*)          key
     end subroutine Util_TableSetIntArray

#ifdef CCTK_INT1
     subroutine Util_TableSetInt1Array (info, handle, N_elements, array, key)
       implicit none
       integer               info
       integer               handle
       integer               N_elements
       CCTK_INT1             array(*)
       character(*)          key
     end subroutine Util_TableSetInt1Array
#endif

#ifdef CCTK_INT2
     subroutine Util_TableSetInt2Array (info, handle, N_elements, array, key)
       implicit none
       integer               info
       integer               handle
       integer               N_elements
       CCTK_INT2             array(*)
       character(*)          key
     end subroutine Util_TableSetInt2Array
#endif

#ifdef CCTK_INT4
     subroutine Util_TableSetInt4Array (info, handle, N_elements, array, key)
       implicit none
       integer               info
       integer               handle
       integer               N_elements
       CCTK_INT4             array(*)
       character(*)          key
     end subroutine Util_TableSetInt4Array
#endif

#ifdef CCTK_INT8
     subroutine Util_TableSetInt8Array (info, handle, N_elements, array, key)
       implicit none
       integer               info
       integer               handle
       integer               N_elements
       CCTK_INT8             array(*)
       character(*)          key
     end subroutine Util_TableSetInt8Array
#endif

     ! arrays of real numbers

     subroutine Util_TableSetRealArray (info, handle, N_elements, array, key)
       implicit none
       integer               info
       integer               handle
       integer               N_elements
       CCTK_REAL             array(*)
       character(*)          key
     end subroutine Util_TableSetRealArray

#ifdef HAVE_CCTK_REAL4
     subroutine Util_TableSetReal4Array (info, handle, N_elements, array, key)
       implicit none
       integer               info
       integer               handle
       integer               N_elements
       CCTK_REAL4            array(*)
       character(*)          key
     end subroutine Util_TableSetReal4Array
#endif

#ifdef HAVE_CCTK_REAL8
     subroutine Util_TableSetReal8Array (info, handle, N_elements, array, key)
       implicit none
       integer               info
       integer               handle
       integer               N_elements
       CCTK_REAL8            array(*)
       character(*)          key
     end subroutine Util_TableSetReal8Array
#endif

#ifdef HAVE_CCTK_REAL16
     subroutine Util_TableSetReal16Array (info, handle, N_elements, array, key)
       implicit none
       integer               info
       integer               handle
       integer               N_elements
       CCTK_REAL16           array(*)
       character(*)          key
     end subroutine Util_TableSetReal16Array
#endif

     ! arrays of complex numbers

     subroutine Util_TableSetComplexArray &
          (info, handle, N_elements, array, key)
       implicit none
       integer               info
       integer               handle
       integer               N_elements
       CCTK_COMPLEX          array(*)
       character(*)          key
     end subroutine Util_TableSetComplexArray

#ifdef HAVE_CCTK_REAL4
     subroutine Util_TableSetComplex8Array &
          (info, handle, N_elements, array, key)
       implicit none
       integer               info
       integer               handle
       integer               N_elements
       CCTK_COMPLEX8         array(*)
       character(*)          key
     end subroutine Util_TableSetComplex8Array
#endif

#ifdef HAVE_CCTK_REAL8
     subroutine Util_TableSetComplex16Array &
          (info, handle, N_elements, array, key)
       implicit none
       integer               info
       integer               handle
       integer               N_elements
       CCTK_COMPLEX16         array(*)
       character(*)          key
     end subroutine Util_TableSetComplex16Array
#endif

#ifdef HAVE_CCTK_REAL16
     subroutine Util_TableSetComplex32Array &
          (info, handle, N_elements, array, key)
       implicit none
       integer               info
       integer               handle
       integer               N_elements
       CCTK_COMPLEX32         array(*)
       character(*)          key
     end subroutine Util_TableSetComplex32Array
#endif

     ! get routines

     ! pointers

     subroutine Util_TableGetPointer (length, handle, value, key)
       implicit none
       integer               length
       integer               handle
       CCTK_POINTER          value
       character(*)          key
     end subroutine Util_TableGetPointer

     subroutine Util_TableGetFPointer (length, handle, value, key)
       implicit none
       integer               length
       integer               handle
       CCTK_FPOINTER         value
       character(*)          key
     end subroutine Util_TableGetFPointer

     ! integers

     subroutine Util_TableGetByte (length, handle, value, key)
       implicit none
       integer               length
       integer               handle
       CCTK_BYTE              value
       character(*)          key
     end subroutine Util_TableGetByte

     subroutine Util_TableGetInt (length, handle, value, key)
       implicit none
       integer               length
       integer               handle
       CCTK_INT              value
       character(*)          key
     end subroutine Util_TableGetInt

#ifdef CCTK_INT1
     subroutine Util_TableGetInt1 (length, handle, value, key)
       implicit none
       integer               length
       integer               handle
       CCTK_INT1             value
       character(*)          key
     end subroutine Util_TableGetInt1
#endif

#ifdef CCTK_INT2
     subroutine Util_TableGetInt2 (length, handle, value, key)
       implicit none
       integer               length
       integer               handle
       CCTK_INT2             value
       character(*)          key
     end subroutine Util_TableGetInt2
#endif

#ifdef CCTK_INT4
     subroutine Util_TableGetInt4 (length, handle, value, key)
       implicit none
       integer               length
       integer               handle
       CCTK_INT4             value
       character(*)          key
     end subroutine Util_TableGetInt4
#endif

#ifdef CCTK_INT8
     subroutine Util_TableGetInt8 (length, handle, value, key)
       implicit none
       integer               length
       integer               handle
       CCTK_INT8             value
       character(*)          key
     end subroutine Util_TableGetInt8
#endif

     ! real numbers

     subroutine Util_TableGetReal (length, handle, value, key)
       implicit none
       integer               length
       integer               handle
       CCTK_REAL             value
       character(*)          key
     end subroutine Util_TableGetReal

#ifdef HAVE_CCTK_REAL4
     subroutine Util_TableGetReal4 (length, handle, value, key)
       implicit none
       integer               length
       integer               handle
       CCTK_REAL4            value
       character(*)          key
     end subroutine Util_TableGetReal4
#endif

#ifdef HAVE_CCTK_REAL8
     subroutine Util_TableGetReal8 (length, handle, value, key)
       implicit none
       integer               length
       integer               handle
       CCTK_REAL8            value
       character(*)          key
     end subroutine Util_TableGetReal8
#endif

#ifdef HAVE_CCTK_REAL16
     subroutine Util_TableGetReal16 (length, handle, value, key)
       implicit none
       integer               length
       integer               handle
       CCTK_REAL16           value
       character(*)          key
     end subroutine Util_TableGetReal16
#endif

     ! complex numbers

     subroutine Util_TableGetComplex (length, handle, value, key)
       implicit none
       integer               length
       integer               handle
       CCTK_COMPLEX          value
       character(*)          key
     end subroutine Util_TableGetComplex

#ifdef HAVE_CCTK_REAL4
     subroutine Util_TableGetComplex8 (length, handle, value, key)
       implicit none
       integer               length
       integer               handle
       CCTK_COMPLEX8         value
       character(*)          key
     end subroutine Util_TableGetComplex8
#endif

#ifdef HAVE_CCTK_REAL8
     subroutine Util_TableGetComplex16 (length, handle, value, key)
       implicit none
       integer               length
       integer               handle
       CCTK_COMPLEX16        value
       character(*)          key
     end subroutine Util_TableGetComplex16
#endif

#ifdef HAVE_CCTK_REAL16
     subroutine Util_TableGetComplex32 (length, handle, value, key)
       implicit none
       integer               length
       integer               handle
       CCTK_COMPLEX32        value
       character(*)          key
     end subroutine Util_TableGetComplex32
#endif

     ! arrays of pointers

     subroutine Util_TableGetPointerArray &
          (length, handle, N_elements, array, key)
       implicit none
       integer               length
       integer               handle
       integer               N_elements
       CCTK_POINTER          array(*)
       character(*)          key
     end subroutine Util_TableGetPointerArray

     subroutine Util_TableGetFPointerArray &
          (length, handle, N_elements, array, key)
       implicit none
       integer               length
       integer               handle
       integer               N_elements
       CCTK_FPOINTER         array(*)
       character(*)          key
     end subroutine Util_TableGetFPointerArray

     ! integers

     subroutine Util_TableGetByteArray (length, handle, N_elements, array, key)
       implicit none
       integer               length
       integer               handle
       integer               N_elements
       CCTK_BYTE              array(*)
       character(*)          key
     end subroutine Util_TableGetByteArray

     subroutine Util_TableGetIntArray (length, handle, N_elements, array, key)
       implicit none
       integer               length
       integer               handle
       integer               N_elements
       CCTK_INT              array(*)
       character(*)          key
     end subroutine Util_TableGetIntArray

#ifdef CCTK_INT1
     subroutine Util_TableGetInt1Array (length, handle, N_elements, array, key)
       implicit none
       integer               length
       integer               handle
       integer               N_elements
       CCTK_INT1             array(*)
       character(*)          key
     end subroutine Util_TableGetInt1Array
#endif

#ifdef CCTK_INT2
     subroutine Util_TableGetInt2Array (length, handle, N_elements, array, key)
       implicit none
       integer               length
       integer               handle
       integer               N_elements
       CCTK_INT2             array(*)
       character(*)          key
     end subroutine Util_TableGetInt2Array
#endif

#ifdef CCTK_INT4
     subroutine Util_TableGetInt4Array (length, handle, N_elements, array, key)
       implicit none
       integer               length
       integer               handle
       integer               N_elements
       CCTK_INT4             array(*)
       character(*)          key
     end subroutine Util_TableGetInt4Array
#endif

#ifdef CCTK_INT8
     subroutine Util_TableGetInt8Array (length, handle, N_elements, array, key)
       implicit none
       integer               length
       integer               handle
       integer               N_elements
       CCTK_INT8             array(*)
       character(*)          key
     end subroutine Util_TableGetInt8Array
#endif

     ! real numbers

     subroutine Util_TableGetRealArray (length, handle, N_elements, array, key)
       implicit none
       integer               length
       integer               handle
       integer               N_elements
       CCTK_REAL             array(*)
       character(*)          key
     end subroutine Util_TableGetRealArray

#ifdef HAVE_CCTK_REAL4
     subroutine Util_TableGetReal4Array &
          (length, handle, N_elements, array, key)
       implicit none
       integer               length
       integer               handle
       integer               N_elements
       CCTK_REAL4            array(*)
       character(*)          key
     end subroutine Util_TableGetReal4Array
#endif

#ifdef HAVE_CCTK_REAL8
     subroutine Util_TableGetReal8Array &
          (length, handle, N_elements, array, key)
       implicit none
       integer               length
       integer               handle
       integer               N_elements
       CCTK_REAL8            array(*)
       character(*)          key
     end subroutine Util_TableGetReal8Array
#endif

#ifdef HAVE_CCTK_REAL16
     subroutine Util_TableGetReal16Array &
          (length, handle, N_elements, array, key)
       implicit none
       integer               length
       integer               handle
       integer               N_elements
       CCTK_REAL16           array(*)
       character(*)          key
     end subroutine Util_TableGetReal16Array
#endif

     ! complex numbers

     subroutine Util_TableGetComplexArray &
          (length, handle, N_elements, array, key)
       implicit none
       integer               length
       integer               handle
       integer               N_elements
       CCTK_COMPLEX          array(*)
       character(*)          key
     end subroutine Util_TableGetComplexArray

#ifdef HAVE_CCTK_REAL4
     subroutine Util_TableGetComplex8Array &
          (length, handle, N_elements, array, key)
       implicit none
       integer               length
       integer               handle
       integer               N_elements
       CCTK_COMPLEX8         array(*)
       character(*)          key
     end subroutine Util_TableGetComplex8Array
#endif

#ifdef HAVE_CCTK_REAL8
     subroutine Util_TableGetComplex16Array &
          (length, handle, N_elements, array, key)
       implicit none
       integer               length
       integer               handle
       integer               N_elements
       CCTK_COMPLEX16        array(*)
       character(*)          key
     end subroutine Util_TableGetComplex16Array
#endif

#ifdef HAVE_CCTK_REAL16
     subroutine Util_TableGetComplex32Array &
          (length, handle, N_elements, array, key)
       implicit none
       integer               length
       integer               handle
       integer               N_elements
       CCTK_COMPLEX32        array(*)
       character(*)          key
     end subroutine Util_TableGetComplex32Array
#endif

     ! Table Iterator API

     ! create/destroy

     subroutine Util_TableItCreate (ihandle, handle)
       implicit none
       integer               ihandle
       integer               handle
     end subroutine Util_TableItCreate

     subroutine Util_TableItClone (newihandle, ihandle)
       implicit none
       integer               newihandle
       integer               ihandle
     end subroutine Util_TableItClone

     subroutine Util_TableItDestroy (ierr, ihandle)
       implicit none
       integer               ierr
       integer               ihandle
     end subroutine Util_TableItDestroy

     ! test for "null-pointer" state

     subroutine Util_TableItQueryIsNull (info, ihandle)
       implicit none
       integer               info
       integer               ihandle
     end subroutine Util_TableItQueryIsNull

     subroutine Util_TableItQueryIsNonNull (info, ihandle)
       implicit none
       integer               info
       integer               ihandle
     end subroutine Util_TableItQueryIsNonNull

     ! query what the iterator points to

     subroutine Util_TableItQueryTableHandle (handle, ihandle)
       implicit none
       integer               handle
       integer               ihandle
     end subroutine Util_TableItQueryTableHandle

     subroutine Util_TableItQueryKeyValueInfo &
          (length, ihandle, key_buffer, type_code, N_elements)
       implicit none
       integer               length
       integer               ihandle
       character(*)          key_buffer
       CCTK_INT              type_code
       CCTK_INT              N_elements
     end subroutine Util_TableItQueryKeyValueInfo

     ! change value of iterator

     subroutine Util_TableItAdvance (info, ihandle)
       implicit none
       integer               info
       integer               ihandle
     end subroutine Util_TableItAdvance

     subroutine Util_TableItResetToStart (info, ihandle)
       implicit none
       integer               info
       integer               ihandle
     end subroutine Util_TableItResetToStart

     subroutine Util_TableItSetToNull (ierr, ihandle)
       implicit none
       integer               ierr
       integer               ihandle
     end subroutine Util_TableItSetToNull

     subroutine Util_TableItSetToKey (ierr, ihandle, key)
       implicit none
       integer               ierr
       integer               ihandle
       character(*)          key
     end subroutine Util_TableItSetToKey

  end interface

end module util_Table
