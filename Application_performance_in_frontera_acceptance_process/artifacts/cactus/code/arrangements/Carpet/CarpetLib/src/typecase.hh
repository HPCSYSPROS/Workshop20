// Instantiate type cases for all available types    -*-C++-*-
// (C) 2001 Erik Schnetter <schnetter@uni-tuebingen.de>

// Usage:
// Define the macro TYPECASE(N,T) to be a typecase for the type T with name N,
// then include this file,
// then undefine the macro TYPECASE.

// Decide which types to typecase

// Should all types be used?
#ifdef CARPET_ALL
#undef CARPET_BYTE
#undef CARPET_ALL_INT
#undef CARPET_ALL_REAL
#undef CARPET_ALL_COMPLEX
#define CARPET_BYTE
#define CARPET_ALL_INT
#define CARPET_ALL_REAL
#define CARPET_ALL_COMPLEX
#endif

// Should all integer/real/complex types be used?
#ifdef CARPET_ALL_INT
#undef CARPET_INT1
#undef CARPET_INT2
#undef CARPET_INT4
#undef CARPET_INT8
#undef CARPET_INT16
#define CARPET_INT1
#define CARPET_INT2
#define CARPET_INT4
#define CARPET_INT8
#define CARPET_INT16
#endif

#ifdef CARPET_ALL_REAL
#undef CARPET_REAL4
#undef CARPET_REAL8
#undef CARPET_REAL16
#define CARPET_REAL4
#define CARPET_REAL8
#define CARPET_REAL16
#endif

#ifdef CARPET_ALL_COMPLEX
#undef CARPET_COMPLEX8
#undef CARPET_COMPLEX16
#undef CARPET_COMPLEX32
#define CARPET_COMPLEX8
#define CARPET_COMPLEX16
#define CARPET_COMPLEX32
#endif

// If no types are specified, use a sensible default
#if !defined(CARPET_BYTE) && !defined(CARPET_INT) && !defined(CARPET_INT1) &&  \
    !defined(CARPET_INT2) && !defined(CARPET_INT4) && !defined(CARPET_INT8) && \
    !defined(CARPET_INT16) && !defined(CARPET_REAL) &&                         \
    !defined(CARPET_REAL4) && !defined(CARPET_REAL8) &&                        \
    !defined(CARPET_REAL16) && !defined(CARPET_COMPLEX) &&                     \
    !defined(CARPET_COMPLEX8) && !defined(CARPET_COMPLEX16) &&                 \
    !defined(CARPET_COMPLEX32)
// Assume the user just wants INT, REAL, and COMPLEX
#define CARPET_INT
#define CARPET_REAL
#define CARPET_COMPLEX
#endif

// Translate the default types to their specific counterparts
#ifdef CARPET_INT
#ifdef CCTK_INTEGER_PRECISION_1
#undef CARPET_INT1
#define CARPET_INT1
#endif
#ifdef CCTK_INTEGER_PRECISION_2
#undef CARPET_INT2
#define CARPET_INT2
#endif
#ifdef CCTK_INTEGER_PRECISION_4
#undef CARPET_INT4
#define CARPET_INT4
#endif
#ifdef CCTK_INTEGER_PRECISION_8
#undef CARPET_INT8
#define CARPET_INT8
#endif
#ifdef CCTK_INTEGER_PRECISION_16
#undef CARPET_INT16
#define CARPET_INT16
#endif
#endif
#ifdef CARPET_REAL
#ifdef CCTK_REAL_PRECISION_4
#undef CARPET_REAL4
#define CARPET_REAL4
#endif
#ifdef CCTK_REAL_PRECISION_8
#undef CARPET_REAL8
#define CARPET_REAL8
#endif
#ifdef CCTK_REAL_PRECISION_16
#undef CARPET_REAL16
#define CARPET_REAL16
#endif
#endif
#ifdef CARPET_COMPLEX
#ifdef CCTK_REAL_PRECISION_4
#undef CARPET_COMPLEX8
#define CARPET_COMPLEX8
#endif
#ifdef CCTK_REAL_PRECISION_8
#undef CARPET_COMPLEX16
#define CARPET_COMPLEX16
#endif
#ifdef CCTK_REAL_PRECISION_16
#undef CARPET_COMPLEX32
#define CARPET_COMPLEX32
#endif
#endif

// Typecase the desired types

#ifndef CARPET_NO_INT

#ifdef CARPET_BYTE
TYPECASE(CCTK_VARIABLE_BYTE, CCTK_BYTE)
#endif

// #ifdef CARPET_INT
// TYPECASE(CCTK_VARIABLE_INT, CCTK_INT)
// #endif
#ifdef CARPET_INT1
#ifdef HAVE_CCTK_INT1
TYPECASE(CCTK_VARIABLE_INT1, CCTK_INT1)
#endif
#endif
#ifdef CARPET_INT2
#ifdef HAVE_CCTK_INT2
TYPECASE(CCTK_VARIABLE_INT2, CCTK_INT2)
#endif
#endif
#ifdef CARPET_INT4
#ifdef HAVE_CCTK_INT4
TYPECASE(CCTK_VARIABLE_INT4, CCTK_INT4)
#endif
#endif
#ifdef CARPET_INT8
#ifdef HAVE_CCTK_INT8
TYPECASE(CCTK_VARIABLE_INT8, CCTK_INT8)
#endif
#endif
#ifdef CARPET_INT16
#ifdef HAVE_CCTK_INT16
TYPECASE(CCTK_VARIABLE_INT16, CCTK_INT16)
#endif
#endif

#endif

#ifndef CARPET_NO_REAL

// #ifdef CARPET_REAL
// TYPECASE(CCTK_VARIABLE_REAL, CCTK_REAL)
// #endif
#ifdef CARPET_REAL4
#ifdef HAVE_CCTK_REAL4
TYPECASE(CCTK_VARIABLE_REAL4, CCTK_REAL4)
#endif
#endif
#ifdef CARPET_REAL8
#ifdef HAVE_CCTK_REAL8
TYPECASE(CCTK_VARIABLE_REAL8, CCTK_REAL8)
#endif
#endif
#ifdef CARPET_REAL16
#ifdef HAVE_CCTK_REAL16
TYPECASE(CCTK_VARIABLE_REAL16, CCTK_REAL16)
#endif
#endif

#endif

#ifndef CARPET_NO_COMPLEX

// #  ifdef CARPET_COMPLEX
// TYPECASE(CCTK_VARIABLE_COMPLEX, CCTK_COMPLEX)
// #  endif
#ifdef CARPET_COMPLEX8
#ifdef HAVE_CCTK_COMPLEX8
TYPECASE(CCTK_VARIABLE_COMPLEX8, CCTK_COMPLEX8)
#endif
#endif
#ifdef CARPET_COMPLEX16
#ifdef HAVE_CCTK_COMPLEX16
TYPECASE(CCTK_VARIABLE_COMPLEX16, CCTK_COMPLEX16)
#endif
#endif
#ifdef CARPET_COMPLEX32
#ifdef HAVE_CCTK_COMPLEX32
TYPECASE(CCTK_VARIABLE_COMPLEX32, CCTK_COMPLEX32)
#endif
#endif

#endif

// Unset all macros
#undef CARPET_ALL
#undef CARPET_BYTE
#undef CARPET_ALL_INT
#undef CARPET_ALL_REAL
#undef CARPET_ALL_COMPLEX

#undef CARPET_NO_INT
#undef CARPET_NO_REAL
#undef CARPET_NO_COMPLEX

#undef CARPET_BYTE
#undef CARPET_INT
#undef CARPET_INT1
#undef CARPET_INT2
#undef CARPET_INT4
#undef CARPET_INT8
#undef CARPET_INT16
#undef CARPET_REAL
#undef CARPET_REAL4
#undef CARPET_REAL8
#undef CARPET_REAL16
#undef CARPET_COMPLEX
#undef CARPET_COMPLEX8
#undef CARPET_COMPLEX16
#undef CARPET_COMPLEX32
