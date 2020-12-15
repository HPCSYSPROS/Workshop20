/* cctk_Config.h.  Generated automatically by configure.  */
#ifndef _CCTK_CONFIG_H_
#define _CCTK_CONFIG_H_

/* Define if your processor stores words with the most significant
   byte first (like Motorola and SPARC, unlike Intel and VAX).  */
/* #undef WORDS_BIGENDIAN */

/* Define if standard C headers are available. */
#define STDC_HEADERS 1

/* Modifier for Fortran function definitions. */
#define CCTK_FCALL 

/* What debugging options to use */
/* #undef CCTK_DEBUG */
/* #undef CCTK_TRACEMEMORY */

/* Various library functions */
#define HAVE_GETHOSTBYNAME 1
#define HAVE_GETOPT_LONG_ONLY 1
#define HAVE___CXA_DEMANGLE 1
#define HAVE_DLADDR 1
#define HAVE_BACKTRACE 1
#define HAVE_BACKTRACE_SYMBOLS 1
#define HAVE_CRYPT 1
#define HAVE_FINITE 1
#define HAVE_COPYSIGN 1
#define HAVE_FPCLASSIFY 1
#define HAVE_ISFINITE 1
#define HAVE_ISINF 1
#define HAVE_ISNAN 1
#define HAVE_ISNORMAL 1
#define HAVE_SIGNBIT 1
#define HAVE_MKSTEMP 1
#define HAVE_VA_COPY 1

/* Do we have mode_t ? */
#define HAVE_MODE_T 1

/* Do we have SOCKET ? */
/* #undef HAVE_SOCKET_TYPE */

/* Do we have socklen_t ? Default to 'int' if not. */
#define HAVE_SOCKLEN_T 1
#ifdef HAVE_SOCKLEN_T
#  define CCTK_SOCKLEN_T socklen_t
#else
#  define CCTK_SOCKLEN_T int
#endif

/* Do we have hrtime_t ? */
/* #undef HAVE_HRTIME_T */

/* Some timing functions */
/* #undef HAVE_GETHRTIME */
/* #undef HAVE_READ_REAL_TIME */
/* #undef HAVE_TIME_BASE_TO_TIME */
#define HAVE_CLOCK_GETTIME 1
/* #undef HAVE_MACH_ABSOLUTE_TIME */

/* Cray UNICOS _rtc() (real-time clock) intrinsic */
/* #undef HAVE__RTC */

/* Some include things */
#define HAVE_TIME_H 1
/* #undef HAVE_SYS_FILIO_H */
#define HAVE_SYS_IOCTL_H 1
#define HAVE_SYS_SOCKET_H 1
#define HAVE_SYS_TIME_H 1
#define HAVE_SYS_TYPES_H 1
#define HAVE_UNISTD_H 1
#define HAVE_STRING_H 1
#define HAVE_ASSERT_H 1
#define HAVE_TGMATH_H 1
#define HAVE_SYS_STAT_H 1
#define HAVE_GETOPT_H 1
#define HAVE_REGEX_H 1
#define HAVE_NETINET_IN_H 1
#define HAVE_NETDB_H 1
#define HAVE_ARPA_INET_H 1
/* #undef HAVE_WINSOCK2_H */
#define HAVE_CRYPT_H 1
#define HAVE_DIRENT_H 1
/* #undef HAVE_C_ASM_H */
/* #undef HAVE_INTRINSICS_H */
/* #undef HAVE_MACH_MACH_TIME_H */
#define HAVE_SIGNAL_H 1
#define HAVE_MALLOC_H 1
#define HAVE_MALLINFO 1
#define HAVE_MALLOPT 1
#define HAVE_M_MMAP_THRESHOLD_VALUE 1
#define HAVE_SCHED_H 1
#define HAVE_EXECINFO_H 1
#define HAVE_SCHED_GETAFFINITY 1
#define HAVE_GETPID 1

#define TIME_WITH_SYS_TIME 1

/* Which format does the C++ STL on this machine provide. */
#define HAVE_VECTOR 1
/* #undef HAVE_VECTOR_H */

/* Timing stuff */
#define HAVE_TIME_GETTIMEOFDAY 1
#define GETTIMEOFDAY_NEEDS_TIMEZONE 1
#define HAVE_TIME_GETRUSAGE 1
/* #undef HAVE_TIME__FTIME */

/* The name of the NULL device for redirecting things to */
#define NULL_DEVICE "/dev/null"

/* Define the machine architecture for the build */
#define CCTK_BUILD_OS "linux-gnu"
#define CCTK_BUILD_CPU "x86_64"
#define CCTK_BUILD_VENDOR "unknown"

/******************************************************************************/

/* Size info for various datatypes */
#define SIZEOF_LONG_LONG 8
#define SIZEOF_LONG_INT 8
#define SIZEOF_INT 4
#define SIZEOF_SHORT_INT 2
#define SIZEOF_LONG_DOUBLE 16
#define SIZEOF_DOUBLE 8
#define SIZEOF_FLOAT 4
#define SIZEOF_CHAR_P 8

/* The chosen CCTK precision */

/* Floating point precision */
/* #undef CCTK_REAL_PRECISION_16 */
#define CCTK_REAL_PRECISION_8 1
/* #undef CCTK_REAL_PRECISION_4 */

/* Integer precision */
/* #undef CCTK_INTEGER_PRECISION_16 */
/* #undef CCTK_INTEGER_PRECISION_8 */
#define CCTK_INTEGER_PRECISION_4 1
/* #undef CCTK_INTEGER_PRECISION_2 */
/* #undef CCTK_INTEGER_PRECISION_1 */

/* Integer sizes */
#define HAVE_CCTK_INT16 1
#define HAVE_CCTK_INT8 1
#define HAVE_CCTK_INT4 1
#define HAVE_CCTK_INT2 1
#define HAVE_CCTK_INT1 1

/* Float sizes */
#define HAVE_CCTK_REAL16 1
#define HAVE_CCTK_REAL8 1
#define HAVE_CCTK_REAL4 1

/******************************************************************************/

#ifdef CCODE

/* CCTK C/C++ Integer datatypes */
#define CCTK_INT16_TYPE __int128
#define CCTK_INT8_TYPE long int
#define CCTK_INT4_TYPE int
#define CCTK_INT2_TYPE short int
#define CCTK_INT1_TYPE signed char

/* CCTK C/C++ Float datatypes */
#define CCTK_REAL16_TYPE long double
#define CCTK_REAL8_TYPE double
#define CCTK_REAL4_TYPE float

/* Disable 'restrict' for compiler versions known to be buggy */
/* We know that 20120731 fails and 20130728 passes. */
#if (defined __INTEL_COMPILER                &&                 \
     __INTEL_COMPILER_BUILD_DATE >= 20120731 &&                 \
     __INTEL_COMPILER_BUILD_DATE <  20130728 &&                 \
     !defined CCTK_INTEL_COMPILER_DONT_DISABLE_RESTRICT)
#  define CCTK_DISABLE_RESTRICT 1
#endif

/****************************************************************************/
/* C specific stuff */
/****************************************************************************/
#ifndef __cplusplus

/* Define to empty if the 'inline' keyword does not work. */
#define HAVE_CCTK_C_INLINE 1
#define HAVE_STANDARD_CCTK_C_INLINE 1
#define CCTK_C_INLINE inline

#ifdef HAVE_CCTK_C_INLINE
#  ifdef HAVE_STANDARD_CCTK_C_INLINE
/* 'inline' works -- do nothing (CCTK_C_INLINE is 'inline') */
#  else
/* need a macro to use 'inline' */
#    define inline CCTK_C_INLINE
#  endif
#else
/* 'inline' does not work (CCTK_C_INLINE is empty ) */
#  define inline
#endif

/* Define to 'static' if the 'static inline' keyword combination does
   not work. */
#define CCTK_C_STATIC_INLINE static inline
#define CCTK_STATIC_INLINE CCTK_C_STATIC_INLINE

/* Define to 'extern' if the 'extern inline' keyword combination does
   not work. */
/* #undef CCTK_C_EXTERN_INLINE */
#define CCTK_EXTERN_INLINE CCTK_C_EXTERN_INLINE

/* Define to empty if the 'const' keyword does not work. */
/* #undef const */

/* Define to empty if the 'restrict' keyword does not work. */
#ifdef CCTK_DISABLE_RESTRICT
#  define CCTK_C_RESTRICT
#else
/* #  undef HAVE_CCTK_C_RESTRICT */
/* #  undef CCTK_C_RESTRICT */
#endif

#ifdef CCTK_C_RESTRICT
#  define restrict CCTK_C_RESTRICT
#endif

/* Allow the use of CCTK_RESTRICT as a qualifier always. */
#ifdef CCTK_C_RESTRICT
#  define CCTK_RESTRICT CCTK_C_RESTRICT
#else
#  define CCTK_RESTRICT restrict
#endif

/* Disable _Pragma if unsupported */
#define HAVE_CCTK_C__PRAGMA 1
#ifndef HAVE_CCTK_C__PRAGMA
#  define _Pragma(x)
#endif

/* Whether copysign exists, and how it should be called */
#ifdef HAVE_COPYSIGN
#  define HAVE_CCTK_C_COPYSIGN HAVE_COPYSIGN
#  define CCTK_C_COPYSIGN copysign
#  define HAVE_CCTK_COPYSIGN HAVE_CCTK_C_COPYSIGN
#  define CCTK_COPYSIGN CCTK_C_COPYSIGN
#endif

/* Whether fpclassify exists, and how it should be called */
#ifdef HAVE_FPCLASSIFY
#  define HAVE_CCTK_C_FPCLASSIFY HAVE_FPCLASSIFY
#  define CCTK_C_FPCLASSIFY fpclassify
#  define HAVE_CCTK_FPCLASSIFY HAVE_CCTK_C_FPCLASSIFY
#  define CCTK_FPCLASSIFY CCTK_C_FPCLASSIFY
#endif

/* Whether isfinite exists, and how it should be called */
#ifdef HAVE_ISFINITE
#  define HAVE_CCTK_C_ISFINITE HAVE_ISFINITE
#  define CCTK_C_ISFINITE isfinite
#  define HAVE_CCTK_ISFINITE HAVE_CCTK_C_ISFINITE
#  define CCTK_ISFINITE CCTK_C_ISFINITE
#endif

/* Whether isinf exists, and how it should be called */
#ifdef HAVE_ISINF
#  define HAVE_CCTK_C_ISINF HAVE_ISINF
#  define CCTK_C_ISINF isinf
#  define HAVE_CCTK_ISINF HAVE_CCTK_C_ISINF
#  define CCTK_ISINF CCTK_C_ISINF
#endif

/* Whether isnan exists, and how it should be called */
#ifdef HAVE_ISNAN
#  define HAVE_CCTK_C_ISNAN HAVE_ISNAN
#  define CCTK_C_ISNAN isnan
#  define HAVE_CCTK_ISNAN HAVE_CCTK_C_ISNAN
#  define CCTK_ISNAN CCTK_C_ISNAN
#endif

/* Whether isnormal exists, and how it should be called */
#ifdef HAVE_ISNORMAL
#  define HAVE_CCTK_C_ISNORMAL HAVE_ISNORMAL
#  define CCTK_C_ISNORMAL isnormal
#  define HAVE_CCTK_ISNORMAL HAVE_CCTK_C_ISNORMAL
#  define CCTK_ISNORMAL CCTK_C_ISNORMAL
#endif

/* Whether signbit exists, and how it should be called */
#ifdef HAVE_SIGNBIT
#  define HAVE_CCTK_C_SIGNBIT HAVE_SIGNBIT
#  define CCTK_C_SIGNBIT signbit
#  define HAVE_CCTK_SIGNBIT HAVE_CCTK_C_SIGNBIT
#  define CCTK_SIGNBIT CCTK_C_SIGNBIT
#endif

/* Whether __attribute__((const)) exists. */
/* #undef HAVE_CCTK_C_ATTRIBUTE_CONST */
#ifdef HAVE_CCTK_C_ATTRIBUTE_CONST
#  define CCTK_ATTRIBUTE_CONST __attribute__((__const__))
#else
#  define CCTK_ATTRIBUTE_CONST
#endif

/* Whether __attribute__((pure)) exists. */
/* #undef HAVE_CCTK_C_ATTRIBUTE_PURE */
#ifdef HAVE_CCTK_C_ATTRIBUTE_PURE
#  define CCTK_ATTRIBUTE_PURE __attribute__((__pure__))
#else
#  define CCTK_ATTRIBUTE_PURE
#endif

/* Whether __attribute__((noinline)) exists. */
#define HAVE_CCTK_C_ATTRIBUTE_NOINLINE 1
#ifdef HAVE_CCTK_C_ATTRIBUTE_NOINLINE
#  define CCTK_ATTRIBUTE_NOINLINE __attribute__((__noinline__))
#else
#  define CCTK_ATTRIBUTE_NOINLINE
#endif

/* Whether __attribute__((always_inline)) exists. */
#define HAVE_CCTK_C_ATTRIBUTE_ALWAYS_INLINE 1
#ifdef HAVE_CCTK_C_ATTRIBUTE_ALWAYS_INLINE
#  define CCTK_ATTRIBUTE_ALWAYS_INLINE __attribute__((__always_inline__))
#else
#  define CCTK_ATTRIBUTE_ALWAYS_INLINE
#endif

/* Whether __attribute__((unused)) exists. */
#define HAVE_CCTK_C_ATTRIBUTE_UNUSED 1
#ifdef HAVE_CCTK_C_ATTRIBUTE_UNUSED
#  define CCTK_ATTRIBUTE_UNUSED __attribute__((__unused__))
#else
#  define CCTK_ATTRIBUTE_UNUSED
#endif

/* Whether __attribute__((aligned(...))) exists. */
#define HAVE_CCTK_C_ATTRIBUTE_ALIGNED 1
#ifdef HAVE_CCTK_C_ATTRIBUTE_ALIGNED
#  define CCTK_ATTRIBUTE_ALIGNED(x) __attribute__((__aligned__(x)))
#else
#  define CCTK_ATTRIBUTE_ALIGNED(x)
#endif

/* Whether __attribute__((cold)) exists. */
#define HAVE_CCTK_C_ATTRIBUTE_COLD 1
#ifdef HAVE_CCTK_C_ATTRIBUTE_COLD
#  define CCTK_ATTRIBUTE_COLD __attribute__((__cold__))
#else
#  define CCTK_ATTRIBUTE_COLD
#endif

/* Whether __attribute__((hot)) exists. */
#define HAVE_CCTK_C_ATTRIBUTE_HOT 1
#ifdef HAVE_CCTK_C_ATTRIBUTE_HOT
#  define CCTK_ATTRIBUTE_HOT __attribute__((__hot__))
#else
#  define CCTK_ATTRIBUTE_HOT
#endif

/* Whether __attribute__((format(...))) exists. */
#define HAVE_CCTK_C_ATTRIBUTE_FORMAT 1
#ifdef HAVE_CCTK_C_ATTRIBUTE_FORMAT
#  define CCTK_ATTRIBUTE_FORMAT(archetype, format, firstarg) __attribute__((__format__(archetype, format, firstarg)))
#else
#  define CCTK_ATTRIBUTE_FORMAT(archetype, format, firstarg)
#endif

/* Whether __attribute__((noreturn)) exists. */
#define HAVE_CCTK_C_ATTRIBUTE_NORETURN 1
#ifdef HAVE_CCTK_C_ATTRIBUTE_NORETURN
#  define CCTK_ATTRIBUTE_NORETURN __attribute__((__noreturn__))
#else
#  define CCTK_ATTRIBUTE_NORETURN
#endif

/* Whether __attribute__((nonnull)) exists. */
#define HAVE_CCTK_C_ATTRIBUTE_NONNULL 1
#ifdef HAVE_CCTK_C_ATTRIBUTE_NONNULL
#  define CCTK_ATTRIBUTE_NONNULL(m) __attribute__((__nonnull__(m)))
#else
#  define CCTK_ATTRIBUTE_NONNULL(m)
#endif

/* Whether __attribute__((returns_nonnull)) exists. */
#define HAVE_CCTK_C_ATTRIBUTE_RETURNS_NONNULL 1
#ifdef HAVE_CCTK_C_ATTRIBUTE_RETURNS_NONNULL
#  define CCTK_ATTRIBUTE_RETURNS_NONNULL __attribute__((__returns_nonnull__))
#else
#  define CCTK_ATTRIBUTE_RETURNS_NONNULL
#endif

/* Whether __builtin_expect exists. */
#define HAVE_CCTK_C_BUILTIN_EXPECT 1
#ifdef HAVE_CCTK_C_BUILTIN_EXPECT
#  define CCTK_BUILTIN_EXPECT __builtin_expect
#else
#  define CCTK_BUILTIN_EXPECT(x,y) (x)
#endif

/* Whether __builtin_unreachable exists. */
#define HAVE_CCTK_C_BUILTIN_UNREACHABLE 1
#ifdef HAVE_CCTK_C_BUILTIN_UNREACHABLE
#  define CCTK_BUILTIN_UNREACHABLE __builtin_unreachable
#else
#  define CCTK_BUILTIN_UNREACHABLE() CCTK_Abort(0, 1)
#endif

/* Whether __builtin_assume_aligned exists. */
#define HAVE_CCTK_C_BUILTIN_ASSUME_ALIGNED 1
#ifdef HAVE_CCTK_C_BUILTIN_ASSUME_ALIGNED
#  define CCTK_BUILTIN_ASSUME_ALIGNED __builtin_assume_aligned
#else
#  define CCTK_BUILTIN_ASSUME_ALIGNED(expr, ...) (expr)
#endif

/* OpenMP collapse clause */
#if (defined CCTK_DISABLE_OMP_COLLAPSE ||                               \
     (defined __IBMC__ && defined _ARCH_450D) ||                        \
     (defined __INTEL_COMPILER && __INTEL_COMPILER_BUILD_DATE < 20100801))
/* see http://software.intel.com/en-us/articles/intel-professional-edition-compilers-111-fixes-list/ */
#  define collapse(N)
#  ifndef CCTK_DISABLE_OMP_COLLAPSE
#    error "OpenMP collapse directive disabled for C, but enabled for Fortran -- likely an error in the option list. Please add -DCCTK_DISABLE_OMP_COLLAPSE to FPPFLAGS and CPPFLAGS."
#  endif
#else
/* #  undef collapse */
#endif

#endif /* ! defined __cplusplus */
/****************************************************************************/

/****************************************************************************/
/* C++ specific stuff */
/****************************************************************************/
#ifdef __cplusplus

#define CCTK_STATIC_INLINE static inline
#define CCTK_EXTERN_INLINE extern

/* Whether copysign exists, and how it should be called */
#define HAVE_CCTK_CXX_COPYSIGN 1
#define CCTK_CXX_COPYSIGN std::copysign
#ifdef HAVE_CCTK_CXX_COPYSIGN
#  define HAVE_CCTK_COPYSIGN HAVE_CCTK_CXX_COPYSIGN
#  define CCTK_COPYSIGN CCTK_CXX_COPYSIGN
#endif

/* Whether fpclassify exists, and how it should be called */
#define HAVE_CCTK_CXX_FPCLASSIFY 1
#define CCTK_CXX_FPCLASSIFY std::fpclassify
#ifdef HAVE_CCTK_CXX_FPCLASSIFY
#  define HAVE_CCTK_FPCLASSIFY HAVE_CCTK_CXX_FPCLASSIFY
#  define CCTK_FPCLASSIFY CCTK_CXX_FPCLASSIFY
#endif

/* Whether isinf exists, and how it should be called */
#define HAVE_CCTK_CXX_ISINF 1
#define CCTK_CXX_ISINF std::isinf
#ifdef HAVE_CCTK_CXX_ISINF
#  define HAVE_CCTK_ISINF HAVE_CCTK_CXX_ISINF
#  define CCTK_ISINF CCTK_CXX_ISINF
#endif

/* Whether isinf exists, and how it should be called */
#define HAVE_CCTK_CXX_ISINF 1
#define CCTK_CXX_ISINF std::isinf
#ifdef HAVE_CCTK_CXX_ISINF
#  define HAVE_CCTK_ISINF HAVE_CCTK_CXX_ISINF
#  define CCTK_ISINF CCTK_CXX_ISINF
#endif

/* Whether isnan exists, and how it should be called */
#define HAVE_CCTK_CXX_ISNAN 1
#define CCTK_CXX_ISNAN std::isnan
#ifdef HAVE_CCTK_CXX_ISNAN
#  define HAVE_CCTK_ISNAN HAVE_CCTK_CXX_ISNAN
#  define CCTK_ISNAN CCTK_CXX_ISNAN
#endif

/* Whether isnormal exists, and how it should be called */
#define HAVE_CCTK_CXX_ISNORMAL 1
#define CCTK_CXX_ISNORMAL std::isnormal
#ifdef HAVE_CCTK_CXX_ISNORMAL
#  define HAVE_CCTK_ISNORMAL HAVE_CCTK_CXX_ISNORMAL
#  define CCTK_ISNORMAL CCTK_CXX_ISNORMAL
#endif

/* Whether signbit exists, and how it should be called */
#define HAVE_CCTK_CXX_SIGNBIT 1
#define CCTK_CXX_SIGNBIT std::signbit
#ifdef HAVE_CCTK_CXX_SIGNBIT
#  define HAVE_CCTK_SIGNBIT HAVE_CCTK_CXX_SIGNBIT
#  define CCTK_SIGNBIT CCTK_CXX_SIGNBIT
#endif

/* Whether __attribute__((const)) exists. */
/* #undef HAVE_CCTK_CXX_ATTRIBUTE_CONST */
#ifdef HAVE_CCTK_CXX_ATTRIBUTE_CONST
#  define CCTK_ATTRIBUTE_CONST __attribute__((__const__))
#else
#  define CCTK_ATTRIBUTE_CONST
#endif
/* #undef HAVE_CCTK_CXX_MEMBER_ATTRIBUTE_CONST */
#ifdef HAVE_CCTK_CXX_MEMBER_ATTRIBUTE_CONST
#  define CCTK_MEMBER_ATTRIBUTE_CONST __attribute__((__const__))
#else
#  define CCTK_MEMBER_ATTRIBUTE_CONST
#endif

/* Whether __attribute__((pure)) exists. */
/* #undef HAVE_CCTK_CXX_ATTRIBUTE_PURE */
#ifdef HAVE_CCTK_CXX_ATTRIBUTE_PURE
#  define CCTK_ATTRIBUTE_PURE __attribute__((__pure__))
#else
#  define CCTK_ATTRIBUTE_PURE
#endif
/* #undef HAVE_CCTK_CXX_MEMBER_ATTRIBUTE_PURE */
#ifdef HAVE_CCTK_CXX_MEMBER_ATTRIBUTE_PURE
#  define CCTK_MEMBER_ATTRIBUTE_PURE __attribute__((__pure__))
#else
#  define CCTK_MEMBER_ATTRIBUTE_PURE
#endif

/* Whether __attribute__((noinline)) exists. */
#define HAVE_CCTK_CXX_ATTRIBUTE_NOINLINE 1
#ifdef HAVE_CCTK_CXX_ATTRIBUTE_NOINLINE
#  define CCTK_ATTRIBUTE_NOINLINE __attribute__((__noinline__))
#else
#  define CCTK_ATTRIBUTE_NOINLINE
#endif
#define HAVE_CCTK_CXX_MEMBER_ATTRIBUTE_NOINLINE 1
#ifdef HAVE_CCTK_CXX_MEMBER_ATTRIBUTE_NOINLINE
#  define CCTK_MEMBER_ATTRIBUTE_NOINLINE __attribute__((__noinline__))
#else
#  define CCTK_MEMBER_ATTRIBUTE_NOINLINE
#endif

/* Whether __attribute__((always_inline)) exists. */
#define HAVE_CCTK_CXX_ATTRIBUTE_ALWAYS_INLINE 1
#ifdef HAVE_CCTK_CXX_ATTRIBUTE_ALWAYS_INLINE
#  define CCTK_ATTRIBUTE_ALWAYS_INLINE __attribute__((__always_inline__))
#else
#  define CCTK_ATTRIBUTE_ALWAYS_INLINE
#endif
#define HAVE_CCTK_CXX_MEMBER_ATTRIBUTE_ALWAYS_INLINE 1
#ifdef HAVE_CCTK_CXX_MEMBER_ATTRIBUTE_ALWAYS_INLINE
#  define CCTK_MEMBER_ATTRIBUTE_ALWAYS_INLINE __attribute__((__always_inline__))
#else
#  define CCTK_MEMBER_ATTRIBUTE_ALWAYS_INLINE
#endif

/* Whether __attribute__((unused)) exists. */
#define HAVE_CCTK_CXX_ATTRIBUTE_UNUSED 1
#ifdef HAVE_CCTK_CXX_ATTRIBUTE_UNUSED
#  define CCTK_ATTRIBUTE_UNUSED __attribute__((__unused__))
#else
#  define CCTK_ATTRIBUTE_UNUSED
#endif

/* Whether __attribute__((aligned(...))) exists. */
#define HAVE_CCTK_CXX_ATTRIBUTE_ALIGNED 1
#ifdef HAVE_CCTK_CXX_ATTRIBUTE_ALIGNED
#  define CCTK_ATTRIBUTE_ALIGNED(x) __attribute__((__aligned__(x)))
#else
#  define CCTK_ATTRIBUTE_ALIGNED(x)
#endif

/* Whether __attribute__((cold)) exists. */
#define HAVE_CCTK_CXX_ATTRIBUTE_COLD 1
#ifdef HAVE_CCTK_CXX_ATTRIBUTE_COLD
#  define CCTK_ATTRIBUTE_COLD __attribute__((__cold__))
#else
#  define CCTK_ATTRIBUTE_COLD
#endif

/* Whether __attribute__((hot)) exists. */
#define HAVE_CCTK_CXX_ATTRIBUTE_HOT 1
#ifdef HAVE_CCTK_CXX_ATTRIBUTE_HOT
#  define CCTK_ATTRIBUTE_HOT __attribute__((__hot__))
#else
#  define CCTK_ATTRIBUTE_HOT
#endif

/* Whether __attribute__((format(...))) exists. */
#define HAVE_CCTK_CXX_ATTRIBUTE_FORMAT 1
#ifdef HAVE_CCTK_CXX_ATTRIBUTE_FORMAT
#  define CCTK_ATTRIBUTE_FORMAT(archetype, format, firstarg) __attribute__((__format__(archetype, format, firstarg)))
#else
#  define CCTK_ATTRIBUTE_FORMAT(archetype, format, firstarg)
#endif

/* Whether __attribute__((noreturn)) exists. */
#define HAVE_CCTK_CXX_ATTRIBUTE_NORETURN 1
#ifdef HAVE_CCTK_CXX_ATTRIBUTE_NORETURN
#  define CCTK_ATTRIBUTE_NORETURN __attribute__((__noreturn__))
#else
#  define CCTK_ATTRIBUTE_NORETURN
#endif

/* Whether __attribute__((nonnull)) exists. */
#define HAVE_CCTK_CXX_ATTRIBUTE_NONNULL 1
#ifdef HAVE_CCTK_CXX_ATTRIBUTE_NONNULL
#  define CCTK_ATTRIBUTE_NONNULL(m) __attribute__((__nonnull__(m)))
#else
#  define CCTK_ATTRIBUTE_NONNULL(m)
#endif

/* Whether __attribute__((returns_nonnull)) exists. */
#define HAVE_CCTK_CXX_ATTRIBUTE_RETURNS_NONNULL 1
#ifdef HAVE_CCTK_CXX_ATTRIBUTE_RETURNS_NONNULL
#  define CCTK_ATTRIBUTE_RETURNS_NONNULL __attribute__((__returns_nonnull__))
#else
#  define CCTK_ATTRIBUTE_RETURNS_NONNULL
#endif

/* Whether __builtin_expect exists. */
#define HAVE_CCTK_CXX_BUILTIN_EXPECT 1
#ifdef HAVE_CCTK_CXX_BUILTIN_EXPECT
#  define CCTK_BUILTIN_EXPECT __builtin_expect
#else
#  define CCTK_BUILTIN_EXPECT(x,y) (x)
#endif

/* Whether __builtin_unreachable exists. */
#define HAVE_CCTK_CXX_BUILTIN_UNREACHABLE 1
#ifdef __CUDACC__
#  define CCTK_BUILTIN_UNREACHABLE() asm("trap;")
#else
#  ifdef HAVE_CCTK_CXX_BUILTIN_UNREACHABLE
#    define CCTK_BUILTIN_UNREACHABLE __builtin_unreachable
#  else
#    define CCTK_BUILTIN_UNREACHABLE() CCTK_Abort(0, 1)
#  endif
#endif

/* Whether __builtin_assume_aligned exists. */
#define HAVE_CCTK_CXX_BUILTIN_ASSUME_ALIGNED 1
#ifdef HAVE_CCTK_CXX_BUILTIN_ASSUME_ALIGNED
#  define CCTK_BUILTIN_ASSUME_ALIGNED __builtin_assume_aligned
#else
#  define CCTK_BUILTIN_ASSUME_ALIGNED(expr, ...) (expr)
#endif

/* Whether static_assert exists. */
#define HAVE_CCTK_CXX_STATIC_ASSERT 1
#ifdef HAVE_CCTK_CXX_STATIC_ASSERT
#  define CCTK_STATIC_ASSERT(cond,msg) static_assert(cond, msg)
#else
#  define CCTK_STATIC_ASSERT_NAME1(x, y) x##y
#  define CCTK_STATIC_ASSERT_NAME2(x, y) CCTK_STATIC_ASSERT_NAME1(x, y)
#  define CCTK_STATIC_ASSERT(cond, msg) typedef int CCTK_STATIC_ASSERT_NAME2(cctk_sa_, __LINE__)[(cond) ? 1 : -1] CCTK_ATTRIBUTE_UNUSED
#  define static_assert(cond, msg) CCTK_STATIC_ASSERT(cond, msg)
#endif

/* Whether C++11 is supported. */
#define HAVE_CCTK_CXX_AUTO_SPECIFIER 1
#define HAVE_CCTK_CXX_LAMBDA 1
#define HAVE_CCTK_CXX_RANGE_BASED_FOR 1

/* Some C++ compilers recognise the restrict keyword */
/* Define to empty if the keyword does not work. */
#ifdef CCTK_DISABLE_RESTRICT
#  define CCTK_CXX_RESTRICT
#else
#  define HAVE_CCTK_CXX_RESTRICT 1
#  define CCTK_CXX_RESTRICT __restrict__
#endif

#ifdef CCTK_CXX_RESTRICT
#  define restrict CCTK_CXX_RESTRICT
#endif

/* Allow the use of CCTK_RESTRICT as a qualifier always. */
#ifdef CCTK_CXX_RESTRICT
#  define CCTK_RESTRICT CCTK_CXX_RESTRICT
#else
#  define CCTK_RESTRICT restrict
#endif

/* OpenMP collapse clause */
#if (defined CCTK_DISABLE_OMP_COLLAPSE ||                               \
     (defined __IBMCPP__ && defined _ARCH_450D) ||                      \
     ( defined __INTEL_COMPILER && __INTEL_COMPILER_BUILD_DATE < 20100801))
/* see http://software.intel.com/en-us/articles/intel-professional-edition-compilers-111-fixes-list/ */
#  define collapse(N)
#else
/* #  undef collapse */
#endif

#endif /* __cplusplus */

/****************************************************************************/

#endif /* CCODE */

#ifdef FCODE

#define HAVE_CCTK_FORTRAN_REAL4 1
#define HAVE_CCTK_FORTRAN_REAL8 1
#define HAVE_CCTK_FORTRAN_REAL16 1

#define HAVE_CCTK_FORTRAN_COMPLEX8 1
#define HAVE_CCTK_FORTRAN_COMPLEX16 1
#define HAVE_CCTK_FORTRAN_COMPLEX32 1

#define CCTK_REAL16_KIND 16
#define CCTK_COMPLEX32_KIND 32

/* OpenMP collapse clause */
#ifdef CCTK_DISABLE_OMP_COLLAPSE
/* see http://software.intel.com/en-us/articles/intel-professional-edition-compilers-111-fixes-list/ */
#  define collapse(N)
#  define COLLAPSE(N)
#else
/* #  undef collapse */
/* #  undef COLLAPSE */
#endif

#endif /* FCODE */

/* Now include the code to pick an appropriate precison for reals and ints */
#include "cctk_Types.h"

/* Include any other stuff which is specific to this architecture */
#include "cctk_Archdefs.h"

/* Include any extra stuff from optional extra packages. */
#include "cctk_Extradefs.h"

#endif /* _CCTK_CONFIG_H_ */
