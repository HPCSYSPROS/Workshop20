dnl /*@@
dnl   @file      aclocal.m4
dnl   @date      Thu Oct 21 00:05:58 1999
dnl   @author    Tom Goodale
dnl   @desc
dnl   Local Cactus macros
dnl   @enddesc
dnl   @version $Header$
dnl @@*/


dnl  These are copies of the standard autoconf macros, except they
dnl  use AC_TRY_COMPILE rather than AC_TRY_CPP to check for headers.
dnl  This gets round the problem on cygwin where the gnu cpp finds
dnl  the gcc headers and not the ones for the actual compiler.

dnl CCTK_CHECK_HEADER(HEADER-FILE, [ADDITIONAL_CODE [, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
AC_DEFUN(CCTK_CHECK_HEADER,
[dnl Do the transliteration at runtime so arg 1 can be a shell variable.
cctk_safe=`echo "$1" | sed 'y%./+-%__p_%'`
AC_MSG_CHECKING([for $1])
AC_CACHE_VAL(cctk_cv_header_$cctk_safe,
[AC_TRY_COMPILE([$2
#include <$1>], [ ], eval "cctk_cv_header_$cctk_safe=yes",
  eval "cctk_cv_header_$cctk_safe=no")])dnl
if eval "test \"`echo '$cctk_cv_header_'$cctk_safe`\" = yes"; then
  AC_MSG_RESULT(yes)
  ifelse([$3], , :, [$3])
else
  AC_MSG_RESULT(no)
ifelse([$4], , , [$4
])dnl
fi
])

dnl CCTK_CHECK_HEADERS(HEADER-FILE... [, ADDITIONAL_CODE [, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
AC_DEFUN(CCTK_CHECK_HEADERS,
[for cctk_hdr in $1
do
CCTK_CHECK_HEADER($cctk_hdr,
[$2 ],
[changequote(, )dnl
  cctk_tr_hdr=HAVE_`echo $cctk_hdr | sed 'y%abcdefghijklmnopqrstuvwxyz./-%ABCDEFGHIJKLMNOPQRSTUVWXYZ___%'`
changequote([, ])dnl
  AC_DEFINE_UNQUOTED($cctk_tr_hdr) $3], [$4])dnl
done
])

dnl A simplified version of AC_CHECK_MEMBER from autoconf 2.50 and up
dnl CCTK_CHECK_MEMBER(AGGREGATE, MEMBER,
dnl                 [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND],
dnl                 [INCLUDES])
dnl ---------------------------------------------------------
dnl AGGREGATE, MEMBER is for instance `struct passwd pw_gecos', shell
dnl variables are not a valid argument.
AC_DEFUN([CCTK_CHECK_MEMBER],
[dnl Extract the aggregate name, and the member name
ac_member_var=`echo $1['_']$2 | sed 'y% %_%'`
AC_MSG_CHECKING([for $1.$2])
AC_CACHE_VAL(ac_cv_member_$ac_member_var,
[dnl
AC_TRY_COMPILE([$5],
[static $1 ac_aggr;
if (ac_aggr.$2)
return 0;],
		eval "ac_cv_member_$ac_member_var=yes",
		eval "ac_cv_member_$ac_member_var=no"dnl
)
if eval "test \"`echo '$''{'ac_cv_member_$ac_member_var'}'`\" = no"; then
AC_TRY_COMPILE([$5],
[static $1 ac_aggr;
if (sizeof ac_aggr.$2)
return 0;],
		eval "ac_cv_member_$ac_member_var=yes",
		eval "ac_cv_member_$ac_member_var=no"dnl
)
fi dnl
])dnl
if eval "test \"`echo '$ac_cv_member_'$ac_member_var`\" = yes"; then
  AC_MSG_RESULT(yes)
  ifelse([$3], ,
[ac_tr_member=HAVE_MEMBER`echo $1 | sed -e 'y% %_%' \
    -e 'y/abcdefghijklmnopqrstuvwxyz/ABCDEFGHIJKLMNOPQRSTUVWXYZ/'`
  AC_DEFINE_UNQUOTED($ac_tr_member)
], [$3])
else
  AC_MSG_RESULT(no)
ifelse([$4], , , [$4
])dnl
fi 
])


dnl A version of AC_TRY_COMPILER(TEST-PROGRAM, WORKING-VAR, CROSS-VAR) which,
dnl if the TEST-PROGRAM could not be executed, does not throw away the stderr
dnl but redirects it into the config.log logfile instead.
dnl
dnl Try to compile, link and execute TEST-PROGRAM.  Set WORKING-VAR to
dnl `yes' if the current compiler works, otherwise set it ti `no'.  Set
dnl CROSS-VAR to `yes' if the compiler and linker produce non-native
dnl executables, otherwise set it to `no'.  Before calling
dnl `AC_TRY_COMPILER()', call `AC_LANG_*' to set-up for the right
dnl language.
dnl
dnl CCTK_TRY_COMPILER(TEST-PROGRAM, WORKING-VAR, CROSS-VAR)
AC_DEFUN(CCTK_TRY_COMPILER,
[cat > conftest.$ac_ext << EOF
ifelse(AC_LANG, [FORTRAN77], ,
[
[#]line __oline__ "configure"
#include "confdefs.h"
])
[$1]
EOF
if AC_TRY_EVAL(ac_link) && test -s conftest${ac_exeext}; then
  [$2]=yes
  # If we can't run a trivial program, we are probably using a cross compiler.
  if (./conftest; exit) 2>&AC_FD_CC; then
    [$3]=no
  else
    [$3]=yes
  fi
else
  echo "configure: failed program was:" >&AC_FD_CC
  cat conftest.$ac_ext >&AC_FD_CC
  [$2]=no
fi
rm -fr conftest*])


dnl CCTK_TRY_LINK_2(INCLUDES, FUNCTION-BODY, OTHER-FUNCTION_BODY,
dnl                [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
AC_DEFUN(CCTK_TRY_LINK_2,
[cat > conftest.$ac_ext <<EOF
ifelse(AC_LANG, [FORTRAN77],
[
      program main
      call [$2]
      call main2
      end
],
[dnl This sometimes fails to find confdefs.h, for some reason.
dnl [#]line __oline__ "[$]0"
[#]line __oline__ "configure"
#include "confdefs.h"
[$1]
int main() {
[$2]
; return 0; }
])EOF
cat > conftest2.$ac_ext <<EOF
ifelse(AC_LANG, [FORTRAN77],
[
      subroutine main2
      call [$3]
      end
],
[dnl This sometimes fails to find confdefs.h, for some reason.
dnl [#]line __oline__ "[$]0"
[#]line __oline__ "configure"
#include "confdefs.h"
[$1]
int main2() {
[$3]
; return 0; }
])EOF
if AC_TRY_EVAL(ac_link conftest2.$ac_ext) && test -s conftest${ac_exeext}; then
  ifelse([$4], , :, [rm -rf conftest*
  $4])
else
  echo "configure: failed program was:" >&AC_FD_CC
  cat conftest.$ac_ext >&AC_FD_CC
ifelse([$5], , , [  rm -rf conftest*
  $5
])dnl
fi
rm -f conftest*])


dnl CCTK_FIND_NULLDEVICE
dnl Have to do it in this rather bizarre way
dnl as cygwin emulates /dev/null
AC_DEFUN(CCTK_FIND_NULLDEVICE,
[AC_MSG_CHECKING([for the null device])
AC_CACHE_VAL(cctk_cv_nulldevice,
[
if test -d /dev ; then
  eval "cctk_cv_nulldevice=/dev/null"
else
  cat > NUL <<EOF
test
EOF
  if eval "`cat NUL > /dev/null 2>/dev/null`" ; then
    eval "cctk_cv_nulldevice=NUL"
  fi
fi
])
if eval "test -n \"$cctk_cv_nulldevice\"" ; then
  AC_MSG_RESULT($cctk_cv_nulldevice)
  AC_DEFINE_UNQUOTED(NULL_DEVICE, "$cctk_cv_nulldevice")
else
  AC_MSG_RESULT("not found")
fi
])

AC_DEFUN(CCTK_TIME__FTIME,
[AC_MSG_CHECKING([for availability of _ftime timing])
AC_CACHE_VAL(cctk_cv_time_ftime,
[AC_TRY_LINK([#include <stdio.h>
#include <time.h>
#include <sys/types.h>
#include <sys/timeb.h>],
[  struct _timeb timebs;
  _ftime(&timebs);
  printf("%f\n",(double)(timebs.time + timebs.millitm/1000.0));
  return 0;], eval "cctk_cv_time_ftime=yes",
  eval "cctk_cv_time_ftime=no")])dnl
if eval "test \"`echo '$cctk_cv_time_ftime'`\" = yes"; then
  AC_MSG_RESULT(yes)
  AC_DEFINE_UNQUOTED(HAVE_TIME__FTIME)
else
  AC_MSG_RESULT(no)
fi
])dnl

AC_DEFUN(CCTK_TIME_GETRUSAGE,
[AC_MSG_CHECKING([for availability of getrusage timing])
AC_CACHE_VAL(cctk_cv_time_getrusage,
[AC_TRY_LINK([#include <stdio.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>],
[struct rusage ru;
 getrusage(RUSAGE_SELF, &ru);
 printf("%f\n",(double)(ru.ru_utime.tv_sec + (double)ru.ru_utime.tv_usec/1000000.0));
 return 0;], eval "cctk_cv_time_getrusage=yes",
  eval "cctk_cv_time_getrusage=no")])dnl
if eval "test \"`echo '$cctk_cv_time_getrusage'`\" = yes"; then
  AC_MSG_RESULT(yes)
  AC_DEFINE_UNQUOTED(HAVE_TIME_GETRUSAGE)
else
  AC_MSG_RESULT(no)
fi
])dnl

AC_DEFUN(CCTK_TIME_GETTIMEOFDAY,
[AC_MSG_CHECKING([for availability of gettimeofday timing])
AC_CACHE_VAL(cctk_cv_time_gettimeofday,
[AC_TRY_LINK([],
[gettimeofday(0, 0);
 return 0;], eval "cctk_cv_time_gettimeofday=yes",
  eval "cctk_cv_time_gettimeofday=no")])dnl
if eval "test \"`echo '$cctk_cv_time_gettimeofday'`\" = yes"; then
  AC_MSG_RESULT(yes)
  AC_DEFINE_UNQUOTED(HAVE_TIME_GETTIMEOFDAY)
else
  AC_MSG_RESULT(no)
fi
]dnl
if eval "test \"`echo '$cctk_cv_time_gettimeofday'`\" = yes"; then
[AC_MSG_CHECKING([if gettimeofday needs timezone])
AC_CACHE_VAL(cctk_cv_time_gettimeofday_timezone,
[AC_TRY_LINK([#include <stdio.h>
#include <sys/time.h>
#include <unistd.h>],
[struct timeval tp;
 struct timezone tzp;
 gettimeofday(&tp, &tzp);
 printf("%f\n", (double)(tp.tv_sec + (double)tp.tv_usec/1000000.0));
 return 0;], eval "cctk_cv_time_gettimeofday_timezone=yes",
  eval "cctk_cv_time_gettimeofday_timezone=no")])dnl
if eval "test \"`echo '$cctk_cv_time_gettimeofday_timezone'`\" = yes"; then
  AC_MSG_RESULT(yes)
  AC_DEFINE_UNQUOTED(GETTIMEOFDAY_NEEDS_TIMEZONE)
else
  AC_MSG_RESULT(no)
fi
]dnl
fi
)dnl

AC_DEFUN(CCTK_PROG_CC_WORKS,
[AC_MSG_CHECKING([whether the C compiler ($CC $CFLAGS $LDFLAGS) works])
AC_LANG_SAVE
AC_LANG_C
rm -fr conftest*
CCTK_TRY_COMPILER([main(){return(0);} int PilotMain(){return(0);}], ac_cv_prog_cc_works, ac_cv_prog_cc_cross)
AC_LANG_RESTORE
AC_MSG_RESULT($ac_cv_prog_cc_works)
if test $ac_cv_prog_cc_works = no; then
  AC_MSG_ERROR([installation or configuration problem: C compiler cannot create executables (see configs/<configname>/config-data/config.log for details).])
fi
changequote({, })
CROSS_COMPILE=`echo $CROSS_COMPILE | tr '[:upper:]' '[:lower:]'`
changequote([, ])
AC_MSG_CHECKING([whether the C compiler ($CC $CFLAGS $LDFLAGS) is a cross-compiler])
AC_MSG_RESULT($ac_cv_prog_cc_cross)
if test $ac_cv_prog_cc_cross = yes -a "x$CROSS_COMPILE" != xyes; then
  AC_MSG_ERROR([Could not run executable generated by C compiler (see configs/<configname>/config-data/config.log for details). If this is a cross configuration please set CROSS_COMPILE=yes.])
fi
cross_compiling=$ac_cv_prog_cc_cross
])

AC_DEFUN(CCTK_PROG_CXX_WORKS,
[AC_MSG_CHECKING([whether the C++ compiler ($CXX $CXXFLAGS $LDFLAGS) works])
AC_LANG_SAVE
AC_LANG_CPLUSPLUS
rm -fr conftest*
CCTK_TRY_COMPILER([int main(){return(0);} extern "C" int PilotMain(){return(0);}], ac_cv_prog_cxx_works, ac_cv_prog_cxx_cross)
AC_LANG_RESTORE
AC_MSG_RESULT($ac_cv_prog_cxx_works)
if test $ac_cv_prog_cxx_works = no; then
  AC_MSG_ERROR([installation or configuration problem: C++ compiler cannot create executables (see configs/<configname>/config-data/config.log for details).])
fi
AC_MSG_CHECKING([whether the C++ compiler ($CXX $CXXFLAGS $LDFLAGS) is a cross-compiler])
AC_MSG_RESULT($ac_cv_prog_cxx_cross)
if test $ac_cv_prog_cxx_cross = yes -a "x$CROSS_COMPILE" != xyes; then
  AC_MSG_ERROR([Could not run executable generated by C++ compiler (see configs/<configname>/config-data/config.log for details). If this is a cross configuration please set CROSS_COMPILE=yes.])
fi
cross_compiling=$ac_cv_prog_cxx_cross
])

AC_DEFUN(CCTK_HEADER_REGEX,
[AC_MSG_CHECKING([for regex.h])
AC_CACHE_VAL(cctk_cv_header_regex_h,
[AC_TRY_COMPILE([#include <stdio.h>
#include <regex.h>], [return 0;], eval "cctk_cv_header_regex_h=yes",
  eval "cctk_cv_header_regex_h=no")])dnl
if eval "test \"`echo '$cctk_cv_header_regex_h'`\" = yes"; then
  AC_MSG_RESULT(yes)
  AC_DEFINE_UNQUOTED(HAVE_REGEX_H)
else
  AC_MSG_RESULT(no)
fi
])

# CCTK_CHECK_FUNCS(FUNCTION..., [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
# ---------------------------------------------------------------------
AC_DEFUN([CCTK_CHECK_FUNCS],
[ac_link='${CC-cc} -o conftest$ac_exeext $CFLAGS $CPPFLAGS $LDFLAGS conftest.$ac_ext `CCTK_Wrap "$LIBDIR_PREFIX" "$LIBDIR_SUFFIX" "$LIBDIRS"` `CCTK_Wrap "$LIBLINK_PREFIX" "$LIBLINK_SUFFIX" "$LIBS"` >&5'
dnl AC_CHECK_FUNCS does not properly quote its last argument
AC_CHECK_FUNCS([$1],[$2],[[$3]])
])

# CCTK_CHECK_FUNC(FUNCTION, [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
# ---------------------------------------------------------------------
AC_DEFUN([CCTK_CHECK_FUNC],
[ac_link='${CC-cc} -o conftest$ac_exeext $CFLAGS $CPPFLAGS $LDFLAGS conftest.$ac_ext `CCTK_Wrap "$LIBDIR_PREFIX" "$LIBDIR_SUFFIX" "$LIBDIRS"` `CCTK_Wrap "$LIBLINK_PREFIX" "$LIBLINK_SUFFIX" "$LIBS"` >&5'
AC_CHECK_FUNC([$1],[$2],[$3])
])


# CCTK_CHECK_LIB(LIBRARY, FUNCTION,
#              [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND],
#              [OTHER-LIBRARIES])
# ------------------------------------------------------
AC_DEFUN(CCTK_CHECK_LIB,
[AC_MSG_CHECKING([for $2 in library $1])
dnl Use a cache variable name containing both the library and function name,
dnl because the test really is for library $1 defining function $2, not
dnl just for library $1.  Separate tests with the same $1 and different $2s
dnl may have different results.
ac_lib_var=`echo $1['_']$2 | sed 'y%./+-%__p_%'`
AC_CACHE_VAL(ac_cv_lib_$ac_lib_var,
[ac_link='${CC-cc} -o conftest$ac_exeext $CFLAGS $CPPFLAGS $LDFLAGS conftest.$ac_ext `CCTK_Wrap "$LIBDIR_PREFIX" "$LIBDIR_SUFFIX" "$LIBDIRS"` `CCTK_Wrap "$LIBLINK_PREFIX" "$LIBLINK_SUFFIX" "$LIBS"` >&5'
ac_save_LIBS="$LIBS"
LIBS="$1 $5 $LIBS"
AC_TRY_LINK(dnl
ifelse(AC_LANG, [FORTRAN77], ,
ifelse([$2], [main], , dnl Avoid conflicting decl of main.
[/* Override any gcc2 internal prototype to avoid an error.  */
]ifelse(AC_LANG, CPLUSPLUS, [#ifdef __cplusplus
extern "C"
#endif
])dnl
[/* We use char because int might match the return type of a gcc2
    builtin and then its argument prototype would still apply.  */
dnl We need a space before the parentheses to avoid accidentally calling
dnl an m4 function, if $2 expands to an m4 function.
char $2 ();
])),
            [$2()],
            eval "ac_cv_lib_$ac_lib_var=yes",
            eval "ac_cv_lib_$ac_lib_var=no")
LIBS="$ac_save_LIBS"
])dnl
if eval "test \"`echo '$ac_cv_lib_'$ac_lib_var`\" = yes"; then
  AC_MSG_RESULT(yes)
  ifelse([$3], ,
[changequote(, )dnl
  ac_tr_lib=HAVE_LIB`echo $1 | sed -e 's/[^a-zA-Z0-9_]/_/g' \
    -e 'y/abcdefghijklmnopqrstuvwxyz/ABCDEFGHIJKLMNOPQRSTUVWXYZ/'`
changequote([, ])dnl
  AC_DEFINE_UNQUOTED($ac_tr_lib)
  LIBS="$1 $LIBS"
], [$3])
else
  AC_MSG_RESULT(no)
ifelse([$4], , , [$4
])dnl
fi
])

AC_DEFUN(CCTK_CHECK_LIB_FUNC,
[CCTK_CHECK_LIB([$1], [$2],
[ifelse([$3], , [changequote(, )dnl
  cctk_tr_lib=HAVE_LIB`echo $1 | sed -e 's/[^a-zA-Z0-9_]/_/g' \
    -e 'y/abcdefghijklmnopqrstuvwxyz/ABCDEFGHIJKLMNOPQRSTUVWXYZ/'`
  cctk_tr_func=HAVE_`echo $2 | sed -e 's/[^a-zA-Z0-9_]/_/g' \
    -e 'y/abcdefghijklmnopqrstuvwxyz/ABCDEFGHIJKLMNOPQRSTUVWXYZ/'`
changequote([, ])dnl
  AC_DEFINE_UNQUOTED($cctk_tr_lib)
  AC_DEFINE_UNQUOTED($cctk_tr_func)
  LIBS="$1 $LIBS"
], [$3])],dnl
[ifelse([$4], , , [$4
])])dnl
])



# CCTK_CHECK_HEADER_LIB(HEADER, LIBRARY, FUNCTION, ARGUMENTS,
#                       [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND],
#                       [OTHER-LIBRARIES])
# ------------------------------------------------------
AC_DEFUN(CCTK_CHECK_HEADER_LIB,
[AC_MSG_CHECKING([for $3 in header $1 and library $2])
dnl Use a cache variable name containing the header, library, and function
dnl name, because the test really is for header $1 and library $2 defining
dnl function $3, not just for header $1 and library $2.  Separate tests with
dnl the same $1 or $2 and different $3s may have different results.
ac_lib_var=`echo $1['_']$2['_']$3 | sed 'y%./+-%__p_%'`
AC_CACHE_VAL(ac_cv_lib_$ac_lib_var,
[ac_link='${CC-cc} -o conftest$ac_exeext $CFLAGS $CPPFLAGS $LDFLAGS conftest.$ac_ext `CCTK_Wrap "$LIBDIR_PREFIX" "$LIBDIR_SUFFIX" "$LIBDIRS"` `CCTK_Wrap "$LIBLINK_PREFIX" "$LIBLINK_SUFFIX" "$LIBS"` >&5'
ac_save_LIBS="$LIBS"
LIBS="$2 $7 $LIBS"
AC_TRY_LINK(dnl
ifelse(AC_LANG, [FORTRAN77], ,
ifelse([$3], [main], , dnl Avoid conflicting decl of main.
[]ifelse(AC_LANG, CPLUSPLUS, [#ifdef __cplusplus
extern "C"
#endif
])dnl
[#include <$1>
])),
            [$3 $4],
            eval "ac_cv_lib_$ac_lib_var=yes",
            eval "ac_cv_lib_$ac_lib_var=no")
LIBS="$ac_save_LIBS"
])dnl
if eval "test \"`echo '$ac_cv_lib_'$ac_lib_var`\" = yes"; then
  AC_MSG_RESULT(yes)
  ifelse([$5], ,
[changequote(, )dnl
  ac_tr_lib=HAVE_LIB`echo $2 | sed -e 's/[^a-zA-Z0-9_]/_/g' \
    -e 'y/abcdefghijklmnopqrstuvwxyz/ABCDEFGHIJKLMNOPQRSTUVWXYZ/'`
changequote([, ])dnl
  AC_DEFINE_UNQUOTED($ac_tr_lib)
  LIBS="$2 $LIBS"
], [$5])
else
  AC_MSG_RESULT(no)
ifelse([$6], , , [$6
])dnl
fi
])

AC_DEFUN(CCTK_CHECK_HEADER_LIB_FUNC,
[CCTK_CHECK_HEADER_LIB([$1], [$2], [$3], [$4],
[ifelse([$5], , [changequote(, )dnl
  cctk_tr_header=HAVE_`echo $1 | sed -e 's/[^a-zA-Z0-9_]/_/g' \
    -e 'y/abcdefghijklmnopqrstuvwxyz/ABCDEFGHIJKLMNOPQRSTUVWXYZ/'`
  cctk_tr_lib=HAVE_LIB`echo $2 | sed -e 's/[^a-zA-Z0-9_]/_/g' \
    -e 'y/abcdefghijklmnopqrstuvwxyz/ABCDEFGHIJKLMNOPQRSTUVWXYZ/'`
  cctk_tr_func=HAVE_`echo $3 | sed -e 's/[^a-zA-Z0-9_]/_/g' \
    -e 'y/abcdefghijklmnopqrstuvwxyz/ABCDEFGHIJKLMNOPQRSTUVWXYZ/'`
changequote([, ])dnl
  AC_DEFINE_UNQUOTED($cctk_tr_header)
  AC_DEFINE_UNQUOTED($cctk_tr_lib)
  AC_DEFINE_UNQUOTED($cctk_tr_func)
  LIBS="$2 $LIBS"
], [$5])]dnl
)
ifelse([$6], , , [$6
])dnl
])



dnl Do nothing if the compiler accepts the restrict keyword.
dnl Otherwise define restrict to __restrict__ or __restrict if one of
dnl those work, otherwise define restrict to be empty.
AC_DEFUN(CCTK_CHECK_C_RESTRICT,
[AC_CACHE_CHECK([for C restrict], cctk_cv_c_restrict,
[cctk_cv_c_restrict=no
for ac_kw in restrict __restrict__ __restrict; do
  AC_TRY_COMPILE([
double * $ac_kw p1;
double * $ac_kw p2[3];
struct s1 { char * $ac_kw arr; };
struct s2 { char * $ac_kw arr[3]; };
void f1 (void * $ac_kw p);
void f2 (void * $ac_kw p[]);
void f3 (void * $ac_kw p[3]);
/* void f4 (void * $ac_kw p[$ac_kw]); */
/* void f5 (void * $ac_kw p[$ac_kw 3]); */
void f1 (void * $ac_kw p) { }
void f2 (void * $ac_kw p[]) { }
void f3 (void * $ac_kw p[3]) { }
/* void f4 (void * $ac_kw p[$ac_kw]) { } */
/* void f5 (void * $ac_kw p[$ac_kw 3]) { } */
], [
double * $ac_kw v1;
double * $ac_kw v2[3];
], [cctk_cv_c_restrict=$ac_kw; break])
done
])
case "$cctk_cv_c_restrict" in
  restrict | yes) ;;
  no) AC_DEFINE(CCTK_C_RESTRICT, ) ;;
  *)  AC_DEFINE_UNQUOTED(CCTK_C_RESTRICT, $cctk_cv_c_restrict)
      AC_DEFINE(HAVE_CCTK_C_RESTRICT, 1) ;;
esac
])

AC_DEFUN(CCTK_CHECK_CXX_RESTRICT,
[AC_CACHE_CHECK([for C++ restrict], cctk_cv_cxx_restrict,
[cctk_cv_cxx_restrict=no
AC_LANG_SAVE
AC_LANG_CPLUSPLUS
for ac_kw in restrict __restrict__ __restrict; do
  AC_TRY_COMPILE([
double * $ac_kw p1;
double * $ac_kw p2[3];
struct s1 { char * $ac_kw arr; };
struct s2 { char * $ac_kw arr[3]; };
void f1 (void * $ac_kw p);
void f2 (void * $ac_kw p[]);
void f3 (void * $ac_kw p[3]);
// void f4 (void * $ac_kw p[$ac_kw]);
// void f5 (void * $ac_kw p[$ac_kw 3]);
void f1 (void * $ac_kw p) { }
void f2 (void * $ac_kw p[]) { }
void f3 (void * $ac_kw p[3]) { }
// void f4 (void * $ac_kw p[$ac_kw]) { }
// void f5 (void * $ac_kw p[$ac_kw 3]) { }
], [
double * $ac_kw v1;
double * $ac_kw v2[3];
], [cctk_cv_cxx_restrict=$ac_kw; break])
done
AC_LANG_RESTORE
])
case "$cctk_cv_cxx_restrict" in
  restrict | yes) ;;
  no) AC_DEFINE(CCTK_CXX_RESTRICT, ) ;;
  *)  AC_DEFINE_UNQUOTED(CCTK_CXX_RESTRICT, $cctk_cv_cxx_restrict)
      AC_DEFINE(HAVE_CCTK_CXX_RESTRICT, 1) ;;
esac
])

dnl Do nothing if the compiler accepts the inline keyword.  Otherwise
dnl define inline to __inline__ or __inline if one of those work,
dnl otherwise define inline to be empty.
AC_DEFUN(CCTK_CHECK_C_INLINE,
[AC_CACHE_CHECK([for C inline], cctk_cv_c_inline,
[cctk_cv_c_inline=no
for ac_kw in inline __inline__ __inline; do
  AC_TRY_COMPILE(, [} $ac_kw int foo() {], [cctk_cv_c_inline=$ac_kw; break])
done
])
case "$cctk_cv_c_inline" in
  inline) 
      AC_DEFINE_UNQUOTED(CCTK_C_INLINE, inline)
      AC_DEFINE(HAVE_CCTK_C_INLINE, 1)
      AC_DEFINE(HAVE_STANDARD_CCTK_C_INLINE, 1) ;;
  no) AC_DEFINE_UNQUOTED(CCTK_C_INLINE, )
      AC_DEFINE(CCTK_C_INLINE, ) ;;
  *)  AC_DEFINE_UNQUOTED(CCTK_C_INLINE, $cctk_cv_c_inline)
      AC_DEFINE(HAVE_CCTK_C_INLINE, 1) ;;
esac
])

dnl Define the macro EXTERN_INLINE to the keywords that the compiler needs
dnl to obtain what is obtained by "extern inline" in the C99 standard.  If
dnl the compiler does not support the "inline" keyword, define the macro
dnl to "extern".
AC_DEFUN(CCTK_CHECK_C_EXTERN_INLINE,
[AC_CACHE_CHECK([for C extern inline], cctk_cv_c_extern_inline,
[cctk_cv_c_extern_inline=no
for ac_kw in 'extern inline' 'extern __inline__' 'extern __inline' extern; do
  CCTK_TRY_LINK_2(, [foo();], [;} $ac_kw foo() {], [cctk_cv_c_extern_inline=$ac_kw; break])
done
])
case "$cctk_cv_c_extern_inline" in
  no) AC_DEFINE(CCTK_C_EXTERN_INLINE, extern) ;;
  *)  AC_DEFINE_UNQUOTED(CCTK_C_EXTERN_INLINE, $cctk_cv_c_extern_inline) ;;
esac
])

dnl Define the macro STATIC_INLINE to the keywords that the compiler needs
dnl to obtain what is obtained by "static inline" in the C99 standard.  If
dnl the compiler does not support the "inline" keyword, define the macro
dnl to "static".
AC_DEFUN(CCTK_CHECK_C_STATIC_INLINE,
[AC_CACHE_CHECK([for C static inline], cctk_cv_c_static_inline,
[cctk_cv_c_static_inline=no
for ac_kw in 'static inline' 'static __inline__' 'static __inline' static; do
  CCTK_TRY_LINK_2(, [;} $ac_kw ifoo(){} foo(){ifoo();], [;} $ac_kw ifoo(){} foo2(){ifoo();], [cctk_cv_c_static_inline=$ac_kw; break])
done
])
case "$cctk_cv_c_static_inline" in
  no) AC_DEFINE(CCTK_C_STATIC_INLINE, static) ;;
  *)  AC_DEFINE_UNQUOTED(CCTK_C_STATIC_INLINE, $cctk_cv_c_static_inline) ;;
esac
])

dnl Do nothing if the compiler accepts the _Pragma keyword.
dnl Otherwise define _Pragma to be empty.
AC_DEFUN(CCTK_C__PRAGMA,
[AC_CACHE_CHECK([for C _Pragma], cctk_cv_have_c__Pragma,
[cctk_cv_have_c__Pragma=no
AC_TRY_COMPILE([
#define LOOP(i) _Pragma("omp for") for (int i=0; i<10; ++i)
],[
  int s=0;
#pragma omp parallel reduction(+: s)
  LOOP(i) {
    s+=i;
  }
], cctk_cv_have_c__Pragma=yes, cctk_cv_have_c__Pragma=no)
])
if test "$cctk_cv_have_c__Pragma" = "yes" ; then
   AC_DEFINE(HAVE_CCTK_C__PRAGMA)
fi
])

dnl The autoconf 2.13 function AC_TRY_COMPILE does not work for Fortran.
dnl This version is corrected and should work for both C and Fortran.
dnl CCTK_TRY_COMPILE(INCLUDES, FUNCTION-BODY,
dnl             [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
AC_DEFUN(CCTK_TRY_COMPILE,
[cat > conftest.$ac_ext <<EOF
ifelse(AC_LANG, [FORTRAN77],
[[$1]
      program main
[$2]
      end
],
[dnl This sometimes fails to find confdefs.h, for some reason.
dnl [#]line __oline__ "[$]0"
[#]line __oline__ "configure"
#include "confdefs.h"
[$1]
int main() {
[$2]
; return 0; }
])EOF
if AC_TRY_EVAL(ac_compile); then
  ifelse([$3], , :, [rm -rf conftest*
  $3])
else
  echo "configure: failed program was:" >&AC_FD_CC
  cat conftest.$ac_ext >&AC_FD_CC
ifelse([$4], , , [  rm -rf conftest*
  $4
])dnl
fi
rm -f conftest*])

AC_DEFUN(CCTK_FORTRAN_REAL4,
[AC_CACHE_CHECK([for Fortran REAL*4], cctk_cv_have_fortran_real4,
[cctk_cv_have_fortran_real4=no
AC_LANG_SAVE
AC_LANG_FORTRAN77
CCTK_TRY_COMPILE(,[      REAL*4 a], cctk_cv_have_fortran_real4=yes, cctk_cv_have_fortran_real4=no)
AC_LANG_RESTORE
])
if test "$cctk_cv_have_fortran_real4" = "yes" ; then
   AC_DEFINE(HAVE_CCTK_FORTRAN_REAL4)
fi
])

AC_DEFUN(CCTK_FORTRAN_REAL8,
[AC_CACHE_CHECK([for Fortran REAL*8], cctk_cv_have_fortran_real8,
[cctk_cv_have_fortran_real8=no
AC_LANG_SAVE
AC_LANG_FORTRAN77
CCTK_TRY_COMPILE(,[      REAL*8 a], cctk_cv_have_fortran_real8=yes, cctk_cv_have_fortran_real8=no)
AC_LANG_RESTORE
])
if test "$cctk_cv_have_fortran_real8" = "yes" ; then
   AC_DEFINE(HAVE_CCTK_FORTRAN_REAL8)
fi
])

AC_DEFUN(CCTK_FORTRAN_REAL16,
[AC_CACHE_CHECK([for Fortran REAL*16], cctk_cv_have_fortran_real16,
[cctk_cv_have_fortran_real16=no
AC_LANG_SAVE
AC_LANG_FORTRAN77
CCTK_TRY_COMPILE(,[      REAL*$REAL16_KIND a], cctk_cv_have_fortran_real16=yes, cctk_cv_have_fortran_real16=no)
AC_LANG_RESTORE
])
if test "$cctk_cv_have_fortran_real16" = "yes" ; then
   AC_DEFINE(HAVE_CCTK_FORTRAN_REAL16)
fi
])

AC_DEFUN(CCTK_FORTRAN_COMPLEX8,
[AC_CACHE_CHECK([for Fortran COMPLEX*8], cctk_cv_have_fortran_complex8,
[cctk_cv_have_fortran_complex8=no
AC_LANG_SAVE
AC_LANG_FORTRAN77
CCTK_TRY_COMPILE(,[      COMPLEX*8 a], cctk_cv_have_fortran_complex8=yes, cctk_cv_have_fortran_complex8=no)
AC_LANG_RESTORE
])
if test "$cctk_cv_have_fortran_complex8" = "yes" ; then
   AC_DEFINE(HAVE_CCTK_FORTRAN_COMPLEX8)
fi
])

AC_DEFUN(CCTK_FORTRAN_COMPLEX16,
[AC_CACHE_CHECK([for Fortran COMPLEX*16], cctk_cv_have_fortran_complex16,
[cctk_cv_have_fortran_complex16=no
AC_LANG_SAVE
AC_LANG_FORTRAN77
CCTK_TRY_COMPILE(,[      COMPLEX*16 a], cctk_cv_have_fortran_complex16=yes, cctk_cv_have_fortran_complex16=no)
AC_LANG_RESTORE
])
if test "$cctk_cv_have_fortran_complex16" = "yes" ; then
   AC_DEFINE(HAVE_CCTK_FORTRAN_COMPLEX16)
fi
])

AC_DEFUN(CCTK_FORTRAN_COMPLEX32,
[AC_CACHE_CHECK([for Fortran COMPLEX*32], cctk_cv_have_fortran_complex32,
[cctk_cv_have_fortran_complex32=no
AC_LANG_SAVE
AC_LANG_FORTRAN77
CCTK_TRY_COMPILE(,[      COMPLEX*$COMPLEX32_KIND a], cctk_cv_have_fortran_complex32=yes, cctk_cv_have_fortran_complex32=no)
AC_LANG_RESTORE
])
if test "$cctk_cv_have_fortran_complex32" = "yes" ; then
   AC_DEFINE(HAVE_CCTK_FORTRAN_COMPLEX32)
fi
])



AC_DEFUN(CCTK_C_ATTRIBUTE_CONST,
[AC_CACHE_CHECK([for C function __attribute__((__const__))], cctk_cv_have_c_attribute_const,
[cctk_cv_have_c_attribute_const=no
AC_TRY_COMPILE(, double foo (double) __attribute__((__const__));, cctk_cv_have_c_attribute_const=yes, cctk_cv_have_c_attribute_const=no)
])
if test "$cctk_cv_have_c_attribute_const" = "yes" ; then
   AC_DEFINE(HAVE_CCTK_C_ATTRIBUTE_CONST)
fi
])

AC_DEFUN(CCTK_CXX_ATTRIBUTE_CONST,
[AC_CACHE_CHECK([for C++ function __attribute__((__const__))], cctk_cv_have_cxx_attribute_const,
[cctk_cv_have_cxx_attribute_const=no
AC_LANG_SAVE
AC_LANG_CPLUSPLUS
AC_TRY_COMPILE(, double foo (double) __attribute__((__const__));, cctk_cv_have_cxx_attribute_const=yes, cctk_cv_have_cxx_attribute_const=no)
AC_LANG_RESTORE
])
if test "$cctk_cv_have_cxx_attribute_const" = "yes" ; then
   AC_DEFINE(HAVE_CCTK_CXX_ATTRIBUTE_CONST)
fi
])

AC_DEFUN(CCTK_CXX_MEMBER_ATTRIBUTE_CONST,
[AC_CACHE_CHECK([for C++ member function __attribute__((__const__))], cctk_cv_have_cxx_member_attribute_const,
[cctk_cv_have_cxx_member_attribute_const=no
AC_LANG_SAVE
AC_LANG_CPLUSPLUS
AC_TRY_COMPILE(, struct bar { double foo (double) __attribute__((__const__)); };, cctk_cv_have_cxx_member_attribute_const=yes, cctk_cv_have_cxx_member_attribute_const=no)
AC_LANG_RESTORE
])
if test "$cctk_cv_have_cxx_member_attribute_const" = "yes" ; then
   AC_DEFINE(HAVE_CCTK_CXX_MEMBER_ATTRIBUTE_CONST)
fi
])



AC_DEFUN(CCTK_C_ATTRIBUTE_PURE,
[AC_CACHE_CHECK([for C function __attribute__((__pure__))], cctk_cv_have_c_attribute_pure,
[cctk_cv_have_c_attribute_pure=no
AC_TRY_COMPILE(, double foo (double) __attribute__((__pure__));, cctk_cv_have_c_attribute_pure=yes, cctk_cv_have_c_attribute_pure=no)
])
if test "$cctk_cv_have_c_attribute_pure" = "yes" ; then
   AC_DEFINE(HAVE_CCTK_C_ATTRIBUTE_PURE)
fi
])

AC_DEFUN(CCTK_CXX_ATTRIBUTE_PURE,
[AC_CACHE_CHECK([for C++ function __attribute__((__pure__))], cctk_cv_have_cxx_attribute_pure,
[cctk_cv_have_cxx_attribute_pure=no
AC_LANG_SAVE
AC_LANG_CPLUSPLUS
AC_TRY_COMPILE(, double foo (double) __attribute__((__pure__));, cctk_cv_have_cxx_attribute_pure=yes, cctk_cv_have_cxx_attribute_pure=no)
AC_LANG_RESTORE
])
if test "$cctk_cv_have_cxx_attribute_pure" = "yes" ; then
   AC_DEFINE(HAVE_CCTK_CXX_ATTRIBUTE_PURE)
fi
])

AC_DEFUN(CCTK_CXX_MEMBER_ATTRIBUTE_PURE,
[AC_CACHE_CHECK([for C++ member function __attribute__((__pure__))], cctk_cv_have_cxx_member_attribute_pure,
[cctk_cv_have_cxx_member_attribute_pure=no
AC_LANG_SAVE
AC_LANG_CPLUSPLUS
AC_TRY_COMPILE(, struct bar { double foo (double) __attribute__((__pure__)); };, cctk_cv_have_cxx_member_attribute_pure=yes, cctk_cv_have_cxx_member_attribute_pure=no)
AC_LANG_RESTORE
])
if test "$cctk_cv_have_cxx_member_attribute_pure" = "yes" ; then
   AC_DEFINE(HAVE_CCTK_CXX_MEMBER_ATTRIBUTE_PURE)
fi
])



AC_DEFUN(CCTK_C_ATTRIBUTE_NOINLINE,
[AC_CACHE_CHECK([for C function __attribute__((__noinline__))], cctk_cv_have_c_attribute_noinline,
[cctk_cv_have_c_attribute_noinline=no
AC_TRY_COMPILE(, double foo (double) __attribute__((__noinline__));, cctk_cv_have_c_attribute_noinline=yes, cctk_cv_have_c_attribute_noinline=no)
])
if test "$cctk_cv_have_c_attribute_noinline" = "yes" ; then
   AC_DEFINE(HAVE_CCTK_C_ATTRIBUTE_NOINLINE)
fi
])

AC_DEFUN(CCTK_CXX_ATTRIBUTE_NOINLINE,
[AC_CACHE_CHECK([for C++ function __attribute__((__noinline__))], cctk_cv_have_cxx_attribute_noinline,
[cctk_cv_have_cxx_attribute_noinline=no
AC_LANG_SAVE
AC_LANG_CPLUSPLUS
AC_TRY_COMPILE(, double foo (double) __attribute__((__noinline__));, cctk_cv_have_cxx_attribute_noinline=yes, cctk_cv_have_cxx_attribute_noinline=no)
AC_LANG_RESTORE
])
if test "$cctk_cv_have_cxx_attribute_noinline" = "yes" ; then
   AC_DEFINE(HAVE_CCTK_CXX_ATTRIBUTE_NOINLINE)
fi
])

AC_DEFUN(CCTK_CXX_MEMBER_ATTRIBUTE_NOINLINE,
[AC_CACHE_CHECK([for C++ member function __attribute__((__noinline__))], cctk_cv_have_cxx_member_attribute_noinline,
[cctk_cv_have_cxx_member_attribute_noinline=no
AC_LANG_SAVE
AC_LANG_CPLUSPLUS
AC_TRY_COMPILE(, struct bar { double foo (double) __attribute__((__noinline__)); };, cctk_cv_have_cxx_member_attribute_noinline=yes, cctk_cv_have_cxx_member_attribute_noinline=no)
AC_LANG_RESTORE
])
if test "$cctk_cv_have_cxx_member_attribute_noinline" = "yes" ; then
   AC_DEFINE(HAVE_CCTK_CXX_MEMBER_ATTRIBUTE_NOINLINE)
fi
])



AC_DEFUN(CCTK_C_ATTRIBUTE_ALWAYS_INLINE,
[AC_CACHE_CHECK([for C function __attribute__((__always_inline__))], cctk_cv_have_c_attribute_always_inline,
[cctk_cv_have_c_attribute_always_inline=no
AC_TRY_COMPILE(, double foo (double) __attribute__((__always_inline__));, cctk_cv_have_c_attribute_always_inline=yes, cctk_cv_have_c_attribute_always_inline=no)
])
if test "$cctk_cv_have_c_attribute_always_inline" = "yes" ; then
   AC_DEFINE(HAVE_CCTK_C_ATTRIBUTE_ALWAYS_INLINE)
fi
])

AC_DEFUN(CCTK_CXX_ATTRIBUTE_ALWAYS_INLINE,
[AC_CACHE_CHECK([for C++ function __attribute__((__always_inline__))], cctk_cv_have_cxx_attribute_always_inline,
[cctk_cv_have_cxx_attribute_always_inline=no
AC_LANG_SAVE
AC_LANG_CPLUSPLUS
AC_TRY_COMPILE(, double foo (double) __attribute__((__always_inline__));, cctk_cv_have_cxx_attribute_always_inline=yes, cctk_cv_have_cxx_attribute_always_inline=no)
AC_LANG_RESTORE
])
if test "$cctk_cv_have_cxx_attribute_always_inline" = "yes" ; then
   AC_DEFINE(HAVE_CCTK_CXX_ATTRIBUTE_ALWAYS_INLINE)
fi
])

AC_DEFUN(CCTK_CXX_MEMBER_ATTRIBUTE_ALWAYS_INLINE,
[AC_CACHE_CHECK([for C++ member function __attribute__((__always_inline__))], cctk_cv_have_cxx_member_attribute_always_inline,
[cctk_cv_have_cxx_member_attribute_always_inline=no
AC_LANG_SAVE
AC_LANG_CPLUSPLUS
AC_TRY_COMPILE(, struct bar { double foo (double) __attribute__((__always_inline__)); };, cctk_cv_have_cxx_member_attribute_always_inline=yes, cctk_cv_have_cxx_member_attribute_always_inline=no)
AC_LANG_RESTORE
])
if test "$cctk_cv_have_cxx_member_attribute_always_inline" = "yes" ; then
   AC_DEFINE(HAVE_CCTK_CXX_MEMBER_ATTRIBUTE_ALWAYS_INLINE)
fi
])



AC_DEFUN(CCTK_C_ATTRIBUTE_UNUSED,
[AC_CACHE_CHECK([for C __attribute__((__unused__))], cctk_cv_have_c_attribute_unused,
[cctk_cv_have_c_attribute_unused=no
AC_TRY_COMPILE(, double * foo __attribute__((__unused__));, cctk_cv_have_c_attribute_unused=yes, cctk_cv_have_c_attribute_unused=no)
])
if test "$cctk_cv_have_c_attribute_unused" = "yes" ; then
   AC_DEFINE(HAVE_CCTK_C_ATTRIBUTE_UNUSED)
fi
])

AC_DEFUN(CCTK_CXX_ATTRIBUTE_UNUSED,
[AC_CACHE_CHECK([for C++ __attribute__((__unused__))], cctk_cv_have_cxx_attribute_unused,
[cctk_cv_have_cxx_attribute_unused=no
AC_LANG_SAVE
AC_LANG_CPLUSPLUS
AC_TRY_COMPILE(, double * foo __attribute__((__unused__));, cctk_cv_have_cxx_attribute_unused=yes, cctk_cv_have_cxx_attribute_unused=no)
AC_LANG_RESTORE
])
if test "$cctk_cv_have_cxx_attribute_unused" = "yes" ; then
   AC_DEFINE(HAVE_CCTK_CXX_ATTRIBUTE_UNUSED)
fi
])



AC_DEFUN(CCTK_C_ATTRIBUTE_ALIGNED,
[AC_CACHE_CHECK([for C __attribute__((__aligned__(...)))], cctk_cv_have_c_attribute_aligned,
[cctk_cv_have_c_attribute_aligned=no
AC_TRY_COMPILE(, double * foo __attribute__((__aligned__(16)));, cctk_cv_have_c_attribute_aligned=yes, cctk_cv_have_c_attribute_aligned=no)
])
if test "$cctk_cv_have_c_attribute_aligned" = "yes" ; then
   AC_DEFINE(HAVE_CCTK_C_ATTRIBUTE_ALIGNED)
fi
])

AC_DEFUN(CCTK_CXX_ATTRIBUTE_ALIGNED,
[AC_CACHE_CHECK([for C++ __attribute__((__aligned__(...)))], cctk_cv_have_cxx_attribute_aligned,
[cctk_cv_have_cxx_attribute_aligned=no
AC_LANG_SAVE
AC_LANG_CPLUSPLUS
AC_TRY_COMPILE(, double * foo __attribute__((__aligned__(16)));, cctk_cv_have_cxx_attribute_aligned=yes, cctk_cv_have_cxx_attribute_aligned=no)
AC_LANG_RESTORE
])
if test "$cctk_cv_have_cxx_attribute_aligned" = "yes" ; then
   AC_DEFINE(HAVE_CCTK_CXX_ATTRIBUTE_ALIGNED)
fi
])



AC_DEFUN(CCTK_C_ATTRIBUTE_COLD,
[AC_CACHE_CHECK([for C __attribute__((__cold__))], cctk_cv_have_c_attribute_cold,
[cctk_cv_have_c_attribute_cold=no
AC_TRY_COMPILE(, double * foo __attribute__((__cold__));, cctk_cv_have_c_attribute_cold=yes, cctk_cv_have_c_attribute_cold=no)
])
if test "$cctk_cv_have_c_attribute_cold" = "yes" ; then
   AC_DEFINE(HAVE_CCTK_C_ATTRIBUTE_COLD)
fi
])

AC_DEFUN(CCTK_CXX_ATTRIBUTE_COLD,
[AC_CACHE_CHECK([for C++ __attribute__((__cold__))], cctk_cv_have_cxx_attribute_cold,
[cctk_cv_have_cxx_attribute_cold=no
AC_LANG_SAVE
AC_LANG_CPLUSPLUS
AC_TRY_COMPILE(, double * foo __attribute__((__cold__));, cctk_cv_have_cxx_attribute_cold=yes, cctk_cv_have_cxx_attribute_cold=no)
AC_LANG_RESTORE
])
if test "$cctk_cv_have_cxx_attribute_cold" = "yes" ; then
   AC_DEFINE(HAVE_CCTK_CXX_ATTRIBUTE_COLD)
fi
])



AC_DEFUN(CCTK_C_ATTRIBUTE_HOT,
[AC_CACHE_CHECK([for C __attribute__((__hot__))], cctk_cv_have_c_attribute_hot,
[cctk_cv_have_c_attribute_hot=no
AC_TRY_COMPILE(, double * foo __attribute__((__hot__));, cctk_cv_have_c_attribute_hot=yes, cctk_cv_have_c_attribute_hot=no)
])
if test "$cctk_cv_have_c_attribute_hot" = "yes" ; then
   AC_DEFINE(HAVE_CCTK_C_ATTRIBUTE_HOT)
fi
])

AC_DEFUN(CCTK_CXX_ATTRIBUTE_HOT,
[AC_CACHE_CHECK([for C++ __attribute__((__hot__))], cctk_cv_have_cxx_attribute_hot,
[cctk_cv_have_cxx_attribute_hot=no
AC_LANG_SAVE
AC_LANG_CPLUSPLUS
AC_TRY_COMPILE(, double * foo __attribute__((__hot__));, cctk_cv_have_cxx_attribute_hot=yes, cctk_cv_have_cxx_attribute_hot=no)
AC_LANG_RESTORE
])
if test "$cctk_cv_have_cxx_attribute_hot" = "yes" ; then
   AC_DEFINE(HAVE_CCTK_CXX_ATTRIBUTE_HOT)
fi
])



AC_DEFUN(CCTK_C_ATTRIBUTE_FORMAT,
[AC_CACHE_CHECK([for C __attribute__((__format__(printf, 1, 2)))], cctk_cv_have_c_attribute_format,
[cctk_cv_have_c_attribute_format=no
AC_TRY_COMPILE(void xyzzy(const char*, ...) __attribute__((__format__(printf, 1, 2)));, xyzzy("%d",42);, cctk_cv_have_c_attribute_format=yes, cctk_cv_have_c_attribute_format=no)
])
if test "$cctk_cv_have_c_attribute_format" = "yes" ; then
   AC_DEFINE(HAVE_CCTK_C_ATTRIBUTE_FORMAT)
fi
])

AC_DEFUN(CCTK_CXX_ATTRIBUTE_FORMAT,
[AC_CACHE_CHECK([for C++ __attribute__((__format__(printf, 1, 2)))], cctk_cv_have_cxx_attribute_format,
[cctk_cv_have_cxx_attribute_format=no
AC_LANG_SAVE
AC_LANG_CPLUSPLUS
AC_TRY_COMPILE(void xyzzy(const char*, ...) __attribute__((__format__(printf, 1, 2)));, xyzzy("%d",42);, cctk_cv_have_cxx_attribute_format=yes, cctk_cv_have_cxx_attribute_format=no)
AC_LANG_RESTORE
])
if test "$cctk_cv_have_cxx_attribute_format" = "yes" ; then
   AC_DEFINE(HAVE_CCTK_CXX_ATTRIBUTE_FORMAT)
fi
])



AC_DEFUN(CCTK_C_ATTRIBUTE_NORETURN,
[AC_CACHE_CHECK([for C __attribute__((__noreturn__))], cctk_cv_have_c_attribute_noreturn,
[cctk_cv_have_c_attribute_noreturn=no
AC_TRY_COMPILE(void xyzzy(void) __attribute__((__noreturn__));, xyzzy(), cctk_cv_have_c_attribute_noreturn=yes, cctk_cv_have_c_attribute_noreturn=no)
])
if test "$cctk_cv_have_c_attribute_noreturn" = "yes" ; then
   AC_DEFINE(HAVE_CCTK_C_ATTRIBUTE_NORETURN)
fi
])

AC_DEFUN(CCTK_CXX_ATTRIBUTE_NORETURN,
[AC_CACHE_CHECK([for C++ __attribute__((__noreturn__))], cctk_cv_have_cxx_attribute_noreturn,
[cctk_cv_have_cxx_attribute_noreturn=no
AC_LANG_SAVE
AC_LANG_CPLUSPLUS
AC_TRY_COMPILE(void xyzzy(void) __attribute__((__noreturn__));, xyzzy(), cctk_cv_have_cxx_attribute_noreturn=yes, cctk_cv_have_cxx_attribute_noreturn=no)
AC_LANG_RESTORE
])
if test "$cctk_cv_have_cxx_attribute_noreturn" = "yes" ; then
   AC_DEFINE(HAVE_CCTK_CXX_ATTRIBUTE_NORETURN)
fi
])



AC_DEFUN(CCTK_C_ATTRIBUTE_NONNULL,
[AC_CACHE_CHECK([for C __attribute__((__nonnull__))], cctk_cv_have_c_attribute_nonnull,
[cctk_cv_have_c_attribute_nonnull=no
AC_TRY_COMPILE([
    void xyzzy1(void *dest, const void *src, int len)
      __attribute__((__nonnull__ (1,2)));
    void xyzzy1(void *dest, const void *src, int len)
    {
    }
    void xyzzy2(void *dest, const void *src, int len)
      __attribute__((__nonnull__));
    void xyzzy2(void *dest, const void *src, int len)
    {
    }
  ], [
    int a, b;
    xyzzy1(&a, &b, 1);
    xyzzy2(&a, &b, 1);
  ],
  cctk_cv_have_c_attribute_nonnull=yes,
  cctk_cv_have_c_attribute_nonnull=no)
])
if test "$cctk_cv_have_c_attribute_nonnull" = "yes" ; then
   AC_DEFINE(HAVE_CCTK_C_ATTRIBUTE_NONNULL)
fi
])

AC_DEFUN(CCTK_CXX_ATTRIBUTE_NONNULL,
[AC_CACHE_CHECK([for C++ __attribute__((__nonnull__))], cctk_cv_have_cxx_attribute_nonnull,
[cctk_cv_have_cxx_attribute_nonnull=no
AC_LANG_SAVE
AC_LANG_CPLUSPLUS
AC_TRY_COMPILE([
    void xyzzy1(void *dest, const void *src, int len)
      __attribute__((__nonnull__ (1,2)));
    void xyzzy1(void *dest, const void *src, int len)
    {
    }
    void xyzzy2(void *dest, const void *src, int len)
      __attribute__((__nonnull__));
    void xyzzy2(void *dest, const void *src, int len)
    {
    }
  ], [
    int a, b;
    xyzzy1(&a, &b, 1);
    xyzzy2(&a, &b, 1);
  ],
  cctk_cv_have_cxx_attribute_nonnull=yes,
  cctk_cv_have_cxx_attribute_nonnull=no)
AC_LANG_RESTORE
])
if test "$cctk_cv_have_cxx_attribute_nonnull" = "yes" ; then
   AC_DEFINE(HAVE_CCTK_CXX_ATTRIBUTE_NONNULL)
fi
])



AC_DEFUN(CCTK_C_ATTRIBUTE_RETURNS_NONNULL,
[AC_CACHE_CHECK([for C __attribute__((__returns_nonnull__))], cctk_cv_have_c_attribute_returns_nonnull,
[cctk_cv_have_c_attribute_returns_nonnull=no
AC_TRY_COMPILE([
    void* xyzzy(void) __attribute__((__returns_nonnull__));
    void* xyzzy(void)
    {
      static int a;
      return &a;
    }
  ], [
    void* a = xyzzy();
  ],
  cctk_cv_have_c_attribute_returns_nonnull=yes,
  cctk_cv_have_c_attribute_returns_nonnull=no)
])
if test "$cctk_cv_have_c_attribute_returns_nonnull" = "yes" ; then
   AC_DEFINE(HAVE_CCTK_C_ATTRIBUTE_RETURNS_NONNULL)
fi
])

AC_DEFUN(CCTK_CXX_ATTRIBUTE_RETURNS_NONNULL,
[AC_CACHE_CHECK([for C++ __attribute__((__returns_nonnull__))], cctk_cv_have_cxx_attribute_returns_nonnull,
[cctk_cv_have_cxx_attribute_returns_nonnull=no
AC_LANG_SAVE
AC_LANG_CPLUSPLUS
AC_TRY_COMPILE([
    void* xyzzy(void) __attribute__((__returns_nonnull__));
    void* xyzzy(void)
    {
      static int a;
      return &a;
    }
  ], [
    void* a = xyzzy();
  ],
  cctk_cv_have_cxx_attribute_returns_nonnull=yes,
  cctk_cv_have_cxx_attribute_returns_nonnull=no)
AC_LANG_RESTORE
])
if test "$cctk_cv_have_cxx_attribute_returns_nonnull" = "yes" ; then
   AC_DEFINE(HAVE_CCTK_CXX_ATTRIBUTE_RETURNS_NONNULL)
fi
])



AC_DEFUN(CCTK_C_BUILTIN_EXPECT,
[AC_CACHE_CHECK([for C __builtin_expect], cctk_cv_have_c_builtin_expect,
[cctk_cv_have_c_builtin_expect=no
AC_TRY_LINK(, __builtin_expect(0,0);, cctk_cv_have_c_builtin_expect=yes, cctk_cv_have_c_builtin_expect=no)
])
if test "$cctk_cv_have_c_builtin_expect" = "yes" ; then
   AC_DEFINE(HAVE_CCTK_C_BUILTIN_EXPECT)
fi
])

AC_DEFUN(CCTK_CXX_BUILTIN_EXPECT,
[AC_CACHE_CHECK([for C++ __builtin_expect], cctk_cv_have_cxx_builtin_expect,
[cctk_cv_have_cxx_builtin_expect=no
AC_LANG_SAVE
AC_LANG_CPLUSPLUS
AC_TRY_LINK(, __builtin_expect(0,0);, cctk_cv_have_cxx_builtin_expect=yes, cctk_cv_have_cxx_builtin_expect=no)
AC_LANG_RESTORE
])
if test "$cctk_cv_have_cxx_builtin_expect" = "yes" ; then
   AC_DEFINE(HAVE_CCTK_CXX_BUILTIN_EXPECT)
fi
])



AC_DEFUN(CCTK_C_BUILTIN_UNREACHABLE,
[AC_CACHE_CHECK([for C __builtin_unreachable], cctk_cv_have_c_builtin_unreachable,
[cctk_cv_have_c_builtin_unreachable=no
AC_TRY_LINK(, __builtin_unreachable();, cctk_cv_have_c_builtin_unreachable=yes, cctk_cv_have_c_builtin_unreachable=no)
])
if test "$cctk_cv_have_c_builtin_unreachable" = "yes" ; then
   AC_DEFINE(HAVE_CCTK_C_BUILTIN_UNREACHABLE)
fi
])

AC_DEFUN(CCTK_CXX_BUILTIN_UNREACHABLE,
[AC_CACHE_CHECK([for C++ __builtin_unreachable], cctk_cv_have_cxx_builtin_unreachable,
[cctk_cv_have_cxx_builtin_unreachable=no
AC_LANG_SAVE
AC_LANG_CPLUSPLUS
AC_TRY_LINK(, __builtin_unreachable();, cctk_cv_have_cxx_builtin_unreachable=yes, cctk_cv_have_cxx_builtin_unreachable=no)
AC_LANG_RESTORE
])
if test "$cctk_cv_have_cxx_builtin_unreachable" = "yes" ; then
   AC_DEFINE(HAVE_CCTK_CXX_BUILTIN_UNREACHABLE)
fi
])



AC_DEFUN(CCTK_C_BUILTIN_ASSUME_ALIGNED,
[AC_CACHE_CHECK([for C __builtin_assume_aligned], cctk_cv_have_c_builtin_assume_aligned,
[cctk_cv_have_c_builtin_assume_aligned=no
AC_TRY_LINK(,
__builtin_assume_aligned((void*)1000, 10);
__builtin_assume_aligned((void*)1001, 10, 1);
, cctk_cv_have_c_builtin_assume_aligned=yes, cctk_cv_have_c_builtin_assume_aligned=no)
])
if test "$cctk_cv_have_c_builtin_assume_aligned" = "yes" ; then
   AC_DEFINE(HAVE_CCTK_C_BUILTIN_ASSUME_ALIGNED)
fi
])

AC_DEFUN(CCTK_CXX_BUILTIN_ASSUME_ALIGNED,
[AC_CACHE_CHECK([for C++ __builtin_assume_aligned], cctk_cv_have_cxx_builtin_assume_aligned,
[cctk_cv_have_cxx_builtin_assume_aligned=no
AC_LANG_SAVE
AC_LANG_CPLUSPLUS
AC_TRY_LINK(,
__builtin_assume_aligned((void*)1000, 10);
__builtin_assume_aligned((void*)1001, 10, 1);
, cctk_cv_have_cxx_builtin_assume_aligned=yes, cctk_cv_have_cxx_builtin_assume_aligned=no)
AC_LANG_RESTORE
])
if test "$cctk_cv_have_cxx_builtin_assume_aligned" = "yes" ; then
   AC_DEFINE(HAVE_CCTK_CXX_BUILTIN_ASSUME_ALIGNED)
fi
])



AC_DEFUN(CCTK_CXX_STATIC_ASSERT,
[AC_CACHE_CHECK([for C++ static_assert], cctk_cv_have_cxx_static_assert,
[cctk_cv_have_cxx_static_assert=no
AC_LANG_SAVE
AC_LANG_CPLUSPLUS
AC_TRY_LINK(, static_assert(1, "good");, cctk_cv_have_cxx_static_assert=yes, cctk_cv_have_cxx_static_assert=no)
AC_LANG_RESTORE
])
if test "$cctk_cv_have_cxx_static_assert" = "yes" ; then
   AC_DEFINE(HAVE_CCTK_CXX_STATIC_ASSERT)
fi
])



dnl Check for a function that may be provided by cmath or math.h
AC_DEFUN(CCTK_CHECK_CXX_STDMATHFUNC,
[cctk_func=`echo $1 | sed 'y%./+-%__p_%'`
 cctk_tr_func=`echo $cctk_func |
   sed -e 'y/abcdefghijklmnopqrstuvwxyz/ABCDEFGHIJKLMNOPQRSTUVWXYZ/'`
AC_MSG_CHECKING([for C++ $1])
AC_CACHE_VAL(cctk_cv_cxx_$cctk_func,
[AC_LANG_SAVE
AC_LANG_CPLUSPLUS
cctk_cv_cxx_func=no
for ac_func in "std::$cctk_func" "$cctk_func" "::$cctk_func"; do
for ac_nargs in 1 2; do
  case $ac_nargs in
    1) ac_args='(1.0)'; ac_argsf='(1.0f)' ;;
    2) ac_args='(1.0, 1.0)'; ac_argsf='(1.0f, 1.0f)' ;;
  esac
  AC_TRY_COMPILE([
/* See note in cctk_Math.h regarding these include statements */
#include <cmath>
#include <math.h>
], [
{
  $ac_func $ac_argsf;
  $ac_func $ac_args;
}
using namespace std;
{
  $ac_func $ac_argsf;
  $ac_func $ac_args;
}
], [cctk_cv_cxx_func="$ac_func"; break 2])
done
done
AC_LANG_RESTORE
eval cctk_cv_cxx_$cctk_func=\$cctk_cv_cxx_func
])
AC_MSG_RESULT($cctk_cv_cxx_func)
case "$cctk_cv_cxx_func" in
  no) : ;;
  *)  AC_DEFINE_UNQUOTED(CCTK_CXX_$cctk_tr_func, $cctk_cv_cxx_func)
      AC_DEFINE_UNQUOTED(HAVE_CCTK_CXX_$cctk_tr_func, 1) ;;
esac
])



AC_DEFUN(CCTK_CXX_AUTO_SPECIFIER,
[AC_CACHE_CHECK([for C++ auto specifier], cctk_cv_have_cxx_auto_specifier,
[cctk_cv_have_cxx_auto_specifier=no
AC_LANG_SAVE
AC_LANG_CPLUSPLUS
AC_TRY_LINK(,
  int x; auto y = x,
  cctk_cv_have_cxx_auto_specifier=yes,
  cctk_cv_have_cxx_auto_specifier=no)
AC_LANG_RESTORE
])
if test "$cctk_cv_have_cxx_auto_specifier" = "yes" ; then
   AC_DEFINE(HAVE_CCTK_CXX_AUTO_SPECIFIER)
fi
])



AC_DEFUN(CCTK_CXX_LAMBDA,
[AC_CACHE_CHECK([for C++ lambda expressions], cctk_cv_have_cxx_lambda,
[cctk_cv_have_cxx_lambda=no
AC_LANG_SAVE
AC_LANG_CPLUSPLUS
AC_TRY_LINK(,
  [int x; [x](int y) { return x+y; };],
  cctk_cv_have_cxx_lambda=yes,
  cctk_cv_have_cxx_lambda=no)
AC_LANG_RESTORE
])
if test "$cctk_cv_have_cxx_lambda" = "yes" ; then
   AC_DEFINE(HAVE_CCTK_CXX_LAMBDA)
fi
])



AC_DEFUN(CCTK_CXX_RANGE_BASED_FOR,
[AC_CACHE_CHECK([for C++ range-based for statements], cctk_cv_have_cxx_range_based_for,
[cctk_cv_have_cxx_range_based_for=no
AC_LANG_SAVE
AC_LANG_CPLUSPLUS
AC_TRY_LINK(
  [#include <vector>],
  [std::vector<int> xs(10);
   for (int& x: xs) x = 42;],
  cctk_cv_have_cxx_range_based_for=yes,
  cctk_cv_have_cxx_range_based_for=no)
AC_LANG_RESTORE
])
if test "$cctk_cv_have_cxx_range_based_for" = "yes" ; then
   AC_DEFINE(HAVE_CCTK_CXX_RANGE_BASED_FOR)
fi
])



AC_DEFUN(CCTK_CHECK_C99,
[AC_CACHE_CHECK([for C99 features], cctk_cv_have_c99,
[cctk_cv_have_c99=no
AC_TRY_COMPILE(,
// Test C99 comments
// Test variable declarations in the middle of blocks
0;
int x = 1;
for (int i=0; i<10; ++i) x+=i;
, cctk_cv_have_c99=yes, cctk_cv_have_c99=no)
])
if test "$cctk_cv_have_c99" = "yes" ; then
   AC_DEFINE(HAVE_CCTK_C99)
fi
])



dnl CCTK_CHECK_DEFINED(DEFINED, [[ADDITIONAL_CODE [, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
AC_DEFUN(CCTK_CHECK_DEFINED,
[dnl Do the transliteration at runtime so arg 1 can be a shell variable.
cctk_safe=`echo "$1" | sed 'y%./+-%__p_%'`
AC_MSG_CHECKING([for $1])
AC_CACHE_VAL(cctk_cv_defined_$cctk_safe,
[AC_TRY_COMPILE([$2
#ifndef $1
#error "$1 not defined"
#endif], [ ], eval "cctk_cv_defined_$cctk_safe=yes",
  eval "cctk_cv_defined_$cctk_safe=no")])dnl
if eval "test \"`echo '$cctk_cv_defined_'$cctk_safe`\" = yes"; then
  AC_MSG_RESULT(yes)
  ifelse([$3], , :, [$3])
else
  AC_MSG_RESULT(no)
ifelse([$4], , , [$4])
fi
])



dnl CCTK_CHECK_HAS_PROTOTYPE(FUNCTION, [FEATURE-DESCRIPTION [, ADDITIONAL_CODE [, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]]])
AC_DEFUN(CCTK_CHECK_HAS_PROTOTYPE,
[dnl Do the transliteration at runtime so arg 1 can be a shell variable.
cctk_safe=`echo "$1" | sed 'y%./+-%__p_%'`
cctk_lang=AC_LANG
AC_MSG_CHECKING([ifelse([$2], , [for $1], [$2])])
AC_CACHE_VAL(cctk_cv_has_${cctk_lang}_prototype_$cctk_safe,
[AC_TRY_COMPILE([$3
], [char *p = (char*) $1;], eval "cctk_cv_has_${cctk_lang}_prototype_$cctk_safe=yes",
  eval "cctk_cv_has_${cctk_lang}_prototype_$cctk_safe=no")])dnl
if eval "test \"`echo '$cctk_cv_has_'$cctk_lang'_prototype_'$cctk_safe`\" = yes"; then
  AC_MSG_RESULT(yes)
  ifelse([$4], , :, [$4])
else
  AC_MSG_RESULT(no)
ifelse([$5], , , [$5])
fi
])



dnl CCTK_C_CHECK_HAS_PROTOTYPE(FUNCTION, [FEATURE-DESCRIPTION [, ADDITIONAL_CODE [, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]]])
AC_DEFUN(CCTK_C_CHECK_HAS_PROTOTYPE,
[dnl short-hand for C protype check
CCTK_CHECK_HAS_PROTOTYPE([$1], [ifelse([$2], , [for C $1], [$2])],
[$3], [$4], [$5])dnl
])



dnl CCTK_CXX_CHECK_HAS_PROTOTYPE(FUNCTION, [FEATURE-DESCRIPTION [, ADDITIONAL_CODE [, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]]])
AC_DEFUN(CCTK_CXX_CHECK_HAS_PROTOTYPE,
[dnl short-hand for C++ protype check
AC_LANG_SAVE
AC_LANG_CPLUSPLUS
CCTK_CHECK_HAS_PROTOTYPE([$1], [ifelse([$2], , [for C++ $1], [$2])],
[$3], [$4], [$5])dnl
AC_LANG_RESTORE
])



AC_DEFUN(CCTK_FORTRAN_CRAY_POINTERS,
[AC_CACHE_CHECK([for Fortran Cray pointers], cctk_cv_have_fortran_cray_pointers,
[cctk_cv_have_fortran_cray_pointers=no
AC_LANG_SAVE
AC_LANG_FORTRAN77
CCTK_TRY_COMPILE(
[
      subroutine sub(pointers, n)
      implicit none

C     Find the integer type used for pointers
      integer dummy
      pointer (pdummy, dummy)
      integer, parameter :: pk = kind(pdummy)

C     An array of pointers that is passed in
      integer(pk) pointers(3)

C     The array size
      integer n

C     Declare local variables which are pointers, using the Cray pointer
C     extension

C     Explanation: The variables "a" and "pa" belong together. Only "pa"
C     is a variable. Whenever "a" is used, the pointer which is stored
C     in "pa" is dereferenced. In C, one would write "*pa" instead of
C     "a".
      double precision a(n,n), b(n,n), c(n,n)
      pointer (pa, a)
      pointer (pb, b)
      pointer (pc, c)

C     Local loop indices
      integer i, j, k

C     Set the pointers from the array which is passed in
      pa = pointers(1)
      pb = pointers(2)
      pc = pointers(3)

C     Do some work on the arrays, as if they were not pointers
      do i = 1, n
         do j = 1, n
            a(i,j) = 0
            do k = 1, n
               a(i,j) = a(i,j) + b(i,k) * c(k,j)
            end do
         end do
      end do

      end subroutine sub
],
[
      implicit none

      integer, parameter :: n = 10
      double precision a(n,n), b(n,n), c(n,n)

C     Find the integer type used for pointers
      integer dummy
      pointer (pdummy, dummy)
      integer, parameter :: pk = kind(pdummy)

      integer(pk) pointers(3)
      pointers(1) = %loc(a)
      pointers(2) = %loc(b)
      pointers(3) = %loc(b)
      call sub(pointers, n)
],
cctk_cv_have_fortran_cray_pointers=yes,
cctk_cv_have_fortran_cray_pointers=no)
AC_LANG_RESTORE
])
if test "$cctk_cv_have_fortran_cray_pointers" = "yes" ; then
   AC_DEFINE(HAVE_CCTK_FORTRAN_CRAY_POINTERS)
fi
])
