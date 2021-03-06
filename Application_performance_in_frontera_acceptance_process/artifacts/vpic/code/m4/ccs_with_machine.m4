dnl ----------------------------------------------------------------------------
dnl CCS_WITH_MACHINE
dnl
dnl ARG1: machine file directory ($1)
dnl ----------------------------------------------------------------------------
AC_DEFUN([CCS_WITH_MACHINE], [
    AC_ARG_VAR([MACHINE_FILE], [Floating Point Type e.g. double, float, ...])

    AC_ARG_WITH([machine],
        [AC_HELP_STRING([--with-machine],
        [Specify machine file for configuration])],
        [
            if test -n "$MACHINE_FILE" ; then
                AC_SUBST(MACHINE_FILE, $MACHINE_FILE)
                AC_SUBST(HAS_MACHINE_FILE, [1])
            elif test "$withval" != no ; then
                AC_SUBST(MACHINE_FILE, $withval)
                AC_SUBST(HAS_MACHINE_FILE, [1])
            fi
        ],
        [
            if test -n "$MACHINE_FILE" ; then
                AC_SUBST(MACHINE_FILE, $MACHINE_FILE)
                AC_SUBST(HAS_MACHINE_FILE, [1])
            else
                AC_SUBST(HAS_MACHINE_FILE, [0])
            fi
    ])

    AC_MSG_CHECKING(machine specification)
    if test "$HAS_MACHINE_FILE" = "1" ; then
        if test -f "$1/$MACHINE_FILE" ; then
            . "$1/$MACHINE_FILE"
            AC_MSG_RESULT($MACHINE_FILE)
        elif test -f "$MACHINE_FILE" ; then
            . "$MACHINE_FILE"
            AC_MSG_RESULT($MACHINE_FILE)
        else
            AC_MSG_RESULT(using config options)
        fi
    else
        AC_MSG_RESULT(using config options)
    fi
])
