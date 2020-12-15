#! /bin/sh

# Set up shell
if [ "$(echo ${VERBOSE} | tr '[:upper:]' '[:lower:]')" = 'yes' ]; then
    set -x                      # Output commands
fi


if [ "${F77}" = "none" ] && [ "${F90}" = "none" ]; then
    echo 'BEGIN ERROR'
    echo "Fortran thorn requires that a Fortran compiler is defined, but F77 = '$F77' and F90 = '$F90'.  Aborting."
    echo 'END ERROR'
    exit 1
fi

# Try to use ## for token concatenation.
# If this works, we likely have an ANSI cpp.

rm config.tmp config.out 2> /dev/null

if test -z "${FPP}"; then
    # FPP is not defined; try to guess which fpp will be chosen later
    FPP=cpp
fi

cat > config.tmp <<EOF
#define CONCAT(a,b) a##b
CONCAT(hello,world)
EOF
${FPP} ${FPPFLAGS} config.tmp > config.out
grep 'helloworld' config.out > /dev/null 2> /dev/null
# grep returns 0 for success, non-zero for failure
cpp_ansi=$?
rm config.tmp config.out 2> /dev/null

cat > config.tmp <<EOF
#define CONCAT(a,b) a/**/b
CONCAT(hello,world)
EOF
${FPP} ${FPPFLAGS} config.tmp > config.out
grep 'helloworld' config.out > /dev/null 2> /dev/null
# grep returns 0 for success, non-zero for failure
cpp_traditional=$?
rm config.tmp config.out 2> /dev/null



if test ${cpp_ansi} = 0; then
    echo 'BEGIN MESSAGE'
    echo 'Found an ANSI-like Fortran cpp'
    echo 'END MESSAGE'
    echo 'BEGIN DEFINE'
    echo 'FORTRAN_CPP_ANSI 1'
    echo 'END DEFINE'
elif test ${cpp_traditional} = 0; then
    echo 'BEGIN MESSAGE'
    echo 'Found a traditional Fortran cpp'
    echo 'END MESSAGE'
    echo 'BEGIN DEFINE'
    echo 'FORTRAN_CPP_ANSI 0'
    echo 'END DEFINE'
else
    echo 'BEGIN ERROR'
    echo 'No Fortran preprocessor defined'
    echo 'END ERROR'
fi
