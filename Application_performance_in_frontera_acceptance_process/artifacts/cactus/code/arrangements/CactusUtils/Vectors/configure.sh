#! /bin/bash

# Set up shell
if [ "$(echo ${VERBOSE} | tr '[:upper:]' '[:lower:]')" = 'yes' ]; then
    set -x                      # Output commands
fi
set -e                          # Abort on errors



################################################################################
# Examine settings
################################################################################



case $(echo "x$VECTORISE" | tr '[:upper:]' '[:lower:]') in
    (xyes) VECTORISE=1 ;;
    (xno)  VECTORISE=0 ;;
    (x)    VECTORISE=0 ;;        # default
    (*)    echo "BEGIN ERROR"
           echo "Illegal value of option VECTORISE"
           echo "END ERROR"
           exit 1;;
esac

# Disable vectorisation when optimisation is disabled. (The
# vectorisation macros depend on optimisation for efficient code;
# without optimisation, the code is most likely much slower than
# usual.)
case $(echo "x$OPTIMISE$OPTIMIZE" | tr '[:upper:]' '[:lower:]') in
    (x)    ;;                   # treat as 'yes'
    (xyes) ;;                   # do nothing
    (xno)  VECTORISE=0 ;;       # disable vectorisation
    (*)    echo "BEGIN ERROR"
           echo "Illegal value of option OPTIMISE or OPTIMIZE"
           echo "END ERROR"
           exit 1
esac

case $(echo "x$VECTORISE_ALIGNED_ARRAYS" | tr '[:upper:]' '[:lower:]') in
    (xyes) VECTORISE_ALIGNED_ARRAYS=1 ;;
    (xno)  VECTORISE_ALIGNED_ARRAYS=0 ;;
    (x)    VECTORISE_ALIGNED_ARRAYS=0 ;; # default
    (*)    echo "BEGIN ERROR"
           echo "Illegal value of option VECTORISE_ALIGNED_ARRAYS"
           echo "END ERROR"
           exit 1
esac

case $(echo "x$VECTORISE_ALWAYS_USE_UNALIGNED_LOADS" | tr '[:upper:]' '[:lower:]') in
    (xyes) VECTORISE_ALWAYS_USE_UNALIGNED_LOADS=1 ;;
    (xno)  VECTORISE_ALWAYS_USE_UNALIGNED_LOADS=0 ;;
    (x)    VECTORISE_ALWAYS_USE_UNALIGNED_LOADS=0 ;; # default
    (*)    echo "BEGIN ERROR"
           echo "Illegal value of option VECTORISE_ALWAYS_USE_UNALIGNED_LOADS"
           echo "END ERROR"
           exit 1
esac

case $(echo "x$VECTORISE_ALWAYS_USE_ALIGNED_LOADS" | tr '[:upper:]' '[:lower:]') in
    (xyes) VECTORISE_ALWAYS_USE_ALIGNED_LOADS=1 ;;
    (xno)  VECTORISE_ALWAYS_USE_ALIGNED_LOADS=0 ;;
    (x)    VECTORISE_ALWAYS_USE_ALIGNED_LOADS=0 ;; # default
    (*)    echo "BEGIN ERROR"
           echo "Illegal value of option VECTORISE_ALWAYS_USE_ALIGNED_LOADS"
           echo "END ERROR"
           exit 1
esac

case $(echo "x$VECTORISE_STREAMING_STORES" | tr '[:upper:]' '[:lower:]') in
    (xyes) VECTORISE_STREAMING_STORES=1 ;;
    (xno)  VECTORISE_STREAMING_STORES=0 ;;
    (x)    VECTORISE_STREAMING_STORES=0 ;; # default
    (*)    echo "BEGIN ERROR"
           echo "Illegal value of option VECTORISE_STREAMING_STORES"
           echo "END ERROR"
           exit 1
esac

case $(echo "x$VECTORISE_INLINE" | tr '[:upper:]' '[:lower:]') in
    (xyes) VECTORISE_INLINE=1 ;;
    (xno)  VECTORISE_INLINE=0 ;;
    (x)    VECTORISE_INLINE=0 ;; # default
    (*)    echo "BEGIN ERROR"
           echo "Illegal value of option VECTORISE_INLINE"
           echo "END ERROR"
           exit 1
esac



################################################################################
# Configure Cactus
################################################################################

# Pass options to Cactus
echo "BEGIN DEFINE"
echo "VECTORISE                            $VECTORISE"
echo "VECTORISE_ALIGNED_ARRAYS             $VECTORISE_ALIGNED_ARRAYS"
echo "VECTORISE_ALWAYS_USE_UNALIGNED_LOADS $VECTORISE_ALWAYS_USE_UNALIGNED_LOADS"
echo "VECTORISE_ALWAYS_USE_ALIGNED_LOADS   $VECTORISE_ALWAYS_USE_ALIGNED_LOADS"
echo "VECTORISE_INLINE                     $VECTORISE_INLINE"
echo "VECTORISE_STREAMING_STORES           $VECTORISE_STREAMING_STORES"
echo "END DEFINE"

echo "BEGIN MAKE_DEFINITION"
echo "VECTORISE                            = $VECTORISE"
echo "VECTORISE_ALIGNED_ARRAYS             = $VECTORISE_ALIGNED_ARRAYS"
echo "VECTORISE_ALWAYS_USE_UNALIGNED_LOADS = $VECTORISE_ALWAYS_USE_UNALIGNED_LOADS"
echo "VECTORISE_ALWAYS_USE_ALIGNED_LOADS   = $VECTORISE_ALWAYS_USE_ALIGNED_LOADS"
echo "VECTORISE_INLINE                     = $VECTORISE_INLINE"
echo "VECTORISE_STREAMING_STORES           = $VECTORISE_STREAMING_STORES"
echo "END MAKE_DEFINITION"
