#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

# For help, see :
# - info:/automake/General%20Index
# - info:/autoconf/Autoconf%20Macro%20Index

AC_PREREQ(2.59)
#See info:/autoconf/Initializing%20configure
AC_INIT(tools, 2.8, florent.lyard@legos.obs-mip.fr)
#See info:/autoconf/Input
AC_CONFIG_SRCDIR([src/archive.h])
AC_CONFIG_HEADER([config.h])
# AC_CONFIG_MACRO_DIR([m4])


echo "############################################"
echo "#"
echo "#            Standard checks"
echo "#"
echo "############################################"

# see info:/automake/Public%20Macros
AM_INIT_AUTOMAKE(foreign) # info:/automake/Options

# for MAC users
AC_PROG_SED
AC_CHECK_PROGS(READLINK,[greadlink readlink])

echo "################"
echo "# compile mode #"
echo "################"
# see info:/autoconf/External%20Software for help on AC_ARG_WITH
AC_ARG_WITH(debug,[AS_HELP_STRING([--with-debug],[make debugging easier by lowering optimisation. You can also set CFLAGS. Default: "-O3", i.e. aggressive optimisation without debugging. CXXFLAGS will be overwritten to the value of CFLAGS.])])
if test "$with_debug" = yes;then
  # even -O1 is annoying ...
  CFLAGS="-g"
elif test "x$CFLAGS" = x;then
  CFLAGS="-O3"
  fi

#CXXFLAGS="$CFLAGS" # done in objects/Makefile.am for many good reasons

# info:/autoconf/C++%20Compiler info:/autoconf/C%20Compiler and info:/autoconf/Fortran%20Compiler
# say that these 3 macros set CXXFLAGS=CFLAGS=FCFLAGS="-g -O2" unless already set
AC_PROG_CXX
AC_PROG_CC
AC_PROG_FC
echo "CXX=\"$CXX\" CFLAGS=\"$CFLAGS\""
echo "CPPFLAGS=\"$CPPFLAGS\""

AM_PROG_LEX
AC_PROG_YACC

AC_MSG_CHECKING([whether automake produces .h or .hpp from .ypp files])
AC_DEFINE(AM_YPP_TARGET_IS_HPP,[],[whether automake actually produces .hpp from .ypp files])
# info:/automake/Yacc and Lex
# the test is crude, but should work for most distributions, including Debian
if test $am__api_version = 1.10 -o $am__api_version = 1.11 ;then
  AM_YPP_TARGET=h
  AC_DEFINE([AM_YPP_TARGET_IS_HPP],[0])
else
  AM_YPP_TARGET=hpp
  AC_DEFINE([AM_YPP_TARGET_IS_HPP],[1])
  fi
AC_MSG_RESULT([.$AM_YPP_TARGET])
AC_SUBST(AM_YPP_TARGET,[$AM_YPP_TARGET])

# EITHER
dnl AC_CONFIG_MACRO_DIR([m4])
dnl LT_INIT
# OR
AC_PROG_RANLIB
# EO EITHER/OR

AC_PROG_MAKE_SET


echo "######################"
echo "#   symbolic links   #"
echo "######################"
# SEPARATE BUILD DIRECTORY TO ALLOW SEVERAL PEOPLE TO SHARE THE SAME SOURCE
rm objects/{cleanCmdline.py,ttb-test.sh,ncRmDim.py,plotCD.py,topo_merge}
mkdir -p objects
ln -vs `$READLINK -f $srcdir`/src/{cleanCmdline.py,ttb-test.sh,ncRmDim.py,plotCD.py} objects/
ln -vs `$READLINK -f $srcdir`/src/poc-formula-{scan.lpp,parse.ypp} objects/
ln -vs topo-merge objects/topo_merge

echo "##################################################"
echo "# Checks for miscellaneous libraries and headers"
echo "##################################################"

AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h string.h strings.h unistd.h shapefil.h libshp/shapefil.h])
AC_TYPE_UINT64_T

AC_CHECK_LIB([shp],[main], ,)
AC_CHECK_LIB([m],[main], ,)
AC_CHECK_LIB([rt],[main], ,)

AC_CHECK_LIB([proj],[pj_init_plus_ctx], ,[AC_MSG_ERROR([Can not find proj: install proj 4.8 to 5 (preferably proj 4.9.0b2 to 5.1.0) or give LDFLAGS="\$LDFLAGS -L/proj_PATH/lib" option to configure])])
AC_CHECK_DECLS([pj_init_plus_ctx],[],[AC_MSG_ERROR([Can not find pj_init_plus_ctx in proj_api.h: install proj 4.8 to 5 (preferably proj 4.9.0b2 to 5.1.0) or give CPPFLAGS="\$CPPFLAGS -I/proj_PATH/include" option to configure])],[#include <proj_api.h>])

AC_CHECK_LIB([gslcblas],[main], ,[AC_MSG_ERROR([Can not find gslcblas: install it or give LDFLAGS="\$LDFLAGS -L/gsl_PATH/lib" option to configure])])
AC_CHECK_LIB([gsl],[main], ,[AC_MSG_ERROR([Can not find gsl: install it or give LDFLAGS="\$LDFLAGS -L/gsl_PATH/lib" option to configure])])
AC_CHECK_HEADERS([gsl/gsl_cblas.h],[],[AC_MSG_ERROR([Can not find gsl/gsl_cblas.h: install it or give CPPFLAGS="\$CPPFLAGS -I/gsl_PATH/include" option to configure])])


echo "############################################"
echo "# Checking netcdf library and dependencies"
echo "############################################"
AC_CHECK_HEADERS([netcdf.h],[],[AC_MSG_ERROR([Can not find netcdf.h: install it or give CPPFLAGS="\$CPPFLAGS -I/netcdf_PATH/include" option to configure])])
AC_CHECK_LIB([netcdf],[main], ,[AC_MSG_ERROR([Can not find netcdf: install it or give LDFLAGS="\$LDFLAGS -L/netcdf_PATH/lib" option to configure])])
AC_CHECK_LIB([netcdff],[main])

###########################
# Checking whether netcdf library is thread safe (IT IS NOT SO FAR)
# netcdf-4.1.1/libsrc/nc.c:23 shows NetCDF is not thread safe
# netcdf-4.1.3/libsrc/nclistmgr.c also shows NetCDF is not thread safe
# info:/autoconf/Shell Substitutions
if false;then
NETCDF_VERSION=`nc-config --version | $SED -e 's/.* //'`
# info:/autoconf/Defining Symbols
AC_DEFINE_UNQUOTED([NETCDF_VERSION],"$NETCDF_VERSION",[netcdf version number])
fi
# info:/autoconf/Common Shell Constructs
AC_DEFINE([NETCDF_CAN_PARALLEL_IO], [0], [whether netcdf library is thread safe])
# the day it is, if ever, you can replace the line above by the 2 lines below, putting the right version number
# AC_DEFINE([NETCDF_CAN_PARALLEL_IO], [], [whether netcdf library is thread safe])
# AS_VERSION_COMPARE($NETCDF_VERSION,[9.9.9],AC_DEFINE([NETCDF_CAN_PARALLEL_IO],[0]),AC_DEFINE([NETCDF_CAN_PARALLEL_IO],[1]),AC_DEFINE([NETCDF_CAN_PARALLEL_IO],[1]))
###########################

echo "#####################"
echo "# Checking GRIB_API #"
echo "#####################"
AC_CHECK_LIB([grib_api],[grib_get_data], ,[AC_MSG_WARN([Can not find grib_get_data() in grib_api: install it or give LDFLAGS="\$LDFLAGS -L/GRIB_API_PATH/lib" option to configure to get GRIB support])])
if test $ac_cv_lib_grib_api_grib_get_data = yes ; then
  AC_CHECK_HEADERS([grib_api.h], ,[AC_MSG_ERROR([Can not find grib_api.h: check readline installation or give CPPFLAGS="\$CPPFLAGS -I/GRIB_API_PATH/include" option to configure])])
  fi


###########################
# Checks for poc-solver
###########################

echo "############################################"
echo "#  checking for dependencies directory"
echo "############################################"
#see info:/autoconf/Autoconf%20Language
#Note : info:/autoconf/Pretty%20Help%20Strings does not work !
AC_ARG_VAR(SOFTSDIR,[The dependencies directory where to automatically find the poc-solver directory.
IT DOES NOT CHECK FOR /home/softs BECAUSE IT IS WRONG TO DO SO ON SERVERS!
If not set, it is set to $srcdir/../])
if test "$SOFTSDIR";then
  echo -n "SOFTSDIR already"
else
  d=`$READLINK -f $srcdir/..`
  echo -n "SOFTSDIR"
  fi
AC_SUBST(SOFTSDIR,[${SOFTSDIR:="$d"}])
echo " set to $SOFTSDIR"

echo "############################################"
echo "#"
echo "# checking for poc-solver includes and library"
echo "#"
echo "############################################"

AC_ARG_VAR(POCSOLVERDIR,[The poc-solver directory. If not set, it automatically checks for, in that order :
$SOFTSDIR/poc-solvers $srcdir/../poc-solvers $srcdir/../../poc-solvers
])
test "$POCSOLVERBUILDDIR" ||
if test "$POCSOLVERDIR";then
  echo -n "POCSOLVERDIR already"
else
  for d in $SOFTSDIR/poc-solvers* $srcdir/../poc-solvers* $srcdir/../../poc-solvers*;do
    AC_MSG_CHECKING($d )
    if test -d "$d";then
      if echo "$d" | grep -qE "~$" ;then
        AC_MSG_RESULT(skipped)
        continue
        fi
      if test `$d/cmade/bin/testlib --version | sed -re "s/\..*//"` != 1 ;then
        AC_MSG_RESULT(wrong version)
        continue
        fi
      AC_MSG_RESULT(found)
      d=`$READLINK -f $d`
      break
    else
      AC_MSG_RESULT(does not exist)
      d=
      fi
    done
  if test -z "$d";then #not found
    echo "###########################################"
    echo "# FAILED TO FIND POC-SOLVER DIRECTORY !!! #"
    echo "###########################################"
    AC_MSG_ERROR([On this machine, you MUST give POCSOLVERDIR=... or SOFTSDIR=... option to configure])
    fi
  echo -n "POCSOLVERDIR"
  fi
AC_SUBST(POCSOLVERDIR,[${POCSOLVERDIR:="$d"}])
echo " set to $POCSOLVERDIR"

AC_ARG_VAR(POCSOLVERBUILDDIR,[The poc-solver build directory. If not set, it automatically checks for, in that order :
$POCSOLVERDIR/cmade $POCSOLVERDIR/sequential $POCSOLVERDIR/parallel
])
if test "$POCSOLVERBUILDDIR";then
  echo -n "POCSOLVERBUILDDIR already"
else
  for d in $POCSOLVERDIR/cmade $POCSOLVERDIR/sequential $POCSOLVERDIR/parallel;do
    AC_MSG_CHECKING($d )
    # setting d anyway to make sure POCSOLVERBUILDDIR is set to a non-existant directory
    # rather than "" to avoid errors with /include and /lib directories later
    d=`$READLINK -f $d`
    if test -d "$d";then
      AC_MSG_RESULT(found)
      break
    else
      AC_MSG_RESULT(does not exist)
      fi
    done
  if test \! -d "$d";then #not found
    echo "#################################################"
    echo "# FAILED TO FIND POC-SOLVER BUILD DIRECTORY !!! #"
    echo "#################################################"
    AC_MSG_WARN([On this machine, you may need to give POCSOLVERBUILDDIR=... option to configure])
    fi
  echo -n "POCSOLVERBUILDDIR"
  fi
AC_SUBST(POCSOLVERBUILDDIR,[${POCSOLVERBUILDDIR:="$d"}])
echo " set to $POCSOLVERBUILDDIR"

AC_ARG_VAR(POCSOLVERINC,[The poc-solver include directory. If not set, it automatically checks for, in that order :
$POCSOLVERBUILDDIR/include $POCSOLVERDIR/include $POCSOLVERDIR/src
])
if test "$POCSOLVERINC";then
  echo -n "POCSOLVERINC already"
else
  for d in $POCSOLVERBUILDDIR/include $POCSOLVERDIR/include $POCSOLVERDIR/src;do
    AC_MSG_CHECKING($d )
    if test -d "$d";then
      AC_MSG_RESULT(found)
      d=`$READLINK -f $d`
      break
    else
      AC_MSG_RESULT(does not exist)
      d=
      fi
    done
  if test -z "$d";then #not found
    echo "###################################################"
    echo "# FAILED TO FIND POC-SOLVER INCLUDE DIRECTORY !!! #"
    echo "###################################################"
    AC_MSG_ERROR([On this machine, you MUST give POCSOLVERBUILDDIR=... or POCSOLVERINC=... option to configure])
    fi
  echo -n "POCSOLVERINC"
  fi
AC_SUBST(POCSOLVERINC,[${POCSOLVERINC:="$d"}])
echo " set to $POCSOLVERINC"

AC_ARG_VAR(POCSOLVERLIB,[Where to find poc-solver library. If not set, it automatically checks the following list of directories, in that order :
$POCSOLVERBUILDDIR/lib $POCSOLVERDIR/lib
])
if test "$POCSOLVERLIB";then
  echo -n "POCSOLVERLIB already"
else
  for d in $POCSOLVERBUILDDIR/lib $POCSOLVERDIR/lib;do
    AC_MSG_CHECKING(poc-solvers.LDFLAGS in $d )
    if test -f "$d/poc-solvers.LDFLAGS";then
      AC_MSG_RESULT(found)
      d=`$READLINK -f $d`
      break
    else
      AC_MSG_RESULT(does not exist)
      d=
      fi
    done
  if test -z "$d";then #not found
    echo "#########################################"
    echo "# FAILED TO FIND POC-SOLVER LIBRARY !!! #"
    echo "#########################################"
    AC_MSG_ERROR([On this machine, you MUST give POCSOLVERBUILDDIR=... or POCSOLVERLIB=... option to configure])
    fi
  echo -n "POCSOLVERLIB"
  fi
AC_SUBST(POCSOLVERLIB,[${POCSOLVERLIB:="$d"}])
echo " set to $POCSOLVERLIB"
AC_SUBST(POCSOLVER_LDFLAGS,[${POCSOLVER_LDFLAGS:="-L$POCSOLVERLIB -lpoc-solver `cat $POCSOLVERLIB/poc-solvers.LDFLAGS`"}])
LIBS="$POCSOLVER_LDFLAGS $LIBS"


###################
# LAPACK detection
###################
#info:/autoconf/Generic%20Functions
missing=
while true;do
  AC_CHECK_FUNC([dgetrf_],
      [
      AC_DEFINE(LAPACK, 1, [take 1 value if LAPACK is used])
      break
      ]
    )
  AC_CHECK_FUNC([dgetrf],
      [
      AC_DEFINE(LAPACKC, 1, [take 1 value if LAPACK is used])
      break
      ]
    )
  #
  # LINPACK is deprecated ...
  #
  if test $missing;then AC_MSG_ERROR([No proper solver found]);break;fi
  missing=1
  unset ac_cv_func_dgetrf_
  unset ac_cv_func_dgetrf
  # See info:/autoconf/Libraries
  # Running it first without the ACTION-IF-FOUND for AC_CHECK_LIB to update LIBS
  AC_CHECK_LIB([lapack],[main])
  AC_CHECK_LIB([lapack],[main],
    [continue],
    [AC_MSG_WARN([Can not find lapack: install it or give LDFLAGS="\$LDFLAGS -L/lapack_PATH/lib" option to configure])]
    )
  AC_MSG_ERROR([No solver found])
  done


echo "############################"
echo "# checking compile options #"
echo "############################"

AC_ARG_WITH(sts,[AS_HELP_STRING([--with-sts],[sts, the spectral tidal solver. Enabled by default])])
AC_MSG_CHECKING([whether to compile sts])
if test "$with_sts" != no;then
  with_sts=yes
  fi
AM_CONDITIONAL([STS], [test "$with_sts" = yes])
AC_MSG_RESULT([$with_sts])

AC_ARG_WITH(ERROR_ON_EXPIRY,[AS_HELP_STRING([--with-ERROR_ON_EXPIRY],[enable ERROR_ON_EXPIRY to help the removal of expired code.])])
AC_MSG_CHECKING([whether to enable ERROR_ON_EXPIRY to help the removal of expired code])
AC_DEFINE([ERROR_ON_EXPIRY], [], [See expire macro for more information])
if test "$with_ERROR_ON_EXPIRY" != no;then
  with_ERROR_ON_EXPIRY=yes
  fi
if test "$with_ERROR_ON_EXPIRY" = yes;then
  AC_DEFINE([ERROR_ON_EXPIRY], [1])
else
  AC_DEFINE([ERROR_ON_EXPIRY], [0])
  fi
AC_MSG_RESULT([$with_ERROR_ON_EXPIRY])

AC_ARG_WITH(experimental,[AS_HELP_STRING([--with-experimental],[Experimental or ill-documented utilities. Disabled by default])])
AC_MSG_CHECKING([whether to compile experimental or ill-documented utilities])
if test "$with_experimental" != yes;then
  with_experimental=no
  fi
AM_CONDITIONAL([EXPERIMENTAL], [test "$with_experimental" = yes])
AC_MSG_RESULT([$with_experimental])


echo "#######################"
echo "# Checking endianness #"
echo "#######################"
AC_C_BIGENDIAN
case "$ac_cv_c_bigendian" in
no)
  need_swap=1
  need_ieeefp=0
  ;;
yes)
  need_swap=0
  need_ieeefp=1
  AC_CHECK_HEADER(ieeefp.h)
  ;;
*)
  ;;
  esac
AC_DEFINE_UNQUOTED(NEED_SWAP, [$need_swap], [take 1 value if swap is necessary])
AC_DEFINE_UNQUOTED(NEED_IEEEFP, [$need_ieeefp], [take 1 value if ieeefp.h is necessary])


echo "############################################"
echo "#"
echo "#            Final checks"
echo "#"
echo "############################################"

echo "Checks for typedefs, structures, and compiler characteristics."
AC_C_CONST
AC_TYPE_SIZE_T
AC_CHECK_TYPES([ptrdiff_t])

echo "Checks for library functions."
# ENABLING THE LINES BELOW WILL INDUCE ERRORS WITH gcc 4.2 !!!
#AC_FUNC_MALLOC
#AC_FUNC_REALLOC
AC_FUNC_STAT
AC_CHECK_FUNCS([floor sqrt strcasecmp strdup strrchr strstr])

echo ""
echo "LIBS='$LIBS'"
echo "POCSOLVER_LDFLAGS='$POCSOLVER_LDFLAGS'"


echo "##########################"
echo "# Creating the Makefiles #"
echo "##########################"
AC_CONFIG_FILES([Makefile
                objects/Makefile
                data/Makefile
                ])
AC_OUTPUT
