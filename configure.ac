dnl --------------------------------
dnl Initialization macros.
dnl --------------------------------

AC_INIT([Wedge],[0.3])
AC_CONFIG_SRCDIR([wedge/wedgebase.h])
AC_CONFIG_HEADERS([wedge/config.h])
LT_INIT

AC_ARG_ENABLE([debug],
  [  --enable-debug    Turn on debugging (only for test program)],
  [case "${enableval}" in
       yes) debug=true ;;
       no)  debug=false ;;
       *) AC_MSG_ERROR([bad value ${enableval} for --enable-debug]) ;;
    esac],[debug=false])
AM_CONDITIONAL([DEBUG], [test x$debug = xtrue])

AC_ARG_ENABLE([arageli],
  [  --enable-arageli    Turn on support for Arageli],
  [case "${enableval}" in
       yes) arageli=true ;;
       no)  arageli=false ;;
       *) AC_MSG_ERROR([bad value ${enableval} for --enable-arageli]) ;;
    esac],[arageli=false])

dnl -----------------------------------------------
dnl Definition of the GiNaC version to use
dnl -----------------------------------------------

GINAC_MAJOR_VERSION=1
GINAC_MINOR_VERSION=6
GINAC_MICRO_VERSION=7

GINAC=ginac-$GINAC_MAJOR_VERSION.$GINAC_MINOR_VERSION.$GINAC_MICRO_VERSION
AC_SUBST(GINAC)

AC_MSG_CHECKING([rtti version])
if test $GINAC_MAJOR_VERSION -le 1 -a $GINAC_MINOR_VERSION -le 3 ; then
	AC_DEFINE([RTTI_METHOD], [0], [0 for integer tinfo, 1 for pointer tinfo, 2 for no tinfo])
	AC_MSG_RESULT([integer tinfo])
elif test $GINAC_MAJOR_VERSION -le 1 -a $GINAC_MINOR_VERSION -le 4 ; then
	AC_DEFINE([RTTI_METHOD], [1], [0 for integer tinfo, 1 for pointer tinfo, 2 for no tinfo])
	AC_MSG_RESULT([pointer tinfo])
else
	AC_DEFINE([RTTI_METHOD], [2], [0 for integer tinfo, 1 for pointer tinfo, 2 for no tinfo])
	AC_MSG_RESULT([no tinfo]) 
fi

dnl -----------------------------------------------
dnl Check whether "patch" is installed (needed by getginac script)
dnl -----------------------------------------------

AC_MSG_CHECKING([for patch])
AC_CHECK_PROG([PATCHPRESENT], [patch], [1], [0])
if test "$PATCHPRESENT" = "1" ; then 
	AC_MSG_RESULT(yes)
else
	AC_MSG_RESULT(no)
	AC_MSG_ERROR([Unable to run patch. Install patch (the program so called, not a specific patch) and try again.])
fi


dnl -----------------------------------------------
dnl Run script getginac.sh
dnl -----------------------------------------------

if ! test -f $GINAC/ginac/ncmul.cpp ; then
	AC_MSG_NOTICE([$GINAC not found; invoking script getginac.sh])
	$srcdir/getginac.sh $GINAC $srcdir
	if ! test $? -eq 0 ; then
		echo "Script getginac.sh failed. See README for instructions"
		exit 1;
	fi
fi


dnl -----------------------------------------------
dnl Package name and version number (user defined)
dnl -----------------------------------------------

WEDGE_LIBRARY_NAME=wedge

#release versioning
WEDGE_MAJOR_VERSION=0
WEDGE_MINOR_VERSION=3
WEDGE_MICRO_VERSION=0

#API version (often = WEDGE_MAJOR_VERSION.WEDGE_MINOR_VERSION)
WEDGE_API_VERSION=$WEDGE_MAJOR_VERSION.$WEDGE_MINOR_VERSION
AC_SUBST(WEDGE_API_VERSION)

#shared library versioning
WEDGE_LIBRARY_VERSION=5:0:0
#                       | | |
#                +------+ | +---+
#                |        |     |
#             current:revision:age
#                |        |     |
#                |        |     +- increment if interfaces have been added
#                |        |        set to zero if interfaces have been removed
#                                  or changed
#                |        +- increment if source code has changed
#                |           set to zero if current is incremented
#                +- increment if interfaces have been added, removed or changed


dnl --------------------------------
dnl Package name and version number
dnl --------------------------------

AC_SUBST(WEDGE_LIBRARY_VERSION)

PACKAGE=$WEDGE_LIBRARY_NAME
AC_SUBST(WEDGE_LIBRARY_NAME)

WEDGE_VERSION=$WEDGE_MAJOR_VERSION.$WEDGE_MINOR_VERSION.$WEDGE_MICRO_VERSION
WEDGE_RELEASE=$WEDGE_MAJOR_VERSION.$WEDGE_MINOR_VERSION
AC_SUBST(WEDGE_RELEASE)
AC_SUBST(WEDGE_VERSION)

VERSION=$WEDGE_VERSION

AM_INIT_AUTOMAKE($PACKAGE, $VERSION, no-define)

dnl -----------------------------------------------
dnl Checks for programs.
dnl -----------------------------------------------

AC_PROG_CC
AC_PROG_CXX
AM_PROG_LIBTOOL
AM_SANITY_CHECK
AC_CONFIG_SUBDIRS([$GINAC])
AC_LANG([C++])
AC_PATH_PROG(DOXYGEN, doxygen, "")
AM_CONDITIONAL(CONFIG_DOXYGEN, [test ! -z "$DOXYGEN"])

AC_CHECK_SIZEOF(int)
AC_CHECK_SIZEOF(long)
AC_CHECK_SIZEOF(long long)
AC_CHECK_SIZEOF(void *)

AC_COMPILE_IFELSE([AC_LANG_PROGRAM([#include <boost/spirit.hpp>],[])],[],AC_MSG_FAILURE([Boost C++ Libraries not installed]))
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([#include <boost/date_time/posix_time/posix_time_types.hpp>],[])],[],AC_MSG_FAILURE([Boost date-time C++ Libraries not installed]))

AC_MSG_CHECKING([for CoCoA])
LIBS_SAVED="$LIBS"
LIBS="$LIBS -lcocoa -lgmp -lcln -lboost_system"
export LIBS
AC_LINK_IFELSE([AC_LANG_PROGRAM([#include <CoCoA/library.H>],[CoCoA::GlobalManager CoCoAFoundations;])],[COCOA=true],[COCOA=false])
if test "$COCOA" = "true"; then
	AC_DEFINE([HAVE_COCOA], [], [Define to 1 if CoCoAlib is installed])
	AC_MSG_RESULT(yes)
else
	AC_MSG_RESULT(no)
fi
AM_CONDITIONAL(CONFIG_COCOA, [test x$COCOA = xtrue])
LIBS="$LIBS_SAVED"

if test "$arageli" = "true"; then
	AC_MSG_CHECKING([for Arageli])
	LIBS_SAVED="$LIBS"
	LIBS="$LIBS -larageli"
	export LIBS
	AC_LINK_IFELSE([AC_LANG_PROGRAM([#include <arageli/arageli.hpp>],[Arageli::matrix<float> m;])],[ARAGELI=true],[ARAGELI=false])
	if test "$ARAGELI" = "true"; then
		AC_DEFINE([HAVE_ARAGELI], [], [Define to 1 if Arageli is installed])
		AC_MSG_RESULT(yes)
	else
		AC_MSG_RESULT(no)
	fi
	LIBS="$LIBS_SAVED"
fi
AM_CONDITIONAL(CONFIG_ARAGELI, [test x$ARAGELI = xtrue])

AM_CXXFLAGS="-ffor-scope -Wctor-dtor-privacy -Wreorder -Wold-style-cast -Wsign-promo -Wchar-subscripts -Winit-self -Wmissing-braces -Wparentheses -Wreturn-type -Wswitch -Wtrigraphs -Wextra -Wno-sign-compare -DIN_GINAC -std=c++14 -Wno-deprecated-declarations"
AC_SUBST(AM_CXXFLAGS)

dnl -----------------------------------------------
dnl Generates Makefile's, configuration files and scripts
dnl -----------------------------------------------

AC_CONFIG_FILES(Makefile wedge-0.3.pc wedge/Makefile test/Makefile doc/Makefile doc/doxyfile examples/Makefile dictionary/Makefile dictionary/doc/Makefile)
AC_OUTPUT

