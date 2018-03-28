#!/bin/sh
if [ -z "$2" ]; then
	cd `dirname $0`
else
	cd $2
fi
if [ -z "$1" ]; then 
	GINAC=ginac-1.5.1
else
	GINAC=$1
fi
if [ ! -f $GINAC.tar.bz2 ]; then
	echo Cannot find $GINAC.tar.bz2; terminating
	exit 1
fi

bzip2 -d $GINAC.tar.bz2 -c | tar xf -
if [ ! $? -eq 0 ]; then
	exit 1
fi
patch -i Doxyfile.patch $GINAC/doc/reference/DoxyfileHTML.in


echo $GINAC | egrep "ginac-1.7" >/dev/null # check if GINAC 1.7.0 or later
if [ $? -eq 0 ]; then
  PRINT_PATCH=print11.patch
else
  PRINT_PATCH=print.patch
fi
#patch operators.cpp to make set_print_context and get_print_context available from Wedge
patch -i operators.patch $GINAC/ginac/operators.cpp
#patch print.h to make set_print_context and get_print_context available from Wedge
patch -i $PRINT_PATCH $GINAC/ginac/print.h

#patch function.hppy to make the latex name of a function available in Wedge
patch -i function_hppy.patch $GINAC/ginac/function.hppy
#patch function.cppy to make it possible to customize function output
patch -i function_cppy.patch $GINAC/ginac/function.cppy

#patch ncmul.cpp to make the inherited template class Wedge::Lambda work correctly.
#undo Chris Dams's "Made also ncmuls rename dummy indices" patch for a considerable performance gain.
echo $GINAC | egrep "ginac-1.7" >/dev/null # check if c++11, i.e. GINAC 1.6.7 or later
if [ $? -eq 0 ]; then
	patch -i ncmul11.patch $GINAC/ginac/ncmul.cpp
else
	patch -i ncmul.patch $GINAC/ginac/ncmul.cpp
  patch -i ncmul2.patch $GINAC/ginac/ncmul.cpp
fi

if [ ! $? -eq 0 ]; then
	rm $GINAC/ginac/ncmul.cpp
	exit 1
else
	exit 0
fi

