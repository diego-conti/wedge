## top directory

#Build in these directories:
SUBDIRS= $(GINAC) wedge test doc dictionary examples 

#Distribute these directories:
DIST_SUBDIRS = wedge test doc dictionary examples 

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = wedge-0.2.pc
(pkgconfig_DATA): config.status

EXTRA_DIST= wedge-0.2.pc getginac.sh ncmul.patch ncmul2.patch ncmul11.patch matrix.patch Doxyfile.patch cxxtest function_hppy.patch function_cppy.patch operators.patch print.patch $(GINAC).tar.bz2

AM_CXXFLAGS= -ffor-scope -Wctor-dtor-privacy -Wreorder -Wold-style-cast -Wsign-promo -Wchar-subscripts -Winit-self -Wmissing-braces -Wmissing-include-dirs -Wparentheses -Wsequence-point -Wreturn-type -Wswitch -Wtrigraphs -Wextra -Wno-sign-compare -Wnon-virtual-dtor -fno-var-tracking-assignments -std=c++11 -Wno-deprecated-declarations

