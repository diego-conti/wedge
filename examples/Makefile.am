#Build in these directories:

#Distribute these directories:

noinst_PROGRAMS = clifford gkspinors curvature groups symplectic cartan almostcomplex estable so3 performance halfflat 

if CONFIG_COCOA
noinst_PROGRAMS +=  einsteinsasaki specialgeometries ricci dictionary
LDADD = ../$(GINAC)/ginac/libginac.la ../wedge/libwedge.la -lcocoa -lcln -lgmp -lboost_system
ricci_LDADD = ../$(GINAC)/ginac/libginac.la ../wedge/libwedge.la  ../dictionary/libdictionary.la  -lcln -lgmp -lboost_system
dictionary_LDADD = ../$(GINAC)/ginac/libginac.la ../wedge/libwedge.la  ../dictionary/libdictionary.la  -lcln -lgmp -lboost_system
else 
LDADD =  ../$(GINAC)/ginac/libginac.la ../wedge/libwedge.la  -lcln -lgmp
endif

so3_SOURCES = so3.cpp
estable_SOURCES = estable.cpp
cartan_SOURCES = cartan.cpp
symplectic_SOURCES = symplectic.cpp
clifford_SOURCES = clifford.cpp
gkspinors_SOURCES = gkspinors.cpp
almostcomplex_SOURCES = almostcomplex.cpp
einsteinsasaki_SOURCES = einsteinsasaki.cpp
curvature_SOURCES = curvature.cpp
groups_SOURCES = groups.cpp
specialgeometries_SOURCES = specialgeometries.cpp
performance_SOURCES = performance.cpp
halfflat_SOURCES= halfflat.cpp coherentsplittings.h
ricci_SOURCES = ricci.cpp parambasis.h
dictionary_SOURCES = dictionary.cpp

AM_CPPFLAGS = -I$(top_srcdir)/$(GINAC) -I$(top_srcdir)/$(GINAC)/ginac -I$(top_srcdir)


