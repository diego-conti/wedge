TESTS = testalgebra testmanifolds testbase testalgebraic testconnections testrepresentations testfunctions

if CONFIG_COCOA
LDADD = ../$(GINAC)/ginac/libginac.la ../wedge/libwedge.la -lcocoa -lcln -lgmp -lboost_system
else
LDADD = ../$(GINAC)/ginac/libginac.la ../wedge/libwedge.la -lcln -lgmp
endif


AM_CXXFLAGS=-ffor-scope -Wctor-dtor-privacy -Wreorder -Wsign-promo -Wchar-subscripts -Winit-self -Wmissing-braces -Wparentheses -Wreturn-type -Wswitch -Wtrigraphs -Wextra -Wno-sign-compare -DIN_GINAC -std=c++11 -Wno-deprecated-declarations


if DEBUG
cxxtestflags=--no-eh
conditional_recompile=debug_mode
AM_CPPFLAGS= -DDEBUG -I$(top_srcdir)/$(GINAC) -I$(top_srcdir)/$(GINAC)/ginac -I$(top_srcdir)/cxxtest -I$(top_srcdir)
else
conditional_recompile=normal_mode
AM_CPPFLAGS= -I$(top_srcdir)/$(GINAC) -I$(top_srcdir)/$(GINAC)/ginac -I$(top_srcdir)/cxxtest -I$(top_srcdir)
endif

if CONFIG_ARAGELI
TESTS+=aragelispeed
aragelispeed_SOURCES= aragelispeed.cpp 
aragelispeed_LDADD = ../$(GINAC)/ginac/libginac.la ../wedge/libwedge.la  -larageli
endif

testpoly_SOURCES= poly.h test.h
testalgebra_SOURCES= algebra.h test.h
testmanifolds_SOURCES= manifolds.h test.h
testbase_SOURCES= base.h test.h
testalgebraic_SOURCES= algebraic.h test.h
testconnections_SOURCES= connections.h test.h
testrepresentations_SOURCES = representations.h test.h
testfunctions_SOURCES = functions.h test.h

nodist_testalgebra_SOURCES=runneralgebra.cpp
nodist_testmanifolds_SOURCES=runnermanifolds.cpp
nodist_testbase_SOURCES=runnerbase.cpp
nodist_testalgebraic_SOURCES=runneralgebraic.cpp
nodist_testconnections_SOURCES=runnerconnections.cpp
nodist_testrepresentations_SOURCES = runnerrepresentations.cpp
nodist_testpoly_SOURCES=runnerpoly.cpp
nodist_testfunctions_SOURCES = runnerfunctions.cpp

testdictionary_LDADD = ../$(GINAC)/ginac/libginac.la ../wedge/libwedge.la  -lcln -lgmp $(top_srcdir)/dictionary/libdictionary.la

if CONFIG_COCOA
TESTS+=testpoly #testdictionary
testdictionary_SOURCES = dictionary.h test.h
nodist_testdictionary_SOURCES = runnerdictionary.cpp
check_PROGRAMS= $(TESTS)
endif

runneralgebra.cpp : algebra.h test.h $(conditional_recompile)
	$(top_srcdir)/cxxtest/cxxtestgen.pl  --error-printer -o runneralgebra.cpp algebra.h test.h $(cxxtestflags)

runnermanifolds.cpp : manifolds.h test.h $(conditional_recompile)
	$(top_srcdir)/cxxtest/cxxtestgen.pl  --error-printer -o runnermanifolds.cpp manifolds.h test.h $(cxxtestflags)

runnerbase.cpp : base.h test.h $(conditional_recompile)
	$(top_srcdir)/cxxtest/cxxtestgen.pl  --error-printer -o runnerbase.cpp base.h test.h $(cxxtestflags)

runneralgebraic.cpp : algebraic.h test.h $(conditional_recompile)
	$(top_srcdir)/cxxtest/cxxtestgen.pl  --error-printer -o runneralgebraic.cpp algebraic.h test.h $(cxxtestflags)

runnerconnections.cpp : connections.h test.h $(conditional_recompile)
	$(top_srcdir)/cxxtest/cxxtestgen.pl  --error-printer -o runnerconnections.cpp connections.h test.h $(cxxtestflags)

runnerrepresentations.cpp : representations.h test.h $(conditional_recompile)
	$(top_srcdir)/cxxtest/cxxtestgen.pl  --error-printer -o runnerrepresentations.cpp representations.h test.h $(cxxtestflags)

runnerpoly.cpp : poly.h test.h $(conditional_recompile)
	$(top_srcdir)/cxxtest/cxxtestgen.pl  --error-printer -o runnerpoly.cpp $(testpoly_SOURCES) $(cxxtestflags)
		
runnerdictionary.cpp : dictionary.h test.h $(conditional_recompile)
	$(top_srcdir)/cxxtest/cxxtestgen.pl  --error-printer -o runnerdictionary.cpp dictionary.h test.h $(cxxtestflags)

runnerfunctions.cpp : functions.h test.h $(conditional_recompile)
	$(top_srcdir)/cxxtest/cxxtestgen.pl  --error-printer -o runnerfunctions.cpp $(testfunctions_SOURCES) $(cxxtestflags)
	
debug_mode :
	touch debug_mode
	rm normal_mode -f

normal_mode :
	touch normal_mode
	rm debug_mode -f
	
