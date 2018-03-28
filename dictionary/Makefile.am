AUTOMAKE_OPTIONS = subdir-objects
if CONFIG_COCOA 

h_sources =  dictionary.h  invariantforms.h  parsedictionary.h vvaluedforms.h basedictionary.h ode.h principal.h productmanifold.h

nobase_h_sources = dictionaries/SU3SO3Dictionary.h dictionaries/SU3U2Dictionary.h dictionaries/SU3SU2Dictionary.h dictionaries/SU3U2DictionaryWithWeight.h groups/SU3.h representations/R3.h representations/C2.h representations/C2_U2.h

cc_sources = dictionary.cpp invariantforms.cpp groups/groups.cpp	

AM_CPPFLAGS = -I$(top_srcdir)/$(GINAC) -I$(top_srcdir)/$(GINAC)/ginac -I$(top_srcdir) -fno-var-tracking-assignments

EXTRA_DIST = README doc/dictionary.doxyfile COPYING

SUBDIRS = doc

libdictionary_la_LIBADD = ../$(GINAC)/ginac/libginac.la ../wedge/libwedge.la -lcocoa -lcln -lgmp
libdictionary_la_LDFLAGS= -version-info 1:0:0 -release 3.0.0

lib_LTLIBRARIES= libdictionary.la
libdictionary_la_SOURCES= $(h_sources) $(cc_sources)
library_include_HEADERS = $(h_sources)
nobase_library_include_HEADERS = $(nobase_h_sources)
library_includedir=$(includedir)/$(WEDGE_LIBRARY_NAME)-$(WEDGE_API_VERSION)/dictionary

endif





