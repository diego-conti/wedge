## Source directory

h_sources = anylinalg.h adjoint.h expressions.h linear.h spinor.h aragelilinalg.h function.h logging.h structures.h basis.h ginaclinalg.h manifold.h tensor.h bilinear.h gl.h  manifoldwith.h tensorlambda.h cocoapolyalg.h gstructure.h parameters.h torsionfreeconnection.h concretemanifold.h lambda.h parse.h     utilities.h connection.h liegroupextension.h pforms.h vectorspace.h coordinates.h liegroup.h polybasis.h wedgealgebraic.h derivation.h liegroupstructures.h riemannianconnection.h wedgebase.h differentialform.h  liesubgroup.h riemannianstructure.h wexception.h liederivative.h linearcombinations.h affinebasis.h leibniz.h submersion.h config.h repgl.h repsl2.h repso.h so.h linearaction.h stabilizer.h su.h latex.h submersionwith.h transversestructure.h transverseconnection.h named.h printcontext.h fderivative.h omitfunctionargument.h pseudoriemannianstructure.h bilinearform.h pseudolevicivita.h simplifier.h normalform.h spiritsemanticactions.hpp

cc_sources = connection.cpp liegroup.cpp riemannianstructure.cpp wedgealgebraic.cpp differentialform.cpp logging.cpp spinor.cpp wexception.cpp ginaclinalg.cpp manifold.cpp structures.cpp parse.cpp utilities.cpp liesubgroup.cpp function.cpp submersion.cpp transverseconnection.cpp transversestructure.cpp linearaction.cpp latex.cpp omitfunctionargument.cpp fderivative.cpp simplifier.cpp normalform.cpp liederivative.cpp bilinearform.cpp pseudoriemannianstructure.cpp polybasis.cpp spiritgrammars.hpp

library_includedir=$(includedir)/$(WEDGE_LIBRARY_NAME)-$(WEDGE_API_VERSION)/$(WEDGE_LIBRARY_NAME)
library_include_HEADERS = $(h_sources)

AM_CPPFLAGS = -I$(top_srcdir)/$(GINAC) -I$(top_srcdir)/$(GINAC)/ginac -I$(top_srcdir)

lib_LTLIBRARIES= libwedge.la
libwedge_la_SOURCES= $(h_sources) $(cc_sources)
libwedge_la_LDFLAGS= -version-info $(WEDGE_LIBRARY_VERSION) -release $(WEDGE_RELEASE)

libwedge_la_LIBADD= ../$(GINAC)/ginac/libginac.la -lcln -lgmp

