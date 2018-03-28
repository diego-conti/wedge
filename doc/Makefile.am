## Process this file with automake to produce Makefile.in

# Targets
if CONFIG_DOXYGEN
doc_DATA = html/index.html

html::
	@[ -d html ] || mkdir html >/dev/null 2>&1

html/index.html: 
	  @echo "Running ${DOXYGEN} doxyfile..."; \
	  ${DOXYGEN} doxyfile
endif
## CONFIG_DOXYGEN

install::

CLEANFILES = html/* 

EXTRA_DIST = doxyfile.in
