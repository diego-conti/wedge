## Process this file with automake to produce Makefile.in

# Targets
if CONFIG_DOXYGEN
doc_DATA = html/index.html

html::
	@[ -d html ] || mkdir html >/dev/null 2>&1

html/index.html:
	  @echo "Running ${DOXYGEN} dictionary.doxyfile..."; \
	  ${DOXYGEN} dictionary.doxyfile
endif
## CONFIG_DOXYGEN

CLEANFILES = html/* latex/*

EXTRA_DIST = dictionary.doxyfile
