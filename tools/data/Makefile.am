#*******************************************************************************
# TeX FILES WITH THE OUTPUT OF <prog> --help
# grep's main source to search for commands
# i.e. for the argument of all the calls to \cmdhelp in ttb.tex
CMDS := $(shell $(SED) -re 's|.*\\cmdhelp\{([^}]*)\}.*|\1|;t;d' $(srcdir)/ttb.tex)
HELPTEXSUF := _-h.tex
EXECPATH= ../objects
# work out the list of such TeX files
# see info:/make/Text Functions
HELPTEXS := $(CMDS:=$(HELPTEXSUF))
# produce TeX files
%$(HELPTEXSUF): $(EXECPATH)/% Makefile
	PATH=$(EXECPATH)/ $* --help | $(SED) -e '1 {/^ *$$/ d}' >$@

#*******************************************************************************
# TeX FILES WITH THE OUTPUT OF <prog>
# grep's main source to search for the argument of all the calls to \cmdoutput in ttb.tex
# so far, only showarg is concerned
OUTPUTS := $(shell $(SED) -re 's|.*\\cmdoutput\{([^}]*)\}.*|\1|;t;d' $(srcdir)/ttb.tex)
OUTPUTTEXSUF := _.tex
# work out the list of such TeX files
# see info:/make/Text Functions
OUTPUTTEXS := $(OUTPUTS:=$(OUTPUTTEXSUF))
# produce TeX files
%$(OUTPUTTEXSUF): $(EXECPATH)/% Makefile
	PATH=$(EXECPATH)/ $* >$@

#*******************************************************************************
# TeX FILE WITH THE OUTPUT OF THE DESCRIPTION SECTION OF THE INLINE HELP
DESCTEX := descriptions.tex
$(DESCTEX): Makefile
# 	echo $(srcdir)/../src/*.cpp >$(DESCTEX)
	set -x ; for f in `grep -l print_help $(srcdir)/../src/*.cpp`; do \
	  echo '---------------------------------------------------------------------' ;\
	  echo $$f | $(SED) -re 's/.*\///;s/.cpp$$//' ;\
	  $(SED) -re '/DESCRIPTION\\n/,/\);|"(\\n)+"/! d;s/^ *"//;s/(DESCRIPTION)?\\n"(\);)?$$//;t;d' $$f ;\
	  done >$(DESCTEX)

#*******************************************************************************
ROOTTEXS = ttb.tex
$(abs_builddir)/literature.bib: $(srcdir)/literature.bib
	cp $(srcdir)/literature.bib $(abs_builddir)
# see info:/make/Text Functions
PDFS = $(ROOTTEXS:.tex=.pdf)
%.pdf: %.tex $(HELPTEXS) $(DESCTEX) $(abs_builddir)/literature.bib
	@echo commands: $(CMDS)
	@echo outputs of: $(OUTPUTS)
	set -x ; latex_count=5 ;\
	while [ $$latex_count -gt 0 ];do \
	  latex_count=$$[$$latex_count-1] ;\
	  pdflatex -synctex=1 -interaction=nonstopmode $(srcdir)/$* >/dev/null ;\
	  if grep "undefined citations" $*.log;then \
	    bibtex $*; \
	    pdflatex -synctex=1 -interaction=nonstopmode $(srcdir)/$* >/dev/null ;\
	    continue;fi;\
	  grep "Rerun to" $*.log || break ;\
	  done ;\
	touch $@

#*******************************************************************************
BUILT_SOURCES = $(HELPTEXS) $(OUTPUTTEXS)
dist_noinst_DATA = $(ROOTTEXS) literature.bib
pkgdata_DATA = $(PDFS)

#*******************************************************************************
clean:
	$(RM) cmdline.* *$(HELPTEXSUF) *$(OUTPUTTEXSUF)
	$(SHELL) -x -O extglob $(.SHELLFLAGS) '$(RM) ttb.!(tex|pdf)' |true
