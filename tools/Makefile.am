
###########################################################
# USUAL STUFF

# Libtool stuff
# ACLOCAL_AMFLAGS= -I m4

#SUBDIRS = objects data
SUBDIRS = objects

# see info:/automake/Data
dist_noinst_DATA = Doxyfile

#overwriting DIST_COMMON to avoid having : decomp install-sh missing
DIST_COMMON = README $(am__configure_deps) $(dist_noinst_DATA) \
	$(srcdir)/Makefile.am COPYING INSTALL

###########################################################
# DEBUG THIS FILE. You can also use the -n make option.
debugMake:
	@printf "%s\n" $(foreach var,$(filter-out .VARIABLES am__remove_distdir,$(.VARIABLES)),'$(var)=$($(var))')
	printenv

###########################################################
# ENABLE STUFF LIKE make comodo-detidor FROM THIS DIRECTORY
BPF = bin_PROGRAMS.mk

objects/$(BPF) : objects/Makefile
	make -C objects $(BPF)

# using -include because it has a fail safe. See info:/make/Include
-include objects/$(BPF)

###########################################################
# PROMPTS FOR THE TARGET: VERY USEFULL WITH KDEVELOP

# using make rather than $(MAKE) to avoid propagating -j option
# unsetting MAKEFLAGS in the environment to avoid propagating -k option
ask:
	MAKEFLAGS=;f=.targetList~;t="$$(cat $$f)";set -x;\
	kdialog --inputbox "target ($$t)" -- "$$t">$$f~ && \
	  mv $$f~ $$f && eval "make -C objects `cat $$f`"

###########################################################
# make doxygen THAT CAUSES LESS PROBLEMS WITH KDEVELOP
DOD=../$(PACKAGE)-Doxygen

doxygen:
	{ cat Doxyfile;echo OUTPUT_DIRECTORY=$(DOD); } | doxygen -
	@if grep -qre '^ *GENERATE_HTML *= *YES' Doxyfile ;then\
	  echo -e "To browse, run:\n " xdg-open $$PWD/$(DOD)/html/index.html ;\
	else \
	  echo "Browse documentation from :" `$(READLINK) -f $(DOD)` ;\
	  fi

###########################################################
# ENABLE make bundles HGBD=../customBundleDir
# see the source for more information

# using -include to avoid automake's include. See info:/automake/Include
-include Makefile.hgftp
