
################################################################################
# FLAGS

# What is beween `@' is replaced by configure when producing Makefile's
POCSOLVERINC= @POCSOLVERINC@

POCSOLVER_LDFLAGS=@POCSOLVER_LDFLAGS@

#overwriting DIST_COMMON to avoid having : decomp
DIST_COMMON = $(include_HEADERS) $(noinst_HEADERS) \
        $(srcdir)/Makefile.am

CFLAGS += -fopenmp
CXXFLAGS = $(CFLAGS)
FCFLAGS = $(CFLAGS)
# see info:/automake/Yacc and Lex
AM_YFLAGS = -d

################################################################################
# HEADER FILES

# SEPARATE BUILD DIRECTORY TO ALLOW SEVERAL PEOPLE TO SHARE THE SAME SOURCE
SRC=$(top_srcdir)/src
abs_SRC=$(abs_top_srcdir)/src

# NOTE : DEFAULT_INCLUDES contains .
REVISION_FILE=revision.def
FLAGS_FILE=flags.def
BPF = bin_PROGRAMS.mk

# DO NOT MAKE A SHORTCUT VARIABLE TO `tools-fortran-sizes'
# WITHOUT FULLY TESTING IT FIRST WITH
#   make distclean;autoreconf -si && ./configure && make -j...
# AS IT SHOULD MESS UP THE DEPENDENCIES CALCULATIONS OF AUTOMAKE IN MOST PLACES!
FORTRAN_SIZES=$(SRC)/tools-fortran-sizes.def

#TODO : some sorting according to header type would be nice
#NOTE : any file listed in include_HEADERS will be installed, even if they are also listed in noinst_HEADERS!
include_HEADERS = \
$(SRC)/swap.h $(SRC)/grd.h $(SRC)/xyz.h $(SRC)/map.h $(SRC)/archive.h $(SRC)/xtrack-io.h \
  $(SRC)/poc-time.h $(SRC)/tides.h $(SRC)/admittance.h $(SRC)/admittance-mts.h $(SRC)/spectrum.h\
  $(SRC)/fe.h $(SRC)/fe.def $(SRC)/fe-classes.h $(SRC)/fe-proto.h $(SRC)/fe-integration.h \
  $(SRC)/mass.h $(SRC)/statistic.h $(SRC)/altimetry.h \
  $(SRC)/polygones.h $(SRC)/polygons.hpp $(SRC)/polygons_implementation.hpp $(SRC)/polygons_exceptions.hpp $(SRC)/constants.h $(SRC)/exceptions.hpp $(SRC)/maths.h\
  $(SRC)/filter.h $(SRC)/fourier.h $(SRC)/loess.h $(SRC)/lanczos.h \
  $(SRC)/ascii.h $(SRC)/rutin.h $(SRC)/topo.h \
  $(SRC)/tools-structures.h $(SRC)/tools-define.h $(SRC)/poc-assertions.h \
  $(SRC)/bmg.h $(SRC)/geo.h $(SRC)/gshhs.h $(SRC)/sym-io.h\
  $(SRC)/mgr.h $(SRC)/mgr-converter.h   $(SRC)/list.h $(SRC)/zapper.h\
  $(SRC)/poc-netcdf-assertions.h $(SRC)/netcdf-classes.h $(SRC)/netcdf-proto.h \
  $(SRC)/poc-netcdf.hpp $(SRC)/poc-netcdf-io.h $(SRC)/poc-netcdf-data.hpp \
  $(SRC)/poc-data-operators.hpp \
  $(SRC)/poc-grib.h \
  $(SRC)/legend.h $(SRC)/functions.h $(SRC)/poc-list.hpp \
  $(SRC)/discretisation.h $(SRC)/triangle.h \
  $(SRC)/cefmo.h $(SRC)/quadrangle.h $(SRC)/gauss.h $(SRC)/periodic.h \
  $(SRC)/parallel.h $(SRC)/datastream.h $(SRC)/harmonic.h $(SRC)/academic.h \
  $(SRC)/meteo.h \
  $(FORTRAN_SIZES)

# see info:/automake/Built%20Sources%20Example
BUILT_SOURCES = $(REVISION_FILE) $(FLAGS_FILE) $(BPF) $(FORTRAN_SIZES) \
  poc-formula-parse.@AM_YPP_TARGET@
# see info:/make/Force%20Targets
REVISION_FILE_PREREQUISITES = $(shell ls $(abs_top_srcdir)/.hg/branch || echo FORCE)
$(REVISION_FILE): $(REVISION_FILE_PREREQUISITES)
# creating a temporary file
#  so that the file is only updated if and only if necessary
#   so that targets are remade only when is necessary
	which hg \
	  && hg parent -R $(abs_top_srcdir) --template '#define HG_REV "{rev}:{node|short}"\n#define REVISION "Mercurial revision " HG_REV " of {date|isodate}"\n' > $(REVISION_FILE)~ \
	  || touch $(REVISION_FILE)~
# do not overwrite a non-empty target with an empty temporary file
# and do not touch the target if there is no difference with the temporary file
	test -s $(REVISION_FILE) -a "!" -s $(REVISION_FILE)~ \
	  || diff -q $(REVISION_FILE)~ $(REVISION_FILE) \
	  || mv -v $(REVISION_FILE)~ $(REVISION_FILE)
FORCE:

#this is most necessary to avoid crashes in parallel mode
tools-fortran-sizes.$(OBJEXT): $(REVISION_FILE)

$(FORTRAN_SIZES): tools-fortran-sizes
	./tools-fortran-sizes >$(FORTRAN_SIZES)

#generation des fichiers d'include prive
#######################################
# *** USEFUL FOR THE TARBALLS !!! *** #
#######################################
noinst_HEADERS = $(SRC)/archive.def \
  $(REVISION_FILE) $(SRC)/version-macros.def \
                 $(SRC)/loess.def $(SRC)/legend.def $(SRC)/poc-netcdf.def $(SRC)/fe.def $(SRC)/map.def \
                 $(SRC)/tides.def $(SRC)/sym-io.def $(SRC)/matrix.h $(SRC)/mss.v2.h

################################################################################
# AM FLAGS

#including $(SRC) for (inexisting) examples and dependencies via tools-config
AM_CPPFLAGS := -I$(SRC) -I$(POCSOLVERINC)
# TO CREATE PORTABLE CODE, ENABLE THE OPTIONS BELOW 
# AM_CFLAGS = -ansi -pedantic -Wall
# AM_CLAGS = -std=c99 -pedantic -W -Wall

# Compilers whose version do not recognise an option will ONLY output a warning.
# PLEASE NOTE THAT A WARNING IS NOT AN ERROR !!!!
AM_CFLAGS = -Wuninitialized -Wno-unused-but-set-variable
# To also enable warnings for unused variables, add the argument on the line below to the make command.
# "AM_CFLAGS = -Wunused -Wuninitialized -Wno-unused-but-set-variable"
# -Wall => -Wunused 
AM_CXXFLAGS = $(AM_CFLAGS)

################################################################################
# LIBRARIES

# lib_LIBRARIES = libgeo.a libtime.a libmap.a librutin.a \
#   libpoc-netcdf.a libfilter.a libtide.a libgrd.a liblist.a\
#   libplg.a libarchive.a libfe.a libsympho.a \
#   libmgr.a liblegend.a\
#   libsparse.a

libgeo_a_SOURCES   = $(SRC)/geo.cpp $(SRC)/geo-02.cpp
libtime_a_SOURCES  = $(SRC)/poc-time.cpp
libmap_a_SOURCES   = $(SRC)/map-utilities.cpp $(SRC)/map-interpolation.cpp \
  $(SRC)/map-indint.cpp $(SRC)/map-operators.cpp $(SRC)/zones.cpp \
  $(SRC)/map-slices.cpp $(SRC)/map03.cpp $(SRC)/map04.cpp $(SRC)/map05.cpp $(SRC)/map06.cpp \
  $(SRC)/maths00.cpp $(SRC)/maths01.cpp $(SRC)/bmg01.cpp $(SRC)/ascii.cpp $(SRC)/topo01.cpp $(SRC)/topo-filter.cpp

librutin_a_SOURCES = \
  $(SRC)/rutin.cpp $(SRC)/functions.cpp $(SRC)/poc-assertions.cpp \
  $(SRC)/poc-procs.cpp $(SRC)/poc-string.cpp $(SRC)/poc-list.cpp \
  $(SRC)/swap.cpp $(SRC)/statistic.cpp

libold_netcdf_a_SOURCES = $(SRC)/netcdf-purgatory.cpp $(SRC)/netcdf-parse.cpp $(SRC)/netcdf-utilities.cpp $(SRC)/netcdf-vio.cpp \
                          $(SRC)/netcdf-standards-01.cpp $(SRC)/netcdf-standards-02.cpp $(SRC)/netcdf-structured.cpp $(SRC)/netcdf-unstructured.cpp $(SRC)/netcdf-parallel.cpp\
                          $(SRC)/netcdf-time.cpp         $(SRC)/netcdf-io.cpp      $(SRC)/netcdf-genesis.cpp
libpoc_netcdf_a_SOURCES = \
  $(SRC)/poc-netcdf-input.cpp $(SRC)/poc-netcdf-output.cpp \
  $(SRC)/poc-netcdf-iio.cpp $(SRC)/poc-netcdf-comodo-parsing.cpp $(SRC)/poc-netcdf-comodo.cpp $(SRC)/poc-netcdf-comodo-iio.cpp
poc_grib_SOURCES = \
  $(SRC)/poc-grib.cpp

libfilter_a_SOURCES =  $(SRC)/fourier.cpp $(SRC)/loess.cpp $(SRC)/lanczos.cpp \
  $(SRC)/filter.cpp $(SRC)/filter-02.cpp $(SRC)/harmonic-filter.cpp

libtide_a_SOURCES   =  $(SRC)/admittance.cpp $(SRC)/admittance-mts.cpp \
  $(SRC)/tides.cpp $(SRC)/analysis.cpp \
  $(SRC)/tides-astronomic.cpp  $(SRC)/tides-io.cpp $(SRC)/tides01.cpp $(SRC)/tides-friction.cpp \
  $(SRC)/harmonic_prediction.cpp $(SRC)/woce.cpp $(SRC)/harmonic_t.cpp $(SRC)/harmonic-01.cpp $(SRC)/harmonic-02.cpp $(SRC)/spectrum.cpp

libgrd_a_SOURCES     = $(SRC)/grd.cpp $(SRC)/xyz01.cpp $(SRC)/xyz02.cpp
liblist_a_SOURCES    = $(SRC)/list.cpp
libplg_a_SOURCES     = $(SRC)/polygones-init.cpp $(SRC)/polygones-io.cpp \
  $(SRC)/polygones-search.cpp $(SRC)/polygones-intersect.cpp \
  $(SRC)/polygones-contours.cpp $(SRC)/polygones-sampling.cpp \
  $(SRC)/polygones-extraction.cpp $(SRC)/polygones-edition.cpp \
  $(SRC)/polygones-processing.cpp \
  $(SRC)/polygones-03.cpp $(SRC)/polygones-04.cpp $(SRC)/polygones-05.cpp
libarchive_a_SOURCES = $(SRC)/archive.cpp $(SRC)/archive-clx.cpp $(SRC)/archive-cdf.cpp $(SRC)/io-obsolete.cpp

libfe_a_SOURCES = $(SRC)/fe-integration-optimisation.cpp $(SRC)/fe-integration-triangle.cpp \
  $(SRC)/fe-integration-gauss.cpp $(SRC)/fe-integration-quadrangle.cpp \
  $(SRC)/fe01.cpp $(SRC)/fe02.cpp $(SRC)/fe-search.cpp $(SRC)/fe-slices.cpp \
  $(SRC)/fe-tables.cpp $(SRC)/fe-topology.cpp $(SRC)/fe-edition.cpp $(SRC)/fe-edition-vertices.cpp \
  $(SRC)/fe05.cpp $(SRC)/fe06.cpp $(SRC)/fe-partition.cpp $(SRC)/fe08.cpp\
  $(SRC)/fe-netcdf-01.cpp $(SRC)/fe-netcdf-02.cpp $(SRC)/fe-netcdf3D.cpp \
  $(SRC)/fe-io.cpp $(SRC)/fe-save.cpp \
  $(SRC)/fe-matrix-LGP0.cpp $(SRC)/fe-matrix-LGP1.cpp            $(SRC)/fe-matrix-DGP1.cpp            $(SRC)/fe-matrix-NCP1.cpp\
  $(SRC)/fe-matrix-LGP2.cpp $(SRC)/fe-matrix-DGP2.cpp            $(SRC)/fe-matrix-QP0.cpp\
  $(SRC)/triangle.cpp       $(SRC)/discretisation-initialise.cpp $(SRC)/discretisation-utilities.cpp\
  $(SRC)/mass.cpp           $(SRC)/fe-delauney.cpp               $(SRC)/fe-meshing.cpp                $(SRC)/fe-critere.cpp                $(SRC)/fe-remapping.cpp \
  $(SRC)/fe-gradients.cpp \
  $(SRC)/quadrangle-00.cpp  $(SRC)/quadrangle-01.cpp $(SRC)/quadrangle-io.cpp $(SRC)/periodic.cpp

libsympho_a_SOURCES =  $(SRC)/sym-io.cpp

libmgr_a_SOURCES = \
  $(SRC)/mgr.cpp $(SRC)/mgr02.cpp $(SRC)/mgr-io.cpp $(SRC)/mgr-woce.cpp \
  $(SRC)/mgr-converter01.cpp $(SRC)/mgr-loaders.cpp \
  $(SRC)/mgr-moorings-stations.cpp

liblegend_a_SOURCES =  $(SRC)/legend01.cpp $(SRC)/legend03.cpp

libsparse_a_SOURCES =  $(SRC)/solvers.cpp $(SRC)/linear.cpp $(SRC)/linear-z.cpp $(SRC)/linear-utilities.cpp $(SRC)/matrix.cpp

LIBSOURCES = $(SRC)/fe-ascii.cpp $(libgeo_a_SOURCES) $(libtime_a_SOURCES)\
  $(SRC)/map-interpolation-1D.cpp \
             $(libmap_a_SOURCES)\
             $(librutin_a_SOURCES)\
  $(libarchive_a_SOURCES) $(libold_netcdf_a_SOURCES) $(libpoc_netcdf_a_SOURCES)  \
  $(poc_grib_SOURCES) \
  $(liblegend_a_SOURCES)\
             $(libfilter_a_SOURCES)\
             $(libtide_a_SOURCES)\
             $(SRC)/sts_fortran.cpp\
             $(libgrd_a_SOURCES)\
             $(liblist_a_SOURCES) $(libplg_a_SOURCES) $(libfe_a_SOURCES) $(libsparse_a_SOURCES)\
             $(libsympho_a_SOURCES) $(libmgr_a_SOURCES) $(libmgr_converter_a_SOURCES) $(libmgr_moorings_stations_a_SOURCES) \
  $(SRC)/datastream.cpp \
  $(SRC)/poc-data-operators.cpp \
  poc-formula-scan.lpp poc-formula-parse.ypp \
  $(SRC)/poc-formula.cpp $(SRC)/poc-formula-functions.cpp $(SRC)/map-editors.cpp \
  $(SRC)/poc-formula-vio.cpp

# Last paragraph of info:/automake/Linking
# We recommend that you avoid using `-l' options in `LDADD' or `PROG_LDADD' when referring to libraries built by your package. Instead, write the file name of the library explicitly (...). Use `-l' only to list third-party libraries. If you follow this rule, the default value of `PROG_DEPENDENCIES' will list all your local libraries and omit the other ones.

# EITHER
# # FAILS BECAUSE libpoc-solver IS STATIC
# # see info:/libtool/Using Automake
# lib_LTLIBRARIES = libtools.la
# libtools_la_SOURCES = $(LIBSOURCES)
# libtools_la_DEPENDENCIES = $(REVISION_FILE)
# LDADD = libtools.la
# OR
lib_LIBRARIES = libtools.a lib4tugo.a

libtools_a_SOURCES = $(LIBSOURCES)
libtools_a_DEPENDENCIES = $(REVISION_FILE)
LDADD = libtools.a

lib4tugo_a_SOURCES = $(librutin_a_SOURCES) $(libpoc_netcdf_a_SOURCES) \
  $(SRC)/map-utilities.cpp $(SRC)/map-operators.cpp \
  $(SRC)/map-interpolation-1D.cpp $(SRC)/map-indint.cpp $(SRC)/map-interpolation.cpp \
  $(SRC)/geo.cpp $(SRC)/filter.cpp $(SRC)/lanczos.cpp $(SRC)/loess.cpp $(SRC)/fourier.cpp \
  $(SRC)/poc-data-operators.cpp $(SRC)/poc-formula-scan.lpp $(SRC)/poc-formula-parse.ypp \
  $(SRC)/poc-formula.cpp $(SRC)/poc-formula-functions.cpp $(SRC)/map-editors.cpp \
  $(SRC)/poc-formula-vio.cpp
lib4tugo_a_DEPENDENCIES = Makefile
# EO EITHER/OR


################################################################################
# EXECUTABLES

bin_PROGRAMS = tools-config\
  delauney fe-format fe-operation\
  geo-toponyms topo-create topo-merge topo-ensemble topo-assembly topo-compare topo-antarctic topo-format topo_smooth topo-operation\
  mesh-partition mesh-assembly mesh-refine mesh-upgrade mesh-split mesh-topo mesh_correlation mesh-editor\
  mesh-doctor mesh-diagnostics mesh-format mesh-renum mesh-academic mesh-generator mesh-limits mesh-zmin mesh-reshape\
  mesh-continents\
  mgr-converter mgr-format mgr-editing mgr-predictor\
  variability extract-nf-v1 extract-nf-v2 extract-nf-v3 extract-sections extract-vprofile turbulence\
  kelvin tides-energy tides-compare tides-merge tides-sample tides-converter tides-validate tides-diagnostics vertical-eigenmodes\
  rectify mesh2mesh mesh2node mesh2restart2d gridit-cdf gridit-cls gridit-versatile gridit-surfref search_nodes \
  shom-grd2nc cm93toXYZ \
  synthetic operation sectionxz sectionxy regular2node3d \
  convert3d statexyz analysis-tugo analysis-ecoop timeseries-validate haw2gnu format symsolve \
  sample regular2node mog3drestart reduce-v1 reduce-v2 meteo-filter \
  convert-nf grid-archives excalibur gridit getharmo detidor adcp SBprofiler clock spectral\
  symtools symharmo symtides grid-config grid-parameter\
  merge serie interpolation-slice interpolator interpolator-cdf fix beta mssV2\
  ww3tools\
  shorelines-doctor shorelines-converter shorelines-smoother shorelines-extractor shorelines-merge shorelines-distance\
  magic cto statistic1-nf statistic2-nf statistic3-nf topex detide-nf restart-nf\
  showarg gdr.v2.1 mss stream-2cdf\
  fe2obc mgr2obc atlas2obc fes-assembly fes-converter symmic xover\
  mapxyz-v0 mapxyz-v1 mapxyz-unstructured\
  gmt-converter altimetry-error altimetry-decimate altimetry-filter\
  penalty xtrack-format altimetry-detidor mercator-detidor comodo-detidor radar-detidor\
  comodo-filter metocean comodo-control comodo-compliance comodo-sfilter concatenator\
  RDI-converter\
  fft-test\
  ww3-U2S academics predictor xtt grid-track metocean-predictor topo-detidor \
  comodo-admittance comodo-energy comodo-calculator ellipse ensemble-extract\
  loading omptest\
  tools-fortran-sizes


#-------------------------------------------------------------------------------
# compilation-related

$(BPF): Makefile
	printf "%s\n" $(bin_PROGRAMS)|$(SED) -re 's|^.+$$|\0:\n\t$$(MAKE) -C objects $(REVISION_FILE) $(FLAGS_FILE) \0\n|' >$(BPF)

# Using an executable rather than a pkg-config file because it can detect
# whether it is installed or not and output accordingly

$(FLAGS_FILE): Makefile
	printf "%s\n" \
	  '#define TOOLS_BUILDDIR "'"$(abs_builddir)"'"' \
	  '#define TOOLS_BINDIR "'"$(bindir)"'"' \
	  '' \
	  '#define TOOLS_LDFLAGS "'"$(POCSOLVER_LDFLAGS)"'"' \
	  '' \
	  '#define TOOLS_INCFLAGS "'"$(AM_CPPFLAGS:-I$(SRC)=-I$(abs_SRC))"'"' \
	  > $(FLAGS_FILE)

tools_config_SOURCES = $(SRC)/tools-config.cpp

tools_fortran_sizes_SOURCES = $(FORTRAN_SIZES:def=cpp)

#-------------------------------------------------------------------------------

delauney_SOURCES = $(SRC)/delauney.cpp

xover_SOURCES = $(SRC)/xover.cpp

altimetry_error_SOURCES = $(SRC)/altimetry.cpp $(SRC)/altimetry-error.cpp

altimetry_decimate_SOURCES = $(SRC)/altimetry.cpp $(SRC)/altimetry-decimate.cpp

altimetry_filter_SOURCES = $(SRC)/altimetry.cpp $(SRC)/altimetry-filter.cpp

penalty_SOURCES = $(SRC)/penalty.cpp

##-----------------------------------------------------------------------------
# unstructured grid utilities
mesh_split_SOURCES = $(SRC)/mesh-split.cpp

mesh_upgrade_SOURCES = $(SRC)/mesh-upgrade.cpp

mesh_reshape_SOURCES = $(SRC)/mesh-reshape.cpp

mesh_refine_SOURCES = $(SRC)/mesh-refine.cpp

mesh_editor_SOURCES = $(SRC)/mesh-editor.cpp

mesh_assembly_SOURCES = $(SRC)/mesh-assembly.cpp

mesh_partition_SOURCES = $(SRC)/mesh-partition.cpp

mesh_topo_SOURCES = $(SRC)/mesh-topo.cpp

mesh_zmin_SOURCES = $(SRC)/mesh-zmin.cpp

mesh_correlation_SOURCES = $(SRC)/mesh-correlation.cpp

mesh_limits_SOURCES = $(SRC)/mesh-limits.cpp

mesh_generator_SOURCES = $(SRC)/mesh-generator.cpp

mesh_continents_SOURCES = $(SRC)/mesh-continents.cpp

mesh_academic_SOURCES =  $(SRC)/academic-estuary.cpp $(SRC)/mesh-academic.cpp 

mesh_diagnostics_SOURCES = $(SRC)/mesh-diagnostics.cpp

mesh_doctor_SOURCES = $(SRC)/mesh-doctor.cpp

mesh_renum_SOURCES = $(SRC)/mesh-renum.cpp

mesh_format_SOURCES = $(SRC)/mesh-format.cpp

fes_assembly_SOURCES = $(SRC)/tugo-obc.cpp $(SRC)/fes-assembly.cpp

fes_converter_SOURCES = $(SRC)/fes-converter.cpp

mesh2node_SOURCES = $(SRC)/mesh2node.cpp

mesh2mesh_SOURCES = $(SRC)/mesh2mesh.cpp

mesh2restart2d_SOURCES = $(SRC)/mesh2restart2d.cpp

fe_format_SOURCES = $(SRC)/surfref-io.cpp $(SRC)/fe-format.cpp

fe_operation_SOURCES = $(SRC)/surfref-io.cpp $(SRC)/fe-operation.cpp

##-----------------------------------------------------------------------------
# tides utilities
tides_converter_SOURCES = $(SRC)/tides-converter.cpp

tides_compare_SOURCES = $(SRC)/tides-compare.cpp

tides_merge_SOURCES = $(SRC)/tides-merge.cpp

tides_energy_SOURCES = $(SRC)/tides03.cpp $(SRC)/spectral-inversion.cpp $(SRC)/spectral-energy01.cpp $(SRC)/spectral-energy02.cpp $(SRC)/tides-energy.cpp

extract_vprofile_SOURCES =$(SRC)/extract-vprofile.cpp

turbulence_SOURCES =$(SRC)/turbulence.cpp

kelvin_SOURCES = $(SRC)/tugo-obc.cpp $(SRC)/kelvin.cpp

tides_sample_SOURCES = $(SRC)/tides-sample.cpp

altimetry_detidor_SOURCES = $(SRC)/altimetry-detidor.cpp $(SRC)/xtrack-io.cpp $(SRC)/xtrack-io.obsolete.cpp 

mercator_detidor_SOURCES = $(SRC)/mercator-detidor.cpp

comodo_detidor_SOURCES = $(SRC)/comodo-detidor.cpp

radar_detidor_SOURCES = $(SRC)/radar-detidor.cpp

comodo_control_SOURCES = $(SRC)/comodo-control.cpp

comodo_compliance_SOURCES = $(SRC)/comodo-compliance.cpp

fft_test_SOURCES = $(SRC)/fft-test.cpp

comodo_filter_SOURCES = $(SRC)/comodo-filter.cpp

metocean_SOURCES = $(SRC)/metocean.cpp

loading_SOURCES = $(SRC)/loading.cpp
# check with -fopt-info-optimized-missed=loading.opt_info, then
# > grep VECTORIZED loading.opt_info|sort -t : -k 2 -gu
# Making sure this is ALWAYS optimised
#loading_CXXFLAGS = $(CXXFLAGS) -O3

omptest_SOURCES = $(SRC)/omptest.cpp
omptest_CXXFLAGS = $(CXXFLAGS) -O3

predictor_SOURCES = $(SRC)/predictor.cpp

xtt_SOURCES = $(SRC)/xtt.cpp

grid_track_SOURCES = $(SRC)/grid-track.cpp

mgr_predictor_SOURCES = $(SRC)/mgr-predictor.cpp

metocean_predictor_SOURCES = $(SRC)/metocean-predictor.cpp

topo_detidor_SOURCES = $(SRC)/topo-detidor.cpp

comodo_admittance_SOURCES = $(SRC)/comodo-admittance.cpp

comodo_energy_SOURCES = $(SRC)/comodo-energy.cpp

comodo_calculator_SOURCES = $(SRC)/comodo-calculator.cpp

ellipse_SOURCES = $(SRC)/ellipse.cpp

comodo_sfilter_SOURCES = $(SRC)/comodo-sfilter.cpp

concatenator_SOURCES = $(SRC)/concatenator.cpp

fe2obc_SOURCES = $(SRC)/tugo-obc.cpp $(SRC)/fe2obc.cpp

mgr2obc_SOURCES = $(SRC)/tugo-obc.cpp $(SRC)/mgr2obc.cpp

atlas2obc_SOURCES = $(SRC)/tugo-obc.cpp $(SRC)/atlas2obc.cpp

tides_validate_SOURCES = $(SRC)/tides-validate.cpp

ensemble_extract_SOURCES = $(SRC)/ensemble-extract.cpp

tides_diagnostics_SOURCES = $(SRC)/crozet.cpp $(SRC)/eos.cpp $(SRC)/density.cpp $(SRC)/sturm-liouville.cpp $(SRC)/tides-diagnostics.cpp

vertical_eigenmodes_SOURCES = $(SRC)/eos.cpp $(SRC)/nemo-api.cpp $(SRC)/density.cpp $(SRC)/sturm-liouville.cpp $(SRC)/vertical-modes-io.cpp $(SRC)/vertical-decomposition.cpp $(SRC)/vertical-energy.cpp $(SRC)/vertical-eigenmodes.cpp

rectify_SOURCES = $(SRC)/rectify.cpp

mapxyz_v0_SOURCES = $(SRC)/mapxyz-v0.cpp
mapxyz_v1_SOURCES =  $(SRC)/mapxyz-v1.cpp

mapxyz_unstructured_SOURCES =  $(SRC)/mapxyz-unstructured.cpp

academics_SOURCES = $(SRC)/academics.cpp $(SRC)/xtrack-io.cpp

##-----------------------------------------------------------------------------
# topography utilities
geo_toponyms_SOURCES = $(SRC)/geo-toponyms.cpp

topo_create_SOURCES = $(SRC)/topo-create.cpp

topo_merge_SOURCES = $(SRC)/topo-merge.cpp

topo_ensemble_SOURCES = $(SRC)/topo-ensemble.cpp

topo_assembly_SOURCES = $(SRC)/topo-assembly.cpp

topo_antarctic_SOURCES = $(SRC)/topo-antarctic.cpp

topo_format_SOURCES = $(SRC)/topo-format.cpp

topo_smooth_SOURCES = $(SRC)/topo-smooth.cpp

topo_compare_SOURCES = $(SRC)/topo-compare.cpp

topo_operation_SOURCES = $(SRC)/topo-operation.cpp

gmt_converter_SOURCES = $(SRC)/gmt-converter.cpp

xtrack_format_SOURCES = $(SRC)/xtrack-io.cpp $(SRC)/xtrack-format.cpp

##-----------------------------------------------------------------------------
# UG to SG utilities

gridit_SOURCES = $(SRC)/gridit.cpp

gridit_surfref_SOURCES = $(SRC)/surfref-io.cpp $(SRC)/gridit-surfref.cpp

gridit_cdf_SOURCES = $(SRC)/gridit-cdf.cpp

gridit_versatile_SOURCES = $(SRC)/gridit-versatile.cpp

gridit_cls_SOURCES = $(SRC)/gridit-cls.cpp

shom_grd2nc_SOURCES = $(SRC)/shom-grd2nc.cpp

cm93toXYZ_SOURCES = $(SRC)/cm93toXYZ.cpp

search_nodes_SOURCES = $(SRC)/search_nodes.cpp

synthetic_SOURCES = $(SRC)/sigma.cpp $(SRC)/quaker.cpp $(SRC)/synthetic.cpp

operation_SOURCES = $(SRC)/operation.cpp

sectionxz_SOURCES = $(SRC)/sectionxz.cpp

sectionxy_SOURCES = $(SRC)/sectionxy.cpp

regular2node_SOURCES = $(SRC)/regular2node.cpp

regular2node3d_SOURCES = $(SRC)/regular2node-3d.cpp

convert3d_SOURCES = $(SRC)/convert3d.cpp

statexyz_SOURCES = $(SRC)/statexyz.cpp

mss_SOURCES = $(SRC)/mss-io.cpp $(SRC)/mss.v1.cpp

mssV2_SOURCES = $(SRC)/mss-comp.v2.cpp $(SRC)/mss-io.v2.cpp $(SRC)/mss-main.v2.cpp

sample_SOURCES = $(SRC)/sample.cpp

format_SOURCES = $(SRC)/format.cpp

analysis_tugo_SOURCES = $(SRC)/analysis-tugo.cpp

analysis_ecoop_SOURCES = $(SRC)/analysis-ecoop.cpp

haw2gnu_SOURCES = $(SRC)/haw2gnu.cpp

mog3drestart_SOURCES = $(SRC)/mog3d-restart.cpp

reduce_v1_SOURCES = $(SRC)/reduce-v1.cpp

reduce_v2_SOURCES = $(SRC)/reduce-v2.cpp

meteo_filter_SOURCES = $(SRC)/meteo-filter.cpp

convert_nf_SOURCES = $(SRC)/convert-nf.cpp

grid_archives_SOURCES = $(SRC)/archive-clx.cpp $(SRC)/grid-archives.cpp

##-----------------------------------------------------------------------------
# symphonie utilities
symtools_SOURCES = $(SRC)/symtools.cpp

symtides_SOURCES = $(SRC)/symtides.cpp

symharmo_SOURCES = $(SRC)/symharmo.cpp

symsolve_SOURCES = $(SRC)/symsolve.cpp

symmic_SOURCES = $(SRC)/quaker.cpp $(SRC)/symmic.cpp 

ww3_U2S_SOURCES = $(SRC)/quaker.cpp $(SRC)/ww3-U2S.cpp 

grid_parameter_SOURCES = $(SRC)/grid-parameter.cpp

grid_config_SOURCES = $(SRC)/grid-config.cpp

##-----------------------------------------------------------------------------
# WW3 utilities

ww3tools_SOURCES = $(SRC)/ww3-io.cpp $(SRC)/ww3tools.cpp

#reswap il manque la fonction setswap ???
#reswap_SOURCES = $(SRC)/reswap.cpp

gdr_v2_1_SOURCES = $(SRC)/gdr.v2.1.cpp

excalibur_SOURCES = $(SRC)/archive-clx.cpp $(SRC)/excalibur.cpp

cto_SOURCES = $(SRC)/cto.cpp

getharmo_SOURCES = $(SRC)/getharmo.cpp

detidor_SOURCES = $(SRC)/detidor.cpp

mgr_converter_SOURCES = $(SRC)/mgr-reconstruct.cpp $(SRC)/mgr-converter.cpp

timeseries_validate_SOURCES = $(SRC)/timeseries-validate.cpp

mgr_format_SOURCES = $(SRC)/mgr-format.cpp

mgr_editing_SOURCES = $(SRC)/mgr-editing.cpp

shorelines_doctor_SOURCES = $(SRC)/shorelines-doctor.cpp

shorelines_converter_SOURCES = $(SRC)/shorelines-converter.cpp

shorelines_smoother_SOURCES = $(SRC)/shorelines-smoother.cpp

shorelines_extractor_SOURCES = $(SRC)/fe-delauney $(SRC)/shorelines-extractor.cpp

shorelines_merge_SOURCES = $(SRC)/shorelines-merge.cpp

shorelines_distance_SOURCES = $(SRC)/shorelines-distance.cpp

SBprofiler_SOURCES = $(SRC)/SBprofiler.cpp $(SRC)/SBprofiler-io.cpp

adcp_SOURCES = $(SRC)/adcp.cpp $(SRC)/adcp-io.cpp

clock_SOURCES = $(SRC)/clock.cpp

spectral_SOURCES = $(SRC)/spectral.cpp

merge_SOURCES = $(SRC)/merge.cpp

serie_SOURCES = $(SRC)/serie.cpp

interpolation_slice_SOURCES = $(SRC)/interpolation-slice.cpp

interpolator_SOURCES = $(SRC)/interpolator.cpp

interpolator_cdf_SOURCES = $(SRC)/interpolator-cdf.cpp

fix_SOURCES = $(SRC)/fix.cpp

beta_SOURCES = $(SRC)/beta.cpp

magic_SOURCES = $(SRC)/magic.cpp

statistic1_nf_SOURCES = $(SRC)/statistic1-nf.cpp

statistic2_nf_SOURCES = $(SRC)/bmg_fe.cpp $(SRC)/statistic2-nf.cpp

statistic3_nf_SOURCES = $(SRC)/statistic3-nf.cpp

topex_SOURCES = $(SRC)/topex.cpp

detide_nf_SOURCES =$(SRC)/detide-nf.cpp

restart_nf_SOURCES = $(SRC)/restart-nf.cpp

extract_nf_v1_SOURCES = $(SRC)/extract-nf-v1.cpp

extract_nf_v2_SOURCES = $(SRC)/metrics.cpp $(SRC)/extract-nf-v2.cpp

extract_nf_v3_SOURCES = $(SRC)/extract-nf-v3.cpp

extract_sections_SOURCES = $(SRC)/metrics.cpp $(SRC)/extract-sections.cpp

showarg_SOURCES = $(SRC)/showarg.cpp

stream_2cdf_SOURCES = $(SRC)/stream-2cdf.cpp 

RDI_converter_SOURCES = $(SRC)/RDI-converter.cpp

variability_SOURCES = $(SRC)/variability.cpp

#-------------------------------------------------------------------------------
if STS
bin_PROGRAMS += stugo
STSRC=$(SRC)/sts
stugo_SOURCES = $(STSRC)/main.cpp $(STSRC)/mog2d.cpp $(STSRC)/mesh.cpp $(STSRC)/topography.cpp \
  $(STSRC)/parameters.cpp $(STSRC)/harmonic-initialise.cpp $(STSRC)/initialise2D.cpp $(STSRC)/misc.cpp $(STSRC)/regular-tides.cpp \
  $(STSRC)/boundary-init.cpp $(STSRC)/boundary2d-generic.cpp $(STSRC)/boundary2d-obsolete.cpp \
  $(STSRC)/tide-obc.cpp $(STSRC)/BC_det.cpp $(STSRC)/forcing.cpp $(STSRC)/operators.cpp \
  $(STSRC)/cefmo.cpp $(STSRC)/cefmo2D.cpp $(STSRC)/cefmo-friction.cpp $(STSRC)/cefmo2D-Q1xQ0.cpp \
  $(STSRC)/IWdrag.cpp $(STSRC)/friction.cpp $(STSRC)/friction-additional.cpp \
  $(STSRC)/configure-hugo.cpp $(STSRC)/init_config_constructor.cpp $(STSRC)/init_config_read.cpp \
  $(STSRC)/momentum2d.cpp
endif

#-------------------------------------------------------------------------------
if EXPERIMENTAL
bin_PROGRAMS += spectrum-test interp-test sts-test

spectrum_test_SOURCES = $(SRC)/spectrum-test.F90
spectrum_test_LDADD = $(LDADD) -lstdc++

interp_test_SOURCES = $(SRC)/interp-test.F90 $(SRC)/ttb-fortran.F90
interp_test_LDADD = $(LDADD) -lstdc++

sts_test_SOURCES = $(SRC)/sts-test.F90 $(SRC)/ttb-fortran.F90
sts_test_LDADD = $(LDADD) -lstdc++
sts_test_DEPENDENCIES = stugo
endif
