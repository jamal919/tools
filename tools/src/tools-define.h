
/**************************************************************************

  T-UGO tools, 2006-2011

  Unstructured Ocean Grid initiative

***************************************************************************/
/** \file

\author  Florent Lyard      LEGOS/CNRS, Toulouse, France. florent.lyard@legos.obs-mip.fr
\author  Laurent Roblou     LEGOS/CNRS, Toulouse, France
\author  Damien Allain      LEGOS/CNRS, Toulouse, France
\author  Yves Soufflet      LEGOS/CNRS, Toulouse, France
\author  Clément Mayet      LEGOS, Toulouse, France (PhD)
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada

\brief global defines and a few constants. Also contains the \ref index "main page".

*/
/*----------------------------------------------------------------------------*/

/** \mainpage TUGOm tools

Version and revision numbers available respectively from #VERSION and #REVISION.

\section shortcuts Direct links to usefull pages
<table><tr>
<td><a href=globals_func.html><b>Alphabetical list of functions</b></a>
<td><a href=todo.html><b>Todo list</b></a>
<td><a href=bug.html><b>Bug list</b></a>
</table>

\section tipstoc Coding tips
- \ref oneLineRead
- \ref oneLineSave

\section budget energy budget
See tides-energy.cpp comodo-energy.cpp

\section detiding
See altimetry-detidor.cpp comodo-detidor.cpp detide-nf.cpp detidor.cpp mercator-detidor.cpp comodo-filter.cpp showarg.cpp

\subsection prediction prediction and atlases
See academics.cpp predictor.cpp comodo-admittance.cpp

\subsection filtering
See comodo-sfilter.cpp filter-02.cpp filter.cpp

\section libraries
TUGOm :
- uses bathymetries and meshes as input
- and produces tidal atlases as output.
\dot
digraph TUGOm {
  rankdir=LR;
  { bathymetry mesh } -> TUGOm -> model;
  model [label="tidal models"];
  TUGOm [shape=box];
  }
\enddot

The tools libraries must therefore provide functionalities to manipulate :
- bathymetries
- meshes
- tidal elevation maps
- tidal gauges data.

This implies interpolation and, if there is a solver available, projection and spectral analysis.
\dot
digraph tools {
  node [shape=box];
  netcdf -> time; netcdf [label="poc-netcdf"];time [label="poc-time"];
  gauges -> spectral;
  { projection spectral } -> solver; spectral [label="spectral analysis"];
  grid -> interpolation -> beta; grid [label="grid conversion"];
  }
\enddot

\section notes General notes on the sources

\todo remove \c malloc
\todo replace
\code
  exitIfNull(
    depths=(double *) malloc(grid.nz*sizeof(double))
    );
\endcode
with
\code
  exitIfNull(
    local=new T[grid.nz]
    );
\endcode

\note perlReplace.pl is usefull.
So is :
\code
l='/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*\/'
l+'\n\n'+txt.value+'\n'+l+'\n{return '+txt.value.replace(/([\(,]) *.*?[ \*]([a-z][a-z0-9]*)/g,'$1$2')+';}'+'\n\n\n'
\endcode

\note keep tags up to date :
\code
hg grep --all AC_INIT configure.ac|sed -re 's/.*:([0-9]+):\+.*, ([0-9][^,]*),.*./hg tag -r \1 \2/;t;d'
\endcode

\note avoid spagetti
\code
[ $(hg st -mardn `hg in -q --template '{files} '` ) ] || hg pull -u
\endcode

\bug Some comments are with the wrong encoding. Look for '�' in the code

\bug Some array of classes may be \c malloc 'ed. This is wrong : they must be \c new 'ed instead.
\todo Find \c malloc 'ed array of classes with :
\code
classesList=`sed -re 's/^ *class +([0-9A-Za-z_]+).*$/\1/;t;d' src/*.* |sort -u`
classesList=`echo $classesList |sed -e 's/ /|/g'`
for class in `sed -re "s/.*\b($classesList)\b.+\bmalloc\b.*$/\1/;t;d" src/*.cpp|sort -u`;do echo "$class\b.+\bmalloc";done
\endcode

\note keep includes clean, in particular DO NOT USE malloc.h USE stdlib.h INSTEAD
\code
for f in !(cfortran).h{,pp};do echo "#include \"$f\""|g++ -E -x c++ - &>/dev/null || echo $f;done
\endcode

\todo Make sure all main() call fct_echo() with :
\code
grep -L fct_echo `grep -lE '\bmain\b *\(' *.cpp`
\endcode

\todo 2012-01-24 Damien Allain : better sort functions

\todo Find obsoleted functions and variable with :
\code
v=cdfvar_t::dim[^’]*
make clean;true>nohup.out;nohup make -kj8
sed -e 's/: warning: ‘'"$v"'’ is deprecated.*$//;t;d' nohup.out|sort -t : -k 1,1 -k 2,2n -u
\endcode

\note Find values of PI with :
\code
sed -re 's/.*[^0-9](3\.14[0-9]*|acos\(-1\.?0*\)).*?/\1/;t;d' *.*|sort -u
\endcode
and values of anything with :
\code
sed -re 's/^.*[^\.e0-9](([0-9]+(\.[0-9]*)?|[0-9]*\.[0-9]+)([eE][-+]?[0-9]+)?).*?$/\1/;t;d' *.cpp *.h *.def |sort -ug
\endcode

\note list options with:
\code
sed -re 's|.*keyword,"([^"]+)".*|    "  \1 : \\n"|;t;s|^\s+case '"'"'(.)'"'"'\s*:\s*|    "  -\1 : \\n"|;t;d' ....cpp
\endcode

\note check for not so usefull files with :
\code
for f in `hg st -nmc !(*.py|*.sh|perlReplace.pl|toDelete*)`;do grep -qE "^[^#]*$f\b" ../objects/Makefile.am || echo $f;done
for f in *.h;do grep -q $f *.cpp !(${f/%.h}).h || echo $f;done
\endcode

\note check for files to clean-up with :
\code
make clean
nohup make -kj12 "AM_CFLAGS = -Wunused -Wuninitialized -Wno-unused-but-set-variable" >nohup.out
perl -e 'while(<>){s/(:\d+){2}: warning: .*$// or next;chomp;$c{$_}++}foreach $f(keys %c){print "$c{$f} $f\n"}' nohup.out|sort -g
\endcode

\note
- ellipse.cpp and tides-sample.cpp and ellipse_Madmp() and ellipse_parameter() may be similar
- SBprofiler.cpp and adcp.cpp are identical !
- SBprofiler-io.cpp and adcp-io.cpp have many similarities
- topo-format.cpp is almost identical to topo-merge.cpp, but both only have a main()
- MSS : mean sea surface stuff, may be used by CTOH
- mss-io.cpp and mss-io.v2.cpp have too many similar or identical functions
- mss-comp.v2.cpp and mss.v1.cpp have many similarities
- fe2obc.cpp and atlas2obc.cpp have 2 functions, that are both similar
- array_index() and pos() are NOT identical : the argument order is different

\todo Make a wiki to find out what these programmes are for :
adcp.cpp
SBprofiler.cpp
altimetry-detidor.cpp
altimetry-error.cpp
analysis-ecoop.cpp
analysis-tugo.cpp
atlas2obc.cpp
beta.cpp
clock.cpp
convert3d.cpp
cto.cpp
detide-nf.cpp
excalibur.cpp
extract-nf-v1.cpp
extract-nf-v2.cpp
extract-nf-v3.cpp
fe2obc.cpp
fe-format.cpp
fes-assembly.cpp
fix.cpp
format.cpp
gdr.v2.1.cpp
getharmo.cpp
gmt-converter.cpp
gridit-cdf.cpp
gridit-cls.cpp
gridit.cpp
gridit-surfref.cpp
haw2gnu.cpp
interpolator-cdf.cpp
interpolator.cpp
kelvin.cpp
magic.cpp
mapxyz-v0.cpp
mapxyz-v1.cpp
mercator-detidor.cpp
mesh2mesh.cpp
mesh2node.cpp
mesh2restart2d.cpp
mesh-assembly.cpp
mesh-correlation.cpp
mesh-critere.cpp
mesh-doctor.cpp
mesh-generator.cpp
mesh-partition.cpp
mesh-refine.cpp
mesh-renum.cpp
mesh-split.cpp
mesh-topo.cpp
mesh-upgrade.cpp
mesh-zmin.cpp
mgr2obc.cpp
mog3d-restart.cpp
mss-main.v2.cpp
mss.v1.cpp
operation.cpp
penalty.cpp
RDI-converter.cpp
rectify.cpp
regular2node-3d.cpp
regular2node.cpp
restart-nf.cpp
sample.cpp
search_nodes.cpp
sectionxy.cpp
sectionxz.cpp
serie.cpp
spectral.cpp
statexyz.cpp
statistic1-nf.cpp
statistic2-nf.cpp
statistic3-nf.cpp
stream-2cdf.cpp
symharmo.cpp
symmic.cpp
symsolve.cpp
symtides.cpp
symtools.cpp
synthetic.cpp
tides-compare.cpp
tides-sample.cpp
topex.cpp
topo-create.cpp
topo-format.cpp
topo-merge.cpp
topo-smooth.cpp
triangle.cpp
variability.cpp
ww3tools.cpp
xover.cpp
xtrack-format.cpp

\todo 2012-04-04,2013-07-29 Damien Allain : get rid of ATLAS and of LINPACK*, stick these solvers paragraphs in a proper function in poc-solvers
\note 2012-04-04 Damien Allain : think about getting rid of comodo patches

 */

 
#ifndef TOOLS_DEFINE
#define TOOLS_DEFINE 1

//#define  Boolean char //already in tools-structures.h
#define True  1
#define False 0

#define TIMESRIES_FORMAT_MATROOS 0
#define TIMESRIES_FORMAT_CLASSIC 1


/*----------------------------------------------------------------------------*/
/// version of cpp/gcc
/**
copied from info:/cpp/Common Predefined Macros
*/
/*----------------------------------------------------------------------------*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
#define GCC_VERSION (__GNUC__*10000 + __GNUC_MINOR__*100 + __GNUC_PATCHLEVEL__)
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/


/*----------------------------------------------------------------------------*/
/// safe deprecated attribute, with, if your version of gcc supports it, a descriptive message
/**
See info:gcc/Function%20Attributes
and http://gcc.gnu.org/gcc-4.5/changes.html
*/
/*----------------------------------------------------------------------------*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
#if GCC_VERSION >= 40500  /* Test for GCC >=4.5.0 */
#define safe_deprecated(x) deprecated(x)
#else
#define safe_deprecated(x) deprecated
#endif
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/


/*----------------------------------------------------------------------------*/
/// used only by #TOSTRING
/**
See info:/cpp/Stringification
*/
/*----------------------------------------------------------------------------*/
#define STRINGIFY(x) #x
/*----------------------------------------------------------------------------*/
/// stringifies
/**
Useful for assertions. See info:/cpp/Stringification and #STRINGIFY
*/
/*----------------------------------------------------------------------------*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
#define TOSTRING(x) STRINGIFY(x)
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

#define MAX(a,b)  (((a) > (b)) ? (a) : (b))
#define MIN(a,b)  (((a) < (b)) ? (a) : (b))

#endif
