poc-solvers dependencies installation when you are not root
###########################################################

Summary
=======

You must have sed g++, automake, cmake, flex, bison and mercurial already installed somewhere.

The poc-solvers need SuiteSparse and LAPACK. They may install PaStiX and HIPS.
SuiteSparse needs LAPACK and BLAS. LAPACK needs BLAS.
PaStiX needs SCOTCH and hwloc if it can not have it with MPI.
HIPS needs SCOTCH.
SCOTCH needs zlib.
MUMPS 5 needs BLAS, BLACS, LAPACK, ScaLAPACK, SCOTCH 6.0.1 or above and METIS.
ScaLAPACK now contains BLACS but still needs BLAS and LAPACK.

TUGOm needs METIS when running with MPI.

Shell variables
===============

Set $ncpu to the number of cores on your machine::
  
  ncpu=`sed -re '/cpu[0-9]+ /! d' /proc/stat |wc -l`

also set $li to the directory where you can and want to install everything::
  
  li=~/.local

Note you can clean-up everything with::
  
  # rm -r $li/{lib,include,bin,share/man}/*

Make sure $li/{lib,include} are created::
  
  mkdir -p $li/{lib,include}

Standard install
================

Configure::
  
  ./configure --prefix=$li

Compile::
  
  make -j$ncpu

Install::
  
  make install

In case you want to uninstall::
  
  make uninstall

Dependencies
============

BLAS (Basic Linear Algebra Subprograms)
---------------------------------------
- Web page: http://www.netlib.org/blas/
- Tarball: http://www.netlib.org/blas/blas.tgz

Compile::
  
  make -j$ncpu

Install::

  cp -v *.a $li/lib/libblas.a

LAPACK (Linear Algebra PACKage)
-------------------------------
- Web page: http://www.netlib.org/lapack/
- Tarball: http://www.netlib.org/lapack/lapack-3.5.0.tgz (later versions may take ages)
- Documentation in README

Configure::
  
  sed -re "s|^( *BLASLIB *= *).*|\1$li/lib/libblas.a|" make.inc.example >make.inc

Compile::
  
  make -kj$ncpu;make -kj$ncpu;make

Install::
  
  cp -v liblapack.a $li/lib

PETSc
-----
Only in newer version of poc-solvers!!!
- Web page: https://www.mcs.anl.gov/petsc/
- Download page: https://www.mcs.anl.gov/petsc/download/index.html
- Tarball: http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.11.0.tar.gz
- Poor documentation in docs/installation.html

Configure, making sure the temporary directory used by Python's tempfile
is mounted with exec permissions ::
  
  # DO SEE tempfile.tempdir DOCUMENTATION FOR MORE DETAILS
  if findmnt -O noexec /tmp;then
    export TMPDIR=~/tmp
    mkdir $TMPDIR
    fi
  
  # Log file is arch-linux2-c-debug/lib/petsc/conf/configure.log
  ./configure --prefix=$li --with-shared-libraries=0

Then see `Standard install`_.

In case you want to uninstall::
  
  rm -rf $li/lib/*petsc*

ScaLAPACK (Scalable Linear Algebra PACKage)
-------------------------------------------
- Web page: http://www.netlib.org/scalapack/
- Tarball: http://www.netlib.org/scalapack/scalapack-2.0.2.tgz
- Documentation in README

Configure::
  
  cat SLmake.inc.example >SLmake.inc

Compile::
  
  make -j$ncpu;rm libscalapack.a;make

Install::
  
  cp -v libscalapack.a $li/lib

SuiteSparse
-----------
It is a suite of sparse matrix algorithms.

- Web page: http://www.SuiteSparse.com
- Tarball: http://faculty.cse.tamu.edu/davis/SuiteSparse/SuiteSparse-4.4.3.tar.gz (later versions may install improperly)
- Documentation in README.txt

Configure::
  
  cd SuiteSparse_config
  sed -re "s|^( *INSTALL_[^ ]+ *= *).*/([^/]+$)|\1$li/\2|;\
    s|^( *(LAPACK\|BLAS) *= *)|\1-L$li/lib |;" \
    SuiteSparse_config_linux.mk >SuiteSparse_config.mk
  cd ..

Compile::
  
  make -j$ncpu library

Install::
  
  make install

zlib
----
- Web page: http://www.zlib.net/
- Tarball: http://prdownloads.sourceforge.net/libpng/zlib-1.2.8.tar.gz?download
- Documentation in README

See `Standard install`_.

SCOTCH
------
- Web page: http://www.labri.fr/perso/pelegrin/scotch/
- Download page: http://gforge.inria.fr/frs/?group_id=248
- Tarball: http://gforge.inria.fr/frs/download.php/file/34618/scotch_6.0.4.tar.gz
- Documentation in INSTALL.txt

Go to the build directory::
  
  cd src/

Configure (for Mac, replace ``Makefile.inc.x86-64_pc_linux2`` by ``Makefile.inc.i686_mac_darwin10``)::
  
  sed -re "s|^( *CLIBFLAGS\t*=).*|\1 -I$li/include|;\
      s|^( *LDFLAGS\t*=)|\1 -L$li/lib|" \
      Make.inc/Makefile.inc.x86-64_pc_linux2 >Makefile.inc

Compile::
  
  make -j$ncpu ptscotch
  make -j$ncpu esmumps

Install::
  
  make prefix=$li install
  cp -v ../lib/libesmumps.a $li/lib

Uninstall with::
  
  rm $li/lib/*scotch*.a $li/lib/*scotch*.h

hwloc
-----
You should have hwloc through MPI,
but on some plateform you may need to install it for PaStiX to compile.
- Web page: https://www.open-mpi.org/projects/hwloc/
- Download page: https://www.open-mpi.org/software/hwloc/v1.11/
- Tarball: https://www.open-mpi.org/software/hwloc/v1.11/downloads/hwloc-1.11.7.tar.bz2

See `Standard install`_.

METIS
-----
- Web page: http://glaros.dtc.umn.edu/gkhome/metis/metis/overview
- Tarball: http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz
- Documentation in BUILD.txt

Configure::
  
  make config prefix=$li

Compile::
  
  make -j$ncpu

Install::
  
  make install

MUMPS: a MUltifrontal Massively Parallel sparse direct Solver
-------------------------------------------------------------
- Web and request page: http://mumps.enseeiht.fr/
- Tarball by e-mail: MUMPS_5.1.0.tar.gz
- Documentation in INSTALL

Configure::
  
  sed -re "s|^( *I[^ ]+ *= *-I)/usr/include/.+|\1$li/include|;\
      s|^( *L[^ ]+DIR *= *)/usr/lib *$|\1$li/lib|;\
      s|^( *INCPAR *=).*|\1|;\
      s|^( *CC *=).*|\1 mpicc|;\
      s|^( *F. *=).*|\1 mpifort|;\
      s|^( *SCALAP *=).*|\1 -lscalapack|;\
      "'s|^( *LIBPAR *=).*|\1 $(SCALAP) $(LAPACK)|' \
      Make.inc/Makefile.debian.PAR >Makefile.inc

Compile::
  
  make all -j$ncpu

Install::
  
  cp -v lib/*.a $li/lib
  cp -v include/*.h $li/include/

poc-solvers
===========
Configure and compile::
  
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$li/lib
  export CMAKE_INCLUDE_PATH=$li/include
  ./INSTALL.sh -s $li

You can also set the environment with modules::
  
  module use -a $PWD
  module load modulefile

To make a small distribution::
  
  tar czvf poc-solvers.tar.gz --exclude=\*.kdev4 DEPENDENCIES_NOT_ROOT.rst README INSTALL.sh src/
