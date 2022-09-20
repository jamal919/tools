FROM ubuntu:18.04

ENV TZ=Europe/Paris
RUN ln -snf /usr/share/zoneinfo/${TZ} /etc/localtime && echo ${TZ} > /etc/timezone

RUN apt-get update \
    && apt-get install -y wget mercurial make automake cmake gcc g++ gfortran flex bison lbzip2 pkg-config \
    && apt-get install -y python3 python3-distutils netcdf-bin libnetcdf-dev libnetcdff-dev openmpi-bin libopenmpi-dev

RUN mkdir /work
ENV PATH=${PATH}:/usr/local/bin
ENV LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/usr/local/lib

RUN cd /work \
    && wget http://www.netlib.org/blas/blas-3.10.0.tgz \
    && tar -xvzf blas-3.10.0.tgz \
    && cd BLAS-3.10.0 \
    && make \
    && cp -v *.a /usr/local/lib/libblas.a \
    && rm -rf /work/blas-3.10.0.tgz /work/BLAS-3.10.0

RUN cd /work \
    && wget -c https://github.com/Reference-LAPACK/lapack/archive/refs/tags/v3.10.1.tar.gz \
    && tar -xvzf v3.10.1.tar.gz \
    && cd lapack-3.10.1 \
    && sed -re "s|^( *BLASLIB *= *).*|\1/usr/local/lib/libblas.a|" make.inc.example > make.inc \
    && make \
    && cp -v lib*.a /usr/local/lib \
    && rm -rf /work/v3.10.1.tar.gz /work/lapack-3.10.1

RUN cd /work \
    && ln -sf /usr/bin/python3 /usr/bin/python \
    && wget -c http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.11.0.tar.gz \
    && tar -xvzf petsc-3.11.0.tar.gz \
    && cd petsc-3.11.0 \
    && ./configure --with-shared-libraries=0 --with-debugging=0 CPPFLAGS="-I/usr/local/include" COPTFLAGS="-O3" CXXOPTFLAGS="-O3" FOPTFLAGS="-O3" --prefix=/usr/local \
    && make PETSC_DIR=/work/petsc-3.11.0 PETSC_ARCH=arch-linux-c-opt all \
    && make PETSC_DIR=/work/petsc-3.11.0 PETSC_ARCH=arch-linux-c-opt install \
    && rm -rf /work/petsc-3.11.0.tar.gz /work/petsc-3.11.0

RUN cd /work \
    && wget -c  http://www.netlib.org/scalapack/scalapack-2.0.2.tgz \
    && tar -xvzf scalapack-2.0.2.tgz \
    && cd scalapack-2.0.2 \
    && mkdir build; cd build; cmake ../ \
    && make \
    && cp -v lib/libscalapack.a /usr/local/lib \
    && rm -rf /work/scalapack-2.0.2.tgz /work/scalapack-2.0.2

RUN cd /work \
    && wget -c http://faculty.cse.tamu.edu/davis/SuiteSparse/SuiteSparse-4.4.3.tar.gz \
    && tar -xvzf SuiteSparse-4.4.3.tar.gz \
    && cd SuiteSparse/SuiteSparse_config \
    && sed -re "s|^( *INSTALL_[^ ]+ *= *).*/([^/]+$)|\1/usr/local/\2|;s|^( *(LAPACK\|BLAS) *= *)|\1-L/usr/local/lib |;" SuiteSparse_config_linux.mk > SuiteSparse_config.mk \
    && cd ../ \
    && make library \
    && make install \
    && rm -rf /work/SuiteSparse-4.4.3.tar.gz /work/SuiteSparse

# https://sources.easybuild.io/s/SCOTCH/scotch_6.0.4.tar.gz
COPY ./scotch_6.0.4 /work/scotch_6.0.4
RUN cd /work/scotch_6.0.4/src \
    && sed -re "s|^( *CCD\t*=).*|\1 mpicc|; s|^( *CLIBFLAGS\t*=).*|\1 -I/usr/local/include|; s|^( *LDFLAGS\t*=)|\1 -L/usr/local/lib|" Make.inc/Makefile.inc.x86-64_pc_linux2 > Makefile.inc \
    && make ptscotch \
    && make esmumps \
    && make prefix=/usr/local install \
    && cp -v ../lib/libesmumps.a /usr/local/lib \
    && rm -rf /work/scotch_6.0.4.tar.gz /work/scotch_6.0.4

RUN cd /work \
    && wget -c http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz \
    && tar -xvzf metis-5.1.0.tar.gz \
    && cd metis-5.1.0 \
    && make config prefix=/usr/local \
    && make \
    && make install \
    && rm -rf /work/metis-5.1.0.tar.gz /work/metis-5.1.0

COPY ./MUMPS_5.1.1 /work/MUMPS_5.1.1
RUN cd /work/MUMPS_5.1.1 \
    && sed -re "s|^( *I[^ ]+ *= *-I)/usr/include/.+|\1/usr/local/include|; \
    s|^( *L[^ ]+DIR *= *)/usr/lib *$|\1/usr/local/lib|; \ 
    s|^( *INCPAR *=).*|\1|; \
    s|^( *CC *=).*|\1 mpicc|; \
    s|^( *F. *=).*|\1 mpifort|; \
    s|^( *SCALAP *=).*|\1 -lscalapack|; \
    s|^( *LIBPAR *=).*|\1 \$(SCALAP) \$(LAPACK)|" \
    Make.inc/Makefile.debian.PAR > Makefile.inc \
    && make \
    && cp -v lib/*.a /usr/local/lib \
    && cp -v include/*.h /usr/local/include \
    && rm -rf /work/MUMPS_5.1.1

# ftp://ftp.legos.obs-mip.fr/pub/ecola/poc-solvers/poc-solvers.hg_bundle - rev.388
# https://sources.easybuild.io/p/PaStiX/pastix_5.2.2.22.tar.bz2
# sed -i 's,wget https://gforge.inria.fr/frs/download.php/file/$remodeDirName/$baseName.tar.bz2,wget --no-check-certificate https://sources.easybuild.io/p/PaStiX/pastix_5.2.2.22.tar.bz2,g' INSTALL.sh
COPY ./poc-solvers /work/poc-solvers
COPY ./pastix_5.2.2.22.tar.bz2 /work/pastix_5.2.2.22.tar.bz2

RUN ln -sf /usr/lib/x86_64-linux-gnu/libX11.so.6 /usr/local/lib/libX11.so \
    && cd /work/poc-solvers \
    && export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib:/usr/lib64 \
    && export CMAKE_INCLUDE_PATH=/usr/local/include \
    && ./INSTALL.sh -s /usr/local \
    && rm -rf /work/pastix_5.2.2.22.tar.bz2
    
RUN cd /work \
    && wget -c https://download.osgeo.org/proj/proj-5.1.0.tar.gz \
    && tar -xvzf proj-5.1.0.tar.gz \
    && cd proj-5.1.0 \
    && ./configure --prefix=/usr/local \
    && make \
    && make install \
    && rm -rf /work/proj-5.1.0.tar.gz /work/proj-5.1.0

RUN cd /work \
    && wget -c https://ftp.gnu.org/gnu/gsl/gsl-1.16.tar.gz \
    && tar -xvzf gsl-1.16.tar.gz \
    && cd gsl-1.16 \
    && ./configure --disable-shared --prefix=/usr/local \
    && make \
    && make install \
    && rm -rf /work/gsl-1.16.tar.gz /work/gsl-1.16

RUN cd /work \
    && wget -c http://download.osgeo.org/shapelib/shapelib-1.5.0.tar.gz \
    && tar -xvzf shapelib-1.5.0.tar.gz \
    && cd shapelib-1.5.0 \
    && ./configure --prefix=/usr/local \
    && make \
    && make install \
    && rm -rf /work/shapelib-1.5.0.tar.gz /work/shapelib-1.5.0

# ftp://ftp.legos.obs-mip.fr/pub/ecola/tools/tools.hg_bundle - rev.3873
# Removed reference to "undefined" SOLVER_ID_GMRES, and removes pdf making (/data) in Makefile.am
COPY ./tools /work/tools
RUN cd /work/tools \
    && autoreconf -si \
    && ./configure --prefix=/usr/local POCSOLVERDIR="../poc-solvers" LDFLAGS="$LDFLAGS -L/usr/local/lib" CPPFLAGS="$CPPFLAGS -I/usr/local/include" CFLAGS="-fpermissive -O2" \
    && make -k -j8 \
    && make \
    && make install \
    && rm -rf /work

RUN apt-get clean

WORKDIR /mnt

CMD ["/bin/bash"]
