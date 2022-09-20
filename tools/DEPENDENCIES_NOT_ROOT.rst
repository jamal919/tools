Sirocco tools dependencies installation when you are not root
#############################################################

Summary
=======

The SIROCCO Tools need at least proj-4.9.1, GSL, NetCDF and LAPACK.
They will get LAPACK through the poc-solvers.
NetCDF needs HDF5.

The SIROCCO Tools optionally need Shapelib and at least grib_api-1.20.0 .

Shell variables
===============

Set $ncpu to the number of cores on your machine::
  
  ncpu=`sed -re '/cpu[0-9]+ /! d' /proc/stat |wc -l`

also set $li to the directory where you can and want to install everything::
  
  li=~/.local

Avoid lib{,64} confusion::
  
  ln -s lib $li/lib64

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

PROJ.4 (Cartographic Projections Library)
-----------------------------------------
- Web page: https://trac.osgeo.org/proj/
- Download page: https://proj4.org/download.html#download
- Tarball: http://download.osgeo.org/proj/proj-5.1.0.tar.gz

See `Standard install`_.

GSL (GNU Scientific Library)
----------------------------
- Web page: http://www.gnu.org/software/gsl/
- Tarball: http://ftpmirror.gnu.org/gsl/gsl-1.16.tar.gz

See `Standard install`_.

HDF5
----
- Web page: http://www.hdfgroup.org/HDF5/
- Download page: http://www.hdfgroup.org/HDF5/release/obtainsrc.html
- Tarball: http://www.hdfgroup.org/ftp/HDF5/current/src/hdf5-1.8.14.tar.gz
- Documentation in release_docs/INSTALL

See `Standard install`_.

Network Common Data Form (NetCDF)
---------------------------------
- Web page: http://www.unidata.ucar.edu/software/netcdf/
- Download page: http://www.unidata.ucar.edu/downloads/netcdf/
- Tarball: ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4.3.3.1.tar.gz

Configure::
  
  ./configure --prefix=$li LDFLAGS="$LDFLAGS -L$li/lib" CPPFLAGS="$CPPFLAGS -I$li/include"

Then see `Standard install`_.

GRIB_API
--------
- Web page: https://software.ecmwf.int/wiki/display/GRIB/Home
- Download page: https://software.ecmwf.int/wiki/display/GRIB/Releases
- Tarball: https://software.ecmwf.int/wiki/download/attachments/3473437/grib_api-1.20.0-Source.tar.gz
- Documentation: https://software.ecmwf.int/wiki/display/GRIB/GRIB+API+CMake+installation

Configure::
  
  mkdir build
  cd build
  CMAKE_LIBRARY_PATH=$LD_LIBRARY_PATH cmake .. -DCMAKE_INSTALL_PREFIX=$li

Then see `Standard install`_.

Shapelib
--------
- Web page: http://shapelib.maptools.org/
- Download page: http://download.osgeo.org/shapelib/
- Tarball: http://download.osgeo.org/shapelib/shapelib-1.5.0.tar.gz

See `Standard install`_.

Tools
=====

Configure::
  
  autoreconf -si
  # Mac users may need to add ``CXX=g++-7`` to the line below
  ./configure LDFLAGS="$LDFLAGS -L$li/lib" CPPFLAGS="$CPPFLAGS -I$li/include"

Compile::
  
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$li/lib
  make -j$ncpu -C objects/ ../src/tools-fortran-sizes.def
  make -j$ncpu
